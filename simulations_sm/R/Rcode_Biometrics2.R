


#====================================================#
#=== step 1:  the main code =========================#
#====================================================#
#This code very ggood results for n=200.  
##remove.packages("hal9001")
#devtools::install_github("tlverse/hal9001", build_vignettes = FALSE)
library(hal9001)
library(SuperLearner)
library(randomForest)

#rm(list=ls())
#======================================================================#
#======== test the effect of undersmooth HAL in one stage design ======#
#======================================================================#

### data generation

### ps = 6, Comparison with augmented inverse probability weighted estimators Section S3.4



dat_hal <- function(n, ps){
  
  if(ps == 6){### 
    x1 = runif(n,-2,2)
    x2 = rnorm(n,0,0.5)
    p = exp((2+2*x1-2*x2+x1*x2)/4) / (1 + exp((2+2*x1-2*x2+x1*x2)/4)) ; summary(p)
    A = rbinom(n, size = 1, prob = p)
    Y = +0*A +x1+1*x1^3 + rnorm(n, mean = 0, sd = 0.1)
    
  }

    Q0<-1
  return(data.frame(Y, A,Q0,p,x2, x1))
}




onestage_under_hal <- function(dat, rate = 0.8,nfolds=5, ss=2,trunc=0.05){
  
  library(hal9001)
  if(ss==1) {true=0.95}
  if(ss==2) {true=1.4}
  if(ss==3) {true=0.35}
  if(ss==4) {true=0}
  if(ss==5) {true=0}
  if(ss==6) {true=0}
  n=nrow(dat)
  #nfolds=10
  folds <- cut(seq(1,nrow(dat)),breaks=nfolds,labels=FALSE)
  preds_rf<-dat$A
  coef_cv<-vector(mode = "list", length = nfolds)
  
  for(i in 1:nfolds){
    testIndexes <- which(folds==i,arr.ind=TRUE)
    test <- dat[testIndexes,]
    train <- dat[-testIndexes,]
    mod_cv = fit_hal(X = train[,-c(1:4)], Y = train$A, family = "binomial", yolo = FALSE, return_x_basis = T,max_degree = 3)
    preds_rf[testIndexes] = predict(mod_cv, new_data = test[,-c(1:4)])
    
  }
  preds_rf = preds_rf
  preds_rf[preds_rf > 0.98] = 0.98
  preds_rf[preds_rf < 0.02] = 0.02
  #======================================#
  #====== step 1: baseline model ========#
  #======================================#
  
  ### logistic: remember to change the covariate if you change the data generation function
  mod_lr = glm(A~x2, data = dat, family = "binomial")
  preds_lr = predict(mod_lr, new_data = dat[,-c(1:4)], type = "response")
  
 
  ### ipw estimator
  rf_ipw = sum(dat$A*dat$Y/preds_rf)/sum(dat$A/preds_rf)
  lr_ipw = sum(dat$A*dat$Y/preds_lr)/sum(dat$A/preds_lr)
  unadj_est = mean(dat$Y[dat$A==1])
  
  ### se estimation for hal
  #cv_se_half1 = sqrt(sum((dat$A[s]*(dat$Y[s] - cv_ipw)/preds_cv[s] /sum(dat$A[s]/preds_cv[s]) * n)^2)/n^2); cv_se_half1
  lr_se = sqrt(sum((dat$A*(dat$Y - lr_ipw)/preds_lr /sum(dat$A/preds_lr) * n)^2)/n^2)
  
  
  Q1<-Q0<-dat$A
  for(i in 1:nfolds){
    testIndexes <- which(folds==i,arr.ind=TRUE)
    test <- dat[testIndexes,]
    train <- dat[-testIndexes,]
    Ymod_cv = tryCatch(fit_hal(X = train[,c("A","x1","x2")], Y = train$Y, yolo = FALSE, max_degree = 3), error=function(e) e, warning=function(w) w)
    ntest<-test[,c("A","x1","x2")]
    ntest$A<-1
    Q1[testIndexes] = predict(Ymod_cv, new_data = ntest)
    
 
  }

  ndat<-dat[,c("A","x1")]
  ndat$A<-1
  Q.lm.t<-(lm(Y~A+x1,data=dat))
  # Q.lm.t<-(lm(Y~A+I(x2^2)+x1*x2,data=dat))
  Q1.ymiss<-predict(Q.lm.t, new_data = ndat)
  
  true.ps<-dat$p
  dr.est<-mean(dat$A*dat$Y/preds_rf 
                     - (dat$A - preds_rf)/preds_rf *Q1)
  
  dr.est.dmiss<-mean(dat$A*dat$Y/preds_lr 
               - (dat$A - preds_lr)/preds_lr *Q1)
  
  dr.est.ymiss<-mean(dat$A*dat$Y/preds_rf 
                     - (dat$A - preds_rf)/preds_rf *Q1.ymiss)
  
  sigma.star<- sqrt(mean((dat$A*dat$Y/preds_rf 
                          - (dat$A - preds_rf)/preds_rf *Q1 - dr.est)^2)/n)  
  sigma.star.dmiss<- sqrt(mean((dat$A*dat$Y/preds_lr 
                                - (dat$A - preds_lr)/preds_lr *Q1 - dr.est.dmiss)^2)/n)  
  sigma.star.ymiss<- sqrt(mean((dat$A*dat$Y/preds_rf 
                                - (dat$A - preds_rf)/preds_rf *Q1.ymiss - dr.est.ymiss)^2)/n)  
  
  D_CAR_dr<- abs(t(Q1) %*% as.matrix((dat$A - preds_lr)/preds_lr)/n)
  
  
  cp_dr = mean(dr.est - 1.96*sigma.star <= true & dr.est + 1.96*sigma.star >= true)
  cp_dr_dmiss = mean(dr.est.dmiss - 1.96*sigma.star <= true & dr.est.dmiss + 1.96*sigma.star >= true)
  cp_dr_ymiss = mean(dr.est.ymiss - 1.96*sigma.star <= true & dr.est.ymiss + 1.96*sigma.star >= true)
  cp_lr = mean(lr_ipw - 1.96*lr_se <= true & lr_ipw + 1.96*lr_se >= true)
  
  
  dr_est_vec<-c(dr.est,dr.est.dmiss,dr.est.ymiss)
  names(dr_est_vec)<-c("dr_hal","dr_dmiss","dr_ymiss")
  bias<-( dr_est_vec-true)
  rnbias<-sqrt(n) * ( bias)
  names(rnbias)<-c("rnbias_dr_hal","rnbias_dr_dmiss","rnbias_dr_ymiss")
  sigma.star.all<-c(sigma.star,sigma.star.dmiss,sigma.star.ymiss)
  names(sigma.star.all)<-c("sigma.star","sigma.star.dmiss","sigma.star.ymiss")
  
  
  rnmse<-sqrt(n)*(bias^2+sigma.star.all^2)
  names(rnmse)<-c("rnmse_dr_hal","rnmse_dr_dmiss","rnmse_dr_ymiss")
  rndcar_dr<-sqrt(n)*D_CAR_dr
  names(rndcar_dr)<-"rndcar_dr"
  
  cp.drs<-c(cp_dr,cp_dr_dmiss,cp_dr_ymiss)
  names(cp.drs)<-c("cp_dr_hal","cp_dr_dmiss","cp_dr_ymiss")
  
  result = c(unadj_est=unadj_est,lr_ipw=lr_ipw,lr_se=lr_se,cp_lr=cp_lr,  
             dr_est_vec,cp.drs, 
             sigma.star.all=sigma.star.all, 
             rnbias,rnmse,rndcar_dr)
  return(result)
}


nsim=48*5
sim.1 <- vector("list", length = nsim)
rate=0.95
nfolds=5
trunc=0.00

set.seed(123)
n = 5000
ss=6
for(i in 1:nsim){
  dat<-dat_hal(n = n, ps = ss )
  sim.1[[i]]<-dat
}
### example
#onestage_under_hal(dat = sim.1[[2]], rate = 0.95,nfolds=nfolds,ss=ss)


library(doParallel)#Load parallel library
no_cores <- detectCores()-0 #Set number of cores
cl <- makeCluster(no_cores) 
registerDoParallel((cl))

results.list <- foreach(k =1:length(sim.1), .packages = c("randomForest","SuperLearner","hal9001","MASS","sandwich")) %dopar%  
  tryCatch(onestage_under_hal(dat = sim.1[[k]], rate = rate,nfolds=nfolds,ss=ss,trunc=trunc), error = function(e) {print(e); NA})


stopCluster(cl)
round(apply(do.call(rbind,results.list),2,function(x) mean(x,na.rm=TRUE)),3)



