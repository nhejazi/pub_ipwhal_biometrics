#====================================================#
#=== step 1:  the main code =========================#
#====================================================#
##remove.packages("hal9001")
#devtools::install_github("tlverse/hal9001", build_vignettes = FALSE)
library(hal9001)
#rm(list=ls())
#======================================================================#
#======== test the effect of undersmooth HAL in one stage design ======#
#======================================================================#

### data generation

### ps = 1, randomized design Section S3.1
### ps = 2, Observational study without positivity violations Section S3.2
### ps = 3, Observational study with near-positivity violations Section S3.3
### ps = 4, Scenario 2 Section 5
### ps = 5, Scenario 1 Section 5
### ps = 6, Comparison with augmented inverse probability weighted estimators Section S3.4
### ps = 7, Observational study with 10 covariates Section S3.5

dat_hal <- function(n, ps){
  if(ps == 1){### randomized trial
    x1 = runif(n,0.2,0.8)
    x2 = rbinom(n,1,0.3)
    p = 0.5
    A = rbinom(n, size = 1, prob = p)
    Y = A*(x1+x2+x1*x2)+(1-A)*x1 + rnorm(n, mean = 0, sd = 0.1)
  }
  if(ps == 2){###
    x1 = runif(n,0.2,0.8)
    x2 = rbinom(n,1,0.6)
    p = exp(2*x1-x2-x1*x2) / (1 + exp(2*x1-x2-x1*x2)) ; summary(p)
    A = rbinom(n, size = 1, prob = p)
    Y = A*(x1+x2+x1*x2)+(1-A)*x1 + rnorm(n, mean = 0, sd = 0.1)
  }
  if(ps == 3){###
    x1 = runif(n,0.0,0.6)
    x2 = rbinom(n,1,0.05)
    p = exp(2*x1-4*x2+x1*x2) / (1 + exp(2*x1-4*x2+x1*x2)) ; summary(p)
    A = rbinom(n, size = 1, prob = p)
    Y = A*(x1+x2)+(1-A)*x1 + rnorm(n, mean = 0, sd = 0.1)
  }
  if(ps == 4){###
    x1 = runif(n,-2,2)
    x2 = rnorm(n,0,0.5)
    p = exp((1*x2^2-exp(x1/2)-0*x1*x2)/2) / (1 + exp((1*x2^2-exp(x1/2)-0*x1*x2)/2)) ; summary(p)
    A = rbinom(n, size = 1, prob = p)
    Y = +0*A +(-2*x2^2+2*x1+2*mean(x2^2)+x2+x1*x2)/1 + rnorm(n, mean = 0, sd = 0.1)
  }
  if(ps == 5){###
    x1 = runif(n,-2,2)
    x2 = rnorm(n,0,0.5)
    p = exp((1.5*x2-x1)/2) / (1 + exp((1.5*x2-x1)/2)) ; summary(p)
    A = rbinom(n, size = 1, prob = p)
    Y = +0*A -x2/1.5+x1/2 + rnorm(n, mean = 0, sd = 0.1)
  }
  if(ps == 6){### 
    x1 = runif(n,-2,2)
    x2 = rnorm(n,0,0.5)
    p = exp((2+2*x1-2*x2+x1*x2)/4) / (1 + exp((2+2*x1-2*x2+x1*x2)/4)) ; summary(p)
    A = rbinom(n, size = 1, prob = p)
    Y = +0*A +x1+1*x1^3 + rnorm(n, mean = 0, sd = 0.1)
    
  }
  if(ps == 7){### 
    x1 = runif(n,-2,2)
    x2 = runif(n,-2,2)
    x3 = runif(n,-2,2)
    x4 = runif(n,-2,2)
    x5 = runif(n,-2,2)
    x6 = rbinom(n,1,0.6)
    x7 = rbinom(n,1,0.6)
    x8 = rbinom(n,1,0.6)
    x9 = rbinom(n,1,0.6)
    x10 = rbinom(n,1,0.6)
    
  
    p = exp((1*x2^2-exp(x1/2)-x3+x4-exp(x5/2)+x6+x7)/2) / (1 + exp((1*x2^2-exp(x1/2)-x3+x4-exp(x5/2)+x6+x7)/2)) ; summary(p)
    A = rbinom(n, size = 1, prob = p)
    Y = +0*A +(-2*x2^2+2*x1+2*mean(x2^2)+x2+x1*x2+x3+x4+2*x5^2-2*mean(x5^2))/1 + rnorm(n, mean = 0, sd = 0.1)
  }
  Q0<-1
  if(ps != 7) {return(data.frame(Y, A,Q0,p,x2, x1))} else {
    return(data.frame(Y, A,Q0,p,x1, x2, x3, x4, x5, x6, x7, x8, x9, x10))}
}



#====================== input ==================================#
#=== this function require HAL 0.2.2 ===========================#
# n : the sample size                                           #
# ps : the ps model type in data generation process             #
# seed : the random seed                                        #
# rate : how fast we decrease the lambda                        #
#===============================================================#

#===================== output ==================================#
# lr_ipw : the ipw estimate use logistic regression             #
# cv_ipw : the ipw estimate use cv hal                          #
# cv_ipw_trun:  the truncated ipw estimate use cv hal           #
# cv_se: the se estimation for cv hal                           #
# cv_se_trun : the se estimation for truncated cv hal           #
# preds_lr: fitted prob for logistic regression                 #
# preds_cv : fitted prob for cv hal                             #
# preds_cv_trun : truncated fitted prob for cv hal              #
# rnbias_DCAR1 : the root n bias of ipw estimate use under hal  #
#                                       with DCAR criterion     #
# rnbias_SCORE:  the root n bias of ipw estimate use under hal  #
#                                       with SCORE criterion    #
# rnmse: root n MSE of different estimators                     #
# sigma.star: variane of efficient estimators based on          #
#                               different undersmooting levels  #
# cp_under: coverage of confidence intervals                    #
# preds_under : fitted prob for under hal                       #
# preds_under_trun : truncated fitted prob for under hal        #
# flag : 0, normal under lambda; 1, cv_lambda is already below  #
#        the cutoff; 2, some bad fit and set under_lambda =     #
#        cv_lambda; 3, something else                           #
# cv_lambda : omit                                              #
# under_lambda : omit                                           #
#===============================================================#

onestage_under_hal <- function(dat, rate = 0.8,nfolds=5, ss=2,trunc=0.05){
  
  library(hal9001)
  if(ss==1) {true=0.95}
  if(ss==2) {true=1.4}
  if(ss==3) {true=0.35}
  if(ss==4) {true=0}
  if(ss==5) {true=0}
  if(ss==6) {true=0}
  if(ss==7) {true=0}
  n=nrow(dat)
  folds <- cut(seq(1,nrow(dat)),breaks=nfolds,labels=FALSE)
  
  #======================================#
  #====== step 1: baseline model ========#
  #======================================#
  
  ### logistic: remember to change the covariate if you change the data generation function
  mod_lr = glm(A~x1+x2, data = dat, family = "binomial")
  preds_lr = predict(mod_lr, new_data = dat[,-c(1:4)], type = "response")
  
  cv_lambda<-NULL
  preds_cv<-dat$A
  coef_cv<-vector(mode = "list", length = nfolds)
  
  for(i in 1:nfolds){
    testIndexes <- which(folds==i,arr.ind=TRUE)
    test <- dat[testIndexes,]
    train <- dat[-testIndexes,]
    mod_cv = fit_hal(X = train[,-c(1:4)], Y = train$A, family = "binomial", yolo = FALSE, return_x_basis = T,max_degree = 2)
    cv_lambda[i] = mod_cv$lambda_star
    preds_cv[testIndexes] = predict(mod_cv, new_data = test[,-c(1:4)])
    coef_cv[[i]]<-mod_cv$coefs
    
  }
  l1cvA<-mean(sapply((coef_cv),function(x) sum(abs(x))))
  
  ### ipw estimator
  lr_ipw = sum(dat$A*dat$Y/preds_lr)/sum(dat$A/preds_lr)
  cv_ipw = sum(dat$A*dat$Y/preds_cv)/sum(dat$A/preds_cv); cv_ipw
  unadj_est = mean(dat$Y[dat$A==1])
  
  ### se estimation for hal
  cv_se = sqrt(sum((dat$A*(dat$Y - cv_ipw)/preds_cv /sum(dat$A/preds_cv) * n)^2)/n^2); cv_se
  #cv_se_half1 = sqrt(sum((dat$A[s]*(dat$Y[s] - cv_ipw)/preds_cv[s] /sum(dat$A[s]/preds_cv[s]) * n)^2)/n^2); cv_se_half1
  lr_se = sqrt(sum((dat$A*(dat$Y - lr_ipw)/preds_lr /sum(dat$A/preds_lr) * n)^2)/n^2)
  
  
  
  Ymod_cv = tryCatch(fit_hal(X = dat[,c("A","x1","x2")], Y = dat$Y, yolo = FALSE, max_degree = 2), error=function(e) e, warning=function(w) w)
  ndat<-dat[,c("A","x1","x2")]
  ndat$A<-1
  Q<-predict(Ymod_cv, new_data = ndat)
  coef_Y<-Ymod_cv$coefs
  l1Y<-sum(abs(coef_Y))
  
  Q0<-dat$Q0
  
  sigma.star<- sqrt(mean((dat$A*dat$Y/preds_lr - (dat$A - Q0)/preds_lr *Q0 - lr_ipw)^2)/n)
  
  
  
  
  
  preds_grid<-dat$A
  D_CAR_under<-score1<-score2<-score3<-temp1<-temp2<-temp5<-score5<-score6<-temp6<-ipw_loop<-temp3<-score3<-temp7<-score7<-
    temp8<-score8<-temp9<-score9<-Ahat.sd<-temp10<-score10<-NULL
  
  preds_under<-matrix(0,ncol=20,nrow=n)
  sum_grid<-l1bound<-NULL
  coef_grid<-basis_func_grid<-vector(mode = "list", length = nfolds)
  for(l in 1:20){
    under_lambda = cv_lambda*rate^l
    for(i in 1:nfolds){
      testIndexes <- which(folds==i,arr.ind=TRUE)
      test <- dat[testIndexes,]
      train <- dat[-testIndexes,]
      
      mod_grid = fit_hal(X = train[,-(1:4)], Y = train$A, family = "binomial", yolo = FALSE, return_x_basis = T, lambda = under_lambda[i], 
                         cv_select = F,max_degree = 2)
      preds_grid[testIndexes] = predict(mod_grid, new_data = test[,-c(1:4)])
      basis_func_grid[[i]] = as.matrix(cbind(1, mod_grid$x_basis)) ### n by p, the row 1 is for intercept term
      coef_grid[[i]]<-mod_grid$coefs
      sum_grid[i]<-sum(abs(mod_grid$coefs))
      
      
    }
    
    A<-dat$A
    A_hat = preds_grid
    A_hat[A_hat > (1-trunc)] = (1-trunc)
    A_hat[A_hat < trunc] = trunc
    
    for(i in 1:nfolds){
      testIndexes <- which(folds==i,arr.ind=TRUE)
      A<-dat$A
      #   temp1[i]<-min( (abs(t(basis_func_grid[[i]]) %*% as.matrix((A[-testIndexes] - preds_grid[-testIndexes])/preds_grid[-testIndexes])) /dim(basis_func_grid[[i]])[1])[abs(coef_grid[[i]])>0])
      
      temp1[i]<-sum(abs(coef_grid[[i]]) * abs(t((basis_func_grid[[i]])) %*% 
                                                as.matrix((A[-testIndexes] - A_hat[-testIndexes]))/dim(basis_func_grid[[i]])[1]))/sum(abs(coef_grid[[i]]))
      
      
      temp2[i]<-sum(abs(coef_cv[[i]]) * abs(t((basis_func_grid[[i]])) %*% 
                                              as.matrix((A[-testIndexes] - A_hat[-testIndexes]))/dim(basis_func_grid[[i]])[1]))/sum(abs(coef_grid[[i]]))
      
      xx<-(  abs(t((basis_func_grid[[i]])) %*% 
                   as.matrix((A[-testIndexes] - A_hat[-testIndexes]))/dim(basis_func_grid[[i]])[1]))
      
      temp3[i]<-quantile((xx[abs(coef_grid[[i]])>0]),prob=c(0.02,0.25,0.5,0.75))[1]
      temp5[i]<-sum(abs(coef_grid[[i]]!=0) * abs(t((basis_func_grid[[i]])) %*% as.matrix((A[-testIndexes] - A_hat[-testIndexes])/
                                                                                           A_hat[-testIndexes]/(1-A_hat[-testIndexes]))/ dim(basis_func_grid[[i]])[1]))/sum(abs(coef_grid[[i]]))
      
      temp6[i]<-sum(abs(coef_grid[[i]]!=0) * abs(t((basis_func_grid[[i]])) %*% as.matrix((A[-testIndexes] - A_hat[-testIndexes])/
                                                                                           A_hat[-testIndexes])/ dim(basis_func_grid[[i]])[1]))/sum(abs(coef_grid[[i]]))
      
      temp7[i]<-(max(abs(xx[abs(coef_grid[[i]])>0])))/sum(abs(coef_grid[[i]]))
      
      xx<-( coef_grid[[i]]* (t((basis_func_grid[[i]])) %*% 
                               as.matrix((A[-testIndexes] - A_hat[-testIndexes]))/dim(basis_func_grid[[i]])[1]))
      temp8[i]<-(quantile(abs(xx[abs(coef_grid[[i]])>0]),prob=c(0.01)))
      
      xx<-(  abs(t((basis_func_grid[[i]])) %*% 
                   as.matrix((A[-testIndexes] - A_hat[-testIndexes]))/dim(basis_func_grid[[i]])[1]))
      temp10[i]<-quantile((xx[abs(coef_grid[[i]])>0]),prob=c(0.0,0.25,0.5,0.95))[1]
      
    }
    
    preds_under[,l]<-preds_grid
    l1bound[l]<-mean(abs(sum_grid))
    D_CAR_under[l]<-abs(t(Q) %*% as.matrix((A - A_hat)/A_hat))/n
    ipw_loop[l] = sum(dat$A * dat$Y/A_hat)/sum(A/A_hat)
    
    score1[l]<-mean(temp1)
    
    score2[l] <- mean(temp2)
    score3[l] <- mean(temp3)
    score5[l] <- mean(temp5)
    score6[l] <- mean(temp6)
    score7[l] <- mean(temp7)
    score8[l] <- mean(temp8,na.rm=T)
    score10[l] <- mean(temp10,na.rm=T)
    
  } 
  
  
  D_CAR_under; which.min(D_CAR_under)

  under_lambda_DCAR1<- cv_lambda*rate^(which.min(D_CAR_under))     
  under_lambda_SCORE<-cv_lambda*rate^which.min(score6)
  under_lambda_1<-cv_lambda*rate^0
  under_lambda_1<-cv_lambda*rate^5
  under_lambda_2<-cv_lambda*rate^10
  under_lambda_3<-cv_lambda*rate^15  
  under_lambda_4<-cv_lambda*rate^20  
  under_lambda_5<-cv_lambda*rate^20  
  under_lambda_6<-cv_lambda*rate^20 

  count.DCAR1<-which.min(D_CAR_under)
  count.score<-which.min(score6)
  count.1<-0
  count.1<-5 
  count.2<-10
  count.3<-15
  count.4<-20 
  count.5<-20 
  count.6<-20 
  
  counts<-c(count.DCAR1,count.score,count.1,count.2,count.3,count.4,count.5,count.6)
  
  #======================================#
  #=== step 2: re-estimation ============#
  #======================================#
  lambda_vec<-rbind(under_lambda_DCAR1,under_lambda_SCORE,under_lambda_1,under_lambda_2,under_lambda_3,under_lambda_4,
                    under_lambda_5,under_lambda_6)
  
  
  
  Q.lm.t<-(lm(Y~A+I(x1^2)+x1*x2,data=dat))
  Q.lm<-predict(Q.lm.t, new_data = ndat)
  
  D_CAR_cv<-abs(t(Q)%*% as.matrix((A - preds_cv)/preds_cv)/n);D_CAR_cv
  D_CAR_lr<- abs(t(Q.lm) %*% as.matrix((A - preds_lr)/preds_lr)/n)
  
  
  preds_under_trun = preds_under ; summary(preds_under)
  preds_under_trun[preds_under_trun > (1-trunc)] = (1-trunc)
  preds_under_trun[preds_under_trun < trunc] = trunc
  
  D_CAR_under<-under_ipw<-sigma.eif<-cp_under<-NULL
  
  for(i in 1:nrow(lambda_vec)) {
    temp.sel<-counts[i]
    D_CAR_under[i]<-abs(t(Q) %*% as.matrix((A - preds_under_trun[,temp.sel])/preds_under_trun[,temp.sel])/n)
    under_ipw[i]<-ipw_loop[temp.sel]
    sigma.eif[i]<-sqrt(mean((dat$A*dat$Y/preds_under_trun[,temp.sel] - (dat$A - preds_under_trun[,temp.sel])/preds_under_trun[,temp.sel] *Q - 
                               under_ipw[i])^2)/n)
    cp_under[i]<-mean(under_ipw[i] - 1.96* sigma.eif[1] <= true & under_ipw[i] + 1.96* sigma.eif[1] >= true)
  }
  
  
  sigma.star_DCAR1<- sigma.eif[1]
 
  cp_cv = mean(cv_ipw - 1.96*sigma.star_DCAR1 <= true & cv_ipw + 1.96*sigma.star_DCAR1 >= true)
  cp_lr = mean(lr_ipw - 1.96*lr_se <= true & lr_ipw + 1.96*lr_se >= true)
  
  bias<-( under_ipw-true)
  rnbias<-sqrt(n) * ( bias)
  names(rnbias)<-c("rnbias_DCAR1","rnbias_SCORE","rnbias_1","rnbias_2","rnbias_3","rnbias_4","rnbias_5","rnbias_6")
  rnmse<-sqrt(n)*(bias^2+sigma.eif^2)
  names(rnmse)<-c("rnmse_DCAR1","rnmse_SCORE","rnmse_1","rnmse_2","rnmse_3","rnmse_4","rnmse_5","rnmse_6")
  
  names(cp_under)<-c("cp_DCAR1","cp_SCORE","cp_1","cp_2","cp_3","cp_4","cp_5","cp_6")
  
#  result = c(unadj_est=unadj_est,lr_ipw=lr_ipw,lr_se=lr_se,cp_lr=cp_lr, cv_ipw=cv_ipw,cp_cv=cp_cv, cv_se=cv_se,  
#             under_ipw=under_ipw,cp_under=cp_under,
#             cv_lambda=cv_lambda, under_lambda = apply(lambda_vec,1,mean), D_CAR_under=D_CAR_under,
#             D_CAR_cv=D_CAR_cv,D_CAR_lr=D_CAR_lr, 
#             sigma.star=sigma.eif, 
#             l1bound=l1bound,counts=counts,rnbias,rnmse)
  result = c(unadj_est=unadj_est,lr_ipw=lr_ipw,cp_lr=cp_lr, cv_ipw=cv_ipw,cp_cv=cp_cv,  
             under_ipw=under_ipw,cp_under=cp_under,
             cv_lambda=cv_lambda, under_lambda = apply(lambda_vec,1,mean), D_CAR_under=D_CAR_under,
             D_CAR_cv=D_CAR_cv,D_CAR_lr=D_CAR_lr, 
             sigma.star=sigma.eif, 
             rnbias,rnmse)
  return(result)
}



nsim=48*1
sim.1 <- vector("list", length = nsim)
rate=0.95 # rate of change of lambda on a grid
nfolds=10 # number of cross-fitting folds
trunc=0.00 # truncation level, e.g., 0, 0.02,0.05.

set.seed(123)
n = 200
ss=5
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

results.list <- foreach(k =1:length(sim.1), .packages = c("hal9001","MASS","sandwich")) %dopar%  
  tryCatch(onestage_under_hal(dat = sim.1[[k]], rate = rate,nfolds=nfolds,ss=ss,trunc=trunc), error = function(e) {print(e); NA})

stopCluster(cl)
round(apply(do.call(rbind,results.list),2,function(x) mean(x,na.rm=TRUE)),3)

