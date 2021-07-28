if (grepl("savio2", Sys.info()["nodename"])) {
  .libPaths("/global/scratch/nhejazi/R")
  Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")
}

# set CRAN mirror
options(repos = structure(c(CRAN = "https://cran.rstudio.com/")))

# lie to pkgbuild, as per Jeremy
pkgbuild:::cache_set("has_compiler", TRUE)

# from CRAN
install.packages(c("R.utils", "Rcpp", "here", "remotes", "devtools",
                   "data.table", "tidyverse", "origami",
                   "foreach", "doRNG", "future", "future.apply", "doFuture"),
                 lib = "/global/scratch/nhejazi/R")

# use remotes to install from GitHub
remotes::install_github("tlverse/hal9001@master",
                        lib = "/global/scratch/nhejazi/R")

# install helper package locally
remotes::install_local(here::here("..", "uhalipw"),
                       force = TRUE, upgrade = FALSE,
                       lib = "/global/scratch/nhejazi/R")

# update all packages
update.packages(ask = FALSE, lib.loc = "/global/scratch/nhejazi/R")
