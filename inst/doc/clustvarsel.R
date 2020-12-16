## ----setup, include=FALSE-----------------------------------------------------
library(knitr)
opts_chunk$set(fig.align = "center", 
               out.width = "90%",
               fig.width = 6, 
               fig.height = 5.5,
               dev.args = list(pointsize=10),
               par = TRUE, # needed for setting hook 
               collapse = TRUE, # collapse input & ouput code in chunks
               cache = FALSE,
               warning = FALSE, message = FALSE)
knit_hooks$set(par = function(before, options, envir)
  { if(before && options$fig.show != "none") 
       par(mar=c(4.1,4.1,1.1,1.1), mgp=c(3,1,0), tcl=-0.5)
})

## -----------------------------------------------------------------------------
library(mclust)
library(clustvarsel)

## -----------------------------------------------------------------------------
n <- 200      # sample size
pro <- 0.5    # mixing proportion
mu1 <- c(0,0) # mean vector for the first cluster
mu2 <- c(3,3) # mean vector for the second cluster
sigma1 <- matrix(c(1,0.5,0.5,1),2,2)       # covar matrix for the first cluster
sigma2 <- matrix(c(1.5,-0.7,-0.7,1.5),2,2) # covar matrix for the second cluster
X <- matrix(0, n, 5, dimnames = list(NULL, paste0("X", 1:5)))
set.seed(1234) # for replication
u <- runif(n)
Class <- ifelse(u < pro, 1, 2)
X[u < pro, 1:2]  <- MASS::mvrnorm(sum(u < pro), mu = mu1, Sigma = sigma1)
X[u >= pro, 1:2] <- MASS::mvrnorm(sum(u >= pro), mu = mu2, Sigma = sigma2)
X[, 3] <- X[, 1] + rnorm(n)
X[, 4] <- rnorm(n, mean = 1.5, sd = 2)
X[, 5] <- rnorm(n, mean = 2, sd = 1)
clPairs(X, Class)

## ----eval=FALSE---------------------------------------------------------------
#  out <- clustvarsel(X)
#  ## iter 1
#  ## + adding step
#  ##   Var  BICdiff Step Decision
#  ## 1  X2 11.17931  Add Accepted
#  ## iter 2
#  ## + adding step
#  ##   Var BICdiff Step Decision
#  ## 2  X1 85.1953  Add Accepted
#  ## iter 3
#  ## + adding step
#  ## - removing step
#  ##   Var   BICdiff   Step Decision
#  ## 3  X3 -14.91130    Add Rejected
#  ## 4  X1  85.19104 Remove Rejected
#  ## final iter
#  ## * fitting model on selected subset
#  
#  out
#  ## ------------------------------------------------------
#  ## Variable selection for Gaussian model-based clustering
#  ## Stepwise (forward/backward) greedy search
#  ## ------------------------------------------------------
#  ##
#  ##  Variable proposed Type of step   BICclust Model G   BICdiff Decision
#  ##                 X2          Add  -822.6398     E 2  11.17931 Accepted
#  ##                 X1          Add -1482.8408   VEV 2  85.19530 Accepted
#  ##                 X3          Add -2047.7064   EEV 2 -14.91130 Rejected
#  ##                 X1       Remove  -822.6355     E 2  85.19104 Rejected
#  ##
#  ## Selected subset: X2, X1
#  
#  summary(out$model)
#  ## ----------------------------------------------------
#  ## Gaussian finite mixture model fitted by EM algorithm
#  ## ----------------------------------------------------
#  ##
#  ## Mclust VEV (ellipsoidal, equal shape) model with 2 components:
#  ##
#  ##  log.likelihood   n df       BIC       ICL
#  ##       -714.9288 200 10 -1482.841 -1493.898
#  ##
#  ## Clustering table:
#  ##   1   2
#  ##  99 101

## ----eval=FALSE---------------------------------------------------------------
#  out <- clustvarsel(X, direction = "backward", fit = FALSE)
#  ## iter 1
#  ## - removing step
#  ##   Var   BICdiff   Step Decision
#  ## 1  X3 -38.09195 Remove Accepted
#  ## iter 2
#  ## - removing step
#  ##   Var   BICdiff   Step Decision
#  ## 2  X4 -23.50212 Remove Accepted
#  ## iter 3
#  ## - removing step
#  ## + adding step
#  ##   Var   BICdiff   Step Decision
#  ## 3  X5 -16.25831 Remove Accepted
#  ## 4  X3 -14.91130    Add Rejected
#  ## iter 4
#  ## - removing step
#  ## + adding step
#  ##   Var   BICdiff   Step Decision
#  ## 5  X1  95.55735 Remove Rejected
#  ## 6  X3 -14.91130    Add Rejected
#  
#  out
#  ## ------------------------------------------------------
#  ## Variable selection for Gaussian model-based clustering
#  ## Stepwise (backward/forward) greedy search
#  ## ------------------------------------------------------
#  ##
#  ##  Variable proposed Type of step   BICclust Model G   BICdiff Decision
#  ##                 X3       Remove -2925.4851   EVE 2 -38.09195 Accepted
#  ##                 X4       Remove -2067.0668   EVE 2 -23.50212 Accepted
#  ##                 X5       Remove -1482.8408   VEV 2 -16.25831 Accepted
#  ##                 X3          Add -2047.7064   EEV 2 -14.91130 Rejected
#  ##                 X1       Remove  -833.0019     V 2  95.55735 Rejected
#  ##                 X3          Add -2047.7064   EEV 2 -14.91130 Rejected
#  ##
#  ## Selected subset: X1, X2

## ----eval=FALSE---------------------------------------------------------------
#  out <- clustvarsel(X, search = "headlong")
#  ## iter 1
#  ## + adding step
#  ##   Var  BICdiff Step Decision
#  ## 1  X2 11.17931  Add Accepted
#  ## iter 2
#  ## + adding step
#  ##   Var BICdiff Step Decision
#  ## 2  X1 85.1953  Add Accepted
#  ## iter 3
#  ## + adding step
#  ## - removing step
#  ##   Var  BICdiff   Step Decision
#  ## 3  X3 -14.9113    Add Rejected
#  ## 4  X1  85.1953 Remove Rejected
#  ## final iter
#  ## * fitting model on selected subset
#  
#  out
#  ## ------------------------------------------------------
#  ## Variable selection for Gaussian model-based clustering
#  ## Headlong (forward/backward) search
#  ## ------------------------------------------------------
#  ##
#  ##  Variable proposed Type of step   BICclust Model G   BICdiff Decision
#  ##                 X2          Add  -822.6398     E 2  11.17931 Accepted
#  ##                 X1          Add -1482.8408   VEV 2  85.19530 Accepted
#  ##                 X3          Add -2047.7064   EEV 2 -14.91130 Rejected
#  ##                 X1       Remove  -822.6398     E 2  85.19530 Rejected
#  ##
#  ## Selected subset: X2, X1

## -----------------------------------------------------------------------------
n <- 200
p <- 10
mu <- rep(0,p)
sigma1 <- matrix(c(1,0.5,0.5,1),2,2)
sigma2 <- matrix(c(1.5,-0.7,-0.7,1.5),2,2)
sigma <- Matrix::bdiag(sigma1, sigma2, diag(6))
set.seed(12345)
X <- MASS::mvrnorm(n, mu, sigma)
colnames(X) <- paste0("X", 1:p)
clPairs(X)

## ----eval=FALSE---------------------------------------------------------------
#  mod <- Mclust(X)
#  summary(mod$BIC)
#  ## Best BIC values:
#  ##              EII,2       VII,2       EII,3
#  ## BIC      -5899.073 -5901.88582 -5922.53452
#  ## BIC diff     0.000    -2.81261   -23.46131
#  summary(mod)
#  ## ----------------------------------------------------
#  ## Gaussian finite mixture model fitted by EM algorithm
#  ## ----------------------------------------------------
#  ##
#  ## Mclust EII (spherical, equal volume) model with 2 components:
#  ##
#  ##  log.likelihood   n df       BIC       ICL
#  ##       -2891.255 200 22 -5899.073 -5953.953
#  ##
#  ## Clustering table:
#  ##   1   2
#  ##  90 110

## ----eval=FALSE---------------------------------------------------------------
#  (out1 <- clustvarsel(X, verbose = FALSE))
#  ## ------------------------------------------------------
#  ## Variable selection for Gaussian model-based clustering
#  ## Stepwise (forward/backward) greedy search
#  ## ------------------------------------------------------
#  ##
#  ##  Variable proposed Type of step   BICclust Model G    BICdiff Decision
#  ##                 X3          Add  -666.9232     E 2 -4.8472740 Accepted
#  ##                 X7          Add -1242.8998   EII 2  1.5678544 Accepted
#  ##                X10          Add -1814.8848   EII 2  1.5172709 Accepted
#  ##                X10       Remove -1242.8998   EII 2  1.5172709 Rejected
#  ##                 X6          Add -2390.3224   EII 2  0.6539224 Accepted
#  ##                 X6       Remove -1814.8848   EII 2  0.6539224 Rejected
#  ##                 X2          Add -2942.2744   EII 2  0.1214452 Accepted
#  ##                 X2       Remove -2390.3224   EII 2  0.1214452 Rejected
#  ##                 X8          Add -3495.9758   EII 2  0.1104619 Accepted
#  ##                 X8       Remove -2942.2744   EII 2  0.1104619 Rejected
#  ##                 X9          Add -4081.9212   EII 2 -2.0296117 Rejected
#  ##                 X8       Remove -2942.2744   EII 2  0.1104619 Rejected
#  ##
#  ## Selected subset: X3, X7, X10, X6, X2, X8
#  
#  summary(out1$model)
#  # or
#  # mod1 <- Mclust(X[,out1$subset])
#  # summary(mod1)
#  ## ----------------------------------------------------
#  ## Gaussian finite mixture model fitted by EM algorithm
#  ## ----------------------------------------------------
#  ##
#  ## Mclust XII (spherical multivariate normal) model with 1 component:
#  ##
#  ##  log.likelihood   n df       BIC       ICL
#  ##       -1727.188 200  7 -3491.465 -3491.465
#  ##
#  ## Clustering table:
#  ##   1
#  ## 200

## ----eval=FALSE---------------------------------------------------------------
#  (out2 <- clustvarsel(X, direction = "backward", verbose = FALSE))
#  ## ------------------------------------------------------
#  ## Variable selection for Gaussian model-based clustering
#  ## Stepwise (backward/forward) greedy search
#  ## ------------------------------------------------------
#  ##
#  ##  Variable proposed Type of step  BICclust Model G     BICdiff Decision
#  ##                 X2       Remove -5346.204   EII 2 -48.0289667 Accepted
#  ##                 X4       Remove -4696.525   EII 2  -7.1909995 Accepted
#  ##                 X3       Remove -4032.114   EII 2  -2.3352781 Accepted
#  ##                 X3          Add -4696.525   EII 2  -2.3352781 Rejected
#  ##                 X7       Remove -3453.786   EII 2  -1.1967844 Accepted
#  ##                 X7          Add -4032.114   EII 2  -1.1967844 Rejected
#  ##                 X8       Remove -2899.545   EII 2  -0.4288443 Accepted
#  ##                 X8          Add -3453.786   EII 2  -0.4288443 Rejected
#  ##                 X5       Remove -2308.785   EII 2   1.0712786 Rejected
#  ##                 X8          Add -3453.786   EII 2  -0.4288443 Rejected
#  ##
#  ## Selected subset: X1, X5, X6, X9, X10
#  
#  summary(out2$model)
#  ## ----------------------------------------------------
#  ## Gaussian finite mixture model fitted by EM algorithm
#  ## ----------------------------------------------------
#  ##
#  ## Mclust XII (spherical multivariate normal) model with 1 component:
#  ##
#  ##  log.likelihood   n df       BIC       ICL
#  ##       -1424.072 200  6 -2879.934 -2879.934
#  ##
#  ## Clustering table:
#  ##   1
#  ## 200

## -----------------------------------------------------------------------------
sessionInfo()

