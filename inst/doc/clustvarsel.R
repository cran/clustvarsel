## ----setup, include=FALSE------------------------------------------------
library(knitr)
opts_chunk$set(fig.align = "center", out.width = "100%",
               fig.width = 5, fig.height = 5,
               dev.args = list(pointsize=10),
               par = TRUE, # needed for setting hook 
               collapse = TRUE, # collapse input & ouput code in chunks
               cache = FALSE,
               warning = FALSE, message = FALSE)
knit_hooks$set(par = function(before, options, envir)
  { if(before && options$fig.show != "none") 
       par(mar=c(4.1,4.1,1.1,1.1), mgp=c(3,1,0), tcl=-0.5)
})

## ------------------------------------------------------------------------
library(clustvarsel)

## ---- echo=-1------------------------------------------------------------
set.seed(20170224)
n = 200      # sample size
pro = 0.5    # mixing proportion
mu1 = c(0,0) # mean vector for the first cluster
mu2 = c(3,3) # mean vector for the second cluster
sigma1 = matrix(c(1,0.5,0.5,1),2,2)       # covar matrix for the first cluster
sigma2 = matrix(c(1.5,-0.7,-0.7,1.5),2,2) # covar matrix for the second cluster
X = matrix(0, n, 5, dimnames = list(NULL, paste0("X", 1:5)))
set.seed(1234) # for replication
u = runif(n)
Class = ifelse(u < pro, 1, 2)
X[u < pro, 1:2]  = MASS::mvrnorm(sum(u < pro), mu = mu1, Sigma = sigma1)
X[u >= pro, 1:2] = MASS::mvrnorm(sum(u >= pro), mu = mu2, Sigma = sigma2)
X[, 3] = X[, 1] + rnorm(n)
X[, 4] = rnorm(n, mean = 1.5, sd = 2)
X[, 5] = rnorm(n, mean = 2, sd = 1)
clPairs(X, Class)

## ------------------------------------------------------------------------
out = clustvarsel(X, verbose = TRUE)
out

out = clustvarsel(X, direction = "backward", verbose = TRUE)
out

out = clustvarsel(X, search = "headlong", verbose = TRUE)
out

## ---- echo=-1------------------------------------------------------------
set.seed(20170225)
n = 200
p = 10
mu = rep(0,p)
sigma1 = matrix(c(1,0.5,0.5,1),2,2)
sigma2 = matrix(c(1.5,-0.7,-0.7,1.5),2,2)
sigma = Matrix::bdiag(sigma1, sigma2, diag(6))
set.seed(12345)
X = MASS::mvrnorm(n, mu, sigma)
colnames(X) = paste0("X", 1:p)
pairs(X, gap = 0)

## ------------------------------------------------------------------------
mod = Mclust(X)
summary(mod$BIC)
summary(mod)
plot(mod, what = "classification")

## ------------------------------------------------------------------------
(out1 = clustvarsel(X))
mod1 = Mclust(X[,out1$subset])
summary(mod1)
plot(mod1, what = "classification")

## ------------------------------------------------------------------------
(out2 = clustvarsel(X, direction = "backward"))
mod2 = Mclust(X[,out2$subset])
summary(mod2)
plot(mod2, what = "classification")

## ------------------------------------------------------------------------
sessionInfo()

