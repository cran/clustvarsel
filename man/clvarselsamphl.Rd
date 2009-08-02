\name{clvarselsamphl}
\alias{clvarselsamphl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Headlong Search Variable Selection for Model-Based Clustering with hierarchical clustering sub-sampling}
\description{A function which uses a headlong search, with sub-sampling
  at the hierarchical clustering stage of Mclust, to find the (locally)
  optimal subset of variables in a dataset that have group/cluster
  information. This function is called by the clustvarsel function when the option `samp' is set to TRUE and `search' is set to ``headlong''.
}
\usage{
clvarselsamphl(X, G, emModels1 = c("E","V"), emModels2 = c("EII","VII","EEI",
               "VEI","EVI","VVI","EEE","EEV","VEV","VVV"), sampsize=2000,
               allow.EEE=TRUE, forcetwo=TRUE, upper=0, lower=-10, itermax=100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{A matrix of data with rows corresponding to observations and
    columns (at least 2) corresponding to variables. Categorical
    variables are not permitted.}
  \item{G}{A scalar specifying the maximum number of clusters believed to be present in X.}
  \item{emModels1}{A vector of character strings indicating the models to be fitted in the EM phase of univariate clustering. Possible models:

``E'' for spherical, equal variance

``V'' for spherical, variable variance

The default is all of the above.}
  \item{emModels2}{A vector of character strings indicating the models to be fitted in the EM phase of multivariate clustering. Possible models:

  ``EII'': spherical, equal volume 

           ``VII'': spherical, unequal volume 

           ``EEI'': diagonal, equal volume, equal shape 

           ``VEI'': diagonal, varying volume, equal shape 

           ``EVI'': diagonal, equal volume, varying shape 

           ``VVI'': diagonal, varying volume, varying shape 

           ``EEE'': ellipsoidal, equal volume, shape, and orientation 

           ``EEV'': ellipsoidal, equal volume and equal shape

           ``VEV'': ellipsoidal, equal shape 

           ``VVV'': ellipsoidal, varying volume, shape, and orientation 


           The default is all of the above.}
  \item{sampsize}{The number of observations to be used in the hierarchical clustering subset.}
  \item{allow.EEE}{A logical value indicating whether a new clustering will be run with equal variance hierarchical clustering starting values if the clusterings with variable variance hierarchical clustering starting values do not produce any viable BIC values.}
  \item{forcetwo}{A logical value indicating whether at least two variables will be forced to be selected initially (regardless of whether BIC evidence suggests bivariate clustering or not).}
  \item{upper}{A scalar value indicating the minimum BIC difference
    between clustering and no clustering used to select a clustering variable. Default is 0.}
  \item{lower}{A scalar value indicating the level of BIC difference between clustering and no clustering below which a variable will be removed from consideration in the algorithm. Default is -10.}
  \item{itermax}{A scalar value giving the maximum number of iterations (of addition and removal steps) the algorithm is allowed to run for.}
}
\details{
This function is called by `clustvarsel' when the option `samp' is set
to TRUE and `search' is set to ``headlong''.

The default value for `forcetwo' is TRUE because often in practice there will be little evidence of clustering on the univariate or bivariate level although there is multivariate clustering present and these variables are used as starting points to attempt to find this clustering, if necessary being removed later in the algorithm.

The default value for `allow.EEE' is TRUE but if necessary to speed up
the algorithm it can be set to FALSE. Other speeding-up restrictions
include reducing the `emModels1' (to ``E'', say) and the `emModels2' to
a smaller set of covariance parameterizations. Reducing the maximum
possible number of clusters present in the data will also increase the
speed of the algorithm. The headlong search may be quicker than the
greedy search in larger data sets (depending on the values of the upper
and lower bounds chosen for the BIC difference). Lower values of
`upper' and higher values of `lower' will possibly speed up the search
(although they may make the solution found less optimal).

The defaults for the `eps', `tol' and `itmax' options for the Mclust
steps run in the algorithm can be changed by setting the variables
.Mclust\$eps, .Mclust\$tol and .Mclust\$itmax respectively to new values.

}
\value{A list giving:
  \item{sel.var}{The matrix of selected variables.}
  \item{steps.info}{A matrix with a row for each step of the algorithm giving:

 the name of the best variable proposed, 

 the BIC of the clustering variables' model at the end of the step,

 the BIC difference between clustering and not clustering for the variable,

 the type of step (addition/removal),

 the decision for the variable.}


}
\references{A. E. Raftery and N. Dean (2006). Variable Selection for Model-Based Clustering, Journal of the American Statistical Association, Volume 101, no. 473, pp. 168-178
\url{http://www.stat.washington.edu/www/research/reports/2004/tr452.pdf}

J. H. Badsberg (1992). Model search in contingency tables by CoCo. In Y. Dodge and J. Whittaker (Eds.), Computational Statistics, Volume 1, pp. 251-256
 }
\author{N. Dean and A. E. Raftery}
\seealso{\code{\link{clustvarsel}}, \code{\link{clvarselnosamphl}},
  \code{\link{clvarselsampgr}}, \code{\link{clvarselnosampgr}}, \code{\link[mclust]{Mclust}}}
\examples{
#Create 3-d data with 2 clusters in the first two variables and no
#clustering in the rest
X<-matrix(0,200,3)
colnames(X)<-1:3
#clusters have mixing proportion pro, means mu1 and mu2 and variances
#sigma1 and sigma2
pro<-0.5
mu1<-c(0,0)
mu2<-c(3,3)
sigma1<-matrix(c(1,0.5,0.5,1),2,2,byrow=TRUE)
sigma2<-matrix(c(1.5,-0.7,-0.7,1.5),2,2,byrow=TRUE)
u<-runif(200)
library(MASS)
for(i in 1:200)
{
ifelse(u[i]<pro,X[i,1:2]<-mvrnorm(1,mu1,sigma1),X[i,1:2]<-mvrnorm(1,mu2,sigma2))
X[i,3]<-rnorm(1,1.5,2)
}
#look at the data
#pairs(X)
#Find the clustering variables
m<-clvarselsamphl(X,G=3,sampsize=100)
#Look at the names of the variables selected
colnames(m$sel.var)
m$steps.info
#look at the clustering produced by the variables selected
result<-Mclust(m$sel.var,1:3)
result
}
\keyword{cluster}

