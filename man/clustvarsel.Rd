\name{clustvarsel}
\alias{clustvarsel}
\alias{print.clustvarsel}
  
\title{Variable Selection for Gaussian Model-Based Clustering}

\description{A function which implements variable selection methodology for model-based clustering which allows to find the (locally) optimal subset of variables in a dataset that have group/cluster information. 
A greedy or headlong search can be used, either in a forward-backward or backward-forward direction, with or without sub-sampling at the hierarchical clustering stage for starting  \code{\link[mclust]{mclust}} models. By default the algorithm uses a sequential search, but parallelisation is also available.}

\usage{
clustvarsel(data, 
            G = 1:9, 
            search = c("greedy", "headlong"),
            direction = c("forward", "backward"),
            emModels1 = c("E", "V"), 
            emModels2 = mclust.options("emModelNames"),
            samp = FALSE, 
            sampsize = round(nrow(data)/2), 
            hcModel = "VVV", 
            allow.EEE = TRUE, 
            forcetwo = TRUE,
            BIC.diff = 0, 
            BIC.upper = 0, 
            BIC.lower = -10, 
            itermax = 100, 
            parallel = FALSE,
            fit = TRUE,
            verbose = interactive())
}
                        
\arguments{

\item{data}{A numeric matrix or data frame where rows correspond to observations and columns correspond to variables. Categorical variables are not allowed.}

\item{G}{An integer vector specifying the numbers of mixture components (clusters) for which the BIC is to be calculated. The default is \code{G = 1:9}.}

\item{search}{A character vector indicating whether a \code{"greedy"} or, potentially quicker but less optimal, \code{"headlong"} algorithm is to be used in the search for clustering variables.}

\item{direction}{A character vector indicating the type of search: \code{"forward"} starts from the empty model and at each step of the algorithm adds/removes a variable until the stopping criterion is satisfied; \code{"backward"} starts from the model with all the available variables and at each step of the algorithm removes/adds a variable until the stopping criterion is satisfied.}

\item{emModels1}{A vector of character strings indicating the models to be fitted in the EM phase of univariate clustering. Possible models are \code{"E"} and \code{"V"}, described in \code{\link{mclustModelNames}}.}

\item{emModels2}{A vector of character strings indicating the models to be fitted in the EM phase of multivariate clustering. Possible models are those described in \code{\link{mclustModelNames}}.}

\item{samp}{A logical value indicating whether or not a subset of observations is to be used in the hierarchical clustering phase used to get starting values for the EM algorithm.}

\item{sampsize}{The number of observations to be used in the hierarchical clustering subset. By default, a random sample of approximately half of the sample size is used.}

\item{hcModel}{A character string specifying the model to be used in hierarchical clustering for choosing the starting values used by the EM algorithm. By default, the \code{"VVV"} covariance structure is used (see \code{\link{hc}}).}

\item{allow.EEE}{A logical value indicating whether a new clustering will be run with equal within-cluster covariance for hierarchical clustering to get starting values, if the clusterings with variable within-cluster covariance for hierarchical clustering do not produce any viable BIC values.}

\item{forcetwo}{A logical value indicating whether at least two variables will be forced to be selected initially, regardless of whether BIC evidence suggests bivariate clustering or not.}

\item{BIC.diff}{A numerical value indicating the minimum BIC difference between clustering and no clustering used to accept the inclusion of a variable in the set of clustering variables in a forward step of the greedy search algorithm. Furthermore, minus \code{BIC.diff} is used to accept the exclusion of a selected variable from the set of clustering variable in a backward step of the greedy search algorithm. Default is 0.}

\item{BIC.upper}{A numerical value indicating the minimum BIC difference between clustering and no clustering used to select a clustering variable in the headlong search. Default is 0.}

\item{BIC.lower}{A numerical value indicating the level of BIC difference between clustering and no clustering below which a variable will be removed from consideration in the headlong algorithm. Default is -10.}

\item{itermax}{An integer value giving the maximum number of iterations (of addition and removal steps) the selected algorithm is allowed to run for.}

\item{parallel}{An optional argument which allows to specify if the selected \code{"greedy"} algorithm should be run sequentially or in parallel. 

For a single machine with multiple cores, possible values are: 
\itemize{
  \item a logical value specifying if parallel computing should be used (\code{TRUE}) or not (\code{FALSE}, default) for running the algorithm;
  \item  a numerical value which gives the number of cores to employ. By default, this is obtained from the function \code{\link[parallel]{detectCores}}; 
  \item a character string specifying the type of parallelisation to use. This depends on system OS: on Windows OS only \code{"snow"} type functionality is available, while on Unix/Linux/Mac OSX both \code{"snow"} and \code{"multicore"} (default) functionalities are available.
}
In all the cases described above, at the end of the search the cluster is automatically stopped by shutting down the workers.

If a cluster of multiple machines is available, the algorithm can be run in parallel using all, or a subset of, the cores available to the machines belonging to the cluster. However, this option requires more work from the user, who needs to set up and register a parallel back end. 
In this case, the cluster must be explicitly stopped with \code{\link[parallel]{stopCluster}}.
}

\item{fit}{A logical specifying if the model with the selected "best" subset should be fitted at the end of the procedure. By default is set to \code{TRUE}.}

\item{verbose}{A logical indicating if info must be provided at each step of the algorithm. By default is set to \code{TRUE} during interactive sessions, and \code{FALSE} otherwise.}
}

\details{
This function implements variable selection methodology for model-based clustering. 
The selection methods available are \code{"greedy"} search and \code{"headlong"} search (see argument \code{search}).
The \emph{greedy search} at each step either checks all the variables not currently included in the set of clustering variables singly for inclusion into the set, or checks all variables in the set of clustering variables singly for exclusion. 
The \emph{headlong search} only checks until a variable is included or excluded (i.e., does not necessarily check all possible variables for inclusion/exclusion at each step) and any variable with evidence of clustering below a certain level at any stage is removed from consideration for the remainder of the algorithm. 
Greedy search can be performed \emph{forward-backward} starting from the empty model, or \emph{backward-forward} starting from the model with all the variables included (see argument \code{direction}). Currently, headlong search can only run \emph{forward-backward}.
% Both greedy and headlong searches can be performed forward-backward starting from the empty model, or backward-forward from the model with all the variables included (see argument \code{direction}). 

The criterion to assess each variable's evidence for being useful to the clustering is given by the difference between the BIC for the clustering model using the set of clustering variables already included and the variable being checked, and the sum of BICs for the model with clustering using the set of already selected clustering variables without the variable being checked and the model for the variable being checked being conditionally independent of the clustering given the other clustering variables. The latter is modeled as a regression of the variable being checked on the other clustering variables. A subset selection procedure is also performed on the regression step.

Clustering models are fitted using \code{\link[mclust]{Mclust}}, where each model is allowed to vary over the specified number of components \code{G} and any of the different covariance parameterisations set in \code{emModels1} and \code{emModels2}.

The default value for \code{forcetwo} is \code{TRUE} because often in practice there will be little evidence of clustering on the univariate or bivariate level although there is multivariate clustering present and these variables are used as starting points to attempt to find this clustering, if necessary being removed later in the algorithm.

The default value for \code{allow.EEE} is \code{TRUE} but if necessary to speed up the algorithm it can be set to \code{FALSE}. Other speeding-up restrictions include reducing the \code{emModels1} (to \code{"E"}, say) and the \code{emModels2} to a smaller set of covariance parameterisations. 
Reducing the maximum possible number of clusters present in the data will also increase the speed of the algorithm.
Another time-saving device is the \code{samp} option which uses the same algorithm but uses only a subset of the observations in the expensive hierarchical phase for initialisation of EM algorithm in \code{\link[mclust]{Mclust}}. 
Finally, the headlong search may be quicker than the greedy search option in data sets with large numbers of variables (depending on the values of the upper and lower bounds chosen for the BIC difference).

The defaults for the \code{eps}, \code{tol} and \code{itmax} options for the Mclust
steps run in the algorithm can be changed by using the \code{\link[mclust]{mclust.options}} function.
}

\value{An object of class \code{'clustvarsel'} with the following components:

\item{variables}{The name of input variables.}

\item{subset}{A vector of values specifying for the selected variables the columns position in the input data.frame or matrix.}

\item{steps.info}{A matrix with a row for each step of the algorithm providing:
\tabular{ll}{
- \tab the name of the variable proposed; \cr
- \tab the BIC of the clustering variables' model at the end of the step; \cr
- \tab the BIC difference between clustering and not clustering for the variable; \cr
- \tab the type of step (Add/Remove); \cr
- \tab the decision for the proposed step (Accepted/Rejected).
} 
}

\item{search}{A string specifying the type of search employed.}

\item{direction}{A string specifying the direction of search employed.}

\item{model}{If input argument \code{fit = TRUE}, a \code{'Mclust'} object containing the final model fitted on the "best" subset selected.}
}


\references{
Raftery, A. E. and Dean, N. (2006) Variable Selection for Model-Based Clustering. \emph{Journal of the American Statistical Association}, 101(473), 168-178.
\url{http://www.stat.washington.edu/www/research/reports/2004/tr452.pdf}

Badsberg, J. H. (1992) Model search in contingency tables by CoCo. In Dodge, Y. and Whittaker, J. (Eds.), \emph{Computational Statistics}, Volume 1, pp. 251-256

Maugis, C., Celeux, G., Martin-Magniette M. (2009) Variable Selection for Clustering With Gaussian Mixture Models. \emph{Biometrics}, 65(3), 701-709.

Scrucca, L. and Raftery, A. E. (2018) clustvarsel: A Package Implementing Variable Selection for Gaussian Model-based Clustering in R. \emph{Journal of Statistical Software}, 84(1), pp. 1-28.
}

\author{N. Dean, A. E. Raftery, and L. Scrucca.}

\seealso{\code{\link[mclust]{Mclust}}}

\examples{
# Simulate data with 2 clusters in the first two variables and no 
# clustering in the rest. Clusters have mixing proportion pro, means 
# mu1 and mu2 and variances sigma1 and sigma2
require(MASS)
n <- 200
pro <- 0.5
mu1 <- c(0,0)
mu2 <- c(3,3)
sigma1 <- matrix(c(1,0.5,0.5,1),2,2,byrow=TRUE)
sigma2 <- matrix(c(1.5,-0.7,-0.7,1.5),2,2,byrow=TRUE)
X <- matrix(0, n, 5)
colnames(X) <- paste("X", 1:ncol(X), sep ="")
# generate the clustering variables
u <- runif(n)
Class <- ifelse(u < pro, 1, 2)
X[u < pro, 1:2]  <- mvrnorm(sum(u < pro), mu = mu1, Sigma = sigma1)
X[u >= pro, 1:2] <- mvrnorm(sum(u >= pro), mu = mu2, Sigma = sigma2)
# generate the non-clustering variables
X[,3] <- X[,1] + rnorm(n)
X[,4] <- rnorm(n, mean = 1.5, sd = 2)
X[,5] <- rnorm(n, mean = 2, sd = 1)
# plot the data
clPairs(X, Class, gap = 0)

# sequential forward greedy search (default)
out <- clustvarsel(X, G = 1:5)
out
summary(out$model)
table(Class, out$mod$classification)

\dontrun{
# sequential backward greedy search
clustvarsel(X, G = 1:5, direction = "backward")

# sequential backward greedy search with subsampling at hierarchical 
# intialisation stage
clustvarsel(X, G = 1:5, direction = "backward", 
            samp = TRUE, sampsize = 50)

# parallel backward greedy search 
clustvarsel(X, G = 1:5, direction = "backward", parallel = TRUE)

# headlong search 
clustvarsel(X, G = 1:5, search = "headlong")
}
}

\keyword{cluster}
