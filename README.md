# clustvarsel

[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/clustvarsel)](https://cran.r-project.org/package=clustvarsel)
[![CRAN\_MonthlyDownloads](http://cranlogs.r-pkg.org/badges/clustvarsel)](https://cran.r-project.org/package=clustvarsel)

An [R](https://www.r-project.org/) package implementing *Variable Selection for Gaussian Model-Based Clustering*.

Variable selection for Gaussian model-based clustering as implemented in the **mclust** package. The methodology allows to find the (locally) optimal subset of variables in a data set that have group/cluster information. A greedy or headlong search can be used, either in a forward-backward or backward-forward direction, with or without sub-sampling at the hierarchical clustering stage for starting mclust models. By default the algorithm uses a sequential search, but parallelisation is also available.

## Installation

You can install the released version of **clustvarsel** from CRAN using:

```
install.packages("clustvarsel")
```

## Usage

Usage of the main functions and several examples are included in the
papers shown in the references section below.

For an intro see the vignette **A quick tour of clustvarsel**, which is available
as

```
vignette("clustvarsel")
```

The vignette is also available in the *Vignette* section on the navigation bar on top of the package's web page.

## References

Raftery, A. E. and Dean, N. (2006) Variable Selection for Model-Based Clustering. *Journal of the American Statistical Association*, 101(473), 168-178.

Maugis, C., Celeux, G., Martin-Magniette M. (2009) Variable Selection for Clustering With Gaussian Mixture Models. *Biometrics*, 65(3), 701-709.

Scrucca, L. and Raftery, A. E. (2018) clustvarsel: A Package Implementing Variable Selection for Gaussian Model-based Clustering in R. *Journal of Statistical Software*, 84(1), pp. 1-28.

