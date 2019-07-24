# Overview

This R package implements the multi-view Banded Spectral Clustering (mvBSC) algorithm developed by Luwan Zhang, Katherine Liao, Issac Kohane, and Tianxi Cai. The technical details of the mvBSC algorithm can be found [here](https://arxiv.org/abs/1804.02097).  

# Installation

If `devtools` is not installed, uncomment the code below and install it from CRAN.

``` r
# install.packages("devtools")
```

Install development version from GitHub:

``` r
devtools::install_github("celehs/mvBSC")
```

Load package `mvBSC` into R:

``` r
library(mvBSC)
```

# Main Steps

The main fitting function is `mvbsc`, which performs multi-view banded spectral clustering for a given set of parameters. To tune parameters, first apply the function `expand.grid` to create a data frame that contains possible combinations of candidate paramters, then run `mvbsc` for each parameter combination and record the evaluation metrics such as ARI, NMI, and F_1 measure. Choose the optimal setting which yields the best model performence (e.g., maximum ARI + NMI) and save the corresponding clustering results for printing and visualization.

# References

L. Zhang, K. Liao, I. Kohane, T. Cai. Multi-view Banded Spectral Clustering with Application to ICD9 Clustering. <https://arxiv.org/abs/1804.02097>
