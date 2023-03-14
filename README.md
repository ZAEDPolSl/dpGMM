# dpGMM
R package dpGMM for dynamic programming based Gaussian mixture modelling clustering for 1D and 2D data.

## Description
Package functionality:
1) Variety of criteria for component selection (AIC, AICc, BIC, ICL-BIC and LR),
2) Merging the components with a small deviation,
3) Quick stop if LR test does not show improvement. Possibility to control significance level,
4) Analysis of single measurements/vector as well as binned data,
5) Distribution plot with selected components in ggplot,
6) QQ-plot of fitted distribution and standard normal distribution,
7) Gaussian Mixture Modeling for 2D data.

## Installation
You can install the package from [GitHub](https://github.com/) with:
``` r
# install.packages("devtools")
devtools::install_github("ZAEDPolSl/rGMM")
```

## Manual
Tu damy tutorial jak używać
``` r
library(rGMM)
```

## References
[1] Polanski, Andrzej, et al. "Initializing the EM algorithm for univariate Gaussian, multi-component, heteroscedastic mixture models by dynamic programming partitions." International Journal of Computational Methods 15.03 (2018): 1850012.
