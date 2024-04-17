# dpGMM
R package dpGMM for dynamic programming based Gaussian mixture modelling clustering for 1D and 2D data.

## Description
Package functionality:
1) Variety of criteria for selection component number (AIC, AICc, BIC, ICL-BIC and LR),
2) Posibility of merging the components with a small satndard deviation within each other,
3) Control of minimum allow variance of component to avoid picks,
4) Quick stop if LR test does not show improvement (off/on) accompanien by possibility to control significance level,
5) Analysis of single measurements/vector as well as binned data,
6) Distribution plot with selected components in ggplot,
7) QQ-plot of fitted distribution and standard normal distribution,
8) Gaussian Mixture Modeling for 2D data.

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

data(example)
mix_test <- runGMM(example$Dist)

```

## References
[1] Polanski, Andrzej, et al. "Initializing the EM algorithm for univariate Gaussian, multi-component, heteroscedastic mixture models by dynamic programming partitions." International Journal of Computational Methods 15.03 (2018): 1850012.
