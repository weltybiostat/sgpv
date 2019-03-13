sgpv
========

The `sgpv` package contains functions to calculate the second-generation p-value in R, its associated delta-gap, power functions, and the false discovery risk or false confirmation risk for an alternative or null finding (sgpv = 0 or sgpv = 1) when assumptions are made about the distributions over the null and alternative spaces.

News
----
Version 0.0.1.0000 now includes a function to plot the intervals colored by SGPV status, as found in Blume et al. (2018) and Blume et al. (2019)

Installation
------------

``` r
# install.packages("devtools")
devtools::install_github("weltybiostat/sgpv")
```

Example
-------

The `sgpvalue()` function calculates the second-generation p-value and delta-gap (if applicable) for uncertainty intervals with lower bounds `est.lo` and upper bounds `est.hi` and an indifference zone (i.e. interval null hypothesis) of (`null.lo`, `null.hi`).  Note that this example is in terms of odds ratios, and the second-generation p-value should be calculated on the "symmetric" scale, i.e. log odds ratios in this case.

``` r
library(sgpv)
lb = log(c(1.05, 1.3, 0.97))
ub = log(c(1.8, 1.8, 1.02))
sgpvalue(est.lo = lb, est.hi = ub, null.lo = log(1/1.1), null.hi = log(1.1))

# $p.delta
# [1] 0.1220227 0.0000000 1.0000000

# $delta.gap
# [1]       NA 1.752741       NA
```

References
----------

Paper appearing in the upcoming special issue of The American Statisician:

Blume JD, Greevy RA Jr., Welty VF, Smith JR, Dupont WD (2019). An Introduction to Second-generation p-values. The American Statistician. In press. https://doi.org/10.1080/00031305.2018.1537893

Original proposal appearing in PLoS ONE:

Blume JD, D’Agostino McGowan L, Dupont WD, Greevy RA Jr. (2018). Second-generation p-values: Improved rigor, reproducibility, & transparency in statistical analyses. PLoS ONE 13(3): e0188299. https://doi.org/10.1371/journal.pone.0188299

