sgpv
========

The `sgpv` package contains functions to calculate second-generation p-values, their associated delta-gaps, and the false discovery risk or false confirmation risk for an alternative or null finding (SGPV = 0 or SGPV = 1) when assumptions are made about the distributions over the null and alternative spaces. It also contains several functions for a variety of plotting types relevant to SGPV usage.

News
----
Version 1.1.0 updated November 2020, with newly added functions `plotman` and `plotsgpower`.


Installation
------------

From CRAN:

``` r
install.packages("sgpv")
```

From GitHub: 

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

Introductory paper appearing in the special issue of The American Statisician:

Jeffrey D. Blume, Robert A. Greevy, Valerie F. Welty, Jeffrey R. Smith & William D. Dupont (2019) An Introduction to Second-Generation p-Values, The American Statistician, 73:sup1, 157-167, https://doi.org/10.1080/00031305.2018.1537893 

Original proposal appearing in PLoS ONE:

Blume JD, D’Agostino McGowan L, Dupont WD, Greevy RA Jr. (2018). Second-generation p-values: Improved rigor, reproducibility, & transparency in statistical analyses. PLoS ONE 13(3): e0188299. https://doi.org/10.1371/journal.pone.0188299

