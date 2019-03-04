sgpv
========

The `sgpv` package contains functions to calculate the second-generation p-value, its associated delta-gap, and the false discovery risk or false confirmation risk for an alternative or null finding (sgpv = 0 or sgpv = 1) when assumptions are made about the distributions over the null and alternative spaces.

News
----
Version 0.0.1.0000 coming Spring 2019!


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
