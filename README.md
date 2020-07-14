FDRestimation
========

The `FDRestimation` package contains functions to calculate method adjusted FDRs and p-values. 

News
----
Version 1.0.0

Installation
------------

``` r
# install.packages("devtools")
devtools::install_github("murraymegan/FDRestimation")
```

Example
-------

The `p.fdr()` function calculates ?

``` r
library(FDRestimation)

pi0 <- 0.8
pi1 <- 1-pi0
n <- 10000
n.0 <- ceiling(n*pi0)
n.1 <- n-n.0

sim.data <- c(rnorm(n.1,5,1),rnorm(n.0,0,1))
sim.data.p <- 2*pnorm(-abs(sim.data))

fdr.output = p.fdr(pvalues=sim.data.p, adjust.method="BH")
```

References
----------

Paper appearing in ?

