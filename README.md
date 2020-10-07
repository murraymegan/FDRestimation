FDRestimation
========

The `FDRestimation` package contains functions to calculate method adjusted FDRs and p-values. 

Authors
-------
Megan Hollister Murray     
Vanderbilt University  
PhD Student, Department of Biostatistics  
Nashville, TN  
<i class="fas fa-envelope"></i>  megan.c.hollister@vanderbilt.edu  
  
Jeffrey D. Blume
Vanderbilt University  
Professor of Biostatistics, Biomedical Informatics and Biochemistry  
Vice Chair for Education, Biostatistics  
Director of Graduate Education, Data Science Institute  Vanderbilt University  
Nashville, TN  
<i class="fas fa-envelope"></i>  j.blume@vanderbilt.edu  

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

A corresponding paper explaining and illustarting this package is in the process of being submitted to The R Journal.

