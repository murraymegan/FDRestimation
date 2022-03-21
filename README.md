FDRestimation
========

`FDRestimation` is a user-friendly R package that directly computes and displays false discovery rates from p-values or z-scores under a variety of assumptions. 

Authors
-------
Megan Hollister Murray     
Vanderbilt University  
PhD Student, Department of Biostatistics  
<i class="fas fa-envelope"></i>  megan.c.hollister@vanderbilt.edu  
  
Jeffrey D. Blume  
Vanderbilt University  
Professor of Biostatistics, Biomedical Informatics and Biochemistry  
Vice Chair for Education, Biostatistics  
Director of Graduate Education, Data Science Institute  
<i class="fas fa-envelope"></i>  j.blume@vanderbilt.edu  

News
----
Version 1.0.1

Installation
------------

This package is availble on CRAN: https://CRAN.R-project.org/package=FDRestimation

``` r
install.packages("FDRestimation")
```


Or directly from GitHub:

``` r
# install.packages("devtools")
devtools::install_github("murraymegan/FDRestimation")
```

Example
-------

The `p.fdr()` function is used to compute FDRs and multiple-comparison adjusted p-values from a vector of raw p-values. This function allows for the following adjustment methods: Benjamini-Hochberg, Benjamini-Yeukateili (with both positive and negative correlation), Bonferroni, Holm, Hochberg, and Sidak. It also allows the user to specify the threshold for important findings, the assumed $pi_0$ value, the desired $pi_0$ estimation method, whether to sort the results, and whether to remove NAs in the imputed raw p-value vector count. 

The underlying methods for estimating the null proportion can be set by using the `estim.method` and `set.pi0` arguments. The default value of `set.pi0` is 1, meaning it assumes that all features are null features. Accordingly, this approach will yield conservative estimates of the FDR. Alternatively, and less conservatively, one can attempt to estimate the null proportion from the data. To do this, we recommend using "last.hist", as it was the simplest routine and one of the most accurate in our simulations. 

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

This `plot.p.fdr()` function is used to plot the results of `p.fdr`. By default, the adjusted FDRs, adjusted p-values and raw p-values are plotted along with two threshold lines to help contextualize the points. Any combination of p-values and thresholds can be removed from the plot. The user can set the axis limits, the location of the legend, the title of the plot and the plotting symbols and colors.

``` r
plot(fdr.output)
```

The `get.pi0()` function is used to estimate the null proportion from the raw p-values. The user can choose one of 6 different methods included in our function: Last Histogram Height, Storey, Meinshausen, Jiang, Nettleton, and Pounds. The user may also change the methods of determining the number of histogram breaks, which is an essential component for many of the methods implemented here. The histogram breaks method defaults to "scott" because this was the most accurate and reliable method in our simulatons. 

``` r 
get.pi0(pvalues=sim.data.p, estim.method="last.hist")
```

Help
----

Use the `help()` function to see full documentation on any of our functions. 

``` r 
help(p.fdr)
```

References
----------

A corresponding paper explaining and illustrating this package is in F1000 research.  
https://f1000research.com/articles/10-441


