pkgname <- "FDRestimation"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "FDRestimation-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('FDRestimation')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("get.pi0")
### * get.pi0

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: get.pi0
### Title: pi0 Estimation
### Aliases: get.pi0
### Keywords: FDR p-values

### ** Examples


# Example 1
pi0 = 0.8
pi1 = 1-pi0
n = 10000
n.0 = ceiling(n*pi0)
n.1 = n-n.0

sim.data = c(rnorm(n.1,3,1),rnorm(n.0,0,1))
sim.data.p = 2*pnorm(-abs(sim.data))

get.pi0(sim.data.p, estim.method = "last.hist")
get.pi0(sim.data.p, estim.method = "storey")
get.pi0(sim.data.p, estim.method = "set.pi0")




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("get.pi0", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("p.fdr")
### * p.fdr

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: p.fdr
### Title: FDR Computation
### Aliases: p.fdr
### Keywords: FDR p-values

### ** Examples


# Example 1
pi0 = 0.8
pi1 = 1-pi0
n = 10000
n.0 = ceiling(n*pi0)
n.1 = n-n.0

sim.data = c(rnorm(n.1,3,1),rnorm(n.0,0,1))
sim.data.p = 2*pnorm(-abs(sim.data))

fdr.output = p.fdr(pvalues=sim.data.p, adjust.method="BH")

fdr.output$fdrs
fdr.output$pi0

# Example 2

sim.data.p = output = c(runif(800),runif(200, min=0, max=0.01))
fdr.output = p.fdr(pvalues=sim.data.p, adjust.method="Holm", sort.results = TRUE)

fdr.output$`Results Matrix`




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("p.fdr", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plot.p.fdr")
### * plot.p.fdr

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plot.p.fdr
### Title: FDR plotting
### Aliases: plot.p.fdr
### Keywords: FDR p-values plot

### ** Examples


# Example 1

sim.data.p = c(runif(80),runif(20, min=0, max=0.01))
fdr.output = p.fdr(pvalues=sim.data.p)

plot(fdr.output)
plot(fdr.output, x.axis="Zvalues")





base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plot.p.fdr", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("print.summary.p.fdr")
### * print.summary.p.fdr

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: print.summary.p.fdr
### Title: Print the summary of p.fdr.object
### Aliases: print.summary.p.fdr
### Keywords: FDR p-values summary

### ** Examples


# Example 1
pi0 = 0.8
pi1 = 1-pi0
n = 10
n.0 = ceiling(n*pi0)
n.1 = n-n.0

sim.data = c(rnorm(n.1,5,1),rnorm(n.0,0,1))
sim.data.p = 2*pnorm(-abs(sim.data))

fdr.output = p.fdr(pvalues=sim.data.p, adjust.method="BH")

summary(fdr.output)





base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("print.summary.p.fdr", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("summary.p.fdr")
### * summary.p.fdr

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: summary.p.fdr
### Title: Summary of p.fdr.object
### Aliases: summary.p.fdr
### Keywords: FDR summary

### ** Examples


# Example 1
pi0 = 0.8
pi1 = 1-pi0
n = 10
n.0 = ceiling(n*pi0)
n.1 = n-n.0

sim.data = c(rnorm(n.1,5,1),rnorm(n.0,0,1))
sim.data.p = 2*pnorm(-abs(sim.data))

fdr.output = p.fdr(pvalues=sim.data.p, adjust.method="BH")

summary(fdr.output)





base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("summary.p.fdr", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
