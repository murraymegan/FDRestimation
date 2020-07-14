################################################################
##	Purpose: 	Estimate Null Proportion of Data (pi0)
##
##	Function:	get.pi0
##	Version:	1.0
##
##	Author:		Megan H. Murray and Jeffrey D. Blume
##	Date:		  June 29, 2020
################################################################
#
#' pi0 Estimation
#'
#' @description This function estimates the null proportion of data or pi0 value.
#'
#' @param pvalues A numeric vector of raw p-values.
#' @param estim.method A string used to determine which method is used to estimate the pi0 value. Defaults to "last.hist".
#' @param zvalues A numeric vector of z-values to be used in pi0 estimation or a string with options "two.sided", "greater" or "less". Defaults to "two.sided".
#' @param threshold A numeric value in the interval \code{[0,1]} used in a multiple comparison hypothesis tests to determine significance from the null. Defaults to 0.05.
#' @param default.odds A numeric value determining the ratio of pi1/pi0 used in the computation of one FDR. Defaults to 1.
#' @param hist.breaks A numeric or string variable representing how many breaks in the pi0 estimation histogram methods. Defaults to "scott".
#'
#' @details We run into errors or warnings when
#'
#' @return \item{pi0} A numeric value represnting the proportion of the given data that come from the null distribution. A value in the interval \code{[0,1]}.
#'
#' @seealso \code{\link{plot.p.fdr}, \link{p.fdr}, \link{summary.p.fdr}}
#' @keywords
#' @export
#' @examples
#'
#' # Example 1
#' pi0 = 0.8
#' pi1 = 1-pi0
#' n = 10
#' n.0 = ceiling(n*pi0)
#' n.1 = n-n.0
#'
#' sim.data = c(rnorm(n.1,5,1),rnorm(n.0,0,1))
#' sim.data.p = 2*pnorm(-abs(sim.data))
#'
#' get.pi0(sim.data.p, estim.method = "last.hist")
#'
#' @references
#' R Journal 2020?
#'

get.pi0 = function(pvalues,
                   estim.method = "last.hist",
                   zvalues = "two.sided",
                   threshold=0.05,
                   default.odds=1,
                   hist.breaks="scott"){

  n=length(pvalues)

  #Zvalues method
  if(is.character(zvalues)){
    if(zvalues=="greater"){
      zvalues = qnorm(pvalues, lower.tail = FALSE)
    }else if(zvalues=="two.sided"){
      zvalues = qnorm(pvalues/2, lower.tail = FALSE)
    }else if(zvalues=="less"){
      zvalues = qnorm(pvalues, lower.tail = TRUE)
    }
  }

  #Null Proportion Estimation
  if(estim.method=="last.hist"){
    try.hist = hist(pvalues, breaks=hist.breaks, plot=FALSE)
    try.mids = try.hist$mids
    try.count = try.hist$counts

    if(tail(try.mids,1)<0.5){
      pi0=0
    }else{
      pi0 = min(tail(try.count,1)*length(try.mids)/sum(try.count),1)
    }
  }else if(estim.method=="Meinshausen"){
    ord <- order(pvalues)
    pvalues <- pvalues[ord]
    cutoff = threshold/n
    pvalues[pvalues < cutoff] <- cutoff

    find.beta <-function(m,alpha){
      quant.ecd <- -log(0.5*log(1/(1-alpha)))
      c.n <- 2*log(log(m))+0.5*log(log(log(m)))-0.5*log(4*pi)
      b.n <- sqrt(2*log(log(m)))
      beta <- 1/sqrt(m)*(quant.ecd+c.n)/b.n
    }

    get.boundingfunction.independent <- function(m,alpha,at=(1:1000)/1000,method="asymptotic")
    {
      if(method!="asymptotic") stop(paste("Method ", method, " not available"))
      beta <- find.beta(m,alpha)

      if( max(at)>1 | min(at)<0 ) stop( "bounding-function can only be evaluated within [0,1]")
      if(length(at)<1) stop(" bounding-function must be evaluated at least at one point ")

      boundingfunction <-   m*(  at + beta*sqrt(at*(1-at))  )
      return(boundingfunction)

    }

    boundingfunction <- get.boundingfunction.independent(n, threshold,
                                                         pvalues)
    lowerbound <- numeric(n)
    cummax <- 0
    for (p in 1:n) {
      cummax <- max(cummax, floor((p - floor(boundingfunction[p]))/max(0.2,
                                                                       (1 - pvalues[p]))))
      lowerbound[p] <- cummax
    }
    pi0=min(tail(lowerbound,1)/n,1)

  }else if(estim.method=="Storey"){
    #Default
    #lambda = seq(0.05, 0.95, by=0.05)
    #Change?
    #lambda = c(0.001,lambda, 0.999)

    lambda = seq(0, 0.95, by=0.05)
    lambda <- sort(lambda)
    ll <- length(lambda)

    if (ll == 1) {
      pi0 <- mean(pvalues >= lambda)/(1 - lambda)
      pi0.lambda <- pi0
      pi0 <- min(pi0, 1)
      pi0Smooth <- NULL
    }else{
      ind <- length(lambda):1

      # We add in (nbins=ll)
      pi0 <- cumsum(tabulate(findInterval(pvalues , vec = lambda), nbins=ll)[ind])/(length(pvalues ) *
                                                                                      (1 - lambda[ind]))
      pi0 <- pi0[ind]
      pi0.lambda <- pi0

      spi0 <- smooth.spline(lambda, pi0, df = 3)
      #Original
      #pi0Smooth <- predict(spi0, x = lambda)$y

      pi0Smooth <- predict(spi0, x = c(lambda,1))$y

      pi0 <- max(min(pi0Smooth[ll+1], 1),0)
    }
  }else if(estim.method=="Pounds"){
    pi0 = min(1, 2 * mean(pvalues))
  }else if(estim.method=="Jiang"){
    pi0.jiang = function(p, nbin) {
      m = length(p)
      t = seq(0, 1, length = nbin + 1)
      NB = rep(0, nbin)
      NBaverage = rep(0, nbin)
      NS = rep(0, nbin)
      pi = rep(0, nbin)
      for (i in 1:nbin) {
        NB[i] = length(p[p >= t[i]])
        NBaverage[i] = NB[i]/(nbin - (i - 1))
        NS[i] = length(p[p >= t[i]]) - length(p[p >=
                                                  t[i + 1]])
        pi[i] = NB[i]/(1 - t[i])/m
      }
      i = min(which(NS <= NBaverage))
      pi0 = min(1, mean(pi[(i - 1):nbin]))
      return(pi0)
    }
    histo=hist(pvalues,breaks=hist.breaks, plot=FALSE)
    nbins=length(histo$mids)
    pi0 = pi0.jiang(pvalues, nbins)
  }else if(estim.method=="Nettleton"){
    pi0.histo = function(p, nbin) {
      bin = c(-0.1, (1:nbin)/nbin)
      bin.counts = tabulate(cut(p, bin), nbins=nbin)
      tail.means = rev(cumsum(rev(bin.counts))/(1:nbin))
      index = which(tail.means >= bin.counts)[1]
      tail.means[index]/tail.means[1]
    }
    histo=hist(pvalues,breaks=hist.breaks, plot=FALSE)
    nbins=length(histo$mids)
    pi0 = min(pi0.histo(pvalues, nbins),1)
  }
    return(pi0)
}

###
##
#