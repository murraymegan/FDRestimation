################################################################
##	Purpose: 	Compute FDRs and Adjusted P-Values
##
##	Function:	p.fdr
##	Version:	1.0
##
##	Author:		Megan H. Murray and Jeffrey D. Blume
##	Date:		  March 11, 2020
################################################################
#
#' FDR Computation
#'
#' @description This function computes FDRs and Method Adjusted P-Values.
#'
#' @param
#'
#' @details
#'
#' @return A list containing the following components:
#' \describe{
#' \item{\code{fdr}}{Vector of method adjusted fdrs.}
#'
#' }
#' @seealso \code{\link{plotfdr}}
#' @keywords szjbfajbfja
#' @export
#' @examples
#'
#' ## Example
#' pi0 <- 0.8
#' pi1 <- 1-pi0
#' n <- 10000
#' n.0 <- ceiling(n*pi0)
#' n.1 <- n-n.0
#'
#' data <- c(rnorm(n.1,5,1),rnorm(n.0,0,1))
#' sim.data.p <- 2*pnorm(-abs(sim.data))
#'
#' fdr.output = p.fdr(pvalues=sim.data.p, adjust.method="BH")
#'
#' fdr.output$fdr
#' [1] ?
#' sgpv$pi0
#' [1]    ?
#'
#'
#' @references
#'
#'

p.fdr = function(pvalues,
                 zvalues = NA,
                 threshold=0.05,
                 adjust.method=NA,
                 just.fdr=FALSE,
                 pi0.estim = "one",
                 set.pi0 = NA,
                 na.rm=TRUE){

  library(splines)
  n=length(pvalues)

  if(na.rm){
    pvalues = pvalues[!is.na(pvalues)]
    n = length(pvalues)
  }

  #Null Proportion Estimation
  if(pi0.estim=="one"){
    pi0=1
  }
  else if(pi0.estim=="last.hist"){
    try.hist <- hist(pvalues, breaks=min(n/10,100), plot=FALSE)
    try.mids <- try.hist$mids
    try.count <- try.hist$counts
    pi0 = tail(try.count,1)*length(try.mids)/sum(try.count)
  }
  else if(pi0.estim=="loess"){
    try.hist <- hist(pvalues, breaks=n/10, plot=FALSE)
    try.mids <- try.hist$mids
    try.count <- try.hist$counts
    loe <- loess(try.count~try.mids)
    k <- 0.3
    x.loe <- try.mids[try.mids>=k]
    y.loe <- loe$fitted[try.mids>=k]
    avg.loe <- mean(y.loe)
    y.est <- c(rep(avg.loe,(length(try.mids)-length(y.loe))),y.loe)

    pi0 = sum(y.est)/sum(try.count)
  }
  else if(pi0.estim=="lindseys"){
    if(sum(is.na(zvalues))==1){
      stop("No argument for `zvalues` provided")
    }
    sim.hist <-hist(zvalues, plot=FALSE, breaks=min(n/10,100))
    bin <- length(sim.hist$mids)  # number of bins
    df <- 7 # J= degrees of freedom
    breaks <- sim.hist$breaks
    y <- sim.hist$counts
    x <- sim.hist$mids
    K <- length(y)
    k <- seq(K)
    yhat <- glm(y ~ ns(x,df=df), poisson)$fit
    bw <- x[2]-x[1]
    fhat <- yhat/(n*bw)
    #Lindseys
    fdra <- dnorm(x,0,1)/fhat

    plot(x, pmin(fdra,1), type="l", col="blue",
         main="Lindsey's adjust.method FDR vs. z-values")
    abline(h=0.05, col="red", lty=2)

    pi0 = min(quantile(fhat/dnorm(x,0,1),0.05),1)
  }

  #Always calc the individual FDRs
  fdr.bh = pmin(1,pi0*pvalues*n/rank(pvalues))

  if(is.na(adjust.method)| adjust.method == "BH"){
    adj.pvalues = cummin(pmin(1,pvalues*n/rank(pvalues)))

    adj.fdrs = pmin(1,pi0*pvalues*n/rank(pvalues))
  }
  else if(adjust.method == "BY"){
    dep = cumsum(1/(1:n))

    adj.pvalues = cummin(pmin(1,dep*pvalues*n/rank(pvalues)))

    adj.fdrs = pmin(1,pi0*dep*pvalues*n/rank(pvalues))
  }
  else if(adjust.method == "Bon"){
    adj.pvalues = pmin(1, n*pvalues)

    adj.fdrs = pmin(1, pi0*n*pvalues)
  }
  else if(adjust.method == "Holm"){
    adj.pvalues = cummax(pmin(1,(n - rank(pvalues) + 1)*pvalues))

    adj.fdrs = pmin(1, (n - rank(pvalues) + 1)*pvalues)
  }

  out <- as.data.frame(cbind("Method Adjusted FDRs" = adj.fdrs,
                             "Method Adjusted P-Values" = adj.pvalues,
                             "Unadjusted P-Values" = pvalues,
                             "Individual BH FDRs" = fdr.bh))

  return(list("fdrs"=adj.fdrs,
              "Results Matrix"=out,
              "Reject Vector"=ifelse(adj.pvalues<=threshold,
                                     "Reject.H0", "FTR.H0"),
              "pi0"=pi0,
              "threshold"=threshold))
}


###
##
#
