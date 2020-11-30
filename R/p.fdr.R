################################################################
##	Purpose: 	Compute FDRs and Adjusted p-values
##
##	Function:	p.fdr
##	Version:	1.0
##
##	Author:		Megan H. Murray and Jeffrey D. Blume
##	Date:		  November 5th, 2020
################################################################
#
#' FDR Computation
#'
#' @description This function computes FDRs and Method Adjusted p-values.
#'
#' @param pvalues A numeric vector of raw p-values.
#' @param zvalues A numeric vector of z-values to be used in pi0 estimation or a string with options "two.sided", "greater" or "less". Defaults to "two.sided".
#' @param threshold A numeric value in the interval \code{[0,1]} used in a multiple comparison hypothesis tests to determine significance from the null. Defaults to 0.05.
#' @param adjust.method A string used to identify the p-value and false discovery rate adjustment method. Defaults to \code{BH}. Options are \code{BH}, \code{BY}, code{Bon},\code{Holm}, \code{Hoch}, and \code{Sidak}.
#' @param BY.corr A string of either "positive" or "negative" to determine which correlation is used in the BY method. Defaults to \code{positive}.
#' @param just.fdr A Boolean TRUE or FALSE value which output only the FDR vector instead of the list output. Defaults to FALSE.
#' @param default.odds A numeric value determining the ratio of pi1/pi0 used in the computation of one FDR. Defaults to 1.
#' @param estim.method A string used to determine which method is used to estimate the null proportion or pi0 value. Defaults to \code{set.pi0}.
#' @param set.pi0 A numeric value to specify a known or assumed pi0 value in the interval \code{[0,1]}. Defaults to 1. Which means the assumption is that all inputted raw p-values come from the null distribution.
#' @param hist.breaks A numeric or string variable representing how many breaks are used in the pi0 estimation histogram methods. Defaults to "scott".
#' @param ties.method A string a character string specifying how ties are treated. Options are "first", "last", "average", "min", "max", or "random". Defaults to "random".
#' @param sort.results A Boolean TRUE or FALSE value which sorts the output in either increasing or non-increasing order dependent on the FDR vector. Defaults to FALSE.
#' @param na.rm A Boolean TRUE or FALSE value indicating whether NA's should be removed from the inputted raw p-value vector before further computation. Defaults to TRUE.
#'
#' @details We run into errors or warnings when pvalues, zvalues, threshold, set.pi0, BY.corr, or default.odds are not inputted correctly.
#'
#' @return A list containing the following components:
#' @return \item{fdrs}{A numeric vector of method adjusted FDRs.}
#' @return \item{Results Matrix}{A numeric matrix of method adjusted FDRs, method adjusted p-values, and raw p-values.}
#' @return \item{Reject Vector}{A vector containing Reject.H0 and/or FTR.H0 based off of the threshold value and hypothesis test on the adjusted p-values.}
#' @return \item{pi0}{A numeric value for the pi0 value used in the computations. }
#' @return \item{threshold}{A numeric value for the threshold value used in the hypothesis tests.}
#' @return \item{Adjustment Method}{The string with the method name used in computation(needed for the plot.fdr function).}
#'
#' @seealso \code{\link{plot.p.fdr}, \link{summary.p.fdr}, \link{get.pi0}}
#' @keywords  FDR p-values
#' @concept FDR adjusted p-values null proportion
#' @importFrom stats qnorm
#' @importFrom utils tail
#' @importFrom graphics hist
#' @importFrom stats smooth.spline
#' @importFrom stats predict
#' @importFrom graphics axis
#' @importFrom graphics points
#' @importFrom graphics abline
#' @importFrom graphics lines
#' @importFrom graphics legend
#' @importFrom graphics abline
#' @importFrom Rdpack reprompt
#' @export
#' @examples
#'
#' # Example 1
#' pi0 = 0.8
#' pi1 = 1-pi0
#' n = 10000
#' n.0 = ceiling(n*pi0)
#' n.1 = n-n.0
#'
#' sim.data = c(rnorm(n.1,3,1),rnorm(n.0,0,1))
#' sim.data.p = 2*pnorm(-abs(sim.data))
#'
#' fdr.output = p.fdr(pvalues=sim.data.p, adjust.method="BH")
#'
#' fdr.output$fdrs
#' fdr.output$pi0
#'
#' # Example 2
#'
#' sim.data.p = output = c(runif(800),runif(200, min=0, max=0.01))
#' fdr.output = p.fdr(pvalues=sim.data.p, adjust.method="Holm", sort.results = TRUE)
#'
#' fdr.output$`Results Matrix`
#'
#' @references
#' \insertRef{Rpack:bibtex}{Rdpack}
#'
#' \insertRef{R}{FDRestimation}
#'
#' \insertRef{efron:2013}{FDRestimation}
#'
#' \insertRef{bh:1995}{FDRestimation}
#'
#' \insertRef{shaffer:1995}{FDRestimation}
#'
#' \insertRef{wach:2004}{FDRestimation}
#'
#' \insertRef{storey:2003}{FDRestimation}
#'
#' \insertRef{by:2001}{FDRestimation}
#'
#' \insertRef{mein:2006}{FDRestimation}
#'
#' \insertRef{jiang:2008}{FDRestimation}
#'
#' \insertRef{nett:2006}{FDRestimation}
#'
#' \insertRef{pounds:2003}{FDRestimation}
#'
#' \insertRef{holm:1979}{FDRestimation}
#'
#' \insertRef{bon:1936}{FDRestimation}
#'
#' \insertRef{hoch:1988}{FDRestimation}
#'
#' \insertRef{sidak:1967}{FDRestimation}
#'
#' \insertRef{murray2020false}{FDRestimation}


p.fdr = function(pvalues,
                 zvalues = "two.sided",
                 threshold=0.05,
                 adjust.method="BH",
                 BY.corr="positive",
                 just.fdr=FALSE,
                 default.odds=1,
                 estim.method = "set.pi0",
                 set.pi0 = 1,
                 hist.breaks="scott",
                 ties.method = "random",
                 sort.results=FALSE,
                 na.rm=TRUE){

  requireNamespace(c("graphics", "stats", "utils"), quietly=TRUE)

  cl <- match.call()

  # Error Checking
  if(TRUE %in% (pvalues>1|pvalues<0)){
    stop("'pvalues' has value outside acceptable [0,1] range")
  }
  if(threshold>=1|threshold<=0){
    stop("'threshold' has value outside acceptable (0,1) range")
  }
  if(set.pi0>1|set.pi0<0){
    stop("'set.pi0' has value outside acceptable [0,1] range")
  }
  if(default.odds<0){
    stop("'default.odds' has a negative value which is outside its acceptable range")
  }
  if(!(BY.corr %in% c("positive", "negative"))){
    stop("'BY.corr' input not acceptable. Must be either 'positive' or 'negative'")
  }

  n=length(pvalues)

  #Remove NA inputted pvalues, zvalues
  if(na.rm){
    pvalues = pvalues[!is.na(pvalues)]
    zvalues = zvalues[!is.na(zvalues)]
    n = length(pvalues)
  }

  #Zvalues method
  if(is.character(zvalues)){
    if(zvalues=="greater"){
      zvalues = qnorm(pvalues, lower.tail = FALSE)
    }else if(zvalues=="two.sided"){
      zvalues = qnorm(pvalues/2, lower.tail = FALSE)
    }else if(zvalues=="less"){
      zvalues = qnorm(pvalues, lower.tail = TRUE)
    }
  }else{
    if(length(zvalues)!=length(pvalues)){
      stop("'zvalues' is different length than 'pvalues'")
    }
  }

  pi0 = get.pi0(pvalues=pvalues,
                set.pi0=set.pi0,
                estim.method=estim.method,
                zvalues=zvalues ,
                threshold=threshold,
                default.odds=default.odds,
                hist.breaks=hist.breaks)

  #Always calc the individual BH FDRs
  fdr.bh = pmin(1,pi0*pvalues*n/rank(pvalues,ties.method = ties.method))

  #Different Adjustment Methods
  if(adjust.method == "BH"){
    o = order(pvalues, decreasing = TRUE)
    ro = order(o)
    adj.pvalues = cummin(pmin(1,pvalues*n/rank(pvalues,ties.method = ties.method))[o])[ro]

    adj.fdrs = pmin(1,pi0*pvalues*n/rank(pvalues,ties.method = ties.method))
  }else if(adjust.method == "BY"&BY.corr=="positive"){
    dep = cumsum(1/(1:n))
    o = order(pvalues, decreasing = TRUE)
    ro = order(o)
    adj.pvalues = cummin(pmin(1,dep*pvalues*n/rank(pvalues,ties.method = ties.method))[o])[ro]

    adj.fdrs = pmin(1,pi0*dep*pvalues*n/rank(pvalues,ties.method = ties.method))
  }else if(adjust.method == "BY"&BY.corr=="negative"){
    gamma= -1*log(n)+sum(1/(1:n))
    dep = log(1:n)+gamma+1/(2*(1:n))
    o = order(pvalues, decreasing = TRUE)
    ro = order(o)
    adj.pvalues = cummin(pmin(1,dep*pvalues*n/rank(pvalues,ties.method = ties.method))[o])[ro]

    adj.fdrs = pmin(1,pi0*dep*pvalues*n/rank(pvalues,ties.method = ties.method))
  }else if(adjust.method == "Bon"){
    adj.pvalues = pmin(1, n*pvalues)

    adj.fdrs = pmin(1, pi0*n*pvalues)
  }else if(adjust.method == "Holm"){
    o = order(pvalues, decreasing = TRUE)
    ro = order(o)
    adj.pvalues = cummax(pmin(1,(n + 1 - rank(pvalues,ties.method = ties.method))*pvalues)[o])[ro]

    adj.fdrs = pmin(1, pi0*(n + 1 - rank(pvalues,ties.method = ties.method)*pvalues))
  }else if(adjust.method == "Hoch"){
    o = order(pvalues, decreasing = TRUE)
    ro = order(o)
    adj.pvalues = cummin(pmin(1,(n - rank(pvalues,ties.method = ties.method) + 1)*pvalues)[o])[ro]

    adj.fdrs = pmin(1, pi0*(n - rank(pvalues,ties.method = ties.method) + 1)*pvalues)
  }else if(adjust.method == "Sidak"){
    adj.pvalues = 1 - (1 - pvalues)^(n)
    adj.fdrs = pi0*(1 - (1 - pvalues)^(n))
  }

  #Sort results matrix in decreasing order by adj.fdrs
  if(sort.results){
    new.order = order(adj.fdrs, decreasing = FALSE)
  }else{
    new.order = c(1:n)
  }

  #Output
  if(length(pvalues)==1){
    one.fdr=(1+exp(zvalues^2/2)*default.odds)^-1

    class(one.fdr) = "one.fdr"
    return(one.fdr)

  }else if(just.fdr){
    output.fdrs = adj.fdrs[new.order]

    class(output.fdrs) = "fdrs"
    return(output.fdrs)

  }else{
    #Only report 4th column if different from first column
    if(adjust.method == "BH"){
      out.mat = as.data.frame(cbind("BH FDRs" = adj.fdrs[new.order],
                                 "Adjusted p-values" = adj.pvalues[new.order],
                                 "Raw p-values" = pvalues[new.order]))
    }else{
      out.mat = as.data.frame(cbind("FDRs" = adj.fdrs[new.order],
                                    "Adjusted p-values" = adj.pvalues[new.order],
                                 "Raw p-values" = pvalues[new.order],
                                 "BH FDRs" = fdr.bh[new.order]))

      colnames(out.mat)[1:2]= c(paste0(adjust.method," FDRs"),paste0(adjust.method," Adjusted p-values"))
    }

    #Return a list of fdrs, results matrix, pi0, threshold, adjustment method
    output = list("fdrs"=adj.fdrs[new.order],
         "Results Matrix"=out.mat,
         "Reject Vector"=ifelse(adj.pvalues[new.order]<=threshold,
                                "Reject.H0", "FTR.H0"),
         "pi0"=pi0,
         "threshold"=threshold,
         "Adjustment Method"= adjust.method,
         "Call" = cl)

    class(output) = "p.fdr"
    return(output)
  }
}

###
##
#
