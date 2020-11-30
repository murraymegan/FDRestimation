################################################################
##	Purpose: 	Summary of p.fdr.object
##
##	Function:	summary.p.fdr
##	Version:	1.0
##
##	Author:		Megan H. Murray and Jeffrey D. Blume
##	Date:		  November 5th, 2020
################################################################
#
#' Summary of p.fdr.object
#'
#' @description This function summarizes a p.fdr object.
#'
#' @param object A list of output from the p.fdr function.
#' @param digits A numeric value for the number of desired digits in the summary output. Defaults to 3.
#' @param ... Additional arguments affecting the summary produced.
#'
#' @details We run into errors or warnings when
#'
#' @return A list containing the following components:
#' @return \item{Range}{The range on the false discovery rates.}
#' @return \item{Significant Findings}{The number of significant findings. Found using the adjusted p-values and the given threshold. This is also the number of times we decide to reject the null hypothesis that the data is generated from a standard normal distribution.}
#' @return \item{Inconclusive Findings}{The number of inconclusive findings. Found using the adjusted p-values and the given threshold. This is also the number of times we fail to reject the null hypothesis that the data is generated from a standard normal distribution.}
#' @return \item{Assumed/Estimated pi0}{the assumed or estimated pi0 value depending on how the p.fdr function was run.}
#' @return \item{Number of Tests}{The total number of multiple comparison tests completed.}
#' @return \item{Adjustment Method}{The adjustment method used in the p.fdr function. }
#'
#' @seealso \code{\link{plot.p.fdr}, \link{p.fdr}, \link{get.pi0}}
#' @keywords summary FDR
#' @concept summary FDR adjusted p-values
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
#' fdr.output = p.fdr(pvalues=sim.data.p, adjust.method="BH")
#'
#' summary(fdr.output)
#'
#'
#' @references
#' \insertRef{Rpack:bibtex}{Rdpack}
#'
#' \insertRef{R}{FDRestimation}
#'
#' \insertRef{murray2020false}{FDRestimation}
#'
#' @rdname summary.p.fdr
#' @export

summary.p.fdr = function(object,
                         digits = 5,
                         ...){
  p.fdr.object=object
  if(is.null(p.fdr.object$fdrs)){
    stop("Invalid 'p.fdr' object: No 'fdrs' component. ")
  }

  ans = list("FDR Range" = round(range(p.fdr.object$fdrs),digits=digits),
             "Raw p-value Range" = round(range(p.fdr.object$`Results Matrix`$`Raw p-values`), digits=digits),
             "Significant Findings" = sum(p.fdr.object$`Results Matrix`$`Adjusted p-values` <= p.fdr.object$threshold),
             "Inconclusive Findings" = sum(p.fdr.object$`Results Matrix`$`Adjusted p-values` > p.fdr.object$threshold),
             "Assumed/Estimated pi0" = p.fdr.object$pi0,
             "Number of Tests" = length(p.fdr.object$fdrs),
             "Adjustment Method" = p.fdr.object$`Adjustment Method`,
             "Threshold" = p.fdr.object$threshold,
             "Call" =  p.fdr.object$Call)

  class(ans) <- "summary.p.fdr"
  return(ans)
}

###
##
#
