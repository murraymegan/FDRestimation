################################################################
##	Purpose: 	Print summary of p.fdr.object
##
##	Function:	print.summary.p.fdr
##	Version:	1.0
##
##	Author:		Megan H. Murray and Jeffrey D. Blume
##	Date:		  November 5th, 2020
################################################################
#
#' Print the summary of p.fdr.object
#'
#' @description This function prints the summary a p.fdr.object.
#'
#' @param x A list of output from the summary.p.fdr function.
#' @param digits A numeric value for the number of desired digits in the summary output. Defaults to 3.
#' @param ... Further arguments passed to or from other methods.
#'
#' @details We run into errors or warnings when
#'
#' @seealso \code{\link{plot.p.fdr}, \link{p.fdr}, \link{get.pi0}}
#' @keywords  summary FDR p-values
#' @concept summary FDR adjusted p-values null proportion
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
#' @rdname print.summary.p.fdr
#' @export

print.summary.p.fdr = function(x,
                               digits = 3,
                               ...){

  summary.p.fdr.object = x
  if (is.null(summary.p.fdr.object)){
    stop("invalid 'summary.p.fdr' object")
  }

  cat("\nCall:\n")
  print(summary.p.fdr.object$Call)

  cat("\n","Number of tests: ", summary.p.fdr.object$`Number of Tests`, "\n", sep="")
  cat("Raw p-value Range: [", summary.p.fdr.object$`Raw p-value Range`[1], ", ",summary.p.fdr.object$`Raw p-value Range`[2],"]", "\n", sep="")
  cat("\n", "Adjustment Method: ", summary.p.fdr.object$`Adjustment Method`, "\n", sep="")
  cat("False Discovery Rate Range: [", summary.p.fdr.object$`FDR Range`[1], ", ",summary.p.fdr.object$`FDR Range`[2],"]", "\n", sep="")
  cat("Findings at ", summary.p.fdr.object$`Threshold`, " level:", "\n", sep="")
  cat("   Significant (Reject): ", summary.p.fdr.object$`Significant Findings`, "\n", sep="")
  cat("   Inconclusive (Fail to Reject): ", summary.p.fdr.object$`Inconclusive Findings`, "\n", sep="")
  cat("---\n")
  cat("Estimated/Assumed Null Proportion (pi0): ", summary.p.fdr.object$`Assumed/Estimated pi0` , "\n",sep="")
  cat("\n")


}

###
##
#
