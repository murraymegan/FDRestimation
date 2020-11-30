################################################################
##	Purpose: 	Plot FDRs
##
##	Function:	plot.p.fdr
##	Version:	1.0
##
##	Author:		Megan H. Murray and Jeffrey D. Blume
##	Date:		  November 5th, 2020
################################################################
#
#' FDR plotting
#'
#' @description This function creates a plot using a x (p.fdr.object).
#'
#' @param x A p.fdr object that contains the list of output.
#' @param raw.pvalues A Boolean TRUE or FALSE value to indicate whether or not to plot the raw p-value points. Defaults to TRUE.
#' @param adj.pvalues A Boolean TRUE or FALSE value to indicate whether or not to plot the adjusted p-value points. Defaults to TRUE.
#' @param sig.line A Boolean TRUE or FALSE value to indicate whether or not to plot the raw p-value significance line. Defaults to TRUE.
#' @param adj.sig.line A Boolean TRUE or FALSE value to indicate whether or not to plot the adjusted significance threshold. Defaults to TRUE.
#' @param threshold A numeric value to determine the threshold at which we plot significance. Defaults to value used in the p.fdr.object.
#' @param x.axis A string variable to indicate what to plot on the x-axis. Can either be "Rank" or "Zvalues". Defaults to "Rank".
#' @param xlim A numeric interval for x-axis limits.
#' @param ylim A numeric interval for y-axis limits. Defaults to c(0,1).
#' @param zvalues A numeric vector of z-values to be used in pi0 estimation or a string with options "two.sided", "greater" or "less". Defaults to "two.sided".
#' @param legend.where A string "bottomright", "bottomleft", "topleft", "topright". Defaults to "topleft" is x.axis="Rank" and "topright" if x.axis="Zvalues".
#' @param legend.on A Boolean TRUE or FALSE value to indicate whether or not to print the legend.
#' @param main A string variable for the title of the plot.
#' @param pch.adj.p A plotting "character’, or symbol to use for the adjusted p-value points. This can either be a single character or an integer code for one of a set of graphics symbols. Defaults to 17.
#' @param pch.raw.p A plotting "character’, or symbol to use for the raw p-value points. This can either be a single character or an integer code for one of a set of graphics symbols. Defaults to 20.
#' @param pch.adj.fdr A plotting "character’, or symbol to use for the adjusted FDR points. This can either be a single character or an integer code for one of a set of graphics symbols. Defaults to 20.
#' @param col A vector of colors for the points and lines in the plot. If the input has 1 value all points and lines will be that same color. If the input has length of 3 then col.adj.fdr will be the first value, col.adj.p will be the second, and col.raw.p is the third. Defaults to c("dodgerblue","firebrick2", "black").
#' @param ... Graphical parameters. Any argument that can be passed to image.plot and to base plot, such as axes=FALSE, main='title', ylab='latitude'
#'
#' @details We run into errors or warnings when zvalues or col are inputted incorrectly.
#'
#' @seealso \code{\link{summary.p.fdr}, \link{p.fdr}, \link{get.pi0}}
#' @keywords plot FDR p-values
#' @concept  plot FDR adjusted p-values
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
#' @examples
#'
#' # Example 1
#'
#' sim.data.p = c(runif(80),runif(20, min=0, max=0.01))
#' fdr.output = p.fdr(pvalues=sim.data.p)
#'
#' plot(fdr.output)
#' plot(fdr.output, x.axis="Zvalues")
#'
#'
#' @references
#' \insertRef{Rpack:bibtex}{Rdpack}
#'
#' \insertRef{R}{FDRestimation}
#'
#' \insertRef{bh:1995}{FDRestimation}
#'
#' \insertRef{by:2001}{FDRestimation}
#'
#' \insertRef{holm:1979}{FDRestimation}
#'
#' \insertRef{hoch:1988}{FDRestimation}
#'
#' \insertRef{sidak:1967}{FDRestimation}
#'
#' \insertRef{bon:1936}{FDRestimation}
#'
#' \insertRef{murray2020false}{FDRestimation}
#'
#' @rdname plot.p.fdr
#' @export


plot.p.fdr = function(x,
                      raw.pvalues=TRUE,
                      adj.pvalues=TRUE,
                      sig.line=TRUE,
                      adj.sig.line=TRUE,
                      threshold = NA,
                      x.axis="Rank",
                      xlim=NA,
                      ylim=c(0,1),
                      zvalues="two.sided",
                      legend.where=NA,
                      legend.on=TRUE,
                      main=NA,
                      pch.adj.p=17,
                      pch.raw.p=20,
                      pch.adj.fdr=20,
                      col=c("dodgerblue","firebrick2", "black"),
                      ...){

  requireNamespace(c("graphics", "stats", "utils"), quietly=TRUE)

  p.fdr.object=x
  n=length(p.fdr.object$fdrs)

  if(is.na(threshold)){
    threshold = p.fdr.object$threshold
  }

  #Zvalues method
  if(is.character(zvalues)){
    if(zvalues=="greater"){
      zvalues = qnorm(p.fdr.object$`Results Matrix`[,3], lower.tail = FALSE)
    }else if(zvalues=="two.sided"){
      zvalues = qnorm(p.fdr.object$`Results Matrix`[,3]/2, lower.tail = FALSE)
    }else if(zvalues=="less"){
      zvalues = qnorm(p.fdr.object$`Results Matrix`[,3], lower.tail = TRUE)
    }
  }else{
    if(length(zvalues)!=length(p.fdr.object$`Results Matrix`[,3])){
      stop("'zvalues' is different length than 'p.fdr object'")
    }
  }

  if(length(col)==3){
    col.adj.fdr = col[1]
    col.adj.p = col[2]
    col.raw.p = col[3]
  }else if(length(col)==1){
    col.adj.fdr=col.adj.p=col.raw.p = col
  }else{
    stop("Length of 'col' input is not 1 or 3.")
  }

  if(x.axis=="Rank"){
    if(is.na(legend.where)){
      legend.where="topleft"
    }

    if(length(xlim)!=2){
      xlim = c(1,n)
    }
    if(p.fdr.object$`Adjustment Method`=="BH"){
      x.line=0:(n)
      y.line=threshold*(0:(n))/n
    }
    else if(p.fdr.object$`Adjustment Method`=="BY"){
      x.line=0:n
      y.line=0.5*(0:n)/(n*cumsum(1/(0:n)))
    }
    else if(p.fdr.object$`Adjustment Method`=="Hoch"){
      x.line=0:n
      y.line=threshold/(n+1-(0:n))
    }
    else if(p.fdr.object$`Adjustment Method`=="Bon"){
      x.line=0:n
      y.line=rep(threshold/n,n+11)
    }
    else if(p.fdr.object$`Adjustment Method`=="Holm"){
      x.line=0:n
      y.line=threshold/(n+1-(0:n))
    }else{
      x.line=NA
      y.line=NA
    }

    if(is.na(main)){
      main= paste0(p.fdr.object$`Adjustment Method`,
                    " p.fdr Object Plot")
    }
    if(adj.pvalues){
      plot(rank(p.fdr.object$`Results Matrix`[,3],ties.method = "random"),
           p.fdr.object$`Results Matrix`[,2],
           col=col.adj.fdr,
           ylim=ylim,
           xlim=xlim,
           xlab="Ranking of Raw p-values",
           ylab=" ",
           main=main,
           pch=pch.adj.p,
           cex.main=1.2,
           axes = FALSE,
           cex.axis=0.5)
      if(adj.sig.line){
        axis(side = 1)
        axis(side = 2,
             at = c(threshold),
             las=1,
             col.axis = col.adj.fdr,
             cex.axis=0.75)
        axis(side = 2,
             cex.axis=0.75,
             las=1)
        points(rank(p.fdr.object$`Results Matrix`[,3],ties.method = "random"),
               p.fdr.object$`Results Matrix`[,1],
               col=col.adj.p,
               pch=pch.adj.fdr)
        abline(h=threshold,
               col=col.adj.fdr,
               lwd=2)
      }else{
        axis(side = 1)
        axis(side = 2,
             las=1)
        points(rank(p.fdr.object$`Results Matrix`[,3],ties.method = "random"),
               p.fdr.object$`Results Matrix`[,1],
               col=col.adj.p,
               pch=pch.adj.fdr)
      }
    }else{
      plot(rank(p.fdr.object$`Results Matrix`[,3],ties.method = "random"),
           p.fdr.object$`Results Matrix`[,1],
           col=col.adj.p,
           ylim=ylim,
           xlim=xlim,
           xlab="Ranking of Raw p-values",
           ylab=" ",
           main=main,
           pch=pch.adj.fdr,
           las=1)
    }
    if(raw.pvalues){
      points(rank(p.fdr.object$`Results Matrix`[,3],ties.method = "random"),
             p.fdr.object$`Results Matrix`[,3],
             col=col.raw.p,
             pch=pch.raw.p)
      if(sig.line){
        lines(x.line, y.line,
              col=col.raw.p,
              lwd=2)
      }
    }

    #Legends
    if(raw.pvalues&adj.pvalues&sig.line&adj.sig.line&legend.on){
      legend(legend.where,
             legend=c(paste0("Adj p-values (",p.fdr.object$`Adjustment Method`,")"),
                      paste0("FDRs (",p.fdr.object$`Adjustment Method`,")"),
                      "Raw p-values",
                      paste0(p.fdr.object$`Adjustment Method`," Threshold (Raw p-values)"),
                      paste0(p.fdr.object$`Adjustment Method`," Threshold (Adj p-values)")),
             pch=c(pch.adj.p,pch.adj.fdr,pch.raw.p, NA, NA),
             lty=c(NA,NA,NA,1,1),
             lwd=c(NA,NA,NA,2,2),
             col=c(col.adj.fdr,col.adj.p, col.raw.p, col.raw.p,col.adj.fdr),
             cex=0.7,
             bty = "n")
    }
    else if(raw.pvalues&!adj.pvalues&sig.line&!adj.sig.line&legend.on){
      legend(legend.where,
             legend=c(paste0("FDRs (", p.fdr.object$`Adjustment Method`,")"),
                      "Raw p-values",
                      paste0(p.fdr.object$`Adjustment Method`," Threshold (Raw p-values)")),
             pch=c(pch.adj.fdr,pch.raw.p, NA),
             lty=c(NA,NA,1),
             lwd=c(NA,NA,2),
             col=c(col.adj.p, col.raw.p, col.raw.p),
             cex=0.7,
             bty = "n")
    }
    else if(!raw.pvalues&adj.pvalues&!sig.line&adj.sig.line&legend.on){
      legend(legend.where,
             legend=c(paste0("Adj p-values (",p.fdr.object$`Adjustment Method`,")"),
                      paste0("FDRs (", p.fdr.object$`Adjustment Method`,")"),
                      paste0(p.fdr.object$`Adjustment Method`," Threshold (Adj p-values)")),
             pch=c(pch.adj.p,pch.adj.fdr, NA),
             lty=c(NA,NA,1),
             lwd=c(NA,NA,2),
             col=c(col.adj.fdr,col.adj.p, col.adj.fdr),
             cex=0.7,
             bty = "n")
    }
    else if(sig.line&legend.on){
      legend(legend.where,
             legend=c(paste0("FDRs (",p.fdr.object$`Adjustment Method`,")")),
             pch=c(pch.adj.fdr),
             col=c(col.adj.p),
             cex=0.7,
             bty = "n")

    }

  }else if(x.axis=="Zvalues"){
      if(is.na(legend.where)){
        legend.where="topright"
      }

      if(is.na(main)){
        main= paste0(p.fdr.object$`Adjustment Method`,
                      " p.fdr Z-values Plot")
      }

      if(adj.pvalues){
        plot(zvalues,
             p.fdr.object$`Results Matrix`[,2],
             col=col.adj.fdr,
             ylim=ylim,
             xlab="Z-values",
             ylab=" ",
             main=main,
             pch=pch.adj.p,
             cex.main=1.2,
             axes=FALSE,
             cex.axis=0.5,
             las=1)
        axis(side = 1)
        axis(side = 2,
             at = c(threshold),
             las=1,
             col.axis=col.adj.fdr,
             cex.axis=0.75)
        axis(side = 2,
             cex.axis=0.75,
             las=1)
        points(zvalues,
               p.fdr.object$`Results Matrix`[,1],
               col=col.adj.p,
               pch=pch.adj.fdr)
        abline(h=threshold,
               col=col.adj.fdr,
               lwd=2)
      }else{
        plot(zvalues,
             p.fdr.object$`Results Matrix`[,1],
             col=col.adj.p,
             ylim=ylim,
             xlab="Z-values",
             ylab=" ",
             main=main,
             pch=pch.adj.fdr,
             las=1)
      }

      if(raw.pvalues){
        points(zvalues,
               p.fdr.object$`Results Matrix`[,3],
               col=col.raw.p,
               pch=pch.raw.p)
      }

    #Legends
    if(raw.pvalues&adj.pvalues&sig.line&adj.sig.line){
      legend(legend.where,
             legend=c(paste0("Adj p-values (",p.fdr.object$`Adjustment Method`,")"),
                      paste0("FDRs (",p.fdr.object$`Adjustment Method`,")"),
                      "Raw p-values",
                      paste0(p.fdr.object$`Adjustment Method`," Threshold (Adj p-values)")),
             pch=c(pch.adj.p,pch.adj.fdr,pch.raw.p, NA),
             lty=c(NA,NA,NA,1),
             lwd=c(NA,NA,NA,2),
             col=c(col.adj.fdr,col.adj.p, col.raw.p,col.adj.fdr),
             cex=0.7,
             bty = "n")
    }
    else if(raw.pvalues&!adj.pvalues&sig.line&!adj.sig.line){
      legend(legend.where,
             legend=c(paste0("FDRs (", p.fdr.object$`Adjustment Method`,")"),
                      "Raw p-values"),
             pch=c(pch.adj.fdr,pch.raw.p),
             col=c(col.adj.p, col.raw.p),
             cex=0.7,
             bty = "n")
    }
    else if(!raw.pvalues&adj.pvalues&!sig.line&adj.sig.line){
      legend(legend.where,
             legend=c(paste0("Adj p-values (",p.fdr.object$`Adjustment Method`,")"),
                      paste0("FDRs (", p.fdr.object$`Adjustment Method`,")"),
                      paste0(p.fdr.object$`Adjustment Method`," Threshold (Raw p-values)")),
             pch=c(pch.adj.p,pch.adj.fdr, NA),
             lty=c(NA,NA,1),
             lwd=c(NA,NA,2),
             col=c(col.adj.fdr,col.adj.p, col.adj.fdr),
             cex=0.7,
             bty = "n")
    }
    else{
      legend(legend.where,
             legend=c(paste0("FDRs (",p.fdr.object$`Adjustment Method`,")")),
             pch=c(pch.adj.fdr),
             col=c(col.adj.p),
             cex=0.7,
             bty = "n")
    }
  }

}



###
##
#
