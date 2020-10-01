################################################################
##	Purpose: 	Plot FDRs
##
##	Function:	plot.p.fdr
##	Version:	1.0
##
##	Author:		Megan H. Murray and Jeffrey D. Blume
##	Date:		  June 29, 2020
################################################################
#
#' FDR plotting
#'
#' @description This function creates a plot using a p.fdr.object.
#'
#' @param p.fdr.object A p.fdr object that contains the list of output.
#' @param raw.pvalues A Boolean TRUE or FALSE value to indicate whether or not to plot the raw p-value points. Defaults to TRUE.
#' @param adj.pvalues A Boolean TRUE or FALSE value to indicate whether or not to plot the adjusted p-value points. Defaults to TRUE.
#' @param sig.line A Boolean TRUE or FALSE value to indicate whether or not to plot the raw p-value significance line. Defaults to TRUE.
#' @param adj.sig.line A Boolean TRUE or FALSE value to indicate whether or not to plot the adjusted significance threshold. Defaults to TRUE.
#' @param threshold A numeric value to determine the threshold at which we plot significance. Defaults to value used in the p.fdr.object.
#' @param x.axis A string variable to indicate what to plot on the x-axis. Can either be "Rank" or "Zvalues". Defaults to "Rank".
#' @param x.lim A numeric interval for x-axis limits.
#' @param y.lim A numeric interval for y-axis limits. Defaults to c(0,1).
#' @param zvalues A numeric vector of z-values to be used in pi0 estimation or a string with options "two.sided", "greater" or "less". Defaults to "two.sided".
#' @param legend.where A string "bottomright", "bottomleft", "topleft", "topright". Defaults to "topleft" is x.axis="Rank" and "topright" if x.axis="Zvalues".
#' @param title A string variable for the title of the plot.
#' @param pch.star A plotting ‘character’, or symbol to use for the adjusted p-value points. This can either be a single character or an integer code for one of a set of graphics symbols. Defaults to 17.
#'
#' @details We run into errors or warnings when
#'
#' @seealso \code{\link{summary.p.fdr}, \link{p.fdr}, \link{get.pi0}}
#' @keywords
#' @export
#' @examples
#'
#' # Example 1
#'
#' sim.data.p = c(runif(80),runif(20, min=0, max=0.01))
#' fdr.output = p.fdr(pvalues=sim.data.p, sort.results = TRUE)
#'
#' plot(fdr.output)
#' plot(fdr.output, x.axis="Zvalues")
#'
#' @references
#'

plot.p.fdr = function(p.fdr.object,
                      raw.pvalues=TRUE,
                      adj.pvalues=TRUE,
                      sig.line=TRUE,
                      adj.sig.line=TRUE,
                      threshold = NA,
                      x.axis="Rank",
                      x.lim=NA,
                      y.lim=c(0,1),
                      zvalues="two.sided",
                      legend.where=NA,
                      title=NA,
                      pch.star = 17){

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

  if(x.axis=="Rank"){
    if(is.na(legend.where)){
      legend.where="topleft"
    }

    if(length(x.lim)!=2){
      x.lim = c(1,n)
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

    if(is.na(title)){
      title= paste0(p.fdr.object$`Adjustment Method`,
                    " p.fdr Object Plot")
    }
    if(adj.pvalues){
      plot(rank(p.fdr.object$`Results Matrix`[,3],ties.method = "random"),
           p.fdr.object$`Results Matrix`[,2],
           col="dodgerblue2",
           ylim=y.lim,
           xlim=x.lim,
           xlab="Ranking of Raw p-values",
           ylab=" ",
           main=title,
           pch=pch.star,
           cex.main=1.2,
           axes = FALSE,
           cex.axis=0.5)
      if(adj.sig.line){
        axis(side = 1)
        axis(side = 2,
             at = c(threshold),
             las=1,
             col.axis = "dodgerblue2",
             cex.axis=0.75)
        axis(side = 2,
             cex.axis=0.75,
             las=1)
        points(rank(p.fdr.object$`Results Matrix`[,3],ties.method = "random"),
               p.fdr.object$`Results Matrix`[,1],
               col="firebrick2",
               pch=20)
        abline(h=threshold,
               col="dodgerblue2",
               lwd=2)
      }else{
        axis(side = 1)
        axis(side = 2,
             las=1)
        points(rank(p.fdr.object$`Results Matrix`[,3],ties.method = "random"),
               p.fdr.object$`Results Matrix`[,1],
               col="firebrick2",
               pch=20)
      }
    }else{
      plot(rank(p.fdr.object$`Results Matrix`[,3],ties.method = "random"),
           p.fdr.object$`Results Matrix`[,1],
           col="firebrick2",
           ylim=y.lim,
           xlim=x.lim,
           xlab="Ranking of Raw p-values",
           ylab=" ",
           main=title,
           pch=20,
           las=1)
    }
    if(raw.pvalues){
      points(rank(p.fdr.object$`Results Matrix`[,3],ties.method = "random"),
             p.fdr.object$`Results Matrix`[,3],
             col="black",
             pch=20)
      if(sig.line){
        lines(x.line, y.line,
              col="black",
              lwd=2)
      }
    }

    #Legends
    if(raw.pvalues&adj.pvalues&sig.line&adj.sig.line){
      legend(legend.where,
             legend=c(paste0("Adj p-values (",p.fdr.object$`Adjustment Method`,")"),
                      paste0("FDRs (",p.fdr.object$`Adjustment Method`,")"),
                      "Raw p-values",
                      paste0(p.fdr.object$`Adjustment Method`," Threshold (Raw p-values)"),
                      paste0(p.fdr.object$`Adjustment Method`," Threshold (Adj p-values)")),
             pch=c(pch.star,20,20, NA, NA),
             lty=c(NA,NA,NA,1,1),
             lwd=c(NA,NA,NA,2,2),
             col=c("dodgerblue2","firebrick2", "black", "black","dodgerblue2"),
             cex=0.7,
             bty = "n")
    }
    else if(raw.pvalues&!adj.pvalues&sig.line&!adj.sig.line){
      legend(legend.where,
             legend=c(paste0("FDRs (", p.fdr.object$`Adjustment Method`,")"),
                      "Raw p-values",
                      paste0(p.fdr.object$`Adjustment Method`," Threshold (Raw p-values)")),
             pch=c(20,20, NA),
             lty=c(NA,NA,1),
             lwd=c(NA,NA,2),
             col=c("firebrick2", "black", "black"),
             cex=0.7,
             bty = "n")
    }
    else if(!raw.pvalues&adj.pvalues&!sig.line&adj.sig.line){
      legend(legend.where,
             legend=c(paste0("Adj p-values (",p.fdr.object$`Adjustment Method`,")"),
                      paste0("FDRs (", p.fdr.object$`Adjustment Method`,")"),
                      paste0(p.fdr.object$`Adjustment Method`," Threshold (Adj p-values)")),
             pch=c(pch.star,20, NA),
             lty=c(NA,NA,1),
             lwd=c(NA,NA,2),
             col=c("dodgerblue2","firebrick2", "dodgerblue2"),
             cex=0.7,
             bty = "n")
    }
    else{
      legend(legend.where,
             legend=c(paste0("FDRs (",p.fdr.object$`Adjustment Method`,")")),
             pch=c(20),
             col=c("firebrick2"),
             cex=0.7,
             bty = "n")

    }

  }else if(x.axis=="Zvalues"){
      if(is.na(legend.where)){
        legend.where="topright"
      }

      if(is.na(title)){
        title= paste0(p.fdr.object$`Adjustment Method`,
                      " p.fdr Z-values Plot")
      }

      if(adj.pvalues){
        plot(zvalues,
             p.fdr.object$`Results Matrix`[,2],
             col="dodgerblue2",
             ylim=y.lim,
             xlab="Z-values",
             ylab=" ",
             main=title,
             pch=pch.star,
             cex.main=1.2,
             axes=FALSE,
             cex.axis=0.5,
             las=1)
        axis(side = 1)
        axis(side = 2,
             at = c(threshold),
             las=1,
             col.axis="dodgerblue2",
             cex.axis=0.75)
        axis(side = 2,
             cex.axis=0.75,
             las=1)
        points(zvalues,
               p.fdr.object$`Results Matrix`[,1],
               col="firebrick2",
               pch=20)
        abline(h=threshold,
               col="dodgerblue2",
               lwd=2)
      }else{
        plot(zvalues,
             p.fdr.object$`Results Matrix`[,1],
             col="firebrick2",
             ylim=y.lim,
             xlab="Z-values",
             ylab=" ",
             main=title,
             pch=20,
             las=1)
      }

      if(raw.pvalues){
        points(zvalues,
               p.fdr.object$`Results Matrix`[,3],
               col="black",
               pch=20)
      }

    #Legends
    if(raw.pvalues&adj.pvalues&sig.line&adj.sig.line){
      legend(legend.where,
             legend=c(paste0("Adj p-values (",p.fdr.object$`Adjustment Method`,")"),
                      paste0("FDRs (",p.fdr.object$`Adjustment Method`,")"),
                      "Raw p-values",
                      paste0(p.fdr.object$`Adjustment Method`," Threshold (Adj p-values)")),
             pch=c(pch.star,20,20, NA),
             lty=c(NA,NA,NA,1),
             lwd=c(NA,NA,NA,2),
             col=c("dodgerblue2","firebrick2", "black","dodgerblue2"),
             cex=0.7,
             bty = "n")
    }
    else if(raw.pvalues&!adj.pvalues&sig.line&!adj.sig.line){
      legend(legend.where,
             legend=c(paste0("FDRs (", p.fdr.object$`Adjustment Method`,")"),
                      "Raw p-values"),
             pch=c(20,20),
             col=c("firebrick2", "black"),
             cex=0.7,
             bty = "n")
    }
    else if(!raw.pvalues&adj.pvalues&!sig.line&adj.sig.line){
      legend(legend.where,
             legend=c(paste0("Adj p-values (",p.fdr.object$`Adjustment Method`,")"),
                      paste0("FDRs (", p.fdr.object$`Adjustment Method`,")"),
                      paste0(p.fdr.object$`Adjustment Method`," Threshold (Raw p-values)")),
             pch=c(pch.star,20, NA),
             lty=c(NA,NA,1),
             lwd=c(NA,NA,2),
             col=c("dodgerblue2","firebrick2", "dodgerblue2"),
             cex=0.7,
             bty = "n")
    }
    else{
      legend(legend.where,
             legend=c(paste0("FDRs (",p.fdr.object$`Adjustment Method`,")")),
             pch=c(20),
             col=c("firebrick2"),
             cex=0.7,
             bty = "n")
    }
  }

}



###
##
#
