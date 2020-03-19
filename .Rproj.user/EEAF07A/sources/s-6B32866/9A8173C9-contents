\name{p.fdr}
\alias{p.fdr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{FDR Computation}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
p.fdr(pvalues, zvalues = NA, threshold = 0.05, adjust.method = NA, just.fdr = FALSE, pi0.estim = "one", set.pi0 = NA, na.rm = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pvalues}{
%%     ~~Describe \code{pvalues} here~~
}
  \item{zvalues}{
%%     ~~Describe \code{zvalues} here~~
}
  \item{threshold}{
%%     ~~Describe \code{threshold} here~~
}
  \item{adjust.method}{
%%     ~~Describe \code{adjust.method} here~~
}
  \item{just.fdr}{
%%     ~~Describe \code{just.fdr} here~~
}
  \item{pi0.estim}{
%%     ~~Describe \code{pi0.estim} here~~
}
  \item{set.pi0}{
%%     ~~Describe \code{set.pi0} here~~
}
  \item{na.rm}{
%%     ~~Describe \code{na.rm} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (pvalues, zvalues = NA, threshold = 0.05, adjust.method = NA,
    just.fdr = FALSE, pi0.estim = "one", set.pi0 = NA, na.rm = TRUE)
{
    library(splines)
    n = length(pvalues)
    if (na.rm) {
        pvalues = pvalues[!is.na(pvalues)]
        n = length(pvalues)
    }
    if (null.estim == "one") {
        pi0 = 1
    }
    else if (null.estim == "last.hist") {
        try.hist <- hist(pvalues, breaks = min(n/10, 100), plot = FALSE)
        try.mids <- try.hist$mids
        try.count <- try.hist$counts
        pi0 = tail(try.count, 1) * length(try.mids)/sum(try.count)
    }
    else if (null.estim == "loess") {
        try.hist <- hist(pvalues, breaks = n/10, plot = FALSE)
        try.mids <- try.hist$mids
        try.count <- try.hist$counts
        loe <- loess(try.count ~ try.mids)
        k <- 0.3
        x.loe <- try.mids[try.mids >= k]
        y.loe <- loe$fitted[try.mids >= k]
        avg.loe <- mean(y.loe)
        y.est <- c(rep(avg.loe, (length(try.mids) - length(y.loe))),
            y.loe)
        pi0 = sum(y.est)/sum(try.count)
    }
    else if (null.estim == "lindseys") {
        if (sum(is.na(zvalues)) == 1) {
            stop("No argument for `zvalues` provided")
        }
        sim.hist <- hist(zvalues, plot = FALSE, breaks = min(n/10,
            100))
        bin <- length(sim.hist$mids)
        df <- 7
        breaks <- sim.hist$breaks
        y <- sim.hist$counts
        x <- sim.hist$mids
        K <- length(y)
        k <- seq(K)
        yhat <- glm(y ~ ns(x, df = df), poisson)$fit
        bw <- x[2] - x[1]
        fhat <- yhat/(n * bw)
        fdra <- dnorm(x, 0, 1)/fhat
        plot(x, pmin(fdra, 1), type = "l", col = "blue", main = "Lindsey's adjust.method FDR vs. z-values")
        abline(h = 0.05, col = "red", lty = 2)
        pi0 = min(quantile(fhat/dnorm(x, 0, 1), 0.05), 1)
    }
    fdr.bh = pmin(1, pi0 * pvalues * n/rank(pvalues))
    if (is.na(adjust.method) | adjust.method == "BH") {
        adj.pvalues = cummin(pmin(1, pvalues * n/rank(pvalues)))
        adj.fdrs = pmin(1, pi0 * pvalues * n/rank(pvalues))
    }
    else if (adjust.method == "BY") {
        dep = cumsum(1/(1:n))
        adj.pvalues = cummin(pmin(1, dep * pvalues * n/rank(pvalues)))
        adj.fdrs = pmin(1, pi0 * dep * pvalues * n/rank(pvalues))
    }
    else if (adjust.method == "Bon") {
        adj.pvalues = pmin(1, n * pvalues)
        adj.fdrs = pmin(1, pi0 * n * pvalues)
    }
    else if (adjust.method == "Holm") {
        adj.pvalues = cummax(pmin(1, (n - rank(pvalues) + 1) *
            pvalues))
        adj.fdrs = pmin(1, (n - rank(pvalues) + 1) * pvalues)
    }
    out <- as.data.frame(cbind(`Method Adjusted FDRs` = adj.fdrs,
        `Method Adjusted P-Values` = adj.pvalues, `Unadjusted P-Values` = pvalues,
        `Individual BH FDRs` = fdr.bh))
    return(list(fdrs = adj.fdrs, `Results Matrix` = out, `Reject Vector` = ifelse(adj.pvalues <=
        threshold, "Reject.H0", "FTR.H0"), pi0 = pi0, threshold = threshold))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
