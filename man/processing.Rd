\name{estquantiles}
\alias{estquantiles}
\alias{fitsimple}
\alias{qPCA}
\alias{qPCA2quantile}
\alias{expMC}
\alias{qplotMC}
\alias{mapofstations}
\alias{qcat}
\alias{qsubset}
\alias{qweed}
\alias{cleanduplicates}
\alias{dist2coast}
\title{Calculations  and processing of quantiles.}
\description{
  The set of functions are used to analyse precipitation using an
  extremely quick and simple approach. For an exponential distribution,
  the mean (mu) specifies any quantile (q_p) according to \deqn{q_p = -ln(1-p)
    \mu}{q_p = -ln(1-p)*mu}. Using this equation, the quantiles of 24-hr
  accumulated precipitation records are plotted against the empirical
  value 'quantile(X,p)', where p is the probability (\deqn{p \in
    [0,1]}{0 .le. p .l.e 1}).
  
  The motivation for this analysis goes back to a couple of papers
  (Benestad, 2007, 2010) which suggested that the 24-hr precipitation
  was close to exponential and that one single parameter - the slope
  of the curve m - could be used to specify any quantile.

  These functions used an improved description of the quantiles, basing
  the estimates on the mean (mu) rather than the slope (m) that is less
  well estimated.
  
   [ref: Benestad, R.E. (2007), 'Novel Methods for Inferring Future Changes in
   Extreme Rainfall over Northern Europe' , \emph{Climate Research},
   \bold{CR34}:195-210, doi: 10.3354/cr00693,
   \url{http://www.google.com/url?q=http\%3A\%2F\%2Fwww.int-res.com\%2Farticles\%2Fcr_oa\%2Fc034p195.pdf&sa=D&sntz=1&usg=AFQjCNEwNhTO47jhNOqQJ9OgoneQb7KKFA};
   
     Benestad, R.E. (2010) 'Downscaling Precipitation Extremes:
   Correction of Analog Models through PDF Predictions', \emph{Theor. \&
   Appl. Clim.}, \bold{Volume 100}, Issue 1, DOI:
 10.1007/s00704-009-0158-1. \url{http://www.google.com/url?q=http\%3A\%2F\%2Fwww.springerlink.com\%2Fcontent\%2Fy65j54125146p95g\%2F&sa=D&sntz=1&usg=AFQjCNG1zv8FinCFtwlE3rtvF0d_bm_FfQ};

Benestad, R.E., D. Nychka and L.O. Mearns, 'Specification of wet-day daily rainfall quantiles from the mean value', Tellus A, 64, 14981, DOI: 10.3402/tellusa.v64i0.14981 \url{http://www.tellusa.net/index.php/tellusa/article/view/14981}
 ]

   The functions in this package are as follows:
   

\code{expMC} and \code{qplotMC} - functions to perform Monte-Carlo
simulations based on the function \code{rexp()}. Used for estimating
confidence intervals.

\code{mapofstations} - function to show the geographical distribution of
the stations, plotted as colur-coded symbols. The critereon for colour
coding is flexiple; quantiles, mean precip, or scenarios for the future.

\code{qcat} and \code{qweed} - functions to post-process the data: one
combines two results and the other is used for weeding out stations,
e.g. with a low number of wet days. \code{cleanduplicates} removes
duplicated locations, keeping the stations with the highest number of
wet days (data). \code{qsubset} extracts a subset of the percentiles.

\code{estquantiles} and \code{fitsimple} - the main functions for
computing quantiles, setting the variable \code{x1} to hold \deqn{q_p =
  -ln(1-p)\mu}{q_p = -ln(1-p)*mu} and and \code{x2} to hold
\code{quantile(X,p)}, with e.g. \code{p=seq(0.5,0.95,by=0.05)} as default. 

\code{qPCA} - performs a principal component analysis (actually SVD) on the
quantiles. For each data record, a vector x is constructed by merging
the quantiles held in x1 with x2 (from
\code{fitsimple}). \code{qPCA2quantile} reconstructs the percentiles,
given a \code{qPCA} object and weights for the principal components
(PCs). See \code{link{testqPCA2quantile}}.

\code{dis2coast} estimates the distance to the coast in km (this
information is added as an attribute 'dist2coast.km'.
}

\usage{
estquantiles(mu.qp,
             prs=c(seq(0.50,0.95,by=0.05),0.96,0.97,0.98,0.99))
fitsimple(X,x.0=1,prs=c(seq(0.50,0.95,by=0.05),0.96,0.97,0.98,0.99))
expMC(mu,p,N,N.runs=1000)
qplotMC(mu=10,N=90,prs=c(seq(0.50,0.95,by=0.05),0.96,0.97,0.98,0.99),
        N.runs=1000,plot=TRUE)
qPCA(mu.qp,plot=TRUE)
qPCA2quantile(x,pca,p=0.95,silent=FALSE)
mapofstations(mu.qp,colourcoding="0.95",add=FALSE,silent=FALSE,
              googleearth=FALSE,kmz.file="~/Desktop/qqplotter.kmz")
qcat(mu.qp1,mu.qp2)
qsubset(mu.qp,crit="attr(mu.qp,'probabilities') <= 0.99")
qweed(mu.qp,crit="attr(mu.qp,'n.wet')> 1000")
cleanduplicates(mu.qp,silent=TRUE,plot=FALSE)
dist2coast(mu.qp)
}

\arguments{
\item{X}{A vector containing the 24-hr accumulated precipitation amounts
(mm/day).}
\item{x.0}{Threshold value for defining a wet-day.}
\item{prs}{Probability levels for defining the set of quantiles.}
\item{mu.qp}{A \code{qqplotter} (matrix) object holding the quantiles. See \code{link{mu.qp}}.}
\item{mu}{Mean wet-day precipitation}
\item{p}{Probability level}
\item{N}{Number of data points (mimicking number of wet days)}
\item{N.runs}{Number of runs with Monte-Carlo simulations}
\item{plot}{TRUE: plot results}
\item{silent}{FALSE: output transcript}
\item{colourcoding}{Criterion for colour coding - an expression}
\item{add}{TRUE - plot on top of old figure}
\item{googleearth}{TRUE: produce a KML-file for GoogleEarth}
\item{kmz.file}{Name of KML-file}
\item{mu.qp1}{a \code{qqplotter} object. See \code{link{mu.qp}}}
\item{mu.qp2}{a \code{qqplotter} object. See \code{link{mu.qp}}}
\item{crit}{The critereon for weeding out data - an expression.}
\item{pca}{A  \code{qPCA} object - the output from a \code{qPCA} call.}
\item{x}{A matrix, vector or a single number describing PCs of \code{pca}.}
}

\value{List objects containing matrices \code{x1} and \code{x2} with
  dimensions  \code{[N.p,N.l]}, where \code{N.p} is the number of
quantiles for each record and \code{N.l} the number of locations.
}
\author{R.E. Benestad}
\examples{
require(PrecipStat)
data(mu.qp)
names(attributes(mu.qp))

mu.qp.wet <- qweed(mu.qp,crit="attr(mu.qp,'mean_precip')> 50")

mu.qpx <- qweed(mu.qp)
mapofstations(mu.qpx,colourcoding="rmse")

# PCA of the quantiles:
pca <- qPCA(mu.qp)
x <- pca$v[1:1000,1:2]
qPCA2quantile(x,pca,p=0.95)

testqPCA2quantile()
}
\keyword{models}
\keyword{manip}
\keyword{methods}
