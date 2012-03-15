\name{qExtrapolation}
\alias{extract.qdata}
\alias{qextrap.test}
\alias{predqPCA}
\alias{synthPCA.qp}
\alias{qpextrapPCA}
\alias{qextrap.eval}
\alias{fqextr}
\alias{annual.max}
\alias{mu2q}
\alias{mu2extr}
\alias{clim2qextr}
\alias{qrl}
\title{Estimating very high percentiles through regression analysis and
  an extrapolation, and an evaluation of this extrapolation.}
\description{
  One question is whether it is possible to extrapolate the shapes of
  the qPCA EOFs to describe more extreme percentiles. This set of
  functions are written to answer this question. The idea is to
  subsample the GDCN data set and keep the station with sufficiently
  long series to provide good estimates of these high percentiles. Then
  the data is split into two samples: one consisting of percentiles
  between 0.50 and 0.99, and another for percentiles greater than
  0.99. \code{\link{qPCA}} is performed on the low-percentile sample,
  and a regression is used to extend the shape of the two leading
  EOFs. These extended shapes are then compared with correpsonding EOFs
  estimated from the complete sample. The functions also predict the
  higher percentiles, which are then compared with those from the data.
  
  \code{extract.qdata} extract the longest series in the GDCN data set
  in order to estimate a number of percentiles at levels in p=[0.99 to
  0.9995].

  \code{qextrap.test} is used to see if a qPCA can be extrapolated
  through the means of regression analysis. A simple and smooth shape in
  the EOFs gives hope.
  
  \code{qextrap.test} shows the extrapolated qPCA shapes and compare
  with those estimated from the complete sample.
  
  \code{synthPCA.qp} generate 'synthetic' qPCA object by extrapolating
  their shapes towards higher percentiles. Involves regression analysis.

  \code{qextrap.eval} evaluates the extrapolated results, estimating
  errors.

  \code{fqextr} estimates precentiles and return value statistics in a
  similar fashion as \code{\link[evd]{fgev}}.
}

\usage{
extract.qdata(do.recompute=FALSE,
              disk2.path="GDCN-disk2/",plot=FALSE)
qextrap.test(mu.qp.all,pdiv=0.99)
predqPCA(pca,pca.all)
synthPCA.qp(mu.qp.all,pdiv=0.99,plot=TRUE)
qpextrapPCA(pcaextrap,cpx="max(probs)",neofs=2,plot=TRUE)
qextrap.eval(results,plot=TRUE)
annual.max(x,n.min=300)
fqextr(x=NULL,mu.qp=NULL,x.0=1,neofs=2,plot=TRUE)
mu2q(x=NULL,plot=TRUE,x.0=1)
mu2extr(mu,neofs=2,pca=NULL,plot=TRUE)
clim2qextr(x,p=0.95,x.0=1,modelq=NULL)
qrl(x=NULL,x.0=1)
}

\arguments{
  \item{do.recompute}{TRUE: do the original extraction from the GCDN.}
  \item{disk2.path}{path to where disc2 data are stored}
  \item{mu.qp.all}{qqplotter object - all the data (all percentiles)}
  \item{pdiv}{probability threshold for dividing samples of low and
  hifgh percentiles}
  \item{pcaextrap}{qPCA-object extrapolated to higher percentiles.}
  \item{cpx}{The highest percentile}
  \item{neofs}{number of EOFs}
  \item{results}{the output from \code{qpextrapPCA}}
  \item{plot}{TRUE: produce graphics}
  \item{pca}{qPCA-object}
  \item{pca.all}{qPCA-object for all percentiles.}
  \item{x}{time series - vector with an attribute "year" for time stamp.}
  \item{x.0}{threshold for wet-day}
  \item{n.min}{smallest number of valid data in a year for acceptable
  maximum estimate.}
  \item{mu.qp}{qqplotter object}
  \item{mu}{mean wet-day value}
  \item{p}{probability}
  \item{modelq}{model derived from \code{\link{modelqPCs}}}
}

\value{List objects containing matrices.
}
\author{R.E. Benestad}
\examples{
\dontrun{
extract.qdata() -> mu.qp.all
synthPCA.qp(mu.qp.all) -> pcaextrap
qpextrapPCA(pcaextrap) -> results
qextrap.eval(results)
}
}
\keyword{models}
\keyword{manip}
\keyword{methods}
