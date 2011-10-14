\name{makemap}
\alias{makemap}
\alias{modelqPCs}
\alias{qqatribution}
\alias{qqscenario}
\alias{qpAnywhere}
\alias{rqstats}
\alias{exploreqmodel}
\title{Mapping qqplotter results for North America.}
\description{
  The function \code{makemaps} uses the kriging methods from the
  \code{LatticeKrig} (\url{http://www.image.ucar.edu/~nychka/}) package
  to produce maps. The default is for North America, but the function
  can also be applied to other regions.

  The function \code{modelqPCs} seeks a relationship between the two
  leading qPCs from \code{\link{qPCA}} and the mean total precipitation,
  elevation and distance to the coast.

  \code{qpAnywhere} is a quick way of extracting percentiles for the
  station closest to a given coordinate.

  \code{qqatribution} applies the qPC-model to the time dimensions of
  the station with the longest data record, splitting the time series
  into 5-yr segments, and use wet-day mean and 
  frequency to predict 95th-percentils. \code{rqstats} applies the
  same temporal evaluation to all stations in \code{mu.qp} and returns a
  distribution for the correlation coefficients.

  \code{exploreqmodel} compares the effect of varying one of the
  covariates (wet-day mean, wet-day frequency, elevation, and distance
  from the coast) on the prediction of the 95th percentile.

  Some functions can take while, but they will save temporary files
  underways, so that jobs can be interrupted and stard from near the
  point where they stopped. When finished, the results taking time to
  compute will be stored in local files with same name as the function
  so that a subsequent call will run faster.

  \code{qqscenario} is merely a demonstration of how the qqplotter frame
  work can be used as a means of downscaling future projections for the
  95th percentile. These results should not be interpreted as real
  projections; The two time slices from the NCAR CSM3.0 model was
  deliberately chosen from different runs to avoid the use of  these
  results as projections for the future - for more information see
  \code{\link[PrecipStat]{mu.1990.1999}}  
}

\usage{
makemap(mu.qp=NULL,statistic="ratio",keep.longest=TRUE,
             lon.rng=c(-180,-50),lat.rng=c(25,70),stnr.given=TRUE,
            nx=200, ny=100, NC=64,zlim=NULL)
modelqPCs(mu.qp=NULL,keep.longest=TRUE,
             lon.rng=c(-180,-50),lat.rng=c(25,70),
                      covariates=c("mu","altitude","dist2coast","wet.freq"),
             Eval.Region=c("NAmerica","WestEurope","EastEurope","Brazil",
                           "S.E.Asia"))
qqscenario(result=NULL,type=c("dq95dw/q95dw","q95dw.sce"))
qqatribution(mu.qp=NULL,result=NULL,
             path="~/GDCN/",x.0=1,months=NULL)
qpAnywhere(lon,lat,mu.qp=NULL,p=0.95)
rqstats(mu.qp=NULL,result=NULL,
                    path="~/GDCN/",x.0=1,months=NULL)
exploreqmodel(mu.qp=NULL,result=NULL)
}

\arguments{
\item{mu.qp}{A qqplotter object: \code{\link[PrecipStat]{mu.qp}}}  
\item{statistic}{The statistic to map: if can either be a character
  string describing a stastistic: c("ratio","mu","mean_precip",
           "n.wet","wetfreq","altitude","error",
           "sigma","q95","PC1","PC2","q95.pred","dist2coast","rmse"), or
	   a vector of numerics with attributes c("name","longitude","latitude","station_number").}
\item{keep.longest}{TRUE: for stations with same coordinate, keep the
  one with longest record.}
\item{lon.rng}{Longitude range of region of interest}
\item{lat.rng}{Latitude  range of region of interest}
\item{stnr.given}{FALSE will accept station numbers set to NA}
\item{nx}{parameter passed on to \code{\link[fields]{fields.x.to.grid}}}
\item{ny}{parameter passed on to \code{\link[fields]{fields.x.to.grid}}}
\item{NC}{parameter passed on to \code{\link[fields]{fields.x.to.grid}}}
\item{zlim}{parameter passed on to \code{\link[fields]{image.plot}}}
\item{covariates}{A set of co-variates for the regression analysis of
  the PCs from
  \code{\link{qPCA}}-call. c("mP","mu","n.wet","n.dry","rmse",
  "sigma","country","wet.freq","altitude","dist2coast","sin.lat","cos.lat",
  "sin.lon","cos.lon")}.
\item{Eval.Region}{a vector of character strings for regions for which
  to make maps of predicted q95. Fewer numer of regions gives faster analysis.}
\item{result}{The results from \code{modelqPCs}.}
\item{path}{Path to the GDCN data}
\item{x.0}{Threshold deciding wet-day}
\item{months}{Months to include in the analysis}
\item{lon}{longitude}
\item{lat}{latitude}
\item{p}{(100*p)-percentile}
\item{type}{Type of map: The codes distinquish between wet-day
  quantities (e.g. \code{q95} and \code{dq95}) and dry+wet day
  quantities (e.g. \code{q95dw} and \code{dq95dw}). The subscript 'd'
  means change, and suffix '.sce' scenario. Maps can be made for ratios
  and differences. The possible choices are:
  c("dq95","q95.sce","dq95/q95","q95dw","q95dw.sce",
  "dq95dw","dq95dw/q95dw","q95dw.sce/q95dw")}.
}

\value{The output from \code{\link[fields]{predict.surface}} in \code{\link{fields}}.
}
\author{D. Nychka and modified by R.E. Benestad}
\examples{
\dontrun{

# Map the various statistics
statistics <-  c("ratio","mu","mean_precip","n.wet","wetfreq","altitude",
                 "sigma","q95","PC1","PC2","q95.pred","dist2coast","error",
                 "rmse")
for (i in 1:length(statistics)) {
  makemap(statistic=statistics[i]) -> output
  dev.copy2eps(file=paste("makemap_",statistic,".eps",sep=""))
}

# Predict PCs1-2 and the 95th percentile
modelqPCs() 
}
}
\keyword{models}
\keyword{manip}
\keyword{methods}
