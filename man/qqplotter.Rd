\name{qqplotter}
\alias{qqplotter}
\alias{qplotGDCN}
\alias{qplotESCN}
\alias{qplotRCM}
\alias{qplotres}
\alias{niceqqplot}
\alias{compare1q}
\title{Quantile-quantile plots for 24-hr precip assuming exponential distribution}
\description{
   
\code{qqplotter} - carries out diverse jobs, such as box qq-plots, calls
qplorGDCN, \code{\link{qPCA}}, niceqqplot, etc. This function generates
a number of plots and statistics, as well as performing various analysis
in terms of PCA, etc.
   
\code{qplotGDCN} - the function producing qq-plots for the GDCN data. It
assumes that these data are located in a directory (or link to one)
'~/GDCN/' (in the home directory), and that these data are unzipped ASCII-files with the suffix
'.dly'.

\code{qplotESCN} - same as \code{qplotGDCN} but analysing the ECA\&D
data for Europa as opposed to the U.S. GDCN data. The function assumes
that the data are stored in a directory (or link to one) called
'data.eca/' in the local path. It assumes that the data are unzipped
ASCII-files with file names starting with 'RR'.

\code{qplotRCM} - plots corresponding qq-plots for RCMs, contrasting
RCMs from the NARCCAP project with that of the ENSEMBLES project. The
NARCCAP data are assumed to be in R-data-format, whereas the ENSEMBLES
data are assumed to be in (the original) netCDF format (\url{http://ensemblesrt3.dmi.dk/}). 

\code{qplotres} - similar as \code{qplotRCM}, but contrasts qq-plots for
HadRM3 with different spatial resolutions.


\code{niceqqplot} - produces a colourful qq-plot from pre-calculated
results.

  Some functions can take while, but they will save temporary files
  underways, so that jobs can be interrupted and stard from near the
  point where they stopped. When finished, the results taking time to
  compute will be stored in local files with same name as the function
  so that a subsequent call will run faster.
}

\usage{
qplotGDCN(x.0=1,Pmax=180,prs=c(seq(0.50,0.95,by=0.05),0.96,0.97,0.98,0.99),
          max.stations=NULL,xlim=c(0,200),ylim=c(0,200),
          months=NULL,years=NULL,N.min=NULL,path="~/GDCN/")
qplotESCN(x.0=1,Pmax=180,prs=c(seq(0.50,0.95,by=0.05),0.96,0.97,0.98,0.99),
          max.stations=NULL,xlim=c(0,200),ylim=c(0,200),
          months=NULL,years=NULL,N.min=NULL)
qplotRCM(mu.qp=NULL,x.0=1,sim.rng=TRUE,
                     prs=c(seq(0.50,0.95,by=0.05),0.96,0.97,0.98,0.99),
                     do.narccap=FALSE, do.ensembles=TRUE,Max.pts=153569,
                     xlim=c(0,200),ylim=c(0,200),months=NULL,
                     path.narccap="~/NARCCAP/data/",
                     path.ensembles="~/ENSEMBLES.daily.50km")
qplotres(mu.qp=NULL,x.0=1,sim.rng=TRUE,
                     prs=c(seq(0.50,0.95,by=0.05),0.96,0.97,0.98,0.99),
                     xlim=c(0,100),ylim=c(0,100),
                     Max.pts=90000,path.25km="~/ENSEMBLES.daily",
                     path.50km="~/ENSEMBLES.daily.50km")
niceqqplot(mu.qp,col.axis="black",addcont=TRUE,addfit=TRUE,
           xlim=c(0,200),ylim=c(0,200),colleg=TRUE,
           x.0=1,colourcoding="mean_precip",
           CI=TRUE,N.min=NULL,p=NULL,add=FALSE,
           col=NULL,pch=19,cex=1)

qqplotter(x.0=1,Pmax=180,do.gdcn=FALSE,do.escn=FALSE,
          do.rcm=FALSE,do.res=FALSE,do.narccap=TRUE,
          do.ensembles=TRUE,max.stations=NULL,europe2=TRUE,
          prs=c(seq(0.50,0.95,by=0.05),0.96,0.97,0.98,0.99))
compare1q(mu.qp1,mu.qp2,p=0.95,N.min=500,type="rcm",
          xlim=c(0,100),ylim=c(0,100),
          col=c("blue","red","darkgreen","steelblue",
                "darkred","green","magenta","cyan","grey30",
                "lightblue","pink","grey80","wheat","brown"))

}

\arguments{
\item{x.0}{Threshold value for defining a wet-day.}
\item{prs}{Probability levels for defining the set of quantiles.}
\item{mu.qp}{A \code{qqplotter} (matrix) object holding the quantiles. See \code{link{mu.qp}}.}
\item{Pmax}{Expected maximum (wet+dry) mean precipitation. Used for
  colour coding the points in the qqplots.}
\item{max.stations}{The maximum number of stations to plot (used for
  quick testing).}
\item{sim.rng}{TRUE: compare the RCM results against observed
  precipitation for which the quantiles have similar range.}
\item{do.narccap}{Unfinished - should be FALSE.}
\item{do.ensembles}{TRUE: analyse the ENSEMBLES RCMs, other wise skip ENSEMBLES.}
\item{Max.pts}{Maximum number of RCM grid-points to include in the analysis.}
\item{path}{Path pointing to the GDCN results.}
\item{path.narccap}{Path pointing to the NARCCAP RCM results.}
\item{path.ensembles}{Path pointing to the ENSEMBLES RCM results.}
\item{path.25km}{Path pointing to high-resolution ENSEMBLES RCM results.}
\item{path.50km}{Path pointing to low-resolution ENSEMBLES RCM results.}
\item{do.gdcn}{TRUE: carry out the analysis for GDCN from scratch,
  otherwise use pre-calculated results.}
\item{do.escn}{TRUE: carry out the analysis for ESCN from scratch,
  otherwise use pre-calculated results.}
\item{do.rcm}{TRUE: carry out the analysis for RCMs from scratch,
  otherwise use pre-calculated results.}
\item{do.res}{TRUE: carry out the analysis for RCMs with different
  resolution from scratch, otherwise use pre-calculated results.}
\item{col.axis}{Colour of axis. "White" for fancy presentation in PowerPoint.}
\item{addcont}{TRUE: add contours to cloud of points.}
\item{addfit}{TRUE: add a bestfit by lm(x2 ~ x1).}
\item{colleg}{TRUE: add a colourbar legend.}
\item{xlim}{see \code{\link{plot}}}
\item{ylim}{see \code{\link{plot}}}
\item{p}{Probability level}
\item{col}{Defines the colours for the colourbar legend.}
\item{cex}{Point size. See \code{\link{plot}}.}
\item{pch}{Point symbol. See \code{\link{plot}}.}
\item{type}{Type of criterion for distinguishing between RCMs/data sets.}
\item{colourcoding}{Criterion for colour coding - an expression}
\item{add}{TRUE - plot on top of old figure}
\item{mu.qp1}{a \code{qqplotter} object.}
\item{mu.qp2}{a \code{qqplotter} object.}
\item{months}{Months to include in the analysis.}
\item{years}{Years to include in the analysis.}
\item{N.min}{Minimum number of wet-days that is acceptable.}
\item{CI}{TRUE: plot confidence intervals.}
\item{europe2}{Use European stations too.}
}

\value{\code{qqplotter} objects, (matrix) object holding the quantiles. See \code{link{mu.qp}}.
}
\author{R.E. Benestad}
\examples{
\dontrun{
data(mu.qp)
names(attributes(mu.qp))
niceqqplot(mu.qp,N.min=1000)

mu.qp <- qplotGDCN()
prs.extreme <- 0.9 + 0.09*seq(0,1,length=20)^0.25
mu.qp.extr <- qplotGDCN(prs=prs.extreme,xlim=c(0,2000),ylim=c(0,2000))
attr(mu.qp.extr,"probabilities") <- prs.extreme
save(file="mu.qp.extr.rda",mu.qp.extr)

new.dev()
qqplotter()
}
}
\keyword{datagen}
\keyword{methods}
\keyword{hplot}
