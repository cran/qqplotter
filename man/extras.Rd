\name{colorbar}
\alias{colorbar}
\alias{pointdensity}
\alias{Sum.is.zero} 
\alias{Sum.gt.zero}
\alias{testqPCA2quantile}
\alias{Points}
\alias{Lines}
\alias{replace.char}
\title{Supporting functions for qqplotter.}
\description{

\code{colorbar} - a function to make a colour bare legend for some of
the plots.

\code{pointdensity} - a function that produces a density surface of the
a set of points. 
}

\usage{
Points(x,y)
Lines(t,y,col="black",
      tshadeoffset=0.03,yshadeoffset=0.03)
colorbar(ckey,fig=c(0.33,0.66,0.30,0.34))
pointdensity(x1,x2,resol=10)
Sum.is.zero(X) 
Sum.gt.zero(X)
testqPCA2quantile(mu.qp=NULL,subset=1:1000,n.eofs=2)
replace.char(c, s, ny.c)
}

\arguments{
\item{X}{A vector containing the 24-hr accumulated precipitation amounts
(mm/day).}
\item{ckey}{Colour key for colourbar legend.}
\item{fig}{Defines the size and location of the colourbar legend.}
\item{x1}{vector holding quantiles estimated from mu.}
\item{x2}{vector holding quantiles estimated using 'quantile(X,p)'.}
\item{resol}{resolution of the grid boxes used for estimating the
  density of the points.}
\item{mu.qp}{\code{qqplotter} (matrix) object holding the quantiles. See \code{link{mu.qp}}}
\item{subset}{Subset of locations to test}
\item{n.eofs}{Number of EOFs to include in the recnstruction}
\item{t}{time coordinates}
\item{x}{x-coordinates}
\item{y}{y-coordinates}
\item{col}{col as in \code{\link{plot}}.}
\item{tshadeoffset}{for shade effect - x-axis}
\item{yshadeoffset}{for shade effect - y-axis}
\item{c}{Character to be replaced}
\item{s}{String in which one character is to be replaced}
\item{ny.c}{The new character substituting the old one.}
}

\value{
}
\author{R.E. Benestad}
\examples{

}
\keyword{hplot}
\keyword{manip}
\keyword{misc}
