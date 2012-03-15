\name{Data}
\alias{pca.qextr}
\alias{samnanger}
\alias{qerrors}
\title{Pre-calculated results for estimating very high percentiles.}
\description{
\code{pca.qextr} contains qPCA results where the shape has been
extrapolated to higher percentiles through the \code{\link{synthPCA.qp}}
function.

\code{qerrors} are error estimates from \code{\link{qextrap.eval}}.

\code{samnanger} is a sample data set: precipitation from Samnanger in
Norway. See the data attributes for more details.
}
\usage{
data(pca.qextr)
data(qerrors)
data(samnanger)
}


\value{
  A \code{\link{qPCA}} object.
}
\author{R.E. Benestad}
\examples{
# extract.qdata() -> mu.qp.all
# synthPCA.qp(mu.qp.all) -> pcaextrap
# pca.qextr <- pcaextrap$pca.syn
# attr(pca.qextr,'pcafit1') <- pcaextrap$pcafit1
# attr(pca.qextr,'pcafit2') <- pcaextrap$pcafit2
# attr(pca.qextr,'pdiv') <- pcaextrrap$pdiv

data(pca.qextr)
}
\keyword{datasets}
