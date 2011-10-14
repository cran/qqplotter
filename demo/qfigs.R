# Script for making the plot in the publication

# Retrieve daily rain data from the U.S., Europe, ++

require(PrecipStat)

data(mu.qp)
data(mu.qp.escn)

# Only include stations with more than 1000 valid wet-day measurements
mu.qp <- qweed(mu.qp)
mu.qp.escn <- qweed(mu.qp.escn)

MU.qp <- qcat(mu.qp,mu.qp.escn)
niceqqplot(MU.qp,p=0.95,addcont=FALSE,addfit=FALSE) -> a
dev.new()

niceqqplot(MU.qp,p=0.99,addcont=FALSE) -> b

dev.new()
mapofstations(MU.qp)

data(mu.qp.rcm)
dev.new()
compare1q(mu.qp.escn,mu.qp.rcm) 

dev.new()
qqplotter()




