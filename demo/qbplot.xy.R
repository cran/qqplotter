require(PrecipStat)

data(mu.qp)
dev.new()
mu.qp$x1[mu.qp$x1>1000] <- NA
mu.qp$x2[mu.qp$x2>1000] <- NA
bplot.xy(x=c(mu.qp$x1), y=c(mu.qp$x2), N = 10,
           main="Wet-day 24-hr precipitation from GDCN",
           xlab="-log(1-p)*mu [mm/day]",ylab="quantile(X,p) [mm/day]",lwd=3,
           sub=paste("thresh.= ",x.0,"mm/day;",
           " [qqplotter.R]"))
pd <- pointdensity(mu.qp$x1,mu.qp$x2)
contour(pd$x,pd$y,pd$z,add=TRUE,col="grey",lev=seq(25,20000,by=250))
lines(c(0,900),c(0,900),col="grey",lty=2)
grid()
