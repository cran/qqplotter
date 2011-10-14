require(PrecipStat)

data(mu.qp.world)
mapofstations(mu.qp.world)

pca <- qPCA(mu.qp.world)

dev.new()
plot(c(0,200),c(0,200),type="l",col="grey",
      main="Variations in PCA modes #1 & 2",
      xlab="-log(1-p)*mu [mm/day]",ylab="quantile(X,p) [mm/day]",lwd=3,
            sub=paste("thresh.= ",attr(mu.qp.world,"x.0"),"mm/day;",
            " [qqplotter.R]"))

N <- length(pca$v[,1])
for (i in 1:N) lines(pca$x1[,1]*pca$d[1]*pca$v[i,1]+
                     pca$x1[,2]*pca$d[2]*pca$v[i,2],
                     pca$qp[,1]*pca$d[1]*pca$v[i,1]+
                     pca$qp[,2]*pca$d[2]*pca$v[i,2],col="blue")
for (i in 1:N) lines(pca$x1[,1]*pca$d[1]*pca$v[i,1],
                     pca$qp[,1]*pca$d[1]*pca$v[i,1],col="red")

lines(c(0,200),c(0,200),col="grey")
legend(150,30,c("Data","mode 1","mode 1+2"),
        col=c("grey70","red","blue"),pch=c(19,26,26),
        lty=c(0,1,1),bg="grey95")

testqPCA2quantile()

# Carry out the test with all EOFs and for the whole set:
testqPCA2quantile(mu.qp.world,subset=NULL,n.eofs=NULL) -> Scatter

# Analyse the test results:

# Contours of the scatter plot:
xydensity <- table(3*round(Scatter$X/3),3*round(Scatter$Y/3))
contour(as.numeric(rownames(xydensity)),as.numeric(colnames(xydensity)),
        xydensity,xlim=c(0,200),ylim=c(0,200),
        main="Density of scatter between reconstruction and actual q95",
        xlab="Original (mm/day)",ylab="Reconstructed (mm/day)")
grid()
lines(c(0,200),c(0,200),lty=2,col="grey")

# Examine the locations with largest scatter:
sstats <- round(100*abs(Scatter$X - Scatter$Y)/Scatter$Y)
print(summary(sstats))
breaks <- seq(-5,max(sstats,na.rm=TRUE)+5,by=5)
h <- hist(sstats,breaks=breaks)
plot(h$mids,h$count,log="xy",type="l",lwd=3,
     main="Percentage error")


print(paste(round(100*sum(sstats > 10,na.rm=TRUE)/sum(is.finite(sstats))),
            "% of the locations deviate with more than 10%"))
print(paste(round(100*sum(sstats > 25,na.rm=TRUE)/sum(is.finite(sstats))),
            "% of the locations deviate with more than 25%"))
print(paste(round(100*sum(sstats > 50,na.rm=TRUE)/sum(is.finite(sstats))),
            "% of the locations deviate with more than 50%"))
print(paste(round(100*sum(sstats > 100,na.rm=TRUE)/sum(is.finite(sstats))),
            "% of the locations deviate with more than 100%"))

dev.new()
xydensity <- table(3*round(Scatter$X/3),3*round(Scatter$Y/3))
contour(as.numeric(rownames(xydensity)),as.numeric(colnames(xydensity)),
        xydensity,xlim=c(0,200),ylim=c(0,200),
        main="Density of scatter between reconstruction and actual q95",
        xlab="Original (mm/day)",ylab="Reconstructed (mm/day)")
grid()
lines(c(0,200),c(0,200),lty=2,col="grey")

# Plot the locations where the largest discrepancies are found:
idodgy <- sstats > 100
print(table(attr(mu.qp.world,"country")[idodgy]))

dev.new()
h0 <- hist(attr(mu.qp.world,'n.wet'))
h1 <- hist(attr(mu.qp.world,'n.wet')[idodgy])
plot(h0$mids,h0$density,log="xy",type="l",lwd=3,
     main="Number of wet days",xlab="Number of wet days")
grid()
lines(h1$mids,h1$density,type="l",lwd=3,col="red")

dev.new()
h0 <- hist(attr(mu.qp.world,'mu'))
h1 <- hist(attr(mu.qp.world,'mu')[idodgy])
plot(h0$mids,h0$density,log="xy",type="l",lwd=3,
     main="Wet-day mean",xlab="Wet-day mean (mm/day)")
grid()
lines(h1$mids,h1$density,type="l",lwd=3,col="red")



dev.new()
plot(attr(mu.qp.world,"longitude")[idodgy],
       attr(mu.qp.world,"latitude")[idodgy],
       pch=19,col="red",main="Questionable rain guage records")
addland()

attr(mu.qp.world,'dodgy') <- idodgy
mu.qp.dodgy <- qweed(mu.qp.world,crit="attr(mu.qp.world,'dodgy')")
Xn <- estquantiles(mu.qp.dodgy)

dev.new()
plot(Xn,mu.qp.dodgy,pch=19,col="grey75",main="Questionable data",
     xlab="-log(1-p)*mu [mm/day]",ylab="quantile(X,p) [mm/day]")
v.dodgy <- pca$v[idodgy,]
N.dodgy <- length(pca$v[idodgy,1])
for (i in 1:N.dodgy) lines(pca$x1[,1]*pca$d[1]*pca$v[i,1]+
                            pca$x1[,2]*pca$d[2]*v.dodgy[i,2],
                            pca$qp[,1]*pca$d[1]*v.dodgy[i,1]+
                            pca$qp[,2]*pca$d[2]*v.dodgy[i,2],col="blue")
for (i in 1:N.dodgy) lines(pca$x1[,1]*pca$d[1]*v.dodgy[i,1],
                           pca$qp[,1]*pca$d[1]*v.dodgy[i,1],col="red")


lines(c(0,800),c(0,800),col="grey")
legend(600,200,c("Data","mode 1","mode 1+2"),
       col=c("grey70","red","blue"),pch=c(19,26,26),
       lty=c(0,1,1),bg="grey95")


# Regression analysis against distance from coast, mean precip, latitude:
print("Regression against geographical parameters")
print("x1 = mean precip; x2= altitude; x3=rnorm(N);")
print("x4-5 = latitude; x5-7=longitude;")
N <- length(pca$v[,1])

# If there is information about distance to coast, use it,
# otherwise use random noise
if (length(attr(MU.qp,"dist2coast.km"))>0) {
  print("x3 = distance to coast")
  d <- attr(MU.qp,"dist2coast.km")
} else {
  print("x3 = random noise")
  d<-rnorm(N)
}

exclude <- (1:dim(mu.qp.world)[2])[!is.finite(colMeans(mu.qp.world))]
print("Regression results for leading PC:")
calibr1 <- data.frame(y=pca$v[,1],x1=attr(mu.qp.world,"mean_precip")[-exclude],
                      x2=attr(mu.qp.world,"altitude")[-exclude],x3=d,
                      x4=sin(pi/180*attr(mu.qp.world,"latitude"))[-exclude],
                      x5=cos(pi/180*attr(mu.qp.world,"latitude"))[-exclude],
                      x6=sin(pi/180*attr(mu.qp.world,"longitude")[-exclude]),
                      x7=cos(pi/180*attr(mu.qp.world,"longitude"))[-exclude])
print(summary(lm(y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7,data=calibr1)))

print("Regression results for second PC:")
calibr2 <- data.frame(y=pca$v[,2],x1=attr(mu.qp.world,"mean_precip")[-exclude],
                      x2=attr(mu.qp.world,"altitude")[-exclude],x3=d,
                      x4=sin(pi/180*attr(mu.qp.world,"latitude"))[-exclude],
                      x5=cos(pi/180*attr(mu.qp.world,"latitude"))[-exclude],
                      x6=sin(pi/180*attr(mu.qp.world,"longitude"))[-exclude],
                      x7=cos(pi/180*attr(mu.qp.world,"longitude"))[-exclude])
print(summary(lm(y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7,data=calibr2)))

# Produce maps

statistics <-  c("q95","PC1","PC2","q95.pred","ratio",
                 "error","mu","mean_precip","n.wet","wetfreq",
                 "altitude","sigma","dist2coast")

for (i in 1:length(statistics)) {
  NorthAmerica(statistic=statistics[i]) -> output
  dev.copy2eps(file=paste("NorthAmerica_",statistics[i],".eps",sep=""))
}
 
