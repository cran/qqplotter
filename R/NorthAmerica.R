
Points <- function(x,y) {
  cols <- paste("grey",seq(0,100,length=101),sep="")
  cexs <- log(seq(5,1.25,length=101))
  dx <- ( max(x,na.rm=TRUE)-min(x,na.rm=TRUE) )/50000
  dy <- ( max(y,na.rm=TRUE)-min(y,na.rm=TRUE) )/50000
  for (i in 1:101) points(x+(i-1)*dx,y+(i-1)*dy,cex=cexs[i],col=cols[i]) 
}


Lines <- function(t,y,col="black",
                        tshadeoffset=0.03,yshadeoffset=0.03) {
  greys <- rgb( (seq(1,0.06,length=100)^0.1),
                (seq(1,0.06,length=100)^0.1),
                (seq(1,0.06,length=100)^0.1) )
  toffs <- tshadeoffset*(max(t,na.rm=TRUE)- min(t,na.rm=TRUE))*0.25
  yoffs <- yshadeoffset*(max(y,na.rm=TRUE)- min(y,na.rm=TRUE))
  for (i in 1:100) {
    ii <- ( (i - 100)/100 )
    lines(t-toffs+ii*toffs,y-yoffs+ii*yoffs,lwd=2,col=greys[i])
    lines(t-toffs-ii*toffs,y-yoffs-ii*yoffs,lwd=2,col=greys[i])
  }
  lines(t,y,lwd=5,col=col)
}

replace.char <- function (c, s, ny.c)  {
    while (instring(c, s)[1] > 0) {
        ii <- instring(c, s)
        if (ii > 1) 
            s <- paste(substr(s, 1, ii - 1), ny.c, substr(s, 
                ii + 1, nchar(s)), sep = "")
        else s <- paste(ny.c, substr(s, ii + 1, nchar(s)), sep = "")
    }
    s
}


makemap <- function(mu.qp=NULL,statistic="ratio",keep.longest=TRUE,
             lon.rng=c(-180,-50),lat.rng=c(25,70),stnr.given=TRUE,
            nx=200, ny=100, NC=64,zlim=NULL) {

# This code was provided by Doug Nychka, and it has been modified to
# work for a slightly larger region.

  require( LatticeKrig)
  require( PrecipStat)

# REB: use the greater data set, but weed out locations outside North America
# and locations with number of wet days fewer than 1000.

  if (is.null(mu.qp)) {
    data(mu.qp.world.jjas,envir = environment())
    mu.qp <- mu.qp.world.jjas; rm(mu.qp.world.jjas)
  }

  mu.qp<-  qweed(mu.qp,crit=paste("attr(mu.qp,'longitude')<",lon.rng[2]))
  mu.qp <- qweed(mu.qp,crit=paste("attr(mu.qp,'longitude')>",lon.rng[1]))
  mu.qp <- qweed(mu.qp,crit=paste("attr(mu.qp,'latitude')>",lat.rng[1]))
  mu.qp <- qweed(mu.qp,crit=paste("attr(mu.qp,'latitude')<",lat.rng[2]))
  mu.qp <- qweed(mu.qp,crit="is.finite(colSums(mu.qp))")

# Reve duplicated record, but keet the one with the most data:
  if (keep.longest) mu.qp <- cleanduplicates(mu.qp,silent=TRUE,plot=TRUE)

#mapofstations(mu.qp); dev.new(); stop()

  hold<- attributes(mu.qp)

# locations
  x<- cbind(hold$longitude, hold$latitude)
  #print("HERE"); print(dim(x))
  elev<- hold$altitude

#quilt.plot( x, elev, nx=256, ny=128)

# checking quantile data
  Q<- mu.qp
  mu<- hold$mu
  Pr<- hold$probabilities
  II<-Pr==.95
  Q95<- Q[II,]

  if (is.element(statistic,c("PC1","PC2","q95.pred","error"))) {
   pca <- qPCA(mu.qp)
   if ( (statistic=="q95.pred") | (statistic=="error") )
     q95.pred <- qPCA2quantile(pca$v[,1:2],pca,p=0.95)
  }

### fitting using Lattice Krig
# the data:

#print(summary(hold)); print(summary(hold$rmse));stop()

# Pick the desired statistic
  if (is.character(statistic)) y<- switch(statistic,
           "ratio"=Q95/( -log( 1-.95)*mu),
           "mu"=mu,
           "mean_precip"=hold$"mean_precip",
           "dist2coast"=hold$dist2coast.km,
           "n.wet"=hold$n.wet,
           "wetfreq"=round(100*hold$n.wet/
             (hold$n.wet+hold$n.dry)),
           "altitude"=hold$altitude,
           "rmse"=hold$rmse,
           "sigma"=hold$sigma,
           "q95"=Q95,
           "PC1"=pca$v[,1],
           "PC2"=pca$v[,2],
           "q95.pred"=q95.pred,
           "error"=Q95 - q95.pred) else
  if (is.numeric(statistic)) {
# If y is numeric/vector passed on - need to ensure
# that the locations held are the same as those yn mu.qp    
    Y <- statistic
    X <- cbind(attr(Y,'longitude'),attr(Y,'latitude'))
    statistic <- attr(statistic,'name')
    print(paste("y is a numeric passed on ... name=",statistic))
    print(c(length(Y),length(mu.qp[1,]))); print(dim(x))
    rep.info<- cat.matrix( X )
    good<- !duplicated(rep.info)
    y <- Y[good]; x2 <- X[good,]
    print(paste(length(y),"data points after check for duplicates and ",
          sum(is.finite(y)),"valid data;",
          " The dimensions of x2 are",dim(x2)[1],"&",dim(x2)[2]))

    if (stnr.given) {
    # Match station numbers by station number
      print("Find matching station numbers")
      matching.stnr <- is.element(attr(Y,'station_number')[good],
                                  attr(mu.qp,'station_number'))
      y <- y[matching.stnr]
      attr(y,'station_number') <- attr(Y,'station_number')[good][matching.stnr]
      # Check for duplicates in station numbers in y:
      stnr.check <- table(attr(y,'station_number'))
      stnr.dupl <- as.numeric(stnr.check)>1
      if (sum(stnr.dupl>0)) {
        print("Found dubplicatted station numbers in y")
        Stnr.Dup <- as.numeric(rownames(stnr.check)[stnr.dupl])
        print(Stnr.Dup)
        N.dup <- length(Stnr.Dup)
        for (II in 1:N.dup) {
          idup <- is.element(attr(y,'station_number'),Stnr.Dup[II])
          iII <- (1:length(y))[idup]
          print(paste(Stnr.Dup[II],sum(idup),iII[1],iII[2]))
          y[iII[2]] <- NA
        }
        good.stnrs <- attr(y,'station_number')[is.finite(y)]
        y <- y[is.finite(y)]
        attr(y,'station_number') <- good.stnrs
      }
      
      no.match.stnr <- !is.element(attr(mu.qp,'station_number'),
                                   attr(y,'station_number'))
      
      if (sum(!no.match.stnr)!= length(y)) {
        print("Check duplicate station numbers in mu.qp")
        stnr.check <- table(attr(mu.qp,'station_number'))
        stnr.dupl <- as.numeric(stnr.check)>1
        Stnr.Dup <- as.numeric(rownames(stnr.check)[stnr.dupl])
        print(Stnr.Dup)
        N.dup <- length(Stnr.Dup)
        for (II in 1:N.dup) {
          idup <- is.element(attr(mu.qp,'station_number'),Stnr.Dup[II])
          iII <- (1:length(y))[idup]
          print(paste(Stnr.Dup[II],sum(idup),iII[1],iII[2]))
          attr(mu.qp,'station_number')[iII[2]] <- NA
          x <- x[-iII[2],]
        }
        mu.qp <- qweed(mu.qp,crit="is.finite(attr(mu.qp,'station_number'))")
      }
      
      if (sum(no.match.stnr)>0) {
        # If some locations in mu.qp are not included in this region
        # remove these and their coordinates in x.
        print("mu.qp contains surplus stations")
        attr(mu.qp,'station_number')[no.match.stnr] <- NA
        mu.qp <- qweed(mu.qp,crit="is.finite(attr(mu.qp,'station_number'))")
        x <- x[attr(mu.qp,'weeded'),]
      }
      
    } else {
    # Match station numbers by longitudexlatitude
      print("Select data according to longitude/latitude information")
      #print(dim(mu.qp)); print(length(Y)); print(length(attr(Y,'longitude'))) 
      #print(length(good))
      y.lon.lat <- paste(round(attr(Y,'longitude')[good],1),
                         round(attr(Y,'latitude')[good],1),sep='x')
      mu.qp.lon.lat <- paste(round(attr(mu.qp,'longitude'),1),
                             round(attr(mu.qp,'latitude'),1),sep='x')
      #print(mu.qp.lon.lat[!is.na(mu.qp.lon.lat)])
      #print(paste("N good lons=",sum(is.finite(attr(mu.qp,'longitude')))))
      points(attr(mu.qp,'longitude'),attr(mu.qp,'latitude'),
             pch="x",col="darkgreen")
      #print(summary(round(attr(mu.qp,'longitude'),1)))
      #print(summary(round(attr(mu.qp,'latitude'),1)))
      matching.stnr <- is.element(y.lon.lat,mu.qp.lon.lat)
      print(paste("Found ",sum(matching.stnr),"lon/lat matches"))
      y <- y[matching.stnr]
    }

    attr(y,'name') <- statistic
    attr(y,'station_number') <- attr(Y,'station_number')[good][matching.stnr]
    attr(y,'longitude') <- attr(Y,'longitude')[good][matching.stnr]
    attr(y,'latitude') <- attr(Y,'latitude')[good][matching.stnr]

    print(paste("Number of stations extracted:",length(y),"=",
                length(mu.qp[1,]),"...?"))    
  }

# quick plot
  print("quilt plot"); print(dim(x)); print(length(y))
  print(summary(x)); print(summary(y))
  quilt.plot( x, y, nx=256, ny=128, main=statistic)
  dev2bitmap(file=paste("NAmerica_",statistic,"_quilt-plot.png",sep=""))

# find repeated stations
  if (!keep.longest) {
    rep.info<- cat.matrix( x)
    good<- !duplicated(rep.info)
    xM<- x[good,]
    yM<- y[good]
  } else {
    xM <- x; yM <- y
  }
# remove missing values
  good2<-  !is.na(yM)
  xM<- xM[good2,]
  yM<- yM[good2]

#################### example fit 
#data( elevation)
  print("Using 5-min topography data")
  data(etopo5,envir = environment()); elevation <- etopo5; rm(etopo5)
# logitude range 180W - 180E
  elevation$x[elevation$x>180] <-
  elevation$x[elevation$x>180]-360
  NAmerica.x <- (elevation$x > lon.rng[1]) &
                (elevation$x < lon.rng[2])
  NAmerica.y <- (elevation$y > lat.rng[1]) &
                (elevation$y < lat.rng[2])
  elevation$z <- elevation$z[NAmerica.x,NAmerica.y]
  elevation$x <- elevation$x[NAmerica.x]
  elevation$y <- elevation$y[NAmerica.y]
  elevation$z[elevation$z<0] <- 0
#image(elevation$x,elevation$y,elevation$z)

  print("create grid")
# create a 400X200 grid  
  gl0<- fields.x.to.grid( xM, nx=nx, ny=ny)
  names(gl0)<- c("x", "y")
# interpolate PRISM elevations to a grid
  Z0<- interp.surface.grid( elevation,gl0)
# elevations at station locations
  ZM<- interp.surface( elevation,xM)

# the range for lambda might need some adjustment for each problem
  lambda<- c(seq(-10,-5.5,,4), seq( 5.5,20,,4), seq(-5,5,,10) )
  lambda<- exp(sort(lambda))
  phi<-  seq( 4.5, 20,,15)
  a.wght<-  (4 + 1/phi^2)
  print(summary(xM)); print(summary(yM)); print(summary(c(ZM)))
  ZM[!is.finite(ZM)] <- 0

  out4<-MLE.LKrig( xM, yM,Z=ZM, NC=NC,alpha=1, a.wght=a.wght,
                  lambda=lambda, verbose=TRUE)
# thin plate spline like fit with a=4 
# out4a<-MLE.LKrig( xM, yM,Z=ZM,normalize=FALSE, NC=64,alpha=1, a.wght=4, lambda=lambda, verbose=TRUE)

# Quick look at likelihood values multiple points are due to collapsing plot on the a (range) parameter
#

  plot( out4[,c(5,9)], ylab="log likelihood", xlab="effective number of parameters", cex=.5)
  dev2bitmap(file=paste("NAmerica_",statistic,"_log-likelihood.png",sep=""))

# find maxmimum values of likelihood for this grid search
  I<- which.max(out4[,"lnProfileLike"])
  a4<-out4[I,2]
  lam4<- out4[I,3]
# fit with these parameters
#  fit4<- LKrig( xM, yM,Z=ZM, NC=64,  a.wght=a4, lambda=lam0)
  fit4<- LKrig( xM, yM,Z=ZM, NC=NC,  a.wght=a4, lambda=lam4)
# evaulating surface w/o elevation on the grid
  out.p<- predict.surface( fit4, grid.list=gl0, drop.Z=TRUE)
# now add in elevation part as a linear adjustment
  out.p$z<- out.p$z + Z0$z*c(fit4$d.coef[4] )
# set ocean grid points to NA
   mask<- in.land.grid( out.p)
  out.p$z[!mask] <- NA

#  image.plot( out.p, axes=FALSE, zlim=c(.5,2))

# "n.wet"

  print(paste("Select appropriate zlim for",statistic))
  if (is.null(zlim)) zlim <- switch(statistic,
                 "q95"=c(0,60),"q95.pred"=c(0,60),"error"=c(-3,3),
                 "q95.hat"=c(0,60),"q95.ext"=c(0,60),"dq95"=c(-10,10),
                 "q95.sce"=c(0,60),"q95dw"=c(0,40),"q95dw.sce"=c(0,40),
                 "ratio"=c(0.75,1.5),"sigma"=c(0,5000),"rmse"=c(0,500),
                 "mu"=c(0,15),"altitude"=c(0,3000),"wetfreq"=c(0,100),
                 "mean_precip"=c(0,8),"dist2coast"=c(0,1500),
                 "PC1"=c(-0.02,0.02),"PC2"=c(-0.03,0.03),"dq95dw"=c(-5,5),
                 "residual_PC1"=c(-0.02,0.02),"residual_PC1"=c(-0.03,0.03),
                 "n.wet"=c(0,40000),n.dry=c(0,40000),"dq95/q95"=c(-50,50),
                 "q95.sce/q95"=c(0.75,1.5),"dq95dw/q95dw"=c(-50,50),
                 "q95dw.sce/q95dw"=c(50,150),otherwise=NULL)
  col=two.colors(n=256, start="darkred", end="darkblue", middle="white")
  image.plot( out.p, axes=FALSE, col=col,main=statistic, zlim=zlim)
  world( add=TRUE, lwd=3)
#  contour(out.p,levels=c(1.1,1.2,1.3), col="grey", labex=0, add=TRUE)
  contour(out.p, col="grey", labex=0, add=TRUE)
  invisible(out.p)
  region <- paste("_",lon.rng[1],"E-",lon.rng[2],"E.",
                  lat.rng[1],"N-",lat.rng[2],"N",sep="")

  statistic <- replace.char("/",statistic,'.')
  dev2bitmap(file=paste("NAmerica_",statistic,region,".png",sep=""))
  dev2bitmap(file=paste("NAmerica_",statistic,region,".pdf",sep=""),
             type="pdfwrite")
  dev.copy2eps(file=paste("NAmerica_",statistic,region,".eps",sep=""))
}





modelqPCs <- function(mu.qp=NULL,keep.longest=TRUE,
             lon.rng=c(-180,-50),lat.rng=c(25,70),
             covariates=c("mu","altitude","dist2coast","wet.freq"),
             Eval.Region=c("NAmerica","WestEurope","EastEurope","Brazil",
                           "S.E.Asia")) {
  require( LatticeKrig)
  require( PrecipStat)

# REB: use the greater data set, but weed out locations outside North America
# and locations with number of wet days fewer than 1000.

  eval.region <- lower.case(Eval.Region)
  
  if (is.null(mu.qp)) {
    data(mu.qp.world.jjas,envir = environment())
    mu.qp <- mu.qp.world.jjas
  } else mu.qp.world.jjas <- mu.qp

# Complete set of validation purposes:
  pca.all <- qPCA(mu.qp,plot=FALSE); exclude <- attr(pca.all,'exclude')
  idep <- (attr(pca.all,'longitude') > min(lon.rng)) &
          (attr(pca.all,'longitude') < max(lon.rng)) &
          (attr(pca.all,'latitude') > min(lat.rng)) &
          (attr(pca.all,'latitude') < max(lat.rng)) 

#  The dependent data for treiining the global analysis

  calibr.dep <- data.frame(
                 pc1=pca.all$v[idep,1],
                 pc2=pca.all$v[idep,2],
                 mP=attr(mu.qp,"mean_precip")[-exclude][idep],
                 mu=attr(mu.qp,"mu")[-exclude][idep],
                 n.wet=attr(mu.qp,"n.wet")[-exclude][idep],
                 n.dry=attr(mu.qp,"n.dry")[-exclude][idep],
                 rmse=attr(mu.qp,"rmse")[-exclude][idep],
                 sigma=attr(mu.qp,"sigma")[-exclude][idep],
                 country=as.factor(attr(mu.qp,"country"))[-exclude][idep],
                 wet.freq=(round(100*attr(mu.qp,"n.wet")/
              (attr(mu.qp,"n.wet")+attr(mu.qp,"n.dry")),1))[-exclude][idep],
                 altitude=attr(mu.qp,"altitude")[-exclude][idep],
                 dist2coast=attr(mu.qp,"dist2coast.km")[-exclude][idep],
                 sin.lat=sin(pi/180*attr(mu.qp,"latitude"))[-exclude][idep],
                 cos.lat=cos(pi/180*attr(mu.qp,"latitude"))[-exclude][idep],
                 sin.lon=sin(pi/180*attr(mu.qp,"longitude"))[-exclude][idep],
                 cos.lon=cos(pi/180*attr(mu.qp,"longitude"))[-exclude][idep]  )

  # Independent data for evaluation
  calibr.ind <- data.frame(
                 pc1=pca.all$v[!idep,1],
                 pc2=pca.all$v[!idep,2],
                 mP=attr(mu.qp,"mean_precip")[-exclude][!idep],
                 mu=attr(mu.qp,"mu")[-exclude][!idep],
                 n.wet=attr(mu.qp,"n.wet")[-exclude][!idep],
                 n.dry=attr(mu.qp,"n.dry")[-exclude][!idep],
                 rmse=attr(mu.qp,"rmse")[-exclude][!idep],
                 sigma=attr(mu.qp,"sigma")[-exclude][!idep],
                 country=as.factor(attr(mu.qp,"country"))[-exclude][!idep],
                 wet.freq=(round(100*attr(mu.qp,"n.wet")/
              (attr(mu.qp,"n.wet")+attr(mu.qp,"n.dry")),1))[-exclude][!idep],
                 altitude=attr(mu.qp,"altitude")[-exclude][!idep],
                 dist2coast=attr(mu.qp,"dist2coast.km")[-exclude][!idep],
                 sin.lat=sin(pi/180*attr(mu.qp,"latitude"))[-exclude][!idep],
                 cos.lat=cos(pi/180*attr(mu.qp,"latitude"))[-exclude][!idep],
                 sin.lon=sin(pi/180*attr(mu.qp,"longitude"))[-exclude][!idep],
                 cos.lon=cos(pi/180*attr(mu.qp,"longitude"))[-exclude][!idep]  )

  
  mu.qp<-  qweed(mu.qp,crit=paste("attr(mu.qp,'longitude')<",max(lon.rng)))
  mu.qp <- qweed(mu.qp,crit=paste("attr(mu.qp,'longitude')>",min(lon.rng)))
  mu.qp <- qweed(mu.qp,crit=paste("attr(mu.qp,'latitude')>",min(lat.rng)))
  mu.qp <- qweed(mu.qp,crit=paste("attr(mu.qp,'latitude')<",max(lat.rng)))

# Reve duplicated record, but keet the one with the most data:
  if (keep.longest) mu.qp <- cleanduplicates(mu.qp,silent=TRUE,plot=TRUE)
#mapofstations(mu.qp); dev.new(); stop()

# Store the coordinates in 'x'
  x<- cbind(attr(mu.qp,'longitude'), attr(mu.qp,'latitude'))
  
  pca.cal <- qPCA(mu.qp,plot=FALSE)
  ii <- is.element(attr(mu.qp.world.jjas,'probabilities'),0.95)
  excl <- attr(pca.cal,'exclude')
  if (length(excl)==0) q95 <- c(mu.qp[ii,]) else
                       q95 <- c(mu.qp[ii,-excl])

  #HHH
  
  # Compare the North American and the World-wide EOFs:
#  dev.new()
#  plot(pca.cal$u[,1]*pca.cal$d[1],pca.all$u[,1]*pca.all$d[1],
#       xlab="PCA.cal",ylab="PCA.all",main="Leading EOF")
#  lines(c(0,800),c(0,800),col="grey")
#  dev2bitmap("pca-EOF1-comparison.png")

#  dev.new()
#  h1 <- hist(pca.cal$v[,1])
#  h2 <- hist(pca.all$v[,1])
#  plot(h1$mids,h1$density,type="l",lwd=2,
#       xlab="PC#1",ylab="Freq.",main="PC#1 - cal v.s. all",
#       sub="black: cal; red: all",xlim=range(c(h1$mids,h2$mids)))
#  lines(h2$mids,h2$density,lwd=2,col="red")
#  grid()
#  dev2bitmap("pca-PC1-comparison.png")

  N <- length(pca.cal$v[,1])

  print("quilt plot")
  dev.new()
  quilt.plot( x, pca.cal$v[,1], nx=256, ny=128, main="original PC#1")
  dev2bitmap("pca-PC1-original.png")
### fitting using Lattice Krig


# North American data for best-fit: dependent data
  calibr <- data.frame(
                      pc1=pca.cal$v[,1],
                      pc2=pca.cal$v[,2],
                      mP=attr(mu.qp,"mean_precip"),
                      mu=attr(mu.qp,"mu"),
                      n.wet=attr(mu.qp,"n.wet"),
                      n.dry=attr(mu.qp,"n.dry"),
                      rmse=attr(mu.qp,"rmse"),
                      sigma=attr(mu.qp,"sigma"),
                      country=as.factor(attr(mu.qp,"country")),
                      wet.freq=round(100*attr(mu.qp,"n.wet")/
                           (attr(mu.qp,"n.wet")+attr(mu.qp,"n.dry")),1),
                      altitude=attr(mu.qp,"altitude"),
                      dist2coast=attr(mu.qp,"dist2coast.km"),
                      sin.lat=sin(pi/180*attr(mu.qp,"latitude")),
                      cos.lat=cos(pi/180*attr(mu.qp,"latitude")),
                      sin.lon=sin(pi/180*attr(mu.qp,"longitude")),
                      cos.lon=cos(pi/180*attr(mu.qp,"longitude"))  )
  
  print("Regression results for two leading PCs:")
  modelexpr1 <- "pc1model <- lm(pc1 ~ "
  modelexpr2 <- "pc2model <- lm(pc2 ~ "
  for (i in 1:length(covariates))
    {
    if (i < length(covariates)) {
      modelexpr1 <- paste(modelexpr1,covariates[i],"+")
      modelexpr2 <- paste(modelexpr2,covariates[i],"+")
    } else {
      modelexpr1 <- paste(modelexpr1,covariates[i],",data=calibr)")
      modelexpr2 <- paste(modelexpr2,covariates[i],",data=calibr)")
    }
  }
  print(modelexpr1)
  print("Regression results for first PC:")
  eval(parse(text=modelexpr1))
  print(summary(pc1model))
  print(modelexpr2)
  print("Regression results for second PC:")
  eval(parse(text=modelexpr2))
  print(summary(pc2model))

# The model prediction based on dependent data: best-fit
  pc1.hat <- predict(pc1model,newdata=calibr)
  residual1 <- pca.cal$v[,1] - pc1.hat
  pc2.hat <- predict(pc2model,newdata=calibr)
  residual2 <- pca.cal$v[,2] - pc2.hat

  dev.new()  
  quilt.plot( x, pc1.hat, nx=256, ny=128, main="predicted PC#1")
  dev2bitmap("pca-PC1-predicted.png")

  dev.new()  
  quilt.plot( x, residual1, nx=256, ny=128, main="residual of PC#1")
  dev2bitmap("pca-PC1-residual.png")

  print("Predictions")  
  q95.hat <- qPCA2quantile(cbind(pc1.hat,pc2.hat),pca.cal,p=0.95)
  q95.rcs <- qPCA2quantile(pca.cal$v[,1:2],pca.cal,p=0.95)

  # Meta-data for passing vector as statistic to makemap()
  attr(q95.hat,'name') <- "q95.pred"
  attr(q95.hat,'longitude') <- attr(mu.qp,'longitude')
  attr(q95.hat,'latitude') <- attr(mu.qp,'latitude')
  attr(q95.hat,'station_number') <- attr(mu.qp,'station_number')

  attr(q95,'name') <- "q95"
  attr(q95,'longitude') <- attr(mu.qp,'longitude')
  attr(q95,'latitude') <- attr(mu.qp,'latitude')
  attr(q95,'station_number') <- attr(mu.qp,'station_number')

  attr(residual1,'name') <- "residual_PC1"
  attr(residual1,'longitude') <- attr(mu.qp,'longitude')
  attr(residual1,'latitude') <- attr(mu.qp,'latitude')
  attr(residual1,'station_number') <- attr(mu.qp,'station_number')
  
  dev.new()  

  plot(q95.hat,q95.rcs,pch=19,
     col="red",main="q95 from predicted PCs #1-2",
     xlab="From predicted PCs 1-2",ylab="From original PCs #1-2")
  lines(c(0,200),c(0,200))
  grid()
  Stats <- cor.test(c(q95.hat),c(q95.rcs))
  print(Stats)
  text(5,130,paste("r=",round(Stats$estimate,2),"(",
                   round(Stats$conf.int[1],2),"  -  ",
                   round(Stats$conf.int[2],2),")"),pos=4)
  text(5,125,paste("p-value=",round(100*Stats$p.value,2),"%"),pos=4)
  dev2bitmap("q95-predicted-from-PCs1-2.png")

  if (is.element("namerica",eval.region)) {
  # North America
    dev.new()
    makemap(mu.qp=mu.qp,statistic=q95.hat,keep.longest=keep.longest,
                 lon.rng=lon.rng,lat.rng=lat.rng)
    dev.new()
    makemap(mu.qp=mu.qp,statistic=q95,keep.longest=keep.longest,
                 lon.rng=lon.rng,lat.rng=lat.rng)
    dev.new()
    makemap(mu.qp=mu.qp,statistic=residual1,keep.longest=keep.longest,
                 lon.rng=lon.rng,lat.rng=lat.rng)
  }
  # Plot the locations of dependent data from PC1:
  dev.new()
  quilt.plot( cbind(attr(pca.all,'longitude')[idep],
                  attr(pca.all,'latitude')[idep]),
           pca.all$v[idep,1], nx=256, ny=128, main="PC#1: Dependent data")
  addland()
  dev2bitmap("pca-PC1-world-dependent.png")
  
  dev.new()
  quilt.plot( cbind(attr(pca.all,'longitude')[!idep],
                  attr(pca.all,'latitude')[!idep]),
           pca.all$v[!idep,1], nx=256, ny=128, main="PC#1: Independent data")
  addland()
  dev2bitmap("pca-PC1-world-independent.png")

  print("Global analysis:")
# Make a model based on the PCA for the whole world, but only use the
# subset of dependent variables: Capital 1st letter for global stuff
  Modelexpr1 <- paste("P",substr(modelexpr1,2,nchar(modelexpr1)-1),
                      ".dep)",sep="")
  eval(parse(text=Modelexpr1))
  Pc1.hat <- predict(Pc1model,newdata=calibr.ind)
  Modelexpr2 <- paste("P",substr(modelexpr2,2,nchar(modelexpr2)-1),
                      ".dep)",sep="")
  eval(parse(text=Modelexpr2))
  Pc2.hat <- predict(Pc2model,newdata=calibr.ind)

# checking quantile data
  Pr<- attr(mu.qp.world.jjas,'probabilities'); II<-Pr==.95
  Q95.all<- c(mu.qp.world.jjas[II,-exclude])
  Pr<- attr(mu.qp,'probabilities'); ii<-Pr==.95
  excl <- attr(mu.qp,'exclude')
  if (length(excl)==0) q95 <- c(mu.qp[ii,]) else
                       q95 <- c(mu.qp[ii,-excl])
  
  q95.ext <- qPCA2quantile(cbind(Pc1.hat,Pc2.hat),pca.all,p=0.95)
#  q95.ind <- qPCA2quantile(pca.all$v[!idep,1:2],pca.all,p=0.95)
  q95.ind <- Q95.all[!idep]

  print(c(length(q95.hat),length(q95)))
  
  #Evaluate the prediction of out-of-sample q95:
  dev.new()
  plot(q95.ext,Q95.all[!idep],
       ylab="Actual q95",xlab="predicted q95",
       main="Evaluation of q95-prediction",
       sub=paste(sub="modelqPCs(); N(ind)=",length(q95.ext),
             " N(dep)=",length(q95)),
       xlim=c(0,200),ylim=c(0,200),pch=19,col="blue")
  points(q95.hat,q95,pch=19,col="red")
  lines(c(0,600),c(0,600))
  stats <- cor.test(c(q95.ext),c(Q95.all[!idep]))
  print(stats)
  text(10,190,paste("r=",round(stats$estimate,2),"(",
                   round(stats$conf.int[1],2),"  -  ",
                   round(stats$conf.int[2],2),")"),pos=4)
  text(10,175,paste("p-value=",round(100*stats$p.value,2),"%"),pos=4)
  legend(165,20,c("Independent","Dependent"),pch=19,col=c("blue","red"),
                  cex=0.7,bg="grey95")
  dev2bitmap("modelqPCs-evaluation.png")
  
  # Meta-data for passing vector as statistic to makemap()
  # Insert [-exclude] to acount for NA-data removed in qPCA.
  attr(q95.ind,'name') <- "q95"
  attr(q95.ind,'longitude') <-
    attr(mu.qp.world.jjas,'longitude')[-exclude][!idep]
  attr(q95.ind,'latitude') <-
    attr(mu.qp.world.jjas,'latitude')[-exclude][!idep]
  attr(q95.ind,'station_number') <-
      attr(mu.qp.world.jjas,'station_number')[-exclude][!idep]

  attr(q95.ext,'name') <- "q95.pred"
  attr(q95.ext,'longitude') <-
    attr(mu.qp.world.jjas,'longitude')[-exclude][!idep]
  attr(q95.ext,'latitude') <-
    attr(mu.qp.world.jjas,'latitude')[-exclude][!idep]
  attr(q95.ext,'station_number') <-
      attr(mu.qp.world.jjas,'station_number')[-exclude][!idep]



  if (is.element("westeurope",eval.region)) {  
  # Evaluation: over W.Europe
    dev.new()
    print("Evaluate over western Europe:")
    makemap(mu.qp=mu.qp.world.jjas,statistic=q95.ind,
                 keep.longest=keep.longest,
                 lon.rng=c(-20,30),lat.rng=c(30,65))
    makemap(mu.qp=mu.qp.world.jjas,statistic=q95.ext,
                 keep.longest=keep.longest,
                 lon.rng=c(-20,30),lat.rng=c(30,65))
  }

  if (is.element("easteurope",eval.region)) {  
  
    # Evaluation: over E.Europe
    print("Evaluate over Eastern Europe:")
    makemap(mu.qp=mu.qp.world.jjas,statistic=q95.ind,
                 keep.longest=keep.longest,
                 lon.rng=c(20,70),lat.rng=c(30,90))
    makemap(mu.qp=mu.qp.world.jjas,statistic=q95.ext,
                 keep.longest=keep.longest,
                 lon.rng=c(20,70),lat.rng=c(30,90))
  }

  if (is.element("brazil",eval.region)) {  
    # Evaluation: over Brazil
    print("Evaluate over Brazil:")
    makemap(mu.qp=mu.qp.world.jjas,statistic=q95.ind,
                 keep.longest=keep.longest,zlim=c(0,80),
                 lon.rng=c(-100,-30),lat.rng=c(-60,10))
    makemap(mu.qp=mu.qp.world.jjas,statistic=q95.ext,
                 keep.longest=keep.longest,zlim=c(0,80),
                 lon.rng=c(-100,-30),lat.rng=c(-60,10))
}

  if (is.element("s.e.asia",eval.region)) {    
    # Evaluation: over South-East Asia
    print("Evaluate over South-East Asia:")
    makemap(mu.qp=mu.qp.world.jjas,statistic=q95.ind,
                 keep.longest=keep.longest,zlim=c(0,120),
                 lon.rng=c(70,150),lat.rng=c(0,50))
    makemap(mu.qp=mu.qp.world.jjas,statistic=q95.ext,
                 keep.longest=keep.longest,zlim=c(0,120),
                 lon.rng=c(70,150),lat.rng=c(0,50))
  }


  if (is.element("australia",eval.region)) {    
      # Evaluation: over Australia
    print("Evaluate over Australia:")
    #print(dim(mu.qp)); print(length(attr(mu.qp,'longitude')))
    makemap(mu.qp=mu.qp.world.jjas,statistic=q95.ind,
                 keep.longest=keep.longest,stnr.given=FALSE,
                 lon.rng=c(100,180),lat.rng=c(-50,-10))
    makemap(mu.qp=mu.qp.world.jjas,statistic=q95.ext,
                 keep.longest=keep.longest,stnr.given=FALSE,
                 lon.rng=c(100,180),lat.rng=c(-50,-10))
  }

  if (is.element("japan",eval.region)) {    
      # Evaluation: over Australia
    print("Evaluate over Japan:")
    makemap(mu.qp=mu.qp.world.jjas,statistic=q95.ind,
                 keep.longest=keep.longest,stnr.given=FALSE,
                 lon.rng=c(120,150),lat.rng=c(30,50))
    makemap(mu.qp=mu.qp.world.jjas,statistic=q95.ext,
                 keep.longest=keep.longest,stnr.given=FALSE,
                 lon.rng=c(120,150),lat.rng=c(30,50))
  }

  if (is.element("norway",eval.region)) {    
      # Evaluation: over Norway
    print("Evaluate over Norway:")
    makemap(mu.qp=mu.qp.world.jjas,statistic=q95.ind,
                 keep.longest=keep.longest,stnr.given=TRUE,
                 lon.rng=c(0,40),lat.rng=c(55,70))
    makemap(mu.qp=mu.qp.world.jjas,statistic=q95.ext,
                 keep.longest=keep.longest,stnr.given=TRUE,
                 lon.rng=c(0,40),lat.rng=c(55,70))
  }

  if (is.element("baltic",eval.region)) {    
      # Evaluation: over the Baltic
    print("Evaluate over Baltic:")
    makemap(mu.qp=mu.qp.world.jjas,statistic=q95.ind,
                 keep.longest=keep.longest,stnr.given=TRUE,
                 lon.rng=c(-10,50),lat.rng=c(45,65))
    makemap(mu.qp=mu.qp.world.jjas,statistic=q95.ext,
                 keep.longest=keep.longest,stnr.given=TRUE,
                 lon.rng=c(-10,50),lat.rng=c(45,65))
  }

  if (is.element("safrica",eval.region)) {    
      # Evaluation: over the South Africa
    print("Evaluate over South Africa:")
    makemap(mu.qp=mu.qp.world.jjas,statistic=q95.ind,
                 keep.longest=keep.longest,
                 lon.rng=c(0,40),lat.rng=c(-45,-20))
    makemap(mu.qp=mu.qp.world.jjas,statistic=q95.ext,
                 keep.longest=keep.longest,
                 lon.rng=c(0,40),lat.rng=c(-45,-20))
  }
  
  
  print("Save the results:")
  result <- list(pc1model=pc1model,pc2model=pc2model,pca=pca.cal,
                 Pc1model,Pc2model,pca.all=pca.all,idep=idep,mu.qp=mu.qp)

  save(file="modelqPCs.rda",result)
  invisible(result)
}

qqscenario <- function(result=NULL,type=c("dq95dw/q95dw","q95dw.sce")) {

  # This function is a crude and dirty way to demonstrate the
  # principles of using qqplotter in downscaling.
  require( PrecipStat)

  print("A simple & crude demonstration of how change in mu")
  print("and f.wet can alter the 95-th quantile")
  
  if (is.null(result)) load("modelqPCs.rda")
  mu.qp <- result$mu.qp
  
  data(f.wet.1990.1999,envir = environment())
  data(mu.1990.1999,envir = environment())
  data(f.wet.2080.2099,envir = environment())
  data(mu.2080.2099,envir = environment())
  df.wet <- round(100*f.wet.2080.2099,1) -
            round(100*f.wet.1990.1999,1)
  dmu <- mu.2080.2099 - mu.1990.1999
  nx <- length(attr(f.wet.1990.1999,'longitude'))
  ny <- length(attr(f.wet.1990.1999,'latitude'))
  dev.new()

  # Colour scheme for the graphics:
  redblue <- rgb( c(seq(1,0.1,length=30),rep(0,20)),
               rep(0,50),
               c(rep(0,20),seq(0.1,1,length=30)) )

  image(attr(f.wet.1990.1999,'longitude'),
        attr(f.wet.1990.1999,'latitude'),df.wet,
        main="d f_wet",col=redblue,zlim=c(-30,30)); addland()
  contour(attr(f.wet.1990.1999,'longitude'),
        attr(f.wet.1990.1999,'latitude'),df.wet,
        main="d f_wet",add=TRUE,col="white",lwd=1)
  dev2bitmap("qq-dfwet.png",res=150) 
  
  dev.new()
  image(attr(mu.1990.1999,'longitude'),
        attr(mu.1990.1999,'latitude'),dmu,
        main="d mu",col=redblue,zlim=c(-7,7)); addland()
  contour(attr(mu.1990.1999,'longitude'),
        attr(mu.1990.1999,'latitude'),dmu,
        main="d mu",add=TRUE,col="white",lwd=1)
  dev2bitmap("qq-dmu.png",res=150) 

  # Scale the GCM results so that the wet-day mean and wet-day frequency
  # for the first time-slice corresponds with the corresponding values
  # from the stations in the same region:

#  mu.gcm <- mean(c(mu.1990.1999),na.rm=TRUE)
#  fwet.gcm <- mean(c(f.wet.1990.1999),na.rm=TRUE)
#  scale.mu <- round(mean(attr(mu.qp,'mu'),na.rm=TRUE)/mu.gcm,4)
#  wet.freq <- mean(attr(mu.qp,"n.wet")/
#             (attr(mu.qp,"n.wet")+attr(mu.qp,"n.dry")),na.rm=TRUE)
#  scale.fwet <- round(wet.freq/fwet.gcm,4)
  mu.gcm <- c(mu.1990.1999)
  fwet.gcm <- c(f.wet.1990.1999)
  scale.mu <- round(attr(mu.qp,'mu')/mu.gcm,4)
  wet.freq <- attr(mu.qp,"n.wet")/
             (attr(mu.qp,"n.wet")+attr(mu.qp,"n.dry"))
  scale.fwet <- round(wet.freq/fwet.gcm,4)

  print(paste("Mean mu from GDCN=",round(mean(attr(mu.qp,'mu'),na.rm=TRUE),1),
              "mean mu from GCM=",round(mean(mu.gcm),1),"mean f.wet from GDCN=",
              round(100*wet.freq),"% mean f.wet from GCM=",round(100*mean(fwet.gcm)),
              " -> mu scaling factor=",mean(scale.mu)," f.wet scaling factor=",
              mean(scale.fwet)))
  
  N <- length(result$pca$v[,1])

  if (file.exists("qqscenario.rda")) {
    load("qqscenario.rda")
    if (length(Dmu)==N) interpolate <- FALSE
  } else interpolate <- TRUE
 
  if (interpolate) {
    lonxy <- rep(attr(mu.1990.1999,'longitude'),ny)
    latxy <- sort(rep(attr(mu.1990.1999,'latitude'),nx))
    if (file.exists("qqscenario.temp.rda")) {
      print("Starting from previous unfinished results")
      load("qqscenario.temp.rda")
      print(paste("i1=",i1,"N=",N))
    } else {
      Dmu <- rep(NA,N); Df.wet <- Dmu
      i1 <- 1
    }

    print(range(lonxy))
    print(range(latxy))
    
    for (i in i1:N) {
      Dmu[i] <- interp(lonxy,latxy,dmu,
                       attr(result$pca,'longitude')[i],
                       attr(result$pca,'latitude')[i])$z
      Df.wet[i] <- interp(lonxy,latxy,df.wet,
                       attr(result$pca,'longitude')[i],
                       attr(result$pca,'latitude')[i])$z
     if (mod(i,10)==0) cat(".")
     if (mod(i,100)==0) {
       cat(Dmu[i]); cat(":"); cat(attr(result$pca,'longitude')[i])
       cat("x"); cat(attr(result$pca,'latitude')[i])
     }
     if (mod(i,1000)==0) {
         cat("o")
         i1 <- i + 1
         save(file="qqscenario.temp.rda",i1,Dmu,Df.wet)
       }
    }
    file.remove("qqscenario.temp.rda")
    save(file="qqscenario.rda",Dmu,Df.wet)
  }
  Dmu <- Dmu * scale.mu; Df.wet <- Df.wet*scale.fwet

  print(c( length(attr(mu.qp,'mu')),length(Dmu),
           length(attr(result$pca,'exclude')) ))
  
  print("Quilt plot of mu & f.wet:")  
  quilt.plot( cbind(attr(result$pca,'longitude'),
                   attr(result$pca,'latitude')), Dmu, nx=256, ny=128,
             col=two.colors(start="darkred", end="darkblue", middle="white"))
  addland()
  dev2bitmap("Dmu-quilt.png")
  quilt.plot( cbind(attr(result$pca,'longitude'),
                   attr(result$pca,'latitude')), Df.wet, nx=256, ny=128,
             col=two.colors(start="darkred", end="darkblue", middle="white"))
  addland()
  dev2bitmap("Df.wet-quilt.png")
  
  predictors0 <- data.frame(mu=attr(mu.qp,'mu'),
            wet.freq=round(100*attr(mu.qp,"n.wet")/
         (attr(mu.qp,"n.wet")+attr(mu.qp,"n.dry")),1),
            altitude=attr(mu.qp,'altitude'),
            dist2coast=attr(mu.qp,'dist2coast.km'))
  # Predict PCs 1 & 2 from variations in mu and wet.freq:
  pc1.0 <- predict(result$pc1model,newdata=predictors0)
  pc2.0 <- predict(result$pc2model,newdata=predictors0)

  # Use predicted PCs to reconstruct q95, given the other PCA products
  q95.0 <- qPCA2quantile(cbind(pc1.0,pc2.0),result$pca,p=0.95)

  predictors1 <- data.frame(mu=attr(mu.qp,'mu') + Dmu,
            wet.freq=round(100*attr(mu.qp,"n.wet")/
         (attr(mu.qp,"n.wet")+attr(mu.qp,"n.dry")),1) + Df.wet,
            altitude=attr(mu.qp,'altitude'),
            dist2coast=attr(mu.qp,'dist2coast.km'))
  # Predict PCs 1 & 2 from variations in mu and wet.freq:
  pc1 <- predict(result$pc1model,newdata=predictors1)
  pc2 <- predict(result$pc2model,newdata=predictors1)

  # Use predicted PCs to reconstruct q95, given the other PCA products
  q95.sce <- qPCA2quantile(cbind(pc1,pc2),result$pca,p=0.95)
  dq95 <- q95.sce - q95.0

  dev.new()
  attr(dq95,'name') <- "dq95"
  attr(dq95,'longitude') <- attr(result$mu.qp,'longitude')
  attr(dq95,'latitude') <- attr(result$mu.qp,'latitude')
  attr(dq95,'station_number') <-
      attr(result$mu.qp,'station_number')
  #print("make maps: dq95")
  if (sum(is.element(type,"dq95"))>0)
    makemap(mu.qp=mu.qp,statistic=dq95)

  attr(q95.sce,'name') <- "q95.sce"
  attr(q95.sce,'longitude') <- attr(result$mu.qp,'longitude')
  attr(q95.sce,'latitude') <- attr(result$mu.qp,'latitude')
  attr(q95.sce,'station_number') <-
      attr(result$mu.qp,'station_number')
  #print("make maps: q95.sce")  
  if (sum(is.element(type,"q95.sce"))>0)
    makemap(mu.qp=mu.qp,statistic=q95.sce)

  dq95.q95 <- 100*(q95.sce - q95.0)/q95.0
  attr(dq95.q95,'name') <- "dq95/q95"
  attr(dq95.q95,'longitude') <- attr(result$mu.qp,'longitude')
  attr(dq95.q95,'latitude') <- attr(result$mu.qp,'latitude')
  attr(dq95.q95,'station_number') <-
     attr(result$mu.qp,'station_number')
  if (sum(is.element(type,"dq95/q95"))>0)
    makemap(mu.qp=mu.qp,statistic=dq95.q95)
  
  print("make maps for wet + dry days")
  f.wet <- attr(mu.qp,"n.wet")/
         (attr(mu.qp,"n.wet")+attr(mu.qp,"n.dry"))
  q95dw <- q95.0; q95dw[] <- NA; q95dw.sce <- q95dw
  
  # The 95th percentiles considering all days: wet + dry
  # p is Pr(A|B) where Pr(A) = Pr(X>x) and P(B) = f.wet.
  # Pr(A) = Pr(A|B)Pr(B)/Pr(B|A) (Baye's theorem). Pr(B|A)=1
  # Denote Pr(A|B) with '1-p' and Pr(A) with '1-P'.
  # 1-P = (1-p) f.wet
  print("Estimating q95 for dry+wet days")
  #print(attr(mu.qp,'probabilities'))
  for (i in 1:N) {
      p.all <- 1 - (1-0.95)/f.wet[i]
      #print(round(c(i,f.wet[i],p.all,range(mu.qp[,i]),pc1[i],pc2[i]),2))
      if ( (is.finite(p.all)) & (sum(is.finite(mu.qp[,i]))>2) &
           (is.finite(pc1[i])) & (is.finite(pc2[i])) ) {
        q95dw[i] <- approx(attr(mu.qp,'probabilities'),c(mu.qp[,i]),p.all)$y
        q95dw.sce[i] <- qPCA2quantile(c(pc1[i],pc2[i]),
                                    result$pca,p=p.all,silent=TRUE)
        if (mod(i,100)==0) cat(".")
        if (mod(i,1000)==0) cat(i)
      }
  }
  attr(q95dw,'name') <- "q95dw"
  attr(q95dw,'longitude') <- attr(result$mu.qp,'longitude')
  attr(q95dw,'latitude') <- attr(result$mu.qp,'latitude')
  attr(q95dw,'station_number') <-
      attr(result$mu.qp,'station_number')
  if (sum(is.element(type,"q95dw"))>0)
    makemap(mu.qp=mu.qp,statistic=q95dw)

  attr(q95dw.sce,'name') <- "q95dw.sce"
  attr(q95dw.sce,'longitude') <- attr(result$mu.qp,'longitude')
  attr(q95dw.sce,'latitude') <- attr(result$mu.qp,'latitude')
  attr(q95dw.sce,'station_number') <-
      attr(result$mu.qp,'station_number')
  if (sum(is.element(type,"q95dw.sce"))>0)
    makemap(mu.qp=mu.qp,statistic=q95dw.sce)

  dq95dw <- q95dw.sce - q95dw 
  attr(dq95dw,'name') <- "dq95dw"
  attr(dq95dw,'longitude') <- attr(result$mu.qp,'longitude')
  attr(dq95dw,'latitude') <- attr(result$mu.qp,'latitude')
  attr(dq95dw,'station_number') <-
  attr(result$mu.qp,'station_number')
  if (sum(is.element(type,"dq95dw"))>0)
    makemap(mu.qp=mu.qp,statistic=dq95dw)

  dq95dw.dq95 <- 100*(q95dw.sce - q95dw)/q95dw 
  attr(dq95dw.dq95,'name') <- "dq95dw/q95dw"
  attr(dq95dw.dq95,'longitude') <- attr(result$mu.qp,'longitude')
  attr(dq95dw.dq95,'latitude') <- attr(result$mu.qp,'latitude')
  attr(dq95dw.dq95,'station_number') <-
     attr(result$mu.qp,'station_number')
  if (sum(is.element(type,"dq95dw/q95dw"))>0)
    makemap(mu.qp=mu.qp,statistic=dq95dw.dq95)

  q95dw.dq95 <- 100*q95dw.sce/q95dw 
  attr(q95dw.dq95,'name') <- "q95dw.sce/q95dw"
  attr(q95dw.dq95,'longitude') <- attr(result$mu.qp,'longitude')
  attr(q95dw.dq95,'latitude') <- attr(result$mu.qp,'latitude')
  attr(q95dw.dq95,'station_number') <-
     attr(result$mu.qp,'station_number')
  if (sum(is.element(type,"q95dw.sce/q95dw"))>0)
    makemap(mu.qp=mu.qp,statistic=q95dw.dq95)  

}



qqatribution <- function(mu.qp=NULL,result=NULL,
                         path="~/GDCN/",x.0=1,months=NULL) {
  require( PrecipStat)

  if (is.null(mu.qp)) data(mu.qp,envir = environment())
  if (is.null(result)) load("modelqPCs.rda")
  # This function explores whether the model developed in modelqPCs()
  # is able to predict changes in time for a random location. 

  # Find the stations with longest data record:
  ndays <- attr(mu.qp,'n.wet') + attr(mu.qp,'n.dry')
  longestrecord <- is.element(ndays,max(ndays,na.rm=TRUE))
  stnr <- attr(mu.qp,'station_number')[longestrecord]
  print(paste("Station with the longest record: ndays=",max(ndays,na.rm=TRUE),
              "station number=",stnr))
  d2c <- attr(mu.qp,"dist2coast.km")[longestrecord]
  z <- attr(mu.qp,"altitude")[longestrecord]
  
  # Read the time series for the longest station record:
  list <- list.files(path,pattern=".dly",full.names = TRUE)
  data(gdcn.inv,envir = environment())
  filename <- list[grep(stnr,list)]
  X <- readGDCN(filename)
#  TX <- readGDCN(filename,type="TMAX")
#  TN <- readGDCN(filename,type="TMIN")
#  TM <- 0.5*(TX+TN)
  meanP <- round(mean(X,na.rm=TRUE),2)
  n.dry <- sum(X < x.0,na.rm=TRUE)
  n.wet <- sum(X >= x.0,na.rm=TRUE)
  sigma <- round(var(c(X[X >= x.0]),na.rm=TRUE),2)
  imatch <- is.element(gdcn.inv$stnr,attr(X,"Station_number"))
  if (sum(imatch)>0) {
    lat <- gdcn.inv$lat[imatch]
    lon <- gdcn.inv$lon[imatch]
    alt <- gdcn.inv$alt[imatch]
    stnr <- attr(X,"Station_number")
    country <- attr(X,"Country_code")
  }
  if (!is.null(months)) X <- X[is.element(attr(X,"month"),months),]

  yrx5 <- 5*trunc(attr(X,'year')/5)
  pentads <- as.numeric( rownames(table(yrx5)) ); N <- length(pentads)
  mu <- rep(NA,N); f.wet <- mu; q95 <- mu
  for ( i in 1:N) {
    ipentad <- is.element(yrx5,pentads[i])
    x <- X[ipentad,];
#    tas <- mean(TM[ipentad,],na.rm=TRUE)
    fwet <- sum(x >= x.0,na.rm=TRUE)/sum(is.finite(x))
    x <- x[x>=x.0]
    mu[i] <- mean(x,na.rm=TRUE)
    q95[i] <- quantile(x,0.95,na.rm=TRUE)
  }

#  tas <- tas - mean(tas,na.rm=TRUE)
  predictors <- data.frame(mu=mu,wet.freq=fwet,
                           altitude=rep(z,N),
                           dist2coast=rep(d2c,N))
  # Predict PCs 1 & 2 from variations in mu and wet.freq:
  pc1 <- predict(result$pc1model,newdata=predictors)
  pc2 <- predict(result$pc2model,newdata=predictors)

  # Use predicted PCs to reconstruct q95, given the other PCA products
  q95.hat <- qPCA2quantile(cbind(pc1,pc2),result$pca,p=0.95)

  # Plot the comparison
  plot(pentads,q95,pch=19,type="n",
       ylab="wet-day q95 (mm/day)",xlab="year",
       main=paste("Station #",stnr,"prediction of q95(t)"),
       sub="prediction of PCs 1&2; reconstruction of q95")
  Lines(pentads,q95.hat)
  Points(pentads,q95)
  
  stats <- cor.test(c(q95),c(q95.hat))
  print(stats)
  text(min(pentads),max(q95),paste("r=",round(stats$estimate,2),"(",
                   round(stats$conf.int[1],2),"  -  ",
                   round(stats$conf.int[2],2),")"),pos=4)
  text(min(pentads),max(q95)*0.975,
       paste("p-value=",round(100*stats$p.value,2),"%"),pos=4)
  legend(max(pentads)-20,1.05*min(q95,na.rm=TRUE),
         c("predicted","observed"),pch=c(26,19),
         col=c("black","grey"),lty=c(1,0),lwd=3,cex=0.7,bg="grey95")
  dev2bitmap("qqatribution-evaluation.png",res=150)
}




qpAnywhere <- function(lon,lat,mu.qp=NULL,p=0.95) {
  require( PrecipStat)

  if (is.null(mu.qp)) {
    data(mu.qp.world,envir = environment())
    mu.qp <- mu.qp.world; rm(mu.qp.world)
  }
  dist <- round(distAB(lon,lat,attr(mu.qp,'longitude'),
                       attr(mu.qp,'latitude'))/1000)
  ii <- is.element(dist,min(dist,na.rm=TRUE))
  f.wet <- ( attr(mu.qp,"n.wet")/
         ( attr(mu.qp,"n.wet")+attr(mu.qp,"n.dry")) )[ii]
  print(paste("Station #",attr(mu.qp,'station_number')[ii]," at ",
        attr(mu.qp,'longitude')[ii],"E/",attr(mu.qp,'latitude')[ii],
        "N is ",dist[ii]," km away. p=",p," f.wet=",round(100*f.wet,1),
              "%",sep=""))
              
  ip <- is.element(attr(mu.qp,'probabilities'),p)
  if (sum(ip)>0) qp <- mu.qp[ip,ii] else {
    qp <- approx(attr(mu.qp,'probabilities'),c(mu.qp[,ii]),p)$y
  }

# p is Pr(A|B) where Pr(A) = Pr(X>x) and P(B) = f.wet.
# Pr(A) = Pr(A|B)Pr(B)/Pr(B|A) (Baye's theorem). Pr(B|A)=1
# Denote Pr(A|B) with '1-p' and Pr(A) with '1-P'.
# 1-P = (1-p) f.wet
  p.all <- 1 - (1-p)/f.wet
  qp.all <- approx(attr(mu.qp,'probabilities'),c(mu.qp[,ii]),p.all)$y
  print(paste("The ",p*100,"-percentile for wet-days is ",qp,"mm/day",
              "and for all days, q_p is ",round(qp.all,1),
              "mm/day; the wet-day mean is ",attr(mu.qp,'mu')[ii],
              " and the mean for all days is ",
              attr(mu.qp,'mean_precip')[ii],sep=""))
  result <- list(lon=attr(mu.qp,'longitude')[ii],lat=attr(mu.qp,'latitude')[ii],
                 stnr=attr(mu.qp,'station_number')[ii],dist=dist[ii],
                 qp=qp,qp.all=qp.all,f.wet=f.wet,mu=attr(mu.qp,'mu')[ii],
                 mean.precip=attr(mu.qp,'mean_precip')[ii],
                 elevation=attr(mu.qp,'altitude')[ii])
  invisible(result)
}




exploreqmodel <- function(mu.qp=NULL,result=NULL) {
  require( PrecipStat)

  if (is.null(mu.qp)) data(mu.qp,envir = environment())
  if (is.null(result)) load("modelqPCs.rda")

  par(mfcol=c(2,2))
  mu <- mean(attr(mu.qp,'mu'),na.rm=TRUE)
  fwet <- mean( attr(mu.qp,"n.wet")/
         ( attr(mu.qp,"n.wet")+attr(mu.qp,"n.dry")),na.rm=TRUE )
  z <- 1000; d2c <- 1000

  N <- 100
  
  predictors <- data.frame(mu=mu*seq(0.5,2,length=N),
                           wet.freq=rep(fwet,N),
                           altitude=rep(z,N),
                           dist2coast=rep(d2c,N))
  # Predict PCs 1 & 2 from variations in mu and wet.freq:
  pc1 <- predict(result$pc1model,newdata=predictors)
  pc2 <- predict(result$pc2model,newdata=predictors)

  # Use predicted PCs to reconstruct q95, given the other PCA products
  q95.hat <- qPCA2quantile(cbind(pc1,pc2),result$pca,p=0.95,silent=TRUE)
  plot(mu*seq(0.5,2,length=N),q95.hat,xlab="mu",ylab="q95",
       main="q95=f(mu)")
  lines(mu*seq(0.5,2,length=N),-log(0.05)*mu*seq(0.5,2,length=N),
        col="red",lty=2)

  
  predictors <- data.frame(mu=rep(mu,N),
                           wet.freq=fwet*seq(0.5,2,length=N),
                           altitude=rep(z,N),
                           dist2coast=rep(d2c,N))
  # Predict PCs 1 & 2 from variations in mu and wet.freq:
  pc1 <- predict(result$pc1model,newdata=predictors)
  pc2 <- predict(result$pc2model,newdata=predictors)

  # Use predicted PCs to reconstruct q95, given the other PCA products
  q95.hat <- qPCA2quantile(cbind(pc1,pc2),result$pca,p=0.95,silent=TRUE)
  plot(fwet*seq(0.5,2,length=N),q95.hat,xlab="f.wet",ylab="q95",
       main="q95=f(f.wet)")


  
  predictors <- data.frame(mu=rep(mu,N),
                           wet.freq=rep(fwet,N),
                           altitude=z*seq(0.5,2,length=N),
                           dist2coast=rep(d2c,N))
  # Predict PCs 1 & 2 from variations in mu and wet.freq:
  pc1 <- predict(result$pc1model,newdata=predictors)
  pc2 <- predict(result$pc2model,newdata=predictors)

  # Use predicted PCs to reconstruct q95, given the other PCA products
  q95.hat <- qPCA2quantile(cbind(pc1,pc2),result$pca,p=0.95,silent=TRUE)
  plot(z*seq(0.5,2,length=N),q95.hat,xlab="z",ylab="q95",
       main="q95=f(z)")

  predictors <- data.frame(mu=rep(mu,N),
                           wet.freq=rep(fwet,N),
                           altitude=rep(z,N),
                           dist2coast=d2c*seq(0.5,2,length=N))
  # Predict PCs 1 & 2 from variations in mu and wet.freq:
  pc1 <- predict(result$pc1model,newdata=predictors)
  pc2 <- predict(result$pc2model,newdata=predictors)

  # Use predicted PCs to reconstruct q95, given the other PCA products
  q95.hat <- qPCA2quantile(cbind(pc1,pc2),result$pca,p=0.95,silent=TRUE)
  plot(d2c*seq(0.5,2,length=N),q95.hat,xlab="d2c",ylab="q95",
       main="q95=f(d2c)")

  dev2bitmap("exploreqmodel.png",res=150)      
}



rqstats <- function(mu.qp=NULL,result=NULL,
                    path="~/GDCN/",x.0=1,months=NULL) {
  require( PrecipStat)

  if (is.null(mu.qp)) data(mu.qp,envir = environment())
  if (is.null(result)) load("modelqPCs.rda")
  # This function explores whether the model developed in modelqPCs()
  # is able to predict changes in time for all stations

  if (!file.exists("rqstats.rda")) {
    N <- length(mu.qp[1,])
    ndays <- attr(mu.qp,'n.wet') + attr(mu.qp,'n.dry')
    d2c <- attr(mu.qp,"dist2coast.km")
    z <- attr(mu.qp,"altitude")
    list <- list.files(path,pattern=".dly",full.names = TRUE)
    data(gdcn.inv,envir = environment())
    rq <- rep(NA,N); rn <- rq

    if (file.exists("rqstats.temp.rda")) {
      load("rqstats.temp.rda")
      print("continuing on previously unfinished results")
      print(paste("i1=",i1,"N=",N))
    } else i1 <- 1
          
    for (i in i1:N) {
      stnr <- attr(mu.qp,'station_number')[i]
  
  # Read the time series 
      filename <- list[grep(stnr,list)]
      X <- readGDCN(filename)
      meanP <- round(mean(X,na.rm=TRUE),2)
      n.dry <- sum(X < x.0,na.rm=TRUE)
      n.wet <- sum(X >= x.0,na.rm=TRUE)
      sigma <- round(var(c(X[X >= x.0]),na.rm=TRUE),2)
      imatch <- is.element(gdcn.inv$stnr,attr(X,"Station_number"))
      if (sum(imatch)>0) {
        lat <- gdcn.inv$lat[imatch]
        lon <- gdcn.inv$lon[imatch]
        alt <- gdcn.inv$alt[imatch]
        stnr <- attr(X,"Station_number")
        country <- attr(X,"Country_code")
      }
      if (!is.null(months)) X <- X[is.element(attr(X,"month"),months),]

      yrx5 <- 5*trunc(attr(X,'year')/5)
     pentads <- as.numeric( rownames(table(yrx5)) ); n <- length(pentads)
     if (n > 10) {
        mu <- rep(NA,n); f.wet <- mu; q95 <- mu
        for ( ii in 1:n) {
          ipentad <- is.element(yrx5,pentads[ii])
          x <- X[ipentad,];
          fwet <- sum(x >= x.0,na.rm=TRUE)/sum(is.finite(x))
          x <- x[x>=x.0]
          mu[ii] <- mean(x,na.rm=TRUE)
          q95[ii] <- quantile(x,0.95,na.rm=TRUE)
        }

        predictors <- data.frame(mu=mu,wet.freq=fwet,
                                 altitude=rep(z[i],n),
                                 dist2coast=rep(d2c[i],n))
  # Predict PCs 1 & 2 from variations in mu and wet.freq:
        pc1 <- predict(result$pc1model,newdata=predictors)
        pc2 <- predict(result$pc2model,newdata=predictors)

  # Use predicted PCs to reconstruct q95, given the other PCA products
        q95.hat <- qPCA2quantile(cbind(pc1,pc2),result$pca,p=0.95,silent=TRUE)

        rq[i] <- cor(c(q95),c(q95.hat))
        rn[i] <- cor(rnorm(n),rnorm(n))
        if (mod(i,100)==0) {
          plot(rq,ylim=c(-1,1),pch=19); points(rn,col="red")
          i1 <- i+1
          save(file="rqstats.temp.rda",i1,rq,rn)
          cat("o")
        }
        if (mod(i,10)==0) cat(".")
      }
    }
    file.remove("rqstats.temp.rda")
    save(file="rqstats.rda",rq,rn)
  } else load("rqstats.rda")
  dev2bitmap("qqtimeseries-correlations.png",res=150)
  
  h0 <- hist(rn); h1 <- hist(rq)
  
  # Plot the comparison
  plot(range(c(h0$mids,h1$mids)),range(c(h0$counts,h1$counts)),type="n",
       ylab="Counts",xlab="Correlation",
       main="Skill of prediction of q95(t)",
       sub=paste(N,"stations"))
  Lines(h1$mids,h1$counts)
  lines(h0$mids,h0$counts,lty=2,col="red")
  
  dev2bitmap("qqtimeseries-evaluation.png",res=150)
}

