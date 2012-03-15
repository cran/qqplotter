# R.E. Benestad (rasmus.benestad@physics.org)
# R-script for reading CDCN daily precipitation data
# NCAR, Mesa Lab, Boulder March 31, 2011

# Look for similarity between exponential distribution and 24-hr precip world wide
# Look for systematic relationship between mu, meanP, and meanT, hence a dependence of shape of p.d.f. on
# local mean T(2m) and precip.


# Reconstruct the quantiles from the formulae
estquantiles <- function(mu.qp,
                         prs=c(seq(0.50,0.95,by=0.05),0.96,0.97,0.98,0.99)) {

# The formula depends on the type of input: Matrix or a single number/vector
  if ( (class(mu.qp)!="qqplotter") & ((class(mu.qp)!="matrix")) ) {
    Xn <- -log(1-prs)*mu.qp
  } else {
    #print("estquantiles - detected a matrix")
    mu <- attr(mu.qp,"mu")
    if (length(attr(mu.qp,'probabilities'))>0)
      prs <- attr(mu.qp,"probabilities")
    Xn <- -log(1-prs)%o%mu   # matrix outer product
    #print(dim(mu.qp)); print(c(length(mu),length(prs))); print(dim(Xn))
  }
  attr(Xn,"probabilities") <- prs
  Xn
}






fitsimple <- function(X,x.0=1,
                      prs=c(seq(0.50,0.95,by=0.05),0.96,0.97,0.98,0.99)) {
  # This function does the actual calculation/estimation of quantiles
  if (is.matrix(X)) {
    print("fitsimple: detected a matrix")
    print(dim(X))
    X[X<x.0] <- NA
    # Assuming matrices of [space,time] dimensions
    mu <- apply(X,1,mean,na.rm=TRUE)
    x2 <- apply(X,1,quantile,prs,na.rm=TRUE)
    #print(c(length(mu),length(prs),length(x2)))
    x1 <- estquantiles(mu,prs)
  } else {
    wdp <- X[X > x.0]
    mu <- mean(wdp,na.rm=TRUE)
    x1 <- estquantiles(mu,prs)
    x2 <- quantile(wdp,prs,na.rm=TRUE)
  }
  #print(paste("fitsimple sizes - x1:",length(x1)," x2:",length(x2)))
  attr(x1,"description") <- "-log(1-prs)*mu"
  attr(x2,"description") <- "quantile(X,p)"
  results <- data.frame(x=x1,y=x2)
  attr(results,"mu") <- mu
  attr(results,"x.0") <- x.0
  invisible(results)
}





expMC <- function(mu,p,N,N.runs=1000) {
  # Carries out Monte-Carlo simulations that describes exponentially
  # distributed data.
  q.sim <- rep(NA,N.runs)
  for (i in 1:N.runs) q.sim[i] <- quantile(rexp(N,1/mu),p)
  attr(q.sim,"description" <- "Monte-Carlo simulation")
  attr(q.sim,"probability") <- p
  attr(q.sim,"record_length") <- N
  attr(q.sim,"mu") <- mu
  invisible(q.sim)
}






qplotMC <- function(mu=10,N=90,
                    prs=c(seq(0.50,0.95,by=0.05),0.96,0.97,0.98,0.99),
                    N.runs=1000,plot=TRUE) {
  # Makes qqplots of the Monte-Carlo simulations and uses these to
  # estimate 5-95% confidence interval.
  np <- length(prs)
  q.sim <- matrix(rep(NA,N.runs*np),np,N.runs)
  for (i in 1:np) {
    q.sim[i,] <- expMC(mu,prs[i],N,N.runs=N.runs)
  }
  Xn <- estquantiles(mu,prs=prs) 
  lower <- apply(q.sim,1,quantile,0.05)
  upper <- apply(q.sim,1,quantile,0.95)
  if (plot) {
    plot(range(q.sim),range(q.sim),type="n",
         main="Monte-Carlo simulation")
    for (i in 1:N.runs) lines(Xn,q.sim[,i],col="pink")
    lines(c(0,1000),c(0,1000),col="grey")
    lines(Xn,lower,lty=2)
    lines(Xn,upper,lty=2)
  }
  print(dim(q.sim))
  CI <- list(Xn=Xn,lower=lower,upper=upper)
  invisible(CI)
}





qPCA <- function(mu.qp,plot=TRUE) {
  # Carries out PCA (EOF) analysis of the cloud of points.
  if (class(mu.qp)!="qqplotter") stop("qPCA: need 'qqplotter' object")
  Xn <- estquantiles(mu.qp)
  X <- rbind(Xn,mu.qp)
  d <- dim(X); n1 <- dim(mu.qp)[1]
  exclude <- (1:d[2])[!is.finite(colMeans(X))]
  if (length(exclude)>0) X <- X[,-exclude]
  pca <- svd(X)
  if (plot) {
    plot(100*pca$d^2/sum(pca$d^2),type="b",lty=2,
         main="PCA variance",xlab="mode #",ylab="variance (%)",
         sub="[qqplotter - qPCA() - pca$var]")
    grid()
  #dev.copy2eps(file="qqplotter5a.eps")
    dev.new()
  }
  # The vectors consist of percentiles according the exponential
  # distribution and the actual percentiles.
  x1<- pca$u[1:n1,]
  x2 <-pca$u[(n1+1):d[1],]
  attr(x2,"description") <- "-log(1-prs)*mu"
  attr(x1,"description") <- "quantile(X,p)"
  pca$x1 <- x1
  pca$qp <- x2
  pca$var <- 100*pca$d^2/sum(pca$d^2)
  if (plot) {
    plot(abs(x1[,1]),abs(x2[,1]),
         main="Leading PCA mode",
         xlab="-log(1-p)*mu [mm/day]",ylab="quantile(X,p) [mm/day]",
         sub="[qqplotter - qPCA()]")
    lines(c(0,1),c(0,1),col="grey",lty=2)
    grid()
  }
  #dev.copy2eps(file="qqplotter5b.eps")
  if (length(exclude)>1) incl <- -exclude else
                         incl <- 1:length(attr(mu.qp,"mu"))
  attr(pca,"description") <- attr(mu.qp,"description")
  attr(pca,"mean_precip") <- attr(mu.qp,"mean_precip")[incl]
  attr(pca,"mu") <- attr(mu.qp,"mu")[incl]
  attr(pca,"longitude") <- attr(mu.qp,"longitude")[incl]
  attr(pca,"latitude") <- attr(mu.qp,"latitude")[incl]
  attr(pca,"station_number") <- attr(mu.qp,"station_number")[incl]
  attr(pca,"altitude") <- attr(mu.qp,"altitude")[incl]
  attr(pca,"n.wet") <- attr(mu.qp,"n.wet")[incl]
  attr(pca,"n.dry") <- attr(mu.qp,"n.dry")[incl]
  attr(pca,"sigma") <- attr(mu.qp,"sigma")[incl]
  attr(pca,"rmse") <- attr(mu.qp,"rmse")[incl]
  attr(pca,"country") <- attr(mu.qp,"country")[incl]
  attr(pca,"history") <- paste(attr(mu.qp,"history")," & qPCA")
  attr(pca,"probabilities") <- attr(mu.qp,"probabilities")
  attr(pca,"formulae") <- attr(mu.qp,"formulae")
  attr(pca,"x.0") <- attr(mu.qp,"x.0")
  attr(pca,"source") <- attr(mu.qp,"source")
  attr(pca,"n1") <- n1
  attr(pca,"exclude") <- exclude
  if (length(attr(mu.qp,"rcm"))>0)
    attr(pca,"rcm") <- attr(mu.qp,"rcm")[incl]
  if (length(attr(mu.qp,"dist2coast.km"))>0)
    attr(pca,"dist2coast.km") <- attr(mu.qp,"dist2coast.km")[incl]
  class(pca) <- "qPCA"
  invisible(pca)
}


qPCA2quantile <- function(x,pca,p=0.95,silent=FALSE) {
  # Reconstruct the original data from PCA results.
  if (!silent) print("qPCA2quantile:")
  if (class(pca)!="qPCA") stop("The pca-argument must be a qPCA object")
  attr(pca,"n1") -> n1; d <- dim(pca$u)
  
  if (class(x)=="matrix") {
    D <- dim(x); n.eofs <- D[2];
    if (!silent) print(paste("Taking x to be a matrix with",n.eofs,"modes"))    
  # Reconstruct the data from the EOFs, replacing the original PCs with x
    X <- pca$u[,1:n.eofs] %*% diag(pca$d[1:n.eofs]) %*% t(x)
  # Extract the part of the data describing the 'actual' quantiles:
    if (!silent) {print(dim(X));print(paste("n1+1=",(n1+1),"d[1]=",d[1]))}
    x2 <-X[(n1+1):d[1],]
  } else if (class(x)=="numeric") {
    if (length(x) <=2) {
      if (!silent)
        print("Assuming the vector describes 2 leading modes from one location")
      n.eofs <- length(x)
      #print(dim(pca$u[,1:n.eofs] %*% diag(pca$d[1:n.eofs])))
      X <- ( pca$u[,1:n.eofs] %*% diag(pca$d[1:n.eofs]) ) %o% t(x)
   # Extract the part of the data describing the 'actual' quantiles:
      x2 <-X[(n1+1):d[1]]
    } else {
       if (!silent)
    print("Assuming the vector describes the leading mode from many location")
      X <- c(pca$u[,1] * pca$d[1]) %o% c(x)
   # Extract the part of the data describing the 'actual' quantiles:
       if (!silent) print(dim(X))
      x2 <-X[(n1+1):d[1],]
    }   
  }
  ip <- is.element(attr(pca,"probabilities"),p)
  #print("HERE");print(length(ip)); print(dim(x2))
  if (sum(ip)==1) {
    if (is.matrix(x2)) percentile <- x2[ip,] else
                       percentile <- x2[ip]
  } else {
    # Interpolate
    if (is.matrix(x2)) {
      print("Interpolate")
      percentile <- rep(NA,dim(x2)[2])
      for (i in 1:dim(x2)[2])
        percentile[i] <- approx(attr(pca,"probabilities"),x2,p)$y
    } else {
      percentile <- approx(attr(pca,"probabilities"),x2,p)$y
    }
  }
  attr(percentile,"probability") <- p
  attr(percentile,"method") <- "qPCA2quantile"
  invisible(percentile)
}


testqPCA2quantile <- function(mu.qp=NULL,subset=1:1000,n.eofs=2) {
  # Test the reconstruction of data from PCA.
  require(PrecipStat)
  print("testqPCA2quantile:")
  if (is.null(mu.qp)) data(mu.qp,envir = environment())
  pca <- qPCA(mu.qp)
  N <- length(pca$v[,1])
  x <- rep(NA,N); y <- x

  # If set to NULL, use the comlete PCA decomposition.
  if (is.null(subset)) subset <- 1:N
  if (is.null(n.eofs)) n.eofs <- length(pca$u[,1])
  
  dev.new()
  plot(c(0,100),c(0,100),type="l",lwd=2,col="grey",
       main="95th-quantile reconstruction",xlab="Original",
       ylab="Reconstructed from qPCA object",
       sub=paste("testqPCA2quantile: the first",length(subset),
         "locations"))
  legend(0,100,c("Single site","Multiple sites/EOFs",
                 "Multiple sites - one EOF"),pch=c(19,21,4),
         col=c("blue","red","grey"),bg="grey95")
  grid()
  
  
  print("Examine 95th percentiles for one site at the time")
  print("Only using the 2 leading EOFs")
  for (i in subset[1:min(1000,length(subset))]) {
    x[i] <- mu.qp[is.element(attr(mu.qp,"probabilities"),0.95),i]
    qPCA2quantile(pca$v[i,1:n.eofs],pca,p=0.95) -> y[i]
    points(x[i],y[i],pch=19,col="blue",cex=0.75)
  }

  print(paste("Examine 95th percentiles for ",length(subset),
              " sites simulataneously"))
  print(paste("Using the ",n.eofs," leading EOFs"))
  qPCA2quantile(pca$v[subset,1:n.eofs],pca,p=0.95) -> Y
  X <- mu.qp[is.element(attr(mu.qp,"probabilities"),0.95),subset]
  print(c(length(X),length(Y)))
  points(X,Y,pch=21,col="red",cex=0.75)

  print(paste("Examine 95th percentiles for ",length(subset),
              "sites simulataneously"))
  print("Only using the leading EOF")
  qPCA2quantile(c(pca$v[subset,1]),pca,p=0.95) -> Y.1
  X.1 <- mu.qp[is.element(attr(mu.qp,"probabilities"),0.95),subset]
  points(X.1,Y.1,pch=4,col="grey",cex=0.75)

  results <- list(x=x,y=y,X=X,Y=Y,X.1=X.1,Y.1=Y.1)
  invisible(results)
}




colorbar <-  function(ckey,fig=c(0.33,0.66,0.30,0.34)) {
  # Graphics - Produces a colour bar
  print("colorbar")
  fig.old <- c(0,1,0,1)
  mar0 <- par()$mar
  par(fig=fig,new=TRUE,mar=rep(0,4),xaxt="n",yaxt="n",cex.axis=0.5)
  breaks <- seq(0,1,length=length(ckey)+1)
  colbar <- rbind(c(breaks),c(breaks)); d <- dim(colbar)
  image(colbar,breaks=breaks,col=ckey)
  par(fig=fig.old,mar=mar0,new=TRUE)
}





mapofstations <- function(mu.qp,colourcoding="0.95",add=FALSE,silent=FALSE,
                 googleearth=FALSE,kmz.file="~/Desktop/qqplotter.kmz") {
  # Plots the location of the stations in the qqplotter object.
  if (class(mu.qp)!="qqplotter") stop("mapofstations: need 'qqplotter' object")
  require(clim.pact)
  if (substr(colourcoding,1,2)=="0.") {
    mu <- attr(mu.qp,"mu")
    crit <- estquantiles(mu.qp)
    colourcoding=paste(colourcoding,"quantile")
    extrainfo <- paste("calculated from -log(1-p) mu, ",
      "assuming exponential distribution.")
  } else {
    attr(mu.qp,colourcoding) -> crit
    if (colourcoding=="rmse") extrainfo <- paste("measure of difference",
          "to an exponential distribution") else
                              extrainfo <- "derived from observations"
  }
  crit <- crit - min(crit,na.rm=TRUE)
  crit[crit > quantile(crit,0.99,na.rm=TRUE)] <- quantile(crit,0.99,na.rm=TRUE)

  Pmax <- max(crit,na.rm=TRUE) - min(crit,na.rm=TRUE)
  N <- length(crit)
  cols <- rgb(rep(0.5,N),rep(0.5,N),rep(0.5,N))
  ok <- is.finite(crit)
  cols[ok] <- rgb((1-(crit[ok]/Pmax))^0.5,rep(0,sum(ok)),(crit[ok]/Pmax)^0.5)

  if (silent) colourcoding=""
  if (!add) {
    par(las=1)
    plot(attr(mu.qp,"longitude"),attr(mu.qp,"latitude"),
                type="n",main=colourcoding,xlab="",ylab="",ylim=c(-90,90),
                xlim=range(attr(mu.qp,"longitude"),na.rm=TRUE)+c(-45,0))
  }
  points(attr(mu.qp,"longitude"),attr(mu.qp,"latitude"),
         col=cols,pch=19,cex=0.5)
  addland()
  s <- seq(0,1,length=100)
  colkey1 <- rgb( (1-s)^0.5,rep(0,100),s^0.5)
  #print(Pmax); print(summary(crit))
  text(-205,-73,colourcoding,cex=0.75,srt=90,pos=4)
  text(-215,-52,round(max(crit,na.rm=TRUE),1),cex=0.75,col="grey30")
  text(-215,-82,round(min(crit,na.rm=TRUE),1),cex=0.75,col="grey30")
  colorbar(colkey1,c(0.15,0.18,0.22,0.30))

  if (googleearth)  {
    # Makes a GoogleEarth file, so the locations can be shown in
    # GoogleEarth.
    kmz.cont <- c('<?xml version="1.0" encoding="UTF-8"?>',
             '<kml xmlns="http://www.opengis.net/kml/2.2">',
             '<Document>')
    for (i in 1:N) {
      pid <- paste(attr(mu.qp,"location")[i])
      info <- paste("Station number=", attr(mu.qp,"station_number")[i],
      "; Statistics: mean precip (all days)=",
                    attr(mu.qp,"mean_precip")[i],
      "mm/day; mean wet-day precip=",round(attr(mu.qp,"mean_precip"),1)[i],
      "mm/day; number of wet days=",attr(mu.qp,"n.wet")[i],
      "; number  of dry days=",attr(mu.qp,"n.dry")[i],
      "; wet-day standard deviation=",round(attr(mu.qp,"sigma"),1)[i],
      "; RMSE of quantiles against exponential distribution",
       round(attr(mu.qp,"rmse")[i],2))
      b <- c('<Folder>',  # 11 elements as common header
           paste('<Placemark id="',pid,'">',sep=""),
           paste('<name>',attr(mu.qp,"location")[i],'</name>'),
           '<description>',
           '<![CDATA[',paste('qqplotter results:',
                             attr(mu.qp,"location")[i],colourcoding,
                             round(crit[i],1),"mm/day",extrainfo,
                             info),
           ']]>','</description>','<Point>',
           paste('<coordinates>',attr(mu.qp,"longitude")[i],',',
                 attr(mu.qp,"latitude")[i],',',
                 attr(mu.qp,"altitude")[i],'</coordinates>',sep=""),
           '</Point>','</Placemark>','</Folder>')
      kmz.cont <- c(kmz.cont,b)  
    }
    kmz.cont <- c(kmz.cont,'</Document>','</kml>')
    writeLines(kmz.cont, kmz.file)
  }
  invisible(cols)
}





qcat <- function(mu.qp1,mu.qp2) {
  # A tool for concatinating two qqplotter objects.
  if (class(mu.qp1)!="qqplotter" | class(mu.qp2)!="qqplotter")
    stop("qcat: need 'qqplotter' object")
  prs1 <- attr(mu.qp1,"probabilities")
  prs2 <- attr(mu.qp2,"probabilities")
  if (!identical(prs1,prs2))
    stop("qcat: objects need same probability levels")
  months1 <- attr(mu.qp1,"months")
  months2 <- attr(mu.qp2,"months")
  if (!identical(months1,months2))
    stop("qcat: objects need to include same months")
  
  mu.qp <- cbind(mu.qp1,mu.qp2)
  attr(mu.qp,"description") <- c(attr(mu.qp1,"description"),
                                 attr(mu.qp2,"description"))
  attr(mu.qp,"mean_precip") <- c(attr(mu.qp1,"mean_precip"),
                                 attr(mu.qp2,"mean_precip"))
  attr(mu.qp,"mu") <- c(attr(mu.qp1,"mu"),
                        attr(mu.qp2,"mu"))
  attr(mu.qp,"longitude") <- c(attr(mu.qp1,"longitude"),
                               attr(mu.qp2,"longitude"))
  attr(mu.qp,"latitude") <- c(attr(mu.qp1,"latitude"),
                              attr(mu.qp2,"latitude"))
  attr(mu.qp,"station_number") <- c(attr(mu.qp1,"station_number"),
                                    attr(mu.qp2,"station_number"))
  attr(mu.qp,"altitude") <- c(attr(mu.qp1,"altitude"),
                              attr(mu.qp2,"altitude"))
  attr(mu.qp,"n.wet") <- c(attr(mu.qp1,"n.wet"),
                           attr(mu.qp2,"n.wet"))
  attr(mu.qp,"n.dry") <- c(attr(mu.qp1,"n.dry"),
                           attr(mu.qp2,"n.dry"))
  attr(mu.qp,"sigma") <- c(attr(mu.qp1,"sigma"),
                           attr(mu.qp2,"sigma"))
  attr(mu.qp,"rmse") <- c(attr(mu.qp1,"rmse"),
                          attr(mu.qp2,"rmse"))
  attr(mu.qp,"country") <- c(attr(mu.qp1,"country"),
                             attr(mu.qp2,"country"))
  attr(mu.qp,"history") <- c(attr(mu.qp1,"history"),
                             attr(mu.qp2,"history"),"qcat()")
  if ( (length(attr(mu.qp1,"rcm"))>0) &
       (length(attr(mu.qp2,"rcm"))>0) )
    attr(mu.qp,"rcm") <- c(attr(mu.qp1,"rcm"),
                           attr(mu.qp2,"rcm"))
  if ( (length(attr(mu.qp1,"dist2coast.km"))>0) &
       (length(attr(mu.qp2,"dist2coast.km"))>0) )
    attr(mu.qp,"dist2coast.km") <- c(attr(mu.qp1,"dist2coast.km"),
                                     attr(mu.qp2,"dist2coast.km"))
  attr(mu.qp,"probabilities") <- prs1
  attr(mu.qp,"months") <- months1
  attr(mu.qp,"formulae") <- attr(mu.qp1,"formulae")
  class(mu.qp) <- "qqplotter"
  invisible(mu.qp)
}

qsubset <- function(mu.qp,crit="attr(mu.qp,'probabilities') <= 0.99") {
  print(paste("qsubset(select='",crit,"')"))
  eval(parse(text=paste("keep=",crit)))
  print(paste("Keeping",sum(keep,na.rm=TRUE),"percentiles"))
  mu.qp <- mu.qp.all[keep,]
  attr(mu.qp,"description") <- attr(mu.qp.all,"description")
  attr(mu.qp,"mean_precip") <- attr(mu.qp.all,"mean_precip")
  attr(mu.qp,"mu") <- attr(mu.qp.all,"mu")
  attr(mu.qp,"longitude") <- attr(mu.qp.all,"longitude")
  attr(mu.qp,"latitude") <- attr(mu.qp.all,"latitude")
  attr(mu.qp,"station_number") <- attr(mu.qp.all,"station_number")
  attr(mu.qp,"altitude") <- attr(mu.qp.all,"altitude")
  attr(mu.qp,"n.wet") <- attr(mu.qp.all,"n.wet")
  attr(mu.qp,"n.dry") <- attr(mu.qp.all,"n.dry")
  attr(mu.qp,"sigma") <- attr(mu.qp.all,"sigma")
  attr(mu.qp,"rmse") <- attr(mu.qp.all,"rmse")
  attr(mu.qp,"country") <- attr(mu.qp.all,"country")
  attr(mu.qp,"history") <- attr(mu.qp.all,"history")
  attr(mu.qp,"probabilities") <- attr(mu.qp.all,"probabilities")
  attr(mu.qp,"formulae") <- attr(mu.qp.all,"formulae")
  attr(mu.qp,"dist2coast.km") <- attr(mu.qp.all,"dist2coast.km")
  attr(mu.qp,"probabilities") <- attr(mu.qp.all,"probabilities")[keep]
  class(mu.qp) <- class(mu.qp.all)
  invisible(mu.qp)
}


qweed <- function(mu.qp,crit="attr(mu.qp,'n.wet')> 1000") {
  # A tool to weed out stations according to a given criteria.
  if (class(mu.qp)!="qqplotter") stop("qweed: need 'qqplotter' object")
  print(paste("qweed(weed='",crit,"')"))
  eval(parse(text=paste("iweed=",crit)))
  print(paste("Keeping",sum(iweed,na.rm=TRUE),"data points"))
  print(c(length(iweed),dim(mu.qp)))
  weeded <- (1:length(mu.qp[1,]))[iweed]
  mu.qpx <- mu.qp[,iweed]
  attr(mu.qpx,"description") <- attr(mu.qp,"description")
  attr(mu.qpx,"mean_precip") <- attr(mu.qp,"mean_precip")[iweed]
  attr(mu.qpx,"mu") <- attr(mu.qp,"mu")[iweed]
  attr(mu.qpx,"longitude") <- attr(mu.qp,"longitude")[iweed]
  attr(mu.qpx,"latitude") <- attr(mu.qp,"latitude")[iweed]
  attr(mu.qpx,"station_number") <- attr(mu.qp,"station_number")[iweed]
  attr(mu.qpx,"altitude") <- attr(mu.qp,"altitude")[iweed]
  attr(mu.qpx,"n.wet") <- attr(mu.qp,"n.wet")[iweed]
  attr(mu.qpx,"n.dry") <- attr(mu.qp,"n.dry")[iweed]
  attr(mu.qpx,"sigma") <- attr(mu.qp,"sigma")[iweed]
  attr(mu.qpx,"rmse") <- attr(mu.qp,"rmse")[iweed]
  attr(mu.qpx,"country") <- attr(mu.qp,"country")[iweed]
  attr(mu.qpx,"history") <- attr(mu.qp,"history")
  attr(mu.qpx,"probabilities") <- attr(mu.qp,"probabilities")
  attr(mu.qpx,"months") <- attr(mu.qp,"months")
  attr(mu.qpx,"formulae") <- attr(mu.qp,"formulae")
  attr(mu.qpx,"weeded") <- weeded
  if (length(attr(mu.qp,"rcm"))>0)
    attr(mu.qpx,"rcm") <- attr(mu.qp,"rcm")[iweed]
  if (length(attr(mu.qp,"dist2coast.km"))>0)
    attr(mu.qpx,"dist2coast.km") <- attr(mu.qp,"dist2coast.km")[iweed]
  class(mu.qpx) <- "qqplotter"
  invisible(mu.qpx)
}

cleanduplicates <- function(mu.qp,silent=TRUE,plot=FALSE) {
  mu.qp <- qweed(mu.qp,crit="is.finite(attr(mu.qp,'longitude'))")
  mu.qp <- qweed(mu.qp,crit="is.finite(attr(mu.qp,'latitude'))")
  lons <- round(attr(mu.qp,"longitude"),4)
  lats <- round(attr(mu.qp,"latitude"),4)
  if (plot) {plot(lons,lats,cex=0.5); addland()}
  N <- length(lons)
  lonlats <- paste(round(lons,4),"E/",round(lats,4),"N",sep="")
  lonlat.tab <- table(lonlats)
  dupl <- rownames(lonlat.tab)[lonlat.tab > 1 ]
  N.dupl <- length(dupl)
  if (N.dupl>0) {
    if (!silent) print(paste(N.dupl,"duplicates"))
    for (i in 1:N.dupl) {
      ii <- (1:N)[is.element(lonlats,dupl[i])]
      if ( (length(ii)>1) & (length(ii)<10) ){
        if (plot) points(lons[ii],lats[ii],pch=19,col="red",cex=0.5)
        if (!silent) print(c(i,sum(ii),NA,
                attr(mu.qp,"longitude")[ii],NA,attr(mu.qp,"latitude")[ii],
                NA,attr(mu.qp,"n.wet")[ii]))
        srt <- order(attr(mu.qp,"n.wet")[ii], decreasing = TRUE)
        attr(mu.qp,"n.wet")[ii][srt][2:sum(ii)] <- 0
        #print(attr(mu.qp,"n.wet")[ii])
      } else if (length(ii)>=10) {
        print(paste("i=",i,"Something suspitious happened, as there were ",
                    length(ii),"stations with the same coordinates:",
                    dupl[i]))
        print(paste("sum(lonlat.tab[lonlat.tab > 1 ])=",
                    sum(lonlat.tab[lonlat.tab > 1 ])))
        if (plot) {
          dev.new()
          plot(attr(mu.qp,"longitude"),attr(mu.qp,"latitude"),
             pch=19,col="grey",xlim=c(-180,180),ylim=c(-90,90)); addland()
          points(attr(mu.qp,"longitude")[ii],attr(mu.qp,"latitude")[ii],
             pch=19,col="red")
        }
        print(attr(mu.qp,"longitude")[ii])
        print(attr(mu.qp,"latitude")[ii])
        print(attr(mu.qp,"n.wet")[ii])
        print(lonlat.tab[lonlat.tab > 1 ])
        stop("Halt in cleanduplicates - unexpected results")      
      }
    }
    if (!silent) print(paste("Number of locations BEFORE weeding=",
                             length(mu.qp)," and ",
                             sum(attr(mu.qp,'n.wet')==0),"locations had",
           "'n.wet' set to zero"))
    mu.qp <- qweed(mu.qp,crit="attr(mu.qp,'n.wet')> 1000")   
  }
  if (!silent) print(paste("Number of locations AFTER weeding=",length(mu.qp),
               " and ",sum(attr(mu.qp,'n.wet')==0),"locations had",
         "'n.wet' set to zero"))

  # Last check:
  lons <- round(attr(mu.qp,"longitude"),4)
  lats <- round(attr(mu.qp,"latitude"),4)
  if (!silent) print(paste("Check: should be zero: ",sum(table(lons,lats)>1)))

  invisible(mu.qp)
}

dist2coast <- function(mu.qp) {
  print("dis2coast takes a bit of time to complete, but the")
  print("process will save intermediate results in a temporary")
  print("data file, so that the task can be carried out in")
  N <- length(attr(mu.qp,'mu')); d <- rep(NA,N)
  print(paste("several steps. N=",N))
  if (file.exists("dist2coast.temp.rda")) {
    load("dist2coast.temp.rda")
    print(paste("Using the unfinished results from previous call: i1=",i1))
  } else {
    i1 <- 1
  }
  data(addland1,envir = environment())
  for (i in i1:N) {
     d[i]<- round(min(distAB(attr(mu.qp,"longitude")[i],
            attr(mu.qp,"latitude")[i],lon.cont,lat.cont),na.rm=TRUE)/1000)
     if (!is.finite(d[i])) d[i] <- NA
     if (mod(i,10)==0) cat(".")
     if (mod(i,100)==0) cat(d[i])
     if (mod(i,1000)==0) {
       cat("o")
       i1 <- i + 1
       save(file="dist2coast.temp.rda",i1,d)
     }
  }
  file.remove("dist2coast.temp.rda")
  attr(mu.qp,"dist2coast.km") <- d
  invisible(mu.qp)
}


pointdensity <- function(x1,x2,resol=NULL) {
  # Creates a density map of from a scatter plot.

  #require( LatticeKrig)
  good <- is.finite(x1) & is.finite(x2)
  x1 <- x1[good]; x2 <- x2[good]
  n <- length(x1); nn <- max(c(10,round(n/100)))
  
  if (is.null(resol)) {
    resol <- 0.5* ( max(x1) - min(x1) )/nn +
             0.5* ( max(x2) - min(x2) )/nn
  }
  x1 <- round(x1/resol)
  x2 <- round(x2/resol)
  freq <- table(x1,x2); Y <- as.matrix(freq)
  x <- as.numeric(rownames(freq))*resol
  y <- as.numeric(colnames(freq))*resol
  X1 <- min(x); X2 <- max(x)
  Y1 <- min(y); Y2 <- max(y)
  
  Y[Y < 0] <- 0
  points <- data.frame(x=rep(x,length(y)),y=sort(rep(y,length(x))),
                           z=c(Y))
  
  #x11(); plot(points$x,points$y,pch=".")
  #text(points$x,points$y,round(points$z/100,1),cex=0.6,col="grey")
  #contour(x,y,Y,add=TRUE,col="red")
  
  pdens <- with(points, interp(x, y, z,
                               xo=seq(X1,X2, length=1000),
                               yo=seq(Y1,Y2, length=1000)))
    
  #contour(pdens$x,pdens$y,pdens$z,add=TRUE)

  invisible(pdens)
}





qplotGDCN <- function(x.0=1,Pmax=180,
                      prs=c(seq(0.50,0.95,by=0.05),0.96,0.97,0.98,0.99),
                      max.stations=NULL,xlim=c(0,200),ylim=c(0,200),
                      months=NULL,years=NULL,N.min=NULL,path="~/GDCN/") {
  # Main routine reading in, analysing, plotting, and storing the
  # quantiles from the GDCN data set.
  require(PrecipStat)
  
  list <- list.files(path,pattern=".dly",full.names = TRUE)
  data(gdcn.inv,envir = environment())
  
  N <- length(list); np <- length(prs)
  if (!is.null(max.stations)) N <- max.stations
  meanP <- rep(NA,N); mu <- meanP; country <- mu
  n.dry <- meanP; n.wet <- meanP; sigma <- meanP
  lat <- mu; lon <- mu; alt <- mu; stnr <- mu; rmse <- mu
  mu.qp <- matrix(rep(NA,N*np),np,N)
  
  plot(xlim,ylim,type="l",
     main="Wet-day 24-hr precipitation from GDCN",col="grey",
     xlab="-log(1-p)*mu [mm/day]",ylab="quantile(X,p) [mm/day]",lwd=3,
     sub=paste("thresh.= ",x.0,"mm/day; #stations=",N,
       " [qqplotter.R]"))
  text(0,ylim[2]*0.975,"Similarity with exponential distribution &",pos=4)
  text(0,ylim[2]*0.950,"connection between mean (mu) and percentiles (q_p)",
       pos=4,cex=0.75)
  grid()

  # Check if there are unfinished results
  if (file.exists("readGDCN.temp.rda")) {
    print("read intermediate results from previous run")
    load("readGDCN.temp.rda")
  } else {
    i1 <- 1
  }
  
  # Loop through all station files

  for (i in i1:N) {
    X <- readGDCN(list[i])
    imatch <- (sum(is.element(gdcn.inv$stnr,attr(X,"Station_number"))>0)) &
              (sum(X >= x.0,na.rm=TRUE) > 100)
    if (imatch) {
      meanP[i] <- round(mean(X,na.rm=TRUE),2)
      n.dry[i] <- sum(X < x.0,na.rm=TRUE)
      n.wet[i] <- sum(X >= x.0,na.rm=TRUE)
      if (n.wet[i]>100)
        sigma[i] <- round(var(c(X[X >= x.0]),na.rm=TRUE),2)
      lat[i] <- gdcn.inv$lat[imatch]
      lon[i] <- gdcn.inv$lon[imatch]
      alt[i] <- gdcn.inv$alt[imatch]
      stnr[i] <- attr(X,"Station_number")
      country[i] <- attr(X,"Country_code")
    
      # If specified, only include given months
      # N.B. There was a bug in versions 0.9-1.03 where the routine
      # selected wrong months. Bug corrected in v. 1.04 and mu.qp.jjas
      # and mu.qp.world.jjas have been re-computed.
      if (!is.null(months)) X <- X[is.element(attr(X,"month"),months),]
      if (!is.null(years)) X <- X[is.element(attr(X,"year"),years),]

      # If specified, only include stations with minimum valid data
      doit <- TRUE
      if (!is.null(N.min))
        if (sum(is.finite(X)) < N.min) doit <- FALSE
    } else doit <- FALSE
    if (doit) {
      # Apply the analysis:
      quants <- fitsimple(c(X),x.0,prs)
      mu[i] <- round(attr(quants,"mu"))
      rmse[i] <- sqrt(sum((quants$x2 - quants$x1)^2))/length(quants$x2)
      # Plot the results:
      col <- rgb(max(0,1-(meanP[i]/Pmax)^0.5),0, min(1,(meanP[i]/Pmax)^0.5))
      points(quants,pch=19,cex=0.6,col=col)
      print(paste(i,length(list),list[i],country[i],stnr[i],lat[i],
                  lon[i],alt[i],mu[i],meanP[i]))
        # Store the main results
      mu.qp[,i] <- quants$y
    } else stnr[i] <- NA

    if (mod(i,1000)==0) {
      # Intermediate save:
      i1 <- i+1
      save(file="readGDCN.temp.rda",i1,mu,rmse,mu.qp,stnr,meanP,
           n.dry,n.wet,sigma,lat,lon,alt,country)
    }
  } # End-of-loop: station files
  file.remove("readGDCN.temp.rda")
      
  attr(mu.qp,"description") <- "quantile(X,p) [mm/day]"
  attr(mu.qp,"mean_precip") <- meanP
  attr(mu.qp,"mu") <- mu
  attr(mu.qp,"longitude") <- lon
  attr(mu.qp,"latitude") <- lat
  attr(mu.qp,"station_number") <- stnr
  attr(mu.qp,"altitude") <- alt
  attr(mu.qp,"n.wet") <- n.wet
  attr(mu.qp,"n.dry") <- n.dry
  attr(mu.qp,"sigma") <- sigma
  attr(mu.qp,"rmse") <- rmse
  attr(mu.qp,"country") <- country
  attr(mu.qp,"history") <- "qqplotter.R - qplotGDCN()"
  attr(mu.qp,"probabilities") <- prs
  attr(mu.qp,"months") <- months
  attr(mu.qp,"formulae") <- "-log(1-p)*mu [mm/day]"
  attr(mu.qp,"x.0") <- x.0
  attr(mu.qp,"source") <- "GDCN"
  class(mu.qp) <- "qqplotter"
  #Save the main results:
  save(file="mu.qp.rda",mu.qp)

  Xn <- estquantiles(mu.qp)
  pd <- pointdensity(Xn,mu.qp)
  contour(pd$x,pd$y,pd$z,add=TRUE,col="steelblue",lev=seq(25,20000,by=25))

  qqfit <- lm(c(mu.qp) ~ c(Xn))
  print(summary(qqfit))
  coef <- round(summary(qqfit)$coefficients,3)
  text(550,50,paste("y = ",coef[1],"(+-",coef[3],") + ",coef[2],"(+-",
                    coef[4],") x"),cex=0.75,pos=2)
  abline(qqfit,col="red")
  lines(c(0,500),c(0,500),col="grey")

  # Add colour bar and legend
  colkey1 <- rgb(1-(0:Pmax/Pmax)^05,rep(0,Pmax+1),(0:Pmax/Pmax)^0.5)
  text(45,550,"mean P",cex=0.75)
  text(45,520,"wet",cex=0.5,col="grey")
  text(45,360,"dry",cex=0.5,col="grey")
  colorbar(colkey1,c(0.17,0.23,0.6,0.75))

  # Saved the graphics"
  #dev2bitmap(file="qqplotter.pdf",type="pdfwrite")
  #dev2bitmap(file="qqplotter.png",res=150)
  #dev.copy2eps(file="qqplotter.eps")

  invisible(mu.qp)
}






qplotESCN <- function(x.0=1,Pmax=180,
                      prs=c(seq(0.50,0.95,by=0.05),0.96,0.97,0.98,0.99),
                      max.stations=NULL,xlim=c(0,200),ylim=c(0,200),
                      months=NULL,years=NULL,N.min=NULL) {
  # Main routine reading in, analysing, plotting, and storing the
  # quantiles from the ESCN ECA&D data set.
  locs <- read.table("data.eca/stations.txt", skip=17,sep = ",",
                     header = TRUE)
  locs$LAT <- as.character(locs$LAT)
  locs$LON <- as.character(locs$LON)
  locs$LAT <- as.numeric(substr(locs$LAT,1,3)) +
              as.numeric(substr(locs$LAT,5,6))/60 +
              as.numeric(substr(locs$LAT,8,9))/3600
  locs$LON <- as.numeric(substr(locs$LON,1,3)) +
              as.numeric(substr(locs$LON,5,6))/60 +
              as.numeric(substr(locs$LON,8,9))/3600  
  list <- list.files(path="data.eca",pattern="RR_",full.names=TRUE)
  N <- length(list); np <- length(prs)
  if (!is.null(max.stations)) N <- max.stations
  meanP <- rep(NA,N); mu <- meanP; country <- mu
  n.dry <- meanP; n.wet <- meanP; sigma <- meanP
  lat <- mu; lon <- mu; alt <- mu; stnr <- mu; rmse <- mu
  mu.qp <- matrix(rep(NA,N*np),np,N)
  
  plot(xlim,ylim,type="l",
     main="Wet-day 24-hr precipitation from ECA&D",col="grey",
     xlab="-log(1-p)*mu [mm/day]",ylab="quantile(X,p) [mm/day]",lwd=3,
     sub=paste("thresh.= ",x.0,"mm/day; #stations=",N,
       " [qqplotter.R]"))
  text(0,ylim[2]*0.975,"Similarity with exponential distribution &",pos=4)
  text(0,ylim[2]*0.950,"connection between mean (mu) and percentiles (q_p)",
       pos=4,cex=0.75)
  grid()

  # Check if there are unfinished results
  if (file.exists("readESCN.temp.rda")) {
    print("read intermediate results from previous run")
    load("readGDCN.temp.rda")
  } else {
    i1 <- 1
  }
  # Loop through all station files

  for (i in i1:N) {
    escn <- read.table(list[i], skip=20,sep = ",",header=TRUE)
    X <- 0.1*escn$RR
    uscr <- instring("D",list[i]); dot <- instring(".",list[i])
    #print(c(uscr,dot[length(dot)],nchar(list[i])))
    #print(substr(list[i],uscr+1,dot[length(dot)]-1))
    #print(list[i])
    stnr[i] <- as.numeric(substr(list[i],uscr+1,dot[length(dot)]-1))
    imatch <- is.element(locs$STAID,stnr[i])
    if (sum(imatch)>0) {
      lat[i] <- round(locs$LAT[imatch],4)
      lon[i] <- round(locs$LON[imatch],4)
      alt[i] <- locs$HGHT[imatch]
      country[i] <- locs$CN[imatch]
    }
    escn$month <- as.numeric(substr(as.character(escn$DATE),5,6))
    X[X<0] <- NA
    meanP[i] <- round(mean(X,na.rm=TRUE),2)
    n.dry[i] <- sum(X < x.0,na.rm=TRUE)
    n.wet[i] <- sum(X >= x.0,na.rm=TRUE)
    if (n.wet[i]>100) sigma[i] <- round(var(X[X >= x.0],na.rm=TRUE),2)

    # If specified, only include given months
    if (!is.null(months)) X <- X[is.element(escn$month,months)]
    if (!is.null(years)) X <- X[is.element(escn$year,years)]

    # If specified, only include stations with minimum valid data
    doit <- TRUE
    if (!is.null(N.min))
      if (sum(is.finite(X)) < N.min) doit <- FALSE
    if (doit) {
      # Apply the analysis:
      quants <- fitsimple(X,x.0,prs)
      mu[i] <- round(attr(quants,"mu"),4)
      rmse[i] <- sqrt(sum(quants$x2 - quants$x1))
      
      # Plot the results:
      col <- rgb(max(0,1-(meanP[i]/Pmax)^0.5),0, min(1,(meanP[i]/Pmax)^0.5))
      points(quants,pch=19,cex=0.6,col=col)
      print(paste(i,length(list),list[i],country[i],stnr[i],lat[i],
                  lon[i],alt[i],mu[i],meanP[i]))
      # Store the main results
      mu.qp[,i] <- quants$y

      if (mod(i,1000)==0) {
      # Intermediate save:
      i1 <- i+1
      save(file="readESCN.temp.rda",i1,mu,rmse,mu.qp,stnr,
           meanP,n.dry,n.wet,sigma,lat,lon,alt,country)
    }
    }
  } # End-of-loop: station files
  file.remove("readGDCN.temp.rda")
      
  attr(mu.qp,"description") <- "quantile(X,p) [mm/day]"
  attr(mu.qp,"mean_precip") <- meanP
  attr(mu.qp,"mu") <- mu
  attr(mu.qp,"n.wet") <- n.wet
  attr(mu.qp,"n.dry") <- n.dry
  attr(mu.qp,"sigma") <- sigma
  attr(mu.qp,"rmse") <- rmse
  attr(mu.qp,"longitude") <- lon
  attr(mu.qp,"latitude") <- lat
  attr(mu.qp,"station_number") <- stnr
  attr(mu.qp,"altitude") <- alt
  attr(mu.qp,"country") <- country
  attr(mu.qp,"history") <- "qqplotter.R - qplotESCN()"
  attr(mu.qp,"probabilities") <- prs
  attr(mu.qp,"months") <- months
  attr(mu.qp,"formulae") <- "-log(1-p)*mu [mm/day]"
  attr(mu.qp,"x.0") <- x.0
  attr(mu.qp,"source") <- "ECA&D"
  class(mu.qp) <- "qqplotter"
  #Save the main results:
  mu.qp.escn <- mu.qp
  save(file="mu.qp.escn.rda",mu.qp.escn)

  Xn <- estquantiles(mu.qp.escn)
  pd <- pointdensity(Xn,mu.qp.escn)
  contour(pd$x,pd$y,pd$z,add=TRUE,col="steelblue",lev=seq(25,20000,by=25))
  
  qqfit <- lm(c(mu.qp.escn) ~ c(Xn))
  #qqfit <- lm(x2 ~ x1 + I(x1^2))
  print(summary(qqfit))
  coef <- round(summary(qqfit)$coefficients,3)
  text(550,50,paste("y = ",coef[1],"(+-",coef[3],") + ",coef[2],"(+-",
                    coef[4],") x"),cex=0.75,pos=2)
  abline(qqfit,col="red")
  lines(c(0,500),c(0,500),col="grey")

  # Add colour bar and legend
  colkey1 <- rgb(1-(0:Pmax/Pmax)^05,rep(0,Pmax+1),(0:Pmax/Pmax)^0.5)
  text(45,550,"mean P",cex=0.75)
  text(45,520,"wet",cex=0.5,col="grey")
  text(45,360,"dry",cex=0.5,col="grey")
  colorbar(colkey1,c(0.17,0.23,0.6,0.75))

  # Saved the graphics"
  #dev2bitmap(file="qqplotter4.pdf",type="pdfwrite")
  #dev2bitmap(file="qqplotter4.png",res=150)
  #dev.copy2eps(file="qqplotter4.eps")

  invisible(mu.qp.escn)
}



Sum.is.zero <- function(X) sum(X==0,na.rm=TRUE)


Sum.gt.zero <- function(X) sum(X>0,na.rm=TRUE)


qplotRCM <- function(mu.qp=NULL,x.0=1,sim.rng=TRUE,
                     prs=c(seq(0.50,0.95,by=0.05),0.96,0.97,0.98,0.99),
                     do.narccap=FALSE, do.ensembles=TRUE,Max.pts=153569,
                     xlim=c(0,200),ylim=c(0,200),months=NULL,
                     path.narccap="~/NARCCAP/data/",
                     path.ensembles="~/ENSEMBLES.daily.50km") {
    # Main routine reading in, analysing, plotting, and storing the
    # quantiles from RCMs.
    #ENSEMBLES & NARCCAP
  require(clim.pact)
  require(ncdf)
  require(PrecipStat)

  if (is.null(mu.qp)) data(mu.qp,envir = environment())
  print(paste("qplotRCM(mu.qp,",x.0,",",sim.rng,", prs,",do.narccap,",",
              do.ensembles,",",Max.pts,")"))
  # Set up plot with main results from observations (qplotGDCN)
  if (is.null(mu.qp)) {
    print("Results from GCDN not available")
    mu.qp <- list(x1=c(0,600),x2=c(0,600)); Xn <- mu.qp
    type <- "l"
    N <- 0
  } else {
    if (class(mu.qp)!="qqplotter") stop("qplotRCM: need 'qqplotter' object")
    type="p"
    # sim.rng = TRUE for comparing with stations with similar range of
    # 24-hr precipitation
    D <- dim(mu.qp)
    Xn <- estquantiles(mu.qp)
    print(D)
    N <- D[2]
    if (sim.rng) {
      print("Only include GCDN data with 95-percentlie less than 200mm/day")
      ikeep <- mu.qp[D[1],] < 200
      ikeep[is.na(ikeep)] <- FALSE
      Xn <- Xn[,ikeep]
      mu.qp <- mu.qp[,ikeep]
      print(dim(mu.qp))
      N <- sum(ikeep)
      print(paste("Keeping",N,"stations"))
    }
  }

  dev.new()
  par(las=1)
  plot(c(Xn),c(mu.qp),type=type,pch=19,col="grey75",
       main="Wet-day 24-hr precipitation from NARCCAP/ENSEMBLES",
       xlab="-log(1-p)*mu [mm/day]",ylab="quantile(X,p) [mm/day]",lwd=3,
       sub=paste("thresh.= ",x.0,"mm/day; [qqplotter.R]"),
       xlim=xlim,ylim=ylim)
  lines(c(0,1600),c(0,1600),col="grey40")
  text(0,ylim[2]*0.975,"Similarity with exponential distribution &",pos=4)
  text(0,ylim[2]*0.950,"connection between mean (mu) and percentiles (q_p)",
       pos=4,cex=0.75)
  grid()

  N <- Max.pts; np <- length(prs)
  meanP <- rep(NA,N); mu <- meanP
  n.dry <- meanP; n.wet <- meanP; sigma <- meanP
  lat <- mu; lon <- mu; alt <- mu; rmse <- mu
  expnt <- rep("NA",N); rcmnm <- expnt
  mu.qp.rcm<- matrix(rep(NA,N*np),np,N)
  iv <- 0
  
  if (do.narccap) {
    
  # Search for NARCCAP files:
  narccap <- list.files(path=path.narccap,
                       pattern="_ncep",full.names = TRUE)
  rcm.name <- list.files(path=path.narccap,
                     pattern="_ncep")
  rcm.name <- substr(rcm.name,1,nchar(rcm.name)-5)

  # Loop through NARCCAP files:
  print("NARCCAP results:")
  print("This part is not yet ready - The NARCCAP data was stored as monthly")
  n2 <- length(narccap)
  for (i in 1:n2) {
    print(narccap[i])
    load(narccap[i])
    eval(parse(text=paste("rcm <- ",rcm.name[i])))
    eval(parse(text=paste("rm(",rcm.name[i],")")));  gc(reset=TRUE)        

    if (!exists("rcm$months")) rcm$months <- rcm$months.c
    if (!is.null(months)) it <- is.element(rcm$months,months) else
                          it <- is.finite(rcm$months)
      
    # The fields do not have identical names...
    if (length(grep("pr.c",names(rcm)))>0)
       d <- dim(rcm$pr.c) else
    if (length(grep("pr",names(rcm)))>0)
       d <- dim(rcm$pr)
   
    if (length(grep("pr.c",names(rcm)))>0)
          X <- rcm$pr.c else
        if (length(grep("pr",names(rcm)))>0)
          X <- rcm$pr else X <- NA
    rm("rcm"); gc(reset=TRUE)
    print("Extracting some, removing the rest - clearing space")
    dim(X) <- c(d[1]*d[2],d[3])
    print(dim(X))
    print("Call 'fitsimple'")
    quants.rcm <- fitsimple(X,x.0,prs)
    points(quants.rcm,pch=19,cex=0.5,col="red")
    iv <- max(iv) + 1:(d[1]*d[2])
    print("Get and keep the stats.")
    mu.qp.rcm[,iv] <- quants.rcm$y
    meanP[iv] <- round( apply(X,1,mean,na.rm=TRUE),2 )
    X[X < x.0] <- 0
    n.wet[iv] <- apply(XX,1,Sum.gt.zero)
    n.dry[iv] <- apply(XX,1,Sum.is.zero)
    X[X < x.0] <- NA
    mu[iv] <-    round( apply(X,1,mean,na.rm=TRUE),4 )
    sigma[iv] <- round( apply(X,1,var,na.rm=TRUE),2 )
    expnt[iv] <- "NARCCAP"
    rcmnm[iv] <- rcm.name[i]
    rmse[iv] <- sqrt(sum((quants.rcm$x2-quants.rcm$x1)^2,na.rm=TRUE))/
                          length(quants.rcm$x2)
    print("Do the next RCM")
  } # end-loop NARCCAP files

  }

  # ENSEMBLES
  if (do.ensembles) {
    
  print("Read the ENSEMBLES data")
  ensembles <- list.files(path=path.ensembles,pattern="pr",
                          full.names = TRUE)
  rcm.name <- list.files(path=path.ensembles,pattern="pr")
  # Loop through ENSEMBLES files:
  n3 <- length(ensembles)
  for (i in 1:n3) {
    print(ensembles[i])
    uscr <- instring("_",rcm.name[i])
    rcm.name[i] <- substr(rcm.name[i],1,uscr-1)
    rcm <- readENSEMBLES(ensembles[i])

    if (!is.null(months))
      it <- is.element(attr(rcm,"date")$month,months) else
      it <- is.finite(attr(rcm,"date")$month)
      
    d <- dim(rcm)
    # loop through each grid-point. The unit is "kg m-2 s-1".
    for (ix in 1:d[1]) {
      for (iy in 1:d[2]) {
        X <- rcm[ix,iy,it]*3600*24
        # Apply analysis
        quants.rcm <- fitsimple(X,x.0,prs)
        # Plot results
        points(quants.rcm,pch=4,cex=0.6,col="blue")
        # Store the main results
        iv <- iv+1
        if (iv > Max.pts)
          print(paste("Max.pts=",Max.pts,"iv=",iv)) else {
          mu.qp.rcm[,iv] <- quants.rcm$y
          meanP[iv] <- round(mean(X,na.rm=TRUE),2)
          n.dry[iv] <- sum(X < x.0,na.rm=TRUE)
          n.wet[iv] <- sum(X >= x.0,na.rm=TRUE)
          sigma[iv] <- round(var(X[X >= x.0],na.rm=TRUE),2)
          expnt[iv] <- "ENSEMBLES"
          mu[iv] <- round(attr(quants.rcm,"mu"),4)
          rcmnm[iv] <- rcm.name[i]
          rmse[iv] <- sqrt(sum((quants.rcm$x2-quants.rcm$x1)^2,na.rm=TRUE))/
                           length(quants.rcm$x2)

        }
      }
    }
    # Clear the space
  } # end-loop ENSEMBLES files
  }

  mu.qp.rcm <- mu.qp.rcm[,1:iv]

  attr(mu.qp.rcm,"mean_precip") <- meanP[1:iv]
  attr(mu.qp.rcm,"mu") <- mu[1:iv]
  attr(mu.qp.rcm,"n.wet") <- n.wet[1:iv]
  attr(mu.qp.rcm,"n.dry") <- n.dry[1:iv]
  attr(mu.qp.rcm,"sigma") <- sigma[1:iv]
  attr(mu.qp.rcm,"history") <- "qqplotter.R - qplotRCM()"
  attr(mu.qp.rcm,"probabilities") <- prs
  attr(mu.qp.rcm,"formulae") <- "-log(1-p)*mu [mm/day]"
  attr(mu.qp,"x.0") <- x.0
  attr(mu.qp.rcm,"description") <- "quantile(X,p) [mm/day]"
  attr(mu.qp.rcm,"experiment") <- expnt[1:iv]
  attr(mu.qp.rcm,"rcm") <- rcmnm[1:iv]
  attr(mu.qp.rcm,"source") <-
    c("NARCCAP","ENSEMBLES")[c(do.narccap,do.ensembles)]
  attr(mu.qp.rcm,"rmse") <- rmse[1:iv]
  attr(mu.qp.rcm,"history") <- "qqplotter.R - qplotRCM()"
  class(mu.qp.rcm) <- "qqplotter"
  #Save the main results:
  save(file="mu.qp.rcm.rda",mu.qp.rcm)

  print(table(attr(mu.qp.rcm,"experiment")))
  narccap <- is.element(as.character(attr(mu.qp.rcm,"experiment")),
                                     "NARCCAP")
  ensembles <- is.element(as.character(attr(mu.qp.rcm,"experiment")),
                                     "ENSEMBLES")
  Xn <- estquantiles(mu.qp.rcm)
  if ( (do.narccap) & (sum(narccap)>100) ) {
    pd1 <- pointdensity(c(Xn[,narccap]),c(mu.qp.rcm[,narccap]),resol=1)
    contour(pd1$x,pd1$y,pd1$z,add=TRUE,col="darkred",
            lev=seq(25,20000,by=5000))
  }
  if ( (do.ensembles) & (sum(ensembles)>100) ) {
    print(summary(c(Xn[,ensembles])))
    print(summary(c(mu.qp.rcm[,ensembles])))
    pd2 <- pointdensity(c(Xn[,ensembles]),c(mu.qp.rcm[,ensembles]),resol=1)
    contour(pd2$x,pd2$y,pd2$z,add=TRUE,col="lightblue",
            lev=seq(25,20000,by=5000))
  }
  lines(c(0,900),c(0,900),col="grey",lty=2)
  grid()
  
  # Saved the graphics"
  #dev2bitmap(file="qqplotter2.pdf",type="pdfwrite")
  #dev2bitmap(file="qqplotter2.png",res=150)
  #dev.copy2eps(file="qqplotter2.eps")

  invisible(mu.qp.rcm)
}





qplotres <- function(mu.qp=NULL,x.0=1,sim.rng=TRUE,
                     prs=c(seq(0.50,0.95,by=0.05),0.96,0.97,0.98,0.99),
                     xlim=c(0,100),ylim=c(0,100),
                     Max.pts=90000,path.25km="~/ENSEMBLES.daily",
                     path.50km="~/ENSEMBLES.daily.50km") {
  # Main routine reading in, analysing, plotting, and storing the
  # quantiles from RCMs with different spatial resolution.
  
  require(ncdf)
  require(PrecipStat)
  
  if (is.null(mu.qp)) data(mu.qp,envir = environment())

  # Set up plot with main results from observations (qplotGDCN)
  if (is.null(mu.qp)) {
    print("Results from GCDN not available")
    mu.qp <- list(x1=c(0,600),x2=c(0,600))
    type <- "l"
    N <- 0
  } else {
    if (class(mu.qp)!="qqplotter") stop("qplotres: need 'qqplotter' object")
    type="p"
    # sim.rng = TRUE for comparing with stations with similar range of
    # 24-hr precipitation
    D <- dim(mu.qp)
    Xn <- estquantiles(mu.qp)
    #print(D)
    N <- D[2]
    if (sim.rng) {
      print("Only include GCDN data with 95-percentlie less than 300mm/day")
      mu.qp <- qweed(mu.qp,"mu.qp[14,] < 300")
      print(dim(mu.qp))
      N <- dim(mu.qp)[2]
      print(paste("Keeping",N,"stations"))
    }
  }
  Xn <- estquantiles(mu.qp)
  dev.new()
  plot(c(Xn),c(mu.qp),type=type,pch=19,col="grey75",
       main="Wet-day 24-hr precipitation from ENSEMBLES 25-50km resolution",
       xlab="-log(1-p)*mu [mm/day]",ylab="quantile(X,p) [mm/day]",lwd=3,
       sub=paste("thresh.= ",x.0,"mm/day; #stations=",N,
           " [qqplotter.R]"),xlim=xlim,ylim=ylim)
  lines(c(0,600),c(0,600),col="grey40")
  text(0,ylim[2]*0.975,"Similarity with exponential distribution &",pos=4)
  text(0,ylim[2]*0.950,"connection between mean (mu) and percentiles (q_p)",
       pos=4,cex=0.75)
  grid()

  N <- Max.pts; np <- length(prs)
  meanP <- rep(NA,N); mu <- meanP;
  n.dry <- meanP; n.wet <- meanP; sigma <- meanP
  lat <- mu; lon <- mu; alt <- mu; rmse <- mu
  mu.qp.res <- matrix(rep(NA,N*np),np,N); expnt <- rep("NA",N)
  iv <- 0
  
  # ENSEMBLES
    
  print("Read the 50-km ENSEMBLES data")
  ensembles <- list.files(path=path.50km,pattern="HadRM",
                          full.names = TRUE)
  # Loop through ENSEMBLES files:
  n3 <- length(ensembles)
  for (i in 1:n3) {
    print(ensembles[i])
    ncid <- open.ncdf(ensembles[i])
    lons <- get.var.ncdf(ncid,"rlon"); nx <- length(lons)
    lats <- get.var.ncdf(ncid,"rlat"); ny <- length(lats)
    time <- get.var.ncdf(ncid,"time"); nt <- length(time)
    #The total interval is 1950:2100
    rcm <- get.var.ncdf(ncid,"pr",
                        start=c(round(nx/4),round(ny/4),1),
                        count=c(round(3*nx/4),round(3*ny/4),nt))
    d <- dim(rcm)
    # loop through each grid-point
    for (ix in 1:d[1]) {
      for (iy in 1:d[2]) {
        X <- rcm[ix,iy,]*3600*24
        # Apply analysis
        quants.rcm <- fitsimple(X,x.0,prs)
        # Plot results
        points(quants.rcm,pch=19,cex=0.6,col="darkblue")
        # Store the main results
        iv <- iv+1
        mu.qp.res[,iv] <- quants.rcm$y
        meanP[i] <- round(mean(X,na.rm=TRUE),2)
        n.dry[i] <- sum(X < x.0,na.rm=TRUE)
        n.wet[i] <- sum(X >= x.0,na.rm=TRUE)
        sigma[i] <- round(var(X[X >= x.0],na.rm=TRUE),2)
        expnt[iv] <- "lowres"
        mu[iv] <- round(attr(quants.rcm,"mu"))
        rmse[iv] <- sqrt(sum((quants.rcm$x2-quants.rcm$x1)^2,na.rm=TRUE))/
                        length(quants.rcm$x2)

      }
    }
    # Clear the space
  } # end-loop ENSEMBLES files

  print("Read the 25-km ENSEMBLES data")
  ensembles <- list.files(path=path.25km,pattern="HadRM",
                          full.names = TRUE)
  # Loop through ENSEMBLES files:
  n3 <- length(ensembles)
  for (i in 1:n3) {
    print(ensembles[i])
    ncid <- open.ncdf(ensembles[i])
    lons <- get.var.ncdf(ncid,"rlon"); nx <- length(lons)
    lats <- get.var.ncdf(ncid,"rlat"); ny <- length(lats)
    time <- get.var.ncdf(ncid,"time"); nt <- length(time)
    #The total interval is 1950:2100
    rcm <- get.var.ncdf(ncid,"pr",
                        start=c(round(nx/4),round(ny/4),1),
                        count=c(round(3*nx/4),round(3*ny/4),nt))
    d <- dim(rcm)
    # loop through each grid-point
    
    for (ix in 1:d[1]) {
      for (iy in 1:d[2]) {
        X <- rcm[ix,iy,]*3600*24
        # Apply analysis
        quants.rcm <- fitsimple(X,x.0,prs)
        # Plot results
        points(quants.rcm,pch=4,cex=0.5,col="lightblue")
        # Store the main results
        iv <- iv+1
        mu.qp.res[,iv] <- quants.rcm$y
        meanP[i] <- round(mean(X,na.rm=TRUE),2)
        n.dry[i] <- sum(X < x.0,na.rm=TRUE)
        n.wet[i] <- sum(X >= x.0,na.rm=TRUE)
        sigma[i] <- round(var(X[X >= x.0],na.rm=TRUE),2)
        expnt[iv] <- "highres"
        mu[iv] <- round(attr(quants.rcm,"mu"))
        rmse[iv] <- sqrt(sum((quants.rcm$x2-quants.rcm$x1)^2,na.rm=TRUE))/
                        length(quants.rcm$x2)
      }
    }
    # Clear the space
  } # end-loop ENSEMBLES files
  print("Save the data")
  attr(mu.qp.res,"mean_precip") <- meanP
  attr(mu.qp.res,"mu") <- mu
  attr(mu.qp.res,"n.wet") <- n.wet
  attr(mu.qp.res,"n.dry") <- n.dry
  attr(mu.qp.res,"sigma") <- sigma
  attr(mu.qp.res,"rmse") <- rmse[1:iv]
  attr(mu.qp.res,"history") <- "qqplotter.R - qplotres()"
  attr(mu.qp.res,"probabilities") <- prs
  attr(mu.qp.res,"formulae") <- "-log(1-p)*mu [mm/day]"
  attr(mu.qp,"x.0") <- x.0
  attr(mu.qp.res,"description") <- "quantile(X,p) [mm/day]"
  attr(mu.qp.res,"experiment") <- expnt[1:iv]
  attr(mu.qp.res,"history") <- "qqplotter.R - qplotres()"
  class(mu.qp.res) <- "qqplotter"
  #Save the main results:
  save(file="mu.qp.res.rda",mu.qp.res)

  lowres <- is.element(as.character(attr(mu.qp.res,"experiment")),
                                     "lowres")
  Xn <- estquantiles(mu.qp)
  pd0 <- pointdensity(Xn,mu.qp,resol=3)
  contour(pd0$x,pd0$y,pd0$z,add=TRUE,lev=seq(25,20000,by=500))
  Xn <- estquantiles(mu.qp.res)
  pd1 <- pointdensity(Xn[,lowres],mu.qp.res[,lowres],resol=1)
  contour(pd1$x,pd1$y,pd1$z,add=TRUE,col="darkblue",lev=seq(25,20000,by=500))
  pd2 <- pointdensity(Xn[,!lowres],mu.qp.res[,!lowres],resol=1)
  contour(pd2$x,pd2$y,pd2$z,add=TRUE,col="steelblue",lev=seq(25,20000,by=500))
  lines(c(0,900),c(0,900),col="grey",lty=2)
  grid()
  
  # Saved the graphics"
  #dev2bitmap(file="qqplotter3.pdf",type="pdfwrite")
  #dev2bitmap(file="qqplotter3.png",res=150)
  #dev.copy2eps(file="qqplotter3.eps")

  invisible(mu.qp.res)
}






niceqqplot <- function(mu.qp,col.axis="black",addcont=TRUE,addfit=TRUE,
                       xlim=c(0,200),ylim=c(0,200),colleg=TRUE,
                       x.0=1,colourcoding="mean_precip",
                       CI=TRUE,N.min=NULL,p=NULL,add=FALSE,
                       col=NULL,pch=19,cex=1) {
  # Plot qqplots, given a qqplotter object mu.qp
  if (class(mu.qp)!="qqplotter") stop("niceqqplot: need 'qqplotter' object")
  if (!is.null(N.min)) mu.qp <- qweed(mu.qp,
                       crit=paste("attr(mu.qp,'n.wet')> ",N.min))
  N <- dim(mu.qp)[2]
  attr(mu.qp,colourcoding) -> crit
  crit[!is.finite(crit)] <- 0
  crit[crit > quantile(crit,0.99,na.rm=TRUE)] <- quantile(crit,0.99,na.rm=TRUE)
  attr(mu.qp,"mu") ->  mu
  if (!is.null(p))
    ip <- max( (1:dim(mu.qp)[1])[attr(mu.qp,"probabilities")<=p]) else
    ip <- (1:dim(mu.qp)[1])
  p <- attr(mu.qp,"probabilities")[ip]
  #print(attr(mu.qp,"probabilities"))

  src <- rownames(table(attr(mu.qp,"source")))
  if (length(src)>1) src <- paste(src[1],"&",src[2])
  if (!add) {
    if (col.axis=="white") {
      xaxt <- "n"; yaxt <- "n"
    } else {
      xaxt <- "s"; yaxt <- "s"
    }
    par(col.axis=col.axis,col.lab=col.axis,xaxt=xaxt,yaxt=yaxt,las=1)
    plot(xlim,ylim,type="l",
       main="Wet-day 24-hr precipitation",col="grey40",
       xlab="-log(1-p)*mu [mm/day]",ylab="quantile(X,p) [mm/day]",lwd=3,
       sub=paste(src,"thresh.= ",x.0,"mm/day; #stations=",N,
         " [qqplotter.R]"))
    text(0,ylim[2]*0.975,"Similarity with exponential distribution &",pos=4)
    text(0,ylim[2]*0.950,"connection between mean (mu) and percentiles (q_p)",
         pos=4,cex=0.75)
    if (!is.null(N.min)) {
      text(mean(xlim)*1.2,ylim[1]*1.2,paste("Stations with N.wet > ",
                                            N.min,"days"))
    }
    if ( (!is.null(p)) & (length(p)==1) ) {
      text(mean(xlim)*1.2,ylim[1]+10,paste(p," quantile"))
    }
    if (col.axis=="black") grid()
  }
  Xn <- estquantiles(mu.qp)
  #print(dim(Xn)); print(dim(mu.qp)); print(N)
  if (is.null(col)) colcode <- TRUE else
                    colcode <- FALSE
  ii <- 1:N
  # Plot the results:
  # print(summary(crit[ii]))
  if (colcode) {
      Pmax <- max(crit[ii],na.rm=TRUE) 
      col <- rgb( rep(1,N)-(crit[ii]/Pmax)^0.5,
                  rep(0,N),
                  (crit[ii]/Pmax)^0.5 )
      #if (sum(ip)>1) col <- t(matrix(rep(col,sum(ip)),N,sum(ip)))
      #print(rbind(crit[1:10],col[1:10]))
    }
  #print(c(length(Xn[ip,ii]),length(mu.qp[ip,ii]),length(col)))
  points(Xn[ip,ii],mu.qp[ip,ii],pch=pch,col=col,cex=cex)
  rmse <- sqrt(sum((Xn[,ii] - mu.qp[,ii])^2,na.rm=TRUE))/length(Xn[,ii])
  attr(mu.qp,"rmse") <- rmse

  # End-of-loop: station files

  if (addcont) {
    # Add contour of point density
    pd <- pointdensity(Xn[ip,ii],mu.qp[ip,ii])
    contour(pd$x,pd$y,pd$z,add=TRUE,col="grey80")
  }

  if (addfit) {
    # Add linear model fit to cloud of points.
    qqfit <- lm(c(mu.qp[ip,ii]) ~ c(Xn[ip,ii]))
    print(summary(qqfit))
    coef <- round(summary(qqfit)$coefficients,3)
    text(550,50,paste("y = ",coef[1],"(+-",coef[3],") + ",coef[2],"(+-",
                      coef[4],") x"),cex=0.75,pos=2)
    abline(qqfit,col="red")
  }

  if (CI) {
    # Plot confidence interval if sufficient number of data points:
    N <- quantile(attr(mu.qp,"n.wet"),0.1,na.rm=TRUE)
    if (N > 100) {
    mu.qp <- qweed(mu.qp,crit="attr(mu.qp,'n.wet')>100")
      print(paste("Monte-Carlo -> CI; mu=",
                  quantile(attr(mu.qp,"mu"),0.95,na.rm=TRUE)))
    
      ci <- qplotMC(mu=max(attr(mu.qp,"mu"),na.rm=TRUE),
                    N=N,prs=attr(mu.qp,"probabilities"),plot=FALSE)
      lines(c(0,ci$Xn),c(0,ci$lower),lty=2)
      lines(c(0,ci$Xn),c(0,ci$upper),lty=2)
    }
  }
  
  lines(c(0,500),c(0,500),col="grey40",lwd=3)

  if (colleg & colcode) {
    # Add colour bar and legend
    colkey1 <- rgb(1-(0:Pmax/Pmax)^05,rep(0,Pmax+1),(0:Pmax/Pmax)^0.5)
    text(0,145,colourcoding,cex=0.75,srt=90)
    text(15,175,round(max(crit[ii],na.rm=TRUE)),cex=0.75,col="grey20")
    text(15,120,round(min(crit[ii],na.rm=TRUE)),cex=0.75,col="grey20")
    colorbar(colkey1,c(0.17,0.23,0.6,0.75))
  }
  results <- list(mu.qp=mu.qp,Xn=Xn)
  invisible(results)
}






compare1q <- function(mu.qp1,mu.qp2,p=0.95,N.min=500,type="rcm",
                      xlim=c(0,100),ylim=c(0,100),
                      col=c("blue","red","darkgreen","steelblue",
                        "darkred","green","magenta","cyan","grey30",
                        "lightblue","pink","grey80","wheat","brown")) {
  # Plot qqplots for several qqplotter objects for comparison
  niceqqplot(mu.qp1,colleg=FALSE,p=p,addcont=FALSE,CI=FALSE,
             addfit=FALSE,N.min=N.min,xlim=xlim,ylim=ylim,
             col="grey30",cex=1.25)
  exptyp <- rownames(table(attr(mu.qp2,type)))
  print(exptyp)
  for (i in 1:length(exptyp)) {
    crit <- paste("attr(mu.qp,'",type,"')=='",exptyp[i],"'",sep="")
    mu.qp.1 <- qweed(mu.qp2,crit=crit)
    if (sum(attr(mu.qp.1,'n.wet')>N.min) > 100) {
      niceqqplot(mu.qp.1,colleg=FALSE,p=p,add=TRUE,col=col[i],pch=4,
                 addfit=FALSE,CI=FALSE,addcont=FALSE,N.min=N.min,cex=0.5)
    } else {
      print(summary(attr(mu.qp.1,'n.wet')))
      exptyp[i] <- "NA"
    }
  }
  col <- col[-grep("NA",exptyp)]
  exptyp <- exptyp[-grep("NA",exptyp)]
  legend(0.8*xlim[2],0.5*ylim[2],exptyp,col[1:length(exptyp)],
         pch=4,cex=0.7,bg="grey95")
}





qqplotter <- function(x.0=1,Pmax=180,do.gdcn=FALSE,do.escn=FALSE,
                      do.rcm=FALSE,do.res=FALSE,do.narccap=TRUE,
                      do.ensembles=TRUE,max.stations=NULL,
                      europe2=TRUE,
                      prs=c(seq(0.50,0.95,by=0.05),0.96,0.97,0.98,0.99)) {
# A 'master' routine that uses the above functions to carry out the
# qqplotter analysis from scratch.
require(fields)
require(clim.pact)
require(PrecipStat)

# Do the time-consuming analysis for GDCN?

print(paste("qqploter: wet-day threshold=",x.0,
      "max expected mean precip=",Pmax," for (colouring the plots)",
      "do GDCN=",do.gdcn))
print("For the GDCN data see http://www.ncdc.noaa.gov/oa/climate/research/gdcn/gdcn.html")
print("For the ENSEMBLES results see http://ensemblesrt3.dmi.dk/")


# Carry out the actual operations:
if (do.gdcn) qplotGDCN(x.0,Pmax,prs,max.stations) -> mu.qp

# Retrieve current version if exists, otherwise use pre-calculated version:
if (!exists("mu.qp") & file.exists("mu.qp.rda")) load("mu.qp.rda") else
                                         data(mu.qp,envir = environment())

if (europe2) {
  print("Include ESCN-stations: qcat GDCN & ECA&D:")
  data(mu.qp.escn,envir = environment())
  MU.qp <- qcat(mu.qp,mu.qp.escn)
} else MU.qp <- mu.qp

dev.new()
par(las=1)
Xn0 <- estquantiles(MU.qp)
# QQ-box-plot:
bplot.xy(x=c(Xn0), y=c(MU.qp), N = 10,
           main="Wet-day 24-hr precipitation from GDCN",
           xlab="-log(1-p)*mu [mm/day]",ylab="quantile(X,p) [mm/day]",lwd=3,
           sub=paste("thresh.= ",x.0,"mm/day;",
           " [qqplotter.R]"))
pd <- pointdensity(Xn0,MU.qp)
contour(pd$x,pd$y,log(pd$z)/log(10),add=TRUE,col="grey")
lines(c(0,900),c(0,900),col="grey",lty=2)
grid()
#dev.copy2eps(file="qqplotter1b.eps")

# Make a table of rain guages per country
if (!europe2) {
  print("Table of countries included:")
  country <- table(as.character(attr(MU.qp,"country")))
  country.codes <- read.fwf("http://www.ncdc.noaa.gov/oa/climate/research/gdcn/Appendix_A_V1_0.txt",widths=c(3,3,24),col.names=c("code","dash","country"))
  categories <- rownames(country)
  contries <- rep("NA",length(categories))
  for (i in 1:length(categories)) {
    ii <- grep(categories[i],as.character(country.codes$code))
    contries[i] <- as.character(country.codes$country[ii])
  }
  rownames(country) <- contries
  print(country)
}

# PCA of the cloud of points in the qq-plot:
print("Perform PCA on the cloud of quantiles:")
dev.new()
par(las=1)
qPCA(MU.qp) -> pca
exclude <- attr(pca,'exclude')
print(round(pca$var[1:10],1))
niceqqplot(MU.qp,col="grey",addcont=TRUE,addfit=FALSE,CI=FALSE)
mode1 <- data.frame(x=abs(pca$x1[,1]*pca$d[1]*0.02),
                    y=abs(pca$qp[,1]*pca$d[1]*0.02))
mode2 <- data.frame(x=abs(pca$x1[,2]*pca$d[2]*0.02),
                    y=abs(pca$qp[,2]*pca$d[2]*0.02))
lines(mode1,lwd=2)
lines(mode1+mode2,lwd=2,lty=2)
lines(mode1-mode2,lwd=2,lty=2)
print(summary(mode1))
#dev.copy2eps(file="qqplotter5b.eps")


dev.new()
par(las=1)
attr(MU.qp,"pc1") <- pca$v[,1]
# Remove the outliers for a clearer colour coding
attr(MU.qp,"pc1")[attr(MU.qp,"pc1")>quantile(attr(MU.qp,"pc1"),0.99)]<-NA
attr(MU.qp,"pc1")[attr(MU.qp,"pc1")<quantile(attr(MU.qp,"pc1"),0.01,
                                             na.rm=TRUE)]<-NA
mapofstations(MU.qp,colourcoding="pc1")

dev.new()
attr(MU.qp,"pc2") <- pca$v[,2]
attr(MU.qp,"pc2")[attr(MU.qp,"pc2")>quantile(attr(MU.qp,"pc2"),0.99)]<-NA
attr(MU.qp,"pc2")[attr(MU.qp,"pc2")<quantile(attr(MU.qp,"pc2"),0.01,
                                             na.rm=TRUE)]<-NA
mapofstations(MU.qp,colourcoding="pc2")

# Compare distributions of PCs:
print("Compare distributions of PCs")
dev.new()
par(las=1)
hist(pca$v[,1],breaks=seq(-0.3,0.3,length=100)) -> h1
hist(pca$v[,2],breaks=seq(-0.5,0.5,length=100)) -> h2
plot(h1$mid,h1$counts,xlim=c(-0.1,0.1),type="l",lwd=2,log="y",
     main="Histograms of PC loadings",
     sub="black=mode#1; blue= mode#2")
lines(h2$mid,h2$counts,col="blue")
grid()

# Regression analysis against distance from coast, mean precip, latitude:
print("Regression against geographical parameters")
print("x1 = mean precip; x2= altitude; x3=rnorm(N);")
print("x4-5 = latitude; x5-7=longitude;")
N <- length(pca$v[,1])
#d <- rep(NA,N)

# If there is information about distance to coast, use it,
# otherwise use random noise
if (length(attr(MU.qp,"dist2coast.km"))>0) {
  print("x3 = distance to coast")
  d <- attr(MU.qp,"dist2coast.km")[-exclude]
} else {
  print("x3 = random noise")
  d<-rnorm(N)
}

print(length(attr(MU.qp,"mean_precip")[-exclude]))

exclude <- (1:dim(MU.qp)[2])[!is.finite(colMeans(MU.qp))]
print("Regression results for leading PC:")
calibr1 <- data.frame(y=pca$v[,1],x1=attr(MU.qp,"mean_precip")[-exclude],
                      x2=attr(MU.qp,"altitude")[-exclude],x3=d,
                      x4=sin(pi/180*attr(MU.qp,"latitude"))[-exclude],
                      x5=cos(pi/180*attr(MU.qp,"latitude"))[-exclude],
                      x6=sin(pi/180*attr(MU.qp,"longitude")[-exclude]),
                      x7=cos(pi/180*attr(MU.qp,"longitude"))[-exclude])
print(summary(lm(y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7,data=calibr1)))

print("Regression results for second PC:")
calibr2 <- data.frame(y=pca$v[,2],x1=attr(MU.qp,"mean_precip")[-exclude],
                      x2=attr(MU.qp,"altitude")[-exclude],x3=d,
                      x4=sin(pi/180*attr(MU.qp,"latitude"))[-exclude],
                      x5=cos(pi/180*attr(MU.qp,"latitude"))[-exclude],
                      x6=sin(pi/180*attr(MU.qp,"longitude"))[-exclude],
                      x7=cos(pi/180*attr(MU.qp,"longitude"))[-exclude])
print(summary(lm(y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7,data=calibr2)))

# Demonstrate that most of the scatter can be described by
# the two leading modes:
print("Reconstruct the quantiles from leading PCS:")
dev.new()
par(las=1)
plot(c(0,200),c(0,200),type="l",col="grey",
     main="Variations in PCA modes #1 & 2",
     xlab="-log(1-p)*mu [mm/day]",ylab="quantile(X,p) [mm/day]",lwd=3,
           sub=paste("thresh.= ",x.0,"mm/day;",
           " [qqplotter.R]"))
grid()
Xn <- estquantiles(MU.qp)
points(Xn,MU.qp,pch=19,col="grey75")

expfit <- lm(c(MU.qp) ~ c(Xn))
print("Regression analysis percentiles for exponential distribution")
print(dim(Xn))
print(summary(expfit))

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

qpca1 <- data.frame(y=c(MU.qp[,-exclude]),
                    x=c(outer(pca$qp[,1]*pca$d[1],pca$v[,1])))
pcafit1 <- lm(y ~ x,data=qpca1)
print("PCA 1 fit to the cloud of points in the figure:")
print(summary(pcafit1))

qpca2 <- data.frame(y=c(MU.qp[,-exclude]),
                    x=c(outer(pca$qp[,1]*pca$d[1],pca$v[,1]))+
                      c(outer(pca$qp[,2]*pca$d[2],pca$v[,2])))
pcafit2 <- lm(y ~ x,data=qpca2)
print("PCA 1+2 fit to the cloud of points in the figure:")
print(summary(pcafit2))

i95 <- is.element(attr(MU.qp,"probabilities"),0.95)
q95pca1 <- data.frame(y=MU.qp[i95,-exclude],
                      x=pca$qp[i95,1]*pca$d[1]*pca$v[,1])
q95pca2 <- data.frame(y=MU.qp[i95,-exclude],
                      x=pca$qp[i95,1]*pca$d[1]*pca$v[,1]+
                        pca$qp[i95,2]*pca$d[2]*pca$v[,2])
print("Variance for 95th percentile for mode1:")
print(summary(lm(y ~x, data=q95pca1)))
print("Variance for 95th percentile for modes1+2:")
print(summary(lm(y ~x, data=q95pca2)))

i99 <- is.element(attr(MU.qp,"probabilities"),0.99)
q99pca1 <- data.frame(y=MU.qp[i99,-exclude],
                      x=pca$qp[i99,1]*pca$d[1]*pca$v[,1])
q99pca2 <- data.frame(y=MU.qp[i99,-exclude],
                      x=pca$qp[i99,1]*pca$d[1]*pca$v[,1]+
                        pca$qp[i99,2]*pca$d[2]*pca$v[,2])
print("Variance for 99th percentile for mode1:")
print(summary(lm(y ~x, data=q99pca1)))
print("Variance for 99th percentile for modes1+2:")
print(summary(lm(y ~x, data=q99pca2)))
stat99 <- summary(lm(y ~x, data=q99pca2))

dev.new()
par(las=1)
plot(q99pca2,
     main="99th percentile",
     xlab="f(PCA modes 1 + 2)",ylab="observed",
           sub=paste("thresh.= ",x.0,"mm/day;",
           " [qqplotter.R]"))
lines(c(0,200),c(0,200),lwd=2)
abline(lm(y ~x, data=q99pca2),col="red",lty=2)
text(20,600,paste("y=",round(stat99$coefficients[1],2),
                  "[+-",round(stat99$coefficients[3],2),"] + ",
                  round(stat99$coefficients[3],2),
                  "[+-",round(stat99$coefficients[4],2),"] x"),pos=4)
text(20,570,paste("R2=",100*round(stat99$r.squared,2),"%"),pos=4)

             
# Measure of difference:
print("Estimate scores describing bias from exponential:")
Xn <- estquantiles(MU.qp)
ip <- is.element(attr(MU.qp,"probabilities"),c(0.90,0.95,0.99))
ratio.score <- apply(Xn[ip,],2,mean,na.rm=TRUE)/
               apply(MU.qp[ip,],2,mean,na.rm=TRUE)
ratio.score[ratio.score>quantile(ratio.score,0.99,na.rm=TRUE)] <- NA
ratio.score[ratio.score<quantile(ratio.score,0.01,na.rm=TRUE)] <- NA
attr(MU.qp,"ratio.score") <-ratio.score

dev.new()
par(las=1)
niceqqplot(MU.qp,p=0.95,addcont=TRUE,CI=TRUE)

dev.new()
par(las=1)
niceqqplot(MU.qp,p=0.99,addcont=TRUE,CI=TRUE)

dev.new()
mapofstations(MU.qp)

# Test Wilson & Toumi (2005)
print("Monte-Carlo test for stretched exponential tail assumption:")
a1 <- rnorm(100000,mean=1)
a2 <- rnorm(100000,mean=1)
a3 <- rnorm(100000,mean=1)
a123 <- a1*a2*a3
mu <- mean(a123)
prs <- seq(0.05,1,by=0.05)

dev.new()
par(las=1)
plot(-log(1-prs)*mu,quantile(a123,prs),pch=19,
     main="Product of 3 normal variables",
     xlab="-log(1-p)*mu",ylab="L-moment estimator",
     sub="Monte-Carlo simulation [N=100,000]")
lines(c(0,3),c(0,3),col="grey")

par(las=1)
compare1q(mu.qp,mu.qp.rcm)

}





