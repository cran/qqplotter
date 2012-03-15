# Test if an extrapolation of the PCA from qqplotter works for
# predicting higher quantiles: 99.9% and higher. Compare with
# return levels based on GEV and GPD.

# Note - this script only runs if the GDCN and ESCN data are located in
# the specified paths, or the extracts are present as rda-files
# (which may be available at ftp://ftp.met.no/users/rasmusb/).
# Or contact: rasmus.benestad@met.no

#pdiv <- 0.99
#cpx <- "max(probs)"
#cpx=pdiv

extract.qdata <- function(do.recompute=FALSE,
                          disk2.path="GDCN-disk2/",plot=FALSE) {
  do.gdcn1 <- !file.exists("mu.qp.extr1.rda")
  do.gdcn2 <- !file.exists("mu.qp.extr2.rda")
  do.escn3 <- !file.exists("mu.qp.extr3.rda")

  if (do.recompute) {
    prs=c( seq(0.50,0.95,by=0.05),c(0.96,0.97,0.98,0.99),
           seq(0.9905,0.9995,by=0.0005) )
    x.0=1
    Pmax=300
    max.stations=NULL

    if (do.gdcn1) {
      qplotGDCN(x.0,Pmax,prs,max.stations,xlim=c(0,400),
                ylim=c(0,400)) -> mu.qp.extr1
      save(file="mu.qp.extr1.rda",mu.qp.extr1)
    } else load("mu.qp.extr1.rda")

    if (do.gdcn2) {
      qplotGDCN(x.0,Pmax,prs,max.stations,xlim=c(0,400),ylim=c(0,400),
                path=disk2.path) -> mu.qp.extr2
      save(file="mu.qp.extr2.rda",mu.qp.extr2)
    } else load("mu.qp.extr2.rda")
      mu.qp.extr2a <- qweed(mu.qp.extr2,
         crit="attr(mu.qp,'latitude')> 50")
      mu.qp.extr2b <- qweed(mu.qp.extr2,
       crit="abs(attr(mu.qp,'longitude'))> 50")
    mu.qp.extr2 <- qcat(mu.qp.extr2a,mu.qp.extr2b)
    cleanduplicates(mu.qp.extr2,silent=TRUE,plot=plot)
    if (plot) dev.new()
  
    if (do.escn3) {
      qplotESCN(x.0,Pmax,prs,max.stations,xlim=c(0,400),
                ylim=c(0,400)) -> mu.qp.extr3
      save(file="mu.qp.extr3.rda",mu.qp.extr3)
    } else load("mu.qp.extr3.rda")

    # 3-rd largest for q_0.9995 and 6th largest for q_0.999
    mu.qp.extr1 <- qweed(mu.qp.extr1,crit="attr(mu.qp,'n.wet')> 6000")
    mu.qp.extr2 <- qweed(mu.qp.extr2,crit="attr(mu.qp,'n.wet')> 6000")
    mu.qp.extr3 <- qweed(mu.qp.extr3,crit="attr(mu.qp,'n.wet')> 6000")
    mu.qp.extr12 <- qcat(mu.qp.extr1,mu.qp.extr2)
    mu.qp.all <- qcat(mu.qp.extr12,mu.qp.extr3)
    #mapofstations(mu.qp.all)
    save(file="mu.qp.all.rda",mu.qp.all)
  } else load("mu.qp.all.rda")
  mu.qp.all <- cleanduplicates(mu.qp.all,silent=TRUE,plot=plot)
  invisible(mu.qp.all)
}


qextrap.test <- function(mu.qp.all,pdiv=0.99) {
# Cut off the high-percentiles:

  mu.qp.lo <- qsubset(mu.qp.all,
                      crit=paste("attr(mu.qp,'probabilities') <=",pdiv)) 
  mu.qp.up <- qsubset(mu.qp.all,
                      crit=paste("attr(mu.qp,'probabilities') >",pdiv)) 

# Carry out the PCA: X = U %*% W %*% t(V)

  pca.all <- qPCA(mu.qp.all,plot=FALSE)
  pca.lo <- qPCA(mu.qp.lo,plot=FALSE)
  pca.up <- qPCA(mu.qp.up,plot=FALSE)

# Compare the PCA with the extended PCA:
# The sign of the EOF may be aritrary: show absolute values
# |U| = I, but the dimensions of U differ for different number of
# percentiles -> scale the EOFs to compensate for this.
  W1.lo <- pca.lo$d[1]*sd(pca.lo$v[,1])
  W1.up <- pca.up$d[1]*sd(pca.up$v[,1])
  W1.all <- pca.all$d[1]*sd(pca.all$v[,1])
  W2.all <- pca.all$d[2]*sd(pca.all$v[,2])

# Calibrate model: for 1st EOF 
  calibrate1 <- data.frame(y=abs(W1.lo*pca.lo$x1[,1]),
                           x=abs(W1.lo*pca.lo$qp[,1]))
  extension1 <- data.frame(x=c(abs(W1.lo*pca.lo$qp[,1]),
                               abs(W1.up*pca.up$qp[,1])))
  pcafit1 <- lm(y ~ x + I(x^2),data=calibrate1)

  dev.new()
  par(bty="n",col.sub="grey")                         
                         
  plot(c(0,80),c(0,80),type="l",col="grey70",lty=2,
       main = "Leading PCA mode", 
       xlab = "-log(1-p)*mu [mm/day]", ylab = "quantile(X,p) [mm/day]", 
       sub = "[qextrapolation.R]")
  grid()
  points(abs(W1.lo*pca.lo$qp[,1]),abs(W1.lo*pca.lo$x1[,1]),pch=19,type="b")
  points(abs(W1.up*pca.up$qp[,1]),abs(W1.up*pca.up$x1[,1]),pch=19,
         col="grey",type="b")
  lines(abs(W1.lo*pca.lo$qp[,1]),abs(predict(pcafit1)),lty=1,lwd=2,
        col="darkred")
  lines(extension1$x,abs(predict(pcafit1,newdata=extension1)),
        lty=2,col="red")
  legend(0,80,c("original","extended","fit"),
         col=c("black","grey","red"),pch=c(19,21,NA),lty=c(0,0,2),bty="n")

  # Second PCA mode:
  dev.new()
  par(bty="n",col.sub="grey")
  
  W2.lo <- pca.lo$d[2]*sd(pca.lo$v[,2])
  W2.up <- pca.up$d[2]*sd(pca.up$v[,2])
                         
  calibrate2 <- data.frame(y=abs(W2.lo*pca.lo$x1[,2]),
                           x=abs(W2.lo*pca.lo$qp[,2]) )
  extension2 <- data.frame(x=c(abs(W2.lo*pca.lo$qp[,2]),
                               abs(W2.up*pca.up$qp[,2])))
  pcafit2 <- lm(y ~ x,data=calibrate2)

  plot(c(0,30),c(0,30),type="l",col="grey70",lty=2,
       main = "Second PCA mode", 
       xlab = "-log(1-p)*mu [mm/day]", ylab = "quantile(X,p) [mm/day]")
  grid()
  points(abs(W2.lo*pca.lo$qp[,2]),abs(W2.lo*pca.lo$x1[,2]),pch=19,type="b")
  points(abs(W2.up*pca.up$qp[,2]),abs(W2.up*pca.up$x1[,2]),pch=19,
         col="grey",type="b")
  lines(abs(W2.lo*pca.lo$qp[,2]),abs(predict(pcafit2)),lty=1,lwd=2,
        col="darkred")
  lines(extension2$x,abs(predict(pcafit2,newdata=extension2)),
     lty=2,col="red")
  legend(0,30,c("original","extended","fit"),
         col=c("black","grey","red"),pch=c(19,21,NA),lty=c(0,0,2),bty="n")

}



predqPCA <- function(pca,pca.all) {
    # Calibrate model: for 1st EOF

  calibrate1 <- data.frame(y=pca$x1[,1],
                           x=pca$qp[,1])
  extension1 <- data.frame(x=pca.all$qp[,1])
  pcafit1 <- lm(y ~ x + I(x^2),data=calibrate1)

  # Second PCA mode:                         
  calibrate2 <- data.frame(y=pca$x1[,2],
                           x=pca$qp[,2] )
  extension2 <- data.frame(x=pca.all$qp[,2])
  pcafit2 <- lm(y ~ x,data=calibrate2)

  pca.syn <- pca
  attr(pca.syn,"n1") <- attr(pca.all,"n1")
  attr(pca.syn,'probabilities') <- attr(pca.all,'probabilities')
  print(length(extension1$x))
  pca.syn$u <- cbind( c(predict(pcafit1,newdata=extension1),extension1$x),
                      c(predict(pcafit2,newdata=extension2),extension2$x) )
  print(dim(pca.syn$u))
  pca.syn$x1 <- cbind(predict(pcafit1,newdata=extension1),
                      predict(pcafit2,newdata=extension2))
  pca.syn$qp <- cbind(extension1$x,extension2$x)
  I <- det(t(pca.syn$u) %*% pca.syn$u)
  print(paste("|UU^T|=",I,"length(pca$u[,1])=",length(pca$u[,1]),
              "length(pca.all$u[,1])=",length(pca.all$u[,1]),
              "length(pca.all$u[,1])/length(pca$u[,1])=",
              round(length(pca.all$u[,1])/length(pca$u[,1]),4)))
  print(paste("pca.all$d[1]=",pca.all$d[1],"pca$d[1]=",pca$d[1],
              "pca.all$d[1]/pca$d[1]=",round(pca.all$d[1]/pca$d[1],4)))
  

  # Some kind of scaling is needed as d is different for PCAs with different
  # dimensions (SVD -> |UU^T|=1).
  
#  pca.syn$d <- pca.all$d[1:2]
  pca.syn$d <- pca.syn$d[1:2]*I^2
  pca.syn$v <- pca.syn$v[,1:2]
  pca.syn$var <- pca.all$var[1:2]
  results <- list(pca.syn=pca.syn,pcafit1=pcafit1,pcafit2=pcafit2,
                  extension1=extension1,extension2=extension2)
  invisible(results)
}


synthPCA.qp <- function(mu.qp.all,pdiv=0.99,plot=TRUE) {
  # Extrapolate the PCA shaps given a set of quantiles, which are 
  # estimated according to qp = -ln(1-p) mu.
  # Hence, mu, p, and two PC-loadings are the only parameters requied
  # to estimate the actual percentiles for any station.
  
  # Cut off the high-percentiles for testing purposes:

  mu.qp.lo <- qsubset(mu.qp.all,
                      crit=paste("attr(mu.qp,'probabilities') <=",pdiv)) 
  mu.qp.up <- qsubset(mu.qp.all,
                      crit=paste("attr(mu.qp,'probabilities') >",pdiv)) 

# Carry out the PCA for the different samples: X = U %*% W %*% t(V)

  pca.all <- qPCA(mu.qp.all,plot=FALSE)
  pca.lo <- qPCA(mu.qp.lo,plot=FALSE)
  pca.up <- qPCA(mu.qp.up,plot=FALSE)


  qmodels <- predqPCA(pca.lo,pca.all)
  attach(qmodels)

  # In-sample reference: Calibrate model: for 1st EOF
  # Reference level: using dependent data:

  Qmodels <- predqPCA(pca.all,pca.all)
  Pcafit1 <- Qmodels$pcafit1
  Pcafit2 <- Qmodels$pcafit2
  pca.ref <- Qmodels$pca.syn

  i1 <- 1:length(extension1$x)
  i2 <- i1 + length(extension1$x)

  W1.lo <- pca.lo$d[1]*sd(pca.lo$v[,1])
  W2.lo <- pca.lo$d[2]*sd(pca.lo$v[,2])
  W1.all <- pca.all$d[1]*sd(pca.all$v[,1])
  W2.all <- pca.all$d[2]*sd(pca.all$v[,2])

  # Compare the PCA with the extended PCA:

  if (plot) {
    dev.new()
    par(bty="n",col.sub="grey")
    plot(W1.all*pca.all$u[i2,1],W1.all*pca.all$u[i1,1],type="l",lwd=4,
       main="Test leading PCA shapes from qp:",col="grey",
       xlim=c(-60,0),ylim=c(-60,2),
       xlab = "-log(1-p)*mu [mm/day]", ylab = "quantile(X,p) [mm/day]",
       sub = "[qextrapolation.R]")
    points(W1.all*pca.all$qp[,1],W1.all*pca.all$x1[,1],type="b",col="grey")
    points(W1.lo*pca.lo$qp[,1],W1.lo*pca.lo$x1[,1],type="b")
    lines(W1.all*pca.ref$u[i2,1],W1.all*pca.ref$u[i1,1],col="pink",lwd=3)
    lines(W1.lo*pca.syn$u[i2,1],W1.lo*pca.syn$u[i1,1],lty=2,col="red",lwd=3)
    legend(-60,0,c("original all","dependent all","out-of-sample",
                 "low-part predict","all predict"),
         col=c("grey","pink","red","black","grey"),
         pch=c(NA,NA,NA,21,21),lwd=c(3,3,3,1,1),lty=c(1,1,2,1,1),bty="n")

    dev.new()
    par(bty="n",col.sub="grey")
    plot(W2.all*pca.all$u[i2,2],W2.all*pca.all$u[i1,2],type="l",lwd=4,
       main="Test second PCA shapes from qp:",col="grey",
       xlim=c(-15,5),ylim=c(-30,20),
       xlab = "-log(1-p)*mu [mm/day]", ylab = "quantile(X,p) [mm/day]",
       sub = "[qextrapolation.R]")
    points(W2.all*pca.all$qp[,2],W2.all*pca.all$x1[,2],type="b",col="grey")
    points(W2.lo*pca.lo$qp[,2],W2.lo*pca.lo$x1[,2],type="b")
    lines(W2.all*pca.ref$u[i2,2],W2.all*pca.ref$u[i1,2],col="pink",lwd=3)
    lines(W2.lo*pca.syn$u[i2,2],W2.lo*pca.syn$u[i1,2],lty=2,col="red",lwd=3)
    legend(-15,-20,c("original","dependent","out-of-sample",
                 "low-part predict","all predict"),
         col=c("grey","pink","red","black","grey"),
         pch=c(NA,NA,NA,21,21),lwd=c(3,3,3,1,1),lty=c(1,1,2,1,1),bty="n")
  }
  
  results <- list(pca.all=pca.all,pca.lo=pca.lo,pca.syn=pca.syn,
                  pca.ref=pca.ref,pdiv=pdiv,
                  pcafit1=pcafit1,pcafit2=pcafit2)
  detach(qmodels)
  invisible(results)

}




qpextrapPCA <- function(pcaextrap,cpx="max(probs)",neofs=2,plot=TRUE) {
  attach(pcaextrap)
  ns <- length(pca.syn$v[,1])
  probs <- attr(mu.qp.all,"probabilities")
  np <- length(probs)

  qp <- rep(NA,ns*np);  dim(qp) <- c(np,ns); qp.ref <- qp; qp.all <- qp
  for (ip in 1:np) {
    qPCA2quantile(pca.syn$v[,1:neofs],pca.syn,p=probs[ip],
                  silent=TRUE) -> qp[ip,]
    qPCA2quantile(pca.all$v[,1:neofs],pca.all,p=probs[ip],
                  silent=TRUE) -> qp.all[ip,]
    qPCA2quantile(pca.ref$v[,1:neofs],pca.ref,p=probs[ip],
                  silent=TRUE) -> qp.ref[ip,]
  }

#print(dim(mu.qp.all)); print(dim(qp))

  ilo <- probs <= pdiv
  eval(parse(text=paste("px <-",cpx)))

  if (plot) {
    dev.new()
    par(bty="n",col.sub="grey")
    plot(c(0,500),c(0,500),type="l",
       main = paste("All percentiles for prob below",pdiv), 
       xlab = "quantile(X,p) [mm/day]", ylab ="f(PCA) [mm/day]", 
       sub = "[qextrapolation.R]")
    points(mu.qp.all[ilo,],qp.ref[ilo,],pch=19,col="pink",cex=0.7)
    points(mu.qp.all[ilo,],qp.all[ilo,],pch=19,col="grey",cex=0.7)
    points(mu.qp.all[ilo,],qp[ilo,],pch=19,col="red",cex=0.7)
    legend(0,500,c("no extrapolation","in-sample","out-of-sample"),
         col=c("grey","pink","red"),pch=19,bty="n")

    dev.new()
    par(bty="n",col.sub="grey")
    plot(c(0,500),c(0,500),type="l",
       main = paste("All percentiles for prob above",pdiv), 
       xlab = "quantile(X,p) [mm/day]", ylab ="f(PCA) [mm/day]", 
       sub = "[qextrapolation.R]")
    points(mu.qp.all[!ilo,],qp.ref[!ilo,],pch=19,col="pink",cex=0.7)
    points(mu.qp.all[!ilo,],qp.all[!ilo,],pch=19,col="grey",cex=0.7)
    points(mu.qp.all[!ilo,],qp[!ilo,],pch=19,col="red",cex=0.7)
    legend(0,500,c("no extrapolation","in-sample","out-of-sample"),
         col=c("grey","pink","red"),pch=19,bty="n")

    dev.new()
    par(bty="n",col.sub="grey")
    plot(c(0,200),c(0,200),type="l",
       main = paste(px," percentiles"), 
       xlab = "quantile(X,p) [mm/day]", ylab = "f(PCA) [mm/day]", 
       sub = "[qextrapolation.R]")
    ipx <- is.element(probs,px)
    points(mu.qp.all[ipx,],qp.ref[ipx,],pch=19,col="pink")
    points(mu.qp.all[ipx,],qp.all[ipx,],pch=19,col="grey")
    points(mu.qp.all[ipx,],qp[ipx,],pch=19,col="red")
    legend(60,20,c("no extrapolation","in-sample","out-of-sample"),
         col=c("grey","pink","red"),pch=19,bty="n")
  }
  
  detach(pcaextrap)
  results <- list(mu.qp.all=mu.qp.all,qp.ref=qp.ref,qp.all=qp.all,qp=qp,
                  px=px)
  invisible(results)
}


qextrap.eval <- function(results,plot=TRUE) {
  attach(results)
  probs <- attr(mu.qp.all,"probabilities")
  np <- length(probs)
  rmse <- rep(NA,np)
  meanprop <- rmse
  ci90prop <- rmse
  RP <- mu.qp.all

  f.wet <- attr(mu.qp.all,'n.wet')/
     (attr(mu.qp.all,'n.wet')+attr(mu.qp.all,'n.dry'))
  for (ip in 1:np) {
    rmse[ip] <- sqrt( sum( (qp[ip,] - mu.qp.all[ip,])^2 ) )/length(qp[ip,])
    meanprop[ip] <- mean( 100*(abs(qp[ip,] - mu.qp.all[ip,]))/mu.qp.all[ip,] )
    ci90prop[ip] <- quantile(
                     100*(abs(qp[ip,] - mu.qp.all[ip,]))/mu.qp.all[ip,],0.90)
    RP[ip,] <- 1/(365.25 * (1-probs[ip])*f.wet)
  }

  if (plot) {
    dev.new()
    par(bty="n",col.sub="grey")
    hist(c(RP[is.element(probs,px),]),
       main="Return period (wet+dry days)",
       xlab="years",ylab="cases",
       sub=paste(px," percentiles"))
    grid()

    dev.new()
    par(bty="n",col.sub="grey",col.axis="white")
    plot(exp(range(100*probs)),c(0,100),type="n",
     main="Proportional error estimates of qextrapolation",
     xlab="probability level (%)",ylab="%", 
     sub = "[qextrapolation.R]")
    grid()
    par(col.axis="black")
    axis(1,at=exp((100*probs)),labels=100*probs)
    axis(2)
    lines(exp((100*probs)),meanprop,lwd=3)
    lines(exp((100*probs)),ci90prop,lwd=2,lty=2)

    dev.new()
    par(bty="n",col.sub="grey",col.axis="white")
    plot(exp(100*probs),rmse,type="l",lwd=3,
     main="RMSE of qextrapolation",
     xlab="probability level (%)",ylab="mm/day", 
     sub = "[qextrapolation.R]")
    par(col.axis="black")
    axis(1,at=exp((100*probs)),labels=100*probs)
    axis(2)
    grid()
  }

  qerrors <- list(probs=probs,meanprop=meanprop,ci90prop=ci90prop,rmse=rmse)
  invisible(qerrors)
}

annual.max <-  function(x,n.min=300) {
  years <- as.numeric(rownames(table(attr(x,'year'))))
  annual.max <- rep(NA,length(years)); ndays <- annual.max
  for (i in 1:length(years)) {
    ii <- is.element(attr(x,'year'),years[i]) &
          is.finite(x)
    if (sum(ii)> 300) annual.max[i] <- max(x[ii])
    ndays[i] <- sum(ii)
  }
  attr(annual.max,'year') <- years
  attr(annual.max,'Stnr') <- attr(x,'Stnr')
  attr(annual.max,'location') <- attr(x,'location')
  attr(annual.max,'lon') <- attr(x,'lon')
  attr(annual.max,'lat') <- attr(x,'lat')
  attr(annual.max,'alt') <- attr(x,'alt')
  attr(annual.max,'km.from.coast') <- attr(x,'km.from.coast')
  attr(annual.max,'n.days') <- ndays
  attr(annual.max,'n.min') <- n.min
  invisible(annual.max)
}


fqextr <- function(x=NULL,mu.qp=NULL,x.0=1,neofs=2,plot=TRUE) {
  # Add the station data statistics to the data mu.qo. Perform a qPCA
  # and then extraploate the EOF shape to higher percentiles.

  data(qerrors,envir=environment())
  data(pca.qextr,envir=environment())
  if (is.null(x)) {
    data(samnanger,envir=environment())
    x <- samnanger
  } 
  if (is.null(mu.qp)) data(mu.qp,envir=environment())
  
  dry <- (x < x.0) & is.finite(x)
  wet <- (x >= x.0) & is.finite(x)
  X <- x[wet]
  mu <- mean(X); ave <- mean(x,na.rm=TRUE)
  sigma <- round(var(c(X[X >= x.0]),na.rm=TRUE),2)
  x1 <- quantile(X,attr(mu.qp,'probabilities'))
  quants <- fitsimple(c(x),x.0,attr(mu.qp,"probabilities"))
  rmse <- sqrt(sum((quants$x2 - quants$x1)^2))/length(quants$x2)
  MU.qp <- cbind(x1,mu.qp)
  attr(MU.qp,"description") <- attr(mu.qp,"description")
  attr(MU.qp,"mean_precip") <- c(ave,attr(mu.qp,"mean_precip"))
  attr(MU.qp,"mu") <- c(mu,attr(mu.qp,"mu"))
  attr(MU.qp,"longitude") <- c(attr(x,"lon"),attr(mu.qp,"longitude"))
  attr(MU.qp,"latitude") <- c(attr(x,"lat"),attr(mu.qp,"latitude"))
  attr(MU.qp,"dist2coast.km") <- c(attr(x,"km.from.coast"),
                                   attr(mu.qp,"dist2coast.km"))
  attr(MU.qp,"station_number") <- c(attr(x,"Stnr"),attr(mu.qp,"station_number"))
  attr(MU.qp,"altitude") <- c(attr(x,"altitude"),attr(mu.qp,"altitude"))
  attr(MU.qp,"n.wet") <- c(sum(wet),attr(mu.qp,"n.wet"))
  attr(MU.qp,"n.dry") <- c(sum(dry),attr(mu.qp,"n.dry"))
  attr(MU.qp,"sigma") <- c(sigma,attr(mu.qp,"sigma"))
  attr(MU.qp,"rmse") <- c(rmse,attr(mu.qp,"rmse"))
  attr(MU.qp,"country") <- c(attr(x,"country"),attr(MU.qp,"country"))
  attr(MU.qp,"history") <- paste(attr(mu.qp,"history"),"added one station")
  attr(MU.qp,"probabilities") <- attr(mu.qp, "probabilities")
  attr(MU.qp,"formulae") <- attr(mu.qp,"formulae")
  attr(MU.qp,"class") <- attr(mu.qp,"class")
                          
  pca <- qPCA(MU.qp,plot=FALSE)
  synthesise <- predqPCA(pca,pca.qextr)
  attach(synthesise)
  probs <- attr(pca.qextr,'probabilities'); np <- length(probs)
  qp <- rep(NA,np)
  for (ip in 1:np) qp[ip] <- 
      qPCA2quantile(pca.syn$v[1,1:neofs],pca.syn,p=probs[ip],
                  silent=TRUE)
  x1 <- quantile(X,probs)
  qerr <- approx(qerrors$probs,qerrors$ci90prop/100,probs)$y
  
  if (plot) {
    dev.new()
    par(bty="n",col.sub="grey")
    plot(c(0,200),c(0,200),type="l",
     main=paste(attr(x,"location"),"24-hr wet-day precipitation"),
     xlab="observed quantile (mm/day)",ylab="predicted quantile (mm/day)", 
     sub = "[fqextr]")
    grid()
    points(x1,qp,pch=19,type="b",col="red")
    lines(x1,qp+qp*qerr,lty=2,col="red")
    lines(x1,qp-qp*qerr,lty=2,col="red")
  }

  attr(qp,'location') <- attr(x,'location')
  attr(qp,'Stnr') <- attr(x,'Stnr')
  attr(qp,'lon') <- attr(x,'lon')
  attr(qp,'lat') <- attr(x,'lat')
  attr(qp,'probabilities') <- probs
  attr(qp,'mu') <- mu
  attr(qp,'f.wet') <- sum(wet)/(sum(wet)+sum(dry))
  attr(qp,'n.data') <- sum(wet)+sum(dry)
  attr(qp,'years') <- range(attr(x,'year'))  
  attr(qp,"km.from.coast") <- attr(x,"km.from.coast")
  attr(qp,"altitude") <- attr(x,"altitude")
  
  detach(synthesise)
  invisible(qp)  
}




mu2q <- function(x=NULL,plot=TRUE,x.0=1) {

  # p.d.f. describing a fitted exponential distribution
  f <- function(x,mu) f <- 1/mu* exp(-1*x/mu)
  
  # Estimates the variance from the second moment
  secondmoment <- function(mu,X=seq(0,1000,by=0.1)) {
    dX <- X[2]-X[1]
    s <- sum(X^2*f(X,mu)*dX)
    s
  }
  
  data(qerrors,envir=environment())
  if (is.null(x)) {
    data(samnanger,envir=environment())
    x <- samnanger
  }
  n <- length(x)
  mu <- rep(NA,n); sigma <- mu
  f.wet <- mu
  print(paste("n=",n))
  for (i in 1:n) {
    if (mod(i,1000)==0) plot(1:n,mu,type="l",main="wet-day mean")   
    y <- x[1:i]
    wet <- (y >= x.0) & is.finite(y)
    dry <- (y < x.0) & is.finite(y)
    Y <- y[wet]
    mu[i] <- sum(Y)/sum(is.finite(y))
    f.wet[i] <- sum(wet)/(sum(wet) + sum(dry))
    sigma[i] <- sqrt( secondmoment(mu[i]) )/sqrt(i)
  }

  par(bty="n",col.sub="grey")
  plot(1:n,mu,type="l",lwd=2,
       main=paste(attr(x,'location'),"wet-day mean"),
       ylab="mean wet-day precipitation (mm/day)",xlab="record length")
  grid()
  lines(1:n,mu + 2*sigma,lty=2,col="grey")
  lines(1:n,mu - 2*sigma,lty=2,col="grey")

  dev.new()
  plot(1:n,f.wet,type="l",lwd=2,
       main=paste(attr(x,'location'),"wet-day fraction"),
       ylab="wet-day frequency (fraction)",xlab="record length")
  grid()
}

# Reconstruct the percentiles from wet-day mean and the qPCA:

mu2extr <- function(mu,neofs=2,pca=NULL,plot=TRUE) {
  if (is.null(pca)) {
    data(pca.qextr,envir=environment()); pca <- pca.qextr
  }
  p <- attr(pca,'probabilities')
  qp <- -log(1-p)*mu
  x1 <- rep(0,length(p))
  mudiff <- abs(mu - attr(pca,'mu'))
  imatch <- (1:length(mudiff))[is.element(mudiff,min(mudiff,na.rm=TRUE))][1]
  pc <- abs(pca$v[imatch,1:neofs])
  print(pc)
  
  # find the right scaling between the percentiles and the values in the EOF:
  # Use the product between scaling and EOF to predict real percentiles:
  for (i in 1:neofs) {
    W <- pca$var[i]/100 
    scale <- data.frame(x=abs(pca$qp[,i]),y=qp)
    scalefit <- lm(y ~ x,data=scale)
    print(summary(scalefit))
    x1 <- x1 + abs(pca$x1[,i])*scalefit$coefficients[2]
  }

  if (plot) {
    par(bty="n",col.sub="grey")
    plot(qp,x1,pch=19,type="b")
    lines(c(0,200),c(0,200),col="grey")
  }

  attr(x1,'qp') <- qp
  attr(x1,'probabilities') <- p
  attr(x1,'mu') <- mu
  attr(x1,'neofs') <- neofs
  attr(x1,'pca') <- pca
  invisible(x1)
}

clim2qextr <- function(x,p=0.95,x.0=1,modelq=NULL) {
  print("clim2qextr: estimating quantiles through modelling PCs.")
  if (is.null(modelq) & file.exists("modelqPCs.rda")) {
    load("modelqPCs.rda"); modelq <- result
  } else modelq <- modelqPCs()

  dry <- (x < x.0) & is.finite(x)
  wet <- (x >= x.0) & is.finite(x)
  f.wet <- sum(wet)/(sum(wet)+sum(dry))
  if (!is.null(attr(x,'km.from.coast')) & !is.null(attr(x,'alt'))) {
    d2c <- attr(x,'km.from.coast')
    altitude <- attr(x,'alt')
  } else {
      data(mu.qp.world)
      data(addland1,envir = environment())
      i <- is.element(attr(mu.qp.world,"station_number"),
                      attr(x,"Station_number"))
      if (sum(i)==1) {
        d2c<- round(min(distAB(attr(mu.qp.world,"longitude")[i],
                                    attr(mu.qp.world,"latitude")[i],
                                    lon.cont,lat.cont),na.rm=TRUE)/1000)
        altitude <-attr(mu.qp.world,"altitude")[i]
      } else  stop(paste("Did not find the station metadata for",
                         attr(x,"Station_number")))
    }
    
    
  predictor <- data.frame(mu=mean(x[wet]),
                          wet.freq=f.wet,
                          dist2coast=d2c,
                          altitude=altitude)
  pc1 <- predict(modelq$pc1model,newdata=predictor)
  pc2 <- predict(modelq$pc2model,newdata=predictor)
  np <- length(p)
  qp <- rep(NA,np)
  for (ip in 1:np)
    qp[ip] <- qPCA2quantile(cbind(pc1,pc2),result$pca,p=p[ip],silent=TRUE)
  attr(qp,'probabilties') <- p
  attr(qp,'f.wet') <- f.wet
  attr(qp,'mu') <- predictor$mu
  attr(qp,'dist2coast') <- predictor$dist2coast
  attr(qp,'altitude') <- predictor$altitude
  invisible(qp)
}


qrl <- function(x=NULL,x.0=1) {

  require(evd)
  data(mu.qp)
  
  if (is.null(x)) {
    data(samnanger,envir=environment())
    x <- samnanger
  } 
  data(qerrors,envir=environment())
  qp <- fqextr(x,plot=FALSE)
  p <- attr(qp,'probabilities')

  X <- x[(x >= x.0) & is.finite(x)]
  mu <- mean(X)
  qp.mu <- mu2extr(mu,plot=FALSE)
  p.clim <- attr(mu.qp,'probabilities')
  qp.clim <- clim2qextr(x,p=p.clim)
  
  #print(qpercs)
  
  RP <- 1/(365.25 * (1-p)*attr(qp,'f.wet'))
  RP.clim <- 1/(365.25 * (1-p.clim)*attr(qp.clim,'f.wet'))
  qerr <- approx(qerrors$probs,qerrors$ci90prop/100,p)$y

  # GEV fit:
  Xyrmax <- annual.max(x)
  extremes <- fgev(Xyrmax)
  par(bty="n",col.sub="grey")
  RLgev <- rl(extremes,main="Return Level: GEV & qPCA",
              xlim=c(1/12,200),ylim=c(10,200),
              sub=attr(x,'location'))
  grid()
  
  points(RP,qp,type="b",pch=19,col="red") 
  lines(RP,qp+qp*qerr,lty=2,col="red")
  lines(RP,qp-qp*qerr,lty=2,col="red")

  points(RP,quantile(x[(x > x.0) & is.finite(x)],p),col="grey")
  lines(rep(1/(365.25 * (1-0.99)*attr(qp,'f.wet')),2),c(0,200),
        lty=3,col="red")

  points(RP,qp.mu,type="b",pch=19,col="blue") 
  points(RP.clim,qp.clim,type="b",pch=19,col="darkgreen") 
  
  legend(20,quantile(qp,0.35),c("qPCA","GEV","mu","local","quantile"),
         col=c("red","black","blue","darkgreen","grey"),
         lty=c(1,1,1,1,0),pch=c(19,4,19,19,21),bty="n")

  arrows(1/12,15,1/(365.25 * (1-0.99)*attr(qp,'f.wet')),15,
         col="red")
  text(1/15,20,"calibration",pos=4,col="red")
}

