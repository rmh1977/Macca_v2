######################################################
# general purpose fitting plots for Macca ############
######################################################
# R Hillary CSIRO 2019 ###############################
######################################################

library(ggplot2)

# useful functions

logit <- function(x){return(log(x/(1-x)))}
ilogit <- function(x){return(1/(1+exp(-x)))}

# plot length data
# args: report object, fishery, years, length bins

plot.lenfreq <- function(lfdat,lfcov,rep,ff,phiLF) {

  # combine across sexes

  rx <- rf[ff]+1
  fx <- which(ffnm==ff)
  idx <- which(lfcov[,2] == fx-1)
  yx <- sort(lfcov[idx,1])
  sx <- rep$sexrat[yx+1,,,rx]
  px <- rep$ply[yx+1,,,fx]
  phat <- apply(sx*px,c(1,2),sum)
  phat <- t(apply(phat,1,function(x){x <- x/sum(x)}))
  neff <- lfdat[idx,1]/phiLF
  nhat <- apply(phat,2,function(x,neff){x <- x*neff},neff)
  nobs <- lfdat[idx,-1]/phiLF

  # translate to dataframe

  yz <- yx+min(yrs)
  df <- expand.grid(obs=NA,pred=NA,y=range(yz)[1]:range(yz)[2],l=mulbins)
  for(y in 1:length(yz)) {

    df[df$y == yz[y],'obs'] <- nobs[y,]
    df[df$y == yz[y],'pred'] <- nhat[y,]
 
  }

  pp <- ggplot(df,aes(x=l,y=obs))+geom_point(colour='magenta',size=1)+facet_wrap(~y)+geom_line(aes(x=l,y=pred),colour='blue')+xlab("length")+ylab("frequency")
  pp

}

# length-conditional age data
# report object, fishery, sex, years, lengths, ages

plot.agelendat <- function(alf,alfcov,rep,ff,sex) {

  # plot mean age (obs and pred) for each length class


  fx <- which(ffnm==ff)
  idx <- which(alfcov[,2] == fx-1)
  alfcovx <- alfcov[idx,]
  alfx <- alf[idx,]
  sx <- ifelse(sex=='female',0,1)
  idx2 <- which(alfcovx[,3] == sx)
  alfcovx <- alfcovx[idx2,]
  alfx <- alfx[idx2,]
  yx <- sort(alfcovx[,1])
  yz <- min(yrs)+unique(yx)
  
  df <- expand.grid(obs=NA,mu=NA,lq=NA,uq=NA,y=range(yz)[1]:range(yz)[2],l=mulbins)

  for(i in 1:dim(alfcovx)[1]) {

    # observed mean age (if there are data in that length bin)     

    nobs <- alfx[i,-1]
    pobs <- nobs/sum(nobs)
    yy <- alfcovx[i,1]+1
    ll <- alfcovx[i,4]+1
    muaobs <- sum(ages*pobs)  
    df[df$y == yrs[yy] & df$l == mulbins[ll],'obs'] <- muaobs 
    phat <- rep$paaly[yy,,ll,sx+1,fx]  
    muahat <- sum(phat*ages)
    sdahat <- sqrt(sum(phat*(ages-muahat)^2))
    df[df$y == yrs[yy] & df$l == mulbins[ll],'mu'] <- muahat
    df[df$y == yrs[yy] & df$l == mulbins[ll],'lq'] <- ifelse(muahat-1.96*sdahat > 0,muahat-1.96*sdahat,0)
    df[df$y == yrs[yy] & df$l == mulbins[ll],'uq'] <- muahat+1.96*sdahat 
    
  }

   # plot it

  pp <- ggplot(df,aes(x=l,y=obs))+geom_point(colour='magenta',size=1)+facet_wrap(~y)+geom_line(aes(x=l,y=mu),colour='blue')+geom_line(aes(x=l,y=uq),linetype='dashed',colour='blue')+geom_line(aes(x=l,y=lq),linetype='dashed',colour='blue')+xlab("length")+ylab("mean age")
  pp 

}

##############################
# plot the tagging data fits #
##############################

# by release event (across length-class and recap area)

plot.tagfits.relyr <- function(rep,tagrel,tagnrec,tagcov,tagrec,nrecmax) {

  nt <- length(tagrel)
  yrng <- range(unique(tagcov[,3]+yrs[1]))
  yrel <- yrng[1]:yrng[2]
  xdf <- expand.grid(yrel=yrel,yrec=yrel+1,obs=NA,pred=NA)
  xdf <- subset(xdf,yrec>yrel)
  xdf <- subset(xdf,yrec-yrel<=nrecmax)
  totrec <- apply(tagrec,c(1,2),sum)
  prhat <- rep$prhat[,-1,]
  zz <- apply(prhat,c(2,3),function(x,tagrel){x*tagrel},tagrel)
  totrhat <- apply(zz,c(1,2),sum)

  for(y in yrel) {

    iy <- y-yrs[1]
    idx <- which(tagcov[,3]==iy)
    nT <- sum(tagrel[idx])
    if(nT > 0) {
   
      Robs <- apply(totrec[idx,],2,sum)
      Rhat <- apply(totrhat[idx,],2,sum)
      xdf[xdf$yrel == y,'obs'] <- Robs[1:tagnrec[idx][1]]
      xdf[xdf$yrel == y,'pred'] <- Rhat[1:tagnrec[idx][1]]

    }

  }

  pp <- ggplot(xdf,aes(x=yrec,y=obs))+geom_point(colour='blue')+geom_point(aes(x=yrec,y=pred),colour='magenta',shape=17)+facet_wrap(~yrel)
  pp

}

# total by rel year (across rec year, length-class and recap area)

plot.tagfits.relyr.tot <- function(rep,tagrel,tagnrec,tagcov,tagrec,nrecmax) {

  nt <- length(tagrel)
  yrng <- range(unique(tagcov[,3]+yrs[1]))
  yrel <- yrng[1]:yrng[2]
  xdf <- expand.grid(yrel=yrel,yrec=yrel+1,obs=NA,pred=NA)
  xdf <- subset(xdf,yrec>yrel)
  xdf <- subset(xdf,yrec-yrel<=nrecmax)
  totrec <- apply(tagrec,c(1,2),sum)
  prhat <- rep$prhat[,-1,]
  zz <- apply(prhat,c(2,3),function(x,tagrel){x*tagrel},tagrel)
  totrhat <- apply(zz,c(1,2),sum)

  for(y in yrel) {

    iy <- y-yrs[1]
    idx <- which(tagcov[,3]==iy)
    nT <- sum(tagrel[idx])
    if(nT > 0) {
   
      Robs <- apply(totrec[idx,],2,sum)
      Rhat <- apply(totrhat[idx,],2,sum)
      xdf[xdf$yrel == y,'obs'] <- Robs[1:tagnrec[idx][1]]
      xdf[xdf$yrel == y,'pred'] <- Rhat[1:tagnrec[idx][1]]

    }

  }   

  ydf <- aggregate(xdf$obs,by=list(yrel=xdf$yrel),FUN=sum,na.rm=T) 
  zdf <- aggregate(xdf$pred,by=list(yrel=xdf$yrel),FUN=sum,na.rm=T)
  wdf <- data.frame(yrel=ydf$yrel,obs=ydf$x,pred=zdf$x)

  pp <- ggplot(wdf,aes(x=yrel,y=obs))+geom_point(colour='blue')+geom_point(aes(x=yrel,y=pred),colour='magenta',shape=17)
  pp

}

# total by rec year (across rel year, length-class and recap area)

plot.tagfits.recyr.tot <- function(rep,tagrel,tagnrec,tagcov,tagrec,nrecmax) {

  nt <- length(tagrel)
  yrng <- range(unique(tagcov[,3]+yrs[1]))
  yrel <- yrng[1]:yrng[2]
  xdf <- expand.grid(yrel=yrel,yrec=yrel+1,obs=NA,pred=NA)
  xdf <- subset(xdf,yrec>yrel)
  xdf <- subset(xdf,yrec-yrel<=nrecmax)
  totrec <- apply(tagrec,c(1,2),sum)
  prhat <- rep$prhat[,-1,]
  zz <- apply(prhat,c(2,3),function(x,tagrel){x*tagrel},tagrel)
  totrhat <- apply(zz,c(1,2),sum)

  for(y in yrel) {

    iy <- y-yrs[1]
    idx <- which(tagcov[,3]==iy)
    nT <- sum(tagrel[idx])
    if(nT > 0) {
   
      Robs <- apply(totrec[idx,],2,sum)
      Rhat <- apply(totrhat[idx,],2,sum)
      xdf[xdf$yrel == y,'obs'] <- Robs[1:tagnrec[idx][1]]
      xdf[xdf$yrel == y,'pred'] <- Rhat[1:tagnrec[idx][1]]

    }

  }

  ydf <- aggregate(xdf$obs,by=list(yrec=xdf$yrec),FUN=sum,na.rm=T) 
  zdf <- aggregate(xdf$pred,by=list(yrec=xdf$yrec),FUN=sum,na.rm=T)
  wdf <- data.frame(yrec=ydf$yrec,obs=ydf$x,pred=zdf$x)

  pp <- ggplot(wdf,aes(x=yrec,y=obs))+geom_point(colour='blue')+geom_point(aes(x=yrec,y=pred),colour='magenta',shape=17)
  pp

}

plot.tagfits.recyrarea <- function(rep,tagrel,tagnrec,tagcov,tagrec,nrecmax) {

  nt <- length(tagrel)
  yrng <- range(unique(tagcov[,3]+yrs[1]))
  yrel <- yrng[1]:yrng[2] 
  xdf <- expand.grid(yrel=yrel,yrec=yrel+1,rrec=c('North','South'),obs=NA,pred=NA)
  xdf <- subset(xdf,yrec>yrel)
  xdf <- subset(xdf,yrec-yrel<=nrecmax)
  robs <- tagrec
  rhat <- apply(rep$prhat[,-1,],c(2,3),function(x,tagrel){x <- x*tagrel},tagrel)

  rnm <- c("North","South")
  for(y in yrel) {

    iy <- y-yrs[1]
    idx <- which(tagcov[,3]==iy)
    nT <- sum(tagrel[idx])
    if(nT > 0) {
   
      for(r in 1:nr) {

        Robs <- apply(robs[idx,,r],2,sum)
        Rhat <- apply(rhat[idx,,r],2,sum)
        xdf[xdf$yrel == y & xdf$rrec == rnm[r],'obs'] <- Robs[1:tagnrec[idx][1]]
        xdf[xdf$yrel == y & xdf$rrec == rnm[r],'pred'] <- Rhat[1:tagnrec[idx][1]]

      }
    }

  } 

  ydf <- aggregate(xdf$obs,by=list(yrec=xdf$yrec,rrec=xdf$rrec),FUN=sum,na.rm=T) 
  zdf <- aggregate(xdf$pred,by=list(yrec=xdf$yrec,rrec=xdf$rrec),FUN=sum,na.rm=T)
  wdf <- data.frame(yrec=ydf$yrec,rrec=ydf$rrec,obs=ydf$x,pred=zdf$x)

  pp <- ggplot(wdf,aes(x=yrec,y=obs))+geom_point(colour='blue')+geom_point(aes(x=yrec,y=pred),colour='magenta',shape=17)+facet_wrap(~rrec,scales='free_y')
  pp


}
