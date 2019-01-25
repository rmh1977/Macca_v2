######################################################
# general purpose fitting plots for Macca ############
######################################################
# R Hillary CSIRO 2019 ###############################
######################################################

library(ggplot2)

# plot length data
# args: report object, fishery, years, length bins

plot.lenfreq <- function(rep,ff,yrs,ll) {

  # combine across sexes

  fr <- rf[ff]+1
  sx <- rep$sexrat[yrs,,,fr]
  px <- rep$ply[yrs,,,fr]
  phat <- apply(sx*px,c(1,2),sum)
  pobs <- t(apply(LF[yrs,,ff],1,function(x){x <- x/sum(x)}))

  # translate to dataframe

  df <- expand.grid(obs=NA,pred=NA,y=yrs,l=ll)
  for(y in 1:length(yrs)) {

    df[df$y == yrs[y],'obs'] <- pobs[y,]
    df[df$y == yrs[y],'pred'] <- phat[y,]
 
  }

  pp <- ggplot(df,aes(x=l,y=obs))+geom_point(colour='magenta',size=1)+facet_wrap(~y,ncol=7)+geom_line(aes(x=l,y=pred),colour='blue')+xlab("length")+ylab("frequency")
  pp

}

# length-conditional age data
# report object, fishery, sex, years, lengths, ages

plot.agelendat <- function(rep,ff,sex,yrs,ll,ages) {

  # plot mean age (obs and pred) for each length class

  fr <- rf[ff]+1
  if(sex == 'male') {
    
    dat <- ALFm[yrs,,,ff]
    px <- rep$paaly[yrs,,,2,ff] 

  }

  if(sex == 'female') {
    
    dat <- ALFf[yrs,,,ff] 
    px <- rep$paaly[yrs,,,1,ff]  

  }
  
  df <- expand.grid(obs=NA,mu=NA,lq=NA,uq=NA,y=yrs,l=ll)

  for(y in 1:length(yrs)) {
    for(l in 1:length(ll)) {
  
      # observed mean age (if there are data in that length bin)     

      neff <- sum(dat[y,l,])
      if(neff > 0) {

        pobs <- dat[y,l,]/neff
        muaobs <- sum(pobs*ages)
        df[df$y == yrs[y] & df$l == ll[l],'obs'] <- muaobs 
      }
        
      muahat <- sum(px[y,,l]*ages)
      sdahat <- sqrt(sum(px[y,,l]*(ages-muahat)^2))
        
      df[df$y == yrs[y] & df$l == ll[l],'mu'] <- muahat
      df[df$y == yrs[y] & df$l == ll[l],'lq'] <- ifelse(muahat-1.96*sdahat > 0,muahat-1.96*sdahat,0)
      df[df$y == yrs[y] & df$l == ll[l],'uq'] <- muahat+1.96*sdahat 

    }
  }

   # plot it

  pp <- ggplot(df,aes(x=l,y=obs))+geom_point(colour='magenta',size=1)+facet_wrap(~y,ncol=7)+geom_line(aes(x=l,y=mu),colour='blue')+geom_line(aes(x=l,y=uq),linetype='dashed',colour='blue')+geom_line(aes(x=l,y=lq),linetype='dashed',colour='blue')+xlab("length")+ylab("mean age")
  pp 

}

##############################
# plot the tagging data fits #
##############################

# by release event (across length-class and recap area)

plot.tagfits.relyr <- function(rep,tmin,sex,nrec) {

  nt <- length(tagrel)
  xdf <- expand.grid(yrel=NA,tev=NA,nrel=NA,lrel=NA,srel=NA,rrel=NA,yrec=rep(NA,nrec),rrec=1:nr,obs=NA,mu=NA,lq=NA,uq=NA)

  for(t in 1:nt) {

    px <- as.vector(rep$prhat[t,2:(nrec+1),])
    yx <- tagcov[t,4]+1
    obs <- as.vector(tagrec[t,1:nrec,])
    nT <- tagrel[t]
    rhat <- as.vector(nT*px)
    sdrhat <- as.vector(sqrt(nT*px*(1-px)))
    yrec <- (tagcov[t,4]+2):(tagcov[t,4]+1+nrec)
    lx <- mulbins[tagcov[t,2]+1]
    if(tagcov[t,1]==0) sx <- 'female'
    if(tagcov[t,1]==1) sx <- 'male'
    rx <- tagcov[t,3]+1
    xdf[1:(nr*nrec),'tev'] <- t
    xdf[1:(nr*nrec),'nrel'] <- nT
    xdf[1:(nr*nrec),'lrel'] <- lx
    xdf[1:(nr*nrec),'srel'] <- sx 
    xdf[1:(nr*nrec),'rrel'] <- rx 
    xdf[1:(nr*nrec),'yrel'] <- yx
    xdf[1:(nr*nrec),'yrec'] <- c((yx+1):(yx+nrec),(yx+1):(yx+nrec))
    xdf[1:(nr*nrec),'obs'] <- obs
    xdf[1:(nr*nrec),'mu'] <- rhat
    lqx <- rhat-1.96*sdrhat
    lqx[lqx < 0] <- 0
    xdf[1:(nr*nrec),'lq'] <- lqx
    xdf[1:(nr*nrec),'uq'] <- rhat+1.96*sdrhat 
    if(t == 1) df <- xdf
    if(t > 1) df <- rbind(df,xdf[1:(nr*nrec),])

  }

  dff <- subset(df,nrel >= tmin)
  pp <- ggplot(dff,aes(x=yrec,y=obs))+geom_point(size=1,colour='magenta')+geom_line(aes(x=yrec,y=mu),colour='blue')+geom_line(aes(x=yrec,y=lq),colour='blue',linetype='dashed')+facet_grid(tev~rrec) 
  pp

}

# by release area

plot.tagfits.relarea <- function(rep,tmin,sex,nrec) {

  nt <- length(tagrel)
  xdf <- expand.grid(yrel=NA,tev=NA,nrel=NA,lrel=NA,srel=NA,rrel=NA,yrec=rep(NA,nrec),rrec=1:nr,obs=NA,mu=NA,lq=NA,uq=NA)

  for(t in 1:nt) {

    px <- as.vector(rep$prhat[t,2:(nrec+1),])
    yx <- tagcov[t,4]+1
    obs <- as.vector(tagrec[t,1:nrec,])
    nT <- tagrel[t]
    rhat <- as.vector(nT*px)
    sdrhat <- as.vector(sqrt(nT*px*(1-px)))
    yrec <- (tagcov[t,4]+2):(tagcov[t,4]+1+nrec)
    lx <- mulbins[tagcov[t,2]+1]
    if(tagcov[t,1]==0) sx <- 'female'
    if(tagcov[t,1]==1) sx <- 'male'
    rx <- tagcov[t,3]+1
    xdf[1:(nr*nrec),'tev'] <- t
    xdf[1:(nr*nrec),'nrel'] <- nT
    xdf[1:(nr*nrec),'lrel'] <- lx
    xdf[1:(nr*nrec),'srel'] <- sx 
    xdf[1:(nr*nrec),'rrel'] <- rx 
    xdf[1:(nr*nrec),'yrel'] <- yx
    xdf[1:(nr*nrec),'yrec'] <- c((yx+1):(yx+nrec),(yx+1):(yx+nrec))
    xdf[1:(nr*nrec),'obs'] <- obs
    xdf[1:(nr*nrec),'mu'] <- rhat
    lqx <- rhat-1.96*sdrhat
    lqx[lqx < 0] <- 0
    xdf[1:(nr*nrec),'lq'] <- lqx
    xdf[1:(nr*nrec),'uq'] <- rhat+1.96*sdrhat 
    if(t == 1) df <- xdf
    if(t > 1) df <- rbind(df,xdf[1:(nr*nrec),])

  }

  dff <- subset(df,nrel >= tmin)
  pp <- ggplot(dff,aes(x=yrec,y=obs))+geom_point(size=1,colour='magenta')+geom_line(aes(x=yrec,y=mu),colour='blue')+geom_line(aes(x=yrec,y=lq),colour='blue',linetype='dashed')+facet_grid(tev~rrec) 
  pp

}

