######################################################
# testing population dynamics calculations ###########
######################################################
# R Hillary CSIRO 2018 ###############################
######################################################

ages <- 1:52
nages <- length(ages)
R <- 2
ns <- 2
R0 <- 1e+6
M <- rep(0.13,nages)
zeta <- c(0.5,0.5) # sex ratio
eta <- c(0.45,0.55) # spatial recruitment ratio
pim <- array(dim=c(R,R,nages)) # movement matrices
sm50 <- 139
sm95 <- 185
a.wl <- 4.4e-6
b.wl <- 3.14
Linf <- c(165,165)
k <- c(0.057,0.052)
t0 <- c(-0.19,-0.34)
mula <- m <- w <- array(dim=c(ns,nages))
mula[1,] <- Linf[1]*(1-exp(-k[1]*(ages-t0[1])))
mula[2,] <- Linf[2]*(1-exp(-k[2]*(ages-t0[2])))
cvla <- 0.15
for(a in 1:nages) {

  lref1 <- seq(mula[1,a]*(1-1.96*cvla),mula[1,a]*(1+1.96*cvla),length=50)
  dl1 <- dlnorm(lref1,log(mula[1,a]),sqrt(log(1+cvla^2)))
  dl1 <- dl1/sum(dl1)
  lref2 <- seq(mula[2,a]*(1-1.96*cvla),mula[2,a]*(1+1.96*cvla),length=50)
  dl2 <- dlnorm(lref2,log(mula[2,a]),sqrt(log(1+cvla^2)))
  dl2 <- dl2/sum(dl2) 
  m[1,a] <- sum(dl1*(1/(1+19^{-(lref1-sm50)/(sm95-sm50)})))
  m[2,a] <- sum(dl2*(1/(1+19^{-(lref2-sm50)/(sm95-sm50)})))
  w[1,a] <- sum(dl1*a.wl*lref1^b.wl)/1e+3
  w[2,a] <- sum(dl2*a.wl*lref2^b.wl)/1e+3

}

# define the movement model

pi12 <- 0.01
pi21 <- 0.05
pim[] <- matrix(c(1-pi12,pi12,pi21,1-pi21),ncol=2,byrow=T)

########################################
# 1. long-term eqm unfished population #
########################################

ny <- 500
Nbar <- array(dim=c(ny,nages,ns,R))
for(s in 1:ns) 
  for(r in 1:R) Nbar[1,,s,r] <- eta[r]*zeta[s]*R0*exp(-ages*M[1]) 
for(y in 2:ny) {

  for(s in 1:ns) 
    for(r in 1:R) Nbar[y,1,s,r] <- R0*eta[r]*zeta[s]
  for(s in 1:ns)
    for(a in 2:(nages-1)) 
      for(r in 1:R) Nbar[y,a,s,r] <- sum(pim[,r,a-1]*Nbar[y-1,a-1,s,]*exp(-M[a-1]))

  for(s in 1:ns)
    for(r in 1:R) Nbar[y,nages,s,r] <- sum(pim[,r,nages-1]*Nbar[y-1,nages-1,s,]*exp(-M[nages-1]))+sum(pim[,r,nages]*Nbar[y-1,nages,s,]*exp(-M[nages]))

}
  
############################################
# 2. Eqm short cut using simultaneous eqns #
############################################

Neqm <- array(dim=c(nages,ns,R))
for(s in 1:ns) Neqm[1,s,] <- R0*eta*zeta[s]
for(a in 2:(nages-1))
  for(s in 1:ns)
    for(r in 1:R) Neqm[a,s,r] <- sum(pim[,r,a-1]*Neqm[a-1,s,]*exp(-M[a-1]))

# solve simulatenous eqns for N in the plus group

a <- rep(NA,R)
B <- matrix(nrow=R,ncol=R)
for(s in 1:ns) {

  for(r in 1:R) a[r] <- sum(pim[,r,nages-1]*Neqm[nages-1,s,]*exp(-M[nages-1]))
  for(r in 1:R)
    for(ss in 1:R) B[ss,r] <- pim[ss,r,nages]*exp(-M[nages])
  npg <- t(a) %*% solve(diag(1,nrow=R)-B)
  Neqm[nages,s,] <- npg

}

#############################
# testing the popn dynamics #
#############################

B0 <- sum(apply(Neqm[,1,],1,sum)*w[1,]*m[1,])

# selectivity

nf <- 3 # NMR, SMR, AT longlines
rf <- c(1,2,2) # areas fisheries operate in

sl50 <- c(65,80,55)
sl95 <- c(75,110,65)

sel <- array(dim=c(nf,ns,nages))
for(f in 1:nf) {
  for(a in 1:nages) {

    lref1 <- seq(mula[1,a]*(1-1.96*cvla),mula[1,a]*(1+1.96*cvla),length=50)
    dl1 <- dlnorm(lref1,log(mula[1,a]),sqrt(log(1+cvla^2)))
    dl1 <- dl1/sum(dl1)
    lref2 <- seq(mula[2,a]*(1-1.96*cvla),mula[2,a]*(1+1.96*cvla),length=50)
    dl2 <- dlnorm(lref2,log(mula[2,a]),sqrt(log(1+cvla^2)))
    dl2 <- dl2/sum(dl2) 
    sel[f,1,a] <- sum(dl1*(1/(1+19^{-(lref1-sl50[f])/(sl95[f]-sl50[f])})))
    sel[f,2,a] <- sum(dl2*(1/(1+19^{-(lref2-sl50[f])/(sl95[f]-sl50[f])})))
  
  }
}

# total catch

Ctot <- 440
Cf <- c(100,140,200)

# SRR params

hh <- 0.75
rho <- B0/R0
alp <- 4*hh/(rho*(1-hh))
bet <- (5*hh-1)/(B0*(1-hh))

# other stuff

ny <- 35
N <- hya <- Cya <- array(dim=c(ny,nages,ns,R))
Cyfa <- hyfa <- array(dim=c(ny,nages,ns,nf))
tauf <- 0.5
h <- array(dim=c(ny,nf,ns))
phi <- m[1,]*w[1,]
S <- rep(NA,ny)
S[1] <- B0
N[1,,,] <- Neqm
Cyfa[1,,,] <- Cya[1,,,] <- hyfa[1,,,] <- 0
h[1,,] <- hya[1,,,] <- 0
hmax <- 0.8

for(y in 2:ny) {

  # recruitment

  Rtot <- alp*S[y-1]/(1+bet*S[y-1])
  for(s in 1:ns) 
    for(r in 1:R) N[y,1,s,r] <- Rtot*eta[r]*zeta[s]

  # population dynamics

 
    # from 2 to A-1

  for(s in 1:ns)
    for(a in 2:(nages-1))
      for(r in 1:R)   
        N[y,a,s,r] <- sum(pim[,r,a-1]*N[y-1,a-1,s,]*exp(-M[a-1])*(1-hya[y-1,a-1,s,])) 

    # plus group

  for(s in 1:ns)
    for(r in 1:R) 
        N[y,nages,s,r] <- sum(pim[,r,nages-1]*N[y-1,nages-1,s,]*exp(-M[nages-1])*(1-hya[y-1,nages-1,s,]))+sum(pim[,r,nages]*N[y-1,nages,s,]*exp(-M[nages])*(1-hya[y-1,nages,s,])) 

  # calculate harvest rates

  for(f in 1:nf) 
    for(s in 1:ns) 
      h[y,f,s] <- Cf[f]*zeta[s]/sum(N[y,,s,rf[f]]*exp(-tauf*M)*sel[f,s,]*w[s,])

  # calculate catch-at-age (factoring in hmax)

  for(s in 1:ns) {
    for(a in 1:nages) {
      for(f in 1:nf) hyfa[y,,s,f] <- h[y,f,s]*sel[f,s,]
       
      # check maxmim total harvest rate not exceeded
      # this is enforce within a given region

      for(r in 1:R) {

        # multiple fisheries in a given area?

        if(length(rf[rf == r]) > 1) {

          idx <- grep(r,rf) 
          if(sum(hyfa[y,a,s,idx]) > hmax) hyfa[y,a,s,idx] <- hmax * Cf[idx]/sum(Cf[idx])

          hya[y,a,s,r] <- sum(hyfa[y,a,s,idx])
          Cyfa[y,a,s,idx] <- hyfa[y,a,s,idx] * N[y,a,s,r]*exp(-tauf*M[a])
          Cya[y,a,s,r] <- sum(Cyfa[y,a,s,idx]) 

        } else {

          idx <- grep(r,rf)  
          if(hyfa[y,a,s,rf[r]] > hmax) hyfa[y,a,s,rf[r]] <- hmax
          hya[y,a,s,r] <- hyfa[y,a,s,rf[r]]
          Cyfa[y,a,s,idx] <- hyfa[y,a,s,idx] * N[y,a,s,r]*exp(-tauf*M[a]) 
          Cya[y,a,s,r] <- Cyfa[y,a,s,idx]

        }
      }
    }
  }

  # calculate female SSB across both areas

  S[y] <- sum(apply(N[y,,1,],1,sum)*phi)

}

ts.plot(S/B0,col='blue',lwd=2,ylim=c(0,1))
abline(h=0.5,col='magenta',lty=2)

##############################################
# generate some pseudo-data for use later on #
##############################################

# length frequency data for each fishery
# 1. create length partition
# 2. calculate P(l | a)
# 3. define effective sample sizes
# 4. simulate the length data 

# 1. length partition

lbins <- c(seq(0,20,by=10),seq(25,105,by=5),seq(110,190,by=10))
nbins <- length(lbins)-1
mulbins <- 0.5*(lbins[-c(1)]+lbins[-c(nbins+1)])

# 2. calculate P(l | a)

pla <- array(NA,dim=c(ns,nbins,nages))
dimnames(pla)[[1]] <- c('female','male')
dimnames(pla)[[3]] <- ages
dimnames(pla)[[2]] <- lbins[1:nbins]

sdla <- sqrt(log(1+cvla^2))
for(s in 1:ns) {
  for(a in 1:nages) {

    pl <- dlnorm(mulbins,log(mula[s,a]),sdla)
    pl <- pl/sum(pl)
    pla[s,,a] <- pl
  }
}

# sum to 1 check along length dimension

apply(pla,c(1,3),sum)

# faaaaancy ggplot of P(l | a) cos you know...

library(ggplot2)
df <- expand.grid(age=ages,length=lbins[1:nbins],sex=c('female','male'),p=NA)
for(a in ages) {
  
  df[df$age == a & df$sex == 'female','p'] <- pla['female',,a]
  df[df$age == a & df$sex == 'male','p'] <- pla['male',,a] 

}

ggplot(df,aes(age,as.factor(length),z=p))+geom_tile(aes(fill=p))+facet_wrap(~sex)+scale_fill_gradient(low="purple", high="cyan")+theme_bw()+ylab("Length")

# lattice version

library(lattice)
levelplot(p~age*length|sex,data=df,xlab='Age',ylab='Length',main='P(length | age)',scales=list(y=list(rot=45), x=list(rot=45)),cex=0.5)

# get length-conditional age probabilities for conversions

A <- matrix(scan("aerr_mat.dat"),ncol=nages,byrow=T)
adf <- expand.grid(a=ages,aa=ages,p=NA)
for(a in 1:nages) adf[adf$a == a,'p'] <- A[,a]
ggplot(adf,aes(a,aa,z=p))+geom_tile(aes(fill=p))+scale_fill_gradient(low="purple", high="cyan")+theme_bw()+xlab("true age")+ylab("observed age")

# steps:
# 1. Get "true" P(a | l,y,f) for each fishery and population overall
# 2. Get "observed" P(a | l,y,f) factoring in ageing error

# 1. true P(a | l,y,f) for each fishery ~ pla * pay

paly <- paaly <- array(dim=c(ny,nages,nbins,ns,nf))
palyt <- array(dim=c(ny,nages,nbins,ns,R))
pay <- array(dim=c(ny,nages,ns,nf)) # prior age-distro
ply <- array(dim=c(ny,nbins,ns,nf)) # length distro

# pay

for(s in 1:ns) 
  for(y in 1:ny) 
    for(f in 1:nf) {
      for(a in 1:nages) pay[y,a,s,f] <- N[y,a,s,rf[f]] * sel[f,s,a]

      pay[y,,s,f] <- pay[y,,s,f]/sum(pay[y,,s,f])

    }

# paly (and ply)

for(y in 1:ny)
  for(s in 1:ns)
    for(l in 1:nbins)
      for(f in 1:nf) {
        for(a in 1:nages) paly[y,a,l,s,f] <- pla[s,l,a] * pay[y,a,s,f]
        ply[y,l,s,f] <- sum(paly[y,,l,s,f])
        paly[y,,l,s,f] <- paly[y,,l,s,f] / ply[y,l,s,f]
      }

# get paaly (paly adjusted for ageing error)

for(y in 1:ny)
  for(s in 1:ns)
    for(l in 1:nbins)
      for(f in 1:nf) 
        for(a in 1:nages) paaly[y,a,l,s,f] <- sum(paly[y,,l,s,f] * A[a,])

# get palyt (true P(a | l) in actual popn)

for(y in 1:ny)
  for(s in 1:ns)
    for(r in 1:R) {

      pa <- N[y,,s,r]/sum(N[y,,s,r]) 
      for(l in 1:nbins) {
        
        palyt[y,,l,s,r] <- pla[s,l,] * pa
        pltmp <- sum(palyt[y,,l,s,r])
        palyt[y,,l,s,r] <- palyt[y,,l,s,r]/pltmp
      }
    }

# get predicted catch-at-length given catch-at-age

Cyfl <- array(dim=c(ny,nbins,ns,nf))
for(s in 1:ns)
  for(y in 1:ny)
    for(f in 1:nf)
      for(l in 1:nbins)
        Cyfl[y,l,s,f] <- sum(Cyfa[y,,s,f]*paly[y,,l,s,f])

# check total catches (numbers & biomass)

sum(apply(Cya,c(1,3),sum)-apply(Cyfl,c(1,3),sum))
Cfy <- array(dim=c(ny,nf,ns))
for(s in 1:ns) {
  for(y in 1:ny) {
    for(f in 1:nf) {
      Cfy[y,f,s] <- sum(Cyfa[y,,s,f]*w[s,])
    }
  }
}

apply(Cfy,1,sum)
Ctot

# 3. effective sample size-by-fishery

neff <- rep(250,nf)

# 4. Simulate LF data for each fishery

fnm <- c("NMR","SMR","AT")
nsim <- 1
lfdat <- expand.grid(y=1:ny,l=lbins[1:nbins],sex=c('female','male'),f=fnm,iter=1:nsim,c=NA)
for(n in 1:nsim)
  for(y in 1:ny)
    for(f in 1:nf) {

      if(y > 1) {
      ptmp <- Cyfl[y,,1,f]
      ptmp <- ptmp/sum(ptmp)
      lfdat[lfdat$y == y & lfdat$f == fnm[f] & lfdat$sex == 'female' & lfdat$iter == n,'c'] <- rmultinom(1,neff[f],ptmp)
      ptmp <- Cyfl[y,,2,f]
      ptmp <- ptmp/sum(ptmp) 
      lfdat[lfdat$y == y & lfdat$f == fnm[f] & lfdat$sex == 'male' & lfdat$iter == n,'c'] <- rmultinom(1,neff[f],ptmp)
      } else {

        lfdat[lfdat$y == y & lfdat$f == fnm[f] & lfdat$iter == n,'c'] <- 0

      }
    }


ggplot(subset(lfdat,sex=='female' & y >= ny-10),aes(x=l,y=c))+geom_line(colour='blue')+facet_grid(f~y)+ylab("numbers")+xlab("length")+theme(axis.text.x = element_text(angle = 60, hjust = 1))

ggplot(subset(lfdat,sex=='male' & y >= ny-10),aes(x=l,y=c))+geom_line(colour='blue')+facet_grid(f~y)+ylab("numbers")+xlab("length")+theme(axis.text.x = element_text(angle = 60, hjust = 1))

###################################################
# simulate age given length data for each fishery #
###################################################

# number of fish aged in total

nage <- c(200,100) # female, male
nasf <- matrix(nrow=ns,ncol=nf)
for(s in 1:ns) nasf[s,] <- round(nage[s] * Cf/Ctot)

# simulate some length-conditional age data for each fishery

alfdat <- expand.grid(y=1:ny,a=ages,l=lbins[1:nbins],sex=c("female","male"),f=fnm,iter=1:nsim,n=0)
zdf <- expand.grid(y=1:ny,l=lbins[1:nbins],sex=c("female","male"),f=fnm,iter=1:nsim,mua=NA)
lref <- 1:nbins

for(n in 1:nsim)
  for(y in 1:ny) {
    for(f in 1:nf) 

        if(y > 1) {
        
          # length samples first

          nn <- nasf[1,f]
          nl <- rmultinom(1,nn,ply[y,,1,f])
          idx <- lref[-c(grep(0,nl))]
          for(l in idx) {

            na <- rmultinom(1,nl[l],paly[y,,l,1,f])
            alfdat[alfdat$y == y & alfdat$l == lbins[l] & alfdat$f == fnm[f] & alfdat$sex == 'female' & alfdat$iter == n,'n'] <- na
            zdf[zdf$y == y & zdf$l == lbins[l] & zdf$f == fnm[f] & zdf$sex == 'female' & zdf$iter == n,'mua'] <- sum(na*ages/sum(na))
          }


          nn <- nasf[2,f]
          nl <- rmultinom(1,nn,ply[y,,2,f])
          idx <- lref[-c(grep(0,nl))] 
          for(l in idx) {

            na <- rmultinom(1,nl[l],paly[y,,l,2,f]) 
            alfdat[alfdat$y == y & alfdat$l == lbins[l] & alfdat$f == fnm[f] & alfdat$sex == 'male' & alfdat$iter == n,'n'] <- na
            zdf[zdf$y == y & zdf$l == lbins[l] & zdf$f == fnm[f] & zdf$sex == 'male' & zdf$iter == n,'mua'] <- sum(na*ages/sum(na)) 
          }
      
        } else {

          alfdat[alfdat$y == y & alfdat$l == lbins[l] & alfdat$f == fnm[f] & alfdat$iter == n,'n'] <- 0

        }

    cat("Year",y,"of",ny,"\n")
  }

ggplot(subset(zdf,sex=='female' & y >= ny-10),aes(x=l,y=mua))+geom_point(colour='blue',size=0.5)+facet_grid(f~y)+ylab("Mean age")+xlab("length")+theme(axis.text.x = element_text(angle = 60, hjust = 1))

ggplot(subset(zdf,sex=='male' & y >= ny-10),aes(x=l,y=mua))+geom_point(colour='blue',size=0.5)+facet_grid(f~y)+ylab("Mean age")+xlab("length")+theme(axis.text.x=element_text(angle=45))

#####################################################
# simulating tagging data via spatial Brownie model #
#####################################################

# size transition matrix

lvbinc <- function(lrel,tau,k,Linf) {

  return(max((Linf-lrel)*(1-exp(-k*tau)),0))
}

T <- array(dim=c(ns,nbins,nbins)) 
tau <- 1

for(s in 1:ns) {
  for(i in 1:nbins) {

    lx <- lbins[i]
    ly <- lbins[i+1]

    for(j in 1:nbins) {

      llj <- lbins[j]
      luj <- lbins[j+1]
          
      lli <- lx + lvbinc(lx,tau,k[s],Linf[s]) 
      lui <- ly + lvbinc(ly,tau,k[s],Linf[s]) 
      
      # need to work out Lebesgue measure of intersection 
      # of image and actual length bin / length bin

      if(lli > llj & lui < luj) {
        ptmp <- 1
      } else {
        tmp <- c(max(llj,lli),min(luj,lui))
        mu <- ifelse(tmp[1] < tmp[2],tmp[2]-tmp[1],0)
        nu <- lui-lli
        ptmp <- mu/nu 
      }
        
      T[s,i,j] <- ptmp
    }
  }
}

# size-based harvest & movement rates

hyfl <- array(dim=c(ny,nbins,ns,f))
for(s in 1:ns)
  for(y in 1:ny)
    for(f in 1:nf)
      for(l in 1:nbins)
        hyfl[y,l,s,f] <- sum(hyfa[y,,s,f]*paly[y,,l,s,f])

hyl <- array(dim=c(ny,nbins,ns,R))
for(s in 1:ns)
  for(y in 1:ny)
    for(l in 1:nbins)
      for(r in 1:R) {

        idx <- grep(r,rf)  
        if(length(rf[rf == r]) > 1) hyl[y,l,s,r] <- sum(hyfl[y,l,s,rf[idx]])
        else hyl[y,l,s,r] <- hyfl[y,l,s,rf[idx]]

      }

pil <- array(dim=c(R,R,nbins))
for(l in 1:nbins) pil[,,l] <- pim[,,1]


# number of release events (per year and per area

nrr <- 1 # per area
nry <- 1 # per year
nrel <- nrr*R*nry

# number of recapture events per release event

nrec <- 5

# years of releases

yrel <- ny-10+1:nry-1

# tag loss parameters

ptg <- 0.95 # instantanous tag loss
nutg <- 0.05 # rate of tag loss over time

# reporting rate

prep <- 0.94

# number of releases per tonne (2)

relpt <- 2
Trel <- Ctot*2

# sex ratio of the releases 

sext <- c(0.5,0.5)

# releases per sex per fishery

Tinit <- array(dim=c(ns,nf))
for(s in 1:ns) Tinit[s,] <- Trel * sext[s] * Cf/Ctot

tres <- list() # dummy list to bung everything into

Tyl <- omega <- array(0,dim=c(nrec+1,nbins,ns,R))
Ml <- rep(0.13,nbins)
nl <- rep(0,nbins)

for(n in 1:nsim) {

  # simulate length distro of releases

  kk <- 1
  for(rr in 1:nry) {
    for(ra in 1:R) { # release area
      for(s in 1:ns) {

        # loop over all fisheries in a given area

        idx <- grep(ra,rf)
        nl[] <- 0
        for(f in idx) nl[] <- nl[] + as.vector(rmultinom(1,Tinit[s,f],ply[yrel[rr],,s,f]))
        
        # simulate tag dynamics for cases where there are tags released

        for(lr in 1:nbins) {

          if(nl[lr] > 0) {

            Tyl[,,,] <- 0

            # set up objects for tag dynamics
        
            Tyl[1,lr,s,ra] <- nl[lr] 
            # tag dynamics for nrec+1 years
            for(tt in 2:(nrec+1)) {
        
              for(l in 1:nbins) 
              for(r in 1:R) {
                
                 ttmp <- 0
                if(tt == 2) {

                  for(rs in 1:R) ttmp <- ttmp+sum(pil[rs,r,]*Tyl[tt-1,,s,rs]*T[s,,l]*ptg*exp(-Ml-nutg^2)*(1-hyl[yrel+tt-2,,s,rs]))
                  Tyl[tt,l,s,r] <- ttmp

                } else {

                  for(rs in 1:R) ttmp <- ttmp+sum(pil[rs,r,]*Tyl[tt-1,,s,rs]*T[s,,l]*exp(-Ml-nutg^2)*(1-hyl[yrel+tt-2,,s,rs]))
                  Tyl[tt,l,s,r] <- ttmp 

                }
              } 
            }

            # store the useful stuff

            tres[[kk]] <- list()
            tres[[kk]][['rel_yr']] <- yrel[rr]
            tres[[kk]][['rel_len']] <- lr
            tres[[kk]][['rel_area']] <- ra
            tres[[kk]][['rel_sex']] <- s
            tres[[kk]][['nl']] <- nl
            tres[[kk]][['Tn']] <- Tyl

            kk <- kk+1

          }
        }
      }
    }
  }
}

# now generate some recaptures

prec <- array(dim=c(nrec,R))
psurv <- rep(NA,nrec+1)
ntag <- length(tres)

tdf <- expand.grid(sex=NA,yrel=yrel,lrel=NA,rrel=NA,rrec=1:R,recev=1:nrec,relev=1:ntag,ps=NA,prec=NA,nrel=0,nrectot=0,R=NA)

ntot <- 0
for(i in 1:ntag) {

  sr <- tres[[i]]$rel_sex
  lr <- tres[[i]]$rel_len
  rr <- tres[[i]]$rel_area
  yr <- tres[[i]]$rel_yr
  psurv[1] <- 1 
  Tx <- wx <- tres[[i]]$Tn
  wx[1,,sr,] <- Tx[1,,sr,]/sum(Tx[1,,sr,])

  for(j in 2:(nrec+1)) {
          
    wx[j,,sr,] <- Tx[j,,sr,]/sum(Tx[j,,sr,]) 
          
    # survival

    ptmp <- 0
    for(rs in 1:R) ptmp <- ptmp+sum(wx[j-1,,sr,rs]*exp(-Ml-nutg^2)*(1-hyl[j+yr-2,,sr,rs]))
    if(j == 2) {
      
      psurv[j] <- psurv[j-1]*ptmp*ptg

    } else {

      psurv[j] <- psurv[j-1]*ptmp 

    }

    # probability of recapture
  
    for(rs in 1:R) prec[j-1,rs] <- sum(wx[j,,sr,rs]*hyl[j+yr-1,,sr,rs])
    prec[j-1,] <- prec[j-1,] * psurv[j] * prep

  }

  # multinomial likelihood for recaptures

  ptmp <- as.vector(prec)
  psum <- sum(ptmp)
  ptmp <- c(ptmp,1-psum)
  nn <- Tx[1,lr,sr,rr]
  ntot <- ntot+nn
  Rec <- as.vector(rmultinom(1,nn,ptmp))
  Rmat <- prec
  for(r in 1:R)
    for(j in 1:nrec) Rmat[j,r] <- Rec[j+(r-1)*nrec]

  # write to the data frame

  sx <- ifelse(sr == 1,'female','male')
  tdf[tdf$relev == i,'sex'] <- sx
  tdf[tdf$relev == i,'nrel'] <- nn
  tdf[tdf$relev == i,'lrel'] <- lr 
  tdf[tdf$relev == i,'yrel'] <- yr
  tdf[tdf$relev == i,'rrel'] <- rr
  for(r in 1:R) {

    tdf[tdf$relev == i & tdf$rrec == r,'ps'] <- psurv[-c(1)]
    tdf[tdf$relev == i & tdf$rrec == r,'prec'] <- prec[,r]
    tdf[tdf$relev == i & tdf$rrec == r,'R'] <- Rmat[,r]

  }
  tdf[tdf$relev == i,'nrectot'] <- sum(Rmat)

  
  cat("Tagging event",i,"of",ntag,"\n")
}

# overall rate of return

sum(tdf$R)/ntot

# storage format
# matrix with dims: ntag, 5+R*rrec

Tdata <- matrix(nrow=length(unique(tdf$relev)),ncol=5+R*nrec)
colnames(Tdata) <- c("nrel","sex","lrel","rrel","yrel",rep("R",R*nrec))

ii <- 1
for(i in unique(tdf$relev)) {

  df <- subset(tdf,relev==i)
  Tdata[ii,1] <- df$nrel[1]
  Tdata[ii,2] <- ifelse(df$sex[1] == 'female',1,2)
  Tdata[ii,3] <- df$lrel[1]
  Tdata[ii,4] <- df$rrel[1]
  Tdata[ii,5] <- df$yrel[1]
  Tdata[ii,6:(nrec+5)] <- subset(df,rrec==1)$R
  Tdata[ii,(nrec+6):(nrec+10)] <- subset(df,rrec==2)$R
  ii <- ii+1

}

# write the four data sets:
# 1. Total catch biomass by fishery
# 2. LF by fishery
# 3. ALF by fishery
# 4. Tagging data

# 1.

ctab <- data.frame(year=1:ny,'NMR'=rep(Cf[1],ny),'SMR'=rep(Cf[2],ny),'AT'=rep(Cf[3],ny))
write.table(ctab,file='catch_biomass.dat',row.names=FALSE)

# 2. LF data (combined by sex)

lf <- array(dim=c(ny,nbins,nf))
lfdat$f <- as.character(lfdat$f)

for(yy in 1:ny)
  for(ff in 1:nf) {

    cdf <- subset(lfdat,yy==y & f == fnm[ff])
    zz <- subset(cdf,sex=='female')$c+subset(cdf,sex=='female')$c
    lf[yy,,ff] <- zz
  }

for(ff in 1:nf) {

  ttl <- paste(paste("lf",fnm[ff],sep="_"),"dat",sep=".")
  write(t(lf[,,ff]),file=ttl,ncolumns=nbins)
}

# 3. ALF data (by fishery and sex)

almat <- matrix(nrow=nbins,ncol=nages)
for(ss in 1:ns) {

  sx <- ifelse(ss == 1,'female','male') 

  for(ff in 1:nf) {

    fx <- fnm[ff]

    ttl <- paste(paste(paste("alf",fx,sep="_"),sx,sep="_"),'dat',sep=".")
    for(yy in 1:ny) {

      
      df <- subset(alfdat,y == yy & sex == sx & f == fx)
      for(ll in 1:nbins) almat[ll,] <- subset(df,l == lbins[ll])$n
      if(yy == 1) {
        
        write(t(almat),file=ttl,ncolumns=nages)
        
      } else {

        write(t(almat),file=ttl,ncolumns=nages,append=TRUE)

      }
    }
  }
}

# 4. tagging data

write.table(Tdata,file='sim_tag_data.dat',row.names=FALSE)

# save it all

save.image("test_pdyn_datagen.rda")

