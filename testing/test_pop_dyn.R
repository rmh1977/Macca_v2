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
N <- array(dim=c(ny,nages,ns,R))
C <- array(dim=c(ny,nages,ns,nf))
tauf <- 0.5
h <- array(dim=c(ny,nf,ns))
phi <- m[1,]*w[1,]
S <- rep(NA,ny)
S[1] <- B0
N[1,,,] <- Neqm
C[1,,,] <- 0
h[1,,] <- 0
hmax <- 0.8

for(y in 2:ny) {

  # recruitment

  Rtot <- alp*S[y-1]/(1+bet*S[y-1])
  for(s in 1:ns) 
    for(r in 1:R) N[y,1,s,r] <- Rtot*eta[r]*zeta[s]

  # population dynamics

  for(s in 1:ns) {
    for(a in 2:(nages-1))
      for(r in 1:R) {
        
        hsum <- 0
        for(f in 1:nf) {
          if(rf[f] == r) hsum <- hsum+h[y-1,f,s]*sel[f,s,a-1]
        }
        hsum <- min(hmax,hsum)
        N[y,a,s,r] <- sum(pim[,r,a-1]*N[y-1,a-1,s,]*exp(-M[a-1])*(1-hsum)) 

      }

    # plus group

    for(r in 1:R) {
        
        hsum1 <- hsum2 <- 0
        for(f in 1:nf) {
          if(rf[f] == r) {
            
            hsum1 <- hsum1+h[y-1,f,s]*sel[f,s,nages-1]
            hsum2 <- hsum2+h[y-1,f,s]*sel[f,s,nages] 

          }
        }
        hsum1 <- min(hmax,hsum1)
        hsum2 <- min(hmax,hsum2)
        N[y,nages,s,r] <- sum(pim[,r,nages-1]*N[y-1,nages-1,s,]*exp(-M[nages-1])*(1-hsum1))+sum(pim[,r,nages]*N[y-1,nages,s,]*exp(-M[nages])*(1-hsum2)) 

      } 

  }

  # calculate harvest rates

  for(f in 1:nf) 
    for(s in 1:ns) 
      h[y,f,s] <- Cf[f]*zeta[s]/sum(N[y,,s,rf[f]]*exp(-tauf*M)*sel[f,s,]*w[s,])

  # calculate catch-at-age (factoring in hmax)

  for(s in 1:ns) {
    for(a in 1:nages) {
      for(f in 1:nf) {
        htmp <- h[y,f,s]*sel[f,s,a]
        htmp <- min(hmax,htmp)
        C[y,a,s,f] <- htmp * N[y,a,s,rf[f]]*exp(-tauf*M[a])
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

# get predicted catch-at-length given catch-at-age

Cl <- array(dim=c(ny,nbins,ns,nf))
for(s in 1:ns)
  for(y in 1:ny)
    for(f in 1:nf)
      for(l in 1:nbins)
        Cl[y,l,s,f] <- sum(C[y,,s,f]*pla[s,l,])

# check total catches (numbers & biomass)

sum(apply(C,c(1,3),sum)-apply(Cl,c(1,3),sum))
Cfy <- array(dim=c(ny,nf,ns))
for(s in 1:ns) {
  for(y in 1:ny) {
    for(f in 1:nf) {
      Cfy[y,f,s] <- sum(C[y,,s,f]*w[s,])
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
      ptmp <- Cl[y,,1,f]
      ptmp <- ptmp/sum(ptmp)
      lfdat[lfdat$y == y & lfdat$f == fnm[f] & lfdat$sex == 'female' & lfdat$iter == n,'c'] <- rmultinom(1,neff[f],ptmp)
      ptmp <- Cl[y,,2,f]
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

A <- matrix(scan("aerr_mat.dat"),ncol=nages,byrow=T)
adf <- expand.grid(a=ages,aa=ages,p=NA)
for(a in 1:nages) adf[adf$a == a,'p'] <- A[,a]
ggplot(adf,aes(a,aa,z=p))+geom_tile(aes(fill=p))+scale_fill_gradient(low="purple", high="cyan")+theme_bw()+xlab("true age")+ylab("observed age")

# steps:
# 1. Get "true" P(a | l,y,f) for each fishery
# 2. Get "observed" P(a | l,y,f) factoring in ageing error
# 3. Generate samples from the correct distro

# 1. true P(a | l,y,f) for each fishery ~ pla * pay

paly <- paaly <- array(dim=c(ny,nages,nbins,ns,nf))
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

# 2. get paaly (paly adjusted for ageing error)

for(y in 1:ny)
  for(s in 1:ns)
    for(l in 1:nbins)
      for(f in 1:nf) 
        for(a in 1:nages) paaly[y,a,l,s,f] <- sum(paly[y,,l,s,f] * A[a,])

# 3. simulate some length-conditional age data for each fishery

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

ggplot(subset(zdf,sex=='male'),aes(x=l,y=mua))+geom_point(colour='blue',size=0.5)+facet_grid(f~y)+ylab("Mean age")+xlab("length")+theme(axis.text.x=element_text(angle=45))

save.image("test_pdyn_datagen.rda")

