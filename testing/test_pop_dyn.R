######################################################
# testing population dynamics calculations ###########
######################################################
# R Hillary CSIRO 2018 ###############################
######################################################

ages <- 0:50
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

B0 <- sum(Neqm[,1,1]*w[1,]*m[1,])+sum(Neqm[,1,2]*w[1,]*m[1,])

# selectivity

nf <- 3 # NMR, SMR, AT longlines
rf <- c(1,2,2) # areas fisheries operate in

sl50 <- c(65,95,55)
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
N <- C <- array(dim=c(ny,nages,ns,R))
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

