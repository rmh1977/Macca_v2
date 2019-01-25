######################################################
# run toothfish model on simulated data set ##########
######################################################
# R. Hillary CSIRO O & A 2018 ########################
######################################################

require(TMB)

compile("mtest_schnute.cpp")
dyn.load("mtest_schnute.so")

##########################
# set up the data object #
##########################

# dimensions, age and length partitions

ny <- 35
ages <- 1:52
na <- length(ages)
nr <- 2 # area 1 (N), area 2 (S)
ns <- 2
lbins <- c(seq(0,20,by=10),seq(25,105,by=5),seq(110,190,by=10))
nl <- length(lbins)-1
mulbins <- 0.5*(lbins[-c(1)]+lbins[-c(nl+1)])
nf <- 3 # NMR, SMR, AT longline
rf <- c(1,2,2)-1 # areas fisheries operate in C++ offset to zero
nrecmax <- 5 # maxium five years of recaptures post release year

# catch biomass

C <- read.table("catch_biomass.dat",header=T)
cb <- as.matrix(C[,-1])

# fraction of year gone before catch taken (Jan 1 being year 1)

tau <- matrix(0.5,nrow=ny,ncol=nf)

# length frequency data per fishery and sex-aggregated

LF <- array(dim=c(ny,nl,nf))
lfAT <- matrix(scan("lf_AT.dat"),nrow=ny,byrow=T)
lfSMR <- matrix(scan("lf_SMR.dat"),nrow=ny,byrow=T)
lfNMR <- matrix(scan("lf_NMR.dat"),nrow=ny,byrow=T)
LF[,,1] <- lfNMR
LF[,,2] <- lfSMR
LF[,,3] <- lfAT

# flag for years which don't have data (0) and years which do (1)

lfflag <- matrix(nrow=ny,ncol=nf)
nlftot <- apply(LF,c(1,3),sum)
for(y in 1:ny)
  for(f in 1:nf) lfflag[y,f] <- ifelse(nlftot[y,f] > 0,1,0)

# age-given-length data - females and males separately

ALFf <- ALFm <- array(dim=c(ny,nl,na,nf))
alfATf <- matrix(scan("alf_AT_female.dat"),ncol=na,byrow=T)
alfSMRf <- matrix(scan("alf_SMR_female.dat"),ncol=na,byrow=T)
alfNMRf <- matrix(scan("alf_NMR_female.dat"),ncol=na,byrow=T)
for(y in 1:ny) {

  ALFf[y,,,1] <- alfNMRf[((y-1)*nl+1):(y*nl),]
  ALFf[y,,,2] <- alfSMRf[((y-1)*nl+1):(y*nl),]
  ALFf[y,,,3] <- alfATf[((y-1)*nl+1):(y*nl),] 

}

alfATm <- matrix(scan("alf_AT_male.dat"),ncol=na,byrow=T)
alfSMRm <- matrix(scan("alf_SMR_male.dat"),ncol=na,byrow=T)
alfNMRm <- matrix(scan("alf_NMR_male.dat"),ncol=na,byrow=T)
for(y in 1:ny) {

  ALFm[y,,,1] <- alfNMRm[((y-1)*nl+1):(y*nl),]
  ALFm[y,,,2] <- alfSMRm[((y-1)*nl+1):(y*nl),]
  ALFm[y,,,3] <- alfATm[((y-1)*nl+1):(y*nl),] 

}

alfflag <- almflag <- matrix(nrow=ny,ncol=nf)
nalftotf <- apply(ALFf,c(1,4),sum)
nalftotm <- apply(ALFm,c(1,4),sum)
for(y in 1:ny)
  for(f in 1:nf) {
    
    alfflag[y,f] <- ifelse(nalftotf[y,f] > 0,1,0)
    almflag[y,f] <- ifelse(nalftotm[y,f] > 0,1,0) 
  }

# tagging data

tagdat <- read.table("sim_tag_data.dat",header=T)
tagrel <- tagdat[,1]
tagcov <- as.matrix(tagdat[,2:5])
colnames(tagcov) <- NULL
tagcov[,] <- tagcov[,]-1 # zero offset in C++ brah
nt <- dim(tagcov)[1]
tagrec <- array(dim=c(nt,nrecmax,nr))
for(t in 1:nt) {

  rtmp <- as.numeric(tagdat[t,-c(1:5)])
  for(r in 1:nr) tagrec[t,,r] <- rtmp[(1+(r-1)*nrecmax):(r*nrecmax)]

}

# ageing error matrix

A <- matrix(scan("aerr_mat.dat"),nrow=na,byrow=T)

# dimensions vector

dms <- c(ny,na,nl,nf,nr,nt,nrecmax)

# fixed quantities

hh <- 0.75
sigmaR <- 0.05
M <- 0.13
zeta <- c(0.5,0.5)
hmax <- 0.8
kap <- 25
psi <- 0.95
nu <- 0.05
prep <- 0.94
sdl <- 0.05

# load up save one

load("mtest.rda")

# bring together the data, dimensions and fixed quantities

a1 <- 5 # young ref. age for growth
a2 <- 20 # old ref. age for growth

data <- list(C = cb,
             tau = tau,
             LF = LF,
             lfflag = lfflag,
             ALFf = ALFf,
             alfflag = alfflag,
             ALFm = ALFm,
             almflag = almflag,
             tagrel = tagrel,
             tagcov = tagcov,
             tagrec = tagrec,
             A = A,
             aref = c(a1,a2),
             dms = dms,
             rf = rf,
             mulbins = mulbins,
             lbins = lbins,
             hh = hh,
             sigmaR = sigmaR,
             M = M,
             zeta = zeta,
             hmax = hmax,
             psi = psi,
             nu = nu,
             prep = prep,
             sdl = sdl)

###########################################
# prepare the initial parameter estimates #
###########################################

Linf <- c(165,165)
k <- c(0.057,0.052)
t0 <- c(-0.19,-0.34)
l1 <- Linf*(1-exp(-k*(a1-t0)))
l2 <- Linf*(1-exp(-k*(a2-t0)))

pars <- list(logR0 = log(1e+6),
             epsR = rep(1e-6,ny),
             logeta = log(c(0.45,0.55)),
             logpim = log(c(1-0.01,0.01,0.05,1-0.05)),
             selpars = log(c(65,10,80,30,55,10)),
             logl1 = log(l1),
             logl2 = log(l2),
             logk = log(c(0.057,0.052)),
             logsdla = log(c(sqrt(log(1+0.15^2)),sqrt(log(1+0.15^2)))))

#####################################
# phasing process: the map function #
#####################################

# phase 1: lnR0 only

map1 <- list(epsR = rep(factor(NA),ny),
             logeta = rep(factor(NA),nr),
             logpim = rep(factor(NA),nr*nr),
             selpars = rep(factor(NA),nf*2),
             logl1 = rep(factor(NA),ns),
             logl2 = rep(factor(NA),ns),
             logk = rep(factor(NA),ns),
             logsdla = rep(factor(NA),ns))

obj1 <- MakeADFun(data=data,parameters=pars,map=map1,DLL="mtest_schnute")

rep <- obj1$rep()

res1 <- do.call("optim",obj1)
rep1 <- obj1$rep()
sdrep1 <- sdreport(obj1)

# phase 2: spatial recruitment fraction

map2 <- list(epsR = rep(factor(NA),ny),
             logpim = rep(factor(NA),nr*nr),
             selpars = rep(factor(NA),nf*2),
             logl1 = rep(factor(NA),ns),
             logl2 = rep(factor(NA),ns), 
             logk = rep(factor(NA),ns),
             logsdla = rep(factor(NA),ns))

pars2 <- pars
pars2$logeta <- jitter(pars2$logeta)

obj2 <- MakeADFun(data=data,parameters=pars2,map=map2,DLL="mtest_schnute")
 
res2 <- do.call("optim",obj2)
rep2 <- obj2$rep()
sdrep2 <- sdreport(obj2)

# phase 3: fix spatial recruitment, estimate movement

map3 <- list(epsR = rep(factor(NA),ny),
             logeta = rep(factor(NA),nr),
             selpars = rep(factor(NA),nf*2),
             logl1 = rep(factor(NA),ns),
             logl2 = rep(factor(NA),ns), 
             logk = rep(factor(NA),ns),
             logsdla = rep(factor(NA),ns))

pars3 <- pars

obj3 <- MakeADFun(data=data,parameters=pars3,map=map3,DLL="mtest_schnute")
 
res3 <- do.call("optim",obj3)
rep3 <- obj3$rep()
sdrep3 <- sdreport(obj3)

# phase 4: estimate selectivity, R0, fixed everything else

map4 <- list(epsR = rep(factor(NA),ny),
             logeta = rep(factor(NA),nr),
             logpim = rep(factor(NA),nr*nr), 
             logl1 = rep(factor(NA),ns),
             logl2 = rep(factor(NA),ns), 
             logk = rep(factor(NA),ns),
             logsdla = rep(factor(NA),ns))

pars4 <- pars
  
obj4 <- MakeADFun(data=data,parameters=pars4,map=map4,DLL="mtest_schnute")
 
res4 <- do.call("optim",obj4)
rep4 <- obj4$rep()
sdrep4 <- sdreport(obj4)

# phase 5: estimate growth (k & variance only)

map5 <- list(epsR = rep(factor(NA),ny),
             logeta = rep(factor(NA),nr),
             logpim = rep(factor(NA),nr*nr), 
             selpars = rep(factor(NA),nf*2),
             logl1 = rep(factor(NA),ns),
             logl2 = rep(factor(NA),ns))


pars5 <- pars
  
obj5 <- MakeADFun(data=data,parameters=pars5,map=map5,DLL="mtest_schnute")
 
res5 <- do.call("optim",obj5)
rep5 <- obj5$rep()
sdrep5 <- sdreport(obj5)




