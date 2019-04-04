#####################################################
# m1c: 4th iteration of new Macca assessment ########
#####################################################
# R. Hillary, J. Day CSIRO 2019 #####################
#####################################################

require(TMB)

compile("m1c.cpp")
dyn.load("m1c.so")

source("m1_utils.R")

load("MIdata_maxrecev_7.rda")

##########################################
# set up the data and fixed input object #
##########################################

A <- matrix(scan("aerr_mat.dat"),ncol=na,byrow=T)
aref <- c(5,20) # reference ages for Schnute growth parameterisation
dms <- c(ny,na,nl,nf,nr) # key dims: y, a, l, f, r
datflag <- c(1,1,1) # LF, ALF, tags
recyr <- c(1,28) # min and max years for recruitment estimation
nobs <- c(dim(LF)[1],dim(ALF)[1],dim(tagrec)[1],nrecmax) # number of obs.
phiLF <- rep(2,nf) # fisheries-specifc OD factor for LF data
phiT <- 1.5 # OD factor for tagging data
hh <- 0.75 # steepness
sigmaR <- 0.25 # stk-rec variability
M <- 0.13 # natural mortality
zeta <- c(0.5,0.5) # birth sex-ratio
hmax <- 0.8 # maximum harvest-rate per age-class and area
prtotmax <- 0.9 # maximum total chance of recovering a tag release
psi <- 1 # immediate tag shedding probability
nu <- 0 # continuous rate of tag shedding
sdl <- 0.05 # random variability (SD) in growth increments

data <- list(C=unname(C),
             tau=unname(tau),
             LFcov=LFcov,
             LF=LF,
             ALFcov=ALFcov,
             ALF=ALF,
             tagrel=tagrel,
             tagnrec=tagnrec,
             tagcov=tagcov,
             tagrec=tagrec,
             A=A,
             aref=aref,
             dms=dms,
             datflag=datflag,
             recyr=recyr,
             nobs=nobs,
             rf=unname(rf),
             mulbins=mulbins,
             lbins=lbins,
             phiLF=phiLF,
             phiT=phiT,
             hh=hh,
             sigmaR=sigmaR,
             M=M,
             zeta=zeta,
             hmax=hmax,
             prtotmax=prtotmax,
             psi=psi,
             nu=nu,
             prep=rrdf$RR,
             sdl=sdl)

######################################
# set up the initial parameter guess #
######################################

R0 <- 8e+5
ilogeta <- logit(0.3)
ilogpim <- logit(c(0.99,0.95))
selparsATT <- c(55,5,30) # dr95 is fixed at 3
selparsNVT <- c(20,0.2857)
selparsATL <- c(60,10)
selparsNMRL <- c(70,20)
selparsSMRL <- c(90,30)
l1 <- c(49,49) # L(a1) female and male
l2 <- c(116,102) # L(a2) female and male
k <- c(0.048,0.071) # VB-k female and male
sdla <- c(0.15,0.14) # SD (mula) female and male
nyrec <- length(recyr[1]:recyr[2])

pars1 <- list(logR0 = log(R0),
             epsR = rep(0,nyrec),
             ilogeta = ilogeta,
             ilogpim = ilogpim,
             selparsATT = log(selparsATT),
             selparsNVT = log(selparsNVT), 
             selparsATL = log(selparsATL), 
             selparsNMRL = log(selparsNMRL), 
             selparsSMRL = log(selparsSMRL), 
             logl1 = log(l1),
             logl2 = log(l2),
             logk = log(k),
             logsdla = log(sdla))

###############################
# start the estimation phases #
###############################

# phase 1: R0

map1 <- list(epsR = rep(factor(NA),nyrec),
             ilogeta = factor(NA),
             ilogpim = rep(factor(NA),nr),
             selparsATT = rep(factor(NA),3), 
             selparsNVT = rep(factor(NA),2),
             selparsATL = rep(factor(NA),2), 
             selparsNMRL = rep(factor(NA),2), 
             selparsSMRL = rep(factor(NA),2), 
             logl1 = rep(factor(NA),ns),
             logl2 = rep(factor(NA),ns),
             logk = rep(factor(NA),ns),
             logsdla = rep(factor(NA),ns))

obj1 <- MakeADFun(data=data,parameters=pars1,map=map1,DLL="m1c")
res1 <- do.call("optim",obj1)
rep1 <- obj1$rep()
sdrep1 <- sdreport(obj1)
summary(sdrep1)

# phase 2: selectivity

map2 <- list(epsR = rep(factor(NA),nyrec),
             ilogeta = factor(NA),
             ilogpim = rep(factor(NA),nr), 
             logl1 = rep(factor(NA),ns),
             logl2 = rep(factor(NA),ns),
             logk = rep(factor(NA),ns),
             logsdla = rep(factor(NA),ns))

pars2 <- list(logR0 = res1$par[['logR0']],
             epsR = rep(0,nyrec),
             ilogeta = ilogeta,
             ilogpim = ilogpim,
             selparsATT = log(selparsATT),
             selparsNVT = log(selparsNVT), 
             selparsATL = log(selparsATL), 
             selparsNMRL = log(selparsNMRL), 
             selparsSMRL = log(selparsSMRL), 
             logl1 = log(l1),
             logl2 = log(l2),
             logk = log(k),
             logsdla = log(sdla))

obj2 <- MakeADFun(data=data,parameters=pars2,map=map2,DLL="m1c")
res2 <- do.call("optim",obj2)
rep2 <- obj2$rep()
sdrep2 <- sdreport(obj2)
summary(sdrep2)

# phase 3: spatial recruitment 

map3 <- list(epsR = rep(factor(NA),nyrec),
             ilogpim = rep(factor(NA),nr), 
             logl1 = rep(factor(NA),ns),
             logl2 = rep(factor(NA),ns),
             logk = rep(factor(NA),ns),
             logsdla = rep(factor(NA),ns))

pars3 <- list(logR0 = res2$par[['logR0']],
             epsR = rep(0,nyrec),
             ilogeta = ilogeta,
             ilogpim = ilogpim,
             selparsATT = unname(res2$par[grep('selparsATT',names(res2$par))]),
             selparsNVT = unname(res2$par[grep('selparsNVT',names(res2$par))]), 
             selparsATL = unname(res2$par[grep('selparsATL',names(res2$par))]), 
             selparsNMRL = unname(res2$par[grep('selparsNMRL',names(res2$par))]), 
             selparsSMRL = unname(res2$par[grep('selparsSMRL',names(res2$par))]), 
             logl1 = log(l1),
             logl2 = log(l2),
             logk = log(k),
             logsdla = log(sdla))

obj3 <- MakeADFun(data=data,parameters=pars3,map=map3,DLL="m1c")
res3 <- do.call("optim",obj3)
rep3 <- obj3$rep()
sdrep3 <- sdreport(obj3)
summary(sdrep3)

# phase 4: movement 

map4 <- list(epsR = rep(factor(NA),nyrec),
             logl1 = rep(factor(NA),ns),
             logl2 = rep(factor(NA),ns),
             logk = rep(factor(NA),ns),
             logsdla = rep(factor(NA),ns))

pars4 <- list(logR0 = res3$par[['logR0']],
             epsR = rep(0,nyrec),
             ilogeta = res3$par[['ilogeta']],
             ilogpim = ilogpim,
             selparsATT = unname(res3$par[grep('selparsATT',names(res3$par))]),
             selparsNVT = unname(res3$par[grep('selparsNVT',names(res3$par))]), 
             selparsATL = unname(res3$par[grep('selparsATL',names(res3$par))]), 
             selparsNMRL = unname(res3$par[grep('selparsNMRL',names(res3$par))]), 
             selparsSMRL = unname(res3$par[grep('selparsSMRL',names(res3$par))]), 
             logl1 = log(l1),
             logl2 = log(l2),
             logk = log(k),
             logsdla = log(sdla))

obj4 <- MakeADFun(data=data,parameters=pars4,map=map4,DLL="m1c")
res4 <- do.call("optim",obj4)
rep4 <- obj4$rep()
sdrep4 <- sdreport(obj4)
summary(sdrep4)

# phase 5: recruitment 

map5 <- list(logl1 = rep(factor(NA),ns),
             logl2 = rep(factor(NA),ns),
             logk = rep(factor(NA),ns),
             logsdla = rep(factor(NA),ns))

pars5 <- list(logR0 = res4$par[['logR0']],
             epsR = rep(0,nyrec),
             ilogeta = res4$par[['ilogeta']],
             ilogpim =unname(res4$par[grep('ilogpim',names(res4$par))]),
             selparsATT = unname(res4$par[grep('selparsATT',names(res4$par))]),
             selparsNVT = unname(res4$par[grep('selparsNVT',names(res4$par))]), 
             selparsATL = unname(res4$par[grep('selparsATL',names(res4$par))]), 
             selparsNMRL = unname(res4$par[grep('selparsNMRL',names(res4$par))]), 
             selparsSMRL = unname(res4$par[grep('selparsSMRL',names(res4$par))]), 
             logl1 = log(l1),
             logl2 = log(l2),
             logk = log(k),
             logsdla = log(sdla))

obj5 <- MakeADFun(data=data,parameters=pars5,map=map5,DLL="m1c")
res5 <- do.call("optim",obj5)
rep5 <- obj5$rep()
sdrep5 <- sdreport(obj5)
summary(sdrep5)

# plots

plot.lenfreq(LF,LFcov,rep5,'SMRL',phiLF[3])
plot.agelendat(ALF,ALFcov,rep5,'SMRL','male')
plot.tagfits.relyr(rep5,tagrel,tagnrec,tagcov,tagrec,nrecmax)
plot.tagfits.relyr.tot(rep5,tagrel,tagnrec,tagcov,tagrec,nrecmax)
plot.tagfits.recyr.tot(rep5,tagrel,tagnrec,tagcov,tagrec,nrecmax)
plot.tagfits.recyrarea(rep5,tagrel,tagnrec,tagcov,tagrec,nrecmax)

######################
# MCMC using tmbstan #
######################

library(tmbstan)

mcmc <- tmbstan(obj5, chains=1,control=list(max_treedepth=12))
traceplot(mcmc, pars=names(obj5$par), inc_warmup=TRUE)
traceplot(mcmc, pars=names(obj5$par), inc_warmup=FALSE)
pairs(mcmc, pars=names(obj5$par)[1])
parlist <- extract(mcmc)

