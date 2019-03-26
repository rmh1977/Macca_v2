#####################################################
# m1a: 2nd iteration of new Macca assessment ########
#####################################################
# R. Hillary, J. Day CSIRO 2019 #####################
#####################################################

require(TMB)

compile("m1b.cpp")
dyn.load("m1b.so")

################
# set up stuff #
################

yrs <- 1985:2018 # back 10 years before catches begin
ny <- length(yrs)
iyrs <- 0:(ny-1) # TMB model index years
names(iyrs) <- as.character(yrs)
ages <- 1:52
na <- length(ages)
nr <- 2 # area 1 (N), area 2 (S)
ns <- 2
lbins <- c(0,10,20,30,seq(35,140,5),seq(150,190,10))
nbins <- length(lbins)-1
nl <- length(lbins)-1
mulbins <- 0.5*(lbins[-c(1)]+lbins[-c(nl+1)])
nf <- 5 # ATT, NVT, ATL, NMRL, SMRL longline
ffnm <- c("ATT","NVT","ATL","NMRL","SMRL")
rf <- c(2,1,2,1,2)-1 # areas fisheries operate in C++ offset to zero
names(rf) <- ffnm
nrecmax <- 7 # maxium years of recaptures post release year

######################################
# set up the data: ###################
# 1. catch biomass ###################
# 2. catch timing (fraction of year) #
# 3. length frequency data ###########
# 4. age-length frequency data #######
# 5. tagging data ####################
######################################

# catch biomass and timing by fishery

datdir <- "../../../data/toothfish2019/Data/2019/processed/"
fnm <- paste(datdir,"CatchBiomass.RData",sep="")
load(fnm)
gnm <- paste(datdir,"CatchTiming.RData",sep="")
load(gnm)

C <- tau <- matrix(0,nrow=ny,ncol=nf)
rownames(C) <- rownames(tau) <- as.character(yrs)
colnames(C) <- colnames(tau) <- ffnm

for(i in 1:dim(cbdf)[1]) {

  yy <- as.character(cbdf$Year[i])
  ff <- as.character(cbdf$Fleet[i])
  cc <- cbdf$Catch_t[i]
  tt <- taudf$tau[i]
  C[yy,ff] <- cc
  tau[yy,ff] <- tt

}

# modify 2001 and 2010 missing catch to avoid zero probs in tag likelihood

C['2001',] <- rep(0.01,nf)
C['2010','NVT'] <- C['2010','NMRL'] <- 0.01  

# fisheries with logistic (0) or double-logistic (1) selectivty functions

self <- c(1,1,0,0,0)
names(self) <- ffnm
nlg <- length(self[self==0])
ndl <- length(self[self==1])
flg <- c(2,3,4)
fdl <- c(0,1)

# length frequency data (using nhauls as initial ESS)

fnm <- paste(datdir,"LF_nhaul.dat",sep="")
lfdat <- matrix(scan(fnm),ncol=nbins+3,byrow=T)

iy <- unname(iyrs[as.character(lfdat[,1])])
iff <- lfdat[,2]-1
LFcov <- unname(cbind(iy,iff))
LF <- lfdat[,-c(1,2)]

# age-given-length frequency data

fnm <- paste(datdir,"ALF.dat",sep="")
alfdat <- matrix(scan(fnm),ncol=na+5,byrow=T)
iy <- unname(iyrs[as.character(alfdat[,1])])
iff <- alfdat[,2]-1
isx <- alfdat[,3]-1
il <- alfdat[,4]-1
ALFcov <- unname(cbind(iy,iff,isx,il))
ALF <- alfdat[,-c(1:4)]

# tagging data

fnm <- paste(datdir,"tags.dat",sep="")
tdf <- read.table(fnm,header=T)
tmin <- 5 # minimum number of releases
tdf <- subset(tdf,T >= tmin)
relyrs <- sort(unique(tdf$relyr))

tlist <- list()
ii <- 1 
for(y in relyrs) {

  yrecmax <- min(max(yrs),y+nrecmax)
  xdf <- subset(tdf,relyr == y & recyr <= yrecmax)
  ll <- sort(unique(xdf$lrel))
  for(l in ll) {

    ydf <- subset(xdf,lrel == l)

    rrel <- unique(ydf$relarea)
    if(length(rrel) == 1) {

      rrel <- ydf$relarea[1] 
      Rtmp1 <- subset(ydf,recarea==1)$R
      Rtmp2 <- subset(ydf,recarea==2)$R 
      tlist[[ii]] <- list()
      tlist[[ii]]$T <- ydf$T[1]
      tlist[[ii]]$yrel <- y
      tlist[[ii]]$l <- l 
      tlist[[ii]]$rrel <- rrel
      tlist[[ii]]$nrec <- yrecmax-y
      tlist[[ii]]$R <- unname(cbind(Rtmp1,Rtmp2))
      ii <-ii+1

    } else {

      for(rr in rrel) {

        zdf <- subset(ydf,relarea == rr)
        Rtmp1 <- subset(zdf,recarea==1)$R
        Rtmp2 <- subset(zdf,recarea==2)$R 
        tlist[[ii]] <- list()
        tlist[[ii]]$T <- zdf$T[1]
        tlist[[ii]]$yrel <- y
        tlist[[ii]]$l <- l 
        tlist[[ii]]$rrel <- rr
        tlist[[ii]]$nrec <- yrecmax-y
        tlist[[ii]]$l <- l 
        tlist[[ii]]$R <- unname(cbind(Rtmp1,Rtmp2))
        ii <-ii+1 

      }
    }
  }
}

ntagev <- length(tlist)
tagrel <- tagnrec <- rep(NA,ntagev)
tagcov <- matrix(ncol=3,nrow=ntagev)
tagrec <- array(0,dim=c(ntagev,nrecmax,nr))
for(i in 1:ntagev) {

  tagrel[i] <- tlist[[i]]$T
  tagnrec[i] <- tlist[[i]]$nrec
  tagcov[i,] <- c(tlist[[i]]$l-1,tlist[[i]]$rrel-1,unname(iyrs[as.character(tlist[[i]]$yrel)]))
  nrectrue <- dim(tlist[[i]]$R)[1]
  tagrec[i,1:nrectrue,] <- tlist[[i]]$R

}

# reporting rates over time

rrdf <- read.table("tagrep.dat",header=T)

# save the data

fnm <- paste(paste("MIdata_maxrecev",nrecmax,sep="_"),"rda",sep=".")
save.image(fnm)

##########################################
# set up the data and fixed input object #
##########################################

A <- matrix(scan("aerr_mat.dat"),ncol=na,byrow=T)
aref <- c(5,20) # reference ages for Schnute growth parameterisation
dms <- c(ny,na,nl,nf,nr) # key dims: y, a, l, f, r
nobs <- c(dim(LF)[1],dim(ALF)[1],dim(tagrec)[1],nrecmax) # number of obs.
dmsel <- c(nlg,ndl) # number of selectivity types
phiLF <- rep(1.1,nf) # fisheries-specifc OD factor for LF data
phiT <- 1.5 # OD factor for tagging data
hh <- 0.75 # steepness
sigmaR <- 0.4 # stk-rec variability
M <- 0.13 # natural mortality
zeta <- c(0.5,0.5) # birth sex-ratio
hmax <- 0.8 # maximum harvest-rate per age-class and area
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
             nobs=nobs,
             rf=unname(rf),
             dmsel=dmsel,
             flg=flg,
             fdl=fdl,
             mulbins=mulbins,
             lbins=lbins,
             phiLF=phiLF,
             phiT=phiT,
             hh=hh,
             sigmaR=sigmaR,
             M=M,
             zeta=zeta,
             hmax=hmax,
             psi=psi,
             nu=nu,
             prep=rrdf$RR,
             sdl=sdl)

######################################
# set up the initial parameter guess #
######################################

R0 <- 1e+6
eta <- c(0.45,0.55)
pim <- c(1-0.01,0.01,0.05,1-0.05)
selparslg <- c(60,10,70,20,90,30) # ATL, NMRL, SMRL
selparsdl <- c(50,5,5,55,55,8,8,35) # ATT, NVT
l1 <- c(52,51) # L(a1) female and male
l2 <- c(119,102) # L(a2) female and male
k <- c(0.08,0.05) # VB-k female and male
sdla <- c(0.15,0.15) # SD (mula) female and male

pars <- list(logR0 = log(R0),
             epsR = rep(0,ny),
             logeta = log(eta),
             logpim = log(pim),
             selparslg = log(selparslg),
             selparsdl = log(selparsdl),
             logl1 = log(l1),
             logl2 = log(l2),
             logk = log(k),
             logsdla = log(sdla))

###############################
# start the estimation phases #
###############################

# phase 1: R0 

map1 <- list(epsR = rep(factor(NA),ny),
             logeta = rep(factor(NA),nr),
             logpim = rep(factor(NA),nr*nr),
             selparslg = rep(factor(NA),nlg*2),
             selparsdl = rep(factor(NA),ndl*4),
             logl1 = rep(factor(NA),ns),
             logl2 = rep(factor(NA),ns),
             logk = rep(factor(NA),ns),
             logsdla = rep(factor(NA),ns))

obj1 <- MakeADFun(data=data,parameters=pars,map=map1,DLL="m1b")

rep <- obj1$rep()



