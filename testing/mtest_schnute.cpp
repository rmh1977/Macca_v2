//////////////////////////////////////////////////////
// test model for revised MITF assessment ////////////
// 1. Catch biomass by fishery ///////////////////////
// 2. LF data by fishery (sex aggregated) ////////////
// 3. ALF data by sex and fishery ////////////////////
// 4. Tag data by sex and length /////////////////////
//////////////////////////////////////////////////////
// R. Hillary, J. Day CSIRO 2018 /////////////////////
//////////////////////////////////////////////////////

#include <TMB.hpp>

using namespace density;

// square function
template <class Type>
Type square(Type x){return x*x;}

// inverse-logistic function
template <class Type>
Type ilogit(Type x){return Type(1)/(Type(1)+exp(-x));}

// useful function for transition matrix 
template <class Type>
Type ptrans(Type llxi,Type llxj,Type luxi,Type luxj) {
  
  Type tmpl,tmpu,ru,rl;
  tmpl = CppAD::CondExpGt(llxi,llxj,llxi,llxj);
  tmpu = CppAD::CondExpLt(luxi,luxj,luxi,luxj);
  ru = CppAD::CondExpLt(tmpl,tmpu,tmpu-tmpl,Type(0));
  rl = luxi-llxi;
  return(ru/rl);

}

template<class Type>
Type objective_function<Type>::operator() ()
{

  // data call ins & dimensions

  DATA_ARRAY(C); // catch biomass by fishery (y,f)
  DATA_ARRAY(tau); // timing of catches through "year" (y,f)
  DATA_ARRAY(LF); // LF by fishery (y, l  f)
  DATA_ARRAY(lfflag); // flag for years each fishery has LF data (y,f)
  DATA_ARRAY(ALFf); // female ALF by fishery (y, l, a, f) 
  DATA_ARRAY(alfflag); // flag for years fishery has female ALF data (y,f)
  DATA_ARRAY(ALFm); // male ALF by fishery (y, l, a, f)
  DATA_ARRAY(almflag); // flag for years fishery has male ALF data (y,f)
  DATA_VECTOR(tagrel); // vector of tag release numbers
  DATA_IMATRIX(tagcov); // tag release event covariates
  DATA_ARRAY(tagrec); // tag recaptures (nt,nrec,R)
  DATA_ARRAY(A); // ageing error matrix (a,aa)
  DATA_VECTOR(aref); // a1 and a2 for Schnute parameterisation

  // dimensions

  DATA_IVECTOR(dms);
  int ny = dms[0];
  int na = dms[1];
  int nl = dms[2];
  int nf = dms[3];
  int nr = dms[4];
  int nt = dms[5];
  int nrec = dms[6];
  int ns = 2; // 0 = female, 1 = male
  DATA_IVECTOR(rf); // area for each fishery

  ////////////////////////////////
  // fixed variables and priors //
  ////////////////////////////////

  DATA_VECTOR(mulbins); // mean length per length bin
  DATA_VECTOR(lbins); // boundaries of the length partition
  DATA_SCALAR(hh); // steepness
  DATA_SCALAR(sigmaR); // S-R variability
  DATA_SCALAR(M); // natural mortality
  DATA_VECTOR(zeta); // sex ratio at birth
  DATA_SCALAR(hmax); // maximum hrate per age-class & area
  DATA_SCALAR(psi); // immediate tag shedding/mortality probability
  DATA_SCALAR(nu); // continuous rate of tag shedding
  DATA_SCALAR(prep); // tag reporting rate
  DATA_SCALAR(sdl); // random variability (SD) in growth increments

  // parameters

  PARAMETER(logR0); // mean unfished total recruitment
  PARAMETER_VECTOR(epsR); // recruitment deviations
  PARAMETER_VECTOR(logeta); // spatial recruitment parameter
  PARAMETER_VECTOR(logpim); // movement matrix parameters
  PARAMETER_VECTOR(selpars); // selectivity parameters for each fleet
  PARAMETER_VECTOR(logl1); // log-L1(a1)
  PARAMETER_VECTOR(logl2); // log-L2(a2) 
  PARAMETER_VECTOR(logk); // log-k
  PARAMETER_VECTOR(logsdla); // log-SD in mean length-at-age

  // weight-at-length & female maturity-at-length

  Type awl = Type(4.4e-6);
  Type bwl = Type(3.14);
  Type ml50 = Type(139);
  Type ml95 = Type(185);

  // variable declarations

  Type objf,nlogl,nlogp;
  Type Rtot,alp,bet,B0;

  // population dynamics stuff

  array<Type> N(ny,na,ns,nr);
  array<Type> sexrat(ny,nl,ns,nr); 
  matrix<Type> veqm(nr,1);
  matrix<Type> npgeqm(1,nr);
  matrix<Type> Meqm(nr,nr);
  array<Type> SSB(ny,ns,nr);
  vector<Type> SSBftot(ny);
  array<Type> Cya(ny,na,ns,nr);
  array<Type> Cyfa(ny,na,ns,nf);
  array<Type> hyf(ny,ns,nf);
  array<Type> hyfa(ny,na,ns,nf);
  array<Type> hya(ny,na,ns,nr);
  array<Type> hylf(ny,nl,ns,nf);
  array<Type> hyl(ny,nl,ns,nr);
  array<Type> sel(na,ns,nf); 
  matrix<Type> m(na,ns);
  matrix<Type> w(na,ns);
  matrix<Type> mula(na,ns); 
  matrix<Type> pim(nr,nr);
  vector<Type> eta(nr);
  vector<Type> l1(ns);
  vector<Type> l2(ns);
  vector<Type> Linf(ns);
  vector<Type> k(ns);
  vector<Type> t0(ns);
  vector<Type> psig(ns);
  vector<Type> omegag(ns);
  vector<Type> ldelg(ns);

  // length and age stuff
  
  array<Type> pla(ns,nl,na); // P(l | a,s) - length-at-age distro
  array<Type> pay(ny,na,ns,nf); // P(a | f,s) - prior age distro
  array<Type> paly(ny,na,nl,ns,nf); // P(a | l,f,s) - truq age-at-length distro
  array<Type> paaly(ny,na,nl,ns,nf); // P(a | l,f,s,a) - adjusted for age error
  array<Type> ply(ny,nl,ns,nf); // P(l | f,s) - length distro
  array<Type> Gamma(ns,nl,nl); // growth transition matrix

  /////////////////////////
  // population dynamics //
  /////////////////////////

  // size-at-age (Schnute formulation)

  for(int s=0;s<ns;s++) {

    l1(s) = exp(logl1(s));
    l2(s) = exp(logl2)(s);
    ldelg(s) = l2(s)-l1(s);
    k(s) = exp(logk(s));
    omegag(s) = Type(1)-exp(-k(s)*(aref(1)-aref(0)));
    psig(s) = (l1(s)*omegag(s))/ldelg(s);
    Linf(s) = l1(s)+ldelg(s)/omegag(s);
    t0(s) = aref(0)-log(psig(s)+Type(1))/k(s);

  }

  for(int s=0;s<ns;s++) 
    for(int a=0;a<na;a++) 
      mula(a,s) = l1(s)+(l2(s)-l1(s))*(Type(1)-exp(-k(s)*(double(a+1)-aref(0))))/(Type(1)-exp(-k(s)*(aref(1)-aref(0))));
  
  // maturity, selectivity and weight-at-age

  Type ll,lu,lref,ldel,psum,summ,sumw,sumf,s50,s95;
  vector<Type> dl(25);
  vector<Type> cvla(ns);
  vector<Type> sdla(ns);
  sdla = exp(logsdla);
  cvla = sqrt(exp(square(sdla))-Type(1));

  // P(l | a,s)

  for(int s=0;s<ns;s++) {
    for(int a=0;a<na;a++) { 
      
      psum = Type(0);
      for(int l=0;l<nl;l++) psum += dnorm(log(mulbins(l)),log(mula(a,s)),sdla(s),false);
      for(int l=0;l<nl;l++) pla(s,l,a) = dnorm(log(mulbins(l)),log(mula(a,s)),sdla(s),false)/psum;
    }
  }

  for(int s=0;s<ns;s++) {
    for(int a=0;a<na;a++) {
          
      ll = mula(a,s)-Type(1.96)*cvla(s)*mula(a,s);
      //if(ll < Type(0)) ll = mula(a,s)*Type(0.01);
      ll = CppAD::CondExpLt(ll,Type(0),mula(a,s)*Type(0.01),ll);
      lu = mula(a,s)+Type(1.96)*cvla(s)*mula(a,s);
      ldel = (lu-ll)/Type(25);
      for(int i=0;i<25;i++) {

          lref = ll+double(i)*ldel;
          dl(i) = dnorm(log(lref),log(mula(a,s)),exp(logsdla(s)),false);
         
      }

      psum = dl.sum();
      for(int i=0;i<25;i++) dl(i) /= psum; 

      // maturity and weight-at-age

      summ = sumw = Type(0);
      for(int i=0;i<25;i++) {
          
        lref = ll+double(i)*ldel;
        summ += dl(i)/(Type(1)+pow(Type(19),-(lref-ml50)/(ml95-ml50)));
        sumw += dl(i)*awl*pow(lref,bwl);

      }

      m(a,s) = summ;
      w(a,s) = sumw/Type(1000);

      // selectivity-at-age

      for(int f=0;f<nf;f++) {

        s50 = exp(selpars(0+f*2));
        s95 = s50+exp(selpars(1+f*2)); // selpars(1,f) = log(s95-s50) ensuring +ve

        sumf = Type(0);
        for(int i=0;i<25;i++) {
          
          lref = ll+double(i)*ldel;
          sumf += dl(i)/(Type(1)+pow(Type(19),-(lref-s50)/(s95-s50)));

        }

        sel(a,s,f) = sumf; 
      }
    }
  }

  // initial numbers-at-age

  eta = exp(logeta);
  summ = eta.sum();
  for(int r=0;r<nr;r++) eta(r) /= summ; // sum to 1 condition for spatial rec
  for(int rr=0;rr<nr;rr++) 
    for(int r=0;r<nr;r++) pim(rr,r) = exp(logpim(r+rr*nr));

  for(int rr=0;rr<nr;rr++) {

    psum = Type(0);
    for(int r=0;r<nr;r++) psum += pim(rr,r);
    for(int r=0;r<nr;r++) pim(rr,r) /= psum; // rows of P(i -> j) sum to 1 

  }

  Type R0 = exp(logR0);
  for(int s=0;s<ns;s++) 
    for(int r=0;r<nr;r++) N(0,0,s,r) = R0 * eta(r) * zeta(s); // recruitment
  for(int a=1;a<na-1;a++) {
    for(int s=0;s<ns;s++) {
      for(int r=0;r<nr;r++) {
        
        N(0,a,s,r) = Type(0);
        for(int rr=0;rr<nr;rr++) N(0,a,s,r) += pim(rr,r) * N(0,a-1,s,rr) * exp(-M);
      }
    }
  }

  matrix<Type> Minv(nr,nr);
  matrix<Type> veqmt(1,nr);
  for(int s=0;s<ns;s++) {
    for(int r=0;r<nr;r++) {

      veqm(r,0) = Type(0);
      for(int rr=0;rr<nr;rr++) {
        
        veqm(r,0) += pim(rr,r) * N(0,na-2,s,rr) * exp(-M);
        if(r == rr) Meqm(rr,r) = Type(1)-pim(rr,r)*exp(-M);
        if(r != rr) Meqm(rr,r) = -pim(rr,r)*exp(-M);

      }
    }

    Minv = Meqm.inverse();
    veqmt = veqm.transpose();
    npgeqm = veqmt*Minv; // solution to simultaneous eqns
    for(int r=0;r<nr;r++) N(0,na-1,s,r) = npgeqm(0,r);
  }

  // unfished (female) SSB for S-R relationship

  B0 = Type(0);
  for(int a=0;a<na;a++) 
    for(int r=0;r<nr;r++) B0 += N(0,a,0,r) * m(a,0) * w(a,0);
  SSBftot(0) = B0;

  for(int s=0;s<ns;s++) {
    for(int r=0;r<nr;r++) {

      SSB(0,s,r) = Type(0);
      for(int a=0;a<na;a++) SSB(0,s,r) += N(0,a,s,r) * m(a,s) * w(a,s);

    }
  }

  // other S-R parameters

  Type rho = B0/R0;
  alp = Type(4)*hh/(rho*(Type(1)-hh));
  bet = (Type(5)*hh-Type(1))/(B0*(Type(1)-hh));

  // initial harvest rates

  Type Xtmp,csum,htmp,hsum;
  int fcnt;

  // calculate overall harvest rates
    
  for(int f=0;f<nf;f++) {
    for(int s=0;s<ns;s++) {

      Xtmp = Type(0);
      for(int a=0;a<na;a++) Xtmp += N(0,a,s,rf(f))*exp(-tau(0,f)*M)*sel(a,s,f)*w(a,s);
      hyf(0,s,f) = C(0,f)*zeta(s)/Xtmp;

    }
  } 

  for(int s=0;s<ns;s++) {
    for(int a=0;a<na;a++) {

      // hmax is max. TOTAL harvest rate per age-class 
      // enforeced within a given region

      for(int r=0;r<nr;r++) {
        
        fcnt = 0; 
        for(int f=0;f<nf;f++) if(rf(f) == r) fcnt++; // multiple fisheries in area?
        if(fcnt > 1) {
            
          hsum = csum = Type(0);
          for(int f=0;f<nf;f++) {

            if(rf(f) == r) {
                
              hsum += hyf(0,s,f) * sel(a,s,f);
              csum += C(0,f);

            }
          }

          for(int f=0;f<nf;f++) {
                 
            if(rf(f) == r) {

              htmp = hyf(0,s,f) * sel(a,s,f);
              hyfa(0,a,s,f) = CppAD::CondExpLt(hsum,hmax,htmp,hmax*C(0,f)/csum);

            }   
          }
          

        } else {

          htmp = hyf(0,s,rf(r)) * sel(a,s,rf(r));
          hyfa(0,a,s,rf(r)) = CppAD::CondExpLt(htmp,hmax,htmp,hmax);

        }
      }

      // catch-at-age & total harvest rates per area

      for(int f=0;f<nf;f++) Cyfa(0,a,s,f) = hyfa(0,a,s,f)*exp(-tau(0,f)*M)*N(0,a,s,rf(f));
              
      for(int r=0;r<nr;r++) {

        hsum = csum = Type(0);
        for(int f=0;f<nf;f++) {

          if(rf(f) == r) {
              
            hsum += hyfa(0,a,s,f);
            csum += Cyfa(0,a,s,f);

          }
        }

        hya(0,a,s,r) = hsum;
        Cya(0,a,s,r) = csum;
      }
    }
  } 

  
  // annual loop
 
  for(int y=1;y<ny;y++) {

    // recruitment

    Rtot = alp*SSBftot(y-1)/(Type(1)+bet*SSBftot(y-1)); // expected value
    for(int s=0;s<ns;s++)
      for(int r=0;r<nr;r++)
        N(y,0,s,r) = Rtot * eta(r) * zeta(s) * exp(epsR(y)-square(sigmaR)/Type(2));

    for(int s=0;s<ns;s++) {

      // ages 2 to A-1

      for(int a=1;a<na-1;a++) {
        for(int r=0;r<nr;r++) {

          N(y,a,s,r) = Type(0);
          for(int rr=0;rr<nr;rr++) N(y,a,s,r) += pim(rr,r)*N(y-1,a-1,s,rr)*exp(-M)*(Type(1)-hya(y-1,a-1,s,rr));

        }
      }

      // plus group

      for(int r=0;r<nr;r++) {

          N(y,na-1,s,r) = Type(0);
          for(int rr=0;rr<nr;rr++) N(y,na-1,s,r) += pim(rr,r)*N(y-1,na-1,s,rr)*exp(-M)*(Type(1)-hya(y-1,na-1,s,rr))+pim(rr,r)*N(y-1,na-2,s,rr)*exp(-M)*(Type(1)-hya(y-1,na-2,s,rr));

        } 

      // SSB
 
      for(int r=0;r<nr;r++) {

        SSB(y,s,r) = Type(0);
        for(int a=0;a<na;a++) SSB(y,s,r) += N(y,a,s,r) * m(a,s) * w(a,s);

      }
    }

    SSBftot(y) = Type(0);
    for(int r=0;r<nr;r++) SSBftot(y) += SSB(y,0,r);

    // calculate overall harvest rates
    
    for(int f=0;f<nf;f++) {
      for(int s=0;s<ns;s++) {

        Xtmp = Type(0);
        for(int a=0;a<na;a++) Xtmp += N(y,a,s,rf(f))*exp(-tau(y,f)*M)*sel(a,s,f)*w(a,s);
        hyf(y,s,f) = C(y,f)*zeta(s)/Xtmp;

      }
    }

    // harvest rates & catch-at-age (modified if needed given hmax)

    for(int s=0;s<ns;s++) {
      for(int a=0;a<na;a++) {

        // hmax is max. TOTAL harvest rate per age-class 
        // enforced within a given region

        for(int r=0;r<nr;r++) {
        
          fcnt = 0; 
          for(int f=0;f<nf;f++) if(rf(f) == r) fcnt++; // multiple fisheries in area?
          if(fcnt > 1) {
            
            hsum = csum = Type(0);
            for(int f=0;f<nf;f++) {

              if(rf(f) == r) {
                
                hsum += hyf(y,s,f) * sel(a,s,f);
                csum += C(y,f);

             }
           }

            for(int f=0;f<nf;f++) {
                 
              if(rf(f) == r) {

                htmp = hyf(y,s,f) * sel(a,s,f);
                hyfa(y,a,s,f) = CppAD::CondExpLt(hsum,hmax,htmp,hmax*C(0,f)/csum);

              }   
            }
          

          } else {

            htmp = hyf(y,s,rf(r)) * sel(a,s,rf(r));
            hyfa(y,a,s,rf(r)) = CppAD::CondExpLt(htmp,hmax,htmp,hmax);

          }
        }

        // catch-at-age & total harvest rates per area

        for(int f=0;f<nf;f++) Cyfa(y,a,s,f) = hyfa(y,a,s,f)*exp(-tau(y,f)*M)*N(y,a,s,rf(f));
              
        for(int r=0;r<nr;r++) {

          hsum = csum = Type(0);
          for(int f=0;f<nf;f++) {

            if(rf(f) == r) {
              
              hsum += hyfa(y,a,s,f);
              csum += Cyfa(y,a,s,f);

            }
          }

          hya(y,a,s,r) = hsum;
          Cya(y,a,s,r) = csum;
        }
      }
    } 

  }

  // length-based sex ratio (useful for later on)

  Type nsum;
  vector<Type> nprior(na);
  vector<Type> paltmp(na);
  for(int y=0;y<ny;y++) {
    for(int l=0;l<nl;l++) {
      for(int r=0;r<nr;r++) {

        // calculate "true" P(a | l) in population (not fishery)

        for(int s=0;s<ns;s++) {

          nsum = Type(0);
          for(int a=0;a<na;a++) nsum += N(y,a,s,r);
          for(int a=0;a<na;a++) nprior(a) = N(y,a,s,r)/nsum;
          psum = Type(0);
          for(int a=0;a<na;a++) psum += pla(s,l,a)*nprior(a);
          for(int a=0;a<na;a++) paltmp(a) = pla(s,l,a)*nprior(a)/psum;
          nsum = Type(0);
          for(int a=0;a<na;a++) nsum += N(y,a,s,r) * paltmp(a);
          sexrat(y,l,s,r) = nsum;

        }

        nsum = Type(0); 
        for(int s=0;s<ns;s++) nsum += sexrat(y,l,s,r);
        for(int s=0;s<ns;s++) sexrat(y,l,s,r) = sexrat(y,l,s,r)/nsum;

      }
    }
  }

  ///////////////////////////////////////////////////
  // age-given length and length frequency distros //
  ///////////////////////////////////////////////////

  for(int y=0;y<ny;y++) {
    for(int s=0;s<ns;s++) {
      for(int f=0;f<nf;f++) {
        
        // 1. prior age distro ~ N(a,y) * sel(a)

        psum = Type(0);
        for(int a=0;a<na;a++) psum += N(y,a,s,rf(f))*sel(a,s,f);
        for(int a=0;a<na;a++) pay(y,a,s,f) = N(y,a,s,rf(f))*sel(a,s,f)/psum; 
        
        // 2. age-given-length distros
        
        for(int l=0;l<nl;l++) {

          // "true" age-given-length

          psum = Type(0);
          for(int a=0;a<na;a++) psum += pla(s,l,a)*pay(y,a,s,f);
          ply(y,l,s,f) = psum;
          for(int a=0;a<na;a++) paly(y,a,l,s,f) = pla(s,l,a)*pay(y,a,s,f)/psum;

          // adjusted for ageing error
        
          for(int a=0;a<na;a++) {
          
            psum = Type(0);
            for(int aa=0;aa<na;aa++) psum += paly(y,aa,l,s,f)*A(a,aa);
            paaly(y,a,l,s,f) = psum;

          }
        }
      }
    }
  }

  //////////////////////////////////////////////////////////////
  // length-based harvest rates needed for tagging likelihood //
  //////////////////////////////////////////////////////////////

  // fleet specific first (as P(a | l) depends on fleet)

  for(int s=0;s<ns;s++) {
    for(int y=0;y<ny;y++) {
      for(int f=0;f<nf;f++) {
        for(int l=0;l<nl;l++) {

          htmp = Type(0);
          for(int a=0;a<na;a++) htmp += hyfa(y,a,s,f)*paly(y,a,l,s,f);
          hylf(y,l,s,f) = htmp;

        }
      }
    }
  }

  // total harvest rate-at-length for each area

  for(int s=0;s<ns;s++) {
    for(int y=0;y<ny;y++) {
      for(int l=0;l<nl;l++) {
        for(int r=0;r<nr;r++) {
          
          htmp =  Type(0);
          for(int f=0;f<nf;f++) {

            if(rf(f) == r) htmp += hylf(y,l,s,f);

          }

          hyl(y,l,s,r) = htmp;

        }
      }
    }
  }

  //////////////////////////////
  // growth transition matrix //
  //////////////////////////////

  vector<Type> epsl(25);
  vector<Type> pl(25);
  Type el,eu,gdel,gtmp,ptmp,lx,ly,llj,luj,lli,lui,xl,xu,xt;
  el = Type(0)-1.96*sdl;
  eu = Type(0)+1.96*sdl;

  for(int i=0;i<25;i++) {
    
    epsl(i) = el+double(i)*(eu-el)/Type(25);
    dl(i) = dnorm(epsl(i),Type(0),sdl,false);

  }

  psum = dl.sum();
  for(int i=0;i<25;i++) dl(i) /= psum;

  for(int s=0;s<ns;s++) {
    for(int i=0;i<nl;i++) {

      lx = lbins(i);
      ly = lbins(i+1);

      for(int j=0;j<nl;j++) {
      
        llj = lbins(j);
        luj = lbins(j+1);

        psum = Type(0);
        for(int n=0;n<25;n++) {
          
          gtmp = (Linf(s)-lx)*(Type(1)-exp(-k(s)*Type(1)));
          gdel = CppAD::CondExpLt(gtmp,Type(0),Type(0),gtmp); 
          lli = lx+gdel*exp(epsl(n));
          gtmp = (Linf(s)-ly)*(Type(1)-exp(-k(s)*Type(1)));
          gdel = CppAD::CondExpLt(gtmp,Type(0),Type(0),gtmp); 
          lui = ly+gdel*exp(epsl(n));

          // looks funny but needed for TMB/CppAD handling of if/else/max/min shite

          xl = CppAD::CondExpGt(lli,llj,Type(1),Type(0));
          xu = CppAD::CondExpLt(lui,luj,Type(1),Type(0));
          xt = xl+xu;
          ptmp = CppAD::CondExpEq(xt,Type(2),Type(1),ptrans(lli,llj,lui,luj));
          psum += ptmp*dl(n);

        }
  
        Gamma(s,i,j) = psum;
      }
    }
  }

  ////////////////
  // likelihood //
  ////////////////

  // length frequency data (multinomial)
  // not split by sex 
  // compute mean LF across both sexes
  // uses sex-ratio by length

  vector<Type> nlf(nf);
  Type neff,plhat;
  nlf.setZero();

  for(int f=0;f<nf;f++) {
    for(int y=0;y<ny;y++) {

      if(lfflag(y,f) == 1) {
 
        for(int l=0;l<nl;l++) {

          plhat = Type(0);
          for(int s=0;s<ns;s++) plhat += ply(y,l,s,f) * sexrat(y,l,s,rf(f));
          nlf(f) -= LF(y,l,f) * log(plhat);
        }

      }
    }
  }

  // age-given-length data (multinomial)
  
  vector<Type> nalff(nf);
  vector<Type> nalfm(nf); 
  nalff.setZero();
  nalfm.setZero(); 

  for(int f=0;f<nf;f++) {
    for(int y=0;y<ny;y++) { 
      for(int s=0;s<ns;s++) {
        for(int l=0;l<nl;l++) {

          if(s == 0) { // females

            if(alfflag(y,f) == 1) {
 
              for(int a=0;a<na;a++) nalff(f) -= ALFf(y,l,a,f)*log(paaly(y,a,l,s,f));

            }

          } 

          if(s == 1) { // males

            if(almflag(y,f) == 1) {
 
              for(int a=0;a<na;a++) nalfm(f) -= ALFm(y,l,a,f)*log(paaly(y,a,l,s,f));

            } 

          }

        }
      }
    }
  }

  // tagging data (Brownie - multinomial)

  vector<Type> ntag(nt);
  ntag.setZero();
  Type nT,nR;
  array<Type> om(nrec+1,nl,nr);
  matrix<Type> Rec(nrec+1,nr);
  matrix<Type> prec(nrec+1,nr);
  array<Type> prhat(nt,nrec+1,nr);
  prhat.setZero();
  vector<Type> ps(nrec+1);
  int stag,ltag,rtag,ytag;

  for(int t=0;t<nt;t++) {

    nT = tagrel(t); // number of releases
    stag = tagcov(t,0); // sex of release 
    ltag = tagcov(t,1); // length bin of release
    rtag = tagcov(t,2); // region of release
    ytag = tagcov(t,3); // year of release
    for(int y=0;y<nrec+1;y++)
      for(int r=0;r<nr;r++) Rec(y,r) = y == 0 ? Type(0) : tagrec(t,y-1,r);
    nR = Type(0);
    for(int y=0;y<nrec+1;y++)
      for(int r=0;r<nr;r++) nR += Rec(y,r); 

    // tag dynamics

    om.setZero();
    om(0,ltag,rtag) = nT*psi; // releases in right length bin and region
    for(int y=1;y<nrec+1;y++) {
      for(int r=0;r<nr;r++) {
        for(int l=0;l<nl;l++) {

            om(y,l,r) = Type(0);
            for(int rr=0;rr<nr;rr++)
              for(int kk=0;kk<nl;kk++) om(y,l,r) += om(y-1,kk,rr)*Gamma(stag,kk,l)*exp(-M-square(nu))*(Type(1)-hyl(ytag+y-1,kk,stag,rr))*pim(rr,r);
          
        }
      }
    }

    for(int y=0;y<nrec+1;y++) {

      nsum = Type(0);
      for(int r=0;r<nr;r++) 
        for(int l=0;l<nl;l++) nsum += om(y,l,r);

      for(int r=0;r<nr;r++) 
        for(int l=0;l<nl;l++) om(y,l,r) /= nsum; // sum to 1 across size and regions

    }
         
    // tag survival probabilty

    
    ps(0) = Type(1);
    for(int y=1;y<nrec+1;y++) {
    
      ptmp = Type(0);
      for(int r=0;r<nr;r++) 
        for(int l=0;l<nl;l++) 
          ptmp += om(y-1,l,r)*exp(-M-square(nu))*(Type(1)-hyl(ytag+y-1,l,stag,r));

      ps(y) = y == 1 ? ps(y-1)*psi*ptmp : ps(y-1)*ptmp; // instant tag loss/mortality bit 

    }

    // tag recapture probability

    psum = Type(0);
    for(int y=0;y<nrec+1;y++) {
      for(int r=0;r<nr;r++) {
    
        ptmp = Type(0);
        for(int l=0;l<nl;l++) ptmp += om(y,l,r)*hyl(ytag+y,l,stag,r);
        prec(y,r) = ps(y) * prep * ptmp;
        prhat(t,y,r) = prec(y,r);
        if(y > 0) psum += prec(y,r);

      }
    }

    // neg. log-likelihood of the recapture history of tagging "event"
  
    for(int y=1;y<nrec+1;y++) 
      for(int r=0;r<nr;r++) ntag(t) -= Rec(y,r)*log(prec(y,r));
      
    ntag(t) -= (nT-nR)*log(Type(1)-psum); // P(never recapturing a tag)
  
  }

  ////////////
  // priors //
  ////////////
  
  // recruitment deviation prior

  nlogp = -dnorm(epsR,-square(sigmaR)/Type(2),sigmaR,true).sum();

  //////////////////////////////
  // total objective function //
  //////////////////////////////

  // length data + age-given-length data + tags + priors

  nlogl = nlf.sum()+nalff.sum()+nalfm.sum()+ntag.sum();

  //objf = Type(0);
  objf = nlogl+nlogp;

  // ADREPORT section

  // REPORT section
 
  REPORT(paly);
  REPORT(paaly);
  REPORT(ply);
  REPORT(pla);
  REPORT(prhat);
  REPORT(sexrat);
  REPORT(Gamma); 
  REPORT(sel);
  REPORT(Linf);
  REPORT(t0);
  REPORT(mula);
  REPORT(pim);
  REPORT(B0);
  REPORT(R0); 
  REPORT(N);
  REPORT(nlf);
  REPORT(nalff);
  REPORT(nalfm);
  REPORT(ntag);
  REPORT(nlogp);
  REPORT(nlogl);
  REPORT(SSBftot);

  return(objf);
}

