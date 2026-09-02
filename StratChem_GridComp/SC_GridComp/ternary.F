      SUBROUTINE ternary(tin,ptot,qh2o,qhno3,qh2so4,ws,wn,v,parthno3,
     &                   err,denssts,pi)

c... subroutine:   ternary
c
c... modified by:  David B. Considine
c...               
c...               NASA Goddard Space Flight Center
c...               Greenbelt, MD 20771
c... email:        dbc@welkin.gsfc.nasa.gov
c... phone:        (301) 286-4299
c... 
c... date:         12/03/96
c
c... purpose:
c...  
c...     This subroutine is written for the Global Modeling Initiative of
c...     the Atmospheric Effects of Aviation program, to be used as part
c...     of the heterogeneous chemistry package for the GMI core model.
c...     Modifications have been made to the original code written by
c...     Dr. Ken Carslaw (see documentation below) to fit as a subroutine
c...     into the structure of the GMI model. As few changes as possible
c...     have been made to the original code.  All changes have been documented
c...     below.
c
c... input variable descriptions:
c
c... t - temperature
c... ptot - total pressure (input in mb)
c... qh2o - H2O mixing ratio (ppmv)
c... qhno3 - HNO3 mixing ratio (ppbv)
c... qh2so4 - mixing ratio of h2so4 (ppbv)
c... qhcl - mixing ratio of hcl (ppbv)
c... qhocl - mixing ratio of hocl (ppbv)
c... ==========================================================================================
c...          fortran77 code to calculate composition of aqueous hno3/h2so4/hcl/hocl
c...         stratospheric aerosols (carslaw, luo, peter - geophys. res. lett., 1995)
c...                           carslaw@nike.mpch-mainz.mpg.de
c
c...                           distribution version 1 (2 aug 1995)
c...         
c...         dr ken carslaw
c...         max-planck-institut fuer chemie
c...         postfach 3060
c...         mainz 55020
c...         germany
c...         tel: (49) (0)6131 305 333
c...         fax: (49) (0)6131 305 328
c
c... ==========================================================================================
c... hno3/h2so4 composition based on thermodynamic model of the system
c... hcl/hbr/hno3/h2so4/h2o {carslaw et al, j. phys. chem., 1995 - a thermodynamic model
c... of the system hcl-hno3-h2so4-h2o, including solubilities of hbr, from <200k to 328k. 
c... (also available in fortran77)}
c... hcl solubility parametrisation from luo et al., grl, 1995
c... hocl solubility parametrisation from huthwelker et al., j. at. sci., 1995
c
c... ****************************************************************************
c... the model is valid for 2e-5 mb<pw<2e-3 mb (pw is the water partial pressure)
c... the upper temperature limit is 240k
c... the lower temperature limit is 185k, or tice-3k, which ever is the higher
c... hno3 solubilities are calculated to a maximum temperature of 215k
c... solubilities are calculated on a molality basis (moles of solute per kg of water)
c... ****************************************************************************
c
c... the solubilities of hcl and hocl are assumed not to affect the uptake of hno3 or h2o.
c... this introduces only small errors at the very lowest temperatures, but for a full
c... calculation where interactions between all species in solution are considered, use
c... the model of carslaw et al., given above.
c... the program has been written so that the different steps are clear, and is therefore
c... not coded for optimum calculation speed.  this is left for the user to do if desired.
c... ==========================================================================================

c... dbc - set implicit none and declare all variables
      IMPLICIT NONE
      INTEGER, PARAMETER :: DBL = KIND(0.00D+00)
 
c... dbc - change to real*8 variable definition
      REAL(KIND=DBL) :: ns,msb,mnb,ms,mn,mcl,mhocl
      REAL, SAVE :: qn(10),qs(10),kn(7),ks(7)

c... dbc - declare parameter r real*8
      REAL(KIND=DBL), PARAMETER :: r=8.205e-5

c... dbc - declare logical for error
      LOGICAL :: err

c... dbc - input variable declaration
      REAL :: t,tin,ptot,qhcl,qhocl,qh2so4
      REAL(KIND=DBL) :: qh2o,qhno3,pi

c... dbc - declare currently undeclared internal variables
      REAL :: pw,pn0,phcl0,phocl0,tt,tr,pr,tice,xsb,hsb,xnb
      REAL :: hnb,a,b,c,phi,pn,ws,wn,parthno3,v,w1,w2,hhcl
      REAL :: hhocl,whocl,denssts,density

c... ======================================================================
      data qn/14.5734,0.0615994,-1.14895,0.691693,-0.098863,
     &  0.0051579,0.123472,-0.115574,0.0110113,0.0097914/
      data qs/14.4700,0.0638795,-3.29597,1.778224,-0.223244,
     &  0.0086486,0.536695,-0.335164,0.0265153,0.0157550/
      data kn/-39.136,6358.4,83.29,-17650.0,198.53,
     &  -11948.,-28.469/
      data ks/-21.661,2724.2,51.81,-15732.0,47.004,
     &       -6969.0,-4.6183/
c... ======================================================================

c... dbc - initially set err to false.  subroutine can change this if a problem
c...       occurs
c     err = .false.
c... initially set t to tin
      t = tin

      pw = qh2o*1.e-6*ptot/1013.0

      pn0 = ptot/1013.0*qhno3*1.e-9

c... dbc ns is the number of moles of h2so4 per meter**3
c...   the factor of 100 convertss from millibar to pascal
      ns = ptot*100.0*qh2so4*1.e-9/8.314409/t

c...       ns = qh2so4*1.e6/6.02e23

c      qhcl = 0.0
c      phcl0 = ptot/1013.0*qhcl*1.e-9
c
c      qhocl = 0.0
c      phocl0 = ptot/1013.0*qhocl*1.e-9
c... ========================================================================
cdbc      tt=r*t*ns
cdbc      tr=1.0e4/t-43.4782608
      pr = log(pw)+18.4
      tice = 2668.70/(10.4310-(log(pw)+log(760.0))/log(10.0))
c
c... code amended by d. considine: carslaw's routine only valid for
c... temperatures greater 185 which are also greater than 3 degrees 
c... below frostpoint. so, when temps are less than this, set temps
c... to minimum values.

      if(t.lt.tice-3.0) t = tice-3.0
      if(t.lt.185.0) t = 185.0

      tt = r*t*ns
      tr = 1.0e4/t-43.4782608

      if((t.ge.tice-3.0).and.(t.ge.185.0))then
c... the h2so4/h2o pure solution concentration
        xsb = 1.0/(2.0*(ks(3)+ks(4)/t))*( -ks(1)-ks(2)/t
     &    -sqrt((ks(1)+ks(2)/t)**2-4.0*(ks(3)+ks(4)/t)*(ks(5)+ks(6)/t
     &    +ks(7)*log(t)-log(pw)))) 
        msb = 55.51*xsb/(1.0-xsb)

        if((t.le.215.0).and.(pn0.gt.0.0))then
c... the hno3/h2so4/h2o solution composition
          hsb = qs(1)+qs(2)*tr**2+(qs(3)+qs(4)*tr+qs(5)*tr**2
     &      +qs(6)*tr**3)*pr + (qs(7)+qs(8)*tr+qs(9)*tr**2)
     &      *pr**2+qs(10)*tr*pr**3
          hsb = exp(hsb)
          xnb = 1.0/(2.0*(kn(3)+kn(4)/t))*( -kn(1)-kn(2)/t-
     &      sqrt((kn(1)+kn(2)/t)**2-4.0*(kn(3)+kn(4)/t)*(kn(5)+
     &      kn(6)/t+kn(7)*log(t)-log(pw))))
          mnb = 55.51*xnb/(1.0-xnb)
          hnb = qn(1)+qn(2)*tr**2+(qn(3)+qn(4)*tr+qn(5)*tr**2
     &      +qn(6)*tr**3)*pr + (qn(7)+qn(8)*tr+qn(9)*tr**2)
     &      *pr**2+qn(10)*tr*pr**3
          hnb = exp(hnb)
          a = (tt*hnb*mnb**2 - tt*hsb*mnb*msb - 2.0*mnb**2*msb 
     &      + mnb*msb**2 + hnb*mnb*msb*pn0 - hsb*msb**2*pn0)/
     &      (mnb**2 - mnb*msb)
          b = msb*(-2.0*tt*hnb*mnb+tt*hsb*msb+mnb*msb-
     &      hnb*msb*pn0)/(mnb-msb)
          c = (tt*hnb*mnb*msb**2)/(mnb - msb) 
          phi = atan(sqrt(4.0*(a**2-3.0*b)**3-(-2.0*a**3
     &      +9.0*a*b-27.0*c)**2)/(-2.0*a**3+9.0*a*b
     &      -27.0*c) ) 
          if(phi.lt.0) phi = phi+pi
            ms = -1.0/3.0 
     &      *(a+2.0*sqrt(a**2-3.0*b)*cos((pi+phi)/3.0))
            mn = mnb*(1.0-ms/msb)
            pn = mn/(hnb*mn/(mn+ms)+hsb*ms/(mn+ms))
            ws = ms*0.098076/(1.0+ms*0.098076+mn*0.063012)
c            wn = mn*0.063012/(1.0+ms*0.098076+mn*0.063012)
            parthno3 = (1.0-(pn0-pn)/pn0)
           else
c... assume solution is pure h2so4/h2o
                 ms = msb
                 mn = 0.0
                 ws = msb*0.098076/(1.0+msb*0.098076)
c                 wn = 0.0
c                 pn = pn0
                 parthno3 = 1.0
               endif
           v = ns*98.076/ws/density(ms,mn,t)*1.e6
c... total liquid aerosol volume in um3/cm3
c... +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c... dbc introduce variable to return density

           denssts = density(ms,mn,t)

c... ===============================================================================
c... the solubility of hcl
c... h* (mol/kg/atm) adapted from luo et al., grl, 1995.
c... calculated concentrations assume that hcl is trace component of the aerosol.
c              if(phcl0.gt.0.0)then
c                 w1 = wn
c                 w2 = ws
c                 hhcl = exp(-(21.0+46.610*w1+4.0690*w2
c     &           -4.8370*sqrt(w1)
c     &           +2.1860*sqrt(w2)-63.00*w1**2-40.170*w1*w2
c     &           -1.5710*w2**2)
c     &           -1.0/t*(-7437.0-8327.80*w1+1300.90*w2
c     &           +1087.20*sqrt(w1)
c     &           -242.710*sqrt(w2)+18749.0*w1**2+18500.0*w1*w2
c     &           +5632.0*w2**2)-log(w1+0.610*w2)
c     &           -log(36.461/(1000.0+98.076*ms+63.012*mn)))
c     &           *1.013e3
c                 mcl = (1.0/r/t*phcl0)/(ns/ms + 1.0/r/t/hhcl)
ccc                 wcl = mcl*36.461/(1000.0+98.076*ms+63.012*mn)
ccc                 phcl = mcl/hhcl
ccc                 parthcl = 1.0-(phcl0-phcl)/phcl0
c              endif
c... ===============================================================================
c... the solubility of hocl 
c... h* (mol/kg/atm) from huthwelker et al., j. at. sci, 1995
c... as an approximation, assume h* depends upon total molality
ccc              if(phocl0.gt.0)then
c              if(phocl0.gt.0.) then
c                 hhocl = exp(6.49460-(-0.041070+54.56/t)*(ms+mn)
c     &           -5862.0*(1.0/298.150-1.0/t))
c                 mhocl = (1.0/r/t*phocl0)/(ns/ms + 1.0/r/t/hhocl)
c                 whocl = mhocl*52.46/(1000.0+98.076*ms+63.012*mn)
c... solubility of hocl is low enough to ignore gas phase removal
c              endif
        else

c... dbc delete printouts and include an error flag to determine if any
c...     conditions have been violated
c        err = .true.

        endif
c
        RETURN
        END SUBROUTINE ternary
