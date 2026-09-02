      SUBROUTINE calcsulf(tin,ptot,qh2o,sasulf,h2so4gas,rlbs,siglbs,
     &                    nlbs,constantn,mdens,pi)

c... subroutine: calcsulf
c
c... written by: David B. Considine
c...             Earth System Science Interdisciplinary Center and
c...             Department of Meteorology
c...             University of Maryland
c...             College Park, MD 20742
c... email:      dbc@welkin.gsfc.nasa.gov
c... phone:      (301) 405-4282
c... 
c... date:       3/5/99
c
c... purpose:
c...     
c...    This subroutine is written for the GSFC 3D CTM.  It is a modified
c...    version of a subroutine originally written for the Global Modeling
c...    Initiative of the Atmospheric Effects of Aviation
c...    Project.  It is part of the psc parameterization for the CTM.
c
c... description:  
c...    
c...    Subroutine calcsulf takes as input the surface area
c...    density for the background sulfate aerosol distribution.
c...    It calculates the mixing ratio of condensed phase
c...    H2SO4 in ppbv consistent with the input surface
c...    area density.  It does this using an assumed lognormal
c...    size distribution for the sulfate aerosols to converts
c...    from surface area density to volume density.  From the
c...    volume density and the weight percent of H2SO4 calculated
c...    using a binary H2SO4/H2O aerosol composition assumption,
c...    the volume mixing ratio is calculated.

c... set implicit none
      IMPLICIT NONE
      INTEGER, PARAMETER :: DBL = KIND(0.00D+00)

!    Input
      REAL,INTENT(IN)  :: tin,ptot,qh2o,sasulf,rlbs
      REAL,INTENT(IN)  :: siglbs,nlbs,mdens
      LOGICAL,INTENT(IN) :: constantn
      REAL(KIND=DBL), INTENT(IN) :: pi

!    Output
      REAL,INTENT(OUT) :: h2so4gas

!    Local
      REAL(KIND=DBL) :: xsb,t,pw,msb,denss,wsb,tice
      REAL(KIND=DBL) :: logsiglbssq,rlbscm,v,dens

c... set up coefficient array
      REAL, SAVE :: ks(7) = (/ -21.661,  2724.2, 51.81, -15732.0,
     &                          47.004, -6969.0, -4.6183/)

c... mass h2so4 in amu
      REAL, PARAMETER :: h2so4amu=98.08
c... grams per amu
      REAL, PARAMETER :: gpamu=1.66e-24

c... converts from input surface area density to volume density
c... use constant liquid binary sulfate particle number density 
c... assumption if flag set
      logsiglbssq = log(siglbs)*log(siglbs)
      rlbscm = rlbs*1.e-4

      if(constantn)then
         v = sasulf**1.5/3.0/(sqrt(4.0*pi*nlbs))*
     &     exp(3.0/2.0*logsiglbssq)

ccc         rmed=(3.0*v/(4.0*pi*nlbs))**(1.0/3.0)
ccc            *exp(-3.0/2.0*logsiglbssq)
       else
         v = rlbscm*sasulf/3.0*
     &      *exp(5.0*logsiglbssq/2.0)
       endif

c... The code convertss from volume density (cm**3 aerosol/cm**3 air)
c... to mixing ratio of condensed phase H2SO4. First step: calculate
c... the density of a binary H2SO4/H2O solution using Carslaw's method.

c... first calculate the mole fraction of H2SO4/H2O solution - this
c... is taken straight from Carslaw's ternary code

      pw = qh2o*ptot/1013.0
      tice = 2668.70/(10.4310-(log(pw)+log(760.0))/log(10.0))

      t = tin
      if(t.lt.tice-3.0) t = tice-3.0
      if(t.lt.185.0) t = 185.0

      xsb = 1.0/(2.0*(ks(3)+ks(4)/t))*( -ks(1)-ks(2)/t-
     &  sqrt((ks(1)+ks(2)/t)**2-4.0*(ks(3)+ks(4)/t)*(ks(5)+ks(6)/t
     &  +ks(7)*log(t)-log(pw)))) 

c...  get molality (moles/mass) of binary solution
      msb = 55.51*xsb/(1.0-xsb)

c... now calculate density - again this is due to carslaw
      denss = (1000.0+123.64*msb-5.6e-4*msb*t**2 
     &       -29.54*msb**1.5 + 1.814e-4*msb**1.5*t**2
     &       + 2.343*msb**2  -1.487e-3*msb**2*t 
     &       -1.324e-5*msb**2*t**2)*1.e-3

c... now calculate weight fraction of H2SO4 in this binary solution:
      wsb = msb*0.098076/(1.0+msb*0.098076)

c... now I can calculate the concentration of H2SO4 in molecules/cm3 air
      h2so4gas = v*denss*wsb/h2so4amu/gpamu

c... now converts to mixing ratio
      h2so4gas = h2so4gas/mdens

      RETURN
      END SUBROUTINE calcsulf
