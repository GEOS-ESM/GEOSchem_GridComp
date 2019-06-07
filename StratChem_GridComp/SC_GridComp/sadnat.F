      SUBROUTINE sadnat(temp,press,pi,h2ogas,h2ocond,hno3gas,hno3cond,
     &                  sanat,rmed,constantn,satratnat,nnat,densnat,
     &                  signat,rnat,mdens)

!    subroutine:   sadnat

!    written by:   David B. Considine
!                  Earth System Science Interdisciplinary Center and
!                  Department of Meteorology
!                  University of Maryland
!                  College Park, MD 20742

!    email:        dbc@welkin.gsfc.nasa.gov
!    phone:        (301) 405-5282
!    
!    date:         3/5/99
!    
!    purpose:
!        
!       This subroutine is written for the GSFC 3D CTM as part
!       of its polar stratospheric cloud parameterization.  It is
!       a modified version of a subroutine written for the Global Modeling
!       Initiative of the Atmospheric Effects of Aviation
!       Program. 

!    description:  

!       Subroutine sadnat calculates the condensed phase
!       number density of HNO3 and a surface area density
!       consistent with the condensed phase number density.
!       The routine assumes a lognormal size distribution.
!       The user must specify the width of the size distribution
!       in signat.  The user must also specify either the 
!       median radius of the distribution rnat or the total
!       number of nat particles per cubic centimeter, nnat.
!       The logical constantn is a switch to determine whether
!       the surface area density is calculated by assuming a
!       constant number of nat particles and a variable median
!       radius (constantn = true) or by assuming a constant median 
!       radius and a variable number of particles per unit volume.
!       (constantn = false).  The equilibrium vapor pressure of
!       HNO3 over nat is calculated according to Hanson and 
!       Mauersberger, Geophys. Res. Lett., 15,855-858, 1988.
!       If condensed phase H2O is non-zero the code sets the
!       calculated NAT surface area density to zero.  Thus the
!       code can be used to calculate the amount of HNO3 condensing
!       in Type II PSCs, under the assumption that a Type II PSC
!       consists of a cocondensed ice and NAT phase.

!    input variables:

!       temp - local temperature (Kelvin)
!       press - local pressure (millibar)
!       mdens - local density (molecules/cm**3)

! Since the majority of the calculations are in pressure, the following
! must be passed in mole fraction, not number density as indicated.

!       hno3gas - gas phase hno3 concentration (molecules/cm**3)
!       h2ogas - gas phase h2o concentration (molecules/cm**3)
!       hno3cond - condensed phase hno3 concentration (molecules/cm**3)
!       h2ocond - condensed phase h2o concentration (molecules/cm**3)

!    output variables:

!       hno3gas - see above
!       hno3cond - number density of condensed phase hno3 (molecules/cm**3)
!       sanat  - surface area density of condensed nat (cm**(-1))

!    internal variables

!       hno3amb - ambient hno3 concentration (sum of gas and condensed)
!       presstorr - pressure convertsed to Torr (Torr)
!       hno3pp - partial pressure of ambient hno3 concentration (Torr)
!       hno3eq - equilibrium partial pressure of hno3 vapor over nat
!       tsathno3 - saturation temperature (Kelvin)
!       rnatcm - median radius of nat particles in cm

!    parameters

!       constantn - logical to determine strategy for calculating
!                   surface area density.  If true, assumes number
!                   of nat particles remains constant and median
!                   radius changes.  If false, median radius is
!                   assumed constant and number of nat particles adjusts

!       nnat - number of nat particles/volume (particles/cm**3)
!       densnat - mass density of nat (grams/cm**3)
!       masshno3 - mass of a molecule of hno3 (grams)
!       ndensnat - number density of condensed phase nat (molecules/cm**3)
!       signat - 
!       rnat = median radius of nat particles, in microns

!    code starts here:

!    declare all variables
      IMPLICIT NONE
      INTEGER, PARAMETER :: DBL = KIND(0.00D+00)

!    declare logicals
      LOGICAL :: constantn

!    declare real input variables
      REAL :: temp,press,mdens,pi
      REAL :: nnat,densnat,signat,rnat,satratnat
      REAL(KIND=DBL) :: hno3gas,h2ogas,h2ocond

!    delare real output variables
      REAL :: rmed
      REAL(KIND=DBL) :: sanat,hno3cond

!    declare real internal variables
      REAL(KIND=DBL) :: presstorr,rnatcm
      REAL(KIND=DBL) :: b,c,hno3amb,hno3pp,hno3eq,h2opp
      REAL(KIND=DBL) :: loghno3eq,tnuchno3,tsathno3 

!    declare real parameters
      REAL(KIND=DBL) :: logsignatsq
      REAL(KIND=DBL) :: ndensnat,oneThird,twoThirds

!    set internal parameters
      REAL(KIND=DBL), PARAMETER :: ahm = -2.7836
      REAL(KIND=DBL), PARAMETER :: bhm = -0.00088
      REAL(KIND=DBL), PARAMETER :: chm = 38.9855
      REAL(KIND=DBL), PARAMETER :: dhm = -11397.0
      REAL(KIND=DBL), PARAMETER :: ehm = 0.009179
      REAL(KIND=DBL), PARAMETER :: masshno3 = 1.943e-22

!    Combine gas and condensed phase hno3.  This is done because
!    this is an equilibrium routine - given the total number of
!    HNO3 molecules per unit volume, the routine calculates the
!    fraction that should be in NAT and and gas phase.
      hno3amb = hno3gas+hno3cond

!    modification added 10/28/00.  If hno3gas and hno3cond both = 0,
!    code, unphysically large hno3gas results. (problem occured during
!    run with instantaneous hno3cond removal. Don't allow this to
!    happen by including small hno3amb always.
      if(hno3amb.le.0) hno3amb = 1.e-15

!    get partial pressure of hno3 and h2o in Torr
      presstorr = 760./1013.*press
      hno3pp = presstorr*hno3amb
      h2opp = presstorr*h2ogas

!    saturation temperature
      b = (ahm*log10(h2opp)-log10(hno3pp)+chm)/
     &  (ehm+bhm*log10(h2opp))

      c = dhm/(ehm+bhm*log10(h2opp))

      tsathno3 = (-b+sqrt(b*b-4.0*c))/2.0

!    nucleation temperature
!    Need to recalculate b only if satratnat is not 1

      IF(satratnat /= 1.00) THEN
      b = (ahm*log10(h2opp)-log10(hno3pp/satratnat)+chm)/
     &  (ehm+bhm*log10(h2opp))
      END IF

      tnuchno3 = (-b+sqrt(b*b-4.0*c))/2.0

!    if temp is below saturation temperature, calculate condensed hno3

      if(temp.gt.tsathno3) then

!    in this case the temperature is too high for NAT aerosols.

         hno3cond = 0.0
         sanat = 0.0
         hno3gas = hno3amb

       else if(temp.gt.tnuchno3.and.hno3cond.eq.0.0
     &        .and.h2ocond.eq.0.0) then

!    dbc 20000326: added h2ocond to above condition.  If h2ocond
!    doesnt = 0, then ice exists in gridbox and NAT is cocondensing
!    on the ice aerosols.  No nucleation barrier for this process.
  
!    if no HNO3 is condensed already, none can form
         hno3cond = 0.0
         sanat = 0.0
         hno3gas = hno3amb
       else

!    in this case, preexisting NAT aerosols can exist. If no HNO3
!    is condensed already, however, none can form.  Note that in
!    the GSFC 3d CTM we are transporting hno3cond and hno3gas
!    separately.  This means that in the CTM condensed phase HNO3
!    can advect into a region above tnuchno3 and below tsathno3
!    and survive.  In the GMI model (current 3/5/99) hno3gas and
!    hno3cond are transported together.  Because of this, NAT
!    aerosols can only exist at some location if that location
!    has dropped below the nucleation temperature sometime in the
!    past.

         loghno3eq = (ahm+bhm*temp)*log10(h2opp)+chm+dhm/temp+ehm*temp
         hno3eq = 10.0**loghno3eq
         hno3cond = hno3pp-hno3eq

!    converts from torr to mixing ratio

         hno3cond = hno3cond/presstorr
         hno3gas = hno3eq/presstorr

!    is there any condensed phase h2o? If so, don't calculate SAD.
!    If condensed h2o exists then the code assumes that the NAT
!    is cocondensed in a type 2 PSC.
         if(h2ocond.gt.0.0) then
           sanat = 0.0
          else

!    converts hno3 condensed mixing ratio to surface area density
!    use constant NAT particle number density assumption if flag set
!    first converts mixing ratio to number density

            hno3cond = hno3cond*mdens

            logsignatsq = log(signat)*log(signat)
            ndensnat = densnat/masshno3
            rnatcm = rnat*1.e-4
	    
!    represent condensed nitric acid in default REAL

	    oneThird = 1.00E+00/3.00E+00
	    twoThirds = 2.00E+00/3.00E+00

            if(constantn) then
               sanat = (4.0*pi*nnat)**(oneThird)*
     &            (3.0*hno3cond/ndensnat)**(twoThirds)
     &             *exp(-1.0*logsignatsq)
               rmed = (3.0*hno3cond/(ndensnat*4.0*pi*nnat))**(oneThird)
     &            *exp(-3.0/2.0*logsignatsq)
             else
               sanat=3.0/rnatcm*(hno3cond/ndensnat)
     &            *exp(-5.0*logsignatsq/2.0)
             endif

!    converts HNO3 back to mixing ratio

            hno3cond = hno3cond/mdens

          endif
       endif

      RETURN
      END SUBROUTINE sadnat
