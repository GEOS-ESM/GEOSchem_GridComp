      SUBROUTINE sadice(temp,press,pi,h2ogas,h2ocond,saice,
     &                  rmedice,constantn,satratice,nice,
     &                  densice,sigice,rice,mdens)

!    subroutine:   sadice

!    written by:   David B. Considine
!                  Earth System Science Interdisciplinary Center
!                  and Department of Meteorology
!                  3439 Computer and Space Sciences Building
!                  University of Maryland
!                  College Park, MD 20742
!    email:        dbc@welkin.gsfc.nasa.gov
!    phone:        (301) 405-5282
!    
!    date:         3/4/99
!    
!    purpose:
!        
!       This subroutine is written for the Goddard Space Flight
!       Center 3D CTM.  It is a modified version of a subroutine
!       written for the 3D model of the Global Modeling
!       Initiative of the Atmospheric Effects of Aviation
!       Program.  

!    description:  

!       Subroutine sadice calculates the condensed phase
!       number density of H2O and a surface area density
!       consistent with the condensed phase number density.
!       The routine assumes a lognormal size distribution.
!       The user must specify the width of the size distribution
!       in sigice.  The user must also specify either the 
!       median radius of the distribution rice or the total
!       number of ice particles per cubic centimeter, nice.
!       The logical constantn is a switch to determine whether
!       the surface area density is calculated by assuming a
!       constant number of ice particles and a variable median
!       radius (constantn = true) or by assuming a constant median 
!       radius and a variable number of particles per unit volume.
!       (constantn = false).  The equilibrium vapor pressure of
!       H2O over ice is calculated according to Marti and 
!       Mauersberger, Geophys. Res. Lett., 20,363-366, 1993

!    input variables:
!       
!       temp - local temperature (Kelvin)
!       press - local pressure (millibar)
!       h2ogas - gas phase h2o concentration (mixing ratio units)
!       h2ocond - initial condensed phase h2o component (mixing ratio units)

!    output variables:

!       h2ocond - condensed phase h2o mixing ratio
!       saice  - surface area density of condensed ice (cm**(-1))

!    internal variables

!       presspasc - pressure convertsed to Pascal (Pascal)
!       h2opp - partial pressure of ambient h2o concentration (Pascal)
!       h2oeq - equilibrium partial pressure of h2o vapor over ice
!       tsath2o - saturation temperature (Kelvin)
!       ricecm - median radius of ice particles in cm
!       dens - local density (molecules/cm**3)

!    parameters

!       constantn - logical to determine strategy for calculating
!                   surface area density.  If true, assumes number
!                   of ice particles remains constant and median
!                   radius changes.  If false, median radius is
!                   assumed constant and number of ice particles adjusts

!       nice - number of ice particles/volume (particles/cm**3)
!       densice - mass density of ice (grams/cm**3)
!       massh2o - mass of a molecule of h2o (grams)
!       ndensice - number density of condensed phase ice (molecules/cm**3)
!       sigice - 
!       rice = median radius of ice particles, in microns

!    code starts here:

!    declare all variables
      IMPLICIT NONE
      INTEGER, PARAMETER :: DBL = KIND(0.00D+00)

!    declare logicals
      LOGICAL :: constantn

!    declare real input variables
      REAL :: temp,press,mdens,densice,sigice
      REAL :: satratice,nice,rice,pi
      REAL(KIND=DBL) :: h2ogas,h2ocond

!    delare real output variables
      REAL :: rmedice
      REAL(KIND=DBL) :: saice

!    declare real internal variables
      REAL(KIND=DBL) :: presspasc,ricecm
      REAL(KIND=DBL) :: h2opp,h2oeq,tsath2o,tnuch2o,h2oamb
      REAL(KIND=DBL) :: oneThird,twoThirds

!    declare real parameters
      REAL(KIND=DBL) :: logsigicesq,ndensice

      REAL(KIND=DBL), PARAMETER :: amm = -2663.5
      REAL(KIND=DBL), PARAMETER :: bmm = 12.537
      REAL(KIND=DBL), PARAMETER :: massh2o = 2.991e-23

! In GEOS-5 we take H2OCOND from the AGCM as QCTOT.  Find the
! surface area density and effective radius only if H2OCOND exists.
! -----------------------------------------------------------------
      saice = 0.0
      rmedice = 0.0

      HaveIce: IF(h2ocond > 0.00) THEN

! Number density
! --------------
       h2ocond = mdens*h2ocond

! Convert h2o condensed phase number density to surface area density
! ------------------------------------------------------------------
       logsigicesq = log(sigice)*log(sigice)
       ndensice = densice/massh2o
       ricecm = rice*1.e-4

       oneThird=1.00E+00/3.00E+00
       twoThirds=2.00E+00/3.00E+00

       IF(constantn) THEN
        saice = (4.0*pi*nice)**(oneThird)*
     &          (3.0*h2ocond/ndensice)**(twoThirds)*
     &          EXP(-1.0*logsigicesq)
        rmedice = (3.0*h2ocond/(ndensice*4.0*pi*nice))**(oneThird)*
     &            EXP(-3.0/2.0*logsigicesq)
       ELSE
        saice = 3.0/ricecm*(h2ocond/ndensice)*
     &          EXP(-5.0*logsigicesq/2.0)
       END IF

! Convert h2o number density back to mixing ratio
! -----------------------------------------------
       h2ocond = h2ocond/mdens

      END IF HaveIce
 
      RETURN
      END SUBROUTINE sadice
