
!=============================================================================
!
! $Id$
!
! CODE DEVELOPERS
!   Condense originally written by:
!     David B. Considine
!     Code 916
!     NASA Goddard Space Flight Center
!     Greenbelt, MD  20771
!     dbc@welkin.gsfc.nasa.gov
!     301-286-4299
!   LLNL modifications:
!     John Tannahill, LLNL
!     jrt@llnl.gov
!
! FILE
!   condense.F
!
!   Aerosol surface area density and condensed phase mixing ratio module
!   for the global modeling initiative three dimensional assessment model.
!
!   Version:  1.1
!   Date:     4/7/97
!
!   SEE condense.doc FOR ADDITIONAL DOCUMENTATION.
!
! ROUTINES
!   Condense
!   Sad_Ice*
!   Sad_NAT*
!   Sad_STS*
!   Calcsulf_DBC*
!   Ternary_DBC*
!   Density_DBC*
!   Sediment_DBC*
!   Fall_Vel*
!
! *These names changed from the originals so GMICHEM and STRATCHEM can be
!  compiled without name duplication in GEOS-5.  Eric Nielsen 7Nov2007
!
!=============================================================================


!-----------------------------------------------------------------------------
!
! ROUTINE
!   Condense
!
! DESCRIPTION
!   This subroutine is written for the Global Modeling Initiative of the
!   Atmospheric Effects of Aviation Program.  It is part of the heterogeneous
!   chemistry package to be used in the GMI core model.
!
!   Update_Condphase calculates surface area densities for ICE, NAT and STS
!   (supercooled ternary sulfate) aerosols given appropriate input values from
!   the GMI core model.  It is set up to loop over the spatial indices of the
!   three dimensional arrays passed to the subroutine as arguments.  Only the
!   subset of the arrays specified by the dimensions i1, i2, ju1, j2, k1, &
!   k2 is considered.  The subroutine calls the subroutines Sadice, Sadnat,
!   and Sadsts, which calculate the ice, nat, and sts surface area densities.
!   The various parameters for these subroutines which control the assumed
!   particle size distributions for each of the aerosol types are set in the
!   parameter section of Update_Condphase, and then are passed to the other
!   subroutines.
!
!   HNO3 and H2O arrays:  It is extremely important to understand the various
!   arrays for HNO3 and H2O to use this subroutine correctly.  There are two
!   arrays for HNO3:  hno3gas and hno3cond.  The code assumes that both of
!   these arrays are transported species and will adjust the values in both of
!   the arrays to account for sequestration of HNO3 into condensed phase from
!   gas phase and sedimentation of HNO3.  The hno3cond array stores the amount
!   of HNO3 that is in condensed phase at a gridpoint but does not distinguish
!   whether the condensed phase HNO3 is in the form of STS, NAT, or Ice
!   aerosol.  There are three arrays for H2O:  h2oback, dehyd, and h2ocond.
!   The code assumes that h2oback is a fixed H2O distribution from
!   observations and does not change its value.  The h2ocond array should not
!   be transported, and records the mixing ratio of condensed phase H2O that
!   exists at any point in the model.  The dehyd array should be transported
!   and is used by the code to account for dehydration of parcels due to
!   sedimentation.  In other words, the total number of H2O molecules at a
!   particular location can be found by subtracting dehyd from h2oback.
!   Because negative mixing ratios typically blow up transport algorithms,
!   dehyd is a positive array like the others and cannot account for
!   evaporation and moistening of layers below a sedimentation region.  If
!   aircraft H2O is being transported, it needs to be added to h2oback before
!   calling the Update_Condphase routine.
!
! ARGUMENTS
!   INPUT:
!     dt         : time step (s)
!     tropp      : Tropopause pressure (hPa)
!     pres3c     : pressure (mb)
!     temp3      : temperature (degK)
!     dehyd      : transported dehydration (mixing ratio)
!     dz         : gridbox height (cm)
!     h2oback    : background h2o array (mixing ratio)
!     hno3cond   : condensed phase hno3 array (mixing ratio)
!     hno3gas    : gas phase hno3 array (mixing ratio)
!     lbssad     : liquid binary sulfate background surface area density
!                  (cm^-1)
!     dehyd_opt  : 0=Disabled. 1=Enabled.
!
!   OUTPUT:
!     denssts    : density of sts aerosol calculated by Ternary routine
!                  (g/cm^3)
!     h2ocond    : condensed phase h2o array (mixing ratio)
!     h2so4gas   : h2so4 inferred from lbssad; actually a misnomer,
!                  since it is assumed that all inferred h2so4 is in
!                  condensed phase (mixing ratio)
!     icesad     : surface area density of ICE aerosols (cm^-1)
!     natsad     : surface area density of NAT aerosols (cm^-1)
!     stssad     : surface area density of STS aerosols (cm^-1)
!     reffice    : effective radius of ICE aerosols  (cm)
!     reffnat    : effective radius of NAT aerosols  (cm)
!     reffsts    : effective radius of STS aerosols  (cm)
!     vfall      : effective aerosol fall velocities (cm/s)
!
! PARAMETERS CONTROLLING SUBROUTINE CALCULATIONS
!   CALCSTS      : if true, calculate STS and not NAT, otherwise NAT only
!   DOSEDIMENT   : if true, do denitrification step
!   PMAX         : maximum pressure
!   PMIN         : minimum pressure
!   TMAX         : maximum local temperature -> only calculate if below tmax
!
!   CONSTANTNICE : if true, use NICE and SIGICE; if false, use RICE and SIGICE
!   DENSICE      : mass density of ICE (g/cm^3)
!   NICE         : number of ICE particles (cm^-3)
!   RICE         : median radius of ICE particle size distribution if used
!                  (microns)
!   SIGICE       : width of lognormal particle size distribution
!   SATRATICE    : saturation ratio for ICE aerosols
!
!   CONSTANTNNAT : if true, use NNAT and SIGNAT to determine SAD DENSNAT
!                  density of NAT aerosols (g/cm^3)
!   NNAT         : number of NAT particles (cm^-3)
!   RNAT         : if  CONSTANTNNAT is false, RNAT is median radius of
!                  NAT particle size distribution
!   SIGNAT       : width of lognormal particle size distribution for NAT
!                  aerosols
!   SATRATNAT    : saturation ratio for NAT aerosols
!
!   CONSTANTNSTS : if true, use NSTS instead of RSTS
!   NSTS         : number of STS particles (cm^-3)
!   RSTS         : median radias of STS particle size distribution (microns)
!   SIGSTS       : controls width of STS particle size distribution
!
!-----------------------------------------------------------------------------

      SUBROUTINE Condense  &
     &  (dt, tropp, pres3c, temp3, dehyd, dehyd_opt, dz,  &
     &   h2oback, hno3cond, hno3gas, lbssad,  &
     &   denssts, h2ocond, h2so4gas, icesad, natsad, stssad,  &
     &   reffice, reffnat, reffsts, vfall, &
     &   pr_diag, loc_proc, londeg, latdeg, NoPSCZone, PMAX, &
     &   ilo, ihi, julo, jhi, i1, i2, ju1, j2, k1, k2)

      implicit none

!     ----------------------
!     Argument declarations.
!     ----------------------

      logical, intent(in) :: pr_diag
      integer, intent(in) :: loc_proc
      integer, intent(in) :: ilo, ihi, julo, jhi
      integer, intent(in) :: i1, i2, ju1, j2, k1, k2
      integer, intent(in) :: dehyd_opt
      real*8,  intent(in) :: tropp(i1:i2, ju1:j2)

      integer, intent(in) :: NoPSCZone
      integer, intent(in) :: PMAX         ! hPa, replaces PARAMETER PMAX
      real*8,  intent(in) :: londeg(i1:i2, ju1:j2), latdeg(i1:i2, ju1:j2)

      real*8  :: dt
      real*8  :: pres3c  (ilo:ihi, julo:jhi, k1:k2)
      real*8  :: temp3   (ilo:ihi, julo:jhi, k1:k2)
      real*8  :: dehyd   (i1:i2,   ju1:j2,   k1:k2)
      real*8  :: dz      (i1:i2,   ju1:j2,   k1:k2)
      real*8  :: h2oback (i1:i2,   ju1:j2,   k1:k2)
      real*8  :: hno3cond(i1:i2,   ju1:j2,   k1:k2)
      real*8  :: hno3gas (i1:i2,   ju1:j2,   k1:k2)
      real*8  :: lbssad  (i1:i2,   ju1:j2,   k1:k2)
      real*8  :: denssts (i1:i2,   ju1:j2,   k1:k2)
      real*8  :: h2ocond (i1:i2,   ju1:j2,   k1:k2)
      real*8  :: h2so4gas(i1:i2,   ju1:j2,   k1:k2)
      real*8  :: icesad  (i1:i2,   ju1:j2,   k1:k2)
      real*8  :: natsad  (i1:i2,   ju1:j2,   k1:k2)
      real*8  :: stssad  (i1:i2,   ju1:j2,   k1:k2)
      real*8  :: reffice (i1:i2,   ju1:j2,   k1:k2)
      real*8  :: reffnat (i1:i2,   ju1:j2,   k1:k2)
      real*8  :: reffsts (i1:i2,   ju1:j2,   k1:k2)
      real*8  :: vfall   (i1:i2,   ju1:j2,   k1:k2)

!     -----------------------
!     Parameter declarations.
!     -----------------------

!     ----------------------------------------
!     Parameters controlling overall behavior:
!     ----------------------------------------

      logical, parameter :: CALCSTS    = .false.
      logical, parameter :: DOSEDIMENT = .true.

!     ------------------------------------------------------
!     If PSCMaxP <= 0, the tropopuase pressure is the upper
!     limit to the range of application. Otherwise, the
!     lesser of PSCMaxP and the tropopause pressure is used.
!     ------------------------------------------------------

      real*8,  parameter :: PMIN =  10.0d0
      real*8,  parameter :: TMAX = 240.0d0

!     ------------------------------------
!     Parameters controlling ice aerosols:
!     ------------------------------------

      logical, parameter :: CONSTANTNICE = .true.

      real*8,  parameter :: DENSICE   =  1.0d0
      real*8,  parameter :: NICE      =  1.0d-2
      real*8,  parameter :: RICE      = 10.0d0
      real*8,  parameter :: SATRATICE =  1.0d0
      real*8,  parameter :: SIGICE    =  1.6d0

!     ------------------------------------
!     Parameters controlling nat aerosols:
!     ------------------------------------

      logical, parameter :: CONSTANTNNAT = .true.

      real*8,  parameter :: DENSNAT   = 1.6d0
      real*8,  parameter :: NNAT      = 0.1d0
      real*8,  parameter :: RNAT      = 0.4d0
      real*8,  parameter :: SIGNAT    = 1.6d0
      real*8,  parameter :: SATRATNAT = 1.0d0

!     ------------------------------------
!     Parameters controlling sts aerosols:
!     ------------------------------------

      logical, parameter :: CONSTANTNSTS = .true.

      real*8,  parameter :: NSTS   = 10.0d0
      real*8,  parameter :: RSTS   =  4.0d0
      real*8,  parameter :: SIGSTS =  1.6d0

!     ------------------------------------
!     Parameters controlling lbs aerosols:
!     ------------------------------------

      logical, parameter :: CONSTANTNLBS = .true.

      real*8,  parameter :: NLBS   = 10.0d0
      real*8,  parameter :: RLBS   =  0.1d0
      real*8,  parameter :: SIGLBS =  1.6d0

!     ----------------------
!     Variable declarations.
!     ----------------------

      integer :: il, ij, ik

      real*8  :: h2ogas
      real*8  :: pscPMAX
      real*8  :: templat

!     ----------------
!     Begin execution.
!     ----------------

!     ==================
      do ik = k1, k2
        do ij = ju1, j2
          do il = i1, i2
!     ==================

!           ------------------------------------------------------
!           Range check.  Use tropopause pressure if PSCMaxP <= 0.
!              Otherwise, choose the MIN of PSCMaxP and tropp.
!           ------------------------------------------------------
            IF(PMAX <= 0) THEN
             pscPMAX = tropp(il,ij)
            ELSE
             pscPMAX = MIN(PMAX*1.0,tropp(il,ij))
            END IF

!           ---------------------------------------------------------
!           Some time can be saved if calculations are not done when
!           the temperature is too warm for condensed phases to form.
!           Also, we can arbitrarily limit the pressures to those in
!           the stratosphere where condensation is likely (i.e., the
!           lower stratosphere).
!           ---------------------------------------------------------

            if ((temp3 (il,ij,ik) <= TMAX)    .and.  &
     &          (pres3c(il,ij,ik) <= pscPMAX) .and.  &
     &          (pres3c(il,ij,ik) >= PMIN)    .and. ABS(latdeg(il,ij)) >= NoPSCZone) then

!              ===========
               call Sad_Ice  &
!              ===========
     &           (temp3(il,ij,ik), pres3c(il,ij,ik),  &
     &            h2oback(il,ij,ik), dehyd(il,ij,ik), dehyd_opt, &
     &            h2ocond(il,ij,ik), icesad(il,ij,ik),  &
     &            reffice(il,ij,ik), CONSTANTNICE, SATRATICE,  &
     &            NICE, DENSICE, SIGICE, RICE)

              if (CALCSTS .and.  &
     &            (h2ocond(il,ij,ik) == 0.0d0)) then

	        h2ogas = h2oback(il,ij,ik)
	        IF(dehyd_opt == 1) h2ogas = h2ogas - dehyd(il,ij,ik)

!               ===========
                call Sad_STS  &
!               ===========
     &            (temp3(il,ij,ik), pres3c(il,ij,ik),  &
     &             h2ogas,  &
     &             hno3gas(il,ij,ik),  &
     &             hno3cond(il,ij,ik), stssad(il,ij,ik),  &
     &             reffsts(il,ij,ik), denssts(il,ij,ik),  &
     &             CONSTANTNSTS, NSTS, SIGSTS, RSTS,  &
     &             lbssad(il,ij,ik), h2so4gas(il,ij,ik),  &
     &             RLBS, SIGLBS, NLBS, CONSTANTNLBS,  &
     &             h2oback(il,ij,ik))

              else

                stssad(il,ij,ik) = 0.d0

!               -----------------------------------------------------------
!               In the nat/ice calculation (CALCSTS = false), this code
!               is set up so that sadnat will calculate the amount of HNO3
!               removed from gas phase even if ice has formed at this
!               location.  Subroutine sets NAT surface area density to zero
!               if ice exists at gridpoint.
!               -----------------------------------------------------------

	        h2ogas = h2oback(il,ij,ik)
	        IF(dehyd_opt == 1) h2ogas = h2ogas - dehyd(il,ij,ik)

!               ===========
                call Sad_NAT  &
!               ===========
     &            (temp3(il,ij,ik), pres3c(il,ij,ik), h2ogas,  &
     &             h2ocond(il,ij,ik), hno3gas(il,ij,ik),  &
     &             hno3cond(il,ij,ik), natsad(il,ij,ik),  &
     &             reffnat(il,ij,ik), CONSTANTNNAT, SATRATNAT, NNAT,  &
     &             DENSNAT, SIGNAT, RNAT)

              end if

            else

!             ---------------------------------------------------------
!             Temperature and pressure ranges are not correct.  In this
!             case, we want to make sure that any values which on a
!             previous step might have been assigned a value and which
!             might be used elsewhere are now set to 0.
!             ---------------------------------------------------------

              hno3gas(il,ij,ik)  =  &
     &          hno3gas(il,ij,ik) + hno3cond(il,ij,ik)
              hno3cond(il,ij,ik) = 0.0d0

              dehyd(il,ij,ik)    = 0.0d0
              h2ocond(il,ij,ik)  = 0.0d0

              icesad(il,ij,ik)   = 0.0d0
              natsad(il,ij,ik)   = 0.0d0
              stssad(il,ij,ik)   = 0.0d0

              reffice(il,ij,ik)  = 0.0d0
              reffnat(il,ij,ik)  = 0.0d0
              reffsts(il,ij,ik)  = 0.0d0

            end if

!     ==========
          end do
        end do
      end do
!     ==========


      if (DOSEDIMENT) then
!       =============
        call Sediment_DBC  &
!       =============
     &    (dz, dt, temp3, pres3c,  &
     &     hno3cond, h2ocond, dehyd, dehyd_opt, &
     &     RSTS, RNAT, RICE, reffsts, reffnat, reffice, SIGSTS,  &
     &     SIGNAT, SIGICE, CONSTANTNSTS, CONSTANTNNAT, CONSTANTNICE,  &
     &     DENSNAT, DENSICE, vfall, CALCSTS, denssts, &
     &     ilo, ihi, julo, jhi, i1, i2, ju1, j2, k1, k2)
      end if

      return
      END SUBROUTINE Condense


!     =================

      SUBROUTINE Sad_Ice(temp3,pres3c,h2oback,dehyd,dehyd_opt,h2ocond,icesad,  &
     &           reffice,constantn,satratice,nice,densice,sigice,rice)

! subroutine:   sadice
!
! written by:   David B. Considine
!               Code 916
!               NASA Goddard Space Flight Center
!               Greenbelt, MD 20771
! email:        dbc@welkin.gsfc.nasa.gov
! phone:        (301) 286-4299
!
! date:         11/25/96
!
! purpose:
!
!    This subroutine is written for the Global Modeling
!    Initiative of the Atmospheric Effects of Aviation
!    Program.  It is part of the heterogeneous chemistry
!    package to be used in the GMI core model.
!
! description:
!
!    Subroutine sadice calculates the condensed phase
!    number density of H2O and a surface area density
!    consistent with the condensed phase number density.
!    The routine assumes a lognormal size distribution.
!    The user must specify the width of the size distribution
!    in sigice.  The user must also specify either the
!    median radius of the distribution rice or the total
!    number of ice particles per cubic centimeter, nice.
!    The logical constantn is a switch to determine whether
!    the surface area density is calculated by assuming a
!    constant number of ice particles and a variable median
!    radius (constantn = true) or by assuming a constant median
!    radius and a variable number of particles per unit volume.
!    (constantn = false).  The equilibrium vapor pressure of
!    H2O over ice is calculated according to Marti and
!    Mauersberger, Geophys. Res. Lett., 20,363-366, 1993
!
! input variables:
!
!    temp3 - local temperature (degK)
!    pres3c - local pressure (mb)
!    dens - local density (molecules/cm^3)
!    h2oback - fixed background h2o concentration (mixing ratio)
!    dehyd - dehydrated h2o component (mixing ratio)
!    dehyd_opt - Use dehyd? 0=disabled, 1=enabled.
!
! output variables:
!
!    h2ocond - condensed phase h2o mixing ratio
!    icesad  - surface area density of condensed ice (cm^-1)
!
! internal variables:
!
!    presspasc - pressure converted to Pascal (Pascal)
!    h2opp - partial pressure of ambient h2o concentration (Pascal)
!    h2oeq - equilibrium partial pressure of h2o vapor over ice
!    tsath2o - saturation temperature (degK)
!    ricecm - median radius of ice particles in cm
!
! parameters
!
!    constantn - logical to determine strategy for calculating
!                surface area density.  If true, assumes number
!                of ice particles remains constant and median
!                radius changes.  If false, median radius is
!                assumed constant and number of ice particles adjusts
!
!    nice - number of ice particles/volume (particles/cm^3)
!    densice - mass density of ice (g/cm^3)
!    massh2o - mass of a molecule of h2o (g)
!    ndensice - number density of condensed phase ice (molecules/cm^3)
!    sigice -
!    rice = median radius of ice particles, in microns
!
! code starts here:

! declare all variables

      implicit none

#     include "gmi_phys_constants.h"

! declare logicals

      logical constantn

! declare integers

      INTEGER, INTENT(IN) :: dehyd_opt

! declare real input variables

      real*8 temp3,pres3c,dens,h2oback,dehyd,h2ocond

! delare real output variables

      real*8 icesad,reffice

! declare real internal variables

      real*8 presspasc,h2opp,tsath2o,tnuch2o,h2oeq,ricecm,h2oamb

! declare real constants

      real*8 amm,bmm,pi,massh2o

! declare real parameters

      real*8 nice,densice,ndensice,sigice,logsigicesq,rice,  &
     &       satratice

      parameter(amm = -2663.5d0)
      parameter(bmm = 12.537d0)
      parameter(pi = GMI_PI)
      parameter(massh2o = 2.991d-23)

! In GEOS-5 we take H2OCOND from the AGCM as QCTOT.  Find the
! surface area density and effective radius only if H2OCOND exists.
! -----------------------------------------------------------------
      icesad = 0.00D+00
      reffice = 0.00D+00
      
      HaveIce: IF(h2ocond > 0.00D+00) THEN

! Air number density
! ------------------
       presspasc=100.00D+00*pres3c
       dens = AVOGAD*presspasc/GAS_CONST_J/temp3/1.00D+06

! Convert to number density (molecules/cm^3 air)
! ---------------------------------------------
       h2ocond = h2ocond*dens

! Convert h2o condensed phase number density to surface area density
! ------------------------------------------------------------------
       logsigicesq = LOG(sigice)*LOG(sigice)
       ndensice = densice/massh2o
       ricecm = rice*1.00D-04

       IF(constantn) THEN

         icesad = (4.00D+00*pi*nice)**(1.00D+00/3.00D+00)* &
	          (3.00D+00*h2ocond/ndensice)**(2.00D+00/3.00D+00)*  &
                  EXP(-1.00D+00*logsigicesq)

         reffice = (3.00D+00*h2ocond/(ndensice*4.00D+00*pi*nice))**(1.00D+00/3.00D+00)*  &
                   EXP(-3.00D+00/2.00D+00*logsigicesq)

       ELSE

         icesad = 3.00D+00/ricecm*(h2ocond/ndensice)*  &
                  EXP(-5.00D+00*logsigicesq/2.00D+00)

       END IF

! Convert h2ocond number density back to mixing ratio
! ---------------------------------------------------
       h2ocond = h2ocond/dens

      END IF HaveIce
 
      RETURN
      END SUBROUTINE Sad_Ice


!     =================

      SUBROUTINE Sad_NAT(temp3,pres3c,h2ogas,h2ocond,hno3gas,hno3cond,  &
     &           natsad,reff,constantn,satratnat,nnat,densnat,signat,  &
     &           rnat)

! subroutine:   sadnat
!
! written by:   David B. Considine
!               Code 916
!               NASA Goddard Space Flight Center
!               Greenbelt, MD 20771
! email:        dbc@welkin.gsfc.nasa.gov
! phone:        (301) 286-4299
!
! date:         11/25/96
!
! purpose:
!
!    This subroutine is written for the Global Modeling
!    Initiative of the Atmospheric Effects of Aviation
!    Program.  It is part of the heterogeneous chemistry
!    package to be used in the GMI core model.
!
! description:
!
!    Subroutine sadnat calculates the condensed phase
!    number density of HNO3 and a surface area density
!    consistent with the condensed phase number density.
!    The routine assumes a lognormal size distribution.
!    The user must specify the width of the size distribution
!    in signat.  The user must also specify either the
!    median radius of the distribution rnat or the total
!    number of nat particles per cubic centimeter, nnat.
!    The logical constantn is a switch to determine whether
!    the surface area density is calculated by assuming a
!    constant number of nat particles and a variable median
!    radius (constantn = true) or by assuming a constant median
!    radius and a variable number of particles per unit volume.
!    (constantn = false).  The equilibrium vapor pressure of
!    HNO3 over nat is calculated according to Hanson and
!    Mauersberger, Geophys. Res. Lett., 15,855-858, 1988.
!
! input variables:
!
!    temp3 - local temperature (degK)
!    pres3c - local pressure (mb)
!    dens - local density (molecules/cm^3)
!    hno3gas - gas phase hno3 concentration (molecules/cm^3)
!    h2ogas - gas phase h2o concentration (molecules/cm^3)
!    hno3cond - condensed phase hno3 concentration (molecules/cm^3)
!    h2ocond - condensed phase h2o concentration (molecules/cm^3)
!
! output variables:
!
!    hno3gas - see above
!    hno3cond - number density of condensed phase hno3 (molecules/cm^3)
!    natsad  - surface area density of condensed nat (cm^-1)
!
! internal variables:
!
!    hno3amb - ambient hno3 concentration (sum of gas and condensed)
!    presstorr - pressure converted to Torr (Torr)
!    hno3pp - partial pressure of ambient hno3 concentration (Torr)
!    hno3eq - equilibrium partial pressure of hno3 vapor over nat
!    tsathno3 - saturation temperature (degK)
!    rnatcm - median radius of nat particles in cm
!
! parameters
!
!    constantn - logical to determine strategy for calculating
!                surface area density.  If true, assumes number
!                of nat particles remains constant and median
!                radius changes.  If false, median radius is
!                assumed constant and number of nat particles adjusts
!
!    nnat - number of nat particles/volume (particles/cm^3)
!    densnat - mass density of nat (g/cm^3)
!    masshno3 - mass of a molecule of hno3 (g)
!    ndensnat - number density of condensed phase nat (molecules/cm^3)
!    signat -
!    rnat = median radius of nat particles, in microns
!
! code starts here:

! declare all variables

      implicit none

#     include "gmi_phys_constants.h"

! declare logicals

      logical constantn

! declare real input variables

      real*8 temp3,pres3c,dens,hno3gas,h2ogas,h2ocond

! delare real output variables

      real*8 hno3cond,natsad,reff

! declare real internal variables

      real*8 presstorr,hno3pp,tsathno3,tnuchno3,hno3eq,rnatcm,  &
     &       hno3amb,h2opp,b,c,loghno3eq

! declare real constants

      real*8 ahm,bhm,chm,dhm,ehm,pi,masshno3

! declare real parameters

      real*8 nnat,densnat,ndensnat,signat,logsignatsq,rnat,  &
     &       satratnat

! set internal parameters

      parameter(ahm = -2.7836d0)
      parameter(bhm = -0.00088d0)
      parameter(chm = 38.9855d0)
      parameter(dhm = -11397.0d0)
      parameter(ehm = 0.009179d0)
      parameter(pi = GMI_PI)
      parameter(masshno3=1.943d-22)

! Combine gas and condensed phase hno3.  This is done because
! this is an equilibrium routine - given the total number of
! HNO3 molecules per unit volume, the routine calculates the
! fraction that should be in NAT and and gas phase.

      hno3amb=hno3gas+hno3cond

! Modification added 11/21/01.  If hno3gas and hno3cond both = 0,
! unphysically large hno3gas results (in GSFC 3D CTM with this
! PSC param installed). (Problem occured during run with
! instantaneous hno3cond removal.) Don't allow this to happen by
! including small hno3amb always.

      if(hno3amb.le.0.0d0)hno3amb=1.0d-15

! calculate local density using ideal gas law

      dens=AVOGAD*pres3c*100.d0/GAS_CONST_J/temp3/1.d6

! get partial pressure of hno3 and h2o in Torr

      presstorr=760.d0/1013.d0*pres3c
      hno3pp=presstorr*hno3amb

      h2opp=presstorr * (Max (h2ogas, 0.5d-06))

! saturation temperature

      b=(ahm*log10(h2opp)-log10(hno3pp)+chm)/  &
     &  (ehm+bhm*log10(h2opp))

      c=dhm/(ehm+bhm*log10(h2opp))

      tsathno3=(-b+sqrt(b*b-4.0d0*c))/2.d0

! nucleation temperature

      b=(ahm*log10(h2opp)-log10(hno3pp/satratnat)+chm)/  &
     &  (ehm+bhm*log10(h2opp))

      tnuchno3=(-b+sqrt(b*b-4.0d0*c))/2.d0

! if temp3 is below saturation temperature, calculate condensed hno3

      if(temp3.gt.tsathno3)then

         hno3cond=0.d0
         natsad=0.d0
         hno3gas=hno3amb

      else if(temp3.gt.tnuchno3.and.hno3cond.eq.0.0d0.and.  &
     &        h2ocond.eq.0.0d0)then

! Change:  11/21/01: added h2ocond to above condition.  If h2ocond
! doesnt = 0, then ice exists in gridbox and NAT is cocondensing
! on the ice aerosols.  No nucleation barrier for this process.

         hno3cond=0.d0
         natsad=0.d0
         hno3gas=hno3amb

      else

         loghno3eq=(ahm+bhm*temp3)*log10(h2opp)+chm+dhm/temp3+ehm*temp3
         hno3eq=10.d0**loghno3eq
         hno3cond=hno3pp-hno3eq

! convert from torr to mixing ratio

         hno3cond=hno3cond/presstorr
         hno3gas=hno3eq/presstorr

! is there any condensed phase h2o? If so, don't calculate SAD:

         if(h2ocond.gt.0.d0)then

            natsad=0.d0

         else

! convert hno3 condensed mixing ratio to surface area density
! use constant NAT particle number density assumption if flag set
! first convert mixing ratio to number density

         hno3cond=hno3cond*dens


         logsignatsq=log(signat)*log(signat)
         ndensnat=densnat/masshno3
         rnatcm=rnat*1.d-4


            if(constantn)then

               natsad=(4.d0*pi*nnat)**(1.d0/3.d0)*  &
     &            (3.d0*hno3cond/ndensnat)**(2.d0/3.d0)  &
     &             *exp(-1.d0*logsignatsq)

               reff=(3.d0*hno3cond/(ndensnat*4.d0*pi*nnat))**(1.d0/3.d0)  &
     &            *exp(-3.d0/2.d0*logsignatsq)

            else

               natsad=3.d0/rnatcm*(hno3cond/ndensnat)  &
     &            *exp(-5.d0*logsignatsq/2.d0)

            endif

! convert HNO3 back to mixing ratio

         hno3cond=hno3cond/dens

         endif

      endif


      return
      END SUBROUTINE Sad_NAT


!     =================

      SUBROUTINE Sad_STS  &
     &          (temp3,pres3c,h2ogas,hno3gas,hno3cond,  &
     &           stssad,reff,denssts,constantn,nsts,sigsts,rsts,  &
     &           lbssad,h2so4gas,rlbs,siglbs,nlbs,constantnlbs,h2oback)

! subroutine:   sadsts
!
! written by: David B. Considine
!               Code 916
!               NASA Goddard Space Flight Center
!               Greenbelt, MD 20771
! email:        dbc@welkin.gsfc.nasa.gov
! phone:        (301) 286-4299
!
! date:         12/03/96
!
! purpose:
!
!    This subroutine is written for the Global Modeling
!    Initiative of the Atmospheric Effects of Aviation
!    Program.  It is part of the heterogeneous chemistry
!    package to be used in the GMI core model.
!
! description:
!
!    Subroutine sadsts calculates the surface area density
!    of supercooled ternary sulfate aerosols when called.
!    It uses the equilibrium model of Ken Carslaw to calculate
!    the total volume density of the STS aerosol and
!    converts the volume density to surface area density
!    based on size distribution parameters passed to the
!    code from the calling subroutine.  The code takes as
!    input an observed sulfate aerosol surface area density
!    from observations.  It uses a sulfate aerosol size
!    distribution assumption to convert from surface area
!    density to volume density.  It then takes this volume density
!    field and uses binary H2SO4/H2O assumptions to calculate
!    the condensed phase H2SO4 amount consistent with the
!    sulfate aerosol surface area density.  This H2SO4 amount
!    is then used along with HNO3 and H2O to calculate the
!    surface area consistent with a ternary composition assumption.
!    NOTE: the field stssad is the DIFFERENCE between the calculated
!    surface area density assuming sts composition and the background
!    sulfate surface area density, when the difference produces a
!    positive result.  This allows the user of the code to simply
!    sum the background and sts surface area densities to get a
!    valid total.

! variable declarations

      implicit none

#     include "gmi_phys_constants.h"

! declare input variables

      logical constantn,constantnlbs

      real*8 temp3,pres3c,hno3gas,h2ogas,  &
     & h2so4gas,lbssad,hno3cond,nsts,sigsts,rsts,  &
     & rlbs,siglbs,nlbs,h2oback

! declare output variables

      real*8 stssad,reff

! declare internal variables

      logical saderr

      real*8 h2oppmv,hno3ppbv,h2so4ppbv,  &
     &       ws,sadwn,v,parthno3,hno3amb,  &
     &       logsigstssq,rstscm,pi,denssts

! set internal parameters

      parameter(pi = GMI_PI)

! set up call to Carslaw's code (subroutine ternary)

! get total (gas + condensed phase HNO3)

      hno3amb=hno3gas+hno3cond

! convert h2o mixing ratio to ppmv:

      h2oppmv= (Max (h2ogas, 0.5d-06)) * 1.d6

! Convert hno3 concentration to ppbv.  Since we are doing an equilibrium
! calculation here, pass the sum of gas and condensed phase HNO3
! to Carslaw's code.  Carslaw's code will decide partitioning between
! gas and condensed phases.

      hno3ppbv=hno3amb*1.d9

! get h2so4 concentration in ppbv from liquid binary sulfate
! surface area density

      call calcsulf_DBC  &
     &  (temp3,pres3c,h2oback*1.d6,lbssad,h2so4gas,rlbs,siglbs,  &
     &   nlbs,constantnlbs)

      h2so4ppbv=h2so4gas*1.d9

! now call subroutine ternary:

      call ternary_DBC(temp3,pres3c,h2oppmv,hno3ppbv,h2so4ppbv,ws,sadwn,v,  &
     &     parthno3,saderr,denssts)

! parthno3 is the fraction of the hno3 molecules that are in gas phase.
! Use this to calculate gas and condensed phase HNO3 amounts
! I have noticed that parthno3 is sometimes slightly greater than 1,
! which could create negative hno3cond values.  Trap for this
! condition

      if(parthno3.gt.1.d0)parthno3=1.d0
      hno3gas=hno3amb*parthno3
      hno3cond=hno3amb-hno3gas

! convert sts volume to surface area density
! use constant STS particle number density assumption if flag set

         logsigstssq=log(sigsts)*log(sigsts)
         rstscm=rsts*1.d-4

! convert volume units from microns cubed/cm^3 to cm^3/cm^3

         v=v*1.d-12

            if(constantn)then

               stssad=(4.d0*pi*nsts)**(1.d0/3.d0)*  &
     &            (3.d0*v)**(2.d0/3.d0)  &
     &             *exp(-1.d0*logsigstssq)

               reff=(3.d0*v/(4.d0*pi*nsts))**(1.d0/3.d0)  &
     &            *exp(-3.d0/2.d0*logsigstssq)

            else

               stssad=3.d0/rstscm*v  &
     &            *exp(-5.d0*logsigstssq/2.d0)

            endif

! return the difference between the calculated sts sad and lbssad,
! providing stssad > lbssad

      stssad=stssad-lbssad

      if(stssad.lt.0.d0)then
         stssad=0.d0
         reff=0.d0
         hno3cond=0.d0
         hno3gas=hno3amb
      endif

      return
      END SUBROUTINE Sad_STS


!     ===================

      SUBROUTINE calcsulf_DBC(t,ptot,qh2o,lbssad,h2so4gas,rlbs,siglbs,  &
     &                     nlbs,constantn)

! subroutine: calcsulf
!
! written by: David B. Considine
!             Code 916
!             NASA Goddard Space Flight Center
!             Greenbelt, MD 20771
! email:      dbc@welkin.gsfc.nasa.gov
! phone:      (301) 286-4299
!
! date:       12/11/96
!
! purpose:
!
!    This subroutine is written for the Global Modeling
!    Initiative of the Atmospheric Effects of Aviation
!    Program.  It is part of the heterogeneous chemistry
!    package to be used in the GMI core model.
!
! description:
!
!    Subroutine calcsulf takes as input the surface area
!    density for the background sulfate aerosol distribution.
!    It calculates the mixing ratio of condensed phase
!    H2SO4 in ppbv consistent with the input surface
!    area density.  It does this using an assumed lognormal
!    size distribution for the sulfate aerosols to convert
!    from surface area density to volume density.  From the
!    volume density and the weight percent of H2SO4 calculated
!    using a binary H2SO4/H2O aerosol composition assumption,
!    the volume mixing ratio is calculated.

! set implicit none

      implicit none

#     include "gmi_phys_constants.h"

! declare variables:

      logical constantn

      real*8 xsb,ks(7),t,ptot,qh2o,pw,msb,denss,wsb,h2so4amu,gpamu,  &
     &       lbssad,h2so4gas,logsiglbssq,siglbs,rlbs,rlbscm,nlbs,  &
     &       v,pi,dens,tt,tice

! set up coefficient arrays

        DATA KS/-21.661d0,2724.2d0,51.81d0,-15732.0d0,47.004d0,  &
     &          -6969.0d0,-4.6183d0/

      save ks

! set parameters

! mass h2so4 in amu

      parameter(h2so4amu=98.08d0)

! grams per amu

      parameter(gpamu=1.66d-24)

! pi:

      parameter(pi=GMI_PI)

! convert from input surface area density to volume density
! use constant liquid binary sulfate particle number density
! assumption if flag set

         logsiglbssq=log(siglbs)*log(siglbs)
         rlbscm=rlbs*1.d-4

            if(constantn)then

               v=lbssad**1.5d0/3.d0/((4.d0*pi*nlbs)**0.5d0)*  &
     &           exp(3.d0/2.d0*logsiglbssq)


            else

               v=rlbscm*lbssad/3.d0*  &
     &            exp(5.d0*logsiglbssq/2.d0)

            endif

! The code converts from volume density (cm^3 aerosol/cm^3 air)
! to mixing ratio of condensed phase H2SO4. First step: calculate
! the density of a binary H2SO4/H2O solution using Carslaw's method.

! first calculate the mole fraction of H2SO4/H2O solution - this
! is taken straight from Carslaw's ternary code

! convert h2o mixing ratio to partial pressure of h2o

        PW=QH2O*1.d-6*PTOT/1013.0d0

        TICE=2668.70d0/(10.4310d0-(LOG(PW)+LOG(760.0d0))/LOG(10.0d0))

        tt = Max (t, tice-3.0d0, 185.0d0)

        XSB = 1.0d0/(2.0d0*(KS(3)+KS(4)/tt))*( -KS(1)-KS(2)/tt-  &
     &  ((KS(1)+KS(2)/tt)**2-4.0d0*(KS(3)+KS(4)/tt)*(KS(5)+KS(6)/tt  &
     &  +KS(7)*LOG(tt)-LOG(PW)))**0.5d0)

!  get molality (moles/mass) of binary solution

        MSB=55.51d0*XSB/(1.0d0-XSB)

! now calculate density - again this is due to Carslaw

        DENSS=(1000.0d0+123.64d0*MSB-5.6d-4*MSB*tt**2  &
     &       -29.54d0*MSB**1.5d0 + 1.814d-4*MSB**1.5d0*tt**2  &
     &       + 2.343d0*MSB**2  -1.487d-3*MSB**2*tt  &
     &       -1.324d-5*MSB**2*tt**2)*1.d-3

! now calculate weight fraction of H2SO4 in this binary solution:

         WSB = MSB*0.098076d0/(1.0d0+MSB*0.098076d0)

! now I can calculate the concentration of H2SO4 in molecules/cm3 air:

      h2so4gas = v*denss*wsb/h2so4amu/gpamu

! now convert to mixing ratio
! first calculate density in total molecules/cm^3

      dens=AVOGAD*ptot*100.d0/GAS_CONST_J/tt/1.d6
      h2so4gas = h2so4gas/dens

      return
      END SUBROUTINE calcsulf_DBC


!     ==================

      SUBROUTINE ternary_DBC(tin,ptot,qh2o,qhno3,qh2so4,ws,sadwn,v,parthno3,  &
     &           saderr,denssts)

! subroutine:   ternary
!
! modified by:  David B. Considine
!               Code 916
!               NASA Goddard Space Flight Center
!               Greenbelt, MD 20771
! email:        dbc@welkin.gsfc.nasa.gov
! phone:        (301) 286-4299
!
! date:         12/03/96
!
! purpose:
!
!     This subroutine is written for the Global Modeling Initiative of
!     the Atmospheric Effects of Aviation program, to be used as part
!     of the heterogeneous chemistry package for the GMI core model.
!     Modifications have been made to the original code written by
!     Dr. Ken Carslaw (see documentation below) to fit as a subroutine
!     into the structure of the GMI model. As few changes as possible
!     have been made to the original code.  All changes have been documented
!     below.
!
! input variables:
!
! t - temperature
! ptot - total pressure (input in mb)
! qh2o - H2O mixing ratio (ppmv)
! qhno3 - HNO3 mixing ratio (ppbv)
! qh2so4 - mixing ratio of h2so4 (ppbv)
! qhcl - mixing ratio of hcl (ppbv)
! qhocl - mixing ratio of hocl (ppbv)
! ==========================================================================================
!          FORTRAN77 CODE TO CALCULATE COMPOSITION OF AQUEOUS HNO3/H2SO4/HCL/HOCL
!         STRATOSPHERIC AEROSOLS (CARSLAW, LUO, PETER - GEOPHYS. RES. LETT., 1995)
!                           CARSLAW@NIKE.MPCH-MAINZ.MPG.DE
!
!                           DISTRIBUTION VERSION 1 (2 AUG 1995)
!
!         DR KEN CARSLAW
!         MAX-PLANCK-INSTITUT FUER CHEMIE
!         POSTFACH 3060
!         MAINZ 55020
!         GERMANY
!         TEL: (49) (0)6131 305 333
!         FAX: (49) (0)6131 305 328
!
! ==========================================================================================
! HNO3/H2SO4 COMPOSITION BASED ON THERMODYNAMIC MODEL OF THE SYSTEM
! HCL/HBR/HNO3/H2SO4/H2O {CARSLAW ET AL, J. PHYS. CHEM., 1995 - A THERMODYNAMIC MODEL
! OF THE SYSTEM HCL-HNO3-H2SO4-H2O, INCLUDING SOLUBILITIES OF HBR, FROM <200K TO 328K.
! (ALSO AVAILABLE IN FORTRAN77)}
! HCL SOLUBILITY PARAMETRISATION FROM LUO ET AL., GRL, 1995
! HOCL SOLUBILITY PARAMETRISATION FROM HUTHWELKER ET AL., J. AT. SCI., 1995
!
! ****************************************************************************
! THE MODEL IS VALID FOR 2E-5 MB<PW<2E-3 MB (PW IS THE WATER PARTIAL PRESSURE)
! THE UPPER TEMPERATURE LIMIT IS 240K
! THE LOWER TEMPERATURE LIMIT IS 185K, OR TICE-3K, WHICH EVER IS THE HIGHER
! HNO3 SOLUBILITIES ARE CALCULATED TO A MAXIMUM TEMPERATURE OF 215K
! SOLUBILITIES ARE CALCULATED ON A MOLALITY BASIS (MOLES OF SOLUTE PER KG OF WATER)
! ****************************************************************************
!
! THE SOLUBILITIES OF HCL AND HOCL ARE ASSUMED NOT TO AFFECT THE UPTAKE OF HNO3 OR H2O.
! THIS INTRODUCES ONLY SMALL ERRORS AT THE VERY LOWEST TEMPERATURES, BUT FOR A FULL
! CALCULATION WHERE INTERACTIONS BETWEEN ALL SPECIES IN SOLUTION ARE CONSIDERED, USE
! THE MODEL OF CARSLAW ET AL., GIVEN ABOVE.
! THE PROGRAM HAS BEEN WRITTEN SO THAT THE DIFFERENT STEPS ARE CLEAR, AND IS THEREFORE
! NOT CODED FOR OPTIMUM CALCULATION SPEED.  THIS IS LEFT FOR THE USER TO DO IF DESIRED.
! ==========================================================================================
!

! dbc - remove original program name
!        PROGRAM ANALYTIC

! dbc - set implicit none and declare all variables

        implicit none

#       include "gmi_phys_constants.h"

! dbc - change to real*8 variable definition

!        REAL NS,MSB,MNB,MS,MN,MCL,MHOCL
!        REAL QN(10),QS(10),KN(7),KS(7)

        real*8 ns,msb,mnb,ms,mn
        real*8 qn(10),qs(10),kn(7),ks(7)

! dbc - declare parameter r real*8

        real*8 r
        PARAMETER (R=8.205d-5)

! dbc - declare logical for error

        logical saderr

! dbc - input variable declaration

        real*8 t,tin,ptot,qh2o,qhno3,qh2so4

! dbc - declare currently undeclared internal variables

        real*8 pw,pn0,tt,tr,pr,tice,xsb,hsb,xnb,  &
     &         hnb,a,b,c,phi,pi,pn,ws,sadwn,parthno3,v,  &
     &         denssts,density_dbc


! ======================================================================
        DATA QN/14.5734d0,0.0615994d0,-1.14895d0,0.691693d0,-0.098863d0,  &
     &  0.0051579d0,0.123472d0,-0.115574d0,0.0110113d0,0.0097914d0/
        DATA QS/14.4700d0,0.0638795d0,-3.29597d0,1.778224d0,-0.223244d0,  &
     &  0.0086486d0,0.536695d0,-0.335164d0,0.0265153d0,0.0157550d0/
        DATA KN/-39.136d0,6358.4d0,83.29d0,-17650.0d0,198.53d0,  &
     &  -11948.d0,-28.469d0/
        DATA KS/-21.661d0,2724.2d0,51.81d0,-15732.0d0,47.004d0,  &
     &       -6969.0d0,-4.6183d0/
! ======================================================================

        save qn, qs, kn, ks

! dbc - Initially set saderr to false.  Subroutine can change this if a problem
!       occurs

        saderr=.false.

! initially set t to tin

        t=tin

        PW=QH2O*1.d-6*PTOT/1013.0d0

        PN0=PTOT/1013.0d0*QHNO3*1.d-9

! dbc ns is the number of moles of H2SO4 per m^3
!     the factor of 100 converts from millibar to Pascal

        NS=PTOT*100.0d0*QH2SO4*1.d-9/8.314d0/T

        PR=LOG(PW)+18.4d0
        TICE=2668.70d0/(10.4310d0-(LOG(PW)+LOG(760.0d0))/LOG(10.0d0))
!
! code amended by D. Considine: Carslaw's routine only valid for
! temperatures greater 185 which are also greater than 3 degrees
! below frostpoint. So, when temps are less than this, set temps
! to minimum values.

        if(t.lt.tice-3.d0)t=tice-3.d0
        if(t.lt.185.d0)t=185.d0

        TT=R*T*NS
        TR=1.0d4/T-43.4782608d0

        IF((T.GE.TICE-3.0d0).AND.(T.GE.185.0d0))THEN
! THE H2SO4/H2O PURE SOLUTION CONCENTRATION
           XSB = 1.0d0/(2.0d0*(KS(3)+KS(4)/T))*( -KS(1)-KS(2)/T-  &
     &     ((KS(1)+KS(2)/T)**2-4.0d0*(KS(3)+KS(4)/T)*(KS(5)+KS(6)/T  &
     &     +KS(7)*LOG(T)-LOG(PW)))**0.5d0)
           MSB=55.51d0*XSB/(1.0d0-XSB)
!
              IF((T.LE.215.0d0).AND.(PN0.GT.0.0d0))THEN
! THE HNO3/H2SO4/H2O SOLUTION COMPOSITION
                 HSB= QS(1)+QS(2)*TR**2+(QS(3)+QS(4)*TR+QS(5)*TR**2  &
     &           +QS(6)*TR**3)*PR + (QS(7)+QS(8)*TR+QS(9)*TR**2)  &
     &           *PR**2+QS(10)*TR*PR**3
                 HSB=EXP(HSB)
                 XNB = 1.0d0/(2.0d0*(KN(3)+KN(4)/T))*( -KN(1)-KN(2)/T-  &
     &           ((KN(1)+KN(2)/T)**2-4.0d0*(KN(3)+KN(4)/T)*(KN(5)+  &
     &           KN(6)/T+KN(7)*LOG(T)-LOG(PW)))**0.5d0)
                 MNB=55.51d0*XNB/(1.0d0-XNB)
                 HNB= QN(1)+QN(2)*TR**2+(QN(3)+QN(4)*TR+QN(5)*TR**2  &
     &           +QN(6)*TR**3)*PR + (QN(7)+QN(8)*TR+QN(9)*TR**2)  &
     &           *PR**2+QN(10)*TR*PR**3
                 HNB=EXP(HNB)
                 A=(TT*HNB*MNB**2 - TT*HSB*MNB*MSB - 2.0d0*MNB**2*MSB  &
     &           + MNB*MSB**2 + HNB*MNB*MSB*PN0 - HSB*MSB**2*PN0)/  &
     &           (MNB**2 - MNB*MSB)
                 B=MSB*(-2.0d0*TT*HNB*MNB+TT*HSB*MSB+MNB*MSB-  &
     &           HNB*MSB*PN0)/(MNB-MSB)
                 C=(TT*HNB*MNB*MSB**2)/(MNB - MSB)
                 PHI= ATAN(SQRT(4.0d0*(A**2-3.0d0*B)**3-(-2.0d0*A**3  &
     &           +9.0d0*A*B-27.0d0*C)**2)/(-2.0d0*A**3+9.0d0*A*B  &
     &           -27.0d0*C) )
                 PI = GMI_PI
              IF(PHI.LT.0d0) PHI=PHI+PI
                 MS=-1.0d0/3.0d0  &
     &           *(A+2.0d0*SQRT(A**2-3.0d0*B)*COS((PI+PHI)/3.0d0))
                 MN = MNB*(1.0d0-MS/MSB)
                 PN = MN/(HNB*MN/(MN+MS)+HSB*MS/(MN+MS))
                 WS = MS*0.098076d0/(1.0d0+MS*0.098076d0+MN*0.063012d0)
                 sadwn =  &
     &             MN*0.063012d0/(1.0d0+MS*0.098076d0+MN*0.063012d0)
                 PARTHNO3=(1.0d0-(PN0-PN)/PN0)

                 mn = Max (mn, 0.0d0)

              ELSE
! ASSUME SOLUTION IS PURE H2SO4/H2O
                 MS = MSB
                 MN = 0.0d0
                 WS = MSB*0.098076d0/(1.0d0+MSB*0.098076d0)
                 sadwn = 0.0d0
                 PN = PN0
                 PARTHNO3=1.0d0
              ENDIF
           V=NS*98.076d0/WS/DENSITY_DBC(MS,MN,T)*1.d6
! TOTAL LIQUID AEROSOL VOLUME IN UM3/CM3
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! dbc introduce variable to return density

           denssts=DENSITY_DBC(MS,MN,T)

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ===============================================================================
! THE SOLUBILITY OF HCL
! H* (MOL/KG/ATM) ADAPTED FROM LUO ET AL., GRL, 1995.
! CALCULATED CONCENTRATIONS ASSUME THAT HCL IS TRACE COMPONENT OF THE AEROSOL.
!c            IF(PHCL0.GT.0.0d0)THEN
!c               W1 = sadwn
!c               W2 = WS
!c               HHCL = EXP(-(21.0d0+46.610d0*W1+4.0690d0*W2
!c   c           -4.8370d0*SQRT(W1)
!c   >           +2.1860d0*SQRT(W2)-63.00d0*W1**2-40.170d0*W1*W2
!c   c           -1.5710d0*W2**2)
!c   >           -1.0/T*(-7437.0d0-8327.80d0*W1+1300.90d0*W2
!c   >           -1.0d0/T*(-7437.0d0-8327.80d0*W1+1300.90d0*W2
!c   c           +1087.20d0*SQRT(W1)
!c   >           -242.710d0*SQRT(W2)+18749.0d0*W1**2+18500.0d0*W1*W2
!c   >           +5632.0d0*W2**2)-LOG(W1+0.610d0*W2)
!c   >           -LOG(36.461d0/(1000.0d0+98.076d0*MS+63.012d0*MN)))
!c   c           *1.013d3
!c               MCL=(1.0d0/R/T*PHCL0)/(NS/MS + 1.0d0/R/T/HHCL)
!c               WCL=MCL*36.461d0/(1000.0d0+98.076d0*MS+63.012d0*MN)
!c               PHCL=MCL/HHCL
!c               PARTHCL=1.0d0-(PHCL0-PHCL)/PHCL0
!c            ENDIF
! ===============================================================================
! THE SOLUBILITY OF HOCL
! H* (MOL/KG/ATM) FROM HUTHWELKER ET AL., J. AT. SCI, 1995
! AS AN APPROXIMATION, ASSUME H* DEPENDS UPON TOTAL MOLALITY
!c            IF(PHOCL0.GT.0d0)THEN
!c            IF(PHOCL0.GT.0.d0)THEN
!c               HHOCL=EXP(6.49460d0-(-0.041070d0+54.56d0/T)*(MS+MN)
!c   >           -5862.0d0*(1.0d0/298.150d0-1.0d0/T))
!c               MHOCL=(1.0d0/R/T*PHOCL0)/(NS/MS + 1.0d0/R/T/HHOCL)
!c               WHOCL=MHOCL*52.46d0/(1000.0d0+98.076d0*MS+63.012d0*MN)
! SOLUBILITY OF HOCL IS LOW ENOUGH TO IGNORE GAS PHASE REMOVAL
!c            ENDIF
! ===============================================================================

        ELSE

! dbc: error flag to determine if any conditions have been violated

         saderr=.true.

        ENDIF

        END SUBROUTINE ternary_DBC
!
! ============================================================================

        FUNCTION DENSITY_DBC(CS,CN,T)

! - dbc declare variables

        implicit none

        real*8 cs,cn,t,denss,densn,density_dbc

! LIQUID SOLUTION DENSITY G/CM3

        DENSS=1000.0d0+123.64d0*CS-5.6d-4*CS*T**2  &
     &       -29.54d0*CS**1.5d0 + 1.814d-4*CS**1.5d0*T**2  &
     &       + 2.343d0*CS**2  -1.487d-3*CS**2*T  &
     &       -1.324d-5*CS**2*T**2

        DENSN=1000.0d0+85.107d0*CN-5.043d-4*CN*T**2  &
     &       -18.96d0*CN**1.5d0 + 1.427d-4*CN**1.5d0*T**2  &
     &       + 1.458d0*CN**2  -1.198d-3*CN**2*T  &
     &       -9.703d-6*CN**2*T**2

        DENSITY_DBC=0.001d0/((1.0d0/DENSS*CS/(CS+CN)  &
     &  +1.0d0/DENSN*CN/(CS+CN)))

        RETURN
        END FUNCTION DENSITY_DBC


!     ===================

      SUBROUTINE Sediment_DBC  &
     &  (dz, dt, temp3, pres3c,  &
     &   hno3cond, h2ocond, dehyd, dehyd_opt, &
     &   rsts, rnat, rice, reffsts, reffnat, reffice, sigsts,  &
     &   signat, sigice, constantnsts, constantnnat, constantnice,  &
     &   densnat, densice, vfall, calcsts, denssts, &
     &   ilo, ihi, julo, jhi, i1, i2, ju1, j2, k1, k2)

      implicit none

#     include "gmi_phys_constants.h"


      integer, intent(in) :: ilo, ihi, julo, jhi, dehyd_opt
      integer, intent(in) :: i1, i2, ju1, j2, k1, k2
      real*8   dz      (i1:i2,   ju1:j2,   k1:k2)
      real*8   dt
      real*8   temp3   (ilo:ihi, julo:jhi, k1:k2)
      real*8   pres3c  (ilo:ihi, julo:jhi, k1:k2)
      real*8   hno3cond(i1:i2,   ju1:j2,   k1:k2)
      real*8   h2ocond (i1:i2,   ju1:j2,   k1:k2)
      real*8   dehyd   (i1:i2,   ju1:j2,   k1:k2)
      real*8   rsts, rnat
      real*8   rice
      real*8   reffsts (i1:i2,   ju1:j2,   k1:k2)
      real*8   reffnat (i1:i2,   ju1:j2,   k1:k2)
      real*8   reffice (i1:i2,   ju1:j2,   k1:k2)
      real*8   sigsts, signat, sigice
      logical  constantnsts, constantnnat, constantnice
      real*8   densnat, densice
      real*8   vfall   (i1:i2,   ju1:j2,   k1:k2)
      logical  calcsts
      real*8   denssts (i1:i2,   ju1:j2,   k1:k2)


      logical  go

      integer  idtsed
      integer  ifac
      integer  il, ij, ik
      integer  itime
      integer  maxht, minht

      real*8   dens, densp
      real*8   dtsed
      real*8   dzmin
      real*8   fac
      real*8   Fall_Vel
      real*8   fluxcorr
      real*8   loss
      real*8   prod
      real*8   vfallmax


!     -----------------------------------------
!     Loop over longitude and latitude indices.
!     -----------------------------------------

!     ===============
      do il = i1, i2
      do ij = ju1, j2
!     ===============

! at each lat and lon, first decide if sedimentation needs
! to be calculated

      go = .false.

      do 11 ik = k1,k2

         if(hno3cond(il,ij,ik).gt.0.d0)then
            go = .true.
            goto 12
         endif

 11   continue

 12   if(go)then

! determine fall velocity field and max fall velocity,
! and upper alt of condensed phase

      maxht=k1
      minht=k2

      vfallmax=0.d0

      do 13 ik = k1,k2

         if(hno3cond(il,ij,ik).gt.0.d0  &
     &      .or.h2ocond(il,ij,ik).gt.0.d0)then

            if(ik.gt.maxht)maxht=ik
            if(ik.lt.minht)minht=ik

         endif

 13   continue

      if(maxht.eq.k2)maxht=k2-1
      if(minht.eq.k1)minht=k1+1

! initialize dzmin with an arbitrary large number

      dzmin=1.d7

      do 14 ik = minht,maxht

! find the minimum gridbox height

      if(dz(il,ij,ik).lt.dzmin)dzmin=dz(il,ij,ik)


! If there is condensed phase h2o at the gridpoint, then
! code assumes all condensed hno3 is contained in h2o particles.
! Therefore, sediment hno3 using h2o particle size distribution.

         if(h2ocond(il,ij,ik).gt.0.d0)then

            fluxcorr=exp(8.d0*log(sigice)*log(sigice))

         else

            if(calcsts)then

               fluxcorr=exp(8.d0*log(sigsts)*log(sigsts))

            else

               fluxcorr=exp(8.d0*log(signat)*log(signat))

            endif

         endif

! calculate local density using ideal gas law

         dens=AVOGAD*pres3c(il,ij,ik)*100.d0/GAS_CONST_J/  &
     &        temp3(il,ij,ik)/1.d6

         if(h2ocond(il,ij,ik).gt.0.d0)then

! calculate fall velocity using ice particle size distribution when there is
! any condensed h2o.

            if(constantnice)then

               vfall(il,ij,ik)=fluxcorr*  &
     &         Fall_Vel(reffice(il,ij,ik)*1.d4,densice,  &
     &                 dens,temp3(il,ij,ik))

            else

               vfall(il,ij,ik)=fluxcorr*  &
     &         Fall_Vel(rice,densice,dens,temp3(il,ij,ik))

            endif

         else

            if(calcsts)then

               if(constantnsts)then

                  vfall(il,ij,ik)=fluxcorr*  &
     &            Fall_Vel(reffsts(il,ij,ik)*1.d4,  &
     &            denssts(il,ij,ik),  &
     &            dens,temp3(il,ij,ik))

               else

                  vfall(il,ij,ik)=fluxcorr*  &
     &            Fall_Vel(rsts,denssts(il,ij,ik),  &
     &            dens,temp3(il,ij,ik))

               endif

            else

               if(constantnnat)then

                  vfall(il,ij,ik)=fluxcorr*  &
     &            Fall_Vel(reffnat(il,ij,ik)*1.d4,densnat,  &
     &            dens,temp3(il,ij,ik))

               else

                  vfall(il,ij,ik)=fluxcorr*  &
     &            Fall_Vel(rnat,densnat,dens,temp3(il,ij,ik))

               endif
            endif
        endif

! vfallmax is recalculated at each latitude and longitude

         if(vfall(il,ij,ik).gt.vfallmax)vfallmax=  &
     &   vfall(il,ij,ik)

 14   continue

! calculate time step for sedimentation.  The time step is calculated
! so that no gridbox can fall through the gridbox below it in one
! sedimentation time step.

      dtsed=dzmin/vfallmax

      if(dt/dtsed.le.1.d0)then

         idtsed=1
         dtsed=dt             ! no need to change timestep

      else

         fac    = dt / dtsed
         ifac   = fac
         idtsed = ifac + 1      ! number of iterations necessary

         if(idtsed.gt.10)idtsed=10 ! limit number of iterations
         dtsed=dt/(idtsed*1.d0)
         vfallmax=dzmin/dtsed

      endif


!     -----------------------------------------------------------------
!     Do sedimentation.  Loop idtsedal times, use time step of dtsedal:
!     -----------------------------------------------------------------

      do 16 itime = 1, idtsed
        do 15 ik = 1, maxht

!         -----------------------------------------------------------------
!         Calculate production and loss terms for sedimentation of hno3;
!         these terms should be calculated using the fall velocity
!         interpolated to the gridbox edges.  This is not being done simply
!         for ease of coding.
!         -----------------------------------------------------------------

          dens  = AVOGAD * pres3c(il,ij,ik)   * 100.0d0 / GAS_CONST_J /  &
     &            temp3(il,ij,ik)   / 1.0d6

          densp = AVOGAD * pres3c(il,ij,ik+1) * 100.0d0 / GAS_CONST_J /  &
     &            temp3(il,ij,ik+1) / 1.0d6

          prod = dtsed / dz(il,ij,ik+1) *  &
     &           Min (vfall(il,ij,ik+1), vfallmax) *  &
     &           hno3cond(il,ij,ik+1) * densp

          loss = dtsed / dz(il,ij,ik) *  &
     &           Min (vfall(il,ij,ik), vfallmax) *  &
     &           hno3cond(il,ij,ik) * dens

!         ---------------------------
!         Actually sediment the hno3.
!         ---------------------------

          hno3cond(il,ij,ik) =  &
     &      hno3cond(il,ij,ik) + ((prod - loss) / dens)

          if (hno3cond(il,ij,ik) .lt. 0.0d0) then
            hno3cond(il,ij,ik) = 0.0d0
          end if

          if ((h2ocond(il,ij,ik)   .gt. 0.0d0) .or.  &
     &        (h2ocond(il,ij,ik+1) .gt. 0.0d0)) then

!           -----------------------------------------------------------
!           Calculate the changes to the dehyd array.  Sedimentation of
!           h2o into a gridbox constitutes a loss of dehydration.
!           Sedimentation out of a gridbox is a production term.
!           -----------------------------------------------------------

            prod  = dtsed / dz(il,ij,ik) *  &
     &              Min (vfall(il,ij,ik), vfallmax) *  &
     &              h2ocond(il,ij,ik)

            loss  = dtsed / dz(il,ij,ik+1) *  &
     &              Min (vfall(il,ij,ik+1), vfallmax) *  &
     &              h2ocond(il,ij,ik+1) * densp

            IF(dehyd_opt == 1) THEN
	     dehyd(il,ij,ik) = dehyd(il,ij,ik) + prod - (loss / dens)
	    END IF

            h2ocond(il,ij,ik) =  &
     &        h2ocond(il,ij,ik) + (loss / dens) - prod

            if (h2ocond(il,ij,ik) .lt. 0.0d0) then
              h2ocond(il,ij,ik) = 0.0d0
            end if

          end if

 15     continue
 16   continue


      end if ! matches if (go)

!     ======
      end do
      end do
!     ======


      return
      END SUBROUTINE Sediment_DBC


!     ================

      FUNCTION Fall_Vel(rad,density,ndens,tem)

! Function fallvel ---------- David B. Considine, 11/11/96

! Code to calculate the fall velocity for a particle of radius r,
! which is passed to the function.  The fall velocity is calculated
! according to Kasten [1968]. (Kasten, F., Falling Speed of Aerosol
! Particles, J. Appl. Met., Oct, 1968
!
! output
!
! fallvel - speed of falling particle in cm/s
!
! variables
!
! rad - input radius of falling particle (microns)
! density - input density of particle (g/cm^3)
! ndens - atmospheric number density (molecules/cm^3)
! dynvis - dynamic viscosity
! tem - temperature (degK)
! mfp - mean free path
!
! parameters
!
! sigsq - square of effective collision diameter
! bet -
! s - sutherland's constant
! a,b,cc - constants used in Kasten's formula
! grav - force of gravity (m/s^2)

! declare variables

      implicit none

#     include "gmi_phys_constants.h"

      real*8 radius,rho,ndens,tem,sigsq,mfp,bet,s,dynvis,a,b,cc,  &
     &     grav,Fall_Vel,rad,density

! set internal parameters

      parameter(sigsq=1.3323d-19)
      parameter(bet=1.458d-6)
      parameter(s=110.4d0)
      parameter(a=1.249d0)
      parameter(b=0.42d0)
      parameter(cc=0.87d0)
      parameter(grav=GMI_G)

! don't cause problems if subroutine called with rad or density
! equal to 0:

      if(rad.eq.0.d0.or.density.eq.0.d0)then

         Fall_Vel=0.d0
         return

      else

! code assumes radius is in microns and rho is in g/cm^3.
! convert radius to meters and rho to kg/m^3:

         radius=rad*1.d-6
         rho=density*1.d3

! Calculate mean free path. According to the U.S. Standard Atmosphere,
! 1976, the mean free path is given by the relationship: mfp =
! sqrt(2)/(2*pi*sig^2*N), where sig is the effective collision diameter
! (3.65e-10 meters) and N is the number density (units of #/m^3)
! mfp units are 1/meters


         mfp=.22508d0/(sigsq*ndens*1.d6)

! Calculate dynamic viscosity. The formula is:
! dynvis=bet*T^(3/2)/(T+S), where bet is ..., s is Sutherland's
! constant (110.4 degK).  The units are Newtons/m^2*s.

         dynvis=bet*tem**1.5d0/(tem+s)

! Calculate fall velocity.  The factor of 100 results in cm/s
! output.

         Fall_Vel=0.2222d0*rho*radius*radius*grav/dynvis*100.d0*  &
     &                    (1.+ mfp/radius*(a+b*exp(-cc*radius/mfp)))


         return

      endif

      END FUNCTION Fall_Vel

