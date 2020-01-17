
!  $Id: Lightning_mod.F90,v 0.1 2019/10/07 15:52:00 mrdamon Exp $

module Lightning_mod

#include "MAPL_Generic.h"

!----------------------------------------------------------------------
!BOP

! !MODULE: 

!    Lightning_mod -- Container for the GEOS lightning utilities.

! !USES:

use MAPL_Mod
use Lightning_Toolbox_Mod, only : CalcFlashRate 

   implicit none
   private

  
! PUBLIC MEMBER FUNCTIONS:
!

   PUBLIC :: getLightning
   PUBLIC :: computeFlashRate
   PUBLIC :: flashfit_flipped ! fields passed are vertially flipped from GEOS orientation 
   PUBLIC :: hemcoFlashrate
   PUBLIC :: flash_rate_Lopez
   PUBLIC :: emiss_lightning
   PUBLIC :: partition
   PUBLIC :: computeCAPE
   PUBLIC :: readLightRatioGlobalData
   PUBLIC :: BUOYANCY


! !PARAMETERS:

   real*8, parameter :: AVOGAD  = 6.0221367d+23 ! Avogadro number (mole^-1)
   real*8, parameter :: CONVFAC = 2.33e-23      ! 14g of N per mole of NO/ 6.02e23 molecules of NO/mole of NO.
   real*8, parameter :: PRODFAC = 1.0e26        ! (joules per flash) * molecules of NO per joule
   real*8, parameter :: Pa2hPa  = 0.01
   integer,parameter  :: hp = KIND( REAL( 0.0, 8 ) ) ! HEMCO type


! !DESCRIPTION:

!  This module contains parametereizations and utilities pertaining to the calculation of
!  flashrate (or strokerate) and NOx from lightning in the Earth's atmosphere.
!  The routines for calculation flashrate / stroke rate are:
! 
!    flashfit_flipped - adapted from the offline GMI-CTM model (Dale Allen)
!    computeFlashRate - adapted from MOIST / CTM cinderalla component
!    hemcoFlashrate - GEOS-Chem's flashrate calculation 
!    flash_rate_Lopez (WARNING: this routine is producing invalid results!)
!   
!  Please review the documentation below for each routine carefully. 
!
!
!
! !REVISION HISTORY:
!
! March 17 2018 - Megan Damon. Integrating code from the GMI-CTM.
! October 15 2018 - Megan Damon. Integrating code from GEOS-5 for computeFlashRate
! November 17 2018 - Megan Damon. Bringing into GEOS from standalone development
! Ocotober 7 2019 - Megan Damon. Preparing delivery to Christoph Keller
!
!EOP
!-------------------------------------------------------------------------


contains


!-----------------------------------------------------------------------
!BOP
! !IROUTINE: getLightning
!
! !DESCRIPTION:
!  
!
!
!  ORIGIN AND CONTACT
!
! !REVISION HISTORY:
! 7 Nov 2019 Damon       First crack
!EOP
!-----------------------------------------------------------------------

subroutine getLightning (flashSource, Q, TH, PLE, CAPE, BYNCY)

  character(len=*), intent(in) :: flashSource ! MOIST, FIT, HEMCO, or Lopez

  real, dimension(:,:,:),   intent(in)  :: Q     ! specific humidity (kg kg-1)
  real, dimension(:,:,:),   intent(in)  :: TH    ! potential temperature (K)
  real, dimension(:,:,:),   intent(in)  :: PLE   ! edge pressures (Pa)

  real, dimension(:,:),     intent(out)  :: CAPE       ! 
  real, dimension(:,:,:),   intent(out)  :: BYNCY

  logical, save :: first = .true.
  integer, save :: IM, JM, LM, KM
  
  integer :: STATUS

  integer :: shape3d (3)


  print*, "I am getLightning, flash source: ", flashSource

  if (first) then

     print*, flashSource(1:5), " flash rate calculation"

!     shape3d = shape(TH)

!     IM = shape3d(1)
!     JM = shape3d(2)
!     LM = shape3d(3)
!     KM = LM + 1

     first = .False.
  endif

  !-----------------------------------------------------------------------
  ! MOIST-derived flashrate calculation (Dale Allen, old)
  !-----------------------------------------------------------------------
  if (flashSource == "MOIST") then

     print*, "I will call: computeCAPE, computeFlashRate"

     !call computeCAPE (TH, Q, PLE, CAPE, BYNCY, IM, JM, LM)

  !-----------------------------------------------------------------------
  ! Dale Allen flashrate calculation adapted from the GMI-CTM
  !-----------------------------------------------------------------------
  else if (flashSource == "FIT") then

     print*, "I will call ExtData?"
     print*, "I will call flashfit"

  
  !-----------------------------------------------------------------------
  ! HEMCO flash rate
  !-----------------------------------------------------------------------
  else if (flashSource == "HEMCO") then

     print*, "I will call computeCAPE"
     print*, "I will call the wrapper routine: hemcoFlashrate"


  !-----------------------------------------------------------------------
  ! Lopez (2016?) flash rate from ECMWF
  !-----------------------------------------------------------------------         
  else if (flashSource == "LOPEZ") then

     print*, "I will call BUOYANCY"
     print*, "I will call flash_rate_Lopez"

       
  else 
     print*, ""
     print*, "Flashrate source not supported!"
     print*, ""
     stop
  endif

end subroutine getLightning

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: BUOYANCY
!
! !DESCRIPTION:
!  Calculates buoyancy, CAPE (Joule/kg), and the inhibition.
!
!
!  ORIGIN AND CONTACT
!  Extracted from GEOS_MoistGridCompMod
!  Guidance from Saulo Freitas:
!     I’d recommend you stick with the Max’s version (1.135). 
!     However, to get units for CAPE as Joule/kilogram you need to divide 
!     “HC” (line 9122) and “HSe” (line 9125) by MAPL_CP.
!
! !REVISION HISTORY:
! 1 Nov 2019 Damon       Bringing routine to GEOS_Shared
!EOP
!-----------------------------------------------------------------------


  subroutine BUOYANCY( T, QLM, QS, DQS, DZ, ZLO,    BUOY, CAPE, INHB )



    ! !DESCRIPTION: Computes the buoyancy $ g \frac{T_c-T_e}{T_e} $ at each level
    !  for a parcel raised from the surface. $Hc$ is the MSE of
    !  the parcel and $HSe$ is the MSE of the environment.

!    integer,                  intent(in)  :: IM,JM,LM
    real, dimension(:,:,:),   intent(in)  :: T, QS, DQS, DZ, ZLO
    real, dimension(:,:  ),   intent(in)  :: QLM
    real, dimension(:,:,:),   intent(out) :: BUOY
    real, dimension(:,:),     intent(out) :: CAPE, INHB

    integer   :: I,J,L, LM
    real      :: HC, HSe, DelT
    logical   :: UNSTABLE

    print*, "min/max T: ", minval(T), maxval(T)
    print*, "min/max QLM: ", minval(QLM), maxval(QLM)
    print*, "min/max QS: ", minval(QS), maxval(QS)
    print*, "min/max DQS: ", minval(DQS), maxval(DQS)
    print*, "min/max DZ: ", minval(DZ), maxval(DZ)
    print*, "min/max ZLO: ", minval(ZLO), maxval(ZLO)



    LM = size(BUOY,3)

    BUOY = MAPL_UNDEF

    do J=1,size(BUOY,2)
       do I=1,size(BUOY,1)

          CAPE(I,J) = 0.0
          INHB(I,J) = 0.0

!          print*, "I,J,LM: ", I,J,LM
!          print*, "CAPE, INHB is zero"

          UNSTABLE  = .false.
          HC        = MAPL_CP*T(I,J,LM) + MAPL_GRAV*ZLO(I,J,LM) + MAPL_ALHL*QLM(I,J)
          HC = HC / MAPL_CP

!          print*, "Before loop: ', HC: ", Hc
          do L=LM-1,1,-1
             HSe    = MAPL_CP*T(I,J,L ) + MAPL_GRAV*ZLO(I,J,L ) + MAPL_ALHL*QS(I,J,L)
             HSe = HSe / MAPL_CP

!             print*, "After loop HC, HSe: ", HC, HSe

             if(HC > HSe) then
!                print*, "HC > HSe"
                DelT        = (HC - HSe) / (1.+ (MAPL_ALHL/MAPL_CP)*DQS(I,J,L))
                BUOY(I,J,L) = MAPL_GRAV * (DelT/T(I,J,L))
                CAPE(I,J)   = CAPE(I,J) + BUOY(I,J,L)*DZ(I,J,L)
                UNSTABLE    = .true.
             else
                if(.not.UNSTABLE) then
                   if(ZLO(I,J,L)>1000.0) then
!                      print*, "ZLO > 1000"
                      exit  ! for economy, since CAPE is non-sense if there is this much inhibition
                   else
!                      print*, "Not exiting"
                      DelT        = (HC - HSe) / (1.+ (MAPL_ALHL/MAPL_CP)*DQS(I,J,L))
                      BUOY(I,J,L) = MAPL_GRAV * (DelT/T(I,J,L))
                      INHB(I,J)   = CAPE(I,J) + BUOY(I,J,L)*DZ(I,J,L)
                   end if
                else
 !                  print*, "UNSTABLE"
                   exit ! if UNSTABLE and delt is negative (i.e., Stop at first LNB)
                end if
             end if
          end do

          if(CAPE(I,J) == 0.0) then
             CAPE(I,J) = MAPL_UNDEF
             INHB(I,J) = MAPL_UNDEF
          end if

       end do
    end do

    print*, "min/max BUOY: ", minval(BUOY), maxval(BUOY)
    print*, "min/max CAPE: ", minval(CAPE), maxval(CAPE)
    print*, "min/max INHB: ", minval(INHB), maxval(INHB)

  end subroutine BUOYANCY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!-----------------------------------------------------------------------
!BOP
! !IROUTINE: computeFlashRate
!
! !DESCRIPTION:
!  Generate lightning flash rates [km$^{-2}$ s$^{-1}$] using a six-variable polynomial fit.\\
!
!
!  ORIGIN AND CONTACT
!  Dr. Dale Allen, Associate Research Scientist
!  Dept. of Atmospheric and Oceanic Science
!  University of Maryland
!  College Park, MD 20742
!  301-405-7629 (ph); 301-314-9482 (fax)
!  http://www.meto.umd.edu/~allen
!  
!
!  !FORMULATION NOTES
!  Predictor variables are set to zero where CN_PRCP is zero or where the 
!   optical depth cloud top height is less than 5.5 km.
!  The fit returns flash rates in units km$^{-2}$ day$^{-1}$.  Convert to 
!   km$^{-2}$ s$^{-1}$ for the export state.
!
!
!  OTHER NOTES OF INTEREST
!  MOIST sets CNV_TOPP to zero if there is an absence of convection.
!
! !REVISION HISTORY:
! 30 Nov 2011 Nielsen     First crack
! 29 Feb 2012 Nielsen     Accomodate CNV\_TOPP MAPL\_UNDEF for and after Fortuna-2\_5\_p4
! 04 Nov 2014 Kouatchou   Adapted the subroutine for GEOSctm
! 18 Nov 2018 Damon       Bringing routine to GEOS_Shared
!EOP
!-----------------------------------------------------------------------

  subroutine computeFlashRate (STATE, nc, lm, TS, CCTP, FROCEAN, CN_PRCP, CAPE, &
       CNV_MFC, TH, PLE, ZLE, strokeRate, RC)

! ! USES
    implicit none
!
! !INPUT PARAMETERS:
    integer, intent(in) :: nc ! Number of cells
    integer, intent(in) :: lm ! Number of layers

    real, intent(in), dimension(nc) :: TS          ! Surface temperature [K]
    real, intent(in), dimension(nc) :: CCTP        ! Convective cloud top pressure [Pa] with MAPL_UNDEFS
    real, intent(in), dimension(nc) :: FROCEAN     ! Areal ocean fraction
    real, intent(in), dimension(nc) :: CN_PRCP     ! Convective precipation [kg m^{-2} s^{-1}]
    REAL, INTENT(IN), DIMENSION(nc) :: CAPE     ! Convective available potential energy [J m^{-2}]

    real, intent(in), dimension(nc,lm) :: TH        ! Potential temperature [K]
    real, intent(in), dimension(nc,0:lm) :: CNV_MFC ! Convective mass flux [kg m^{-2} s^{-1}]
    real, intent(in), dimension(nc,0:lm) :: PLE     ! Layer interface pressures  [Pa]
    real, intent(in), dimension(nc,0:lm) :: ZLE     ! Layer depths [m]

!
! !OUTPUT PARAMETERS:
    real, intent(out), dimension(nc) :: strokeRate ! Flashes per second
    integer, optional, intent(out) :: RC 
    
!
! !INPUT/OUTPUT PARAMETERS:
    TYPE(MAPL_MetaComp), POINTER :: STATE ! Internal MAPL_Generic state

!BOC
!
! !LOCAL VARIABLES:


    REAL, DIMENSION(nc) :: A1X1
    REAL, DIMENSION(nc) :: A2X2
    REAL, DIMENSION(nc) :: A3X3
    REAL, DIMENSION(nc) :: A4X4
    REAL, DIMENSION(nc) :: A5X5

    real :: a0c,a0m         ! Coefficients at continental and marine locations
    REAL :: a1c,a1m
    REAL :: a2c,a2m
    REAL :: a3c,a3m
    REAL :: a4c,a4m
    REAL :: a5c,a5m

    REAL :: x1Divisor       ! Divisors for x1-x5.
    REAL :: x2Divisor
    REAL :: x3Divisor
    REAL :: x4Divisor
    REAL :: x5Divisor

    real :: x5Power ! Exponent for the surface temperature deviation predictor


    REAL :: sfcTLimit       ! Temperature thresholds
    REAL :: airTLimit

    REAL :: hPaCldTop       ! Cloud top limiter for weak/no convection

    REAL, ALLOCATABLE, DIMENSION(:) :: x1         ! Five independent variables
    REAL, ALLOCATABLE, DIMENSION(:) :: x2
    REAL, ALLOCATABLE, DIMENSION(:) :: x3
    REAL, ALLOCATABLE, DIMENSION(:) :: x4
    REAL, ALLOCATABLE, DIMENSION(:) :: x5


    REAL, ALLOCATABLE, DIMENSION(:) :: cloudTopAG ! Cloud top height above ground
    REAL, ALLOCATABLE, DIMENSION(:) :: cnv_topp   ! Convective cloud top pressure with MAPL_UNDEFs
                                                  ! changed to zero
      

    REAL, ALLOCATABLE, DIMENSION(:,:) :: dZ       ! Layer depths [m]
    REAL, ALLOCATABLE, DIMENSION(:,:) :: p        ! Pressure at middle of layer [Pa]
    REAL, ALLOCATABLE, DIMENSION(:,:) :: T        ! Air temperature at middle of layer [K]

    INTEGER, ALLOCATABLE, DIMENSION(:)   :: weakCnvMask   ! Weak or no convection mask
    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: mask          ! Working mask
    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: cloudTopMask  ! Mask is 1 below cloud top


    integer                    :: STATUS, k, i, n
    CHARACTER(LEN=*), PARAMETER :: Iam = "computeFlashRate"





    ! Preliminaries
    ! -------------
    RC = 0
    strokeRate(:) = 0.0
    

    print*, "after strokeRate"

    ! Coefficients of the predictors, marine locations
    ! ------------------------------------------------
    CALL MAPL_GetResource(STATE,a0m,'MARINE_A0:',DEFAULT= 0.0139868,RC=STATUS)
    VERIFY_(STATUS)
    CALL MAPL_GetResource(STATE,a1m,'MARINE_A1:',DEFAULT= 0.0358764,RC=STATUS)
    VERIFY_(STATUS)
    CALL MAPL_GetResource(STATE,a2m,'MARINE_A2:',DEFAULT=-0.0610214,RC=STATUS)
    VERIFY_(STATUS)
    CALL MAPL_GetResource(STATE,a3m,'MARINE_A3:',DEFAULT=-0.0102320,RC=STATUS)
    VERIFY_(STATUS)
    CALL MAPL_GetResource(STATE,a4m,'MARINE_A4:',DEFAULT= 0.0031352,RC=STATUS)
    VERIFY_(STATUS)
    CALL MAPL_GetResource(STATE,a5m,'MARINE_A5:',DEFAULT= 0.0346241,RC=STATUS)
    VERIFY_(STATUS)


    ! Coefficients of the predictors, continental locations
    ! -----------------------------------------------------
    CALL MAPL_GetResource(STATE,a0c,'CONTINENT_A0:',DEFAULT=-0.0183172,RC=STATUS)
    VERIFY_(STATUS)
    CALL MAPL_GetResource(STATE,a1c,'CONTINENT_A1:',DEFAULT=-0.0562338,RC=STATUS)
    VERIFY_(STATUS)
    CALL MAPL_GetResource(STATE,a2c,'CONTINENT_A2:',DEFAULT= 0.1862740,RC=STATUS)
    VERIFY_(STATUS)
    CALL MAPL_GetResource(STATE,a3c,'CONTINENT_A3:',DEFAULT=-0.0023363,RC=STATUS)
    VERIFY_(STATUS)
    CALL MAPL_GetResource(STATE,a4c,'CONTINENT_A4:',DEFAULT=-0.0013838,RC=STATUS)
    VERIFY_(STATUS)
    CALL MAPL_GetResource(STATE,a5c,'CONTINENT_A5:',DEFAULT= 0.0114759,RC=STATUS)
    VERIFY_(STATUS)

    ! Divisors for nondimensionalization of the predictors
    ! ----------------------------------------------------
    CALL MAPL_GetResource(STATE,x1Divisor,'X1_DIVISOR:',DEFAULT=4.36,RC=STATUS)
    VERIFY_(STATUS)
    CALL MAPL_GetResource(STATE,x2Divisor,'X2_DIVISOR:',DEFAULT=9.27,RC=STATUS)
    VERIFY_(STATUS)
    CALL MAPL_GetResource(STATE,x3Divisor,'X3_DIVISOR:',DEFAULT=34.4,RC=STATUS)
    VERIFY_(STATUS)
    CALL MAPL_GetResource(STATE,x4Divisor,'X4_DIVISOR:',DEFAULT=21.4,RC=STATUS)
    VERIFY_(STATUS)
    CALL MAPL_GetResource(STATE,x5Divisor,'X5_DIVISOR:',DEFAULT=14600.,RC=STATUS)
    VERIFY_(STATUS)


    ! Exponent for the surface temperature deviation predictor
    ! --------------------------------------------------------
    CALL MAPL_GetResource(STATE,x5Power,'X5_EXPONENT:',DEFAULT=3.00,RC=STATUS)
    VERIFY_(STATUS)

    ! Threshold temperatures
    ! ----------------------
    CALL MAPL_GetResource(STATE,sfcTLimit,'SFC_T_LIMIT:',DEFAULT=273.0,RC=STATUS)
    VERIFY_(STATUS)
    CALL MAPL_GetResource(STATE,airTLimit,'AIR_T_LIMIT:',DEFAULT=263.0,RC=STATUS)
    VERIFY_(STATUS)
    
    ! Cloud-top pressure limiter
    ! --------------------------
    CALL MAPL_GetResource(STATE,hPaCldTop,'CLOUD_TOP_LIMIT:',DEFAULT=500.,RC=STATUS)
    VERIFY_(STATUS)

    print*, "After Resources"

    ! Layer depths [m]
    ! ----------------
    ALLOCATE(dZ(nc,lm),STAT=STATUS)
    VERIFY_(STATUS)
    dZ = zle(:,0:lm-1)-zle(:,1:lm)


    ! Pressure at mid-layer [Pa]
    ! --------------------------
    ALLOCATE(p(nc,lm),STAT=STATUS)
    VERIFY_(STATUS)
    p = (ple(:,1:lm)+ple(:,0:lm-1))*0.50


    ! Temperature at mid-layer [K]
      ! ----------------------------
    ALLOCATE(T(nc,lm),STAT=STATUS)
    VERIFY_(STATUS)
    T = TH*((p*1.00E-05)**(MAPL_RGAS/MAPL_CP))

    ! Reset CNV_TOPP's MAPL_UNDEFs to zeroes
    ! --------------------------------------
    ALLOCATE(cnv_topp(nc),STAT=STATUS)
    WHERE(CCTP == MAPL_UNDEF)
       cnv_topp = 0.00
    ELSEWHERE
       cnv_topp = CCTP
    END WHERE


    ! Set weak/no convection mask
    ! ---------------------------
    ALLOCATE(weakCnvMask(nc),STAT=STATUS)
    VERIFY_(STATUS)
    weakCnvMask = 0
    WHERE(cn_prcp == 0.00 .OR. cnv_topp >= hPaCldTop*100.00 .OR. CAPE >= MAPL_UNDEF) weakCnvMask = 1

    ! Convective cloud top mask
    ! -------------------------
    ALLOCATE(cloudTopMask(nc,lm),STAT=STATUS)
    VERIFY_(STATUS)
    cloudTopMask = 0
    DO k = 1,lm
       WHERE(ple(1:nc,k) > cnv_topp(1:nc) .AND. cnv_topp(1:nc) > 0.00) cloudTopMask(1:nc,k) = 1
    END DO

    print*, "Before cloud top distance"

    ! Cloud top distance above ground [m]
    ! -----------------------------------
    ALLOCATE(cloudTopAG(nc),STAT=STATUS)
    VERIFY_(STATUS)
    cloudTopAG = 0.00
    DO i = 1,nc
       n = SUM(cloudTopMask(i,1:lm))
       IF(n > 0) cloudTopAG(i) = SUM(dZ(i,lm-n+1:lm))
    END DO


    ! X1: Cold cloud depth: Vertical extent [km] where T < airTLimit and p > cnv_topp
    ! -------------------------------------------------------------------------------
    ALLOCATE(x1(nc),STAT=STATUS)
    VERIFY_(STATUS)
    ALLOCATE(mask(nc,lm),STAT=STATUS)
    VERIFY_(STATUS)

    mask = 0
    WHERE(T < airTLimit .AND. cloudTopMask == 1) mask = 1

    x1 = 0.00
    DO i = 1,nc
       DO k = 1,lm
          IF(mask(i,k) == 1) x1(i) = x1(i)+dZ(i,k)*0.001
       END DO
    END DO
    WHERE(weakCnvMask == 1) x1 = 0.00
    x1 = x1/x1Divisor

    
    print*, "Before X4"

    ! X4: Integrated convective mass flux
    ! -----------------------------------
    ALLOCATE(x4(nc),STAT=STATUS)
    VERIFY_(STATUS)
    x4 = 0.00
    DO i = 1,nc
       DO k = 1,lm
          IF(mask(i,k) == 1) x4(i) = x4(i)+cnv_mfc(i,k)*dZ(i,k)
       END DO
    END DO
    WHERE(weakCnvMask == 1) x4 = 0.00
    x4 = x4/x4Divisor


    ! X5: Surface temperature deviation from sfcTLimit, positive only.
    ! Note: UNDEF TS test retains the ability to boot-strap moist_import_rst.
    ! -----------------------------------------------------------------------
    ALLOCATE(x5(nc),STAT=STATUS)
    VERIFY_(STATUS)
    WHERE(TS == MAPL_UNDEF)
       x5 = 0.00
    ELSEWHERE
       x5 = TS-sfcTLimit
    END WHERE
    WHERE(weakCnvMask == 1) x5 = 0.00
    WHERE(x5 < 0.00) x5 = 0.00
    x5 = x5**x5Power/x5Divisor


    ! X2: Total cloud depth [km]
    ! --------------------------
    ALLOCATE(x2(nc),STAT=STATUS)
    VERIFY_(STATUS)
    x2 = cloudTopAG*0.001
    WHERE(weakCnvMask == 1) x2 = 0.00
    x2 = x2/x2Divisor


    print*, "Before CAPE"

    ! X3: CAPE
    ! --------
    ALLOCATE(x3(nc),STAT=STATUS)
    VERIFY_(STATUS)
    x3 = CAPE
    WHERE(weakCnvMask == 1) x3 = 0.00
    x3 = x3/x3Divisor

    ! Polynomial fit [units: km^{-2} s^{-1}] and individual
    ! terms including marine and continental discrimination
    ! -----------------------------------------------------
    WHERE(frOcean >= 0.01)
       strokeRate = (a0m + a1m*x1 + a2m*x2 + a3m*x3 + a4m*x4 + a5m*x5)/86400.00
       A1X1 = a1m*x1/86400.00
       A2X2 = a2m*x2/86400.00
       A3X3 = a3m*x3/86400.00
       A4X4 = a4m*x4/86400.00
       A5X5 = a5m*x5/86400.00
    ELSEWHERE
       strokeRate = (a0c + a1c*x1 + a2c*x2 + a3c*x3 + a4c*x4 + a5c*x5)/86400.00
       A1X1 = a1c*x1/86400.00
       A2X2 = a2c*x2/86400.00
       A3X3 = a3c*x3/86400.00
       A4X4 = a4c*x4/86400.00
       A5X5 = a5c*x5/86400.00
    END WHERE


    print*, "Before elim negatives"

    ! Eliminate negatives
    ! -------------------
    WHERE(strokeRate < 0.00) strokeRate = 0.00
    
    ! Set rate to zero where any of x1 through x5 are zero
    ! ----------------------------------------------------
    WHERE(x1 == 0.00) strokeRate = 0.00
    WHERE(x2 == 0.00) strokeRate = 0.00
    WHERE(x3 == 0.00) strokeRate = 0.00
    WHERE(x4 == 0.00) strokeRate = 0.00
    WHERE(x5 == 0.00) strokeRate = 0.00

!     print*, "nc = ", nc
!     print*, "min/max TS: ", minval(TS), maxval(TS)
!     print*, "min/max PLE: ", minval(PLE), maxval(PLE)
!     print*, "min/max ZLE: ", minval(ZLE), maxval(ZLE)
!     print*, "min/max CCTP: ", minval(CCTP), maxval(CCTP)
!     print*, "min/max FROCEAN: ", minval(FROCEAN), maxval(FROCEAN)
!     print*, "min/max CN_PRCP: ", minval(CN_PRCP), maxval(CN_PRCP)
!     print*, "a0m-5m: ", a0m, a1m, a2m, a3m, a4m, a5m
!     print*, "a0m-5m: ", a0c, a1c, a2c, a3c, a4c, a5c
!     print*, "x1Divisor-x5: ", x1Divisor, x2Divisor, x3Divisor, x4Divisor, &
!          x5Divisor
!     print*, "Threshold temperatures: ", sfcTLimit, airTLimit
!     print*, "Cloud-top pressure limiter: ", hPaCldTop
!     print*, "Size of dZ: ", shape(dZ)
!     print*, "Layer depths 1, lm: ", dZ(1,1), dZ(1,lm)
!     print*, "Pressure [Pa]  at mid-layer 1, lm: ", p(1,1), p(1,lm)
!     print*, "Temperature [K] at mid-layer 1, lm: ", T(1,1), T(1,lm)
!     print*, "min/max cnv_topp: ", minval(cnv_topp), maxval(cnv_topp)
!     print*, "min/max strokeRate:", minval(strokeRate), maxval(strokeRate)
!     print*, "min/max weakCnvMask: ", minval(weakCnvMask), maxval(weakCnvMask)
!     print*, "min/max cloudTopMask: ", minval(cloudTopMask), maxval(cloudTopMask)
!     print*, "min/max cloudTopAG: ", minval(cloudTopAG), maxval(cloudTopAG)
!     print*, "x1(1): ", x1(1)
!     print*, "x2(1): ", x2(1)
!     print*, "x3(1): ", x3(1)
!     print*, "x4(1): ", x4(1)
!     print*, "x5(1): ", x5(1)
!    print*, "computeFlashRate: A1X1-5(1): ", A1X1(1), A2X2(1), A3X3(1), A4X4(1), A5X5(1)


    DEALLOCATE(x1,STAT=STATUS)
    VERIFY_(STATUS)
    DEALLOCATE(x2,STAT=STATUS)
    VERIFY_(STATUS)
    DEALLOCATE(x3,STAT=STATUS)
    VERIFY_(STATUS)
    DEALLOCATE(x4,STAT=STATUS)
    VERIFY_(STATUS)
    DEALLOCATE(x5,STAT=STATUS)
    VERIFY_(STATUS)
    DEALLOCATE(cnv_topp,STAT=STATUS)
    VERIFY_(STATUS)
    DEALLOCATE(dZ,STAT=STATUS)
    VERIFY_(STATUS)
    DEALLOCATE(p,STAT=STATUS)
    VERIFY_(STATUS)
    DEALLOCATE(T,STAT=STATUS)
    VERIFY_(STATUS)
    DEALLOCATE(cloudTopAG,STAT=STATUS)
    VERIFY_(STATUS)
    DEALLOCATE(mask,STAT=STATUS)
    VERIFY_(STATUS)
    DEALLOCATE(cloudTopMask,STAT=STATUS)
    VERIFY_(STATUS)
    DEALLOCATE(weakCnvMask,STAT=STATUS)
    VERIFY_(STATUS)
    
    return

  end subroutine computeFlashRate

!-------------------------------------------------------------------------
!EOP


!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: flashfit_flipped

!
! !INTERFACE:
!
    subroutine flashfit_flipped (cldmas, threshold, ratio_local, ratio_global, midlatAdj, &
         desired_g_N_prod_rate, flashrate)

! !USES
      implicit none

! !INPUT PARAMETERS:

      !     Convective mass flux at desired layer (kg m-2 min-1)
      !     Warning: Make sure units are correct and that you've accessed the proper vertical layer. 
      real*8, intent(in) :: cldmas(:, :)    
      
      !     Desired mean lightning NO production rate (Tg yr-1) (specified in namelist)       	
      real*8, intent(in) :: desired_g_N_prod_rate
      
      !     Mass flux threshold (kg m-2 min-1) below which flash rate is assumed to be zero (Read in from file) 
      real*8, intent(in) :: threshold            
      
      !     Adjustment factor local flash rates must be multiplied so that monthly average local flash rates (after
      !     "ratio_global" adjustment) match monthly average local v2.2 climatological OTD/LIS flash rates (Read in from file).  
      !     Note: ratio_local varies monthly. 
      real*8, intent(inout) :: ratio_local(:,:) 

      real*8, intent(in) :: midLatAdj(:,:) 
              
      !     Adjustment factor local flash rates must be multiplied so that globally averaged 
      !     flash rate matches v2.2 OTD/LIS climatological globally averaged flash rate (Read in from file) 
      !     Note: ratio_global varies monthly. 
      real*8, intent(in) :: ratio_global             
            
      
! !OUTPUT PARAMETERS:

      !     Total (CG+IC) flash rate (flashes per grid box per s)
      real*8, intent(inout) :: flashrate(:,:)

! ! DESCRIPTION:
!EOP
!-------------------------------------------------------------------------


      ! !Local variables
      real*8, allocatable :: cldmas_local(:,:)
      integer :: shapeCldmas(2)
      integer :: RC
      
      shapeCldmas(:) = shape(cldmas)
      
      allocate(cldmas_local(shapeCldmas(1),shapeCldmas(2)), STAT=RC)

      ! If deep convection exists
      cldmas_local(:,:) = cldmas(:,:)
      where(cldmas_local > 0.025) cldmas_local = 0.025 
          
      where(cldmas_local > threshold) ! comparision is made in s-1 
         ! calculation is made in min-1
         flashrate = ratio_global*(cldmas_local(:,:)*60.-threshold*60.)*ratio_local(:,:)/60.
      elsewhere
         flashrate = 0.  
      end where

      ! Calculate N production rate (g s-1).
      ! For specified flashrate and PRODFAC, global production rate = 3.41 Tg N yr-1.
      ! Assume IC and CG flash rates produce the same amount of N / flash.
      !      pnox = flashrate(:,:)*PRODFAC*CONVFAC
      ! pnox = midLatAdj(:,:)*flashrate(:,:)*PRODFAC*CONVFAC 
      
      !adjust for desired global Nitrogen production rate
      !pnox = pnox * desired_g_N_prod_rate /3.41        
    
      deallocate(cldmas_local)
      return

    end subroutine flashfit_flipped
!EOP

!==========================================================================

!BOP

! !IROUTINE: hemcoFlashrate

! !DESCRIPTION: Wrapper routine for GEOS Chem's flashrate routine.
!
! !INTERFACE:
!
    subroutine hemcoFlashrate (cellArea, lwi, lonslocal, latslocal, airTemp, ple, geoPotHeight, &
         cnvMfc, otdLisScale, flashRate, RC )

! !USES
      implicit none

! !ARGUMENTS

      !     Grid box cell area (m2)
      real, intent(in) :: cellArea(:,:)    
      
      !     Land water ice flag 
      !     need to translate to (0=land; 1=water; 2=ice)
      real, intent(in) :: lwi(:,:)

      ! Air temperature in K
      real, intent(in) :: airTemp(:,:,:)

      ! Pressure edges (Pa)
      real, intent(in) :: ple(:,:,:)

      !     Geo potential height (m)
      real, intent(in) :: geoPotHeight(:,:,:)    

      ! Cumulative mass flux (kg/m2/s)
      real, intent(in) :: cnvMfc(:,:,:)

!      ! Buoyancy (not sure about units. coming from computeCAPE)
!      real, intent(in) :: buoyancy(:,:,:)

      ! otdLisScale = 0.355 for c48; 0.1 for c90
      real, intent(in) :: otdLisScale

      ! Flash rate in flashes km-2 s-1
      real, intent(out) :: flashRate(:,:) 

      ! Longitude and latitude (local values) in degrees
      real, intent(in) :: lonslocal(:,:)
      real, intent(in) :: latslocal(:,:)

      ! return code
      integer, intent(inout) :: RC

!EOP

! Locals


      
      real(hp) :: otdLisScaleHp

      real(hp), pointer, dimension(:,:)   ::  cellAreaHemco => null()
      real(hp), pointer, dimension(:,:)   ::  otdLisHemco => null()
      real(hp), pointer, dimension(:,:)   ::  lonsHemco => null()
      real(hp), pointer, dimension(:,:)   ::  latsHemco => null()
      real(hp), pointer, dimension(:,:)   ::  convFracHemco => null()
      real(hp), pointer, dimension(:,:,:)   ::  airTempHemco => null()
      real(hp), pointer, dimension(:,:,:)   ::  pleHemco => null()
      real(hp), pointer, dimension(:,:,:)   ::  gridBoxHeightHemco => null()
      real(hp), pointer, dimension(:,:,:)   ::  gridBoxHeightHemco2 => null()
      real(hp), pointer, dimension(:,:,:)   ::  cnvMfcHemco => null()
      real(hp), pointer, dimension(:,:,:)   ::  buoyancyHemco => null()

      integer,  pointer, dimension(:,:)  ::   lwiHemco => null()
      
      ! outputs
      real(hp), pointer, dimension(:,:) :: lNox2d => null() ! Total lightning NOx molecules per 6 hrs
      real(hp), pointer, dimension(:,:) :: lfr2d => null()  ! flashrate (flashes/min/km2)
      real(hp), pointer, dimension(:,:) :: ic2d => null()   ! Intra-cloud flash rate (flashes/min/km2)
      real(hp), pointer, dimension(:,:) :: cg2d => null()   ! Cloud-to-ground flash rate (flashes/min/km2)
      real(hp), pointer, dimension(:,:) :: h02d => null()   ! Convective cloud top height (m)
      integer,  pointer, dimension(:,:) :: ltop2d => null() ! Lightning top level


      integer :: shape3d (3)
      integer :: IM, JM, LM, KM
      integer :: levCount, kR

      real, parameter :: MIN_TO_SEC = 1.0/60.0

!EOP
!-------------------------------------------------------------------------


      shape3d = shape(airTemp)
      IM = shape3d(1)
      JM = shape3d(2)
      LM = shape3d(3)
      KM = LM + 1

      ALLOCATE(pleHemco     (1:IM, 1:JM, 1:LM+1),   STAT=RC)
      ALLOCATE(cnvMfcHemco  (1:IM, 1:JM, 1:LM+1),   STAT=RC)
      ALLOCATE(airTempHemco (1:IM, 1:JM, 1:LM),   STAT=RC)
      ALLOCATE(buoyancyHemco(1:IM, 1:JM, 1:LM),   STAT=RC)
      ALLOCATE(gridBoxHeightHemco (1:IM, 1:JM, 1:LM),   STAT=RC)
      ALLOCATE(gridBoxHeightHemco2 (1:IM, 1:JM, 1:LM),   STAT=RC)
      ALLOCATE(cellAreaHemco(1:IM, 1:JM),   STAT=RC)
      ALLOCATE(convFracHemco(1:IM, 1:JM),   STAT=RC)
      ALLOCATE(lwiHemco     (1:IM, 1:JM),   STAT=RC)
      ALLOCATE(otdLisHemco  (1:IM, 1:JM),   STAT=RC)
      ALLOCATE(lonsHemco    (1:IM, 1:JM),   STAT=RC)
      ALLOCATE(latsHemco    (1:IM, 1:JM),   STAT=RC)


      ALLOCATE(lNox2d    (1:IM, 1:JM),   STAT=RC)
      ALLOCATE(lfr2d     (1:IM, 1:JM),   STAT=RC)
      ALLOCATE(ic2d   (1:IM, 1:JM),   STAT=RC)
      ALLOCATE(cg2d   (1:IM, 1:JM),   STAT=RC)
      ALLOCATE(h02d   (1:IM, 1:JM),   STAT=RC)
      ALLOCATE(ltop2d (1:IM, 1:JM),   STAT=RC)

      ! edge layer fields
      pleHemco(:,:,1:KM) = ple(:,:,KM:1:-1)
      cnvMfcHemco(:,:,1:KM) = cnvMfc(:,:,KM:1:-1)

      ! fields at grid box center
      airTempHemco(:,:,1:LM) = airTemp(:,:,LM:1:-1)

      ! buoyancy is only needed if convFrac is provided and if the corresponding toggle is 
      ! turned on in CalcFlashRate below - but it's hardcoded to false.
      buoyancyHemco(:,:,1:LM) = real(-1) !buoyancy(:,:,LM:1:-1)


      do levCount=1,LM
         gridBoxHeightHemco(:,:,levCount) = geoPotHeight(:,:,levCount) - geoPotHeight(:,:,levCount+1)
      enddo
      gridBoxHeightHemco2(:,:,1:LM) = gridBoxHeightHemco(:,:,LM:1:-1)


      cellAreaHemco(:,:) = cellArea(:,:)
      otdLisHemco(:,:) = real(-1) ! per instruction from Christoph Keller (C.K.)
      convFracHemco(:,:) = real(-1) ! dummy value per C.K.
      lonsHemco(:,:) = lonslocal(:,:)
      latsHemco(:,:) = latslocal(:,:)

      lNox2d(:,:) = real(0)
      lfr2d(:,:) = real(0)
      ic2d(:,:) = real(0)
      cg2d(:,:) = real(0)
      h02d(:,:) = real(0)
      ltop2d(:,:) = real(0)


      where (lwi == real(0)) 
         lwiHemco = 1
      elsewhere (lwi == real(1))
         lwiHemco = 0
      elsewhere
         lwiHemco = 2
      end where

      otdLisScaleHp = otdLisScale

      call CalcFlashRate (MAPL_AM_I_ROOT(), IM, JM, LM, cellAreaHemco, lwiHemco, &
           otdLisHemco, lonsHemco, latsHemco, airTempHemco, pleHemco, gridBoxHeightHemco2, cnvMfcHemco, &
           buoyancyHemco, convFracHemco, .False., otdLisScaleHp, lNox2d, lfr2d, ic2d, cg2d, ltop2d, h02d, RC)

      ! Convert to km-2 s-1
      IF ( RC==0 ) THEN
         flashRate (:,:) = real(lfr2d(:,:)) * MIN_TO_SEC
      ELSE
         flashRate (:,:) = real(0.0)
      ENDIF

      DEALLOCATE(lwiHemco)
      DEALLOCATE(cellAreaHemco)
      DEALLOCATE(convFracHemco)
      DEALLOCATE(otdLisHemco)
      DEALLOCATE(airTempHemco)
      DEALLOCATE(pleHemco)
      DEALLOCATE(gridBoxHeightHemco)
      DEALLOCATE(gridBoxHeightHemco2)
      DEALLOCATE(cnvMfcHemco)
      DEALLOCATE(buoyancyHemco)
      DEALLOCATE(lonsHemco)
      DEALLOCATE(latsHemco)

      DEALLOCATE(lNox2d)
      DEALLOCATE(lfr2d)
      DEALLOCATE(ic2d)
      DEALLOCATE(cg2d)
      DEALLOCATE(h02d)
      DEALLOCATE(ltop2d)

      return

    end subroutine hemcoFlashrate


!-----------------------------------------------------------------------
!BOP
! !IROUTINE: flash_rate_Lopez
!
! !DESCRIPTION:
! Lightning parameterization based on 
! "A Lightning Parameterization for the ECMWF Integrated Forecasting System"
! P. Lopez 2016 MWR
! !REVISION HISTORY:!
!  - coded/adapted to the GF scheme by Saulo Freitas (10-Aug-2019)
!  - code brought to this module by Megan Damon 
!=====================================================================================

!EOP
!-----------------------------------------------------------------------


    subroutine flash_rate_Lopez(IM, JM, LM, FRLAND, ZKBCON, CAPE,  &
                               ZLE, PFI_CN_GF, CNV_QC, TE, PLO,  &
                               LFR_Lopez)
     implicit none 
     integer ,intent(in)  :: im,jm,lm
     real    ,intent(in), dimension(im,jm,0:lm) ::  &
                              ZLE        & ! layer depths [m]
                             ,PFI_CN_GF       ! 3d_ice_precipitation_flux_GF [kg/m2/s]

     real    ,intent(in), dimension(im,jm,lm) ::  &
                             TE             &  ! air temp [K]
                             ,PLO           &  ! press (hPa)
                             ,CNV_QC        ! grid_mean_convective_condensate [kg/kg]


     real    ,intent(in), dimension(im,jm) :: FRLAND    & ! fraction of land
          , cape   &  ! Convective available potential energy [J/kg]
          , zkbcon    ! cloud_base_height_deep_GF [m]



     real    ,intent(out), dimension(im,jm) ::  &
                              LFR_Lopez           ! lightning flash density rate (units: 1/km2/day)
                         
    !-- locals
    real, parameter :: V_graup     = 3.0  ! m/s
    real, parameter :: V_snow      = 0.5  ! m/s
    real, parameter :: beta_land   = 0.70 ! 1
    real, parameter :: beta_ocean  = 0.45 ! 1
    real, parameter :: alpha       = 37.5 ! 1
    real, parameter :: t_initial   =  0.0 + 273.15 ! K
    real, parameter :: t_final     = -25. + 273.15 ! K
    
    integer :: i,j,k, k_initial, k_final                          
    real    :: Q_R, z_base,beta,prec_flx_fr,dz
    real,    dimension(1:lm) :: q_graup,q_snow,rho



     print*, ""
     print*, "Hello from flash_rate_Lopez"
     print*, ""

     print*, "im, jm, lm: ", im, jm, lm
     print*, "FRLAND: ", minval(FRLAND), maxval(FRLAND)
     print*, "CAPE: ", minval(CAPE), maxval(CAPE)
     print*, "ZLE: ", minval(ZLE), maxval(ZLE)
     print*, "PFI_CN_GF ", minval(PFI_CN_GF), maxval(PFI_CN_GF)
     print*, "TE: ", minval(TE), maxval(TE)
     print*, "PLO: ", minval(PLO), maxval(PLO)
     print*, "CNV_QC: ", minval(CNV_QC), maxval(CNV_QC)
     print*, "zkbcon: ", minval(zkbcon), maxval(zkbcon)

    DO j = 1, jm
        DO i = 1, im
           LFR_Lopez(i,j) = 0.0
           
           if(ZKBCON(i,j) <= 0. .or. CAPE(i,j) == MAPL_UNDEF) cycle !-> no convection 

           beta= frland(i,j)*beta_land + (1.-frland(i,j))*beta_ocean
           q_graup(:) = 0. 
           q_snow (:) = 0. 

           do k=lm,1,-1
                rho(k) =  100.*PLO(i,j,k) / (MAPL_RGAS*TE(i,j,k) )
                prec_flx_fr = PFI_CN_GF(i,j,k) / rho(k)

                q_graup(k) =      beta *prec_flx_fr/V_graup ! - graupel mixing ratio (kg/kg)
                q_snow (k) =  (1.-beta)*prec_flx_fr/V_snow  ! - snow    mixing ratio (kg/kg)

          enddo
          k_initial = minloc(abs(te(i,j,:)-t_initial),1)  
          k_final   = minloc(abs(te(i,j,:)-t_final  ),1)  

          Q_R = 0.0
          do k = k_final, k_initial
             
             dz  = - (zle(i,j,k)-zle(i,j,k-1))  

             Q_R = Q_R + dz*rho(k)*(q_graup(k)*(CNV_QC(i,j,k)+q_snow(k)))
         enddo

         z_base = ZKBCON(i,j)/1000. !- convert to [km]

        !--- lightning flash density (units: number of flashes/km2/day) - equation 5
        !--- (to compare with Lopez 2016's results, convert to per year: LFR_Lopez*365)
        !
        LFR_Lopez(i,j) = alpha * Q_R *sqrt(max(0.,cape(i,j))) * min(z_base,1.8)**2 

     enddo
   enddo
        
        
if (maxval(LFR_Lopez)>1000.) then
       print*,"1================================================= LOPEZ"
       print*," LFR=",maxval(LFR_Lopez),maxval(cape),maxval(CNV_QC),maxval(PFI_CN_GF)
       print*,"2================================================= LOPEZ"
       call flush(6)
endif        
  end subroutine flash_rate_Lopez
!========================================================================================


!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: flashfit
! NOTE: The flashfit routine was written for GMI-CTM vertical orientation
! The surface is level 1, to the top of atomsphere
! This routine will flip GEOS-5 oriented arrays and call the 
! GMI-CTM routine
!
! !INTERFACE:

    subroutine flashfit ()

      print*, "I will call flashfit_flipped after I flip GEOS-5 arrays"

    end subroutine flashfit

!-------------------------------------------------------------------------
!BOP
! ! ROUTINE: computeCAPE
!
! !INTERFACE: 
!
  subroutine computeCAPE (TH, Q, PLE, CAPE, BUOY, IM, JM, LM)
    
! !DESCRIPTION: 
!
! !USES:
    use GEOS_UtilsMod    
    implicit none


! !INPUT PARMETERS:

      integer,                     intent(in)  :: IM,JM,LM
      real, dimension(IM,JM,LM),   intent(in)  :: TH  ! potential temperature
      real, dimension(IM,JM,LM),   intent(in)  :: Q   ! specific humidity
      real, dimension(IM,JM,0:LM), intent(in)  :: PLE   ! pressure

! !OUTPUT PARAMETERS:
      real, dimension(IM,JM),      intent(out) :: CAPE
      real, dimension(IM,JM,LM),   intent(out) :: BUOY

! !DESCRIPTION:

! !REVISION HISTORY:
!
! 17Nov18 - Megan Damon.
!
!EOP
!-------------------------------------------------------------------------

! !Local Variables
      integer                         :: L
      real,    dimension(IM,JM,  LM)  :: DQS, QSS, PLO, TEMP, PK, DM, DP
      real,    dimension(IM,JM,  LM)  :: ZLO
      real,    dimension(IM,JM,0:LM)  :: ZLE
      real,    dimension(IM,JM,0:LM)    :: CNV_PLE
      real,    dimension(IM,JM,0:LM)    :: PKE
      real,    dimension(IM,JM   )  :: HC
      logical, dimension(IM,JM   )  :: UNSTABLE

! Initialize local variables
      CNV_PLE  = PLE*.01
      PLO      = 0.5*(CNV_PLE(:,:,0:LM-1) +  CNV_PLE(:,:,1:LM  ) )
      PKE      = (CNV_PLE/1000.)**(MAPL_RGAS/MAPL_CP)
      DP       = ( PLE(:,:,1:LM)-PLE(:,:,0:LM-1) )
      PK       = (PLO/1000.)**(MAPL_RGAS/MAPL_CP)
      DM       = DP*(1./MAPL_GRAV)
      TEMP     = TH*PK
      DQS      = GEOS_DQSAT(TEMP, PLO, qsat=QSS)

      ZLE(:,:,LM) = 0.
      do L=LM,1,-1
         ZLE(:,:,L-1) = TH (:,:,L) * (1.+MAPL_VIREPS*Q(:,:,L))
         ZLO(:,:,L  ) = ZLE(:,:,L) + (MAPL_CP/MAPL_GRAV)*( PKE(:,:,L)-PK (:,:,L  ) ) * ZLE(:,:,L-1)
         ZLE(:,:,L-1) = ZLO(:,:,L) + (MAPL_CP/MAPL_GRAV)*( PK (:,:,L)-PKE(:,:,L-1) ) * ZLE(:,:,L-1)
      end do

      ! From BUYOANCY

       HC  =  TEMP(:,:,LM) + (MAPL_GRAV/MAPL_CP)*ZLO(:,:,LM) + (MAPL_ALHL/MAPL_CP)*Q(:,:,LM)

       do L=LM-1,1,-1
          BUOY(:,:,L) = HC - (TEMP(:,:,L) + (MAPL_GRAV/MAPL_CP)*ZLO(:,:,L) + (MAPL_ALHL/MAPL_CP)*QSS(:,:,L))
          BUOY(:,:,L) = BUOY(:,:,L) / ( (1.+ (MAPL_ALHL/MAPL_CP)*DQS(:,:,L))*TEMP(:,:,L) )
       enddo

       BUOY(:,:,LM) = 0.0

       UNSTABLE = .false.

       CAPE = 0.

! New formulation
       do L=1,LM-1
         where(BUOY(:,:,L)>0.)
            CAPE = CAPE + BUOY(:,:,L)*DM(:,:,L)
         end where
      end do
! Old formulation
!       do L=1,LM-1
!          where(BUOY(:,:,L)>0.) UNSTABLE=.true.
!          where(UNSTABLE)
!             CAPE = CAPE + BUOY(:,:,L)*DM(:,:,L)
!          end where
!       end do

       UNSTABLE = CAPE > 0.0

       where(.not.UNSTABLE)
          CAPE=MAPL_UNDEF
       end where

      return

      end subroutine computeCAPE
!EOC

!-----------------------------------------------------------------------------
! !ROUTINE
!   emiss_lightning
!
! !DESCRIPTION
!   
!  Generate NOx production rates from parameterized lightning and distribute vertically
!  using profiles from Pickering (2007).
!  This routine came from GmiEmiss_lightning_mod in the GEOS-CTM
!
! !REVISION HISTORY:
!   July 28, 2008 - Dale Allen.  First version, obsolete.  
!   December 30, 2011 - Eric Nielsen: Simplified for GEOS-5 with flash rates imported from MOIST.
!                       Target: ~4 Tg N yr^{-1}. Best advice: ~250 moles NO per flash, which 
!                       translates to numberNOperFlash = 1.5E+26.
!   January 24, 2013 - Eric Nielsen: kgNOx3D changed from N to NO production rate for export 
!                       EM_LGTNO. Purpose is to allow direct comparison with GMI NO_lgt.
!-----------------------------------------------------------------------------
 SUBROUTINE emiss_lightning (i1, i2, j1, j2, k1, k2, minDeepCloudTop, ampFactor, numberNOperFlash, &
                            lwi, flashrate, cellDepth, dtrn, pNOx3D, kgNOx3D, rc)

  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: i1, i2, j1, j2, k1, k2         ! Index ranges on this processor
  REAL,    INTENT(IN)  :: minDeepCloudTop                ! Minimum cloud top [km] for selecting deep convection profiles
  REAL,    INTENT(IN)  :: ampFactor                      ! > 0, for targeting the observed nitrogen production rate [3.41 Tg yr^{-1}]
  REAL,    INTENT(IN)  :: numberNOperFlash               ! NO molecules generated by each flash
  INTEGER, INTENT(IN)  :: lwi(i1:i2, j1:j2)              ! Flag: 1=water 2=land 3=ice
  REAL*8,  INTENT(IN)  :: flashrate(i1:i2, j1:j2)        ! Flash rate [km^{-2} s^{-1}]
  REAL*8,  INTENT(IN)  :: dtrn(i1:i2, j1:j2, k1:k2)      ! Detrainment [kg m^{-2} s^{-1}]
  REAL*8,  INTENT(IN)  :: cellDepth(i1:i2, j1:j2, k1:k2) ! Grid cell depth [m]
  
  REAL*8, INTENT(OUT)  :: pNOx3D(i1:i2, j1:j2, k1:k2)    ! Lightning NO production rate [m^{-3} s^{-1}]
  REAL*8, INTENT(OUT)  :: kgNOx3D(i1:i2, j1:j2, k1:k2)   ! NO production rate [kg m^{-3} s^{-1}]

! Local
! -----
  INTEGER :: k
  INTEGER :: status, rc
  REAL*8, ALLOCATABLE  :: pNOx2D(:,:)                    ! Lightning NO production [molecules NO m^{-2} s^{-1}]
  CHARACTER(LEN=*), PARAMETER :: Iam = "emiss_lightning"
  rc = 0
  status = 0


! Validate ranges for ampFactor and numberNOperFlash
! --------------------------------------------------
  IF(ampFactor <= 0.00) THEN
   IF(MAPL_AM_I_ROOT()) PRINT *,TRIM(IAm)//": Invalid ampFactor ",ampFactor
   status = 1
   VERIFY_(status)
  END IF
  IF(numberNOperFlash <= 0.00) THEN
   IF(MAPL_AM_I_ROOT()) PRINT *,TRIM(IAm)//": Invalid numberNOperFlash ",numberNOperFlash
   status = 1
   VERIFY_(status)
  END IF

! Grab some memory
! ----------------
  ALLOCATE(pNOx2D(i1:i2, j1:j2),STAT=status)
  VERIFY_(status)
  pNOx2D(:,:) = 0.00

! Calculate the NOx produdction rate [molecules NO m^{-2} s^{-1}]
! ---------------------------------------------------------------
  pNOx2D(:,:) = 1.00E-06*flashrate(:,:)*numberNOPerFlash

! Amplification/suppression factor: > 0
! -------------------------------------
  pNOx2D(:,:) = pNOx2D(:,:)*ampFactor
     
! Partition vertically without changing units
! -------------------------------------------

  CALL partition(i1, i2, j1, j2, k1, k2, pNOx2D,dtrn,cellDepth,minDeepCloudTop,lwi,pNOx3D)

! Place output in useful units
! ----------------------------
  DO k = k1,k2

! Number density tendency [m^{-3} s^{-1}]
! ---------------------------------------
   pNOx3D(i1:i2,j1:j2,k) = pNOx3D(i1:i2,j1:j2,k)/cellDepth(i1:i2,j1:j2,k)

! NO density tendency [kg N m^{-3} s^{-1}]
! ----------------------------------------
   kgNOx3D(i1:i2,j1:j2,k) = pNOx3D(i1:i2,j1:j2,k)*30.0064/(1000.00*AVOGAD)

  END DO

! Clean up
! --------
  DEALLOCATE(pNOx2D,STAT=status)
  VERIFY_(status)

  RETURN
END SUBROUTINE emiss_lightning





!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: partition
! This routine was "parition" in GmiEmissionLightning_mod in GEOS-CTM
!
! INTERFACE:

 subroutine partition (i1,i2,j1,j2,k1,k2,pNOx2D,dtrn,cellDepth,minDeepCloudTop,lwi,pNOx3D)

! !INPUT PARAMETERS:
   integer :: i1, i2, j1, j2, k1, k2
   real*8, intent(in) :: pNOx2D(:,:) 
   real*8, intent(in) :: dtrn(:,:,:) ! Detrainment [kg m^{-2}s^{-1}]
   real*8, intent(in)  :: cellDepth(:,:,:) ! Grid cell depth [m]
   real, intent(in)  :: minDeepCloudTop  ! Minimum cloud top [km] for selecting deep convection profiles
   integer, intent(in) :: lwi(:,:)

! !OUTPUT PARAMETERS:   
   REAL*8, intent(out) :: pNOx3D(i1:i2, j1:j2, k1:k2) ! Scaled production rate (no units conversion here)

! !Local Variables
   character(len=*), parameter :: Iam = "partition"
   
   integer, parameter :: numKm = 17  ! Number of elements (kilometers) in each specified profile

   real, allocatable :: r(:,:)       ! Specified NOx distribution profiles
   real, allocatable :: w(:)         ! Weights applied to scaled cloud layers
   real, allocatable :: z(:)         ! Layer edge heights above ground
   real, allocatable :: zScaled(:)   ! Scaled layer edge heights

   integer :: status, i, j, k, cl
   integer :: nTop, profileNumber
   real :: zLower, zUpper
   
   status = 0

   !print*, "In subroutine ", trim(Iam)

   ! Specify the percentage NOx distributions in each km for a numKm-depth cloud.
   ! Deep convection is arbitrarily assigned when the cloud top is greater than 7 km.
   ! --------------------------------------------------------------------------------
   ALLOCATE(r(numKm,3))


   ! Deep convection, continental
   ! ----------------------------
   r(1:numKm,1) = (/ 0.23, 0.47, 0.56, 1.40, 2.70, 4.00, 5.03, 6.24, &
        8.60,10.28,11.62,12.34,12.70,12.34, 7.63, 3.02, 0.84 /)

   ! Deep convection, marine
   ! -----------------------
   r(1:numKm,2) = (/ 0.60, 1.50, 2.90, 4.30, 5.40, 6.70, 7.70, 8.50, &
        9.60,10.20,10.50,10.20, 8.20, 6.50, 4.50, 2.20, 0.50 /)

   ! Other
   ! -----
   r(1:numKm,3) = (/ 2.40, 5.00, 7.40, 9.30,10.60,11.40,11.50,11.00, &
        9.90, 8.30, 6.30, 4.20, 2.20, 0.50, 0.00, 0.00, 0.00 /)

   ALLOCATE(z(k1:k2))
   ALLOCATE(zScaled(k1:k2))


   ! Work in each column
   ! -------------------
   DO j = j1,j2
      DO i = i1,i2
         
         SeeingLightning : IF(pNOx2D(i,j) > 0.00) THEN

            
            ! Define cloud top to be highest layer with dtrn > 0, but at least 2.
            ! -------------------------------------------------------------------
            DO k = k2,1,-1
               IF(dtrn(i,j,k) > 0) EXIT
            END DO
            nTop = k
            IF(nTop < 2) nTop = 2

            ! Sum grid box thicknesses (m) to obtain layer edge heights
            ! ---------------------------------------------------------
            DO k = 1,k2
               z(k) = SUM(cellDepth(i,j,1:k))
            END DO
            
            ! Select NOx distribution profile. LWI flag is: 1=water 2=land 3=ice
            ! ------------------------------------------------------------------
            IF(z(nTop) > minDeepCloudTop*1000.00) THEN
               IF(lwi(i,j) == 2) THEN
                  profileNumber = 1
               ELSE
                  profileNumber = 2
               END IF
            ELSE
               profileNumber = 3
            ENDIF
            
            ! Scale factor, 0 at ground, numKm at cloud top
            ! ---------------------------------------------
            zScaled(1:k2) = z(1:k2)*numKm/z(nTop)
            
            ! Grab a little work space and intialize
            ! --------------------------------------
            ALLOCATE(w(nTop))
            w(1:nTop) = 0.00
            cl = 1
            zLower = 0.00

            !print*, "Seeing lightning: ", i, j, pNOx2D(i,j), nTop
            !print*, "layer edge heights: ", z(:)
            !print*, "profileNumber: ", profileNumber
            !print*, "zscaled: ", zScaled(:)
            !print*, "cl: ", cl


            ! Work through each km in the specified distribution
            ! --------------------------------------------------
            Kilometers: DO k = 1,numKm
               
               ! ... segment-by-segment
               ! ----------------------
               Segment: DO
                  
                  ! Push up to the lesser of scaled cloud height or next km
                  ! -------------------------------------------------------
                  zUpper = MIN(zScaled(cl),k*1.)
                  IF(zScaled(cl) > numKm) EXIT
                  
                  
                  ! Add increment to scaled weighting for the current cloud layer
                  ! -------------------------------------------------------------
                  w(cl) = w(cl) + (zUpper-zLower) * r(k,profileNumber) * 0.01
                  
                  ! Advance to next cloud layer if any of it lies within this km
                  ! ------------------------------------------------------------
                  IF(zUpper == zScaled(cl)) cl = cl+1
                  
                  ! Shift bounds before working on the next segment
                  ! -----------------------------------------------
                  zLower = zUpper
                  
                  ! At top of this km. Advance to the next one
                  ! ------------------------------------------
                  IF(zUpper == k) EXIT
                  
                  
               END DO Segment

            END DO Kilometers


            ! Finalize vertical distribution and clean up
            ! -------------------------------------------
            pNOx3D(i,j,1:nTop) = w(1:nTop)*pNOx2D(i,j)
            
            
            DEALLOCATE(w)

         endif SeeingLightning
           
      enddo
   enddo


   DEALLOCATE(r)
   DEALLOCATE(z)
   DEALLOCATE(zScaled)
   
   
 end subroutine partition

!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: partitionnox
!
! !INTERFACE:

 subroutine partitionnox (lprslay, imonth, dtrn, mass, pnox2d, pnox3d)
! !NOTE: This routine is NOT fully implemented.
   
! !USES
   implicit none

! !INPUT PARAMETERS:
   real*8, allocatable, intent(in) :: lprslay(:,:,:)
   integer, intent(in)  :: imonth
   real*8,  intent(in)  :: dtrn(:,:,:)
   real*8,  intent(in)  :: mass(:,:,:)
   !      integer, intent(in)  :: cmi_flags(i1:i2, ju1:j2)
   real*8,  intent(in)  :: pnox2d(:,:)

! !OUTPUT PARAMETERS:
   real*8,  intent(inout)  :: pnox3d(:,:,:)


! ! DESCRIPTION:
!     Modeling studies of Pickering (2007) are used as a basis
!      to partition lightning NO in the vertical.
!     Update 4/4/08 (cont. tropical partitioning changed) 
!     Update 6/16/09 (Midlat model cloud profile shifted upward) 
!     INPUTS
!     lprslay:  3-d: Pressure at grid box edges (hPa)
!     imonth: month integer
!     dtrn:     3-d: Detrainment rate (kg m-2 s-1) (units not important)
!     dlatcen:  2-d: Latitude at grid box centers (degrees N)
!     pnox2d:   2-d: lightning NOx production rate in column (g N/s)
!     mass:     3-d: Background (atmospheric) mass in each grid volume (kg)
!     OUTPUTS
!     pnox3d:  lightning NOx production rate in each grid volume (ppv/s)


   ! need to derive/replace: 
   !     cmi_flags 2-d: continental(2), marine(1), or ice flag.



   ! initial arguments
   !     i1,i2, ju1,j2, ilo,ihi, julo,jhi, k1,k2, &
   !     & dlatcen, lprslay, mass, dtrn, cmi_flags,  &
  !     & pnox2d, pnox3d,imonth,flashrate, i2_gl, j2_gl)



! !Local Variables:
   integer :: il,ij,ik,iktop        ! old indices / counters
   integer :: IM, JM, LME, LM       ! new indices using shape array information

   real*8 :: ppmtoppv                  ! conversion from parts per mass to parts per vol.
   
   real*8, allocatable :: ztop(:,:,:)  ! cloud layer tops
   real*8, allocatable :: htedge2(:)   ! Top edge heights after adjusting to match model-cloud
   real*8, allocatable :: fd0(:)       ! Used to calc cont of each "model-cloud"
   real*8, allocatable :: yout(:)      ! Layer to each CTM cloud layer

   integer, allocatable :: inox1(:,:)  ! "Model-cloud" type indices
   integer, allocatable :: ntop(:,:)   ! Top cloud layer indices
   integer :: ikmm                     ! Number of vertical layers in model
   
   integer,parameter :: pmax(4) = (/ 17, 17, 16, 15 /)
   integer ::  shapeArray (3)
   real*8 :: r0(17,4)
   
   ikmm = 40 !? What is this, and does it work for layers other than GEOS (72)?
  
   !     r0(*,i): % of lightning NOy mass dep. into each layer for cld type i
   !     "Model cloud" has 17 1-km layers. Percent of total NOx mass
   !     deposited into each layer is given for 1) tropical continental,  
   !     2) tropical marine, 3) subtropical, and 4) midlatitude "model clouds". 

   data r0/  0.23, 0.47, 0.56, 1.40, 2.70, 4.00, 5.03, 6.24, &
        &          8.60,10.28,11.62,12.34,12.70,12.34, 7.63, 3.02, 0.84, &        
        &          0.60, 1.50, 2.90, 4.30, 5.40, 6.70, 7.70, 8.50, &
        &          9.60,10.20,10.50,10.20, 8.20, 6.50, 4.50, 2.20, 0.50, & 
        &          1.00, 2.10, 3.90, 5.80, 7.70, 9.30,10.50,11.00, &
        &         11.00,10.40, 9.20, 7.50, 5.50, 3.40, 1.50, 0.20, 0.00, &
        &          1.00, 2.37, 4.95, 7.32, 9.21,10.49,11.29,11.39, &
        &         10.89, 9.80, 8.22, 6.24, 4.16, 2.18, 0.49, 0.00, 0.00/ 

   shapeArray = shape(lprslay)
   IM = shapeArray(1)
   JM = shapeArray(2)
   LME = shapeArray(3)
   LM = LME-1

   ALLOCATE (ztop(1:IM,1:JM,1:LME-1))

   !     Assuming a scale height of 8 km, estimate heights at layer tops.
   do ik=1,LME-1
      ztop(:,:,ik) = -8. * dlog(lprslay(:,:,ik)/lprslay(:,:,0))
   enddo

   ! Old way : 
   !     Assuming a scale height of 8 km, estimate heights at layer tops.
   !      do ij = ju1,j2
   !         do il = i1,i2
   !            do ik= k1,k2
   !            ztop(il,ij,ik) = -8.* dlog(lprslay(il,ij,ik)/lprslay(il,ij,0))
   !            !ztop(il,ij,ik) = -8.* dlog(lprslay(il,ij,ik)/pctm2(il,ij))
   !            end do
   !         end do
   !      end do
   
   allocate(inox1(1:IM,1:JM))

   ! !     Determine which "model cloud" to use at each grid point. Decision
   ! !     is based on whether pts are tropical continental, tropical marine, 
   ! !     subtropical, or midlatitude. 
   !       select case (imonth)
   !       case (1,2,3,12)
   ! !     Southern Hemisphere Summer      
   !          where     ((abs(dlatcen) <= 15.).and.(cmi_flags == 2))
   !             inox1 = 1      ! Tropical continental
   !          elsewhere     ((abs(dlatcen) <= 15.).and.(cmi_flags /= 2))
   !             inox1 = 2      ! Tropical marine
   !          elsewhere ((dlatcen > 15.).and.(dlatcen <= 30.))
   !             inox1 = 3      ! Subtropical 
   !          elsewhere ((dlatcen >= -40.).and.(dlatcen < -15.))
   !             inox1 = 3      ! Subtropical 
   !          elsewhere 
   !             inox1 = 4      ! Midlatitude 
   !          end where
   ! !     In-between months     
   !       case (4,5,10,11)
   !          where     ((abs(dlatcen) <= 15.).and.(cmi_flags == 2))
   !             inox1 = 1      ! Tropical continental
   !          elsewhere     ((abs(dlatcen) <= 15.).and.(cmi_flags /= 2))
   !             inox1 = 2      ! Tropical marine
   !          elsewhere ((dlatcen > 15.).and.(dlatcen <= 30.))
   !             inox1 = 3      ! Subtropical 
   !          elsewhere ((dlatcen >= -30.).and.(dlatcen < -15.))
   !             inox1 = 3      ! Subtropical 
   !          elsewhere 
   !             inox1 = 4      ! Midlatitude 
   !          end where
   ! !     Northern Hemisphere Summer      
   !       case (6,7,8,9)
   !          where     ((abs(dlatcen) <= 15.).and.(cmi_flags == 2))
   !             inox1 = 1      ! Tropical continental
   !          elsewhere     ((abs(dlatcen) <= 15.).and.(cmi_flags /= 2))
   !             inox1 = 2      ! Tropical marine
   !          elsewhere ((dlatcen > 15.).and.(dlatcen <= 40.))
   !             inox1 = 3      ! Subtropical 
   !          elsewhere ((dlatcen >= -30.).and.(dlatcen < -15.))
   !             inox1 = 3      ! Subtropical 
   !          elsewhere 
   !             inox1 = 4      ! Midlatitude 
   !          end where
   !       case default
   !          print*, "CASE default in partitionnox"
   !       end select 	 

   !     Define cloud top to be highest layer with dtrn > 0.
   !     Caution: Code assumes cloud does not extend above model layer 30.
   !     Warning.  This statement should be examined as vertical res is increased.
   iktop = min(40,LM)
   
   ALLOCATE (ntop(1:IM,1:JM))
   
   !     Determine model layer of cloud top. new way:
   do ij=1,JM
      do il=1,IM
         do ik=LM,2,-1
            if (dtrn(il,ij,ik).gt.0) then
               ntop(il,ij) = ik
               goto 213
            endif
            ntop(il,ij) = 2
         enddo
213 continue
         enddo
      enddo
      
      !old way - 
      !      do ij=ju1,j2
      !      do il=i1,i2
      !         do ik=iktop,2,-1
      !            if (dtrn(il,ij,ik).gt.0) then
      !               ntop(il,ij) = ik
      !               goto 213
      !            end if
      !            ntop(il,ij) = 2
      !         end do
      ! 213     continue
      !      end do
      !      end do

      !     note: ppmtoppv converts from ppm to ppv.  N(MW) = 14.
      pnox3d = 0.
      

      ! Megan Damon
      ! adding allocating of htedge2
      allocate(htedge2(1:LM))
      allocate(fd0(0:LM))
      allocate(yout(1:LM))
      
      ppmtoppv = 28.97 / 14.
      
      !       do ij=ju1,j2
      !       do il=i1,i2
      !          if (pnox2d(il,ij).gt.0.0) then
      !             iktop = ntop(il,ij)                 !Index of top cloud layer
      !             ii = inox1(il,ij)                   !"Model-cloud index"
      
      ! !           CTM "cloud" is expanded or contracted so that its depth = 17 km, which is the
      ! !           assumed depth of the "model-cloud".
      !             htedge2 =  ztop(il,ij,:)*real(pmax(ii))/ztop(il,ij,iktop)
      
      !             yout = 0.
      !             fd0 = 0.
      ! 	    do ik=1,pmax(ii) 
      !                do ikk=1,ikmm
      !                   fd0(ikk) = max(min(htedge2(ikk)-(ik-1),1.),0.)
!                end do
      !                do ikk=1,ikmm-1
      !                   yout(ikk) = yout(ikk)+r0(ik,ii)*(fd0(ikk)-fd0(ikk-1))
      !                end do
      !             enddo
      !             yout = yout * 1. / sum(yout)
      
      !             do ikk=1,iktop
      !                pnox3d(il,ij,ikk)=yout(ikk)*pnox2d(il,ij)
      !             end do
      
      !          end if
      !       end do
      !       end do
      
      
      !     Convert from g N s-1 to ppv s-1.
      pnox3d = 1.0d0 ! remove this once there are actual values from above loop
      pnox3d = 0.001 * pnox3d * ppmtoppv / mass

      deallocate(ztop)
      deallocate(htedge2)
      deallocate(fd0)
      deallocate(yout)
      deallocate(inox1)
      deallocate(ntop)
      
      return
    end subroutine partitionnox


!=============================================================================



!-------------------------------------------------------------------------
!BOP
!
! !INTERFACE:
    subroutine calculateProductionNox (midLatAdj, flashrate, desired_g_N_prod_rate, pnox)

! !USES:
      implicit none

! !INPUT PARAMETERS:
      real*8, intent(in) :: midLatAdj(:,:) 

      !     Total (CG+IC) flash rate (flashes per grid box per s)
      real*8, intent(in) :: flashrate(:,:)


      !     Desired mean lightning NO production rate (Tg yr-1) (specified in namelist)       	
      real*8, intent(in) :: desired_g_N_prod_rate

! !OUTPUT PARAMETERS:
        !     Lightning NO production rate (g N per grid box per s)
        real*8, intent(out) :: pnox(:,:)

! !DESCRIPTION
! Calculate N production rate (g s-1).
! For specified flashrate and PRODFAC, global production rate = 3.41 Tg N yr-1.
! Assume IC and CG flash rates produce the same amount of N / flash.
!      pnox = flashrate(:,:)*PRODFAC*CONVFAC
!EOP
!-------------------------------------------------------------------------

        pnox = midLatAdj(:,:)*flashrate(:,:)*PRODFAC*CONVFAC 
        
        !adjust for desired global Nitrogen production rate
        pnox = pnox * desired_g_N_prod_rate /3.41        
        
      end subroutine calculateProductionNox

      subroutine readLightRatioGlobalData &
           &  (light_ratioGlobal_infile_name, ratioGlobalLight, nym, pr_diag, rootProc)


        character (len=*) :: light_ratioGlobal_infile_name
        real*8, intent(out) :: ratioGlobalLight
        integer, intent(in) :: nym
        logical, intent(in) :: pr_diag, rootProc

        integer :: ii, mm, asc_lun, ierr
        real*8  :: readRatio


        character (len=256) :: err_msg

        integer :: num_reads = 12 * 37 ! roughly 1 entry for each month in satellite era

        asc_lun = 9

        OPEN(UNIT=asc_lun, FILE=TRIM(light_ratioGlobal_infile_name), STATUS='old', ACTION='read', &
             FORM='formatted', ACCESS='sequential', IOSTAT=ierr)

        if (ierr /= 0) then
           err_msg = 'Failed to OPEN '//TRIM(light_ratioGlobal_infile_name)
           stop
        end if

        ratioGlobalLight = 0.0

        do mm = 1, num_reads
           Read (asc_lun, *) ii, readRatio
           if (ii == nym) then
              ratioGlobalLight = readRatio
            exit
         end if
      end do

      Close (asc_lun)

      IF(rootProc) THEN
         WRITE(6,*) 'Ratio global for lightning read from', ratioGlobalLight
      END IF

      return

    end subroutine readLightRatioGlobalData

!========================================================================================



      
  end module Lightning_mod
