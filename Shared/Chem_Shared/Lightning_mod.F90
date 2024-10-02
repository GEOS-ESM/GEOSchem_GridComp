!  $Id: Lightning_mod.F90,v 1.1.4.1.2.1.2.1 2021/08/06 16:16:34 mmanyin Exp $

module Lightning_mod

#include "MAPL_Generic.h"

!----------------------------------------------------------------------
!BOP

! !MODULE: 

!    Lightning_mod -- Container for the GEOS lightning utilities.

! !USES:

use ESMF
use MAPL
use GEOS_UtilsMod
use Lightning_Toolbox_Mod, only : CalcFlashRate 
use m_set_eta, only: set_eta
use Chem_UtilMod, only: Chem_UtilResVal

use, intrinsic :: iso_fortran_env, only: REAL64



   implicit none
   private

  
! PUBLIC MEMBER FUNCTIONS:
!

   PUBLIC :: getLightning
   PUBLIC :: read_flash_source
   PUBLIC :: read_lightning_config
   PUBLIC :: update_lightning_ratio
   PUBLIC :: HEMCO_FlashRate       ! Had been called from MOIST to provide a secondary version of LFR

   ! Enumerated values, along with matching strings
   integer, parameter, public    :: FLASH_SOURCE_MOIST      = 1   ! values also serve as indices into flashSourceNames
   integer, parameter, public    :: FLASH_SOURCE_FIT        = 2
   integer, parameter, public    :: FLASH_SOURCE_HEMCO      = 3
   integer, parameter, public    :: FLASH_SOURCE_LOPEZ      = 4
   integer, parameter, public    :: FLASH_SOURCE_UNDEFINED  = 5

   integer, parameter, public    :: FLASH_SOURCE_count      = FLASH_SOURCE_UNDEFINED-1

   character (len=9),  public    :: flashSourceNames(5) = (/ 'MOIST    ', 'FIT      ', 'HEMCO    ', 'LOPEZ    ', 'UNDEFINED' /)


   PRIVATE :: MOIST_FlashRate
   PRIVATE ::   FIT_FlashRate ! fields passed are vertically flipped from GEOS orientation 
   PRIVATE :: LOPEZ_FlashRate
   PRIVATE :: partition
   PRIVATE :: readLightRatioGlobalData
   PRIVATE :: BUOYANCY
   PRIVATE :: identify_flash_source
   PRIVATE :: emiss_lightning


! !PARAMETERS:

! #include "gmi_phys_constants.h"    (Had been included in GmiEmissionLightning_mod)

   real,         parameter :: MWT_NO       = 30.0064                  ! molecular weight of NO
   real(REAL64), parameter :: CONVFAC      = 2.33e-23                 ! 14g of N per mole of NO/ 6.02e23 molecules of NO/mole of NO.
   real(REAL64), parameter :: PRODFAC      = 1.0e26                   ! (joules per flash) * molecules of NO per joule
   real(REAL64), parameter :: Pa2hPa       = 0.01
   real(REAL64), parameter :: REF_PRESS    = 412.5d0                  ! Dale Allen's new FIT routine mb/hPa
   real(REAL64), parameter :: SEC_PER_MIN  = 60.0
   real(REAL64), parameter :: MIN_PER_SEC  = (1.0/SEC_PER_MIN)
   real(REAL64), parameter :: SEC_PER_DAY  = 86400.0   ! Seconds per day
   real(REAL64), parameter :: DAY_PER_SEC  = (1.0/SEC_PER_DAY)

   integer,      parameter :: hp = KIND( REAL( 0.0, 8 ) ) ! HEMCO type, as in the Lightning Toolbox
   integer,      parameter :: r8 = 8
   integer,      parameter :: r4 = 4



! !DESCRIPTION:

!  This module contains parameterizations and utilities pertaining to the calculation of
!  flashrate (or strokerate) and NOx from lightning in the Earth's atmosphere.
!  The routines for calculation flashrate / stroke rate are:
! 
!      FIT_FlashRate - adapted from the offline GMI-CTM model (Dale Allen)
!    MOIST_FlashRate - adapted from MOIST / CTM cinderalla component
!    HEMCO_FlashRate - GEOS-Chem's flashrate calculation 
!    LOPEZ_FlashRate (WARNING: this routine is producing invalid results!)
!   
!  Please review the documentation below for each routine carefully. 
!
!
!
! !REVISION HISTORY:
!
! March 17 2018 - Megan Damon. Integrating code from the GMI-CTM.
! October 15 2018 - Megan Damon. Integrating code from GEOS-5 for MOIST_FlashRate
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
!  ORIGIN AND CONTACT
!
! !REVISION HISTORY:
! 7  Nov 2019 Damon      First crack
! 18 Dec 2019 Damon      Code review with Tom Clune
!    Sep 2020 Manyin     Refactored
!EOP
!-----------------------------------------------------------------------

subroutine getLightning (GC, ggState, CLOCK, &
     flash_source_enum, ratioGlobalLight, &
     DT, &
     LONSLOCAL, LATSLOCAL, & 
     CNV_PRCP, FRLAND, FROCEAN, LWI, PBLH, mcor, cellArea, midLatAdj, ratioLocal, TS, &
     CNV_MFC, CNV_QC, T, TH, &
     PFICU, PLE, Q, ZLE, &
     minDeepCloudTop, lightNOampFactor, numberNOperFlash, &
     MOIST_flashFactor, FIT_flashFactor, HEMCO_flashFactor, LOPEZ_flashFactor, &
     CNV_MFD, usePreconCape, CAPE_PRECON, INHB_PRECON, BYNCY_PRECON, &
     CAPE, BYNCY, flashRate, light_NO_prod, PHIS,  &
     RC) 

  type(ESMF_GridComp),           intent(inout)    :: GC     ! Gridded component

  type (MAPL_MetaComp), pointer, intent(in)       :: ggState
  type (ESMF_Clock),             intent(inout)    :: CLOCK 

  integer,                       intent(in)       :: flash_source_enum  ! indicating MOIST, FIT, HEMCO, or LOPEZ
  real,                          intent(in)       :: ratioGlobalLight   ! only needed for FIT
  real,                          intent(in)       :: MOIST_flashFactor
  real,                          intent(in)       :: FIT_flashFactor
  real,                          intent(in)       :: HEMCO_flashFactor
  real,                          intent(in)       :: LOPEZ_flashFactor

  real(r8) :: DT ! model heart beat (seconds)      

  real, intent(in) :: LONSLOCAL(:,:), LATSLOCAL(:,:) ! degrees of lat/lon coordinates 
                                                     ! MRD: these are required for HEMCO flash rate
                                                     ! tropics, midlats differentiated  
                                                     ! Eurasian and America flashrates also differentiated
   

  real, dimension(:,:),     intent(in)  :: CNV_PRCP ! convective precip  (kg m-2 s-1)
  real, dimension(:,:),     intent(in)  :: FRLAND   ! land fraction (none)
  real, dimension(:,:),     intent(in)  :: FROCEAN  ! land fraction (none)
  real, dimension(:,:),     intent(in)  :: LWI      ! land water ice flag  (0=water; 1=land; 2=ice)
  real, dimension(:,:),     intent(in)  :: PBLH     ! planetary boundary layer height (m)

  ! Only associated when (flash_source_enum == FLASH_SOURCE_FIT)
  real, pointer, dimension(:,:),     intent(in)  :: mcor       ! cellArea for GMI (m2)
  real, pointer, dimension(:,:),     intent(in)  :: midLatAdj  ! input field for FIT (1)
  real, pointer, dimension(:,:),     intent(in)  :: ratioLocal ! input field for FIT (1) 

  real, dimension(:,:),     intent(in)  :: cellArea   ! cellArea for this model (m2)
  real, dimension(:,:),     intent(in)  :: TS         ! surface temperature (K)
  real, dimension(:,:),     intent(in)  :: PHIS       ! geopotential height at the surface

  real, dimension(:,:,:),   intent(in)  :: CNV_QC   ! grid mean convective condensate (?)

  ! edge vars  (vertical indices: 1 to LM+1)
  real, dimension(:,:,:),   intent(in)  :: CNV_MFC  ! cumulative mass flux (kg m-2 s-1)   [top-down]
  real, dimension(:,:,:),   intent(in)  :: PLE      ! edge pressures (Pa)
  real, dimension(:,:,:),   intent(in)  :: ZLE      ! geopotential height (m)             [top-down]

  real, dimension(:,:,:),   intent(in)  :: PFICU    ! flux of ice convective precip (kg m-2 s-1)
  real, dimension(:,:,:),   intent(in)  :: T        ! air temperature (K)                 [top-down]
  real, dimension(:,:,:),   intent(in)  :: TH       ! potential temperature (K)           [top-down]
  real, dimension(:,:,:),   intent(in)  :: Q        ! specific humidity (kg kg-1)         [top-down]

  real,                     intent(in)  :: minDeepCloudTop    ! Minimum cloud top [km] for selecting deep convection profiles
  real,                     intent(in)  :: lightNOampFactor   ! > 0, for targeting the observed nitrogen production rate [3.41 Tg yr^{-1}]
  real,                     intent(in)  :: numberNOperFlash   ! NO molecules generated by each flash
  real, dimension(:,:,:),   intent(in)  :: CNV_MFD            ! detraining_mass_flux  (kg m-2 s-1)     [top-down]

  logical,                  intent(in)  :: usePreconCape   ! Whether to use CAPE, INHB and BYNCY from MOIST
                                                           ! (only affects LOPEZ and MOIST)
                                                           ! LOPEZ only uses CAPE, and should benefit from PRECON
                                                           ! MOIST prefers CAPE computed using the old MERRA2 approach
                                                           !       which needs BYNCY, so may benefit from PRECON
  real, dimension(:,:),     intent(in)  ::  CAPE_PRECON
  real, dimension(:,:),     intent(in)  ::  INHB_PRECON
  real, dimension(:,:,:),   intent(in)  :: BYNCY_PRECON

  real, dimension(:,:),     intent(out)  :: CAPE           ! can remove these eventually
  real, dimension(:,:,:),   intent(out)  :: BYNCY          ! can remove these eventually               [top-down]
  real, dimension(:,:),     intent(out)  :: flashRate      ! (flashes/km2/s)

  ! May be NULL:
  real*4, pointer, dimension(:,:,:),   intent(inout)  :: light_NO_prod  ! Production of NO from lightning (m-3 s-1) [top-down]

  integer, optional, intent(out) :: RC

! type (ESMF_Grid)  :: esmfGrid

  ! Lopez calculation; TODO these variables can be contained within a contains block
  real, allocatable :: PLO(:,:,:)     !  pressure at middle of gridbox [hPa]
  real, allocatable :: PKE(:,:,:)
  real, allocatable :: PK(:,:,:)      
  real, allocatable :: QSS(:,:,:)     !  saturation specific humidity
  real, allocatable :: DQS(:,:,:)     !  derivative of satuation specific humidity wrt temperature
  real, allocatable :: DZ(:,:,:)
  real, allocatable :: DM(:,:,:)
  real, allocatable :: THV(:,:,:)
  real, allocatable :: ZLO(:,:,:)
  real, allocatable :: ZLE2(:,:,:)
  real, allocatable :: CAPE_MERRA2(:,:)
  real, allocatable :: INHB(:,:)

  integer, allocatable :: LWI_INT(:,:)
  real*8, allocatable :: LFR_R8(:,:)
  real*4, allocatable :: gridBoxThickness_bottomup_R4(:,:,:)
  real*4, allocatable :: DTRN_bottomup_R4(:,:,:)


  real, allocatable :: cldmas0(:,:) ! Dale Allen's Fit routine

  logical           :: need_cape
      
  real,       pointer,      dimension(:,:)   ::  CNV_TOPP => null()

  real(REAL64),  allocatable   :: flashRateDale  (:,:)
  real,          allocatable   :: flashRateLopez (:,:)

  real(REAL64),  allocatable   :: AK(:)
  real(REAL64),  allocatable   :: BK(:)
  real(REAL64),  allocatable   :: PREF(:)

  logical, save :: first = .TRUE.
  integer       :: IM, JM, LM, K0, KM, nc
  integer       :: i, j, L, ls
  integer       :: modelLevel ! only used for FIT
  real(REAL64)  :: ptop, pint

  integer :: STATUS
  integer :: shape3d (3)

  integer :: mon, year, nym
  type(ESMF_Time) :: curr_date             ! temporary date used in logic
  character(len=80) :: cYear
  character(len=80) :: cMonth
  character(len=80) :: cNym
  character(len=3)  :: binName

  character(len=*), parameter :: Iam = "getLightning" ! can be moduleName // getLightning TODO

  if (PRESENT(RC)) RC = ESMF_SUCCESS

  shape3d = shape(TH)

  IM = shape3d(1)
  JM = shape3d(2)
  LM = shape3d(3)

  K0 = 1           ! = LBOUND of deferred shape EDGE array
  KM = LM + 1      ! = UBOUND of deferred shape EDGE array

  nc = IM * JM

  if (first) then
    if( MAPL_AM_I_ROOT() ) then
      print*, "I am getLightning, flash source: ", TRIM(flashSourceNames(flash_source_enum)) ! TODO every process is problematic (writeParallel); maybe with underscore
      print*, TRIM(flashSourceNames(flash_source_enum)), " flash rate calculation"
    endif
    first = .FALSE.
  endif
     
  need_cape = ((flash_source_enum == FLASH_SOURCE_LOPEZ) .OR.  &
               (flash_source_enum == FLASH_SOURCE_HEMCO) .OR.  &
               (flash_source_enum == FLASH_SOURCE_MOIST))

  need_cape = .TRUE.

  IF ( need_cape ) THEN

     ALLOCATE(        INHB(IM, JM),         __STAT__ )

     ! on gridbox centers:
     ALLOCATE(         PLO(IM, JM, 1:LM),   __STAT__ )
     ALLOCATE(          PK(IM, JM, 1:LM),   __STAT__ )
     ALLOCATE(         DQS(IM, JM, 1:LM),   __STAT__ )
     ALLOCATE(         QSS(IM, JM, 1:LM),   __STAT__ )
     ALLOCATE(          DZ(IM, JM, 1:LM),   __STAT__ )
     ALLOCATE(         THV(IM, JM, 1:LM),   __STAT__ )
     ALLOCATE(         ZLO(IM, JM, 1:LM),   __STAT__ )

     ! on gridbox edges:
     ALLOCATE(         PKE(IM, JM, 0:LM),   __STAT__ )
     ALLOCATE(        ZLE2(IM, JM, 0:LM),   __STAT__ )


     PLO = 0.5*(PLE(:,:,K0:KM-1) +  PLE(:,:,K0+1:KM  ) )*0.01
     PKE(:,:,0:LM) = (PLE(:,:,K0:KM)*(1.0/MAPL_P00))**(MAPL_RGAS/MAPL_CP)

     PK       = (PLO*(100.0/MAPL_P00))**(MAPL_RGAS/MAPL_CP)

     DQS = GEOS_DQSAT (TH*PK, PLO, qsat=QSS)

     THV(:,:,:) = TH(:,:,:) * (1.+MAPL_VIREPS*Q(:,:,:))

     ZLE2(:,:,LM) = 0.               !   ORIGINAL FORMULATION - does Saulo need this?
     ZLE2(:,:,LM) = PHIS/MAPL_GRAV   ! Better match for archived ZLE
     do L=LM,1,-1
         ZLO(:,:,L  ) = ZLE2(:,:,L) + (MAPL_CP/MAPL_GRAV)*( PKE(:,:,L)-PK (:,:,L  ) ) * THV(:,:,L)
         DZ (:,:,L  ) =               (MAPL_CP/MAPL_GRAV)*( PKE(:,:,L)-PKE(:,:,L-1) ) * THV(:,:,L)
        ZLE2(:,:,L-1) = ZLE2(:,:,L) + DZ(:,:,L)
!       ZLE2(:,:,L-1) = ZLE2(:,:,L) + (MAPL_CP/MAPL_GRAV)*( PKE(:,:,L)-PKE(:,:,L-1) ) * THV(:,:,L)
     end do

! TODO - decide if we can remove ZLE or ZLE2 since they are nearly duplicates
!        although the indices are mismatched (ZLE starts at 1, ZLE2 starts at 0)
!  IF(MAPL_AM_I_ROOT()) THEN
!    DO i=K0,KM
!      PRINT *,'ZLE, ZLE2:', i, ZLE(1,1,i), ZLE2(1,1,i-1), ZLE(1,1,i)-ZLE2(1,1,i-1)
!    ENDDO
!  ENDIF

!      IF(MAPL_AM_I_ROOT()) THEN
!        PRINT*,' ZLE vert: ', LBOUND(ZLE, 3), UBOUND(ZLE, 3)
!        PRINT*,'ZLE2 vert: ', LBOUND(ZLE2,3), UBOUND(ZLE2,3)
!        PRINT*,' PLE vert: ', LBOUND(PLE, 3), UBOUND(PLE, 3)
!   
!   !  ZLE vert:            1          73
!   ! ZLE2 vert:            0          72
!   !  PLE vert:            1          73
!   
!      ENDIF

     ! ALL SET to call BUOYANCY (T, Q, QSS, DQS, DZ, ZLO, BYNCY, CAPE, INHB)

     ! DEALLOCATE vars at the end

  END IF   ! need_cape



  !-----------------------------------------------------------------------
  ! Lopez (2016?) flash rate from ECMWF
  !-----------------------------------------------------------------------         
  if (flash_source_enum == FLASH_SOURCE_LOPEZ) then

     ! callLopezCalcuations() put above end subroutine getLightning

     ! MEM- Recommend PRECON, but leaving the other option available
     if (usePreconCape) then
       CAPE  =  CAPE_PRECON
       INHB  =  INHB_PRECON
       BYNCY = BYNCY_PRECON
     else
       call BUOYANCY (T, Q, QSS, DQS, DZ, ZLO, BYNCY, CAPE, INHB)
     end if

     ALLOCATE( flashRateLopez(IM, JM),   STAT=STATUS); VERIFY_(STATUS)
     flashRateLopez = real(0)
!    flashRateLopez(:,:) = 0.01

      call LOPEZ_FlashRate(IM, JM, LM, FRLAND, PBLH, CAPE, ZLE2, PFICU, &
           CNV_QC, T, PLO, LOPEZ_flashFactor, flashRateLopez)

!print*, "min/max LFR LOPEZ: ", minval(flashRateLopez), maxval(flashRateLopez)
           
     ! flashrates come back as flashes/km2/day
     ! LFR should be in flashes/km2/s
     flashRate(:,:) = flashRateLopez(:,:) * DAY_PER_SEC 

! print*,'Lopez - minmax Flash Rate: ', MINVAL(flashrate), MAXVAL(flashrate)
! print*,'Lopez - minmax BUOY: ', MINVAL(BYNCY), MAXVAL(BYNCY)
! print*,'Lopez - minmax CAPE: ', MINVAL(CAPE), MAXVAL(CAPE)
! print*,'Lopez - max ZLE, ZLE2: ', MAXVAL(ZLE), MAXVAL(ZLE2)

     DEALLOCATE(flashRateLopez, __STAT__)


  !-----------------------------------------------------------------------
  ! HEMCO flash rate
  !-----------------------------------------------------------------------
  else if (flash_source_enum == FLASH_SOURCE_HEMCO) then

!    The HEMCO flash routine may use BYNCY in the future
!    if (usePreconCape) then
!      CAPE  =  CAPE_PRECON
!      INHB  =  INHB_PRECON
!      BYNCY = BYNCY_PRECON
!    else
!      call BUOYANCY (T, Q, QSS, DQS, DZ, ZLO, BYNCY, CAPE, INHB)
!    end if


     call HEMCO_FlashRate (cellArea, LWI, LONSLOCAL, LATSLOCAL, T, PLE, &
          ZLE, CNV_MFC, HEMCO_flashFactor, flashRate, __RC__ )

!print*, "min/max LFR HEMCO: ", minval(flashRate), maxval(flashRate)
     ! HEMCO flashrates returned as: flashes/km2/s


  !-----------------------------------------------------------------------
  ! MOIST-derived flashrate calculation (Dale Allen, old)
  !-----------------------------------------------------------------------
  else if (flash_source_enum == FLASH_SOURCE_MOIST) then

     ALLOCATE(          DM(IM, JM, 1:LM),   __STAT__ )
     ALLOCATE( CAPE_MERRA2(IM, JM),         __STAT__ )
     ALLOCATE(    CNV_TOPP(IM, JM),         __STAT__ )

     DM(:,:,1:LM) = ( PLE(:,:,K0+1:KM)-PLE(:,:,K0:KM-1) ) * (1./MAPL_GRAV)  ! DELP / g

     ! (Don't use CAPE_PRECON since the units are different from CAPE_MERRA2)
     ! The MOIST approach wants CAPE computed in the old MERRA2 way
     ! which uses BYNCY, so usePreconCape really means usePreconBuoyancy in this case
     ! Probably best to use PRECON
     if (usePreconCape) then

       CAPE  =  CAPE_PRECON
       INHB  =  INHB_PRECON
       BYNCY = BYNCY_PRECON

       CAPE_MERRA2(:,:) = 0.

       do L=1,LM-1
          where(BYNCY(:,:,L)>0.)
             CAPE_MERRA2 = CAPE_MERRA2 + BYNCY(:,:,L)*DM(:,:,L)
          end where
       end do

       where(CAPE_MERRA2 <= 0.0)
          CAPE_MERRA2=MAPL_UNDEF
       end where

     else
       call BUOYANCY (T, Q, QSS, DQS, DZ, ZLO, BYNCY, CAPE, INHB, &
                      DM=DM, CAPE_MERRA2=CAPE_MERRA2)
     end if




     ! Determine the pressure at convective cloud top

     ! WHERE CNV_MFC != 0 ... TODO
     CNV_TOPP(:,:) = MAPL_UNDEF
     do j=1, JM
        do i=1, IM
           do l=1, LM
              if (CNV_MFC(i,j,l)/=0.0) then
                  CNV_TOPP(i,j) = PLE(i,j,l)
                  EXIT    ! exit the top-down loop
              endif
           enddo
        enddo
     enddo
!print*, 'count of TOPP >= 50000:', COUNT(CNV_TOPP >= 50000.)

     call MOIST_FlashRate(ggState, nc, LM, TS, CNV_TOPP, FROCEAN, &
          CNV_PRCP, CAPE_MERRA2, CNV_MFC, TH, PLE, ZLE, MOIST_flashFactor, &
          flashRate, __RC__ ) !TODO signify Dale Allen Old

!print*, "min/max LFR MOIST: ", minval(flashRate), maxval(flashRate)

     DEALLOCATE( DM, CAPE_MERRA2, CNV_TOPP, __STAT__ )


  !-----------------------------------------------------------------------
  ! Dale Allen flashrate calculation adapted from the GMI-CTM
  !-----------------------------------------------------------------------
  else if (flash_source_enum == FLASH_SOURCE_FIT) then

     ALLOCATE ( AK(1:LM+1), __STAT__ )
     ALLOCATE ( BK(1:LM+1), __STAT__ )
     ALLOCATE ( PREF(0:LM), __STAT__ )

!    call ESMF_GridCompGet ( GC, Grid=esmfGrid, __RC__ )
!    call ESMF_AttributeGet(esmfGrid,name="GridAK",valuelist=AK, __RC__ )
!    call ESMF_AttributeGet(esmfGrid,name="GridBK",valuelist=BK, __RC__ )

     call set_eta(LM,ls,ptop,pint,AK,BK)

     PREF(0:LM) = AK(1:LM+1) + BK(1:LM+1) * MAPL_P00


     ! PREF is in Pa, while pressure is in mb
     modelLevel  = max(1,count(PREF <= REF_PRESS * 100.0 )) ! TODO no magic number 
     
     ALLOCATE ( cldmas0(IM,JM), stat=STATUS)

     ! PREF is on 0-LM; CNV_MFC is on 1-LM
     cldmas0 (:,:) = CNV_MFC(:,:,modelLevel+1)


     ALLOCATE (flashRateDale (IM,JM), stat=STATUS)


     call FIT_FlashRate (cldmas0, 0.0d0, ratioLocal, ratioGlobalLight, midLatAdj, &
          5.0d0, FIT_flashFactor, flashRateDale)

     ! flashRateDale has units  [flashes / gridbox / sec]
     ! flashrates from Dale's routine are assumed on a particular grid, that is
     ! non-CS-based. This adjustment scales the flashrates
     ! from the expected (mcor-GMI) to the actual grid (CS-assumed)

!    flashRate (:,:) = flashRateDale(:,:) / ( mcor(:,:)/cellArea(:,:) )
!    flashRate (:,:) = flashRate(:,:) / (cellArea(:,:) / 1e6)

     ! MEM - Probably should use callArea instead of mcor here:
     flashRate (:,:) = 1e6 * flashRateDale(:,:) /  mcor(:,:)

!print*, "min/max LFR FIT: ", minval(flashRateDale), maxval(flashRateDale)
     BYNCY (:,:,:) = real(0) ! MAPLE_UNDEF
     CAPE  (:,:)   = real(0) ! MAPLE_UNDEF

     DEALLOCATE (flashRateDale, AK, BK, PREF, cldmas0)
 
       
  else 
     print*, ""
     print*, "Flashrate source not supported!"
     print*, ""
     STATUS=88
     VERIFY_(STATUS)
  endif

  IF ( need_cape ) THEN
     DEALLOCATE( INHB, PLO, PK, DQS, QSS, DZ, THV, ZLO, PKE, ZLE2,  __STAT__ )
  END IF

!print*, "min/max Final LFR: ", minval(flashRate), maxval(flashRate)


! CALL EMISS

  if ( associated(light_NO_prod) ) then

    ALLOCATE (                     LWI_INT(IM,JM   ), __STAT__)
    ALLOCATE (                      LFR_R8(IM,JM   ), __STAT__)
    ALLOCATE (gridBoxThickness_bottomup_R4(IM,JM,LM), __STAT__)
    ALLOCATE (            DTRN_bottomup_R4(IM,JM,LM), __STAT__)

    LFR_R8 = DBLE(flashRate)
!   LWI_INT(:,:) = -99
!   WHERE(LWI > -0.5 .AND. LWI < 0.5 ) LWI_INT = 0  + 1
!   WHERE(LWI >  0.5 .AND. LWI < 1.5 ) LWI_INT = 1  + 1
!   WHERE(LWI >  1.5 .AND. LWI < 2.5 ) LWI_INT = 2  + 1

    LWI_INT = FLOOR( LWI+0.1 )

    DO L=1,LM
      gridBoxThickness_bottomup_R4(:,:,(LM+1-L)) =     zle(:,:,L)-zle(:,:,L+1)
                  DTRN_bottomup_R4(:,:,(LM+1-L)) = CNV_MFD(:,:,L)
    END DO

    CALL emiss_lightning (1, IM, 1, JM, 1, LM, &
                          minDeepCloudTop, lightNOampFactor, numberNOperFlash, &
                          LWI_INT, LFR_R8, gridBoxThickness_bottomup_R4, DTRN_bottomup_R4, &
                          light_NO_prod, __RC__)

    ! reverse the vertical, making light_NO_prod  top-down:
    ! (use DTRN as a temp field)
    DTRN_bottomup_R4(:,:,1:LM) = light_NO_prod(:,:,LM:1:-1)
    light_NO_prod = DTRN_bottomup_R4

    DEALLOCATE (LWI_INT, LFR_R8, gridBoxThickness_bottomup_R4, DTRN_bottomup_R4, __STAT__)

  end if

end subroutine getLightning

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: read_flash_source
!
! !DESCRIPTION: Read lightning flash source from the resource file
!  
! !ORIGIN AND CONTACT:  Michael Manyin
!
! !REVISION HISTORY:
! 16 Apr 2021 Manyin     First crack
!EOP
!-----------------------------------------------------------------------

subroutine read_flash_source ( rcfilen, flash_source_enum, RC )

  ! IN:
  character(len=ESMF_MAXSTR),  intent(in)   :: rcfilen

  ! OUT:
  integer,                     intent(out)  :: flash_source_enum
  integer, optional,           intent(out)  :: RC

  ! Local vars:
  character(len=ESMF_MAXSTR)   :: IAm
  integer                      :: STATUS
  type (ESMF_Config)           :: esmfConfig
  character(len=ESMF_MAXSTR)   :: flashSource

  IAm = "read_flash_source"

  esmfConfig = ESMF_ConfigCreate(__RC__)

  call ESMF_ConfigLoadFile(esmfConfig, TRIM(rcfilen), __RC__ )

  ! How to calculate flashrate
  call ESMF_ConfigGetAttribute( esmfConfig, flashSource,    &
                                Label    = "flashSource:",  &
                                Default  = 'MOIST',  __RC__ )

  call identify_flash_source( flashSource, flash_source_enum )
  IF ( flash_source_enum == FLASH_SOURCE_UNDEFINED ) THEN
    print*,'Invalid Flash Source: '//TRIM(flashSource)//' from '//TRIM(rcfilen)
    STATUS = 101
    VERIFY_(STATUS)
  END IF

  call ESMF_ConfigDestroy(esmfConfig, __RC__)

end subroutine read_flash_source

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: read_lightning_config
!
! !DESCRIPTION: Retrieve values from the resource file
!  
! !ORIGIN AND CONTACT:  Michael Manyin
!
! !REVISION HISTORY:
! 11 Sep 2020 Manyin     First crack
! 16 Apr 2021 Manyin     Now use resolution vectors
!EOP
!-----------------------------------------------------------------------

subroutine read_lightning_config ( im, jm, rcfilen, flash_source_enum,   &
                                   SimType,                              &
                                   ratioGlobalFile, minDeepCloudTop,     &
                                   lightNOampFactor, numberNOperFlash,   &
                                   MOIST_flashFactor, FIT_flashFactor,   &
                                   HEMCO_flashFactor, LOPEZ_flashFactor, &
                                   usePreconCape,                        &
                                   RC )

  ! IN:
  integer,                     intent(in)   :: im,jm
  character(len=ESMF_MAXSTR),  intent(in)   :: rcfilen
  integer,                     intent(in)   :: flash_source_enum
  character(len=*),            intent(in)   :: SimType   ! CTM, FREE or REPLAY

  ! OUT:
  character(len=ESMF_MAXSTR),  intent(out)  :: ratioGlobalFile
  real,                        intent(out)  :: minDeepCloudTop
  real,                        intent(out)  :: lightNOampFactor
  real,                        intent(out)  :: numberNOperFlash
  real,                        intent(out)  :: MOIST_flashFactor
  real,                        intent(out)  :: FIT_flashFactor
  real,                        intent(out)  :: HEMCO_flashFactor
  real,                        intent(out)  :: LOPEZ_flashFactor
  logical,                     intent(out)  :: usePreconCape
  integer, optional,           intent(out)  :: RC

  ! Local vars:
  integer, parameter           :: NHRES = 6
  character(len=ESMF_MAXSTR)   :: IAm
  integer                      :: STATUS
  type (ESMF_Config)           :: esmfConfig
  real                         :: resolution_vector(NHRES)

  IAm = "read_lightning_config"

  _ASSERT(  TRIM(SimType) == 'CTM' .OR. TRIM(SimType) == 'FREE' .OR. TRIM(SimType) == 'REPLAY', 'read_lightning_config: SimType must be CTM, FREE or REPLAY' )

  ! default values that should never be used:
  ratioGlobalFile   = "default_filename"
  MOIST_flashFactor = MAPL_UNDEF
    FIT_flashFactor = MAPL_UNDEF
  HEMCO_flashFactor = MAPL_UNDEF
  LOPEZ_flashFactor = MAPL_UNDEF

  esmfConfig = ESMF_ConfigCreate(__RC__)

  call ESMF_ConfigLoadFile(esmfConfig, TRIM(rcfilen), __RC__ )

  IF ( flash_source_enum == FLASH_SOURCE_MOIST ) THEN

     call ESMF_ConfigGetAttribute( esmfConfig, resolution_vector,                  &
                                   Label    = "MOIST_flashFactor_resvec_"//TRIM(SimType)//":",  __RC__ )

     MOIST_flashFactor = Chem_UtilResVal(im, jm, resolution_vector, __RC__)

  END IF

  IF ( flash_source_enum == FLASH_SOURCE_FIT ) THEN

     call ESMF_ConfigGetAttribute( esmfConfig, resolution_vector,                &
                                   Label    = "FIT_flashFactor_resvec_"//TRIM(SimType)//":",  __RC__ )

     FIT_flashFactor = Chem_UtilResVal(im, jm, resolution_vector, __RC__)

     call ESMF_ConfigGetAttribute( esmfConfig, ratioGlobalFile,                             &
                                   Label    = "ratioGlobalFile:",                           &
                                   Default  = '/home/mrdamon/Files/RatioGlobal.asc', __RC__ )

  END IF

  IF ( flash_source_enum == FLASH_SOURCE_HEMCO ) THEN

     ! NOTE: HEMCO_flashFactor is an alias for otdLisScale
     call ESMF_ConfigGetAttribute( esmfConfig, resolution_vector,                  &
                                   Label    = "HEMCO_flashFactor_resvec_"//TRIM(SimType)//":",  __RC__ )

     HEMCO_flashFactor = Chem_UtilResVal(im, jm, resolution_vector, __RC__)

  END IF

  IF ( flash_source_enum == FLASH_SOURCE_LOPEZ ) THEN

     ! NOTE: LOPEZ_flashFactor was originally called 'alpha'
     call ESMF_ConfigGetAttribute( esmfConfig, resolution_vector,                  &
                                   Label    = "LOPEZ_flashFactor_resvec_"//TRIM(SimType)//":",  __RC__ )

     LOPEZ_flashFactor = Chem_UtilResVal(im, jm, resolution_vector, __RC__)

  END IF


  ! minDeepCloudTop  = 7.0 --> Minimum cloud top [km] for selecting deep convection profiles
  ! lightNOampFactor = 1.0 --> NO production amplification/suppression factor, > 0
  ! numberNOperFlash = 1.50E+26 --> NO molecules generated by each flash


  call ESMF_ConfigGetAttribute(esmfConfig, minDeepCloudTop,    &
                               Label    = "minDeepCloudTop:",  &
                               Default  =  7.0,         __RC__ )

  call ESMF_ConfigGetAttribute(esmfConfig, lightNOampFactor,   &
                               Label    = "lightNOampFactor:", &
                               Default  =  1.0,         __RC__ )

  call ESMF_ConfigGetAttribute(esmfConfig, numberNOperFlash,   &
                               Label    = "numberNOperFlash:", &
                               Default  =  1.50E+26,    __RC__ )

  call ESMF_ConfigGetAttribute(esmfConfig, usePreconCape,      &
                               Label    = "usePreconCape:",    &
                               Default  =  .TRUE.,      __RC__ )

  call ESMF_ConfigDestroy(esmfConfig, __RC__)

end subroutine read_lightning_config

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
!  Manyin - For now we use version 1.137 to match Icarus-3_2_p9  CCM
!         - Added DM and CAPE_MERRA2 to return CAPE suitable for MOIST flash calc
!
! !REVISION HISTORY:
! 1 Nov 2019 Damon       Bringing routine to GEOS_Shared
!EOP
!-----------------------------------------------------------------------
  subroutine BUOYANCY( T, Q, QS, DQS, DZ, ZLO, BUOY, CAPE, INHB, DM, CAPE_MERRA2)


    ! !DESCRIPTION: Computes the buoyancy $ g \frac{T_c-T_e}{T_e} $ at each level
    !  for a parcel raised from the surface. $T_c$ is the virtual temperature of
    !  the parcel and $T_e$ is the virtual temperature of the environment.

    real, dimension(:,:,:),   intent(in)            :: T        ! air temperature (K)                 [top-down]
    real, dimension(:,:,:),   intent(in)            :: Q        ! specific humidity (kg kg-1)         [top-down]
    real, dimension(:,:,:),   intent(in)            :: QS       !  saturation specific humidity
    real, dimension(:,:,:),   intent(in)            :: DQS      !  derivative of satuation specific humidity wrt temperature
    real, dimension(:,:,:),   intent(in)            :: DZ       !  gridbox height
    real, dimension(:,:,:),   intent(in)            :: ZLO      !  geopotential height at center of gridbox

    real, dimension(:,:,:),   intent(out)           :: BUOY         !  [m/s/s]
    real, dimension(:,:),     intent(out)           :: CAPE         !  [J/kg]
    real, dimension(:,:),     intent(out)           :: INHB         !  [J/kg]

    real, dimension(:,:,:),   intent(in),  optional :: DM           !  [kg/m^2]  DELP/g
    real, dimension(:,:),     intent(out), optional :: CAPE_MERRA2  !  [kg/(m*s^2)]  Only use w/ Dale's orig lightning alg

    integer :: L, LM

    LM = size(T,3)

    BUOY(:,:,LM) =  T(:,:,LM) + (MAPL_GRAV/MAPL_CP)*ZLO(:,:,LM) + (MAPL_ALHL/MAPL_CP)*Q(:,:,LM)

    do L=LM-1,1,-1
       BUOY(:,:,L) = BUOY(:,:,LM) - (T(:,:,L) + (MAPL_GRAV/MAPL_CP)*ZLO(:,:,L) + (MAPL_ALHL/MAPL_CP)*QS(:,:,L))
       BUOY(:,:,L) = MAPL_GRAV*BUOY(:,:,L) / ( (1.+ (MAPL_ALHL/MAPL_CP)*DQS(:,:,L))*T(:,:,L) )
    enddo

    BUOY(:,:,LM) = 0.0

    CAPE(:,:) = 0.
    INHB(:,:) = 0.

    do L=1,LM-1
       where(BUOY(:,:,L)>0.)
          CAPE = CAPE + BUOY(:,:,L)*DZ(:,:,L)
       end where
       where(BUOY(:,:,L)<0.)
          INHB = INHB - BUOY(:,:,L)*DZ(:,:,L)
       end where
    end do

    where(CAPE <= 0.0)
       CAPE=MAPL_UNDEF
       INHB=MAPL_UNDEF
    end where

    if ( present(CAPE_MERRA2) ) then

      if ( present(DM) ) then

        CAPE_MERRA2(:,:) = 0.

        do L=1,LM-1
           where(BUOY(:,:,L)>0.)
              CAPE_MERRA2 = CAPE_MERRA2 + BUOY(:,:,L)*DM(:,:,L)
           end where
        end do

        where(CAPE_MERRA2 <= 0.0)
           CAPE_MERRA2=MAPL_UNDEF
        end where

      else
        print*,'BUOYANCY: Cannot compute old CAPE without DM'
      endif

    endif

!   print*, "min/max BUOY: ", minval(BUOY), maxval(BUOY)
!   print*, "min/max CAPE: ", minval(CAPE), maxval(CAPE)
!   print*, "min/max INHB: ", minval(INHB), maxval(INHB)

  end subroutine BUOYANCY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!-----------------------------------------------------------------------
!BOP
! !IROUTINE: MOIST_FlashRate
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
! 09 Sep 2020 Manyin      renamed computeFlashRate to be MOIST_FlashRate
!EOP
!-----------------------------------------------------------------------

  subroutine MOIST_FlashRate (STATE, nc, lm, TS, CCTP, FROCEAN, CN_PRCP, CAPE, &
       CNV_MFC, TH, PLE, ZLE, MOIST_flashFactor, strokeRate, RC)

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
    REAL, INTENT(IN), DIMENSION(nc) :: CAPE        ! Convective available potential energy [J m^{-2}]

    real, intent(in), dimension(nc,lm) :: TH        ! Potential temperature [K]
    real, intent(in), dimension(nc,0:lm) :: CNV_MFC ! Convective mass flux [kg m^{-2} s^{-1}]
    real, intent(in), dimension(nc,0:lm) :: PLE     ! Layer interface pressures  [Pa]
    real, intent(in), dimension(nc,0:lm) :: ZLE     ! Layer depths [m]

    real, intent(in)           :: MOIST_flashFactor ! Global scaling term
!
! !OUTPUT PARAMETERS:
    real, intent(out), dimension(nc) :: strokeRate ! Flashes per km^2 per second
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
    CHARACTER(LEN=*), PARAMETER :: Iam = "MOIST_FlashRate"





    ! Preliminaries
    ! -------------
    RC = ESMF_SUCCESS
    strokeRate(:) = 0.0
    

    ! Coefficients of the predictors, marine locations
    ! ------------------------------------------------
    CALL MAPL_GetResource(STATE,a0m,'MARINE_A0:',DEFAULT= 0.0139868, __RC__)
    CALL MAPL_GetResource(STATE,a1m,'MARINE_A1:',DEFAULT= 0.0358764, __RC__)
    CALL MAPL_GetResource(STATE,a2m,'MARINE_A2:',DEFAULT=-0.0610214, __RC__)
    CALL MAPL_GetResource(STATE,a3m,'MARINE_A3:',DEFAULT=-0.0102320, __RC__)
    CALL MAPL_GetResource(STATE,a4m,'MARINE_A4:',DEFAULT= 0.0031352, __RC__)
    CALL MAPL_GetResource(STATE,a5m,'MARINE_A5:',DEFAULT= 0.0346241, __RC__)


    ! Coefficients of the predictors, continental locations
    ! -----------------------------------------------------
    CALL MAPL_GetResource(STATE,a0c,'CONTINENT_A0:',DEFAULT=-0.0183172, __RC__)
    CALL MAPL_GetResource(STATE,a1c,'CONTINENT_A1:',DEFAULT=-0.0562338, __RC__)
    CALL MAPL_GetResource(STATE,a2c,'CONTINENT_A2:',DEFAULT= 0.1862740, __RC__)
    CALL MAPL_GetResource(STATE,a3c,'CONTINENT_A3:',DEFAULT=-0.0023363, __RC__)
    CALL MAPL_GetResource(STATE,a4c,'CONTINENT_A4:',DEFAULT=-0.0013838, __RC__)
    CALL MAPL_GetResource(STATE,a5c,'CONTINENT_A5:',DEFAULT= 0.0114759, __RC__)

    ! Divisors for nondimensionalization of the predictors
    ! ----------------------------------------------------
    CALL MAPL_GetResource(STATE,x1Divisor,'X1_DIVISOR:',DEFAULT=4.36, __RC__)
    CALL MAPL_GetResource(STATE,x2Divisor,'X2_DIVISOR:',DEFAULT=9.27, __RC__)
    CALL MAPL_GetResource(STATE,x3Divisor,'X3_DIVISOR:',DEFAULT=34.4, __RC__)
    CALL MAPL_GetResource(STATE,x4Divisor,'X4_DIVISOR:',DEFAULT=21.4, __RC__)
    CALL MAPL_GetResource(STATE,x5Divisor,'X5_DIVISOR:',DEFAULT=14600., __RC__)


    ! Exponent for the surface temperature deviation predictor
    ! --------------------------------------------------------
    CALL MAPL_GetResource(STATE,x5Power,'X5_EXPONENT:',DEFAULT=3.00, __RC__)

    ! Threshold temperatures
    ! ----------------------
    CALL MAPL_GetResource(STATE,sfcTLimit,'SFC_T_LIMIT:',DEFAULT=273.0, __RC__)
    CALL MAPL_GetResource(STATE,airTLimit,'AIR_T_LIMIT:',DEFAULT=263.0, __RC__)
    
    ! Cloud-top pressure limiter
    ! --------------------------
    CALL MAPL_GetResource(STATE,hPaCldTop,'CLOUD_TOP_LIMIT:',DEFAULT=500., __RC__)

    ! Layer depths [m]
    ! ----------------
    ALLOCATE(dZ(nc,lm), __STAT__)
    dZ = zle(:,0:lm-1)-zle(:,1:lm)


    ! Pressure at mid-layer [Pa]
    ! --------------------------
    ALLOCATE(p(nc,lm), __STAT__)
    p = (ple(:,1:lm)+ple(:,0:lm-1))*0.50


    ! Temperature at mid-layer [K]
      ! ----------------------------
    ALLOCATE(T(nc,lm), __STAT__)
    T = TH*((p*(1.0/MAPL_P00))**(MAPL_RGAS/MAPL_CP))

    ! Reset CNV_TOPP's MAPL_UNDEFs to zeroes
    ! --------------------------------------
    ALLOCATE(cnv_topp(nc), __STAT__)
    WHERE(CCTP == MAPL_UNDEF)
       cnv_topp = 0.00
    ELSEWHERE
       cnv_topp = CCTP
    END WHERE


    ! Set weak/no convection mask
    ! ---------------------------
    ALLOCATE(weakCnvMask(nc), __STAT__)
    weakCnvMask = 0
    WHERE(cn_prcp == 0.00 .OR. cnv_topp >= hPaCldTop*100.00 .OR. CAPE >= MAPL_UNDEF) weakCnvMask = 1

    ! Convective cloud top mask
    ! -------------------------
    ALLOCATE(cloudTopMask(nc,lm), __STAT__)
    cloudTopMask = 0
    DO k = 1,lm
       WHERE(ple(1:nc,k) > cnv_topp(1:nc) .AND. cnv_topp(1:nc) > 0.00) cloudTopMask(1:nc,k) = 1
    END DO

    ! Cloud top distance above ground [m]
    ! -----------------------------------
    ALLOCATE(cloudTopAG(nc), __STAT__)
    cloudTopAG = 0.00
    DO i = 1,nc
       n = SUM(cloudTopMask(i,1:lm))
       IF(n > 0) cloudTopAG(i) = SUM(dZ(i,lm-n+1:lm))
    END DO


    ! X1: Cold cloud depth: Vertical extent [km] where T < airTLimit and p > cnv_topp
    ! -------------------------------------------------------------------------------
    ALLOCATE(x1(nc), __STAT__)
    ALLOCATE(mask(nc,lm), __STAT__)

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

! print*,'MOIST - Count of weakCnvMask==0: ', COUNT(weakCnvMask == 0)
    
    ! X4: Integrated convective mass flux
    ! -----------------------------------
    ALLOCATE(x4(nc), __STAT__)
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
    ALLOCATE(x5(nc), __STAT__)
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
    ALLOCATE(x2(nc), __STAT__)
    x2 = cloudTopAG*0.001
    WHERE(weakCnvMask == 1) x2 = 0.00
    x2 = x2/x2Divisor


    ! X3: CAPE
    ! --------
    ALLOCATE(x3(nc), __STAT__)
    x3 = CAPE
    WHERE(weakCnvMask == 1) x3 = 0.00
    x3 = x3/x3Divisor

    ! Polynomial fit [units: km^{-2} s^{-1}] and individual
    ! terms including marine and continental discrimination
    ! -----------------------------------------------------
    WHERE(frOcean >= 0.01)
       strokeRate = (a0m + a1m*x1 + a2m*x2 + a3m*x3 + a4m*x4 + a5m*x5)/SEC_PER_DAY
       A1X1 = a1m*x1/SEC_PER_DAY
       A2X2 = a2m*x2/SEC_PER_DAY
       A3X3 = a3m*x3/SEC_PER_DAY
       A4X4 = a4m*x4/SEC_PER_DAY
       A5X5 = a5m*x5/SEC_PER_DAY
    ELSEWHERE
       strokeRate = (a0c + a1c*x1 + a2c*x2 + a3c*x3 + a4c*x4 + a5c*x5)/SEC_PER_DAY
       A1X1 = a1c*x1/SEC_PER_DAY
       A2X2 = a2c*x2/SEC_PER_DAY
       A3X3 = a3c*x3/SEC_PER_DAY
       A4X4 = a4c*x4/SEC_PER_DAY
       A5X5 = a5c*x5/SEC_PER_DAY
    END WHERE


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
!    print*, "MOIST_FlashRate: A1X1-5(1): ", A1X1(1), A2X2(1), A3X3(1), A4X4(1), A5X5(1)

    strokeRate = strokeRate * MOIST_flashFactor

    DEALLOCATE(x1, __STAT__)
    DEALLOCATE(x2, __STAT__)
    DEALLOCATE(x3, __STAT__)
    DEALLOCATE(x4, __STAT__)
    DEALLOCATE(x5, __STAT__)
    DEALLOCATE(cnv_topp, __STAT__)
    DEALLOCATE(dZ, __STAT__)
    DEALLOCATE(p, __STAT__)
    DEALLOCATE(T, __STAT__)
    DEALLOCATE(cloudTopAG, __STAT__)
    DEALLOCATE(mask, __STAT__)
    DEALLOCATE(cloudTopMask, __STAT__)
    DEALLOCATE(weakCnvMask, __STAT__)
    
    return

  end subroutine MOIST_FlashRate

!-------------------------------------------------------------------------
!EOP


!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: FIT_FlashRate   (previously called DALE_FlashRate)

!
! !INTERFACE:
!
    subroutine FIT_FlashRate (cldmas, threshold, ratio_local, ratio_global, midlatAdj, &
         desired_g_N_prod_rate, FIT_flashFactor, flashrate)

! !USES
      implicit none

! !INPUT PARAMETERS:

      !     Convective mass flux at desired layer (kg m-2 min-1)
      !     Warning: Make sure units are correct and that you've accessed the proper vertical layer. 
      real, intent(in) :: cldmas(:, :)    
      
      !     Desired mean lightning NO production rate (Tg yr-1) (specified in namelist)       	
      real*8, intent(in) :: desired_g_N_prod_rate
      
      !     Mass flux threshold (kg m-2 min-1) below which flash rate is assumed to be zero (Read in from file) 
      real*8, intent(in) :: threshold            
      
      !     Adjustment factor local flash rates must be multiplied so that monthly average local flash rates (after
      !     "ratio_global" adjustment) match monthly average local v2.2 climatological OTD/LIS flash rates (Read in from file).  
      !     Note: ratio_local varies monthly. 
      !real*8, intent(in) :: ratio_local(:,:) 
      real, pointer, intent(in) :: ratio_local(:,:) 

      !real*8, intent(in) :: midLatAdj(:,:) 
      real, pointer, intent(in) :: midLatAdj(:,:) 
              
      !     Adjustment factor local flash rates must be multiplied so that globally averaged 
      !     flash rate matches v2.2 OTD/LIS climatological globally averaged flash rate (Read in from file) 
      !     Note: ratio_global varies monthly. 
      real, intent(in) :: ratio_global             

      !     Global scaling term
      real, intent(in) :: FIT_flashFactor             
            
      
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

      flashrate = flashrate * FIT_flashFactor

      ! Calculate N production rate (g s-1).
      ! For specified flashrate and PRODFAC, global production rate = 3.41 Tg N yr-1.
      ! Assume IC and CG flash rates produce the same amount of N / flash.
      !      pnox = flashrate(:,:)*PRODFAC*CONVFAC
      ! pnox = midLatAdj(:,:)*flashrate(:,:)*PRODFAC*CONVFAC 
      
      !adjust for desired global Nitrogen production rate
      !pnox = pnox * desired_g_N_prod_rate /3.41        
    
      deallocate(cldmas_local)
      return

    end subroutine FIT_FlashRate
!EOP

!==========================================================================

!BOP

! !IROUTINE: HEMCO_FlashRate

! !DESCRIPTION: Wrapper routine for GEOS Chem's flashrate routine.
!
! !INTERFACE:
!
    subroutine HEMCO_FlashRate (cellArea, lwi, lonslocal, latslocal, &
         airTemp, ple, geoPotHeight, cnvMfc, otdLisScale, flashRate, RC )

! !USES
      implicit none

! !ARGUMENTS

    !! IN
      !     Grid box cell area (m2)
      real, intent(in) :: cellArea(:,:)    
      
      !     Land water ice flag  (0=water; 1=land; 2=ice)
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

!     ! Buoyancy (m/s/s) - not needed when LCNVFRC == FALSE
!     real, intent(in) :: buoyancy(:,:,:)

      ! otdLisScale = 0.355 for c48; 0.1 for c90
      real, intent(in) :: otdLisScale

      ! Longitude and latitude (local values) in degrees
      real, intent(in) :: lonslocal(:,:)
      real, intent(in) :: latslocal(:,:)

    !! OUT
      ! Flash rate in flashes km-2 s-1
      real, intent(out) :: flashRate(:,:)

      ! return code
      integer, optional, intent(out) :: RC

!EOP

! Locals

      integer :: STATUS
      character(len=*), parameter :: Iam = "HEMCO_FlashRate" ! can be moduleName // getLightning TODO
      
      real(hp) :: otdLisScaleHp

      real(hp), pointer, dimension(:,:)     ::       cellAreaHemco  => null()
      real(hp), pointer, dimension(:,:)     ::         otdLisHemco  => null()
      real(hp), pointer, dimension(:,:)     ::           lonsHemco  => null()
      real(hp), pointer, dimension(:,:)     ::           latsHemco  => null()
      real(hp), pointer, dimension(:,:)     ::       convFracHemco  => null()
      real(hp), pointer, dimension(:,:,:)   ::        airTempHemco  => null()
      real(hp), pointer, dimension(:,:,:)   ::            pleHemco  => null()
      real(hp), pointer, dimension(:,:,:)   ::  gridBoxHeightHemco  => null()
      real(hp), pointer, dimension(:,:,:)   ::  gridBoxHeightHemco2 => null()
      real(hp), pointer, dimension(:,:,:)   ::         cnvMfcHemco  => null()
      real(hp), pointer, dimension(:,:,:)   ::       buoyancyHemco  => null()

      integer,  pointer, dimension(:,:)     ::            lwiHemco  => null()
      
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

!EOP
!-------------------------------------------------------------------------


      shape3d = shape(airTemp)
      IM = shape3d(1)
      JM = shape3d(2)
      LM = shape3d(3)
      KM = LM + 1

      ALLOCATE(          pleHemco (1:IM, 1:JM, 1:LM+1), __STAT__ )
      ALLOCATE(       cnvMfcHemco (1:IM, 1:JM, 1:LM+1), __STAT__ )
      ALLOCATE(      airTempHemco (1:IM, 1:JM, 1:LM),   __STAT__ )
      ALLOCATE(     buoyancyHemco (1:IM, 1:JM, 1:LM),   __STAT__ )
      ALLOCATE(gridBoxHeightHemco (1:IM, 1:JM, 1:LM),   __STAT__ )
      ALLOCATE(gridBoxHeightHemco2(1:IM, 1:JM, 1:LM),   __STAT__ )
      ALLOCATE(     cellAreaHemco (1:IM, 1:JM),         __STAT__ )
      ALLOCATE(     convFracHemco (1:IM, 1:JM),         __STAT__ )
      ALLOCATE(          lwiHemco (1:IM, 1:JM),         __STAT__ )
      ALLOCATE(       otdLisHemco (1:IM, 1:JM),         __STAT__ )
      ALLOCATE(         lonsHemco (1:IM, 1:JM),         __STAT__ )
      ALLOCATE(         latsHemco (1:IM, 1:JM),         __STAT__ )

      ALLOCATE(lNox2d (1:IM, 1:JM),   __STAT__ )
      ALLOCATE( lfr2d (1:IM, 1:JM),   __STAT__ )
      ALLOCATE(  ic2d (1:IM, 1:JM),   __STAT__ )
      ALLOCATE(  cg2d (1:IM, 1:JM),   __STAT__ )
      ALLOCATE(  h02d (1:IM, 1:JM),   __STAT__ )
      ALLOCATE(ltop2d (1:IM, 1:JM),   __STAT__ )

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
           buoyancyHemco, convFracHemco, .False., otdLisScaleHp, lNox2d, lfr2d, ic2d, cg2d, ltop2d, h02d, __RC__)

      ! Convert to km-2 s-1
      flashRate (:,:) = real(lfr2d(:,:)) * MIN_PER_SEC

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

    end subroutine HEMCO_FlashRate


!-----------------------------------------------------------------------
!BOP
! !IROUTINE: LOPEZ_FlashRate
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


    subroutine LOPEZ_FlashRate(IM, JM, LM, FRLAND, ZKBCON, CAPE,  &
                               ZLE, PFI_CN_GF, CNV_QC, TE, PLO, LOPEZ_flashFactor, &
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

     real    ,intent(in)                   :: LOPEZ_flashFactor   ! global scaling term
                                                                  ! originally called 'alpha'

     real    ,intent(out), dimension(im,jm) ::  &
                              LFR_Lopez           ! lightning flash density rate (units: flashes/km2/day)
                         
    !-- locals
    real, parameter :: V_graup     = 3.0  ! m/s
    real, parameter :: V_snow      = 0.5  ! m/s
    real, parameter :: beta_land   = 0.70 ! 1
    real, parameter :: beta_ocean  = 0.45 ! 1
    real, parameter :: t_initial   =  0.0 + 273.15 ! K
    real, parameter :: t_final     = -25. + 273.15 ! K
    
    integer :: i,j,k, k_initial, k_final                          
    real    :: tdif, td2
    real    :: Q_R, z_base,beta,prec_flx_fr,dz
    real,    dimension(1:lm) :: q_graup,q_snow,rho


!!! DEBUGGING
!    print*, ""
!    print*, "Hello from LOPEZ_FlashRate"
!    print*, ""

!    print*, "im, jm, lm: ", im, jm, lm
!    print*, "FRLAND: ", minval(FRLAND), maxval(FRLAND)
!    print*, "CAPE: ", minval(CAPE), maxval(CAPE)
!    print*, "ZLE: ", minval(ZLE), maxval(ZLE)
!    print*, "PFI_CN_GF ", minval(PFI_CN_GF), maxval(PFI_CN_GF)
!    print*, "TE: ", minval(TE), maxval(TE)
!    print*, "PLO: ", minval(PLO), maxval(PLO)
!    print*, "CNV_QC: ", minval(CNV_QC), maxval(CNV_QC)
!    print*, "zkbcon: ", minval(zkbcon), maxval(zkbcon)
!!!

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

           k_initial = lm
           tdif =  ABS(t_initial - te(i,j,lm))
           do k=lm-1,1,-1
             td2 = ABS(t_initial - te(i,j,k))
             if ( td2 < tdif ) then
               k_initial = k
               tdif = td2
             endif
             if ( te(i,j,k) < t_initial ) exit
           enddo

           k_final = lm
           tdif =  ABS(t_final - te(i,j,lm))
           do k=lm-1,1,-1
             td2 = ABS(t_final - te(i,j,k))
             if ( td2 < tdif ) then
               k_final = k
               tdif = td2
             endif
             if ( te(i,j,k) < t_final ) exit
           enddo
! print*,'K final,  initial: ', k_final, k_initial

          Q_R = 0.0
          do k = k_final, k_initial
             
             dz  = zle(i,j,k-1)-zle(i,j,k)  

             Q_R = Q_R + dz*rho(k)*(q_graup(k)*(CNV_QC(i,j,k)+q_snow(k)))
         enddo

         z_base = ZKBCON(i,j)/1000. !- convert to [km]

        !--- lightning flash density (units: number of flashes/km2/day) - equation 5
        !--- (to compare with Lopez 2016's results, convert to per year: LFR_Lopez*365)
        !
        LFR_Lopez(i,j) = LOPEZ_flashFactor * Q_R *sqrt(max(0.,cape(i,j))) * min(z_base,1.8)**2 

     enddo
   enddo
        
        
if (maxval(LFR_Lopez)>1000.) then
       print*,"1================================================= LOPEZ"
       print*," LFR=",maxval(LFR_Lopez),maxval(cape),maxval(CNV_QC),maxval(PFI_CN_GF)
       print*,"2================================================= LOPEZ"
       call flush(6)
endif        
  end subroutine LOPEZ_FlashRate
!========================================================================================

!-----------------------------------------------------------------------------
! !ROUTINE
!   emiss_lightning
!
! !DESCRIPTION
!   
!  Generate NOx production rates from parameterized lightning and distribute vertically
!  using profiles from Pickering (2007).
!
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
                            lwi, flashrate, cellDepth, dtrn, pNOx3D, rc, kgNOx3D)

  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: i1, i2, j1, j2, k1, k2         ! Index ranges on this processor
  REAL,    INTENT(IN)  :: minDeepCloudTop                ! Minimum cloud top [km] for selecting deep convection profiles
  REAL,    INTENT(IN)  :: ampFactor                      ! > 0, for targeting the observed nitrogen production rate [3.41 Tg yr^{-1}]
  REAL,    INTENT(IN)  :: numberNOperFlash               ! NO molecules generated by each flash


  INTEGER, INTENT(IN)  :: lwi(i1:i2, j1:j2)              ! Flag: 0=water 1=land 2=ice

  REAL*8,  INTENT(IN)  :: flashrate(i1:i2, j1:j2)        ! Flash rate [km^{-2} s^{-1}]
  REAL*4,  INTENT(IN)  :: cellDepth(i1:i2, j1:j2, k1:k2) ! Grid cell depth [m]             (bottom-up)
  REAL*4,  INTENT(IN)  ::      dtrn(i1:i2, j1:j2, k1:k2) ! Detrainment [kg m^{-2} s^{-1}]  (bottom-up)
  
  REAL*4,  INTENT(INOUT), POINTER  :: pNOx3D(:,:,:)      ! Lightning NO production rate [m^{-3} s^{-1}]  (bottom-up)
                                                         ! Must be ASSOCIATED
                                                         ! (i1:i2, j1:j2, k1:k2)

  REAL,    INTENT(OUT), OPTIONAL :: kgNOx3D(i1:i2, j1:j2, k1:k2)   ! NO production rate [kg m^{-3} s^{-1}]  (bottom-up)
  INTEGER, INTENT(OUT)  :: rc                             ! Return code - TODO: make this optional!

! Local
! -----
  INTEGER :: k
  INTEGER :: status
  REAL*8, ALLOCATABLE  :: pNOx2D(:,:)                    ! Lightning NO production [molecules NO m^{-2} s^{-1}]
  CHARACTER(LEN=*), PARAMETER :: Iam = "emiss_lightning"
  rc = 0
  status = 0

  _ASSERT( LBOUND(pNOx3D, 1) == i1, 'bad LBOUND 1 for pNOX3D' )
  _ASSERT( UBOUND(pNOx3D, 1) == i2, 'bad UBOUND 1 for pNOX3D' )
  _ASSERT( LBOUND(pNOx3D, 2) == j1, 'bad LBOUND 2 for pNOX3D' )
  _ASSERT( UBOUND(pNOx3D, 2) == j2, 'bad UBOUND 2 for pNOX3D' )
  _ASSERT( LBOUND(pNOx3D, 3) == k1, 'bad LBOUND 3 for pNOX3D' )
  _ASSERT( UBOUND(pNOx3D, 3) == k2, 'bad UBOUND 3 for pNOX3D' )

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
  CALL partition(i1, i2, j1, j2, k1, k2, minDeepCloudTop, lwi, pNOx2D, dtrn, cellDepth, pNOx3D, rc)

! Place output in useful units
! ----------------------------

! Number density tendency [m^{-3} s^{-1}]
! ---------------------------------------
  pNOx3D(:,:,:) = pNOx3D(:,:,:)/cellDepth(:,:,:)

! NO density tendency [kg NO m^{-3} s^{-1}]
! -----------------------------------------
  IF ( PRESENT(kgNOx3D)) &
               kgNOx3D(:,:,:) = pNOx3D(:,:,:)*MWT_NO/MAPL_AVOGAD

! Clean up
! --------
  DEALLOCATE(pNOx2D,STAT=status)
  VERIFY_(status)

  RETURN
 END SUBROUTINE emiss_lightning

!-----------------------------------------------------------------------------
!BOP
!
! !ROUTINE
!   partition
!
! !DESCRIPTION
!  Vertically distribute NOx
!
!  This routine came from GmiEmiss_lightning_mod in the GEOS-CTM
!-----------------------------------------------------------------------------
 SUBROUTINE partition(i1, i2, j1, j2, k1, k2, minDeepCloudTop, lwi, pNOx2D, dtrn, cellDepth, pNOx3D, rc)

  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: i1, i2, j1, j2, k1, k2    ! Index ranges on this processor
  REAL,    INTENT(IN)  :: minDeepCloudTop           ! Minimum cloud top [km] for selecting deep convection profiles


  INTEGER, INTENT(IN)  ::       lwi(i1:i2, j1:j2)        ! Flag: 0=water 1=land 2=ice

  REAL*8,  INTENT(IN)  ::    pNOx2D(i1:i2, j1:j2)        ! Lightning NO production rate [molecules NO m^{-2} s^{-1}]
  REAL*4,  INTENT(IN)  ::      dtrn(i1:i2, j1:j2, k1:k2) ! Detrainment [kg m^{-2}s^{-1}]    (bottom-up)
  REAL*4,  INTENT(IN)  :: cellDepth(i1:i2, j1:j2, k1:k2) ! Grid cell depth [m]              (bottom-up)

  REAL*4,  INTENT(INOUT), POINTER :: pNOx3D(:,:,:)       ! Scaled production rate (no units conversion here)    (bottom-up)
                                                         ! (i1:i2, j1:j2, k1:k2)

! Local variables
! ---------------
  CHARACTER(LEN=*), PARAMETER :: Iam = "partition"

  INTEGER :: i,j,k
  INTEGER :: cl                     ! vertical index
  INTEGER :: nTop                   ! Top model level at which detrainment is non-zero
  INTEGER :: profileNumber
  INTEGER :: rc, status
  INTEGER, PARAMETER :: numKm = 17  ! Number of elements (kilometers) in each specified profile

  REAL :: zLower, zUpper

  REAL, ALLOCATABLE :: r(:,:)       ! Specified NOx distribution profiles

  REAL ::       w(k1:k2)     ! Weights applied to scaled cloud layers
  REAL ::       z(k1:k2)     ! Height above ground for top edge of grid-box
  REAL :: zScaled(k1:k2)     ! Scaled layer edge heights
   
  rc = 0
  status = 0

  _ASSERT( LBOUND(pNOx3D, 1) == i1, 'bad LBOUND 1 for pNOx3D' )
  _ASSERT( UBOUND(pNOx3D, 1) == i2, 'bad UBOUND 1 for pNOx3D' )
  _ASSERT( LBOUND(pNOx3D, 2) == j1, 'bad LBOUND 2 for pNOx3D' )
  _ASSERT( UBOUND(pNOx3D, 2) == j2, 'bad UBOUND 2 for pNOx3D' )
  _ASSERT( LBOUND(pNOx3D, 3) == k1, 'bad LBOUND 3 for pNOx3D' )
  _ASSERT( UBOUND(pNOx3D, 3) == k2, 'bad UBOUND 3 for pNOx3D' )

! Specify the percentage NOx distributions in each km for a numKm-depth cloud.
! Deep convection is arbitrarily assigned when the cloud top is greater than 7 km.
! --------------------------------------------------------------------------------
   ALLOCATE(r(numKm,3),STAT=status)
   VERIFY_(status)

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

       IF(lwi(i,j) == 1      ) THEN                 !  LAND
        profileNumber = 1
!  PRINT*,'PROFILE_1  MAX WGT ', z(nTop)*13./17.
       ELSE
        profileNumber = 2
!  PRINT*,'PROFILE_2  MAX WGT ', z(nTop)*11./17.
       END IF
      ELSE
       profileNumber = 3
!  PRINT*,'PROFILE_3  MAX WGT ', z(nTop)* 7./17.
      ENDIF

! Scale factor, 0 at ground, numKm at cloud top
! The distance from ground through cloud top is artificially
! scaled to be 0 through 17 km.
!
! The z quotient is <  1 below nTop,  1 at nTop, >  1 above nTop
! Scale factor   is < 17 below nTop, 17 at nTop, > 17 above nTop
! zScaled(i) is the top of the model level i (km)
! ---------------------------------------------
    ! zScaled(1:k2) = numKm*(z(1:k2)/z(nTop))
      zScaled(1:k2) = z(1:k2)*numKm/z(nTop)

! Intialize
! ---------
      w(:) = 0.00    ! weights - will only use indices  1:nTop
      cl = 1         ! model level index - start at the bottom
      zLower = 0.00  ! edge of gridbox or kilometer boundary (km)

! Compute the weight (w) to be applied at each level
! --------------------------------------------------
      Kilometers: DO k = 1,numKm      ! k = index into profile
                                      ! k is also the height of top edge (km)

! ... segment-by-segment
! A segment extends from zLower to zUpper
! Each of those can be the edge of a gridbox or a kilometer boundary
! ------------------------------------------------------------------
       Segment: DO

! Push up to the lesser of scaled cloud height or next km
! -------------------------------------------------------
        zUpper = MIN(zScaled(cl),k*1.)
        IF(zScaled(cl) > numKm) EXIT

! Add increment to scaled weighting for the current cloud layer
! Accumulate the weight
! Convert profile value (r) from 0-100 range, to 0-1 range
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
        IF(zUpper == k*1.) EXIT

       END DO Segment

      END DO Kilometers

! Finalize vertical distribution and clean up
! -------------------------------------------
      pNOx3D(i,j,1:nTop) = REAL( w(1:nTop)*pNOx2D(i,j) )

! PRINT*,'TOTAL for w is ', SUM(w(1:nTop))  ->  this is 1.0

     END IF SeeingLightning

! Next column
! -----------
    END DO
   END DO

! Clean up
! --------
   DEALLOCATE(r,STAT=status)
   VERIFY_(status)

  RETURN
 END SUBROUTINE partition

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


    subroutine identify_flash_source( flash_source_str, flash_source_enum )

      character(len=ESMF_MAXSTR), intent(in)  :: flash_source_str
      integer,                    intent(out) :: flash_source_enum

      integer :: i

      flash_source_enum = FLASH_SOURCE_UNDEFINED
      do i= 1, FLASH_SOURCE_count
        if ( TRIM(flash_source_str) == TRIM(flashSourceNames(i)) ) flash_source_enum = i
      end do

!     print*, "IN identify_flash_source, flash source: ", TRIM(flashSourceNames(flash_source_enum))
      return

    end subroutine identify_flash_source


    ! Only needed for FIT
    subroutine update_lightning_ratio( ratioGlobalFile, CLOCK, year, month, ratioGlobalLight, RC )
      character(len=ESMF_MAXSTR),    intent(in)       :: ratioGlobalFile
      type (ESMF_Clock),             intent(inout)    :: CLOCK 
      integer,                       intent(inout)    :: year      ! in: previous year  processed;   out: year  from CLOCK
      integer,                       intent(inout)    :: month     ! in: previous month processed;   out: month from CLOCK
      real,                          intent(inout)    :: ratioGlobalLight  ! in: previous answer;  out: may be updated
      integer, optional,             intent(  out)    :: RC        ! Error code

    ! ErrLog Variables
      character(len=ESMF_MAXSTR)   :: IAm
      integer                      :: STATUS

      ! temp vars
      type(ESMF_Time) :: curr_date
      integer         :: mm, yy

      IAm = 'update_lightning_ratio'

      call ESMF_ClockGet(CLOCK, currTime=curr_date, __RC__ )
      call ESMF_TimeGet(curr_date,mm=mm,yy=yy,__RC__ )

      if ( year /= yy .OR. month /= mm ) then
        call readLightRatioGlobalData( CLOCK, ratioGlobalFile, ratioGlobalLight, MAPL_AM_I_ROOT(), __RC__ )
        _ASSERT(ratioGlobalLight > 0.0, 'ratioGlobalLight must be positive')
        year  = yy
        month = mm
      end if

    end subroutine update_lightning_ratio

      subroutine readLightRatioGlobalData &
           &  (CLOCK, light_ratioGlobal_infile_name, ratioGlobalLight, rootProc, RC)


        type (ESMF_Clock),             intent(inout)    :: CLOCK 
        character (len=*) :: light_ratioGlobal_infile_name
        real, intent(out) :: ratioGlobalLight
        logical, intent(in) :: rootProc
        integer, optional, intent(out)    :: RC        ! Error code

      ! ErrLog Variables
        character(len=ESMF_MAXSTR)   :: IAm
        integer                      :: STATUS

        integer :: ii, mm, asc_lun, ierr
        real    :: readRatio

        integer :: mon, year, nym
        type(ESMF_Time) :: curr_date             ! temporary date used in logic
        character(len=80) :: cYear
        character(len=80) :: cMonth
        character(len=80) :: cNym
        character(len=3)  :: binName


        character (len=256) :: err_msg

        integer :: num_reads = 12 * 37 ! roughly 1 entry for each month in satellite era

        IAm = "readLightRatioGlobalData"

        call ESMF_ClockGet(CLOCK, currTime=curr_date, __RC__ )
        call ESMF_TimeGet(curr_date,mm=mon,yy=year,__RC__ )

        write(unit=cYear,fmt=*) year
        write(unit=cMonth,fmt=*) mon

        if (mon >= 10) then
           cNym = trim(adjustl(cYear))//trim(adjustl(cMonth))
        else
           cNym = trim(adjustl(cYear))//'0'//trim(adjustl(cMonth))
        endif
        read(unit=cNym,fmt=*), nym

        asc_lun = 9

        OPEN(UNIT=asc_lun, FILE=TRIM(light_ratioGlobal_infile_name), STATUS='old', ACTION='read', &
             FORM='formatted', ACCESS='sequential', IOSTAT=ierr)

        if (ierr /= 0) then
           err_msg = 'Failed to OPEN '//TRIM(light_ratioGlobal_infile_name)
           STATUS=99
           VERIFY_(STATUS)
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
         WRITE(6,*) 'Ratio global for lightning ', year, mon, ratioGlobalLight
      END IF

      return

    end subroutine readLightRatioGlobalData


!========================================================================================



      
  end module Lightning_mod
