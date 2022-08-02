#include "MAPL_Generic.h"

!=============================================================================
!BOP

! !MODULE: GEOS_ChemEnvGridCompMod -- Prepares Environment for GEOSchem

! !INTERFACE:

module GEOS_ChemEnvGridCompMod

! !USES:

  use ESMF
  use MAPL

  use OVP,          only:  OVP_init, OVP_end_of_timestep_hms, OVP_mask, OVP_apply_mask
  use Lightning_mod

  implicit none
  private

  INTEGER, SAVE, ALLOCATABLE :: MASK_10AM(:,:)
  INTEGER, SAVE, ALLOCATABLE :: MASK_2PM(:,:)
  INTEGER, SAVE              :: OVP_FIRST_HMS
  INTEGER, SAVE              :: OVP_RUN_DT
  INTEGER, SAVE              :: OVP_GC_DT
  INTEGER, SAVE              :: OVP_MASK_DT
  LOGICAL                    :: OVP_setup_done

  integer,  parameter        :: r8 = 8
  character(len=ESMF_MAXSTR) :: rcfilen = 'ChemEnv.rc'

  ! Retrieve from rc file:
  integer                    :: flash_source_enum = FLASH_SOURCE_UNDEFINED
  character(len=ESMF_MAXSTR) :: ratioGlobalFile     ! only needed for FIT
  real                       :: minDeepCloudTop     ! needed for emiss_lightning
  real                       :: lightNOampFactor    ! needed for emiss_lightning
  real                       :: numberNOperFlash    ! needed for emiss_lightning
  real                       :: MOIST_flashFactor
  real                       :: FIT_flashFactor
  real                       :: HEMCO_flashFactor
  real                       :: LOPEZ_flashFactor

  ! May change during the course of the run:
  integer                    :: year_for_ratio = 0
  integer                    :: month_for_ratio = 0
  real                       :: ratioGlobalLight    ! only needed for FIT

! !PUBLIC MEMBER FUNCTIONS:

  public SetServices

! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE                    :: Initialize_    ! Init method
  PRIVATE                    :: Run1           ! Run method
  PRIVATE                    :: Run2           ! Run method

  PRIVATE                    :: Airdens        ! Compute Air Density

! !DESCRIPTION: This is a Cinderella gridded component (GC) 
!EOP

contains

!BOP

! !IROUTINE: SetServices -- Sets ESMF services for this component

! !INTERFACE:

    subroutine SetServices ( GC, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
    integer, intent(OUT)               :: RC  ! return code

! !DESCRIPTION:  The SetServices for the Chemistry Env GC needs to register its
!   Initialize and Run.  It uses the MAPL\_Generic construct for defining 
!   state specs and couplings among its children.  In addition, it creates the   
!   children GCs and runs their respective SetServices.

!EOP

!=============================================================================
!
! ErrLog Variables

    character(len=ESMF_MAXSTR)              :: IAm
    integer                                 :: STATUS
    character(len=ESMF_MAXSTR)              :: COMP_NAME

! Locals

    type (ESMF_Config)         :: CF

!=============================================================================

! Begin...

! Get my name and set-up traceback handle
! ---------------------------------------

    Iam = 'ChemEnv::SetServices'
    call ESMF_GridCompGet( GC, NAME=COMP_NAME, CONFIG=CF, __RC__ )
    Iam = trim(COMP_NAME) // 'SetServices'

    OVP_setup_done = .FALSE.

! Register services for this component
! ------------------------------------
    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_INITIALIZE, Initialize_, __RC__ )
    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,        Run1,        __RC__ ) 
    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,        Run2,        __RC__ ) 

!BOS

    call read_flash_source ( rcfilen, flash_source_enum, __RC__ )

! !IMPORT STATE:

    call MAPL_AddImportSpec(GC,                                    &
         SHORT_NAME = 'PLE',                                       &
         LONG_NAME  = 'air_pressure',                              &
         UNITS      = 'Pa',                                        &
         DIMS       =  MAPL_DimsHorzVert,                          &
         VLOCATION  =  MAPL_VLocationEdge,                  __RC__ )

    call MAPL_AddImportSpec(GC,                                    &
         SHORT_NAME = 'TH',                                        &
         LONG_NAME  = 'potential_temperature',                     &
         UNITS      = 'K',                                         &
         DIMS       =  MAPL_DimsHorzVert,                          &
         VLOCATION  =  MAPL_VLocationCenter,                __RC__ )

     call MAPL_AddImportSpec(GC,                                   &
        SHORT_NAME = 'Q',                                          &
        LONG_NAME  = 'specific_humidity',                          &
        UNITS      = 'kg kg-1',                                    &
        DIMS       = MAPL_DimsHorzVert,                            &
        VLOCATION  = MAPL_VLocationCenter,                  __RC__ )

!   Geopotential Height
!   -------------------
    call MAPL_AddImportSpec(GC,                                    &
         SHORT_NAME='PHIS',                                        &
         LONG_NAME ='Geopotential Height at Surface',              &
         UNITS     ='m+2 s-2',                                     &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                    __RC__ )

!   Convective precip
!   -----------------
    call MAPL_AddImportSpec(GC,                                    &
         SHORT_NAME='CN_PRCP',                                     &
         LONG_NAME ='Surface Conv. rain flux needed by land',      &
         UNITS     ='kg m-2 s-1',                                  &
         DEFAULT   = MAPL_UNDEF,                                   &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                    __RC__ )

!   Total precip
!   ------------
    call MAPL_AddImportSpec(GC,                                    &
         SHORT_NAME='TPREC',                                       &
         LONG_NAME ='total_precipitation',                         &
         UNITS     ='kg m-2 s-1',                                  &
         DEFAULT   = MAPL_UNDEF,                                   &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                    __RC__ )

!   call MAPL_AddImportSpec(GC,                                    &
!      SHORT_NAME  = 'FRLANDICE',                                  &
!      LONG_NAME  = 'fraction_of_land_ice',                        &
!      UNITS  = '1',                                               &
!      DIMS  = MAPL_DimsHorzOnly,                                  &
!      VLOCATION  = MAPL_VLocationNone,                            &
!         RC=STATUS  )
!   VERIFY_(STATUS)

!!  NEEDED?:
!    call MAPL_AddImportSpec(GC,				   &
!         SHORT_NAME	  = 'FRACI',  			           &
!         LONG_NAME          = 'ice_covered_fraction_of_tile',     &
!         UNITS              = '1',                                &
!         DIMS               = MAPL_DimsHorzOnly,                  &
!         VLOCATION          = MAPL_VLocationNone,                 &
!    						          __RC__ )

    if (flash_source_enum == FLASH_SOURCE_FIT) then
      call MAPL_AddImportSpec(GC,                                    &
               SHORT_NAME  = 'MCOR',                                 &
               LONG_NAME   = 'agrid_cell_area',                      &
               UNITS       = 'm+2',                                  &
               DIMS        = MAPL_DimsHorzOnly,                      &
               VLOCATION   = MAPL_VLocationNone,              __RC__ )

      call MAPL_AddImportSpec(GC,                                    &
               SHORT_NAME  = 'RATIO_LOCAL',                          &
               LONG_NAME   = 'local_ratios_lightning',               &
               UNITS       = '1',                                    &
               DIMS        = MAPL_DimsHorzOnly,                      &
               VLOCATION   = MAPL_VLocationNone,              __RC__ )
    
      call MAPL_AddImportSpec(GC,                                    &
               SHORT_NAME  = 'MIDLAT_ADJ',                           &
               LONG_NAME   = 'midlat_adjustment_lightning',          &
               UNITS       = '1',                                    &
               DIMS        = MAPL_DimsHorzOnly,                      &
               VLOCATION   = MAPL_VLocationNone,              __RC__ )
    endif

    call MAPL_AddImportSpec(GC,                                    &
             SHORT_NAME = 'CNV_MFC',                               &
             LONG_NAME  = 'cumulative_mass_flux',                  &
             UNITS      = 'kg m-2 s-1',                            &
             DIMS       = MAPL_DimsHorzVert,                       &
             VLOCATION  = MAPL_VLocationEdge,               __RC__ )

    call MAPL_AddImportSpec(GC,                                    &
             SHORT_NAME = 'CNV_MFD',                               &
             LONG_NAME  = 'detraining_mass_flux',                  &
             UNITS      = 'kg m-2 s-1',                            &
             DIMS       = MAPL_DimsHorzVert,                       &
             VLOCATION  = MAPL_VLocationCenter,             __RC__ )

    call MAPL_AddImportSpec(GC,                                    &
             SHORT_NAME = 'PFI_CN',                                &
             LONG_NAME  = 'ice_convective_precipitation',          &
             UNITS      = 'kg m-2 s-1',                            &
             DIMS       = MAPL_DimsHorzVert,                       &
             VLOCATION  = MAPL_VLocationEdge,               __RC__ )

    call MAPL_AddImportSpec ( gc,                                  &
             SHORT_NAME = 'T',                                     &
             LONG_NAME  = 'air_temperature',                       &
             UNITS      = 'K',                                     &
             DIMS       = MAPL_DimsHorzVert,                       &
             VLOCATION  = MAPL_VLocationCenter,             __RC__ )

    call MAPL_AddImportSpec(GC,                                    &
             SHORT_NAME = 'TS',                                    &
             LONG_NAME  = 'surface_skin_temperature',              &
             UNITS      = 'K',                                     &
             DIMS       = MAPL_DimsHorzOnly,                       &
             VLOCATION  = MAPL_VLocationNone,               __RC__ )

    call MAPL_AddImportSpec(GC,                                    &
             SHORT_NAME = 'FROCEAN',                               &
             LONG_NAME  = 'fraction_of_ocean',                     &
             UNITS      = '1',                                     &
             DIMS       = MAPL_DimsHorzOnly,                       &
             VLOCATION  = MAPL_VLocationNone,               __RC__ )

    call MAPL_AddImportSpec(GC,                                    &
             SHORT_NAME = 'FRLAND',                                &
             LONG_NAME  = 'fraction_of_land',                      &
             UNITS      = '1',                                     &
             DIMS       = MAPL_DimsHorzOnly,                       &
             VLOCATION  = MAPL_VLocationNone,               __RC__ )

    call MAPL_AddImportSpec(GC,                                    &
             SHORT_NAME         = 'ZLE',                           &
             LONG_NAME          = 'geopotential_height',           &
             UNITS              = 'm',                             &
             DIMS               = MAPL_DimsHorzVert,               &
             VLOCATION          = MAPL_VLocationEdge,       __RC__ )

    call MAPL_AddImportSpec(GC,                                    &
             SHORT_NAME         = 'LWI',                           &
             LONG_NAME          = 'land(1)_water(0)_ice(2)_flag',  &
             UNITS              = '1',                             &
             DIMS               = MAPL_DimsHorzOnly,               &
             VLOCATION          = MAPL_VLocationNone,       __RC__ )

    call MAPL_AddImportSpec(GC,                                    &
             SHORT_NAME     = 'ZPBL',                              &
             LONG_NAME      = 'planetary_boundary_layer_height',   &
             UNITS          = 'm',                                 &
             DIMS           = MAPL_DimsHorzOnly,                   &
             VLOCATION      = MAPL_VLocationNone,           __RC__ )

    call MAPL_AddImportSpec(GC,                                    &
             SHORT_NAME        = 'AREA',                           &
             LONG_NAME         = 'agrid_cell_area',                &
             UNITS             = 'm+2',                            &
             DIMS              = MAPL_DimsHorzOnly,                &
             VLOCATION         = MAPL_VLocationNone,        __RC__ )

    call MAPL_AddImportSpec(GC,                                    &
             SHORT_NAME = 'CNV_QC',                                &
             LONG_NAME  = 'grid_mean_convective_condensate',       &
             UNITS      = 'kg kg-1',                               &
             DIMS       = MAPL_DimsHorzVert,                       &
             VLOCATION  = MAPL_VLocationCenter,             __RC__ )



! !EXPORT STATE:

!    AIRDENS: Provided for Children
!    ------------------------------
     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'AIRDENS',                      &
        LONG_NAME          = 'moist_air_density',            &
        UNITS              = 'kg m-3',                       &
        DIMS               = MAPL_DimsHorzVert,              &
        VLOCATION          = MAPL_VLocationCenter,    __RC__ )

!    Density of dry air
!    ------------------
     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'AIRDENS_DRYP',                 &
        LONG_NAME          = 'partial_dry_air_density',      &
        UNITS              = 'kg dry m-3 tot',               &
        DIMS               = MAPL_DimsHorzVert,              &
        VLOCATION          = MAPL_VLocationCenter,    __RC__ )

!    DELP (This should be wired from DYN)
!    ------------------------------------
     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'DELP',                         &
        LONG_NAME          = 'pressure_thickness',           &
        UNITS              = 'Pa',                           &
        DIMS               = MAPL_DimsHorzVert,              &
        VLOCATION          = MAPL_VLocationCenter,    __RC__ )

!   Total precip
!   ------------
    call MAPL_AddExportSpec(GC,                              &
         SHORT_NAME        = 'TPREC',                        &
         LONG_NAME         = 'total_precipitation',          &
         UNITS             = 'kg m-2 s-1',                   &
         DIMS              = MAPL_DimsHorzOnly,              &
         VLOCATION         = MAPL_VLocationNone,      __RC__ )

!    Convective precip
!    -----------------
     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'CN_PRCP',                      &
        LONG_NAME          = 'Convective precipitation',     &
        UNITS              = 'kg m-2 s-1',                   &
        DIMS               = MAPL_DimsHorzOnly,              &
        VLOCATION          = MAPL_VLocationNone,      __RC__ )

!    Non-convective precip
!    ---------------------
     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'NCN_PRCP',                     &
        LONG_NAME          = 'Non-convective precipitation', &
        UNITS              = 'kg m-2 s-1',                   &
        DIMS               = MAPL_DimsHorzOnly,              &
        VLOCATION          = MAPL_VLocationNone,      __RC__ )

!    10am overpass AIRDENS
!    ---------------------
     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'OVP10_AIRDENS',                &
        LONG_NAME          = 'moist_air_density_10am_local', &
        UNITS              = 'kg m-3',                       &
        DIMS               = MAPL_DimsHorzVert,              &
        VLOCATION          = MAPL_VLocationCenter,  RC=STATUS)
     VERIFY_(STATUS)

!    2pm  overpass AIRDENS
!    ---------------------
     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'OVP14_AIRDENS',                &
        LONG_NAME          = 'moist_air_density_2pm_local',  &
        UNITS              = 'kg m-3',                       &
        DIMS               = MAPL_DimsHorzVert,              &
        VLOCATION          = MAPL_VLocationCenter,  RC=STATUS)
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'LFR',                          &
        LONG_NAME          = 'lightning_flash_rate',         &
        UNITS              = 'km-2 s-1',                     &
        DIMS               = MAPL_DimsHorzOnly,              &
        VLOCATION          = MAPL_VLocationNone,      __RC__ )

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'BYNCY',                        &
        LONG_NAME          = 'buoyancy_of surface_parcel',   &
        UNITS              = 'm s-2',                        &
        DIMS               =  MAPL_DimsHorzVert,             &
        VLOCATION          =  MAPL_VLocationCenter,   __RC__ )

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'CAPE',                         &
        LONG_NAME          = 'convective_avail_pot_energy',  &
        UNITS              = 'J m^{-2} <check this!>',       &
        DIMS               = MAPL_DimsHorzOnly,              &
        VLOCATION          = MAPL_VLocationNone,      __RC__ )

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'LIGHT_NO_PROD',                &
        LONG_NAME          = 'lightning_NO_prod_rate',       &
        UNITS              = 'm-3 s-1',                      &
        DIMS               =  MAPL_DimsHorzVert,             &
        PRECISION          =  ESMF_KIND_R4,                  &
        VLOCATION          =  MAPL_VLocationCenter,   __RC__ )


!EOS


! Create children's gridded components and invoke their SetServices
! -----------------------------------------------------------------
    call MAPL_GenericSetServices    ( GC, __RC__ )

    RETURN_(ESMF_SUCCESS)
  
  end subroutine SetServices


! !IROUTINE: Initialize_ -- Finish setting up

! !INTERFACE:

   subroutine Initialize_ ( GC, impChem, expChem, clock, RC )

! !USES:

   implicit NONE

! !INPUT PARAMETERS:

   type(ESMF_Clock),  intent(inout) :: clock      ! The clock

! !OUTPUT PARAMETERS:

   type(ESMF_GridComp), intent(inout)  :: GC      ! Grid Component
   type(ESMF_State), intent(inout) :: impChem     ! Import State
   type(ESMF_State), intent(inout) :: expChem     ! Export State
   integer, intent(out) ::  RC                    ! Error return code:
                                                  !  0 - all is well
                                                  !  1 - 

! !DESCRIPTION:  Read the resource file

!EOP

!=============================================================================
!
! ErrLog Variables

    character(len=ESMF_MAXSTR)              :: IAm
    integer                                 :: STATUS
    character(len=ESMF_MAXSTR)              :: COMP_NAME

! Locals

    type (ESMF_Grid)           :: esmfGrid
    integer                    :: dims(3), im, jm
    type (ESMF_Config)         :: CF
    character(len=ESMF_MAXSTR) :: rc_string
    character(len=10)          :: SimType      ! CTM, FREE or REPLAY

!=============================================================================

! Begin...

! Get my name and set-up traceback handle
! ---------------------------------------

    Iam = 'ChemEnv::Initialize_'
    call ESMF_GridCompGet( GC, NAME=COMP_NAME, GRID=esmfGrid, __RC__ )
    Iam = trim(COMP_NAME) // 'Initialize_'


!   Call GenericInitialize (for every Child?)
!   -----------------------------------------
    call MAPL_GenericInitialize ( GC, impChem, expChem, clock, __RC__ )

    call MAPL_GridGet ( esmfGrid, globalCellCountPerDim=dims, __RC__ )
    im = dims(1)
    jm = dims(2)

    ! SimType is added by ctm_setup
    call ESMF_GridCompGet ( GC, CONFIG=CF, __RC__ )
    call ESMF_ConfigGetAttribute ( CF, rc_string, Label="SimType:", Default="undefined", __RC__ )
1   IF(MAPL_AM_I_ROOT()) PRINT*,'The SimType resource is '//TRIM(rc_string)

    IF ( TRIM(rc_string) == "CTM" ) THEN
      SimType = "CTM"
    ELSE 
      call ESMF_ConfigGetAttribute ( CF, rc_string, Label="REPLAY_MODE:", Default="undefined", __RC__ )
!     IF(MAPL_AM_I_ROOT()) PRINT*,'The REPLAY_MODE is '//TRIM(rc_string)

      IF ( TRIM(rc_string) == "undefined" ) THEN
        SimType = "FREE"
      ELSE 
        SimType = "REPLAY"
      END IF 
    END IF 
!   IF(MAPL_AM_I_ROOT()) PRINT*,'The SimType is '//SimType


    call read_lightning_config ( im, jm, rcfilen, flash_source_enum,   &
                                 SimType, ratioGlobalFile,             &
                                 minDeepCloudTop, lightNOampFactor,    &
                                 numberNOperFlash,                     &
                                 MOIST_flashFactor, FIT_flashFactor,   &
                                 HEMCO_flashFactor, LOPEZ_flashFactor, &
                                 __RC__ )

    IF(MAPL_AM_I_ROOT()) THEN
      if ( flash_source_enum == FLASH_SOURCE_MOIST ) PRINT*,'MOIST_flashFactor is ',MOIST_flashFactor
      if ( flash_source_enum == FLASH_SOURCE_FIT   ) PRINT*, ' FIT_flashFactor is ',  FIT_flashFactor
      if ( flash_source_enum == FLASH_SOURCE_HEMCO ) PRINT*,'HEMCO_flashFactor is ',HEMCO_flashFactor
      if ( flash_source_enum == FLASH_SOURCE_LOPEZ ) PRINT*,'LOPEZ_flashFactor is ',LOPEZ_flashFactor
    ENDIF

    RETURN_(ESMF_SUCCESS)
  
  end subroutine Initialize_

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !IROUTINE: Run1 -- Phase 1 run method

! !INTERFACE:

  subroutine Run1 ( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:

  type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
  type(ESMF_State),    intent(inout) :: IMPORT ! Import state
  type(ESMF_State),    intent(inout) :: EXPORT ! Export state
  type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
  integer, optional,   intent(  out) :: RC     ! Error code

! !DESCRIPTION: The Run1 method of the Chemistry Environment Gridded 
! Component. It calculates air density used by chemistry and emissions.

!EOP

! ErrLog Variables
  character(len=ESMF_MAXSTR)           :: IAm
  integer                              :: STATUS
  character(len=ESMF_MAXSTR)           :: COMP_NAME

  type (MAPL_MetaComp), pointer        :: ggState
  integer                    :: ndt
  real(r8)                   :: DT

  real(r8), pointer          :: R8D2(:,:)
  real, pointer              :: LONSLOCAL(:,:), LATSLOCAL(:,:)
  type(ESMF_Grid)            :: esmfGrid

!!!
!!!

! Imports
  real, pointer, dimension(:,:,:) ::         PLE => null()
  real, pointer, dimension(:,:)   :: RATIO_LOCAL => null()
  real, pointer, dimension(:,:)   ::  MIDLAT_ADJ => null()

  real, pointer, dimension(:,:,:) ::           T => null()
  real, pointer, dimension(:,:,:) ::          TH => null()
  real, pointer, dimension(:,:,:) ::           Q => null()
  real, pointer, dimension(:,:,:) ::         ZLE => null()
  real, pointer, dimension(:,:)   ::    cellArea => null()
  real, pointer, dimension(:,:)   ::        mcor => null()
  real, pointer, dimension(:,:)   ::         LWI => null()
  real, pointer, dimension(:,:)   ::        PBLH => null()

  real, pointer, dimension(:,:)   ::          TS => null()
  real, pointer, dimension(:,:)   ::     FROCEAN => null()
  real, pointer, dimension(:,:)   ::      FRLAND => null()
  real, pointer, dimension(:,:)   ::     CN_PRCP => null()
  real, pointer, dimension(:,:)   ::        PHIS => null()
  real, pointer, dimension(:,:,:) ::     CNV_MFC => null()
  real, pointer, dimension(:,:,:) ::     CNV_MFD => null()
  real, pointer, dimension(:,:,:) ::      PFI_CN => null()
  real, pointer, dimension(:,:,:) ::      CNV_QC => null()

! Exports
  real,   pointer, dimension(:,:,:)  ::           delp => null()
  real,   pointer, dimension(:,:)    ::            LFR => null()
  real,   pointer, dimension(:,:,:)  ::          BYNCY => null()
  real,   pointer, dimension(:,:)    ::           CAPE => null()
  real*4, pointer, dimension(:,:,:)  ::  LIGHT_NO_PROD => null()

!=============================================================================
 
    integer                            :: k, k0

! Begin... 

!   Get the target components name and set-up traceback handle.
!   -----------------------------------------------------------
    Iam = "Run1"
    call ESMF_GridCompGet ( GC, name=COMP_NAME, __RC__ )
    Iam = trim(COMP_NAME) // "Run1"

    call MAPL_GetObjectFromGC ( GC, ggState, __RC__ )

    ! Get the time-step
    ! -----------------------
    call MAPL_GetResource( ggState, ndt, 'RUN_DT:', default=0, __RC__ )
    DT = ndt

    ! Get the grid related information
    !---------------------------------
    call ESMF_GridCompGet ( GC, GRID=esmfGrid, __RC__ )

    ! Get the local longitudes 
    !---------------------------------
    call ESMF_GridGetCoord(esmfGrid, localDE=0, coordDim=1, &
         staggerloc=ESMF_STAGGERLOC_CENTER, &
         datacopyFlag = ESMF_DATACOPY_REFERENCE,       &
         farrayPtr=R8D2, __RC__ )

    allocate(LONSLOCAL(size(R8D2,1),size(R8D2,2)), __STAT__ )

    LONSLOCAL = R8D2*(180/MAPL_PI)

    ! Get the local latitudes
    !---------------------------------
    call ESMF_GridGetCoord(esmfGrid, localDE=0, coordDim=2, &
         staggerloc=ESMF_STAGGERLOC_CENTER, &
         datacopyFlag = ESMF_DATACOPY_REFERENCE,       &
         farrayPtr=R8D2, __RC__ )

    allocate(LATSLOCAL(size(R8D2,1),size(R8D2,2)), __STAT__ )

    LATSLOCAL = R8D2*(180/MAPL_PI)

!   Every time, double check that the ratio has not changed
!   We pass in the year and month last used, to save reading again
!   --------------------------------------------------------------
    if ( flash_source_enum == FLASH_SOURCE_FIT ) then
      call update_lightning_ratio( ratioGlobalFile, CLOCK, &
                                   year_for_ratio, month_for_ratio, &
                                   ratioGlobalLight, __RC__ )
    end if

!   Get the imports...
!   ------------------
    call MAPL_GetPointer ( IMPORT,  PLE,  'PLE', __RC__ )

!   Get the exports...
!   ------------------
    call MAPL_GetPointer ( EXPORT, delp,   'DELP',   ALLOC=.TRUE.,     __RC__ )

!   Compute moist (rho) and dry (rhoDry) air density
!   ------------------------------------------------
    CALL AirDens ( GC, IMPORT, EXPORT, __RC__ )

    if ( associated(delp) ) then
       k0 = 1-lbound(PLE,3)

       do k = 1, size(delp,3)
          delp(:,:,k) = PLE(:,:,k-k0+1) - PLE(:,:,k-k0)
       end do
    end if


    !------------------------------------------------
    ! Flash Rate (LFR) for Lighting Parameterization
    !------------------------------------------------
    call MAPL_GetPointer ( EXPORT, LFR,            'LFR',             ALLOC=.TRUE., __RC__ )     !   Can we get rid of ALLOC=T ?
    call MAPL_GetPointer ( EXPORT, BYNCY,          'BYNCY',           ALLOC=.TRUE., __RC__ )
    call MAPL_GetPointer ( EXPORT, CAPE,           'CAPE',            ALLOC=.TRUE., __RC__ )
    call MAPL_GetPointer ( EXPORT, LIGHT_NO_PROD,  'LIGHT_NO_PROD',                 __RC__ )     ! If not requested, do not compute

    ! for FIT flashrate option
    if (flash_source_enum == FLASH_SOURCE_FIT) then
      call MAPL_GetPointer ( IMPORT, mcor,        'MCOR',         __RC__ )
      call MAPL_GetPointer ( IMPORT, RATIO_LOCAL, 'RATIO_LOCAL',  __RC__ )
      call MAPL_GetPointer ( IMPORT, MIDLAT_ADJ,  'MIDLAT_ADJ',   __RC__ )
    end if

    call MAPL_GetPointer ( IMPORT, CNV_MFC,     'CNV_MFC', __RC__ )
    call MAPL_GetPointer ( IMPORT, CNV_MFD,     'CNV_MFD', __RC__ )
    call MAPL_GetPointer ( IMPORT, CN_PRCP,     'CN_PRCP', __RC__ )
    call MAPL_GetPointer ( IMPORT, PHIS,        'PHIS',    __RC__ )
    call MAPL_GetPointer ( IMPORT, PFI_CN,      'PFI_CN', ALLOC=.TRUE., __RC__  )     !  ??

    call MAPL_GetPointer ( IMPORT, T,           'T',        __RC__ )
    call MAPL_GetPointer ( IMPORT, TH,          'TH',       __RC__ )
    call MAPL_GetPointer ( IMPORT, TS,          'TS',       __RC__ )
    call MAPL_GetPointer ( IMPORT, LWI,         'LWI',      __RC__ )
    call MAPL_GetPointer ( IMPORT, PBLH,        'ZPBL',     __RC__ )
    call MAPL_GetPointer ( IMPORT, FROCEAN,     'FROCEAN',  __RC__ )
    call MAPL_GetPointer ( IMPORT, FRLAND,      'FRLAND',   __RC__ )
    call MAPL_GetPointer ( IMPORT, Q,           'Q',        __RC__ )
    call MAPL_GetPointer ( IMPORT, ZLE,         'ZLE',      __RC__ )
    call MAPL_GetPointer ( IMPORT, CNV_QC,      'CNV_QC',   __RC__ )
    call MAPL_GetPointer ( IMPORT, cellArea,    'AREA',     __RC__ )

              BYNCY(:,:,:) = real(0)
               CAPE(:,:)   = real(0)
                LFR(:,:)   = real(0)

    if ( associated(LIGHT_NO_PROD) ) then
      LIGHT_NO_PROD(:,:,:) = REAL(0)
    end if

! print*,'ENV MOIST CN_PRCP ', MINVAL(CN_PRCP), MAXVAL(CN_PRCP)
! print*,'ENV MOIST FROCEAN ', MINVAL(FROCEAN), MAXVAL(FROCEAN)
! print*,'ENV MOIST TS      ', MINVAL(TS     ), MAXVAL(TS     )
! print*,'ENV MOIST CNV_MFC ', MINVAL(CNV_MFC), MAXVAL(CNV_MFC)
! print*,'ENV MOIST PLE     ', MINVAL(PLE    ), MAXVAL(PLE    )
! print*,'ENV MOIST TH      ', MINVAL(TH     ), MAXVAL(TH     )
! print*,'ENV MOIST Q       ', MINVAL(Q      ), MAXVAL(Q      )
! print*,'ENV MOIST ZLE     ', MINVAL(ZLE    ), MAXVAL(ZLE    )


!!  NOTE: CNV_MFD is dtrn in GMI

    call getLightning ( GC, ggState, CLOCK, &
                flash_source_enum, ratioGlobalLight, DT, &
                LONSLOCAL, LATSLOCAL, CN_PRCP, FRLAND, FROCEAN, LWI,   &
                PBLH, mcor, cellArea, MIDLAT_ADJ, RATIO_LOCAL,   &
                TS, CNV_MFC, CNV_QC, T, TH, PFI_CN, PLE, Q, ZLE,   &
                minDeepCloudTop, lightNOampFactor, numberNOperFlash, &
                MOIST_flashFactor, FIT_flashFactor, HEMCO_flashFactor, LOPEZ_flashFactor, &
                CNV_MFD, CAPE, BYNCY, LFR, LIGHT_NO_PROD, PHIS, &
                __RC__)

!   call pmaxmin('LFR', LFR, flashRateMin, flashRateMax, nc, 1, 1.)

    DEALLOCATE(LATSLOCAL, LONSLOCAL, __STAT__ )


!   All Done
!   --------
    RETURN_(ESMF_SUCCESS)

 end subroutine Run1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !IROUTINE: Run2 -- Phase 2 run method

! !INTERFACE:

  subroutine Run2 ( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:

  type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
  type(ESMF_State),    intent(inout) :: IMPORT ! Import state
  type(ESMF_State),    intent(inout) :: EXPORT ! Export state
  type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
  integer, optional,   intent(  out) :: RC     ! Error code

! !DESCRIPTION: The Run method of the Chemistry Composite Gridded Component.
!  It acts as a driver for the initializtion of the children.

!EOP

! ErrLog Variables

  character(len=ESMF_MAXSTR)           :: IAm 
  integer                              :: STATUS
  character(len=ESMF_MAXSTR)           :: COMP_NAME

! Imports
  real, pointer, dimension(:,:,:)      :: PLE => null()

! Exports
  real, pointer, dimension(:,:,:)      :: delp => null()
  real, pointer, dimension(:,:)        :: tprec => null()
  real, pointer, dimension(:,:)        :: cn_prcp => null()
  real, pointer, dimension(:,:)        :: ncn_prcp => null()

!=============================================================================
 
    type (MAPL_MetaComp), pointer      :: MAPL
    type (ESMF_FieldBundle)            :: Bundle
    type (ESMF_Grid)                   :: GRID
    type (ESMF_Time)                   :: CurrentTime
    character(len=ESMF_MAXSTR)         :: PRECIP_FILE
    integer                            :: year, month, day, hr, mn, se
  
    real, pointer, dimension(:,:)      :: pr_total
    real, pointer, dimension(:,:)      :: pr_conv
    real, pointer, dimension(:,:)      :: pr_snow

    logical                            :: observed_precip

    integer                            :: k, k0

    REAL, POINTER, DIMENSION(:,:,:)    :: DATA_FOR_OVP_3D => NULL()
    REAL, POINTER, DIMENSION(:,:,:)    :: OVP10_OUTPUT_3D => NULL()
    REAL, POINTER, DIMENSION(:,:,:)    :: OVP14_OUTPUT_3D => NULL()

    INTEGER                            :: CURRENT_HMS  !  for the end of the timestep

    real(ESMF_KIND_R4), pointer, dimension(:,:) :: LONS


! Begin... 

!   Get the target components name and set-up traceback handle.
!   -----------------------------------------------------------
    call ESMF_GridCompGet ( GC, name=COMP_NAME, Grid=GRID, __RC__ )
    Iam = trim(COMP_NAME) // "Run2"

! Get my internal MAPL_Generic state
!-----------------------------------
    call MAPL_GetObjectFromGC ( GC, MAPL, __RC__ )

!   Get to the imports...
!   ---------------------
    call MAPL_GetPointer ( IMPORT,  PLE,  'PLE', __RC__ )

!   Get to the exports...
!   ---------------------
    call MAPL_GetPointer ( EXPORT, delp,   'DELP',        __RC__ )

!   Compute moist (rho) and dry (rhoDry) air density
!   ------------------------------------------------
    CALL AirDens ( GC, IMPORT, EXPORT, __RC__ )

    if ( associated(delp) ) then
       k0 = 1-lbound(PLE,3)

       do k = 1, size(delp,3)
          delp(:,:,k) = PLE(:,:,k-k0+1) - PLE(:,:,k-k0)
       end do
    end if

! Import precip from MOIST or read observed precip from a file. 
! Export total, convective and non-convective precipitation.
!--------------------------------------------------------------
    call MAPL_GetPointer ( EXPORT, tprec,    'TPREC',     __RC__ )
    call MAPL_GetPointer ( EXPORT, cn_prcp,  'CN_PRCP',   __RC__ )
    call MAPL_GetPointer ( EXPORT, ncn_prcp, 'NCN_PRCP',  __RC__ )


    call MAPL_GetResource(MAPL, PRECIP_FILE, LABEL="PRECIP_FILE:", default="null", __RC__)

    if (PRECIP_FILE /= "null") then
        observed_precip = .true.        ! export observed precip
    else
        observed_precip = .false.       ! export modeled precip
    endif


    if (observed_precip) then
       bundle = ESMF_FieldBundleCreate (NAME='PRECIP', __RC__)
       call ESMF_FieldBundleSet(bundle, GRID=GRID, __RC__)

       call ESMF_ClockGet(CLOCK, currTime=CurrentTime, __RC__)

       call ESMF_TimeGet (CurrentTime, yy=year, mm=month, dd=day, h=hr, m=mn, s=se, __RC__)
       call ESMF_TimeSet (CurrentTime, yy=year, mm=month, dd=day, h=hr, m=30, s=0,  __RC__)

       call MAPL_CFIORead ( PRECIP_FILE, CurrentTime, Bundle, NOREAD=.true., __RC__ )
       call MAPL_CFIORead ( PRECIP_FILE, CurrentTime, Bundle, __RC__ )

       call ESMFL_BundleGetPointerToData ( Bundle, 'PRECTOT', pr_total, __RC__ )
       call ESMFL_BundleGetPointerToData ( Bundle, 'PRECCON', pr_conv,  __RC__ )
       call ESMFL_BundleGetPointerToData ( Bundle, 'PRECSNO', pr_snow,  __RC__ )
    else
       call MAPL_GetPointer ( IMPORT, pr_total, 'TPREC',   __RC__ )
       call MAPL_GetPointer ( IMPORT, pr_conv,  'CN_PRCP', __RC__ )
    end if


    if (associated(tprec))       tprec = pr_total
    if (associated(cn_prcp))   cn_prcp = pr_conv
    if (associated(ncn_prcp)) ncn_prcp = (pr_total - pr_conv)


    if (observed_precip) then
       call ESMF_FieldBundleDestroy(bundle, __RC__) 

       if (associated(pr_total)) deallocate(pr_total)
       if (associated(pr_conv))  deallocate(pr_conv)
       if (associated(pr_snow))  deallocate(pr_snow)
    end if

    IF ( OVP_setup_done .eqv. .FALSE. ) THEN

!     Set up Overpass Masks
!     --------------------
      CALL OVP_init ( GC, "CHEM_DT:", LONS, OVP_RUN_DT, OVP_GC_DT, __RC__ ) !  Get LONS, timesteps

!     IF(MAPL_AM_I_ROOT()) PRINT*,'in CHEMENV the RUN_DT and CHEM_DT values are: ', OVP_RUN_DT, OVP_GC_DT

      ! In this case we update the Exports only after each CHEMENV timestep:
      OVP_MASK_DT = OVP_GC_DT

      OVP_FIRST_HMS = OVP_end_of_timestep_hms( CLOCK, OVP_MASK_DT )
!     IF(MAPL_AM_I_ROOT()) PRINT*,'CHEMENV FIRST_HMS =',OVP_FIRST_HMS

      CALL OVP_mask ( LONS=LONS, DELTA_TIME=OVP_MASK_DT, OVERPASS_HOUR=10, MASK=MASK_10AM )
      CALL OVP_mask ( LONS=LONS, DELTA_TIME=OVP_MASK_DT, OVERPASS_HOUR=14, MASK=MASK_2PM  )

      OVP_setup_done = .TRUE.

    END IF


!  Record the Overpass values
!  -------------------------------------------------------------------

   CURRENT_HMS = OVP_end_of_timestep_hms( CLOCK, OVP_RUN_DT )
!  IF(MAPL_AM_I_ROOT()) PRINT*,'AGCM CURRENT_HMS =',CURRENT_HMS

! AIRDENS overpass

   CALL MAPL_GetPointer(EXPORT,  DATA_FOR_OVP_3D,         'AIRDENS', __RC__)
   CALL MAPL_GetPointer(EXPORT,  OVP10_OUTPUT_3D,   'OVP10_AIRDENS', __RC__)
   CALL MAPL_GetPointer(EXPORT,  OVP14_OUTPUT_3D,   'OVP14_AIRDENS', __RC__)

   CALL OVP_apply_mask( DATA_FOR_OVP_3D, OVP10_OUTPUT_3D, MASK_10AM, OVP_FIRST_HMS, CURRENT_HMS, K_EDGES=.FALSE., __RC__ )
   CALL OVP_apply_mask( DATA_FOR_OVP_3D, OVP14_OUTPUT_3D, MASK_2PM,  OVP_FIRST_HMS, CURRENT_HMS, K_EDGES=.FALSE., __RC__ )


!   All Done
!   --------
    RETURN_(ESMF_SUCCESS)

 end subroutine Run2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!-----------------------------------
! Compute Density for chemistry
!-----------------------------------

! !INTERFACE:

 subroutine Airdens ( GC, IMPORT, EXPORT, RC )

! !ARGUMENTS:

  type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
  type(ESMF_State),    intent(inout) :: IMPORT ! Import state
  type(ESMF_State),    intent(inout) :: EXPORT ! Export state
  integer, optional,   intent(  out) :: RC     ! Error code

! !DESCRIPTION: 
! Calculates air density used by chemistry and emissions.

!EOP

! Imports 
  real, pointer, dimension(:,:,:)      :: PLE => null()
  real, pointer, dimension(:,:,:)      :: th => null()
  real, pointer, dimension(:,:,:)      :: q  => null()

! Exports
  real, pointer, dimension(:,:,:)      :: rho => null()
  real, pointer, dimension(:,:,:)      :: rhoDry => null()

! Error handling
  character(len=ESMF_MAXSTR)           :: IAm = 'Airdens'
  integer                              :: STATUS
  character(len=ESMF_MAXSTR)           :: COMP_NAME

! Local variables 
  integer           :: k, k0, iml, jml, nl  ! dimensions
  real              :: eps
  real, allocatable :: npk(:,:,:)           ! normalized pk = (pe/p0)^kappa

!=============================================================================
! Begin... 

!   Get the target components name and set-up traceback handle.
!   -----------------------------------------------------------
    call ESMF_GridCompGet ( GC, name=COMP_NAME, __RC__ )
    Iam = trim(COMP_NAME) // "Airdens"

!   Get to the imports...
!   ---------------------
    call MAPL_GetPointer ( IMPORT,  PLE,  'PLE', __RC__ )
    call MAPL_GetPointer ( IMPORT,  th,  'TH',  __RC__ )
    call MAPL_GetPointer ( IMPORT,   q,  'Q',   __RC__ )

!   Get to the exports...
!   ---------------------
    call MAPL_GetPointer ( EXPORT, rho,    'AIRDENS',      __RC__ )
    call MAPL_GetPointer ( EXPORT, rhoDry, 'AIRDENS_DRYP', __RC__ )

!   Compute air density
!   -------------------
    iml = size(q,1)
    jml = size(q,2)
    nl  = size(q,3)
    k0  = 1-lbound(PLE,3)

    allocate(npk(iml,jml,nl+1),stat=STATUS) ! work space
    VERIFY_(STATUS)

    eps = MAPL_RVAP / MAPL_RGAS - 1.0

  !   Compute normalized PLE**Kappa
  !   ----------------------------
    npk = (PLE/MAPL_P00)**MAPL_KAPPA

  !   Compute rho from hydrostatic equation
  !   -------------------------------------
    DO k = 1, nl

        IF(ASSOCIATED(rho)) THEN

         rho(:,:,k) =       ( PLE(:,:,k-k0+1) - PLE(:,:,k-k0) ) /      &
                      ( MAPL_CP * ( th(:,:,k)*(1. + eps*q(:,:,k) ) ) &
                              * ( npk(:,:,k+1) - npk(:,:,k) ) )

        END IF

        IF(ASSOCIATED(rhoDry)) THEN

         rhoDry(:,:,k) =    ( 1.0 - q(:,:,k)                ) *      &
                            ( PLE(:,:,k-k0+1) - PLE(:,:,k-k0) ) /      &
                      ( MAPL_CP * ( th(:,:,k)*(1. + eps*q(:,:,k) ) ) &
                              * ( npk(:,:,k+1) - npk(:,:,k) ) )

        END IF

    END DO

    deallocate(npk)

!   All Done
!   --------
    RETURN_(ESMF_SUCCESS)

 end subroutine Airdens

end module GEOS_ChemEnvGridCompMod
