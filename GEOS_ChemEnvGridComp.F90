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

      !Derived types for the internal state
      type T_ChemEnv_STATE
         private
         INTEGER, ALLOCATABLE       :: MASK_10AM(:,:)
         INTEGER, ALLOCATABLE       :: MASK_2PM(:,:)
         INTEGER                    :: OVP_FIRST_HMS
         INTEGER                    :: OVP_RUN_DT
         INTEGER                    :: OVP_GC_DT
         INTEGER                    :: OVP_MASK_DT
         LOGICAL                    :: OVP_setup_done = .FALSE.
      end type T_ChemEnv_STATE

      type ChemEnv_WRAP
         type (T_ChemEnv_STATE), pointer :: PTR
      end type ChemEnv_WRAP

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
  logical                    :: usePreconCape       ! use CAPE, INHB and BYNCY from MOIST

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

    type (ESMF_Config)              :: CF
    type (T_ChemEnv_STATE), pointer :: state
    type (ChemEnv_WRAP)             :: wrap

!=============================================================================

! Begin...

! Get my name and set-up traceback handle
! ---------------------------------------

    Iam = 'ChemEnv::SetServices'
    call ESMF_GridCompGet( GC, NAME=COMP_NAME, CONFIG=CF, __RC__ )
    Iam = trim(COMP_NAME) // 'SetServices'

    ! Wrap internal state for storing in GC; rename legacyState
    ! -------------------------------------
    allocate ( state, stat=STATUS )
    VERIFY_(STATUS)
    wrap%ptr => state


! Register services for this component
! ------------------------------------
    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_INITIALIZE, Initialize_, __RC__ )
    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,        Run1,        __RC__ ) 
    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,        Run2,        __RC__ ) 

    ! Save pointer to the wrapped internal state in the GC
    !-----------------------------------------------------

    call ESMF_UserCompSetInternalState ( GC, 'ChemEnv', wrap, STATUS )
    VERIFY_(STATUS)
    state%OVP_setup_done = .FALSE.

!BOS

    call read_flash_source ( rcfilen, flash_source_enum, __RC__ )

! !IMPORT STATE:

    call MAPL_AddImportSpec(GC,                                         &
         SHORT_NAME = 'PLE',                                            &
         LONG_NAME  = 'air_pressure',                                   &
         UNITS      = 'Pa',                                             &
         DIMS       =  MAPL_DimsHorzVert,                               &
         VLOCATION  =  MAPL_VLocationEdge,                       __RC__ )

    call MAPL_AddImportSpec(GC,                                         &
         SHORT_NAME = 'TH',                                             &
         LONG_NAME  = 'potential_temperature',                          &
         UNITS      = 'K',                                              &
         DIMS       =  MAPL_DimsHorzVert,                               &
         VLOCATION  =  MAPL_VLocationCenter,                     __RC__ )

    call MAPL_AddImportSpec(GC,                                         &
         SHORT_NAME = 'Q',                                              &
         LONG_NAME  = 'specific_humidity',                              &
         UNITS      = 'kg kg-1',                                        &
         DIMS       = MAPL_DimsHorzVert,                                &
         VLOCATION  = MAPL_VLocationCenter,                      __RC__ )

!   Geopotential Height
!   -------------------
    call MAPL_AddImportSpec(GC,                                         &
         SHORT_NAME = 'PHIS',                                           &
         LONG_NAME  = 'Geopotential Height at Surface',                 &
         UNITS      = 'm+2 s-2',                                        &
         DIMS       =  MAPL_DimsHorzOnly,                               &
         VLOCATION  =  MAPL_VLocationNone,                       __RC__ )

!   Total precip
!   ------------
    call MAPL_AddImportSpec(GC,                                         &
         SHORT_NAME = 'PRECTOT',                                        &
         LONG_NAME  = 'total_precipitation',                            &
         UNITS      = 'kg m-2 s-1',                                     &
         DEFAULT    =  MAPL_UNDEF,                                      &
         DIMS       =  MAPL_DimsHorzOnly,                               &
         VLOCATION  =  MAPL_VLocationNone,                       __RC__ )

!   Convective precip
!   -----------------
    call MAPL_AddImportSpec(GC,                                         &
         SHORT_NAME = 'CN_PRCP',                                        &
         LONG_NAME  = 'convective_precipitation',                       &
         UNITS      = 'kg m-2 s-1',                                     &
         DEFAULT    =  MAPL_UNDEF,                                      &
         DIMS       =  MAPL_DimsHorzOnly,                               &
         VLOCATION  =  MAPL_VLocationNone,                       __RC__ )

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
      call MAPL_AddImportSpec(GC,                                       &
           SHORT_NAME = 'MCOR',                                         &
           LONG_NAME  = 'agrid_cell_area',                              &
           UNITS      = 'm+2',                                          &
           DIMS       = MAPL_DimsHorzOnly,                              &
           VLOCATION  = MAPL_VLocationNone,                      __RC__ )

      call MAPL_AddImportSpec(GC,                                       &
           SHORT_NAME = 'RATIO_LOCAL',                                  &
           LONG_NAME  = 'local_ratios_lightning',                       &
           UNITS      = '1',                                            &
           DIMS       = MAPL_DimsHorzOnly,                              &
           VLOCATION  = MAPL_VLocationNone,                      __RC__ )
    
      call MAPL_AddImportSpec(GC,                                       &
           SHORT_NAME = 'MIDLAT_ADJ',                                   &
           LONG_NAME  = 'midlat_adjustment_lightning',                  &
           UNITS      = '1',                                            &
           DIMS       = MAPL_DimsHorzOnly,                              &
           VLOCATION  = MAPL_VLocationNone,                      __RC__ )
    endif

    call MAPL_AddImportSpec(GC,                                         &
         SHORT_NAME = 'CNV_MFC',                                        &
         LONG_NAME  = 'cumulative_mass_flux',                           &
         UNITS      = 'kg m-2 s-1',                                     &
         DIMS       = MAPL_DimsHorzVert,                                &
         VLOCATION  = MAPL_VLocationEdge,                        __RC__ )

    call MAPL_AddImportSpec(GC,                                         &
         SHORT_NAME = 'CNV_MFD',                                        &
         LONG_NAME  = 'detraining_mass_flux',                           &
         UNITS      = 'kg m-2 s-1',                                     &
         DIMS       = MAPL_DimsHorzVert,                                &
         VLOCATION  = MAPL_VLocationCenter,                      __RC__ )

    call MAPL_AddImportSpec(GC,                                         &
         SHORT_NAME = 'PFI_CN',                                         &
         LONG_NAME  = '3D_flux_of_ice_convective_precipitation',                   &
         UNITS      = 'kg m-2 s-1',                                     &
         DIMS       = MAPL_DimsHorzVert,                                &
         VLOCATION  = MAPL_VLocationEdge,                        __RC__ )

    call MAPL_AddImportSpec ( gc,                                       &
         SHORT_NAME = 'T',                                              &
         LONG_NAME  = 'air_temperature',                                &
         UNITS      = 'K',                                              &
         DIMS       = MAPL_DimsHorzVert,                                &
         VLOCATION  = MAPL_VLocationCenter,                      __RC__ )

    call MAPL_AddImportSpec(GC,                                         &
         SHORT_NAME = 'TS',                                             &
         LONG_NAME  = 'surface_skin_temperature',                       &
         UNITS      = 'K',                                              &
         DIMS       = MAPL_DimsHorzOnly,                                &
         VLOCATION  = MAPL_VLocationNone,                        __RC__ )

    call MAPL_AddImportSpec(GC,                                         &
         SHORT_NAME = 'FROCEAN',                                        &
         LONG_NAME  = 'fraction_of_ocean',                              &
         UNITS      = '1',                                              &
         DIMS       = MAPL_DimsHorzOnly,                                &
         VLOCATION  = MAPL_VLocationNone,                        __RC__ )

    call MAPL_AddImportSpec(GC,                                         &
         SHORT_NAME = 'FRLAND',                                         &
         LONG_NAME  = 'fraction_of_land',                               &
         UNITS      = '1',                                              &
         DIMS       = MAPL_DimsHorzOnly,                                &
         VLOCATION  = MAPL_VLocationNone,                        __RC__ )

    call MAPL_AddImportSpec(GC,                                         &
         SHORT_NAME = 'ZLE',                                            &
         LONG_NAME  = 'geopotential_height',                            &
         UNITS      = 'm',                                              &
         DIMS       = MAPL_DimsHorzVert,                                &
         VLOCATION  = MAPL_VLocationEdge,                        __RC__ )

    call MAPL_AddImportSpec(GC,                                         &
         SHORT_NAME = 'LWI',                                            &
         LONG_NAME  = 'land(1)_water(0)_ice(2)_flag',                   &
         UNITS      = '1',                                              &
         DIMS       = MAPL_DimsHorzOnly,                                &
         VLOCATION  = MAPL_VLocationNone,                        __RC__ )

    call MAPL_AddImportSpec(GC,                                         &
         SHORT_NAME = 'ZPBL',                                           &
         LONG_NAME  = 'planetary_boundary_layer_height',                &
         UNITS      = 'm',                                              &
         DIMS       = MAPL_DimsHorzOnly,                                &
         VLOCATION  = MAPL_VLocationNone,                        __RC__ )

    call MAPL_AddImportSpec(GC,                                         &
         SHORT_NAME = 'AREA',                                           &
         LONG_NAME  = 'agrid_cell_area',                                &
         UNITS      = 'm+2',                                            &
         DIMS       = MAPL_DimsHorzOnly,                                &
         VLOCATION  = MAPL_VLocationNone,                        __RC__ )

    call MAPL_AddImportSpec(GC,                                         &
         SHORT_NAME = 'CNV_QC',                                         &
         LONG_NAME  = 'grid_mean_convective_condensate',                &
         UNITS      = 'kg kg-1',                                        &
         DIMS       = MAPL_DimsHorzVert,                                &
         VLOCATION  = MAPL_VLocationCenter,                      __RC__ )

    call MAPL_AddImportSpec ( gc,                                       &
         SHORT_NAME = 'QLTOT',                                          &
         LONG_NAME  = 'mass_fraction_of_cloud_liquid_water',            &
         UNITS      = 'kg kg-1',                                        &
         DIMS       = MAPL_DimsHorzVert,                                &
         VLOCATION  = MAPL_VLocationCenter,                      __RC__ )

    call MAPL_AddImportSpec(GC,                                         &
         SHORT_NAME = 'PS',                                             &
         LONG_NAME  = 'surface_pressure',                               &
         UNITS      = 'Pa',                                             &
         DIMS       = MAPL_DimsHorzOnly,                                &
         VLOCATION  = MAPL_VLocationNone,                         __RC__)

    call MAPL_AddImportSpec(GC,                                         &
         SHORT_NAME = 'U10M',                                           &
         LONG_NAME  = '10-meter_eastward_wind',                         &
         UNITS      = 'm s-1',                                          &
         DIMS       = MAPL_DimsHorzOnly,                                &
         VLOCATION  = MAPL_VLocationNone,                         __RC__)

    call MAPL_AddImportSpec(GC,                                         &
         SHORT_NAME = 'V10M',                                           &
         LONG_NAME  = '10-meter_northward_wind',                        &
         UNITS      = 'm s-1',                                          &
         DIMS       = MAPL_DimsHorzOnly,                                &
         VLOCATION  = MAPL_VLocationNone,                         __RC__)

!    call MAPL_AddImportSpec ( gc,                                       &
!         SHORT_NAME = 'EPV',                                            &
!         LONG_NAME  = 'ertels_potential_vorticity',                     &
!         UNITS      = 'K m+2 kg-1 s-1',                                 &
!         RESTART    = MAPL_RestartSkip,                                 &
!         DIMS       = MAPL_DimsHorzVert,                                &
!         VLOCATION  = MAPL_VLocationCenter,                       __RC__)

    call MAPL_AddImportSpec(GC,                                         &
         SHORT_NAME = 'PPBL',                                           &
         LONG_NAME  = 'pbltop_pressure',                                &
         UNITS      = '1',                                              &
         RESTART    = MAPL_RestartSkip,                                 &
         DIMS       = MAPL_DimsHorzOnly,                                &
         VLOCATION  = MAPL_VLocationNone,                         __RC__)

    call MAPL_AddImportSpec(GC,                                         &
         SHORT_NAME = 'TROPP',                                          &
         LONG_NAME  = 'tropopause_pressure_based_on_blended_estimate',  &
         UNITS      = 'Pa',                                             &
         DIMS       = MAPL_DimsHorzOnly,                                &
         VLOCATION  = MAPL_VLocationNone,                         __RC__)

    call MAPL_AddImportSpec(GC,                                         &
         SHORT_NAME ='BYNCY',                                           &
         LONG_NAME  ='buoyancy_of surface_parcel',                      &
         UNITS      ='m s-2',                                           &
         RESTART    = MAPL_RestartSkip,                                 &
         DIMS       = MAPL_DimsHorzVert,                                &
         VLOCATION  = MAPL_VLocationCenter,                       __RC__)

    call MAPL_AddImportSpec(GC,                                         &
         SHORT_NAME ='CAPE',                                            &
         LONG_NAME  ='cape_for_surface_parcel',                         &
         UNITS      ='J kg-1',                                          &
         RESTART    = MAPL_RestartSkip,                                 &
         DIMS       = MAPL_DimsHorzOnly,                                &
         VLOCATION  = MAPL_VLocationNone,                         __RC__)
       
    call MAPL_AddImportSpec(GC,                                         &
         SHORT_NAME ='INHB',                                            &
         LONG_NAME  ='inhibition_for_surface_parcel',                   &
         UNITS      ='J kg-1',                                          &
         RESTART    = MAPL_RestartSkip,                                 &
         DIMS       = MAPL_DimsHorzOnly,                                & 
         VLOCATION  = MAPL_VLocationNone,                         __RC__)


    call MAPL_AddImportSpec(GC,                               &
         SHORT_NAME='ZLCL',                                        &
         LONG_NAME ='lifting_condensation_level',                  &
         UNITS     ='m'  ,                                         &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)
  
    call MAPL_AddImportSpec(GC,                               &
         SHORT_NAME='ZLFC',                                        &
         LONG_NAME ='level_of_free_convection',                    &
         UNITS     ='m'  ,                                         &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)



! !EXPORT STATE:

!    AIRDENS: Provided for Children
!    ------------------------------
    call MAPL_AddExportSpec(GC,                                         &
         SHORT_NAME = 'AIRDENS',                                        &
         LONG_NAME  = 'moist_air_density',                              &
         UNITS      = 'kg m-3',                                         &
         DIMS       = MAPL_DimsHorzVert,                                &
         VLOCATION  = MAPL_VLocationCenter,                      __RC__ )

!    Density of dry air
!    ------------------
    call MAPL_AddExportSpec(GC,                                         &
         SHORT_NAME = 'AIRDENS_DRYP',                                   &
         LONG_NAME  = 'partial_dry_air_density',                        &
         UNITS      = 'kg dry m-3 tot',                                 &
         DIMS       = MAPL_DimsHorzVert,                                &
         VLOCATION  = MAPL_VLocationCenter,                      __RC__ )

!    DELP (This should be wired from DYN)
!    ------------------------------------
    call MAPL_AddExportSpec(GC,                                         &
         SHORT_NAME = 'DELP',                                           &
         LONG_NAME  = 'pressure_thickness',                             &
         UNITS      = 'Pa',                                             &
         DIMS       = MAPL_DimsHorzVert,                                &
         VLOCATION  = MAPL_VLocationCenter,                      __RC__ )

!   Total precip
!   ------------
    call MAPL_AddExportSpec(GC,                                         &
         SHORT_NAME = 'TPREC',                                          &
         LONG_NAME  = 'total_precipitation',                            &
         UNITS      = 'kg m-2 s-1',                                     &
         DIMS       = MAPL_DimsHorzOnly,                                &
         VLOCATION  = MAPL_VLocationNone,                        __RC__ )

!    Convective precip
!    -----------------
    call MAPL_AddExportSpec(GC,                                         &
         SHORT_NAME = 'CN_PRCP',                                        &
         LONG_NAME  = 'Convective precipitation',                       &
         UNITS      = 'kg m-2 s-1',                                     &
         DIMS       = MAPL_DimsHorzOnly,                                &
         VLOCATION  = MAPL_VLocationNone,                        __RC__ )

!    Non-convective precip
!    ---------------------
    call MAPL_AddExportSpec(GC,                                         &
         SHORT_NAME = 'NCN_PRCP',                                       &
         LONG_NAME  = 'Non-convective precipitation',                   &
         UNITS      = 'kg m-2 s-1',                                     &
         DIMS       = MAPL_DimsHorzOnly,                                &
         VLOCATION  = MAPL_VLocationNone,                        __RC__ )


    call MAPL_AddExportSpec(GC,                                         &
         SHORT_NAME = 'LFR',                                            &
         LONG_NAME  = 'lightning_flash_rate',                           &
         UNITS      = 'km-2 s-1',                                       &
         DIMS       = MAPL_DimsHorzOnly,                                &
         VLOCATION  = MAPL_VLocationNone,                        __RC__ )

    call MAPL_AddExportSpec(GC,                                         &
         SHORT_NAME = 'BYNCY',                                          &
         LONG_NAME  = 'buoyancy_of surface_parcel',                     &
         UNITS      = 'm s-2',                                          &
         DIMS       =  MAPL_DimsHorzVert,                               &
         VLOCATION  =  MAPL_VLocationCenter,                     __RC__ )

    call MAPL_AddExportSpec(GC,                                         &
         SHORT_NAME = 'CAPE',                                           &
         LONG_NAME  = 'convective_avail_pot_energy',                    &
         UNITS      = 'J m^{-2} <check this!>',                         &
         DIMS       = MAPL_DimsHorzOnly,                                &
         VLOCATION  = MAPL_VLocationNone,                        __RC__ )

    call MAPL_AddExportSpec(GC,                                         &
         SHORT_NAME = 'LIGHT_NO_PROD',                                  &
         LONG_NAME  = 'lightning_NO_prod_rate',                         &
         UNITS      = 'm-3 s-1',                                        &
         DIMS       =  MAPL_DimsHorzVert,                               &
         PRECISION  =  ESMF_KIND_R4,                                    &
         VLOCATION  =  MAPL_VLocationCenter,                     __RC__ )

!!!! >>>>>>>>>>>>>>>>  OVP

!    10am overpass AIRDENS
!    ---------------------
    call MAPL_AddExportSpec(GC,                                         &
         SHORT_NAME = 'OVP10_AIRDENS',                                  &
         LONG_NAME  = 'moist_air_density_10am_local',                   &
         UNITS      = 'kg m-3',                                         &
         DIMS       = MAPL_DimsHorzVert,                                &
         VLOCATION  = MAPL_VLocationCenter,                      __RC__ )

!    2pm  overpass AIRDENS
!    ---------------------
    call MAPL_AddExportSpec(GC,                                         &
         SHORT_NAME = 'OVP14_AIRDENS',                                  &
         LONG_NAME  = 'moist_air_density_2pm_local',                    &
         UNITS      = 'kg m-3',                                         &
         DIMS       = MAPL_DimsHorzVert,                                &
         VLOCATION  = MAPL_VLocationCenter,                      __RC__ )

    CALL MAPL_AddExportSpec(GC,                                         &
         SHORT_NAME = 'OVP10_T',                                        &
         LONG_NAME  = 'air_temperature_10am_local',                     &
         UNITS      = 'K',                                              &
         DIMS       = MAPL_DimsHorzVert,                                &
         VLOCATION  = MAPL_VLocationCenter,                      __RC__ )

    CALL MAPL_AddExportSpec(GC,                                         &
         SHORT_NAME = 'OVP14_T',                                        &
         LONG_NAME  = 'air_temperature_2pm_local',                      &
         UNITS      = 'K',                                              &
         DIMS       = MAPL_DimsHorzVert,                                &
         VLOCATION  = MAPL_VLocationCenter,                      __RC__ )

    CALL MAPL_AddExportSpec(GC,                                         &
         SHORT_NAME = 'OVP10_PL',                                       &
         LONG_NAME  = 'mid_level_pressure_10am_local',                  &
         UNITS      = 'Pa',                                             &
         DIMS       = MAPL_DimsHorzVert,                                &
         VLOCATION  = MAPL_VLocationCenter,                      __RC__ )

    CALL MAPL_AddExportSpec(GC,                                         &
         SHORT_NAME = 'OVP14_PL',                                       &
         LONG_NAME  = 'mid_level_pressure_2pm_local',                   &
         UNITS      = 'Pa',                                             &
         DIMS       = MAPL_DimsHorzVert,                                &
         VLOCATION  = MAPL_VLocationCenter,                      __RC__ )

    CALL MAPL_AddExportSpec(GC,                                         &
         SHORT_NAME = 'OVP10_PLE',                                      &
         LONG_NAME  = 'edge_pressure_10am_local',                       &
         UNITS      = 'Pa',                                             &
         DIMS       = MAPL_DimsHorzVert,                                &
         VLOCATION  = MAPL_VLocationEdge,                        __RC__ )

    CALL MAPL_AddExportSpec(GC,                                         &
         SHORT_NAME = 'OVP14_PLE',                                      &
         LONG_NAME  = 'edge_pressure_2pm_local',                        &
         UNITS      = 'Pa',                                             &
         DIMS       = MAPL_DimsHorzVert,                                &
         VLOCATION  = MAPL_VLocationEdge,                        __RC__ )

    CALL MAPL_AddExportSpec(GC,                                         &
         SHORT_NAME = 'OVP10_QV_VMR',                                   &
         LONG_NAME  = 'water_vapor_10am_local',                         &
         UNITS      = 'mol mol-1',                                      &
         DIMS       = MAPL_DimsHorzVert,                                &
         VLOCATION  = MAPL_VLocationCenter,                      __RC__ )

    CALL MAPL_AddExportSpec(GC,                                         &
         SHORT_NAME = 'OVP14_QV_VMR',                                   &
         LONG_NAME  = 'water_vapor_2pm_local',                          &
         UNITS      = 'mol mol-1',                                      &
         DIMS       = MAPL_DimsHorzVert,                                &
         VLOCATION  = MAPL_VLocationCenter,                      __RC__ )

    CALL MAPL_AddExportSpec(GC,                                         &
         SHORT_NAME = 'OVP10_QLTOT',                                    &
         LONG_NAME  = 'mass_fraction_of_cloud_liquid_water_10am_local', &
         UNITS      = 'kg kg-1',                                        &
         DIMS       = MAPL_DimsHorzVert,                                &
         VLOCATION  = MAPL_VLocationCenter,                      __RC__ )

    CALL MAPL_AddExportSpec(GC,                                         &
         SHORT_NAME = 'OVP14_QLTOT',                                    &
         LONG_NAME  = 'mass_fraction_of_cloud_liquid_water_2pm_local',  &
         UNITS      = 'kg kg-1',                                        &
         DIMS       = MAPL_DimsHorzVert,                                &
         VLOCATION  = MAPL_VLocationCenter,                      __RC__ )

    CALL MAPL_AddExportSpec(GC,                                         &
         SHORT_NAME = 'OVP10_PS',                                       &
         LONG_NAME  = 'surface_pressure_10am_local',                    &
         UNITS      = 'Pa',                                             &
         DIMS       = MAPL_DimsHorzOnly,                                &
         VLOCATION  = MAPL_VLocationNone,                        __RC__ )

    CALL MAPL_AddExportSpec(GC,                                         &
         SHORT_NAME = 'OVP14_PS',                                       &
         LONG_NAME  = 'surface_pressure_2pm_local',                     &
         UNITS      = 'Pa',                                             &
         DIMS       = MAPL_DimsHorzOnly,                                &
         VLOCATION  = MAPL_VLocationNone,                        __RC__ )

    CALL MAPL_AddExportSpec(GC,                                         &
         SHORT_NAME = 'OVP10_PPBL',                                     &
         LONG_NAME  = 'pbltop_pressure_10am_local',                     &
         UNITS      = 'Pa',                                             &
         DIMS       = MAPL_DimsHorzOnly,                                &
         VLOCATION  = MAPL_VLocationNone,                        __RC__ )

    CALL MAPL_AddExportSpec(GC,                                         &
         SHORT_NAME = 'OVP14_PPBL',                                     &
         LONG_NAME  = 'pbltop_pressure_2pm_local',                      &
         UNITS      = 'Pa',                                             &
         DIMS       = MAPL_DimsHorzOnly,                                &
         VLOCATION  = MAPL_VLocationNone,                        __RC__ )

    CALL MAPL_AddExportSpec(GC,                                         &
         SHORT_NAME = 'OVP10_TROPP',                                    &
         LONG_NAME  = 'tropopause_pressure_10am_local',                 &
         UNITS      = 'Pa',                                             &
         DIMS       = MAPL_DimsHorzOnly,                                &
         VLOCATION  = MAPL_VLocationNone,                        __RC__ )

    CALL MAPL_AddExportSpec(GC,                                         &
         SHORT_NAME = 'OVP14_TROPP',                                    &
         LONG_NAME  = 'tropopause_pressure_2pm_local',                  &
         UNITS      = 'Pa',                                             &
         DIMS       = MAPL_DimsHorzOnly,                                &
         VLOCATION  = MAPL_VLocationNone,                        __RC__ )

    CALL MAPL_AddExportSpec(GC,                                         &
         SHORT_NAME = 'OVP10_U10M',                                     &
         LONG_NAME  = 'eastward_10m_wind_speed_10am_local',             &
         UNITS      = 'm s-1',                                          &
         DIMS       = MAPL_DimsHorzOnly,                                &
         VLOCATION  = MAPL_VLocationNone,                        __RC__ )

    CALL MAPL_AddExportSpec(GC,                                         &
         SHORT_NAME = 'OVP14_U10M',                                     &
         LONG_NAME  = 'eastward_10m_wind_speed_2pm_local',              &
         UNITS      = 'm s-1',                                          &
         DIMS       = MAPL_DimsHorzOnly,                                &
         VLOCATION  = MAPL_VLocationNone,                        __RC__ )

    CALL MAPL_AddExportSpec(GC,                                         &
         SHORT_NAME = 'OVP10_V10M',                                     &
         LONG_NAME  = 'northward_10m_wind_speed_10am_local',            &
         UNITS      = 'm s-1',                                          &
         DIMS       = MAPL_DimsHorzOnly,                                &
         VLOCATION  = MAPL_VLocationNone,                        __RC__ )

    CALL MAPL_AddExportSpec(GC,                                         &
         SHORT_NAME = 'OVP14_V10M',                                     &
         LONG_NAME  = 'northward_10m_wind_speed_2pm_local',             &
         UNITS      = 'm s-1',                                          &
         DIMS       = MAPL_DimsHorzOnly,                                &
         VLOCATION  = MAPL_VLocationNone,                        __RC__ )

!    CALL MAPL_AddExportSpec(GC,                                         &
!         SHORT_NAME = 'OVP10_EPV',                                      &
!         LONG_NAME  = 'ertels_potential_vorticity_10am_local',          &
!         UNITS      = 'K m+2 kg-1 s-1',                                 &
!         DIMS       = MAPL_DimsHorzVert,                                &
!         VLOCATION  = MAPL_VLocationCenter,                      __RC__ )

!    CALL MAPL_AddExportSpec(GC,                                         &
!         SHORT_NAME = 'OVP14_EPV',                                      &
!         LONG_NAME  = 'ertels_potential_vorticity_2pm_local',           &
!         UNITS      = 'K m+2 kg-1 s-1',                                 &
!         DIMS       = MAPL_DimsHorzVert,                                &
!         VLOCATION  = MAPL_VLocationCenter,                      __RC__ )

!!!! <<<<<<<<<<<<<<<<  OVP


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
!   IF(MAPL_AM_I_ROOT()) PRINT*,'The SimType resource is '//TRIM(rc_string)

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
                                 usePreconCape,                        &
                                 __RC__ )

    IF(MAPL_AM_I_ROOT()) THEN
      if ( flash_source_enum == FLASH_SOURCE_MOIST ) PRINT*,'MOIST_flashFactor is ',MOIST_flashFactor
      if ( flash_source_enum == FLASH_SOURCE_FIT   ) PRINT*,'  FIT_flashFactor is ',  FIT_flashFactor
      if ( flash_source_enum == FLASH_SOURCE_HEMCO ) PRINT*,'HEMCO_flashFactor is ',HEMCO_flashFactor
      if ( flash_source_enum == FLASH_SOURCE_LOPEZ ) PRINT*,'LOPEZ_flashFactor is ',LOPEZ_flashFactor

                                                     PRINT*,'usePreconCape = ', usePreconCape
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

  real, pointer, dimension(:,:)   ::          ZLFC => null()
  real, pointer, dimension(:,:)   ::          ZLCL => null()

  real, pointer, dimension(:,:)   ::          TS => null()
  real, pointer, dimension(:,:)   ::     FROCEAN => null()
  real, pointer, dimension(:,:)   ::      FRLAND => null()
  real, pointer, dimension(:,:)   ::     CN_PRCP => null()
  real, pointer, dimension(:,:)   ::        PHIS => null()
  real, pointer, dimension(:,:,:) ::     CNV_MFC => null()
  real, pointer, dimension(:,:,:) ::     CNV_MFD => null()
  real, pointer, dimension(:,:,:) ::      PFI_CN => null()
  real, pointer, dimension(:,:,:) ::      CNV_QC => null()

  real, pointer, dimension(:,:)   ::       CAPE_PRECON => null()
  real, pointer, dimension(:,:)   ::       INHB_PRECON => null()
  real, pointer, dimension(:,:,:) ::      BYNCY_PRECON => null()

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
    call MAPL_GetPointer ( IMPORT, PFI_CN,      'PFI_CN',   __RC__ )

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

    call MAPL_GetPointer ( IMPORT,  CAPE_PRECON, 'CAPE',    __RC__ )
    call MAPL_GetPointer ( IMPORT,  INHB_PRECON, 'INHB',    __RC__ )
    call MAPL_GetPointer ( IMPORT, BYNCY_PRECON, 'BYNCY',   __RC__ )

    call MAPL_GetPointer ( IMPORT, ZLFC, 'ZLFC',   __RC__ )
    call MAPL_GetPointer ( IMPORT, ZLCL, 'ZLCL',   __RC__ )


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
                CNV_MFD, usePreconCape, CAPE_PRECON, INHB_PRECON, BYNCY_PRECON, &
                CAPE, BYNCY, ZLFC, ZLCL, LFR, LIGHT_NO_PROD, PHIS, &
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

  type (T_ChemEnv_STATE), pointer  :: ChemEnv_STATE
  type (ChemEnv_wrap)              :: WRAP

! Imports
  real, pointer, dimension(:,:,:)      :: PLE => null()

! Exports
  real, pointer, dimension(:,:,:)      :: delp => null()
  real, pointer, dimension(:,:)        ::    tprec => null()
  real, pointer, dimension(:,:)        ::  cn_prcp => null()
  real, pointer, dimension(:,:)        :: ncn_prcp => null()

!=============================================================================
 
    type (MAPL_MetaComp), pointer      :: MAPL
    type (ESMF_Grid)                   :: GRID
    type (ESMF_Time)                   :: CurrentTime
    integer                            :: year, month, day, hr, mn, se
  
    real, pointer, dimension(:,:)      :: pr_total
    real, pointer, dimension(:,:)      :: pr_conv

    integer                            :: k, k0

    REAL, POINTER, DIMENSION(:,:,:)    :: DATA_FOR_OVP_3D => NULL()
    REAL, POINTER, DIMENSION(:,:,:)    :: OVP10_OUTPUT_3D => NULL()
    REAL, POINTER, DIMENSION(:,:,:)    :: OVP14_OUTPUT_3D => NULL()
    REAL, POINTER, DIMENSION(:,:)      :: DATA_FOR_OVP_2D => NULL()
    REAL, POINTER, DIMENSION(:,:)      :: OVP10_OUTPUT_2D => NULL()
    REAL, POINTER, DIMENSION(:,:)      :: OVP14_OUTPUT_2D => NULL()

    INTEGER                            :: CURRENT_HMS  !  for the end of the timestep

    real(ESMF_KIND_R4), pointer, dimension(:,:) :: LONS
    REAL, ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: scratch_3d
    INTEGER                                     :: idim, jdim, kdim


! Begin... 

!   Get the target components name and set-up traceback handle.
!   -----------------------------------------------------------
    call ESMF_GridCompGet ( GC, name=COMP_NAME, Grid=GRID, __RC__ )
    Iam = trim(COMP_NAME) // "Run2"

! Get my internal MAPL_Generic state
!-----------------------------------
    call MAPL_GetObjectFromGC ( GC, MAPL, __RC__ )

    ! Get my private state from the component
    !----------------------------------------
    call ESMF_UserCompGetInternalState(gc, 'ChemEnv', WRAP, STATUS)
    VERIFY_(STATUS)

    ChemEnv_STATE => WRAP%PTR

!   Get to the imports...
!   ---------------------
    call MAPL_GetPointer ( IMPORT,  PLE,  'PLE', __RC__ )

!   Dimensions of fields with VLOC=center
!   -------------------------------------
    idim = size(PLE,1)
    jdim = size(PLE,2)
    kdim = size(PLE,3) - 1

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

! Import total precip from SURFACE (either model or observed precip)
! Export total, convective and non-convective precipitation.
!--------------------------------------------------------------
    call MAPL_GetPointer ( EXPORT,    tprec,    'TPREC',  __RC__ )
    call MAPL_GetPointer ( EXPORT,  cn_prcp,  'CN_PRCP',  __RC__ )
    call MAPL_GetPointer ( EXPORT, ncn_prcp, 'NCN_PRCP',  __RC__ )


    call MAPL_GetPointer ( IMPORT, pr_total, 'PRECTOT', __RC__ )
    call MAPL_GetPointer ( IMPORT, pr_conv,  'CN_PRCP', __RC__ )
 

    if (associated(   tprec))    tprec =  pr_total
    if (associated( cn_prcp))  cn_prcp =  pr_conv
    if (associated(ncn_prcp)) ncn_prcp = (pr_total - pr_conv)


    IF ( ChemEnv_STATE%OVP_setup_done .eqv. .FALSE. ) THEN

!     Set up Overpass Masks
!     --------------------
      CALL OVP_init ( GC, "CHEM_DT:", LONS, ChemEnv_STATE%OVP_RUN_DT, ChemEnv_STATE%OVP_GC_DT, __RC__ ) !  Get LONS, timesteps

!     IF(MAPL_AM_I_ROOT()) PRINT*,'in CHEMENV the RUN_DT and CHEM_DT values are: ', ChemEnv_STATE%OVP_RUN_DT, ChemEnv_STATE%OVP_GC_DT

      ! In this case we update the Exports only after each CHEMENV timestep:
      ChemEnv_STATE%OVP_MASK_DT = ChemEnv_STATE%OVP_GC_DT

      ChemEnv_STATE%OVP_FIRST_HMS = OVP_end_of_timestep_hms( CLOCK, ChemEnv_STATE%OVP_MASK_DT )
!     IF(MAPL_AM_I_ROOT()) PRINT*,'CHEMENV FIRST_HMS =',ChemEnv_STATE%OVP_FIRST_HMS

      CALL OVP_mask ( LONS=LONS, DELTA_TIME=ChemEnv_STATE%OVP_MASK_DT, &
                      OVERPASS_HOUR=10, MASK=ChemEnv_STATE%MASK_10AM )
      CALL OVP_mask ( LONS=LONS, DELTA_TIME=ChemEnv_STATE%OVP_MASK_DT, &
                      OVERPASS_HOUR=14, MASK=ChemEnv_STATE%MASK_2PM  )

      ChemEnv_STATE%OVP_setup_done = .TRUE.

    END IF


!  Record the Overpass values
!  -------------------------------------------------------------------

   CURRENT_HMS = OVP_end_of_timestep_hms( CLOCK, ChemEnv_STATE%OVP_RUN_DT )
!  IF(MAPL_AM_I_ROOT()) PRINT*,'AGCM CURRENT_HMS =',CURRENT_HMS
!
! AIRDENS overpass
!
   CALL MAPL_GetPointer(EXPORT,  DATA_FOR_OVP_3D,         'AIRDENS', __RC__)
   CALL MAPL_GetPointer(EXPORT,  OVP10_OUTPUT_3D,   'OVP10_AIRDENS', __RC__)
   CALL MAPL_GetPointer(EXPORT,  OVP14_OUTPUT_3D,   'OVP14_AIRDENS', __RC__)

   CALL OVP_apply_mask( DATA_FOR_OVP_3D, OVP10_OUTPUT_3D, ChemEnv_STATE%MASK_10AM, &
                        ChemEnv_STATE%OVP_FIRST_HMS, CURRENT_HMS, K_EDGES=.FALSE., __RC__ )
   CALL OVP_apply_mask( DATA_FOR_OVP_3D, OVP14_OUTPUT_3D, ChemEnv_STATE%MASK_2PM,  &
                        ChemEnv_STATE%OVP_FIRST_HMS, CURRENT_HMS, K_EDGES=.FALSE., __RC__ )
!
! T overpass
!
   CALL MAPL_GetPointer(IMPORT,  DATA_FOR_OVP_3D,         'T', __RC__)
   CALL MAPL_GetPointer(EXPORT,  OVP10_OUTPUT_3D,   'OVP10_T', __RC__)
   CALL MAPL_GetPointer(EXPORT,  OVP14_OUTPUT_3D,   'OVP14_T', __RC__)

   CALL OVP_apply_mask( DATA_FOR_OVP_3D, OVP10_OUTPUT_3D, ChemEnv_STATE%MASK_10AM, &
                        ChemEnv_STATE%OVP_FIRST_HMS, CURRENT_HMS, K_EDGES=.FALSE., __RC__ )
   CALL OVP_apply_mask( DATA_FOR_OVP_3D, OVP14_OUTPUT_3D, ChemEnv_STATE%MASK_2PM,  &
                        ChemEnv_STATE%OVP_FIRST_HMS, CURRENT_HMS, K_EDGES=.FALSE., __RC__ )
!
! PL overpass
!
   allocate(scratch_3d(idim,jdim,kdim), __STAT__ )

   ! We already have PLE
   scratch_3d = 0.5*(PLE(:,:,1:kdim)+PLE(:,:,0:kdim-1))
   DATA_FOR_OVP_3D => scratch_3d

   CALL MAPL_GetPointer(EXPORT,  OVP10_OUTPUT_3D,   'OVP10_PL', __RC__)
   CALL MAPL_GetPointer(EXPORT,  OVP14_OUTPUT_3D,   'OVP14_PL', __RC__)

   CALL OVP_apply_mask( DATA_FOR_OVP_3D, OVP10_OUTPUT_3D, ChemEnv_STATE%MASK_10AM, &
                        ChemEnv_STATE%OVP_FIRST_HMS, CURRENT_HMS, K_EDGES=.FALSE., __RC__ )
   CALL OVP_apply_mask( DATA_FOR_OVP_3D, OVP14_OUTPUT_3D, ChemEnv_STATE%MASK_2PM,  &
                        ChemEnv_STATE%OVP_FIRST_HMS, CURRENT_HMS, K_EDGES=.FALSE., __RC__ )

   deallocate(scratch_3d, __STAT__ )
!
! PLE overpass
!
   CALL MAPL_GetPointer(IMPORT,  DATA_FOR_OVP_3D,         'PLE', __RC__)
   CALL MAPL_GetPointer(EXPORT,  OVP10_OUTPUT_3D,   'OVP10_PLE', __RC__)
   CALL MAPL_GetPointer(EXPORT,  OVP14_OUTPUT_3D,   'OVP14_PLE', __RC__)

   CALL OVP_apply_mask( DATA_FOR_OVP_3D, OVP10_OUTPUT_3D, ChemEnv_STATE%MASK_10AM, &
                        ChemEnv_STATE%OVP_FIRST_HMS, CURRENT_HMS, K_EDGES=.TRUE., __RC__ )
   CALL OVP_apply_mask( DATA_FOR_OVP_3D, OVP14_OUTPUT_3D, ChemEnv_STATE%MASK_2PM,  &
                        ChemEnv_STATE%OVP_FIRST_HMS, CURRENT_HMS, K_EDGES=.TRUE., __RC__ )
!
! QV_VMR overpass
!
   allocate(scratch_3d(idim,jdim,kdim), __STAT__ )

   CALL MAPL_GetPointer(IMPORT,  DATA_FOR_OVP_3D,              'Q', __RC__)
   CALL MAPL_GetPointer(EXPORT,  OVP10_OUTPUT_3D,   'OVP10_QV_VMR', __RC__)
   CALL MAPL_GetPointer(EXPORT,  OVP14_OUTPUT_3D,   'OVP14_QV_VMR', __RC__)

   scratch_3d = DATA_FOR_OVP_3D/(1.0 - DATA_FOR_OVP_3D)  ! mass mixing ratio wrt dry air
   DATA_FOR_OVP_3D => scratch_3d

   CALL OVP_apply_mask( DATA_FOR_OVP_3D, OVP10_OUTPUT_3D, ChemEnv_STATE%MASK_10AM, &
                        ChemEnv_STATE%OVP_FIRST_HMS, CURRENT_HMS, K_EDGES=.FALSE., &
                        SCALE=(MAPL_AIRMW/MAPL_H2OMW), __RC__ )
   CALL OVP_apply_mask( DATA_FOR_OVP_3D, OVP14_OUTPUT_3D, ChemEnv_STATE%MASK_2PM,  &
                        ChemEnv_STATE%OVP_FIRST_HMS, CURRENT_HMS, K_EDGES=.FALSE., &
                        SCALE=(MAPL_AIRMW/MAPL_H2OMW), __RC__ )

   deallocate(scratch_3d, __STAT__ )
!
! QLTOT overpass
!
   CALL MAPL_GetPointer(IMPORT,  DATA_FOR_OVP_3D,         'QLTOT', __RC__)
   CALL MAPL_GetPointer(EXPORT,  OVP10_OUTPUT_3D,   'OVP10_QLTOT', __RC__)
   CALL MAPL_GetPointer(EXPORT,  OVP14_OUTPUT_3D,   'OVP14_QLTOT', __RC__)

   CALL OVP_apply_mask( DATA_FOR_OVP_3D, OVP10_OUTPUT_3D, ChemEnv_STATE%MASK_10AM, &
                        ChemEnv_STATE%OVP_FIRST_HMS, CURRENT_HMS, K_EDGES=.FALSE., __RC__ )
   CALL OVP_apply_mask( DATA_FOR_OVP_3D, OVP14_OUTPUT_3D, ChemEnv_STATE%MASK_2PM,  &
                        ChemEnv_STATE%OVP_FIRST_HMS, CURRENT_HMS, K_EDGES=.FALSE., __RC__ )
!
! PS overpass
!
   CALL MAPL_GetPointer(IMPORT,  DATA_FOR_OVP_2D,         'PS', __RC__)
   CALL MAPL_GetPointer(EXPORT,  OVP10_OUTPUT_2D,   'OVP10_PS', __RC__)
   CALL MAPL_GetPointer(EXPORT,  OVP14_OUTPUT_2D,   'OVP14_PS', __RC__)

   CALL OVP_apply_mask( DATA_FOR_OVP_2D, OVP10_OUTPUT_2D, ChemEnv_STATE%MASK_10AM, &
                        ChemEnv_STATE%OVP_FIRST_HMS, CURRENT_HMS,  __RC__ )
   CALL OVP_apply_mask( DATA_FOR_OVP_2D, OVP14_OUTPUT_2D, ChemEnv_STATE%MASK_2PM,  &
                        ChemEnv_STATE%OVP_FIRST_HMS, CURRENT_HMS,  __RC__ )
!
! PPBL overpass
!
   CALL MAPL_GetPointer(IMPORT,  DATA_FOR_OVP_2D,         'PPBL', __RC__)
   CALL MAPL_GetPointer(EXPORT,  OVP10_OUTPUT_2D,   'OVP10_PPBL', __RC__)
   CALL MAPL_GetPointer(EXPORT,  OVP14_OUTPUT_2D,   'OVP14_PPBL', __RC__)

   CALL OVP_apply_mask( DATA_FOR_OVP_2D, OVP10_OUTPUT_2D, ChemEnv_STATE%MASK_10AM, &
                        ChemEnv_STATE%OVP_FIRST_HMS, CURRENT_HMS, __RC__ )
   CALL OVP_apply_mask( DATA_FOR_OVP_2D, OVP14_OUTPUT_2D, ChemEnv_STATE%MASK_2PM,  &
                        ChemEnv_STATE%OVP_FIRST_HMS, CURRENT_HMS, __RC__ )
!
! TROPP overpass
!
   CALL MAPL_GetPointer(IMPORT,  DATA_FOR_OVP_2D,         'TROPP', __RC__)
   CALL MAPL_GetPointer(EXPORT,  OVP10_OUTPUT_2D,   'OVP10_TROPP', __RC__)
   CALL MAPL_GetPointer(EXPORT,  OVP14_OUTPUT_2D,   'OVP14_TROPP', __RC__)

   CALL OVP_apply_mask( DATA_FOR_OVP_2D, OVP10_OUTPUT_2D, ChemEnv_STATE%MASK_10AM, &
                        ChemEnv_STATE%OVP_FIRST_HMS, CURRENT_HMS, __RC__ )
   CALL OVP_apply_mask( DATA_FOR_OVP_2D, OVP14_OUTPUT_2D, ChemEnv_STATE%MASK_2PM,  &
                        ChemEnv_STATE%OVP_FIRST_HMS, CURRENT_HMS, __RC__ )
!
! U10M overpass
!
   CALL MAPL_GetPointer(IMPORT,  DATA_FOR_OVP_2D,         'U10M', __RC__)
   CALL MAPL_GetPointer(EXPORT,  OVP10_OUTPUT_2D,   'OVP10_U10M', __RC__)
   CALL MAPL_GetPointer(EXPORT,  OVP14_OUTPUT_2D,   'OVP14_U10M', __RC__)

   CALL OVP_apply_mask( DATA_FOR_OVP_2D, OVP10_OUTPUT_2D, ChemEnv_STATE%MASK_10AM, &
                        ChemEnv_STATE%OVP_FIRST_HMS, CURRENT_HMS, __RC__ )
   CALL OVP_apply_mask( DATA_FOR_OVP_2D, OVP14_OUTPUT_2D, ChemEnv_STATE%MASK_2PM,  &
                        ChemEnv_STATE%OVP_FIRST_HMS, CURRENT_HMS, __RC__ )
!
! V10M overpass
!
   CALL MAPL_GetPointer(IMPORT,  DATA_FOR_OVP_2D,         'V10M', __RC__)
   CALL MAPL_GetPointer(EXPORT,  OVP10_OUTPUT_2D,   'OVP10_V10M', __RC__)
   CALL MAPL_GetPointer(EXPORT,  OVP14_OUTPUT_2D,   'OVP14_V10M', __RC__)

   CALL OVP_apply_mask( DATA_FOR_OVP_2D, OVP10_OUTPUT_2D, ChemEnv_STATE%MASK_10AM, &
                        ChemEnv_STATE%OVP_FIRST_HMS, CURRENT_HMS, __RC__ )
   CALL OVP_apply_mask( DATA_FOR_OVP_2D, OVP14_OUTPUT_2D, ChemEnv_STATE%MASK_2PM,  &
                        ChemEnv_STATE%OVP_FIRST_HMS, CURRENT_HMS, __RC__ )
!
! EPV overpass
!
!   CALL MAPL_GetPointer(IMPORT,  DATA_FOR_OVP_3D,         'EPV', __RC__)
!   CALL MAPL_GetPointer(EXPORT,  OVP10_OUTPUT_3D,   'OVP10_EPV', __RC__)
!   CALL MAPL_GetPointer(EXPORT,  OVP14_OUTPUT_3D,   'OVP14_EPV', __RC__)
!
!   CALL OVP_apply_mask( DATA_FOR_OVP_3D, OVP10_OUTPUT_3D, ChemEnv_STATE%MASK_10AM, &
!                        ChemEnv_STATE%OVP_FIRST_HMS, CURRENT_HMS, K_EDGES=.FALSE., __RC__ )
!   CALL OVP_apply_mask( DATA_FOR_OVP_3D, OVP14_OUTPUT_3D, ChemEnv_STATE%MASK_2PM,  &
!                        ChemEnv_STATE%OVP_FIRST_HMS, CURRENT_HMS, K_EDGES=.FALSE., __RC__ )


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
