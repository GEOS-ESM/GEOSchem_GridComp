#include "MAPL_Generic.h"

!=============================================================================
!BOP

! !MODULE: GEOS_ChemEnvGridCompMod -- Prepares Environment for GEOSchem

! !INTERFACE:

module GEOS_ChemEnvGridCompMod

! !USES:

  use ESMF
  use MAPL_Mod

  implicit none
  private

! !PUBLIC MEMBER FUNCTIONS:

  public SetServices

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

    Iam = 'SetServices'
    call ESMF_GridCompGet( GC, NAME=COMP_NAME, CONFIG=CF, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // 'SetServices'

! Register services for this component
! ------------------------------------
    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN, Run1, __RC__ ) 
    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN, Run2, __RC__ ) 

!BOS

! !IMPORT STATE:

    call MAPL_AddImportSpec(GC,                               &
         SHORT_NAME = 'PLE',                                       &
         LONG_NAME  = 'air_pressure',                              &
         UNITS      = 'Pa',                                        &
         DIMS       =  MAPL_DimsHorzVert,                          &
         VLOCATION  =  MAPL_VLocationEdge,                         &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                               &
         SHORT_NAME = 'TH',                                        &
         LONG_NAME  = 'potential_temperature',                     &
         UNITS      = 'K',                                         &
         DIMS       =  MAPL_DimsHorzVert,                          &
         VLOCATION  =  MAPL_VLocationCenter,                       &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME = 'Q',                                         &
        LONG_NAME  = 'specific_humidity',                         &
        UNITS      = 'kg kg-1',                                   &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,                        &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

!   Convective precip
!   -----------------
    call MAPL_AddImportSpec(GC,                               &
         SHORT_NAME='CN_PRCP',                                     &
         LONG_NAME ='Surface Conv. rain flux needed by land',      &
         UNITS     ='kg m-2 s-1',                                  &
         DEFAULT   = MAPL_UNDEF,                                   &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

!   Total precip
!   ------------
    call MAPL_AddImportSpec(GC,                               &
         SHORT_NAME='TPREC',                                       &
         LONG_NAME ='total_precipitation',                         &
         UNITS     ='kg m-2 s-1',                                  &
         DEFAULT   = MAPL_UNDEF,                                   &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,				 &
       SHORT_NAME	  = 'FRLAND',				 &
       LONG_NAME	  = 'fraction_of_land', 		 &
       UNITS		  = '1',				 &
       DIMS		  = MAPL_DimsHorzOnly,		     &
       VLOCATION	  = MAPL_VLocationNone,		     &
    						      RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,				 &
       SHORT_NAME	  = 'FRLANDICE',			 &
       LONG_NAME	  = 'fraction_of_land_ice',		 &
       UNITS		  = '1',				 &
       DIMS		  = MAPL_DimsHorzOnly,		     &
       VLOCATION	  = MAPL_VLocationNone,		     &
    						      RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,				 &
       SHORT_NAME	  = 'FROCEAN',  			 &
       LONG_NAME	  = 'fraction_of_ocean',		 &
       UNITS		  = '1',				 &
       DIMS		  = MAPL_DimsHorzOnly,		     &
       VLOCATION	  = MAPL_VLocationNone,		     &
    						      RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,				 &
       SHORT_NAME	  = 'FRACI',  			 &
       LONG_NAME          = 'ice_covered_fraction_of_tile',      &
       UNITS              = '1',                                 &
       DIMS               = MAPL_DimsHorzOnly,                   &
       VLOCATION          = MAPL_VLocationNone,                  &
    						      RC=STATUS  )
    VERIFY_(STATUS)


    call MAPL_AddImportSpec(GC,                 &
       SHORT_NAME     = 'TS',                   &
       LONG_NAME      = 'TS',                   &
       UNITS          = 'K',                    &
       DIMS           = MAPL_DimsHorzOnly,      &
       VLOCATION      = MAPL_VLocationNone,     RC=STATUS  )
    VERIFY_(STATUS)

! !EXPORT STATE:

!    AIRDENS: Provided for Children
!    ------------------------------
     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'AIRDENS',                      &
        LONG_NAME          = 'moist_air_density',            &
        UNITS              = 'kg m-3',                       &
        DIMS               = MAPL_DimsHorzVert,              &
        VLOCATION          = MAPL_VLocationCenter,  RC=STATUS)
     VERIFY_(STATUS)

!    Density of dry air
!    ------------------
     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'AIRDENS_DRYP',                 &
        LONG_NAME          = 'partial_dry_air_density',      &
        UNITS              = 'kg dry m-3 tot',               &
        DIMS               = MAPL_DimsHorzVert,              &
        VLOCATION          = MAPL_VLocationCenter,  RC=STATUS)
     VERIFY_(STATUS)

!    DELP (This should be wired from DYN)
!    ------------------------------------
     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'DELP',                         &
        LONG_NAME          = 'pressure_thickness',           &
        UNITS              = 'Pa',                           &
        DIMS               = MAPL_DimsHorzVert,              &
        VLOCATION          = MAPL_VLocationCenter,  RC=STATUS)
     VERIFY_(STATUS)

!   Total precip
!   ------------
    call MAPL_AddExportSpec(GC,                              &
         SHORT_NAME        = 'TPREC',                        &
         LONG_NAME         = 'total_precipitation',          &
         UNITS             = 'kg m-2 s-1',                   &
         DIMS              = MAPL_DimsHorzOnly,              &
         VLOCATION         = MAPL_VLocationNone,    RC=STATUS)
    VERIFY_(STATUS)

!    Convective precip
!    -----------------
     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'CN_PRCP',                      &
        LONG_NAME          = 'Convective precipitation',     &
        UNITS              = 'kg m-2 s-1',                   &
        DIMS               = MAPL_DimsHorzOnly,              &
        VLOCATION          = MAPL_VLocationNone,    RC=STATUS)
     VERIFY_(STATUS)

!    Non-convective precip
!    ---------------------
     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'NCN_PRCP',                     &
        LONG_NAME          = 'Non-convective precipitation', &
        UNITS              = 'kg m-2 s-1',                   &
        DIMS               = MAPL_DimsHorzOnly,              &
        VLOCATION          = MAPL_VLocationNone,    RC=STATUS)
     VERIFY_(STATUS)
!EOS


! Create children's gridded components and invoke their SetServices
! -----------------------------------------------------------------
    call MAPL_GenericSetServices    ( GC, RC=STATUS )
    VERIFY_(STATUS)

    RETURN_(ESMF_SUCCESS)
  
  end subroutine SetServices

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

! Imports
  real, pointer, dimension(:,:,:)      :: pe => null()

! Exports
  real, pointer, dimension(:,:,:)      :: delp => null()

!=============================================================================
 
    integer                            :: k, k0

! Begin... 

!   Get the target components name and set-up traceback handle.
!   -----------------------------------------------------------
    call ESMF_GridCompGet ( GC, name=COMP_NAME, __RC__ )
    Iam = trim(COMP_NAME) // "Run1"

!   Get to the imports...
!   ---------------------
    call MAPL_GetPointer ( IMPORT,  pe,  'PLE', RC=STATUS );  VERIFY_(STATUS)

!   Get to the exports...
!   ---------------------
    call MAPL_GetPointer ( EXPORT, delp,   'DELP',        RC=STATUS ); VERIFY_(STATUS)

!   Compute moist (rho) and dry (rhoDry) air density
!   ------------------------------------------------
    CALL AirDens ( GC, IMPORT, EXPORT, __RC__ )

    if ( associated(delp) ) then
       k0 = 1-lbound(pe,3)

       do k = 1, size(delp,3)
          delp(:,:,k) = pe(:,:,k-k0+1) - pe(:,:,k-k0)
       end do
    end if


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
  real, pointer, dimension(:,:,:)      :: pe => null()

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


! Begin... 

!   Get the target components name and set-up traceback handle.
!   -----------------------------------------------------------
    call ESMF_GridCompGet ( GC, name=COMP_NAME, Grid=GRID, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // "Run2"

! Get my internal MAPL_Generic state
!-----------------------------------
    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS )
    VERIFY_(STATUS)

!   Get to the imports...
!   ---------------------
    call MAPL_GetPointer ( IMPORT,  pe,  'PLE', RC=STATUS );  VERIFY_(STATUS)

!   Get to the exports...
!   ---------------------
    call MAPL_GetPointer ( EXPORT, delp,   'DELP',        RC=STATUS ); VERIFY_(STATUS)

!   Compute moist (rho) and dry (rhoDry) air density
!   ------------------------------------------------
    CALL AirDens ( GC, IMPORT, EXPORT, __RC__ )

    if ( associated(delp) ) then
       k0 = 1-lbound(pe,3)

       do k = 1, size(delp,3)
          delp(:,:,k) = pe(:,:,k-k0+1) - pe(:,:,k-k0)
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
  real, pointer, dimension(:,:,:)      :: pe => null()
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
    call MAPL_GetPointer ( IMPORT,  pe,  'PLE', __RC__ )
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
    k0  = 1-lbound(pe,3)

    allocate(npk(iml,jml,nl+1),stat=STATUS) ! work space
    VERIFY_(STATUS)

    eps = MAPL_RVAP / MAPL_RGAS - 1.0

  !   Compute normalized pe**Kappa
  !   ----------------------------
    npk = (pe/MAPL_P00)**MAPL_KAPPA

  !   Compute rho from hydrostatic equation
  !   -------------------------------------
    DO k = 1, nl

        IF(ASSOCIATED(rho)) THEN

         rho(:,:,k) =       ( pe(:,:,k-k0+1) - pe(:,:,k-k0) ) /      &
                      ( MAPL_CP * ( th(:,:,k)*(1. + eps*q(:,:,k) ) ) &
                              * ( npk(:,:,k+1) - npk(:,:,k) ) )

        END IF

        IF(ASSOCIATED(rhoDry)) THEN

         rhoDry(:,:,k) =    ( 1.0 - q(:,:,k)                ) *      &
                            ( pe(:,:,k-k0+1) - pe(:,:,k-k0) ) /      &
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
