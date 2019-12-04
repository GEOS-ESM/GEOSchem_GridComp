#include "MAPL_Generic.h"

!=============================================================================
!BOP

! !MODULE: SSng_GridCompMod - GOCART refactoring of the SS gridded component 

! !INTERFACE:

module SSng_GridCompMod

! !USES:
   USE ESMF
   USE MAPL_Mod

   USE Chem_MieMod

   implicit none
   private

   type(Chem_Mie), dimension(2), save :: gocartngMieTable
   integer, parameter :: instanceComputational = 1
   integer, parameter :: instanceData          = 2

! !PUBLIC MEMBER FUNCTIONS:
   PUBLIC  SetServices
!   PUBLIC  Initialize



! !DESCRIPTION: This module implements GOCARTS' Sea Salt (SS) Gridded Component.

! !REVISION HISTORY:
! 24Oct2019  E.Sherman  First attempt at refactoring.

!EOP
!===========================================================================

contains


!BOP

! !IROUTINE: SetServices 

! !INTERFACE:

  subroutine SetServices ( GC, RC )

! !ARGUMENTS:

    type (ESMF_GridComp), intent(INOUT)   :: GC  ! gridded component
    integer,              intent(  OUT)   :: RC  ! return code


! !DESCRIPTION: This version uses the {\tt GEOS\_GenericSetServices}, which sets
!   the Initialize and Finalize services to generic versions. It also
!   allocates our instance of a generic state and puts it in the 
!   gridded component (GC). Here we only set the two-stage run method and
!   declare the data services.
! \newline

! !REVISION HISTORY: 
! 24oct2019   E.Sherman  First attempt at refactoring

!EOP

!****************************************************************************
!
! ErrLog Variables

    character (len=ESMF_MAXSTR)                 :: IAm
    integer                                     :: STATUS
    character (len=ESMF_MAXSTR)                 :: COMP_NAME


    type (ESMF_Config)                          :: cfg

    character (len=ESMF_MAXSTR)                 :: field_name

    integer                                     :: n, i, n_bins
    real                                        :: DEFVAL
    logical                                     :: data_driven=.true.

    !development testing variables - to be deleted
    real, dimension(:,:), pointer       :: ptr_test

!****************************************************************************
! Begin...

! Get my name and set-up traceback handle
! ---------------------------------------
    Iam = 'SetServices'
    call ESMF_GridCompGet (GC, NAME=COMP_NAME, __RC__)
    Iam = trim(COMP_NAME) // Iam

if (mapl_am_i_root()) print*,'GOCARTng SSng SetServices BEGIN'           ! for testing - to be deleted
if (mapl_am_i_root()) print*,'GOCARTng SS COMP_NAME = ', trim(COMP_NAME) ! for testing - to be deleted


!  Load resource file and get number of bins 
!  -------------------------------------------
   cfg = ESMF_ConfigCreate (__RC__)
   call ESMF_ConfigLoadFile (cfg, 'SSng_GridComp_'//trim(COMP_NAME)//'.rc', RC=STATUS)
   if (STATUS /= 0) then
     if (mapl_am_i_root()) print*,'SSngGridComp_'//trim(COMP_NAME)//'.rc does not exist! loading SSng_GridComp_SS.data.rc instead'
     call ESMF_ConfigLoadFile (cfg, 'SSng_GridComp_SS.rc', __RC__)
   end if

   call ESMF_ConfigGetAttribute (cfg, n_bins, label='bins:', __RC__)

!   Set entry points
!   ------------------------
   call MAPL_GridCompSetEntryPoint (GC, ESMF_METHOD_INITIALIZE,  Initialize, __RC__)
   call MAPL_GridCompSetEntryPoint (GC, ESMF_METHOD_RUN, Run1, __RC__)

!   Is SS data driven?
!   ------------------
    call data_driven_ (COMP_NAME, data_driven, __RC__)


!   INTERNAL STATE
!   ---------------
!   Default INTERNAL state values
!   -----------------------------
    DEFVAL = 0.0

!   Aerosol Tracers to be transported
!   ---------------------------------
    do i = 1, n_bins
        write(field_name, '(A, I0.3)') '::ss', i

        call MAPL_AddInternalSpec(GC,                                        &
          SHORT_NAME = trim(COMP_NAME)//trim(field_name),                    &
          LONG_NAME  = 'Sea Salt Mixing Ratio (bin '//trim(field_name)//')', &
          UNITS      = 'kg kg-1',                                            &
          RESTART    = MAPL_RestartOptional,                                 &
          DEFAULT    = DEFVAL,                                               &
          DIMS       = MAPL_DimsHorzVert,                                    &
          VLOCATION  = MAPL_VLocationCenter, __RC__)
    end do


!   IMPORT STATE
!   -------------
    if (data_driven) then
        do i = 1, n_bins
            write(field_name, '(A, I0.3)') '', i

            call MAPL_AddImportSpec(GC,                                           &
              SHORT_NAME = 'climss'//trim(field_name),                            &
              LONG_NAME  = 'Sea Salt Mixing Ratio (bin '//trim(field_name)//')',  &
              UNITS      = 'kg kg-1',                                             &
              RESTART    = MAPL_RestartSkip,                                      &
              DIMS       = MAPL_DimsHorzVert,                                     &
              VLOCATION  = MAPL_VLocationCenter, __RC__)

!           ! dry deposition
            call MAPL_AddImportSpec(GC,                                           &
              SHORT_NAME = 'climSSDP'//trim(field_name),                          &
              LONG_NAME  = 'Sea Salt Mixing Ratio (bin '//trim(field_name)//')',  &
              UNITS      = 'kg kg-1',                                             &
              DIMS       = MAPL_DimsHorzOnly,                                     &
              VLOCATION  = MAPL_VLocationCenter,                                  &
              RESTART    = MAPL_RestartSkip, __RC__)

!           ! wet deposition    
            call MAPL_AddImportSpec(GC,                                           &
               SHORT_NAME = 'climSSWT'//trim(field_name),                         &
               LONG_NAME  = 'Sea Salt Mixing Ratio (bin '//trim(field_name)//')', &
               UNITS      = 'kg kg-1',                                            &
               DIMS       = MAPL_DimsHorzOnly,                                    &
               VLOCATION  = MAPL_VLocationCenter,                                 &
               RESTART    = MAPL_RestartSkip, __RC__)

!           ! gravitational settling
            call MAPL_AddImportSpec(GC,                                           &
               SHORT_NAME = 'climSSSD'//trim(field_name),                         &
               LONG_NAME  = 'Sea Salt Mixing Ratio (bin '//trim(field_name)//')', &
               UNITS      = 'kg kg-1',                                            &
               DIMS       = MAPL_DimsHorzOnly,                                    &
               VLOCATION  = MAPL_VLocationCenter,                                 &
               RESTART    = MAPL_RestartSkip, __RC__)

!        ! convective scavenging
            call MAPL_AddImportSpec(GC,                                           &
               SHORT_NAME = 'climSSSV'//trim(field_name),                         &
               LONG_NAME  = 'Sea Salt Mixing Ratio (bin '//trim(field_name)//')', &
               UNITS      = 'kg kg-1',                                            &
               DIMS       = MAPL_DimsHorzOnly,                                    &
               VLOCATION  = MAPL_VLocationCenter,                                 &
               RESTART    = MAPL_RestartSkip, __RC__)
        end do
    end if


!   EXPORT STATE
!   -------------
    if (.not. data_driven) then
#       include "SSng_ExportSpec___.h"
    end if

!   This state holds fields needed by radiation
!   ---------------------------------------------
    call MAPL_AddExportSpec(GC,                                 &
      SHORT_NAME = trim(COMP_NAME)//'_AERO',                   &
!       SHORT_NAME = 'AERO',                   &
       LONG_NAME  = 'aerosols_from_'//trim(COMP_NAME),  &
       UNITS      = 'kg kg-1',                                  &
       DIMS       = MAPL_DimsHorzVert,                          &
       VLOCATION  = MAPL_VLocationCenter,                       &
       DATATYPE   = MAPL_StateItem, __RC__)

!   This state is needed by MOIST - It will contain aerosols
!   ----------------------------------------------------------
    call MAPL_AddExportSpec(GC,                                                  &
       SHORT_NAME = trim(COMP_NAME)//'_AERO_ACI',                                &
       LONG_NAME  = 'aerosol_cloud_interaction_aerosols_from_'//trim(COMP_NAME),  &
       UNITS      = 'kg kg-1',                                                   &
       DIMS       = MAPL_DimsHorzVert,                                           &
       VLOCATION  = MAPL_VLocationCenter,                                        &
       DATATYPE   = MAPL_StateItem, __RC__)

!   This bundle is needed by surface for snow albedo modification
!   by aerosol settling and deposition
!   DEVELOPMENT NOTE - Change to StateItem in future
!   ---------------------------------------------------------------
    call MAPL_AddExportSpec(GC,                                   &
       SHORT_NAME = 'AERO_DP',                  &
       LONG_NAME  = 'aerosol_deposition_from_'//trim(COMP_NAME),  &
       UNITS      = 'kg m-2 s-1',                                 &
       DIMS       = MAPL_DimsHorzOnly,                            &
       DATATYPE   = MAPL_BundleItem, __RC__)


if (mapl_am_i_root()) print*,'GOCARTng SSng SetServices END' ! for testing - to be deleted

!   Set generic services
!   ----------------------------------
    call MAPL_GenericSetServices (GC, __RC__)

    RETURN_(ESMF_SUCCESS)

  end subroutine SetServices

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!BOP

! !IROUTINE: Initialize 

! !INTERFACE:

  subroutine Initialize (GC, IMPORT, EXPORT, CLOCK, RC)

! !ARGUMENTS:

    type (ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
    type (ESMF_State),    intent(inout) :: IMPORT ! Import state
    type (ESMF_State),    intent(inout) :: EXPORT ! Export state
    type (ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
    integer, optional,    intent(  out) :: RC     ! Error code

! !DESCRIPTION: This initializes SS' Grid Component. It primaryily fills 
!               GOCART's AERO states with its sea salt fields. 

! !REVISION HISTORY: 
! 24oct2019   E.Sherman  First attempt at refactoring

!EOP

!****************************************************************************
! ErrLog Variables

    character (len=ESMF_MAXSTR)          :: IAm
    integer                              :: STATUS
    character (len=ESMF_MAXSTR)          :: COMP_NAME

! Local derived type aliases

    type (MAPL_MetaComp),      pointer  :: MAPL
    type (ESMF_Grid)                    :: grid
    type (ESMF_State)                   :: INTERNAL
    type (ESMF_State)                   :: AERO, AERO_ACI
    type (ESMF_State)                   :: providerState
    type (ESMF_Config)                  :: cfg
    type (ESMF_FieldBundle)             :: Bundle_DP

    integer                             :: i, n_bins, nbins_test
    integer                              :: instance
    type (ESMF_Field)                   :: field, fld
    character (len=ESMF_MAXSTR)         :: field_name, prefix

    logical                             :: data_driven


    !development testing variables - to be deleted
    real, dimension(:,:,:), pointer       :: ptr_test
    real, dimension(:,:), pointer       :: ptr_test2d

!****************************************************************************

! Begin... 

!   Get the target components name and set-up traceback handle.
!   -----------------------------------------------------------
    Iam = "Initialize"
    call ESMF_GridCompGet (GC, grid=grid, name=COMP_NAME, __RC__)
    Iam = trim(COMP_NAME) // trim(Iam)


if (mapl_am_i_root()) print*,'GOCART SSng Initialize BEGIN' ! for testing - to be deleted

!   Get my internal MAPL_Generic state
!   -----------------------------------
    call MAPL_GetObjectFromGC (GC, MAPL, __RC__)

!  Load resource file and get number of bins 
!  -------------------------------------------
    cfg = ESMF_ConfigCreate (__RC__)
    call ESMF_ConfigLoadFile (cfg, 'SSng_GridComp_'//trim(COMP_NAME)//'.rc', RC=STATUS)
    if (STATUS /= 0) then
      if (mapl_am_i_root()) print*,'SSng_GridComp_'//trim(COMP_NAME)//'.rc does not exist! loading SSng_GridComp_SS.rc instead'
      call ESMF_ConfigLoadFile( cfg, 'SSng_GridComp_SS.rc', __RC__)
    end if

    call ESMF_ConfigGetAttribute (cfg, n_bins, label='bins:', __RC__)

!   Call Generic Initialize 
!   ----------------------------------------
    call MAPL_GenericInitialize (GC, IMPORT, EXPORT, CLOCK, __RC__)

!   Get parameters from generic state.
!   -----------------------------------
    call MAPL_Get (MAPL, INTERNAL_ESMF_STATE=INTERNAL, __RC__)

!   Is SS data driven?
!   ------------------
    call data_driven_ (COMP_NAME, data_driven, __RC__)


!   If using GOCART.data, the data is provided in the import
!   state via ExtData versus the actual GOCART children
!   -------------------------------------------------------- 
    if ( data_driven ) then
       providerState = IMPORT
       prefix = 'clim'
    else
       providerState = EXPORT
       prefix = ''
    end if


!   Fill AERO State with sea salt fields
!   ----------------------------------------
    call ESMF_StateGet (EXPORT, trim(COMP_NAME)//'_AERO'    , AERO    , __RC__)
!    call ESMF_StateGet (EXPORT, 'AERO'    , AERO    , __RC__)
    call ESMF_StateGet (EXPORT, trim(COMP_NAME)//'_AERO_ACI', AERO_ACI, __RC__)
    call ESMF_StateGet (EXPORT, 'AERO_DP' , Bundle_DP, __RC__)


!   Set number of bins for this instance
    call ESMF_AttributeSet (AERO, name='n_bins', value=n_bins, __RC__)

!   Fill EXPORT States and Bundle
!   ------------------------------
    do i = 1, n_bins
        write(field_name, '(A, I0.3)') '', i
        call ESMF_StateGet (INTERNAL, trim(COMP_NAME)//'::ss'//trim(field_name), field, __RC__)
        fld = MAPL_FieldCreate (field, name='ss'//trim(field_name), __RC__)
        call MAPL_StateAdd (AERO,  fld, __RC__)
        call MAPL_StateAdd (AERO_ACI, fld, __RC__)

!       Dry deposition
!       ---------------
        call ESMF_StateGet (providerState, trim(prefix)//'SSDP'//trim(field_name), field, __RC__)
        call MAPL_AllocateCoupling (field, __RC__)
        call MAPL_FieldBundleAdd (Bundle_DP, field, __RC__)

!       Wet deposition (Convective scavenging)
!       --------------------------------------
        call ESMF_StateGet (providerState, trim(prefix)//'SSSV'//trim(field_name), field, __RC__)
        call MAPL_AllocateCoupling (field, __RC__)
        call MAPL_FieldBundleAdd (Bundle_DP, field, __RC__)

!       Wet deposition
!       ---------------
        call ESMF_StateGet (providerState, trim(prefix)//'SSWT'//trim(field_name), field, __RC__)
        call MAPL_AllocateCoupling (field, __RC__)
        call MAPL_FieldBundleAdd (Bundle_DP, field, __RC__)

!       Gravitational Settling
!       ----------------------
        call ESMF_StateGet (providerState, trim(prefix)//'SSSD'//trim(field_name), field, __RC__)
        call MAPL_AllocateCoupling (field, __RC__)
        call MAPL_FieldBundleAdd (Bundle_DP, field, __RC__)
    end do

!   Set AERO States' attributes
!   ----------------------------
    if (data_driven) then
        instance = instanceData
    else
        instance = instanceComputational
    end if


    gocartngMieTable(instance) = Chem_MieCreateng(cfg, COMP_NAME, __RC__)


    ! Mie Table instance/index
    call ESMF_AttributeSet(AERO, name='mie_table_instance', value=instance, __RC__)

    ! state of the atmosphere
    call ESMF_AttributeSet(AERO, name='air_pressure_for_aerosol_optics',             value='PLE', __RC__)
    call ESMF_AttributeSet(AERO, name='relative_humidity_for_aerosol_optics',        value='RH',  __RC__)

    ! aerosol optics
    call ESMF_AttributeSet(AERO, name='band_for_aerosol_optics',                     value=0,     __RC__)
    call ESMF_AttributeSet(AERO, name='extinction_in_air_due_to_ambient_aerosol',    value='EXT', __RC__)
    call ESMF_AttributeSet(AERO, name='single_scattering_albedo_of_ambient_aerosol', value='SSA', __RC__)
    call ESMF_AttributeSet(AERO, name='asymmetry_parameter_of_ambient_aerosol',      value='ASY', __RC__)

    ! add PLE to aero state
    call ESMF_AttributeGet(AERO, name='air_pressure_for_aerosol_optics', value=field_name, __RC__)
    if (field_name /= '') then
        field = MAPL_FieldCreateEmpty(trim(field_name), grid, __RC__)

        call MAPL_FieldAllocCommit(field, dims=MAPL_DimsHorzVert, location=MAPL_VLocationEdge, typekind=MAPL_R4, hw=0, __RC__)
        call MAPL_StateAdd(AERO, field, __RC__)
    end if

    ! add RH to Aero state
    call ESMF_AttributeGet(AERO, name='relative_humidity_for_aerosol_optics', value=field_name, __RC__)
    if (field_name /= '') then
        field = MAPL_FieldCreateEmpty(trim(field_name), grid, __RC__)

        call MAPL_FieldAllocCommit(field, dims=MAPL_DimsHorzVert, location=MAPL_VLocationCenter, typekind=MAPL_R4, hw=0, __RC__)
        call MAPL_StateAdd(AERO, field, __RC__)
    end if

    ! add EXT to aero state
    call ESMF_AttributeGet(AERO, name='extinction_in_air_due_to_ambient_aerosol', value=field_name, __RC__)
    if (field_name /= '') then
        field = MAPL_FieldCreateEmpty(trim(field_name), grid, __RC__)

        call MAPL_FieldAllocCommit(field, dims=MAPL_DimsHorzVert, location=MAPL_VLocationCenter, typekind=MAPL_R4, hw=0, __RC__)
        call MAPL_StateAdd(AERO, field, __RC__)
    end if

    ! add SSA to aero state
    call ESMF_AttributeGet(AERO, name='single_scattering_albedo_of_ambient_aerosol', value=field_name, __RC__)
    if (field_name /= '') then
        field = MAPL_FieldCreateEmpty(trim(field_name), grid, __RC__)

        call MAPL_FieldAllocCommit(field, dims=MAPL_DimsHorzVert, location=MAPL_VLocationCenter, typekind=MAPL_R4, hw=0, __RC__)
        call MAPL_StateAdd(AERO, field, __RC__)
    end if

    ! add ASY to aero state
    call ESMF_AttributeGet(AERO, name='asymmetry_parameter_of_ambient_aerosol', value=field_name, RC=STATUS)
    if (field_name /= '') then
        field = MAPL_FieldCreateEmpty(trim(field_name), grid, __RC__)

        call MAPL_FieldAllocCommit(field, dims=MAPL_DimsHorzVert, location=MAPL_VLocationCenter, typekind=MAPL_R4, hw=0, __RC__)
        call MAPL_StateAdd(AERO, field, __RC__)
    end if

    ! attach the aerosol optics method
    call ESMF_MethodAdd(AERO, label='aerosol_optics', userRoutine=aerosol_optics, __RC__)


if (mapl_am_i_root()) print*,'GOCART SSng Initialize END'  ! for testing - to be deleted



    RETURN_(ESMF_SUCCESS)

  end subroutine Initialize


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!BOP
! !IROUTINE: Run1 

! !INTERFACE:

  subroutine Run1 (GC, IMPORT, EXPORT, CLOCK, RC)

    ! !ARGUMENTS:

    type (ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
    type (ESMF_State),    intent(inout) :: IMPORT ! Import state
    type (ESMF_State),    intent(inout) :: EXPORT ! Export state
    type (ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
    integer, optional,    intent(  out) :: RC     ! Error code:

! !DESCRIPTION: Run method for the Sea Salt Grid Component. Determines whether to run
!               data or computational run method.

!EOP

!****************************************************************************
! ErrLog Variables

    character (len=ESMF_MAXSTR)          :: IAm
    integer                              :: STATUS
    character (len=ESMF_MAXSTR)          :: COMP_NAME


    type (MAPL_MetaComp), pointer        :: MAPL
    type (ESMF_State)                    :: INTERNAL
 
    logical                              :: data_driven

!*****************************************************************************
! Begin... 

!   Get my name and set-up traceback handle
!   ---------------------------------------
    Iam = 'Run1'
    call ESMF_GridCompGet (GC, NAME=COMP_NAME, __RC__)
    Iam = trim(COMP_NAME) // Iam

if (mapl_am_i_root()) print*,'GOCARTng SSng Run1 BEGIN'

!   Get my internal MAPL_Generic state
!   -----------------------------------
    call MAPL_GetObjectFromGC (GC, MAPL, __RC__)

!   Get parameters from generic state.
!   -----------------------------------
    call MAPL_Get (MAPL, INTERNAL_ESMF_STATE=INTERNAL, __RC__)

!   Is DU data driven?
!   ------------------
    call data_driven_ (COMP_NAME, data_driven, __RC__)

!   Update INTERNAL state variables with ExtData
!   ---------------------------------------------
    if (data_driven) then
      call Run_data (GC, IMPORT, EXPORT, INTERNAL, __RC__)
    end if


if (mapl_am_i_root()) print*,'GOCARTng SSng Run1 END'


    RETURN_(ESMF_SUCCESS)

  end subroutine Run1


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!BOP
! !IROUTINE: Run_data -- ExtData Sea Salt Grid Component

! !INTERFACE:

  subroutine Run_data (GC, IMPORT, EXPORT, INTERNAL, RC)

    ! !ARGUMENTS:

    type (ESMF_GridComp), intent(inout) :: GC       ! Gridded component 
    type (ESMF_State),    intent(inout) :: IMPORT   ! Import state
    type (ESMF_State),    intent(inout) :: EXPORT   ! Export state
    type (ESMF_State),    intent(inout) :: INTERNAL ! Interal state
    integer, optional,    intent(  out) :: RC       ! Error code:

! !DESCRIPTION: Updates pointers in Internal state with fields from ExtData. 

!EOP

!****************************************************************************
! ErrLog Variables

    character (len=ESMF_MAXSTR)        :: IAm
    integer                            :: STATUS
    character (len=ESMF_MAXSTR)        :: COMP_NAME


    type (MAPL_MetaComp), pointer      :: MAPL
    type (ESMF_State)                  :: AERO

    integer                            :: i, n_bins
    character (len=ESMF_MAXSTR)        :: field_name

    real, pointer, dimension(:,:,:)    :: ptr3d_int
    real, pointer, dimension(:,:,:)    :: ptr3d_imp


!   test vars - to be deleted
    type (ESMF_Field)                  :: field
    real, pointer, dimension(:,:)      :: ptr_test2d

!*****************************************************************************
! Begin... 

! Get my name and set-up traceback handle
! ---------------------------------------
    Iam = 'Run1'
    call ESMF_GridCompGet (GC, NAME=COMP_NAME, __RC__)
    Iam = trim(COMP_NAME) // Iam


if (mapl_am_i_root()) print*,'SSng ExtData Run_data BEGIN'


    call ESMF_StateGet (EXPORT, trim(COMP_NAME)//'_AERO', AERO, __RC__)
!    call ESMF_StateGet (EXPORT, 'AERO', AERO, __RC__)
    call ESMF_AttributeGet (AERO, name='n_bins', value=n_bins, __RC__)

!   Update interal data pointers with ExtData
!   -----------------------------------------
    do i = 1, n_bins
    write(field_name, '(A, I0.3)') 'ss', i
        call MAPL_GetPointer (INTERNAL, NAME=trim(COMP_NAME)//'::'//trim(field_name), ptr=ptr3d_int, __RC__)
        call MAPL_GetPointer (IMPORT,  NAME='clim'//trim(field_name), ptr=ptr3d_imp, __RC__)

        ptr3d_int = ptr3d_imp
    end do


if (mapl_am_i_root()) print*,'SSng ExtData Run_data END'
 


    RETURN_(ESMF_SUCCESS)

  end subroutine Run_data


!-------------------------------------------------------------------------------------
  subroutine aerosol_optics(state, rc)

    implicit none

!   !ARGUMENTS:
    type (ESMF_State)                                :: state
    integer,            intent(out)                  :: rc

!   !Local
    real, dimension(:,:,:), pointer                  :: ple
    real, dimension(:,:,:), pointer                  :: rh
    real, dimension(:,:,:), pointer                  :: var
    real, dimension(:,:,:), pointer                  :: q
    real, dimension(:,:,:,:), pointer                :: q_4d

    character (len=ESMF_MAXSTR)                      :: fld_name
    type(ESMF_Field)                                 :: fld
    character (len=ESMF_MAXSTR)                      :: COMP_NAME

    real, dimension(:,:,:,:), allocatable            :: ext, ssa, asy  ! (lon:,lat:,lev:,band:)

    integer                                          :: instance
    integer                                          :: n
    integer                                          :: i1, j1, i2, j2, km
    integer                                          :: band, offset
    integer                                          :: n_bins
    integer, parameter                               :: n_bands = 1

    real    :: x
    integer :: i, j, k

    __Iam__('SSng::aerosol_optics')


!   Begin... 


!   Mie Table instance/index
!   ------------------------
    call ESMF_AttributeGet(state, name='mie_table_instance', value=instance, __RC__)

!   Radiation band
!   --------------
    band = 0
    call ESMF_AttributeGet(state, name='band_for_aerosol_optics', value=band, __RC__)
    offset = band - n_bands

!   Pressure at layer edges 
!   ------------------------
    call ESMF_AttributeGet(state, name='air_pressure_for_aerosol_optics', value=fld_name, __RC__)
    call MAPL_GetPointer(state, ple, trim(fld_name), __RC__)

    i1 = lbound(ple, 1); i2 = ubound(ple, 1)
    j1 = lbound(ple, 2); j2 = ubound(ple, 2)
                         km = ubound(ple, 3)

!   Relative humidity
!   -----------------
    call ESMF_AttributeGet(state, name='relative_humidity_for_aerosol_optics', value=fld_name, __RC__)
    call MAPL_GetPointer(state, rh, trim(fld_name), __RC__)

    call ESMF_AttributeGet (state, name='n_bins'   , value=n_bins   , __RC__)


    allocate(ext(i1:i2,j1:j2,km,n_bands), &
             ssa(i1:i2,j1:j2,km,n_bands), &
             asy(i1:i2,j1:j2,km,n_bands), __STAT__)

    allocate(q_4d(i1:i2,j1:j2,km,n_bins), __STAT__)

    do n = 1, n_bins
        write (fld_name, '(A, I0.3)') 'ss', n
        call ESMF_StateGet(state, trim(fld_name), field=fld, __RC__)
        call ESMF_FieldGet(fld, farrayPtr=q, __RC__)

        do k = 1, km
            do j = j1, j2
                do i = i1, i2
                    x = ((PLE(i,j,k) - PLE(i,j,k-1))*0.01)*(100./MAPL_GRAV)
                    q_4d(i,j,k,n) = x * q(i,j,k)
                end do
            end do
        end do
    end do


    call mie_(gocartngMieTable(instance), n_bins, n_bands, offset, q_4d, rh, ext, ssa, asy, __RC__)


    call ESMF_AttributeGet(state, name='extinction_in_air_due_to_ambient_aerosol', value=fld_name, __RC__)
    if (fld_name /= '') then
        call MAPL_GetPointer(state, var, trim(fld_name), __RC__)
        var = ext(:,:,:,1)
    end if

    call ESMF_AttributeGet(state, name='single_scattering_albedo_of_ambient_aerosol', value=fld_name, __RC__)
    if (fld_name /= '') then
        call MAPL_GetPointer(state, var, trim(fld_name), __RC__)
        var = ssa(:,:,:,1)
    end if

   call ESMF_AttributeGet(state, name='asymmetry_parameter_of_ambient_aerosol', value=fld_name, __RC__)
    if (fld_name /= '') then
        call MAPL_GetPointer(state, var, trim(fld_name), __RC__)
        var = asy(:,:,:,1)
    end if

if (mapl_am_i_root()) print*,'SSng aerosol_optics IM RUNNING_ ='



    deallocate(ext, ssa, asy, __STAT__)


    RETURN_(ESMF_SUCCESS)

  contains

    subroutine mie_(mie_table, n_bins, nb, offset, q, rh, ext, ssa, asy, rc)
    implicit none

    type(Chem_Mie),    intent(inout):: mie_table    ! mie table
    integer,           intent(in )  :: n_bins           ! number of bins
    integer,           intent(in )  :: nb               ! number of bands
    integer,           intent(in )  :: offset           ! bands offset 
    real,              intent(in )  :: q(:,:,:,:)       ! aerosol mass mixing ratio, kg kg-1
    real,              intent(in )  :: rh(:,:,:)        ! relative humidity

    real,              intent(out)  :: ext(:,:,:,:)     ! extinction
    real,              intent(out)  :: ssa(:,:,:,:)     ! SSA
    real,              intent(out)  :: asy(:,:,:,:)     ! asymmetry parameter

    integer,           intent(out)  :: rc

    ! local
    integer :: STATUS
    character(len=ESMF_MAXSTR) :: Iam='SSng::aerosol_optics::mie_'
    character (len=ESMF_MAXSTR)                      :: fld_name

    integer :: l, idx, na

     real(kind=8) :: ext_(size(ext,1),size(ext,2),size(ext,3),size(ext,4))
     real(kind=8) :: ssa_(size(ext,1),size(ext,2),size(ext,3),size(ext,4))
     real(kind=8) :: asy_(size(ext,1),size(ext,2),size(ext,3),size(ext,4))


     ASSERT_ (n_bins == size(q,4))

     ext_ = 0.0d0
     ssa_ = 0.0d0
     asy_ = 0.0d0

    do l = 1, n_bins
        write (fld_name, '(A, I0.3)') 'ss', l
        idx = Chem_MieQueryIdx(mie_table, trim(fld_name), __RC__)
        call Chem_MieQueryAllBand4D(mie_table, idx, nb, offset, q(:,:,:,l), rh, ext, ssa, asy, __RC__)

        ext_ = ext_ +          ext     ! total extinction due to dust
        ssa_ = ssa_ +     (ssa*ext)    ! total scattering
        asy_ = asy_ + asy*(ssa*ext)    ! sum of (asy * sca)


!if (l == 2) then
!  if (mapl_am_i_root()) print*,'GOCARTng du001 nb = ', nb
!  if (mapl_am_i_root()) print*,'GOCARTng ss002 offset = ', offset
!  if (mapl_am_i_root()) print*,'GOCARTng ext = ', ext
!   if (mapl_am_i_root()) print*,'GOCARTng q(:,:,:,1) = ', q(:,:,:,1)
!   if (mapl_am_i_root()) print*,'GOCARTng rh = ', rh
!end if

    end do

    ext = ext_
    ssa = ssa_
    asy = asy_


    RETURN_(ESMF_SUCCESS)

    end subroutine mie_


  end subroutine aerosol_optics




!---------------------------------------------------------------------------------------

  subroutine data_driven_(COMP_NAME, data_driven, RC)

    !ARGUMENTS:
    integer, optional,              intent(  out)   :: RC          ! Error code:
    character(len=ESMF_MAXSTR),     intent(in   )   :: COMP_NAME
    logical,                        intent(  out)   :: data_driven
    integer                                        :: i

!   Description: Determines whether Sea Salt is data driven or not.


     __Iam__('data_driven_')

!   Is DU data driven?
!   ------------------
    data_driven = .false.

    i = index(COMP_NAME, 'data')
    if (i > 0) then
      data_driven = .true.
    end if

    RETURN_(ESMF_SUCCESS)
    end subroutine data_driven_

!-----------------------------------------------------------------------------------










end module SSng_GridCompMod

