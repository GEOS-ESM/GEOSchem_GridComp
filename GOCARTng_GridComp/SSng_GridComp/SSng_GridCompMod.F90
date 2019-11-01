#include "MAPL_Generic.h"

!=============================================================================
!BOP

! !MODULE: SSng_GridCompMod - GOCART refactoring of the SS gridded component 

! !INTERFACE:

module SSng_GridCompMod

! !USES:
   USE ESMF
   USE MAPL_Mod

   implicit none
   private

! !PUBLIC MEMBER FUNCTIONS:
   PUBLIC  SetServices
   PUBLIC  Initialize



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
       LONG_NAME  = 'seasalt_aerosols_from_'//trim(COMP_NAME),  &
       UNITS      = 'kg kg-1',                                  &
       DIMS       = MAPL_DimsHorzVert,                          &
       VLOCATION  = MAPL_VLocationCenter,                       &
       DATATYPE   = MAPL_StateItem, __RC__)

!   This state is needed by MOIST - It will contain aerosols
!   ----------------------------------------------------------
    call MAPL_AddExportSpec(GC,                                                  &
       SHORT_NAME = trim(COMP_NAME)//'_AERO_ACI',                                &
       LONG_NAME  = 'aerosol_cloud_interaction_seasalt_from_'//trim(COMP_NAME),  &
       UNITS      = 'kg kg-1',                                                   &
       DIMS       = MAPL_DimsHorzVert,                                           &
       VLOCATION  = MAPL_VLocationCenter,                                        &
       DATATYPE   = MAPL_StateItem, __RC__)

!   This bundle is needed by surface for snow albedo modification
!   by aerosol settling and deposition
!   DEVELOPMENT NOTE - Change to StateItem in future
!   ---------------------------------------------------------------
    call MAPL_AddExportSpec(GC,                                   &
       SHORT_NAME = trim(COMP_NAME)//'_AERO_DP',                  &
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
    type (ESMF_Grid)                    :: GRID
    type (ESMF_State)                   :: INTERNAL
    type (ESMF_State)                   :: SS_State, SS_State_ACI
    type (ESMF_State)                   :: providerState
    type (ESMF_Config)                  :: cfg
    type (ESMF_FieldBundle)             :: SS_Bundle_DP

    integer                             :: i, n_bins, nbins_test
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
    call ESMF_GridCompGet (GC, name=COMP_NAME, __RC__)
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
    call ESMF_StateGet (EXPORT, trim(COMP_NAME)//'_AERO'    , SS_State    , __RC__)
    call ESMF_StateGet (EXPORT, trim(COMP_NAME)//'_AERO_ACI', SS_State_ACI, __RC__)
    call ESMF_StateGet (EXPORT, trim(COMP_NAME)//'_AERO_DP' , SS_Bundle_DP, __RC__)

!   Set number of bins for this instance
    call ESMF_AttributeSet (SS_State, name='n_bins', value=n_bins, __RC__)

    do i = 1, n_bins
        write(field_name, '(A, I0.3)') '', i
        call ESMF_StateGet (INTERNAL, trim(COMP_NAME)//'::ss'//trim(field_name), field, __RC__)
        fld = MAPL_FieldCreate (FIELD, name='ss'//trim(field_name), __RC__)
        call MAPL_StateAdd (SS_State,  field, __RC__)
        call MAPL_StateAdd (SS_State_ACI, field, __RC__)

!       Dry deposition
!       ---------------
        call ESMF_StateGet (providerState, trim(prefix)//'SSDP'//trim(field_name), field, __RC__)
        call MAPL_AllocateCoupling (field, __RC__)
        call MAPL_FieldBundleAdd (SS_Bundle_DP, field, __RC__)

!       Wet deposition (Convective scavenging)
!       --------------------------------------
        call ESMF_StateGet (providerState, trim(prefix)//'SSSV'//trim(field_name), field, __RC__)
        call MAPL_AllocateCoupling (field, __RC__)
        call MAPL_FieldBundleAdd (SS_Bundle_DP, field, __RC__)

!       Wet deposition
!       ---------------
        call ESMF_StateGet (providerState, trim(prefix)//'SSWT'//trim(field_name), field, __RC__)
        call MAPL_AllocateCoupling (field, __RC__)
        call MAPL_FieldBundleAdd (SS_Bundle_DP, field, __RC__)

!       Gravitational Settling
!       ----------------------
        call ESMF_StateGet (providerState, trim(prefix)//'SSSD'//trim(field_name), field, __RC__)
        call MAPL_AllocateCoupling (field, __RC__)
        call MAPL_FieldBundleAdd (SS_Bundle_DP, field, __RC__)
    end do


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
    type (ESMF_State)                  :: SS_State

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


    call ESMF_StateGet (EXPORT, trim(COMP_NAME)//'_AERO', SS_State, __RC__)
    call ESMF_AttributeGet (SS_State, name='n_bins', value=n_bins, __RC__)

!   Update interal data pointers with ExtData
!   -----------------------------------------
    do i = 1, n_bins
    write(field_name, '(A, I0.3)') 'ss', i
        call MAPL_GetPointer (INTERNAL, NAME=trim(COMP_NAME)//'::'//trim(field_name), ptr=ptr3d_int, __RC__)
        call MAPL_GetPointer (IMPORT,  NAME='clim'//trim(field_name), ptr=ptr3d_imp, __RC__)

        ptr3d_int = ptr3d_imp
    end do


!  call ESMF_StateGet( IMPORT, 'climSSDP001', field, __RC__ )
!  call ESMF_FieldGet( field, farrayPtr=ptr_test2d, __RC__ )
!  if (mapl_am_i_root()) print*,'SS.data climSSDP001 = ', ptr_test2d

if (mapl_am_i_root()) print*,'SSng ExtData Run_data END'
 


    RETURN_(ESMF_SUCCESS)

  end subroutine Run_data




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

