#include "MAPL_Generic.h"

!=============================================================================
!BOP

! !MODULE: DU2G_GridCompMod - GOCART refactoring of the DU gridded component 

! !INTERFACE:

module DU2G_GridCompMod

!  !USES:
   use ESMF
   use MAPL_Mod
   use Chem_MieMod2G
   use Chem_AeroGeneric

   implicit none
   private

   type(Chem_Mie), dimension(2), save :: gocart2GMieTable
   integer, parameter :: instanceComputational = 1
   integer, parameter :: instanceData          = 2

! !PUBLIC MEMBER FUNCTIONS:
   public  SetServices


! !DESCRIPTION: This module implements GOCART2G's Dust (DU) Gridded Component.

! !REVISION HISTORY:
! 16Oct2019  Sherman, da Silva, Darmenov, Clune -  First attempt at refactoring.

!EOP
!===========================================================================

contains


!BOP

! !IROUTINE: SetServices 

! !INTERFACE:
  subroutine SetServices (GC, RC)

! !ARGUMENTS:
    type (ESMF_GridComp), intent(INOUT)   :: GC  ! gridded component
    integer,              intent(  OUT)   :: RC  ! return code


! !DESCRIPTION: This version uses MAPL_GenericSetServices, which sets
!   the Initialize and Finalize services to generic versions. It also
!   allocates our instance of a generic state and puts it in the 
!   gridded component (GC). Here we only set the two-stage run method and
!   declare the data services.

! !REVISION HISTORY: 
! 16oct2019   E.Sherman  First attempt at refactoring

!EOP

!****************************************************************************
!
!   Locals
    character (len=ESMF_MAXSTR)                 :: COMP_NAME
    type (ESMF_Config)                          :: cfg

    character (len=ESMF_MAXSTR)                 :: field_name
    character (len=ESMF_MAXSTR), allocatable    :: aerosol_names(:)

    integer                                     :: n, i, nCols, aero_n
    real                                        :: DEFVAL
    logical                                     :: data_driven = .true.
    integer, parameter                          :: bins = 5  ! This should equal the number of bins. Is this how we want to handle this?

    __Iam__('SetServices')

!****************************************************************************
!   Begin...

!   Get my name and set-up traceback handle
!   ---------------------------------------
    call ESMF_GridCompGet (GC, NAME=COMP_NAME, __RC__)
    Iam = trim(COMP_NAME) //'::'// Iam

!   Load resource file  
!   -------------------
    cfg = ESMF_ConfigCreate (__RC__)
    call ESMF_ConfigLoadFile (cfg, 'DU2G_GridComp_'//trim(COMP_NAME)//'.rc', RC=STATUS)
    if (STATUS /= 0) then
        if (mapl_am_i_root()) print*,'DU2G_GridComp_'//trim(COMP_NAME)//'.rc does not exist! Loading DU2G_GridComp_DU.rc instead'
        call ESMF_ConfigLoadFile (cfg, 'DU2G_GridComp_DU.rc', __RC__)
    end if


!   Get names of aerosols and write to aerosol_names to add to AERO State
!   ----------------------------------------------------------------------
    call ESMF_ConfigGetDim (cfg, aero_n, nCols, label=('variable_table::'), __RC__ )
    allocate (aerosol_names(aero_n), __STAT__)
    call ESMF_ConfigFindLabel (cfg, 'variable_table::', __RC__ )

    do i = 1, aero_n
        call ESMF_ConfigNextLine( cfg, __RC__ )
        call ESMF_ConfigGetAttribute( cfg, aerosol_names(i), __RC__ )
    end do


!   Set entry points
!   ------------------------
    call MAPL_GridCompSetEntryPoint (GC, ESMF_Method_Initialize,  Initialize, __RC__)
    call MAPL_GridCompSetEntryPoint (GC, ESMF_Method_Run, Run1, __RC__)
    call MAPL_GridCompSetEntryPoint (GC, ESMF_Method_Run, Run2, __RC__)

!   Is DU data driven?
!   ------------------
    call data_driven_ (COMP_NAME, data_driven, __RC__)


!   INTERNAL STATE
!   ---------------
!   Default internal state values
!   -----------------------------
    DEFVAL = 0.0

!   Aerosol Tracers to be transported
!   ---------------------------------
    do i = 1, aero_n
        call MAPL_AddInternalSpec(GC,                                    &
          short_name = trim(COMP_NAME)//'::'//trim(aerosol_names(i)),    &
          long_name  = 'Dust Mixing Ratio (bin '//trim(field_name)//')', &
          units      = 'kg kg-1',                                        &
          restart    = MAPL_RestartOptional,                             &
          default    = DEFVAL,                                           &
          dims       = MAPL_DimsHorzVert,                                &
          vlocation  = MAPL_VLocationCenter, __RC__)
    end do


!   IMPORT STATE
!   -------------
    if (data_driven) then
       do i = 1, bins 
           write (field_name, '(A, I0.3)') '', i
            call MAPL_AddImportSpec(GC,                                     &
              short_name = 'climdu'//trim(field_name),                        &
              long_name  = 'Dust Mixing Ratio (bin '//trim(field_name)//')',  &
              units      = 'kg kg-1',                                         &
              restart    = MAPL_RestartSkip,                                  &
              dims       = MAPL_DimsHorzVert,                                 &
              vlocation  = MAPL_VLocationCenter, __RC__)

!           ! dry deposition
            call MAPL_AddImportSpec(GC,                                       &
              short_name = 'climDUDP'//trim(field_name),                      &
              long_name  = 'Dust Mixing Ratio (bin '//trim(field_name)//')',  &
              units      = 'kg kg-1',                                         &
              dims       = MAPL_DimsHorzOnly,                                 &
              vlocation  = MAPL_VLocationCenter,                              &
              restart    = MAPL_RestartSkip, __RC__)

!           ! wet deposition    
            call MAPL_AddImportSpec(GC,                                       &
               short_name = 'climDUWT'//trim(field_name),                     &
               long_name  = 'Dust Mixing Ratio (bin '//trim(field_name)//')', &
               units      = 'kg kg-1',                                        &
               dims       = MAPL_DimsHorzOnly,                                &
               vlocation  = MAPL_VLocationCenter,                             &
               restart    = MAPL_RestartSkip, __RC__)

!           ! gravitational settling
            call MAPL_AddImportSpec(GC,                                       &
               short_name = 'climDUSD'//trim(field_name),                     &
               long_name  = 'Dust Mixing Ratio (bin '//trim(field_name)//')', &
               units      = 'kg kg-1',                                        &
               dims       = MAPL_DimsHorzOnly,                                &
               vlocation  = MAPL_VLocationCenter,                             &
               restart    = MAPL_RestartSkip, __RC__)

!           convective scavenging
            call MAPL_AddImportSpec(GC,                                       &
               short_name = 'climDUSV'//trim(field_name),                     &
               long_name  = 'Dust Mixing Ratio (bin '//trim(field_name)//')', &
               units      = 'kg kg-1',                                        &
               dims       = MAPL_DimsHorzOnly,                                &
               vlocation  = MAPL_VLocationCenter,                             &
               restart    = MAPL_RestartSkip, __RC__)
        end do
!  else
!       call MAPL_AddImportSpec(GC,                                       &
!          short_name = 'DU_SRC',             &
!          long_name  = 'erod'  ,                                         &
!          units      = '1',                                              &
!          dims       = MAPL_DimsHorzOnly,                                &
!          vlocation  = MAPL_VLocationNone,                               &
!          restart    = MAPL_RestartSkip, __RC__)
    end if ! (data_driven)


!   EXPORT STATE
!   -------------
    if (.not. data_driven) then
#       include "DU2G_ExportSpec___.h"
    end if

!   This state holds fields needed by radiation
!   ---------------------------------------------
    call MAPL_AddExportSpec (GC,                             &
       short_name = trim(COMP_NAME)//'_AERO',                &
       long_name  = 'aerosols_from_'//trim(COMP_NAME),       &
       units      = 'kg kg-1',                               &
       dims       = MAPL_DimsHorzVert,                       &
       vlocation  = MAPL_VLocationCenter,                    &
       datatype   = MAPL_StateItem, __RC__)

!   This state is needed by MOIST - It will contain aerosols
!   ----------------------------------------------------------
    call MAPL_AddExportSpec (GC,                                                  &
       short_name = trim(COMP_NAME)//'_AERO_ACI',                                 &
       long_name  = 'aerosol_cloud_interaction_aerosols_from_'//trim(COMP_NAME),  &
       units      = 'kg kg-1',                                                    &
       dims       = MAPL_DimsHorzVert,                                            &
       vlocation  = MAPL_VLocationCenter,                                         &
       datatype   = MAPL_StateItem, __RC__)

!   This bundle is needed by surface for snow albedo modification
!   by aerosol settling and deposition
!   ~~~DEVELOPERS NOTE~~~ Change to StateItem when possible
!                         This will require refactoring Radiation
!   ---------------------------------------------------------------
    call MAPL_AddExportSpec (GC,                                  &
       short_name = 'AERO_DP',                                    &
       long_name  = 'aerosol_deposition_from_'//trim(COMP_NAME),  &
       units      = 'kg m-2 s-1',                                 &
       dims       = MAPL_DimsHorzOnly,                            &
       datatype   = MAPL_BundleItem, __RC__)


!   Set generic services
!   ----------------------------------
    call MAPL_GenericSetServices (GC, __RC__)


    RETURN_(ESMF_SUCCESS)

  end subroutine SetServices

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!BOP

! !IROUTINE: Initialize 

! !INTERFACE:

  subroutine Initialize (GC, import, export, clock, RC)

! !ARGUMENTS:

    type (ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
    type (ESMF_State),    intent(inout) :: import ! Import state
    type (ESMF_State),    intent(inout) :: export ! Export state
    type (ESMF_Clock),    intent(inout) :: clock  ! The clock
    integer, optional,    intent(  out) :: RC     ! Error code

! !DESCRIPTION: This initializes DU's Grid Component. It primaryily fills 
!               GOCART's AERO states with its dust fields. 

! !REVISION HISTORY: 
! 16oct2019   E.Sherman  First attempt at refactoring

!EOP

!****************************************************************************
! Locals
    character (len=ESMF_MAXSTR)          :: COMP_NAME
    type (MAPL_MetaComp),       pointer  :: MAPL
    type (ESMF_Grid)                     :: grid
    type (ESMF_State)                    :: internal
    type (ESMF_State)                    :: aero, aero_aci
    type (ESMF_State)                    :: providerState
    type (ESMF_Config)                   :: cfg
    type (ESMF_FieldBundle)              :: Bundle_DP

    integer                              :: i, aero_n, nCols
    integer                              :: instance
    type (ESMF_Field)                    :: field, fld
    character (len=ESMF_MAXSTR)          :: field_name, prefix
    character (len=ESMF_MAXSTR), allocatable           :: aerosol_names(:)

    logical                              :: data_driven

    integer                              :: NUM_BANDS
    __Iam__('Initialize')
!****************************************************************************

!   Begin... 

!   Get the target components name and set-up traceback handle.
!   -----------------------------------------------------------
    call ESMF_GridCompGet (GC, grid=grid, name=COMP_NAME, __RC__)
    Iam = trim(COMP_NAME) // trim(Iam)

!   Get my internal MAPL_Generic state
!   -----------------------------------
    call MAPL_GetObjectFromGC (GC, MAPL, __RC__)

!   Load resource file  
!   -------------------
    cfg = ESMF_ConfigCreate (__RC__)
    call ESMF_ConfigLoadFile (cfg, 'DU2G_GridComp_'//trim(COMP_NAME)//'.rc', RC=STATUS)
    if (STATUS /= 0) then
        if (mapl_am_i_root()) print*,'DU2G_GridComp_'//trim(COMP_NAME)//'.rc does not exist! &
                                      loading DU2G_GridComp_DU.rc instead'
        call ESMF_ConfigLoadFile (cfg, 'DU2G_GridComp_DU.rc', __RC__)
    end if

!   Get names of aerosols and write to list to add to AERO State
!   -------------------------------------------------------------
    call ESMF_ConfigGetDim (cfg, aero_n, nCols, label=('variable_table::'), __RC__ )
    allocate (aerosol_names(aero_n), __STAT__)
    call ESMF_ConfigFindLabel (cfg, 'variable_table::', __RC__ )

    do i = 1, aero_n
        call ESMF_ConfigNextLine (cfg, __RC__ )
        call ESMF_ConfigGetAttribute (cfg, aerosol_names(i), __RC__ )
    end do


!   Call Generic Initialize 
!   ------------------------
    call MAPL_GenericInitialize (GC, import, export, clock, __RC__)

!   Get parameters from generic state.
!   -----------------------------------
    call MAPL_Get (MAPL, INTERNAL_ESMF_STATE=internal, __RC__)

!   Is DU data driven?
!   ------------------
    call data_driven_ (COMP_NAME, data_driven, __RC__)


!   If this is a data component, the data is provided in the import
!   state via ExtData instead of the actual GOCART children
!   ---------------------------------------------------------------- 
    if ( data_driven ) then
        providerState = import
        prefix = 'clim'
    else
        providerState = export
        prefix = ''
    end if


!   Fill AERO States with dust fields
!   ------------------------------------
    call ESMF_StateGet (export, trim(COMP_NAME)//'_AERO'    , aero    , __RC__)
    call ESMF_StateGet (export, trim(COMP_NAME)//'_AERO_ACI', aero_aci, __RC__)
    call ESMF_StateGet (export, 'AERO_DP' , Bundle_DP, __RC__)

!   Add list of aerosol names to AERO state 
    call ESMF_AttributeSet(aero, name='aerosol_names', valueList=aerosol_names, itemCount=size(aerosol_names), __RC__)


!   Fill AERO States and Bundle
!   ------------------------------
    do i = 1, size(aerosol_names)
        write (field_name, '(A, I0.3)') '', i
        call ESMF_StateGet (internal, trim(COMP_NAME)//'::'//trim(aerosol_names(i)), field, __RC__)
        fld = MAPL_FieldCreate (field, name='du'//trim(field_name), __RC__)
        call MAPL_StateAdd (aero    , fld, __RC__)
        call MAPL_StateAdd (aero_aci, fld, __RC__)

!       Dry deposition
!       ---------------
        call ESMF_StateGet (providerState, trim(prefix)//'DUDP'//trim(field_name), field, __RC__)
        call MAPL_AllocateCoupling (field, __RC__)
        call MAPL_FieldBundleAdd (Bundle_DP, field, __RC__)

!        Wet deposition (Convective scavenging)
!        --------------------------------------
        call ESMF_StateGet (providerState, trim(prefix)//'DUSV'//trim(field_name), field, __RC__)
        call MAPL_AllocateCoupling (field, __RC__)
        call MAPL_FieldBundleAdd (Bundle_DP, field, __RC__)

!        Wet deposition
!        ---------------
        call ESMF_StateGet (providerState, trim(prefix)//'DUWT'//trim(field_name), field, __RC__)
        call MAPL_AllocateCoupling (field, __RC__)
        call MAPL_FieldBundleAdd (Bundle_DP, field, __RC__)

!        Gravitational Settling
!        ----------------------
        call ESMF_StateGet (providerState, trim(prefix)//'DUSD'//trim(field_name), field, __RC__)
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


    call MAPL_GetResource (MAPL, NUM_BANDS, 'NUM_BANDS:', __RC__)
    gocart2GMieTable(instance) = Chem_MieCreateng (cfg, NUM_BANDS,__RC__)

    ! Mie Table instance/index
    call ESMF_AttributeSet (aero, name='mie_table_instance', value=instance, __RC__)

    ! Add variables to DU instance's AERO state. This is used in aerosol optics calculations
    call addAero (aero, label='air_pressure_for_aerosol_optics',             label2='PLE', grid=grid, typekind=MAPL_R4, __RC__)
    call addAero (aero, label='relative_humidity_for_aerosol_optics',        label2='RH',  grid=grid, typekind=MAPL_R4, __RC__)
    call addAero (aero, label='extinction_in_air_due_to_ambient_aerosol',    label2='EXT', grid=grid, typekind=MAPL_R8, __RC__)
    call addAero (aero, label='single_scattering_albedo_of_ambient_aerosol', label2='SSA', grid=grid, typekind=MAPL_R8, __RC__)
    call addAero (aero, label='asymmetry_parameter_of_ambient_aerosol',      label2='ASY', grid=grid, typekind=MAPL_R8, __RC__)

    call ESMF_AttributeSet(aero, name='band_for_aerosol_optics',             value=0,     __RC__)

    ! attach the aerosol optics method
    call ESMF_MethodAdd (aero, label='aerosol_optics', userRoutine=aerosol_optics, __RC__)


    RETURN_(ESMF_SUCCESS)

  end subroutine Initialize


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!BOP
! !IROUTINE: Run1 

! !INTERFACE:

  subroutine Run1 (GC, import, export, clock, RC)

    ! !ARGUMENTS:

    type (ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
    type (ESMF_State),    intent(inout) :: import ! Import state
    type (ESMF_State),    intent(inout) :: export ! Export state
    type (ESMF_Clock),    intent(inout) :: clock  ! The clock
    integer, optional,    intent(  out) :: RC     ! Error code:

! !DESCRIPTION: Run method for the Dust Grid Component. Determines whether to run
!               data or computational run method.

!EOP

!****************************************************************************
! Locals
    character (len=ESMF_MAXSTR)       :: COMP_NAME
    type (MAPL_MetaComp), pointer     :: MAPL
    type (ESMF_State)                 :: internal

    logical                           :: data_driven

    __Iam__('Run1')

!*****************************************************************************
!   Begin... 

!   Get my name and set-up traceback handle
!   ---------------------------------------
    call ESMF_GridCompGet (GC, NAME=COMP_NAME, __RC__)
    Iam = trim(COMP_NAME) //'::'// Iam

!   Get my internal MAPL_Generic state
!   -----------------------------------
    call MAPL_GetObjectFromGC (GC, MAPL, __RC__)

!   Get parameters from generic state.
!   -----------------------------------
    call MAPL_Get (MAPL, INTERNAL_ESMF_STATE=internal, __RC__)

!   Is DU data driven?
!   ------------------
    call data_driven_ (COMP_NAME, data_driven, __RC__)

!   Update INTERNAL state variables with ExtData
!   ---------------------------------------------
    if (data_driven) then
      call Run_data (GC, import, export, internal, __RC__)
    else
       call Run2 (GC, import, export, clock, __RC__)
    end if


    RETURN_(ESMF_SUCCESS)

  end subroutine Run1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!BOP
! !IROUTINE: Run2 

! !INTERFACE:

  subroutine Run2 (GC, import, export, clock, RC)

    ! !ARGUMENTS:

    type (ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
    type (ESMF_State),    intent(inout) :: import ! Import state
    type (ESMF_State),    intent(inout) :: export ! Export state
    type (ESMF_Clock),    intent(inout) :: clock  ! The clock
    integer, optional,    intent(  out) :: RC     ! Error code:

! !DESCRIPTION: Run2 method for the Dust Grid Component.

!EOP

!****************************************************************************
! Locals
    character (len=ESMF_MAXSTR)       :: COMP_NAME
    type (MAPL_MetaComp), pointer     :: MAPL
    type (ESMF_State)                 :: internal

    logical                           :: data_driven

    __Iam__('Run2')

!*****************************************************************************
!   Begin... 

!   Get my name and set-up traceback handle
!   ---------------------------------------
    call ESMF_GridCompGet (GC, NAME=COMP_NAME, __RC__)
    Iam = trim(COMP_NAME) // '::' // Iam


!   Get my internal MAPL_Generic state
!   -----------------------------------
    call MAPL_GetObjectFromGC (GC, MAPL, __RC__)

!   Get parameters from generic state.
!   -----------------------------------
!    call MAPL_Get (MAPL, INTERNAL_ESMF_STATE=internal, __RC__)



    RETURN_(ESMF_SUCCESS)

  end subroutine Run2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!BOP
! !IROUTINE: Run_data -- ExtData Dust Grid Component

! !INTERFACE:

  subroutine Run_data (GC, import, export, internal, RC)

    ! !ARGUMENTS:

    type (ESMF_GridComp), intent(inout) :: GC       ! Gridded component 
    type (ESMF_State),    intent(inout) :: import   ! Import state
    type (ESMF_State),    intent(inout) :: export   ! Export state
    type (ESMF_State),    intent(inout) :: internal ! Interal state
    integer, optional,    intent(  out) :: RC       ! Error code:

! !DESCRIPTION: Updates pointers in Internal state with fields from ExtData. 

!EOP

!****************************************************************************
! Locals
    character (len=ESMF_MAXSTR)        :: COMP_NAME
    type (MAPL_MetaComp), pointer      :: MAPL
    type (ESMF_State)                  :: aero 

    integer                            :: i, n
    character (len=ESMF_MAXSTR)        :: field_name

    real, pointer, dimension(:,:,:)    :: ptr3d_int
    real, pointer, dimension(:,:,:)    :: ptr3d_imp

    __Iam__('Run_data')


!*****************************************************************************
!   Begin... 

!   Get my name and set-up traceback handle
!   ---------------------------------------
    call ESMF_GridCompGet (GC, NAME=COMP_NAME, __RC__)

    call ESMF_StateGet (export, trim(COMP_NAME)//'_AERO', aero, __RC__)
    call ESMF_AttributeGet (aero, name='aerosol_names', itemCount=n, __RC__)

!   Update interal data pointers with ExtData
!   -----------------------------------------
    do i = 1, n
    write (field_name, '(A, I0.3)') 'du', i
        call MAPL_GetPointer (internal, NAME=trim(COMP_NAME)//'::'//trim(field_name), ptr=ptr3d_int, __RC__)
        call MAPL_GetPointer (import,  NAME='clim'//trim(field_name), ptr=ptr3d_imp, __RC__)

        ptr3d_int = ptr3d_imp
    end do


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
    real(kind=8), dimension(:,:,:), pointer                  :: var
    real, dimension(:,:,:), pointer                  :: q
    real, dimension(:,:,:,:), pointer                :: q_4d

    character (len=ESMF_MAXSTR)                      :: fld_name
    type(ESMF_Field)                                 :: fld
    character (len=ESMF_MAXSTR)                      :: COMP_NAME
    character (len=ESMF_MAXSTR), allocatable         :: aerosol_names(:)

    real(kind = 8), dimension(:,:,:), allocatable              :: ext, ssa, asy  ! (lon:,lat:,lev:)

    integer                                          :: instance
    integer                                          :: n
    integer                                          :: i1, j1, i2, j2, km
    integer                                          :: band, offset
    integer, parameter                               :: n_bands = 1

    real    :: x
    integer :: i, j, k

    __Iam__('DU2G::aerosol_optics')


!   Begin... 

!   Mie Table instance/index
!   ------------------------
    call ESMF_AttributeGet (state, name='mie_table_instance', value=instance, __RC__)

!   Get aerosol names in this state
!   -------------------------------
    call ESMF_AttributeGet (state, name='aerosol_names', itemCount=n, __RC__)
    allocate (aerosol_names(n), __STAT__)
    call ESMF_AttributeGet (state, name='aerosol_names', valueList=aerosol_names, __RC__)


!   Radiation band
!   --------------
    band = 0
    call ESMF_AttributeGet (state, name='band_for_aerosol_optics', value=band, __RC__)
    offset = band - n_bands


!   Pressure at layer edges 
!   ------------------------
    call ESMF_AttributeGet (state, name='air_pressure_for_aerosol_optics', value=fld_name, __RC__)
    call MAPL_GetPointer (state, ple, trim(fld_name), __RC__)

    i1 = lbound(ple, 1); i2 = ubound(ple, 1)
    j1 = lbound(ple, 2); j2 = ubound(ple, 2)
                         km = ubound(ple, 3)

!   Relative humidity
!   -----------------
    call ESMF_AttributeGet (state, name='relative_humidity_for_aerosol_optics', value=fld_name, __RC__)
    call MAPL_GetPointer (state, rh, trim(fld_name), __RC__)


    allocate(ext(i1:i2, j1:j2, km), &
             ssa(i1:i2, j1:j2, km), &
             asy(i1:i2, j1:j2, km), __STAT__)

    allocate(q_4d(i1:i2, j1:j2, km, size(aerosol_names)), __STAT__)

    do n = 1, size(aerosol_names)
        call ESMF_StateGet (state, trim(aerosol_names(n)), field=fld, __RC__)
        call ESMF_FieldGet (fld, farrayPtr=q, __RC__)

        do k = 1, km
            do j = j1, j2
                do i = i1, i2
                    x = ((PLE(i,j,k) - PLE(i,j,k-1))*0.01)*(100./MAPL_GRAV)
                    q_4d(i,j,k,n) = x * q(i,j,k)
                end do
            end do
        end do
    end do


    call mie_ (gocart2GMieTable(instance), aerosol_names, n_bands, offset, q_4d, rh, ext, ssa, asy, __RC__)

    call ESMF_AttributeGet (state, name='extinction_in_air_due_to_ambient_aerosol', value=fld_name, __RC__)
    if (fld_name /= '') then
        call MAPL_GetPointer (state, var, trim(fld_name), __RC__)
        var = ext(:,:,:)
    end if

    call ESMF_AttributeGet (state, name='single_scattering_albedo_of_ambient_aerosol', value=fld_name, __RC__)
    if (fld_name /= '') then
        call MAPL_GetPointer (state, var, trim(fld_name), __RC__)
        var = ssa(:,:,:)
    end if

    call ESMF_AttributeGet (state, name='asymmetry_parameter_of_ambient_aerosol', value=fld_name, __RC__)
    if (fld_name /= '') then
        call MAPL_GetPointer (state, var, trim(fld_name), __RC__)
        var = asy(:,:,:)
    end if


    deallocate(ext, ssa, asy, __STAT__)


    RETURN_(ESMF_SUCCESS)

  contains

    subroutine mie_(mie_table, aerosol_names, nb, offset, q, rh, bext_, bssa_, gasym_, rc)

    implicit none

    type(Chem_Mie),                intent(inout) :: mie_table        ! mie table
    character(len=ESMF_MAXSTR),    intent(in )   :: aerosol_names(:) ! names of aerosols
    integer,                       intent(in )   :: nb               ! number of bands
    integer,                       intent(in )   :: offset           ! bands offset 
    real,                          intent(in )   :: q(:,:,:,:)       ! aerosol mass mixing ratio, kg kg-1
    real,                          intent(in )   :: rh(:,:,:)        ! relative humidity


    real(kind=8), intent(  out) :: bext_ (size(ext,1),size(ext,2),size(ext,3))
    real(kind=8), intent(  out) :: bssa_ (size(ext,1),size(ext,2),size(ext,3))
    real(kind=8), intent(  out) :: gasym_(size(ext,1),size(ext,2),size(ext,3))

    integer,                       intent(  out) :: rc

    ! local
    integer                               :: STATUS
    character (len=ESMF_MAXSTR)           :: Iam = 'DU2G::aerosol_optics::mie_'
    integer                               :: l, na, idx


    real                           :: bext (size(ext,1),size(ext,2),size(ext,3))     ! extinction
    real                           :: bssa (size(ext,1),size(ext,2),size(ext,3))     ! SSA
    real                           :: gasym(size(ext,1),size(ext,2),size(ext,3))     ! asymmetry parameter

     na = size(aerosol_names)

     ASSERT_ (na == size(q,4))

     bext_  = 0.0d0
     bssa_  = 0.0d0
     gasym_ = 0.0d0

     do l = 1, na
         idx = Chem_MieQueryIdx(mie_table, trim(aerosol_names(l)), __RC__)
         call Chem_MieQuery(mie_table, idx, real(offset+1.), q(:,:,:,l), rh, bext, gasym=gasym, ssa=bssa)

         bext_  = bext_  +             bext     ! total extinction due to dust
         bssa_  = bssa_  +       (bssa*bext)    ! total scattering
         gasym_ = gasym_ + gasym*(bssa*bext)    ! sum of (asy * sca)

!         bext_  = bext_  +             bext     ! total extinction due to dust
!         bssa_  = bssa_  +             bssa     ! total scattering
!         gasym_ = gasym_ +             gasym    ! sum of (asy * sca)
     end do

!     bext  = bext_
!     bssa  = bssa_
!     gasym = gasym_

!if (mapl_am_i_root()) print*,'MIE sum(bext) = ',sum(bext)


     RETURN_(ESMF_SUCCESS)

    end subroutine mie_

  end subroutine aerosol_optics



!---------------------------------------------------------------------------------------
  subroutine data_driven_(COMP_NAME, data_driven, RC)

    !ARGUMENTS:
    integer, optional,               intent(  out)   :: RC          ! Error code:
    character (len=ESMF_MAXSTR),     intent(in   )   :: COMP_NAME
    logical,                         intent(  out)   :: data_driven

    !Local
    integer                                          :: i

!   Description: Determines whether Dust is data driven or not.

     __Iam__('data_driven_')

!   Begin... 

!   Is DU data driven?
!   ------------------
    data_driven = .false.

    i = index(COMP_NAME, 'data')
    if (i > 0) then
      data_driven = .true.
    end if

    RETURN_(ESMF_SUCCESS)

  end subroutine data_driven_

!-------------------------------------------------------------------------------------


end module DU2G_GridCompMod

