#include "MAPL_Generic.h"

!=============================================================================
!BOP

! !MODULE: CA2G_GridCompMod - GOCART Carbonaceous Aerosol gridded component 

! !INTERFACE:
module CA2G_GridCompMod

!  !USES:
   use ESMF
   use MAPL
   use Chem_MieTableMod2G
   use Chem_AeroGeneric
   use iso_c_binding, only: c_loc, c_f_pointer, c_ptr

   use Chem_UtilMod
   use GOCART2G_Process       ! GOCART2G process library

   implicit none
   private

   integer, parameter :: instanceComputational = 1
   integer, parameter :: instanceData          = 2

! !PUBLIC MEMBER FUNCTIONS:
   public  SetServices


! !DESCRIPTION: This module implements GOCART2G's Carbonaceous Aerosol (CA) Gridded Component.

! !REVISION HISTORY:
! 16Oct2019  Sherman, da Silva, Darmenov, Clune -  First attempt at refactoring.

!EOP
!===========================================================================

!  !Carbonaceous aerosol state
   type CA2G_GridComp
       type(Chem_Mie), dimension(2)    :: rad_MieTable, diag_MieTable
       real, allocatable      :: fscav(:)       ! scavenging efficiency

       real                   :: CDT            ! chemistry timestep (secs)
       integer                :: km             ! vertical grid dimension
       integer                :: nbins
       integer                :: instance       ! data or computational instance
       integer                :: myDOW = -1     ! my Day of the week: Sun=1, Mon=2,...,Sat=7

   end type CA2G_GridComp 

   type wrap_
      type (CA2G_GridComp), pointer     :: PTR => null()
   end type wrap_

contains

!============================================================================
!BOP

! !IROUTINE: SetServices 

! !INTERFACE:
  subroutine SetServices (GC, RC)

!   !ARGUMENTS:
    type (ESMF_GridComp), intent(INOUT)   :: GC  ! gridded component
    integer,              intent(  OUT)   :: RC  ! return code

!   !DESCRIPTION: This version uses MAPL_GenericSetServices, which sets
!     the Initialize and Finalize services to generic versions. It also
!     allocates our instance of a generic state and puts it in the 
!     gridded component (GC). Here we only set the two-stage run method
!     and declare the data services.

!   !REVISION HISTORY: 
!   16oct2019   E.Sherman, A.Da Silva, A.Darmenov, T.Clune  First attempt at refactoring

!EOP
!============================================================================

!
!   !Locals
    character (len=ESMF_MAXSTR)                 :: COMP_NAME
    type (ESMF_Config)                          :: cfg
    type (wrap_)                                :: wrap
    type (CA2G_GridComp), pointer               :: self
    type (Chem_Mie)                             :: this

    character (len=ESMF_MAXSTR)                 :: field_name

    integer                                     :: n, i, nCols, nbins
    real                                        :: DEFVAL
    logical                                     :: data_driven = .true.

    __Iam__('SetServices')

!****************************************************************************
!   Begin...

!   Get my name and set-up traceback handle
!   ---------------------------------------
    call ESMF_GridCompGet (GC, NAME=COMP_NAME, __RC__)
    Iam = trim(COMP_NAME) //'::'// Iam

if (mapl_am_I_root()) print*,'CA2G SetServices COMP_NAME = ',trim(COMP_NAME)

!   Wrap internal state for storing in GC
!   -------------------------------------
    allocate (self, __STAT__)
    wrap%ptr => self

!   Load resource file  
!   -------------------
    cfg = ESMF_ConfigCreate (__RC__)
    call ESMF_ConfigLoadFile (cfg, 'CA2G_GridComp_'//trim(COMP_NAME)//'.rc', rc=status)
    if (status /= 0) then
        if (mapl_am_i_root()) print*,'CA2G_GridComp_'//trim(COMP_NAME)//'.rc does not exist! &
                                      Loading CA2G_GridComp_CA.rc instead'
        call ESMF_ConfigLoadFile (cfg, 'CA2G_GridComp_CA.rc', __RC__)
    end if

!   Get nbins from cfg
    call ESMF_ConfigGetAttribute (cfg, self%nbins, label='nbins:', __RC__)
    nbins = self%nbins

if (mapl_am_I_root()) print*,'CA2G TEST1'

!   Parse config file into private internal state
!   ----------------------------------------------
    allocate(self%fscav(nbins), __STAT__)

    call ESMF_ConfigGetAttribute (cfg, self%fscav, label='fscav:', __RC__)
    call ESMF_ConfigGetAttribute (cfg, self%myDOW, label='my_day_of_week:', default=-1, __RC__)

!   Is CA data driven?
!   ------------------
    call determine_data_driven (COMP_NAME, data_driven, __RC__)

!   Set entry points
!   ------------------------
    call MAPL_GridCompSetEntryPoint (GC, ESMF_Method_Initialize,  Initialize, __RC__)
    call MAPL_GridCompSetEntryPoint (GC, ESMF_Method_Run, Run, __RC__)
    if (data_driven /= .true.) then
       call MAPL_GridCompSetEntryPoint (GC, ESMF_Method_Run, Run2, __RC__)
    end if

!   INTERNAL STATE
!   ---------------
!   Default internal state values
!   -----------------------------
    DEFVAL = 0.0

#include "CA2G_Internal___.h"

if (mapl_am_I_root()) print*,'CA2G TEST2'

!   IMPORT STATE
!   -------------
    if (data_driven) then

!      Pressure at layer edges
!      -----------------------
       call MAPL_AddImportSpec(GC,                            &
          SHORT_NAME = 'PLE',                                 &
          LONG_NAME  = 'air_pressure',                        &
          UNITS      = 'Pa',                                  &
          DIMS       = MAPL_DimsHorzVert,                     &
          VLOCATION  = MAPL_VLocationEdge,                    &
          RESTART    = MAPL_RestartSkip,     __RC__)

!      RH: is between 0 and 1
!      ----------------------
       call MAPL_AddImportSpec(GC,                            &
          SHORT_NAME = 'RH2',                                 &
          LONG_NAME  = 'Rel_Hum_after_moist',                 &
          UNITS      = '1',                                   &
          DIMS       = MAPL_DimsHorzVert,                     &
          VLOCATION  = MAPL_VLocationCenter,                  &
          RESTART    = MAPL_RestartSkip,     __RC__)

! UPDATE  OCphilic and OCphobic METADATA!!!!!
          call MAPL_AddImportSpec(GC,                            &
             short_name = 'climOCphobic',                        &
             long_name  = 'Organic Carbon phobic Mixing Ratio',  &
             units      = 'kg kg-1',                             &
             restart    = MAPL_RestartSkip,                      &
             dims       = MAPL_DimsHorzVert,                     &
             vlocation  = MAPL_VLocationCenter, __RC__)

          call MAPL_AddImportSpec(GC,                            &
             short_name = 'climOCphilic',                        &
             long_name  = 'Organic Carbon philic Mixing Ratio',  &
             units      = 'kg kg-1',                             &
             restart    = MAPL_RestartSkip,                      &
             dims       = MAPL_DimsHorzVert,                     &
             vlocation  = MAPL_VLocationCenter, __RC__)

       do i = 1, nbins
          write (field_name, '(A, I0.3)') '', i

!        !dry deposition
          call MAPL_AddImportSpec(GC,                                       &
             short_name = 'climOCDP'//trim(field_name),                      &
             long_name  = 'Organic Carbon Mixing Ratio (bin '//trim(field_name)//')',  &
             units      = 'kg kg-1',                                         &
             dims       = MAPL_DimsHorzOnly,                                 &
             vlocation  = MAPL_VLocationCenter,                              &
             restart    = MAPL_RestartSkip, __RC__)

!        !wet deposition    
          call MAPL_AddImportSpec(GC,                                       &
             short_name = 'climOCWT'//trim(field_name),                     &
             long_name  = 'Organic Carbon Mixing Ratio (bin '//trim(field_name)//')', &
             units      = 'kg kg-1',                                        &
             dims       = MAPL_DimsHorzOnly,                                &
             vlocation  = MAPL_VLocationCenter,                             &
             restart    = MAPL_RestartSkip, __RC__)

!        !gravitational settling
          call MAPL_AddImportSpec(GC,                                       &
             short_name = 'climOCSD'//trim(field_name),                     &
             long_name  = 'Organic Carbon Mixing Ratio (bin '//trim(field_name)//')', &
             units      = 'kg kg-1',                                        &
             dims       = MAPL_DimsHorzOnly,                                &
             vlocation  = MAPL_VLocationCenter,                             &
             restart    = MAPL_RestartSkip, __RC__)

!        !convective scavenging
          call MAPL_AddImportSpec(GC,                                       &
             short_name = 'climOCSV'//trim(field_name),                     &
             long_name  = 'Organic Carbon Mixing Ratio (bin '//trim(field_name)//')', &
             units      = 'kg kg-1',                                        &
             dims       = MAPL_DimsHorzOnly,                                &
             vlocation  = MAPL_VLocationCenter,                             &
             restart    = MAPL_RestartSkip, __RC__)
        end do
    end if ! (data_driven)


!   Computational Import and Export states
!   --------------------------------------
    if (.not. data_driven) then
#include "CA2G_Export___.h"
#include "CA2G_Import___.h"
    end if

if (mapl_am_I_root()) print*,'CA2G TEST3'

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
       short_name = trim(COMP_NAME)//'_AERO_DP',                                    &
       long_name  = 'aerosol_deposition_from_'//trim(COMP_NAME),  &
       units      = 'kg m-2 s-1',                                 &
       dims       = MAPL_DimsHorzOnly,                            &
       datatype   = MAPL_BundleItem, __RC__)


!   Store internal state in GC
!   --------------------------
    call ESMF_UserCompSetInternalState ( GC, 'CA2G_GridComp', wrap, STATUS )
    VERIFY_(STATUS)

!   Set generic services
!   ----------------------------------
    call MAPL_GenericSetServices (GC, __RC__)

    RETURN_(ESMF_SUCCESS)

if (mapl_am_I_root()) print*,'CA2G SetServices END'

  end subroutine SetServices


!============================================================================
!BOP

! !IROUTINE: Initialize 

! !INTERFACE:
  subroutine Initialize (GC, import, export, clock, RC)

!   !ARGUMENTS:
    type (ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
    type (ESMF_State),    intent(inout) :: import ! Import state
    type (ESMF_State),    intent(inout) :: export ! Export state
    type (ESMF_Clock),    intent(inout) :: clock  ! The clock
    integer, optional,    intent(  out) :: RC     ! Error code

! !DESCRIPTION: This initializes CA's Grid Component. It primaryily fills 
!               GOCART's AERO states with its carbonaceous aerosol fields. 

! !REVISION HISTORY: 
! 16oct2019   E.Sherman  First attempt at refactoring

!EOP
!============================================================================
!   !Locals
    character (len=ESMF_MAXSTR)          :: COMP_NAME
    type (MAPL_MetaComp),       pointer  :: MAPL
    type (ESMF_Grid)                     :: grid
    type (ESMF_State)                    :: internal
    type (ESMF_State)                    :: aero, aero_aci
    type (ESMF_State)                    :: providerState
    type (ESMF_Config)                   :: cfg
    type (ESMF_FieldBundle)              :: Bundle_DP
    type (wrap_)                         :: wrap
    type (CA2G_GridComp), pointer        :: self

    integer, allocatable                 :: mieTable_pointer(:)
    integer                              :: i, nbins, nCols, dims(3), km
    integer                              :: instance
    type (ESMF_Field)                    :: field, fld
    character (len=ESMF_MAXSTR)          :: field_name, prefix
    integer                              :: varCount
    character (len=ESMF_MAXSTR), allocatable   :: varList(:)

    real                                 :: CDT         ! chemistry timestep (secs)
    integer                              :: HDT         ! model     timestep (secs)

    logical                              :: data_driven
    integer                              :: NUM_BANDS

    __Iam__('Initialize')

!****************************************************************************

!   Begin... 

if (mapl_am_I_root()) print*,'CA2G Init BEGIN'

!   Get the target components name and set-up traceback handle.
!   -----------------------------------------------------------
    call ESMF_GridCompGet (GC, grid=grid, name=COMP_NAME, __RC__)
    Iam = trim(COMP_NAME) // '::' //trim(Iam)

!   Get my internal MAPL_Generic state
!   -----------------------------------
    call MAPL_GetObjectFromGC (GC, MAPL, __RC__)

!   Get my internal private state
!   -----------------------------
    call ESMF_UserCompGetInternalState(GC, 'CA2G_GridComp', wrap, STATUS)
    VERIFY_(STATUS)
    self => wrap%ptr

!   Get dimensions
!   ---------------
    call MAPL_GridGet (grid, globalCellCountPerDim=dims, __RC__ )
    km = dims(3)
    self%km = km

!   Get DTs
!   -------
    call MAPL_GetResource(mapl, HDT, Label='RUN_DT:', __RC__)
    call MAPL_GetResource(mapl, CDT, Label='GOCART_DT:', default=real(HDT), __RC__)
    self%CDT = CDT

!   Load resource file  
!   -------------------
    cfg = ESMF_ConfigCreate (__RC__)
    call ESMF_ConfigLoadFile (cfg, 'CA2G_GridComp_'//trim(COMP_NAME)//'.rc', rc=status)
    if (status /= 0) then
        if (mapl_am_i_root()) print*,'CA2G_GridComp_'//trim(COMP_NAME)//'.rc does not exist! &
                                      Loading CA2G_GridComp_CA.rc instead'
        call ESMF_ConfigLoadFile (cfg, 'CA2G_GridComp_CA.rc', __RC__)
    end if

!   Call Generic Initialize 
!   ------------------------
    call MAPL_GenericInitialize (GC, import, export, clock, __RC__)

!   Get parameters from generic state.
!   -----------------------------------
    call MAPL_Get (MAPL, INTERNAL_ESMF_STATE=internal, __RC__)

!   Is DU data driven?
!   ------------------
    call determine_data_driven (COMP_NAME, data_driven, __RC__)

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

!   Add attribute information to internal state varaibles
!   -----------------------------------------------------
    call ESMF_StateGet (internal, itemNameList=varList, __RC__)

!   Fill AERO States with dust fields
!   ------------------------------------
    call ESMF_StateGet (export, trim(COMP_NAME)//'_AERO'    , aero    , __RC__)
    call ESMF_StateGet (export, trim(COMP_NAME)//'_AERO_ACI', aero_aci, __RC__)
    call ESMF_StateGet (export, trim(COMP_NAME)//'_AERO_DP' , Bundle_DP, __RC__)

    call ESMF_StateGet (internal, 'OCphilic', field, __RC__)
    fld = MAPL_FieldCreate (field, 'OCphilic', __RC__)
    call MAPL_StateAdd (aero, fld, __RC__)
    call MAPL_StateAdd (aero_aci, fld, __RC__)

    ! ADD OTHER ATTRIBUTE, HENTRY COEFFICIENTS?
    call ESMF_AttributeSet(field, NAME='ScavengingFractionPerKm', VALUE=self%fscav(1), __RC__)

    call ESMF_StateGet (internal, 'OCphobic', field, __RC__)
    fld = MAPL_FieldCreate (field, 'OCpoblic', __RC__)
    call MAPL_StateAdd (aero, fld, __RC__)
    call MAPL_StateAdd (aero_aci, fld, __RC__)

    ! ADD OTHER ATTRIBUTE, HENTRY COEFFICIENTS?
    call ESMF_AttributeSet(field, NAME='ScavengingFractionPerKm', VALUE=self%fscav(2), __RC__)

!   Dry deposition
!   ---------------
    call append_to_bundle('OCDP', providerState, prefix, Bundle_DP, __RC__)

!   Wet deposition (Convective scavenging)
!   --------------------------------------
    call append_to_bundle('OCSV', providerState, prefix, Bundle_DP, __RC__)

!   Wet deposition
!   ---------------
    call append_to_bundle('OCWT', providerState, prefix, Bundle_DP, __RC__)

!   Gravitational Settling
!   ----------------------
    call append_to_bundle('OCSD', providerState, prefix, Bundle_DP, __RC__)

!   Set AERO States' attributes
!   ----------------------------
    if (data_driven) then
       instance = instanceData
    else
       instance = instanceComputational
    end if

    self%instance = instance

!   Create Radiation Mie Table
!   --------------------------
    call MAPL_GetResource (MAPL, NUM_BANDS, 'NUM_BANDS:', __RC__)

!   Get file names for the optical tables
    call ESMF_ConfigGetAttribute (cfg, self%rad_MieTable(instance)%optics_file, &
                                  label="aerosol_radBands_optics_file:", __RC__ )

    allocate (self%rad_MieTable(instance)%channels(NUM_BANDS), __STAT__ )

    call ESMF_ConfigGetAttribute (cfg, self%rad_MieTable(instance)%channels, label= "BANDS:", &
                                 count=self%rad_MieTable(instance)%nch, rc=status)

    if (rc /= 0) then
       do i = 1, NUM_BANDS
          self%rad_MieTable(instance)%channels(i) = i
       end do
    end if

    allocate (self%rad_MieTable(instance)%mie_aerosol, __STAT__)
    self%rad_MieTable(instance)%mie_aerosol = Chem_MieTableCreate (self%rad_MieTable(instance)%optics_file, rc)
    call Chem_MieTableRead (self%rad_MieTable(instance)%mie_aerosol, NUM_BANDS, self%rad_MieTable(instance)%channels, rc)

!   Create Diagnostics Mie Table
!   -----------------------------
!   Get file names for the optical tables
    call ESMF_ConfigGetAttribute (cfg, self%diag_MieTable(instance)%optics_file, &
                                  label="aerosol_monochromatic_optics_file:", __RC__ )
    call ESMF_ConfigGetAttribute (cfg, self%diag_MieTable(instance)%nch, label="n_channels:", __RC__)
    call ESMF_ConfigGetAttribute (cfg, self%diag_MieTable(instance)%nmom, label="n_moments:", default=0,  __RC__)
    allocate (self%diag_MieTable(instance)%channels(self%diag_MieTable(instance)%nch), __STAT__ )
    call ESMF_ConfigGetAttribute (cfg, self%diag_MieTable(instance)%channels, label= "r_channels:", __RC__)

    allocate (self%diag_MieTable(instance)%mie_aerosol, __STAT__)
    self%diag_MieTable(instance)%mie_aerosol = Chem_MieTableCreate (self%diag_MieTable(instance)%optics_file, __RC__ )
    call Chem_MieTableRead (self%diag_MieTable(instance)%mie_aerosol, self%diag_MieTable(instance)%nch, &
                            self%diag_MieTable(instance)%channels, rc, nmom=self%diag_MieTable(instance)%nmom)

    ! Mie Table instance/index
    call ESMF_AttributeSet(aero, name='mie_table_instance', value=instance, __RC__)

    ! Add variables to SS instance's aero state. This is used in aerosol optics calculations
    call add_aero (aero, label='air_pressure_for_aerosol_optics',             label2='PLE', grid=grid, typekind=MAPL_R4, __RC__)
    call add_aero (aero, label='relative_humidity_for_aerosol_optics',        label2='RH',  grid=grid, typekind=MAPL_R4,__RC__)
!   call ESMF_StateGet (import, 'PLE', field, __RC__)
!   call MAPL_StateAdd (aero, field, __RC__)
!   call ESMF_StateGet (import, 'RH2', field, __RC__)
!   call MAPL_StateAdd (aero, field, __RC__)

    call add_aero (aero, label='extinction_in_air_due_to_ambient_aerosol',    label2='EXT', grid=grid, typekind=MAPL_R8,__RC__)
    call add_aero (aero, label='single_scattering_albedo_of_ambient_aerosol', label2='SSA', grid=grid, typekind=MAPL_R8,__RC__)
    call add_aero (aero, label='asymmetry_parameter_of_ambient_aerosol',      label2='ASY', grid=grid, typekind=MAPL_R8,__RC__)

    call ESMF_AttributeSet(aero, name='band_for_aerosol_optics',             value=0,     __RC__)

    mieTable_pointer = transfer(c_loc(self), [1])
    call ESMF_AttributeSet(aero, name='mieTable_pointer', valueList=mieTable_pointer, itemCount=size(mieTable_pointer), __RC__)

! ADD FOR OCphilic
    call ESMF_AttributeSet(aero, name='internal_varaible_name', value='OCphobic', __RC__)

!    call ESMF_MethodAdd(AERO, label='aerosol_optics', userRoutine=aerosol_optics, __RC__)

    RETURN_(ESMF_SUCCESS)

if (mapl_am_I_root()) print*,'CA2G Init END'

  end subroutine Initialize

!============================================================================

!BOP
! !IROUTINE: Run 

! !INTERFACE:
  subroutine Run (GC, import, export, clock, rc)

!   !ARGUMENTS:
    type (ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
    type (ESMF_State),    intent(inout) :: import ! Import state
    type (ESMF_State),    intent(inout) :: export ! Export state
    type (ESMF_Clock),    intent(inout) :: clock  ! The clock
    integer, optional,    intent(  out) :: rc     ! Error code:

! !DESCRIPTION: Run method for the Sea Salt Grid Component. Determines whether to run
!               data or computational run method.

!EOP
!============================================================================
!   !Locals
    character (len=ESMF_MAXSTR)       :: COMP_NAME
    type (MAPL_MetaComp), pointer     :: MAPL
    type (ESMF_State)                 :: internal

    logical                           :: data_driven

    __Iam__('Run')

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

!   Is SS data driven?
!   ------------------
    call determine_data_driven (COMP_NAME, data_driven, __RC__)

!   Update INTERNAL state variables with ExtData
!   ---------------------------------------------
    if (data_driven) then
       call Run_data (GC, import, export, internal, __RC__)
    else
       call Run1 (GC, import, export, clock, __RC__)
    end if

    RETURN_(ESMF_SUCCESS)

  end subroutine Run


!============================================================================
!BOP
! !IROUTINE: Run1 

! !INTERFACE:
  subroutine Run1 (GC, import, export, clock, RC)

!   !ARGUMENTS:
    type (ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
    type (ESMF_State),    intent(inout) :: import ! Import state
    type (ESMF_State),    intent(inout) :: export ! Export state
    type (ESMF_Clock),    intent(inout) :: clock  ! The clock
    integer, optional,    intent(  out) :: RC     ! Error code:

! !DESCRIPTION:  Computes emissions/sources for Sea Salt

!EOP
!============================================================================
!   !Locals
    character (len=ESMF_MAXSTR)       :: COMP_NAME
    type (MAPL_MetaComp), pointer     :: mapl
    type (ESMF_State)                 :: internal
    type (ESMF_Grid)                  :: grid
    type (wrap_)                      :: wrap
    type (CA2G_GridComp), pointer     :: self
    type(ESMF_Time)                   :: time

    character(len=ESMF_MAXSTR), dimension(:), allocatable :: itemNameList
    integer          :: idow
    character(len=3) :: cdow
    integer         :: nymd, nhms, iyr, imm, idd, ihr, imn, isc


#include "CA2G_DeclarePointer___.h"

   __Iam__('Run1')

!*****************************************************************************
!   Begin... 

if(mapl_am_i_root()) print*,'CA2G Run1 BEGIN'

!   Get my name and set-up traceback handle
!   ---------------------------------------
    call ESMF_GridCompGet (GC, grid=grid, NAME=COMP_NAME, __RC__)
    Iam = trim(COMP_NAME) //'::'// Iam

!   Get my internal MAPL_Generic state
!   -----------------------------------
    call MAPL_GetObjectFromGC (GC, mapl, __RC__)

!   Get parameters from generic state.
!   -----------------------------------
    call MAPL_Get (mapl, INTERNAL_ESMF_STATE=internal, __RC__)

#include "CA2G_GetPointer___.h"

!   Get my private internal state
!   ------------------------------
    call ESMF_UserCompGetInternalState(GC, 'CA2G_GridComp', wrap, STATUS)
    VERIFY_(STATUS)
    self => wrap%ptr

    allocate(itemNameList(self%nbins), __STAT__)
    call ESMF_StateGet(internal, itemNameList=itemNameList, __RC__)
if(mapl_am_i_root()) print*,'CA2G itemNameList = ',itemNameList

!   Extract nymd(yyyymmdd) from clock
!   ---------------------------------
    call ESMF_ClockGet (clock, currTime=time, __RC__)
    call ESMF_TimeGet (time ,YY=iyr, MM=imm, DD=idd, H=ihr, M=imn, S=isc, __RC__)
    call MAPL_PackTime (nymd, iyr, imm , idd)
    call MAPL_PackTime (nhms, ihr, imn, isc)

!   Reset tracer to zero at 0Z on specific day of week
!   --------------------------------------------------
!#if 0
    idow = Chem_UtilIdow(nymd)
    if ( (nhms==0) .and. (idow == self%myDOW) ) then
       cdow = Chem_UtilCdow(nymd)
       OCphobic = tiny(1.) ! avoid division by zero
       OCphilic = tiny(1.) ! avoid division by zero
       if ( MAPL_AM_I_ROOT() ) then
          print *, '<> OC '//cdow//' tracer being set to zero on ', nymd, nhms
       end if
    end if
!#endif



    RETURN_(ESMF_SUCCESS)

if(mapl_am_i_root()) print*,'CA2G Run1 END'

  end subroutine Run1


!============================================================================
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

! !DESCRIPTION: Run2 method for the CA Grid Component.

!EOP
!============================================================================
! Locals
    character (len=ESMF_MAXSTR)       :: COMP_NAME
    type (MAPL_MetaComp), pointer     :: MAPL
    type (ESMF_State)                 :: internal
    type (wrap_)                      :: wrap
    type (CA2G_GridComp), pointer     :: self

    integer                           :: import_shape(2), i2, j2
    integer                           :: n
    real, allocatable, dimension(:,:) :: drydepositionfrequency, dqa
    real                              :: fwet
    logical                           :: KIN

    real, parameter ::  cpd    = 1004.16


#include "CA2G_DeclarePointer___.h"

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
    call MAPL_Get (MAPL, INTERNAL_ESMF_STATE=internal, __RC__)

#include "CA2G_GetPointer___.h"

!   Get my private internal state
!   ------------------------------
    call ESMF_UserCompGetInternalState(GC, 'CA2G_GridComp', wrap, STATUS)
    VERIFY_(STATUS)
    self => wrap%ptr





    RETURN_(ESMF_SUCCESS)

  end subroutine Run2

!============================================================================
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
!============================================================================
! Locals
    character (len=ESMF_MAXSTR)        :: COMP_NAME
    type (MAPL_MetaComp), pointer      :: MAPL
    type (ESMF_State)                  :: aero
    type (wrap_)                      :: wrap
    type (CA2G_GridComp), pointer     :: self

    integer                            :: i, n
    character (len=ESMF_MAXSTR)        :: field_name

    real, pointer, dimension(:,:,:)    :: ptr3d_int
    real, pointer, dimension(:,:,:)    :: ptr3d_imp

    __Iam__('Run_data')

!*****************************************************************************
! Begin... 

! Get my name and set-up traceback handle
! ---------------------------------------
    call ESMF_GridCompGet (GC, NAME=COMP_NAME, __RC__)
    Iam = trim(COMP_NAME) //'::'//Iam

if (mapl_am_I_root()) print*,trim(comp_name),' Run_data BEGIN'

!   Update interal data pointers with ExtData
!   -----------------------------------------
    call MAPL_GetPointer (internal, NAME='OCphobic', ptr=ptr3d_int, __RC__)
    call MAPL_GetPointer (import, NAME='climOCphobic', ptr=ptr3d_imp, __RC__)
    ptr3d_int = ptr3d_imp

    call MAPL_GetPointer (internal, NAME='Cphilic', ptr=ptr3d_int, __RC__)
    call MAPL_GetPointer (import, NAME='climOCphilic', ptr=ptr3d_imp, __RC__)
    ptr3d_int = ptr3d_imp

if (mapl_am_I_root()) print*,trim(comp_name),' Run_data END'

    RETURN_(ESMF_SUCCESS)

  end subroutine Run_data



!---------------------------------------------------------------------------------------

end module CA2G_GridCompMod


