#include "MAPL_Generic.h"

!=============================================================================
!BOP

! !MODULE: NI2G_GridCompMod - GOCART Nitrate gridded component 

! !INTERFACE:
module NI2G_GridCompMod

! !USES:
   use ESMF
   use MAPL
   use Chem_MieTableMod2G
   use Chem_AeroGeneric
   use iso_c_binding, only: c_loc, c_f_pointer, c_ptr

   use Chem_UtilMod
   use GOCART2G_Process       ! GOCART2G process library
   use GA_GridCompMod

   implicit none
   private

   integer, parameter :: instanceComputational = 1
   integer, parameter :: instanceData          = 2
   real, parameter :: OCEAN=0.0, LAND = 1.0, SEA_ICE = 2.0
   integer, parameter     :: DP=kind(1.0d0)

! !PUBLIC MEMBER FUNCTIONS:
   PUBLIC  SetServices

real, parameter ::  chemgrav   = 9.80616

! !DESCRIPTION: This module implements GOCART's Nitrate (NI) Gridded Component.

! !REVISION HISTORY:
! 01July2020  Sherman, da Silva, Darmenov, Clune -  First attempt at refactoring.

!EOP
!===========================================================================

!  !Nitrate state
   type, extends(GA_GridComp) :: NI2G_GridComp
       logical           :: first
       logical           :: recycle_HNO3 = .false.
       real, allocatable :: xhno3(:,:,:)   ! buffer for NITRATE_HNO3 [kg/(m^2 sec)]
   end type NI2G_GridComp

   type wrap_
      type (NI2G_GridComp), pointer     :: PTR => null()
   end type wrap_

contains

!============================================================================
!BOP

! !IROUTINE: SetServices 

! !INTERFACE:
  subroutine SetServices ( GC, RC )

!   !ARGUMENTS:
    type (ESMF_GridComp), intent(INOUT)   :: GC  ! gridded component
    integer,              intent(  OUT)   :: RC  ! return code

!    DESCRIPTION: This version uses MAPL_GenericSetServices, which sets
!     the Initialize and Finalize services to generic versions. It also
!     allocates our instance of a generic state and puts it in the 
!     gridded component (GC). Here we only set the two-stage run method
!     and declare the data services.

!   !REVISION HISTORY:
!   24oct2019   E.Sherman, A.Da Silva, A.Darmenov, T.Clune  First attempt at refactoring

!EOP
!============================================================================
!

!   !Locals
    character (len=ESMF_MAXSTR)                 :: COMP_NAME
    type (ESMF_Config)                          :: cfg
    type (wrap_)                                :: wrap
    type (NI2G_GridComp), pointer               :: self
    type (Chem_Mie)                             :: this

    character (len=ESMF_MAXSTR)                 :: field_name
    character (len=ESMF_MAXSTR), allocatable    :: aerosol_names(:)

    integer                                     :: n, i, nCols, nbins
    real                                        :: DEFVAL
    logical                                     :: data_driven=.true.

    !development testing variables - to be deleted
    real, dimension(:,:), pointer       :: ptr_test

    __Iam__('SetServices')

!****************************************************************************
!   Begin...

!   Get my name and set-up traceback handle
!   ---------------------------------------
    call ESMF_GridCompGet (GC, NAME=COMP_NAME, __RC__)
    Iam = trim(COMP_NAME) // '::' // Iam

if(mapl_am_i_root()) print*,trim(comp_name),' SetServices BEGIN'

!   Wrap internal state for storing in GC
!   -------------------------------------
    allocate (self, __STAT__)
    wrap%ptr => self

!   Load resource file 
!   -------------------
    cfg = ESMF_ConfigCreate (__RC__)
    call ESMF_ConfigLoadFile (cfg, 'NI2G_GridComp_'//trim(COMP_NAME)//'.rc', rc=status)
    if (status /= 0) then
       if (mapl_am_i_root()) print*,'NI2G_GridComp_'//trim(COMP_NAME)//'.rc does not exist! Loading NI2G_GridComp_NI.rc instead'
       call ESMF_ConfigLoadFile (cfg, 'NI2G_GridComp_NI.rc', __RC__)
    end if

    ! process generic config items
    call self%GA_GridComp%load_from_config( cfg, __RC__)

!   Is NI data driven?
!   ------------------
    call determine_data_driven (COMP_NAME, data_driven, __RC__)

!   Set entry points
!   ------------------------
    call MAPL_GridCompSetEntryPoint (GC, ESMF_METHOD_INITIALIZE,  Initialize, __RC__)
    call MAPL_GridCompSetEntryPoint (GC, ESMF_METHOD_RUN, Run, __RC__)
    if (data_driven /= .true.) then
       call MAPL_GridCompSetEntryPoint (GC, ESMF_Method_Run, Run2, __RC__)
    end if

    DEFVAL = 0.0

!   Import and Internal states if data instance 
!   -------------------------------------------

!    if (data_driven) then

!    do mapl_add*spec

!    end if ! (data_driven)


!   Import, Export, Internal states for computational instance 
!   ----------------------------------------------------------
    if (.not. data_driven) then
#include "NI2G_Export___.h"
#include "NI2G_Import___.h"
#include "NI2G_Internal___.h"
    end if


!   This state holds fields needed by radiation
!   ---------------------------------------------
    call MAPL_AddExportSpec(GC,                                 &
      SHORT_NAME = trim(COMP_NAME)//'_AERO',                   &
       LONG_NAME  = 'aerosols_from_'//trim(COMP_NAME),  &
       UNITS      = 'kg kg-1',                                  &
       DIMS       = MAPL_DimsHorzVert,                          &
       VLOCATION  = MAPL_VLocationCenter,                       &
       DATATYPE   = MAPL_StateItem, __RC__)

!   This state is needed by MOIST - It will contain aerosols
!   ----------------------------------------------------------
!    call MAPL_AddExportSpec(GC,                                                  &
!       SHORT_NAME = trim(COMP_NAME)//'_AERO_ACI',                                &
!       LONG_NAME  = 'aerosol_cloud_interaction_aerosols_from_'//trim(COMP_NAME),  &
!       UNITS      = 'kg kg-1',                                                   &
!       DIMS       = MAPL_DimsHorzVert,                                           &
!       VLOCATION  = MAPL_VLocationCenter,                                        &
!       DATATYPE   = MAPL_StateItem, __RC__)

!   This bundle is needed by surface for snow albedo modification
!   by aerosol settling and deposition
!   DEVELOPMENT NOTE - Change to StateItem in future
!   ---------------------------------------------------------------
!    call MAPL_AddExportSpec(GC,                                   &
!       SHORT_NAME = trim(COMP_NAME)//'_AERO_DP',                  &
!       LONG_NAME  = 'aerosol_deposition_from_'//trim(COMP_NAME),  &
!       UNITS      = 'kg m-2 s-1',                                 &
!       DIMS       = MAPL_DimsHorzOnly,                            &
!       DATATYPE   = MAPL_BundleItem, __RC__)


!   Store internal state in GC
!   --------------------------
    call ESMF_UserCompSetInternalState ( GC, 'NI2G_GridComp', wrap, STATUS )
    VERIFY_(STATUS)

!   Set generic services
!   ----------------------------------
    call MAPL_GenericSetServices (GC, __RC__)


if(mapl_am_i_root()) print*,trim(comp_name),' SetServices END'

    RETURN_(ESMF_SUCCESS)

  end subroutine SetServices



!============================================================================
!BOP

! !IROUTINE: Initialize 

! !INTERFACE:
  subroutine Initialize (GC, IMPORT, EXPORT, CLOCK, RC)

!   !ARGUMENTS:
    type (ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
    type (ESMF_State),    intent(inout) :: IMPORT ! Import state
    type (ESMF_State),    intent(inout) :: EXPORT ! Export state
    type (ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
    integer, optional,    intent(  out) :: RC     ! Error code

! !DESCRIPTION: This initializes the Nitrate gridded component. It primaryily 
!               fills GOCART's AERO states with its nitrate fields. 

! !REVISION HISTORY: 
! 30June2020   E.Sherman  First attempt at refactoring

!EOP
!============================================================================
!   !Locals 
    character (len=ESMF_MAXSTR)          :: COMP_NAME
    type (MAPL_MetaComp),      pointer   :: MAPL
    type (ESMF_Grid)                     :: grid
    type (ESMF_State)                    :: internal
    type (ESMF_State)                    :: aero, aero_aci
    type (ESMF_State)                    :: providerState
    type (ESMF_Config)                   :: cfg
    type (ESMF_FieldBundle)              :: Bundle_DP
    type (wrap_)                         :: wrap
    type (NI2G_GridComp), pointer        :: self

    integer, allocatable                 :: mieTable_pointer(:)
    integer                              :: i, j, nbins, nCols, dims(3), km
    integer                              :: instance
    type (ESMF_Field)                    :: field, fld
    character (len=ESMF_MAXSTR)          :: field_name, prefix, ind
    real                                 :: CDT         ! chemistry timestep (secs)
    integer                              :: HDT         ! model     timestep (secs)

    logical                              :: data_driven
    integer                              :: NUM_BANDS

    type(ESMF_Calendar)     :: calendar
    type(ESMF_Time)         :: currentTime
    type(ESMF_Alarm)        :: alarm_HNO3
    type(ESMF_Time)         :: ringTime
    type(ESMF_TimeInterval) :: ringInterval
    integer                 :: year, month, day, hh, mm, ss

real, pointer :: ssptr(:,:,:,:)

    __Iam__('Initialize')

!****************************************************************************
!   Begin... 

!   Get the target components name and set-up traceback handle.
!   -----------------------------------------------------------
    call ESMF_GridCompGet (GC, grid=grid, name=COMP_NAME, __RC__)
    Iam = trim(COMP_NAME) // '::' //trim(Iam)

if(mapl_am_i_root()) print*,trim(comp_name),' Init BEGIN'

!   Get my internal MAPL_Generic state
!   ----------------------------------- 
    call MAPL_GetObjectFromGC (GC, MAPL, __RC__)
    
!   Get my internal private state  
!   -----------------------------
    call ESMF_UserCompGetInternalState(GC, 'NI2G_GridComp', wrap, STATUS)
    VERIFY_(STATUS)
    self => wrap%ptr

!   Get dimensions
!   ---------------
    call MAPL_GridGet (grid, globalCellCountPerDim=dims, __RC__ )
    km = dims(3)
    self%km = km

    allocate(self%xhno3(dims(1),dims(2),dims(3)), __STAT__)
if(mapl_am_i_root())print*,'NI dims(1) = ',dims(1)
if(mapl_am_i_root())print*,'NI dims(2) = ',dims(2)
if(mapl_am_i_root())print*,'NI dims(3) = ',dims(3)

!   Get DTs
!   -------
    call MAPL_GetResource(mapl, HDT, Label='RUN_DT:', __RC__)
    call MAPL_GetResource(mapl, CDT, Label='GOCART_DT:', default=real(HDT), __RC__)
    self%CDT = CDT

!  Load resource file and get number of bins 
!  -------------------------------------------
    cfg = ESMF_ConfigCreate (__RC__)
    call ESMF_ConfigLoadFile (cfg, 'NI2G_GridComp_'//trim(COMP_NAME)//'.rc', rc=status)
    if (status /= 0) then
      if (mapl_am_i_root()) print*,'NI2G_GridComp_'//trim(COMP_NAME)//'.rc does not exist! &
                                    loading NI2G_GridComp_NI.rc instead'
      call ESMF_ConfigLoadFile( cfg, 'NI2G_GridComp_NI.rc', __RC__)
    end if

!   Call Generic Initialize 
!   ----------------------------------------
    call MAPL_GenericInitialize (GC, import, export, clock, __RC__)

!   Get parameters from generic state.
!   -----------------------------------
    call MAPL_Get ( mapl, INTERNAL_ESMF_STATE = internal, __RC__)

!   Is NI data driven?
!   ------------------
    call determine_data_driven (COMP_NAME, data_driven, __RC__)

!   Set HNO3 recycle alarm
    if (.not. data_driven) then
        call ESMF_ClockGet(clock, calendar=calendar, currTime=currentTime, __RC__)
        call ESMF_TimeGet(currentTime, YY=year, MM=month, DD=day, H=hh, M=mm, S=ss, __RC__)
        call ESMF_TimeSet(ringTime, YY=year, MM=month, DD=day, H=0, M=0, S=0, __RC__)
        call ESMF_TimeIntervalSet(ringInterval, H=3, calendar=calendar, __RC__)

        do while (ringTime < currentTime)! DO WE NEED THIS?
            ringTime = currentTime + ringInterval
        end do

        alarm_HNO3 = ESMF_AlarmCreate(Clock        = clock,        &
                                      Name         = 'HNO3_RECYCLE_ALARM', &
                                      RingInterval = ringInterval, &
                                      RingTime     = currentTime,  &
                                      Enabled      = .true.   ,    &
                                      Sticky       = .false.  , __RC__)
    end if



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

    call ESMF_StateGet (internal, 'NH3', field, __RC__)
    call ESMF_AttributeSet(field, NAME='ScavengingFractionPerKm', VALUE=self%fscav(1), __RC__)
    call ESMF_StateGet (internal, 'NH4a', field, __RC__)
    call ESMF_AttributeSet(field, NAME='ScavengingFractionPerKm', VALUE=self%fscav(2), __RC__)

!   Fill AERO State with N03an(1,2,3) fields
!   ----------------------------------------
    call ESMF_StateGet (export, trim(COMP_NAME)//'_AERO'    , aero    , __RC__)

    call ESMF_StateGet (internal, 'NO3an1', field, __RC__)
    call ESMF_AttributeSet(field, NAME='ScavengingFractionPerKm', VALUE=self%fscav(3), __RC__)
    fld = MAPL_FieldCreate (field, 'NO3an1', __RC__)
    call MAPL_StateAdd (aero, fld, __RC__)

    call ESMF_StateGet (internal, 'NO3an2', field, __RC__)
    call ESMF_AttributeSet(field, NAME='ScavengingFractionPerKm', VALUE=self%fscav(4), __RC__)
    fld = MAPL_FieldCreate (field, 'NO3an2', __RC__)
    call MAPL_StateAdd (aero, fld, __RC__)

    call ESMF_StateGet (internal, 'NO3an3', field, __RC__)
    call ESMF_AttributeSet(field, NAME='ScavengingFractionPerKm', VALUE=self%fscav(5), __RC__)
    fld = MAPL_FieldCreate (field, 'NO3an3', __RC__)
    call MAPL_StateAdd (aero, fld, __RC__)


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

!#if 0
!   Create Diagnostics Mie Table
!   -----------------------------
!   Get file names for the optical tables
    call ESMF_ConfigGetAttribute (cfg, self%diag_MieTable(instance)%optics_file, &
                                  label="aerosol_monochromatic_optics_file:", __RC__ )
    call ESMF_ConfigGetAttribute (cfg, self%diag_MieTable(instance)%nch, label="n_channels:", __RC__)
    call ESMF_ConfigGetAttribute (cfg, self%diag_MieTable(instance)%nmom, label="n_moments:", default=0,  __RC__)
    allocate (self%diag_MieTable(instance)%channels(self%diag_MieTable(instance)%nch), __STAT__ )
    call ESMF_ConfigGetAttribute (cfg, self%diag_MieTable(instance)%channels, &
                                  label="aerosol_monochromatic_optics_wavelength:", __RC__)
    allocate (self%diag_MieTable(instance)%mie_aerosol, __STAT__)
    self%diag_MieTable(instance)%mie_aerosol = Chem_MieTableCreate (self%diag_MieTable(instance)%optics_file, __RC__ )
    call Chem_MieTableRead (self%diag_MieTable(instance)%mie_aerosol, self%diag_MieTable(instance)%nch, &
                            self%diag_MieTable(instance)%channels, rc, nmom=self%diag_MieTable(instance)%nmom)

    ! Mie Table instance/index
    call ESMF_AttributeSet(aero, name='mie_table_instance', value=instance, __RC__)

    ! Add variables to NI instance's aero state. This is used in aerosol optics calculations
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

!    call ESMF_AttributeSet(aero, name='internal_varaible_name', value='SS', __RC__)

!    call ESMF_MethodAdd(AERO, label='aerosol_optics', userRoutine=aerosol_optics, __RC__)

!#endif
if(mapl_am_i_root()) print*,trim(comp_name),' Init END'



    RETURN_(ESMF_SUCCESS)

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

!if(mapl_am_i_root()) print*,'NI2G Run BEGIN'

!   Get my internal MAPL_Generic state
!   -----------------------------------
    call MAPL_GetObjectFromGC (GC, MAPL, __RC__)

!   Get parameters from generic state.
!   -----------------------------------
    call MAPL_Get (MAPL, INTERNAL_ESMF_STATE=internal, __RC__)

!   Is NI data driven?
!   ------------------
    call determine_data_driven (COMP_NAME, data_driven, __RC__)

!   Update INTERNAL state variables with ExtData
!   ---------------------------------------------
!    if (data_driven) then
!       call Run_data (GC, import, export, internal, __RC__)
!    else
       call Run1 (GC, import, export, clock, __RC__)
!    end if

!if(mapl_am_i_root()) print*,'NI2G Run END'

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

! !DESCRIPTION:  Computes emissions/sources for Nitrate

!EOP
!============================================================================
!   !Locals
    character (len=ESMF_MAXSTR)       :: COMP_NAME
    type (MAPL_MetaComp), pointer     :: mapl
    type (ESMF_State)                 :: internal
    type (ESMF_Grid)                  :: grid
    type (wrap_)                      :: wrap
    type (NI2G_GridComp), pointer     :: self

#include "NI2G_DeclarePointer___.h"

   __Iam__('Run1')

!*****************************************************************************
!   Begin... 

!   Get my name and set-up traceback handle
!   ---------------------------------------
    call ESMF_GridCompGet (GC, grid=grid, NAME=comp_name, __RC__)
    Iam = trim(comp_name) //'::'// Iam

!if(mapl_am_i_root()) print*,trim(comp_name),'2G Run1 BEGIN'

!   Get my internal MAPL_Generic state
!   -----------------------------------
    call MAPL_GetObjectFromGC (GC, mapl, __RC__)

!   Get parameters from generic state.
!   -----------------------------------
    call MAPL_Get (mapl, INTERNAL_ESMF_STATE=internal, __RC__)

#include "NI2G_GetPointer___.h"

!   Get my private internal state
!   ------------------------------
    call ESMF_UserCompGetInternalState(GC, 'NI2G_GridComp', wrap, STATUS)
    VERIFY_(STATUS)
    self => wrap%ptr

!   NH3 Emissions
!   -------------
    if (associated(NH3EM)) then
       NH3EM = 0.
       if (associated(EMI_NH3_BB)) NH3EM = NH3EM + EMI_NH3_BB
       if (associated(EMI_NH3_AG)) NH3EM = NH3EM + EMI_NH3_AG
       if (associated(EMI_NH3_EN)) NH3EM = NH3EM + EMI_NH3_EN
       if (associated(EMI_NH3_RE)) NH3EM = NH3EM + EMI_NH3_RE
       if (associated(EMI_NH3_TR)) NH3EM = NH3EM + EMI_NH3_TR
       if (associated(EMI_NH3_IN)) NH3EM = NH3EM + EMI_NH3_IN
       if (associated(EMI_NH3_OC)) NH3EM = NH3EM + EMI_NH3_OC
    end if

    if (associated(EMI_NH3_BB)) &
       NH3(:,:,self%km) = NH3(:,:,self%km)+self%cdt*chemgrav/delp(:,:,self%km)*EMI_NH3_BB
    if (associated(EMI_NH3_AG)) &
       NH3(:,:,self%km) = NH3(:,:,self%km)+self%cdt*chemgrav/delp(:,:,self%km)*EMI_NH3_AG
    if (associated(EMI_NH3_EN)) &
       NH3(:,:,self%km) = NH3(:,:,self%km)+self%cdt*chemgrav/delp(:,:,self%km)*EMI_NH3_EN
    if (associated(EMI_NH3_IN)) &
       NH3(:,:,self%km) = NH3(:,:,self%km)+self%cdt*chemgrav/delp(:,:,self%km)*EMI_NH3_IN
    if (associated(EMI_NH3_RE)) &
       NH3(:,:,self%km) = NH3(:,:,self%km)+self%cdt*chemgrav/delp(:,:,self%km)*EMI_NH3_RE
    if (associated(EMI_NH3_TR)) &
       NH3(:,:,self%km) = NH3(:,:,self%km)+self%cdt*chemgrav/delp(:,:,self%km)*EMI_NH3_TR
    if (associated(EMI_NH3_OC)) &
       NH3(:,:,self%km) = NH3(:,:,self%km)+self%cdt*chemgrav/delp(:,:,self%km)*EMI_NH3_OC


!if(mapl_am_i_root()) print*,'NI2G sum(DU)',sum(DU)



!if(mapl_am_i_root()) print*,trim(comp_name),'2G Run1 END'


    RETURN_(ESMF_SUCCESS)

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

! !DESCRIPTION: Run2 method for the Dust Grid Component.

!EOP
!============================================================================
! Locals
    character (len=ESMF_MAXSTR)       :: COMP_NAME
    type (MAPL_MetaComp), pointer     :: MAPL
    type (ESMF_State)                 :: internal
    type (wrap_)                      :: wrap
    type (NI2G_GridComp), pointer     :: self

    integer                           :: n
    real, allocatable, dimension(:,:) :: drydepositionfrequency, dqa
    real                              :: fwet
    logical                           :: KIN

    type (ESMF_ALARM)               :: alarm
    logical                         :: alarm_is_ringing

real :: rmedDU(5), rmedSS(5), fnumDU(5), fnumSS(5)

#include "NI2G_DeclarePointer___.h"

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
    call MAPL_Get (MAPL, INTERNAL_ESMF_STATE=INTERNAL, __RC__)

#include "NI2G_GetPointer___.h"

!if(mapl_am_i_root()) print*,trim(comp_name),'2G Run2 BEGIN'

!do n=1,5
!   if(mapl_am_i_root()) print*,'n = ', n,' : Run2 B SS2G sum(ss00n) = ',sum(SS(:,:,:,n))
!end do

!   Get my private internal state
!   ------------------------------
    call ESMF_UserCompGetInternalState(GC, 'NI2G_GridComp', wrap, STATUS)
    VERIFY_(STATUS)
    self => wrap%ptr

!   check hno3 alarm
    call ESMF_ClockGetAlarm(clock, 'HNO3_RECYCLE_ALARM', alarm, __RC__)
    alarm_is_ringing = ESMF_AlarmIsRinging(alarm, __RC__)

    if (alarm_is_ringing) then
       call ESMF_AlarmRingerOff(alarm, __RC__)
    end if

!   Save local copy of HNO3 for first pass through run method regardless
    if (self%first) then
       self%xhno3 = MAPL_UNDEF
       self%first = .false.
    end if

!   Recycle HNO3 every 3 hours
    if (alarm_is_ringing) then
       self%xhno3 = NITRATE_HNO3
    end if

!if(mapl_am_i_root()) print*,'NI2G sum(hno3) = ',sum(nitrate_hno3)
!if(mapl_am_i_root()) print*,'NI2G sum(self%xhno3) = ',sum(self%xhno3)
!if(mapl_am_i_root()) print*,'NI2G sum(DU) = ',sum(DU)
!if(mapl_am_i_root()) print*,'NI2G sum(SS) = ',sum(SS)

! This could be in incorrect alarm. This alarm is currently for HNO3_RECYCLE_ALARM, but a
! new alarm might need to be created just for this GC, or all of gocart2g?
!RUN_ALARM: if (alarm_is_ringing) then

!if(mapl_am_i_root()) print*,'NI2G - HNO3 ALARM HAS RUNG!!!!'

!FOR TESTING ONLY
rmedDU=(/0.73, 1.4, 2.4, 4.5, 8.0/)
fnumDU=(/2.45e14, 3.28e13, 6.52e12, 9.89e11, 1.76e11/)
rmedSS=(/0.079, 0.316, 1.119, 2.818, 7.772/)
fnumSS=(/3.017e17, 1.085e16, 1.207e14, 9.391e12, 2.922e11/)


    call NIheterogenousChem (NIHT, self%xhno3, MAPL_AVOGAD, MAPL_AIRMW, MAPL_PI, MAPL_RUNIV, &
                             airdens, t, rh2, delp, DU, SS, rmedDU*1.e-6, rmedSS*1.e-6, &
                             fnumDU, fnumSS, 5, 5, self%km, self%cdt, chemgrav, NO3an1, NO3an2, &
                             NO3an3, HNO3CONC, HNO3SMASS, HNO3CMASS, rc)

!if(mapl_am_i_root()) print*,'NI2G sum(NIHT(:,:,1)) = ',sum(NIHT(:,:,1))
!if(mapl_am_i_root()) print*,'NI2G sum(NIHT(:,:,2)) = ',sum(NIHT(:,:,2))
!if(mapl_am_i_root()) print*,'NI2G sum(NIHT(:,:,3)) = ',sum(NIHT(:,:,3))
!if(mapl_am_i_root()) print*,'NI2G sum(NO3an1) = ',sum(NO3an1)
!if(mapl_am_i_root()) print*,'NI2G sum(NO3an2) = ',sum(NO3an2)
!if(mapl_am_i_root()) print*,'NI2G sum(NO3an3) = ',sum(NO3an3)
!if(mapl_am_i_root()) print*,'NI2G sum(HNO3CONC) = ',sum(HNO3CONC)
!if(mapl_am_i_root()) print*,'NI2G sum(HNO3SMASS) = ',sum(HNO3SMASS)
!if(mapl_am_i_root()) print*,'NI2G sum(HNO3CMASS) = ',sum(HNO3CMASS)


!end if RUN_ALARM



!if(mapl_am_i_root()) print*,trim(comp_name),'2G Run2 END'



    RETURN_(ESMF_SUCCESS)
  
  end subroutine Run2




end module NI2G_GridCompMod
