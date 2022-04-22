#include "MAPL_Generic.h"

!=============================================================================
!BOP

! !MODULE: XX_GridCompMod - XX gridded component 

! !INTERFACE:
module XX_GridCompMod

!  !USES:
   use ESMF
   use MAPL
   use QuickChem_Generic
   use iso_c_binding, only: c_loc, c_f_pointer, c_ptr

   use QuickChem_Process       ! QuickChem process library
   use QC_EnvironmentMod
   use MAPL_StringTemplate, only: StrTemplate
   
   implicit none
   private

   integer, parameter :: instanceComputational = 1
   integer, parameter :: instanceData          = 2

! !PUBLIC MEMBER FUNCTIONS:
   public  SetServices

! !DESCRIPTION: This module implements QuickChem's template (XX) Gridded Component.

! !REVISION HISTORY:
! 16Oct2019  Sherman, da Silva, Darmenov, Clune -  First attempt at refactoring

!EOP
!===========================================================================

!  !XX state
   type, extends(QC_Environment) :: XX_GridComp
!      real, allocatable      :: rlow(:)        ! example of contents, was used in Run2
!      character(len=:), allocatable :: emission_scheme     ! See G2G DU for how this is used
!                                                           ! NOTE: it affects XX_Import___.h and XX_GetPointer___.h
   end type XX_GridComp

   type wrap_
      type (XX_GridComp), pointer     :: PTR => null()
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
    character (len=ESMF_MAXSTR)        :: COMP_NAME
    type (ESMF_Config)                 :: cfg
    type (ESMF_Config)                 :: universal_cfg
    type (wrap_)                       :: wrap
    type (XX_GridComp), pointer      :: self

    character (len=ESMF_MAXSTR)        :: field_name
    integer                            :: i
    logical                            :: data_driven = .true.

    __Iam__('SetServices')

!****************************************************************************
!   Begin...

!   Get my name and set-up traceback handle
!   ---------------------------------------
    call ESMF_GridCompGet (GC, NAME=COMP_NAME, config=universal_cfg, __RC__)
    IF (index(Iam,'::') == 0) Iam = trim(COMP_NAME)//'::'//Iam

!   Wrap internal state for storing in GC
!   -------------------------------------
    allocate (self, __STAT__)
    wrap%ptr => self

!   Load resource file  
!   -------------------
    cfg = ESMF_ConfigCreate (__RC__)
    call ESMF_ConfigLoadFile (cfg, 'XX_instance_'//trim(COMP_NAME)//'.rc', rc=status)

    if (status /= 0) then
        if (mapl_am_i_root()) print*,'XX_instance_'//trim(COMP_NAME)//'.rc does not exist! Loading XX_GridComp_XX.rc instead'
        call ESMF_ConfigLoadFile (cfg, 'XX_instance_XX.rc', __RC__)
    end if

    ! process generic config items
    call self%QC_Environment%load_from_config(cfg, universal_cfg, __RC__)

!   allocate( self%rlow(self%nbins),  __STAT__)  !! example

    ! process DU-specific items
!   call ESMF_ConfigGetAttribute (cfg, self%rlow,       label='radius_lower:', __RC__)   !! example of reading a vector of values from XX_instance_XX.rc


!   Is XX data driven?
!   ------------------
    call determine_data_driven (COMP_NAME, data_driven, __RC__)

!   Set entry points
!   ------------------------
    call MAPL_GridCompSetEntryPoint (GC, ESMF_Method_Initialize,  Initialize, __RC__)
    call MAPL_GridCompSetEntryPoint (GC, ESMF_Method_Run, Run, __RC__)
    if (data_driven .neqv. .true.) then
       call MAPL_GridCompSetEntryPoint (GC, ESMF_Method_Run, Run2, __RC__)
    end if


!   IMPORT STATE
!   -------------
    if (data_driven) then

       call MAPL_AddInternalSpec(gc,&
          short_name='XX', &
          long_name='Sample QuickChem species', &
          units='???', &
          dims=MAPL_DimsHorzVert, &
          vlocation=MAPL_VlocationCenter, &
          restart=MAPL_RestartOptional, &
          ungridded_dims=[self%nbins], &
          friendlyto='DYNAMICS:TURBULENCE:MOIST', &
          add2export=.true., __RC__)

!      Pressure at layer edges
!      -----------------------
       call MAPL_AddImportSpec(GC,                            &
          SHORT_NAME = 'PLE',                                 &
          LONG_NAME  = 'air_pressure',                        &
          UNITS      = 'Pa',                                  &
          DIMS       = MAPL_DimsHorzVert,                     &
          VLOCATION  = MAPL_VLocationEdge,                    &
          RESTART    = MAPL_RestartSkip,     __RC__)

!!      RH: is between 0 and 1
!!      ----------------------
!       call MAPL_AddImportSpec(GC,                            &
!          SHORT_NAME = 'RH2',                                 &
!          LONG_NAME  = 'Rel_Hum_after_moist',                 &
!          UNITS      = '1',                                   &
!          DIMS       = MAPL_DimsHorzVert,                     &
!          VLOCATION  = MAPL_VLocationCenter,                  &
!          RESTART    = MAPL_RestartSkip,     __RC__)
 
       do i = 1, self%nbins 
          write (field_name, '(A, I0.3)') '', i
          call MAPL_AddImportSpec(GC,                                        &
             short_name = 'climxx'//trim(field_name),                        &
             long_name  = 'XX Mixing Ratio (bin '//trim(field_name)//')',    &
             units      = 'kg kg-1',                                         &
             restart    = MAPL_RestartSkip,                                  &
             dims       = MAPL_DimsHorzVert,                                 &
             vlocation  = MAPL_VLocationCenter, __RC__)

!  Other imports (in DUST) were:
!            short_name = 'climXXDP'//trim(field_name),                     &
!            short_name = 'climXXWT'//trim(field_name),                     &
!            short_name = 'climXXSD'//trim(field_name),                     &
!            short_name = 'climXXSV'//trim(field_name),                     &

        end do
    end if ! (data_driven)

!   Computational Instance
!   ----------------------
    if (.not. data_driven) then

#include "XX_Export___.h"

#include "XX_Import___.h"

#include "XX_Internal___.h"

    end if

!   Store internal state in GC
!   --------------------------
    call ESMF_UserCompSetInternalState ( GC, 'XX_GridComp', wrap, STATUS )
    VERIFY_(STATUS)

!   Set generic services
!   ----------------------------------
    call MAPL_GenericSetServices (GC, __RC__)

    RETURN_(ESMF_SUCCESS)

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

! !DESCRIPTION: This initializes XX's Grid Component.

! !REVISION HISTORY: 
! 16oct2019  E.Sherman, A.da Silva, T.Clune, A.Darmenov - First attempt at refactoring

!EOP
!============================================================================
!   !Locals
    character (len=ESMF_MAXSTR)          :: COMP_NAME
    type (MAPL_MetaComp),       pointer  :: MAPL
    type (ESMF_Grid)                     :: grid
    type (ESMF_State)                    :: internal
!   type (ESMF_State)                    :: providerState
    type (ESMF_Config)                   :: cfg
    type (ESMF_Config)                   :: universal_cfg
    type (wrap_)                         :: wrap
    type (XX_GridComp), pointer          :: self

    integer                              :: i, dims(3), km
    integer                              :: instance
    type (ESMF_Field)                    :: field, fld
    character (len=ESMF_MAXSTR)          :: bin_index, prefix
    real                                 :: CDT         ! chemistry timestep (secs)
    integer                              :: HDT         ! model     timestep (secs)
    real, pointer, dimension(:,:,:)      :: int_ptr
!   real, pointer, dimension(:,:,:)      :: ple
    logical                              :: data_driven

    __Iam__('Initialize')

!****************************************************************************
!   Begin... 

!   Get the target components name and set-up traceback handle.
!   -----------------------------------------------------------
    call ESMF_GridCompGet (GC, grid=grid, name=COMP_NAME, config=universal_cfg, __RC__)
    IF (index(Iam,'::') == 0) Iam = trim(COMP_NAME)//'::'//Iam

!   Get my internal MAPL_Generic state
!   -----------------------------------
    call MAPL_GetObjectFromGC (GC, MAPL, __RC__)

!   Get my internal private state
!   -----------------------------
    call ESMF_UserCompGetInternalState(GC, 'XX_GridComp', wrap, STATUS)
    VERIFY_(STATUS)
    self => wrap%ptr

    call MAPL_GridGet ( grid, globalCellCountPerDim=dims, __RC__ )

!   Dust emission tuning coefficient [kg s2 m-5]. NOT bin specific.
!   ---------------------------------------------------------------
!   self%Ch_DU = Chem_UtilResVal(dims(1), dims(2), self%Ch_DU_res(:), __RC__)  !! example of using res vec

!   Get dimensions
!   ---------------
    km = dims(3)
    self%km = km

!   Get DTs
!   -------
    call MAPL_GetResource(mapl, HDT, Label='RUN_DT:', __RC__)                        
    call MAPL_GetResource(mapl, CDT, Label='QUICKCHEM_DT:', default=real(HDT), __RC__)
    self%CDT = CDT

!   Load resource file  
!   -------------------
    cfg = ESMF_ConfigCreate (__RC__)
    call ESMF_ConfigLoadFile (cfg, 'XX_instance_'//trim(COMP_NAME)//'.rc', RC=STATUS)
    if (status /= 0) then
        if (mapl_am_i_root()) print*,'XX_instance_'//trim(COMP_NAME)//'.rc does not exist! &
                                      loading XX_instance_XX.rc instead'
        call ESMF_ConfigLoadFile (cfg, 'XX_instance_XX.rc', __RC__)
    end if

!   Call Generic Initialize 
!   ------------------------
    call MAPL_GenericInitialize (GC, import, export, clock, __RC__)

!   Get parameters from generic state.
!   -----------------------------------
    call MAPL_Get (MAPL, INTERNAL_ESMF_STATE=internal, __RC__)

!   Is XX data driven?
!   ------------------
    call determine_data_driven (COMP_NAME, data_driven, __RC__)

!   If this is a data component, the data is provided in the import
!   state via ExtData instead of the actual QuickChem children
!   ---------------------------------------------------------------- 
!   if ( data_driven ) then
!       providerState = import
!       prefix = 'clim'
!   else
!       providerState = export
!       prefix = ''
!   end if

!   Add attribute information for XX export.
!   call ESMF_StateGet (export, 'XX', field, __RC__)
!   call ESMF_AttributeSet(field, NAME='radius', valueList=self%radius, itemCount=self%nbins, __RC__)
!   call ESMF_AttributeSet(field, NAME='fnum', valueList=self%fnum, itemCount=self%nbins, __RC__)

!   Add attribute information to internal state variables
!   -----------------------------------------------------

!   call ESMF_StateGet (internal, 'XX', field, __RC__)
!   call ESMF_AttributeSet(field, NAME='ScavengingFractionPerKm', value=self%fscav(1), __RC__)

!   if (.not. data_driven) then
!
!      Set internal XX values to 0 where above klid
!      call MAPL_GetPointer (internal, int_ptr, 'XX', __RC__)
!      call setZeroKlid4d (self%km, self%klid, int_ptr)         ! see GOCART2G Generic F90
!
!   end if

    if (data_driven) then
       instance = instanceData
    else
       instance = instanceComputational
    end if

    self%instance = instance

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

!   !DESCRIPTION: Run method for the Dust Grid Component. Determines whether to
!                 run data or computational run method.

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
    IF (index(Iam,'::') == 0) Iam = trim(COMP_NAME)//'::'//Iam

!   Get my internal MAPL_Generic state
!   -----------------------------------
    call MAPL_GetObjectFromGC (GC, MAPL, __RC__)

!   Get parameters from generic state.
!   -----------------------------------
    call MAPL_Get (MAPL, INTERNAL_ESMF_STATE=internal, __RC__)

!   Is this species data driven?
!   ----------------------------
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

! !DESCRIPTION:  Computes emissions/sources for this species

!EOP
!============================================================================
!   !Locals
    character (len=ESMF_MAXSTR)       :: COMP_NAME
    type (MAPL_MetaComp), pointer     :: mapl
    type (ESMF_State)                 :: internal
    type (ESMF_Grid)                  :: grid
    type (wrap_)                      :: wrap
    type (XX_GridComp), pointer     :: self
    type(ESMF_Time)                   :: time

    integer  :: nymd, nhms, iyr, imm, idd, ihr, imn, isc
    integer  :: import_shape(2), i2, j2
!    real, dimension(:,:), pointer     :: emissions_surface
!    real, dimension(:,:,:), allocatable     :: emissions_surface
!    real, dimension(:,:,:,:), allocatable :: emissions
!    real, dimension(:,:,:), allocatable   :: emissions_point
!    character (len=ESMF_MAXSTR)  :: fname ! file name for point source emissions
!    integer, pointer, dimension(:)  :: iPoint, jPoint
!    logical :: fileExists
!    real :: qmax, qmin
    integer :: n, ijl
!    real, dimension(:,:), allocatable   :: z_
!    real, dimension(:,:), allocatable   :: ustar_
!    real, dimension(:,:), allocatable   :: ustar_t_
!    real, dimension(:,:), allocatable   :: ustar_ts_
!    real, dimension(:,:), allocatable   :: R_
!    real, dimension(:,:), allocatable   :: H_w_
!    real, dimension(:,:), allocatable   :: f_erod_

#include "XX_DeclarePointer___.h"

   __Iam__('Run1')

!*****************************************************************************
!   Begin... 

!   Get my name and set-up traceback handle
!   ---------------------------------------
    call ESMF_GridCompGet (GC, grid=grid, NAME=COMP_NAME, __RC__)
    IF (index(Iam,'::') == 0) Iam = trim(COMP_NAME)//'::'//Iam

!   Get my internal MAPL_Generic state
!   -----------------------------------
    call MAPL_GetObjectFromGC (GC, mapl, __RC__)

!   Get parameters from generic state.
!   -----------------------------------
    call MAPL_Get (mapl, INTERNAL_ESMF_STATE=internal, __RC__)

!   Get my private internal state
!   ------------------------------
    call ESMF_UserCompGetInternalState(GC, 'XX_GridComp', wrap, STATUS)
    VERIFY_(STATUS)
    self => wrap%ptr

!   Extract nymd(yyyymmdd) from clock
!   ---------------------------------
    call ESMF_ClockGet (clock, currTime=time, __RC__)
    call ESMF_TimeGet (time ,YY=iyr, MM=imm, DD=idd, H=ihr, M=imn, S=isc, __RC__)
    call MAPL_PackTime (nymd, iyr, imm , idd)
    call MAPL_PackTime (nhms, ihr, imn, isc)

#include "XX_GetPointer___.h"

!   Set xx_src to 0 where undefined
!   --------------------------------
    where (1.01*xx_src > MAPL_UNDEF) xx_src = 0.

!   Get dimensions
!   ---------------
    import_shape = shape(lwi)
    i2 = import_shape(1)
    j2 = import_shape(2)
    ijl  = ( i2 - 1 + 1 ) * ( j2 - 1 + 1 )

!   EMISSIONS:
!   Assume xx_src units are 'per second'
!   For more options, including point emissions, see GOCART2G species like DU
!   XX does not support nbins>1, so use 3 indices instead of 4:
    XX(:,:,self%km) = XX(:,:,self%km) + &
                      xx_src * (self%CDT * MAPL_GRAV / delp(:,:,self%km))

!   from here, G2G species called DustEmissionGOCART2G and UpdateAerosolState

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
    type (XX_GridComp), pointer     :: self

#include "XX_DeclarePointer___.h"

    __Iam__('Run2')

!*****************************************************************************
!   Begin... 

!   Get my name and set-up traceback handle
!   ---------------------------------------
    call ESMF_GridCompGet (GC, NAME=COMP_NAME, __RC__)
    IF (index(Iam,'::') == 0) Iam = trim(COMP_NAME)//'::'//Iam

!   Get my internal MAPL_Generic state
!   -----------------------------------
    call MAPL_GetObjectFromGC (GC, MAPL, __RC__)

!   Get parameters from generic state.
!   -----------------------------------
    call MAPL_Get (MAPL, INTERNAL_ESMF_STATE=internal, __RC__)

!   Get my private internal state
!   ------------------------------
    call ESMF_UserCompGetInternalState(GC, 'XX_GridComp', wrap, STATUS)
    VERIFY_(STATUS)
    self => wrap%ptr

#include "XX_GetPointer___.h"

!!  In G2G, here is where they call Chem_Settling, and (for each bin) DryDeposition

!!  In G2G, here is where they call WetRemovalGOCART2G for each bin, and Aero_Compute_Diags

    RETURN_(ESMF_SUCCESS)

  end subroutine Run2

!============================================================================
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

! Locals
    character (len=ESMF_MAXSTR)        :: COMP_NAME
    type (wrap_)                       :: wrap
    type (XX_GridComp), pointer        :: self

    integer                            :: i
    character (len=ESMF_MAXSTR)        :: field_name

    real, pointer, dimension(:,:,:)    :: ptr3d_int
    real, pointer, dimension(:,:,:)    :: ptr3d_imp

    __Iam__('Run_data')

!EOP
!*****************************************************************************
!   Begin... 

!   Get my name and set-up traceback handle
!   ---------------------------------------
    call ESMF_GridCompGet (GC, NAME=COMP_NAME, __RC__)
    IF (index(Iam,'::') == 0) Iam = trim(COMP_NAME)//'::'//Iam

!   Get my private internal state
!   ------------------------------
    call ESMF_UserCompGetInternalState(GC, 'XX_GridComp', wrap, STATUS)
    VERIFY_(STATUS)
    self => wrap%ptr

!   Update internal data pointers with ExtData
!   ------------------------------------------
    call MAPL_GetPointer (internal, NAME='XX', ptr=ptr3d_int, __RC__)

    IF ( self%nbins /= 1 ) THEN
      PRINT*,'expecting only 1 XX bin'
      RETURN_(123)
    ENDIF

    do i = 1, self%nbins
    write (field_name, '(A, I0.3)') 'xx', i
        call MAPL_GetPointer (import,  NAME='clim'//trim(field_name), ptr=ptr3d_imp, __RC__)

!       ptr4d_int(:,:,:,i) = ptr3d_imp
        ptr3d_int(:,:,:  ) = ptr3d_imp
    end do

    RETURN_(ESMF_SUCCESS)

  end subroutine Run_data

!-------------------------------------------------------------------------------------

! INTERESTING THINGS
!
!   _ASSERT(any(self%emission_scheme == [character(len=7) :: 'ginoux','other']), "Error. Unallowed emission scheme: "//trim(self%emission_scheme)//". Allowed: ginoux, other")

!     associate (scheme => self%emission_scheme)
! #include "XX_Import___.h"
!     end associate

!     associate (scheme => self%emission_scheme)
! #include "XX_GetPointer___.h"
    ! end associate

end module XX_GridCompMod
