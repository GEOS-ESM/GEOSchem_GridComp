#include "MAPL_Generic.h"

!=============================================================================
!BOP

! !MODULE: GOCART2G_GridCompMod - The GOCART 2nd Generation Aerosol Grid Component

! !INTERFACE:

module GOCART2G_GridCompMod

! !USES:

   USE ESMF
   USE MAPL_Mod

! !Establish the Childen's SetServices
! !-----------------------------------
   use DU2G_GridCompMod,    only   : DU2GSetServices  => SetServices
   use SSng_GridCompMod,    only   : SSngSetServices  => SetServices
!   use SUng_GridCompMod,    only   : SUngSetServices  => SetServices
!   use BCng_GridCompMod,    only   : BCngSetServices  => SetServices
!   use OCng_GridCompMod,    only   : OCngSetServices  => SetServices
!   use NIng_GridCompMod,    only   : NIngSetServices  => SetServices

   implicit none
   private

! !PUBLIC MEMBER FUNCTIONS:
   PUBLIC  SetServices

! Private State
  type GOCART_State
     private
     character (len=ESMF_MAXSTR), pointer    :: instances_DU(:) => null()
     character (len=ESMF_MAXSTR), pointer    :: instances_SS(:) => null()
     character (len=ESMF_MAXSTR), pointer    :: instances_SU(:) => null()
     character (len=ESMF_MAXSTR), pointer    :: instances_BC(:) => null()
     character (len=ESMF_MAXSTR), pointer    :: instances_OC(:) => null()
     character (len=ESMF_MAXSTR), pointer    :: instances_NI(:) => null()
     integer                                 :: active_instances_DU = 0
     integer                                 :: active_instances_SS = 0
     integer                                 :: active_instances_SU = 0
     integer                                 :: active_instances_BC = 0
     integer                                 :: active_instances_OC = 0
     integer                                 :: active_instances_NI = 0
  end type GOCART_State

  type wrap_
     type (GOCART_State), pointer     :: PTR => null()
  end type wrap_


! !DESCRIPTION:
!
!   {\tt GOCART} is a gridded component from the GOCART model and includes 
!  dust, sea salt, sulfates, nitrate, organic and black carbon. 
!
!
! !REVISION HISTORY:
!
!  25feb2005  da Silva   First crack.
!  19jul2006  da Silva   First separate GOCART component.
!  14Oct2019  E.Sherman, A.Darmenov, da Silva  First attempt at refactoring. 
!
!EOP
!-------------------------------------------------------------------------

  integer ::     DU2G = -1
  integer ::     SSng = -1
  integer ::     SUng = -1
  integer ::     BCng = -1
  integer ::     OCng = -1
  integer ::     NIng = -1


contains

!BOP

! !IROUTINE: SetServices -- Sets ESMF services for this component

! !INTERFACE:

  subroutine SetServices (GC, RC)

! !ARGUMENTS:

    type (ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
    integer, optional                   :: RC  ! return code

! !DESCRIPTION: This version uses MAPL_GenericSetServices, which sets
!   the Initialize and Finalize services to generic versions. It also
!   allocates our instance of a generic state and puts it in the 
!   gridded component (GC). Here we only set the two-stage run method and
!   declare the data services.

! !REVISION HISTORY: 
!  14oct2019  Sherman, da Silva, Darmenov, Clune - First attempt at refactoring for ESMF compatibility


!EOP

!****************************************************************************
!
! ErrLog Variables
    character (len=ESMF_MAXSTR)                   :: IAm
    integer                                       :: STATUS
    character (len=ESMF_MAXSTR)                   :: COMP_NAME

! Locals
    type (ESMF_Config)                            :: myCF
    type (GOCART_State), pointer                  :: self
    type (wrap_)                                  :: wrap

    integer                                       :: i,nq

!****************************************************************************

! Begin...

!   Get my name and set-up traceback handle
!   ---------------------------------------
    call ESMF_GridCompGet (GC, NAME=COMP_NAME, __RC__)
    Iam = trim(COMP_NAME) // '::' //  'SetServices'


!   Wrap internal state for storing in GC
!   -------------------------------------
    allocate (self, __STAT__)
    wrap%ptr => self

!   Set the Initialize, Run, Finalize entry points
!   ------------------------------------------------
    call MAPL_GridCompSetEntryPoint (GC, ESMF_METHOD_INITIALIZE,  Initialize,  __RC__)
    call MAPL_GridCompSetEntryPoint (GC, ESMF_METHOD_RUN,  Run1, __RC__)
    call MAPL_GridCompSetEntryPoint (GC, ESMF_METHOD_RUN,  Run2, __RC__)

!   Store internal state in GC
!   --------------------------
    call ESMF_UserCompSetInternalState (GC, 'GOCART_State', wrap, STATUS)
    VERIFY_(STATUS)

!   Get instances to determine what children will be born
!   -----------------------------------------------------
    myCF = ESMF_ConfigCreate (__RC__)
    call ESMF_ConfigLoadFile (myCF, 'GOCART2G_GridComp.rc', __RC__)

    call getInstances_('DU', myCF, instances=self%instances_DU, active_instances=self%active_instances_DU, __RC__)
    call getInstances_('SS', myCF, instances=self%instances_SS, active_instances=self%active_instances_SS, __RC__)
!    call getInstances_('SU', myCF, instances=self%instances_SU, active_instances=self%active_instances_SU, __RC__)
!    call getInstances_('BC', myCF, instances=self%instances_BC, active_instances=self%active_instances_BC, __RC__)
!    call getInstances_('OC', myCF, instances=self%instances_OC, active_instances=self%active_instances_OC, __RC__)
!    call getInstances_('NI', myCF, instances=self%instances_NI, active_instances=self%active_instances_NI, __RC__)

    call ESMF_ConfigDestroy(myCF, __RC__)

!   Create children`s gridded components and invoke their SetServices
!   -----------------------------------------------------------------
    do i = 1, size(self%instances_DU)
        DU2G = MAPL_AddChild(GC, NAME=trim(self%instances_DU(i)), SS=DU2GSetServices, __RC__)
    end do

    do i = 1, size(self%instances_SS)
        SSng = MAPL_AddChild(GC, NAME=trim(self%instances_SS(i)), SS=SSngSetServices, __RC__)
    end do

!    do i = 1, size(self%instances_SU)
!        SUng = MAPL_AddChild(GC, NAME=trim(self%instances_SU(i)), SS=SUngSetServices, __RC__)
!    end do

!    do i = 1, size(self%instances_BC)
!        BCng = MAPL_AddChild(GC, NAME=trim(self%instances_BC(i)), SS=BCngSetServices, __RC__)
!    end do

!    do i = 1, size(self%instances_OC)
!        OCng = MAPL_AddChild(GC, NAME=trim(self%instances_OC(i)), SS=OCngSetServices, __RC__)
!    end do

!    do i = 1, size(self%instances_NI)
!        NIng = MAPL_AddChild(GC, NAME=trim(self%instances_(i)), SS=NIngSetServices, __RC__)
!    end do


!   Define EXPORT states
!   This state is needed by radiation - It will contain 
!   aerosols and aerosol optics
!   --------------------------------------------------------
    call MAPL_AddExportSpec(GC,                       &
       SHORT_NAME = 'AERO2G',                         &
       LONG_NAME  = 'aerosol_mass_mixing_ratios_ng',  &
       UNITS      = 'kg kg-1',                        &
       DIMS       = MAPL_DimsHorzVert,                &
       VLOCATION  = MAPL_VLocationCenter,             &
       DATATYPE   = MAPL_StateItem, __RC__)

!   This state is needed by MOIST - It will contain
!   aerosols
!   --------------------------------------------------------
    call MAPL_AddExportSpec(GC,                       &
       SHORT_NAME = 'AERO2G_ACI',                     &
       LONG_NAME  = 'aerosol_cloud_interaction_ng',   &
       UNITS      = 'kg kg-1',                        &
       DIMS       = MAPL_DimsHorzVert,                &
       VLOCATION  = MAPL_VLocationCenter,             &
       DATATYPE   = MAPL_StateItem, __RC__)

!   This bundle is needed by surface for snow albedo modification
!   by aerosol settling and deposition
!   ~~~DEVELOPERS NOTE~~~ Change to StateItem when possible
!                         This will require refactoring Radiation
!   --------------------------------------------------------
    call MAPL_AddExportSpec(GC,                       &
       SHORT_NAME = 'AERO2G_DP',                      &
       LONG_NAME  = 'aerosol_deposition_ng',          &
       UNITS      = 'kg m-2 s-1',                     &
       DIMS       = MAPL_DimsHorzOnly,                &
       DATATYPE   = MAPL_BundleItem, __RC__)



!   Set generic services
!   ----------------------------------
    call MAPL_GenericSetServices (GC, __RC__)


    RETURN_(ESMF_SUCCESS)

  end subroutine SetServices


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!BOP

! !IROUTINE: Initialize -- Initialize method for the composite Gridded Component

! !INTERFACE:

  subroutine Initialize (GC, IMPORT, EXPORT, CLOCK, RC)

! !ARGUMENTS:

    type (ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
    type (ESMF_State),    intent(inout) :: IMPORT ! Import state
    type (ESMF_State),    intent(inout) :: EXPORT ! Export state
    type (ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
    integer, optional,    intent(  out) :: RC     ! Error code

! !DESCRIPTION:  This initializes the GOCART Grid Component. It primarily creates
!                its exports and births its children.

! !REVISION HISTORY: 
! 14oct2019   E.Sherman  First attempt at refactoring


!EOP

!****************************************************************************
! ErrLog Variables

    character (len=ESMF_MAXSTR)              :: IAm
    integer                                  :: STATUS
    character (len=ESMF_MAXSTR)              :: COMP_NAME

! Local derived type aliases
    type (MAPL_MetaComp),       pointer      :: MAPL
    type (ESMF_GridComp),       pointer      :: GCS(:)
    type (ESMF_State),          pointer      :: GEX(:)
    type (ESMF_State)                        :: INTERNAL
    type (ESMF_Grid)                         :: grid

    type (ESMF_State)                        :: AERO, AERO_ACI
    type (ESMF_FieldBundle)                  :: AERO_DP, child_bundle
    type (ESMF_State)                        :: child_state
    type (ESMF_Field), allocatable           :: fieldList(:)
    type (ESMF_Field)                        :: field
    character (len=ESMF_MAXSTR)              :: field_name

    type (GOCART_State),      pointer        :: self
    type (wrap_)                             :: wrap

    character (len=ESMF_MAXSTR)              :: CHILD_NAME
    character (len=ESMF_MAXSTR), allocatable :: AEROlist(:)

    integer                                  :: i, j, k, fieldCount

!****************************************************************************

! Begin... 

!   Get the target components name and set-up traceback handle.
!   -----------------------------------------------------------
    call ESMF_GridCompGet (GC, grid=grid, name=COMP_NAME, __RC__)
    Iam = trim(COMP_NAME) // '::' // 'Initialize'

   if (MAPL_AM_I_ROOT()) then
       print *, TRIM(Iam)//': Starting...'
       print *,' '
   end if

!   Get my internal MAPL_Generic state
!   -----------------------------------
    call MAPL_GetObjectFromGC (GC, MAPL, __RC__)

!   Call Generic Initialize
!   ----------------------------------------
    call MAPL_GenericInitialize (GC, IMPORT, EXPORT, CLOCK, __RC__)

!   Get my internal state
!   ---------------------
    call ESMF_UserCompGetInternalState (GC, 'GOCART_State', wrap, STATUS)
    VERIFY_(STATUS)
    self => wrap%ptr

!   Get children and their export states from my generic state
!   -----------------------------------------------------------
    call MAPL_Get (MAPL, GCS=GCS, GEX=GEX, __RC__ )

!   Fill AERO, AERO_ACI, and AERO_DP with the children's states
!   ------------------------------------------------------------
    call ESMF_StateGet (EXPORT, 'AERO2G'     , AERO     , __RC__)
    call ESMF_StateGet (EXPORT, 'AERO2G_ACI' , AERO_ACI , __RC__)
    call ESMF_StateGet (EXPORT, 'AERO2G_DP'  , AERO_DP  , __RC__)


!   Add children's AERO states to GOCART2G's AERO states
!   -----------------------------------------------------
    do i = 1, size(GCS)
        call ESMF_GridCompGet (GCS(i), NAME=CHILD_NAME, __RC__ )
                
        call ESMF_StateGet (GEX(i), trim(CHILD_NAME)//'_AERO', child_state, __RC__)
        call ESMF_StateAdd (AERO, (/child_state/), __RC__)

        call ESMF_StateGet (GEX(i), trim(CHILD_NAME)//'_AERO_ACI', child_state, __RC__)
        call ESMF_StateAdd (AERO_ACI, (/child_state/), __RC__)

        call ESMF_StateGet (GEX(i), 'AERO_DP', child_bundle, __RC__)
        call ESMF_FieldBundleGet (child_bundle, fieldCount=fieldCount, __RC__)
        allocate (fieldList(FieldCount), __STAT__)
        call ESMF_FieldBundleGet (child_bundle, fieldList=fieldList, __RC__)
        call ESMF_FieldBundleAdd (AERO_DP, fieldList, __RC__)
        deallocate(fieldList, __STAT__)
    end do


    ! This attribute indicates if the aerosol optics method is implemented or not. 
    ! Radiation will not call the aerosol optics method unless this attribute is 
    ! explicitly set to true.
    call ESMF_AttributeSet(aero, name='implements_aerosol_optics_method', value=.true., __RC__)

    ! set list of active aerosol instances 
    call getAERO_ (self, AEROlist, __RC__)
    call ESMF_AttributeSet(AERO, name='active_aerosol_instances', valueList=AEROlist, itemCount=size(AEROlist), __RC__)

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
    call ESMF_MethodAdd(AERO, label='run_aerosol_optics', userRoutine=run_aerosol_optics, __RC__)


    RETURN_(ESMF_SUCCESS)

  end subroutine Initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!BOP
! !IROUTINE: RUN -- Run method for the CONVECT component

! !INTERFACE:

  subroutine Run1 (GC, IMPORT, EXPORT, CLOCK, RC)

    ! !ARGUMENTS:

    type (ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
    type (ESMF_State),    intent(inout) :: IMPORT ! Import state
    type (ESMF_State),    intent(inout) :: EXPORT ! Export state
    type (ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
    integer, optional,    intent(  out) :: RC     ! Error code:

! !DESCRIPTION: This version uses the MAPL\_GenericSetServices. This function sets
!                the Initialize and Finalize services, as well as allocating

!EOP

!****************************************************************************
! ErrLog Variables

    character(len=ESMF_MAXSTR)          :: IAm
    integer                             :: STATUS
    character(len=ESMF_MAXSTR)          :: COMP_NAME

! Local derived type aliases

    type (MAPL_MetaComp),      pointer  :: MAPL
    type (ESMF_GridComp),      pointer  :: GCS(:)
    type (ESMF_State),         pointer  :: GIM(:)
    type (ESMF_State),         pointer  :: GEX(:)
    type (ESMF_State)                   :: INTERNAL

    integer                             :: i

!****************************************************************************
! Begin... 

!   Get my name and set-up traceback handle
!   ---------------------------------------
    Iam = 'Run1'
    call ESMF_GridCompGet( GC, NAME=COMP_NAME, __RC__ )
    Iam = trim(COMP_NAME) // Iam

!   Get my internal MAPL_Generic state
!   -----------------------------------
    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS )
    VERIFY_(STATUS)

!   Get parameters from generic state.
!   -----------------------------------
    call MAPL_Get ( MAPL, GCS=GCS, GIM=GIM, GEX=GEX, INTERNAL_ESMF_STATE=INTERNAL, __RC__ )

!   Run the children
!   -----------------
    do i = 1, size(GCS)
      call ESMF_GridCompRun (GCS(i), importState=GIM(i), exportState=GEX(i), clock=CLOCK, __RC__)
    end do


    RETURN_(ESMF_SUCCESS)

  end subroutine Run1


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!BOP
! !IROUTINE: RUN -- Run method for the CONVECT component

! !INTERFACE:

  subroutine Run2 (GC, IMPORT, EXPORT, CLOCK, RC)

    ! !ARGUMENTS:

    type (ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
    type (ESMF_State),    intent(inout) :: IMPORT ! Import state
    type (ESMF_State),    intent(inout) :: EXPORT ! Export state
    type (ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
    integer, optional,    intent(  out) :: RC     ! Error code:

! !DESCRIPTION: This version uses the MAPL\_GenericSetServices. This function sets
!                the Initialize and Finalize services, as well as allocating

!EOP

!****************************************************************************
! ErrLog Variables

    character(len=ESMF_MAXSTR)          :: IAm
    integer                             :: STATUS
    character(len=ESMF_MAXSTR)          :: COMP_NAME

! Local derived type aliases

    type (MAPL_MetaComp),      pointer  :: MAPL
    type (ESMF_GridComp),      pointer  :: GCS(:)
    type (ESMF_State),         pointer  :: GIM(:)
    type (ESMF_State),         pointer  :: GEX(:)
    type (ESMF_State)                   :: INTERNAL

    integer                             :: i

!****************************************************************************
! Begin... 

!   Get my name and set-up traceback handle
!   ---------------------------------------
    Iam = 'Run2'
    call ESMF_GridCompGet( GC, NAME=COMP_NAME, __RC__ )
    Iam = trim(COMP_NAME) // Iam


!   Get my internal MAPL_Generic state
!   -----------------------------------
    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS )
    VERIFY_(STATUS)

!   Get parameters from generic state.
!   -----------------------------------
    call MAPL_Get ( MAPL, GCS=GCS, GIM=GIM, GEX=GEX, INTERNAL_ESMF_STATE=INTERNAL, __RC__ )

!   Run the children
!   -----------------
    do i = 1, size(GCS)
      call ESMF_GridCompRun (GCS(i), importState=GIM(i), exportState=GEX(i), phase=2, clock=CLOCK, __RC__)
    end do

    RETURN_(ESMF_SUCCESS)

  end subroutine Run2


!===============================================================================

  subroutine getInstances_ (aerosol, myCF, instances, active_instances, rc)

!   Description: Fills the GOCART_State (aka, self%instance_XX) with user
!                defined instances from the GOCART2G_GridComp.rc.

    implicit none

    character (len=*),                intent(in   )  :: aerosol
    type (ESMF_Config),               intent(inout)  :: myCF
    integer,                          intent(inout)  :: active_instances
    integer,                          intent(  out)  :: rc
    character (len=ESMF_MAXSTR), pointer             :: instances(:)

!   locals
    integer                                          :: STATUS
    character (len=ESMF_MAXSTR)                      :: Iam = 'GOCART2G::getInstances_'

    integer                                          :: i, n_inst_act, n_inst_pass
    character (len=ESMF_MAXSTR)                      :: inst_name
!--------------------------------------------------------------------------------------

!   Begin...
    n_inst_act  = ESMF_ConfigGetLen (myCF, label='ACTIVE_INSTANCES_'//trim(aerosol)//':', __RC__)
    n_inst_pass = ESMF_ConfigGetLen (myCF, label='PASSIVE_INSTANCES_'//trim(aerosol)//':', __RC__)
    allocate (instances(n_inst_act + n_inst_pass), __STAT__)

    active_instances = n_inst_act

!   !Fill the instances list with active instances first
!   !getAERO_ depends on the active instances being filled first.
    call ESMF_ConfigFindLabel (myCF, 'ACTIVE_INSTANCES_'//trim(aerosol)//':', __RC__)
    do i = 1, n_inst_act
        call ESMF_ConfigGetAttribute (myCF, instances(i), __RC__)
    end do
!   !Now fill instances list with passive instances
    call ESMF_ConfigFindLabel (myCF, 'PASSIVE_INSTANCES_'//trim(aerosol)//':', __RC__)
    do i = n_inst_act+1, n_inst_act+n_inst_pass
        call ESMF_ConfigGetAttribute (myCF, instances(i), __RC__)
    end do

    RETURN_(ESMF_SUCCESS)

  end subroutine getInstances_


!====================================================================================
  
  subroutine getAERO_ (self, AEROlist, rc)

!   Description: Fills the AEROlist with all active instances of each GOCART2G child.
!                Active instances are defined in 'ACTIVE_INSTANCES_XX:' in GOCART2G_GridComp.rc 
!                If additional GOCART2g children are added, this subroutine will need to be updated.


    implicit none

    type (GOCART_State), pointer,                   intent(in   )     :: self
    character (len=ESMF_MAXSTR), allocatable,       intent(inout)     :: AEROlist(:)  !names of active instances 
    integer,                                        intent(  out)     :: rc

!   locals
    integer                                              :: STATUS
    character (len=ESMF_MAXSTR)                          :: Iam = 'GOCART2G::getAERO_'
    integer                                              :: i, ind, tot_active_inst
!--------------------------------------------------------------------------------------

!   Begin...

!   !Get total number of active instances
    tot_active_inst = self%active_instances_DU + self%active_instances_SS + self%active_instances_SU + &
                      self%active_instances_BC + self%active_instances_OC + self%active_instances_NI

    allocate (AEROlist(tot_active_inst), __STAT__)

    ind = 1

!   !Fill the AEROlist with names of the active instances        
    do i = 1, self%active_instances_DU
        AEROlist(ind) = self%instances_DU(i)
        ind = ind + 1
    end do

    do i = 1, self%active_instances_SS
        AEROlist(ind) = self%instances_SS(i)
        ind = ind + 1
    end do

    do i = 1, self%active_instances_SU
        AEROlist(ind) = self%instances_SU(i)
        ind = ind + 1
    end do

    do i = 1, self%active_instances_BC
        AEROlist(ind) = self%instances_BC(i)
        ind = ind + 1
    end do

    do i = 1, self%active_instances_OC
        AEROlist(ind) = self%instances_OC(i)
        ind = ind + 1
    end do

    do i = 1, self%active_instances_NI
        AEROlist(ind) = self%instances_NI(i)
        ind = ind + 1
    end do

    RETURN_(ESMF_SUCCESS)

  end subroutine getAERO_


!====================================================================================

  subroutine run_aerosol_optics(state, rc)

    implicit none

!   !ARGUMENTS:
    type (ESMF_State)                                :: state
    integer,            intent(out)                  :: rc

!   !Local
    real, dimension(:,:,:), pointer                  :: ple
    real, dimension(:,:,:), pointer                  :: rh
    real, dimension(:,:,:), pointer                  :: var

    character (len=ESMF_MAXSTR)                      :: fld_name

    real, dimension(:,:,:),pointer                   :: ext_, ssa_, asy_      ! (lon:,lat:,lev:,band:)
    real, dimension(:,:,:), allocatable              :: ext,  ssa,  asy       ! (lon:,lat:,lev:,band:)

    integer                                          :: i 
    integer                                          :: i1, j1, i2, j2, km
    integer                                          :: band
    integer, parameter                               :: n_bands = 1

    character (len=ESMF_MAXSTR), allocatable         :: AEROlist(:)
    type (ESMF_State)                                :: child_state
    real, pointer,     dimension(:,:,:)              :: AS_PTR_3D


    __Iam__('GOCART2G::run_aerosol_optics')


!   Begin... 

!   Radiation band
!   --------------
    call ESMF_AttributeGet(state, name='band_for_aerosol_optics', value=band, __RC__)

!   Relative humidity
!   -----------------
    call ESMF_AttributeGet(state, name='relative_humidity_for_aerosol_optics', value=fld_name, __RC__)
    call MAPL_GetPointer(state, RH, trim(fld_name), __RC__)

!   Pressure at layer edges 
!   ------------------------
    call ESMF_AttributeGet(state, name='air_pressure_for_aerosol_optics', value=fld_name, __RC__)
    call MAPL_GetPointer(state, PLE, trim(fld_name), __RC__)

    i1 = lbound(ple, 1); i2 = ubound(ple, 1)
    j1 = lbound(ple, 2); j2 = ubound(ple, 2)
                         km = ubound(ple, 3)


    allocate(ext(i1:i2,j1:j2,km),  &
             ssa(i1:i2,j1:j2,km),  &
             asy(i1:i2,j1:j2,km), __STAT__)


    call ESMF_AttributeGet (state, name='active_aerosol_instances', itemCount=i, __RC__)
    allocate (AEROlist(i), __STAT__)
    call ESMF_AttributeGet (state, name='active_aerosol_instances', valueList=AEROlist, __RC__)

    ext = 0.0d0
    ssa = 0.0d0
    asy = 0.0d0


   do i = 1, size(AEROlist)
        call ESMF_StateGet(state, trim(AEROlist(i))//'_AERO', child_state, __RC__)

!       ! set RH for aerosol optics
        call ESMF_AttributeGet(child_state, name='relative_humidity_for_aerosol_optics', value=fld_name, __RC__)

        if (fld_name /= '') then
            call MAPL_GetPointer(child_state, AS_PTR_3D, trim(fld_name), __RC__)
            AS_PTR_3D = RH
        end if

!       ! set PLE for aerosol optics
        call ESMF_AttributeGet(child_state, name='air_pressure_for_aerosol_optics', value=fld_name, __RC__)

        if (fld_name /= '') then
            call MAPL_GetPointer(child_state, AS_PTR_3D, trim(fld_name), __RC__)
            AS_PTR_3D = PLE
        end if

        call ESMF_AttributeSet(child_state, name='band_for_aerosol_optics', value=band, __RC__)

!       ! execute the aero provider's optics method
        call ESMF_MethodExecute(child_state, label="aerosol_optics", __RC__)

!       ! EXT from AERO_PROVIDER
        call ESMF_AttributeGet(child_state, name='extinction_in_air_due_to_ambient_aerosol', value=fld_name, __RC__)
        if (fld_name /= '') then
            call MAPL_GetPointer(child_state, ext_, trim(fld_name), __RC__)
        end if

!       ! SSA from AERO_PROVIDER
        call ESMF_AttributeGet(child_state, name='single_scattering_albedo_of_ambient_aerosol', value=fld_name, __RC__)
        if (fld_name /= '') then
            call MAPL_GetPointer(child_state, ssa_, trim(fld_name), __RC__)
        end if

!       ! ASY from AERO_PROVIDER
        call ESMF_AttributeGet(child_state, name='asymmetry_parameter_of_ambient_aerosol', value=fld_name, __RC__)
        if (fld_name /= '') then
            call MAPL_GetPointer(child_state, asy_, trim(fld_name), __RC__)
        end if

        ext = ext + ext_
        ssa = ssa + ssa_
        asy = asy + asy_
    end do



    call ESMF_AttributeGet(state, name='extinction_in_air_due_to_ambient_aerosol', value=fld_name, __RC__)
    if (fld_name /= '') then
        call MAPL_GetPointer(state, var, trim(fld_name), __RC__)
        var = ext(:,:,:)
    end if

    call ESMF_AttributeGet(state, name='single_scattering_albedo_of_ambient_aerosol', value=fld_name, __RC__)
    if (fld_name /= '') then
        call MAPL_GetPointer(state, var, trim(fld_name), __RC__)
        var = ssa(:,:,:)
    end if

    call ESMF_AttributeGet(state, name='asymmetry_parameter_of_ambient_aerosol', value=fld_name, __RC__)
    if (fld_name /= '') then
        call MAPL_GetPointer(state, var, trim(fld_name), __RC__)
        var = asy(:,:,:)
    end if

    deallocate(ext, ssa, asy, __STAT__)



   RETURN_(ESMF_SUCCESS)

  end subroutine run_aerosol_optics




end module GOCART2G_GridCompMod






