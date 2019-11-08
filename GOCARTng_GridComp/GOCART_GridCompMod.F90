#include "MAPL_Generic.h"

!=============================================================================
!BOP

! !MODULE: GOCARTng_GridCompMod - The GOCART Next Generation Aerosol Grid Component

! !INTERFACE:

module GOCARTng_GridCompMod

! !USES:

   USE ESMF
   USE MAPL_Mod

! !Establish the Childen's SetServices
! !-----------------------------------
   use DUng_GridCompMod,    only   : DUngSetServices  => SetServices
   use SSng_GridCompMod,    only   : SSngSetServices  => SetServices


   implicit none
   private

! !PUBLIC MEMBER FUNCTIONS:
   PUBLIC  SetServices
!   PUBLIC  Initialize

! Private State
  type GOCARTng_State
     private
     character (len=ESMF_MAXSTR), pointer    :: instances_DU(:) => null()
     character (len=ESMF_MAXSTR), pointer    :: instances_SS(:) => null()
     character (len=ESMF_MAXSTR), pointer    :: instances_SU(:) => null()
     character (len=ESMF_MAXSTR), pointer    :: instances_BC(:) => null()
     character (len=ESMF_MAXSTR), pointer    :: instances_OC(:) => null()
     character (len=ESMF_MAXSTR), pointer    :: instances_NI(:) => null()
!     character (len=ESMF_MAXSTR), pointer    :: active_instances_DU = 0
!     character (len=ESMF_MAXSTR), pointer    :: active_instances_SS = 0
!     character (len=ESMF_MAXSTR), pointer    :: active_instances_SU = 0
!     character (len=ESMF_MAXSTR), pointer    :: active_instances_BC = 0
!     character (len=ESMF_MAXSTR), pointer    :: active_instances_OC = 0
!     character (len=ESMF_MAXSTR), pointer    :: active_instances_NI = 0
     integer                                 :: active_instances_DU = 0
     integer                                 :: active_instances_SS = 0
     integer                                 :: active_instances_SU = 0
     integer                                 :: active_instances_BC = 0
     integer                                 :: active_instances_OC = 0
     integer                                 :: active_instances_NI = 0
  end type GOCARTng_State

  type wrap_
     type (GOCARTng_State), pointer     :: PTR => null()
  end type wrap_


! !DESCRIPTION:
!
!   {\tt GOCART} is a gridded component from the GOCART model and includes 
!  dust, sea salt, sulfates, organic and black carbon. In addition, we
!  also include closely related components for CO and CO2 with relatively
!  simple parameterization of the chemical processes, but sharing
!  consistent emissions with the aerosols.
!
!
! !REVISION HISTORY:
!
!  25feb2005  da Silva   First crack.
!  19jul2006  da Silva   First separate GOCART component.
!  14Oct2019  E.Sherman  First attempt at refactoring for ESMF compatibility
!
!EOP
!-------------------------------------------------------------------------


  integer ::     DUng = -1
  integer ::     SSng = -1

contains

!BOP

! !IROUTINE: SetServices -- Sets ESMF services for this component

! !INTERFACE:

  subroutine SetServices (GC, RC)

! !ARGUMENTS:

    type (ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
    integer, optional                   :: RC  ! return code

! !DESCRIPTION: This version uses the {\tt GEOS\_GenericSetServices}, which sets
!   the Initialize and Finalize services to generic versions. It also
!   allocates our instance of a generic state and puts it in the 
!   gridded component (GC). Here we only set the two-stage run method and
!   declare the data services.
! \newline

! !REVISION HISTORY: 
!  14oct2019  E.Sherman  First attempt at refactoring for ESMF compatibility


!EOP

!****************************************************************************
!
! ErrLog Variables

    character (len=ESMF_MAXSTR)                   :: IAm
    integer                                       :: STATUS
    character (len=ESMF_MAXSTR)                   :: COMP_NAME

! Locals

    type (ESMF_Config)                            :: myCF
    type (GOCARTng_State), pointer                :: self
    type (wrap_)                                  :: wrap

    integer                                       :: i,nq

!****************************************************************************

! Begin...

!   Get my name and set-up traceback handle
!   ---------------------------------------
    call ESMF_GridCompGet (GC, NAME=COMP_NAME, __RC__)
    Iam = trim(COMP_NAME) // '::' //  'SetServices'


    if (mapl_am_i_root()) print*,'GOCARTng SetServices BEGIN'

!   Wrap internal state for storing in GC
!   -------------------------------------
    allocate (self, __STAT__)
    wrap%ptr => self

!   Set the Initialize, Run, Finalize entry points
!   ------------------------------------------------
    call MAPL_GridCompSetEntryPoint (GC, ESMF_METHOD_INITIALIZE,  Initialize,  __RC__)
    call MAPL_GridCompSetEntryPoint (GC, ESMF_METHOD_RUN,  Run1, __RC__)

!   Store internal state in GC
!   --------------------------
    call ESMF_UserCompSetInternalState (GC, 'GOCARTng_state', wrap, STATUS)
    VERIFY_(STATUS)

!   What children will be born?
!   ----------------------------
    myCF = ESMF_ConfigCreate (__RC__)
    call ESMF_ConfigLoadFile (myCF, 'GOCARTng_GridComp.rc', __RC__)

    call getInstances_('DU', myCF, instances=self%instances_DU, active_instances=self%active_instances_DU, __RC__)
    call getInstances_('SS', myCF, instances=self%instances_SS, active_instances=self%active_instances_SS, __RC__)
    call getInstances_('SU', myCF, instances=self%instances_SU, active_instances=self%active_instances_SU, __RC__)
    call getInstances_('BC', myCF, instances=self%instances_BC, active_instances=self%active_instances_BC, __RC__)
    call getInstances_('OC', myCF, instances=self%instances_OC, active_instances=self%active_instances_OC, __RC__)
    call getInstances_('NI', myCF, instances=self%instances_NI, active_instances=self%active_instances_NI, __RC__)

    call ESMF_ConfigDestroy(myCF, __RC__)

!   Create children`s gridded components and invoke their SetServices
!   -----------------------------------------------------------------
    do i = 1, size(self%instances_DU)
        DUng = MAPL_AddChild(GC, NAME=trim(self%instances_DU(i)), SS=DUngSetServices, __RC__)
    end do

    do i = 1, size(self%instances_SS)
        SSng = MAPL_AddChild(GC, NAME=trim(self%instances_SS(i)), SS=SSngSetServices, __RC__)
    end do


!   Define EXPORT states
!   This state is needed by radiation - It will contain 
!   aerosols and aerosol optics
!   --------------------------------------------------------
    call MAPL_AddExportSpec(GC,                       &
       SHORT_NAME = 'AEROng',                           &
       LONG_NAME  = 'aerosol_mass_mixing_ratios_ng',  &
       UNITS      = 'kg kg-1',                        &
       DIMS       = MAPL_DimsHorzVert,                &
       VLOCATION  = MAPL_VLocationCenter,             &
       DATATYPE   = MAPL_StateItem, __RC__)

!   This state is needed by MOIST - It will contain
!   aerosols
!   --------------------------------------------------------
    call MAPL_AddExportSpec(GC,                       &
       SHORT_NAME = 'AERO_ACI',                       &
       LONG_NAME  = 'aerosol_cloud_interaction_ng',   &
       UNITS      = 'kg kg-1',                        &
       DIMS       = MAPL_DimsHorzVert,                &
       VLOCATION  = MAPL_VLocationCenter,             &
       DATATYPE   = MAPL_StateItem, __RC__)

!   This bundle is needed by surface for snow albedo modification
!   by aerosol settling and deposition
!   DEVELOPMENT NOTE - Change to StateItem in future
!   --------------------------------------------------------
    call MAPL_AddExportSpec(GC,                       &
       SHORT_NAME = 'AERO_DP',                        &
       LONG_NAME  = 'aerosol_deposition_ng',          &
       UNITS      = 'kg m-2 s-1',                     &
       DIMS       = MAPL_DimsHorzOnly,                &
       DATATYPE   = MAPL_BundleItem, __RC__)




!   Set generic services
!   ----------------------------------
    call MAPL_GenericSetServices (GC, __RC__)



    if (mapl_am_i_root()) print*,'GOCARTng SetServices END'

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

! !DESCRIPTION:  This initializes the GOCART Grid Component. It primary creates
!                its exports and births its children.

! !REVISION HISTORY: 
! 14oct2019   E.Sherman  First attempt at refactoring


!EOP

!****************************************************************************
! ErrLog Variables

    character (len=ESMF_MAXSTR)          :: IAm
    integer                              :: STATUS
    character (len=ESMF_MAXSTR)          :: COMP_NAME

! Local derived type aliases
    type (MAPL_MetaComp),       pointer  :: MAPL
    type (ESMF_GridComp),       pointer  :: GCS(:)
    type (ESMF_State),          pointer  :: GEX(:)
    type (ESMF_State)                    :: INTERNAL
    type (ESMF_Grid)                     :: grid

    type (ESMF_State)                    :: AERO, AERO_ACI
    type (ESMF_FieldBundle)              :: AERO_DP, child_bundle
    type (ESMF_State)                    :: child_state
    type (ESMF_Field), allocatable       :: fieldList(:)
    type (ESMF_Field)                    :: field
    character (len=ESMF_MAXSTR)          :: field_name

    type (GOCARTng_state),      pointer  :: self
    type (wrap_)                         :: wrap

    character (len=ESMF_MAXSTR)          :: CHILD_NAME
    character (len=ESMF_MAXSTR), allocatable :: AEROlist(:)

    integer                              :: i, j, k, fieldCount

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
    call ESMF_UserCompGetInternalState (GC, 'GOCARTng_state', wrap, STATUS)
    VERIFY_(STATUS)
    self => wrap%ptr
!   if ( present(state) ) then
!        state => wrap%ptr
!   end if

!   Get children and their export states from my generic state
!   -----------------------------------------------------------
    call MAPL_Get (MAPL, GCS=GCS, GEX=GEX, __RC__ )

!   Fill AERO, AERO_ACI, and AERO_DP with the analogous children's states
!   ----------------------------------------------------------------------
    call ESMF_StateGet (EXPORT, 'AEROng'  , AERO     , __RC__)
    call ESMF_StateGet (EXPORT, 'AERO_ACI', AERO_ACI , __RC__)
    call ESMF_StateGet (EXPORT, 'AERO_DP' , AERO_DP  , __RC__)

!   Get list of AERO_PROVIDERs 
!   ---------------------------
    call getAERO_ (self, AEROlist, __RC__)

if (mapl_am_i_root()) print*,'GOCARTng AEROlist = ', AEROlist

    call ESMF_AttributeSet(AERO, name='active_aerosol_instances', valueList=AEROlist, itemCount=size(AEROlist), __RC__)


!   Add AERO_PROVIDER's states to GOCART's AERO states
!   ---------------------------------------------------
    do i = 1, size(GCS)
        call ESMF_GridCompGet (GCS(i), NAME=CHILD_NAME, __RC__ )
        do j = 1, size(AEROlist)
            if (trim(AEROlist(j)) == trim(CHILD_NAME)) then
                call ESMF_StateGet (GEX(i), trim(CHILD_NAME)//'_AERO', child_state, __RC__)
                call ESMF_StateAdd (AERO, (/child_state/), __RC__)

                call ESMF_StateGet (GEX(i), trim(CHILD_NAME)//'_AERO_ACI', child_state, __RC__)
                call ESMF_StateAdd (AERO_ACI, (/child_state/), __RC__)

                call ESMF_StateGet (GEX(i), trim(CHILD_NAME)//'_AERO_DP', child_bundle, __RC__)
                call ESMF_FieldBundleGet (child_bundle, fieldCount=fieldCount, __RC__)
                allocate (fieldList(FieldCount), __STAT__)
                call ESMF_FieldBundleGet (child_bundle, fieldList=fieldList, __RC__)
                call ESMF_FieldBundleAdd (AERO_DP, fieldList, __RC__)
                deallocate(fieldList, __STAT__)
            end if
        end do     
    end do


!   Verify that childen's states are properly added - for testing to be deleted
    if(mapl_am_i_root()) print*,'GOCARTng AERO print state = '
    if(mapl_am_i_root()) then
        call esmf_stateprint(AERO, __RC__)
    end if

!   Verify that childen's states are properly added - for testing to be deleted
!    if(mapl_am_i_root()) print*,'GOCARTng AERO_ACI print state = '
!    if(mapl_am_i_root()) then
!        call esmf_stateprint(AERO_ACI, __RC__)
!    end if

!   Verify that childen's states are properly added - for testing to be deleted
    if(mapl_am_i_root()) print*,'GOCARTng AERO_DP print bundle = '
    if(mapl_am_i_root()) then
        call esmf_fieldbundleprint(AERO_DP, __RC__)
    end if
  
!   Verify contents of the sea salt state - for testing to be deleted
!  if(mapl_am_i_root()) print*,'last child_state print state = '
!  if(mapl_am_i_root()) then
!      call esmf_stateprint(child_state, __RC__)
!  end if

!if (mapl_am_i_root()) print*,'GOCARTng Initialize END'


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

if (mapl_am_i_root()) print*,'GOCARTng Run1 BEGIN'


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


if (mapl_am_i_root()) print*,'GOCARTng Run1 END'

    RETURN_(ESMF_SUCCESS)

  end subroutine Run1

!===============================================================================

  subroutine getInstances_ (aerosol, myCF, instances, active_instances, rc)

!   Description: Fills the GOCARTng_State (aka, self%instance_XX) with user
!                defined instances from the GOCARTng_GridComp.rc.

    implicit none

    character (len=*),                intent(in   )  :: aerosol
    type (ESMF_Config),               intent(inout)  :: myCF
    integer,                          intent(  out)  :: rc
    character (len=ESMF_MAXSTR), pointer             :: instances(:)
    integer,                          intent(inout)  :: active_instances

!   locals
    integer                                          :: STATUS
    character (len=ESMF_MAXSTR)                      :: Iam = 'GOCARTng::getInstances_'

    integer                                          :: i, n_inst_act, n_inst_pass
    character (len=ESMF_MAXSTR)                      :: inst_name
!--------------------------------------------------------------------------------------

!   Begin...
    n_inst_act  = ESMF_ConfigGetLen (myCF, label='ACTIVE_INSTANCES_'//trim(aerosol)//':', __RC__)
    n_inst_pass = ESMF_ConfigGetLen (myCF, label='PASSIVE_INSTANCES_'//trim(aerosol)//':', __RC__)
    allocate (instances(n_inst_act + n_inst_pass), __STAT__)

    active_instances = n_inst_act

    call ESMF_ConfigFindLabel (myCF, 'ACTIVE_INSTANCES_'//trim(aerosol)//':', __RC__)
    do i = 1, n_inst_act
        call ESMF_ConfigGetAttribute (myCF, instances(i), __RC__)
    end do

    call ESMF_ConfigFindLabel (myCF, 'PASSIVE_INSTANCES_'//trim(aerosol)//':', __RC__)
    do i = n_inst_act+1, n_inst_act+n_inst_pass
        call ESMF_ConfigGetAttribute (myCF, instances(i), __RC__)
    end do

    RETURN_(ESMF_SUCCESS)

  end subroutine getInstances_


!====================================================================================
  
  subroutine getAERO_ (self, AEROlist, rc)

!   Description: Fills the AEROlist with the first instance of each GOCARTng child.
!                The first defined instance (ACTIVE_INSTANCES_XX:) in GOCARTng_GridComp.rc 
!                is the AERO_PROVIDER.
!                If additional aerosol gridded components are added, this subroutine
!                will need to be updated.
  
    implicit none

    type (GOCARTng_state), pointer,                 intent(in   )     :: self
    character (len=ESMF_MAXSTR), allocatable,       intent(inout)     :: AEROlist(:) 
    integer,                                        intent(  out)     :: rc

!   locals
    integer                                              :: STATUS
    character (len=ESMF_MAXSTR)                          :: Iam = 'GOCARTng::getAERO_'
    integer                                              :: i, ind, tot_active_inst
!--------------------------------------------------------------------------------------

!   Begin...
!    AEROlist(1) = trim(self%instances_DU(1))
!    AEROlist(2) = trim(self%instances_SS(1))
!    AEROlist(3) = trim(self%instances_SU(1))
!    AEROlist(4) = trim(self%instances_BC(1))
!    AEROlist(5) = trim(self%instances_OC(1))
!    AEROlist(6) = trim(self%instances_NI(1))

    tot_active_inst = self%active_instances_DU + self%active_instances_SS + self%active_instances_SU + &
                      self%active_instances_BC + self%active_instances_OC + self%active_instances_NI

    allocate (AEROlist(tot_active_inst), __STAT__)

    ind = 1
        
    do i = 1, self%active_instances_DU
        AEROlist(ind) = self%instances_DU(i)
        ind = ind + 1
    end do

    do i = 1, self%active_instances_SS
        AEROlist(ind) = self%instances_SS(i)
        ind = ind + 1
    end do






    RETURN_(ESMF_SUCCESS)

  end subroutine getAERO_


!====================================================================================

  subroutine get_optics (state, rc)

    implicit none

!   !ARGUMENTS:
    type (ESMF_State)                                :: state
    integer,            intent(out)                  :: rc

!   !Local
    type (ESMF_State)                                :: child_state
    integer                                          :: i
!   character (len=ESMF_MAXSTR), allocatable         :: AEROlist(:)



    __Iam__('GOCARTng::getOptics_')



!  call ESMF_AttributeGet (state, name='active_aerosol_instances', itemCount=i, __RC__)
!  allocate(AEROlist(i), __STAT__)
!  call ESMF_AttributeGet (state, name='active_aerosol_instances', valueList=AEROlisttest, __RC__)

!f (mapl_am_i_root()) print*,'AEROlisttest = ', AEROlisttest

!   execute childrens' aerosol_optics method and sum output  to get totals
!   do i = 1 , size(AEROlist)
!       call ESMF_StateGet(state, AEROlist(i), child_state, __RC__)

!     call ESMF_AttributeGet(child_state, name='relative_humidity_for_aerosol_optics', value=AS_FIELD_NAME, RC=STATUS)
!     VERIFY_(STATUS)

!     if (AS_FIELD_NAME /= '') then
!        call MAPL_GetPointer(AEROng, AS_PTR_3D, trim(AS_FIELD_NAME), RC=STATUS)
!        VERIFY_(STATUS)

!        AS_PTR_3D = RH
!     end if








!       call ESMF_MethodExecute(AEROlist(i), label="aerosol_optics", __RC__)
!    end do



  end subroutine get_optics




end module GOCARTng_GridCompMod






