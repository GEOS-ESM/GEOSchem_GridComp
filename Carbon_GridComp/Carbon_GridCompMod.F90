#include "MAPL_Generic.h"

!=============================================================================
!BOP

! !MODULE: Carbon_GridCompMod - Carbon Component. Refactored from Legacy GOCART.

! !INTERFACE:

module Carbon_GridCompMod

! !USES:

   use ESMF
   use MAPL

! !Establish the Childen's SetServices
 !-----------------------------------
!   use CO_GridCompMod,     only   : CO_setServices   => SetServices
!   use CO2_GridCompMod,    only   : CO2_setServices  => SetServices
!   use OH_GridCompMod,     only   : OH_setServices   => SetServices

   implicit none
   private

! !PUBLIC MEMBER FUNCTIONS:
   public  SetServices


  ! Private State
  type :: Instance
     integer :: id = -1
!     logical :: is_active
     character(:), allocatable :: name
  end type Instance

  type Constituent
     type(Instance), allocatable :: instances(:)
 !    integer :: n_active
  end type Constituent

  type Carbon_State
     private
     type(Constituent) :: CO
     type(Constituent) :: CO2
     type(Constituent) :: OH
  end type Carbon_State

  type wrap_
     type (Carbon_State), pointer  :: PTR => null()
  end type wrap_


! !DESCRIPTION:
!
!  Carbon is a gridded component refactored from the GOCART model and 
!  includes CO, CO2, and OH.

!
! !REVISION HISTORY:
!
!  25feb2005   da Silva   First crack.
!  19jul2006   da Silva   First separate GOCART component.
!  15June2021  E.Sherman  First attempt at refactoring. 

!
!EOP
!============================================================================

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
!============================================================================
!
!   Locals
    character (len=ESMF_MAXSTR)                   :: comp_name
    type (ESMF_Config)                            :: myCF      ! Carbon_GridComp.rc
    type (ESMF_Config)                            :: cf        ! universal config
    type (GOCART_State), pointer                  :: self
    type (wrap_)                                  :: wrap

    __Iam__('SetServices')

!****************************************************************************
! Begin...

!   Get my name and set-up traceback handle
!   ---------------------------------------
    call ESMF_GridCompGet (GC, name=comp_name, config=cf, __RC__)
    Iam = trim(comp_name)//'::'//'SetServices'

!   Wrap internal state for storing in GC
!   -------------------------------------
    allocate (self, __STAT__)
    wrap%ptr => self

!   Set the Initialize, Run, Finalize entry points
!   ------------------------------------------------
    call MAPL_GridCompSetEntryPoint (GC, ESMF_Method_Initialize,  Initialize,  __RC__)
    call MAPL_GridCompSetEntryPoint (GC, ESMF_Method_Run,  Run1, __RC__)
    call MAPL_GridCompSetEntryPoint (GC, ESMF_Method_Run,  Run2, __RC__)

!   Store internal state in GC
!   --------------------------
    call ESMF_UserCompSetInternalState (GC, 'Carbon_State', wrap, STATUS)
    VERIFY_(STATUS)

!   Get instances to determine what grid comps will be instantiated
!   ---------------------------------------------------------------
    call getInstances_('CO',  myCF, species=self%CO, __RC__)
    call getInstances_('CO2', myCF, species=self%CO2, __RC__)
    call getInstances_('OH',  myCF, species=self%OH, __RC__)

!   Create children's gridded components and invoke their SetServices
!   -----------------------------------------------------------------
    call createInstances_(self, GC, __RC__)

if(mapl_am_i_root()) print*,'CARBON COMP SET SERVICES'

!   Set generic services
!   ----------------------------------
    call MAPL_GenericSetServices (GC, __RC__)

    RETURN_(ESMF_SUCCESS)

  end subroutine SetServices







!===============================================================================

  subroutine getInstances_ (aerosol, myCF, species, rc)

!   Description: Fills the GOCART_State (aka, self%instance_XX) with user
!                defined instances from the GOCART2G_GridComp.rc.

    implicit none

    character (len=*),                intent(in   )  :: carbonType
    type (ESMF_Config),               intent(inout)  :: myCF
    type(Constituent),                intent(inout)  :: species
    integer,                          intent(  out)  :: rc


!   locals
    integer                                          :: i
    integer                                          :: n_instances
    character (len=ESMF_MAXSTR)                      :: inst_name

    __Iam__('Carbon::getInstances_')

!--------------------------------------------------------------------------------------

!   Begin...
    n_instances  = ESMF_ConfigGetLen (myCF, label='INSTANCES_'//trim(carbonType)//':', __RC__)
    allocate (species%instances(n_instances), __STAT__)

!   !Fill the instances list with active instances first
    call ESMF_ConfigFindLabel (myCF, 'INSTANCES_'//trim(carbonType)//':', __RC__)
    do i = 1, n_instances
       call ESMF_ConfigGetAttribute (myCF, inst_name, __RC__)
       species%instances(i)%name = inst_name
    end do
    species%n_active = n_active



    RETURN_(ESMF_SUCCESS)

  end subroutine getInstances_


!====================================================================================
  subroutine createInstances_ (self, GC, rc)

!   Description: Creates GOCART2G children. Active instances must be created first. If
!     additional GOCART2G children are added, this subroutine will need to be updated.

    implicit none

    type (GOCART_State), pointer,            intent(in   )     :: self
    type (ESMF_GridComp),                    intent(inout)     :: GC
    integer,                                 intent(  out)     :: rc

    ! locals
    integer                                                    :: i

    __Iam__('GOCART2G::createInstances_')

!-----------------------------------------------------------------------------------
!   Begin...

!   Active instances must be created first! This ordering is necessary for
!   filing the AERO states that are passed to radiation.
!   This is achieved by arranging the names of the active instances first.

    call addChildren__ (gc, self%CO,  setServices=CO_setServices, __RC__)
    call addChildren__ (gc, self%CO2, setServices=CO2_setServices, __RC__)
    call addChildren__ (gc, self%OH,  setServices=OH_setServices, __RC__)

    RETURN_(ESMF_SUCCESS)

    contains

        subroutine addChildren__ (gc, species, setServices, rc)

          type (ESMF_GridComp),            intent(inout)     :: gc
          type(Constituent),               intent(inout)     :: species
          external                                           :: setServices
          integer,                         intent(  out)     :: rc

          ! local
          integer  :: n

          __Iam__('GOCART2G::createInstances_::addChildren__')

          n=size(species%instances)

          do i = 1, n
             species%instances(i)%id = MAPL_AddChild(gc, name=species%instances(i)%name, SS=SetServices, __RC__)
          end do

        RETURN_(ESMF_SUCCESS)

     end subroutine addChildren__

  end subroutine createInstances_




end module Carbon_GridCompMod

