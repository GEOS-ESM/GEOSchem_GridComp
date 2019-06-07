#include "MAPL_Exceptions.h"
#include "MAPL_Generic.h"


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: MAML_DryRemovalMod - Gravitational sedimentation/settling and dry
!                               deposition of aerosol particles and gases
!
! !INTERFACE:
!
   module MAML_DryRemovalMod
!
! !USES:
!
   use MAPL_Mod

   use MAPL_ConstantsMod


   implicit NONE
   private

!
! !PUBLIC MEMBER FUNCTIONS:

   public MAML_DryRemoval

!
! !PUBLIC PARAMETERS:

!
! !PRIVATE PARAMETERS:

   real, private, parameter :: g_E    = MAPL_GRAV     ! standard gravity,        'm s-2'

!
! !DESCRIPTION: 
!
!  {\tt MAML\_SettlingMod} provides a collection of methods for 
!  modeling graviational sedimentation/settling of aerosol particles
!
!
! !REVISION HISTORY:
!
!  29Nov2011  A. Darmenov  Initial version
!
!EOP
!-------------------------------------------------------------------------


   interface MAML_DryRemoval
       module procedure MAML_DryRemovalAerosol
       module procedure MAML_DryRemovalGas
   end interface MAML_DryRemoval

   contains

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: MAML_DryRemovalTendencySolverAerosol --- 
!
! !INTERFACE:

   subroutine MAML_DryRemovalTendencySolverAerosol(dqdt, q, delp, dz, v_t, v_d)

! !USES:

   implicit None


! !INPUT/OUTPUT PARAMETERS:
   real, dimension(:), intent(inout) :: dqdt    ! dq/dt - mixing ratio tendency due to
                                                !         gravitational sedimentation

! !INPUT PARAMETERS:

   real, dimension(:), intent(in)    :: q       ! mixing ratio,                 'kg kg-1' or '# kg-1'
   real, dimension(:), intent(in)    :: delp    ! pressure thickness of levels, 'Pa'
   real, dimension(:), intent(in)    :: dz      ! thickness of levels,          'm'
   real, dimension(:), intent(in)    :: v_t     ! settling velocity,            'm s-1'
   real, intent(in)                  :: v_d     ! deposition velocity,          'm s-1'


! !OUTPUT PARAMETERS:


! !DESCRIPTION: Calculates mixing ratio tendency (dq/dt) due to 
!               gravitational sedimentation/settling and dry deposition.
!
! !REVISION HISTORY:
!
!  19Nov2011  A. Darmenov   First crack.
!  20Dec2011  A. Darmenov   Gravitational sedimentation and dry deposition 
!                           are done in parallel
!
!EOP
!-------------------------------------------------------------------------
   
                   __Iam__('MAML_DryRemovalTendencySolverAerosol')

   real, allocatable, dimension(:) :: v
   integer :: k1, k2, km
   integer :: rc

   k1 = lbound(q, 1)
   km = ubound(q, 1)

   allocate(v(k1:km), __STAT__)
   
   v(k1:km-1) = v_t(k1:km-1)
   v(km) = v_d
   
   k2   = k1 + 1
   dqdt = 0.0

   ! zero flux in from top of the atmosphere (k=k1)
   dqdt(k1) = 0.0 - q(k1) * v(k1) / dz(k1)

   ! levels k2:km --- flux in from the level above and flux out into the level below
   dqdt(k2:km) = (q(k1:km-1) * (v(k1:km-1) / dz(k1:km-1))) * (delp(k1:km-1)/delp(k2:km)) - &
                 (q(k2:km  ) * (v(k2:km  ) / dz(k2:km  )))

   deallocate(v, __STAT__)

   end subroutine MAML_DryRemovalTendencySolverAerosol


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: MAML_DryRemovalTendencySolverGas --- 
!
! !INTERFACE:

   subroutine MAML_DryRemovalTendencySolverGas(dqdt, q, dz, v_d, dt)

! !USES:

   implicit None


! !INPUT/OUTPUT PARAMETERS:
   real, intent(inout) :: dqdt    ! dq/dt - mixing ratio tendency due to
                                  !         gravitational sedimentation

! !INPUT PARAMETERS:

   real, intent(in)    :: q       ! mixing ratio,                    'kg kg-1' or '# kg-1'
   real, intent(in)    :: dz      ! thickness of the surface level,  'm'
   real, intent(in)    :: v_d     ! deposition velocity,             'm s-1'
   real, intent(in)    :: dt      ! time step,                       's'


! !OUTPUT PARAMETERS:


! !DESCRIPTION: Calculates mixing ratio tendency (dq/dt) due to dry deposition.
!               Follows Kerkweg et al., 2006 approach to avoid depletion (q < 0) 
!               of the species by introducing effective deposition velocity.
!
! !REVISION HISTORY:
!
!  20Dec2011  A. Darmenov   First crack.
!
!EOP
!-------------------------------------------------------------------------
   
                   __Iam__('MAML_DryRemovalTendencySolverGas')

  
   real :: v_eff    ! effective deposition velocity: derived by 
                    ! integrating the equation:
                    !
                    !         dq/dt = - q * (v_d / dz)


   v_eff = dz/dt * (1 - exp(-v_d * dt/dz))

   dqdt = - q * v_eff / dz

   end subroutine MAML_DryRemovalTendencySolverGas


   
!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: MAML_DryRemovalAerosol --- 
!
! !INTERFACE:

   subroutine MAML_DryRemovalAerosol(q, delp, dz, v_t, v_d, dt, flux)

! !USES:

   implicit None


! !INPUT/OUTPUT PARAMETERS:
   real, dimension(:), intent(inout) :: q         ! mixing ratio,                 'kg kg-1' or '# kg-1'
   real, optional,     intent(inout) :: flux      ! deposition flux               'kg m-2 s-1' or '# m-2 s-1'

! !INPUT PARAMETERS:
   real, dimension(:), intent(in)    :: delp      ! pressure thickness of levels, 'Pa'
   real, dimension(:), intent(in)    :: dz        ! thickness of levels,          'm'
   real, dimension(:), intent(in)    :: v_t       ! settling velocity,            'm s-1'
   real, intent(in)                  :: v_d       ! deposition velocity,          'm s-1'
   real, intent(in)                  :: dt        ! model time step,              's'

! !OUTPUT PARAMETERS:


! !DESCRIPTION: Calculates the changes in the mixing ratio due to 
!               gravitational sedimentation/settling.
!
! !REVISION HISTORY:
!
!  6Nov2011  A. Darmenov
!
!EOP
!-------------------------------------------------------------------------
   
                  __Iam__('MAML_DryRemovalAerosol')

   real, pointer, dimension(:) :: dqdt            ! dq/dt - mixing ratio tendency due to
                                                  !         gravitational sedimentation

   real    :: q_column_initial                    ! initial value of column integrated mixing ratio
   real    :: q_column_final                      ! final   value of column integrated mixing ratio

   real    :: dt_step, dt_cfl                     ! integration time step
   integer :: n_steps                             ! number of time steps

   integer :: k1, km                              ! indexes
   integer :: n                                   ! loop counter

   integer :: rc                                  ! return code


   k1 = lbound(q, 1)
   km = ubound(q, 1)
   
   allocate(dqdt(k1:km), __STAT__)

   ! If there is no time splitting, the flux is simply the settling flux out from the 
   ! surface layer, because of the mass conservation. Instead of integrating in time,
   ! we calculate the flux as the difference in the column integrated mass before and 
   ! after the settling.

   q_column_initial = 0.0
   q_column_final   = 0.0

   if (present(flux)) then
       q_column_initial = sum(q * delp/g_E)
   end if
  
   ! test if the time step is sufficiently small to maintain numerical stability
   dt_cfl = min(dt, minval(dz / v_t))
   dt_cfl = min(dt_cfl, dz(km) / v_d)

   if (dt_cfl < dt) then
      n_steps = ceiling(dt / dt_cfl)
   else
      n_steps = 1
   end if

   ! time integration
   dt_step = dt / n_steps

   do n = 1, n_steps
       dqdt = 0.0

       ! NOTE: In case of polydisperse aerosols, it would be more correct
       !       to update the settling velocities in the innner time loop.
       call MAML_DryRemovalTendencySolverAerosol(dqdt, q, delp, dz, v_t, v_d)

       q = q + (dqdt * dt_step)
   end do

   ! calculate the flux due to deposition
   if (present(flux)) then
       q_column_final = sum(q * delp/g_E)

       flux = -(q_column_final - q_column_initial) / dt
   end if

   deallocate(dqdt, __STAT__)

   end subroutine MAML_DryRemovalAerosol


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: MAML_DryRemovalGas --- 
!
! !INTERFACE:

   subroutine MAML_DryRemovalGas(q, delp, dz, v_d, dt, flux)

! !USES:

   implicit None


! !INPUT/OUTPUT PARAMETERS:
   real, intent(inout)           :: q         ! mixing ratio,          'kg kg-1' or '# kg-1'
   real, optional, intent(inout) :: flux      ! deposition flux        'kg m-2 s-1' or '# m-2 s-1'

! !INPUT PARAMETERS:
   real, intent(in)              :: delp      ! pressure thickness of levels, 'Pa'
   real, intent(in)              :: dz        ! thickness of levels,   'm'
   real, intent(in)              :: v_d       ! deposotion velocity,   'm s-1'
   real, intent(in)              :: dt        ! model time step,       's'

! !OUTPUT PARAMETERS:


! !DESCRIPTION: Calculates the changes in the mixing ratio due to 
!               gravitational sedimentation/settling.
!
! !REVISION HISTORY:
!
!  6Nov2011  A. Darmenov
!
!EOP
!-------------------------------------------------------------------------
   
                  __Iam__('MAML_DryRemovalGas')

   real    :: dqdt                ! dq/dt - mixing ratio tendency due to
                                  !         gravitational sedimentation

   real    :: q_initial           ! initial value of the mixing ratio


   ! save the initial mixing ratio
   q_initial = q 

   call MAML_DryRemovalTendencySolverGas(dqdt, q, dz, v_d, dt)

   q = q + (dqdt * dt)

   ! calculate the flux due to deposition
   if (present(flux)) then
       flux = -(q - q_initial) * (delp/g_E) / dt
   end if

   end subroutine MAML_DryRemovalGas

   end module MAML_DryRemovalMod

