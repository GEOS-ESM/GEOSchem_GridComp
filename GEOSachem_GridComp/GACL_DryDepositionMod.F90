#include "MAPL_Exceptions.h"
#include "MAPL_Generic.h"


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: GACL_DryDepositionMod - Dry deposition of gases and particles
!
! !INTERFACE:
!
   module GACL_DryDepositionMod
!
! !USES:
!
   use MAPL_Mod

   use GACL_ConstantsMod


   implicit NONE
   private

!
! !PUBLIC MEMBER FUNCTIONS:

   public DepositionVelocity

!
! !PUBLIC PARAMETERS:

!
! !PRIVATE PARAMETERS:


!
! !DESCRIPTION: 
!
!  {\tt GACL\_DryDepositionMod} provides a collection of methods for 
!  modeling dry deposition of gases and aerosol particles.
!
!
! !REVISION HISTORY:
!
!  15Dec2011  A. Darmenov  Initial version
!
!EOP
!-------------------------------------------------------------------------


   interface DepositionVelocity
       module procedure DepositionVelocityAerosol
       module procedure DepositionVelocityGas
   end interface DepositionVelocity


   contains

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: schmidt_number --- calculates the Schmidt's number
!
! !INTERFACE:

   pure function schmidt_number(viscosity, D) result (Sc)
! !USES:

   implicit None

   real :: Sc                            ! Schmidt number, ''

! !INPUT/OUTPUT PARAMETERS:

! !INPUT PARAMETERS:

   real, intent(in) :: viscosity         ! kinematic viscosity, 'm2 s-1'
   real, intent(in) :: D                 ! Brownian diffusivity coefficient, ''

! !OUTPUT PARAMETERS:


! !DESCRIPTION: calculates the Schmidt's number (see Eq. 
!
!
! !REVISION HISTORY:
!
!  29Nov2011  A. Darmenov
!
!EOP
!-------------------------------------------------------------------------

   !                __Iam__('schmidt_number')

   Sc = viscosity / D

   end function schmidt_number


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: stokes_number --- calculates the Sokes number
!
! !INTERFACE:

   pure function stokes_number(settling_velocity, friction_velocity, viscosity) result (St)
! !USES:

   use GACL_ConstantsMod, only : g => g_earth

   implicit None

   real :: St                            ! Stokes number, ''

! !INPUT/OUTPUT PARAMETERS:

! !INPUT PARAMETERS:

   real, intent(in) :: settling_velocity ! settling/sedimentation velocity, 'm s-1'
   real, intent(in) :: friction_velocity ! friction velocity,               'm s-1'
   real, intent(in) :: viscosity         ! kinematic viscosity,             'm2 s-1'

! !OUTPUT PARAMETERS:


! !DESCRIPTION: calculates the Stokes number (see Eq. 8.105, Seinfeld and Pandis)
!
! !REVISION HISTORY:
!
!  15Dec2011  A. Darmenov
!
!EOP
!-------------------------------------------------------------------------

   !                __Iam__('stokes_number')

   St = settling_velocity * friction_velocity**2  / (g * viscosity)

   end function stokes_number
   

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: stokes_number --- calculates the quasi-laminar resistance for
!                              particles
!
! !INTERFACE:

   pure function quasi_laminar_resistance(friction_velocity, Sc, St) result (r_b)
! !USES:

   implicit None

   real :: r_b                           ! Quasi-laminar resistance, 'm-1 s'

! !INPUT/OUTPUT PARAMETERS:

! !INPUT PARAMETERS:

   real, intent(in) :: friction_velocity ! friction velocity, 'm s-1'
   real, intent(in) :: Sc                ! Schmidt number,    ''
   real, intent(in) :: St                ! Stokes number,     ''

! !OUTPUT PARAMETERS:


! !DESCRIPTION: calculates the quasi-laminar resistance (see Eq. XX, Seinfeld and Pandis)
!
! !REVISION HISTORY:
!
!  15Dec2011  A. Darmenov
!
!EOP
!-------------------------------------------------------------------------

   !                __Iam__('quasi_laminar_resistance')

   r_b = 1 / (friction_velocity * (Sc**(-0.5) + 10.0**(-3/St)))

   end function quasi_laminar_resistance


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: DepositionVelosityAerosol --- 
!
! !INTERFACE:

   pure function DepositionVelocityAerosol(v_t, r_a, r_b) result (v_d)
! !USES:

   implicit None

   real :: v_d

! !INPUT/OUTPUT PARAMETERS:

! !INPUT PARAMETERS:

   real, intent(in) :: v_t                    ! settling velocity,        'm s-1'

   real, intent(in) :: r_a                    ! aerodynamic resistance,   'm-1 s'
   real, intent(in) :: r_b                    ! quasi-laminar resistance, 'm-1 s'


! !OUTPUT PARAMETERS:


! !DESCRIPTION: Calculates the deposition velocity of particles following  
!               Venkatram and Pleim, 1999.
!
! !REVISION HISTORY:
!
!  15Dec2011  A. Darmenov
!
!EOP
!-------------------------------------------------------------------------
   
   !          __Iam__('DepositionVelocityAerosol')

   v_d = v_t / (1 - exp(-(r_a + r_b) * v_t))

   end function DepositionVelocityAerosol


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GACL_DepositionVelosityAerosol --- 
!
! !INTERFACE:

   pure function DepositionVelocityGas(r_a, r_b) result (v_d)
! !USES:

   implicit None

   real :: v_d

! !INPUT/OUTPUT PARAMETERS:

! !INPUT PARAMETERS:

   real, intent(in) :: r_a                    ! aerodynamic resistance,   'm-1 s'
   real, intent(in) :: r_b                    ! quasi-laminar resistance, 'm-1 s'


! !OUTPUT PARAMETERS:


! !DESCRIPTION: Calculates the deposition velocity of particles following  
!               Venkatram and Pleim, 1999.
!
! !REVISION HISTORY:
!
!  15Dec2011  A. Darmenov
!
!EOP
!-------------------------------------------------------------------------
   
   !          __Iam__('GACL_DepositionVelocityGas')

   v_d = 1 / (r_a + r_b)

   end function DepositionVelocityGas



   end module GACL_DryDepositionMod

