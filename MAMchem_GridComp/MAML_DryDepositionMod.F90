#include "MAPL_Exceptions.h"
#include "MAPL_Generic.h"


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: MAML_DryDepositionMod - Dry deposition of gases and particles
!
! !INTERFACE:
!
   module MAML_DryDepositionMod
!
! !USES:
!
   use MAPL_Mod

   use MAPL_ConstantsMod


   implicit NONE
   private

!
! !PUBLIC MEMBER FUNCTIONS:

   public MAML_DepositionVelocity

   public schmidt_number
   public stokes_number
   public quasi_laminar_resistance
   public aerodynamic_resistance

!
! !PUBLIC PARAMETERS:

!
! !PRIVATE PARAMETERS:

   real, private, parameter :: g_E        = MAPL_GRAV     ! standard gravity,        'm s-2'
   real, private, parameter :: von_karman = MAPL_KARMAN   ! von Karman's constant,   '1'
   real, private, parameter :: c_p        = MAPL_CP       ! specific heat capacity of dry air, 'J kg-1 K01'

!
! !DESCRIPTION: 
!
!  {\tt MAML\_DryDepositionMod} provides a collection of methods for 
!  modeling dry deposition of gases and aerosol particles.
!
!
! !REVISION HISTORY:
!
!  15Dec2011  A. Darmenov  Initial version
!
!EOP
!-------------------------------------------------------------------------


   interface MAML_DepositionVelocity
       module procedure MAML_DepositionVelocityAerosol
!      module procedure MAML_DepositionVelocityGas
   end interface MAML_DepositionVelocity


   contains

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: schmidt_number --- calculates the Schmidt's number
!
! !INTERFACE:

   function schmidt_number(viscosity, D) result (Sc)
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

                    __Iam__('schmidt_number')

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

   function stokes_number(settling_velocity, friction_velocity, viscosity) result (St)
! !USES:

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

                    __Iam__('stokes_number')

   St = settling_velocity * friction_velocity**2  / (g_E * viscosity)

   end function stokes_number
   

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: quasi_laminar_resistance --- calculates the quasi-laminar 
!            resistance for particles
!
! !INTERFACE:

   function quasi_laminar_resistance(friction_velocity, Sc, St) result (r_b)
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

                    __Iam__('quasi_laminar_resistance')

   r_b = 1 / (friction_velocity * (Sc**(-0.5) + 10.0**(-3/St)))

   end function quasi_laminar_resistance


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: aerodynamic_resistance --- calculates the aerodynamic resistance
!
! !INTERFACE:

   function aerodynamic_resistance(temperature, density_air, flux_sh, &
                                   friction_velocity, dz, z0h) result (r_a)
! !USES:

   implicit None

   real :: r_a                           ! Aerodynamic resistance, 'm-1 s'

! !INPUT/OUTPUT PARAMETERS:

! !INPUT PARAMETERS:

   real, intent(in) :: temperature       ! temperature,       'K'
   real, intent(in) :: density_air       ! density of air,    'kg m-3'
   real, intent(in) :: flux_sh           ! sensible heat flux at the surface, 'W m-2'
   real, intent(in) :: friction_velocity ! friction velocity, 'm s-1'
   real, intent(in) :: dz                ! depth of the surface layer, 'm'
   real, intent(in) :: z0h               ! roughness height for sensible heat, 'm'

! !OUTPUT PARAMETERS:


! !DESCRIPTION: calculates the aerodynamic resistance (see Eq. XX, Seinfeld and Pandis)
!
! !REVISION HISTORY:
!
!  15Dec2011  A. Darmenov
!
!EOP
!-------------------------------------------------------------------------

                    __Iam__('aerodynamic_resistance')

   ! local
   real, parameter :: k = von_karman

   real :: z_ref
   real :: f
   real :: psi_h
   real :: log_f
   real :: z0h_
   real :: eps
   real :: L


   z_ref = 0.5 * dz

   L = monin_obukhov_length(temperature, density_air, flux_sh, friction_velocity)

   f = z_ref / L

   if(f > 1.0) then
       f = 1.0
   end if

   if ( (f > 0.0) .and. (f <= 1.0)) then
       psi_h = -5.0*f
   else if (f < 0.0) then
       eps = min(1.0, -f)
       log_f = log(eps)
       psi_h = exp(0.598 + 0.39*log_f - 0.09*(log_f**2))
   endif

   z0h_ = max(z0h, 1e2 * tiny(1.0))

   r_a = (log(z_ref / z0h_) - psi_h) / (k * friction_velocity)

   end function aerodynamic_resistance


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: monin_obukhov_length --- calculates the Monin-Obukhov length
!
! !INTERFACE:

   function monin_obukhov_length(temperature, density_air, flux_sh, friction_velocity) result (L)
! !USES:

   implicit None

   real :: L                             !  Monin-Obukhov length, 'm'

! !INPUT/OUTPUT PARAMETERS:

! !INPUT PARAMETERS:

   real, intent(in) :: temperature       ! temperature,       'K'
   real, intent(in) :: density_air       ! density of air,    'kg m-3'
   real, intent(in) :: flux_sh           ! sensible heat flux at the surface, 'W m-2'
   real, intent(in) :: friction_velocity ! friction velocity, 'm s-1'

! !OUTPUT PARAMETERS:


! !DESCRIPTION: calculates the calculates the Monin-Obukhov length
!
! !REVISION HISTORY:
!
!  15Dec2011  A. Darmenov
!
!EOP
!-------------------------------------------------------------------------

                    __Iam__('monin_obukhov_length')

   ! local
   real, parameter :: k = von_karman 
   real, parameter :: g = g_E

   if (abs(flux_sh) > 1e3*epsilon(0.0)) then
       L = - density_air * c_p * temperature * friction_velocity**3 / (k * g * flux_sh)
   else
       L = 1/(1e3*epsilon(0.0))
   end if

   end function monin_obukhov_length


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: MAML_DepositionVelosityAerosol --- 
!
! !INTERFACE:

   function MAML_DepositionVelocityAerosol(v_t, r_a, r_b) result (v_d)
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
   
              __Iam__('MAML_DepositionVelocityAerosol')

   v_d = v_t / (1 - exp(-(r_a + r_b) * v_t))

   end function MAML_DepositionVelocityAerosol

   end module MAML_DryDepositionMod

