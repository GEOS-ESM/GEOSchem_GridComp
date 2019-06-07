#include "MAPL_Generic.h"


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: MAML_SettlingMod - Gravitational sedimentation/settling of aerosol 
!                             particles
!
! !INTERFACE:
!
   module MAML_SettlingMod
!
! !USES:
!
   use MAPL_Mod

   use MAPL_ConstantsMod


   implicit NONE
   private

!
! !PUBLIC MEMBER FUNCTIONS:

   public MAML_SettlingVelocity
   public MAML_Settling

   public dynamic_viscosity_air
   public kinematic_viscosity_air
   public free_mean_path_air
   public knudsen_number
   public slip_flow_correction
   public particle_diffusion_coefficient
   public stokes_settling_velocity


!
! !PUBLIC PARAMETERS:

!
! !PRIVATE PARAMETERS:

   real, private, parameter :: pi     = MAPL_PI

   real, private, parameter :: R_univ = MAPL_RUNIV    ! Universal gas constant,  'J K-1 Kmole-1'
   real, private, parameter :: N_avog = MAPL_AVOGAD   ! Avogadro constant,       'Kmole-1'
   real, private, parameter :: k_B    = R_univ/N_avog ! Boltzmann's constant,    'J K-1'
   real, private, parameter :: MW_air = MAPL_AIRMW    ! molecular weight of air, 'kg Kmole-1'
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


   interface MAML_SettlingVelocity
       module procedure MAML_SettlingVelocityMonodisperseAerosol
       module procedure MAML_SettlingVelocityPolydisperseAerosol
   end interface MAML_SettlingVelocity


   contains

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: dynamic_viscosity_air --- 
!
! !INTERFACE:

   function dynamic_viscosity_air(T) result (viscosity)
! !USES:

   implicit None

   real :: viscosity                     ! dynamic viscosity air, 'kg m-1 s-1' 

! !INPUT/OUTPUT PARAMETERS:

! !INPUT PARAMETERS:

   real, intent(in) :: T                 ! temperature, 'K'

! !OUTPUT PARAMETERS:


! !DESCRIPTION: Calculates the dynamic viscosity of air, following 
!               Sutherland's equation (List, 1984)
!
! !REVISION HISTORY:
!
!  29Nov2011  A. Darmenov
!
!EOP
!-------------------------------------------------------------------------

   viscosity = 1.8325e-5 * (416.16 / (T + 120)) * (T/296.16)**1.5

   end function dynamic_viscosity_air


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: kinematic_viscosity_air --- 
!
! !INTERFACE:

   function kinematic_viscosity_air(T, density_air) result (viscosity)
! !USES:

   implicit None

   real :: viscosity                     ! kinematic viscosity air, 'm2 s-1'

! !INPUT/OUTPUT PARAMETERS:

! !INPUT PARAMETERS:

   real, intent(in) :: T                 ! temperature,    'K'
   real, intent(in) :: density_air       ! density of air, 'kg m-3'

! !OUTPUT PARAMETERS:


! !DESCRIPTION: Calculates the kinematic viscosity of air.
!
! !REVISION HISTORY:
!
!  12Dec2011  A. Darmenov
!
!EOP
!-------------------------------------------------------------------------

   viscosity = dynamic_viscosity_air(T) / density_air

   end function kinematic_viscosity_air



!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: free_mean_path_gas --- 
!
! !INTERFACE:

   function free_mean_path_gas(p, T, v, mw) result (path)
! !USES:

   implicit None

   real :: path                          ! free mean path of gas, 'm' 

! !INPUT/OUTPUT PARAMETERS:

! !INPUT PARAMETERS:

   real, intent(in) :: p                 ! pressure,         'Pa'
   real, intent(in) :: T                 ! temperature,      'K'
   real, intent(in) :: v                 ! viscosity,        'kg m-1 s-1'
   real, intent(in) :: mw                ! molecular weight, 'kg Kmole-1'

! !OUTPUT PARAMETERS:


! !DESCRIPTION: Calculates the free mean path of a pure gas, following 
!               Seinfeld and Pandis, equation 8.6
!
! !REVISION HISTORY:
!
!  29Nov2011  A. Darmenov
!
!EOP
!-------------------------------------------------------------------------
      
   path = 2 * v / (p * (8/pi * mw/(R_univ*T))**0.5)

   end function free_mean_path_gas


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: free_mean_path_gas --- 
!
! !INTERFACE:

   function free_mean_path_air(p, T) result (path)
! !USES:

   implicit None

   real :: path                          ! free mean path of air, 'm' 

! !INPUT/OUTPUT PARAMETERS:

! !INPUT PARAMETERS:

   real, intent(in) :: p                 ! pressure,         'Pa'
   real, intent(in) :: T                 ! temperature,      'K'

! !OUTPUT PARAMETERS:


! !DESCRIPTION: Calculates the free mean path of air molecules
!
! !REVISION HISTORY:
!
!  29Nov2011  A. Darmenov
!
!EOP
!-------------------------------------------------------------------------
      
   real :: v_air  ! viscosity,        'kg m-1 s-1'
   
   v_air = dynamic_viscosity_air(T)
   path  = free_mean_path_gas(p, T, v_air, MW_air)

   end function free_mean_path_air


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: knudsen_number --- 
!
! !INTERFACE:

   function knudsen_number(free_mean_path, diameter) result (Kn)
! !USES:

   implicit None

   real :: Kn                            ! free mean path of air, 'm' 

! !INPUT/OUTPUT PARAMETERS:

! !INPUT PARAMETERS:

   real, intent(in) :: free_mean_path    ! free mean path,        'm'
   real, intent(in) :: diameter          ! diameter,              'm'

! !OUTPUT PARAMETERS:


! !DESCRIPTION: Calculates the Knudsen number
!
! !REVISION HISTORY:
!
!  29Nov2011  A. Darmenov
!
!EOP
!-------------------------------------------------------------------------

   Kn = 2 * free_mean_path / diameter

   end function knudsen_number


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: slip_flow_correction --- 
!
! !INTERFACE:

   function slip_flow_correction(Kn) result (Cc)
! !USES:

   implicit None

   real :: Cc                                     ! Cunningham slip-flow correction

! !INPUT/OUTPUT PARAMETERS:

! !INPUT PARAMETERS:

   real, intent(in)              :: Kn            ! Knudsen number

! !OUTPUT PARAMETERS:


! !DESCRIPTION: Calculates the Cunningham slip-flow correction
!
! !REVISION HISTORY:
!
!  29Nov2011  A. Darmenov
!
!EOP
!-------------------------------------------------------------------------
   
   ! following Seinfeld and Pandis, equation 8.34 / Allen and Raabe (1982)
   real, parameter :: A = 1.257
   real, parameter :: B = 0.4
   real, parameter :: C = 1.1

   Cc = 1 + Kn*(A + B*exp(-C/Kn))

   end function slip_flow_correction


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: slip_flow_correction_ --- 
!
! !INTERFACE:

   function slip_flow_correction_(Kn, linearized) result (Cc)
! !USES:

   implicit None

   real :: Cc                                     ! Cunningham slip-flow correction

! !INPUT/OUTPUT PARAMETERS:

! !INPUT PARAMETERS:

   real, intent(in)              :: Kn            ! Knudsen number
   logical, optional, intent(in) :: linearized    ! linearized formulation (default = false) 

! !OUTPUT PARAMETERS:


! !DESCRIPTION: Calculates the Cunningham slip-flow correction
!
! !REVISION HISTORY:
!
!  29Nov2011  A. Darmenov
!
!EOP
!-------------------------------------------------------------------------

   logical         :: linearized_form
   
   ! following Seinfeld and Pandis, equation 8.34 / Allen and Raabe (1982)
   real, parameter :: A = 1.257
   real, parameter :: B = 0.4
   real, parameter :: C = 1.1

   ! linearized form following Binkowski and Shankar (equation A27, 1995)
   real, parameter :: A_lf = 1.246


   if (present(linearized)) then
       linearized_form = linearized
   else
       linearized_form = .false.
   end if

   if (linearized_form) then
       Cc = 1 + Kn*A_lf
   else
       Cc = 1 + Kn*(A + B*exp(-C/Kn))
   end if

   end function slip_flow_correction_


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: particle_diffusion_coefficient --- 
!
! !INTERFACE:

   function particle_diffusion_coefficient(temperature,   &
                                           viscosity_air, &
                                           Cc,            &
                                           diameter) result (Dp)
! !USES:

   implicit None

   real :: Dp

! !INPUT/OUTPUT PARAMETERS:

! !INPUT PARAMETERS:

   real, intent(in) :: temperature       ! temperature,               'K'
   real, intent(in) :: viscosity_air     ! dynamic viscosity of air,  'kg m-1 s-1'

   real, intent(in) :: diameter          ! particle diameter,         'm'

   real, intent(in) :: Cc                ! slip-flow correction factor
   

! !OUTPUT PARAMETERS:


! !DESCRIPTION: Calculates the Brownian particle diffusivity.
!
! !REVISION HISTORY:
!
!  29Nov2011  A. Darmenov
!
!EOP
!-------------------------------------------------------------------------
     
   Dp = (k_B * temperature / (3 * pi * viscosity_air * diameter)) * Cc

   end function particle_diffusion_coefficient


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: particle_diffusion_coefficient --- 
!
! !INTERFACE:

   function stokes_settling_velocity(viscosity_air, &
                                     diameter,      &
                                     density) result (v_s)
! !USES:

   implicit None

   real :: v_s

! !INPUT/OUTPUT PARAMETERS:

! !INPUT PARAMETERS:

   real, intent(in) :: viscosity_air     ! dynamic viscosity of air,  'kg m-1 s-1'

   real, intent(in) :: diameter          ! particle diameter,         'm'
   real, intent(in) :: density           ! particle density,          'm'


! !OUTPUT PARAMETERS:


! !DESCRIPTION: Calculates the Stokes velocity.
!
! !REVISION HISTORY:
!
!  29Nov2011  A. Darmenov
!
!EOP
!-------------------------------------------------------------------------
     
   v_s = (g_E / 18) * (density * diameter**2 / viscosity_air)

   end function stokes_settling_velocity

   

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: MAML_SettlingVelosityMonodisperseAerosol --- 
!
! !INTERFACE:

   function MAML_SettlingVelocityMonodisperseAerosol(pressure,    &
                                                     temperature, &
                                                     diameter,    &
                                                     density) result (v_t)
! !USES:

   implicit None

   real :: v_t

! !INPUT/OUTPUT PARAMETERS:

! !INPUT PARAMETERS:

   real, intent(in) :: pressure               ! pressure,          'Pa'
   real, intent(in) :: temperature            ! temperature,       'K'

   real, intent(in) :: diameter               ! particle diameter, 'm'
   real, intent(in) :: density                ! particle density,  'kg m-3'
   

! !OUTPUT PARAMETERS:


! !DESCRIPTION: Calculates the terminal settling velocity of 
!               monodisperse spherical particles. 
!               (Reference: Seinfeld and Pandis, equation 8.42)
!
! !REVISION HISTORY:
!
!  29Nov2011  A. Darmenov
!
!EOP
!-------------------------------------------------------------------------
   
   real :: viscosity      ! viscosity of air,  'kg m-1 s-1'
   real :: free_mean_path ! free mean path
   real :: Kn             ! Knudsen number
   real :: Cc             ! slip-flow correction
   real :: v_s            ! Stokes velocity


   viscosity = dynamic_viscosity_air(temperature)
   free_mean_path = free_mean_path_air(pressure, temperature)

   Kn = knudsen_number(free_mean_path, diameter)

   Cc = slip_flow_correction(Kn)
   v_s = stokes_settling_velocity(viscosity, diameter, density)

   v_t = v_s * Cc

   end function MAML_SettlingVelocityMonodisperseAerosol


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: MAML_SettlingVelocityPolydisperseAerosol --- 
!
! !INTERFACE:

   function MAML_SettlingVelocityPolydisperseAerosol(pressure,    &
                                                     temperature, &
                                                     diameter,    &
                                                     density,     &
                                                     sigma,       &
                                                     k) result (v_t)
! !USES:

   implicit None

   real :: v_t

! !INPUT/OUTPUT PARAMETERS:

! !INPUT PARAMETERS:

   real, intent(in)    :: pressure            ! pressure,          'Pa'
   real, intent(in)    :: temperature         ! temperature,       'K'

   real, intent(in)    :: diameter            ! geometric mean diameter, 'm'
   real, intent(in)    :: density             ! particle density,  'kg m-3'

   real, intent(in)    :: sigma               ! geometric standard deviation
   integer, intent(in) :: k                   ! k-th moment of the number size distribution
   

! !OUTPUT PARAMETERS:


! !DESCRIPTION: Calculates the terminal settling velocity of 
!               polydisperse spherical particles. Note that the expression 
!               for monodisperse distribution settling velocity is recovered
!               if the geometric standard deviatin is set to 1.
!               References: Binkowski and Shankar, 1995; Tulet et al., 2005 
!
! !REVISION HISTORY:
!
!  29Nov2011  A. Darmenov
!
!EOP
!-------------------------------------------------------------------------
   
   real :: viscosity      ! dynamic viscosity of air,  'kg m-1 s-1'
   real :: free_mean_path ! free mean path
   real :: Kn             ! Knudsen number
   real :: Cc             ! Cunningham-slip-flow correction term
   real :: v_s            ! Stokes terminal velocity

   real :: ln2_sigma      ! log^2(sigma)

   ! linearized form of slip-flow correction following Binkowski and Shankar 
   ! (equation A27, 1995)
   real, parameter :: A = 1.246

   ln2_sigma = (log(sigma))**2

   viscosity = dynamic_viscosity_air(temperature)
   free_mean_path = free_mean_path_air(pressure, temperature)

   Kn = knudsen_number(free_mean_path, diameter)

   Cc = exp((2*k + 2) * ln2_sigma) + A * Kn * exp((k + 0.5) * ln2_sigma)
   v_s = stokes_settling_velocity(viscosity, diameter, density)

   v_t = v_s * Cc

   end function MAML_SettlingVelocityPolydisperseAerosol


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: MAML_SettlingTendencySolver --- 
!
! !INTERFACE:

   subroutine MAML_SettlingTendencySolver(dqdt, q, delp, dz, v_t)

! !USES:

   implicit None


! !INPUT/OUTPUT PARAMETERS:
   real, dimension(:), optional, intent(out) :: dqdt      ! dq/dt - mixing ratio tendency due to
                                                          !         gravitational sedimentation

! !INPUT PARAMETERS:

   real, dimension(:), intent(in)    :: q         ! mixing ratio,                 'kg kg-1' or '# kg-1'
   real, dimension(:), intent(in)    :: delp      ! pressure thickness of levels, 'Pa'
   real, dimension(:), intent(in)    :: dz        ! thickness of levels,          'm'
   real, dimension(:), intent(in)    :: v_t       ! settling velocity,            'm s-1'

! !OUTPUT PARAMETERS:


! !DESCRIPTION: Calculates mixing ratio tendency (dq/dt) due to 
!               gravitational sedimentation/settling.
!
! !REVISION HISTORY:
!
!  29Nov2011  A. Darmenov
!
!EOP
!-------------------------------------------------------------------------
   
   integer :: k1, k2, km

   k1 = lbound(q, 1)
   km = ubound(q, 1)

   k2   = k1 + 1
   dqdt = 0.0

   ! zero flux in from top of the atmosphere (k=k1)
   dqdt(k1) = 0.0 - q(k1) * v_t(k1) / dz(k1)

   ! levels k2:km --- flux in from the level above and flux out into the level below
   dqdt(k2:km) = (q(k1:km-1) * (v_t(k1:km-1) / dz(k1:km-1))) * (delp(k1:km-1)/delp(k2:km)) - &
                 (q(k2:km  ) * (v_t(k2:km  ) / dz(k2:km  )))

   end subroutine MAML_SettlingTendencySolver


   
!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: MAML_Settling --- 
!
! !INTERFACE:

   subroutine MAML_Settling(q, delp, dz, v_t, dt, flux)

! !USES:

   implicit None


! !INPUT/OUTPUT PARAMETERS:
   real, dimension(:), intent(inout) :: q         ! mixing ratio,                 'kg kg-1' or '# kg-1'
   real, optional,     intent(inout) :: flux      ! surface flux                  'kg m-2 s-1' or '# m-2 s-1'

! !INPUT PARAMETERS:
   real, dimension(:), intent(in)    :: delp      ! pressure thickness of levels, 'Pa'
   real, dimension(:), intent(in)    :: dz        ! thickness of levels,          'm'
   real, dimension(:), intent(in)    :: v_t       ! settling velocity,            'm s-1'
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
   
                  __Iam__('MAML_Settling')

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
       call MAML_SettlingTendencySolver(dqdt, q, delp, dz, v_t)

       q = q + (dqdt * dt_step)
   end do

   ! calculate the flux due to settling
   if (present(flux)) then
       q_column_final = sum(q * delp/g_E)

       flux = -(q_column_final - q_column_initial) / dt
   end if

   deallocate(dqdt, __STAT__)

   end subroutine MAML_Settling


   end module MAML_SettlingMod

