!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: GACL_ReactionRatesMod - Reaction rate coefficients.
!
! !INTERFACE:
!
   module GACL_ReactionRatesMod
!
! !USES:
!
   use GACL_ConstantsMod, only : pi


   implicit NONE
   private

!
! !PUBLIC MEMBER FUNCTIONS:

   public henry 

!
! !PUBLIC PARAMETERS:

  ! molar enthalpy of dissolution over R_univ  
  real, public, parameter :: E_R_DMS   = -3500.0      ! K
  real, public, parameter :: E_R_MSA   = -0.0         ! K 
  real, public, parameter :: E_R_SO2   = -3120.0      ! K 
  real, public, parameter :: E_R_H2SO4 = -0.0         ! K
  real, public, parameter :: E_R_NH3   = -4085.0      ! K
  real, public, parameter :: E_R_OH    = -4300.0      ! K
  real, public, parameter :: E_R_H2O2  = -6338.0      ! K
  real, public, parameter :: E_R_O3    = -2560.0      ! K

  ! Henry's Low coefficients at T0 = 298.15K
  real, public, parameter :: H_DMS_298   = 5.4e-1     ! M atm-1
  real, public, parameter :: H_MSA_298   = 1.0e20     ! M atm-1   <-- big as in infinity
  real, public, parameter :: H_SO2_298   = 1.2        ! M atm-1
  real, public, parameter :: H_H2SO4_298 = 1.0e11     ! M atm-1
  real, public, parameter :: H_NH3_298   = 58.0       ! M atm-1
  real, public, parameter :: H_OH_298    = 30.0       ! M atm-1
  real, public, parameter :: H_H2O2_298  = 1.0e5      ! M atm-1
  real, public, parameter :: H_O3_298    = 1.2e-2     ! M atm-1


!
! !PRIVATE PARAMETERS:

!
! !DESCRIPTION: 
!
!  {\tt GACL\_ReactionRatesMod} - Reaction rate coefficients.
!
! !REVISION HISTORY:
!
!  22Oct2012  A. Darmenov  Initial version
!
!EOP
!-------------------------------------------------------------------------


   contains

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Henry - computes Henry;s low coefficents 
!
! !INTERFACE:

   pure real function henry(H0, T0, E_R, T)

! !USES:

   implicit None

! !INPUT/OUTPUT PARAMETERS:

! !INPUT PARAMETERS:
   real, intent(in) :: H0            ! Henry's constant at T=T0, M atm-1
   real, intent(in) :: T0            ! temperature, K

   real, intent(in) :: E_R           ! molar enthalpy of dissolution over R_univ, K
    
   real, intent(in) :: T             ! temperature, K



! !OUTPUT PARAMETERS:

! !DESCRIPTION: 
!
!
! !REVISION HISTORY:
!
!  29Sep2012  A. Darmenov
!
!EOP
!-------------------------------------------------------------------------
   
   henry = H0 * exp(-E_R*(1.0/T - 1.0/T0))
   
   end function henry


#if(0)
!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: dynamic_viscosity_air --- 
!
! !INTERFACE:

   pure function dynamic_viscosity_air(T) result (viscosity)
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

   !                __Iam__('dynamic_viscosity_air')

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

   pure function kinematic_viscosity_air(T, density_air) result (viscosity)
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

   !                __Iam__('kinematic_viscosity_air')

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

   pure function free_mean_path_gas(p, T, v, mw) result (path)
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
      
   !               __Iam__('free_mean_path_gas')

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

   pure function free_mean_path_air(p, T) result (path)
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
      
   !               __Iam__('free_mean_path_air')

   real :: v_air  ! viscosity,        'kg m-1 s-1'
   
   v_air = dynamic_viscosity_air(T)
   path  = free_mean_path_gas(p, T, v_air, MW_air)

   end function free_mean_path_air
#endif

   end module GACL_ReactionRatesMod
