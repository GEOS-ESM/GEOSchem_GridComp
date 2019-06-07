#include "MAPL_Exceptions.h"
#include "MAPL_Generic.h"

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: MAML_SizeMod - methods for calculation of dry and wet sizes.
!
! !INTERFACE:
!
   module MAML_SizeMod
!
! !USES:
!

   use MAPL_Mod
   use MAPL_ConstantsMod, only : MAPL_PI, MAPL_RHOWTR, r8 => MAPL_R8

   use modal_aero_wateruptake, only : modal_aero_kohler


   implicit NONE
   private

!
! !PUBLIC MEMBER FUNCTIONS:

   public MAML_DrySize
   public MAML_WetSize


! !PRIVATE PARAMETERS 
   real, private, parameter :: pi = MAPL_PI
   real, private, parameter :: density_water = MAPL_RHOWTR    ! density of water,  'kg m-3'


!
! !DESCRIPTION: 
!
!  {\tt MAML\_SizeMod} provides a collection of methods to calculate
!  dry and wet size of aerosol particles.
!
!
! !REVISION HISTORY:
!
!  18Nov2011  A. Darmenov  Initial version
!
!EOP
!-------------------------------------------------------------------------



   contains


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: MAML_DrySize --- 
!
! !INTERFACE:

   function MAML_DrySize(q_number,        &
                         q_mass,          &
                         density,         & 
                         sigma,           &
                         Dg_default,      &
                         Dg_min,          &
                         Dg_max,          &
                         vol2num_default, &
                         vol2num_min,     &
                         vol2num_max)     result (Dg_num)
! !USES:

   implicit NONE

   real :: Dg_num

! !INPUT/OUTPUT PARAMETERS:

! !INPUT PARAMETERS:

   real, intent(in)               :: q_number          ! number mixing ratio
   real, dimension(:), intent(in) :: q_mass            ! mass mixing ratio of all components
   real, dimension(:), intent(in) :: density           ! bulk density of all components
   real, intent(in)               :: sigma             ! geometric standard deviation

   real, intent(in)               :: Dg_default        ! default geometric mean number of number size distribution
   real, intent(in)               :: Dg_min            ! minimum geometric mean number of number size distribution
   real, intent(in)               :: Dg_max            ! maximum geometric mean number of number size distribution

   real, intent(in)               :: vol2num_default   ! default volume to number ratio
   real, intent(in)               :: vol2num_min       ! minimum volume to number ratio
   real, intent(in)               :: vol2num_max       ! maximum volume to number ratio

! !OUTPUT PARAMETERS:


! !DESCRIPTION: Calculates the geometric mean diameter of number size distribution.
!
! !REVISION HISTORY:
!
!  17Nov2011  A. Darmenov  First crack.
!
!EOP
!-------------------------------------------------------------------------
                   __Iam__('MAML_DrySize')
       
   ! local variables
   real :: vol         ! volume
   real :: num         ! number mixing ratio
   real :: vol2num     ! mass mixing ratio
   real :: Dg          ! geometric mean diameter

   real :: f           ! factor

   ! set the default values
   Dg = Dg_default
   vol2num = vol2num_default

   ! compute volume mixing ratio
   vol = sum(q_mass/density)

   if (vol > 0) then
       num = max(0.0, q_number)

       if (num <= vol*vol2num_max) then
           Dg = Dg_max
           vol2num = vol2num_max
       else if (num >= vol*vol2num_min) then
           Dg = Dg_min
           vol2num = vol2num_min
       else
           ! lognormal size distribution factor
           f = (pi / 6) * exp(4.5 * log(sigma)**2)

           Dg = (vol / (f * num))**(1/3.)
           vol2num = num / vol
       end if
   end if

   Dg_num = Dg

   end function MAML_DrySize
   

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: MAML_WetSize --- 
!
! !INTERFACE:

   subroutine MAML_WetSize(q_mass,             &
                           Dg_dry,             &
                           density,            &
                           hygroscopicity,     &
                           sigma,              &
                           rh_deliquescence,   &
                           rh_crystallization, &
                           rh,                 &
                           f_cld,              &
                           Dg_wet,             &
                           density_wet,        &
                           q_aerosol_water)
! !USES:

   implicit NONE


! !INPUT/OUTPUT PARAMETERS:
   real, intent(inout)            :: Dg_wet             ! aerosol wet size,     'm'
   real, intent(inout)            :: density_wet        ! aerosol wet density,  'kg m-3'
   real, intent(inout)            :: q_aerosol_water    ! mass mixing ratio of absorbed
                                                        ! by aerosol water,     'kg kg-1'

! !INPUT PARAMETERS:

   real, dimension(:), intent(in) :: q_mass             ! mass mixing ratio of all components
   real, dimension(:), intent(in) :: density            ! bulk density of all components
   real, dimension(:), intent(in) :: hygroscopicity     ! hygroscopicity of all components
   real, intent(in)               :: sigma              ! geometric standard deviation

   real, intent(in)               :: Dg_dry             ! Dry size -- geometric mean diameter of number size distribution

   real, intent(in)               :: rh_deliquescence   ! deliquescence RH point
   real, intent(in)               :: rh_crystallization ! crystallization RH point

   real, intent(in)               :: rh                 ! relative humidity
   real, intent(in)               :: f_cld              ! cloud fraction

! !OUTPUT PARAMETERS:


! !DESCRIPTION: Calculates wet size and density of aerosol particles.
!
! !REVISION HISTORY:
!
!  12Dec2011  A. Darmenov  First crack -- based on modal_aero_wateruptake_sub(),
!                          from CESM-1.0.3
!
!EOP
!-------------------------------------------------------------------------

                   __Iam__('MAML_WetSize')


   real, allocatable, dimension(:) :: q_vol_dry    ! dry volume mixing ratios

   real     :: vol_dry                             ! total dry volume (mixing ratio)
   real     :: mass_dry                            ! total dry mass (mixing ratio)
   real     :: number_dry                          ! total number of particles
   real     :: density_dry                         ! dry density

   real     :: f                                   ! lognormal size distribution factor
   real     :: vol2num                             ! volume to number ratio
   real     :: Dg                                  ! geometric mean diameter

   real     :: f_hysteresis                        ! hysteresis factor

   real     :: particle_number_dry                 ! one particle, that is equal to 1
   real     :: particle_hygroscopicity             !  - volume average hygroscopicity
   real     :: particle_vol_dry                    !  - dry volume
   real     :: particle_mass_dry                   !  - dry mass
   real     :: particle_radius_dry                 !  - dry radius
   real     :: particle_vol_wet                    !  - wet volume 
   real     :: particle_radius_wet                 !  - wet radius
   real     :: particle_vol_water                  !  - volume of aerosol water


   real     :: rh_clr

   real(r8) :: rh_(1)
   real(r8) :: particle_radius_dry_(1)
   real(r8) :: particle_hygroscopicity_(1)
   real(r8) :: particle_radius_wet_(1)

   integer  :: n_species
   integer  :: rc                                 ! return code

   real, parameter :: VMR_DRY_MIN = 1e-30         ! minimum volume mixing ratio
   real, parameter :: MMR_DRY_MIN = 1e-31         ! minimum mass mixing ratio

   real, parameter :: pi_43 = (4/3.0) * pi
   
   real, parameter :: third = (1/3.0)



   ! clear portion RH
   rh_clr = rh
   rh_clr = max(rh_clr, 0.00)
   rh_clr = min(rh_clr, 0.98)

   if (f_cld < 1.0) then
       rh_clr = (rh_clr - 1.0*f_cld) / (1 - f_cld)
   end if

   rh_clr = max(rh_clr, 0.00)


   n_species = size(q_mass)
   allocate(q_vol_dry(n_species), __STAT__)

   ! volume average hygroscopicity
   q_vol_dry = max(0.0, q_mass/density)
   vol_dry   = sum(q_vol_dry)                      ! the total dry volume

   if (vol_dry > VMR_DRY_MIN) then
       particle_hygroscopicity = sum(hygroscopicity * q_vol_dry) / vol_dry
   else
       particle_hygroscopicity = sum(hygroscopicity) / size(hygroscopicity)
   end if

   deallocate(q_vol_dry, __STAT__)

   ! (volume) average density
   mass_dry = sum(max(0.0, q_mass))

   if (mass_dry > 1e-31) then
       density_dry = mass_dry / vol_dry
   else
       density_dry = sum(density) / size(density)
   end if

   
   ! dry volume to number factor
   f = (pi / 6) * exp(4.5 * (log(sigma))**2)
   vol2num = 1 / (f * Dg_dry**3)

   number_dry = vol2num * vol_dry     ! why not pass the number mixing ratio as an input argument?
                                      ! besides, there might be inconsistencies introduced by 
                                      ! using max()/min() range tests


   ! mean (single particle) dry volume and mass
   particle_number_dry = 1.0
   particle_vol_dry    = particle_number_dry / vol2num
   particle_mass_dry   = density_dry * particle_vol_dry
   particle_radius_dry = (particle_vol_dry / pi_43)**third
   
   ! compute the wet radius
   rh_(1) = rh_clr
   particle_radius_dry_(1)     = particle_radius_dry
   particle_hygroscopicity_(1) = particle_hygroscopicity

   call modal_aero_kohler(particle_radius_dry_,     &
                          particle_hygroscopicity_, &
                          rh_,                      &
                          particle_radius_wet_, 1, 1)

   particle_radius_wet = particle_radius_wet_(1)

   particle_radius_wet = max(particle_radius_wet, particle_radius_dry)

   ! swell the particle by applying the ratio of wet to dry size
   Dg_wet = Dg_dry * (particle_radius_wet / particle_radius_dry)

   ! aerosol water volume and mass
   particle_vol_wet = pi_43 * particle_radius_wet**3
   particle_vol_wet = max(particle_vol_wet, particle_vol_dry)

   particle_vol_water = particle_vol_wet - particle_vol_dry
   particle_vol_water = max(0.0, particle_vol_water)


   ! Simple treatment of deliquesence/crystallization hysteresis -- 
   ! for rhcrystal < rh < rhdeliques, aerosol water is a fraction of the 
   ! "upper curve" value, and the fraction is a linear function of RH

   if (rh_clr < rh_crystallization) then
       particle_radius_wet = particle_radius_dry
       particle_vol_wet    = particle_vol_dry
       particle_vol_water  = 0.0
   else if (rh_clr < rh_deliquescence) then
       f_hysteresis = 1.0 / max(1.0e-5, (rh_deliquescence - rh_crystallization))

       particle_vol_water = f_hysteresis * (rh_clr - rh_crystallization) * particle_vol_water
       particle_vol_water = max(0.0, particle_vol_water)

       particle_vol_wet = particle_vol_dry + particle_vol_water
       particle_radius_wet = (particle_vol_wet / pi_43)**third
   end if

 
   ! water absorbed by the aerosols
   q_aerosol_water = density_water * number_dry * particle_vol_water


   ! aerosol wet density (kg/m3)
   if (particle_vol_wet > VMR_DRY_MIN) then
       density_wet = (particle_mass_dry + (density_water * particle_vol_water))/particle_vol_wet
   else
       density_wet = density_dry
   end if
   
   end subroutine MAML_WetSize


   end module MAML_SizeMod
