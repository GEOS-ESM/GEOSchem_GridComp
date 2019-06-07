#include "MAPL_Exceptions.h"
#include "MAPL_Generic.h"


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: MAM_DryRemovalMod - Gravitational sedimentation and deposition 
!                              of aerosol particles
!
! !INTERFACE:
!
   module MAM_DryRemovalMod
!
! !USES:
!
   use ESMF
   use MAPL_Mod


   use MAPL_SimpleBundleMod

   use MAML_SettlingMod
   use MAML_DryDepositionMod
   use MAML_DryRemovalMod

   use MAM_BaseMod

   implicit NONE
   private

!
! !PUBLIC MEMBER FUNCTIONS:

   public MAM_DryRemoval

!
! !PUBLIC PARAMETERS:

!
! !PRIVATE PARAMETERS:
   real, private, parameter :: pi = MAPL_PI
   real, private, parameter :: density_water = MAPL_RHOWTR    ! density of water,  'kg m-3'


!
! !DESCRIPTION: 
!
!  {\tt MAML\_DryRemovalMod} provides a collection of methods for 
!  modeling graviational sedimentation/settling and dry deposition 
!  of aerosol particles
!
!
! !REVISION HISTORY:
!
!  29Nov2011  A. Darmenov  Initial version
!
!EOP
!-------------------------------------------------------------------------


   contains

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: MAM_DryRemoval --- 
!
! !INTERFACE:

   subroutine MAM_DryRemoval(self, import, export, qa, Da, cdt, rc)
! !USES:

   implicit NONE

! !INPUT/OUTPUT PARAMETERS:
   type(MAPL_SimpleBundle), intent(inout) :: qa         ! number/mass mixing ratio
   type(ESMF_State), intent(inout)        :: export     ! export state
   integer, optional, intent(inout)       :: rc         ! return code

! !INPUT PARAMETERS:
   type(MAPL_SimpleBundle), intent(in)    :: Da         ! dry(geometric mean) and wet diameter 
                                                        ! of number size distribution

   type(MAM_Scheme), intent(in)           :: self       ! MAM scheme/configuration 

   type(ESMF_State), intent(inout)        :: import     ! import state

   real,    intent(in)                    :: cdt        ! chemical timestep (secs)

! !OUTPUT PARAMETERS:


! !DESCRIPTION: Gravitational sedimentation and dry deposition of aerosol particles.
!
! !REVISION HISTORY:
!
! 01Dec2011  A. Darmenov  First crack.
! 21Dec2011  A. Darmenov  Dry removal of aerosols.
!
!EOP
!-------------------------------------------------------------------------

              __Iam__('MAM_DryRemoval')

   ! Mode parameters
   ! ---------------
   real                           :: sigma
   character(len=MAM_MAXSTR)      :: mode_name

   integer                        :: n_species
   character(len=MAM_MAXSTR)      :: species_name

   character(len=MAM_MAXSTR)      :: field_name

   integer                        :: m, s
   integer                        :: iq_dgn_dry, iq_dgn_wet
   integer                        :: iq_nmr, iq_mmr
   integer, pointer, dimension(:) :: iq_mmr_species
   real,    pointer, dimension(:) :: density_species, qa_species

   integer                        :: i, im, j, jm, k, km

   real                           :: flux_drydep
   real                           :: density_air, temperature, pressure
   real                           :: diameter_dry, density_dry
   real                           :: diameter_wet, density_wet

   real                           :: f_hg

   real, pointer, dimension(:)    :: vt_nmr, vt_mmr
   real, pointer, dimension(:)    :: dz

   real                           :: vd_nmr, vd_mmr
   real                           :: r_a, r_b
   real                           :: flux_sh
   real                           :: friction_velocity
   real                           :: viscosity_dyn, viscosity_kin
   real                           :: free_mean_path, Kn, Cc, Dp, Sc, St


   !  Input fields from fvGCM
   !  -----------------------
   real, pointer, dimension(:,:,:) :: rhoa
   real, pointer, dimension(:,:,:) :: delp
   real, pointer, dimension(:,:,:) :: ple
   real, pointer, dimension(:,:,:) :: T

   real, pointer, dimension(:,:)   :: SH
   real, pointer, dimension(:,:)   :: z0h
   real, pointer, dimension(:,:)   :: ustar


   !  Exports
   !  -----------------------
   real, pointer, dimension(:,:)   :: flux
   

   integer, parameter :: MOMENT_0 = 0
   integer, parameter :: MOMENT_3 = 3



   !  Get Imports
   !  --------------
   call MAPL_GetPointer(import, rhoa,  'AIRDENS', __RC__)
   call MAPL_GetPointer(import, delp,  'DELP',    __RC__)
   call MAPL_GetPointer(import, ple,   'PLE',     __RC__)
   call MAPL_GetPointer(import, T,     'T',       __RC__)
   call MAPL_GetPointer(import, SH,    'SH',      __RC__)
   call MAPL_GetPointer(import, z0h,   'Z0H',     __RC__)
   call MAPL_GetPointer(import, ustar, 'USTAR',   __RC__)


   !  Local dimensions
   !  ----------------
   im = size(rhoa, 1)
   jm = size(rhoa, 2)
   km = size(rhoa, 3)


   allocate(dz(km),     __STAT__)
   allocate(vt_nmr(km), __STAT__)
   allocate(vt_mmr(km), __STAT__)

   do m = 1, self%n_modes
       call MAM_AerosolModeGet(self%mode(m), name      = mode_name, &
                                             sigma     = sigma,     &
                                             n_species = n_species)

       ! aerosol dry and wet sizes
       field_name = 'DGN_DRY_' // trim(mode_name)
       iq_dgn_dry = MAPL_SimpleBundleGetIndex(Da, field_name, 3, __RC__) 

       field_name = 'DGN_WET_' // trim(mode_name)
       iq_dgn_wet = MAPL_SimpleBundleGetIndex(Da, field_name, 3, __RC__)

       ! number mixing ratio
       field_name = 'NUM_A_' // trim(mode_name)
       iq_nmr = MAPL_SimpleBundleGetIndex(qa, field_name, 3, __RC__)

       ! buffers
       allocate(iq_mmr_species(n_species),  __STAT__)     ! index of the species in the bundle
       allocate(density_species(n_species), __STAT__)     ! density of the species
       allocate(qa_species(n_species),      __STAT__)     ! the mass mixing ratio of the species

       iq_mmr_species  = -1
       qa_species      = 0.0
       density_species = 0.0

       do s = 1, n_species
           species_name = self%mode(m)%species(s)%name
           field_name  = trim(species_name) // '_A_' // trim(mode_name)
           iq_mmr_species(s)  = MAPL_SimpleBundleGetIndex(qa, field_name, 3, __RC__)

           density_species(s) = self%mode(m)%species(s)%component%density
       end do

       do j = 1, jm
           do i = 1, im

               ! calculate the sedimentation velocities in the column
               do k = 1, km
                   ! mid level pressure
                   pressure = (ple(i,j,k-1) + ple(i,j,k)) * 0.5
                   temperature = T(i,j,k)
                   density_air = rhoa(i,j,k)
                
                   ! aerosol size
                   diameter_dry = Da%r3(iq_dgn_dry)%q(i,j,k)
                   diameter_wet = Da%r3(iq_dgn_wet)%q(i,j,k)

                   ! aerosol species mass mixing ratios
                   do s = 1, n_species
                       qa_species(s) = qa%r3(iq_mmr_species(s))%q(i,j,k)
                   end do

                   ! dry density
                   if (any(qa_species > 1.0e-15)) then
                       density_dry = sum(qa_species) / sum(qa_species / density_species)
                   else
                       density_dry = sum(density_species)/n_species
                   end if

                   ! wet density
                   f_hg = max(diameter_wet / diameter_dry, 1.0) ! hygroscopic growth factor
                   density_wet = f_hg**(-3) * density_dry + (1 - f_hg**(-3)) * density_water


                   ! calculate the settling velocity using the wet size and density
                   vt_nmr(k) = MAML_SettlingVelocity(pressure, temperature, &
                                                     diameter_wet, density_wet, sigma, MOMENT_0)
 
                   vt_mmr(k) = MAML_SettlingVelocity(pressure, temperature, &
                                                     diameter_wet, density_wet, sigma, MOMENT_3)
               end do ! k

               dz(:) = delp(i,j,:) / (MAPL_GRAV * rhoa(i,j,:))
               
               ! compute deposition velocity
               temperature       = T(i,j,km)
               density_air       = rhoa(i,j,km)
               flux_sh           = SH(i,j) 
               friction_velocity = ustar(i,j)
               diameter_wet      = Da%r3(iq_dgn_wet)%q(i,j,km) 

               viscosity_dyn = dynamic_viscosity_air(temperature)
               viscosity_kin = kinematic_viscosity_air(temperature, density_air)

               free_mean_path = free_mean_path_air(pressure, temperature)
               Kn = knudsen_number(free_mean_path, diameter_wet)
               Cc = slip_flow_correction(Kn)
            
               Dp = particle_diffusion_coefficient(temperature, viscosity_dyn, Cc, diameter_wet)
               Sc = schmidt_number(viscosity_kin, Dp)
               St = stokes_number(vt_nmr(km), friction_velocity, viscosity_kin)

               r_b = quasi_laminar_resistance(friction_velocity, Sc, St)
               r_a = aerodynamic_resistance(temperature, density_air, flux_sh, friction_velocity, dz(km), z0h(i,j))

               vd_nmr = MAML_DepositionVelocity(vt_nmr(km), r_a, r_b)
               vd_mmr = MAML_DepositionVelocity(vt_mmr(km), r_a, r_b)

               ! update the number mixing ratio in the column
               call MAML_DryRemoval(qa%r3(iq_nmr)%q(i,j,:), delp(i,j,:), &
                                                            dz(:),       &
                                                            vt_nmr(:),   &
                                                            vd_nmr,      &
                                                            cdt)

               ! update the mass mixing ratios in the column and save the mass flux
               do s = 1, n_species
                   species_name = self%mode(m)%species(s)%name
                   iq_mmr = iq_mmr_species(s)

                   ! pointer to the sedimentation flux
                   field_name = 'DP_' // trim(species_name) // '_' // trim(mode_name)

                   call MAPL_GetPointer(export, flux, field_name, __RC__)

                   call MAML_DryRemoval(qa%r3(iq_mmr)%q(i,j,:), delp(i,j,:), &
                                                                dz(:),       &
                                                                vt_mmr(:),   &
                                                                vd_mmr,      &
                                                                cdt,         &
                                                                flux=flux_drydep)
                                                              
                   if (associated(flux)) then
                       flux(i,j) = flux_drydep
                   end if
               end do

           end do ! i
       end do ! j

       deallocate(iq_mmr_species,  __STAT__)
       deallocate(density_species, __STAT__)
       deallocate(qa_species,      __STAT__)
   end do ! m

       
   deallocate(vt_mmr, __STAT__)
   deallocate(vt_nmr, __STAT__)
   deallocate(dz,     __STAT__)

   end subroutine MAM_DryRemoval


end module MAM_DryRemovalMod
