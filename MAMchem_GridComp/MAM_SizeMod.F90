#include "MAPL_Exceptions.h"
#include "MAPL_Generic.h"


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: MAM_SizeMod - Dry and wet size of interstitial aerosols
!
! !INTERFACE:
!
   module MAM_SizeMod
!
! !USES:
!

   use ESMF

   use MAPL_Mod
   use MAPL_SimpleBundleMod

   use MAML_SizeMod
   use MAM_BaseMod


   implicit NONE
   private

!
! !PUBLIC MEMBER FUNCTIONS:

   public MAM_DrySize
   public MAM_WetSize

!
! !PUBLIC PARAMETERS:

!
! !PRIVATE PARAMETERS:

   real, private, parameter :: pi = MAPL_PI


!
! !DESCRIPTION: 
!
!  {\tt MAML\_SizeMod} provides a collection of methods for calculating
!  dry and wet size of aerosol particles.
!
!
! !REVISION HISTORY:
!
!  8Dec2011  A. Darmenov  Initial version
!
!EOP
!-------------------------------------------------------------------------


   contains

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: MAM_DrySize --- calculates the dry size of interstitial aerosols
!
! !INTERFACE:

   subroutine MAM_DrySize(self, import, export, qa, Da, rc)
! !USES:

   implicit NONE

! !INPUT/OUTPUT PARAMETERS:
   type(ESMF_State), intent(inout)     :: import  ! import state
   type(ESMF_State), intent(inout)     :: export  ! export state

   integer, optional, intent(inout)    :: rc      ! return code

! !INPUT PARAMETERS:

   type(MAPL_SimpleBundle), intent(in) :: qa      ! number/mass mixing ratios
   type(MAPL_SimpleBundle), intent(in) :: Da      ! dry sizes

   type(MAM_Scheme), intent(in)        :: self    ! MAM scheme/configuration 

! !OUTPUT PARAMETERS:


! !DESCRIPTION: Calculates the number/volume mean geometric diameter.
!
! !REVISION HISTORY:
!
!  18Nov2011  A. Darmenov  First crack.
!
!EOP
!-------------------------------------------------------------------------

                    __Iam__('MAM_DrySize') 


   real                      :: D           ! number mean geometric diameter
  
   real                      :: sigma
   real                      :: Dg_default, Dg_min, Dg_max
   real                      :: vol2num_default, vol2num_min, vol2num_max
   character(len=MAM_MAXSTR) :: mode_name

   integer                   :: n_species
   character(len=MAM_MAXSTR) :: species_name

   character(len=MAM_MAXSTR) :: field_name

   integer                   :: iq_dgn_dry 
   integer                   :: iq_nmr
   integer, pointer          :: iq_mmr(:)

   integer                   :: m, s
   real, pointer             :: q_mass(:)
   real, pointer             :: density(:)
   real                      :: q_number

   real                      :: f
   integer                   :: i, j, k, im, jm, km


   do m = 1, self%n_modes

       call MAM_AerosolModeGet(self%mode(m), name         = mode_name,  &
                                             sigma        = sigma,      &
                                             size_default = Dg_default, &
                                             size_min     = Dg_min,     &
                                             size_max     = Dg_max,     &
                                             n_species    = n_species)

       ! will need these for the calculations
       f = (pi / 6) * exp(4.5 * log(sigma)**2)      ! volume to number factor
       vol2num_default = 1 / (f * Dg_default**3)
       vol2num_min     = 1 / (f * Dg_min**3)
       vol2num_max     = 1 / (f * Dg_max**3)


       ! size of dry interstitial aerosols
       field_name = 'DGN_DRY_' // trim(mode_name)
       iq_dgn_dry = MAPL_SimpleBundleGetIndex(Da, field_name, 3, __RC__)

       ! number mixing ratio
       field_name = 'NUM_A_' // trim(mode_name)
       iq_nmr = MAPL_SimpleBundleGetIndex(qa, field_name, 3, __RC__)
           
       allocate(q_mass(n_species),  __STAT__)   ! buffer for the species mass mixing ratios
       allocate(density(n_species), __STAT__)   ! buffer for the species densities
       allocate(iq_mmr(n_species),  __STAT__)   ! buffer for the species mass mixing ratios bundle indexes

       density = 0.0
       do s = 1, n_species
           call MAM_AerosolSpeciesGet(self%mode(m)%species(s), name    = species_name, &
                                                               density = density(s))

           field_name = trim(species_name) //'_A_' // trim(mode_name)
           iq_mmr(s) = MAPL_SimpleBundleGetIndex(qa, field_name, 3, __RC__)
       end do

       im = ubound(qa%r3(iq_nmr)%q, 1)
       jm = ubound(qa%r3(iq_nmr)%q, 2)
       km = ubound(qa%r3(iq_nmr)%q, 3)

       do k = 1, km
           do j = 1, jm
               do i = 1, im

                   q_number = qa%r3(iq_nmr)%q(i,j,k)

                   q_mass = 0.0
                   do s = 1, n_species
                       q_mass(s) = qa%r3(iq_mmr(s))%q(i,j,k)     ! OPTIMIZATION_NEEDED: this will probably trash the cache
                   end do

                   D = MAML_DrySize(q_number,        &
                                    q_mass,          &
                                    density,         &
                                    sigma,           &
                                    Dg_default,      &
                                    Dg_min,          &
                                    Dg_max,          &
                                    vol2num_default, &
                                    vol2num_min,     &
                                    vol2num_max)

                   Da%r3(iq_dgn_dry)%q(i,j,k) = D
               end do ! i
           end do ! j
       end do ! k

       deallocate(q_mass, density, iq_mmr, __STAT__)
   end do

   end subroutine MAM_DrySize


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: MAM_WetSize --- calculates the wet size and density of 
!                            interstitial aerosols
!
! !INTERFACE:

   subroutine MAM_WetSize(self, import, export, qa, Da, rc)
! !USES:

       implicit NONE

! !INPUT/OUTPUT PARAMETERS:
       type(ESMF_State), intent(inout)     :: import  ! import state
       type(ESMF_State), intent(inout)     :: export  ! export state

       integer, optional, intent(inout)    :: rc      ! return code

! !INPUT PARAMETERS:

       type(MAPL_SimpleBundle), intent(in) :: qa      ! number/mass mixing ratios
       type(MAPL_SimpleBundle), intent(in) :: Da      ! dry and wet sizes of interstital aerosols

       type(MAM_Scheme), intent(in)        :: self    ! MAM scheme/configuration


! !OUTPUT PARAMETERS:


! !DESCRIPTION: Calculates the number/volume mean geometric diameter.
!
! !REVISION HISTORY:
!
!  18Nov2011  A. Darmenov  First crack.
!
!EOP
!-------------------------------------------------------------------------

                    __Iam__('MAM_WetSize') 


       real :: D           ! number mean geometric diameter
  
       real :: sigma
       real :: rh_deliquescence
       real :: rh_crystallization
       character(len=MAM_MAXSTR) :: mode_name

       integer :: n_species
       character(len=MAM_MAXSTR) :: species_name

       character(len=MAM_MAXSTR) :: field_name

       integer :: iq_mmr_wtr
       integer :: iq_dgn_dry
       integer :: iq_dgn_wet
       integer, allocatable, dimension(:) :: iq_mmr

       integer :: m, s
       real, allocatable, dimension(:) :: q_mass
       real, allocatable, dimension(:) :: density
       real, allocatable, dimension(:) :: hygroscopicity
     
       real :: diameter_dry
       real :: diameter_wet
       real :: density_wet
       real :: q_aerosol_water

       real, pointer, dimension(:,:,:) :: rh
       real, pointer, dimension(:,:,:) :: fcld

       integer :: i, j, k, im, jm, km


       call MAPL_GetPointer(import, rh,   'RH2',  __RC__)
       call MAPL_GetPointer(import, fcld, 'FCLD', __RC__)

       do m = 1, self%n_modes

           call MAM_AerosolModeGet(self%mode(m), name               = mode_name,          &
                                                 sigma              = sigma,              &
                                                 rh_deliquescence   = rh_deliquescence,   &
                                                 rh_crystallization = rh_crystallization, &
                                                 n_species          = n_species)

           ! dry size
           field_name = 'DGN_DRY_' // trim(mode_name)
           iq_dgn_dry = MAPL_SimpleBundleGetIndex(Da, field_name, 3, __RC__)

           ! wet size
           field_name = 'DGN_WET_' // trim(mode_name)
           iq_dgn_wet = MAPL_SimpleBundleGetIndex(Da, field_name, 3, __RC__)

           ! absorbed water  
           field_name = 'WTR_A_' // trim(mode_name)
           iq_mmr_wtr = MAPL_SimpleBundleGetIndex(qa, field_name, 3, __RC__)


           allocate(density(n_species),        __STAT__)   ! buffer for the species densities
           allocate(hygroscopicity(n_species), __STAT__)   ! buffer for the species hygroscopicity

           allocate(q_mass(n_species),         __STAT__)   ! buffer for the species mass mixing ratios
           allocate(iq_mmr(n_species),         __STAT__)   ! buffer for the species mass mixing ratios bundle indexes
                      
           do s = 1, n_species
               call MAM_AerosolSpeciesGet(self%mode(m)%species(s), name           = species_name,    &
                                                                   density        = density(s),      & 
                                                                   hygroscopicity = hygroscopicity(s))

               field_name = trim(species_name) //'_A_' // trim(mode_name)
               iq_mmr(s) = MAPL_SimpleBundleGetIndex(qa, field_name, 3, __RC__)
           end do

           im = size(qa%r3(iq_mmr(1))%q, 1)
           jm = size(qa%r3(iq_mmr(1))%q, 2)
           km = size(qa%r3(iq_mmr(1))%q, 3)

           do k = 1, km
               do j = 1, jm
                   do i = 1, im

                       q_mass = 0.0
                       do s = 1, n_species
                           q_mass(s) = qa%r3(iq_mmr(s))%q(i,j,k)
                       end do

                       diameter_dry = Da%r3(iq_dgn_dry)%q(i,j,k)

                       call MAML_WetSize(q_mass,             &
                                         diameter_dry,       &
                                         density,            &
                                         hygroscopicity,     &
                                         sigma,              &
                                         rh_deliquescence,   &
                                         rh_crystallization, &
                                         rh(i,j,k),          &
                                         fcld(i,j,k),        &
                                         diameter_wet,       &
                                         density_wet,        &
                                         q_aerosol_water)
                      
                       Da%r3(iq_dgn_wet)%q(i,j,k) = diameter_wet

                       qa%r3(iq_mmr_wtr)%q(i,j,k) = q_aerosol_water
                       
                   end do ! i
               end do ! j
           end do ! k

           deallocate(q_mass, density, hygroscopicity, iq_mmr, __STAT__)
       end do

   end subroutine MAM_WetSize


end module MAM_SizeMod
