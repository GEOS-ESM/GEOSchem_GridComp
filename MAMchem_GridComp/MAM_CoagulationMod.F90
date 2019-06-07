#include "MAPL_Exceptions.h"
#include "MAPL_Generic.h"


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: MAM_CoagulationMod - Coagulation of interstitial aerosols
!
! !INTERFACE:
!
   module MAM_CoagulationMod
!
! !USES:
!

   use ESMF

   use MAPL_Mod
   use MAPL_SimpleBundleMod

   use MAML_CoagulationMod

   use MAM3_DataMod 
   use MAM7_DataMod

   use MAM_BaseMod


   implicit NONE
   private

!
! !PUBLIC MEMBER FUNCTIONS:

   public MAM_CoagulationBimodal
#if (1)
   public MAM_Coagulation
#endif   

!
! !PUBLIC PARAMETERS:

!
! !PRIVATE PARAMETERS:


!
! !DESCRIPTION: 
!
!  {\tt MAML\_CoagulationMod} provides methods for computing the changes  
!  of number and mass mixing ratios of MAM interstitial aerosols.
!
!
! !REVISION HISTORY:
!
!  9Jan2012  A. Darmenov  Initial version
!
!EOP
!-------------------------------------------------------------------------


   contains


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: MAM_CoagulationBimodal --- models the coagualtion of interstitial aerosols
!
!
! !INTERFACE:

   subroutine MAM_CoagulationBimodal(self, import, export, qa, Da, cdt, rc)
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


! !DESCRIPTION: Coagulation of interstitial aerosols: 
!               - MAM7 scheme models coagulation by treating intramodal 
!               coagulation of Aitken, accumulation and primary carbon (PCM) 
!               modes as well as intramodal coagulation from Aitken to 
!               accumulation, primary carbon to accumulation and aged Aitken 
!               via primary carbon to accumulation mode.
!               - MAM3 scheme models intra- and intermodal coagulation of 
!               the Aitken and accumulation modes.
!
! !REVISION HISTORY:
!
!  09Jan2012  A. Darmenov  First crack.
!
!EOP
!-------------------------------------------------------------------------

              __Iam__('MAM_CoagulationBimodal')


   ! Mode parameters
   ! ---------------
   real                      :: sigma
   character(len=MAM_MAXSTR) :: mode_name

   integer                   :: n_species
   character(len=MAM_MAXSTR) :: species_name

   character(len=MAM_MAXSTR) :: field_name

   ! indexes 
   integer :: m, s, n 
   integer :: i, im, j, jm, k, km

   integer :: n_coag_modes                                              ! number of coagulation modes
   integer :: coag_species_number_max                                   ! max number of species

   real    :: temperature                                               ! local temperature
   real    :: pressure                                                  ! local mid-level pressure
   real    :: density_air                                               ! local air density

   real    :: f_hg                                                      ! hygroscopic growth factor

   ! buffers
   integer, allocatable, dimension(  :) :: coag_mode_index              ! indexes of the modes included in the coagulation
   integer, allocatable, dimension(  :) :: coag_mode_species_number     ! numbner of species in each mode

   integer, allocatable, dimension(  :) :: coag_intermodal_transfer     ! intermodal coagulation mapping

   integer, allocatable, dimension(  :) :: iq_dgn_dry, iq_dgn_wet       ! dry and wet size bundle indexes

   integer, allocatable, dimension(  :) :: iq_nmr                       ! number mixing ratio bundle indexes
   integer, allocatable, dimension(:,:) :: iq_mmr                       ! mass mixing ratio bundle indexes


   real,    allocatable, dimension(  :) :: qa_number                    ! number mixing ratios
   real,    allocatable, dimension(:,:) :: qa_mass                      ! mass mixing ratios
   real,    allocatable, dimension(:,:) :: density_species              ! species bulk densities

   real,    allocatable, dimension(  :) :: diameter_dry, diameter_wet   ! dry and wet aerosol sizes
   real,    allocatable, dimension(  :) :: density_dry,  density_wet    ! dry and wet aerosol densities
   real,    allocatable, dimension(  :) :: coag_mode_sigma              ! geometric standard deviations

   

   ! other derived variables
   real,    allocatable, dimension(:)   :: dz


   !  Input fields from fvGCM
   !  -----------------------
   real, pointer, dimension(:,:,:) :: rhoa
   real, pointer, dimension(:,:,:) :: delp
   real, pointer, dimension(:,:,:) :: ple
   real, pointer, dimension(:,:,:) :: T

   !  Exports
   !  -----------------------
   real, pointer, dimension(:,:)   :: flux
  

   real, parameter :: density_water = 1000.0 ! kg m-3




   !  Get Imports
   !  --------------
   call MAPL_GetPointer(import, rhoa, 'AIRDENS', __RC__)
   call MAPL_GetPointer(import, delp, 'DELP',    __RC__)
   call MAPL_GetPointer(import, ple,  'PLE',     __RC__)
   call MAPL_GetPointer(import, T,    'T',       __RC__)


   !  Local dimensions
   !  ----------------
   im = size(rhoa, 1)
   jm = size(rhoa, 2)
   km = size(rhoa, 3)


   !  Find the coagulating modes
   !  --------------------------
   n_coag_modes = 2
   allocate(coag_mode_index(n_coag_modes),          __STAT__)
   allocate(coag_mode_species_number(n_coag_modes), __STAT__)

   coag_mode_index          = 0
   coag_mode_species_number = 0

   do m = 1, self%n_modes
       call MAM_AerosolModeGet(self%mode(m), name      = mode_name, &
                                             n_species = n_species  ) 

       if (trim(mode_name) == trim(MAM7_AITKEN_MODE_NAME)) then
           coag_mode_index(1) = m
           coag_mode_species_number(1) = n_species 
       else if (trim(mode_name) == trim(MAM7_ACCUMULATION_MODE_NAME)) then
           coag_mode_index(2) = m
           coag_mode_species_number(2) = n_species
       end if
   end do
    
   ASSERT_(any(coag_mode_index /= 0))
   ASSERT_(any(coag_mode_species_number /= 0))


   !  Allocate memory for bufferes
   !  ----------------------------

   coag_species_number_max = maxval(coag_mode_species_number)                 ! max number of species from the coagulating modes

   allocate(coag_intermodal_transfer(coag_mode_species_number(1)), __STAT__)  ! indexes used for the intermodal transfer

   allocate(coag_mode_sigma(n_coag_modes), __STAT__)                          ! geometric standard deviations of the coagulating modes

   allocate(iq_nmr(n_coag_modes),     __STAT__)                               ! number mixing ratio indexes of the coagulating modes 
   allocate(iq_dgn_dry(n_coag_modes), __STAT__)                               ! dry size indexes of the coagulating modes
   allocate(iq_dgn_wet(n_coag_modes), __STAT__)                               ! wet size indexes of the coagulating modes

   allocate(iq_mmr(coag_species_number_max, n_coag_modes),         __STAT__)  ! mass mixing ratio indexes of the coagulating modes
   allocate(density_species(coag_species_number_max, n_coag_modes), __STAT__) ! species bulk densities of the coagulating modes

   allocate(qa_number(n_coag_modes), __STAT__)                                ! number mixing ratios of the coagulating modes
   allocate(qa_mass(coag_species_number_max, n_coag_modes), __STAT__)         ! mass mixing ratios of the coagulating modes

   allocate(diameter_dry(n_coag_modes), __STAT__)                             ! dry aerosol size
   allocate(diameter_wet(n_coag_modes), __STAT__)                             ! wet aerosol size

   allocate(density_dry(n_coag_modes), __STAT__)                              ! dry aerosol density
   allocate(density_wet(n_coag_modes), __STAT__)                              ! wet aerosol density


   coag_intermodal_transfer(:) = (/1, 2, 3, 6/)

   coag_mode_sigma = 0.0
   density_species = 0.0

   iq_nmr = 0
   iq_mmr = 0

   iq_dgn_dry = 0
   iq_dgn_wet = 0


   do n = 1, n_coag_modes

       m = coag_mode_index(n)
       call MAM_AerosolModeGet(self%mode(m), name      = mode_name, &
                                             sigma     = sigma,     &
                                             n_species = n_species)

       ! geometric standard deviation
       coag_mode_sigma(n) = sigma

       ! aerosol dry and wet sizes
       field_name = 'DGN_DRY_' // trim(mode_name)
       iq_dgn_dry(n) = MAPL_SimpleBundleGetIndex(Da, field_name, 3, __RC__)

       field_name = 'DGN_WET_' // trim(mode_name)
       iq_dgn_wet(n) = MAPL_SimpleBundleGetIndex(Da, field_name, 3, __RC__)

       ! number mixing ratio
       field_name = 'NUM_A_' // trim(mode_name)
       iq_nmr(n) = MAPL_SimpleBundleGetIndex(qa, field_name, 3, __RC__)


       do s = 1, n_species
           species_name = self%mode(m)%species(s)%name
           field_name  = trim(species_name) // '_A_' // trim(mode_name)
           iq_mmr(s, n) = MAPL_SimpleBundleGetIndex(qa, field_name, 3, __RC__)

           density_species(s, n) = self%mode(m)%species(s)%component%density
       end do
   end do
  

   ! coagulation and column integrated number flux due to it
   allocate(dz(km), __STAT__)

   do j = 1, jm
       do i = 1, im

           do k = 1, km

               ! mid level pressure, temperature and air density
               pressure    = 0.5 * (ple(i,j,k-1) + ple(i,j,k))
               temperature = T(i,j,k)
               density_air = rhoa(i,j,k)

               ! dry and wet densities 
               do n = 1, n_coag_modes
                   ! number of species in this mode 
                   n_species = coag_mode_species_number(n)

                   ! aerosol size
                   diameter_dry(n) = Da%r3(iq_dgn_dry(n))%q(i,j,k)
                   diameter_wet(n) = Da%r3(iq_dgn_wet(n))%q(i,j,k)

                   ! aerosol species mass mixing ratios
                   qa_number(n) = 0.0 
                   qa_number(n) = qa%r3(iq_nmr(n))%q(i,j,k)

                   qa_mass(:,n) = 0.0
                   do s = 1, n_species
                       qa_mass(s,n) = qa%r3(iq_mmr(s,n))%q(i,j,k)
                   end do

                   ! dry density
                   if (any(qa_mass(:,n) > 1.0e-15)) then
                       density_dry(n) = sum(qa_mass(1:n_species,n)) / sum(qa_mass(1:n_species,n)/density_species(1:n_species,n))
                   else
                       density_dry(n) = sum(density_species(1:n_species,n))/n_species
                   end if

                   ! wet density
                   f_hg = max(diameter_wet(n) / diameter_dry(n), 1.0) ! hygroscopic growth factor
                   density_wet(n) = f_hg**(-3) * density_dry(n) + (1 - f_hg**(-3)) * density_water
               end do

               ! coagulation 
               call MAML_Coagulation(pressure,                 &
                                     temperature,              &
                                     density_air,              &
                                     qa_number,                &
                                     qa_mass,                  &                                       
                                     diameter_wet,             &
                                     density_wet,              &
                                     coag_mode_sigma,          &
                                     coag_mode_species_number, &
                                     coag_intermodal_transfer, &
                                     cdt)

               ! update the number and mass mixing ratios
               do n = 1, n_coag_modes
                   ! number of species in this mode 
                   n_species = coag_mode_species_number(n)  

                   qa%r3(iq_nmr(n))%q(i,j,k) = qa_number(n)
                                 
                   do s = 1, n_species
                       qa%r3(iq_mmr(s,n))%q(i,j,k) = qa_mass(s,n)
                   end do
               end do 

           end do ! k

           dz(:) = delp(i,j,:) / (MAPL_GRAV * rhoa(i,j,:))

       end do ! i
   end do ! j

   !  Free dynamically allocated memory
   !  ---------------------------------
   deallocate(dz, __STAT__)

   deallocate(density_dry, __STAT__)
   deallocate(density_wet, __STAT__)

   deallocate(diameter_dry, __STAT__)
   deallocate(diameter_wet, __STAT__)

   deallocate(qa_number, __STAT__)
   deallocate(qa_mass, __STAT__)

   deallocate(iq_mmr, __STAT__)
   deallocate(density_species, __STAT__)

   deallocate(iq_nmr,     __STAT__)
   deallocate(iq_dgn_dry, __STAT__)
   deallocate(iq_dgn_wet, __STAT__)

   deallocate(coag_mode_sigma, __STAT__)
   deallocate(coag_intermodal_transfer, __STAT__)

   end subroutine MAM_CoagulationBimodal



#if (1)
!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: MAM_Coagulation --- models the coagualtion of interstitial aerosols
!
!
! !INTERFACE:

   subroutine MAM_Coagulation(self, import, export, qa, Da, cdt, rc)
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


! !DESCRIPTION: Coagulation of interstitial aerosols: 
!               - MAM7 scheme models coagulation by treating intramodal 
!               coagulation of Aitken, accumulation and primary carbon (PCM) 
!               modes as well as intramodal coagulation from Aitken to 
!               accumulation, primary carbon to accumulation and aged Aitken 
!               via primary carbon to accumulation mode.
!               - MAM3 scheme models intra- and intermodal coagulation of 
!               the Aitken and accumulation modes.
!
! !REVISION HISTORY:
!
!  09Jan2012  A. Darmenov  First crack.
!
!EOP
!-------------------------------------------------------------------------

              __Iam__('MAM_Coagulation')


   ! Mode parameters
   ! ---------------
   real                      :: sigma
   character(len=MAM_MAXSTR) :: mode_name

   integer                   :: n_species
   character(len=MAM_MAXSTR) :: species_name

   character(len=MAM_MAXSTR) :: field_name

   ! indexes 
   integer :: m, s, n 
   integer :: i, im, j, jm, k, km

   integer :: n_coag_pairs                                              ! number of pairs of modes
   integer :: n_coag_modes                                              ! number of coagulation modes
   integer :: coag_species_number_max                                   ! max number of species

   real    :: temperature                                               ! local temperature
   real    :: pressure                                                  ! local mid-level pressure
   real    :: density_air                                               ! local air density

   real    :: f_hg                                                      ! hygroscopic growth factor

   ! buffers
   integer, allocatable, dimension(  :) :: coag_mode_index              ! indexes of the modes included in the coagulation
   integer, allocatable, dimension(  :) :: coag_mode_species_number     ! number of species in each mode

   integer, allocatable, dimension(:,:) :: coag_intermodal_transfer     ! intermodal coagulation mapping

   integer, allocatable, dimension(  :) :: iq_dgn_dry, iq_dgn_wet       ! dry and wet size bundle indexes

   integer, allocatable, dimension(  :) :: iq_nmr                       ! number mixing ratio bundle indexes
   integer, allocatable, dimension(:,:) :: iq_mmr                       ! mass mixing ratio bundle indexes


   real,    allocatable, dimension(  :) :: qa_number                    ! number mixing ratios
   real,    allocatable, dimension(:,:) :: qa_mass                      ! mass mixing ratios
   real,    allocatable, dimension(:,:) :: density_species              ! species bulk densities

   real,    allocatable, dimension(  :) :: diameter_dry, diameter_wet   ! dry and wet aerosol sizes
   real,    allocatable, dimension(  :) :: density_dry,  density_wet    ! dry and wet aerosol densities
   real,    allocatable, dimension(  :) :: coag_mode_sigma              ! geometric standard deviations

   real,    allocatable, dimension(  :) :: mass2vol_aitken_age          ! mass to volume factor
   real,    allocatable, dimension(  :) :: mass2vol_pcarbon             ! mass to volume factor
   

   ! other derived variables
   real,    allocatable, dimension(:)   :: dz


   !  Input fields from fvGCM
   !  -----------------------
   real, pointer, dimension(:,:,:) :: rhoa
   real, pointer, dimension(:,:,:) :: delp
   real, pointer, dimension(:,:,:) :: ple
   real, pointer, dimension(:,:,:) :: T

   !  Exports
   !  -----------------------
   real, pointer, dimension(:,:)   :: flux
  

   real, parameter :: density_water = 1000.0 ! kg m-3

   ! MAM7
   integer, parameter :: mam7_mode_ait = 1               ! Aitken mode index
   integer, parameter :: mam7_mode_pcm = 2               ! primary carbon mode index
   integer, parameter :: mam7_mode_acc = 3               ! accumulation mode index

   integer, parameter :: mam7_ait_acc  = 1               ! AIT -> ACC coagulation pair index
   integer, parameter :: mam7_pcm_acc  = 2               ! PCM -> ACC coagulation pair inde
   integer, parameter :: mam7_ait_pcm  = 3               ! AIT -> PCM + 'ageing' -> ACC coagulation pair index

   ! MAM3
   integer, parameter :: mam3_mode_ait = 1               ! Aitken mode index
   integer, parameter :: mam3_mode_acc = 2               ! accumulation mode index

   integer, parameter :: mam3_ait_acc = 1                ! AIT -> ACC coagulation pair index



   !  Get Imports
   !  --------------
   call MAPL_GetPointer(import, rhoa, 'AIRDENS', __RC__)
   call MAPL_GetPointer(import, delp, 'DELP',    __RC__)
   call MAPL_GetPointer(import, ple,  'PLE',     __RC__)
   call MAPL_GetPointer(import, T,    'T',       __RC__)


   !  Local dimensions
   !  ----------------
   im = size(rhoa, 1)
   jm = size(rhoa, 2)
   km = size(rhoa, 3)


   !  Find the coagulating modes
   !  --------------------------
   if (self%id == MAM7_SCHEME) then
       n_coag_modes = 3                    ! AIT, PCM, ACC
       n_coag_pairs = 3                    ! AIT->ACC, PCM->ACC, AIT->PCM
   else if (self%id == MAM3_SCHEME) then
       n_coag_modes = 2                    ! AIT, ACC
       n_coag_pairs = 1                    ! AIT->ACC
   else
       __raise__ (MAM_UNKNOWN_SCHEME_ERROR, "Unsupported MAM model.")
   end if


   allocate(coag_mode_index(n_coag_modes),          __STAT__)
   allocate(coag_mode_species_number(n_coag_modes), __STAT__)

   coag_mode_index          = 0
   coag_mode_species_number = 0

   do m = 1, self%n_modes
       call MAM_AerosolModeGet(self%mode(m), name    = mode_name, &
                                           n_species = n_species  ) 

       if (self%id == MAM7_SCHEME) then
           if (trim(mode_name) == trim(MAM7_AITKEN_MODE_NAME)) then
               coag_mode_index(mam7_mode_ait) = m
               coag_mode_species_number(mam7_mode_ait) = n_species
           else if (trim(mode_name) == trim(MAM7_PRIMARY_CARBON_MODE_NAME)) then
               coag_mode_index(mam7_mode_pcm) = m
               coag_mode_species_number(mam7_mode_pcm) = n_species
           else if (trim(mode_name) == trim(MAM7_ACCUMULATION_MODE_NAME)) then
               coag_mode_index(mam7_mode_acc) = m
               coag_mode_species_number(mam7_mode_acc) = n_species
           end if
       else if (self%id == MAM3_SCHEME) then
           if (trim(mode_name) == trim(MAM3_AITKEN_MODE_NAME)) then
               coag_mode_index(mam3_mode_ait) = m
               coag_mode_species_number(mam3_mode_ait) = n_species
           else if (trim(mode_name) == trim(MAM3_ACCUMULATION_MODE_NAME)) then
               coag_mode_index(mam3_mode_acc) = m
               coag_mode_species_number(mam3_mode_acc) = n_species
           end if
       end if
   end do

   ASSERT_(any(coag_mode_index /= 0))
   ASSERT_(any(coag_mode_species_number /= 0))


   !  Allocate memory for bufferes
   !  ----------------------------

   coag_species_number_max = maxval(coag_mode_species_number)                 ! max number of species from the coagulating modes

   allocate(coag_intermodal_transfer(coag_species_number_max,coag_mode_species_number(1)), __STAT__) ! indexes used for the intermodal transfer

   allocate(coag_mode_sigma(n_coag_modes), __STAT__)                          ! geometric standard deviations of the coagulating modes

   allocate(iq_nmr(n_coag_modes),     __STAT__)                               ! number mixing ratio indexes of the coagulating modes 
   allocate(iq_dgn_dry(n_coag_modes), __STAT__)                               ! dry size indexes of the coagulating modes
   allocate(iq_dgn_wet(n_coag_modes), __STAT__)                               ! wet size indexes of the coagulating modes

   allocate(iq_mmr(coag_species_number_max, n_coag_modes),         __STAT__)  ! mass mixing ratio indexes of the coagulating modes
   allocate(density_species(coag_species_number_max, n_coag_modes), __STAT__) ! species bulk densities of the coagulating modes

   allocate(qa_number(n_coag_modes), __STAT__)                                ! number mixing ratios of the coagulating modes
   allocate(qa_mass(coag_species_number_max, n_coag_modes), __STAT__)         ! mass mixing ratios of the coagulating modes

   allocate(diameter_dry(n_coag_modes), __STAT__)                             ! dry aerosol size
   allocate(diameter_wet(n_coag_modes), __STAT__)                             ! wet aerosol size

   allocate(density_dry(n_coag_modes), __STAT__)                              ! dry aerosol density
   allocate(density_wet(n_coag_modes), __STAT__)                              ! wet aerosol density
   
   allocate(mass2vol_aitken_age(coag_species_number_max), __STAT__)           ! mass to volume factor, needed only for MAM7 coagulation
   allocate(mass2vol_pcarbon(coag_species_number_max),    __STAT__)           ! mass to volume factor, needed only for MAM7 coagulation


   coag_intermodal_transfer(:,:) = 0

   if (self%id == MAM7_SCHEME) then
       n_species = coag_mode_species_number(mam7_mode_ait)
       coag_intermodal_transfer(mam7_ait_acc, 1:n_species) = (/1, 2, 3, 6/)

       n_species = coag_mode_species_number(mam7_mode_pcm)
       coag_intermodal_transfer(mam7_pcm_acc, 1:n_species) = (/4, 5/)

       n_species = coag_mode_species_number(mam7_mode_ait)
       coag_intermodal_transfer(mam7_ait_pcm, 1:n_species) = (/0, 0, 0, 0/)

   else if (self%id == MAM3_SCHEME) then
       n_species = coag_mode_species_number(mam3_mode_ait)
       coag_intermodal_transfer(mam3_ait_acc, 1:n_species) = (/1, 2, 6/)
   end if


   coag_mode_sigma = 0.0
   density_species = 0.0

   iq_nmr = 0
   iq_mmr = 0

   iq_dgn_dry = 0
   iq_dgn_wet = 0


   do n = 1, n_coag_modes

       m = coag_mode_index(n)
       call MAM_AerosolModeGet(self%mode(m), name      = mode_name, &
                                             sigma     = sigma,     &
                                             n_species = n_species)

       ! geometric standard deviation
       coag_mode_sigma(n) = sigma

       ! aerosol dry and wet sizes
       field_name = 'DGN_DRY_' // trim(mode_name)
       iq_dgn_dry(n) = MAPL_SimpleBundleGetIndex(Da, field_name, 3, __RC__)

       field_name = 'DGN_WET_' // trim(mode_name)
       iq_dgn_wet(n) = MAPL_SimpleBundleGetIndex(Da, field_name, 3, __RC__)

       ! number mixing ratio
       field_name = 'NUM_A_' // trim(mode_name)
       iq_nmr(n) = MAPL_SimpleBundleGetIndex(qa, field_name, 3, __RC__)


       do s = 1, n_species
           species_name = self%mode(m)%species(s)%name
           field_name  = trim(species_name) // '_A_' // trim(mode_name)
           iq_mmr(s, n) = MAPL_SimpleBundleGetIndex(qa, field_name, 3, __RC__)

           density_species(s, n) = self%mode(m)%species(s)%component%density
       end do
   end do
  

   ! coagulation calculation and column integrated number flux due to it
   allocate(dz(km), __STAT__)

   do j = 1, jm
       do i = 1, im

           do k = 1, km

               ! mid level pressure, temperature and air density
               pressure    = (ple(i,j,k-1) + ple(i,j,k)) * 0.5
               temperature = T(i,j,k)
               density_air = rhoa(i,j,k)

               ! dry and wet densities 
               do n = 1, n_coag_modes
                   ! number of species in this mode 
                   n_species = coag_mode_species_number(n)

                   ! aerosol size
                   diameter_dry(n) = Da%r3(iq_dgn_dry(n))%q(i,j,k)
                   diameter_wet(n) = Da%r3(iq_dgn_wet(n))%q(i,j,k)

                   ! aerosol species mass mixing ratios
                   qa_number(n) = 0.0 
                   qa_number(n) = qa%r3(iq_nmr(n))%q(i,j,k)

                   qa_mass(:,n) = 0.0
                   do s = 1, n_species
                       qa_mass(s,n) = qa%r3(iq_mmr(s,n))%q(i,j,k)
                   end do

                   ! dry density
                   if (any(qa_mass(:,n) > 1.0e-15)) then
                       density_dry(n) = sum(qa_mass(1:n_species,n)) / &
                                        (sum(qa_mass(1:n_species,n) / density_species(1:n_species,n)))
                   else
                       density_dry(n) = sum(density_species(1:n_species,n)) / n_species
                   end if

                   ! wet density
                   f_hg = max(diameter_wet(n) / diameter_dry(n), 1.0) ! hygroscopic growth factor
                   density_wet(n) = f_hg**(-3) * density_dry(n) + (1 - f_hg**(-3)) * density_water
               end do

               ! coagulation 
               call MAML_Coagulation(pressure,                 &
                                     temperature,              &
                                     density_air,              &
                                     qa_number,                &
                                     qa_mass,                  &
                                     diameter_dry,             &
                                     diameter_wet,             &
                                     density_wet,              &
                                     coag_mode_sigma,          &
                                     coag_mode_species_number, &
                                     coag_intermodal_transfer, &
                                     mass2vol_aitken_age,      &
                                     mass2vol_pcarbon,         &
                                     cdt)

               ! update the number and mass mixing ratios
               do n = 1, n_coag_modes
                   ! number of species in this mode 
                   n_species = coag_mode_species_number(n)  

                   qa%r3(iq_nmr(n))%q(i,j,k) = qa_number(n)
                                 
                   do s = 1, n_species
                       qa%r3(iq_mmr(s,n))%q(i,j,k) = qa_mass(s,n)
                   end do
               end do 

           end do ! k

           dz(:) = delp(i,j,:) / (MAPL_GRAV * rhoa(i,j,:))

       end do ! i
   end do ! j


   !  Free dynamically allocated memory
   !  ---------------------------------
   deallocate(dz, __STAT__)

   deallocate(mass2vol_aitken_age, __STAT__)
   deallocate(mass2vol_pcarbon, __STAT__)

   deallocate(density_dry, __STAT__)
   deallocate(density_wet, __STAT__)

   deallocate(diameter_dry, __STAT__)
   deallocate(diameter_wet, __STAT__)

   deallocate(qa_number, __STAT__)
   deallocate(qa_mass, __STAT__)

   deallocate(iq_mmr, __STAT__)
   deallocate(density_species, __STAT__)

   deallocate(iq_nmr, __STAT__)
   deallocate(iq_dgn_dry, __STAT__)
   deallocate(iq_dgn_wet, __STAT__)

   deallocate(coag_intermodal_transfer, __STAT__)

   deallocate(coag_mode_sigma, __STAT__) 

   deallocate(coag_mode_index, __STAT__)
   deallocate(coag_mode_species_number, __STAT__)

   end subroutine MAM_Coagulation

#endif

end module MAM_CoagulationMod
