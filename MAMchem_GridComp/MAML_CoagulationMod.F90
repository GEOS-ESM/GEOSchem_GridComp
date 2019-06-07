#include "MAPL_Exceptions.h"
#include "MAPL_Generic.h"

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: MAML_CoagulationMod - Coagulation of aerosol particles.
!
! !INTERFACE:
!
   module MAML_CoagulationMod
!
! !USES:
!
   use MAPL_Mod
   use MAPL_ConstantsMod, only : MAPL_PI, MAPL_RHOWTR, r8 => MAPL_R8, r4 => MAPL_R4

   use modal_aero_coag, only : getcoags_wrapper_f

   use MAM_ComponentsDataMod, only : MAM_SOA_COMPONENT_HYGROSCOPICITY, &
                                     MAM_SULFATE_COMPONENT_HYGROSCOPICITY


   implicit NONE
   private

!
! !PUBLIC MEMBER FUNCTIONS:

   public MAML_Coagulation


! !PRIVATE PARAMETERS 
   real, private, parameter :: pi            = MAPL_PI

   real, private, parameter :: R_univ        = MAPL_RUNIV     ! Universal gas constant,  'J K-1 Kmole-1'
   real, private, parameter :: density_water = MAPL_RHOWTR    ! density of water,  'kg m-3'

   
   ! the fSOA_EquivSO4 factor converts an SOA volume to a volume of SO4(+NH4) 
   ! having same hygroscopicity as the SOA
   real, private, parameter :: fSOA_EquivSO4 = (MAM_SOA_COMPONENT_HYGROSCOPICITY /  &
                                                MAM_SULFATE_COMPONENT_HYGROSCOPICITY)

   ! number of SO4(+NH4) monolayers needed to 'age' a carbon particle
   real, private, parameter :: NUMBER_SO4_MONOLAYERS_PCAGE = 3.0


!
! !DESCRIPTION: 
!
!  {\tt MAML\_CoagulationMod} provides a collection of methods to calculate
!  intra- and intermodal coagulation rates.
!
!
! !REVISION HISTORY:
!
!  03Jan2012  A. Darmenov  Initial version -- based on CESM-1.0.3 CAM/MAM
!                                             modal_aero_coag module
!                                            
!
!EOP
!-------------------------------------------------------------------------


   interface MAML_Coagulation
       module procedure MAML_CoagulationBimodal
       module procedure MAML_Coagulation_AIT_PCM_ACC
   end interface MAML_Coagulation


   contains


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: MAML_Coagulation_Bimodal --- 
!
! !INTERFACE:

   subroutine MAML_CoagulationBimodal(pressure,            &
                                      temperature,         &
                                      density_air,         &
                                      q_number,            &
                                      q_mass,              &
                                      Dg_wet,              &
                                      density_wet,         &
                                      sigma,               &
                                      n_species,           &
                                      intermodal_transfer, &
                                      dt)
                                       

! !USES:

   implicit NONE

! !INPUT/OUTPUT PARAMETERS:
   real, dimension(:),   intent(inout) :: q_number             ! number mixing ratios of the aerosol modes
   real, dimension(:,:), intent(inout) :: q_mass               ! mass mixing ratios of the components in the three modes

! !INPUT PARAMETERS:
   real, intent(in)                    :: pressure             ! pressure at mid level
   real, intent(in)                    :: temperature          ! temperature at mid level
   real, intent(in)                    :: density_air          ! air density

   real, intent(in)                    :: dt                   ! time step

   real, dimension(:), intent(in)      :: Dg_wet               ! wet geometric mean diameter of number size distribution

   real, dimension(:), intent(in)      :: density_wet          ! wet density
   real, dimension(:), intent(in)      :: sigma                ! geometric standard deviation

   integer, dimension(:), intent(in)   :: n_species            ! number of species

   integer, dimension(:), intent(in)   :: intermodal_transfer  ! maps the indexes of the source mode species to the 
                                                               ! indexes of the receiving mode species

   

! !OUTPUT PARAMETERS:


! !DESCRIPTION: Bimodal (e.g., Aitken and accumulation modes) coagulation.
!
! !REVISION HISTORY:
!
!  03Jan2011  A. Darmenov  First crack -- based on modal_aero_coag_sub(),
!                          from CESM-1.0.3
!
!EOP
!-------------------------------------------------------------------------
                   __Iam__('MAML_CoagulationBimodal')
       
   ! local parameters
   integer, parameter :: n_coag_modes = 2                ! number of modes considered in the coagulation process
                                                         ! bimodal coagulation, e.g., 1(AIT) + 1(ACC) = 2

   integer, parameter :: mode_i = 1                      ! source (e.g., Aitken) mode index
   integer, parameter :: mode_j = 2                      ! receiving (e.g., accumulation) mode index


   ! local variables
   real(r8) :: P, T                                      ! pressure and temperarute

   real(r8) :: D_wet_i, density_wet_i                    ! wet size and density of the source mode
   real(r8) :: sigma_i, ln_sigma_i                       ! geometric standard deviation and its log

   real(r8) :: D_wet_j, density_wet_j                    ! wet size and density of the receiving mode
   real(r8) :: sigma_j, ln_sigma_j                       ! geometric standard deviation and its log

   real(r8) :: beta_ii0, beta_ii2                        ! intramodal coagulation rates - 'i' mode
   real(r8) :: beta_jj0, beta_jj2                        ! intramodal coagulation rates - 'j' mode
   real(r8) :: beta_ij0, beta_ij2i, beta_ij2j, beta_ij3  ! intermodal coagulation rates
   
   real(r8), dimension(n_coag_modes) :: number_conc      ! initial number concentrations of 'i' and 'j' modes
   real(r8), dimension(n_coag_modes) :: number_conc_new  ! final   number concentrations of 'i' and 'j' modes
   real(r8), dimension(n_coag_modes) :: number_conc_avg  ! average number concentrations
   
   real(r8) :: N_i, N_j                                  ! temporary variables - initial number concentration  
   real(r8) :: tmp_A, tmp_B, tmp_C, tmp_F, tmp_G, tmp_H  ! temporary variables - various terms 

   real(r8) :: frac_transfer_vol, frac_transfer_vol_max  ! fraction of volume that can be transfered between the modes 
   real(r8) :: vol_loss_i, mass_transfer

   integer  :: m, iq, iq_mode_i, iq_mode_j



   do m = 1, n_coag_modes
       number_conc(m) = q_number(m)*density_air
       number_conc(m) = max(0.0, number_conc(m))
   end do


   ! Calculate the coagulation rates -- use double precision
   !                                    where it is required
   ! --------------------------------------------------------
   P = dble(pressure)
   T = dble(temperature)

   D_wet_i = dble(Dg_wet(mode_i))
   D_wet_j = dble(Dg_wet(mode_j))

   density_wet_i = dble(density_wet(mode_i))
   density_wet_j = dble(density_wet(mode_j))

   sigma_i = dble(sigma(mode_i))
   sigma_j = dble(sigma(mode_j))

   ln_sigma_i = log(sigma_i)
   ln_sigma_j = log(sigma_j)

   ! coagulation rates using CMAQ 'fast' method, based on Whitby's 
   ! approximation approach
   call getcoags_wrapper_f(T, P, D_wet_i,       &
                                 D_wet_j,       &
                                 sigma_i,       &
                                 sigma_j,       &
                                 ln_sigma_i,    &
                                 ln_sigma_j,    &
                                 density_wet_i, &
                                 density_wet_j, &
                                 beta_ij0,      &
                                 beta_ij2i,     &
                                 beta_ij2j,     &
                                 beta_ij3,      &
                                 beta_ii0,      &
                                 beta_ii2,      &
                                 beta_jj0,      &
                                 beta_jj2       )


   ! Compute number mixing ratio changes due to 
   ! coagulation between 'i' - source/from mode and 
   ! 'j' - receiving/to mode:
   !
   !                    intramodal          intermodal      
   !                -------------------   ------------------
   !    | dN_i/dt = -beta_ii0 * N_i*N_i - beta_ij0 * N_i*N_j
   !    | dN_j/dt = -beta_jj0 * N_j*N_j
   !
   !    | dS_i/dt = -beta_ii2 * S_i*N_i - beta_ij2i * S_i*N_j
   !    | dS_j/dt = -beta_jj2 * S_j*N_j + beta_ij2j * S_i*N_j
   ! 
   !    | dV_i/dt = -beta_ij3 * V_i*N_j
   !    | dV_j/dt = -dV_i/dt
   !
   ! 
   ! The first system of equations is solved for N_j first, and 
   ! then for N_i, assuming that the coag. coefficients are 
   ! constants during the integration step, and by substituting N_j
   ! with its mean value <N_j> = 1/2 * (N_j(t) + N_j(t+dt)), i.e.
   !
   !      dN_j/dt = -beta_jj0 * N_j*N_j
   !      dN_i/dt = -beta_ii0 * N_i*N_i - (beta_ij0 * <N_j>)*N_j
   ! 
   ! --------------------------------------------------------

   ! update number mixing ratio of the receiving (j) mode
   N_j = number_conc(mode_j)                                        ! N_j(t)
   number_conc_new(mode_j) = N_j / (1.0 + beta_jj0*N_j*dt)          ! N_j(t + dt) 
   number_conc_avg(mode_j) = 0.5 * (number_conc_new(mode_j) + N_j)  ! [N_j(t) + N_j(t + dt)]/2
 
   q_number(mode_j) = number_conc_new(mode_j) / density_air         ! update the input number concentration
   
   ! update number mixing ratio of the source (i) mode
   N_i = number_conc(mode_i)                                        ! N_i(t)

   tmp_A = beta_ij0 * number_conc_avg(mode_j) * dt                  ! recurring terms
   tmp_B = beta_ii0 * dt                                            ! ...
   tmp_C = tmp_A + (tmp_B * N_i)                                    ! ...

   if (abs(tmp_C) < 1e-2) then
       number_conc_new(mode_i) = N_i * exp(-tmp_C)                  ! N_i(t + dt)
   else if (abs(tmp_A) < 1e-3) then
       number_conc_new(mode_i) = exp(-tmp_A) * N_i/(1.0 + tmp_B*N_i)! N_i(t + dt)
   else
       tmp_F = (tmp_B * N_i)/tmp_C                                  ! recurring terms
       tmp_G = exp(-tmp_A)                                          ! ...
       tmp_H = tmp_G*(1.0 - tmp_F)/(1.0 - tmp_G*tmp_F)              ! ...

       number_conc_new(mode_i) = N_i * max(0.0, min(1.0, tmp_H))    ! N_i(t + dt)
   end if

   number_conc_avg(mode_i) = 0.5*(number_conc_new(mode_i) + N_i)    ! [N_i(t) + N_i(t + dt)]/2

   q_number(mode_i) = number_conc_new(mode_i) / density_air


   ! Compute mass mixing ratios changes due to coagulation between 
   ! source and receiving modes
   !  ------------------------------------------------------------

   ! maximum fraction of transfered volume = 1 - eps
   frac_transfer_vol_max = 1.0 - 1.0e1*epsilon(1.0_r8)

   ! first order loss rate for mode 'i' volume
   vol_loss_i = beta_ij3 * number_conc_avg(mode_j)

   ! fraction of 'i' volume transferred to 'j' over time dt
   frac_transfer_vol = 1.0 - exp(-vol_loss_i*dt)

   frac_transfer_vol = min(frac_transfer_vol_max, frac_transfer_vol)
   frac_transfer_vol = max(0.0, frac_transfer_vol)

   do iq = 1, n_species(mode_i)
        iq_mode_i = iq
        iq_mode_j = intermodal_transfer(iq)

       if (iq_mode_j > 0) then
           ! species mass transfered from 'i' to 'j' mode
           mass_transfer = q_mass(iq_mode_i, mode_i)*frac_transfer_vol

           q_mass(iq_mode_i, mode_i) = q_mass(iq_mode_i, mode_i) - mass_transfer
           q_mass(iq_mode_j, mode_j) = q_mass(iq_mode_j, mode_j) + mass_transfer
       end if
   end do
 
   end subroutine MAML_CoagulationBimodal



!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: MAML_Coagulation_AIT_PCM_ACC --- 
!
! !INTERFACE:

   subroutine MAML_Coagulation_AIT_PCM_ACC(pressure,            &
                                           temperature,         &
                                           density_air,         &
                                           q_number,            &
                                           q_mass,              &
                                           Dg_dry,              &
                                           Dg_wet,              &
                                           density_wet,         &
                                           sigma,               &
                                           n_species,           &
                                           intermodal_transfer, &
                                           mass2vol_aitken_age, &
                                           mass2vol_pcarbon,    &
                                           dt)
                                       

! !USES:

   implicit NONE

! !INPUT/OUTPUT PARAMETERS:
   real, dimension(:),   intent(inout) :: q_number             ! number mixing ratios of the  and ACC modes
   real, dimension(:,:), intent(inout) :: q_mass               ! mass mixing ratios of the components in the three modes

! !INPUT PARAMETERS:
   real, intent(in)                    :: pressure             ! pressure at mid level
   real, intent(in)                    :: temperature          ! temperature at mid level
   real, intent(in)                    :: density_air          ! air density

   real, intent(in)                    :: dt                   ! time step

   real, dimension(:), intent(in)      :: Dg_dry               ! dry geometric mean diameter of number size distribution
   real, dimension(:), intent(in)      :: Dg_wet               ! wet geometric mean diameter of number size distribution

   real, dimension(:), intent(in)      :: density_wet          ! wet density
   real, dimension(:), intent(in)      :: sigma                ! geometric standard deviation

   integer, dimension(:), intent(in)   :: n_species            ! number of species

   integer, dimension(:,:), intent(in) :: intermodal_transfer  ! maps the indexes of the source mode species to the 
                                                               ! indexes of the receiving mode species

   real, dimension(:), intent(in)      :: mass2vol_aitken_age  ! conversion factor for the aitken to PCM aging, equal
                                                               ! to 
   real, dimension(:), intent(in)      :: mass2vol_pcarbon     !


! !OUTPUT PARAMETERS:


! !DESCRIPTION: Coagulation of Aitken, primary carbon and accumulation  
!               modes.
!
! !REVISION HISTORY:
!
!  13Jan2011  A. Darmenov  First crack -- based on modal_aero_coag_sub(),
!                          from CESM-1.0.3
!
!EOP
!-------------------------------------------------------------------------
                   __Iam__('MAML_Coagulation_AIT_PCM_ACC')
       
   ! local parameters
   integer, parameter :: n_coag_modes = 3                ! number of modes considered in the coagulation process
                                                         ! 1(AIT) + 1(PCM) + 1(ACC) = 3

   integer, parameter :: n_coag_pairs = 3                ! number of mode pairs that coagulate 
                                                         !    1) AIT -> ACC
                                                         !    2) PCM -> ACC
                                                         !    3) AIT -> PCM + 'ageing' -> ACC

   integer, parameter :: mode_ait = 1                    ! Aitken mode index
   integer, parameter :: mode_pcm = 2                    ! primary carbon mode index
   integer, parameter :: mode_acc = 3                    ! accumulation mode index

   integer, parameter :: ait_acc = 1                     ! AIT -> ACC coagulation pair index
   integer, parameter :: pcm_acc = 2                     ! PCM -> ACC coagulation pair inde
   integer, parameter :: ait_pcm = 3                     ! AIT -> PCM + 'ageing' -> ACC coagulation pair index

   integer, parameter :: coag_mode_source(n_coag_pairs) = (/mode_ait, mode_pcm, mode_ait/)
   integer, parameter :: coag_mode_receiv(n_coag_pairs) = (/mode_acc, mode_acc, mode_pcm/)


   ! local variables 
   real(r8) :: P, T, air_conc                            ! air molar density

   real(r8) :: D_wet_i, density_wet_i                    ! wet size and density of the source mode
   real(r8) :: sigma_i, ln_sigma_i                       ! geometric standard deviation and its log

   real(r8) :: D_wet_j, density_wet_j                    ! wet size and density of the receiving mode
   real(r8) :: sigma_j, ln_sigma_j                       ! geometric standard deviation and its log

   real(r8), dimension(n_coag_pairs) :: beta_ii0         ! intramodal coagulation rates
   real(r8), dimension(n_coag_pairs) :: beta_ii2         ! ...
   real(r8), dimension(n_coag_pairs) :: beta_jj0         ! ...
   real(r8), dimension(n_coag_pairs) :: beta_jj2         ! ...
   real(r8), dimension(n_coag_pairs) :: beta_ij0         ! intermodal coagulation rates
   real(r8), dimension(n_coag_pairs) :: beta_ij2i        ! ...
   real(r8), dimension(n_coag_pairs) :: beta_ij2j        ! ...
   real(r8), dimension(n_coag_pairs) :: beta_ij3         ! ...
   
   real(r8), dimension(n_coag_modes) :: number_conc      ! initial number concentrations of the modes
   real(r8), dimension(n_coag_modes) :: number_conc_new  ! final   number concentrations of the modes
   real(r8), dimension(n_coag_modes) :: number_conc_avg  ! mean (during the time step) number concentrations
   
   real(r8) :: N_0                                       ! temporary variables - initial number concentration
   real(r8) :: tmp_A, tmp_B, tmp_C, tmp_F, tmp_G, tmp_H  ! temporary variables - various terms 

   real(r8) :: frac_transfer_vol, frac_transfer_vol_max  ! fraction of volume that can be transfered between the modes 
   real(r8) :: frac_transfer_vol_pcage                   ! fraction of volume transfered between the modes due to primary carbone aging
   real(r8) :: vol_shell, vol_core                       ! volumes of shell and core 
   real(r8) :: vol_loss, mass_transfer, number_transfer  ! volume, mass and number transfered between the modes

   real(r8) :: dR_so4_monolayers_pcage                   ! change in size(radius) due to number of SO4 monolayers

   real     :: f_vol2sfc_pcarbon                         ! volume to surface factor

   integer  :: m, n
   integer  :: iq, iq_mode_i, iq_mode_j
   integer  :: mode_i, mode_j
   integer  :: coag_pair


   ! air molar density (kmol m-3)
   P = dble(pressure)
   T = dble(temperature)

   do m = 1, n_coag_modes
       number_conc(m) = q_number(m)*density_air
       number_conc(m) = max(0.0, number_conc(m))
   end do


   ! Calculate the coagulation rates -- use double precision
   !                                    where it is required
   ! --------------------------------------------------------
   beta_ij0  = 0.0
   beta_ij2i = 0.0 
   beta_ij2j = 0.0
   beta_ij3  = 0.0 
   beta_ii0  = 0.0
   beta_ii2  = 0.0
   beta_jj0  = 0.0
   beta_jj2  = 0.0

   do n = 1, n_coag_pairs
       mode_i = coag_mode_source(n)
       mode_j = coag_mode_receiv(n)

       D_wet_i = dble(Dg_wet(mode_i))
       D_wet_j = dble(Dg_wet(mode_j))

       density_wet_i = dble(density_wet(mode_i))
       density_wet_j = dble(density_wet(mode_j))

       sigma_i = dble(sigma(mode_i))
       sigma_j = dble(sigma(mode_j))

       ln_sigma_i = log(sigma_i)
       ln_sigma_j = log(sigma_j)

       ! coagulation rates using CMAQ 'fast' method, based on Whitby's 
       ! approximation approach
       call getcoags_wrapper_f(T, P, D_wet_i,       &
                                     D_wet_j,       &
                                     sigma_i,       &
                                     sigma_j,       &
                                     ln_sigma_i,    &
                                     ln_sigma_j,    &
                                     density_wet_i, &
                                     density_wet_j, &
                                     beta_ij0(n),   &
                                     beta_ij2i(n),  &
                                     beta_ij2j(n),  &
                                     beta_ij3(n),   &
                                     beta_ii0(n),   &
                                     beta_ii2(n),   &
                                     beta_jj0(n),   &
                                     beta_jj2(n)    )
   end do


   ! Compute number mixing ratio changes due to 
   ! coagulation between ait, primary carbon and  
   ! accumulation mode
   !
   !                    intramodal          intermodal      
   !                -------------------   ------------------
   !    | dN_ait/dt = -beta_ii0 * N_ait*N_ait - beta_ij0 * N_acc*N_ait - beta_ij0 * N_pcm*N_ait
   !    | dN_pcm/dt = -beta_ii0 * N_pcm*N_pcm - beta_ij0 * N_acc*N_pcm
   !    | dN_acc/dt = -beta_jj0 * N_acc*N_acc
   !
   !    | dV_i/dt = -beta_ij3 * V_i*N_j
   !    | dV_j/dt = -dV_i/dt
   !
   ! 
   ! The first system of equations is solved for N_acc first, and 
   ! then for N_pcm, and finally for N_ait assuming that the coag. 
   ! coefficients are constants during the integration step, and 
   ! by substituting N_ait and N_pcm with their mean values
   ! <N_acc|pcm> = 1/2 * (N_acc|pcm(t) + N_acc|pcm(t+dt)), i.e.
   !
   !      dN_acc/dt = -beta_jj0 * N_acc*N_acc
   !      dN_pcm/dt = -beta_ii0 * N_pcm*N_pcm - (beta_ij0 * <N_acc>)*N_pcm
   !      dN_ait/dt = -beta_jj0 * N_ait*N_ait - (beta_ij0 * <N_acc> + beta_ij0 * <N_pcm>)*N_ait
   !
   !
   ! TODO: General coagulation solver. Assuming that the coagulating modes
   !       are ordered by size (from smaller to larger), calculate the coagulation
   !       rates for every pair and find the intermodal terms
   !       <beta * N>_i = sum {j>i} (beta_ij0 * <N_j>)
   !
   !       Then solve the equations
   !       dN_j/dt = -beta_jj0 * N_j*N_j,                      larges mode
   !       <N_j>   = ...
   !       ...
   !       dN_i/dt = -beta_jj0 * N_i*N_i - <beta * N>_i * N_i, i = j -1
   !       ...
   !       dN_i/dt = -beta_jj0 * N_i*N_i - <beta * N>_i * N_i, i = 1
   !       ...
   ! --------------------------------------------------------

   
   ! update number mixing ratio of the accumulation mode
   N_0 = number_conc(mode_acc)                                            ! N_acc(t)
   number_conc_new(mode_acc) = N_0 / (1.0 + beta_jj0(ait_acc)*N_0*dt)     ! N_acc(t + dt) 
   number_conc_avg(mode_acc) = 0.5 * (number_conc_new(mode_acc) + N_0)    ! [N_acc(t) + N_acc(t + dt)]/2
 
   q_number(mode_acc) = number_conc_new(mode_acc) / air_conc              ! update the input number concentration


   ! update number mixing ratio of the primary carbon mode
   N_0 = number_conc(mode_pcm)                                            ! N_pcm(t)

   tmp_A = beta_ij0(pcm_acc) * number_conc_avg(mode_acc) * dt             ! recurring terms
   tmp_B = beta_ii0(pcm_acc) * dt                                         ! ...
   tmp_C = tmp_A + (tmp_B * N_0)                                          ! ...

   if (abs(tmp_C) < 1e-2) then
       number_conc_new(mode_pcm) = N_0 * exp(-tmp_C)                      ! N_pcm(t + dt)
   else if (abs(tmp_A) < 1e-3) then
       number_conc_new(mode_pcm) = exp(-tmp_A) * N_0/(1.0 + tmp_B*N_0)    ! N_pcm(t + dt)
   else
       tmp_F = (tmp_B * N_0)/tmp_C                                        ! recurring terms
       tmp_G = exp(-tmp_A)                                                ! ...
       tmp_H = tmp_G*(1.0 - tmp_F)/(1.0 - tmp_G*tmp_F)                    ! ...

       number_conc_new(mode_pcm) = N_0 * max(0.0, min(1.0, tmp_H))        ! N_pcm(t + dt)
   end if

   number_conc_avg(mode_pcm) = 0.5*(number_conc_new(mode_pcm) + N_0)      ! [N_pcm(t) + N_pcm(t + dt)]/2

   q_number(mode_pcm) = number_conc_new(mode_pcm) / density_air


   ! update number mixing ratio of the aitken mode
   !
   ! coagulation pair: aitken -> primary carbon
   !
   N_0 = number_conc(mode_ait)                                            ! N_ait(t)

   tmp_A = ( beta_ij0(ait_acc) * number_conc_avg(mode_acc) + &
             beta_ij0(ait_pcm) * number_conc_avg(mode_pcm) ) * dt         ! recurring terms
   tmp_B = beta_ii0(ait_acc) * dt                                         ! ...
   tmp_C = tmp_A + (tmp_B * N_0)                                          ! ...

   if (abs(tmp_C) < 1e-2) then
       number_conc_new(mode_ait) = N_0 * exp(-tmp_C)                      ! N_ait(t + dt)
   else if (abs(tmp_A) < 1e-3) then
       number_conc_new(mode_ait) = exp(-tmp_A) * N_0/(1.0 + tmp_B*N_0)    ! N_ait(t + dt)
   else
       tmp_F = (tmp_B * N_0)/tmp_C                                        ! recurring terms
       tmp_G = exp(-tmp_A)                                                ! ...
       tmp_H = tmp_G*(1.0 - tmp_F)/(1.0 - tmp_G*tmp_F)                    ! ...

       number_conc_new(mode_ait) = N_0 * max(0.0, min(1.0, tmp_H))        ! N_ait(t + dt)
   end if

   number_conc_avg(mode_ait) = 0.5*(number_conc_new(mode_ait) + N_0)      ! [N_ait(t) + N_ait(t + dt)]/2

   q_number(mode_ait) = number_conc_new(mode_ait) / air_conc



   ! Compute mass mixing ratios changes due to coagulation between 
   ! source and receiving modes
   !  ------------------------------------------------------------

   ! maximum fraction of transfered volume = 1 - eps
   frac_transfer_vol_max = 1.0 - 1.0e1*epsilon(1.0_r8)

   ! first order loss rate from aitken to accumulation and primary carbon modes
   vol_loss = (beta_ij3(ait_acc) * number_conc_avg(mode_acc) + &
               beta_ij3(ait_pcm) * number_conc_avg(mode_pcm))

   ! fraction of 'i' volume transferred to 'j' mode
   frac_transfer_vol = 1.0 - exp(-vol_loss*dt)

   frac_transfer_vol = min(frac_transfer_vol_max, frac_transfer_vol)
   frac_transfer_vol = max(0.0, frac_transfer_vol)

   vol_shell = 0.0
   tmp_A = beta_ij3(ait_pcm)*number_conc_avg(mode_pcm)/max(vol_loss, 1.0e-37_r8)

   coag_pair = ait_acc
   mode_i = coag_mode_source(coag_pair)
   mode_j = coag_mode_receiv(coag_pair)

   do iq = 1, n_species(mode_i)
       iq_mode_i = iq
       iq_mode_j = intermodal_transfer(iq, coag_pair)

       if (iq_mode_j > 0) then
           ! species mass transfered from 'i' to 'j' mode
           mass_transfer = q_mass(iq_mode_i, mode_i)*frac_transfer_vol

           q_mass(iq_mode_i, mode_i) = q_mass(iq_mode_i, mode_i) - mass_transfer
           q_mass(iq_mode_j, mode_j) = q_mass(iq_mode_j, mode_j) + mass_transfer

           ! volume of shell material: SO4 and NH4 transfered from aitken to PCM
           vol_shell = vol_shell + (mass_transfer * tmp_A * mass2vol_aitken_age(iq))
       end if
   end do


   ! calculate aging transfer fraction for primary carbon to accumulation: 
   ! duplicates the code in CAM/MAM modal_aero_gasaerexch module

   ! use 1 mol (bi-)sulfate = 65 cm^3 --> 1 molecule = (4.76e-10 m)^3
   dR_so4_monolayers_pcage = NUMBER_SO4_MONOLAYERS_PCAGE * 4.76e-10

   ! volume to surface factor
   f_vol2sfc_pcarbon = exp(2.5 * (log(sigma(mode_pcm)))**2) 

   vol_core = 0.0
   do iq = 1, n_species(mode_pcm)
       vol_core = vol_core + q_mass(iq, mode_pcm) * mass2vol_pcarbon(iq)
   end do

   tmp_A = vol_shell * Dg_dry(mode_pcm) * f_vol2sfc_pcarbon
   tmp_B = 6.0 * dR_so4_monolayers_pcage * vol_core
   tmp_B = max(tmp_B, 0.0)

   if (tmp_A >= tmp_B) then
       frac_transfer_vol_pcage = frac_transfer_vol_max
   else
       frac_transfer_vol_pcage = min(tmp_A/tmp_B, frac_transfer_vol_max)
   end if


   ! calculate mass changes from primary carbon to accumulation by 
   ! direct coagulation and aging
   vol_loss = beta_ij3(pcm_acc)*number_conc_avg(mode_acc)

   frac_transfer_vol = 1.0 - exp(-vol_loss*dt)
   frac_transfer_vol = frac_transfer_vol + frac_transfer_vol_pcage

   frac_transfer_vol = min(frac_transfer_vol_max, frac_transfer_vol)
   frac_transfer_vol = max(0.0, frac_transfer_vol)

   coag_pair = pcm_acc
   mode_i = coag_mode_source(coag_pair)
   mode_j = coag_mode_receiv(coag_pair)
   
   do iq = 1, n_species(mode_i)                                   ! mass mixing ratios
       iq_mode_i = iq
       iq_mode_j = intermodal_transfer(iq, coag_pair)

       if (iq_mode_j > 0) then
           ! species mass transfered from 'i' to 'j' mode
           mass_transfer = q_mass(iq_mode_i, mode_i)*frac_transfer_vol

           q_mass(iq_mode_i, mode_i) = q_mass(iq_mode_i, mode_i) - mass_transfer
           q_mass(iq_mode_j, mode_j) = q_mass(iq_mode_j, mode_j) + mass_transfer
       end if
   end do

   number_transfer  = q_number(mode_i) * frac_transfer_vol_pcage  ! number mixing ratios
   q_number(mode_i) = q_number(mode_i) - number_transfer
   q_number(mode_j) = q_number(mode_j) + number_transfer

   return
   end subroutine MAML_Coagulation_AIT_PCM_ACC


   end module MAML_CoagulationMod
