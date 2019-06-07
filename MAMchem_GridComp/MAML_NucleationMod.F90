#include "MAPL_Exceptions.h"
#include "MAPL_Generic.h"

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: MAML_NucleationMod - Nucleation of aerosol particles.
!
! !INTERFACE:
!
   module MAML_NucleationMod
!
! !USES:
!

   use MAPL_Mod
   use MAPL_ConstantsMod, only : MAPL_PI, r8 => MAPL_R8

   use MAM_ComponentsDataMod
   use modal_aero_newnuc, only : mer07_veh02_nuc_mosaic_1box


   implicit NONE
   private

!
! !PUBLIC MEMBER FUNCTIONS:

   public MAML_Nucleation


! !PRIVATE PARAMETERS 
   real, private, parameter :: pi = MAPL_PI



!
! !DESCRIPTION: 
!
!  {\tt MAML\_NucleationMod} provides a collection of methods to calculate
!  binary and ternary nucleation rates.
!
!
! !REVISION HISTORY:
!
!  26Jan2012  A. Darmenov  Initial version -- based on CESM-1.0.3 CAM/MAM
!                                             modal_aero_newnuc module
!                                            
!
!EOP
!-------------------------------------------------------------------------


   interface MAML_Nucleation
       module procedure MAML_NucleationHomogeneous
   end interface MAML_Nucleation


   contains


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: MAML_NucleationHomogeneous --- 
!
! !INTERFACE:

   subroutine MAML_NucleationHomogeneous(pressure,          &
                                         temperature,       &
                                         density_air,       &
                                         rh,                &
                                         f_cld,             &
                                         z,                 &
                                         pblz,              &
                                         q_number,          &
                                         q_nh4,             &
                                         q_so4,             &
                                         q_h2so4,           &
                                         q_nh3,             &
                                         do_nh3,            &
                                         Dg,                &
                                         Dg_min,            &
                                         Dg_max,            &
                                         density_so4,       &
                                         mw_so4a,           &
                                         mw_nh4a,           &
                                         dq_h2so4_gasprod,  &
                                         dq_h2so4_aeruptk,  &
                                         dt)


! !USES:

   implicit NONE

! !INPUT/OUTPUT PARAMETERS:
   real, intent(inout) :: q_number             ! number mixing ratios of the Aitken mode
   real, intent(inout) :: q_nh4                ! mass mixing ratio of ammonium (NH4) in the Aitken mode
   real, intent(inout) :: q_so4                ! mass mixing ratio of sulfate  (SO4) in the Aitken mode

   real, intent(inout) :: q_h2so4              ! mass mixing ratio of sulfuric acid (H2SO4) 
   real, intent(inout) :: q_nh3                ! mass mixing ratio of ammonia (NH3)


! !INPUT PARAMETERS:
   real, intent(in)    :: pressure             ! pressure at mid level, Pa
   real, intent(in)    :: temperature          ! temperature at mid level, K
   real, intent(in)    :: density_air          ! air density, kg m-3
   real, intent(in)    :: rh                   ! relative humidity
   real, intent(in)    :: f_cld                ! cloud fraction
   real, intent(in)    :: z                    ! mid-layer height above surface, m
   real, intent(in)    :: pblz                 ! PBL height, m

   real, intent(in)    :: Dg                   ! mean diameter of Aitken mode number size distribution
   real, intent(in)    :: Dg_min               ! low limit of the mean diameter
   real, intent(in)    :: Dg_max               ! upper limit of the mean diameter

   real, intent(in)    :: density_so4          ! SO4 bulk density
   real, intent(in)    :: mw_so4a              ! molecular weight of SO4 aerosol
   real, intent(in)    :: mw_nh4a              ! molecular weight of NH4 aerosol


   real, intent(in)    :: dq_h2so4_gasprod     ! H2SO4 gas-phase production change over dt (mol/mol)
   real, intent(in)    :: dq_h2so4_aeruptk     ! H2SO4 gas-phase loss to aerosol over dt (mol/mol)

   real, intent(in)    :: dt                   ! time step

   logical, intent(in) :: do_nh3               ! NH3 flag
   

! !OUTPUT PARAMETERS:


! !DESCRIPTION: Homogeneous nucleation.
!
! !REVISION HISTORY:
!
!  03Jan2011  A. Darmenov  First crack -- based on modal_aero_coag_sub(),
!                          from CESM-1.0.3
!
!EOP
!-------------------------------------------------------------------------
                   __Iam__('MAML_NucleationHomogeneous')


   ! local variables
   real(r8) :: P, T                                      ! pressure and temperature at midlevels
   real(r8) :: zm                                        ! mid-level height
   real(r8) :: pblh                                      ! PBL height

   real(r8) :: d_dry_min, d_dry_max                      ! dry-diameter limits
   real(r8) :: mass_1p, mass_1p_min, mass_1p_max         ! single particle mass and mass limits

   real(r8) :: f                                         ! lognormal size distribution factor

   real(r8) :: mw_so4a_host                              ! molecular weght of sulfate aerosol

   real     :: cld                                       ! cloud fraction in the interval [0, 1]
   real     :: rh_grid                                   ! RH (grid average)
   real(r8) :: rh_non_cld                                ! RH in the non cloudy area of the grid

   real(r8) :: dqdt                                      ! tendency


   real(r8) :: deltat                                    ! time step

   real(r8) :: dplom_mode(1), dphim_mode(1)              ! dry-diameter limits

   real(r8) :: q_h2so4_cur                               ! H2SO4 molar mixing ratio
   real(r8) :: q_h2so4_avg                               ! estimated H2SO4, mol/mol-air
   real(r8) :: q_nh3_cur                                 ! NH3, mol/mol-air

                                                         ! note: aerosol changes are > 0; gas changes are < 0
   real(r8) :: dq_numa                                   ! change to aerosol number mixing ratio,  #/mol-air
   real(r8) :: dq_so4a                                   ! change to aerosol SO4 mixing ratio,   mol/mol-air
   real(r8) :: dq_nh4a                                   ! change to aerosol NH4 mixing ratio,   mol/mol-air
   real(r8) :: dq_h2so4                                  ! change to gas H2SO4 mixing ratio,     mol/mol-air
   real(r8) :: dq_nh3                                    ! change to gas NH3 mixing ratio,       mol/mol-air

   real(r8) :: dens_nh4so4a                              ! dry-density of the new NH4-SO4 aerosol mass, kg m-3

   real(r8) :: dndt_ait, dmdt_ait                        ! number and mass nucleation rates

   real(r8) :: dso4dt_ait                                ! 
   real(r8) :: dnh4dt_ait                                ! 
   real(r8) :: dqdt_numait, dqdt_nh4ait, dqdt_so4ait     !
   real(r8) :: dqdt_h2so4, dqdt_nh3                      !

   real(r8) :: tmp_a, tmp_b, tmp_c, tmp_q2, tmp_q3,   &  ! temporary vars
               tmp_uptake_rate, tmp_frso4
               
   real(r8) :: dndt_aitsv1, dndt_aitsv2, dndt_aitsv3, &  ! temporary values of the nucleation rates
               dmdt_aitsv1, dmdt_aitsv2, dmdt_aitsv3

   integer  :: l_diag_veh02                              ! diagnostics flag, -1 / +1 corresponds to disable / enable

   integer  :: itmp                                      ! size bin of newly formed particles


   ! local parameters
   real(r8), parameter :: q_h2so4_cutoff = 4.0e-16_r8    ! minimal H2SO4 vapor molar mixing ratio for nucleation = 4.0e-16 mol/mol-air,
                                                         ! which corresponds to approximatlelly 1.0e4 molecules/cm3

   integer, parameter  :: nuc_method_flag = 11           ! 1  = merikanto et al (2007) ternary
                                                         ! 2  = vehkamaki et al (2002) binary
                                                         ! 11 = merikanto ternary + first-order boundary layer
                                                         ! 12 = merikanto ternary + second-order boundary layer


   dqdt = 0.0

   ! dry-diameter limits for 'grown' new particles
   d_dry_min = exp(0.67*log(Dg_min) + 0.33*log(Dg))
   d_dry_max = Dg_max



   ! mass_1p_[min|max] = mass (kg) of so4 & nh4 in a single particle of diameter min|max
   !             (assuming same dry density for so4 & nh4)
   ! mass1p_aitlo - dp = dplom_mode(1)
   ! mass1p_aithi - dp = dphim_mode(1)
   f = (pi/6.0) * density_so4
   mass_1p_min = f * (d_dry_min**3)
   mass_1p_max = f * (d_dry_max**3)

   !  mw_so4a_host is molecular weght of sulfate aerosol in host code:
   !  if NH3/NH4 are simulated mw_so4a_host is equal to 96,
   !  - something else when NH3/NH4 are not simulated
   mw_so4a_host = mw_so4a



   ! if completely cloudy all H2SO4 vapor should be cloud-borne
   if (f_cld >= 0.99) &
       return

   ! current H2SO4 mixing ratio (after aerosol uptake)
   q_h2so4_cur = q_h2so4

   ! skip if H2SO4 vapor mixing ratio is less than q_h2so4_cutoff
   if (q_h2so4_cur <= q_h2so4_cutoff) &
       return
    

   tmp_a = max(0.0, dq_h2so4_gasprod)
   tmp_q3 = q_h2so4_cur

   ! tmp_q2 = qh2so4 before aerosol uptake, note that both 
   ! tmp_q3 and tmp_q2 are greater or equal to 0
   tmp_q2 = tmp_q3 + max(0.0, -dq_h2so4_aeruptk)


   ! tmp_b = log(tmp_q2 / tmp_q3) BUT with some checks added
   ! tmp_uptake_rate = tmp_b/dt
   if (tmp_q2 <= tmp_q3) then
       tmp_b = 0.0
   else
       tmp_c = tmp_q2 * exp(-20.0)

       if (tmp_q3 <= tmp_c) then
           tmp_q3 = tmp_c
           tmp_b  = 20.0_r8
       else
           tmp_b = log(tmp_q2/tmp_q3)
       end if
   end if

   ! d[ln(qh2so4)]/dt (1/s) from uptake (condensation) to aerosol
   tmp_uptake_rate = tmp_b/dt

   ! q_h2so4_avg = estimated average q_h2so4 when production and loss are done simultaneously
   if (tmp_b <= 0.1_r8) then
       q_h2so4_avg = tmp_q3*(1.0_r8 + 0.5_r8*tmp_b) - 0.5_r8*tmp_a
   else
       tmp_c = tmp_a/tmp_b
       q_h2so4_avg = ((tmp_q3 - tmp_c)*((exp(tmp_b) - 1.0_r8) / tmp_b)) + tmp_c
   end if

   if (q_h2so4_avg <= q_h2so4_cutoff) &
       return

   if (do_nh3) then
       q_nh3_cur = max(0.0_r8, q_nh3)
   else
       q_nh3_cur = 0.0_r8
   end if


   ! grid average RH
   rh_grid = max(0.0, min(1.0, rh))

   ! non-cloudy area RH
   cld = max(0.0_r8, f_cld)
   rh_non_cld = (rh_grid - cld) / (1.0 - cld)
   rh_non_cld = max(0.0_r8, min(1.0_r8, rh_non_cld))

   ! limit RH to between 0.1% and 99%
   rh_non_cld = max(0.01, min(0.99, rh_non_cld))


   ! call ... routine to get nucleation rates
   l_diag_veh02 = -1    ! diagnostics flag


   ! double precission
   T = temperature
   P = pressure
   zm = z
   deltat = dt

   dplom_mode(1) = d_dry_min
   dphim_mode(1) = d_dry_max

   call mer07_veh02_nuc_mosaic_1box(nuc_method_flag,    &      ! nucleation method
                                    deltat,             &      ! time step, s
                                    T,                  &      ! temperature, K
                                    rh_non_cld,         &      ! relative humidity, as fraction
                                    P,                  &      ! air pressure, Pa 
                                    zm,                 &      ! 
                                    pblh,               &      !
                                    q_h2so4_cur,        &      ! gas h2so4 mixing ratios (mol/mol-air) -- current value (after gas chem and condensation)
                                    q_h2so4_avg,        &      ! -- // --                              -- estimated average value (for simultaneous source/sink calcs)
                                    q_nh3_cur,          &      ! gas nh3 mixing ratios (mol/mol-air)   -- current value 
                                    tmp_uptake_rate,    &
                                    mw_so4a_host,       &
                                    1,                  &      ! ?? nsize  // number of aerosol size bins
                                    1,                  &      ! ??
                                    dplom_mode,         &      ! dry diameter at lower bnd of bin (m)
                                    dphim_mode,         &      ! dry diameter at upper bnd of bin (m)
                                    itmp,               &      ! size bin into which new particles go
                                    dq_numa,            &      ! change to aerosol number mixing ratio (#/mol-air) 
                                    dq_so4a,            &      ! change to aerosol so4 mixing ratio (mol/mol-air) -- aerosol changes are > 0
                                    dq_nh4a,            &      ! change to aerosol nh4 mixing ratio (mol/mol-air)
                                    dq_h2so4,           &      ! change to gas h2so4 mixing ratio (mol/mol-air)   -- gas changes are < 0
                                    dq_nh3,             &      ! change to gas nh3 mixing ratio (mol/mol-air)
                                    dens_nh4so4a,       &      ! dry-density of the new nh4-so4 aerosol mass (kg/m3)
                                    l_diag_veh02)              ! diagnostics


   ! convert dq_numa units from #/mol-air to #/kmol-air
   dq_numa = dq_numa * 1.0e3_r8 

   ! number nucleation rate, #/kmol-air/s
   dndt_ait = dq_numa/deltat

   ! fraction of mass nuc going to SO4
   tmp_a = dq_so4a * mw_so4a
   tmp_b = tmp_a + dq_nh4a*mw_nh4a
   tmp_frso4 = max(tmp_a, 1.0e-35_r8)/max(tmp_b, 1.0e-35_r8)

   ! mass nuc rate (kg/kmol-air/s or g/mol...)
   dmdt_ait = max(0.0_r8, (tmp_b/deltat) )

   dndt_aitsv1 = dndt_ait
   dmdt_aitsv1 = dmdt_ait
   dndt_aitsv2 = 0.0
   dmdt_aitsv2 = 0.0
   dndt_aitsv3 = 0.0
   dmdt_aitsv3 = 0.0

   if (dndt_ait < 1.0e2) then
       ! ignore newnuc if number rate < 100 #/kmol-air/s ~= 0.3 #/mg-air/d
       dndt_ait = 0.0
       dmdt_ait = 0.0
   else
       dndt_aitsv2 = dndt_ait
       dmdt_aitsv2 = dmdt_ait


       ! mirage2 code checked for complete H2SO4 depletion here, 
       ! but this is now done in mer07_veh02_nuc_mosaic_1box
       mass_1p = dmdt_ait/dndt_ait
       dndt_aitsv3 = dndt_ait
       dmdt_aitsv3 = dmdt_ait

       ! apply particle size constraints
       if (mass_1p < mass_1p_min) then
           ! reduce dndt to increase new particle size
           dndt_ait = dmdt_ait/mass_1p_min
       else if (mass_1p > mass_1p_max) then
           ! reduce dmdt to decrease new particle size
           dmdt_ait = dndt_ait*mass_1p_max
       end if
   end if


   ! set tendencies

   ! dso4dt_ait and dnh4dt_ait are in units kmol/kmol-air/s
   dso4dt_ait = dmdt_ait*tmp_frso4/mw_so4a
   dnh4dt_ait = dmdt_ait*(1.0_r8 - tmp_frso4)/mw_nh4a

   dqdt_h2so4 = -dso4dt_ait*(1.0 - cld)
   q_h2so4    = q_h2so4 + dqdt_h2so4*deltat

   dqdt_so4ait = dso4dt_ait*(1.0 - cld)
   q_so4       = q_so4 + dqdt_so4ait*deltat

   dqdt_numait = dndt_ait*(1.0 - cld)
   q_number = q_number + dqdt_numait*deltat

   if (do_nh3 .and. (dnh4dt_ait > 0.0_r8)) then
       dqdt_nh3 = -dnh4dt_ait*(1.0 - cld)
       q_nh3    = q_nh3 + dqdt_nh3*deltat

       dqdt_nh4ait = dnh4dt_ait*(1.0 - cld)
       q_nh4       = q_nh4 + dqdt_nh4ait*deltat
   end if


   return
 
   end subroutine MAML_NucleationHomogeneous


   end module MAML_NucleationMod
