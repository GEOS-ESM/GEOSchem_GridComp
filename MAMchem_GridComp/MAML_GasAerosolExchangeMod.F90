#include "MAPL_Exceptions.h"
#include "MAPL_Generic.h"

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: MAML_GasAerosolExchangeMod - Gas condensation.
!
! !INTERFACE:
!
   module MAML_GasAerosolExchangeMod
!
! !USES:
!
   use shr_kind_mod, only : r8 => shr_kind_r8

   use MAPL_Mod
   use MAPL_ConstantsMod, only : MAPL_PI

   use MAM_ComponentsDataMod
   use modal_aero_gasaerexch, only : modal_aero_soaexch


   implicit NONE
   private

!
! !PUBLIC MEMBER FUNCTIONS:

   public MAML_GasAerosolExchange


! !PRIVATE PARAMETERS 
   real, private, parameter :: pi = MAPL_PI
   real, private, parameter :: mw_air = MAPL_AIRMW  ! molecular weight of dry air, kg/Kmole


!
! !DESCRIPTION: 
!
!  {\tt MAML\_GasAerosolExchangeMod} provides methods to compute
!  gas-aerosol exchange.
!
!
! !REVISION HISTORY:
!
!  25Jun2012  A. Darmenov  Initial version -- based on CESM-1.0.3 CAM/MAM
!                                             modal_aero_gasaerexch module
!                                            
!
!EOP
!-------------------------------------------------------------------------




   contains

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: MAML_GasAerosolExchange --- condensation of H2SO4, NH3 and MSA.
!
!
! !INTERFACE:

   subroutine MAML_GasAerosolExchange(pressure,          &
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
                                      do_nh4g,           & ! << added
                                      do_msag,           & ! << added 
                                      do_soag,           & ! << added
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

   real, intent(inout) :: g_so4                ! gas H2SO4
   real, intent(inout) :: g_nh4                ! gas NH3
   real, intent(inout) :: g_msa                ! gas MSA
   real, intent(inout) :: g_soa                ! gas SOA(SOAG)



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


! !DESCRIPTION: Condensation of H2SO4, NH3, and MSA. Both treated as completely 
!               non-volatile (gas --> aerosol, but no aerosol --> gas): 
!               - gas H2SO4 goes to aerosol SO4
!               - gas MSA (if present) goes to aerosol SO4
!                 aerosol MSA is not distinguished from aerosol SO4
!               - gas NH3 (if present) goes to aerosol NH4
!                 if gas NH3 is not present, then...?

!
! !REVISION HISTORY:
!
!  25Jun2012  A. Darmenov  First crack -- based on modal_aero_gasaerexch_sub(),
!                          from CESM-1.0.3
!
!EOP
!-------------------------------------------------------------------------
                   __Iam__('MAML_GasAerosolExchange')


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



   ! mass to volume factors
   f_m2v_so4 = mw_so4_a / density_so4_a
   f_m2v_nh4 = mw_nh4_a / density_nh4_a
   f_m2v_soa = mw_soa_a / density_soa_a

   f_m2v_pcarbon(:) = 0.0
   n = modeptr_pcarbon
   do l = 1, nspec_amode(n)
       l2 = lspectype_amode(l,n)
       ! fac_m2v converts (kmol-AP/kmol-air) to (m3-AP/kmol-air)
       ! [m3-AP/kmol-AP]  = [kg-AP/kmol-AP]  / [kg-AP/m3-AP]
       fac_m2v_pcarbon(l) = mw_amode(l2) / density_amode(l2)
   end do

   ! volume to surface
   f_v2s_pcarbon = exp(2.5*(alnsg_amode(n)**2))
   xferfrac_max = 1.0 - 10.0*epsilon(1.0)   ! 1 - eps


   ! compute gas-to-aerosol mass transfer rates
   call gas_aer_uptkrates(q, t, pmid, dgncur_awet, uptkrate)


   ! use this for tendency calculations to avoid generating very small negative values
   dt_ = dt * (1.0 + epsilon(1.0))

#if(0)
   jsrf = jsrflx_gaexch

   ! f_gain_so4(n) = fraction of total H2SO4 uptake going to mode n
   ! f_gain_nh4(n) = fraction of total NH3   uptake going to mode n
   sum_uptk_rate_so4 = 0.0
   sum_uprt_nh4 = 0.0
   sum_uprt_soa = 0.0

   do n = 1, ntot_amode
       uptkratebb(n) = uptkrate(n,i,k)

       if (ido_so4a(n) > 0) then
           fgain_so4(n) = uptkratebb(n)
           sum_uprt_so4 = sum_uprt_so4 + fgain_so4(n)

           if (ido_so4a(n) == 1) then
               qold_so4(n) = q(i,k,lptr_so4_a_amode(n)-loffset)
           else
               qold_so4(n) = 0.0
           end if
       else
           fgain_so4(n) = 0.0
           qold_so4(n) = 0.0
       end if

       if (ido_nh4a(n) > 0) then
           ! 2.08 factor is for gas diffusivity (nh3/h2so4)
           ! differences in fuch-sutugin and accom coef ignored

           fgain_nh4(n) = uptkratebb(n)*2.08
           sum_uprt_nh4 = sum_uprt_nh4 + fgain_nh4(n)

           if (ido_nh4a(n) == 1) then
               qold_nh4(n) = q(i,k,lptr_nh4_a_amode(n)-loffset)
           else
               qold_nh4(n) = 0.0
           end if
       else
           fgain_nh4(n) = 0.0
           qold_nh4(n) = 0.0
       end if

       if (ido_soaa(n) > 0) then
           ! 0.81 factor is for gas diffusivity (soa/h2so4)
           ! (differences in fuch-sutugin and accom coef ignored)

           fgain_soa(n) = uptkratebb(n)*0.81
           sum_uprt_soa = sum_uprt_soa + fgain_soa(n)

           if (ido_soaa(n) == 1) then
               qold_soa(n) = q(i,k,lptr_soa_a_amode(n)-loffset)
               l = lptr_pom_a_amode(n)-loffset

               if (l > 0) then
                   qold_poa(n) = q(i,k,l)
               else
                   qold_poa(n) = 0.0
               end if
           else
               qold_soa(n) = 0.0
               qold_poa(n) = 0.0
           end if
       else
           fgain_soa(n) = 0.0
           qold_soa(n) = 0.0
           qold_poa(n) = 0.0
       end if

       uptkrate_soa(n) = fgain_soa(n)
   end do



   if (sum_uprt_so4 > 0.0) then
       do n = 1, ntot_amode
           fgain_so4(n) = fgain_so4(n) / sum_uprt_so4
       end do
   end if

   ! at this point (sum_uprt_so4 <= 0.0) only when all the fgain_so4 are zero
   if (sum_uprt_nh4 > 0.0) then
       do n = 1, ntot_amode
           fgain_nh4(n) = fgain_nh4(n) / sum_uprt_nh4
       end do
   end if

   if (sum_uprt_soa > 0.0) then
       do n = 1, ntot_amode
           fgain_soa(n) = fgain_soa(n) / sum_uprt_soa
       end do
   end if

   ! uptake amount (fraction of gas uptaken) over deltat
   avg_uprt_so4 = (1.0 - exp(-deltatxx*sum_uprt_so4))/deltatxx
   avg_uprt_nh4 = (1.0 - exp(-deltatxx*sum_uprt_nh4))/deltatxx
   avg_uprt_soa = (1.0 - exp(-deltatxx*sum_uprt_soa))/deltatxx
 
   ! sum_dqdt_so4 = so4_a tendency from h2so4 gas uptake (mol/mol/s)
   ! sum_dqdt_msa = msa_a tendency from msa   gas uptake (mol/mol/s)
   ! sum_dqdt_nh4 = nh4_a tendency from nh3   gas uptake (mol/mol/s)
   ! sum_dqdt_soa = soa_a tendency from soa   gas uptake (mol/mol/s)
   sum_dqdt_so4 = q(i,k,l_so4g) * avg_uprt_so4

   if (do_msag) then
       sum_dqdt_msa = q(i,k,l_msag) * avg_uprt_so4
   else
       sum_dqdt_msa = 0.0
   end if

   if (do_nh4g) then
       sum_dqdt_nh4 = q(i,k,l_nh4g) * avg_uprt_nh4
   else
       sum_dqdt_nh4 = 0.0
   end if

   if (do_soag) then
       sum_dqdt_soa = q(i,k,l_soag) * avg_uprt_soa
   else
       sum_dqdt_soa = 0.0
   end if


   ! compute TMR tendencies for so4, nh4, msa interstial aerosol
   ! due to simple gas uptake
   pdel_fac = pdel(i,k)/gravit
   sum_dqdt_nh4_b = 0.0

   do n = 1, ntot_amode
       dqdt_so4(n) = fgain_so4(n)*(sum_dqdt_so4 + sum_dqdt_msa)
 
       if (do_nh4g) then
           dqdt_nh4(n) = fgain_nh4(n)*sum_dqdt_nh4
           qnew_nh4 = qold_nh4(n) + dqdt_nh4(n)*deltat
           qnew_so4 = qold_so4(n) + dqdt_so4(n)*deltat
           qmax_nh4 = 2.0*qnew_so4

           if (qnew_nh4 > qmax_nh4) then
               dqdt_nh4(n) = (qmax_nh4 - qold_nh4(n))/deltatxx
           end if

           sum_dqdt_nh4_b = sum_dqdt_nh4_b + dqdt_nh4(n)
       end if
   end do

   if (( do_soag ) .and. (method_soa > 1)) then
       ! compute TMR tendencies for soag and soa interstial aerosol
       ! using soa parameterization

       niter_max   = 1000
       dqdt_soa(:) = 0.0

       call modal_aero_soaexch( deltat, t(i,k), pmid(i,k),      &
                                niter, niter_max, ntot_soamode, &
                                q(i,k,l_soag), qold_soa, qold_poa, uptkrate_soa, &
                                tmp1, dqdt_soa )

       sum_dqdt_soa = -tmp1
   else if ( do_soag ) then
       ! compute TMR tendencies for soa interstial aerosol
       ! due to simple gas uptake

       do n = 1, ntot_amode
           dqdt_soa(n) = fgain_soa(n)*sum_dqdt_soa
       end do
   else
       dqdt_soa(:) = 0.0
   end if
 
   do n = 1, ntot_amode
       if (ido_so4a(n) == 1) then
           l = lptr_so4_a_amode(n)-loffset
           dqdt(i,k,l) = dqdt_so4(n)
           qsrflx(i,l,jsrf) = qsrflx(i,l,jsrf) + dqdt_so4(n)*pdel_fac
       end if
 
       if (do_nh4g) then
           if (ido_nh4a(n) == 1) then
               l = lptr_nh4_a_amode(n)-loffset
               dqdt(i,k,l) = dqdt_nh4(n)
               qsrflx(i,l,jsrf) = qsrflx(i,l,jsrf) + dqdt_nh4(n)*pdel_fac
           end if
       end if
 
       if (do_soag) then
           if (ido_soaa(n) == 1) then
               l = lptr_soa_a_amode(n)-loffset
               dqdt(i,k,l) = dqdt_soa(n)
               qsrflx(i,l,jsrf) = qsrflx(i,l,jsrf) + dqdt_soa(n)*pdel_fac
           end if
       end if
   end do
 
   ! compute TMR tendencies for h2so4, nh3, and msa gas
   ! due to simple gas uptake
   l = l_so4g
   dqdt(i,k,l) = -sum_dqdt_so4
   qsrflx(i,l,jsrf) = qsrflx(i,l,jsrf) + dqdt(i,k,l)*pdel_fac

   if (do_msag) then
       l = l_msag
       dqdt(i,k,l) = -sum_dqdt_msa
       qsrflx(i,l,jsrf) = qsrflx(i,l,jsrf) + dqdt(i,k,l)*pdel_fac
   end if

   if (do_nh4g) then
       l = l_nh4g
       dqdt(i,k,l) = -sum_dqdt_nh4_b
       qsrflx(i,l,jsrf) = qsrflx(i,l,jsrf) + dqdt(i,k,l)*pdel_fac
   end if
 
   if (do_soag) then
       l = l_soag
       dqdt(i,k,l) = -sum_dqdt_soa
       qsrflx(i,l,jsrf) = qsrflx(i,l,jsrf) + dqdt(i,k,l)*pdel_fac
   end if
 
   ! compute TMR tendencies associated with primary carbon aging
   if (modefrm_pcage > 0) then
       n = modeptr_pcarbon
       vol_shell = deltat * (dqdt_so4(n)*fac_m2v_so4 + dqdt_nh4(n)*fac_m2v_nh4 + &
                             dqdt_soa(n)*fac_m2v_soa*soa_equivso4_factor )

       vol_core = 0.0
       do l = 1, nspec_amode(n)
           vol_core = vol_core + q(i,k,lmassptr_amode(l,n)-loffset)*fac_m2v_pcarbon(l)
       end do

       !   ratio1 = vol_shell/vol_core = 
       !      actual hygroscopic-shell-volume/carbon-core-volume after gas uptake
       !   ratio2 = 6.0_r8*dr_so4_monolayers_pcage/(dgncur_a*fac_volsfc_pcarbon)
       !      = (shell-volume corresponding to n_so4_monolayers_pcage)/core-volume 
       !      The 6.0/(dgncur_a*fac_volsfc_pcarbon) = (mode-surface-area/mode-volume)
       !   Note that vol_shell includes both so4+nh4 AND soa as "equivalent so4",
       !      The soa_equivso4_factor accounts for the lower hygroscopicity of soa.
       !
       !   Define xferfrac_pcage = min( 1.0, ratio1/ratio2)
       !   But ratio1/ratio2 == tmp1/tmp2, and coding below avoids possible overflow 
       !

       tmp1 = vol_shell*dgncur_a(i,k,n)*fac_volsfc_pcarbon
       tmp2 = max( 6.0_r8*dr_so4_monolayers_pcage*vol_core, 0.0_r8 )

       if (tmp1 >= tmp2) then
           xferfrac_pcage = xferfrac_max
       else
           xferfrac_pcage = min( tmp1/tmp2, xferfrac_max )
       end if

       if (xferfrac_pcage > 0.0_r8) then
           do iq = 1, nspecfrm_pcage
               lsfrm = lspecfrm_pcage(iq)-loffset
               lstoo = lspectoo_pcage(iq)-loffset
               xferrate = (xferfrac_pcage/deltat)*q(i,k,lsfrm)
               dqdt(i,k,lsfrm) = dqdt(i,k,lsfrm) - xferrate
               qsrflx(i,lsfrm,jsrf) = qsrflx(i,lsfrm,jsrf) - xferrate*pdel_fac

               if ((lstoo > 0) .and. (lstoo <= pcnst)) then
                   dqdt(i,k,lstoo) = dqdt(i,k,lstoo) + xferrate
                   qsrflx(i,lstoo,jsrf) = qsrflx(i,lstoo,jsrf) + xferrate*pdel_fac
               end if
           end do

           if (ido_so4a(modetoo_pcage) > 0) then
               l = lptr_so4_a_amode(modetoo_pcage)-loffset
               dqdt(i,k,l) = dqdt(i,k,l) + dqdt_so4(modefrm_pcage)
               qsrflx(i,l,jsrf) = qsrflx(i,l,jsrf) + dqdt_so4(modefrm_pcage)*pdel_fac
           end if

           if (ido_nh4a(modetoo_pcage) > 0) then
               l = lptr_nh4_a_amode(modetoo_pcage)-loffset
               dqdt(i,k,l) = dqdt(i,k,l) + dqdt_nh4(modefrm_pcage)
               qsrflx(i,l,jsrf) = qsrflx(i,l,jsrf) + dqdt_nh4(modefrm_pcage)*pdel_fac
           end if

           if (ido_soaa(modetoo_pcage) > 0) then
               l = lptr_soa_a_amode(modetoo_pcage)-loffset
               dqdt(i,k,l) = dqdt(i,k,l) + dqdt_soa(modefrm_pcage)
               qsrflx(i,l,jsrf) = qsrflx(i,l,jsrf) + dqdt_soa(modefrm_pcage)*pdel_fac
           end if

           end if

       end if


! set "temporary testing arrays"
   qold(:,:,:) = q(:,:,:)
   qqcwold(:,:,:) = qqcw(:,:,:)
   dqdtsv1(:,:,:) = dqdt(:,:,:)
   dqqcwdtsv1(:,:,:) = dqqcwdt(:,:,:)


!
! do renaming calcs
!
   dotendrn(:) = .false.
   dotendqqcwrn(:) = .false.
   dorename_atik(1:ncol,:) = .true.
   is_dorename_atik = .true.
   if (ncol >= -13579) then
      call modal_aero_rename_sub(                              &
                       'modal_aero_gasaerexch_sub',            &
                       lchnk,             ncol,      nstep,    &
                       loffset,           deltat,              &
                       latndx,            lonndx,              &
                       pdel,                                   &
                       dotendrn,          q,                   &
                       dqdt,              dqdt_other,          &
                       dotendqqcwrn,      qqcw,                &
                       dqqcwdt,           dqqcwdt_other,       &
                       is_dorename_atik,  dorename_atik,       &
                       jsrflx_rename,     nsrflx,              &
                       qsrflx,            qqcwsrflx            )
   end if


!
!  apply the dqdt to update q (and same for qqcw)
!
   do l = 1, pcnstxx
      if ( dotend(l) .or. dotendrn(l) ) then
         do k = 1, pver
         do i = 1, ncol
            q(i,k,l) = q(i,k,l) + dqdt(i,k,l)*deltat
         end do
         end do
      end if
      if ( dotendqqcw(l) .or. dotendqqcwrn(l) ) then
         do k = 1, pver
         do i = 1, ncol
            qqcw(i,k,l) = qqcw(i,k,l) + dqqcwdt(i,k,l)*deltat
         end do
         end do
      end if
   end do


!   do history file column-tendency fields
   do l = 1, pcnstxx
       lb = l + loffset

       do jsrf = 1, 2

       do jac = 1, 2

           if (jac == 1) then
               if (jsrf == jsrflx_gaexch) then
                   if ( .not. dotend(l) ) cycle

                   fieldname = trim(cnst_name(lb)) // '_sfgaex1'
               else if (jsrf == jsrflx_rename) then
                   if ( .not. dotendrn(l) ) cycle

                   fieldname = trim(cnst_name(lb)) // '_sfgaex2'
               else
                   cycle
               end if

               do i = 1, ncol
                   qsrflx(i,l,jsrf) = qsrflx(i,l,jsrf)*(adv_mass(l)/mwdry)
               end do

               call outfld( fieldname, qsrflx(:,l,jsrf), pcols, lchnk )
           else
               if (jsrf == jsrflx_gaexch) then
                   cycle
               else if (jsrf == jsrflx_rename) then
                   if ( .not. dotendqqcwrn(l) ) cycle
                   fieldname = trim(cnst_name_cw(lb)) // '_sfgaex2'
               else 
                   cycle
               end if

               do i = 1, ncol
                   qqcwsrflx(i,l,jsrf) = qqcwsrflx(i,l,jsrf)*(adv_mass(l)/mwdry)
               end do

               call outfld( fieldname, qqcwsrflx(:,l,jsrf), pcols, lchnk )
           end if

           if (ldiag4 > 0) then
               if (icol_diag > 0) then
                   i = icol_diag

                   if (jac == 1) then
                       tmp1 = qsrflx(i,l,jsrf)
                   else
                       tmp1 = qqcwsrflx(i,l,jsrf)
                   end if

                   write(*,'(a,4i5,2x,a,1p,2e12.4)') &
                             'gasaerexch nstep,lat,lon,l,fieldname,qsrflx,adv_mass',   &
                             nstep, latndx(i), lonndx(i), l, fieldname, tmp1, adv_mass(l)
               end if
           end if

       end do ! jac = ...
       end do ! jsrf = ...
   end do ! l = ...

#endif
   return
 
   end subroutine MAML_GasAerosolExchange


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: MAML_GasAerosolUptake --- Compute uptake rate due to gas 
!                                      condensation.
!
!
! !INTERFACE:

   subroutine MAML_GasAerosolUptake(pressure,          &
                                    temperature,       &
                                    density_air,       &
                                    q_number,          &
                                    Dg_wet,            &
                                    sigma,             &
                                    uptake)


! !USES:

   implicit NONE

! !INPUT/OUTPUT PARAMETERS:
   real, dimension(:), intent(inout) :: q_number ! number mixing ratios of the aerosol modes


   real, dimension(:), intent(inout) :: uptake   ! gas-to-aerosol mass transfer rate


! !INPUT PARAMETERS:
   real, intent(in) :: pressure                  ! pressure at mid level, Pa
   real, intent(in) :: temperature               ! temperature at mid level, K
   real, intent(in) :: density_air               ! air density, kg m-3

   real, dimension(:), intent(in)    :: Dg_wet   ! wet geometric mean diameter of number size distribution
   real, dimension(:), intent(in)    :: sigma    ! geometric standard deviation


! !OUTPUT PARAMETERS:


! !DESCRIPTION: Computes the H2SO4 uptake rate (gas to aerosol phase) for aerosol
!               population with lognormal size distribution N=N(ln(Dp))
!                             /
!               uptake rate = | dx * dN/dx * gas_conden_rate(Dp(x)), where
!                             /
!               
!               Dp    = particle diameter, cm
!               x     = ln(Dp)
!               dN/dx = log-normal particle number density distribution
!                
!               gas_conden_rate(Dp) = 2 * pi * gas_diffusivity * Dp * F(Kn,ac)
!               F(Kn,ac) = Fuchs-Sutugin correction factor
!               Kn = Knudsen number
!               ac = accommodation coefficient = 
!                  = 'number of molecules entering liquid phase' / 
!                    'number of molecular collisions with the surface'
!
!               The uptake rate is computed numerically using Gauss-Hermite 
!               quadrature of order 2.
!

!
! !REVISION HISTORY:
!
!  25Jun2012  A. Darmenov  First crack -- based on gas_aer_uptkrates(),
!                          from CESM-1.0.3
!
!EOP
!-------------------------------------------------------------------------
                   __Iam__('MAML_GasAerosolUptake')


   ! local parameters
   real, parameter    :: beta    = 2.0
   real, parameter    :: sqrt_2  = sqrt(2.0)
   real, parameter    :: sqrt_pi = sqrt(pi)

   integer, parameter :: n_ghq = 2             ! Gauss-Hermite quadrature order, abscissae and weights
   real, dimension(n_ghq), parameter :: x_ghq = (-0.70710678, 0.70710678)
   real, dimension(n_ghq), parameter :: w_ghq = ( 0.88622693, 0.88622693)

!  integer, parameter :: n_ghq = 3
!  real, dimension(n_ghq), parameter :: x_ghq = (-1.22474487, 0.000000000, 1.22474487)
!  real, dimension(n_ghq), parameter :: w_ghq = ( 0.29540897, 1.181635901, 0.29540897)
!
!  integer, parameter :: n_ghq = 4
!  real, dimension(n_ghq), parameter :: x_ghq = (-1.65068012,-0.524647623, 0.524647623, 1.65068012)
!  real, dimension(n_ghq), parameter :: w_ghq = ( 0.08131283, 0.804914090, 0.804914090, 0.08131283)
!
!  integer, parameter :: n_ghq = 5
!  real, dimension(n_ghq), parameter :: x_ghq = (-2.02018287,-0.958572465, 0.000000000, 0.958572465, 2.02018287)
!  real, dimension(n_ghq), parameter :: w_ghq = ( 0.01995324, 0.393619323, 0.945308720, 0.393619323, 0.01995324)

   
   ! local variables
   real :: air_con                             ! dry air molar concentration, kmol-air/m3
   real :: num_a_con                           ! aerosol number molar concentration, kmol m-3

   real :: diffusivity_h2so4                   ! diffusivity of H2SO4(gas), m2/s
   real :: speed_h2so4                         ! mean molecular speed of H2SO4(gas), m/s

   real :: mean_free_path                      ! mean free path, m

   real, dimension(n_ghq) :: ln_dp, dp         ! temporary vars
   real, dimension(n_ghq) :: Kn                ! Knudsen number
   real, dimension(n_ghq) :: fuchs_sutugin     ! Fuchs-Sutugin term

   real :: sum_ghq                             ! Gauss-Hermite quadrature

   real :: C, ln_Dg, ln_sigma                  ! temporary vars

   integer :: n, n_modes
   

   n_modes = size(q_number)


   ! dry air concentration
   air_con = density_air / mw_air

   ! following expressions for H2SO4(gas) are from MOSAIC
   diffusivity_h2so4 = 0.557e-4 * (temperature**1.75) / pressure
   speed_h2so4       = 1.470e1 * sqrt(temperature )

   
   ! Fuchs-Sutugin definition of mean free path
   mean_free_path = 3.0 * diffusivity_h2so4 / speed_h2so4

   do n = 1, n_modes       
       ! concentration of aerosol particles
       num_a_con = q_number(n) * air_con

       ln_Dg    = log(Dg_wet(n))
       ln_sigma = log(sigma(n))

       ! compute function values at gauss-hermite quadrature points
       ln_dp = ln_Dg + beta*ln_sigma**2 + sqrt_2*ln_sigma*x_ghq
       dp = exp(ln_Dp)

       ! knudsen number
       Kn = 2 * mean_free_path/dp

       ! apply accommodation coefficient (ac) = 0.65, after Adams & Seinfeld (JGR, 2002)
       ! fuchs_sutugin(Kn,ac) = (0.75*ac*(1 + Kn)) / (Kn*(1 + Kn + 0.283*ac) + 0.75*ac)
       fuchs_sutugin = (0.4875 * (1 + Kn)) / (Kn*(1.184 + Kn) + 0.4875)

       ! gauss-hermite quadrature 
       sum_ghq = sum(w_ghq * dp * fuchs_sutugin/(dp**beta))

       C = 2*sqrt_pi * num_a_con * exp(beta*ln_Dg + 0.5*(beta*ln_sigma)**2)
       uptake(n) = C * diffusivity_h2so4 * sum_ghq
   end do
 
   return
 
   end subroutine MAML_GasAerosolUptake


   end module MAML_GasAerosolExchangeMod
