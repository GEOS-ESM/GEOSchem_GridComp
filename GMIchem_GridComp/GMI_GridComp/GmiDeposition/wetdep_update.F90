!=============================================================================
!
! $Id$
!
! CODE DEVELOPER
!   Dan Bergmann, LLNL
!   dbergmann@llnl.gov
!
! FILE
!   wetdep_update.F
!
! ROUTINES
!   Update_Wetdep
!   Calc_Wet_Loss_Rate
!   E_Ice
!   Update_Wetdep_Orig_LLNL
!
! HISTORY
!   * March 16, 2005 - Jules Kouatchou
!     In the function E_Ice:
!        - Added two parameters "temp_min" (minimum temperature value) and
!          "temp_max" (maximum temperature value, equal to 273 degK).
!        - Lowered the value of "temp_min" to 100 degK (initially set to 153 degK)
!          to prevent the code from crashing. Such a low value is acceptable
!          because low temperature occur in the lower mesoshpere that does
!          not affect much GMI simulations.
!
!=============================================================================


!-----------------------------------------------------------------------------
!
! ROUTINE
!   Update_Wetdep
!
! DESCRIPTION
!   This routine calculates both in cloud scavenging and below cloud washout,
!   based on a writeup by Daniel Jacob at Harvard.  This is loosely based on
!   Giorgi and Chameides JGR 91 14,367-14376, 1986.
!
! ARGUMENTS
!   pr_wet_depos       : are wet depositions be written to the const output
!                        file?
!   pr_scav            : Is scavenging (3D) written to the history file?
!   chem_opt     : chemistry option
!   ih2o2_num    : const array index for H2O2
!   ihno3_num    : const array index for HNO3
!   tdt          : model time step  (s)
!   mw           : array of species' molecular weights (g/mol)
!   con_precip   : convective precipitation (mm/day)
!   tot_precip   : total      precipitation (mm/day)
!   mcor         : area of grid box (m^2)
!   grid_height  : grid box height  (m)
!   mass         : total mass of the atmosphere within each grid box   (kg)
!   moistq       : moisture changes due to wet processes (g/kg/day)
!   rain         : rainfall across cell edges (mm/day)
!   kel          : temperature (degK)
!   press3c      : atmospheric pressure at the center of each grid box (mb)
!   press3e      : atmospheric pressure at the edge   of each grid box (mb)
!   const        : species concentration, known at zone centers (mixing ratio)
!   wet_depos    : wet deposition accumulated since last output (kg/m^2)
!-micro_aerosol
!   humidity     : specific humidity (g/kg)
!   REL_SCAV_EFF_new
!                : relative scavenging efficiency for aerosols (0-1)
!
!-----------------------------------------------------------------------------

      subroutine Update_Wetdep  &
     &  (pr_wet_depos, pr_scav,  &
     &   chem_opt, ih2o2_num, ihno3_num, tdt, mw, &
     &   con_precip, tot_precip, mcor, grid_height, mass, moistq,  &
     &   rain_cn, rain_ls, &
     &   kel, press3c, press3e, concentration, wet_depos, scav3d, &
#ifdef MICRO_AEROSOL
     &   humidity, REL_SCAV_EFF_new, &
     &   pr_diag, loc_proc, i1, i2, ju1, j2, k1, k2, &
     &   ilo, ihi, julo, jhi, ilong, ivert, num_species)
#else
     &   pr_diag, loc_proc, i1, i2, ju1, j2, k1, k2, &
     &   ilo, ihi, julo, jhi, ilong, ivert, num_species)
#endif

      use GmiArrayBundlePointer_mod, only : t_GmiArrayBundle
      use GmiSpcConcentrationMethod_mod, only : isFixedConcentration

      implicit none

#     include "gmi_phys_constants.h"
#     include "gmi_time_constants.h"

#ifdef MICRO_AEROSOL
#     include "gmi_micro_aerosol.h"
#elif GOCARTaerosol
#     include "gocart_aerosol.h"
#else
#     include "gmi_aerosol.h"
#endif

!c?   Tight coupling to setkin?
#     include "setkin_par.h"
#     include "setkin_depos.h"
#ifdef MICRO_AEROSOL
#     include "umaerosol.h"
#endif


!     ----------------------
!     Argument declarations.
!     ----------------------

      logical           , intent(in   ) :: pr_diag
      integer           , intent(in   ) :: loc_proc
      integer           , intent(in   ) :: i1, i2, ju1, j2, k1, k2, ilong, ivert
      integer           , intent(in   ) :: ilo, ihi, julo, jhi
      integer           , intent(in   ) :: num_species

      logical           , intent(in   ) :: pr_wet_depos, pr_scav
      integer           , intent(in   ) :: chem_opt
      integer           , intent(in   ) :: ih2o2_num
      integer           , intent(in   ) :: ihno3_num
      real*8            , intent(in   ) :: tdt
      real*8            , intent(in   ) :: mw(num_species)
      real*8            , intent(in   ) :: con_precip (i1:i2,   ju1:j2)
      real*8            , intent(in   ) :: tot_precip (i1:i2,   ju1:j2)
      real*8            , intent(in   ) :: mcor       (i1:i2,   ju1:j2)
      real*8            , intent(in   ) :: grid_height(i1:i2,   ju1:j2,   k1:k2)
      real*8            , intent(in   ) :: mass       (i1:i2,   ju1:j2,   k1:k2)
      real*8            , intent(in   ) :: moistq     (i1:i2,   ju1:j2,   k1:k2)
      real*8            , intent(in   ) :: rain_cn    (i1:i2,   ju1:j2,   k1:k2)
      real*8            , intent(in   ) :: rain_ls    (i1:i2,   ju1:j2,   k1:k2)
      real*8            , intent(in   ) :: kel        (ilo:ihi, julo:jhi, k1:k2)
      real*8            , intent(in   ) :: press3c    (ilo:ihi, julo:jhi, k1:k2)
      real*8            , intent(in   ) :: press3e    (ilo:ihi, julo:jhi, k1-1:k2)
      real*8            , intent(inout) :: wet_depos  (i1:i2,   ju1:j2,   1:num_species)
      real*8            , intent(inout) :: scav3d     (i1:i2,   ju1:j2,   k1:k2,1:num_species)
      type (t_GmiArrayBundle), intent(inout) :: concentration(num_species)
#ifdef MICRO_AEROSOL
      real*8            , intent(in   ) :: humidity   (i1:i2,   ju1:j2,   k1:k2)
      real*8            , intent(in   ) :: REL_SCAV_EFF_new(i1:i2, ju1:j2, k1:k2, num_species)
#endif


!     -----------------------
!     Parameter declarations.
!     -----------------------

      real*8, parameter :: CM2PM2    = CMPM * CMPM  ! cm^2/m^2
      real*8, parameter :: CMPMM     = 0.1d0        ! cm/mm
      real*8, parameter :: WATER_DEN = 1.0d0        ! gr/cm^3

!     ----------------------
!     Variable declarations.
!     ----------------------

!     ---------------------------------------------------------
!     rainout : logical array pointing to locations for rainout
!     washout : logical array pointing to locations for washout
!     ---------------------------------------------------------

      logical :: rainout(i1:i2, ju1:j2)
      logical :: washout(i1:i2, ju1:j2)

      logical, save :: first = .true.

      integer :: il, ij, ik, ic

      real*8  :: f_max
      real*8  :: k_con, k_min
      real*8  :: l_plus_w_con, l_plus_w_str
      real*8  :: min_dtotau_1
      real*8  :: tau_con

!     ---------------------------------------------------------------------
!     conv_frac       : fraction of precipitation from convective processes
!     exp_fac         : temporary storage location
!     frac            : temporary storage location
!     frac_area_con   : fraction of grid box and above with convective
!                       precipitation
!     frac_area_str   : fraction of grid box and above with stratiform
!                       precipitation
!     frac_i_liq      : equilibrium fractionation of gas into rainwater
!     frac_not_evap   : fraction of precipitaion NOT evaporated in grid box
!     frac_wash_con   : washout fraction from convective precipitation
!     frac_wash_str   : washout fraction from stratiform precipitation
!     henry_eff       : effective Henry's law constant (mole/atm)
!     kloss           : species and temperature dependent wet scavenging
!                       rate constant (s^-1)
!     kloss_con       : kloss in convective clouds
!     kloss_str       : kloss in stratiform clouds
!     mbot            : mass of scaveged tracer through bottom of box (kg)
!     precip_bottom   : precipitation through bottom of box (g/cm^2/s)
!     precip_top      : precipitation through top    of box (g/cm^2/s)
!     q               : rate of evaporation (negative means condensation)
!                       (g/cm^3/s)
!     --------------------------------------------------------------------

      real*8  :: conv_frac      (i1:i2, ju1:j2)
      real*8  :: exp_fac        (i1:i2, ju1:j2)
      real*8  :: frac           (i1:i2, ju1:j2)
      real*8  :: frac_area_con  (i1:i2, ju1:j2)
      real*8  :: frac_area_str  (i1:i2, ju1:j2)
      real*8  :: frac_i_liq     (i1:i2, ju1:j2)
      real*8  :: frac_not_evap  (i1:i2, ju1:j2)
      real*8  :: frac_wash_con  (i1:i2, ju1:j2)
      real*8  :: frac_wash_str  (i1:i2, ju1:j2)
      real*8  :: henry_eff      (i1:i2, ju1:j2)
      real*8  :: kloss          (i1:i2, ju1:j2)
      real*8  :: kloss_con      (i1:i2, ju1:j2)
      real*8  :: kloss_str      (i1:i2, ju1:j2)
      real*8  :: mbot           (i1:i2, ju1:j2)
      real*8  :: precip_bot     (i1:i2, ju1:j2)
      real*8  :: precip_top     (i1:i2, ju1:j2)
      real*8  :: q              (i1:i2, ju1:j2)

!     --------------------------------------------------------------
!     scav_tracer_con : amount of scavenged tracer due to convective
!                       processes (kg)
!     scav_tracer_str : amount of scavenged tracer due to stratiform
!                       processes (kg)
!     --------------------------------------------------------------

      real*8  :: scav_tracer_con(i1:i2, ju1:j2, 1:num_species)
      real*8  :: scav_tracer_str(i1:i2, ju1:j2, 1:num_species)

#ifdef MICRO_AEROSOL
!-micro_aerosol-------------begin--------------------------------------
!     --------------------
!     adding for umaerosol
!     --------------------

      real*8  :: relhume   (i1:i2, ju1:j2)
      real*8  :: h2osat    (i1:i2, ju1:j2)
      real*8  :: h2ogas    (i1:i2, ju1:j2)
      real*8  :: so4mfrac  (i1:i2, ju1:j2)
      real*8  :: so4dens   (i1:i2, ju1:j2)
      real*8  :: wetmas    (i1:i2, ju1:j2)
      real*8  :: so4aer    (i1:i2, ju1:j2, naer)
      real*8  :: so4radv   (i1:i2, ju1:j2, nso4)
      real*8  :: so4radg   (i1:i2, ju1:j2, nso4)

      real*8  :: WASH_COEFF       (i1:i2, ju1:j2, num_species)
      real*8  :: frac_wash_str_new(i1:i2, ju1:j2)
      real*8  :: frac_wash_con_new(i1:i2, ju1:j2)
!-micro_aerosol-------------end----------------------------------------
#endif

!     ----------------
!     Begin execution.
!     ----------------

      if (pr_diag) then
        Write (6,*) 'Update_Wetdep called by ', loc_proc
      end if

!     ----------------------------------------------
!     In GEOS-5, the following construct is invalid.
!     ----------------------------------------------
!     ==========
!     if (first) then
!     ==========
        first = .false.

        IF(pr_wet_depos) THEN
         wet_depos(i1:i2,ju1:j2,1:num_species) = 0.0d0
        END IF
        IF(pr_scav) THEN
         scav3d(i1:i2,ju1:j2,k1:k2,1:num_species) = 0.0d0
        END IF

!     end if


#ifdef MICRO_AEROSOL
      kloss     (:,:) = 0.0d0
#endif
      henry_eff (:,:) = 0.0d0

      kloss_con (:,:) = 0.0d0
      kloss_str (:,:) = 0.0d0

      precip_bot(:,:) = 0.0d0

      scav_tracer_con(:,:,:) = 0.0d0
      scav_tracer_str(:,:,:) = 0.0d0


!     -----------------------------------------------------------------
!     Find the wet scavenged fraction using a first order rainout
!     parameterization from Giorgi and Chameides
!     (JGR v91 pp14,367-14376 1986).  Use this for both the large-scale
!     and the convective precipitation.
!     -----------------------------------------------------------------

!cc!? For now just hardwire a few values for the radon/lead problem
!c    (chem_opt = 1).  Setkin will eventually provide this information
!c    for all species when doing a full chemistry calculation.

      if (chem_opt == 1) then

        aerosol(1) = 0
        aerosol(2) = 1

        hstar_wet(1)   = 9.3d-3
        hstar_wet(2)   = 0.0d0

        delh_298_over_r_wet(1) = 2600.0d0
        delh_298_over_r_wet(2) = 0.0d0

        retention_eff(1)   = 0.0d0
        retention_eff(2)   = 0.0d0

      else if (chem_opt == 6) then

        aerosol(1) = 15
        aerosol(2) = 15

        hstar_wet(1)   = 0.0d0
        hstar_wet(2)   = 0.0d0

        delh_298_over_r_wet(1) = 0.0d0
        delh_298_over_r_wet(2) = 0.0d0

        retention_eff(1)   = 0.0d0
        retention_eff(2)   = 0.0d0

      else if (chem_opt == 8) then

        hstar_wet(IFSO2) = 600.0d0
        hstar_wet(INSO2) = 600.0d0

      end if

!     --------------------------------------------------------
!     First do the first-order rainout for large-scale clouds.
!     --------------------------------------------------------

!     -----------------------------------------------------
!     Calculate the fractional cloud coverage
!     (negative moistq means condensation).
!     Calculate a value of q for each level (g H2O/cm^3/s).
!     -----------------------------------------------------

      f_max        = 0.3d0

      k_con        = 1.5d-3
      k_min        = 1.0d-4

      l_plus_w_con = 2.0d-6
      l_plus_w_str = 0.5d-6

      tau_con      = 1800.0d0
      min_dtotau_1 = Min ((tdt / tau_con), 1.0d0)

      frac_area_con(:,:) = 1.0d-20
      frac_area_str(:,:) = 1.0d-20


!     =======
      ikloop: do ik = k2, k1, -1
!     =======

!     -----------------------------------------------------
!     When GMAO GEOS-5 [and 4!], conv_frac from the above
!     single-layer calculation is not used.  Instead it is
!     derived here at layer ik from 3D convective and large-
!     scale inputs.  Note that [the old] rain_hk is now 
!     assumed to be included in rain_cn.
!     -----------------------------------------------------

        conv_frac(:,:) = rain_cn(:,:,ik)
 
        WHERE(conv_frac(:,:) + rain_ls(:,:,ik) .NE. 0.0d0)  &
     &        conv_frac(:,:) = conv_frac(:,:) / &
     &        		(conv_frac(:,:) + rain_ls(:,:,ik))

        WHERE(conv_frac(:,:) < 0.0d0) conv_frac(:,:) = 0.0d0

        WHERE(conv_frac(:,:) > 1.0d0) conv_frac(:,:) = 1.0d0


        q(:,:) =  &
     &    -moistq(:,:,ik) *  &
     &    (press3c(i1:i2,ju1:j2,ik) * MWTAIR * MB2CGS) /  &
     &    (GPKG * SECPDY * kel(i1:i2,ju1:j2,ik) * BOLTZMN_E * AVOGAD)

!       ---------------------------------------------------------------------
!       Store the precipitation from the bottom of the previous box into
!       the top of the current box.  Calculate new precipitation value for
!       the bottom of the current box.  These will be used later for washout.
!       ---------------------------------------------------------------------

        precip_top(:,:) = precip_bot(:,:)

        precip_bot(:,:) =  &
     &    Max (0.0d0,  &
     &         precip_bot(:,:) -  &
     &         (moistq(:,:,ik) * mass(:,:,ik) /  &
     &          (SECPDY * mcor(i1:i2,ju1:j2) * CM2PM2)))

!       ---------------------------------------------------------------------
!       Calculate fractional area of convective and stratiform precipitation.
!       This is equation no. 13 in Daniel Jacob's scavenging paper.
!       ---------------------------------------------------------------------

        where (q(:,:) > 0.0d0)

          frac_area_con(:,:) =  &
     &      Max (f_max*q(:,:) * conv_frac(:,:) * min_dtotau_1 /  &
     &           (q(:,:) * conv_frac(:,:) * min_dtotau_1 +  &
     &            f_max*k_con * l_plus_w_con),  &
     &           frac_area_con(:,:))

          frac_area_str(:,:) =  &
     &      Max (q(:,:) * (1.0d0 - conv_frac(:,:)) /  &
     &           ((k_min + q(:,:) * (1.0d0 - conv_frac(:,:)) /  &
     &             l_plus_w_str) *  &
     &            l_plus_w_str),  &
     &           frac_area_str(:,:))

        end where

!       ---------------------------------------------------------
!       Calculate some fractional washout terms which may be used
!       regardless of species.
!       ---------------------------------------------------------

        where (frac_area_str(:,:) > 0.0d0)

          frac_wash_str(:,:) =  &
     &      frac_area_str(:,:) *  &
     &      (1.0d0 - Exp (-precip_bot(:,:) * (1.0d0 - conv_frac(:,:)) *  &
     &       tdt / frac_area_str(:,:)))

        elsewhere

          frac_wash_str(:,:) = 0.0d0

        end where

        where (frac_area_con(:,:) > 0.0d0)

          frac_wash_con(:,:) =  &
     &      frac_area_con(:,:) *  &
     &      (1.0d0 - Exp (-precip_bot(:,:) * conv_frac(:,:) *  &
     &       tdt / frac_area_con(:,:)))

        elsewhere

          frac_wash_con(:,:) = 0.0d0

        end where

!       ------------------------------------------------------
!       Create an array showing where washout occurs.
!       It only occurs where there is no precipitation forming
!       and the temperature is above 268 K.
!       ------------------------------------------------------

        washout(:,:) = ((moistq(:,:,ik)        >=   0.0d0) .and.  &
     &                  (kel (i1:i2,ju1:j2,ik) >= 268.0d0) .and.  &
     &                  (precip_top(:,:)       >    0.0d0))

!       ----------------------------------------------------
!       Create an array showing where rainout occurs.
!       It only occurs where there is precipitation forming.
!       ----------------------------------------------------

        rainout(:,:) = (moistq(:,:,ik) < -1.0d-10)

!       ------------------------------------------------------------------
!       First do the reevaporation by calculating fraction NOT evaporated.
!       ------------------------------------------------------------------

        where (washout(:,:))

          frac_not_evap(:,:) =  &
     &      1.0d0 -  &
     &      ((precip_top(:,:) - precip_bot(:,:)) / precip_top(:,:)) *  &
     &      0.5d0

        elsewhere

          frac_not_evap(:,:) = 1.0d0

        end where

#ifdef MICRO_AEROSOL
!-micro_aerosol-------------begin-------------------------------------------
!c      --------------------
!c      adding for umaerosol
!c      --------------------

!c      compute relative humidity (0-1)

        relhume(:,:) = 1.0d0 - (373.15d0 / kel(i1:i2,ju1:j2,ik))
        relhume(:,:) =  &
     &    1013.25d0 * Exp (13.3185d0 * relhume(:,:)    -  &
     &                      1.9760d0 * relhume(:,:)**2 -  &
     &                      0.6445d0 * relhume(:,:)**3 -  &
     &                      0.1299d0 * relhume(:,:)**4)

        relhume(:,:) =  &
     &    humidity(:,:,ik) * MWTAIR / 18.0d0 /  &
     &    GPKG * press3c(i1:i2,ju1:j2,ik) / relhume(:,:)

        relhume(:,:) = Max (Min (relhume(:,:), 0.95d0), 0.0d0)

!c      compute so4 aerosol size (m)

        so4aer(:,:,1) = concentration(ISO4M1)%pArray3D(:,:,ik)
        so4aer(:,:,2) = concentration(ISO4N1)%pArray3D(:,:,ik)
        so4aer(:,:,3) = concentration(ISO4M2)%pArray3D(:,:,ik)
        so4aer(:,:,4) = concentration(ISO4N2)%pArray3D(:,:,ik)

        do il = i1, i2
          do ij = ju1, j2
            h2osat(il,ij) = h2osat_f(kel(il,ij,ik))
            h2ogas(il,ij) = relhume(il,ij)*h2osat(il,ij)

            so4mfrac(il,ij) = so4mfrac_f(kel(il,ij,ik),  &
     &                                   h2osat(il,ij),  &
     &                                   h2ogas(il,ij))
            so4dens(il,ij) = so4dens_f(kel(il,ij,ik),  &
     &                                 so4mfrac(il,ij))
          end do
        end do

        do iso4 = 1, nso4

          iso4n = iso4 * nmomso4
          iso4m = iso4n - 1

          wetmas(:,:) = max(r2so4min, so4aer(:,:,iso4m)  &
     &                / max(epsilo,   so4aer(:,:,iso4n)))  &
     &                / so4mfrac(:,:)

          so4radv(:,:,iso4) = (r3q * wetmas(:,:)  &
     &                      / (GMI_PI * so4dens(:,:)))**r1td
          so4radv(:,:,iso4) = max(so4radvmin,  &
     &                        min(so4radv(:,:,iso4),so4radvmax))

        end do

!c      calculate washout coefficients (cm^-1) of so4 aerosol

        WASH_COEFF(:,:,:) = 0.0d0
        rg(:) = rg(:) * 1.0d-6         ! um -> meter

        do il = i1, i2
          do ij = ju1, j2

            do iso4 = 1, nso4

              xlnsg = log(sigmod(iso4))
              xlnsg2 = xlnsg * xlnsg

              so4radg(il,ij,iso4) = so4radv(il,ij,iso4)  &
     &                            * exp(-r3h*xlnsg2)

              if (so4radg(il,ij,iso4) <= rg(1)) then

                if (iso4 == 1) then
                  WASH_COEFF(il,ij,ISO4M1) = k1_q3(1)
                  WASH_COEFF(il,ij,ISO4N1) = k0_q0(1)
                else if (iso4 == 2) then
                  WASH_COEFF(il,ij,ISO4M2) = k2_q3(1)
                  WASH_COEFF(il,ij,ISO4N2) = k0_q0(1)
                end if

              else if (so4radg(il,ij,iso4) >= rg(nrg)) then

                if (iso4 == 1) then
                  WASH_COEFF(il,ij,ISO4M1) = k1_q3(nrg)
                  WASH_COEFF(il,ij,ISO4N1) = k0_q0(nrg)
                else if (iso4 == 2) then
                  WASH_COEFF(il,ij,ISO4M2) = k2_q3(nrg)
                  WASH_COEFF(il,ij,ISO4N2) = k0_q0(nrg)
                end if

              else

                irg = 1
                do while((rg(irg) < so4radg(il,ij,iso4))  &
     &                 .and. (irg < nrg))
                  irg = irg + 1
                end do

                if (iso4 == 1) then
                  WASH_COEFF(il,ij,ISO4M1) =  &
     &                             k1_q3(irg-1)  &
     &                           + (k1_q3(irg)-k1_q3(irg-1))  &
     &                           / (rg(irg)-rg(irg-1))  &
     &                           * (so4radg(il,ij,iso4)-rg(irg-1))
                  WASH_COEFF(il,ij,ISO4N1) =  &
     &                             k0_q0(irg-1)  &
     &                           + (k0_q0(irg)-k0_q0(irg-1))  &
     &                           / (rg(irg)-rg(irg-1))  &
     &                           * (so4radg(il,ij,iso4)-rg(irg-1))
                else if (iso4 == 2) then
                  WASH_COEFF(il,ij,ISO4M2) =  &
     &                             k2_q3(irg-1)  &
     &                           + (k2_q3(irg)-k2_q3(irg-1))  &
     &                           / (rg(irg)-rg(irg-1))  &
     &                           * (so4radg(il,ij,iso4)-rg(irg-1))
                  WASH_COEFF(il,ij,ISO4N2) =  &
     &                             k0_q0(irg-1)  &
     &                           + (k0_q0(irg)-k0_q0(irg-1))  &
     &                           / (rg(irg)-rg(irg-1))  &
     &                           * (so4radg(il,ij,iso4)-rg(irg-1))
                end if

              end if

            end do
          end do
        end do

        WASH_COEFF(:,:,:) = 10.0d0 * WASH_COEFF(:,:,:) ! mm-1 -> cm-1

!-micro_aerosol-------------end---------------------------------------------
#endif

!       =======
        icloop: do ic = 1, num_species
        kloss_str(:,:) = 0.0d0
        kloss_con(:,:) = 0.0d0

!       =======

!         ====================================
          if (isFixedConcentration(ic)) cycle icloop
!         ====================================

#ifdef MICRO_AEROSOL
          if (ic == ihno3_num .or. ic == ISO4G) then
#else
          if (ic == ihno3_num) then
#endif

            where (q(:,:) > 0.0d0)
              kloss_con(:,:) = k_con
              kloss_str(:,:) = k_min +  &
     &          q(:,:) * (1.0d0 - conv_frac(:,:)) / l_plus_w_str
            elsewhere
              kloss_con(:,:) = 0.0d0
              kloss_str(:,:) = 0.0d0
            end where

          else if (aerosol(ic) >= 1) then

            where (q(:,:) > 0.0d0)
#ifdef MICRO_AEROSOL
              kloss_con(:,:) = k_con * REL_SCAV_EFF_new(:,:,ik,ic)
#else
              kloss_con(:,:) = k_con * REL_SCAV_EFF(aerosol(ic))
#endif
              kloss_str(:,:) = (k_min +  &
     &          q(:,:) * (1.0d0 - conv_frac(:,:)) / l_plus_w_str) &
#ifdef MICRO_AEROSOL
     &          * REL_SCAV_EFF_new(:,:,ik,ic)
#else
     &          * REL_SCAV_EFF(aerosol(ic))
#endif
            elsewhere
              kloss_con(:,:) = 0.0d0
              kloss_str(:,:) = 0.0d0
            end where

          else if (hstar_wet(ic) > 0.0d0) then

            do ij = ju1, j2

!-micor_aerosol : note gmi version of code is kept, the following is called
!                 twice, in first call kloss_con is used and in second call
!                 kloss_str is used
!             =======================
              call Calc_Wet_Loss_Rate  &
!             =======================
     &          (.true., ic, ih2o2_num, delh_298_over_r_wet(ic), hstar_wet(ic),  &
     &           retention_eff(ic), press3e(:,ij,ik),  &
     &           kel(:,ij,ik), kloss_con(:,ij),  &
     &           i1, i2, ilo, ihi)

!             =======================
              call Calc_Wet_Loss_Rate  &
!             =======================
     &          (.false., ic, ih2o2_num, delh_298_over_r_wet(ic), hstar_wet(ic),  &
     &           retention_eff(ic), press3e(:,ij,ik),  &
     &           kel(:,ij,ik), kloss_str(:,ij),  &
     &           i1, i2, ilo, ihi)

            end do

!-micro_aerosol: note: gmi version is kept, kloss_con and kloss_str is used
!                      in place kloss
            where (q(:,:) > 0.0d0)
              kloss_con(:,:) = k_con * kloss_con(:,:)
              kloss_str(:,:) = (k_min +  &
     &          q(:,:) * (1.0d0 - conv_frac(:,:)) / l_plus_w_str)  &
     &          * kloss_str(:,:)
            elsewhere
              kloss_con(:,:) = 0.0d0
              kloss_str(:,:) = 0.0d0
            end where

#ifdef MICRO_AEROSOL
!-micro_aerosol------------------begin-----------------------------------
          else

            kloss_con(:,:) = 0.0d0
            kloss_str(:,:) = 0.0d0

!-micro_aerosol------------------end-------------------------------------
#endif
          end if

#ifdef MICRO_AEROSOL
!-micro_aerosol------------------begin-----------------------------------
!         ----------------------------------------------------------
!         Calculate some fractional washout terms related to species
!         ----------------------------------------------------------

          if (aerosol(ic) == 1 .or. aerosol(ic) == 2) then

            where (frac_area_str(:,:) > 0.0d0)

              frac_wash_str_new(:,:) =  &
     &          frac_area_str(:,:) *  &
     &          (1.0d0 - Exp (-precip_bot(:,:) *  &
     &          (1.0d0 - conv_frac(:,:)) * WASH_COEFF(:,:,ic)  &
     &           * tdt / frac_area_str(:,:)))

            elsewhere

              frac_wash_str_new(:,:) = 0.0d0

            end where

            where (frac_area_con(:,:) > 0.0d0)

              frac_wash_con_new(:,:) =  &
     &          frac_area_con(:,:) *  &
     &          (1.0d0 - Exp (-precip_bot(:,:)  &
     &          * conv_frac(:,:) * WASH_COEFF(:,:,ic)  &
     &          * tdt / frac_area_con(:,:)))

            elsewhere

              frac_wash_con_new(:,:) = 0.0d0

            end where

          else

            frac_wash_str_new(:,:) = frac_wash_str(:,:)
            frac_wash_con_new(:,:) = frac_wash_con(:,:)

          end if

!-micro_aerosol------------------end-------------------------------------
#endif
!         -----------------------------------------------------------
!         Only apply rainout in grid boxes where precipitation forms.
!         -----------------------------------------------------------

          where (rainout(:,:))

!           --------------------------
!           Do the stratiform rainout.
!           --------------------------

            exp_fac(:,:) = Exp (-kloss_str(:,:) * tdt)

            scav_tracer_str(:,:,ic) =  &
     &        scav_tracer_str(:,:,ic) +  &
     &        frac_area_str(:,:) * (1.0d0 - exp_fac(:,:)) *  &
     &        concentration(ic)%pArray3D(:,:,ik) * mass(:,:,ik) *  &
     &        mw(ic) / MWTAIR

            concentration(ic)%pArray3D(:,:,ik) =  &
     &        concentration(ic)%pArray3D(:,:,ik) *  &
     &        (1.0d0 -  &
     &         (frac_area_str(:,:) * (1.0d0 - exp_fac(:,:))))

!           --------------------------
!           Do the convective rainout.
!           --------------------------

            exp_fac(:,:) = Exp (-kloss_con(:,:) * tdt)

            scav_tracer_con(:,:,ic) =  &
     &        scav_tracer_con(:,:,ic) +  &
     &        frac_area_con(:,:) * (1.0d0 - exp_fac(:,:)) *  &
     &        concentration(ic)%pArray3D(:,:,ik) * mass(:,:,ik) *  &
     &        mw(ic) / MWTAIR

            concentration(ic)%pArray3D(:,:,ik) =  &
     &        concentration(ic)%pArray3D(:,:,ik) *  &
     &        (1.0d0 -  &
     &         (frac_area_con(:,:) * (1.0d0 - exp_fac(:,:))))

          end where

!         -------------------------------------------------------------
!         Do the washout and reevaporation of aerosols and nitric acid.
!         -------------------------------------------------------------

#ifdef MICRO_AEROSOL
!-micro_aerosol------------------begin--------------------------------
!         if ((aerosol(ic) >= 1) .or. (ic == ihno3_num)) then
          if ((aerosol(ic) >= 1) .or. (ic == ihno3_num) .or.  &
     &        (ic == ISO4G)) then

            if ((ic == ISO4G) .or. (ic == ISO4M1) .or.  &
     &          (ic == ISO4M2)) then

!-micro_aerosol------------------end----------------------------------

#else
          if ((aerosol(ic) >= 1) .or. (ic == ihno3_num)) then
#endif

            where (washout(:,:))

#ifdef MICRO_AEROSOL
!-micro_aerosol------------------begin--------------------------------

                concentration(ISO4M2)%pArray3D(:,:,ik) =  &
     &            concentration(ISO4M2)%pArray3D(:,:,ik) + &
!-micro_aerosol------------------end----------------------------------
#else
              concentration(ic)%pArray3D(:,:,ik) =  &
     &          concentration(ic)%pArray3D(:,:,ik) + &
#endif
     &          (1.0d0 - frac_not_evap(:,:)) *  &
     &          (scav_tracer_str(:,:,ic) + scav_tracer_con(:,:,ic)) *  &
     &          (MWTAIR / mw(ic)) / mass(:,:,ik)

!mkch
#ifdef MICRO_AEROSOL
!-micro_aerosol------------------begin--------------------------------

              end where

            else if ((ic == ISO4N1) .or. (ic == ISO4N2)) then

              where (washout(:,:))
!               concentration(ISO4N2)%pArray3D(:,:,ik) =
!    &            concentration(ISO4N2)%pArray3D(:,:,ik) +
!    &            (1.0d0 - frac_not_evap(:,:)) *
!    &            (scav_tracer_str(:,:,ic-1) +
!    &             scav_tracer_con(:,:,ic-1)) *
!    &            (MWTAIR / mw(ic)) / mass(:,:,ik) *
!    &            concentration(ISO4N2)%pArray3D(:,:,ik) / concentration(ISO4M2)%pArray3D(:,:,ik)
                concentration(ISO4N1)%pArray3D(:,:,ik) =  &
     &            concentration(ISO4N1)%pArray3D(:,:,ik)
                concentration(ISO4N2)%pArray3D(:,:,ik) =  &
     &            concentration(ISO4N2)%pArray3D(:,:,ik)
              end where

            else

              where (washout(:,:))
                concentration(ic)%pArray3D(:,:,ik) =  &
     &            concentration(ic)%pArray3D(:,:,ik) +  &
     &            (1.0d0 - frac_not_evap(:,:)) *  &
     &            (scav_tracer_str(:,:,ic) + scav_tracer_con(:,:,ic)) *  &
     &            (MWTAIR / mw(ic)) / mass(:,:,ik)
              end where

            end if

            where (washout(:,:))
!-micro_aerosol------------------end----------------------------------
#endif

              scav_tracer_str(:,:,ic) =  &
     &          scav_tracer_str(:,:,ic) * frac_not_evap(:,:)

              scav_tracer_con(:,:,ic) =  &
     &          scav_tracer_con(:,:,ic) * frac_not_evap(:,:)

!             ------------------------------
!             Now do the stratiform washout.
!             ------------------------------

              scav_tracer_str(:,:,ic) =  &
     &          scav_tracer_str(:,:,ic) + &
#ifdef MICRO_AEROSOL
     &          frac_wash_str_new(:,:) * &
#else
     &          frac_wash_str(:,:) * &
#endif
     &          concentration(ic)%pArray3D(:,:,ik) * mass(:,:,ik) *  &
     &          mw(ic) / MWTAIR

              concentration(ic)%pArray3D(:,:,ik) = &
#ifdef MICRO_AEROSOL
     &          concentration(ic)%pArray3D(:,:,ik) * (1.0d0 - frac_wash_str_new(:,:))
#else
     &          concentration(ic)%pArray3D(:,:,ik) * (1.0d0 - frac_wash_str(:,:))
#endif

!             ------------------------------
!             Now do the convective washout.
!             ------------------------------

              scav_tracer_con(:,:,ic) =  &
     &          scav_tracer_con(:,:,ic) +  &
#ifdef MICRO_AEROSOL
     &          frac_wash_con_new(:,:) * &
#else
     &          frac_wash_con(:,:) * &
#endif
     &          concentration(ic)%pArray3D(:,:,ik) * mass(:,:,ik) *  &
     &          mw(ic) / MWTAIR

              concentration(ic)%pArray3D(:,:,ik) = &
#ifdef MICRO_AEROSOL
     &          concentration(ic)%pArray3D(:,:,ik) * (1.0d0 - frac_wash_con_new(:,:))
#else
     &          concentration(ic)%pArray3D(:,:,ik) * (1.0d0 - frac_wash_con(:,:))
#endif

            end where

          else

!           ---------------------------------------------
!           Do washout and reevaporation of species other
!           than aerosols or nitric acid.
!           ---------------------------------------------

            henry_eff(:,:) =  &
     &        hstar_wet(ic) *  &
     &        Exp (-delh_298_over_r_wet(ic) *  &
     &             ((1.0d0 / kel(i1:i2,ju1:j2,ik)) -  &
     &              (1.0d0 / 298.0d0)))

!           -----------------------------------------------
!           Stratiform or large scale component of washout.
!           -----------------------------------------------

            frac(:,:) =  &
     &        henry_eff(:,:) * GAS_CONST_J / 100.0d0 *  &
     &        kel(i1:i2,ju1:j2,ik) * precip_bot(:,:) *  &
     &        (1.0d0 - conv_frac(:,:)) * tdt /  &
     &        (frac_area_str(:,:) * (grid_height(:,:,ik) * 100.0d0))

            frac_i_liq(:,:) =  &
     &        frac(:,:) / (1.0d0 + frac(:,:))

!           ---------------------------------------
!           Scavenging is limited by mass transfer.
!           ---------------------------------------

            where ((frac_i_liq(:,:) >=  &
     &              (frac_wash_str(:,:) / frac_area_str(:,:))) .and.  &
     &             washout(:,:))

!             ---------------------------
!             First do the reevaporation.
!             ---------------------------

              concentration(ic)%pArray3D(:,:,ik) =  &
     &          concentration(ic)%pArray3D(:,:,ik) +  &
     &          (1.0d0 - frac_not_evap(:,:)) * scav_tracer_str(:,:,ic) *  &
     &          (MWTAIR / mw(ic)) /  mass(:,:,ik)

              scav_tracer_str(:,:,ic) =  &
     &          scav_tracer_str(:,:,ic) * frac_not_evap(:,:)

!             ------------------------------
!             Now do the stratiform washout.
!             ------------------------------

              scav_tracer_str(:,:,ic) =  &
     &          scav_tracer_str(:,:,ic) +  &
     &          frac_wash_str(:,:) *  &
     &          concentration(ic)%pArray3D(:,:,ik) * mass(:,:,ik) *  &
     &          mw(ic) / MWTAIR

              concentration(ic)%pArray3D(:,:,ik) =  &
     &          concentration(ic)%pArray3D(:,:,ik) * (1.0d0 - frac_wash_str(:,:))

            end where

!           -------------------------------------------------
!           Scavenging is limited by Henry's law equilibrium.
!           -------------------------------------------------

            where ((frac_i_liq(:,:) <  &
     &              (frac_wash_str(:,:) / frac_area_str(:,:))) .and.  &
     &             washout(:,:))

              mbot(:,:) =  &
     &          frac_i_liq(:,:) *  &
     &          (frac_area_str(:,:) * concentration(ic)%pArray3D(:,:,ik) *  &
     &           mass(:,:,ik) * mw(ic) / MWTAIR +  &
     &           scav_tracer_str(:,:,ic))

              concentration(ic)%pArray3D(:,:,ik) =  &
     &          concentration(ic)%pArray3D(:,:,ik) -  &
     &          (mbot(:,:) - scav_tracer_str(:,:,ic)) /  &
     &          mass(:,:,ik) * MWTAIR / mw(ic)

              scav_tracer_str(:,:,ic) = mbot(:,:)

            end where

!           --------------------------------
!           Convective component of washout.
!           --------------------------------

            frac(:,:) =  &
     &        henry_eff(:,:) * GAS_CONST_J / 100.0d0 *  &
     &        kel(i1:i2,ju1:j2,ik) * precip_bot(:,:) *  &
     &        conv_frac(:,:) * tdt /  &
     &        (frac_area_con(:,:) * (grid_height(:,:,ik) * 100.0d0))

            frac_i_liq(:,:) =  &
     &        frac(:,:) / (1.0d0 + frac(:,:))

!           ---------------------------------------
!           Scavenging is limited by mass transfer.
!           ---------------------------------------

            where ((frac_i_liq(:,:) >=  &
     &              (frac_wash_con(:,:) / frac_area_con(:,:))) .and.  &
     &             washout(:,:))

!             ---------------------------
!             First do the reevaporation.
!             ---------------------------

              concentration(ic)%pArray3D(:,:,ik) =  &
     &          concentration(ic)%pArray3D(:,:,ik) +  &
     &          (1.0d0 - frac_not_evap(:,:)) * scav_tracer_con(:,:,ic) *  &
     &          (MWTAIR / mw(ic)) /  mass(:,:,ik)

              scav_tracer_con(:,:,ic) =  &
     &          scav_tracer_con(:,:,ic) * frac_not_evap(:,:)

!             ------------------------------
!             Now do the convective washout.
!             ------------------------------

              scav_tracer_con(:,:,ic) =  &
     &          scav_tracer_con(:,:,ic) +  &
     &          frac_wash_con(:,:) *  &
     &          concentration(ic)%pArray3D(:,:,ik) * mass(:,:,ik) *  &
     &          mw(ic) / MWTAIR

              concentration(ic)%pArray3D(:,:,ik) =  &
     &          concentration(ic)%pArray3D(:,:,ik) * (1.0d0 - frac_wash_con(:,:))

            end where

!           -------------------------------------------------
!           Scavenging is limited by Henry's law equilibrium.
!           -------------------------------------------------

            where ((frac_i_liq(:,:) <  &
     &              (frac_wash_con(:,:) / frac_area_con(:,:))) .and.  &
     &             washout(:,:))

              mbot(:,:) =  &
     &          frac_i_liq(:,:) *  &
     &          (frac_area_con(:,:) * concentration(ic)%pArray3D(:,:,ik) *  &
     &           mass(:,:,ik) * mw(ic) / MWTAIR +  &
     &           scav_tracer_con(:,:,ic))

              concentration(ic)%pArray3D(:,:,ik) =  &
     &          concentration(ic)%pArray3D(:,:,ik) -  &
     &          (mbot(:,:) - scav_tracer_con(:,:,ic)) /  &
     &           mass(:,:,ik) * MWTAIR / mw(ic)

              scav_tracer_con(:,:,ic) = mbot(:,:)

            end where

          end if

!         -------------------------------------------------------------
!         If this is the bottom grid box, then put the scavenged amount
!         into the wet_depos array.
!         -------------------------------------------------------------

          if ((ik == k1) .and. (pr_wet_depos)) then

            wet_depos(:,:,ic) =  &
     &        wet_depos(:,:,ic) +  &
     &        (scav_tracer_con(:,:,ic) + scav_tracer_str(:,:,ic)) /  &
     &        mcor(:,:)

          end if

!         -------------------------------------------------------------
!         Save scavenging in units kg m^{-3} s^{-1}.
!         -------------------------------------------------------------
          IF(pr_scav) THEN

            scav3d(:,:,ik,ic) = scav3d(:,:,ik,ic) +  &
     &        (scav_tracer_con(:,:,ic) + scav_tracer_str(:,:,ic)) /  &
     &        (mcor(:,:)*grid_height(:,:,ik)*tdt)

          END IF

!       =============
        end do icloop
!       =============

!     =============
      end do ikloop
!     =============


      return

      end subroutine Update_Wetdep


!-----------------------------------------------------------------------------
!
! ROUTINE
!   Calc_Wet_Loss_Rate
!
! DESCRIPTION
!   This routine calculates the wet loss rate for 1 species at each latitude.
!
! ARGUMENTS
!   convective      : convective clouds? (true or false)
!   ic              : species index for calculating wet loss rate.
!   ih2o2_num       : index of species H2O2 (zero if not present).
!   delh_298_over_r_wet : temperature scaling term for Henry's law.
!   hstar_wet           : Henry's law constant (no temperature dependence).
!   retention_eff   : retention efficiency
!   presse1d        : atmospheric pressure at the edge   of each grid box (mb)
!   kel1d           : temperature (degK)
!   kloss           : wet loss rate returned from routine (s^-1)
!
!-----------------------------------------------------------------------------

      subroutine Calc_Wet_Loss_Rate  &
     &  (convective, ic, ih2o2_num, delh_298_over_r_wet, hstar_wet,  &
     &   retention_eff, presse1d, kel1d, kloss, &
     &   i1, i2, ilo, ihi)

!      use m_PhysicalConstants, only : GAS_CONST_J

      implicit none

#     include "gmi_phys_constants.h"


!     ----------------------
!     Argument declarations.
!     ----------------------

      integer, intent(in) :: i1, i2, ilo, ihi
      logical :: convective
      integer :: ic
      integer :: ih2o2_num
      real*8  :: delh_298_over_r_wet
      real*8  :: hstar_wet
      real*8  :: retention_eff
      real*8  :: presse1d(ilo:ihi)
      real*8  :: kel1d   (ilo:ihi)
      real*8  :: kloss   (i1:i2)


!     ----------------------
!     Function declarations.
!     ----------------------

      real*8, external  :: E_Ice


!     ----------------------
!     Variable declarations.
!     ----------------------

      integer :: il

      real*8  :: sqrt19

!     ---------------------------------------------------------------------
!     ci_ice_over_ci_gas : mixing ratio in ice over mixing ratio in gas
!     ci_liq_over_ci_gas : mixing ratio in liquid over mixing ratio in gas
!     cloud_liq_water    : estimate of cloud liquid water (g/cm^3)
!     henry_eff : effective Henry's law constant (includes temperature
!                 dependence)
!     retention : retention efficiency of tracer in the liquid cloud
!                 condensate as it is converted to precipitation
!                 (retention < 1 accounts for volitilization during riming)
!     ---------------------------------------------------------------------

      real*8  :: ci_ice_over_ci_gas(i1:i2)
      real*8  :: ci_liq_over_ci_gas(i1:i2)
      real*8  :: cloud_liq_water   (i1:i2)
      real*8  :: henry_eff         (i1:i2)
      real*8  :: retention         (i1:i2)


!     ----------------
!     Begin execution.
!     ----------------

!     -------------------------------------------------------------
!     Calculate the wet scavenging species dependent rate constant.
!     -------------------------------------------------------------

      henry_eff(:) =  &
     &  hstar_wet * Exp (-delh_298_over_r_wet *  &
     &  ((1.0d0 / kel1d(i1:i2)) - (1.0d0 / 298.0d0)))

      retention(:) = 0.0d0

      where (kel1d(i1:i2) > 248.0d0)  &
     &  retention(:) = retention_eff

      where (kel1d(i1:i2) > 268.0d0)  &
     &  retention(:) = 1.0d0

      if (convective) then

        cloud_liq_water(:) = 0.0d0

        where (kel1d(i1:i2) >= 248.0d0)  &
     &    cloud_liq_water(:) =  &
     &      2.0d-6 * (kel1d(i1:i2) - 248.0d0) / 20.0d0

        where (kel1d(i1:i2) >= 268.0d0)  &
     &    cloud_liq_water(:) = 2.0d-6

      else

        cloud_liq_water(:) = 0.5d-6

      end if

      ci_liq_over_ci_gas(:) =  &
     &  henry_eff(:) * cloud_liq_water(:) *  &
     &  GAS_CONST_J / 100.0d0 * kel1d(i1:i2)

      ci_ice_over_ci_gas(:) =  0.0d0

      if (ic == ih2o2_num) then

        sqrt19 = Sqrt (1.9d0)

        do il = i1, i2

          if (kel1d(il) < 273.0d0)  &
     &      ci_ice_over_ci_gas(il) =  &
     &      (2.0d-6 - cloud_liq_water(il)) *  &
     &      presse1d(il) /  &
     &      E_Ice (Max (140.0d0, kel1d(il) )) *  &
     &      0.6d0 * sqrt19

        end do

      end if

      kloss(:) =  &
     &  (retention(:) *  &
     &   ci_liq_over_ci_gas(:) /  &
     &   (1.0d0 + ci_liq_over_ci_gas(:) + ci_ice_over_ci_gas(:)) +  &
     &   ci_ice_over_ci_gas(:) /  &
     &   (1.0d0 + ci_liq_over_ci_gas(:) + ci_ice_over_ci_gas(:)))


      return

      end subroutine Calc_Wet_Loss_Rate


!-----------------------------------------------------------------------------
!
! ROUTINE
!   E_Ice
!
! DESCRIPTION
!   This routine computes the saturation vapor pressure of ice at a given
!   Celsius temperature.
!
!   NOTE:
!     The saturation vapor pressure of ice is generated by a
!     formulation from Marti & Mauersberger (GRL, V20, No 5, pp363-366, Mar 1993.
!     Ref says good for 170K <= T <= 273. This formulation fits the old
!     formulation closely down to the previous low limit of 153K (-120C) and
!     is probably more accurate.
!
! CODE DEVELOPER
!   Stephen Steenrod, GSFC (SSAI)
!   steenrod@code613-3.gsfc.nasa.gov
!
! ARGUMENTS
!   tk : ambient temperature (degK)
!
!-----------------------------------------------------------------------------

      function E_Ice (tk)

      use GmiPrintError_mod, only : GmiPrintError

      implicit none


!     ----------------------
!     Argument declarations.
!     ----------------------

      real*8, intent(in)  :: tk


!     ----------------------
!     Function declarations.
!     ----------------------

      real*8  :: E_Ice


!     -----------------------
!     Parameter declarations.
!     -----------------------

      real*8, parameter  :: temp_min = 100.d0 ! minimum temperature (degK)
                                              ! This value can be lowered
                                              ! if needed.
      real*8, parameter  :: temp_max = 273.d0 ! maximum temperature (degK)

!     ---------------------------------------------------
!     Parameters for the fit log P = A/T + B (Marti & Mauersberger, GRL 93).
!     ---------------------------------------------------

      real*8, parameter :: A = -2663.5d0
      real*8, parameter :: B = 12.537d0

      character(len=30) :: err_msg

!     ----------------
!     Begin execution.
!     ----------------

      if( (tk .ge. temp_min) .and. (tk .le. temp_max)) then

        E_ice = 10 ** ( A/tk + B )
!... convert from Pa to hPa
        E_ice = E_ice/100

      else

!       -----------------------------------------
!       Below -120 C, stop with an error message.
!       -----------------------------------------

        Write (6, '(''Temperature: '', f10.5, '' K'')') tk - 273.0d0
        Write (6, '(''T must be in the range'', f10.5, '' <= T <= '',f10.5)')  &
     &                          temp_min, temp_max
        Write (6, '(''STOPPING in E_Ice'')')

        err_msg = ' Problem in E_Ice '
        call GmiPrintError (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)

      end if


      return

      end
