!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!          & Atmospheric Chemistry & Dynamics, Code 614                  !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  WetRemovalGmiMod
!
! !INTERFACE:
!

   module  WetRemovalGmiMod

! !USES:

   implicit none

! !PUBLIC TYPES:
!
   PRIVATE

!
! !PUBLIC MEMBER FUNCTIONS:
!

   PUBLIC  TR_GMI_WetRemoval

!
! !DESCRIPTION:
!
!  This module implements in-cloud scavenging and below-cloud washout;
!   the routines were originally in GMI, in the file wetdep_update.F90 .
!
! !ORIGINAL COMMENTS:
!
!   CODE DEVELOPER
!     Dan Bergmann, LLNL
!     dbergmann@llnl.gov
!   FILE
!     wetdep_update.F
!   HISTORY
!     * March 16, 2005 - Jules Kouatchou
!       In the function E_Ice:
!          - Added two parameters "temp_min" (minimum temperature value) and
!            "temp_max" (maximum temperature value, equal to 273 degK).
!          - Lowered the value of "temp_min" to 100 degK (initially set to 153 degK)
!            to prevent the code from crashing. Such a low value is acceptable
!            because low temperature occur in the lower mesoshpere that does
!            not affect much GMI simulations.
!
! !REVISION HISTORY:
!
!   7Aug2014 - Michael Manyin, first crack
!
!              Derived from GMI code in wetdep_update.F90
!                Deleted the MICRO_AEROSOL (U Mich) compile flags
!                Changed met fields to real*4
!                Changed to top-down k-index
!                Reduced to 1 species of tracer
!                Replaced scheme of indexing into aerosol.h file w/ explicit parameters
!                Made some return parameters optional
!                Eliminated redundant parameters (e.g. julo)
!
!EOP

!-------------------------------------------------------------------------
CONTAINS

!-----------------------------------------------------------------------------
!
! ROUTINE
!   TR_GMI_WetRemoval  (adapted from the routine Update_Wetdep)
!
! DESCRIPTION
!   This routine calculates both in cloud scavenging and below cloud washout,
!   based on a writeup by Daniel Jacob at Harvard.  This is loosely based on
!   Giorgi and Chameides JGR 91 14,367-14376, 1986.
!
! ARGUMENTS
!   h2o2_flag    : if the species should be treated as H2O2
!   hno3_flag    : if the species should be treated as HNO3
!   aero_flag    : treat as an aerosol for the purposes of washout and reevaporation
!   tdt          : model time step  (s)
!   mw           : molecular weight of species (g/mol)
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
!
!-----------------------------------------------------------------------------

      subroutine TR_GMI_WetRemoval                                   &
     &  (i1, i2, j1, j2, k1, k2, tdt, mw,                            &
     &   rel_scav_eff, hstar, delH_298_over_R, retention_eff,        &
     &   aero_flag, h2o2_flag, hno3_flag,                            &
     &   mcor, grid_height, mass, moistq,                            &
     &   rain_cn, rain_ls,                                           &
     &   kel, press3c, press3e, concentration,                       &
     &   wet_depos, scav3d)

      implicit none

#     include "gmi_phys_constants.h"
#     include "gmi_time_constants.h"

!     ----------------------
!     Argument declarations.
!     ----------------------

      integer         , intent(in   ) :: i1, i2, j1, j2, k1, k2

      real*4          , intent(in   ) :: tdt
      real*4          , intent(in   ) :: mw

      ! Parameters which are listed per species in gmi_aerosol.h and setkin_depos.h
      real*4          , intent(in   ) :: rel_scav_eff    ! scavenging efficiency relative to sulfate (unitless)
      real*4          , intent(in   ) :: hstar           ! a parameter for computing Henry's law constant
      real*4          , intent(in   ) :: delH_298_over_R ! a parameter for computing Henry's law constant
      real*4          , intent(in   ) :: retention_eff   ! retention efficiency

      logical         , intent(in   ) :: aero_flag   ! treat as an aerosol for the purposes of washout and reevaporation
      logical         , intent(in   ) :: h2o2_flag   ! if the species should be treated as H2O2
      logical         , intent(in   ) :: hno3_flag   ! if the species should be treated as HNO3

      real*4          , intent(in   ) :: mcor          (i1:i2,  j1:j2)
      real*4          , intent(in   ) :: grid_height   (i1:i2,  j1:j2,  k1:k2)
      real*4          , intent(in   ) :: mass          (i1:i2,  j1:j2,  k1:k2)
      real*4          , intent(in   ) :: moistq        (i1:i2,  j1:j2,  k1:k2)

      real*4          , intent(in   ) :: rain_cn       (i1:i2,  j1:j2,  k1:k2)
      real*4          , intent(in   ) :: rain_ls       (i1:i2,  j1:j2,  k1:k2)

      real*4          , intent(in   ) :: kel           (i1:i2,  j1:j2,  k1:k2)
      real*4          , intent(in   ) :: press3c       (i1:i2,  j1:j2,  k1:k2)
      real*4          , intent(in   ) :: press3e       (i1:i2,  j1:j2,  k1:k2+1) ! edge i is above center i

      real*4          , intent(inout) :: concentration (i1:i2,  j1:j2,  k1:k2)

      real*4, optional, intent(inout) :: wet_depos     (i1:i2,  j1:j2)
      real*4, optional, intent(inout) :: scav3d        (i1:i2,  j1:j2,  k1:k2)



!     -----------------------
!     Parameter declarations.
!     -----------------------

      real*4, parameter :: CM2PM2    = CMPM * CMPM  ! cm^2/m^2
      real*4, parameter :: CMPMM     = 0.1d0        ! cm/mm
      real*4, parameter :: WATER_DEN = 1.0d0        ! gr/cm^3

!     ----------------------
!     Variable declarations.
!     ----------------------

!     ---------------------------------------------------------
!     rainout : logical array pointing to locations for rainout
!     washout : logical array pointing to locations for washout
!     ---------------------------------------------------------

      logical :: rainout(i1:i2, j1:j2)
      logical :: washout(i1:i2, j1:j2)

      integer :: ij, ik

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

      real*8  :: conv_frac      (i1:i2, j1:j2)
      real*8  :: exp_fac        (i1:i2, j1:j2)
      real*8  :: frac           (i1:i2, j1:j2)
      real*8  :: frac_area_con  (i1:i2, j1:j2)
      real*8  :: frac_area_str  (i1:i2, j1:j2)
      real*8  :: frac_i_liq     (i1:i2, j1:j2)
      real*8  :: frac_not_evap  (i1:i2, j1:j2)
      real*8  :: frac_wash_con  (i1:i2, j1:j2)
      real*8  :: frac_wash_str  (i1:i2, j1:j2)
      real*8  :: henry_eff      (i1:i2, j1:j2)
      real*8  :: kloss          (i1:i2, j1:j2)
      real*8  :: kloss_con      (i1:i2, j1:j2)
      real*8  :: kloss_str      (i1:i2, j1:j2)
      real*8  :: mbot           (i1:i2, j1:j2)
      real*8  :: precip_bot     (i1:i2, j1:j2)
      real*8  :: precip_top     (i1:i2, j1:j2)
      real*8  :: q              (i1:i2, j1:j2)

!     --------------------------------------------------------------
!     scav_tracer_con : amount of scavenged tracer due to convective
!                       processes (kg)
!     scav_tracer_str : amount of scavenged tracer due to stratiform
!                       processes (kg)
!     --------------------------------------------------------------

      real*8  :: scav_tracer_con(i1:i2, j1:j2)
      real*8  :: scav_tracer_str(i1:i2, j1:j2)

!     ----------------
!     Begin execution.
!     ----------------

      IF( present(wet_depos) ) THEN
       wet_depos(:,:) = 0.0d0
      END IF
      IF( present(scav3d) ) THEN
       scav3d(:,:,:) = 0.0d0
      END IF

      henry_eff (:,:) = 0.0d0

      kloss_con (:,:) = 0.0d0
      kloss_str (:,:) = 0.0d0

      precip_bot(:,:) = 0.0d0

      scav_tracer_con(:,:) = 0.0d0
      scav_tracer_str(:,:) = 0.0d0


!     -----------------------------------------------------------------
!     Find the wet scavenged fraction using a first order rainout
!     parameterization from Giorgi and Chameides
!     (JGR v91 pp14,367-14376 1986).  Use this for both the large-scale
!     and the convective precipitation.
!     -----------------------------------------------------------------

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
      ikloop: do ik = k1, k2          !  TOP-DOWN
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
     &    (press3c(:,:,ik) * MWTAIR * MB2CGS) /  &
     &    (GPKG * SECPDY * kel(:,:,ik) * BOLTZMN_E * AVOGAD)

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
     &          (SECPDY * mcor(:,:) * CM2PM2)))

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

        washout(:,:) = ((moistq(:,:,ik) >=   0.0d0) .and.  &
     &                    (kel (:,:,ik) >= 268.0d0) .and.  &
     &              (precip_top(:,:)    >    0.0d0))

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

!       =======
!!!     icloop: do ic = 1, num_species
        kloss_str(:,:) = 0.0d0
        kloss_con(:,:) = 0.0d0

!       =======

          if (hno3_flag) then

            where (q(:,:) > 0.0d0)
              kloss_con(:,:) = k_con
              kloss_str(:,:) = k_min +  &
     &          q(:,:) * (1.0d0 - conv_frac(:,:)) / l_plus_w_str
            elsewhere
              kloss_con(:,:) = 0.0d0
              kloss_str(:,:) = 0.0d0
            end where

          else if (rel_scav_eff /= 0.0e0) then

            where (q(:,:) > 0.0d0)
              kloss_con(:,:) = k_con * rel_scav_eff
              kloss_str(:,:) = (k_min +  &
     &          q(:,:) * (1.0d0 - conv_frac(:,:)) / l_plus_w_str) &
     &          * rel_scav_eff
            elsewhere
              kloss_con(:,:) = 0.0d0
              kloss_str(:,:) = 0.0d0
            end where

          else if (hstar > 0.0d0) then

            do ij = j1, j2

!-micro_aerosol : note gmi version of code is kept, the following is called
!                 twice, in first call kloss_con is used and in second call
!                 kloss_str is used
!             =======================
              call Calc_Wet_Loss_Rate  &
!             =======================
     &          (.true., h2o2_flag, delH_298_over_R, hstar,  &
     &           retention_eff, press3e(:,ij,ik),  &
     &           kel(:,ij,ik), kloss_con(:,ij),  &
     &           i1, i2 )

!             =======================
              call Calc_Wet_Loss_Rate  &
!             =======================
     &          (.false., h2o2_flag, delH_298_over_R, hstar,  &
     &           retention_eff, press3e(:,ij,ik),  &
     &           kel(:,ij,ik), kloss_str(:,ij),  &
     &           i1, i2 )

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

          end if

!         -----------------------------------------------------------
!         Only apply rainout in grid boxes where precipitation forms.
!         -----------------------------------------------------------

          where (rainout(:,:))

!           --------------------------
!           Do the stratiform rainout.
!           --------------------------

            exp_fac(:,:) = Exp (-kloss_str(:,:) * tdt)

            scav_tracer_str(:,:) =  &
     &        scav_tracer_str(:,:) +  &
     &        frac_area_str(:,:) * (1.0d0 - exp_fac(:,:)) *  &
     &        concentration(:,:,ik) * mass(:,:,ik) *  &
     &        mw / MWTAIR

            concentration(:,:,ik) =  &
     &        concentration(:,:,ik) *  &
     &        (1.0d0 -  &
     &         (frac_area_str(:,:) * (1.0d0 - exp_fac(:,:))))

!           --------------------------
!           Do the convective rainout.
!           --------------------------

            exp_fac(:,:) = Exp (-kloss_con(:,:) * tdt)

            scav_tracer_con(:,:) =  &
     &        scav_tracer_con(:,:) +  &
     &        frac_area_con(:,:) * (1.0d0 - exp_fac(:,:)) *  &
     &        concentration(:,:,ik) * mass(:,:,ik) *  &
     &        mw / MWTAIR

            concentration(:,:,ik) =  &
     &        concentration(:,:,ik) *  &
     &        (1.0d0 -  &
     &         (frac_area_con(:,:) * (1.0d0 - exp_fac(:,:))))

          end where

!         -------------------------------------------------------------
!         Do the washout and reevaporation of aerosols and nitric acid.
!         -------------------------------------------------------------

          if (aero_flag .or. hno3_flag) then

            where (washout(:,:))

              concentration(:,:,ik) =  &
     &          concentration(:,:,ik) + &
     &          (1.0d0 - frac_not_evap(:,:)) *  &
     &          (scav_tracer_str(:,:) + scav_tracer_con(:,:)) *  &
     &          (MWTAIR / mw) / mass(:,:,ik)

!mkch

              scav_tracer_str(:,:) =  &
     &          scav_tracer_str(:,:) * frac_not_evap(:,:)

              scav_tracer_con(:,:) =  &
     &          scav_tracer_con(:,:) * frac_not_evap(:,:)

!             ------------------------------
!             Now do the stratiform washout.
!             ------------------------------

              scav_tracer_str(:,:) =  &
     &          scav_tracer_str(:,:) + &
     &          frac_wash_str(:,:) * &
     &          concentration(:,:,ik) * mass(:,:,ik) *  &
     &          mw / MWTAIR

              concentration(:,:,ik) = &
     &          concentration(:,:,ik) * (1.0d0 - frac_wash_str(:,:))

!             ------------------------------
!             Now do the convective washout.
!             ------------------------------

              scav_tracer_con(:,:) =  &
     &          scav_tracer_con(:,:) +  &
     &          frac_wash_con(:,:) * &
     &          concentration(:,:,ik) * mass(:,:,ik) *  &
     &          mw / MWTAIR

              concentration(:,:,ik) = &
     &          concentration(:,:,ik) * (1.0d0 - frac_wash_con(:,:))

            end where

          else

!           ---------------------------------------------
!           Do washout and reevaporation of species other
!           than aerosols or nitric acid.
!           ---------------------------------------------

            henry_eff(:,:) =  &
     &        hstar *  &
     &        Exp (-delH_298_over_R *  &
     &             ((1.0d0 / kel(:,:,ik)) -  &
     &              (1.0d0 / 298.0d0)))

!           -----------------------------------------------
!           Stratiform or large scale component of washout.
!           -----------------------------------------------

            frac(:,:) =  &
     &        henry_eff(:,:) * GAS_CONST_J / 100.0d0 *  &
     &        kel(:,:,ik) * precip_bot(:,:) *  &
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

              concentration(:,:,ik) =  &
     &          concentration(:,:,ik) +  &
     &          (1.0d0 - frac_not_evap(:,:)) * scav_tracer_str(:,:) *  &
     &          (MWTAIR / mw) /  mass(:,:,ik)

              scav_tracer_str(:,:) =  &
     &          scav_tracer_str(:,:) * frac_not_evap(:,:)

!             ------------------------------
!             Now do the stratiform washout.
!             ------------------------------

              scav_tracer_str(:,:) =  &
     &          scav_tracer_str(:,:) +  &
     &          frac_wash_str(:,:) *  &
     &          concentration(:,:,ik) * mass(:,:,ik) *  &
     &          mw / MWTAIR

              concentration(:,:,ik) =  &
     &          concentration(:,:,ik) * (1.0d0 - frac_wash_str(:,:))

            end where

!           -------------------------------------------------
!           Scavenging is limited by Henry's law equilibrium.
!           -------------------------------------------------

            where ((frac_i_liq(:,:) <  &
     &              (frac_wash_str(:,:) / frac_area_str(:,:))) .and.  &
     &             washout(:,:))

              mbot(:,:) =  &
     &          frac_i_liq(:,:) *  &
     &          (frac_area_str(:,:) * concentration(:,:,ik) *  &
     &           mass(:,:,ik) * mw / MWTAIR +  &
     &           scav_tracer_str(:,:))

              concentration(:,:,ik) =  &
     &          concentration(:,:,ik) -  &
     &          (mbot(:,:) - scav_tracer_str(:,:)) /  &
     &          mass(:,:,ik) * MWTAIR / mw

              scav_tracer_str(:,:) = mbot(:,:)

            end where

!           --------------------------------
!           Convective component of washout.
!           --------------------------------

            frac(:,:) =  &
     &        henry_eff(:,:) * GAS_CONST_J / 100.0d0 *  &
     &        kel(:,:,ik) * precip_bot(:,:) *  &
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

              concentration(:,:,ik) =  &
     &          concentration(:,:,ik) +  &
     &          (1.0d0 - frac_not_evap(:,:)) * scav_tracer_con(:,:) *  &
     &          (MWTAIR / mw) /  mass(:,:,ik)

              scav_tracer_con(:,:) =  &
     &          scav_tracer_con(:,:) * frac_not_evap(:,:)

!             ------------------------------
!             Now do the convective washout.
!             ------------------------------

              scav_tracer_con(:,:) =  &
     &          scav_tracer_con(:,:) +  &
     &          frac_wash_con(:,:) *  &
     &          concentration(:,:,ik) * mass(:,:,ik) *  &
     &          mw / MWTAIR

              concentration(:,:,ik) =  &
     &          concentration(:,:,ik) * (1.0d0 - frac_wash_con(:,:))

            end where

!           -------------------------------------------------
!           Scavenging is limited by Henry's law equilibrium.
!           -------------------------------------------------

            where ((frac_i_liq(:,:) <  &
     &              (frac_wash_con(:,:) / frac_area_con(:,:))) .and.  &
     &             washout(:,:))

              mbot(:,:) =  &
     &          frac_i_liq(:,:) *  &
     &          (frac_area_con(:,:) * concentration(:,:,ik) *  &
     &           mass(:,:,ik) * mw / MWTAIR +  &
     &           scav_tracer_con(:,:))

              concentration(:,:,ik) =  &
     &          concentration(:,:,ik) -  &
     &          (mbot(:,:) - scav_tracer_con(:,:)) /  &
     &           mass(:,:,ik) * MWTAIR / mw

              scav_tracer_con(:,:) = mbot(:,:)

            end where

          end if

!         -------------------------------------------------------------
!         If this is the bottom grid box, then put the scavenged amount
!         into the wet_depos array.
!         -------------------------------------------------------------

          if ((ik == k2) .and. present(wet_depos)) then

            wet_depos(:,:) =  &
     &        wet_depos(:,:) +  &
     &        (scav_tracer_con(:,:) + scav_tracer_str(:,:)) /  &
     &        mcor(:,:)

          end if

!         -------------------------------------------------------------
!         Save scavenging in units kg m^{-3} s^{-1}.
!         -------------------------------------------------------------
          IF( present(scav3d) ) THEN

            scav3d(:,:,ik) = scav3d(:,:,ik) +  &
     &        (scav_tracer_con(:,:) + scav_tracer_str(:,:)) /  &
     &        (mcor(:,:)*grid_height(:,:,ik)*tdt)

          END IF

!       =============
!!!     end do icloop
!       =============

!     =============
      end do ikloop
!     =============


      return

      end subroutine TR_GMI_WetRemoval

!=============================================================================



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
!   h2o2_flag       : if the species should be treated as H2O2
!   delH_298_over_R : temperature scaling term for Henry's law.
!   hstar           : Henry's law constant (no temperature dependence).
!   retention_eff   : retention efficiency
!   presse1d        : atmospheric pressure at the edge   of each grid box (mb)
!   kel1d           : temperature (degK)
!   kloss           : wet loss rate returned from routine (s^-1)
!
!-----------------------------------------------------------------------------

      subroutine Calc_Wet_Loss_Rate  &
     &  (convective, h2o2_flag, delH_298_over_R, hstar,  &
     &   retention_eff, presse1d, kel1d, kloss, &
     &   i1, i2)

!      use m_PhysicalConstants, only : GAS_CONST_J

      implicit none

#     include "gmi_phys_constants.h"


!     ----------------------
!     Argument declarations.
!     ----------------------

      integer, intent(in) :: i1, i2
      logical :: convective
      logical :: h2o2_flag
      real*4  :: delH_298_over_R
      real*4  :: hstar
      real*4  :: retention_eff
      real*4  :: presse1d(i1:i2)
      real*4  :: kel1d   (i1:i2)
      real*8  :: kloss   (i1:i2)


!     ----------------------
!     Variable declarations.
!     ----------------------

      integer :: il

      real*4  :: sqrt19

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

      real*4  :: ci_ice_over_ci_gas(i1:i2)
      real*4  :: ci_liq_over_ci_gas(i1:i2)
      real*4  :: cloud_liq_water   (i1:i2)
      real*4  :: henry_eff         (i1:i2)
      real*4  :: retention         (i1:i2)


!     ----------------
!     Begin execution.
!     ----------------

!     -------------------------------------------------------------
!     Calculate the wet scavenging species dependent rate constant.
!     -------------------------------------------------------------

      henry_eff(:) =  &
     &  hstar * Exp (-delH_298_over_R *  &
     &  ((1.0d0 / kel1d(:)) - (1.0d0 / 298.0d0)))

      retention(:) = 0.0d0

      where (kel1d(:) > 248.0d0)  &
     &  retention(:) = retention_eff

      where (kel1d(i1:i2) > 268.0d0)  &
     &  retention(:) = 1.0d0

      if (convective) then

        cloud_liq_water(:) = 0.0d0

        where (kel1d(i1:i2) >= 248.0d0)  &
     &    cloud_liq_water(:) =  &
     &      2.0d-6 * (kel1d(:) - 248.0d0) / 20.0d0

        where (kel1d(:) >= 268.0d0)  &
     &    cloud_liq_water(:) = 2.0d-6

      else

        cloud_liq_water(:) = 0.5d-6

      end if

      ci_liq_over_ci_gas(:) =  &
     &  henry_eff(:) * cloud_liq_water(:) *  &
     &  GAS_CONST_J / 100.0d0 * kel1d(:)

      ci_ice_over_ci_gas(:) =  0.0d0

      if (h2o2_flag) then

        sqrt19 = Sqrt (1.9d0)

        do il = i1, i2

          if (kel1d(il) < 273.0d0)  &
     &      ci_ice_over_ci_gas(il) =  &
     &      (2.0d-6 - cloud_liq_water(il)) *  &
     &      presse1d(il) /  &
     &      E_Ice (Max (140.0e0, kel1d(il) )) *  &
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

      real*4, intent(in)  :: tk


!     ----------------------
!     Function declarations.
!     ----------------------

      real*4  :: E_Ice


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

      end function E_Ice

   end module WetRemovalGmiMod
