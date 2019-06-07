!=============================================================================
!
! $Id$
!
! CODE DEVELOPER
!   originated from GRANTOUR programmed by John Walton
!   modified by Xiaohong Liu
!
! FILE
!   sulfchem.F
!
! ROUTINES
!   Do_Sulfchem
!
! HISTORY
!   - December 8, 2005 - Bigyani Das
!     Added two variables do_aerocom, do_dust_emiss to the argument list
!     of Do_Sulfchem, Do_Gas_Sulfchem and Do_Aqu_Sulfchem
!=============================================================================

!-----------------------------------------------------------------------------
!
! ROUTINE
!   Do_Sulfchem
!
! DESCRIPTION
!
!     *****************************************************************
!     * carries out the sulfur chemistry calculation
!     *****************************************************************
!
!     Performs a numerical integration of the full aqueous chemistry
!      equations for the species H2SO4, DMS, SO2 & H2O2 in cloudy air.
!     Performs an analytic calculation of the clear sky part of the
!      chemistry equations for the species H2SO4, DMS, SO2 & H2O2 in
!      clear air.
!     Combines the clear and cloudy concentrations to simulate mixing
!      between clear and cloudy regions on the basis of cloud fraction
!
!     ********** input:
!     *     dtchem  =  model time step in seconds
!     *     dtclc   =  time step for clear air chemistry
!     *     dtclmx  =  max. time step between mixing calculations in seconds
!     *     dtcl    =  time step between mixing calculations in seconds
!     *     dtaq    =  time step used in integration of the aqueous
!                      equations in seconds
!     *     dph     =  assumed pH of the water drops
!
!     *     tgcm       air temperature (K))
!     *     zmair      air density (# cm-3)
!     *     gmair      air density (gm cm-3)
!     *     oh         OH  concentration (cm-3)
!     *     ho2        HO2 concentration (cm-3)
!     *     no3        NO3 concentration (cm-3)
!     *     o3         O3  concentration (cm-3)
!     *     cfac       total cloud fraction (0-1)
!     *     cwac       cloud water mixing ratio for aq chem (g/g) (in-cloud)
!
!     *     qh1p       DMS + OH  reaction rate coefficient
!     *     qh2p       SO2 + OH     "      "        "
!     *     qh3p       HO2 + HO2    "      "        "
!     *     qj4p       H2O2 + hv photolysis rate coef.
!     *     qh5p       H2O2 + OH reaction rate coefficient
!     *     qh6p       SO2(aq) + H2O2(aq)  "        "
!     *     qh7p       SO2(aq) + O3(aq)    "        "
!     *     qh8p       HO2 + HO2 + H2O     "        "
!     *     qh9p       DMS + NO3 reaction rate coefficient
!
!     *     a1->a9     reaction coefficients
!
!
!     ********** output:
!
!     * Here is how I treat separate changes due to the aqueous
!       and clear sky chemical paths.
!
!       I assume that aqueous and clear sky chemistry act on the
!       total concentration, C, giving the concentration changes
!       dCa(C) & dCn(C).
!
!       The new values of the concentration at t+dt are
!
!           C(t+dt) = C(t) + cfp * dCa(C) + ( 1 - cfp ) * dCn(C)
!
!       In the case of H2SO4, these two terms contribute separately
!       to aqeuous and clear sky concentration.
!
!       In the case of all clear air or all cloud, these reduce to
!       what one would expect.
!
!       Definitions
!
!       C           species concentration
!
!       dCa, dCn    Concentration changes for aqueous & clear sky
!                   chemistry.
!
!       cfp         fractional cloud cover.
!
!     *----------------------------------------------------------------
!
! ARGUMENTS
!   itloop    : # of zones (ilong * ilat * ivert)
!   time      : time interval for each chemical calculation (seconds)
!   tgcm      : temperature (K)
!   pgcm      : atmospheric pressure at the center of each grid box (mb)
!   xr        : species concentration, known at zone centers
!               (molecules/cm^3 for gas-species and kg/kg for aerosols)
!   humidity  : specific humidity (g/kg)
!   semiss    : array of sulfur emissions (molecules/cm^3/s)
!   qjgmi     : photolysis rate constants (s^-1)
!   qkgmi     : thermal    rate constants (units vary)
!   cfac      : total cloud fraction [0 - 1]
!   cwac      : in-cloud liquid water in each grid box (g/g)
!   relh      : relative humidity [0 - 1]
!   aqua_infile_name : aquachem input file name
!
!-----------------------------------------------------------------------------

      subroutine Do_Sulfchem  &
     &  (itloop, time, tgcm, pgcm, xr, humidity,  &
     &   semiss, qjgmi, qkgmi, cfac, cwac, relh,  &
     &   aqua_infile_name, pr_diag, loc_proc, num_time_steps,  &
     &   massc, pr_sulf_src, do_aerocom, do_dust_emiss,  &
     &   dms_oh, dms_no3, so2_oh, so2_h2o2, so2_o3)

      implicit none

#     include "gmi_phys_constants.h"
#     include "setkin_par.h"
#     include "sulfchem.h"

!     ----------------------
!     GMIMOD-related values.
!     ----------------------

      character*(*) :: aqua_infile_name

      integer :: itloop
      real*8  :: time
      real*8  :: tgcm    (itloop)
      real*8  :: pgcm    (itloop)
      real*8  :: xr      (itloop, NSP)
      real*8  :: humidity(itloop)
      real*8  :: semiss  (itloop, NSP)
      real*8  :: qjgmi   (itloop, NUM_J)
      real*8  :: qkgmi   (itloop, NUM_K)
      real*8  :: cfac    (itloop)
      real*8  :: cwac    (itloop)
      real*8  :: relh    (itloop)

      logical :: pr_diag
      integer :: loc_proc
      integer :: num_time_steps

      real*8  :: massc   (itloop)
      logical :: pr_sulf_src
      logical :: do_aerocom
      logical :: do_dust_emiss
      real*8  :: dms_oh  (itloop)
      real*8  :: dms_no3 (itloop)
      real*8  :: so2_oh  (itloop)
      real*8  :: so2_h2o2(itloop)
      real*8  :: so2_o3  (itloop)


! ----------------------------------------------------------
!  *** LOCAL  VARIABLES WERE DECLARED HERE.
! ----------------------------------------------------------

      integer :: numcl

      real*8   dtclmx, dtaqset,  &
     &         accset, dph

      parameter (dtclmx = 14400.d0)
      parameter (dtaqset = 600.d0)
      parameter (accset = 0.0001d0)
      parameter (dph = 4.5d0)

      real*8  &
     &  dtaq,    & ! time step used in integration of aqueous equations in seconds
     &  dtchem,  & ! model time step in seconds
     &  dtcl,    & ! time step between mixing calculations in seconds
     &  dtclc  ! time step for clear air chemistry

!     qh6       SO2(aq) + H2O2(aq) reaction rate coefficient
      real*8  qh6

!     air density (molecules/cm3).
      real*8 :: zmair(itloop)

!     air density (g/cm3).
      real*8 :: gmair(itloop)

!     H2O concentration in molecules/cm3.
      real*8 :: h2o(itloop)

!     concentrations of O3, OH, HO2 and NO3 in molecules/cm3.
      real*8 :: o3(itloop)
      real*8 :: oh(itloop)
      real*8 :: ho2(itloop)
      real*8 :: no3(itloop)

!     rtcd (s-1) rate coef. for SO2 -> dust transfer
!     rtcs (s-1) rate coef. for SO2 -> sslt transfer
      real*8 :: rtcd(itloop,4)
      real*8 :: rtcs(itloop,4)

!     cp and dcp have units [molecules/cm3 for gas-species and
!     gm spec./gm air for aerosols]
      real*8 :: cp (itloop,NSP)
      real*8 :: dcp(itloop,NSP)

!     H2O2 photolysi rate coefficient in 1/s
      real*8 :: qj4p(itloop)

!     qh1p       DMS + OH  reaction rate coefficient
!     qh2p       SO2 + OH     "      "        "
!     qh3p       HO2 + HO2    "      "        "
!     qh5p       H2O2 + OH reaction rate coefficient
!     qh7p       SO2(aq) + O3(aq)    "        "
!     qh8p       HO2 + HO2 + H2O     "        "
!     qh9p       DMS + NO3    "      "        "

      real*8 :: qh1p(itloop)
      real*8 :: qh2p(itloop)
      real*8 :: qh3p(itloop)
      real*8 :: qh5p(itloop)
      real*8 :: qh7p(itloop)
      real*8 :: qh8p(itloop)
      real*8 :: qh9p(itloop)

!     henry's law coefficients (hwp, hyp, h3p) for SO2, H2O2, O3
      real*8 :: hwp(itloop)
      real*8 :: hyp(itloop)
      real*8 :: h3p(itloop)

!     array of grid numbers for grids in clouds
      integer :: idxcl(itloop)

!     a1-9 coefficients in aqueous equations (only for cloudy grids)
      real*8 :: acoef(itloop,9)


      zmair = 0.0d0
      gmair = 0.0d0
      h2o   = 0.0d0
      o3    = 0.0d0
      oh    = 0.0d0
      ho2   = 0.0d0
      no3   = 0.0d0
      rtcd  = 0.0d0
      rtcs  = 0.0d0
      cp    = 0.0d0
      dcp   = 0.0d0

      qj4p  = 0.0d0
      qh1p  = 0.0d0
      qh2p  = 0.0d0
      qh3p  = 0.0d0
      qh5p  = 0.0d0
      qh7p  = 0.0d0
      qh8p  = 0.0d0
      qh9p  = 0.0d0
      hwp   = 0.0d0
      hyp   = 0.0d0
      h3p   = 0.0d0
      idxcl = 0.0d0
      acoef = 0.0d0


!     ----------------
!     Begin execution.
!     ----------------

      if (pr_diag) then
        Write (6,*) 'Do_Sulfchem called by ', loc_proc
      end if


!     *****************************************************************
!     * Set chemistry parameters
!     *****************************************************************

      dtchem = time
      dtcl   = dtclmx
      dtaq   = dtaqset

      hpp  = 10**(-dph)
      rhpp = 1./hpp

!     *****************************************************************
!     * Compute air density, "zmair" [#/cm3] and "gmair" [gm/cm3] at
!       model grid points.
!     *****************************************************************

      do ijk = 1, itloop
        zmair(ijk) = xr(ijk, IMGAS)
        gmair(ijk) = ( MWTAIR / AVOGAD ) * zmair(ijk)
      end do

!     *****************************************************************
!     * Calculate water vapor concentration in molecules/cm^3
!     *****************************************************************

      do ijk = 1, itloop
        h2o(ijk) = humidity(ijk) * MWTAIR / MWTH2O / GPKG  &
     &             * xr(ijk, IMGAS)
      end do

      do ijk = 1, itloop
        o3(ijk)  = xr(ijk, IO3)
        oh(ijk)  = xr(ijk, IOH)
        ho2(ijk) = xr(ijk, IHO2)
        no3(ijk) = xr(ijk, INO3)
      end do

!     change unit of pressure at the grid points from millibar to bar
      do ijk = 1, itloop
        pgcm(ijk) = pgcm(ijk) / MBPB
      end do

!     *****************************************************************
!     * Define the H2O2 photolysis rate coefficient
!     *****************************************************************

      do ijk = 1, itloop
        qj4p(ijk) = qjgmi(ijk, 1)
      end do

!     *****************************************************************
!     * Given by GMIMOD concentrations for species changed by
!       sulfur chemistry
!     *****************************************************************

      do ic = 1, NSP
        cp(:, ic) = xr(:, ic)
      end do


!     *****************************************************************
!     * Calculate the SO2 mass transfer rate coefficients for
!       dust and sea salt
!     *****************************************************************

!          =================
!     call Do_So2toPar_Rates
!          =================
!    &  (itloop, pgcm, tgcm,
!    &   relh, aqua_infile_name,
!    &   gmair, cp, rtcd, rtcs)


!     *****************************************************************
!     * Calculate gas-phase reaction rates
!       (qh1p, qh2p, qh3p, qh5p, qh8p, qh9p)
!     *****************************************************************

!          ===============
      call Do_QK_Gas_Rates  &
!          ===============
     &  (itloop, tgcm, zmair,  &
     &   qh1p, qh2p, qh3p, qh5p, qh8p, qh9p)


!     ================================================================
      if (Mod (num_time_steps, (Nint(dtcl) / Nint(dtchem))) == 0) then
!     ================================================================

!     *****************************************************************
!     * Set up for aqueous chemistry
!       (generate numcl, idxcl)
!     *****************************************************************

!          ================
      call Do_Setup_Aquchem  &
!          ================
     &  (itloop, numcl, idxcl, cfac, cwac)


      if (numcl .gt. 0) then

!     *******************************************************************
!     * There are some grid points in clouds. Set up for aqueous chem.
!       calculate the aqueous-phase reaction rates, effective Henry's law
!       coefficients, and aqueous equation coeff. (a1-a9)
!     *******************************************************************

!            ===============
        call Do_QK_Aqu_Rates  &
!            ===============
     &    (itloop, tgcm, hpp, rhpp, numcl, idxcl,  &
     &     oh, ho2, h2o, no3, o3,  &
     &     cwac, zmair,  &
     &     qh1p, qh2p, qh3p, qj4p, qh5p, qh8p, qh9p,  &
     &     qh6, qh7p, hwp, hyp, h3p,  &
     &     acoef)

        dtclc = dtcl

      else

!     ******************************************************************
!     * There are no cloudy grid points, advance the entire dtchem step.
!     * Advance the clear air chemistry only.
!     ******************************************************************

        dtclc = dtcl   ! dtchem

      endif

!     =====
      endif
!     =====

!.... time step over cloud mixing time (dtclc)

!     tcl = 0.0d0

!     ========
! 20  continue
!     ========


!     ***********************************************************
!     * Calculate the species changes due to clear sky chemistry.
!     ***********************************************************

!          ===============
      call Do_Gas_Sulfchem  &
!          ===============
!    &  (dtclc, itloop, cp, dcp, cfac, zmair, semiss,
     &  (dtchem, itloop, cp, dcp, cfac, zmair, semiss,  &
     &   oh, ho2, h2o, no3,  &
     &   massc, pr_sulf_src, do_aerocom, do_dust_emiss,  &
     &   dms_oh, dms_no3, so2_oh,  &
     &   qh1p, qh2p, qh3p, qj4p, qh5p, qh8p, qh9p)


!     =================================================================
      if (Mod (num_time_steps, (Nint(dtcl) / Nint(dtchem))) == 0) then
!     =================================================================

!     *****************************************************************
!     * Integrate aqueous equations over dtclc and add the aqueous
!       change to the clear sky changes already computed.
!     *****************************************************************
!
      if (numcl .gt. 0) then

!            ===============
        call Do_Aqu_Sulfchem  &
!            ===============
     &    (accset, dtaq, dtchem, dtcl, dtclc,  &
     &     itloop, numcl, idxcl, semiss,  &
     &     massc, pr_sulf_src,  do_aerocom, do_dust_emiss,  &
     &     dms_oh, dms_no3, so2_oh, so2_h2o2, so2_o3,  &
     &     cp, dcp, cfac, zmair, acoef, loc_proc)

      end if

!     ======
      end if
!     ======


!     *****************************************************************
!     * Impose the clear plus cloudy sky concentration changes on
!       the grid cells.
!     *****************************************************************

!          ====================
      call Do_Const_Sulf_Update  &
!          ====================
     &  (itloop, cp, dcp)


!     *****************************************************************
!     * Calculate the mass transfer of SO2 vapor to water on the
!       surface of particles.
!     * The SO2 will appear on the particles as SO4.
!     *****************************************************************

!          =================
!     call Do_Const_So2toPar
!          =================
!    &  (dtchem, itloop, zmair,
!    &   rtcd, rtcs, cp)


!.... advance the time

!     tcl = tcl + dtclc
!                                 ========
!     if (tcl .lt. dtchem - 0.01d0) go to 20
!                                 ========

!     *****************************************************************
!     * Give back to GMIMOD data
!     *****************************************************************

      do ic = 1, NSP
        xr(:, ic) = cp(:, ic)
      end do


      return

      end

