!=============================================================================
!
! $Id: emiss_update.F90,v 1.1.1.1.2.2.14.1.8.1.108.1.40.1.188.1.10.1.2.1 2020/12/04 18:48:56 mmanyin Exp $
!
! CODE DEVELOPER
!   John Tannahill, LLNL
!   jrt@llnl.gov
!
! FILE
!   emiss_update.F
!
! ROUTINES
!   Update_Emiss
!   Update_Semiss_Inchem
!
! HISTORY
!   - November 23, 2004 * Jules Kouatchou
!     Modified the routine Update_Emiss by adding the variable
!     "surf_emiss_out2" as argument of Update_Emiss and Add_Emiss_Harvard.
!   - December 8, 2005 * Bigyani Das
!     Added changes to update DMS, dust, sea salt emissions daily introducing
!     medt = 86400.0d0 and changing mdt to medt for the aerocom runs
!     when "do_aerocom" is true.
!   - January 15, 2008 * Eric Nielsen
!     Two-meter air temperature replaces surface air temperature
!=============================================================================
!
! ROUTINE
!   Update_Emiss
!
! DESCRIPTION
!   This routine updates const based on emissions.
!
! ARGUMENTS
!   lwi_flags     : 0=water; 1=land; 2=ice
!   latdeg        : latitude  (deg)
!   mcor          : area of grid box (m^2)
!   emiss_isop    : isoprene    emissions (kg/s)
!   emiss_monot   : monoterpene emissions (kg/s)
!   emiss_nox     : NOx         emissions (kg/s)
!   radswg        : net downward shortwave radiation at ground (W/m^2)
!   TwoMeter_air_temp : 2m air temperature  (degK)
!   surf_rough    : surface roughness   (m)
!   con_precip    : convective precipitation (mm/day)
!   tot_precip    : total precipitation (mm/day)
!   ustar         : ustar (m/s)
!   mass          : total mass of the atmosphere within each grid box (kg)
!   max_cloud     : maximum overlap cloud fraction for LW
!   ran_cloud     : random  overlap cloud fraction for LW
!   kel           : temperature (degK)
!   surf_emiss_out: surface emissions to be written to output (kg/m^2/time)
!   emiss_3d_out: surface emissions to be written to output (kg/m^2/time)
!   surf_emiss_out2: surface emissions to be written to output (kg/m^2/time)
!   const         : species concentration, known at zone centers (mixing ratio)
!   emiss         : array of emissions      (kg/s)
!   emiss_dust_t  : array of dust emissions (kg/s)
!   emiss_dust    : array of dust emissions for 6 hours (kg/s)
!   emiss_aero    : array of aerosol emissions (kg/s)
!   pbl           : boundary layer height (m)
!   humidity      : specific humidity (g/kg)
!   pctm1         : surface pressure at t1 (mb)
!-----------------------------------------------------------------------------

      subroutine Update_Emiss (gmiGrid, idaySoilType, firstBiogenicBase, &
     &   lwi_flags, cosSolarZenithAngle, latdeg, mcor, emiss_isop, emiss_monot,  &
     &   emiss_nox, do_ShipEmission, emiss_hno3, emiss_o3, ihno3_num, io3_num, &
     &   radswg, TwoMeter_air_temp, surf_rough, con_precip, tot_precip, ustar, &
     &   mass, fracCloudCover, kel, surf_emiss_out, surf_emiss_out2, emiss_3d_out, &
     &   aerosolEmiss3D, aerosolSurfEmiss, aerosolSurfEmissMap, concentration, &
     &   emissionArray, emiss_dust_t, emiss_dust,emiss_aero_t, emiss_aero, pbl, &
     &   gridBoxHeight, index_soil, ncon_soil, soil_fert, soil_precip, soil_pulse, &
     &   ireg, iland, iuse, convert_isop, convert_monot, coeff_isop, base_isop, &
     &   base_monot, xlai, IBOC, IBBC, INOC, IFOC, IFBC, ISSLT1, ISSLT2, ISSLT3, &
     &   ISSLT4, IFSO2, INSO2, INDMS, IAN, IMGAS, INO, iisoprene_num, ino_num, &
     &   ico_num, ipropene_num, pr_surf_emiss, pr_emiss_3d, pr_diag, loc_proc,  &
     &   rootProc, met_opt, emiss_opt, chem_opt, trans_opt, emiss_aero_opt, &
     &   emiss_dust_opt, do_aerocom, do_semiss_inchem, do_gcr, do_drydep, emiss_map, &
     &   emiss_map_dust, emiss_map_aero, ndust, nst_dust, nt_dust, naero, nymd, &
     &   num_time_steps, mw, tdt, ndt, emiss_timpyr, num_emiss, isop_scale, &
     &   i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi, i1_gl, i2_gl, ju1_gl, j2_gl, &
     &   ilong, num_species, doMEGANemission, doMEGANviaHEMCO, aefIsop, aefMbo, &
     &   aefMonot, isoLaiPrev, isoLaiCurr, isoLaiNext, pardif, pardir, T_15_AVG, &
     &   emissionSpeciesLayers, exp_fac, mixPBL)

      use GmiArrayBundlePointer_mod, only : t_GmiArrayBundle
      use GmiTimeControl_mod       , only : GmiSplitDateTime
!     use CalcAerosolEmissDiagn_mod, only : calcAerosolEmissDiagn
      use GmiGrid_mod, only : t_gmiGrid

      implicit none

#     include "gmi_emiss_constants.h"

!     ----------------------
!     Argument declarations.
!     ----------------------

      type (t_gmiGrid) , intent(in) :: gmiGrid

      integer, intent(in   ) :: IBOC, IBBC, INOC, IFOC, IFBC, ISSLT1, ISSLT2, ISSLT3, ISSLT4
      integer, intent(in   ) :: IFSO2, INSO2, INDMS, IAN, IMGAS, INO
      logical, intent(in   ) :: pr_surf_emiss, pr_emiss_3d, pr_diag, rootProc
      integer, intent(in   ) :: loc_proc
      integer, intent(in   ) :: i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi
      integer, intent(in   ) :: i1_gl, i2_gl, ju1_gl, j2_gl, ilong
      integer, intent(in   ) :: num_species, emiss_timpyr, num_emiss
      integer, intent(in   ) :: ndust, nst_dust, nt_dust, naero
      integer, intent(in   ) :: nymd, num_time_steps, ndt
      integer, intent(in   ) :: iisoprene_num, ino_num, ico_num, ipropene_num, io3_num, ihno3_num
      real*8 , intent(in   ) :: mw(num_species)
      real*8 , intent(in   ) :: isop_scale(12)
      real*8 , intent(in   ) :: tdt
      integer, intent(in   ) :: met_opt, emiss_opt, chem_opt, trans_opt, emiss_aero_opt, emiss_dust_opt
      integer, intent(in   ) :: emiss_map(num_species), emiss_map_dust(num_species)
      integer, intent(in   ) :: emiss_map_aero(num_species)
      logical, intent(in   ) :: do_semiss_inchem, do_gcr, do_drydep, do_aerocom, do_ShipEmission
      integer, intent(in   ) :: lwi_flags(i1:i2, ju1:j2)    ! 0=water; 1=land; 2=ice
      real*8 , intent(in   ) :: cosSolarZenithAngle(i1:i2, ju1:j2)
      real*8 , intent(in   ) :: latdeg   (i1:i2, ju1:j2)
      real*8 , intent(in   ) :: mcor     (i1:i2, ju1:j2)
      real*8 , intent(inout) :: emiss_isop   (i1:i2, ju1:j2)
      real*8 , intent(  out) :: emiss_monot  (i1:i2, ju1:j2)
      real*8 , intent(  out) :: emiss_nox    (i1:i2, ju1:j2)
      real*8 , intent(in   ) :: emiss_o3     (i1:i2, ju1:j2)
      real*8 , intent(in   ) :: emiss_hno3     (i1:i2, ju1:j2)
      real*8 , intent(in   ) :: radswg       (i1:i2, ju1:j2)
      real*8 , intent(in   ) :: TwoMeter_air_temp(i1:i2, ju1:j2),T_15_AVG(i1:i2, ju1:j2)
      real*8 , intent(in   ) :: surf_rough   (i1:i2, ju1:j2)
      real*8 , intent(in   ) :: con_precip   (i1:i2, ju1:j2)
      real*8 , intent(in   ) :: tot_precip   (i1:i2, ju1:j2)
      real*8 , intent(in   ) :: ustar        (i1:i2, ju1:j2)
      real*8 , intent(in   ) :: mass         (i1:i2, ju1:j2, k1:k2)
      real*8 , intent(in   ) :: kel          (ilo:ihi, julo:jhi, k1:k2)
      real*8 , intent(inout) :: surf_emiss_out (i1:i2, ju1:j2, num_species)
      real*8 , intent(inout) :: surf_emiss_out2 (i1:i2, ju1:j2, 6)
      real*8 , intent(inout) :: emiss_3d_out (i1:i2, ju1:j2, k1:k2, num_species)
      real*8 , intent(inout) :: aerosolSurfEmiss   (i1:i2, ju1:j2, ndust+naero+5)
      real*8 , intent(inout) :: aerosolEmiss3D   (i1:i2, ju1:j2, k1:k2, 5)
      integer, intent(in   ) :: aerosolSurfEmissMap(ndust+naero+5)
      real*8 , intent(in   ) :: emiss_dust_t (i1:i2, ju1:j2, ndust, nst_dust:nst_dust+nt_dust-1)
      real*8 , intent(  out) :: emiss_dust   (i1:i2, ju1:j2, ndust)
      real*8 , intent(  out) :: emiss_aero   (i1:i2, ju1:j2, naero)
      real*8 , intent(in   ) :: emiss_aero_t (i1:i2, ju1:j2, naero, emiss_timpyr)
      real*8 , intent(in   ) :: pbl          (i1:i2, ju1:j2)
      real*8 , intent(in   ) :: gridBoxHeight (i1:i2, ju1:j2, k1:k2) ! for Llnl emission
      real*8 , intent(in   ) :: fracCloudCover (i1:i2, ju1:j2) ! cloud fraction

      logical, intent(in   ) :: doMEGANemission
      logical, intent(in   ) :: doMEGANviaHEMCO
      real*8 , intent(in) :: aefIsop   (i1:i2, ju1:j2)
      real*8 , intent(in) :: aefMbo    (i1:i2, ju1:j2)
      real*8 , intent(in) :: aefMonot  (i1:i2, ju1:j2)
      real*8 , intent(in) :: isoLaiPrev(i1:i2, ju1:j2)
      real*8 , intent(in) :: isoLaiCurr(i1:i2, ju1:j2)
      real*8 , intent(in) :: isoLaiNext(i1:i2, ju1:j2)
      real*8 , intent(in) :: pardif    (i1:i2, ju1:j2)
      real*8 , intent(in) :: pardir    (i1:i2, ju1:j2)
      real*8 , intent(in) :: exp_fac(NPULSE)    ! pulsing decay per time step (day^-1)
      logical, intent(in) :: mixPBL             ! whether to explicitly distribute
                                                ! aerosol emissions within the PBL

      INTEGER, INTENT(IN   ) :: emissionSpeciesLayers(num_species)

      type (t_GmiArrayBundle), intent(inout) :: concentration(num_species)
      type (t_GmiArrayBundle), intent(inout) :: emissionArray(num_emiss)
      
     integer :: iland(i1:i2, ju1:j2, NTYPE), ireg(i1:i2, ju1:j2)
     integer :: iuse (i1:i2, ju1:j2, NTYPE)
     real*8  :: convert_isop (NVEGTYPE), convert_monot(NVEGTYPE)
     real*8  :: base_isop (i1:i2, ju1:j2, NTYPE), base_monot(i1:i2, ju1:j2, NTYPE)
     real*8  :: xlai      (i1:i2, ju1:j2, NTYPE)
     real*8  :: soil_fert(i1:i2, ju1:j2), soil_precip(i1:i2, ju1:j2), soil_pulse(NPULSE+1,i1:i2, ju1:j2)
     integer :: ncon_soil (NVEGTYPE), index_soil(2, i1:i2, ju1:j2)
     real*8  :: coeff_isop   (NPOLY)
     integer, intent(inOut) :: idaySoilType
     logical, intent(inOut) :: firstBiogenicBase
      
!     ----------------------
!     Variable declarations.
!     ----------------------

      integer :: idumyr, month, idumday
      integer :: num_new_emiss_dust
      real*8  :: medt

      real*8  :: tempk(i1:i2, ju1:j2)

!     ----------------
!     Begin execution.
!     ----------------

      if (pr_diag) then
        Write (6,*) 'Update_Emiss called by ', loc_proc
      end if

!     -----------------------------------
!     Update dust emission every 6 hours.
!     -----------------------------------

      if (emiss_dust_opt == 1) then    ! GMI dust emissions

        if (do_aerocom) then
           medt = 86400.0d0
        else
           medt = 21600.0d0
        end if

        num_new_emiss_dust =  &
     &    (num_time_steps / (Nint (medt) / ndt)) + 1

        if (Mod (num_time_steps, (Nint (medt) / ndt)) == 0) then

          emiss_dust(i1:i2,ju1:j2,1:ndust) =  &
     &      emiss_dust_t(i1:i2,ju1:j2,1:ndust,num_new_emiss_dust)

        end if

      end if


      if (btest(emiss_opt,0) .or. btest(emiss_opt,1)) then

!       ===================
        call Add_Emiss_Llnl  &
!       ===================
     &    (do_aerocom, pr_surf_emiss, pr_emiss_3d, mcor, surf_emiss_out, emiss_3d_out,  &
     &     mass, concentration, emissionArray, emiss_dust,  emiss_aero_t, emiss_aero, pbl, &
     &     gridBoxHeight, IBOC, IBBC, INOC, IFOC, IFBC, ISSLT1, ISSLT2, ISSLT3, ISSLT4, &
     &     IFSO2, INSO2, INDMS, pr_diag, loc_proc, chem_opt, trans_opt, emiss_aero_opt, &
     &     emiss_dust_opt, do_semiss_inchem, emiss_map, emiss_map_dust, emiss_map_aero, &
     &     emissionSpeciesLayers, ndust, naero, nymd, mw, tdt, emiss_timpyr, num_emiss, &
     &     i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi, num_species, mixPBL)
 
        ! Compute surface emission diagnostics for aerosols.

!        if ((pr_surf_emiss) .or. (pr_emiss_3d)) then
!           if ((emiss_aero_opt > 0) .or. (emiss_dust_opt > 0)) then
!              call calcAerosolEmissDiagn(aerosolEmiss3D, aerosolSurfEmiss, &
!                       aerosolSurfEmissMap, &
!                       emissionArray, emiss_dust, emiss_aero, tdt, mcor, &
!                       pr_surf_emiss, pr_emiss_3d, &
!                       emiss_aero_opt, emiss_dust_opt, emiss_map, num_species, &
!                       ndust, naero, i1, i2, ju1, j2, k1, k2, num_emiss)
!           end if
!        end if

        if (btest(emiss_opt,1)) then

          IF ( .not. doMEGANviaHEMCO ) THEN 
             emiss_isop (:,:) = 0.0d0 
          END IF
          emiss_monot(:,:) = 0.0d0
          emiss_nox  (:,:) = 0.0d0

          if (met_opt == 3) then
            tempk(i1:i2,ju1:j2) = TwoMeter_air_temp(i1:i2,ju1:j2)
          else
            tempk(i1:i2,ju1:j2) = kel(i1:i2,ju1:j2,1)
          end if

          call Update_Emiss_Harvard (gmiGrid, idaySoilType, firstBiogenicBase, &
     &       iisoprene_num, ino_num, lwi_flags, tdt, mw, cosSolarZenithAngle, latdeg, nymd, &
     &       mcor, tempk, radswg, surf_rough, con_precip, tot_precip, ustar, fracCloudCover,  &
     &       emiss_isop, emiss_monot, emiss_nox, index_soil, ncon_soil, soil_fert, soil_precip, &
     &       soil_pulse, ireg, iland, iuse, convert_isop, convert_monot, coeff_isop,  &
     &       base_isop, base_monot, xlai, pr_diag, loc_proc, rootProc, i1, i2, ju1, j2, &
     &       k1, k2, i1_gl, i2_gl, ju1_gl, j2_gl, ilong, num_species, doMEGANemission, &
     &       doMEGANviaHEMCO, aefIsop, aefMbo, aefMonot, isoLaiPrev, isoLaiCurr, isoLaiNext, &
             pardif, pardir,T_15_AVG, exp_fac)

!         ------------------------------------------------
!         Scale the isoprene emissions on a monthly basis.
!         Default is no scaling (isop_scale = 1.0d0).
!         ------------------------------------------------

          call GmiSplitDateTime (nymd, idumyr, month, idumday)

          emiss_isop(:,:) = emiss_isop(:,:) * isop_scale(month)


          if (.not. do_semiss_inchem) then
!           ======================
            call Add_Emiss_Harvard  &
!           ======================
     &        (pr_surf_emiss, pr_emiss_3d, iisoprene_num, ico_num,  &
     &         ipropene_num,  &
     &         ino_num, tdt, mw, mcor, surf_emiss_out, surf_emiss_out2,  &
     &         emiss_3d_out, mass, concentration, emiss_isop, emiss_monot, emiss_nox, &
     &         do_ShipEmission, emiss_hno3, emiss_o3, ihno3_num, io3_num, &   
     &         pr_diag, loc_proc, i1, i2, ju1, j2, k1, k2, num_species)

          end if


        end if

      end if

!=================================================
!   Do Galactic Cosmic Ray Emissions of N and NO
!=================================================

      if (do_gcr) then
        if (pr_diag) then
          Write (6,*) 'Add_Emiss_Gsfc called by ', loc_proc
        end if
!       ===================
        call Add_Emiss_Gsfc  &
!       ===================
     &    (mass, concentration, nymd, tdt, &
     &     IAN, IMGAS, INO, &
     &     pr_diag, loc_proc, i1, i2, ju1, j2, k1, k2, num_species)
      end if


      return

      end


!-----------------------------------------------------------------------------
!
! ROUTINE
!   Update_Semiss_Inchem
!
! DESCRIPTION
!   This routine updates const based on surface emissions, and does it in the
!   chemistry solver.
!
! ARGUMENTS
!   emiss_isop  : isoprene    emissions (kg/s)
!   emiss_monot : monoterpene emissions (kg/s)
!   emiss_nox   : NOx         emissions (kg/s)
!   mcor        : area of grid box  (m^2)
!   humidity    : specific humidity (g/kg)
!   pctm2       : CTM surface pressure at t1+tdt (mb)
!   kel         : temperature (degK)
!   emissionArray: array of emissions (kg/s)
!
!   surf_emiss  : surface  emissions (units?)
!
!-----------------------------------------------------------------------------

      subroutine Update_Semiss_Inchem  &
     &  (pr_surf_emiss, pr_emiss_3d, emiss_isop, emiss_monot, emiss_nox,  &
     &   do_ShipEmission, emiss_hno3, emiss_o3, ihno3_num, io3_num, &   
     &   mcor, surf_emiss_out, surf_emiss_out2, emiss_3d_out, &
     &   emissionArray, surf_emiss, surf_emiss2, emiss_3d, &
     &   gridBoxHeight, emiss_timpyr, num_emiss, emiss_opt, emiss_map, tdt, nymd, &
     &   ico_num, ino_num, ipropene_num, iisoprene_num, mw, &
     &   pr_diag, loc_proc, i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi, num_species)

      use MAPL_ConstantsMod
      use GmiTimeControl_mod, only : GmiSplitDateTime
      use GmiArrayBundlePointer_mod, only : t_GmiArrayBundle

      implicit none

#     include "gmi_phys_constants.h"
#     include "gmi_time_constants.h"


!     ----------------------
!     Argument declarations.
!     ----------------------

      logical, intent(in) :: pr_diag
      integer, intent(in) :: loc_proc, i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi, num_species
      integer, intent(in) :: num_emiss, emiss_timpyr, emiss_opt
      integer, intent(in) :: emiss_map(num_species)
      real*8 , intent(in) :: tdt, mw(num_species)
      integer, intent(in) :: nymd
      integer, intent(in) :: ico_num, ino_num, ipropene_num, iisoprene_num
      integer, intent(in) :: io3_num, ihno3_num
      logical :: pr_surf_emiss, do_ShipEmission
      logical :: pr_emiss_3d
      real*8 , intent(in   ) :: gridBoxHeight(i1:i2, ju1:j2, k1:k2)
!                             ! height of each grid box (m)
      real*8  :: emiss_isop     (i1:i2, ju1:j2)
      real*8  :: emiss_monot    (i1:i2, ju1:j2)
      real*8  :: emiss_nox      (i1:i2, ju1:j2)
      real*8  :: emiss_hno3     (i1:i2, ju1:j2)
      real*8  :: emiss_o3    (i1:i2, ju1:j2)
      real*8  :: mcor           (i1:i2, ju1:j2)
      real*8  :: surf_emiss_out (i1:i2, ju1:j2, num_species)
      real*8  :: surf_emiss_out2 (i1:i2, ju1:j2, 6)
      real*8  :: emiss_3d_out (i1:i2, ju1:j2, k1:k2,num_species)
      real*8  :: surf_emiss (i1:i2, ju1:j2, num_species)
      real*8  :: surf_emiss2 (i1:i2, ju1:j2, 6)
      real*8  :: emiss_3d(i1:i2, ju1:j2, k1:k2, num_species)

      type(t_GmiArrayBundle), intent(inout) :: emissionArray(num_emiss)

!     -----------------------
!     Parameter declarations.
!     -----------------------

      real*8, parameter :: CMPM3 = CMPM * CMPM * CMPM

!     ----------------------
!     Variable declarations.
!     ----------------------

      integer :: ic, icx
      integer :: ij, il, ik
      integer :: idumday, idumyear
      integer :: imon
      integer :: inum
      integer :: month

      real*8  :: conv_emiss(i1:i2, ju1:j2)


!     ----------------
!     Begin execution.
!     ----------------

      if (pr_diag) then
        Write (6,*) 'Update_Semiss_Inchem called by ', loc_proc
      end if


      if (emiss_timpyr == MONTHS_PER_YEAR) then
        call GmiSplitDateTime (nymd, idumyear, month, idumday)
        imon = month
      else
        imon = 1
      end if


!     --------------------------------------------------------------
!     Calculate conversion factor to go from kg/box/s to molecules/cm^3/s,
!     but leave out the molecular weight term for each species.
!     --------------------------------------------------------------

!     NOTE: MAPL_AVOGAD is [molec/kmol]
      conv_emiss(:,:) =  MAPL_AVOGAD / mcor(:,:) / CMPM3 / gridBoxheight(:,:,k1)

      inum = 0

      surf_emiss(:,:,:) = 0.0d0
      surf_emiss2(:,:,:) = 0.0d0
      emiss_3d(:,:,:,:) = 0.0d0

!     ================================
      SPCLOOP: do icx = 1, num_species
!     ================================

        ic = emiss_map(icx)

        if (ic > 0) then

          inum = inum + 1

          IF(inum > num_emiss) THEN
           PRINT *,"Update_Semiss_Inchem: Asking for more emissions than are available."
           STOP
          END IF

          surf_emiss(:,:,ic) =  &
     &      surf_emiss(:,:,ic) +  &
     &      emissionArray(inum)%pArray3D(:,:,k1) * conv_emiss(:,:) / mw(ic)

           emiss_3d(:,:,k1,ic) =  &
     &      emiss_3d(:,:,k1,ic) +  &
     &      emissionArray(inum)%pArray3D(:,:,k1) * conv_emiss(:,:) / mw(ic)

!                                ============
          if (inum == num_emiss) exit SPCLOOP
!                                ============

        end if

!     ==============
      end do SPCLOOP
!     ==============


      if (emiss_opt == 2) then

!       ==============================
        call Add_Semiss_Harvard_Inchem  &
!       ==============================
     &    (iisoprene_num, ico_num, ipropene_num, ino_num, mw,  &
     &     emiss_isop, emiss_monot, emiss_nox, &
     &     do_ShipEmission, emiss_hno3, emiss_o3, ihno3_num, io3_num, &
     &     conv_emiss, surf_emiss,  &
     &     surf_emiss2, emiss_3d, &
     &     pr_diag, loc_proc, i1, i2, ju1, j2, k1, k2, num_species)

      end if


!     -------------------------------------------------------------
!     Accumulate the diagnostics surf_emiss_out and emiss_3d_out.
!      1. Convert the emission arrays back to kg s^{-1}.
!      2. Divide by the tile surface area.
!      3. Multiply by the time step length.
!     The result is (should be!) kg m^{-2} dt^{-1}.
!     -------------------------------------------------------------

      if (pr_surf_emiss) then

        do ic = 1, num_species

          surf_emiss_out(:,:,ic) =  &
     &      surf_emiss_out(:,:,ic) +  &
     &      surf_emiss(:,:,ic) / conv_emiss(:,:) * mw(ic) *  &
     &      tdt / mcor(:,:)
        end do

          surf_emiss_out2(:,:,1) = surf_emiss_out2(:,:,1) +  &
     &      surf_emiss2(:,:,1) / conv_emiss(:,:) * mw(ico_num) *  &
     &      tdt / mcor(:,:)
          surf_emiss_out2(:,:,2) = surf_emiss_out2(:,:,2) +  &
     &      surf_emiss2(:,:,2) / conv_emiss(:,:) * mw(ico_num) *  &
     &      tdt / mcor(:,:)
          surf_emiss_out2(:,:,3) = surf_emiss_out2(:,:,3) +  &
     &      surf_emiss2(:,:,3) / conv_emiss(:,:) * mw(ipropene_num) *  &
     &      tdt / mcor(:,:)

          surf_emiss_out2(:,:,4) = surf_emiss_out2(:,:,4) +  &
     &      surf_emiss2(:,:,4) / conv_emiss(:,:) * mw(ino_num) *  &
     &      tdt / mcor(:,:)
     
          if (do_ShipEmission) then

             surf_emiss_out2(:,:,5) = surf_emiss_out2(:,:,5) +  &
     &         surf_emiss2(:,:,5) / conv_emiss(:,:) * mw(ihno3_num) *  &
     &         tdt / mcor(:,:)

             surf_emiss_out2(:,:,6) = surf_emiss_out2(:,:,6) +  &
     &         surf_emiss2(:,:,6) / conv_emiss(:,:) * mw(io3_num) *  &
     &         tdt / mcor(:,:)

          end if
       end if

       if (pr_emiss_3d) then
        do ic = 1, num_species

             emiss_3d_out(:,:,k1,ic) =  &
     &          emiss_3d_out(:,:,k1,ic) +  &
     &          emiss_3d(:,:,k1,ic) / conv_emiss(:,:) * mw(ic) *  &
     &          tdt / mcor(:,:)

        end do

      end if

      return

      end

