
!=============================================================================
!
! $Id$
!
! CODE DEVELOPER
!   John Tannahill, LLNL
!   jrt@llnl.gov
!
! FILE
!   emiss_harvard.F
!
! ROUTINES
!   Update_Emiss_Harvard
!   Do_Biogenic_Emiss
!   Do_Soil_Emiss
!   Add_Emiss_Harvard
!   Add_Semiss_Harvard_Inchem
!
! HISTORY
!   - November 23, 2004 * Jules Kouatchou
!     Added the variable "surf_emiss_out2" as argument of the routine
!     Add_Emiss_Harvard. It is used to produce surface emission diagnostics
!     for soil NO, monoterpenes CO, methanol CO and biogenic propene.
!   - August 19, 2005 * Jules Kouatchou
!     Set here the value of NLANDHAR using the function
!     NLANDHAR_expected.
!   - January 15, 2008 * Eric Nielsen
!     Two-meter air temperature replaces surface air temperature
!=============================================================================


!-----------------------------------------------------------------------------
!
! ROUTINE
!   Update_Emiss_Harvard
!
! DESCRIPTION
!   This routine updates the Harvard biogenic (isoprene, CO, & propene) and
!   soil emissions (NOx).
!
! HISTORY
!   * August 12, 2005 * Jules Kouatchou
!     Allocate the variables index_soil, soil_fert, soil_precip, and
!     soil_pulse the first time this routine is called.
!
! ARGUMENTS
!   iisoprene_num : index of isoprene in const
!   ino_num       : index of NO       in const
!   nhms          : hour/minute/sec  (HHMMSS)
!   nymd          : year/month/day   (YYYYMMDD)
!   lwi_flags     : array of flags that indicate land, water, or ice
!   tdt           : model time step  (s)
!   mw            : array of species' molecular weights (g/mol)
!   latdeg        : latitude         (deg)
!   londeg        : longitude        (deg)
!   mcor          : area of grid box (m^2)
!   tempk         : surface or 2m air temperature  (degK)
!   pardif        : Diffuse photosynthetically active radiation (0.35-0.70 um)
!   pardir        : Direct  photosynthetically active radiation (0.35-0.70 um)
!   radswg        : net downward shortwave radiation at ground (W/m^2)
!   surf_rough    : surface roughness        (m)
!   con_precip    : convective precipitation (mm/day)
!   tot_precip    : total precipitation      (mm/day)
!   ustar         : ustar                    (m/s)
!   max_cloud     : maximum overlap cloud fraction for LW
!   ran_cloud     : random  overlap cloud fraction for LW
!   emiss_isop    : isoprene    emissions    (kg/s)
!   emiss_monot   : monoterpene emissions    (kg/s)
!   emiss_nox     : NOx         emissions    (kg/s)
!
!-----------------------------------------------------------------------------

      SUBROUTINE Update_Emiss_Harvard (gmiGrid, idaySoilType, firstBiogenicBase,  &
     &   iisoprene_num, ino_num, lwi_flags, tdt, mw, cosSolarZenithAngle, latdeg, &
         nymd, mcor, tempk, radswg, surf_rough, con_precip, tot_precip, ustar,    &
         fracCloudCover, emiss_isop, emiss_monot, emiss_nox, index_soil,          &
         ncon_soil, soil_fert, soil_precip, soil_pulse, ireg, iland, iuse,        &
         convert_isop, convert_monot, coeff_isop, base_isop, base_monot, xlai,    &
         pr_diag, loc_proc, rootProc, i1, i2, ju1, j2, k1, k2, i1_gl, i2_gl,      &
         ju1_gl, j2_gl, ilong, num_species, doMEGANemission, aefIsop, aefMbo,     &
         aefMonot, isoLaiPrev, isoLaiCurr, isoLaiNext, pardif, pardir,T_15_AVG,   &
         exp_fac)

      USE GmiGrid_mod, ONLY : t_gmiGrid
      USE ReadInputMegan_mod, ONLY : setMEGANisoLAI
      USE GmiEmissionMEGAN_mod, ONLY : calcBiogenicMEGANemission

      implicit none

#     include "gmi_emiss_constants.h"
#     include "gmi_time_constants.h"
#     include "gmi_phys_constants.h"


!     ----------------------
!     Argument declarations.
!     ----------------------

     type (t_gmiGrid) , intent(in) :: gmiGrid
     logical, intent(in) :: pr_diag, rootProc, doMEGANemission
     integer, intent(in) :: loc_proc, nymd
     integer, intent(in) :: i1, i2, ju1, j2, k1, k2
     integer, intent(in) :: i1_gl, i2_gl, ju1_gl, j2_gl, ilong, num_species
     integer, intent(inOut) :: idaySoilType
     logical, intent(inOut) :: firstBiogenicBase
     real*8 , intent(in) :: aefIsop   (i1:i2, ju1:j2)
     real*8 , intent(in) :: aefMbo    (i1:i2, ju1:j2)
     real*8 , intent(in) :: aefMonot  (i1:i2, ju1:j2)
     real*8 , intent(in) :: isoLaiPrev(i1:i2, ju1:j2)
     real*8 , intent(in) :: isoLaiCurr(i1:i2, ju1:j2)
     real*8 , intent(in) :: isoLaiNext(i1:i2, ju1:j2)
     real*8 , intent(in) :: pardif    (i1:i2, ju1:j2)
     real*8 , intent(in) :: pardir    (i1:i2, ju1:j2)
     integer :: iisoprene_num, ino_num
     real*8  :: latdeg     (i1:i2, ju1:j2), cosSolarZenithAngle (i1:i2, ju1:j2)
     real*8, intent(in) :: exp_fac(NPULSE)    ! pulsing decay per time step (day^-1)
     integer :: lwi_flags  (i1:i2, ju1:j2)
     real*8  :: tdt
     real*8  :: mw         (num_species)  , mcor       (i1:i2, ju1:j2)
     real*8  :: tempk      (i1:i2, ju1:j2), radswg     (i1:i2, ju1:j2)
     real*8  :: surf_rough (i1:i2, ju1:j2), con_precip (i1:i2, ju1:j2)
     real*8  :: tot_precip (i1:i2, ju1:j2), ustar      (i1:i2, ju1:j2)
     real*8  :: fracCloudCover  (i1:i2, ju1:j2)
     real*8  :: emiss_isop (i1:i2, ju1:j2), emiss_monot(i1:i2, ju1:j2)
     real*8  :: emiss_nox  (i1:i2, ju1:j2)
     integer :: iland(i1:i2, ju1:j2, NTYPE), ireg(i1:i2, ju1:j2)
     integer :: iuse (i1:i2, ju1:j2, NTYPE)
     real*8  :: convert_isop (NVEGTYPE), convert_monot(NVEGTYPE)
     real*8  :: base_isop (i1:i2, ju1:j2, NTYPE), base_monot(i1:i2, ju1:j2, NTYPE)
     real*8  :: xlai      (i1:i2, ju1:j2, NTYPE)
     real*8  :: soil_fert(i1:i2, ju1:j2), soil_precip(i1:i2, ju1:j2), soil_pulse(NPULSE+1,i1:i2, ju1:j2)
     integer :: ncon_soil (NVEGTYPE), index_soil(2, i1:i2, ju1:j2)
     real*8  :: coeff_isop   (NPOLY)
     integer :: days_btw_m ! days between midmonths in the LAI data
     real*8  :: isoLai(i1:i2, ju1:j2)

     LOGICAL :: PT_15isOK
     REAL*8, INTENT(IN)  :: T_15_AVG(i1:i2, ju1:j2)

!     ----------------------
!     Variable declarations.
!     ----------------------

      integer :: ij
      integer :: nsec_jan1

      real*8  :: decl, tdtinv
      real*8  :: rdistsq
      real*8  :: days, time

!     ----------------
!     Begin execution.
!     ----------------

      if (pr_diag) then
        Write (6,*) 'Update_Emiss_Harvard called by ', loc_proc
      end if

       if (doMEGANemission) then

          days_btw_m = 31
	  isoLai(i1:i2, ju1:j2) = 0.00
	  
          PT_15isOK = .TRUE.

          call setMEGANisoLAI &
     &           (isoLai, isoLaiCurr, isoLaiPrev, isoLaiNext, &
     &            days_btw_m, nymd, i1, i2, ju1, j2, i1_gl, ju1_gl)

!         ======================
          call calcBiogenicMEGANemission  &
!         ======================
     &      (emiss_isop, emiss_monot, days_btw_m, &
     &       aefIsop, aefMbo, aefMonot, isoLai, isoLaiCurr, isoLaiPrev, &
     &       tempk, T_15_AVG, pardif, pardir, cosSolarZenithAngle, &
     &       pr_diag, loc_proc, i1, i2, ju1, j2, PT_15isOK)

          ! Perform unit conversions
          tdtinv = 1.0d0 / tdt

          emiss_isop = emiss_isop * tdtinv / ATOMSC_PER_MOLECISOP *  &
     &                        (mw(iisoprene_num) / AVOGAD) * KGPG

          emiss_monot = emiss_monot * tdtinv / ATOMSC_PER_MOLECMONOT *  &
     &                          (MWTMONOT / AVOGAD) * KGPG

       else
!         ======================
          call Do_Biogenic_Emiss  &
!         ======================
     &      (iisoprene_num, ireg, iland, iuse, tdt, mw, cosSolarZenithAngle, mcor,  &
     &       tempk, fracCloudCover, coeff_isop, convert_isop, convert_monot,  &
     &       xlai, base_isop, base_monot, emiss_isop, emiss_monot, firstBiogenicBase, &
     &       pr_diag, loc_proc, rootProc, i1, i2, ju1, j2, num_species)
      end if

!     ==================
      call Do_Soil_Emiss (gmiGrid, &
!     ==================
     &   ino_num, nymd, lwi_flags, ireg, iland, iuse, index_soil,  &
     &   ncon_soil, tdt, mw, latdeg, cosSolarZenithAngle, mcor, tempk, radswg,  &
     &   surf_rough, con_precip, tot_precip, ustar, fracCloudCover,  &
     &   soil_fert, soil_precip, soil_pulse, xlai, emiss_nox, idaySoilType, &
     &   pr_diag, loc_proc, i1, i2, ju1, j2, ju1_gl, j2_gl, num_species, exp_fac)

      return

      end


!-----------------------------------------------------------------------------
!
! ROUTINE
!   Do_Biogenic_Emiss
!
! DESCRIPTION
!   This routine updates the Harvard biogenic emissions (isoprene, CO, &
!   propene).
!
! ARGUMENTS
!   iisoprene_num : index of isoprene in const
!   ireg   : number of land types in a grid square
!   iland  : land type id in grid square for ireg land types
!   iuse   : fraction of grid box area occupied by land type (mils?)
!   tdt    : model time step  (s)
!   mw     : array of species' molecular weights (g/mol)
!   cossza : cosines of the solar zenith angle
!   mcor   : area of grid box (m^2)
!   tempk  : surface or 2m air temperature (degK)
!   fracCloudCover    : fractional cloud cover
!   coeff_isop    : coefficients used for polynomial fit
!   convert_isop  : isoprene    emissions by landtype   (atomsC/cm^2 leaf/s)
!   convert_monot : monoterpene emissions by landtype   (atomsC/cm^2 leaf/s)
!   xlai          : leaf area index of land type for current month
!   base_isop     : baseline emissions for isoprene     (kgC/box/step?)
!   base_monot    : baseline emissions for monoterpenes (kgC/box/step?)
!   emiss_isop    : isoprene    emissions (kg/s)
!   emiss_monot   : monoterpene emissions (kg/s)
!
!-----------------------------------------------------------------------------

      subroutine Do_Biogenic_Emiss  &
     &  (iisoprene_num, ireg, iland, iuse, tdt, mw, cossza, mcor,  &
     &   tempk, fracCloudCover, coeff_isop, convert_isop, convert_monot,  &
     &   xlai, base_isop, base_monot, emiss_isop, emiss_monot, firstBiogenicBase, &
     &   pr_diag, loc_proc, rootProc, i1, i2, ju1, j2, num_species)

      implicit none

#     include "gmi_phys_constants.h"
#     include "gmi_emiss_constants.h"

!     ----------------------
!     Argument declarations.
!     ----------------------

      logical, intent(in) :: pr_diag, rootProc
      integer, intent(in) :: loc_proc, i1, i2, ju1, j2, num_species
      integer :: iisoprene_num
      integer :: ireg (i1:i2, ju1:j2)
      integer :: iland(i1:i2, ju1:j2, NTYPE)
      integer :: iuse (i1:i2, ju1:j2, NTYPE)
      real*8  :: tdt
      real*8  :: mw   (num_species)
      real*8  :: cossza       (i1:i2, ju1:j2)
      real*8  :: mcor         (i1:i2, ju1:j2)
      real*8  :: tempk        (i1:i2, ju1:j2)
      real*8  :: fracCloudCover   (i1:i2, ju1:j2)
      real*8  :: coeff_isop   (NPOLY)
      real*8  :: convert_isop (NVEGTYPE)
      real*8  :: convert_monot(NVEGTYPE)
      real*8  :: xlai         (i1:i2, ju1:j2, NTYPE)
      real*8  :: base_isop    (i1:i2, ju1:j2, NTYPE)
      real*8  :: base_monot   (i1:i2, ju1:j2, NTYPE)
      real*8  :: emiss_isop   (i1:i2, ju1:j2)
      real*8  :: emiss_monot  (i1:i2, ju1:j2)
      logical, intent(inOut) :: firstBiogenicBase

!     ----------------------
!     Function declarations.
!     ----------------------

      real*8, external  :: Biogenic_Isop
      real*8, external  :: Biogenic_Monot

!     ----------------------
!     Variable declarations.
!     ----------------------

      integer :: il, ij

      integer :: iuse1(NTYPE)

      real*8  :: biemiss, bmemiss
      real*8  :: tdtinv

      real*8  :: base_isop1 (NTYPE)
      real*8  :: base_monot1(NTYPE)
      real*8  :: xlai1      (NTYPE)


!     ----------------
!     Begin execution.
!     ----------------

      if (pr_diag) then
        Write (6,*) 'Do_Biogenic_Emiss called by ', loc_proc
      end if


!     ==========
      if (firstBiogenicBase) then
!     ==========

        firstBiogenicBase = .false.

!       ==================
        call Biogenic_Base  &
!       ==================
     &    (ireg, iland, tdt, convert_isop, convert_monot, mcor,  &
     &     base_isop, base_monot, pr_diag, loc_proc, i1, i2, ju1, j2)

      end if


      tdtinv = 1.0d0 / tdt


      do ij = ju1, j2
        do il = i1, i2

          iuse1(:) = iuse(il,ij,:)

          base_isop1 (:) = base_isop (il,ij,:)
          base_monot1(:) = base_monot(il,ij,:)
          xlai1      (:) = xlai      (il,ij,:)


          biemiss =  &
!           =============
     &      Biogenic_Isop  &
!           =============
     &        (ireg(il,ij), iuse1(:), fracCloudCover(il,ij), cossza(il,ij),  &
     &         tempk(il,ij), coeff_isop(:), base_isop1(:), xlai1(:))

          emiss_isop(il,ij) =  &
     &      biemiss * tdtinv / ATOMSC_PER_MOLECISOP *  &
     &      (mw(iisoprene_num) / AVOGAD) * KGPG


          bmemiss =  &
!           ==============
     &      Biogenic_Monot  &
!           ==============
     &        (ireg(il,ij), iuse1(:), tempk(il,ij), base_monot1(:),  &
     &         xlai1(:))

          emiss_monot(il,ij) =  &
     &      bmemiss * tdtinv / ATOMSC_PER_MOLECMONOT *  &
     &      (MWTMONOT / AVOGAD) * KGPG

        end do
      end do

      return

      end


!-----------------------------------------------------------------------------
!
! ROUTINE
!   Do_Soil_Emiss
!
! DESCRIPTION
!   This routine updates the Harvard soil emissions (NOx)
!
! ARGUMENTS
!   ino_num       : index of NO in const
!   nymd          : year/month/day  (YYYYMMDD)
!   lwi_flags     : array of flags that indicate land, water, or ice
!   ireg          : number of land types in a grid square
!   iland         : land type id in grid square for ireg land types
!   iuse          : fraction of grid box area occupied by land type (mil^-1?)
!   index_soil    : i,j of the grid
!   ncon_soil     : Olson -> soil type
!   tdt           : model time step (s)
!   mw            : array of species' molecular weights (g/mol)
!   latdeg        : latitude (deg)
!   cossza        : cosines of the solar zenith angle
!   mcor          : area of grid box         (m^2)
!   tempk         : surface or 2m air temperature  (degK)
!   radswg        : net downward shortwave radiation at ground (W/m^2)
!   surf_rough    : surface roughness        (m)
!   con_precip    : convective precipitation (mm/day)
!   tot_precip    : total precipitation      (mm/day)
!   ustar         : ustar                    (m/s)
!   fracCloudCover    : fractional cloud cover
!   soil_fert     : fertilizers   (ng N/m^2/s)
!   soil_precip   : two months of observed precipitation (mm/day/box)
!   soil_pulse    : tracking of wet/dry & three types of pulsing (Y&L, 94)
!   xlai          : leaf area index of land type for month #1
!   emiss_nox     : NOx emissions (kg/s)
!
!-----------------------------------------------------------------------------

      subroutine Do_Soil_Emiss (gmiGrid,  &
     &   ino_num, nymd, lwi_flags, ireg, iland, iuse, index_soil,  &
     &   ncon_soil, tdt, mw, latdeg, cossza, mcor, tempk, radswg,  &
     &   surf_rough, con_precip, tot_precip, ustar, fracCloudCover,  &
     &   soil_fert, soil_precip, soil_pulse, xlai, emiss_nox, idaySoilType,   &
     &   pr_diag, loc_proc, i1, i2, ju1, j2, ju1_gl, j2_gl, num_species, exp_fac)

      use GmiGrid_mod, only : t_gmiGrid
      use GmiSoilEmission_mod

      implicit none

#     include "gmi_phys_constants.h"
#     include "gmi_emiss_constants.h"


!     ----------------------
!     Argument declarations.
!     ----------------------

      type (t_gmiGrid) , intent(in) :: gmiGrid
      logical, intent(in) :: pr_diag
      integer, intent(in) :: loc_proc, i1, i2, ju1, j2, ju1_gl, j2_gl, num_species
      integer :: ino_num
      integer :: nymd
      integer :: lwi_flags  (i1:i2, ju1:j2)
      integer :: ireg       (i1:i2, ju1:j2)
      integer :: iland      (i1:i2, ju1:j2, NTYPE)
      integer :: iuse       (i1:i2, ju1:j2, NTYPE)
      integer :: index_soil (2, i1:i2, ju1:j2)
      integer :: ncon_soil  (NVEGTYPE)
      real*8, intent(in) :: exp_fac(NPULSE) ! pulsing decay per time step (day^-1)
      integer, intent(inOut) :: idaySoilType
      real*8  :: tdt
      real*8  :: mw         (num_species)
      real*8  :: latdeg     (i1:i2, ju1:j2)
      real*8  :: cossza     (i1:i2, ju1:j2)
      real*8  :: mcor       (i1:i2, ju1:j2)
      real*8  :: tempk      (i1:i2, ju1:j2)
      real*8  :: radswg     (i1:i2, ju1:j2)
      real*8  :: surf_rough (i1:i2, ju1:j2)
      real*8  :: con_precip (i1:i2, ju1:j2)
      real*8  :: tot_precip (i1:i2, ju1:j2)
      real*8  :: ustar      (i1:i2, ju1:j2)
      real*8  :: fracCloudCover (i1:i2, ju1:j2)
      real*8  :: soil_fert  (i1:i2, ju1:j2)
      real*8  :: soil_precip(i1:i2, ju1:j2)
      real*8  :: soil_pulse (NPULSE+1, i1:i2, ju1:j2)
      real*8  :: xlai       (i1:i2, ju1:j2, NTYPE)
      real*8  :: emiss_nox  (i1:i2, ju1:j2)


!     ----------------------
!     Variable declarations.
!     ----------------------

      integer :: il, ij

      real*8  :: cmpm2
      real*8  :: rfac
      real*8  :: windsp

!     ----------------------------------------------------------------------
!     frac_precip : fraction of grid box undergoing precipitation (unitless)
!     rate_precip : rate of precipitation for grid box (i,j) (mm/day)
!     windsp2     : surface wind speed squared ((m/s)^2)
!     xsoil_nox   : soil NOx (molec/cm^2/s)
!     ----------------------------------------------------------------------

      real*8  :: frac_precip(i1:i2, ju1:j2)
      real*8  :: rate_precip(i1:i2, ju1:j2)
      real*8  :: windsp2    (i1:i2, ju1:j2)
      real*8  :: xsoil_nox  (i1:i2, ju1:j2)

!     ---------------------------------------------
!     canopy_nox : deposition rate constant for NOx
!     ---------------------------------------------

      real*8  :: canopy_nox (i1:i2, ju1:j2, NTYPE)


!     ----------------
!     Begin execution.
!     ----------------

      if (pr_diag) then
        Write (6,*) 'Do_Soil_Emiss called by ', loc_proc
      end if


      frac_precip (:,:) = 0.0d0
      rate_precip (:,:) = 0.0d0
      windsp2     (:,:) = 0.0d0
      xsoil_nox   (:,:) = 0.0d0

      canopy_nox(:,:,:) = 0.0d0


      do ij = ju1, j2
        do il = i1, i2

!         --------------------------------------------------------------
!         Calculate the wind speed (m/s) (see Seinfeld, page 494, 1986).
!         --------------------------------------------------------------

          windsp = (ustar(il,ij) / 0.4d0) *  &
     &             Log (10.0d0 / surf_rough(il,ij))

          windsp2(il,ij) = windsp * windsp

        end do
      end do


!     ====================
      call Calc_Canopy_Nox  &
!     ====================
     &  (ino_num, lwi_flags, ireg, iland, iuse, mw, fracCloudCover,  &
     &   radswg, cossza, tempk, xlai, canopy_nox, &
     &   pr_diag, loc_proc, i1, i2, ju1, j2, num_species)


!     ================
      call Precip_Frac  &
!     ================
     &  (tot_precip, con_precip, frac_precip, rate_precip, &
     &  pr_diag, loc_proc, i1, i2, ju1, j2, j2_gl, latdeg)


      call Soil_Nox ( nymd, ireg, iland, iuse, ncon_soil, tdt, &
     &   frac_precip, rate_precip, radswg, tempk, windsp2, canopy_nox,  &
     &   xlai, soil_fert, soil_precip, soil_pulse, xsoil_nox, idaySoilType, &
     &   pr_diag, loc_proc, i1, i2, ju1, j2, exp_fac)

      cmpm2 = CMPM * CMPM

      do ij = ju1, j2
        do il = i1, i2

          rfac = (mw(ino_num) * KGPG) * (mcor(il,ij) * cmpm2) / AVOGAD

          emiss_nox(il,ij) = xsoil_nox(il,ij) * rfac

        end do
      end do


      return

      end


!-----------------------------------------------------------------------------
!
! ROUTINE
!   Add_Emiss_Harvard
!
! DESCRIPTION
!   This routine adds the Harvard emissions to const:
!
!     1) Take emissions of kg/s and multiply by the time step to get total
!        kg of emission over time step.
!
!     2) Divide by mass of the zone to obtain emissions in terms of mixing
!        ratio; also multiply by the ratio of molecular weight of air to
!        molecular weight of the chemical emission to get volume mixing
!        ratio from mass mixing ratio.
!
!     3) Add emitted mixing ratio amount to existing mixing ratio of const.
!
! ARGUMENTS
!   pr_surf_emiss : should surface emissions be accumulated for output?
!   iisoprene_num : index of isoprene in const
!   ico_num       : index of CO       in const
!   ipropene_num  : index of propene  in const
!   ino_num       : index of NO       in const
!   tdt           : model time step    (s)
!   mw            : array of species' molecular weights (g/mol)
!   mcor          : surface area of grid box (m^2)
!   surf_emiss_out: accumulated surface emissions (kg/m^2/delta t output)
!   mass          : total mass of the atmosphere within each grid box (kg)
!   const         : species concentration, known at zone centers
!                   (mixing ratio)
!   emiss_isop    : isoprene   emissions (kg/s)
!   emiss_monot   : monterpene emissions (kg/s)
!   emiss_nox     : NOx        emissions (kg/s)
!
!-----------------------------------------------------------------------------

      subroutine Add_Emiss_Harvard  &
     &  (pr_surf_emiss, pr_emiss_3d, iisoprene_num, ico_num, ipropene_num,  &
     &   ino_num, tdt, mw, mcor, surf_emiss_out, surf_emiss_out2,  &
     &   emiss_3d_out, mass, concentration, emiss_isop, emiss_monot, emiss_nox, &
     &   do_ShipEmission, emiss_hno3, emiss_o3, ihno3_num, io3_num, &
     &   pr_diag, loc_proc, i1, i2, ju1, j2, k1, k2, num_species)

      use GmiArrayBundlePointer_mod, only : t_GmiArrayBundle

      implicit none

#     include "gmi_phys_constants.h"
#     include "gmi_emiss_constants.h"

!     ----------------------
!     Argument declarations.
!     ----------------------

      logical, intent(in) :: pr_diag
      integer, intent(in) :: loc_proc, i1, i2, ju1, j2, k1, k2, num_species
      logical :: pr_surf_emiss, do_ShipEmission
      logical :: pr_emiss_3d
      integer :: iisoprene_num, ico_num, ipropene_num, ino_num, io3_num, ihno3_num
      real*8  :: tdt
      real*8  :: mw   (num_species)
      real*8  :: mcor           (i1:i2, ju1:j2)
      real*8  :: surf_emiss_out (i1:i2, ju1:j2, num_species)
      real*8  :: surf_emiss_out2(i1:i2, ju1:j2, 6)
      real*8  :: emiss_3d_out   (i1:i2, ju1:j2, k1:k2, num_species)
      real*8  :: mass           (i1:i2, ju1:j2, k1:k2)
      real*8  :: emiss_isop     (i1:i2, ju1:j2)
      real*8  :: emiss_monot    (i1:i2, ju1:j2)
      real*8  :: emiss_nox      (i1:i2, ju1:j2)
      real*8  :: emiss_o3       (i1:i2, ju1:j2)
      real*8  :: emiss_hno3     (i1:i2, ju1:j2)
      type (t_GmiArrayBundle), intent(inout) :: concentration(num_species)

      integer :: ik

!     ----------------------
!     Variable declarations.
!     ----------------------

      real*8  :: emass(i1:i2, ju1:j2)

!     ----------------
!     Begin execution.
!     ----------------

      if (pr_diag) then
        Write (6,*) 'Add_Emiss_Harvard called by ', loc_proc
      end if

      if (iisoprene_num > 0) then

!       ----------------------------
!       Biogenic source of isoprene.
!       ----------------------------

        emass(:,:) =  &
     &    emiss_isop(:,:) * tdt

        concentration(iisoprene_num)%pArray3D(:,:,1) =  &
     &    concentration(iisoprene_num)%pArray3D(:,:,1) +  &
     &    ((emass(:,:) / mass(:,:,1)) *  &
     &     (MWTAIR / mw(iisoprene_num)))

       if (pr_surf_emiss)  &
     &   surf_emiss_out(:,:,iisoprene_num) =  &
     &     surf_emiss_out(:,:,iisoprene_num) +  &
     &     (emass(:,:) / mcor(:,:))


       if (pr_emiss_3d) then
          emiss_3d_out(:,:,1,iisoprene_num) =  &
     &      emiss_3d_out(:,:,1,iisoprene_num) +  &
     &       (emass(:,:) / mcor(:,:))

       end if
      end if

      if (ico_num > 0) then

!       --------------------------------------------------------------------
!       Biogenic source of CO from methanol oxidation, scaled from isoprene.
!       --------------------------------------------------------------------

        emass(:,:) =  &
     &    emiss_isop(:,:) * tdt *  &
     &    ICO_FAC_ISOP

        concentration(ico_num)%pArray3D(:,:,1) =  &
     &    concentration(ico_num)%pArray3D(:,:,1) +  &
     &    ((emass(:,:) / mass(:,:,1)) *  &
     &     (MWTAIR / mw(ico_num)))

       if (pr_surf_emiss)  &
     &   surf_emiss_out(:,:,ico_num) =  &
     &     surf_emiss_out(:,:,ico_num) +  &
     &     (emass(:,:) / mcor(:,:))
!!!!!!!!!!!!!!!!CO_methanol
       if (pr_surf_emiss)  &
     &   surf_emiss_out2(:,:,1) =  &
     &     surf_emiss_out2(:,:,1) +  &
     &     (emass(:,:) / mcor(:,:))
!!!!!!!!!!!!!!!!


       if (pr_emiss_3d) then
          emiss_3d_out(:,:,1,ico_num) =  &
     &      emiss_3d_out(:,:,1,ico_num) +  &
     &       (emass(:,:) / mcor(:,:))

       end if

!       ---------------------------------------------------------------------
!       Biogenic source of CO from monoterpene oxidation.
!       ---------------------------------------------------------------------

        emass(:,:) =  &
     &    emiss_monot(:,:) * tdt *  &
     &    ICO_FAC_MONOT

        concentration(ico_num)%pArray3D(:,:,1) =  &
     &    concentration(ico_num)%pArray3D(:,:,1) +  &
     &    ((emass(:,:) / mass(:,:,1)) *  &
     &     (MWTAIR / mw(ico_num)))

       if (pr_surf_emiss)  &
     &   surf_emiss_out(:,:,ico_num) =  &
     &     surf_emiss_out(:,:,ico_num) +  &
     &     (emass(:,:) / mcor(:,:))

!!!!!!!!!!!!!!!!CO_monoterpene
       if (pr_surf_emiss)  &
     &   surf_emiss_out2(:,:,2) =  &
     &     surf_emiss_out2(:,:,2) +  &
     &     (emass(:,:) / mcor(:,:))
!!!!!!!!!!!!!!!!

       if (pr_emiss_3d) then
          emiss_3d_out(:,:,1,ico_num) =  &
     &      emiss_3d_out(:,:,1,ico_num) +  &
     &       (emass(:,:) / mcor(:,:))
       end if

      end if

!     ---------------------
      if (ipropene_num > 0) then
!     ---------------------

!       -------------------------------------------------
!       Biogenic source of propene, scaled from isoprene.
!       -------------------------------------------------

        emass(:,:) =  &
     &    emiss_isop(:,:) * tdt *  &
     &    BIOSCAL *  &
     &    ((ATOMSC_PER_MOLECISOP / ATOMSC_PER_MOLECPRPE) *  &
     &     (mw(ipropene_num)     / mw(iisoprene_num)))

        concentration(ipropene_num)%pArray3D(:,:,1) =  &
     &    concentration(ipropene_num)%pArray3D(:,:,1) +  &
     &    ((emass(:,:) / mass(:,:,1)) *  &
     &     (MWTAIR / mw(ipropene_num)))

       if (pr_surf_emiss)  &
     &   surf_emiss_out(:,:,ipropene_num) =  &
     &     surf_emiss_out(:,:,ipropene_num) +  &
     &     (emass(:,:) / mcor(:,:))

!!!!!!!!!!!!!!!!Biogenic_propene
       if (pr_surf_emiss)  &
     &   surf_emiss_out2(:,:,3) =  &
     &     surf_emiss_out2(:,:,3) +  &
     &     (emass(:,:) / mcor(:,:))
!!!!!!!!!!!!!!!!

       if (pr_emiss_3d) then
          emiss_3d_out(:,:,1,ipropene_num) =  &
     &      emiss_3d_out(:,:,1,ipropene_num) +  &
     &       (emass(:,:) / mcor(:,:))
       end if
!!!!!!!!!!!!!!!!!!

      end if

!     ----------------
      if (ino_num > 0) then
!     ----------------

!       -------------------
!       Soil source of NOx.
!       -------------------

        emass(:,:) =  &
     &    emiss_nox(:,:) * tdt

        concentration(ino_num)%pArray3D(:,:,1) =  &
     &    concentration(ino_num)%pArray3D(:,:,1) +  &
     &    ((emass(:,:) / mass(:,:,1)) *  &
     &     (MWTAIR / mw(ino_num)))

       if (pr_surf_emiss)  &
     &   surf_emiss_out(:,:,ino_num) =  &
     &     surf_emiss_out(:,:,ino_num) +  &
     &     (emass(:,:) / mcor(:,:))

!!!!!!!!!!!!!!!!Soil_NOx
       if (pr_surf_emiss)  &
     &   surf_emiss_out2(:,:,4) =  &
     &     surf_emiss_out2(:,:,4) +  &
     &     (emass(:,:) / mcor(:,:))
!!!!!!!!!!!!!!!!
       if (pr_emiss_3d) then
          emiss_3d_out(:,:,1,ino_num) =  &
     &      emiss_3d_out(:,:,1,ino_num) +  &
     &       (emass(:,:) / mcor(:,:))

       end if
!!!!!!!!!!!!!!!!!!

      end if

      if (do_ShipEmission) then
         if (ihno3_num > 0) then

            ! --------------
            ! Source of HNO3.
            ! --------------
            emass(:,:) =  emiss_hno3(:,:) * tdt
            if (pr_surf_emiss)  &
     &         surf_emiss_out2(:,:,5) = surf_emiss_out2(:,:,5) +  &
     &                                  (emass(:,:) / mcor(:,:))
         end if

         if (io3_num > 0) then

            ! -------------
            ! Source of O3.
            ! -------------
            emass(:,:) =  emiss_o3(:,:) * tdt
            if (pr_surf_emiss)  &
     &         surf_emiss_out2(:,:,6) = surf_emiss_out2(:,:,6) +  &
     &                                  (emass(:,:) / mcor(:,:))
         end if
      end if

      return

      end


!-----------------------------------------------------------------------------
!
! ROUTINE
!   Add_Semiss_Harvard_Inchem
!
! DESCRIPTION
!   This routine adds in the Harvard emissions, and is used just prior to
!   calling the chemistry solver.
!
! ARGUMENTS
!   iisoprene_num : index of isoprene in const
!   ico_num       : index of CO       in const
!   ipropene_num  : index of propene  in const
!   ino_num       : index of NO       in const
!   mw            : array of species' molecular weights (g/mol)
!   emiss_isop    : isoprene    emissions (kg/s)
!   emiss_monot   : monoterpene emissions (kg/s)
!   emiss_nox     : NOx         emissions (kg/s)
!c?
!   surf_emiss    : surface     emissions (units?)
!   emiss_3d      : 3d          emissions
!
!-----------------------------------------------------------------------------

      subroutine Add_Semiss_Harvard_Inchem  &
     &  (iisoprene_num, ico_num, ipropene_num, ino_num, mw,  &
     &   emiss_isop, emiss_monot, emiss_nox, &
     &   do_ShipEmission, emiss_hno3, emiss_o3, ihno3_num, io3_num, &
     &   conv_emiss, surf_emiss, surf_emiss2, emiss_3d, &
     &   pr_diag, loc_proc, i1, i2, ju1, j2, k1, k2, num_species)

      implicit none

#     include "gmi_emiss_constants.h"


!     ----------------------
!     Argument declarations.
!     ----------------------

      logical, intent(in) :: pr_diag, do_ShipEmission
      integer, intent(in) :: loc_proc, i1, i2, ju1, j2, k1, k2, num_species
      integer :: il, ij, ik
      integer :: iisoprene_num, ico_num, ipropene_num, ino_num, io3_num, ihno3_num
      real*8  :: mw(num_species)
      real*8  :: emiss_isop (i1:i2, ju1:j2)
      real*8  :: emiss_monot(i1:i2, ju1:j2)
      real*8  :: emiss_nox  (i1:i2, ju1:j2)
      real*8  :: emiss_hno3 (i1:i2, ju1:j2)
      real*8  :: emiss_o3(i1:i2, ju1:j2)
      real*8  :: conv_emiss (i1:i2, ju1:j2)
      real*8  :: surf_emiss (i1:i2, ju1:j2, num_species)
      real*8  :: surf_emiss2 (i1:i2, ju1:j2, 6)
      real*8  :: emiss_3d (i1:i2, ju1:j2, k1:k2, num_species)


!     ----------------
!     Begin execution.
!     ----------------

      if (pr_diag) then
        Write (6,*) 'Add_Semiss_Harvard_Inchem called by ', loc_proc
      end if


!     ----------------------
      if (iisoprene_num > 0) then
!     ----------------------

!       ----------------------------
!       Biogenic source of isoprene.
!       ----------------------------

        surf_emiss(:,:,iisoprene_num) =  &
     &    surf_emiss(:,:,iisoprene_num) +  &
     &    (emiss_isop(:,:) * conv_emiss(:,:) / mw(iisoprene_num))

              emiss_3d(:,:,1,iisoprene_num) =  &
     &           emiss_3d(:,:,1,iisoprene_num) +  &
     &           emiss_isop(:,:) * conv_emiss(:,:) / mw(iisoprene_num)

      end if

!     ----------------
      if (ico_num > 0) then
!     ----------------

!       ----------------------------------------------
!       Biogenic source of CO from methanol oxidation,
!       scaled from isoprene.
!       ----------------------------------------------

        surf_emiss(:,:,ico_num) =  &
     &    surf_emiss(:,:,ico_num) +  &
     &    (emiss_isop(:,:) * conv_emiss(:,:) / mw(ico_num)) *  &
     &    ICO_FAC_ISOP

        surf_emiss2(:,:,1) =  &
     &    surf_emiss2(:,:,1) +  &
     &    (emiss_isop(:,:) * conv_emiss(:,:) / mw(ico_num)) *  &
     &    ICO_FAC_ISOP

              emiss_3d(:,:,1,ico_num) =  &
     &           emiss_3d(:,:,1,ico_num) +  &
     &           emiss_isop(:,:) * conv_emiss(:,:) / mw(ico_num) *  &
     &           ICO_FAC_ISOP

!       -------------------------------------------------
!       Biogenic source of CO from monoterpene oxidation.
!       -------------------------------------------------

        surf_emiss(:,:,ico_num) =  &
     &    surf_emiss(:,:,ico_num) +  &
     &    (emiss_monot(:,:) *  conv_emiss(:,:) / mw(ico_num)) *  &
     &    ICO_FAC_MONOT

        surf_emiss2(:,:,2) =  &
     &    surf_emiss2(:,:,2) +  &
     &    (emiss_monot(:,:) *  conv_emiss(:,:) / mw(ico_num)) *  &
     &    ICO_FAC_MONOT

              emiss_3d(:,:,1,ico_num) =  &
     &           emiss_3d(:,:,1,ico_num) +  &
     &           emiss_monot(:,:) * conv_emiss(:,:) / mw(ico_num) *  &
     &           ICO_FAC_MONOT

      end if


!     ---------------------
      if (ipropene_num > 0) then
!     ---------------------

!       -------------------------------------------------
!       Biogenic source of propene, scaled from isoprene.
!       -------------------------------------------------

        surf_emiss(:,:,ipropene_num)  =  &
     &    surf_emiss(:,:,ipropene_num)  +  &
     &    (emiss_isop(:,:) * conv_emiss(:,:) / mw(iisoprene_num)) *  &
     &    BIOSCAL *  &
     &    (ATOMSC_PER_MOLECISOP / ATOMSC_PER_MOLECPRPE)

        surf_emiss2(:,:,3)  =  &
     &    surf_emiss2(:,:,3)  +  &
     &    (emiss_isop(:,:) * conv_emiss(:,:) / mw(iisoprene_num)) *  &
     &    BIOSCAL *  &
     &    (ATOMSC_PER_MOLECISOP / ATOMSC_PER_MOLECPRPE)

              emiss_3d(:,:,1,ipropene_num) =  &
     &           emiss_3d(:,:,1,ipropene_num) +  &
     &           emiss_isop(:,:) * conv_emiss(:,:) / mw(iisoprene_num) *  &
     &           BIOSCAL *  &
     &           (ATOMSC_PER_MOLECISOP / ATOMSC_PER_MOLECPRPE)

      endif
!     ----------------
      if (ino_num > 0) then
!     ----------------

!       -------------------
!       Soil source of NOx.
!       -------------------

        surf_emiss(:,:,ino_num) =  &
     &    surf_emiss(:,:,ino_num) +  &
     &    (emiss_nox(:,:) * conv_emiss(:,:) / mw(ino_num))

        surf_emiss2(:,:,4) =  &
     &    surf_emiss2(:,:,4) +  &
     &    (emiss_nox(:,:) * conv_emiss(:,:) / mw(ino_num))

              emiss_3d(:,:,1,ino_num) =  &
     &           emiss_3d(:,:,1,ino_num) +  &
     &           emiss_nox(:,:) * conv_emiss(:,:) / mw(ino_num)

      end if

      if (do_ShipEmission) then

         if (ihno3_num > 0) then
            ! --------------
            ! Source of HNO3
            ! --------------
            surf_emiss2(:,:,5) =  surf_emiss2(:,:,5) +  &
     &          (emiss_hno3(:,:) * conv_emiss(:,:) / mw(ihno3_num))
         end if

         if (io3_num > 0) then
            ! ------------
            ! Source of O3
            ! ------------
            surf_emiss2(:,:,6) =  surf_emiss2(:,:,6) +  &
     &          (emiss_o3(:,:) * conv_emiss(:,:) / mw(io3_num))
         end if

      end if

      return

      end

