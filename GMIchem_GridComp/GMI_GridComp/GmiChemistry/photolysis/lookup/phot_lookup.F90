
!=============================================================================
!
! $Id$
!
! CODE DEVELOPER
!   Dan Bergmann & John Tannahill, LLNL
!
! FILE
!   phot_lookup.F
!
! ROUTINES
!   Lookup_Qj
!
! HISTORY
!   * September 23, 2004 - Jules Kouatchou
!     Modifications to accomodate for the combined strat/trop chemical
!     mechanism. New variables are: o1d_a, o1d_b, qj_adjust, o1d_T.
!     New calculations are carried out to:
!        - Interpolate the O(1D) correlation factor in solar zenith angle,
!          column ozone, and pressure (see 1d_a, o1d_b, and o1d_T)
!        - Nudge O3 photolysis to represent higher resolution around the
!          O(1D) falloff region (308-328 nm) (see qj_adjust).
!=============================================================================


!-----------------------------------------------------------------------------
!
! ROUTINE
!   Lookup_Qj
!
! DESCRIPTION
!   This is the LLNL implementation of Randy Kawa's photolysis lookup table.
!   It requires a NetCDF file which has the radiative source as a function
!   of wavelength, pressure, column ozone, and solar zenith angle.  This
!   NetCDF file also contains cross sections as a function of species,
!   temperature, and wavelength.
!
! ARGUMENTS
!   do_clear_sky  : do clear sky photolysis?
!   phot_opt      : photolysis option
!   io3_num       : const array index for ozone
!   nymd          : current year/month/day (YYYYMMDD)
!   photintv      : photolysis time step  (s)
!   rsec_jan1     : seconds from Jan. 1st (s)
!   press3c       : atmospheric pressure at the center of each grid box (mb)
!   kel           : temperature (degK)
!   const         : species concentration, known at zone centers (mixing ratio)
!   latdeg        : latitude (deg)
!   mcor          : area of grid box (m^2)
!   mass          : total mass of the atmosphere within each grid box   (kg)
!   fracCloudCover: fractional cloud cover
!   qjgmi         : photolysis rate constants (s^-1)
!
!-----------------------------------------------------------------------------

      subroutine Lookup_Qj (chem_mecha, do_clear_sky, phot_opt, io3_num, nymd, &
                     photintv, rsec_jan1, press3c, kel, concentration,         &
                     solarZenithAngle, mcor, mass, fracCloudCover, qjgmi,      &
                     pr_diag, loc_proc, num_species, num_qjs, num_qjo, &
                     ilo, ihi, julo, jhi, i1, i2, ju1, j2, k1, k2)

      use GmiInterpolation_mod     , only : Interp
      use GmiArrayBundlePointer_mod, only : t_GmiArrayBundle
      use solar_cycle_mod
      use GmiTimeControl_mod  , only : GmiSplitDateTime
      use GmiPrintError_mod, only : GmiPrintError

      implicit none

#     include "gmi_time_constants.h"
#     include "gmi_phys_constants.h"
#     include "phot_lookup_constants.h"
#     include "phot_lookup.h"
#     include "phot_lookup_arrays.h"

!     ----------------------
!     Argument declarations.
!     ----------------------

      CHARACTER(LEN=*), INTENT(IN) :: chem_mecha
      logical, intent(in) :: pr_diag
      integer, intent(in) :: loc_proc
      integer, intent(in) :: num_species, num_qjs, num_qjo
      integer, intent(in) :: ilo, ihi, julo, jhi
      integer, intent(in) :: i1, i2, ju1, j2, k1, k2
      logical :: do_clear_sky
      integer :: phot_opt
      integer :: io3_num
      integer :: nymd
      real*8  :: photintv
      real*8  :: rsec_jan1
      real*8  :: press3c  (ilo:ihi, julo:jhi, k1:k2)
      real*8  :: kel      (ilo:ihi, julo:jhi, k1:k2)
      real*8  :: solarZenithAngle (i1:i2, ju1:j2)
      real*8  :: mcor     (i1:i2, ju1:j2)
      real*8  :: mass     (i1:i2, ju1:j2, k1:k2)
      real*8 , intent(in   ) :: fracCloudCover(i1:i2, ju1:j2)
      type (t_GmiArrayBundle) :: qjgmi(num_qjo)
      type (t_GmiArrayBundle), intent(in) :: concentration(num_species)


!     ----------------------
!     Variable declarations.
!     ----------------------

      character (len=75) :: err_msg

      logical, save :: first = .true.

      integer :: idumday, idumyear
      integer :: iik, ilam, io3
      integer :: isza
      integer :: il, ij, ik, iq
      integer :: month
      integer :: o3_pl_ind, o3_ph_ind
      integer :: o3_pl_indp1, o3_ph_indp1
      integer :: p_ind, p_ind_old, p_indm1
      integer :: sza_ind, t_index
!... variable containing index corresponding to current model year/month,
!...  repeats last 11 years after 2001
      integer :: iscyr
      integer :: inyr, inmon

      real*8  :: days
      real*8  :: decl
      real*8  :: half_the_box
      real*8  :: loc_col_o3
      real*8  :: o3_ph_frac, o3_pl_frac
      real*8  :: p_frac
      real*8  :: rdistsq, ri2_gl, ril
      real*8  :: stime_loc, time
      real*8  :: sza, sza_frac
      real*8  :: total

      real*8,  save :: cos94

      real*8  :: cossza(1)
      real*8  :: loc_col_o3ary(1)
      real*8  :: minary(1)

      real*8  :: rsf(NUMLAM)

      real*8,  save :: lon_dummy(1) = 0.0d0

      real*8  :: o1d_a, o1d_b, qj_adjust
      real*8  :: o1d_T(5)

!     ----------------
!     Begin execution.
!     ----------------

      if (pr_diag) then
        Write (6,*) 'Lookup_Qj called by ', loc_proc
      end if


!     ==========
      if (first) then
!     ==========

        first = .false.

        cos94 = Cos (RADPDEG * 94.0d0)

      end if

!... figure out index for solar cycle array from year and month
      inyr = int(nymd/10000)
      inmon = int(nymd/100)-100*inyr
      if(inyr.ge.70) then
        iscyr = 12*(inyr-47)+inmon
      else
        iscyr = 12*(inyr+100-47)+inmon
      endif
!... recycle last 11 years if beyond end of data
      DO WHILE ((iscyr.gt.NUM_SOLAR_CYC_MON))
        iscyr = iscyr-11*12
      enddo
!... recycle previous 11 years if before beginning of data
      DO WHILE ((iscyr.lt.1))
        iscyr = iscyr+11*12
      enddo
      if (pr_diag) then
        print *,'Solar cycle debug: ',iscyr, loc_proc,s_cycle(3,iscyr)
      end if


      time   = rsec_jan1 + (0.5d0 * photintv)

!      ri2_gl = i2_gl

      call GmiSplitDateTime (nymd, idumyear, month, idumday)


!     ======================
      illoop: do il = i1, i2
!     ======================

!        ril = il
!
!        stime_loc = time + (SECPDY * (ril - 1.0d0) / ri2_gl)
!
!        days = stime_loc / SECPDY
!
!!       ============
!        call Solrpos  &
!!       ============
!     &    (days, decl, rdistsq)


!       =======================
        ijloop: do ij = ju1, j2
!       =======================

!         -----------------------------------------------
!         Calculate solar zenith angle and index into the
!         radiative source function table.
!         -----------------------------------------------

!!         ===========
!          call Solrza  &
!!         ===========
!     &      (days, decl, latdeg(ij), lon_dummy, 1, cossza)

          cossza(1) = cos(solarZenithAngle(il,ij))


!         ---------------------------------------------------------
!         If the solar zenith angle exceeds 94 degrees, set all qjs
!         to zero and exit loop.
!         ---------------------------------------------------------

          if (cossza(1) < cos94) then

            do ik = k1, k2
              do iq = 1, num_qjs
                qjgmi(iq)%pArray3D(il,ij,ik) = 0.0d0
              end do
            end do

!           ============
            cycle ijloop
!           ============

          end if


          sza = Acos (cossza(1)) / RADPDEG

          isza_loop: do isza = 2, NUMSZA

            if (sza <= sza_phot(isza)) then

              sza_ind = isza - 1

!             ==============
              exit isza_loop
!             ==============

            end if

            if (isza == NUMSZA) then
              err_msg =  &
     &          'Error finding solar zenith angle index in Lookup_Qj.'
              call GmiPrintError  &
     &          (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
            end if

          end do isza_loop


          sza_frac = (sza_phot(sza_ind+1) - sza) /  &
     &               (sza_phot(sza_ind+1) - sza_phot(sza_ind))


!         ---------------------------------------------------
!         Initialize pressure index into the radiative source
!         function table.
!         ---------------------------------------------------

          p_ind = NUMPRS

          loc_col_o3 = 0.0d0


!         ==========================
          ikloop: do ik = k2, k1, -1
!         ==========================

            if (phot_opt == 5) then

!             -------------------------------------------------
!             Find the column ozone from the climatology in the
!             lookup table.
!             -------------------------------------------------

              minary(1) = Min (press3c(il,ij,ik), o3_clim_prs(1))

              loc_col_o3ary(1) = loc_col_o3

!             ===========
              call Interp  &
!             ===========
     &          (o3_clim_prs(:), o3_clim(:,il,ij,month), NUM_O3CLIM_PRS,  &
     &           minary(:), loc_col_o3ary(:), 1)

              loc_col_o3 = loc_col_o3ary(1)

            else

!             -----------------------
!             Calculate column ozone.
!             -----------------------

              half_the_box =  &
     &          0.5d0 * concentration(io3_num)%pArray3D(il,ij,ik) * mass(il,ij,ik) *  &
     &          AVOGAD / (mcor(il,ij) * MWTAIR * 10.0d0)

              loc_col_o3 = loc_col_o3 + half_the_box

            end if


!           -------------------------------------
!           Find new pressure index in the table.
!           -------------------------------------

            p_ind_old = p_ind

            iik_loop: do iik = p_ind_old, 1, -1

              p_indm1 = iik

              if (prs_phot(iik) >= press3c(il,ij,ik)) then
!               =============
                exit iik_loop
!               =============
              end if

            end do iik_loop


            if (p_indm1 == 1) then

              p_frac = 0.0d0
              p_ind  = p_indm1

            else if (p_indm1 == NUMPRS) then

              p_frac = 1.0d0
              p_ind  = NUMPRS

            else

              p_ind  = p_indm1 + 1
              p_frac = (prs_phot(p_indm1) - press3c(il,ij,ik)) /  &
     &                 (prs_phot(p_indm1) - prs_phot(p_ind))

            end if


!           -----------------------------------------------------------
!           Find new column ozone indices in the table for both
!           pressure levels surrounding the pressure for this grid box.
!           -----------------------------------------------------------

            io3_loop1: do io3 = 2, NUMO3

              if (loc_col_o3 < col_o3(io3, p_ind)) then

                o3_pl_ind   = io3 - 1
                o3_pl_indp1 = io3

                o3_pl_frac =  &
     &            (col_o3(o3_pl_indp1,p_ind) - loc_col_o3) /  &
     &            (col_o3(o3_pl_indp1,p_ind) - col_o3(o3_pl_ind,p_ind))

                o3_pl_frac = Max (o3_pl_frac, 0.0d0)

!               ==============
                exit io3_loop1
!               ==============

              end if

              if (io3 == NUMO3) then
                o3_pl_ind   = NUMO3
                o3_pl_indp1 = NUMO3
                o3_pl_frac  = 1.0d0
              end if

            end do io3_loop1


            io3_loop2: do io3 = 2, NUMO3

              if (loc_col_o3 < col_o3(io3,p_indm1)) then

                o3_ph_ind   = io3 - 1
                o3_ph_indp1 = io3

                o3_ph_frac =  &
     &            (col_o3(o3_ph_indp1,p_indm1) - loc_col_o3) /  &
     &            (col_o3(o3_ph_indp1,p_indm1) -  &
     &             col_o3(o3_ph_ind,p_indm1))

                o3_ph_frac = Max (o3_ph_frac, 0.0d0)

!               ==============
                exit io3_loop2
!               ==============

              end if

              if (io3 == NUMO3) then
                o3_ph_ind   = NUMO3
                o3_ph_indp1 = NUMO3
                o3_ph_frac  = 1.0d0
              end if

            end do io3_loop2


!           --------------------------------------------------
!           Interpolate the radiative source function in solar
!           zenith angle, column ozone, and pressure.
!           --------------------------------------------------

            do ilam = 1, NUMLAM

              rsf(ilam) =  &
     &         sza_frac *  &
     &          (p_frac *  &
     &           ((o3_pl_frac) *  &
     &            rad_source(ilam,sza_ind,o3_pl_ind,p_ind) +  &
     &            (1.0d0 - o3_pl_frac) *  &
     &            rad_source(ilam,sza_ind,o3_pl_indp1,p_ind)) +  &
     &          (1.0d0 - p_frac) *  &
     &           ((o3_ph_frac) *  &
     &            rad_source(ilam,sza_ind,o3_ph_ind,p_indm1) +  &
     &            (1.0d0 - o3_ph_frac) *  &
     &            rad_source(ilam,sza_ind,o3_ph_indp1,p_indm1))) +  &
     &         (1.0d0 - sza_frac) *  &
     &          (p_frac *  &
     &           ((o3_pl_frac) *  &
     &            rad_source(ilam,sza_ind+1,o3_pl_ind,p_ind) +  &
     &            (1.0d0 - o3_pl_frac) *  &
     &            rad_source(ilam,sza_ind+1,o3_pl_indp1,p_ind)) +  &
     &          (1.0d0 - p_frac) *  &
     &           ((o3_ph_frac) *  &
     &            rad_source(ilam,sza_ind+1,o3_ph_ind,p_indm1) +  &
     &            (1.0d0 - o3_ph_frac) *  &
     &            rad_source(ilam,sza_ind+1,o3_ph_indp1,p_indm1)))

            end do

            if (TRIM(chem_mecha) == 'strat_trop' .OR. TRIM(chem_mecha) == 'strat_trop_aerosol') then

!              --------------------------------------------------
!              Interpolate the O(1D) correlation factor in solar
!              zenith angle, column ozone, and pressure.
!              --------------------------------------------------

               do ilam = 1, 5

                 o1d_a =  &
     &            sza_frac *  &
     &             (p_frac *  &
     &              ((o3_pl_frac) *  &
     &               o1d_coef(ilam,sza_ind,o3_pl_ind,p_ind,1) +  &
     &               (1.0d0 - o3_pl_frac) *  &
     &               o1d_coef(ilam,sza_ind,o3_pl_indp1,p_ind,1)) +  &
     &             (1.0d0 - p_frac) *  &
     &              ((o3_ph_frac) *  &
     &               o1d_coef(ilam,sza_ind,o3_ph_ind,p_indm1,1) +  &
     &               (1.0d0 - o3_ph_frac) *  &
     &               o1d_coef(ilam,sza_ind,o3_ph_indp1,p_indm1,1))) +  &
     &            (1.0d0 - sza_frac) *  &
     &             (p_frac *  &
     &              ((o3_pl_frac) *  &
     &               o1d_coef(ilam,sza_ind+1,o3_pl_ind,p_ind,1) +  &
     &               (1.0d0 - o3_pl_frac) *  &
     &               o1d_coef(ilam,sza_ind+1,o3_pl_indp1,p_ind,1)) +  &
     &             (1.0d0 - p_frac) *  &
     &              ((o3_ph_frac) *  &
     &               o1d_coef(ilam,sza_ind+1,o3_ph_ind,p_indm1,1) +  &
     &               (1.0d0 - o3_ph_frac) *  &
     &               o1d_coef(ilam,sza_ind+1,o3_ph_indp1,p_indm1,1)))

                 o1d_b =  &
     &            sza_frac *  &
     &             (p_frac *  &
     &              ((o3_pl_frac) *  &
     &               o1d_coef(ilam,sza_ind,o3_pl_ind,p_ind,2) +  &
     &               (1.0d0 - o3_pl_frac) *  &
     &               o1d_coef(ilam,sza_ind,o3_pl_indp1,p_ind,2)) +  &
     &             (1.0d0 - p_frac) *  &
     &              ((o3_ph_frac) *  &
     &               o1d_coef(ilam,sza_ind,o3_ph_ind,p_indm1,2) +  &
     &               (1.0d0 - o3_ph_frac) *  &
     &               o1d_coef(ilam,sza_ind,o3_ph_indp1,p_indm1,2))) +  &
     &            (1.0d0 - sza_frac) *  &
     &             (p_frac *  &
     &              ((o3_pl_frac) *  &
     &               o1d_coef(ilam,sza_ind+1,o3_pl_ind,p_ind,2) +  &
     &               (1.0d0 - o3_pl_frac) *  &
     &               o1d_coef(ilam,sza_ind+1,o3_pl_indp1,p_ind,2)) +  &
     &             (1.0d0 - p_frac) *  &
     &              ((o3_ph_frac) *  &
     &               o1d_coef(ilam,sza_ind+1,o3_ph_ind,p_indm1,2) +  &
     &               (1.0d0 - o3_ph_frac) *  &
     &               o1d_coef(ilam,sza_ind+1,o3_ph_indp1,p_indm1,2)))

                 o1d_T(ilam) = o1d_a + (kel(il,ij,ik) - 298.0d0) * o1d_b

               end do
            end if

!           -------------------------------------------
!           Interpolate the O2 photolysis rate in solar
!           zenith angle, column ozone, and pressure.
!           -------------------------------------------
!... weighted by solar_cycle in bin 3 of wavelength bins (primarily 185-220 nm)
            if (num_qj_o2 /= 0) then

              qjgmi(num_qj_o2)%pArray3D(il,ij,ik) =  &
     &         sza_frac * s_cycle(3,iscyr) *  &
     &          (p_frac *  &
     &           ((o3_pl_frac) *  &
     &            o2_qj(sza_ind,o3_pl_ind,p_ind) +  &
     &            (1.0d0 - o3_pl_frac) *  &
     &            o2_qj(sza_ind,o3_pl_indp1,p_ind)) +  &
     &          (1.0d0 - p_frac) *  &
     &           ((o3_ph_frac) *  &
     &            o2_qj(sza_ind,o3_ph_ind,p_indm1) +  &
     &            (1.0d0 - o3_ph_frac) *  &
     &            o2_qj(sza_ind,o3_ph_indp1,p_indm1))) +  &
     &         (1.0d0 - sza_frac) *  &
     &          (p_frac *  &
     &           ((o3_pl_frac) *  &
     &            o2_qj(sza_ind+1,o3_pl_ind,p_ind) +  &
     &            (1.0d0 - o3_pl_frac) *  &
     &            o2_qj(sza_ind+1,o3_pl_indp1,p_ind)) +  &
     &          (1.0d0 - p_frac) *  &
     &           ((o3_ph_frac) *  &
     &           o2_qj(sza_ind+1,o3_ph_ind,p_indm1) +  &
     &           (1.0d0 - o3_ph_frac) *  &
     &           o2_qj(sza_ind+1,o3_ph_indp1,p_indm1)))

            end if


!           -------------------------------------------
!           Interpolate the NO photolysis rate in solar
!           zenith angle, column ozone, and pressure.
!           -------------------------------------------
!... weighted by solar_cycle in bin 6 of wavelength bins (primarily 190.9 and
!...  182.7 nm) in the middle atmosphere, best fit to average of bin 5 and 10.

            if (num_qj_no /= 0) then

              qjgmi(num_qj_no)%pArray3D(il,ij,ik) =  &
     &         sza_frac * s_cycle(6,iscyr) *  &
     &          (p_frac *  &
     &           ((o3_pl_frac) *  &
     &            no_qj(sza_ind,o3_pl_ind,p_ind) +  &
     &            (1.0d0 - o3_pl_frac) *  &
     &            no_qj(sza_ind,o3_pl_indp1,p_ind)) +  &
     &          (1.0d0 - p_frac) *  &
     &           ((o3_ph_frac) *  &
     &            no_qj(sza_ind,o3_ph_ind,p_indm1) +  &
     &            (1.0d0 - o3_ph_frac) *  &
     &            no_qj(sza_ind,  o3_ph_indp1,p_indm1))) +  &
     &         (1.0d0 - sza_frac) *  &
     &          (p_frac *  &
     &           ((o3_pl_frac) *  &
     &            no_qj(sza_ind+1,o3_pl_ind,p_ind) +  &
     &            (1.0d0 - o3_pl_frac) *  &
     &            no_qj(sza_ind+1,o3_pl_indp1,p_ind)) +  &
     &          (1.0d0 - p_frac) *  &
     &           ((o3_ph_frac) *  &
     &            no_qj(sza_ind+1,o3_ph_ind,p_indm1) +  &
     &            (1.0d0 - o3_ph_frac) *  &
     &            no_qj(sza_ind+1,o3_ph_indp1,p_indm1)))

            end if


!           ------------------------------------------------------
!           Calculate the temperature index into the cross section
!           data which lists coss sections for temperatures from
!           150 to 349 degrees K.  Make sure the index is a value
!           between 1 and 200.
!           ------------------------------------------------------

            t_index = kel(il,ij,ik) - 148.5d0

            t_index = Min (200, Max (t_index, 1))


!           ------------------------------------------------------------------
!           For each species, integrate over the wavelengths to
!           calculate a qj.
!
!           There are special cases for:
!             oxygen          : photolytic rxn,  O2   + hv = O  + O
!             nitric oxide    : photolytic rxn,  NO   + hv = N  + O
!             formaldehyde    : photolytic rxn,  CH2O + hv = CO + H2
!             acetone         : photolytic rxn,  ACET + hv = MCO3 + MO2
!             methyl glyoxal  : photolytic rxn,  MGLY + hv = MCO3 + HO2 + CO
!             hydroxy acetone : photolytic rxn,  HACN + hv = MCO3 + HCHO + HO2
!           ------------------------------------------------------------------

            do iq = 1, num_qjs

              if ((iq /= num_qj_o2) .and. (iq /= num_qj_no)) then

                total = 0.0d0

                do ilam = 1, NUMLAM

                  total = total +  &
     &                    (rsf(ilam) * s_cycle(ilam,iscyr) *  &
     &                     cross_section(ilam,t_index,iq))

                end do

                qjgmi(iq)%pArray3D(il,ij,ik) = total

              end if

            end do

            if (TRIM(chem_mecha) == 'strat_trop' .OR. TRIM(chem_mecha) == 'strat_trop_aerosol') then
!              --------------------------------------------------
!              Nudge O3 photolysis to represent higher resolution
!              around the O(1D) falloff region (308-328 nm).
!              --------------------------------------------------

               qj_adjust = sum( rsf(49:53) * s_cycle(49:53,iscyr) *  &
     &                       cross_section(49:53,t_index,3) *  &
     &                       (o1d_T(:) - 1.0d0)  &
     &                      ,DIM = 1)

               qjgmi(2)%pArray3D(il,ij,ik) = qjgmi(2)%pArray3D(il,ij,ik) - qj_adjust
               qjgmi(3)%pArray3D(il,ij,ik) = qjgmi(3)%pArray3D(il,ij,ik) + qj_adjust
            end if

!           --------------------------------------------------
!           Patch up the pressure dependent quantum yields for
!           formaldehyde (CH2O) photolytic reaction.
!           --------------------------------------------------

            if (num_qj_ch2o /= 0) then

              do ilam = 55, 59

                qjgmi(num_qj_ch2o)%pArray3D(il,ij,ik) =  &
     &            qjgmi(num_qj_ch2o)%pArray3D(il,ij,ik) -  &
     &            ((rsf(ilam) * s_cycle(ilam,iscyr) *  &
     &              cross_section(ilam,t_index,num_qj_ch2o)) *  &
     &             (1.0d0 - CH2O_CONST_1(ilam) /  &
     &              (1.0d0 + (press3c(il,ij,ik) / 1013.25d0) *  &
     &               CH2O_CONST_2(ilam))))

              end do

            end if


!           -------------------------------------------------------------
!           Patch up the pressure dependent quantum yields for
!           acetone (ACET) or hydroxy acetone (HACN) photolytic reaction.
!           -------------------------------------------------------------

            if (num_qj_acet /= 0) then

              do ilam = 47, 58

                qjgmi(num_qj_acet)%pArray3D(il,ij,ik) =  &
     &            qjgmi(num_qj_acet)%pArray3D(il,ij,ik) -  &
     &            ((rsf(ilam) * s_cycle(ilam,iscyr) *  &
     &              cross_section(ilam,t_index,num_qj_acet)) *  &
     &             (1.0d0 -  &
     &              (1.0d0 /  &
     &               (ACET_CONST_1(ilam) +  &
     &                ACET_CONST_2(ilam) * press3c(il,ij,ik) /  &
     &                (kel(il,ij,ik) * BOLTZMN_E * BPMB)))))

              end do

            end if

            if (num_qj_hacn /= 0) then

              do ilam = 47, 58

                qjgmi(num_qj_hacn)%pArray3D(il,ij,ik) =  &
     &            qjgmi(num_qj_hacn)%pArray3D(il,ij,ik) -  &
     &            ((rsf(ilam) * s_cycle(ilam,iscyr) *  &
     &              cross_section(ilam,t_index,num_qj_hacn)) *  &
     &             (1.0d0 -  &
     &              (1.0d0 /  &
     &               (ACET_CONST_1(ilam) +  &
     &                ACET_CONST_2(ilam) * press3c(il,ij,ik) /  &
     &                (kel(il,ij,ik) * BOLTZMN_E * BPMB)))))

              end do

            end if


!           --------------------------------------------------
!           Patch up the pressure dependent quantum yields for
!           methyl glyoxal (MGLY) photolytic reaction.
!           --------------------------------------------------

            if (num_qj_mgly /= 0) then

              do ilam = 53, 76

                qjgmi(num_qj_mgly)%pArray3D(il,ij,ik) =  &
     &            qjgmi(num_qj_mgly)%pArray3D(il,ij,ik) -  &
     &            ((rsf(ilam) * s_cycle(ilam,iscyr) *  &
     &              cross_section(ilam,t_index,num_qj_mgly)) *  &
     &             (1.0d0 -  &
     &              (1.0d0 /  &
     &               (MGLY_CONST_1(ilam) +  &
     &                MGLY_CONST_2(ilam) * press3c(il,ij,ik) /  &
     &                (kel(il,ij,ik) * BOLTZMN_E * BPMB)))))

              end do

            end if


!           ------------------------------------------------
!           Add in the bottom half of this box to the column
!           ozone calculation.
!           ------------------------------------------------

            if (phot_opt /= 5) then
              loc_col_o3 = loc_col_o3 + half_the_box
            end if

!         =============
          end do ikloop
!         =============

!       =============
        end do ijloop
!       =============

!     =============
      end do illoop
!     =============


!     -----------------------------------------------------------------------
!     Check to see if cloud fractions should be used to attenuate the
!     photolysis rate constants. If so, calculate a clear sky fraction (csf).
!     -----------------------------------------------------------------------

      if (.not. do_clear_sky) then

        do ik = k1, k2

          do iq = 1, num_qjs

            qjgmi(iq)%pArray3D(:,:,ik) =  &
     &        qjgmi(iq)%pArray3D(:,:,ik) *  &
     &        (0.5d0 * fracCloudCover(:,:) + 0.5d0)

          end do

        end do

      end if


      return

      end

