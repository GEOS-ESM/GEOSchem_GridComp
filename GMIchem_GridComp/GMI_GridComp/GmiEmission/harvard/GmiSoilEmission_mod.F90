!------------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1              !
!------------------------------------------------------------------------------
!BOP
! 
! !MODULE: GmiSoilEmission_mod
!
module GmiSoilEmission_mod
!
! !USES:
      use GmiResistance_mod, only : CanopyResistance, SurfaceResistance
      use GmiTimeControl_mod, only : GmiSplitDateTime
!
      implicit none

#     include "gmi_emiss_constants.h"
#     include "gmi_phys_constants.h"
#     include "gmi_time_constants.h"
!
      private
      public  :: Calc_Canopy_Nox
      public  :: Precip_Frac
      public  :: Soil_Nox
!
! !DECLARED PARAMETERS:
      real*8, parameter :: UNITCONV = 4.3d9  ! ng N/m^2/s -> molec/cm^2/s
!
! !DESCRIPTION:
!  Harvard tropospheric emissions module for 3D applications;
!  by Yuhang Wang, Gerry Gardner, and Prof. Daniel Jacob
!  of Harvard University (Release V2.1)
!
! !HISTORY:
! 8Feb2011 Jules Kouatchou
!     Rewrote Soil_Nox to do all calculations in 2D.
!
! !AUTHOR:
! John Tannahill, LLNL, jrt@llnl.gov
! Jules Kouatchou, NASA/GSFC, Jules.Kouatchou-1@nasa.gov
!
!EOP
!------------------------------------------------------------------------------
      CONTAINS
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Calc_Canopy_Nox
!
! !INTERFACE:
!
      subroutine Calc_Canopy_Nox (ino_num, lsnow, ireg, iland, iuse, mw,       &
                     cfrac, radiat, suncos, tempk, xlai, canopynox, pr_diag,   &
                     loc_proc, i1, i2, ju1, j2, num_species)
!
#     include "gmi_drydep_data.h"
!
! !INPUT PARAMETERS:
      logical, intent(in) :: pr_diag
      integer, intent(in) :: loc_proc
      integer, intent(in) :: i1, i2, ju1, j2
      integer, intent(in) :: num_species
      integer, intent(in) :: ino_num                      ! index of NO in const
      integer, intent(in) :: lsnow (i1:i2, ju1:j2)        ! array of flags that indicate land, water, or ice
                                                          ! (1 => water, 2 => land, 3 => ice)
      integer, intent(in) :: ireg  (i1:i2, ju1:j2)        ! number of land types in a grid square
      integer, intent(in) :: iland (i1:i2, ju1:j2, NTYPE) ! land type id in grid square for ireg land types
      integer, intent(in) :: iuse  (i1:i2, ju1:j2, NTYPE) ! fraction of grid box area occupied by land type (mil^-1?)
      real*8 , intent(in) :: mw    (num_species)          ! species' molecular weights (g/mol)
      real*8 , intent(in) :: cfrac (i1:i2, ju1:j2)        ! fractional cloud cover
      real*8 , intent(in) :: radiat(i1:i2, ju1:j2)        ! solar radiation (W/m^2)
      real*8 , intent(in) :: suncos(i1:i2, ju1:j2)        ! cosines of the solar zenith angle
      real*8 , intent(in) :: tempk (i1:i2, ju1:j2)        ! surface air temperature (degK)
      real*8 , intent(in) :: xlai  (i1:i2, ju1:j2, NTYPE) ! leaf area index of land type for current month
!
! !OUTPUT PARAMETERS:
      real*8, intent(out) :: canopynox(i1:i2, ju1:j2, NTYPE) ! deposition rate constant for NOx
!
! !DECLARED PARAMETERS:
      real*8, parameter :: F0_NO    = 0.0d0  ! reactivity factor for oxidation 
                                             ! of biological substances
      real*8, parameter :: HSTAR_NO = 1.9d-3 ! Henry's Law constant for NO
      real*8, parameter :: PRESS    = 1.5d5
!
! !DESCRIPTION:
!  Calculates canopynox using a resistance-in-series model
!   (Harvard Model / Depvel Version 3.2:  5/27/97).
!
! !LOCAL VARIABLES:
      integer :: idep1, il, ij, iolson, ldt
      real*8  :: rdc, rix, rsurfc, rt, tempc1, tempk1
      real*8  :: rac (NTYPE), rclo(NTYPE), rcls(NTYPE)
      real*8  :: rgso(NTYPE), rgss(NTYPE), ri  (NTYPE)
      real*8  :: rlu (NTYPE)
!EOP
!------------------------------------------------------------------------------
!BOC
      if (pr_diag) Write (6,*) 'Calc_Canopy_Nox called by ', loc_proc

      do ij = ju1, j2
         do il = i1, i2

            tempk1 = tempk(il,ij)
            tempc1 = tempk1 + ABS_ZERO

            ! --------------------------------------------------------------------
            ! Compute bulk surface resistance for gases.  Adjust external surface
            ! resistances for temperature; from Wesely [1989], expression given in
            ! text on p. 1296.  There is no evidence that the resistance continues
            ! to increase at temperatures below -18 C, so at colder temperatures
            ! hold the resistance fixed.
            ! --------------------------------------------------------------------

            if (tempc1 < -18.0d0) then
               rt = 1.2d9
            else
               rt = 1000.0d0 * Exp (-tempc1 - 4.0d0)
            end if

            ! --------------------------------------------------------------------
            ! Get surface resistances; loop over land types, ldt.  The land types
            ! within each grid square are defined using the Olson land-type
            ! database.  Each of the Olson land types is assigned a corresponding
            ! "deposition land type" with characteristic values of surface
            ! resistance components.  There are 74 Olson land-types but only 11
            ! deposition land-types (i.e., many of the Olson land types share the
            ! same deposition characteristics).  Surface resistance components for
            ! the "deposition land types" are from Wesely [1989] except for
            ! tropical forests [Jacob and Wofsy, 1990] and for tundra [Jacob,
            ! et al., 1992].  All surface resistance components are normalized to
            ! a leaf area index of unity.  Olson land types, deposition land
            ! types, and surface resistance components are read from an input data
            ! file; check that file for further details.
            ! --------------------------------------------------------------------

            !================================
            LDTLOOP: do ldt = 1, ireg(il,ij)
            !================================
                                        !=============
               if (iuse(il,ij,ldt) == 0) cycle LDTLOOP
                                        !=============

               if (lsnow(il,ij) == 3) then  ! snow or ice
                  idep1  = 1
               else
                  iolson = iland(il,ij,ldt) + 1
                  idep1  = idep(iolson)
               end if

               call SurfaceResistance (idep1, rac(ldt), rclo(ldt), rcls(ldt),   &
                           rgso(ldt), rgss(ldt), ri(ldt), rlu(ldt), rt, tempc1, &
                           cfrac(il,ij), radiat(il,ij), suncos(il,ij),          &
                           xlai(il,ij,ldt), rix)

               !----------------------------------------------------------------
               ! Compute aerodynamic resistance to lower elements in lower part
               ! of the canopy or structure, assuming level terrain; equation (5)
               ! of Wesely [1989].
               !----------------------------------------------------------------

               rdc = 100.0d0 * (1.0d0 + (1000.0d0 / (radiat(il,ij) + 10.0d0)))

               !---------------------------------------------------------------
               ! Species-dependent corrections to resistances are from equations
               ! (6)-(9) of Wesely [1989].
               !---------------------------------------------------------------

               call CanopyResistance (rdc, rix, PRESS, tempk1, F0_NO, HSTAR_NO,   &
                          mw(ino_num), rac(ldt), rclo(ldt), rcls(ldt), rgso(ldt), &
                          rgss(ldt), rlu(ldt), rsurfc)

               canopynox(il,ij,ldt) =  1.0d0 / rsurfc

            !==============
            end do LDTLOOP
            !==============
         end do
      end do

      return

      end subroutine Calc_Canopy_Nox
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Precip_Frac
!
! !INTERFACE:
!
      subroutine Precip_Frac (preacc, precon, frac, rate, &
                        pr_diag, loc_proc, i1, i2, ju1, j2, j2_gl, latdeg)
!
! !INPUT PARAMETERS:
      logical, intent(in) :: pr_diag
      integer, intent(in) :: loc_proc, i1, i2, ju1, j2, j2_gl
               ! DAO total      precipitation at ground (mm/day)
      real*8 , intent(in) :: preacc(i1:i2, ju1:j2)
               ! DAO convective precipitation at ground (mm/day)
      real*8 , intent(in) :: precon(i1:i2, ju1:j2)
      real*8 , intent(in) :: latdeg(i1:i2, ju1:j2)
!
! !OUTPUT PARAMETERS:
               ! fraction of grid box undergoing precipitation (unitless)
      real*8 , intent(out) :: frac  (i1:i2, ju1:j2)
               ! rate of precipitation for grid box (i,j)      (mm/day)
      real*8 , intent(out) :: rate  (i1:i2, ju1:j2)
!
! !DESCRIPTION:
!   This routine computes the fraction of a grid box that is actually
!   precipitating, along with the precipitation rate.
!
!   This version of Precip_Frac replaces Yuhang Wang's original version, as
!   used in the Harvard code prior to 10/18/99.
!
!  Reference:
!    Liu, H. Y., D. J. Jacob, I. Bey, R. M. Yantosca, and D. M. Koch,
!    Three-dimensional simulation of $210Pb$ and $7Be$ in the Harvard-DAO
!    tropospheric chemistry model, Eos Trans. AGU, 80 (17), S32, 1999a.
!
! !LOCAL VARIABLES:
      integer :: il, ij
      integer :: jstart, jend
      real*8  :: frac_convec, frac_lsprecip
!
!EOP
!------------------------------------------------------------------------------
!BOC
      if (pr_diag) Write (6,*) 'Precip_Frac called by ', loc_proc

      ! ---------------------------------------------------------
      ! For the polar boxes there is no precipitation.
      ! This assumption has been removed for GEOS-5 Ganymed tags.
      ! ---------------------------------------------------------

      jstart = ju1
      jend = j2

      do ij = jstart, jend
         do il = i1, i2

            ! ---------------------------------------------------------------
            ! Large-scale precipitation at (il,ij) = preacc(il,ij) - precon(il,ij).
            ! If there is large-scale precipitation at grid box (il,ij), then
            ! assume that it covers 7% of the area of grid box(il,ij).
            ! ---------------------------------------------------------------

            if ((preacc(il,ij) - precon(il,ij)) > 0.0d0) then
               frac_lsprecip = 7.0d-2
            else
               frac_lsprecip = 0.0d0
            end if

            ! ----------------------------------------------------------------
            ! Convective precipitation at (il,ij) = precon(il,ij); if there is
            ! convective precipitation at (il,ij), then assume that it covers
            ! 0.3% of the area of grid box (il,ij).
            ! ----------------------------------------------------------------

            if (precon(il,ij) > 0.0d0) then
               frac_convec = 3.0d-3
            else
               frac_convec = 0.0d0
            end if

            ! -------------------------------------------------------------------
            ! frac = total fraction of grid box (il,ij) covered by precipitation;
            ! the possible values of frac are:  0.0%, 0.3%, 7.0%, or 7.3%.
            ! -------------------------------------------------------------------

            frac(il,ij) = frac_lsprecip + frac_convec

            ! -----------------------------------------------------------
            ! rate = total precipitation rate in mm/day, adjusted for the
            ! fraction of the grid box that is precipitating.
            ! -----------------------------------------------------------

            if (frac(il,ij) > 0.0d0) then
               rate(il,ij) = preacc(il,ij) / frac(il,ij)
            else
               rate(il,ij) = 0.0d0
            end if

         end do
      end do

      return

      end subroutine Precip_Frac
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Soil_Nox
!
! !INTERFACE:
!
      subroutine Soil_Nox ( nymd, ireg, iland, iuse, nconsoil, tdt, &
                     frac, rate, radiat, tempk, windsqr, canopynox,     &
                     xlai, soilfert, soilprep, soilpuls, xsoilnox, &
                     idaySoilType, pr_diag, loc_proc, i1, i2, ju1, j2, &
                     exp_fac)
!
!INPUT PARAMETERS:
      logical, intent(in) :: pr_diag
      integer, intent(in) :: loc_proc, i1, i2, ju1, j2
      integer, intent(in) :: nymd  ! year/month/day  (YYYYMMDD)
      integer, intent(in) :: ireg(i1:i2, ju1:j2)  ! number of land types in a grid square
      integer, intent(in) :: iland(i1:i2, ju1:j2, NTYPE) ! land type id in grid square for ireg land types
      integer, intent(in) :: iuse(i1:i2, ju1:j2, NTYPE) ! fraction of grid box area occupied by land type (mil^-1?)
      integer, intent(in) :: nconsoil (NVEGTYPE)  ! Olson -> soil type
      real*8, intent(in) :: exp_fac(NPULSE)    ! pulsing decay per time step (day^-1)
      real*8,  intent(in) :: tdt ! model time step (s)
      real*8,  intent(in) :: frac(i1:i2, ju1:j2) ! precipitating fraction of grid square
      real*8,  intent(in) :: rate(i1:i2, ju1:j2) ! local precipitation rate, 
                                                 ! i.e., grid-scale precipitation
                                                 ! rate divided by frac (mm/day)
      real*8,  intent(in) :: radiat(i1:i2, ju1:j2)  ! solar radiation (W/m^2)
      real*8,  intent(in) :: tempk(i1:i2, ju1:j2)   ! temperature     (degK)
      real*8,  intent(in) :: windsqr(i1:i2, ju1:j2) ! surface wind speed squared ((m/s)^2)
      real*8,  intent(in) :: canopynox(i1:i2, ju1:j2, NTYPE) ! deposition rate constant for NOx
      real*8,  intent(in) :: xlai(i1:i2, ju1:j2, NTYPE) ! leaf area index of land type for current month
      real*8,  intent(in) :: soilfert(i1:i2, ju1:j2) ! fertilizers (ng N/m^2/s)
      real*8,  intent(in) :: soilprep(i1:i2, ju1:j2) ! observed precipitation 
                                                     ! (mm/day/box)
!
! !INPUT/OUTPUT PARAMETERS:
      integer, intent(inOut) :: idaySoilType
      real*8, intent(inOut) :: soilpuls(NPULSE+1, i1:i2, ju1:j2) 
               ! tracking of wet/dry & three types of pulsing (Y&L, 94);
      real*8, intent(inOut) :: xsoilnox(i1:i2, ju1:j2) ! soil NOx (molec/cm^2/s)
!
! !DESCRIPTION:
!   This routine calculates the total NOx emission from the soil.
!
!   Based on Yienger and Levy [1995];
!   see Wang et al [1998]: Global Simulation of Tropospheric
!   O3-NOx-hydrocarbon; JGR Vol 103, pages 10713-10725.
!
! \begin{verbatim}
!   soilpuls : tracking of wet/dry & three types of pulsing (Y&L, 94);
!     soilpuls(1,mm)    => flag for wet (-1) or dry (+1) soil; only dry soil
!                          is subject to pulsing
!     soilpuls(nn+1,mm) => fraction of grid square affected by equivalent
!                          fresh pulsing (e.g., sprinkle, shower, heavy rain)
! \end{verbatim}
!
! !LOCAL VARIABLES:
      integer :: iday, imonth
      integer :: ii, ij, it
      integer :: mm, nn
      integer :: ydummy
      real*8  :: riuse
      real*8  :: rpulse
      real*8  :: tempc1

      real*8 :: soil_temp_fac_tmp
      real*8 :: soil_base_tmp
      real*8 :: soil_fert_fac_tmp
      real*8 :: soil_can_red_fac_tmp

!EOP
!------------------------------------------------------------------------------
!BOC
      if (pr_diag) Write (6,*) 'Soil_Nox called by ', loc_proc

!     ====================
      call GmiSplitDateTime (nymd, ydummy, imonth, iday)
!     ====================

      call Soil_Type (iday, idaySoilType, soilprep, soilpuls, pr_diag, loc_proc, &
                      i1, i2, ju1, j2)

      do ij = ju1, j2
         do ii = i1, i2

            ! ---------------------------------------------------
            ! Pulsing factor "function Pulseit ( )";
            ! ECO system dependent;
            ! Temperature factor "function Soil_Temp_Fac ( )";
            ! Base emission with fertilization;
            ! Canopy reduction;
            ! Soil NOx emissions (watch out for trop. evergreen).
            ! ---------------------------------------------------

            rpulse = Pulseit(ii, ij, tdt, frac(ii,ij), rate(ii,ij), soilpuls, &
                             i1, i2, ju1, j2, exp_fac)
  
            tempc1 = tempk(ii,ij) + ABS_ZERO

            do it = 1, ireg(ii,ij)
               nn    = nconsoil(iland(ii,ij,it)+1)
               riuse = iuse(ii,ij,it)

               soil_temp_fac_tmp = Soil_Temp_Fac(nn, tempc1, soilpuls(1,ii,ij))
               soil_base_tmp = Soil_Base(nn, soilprep(ii,ij), &
                                             soilpuls(1,ii,ij), rpulse)
               soil_fert_fac_tmp = Soil_Fert_Fac(nn,soilfert(ii,ij))
               soil_can_red_fac_tmp = Soil_Can_Red_Fac (nn, radiat(ii,ij), &
     &             windsqr(ii,ij), canopynox(ii,ij,it), xlai(ii,ij,it))

               xsoilnox(ii,ij) = xsoilnox(ii,ij) +  &
     &             ((Soil_Temp_Fac_tmp * soil_base_tmp) + Soil_Fert_Fac_tmp) * &
     &             (1.0d0 - soil_can_red_fac_tmp) * riuse / 1000.0d0

            end do
         end do
      end do

      return

      end subroutine Soil_Nox
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Pulseit
!
! !INTERFACE:
!
      function Pulseit (i, j, tdt, frac1, rate1, soilpuls, i1, i2, ju1, j2, &
                        exp_fac)
!
! !INPUT PARAMETERS:
      integer, intent(in) :: i, j, i1, i2, ju1, j2
      real*8, intent(in) :: tdt    ! model time step (s)
      real*8, intent(in) :: frac1  ! precipitating fraction of grid square
      real*8, intent(in) :: rate1  ! local precipitation rate, i.e., grid-scale precipitation
      real*8, intent(in) :: exp_fac(NPULSE)    ! pulsing decay per time step (day^-1)
                        ! rate divided by frac1 (mm/day)
      real*8, intent(inout) :: soilpuls(NPULSE+1, i1:i2, ju1:j2) ! tracking of wet/dry & three types of pulsing (Y&L, 94)
!
! !RETURNED VALUE:
      real*8  :: Pulseit
!
! !DESCRIPTION:
!  This routine calculates the increase (or "pulse") of soil NOx
!  emission due to precipitation falling over a dry grid square and
!  activating dormant microbes.  Based on section 4.1 of Yienger and Levy,
!  JGR 100, 11,447-11,464, 1995.
!
! \begin{verbatim}
!   soilpuls : tracking of wet/dry & three types of pulsing (Y&L, 94);
!     soilpuls(1,mm)    => flag for wet (-1) or dry (+1) soil; only dry soil
!                          is subject to pulsing
!     soilpuls(nn+1,mm) => fraction of grid square affected by equivalent
!                          fresh pulsing (e.g., sprinkle, shower, heavy rain)
! \end{verbatim}
!
! !LOCAL VARIABLES:
      integer :: nn
      real*8  :: area      ! fraction of grid square affected by pulsing
      real*8  :: rfac
      real*8  :: rsecpdy
      real*8  :: tdt_days
!EOP
!------------------------------------------------------------------------------
!BOC
      Pulseit  = 0.0d0
      rsecpdy  = SECPDY
      tdt_days = tdt / rsecpdy

!!     ==========
!      if (first) then
!!     ==========
!         first = .false.
!
!         exp_fac(1) = Exp (-PULSE_DECAY(1) * tdt_days)
!         exp_fac(2) = Exp (-PULSE_DECAY(2) * tdt_days)
!         exp_fac(3) = Exp (-PULSE_DECAY(3) * tdt_days)
!      end if

      if (soilpuls(1,i,j) <= 0.0d0) then  ! wet, no pulsing

         Pulseit = 1.0d0

      else                               ! dry, subject to pulsing

         do nn = 1, NPULSE
            if (soilpuls(nn+1,i,j) < 1.0d-3) then  ! no pulsing, assume evap.
               soilpuls(nn+1,i,j) = 0.0d0
            else        ! pulse from previous time step decays exponentially
               soilpuls(nn+1,i,j) = soilpuls(nn+1,i,j) * exp_fac(nn)
            end if
         end do

         ! -------------------------------------------------------------------
         ! Determine if a new pulse is to be applied to the grid square due
         ! to precipitation over the current time step.  The pulse is applied
         ! to the grid square fraction, frac1, experiencing precipitation.
         ! Assume a characteristic 1-day duration for precipitation in a given
         ! subgrid area of the grid square, so that the full extent of pulsing
         ! (PULSE_FAC) is realized over 24 hours; for a model time step of
         ! tdt seconds, reduce the pulsing by a factor tdt/SECPDY.
         ! -------------------------------------------------------------------

         rfac = frac1 * tdt_days

         if ((rate1 >=  1.0d0) .and. (rate1 <   5.0d0)) then  ! sprinkle
            soilpuls(2,i,j) = soilpuls(2,i,j) + rfac
         else if ((rate1 >=  5.0d0) .and. (rate1 <  15.0d0)) then  ! shower
            soilpuls(3,i,j) = soilpuls(3,i,j) + rfac
         else if (rate1  >= 15.0d0) then   ! heavy rain
            soilpuls(4,i,j) = soilpuls(4,i,j) + rfac
         end if

         ! ------
         ! Scale.
         ! ------

         area = 0.0d0

         ! ---------------------------------------------------------------------
         ! Add up the contributions of the different types of pulses to obtain
         ! the total pulsing multiplicative factor Pulseit; PULSE_FAC is the
         ! multiplicative factor for fresh pulsing of each type.  Also determine
         ! the fractional grid square area, area, affected by pulsing.  Assume
         ! that the area occupied by the different pulses is additive, i.e.,
         ! that successive pulses apply to different areas of the grid square
         ! and that the area co-occupied by a pulse decreases as the pulsing
         ! decays.  If the resulting area is in excess of unity then the pulsing
         ! must be scaled back to the grid square area.  If the area is less
         ! than unity then we have to account for non-pulsing emissions from the
         ! (1-area) non-pulsing fraction of the grid square.
         ! ---------------------------------------------------------------------

         do nn = 1, NPULSE
            Pulseit = Pulseit + (PULSE_FAC(nn) * soilpuls(nn+1,i,j))
            area    = area + soilpuls(nn+1,i,j)
         end do

         if (area < 1.0d0) then
            Pulseit = Pulseit + 1.0d0 - area
         else
            Pulseit = Pulseit / area
            do nn = 1, NPULSE
               soilpuls(nn+1,i,j) = soilpuls(nn+1,i,j) / area
            end do
         end if

      end if

      return

      end function Pulseit
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Soil_Base
!
! !INTERFACE:
!
      function Soil_Base (nn, soilprep_ij, soilpuls1_ij, rpulse) 

!     ----------------------
!     Argument declarations.
!     ----------------------
!
! !INPUT PARAMETERS:
      integer, intent(in) :: nn           ! soil type
      real*8, intent(in) :: soilprep_ij  ! observed precipitation (mm/day/box)
      real*8, intent(in) :: soilpuls1_ij ! tracking of wet/dry & three types of pulsing (Y&L, 94);
      real*8, intent(in) :: rpulse       ! pulsing rate (units?)
!
! !RETURNED VALUE:
      real*8  :: Soil_Base
!
! !DESCRIPTION:
!   This routine updates the land fractions of three types of pulsings;
!   updates only if within the calculated window.
!
! \begin{verbatim}
!   soilpuls : tracking of wet/dry & three types of pulsing (Y&L, 94);
!     soilpuls(1,mm)    => flag for wet (-1) or dry (+1) soil; only dry soil
!                          is subject to pulsing
!     soilpuls(nn+1,mm) => fraction of grid square affected by equivalent
!                          fresh pulsing (e.g., sprinkle, shower, heavy rain)
! \end{verbatim}
!
!EOP
!------------------------------------------------------------------------------
!BOC
      Soil_Base = 0.0d0

      if (nn == 1) then       ! desert

         Soil_Base = 0.0d0

      else if (nn == 2) then  ! tropical rain forest

         if (soilprep_ij > 1.0d0) then  ! wet
            Soil_Base = SOIL_AW(2)
         else                              ! dry
            Soil_Base = SOIL_AD(2)
         end if

      else if ((nn ==8) .or. (nn == 9)) then

         Soil_Base = SOIL_AW(nn)
         if (nn == 9) then
            Soil_Base = Soil_Base / 30.0d0
         end if

      else  ! other

         if (soilpuls1_ij > 0.0d0) then  ! dry
            Soil_Base = SOIL_AD(nn) * rpulse
         else                              ! wet
            Soil_Base = SOIL_AW(nn)
         end if

      end if

!     --------------
!     Convert units.
!     --------------

      Soil_Base = Soil_Base * UNITCONV

      return

      end function Soil_Base
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Soil_Can_Red_Fac
!
! !INTERFACE:
!
      function Soil_Can_Red_Fac (nn, radiat, windsqr, canopynox, xlai)
!
! !INPUT PARAMETERS:
      integer, intent(in) :: nn       ! soil type
      real*8, intent(in) :: radiat    ! solar radiation (W/m^2)
      real*8, intent(in) :: windsqr   ! surface wind speed squared ((m/s)^2)
      real*8, intent(in) :: canopynox ! deposition rate constant for NOx
      real*8, intent(in) :: xlai      ! leaf area index of land type for 
                                      !   current month
!
! !RETURNED VALUE:
      real*8  :: Soil_Can_Red_Fac
!
! !DECLARED PARAMETERS:
!     -----------------------------------------------------------------
!     Coefficient ALPHA (2.8d-2, 5.6d-3) day, night canopy ventilation;
!     time of 1 hour day, 5 hour night;
!     VFDAY, VFNIGHT -> ALPHA scaled.
!     -----------------------------------------------------------------

      real*8, parameter :: VFDAY   = 1.0d-2 ! ventilation velocity in day   (m/s)
      real*8, parameter :: VFNIGHT = 0.2d-2 ! ventilation velocity in night (m/s)
!
! !DESCRIPTION:
!   This routine calculates the canopy reduction factor.
!
!   Wang et al.: [1998] JGR vol. 103 p10713-10725
!
! !LOCAL VARIABLES:
      real*8  ::  vfnew    ! ventilation rate constant for NOx (m/s)
!EOP
!------------------------------------------------------------------------------
!BOC
      Soil_Can_Red_Fac = 0.0d0

      if (radiat > 0.0d0) then
        vfnew = VFDAY
      else
        vfnew = VFNIGHT
      end if

      if ((xlai > 0.0d0) .and. (canopynox > 0.0d0)) then

        vfnew =  vfnew *  &
     &    Sqrt (windsqr / 9.0d0 * 7.0d0 / xlai) *  &
     &    (SOIL_EXT(2) / SOIL_EXT(nn))

        Soil_Can_Red_Fac = canopynox / (canopynox + vfnew)

      end if

      return

      end function Soil_Can_Red_Fac
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Soil_Fert_Fac
!
! !INTERFACE:
!
      function Soil_Fert_Fac (nn, soilfert_ij)
!
! !INPUT PARAMETERS:
      integer, intent(in) :: nn                     ! soil type
      real*8 , intent(in) :: soilfert_ij            ! fertilizers (ng N/m^2/s)
!
! !RETURNED VALUE:
      real*8  :: Soil_Fert_Fac
!
! !DESCRIPTION:
! This routine calculates the NOx emission from fertilizer.
!
!EOP
!------------------------------------------------------------------------------
!BOC
      Soil_Fert_Fac = 0.0d0

!                                    ======
      if ((nn /= 8) .and. (nn /= 9)) return
!                                    ======

!     -------------------------------------------------
!     Agriculture/fertilizer has to be in "ng N/m^2/s";
!     -------------------------------------------------

      Soil_Fert_Fac = soilfert_ij

      if (nn == 9) then
         Soil_Fert_Fac = Soil_Fert_Fac / 30.0d0
      end if

!     --------------
!     Convert units.
!     --------------

      Soil_Fert_Fac = Soil_Fert_Fac * UNITCONV

      return

      end function Soil_Fert_Fac
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Soil_Temp_Fac
!
! !INTERFACE:
!
      function Soil_Temp_Fac (nn, tempc1, soilpuls1_ij)
!
! !INPUT PARAMETERS:
      integer, intent(in) :: nn           ! soil type
      real*8, intent(in)  :: tempc1       ! surface air temperature (degC)
      real*8, intent(in)  :: soilpuls1_ij ! tracking of wet/dry & three types 
                                          ! of pulsing (Y&L, 94)
!
! !RETURNED VALUE:
      real*8  :: Soil_Temp_Fac
!
! !DESCRIPTION:
!   This routine calculates the soil temperature factor.
!   Yienger and Levy [1995] JGR 100, 11447-11464
!
! \begin{verbatim}
!     soilpuls(1,mm)    => flag for wet (-1) or dry (+1) soil; only dry soil
!                          is subject to pulsing
!     soilpuls(nn+1,mm) => fraction of grid square affected by equivalent
!                          fresh pulsing (e.g., sprinkle, shower, heavy rain)
! \end{verbatim}
!
! !LOCAL VARIABLES:
      real*8  :: xtempc1
!
!EOP
!------------------------------------------------------------------------------
!BOC
      Soil_Temp_Fac = 0.0d0

      xtempc1 = tempc1

      if (nn <= 2) then  ! desert and rain forest
         Soil_Temp_Fac = 1.0d0
      else if ((soilpuls1_ij > 0.0d0) .and. (nn /= 8) .and. (nn /= 9)) then

         ! -------------------------------------------------------
         ! Dry:
         !   surface temperature -> soil temperature;
         !   convert the lowest model level air temp to soil temp;
         !   based on observations of Johansson et. al. [1988];
         !   add 5 degrees C to model temperature.
         ! -------------------------------------------------------

         xtempc1 = xtempc1 + 5.0d0
         if (xtempc1 > 30.0d0) then       ! optimal
            Soil_Temp_Fac = 1.0d0
         else if (xtempc1 >  0.0d0) then  ! cold-linear
            Soil_Temp_Fac = xtempc1 / 30.0d0
         else
            Soil_Temp_Fac = 0.0d0
         end if
      else

         ! ---------------------------------------------------------------------
         ! Wet:
         !   surface temperature -> soil temperature;
         !   convert the lowest model level air temp to soil temp;
         !   use the empirical relationships derived by Williams et al. [1992b];
         !   ECO system dependent.
         ! ---------------------------------------------------------------------

         xtempc1 = SOIL_T2(nn) + (SOIL_T1(nn) * xtempc1)
         if (xtempc1 >= 30.0d0) then       ! optimal
            Soil_Temp_Fac = 21.97d0
         else if (xtempc1 >= 10.0d0) then  ! exponential
            Soil_Temp_Fac = Exp (0.103d0 * xtempc1)
         else if (xtempc1 >   0.0d0) then  ! cold-linear
            Soil_Temp_Fac = 0.28d0 * xtempc1
         else
            Soil_Temp_Fac = 0.0d0
         end if

      end if

      return

      end function Soil_Temp_Fac
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Soil_Type
!
! !INTERFACE:
!
      subroutine Soil_Type (iday, idaySoilType, soilprep, soilpuls, &
                            pr_diag, loc_proc, i1, i2, ju1, j2)
!
! !INPUT PARAMETERS:
      logical, intent(in) :: pr_diag
      integer, intent(in) :: loc_proc
      integer, intent(in) :: i1, i2, ju1, j2
      integer, intent(in) :: iday   ! day of the current month
               ! observed precipitation (mm/day/box)
      real*8 , intent(in) :: soilprep(i1:i2, ju1:j2) 
!
! !INPUT/OUTPUT PARAMETERS:
!   soilpuls : tracking of wet/dry & three types of pulsing (Y&L, 94);
!     soilpuls(1,mm)    => flag for wet (-1) or dry (+1) soil; only dry soil
!                          is subject to pulsing
!     soilpuls(nn+1,mm) => fraction of grid square affected by equivalent
!                          fresh pulsing (e.g., sprinkle, shower, heavy rain)
      real*8, intent(inOut) :: soilpuls(NPULSE+1, i1:i2, ju1:j2)
      integer, intent(inOut) :: idaySoilType
!
! !DECLARED PARAMETERS:
      integer, parameter :: IDAYS_TO_TEST = 14      ! number of days for pulse

      real*8,  parameter :: WETSOIL       = 10.0d0  ! criteria for wet soil (mm);
                                                    ! above 10 mm for two weeks
!
! !DESCRIPTION:
!  Determines whether soil is dry or wet for all land grid boxes; updated daily.
!
! !LOCAL VARIABLES:
      integer :: nn, ii, ij
      integer :: ncurr          ! number of days in current  month
      integer :: nprev          ! number of days in previous month

      real*8  :: rain           ! total rain (units?)
      real*8  :: rncurr
      real*8  :: rnprev
!EOP
!------------------------------------------------------------------------------
!BOC
      if (pr_diag) Write (6,*) 'Soil_Type called by ', loc_proc

      if (iday /= idaySoilType) then
         idaySoilType  = iday

         if (iday >= IDAYS_TO_TEST) then
             ncurr = IDAYS_TO_TEST
             nprev = 0
         else
            ncurr = iday
            nprev = IDAYS_TO_TEST - iday
         end if

         rncurr = ncurr
         rnprev = nprev

         do ij = ju1, j2
            do ii = i1, i2
               rain = IDAYS_TO_TEST*soilprep(ii,ij)

               if (rain > WETSOIL) then  ! wet
                  soilpuls(1,ii,ij) = -1.0d0
                  do nn = 1, NPULSE
                     soilpuls(nn+1,ii,ij) = 0.0d0
                  end do
               else                      ! dry
                  soilpuls(1,ii,ij) = 1.0d0
               end if
           end do
        end do

      end if

      return

      end subroutine Soil_Type
!EOC
!------------------------------------------------------------------------------
end module GmiSoilEmission_mod
