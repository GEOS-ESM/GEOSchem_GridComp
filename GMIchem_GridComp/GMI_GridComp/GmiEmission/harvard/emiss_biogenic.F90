
!=============================================================================
!
! $Id$
!
! CODE DEVELOPER
!   John Tannahill, LLNL
!   jrt@llnl.gov
!
!   Original code from:
!     Harvard tropospheric emissions module for 3D applications;
!       by Yuhang Wang, Gerry Gardner, and Prof. Daniel Jacob
!       of Harvard University (Release V1.0)
!
! FILE
!   emiss_biogenic.F
!
! ROUTINES
!   Biogenic_Base
!   Biogenic_Isop
!   Biogenic_Monot
!   Biogenic_Temp_Fac
!
!=============================================================================


!-----------------------------------------------------------------------------
!
! ROUTINE
!   Biogenic_Base
!
! DESCRIPTION
!   This routine calculates the baseline emissions for isoprene and
!   monoterpenes.
!
!   Isoprene is traced in terms of equivalent Carbon atoms.
!
! ARGUMENTS
!   ireg  : number of land types in a grid square
!   iland : land type id in grid square for ireg land types
!   tdt   : model time step (s)
!   convert_isop  : isoprene    emissions by landtype   (atomsC/cm^2 leaf/s)
!   convert_monot : monoterpene emissions by landtype   (atomsC/cm^2 leaf/s)
!   mcor          : grid box area (m^2)
!   base_isop     : baseline emissions for isoprene     (kgC/box/step?)
!   base_monot    : baseline emissions for monoterpenes (kgC/box/step?)
!
!-----------------------------------------------------------------------------

      subroutine Biogenic_Base  &
     &  (ireg, iland, tdt, convert_isop, convert_monot, mcor,  &
     &   base_isop, base_monot, pr_diag, loc_proc, i1, i2, ju1, j2)

!      use m_PhysicalConstants, only : KGPG
!      use m_PhysicalConstants, only : AVOGAD
!      use m_PhysicalConstants, only : CMPM

      implicit none

#     include "gmi_phys_constants.h"
#     include "gmi_emiss_constants.h"

!     ----------------------
!     Argument declarations.
!     ----------------------

      logical, intent(in) :: pr_diag
      integer, intent(in) :: loc_proc, i1, i2, ju1, j2
      integer :: ireg (i1:i2, ju1:j2)
      integer :: iland(i1:i2, ju1:j2, NTYPE)
      real*8  :: tdt
      real*8  :: convert_isop (NVEGTYPE)
      real*8  :: convert_monot(NVEGTYPE)
      real*8  :: mcor      (i1:i2, ju1:j2)
      real*8  :: base_isop (i1:i2, ju1:j2, NTYPE)
      real*8  :: base_monot(i1:i2, ju1:j2, NTYPE)


!     -----------------------
!     Parameter declarations.
!     -----------------------

!     ----------------------
!     Variable declarations.
!     ----------------------

      integer :: il, ij, it

      real*8  :: rfac        ! conversion factor from atomsC/cm^2/s to
                             ! kgC/box/step?
      real*8  :: cmpm2


!     ----------------
!     Begin execution.
!     ----------------

      if (pr_diag) then
        Write (6,*) 'Biogenic_Base called by ', loc_proc
      end if

      cmpm2 = CMPM * CMPM

      do ij = ju1, j2
        do il = i1, i2

          rfac = (MWTCARBON * KGPG) * (mcor(il,ij) * cmpm2) *  &
     &           tdt / AVOGAD

          do it = 1, ireg(il,ij)

!           --------------------------------------
!           ijland+1 is the Olson land type index.
!           --------------------------------------

            base_isop (il,ij,it) =  &
     &        convert_isop (iland(il,ij,it)+1) * rfac

            base_monot(il,ij,it) =  &
     &        convert_monot(iland(il,ij,it)+1) * rfac

          end do

        end do
      end do


      return

      end


!-----------------------------------------------------------------------------
!
! ROUTINE
!   Biogenic_Isop
!
! DESCRIPTION
!   This routine calculates the biogenic isoprene emissions in a grid box.
!
! ARGUMENTS
!   ireg1     : number of Olson land types in the grid box
!   iuse1     : Olson land type fraction   in the grid box (mils?)
!   cfrac1    : fractional cloud cover
!   suncos1   : cosine of the solar zenith angle
!   tempk1    : surface air temperature (degK)
!   sopcoeff  : 2nd order polynomial coefficients for light correction
!   baseisop1 : baseline emissions for isoprene (kgC/box/step)
!   xlai1     : leaf area index of land type for current month
!
!-----------------------------------------------------------------------------

      function Biogenic_Isop  &
     &  (ireg1, iuse1, cfrac1, suncos1, tempk1, sopcoeff,  &
     &   baseisop1, xlai1)

!      use m_PhysicalConstants, only : KGPG
!      use m_PhysicalConstants, only : AVOGAD
!      use m_PhysicalConstants, only : ABS_ZERO

      implicit none

#     include "gmi_phys_constants.h"
#     include "gmi_emiss_constants.h"


!     ----------------------
!     Argument declarations.
!     ----------------------

      integer :: ireg1
      integer :: iuse1(NTYPE)
      real*8  :: cfrac1
      real*8  :: suncos1
      real*8  :: tempk1
      real*8  :: sopcoeff (NPOLY)
      real*8  :: baseisop1(NTYPE)
      real*8  :: xlai1    (NTYPE)


!     ----------------------
!     Function declarations.
!     ----------------------

      real*8  :: Biogenic_Isop

      real*8, external  :: Biofit
      real*8, external  :: Biogenic_Temp_Fac


!     -----------------------
!     Parameter declarations.
!     -----------------------

!     ----------------------------------------------------
!     XNUMOL_C : number of molecules C / kg C (atomsC/kgC)
!     ----------------------------------------------------

      real*8, parameter ::  &
     &  XNUMOL_C = AVOGAD / (MWTCARBON * KGPG)


!     ----------------------
!     Variable declarations.
!     ----------------------

      integer :: inveg

      real*8  :: clight
      real*8  :: embio
      real*8  :: riuse1
      real*8  :: tlai


!     ----------------
!     Begin execution.
!     ----------------

      Biogenic_Isop = 0.0d0


!     ------------------------------------------------------
!     Calculate total of leaf area index * baseline isoprene
!     over all Olson land types that are in this grid box.
!     ------------------------------------------------------

      tlai = 0.0d0

      do inveg = 1, ireg1

        tlai = tlai + (xlai1(inveg) * baseisop1(inveg))

      end do


!     ----------------------------------------------------------------
!     Apply light & temperature corrections to baseline emissions only
!     if it is daytime and if there is nonzero isoprene emission
!     (i.e., xlai1 * baseisop1 > 0).
!     ----------------------------------------------------------------

      if ((suncos1 > 0.0d0) .and. (tlai > 0.0d0)) then

        embio = 0.0d0

!       ------------------------------------------------
!       Loop over all Olson land types in this grid box.
!       ------------------------------------------------

        do inveg = 1, ireg1

          if ((xlai1(inveg) * baseisop1(inveg)) > 0.0d0) then

!           ---------------------------------------------
!           Calculate light correction -- polynomial fit.
!           ---------------------------------------------

            clight =  &
     &        Biofit (cfrac1, suncos1, xlai1(inveg), sopcoeff)

!           ------------------------------------------------------
!           Apply light correction to baseline isoprene emissions.
!           Also multiply by the fraction of the grid box occupied
!           by this Olson land type.
!           ------------------------------------------------------

            riuse1 = iuse1(inveg)

            embio  = embio +  &
     &               ((baseisop1(inveg) * clight * riuse1) / 1000.0d0)

          end if

        end do

!       --------------------------------------------------------------
!       Apply the temperature correction from Gunther et al. (1992) to
!       the isoprene emissions.
!       --------------------------------------------------------------

        if ((tempk1 + ABS_ZERO) > 0.0d0) then
          Biogenic_Isop = Biogenic_Temp_Fac (tempk1) * embio
        else
          Biogenic_Isop = 0.0d0
        end if

      end if


!     ---------------------------------------------
!     Convert from kgC/box/step to atomsC/box/step.
!     ---------------------------------------------

      Biogenic_Isop = Biogenic_Isop * XNUMOL_C


      return

      end


!-----------------------------------------------------------------------------
!
! ROUTINE
!   Biogenic_Monot
!
! DESCRIPTION
!   This routine calculates the biogenic monoterpene emissions in a grid box.
!
! ARGUMENTS
!   ireg1      : number of Olson land types in the grid box
!   iuse1      : Olson land type fraction   in the grid box (mils?)
!   tempk1     : surface air temperature (degK)
!   basemonot1 : baseline emissions for monoterpene (kgC/box/step)
!   xlai1      : leaf area index of land type for current month
!
!-----------------------------------------------------------------------------

      function Biogenic_Monot  &
     &  (ireg1, iuse1, tempk1, basemonot1, xlai1)

!      use m_PhysicalConstants, only : KGPG
!      use m_PhysicalConstants, only : AVOGAD

      implicit none

#     include "gmi_phys_constants.h"
#     include "gmi_emiss_constants.h"


!     ----------------------
!     Argument declarations.
!     ----------------------

      integer :: ireg1
      integer :: iuse1     (NTYPE)
      real*8  :: tempk1
      real*8  :: basemonot1(NTYPE)
      real*8  :: xlai1     (NTYPE)


!     ----------------------
!     Function declarations.
!     ----------------------

      real*8  :: Biogenic_Monot


!     -----------------------
!     Parameter declarations.
!     -----------------------

!     ----------------------------------------------------------------
!     Temperature correction coefficients from Guenther et al. (1995).
!
!     BETA : empirical coefficient given by Guenther (degK^-1)
!     TS   : leaf temperature at standard conditions (degK)
!     ----------------------------------------------------------------

      real*8, parameter ::  &
     &  BETA =   0.09d0,  &
     &  TS   = 303.0d0

!     ----------------------------------------------------
!     XNUMOL_C : number of molecules C / kg C (atomsC/kgC)
!     ----------------------------------------------------

      real*8, parameter ::  &
     &  XNUMOL_C = AVOGAD / (MWTCARBON * KGPG)


!     ----------------------
!     Variable declarations.
!     ----------------------

      integer :: inveg

      real*8  :: embio
      real*8  :: riuse1


!     ----------------
!     Begin execution.
!     ----------------

      Biogenic_Monot = 0.0d0

      embio = 0.0d0


!     ------------------------------------------------
!     Loop over all Olson land types in this grid box.
!     ------------------------------------------------

      do inveg = 1, ireg1

!       ----------------------------------------------------------------
!       Calculate monoterpene emissions for the grid box (kgC/box/step).
!
!       Monoterpenes are now scaled to the leaf area index, and
!       multiplied by the fraction of the grid box occupied by this
!       Olson land type.
!       ----------------------------------------------------------------

        riuse1 = iuse1(inveg)

        embio  =  &
     &    embio +  &
     &    ((basemonot1(inveg) * xlai1(inveg) * riuse1) / 1000.0d0)

      end do


!     --------------------------------------------------------------------
!     Temperature correct monoterpene emissions (Guenther et al., (1995)).
!     Foliar density is accounted for in monoterpene emissions table.
!     --------------------------------------------------------------------

      Biogenic_Monot = embio * Exp (BETA * (tempk1 - TS))


!     ---------------------------------------------
!     Convert from kgC/box/step to atomsC/box/step.
!     ---------------------------------------------

      Biogenic_Monot = Biogenic_Monot * XNUMOL_C


      return

      end


!-----------------------------------------------------------------------------
!
! ROUTINE
!   Biogenic_Temp_Fac
!
! DESCRIPTION
!   This routine corrects the temperature for isoprene emissions,
!   Guenther et al.(92).
!
! ARGUMENTS
!   tempk1 : temperature to be corrected (degK)
!
!-----------------------------------------------------------------------------


      function Biogenic_Temp_Fac (tempk1)

!      use m_PhysicalConstants, only : GAS_CONST_J

      implicit none

#     include "gmi_phys_constants.h"


!     ----------------------
!     Argument declarations.
!     ----------------------

      real*8  :: tempk1


!     ----------------------
!     Function declarations.
!     ----------------------

      real*8  :: Biogenic_Temp_Fac


!     -----------------------
!     Parameter declarations.
!     -----------------------

      real*8, parameter ::  &
     &  CT1 =  95000.0d0,  &
     &  CT2 = 230000.0d0

      real*8, parameter ::  &
     &  T1  =    303.0d0,  &
     &  T3  =    314.0d0


!     ----------------------
!     Variable declarations.
!     ----------------------

      real*8  :: rfac


!     ----------------
!     Begin execution.
!     ----------------

      rfac = GAS_CONST_J * T1 * tempk1


      Biogenic_Temp_Fac =  &
     &           Exp (CT1 / rfac * (tempk1 - T1)) /  &
     &  (1.0d0 + Exp (CT2 / rfac * (tempk1 - T3)))


      return

      end

