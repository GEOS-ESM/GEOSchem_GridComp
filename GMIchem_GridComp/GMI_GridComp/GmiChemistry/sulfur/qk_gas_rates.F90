!=============================================================================
!
! $Id$
!
! CODE DEVELOPER
!   originated from GRANTOUR programmed by John Walton
!   modified by Xiaohong Liu
!
! FILE
!   qk_gas_rates.F
!
! ROUTINES
!   Do_QK_Gas_Rates
!
!=============================================================================

!-----------------------------------------------------------------------------
!
! ROUTINE
!   Do_QK_Gas_Rates
!
! DESCRIPTION
!   This is the routine for the gas chemical reaction rates
!
! ARGUMENTS
!   itloop        : # of zones (ilong * ilat * ivert)
!   tgcm          : temperature (K)
!   zmair         : air density (#/cm3)
!   qh1p, qh2p, qh3p, qh5p, qh8p, qh9p
!                 : gas reaction rate coefficients
!
! ---------------------------------------------------------------------------

      subroutine Do_QK_Gas_Rates  &
     &  (itloop, tgcm, zmair,  &
     &   qh1p, qh2p, qh3p, qh5p, qh8p, qh9p)

      implicit none

!     ----------------------
!     GMIMOD-related values.
!     ----------------------

      integer :: itloop

      real*8 ::  tgcm (itloop)
      real*8 ::  zmair(itloop)
      real*8 ::  qh1p (itloop)
      real*8 ::  qh2p (itloop)
      real*8 ::  qh3p (itloop)
      real*8 ::  qh5p (itloop)
      real*8 ::  qh8p (itloop)
      real*8 ::  qh9p (itloop)


!     ----------------
!     Begin execution.
!     ----------------

!.... reaction rate coefficients in the gas phase
!     (qh1p, qh2p, qh3p, qh5p, qh8p, qh9p)

!          ======
      call chemk1  &
!          ======
     &  (tgcm, itloop, qh1p)

!          ======
      call chemk2  &
!          ======
     &  (tgcm, itloop, qh2p, zmair)

!          ======
      call chemk3  &
!          ======
     &  (tgcm, itloop, qh3p, zmair)

!          ======
      call chemk5  &
!          ======
     &  (tgcm, itloop, qh5p)

!          ======
      call chemk8  &
!          ======
     &  (tgcm, itloop, qh8p, zmair)

!          ======
      call chemk9  &
!          ======
     &  (tgcm, itloop, qh9p)


      return

      end


!-----------------------------------------------------------------------------
!
! ROUTINE
!   chemk1
!
! DESCRIPTION
!
!     *****************************************************************
!     * reaction rate for DMS + OH -> SO2
!     *****************************************************************
!
!     *     tpp  : air temperature (K)
!     *     qh1p : reaction rate coefficient (cm3 molecules-1 s-1)
!
!-----------------------------------------------------------------------------

      subroutine chemk1  &
     &  (tpp, itloop, qh1p)

      implicit none

#     include "sulfchem.h"

!     ----------------------
!     GMIMOD-related values.
!     ----------------------

      integer :: itloop

      real*8  :: tpp (itloop)
      real*8  :: qh1p(itloop)

! -------------------------------------------------------------
!  *** LOCAL  VARIABLES WERE DECLARED HERE.
! -------------------------------------------------------------

      real*8  c2, c3, c4, c5

      data c2 /8.46d-10/, c3 /2.68d-10/, c4 /88.1d0/, c5 /1.04d11/

      real*8 :: t1(itloop)
      real*8 :: t2(itloop)
      real*8 :: t3(itloop)
      real*8 :: t4(itloop)

      t1 = 0.0d0
      t2 = 0.0d0
      t3 = 0.0d0
      t4 = 0.0d0


!     *****************************************************************
!     * DMS + OH --> SO2 rate coef.
!     *****************************************************************

      do 100 ijk = 1, itloop

        t1(ijk) = exp(-234.d0 / tpp(ijk))
        t2(ijk) = exp(7230.d0 / tpp(ijk))
        t3(ijk) = exp(7810.d0 / tpp(ijk))
        t4(ijk) = exp(7460.d0 / tpp(ijk))

  100 continue

      do 200 ijk = 1, itloop

        qh1p(ijk) = (tpp(ijk) * t1(ijk) + c2 * t2(ijk) + c3 * t3(ijk))  &
     &            / (c4 * t4(ijk) + c5 * tpp(ijk))

  200 continue


      return

      end


!-----------------------------------------------------------------------------
!
! ROUTINE
!   chemk2
!
! DESCRIPTION
!
!
!     *****************************************************************
!     * reaction rate for SO2 + OH -> H2SO4
!     *****************************************************************
!
!     *     tpp  : air temperature (K)
!     *     zmp  : air density (molecules cm-3)
!     *     qh2p : reaction rate coefficient (cm3 molecules-1 s-1)
!
!     update reaction rate according to DeMore et al., 1997
!     qhzro300 = 3.0d-31 cm6 molecules-2 s-1
!     qhinf    = 1.5d-12 cm3 molecules-1 s-1
!     fc       = 0.6
!-----------------------------------------------------------------------------

      subroutine chemk2  &
     &  (tpp, itloop, qh2p, zmp)

      implicit none

#     include "sulfchem.h"

!     ----------------------
!     GMIMOD-related values.
!     ----------------------

      integer :: itloop

      real*8  :: tpp (itloop)
      real*8  :: qh2p(itloop)
      real*8  :: zmp (itloop)

! -------------------------------------------------------------
!  *** LOCAL  VARIABLES WERE DECLARED HERE.
! -------------------------------------------------------------

      real*8  &
     &  qhzro300, qhinf, en, fc

!     data qhzro300 /3.16d-31/, qhinf /2.d-12/, en /3.3d0/, fc /0.45d0/
      data qhzro300 /3.0d-31/, qhinf /1.5d-12/, en /3.3d0/, fc /0.6d0/

      real*8 :: t1(itloop)
      real*8 :: t2(itloop)
      real*8 :: t3(itloop)
      real*8 :: t4(itloop)
      real*8 :: t5(itloop)

      t1 = 0.0d0
      t2 = 0.0d0
      t3 = 0.0d0
      t4 = 0.0d0
      t5 = 0.0d0

!     *****************************************************************
!     * SO2 + OH --> H2SO4 rate coef.
!     *****************************************************************

      do 100 ijk = 1, itloop

        t1(ijk) = qhzro300 * (300.d0 / tpp(ijk)) ** en

  100 continue

      do 200 ijk = 1, itloop

        t3(ijk) = t1(ijk) * zmp(ijk) / qhinf

  200 continue

      do 300 ijk = 1, itloop

        t4(ijk) = t1(ijk) * zmp(ijk) / (1.d0 + t3(ijk))
        t5(ijk) = 1.d0 / (1.d0 + log10(t3(ijk))**2)

  300 continue

      do 400 ijk = 1, itloop

        qh2p(ijk) = t4(ijk) * fc ** t5(ijk)

  400 continue


      return

      end


!-----------------------------------------------------------------------------
!
! ROUTINE
!   chemk3
!
! DESCRIPTION
!
!
!     *****************************************************************
!     * reaction rate for HO2 + HO2 -> H2O2 + O2
!     *****************************************************************
!     *
!     *     tpp  : air temperature (K)
!     *     zmp  : air density (molecules cm-3)
!     *     qh3p : reaction rate coefficient (cm3 molecules-1 s-1)
!
!-----------------------------------------------------------------------------

      subroutine chemk3  &
     &  (tpp, itloop, qh3p, zmp)

      implicit none

#     include "sulfchem.h"

!     ----------------------
!     GMIMOD-related values.
!     ----------------------

      integer :: itloop

      real*8  :: tpp (itloop)
      real*8  :: qh3p(itloop)
      real*8  :: zmp (itloop)

!     *****************************************************************
!     * HO2 + HO2 --> H2O2 + O2
!     *****************************************************************

      do 100 ijk = 1, itloop

        qh3p(ijk) = 2.3d-13 * exp(600.d0 / tpp(ijk))  &
     &            + 1.7d-33 * zmp(ijk) * exp(1000.d0 / tpp(ijk))

  100 continue


      return

      end


!-----------------------------------------------------------------------------
!
! ROUTINE
!   chemk5
!
! DESCRIPTION
!
!     *****************************************************************
!     * reaction rate for H2O2 + OH -> H2O + HO2
!     *****************************************************************
!
!     *     tpp  : air temperature (K)
!     *     qh5p : reaction rate coefficient (cm3 molecules-1 s-1)
!
!-----------------------------------------------------------------------------

      subroutine chemk5  &
     &  (tpp, itloop, qh5p)

      implicit none

#     include "sulfchem.h"

!     ----------------------
!     GMIMOD-related values.
!     ----------------------

      integer :: itloop

      real*8  :: tpp (itloop)
      real*8  :: qh5p(itloop)


!     *****************************************************************
!     * H2O2 + OH --> H2O + HO2
!     *****************************************************************

      do 100 ijk = 1, itloop

        qh5p(ijk) = 2.9d-12 * exp(-160.d0 / tpp(ijk))

  100 continue


      return

      end


!-----------------------------------------------------------------------------
!
! ROUTINE
!   chemk8
!
! DESCRIPTION
!
!
!     *****************************************************************
!     * reaction rate for HO2 + HO2 + H2O --> H2O2 + O2 + H2O
!     *****************************************************************
!
!     *     tpp : air temperature (K)
!     *     zmp : air density (molecules cm-3)
!     *     qh8p: reaction rate coefficient (cm6 molecules^-2 s-1)
!
!-----------------------------------------------------------------------------

      subroutine chemk8  &
     &  (tpp, itloop, qh8p, zmp)

      implicit none

#     include "sulfchem.h"

!     ----------------------
!     GMIMOD-related values.
!     ----------------------

      integer :: itloop

      real*8  :: tpp (itloop)
      real*8  :: qh8p(itloop)
      real*8  :: zmp (itloop)

!     *****************************************************************
!     * HO2 + HO2 + H2O --> H2O2 + O2 + H2O
!     *****************************************************************

      do 100 ijk = 1, itloop

        qh8p(ijk) = 3.22d-34 * exp(2800.d0 / tpp(ijk))  &
     &            + 2.38d-54 * exp(3200.d0 / tpp(ijk)) * zmp(ijk)

  100 continue


      return

      end


!-----------------------------------------------------------------------------
!
! ROUTINE
!   chemk9
!
! DESCRIPTION
!
!     *****************************************************************
!     * reaction rate for DMS + NO3 -> SO2 + ...
!     *****************************************************************
!
!     *     qh9p: reaction rate coefficient (cm3 molecules-1 s-1)
!
!-----------------------------------------------------------------------------

      subroutine chemk9  &
     &  (tpp, itloop, qh9p)

      implicit none

#     include "sulfchem.h"

!     ----------------------
!     GMIMOD-related values.
!     ----------------------

      integer :: itloop

      real*8  :: tpp (itloop)
      real*8  :: qh9p(itloop)


!     *****************************************************************
!     * DMS + NO3 --> SO2 + ....
!     *****************************************************************

      do 100 ijk = 1, itloop

!       qh9p(ijk) = 1.9d-13 * exp(520.d0 / tpp(ijk))
!       qh9p(ijk) = 1.9d-13 * exp(500.d0 / tpp(ijk))  ! update 10/2002

!     &################################################

!       Modification on 97/05/14 to turn off the NO3 reaction

        qh9p(ijk) = 0.d0

!     &################################################

  100 continue


      return

      end

