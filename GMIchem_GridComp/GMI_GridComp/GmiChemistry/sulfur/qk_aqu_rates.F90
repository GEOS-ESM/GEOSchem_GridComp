!=============================================================================
!
! $Id$
!
! CODE DEVELOPER
!   originated from GRANTOUR programmed by John Walton
!   modified by Xiaohong Liu
!
! FILE
!   qk_aqu_rates.F
!
! ROUTINES
!   Do_QK_Aqu_Rates
!
!=============================================================================

!-----------------------------------------------------------------------------
!
! ROUTINE
!   Do_QK_Aqu_Rates
!
! DESCRIPTION
!   This is the routine for the aqueous reaction rates,
!     Henry's law constants, and aqueous equation coefficients
!
! ARGUMENTS
!   itloop        : # of zones (ilong * ilat * ivert)
!   tgcm          : temperature (K)
!   hpp           : 10**(-dph), dph: mean aqueous pH (=4.5)
!   rhpp          : 1./hpp
!   numcl         : number of grid points in cloud
!   idxcl         : array of the numbers of the numcl points in cloud
!   oh, ho2, h2o, no3, o3
!                 : gas species concentration (molecules/cm-3)
!   cwac          : in-cloud cloud water mixing ratio (g/g air)
!   zmair         : air density (#/cm3)
!   qh1p, qh2p, qh3p, qh5p, qh8p, qh9p
!                 : gas-phase reaction rate coefficients
!   qj4p          : H2O2 photolysis rate coefficient (1/s)
!   qh6, qh7p     : aqueous reaction rate coefficients (L/mol/s)
!   hwp, hyp, h3p : Henry's law constants (mol/L/atm)
!   acoef         : aqueous equation coefficients (a1-a9)
!
! ---------------------------------------------------------------------------

      subroutine Do_QK_Aqu_Rates  &
     &  (itloop, tgcm, hpp, rhpp, numcl, idxcl,  &
     &   oh, ho2, h2o, no3, o3,  &
     &   cwac, zmair,  &
     &   qh1p, qh2p, qh3p, qj4p, qh5p, qh8p, qh9p,  &
     &   qh6, qh7p, hwp, hyp, h3p,  &
     &   acoef)

      implicit none

#     include "sulfchem.h"

!     ----------------------
!     GMIMOD-related values.
!     ----------------------

      integer :: itloop, numcl

      integer :: idxcl(itloop)

      real*8 :: qh6

      real*8 :: tgcm (itloop)
      real*8 :: oh   (itloop)
      real*8 :: ho2  (itloop)
      real*8 :: h2o  (itloop)
      real*8 :: no3  (itloop)
      real*8 :: o3   (itloop)
      real*8 :: cwac (itloop)
      real*8 :: zmair(itloop)
      real*8 :: qh1p (itloop)
      real*8 :: qh2p (itloop)
      real*8 :: qh3p (itloop)
      real*8 :: qj4p (itloop)
      real*8 :: qh5p (itloop)
      real*8 :: qh7p (itloop)
      real*8 :: qh8p (itloop)
      real*8 :: qh9p (itloop)
      real*8 :: hwp  (itloop)
      real*8 :: hyp  (itloop)
      real*8 :: h3p  (itloop)
      real*8 :: acoef(itloop, 9)

!     ----------------
!     Begin execution.
!     ----------------

!.... aqueous reaction rate (qh6) a function of pH only, therefore
!     constant for all grid points: SO2(aq) + H2O2(aq) ->
!
!          ======
      call chemk6  &
!          ======
     &  (qh6, hpp)


!.... aqueous reaction rates (qh7p): SO2(aq) + O3(aq) ->
!
!          ======
      call chemk7  &
!          ======
     &  (tgcm, itloop, qh7p, hpp)


!.... Henry's law coefficients (hwp, hyp, h3p): H2O2, SO2, O3
!
!          ======
      call chem22  &
!          ======
     &  (tgcm, itloop, numcl, idxcl,  &
     &   hwp, hyp, h3p, rhpp)


!.... aqueous equation coefficients (acoef: a1, ..., a9)
!
!          ======
      call chem23  &
!          ======
     &  (itloop, numcl, idxcl,  &
     &   oh, ho2, h2o, no3, o3, hwp, hyp, h3p,  &
     &   cwac, zmair,  &
     &   qh1p, qh2p, qh3p, qj4p, qh5p, qh6, qh7p, qh8p, qh9p,  &
     &   acoef)


      return

      end


!-----------------------------------------------------------------------------
!
! ROUTINE
!   chemk6
!
! DESCRIPTION
!
!
!     *****************************************************************
!     * reaction rate for SO2(aq) + H2O2(aq) -> H2SO4
!     *****************************************************************
!
!     ********** comments:
!     * reaction rate depends only on pH of the droplets which is
!       assumed constant in all clouds
!
!     *     qh6 : aqueous reaction rate coefficient (L/mol/s)
!     *     hpp : hydrogen ion concentration in drops, 10**(-pH)
!
!-----------------------------------------------------------------------------

      subroutine chemk6  &
     &  (qh6, hpp)

      implicit none

#     include "sulfchem.h"

      real*8  :: qh6

! -------------------------------------------------------------
!  *** LOCAL  VARIABLES WERE DECLARED HERE.
! -------------------------------------------------------------

      real*8 qk1f, qk1rk2, qka1, qka2

      data qk1f /5.2d6/, qk1rk2 /0.1d0/,  &
     &     qka1 /0.016d0/, qka2 /1.d-7/


      qh6 = qk1f * hpp /  &
     &      ((qk1rk2 + hpp) * (1.d0 + hpp/qka1 + qka2/hpp))


      return

      end


!-----------------------------------------------------------------------------
!
! ROUTINE
!   chemk7
!
! DESCRIPTION
!
!
!     *****************************************************************
!     * reaction rate for SO2(aq) + O3(aq) -> H2SO4 + O2
!     *****************************************************************
!
!     *     tpp  : air temperature (K)
!     *     qh7p : aqueous reaction rate coefficient (L/mol/s)
!     *     hpp  : hydrogen ion concentration in drops, 10**(-pH)
!
!     *----------------------------------------------------------------

      subroutine chemk7  &
     &  (tpp, itloop, qh7p, hpp)

      implicit none

#     include "sulfchem.h"

!     ----------------------
!     GMIMOD-related values.
!     ----------------------

      integer :: itloop

      real*8  :: tpp (itloop)
      real*8  :: qh7p(itloop)

! -------------------------------------------------------------
!  *** LOCAL  VARIABLES WERE DECLARED HERE.
! -------------------------------------------------------------

      real*8  :: qk511, qk512, qk521, qk522

      data qk511 / 4.39d11/, qk512 / -4131.d0/,  &
     &     qk521 / 2.56d3/,  qk522 /  -996.d0/

      real*8 :: t1(itloop)
      real*8 :: t2(itloop)

      t1 = 0.0d0
      t2 = 0.0d0

!     *****************************************************************
!     * SO2(aq) + O3(aq) --> H2SO4(aq) + O2
!     *****************************************************************

      do 100 ijk = 1, itloop

        t1(ijk) = qk511 * exp (qk512 / tpp(ijk))
        t2(ijk) = qk521 * exp (qk522 / tpp(ijk))

  100 continue

      do 200 ijk = 1, itloop

        qh7p(ijk) = t1(ijk) + t2(ijk) / hpp

  200 continue


      return

      end


!-----------------------------------------------------------------------------
!
! ROUTINE
!   chem22
!
! DESCRIPTION
!
!     *****************************************************************
!     * calculate the henry's law coefficients for H2O2, SO2, and O3
!     *****************************************************************
!
!     *    tpp   :  air temperature (K)
!     *    numcl :  number of grid points in cloud
!     *    idxcl :  array of the numbers of the numcl points in cloud
!     *    rhpp  :  1./10**(-dph), dph: mean aqueous pH (4.5)
!     *    hwp   :  henry's law constant for H2O2
!     *    hyp   :  henry's law constant for SO2
!     *    h3p   :  henry's law constant for O3
!
!     *----------------------------------------------------------------

      subroutine chem22  &
     &   (tpp, itloop, numcl, idxcl,  &
     &    hwp, hyp, h3p, rhpp)

      implicit none

#     include "sulfchem.h"


!     ----------------------
!     GMIMOD-related values.
!     ----------------------

      integer :: itloop, numcl

      integer :: idxcl(itloop)

      real*8  :: tpp(itloop)
      real*8  :: hwp(itloop)
      real*8  :: hyp(itloop)
      real*8  :: h3p(itloop)

! -------------------------------------------------------------
!  *** LOCAL  VARIABLES WERE DECLARED HERE.
! -------------------------------------------------------------

      real*8  e(5), f(5)

      data e / 1.23d0, 0.01232d0, 6.61d-8, 7.45d4, 0.011d0 /
      data f / 3020.d0, 1960.d0, 1500.d0, 6620.d0, 2300.d0 /

      real*8 :: t1(itloop)
      real*8 :: t2(itloop)
      real*8 :: t3(itloop)
      real*8 :: t4(itloop)

      t1 = 0.0d0
      t2 = 0.0d0
      t3 = 0.0d0
      t4 = 0.0d0

      do 100 n = 1, numcl

        t4(n) = 1.d0 / tpp(idxcl(n)) - 1.d0/298.d0

  100 continue

      do 200 n = 1, numcl

        t1(n)  = e(1) * exp(f(1) * t4(n))
        t2(n)  = e(2) * exp(f(2) * t4(n))
        t3(n)  = e(3) * exp(f(3) * t4(n))

        hwp(n) = e(4) * exp(f(4) * t4(n))
        h3p(n) = e(5) * exp(f(5) * t4(n))

  200 continue

      do 300 n = 1, numcl

        hyp(n) = t1(n) * (1.d0 + t2(n) * rhpp * (1.d0 + t3(n) * rhpp))

  300 continue


      return

      end


!-----------------------------------------------------------------------------
!
! ROUTINE
!   chem23
!
! DESCRIPTION
!
!     *****************************************************************
!     * calculate the coefficients for the aqueous reaction equations
!     *****************************************************************
!
!     *    numcl : number of grid points in cloud
!     *    idxcl : array of the numbers of the numcl points in cloud
!     *    ohp, ho2p, h2op, no3p, o3p
!     *          : gas species concentrations (molecules cm-3)
!     *    hwp, hyp, h3p
!     *          : henry's law constant for H2O2, SO2, O3
!     *    cwp   : in-cloud cloud water mixing ratio (g/g)
!     *    zmp   : air density (#/cm3)
!     *    qh1p  : +----------------------
!     *    qh2p  : |
!     *    qh3p  : |
!     *    qj4p  : |
!     *    qh5p  : |  reaction rates
!     *    qh6   : |
!     *    qh7p  : |
!     *    qh8p  : |
!     *    qh9p  : +----------------------
!     *    a1-a9 : coefficients in aqueous equations
!
!-----------------------------------------------------------------------------

      subroutine chem23  &
     &   (itloop, numcl, idxcl,  &
     &    ohp, ho2p, h2op, no3p, o3p, hwp, hyp, h3p,  &
     &    cwp, zmp,  &
     &    qh1p, qh2p, qh3p, qj4p, qh5p, qh6, qh7p, qh8p, qh9p,  &
     &    a)

      implicit none

#     include "sulfchem.h"

!     ----------------------
!     GMIMOD-related values.
!     ----------------------

      integer :: itloop, numcl

      integer :: idxcl(itloop)

      real*8  :: qh6

      real*8  :: ohp (itloop)
      real*8  :: ho2p(itloop)
      real*8  :: h2op(itloop)
      real*8  :: no3p(itloop)
      real*8  :: o3p (itloop)
      real*8  :: hwp (itloop)
      real*8  :: hyp (itloop)
      real*8  :: h3p (itloop)
      real*8  :: cwp (itloop)
      real*8  :: zmp (itloop)
      real*8  :: qh1p(itloop)
      real*8  :: qh2p(itloop)
      real*8  :: qh3p(itloop)
      real*8  :: qj4p(itloop)
      real*8  :: qh5p(itloop)
      real*8  :: qh7p(itloop)
      real*8  :: qh8p(itloop)
      real*8  :: qh9p(itloop)
      real*8  :: a   (itloop, 9)

! -------------------------------------------------------------
!  *** LOCAL  VARIABLES WERE DECLARED HERE.
! -------------------------------------------------------------

      real*8  rlonum, d00

      real*8 :: t1(numcl)
      real*8 :: t2(numcl)
      real*8 :: t3(numcl)
      real*8 :: t4(numcl)
      real*8 :: t5(numcl)

      t1 = 0.0d0
      t2 = 0.0d0
      t3 = 0.0d0
      t4 = 0.0d0
      t5 = 0.0d0


      rlonum = 1. / alonum
      d00    = wtmair * 1.d-3 * rlonum

      do 100 n = 1, numcl

        i = idxcl(n)
        t1(n) = d00 * zmp(i) * cwp(i)

  100 continue

      do 200 n = 1, numcl

        t2(n)  =    1.d0  / (1.d0 + hyp(n) * t1(n))
        t3(n)  =    1.d0  / (1.d0 + hwp(n) * t1(n))
        t4(n)  = hyp(n)   / (1.d0 + hyp(n) * t1(n))
        t5(n)  = hwp(n)   / (1.d0 + hwp(n) * t1(n))

  200 continue

      do 300 n = 1, numcl

        i = idxcl(n)

        a(i, 1) = qh1p(i) * ohp(i)
        a(i, 2) = qh2p(i) * ohp(i) * t2(n)
        a(i, 3) = qh3p(i) * ho2p(i) * ho2p(i)
        a(i, 4) = qj4p(i) * t3(n)
        a(i, 5) = qh5p(i) * ohp(i) * t3(n)
        a(i, 6) = qh6     * t1(n) * t4(n) * t5(n) * rlonum
        a(i, 7) = qh7p(i) * t1(n) * t4(n) * h3p(n) * o3p(i) * rlonum
        a(i, 8) = qh8p(i) * ho2p(i) * ho2p(i) * h2op(i)
        a(i, 9) = qh9p(i) * no3p(i)

  300 continue


      return

      end

