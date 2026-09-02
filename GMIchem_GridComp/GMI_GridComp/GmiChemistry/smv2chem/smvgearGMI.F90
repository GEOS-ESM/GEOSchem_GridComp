
!=============================================================================
!
! $Id$
!
! CODE DEVELOPER
!   Original code from Mark Z. Jacobson ((C) COPYRIGHT, 1993).
!   LLNL modifications:  John Tannahill
!                        jrt@llnl.gov
!
! FILE
!   smvgear.F
!
! ROUTINES
!   Smvgear
!   Backsub
!   Decomp
!   Pderiv
!   Subfun
!   Update
!
!=============================================================================


!-----------------------------------------------------------------------------
!
! ROUTINE
!   Smvgear
!
! DESCRIPTION
!   This routine is the driver for the Smvgear (Sparse Matrix Vector Gear code)
!   chemistry solver.  It uses a Gear-type integrator that solves first order
!   ordinary differential equations with initial value boundary conditions.
!   Smvgear differs from an original Gear code in that it uses sparse matrix
!   and vectorization techniques to improve its computational speed.
!
!   This version is Smvgear II, 9/96.  It has been modified to include
!   grid-cell reordering prior to each time interval and different chemistry
!   for different atmospheric regions.  The purpose of the reordering is to
!   group cells with stiff equations together and those with non-stiff
!   equations together.  This reordering can save signifcant computer time
!   (e.g., speed the code by a factor of two or more), depending on the
!   variation in stiffness throughout the grid-domain.  When the stiffness is
!   the same throughout the grid-domain (e.g., if all concentrations and rates
!   are the same), then reordering is unnecessary and will not speed solutions.
!
!   This version includes a variable absolute error tolerance.  The absolute
!   tolerance is recalculated every few Gear time steps.  This version also
!   contains different sets of chemistry for different regions of the
!   atmosphere.  Thus, urban, free tropospheric, and stratospheric chemistry
!   can be solved during the same model run.
!
!   References =>
!   ----------
!
!     Jacobson M. Z. (1997) Improvement in Smvgear II through Absolute
!     Error Tolerance Control; in submission.
!
!     Jacobson M. Z. (1995) Computation of Global Photochemistry with Smvgear
!     II, Atmos. Environ., 29a, 2541-2546.
!
!     Jacobson M. Z. (1994) Developing, Coupling, and Applying a Gas, Aerosol,
!     Transport, and Radiation Model to Studying Urban and Regional Air
!     Pollution, PhD thesis, University of California, Los Angeles.
!
!     Jacobson M. Z. and Turco R. P. (1994) Smvgear: A Sparse Matrix,
!     Vectorized Gear Code for Atmospheric Models, Atmos. Environ. 28a,
!     273-284.
!
!     The origins of the Gear integrator used in Smvgear are found in:
!       Gear C. W. (1971) Numerical Initial Value Problems in Ordinary
!       Differential Equations, Prentice-Hall, NJ, pp. 158-166.
!
!     Finally, in subroutine Smvgear, the following ideas originated from
!     Lsodes, the Livermore solver for ordinary differential with sparse
!     matrices (Hindmarsh A. C. and Sherman A. H.):
!       (a) predicting the first time-step;
!       (b) determining corrector convergence differently than in Gear's
!           original code (goc);
!       (c) determining error differently than in goc;
!       (d) summing up the pascal matrix differently than in goc.
!
!     References for the 1987 Lsodes version include:
!
!       Sherman A. H. and Hindmarsh A. C. (1980) Gears: A Package for the
!       Solution of Sparse, Stiff Ordinary Differential Equations, Lawrence
!       Livermore National Laboratory Report UCRL-84102.
!
!       Hindmarsh A. C. (1983) Odepack, A Systematized Collection of ODE
!       Solvers, Scientific Computing, R.S. Stepleman et. al., eds.,
!       North-Holland, Amsterdam, pp. 55-74.
!
! ARGUMENTS
!   do_qqjk_inchem   : if pr_qqjk is on, should qqj's & qqk's be determined
!                      inside the chemistry solver, or outside?
!   do_semiss_inchem : do surface emissions inside the chemistry solver, or
!                      outside?
!   pr_qqjk  : should the periodic qqjk output file be written?
!   pr_smv2  : should the SmvgearII     output file be written
!              (non-parallel mode only)?
!   ifsun    : identifies whether sun is up (=1) or down (=2)
!   ilat     : # of latitudes
!   ilong    : # of longitudes
!   ivert    : # of vertical layers
!   ireord   : 1 => reorder grid-cells and blocks for chemistry
!              2 => solve chemistry
!   itloop   : # of zones (ilong * ilat * ivert)
!   jlooplo  : low ntloop grid-cell - 1 in a grid-block
!   ktloop   : # of grid-cells in a grid-block
!   lunsmv   : logical unit number to write to when pr_smv2 is true
!   nallr    : # of active rxns
!   ncs      : identifies gas chemistry type (1..NCSGAS)
!   nfdh2    : nfdh3 + # of rxns with two   active reactants
!   nfdh3    :         # of rxns with three active reactants
!   nfdl1    : nfdh2 + 1
!   nfdl2    : nfdh3 + 1
!   nfdrep   : nfdh3 + # of rxns with two active reactants that are not
!              followed by a rxn with the same reactants
!   nfdrep1  : nfdrep + 1
!   fracdec  : fraction time step is decreased in Smvgear if convergence
!              test fails
!   hmaxnit  : max time step for night all chem (s)
!   pr_nc_period     : NetCDF output period
!   tdt      : model time step (s)
!   do_cell_chem     : do chemistry for a particular cell?
!   irma,b,c : spc # of each reactant; locates reordered active spc #s
!   jreorder : gives original grid-cell from re-ordered grid-cell
!   jphotrat : tbd
!   ntspec   : # of active + inactive gases
!   inewold  : original spc # of each new jnew spc
!   denair   : density of air (molec/cm^3)
!   corig    : original gas-phase concs used to restart Smvgear if a failure
!              occurs (molec/cm^3)
!   pratk1   : tbd
!   yemis    : surface emissions (units?)
!   smvdm    : amount added to each spc at each grid-cell, for mass balance
!              accounting (# cm^-3 for gas chemistry (?))
!   nfdh1    : nfdh2 + # of rxns with one   active reactant
!   errmx2   : tbd
!   cc2      : array holding values of decomposed matrix
!   cnew     : stores conc (y (estimated)) (molec/cm^3)
!   gloss    : value of first derivatives on output from Subfun; right-side
!              of eqn on input to Backsub; error term (solution from Backsub)
!              on output from Backsub
!   vdiag    : 1 / current diagonal term of the decomposed matrix
!   rrate    : rate constants
!   trate    : rxn rate (moles l^-1-h2o s^-1 or # cm^-3 s^-1 (?))
!   urate    : term of Jacobian (J) = partial derivative
!
!-----------------------------------------------------------------------------

      subroutine SmvgearGMI  &
     &  (savedVars, do_qqjk_inchem, do_semiss_inchem, pr_qqjk, pr_smv2, ifsun,  &
     &   ilat, ilong, ivert, ireord, itloop, jlooplo, ktloop, lunsmv,  &
     &   nallr, ncs, nfdh2, nfdh3, nfdl1, nfdl2, nfdrep, nfdrep1,  &
     &   fracdec, hmaxnit, pr_nc_period, tdt, do_cell_chem, irma, irmb,  &
     &   irmc, jreorder, jphotrat, ntspec, inewold, denair, corig,  &
     &   pratk1, yemis, smvdm, nfdh1, errmx2, cc2, cnew, gloss, vdiag,  &
     &   rrate, trate, urate, &
     &   yda, qqkda, qqjda, qkgmi, qjgmi, &
     &   CTMi1, CTMi2, CTMju1, CTMj2, CTMk1, CTMk2, &
     &   num_qjo, num_qks, num_qjs, num_active, loc_proc)

      use GmiSavedVariables_mod, only : t_ChemistrySaved
      use GmiPrintError_mod, only : GmiPrintError

      implicit none

#     include "smv2chem_par.h"
!#     include "smv2chem2.h"


!     ----------------------
!     Argument declarations.
!     ----------------------

      integer, intent(in) :: CTMi1, CTMi2, CTMju1, CTMj2, CTMk1, CTMk2
      integer, intent(in) :: num_qjo, num_qks, num_qjs, num_active
      real*8 , intent(in   ) :: qjgmi(CTMi1:CTMi2, CTMju1:CTMj2, CTMk1:CTMk2, num_qjo)
      real*8 , intent(in   ) :: qkgmi(CTMi1:CTMi2, CTMju1:CTMj2, CTMk1:CTMk2, num_qks)
      real*8 , intent(inout) :: qqjda(CTMi1:CTMi2, CTMju1:CTMj2, CTMk1:CTMk2, num_qjs)
      real*8 , intent(inout) :: qqkda(CTMi1:CTMi2, CTMju1:CTMj2, CTMk1:CTMk2, num_qks)
      real*8 , intent(inout) :: yda  (CTMi1:CTMi2, CTMju1:CTMj2, CTMk1:CTMk2, num_active)
      logical, intent(in)  :: do_qqjk_inchem
      logical, intent(in)  :: do_semiss_inchem
      logical, intent(in)  :: pr_qqjk
      logical, intent(in)  :: pr_smv2
      integer, intent(in)  :: ifsun
      integer, intent(in)  :: ilat, ilong, ivert
      integer, intent(in)  :: ireord
      integer, intent(in)  :: itloop
      integer, intent(in)  :: jlooplo
      integer, intent(in)  :: ktloop
      integer, intent(in)  :: lunsmv
      integer, intent(in)  :: nallr
      integer, intent(in)  :: ncs
      integer, intent(in)  :: loc_proc
      integer, intent(in)  :: nfdh2,  nfdh3
      integer, intent(in)  :: nfdl1,  nfdl2
      integer, intent(in)  :: nfdrep, nfdrep1
      real*8,  intent(in)  :: fracdec
      real*8,  intent(in)  :: hmaxnit
      real*8,  intent(in)  :: pr_nc_period
      real*8,  intent(in)  :: tdt
      logical, intent(in)  :: do_cell_chem(ilong, ilat, ivert)
      integer, intent(in)  :: irma    (NMTRATE)
      integer, intent(in)  :: irmb    (NMTRATE)
      integer, intent(in)  :: irmc    (NMTRATE)
      integer, intent(in)  :: jreorder(itloop)
      integer, intent(in)  :: jphotrat(ICS)
      integer, intent(in)  :: ntspec  (ICS)
      integer, intent(in)  :: inewold (MXGSAER, ICS)
      real*8,  intent(in)  :: denair  (ktloop)
      real*8,  intent(in)  :: corig   (KBLOOP, MXGSAER)
      real*8,  intent(in)  :: pratk1  (KBLOOP, IPHOT)
      real*8,  intent(in)  :: yemis   (ilat*ilong, IGAS)

      real*8,  intent(inout) :: errmx2(itloop)
      real*8,  intent(inout) :: cc2   (KBLOOP, 0:MXARRAY)
      real*8,  intent(inout) :: cnew  (KBLOOP, MXGSAER)
      real*8,  intent(inout) :: gloss (KBLOOP, MXGSAER)
      real*8,  intent(inout) :: smvdm (KBLOOP, MXGSAER)
      real*8,  intent(inout) :: vdiag (KBLOOP, MXGSAER)
      real*8,  intent(inout) :: rrate (KBLOOP, NMTRATE)
      real*8,  intent(inout) :: trate (KBLOOP, NMTRATE*2)
      real*8,  intent(inout) :: urate (KBLOOP, NMTRATE, 3)

      integer, intent(out) :: nfdh1

      type(t_ChemistrySaved), intent(inOut) :: savedVars

!     ----------------------
!     Variable declarations.
!     ----------------------

!     ------------------------------------------------------------------------
!     idoub     : records # of steps since the last change in step size or
!                 order; it must be at least kstep = nqq+1 before doubling is
!                 allowed
!
!     ifail     : # of times corrector failed to converge while the Jacobian
!                 was old
!     lfail     : # of times accumulated error test failed
!     nfail     : # of times correcter failed to converge after Pderiv was
!                 called
!
!     ifsuccess : identifies whether step is successful (=1) or not (=0)
!     ischan    : # of first-order eqns to solve, = # of spc = order of
!                 original matrix; ischan has a different value for day and
!                 night and for gas- and aqueous-phase chemistry;
!                 # spc with prod or loss terms in Smvgear (?)
!     jeval     :  1 => call Pderiv the next time through the corrector steps;
!                  0 => last step successful and do not need to call Pderiv;
!                 -1 => Pderiv just called, and do not need to call again
!                  until jeval switched to 1
!     jrestar   : counts # of times Smvgear starts over at order 1 because of
!                 excessive failures
!     kstep     : nqq + 1
!
!     npderiv   : total # of times Pderiv is called
!     nsubfun   : total # of times Subfun is called
!
!     nqqisc    : nqq * ischan
!     nqqold    : value of nqq during last time step
!     nqq       : order of integration method; varies between 1 and MAXORD
!     nslp      : last time step # during which Pderiv was called
!     nsteps    : total # of successful time steps taken
!     ------------------------------------------------------------------------

      integer :: i, j, k
      integer :: i1, i2
      integer :: idoub
      integer :: ifail, jfail, lfail, nfail
      integer :: ifsuccess
      integer :: ischan, ischan1
      integer :: jb
      integer :: jeval
      integer :: jg1
      integer :: jgas
      integer :: jnew
      integer :: jrestar
      integer :: jspc
      integer :: k1, k2, k3, k4, k5
      integer :: kloop
      integer :: kstep, kstepisc
      integer :: l3
      integer :: nact
      integer :: ncsp  ! ncs       => for daytime   gas chemistry
                       ! ncs + ICS => for nighttime gas chemistry
      integer :: npderiv, nsubfun
      integer :: nqisc, nqqisc, nqqold
      integer :: nqq
      integer :: nslp
      integer :: nsteps
      integer :: nylowdec

      integer :: ibcb(IGAS)

!     ------------------------------------------------
!     kgrp : counts # of concs above abtol(i), i = 1..
!     ------------------------------------------------

      integer :: kgrp(KBLOOP, 5)

!     ------------------------------------------------------------------------
!     asn1      : value of aset(nqq,1)
!     delt      : current time step (s)
!     drate     : parameter which is used to determine whether convergence
!                 has occurred
!     edwn      : pertst^2*order for one order lower  than current order
!     enqq      : pertst^2*order for current order
!     eup       : pertst^2*order for one order higher than current order
!     hmax      : max time step at a given time (s)
!     hratio    : relative change in delt*aset(1) each change in step or order
!                 when Abs(hratio-1) > hrmax, reset jeval = 1 to call Pderiv
!     hrmax     : max relative change in delt*aset(1) before Pderiv is called
!     order     : floating point value of ischan, the order of # of ODEs
!     rdelmax   : max factor by which delt can be increased in a single step;
!                 as in Lsodes, set it to 1d4 initially to compensate for the
!                 small initial delt, but then set it to 10 after successful
!                 steps and to 2 after unsuccessful steps
!     rdelt     : factor (time step ratio) by which delt is increased or
!                 decreased
!     rdeltdn   : time step ratio at one order lower  than current order
!     rdeltsm   : time step ratio at current order
!     rdeltup   : time step ratio at one order higher than current order
!     rmsrat    : ratio of current to previous rms scaled error; if this
!                 ratio decreases, then convergence is occuring
!     timremain : remaining time in an chem interval (s)
!     tinterval : total chem time interval; same as chemintv (s)
!     told      : stores last value of xelaps in case current step fails
!     xelaps    : elapsed time in chem interval (s)
!     yfac      : = 1 originially, but is decreased if excessive failures
!                 occur in order to reduce absolute error tolerance
!     ------------------------------------------------------------------------

      real*8  :: abtoler1
      real*8  :: asn1, asnqqj
      real*8  :: cnewylow
      real*8  :: cnw
      real*8  :: conp1, conp2, conp3
      real*8  :: consmult
      real*8  :: dcon
      real*8  :: delt, delt1
      real*8  :: der1max, der2max, der3max
      real*8  :: drate
      real*8  :: dtasn1
      real*8  :: edwn, enqq, eup
      real*8  :: errinit, errinit_inv
      real*8  :: errmax_ncs_inv, errymax
      real*8  :: hmax, hrmax
      real*8  :: hmtim
      real*8  :: hratio
      real*8  :: iabove
      real*8  :: order, order_inv
      real*8  :: r1delt, rdelmax, rdelt, rdelta
      real*8  :: rdeltdn, rdeltsm, rdeltup
      real*8  :: real_kstep
      real*8  :: reltol1, reltol2, reltol3
      real*8  :: rmserr, rmserrp, rmsrat, rmstop
      real*8  :: timremain, tinterval
      real*8  :: told
      real*8  :: xelaps, xtimestep
      real*8  :: yfac

!     -------------------------------------------------------------------------
!     dely   : tbd
!     yabst  : absolute error tolerance (molec/cm^-3 for gases)
!
!     cest   : stores value of dtlos when idoub = 1
!     chold  : 1 / (reltol * cnew + abtol); multiply chold by local errors in
!              different error tests
!     dtlos  : an array of length ischan, used for the accumulated corrections;
!              on a successful return; dtlos(kloop,i) contains the estimated
!              one step local error in cnew
!     explic : tbd
!     conc   : an array of length ischan*(MAXORD+1) that carries the
!              derivatives of cnew, scaled by delt^j/factorial(j), where j is
!              the jth derivative; j varies from 1 to nqq; e.g., conc(jspc,2)
!              stores delt*y' (estimated)
!     -------------------------------------------------------------------------

      real*8  :: dely  (KBLOOP)
      real*8  :: yabst (KBLOOP)

      real*8  :: cest  (KBLOOP, MXGSAER)
      real*8  :: chold (KBLOOP, MXGSAER)
      real*8  :: dtlos (KBLOOP, MXGSAER)
      real*8  :: explic(KBLOOP, MXGSAER)

      real*8  :: conc  (KBLOOP, MXGSAER*7)


!     ----------------
!     Begin execution.
!     ----------------

!c    Write (6,*) 'Smvgear called.'


      nact = nnact

!     =======================
#     include "setkin_ibcb.h"
!     =======================

      ifail     = 0
      jfail     = 0
      lfail     = 0
      nfail     = 0
      npderiv   = 0
      nsteps    = 0
      nsubfun   = 0
      nylowdec  = 0

      hrmax     = 0.3d0
      rmserr    = 1.0d0

      ischan    = savedVars%ischang(ncs)
      ischan1   = ischan - 1
      order     = ischan
      order_inv = 1.0d0 / ischan

      tinterval = savedVars%timeintv(ncs)

      ncsp      = (ifsun - 1) * ICS + ncs

      if (ifsun == 1) then
        hmax = savedVars%hmaxday(ncs)
      else
        hmax = hmaxnit
      end if


      yfac   = 1.0d0
      iabove = order * 0.4d0

      errinit     = Min (savedVars%errmax(ncs), 1.0d-03)
      errinit_inv = 1.0d0 / errinit

      errmax_ncs_inv = 1.0d0 / savedVars%errmax(ncs)


!     ----------------------------------------------------
!     Start time interval or re-enter after total failure.
!     ----------------------------------------------------

!     ========
 100  continue
!     ========

      idoub     = 2
      nslp      = MBETWEEN
      jrestar   = 0
      xelaps    = 0.0d0
      told      = 0.0d0
      timremain = tinterval

      reltol1   = yfac * errinit_inv

      reltol2   = yfac * errmax_ncs_inv

      reltol3   = errmax_ncs_inv

      abtoler1  = savedVars%abtol(6,ncs) * reltol1

!     -------------------------------
!     Initialize concentration array.
!     -------------------------------

      do jnew = 1, ischan
        do kloop = 1, ktloop
          cnew(kloop, jnew) = corig(kloop, jnew)
        end do
      end do

!     --------------------------------------------------------------------
!     Re-enter here if total failure or if restarting with new cell block.
!     --------------------------------------------------------------------

!     ========
 150  continue
!     ========

      hratio    = 0.0d0
      asn1      = 1.0d0
      ifsuccess = 1
      rdelmax   = 1.0d4

!     ---------------------
!     Initialize photrates.
!     ---------------------

!DIR$ INLINE
!     ===========
      call UpdateGMI  &
!     ===========
     &  (savedVars, ktloop, nallr, ncs, ncsp, jphotrat, pratk1, rrate, trate)
!DIR$ NOINLINE


!     ----------------------------
!     Initialize first derivative.
!     ----------------------------

!     ===========
      call SubfunGMI  &
!     ===========
     &  (savedVars, ischan, ktloop, nallr, ncsp, nfdh2, nfdh3, nfdl1, nfdl2,  &
     &   nfdrep, nfdrep1, irma, irmb, irmc, cnew, rrate, nsubfun,  &
     &   gloss, trate, nfdh1)


!     ----------------------------------------------------------------
!     Zero first derviatives in surface zones for species with fixed
!     concentration boundary conditions (LLNL addition, PSC, 5/14/99).
!     ----------------------------------------------------------------

      do kloop = 1, ktloop

        if (jreorder(jlooplo+kloop) <= (ilat*ilong)) then

          do jspc = 1, ntspec(ncs)

            jgas = inewold(jspc,1)

            if (ibcb(jgas) == 1) then
              gloss(kloop,jspc) = 0.0d0
            else if (do_semiss_inchem) then
              gloss(kloop,jspc) =  &
     &          gloss(kloop,jspc) +  &
     &          yemis(jreorder(jlooplo+kloop),jgas)
            end if

          end do

        end if

      end do

!     -------------------------------------------
!     Determine initial absolute error tolerance.
!     -------------------------------------------

      do kloop = 1, ktloop
        dely(kloop) = 0.0d0
      end do


!     ==========================
      IREORDIF: if (ireord /= 1) then
!     ==========================

        do k = 1, 5
          do kloop = 1, ktloop
            kgrp(kloop,k) = 0
          end do
        end do

        do jspc = 1, ischan
          do kloop = 1, ktloop

            cnw = cnew(kloop,jspc)

            if (cnw > savedVars%abtol(1,ncs)) then
              kgrp(kloop,1) = kgrp(kloop,1) + 1
            else if (cnw > savedVars%abtol(2,ncs)) then
              kgrp(kloop,2) = kgrp(kloop,2) + 1
            else if (cnw > savedVars%abtol(3,ncs)) then
              kgrp(kloop,3) = kgrp(kloop,3) + 1
            else if (cnw > savedVars%abtol(4,ncs)) then
              kgrp(kloop,4) = kgrp(kloop,4) + 1
            else if (cnw > savedVars%abtol(5,ncs)) then
              kgrp(kloop,5) = kgrp(kloop,5) + 1
            end if

          end do
        end do

        do kloop = 1, ktloop

          k1 = kgrp(kloop,1)
          k2 = kgrp(kloop,2) + k1
          k3 = kgrp(kloop,3) + k2
          k4 = kgrp(kloop,4) + k3
          k5 = kgrp(kloop,5) + k4

          if (k1 > iabove) then
            yabst(kloop) = savedVars%abtol(1,ncs)
          else if (k2 > iabove) then
            yabst(kloop) = savedVars%abtol(2,ncs)
          else if (k3 > iabove) then
            yabst(kloop) = savedVars%abtol(3,ncs)
          else if (k4 > iabove) then
            yabst(kloop) = savedVars%abtol(4,ncs)
          else if (k5 > iabove) then
            yabst(kloop) = savedVars%abtol(5,ncs)
          else
            yabst(kloop) = savedVars%abtol(6,ncs)
          end if

        end do

!c
        do kloop = 1, ktloop
          do jspc = 1, ischan
            cnewylow    = cnew (kloop,jspc) + (yabst(kloop) * reltol1)
            errymax     = gloss(kloop,jspc) / cnewylow
            dely(kloop) = dely (kloop) + (errymax * errymax)
          end do
        end do

!     ====
      else
!     ====

!       ------------------------------------------------------
!       Use lowest absolute error tolerance when reordering.
!       If reordering, set errmx2 then return to Physproc.
!
!       abtoler1 = yfac * abtol(6,ncs) / Min (errmax, 1.0d-03)
!       ------------------------------------------------------

!c
        do kloop = 1, ktloop
          do jspc = 1, ischan
            errymax     = gloss(kloop,jspc) /  &
     &                    (cnew(kloop,jspc) + abtoler1)
            dely(kloop) = dely(kloop) + (errymax * errymax)
          end do
        end do

        do kloop = 1, ktloop
          errmx2(jlooplo+kloop) = dely(kloop)
        end do

!       =========
        go to 700
!       =========

!     ===============
      end if IREORDIF
!     ===============


!     --------------------------------------
!     Calculate initial time step size (s).
!
!     Sqrt (dely / [errinit * order]) =
!       rmsnorm of error scaled to errinit *
!       cnew + abtol / reltol
!     --------------------------------------

      rmstop = 0.0d0

      do kloop = 1, ktloop
        if (dely(kloop) > rmstop) then
          rmstop = dely(kloop)
        end if
      end do

      delt1 = Sqrt (errinit / (savedVars%abst2(ncs) + (rmstop * order_inv)))
      delt  = Max  (Min (delt1, timremain, hmax), HMIN)


!     -----------------------
!     Set initial order to 1.
!     -----------------------

      nqqold = 0
      nqq    = 1
      jeval  = 1
      rdelt  = 1.0d0


!     --------------------------------------------------------------
!     Store initial concentration and first derivatives x time step.
!     --------------------------------------------------------------

      do jspc = 1, ischan

        j = jspc + ischan

        do kloop = 1, ktloop
          conc(kloop,jspc) = cnew(kloop,jspc)
          conc(kloop,j)    = delt * gloss(kloop,jspc)
        end do

      end do


!     ========
 200  continue
!     ========


!     -------------------------------------------------------------------
!     Update coefficients of the order; note that pertst2 is the original
!     pertst^2.
!     -------------------------------------------------------------------

      if (nqq /= nqqold) then

        nqqold = nqq
        kstep  = nqq + 1
        hratio = hratio * savedVars%aset(nqq,1) / asn1
        asn1   = savedVars%aset(nqq,1)
        enqq   = savedVars%pertst2(nqq,1) * order
        eup    = savedVars%pertst2(nqq,2) * order
        edwn   = savedVars%pertst2(nqq,3) * order
        conp3  = 1.4d0 /  (eup**savedVars%enqq3(nqq))
        conp2  = 1.2d0 / (enqq**savedVars%enqq2(nqq))
        conp1  = 1.3d0 / (edwn**savedVars%enqq1(nqq))
        nqqisc = nqq * ischan

      end if


!     ----------------------------------------------------------------
!     Limit size of rdelt, then recalculate new time step and update
!     hratio.  Use hratio to determine whether Pderiv should be called
!     again.
!     ----------------------------------------------------------------

      hmtim  = Min (hmax, timremain)
      rdelt  = Min (rdelt, rdelmax, hmtim/delt)
      delt   = delt   * rdelt
      hratio = hratio * rdelt
      xelaps = xelaps + delt

      if ((Abs (hratio-1.0d0) > hrmax) .or. (nsteps >= nslp)) then
        jeval = 1
      end if


!     ----------------------------------------------------------
!     If time step < HMIN, tighten absoloute error tolerance and
!     restart integration at beginning of time interval.
!     ----------------------------------------------------------

      if (delt < HMIN) then

        if (pr_smv2) then
          Write (lunsmv,950) delt, timremain, yfac, savedVars%errmax(ncs)
        end if

 950    format ('Smvgear:  delt      = ', 1pe9.3, /,  &
     &          '          timremain = ', 1pe9.3, /,  &
     &          '          yfac      = ', 1pe9.3, /,  &
     &          '          errmax    = ', 1pe9.3)

        nylowdec = nylowdec + 1
        yfac     = yfac * 0.01d0

        if (nylowdec == 10) then

          if (pr_smv2) then
            Write (lunsmv,960)
          end if

 960      format ('Smvgear:  too many decreases of yfac.')

        call GmiPrintError ('Problem in Smvgear', .true., 0, 0, 0, 0, 0.0d0, 0.0d0)

        end if

!       =========
        go to 100
!       =========

      end if


!     -------------------------------------------------------------------
!     If the delt is different than during the last step (if rdelt /= 1),
!     then scale the derivatives.
!     -------------------------------------------------------------------

      if (rdelt /= 1.0d0) then

        rdelta = 1.0d0
        i1     = 1

        do j = 2, kstep

          rdelta = rdelta * rdelt
          i1     = i1 + ischan

          do i = i1, i1 + ischan1
            do kloop = 1, ktloop
              conc(kloop,i) = conc(kloop,i) * rdelta
            end do
          end do

        end do

      end if


!     --------------------------------------------------------------
!     If the last step was successful, reset rdelmax = 10 and update
!     the chold array with current values of cnew.
!     --------------------------------------------------------------

!     ================================
      IFSUCCESSIF: if (ifsuccess == 1) then
!     ================================

        rdelmax = 10.0d0

!       ---------------------------------------
!       Determine new absolute error tolerance.
!       ---------------------------------------

        if (Mod (nsteps, 3) == 2) then

          do k = 1, 5
            do kloop = 1, ktloop
              kgrp(kloop,k) = 0
            end do
          end do

          do jspc = 1, ischan
            do kloop = 1, ktloop

              cnw = cnew(kloop,jspc)

              if (cnw > savedVars%abtol(1,ncs)) then
                kgrp(kloop,1) = kgrp(kloop,1) + 1
              else if (cnw > savedVars%abtol(2,ncs)) then
                kgrp(kloop,2) = kgrp(kloop,2) + 1
              else if (cnw > savedVars%abtol(3,ncs)) then
                kgrp(kloop,3) = kgrp(kloop,3) + 1
              else if (cnw > savedVars%abtol(4,ncs)) then
                kgrp(kloop,4) = kgrp(kloop,4) + 1
              else if (cnw > savedVars%abtol(5,ncs)) then
                kgrp(kloop,5) = kgrp(kloop,5) + 1
              end if

            end do
          end do

          do kloop = 1, ktloop

            k1 = kgrp(kloop,1)
            k2 = kgrp(kloop,2) + k1
            k3 = kgrp(kloop,3) + k2
            k4 = kgrp(kloop,4) + k3
            k5 = kgrp(kloop,5) + k4

            if (k1 > iabove) then
              yabst(kloop) = savedVars%abtol(1,ncs)
            else if (k2 > iabove) then
              yabst(kloop) = savedVars%abtol(2,ncs)
            else if (k3 > iabove) then
              yabst(kloop) = savedVars%abtol(3,ncs)
            else if (k4 > iabove) then
              yabst(kloop) = savedVars%abtol(4,ncs)
            else if (k5 > iabove) then
              yabst(kloop) = savedVars%abtol(5,ncs)
            else
              yabst(kloop) = savedVars%abtol(6,ncs)
            end if

          end do

        end if

!c
        do kloop = 1, ktloop
          do jspc = 1, ischan

            chold(kloop,jspc) =  &
     &        reltol3 /  &
     &        (Max (cnew(kloop,jspc), 0.0d0) +  &
     &         (yabst(kloop) * reltol2))

          end do
        end do

!     ==================
      end if IFSUCCESSIF
!     ==================


!     ------------------------------------------------------------------
!     Compute the predicted concentration and derivatives by multiplying
!     previous values by the pascal triangle matrix.
!     ------------------------------------------------------------------

      i1 = nqqisc + 1

      do jb = 1, nqq - 1

        i1 = i1 - ischan

        do i = i1,  nqqisc

          j = i + ischan

          do kloop = 1, ktloop
            conc(kloop,i)  = conc(kloop,i) + conc(kloop,j)
          end do

        end do

      end do

      do jspc = 1,  ischan

        j = jspc + ischan

        do kloop = 1, ktloop
          conc  (kloop,jspc) = conc(kloop,jspc) + conc(kloop,j)
          explic(kloop,jspc) = conc(kloop,j)
        end do

      end do

      do i = ischan + 1, nqqisc

        j = i + ischan

        do kloop = 1, ktloop
          conc(kloop,i) = conc(kloop,i) + conc(kloop,j)
        end do

      end do


!     -------------------------------------------------------------------
!     Correction loop.
!
!     Take up to 3 corrector iterations.  Test convergence by requiring
!     that changes be less than the rms norm weighted by chold.
!     Accumulate the correction in the array dtlos.  It equals the
!     jth derivative of concentration multiplied by delt^kstep /
!     (factorial(kstep-1) * aset(kstep)); thus, it is proportional to the
!     actual errors to the lowest power of delt present (delt^kstep).
!     -------------------------------------------------------------------


!     ========
 250  continue
!     ========


      l3 = 0

      do jspc = 1, ischan
        do kloop = 1, ktloop
          cnew (kloop,jspc) = conc(kloop,jspc)
          dtlos(kloop,jspc) = 0.0d0
        end do
      end do


!     ------------------------------------------------------------------
!     If jeval = 1, re-evaluate predictor matrix P = I - H * aset(1) * J
!     before starting the corrector iteration.  After calling Pderiv,
!     set jeval = -1 to prevent recalling Pderiv unless necessary later.
!     Call Decomp to decompose the matrix.
!     ------------------------------------------------------------------

      if (jeval == 1) then

        r1delt = -asn1 * delt

!DIR$   INLINE
!       ===========
        call PderivGMI  &
!       ===========
     &    (savedVars, ischan, ktloop, ncsp, nfdh2, nfdh3, nfdl1, nfdl2, irma, irmb,  &
     &     irmc, r1delt, cnew, rrate, npderiv, cc2, urate, nfdh1)
!DIR$   NOINLINE

!DIR$   INLINE
!       ===========
        call DecompGMI  &
!       ===========
     &    (savedVars, ischan, ktloop, ncsp, cc2, vdiag)
!DIR$   NOINLINE

        jeval  = -1
        hratio = 1.0d0
        nslp   = nsteps + MBETWEEN
        drate  = 0.7d0

      end if


!     -------------------------------------------------------------
!     Evaluate the first derivative using corrected values of cnew.
!     -------------------------------------------------------------

!     ========
 300  continue
!     ========

!     ===========
      call SubfunGMI  &
!     ===========
     &  (savedVars, ischan, ktloop, nallr, ncsp, nfdh2, nfdh3, nfdl1, nfdl2,  &
     &   nfdrep, nfdrep1, irma, irmb, irmc, cnew, rrate, nsubfun,  &
     &   gloss, trate, nfdh1)


!     ----------------------------------------------------------------
!     Zero first derviatives in surface zones for species with fixed
!     concentration boundary conditions.  Include surface emissions in
!     first derviatives.
!
!     Species with non-zero fluxes (LLNL addition, PSC, 7/24/96).
!     ----------------------------------------------------------------

      do kloop = 1, ktloop

        if (jreorder(jlooplo+kloop) <= (ilat*ilong)) then

          do jspc = 1, ntspec(ncs)

            jgas = inewold(jspc,1)

            if (ibcb(jgas) == 1) then

              gloss(kloop,jspc) = 0.0d0

            else

              if (do_semiss_inchem) then
                gloss(kloop,jspc) = gloss(kloop,jspc) +  &
     &                              yemis(jreorder(jlooplo+kloop),jgas)
              end if

            end if

          end do

        end if

      end do


!     ---------------------------------------------------------------
!     In the case of the chord method, compute error (gloss) from the
!     corrected calculation of the first derivative.
!     ---------------------------------------------------------------

      do jspc = 1, ischan

        j = jspc + ischan

        do kloop = 1, ktloop
          gloss(kloop,jspc) = (delt * gloss(kloop,jspc)) -  &
     &                        (conc(kloop,j) + dtlos(kloop,jspc))
        end do

      end do


!     --------------------------------------------------------------
!     Solve the linear system of equations with the corrector error;
!     Backsub solves backsubstitution over matrix of partial derivs.
!     --------------------------------------------------------------

!     ============
      call BacksubGMI  &
!     ============
     &  (savedVars, ischan, ktloop, ncsp, cc2, vdiag, gloss)


!     ----------------------------------------------------------------
!     Sum up the accumulated error, correct the concentration with the
!     error, and begin to calculate the rmsnorm of the error relative
!     to chold.
!     ----------------------------------------------------------------

      do kloop = 1, ktloop
        dely(kloop) = 0.0d0
      end do

      if (asn1 == 1.0d0) then

        do i = 1, ischan
          do kloop = 1, ktloop

            dtlos(kloop,i) = dtlos(kloop,i) + gloss(kloop,i)
            cnew(kloop,i)  = conc(kloop,i)  + dtlos(kloop,i)
            errymax        = gloss(kloop,i) * chold(kloop,i)
            dely(kloop)    = dely(kloop)    + (errymax * errymax)

          end do
        end do

      else

        do i = 1, ischan
          do kloop = 1, ktloop

            dtlos(kloop,i) = dtlos(kloop,i) + gloss(kloop,i)
            cnew(kloop,i)  = conc(kloop,i)  + (asn1 * dtlos(kloop,i))
            errymax        = gloss(kloop,i) * chold(kloop,i)
            dely(kloop)    = dely(kloop)    + (errymax * errymax)

          end do
        end do

      end if


!     ------------------------------------------------------------------
!     Set the previous rms error and calculate the new rms error.
!     If dcon < 1, then sufficient convergence has occurred.  Otherwise,
!     if the ratio of the current to previous rmserr is decreasing,
!     iterate more.  If it is not, then the convergence test failed.
!     ------------------------------------------------------------------

      rmserrp = rmserr
      der2max = 0.0d0


      do kloop = 1, ktloop

        if (dely(kloop) > der2max) then
          der2max = dely(kloop)
        end if

      end do


      rmserr = Sqrt (der2max * order_inv)


      l3 = l3 + 1

      if (l3 > 1) then
        rmsrat = rmserr / rmserrp
        drate  = Max (0.2d0*drate, rmsrat)
      else
        rmsrat = 1.0d0
      end if

      dcon = rmserr * Min (savedVars%conpst(nqq), savedVars%conp15(nqq)*drate)


!     --------------------------------------------------------
!     If convergence occurs, go on to check accumulated error.
!     --------------------------------------------------------

      if (dcon > 1.0d0) then

!       -------------------------------------------------------------------
!       If nonconvergence after one step, re-evaluate first derivative with
!       new values of cnew.
!       -------------------------------------------------------------------

        if (l3 == 1) then

!         =========
          go to 300
!         =========

!         ----------------------------------------------------------------
!         The corrector iteration failed to converge.
!
!         If the Jacobian matrix is more than one step old, update the
!         Jacobian and try convergence again.  If the Jacobian is current,
!         then reduce the time step, reset the accumulated derivatives to
!         their values before the failed step, and retry with the smaller
!         step.
!         ----------------------------------------------------------------

        else if (jeval == 0) then

          ifail = ifail + 1
          jeval = 1

!         =========
          go to 250
!         =========

        end if

        nfail     = nfail + 1
        rdelmax   = 2.0d0
        jeval     = 1
        ifsuccess = 0
        xelaps    = told
        rdelt     = fracdec

        i1 = nqqisc + 1

        do jb = 1, nqq

          i1 = i1 - ischan

          do i = i1, nqqisc

            j = i + ischan

            do kloop = 1, ktloop
              conc(kloop,i) = conc(kloop,i) - conc(kloop,j)
            end do

          end do

        end do

!       =========
        go to 200
!       =========

      end if


!     -------------------------------------------------------------------
!     The corrector iteration converged.
!
!     Set jeval = 0, so that it does not need to be called the next step.
!     If all else goes well.  Next, test the accumulated error from the
!     convergence process above.
!     -------------------------------------------------------------------

      jeval = 0

      if (l3 > 1) then

        do kloop = 1, ktloop
          dely(kloop) = 0.0d0
        end do

!c
        do kloop = 1, ktloop
          do jspc = 1, ischan
            errymax     = dtlos(kloop,jspc) * chold(kloop,jspc)
            dely(kloop) = dely(kloop) + errymax * errymax
          end do
        end do

        der2max = 0.0d0

        do kloop = 1, ktloop

          if (dely(kloop) > der2max) then
            der2max = dely(kloop)
          end if

        end do

      end if


!     ----------------------------------------------------------------
!     The accumulated error test failed.
!
!     In all cases, reset the derivatives to their values before the
!     last time step.  Next:
!       (a) re-estimate a time step at the same or one lower order and
!           retry the step;
!       (b) if the first attempts fail, retry the step at fracdec the
!           the prior step;
!       (c) iF this fails, reset the order to 1 and go back to the
!           beginning, at order = 1, because errors of the wrong order
!           have accumulated.
!     ----------------------------------------------------------------

!     ==============================
      DER2MAXIF: if (der2max > enqq) then
!     ==============================

        xelaps = told
        lfail  = lfail + 1
        jfail  = jfail  + 1
        i1     = nqqisc + 1

        do jb = 1, nqq

          i1 = i1 - ischan

          do i = i1, nqqisc

            j = i + ischan

            do kloop = 1, ktloop
              conc(kloop,i) = conc(kloop,i) - conc(kloop,j)
            end do

          end do

        end do

        rdelmax = 2.0d0

        if (jfail <= 6) then

          ifsuccess = 0
          rdeltup   = 0.0d0

!         =========
          go to 400
!         =========

        else if (jfail <= 20) then

          ifsuccess = 0
          rdelt     = fracdec

!         =========
          go to 200
!         =========

        else

          delt    = delt * 0.1d0
          rdelt   = 1.0d0
          jfail   = 0
          jrestar = jrestar + 1
          idoub   = 5

          do jspc = 1, ischan
            do kloop = 1, ktloop
              cnew(kloop,jspc) = conc(kloop,jspc)
            end do
          end do

          if (pr_smv2) then
            Write (lunsmv,970) delt, xelaps
          end if

 970      format ('delt dec to ', e13.5, ' at time ', e13.5,  &
     &            ' because of excessive errors.')

          if (jrestar == 100) then

            if (pr_smv2) then
              Write (lunsmv,980)
            end if

 980        format ('Smvgear:  Stopping because of excessive errors.')

            call GmiPrintError ('Problem in Smvgear', .true., 0, 0, 0, 0, 0.0d0, 0.0d0)

          end if

!         =========
          go to 150
!         =========

        end if

!     ====
      else
!     ====

!       -------------------------------------------------------------
!       All successful steps come through here.
!
!       After a successful step, update the concentration and all
!       derivatives, reset told, set ifsuccess = 1, increment nsteps,
!       and reset jfail = 0.
!       -------------------------------------------------------------


        if (pr_qqjk .and. do_qqjk_inchem) then
          xtimestep = xelaps - told

!         =================
          call Do_Smv2_Diag  &
!         =================
     &      (jlooplo, ktloop, pr_nc_period, tdt, told, do_cell_chem,  &
     &       jreorder, inewold, denair, cnew, xtimestep, &
     &       yda, qqkda, qqjda, qkgmi, qjgmi, &
     &       ilong, ilat, ivert, itloop, &
     &       CTMi1, CTMi2, CTMju1, CTMj2, CTMk1, CTMk2, &
     &       num_qjo, num_qks, num_qjs, num_active)
        end if


        jfail     = 0
        ifsuccess = 1
        nsteps    = nsteps + 1
        told      = xelaps

        i1 = 1

        do j = 2, kstep

          i1 = i1 + ischan

          asnqqj = savedVars%aset(nqq,j)

          do jspc = 1, ischan

            i = jspc + i1 - 1

            do kloop = 1, ktloop
              conc(kloop,i) =  &
     &          conc(kloop,i) + (asnqqj * dtlos(kloop,jspc))
            end do

          end do

        end do

!       ------------------------------
!       Update chemistry mass balance.
!       ------------------------------

        if (asn1 == 1.0d0) then

          do i = 1, ischan
            do kloop = 1, ktloop

              smvdm(kloop,i) =  &
     &          smvdm(kloop,i) + dtlos(kloop,i) + explic(kloop,i)

              conc(kloop,i) = conc(kloop,i) + dtlos(kloop,i)

            end do
          end do

        else

          do i = 1, ischan
            do kloop = 1, ktloop

              dtasn1         = asn1 * dtlos(kloop,i)
              smvdm(kloop,i) = smvdm(kloop,i) + dtasn1 + explic(kloop,i)
              conc (kloop,i) = conc (kloop,i) + dtasn1

            end do
          end do

        end if

!       ---------------------------------------------------
!       Exit smvgear if a time interval has been completed.
!       ---------------------------------------------------

        timremain = tinterval - xelaps

        if (timremain <= 1.0d-06) then
!         =========
          go to 700
!         =========
        end if

!       -------------------------------------------------------------------
!       idoub counts the number of successful steps before re-testing the
!       step-size and order:
!         if idoub > 1, decrease idoub and go on to the next time step with
!                       the current step-size and order;
!         if idoub = 1, store the value of the error (dtlos) for the time
!                       step prediction, which will occur when idoub = 0,
!                       but go on to the next step with the current step
!                       size and order;
!         if idoub = 0, test the time step and order for a change.
!       -------------------------------------------------------------------

        if (idoub > 1) then

          idoub = idoub - 1

          if (idoub == 1) then

            do jspc = 1, ischan, 2

              jg1 = jspc + 1

              do kloop = 1, ktloop
                cest(kloop,jspc) = dtlos(kloop,jspc)
                cest(kloop,jg1)  = dtlos(kloop,jg1)
              end do

            end do

          end if

          rdelt = 1.0d0

!         =========
          go to 200
!         =========

        end if

!     ================
      end if DER2MAXIF
!     ================


!     ------------------------------------------------------------------
!     Test whether to change the step-size and order.
!
!     Determine the time step at (a) one order lower than, (b) the same
!     order as, and (c) one order higher than the current order.  In the
!     case of multiple grid-cells in a grid-block, find the minimum
!     step size among all the cells for each of the orders.  Then, in
!     all cases, choose the longest time step among the three steps
!     paired with orders, and choose the order allowing this longest
!     step.
!     ------------------------------------------------------------------

!     ---------------------------------------------------------------
!     Estimate the time step ratio (rdeltup) at one order higher than
!     the current order.  If nqq >= MAXORD, then we do not allow the
!     order to increase.
!     ---------------------------------------------------------------

      if (nqq < MAXORD) then

        do kloop = 1, ktloop
          dely(kloop) = 0.0d0
        end do

        do jspc = 1, ischan
          do kloop = 1, ktloop
            errymax     = (dtlos(kloop,jspc) - cest(kloop,jspc)) *  &
     &                    chold(kloop,jspc)
            dely(kloop) = dely(kloop) + (errymax * errymax)
          end do
        end do

        der3max = 0.0d0

        do kloop = 1, ktloop

          if (dely(kloop) > der3max) then
            der3max = dely(kloop)
          end if

        end do

        rdeltup = 1.0d0 / ((conp3 * der3max**savedVars%enqq3(nqq)) + 1.4d-6)

      else

        rdeltup = 0.0d0

      end if


!     ========
 400  continue
!     ========


!     ------------------------------------------------------------
!     Estimate the time step ratio (rdeltsm) at the current order.
!     der2max was calculated during the error tests earlier.
!     ------------------------------------------------------------

      rdeltsm = 1.0d0 / ((conp2 * der2max**savedVars%enqq2(nqq)) + 1.2d-6)


!     ------------------------------------------------------------------
!     Estimate the time step ratio (rdeltdn) at one order lower than
!     the current order.  if nqq = 1, then we cannot test a lower order.
!     ------------------------------------------------------------------

      if (nqq > 1) then

        do kloop = 1, ktloop
          dely(kloop) = 0.0d0
        end do

        kstepisc = (kstep - 1) * ischan

!c
        do kloop = 1, ktloop
          do jspc = 1, ischan

            i = jspc + kstepisc

            errymax     = conc(kloop,i) * chold(kloop,jspc)
            dely(kloop) = dely(kloop) + (errymax * errymax)

          end do

        end do

        der1max = 0.0d0

        do kloop = 1, ktloop
          if (dely(kloop) > der1max) then
            der1max = dely(kloop)
          end if
        end do

        rdeltdn = 1.0d0 / ((conp1 * der1max**savedVars%enqq1(nqq)) + 1.3d-6)

      else

        rdeltdn = 0.0d0

      end if


!     -----------------------------------------------------------------
!     Find the largest of the predicted time step ratios of each order.
!     -----------------------------------------------------------------

      rdelt = Max (rdeltup, rdeltsm, rdeltdn)


!     ---------------------------------------------------------------
!     If the last step was successful and rdelt is small, keep the
!     current step and order, and allow three successful steps before
!     re-checking the time step and order.
!     ---------------------------------------------------------------

      if ((rdelt < 1.1d0) .and. (ifsuccess == 1)) then

        idoub = 3

!       =========
        go to 200
!       =========

!       --------------------------------------------------------------
!       If the maximum time step ratio is that of one order lower than
!       the current order, decrease the order.  Do not minimize rdelt
!       to <= 1, when ifsuccess = 0 since this is less efficient.
!       --------------------------------------------------------------

      else if (rdelt == rdeltdn) then

        nqq = nqq - 1

!       ---------------------------------------------------------------
!       If the maximum time step ratio is that of one order higher than
!       the current order, increase the order and add a derivative term
!       for the higher order.
!       ---------------------------------------------------------------

      else if (rdelt == rdeltup) then

        real_kstep = kstep
        consmult   = savedVars%aset(nqq,kstep) / real_kstep
        nqq        = kstep
        nqisc      = nqq * ischan

        do jspc = 1, ischan, 2

          jg1 = jspc + 1
          i1  = jspc + nqisc
          i2  = jg1  + nqisc

          do kloop = 1, ktloop
            conc(kloop,i1) = dtlos(kloop,jspc) * consmult
            conc(kloop,i2) = dtlos(kloop,jg1)  * consmult
          end do

        end do

      end if

!     ----------------------------------------------------------------
!     If the last two steps have failed, re-set idoub to the current
!     order + 1.  Do not minimize rdelt if jfail >= 2 since tests show
!     that this merely leads to additional computations.
!     ----------------------------------------------------------------

      idoub = nqq + 1

!     =========
      go to 200
!     =========

!     ========
 700  continue
!     ========

      !PRINT*, "@@@ Jules @@@ ", nsubfun, " by: ", loc_proc

      return

      end subroutine SmvgearGMI


!-----------------------------------------------------------------------------
!
! ROUTINE
!   Backsub
!
! DESCRIPTION
!   This routine performs back-substitutions on the decomposed matrix.  It
!   solves the linear set of equations Ax = B FOR x, the correction vector,
!   where "A" is the L-U decompostion of the original matrix =>
!
!     P = I - H x Bo x J
!
!   I = identity matrix, H = time step, Bo = a coefficient that depends on
!   the order of the integration method, and J is the matrix of partial
!   derivatives.  B is sent from Smvgear as a corrected value of the first
!   derivatives of the ordinary differential equations.  Decomp solved for
!   "A", the decomposed matrix.  See Press, et. al. (1992), Numerical
!   Recipes, Cambridge University Press, for a better description of the
!   back-substitution process.
!
!   This back-substitution process uses sparse matrix techniques,
!   vectorizes around the grid-cell dimension, and uses no partial
!   pivoting.  Tests by Sherman & Hindmarsh (1980), Lawrence Livermore
!   Livermore Laboratory, Rep. UCRL-84102, and by us have confirmed that
!   the removal of partial pivoting has little effect on results.
!
!   Backsub loop # 1 =>
!     First, adjust right side of Ax = B using lower triangular matrix.
!     Sum 1,2,3,4, or 5 terms at a time to improve vectorization.
!
! ARGUMENTS
!   ischan : # of first-order eqns to solve, = # of spc = order of original
!            matrix; ischan has a different value for day and night, and for
!            gas- and aqueous-phase chemistry;
!            # spc with prod or loss terms in Smvgear (?)
!   ktloop : # of grid-cells in a grid-block
!   ncsp   : ncs       => for daytime   gas chemistry
!            ncs + ICS => for nighttime gas chemistry
!   cc2    : array holding values of decomposed matrix
!   vdiag  : 1 / current diagonal term of the decomposed matrix
!   gloss  : first derivative = sum of prod minus loss rates for a spc
!
!-----------------------------------------------------------------------------

      subroutine BacksubGMI  &
     &  (savedVars, ischan, ktloop, ncsp, cc2, vdiag, gloss)

      use GmiSavedVariables_mod, only : t_ChemistrySaved
      use GmiPrintError_mod, only : GmiPrintError

      implicit none

#     include "smv2chem_par.h"
!#     include "smv2chem2.h"


!     ----------------------
!     Argument declarations.
!     ----------------------

      integer, intent(in)  :: ischan
      integer, intent(in)  :: ktloop
      integer, intent(in)  :: ncsp
      real*8,  intent(in)  :: cc2  (KBLOOP, 0:MXARRAY)
      real*8,  intent(in)  :: vdiag(KBLOOP, MXGSAER)

      real*8,  intent(inout) :: gloss(KBLOOP, MXGSAER)
      type(t_ChemistrySaved), intent(inOut) :: savedVars


!     ----------------------
!     Variable declarations.
!     ----------------------

      integer :: i, ij
      integer :: ij0, ij1, ij2, ij3, ij4

      integer :: j0, j1, j2, j3, j4

      integer :: k, kc, kzt
      integer :: kh1, kh2, kh3, kh4, kh5
      integer :: kl1, kl2, kl3, kl4, kl5

      integer :: mc, mzt
      integer :: mh1, mh2, mh3, mh4, mh5
      integer :: ml1, ml2, ml3, ml4, ml5


!     ----------------
!     Begin execution.
!     ----------------

!c    Write (6,*) 'Backsub called.'


      ij = 1


!     ==========================================
      KZTLOOP: do kzt = savedVars%kztlo(ncsp), savedVars%kzthi(ncsp)
!     ==========================================

        i = savedVars%ikztot(kzt)

        kl5 = savedVars%kbl5(kzt)
        kh5 = savedVars%kbh5(kzt)
        kl4 = savedVars%kbl4(kzt)
        kh4 = savedVars%kbh4(kzt)
        kl3 = savedVars%kbl3(kzt)
        kh3 = savedVars%kbh3(kzt)
        kl2 = savedVars%kbl2(kzt)
        kh2 = savedVars%kbh2(kzt)
        kl1 = savedVars%kbl1(kzt)
        kh1 = savedVars%kbh1(kzt)

!       -- Sum 5 terms at a time. --

        do kc = kl5, kh5

          ij0 = ij
          ij1 = ij + 1
          ij2 = ij + 2
          ij3 = ij + 3
          ij4 = ij + 4
          ij  = ij + 5

          j0  = savedVars%kzeroa(kc)
          j1  = savedVars%kzerob(kc)
          j2  = savedVars%kzeroc(kc)
          j3  = savedVars%kzerod(kc)
          j4  = savedVars%kzeroe(kc)

          do k = 1, ktloop
            gloss(k,i) =  &
     &        gloss(k,i) -  &
     &        (cc2(k,ij0) * gloss(k,j0)) -  &
     &        (cc2(k,ij1) * gloss(k,j1)) -  &
     &        (cc2(k,ij2) * gloss(k,j2)) -  &
     &        (cc2(k,ij3) * gloss(k,j3)) -  &
     &        (cc2(k,ij4) * gloss(k,j4))
          end do

        end do

!       -- Sum 4 terms at a time. --

        do kc = kl4, kh4

          ij0 = ij
          ij1 = ij + 1
          ij2 = ij + 2
          ij3 = ij + 3
          ij  = ij + 4

          j0  = savedVars%kzeroa(kc)
          j1  = savedVars%kzerob(kc)
          j2  = savedVars%kzeroc(kc)
          j3  = savedVars%kzerod(kc)

          do k = 1, ktloop
            gloss(k,i) =  &
     &        gloss(k,i) -  &
     &        (cc2(k,ij0) * gloss(k,j0)) -  &
     &        (cc2(k,ij1) * gloss(k,j1)) -  &
     &        (cc2(k,ij2) * gloss(k,j2)) -  &
     &        (cc2(k,ij3) * gloss(k,j3))
          end do

        end do

!       -- Sum 3 terms at a time. --

        do kc = kl3, kh3

          ij0 = ij
          ij1 = ij + 1
          ij2 = ij + 2
          ij  = ij + 3

          j0  = savedVars%kzeroa(kc)
          j1  = savedVars%kzerob(kc)
          j2  = savedVars%kzeroc(kc)

          do k = 1, ktloop
            gloss(k,i) =  &
     &        gloss(k,i) -  &
     &        (cc2(k,ij0) * gloss(k,j0)) -  &
     &        (cc2(k,ij1) * gloss(k,j1)) -  &
     &        (cc2(k,ij2) * gloss(k,j2))
          end do

        end do

!       -- Sum 2 terms at a time. --

        do kc = kl2, kh2

          ij0 = ij
          ij1 = ij + 1
          ij  = ij + 2

          j0  = savedVars%kzeroa(kc)
          j1  = savedVars%kzerob(kc)

          do k = 1, ktloop
            gloss(k,i) =  &
     &        gloss(k,i) -  &
     &        (cc2(k,ij0) * gloss(k,j0)) -  &
     &        (cc2(k,ij1) * gloss(k,j1))
          end do

        end do

!       -- Sum 1 term at a time. --

        do kc = kl1, kh1

          ij0 = ij
          ij  = ij + 1

          j0  = savedVars%kzeroa(kc)

          do k = 1, ktloop
            gloss(k,i) =  &
     &        gloss(k,i) -  &
     &        (cc2(k,ij0) * gloss(k,j0))
          end do

        end do

!     ==============
      end do KZTLOOP
!     ==============


!     ---------------------------------------------------------------
!     Backsub loop # 2.
!
!     Backsubstite with upper triangular matrix to find solution.
!     Again, sum up several terms at a time to improve vectorization.
!     ---------------------------------------------------------------

!     ===========================
      ILOOP: do i = ischan, 1, -1
!     ===========================

        mzt = savedVars%imztot(i,ncsp)

!       ===================
        MZTIF: if (mzt > 0) then
!       ===================

          ml5 = savedVars%mbl5(mzt)
          mh5 = savedVars%mbh5(mzt)
          ml4 = savedVars%mbl4(mzt)
          mh4 = savedVars%mbh4(mzt)
          ml3 = savedVars%mbl3(mzt)
          mh3 = savedVars%mbh3(mzt)
          ml2 = savedVars%mbl2(mzt)
          mh2 = savedVars%mbh2(mzt)
          ml1 = savedVars%mbl1(mzt)
          mh1 = savedVars%mbh1(mzt)

!         -- Sum 5 terms at a time. --

          do mc = ml5, mh5

            ij0 = ij
            ij1 = ij + 1
            ij2 = ij + 2
            ij3 = ij + 3
            ij4 = ij + 4
            ij  = ij + 5

            j0  = savedVars%mzeroa(mc)
            j1  = savedVars%mzerob(mc)
            j2  = savedVars%mzeroc(mc)
            j3  = savedVars%mzerod(mc)
            j4  = savedVars%mzeroe(mc)

            do k = 1, ktloop
              gloss(k,i) =  &
     &          gloss(k,i) -  &
     &          (cc2(k,ij0) * gloss(k,j0)) -  &
     &          (cc2(k,ij1) * gloss(k,j1)) -  &
     &          (cc2(k,ij2) * gloss(k,j2)) -  &
     &          (cc2(k,ij3) * gloss(k,j3)) -  &
     &          (cc2(k,ij4) * gloss(k,j4))
            end do

          end do

!         -- Sum 4 terms at a time. --

          do mc = ml4, mh4

            ij0 = ij
            ij1 = ij + 1
            ij2 = ij + 2
            ij3 = ij + 3
            ij  = ij + 4

            j0  = savedVars%mzeroa(mc)
            j1  = savedVars%mzerob(mc)
            j2  = savedVars%mzeroc(mc)
            j3  = savedVars%mzerod(mc)

            do k = 1, ktloop
              gloss(k,i) =  &
     &          gloss(k,i) -  &
     &          (cc2(k,ij0) * gloss(k,j0)) -  &
     &          (cc2(k,ij1) * gloss(k,j1)) -  &
     &          (cc2(k,ij2) * gloss(k,j2)) -  &
     &          (cc2(k,ij3) * gloss(k,j3))
            end do

          end do

!         -- Sum 3 terms at a time. --

          do mc = ml3, mh3

            ij0 = ij
            ij1 = ij + 1
            ij2 = ij + 2
            ij  = ij + 3

            j0  = savedVars%mzeroa(mc)
            j1  = savedVars%mzerob(mc)
            j2  = savedVars%mzeroc(mc)

            do k = 1, ktloop
              gloss(k,i) =  &
     &          gloss(k,i) -  &
     &          (cc2(k,ij0) * gloss(k,j0)) -  &
     &          (cc2(k,ij1) * gloss(k,j1)) -  &
     &          (cc2(k,ij2) * gloss(k,j2))
            end do

          end do

!         -- Sum 2 terms at a time. --

          do mc = ml2, mh2

            ij0 = ij
            ij1 = ij + 1
            ij  = ij + 2

            j0  = savedVars%mzeroa(mc)
            j1  = savedVars%mzerob(mc)

            do k = 1, ktloop
              gloss(k,i) =  &
     &          gloss(k,i) -  &
     &          (cc2(k,ij0) * gloss(k,j0)) -  &
     &          (cc2(k,ij1) * gloss(k,j1))
            end do

          end do

!         -- Sum 1 term at a time. --

          do mc = ml1, mh1

            ij0 = ij
            ij  = ij + 1

            j0  = savedVars%mzeroa(mc)

            do k = 1, ktloop
              gloss(k,i) =  &
     &          gloss(k,i) -  &
     &          (cc2(k,ij0) * gloss(k,j0))
            end do

          end do

!       ============
        end if MZTIF
!       ============

!       -- Adjust gloss with diagonal element. --

        do k = 1, ktloop
          gloss(k,i) = gloss(k,i) * vdiag(k,i)
        end do

!     ============
      end do ILOOP
!     ============


      return

      end subroutine BacksubGMI


!-----------------------------------------------------------------------------
!
! ROUTINE
!   Decomp
!
! DESCRIPTION
!   This routine decomposes the sparse matrix "P" into the matrix "A" in
!   order to solve the linear set of equations Ax = B for x, which is a
!   correction vector.  Ax = B is solved in Backsub, the original matrix
!   "P" is =>
!
!     P = I - H x Bo x J
!
!   where I = identity matrix, H = time step, Bo = a coefficient that
!   depends on the order of the integration method, and J is the matrix of
!   partial derivatives.  See Press, et. al. (1992), Numerical Recipes,
!   Cambridge University Press, for a better description of the L-U
!   decompostion process.
!
!   This L-U decompostion process uses sparse matrix techniques, vectorizes
!   around the grid-cell dimension, and uses no partial pivoting.  Tests by
!   Sherman & Hindmarsh (1980), Lawrence Livermore National Laboratory,
!   Rep. UCRL-84102, and by us have confirmed that the removal of partial
!   pivoting has little effect on results.
!
! ARGUMENTS
!   ischan : # of first-order eqns to solve, = # of spc = order of original
!            matrix; ischan has a different value for day and night, and for
!            gas- and aqueous-phase chemistry;
!            # spc with prod or loss terms in Smvgear (?)
!   ktloop : # of grid-cells in a grid-block
!   ncsp   : ncs       => for daytime   gas chemistry
!            ncs + ICS => for nighttime gas chemistry
!   cc2    : array of iarray units holding values of each matrix
!            position actually used; originally,
!            cc2 = P = I - delt * aset(nqq,1) * partial_derivatives;
!            however, cc2 is decomposed here
!   vdiag  : 1 / current diagonal term of the decomposed matrix
!
!-----------------------------------------------------------------------------

      subroutine DecompGMI &
     &  (savedVars, ischan, ktloop, ncsp, cc2, vdiag)

      use GmiSavedVariables_mod, only : t_ChemistrySaved
      use GmiPrintError_mod, only : GmiPrintError

      implicit none
#     include "smv2chem_par.h"
!#     include "smv2chem2.h"


!     ----------------------
!     Argument declarations.
!     ----------------------

      integer, intent(in)  :: ischan
      integer, intent(in)  :: ktloop
      integer, intent(in)  :: ncsp

      real*8,  intent(inout) :: cc2  (KBLOOP, 0:MXARRAY)
      real*8,  intent(inout) :: vdiag(KBLOOP, MXGSAER)
      type(t_ChemistrySaved), intent(inOut) :: savedVars


!     ----------------------
!     Variable declarations.
!     ----------------------

      integer :: iar, ic
      integer :: ih1, ih2, ih3, ih4, ih5
      integer :: ij, ija, ijt
      integer :: ik0, ik1, ik2, ik3, ik4
      integer :: il1, il2, il3, il4, il5
      integer :: j, jc, jh, jl
      integer :: k
      integer :: kj0, kj1, kj2, kj3, kj4


!     ----------------
!     Begin execution.
!     ----------------

!c    Write (6,*) 'Decomp called.'


!     -----------------------------------------------------------
!     First loop of L-U decompostion.
!
!     Sum 1,2,3,4, OR 5 terms at a time to improve vectorization.
!     -----------------------------------------------------------

!     =======================
      JLOOP: do j = 1, ischan
!     =======================

!       ==============================================
        IJTLOOP: do ijt = savedVars%ijtlo(j,ncsp), savedVars%ijthi(j,ncsp)
!       ==============================================

          ij  = savedVars%ijval(ijt)
          il5 = savedVars%idl5 (ijt)
          ih5 = savedVars%idh5 (ijt)
          il4 = savedVars%idl4 (ijt)
          ih4 = savedVars%idh4 (ijt)
          il3 = savedVars%idl3 (ijt)
          ih3 = savedVars%idh3 (ijt)
          il2 = savedVars%idl2 (ijt)
          ih2 = savedVars%idh2 (ijt)
          il1 = savedVars%idl1 (ijt)
          ih1 = savedVars%idh1 (ijt)

!         -- Sum 5 terms at a time. --

          do ic = il5, ih5

            ik0 = savedVars%ikdeca(ic)
            ik1 = savedVars%ikdecb(ic)
            ik2 = savedVars%ikdecc(ic)
            ik3 = savedVars%ikdecd(ic)
            ik4 = savedVars%ikdece(ic)

            kj0 = savedVars%kjdeca(ic)
            kj1 = savedVars%kjdecb(ic)
            kj2 = savedVars%kjdecc(ic)
            kj3 = savedVars%kjdecd(ic)
            kj4 = savedVars%kjdece(ic)

            do k = 1, ktloop
              cc2(k,ij) =  &
     &          cc2(k,ij) -  &
     &          (cc2(k,ik0) * cc2(k,kj0)) -  &
     &          (cc2(k,ik1) * cc2(k,kj1)) -  &
     &          (cc2(k,ik2) * cc2(k,kj2)) -  &
     &          (cc2(k,ik3) * cc2(k,kj3)) -  &
     &          (cc2(k,ik4) * cc2(k,kj4))
            end do

          end do

!         -- Sum 4 terms at a time. --

          do ic = il4, ih4

            ik0 = savedVars%ikdeca(ic)
            ik1 = savedVars%ikdecb(ic)
            ik2 = savedVars%ikdecc(ic)
            ik3 = savedVars%ikdecd(ic)

            kj0 = savedVars%kjdeca(ic)
            kj1 = savedVars%kjdecb(ic)
            kj2 = savedVars%kjdecc(ic)
            kj3 = savedVars%kjdecd(ic)

            do k = 1, ktloop
              cc2(k,ij) =  &
     &          cc2(k,ij) -  &
     &          (cc2(k,ik0) * cc2(k,kj0)) -  &
     &          (cc2(k,ik1) * cc2(k,kj1)) -  &
     &          (cc2(k,ik2) * cc2(k,kj2)) -  &
     &          (cc2(k,ik3) * cc2(k,kj3))
            end do

          end do

!         -- Sum 3 terms at a time. --

          do ic = il3, ih3

            ik0 = savedVars%ikdeca(ic)
            ik1 = savedVars%ikdecb(ic)
            ik2 = savedVars%ikdecc(ic)

            kj0 = savedVars%kjdeca(ic)
            kj1 = savedVars%kjdecb(ic)
            kj2 = savedVars%kjdecc(ic)

            do k = 1, ktloop
              cc2(k,ij) =  &
     &          cc2(k,ij) -  &
     &          (cc2(k,ik0) * cc2(k,kj0)) -  &
     &          (cc2(k,ik1) * cc2(k,kj1)) -  &
     &          (cc2(k,ik2) * cc2(k,kj2))
            end do

          end do

!         -- Sum 2 terms at a time. --

          do ic = il2, ih2

            ik0 = savedVars%ikdeca(ic)
            ik1 = savedVars%ikdecb(ic)

            kj0 = savedVars%kjdeca(ic)
            kj1 = savedVars%kjdecb(ic)

            do k = 1, ktloop
              cc2(k,ij) =  &
     &          cc2(k,ij) -  &
     &          (cc2(k,ik0) * cc2(k,kj0)) -  &
     &          (cc2(k,ik1) * cc2(k,kj1))
            end do

          end do

!         -- Sum 1 term  at a time. --

          do ic = il1, ih1

            ik0 = savedVars%ikdeca(ic)

            kj0 = savedVars%kjdeca(ic)

            do k = 1, ktloop
              cc2(k,ij) =  &
     &          cc2(k,ij) -  &
     &          (cc2(k,ik0) * cc2(k,kj0))
            end do

          end do

!       ==============
        end do IJTLOOP
!       ==============

        iar = savedVars%jarrdiag(j,ncsp)

        do k = 1, ktloop
          vdiag(k,j) = 1.0d0 / cc2(k,iar)
        end do

!       ----------------------------
!       Second loop of decompostion.
!       ----------------------------

        jl = savedVars%jloz1(j,ncsp)
        jh = savedVars%jhiz1(j,ncsp)

        do jc = jl, jh

          ija = savedVars%jzeroa(jc)

          do k = 1, ktloop
            cc2(k,ija) = cc2(k,ija) * vdiag(k,j)
          end do

        end do

!     ============
      end do JLOOP
!     ============


      return

      end subroutine DecompGMI


!-----------------------------------------------------------------------------
!
! ROUTINE
!   Pderiv
!
! DESCRIPTION
!   This routine puts the partial derivatives of each ordinary differential
!   equation into a matrix.  The form of the matrix equation is =>
!
!     P = I - H x Bo x J
!
!   where I = identity matrix, H = time step, Bo = coefficient corresponding
!   to the order of the method, and J is the Jacobian matrix of partial
!   derivatives.
!
!   Example of how partial derivatives are placed in an array =>
!
!     Species:          A,   B,   C
!     Concentrations:  (A), (B), (C)
!
!     Reactions:    1) A          --> B      J
!                   2) A  + B     --> C      K1
!                   3) A  + B + C --> D      K2
!
!     First         d(A) / dt = -J(A) - K1(A)(B) - K2(A)(B)(C)
!     Derivatives:  d(B) / dt = +J(A) - K1(A)(B) - K2(A)(B)(C)
!                   d(C) / dt =       + K1(A)(B) - K2(A)(B)(C)
!                   d(D) / dt =                  + K2(A)(B)(C)
!
!     Predictor matrix (P) = I - h * b * J:
!       J = Jacobian matrix of partial derivates
!       I = identity matrix
!       h = time step
!       b = coefficient of method
!       R = h * b = -r1delt
!
!         A                      B                     C                  D
!      _____________________________________________________________________
!     |
!   A | 1-R(-J-K1(B)-K2(B)(C))  -R(-K1(A)-K2(A)(C))   -R(-K2(A)(B))       0
!     |
!   B |  -R(+J-K1(B)-K2(B)(C)) 1-R(-K1(A)-K2(A)(C))   -R(-K2(A)(B))       0
!     |
!   C |  -R(  +K1(B)-K2(B)(C))  -R(+K1(A)-K2(A)(C))  1-R(-K2(A)(B))       0
!     |
!   D |  -R(        +K2(B)(C))  -R(      +K2(A)(C))   -R(+K2(A)(B))       1
!
!
! ARGUMENTS
!   ischan   : # of first-order eqns to solve, = # of spc = order of original
!              matrix; ischan has a different value for day and night, and for
!              gas- and aqueous-phase chemistry;
!              # spc with prod or loss terms in Smvgear (?)
!   ktloop   : # of grid-cells in a grid-block
!   ncsp     : ncs       => for daytime   gas chemistry
!              ncs + ICS => for nighttime gas chemistry
!   nfdh2    : nfdh3 + # of rxns with two   active reactants
!   nfdh3    :         # of rxns with three active reactants
!   nfdl1    : nfdh2 + 1
!   nfdl2    : nfdh3 + 1
!   irma,b,c : spc # of each reactant; locates reordered active spc #s
!   r1delt   : = -aset(nqq,1) * time_step = -coefficient_of_method * dt
!   cnew     : init (and final) spc conc (# cm-3-air or moles l-1-h2o (?))
!   rrate    : rxn rate coefficient =>
!                rates w. 1 reactant:   s^-1
!                rates w. 2 reactants:  l-h2o mole-1 s^-1 or cm^3 #-1 s^-1 (?)
!                rates w. 3 reactants:  l^2-h2o m-2 s^-1  or cm^6 #-2 s^-1 (?)
!   npderiv  : total # of times Pderiv is called
!   cc2      : array of iarray units holding values of each matrix
!              position actually used;
!              cc2 = P = I - delt * aset(nqq,1) * partial_derivatives
!   urate    : term of Jacobian (J) = partial derivative
!   nfdh1    : nfdh2 + # of rxns with one   active reactant
!
!-----------------------------------------------------------------------------

      subroutine PderivGMI  &
     &  (savedVars, ischan, ktloop, ncsp, nfdh2, nfdh3, nfdl1, nfdl2, irma, irmb,  &
     &   irmc, r1delt, cnew, rrate, npderiv, cc2, urate, nfdh1)
   
      use GmiSavedVariables_mod, only : t_ChemistrySaved
      use GmiPrintError_mod, only : GmiPrintError

      implicit none
#     include "smv2chem_par.h"
!#     include "smv2chem2.h"


!     ----------------------
!     Argument declarations.
!     ----------------------

      integer, intent(in)  :: ischan
      integer, intent(in)  :: ktloop
      integer, intent(in)  :: ncsp
      integer, intent(in)  :: nfdh2, nfdh3
      integer, intent(in)  :: nfdl1, nfdl2
      integer, intent(in)  :: irma (NMTRATE)
      integer, intent(in)  :: irmb (NMTRATE)
      integer, intent(in)  :: irmc (NMTRATE)
      real*8,  intent(in)  :: r1delt
      real*8,  intent(in)  :: cnew (KBLOOP, MXGSAER)
      real*8,  intent(in)  :: rrate(KBLOOP, NMTRATE)

      integer, intent(inout) :: npderiv
      real*8,  intent(inout) :: cc2  (KBLOOP, 0:MXARRAY)
      real*8,  intent(inout) :: urate(KBLOOP, NMTRATE, 3)

      integer, intent(out) :: nfdh1
      type(t_ChemistrySaved), intent(inOut) :: savedVars


!     ----------------------
!     Variable declarations.
!     ----------------------

      integer :: ial, iar
      integer :: iarry
      integer :: ja, jb, jc
      integer :: k
      integer :: n, nkn
      integer :: nondiag     ! # of final matrix positions, excluding diagonal
                             ! terms, filled after all matrix processes
      integer :: nondiag1    ! nondiag + 1
      integer :: npdh, npdl

      real*8  :: fracr1


!     ----------------
!     Begin execution.
!     ----------------

!c    Write (6,*) 'Pderiv called.'


!     -----------------------------------------------------------
!     Calculate partial derivatives and sum up partial derivative
!     loss terms.
!     -----------------------------------------------------------

      npderiv  = npderiv + 1
      iarry    = savedVars%iarray(ncsp)
      nondiag  = iarry - ischan
      nondiag1 = nondiag + 1
      nfdh1    = nfdh2 + savedVars%ioner(ncsp)
      npdl     = savedVars%npdlo(ncsp)
      npdh     = savedVars%npdhi(ncsp)


!     -----------------------------------------------------------
!     Partial derivatives for rates with three active loss terms.
!     -----------------------------------------------------------

      do nkn = 1, nfdh3

        ja = irma(nkn)
        jb = irmb(nkn)
        jc = irmc(nkn)

        do k = 1, ktloop
          urate(k,nkn,1) = rrate(k,nkn) * cnew(k,jb) * cnew(k,jc)
          urate(k,nkn,2) = rrate(k,nkn) * cnew(k,ja) * cnew(k,jc)
          urate(k,nkn,3) = rrate(k,nkn) * cnew(k,ja) * cnew(k,jb)
        end do

      end do


!     ---------------------------------------------------------
!     Partial derivatives for rates with two active loss terms.
!     ---------------------------------------------------------

      do nkn = nfdl2, nfdh2

        ja = irma(nkn)
        jb = irmb(nkn)

        do k = 1, ktloop
          urate(k,nkn,1) = rrate(k,nkn) * cnew(k,jb)
          urate(k,nkn,2) = rrate(k,nkn) * cnew(k,ja)
        end do

      end do


!     --------------------------------------------------------
!     Partial derivatives for rates with one active loss term.
!     --------------------------------------------------------

      do nkn = nfdl1, nfdh1
        do k = 1, ktloop
          urate(k,nkn,1) = rrate(k,nkn)
        end do
      end do


!     ------------------------------------------------------------------
!     Put partial derivatives production and loss terms in matrix array.
!     ------------------------------------------------------------------

      do iar = 1, nondiag
        do k = 1, ktloop
          cc2(k,iar) = 0.0d0
        end do
      end do


      do iar = nondiag1, iarry
        do k = 1, ktloop
          cc2(k,iar) = 1.0d0
        end do
      end do


      do n = npdl, npdh

        nkn    = savedVars%nkpdterm(n)
        iar    = savedVars%ipospd  (n)
        ial    = savedVars%iialpd  (n)
        fracr1 = savedVars%fracpl  (n) * r1delt

        do k = 1, ktloop
          cc2(k,iar) = cc2(k,iar) + (fracr1 * urate(k,nkn,ial))
        end do

      end do


      return

      end subroutine PderivGMI


!-----------------------------------------------------------------------------
!
! ROUTINE
!   Subfun
!
! DESCRIPTION
!   This routine evaluates the first derivative of each ordinary
!   differential equation (ODE).  It evaluates derivatives in the special
!   form f = y'(est) = f(x,y,estimated), where f is the right hand side of
!   the differential equation.
!
!   Example =>
!
!     Species:         A,   B,   C
!     Concentrations: (A), (B), (C)
!
!     Reactions:    1) A          --> B      J
!                   2) A  + B     --> C      K1
!                   3) A  + B + C --> D      K2
!
!     First         d(A) / dt = -J(A) - K1(A)(B) - K2(A)(B)(C)
!     Derivatives:  d(B) / dt = +J(A) - K1(A)(B) - K2(A)(B)(C)
!                   d(C) / dt =       + K1(A)(B) - K2(A)(B)(C)
!                   d(D) / dt =                  + K2(A)(B)(C)
!
! ARGUMENTS
!   ischan   : # of first-order eqns to solve, = # of spc = order of original
!              matrix; ischan has a different value for day and night, and for
!              gas- and aqueous-phase chemistry;
!              # spc with prod or loss terms in Smvgear (?)
!   ktloop   : # of grid-cells in a grid-block
!   nallr    : # of active rxns
!   ncsp     : ncs       => for daytime   gas chemistry
!              ncs + ICS => for nighttime gas chemistry
!   nfdh2    : nfdh3 + # of rxns with two   active reactants
!   nfdh3    :         # of rxns with three active reactants
!   nfdl1    : nfdh2 + 1
!   nfdl2    : nfdh3 + 1
!   nfdrep   : nfdh3 + # of rxns with two active reactants that are not
!              followed by a rxn with the same reactants
!   nfdrep1  : nfdrep + 1
!   irma,b,c : spc # of each reactant; locates reordered active spc #s
!   cnew     : init (and final) spc conc (# cm^-3-air or moles l^-1-h2o (?))
!   rrate    : rxn rate coefficient =>
!                rates with 1 reactant:   s^-1
!                rates with 2 reactants:  l-h2o mole^-1 s^-1 or
!                                         cm^3 #-1 s^-1 (?)
!                rates with 3 reactants:  l^2-h2o m-2 s-1  or cm^6 #-2 s-1 (?)
!   nsubfun  : total # of times Subfun is called
!   gloss    : first derivative = sum of prod. minus loss rates for a spc
!   trate    : rxn rate (moles l^-1-h2o s^-1 or # cm^-3 s^-1 (?))
!   nfdh1    : nfdh2 + # of rxns with one   active reactant
!
!-----------------------------------------------------------------------------

      subroutine SubfunGMI  &
     &  (savedVars, ischan, ktloop, nallr, ncsp, nfdh2, nfdh3, nfdl1, nfdl2,  &
     &   nfdrep, nfdrep1, irma, irmb, irmc, cnew, rrate, nsubfun,  &
     &   gloss, trate, nfdh1)

      use GmiSavedVariables_mod, only : t_ChemistrySaved
      use GmiPrintError_mod, only : GmiPrintError

      implicit none
#     include "smv2chem_par.h"
!#     include "smv2chem2.h"


!     ----------------------
!     Argument declarations.
!     ----------------------

      integer, intent(in)  :: ischan
      integer, intent(in)  :: ktloop
      integer, intent(in)  :: nallr
      integer, intent(in)  :: ncsp
      integer, intent(in)  :: nfdh2, nfdh3
      integer, intent(in)  :: nfdl1, nfdl2
      integer, intent(in)  :: nfdrep, nfdrep1
      integer, intent(in)  :: irma (NMTRATE)
      integer, intent(in)  :: irmb (NMTRATE)
      integer, intent(in)  :: irmc (NMTRATE)
      real*8,  intent(in)  :: cnew (KBLOOP, MXGSAER)
      real*8,  intent(in)  :: rrate(KBLOOP, NMTRATE)

      integer, intent(inout) :: nsubfun
      real*8,  intent(inout) :: gloss(KBLOOP, MXGSAER)
      real*8,  intent(inout) :: trate(KBLOOP, NMTRATE*2)

      integer, intent(out) :: nfdh1
      type(t_ChemistrySaved), intent(inOut) :: savedVars


!     ----------------------
!     Variable declarations.
!     ----------------------

      integer ::  ja, jb, jc
      integer ::  jspc
      integer ::  k
      integer ::  n, nc, nh
      integer ::  nh1, nh2, nh3, nh4, nh5
      integer ::  nk0, nk1, nk2
      integer ::  nk3, nk4, nkn
      integer ::  nl1, nl2, nl3, nl4, nl5
      integer ::  npl

!     -----------------------------------------------------------------------
!     concmult : product of concs in a rate; if two consecutive rxns have the
!                same spc reacting (e.g., A + B --> C and A + B --> D + E),
!                then use the same value for both rxns
!     -----------------------------------------------------------------------

      real*8  :: concmult
      real*8  :: fracn


!     ----------------
!     Begin execution.
!     ----------------

!c    Write (6,*) 'Subfun called.'


!     ----------------------
!     Set rates of reaction.
!     ----------------------

      nsubfun = nsubfun + 1
      nfdh1   = nfdh2 + savedVars%ioner(ncsp)


!     ---------------------------------------------------------
!     First derivatives for rates with three active loss terms.
!     ---------------------------------------------------------

      do nkn = 1, nfdh3

        ja = irma(nkn)
        jb = irmb(nkn)
        jc = irmc(nkn)

        nh = nkn + nallr

        do k = 1, ktloop
          trate(k,nkn) =  &
     &      rrate(k,nkn) * cnew(k,ja) * cnew(k,jb) * cnew(k,jc)
          trate(k,nh)  = -trate(k,nkn)
        end do

      end do


!     -------------------------------------------------------
!     First derivatives for rates with two active loss terms.
!     -------------------------------------------------------

      do nkn = nfdl2, nfdrep

        ja = irma(nkn)
        jb = irmb(nkn)

        nh = nkn + nallr

        do k = 1, ktloop
          trate(k,nkn) = rrate(k,nkn) * cnew(k,ja) * cnew(k,jb)
          trate(k,nh)  = -trate(k,nkn)
        end do

      end do


!     -----------------------------------------------------------
!     First derivatives for rates with two active loss terms and
!     where the subsequent reaction has the same reactants, but a
!     different rate.
!     -----------------------------------------------------------

      do nkn = nfdrep1, nfdh2, 2

        ja  = irma(nkn)
        jb  = irmb(nkn)

        nk2 = nkn + 1
        nh  = nkn + nallr
        nh2 = nk2 + nallr

        do k = 1, ktloop
          concmult     =  cnew(k,ja)   * cnew(k,jb)
          trate(k,nkn) =  rrate(k,nkn) * concmult
          trate(k,nk2) =  rrate(k,nk2) * concmult
          trate(k,nh)  = -trate(k,nkn)
          trate(k,nh2) = -trate(k,nk2)
        end do

      end do


!     ------------------------------------------------------
!     First derivatives for rates with one active loss term.
!     ------------------------------------------------------

      do nkn = nfdl1, nfdh1

        ja = irma(nkn)

        nh = nkn + nallr

        do k = 1, ktloop
          trate(k,nkn) =  rrate(k,nkn) * cnew(k,ja)
          trate(k,nh)  = -trate(k,nkn)
        end do

      end do


!     --------------------------------
!     Initialize first derivative = 0.
!     --------------------------------

      do jspc = 1, ischan
        do k = 1, ktloop
          gloss(k,jspc) = 0.0d0
        end do
      end do


!     ---------------------------------------------------------------
!     Sum net (not reproduced) kinetic and photo gains and losses for
!     each species.
!
!     Sum 1,2,3,4, or 5 terms at a time to improve vectorization.
!     ---------------------------------------------------------------

!     ==========================================
      NPLLOOP: do npl = savedVars%npllo(ncsp), savedVars%nplhi(ncsp)
!     ==========================================

        jspc = savedVars%jspnpl(npl)

        nl5 = savedVars%npl5(npl)
        nh5 = savedVars%nph5(npl)
        nl4 = savedVars%npl4(npl)
        nh4 = savedVars%nph4(npl)
        nl3 = savedVars%npl3(npl)
        nh3 = savedVars%nph3(npl)
        nl2 = savedVars%npl2(npl)
        nh2 = savedVars%nph2(npl)
        nl1 = savedVars%npl1(npl)
        nh1 = savedVars%nph1(npl)

!       -- Sum 5 terms at a time. --

        do nc = nl5, nh5

          nk0 = savedVars%lossra(nc)
          nk1 = savedVars%lossrb(nc)
          nk2 = savedVars%lossrc(nc)
          nk3 = savedVars%lossrd(nc)
          nk4 = savedVars%lossre(nc)

          do k = 1, ktloop
            gloss(k,jspc) =  &
     &        gloss(k,jspc) -  &
     &        trate(k,nk0)  - trate(k,nk1) - trate(k,nk2) -  &
     &        trate(k,nk3)  - trate(k,nk4)
          end do

        end do

!       -- Sum 4 terms at a time. --

        do nc = nl4, nh4

          nk0 = savedVars%lossra(nc)
          nk1 = savedVars%lossrb(nc)
          nk2 = savedVars%lossrc(nc)
          nk3 = savedVars%lossrd(nc)

          do k = 1, ktloop
            gloss(k,jspc) =  &
     &        gloss(k,jspc) -  &
     &        trate(k,nk0)  - trate(k,nk1) - trate(k,nk2) -  &
     &        trate(k,nk3)
          end do

        end do

!       -- Sum 3 terms at a time. --

        do nc = nl3, nh3

          nk0 = savedVars%lossra(nc)
          nk1 = savedVars%lossrb(nc)
          nk2 = savedVars%lossrc(nc)

          do k = 1, ktloop
            gloss(k,jspc) =  &
     &        gloss(k,jspc) -  &
     &        trate(k,nk0)  - trate(k,nk1)  - trate(k,nk2)
          end do

        end do

!       -- Sum 2 terms at a time. --

        do nc = nl2, nh2

          nk0 = savedVars%lossra(nc)
          nk1 = savedVars%lossrb(nc)

          do k = 1, ktloop
            gloss(k,jspc) =  &
     &        gloss(k,jspc) -  &
     &        trate(k,nk0)  - trate(k,nk1)
          end do

        end do

!       -- Sum 1 term at a time. --

        do nc = nl1, nh1

          nk0 = savedVars%lossra(nc)

          do k = 1, ktloop
            gloss(k,jspc) =  &
     &        gloss(k,jspc) -  &
     &        trate(k,nk0)
          end do

        end do

!     ==============
      end do NPLLOOP
!     ==============


!     --------------------------------------------------------------
!     Sum production term for reactions where products fractionated.
!     --------------------------------------------------------------

      do n = savedVars%nfrlo(ncsp), savedVars%nfrhi(ncsp)

        jspc  = savedVars%jspcnfr(n)
        nkn   = savedVars%nknfr  (n)
        fracn = savedVars%fracnfr(n)

        do k = 1, ktloop
          gloss(k,jspc) = gloss(k,jspc) + (fracn * trate(k,nkn))
        end do

      end do


      return

      end subroutine SubfunGMI


!-----------------------------------------------------------------------------
!
! ROUTINE
!   Update
!
! DESCRIPTION
!   This routine updates photodissociation rates.
!
!   Photorates are included in first and partial derivative equations.
!
! ARGUMENTS
!   ktloop   : # of grid-cells in a grid-block
!   nallr    : # of active rxns
!   ncs      : identifies gas chemistry type (1..NCSGAS)
!   ncsp     : ncs       => for daytime   gas chemistry
!              ncs + ICS => for nighttime gas chemistry
!   jphotrat : tbd
!   ptratk1  : tbd
!   rrate    : rate constants
!   trate    : rxn rate (moles l^-1-h2o s^-1 or # cm^-3 s^-1 (?))
!
!-----------------------------------------------------------------------------

      subroutine UpdateGMI  &
     &  (savedVars, ktloop, nallr, ncs, ncsp, jphotrat, pratk1, rrate, trate)

      use GmiSavedVariables_mod, only : t_ChemistrySaved
      use GmiPrintError_mod, only : GmiPrintError

      implicit none
#     include "smv2chem_par.h"
!#     include "smv2chem2.h"


!     ----------------------
!     Argument declarations.
!     ----------------------

      integer, intent(in)  :: ktloop
      integer, intent(in)  :: nallr
      integer, intent(in)  :: ncs
      integer, intent(in)  :: ncsp
      integer, intent(in)  :: jphotrat(ICS)
      real*8,  intent(in)  :: pratk1  (KBLOOP, IPHOT)

      real*8,  intent(inout) :: rrate(KBLOOP, NMTRATE)
      real*8,  intent(inout) :: trate(KBLOOP, NMTRATE*2)

      type(t_ChemistrySaved), intent(inOut) :: savedVars

!     ----------------------
!     Variable declarations.
!     ----------------------

      integer :: i, j
      integer :: kloop
      integer :: nh, nk, nkn


!     ----------------
!     Begin execution.
!     ----------------

!c    Write (6,*) 'Update called.'


!     -------------------------------
!     Load photolysis rate constants.
!     -------------------------------

      do j = 1, jphotrat(ncs)

        nkn = savedVars%nknphotrt(j,ncs)

        do kloop = 1, ktloop
          rrate(kloop,nkn) = pratk1(kloop,j)
        end do

      end do


!     ------------------------------------------------------
!     Set rates where photoreaction has no active loss term.
!     ------------------------------------------------------

      do i = 1, savedVars%nolosp(ncsp)

        nk  = savedVars%nknlosp(i,ncs)
        nkn = savedVars%newfold(nk,ncs)
        nh  = nkn + nallr

        do kloop = 1, ktloop
          trate(kloop,nkn) =  rrate(kloop,nkn)
          trate(kloop,nh)  = -trate(kloop,nkn)
        end do

      end do

      return

      end subroutine UpdateGMI

