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
!   jsparseGMI.F90
!
! ROUTINES
!   JsparseGMI
!   KsparseGMI
!
!=============================================================================


!-----------------------------------------------------------------------------
!
! ROUTINE
!   JsparseGMI
!
! DESCRIPTION
!   This routine sets up sparse matrix and other arrays for Smvgear
!   (Sparse Matrix Vectorized Gear code).  It sets arrays for gas-phase,
!   aqueous-phase, and any other type of chemistry.  It also sets arrays
!   for both day and night chemistry of each type.
!
!   Sets up arrays for gas- and aqueous-phase chemistry =>
!
!   Includes arrays for calculating first derivatives, partial derivatives,
!   matrix decompostion, and matrix back-substitution.  First, Jsparse
!   reorders the ordinary differential equations to maximize the number of
!   zeros in the matrix of partial derivatives.  It later sets arrays to
!   eliminate all calculations involving a zero.
!
!   GEOS5: JsparseGMI avoids name comflict with jsparse.F90 in GEOSCHEMchem
!
! ARGUMENTS
!   pr_smv2   : should the SmvgearII output file be written
!               (non-parallel mode only)?
!   lunsmv    : logical unit number to write to when pr_smv2 is true
!   ncs       : identifies gas chemistry type (1..NCSGAS)
!   ih2o      : identifies spc # of water
!   imgas     : array index for air density
!   initrogen : identifies spc # of nitrogen gas
!   ioxygen   : identifies spc # of oxygen   gas
!   jphotrat  : tbd
!   nrates    : # of kinetic rxns (non-photo)
!   ntspec    : # of active + inactive gases
!   inewold   : original spc # of each new jnew spc
!   npphotrat : tbd
!
!-----------------------------------------------------------------------------

      subroutine JsparseGMI(savedVars, pr_smv2, lunsmv, ncs, &
                          jphotrat, nrates, ntspec, inewold, npphotrat,  &
                          ncsp, lzero, jarraypt, &
                          iccount, jccount, kccount, mccount, icnt,     &
                          jcnt, kcnt, mcnt, idecomp, ijtot, kztot,      &
                          mztot, kbsub, mbsub, fkoef, irm, namesp2)

      use GmiPrintError_mod, only : GmiPrintError
      use GmiSavedVariables_mod, only : t_ChemistrySaved

      implicit none

#     include "smv2chem_par.h"
!#     include "smv2chem2.h"

!     ----------------------
!     Argument declarations.
!     ----------------------

      logical, intent(in)  :: pr_smv2
      integer, intent(in)  :: lunsmv
      integer, intent(in)  :: ncs

      integer, intent(out) :: jphotrat (ICS)
      integer, intent(out) :: nrates   (ICS)
      integer, intent(out) :: ntspec   (ICS)
      integer, intent(out) :: inewold  (MXGSAER, ICS)
      integer, intent(out) :: npphotrat(IPHOT,   ICS)

      integer :: ncsp
      integer :: jarraypt(MXGSAER, MXGSAER)
      integer :: lzero   (MXGSAER, MXGSAER)
      integer :: iccount, jccount, kccount, mccount
      integer :: icnt, jcnt, kcnt, mcnt
      integer :: idecomp
      integer :: ijtot, kztot, mztot
      integer :: kbsub, mbsub


      type(t_ChemistrySaved), intent(inOut) :: savedVars


!     ----------------------
!     Variable declarations.
!     ----------------------

      character (len=14) :: namesp2(0:MXGSAER, ICS)

      integer :: i, j, k, l, n
      integer :: i1
      integer :: ial, iap, iar
      integer :: ib, ic
      integer :: icb, icd
      integer :: ifnone
      integer :: ifsun      ! 1 => for daytime, 2 => for nighttime
      integer :: ihireac
      integer :: iminnew, iminold
      integer :: inew, iold
      integer :: iplsithrm
      integer :: iporr
      integer :: ipr, iprod
      integer :: ir, ire
      integer :: ireact
      integer :: isdiff
      integer :: isp
      integer :: ispc1, ispc2, ispc3
      integer :: jal
      integer :: jgas
      integer :: jnew
      integer :: jold      ! mappl(jold) for inactive spc
      integer :: jp, jpr
      integer :: jr, jre
      integer :: jspc, jspcl
      integer :: kdif
      integer :: knumporl
      integer :: kprods
      integer :: lfrac
      integer :: mc
      integer :: minvalu
      integer :: na, nar
      integer :: nallreac  ! total reactant positions in rxn
      integer :: nccount
      integer :: nfrcoun
      integer :: ngfsum, ngtsum, nltsum
      integer :: ngn, nmo, nol, npl
      integer :: ngsum, nlsum
      integer :: nk        ! rxn # of each rxn
      integer :: nk1, nkn
      integer :: nklast
      integer :: nmprod    ! max # of active products  in rxn
      integer :: nmreac    ! max # of active reactants in rxn
      integer :: nochang
      integer :: none
      integer :: npdcoun
      integer :: nplfun, npltot
      integer :: nprcoun
      integer :: nprodhi   ! nallreac+nmprod => highest product position in rxn
      integer :: nprodlo   ! nallreac+1      => lowest  product position in rxn
      integer :: nremain, nrept
      integer :: ntwo, nthree, nfour, nfive

      integer :: ngnfrac  (ICP)
      integer :: nolosrat (ICP)  ! # of active rxns with no loss terms

      integer :: isgaine  (ICS)
      integer :: isgainr  (ICS)  ! # of active spc with >= 1 net chem gain
      integer :: isporl   (ICS)  ! # of active spc with >= 1 net production or
                                 ! loss term
      integer :: nogaine  (ICS)  ! # of active spc with 0 net chem or emiss
                                 ! gains
      integer :: nspcsolv (ICS)  ! # of active spc with >= 1 gross loss
      integer :: nspec    (ICS)  ! total # of active spc

      integer :: newnk    (MAXGL)

      integer :: isaporl  (MXGSAER)  ! counts partial deriv terms for each spc

      integer :: nkphotrat(IPHOT, ICS)

      integer :: ignfrac  (MAXGL, ICS)
      integer :: nkgnfrac (MAXGL, ICS)

!     --------------------------------------------------------------
!     numgain : every net occurence of a production where spc is not
!               lost       in the same rxn (active & inactive spc)
!     numloss : every net occurence of a loss       where spc is not
!               reproduced in the same rxn (active & inactive spc)
!     --------------------------------------------------------------

      integer :: numgain  (MXGSAER, ICP)
      integer :: numloss  (MXGSAER, ICP)

      integer :: numporl  (MXGSAER, ICP)

      integer :: ifemis   (MXGSAER, ICS)
      integer :: iporl    (MXGSAER, ICS)  ! spc # of each isporl spc
      integer :: numgaint (MXGSAER, ICS)  ! every occurence of a production
                                          ! (active & inactive spc)
      integer :: numgfrt  (MXGSAER, ICS)
      integer :: numlost  (MXGSAER, ICS)  ! every occurence of a loss
                                          ! (active & inactive spc)

!     ---------------------------------------------------------------
!     iaprod : # of active products in each nk rxn; used to calculate
!              partial derivatives in Pderiv
!     ---------------------------------------------------------------

      integer :: iaprod   (NMTRATE, ICS)

      integer :: lskip    (NMTRATE, ICS)  ! 0 => rxn included,
                                          ! 1 => rxn not included
      integer :: ncequat  (NMTRATE, ICS)

      integer :: nkoner   (NMTRATE, ICS)  ! rxn # of each ioner rxn
      integer :: nktwor   (NMTRATE, ICS)  ! rxn # of each itwor rxn
      integer :: nkthrr   (NMTRATE, ICS)  ! rxn # of each ithrr rxn

      integer :: nolosrn  (NMTRATE, ICS)  ! rxn # of each rxn with no loss
                                          ! terms

!     ------------------------------------------------------------------
!     nrrep : 0 for each of two rxns where the reactants are identical;
!             more than two rxns have the same reactants, nrrep = 0 for
!             the first two rxns only (?)
!     nruse : 1,2,3 if rxn has 1, 2, or 3 active reactants, respectively
!     ------------------------------------------------------------------

      integer :: nrrep    (NMTRATE, ICS)
      integer :: nruse    (NMTRATE, ICS)

      integer :: isparder (MXGSAER, MXGSAER)

      integer :: jporl    (MXGSAER, MAXGL, ICS)

!     ------------------------------------------------------
!     irm : spc # of each reactant or product in each nk rxn
!     ------------------------------------------------------

      integer :: irm      (NMRPROD, NMTRATE, ICS)

      real*8  :: alfrac, rfrac
      real*8  :: diff

      real*8  :: aporl    (MXGSAER)

      real*8  :: fracp    (MAXGL,   ICS)

      real*8  :: fracgain (MXGSAER, ICS)

!     --------------------------------------------------------------------
!     fk2   : tbd
!     fkoef : 1, 2, fraction, or more => # of a given reactant or products
!             e.g., rxn           A + B  --> 2C + 0.34D  + D
!             value of fkoef      1   1      2    0.34     1
!     --------------------------------------------------------------------

      real*8  :: fk2  (NMRPROD, NMTRATE, ICS)
      real*8  :: fkoef(NMRPROD, NMTRATE, ICS)

      character (len=30) :: err_msg

      integer :: nm3bod (ICS)
      integer :: nmair  (ICS)
      integer :: nmn2   (ICS)
      integer :: nmo2   (ICS)

      integer :: lgas3bod(MAXGL3,  ICS)
      integer :: nreac3b (MAXGL3,  ICS)
      integer :: nreacair(MAXGL3,  ICS)
      integer :: nreacn2 (MAXGL3,  ICS)
      integer :: nreaco2 (MAXGL3,  ICS)


!     ----------------
!     Begin execution.
!     ----------------

      nmreac   = 3
      nallreac = 4
      nmprod   = 12
      nprodlo  = nallreac + 1
      nprodhi  = nallreac + nmprod

!     =======================
!#     include "setkin_smv2.h"
!     =======================

      !savedVars%lgas3bod(:,:) = LGAS3BOD(:,:)
      !savedVars%nreac3b (:,:) = NREAC3B (:,:)
      !savedVars%nreacair(:,:) = NREACAIR(:,:)
      !savedVars%nreaco2 (:,:) = NREACO2 (:,:)
      !savedVars%nreacn2 (:,:) = NREACN2 (:,:)

      !savedVars%nmn2  (:) = NMN2  (:)
      !savedVars%nmo2  (:) = NMO2  (:)
      !savedVars%nmair (:) = NMAIR (:)
      !savedVars%nm3bod(:) = NM3BOD(:)


      savedVars%nmoth (ncs) = 0
      nspec (ncs) = NNACT
      ntspec(ncs) = IGAS - 1

      nrates  (ncs) = ITHERM
      jphotrat(ncs) = IPHOT
      savedVars%ntrates (ncs) = nrates(ncs) + jphotrat(ncs)


      do i = 1, jphotrat(ncs)
        iplsithrm = ITHERM + i

        nkphotrat(i,ncs) = nrates(ncs) + i
        npphotrat(i,ncs) = i

        savedVars%jphotnk(iplsithrm,ncs) = i
      end do


      do i = 1, savedVars%ntrates(ncs)
        lskip(i,ncs) = 0
      end do


      ifnone = 0


      do j = 1, ICS
        do i = 1, NMTRATE

          nruse (i,j) = 0
          iaprod(i,j) = 0

        end do
      end do

      nrept = 0

      savedVars%ioner   (ncs) = 0
      savedVars%itwor   (ncs) = 0
      savedVars%ithrr   (ncs) = 0

      isgaine (ncs) = 0
      isgainr (ncs) = 0
      savedVars%ischang (ncs) = 0
      isporl  (ncs) = 0

      savedVars%nallrat (ncs) = 0
      ngnfrac (ncs) = 0
      nogaine (ncs) = 0
      savedVars%nolosp  (ncs) = 0
      nolosrat(ncs) = 0
      nspcsolv(ncs) = 0


      do i = 1, MXGSAER
        isaporl (i) = 0

        numloss (i,ncs) = 0
        numlost (i,ncs) = 0
        numgain (i,ncs) = 0
        numgaint(i,ncs) = 0
        numporl (i,ncs) = 0
        numgfrt (i,ncs) = 0

        fracgain(i,ncs) = 0.0d0
      end do


      do i = 1, NMTRATE
        i1 = NMTRATE + i

        iaprod  (i,ncs)  = 0
        nruse   (i,ncs)  = 0
        nrrep   (i,ncs)  = 0
        nolosrn (i,ncs)  = 0
        ncequat (i,ncs)  = 0
        savedVars%noldfnew(i,ncs)  = 0
        savedVars%newfold (i,ncs)  = 0
        nkoner  (i,ncs)  = 0
        nktwor  (i,ncs)  = 0
        nkthrr  (i,ncs)  = 0

        savedVars%newfold (i1,ncs) = 0
      end do


      do j = 1, MXGSAER
        do i = 1, MXGSAER
          isparder(i,j) = 0
        end do
      end do


      do nk = 1, savedVars%ntrates(ncs)

        if (lskip(nk,ncs) == 0) then

          savedVars%nallrat(ncs)     = savedVars%nallrat(ncs) + 1
          nar              = savedVars%nallrat(ncs)
          ncequat(nar,ncs) = nk

          if (nk <= nrates(ncs)) savedVars%nallrat(ncsp) = savedVars%nallrat(ncs)

          do k = 1, nmreac

            ireact = irm(k,nk,ncs)

            if ((ireact > 0) .and. (ireact <= nspec(ncs))) then

              do l = 1, nprodhi

                iporr = irm(l,nk,ncs)

                if (((l <= nmreac) .or. (l >= nprodlo)) .and.  &
     &              (iporr > 0)    .and.(iporr <= nspec(ncs))) then
                  isparder(iporr,ireact) = 1
                end if

              end do

            end if

          end do

        end if

      end do


      do ireact = 1, ntspec(ncs)
        do iporr = 1, ntspec(ncs)

          if (isparder(iporr,ireact) == 1) then
            isaporl(iporr) = isaporl(iporr) + 1
          end if

        end do
      end do


!     ------------------------------------------------------------
!     Rearrage species array so that all species with at least one
!     partial derivative term appear first, and those with zero
!     appear last.
!     ------------------------------------------------------------

      nochang = nspec(ncs)

      do jold = 1, ntspec(ncs)

        if (jold > nspec(ncs)) then
          savedVars%mappl(jold,ncs)      = jold
          inewold(jold,ncs)    = jold
        else if (isaporl(jold) > 0) then
          savedVars%ischang(ncs)         = savedVars%ischang(ncs) + 1
          jnew                 = savedVars%ischang(ncs)
          inewold(jnew,ncs)    = jold
          savedVars%mappl(jold,ncs)      = jnew
        else
          inewold(nochang,ncs) = jold
          savedVars%mappl(jold,ncs)      = nochang
          nochang              = nochang - 1
        end if

      end do


!     --------------------------------------------------------------------
!     Rearrage species in ischang array so that species with the fewest
!     partial derivative terms combined are placed first, and those with
!     the most appear last.  However, species with zero partial derivative
!     terms still appear after all ischang species.
!     --------------------------------------------------------------------

      do jnew = 1, savedVars%ischang(ncs)

        jold    = inewold(jnew,ncs)
        minvalu = isaporl(jold)
        iminold = jold
        iminnew = jnew

        do inew = jnew + 1, savedVars%ischang(ncs)

          iold = inewold(inew,ncs)

          if (isaporl(iold) < minvalu) then
            minvalu = isaporl(iold)
            iminold = iold
            iminnew = inew
          end if

        end do

        inewold(iminnew,ncs) = jold
        inewold(jnew,   ncs) = iminold

        savedVars%mappl(jold,   ncs)   = iminnew
        savedVars%mappl(iminold,ncs)   = jnew

      end do


!     -------------------------
!     Count gross and net loss.
!     -------------------------

      savedVars%nolosp(ncsp) = 0
      nklast       = 0

!     ===============================
      NKLOOP: do nk = 1, savedVars%ntrates(ncs)
!     ===============================

!       -----------------------------------------------------------------------
!       Determine occurrences of inactive species in rate equations.
!
!       Corrected treatment of nmoth (error in older version found by I. Plumb)
!       -----------------------------------------------------------------------

        do j = 1, nmreac

          ireact = irm(j,nk,ncs)

          if (ireact > 0) then

            ire = savedVars%mappl(ireact,ncs)

            if (ire > nspec(ncs)) then

              if (nk <= nrates(ncs)) then
                savedVars%nmoth(ncs)        = savedVars%nmoth(ncs) + 1
                nmo               = savedVars%nmoth(ncs)
                savedVars%nreacoth(nmo,ncs) = nk
                savedVars%lgasbino(nmo,ncs) = ireact
              else
                savedVars%nolosp(ncs)       = savedVars%nolosp(ncs) + 1
                nol               = savedVars%nolosp(ncs)
                savedVars%nknlosp (nol,ncs) = nk
                savedVars%losinacp(nol,ncs) = ireact
              end if

            end if

          end if

        end do

!       ================================
        LSKIPIF: if (lskip(nk,ncs) == 0) then
!       ================================

          ial = 0

          do jspc = 1, ntspec(ncs)
            aporl(jspc) = 0.0d0
          end do

!         ------------------------------------------
!         Set array to identify active loss species.
!         ------------------------------------------

          do j = 1, nmreac

            ireact = irm(j,nk,ncs)

            if (ireact > 0) then

              ire              = savedVars%mappl(ireact,ncs)
              aporl(ire)       = aporl(ire) - 1.0d0
              numlost(ire,ncs) = numlost(ire,ncs) + 1

              if (ire <= nspec(ncs)) then
                ial              = ial + 1
                savedVars%irm2(ial,nk,ncs) = ire
              end if

            end if

          end do

!         ---------------------------------------------------------------
!         Set arrays to identify reactions with at least one active loss.
!         ---------------------------------------------------------------

!         ===================
          IALIF: if (ial > 0) then
!         ===================

            nruse(nk,ncs) = ial
            nrrep(nk,ncs) = ial

            if (ial == 1) then
              savedVars%ioner(ncs)             = savedVars%ioner(ncs) + 1
              nkoner(savedVars%ioner(ncs),ncs) = nk
            else if (ial == 2) then
              savedVars%itwor(ncs)             = savedVars%itwor(ncs) + 1
              nktwor(savedVars%itwor(ncs),ncs) = nk
            else if (ial == 3) then
              savedVars%ithrr(ncs)             = savedVars%ithrr(ncs) + 1
              nkthrr(savedVars%ithrr(ncs),ncs) = nk
            end if

!           ----------------------------------------------------------
!           Compare two consecutive reactions; if the species (but not
!           rates) are the same, then save multiplications in Subfun.
!           ----------------------------------------------------------

            if (nklast > 0) then

              if (nruse(nklast,ncs) == ial) then

                isdiff = 0

                do ib = 1, ial
                  jspcl = savedVars%irm2(ib,nklast,ncs)
                  jspc  = savedVars%irm2(ib,nk    ,ncs)
                  if (jspcl /= jspc) isdiff = 1
                end do

                if ((isdiff == 0).and. (nrrep(nklast,ncs) /= 0)) then
                  nrrep(nk,ncs)     = 0
                  nrrep(nklast,ncs) = 0

                  nrept = nrept + 1
                  ispc1 = savedVars%irm2(1,nk,ncs)
                  ispc2 = savedVars%irm2(2,nk,ncs)
                  ispc3 = savedVars%irm2(3,nk,ncs)

                  if (ispc1 > 0) ispc1 = inewold(ispc1,ncs)
                  if (ispc2 > 0) ispc2 = inewold(ispc2,ncs)
                  if (ispc3 > 0) ispc3 = inewold(ispc3,ncs)

                  if (pr_smv2) then
                    Write (lunsmv,800)  &
     &                nrept, nk, namesp2(ispc1,ncs),  &
     &                namesp2(ispc2,ncs), namesp2(ispc3,ncs)
                  end if

 800              format ('Repeat reactants: ', i5, i5, 3(1x,a14))

                end if

              end if

            end if

!         --------------------------------------------------------------
!         Determine the number of reactions with zero active loss terms.
!         --------------------------------------------------------------

!         ==================
          else if (ial == 0) then
!         ==================

            nolosrat(ncs)     = nolosrat(ncs) + 1
            nol               = nolosrat(ncs)
            nolosrn (nol,ncs) = nk

!         ============
          end if IALIF
!         ============

!         ------------------------------------------------------------------
!         Count gross and net production and set a partial derivative array.
!         ------------------------------------------------------------------

          iap = nprodlo - 1

          do k = nprodlo, nprodhi

            iprod = irm(k,nk,ncs)

            if (iprod > 0) then

              ipr    = savedVars%mappl(iprod,ncs)
              rfrac  = fkoef(k,nk,ncs)
              lfrac  = rfrac + SMAL1
              alfrac = lfrac
              diff   = Abs (rfrac-alfrac)

              if (diff > SMAL1) then

!               -- Production term is a fraction. --

                if (ipr <= nspec(ncs)) then
                  ngnfrac (ncs)     = ngnfrac(ncs) + 1
                  ngn               = ngnfrac(ncs)
                  ignfrac (ngn,ncs) = ipr
                  nkgnfrac(ngn,ncs) = nk
                  fracp   (ngn,ncs) = rfrac
                end if

                kprods            = 1
                numgfrt (ipr,ncs) = numgfrt (ipr,ncs) + 1
                fracgain(ipr,ncs) = fracgain(ipr,ncs) + rfrac

              else

!               -- Production term is non-fraction. --

                aporl(ipr)        = aporl(ipr) + rfrac
                kprods            = lfrac
                numgaint(ipr,ncs) = numgaint(ipr,ncs) + lfrac
                fkoef(k,nk,ncs)   = 1.0d0

              end if

              if (ipr <= nspec(ncs)) then

!               -- Identify all production terms. --

                do l = 1, kprods
                  iap              = iap + 1
                  iaprod  (nk,ncs) = iap
                  savedVars%irm2(iap,nk,ncs) = ipr
                  fk2 (iap,nk,ncs) = fkoef(k,nk,ncs)
                end do

              end if

            end if

          end do

!         ---------------------------------------------------------------
!         Find net prod and loss terms for all but fractionated products.
!         ---------------------------------------------------------------

          do jspc = 1, ntspec(ncs)

            if (Abs (aporl(jspc)) < SMAL1) then

              kdif = 0

            else if (aporl(jspc) > 0.0d0) then

              kdif = aporl(jspc) + 0.00001d0

              do l = 1, kdif
                numgain(jspc,ncs)   = numgain(jspc,ncs) + 1
                numporl(jspc,ncs)   = numporl(jspc,ncs) + 1
                npl                 = numporl(jspc,ncs)
                jporl(jspc,npl,ncs) = nk + savedVars%ntrates(ncs)
              end do

            else

              kdif = aporl(jspc) - 0.00001d0
              kdif = -kdif

              do l = 1, kdif
                numloss(jspc,ncs)   = numloss(jspc,ncs) + 1
                numporl(jspc,ncs)   = numporl(jspc,ncs) + 1
                npl                 = numporl(jspc,ncs)
                jporl(jspc,npl,ncs) = nk
              end do

            end if

            if (nk <= nrates(ncs)) then
              numloss(jspc,ncsp) = numloss(jspc,ncs)
              numgain(jspc,ncsp) = numgain(jspc,ncs)
              numporl(jspc,ncsp) = numporl(jspc,ncs)
            end if

          end do

          if (nk <= nrates(ncs)) then
            nolosrat(ncsp) = nolosrat(ncs)
            ngnfrac (ncsp) = ngnfrac (ncs)
            savedVars%ioner   (ncsp) = savedVars%ioner   (ncs)
          end if

          nklast = nk

!       ==============
        end if LSKIPIF
!       ==============

!     =============
      end do NKLOOP
!     =============


!     --------------------------------------------------------------
!     Set array for reordering rates from 3..2..1..0 body reactions.
!     --------------------------------------------------------------

      ic = 0

      do i = 1, savedVars%ithrr(ncs)
        ic  = ic + 1
        nk  = nkthrr(i,ncs)
        nk1 = nk + savedVars%ntrates(ncs)

        savedVars%noldfnew(ic, ncs) = nk
        savedVars%newfold (nk, ncs) = ic
        savedVars%newfold (nk1,ncs) = ic + savedVars%nallrat(ncs)
      end do

      ntwo = savedVars%ithrr(ncs) + savedVars%itwor(ncs)
      icb  = ntwo + 1

      do i = 1, savedVars%itwor(ncs)

        nk  = nktwor(i,ncs)
        nk1 = nk + savedVars%ntrates(ncs)

        if (nrrep(nk,ncs) > 0) then
         ic  = ic + 1
         icd = ic
        else
         icb = icb - 1
         icd = icb
        end if

        savedVars%noldfnew(icd,ncs) = nk
        savedVars%newfold (nk, ncs) = icd
        savedVars%newfold (nk1,ncs) = icd + savedVars%nallrat(ncs)

      end do

      savedVars%inorep(ncs) = ic
      ic          = ntwo

      do i = 1, savedVars%ioner(ncs)
        ic  = ic + 1
        nk  = nkoner(i,ncs)
        nk1 = nk + savedVars%ntrates(ncs)

        savedVars%noldfnew(ic, ncs) = nk
        savedVars%newfold (nk, ncs) = ic
        savedVars%newfold (nk1,ncs) = ic + savedVars%nallrat(ncs)
      end do

      do i = 1, nolosrat(ncs)
        ic  = ic + 1
        nk  = nolosrn(i,ncs)
        nk1 = nk + savedVars%ntrates(ncs)

        savedVars%noldfnew(ic, ncs) = nk
        savedVars%newfold (nk, ncs) = ic
        savedVars%newfold (nk1,ncs) = ic + savedVars%nallrat(ncs)
      end do

      if (ic /= savedVars%nallrat(ncs)) then
        Write (6,810) ic, savedVars%nallrat(ncs)
        err_msg = ' Problem in Jsparse '
        call GmiPrintError (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
      end if

 810  format ('Jsparse: ic /= nallrat = ', 2(i5))


!     ------------------------------------------
!     Set a slightly more efficient photo array.
!     ------------------------------------------

      do j = 1, jphotrat(ncs)
        nk               = nkphotrat(j,ncs)
        nkn              = savedVars%newfold(nk,ncs)
        savedVars%nknphotrt(j,ncs) = nkn
      end do

!     --------------------------------------------------------
!     Determine number of species with gross/net losses/gains.
!     --------------------------------------------------------

      do jold = 1, nspec(ncs)

        jnew = savedVars%mappl(jold,ncs)

        if (numgain(jnew,ncs) > 0) then
         isgainr(ncs) = isgainr(ncs) + 1
        end if

        if (numporl(jnew,ncs) > 0) then
         isporl(ncs)     = isporl(ncs) + 1
         isp             = isporl(ncs)
         iporl (isp,ncs) = jnew
        end if

        if (numlost(jnew,ncs) > 0) then
          nspcsolv(ncs) = nspcsolv(ncs) + 1
        end if

        ifemis(jold, ncs) = 0

        if ((numgain (jnew,ncs) >  0)     .or.  &
     &      (fracgain(jnew,ncs) >  0.0d0) .or.  &
     &      (ifemis  (jold,ncs) == 1)) then
          isgaine(ncs) = isgaine(ncs) + 1
        else if (numloss(jnew,ncs) > 0) then
          nogaine(ncs) = nogaine(ncs) + 1
        end if

      end do


!     -------------------------------------------------
!     Check dimensions resulting from gains and losses.
!     -------------------------------------------------

      ngtsum = 0
      nltsum = 0
      ngsum  = 0
      nlsum  = 0
      ngfsum = 0

      do k = 1, ntspec(ncs)

        j = inewold(k,ncs)

        ngtsum  = ngtsum + numgaint(k,ncs)
        nltsum  = nltsum + numlost (k,ncs)
        ngsum   = ngsum  + numgain (k,ncs)
        nlsum   = nlsum  + numloss (k,ncs)
        ngfsum  = ngfsum + numgfrt (k,ncs)

        if ((numgaint(k,ncs) > MAXGL) .or.  &
     &      (numlost (k,ncs) > MAXGL)) then
          Write (6,820) namesp2(j,ncs), numgaint(k,ncs), numlost(k,ncs)
          err_msg = ' Problem in Jsparse '
          call GmiPrintError (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
        end if

      end do

      if ((savedVars%nmoth  (ncs) > MAXGL2) .or. (savedVars%nolosp(ncs) > MAXGL3) .or.  &
     &    (ngnfrac(ncs) > MAXGL)) then
        Write (6,830) MAXGL2, savedVars%nmoth(ncs), MAXGL3, savedVars%nolosp(ncs),  &
     &                MAXGL, ngnfrac(ncs)
        err_msg = ' Problem in Jsparse '
        call GmiPrintError (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
      end if

      if (pr_smv2) then
        Write (lunsmv,840)  &
     &    savedVars%nallrat(ncsp), savedVars%nallrat(ncs)-savedVars%nallrat(ncsp), &
     &    savedVars%nallrat(ncs)
      end if

 820  format  &
     &  ('Gearset: spec ', a6,' dimens exceeded either numgaint ',  &
     &   'numloss, numgain, or numlost > MAXGL ', 4(i3,1x))

 830  format  &
     &  ('One of the dimensions below is too small:', /,  &
     &   'Dimension: MAXGL2 = ', i4, ' Variable: nmoth   = ', i4, /  &
     &   'Dimension: MAXGL3 = ', i4, ' Variable: nolosp  = ', i4, /  &
     &   'Dimension: MAXGL  = ', i4, ' Variable: ngnfrac = ', i4)

 840  format  &
     &  (/, '# Kinetic reactions: ', i5 ,' Photorates: ', i5,  &
     &   ' Total: ', i5)

 881  format('#irm2(',i3,',',i3,') = ',i5,i7,i7,f12.8,f14.10)

!     ------------------------------------------------
!     Set arrays to take advantage of sparse matrices.
!     ------------------------------------------------

      if (ifnone == 0) then
        ifnone  = 1
        nplfun  = 0
        nfrcoun = 0
        npdcoun = 0
        nprcoun = 0
        npltot  = 0
      end if


!     ==========================
      IFSUNLOOP: do ifsun = 1, 2
!     ==========================

        ncsp = ((ifsun - 1) * ICS) + ncs

        do i = 1, MXGSAER
          do j= 1, MXGSAER
            lzero(j,i) = 0
          end do
          lzero(i,i) = 1
        end do

        do na = 1, savedVars%nallrat(ncsp)

          nk      = ncequat(na,ncs)
          ihireac = nruse  (nk,ncs)

          do ial = 1, ihireac

            ire = savedVars%irm2(ial,nk,ncs)

            do jal = 1, ihireac
              jre            = savedVars%irm2(jal,nk,ncs)
              lzero(jre,ire) = 1
            end do

            do iap = nprodlo, iaprod(nk,ncs)
              jpr            = savedVars%irm2(iap,nk,ncs)
              lzero(jpr,ire) = 1
            end do

          end do

        end do

!       -------------------------------------------------------------
!       Set decomposition and back-substitution sparse matrix arrays.
!       -------------------------------------------------------------

!       ============
        call KsparseGMI  &
!       ============
     &    (savedVars, pr_smv2, lunsmv, ncs, ncsp, lzero, jarraypt, &
                          iccount, jccount, kccount, mccount, icnt,     &
                          jcnt, kcnt, mcnt, idecomp, ijtot, kztot,      &
                          mztot, kbsub, mbsub)

!       -----------------------------------------------------------
!       Set arrays to improve efficiency of first-derivative calcs.
!
!       Set arrays for kinetic and photo production and loss rates.
!       -----------------------------------------------------------

        savedVars%npllo(ncsp) = npltot + 1

        do i = 1, isporl(ncs)

          jspc     = iporl(i,ncs)
          knumporl = numporl(jspc,ncsp)
          nccount  = 0
          npltot   = npltot + 1
          nremain  = knumporl
          nfive    = (nremain + 0.0001d0) / 5
          nremain  =  nremain - nfive   * 5
          nfour    = (nremain + 0.0001d0) / 4
          nremain  =  nremain - nfour   * 4
          nthree   = (nremain + 0.0001d0) / 3
          nremain  =  nremain - nthree  * 3
          ntwo     = (nremain + 0.0001d0) / 2
          nremain  =  nremain - ntwo    * 2
          none     = (nremain + 0.0001d0)
          nremain  =  nremain - none

          savedVars%jspnpl(npltot) = jspc
          savedVars%npl5  (npltot) = nplfun       + 1
          savedVars%nph5  (npltot) = nplfun       + nfive
          savedVars%npl4  (npltot) = savedVars%nph5(npltot) + 1
          savedVars%nph4  (npltot) = savedVars%nph5(npltot) + nfour
          savedVars%npl3  (npltot) = savedVars%nph4(npltot) + 1
          savedVars%nph3  (npltot) = savedVars%nph4(npltot) + nthree
          savedVars%npl2  (npltot) = savedVars%nph3(npltot) + 1
          savedVars%nph2  (npltot) = savedVars%nph3(npltot) + ntwo
          savedVars%npl1  (npltot) = savedVars%nph2(npltot) + 1
          savedVars%nph1  (npltot) = savedVars%nph2(npltot) + none
          nplfun         = savedVars%nph1(npltot)

          do n = 1, knumporl
            nk       = jporl(jspc,n,ncs)
            newnk(n) = savedVars%newfold(nk,ncs)
          end do

          do mc = savedVars%npl5(npltot), savedVars%nph5(npltot)
            savedVars%lossra(mc) = newnk(nccount+1)
            savedVars%lossrb(mc) = newnk(nccount+2)
            savedVars%lossrc(mc) = newnk(nccount+3)
            savedVars%lossrd(mc) = newnk(nccount+4)
            savedVars%lossre(mc) = newnk(nccount+5)
            nccount    = nccount + 5
          end do

          do mc = savedVars%npl4(npltot), savedVars%nph4(npltot)
            savedVars%lossra(mc) = newnk(nccount+1)
            savedVars%lossrb(mc) = newnk(nccount+2)
            savedVars%lossrc(mc) = newnk(nccount+3)
            savedVars%lossrd(mc) = newnk(nccount+4)
            nccount    = nccount + 4
          end do

          do mc = savedVars%npl3(npltot), savedVars%nph3(npltot)
            savedVars%lossra(mc) = newnk(nccount+1)
            savedVars%lossrb(mc) = newnk(nccount+2)
            savedVars%lossrc(mc) = newnk(nccount+3)
            nccount    = nccount + 3
          end do

          do mc = savedVars%npl2(npltot), savedVars%nph2(npltot)
            savedVars%lossra(mc) = newnk(nccount+1)
            savedVars%lossrb(mc) = newnk(nccount+2)
            nccount    = nccount + 2
          end do

          do mc = savedVars%npl1(npltot), savedVars%nph1(npltot)
            savedVars%lossra(mc) = newnk(nccount+1)
            nccount    = nccount + 1
          end do

        end do

        savedVars%nplhi(ncsp) = npltot

!       ------------------------------------
!       Set array for fractionated products.
!       ------------------------------------

        savedVars%nfrlo(ncsp) = nfrcoun + 1

        do i = 1, ngnfrac(ncsp)
          jspc             = ignfrac(i,ncs)
          nfrcoun          = nfrcoun + 1
          savedVars%jspcnfr(nfrcoun) = jspc
          nk               = nkgnfrac(i,ncs)
          savedVars%nknfr  (nfrcoun) = savedVars%newfold(nk,ncs)
          savedVars%fracnfr(nfrcoun) = fracp(i,ncs)
        end do

        savedVars%nfrhi(ncsp) = nfrcoun

!       -------------------------------------------------------------
!       Set arrays to improve efficiency of partial derivative calcs.
!       -------------------------------------------------------------

        savedVars%npdlo(ncsp) = npdcoun + 1

        do na = 1, savedVars%nallrat(ncsp)

          nk      = ncequat(na,ncs)
          ihireac = nruse  (nk,ncs)

          do ial = 1, ihireac

            ir = savedVars%irm2(ial,nk,ncs)

            do jal = 1, ihireac
              jr      = savedVars%irm2(jal,nk,ncs)
              iar     = jarraypt(jr,ir)
              npdcoun = npdcoun + 1

              savedVars%nkpdterm(npdcoun) = savedVars%newfold(nk,ncs)
              savedVars%ipospd  (npdcoun) = iar
              savedVars%iialpd  (npdcoun) = ial
              savedVars%fracpl  (npdcoun) = -1.0d0
            end do

            do iap = nprodlo, iaprod(nk,ncs)
              jp      = savedVars%irm2(iap,nk,ncs)
              iar     = jarraypt(jp,ir)
              npdcoun = npdcoun + 1

              savedVars%nkpdterm(npdcoun) = savedVars%newfold(nk,ncs)
              savedVars%ipospd  (npdcoun) = iar
              savedVars%iialpd  (npdcoun) = ial
              savedVars%fracpl  (npdcoun) = fk2(iap,nk,ncs)
            end do

          end do

        end do

        savedVars%npdhi(ncsp) = npdcoun

!       ---------------------------------------------
!       Check dimensions and print out array savings.
!       ---------------------------------------------

        if ((npltot  > MXCOUNT4)  .or. (nplfun  > MXCOUNT4) .or.  &
     &      (nprcoun > MXCOUNT4)  .or. (nfrcoun > MXCOUNT4) .or.  &
     &      (npdcoun > MXCOUNT2)) then

          Write (6,890) MXCOUNT4, npltot,  MXCOUNT4, nplfun,  &
     &                  MXCOUNT4, nprcoun, MXCOUNT4, nfrcoun,  &
     &                  MXCOUNT2, npdcoun

          err_msg = ' Problem in Jsparse '
          call GmiPrintError (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)

        end if

 890    format  &
     &    ('One of the dimensions below is too small:', /,  &
     &     'Dimension: MXCOUNT4 = ', i4, ' Variable: npltot  = ', i4, /,  &
     &     'Dimension: MXCOUNT4 = ', i4, ' Variable: nplfun  = ', i4, /,  &
     &     'Dimension: MXCOUNT4 = ', i4, ' Variable: nprcoun = ', i4, /,  &
     &     'Dimension: MXCOUNT4 = ', i4, ' Variable: nfrcoun = ', i4, /,  &
     &     'Dimension: MXCOUNT2 = ', i4, ' Variable: npdcoun = ', i4)

!     ================
      end do IFSUNLOOP
!     ================

 891  format('#lzero(',i3,',',i3,') = ',i5,i7,i7)
 892  format('#lzero(',i3,',',i3,') = ',i5,i7,i7,i7)

      return

      end subroutine JsparseGMI


!-----------------------------------------------------------------------------
!
! ROUTINE
!   KsparseGMI
!
! DESCRIPTION
!   This routine sets up sparse matrix and other arrays.  It also sets arrays
!   for gas-phase, aqueous-phase, or any other type of chemistry.  Finally,
!   it sets arrays for both day and night chemistry of each type.
!
!   Sets up arrays for decomposition / back-substitution of sparse matrices
!   by removing all calculations involving a zero.  Sets other arrays to take
!   advantage of sparse matrices.
!
! ARGUMENTS
!   pr_smv2  : should the SmvgearII output file be written
!              (non-parallel mode only)?
!   lunsmv   : logical unit number to write to when pr_smv2 is true
!   ncs      : identifies gas chemistry type (1..NCSGAS)
!   ncsp     : ncs       => for daytime   gas chemistry
!              ncs + ICS => for nighttime gas chemistry
!   lzero    : = 1 if an array spot is filled with a non-zero value; it is
!              updated as we simulate the order of calculations during a
!              practice l-u decomposition (?)
!   jarraypt : identifies 1d array point for each 2d point i,j
!
!-----------------------------------------------------------------------------

      subroutine KsparseGMI(savedVars, pr_smv2, lunsmv, ncs, ncsp, lzero, &
                          jarraypt, iccount, jccount, kccount, &
                          mccount, icnt, jcnt, kcnt, mcnt, idecomp, &
                          ijtot, kztot, mztot, kbsub, mbsub)

      use GmiPrintError_mod, only : GmiPrintError
      use GmiSavedVariables_mod, only : t_ChemistrySaved

      implicit none

#     include "smv2chem_par.h"
!#     include "smv2chem2.h"

!     ----------------------
!     Argument declarations.
!     ----------------------

      logical, intent(in)  :: pr_smv2
      integer, intent(in)  :: lunsmv
      integer, intent(in)  :: ncs
      integer, intent(in)  :: ncsp

      integer, intent(inout) :: lzero (MXGSAER, MXGSAER)

      integer, intent(out) :: jarraypt(MXGSAER, MXGSAER)

      integer :: iccount, jccount, kccount, mccount
      integer :: icnt, jcnt, kcnt, mcnt
      integer :: idecomp
      integer :: ijtot, kztot, mztot
      integer :: kbsub, mbsub

      type(t_ChemistrySaved), intent(inOut) :: savedVars

!     ----------------------
!     Variable declarations.
!     ----------------------

      integer :: i, j, k
      integer :: i1, i2
      integer :: ia, ic

!     -----------------------------------------------------------------------
!     icnta    : # operations in Decomp  loop 1 w/o  sparse matrix reductions
!     jcnta    : # operations in Decomp  loop 2 w/o  sparse matrix reductions
!     kcnta    : # operations in Backsub loop 1 w/o  sparse matrix reductions
!     mcnta    : # operations in Backsub loop 2 w/o  sparse matrix reductions
!     kount0a  : # initial matrix spots  filled w/o  sparse matrix reductions
!     kntarray : # final   matrix spots  filled w/o  sparse matrix reductions
!     -----------------------------------------------------------------------

      integer :: icnta, jcnta, kcnta, mcnta
      integer :: kount0a, kntarray

!     -----------------------------------------------------------------------
!     icntb    : # operations in Decomp  loop 1 with sparse matrix reductions
!     jcntb    : # operations in Decomp  loop 2 with sparse matrix reductions
!     kcntb    : # operations in Backsub loop 1 with sparse matrix reductions
!     mcntb    : # operations in Backsub loop 2 with sparse matrix reductions
!     kount0   : # initial matrix spots  filled with sparse matrix reductions
!     iarray2  : # final   matrix spots  filled with sparse matrix reductions
!     -----------------------------------------------------------------------

      integer :: icntb, jcntb, kcntb, mcntb
      integer :: kount0, iarray2

      integer :: izil
      integer :: j1
      integer :: ka, kb, kc, kd, ke
      integer :: kzil
      integer :: mc
      integer :: mzil
      integer :: none

      integer :: nremain
      integer :: ntwo, nthree, nfour, nfive

!     ------------------------------------------------------------------
!     izerok : used to find 1d array point for each i,k value found in
!              same loop; each occurrence of each izilch calculation (?)
!     ------------------------------------------------------------------

      integer :: izerok(MXCOUNT2)

      integer :: jzero (MXCOUNT3)
      integer :: kzero (MXCOUNT3)
      integer :: mzero (MXCOUNT3)

!     ---------------------------------------------------------------------
!     jzilch : contains non-zeros for lower triangular matrix;
!              # of calculations with non-zero values to fill lower part of
!              decomposed matrix (?)
!     kzilch : contains non-zeros for lower triangular matrix
!     mzilch : contains non-zeros for upper triangular matrix, where back-
!              substitution occurs
!     ---------------------------------------------------------------------

      integer :: jzilch(MXGSAER)
      integer :: kzilch(MXGSAER)
      integer :: mzilch(MXGSAER)

!     --------------------------------------------------------------------
!     izilch : # of calculations with non-zero values during matrix decomp
!     --------------------------------------------------------------------

      integer :: izilch(MXGSAER, MXGSAER)


!     ---------------------------------------------------------------
!     pertst : coefficients used to select step-size and order; thus,
!              only about 1% accuracy needed; see Gear (1971) or
!              Hindmarsh (1973, UCID-30059).
!     ---------------------------------------------------------------

      character(len=30) :: err_msg

!     ----------------
!     Begin execution.
!     ----------------
 
      izerok = 0
      jzero  = 0
      kzero  = 0
      mzero  = 0

      kount0a = 0
      kount0  = 0
      icnta   = 0
      icntb   = 0
      jcnta   = 0
      jcntb   = 0
      kcnta   = 0
      kcntb   = 0
      mcnta   = 0
      mcntb   = 0
      iarray2 = 0


      do k = 1, savedVars%ischang(ncs)
        do j = 1, savedVars%ischang(ncs)

          kount0a = kount0a + 1
          if (lzero(j,k) == 1) kount0 = kount0 + 1

        end do
      end do


      do j = 1, savedVars%ischang(ncs)
        do k = 1, savedVars%ischang(ncs)
          jarraypt(k,j) = 0
        end do
      end do


!     -------------------------
!     Arrays for decomposition.
!     -------------------------

      do j = 1, savedVars%ischang(ncs)

        jzilch(j) = 0
        j1        = j - 1

!    -- First loop of decompostion. --

        do i = 2, savedVars%ischang(ncs)

          izilch(j,i) = 0
          i1          = j1
          if (i <= j1) i1 = i - 1

          do k = 1, i1
            icnta = icnta + 1
            if ((lzero(i,k) == 1) .and.(lzero(k,j) == 1)) then
              izilch(j,i)  = izilch(j,i) + 1
              icnt         = icnt  + 1
              icntb        = icntb + 1
              izerok(icnt) = k
              lzero(i,j)   = 1
            end if
          end do

        end do

!       -- Second loop of decompostion. --

        do i = j + 1, savedVars%ischang(ncs)

          jcnta = jcnta + 1

          if (lzero(i,j) == 1) then
            jzilch(j)   = jzilch(j) + 1
            jcnt        = jcnt  + 1
            jcntb       = jcntb + 1
            jzero(jcnt) = i
          end if

        end do

      end do


!     -----------------------------
!     Arrays for back-substitution.
!     -----------------------------

!     -- First loop of back-substitution. --

      do i = 2, savedVars%ischang(ncs)

        kzilch(i) = 0
        i1        = I - 1

        do j = 1, i1

          kcnta = kcnta + 1

          if (lzero(i,j) == 1) then
            kzilch(i)         = kzilch(i) + 1
            kcntb             = kcntb     + 1
            kcnt              = kcnt      + 1
            iarray2           = iarray2   + 1
            kzero(kcnt)       = j
            jarraypt(i,j)     = iarray2
          end if

        end do

      end do

!     -- Second loop of back-substitution. --

      do i = savedVars%ischang(ncs), 1, -1

        mzilch(i) = 0
        i2        = i + 1

        do j = i2, savedVars%ischang(ncs)

          mcnta = mcnta + 1

          if (lzero(i,j) == 1) then
            mzilch(i)     = mzilch(i) + 1
            mcntb         = mcntb     + 1
            mcnt          = mcnt      + 1
            iarray2       = iarray2   + 1
            mzero(mcnt)   = j
            jarraypt(i,j) = iarray2
          end if

        end do

      end do

 580  format('+lzero(',i3,',',i3,') = ',i5)
 581  format('-jzero(',i3,') = ',i5,i7)
 582  format('*kzero(',i3,') = ',i5,i7)


!     -----------------------------------------------------------
!     Fill jarraypt with remaining array points (along diagonal).
!     -----------------------------------------------------------

      do i = 1, savedVars%ischang(ncs)
        iarray2       = iarray2 + 1
        jarraypt(i,i) = iarray2
      end do

      savedVars%iarray(ncsp) = iarray2
      kntarray     = kcnta + mcnta + savedVars%ischang(ncs)


!     ----------------------------------------------------------
!     Change izero and jzero arrays so their values point to new
!     array positions defined in jarraypt.
!     ----------------------------------------------------------

!     ==============================
      JLOOP1: do j = 1, savedVars%ischang(ncs)
!     ==============================

!       -- First loop of decompostion. --

        savedVars%ijtlo(j,ncsp) = ijtot + 1

!       ==============================
        ILOOP1: do i = 2, savedVars%ischang(ncs)
!       ==============================

          izil = izilch(j,i)

          if (izil > 0) then

            ijtot   = ijtot + 1
            nremain = izil
            nfive   = (nremain + 0.0001d0) / 5
            nremain =  nremain - (nfive  * 5)
            nfour   = (nremain + 0.0001d0) / 4
            nremain =  nremain - (nfour  * 4)
            nthree  = (nremain + 0.0001d0) / 3
            nremain =  nremain - (nthree * 3)
            ntwo    = (nremain + 0.0001d0) / 2
            nremain =  nremain - (ntwo   * 2)
            none    = (nremain + 0.0001d0)
            nremain =  nremain - none

            savedVars%ijval(ijtot) = jarraypt(i,j)
            savedVars%idl5 (ijtot) = idecomp     + 1
            savedVars%idh5 (ijtot) = idecomp     + nfive
            savedVars%idl4 (ijtot) = savedVars%idh5(ijtot) + 1
            savedVars%idh4 (ijtot) = savedVars%idh5(ijtot) + nfour
            savedVars%idl3 (ijtot) = savedVars%idh4(ijtot) + 1
            savedVars%idh3 (ijtot) = savedVars%idh4(ijtot) + nthree
            savedVars%idl2 (ijtot) = savedVars%idh3(ijtot) + 1
            savedVars%idh2 (ijtot) = savedVars%idh3(ijtot) + ntwo
            savedVars%idl1 (ijtot) = savedVars%idh2(ijtot) + 1
            savedVars%idh1 (ijtot) = savedVars%idh2(ijtot) + none
            idecomp      = savedVars%idh1(ijtot)

            do ic = savedVars%idl5(ijtot), savedVars%idh5(ijtot)
              ka         = izerok(iccount+1)
              kb         = izerok(iccount+2)
              kc         = izerok(iccount+3)
              kd         = izerok(iccount+4)
              ke         = izerok(iccount+5)
              iccount    = iccount + 5
              savedVars%ikdeca(ic) = jarraypt(i,ka)
              savedVars%ikdecb(ic) = jarraypt(i,kb)
              savedVars%ikdecc(ic) = jarraypt(i,kc)
              savedVars%ikdecd(ic) = jarraypt(i,kd)
              savedVars%ikdece(ic) = jarraypt(i,ke)
              savedVars%kjdeca(ic) = jarraypt(ka,j)
              savedVars%kjdecb(ic) = jarraypt(kb,j)
              savedVars%kjdecc(ic) = jarraypt(kc,j)
              savedVars%kjdecd(ic) = jarraypt(kd,j)
              savedVars%kjdece(ic) = jarraypt(ke,j)
            end do

            do ic = savedVars%idh5(ijtot) + 1, savedVars%idh4(ijtot)
              ka         = izerok(iccount+1)
              kb         = izerok(iccount+2)
              kc         = izerok(iccount+3)
              kd         = izerok(iccount+4)
              iccount    = iccount + 4
              savedVars%ikdeca(ic) = jarraypt(i,ka)
              savedVars%ikdecb(ic) = jarraypt(i,kb)
              savedVars%ikdecc(ic) = jarraypt(i,kc)
              savedVars%ikdecd(ic) = jarraypt(i,kd)
              savedVars%kjdeca(ic) = jarraypt(ka,j)
              savedVars%kjdecb(ic) = jarraypt(kb,j)
              savedVars%kjdecc(ic) = jarraypt(kc,j)
              savedVars%kjdecd(ic) = jarraypt(kd,j)
            end do

            do ic = savedVars%idh4(ijtot) + 1, savedVars%idh3(ijtot)
              ka         = izerok(iccount+1)
              kb         = izerok(iccount+2)
              kc         = izerok(iccount+3)
              iccount    = iccount + 3
              savedVars%ikdeca(ic) = jarraypt(i,ka)
              savedVars%ikdecb(ic) = jarraypt(i,kb)
              savedVars%ikdecc(ic) = jarraypt(i,kc)
              savedVars%kjdeca(ic) = jarraypt(ka,j)
              savedVars%kjdecb(ic) = jarraypt(kb,j)
              savedVars%kjdecc(ic) = jarraypt(kc,j)
            end do

            do ic = savedVars%idh3(ijtot) + 1, savedVars%idh2(ijtot)
              ka         = izerok(iccount+1)
              kb         = izerok(iccount+2)
              iccount    = iccount + 2
              savedVars%ikdeca(ic) = jarraypt(i,ka)
              savedVars%ikdecb(ic) = jarraypt(i,kb)
              savedVars%kjdeca(ic) = jarraypt(ka,j)
              savedVars%kjdecb(ic) = jarraypt(kb,j)
            end do

            do ic = savedVars%idh2(ijtot) + 1, savedVars%idh1(ijtot)
              ka         = izerok(iccount+1)
              iccount    = iccount + 1
              savedVars%ikdeca(ic) = jarraypt(i,ka)
              savedVars%kjdeca(ic) = jarraypt(ka,j)
            end do

          end if

!       =============
        end do ILOOP1
!       =============

        savedVars%ijthi(j,ncsp) = ijtot

        savedVars%jarrdiag(j,ncsp) = jarraypt(j,j)

!       -- Second loop of decompostion. --

        savedVars%jloz1(j,ncsp) = jccount + 1

        do i = 1, jzilch(j)
          jccount         = jccount + 1
          ia              = jzero(jccount)
          savedVars%jzeroa(jccount) = jarraypt(ia,j)
        end do

        savedVars%jhiz1(j,ncsp) = jccount

!     =============
      end do JLOOP1
!     =============


!     ------------------------------------------------------------
!     Create more back-substitution arrays to increase efficiency.
!     ------------------------------------------------------------

!     -- First loop of back-substitution. --

      savedVars%kztlo(ncsp) = kztot + 1

!     ==============================
      ILOOP2: do i = 2, savedVars%ischang(ncs)
!     ==============================

        kzil = kzilch(i)

        if (kzil > 0) then

          kztot   = kztot + 1
          nremain = kzil
          nfive   = (nremain + 0.0001d0) / 5
          nremain =  nremain - (nfive  * 5)
          nfour   = (nremain + 0.0001d0) / 4
          nremain =  nremain - (nfour  * 4)
          nthree  = (nremain + 0.0001d0) / 3
          nremain =  nremain - (nthree * 3)
          ntwo    = (nremain + 0.0001d0) / 2
          nremain =  nremain - (ntwo   * 2)
          none    = (nremain + 0.0001d0)
          nremain =  nremain - none

          savedVars%ikztot(kztot) = i
          savedVars%kbl5  (kztot) = kbsub       + 1
          savedVars%kbh5  (kztot) = kbsub       + nfive
          savedVars%kbl4  (kztot) = savedVars%kbh5(kztot) + 1
          savedVars%kbh4  (kztot) = savedVars%kbh5(kztot) + nfour
          savedVars%kbl3  (kztot) = savedVars%kbh4(kztot) + 1
          savedVars%kbh3  (kztot) = savedVars%kbh4(kztot) + nthree
          savedVars%kbl2  (kztot) = savedVars%kbh3(kztot) + 1
          savedVars%kbh2  (kztot) = savedVars%kbh3(kztot) + ntwo
          savedVars%kbl1  (kztot) = savedVars%kbh2(kztot) + 1
          savedVars%kbh1  (kztot) = savedVars%kbh2(kztot) + none
          kbsub         = savedVars%kbh1(kztot)

          do kc = savedVars%kbl5(kztot), savedVars%kbh5(kztot)
            savedVars%kzeroa(kc) = kzero(kccount+1)
            savedVars%kzerob(kc) = kzero(kccount+2)
            savedVars%kzeroc(kc) = kzero(kccount+3)
            savedVars%kzerod(kc) = kzero(kccount+4)
            savedVars%kzeroe(kc) = kzero(kccount+5)
            kccount    = kccount + 5
          end do

          do kc = savedVars%kbl4(kztot), savedVars%kbh4(kztot)
            savedVars%kzeroa(kc) = kzero(kccount+1)
            savedVars%kzerob(kc) = kzero(kccount+2)
            savedVars%kzeroc(kc) = kzero(kccount+3)
            savedVars%kzerod(kc) = kzero(kccount+4)
            kccount    = kccount + 4
          end do

          do kc = savedVars%kbl3(kztot), savedVars%kbh3(kztot)
            savedVars%kzeroa(kc) = kzero(kccount+1)
            savedVars%kzerob(kc) = kzero(kccount+2)
            savedVars%kzeroc(kc) = kzero(kccount+3)
            kccount    = kccount + 3
          end do

          do kc = savedVars%kbl2(kztot), savedVars%kbh2(kztot)
            savedVars%kzeroa(kc) = kzero(kccount+1)
            savedVars%kzerob(kc) = kzero(kccount+2)
            kccount    = kccount + 2
          end do

          do kc = savedVars%kbl1(kztot), savedVars%kbh1(kztot)
            savedVars%kzeroa(kc) = kzero(kccount+1)
            kccount    = kccount + 1
          end do

        end if

!     =============
      end do ILOOP2
!     =============


      savedVars%kzthi(ncsp) = kztot


!     -- Second loop of back-substitution. --

      do i = savedVars%ischang(ncs), 1, -1

        mzil = mzilch(i)

        if (mzil > 0) then

          mztot   = mztot + 1
          nremain = mzil
          nfive   = (nremain + 0.0001d0) / 5
          nremain =  nremain - (nfive  * 5)
          nfour   = (nremain + 0.0001d0) / 4
          nremain =  nremain - (nfour  * 4)
          nthree  = (nremain + 0.0001d0) / 3
          nremain =  nremain - (nthree * 3)
          ntwo    = (nremain + 0.0001d0) / 2
          nremain =  nremain - (ntwo   * 2)
          none    = (nremain + 0.0001d0)
          nremain =  nremain - none

          savedVars%imztot(i,ncsp) = mztot
          savedVars%mbl5  (mztot)  = mbsub       + 1
          savedVars%mbh5  (mztot)  = mbsub       + nfive
          savedVars%mbl4  (mztot)  = savedVars%mbh5(mztot) + 1
          savedVars%mbh4  (mztot)  = savedVars%mbh5(mztot) + nfour
          savedVars%mbl3  (mztot)  = savedVars%mbh4(mztot) + 1
          savedVars%mbh3  (mztot)  = savedVars%mbh4(mztot) + nthree
          savedVars%mbl2  (mztot)  = savedVars%mbh3(mztot) + 1
          savedVars%mbh2  (mztot)  = savedVars%mbh3(mztot) + ntwo
          savedVars%mbl1  (mztot)  = savedVars%mbh2(mztot) + 1
          savedVars%mbh1  (mztot)  = savedVars%mbh2(mztot) + none
          mbsub          = savedVars%mbh1(mztot)

          do mc = savedVars%mbl5(mztot), savedVars%mbh5(mztot)
            savedVars%mzeroa(mc) = mzero(mccount+1)
            savedVars%mzerob(mc) = mzero(mccount+2)
            savedVars%mzeroc(mc) = mzero(mccount+3)
            savedVars%mzerod(mc) = mzero(mccount+4)
            savedVars%mzeroe(mc) = mzero(mccount+5)
            mccount    = mccount + 5
          end do

          do mc = savedVars%mbl4(mztot), savedVars%mbh4(mztot)
            savedVars%mzeroa(mc) = mzero(mccount+1)
            savedVars%mzerob(mc) = mzero(mccount+2)
            savedVars%mzeroc(mc) = mzero(mccount+3)
            savedVars%mzerod(mc) = mzero(mccount+4)
            mccount    = mccount + 4
          end do

          do mc = savedVars%mbl3(mztot), savedVars%mbh3(mztot)
            savedVars%mzeroa(mc) = mzero(mccount+1)
            savedVars%mzerob(mc) = mzero(mccount+2)
            savedVars%mzeroc(mc) = mzero(mccount+3)
            mccount    = mccount + 3
          end do

          do mc = savedVars%mbl2(mztot), savedVars%mbh2(mztot)
            savedVars%mzeroa(mc) = mzero(mccount+1)
            savedVars%mzerob(mc) = mzero(mccount+2)
            mccount    = mccount + 2
          end do

          do mc = savedVars%mbl1(mztot), savedVars%mbh1(mztot)
            savedVars%mzeroa(mc) = mzero(mccount+1)
            mccount    = mccount + 1
          end do

        end if

      end do


!     ---------------------------------------------
!     Check dimensions and print out array savings.
!     ---------------------------------------------

      if ((icnt    > MXCOUNT2) .or. (jcnt    > MXCOUNT3) .or.  &
     &    (kcnt    > MXCOUNT3) .or. (mcnt    > MXCOUNT3) .or.  &
     &    (iccount > MXCOUNT2) .or. (jccount > MXCOUNT3) .or.  &
     &    (kccount > MXCOUNT3) .or. (mccount > MXCOUNT3) .or.  &
     &    (ijtot   > MXCOUNT3) .or. (idecomp > MXCOUNT3) .or.  &
     &    (kztot   > MXCOUNT4) .or. (kbsub   > MXCOUNT4) .or.  &
     &    (mztot   > MXCOUNT4) .or. (mbsub   > MXCOUNT4) .or.  &
     &    (iarray2 > MXARRAY)) then

        Write (6,800)  &
     &    MXCOUNT2, icnt,    MXCOUNT3, jcnt,  &
     &    MXCOUNT3, kcnt,    MXCOUNT3, mcnt,  &
     &    MXCOUNT2, iccount, MXCOUNT3, jccount,  &
     &    MXCOUNT3, kccount, MXCOUNT3, mccount,  &
     &    MXCOUNT3, ijtot,   MXCOUNT3, idecomp,  &
     &    MXCOUNT4, kztot,   MXCOUNT4, kbsub,  &
     &    MXCOUNT4, mztot,   MXCOUNT4, mbsub,  &
     &    MXARRAY,  iarray2

        err_msg = ' Problem in Ksparse '
        call GmiPrintError (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)

      end if

 800  format  &
     &  ('One of the dimensions below is too small:', /,  &
     &   'Dimension: MXCOUNT2 = ', i4, ' Variable: icnt    = ', i4, /,  &
     &   'Dimension: MXCOUNT3 = ', i4, ' Variable: jcnt    = ', i4, /,  &
     &   'Dimension: MXCOUNT3 = ', i4, ' Variable: kcnt    = ', i4, /,  &
     &   'Dimension: MXCOUNT3 = ', i4, ' Variable: mcnt    = ', i4, /,  &
     &   'Dimension: MXCOUNT2 = ', i4, ' Variable: iccount = ', i4, /,  &
     &   'Dimension: MXCOUNT3 = ', i4, ' Variable: jccount = ', i4, /,  &
     &   'Dimension: MXCOUNT3 = ', i4, ' Variable: kccount = ', i4, /,  &
     &   'Dimension: MXCOUNT3 = ', i4, ' Variable: mccount = ', i4, /,  &
     &   'Dimension: MXCOUNT3 = ', i4, ' Variable: ijtot   = ', i4, /,  &
     &   'Dimension: MXCOUNT3 = ', i4, ' Variable: idecomp = ', i4, /,  &
     &   'Dimension: MXCOUNT4 = ', i4, ' Variable: kztot   = ', i4, /,  &
     &   'Dimension: MXCOUNT4 = ', i4, ' Variable: kbsub   = ', i4, /,  &
     &   'Dimension: MXCOUNT4 = ', i4, ' Variable: mztot   = ', i4, /,  &
     &   'Dimension: MXCOUNT4 = ', i4, ' Variable: mbsub   = ', i4, /,  &
     &   'Dimension: MXARRAY  = ', i4, ' Variable: iarray2 = ', i4)


      if (pr_smv2) then
        Write (lunsmv,850)  &
     &    ncsp, kount0a, kount0, kntarray, iarray2, icnta, icntb,  &
     &    jcnta, jcntb, kcnta, kcntb, mcnta, mcntb
      end if

 850  format  &
     &  (/, 'Param    Poss Matrix Points -- Nonzeros -- ncsp=', i4, /,  &
     &      'Initmat  ', 4x, i8, 9x, i8, /  &
     &      'Finmat   ', 4x, i8, 9x, i8, /  &
     &      'Decomp1  ', 4x, i8, 9x, i8, /  &
     &      'Decomp2  ', 4x, i8, 9x, i8, /  &
     &      'Backsb1  ', 4x, i8, 9x, i8, /  &
     &      'Backsb2  ', 4x, i8, 9x, i8, /)


      return

      end subroutine KsparseGMI

