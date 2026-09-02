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
!   physproc.F
!
! ROUTINES
!   Physproc
!   Deter_Block_Size
!   Solve_Block
!   Calcrate
!   Reorder_Grid_Cells
!
!=============================================================================


!-----------------------------------------------------------------------------
!
! ROUTINE
!   Physproc
!
! DESCRIPTION
!   This routine solves gas-phase chemical equations.  It divides the grid-
!   domain into grid-blocks, and the code vectorizes around the number of
!   grid-cells in each block;
!     grid-block  => a group of grid-cells in the grid-domain,
!     grid-cell   => a single grid-box,
!     grid-domain => all grid-cells in the model.
!   It sets compact arrays and calls Smvgear (Sparse Matrix Vectorized Gear
!   code).
!
! ARGUMENTS
!   do_qqjk_inchem   : if pr_qqjk is on, should qqj's & qqk's be determined
!                      inside the chemistry solver, or outside?
!   do_semiss_inchem : do surface emissions inside the chemistry solver, or
!                      outside?
!   pr_qqjk   : should the periodic qqjk output file be written?
!   pr_smv2   : should the SmvgearII     output file be written
!               (non-parallel mode only)?
!   ilat      : # of latitudes
!   ilong     : # of longitudes
!   ivert     : # of vertical layers
!   ifreord   : if 1, then reorder grid-cells by stiffness
!   imgas     : array index for air density
!   initrogen : identifies spc # of nitrogen gas
!   ioxygen   : identifies spc # of oxygen   gas
!   itloop    : # of zones (ilong * ilat * ivert)
!   kuloop    : intended # of grid-cells in a grid-block
!   lunsmv    : logical unit number to write to when pr_smv2 is true
!   ncs       : identifies gas chemistry type (1..NCSGAS)
!   fracdec   : fraction time step is decreased in Smvgear if convergence
!               test fails
!   hmaxnit   : max time step for night all chem (s)
!   pr_nc_period     : NetCDF output period
!   tdt       : model time step (s)
!   do_cell_chem     : do chemistry for a particular cell?
!   jphotrat  : tbd
!   nrates    : # of kinetic rxns (non-photo)
!   ntloopncs : tbd
!   ntspec    : # of active + inactive gases
!   inewold   : original spc # of each new jnew spc
!   npphotrat : tbd
!   arate     : thermal    rate constants (units vary)
!   prate     : photolysis rate constants (s^-1)
!   yemis     : surface emissions (molecules/cm^3/sec)
!   jreorder  : gives original grid-cell from re-ordered grid-cell
!   lreorder  : gives original grid-cell from re-ordered cell, except when
!               cell is a virtual boundary cell, then it gives original
!               edge cell from re-ordered virtual boundary cell
!   csuma     : tbd
!   csumc     : tbd
!   errmx2    : tbd
!   cx        : spc conc (molec/cm^3)
!
! HISTORY
!   - July 1, 2004 - Jules Kouatchou
!       o Use the pre-processing option MSG_OPTION to remove MPI
!         statements when MSG_OPTION is set to MSG_NONE.
!-----------------------------------------------------------------------------

      subroutine PhysprocGMI  &
     &  (savedVars, do_qqjk_inchem, do_semiss_inchem, pr_qqjk, pr_smv2, ilat,  &
     &   ilong, ivert, ifreord, imgas, initrogen, ioxygen, itloop,  &
     &   kuloop, lunsmv, ncs, fracdec, hmaxnit, pr_nc_period, tdt,  &
     &   do_cell_chem, jphotrat, nrates, ntloopncs, ntspec, inewold,  &
     &   npphotrat, arate, prate, yemis, jreorder, lreorder, csuma,  &
     &   csumc, errmx2, cx, yda, qqkda, qqjda, qkgmi, qjgmi, &
     &   i1, i2, ju1, j2, k1, k2, num_qjo, num_qks, num_qjs, num_active)

      use GmiSavedVariables_mod, only : t_ChemistrySaved
      use GmiPrintError_mod, only : GmiPrintError

      implicit none

#     include "smv2chem_par.h"

!     ----------------------
!     Argument declarations.
!     ----------------------

      integer, intent(in) :: i1, i2, ju1, j2, k1, k2
      integer, intent(in) :: num_qjo, num_qks, num_qjs, num_active
      real*8 , intent(in   ) :: qjgmi(i1:i2, ju1:j2, k1:k2, num_qjo)
      real*8 , intent(in   ) :: qkgmi(i1:i2, ju1:j2, k1:k2, num_qks)
      real*8 , intent(inout) :: qqjda(i1:i2, ju1:j2, k1:k2, num_qjs)
      real*8 , intent(inout) :: qqkda(i1:i2, ju1:j2, k1:k2, num_qks)
      real*8 , intent(inout) :: yda  (i1:i2, ju1:j2, k1:k2, num_active)
      logical, intent(in)  :: do_qqjk_inchem
      logical, intent(in)  :: do_semiss_inchem
      logical, intent(in)  :: pr_qqjk
      logical, intent(in)  :: pr_smv2
      integer, intent(in)  :: ilat, ilong, ivert
      integer, intent(in)  :: ifreord
      integer, intent(in)  :: imgas
      integer, intent(in)  :: initrogen
      integer, intent(in)  :: ioxygen
      integer, intent(in)  :: itloop
      integer, intent(in)  :: kuloop
      integer, intent(in)  :: lunsmv
      integer, intent(in)  :: ncs
      real*8,  intent(in)  :: fracdec
      real*8,  intent(in)  :: hmaxnit
      real*8,  intent(in)  :: pr_nc_period
      real*8,  intent(in)  :: tdt
      logical, intent(in)  :: do_cell_chem(itloop)
      integer, intent(in)  :: jphotrat (ICS)
      integer, intent(in)  :: nrates   (ICS)
      integer, intent(in)  :: ntloopncs(ICS)
      integer, intent(in)  :: ntspec   (ICS)
      integer, intent(in)  :: inewold  (MXGSAER, ICS)
      integer, intent(in)  :: npphotrat(IPHOT,   ICS)
      real*8,  intent(in)  :: arate    (itloop,  ITHERM)
      real*8,  intent(in)  :: prate    (itloop,  IPHOT)
      real*8,  intent(in)  :: yemis    (ilat*ilong, IGAS)

      type(t_ChemistrySaved), intent(inOut) :: savedVars

      integer, intent(inout) :: jreorder(itloop)
      integer, intent(inout) :: lreorder(itloop)
      real*8,  intent(inout) :: csuma   (itloop)
      real*8,  intent(inout) :: csumc   (itloop)
      real*8,  intent(inout) :: errmx2  (itloop)
      real*8,  intent(inout) :: cx      (itloop, IGAS)


!     ----------------------
!     Variable declarations.
!     ----------------------

      integer :: idaynt
      integer :: ifsun      ! 1 => for daytime, 2 => for nighttime
      integer :: ireord     ! 1 => reorder grid-cells and blocks for chemistry
                            ! 2 => solve chemistry
      integer :: iday
      integer :: jloop, jloopn
      integer :: loreord    ! 1 => if reordering, 2 => if no reordering
      integer :: nblockuse  ! (ireord == loreord) => # of original  blocks
                            ! (ireord /= loreord) => # of reordered blocks
      integer :: nreblock
      integer :: ntloopuse

      integer :: jlowvar(MXBLOCK)
      integer :: ktlpvar(MXBLOCK)
      integer :: myId, myIerr

!     ----------------
!     Begin execution.
!     ----------------


!     --------------------------------------------------------
!     Can the domain be separated into day and night sections?
!     --------------------------------------------------------

      if (prate(1,1) >= 1.0d-80) then
        ifsun = 1
      else
        ifsun = 2
      end if

      idaynt = 1

      if (ifsun == 1) then

        IFSUN1: do jloop = 2, ntloopncs(ncs)
          if (prate(jloop,1) < 1.0d-80) then
            idaynt = 2
            exit IFSUN1
          end if
        end do IFSUN1

      else if (ifsun == 2) then

        IFSUN2: do jloop = 2, ntloopncs(ncs)
          if (prate(jloop,1) > 1.0d-80) then
            idaynt = 2
            exit IFSUN2
          end if
        end do IFSUN2

      end if


!     --------------------------------------------------
!     Reorder cells and blocks then solve chemical odes.
!     --------------------------------------------------

      if ((ifreord == 1) .and. (itloop > 1)) then
        loreord = 1
      else
        loreord = 2
      end if

      do iday = 1, idaynt
        do ireord = loreord, 2

          if (ireord == loreord) then
!           =====================
            call Deter_Block_Size  &
!           =====================
     &        (savedVars, iday, idaynt, itloop, kuloop, ncs, do_cell_chem,  &
     &         ntloopncs, prate, ifsun, nblockuse, ntloopuse, jlowvar,  &
     &         ktlpvar, jreorder)
          else
            nblockuse = nreblock
          end if

          do jloopn = 1, ntloopuse
            lreorder(jloopn) = jreorder(jloopn)
          end do

!         ================
          call Solve_Block  &
!         ================
     &      (savedVars, do_qqjk_inchem, do_semiss_inchem, pr_qqjk, pr_smv2, ifsun,  &
     &       ilat, ilong, ivert, imgas, initrogen, ioxygen, ireord,  &
     &       itloop, lunsmv, nblockuse, ncs, fracdec, hmaxnit,  &
     &       pr_nc_period, tdt, do_cell_chem, jphotrat, nrates, ntspec,  &
     &       jlowvar, ktlpvar, inewold, npphotrat, arate, prate, yemis,  &
     &       jreorder, lreorder, errmx2, cx, yda, qqkda, qqjda, qkgmi, &
             qjgmi, i1, i2, ju1, j2, k1, k2, num_qjo, num_qks, num_qjs,&
             num_active)


          if (ireord == 1) then
!           =======================
            call Reorder_Grid_Cells  &
!           =======================
     &        (itloop, kuloop, ntloopuse, errmx2, jreorder, csuma,  &
     &         nreblock, lreorder, jlowvar, ktlpvar, csumc)
          end if

        end do
      end do

      return

      end subroutine PhysprocGMI

!-----------------------------------------------------------------------------
!
! ROUTINE
!   Deter_Block_Size
!
! DESCRIPTION
!   This routine determines block sizes for chemistry or stiffness.
!
! ARGUMENTS
!   iday      : tbd
!   idaynt    : tbd
!   itloop    : # of zones (ilong * ilat * ivert)
!   kuloop    : intended # of grid-cells in a grid-block
!   ncs       : identifies gas chemistry type (1..NCSGAS)
!   do_cell_chem : do chemistry for a particular cell?
!   ntloopncs : tbd
!   prate     : photolysis rate constants (s^-1)
!   ifsun     : 1 => for daytime, 2 => for nighttime
!   nblockuse : # of original blocks
!   ntloopuse : tbd
!   jlowvar   : tbd
!   ktlpvar   : tbd
!   jreorder  : gives original grid-cell from re-ordered grid-cell
!
!-----------------------------------------------------------------------------

      subroutine Deter_Block_Size  &
     &  (savedVars, iday, idaynt, itloop, kuloop, ncs, do_cell_chem, ntloopncs,  &
     &   prate, ifsun, nblockuse, ntloopuse, jlowvar, ktlpvar, jreorder)


      use GmiSavedVariables_mod, only : t_ChemistrySaved
      use GmiPrintError_mod, only : GmiPrintError

      implicit none

#     include "smv2chem_par.h"

!     ----------------------
!     Argument declarations.
!     ----------------------

      integer, intent(in)  :: iday
      integer, intent(in)  :: idaynt
      integer, intent(in)  :: itloop
      integer, intent(in)  :: kuloop
      integer, intent(in)  :: ncs
      logical, intent(in)  :: do_cell_chem(itloop)
      integer, intent(in)  :: ntloopncs(ICS)
      real*8,  intent(in)  :: prate(itloop, IPHOT)

      integer, intent(out) :: ifsun
      integer, intent(out) :: nblockuse
      integer, intent(out) :: ntloopuse
      integer, intent(out) :: jlowvar (MXBLOCK)
      integer, intent(out) :: ktlpvar (MXBLOCK)
      integer, intent(out) :: jreorder(itloop)
      type(t_ChemistrySaved), intent(inOut) :: savedVars


!     ----------------------
!     Variable declarations.
!     ----------------------

      integer :: iavblok, iavgsize
      integer :: iremain
      integer :: jloop
      integer :: jlooplo  ! low ntloop grid-cell - 1 in a grid-block
      integer :: kblk
      integer :: nblock1

!     ----------------
!     Begin execution.
!     ----------------

      if (idaynt == 1) then

        ntloopuse = 0

        do jloop = 1, ntloopncs(ncs)
          if (do_cell_chem(jloop)) then
            ntloopuse           = ntloopuse + 1
            jreorder(ntloopuse) = jloop
          end if
        end do

      else

        ntloopuse = 0

        if (iday == 1) then

          ifsun = 1

          do jloop = 1, ntloopncs(ncs)
            if ((prate(jloop,1) > 1.0d-80) .and.  &
     &          do_cell_chem(jloop)) then
              ntloopuse           = ntloopuse + 1
              jreorder(ntloopuse) = jloop
            end if
          end do

        else

          ifsun = 2

          do jloop = 1, ntloopncs(ncs)
            if ((prate(jloop,1) < 1.0d-80) .and.  &
     &          do_cell_chem(jloop)) then
              ntloopuse           = ntloopuse + 1
              jreorder(ntloopuse) = jloop
            end if
          end do

        end if

      end if


      nblockuse = 1 + ntloopuse / (kuloop    + 0.0001d0)
      iavblok   = 1 + ntloopuse / (nblockuse + 0.0001d0)
      iavgsize  = Min (iavblok, kuloop)
      jlooplo   = 0
      nblock1   = nblockuse - 1


      if (nblockuse > savedVars%nblockuse_max) then
        savedVars%nblockuse_max = nblockuse
!c      Write (6,*) 'nblockuse_max:  ', savedVars%nblockuse_max
      end if

      if (savedVars%nblockuse_max > MXBLOCK) then
        Write (6,*) 'nblockuse > MXBLOCK', nblockuse
        call GmiPrintError ('Deter_Block_Size', .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
      end if


      do kblk = 1, nblock1
        jlowvar(kblk) = jlooplo
        ktlpvar(kblk) = iavgsize
        jlooplo       = jlooplo + iavgsize
      end do

      iremain = Max (ntloopuse - (nblock1 * iavgsize), 0)

      jlowvar(nblockuse) = jlooplo
      ktlpvar(nblockuse) = iremain

      return

      end subroutine Deter_Block_Size

!-----------------------------------------------------------------------------
!
! ROUTINE
!   Solve_Block
!
! DESCRIPTION
!   This routine solves the chemical odes for each block.
!
! ARGUMENTS
!   do_qqjk_inchem   : if pr_qqjk is on, should qqj's & qqk's be determined
!                      inside the chemistry solver, or outside?
!   do_semiss_inchem : do surface emissions inside the chemistry solver, or
!                      outside?
!   pr_qqjk   : should the periodic qqjk output file be written?
!   pr_smv2   : should the SmvgearII     output file be written
!               (non-parallel mode only)?
!   ifsun     : 1 => for daytime, 2 => for nighttime
!   ilat      : # of latitudes
!   ilong     : # of longitudes
!   ivert     : # of vertical layers
!   imgas     : array index for air density
!   initrogen : identifies spc # of nitrogen gas
!   ioxygen   : identifies spc # of oxygen   gas
!   ireord    : 1 => reorder grid-cells and blocks for chemistry
!               2 => solve chemistry
!   itloop    : # of zones (ilong * ilat * ivert)
!   lunsmv    : logical unit number to write to when pr_smv2 is true
!   nblockuse : (ireord == loreord) => # of original  blocks
!               (ireord /= loreord) => # of reordered blocks
!   ncs       : identifies gas chemistry type (1..NCSGAS)
!   fracdec   : fraction time step is decreased in Smvgear if convergence
!               test fails
!   hmaxnit   : max time step for night all chem (s)
!   pr_nc_period     : NetCDF output period
!   tdt       : model time step (s)
!   do_cell_chem     : do chemistry for a particular cell?
!   jphotrat  : tbd
!   nrates    : # of kinetic rxns (non-photo)
!   ntspec    : # of active + inactive gases
!   jlowvar   : tbd
!   ktlpvar   : tbd
!   inewold   : original spc # of each new jnew spc
!   npphotrat : tbd
!   arate     : thermal    rate constants (units vary)
!   prate     : photolysis rate constants (s^-1)
!   yemis     : surface emissions (molec/cm^3/sec)
!   jreorder  : gives original grid-cell from re-ordered grid-cell
!   lreorder  : gives original grid-cell from re-ordered cell, except when
!               cell is a virtual boundary cell, then it gives original
!               edge cell from re-ordered virtual boundary cell
!   errmx2    : tbd
!   cx        : spc conc (molec/cm^3)
!
!-----------------------------------------------------------------------------

      subroutine Solve_Block  &
     &  (savedVars, do_qqjk_inchem, do_semiss_inchem, pr_qqjk, pr_smv2, ifsun,  &
     &   ilat, ilong, ivert, imgas, initrogen, ioxygen, ireord,  &
     &   itloop, lunsmv, nblockuse, ncs, fracdec, hmaxnit,  &
     &   pr_nc_period, tdt, do_cell_chem, jphotrat, nrates, ntspec,  &
     &   jlowvar, ktlpvar, inewold, npphotrat, arate, prate, yemis,  &
     &   jreorder, lreorder, errmx2, cx, yda, qqkda, qqjda, qkgmi, qjgmi, &
     &   i1, i2, ju1, j2, k1, k2, num_qjo, num_qks, num_qjs, num_active)


      use GmiSavedVariables_mod, only : t_ChemistrySaved
      use GmiPrintError_mod, only : GmiPrintError

      implicit none

#     include "smv2chem_par.h"
!#     include "gem_sys_options.h"

!     ==========================
!#if (MSG_OPTION == MSG_MPI)
!       include "mpif.h"
!#endif
!#     include "gmi_subdomains.h"
!     ==========================


!     ----------------------
!     Argument declarations.
!     ----------------------

      integer, intent(in) :: i1, i2, ju1, j2, k1, k2
      integer, intent(in) :: num_qjo, num_qks, num_qjs, num_active
      real*8 , intent(in   ) :: qjgmi(i1:i2, ju1:j2, k1:k2, num_qjo)
      real*8 , intent(in   ) :: qkgmi(i1:i2, ju1:j2, k1:k2, num_qks)
      real*8 , intent(inout) :: qqjda(i1:i2, ju1:j2, k1:k2, num_qjs)
      real*8 , intent(inout) :: qqkda(i1:i2, ju1:j2, k1:k2, num_qks)
      real*8 , intent(inout) :: yda  (i1:i2, ju1:j2, k1:k2, num_active)
      logical, intent(in)  :: do_qqjk_inchem
      logical, intent(in)  :: do_semiss_inchem
      logical, intent(in)  :: pr_qqjk
      logical, intent(in)  :: pr_smv2
      integer, intent(in)  :: ifsun
      integer, intent(in)  :: ilat, ilong, ivert
      integer, intent(in)  :: imgas
      integer, intent(in)  :: initrogen
      integer, intent(in)  :: ioxygen
      integer, intent(in)  :: ireord
      integer, intent(in)  :: itloop
      integer, intent(in)  :: lunsmv
      integer, intent(in)  :: nblockuse
      integer, intent(in)  :: ncs
      real*8,  intent(in)  :: fracdec
      real*8,  intent(in)  :: hmaxnit
      real*8,  intent(in)  :: pr_nc_period
      real*8,  intent(in)  :: tdt
      logical, intent(in)  :: do_cell_chem(itloop)
      integer, intent(in)  :: jphotrat (ICS)
      integer, intent(in)  :: nrates   (ICS)
      integer, intent(in)  :: ntspec   (ICS)
      integer, intent(in)  :: jlowvar  (MXBLOCK)
      integer, intent(in)  :: ktlpvar  (MXBLOCK)
      integer, intent(in)  :: inewold  (MXGSAER, ICS)
      integer, intent(in)  :: npphotrat(IPHOT,   ICS)
      real*8,  intent(in)  :: arate    (itloop,  ITHERM)
      real*8,  intent(in)  :: prate    (itloop,  IPHOT)
      real*8,  intent(in)  :: yemis    (ilat*ilong, IGAS)

      integer, intent(inout) :: jreorder(itloop)
      integer, intent(inout) :: lreorder(itloop)
      real*8,  intent(inout) :: errmx2  (itloop)
      real*8,  intent(inout) :: cx      (itloop, IGAS)

      type(t_ChemistrySaved), intent(inOut) :: savedVars

!     -----------------------
!     Parameter declarations.
!     -----------------------

      logical, parameter :: DOWRT_SBDIAG = .false.


!     ----------------------
!     Variable declarations.
!     ----------------------

      integer :: ierr
      integer :: j
      integer :: jgas
      integer :: jloop
      integer :: jlooplo  ! low ntloop grid-cell - 1 in a grid-block
      integer :: jnew
      integer :: kblk
      integer :: kloop
      integer :: ktloop   ! # of grid-cells in a grid-block
      integer :: nallr    ! # of active rxns
      integer :: nfdh1    ! nfdh2 + # of rxns with one   active reactant
      integer :: nfdh2    ! nfdh3 + # of rxns with two   active reactants
      integer :: nfdh3    !         # of rxns with three active reactants
      integer :: nfdl1    ! nfdh2 + 1
      integer :: nfdl2    ! nfdh3 + 1
      integer :: nfdrep   ! nfdh3 + # of rxns with two active reactants that
                          ! are not followed by a rxn with the same reactants
      integer :: nfdrep1  ! nfdrep + 1
      integer :: np
      integer :: proc_num

!     ------------------------------------------------------------------
!     irma,b,c : spc # of each reactant; locates reordered active spc #s
!     ------------------------------------------------------------------

      integer :: irma(NMTRATE)
      integer :: irmb(NMTRATE)
      integer :: irmc(NMTRATE)

      real*8  :: denair(MXBLOCK)  ! density of air (molec/cm^3)

!     -----------------------------------------------------------------------
!     cc2    : array holding values of decomposed matrix
!     pratk1 : tbd
!     cblk   : gas-phase concs (molec/cm^3)
!     cnew   : stores conc (y (estimated)) (molec/cm^3)
!     corig  : original gas-phase concs used to restart Smvgear if a failure
!              occurs (molec/cm^3)
!     gloss  : value of first derivatives on output from Subfun; right-side
!              of eqn on input to Backsub; error term (solution from Backsub)
!              on output from Backsub
!     smvdm  : amount added to each spc at each grid-cell, for mass balance
!              accounting (# cm^-3 for gas chemistry (?))
!     vdiag  : 1 / current diagonal term of the decomposed matrix
!     rrate  : rate constants
!              (s^-1, cm^3/molec*s, cm^6/molec*s^2, or cm^9/molec*s^3)
!     trate  : rxn rate (moles l^-1-h2o s^-1 or # cm^-3 s^-1 (?))
!     urate  : term of Jacobian (J) = partial derivative
!     -----------------------------------------------------------------------

      real*8  :: cc2   (KBLOOP, 0:MXARRAY)  = 0.0d0

      real*8  :: pratk1(KBLOOP, IPHOT)      = 0.0d0

      real*8  :: cblk  (KBLOOP, MXGSAER)    = 0.0d0
      real*8  :: cnew  (KBLOOP, MXGSAER)    = 0.0d0
      real*8  :: corig (KBLOOP, MXGSAER)    = 0.0d0
      real*8  :: gloss (KBLOOP, MXGSAER)    = 0.0d0
      real*8  :: smvdm (KBLOOP, MXGSAER)    = 0.0d0
      real*8  :: vdiag (KBLOOP, MXGSAER)    = 0.0d0

      real*8  :: rrate (KBLOOP, NMTRATE)    = 0.0d0
      real*8  :: trate (KBLOOP, NMTRATE*2)  = 0.0d0
      real*8  :: urate (KBLOOP, NMTRATE, 3) = 0.0d0

!     ===============================
!$    integer :: Omp_Get_Max_Threads
!$    integer :: Omp_Get_Num_Threads
!c    integer :: Omp_Get_Thread_Num

!$    integer :: max_thrds, num_thrds
!c    integer :: thrd_num
!     ===============================

!     =========================
!c    real*8  :: wclk1, wclk2
!c    real*8  :: wtime
!c    real*8  :: ctime(MXBLOCK)
!c    real*8  :: sumtime(0:15)
!     =========================
!JK
      integer :: myId, myIerr


!     ----------------
!     Begin execution.
!     ----------------

      if (DOWRT_SBDIAG) then
!       ====================================================================
!#if (MSG_OPTION == MSG_MPI)
!        call Mpi_Comm_Rank (commu_slaves, proc_num, ierr)
!        proc_num = proc_num + 1  ! change proc_num range from 0->N to 1->N+1
!#endif

!$omp   parallel
!$      num_thrds = Omp_Get_Num_Threads ( )
!$omp   end parallel

!$      max_thrds = Omp_Get_Max_Threads ( )
!$      Write (6,900) proc_num, max_thrds, num_thrds, nblockuse
!$ 900   format ('Proc #, Max Thrds, # Thrds, nblockuse:  ', i6, i4, i4, i8)
!       ====================================================================
      end if

!c    call f_hpmstart (1, "chem block loop")

!c    call Mpi_Comm_Rank (commu_slaves, proc_num, ierr)
!c    ctime  (:) = 0.0d0
!c    sumtime(:) = 0.0d0

!$omp   parallel do                                  &
!$omp   default(shared)                              &
!$omp   schedule(runtime)                            &
!$omp   private(j, jgas, jnew, kblk, np)             &
!$omp   private(jloop, jlooplo, kloop, ktloop)       &
!$omp   private(nallr, nfdrep, nfdrep1)              &
!$omp   private(nfdh1, nfdh2, nfdh3, nfdl1, nfdl2)   &
!$omp   private(cblk, cc2, cnew, corig)              &
!$omp   private(denair, gloss)                       &
!$omp   private(irma, irmb, irmc)                    &
!$omp   private(pratk1, smvdm, vdiag)                &
!$omp   private(rrate, trate, urate)

!     ======================
      do kblk = 1, nblockuse
!     ======================

!c      thrd_num  = Omp_Get_Thread_Num  ( )
!c      num_thrds = Omp_Get_Num_Threads ( )
!c      wclk1     = Mpi_Wtime (ierr)

        jlooplo = jlowvar(kblk)
        ktloop  = ktlpvar(kblk)

!       ================
        if (ktloop /= 0) then
!       ================

!         -------------------------------------
!         Set (and rearrange) photofrequencies.
!         -------------------------------------

          do j = 1, jphotrat(ncs)
            np = npphotrat(j,ncs)
            do kloop = 1, ktloop
              jloop           = lreorder(jlooplo+kloop)
              pratk1(kloop,j) = prate(jloop,np)
            end do
          end do

!         -------------------------------------------------------
!         Place large domain gas array (molec/cm^3) into smaller
!         block array.
!
!         The urate=arate loop is where the 2D diurnal averaging
!         factors for the thermal reactions are passed in.  urate
!         later is used for another purpose.
!         -------------------------------------------------------

          do j = 1, nrates(ncs)
            do kloop = 1, ktloop
              jloop            = lreorder(jlooplo+kloop)
              urate(kloop,j,1) = arate(jloop,j)
            end do
          end do

!         --------------------------------------------------------------
!         Calculate rates and solve chemistry.
!
!         ireord = 1 : call Calcrate to find stiffness of each grid-cell
!         ireord = 2 : set chemistry rates and solve equations
!         --------------------------------------------------------------

!         =============
          call CalcrateGMI  &
!         =============
     &      (ifsun, imgas, initrogen, ioxygen, itloop, jlooplo,  &
     &       ktloop, ncs, jreorder, nrates, ntspec, cx, urate,  &
     &       denair, pratk1, cblk, rrate, trate, nallr, nfdh1, nfdh2,  &
     &       nfdh3, nfdl1, nfdl2, nfdrep, nfdrep1, irma, irmb, irmc,  &
     &       corig, smvdm, savedVars)

!         --------------------
!         Solve chemical odes.
!         --------------------
!         ============
          call SmvgearGMI  &
!         ============
     &      (savedVars, do_qqjk_inchem, do_semiss_inchem, pr_qqjk, pr_smv2,  &
     &       ifsun, ilat, ilong, ivert, ireord, itloop, jlooplo,  &
     &       ktloop, lunsmv, nallr, ncs, nfdh2, nfdh3, nfdl1, nfdl2,  &
     &       nfdrep, nfdrep1, fracdec, hmaxnit, pr_nc_period, tdt,  &
     &       do_cell_chem, irma, irmb, irmc, jreorder, jphotrat,  &
     &       ntspec, inewold, denair, corig, pratk1, yemis, smvdm,  &
     &       nfdh1, errmx2, cc2, cnew, gloss, vdiag, rrate, trate,  &
     &       urate, yda, qqkda, qqjda, qkgmi, qjgmi, i1, i2, ju1, j2, &
             k1, k2, num_qjo, num_qks, num_qjs, num_active)

!         -----------------------------------------------------
!         Replace block concentrations (molec/cm^3) into domain
!         concentrations, but only for active species.
!         -----------------------------------------------------

          if (ireord == 2) then

            do jnew = 1, ntspec(ncs)
              jgas = inewold(jnew,ncs)

              if (jgas <= SK_NACT) then

                do kloop = 1, ktloop
                  jloop          = jreorder(jlooplo+kloop)
                  cx(jloop,jgas) = Max (cnew(kloop,jnew), SMAL2)
                end do

              end if

            end do

          end if

!       ======
        end if
!       ======

!c      wclk2       = Mpi_Wtime (ierr)
!c      ctime(kblk) = wclk2 - wclk1
!c      Write (6,910) proc_num, thrd_num, kblk, ctime(kblk)
!910    format ('Proc #, Thrd, Block, Time:  ', i6, i4, i8, e20.8)
!c      sumtime(thrd_num) = sumtime(thrd_num) + ctime(kblk)

!     ======
      end do
!     ======

!c    call f_hpmstop (1)

!c    Write (6,920)
!c   &  proc_num, (sumtime(thrd_num), thrd_num=0,num_thrds-1)
!920  format ('Proc #, sumtime:  ', i6, 16e20.8)

      return

      end subroutine Solve_Block


!-----------------------------------------------------------------------------
!
! ROUTINE
!   Calcrate
!
! DESCRIPTION
!   This routine calculates kinetic reaction and photorates
!   (s^-1, cm^3/s, or cm^6/s^2) and pressure- and temperature- dependence
!   for gas-phase chemical reactions.
!
! ARGUMENTS
!   ifsun     : 1 => for daytime, 2 => for nighttime
!   imgas     : array index for air density
!   initrogen : identifies spc # of nitrogen gas
!   ioxygen   : identifies spc # of oxygen   gas
!   itloop    : # of zones (ilong * ilat * ivert)
!   jlooplo   : low ntloop grid-cell - 1 in a grid-block
!   ktloop    : # of grid-cells in a grid-block
!   ncs       : identifies gas chemistry type (1..NCSGAS)
!   jreorder  : gives original grid-cell from re-ordered grid-cell
!   nrates    : # of kinetic rxns (non-photo)
!   ntspec    : # of active + inactive gases
!   cx        : spc conc (molec/cm^3)
!   urate     : term of Jacobian (J) = partial derivative
!   denair    : density of air  (molec/cm^3)
!   pratk1    : tbd
!   cblk      : gas-phase concs (molec/cm^3)
!   rrate     : rate constants
!               (s^-1, cm^3/molec*s, cm^6/molec*s^2, or cm^9/molec*s^3)
!   trate     : rxn rate (moles l^-1-h2o s^-1 or # cm^-3 s^-1 (?))
!   nallr     : # of active rxns
!   nfdh1     : nfdh2 + # of rxns with one   active reactant
!   nfdh2     : nfdh3 + # of rxns with two   active reactants
!   nfdh3     :         # of rxns with three active reactants
!   nfdl1     : nfdh2 + 1
!   nfdl2     : nfdh3 + 1
!   nfdrep    : nfdh3 + # of rxns with two active reactants that are not
!               followed by a rxn with the same reactants
!   nfdrep1   : nfdrep + 1
!   irma,b,c  : spc # of each reactant; locates reordered active spc #s
!   corig     : original gas-phase concs used to restart Smvgear if a failure
!               occurs (molec/cm^3)
!   smvdm     : amount added to each spc at each grid-cell, for mass balance
!               accounting (# cm^-3 for gas chemistry (?))
!
!-----------------------------------------------------------------------------

      subroutine CalcrateGMI  &
     &  (ifsun, imgas, initrogen, ioxygen, itloop, jlooplo, ktloop,  &
     &   ncs, jreorder, nrates, ntspec, cx, urate, denair, pratk1,  &
     &   cblk, rrate, trate, nallr, nfdh1, nfdh2, nfdh3, nfdl1, nfdl2,  &
     &   nfdrep, nfdrep1, irma, irmb, irmc, corig, smvdm, savedVars)


      use GmiSavedVariables_mod, only : t_ChemistrySaved
      use GmiPrintError_mod, only : GmiPrintError

      implicit none

#     include "smv2chem_par.h"
!#     include "smv2chem2.h"

!     ----------------------
!     Argument declarations.
!     ----------------------

      integer, intent(in)  :: ifsun
      integer, intent(in)  :: imgas
      integer, intent(in)  :: initrogen
      integer, intent(in)  :: ioxygen
      integer, intent(in)  :: itloop
      integer, intent(in)  :: jlooplo
      integer, intent(in)  :: ktloop
      integer, intent(in)  :: ncs
      integer, intent(in)  :: jreorder(itloop)
      integer, intent(in)  :: nrates  (ICS)
      integer, intent(in)  :: ntspec  (ICS)
      real*8,  intent(in)  :: cx      (itloop, IGAS)
      real*8,  intent(in)  :: urate   (KBLOOP, NMTRATE, 3)

      real*8,  intent(inout) :: denair(ktloop)
      real*8,  intent(inout) :: pratk1(KBLOOP, IPHOT)
      real*8,  intent(inout) :: cblk  (KBLOOP, MXGSAER)
      real*8,  intent(inout) :: rrate (KBLOOP, NMTRATE)
      real*8,  intent(inout) :: trate (KBLOOP, NMTRATE*2)

      integer, intent(out) :: nallr
      integer, intent(out) :: nfdh1,  nfdh2, nfdh3
      integer, intent(out) :: nfdl1,  nfdl2
      integer, intent(out) :: nfdrep, nfdrep1
      integer, intent(out) :: irma (NMTRATE)
      integer, intent(out) :: irmb (NMTRATE)
      integer, intent(out) :: irmc (NMTRATE)
      real*8,  intent(out) :: corig(KBLOOP, MXGSAER)
      real*8,  intent(out) :: smvdm(KBLOOP, MXGSAER)

      type(t_ChemistrySaved), intent(inOut) :: savedVars


!     ----------------------
!     Variable declarations.
!     ----------------------

      integer :: i, j
      integer :: jloop
      integer :: jnew
      integer :: jold            ! = mappl(jold) for inactive spc
      integer :: kloop
      integer :: ncsp            ! ncs       => for daytime   gas chemistry
                                 ! ncs + ICS => for nighttime gas chemistry
      integer :: nfdl0           ! nfdh1 + 1
      integer :: nh, nk, nkn

      real*8  :: concn2(ktloop)  ! nitrogen conc  (molec/cm^3)
      real*8  :: conco2(ktloop)  ! oxygen   conc  (molec/cm^3)


!     ----------------
!     Begin execution.
!     ----------------

!     ----------------------------
!     Load kinetic reaction rates.
!     ----------------------------

      do nk = 1, nrates(ncs)
        do kloop = 1, ktloop
          rrate(kloop,nk) = urate(kloop,nk,1)
        end do
      end do


!     -------------------------------------------------------------------
!     Place large domain gas array (molec/cm^3) into smaller block array.
!     -------------------------------------------------------------------

      do jold = 1, ntspec(ncs)

        jnew = savedVars%mappl(jold,ncs)

        do kloop = 1, ktloop

          jloop = jreorder(jlooplo+kloop)

          cblk (kloop,jold) = cx(jloop,jold)
          corig(kloop,jnew) = cx(jloop,jold)
          smvdm(kloop,jnew) = 0.0d0

!         ---------------------------------
!         Load third body number densities.
!         ---------------------------------

          if (initrogen > 0) concn2(kloop) = cx(jloop,initrogen)
          if (ioxygen   > 0) conco2(kloop) = cx(jloop,ioxygen)
          if (imgas     > 0) denair(kloop) = cx(jloop,imgas)

        end do

      end do


!     ------------------------------------------------------
!     Multiply rates by constant species concentrations
!     (either M, O2, N2, or any active or inactive species).
!     ------------------------------------------------------

      do i = 1, savedVars%nmair(ncs)
        nk = savedVars%nreacair(i,ncs)
        do kloop = 1, ktloop
           rrate(kloop,nk) = rrate(kloop,nk) * denair(kloop)
        end do
      end do

      do i = 1, savedVars%nmo2(ncs)
        nk = savedVars%nreaco2(i,ncs)
        do kloop = 1, ktloop
          rrate(kloop,nk) = rrate(kloop,nk) * conco2(kloop)
        end do
      end do

      do i = 1, savedVars%nmn2(ncs)
        nk = savedVars%nreacn2(i,ncs)
        do kloop = 1, ktloop
          rrate(kloop,nk) = rrate(kloop,nk) * concn2(kloop)
        end do
      end do


!     ---------------------------------------------------------------
!     Multiply rate coefficient by any other third body concentration
!     (e.g., H2O); multiply by other inactive concentrations later.
!     ---------------------------------------------------------------

      do i = 1, savedVars%nm3bod(ncs)
        nk   = savedVars%nreac3b (i,ncs)
        jold = savedVars%lgas3bod(i,ncs)
        do kloop = 1, ktloop
          rrate(kloop,nk) = rrate(kloop,nk) * cblk(kloop,jold)
        end do
      end do


!     -----------------------------------------------------------
!     Multiply rate coefficient by other inactive concentrations.
!     This loop must occur after equilibrium reactions.
!     -----------------------------------------------------------

      do i = 1, savedVars%nmoth(ncs)
        nk   = savedVars%nreacoth(i,ncs)
        jold = savedVars%lgasbino(i,ncs)
        do kloop = 1, ktloop
          rrate(kloop,nk) = rrate(kloop,nk) * cblk(kloop,jold)
       end do
      end do


!     --------------------
!     Reorder rrate array.
!     --------------------

      nfdh3   = savedVars%ithrr(ncs)
      nfdl2   = nfdh3  + 1
      nfdrep  = savedVars%inorep(ncs)
      nfdrep1 = nfdrep + 1
      nfdh2   = nfdh3  + savedVars%itwor(ncs)
      nfdl1   = nfdh2  + 1
      nfdh1   = nfdh2  + savedVars%ioner(ncs)
      nfdl0   = nfdh1  + 1
      nallr   = savedVars%nallrat(ncs)

      do nkn = 1, nallr
        nk = savedVars%noldfnew(nkn,ncs)
        irma(nkn) = savedVars%irm2(1,nk,ncs)
        irmb(nkn) = savedVars%irm2(2,nk,ncs)
        irmc(nkn) = savedVars%irm2(3,nk,ncs)
      end do


!     ---------------------------------
!     trate here used as a dummy array.
!     ---------------------------------

      do nk = 1, savedVars%ntrates(ncs)
        do kloop = 1, ktloop
          trate(kloop,nk) = rrate(kloop,nk)
        end do
      end do

      do nkn = 1, nallr
        nk = savedVars%noldfnew(nkn,ncs)
        do kloop = 1, ktloop
          rrate(kloop,nkn) = trate(kloop,nk)
        end do
      end do


!     ----------------------------------------------------------------
!     Set kinetic and photo rates where reaction has no active losses.
!     ----------------------------------------------------------------

      do nkn = nfdl0, nallr
        nh = nkn + nallr
        do kloop = 1, ktloop
          trate(kloop,nkn) =  rrate(kloop,nkn)
          trate(kloop,nh)  = -trate(kloop,nkn)
        end do
      end do


!     ---------------------------------------------------------------
!     Reset ncsp; multiply photorate coef. (s^-1) by concentration
!     (molec/cm^3) when photodissociating species is inactive.  Thus,
!     rate does not need to be multiplied later by concentration.
!     ---------------------------------------------------------------

      ncsp = (ifsun - 1) * ICS + ncs

      do i = 1, savedVars%nolosp(ncsp)
        nk   = savedVars%nknlosp(i,ncs)
        j    = savedVars%jphotnk(nk,ncs)
        jold = savedVars%losinacp(i,ncs)
        do kloop = 1, ktloop
          pratk1(kloop,j)  = pratk1(kloop,j) * cblk(kloop,jold)
        end do
      end do

      return

      end subroutine CalcrateGMI

!-----------------------------------------------------------------------------
!
! ROUTINE
!   Reorder_Grid_Cells
!
! DESCRIPTION
!   This routine reorders the  grid-cells from least to most stiff, with
!   those with similar stiffness grouped together.
!
!   From least to most stiff:
!     smaller errmx2 (csuma) => less stiff
!     (csuma = errmx2)
!
!   Sort using heapsort routine (numerical recipes), an n(logb2)n process.
!   This reordering scheme is very fast, although complicated.  errmx2 from
!   Smvgear:  denotes stiffness (larger value => more stiff).
!
! ARGUMENTS
!   itloop    : # of zones (ilong * ilat * ivert)
!   kuloop    : intended # of grid-cells in a grid-block
!   ntloopuse : tbd
!   errmx2    : tbd
!   jreorder  : gives original grid-cell from re-ordered grid-cell
!   csuma     : tbd
!   nreblock  : tbd
!   lreorder  : gives original grid-cell from re-ordered cell, except when
!               cell is a virtual boundary cell, then it gives original
!               edge cell from re-ordered virtual boundary cell
!   jlowvar   : tbd
!   ktlpvar   : tbd
!   csumc     : tbd
!
!-----------------------------------------------------------------------------

      subroutine Reorder_Grid_Cells  &
     &  (itloop, kuloop, ntloopuse, errmx2, jreorder, csuma,  &
     &   nreblock, lreorder, jlowvar, ktlpvar, csumc)

      implicit none

#     include "smv2chem_par.h"

!     ----------------------
!     Argument declarations.
!     ----------------------

      integer, intent(in)  :: itloop
      integer, intent(in)  :: kuloop
      integer, intent(in)  :: ntloopuse
      real*8,  intent(in)  :: errmx2(itloop)

      integer, intent(inout) :: jreorder(itloop)

      integer, intent(out) :: nreblock
      integer, intent(out) :: lreorder(itloop)
      integer, intent(out) :: jlowvar (MXBLOCK)
      integer, intent(out) :: ktlpvar (MXBLOCK)
      real*8,  intent(out) :: csuma(itloop)
      real*8,  intent(out) :: csumc(itloop)


!     ----------------------
!     Variable declarations.
!     ----------------------

      integer :: iavblok, iavgsize
      integer :: ipar, jpar, jpar1
      integer :: iradd, iradd1
      integer :: iremain
      integer :: irval, lval
      integer :: jllast
      integer :: jloop
      integer :: jlooplo    ! low ntloop grid-cell - 1 in a grid-block
      integer :: jreord
      integer :: kblk
      integer :: nblock1
      integer :: nblockrow  ! # blocks for each reorder group
                            ! (stiffness, sunrise, sunset)
      integer :: ncellrow
      integer :: nnorise

      real*8  :: vallow


!     ----------------
!     Begin execution.
!     ----------------

      do jloop = 1, ntloopuse
        lreorder(jloop) = jreorder(jloop)
        csumc   (jloop) = errmx2  (jloop)
      end do


      jllast = ntloopuse

      do jloop = 1, ntloopuse
        jreorder(jllast) = lreorder(jloop)
        csuma   (jllast) = csumc   (jloop)
        jllast           = jllast - 1
      end do


      nnorise = ntloopuse

      iradd   = 0
      lval    = iradd + (nnorise * 0.5d0) + 1
      irval   = iradd + nnorise
      iradd1  = iradd + 1

!     ===================
      if (irval > iradd1) then
!     ===================

!       =============
        OUTERLOOP: do
!       =============

          if (lval > iradd1) then

            lval            = lval - 1
            vallow          = csuma   (lval)
            jreord          = jreorder(lval)

          else

            vallow          = csuma   (irval)
            jreord          = jreorder(irval)
            csuma   (irval) = csuma   (iradd1)
            jreorder(irval) = jreorder(iradd1)
            irval           = irval - 1

            if (irval == iradd1) then
              csuma   (iradd1) = vallow
              jreorder(iradd1) = jreord
!             ==============
              exit OUTERLOOP
!             ==============
            end if

          end if

          ipar = lval
          jpar = lval + lval - iradd

!         ===================================
          INNERLOOP: do while (jpar <= irval)
!         ===================================

            if (jpar < irval) then
              jpar1 = jpar + 1
              if (csuma(jpar) < csuma(jpar1)) jpar = jpar1
            end if

            if (vallow < csuma(jpar)) then
              csuma   (ipar) = csuma   (jpar)
              jreorder(ipar) = jreorder(jpar)
              ipar           = jpar
              jpar           = jpar + jpar - iradd
!             ===============
              cycle INNERLOOP
!             ===============
            else
!             ==============
              exit INNERLOOP
!             ==============
            end if

!         ================
          end do INNERLOOP
!         ================

          csuma   (ipar) = vallow
          jreorder(ipar) = jreord

!       ================
        end do OUTERLOOP
!       ================

!     ======
      end if
!     ======


!     ---------------------------------------------------------------
!     Determine how many blocks of cells are needed after reordering.
!     ---------------------------------------------------------------

      jlooplo   = 0
      nblockrow = 0
      nreblock  = 0

      ncellrow = nnorise

      if (ncellrow /= 0) then

        nblockrow = 1 + (ncellrow / (kuloop + 0.0001d0))

        iavblok   = 1 + (ncellrow / (nblockrow + 0.0001d0))
        iavgsize  = Min (iavblok, kuloop)
        nblock1   = nblockrow - 1

        do kblk = 1, nblock1
          nreblock          = nreblock + 1
          jlowvar(nreblock) = jlooplo
          ktlpvar(nreblock) = iavgsize
          jlooplo           = jlooplo + iavgsize
        end do

        nreblock = nreblock + 1
        iremain  = Max (ncellrow - (nblock1 * iavgsize), 0)

        jlowvar(nreblock) = jlooplo
        ktlpvar(nreblock) = iremain

      end if

      return

      end subroutine Reorder_Grid_Cells

