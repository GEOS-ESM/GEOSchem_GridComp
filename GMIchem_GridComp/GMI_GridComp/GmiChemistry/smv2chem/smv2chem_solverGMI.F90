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
!   smv2chem_solver.F
!
! ROUTINES
!   Do_Smv2_Solver
!
!=============================================================================


!-----------------------------------------------------------------------------
!
! ROUTINE
!   Do_Smv2_Solver
!
! DESCRIPTION
!   This is the main control routine for the ordinary differential equation
!   solver, "Smvgear II" (Sparse Matrix Vectorized Gear-type code).
!
! ARGUMENTS
!   do_qqjk_inchem   : if pr_qqjk is on, should qqj's & qqk's be determined
!                      inside the chemistry solver, or outside?
!   do_semiss_inchem : do surface emissions inside the chemistry solver, or
!                      outside?
!   pr_diag          : print some diagnostic output to screen?
!   pr_qqjk          : should the periodic qqjk output file be written?
!   pr_smv2          : should the SmvgearII     output file be written
!                      (non-parallel mode only)?
!   loc_proc         : local processor #
!   ilat             : # of latitudes
!   ilong            : # of longitudes
!   ivert            : # of vertical layers
!   itloop           : # of zones (ilong * ilat * ivert)
!   pr_nc_period     : NetCDF output period
!   tdt              : model time step (s)
!   do_cell_chem     : do chemistry for a particular cell?
!   arate            : thermal    rate constants (units vary)
!   prate            : photolysis rate constants (s^-1)
!   yemis            : surface emissions (molec/cm^3/s)
!   cx               : spc conc (molec/cm^3)
!
!-----------------------------------------------------------------------------

      subroutine Do_Smv2_Solver (savedVars, do_qqjk_inchem, &
                    do_semiss_inchem, pr_diag, pr_qqjk, pr_smv2,  &
                    loc_proc, ilat, ilong, ivert, itloop, pr_nc_period,&
                    tdt,  do_cell_chem, arate, prate, yemis, cx, &
                    yda, qqkda, qqjda, qkgmi, qjgmi, i1, i2, ju1, j2, &
                    k1, k2, num_qjo, num_qks, num_qjs, num_active)

      use GmiSavedVariables_mod, only : t_ChemistrySaved

      IMPLICIT NONE

#     include "smv2chem_par.h"
!#     include "smv2chem1.h"

!     ----------------------
!     Argument declarations.
!     ----------------------

      type(t_ChemistrySaved), intent(inOut) :: savedVars
      integer, intent(in) :: i1, i2, ju1, j2, k1, k2
      integer, intent(in) :: num_qjo, num_qks, num_qjs, num_active
      real*8 , intent(in) :: qjgmi(i1:i2, ju1:j2, k1:k2, num_qjo)
      real*8 , intent(in) :: qkgmi(i1:i2, ju1:j2, k1:k2, num_qks)
      logical, intent(in) :: do_qqjk_inchem
      logical, intent(in) :: do_semiss_inchem
      logical, intent(in) :: pr_diag
      logical, intent(in) :: pr_qqjk
      logical, intent(in) :: pr_smv2
      integer, intent(in) :: loc_proc
      integer, intent(in) :: ilat, ilong, ivert
      integer, intent(in) :: itloop
      real*8,  intent(in) :: pr_nc_period
      real*8,  intent(in) :: tdt
      logical, intent(in) :: do_cell_chem(itloop)
      real*8,  intent(in) :: arate(itloop, ITHERM)
      real*8,  intent(in) :: prate(itloop, IPHOT)
      real*8,  intent(in) :: yemis(ilat*ilong, IGAS)

      real*8 , intent(inout) :: qqjda(i1:i2, ju1:j2, k1:k2, num_qjs)
      real*8 , intent(inout) :: qqkda(i1:i2, ju1:j2, k1:k2, num_qks)
      real*8 , intent(inout) :: yda  (i1:i2, ju1:j2, k1:k2, num_active)
      real*8,  intent(inout) :: cx(itloop, IGAS)

!     ----------------------
!     Variable declarations.
!     ----------------------

      integer :: jreorder(itloop)
      integer :: lreorder(itloop)

      real*8  :: errmx2  (itloop)

!     ----------------
!     Begin execution.
!     ----------------

      if (pr_diag) then
        Write (6,*) 'Do_Smv2_Solver called by ', loc_proc
      end if


      jreorder(:) = 0
      lreorder(:) = 0

      errmx2  (:) = 0.0d0

      call PhysprocGMI(savedVars, do_qqjk_inchem, do_semiss_inchem, &
               pr_qqjk, pr_smv2, ilat, ilong, ivert, savedVars%ifreord,&
               savedVars%imgas, savedVars%initrogen, savedVars%ioxygen,&
               itloop, savedVars%kuloop, savedVars%lunsmv, &
               savedVars%ncs, savedVars%fracdec, savedVars%hmaxnit, &
               pr_nc_period, tdt, do_cell_chem, savedVars%jphotrat, &
               savedVars%nrates, savedVars%ntloopncs, savedVars%ntspec, &
               savedVars%inewold, savedVars%npphotrat, arate, prate, &
               yemis, jreorder, lreorder, savedVars%csuma, &
               savedVars%csumc, errmx2, cx, yda, qqkda, qqjda, qkgmi, &
               qjgmi, i1, i2, ju1, j2, k1, k2, num_qjo, num_qks, &
               num_qjs, num_active)

      return

      end subroutine Do_Smv2_Solver
