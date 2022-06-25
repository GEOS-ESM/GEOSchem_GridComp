
!=============================================================================
!
! $Id$
!
! CODE DEVELOPER
!   Original code from Peter Connell, LLNL
!   Gmimod modifications:  John Tannahill
!                          jrt@llnl.gov
!
! FILE
!   smv2chem_diag.F
!
! ROUTINES
!   Do_Smv2_Diag
!
!=============================================================================


!-----------------------------------------------------------------------------
!
! ROUTINE
!   Do_Smv2_Diag
!
! DESCRIPTION
!   This routine collects the Smvgear II chemical diagnostics.  It accumulates
!   averages for species and reaction rates.
!
! ARGUMENTS
!   jlooplo      : low ntloop grid-cell - 1 in a grid-block
!   ktloop       : number of grid-cells     in a grid-block
!   pr_nc_period : NetCDF output period
!   tdt          : model time step (s)
!   told         : stores last value of xelaps in case current step fails
!   do_cell_chem : do chemistry for a particular cell?
!   jreorder     : gives original grid-cell from re-ordered grid-cell
!   inewold      : original spc # of each new jnew spc
!   denair       : density of air (molec/cm^3)
!   cnew         : init (and final) spc conc
!                  (# cm^-3-air or moles l^-1-h2o (?))
!   xtimestep    : xelaps - told
!
!-----------------------------------------------------------------------------

      subroutine Do_Smv2_Diag  &
     &  (jlooplo, ktloop, pr_nc_period, tdt, told, do_cell_chem,  &
     &   jreorder, inewold, denair, cnew, xtimestep, &
     &   yda, qqkda, qqjda, qkgmi, qjgmi, &
     &   ilong, ilat, ivert, itloop, &
     &   i1, i2, ju1, j2, k1, k2, &
     &   num_qjo, num_qks, num_qjs, num_active)

#     include "smv2chem_par.h"


!     ----------------------
!     Argument declarations.
!     ----------------------

      integer, intent(in) :: i1, i2, ju1, j2, k1, k2
      integer, intent(in) :: ilong, ilat, ivert, itloop
      integer, intent(in) :: num_qjo, num_qks, num_qjs, num_active
      real*8 , intent(in   ) :: qjgmi(i1:i2, ju1:j2, k1:k2, num_qjo)
      real*8 , intent(in   ) :: qkgmi(i1:i2, ju1:j2, k1:k2, num_qks)
      real*8 , intent(inout) :: qqjda(i1:i2, ju1:j2, k1:k2, num_qjs)
      real*8 , intent(inout) :: qqkda(i1:i2, ju1:j2, k1:k2, num_qks)
      real*8 , intent(inout) :: yda  (i1:i2, ju1:j2, k1:k2, num_active)

      integer, intent(in)  :: jlooplo
      integer, intent(in)  :: ktloop
      real*8,  intent(in)  :: pr_nc_period
      real*8,  intent(in)  :: tdt
      real*8,  intent(in)  :: told
      logical, intent(in)  :: do_cell_chem(ilong, ilat, ivert)
      integer, intent(in)  :: jreorder(itloop)
      integer, intent(in)  :: inewold (MXGSAER, ICS)
      real*8,  intent(in)  :: denair  (ktloop)
      real*8,  intent(in)  :: cnew    (KBLOOP, MXGSAER)

      real*8,  intent(inout) :: xtimestep



      return

      end

