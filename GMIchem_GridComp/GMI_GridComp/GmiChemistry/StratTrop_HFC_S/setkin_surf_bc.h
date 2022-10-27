!=======================================================================
!
! $Id: setkin_surf_bc.h $
!
! CODE DEVELOPER
!   Peter Connell, LLNL
!   connell2@llnl.gov
!
! FILE
!   setkin_surf_bc.h
!
! DESCRIPTION
!   This include file contains information about surface boundary
!   conditions.
!
!  Input mechanism:        GeosCCM_Combo_2015mechanism.txt
!  Reaction dictionary:    GMI_reactions_JPL15.db
!  Setkin files generated: Tue Mar 30 16:58:02 2021
!
!=======================================================================

!     -----------------------
!     Parameter declarations.
!     -----------------------

!     -------------------------------------------------------
!     NUM_SBC : number of surface bounday conditions to reset
!     K_SBC   : max k to which species will be adjusted
!     -------------------------------------------------------

      integer, parameter :: NUM_SBC = 1

      integer, parameter :: K_SBC   = 2

!     ------------------
!     Integer variables.
!     ------------------

      integer, save :: sbc_map(NUM_SBC) = &
     &  (/ 0 /)

!                                  --^--

