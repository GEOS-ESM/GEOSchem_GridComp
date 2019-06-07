
!=============================================================================
!
! $Id$
!
! CODE DEVELOPER
!   John Tannahill, LLNL (Original code from Mark Jacobson, UCLA/Stanford)
!   jrt@llnl.gov
!
! FILE
!   smv2chem1.h
!
! DESCRIPTION
!   This include file contains the common blocks for SmvgearII for the
!   variables that need to be saved between calls.
!
!   Note that during analysis for adding openmp to the Physproc routine,
!   smv2chem.h was split into smv2chem1.h & smv2chem2.h, to make the handling
!   of data cleaner and clearer.
!
!=============================================================================


!     ----------------
!     Integer commons.
!     ----------------

!     ----------------------------------------------------------------
!     ifreord   : if 1, then reorder grid-cells by stiffness
!     ih2o      : identifies spc # of water vapor
!     imgas     : tbd
!     initrogen : identifies spc # of nitrogen gas
!     ioxygen   : identifies spc # of oxygen   gas
!     kuloop    : intended # of grid-cells in a grid-block
!     lunsmv    : logical unit number to write to when pr_smv2 is true
!     ncs       : identifies gas chemistry type (1..NCSGAS)
!     ----------------------------------------------------------------

      integer :: ifreord
      integer :: ih2o
      integer :: imgas
      integer :: initrogen, ioxygen
      integer :: kuloop
      integer :: lunsmv
      integer :: ncs

!     ====================
      common  / smv21_i1 /  &
!     ====================
     &  ifreord,  &
     &  ih2o,  &
     &  imgas,  &
     &  initrogen, ioxygen,  &
     &  kuloop,  &
     &  lunsmv,  &
     &  ncs


!     -----------------------------------------
!     jphotrat  : tbd
!     nrates    : # of kinetic rxns (non-photo)
!     ntloopncs : tbd
!     ntspec    : # of active + inactive gases
!     -----------------------------------------

      integer :: jphotrat (ICS)
      integer :: nrates   (ICS)
      integer :: ntloopncs(ICS)
      integer :: ntspec   (ICS)

!     ====================
      common  / smv21_i2 /  &
!     ====================
     &  jphotrat,  &
     &  nrates,  &
     &  ntloopncs,  &
     &  ntspec


!     -----------------------------------------------
!     inewold   : original spc # of each new jnew spc
!
!     npphotrat : tbd
!     -----------------------------------------------

      integer :: inewold  (MXGSAER, ICS)

      integer :: npphotrat(IPHOT,   ICS)

!     ====================
      common  / smv21_i3 /  &
!     ====================
     &  inewold,  &
!
     &  npphotrat


!     -------------
!     Real commons.
!     -------------

!     -------------------------------------------------------------------
!     fracdec : fraction time step is decreased in Smvgear if convergence
!               test fails
!     hmaxnit : max time step for night all chem (s)
!     -------------------------------------------------------------------

      real*8  :: fracdec
      real*8  :: hmaxnit

!     ====================
      common  / smv21_r1 /  &
!     ====================
     &  fracdec,  &
     &  hmaxnit

