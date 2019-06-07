
!=============================================================================
!
! $Id$
!
! CODE DEVELOPER
!   Al Franz, LLNL
!   franz2@llnl.gov
!
! FILE
!   gmi_diag_constants_llnl.h
!
! DESCRIPTION
!   This include file contains parameter declarations for the LLNL
!   production and loss diagnostics.
!
!=============================================================================


      integer, parameter :: BEFORE        =  1
      integer, parameter :: AFTER         =  2

      integer, parameter :: ONE_STEP      =  1
      integer, parameter :: CUMULATIVE    =  2

      integer, parameter :: EMISS_OP      =  1
      integer, parameter :: DIFFU_OP      =  2
      integer, parameter :: ADVEC_OP      =  3
      integer, parameter :: SETTLING_OP   =  4
      integer, parameter :: CONVEC_OP     =  5
      integer, parameter :: DRYDEP_OP     =  6
      integer, parameter :: WETDEP_OP     =  7
      integer, parameter :: SIMPDEP_OP    =  8
      integer, parameter :: ADDWAT_OP     =  9
      integer, parameter :: CHEM_OP       = 10
      integer, parameter :: REMWAT_OP     = 11
      integer, parameter :: SYNSPC_OP     = 12

      integer, parameter :: NUM_OPERATORS = 12


      character*16, parameter :: OPERATOR_NAME(NUM_OPERATORS) =  &
     &  (/ "EMISSION        ",  &
     &     "DIFFUSION       ",  &
     &     "ADVECTION       ",  &
     &     "SETTLING        ",  &
     &     "CONVECTION      ",  &
     &     "DRY_DEP         ",  &
     &     "WET_DEP         ",  &
     &     "SIMPLE_DEP      ",  &
     &     "ADD_WATER       ",  &
     &     "CHEMISTRY       ",  &
     &     "REMOVE_WATER    ",  &
     &     "SYNTHETIC SPCS  "/)

