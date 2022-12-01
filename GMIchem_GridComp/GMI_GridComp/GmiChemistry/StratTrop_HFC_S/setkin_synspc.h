!=======================================================================
!
! $Id: setkin_synspc.h $
!
! CODE DEVELOPER
!   Peter Connell, LLNL
!   connell2@llnl.gov
!
! FILE
!   setkin_synspc.h
!
! DESCRIPTION
!   This include file contains information about synthetic stratospheric
!   sources for ozone and nitrogen oxides.
!
!  Input mechanism:        GeosCCM_Combo_2015mechanism.txt
!  Reaction dictionary:    GMI_reactions_JPL15.db
!  Setkin files generated: Tue Mar 30 16:58:02 2021
!
!=======================================================================


      logical, parameter :: USE_SYNOZ = .false.
      logical, parameter :: USE_NODOZ = .false.

      integer, parameter :: MAXNODOZ_ELEM = 10

      integer, parameter :: NOX_ELEM_MAP(MAXNODOZ_ELEM) = &
     &  (/ INO, INO2, INO3, IN2O5, IN2O5,  0,  0,  0,  0,  0 /)

      integer, parameter :: NOY_ELEM_MAP(MAXNODOZ_ELEM) = &
     &  (/  IHNO2,  IHNO3,  IHNO4,  0,  0,  0,  0,  0,  0,  0 /)

!                                  --^--

