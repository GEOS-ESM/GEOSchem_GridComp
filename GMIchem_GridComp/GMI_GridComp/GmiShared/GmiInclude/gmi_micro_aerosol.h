
!=============================================================================
!
! $Id$
!
! CODE DEVELOPER
!   Dan Bergmann, LLNL
!   dbergmann@llnl.gov
!
! FILE
!   imp_aerosol.h
!
! DESCRIPTION
!   This include file contains information about various aerosols
!
!=============================================================================


!     -----------------------
!     Parameter declarations.
!     -----------------------

!     ----------------------------------------------------------------------
!     NUM_AEROSOL : # of different aerosol types having settling velocities:
!                   1) sulfate 1  (0.01 - 0.05 um)
!                   2) sulfate 2  (0.05 - 0.63 um)
!                   3) natural organic
!                   4) fossil fuel OC
!                   5) fossil fuel BC
!                   6) biomass OC
!                   7) biomass BC
!                   8) dust 1     (0.05 - 0.63 um)
!                   9) dust 2     (0.63 - 1.26 um)
!                  10) dust 3     (1.26 - 2.5  um)
!                  11) dust 4     (2.5  - 10   um)
!                  12) sea salt 1 (0.05 - 0.63 um)
!                  13) sea salt 2 (0.63 - 1.26 um)
!                  14) sea salt 3 (1.26 - 2.5  um)
!                  15) sea salt 4 (2.5  - 10   um)
!                  16) stratospheric sulfate (mono-dispersal)
!                  17) aircraft soot
!     ----------------------------------------------------------------------

      integer, parameter :: NUM_AEROSOL = 17
      integer, parameter :: NUM_WAV_AER = 1


!     ----------------------------------------------------
!     AER_DENSITY : density of each aerosol type (kg/m^3).
!     ----------------------------------------------------

      real*8, parameter :: AER_DENSITY(NUM_AEROSOL) =  &
     &  (/ 1769.d0, 1769.d0, 1500.d0, 1500.d0, 1500.d0, 1500.d0,  &
     &     1500.d0, 2500.d0, 2500.d0, 2650.d0, 2650.d0, 2216.d0,  &
     &     2216.d0, 2216.d0, 2216.d0, 1769.d0, 1500.d0 /)

!     --------------------------------------------------------
!     RADIUS_EFF : effective radius of each aerosol type (m).
!     --------------------------------------------------------

      real*8, parameter :: RADIUS_EFF(NUM_AEROSOL) =  &
     &  (/ 0.0405d-6, 0.2695d-6, 0.2064d-6, 0.2130d-6, 0.2130d-6,  &
     &     0.2064d-6, 0.2064d-6, 0.4610d-6, 0.9868d-6, 1.9140d-6,   & ! summer dust
     &     5.8800d-6, 0.4040d-6, 0.9759d-6, 1.9420d-6, 7.1030d-6,  &
     &     0.1000d-6, 0.0444d-6 /)

!     --------------------------------------------------------
!     Constants used in swelling calculation;
!     to turn off swelling for a particular aerosol, set C1=0.
!     --------------------------------------------------------

      real*8, parameter :: C1(NUM_AEROSOL) =  &
     &  (/  0.4809d0,  0.4809d0,  0.2789d0,  0.3926d0,  0.3926d0,  &
     &      0.2789d0,  0.2789d0,  0.0000d0,  0.0000d0,  0.0000d0,  &
     &      0.0000d0,  0.7674d0,  0.7674d0,  0.7674d0,  0.7674d0,  &
     &      0.4809d0,  0.4809d0 /)

      real*8, parameter :: C2(NUM_AEROSOL) =  &
     &  (/  3.082d0,   3.082d0,   3.115d0,   3.101d0,   3.101d0,  &
     &      3.115d0,   3.115d0,   0.000d0,   0.000d0,   0.000d0,  &
     &      0.000d0,   3.079d0,   3.079d0,   3.079d0,   3.079d0,  &
     &      3.082d0,   3.082d0 /)

      real*8, parameter :: C3(NUM_AEROSOL) =  &
     &  (/  3.110d-11, 3.110d-11, 5.415d-11, 4.190d-11, 4.190d-11,  &
     &      5.415d-11, 5.415d-11, 0.000d-00, 0.000d-00, 0.000d-00,  &
     &      0.000d-00, 2.573d-11, 2.573d-11, 2.573d-11, 2.573d-11,  &
     &      3.110d-11, 3.110d-11 /)

      real*8, parameter :: C4(NUM_AEROSOL) =  &
     &  (/ -1.428d0,  -1.428d0,  -1.399d0,  -1.404d0,  -1.404d0,  &
     &     -1.399d0,  -1.399d0,   0.000d0,   0.000d0,   0.000d0,  &
     &      0.000d0,  -1.424d0,  -1.424d0,  -1.424d0,  -1.424d0,  &
     &     -1.428d0,  -1.428d0 /)


!     ----------------------------------------------------------------
!     REL_SCAV_EFF : scavenging efficiency of each aerosol relative to
!                    sulfate (unitless)
!     ----------------------------------------------------------------

      real*8, parameter :: REL_SCAV_EFF(NUM_AEROSOL) =  &
     &  (/ 0.0d0, 1.0d0, 0.4d0, 0.4d0, 0.4d0, 0.4d0, 0.4d0,  &
     &     1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0,  &
     &     1.0d0, 1.0d0, 1.0d0 /)

