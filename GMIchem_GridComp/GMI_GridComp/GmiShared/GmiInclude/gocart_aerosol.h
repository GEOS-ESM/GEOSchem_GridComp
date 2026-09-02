
!=============================================================================
!
! $Id$
!
! CODE DEVELOPER
!   Dan Bergmann, LLNL
!   dbergmann@llnl.gov
!
! FILE
!   gmi_aerosol.h
!
! DESCRIPTION
!   This include file contains information about various aerosols.
!
!=============================================================================


!     -----------------------
!     Parameter declarations.
!     -----------------------

!     ---------------------------------------------------------------------
!     NUM_AEROSOL : # of different aerosol types having settling velocities
!                   and they may have different optical depth
!                   characteristics =>
!                    1) sulfate
!                    2) organic carbon from biomass burning
!                    3) black   carbon from biomass burning
!                    4) organic carbon from fossil fuel
!                    5) black   carbon from fossil fuel
!                    6) organic carbon from natural sources
!                    7) sea salt 1    (sub-micron)
!                    8) sea salt 2    (super-micron)
!                    9) winter dust 1 (sub-micron)
!                   10) winter dust 2 (super-micron)
!                   11) summer dust 1 (sub-micron)
!                   12) summer dust 2 (super-micron)
!                   13) desert dust 1 (sub-micron)
!                   14) desert dust 2 (super-micron)
!                   15) stratospheric sulfate (monodispersed)
!                   16) dust size #1 from GOCART
!                   17) dust size #2 from GOCART
!                   18) dust size #3 from GOCART
!                   19) dust size #4 from GOCART
!                   20) sulfate 1  (0.01 - 0.05 um) UMich aerosol
!                   21) sulfate 2  (0.05 - 0.63 um) UMich aerosol
!                   22) sulfate 3  (0.63 - 1.26 um) UMich aerosol
!                   23) sea salt 1 (0.05 - 0.63 um) UMich aerosol
!                   24) sea salt 2 (0.63 - 1.26 um) UMich aerosol
!                   25) sea salt 3 (1.26 - 2.50 um) UMich aerosol
!                   26) sea salt 4 (2.50 - 10   um) UMich aerosol
!                   27) dust 1     (0.05 - 0.63 um) UMich, (0.1-1   um) GOCART
!                   28) dust 2     (0.63 - 1.26 um) UMich, (1.0-1.8 um) GOCART
!                   29) dust 3     (1.26 - 2.50 um) UMich, (1.8-3.0 um) GOCART
!                   30) dust 4     (2.50 - 10   um) UMich, (3.0-6.0 um) GOCART
!                   31) dust 5     aero for UMich aerosol, (6.0-10. um) GOCART
!                   32) biomass and natural organic UMich aerosol
!                   33) fossil fuel                 UMich aerosol
!     ---------------------------------------------------------------------

      integer, parameter :: NUM_AEROSOL = 33
      integer, parameter :: NUM_WAV_AER = 1


!     ---------------------------------------------------
!     AER_DENSITY : density of each aerosol type (kg/m^3)
!     ---------------------------------------------------

      real*8, parameter :: AER_DENSITY(NUM_AEROSOL) =  &
     &  (/ 1769.d0, 1500.d0, 1500.d0, 1500.d0, 1500.d0,  &
     &     1500.d0, 2216.d0, 2216.d0, 2500.d0, 2500.d0,  &
     &     2500.d0, 2500.d0, 2500.d0, 2500.d0, 1769.d0,  &
     &     2500.d0, 2650.d0, 2650.d0, 2650.d0, 1769.d0,  &
     &     1769.d0, 1769.d0, 2200.d0, 2200.d0, 2200.d0,  &
     &     2200.d0, 2500.d0, 2650.d0, 2650.d0, 2650.d0,  &
     &     2650.d0, 1500.d0, 1500.d0 /)


!     ------------------------------------------------------
!     RADIUS_EFF : effective radius of each aerosol type (m)
!     ------------------------------------------------------

      real*8, parameter :: RADIUS_EFF(NUM_AEROSOL) =  &
     &  (/ 0.3417d-6, 0.3900d-7, 0.2064d-6, 0.8700d-7, 0.2130d-6,  &
     &     0.2064d-6, 0.6660d-6, 6.9770d-6, 0.6755d-6, 3.5530d-6,  &
     &     0.7276d-6, 5.1830d-6, 0.7118d-6, 6.5300d-6, 0.1000d-6,  &
     &     0.7500d-6, 1.5000d-6, 2.5000d-6, 4.0000d-6, 0.0405d-6,  &
     &     0.2695d-6, 0.8270d-6, 0.2600d-6, 1.1900d-6, 2.4300d-6,  &
     &     7.5700d-6, 0.7300d-6, 1.4000d-6, 2.3999d-6, 4.5006d-6,  &
     &     8.0006d-6, 0.2064d-6, 0.2130d-6 /)


!     --------------------------------------------------------
!     Constants used in swelling calculation;
!     to turn off swelling for a particular aerosol, set C1=0.
!     --------------------------------------------------------

      real*8, parameter :: C1(NUM_AEROSOL) =  &
     &  (/  0.4809d0,  0.2789d0,  0.2789d0,  0.3926d0,  0.3926d0,  &
     &      0.2789d0,  0.7674d0,  0.7674d0,  0.0000d0,  0.0000d0,  &
     &      0.0000d0,  0.0000d0,  0.0000d0,  0.0000d0,  0.4809d0,  &
     &      0.0000d0,  0.0000d0,  0.0000d0,  0.0000d0,  0.4809d0,  &
     &      0.4809d0,  0.4809d0,  0.7674d0,  0.7674d0,  0.7674d0,  &
     &      0.76764d0, 0.0000d0,  0.0000d0,  0.0000d0,  0.0000d0,  &
     &      0.0000d0,  0.2789d0,  0.3926d0 /)

      real*8, parameter :: C2(NUM_AEROSOL) =  &
     &  (/  3.082d0,   3.115d0,   3.115d0,   3.101d0,   3.101d0,  &
     &      3.115d0,   3.079d0,   3.079d0,   0.000d0,   0.000d0,  &
     &      0.000d0,   0.000d0,   0.000d0,   0.000d0,   3.082d0,  &
     &      0.000d0,   0.000d0,   0.000d0,   0.000d0,   3.082d0,  &
     &      3.082d0,   3.082d0,   3.079d0,   3.079d0,   3.079d0,  &
     &      3.079d0,   0.000d0,   0.000d0,   0.000d0,   0.000d0,  &
     &      0.000d0,   3.115d0,   3.101d0 /)

      real*8, parameter :: C3(NUM_AEROSOL) =  &
     &  (/  3.110d-11, 5.415d-11, 5.415d-11, 4.190d-11, 4.190d-11,  &
     &      5.415d-11, 2.573d-11, 2.573d-11, 0.000d-00, 0.000d-00,  &
     &      0.000d-00, 0.000d-00, 0.000d-00, 0.000d-00, 3.110d-11,  &
     &      0.000d-00, 0.000d-00, 0.000d-00, 0.000d-00, 3.110d-11,  &
     &      3.110d-11, 3.110d-11, 2.573d-11, 2.573d-11, 2.573d-11,  &
     &      2.573d-11, 0.000d-00, 0.000d-00, 0.000d-00, 0.000d-00,  &
     &      0.000d-00, 5.415d-11, 4.190d-11 /)

      real*8, parameter :: C4(NUM_AEROSOL) =  &
     &  (/ -1.428d0,  -1.399d0,  -1.399d0,  -1.404d0,  -1.404d0,  &
     &     -1.399d0,  -1.424d0,  -1.424d0,   0.000d0,   0.000d0,  &
     &      0.000d0,   0.000d0,   0.000d0,   0.000d0,  -1.428d0,  &
     &      0.000d0,   0.000d0,   0.000d0,   0.000d0,  -1.428d0,  &
     &     -1.428d0,  -1.428d0,  -1.424d0,  -1.424d0,  -1.424d0,  &
     &     -1.424d0,   0.000d0,   0.000d0,   0.000d0,   0.000d0,  &
     &      0.000d0,  -1.399d0,  -1.404d0 /)


!     ----------------------------------------------------------------
!     REL_SCAV_EFF : scavenging efficiency of each aerosol relative to
!                    sulfate; all set to 1.0 for now (unitless)
!     ----------------------------------------------------------------

      real*8, parameter :: REL_SCAV_EFF(NUM_AEROSOL) =  &
     &  (/ 1.0d0, 0.4d0, 1.0d0, 0.4d0, 1.0d0,  &
     &     1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0,  &
     &     1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0,  &
     &     1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0,  &
     &     1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0,  &
     &     1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0,  &
     &     1.0d0, 0.4d0, 0.4d0 /)
