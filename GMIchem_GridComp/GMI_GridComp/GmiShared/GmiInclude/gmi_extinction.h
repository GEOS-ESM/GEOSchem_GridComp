
!=============================================================================
!
! $Id$
!
! CODE DEVELOPER
!   Dan Bergmann, LLNL
!   dbergmann@llnl.gov
!
! FILE
!   gmi_extinction.h
!
! DESCRIPTION
!   This include file contains extinction coefficients for various aerosols.
!
!=============================================================================


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
!                   16) dust type #1 from GOCART (sub-micron)
!                   17) dust type #2 from GOCART (super-micron)
!                   18) dust type #3 from GOCART (super-micron)
!                   19) dust type #4 from GOCART (super-micron)
!                   20) sulfate 1  (0.01 - 0.05 um) UMich aerosol
!                   21) sulfate 2  (0.05 - 0.63 um) UMich aerosol
!                   22) sulfate 3  (0.63 - 1.26 um) UMich aerosol
!                   23) sea salt 1 (0.05 - 0.63 um) UMich aerosol
!                   24) sea salt 2 (0.63 - 1.26 um) UMich aerosol
!                   25) sea salt 3 (1.26 - 2.50 um) UMich aerosol
!                   26) sea salt 4 (2.50 - 10   um) UMich aerosol
!                   27) dust 1     (0.05 - 0.63 um) UMich aerosol
!                   28) dust 2     (0.63 - 1.26 um) UMich aerosol
!                   29) dust 3     (1.26 - 2.50 um) UMich aerosol
!                   30) dust 4     (2.50 - 10   um) UMich aerosol
!                   31) biomass and natural organic UMich aerosol
!                   32) fossil fuel                 UMich aerosol
!     ---------------------------------------------------------------------


!     -------------------------------------------------------------------
!     extinction_coef : extinction coefficients for each aerosol (m^2/gm)
!     -------------------------------------------------------------------

      real*8, parameter :: EXTINCTION_COEF(NUM_AEROSOL,NUM_WAV_AER) =  &
     &  Reshape  &
     &    ((/ 4.9903d0, 2.8203d0, 2.8203d0, 4.5867d0, 4.5867d0,  &
     &        2.8203d0, 1.7933d0, 0.1380d0, 1.6239d0, 0.3232d0,  &
     &        1.3704d0, 0.2086d0, 1.4932d0, 0.1452d0, 4.9903d0,  &
     &        1.4932d0, 0.1452d0, 0.1452d0, 0.1452d0, 0.0000d0,  &
     &        0.0000d0, 0.0000d0, 0.0000d0, 0.0000d0, 0.0000d0,  &
     &        0.0000d0, 0.0000d0, 0.0000d0, 0.0000d0, 0.0000d0,  &
     &        0.0000d0, 0.0000d0 /),  &
     &     (/ NUM_AEROSOL,NUM_WAV_AER /))


!     -----------------------------------------------------------
!     ehc_so4 : extinction humidity coefficients used for sulfate
!     -----------------------------------------------------------

      real*8, parameter :: EHC_SO4(5,NUM_WAV_AER) =  &
     &  Reshape  &
     &    ((/ 5.6881d-01, 3.8235d+00, 9.6641d+00, 2.7364d+02,  &
     &        3.8671d-03 /),  &
     &     (/ 5,NUM_WAV_AER /))


!     ------------------------------------------------------------------
!     ehc_bb : extinction humidity coefficients used for biomass burning
!     ------------------------------------------------------------------

      real*8, parameter :: EHC_BB(13,NUM_WAV_AER) =  &
     &  Reshape  &
     &    ((/ 2.46744d+00, -6.41756d-01,  4.03946d-01, -7.76428d-02,  &
     &       1.98641d-02,  4.37936d+00,  1.28083d-01,  1.43126d+01,  &
     &       8.37816d+02,  7.67010d+00,  2.34498d+01,  7.54765d-04,  &
     &       1.40055d-04 /),  &
     &     (/ 13,NUM_WAV_AER /))


!     --------------------------------------------------------------
!     ehc_ff : extinction humidity coefficients used for fossil fuel
!     --------------------------------------------------------------

      real*8, parameter :: EHC_FF(13,NUM_WAV_AER) =  &
     &  Reshape  &
     &    ((/ 9.27459d-01, -3.92346d-01,  5.06035d-01, -2.07598d-01,  &
     &       5.46221d-02,  3.92589d+00,  3.75507d-01,  1.19650d+01,  &
     &       1.05464d+03,  1.13302d+01,  6.72778d+01,  7.89406d-04,  &
     &      -1.79976d-05 /),  &
     &     (/ 13,NUM_WAV_AER /))


!     ----------------------------------------------------------------
!     ehc_ocnat : extinction humidity coefficients used for natural OC
!     ----------------------------------------------------------------

      real*8, parameter :: EHC_OCNAT(5,NUM_WAV_AER) =  &
     &  Reshape  &
     &    ((/ 4.03946d-01, 4.37936d+00, 1.43126d+01, 8.37816d+02,  &
     &       7.54765d-04 /),  &
     &     (/ 5,NUM_WAV_AER /))


!     ------------------------------------------------------------------
!     ehc_ss1 : extinction humidity coefficients used for small sea salt
!     ------------------------------------------------------------------

      real*8, parameter :: EHC_SS1(5,NUM_WAV_AER) =  &
     &  Reshape  &
     &    ((/ 6.1654d-01, 3.6957d+00, 9.1727d+00, 3.1506d+02,  &
     &        4.9787d-03 /),  &
     &     (/ 5,NUM_WAV_AER /))


!     ------------------------------------------------------------------
!     ehc_ss2 : extinction humidity coefficients used for large sea salt
!     ------------------------------------------------------------------

      real*8, parameter :: EHC_SS2(5,NUM_WAV_AER) =  &
     &  Reshape  &
     &    ((/ 6.1743d-01, 3.6212d+00, 9.1653d+00, 3.2070d+02,  &
     &        5.4119d-03 /),  &
     &     (/ 5,NUM_WAV_AER /))

