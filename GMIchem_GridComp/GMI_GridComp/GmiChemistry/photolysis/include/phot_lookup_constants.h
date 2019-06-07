
!=============================================================================
!
! $Id$
!
! CODE DEVELOPER
!   Dan Bergmann, LLNL
!   dbergmann@llnl.gov
!
! FILE
!   phot_lookup_constants.h
!
! DESCRIPTION
!   This include file contains the declarations of the constants for the
!   photolysis lookup code.
!
!=============================================================================


      integer, parameter ::  &
     &  MAX_NUMQJ_PL = 200    ! maximum number of qj (photolysis) rxns

      integer, parameter ::  &
     &  NUMLAM =  79,    & ! number of "lambdas" (wavelengths) for the radiative
                       ! source function and the photolysis cross-sections
     &  NUMO3  =  12,    & ! number of column ozone       coords. in lookup table
     &  NUMPRS =  35,    & ! number of pressure           coords. in lookup table
     &  NUMSZA =  20,    & ! number of solar zenith angle coords. in lookup table
     &  NUMTMP = 200   ! number of temperature        values  in lookup table

      integer, parameter ::  &
     &  NUM_O3CLIM_MON = 12,    & ! number of months    in the ozone climatology
     &  NUM_O3CLIM_PRS = 30   ! number of pressures in the ozone climatology


      real*8, parameter  :: ACET_CONST_1(47:58) =  &
     &  (/ 1.21363d+00, 1.46791d+00, 1.72995d+00, 2.78158d+00,  &
     &     4.88461d+00, 7.02881d+00, 9.17301d+00, 1.13172d+01,  &
     &     1.34614d+01, 1.34614d+01, 1.34614d+01, 1.34614d+01 /)

      real*8, parameter  :: ACET_CONST_2(47:58) =  &
     &  (/ 8.69881d-20, 1.70934d-19, 3.42893d-19, 5.63761d-19,  &
     &     7.02667d-19, 8.79584d-19, 1.10104d-18, 1.37826d-18,  &
     &     1.72528d-18, 1.72528d-18, 1.72528d-18, 1.72528d-18 /)


      real*8, parameter  :: CH2O_CONST_1(55:59) =  &
     &  (/ 0.9679353d+00, 0.9943720d+00, 0.9990122d+00, 0.9998249d+00,  &
     &     0.9999690d+00 /)

      real*8, parameter  :: CH2O_CONST_2(55:59) =  &
     &  (/ 0.2069d+00, 0.4583d+00, 0.8421d+00, 1.5000d+00, 2.8889d+00 /)


      real*8, parameter  :: MGLY_CONST_1(53:76) =  &
     &  (/ 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00,  &
     &     0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00, 0.000d+00,  &
     &     0.000d+00, 0.000d+00, 1.109d+00, 1.406d+00, 1.772d+00,  &
     &     2.221d+00, 2.767d+00, 3.430d+00, 4.230d+00, 5.190d+00,  &
     &     6.337d+00, 7.702d+00, 9.320d+00, 1.123d+01 /)

      real*8, parameter  :: MGLY_CONST_2(53:76) =  &
     &  (/ 7.485d-21, 1.128d-20, 1.679d-20, 2.469d-20, 3.592d-20,  &
     &     5.170d-20, 7.364d-20, 1.039d-19, 1.452d-19, 2.010d-19,  &
     &     2.760d-19, 3.757d-19, 5.074d-19, 6.800d-19, 9.046d-19,  &
     &     1.195d-18, 1.567d-18, 2.042d-18, 2.645d-18, 3.403d-18,  &
     &     4.354d-18, 5.538d-18, 7.006d-18, 8.815d-18 /)

