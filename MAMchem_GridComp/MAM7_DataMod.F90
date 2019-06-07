!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 710.1     !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: MAM7_DataMod - basic MAM7 parameters and types
!
! !INTERFACE:
!
   MODULE MAM7_DataMod
!
! !USES:
!


   IMPLICIT NONE
   PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:

!
! !DESCRIPTION: 
!
!  {\tt MAM7\_DataMod} provides a collection of parameters and 
!  types used in the MAM code.
!
!
! !REVISION HISTORY:
!
!  14Sep2011  A. Darmenov  Initial version
!
!EOP
!-------------------------------------------------------------------------

! Number of aerosol components
! ----------------------------
  integer, public, parameter :: MAM7_AEROSOL_COMPONENTS = 7


! Number of modes
! ---------------
  integer, public, parameter :: MAM7_MODES = 7


! Mode IDs
! ------------------ 
  integer, public, parameter :: MAM7_AITKEN_MODE_ID         = 1
  integer, public, parameter :: MAM7_ACCUMULATION_MODE_ID   = 2
  integer, public, parameter :: MAM7_PRIMARY_CARBON_MODE_ID = 3
  integer, public, parameter :: MAM7_FINE_SEASALT_MODE_ID   = 4
  integer, public, parameter :: MAM7_FINE_DUST_MODE_ID      = 5
  integer, public, parameter :: MAM7_COARSE_SEASALT_MODE_ID = 6
  integer, public, parameter :: MAM7_COARSE_DUST_MODE_ID    = 7

! Mode Names
! ------------------
  character(len=*), public, parameter :: MAM7_AITKEN_MODE_NAME         = 'AIT'
  character(len=*), public, parameter :: MAM7_ACCUMULATION_MODE_NAME   = 'ACC'
  character(len=*), public, parameter :: MAM7_PRIMARY_CARBON_MODE_NAME = 'PCM'
  character(len=*), public, parameter :: MAM7_FINE_SEASALT_MODE_NAME   = 'FSS'
  character(len=*), public, parameter :: MAM7_FINE_DUST_MODE_NAME      = 'FDU'
  character(len=*), public, parameter :: MAM7_COARSE_SEASALT_MODE_NAME = 'CSS'
  character(len=*), public, parameter :: MAM7_COARSE_DUST_MODE_NAME    = 'CDU'
 
! Number of species
! -----------------
  integer, public, parameter    :: MAM7_AITKEN_MODE_SPECIES         = 4
  integer, public, parameter    :: MAM7_ACCUMULATION_MODE_SPECIES   = 6
  integer, public, parameter    :: MAM7_PRIMARY_CARBON_MODE_SPECIES = 2
  integer, public, parameter    :: MAM7_FINE_SEASALT_MODE_SPECIES   = 3
  integer, public, parameter    :: MAM7_FINE_DUST_MODE_SPECIES      = 3
  integer, public, parameter    :: MAM7_COARSE_SEASALT_MODE_SPECIES = 3
  integer, public, parameter    :: MAM7_COARSE_DUST_MODE_SPECIES    = 3

! Geometric standard deviation
! ------------------------------
  real, public, parameter    :: MAM7_AITKEN_MODE_SIGMA         = 1.6
  real, public, parameter    :: MAM7_ACCUMULATION_MODE_SIGMA   = 1.8
  real, public, parameter    :: MAM7_PRIMARY_CARBON_MODE_SIGMA = 1.6
  real, public, parameter    :: MAM7_FINE_SEASALT_MODE_SIGMA   = 2.0
  real, public, parameter    :: MAM7_FINE_DUST_MODE_SIGMA      = 1.8
  real, public, parameter    :: MAM7_COARSE_SEASALT_MODE_SIGMA = 2.0
  real, public, parameter    :: MAM7_COARSE_DUST_MODE_SIGMA    = 1.8

! Default Size - Median geometric diameter of number size distribution), 'm'
! ------------
  real, public, parameter    :: MAM7_AITKEN_MODE_SIZE         = 0.0260e-6
  real, public, parameter    :: MAM7_ACCUMULATION_MODE_SIZE   = 0.1100e-6
  real, public, parameter    :: MAM7_PRIMARY_CARBON_MODE_SIZE = 0.0500e-6
  real, public, parameter    :: MAM7_FINE_SEASALT_MODE_SIZE   = 0.2000e-6
  real, public, parameter    :: MAM7_FINE_DUST_MODE_SIZE      = 0.1000e-6
  real, public, parameter    :: MAM7_COARSE_SEASALT_MODE_SIZE = 2.0000e-6
  real, public, parameter    :: MAM7_COARSE_DUST_MODE_SIZE    = 1.0000e-6

! Minimum Size - lower limit of the number size distribution), 'm'
! ------------
  real, public, parameter    :: MAM7_AITKEN_MODE_SIZE_MIN         = 0.0087e-6
  real, public, parameter    :: MAM7_ACCUMULATION_MODE_SIZE_MIN   = 0.0535e-6
  real, public, parameter    :: MAM7_PRIMARY_CARBON_MODE_SIZE_MIN = 0.0100e-6
  real, public, parameter    :: MAM7_FINE_SEASALT_MODE_SIZE_MIN   = 0.0500e-6
  real, public, parameter    :: MAM7_FINE_DUST_MODE_SIZE_MIN      = 0.0500e-6
  real, public, parameter    :: MAM7_COARSE_SEASALT_MODE_SIZE_MIN = 1.0000e-6
  real, public, parameter    :: MAM7_COARSE_DUST_MODE_SIZE_MIN    = 0.5000e-6

! Maximum Size - upper limit of the number size distribution), 'm'
! ------------
  real, public, parameter    :: MAM7_AITKEN_MODE_SIZE_MAX         = 0.0520e-6
  real, public, parameter    :: MAM7_ACCUMULATION_MODE_SIZE_MAX   = 0.4400e-6
  real, public, parameter    :: MAM7_PRIMARY_CARBON_MODE_SIZE_MAX = 0.1000e-6
  real, public, parameter    :: MAM7_FINE_SEASALT_MODE_SIZE_MAX   = 1.0000e-6
  real, public, parameter    :: MAM7_FINE_DUST_MODE_SIZE_MAX      = 0.5000e-6
  real, public, parameter    :: MAM7_COARSE_SEASALT_MODE_SIZE_MAX = 4.0000e-6
  real, public, parameter    :: MAM7_COARSE_DUST_MODE_SIZE_MAX    = 2.0000e-6

! Crystallization RH points
! ------------------------
  real, public, parameter    :: MAM7_AITKEN_MODE_RH_CRYSTALLIZATION         = 0.350
  real, public, parameter    :: MAM7_ACCUMULATION_MODE_RH_CRYSTALLIZATION   = 0.350
  real, public, parameter    :: MAM7_PRIMARY_CARBON_MODE_RH_CRYSTALLIZATION = 0.350
  real, public, parameter    :: MAM7_FINE_SEASALT_MODE_RH_CRYSTALLIZATION   = 0.350
  real, public, parameter    :: MAM7_FINE_DUST_MODE_RH_CRYSTALLIZATION      = 0.350
  real, public, parameter    :: MAM7_COARSE_SEASALT_MODE_RH_CRYSTALLIZATION = 0.350
  real, public, parameter    :: MAM7_COARSE_DUST_MODE_RH_CRYSTALLIZATION    = 0.350

! Deliquescence RH points
! ------------------------
  real, public, parameter    :: MAM7_AITKEN_MODE_RH_DELIQUESCENCE         = 0.800
  real, public, parameter    :: MAM7_ACCUMULATION_MODE_RH_DELIQUESCENCE   = 0.800
  real, public, parameter    :: MAM7_PRIMARY_CARBON_MODE_RH_DELIQUESCENCE = 0.800
  real, public, parameter    :: MAM7_FINE_SEASALT_MODE_RH_DELIQUESCENCE   = 0.800
  real, public, parameter    :: MAM7_FINE_DUST_MODE_RH_DELIQUESCENCE      = 0.800
  real, public, parameter    :: MAM7_COARSE_SEASALT_MODE_RH_DELIQUESCENCE = 0.800
  real, public, parameter    :: MAM7_COARSE_DUST_MODE_RH_DELIQUESCENCE    = 0.800

! Sea-salt emission cut-off size ranges in [m]
! ----------------------------------------------
 real, public, parameter    :: MAM7_AIT_SS_D_CUTOFF(2) = (/ 0.02, 0.08 /) * 1.0e-6
 real, public, parameter    :: MAM7_ACC_SS_D_CUTOFF(2) = (/ 0.08, 0.30 /) * 1.0e-6
 real, public, parameter    :: MAM7_FSS_SS_D_CUTOFF(2) = (/ 0.30, 1.00 /) * 1.0e-6
 real, public, parameter    :: MAM7_CSS_SS_D_CUTOFF(2) = (/ 1.00, 10.0 /) * 1.0e-6


! Dust emission cut-off size ranges in units [m]
! ----------------------------------------------
 real, public, parameter    :: MAM7_FDU_DU_D_CUTOFF(2) = (/ 0.10, 2.00 /) * 1.0e-6
 real, public, parameter    :: MAM7_CDU_DU_D_CUTOFF(2) = (/ 2.00, 10.0 /) * 1.0e-6



! ...conveniently packed, use the MAM7_MODE_ID array to inquire data
! -------------------------------------------------------------------

  integer, public, parameter, dimension(MAM7_MODES) ::  &
      MAM7_MODE_ID =    (/ MAM7_AITKEN_MODE_ID,         &
                           MAM7_ACCUMULATION_MODE_ID,   &
                           MAM7_PRIMARY_CARBON_MODE_ID, &
                           MAM7_FINE_SEASALT_MODE_ID,   &
                           MAM7_FINE_DUST_MODE_ID,      &
                           MAM7_COARSE_SEASALT_MODE_ID, &
                           MAM7_COARSE_DUST_MODE_ID /)
  
  character(len=*), public, parameter, dimension(MAM7_MODES) :: &
      MAM7_MODE_NAME  = (/ MAM7_AITKEN_MODE_NAME, &
                           MAM7_ACCUMULATION_MODE_NAME, &
                           MAM7_PRIMARY_CARBON_MODE_NAME, &
                           MAM7_FINE_SEASALT_MODE_NAME, &
                           MAM7_FINE_DUST_MODE_NAME, &
                           MAM7_COARSE_SEASALT_MODE_NAME, &
                           MAM7_COARSE_DUST_MODE_NAME /)

  real, public, parameter, dimension(MAM7_MODES) :: &
      MAM7_MODE_SIGMA = (/ MAM7_AITKEN_MODE_SIGMA, &
                           MAM7_ACCUMULATION_MODE_SIGMA, &
                           MAM7_PRIMARY_CARBON_MODE_SIGMA, &
                           MAM7_FINE_SEASALT_MODE_SIGMA, &
                           MAM7_FINE_DUST_MODE_SIGMA, &
                           MAM7_COARSE_SEASALT_MODE_SIGMA, &
                           MAM7_COARSE_DUST_MODE_SIGMA /)


  integer, public, parameter, dimension(4) :: &
      MAM7_SS_EMISSION_MODE_ID =      (/ MAM7_AITKEN_MODE_ID, &
                                         MAM7_ACCUMULATION_MODE_ID, &
                                         MAM7_FINE_SEASALT_MODE_ID, &
                                         MAM7_COARSE_SEASALT_MODE_ID /)

  real, public, parameter, dimension(4) :: &
      MAM7_SS_EMISSION_D_CUTOFF_LOW = (/ MAM7_AIT_SS_D_CUTOFF(1), &
                                         MAM7_ACC_SS_D_CUTOFF(1), &
                                         MAM7_FSS_SS_D_CUTOFF(1), &
                                         MAM7_CSS_SS_D_CUTOFF(1) /)

  real, public, parameter, dimension(4) :: &
      MAM7_SS_EMISSION_D_CUTOFF_UP  = (/ MAM7_AIT_SS_D_CUTOFF(2), &
                                         MAM7_ACC_SS_D_CUTOFF(2), &
                                         MAM7_FSS_SS_D_CUTOFF(2), &
                                         MAM7_CSS_SS_D_CUTOFF(2) /)

  integer, public, parameter, dimension(2) :: & 
      MAM7_DU_EMISSION_MODE_ID      = (/ MAM7_FINE_DUST_MODE_ID, &
                                         MAM7_COARSE_DUST_MODE_ID /)

  real, public, parameter, dimension(2) :: &
      MAM7_DU_EMISSION_D_CUTOFF_LOW = (/ MAM7_FDU_DU_D_CUTOFF(1), &
                                         MAM7_CDU_DU_D_CUTOFF(1) /)

  real, public, parameter, dimension(2) :: &
      MAM7_DU_EMISSION_D_CUTOFF_UP  = (/ MAM7_FDU_DU_D_CUTOFF(2), &
                                         MAM7_CDU_DU_D_CUTOFF(2) /) 
                                                                              

  END MODULE MAM7_DataMod
