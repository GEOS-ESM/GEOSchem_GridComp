!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: MAM3_DataMod - basic MAM3 parameters and types
!
! !INTERFACE:
!
   MODULE MAM3_DataMod
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
!  {\tt MAM3\_DataMod} provides a collection of parameters and types 
!  used in the MAM3 code.
!
!
! !REVISION HISTORY:
!
!  14Sep2011  A. Darmenov  Initial version
!
!EOP
!-------------------------------------------------------------------------

! Number of modes
! ------------------
  integer, public, parameter :: MAM3_MODES = 3


! Mode IDs
! ------------------
  integer, public, parameter :: MAM3_AITKEN_MODE_ID       = 1
  integer, public, parameter :: MAM3_ACCUMULATION_MODE_ID = 2
  integer, public, parameter :: MAM3_COARSE_MODE_ID       = 3

  integer, public, parameter :: MAM3_AIT_ID = MAM3_AITKEN_MODE_ID
  integer, public, parameter :: MAM3_ACC_ID = MAM3_ACCUMULATION_MODE_ID
  integer, public, parameter :: MAM3_COR_ID = MAM3_COARSE_MODE_ID


! Mode Names
! ------------------
  character(len=*), public, parameter :: MAM3_AITKEN_MODE_NAME       = 'AIT'
  character(len=*), public, parameter :: MAM3_ACCUMULATION_MODE_NAME = 'ACC'
  character(len=*), public, parameter :: MAM3_COARSE_MODE_NAME       = 'COR'
 
  character(len=*), public, parameter :: MAM3_AIT_NAME = MAM3_AITKEN_MODE_NAME
  character(len=*), public, parameter :: MAM3_ACC_NAME = MAM3_ACCUMULATION_MODE_NAME
  character(len=*), public, parameter :: MAM3_COR_NAME = MAM3_COARSE_MODE_NAME


! Geometric standard deviation
! ------------------------------
  real, public, parameter    :: MAM3_AITKEN_MODE_SIGMA       = 1.6
  real, public, parameter    :: MAM3_ACCUMULATION_MODE_SIGMA = 1.8
  real, public, parameter    :: MAM3_COARSE_MODE_SIGMA       = 2.0

  real, public, parameter    :: MAM3_AIT_SIGMA = MAM3_AITKEN_MODE_SIGMA
  real, public, parameter    :: MAM3_ACC_SIGMA = MAM3_ACCUMULATION_MODE_SIGMA
  real, public, parameter    :: MAM3_COR_SIGMA = MAM3_COARSE_MODE_SIGMA


! Sea-salt emission cut-off size ranges in [m]
! ----------------------------------------------
  real, public, parameter    :: MAM3_AIT_SS_D_CUTOFF(2) = (/ 0.02, 0.08 /) * 1.0e-6
  real, public, parameter    :: MAM3_ACC_SS_D_CUTOFF(2) = (/ 0.08, 1.00 /) * 1.0e-6
  real, public, parameter    :: MAM3_COR_SS_D_CUTOFF(2) = (/ 1.00, 10.0 /) * 1.0e-6

! Dust emission cut-off size ranges in units [m]
! ----------------------------------------------
  real, public, parameter    :: MAM3_ACC_DU_D_CUTOFF(2) = (/ 0.10, 1.00 /) * 1.0e-6
  real, public, parameter    :: MAM3_COR_DU_D_CUTOFF(2) = (/ 1.00, 10.0 /) * 1.0e-6



! ...conveniently packed, use the MAM3_MODE_ID array to inquire data
! -------------------------------------------------------------------
  integer, public, parameter, dimension(MAM3_MODES) :: &
      MAM3_MODE_ID = (/ MAM3_AIT_ID, MAM3_ACC_ID, MAM3_COR_ID /)

  character(len=*), public, parameter, dimension(MAM3_MODES) :: &
      MAM3_MODE_NAME = (/ MAM3_AIT_NAME, MAM3_ACC_NAME, MAM3_COR_NAME /)

  real, public, parameter, dimension(MAM3_MODES) :: &
      MAM3_MODE_SIGMA = (/ MAM3_AIT_SIGMA, MAM3_ACC_SIGMA, MAM3_COR_SIGMA /)



  integer, public, parameter, dimension(3) ::                       &
                         MAM3_SS_EMISSION_MODE_ID = (/ MAM3_AIT_ID, &
                                                       MAM3_ACC_ID, &
                                                       MAM3_COR_ID /)

  real, public, parameter, dimension(3) ::                          &
        MAM3_SS_EMISSION_D_CUTOFF_LOW = (/ MAM3_AIT_SS_D_CUTOFF(1), &
                                           MAM3_ACC_SS_D_CUTOFF(1), &
                                           MAM3_COR_SS_D_CUTOFF(1) /)

  real, public, parameter, dimension(3) ::                          &
        MAM3_SS_EMISSION_D_CUTOFF_UP  = (/ MAM3_AIT_SS_D_CUTOFF(2), &
                                           MAM3_ACC_SS_D_CUTOFF(2), &
                                           MAM3_COR_SS_D_CUTOFF(2) /)  

  integer, public, parameter, dimension(2) ::                       &
                         MAM3_DU_EMISSION_MODE_ID = (/ MAM3_ACC_ID, &
                                                       MAM3_COR_ID /)

  real, public, parameter, dimension(2) ::                          &
        MAM3_DU_EMISSION_D_CUTOFF_LOW = (/ MAM3_ACC_DU_D_CUTOFF(1), &
                                           MAM3_COR_DU_D_CUTOFF(1) /)

  real, public, parameter, dimension(2) ::                          &
        MAM3_DU_EMISSION_D_CUTOFF_UP  = (/ MAM3_ACC_DU_D_CUTOFF(2), &
                                           MAM3_COR_DU_D_CUTOFF(2) /)

 END MODULE MAM3_DataMod
