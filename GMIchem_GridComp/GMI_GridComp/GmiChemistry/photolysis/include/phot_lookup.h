
!=============================================================================
!
! $Id$
!
! CODE DEVELOPER
!   Dan Bergmann, LLNL
!   dbergmann@llnl.gov
!
! FILE
!   phot_lookup.h
!
! DESCRIPTION
!   This include file contains the declarations of the photolysis reaction
!   numbers.
!
!=============================================================================


!     ----------------------------------------------------
!     All of the declarations below are "read from table."
!     ----------------------------------------------------

      integer ::  &
     &  num_qj_acet,        & ! photolysis reaction number for ACET
     &  num_qj_ch2o,        & ! photolysis reaction number for CH2O
     &  num_qj_hacn,        & ! photolysis reaction number for HACN
     &  num_qj_mgly,        & ! photolysis reaction number for MGLY
     &  num_qj_no,          & ! photolysis reaction number for NO
     &  num_qj_o2,          & ! photolysis reaction number for O2
     &  n_qj_O3_2OH           ! photolysis reaction number for O3 + hv = 2OH


!     ====================
      common  / gmipl_i1 /  &
!     ====================
     &  num_qj_acet, num_qj_ch2o, num_qj_hacn,     num_qj_mgly,  &
     &  num_qj_no,   num_qj_o2,   n_qj_O3_2OH

