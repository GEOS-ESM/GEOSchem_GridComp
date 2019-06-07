!=======================================================================
!
! $Id$
!
! CODE DEVELOPER
!   Peter Connell, LLNL
!   connell2@llnl.gov
!
! FILE
!   setkin_group.h
!
! DESCRIPTION
!   This include file contains information about grouping species for
!   transport purposes.
!
!  Chemistry input file:    08/2018
!  Reaction dictionary:     GMI_Combo_rxns_119species_SO2_JPL15_OCS.db
!  Setkin files generated:  Tue Sep 18 18:25:41 2018
!
!=======================================================================

      integer, parameter :: NUMGRP      =  2

      integer, parameter :: MAXGRP_ELEM = 10

      character*16, save :: cgrp_nam(NUMGRP)

      data cgrp_nam(1) /"Bry"/
      data cgrp_nam(2) /"Cly"/

      integer :: mge1, ng1

      integer, save :: sgrp_elem_map(MAXGRP_ELEM, NUMGRP)

      data  ((sgrp_elem_map(mge1,ng1), mge1=1,MAXGRP_ELEM), &
     &                                 ng1 =1,NUMGRP) / &
     & 26,27,28,29,30,31, 0, 0, 0, 0, &
     & 32,33,34,35,36,37,38,39, 0, 0 /

      real*8,  save :: sgrp_fac(MAXGRP_ELEM, NUMGRP)

      data  ((sgrp_fac(mge1,ng1), mge1=1,MAXGRP_ELEM), &
     &                            ng1 =1,NUMGRP) / &
     &  1.0D+00,  1.0D+00,  1.0D+00,  1.0D+00,  1.0D+00,  &
     &  1.0D+00,  1.0D+00,  1.0D+00,  1.0D+00,  1.0D+00,  &
     &  1.0D+00,  2.0D+00,  1.0D+00,  2.0D+00,  1.0D+00,  &
     &  1.0D+00,  1.0D+00,  1.0D+00,  1.0D+00,  1.0D+00 /

!      -----------------------------------------------------------------
!      qq2 : grouped species sum
!      -----------------------------------------------------------------

      real*8,  save :: max_bry_adjust
      real*8,  save :: max_cly_adjust
      real*8,  save :: max_nox_adjust

      real*8  :: qqgrp        (i1:i2 ,ju1:j2 ,k1:k2 ,NUMGRP)

      real*8  :: group_adjust (i1:i2 ,ju1:j2 ,k1:k2)
      real*8  :: group_factor (i1:i2 ,ju1:j2 ,k1:k2)
      real*8  :: qq2          (i1:i2 ,ju1:j2 ,k1:k2)

!                                  --^--
