!=======================================================================
!
! $Id: setkin_group.h $
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
!  Input mechanism:        GeosCCM_Combo_2015mechanism.txt
!  Reaction dictionary:    GMI_reactions_JPL15.db
!  Setkin files generated: Tue Mar 30 16:58:02 2021
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
     &  7, 8, 9,10,46,63, 0, 0, 0, 0, &
     & 29,30,31,32,33,50,64,96, 0, 0 /

      real*8,  save :: sgrp_fac(MAXGRP_ELEM, NUMGRP)

      data  ((sgrp_fac(mge1,ng1), mge1=1,MAXGRP_ELEM), &
     &                            ng1 =1,NUMGRP) / &
     &  1.0D+00,  1.0D+00,  1.0D+00,  1.0D+00,  1.0D+00,  &
     &  1.0D+00,  1.0D+00,  1.0D+00,  1.0D+00,  1.0D+00,  &
     &  2.0D+00,  2.0D+00,  1.0D+00,  1.0D+00,  1.0D+00,  &
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
