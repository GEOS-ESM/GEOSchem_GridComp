
!=============================================================================
!
! $Id$
!
! CODE DEVELOPER
!   John Tannahill, LLNL
!   jrt@llnl.gov
!
! FILE
!   gmi_loss_data_18.h
!
! DESCRIPTION
!   This include file contains the tparm variable which can be used to
!   calculate loss frequencies.
!
!=============================================================================


      integer, parameter ::  &
     &  JDIM = 18,  &
     &  KDIM = 20,  &
     &  MDIM = 12


      real*8  :: tparm(KDIM, JDIM, MDIM, 1)

      common  / lossd18 /  &
     &          tparm

