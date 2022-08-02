module cam_logfile

!----------------------------------------------------------------------- 
! 
! Purpose: This module is responsible for managing the logical unit 
!          of CAM's output log
! 
! Author: mvr, Sep 2007
!         asd, Dec, 2014 (port to GEOS5)
! 
!-----------------------------------------------------------------------

#ifndef GEOS5_PORT
   implicit none
   integer :: iulog = 6
#else
   use, intrinsic :: iso_fortran_env, only: iulog => output_unit
   implicit none
#endif

end module cam_logfile
