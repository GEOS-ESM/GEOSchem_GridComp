
!=============================================================================
!
! $Id$
!
! CODE DEVELOPER
!   John Tannahill, LLNL
!   jrt@llnl.gov
!
! FILE
!   gmi_forc_bc.h
!
! DESCRIPTION
!   This include file contains forcing data dimensions as well as an array
!   to hold the forcing data itself.
!
!=============================================================================


      integer, parameter ::  &
     &  FBC_LATDIM = 19,    & ! number of latitude bins (-90,-80,...,89,90)
     &  FBC_MONDIM = 12   ! number of months


!     -----------------------------------------------------
!     forc_bc_data : forcing boundary condition data (ppmv)
!     -----------------------------------------------------

!      real*8, pointer :: forc_bc_data(:,:,:,:)
!
!!     =====================
!      common  / gmifbc_r1 /  &
!!     =====================
!     &  forc_bc_data

