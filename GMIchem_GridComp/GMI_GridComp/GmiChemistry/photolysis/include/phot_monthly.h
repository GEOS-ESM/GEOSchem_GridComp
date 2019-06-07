!=============================================================================
!
! CODE DEVELOPER
!   Huisheng bian, GSFC
!   bian@code916.gsfc.nasa.gov
!
! FILE
!   phot_monthly.h
!
! DESCRIPTION
!   This include file contains the declarations of the arrays for the
!   monthly photolysis rate.
!
! HISTORY
!
!=============================================================================


!     ------------------------------------------------------------------
!     qjmon   : monthly photolysis rate               (phot_opt=2)
!
!     ------------------------------------------------------------------

      real*8, pointer :: qjmon   (:,:,:,:,:)

!     =====================
      common  / gmipr_mon /  &
!     =====================
     &  qjmon

