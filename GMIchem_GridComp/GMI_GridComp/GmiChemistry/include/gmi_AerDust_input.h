!=============================================================================
!
! CODE DEVELOPER
!   Jules Kouatchou
!   Jules.Kouatchou@gsfc.nasa.gov
!
! FILE
!   gmi_AerDust_input.h
!
! DESCRIPTION
!   This include file contains declarations of input variables needed for
!   aerosol/dust optical dept calculations and diagnostics.
!   No variables here is associated with microphysical calculations.
!
! HISTORY
!
!=============================================================================

      ! DAERSL    : Mass density of hydrophobic aerosol (kg/m3)
      ! WAERSL    : Mass density of hydrophilic aerosol (kg/m3)
      ! DUST      : Mass density of Dust                (kg/m3)

      real*8, pointer :: DAERSL  (:,:,:,:)
      real*8, pointer :: WAERSL  (:,:,:,:)
      real*8, pointer :: DUST    (:,:,:,:)

!     =========================
      common  / gmi_MassDensity /  &
!     =========================
     &  DAERSL, WAERSL, DUST

      ! QAA_b     : Aerosol scattering phase functions
      ! RAA_b     : Effective radius associated with aerosol type

      real*8, pointer :: QAA_b(:,:)
      real*8, pointer :: RAA_b(:,:)

      COMMON /aer_vars/ QAA_b, RAA_b
