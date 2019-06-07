
!=============================================================================
!
! $Id$
!
! CODE DEVELOPER
!   Dan Bergmann, LLNL
!   dbergmann@llnl.gov
!
! FILE
!   phot_lookup_arrays.h
!
! DESCRIPTION
!   This include file contains the declarations of the arrays for the
!   photolysis lookup code.
!
! HISTORY
!   * September 23, 2004 - Jules Kouatchou
!     Introduced the variable "o1d_coef" needed for the combined
!     strat/trop chemical mechanism.
!
!=============================================================================


!     ------------------------------------------------------------------
!     o3_clim_prs   : ozone climatology pressure coordinate (phot_opt=5)
!
!     prs_phot      : pressure           values for lookup table
!     sza_phot      : solar zenith angle values for lookup table
!
!     col_o3        : column ozone       values for lookup table
!                     (which vary with pressure)
!
!     no_qj         : photolysis rates for NO
!     o2_qj         : photolysis rates for O2
!
!     cross_section : cross_sections for lookup table
!                     (as a function of temperature and wavelength)
!     ------------------------------------------------------------------

      real*8 :: o3_clim_prs(NUM_O3CLIM_PRS)

      real*8 :: prs_phot(NUMPRS)
      real*8 :: sza_phot(NUMSZA)

      real*8 :: col_o3(NUMO3, NUMPRS)

      real*8 :: no_qj(NUMSZA, NUMO3, NUMPRS)
      real*8 :: o2_qj(NUMSZA, NUMO3, NUMPRS)

      real*8 :: cross_section(NUMLAM, NUMTMP, MAX_NUMQJ_PL)

!     =====================
      common  / gmipla_r1 /  &
!     =====================
     &  o3_clim_prs,  &
     &  prs_phot, sza_phot,  &
     &  col_o3,  &
     &  no_qj, o2_qj,  &
     &  cross_section


!     -------------------------------------------------------
!     o3_clim    : ozone climatology (phot_opt=5)
!     rad_source : radiative source function for lookup table
!     o1d_coef   : table of adjustment factors for O3 + hv
!     -------------------------------------------------------

      real*8, pointer :: o3_clim   (:,:,:,:)
      real*8, pointer :: rad_source(:,:,:,:)
      real*8, pointer :: o1d_coef  (:,:,:,:,:)

!     =====================
      common  / gmipla_r2 /  &
!     =====================
     &  o3_clim, rad_source, o1d_coef

