
!=============================================================================
!
! $Id$
!
! CODE DEVELOPER
!   John Tannahill, LLNL
!   jrt@llnl.gov
!
! FILE
!   gmi_subdomains.h
!
! DESCRIPTION
!   This include file contains domain decomposition and communication
!   related declarations.
!
!=============================================================================


!     ----------
!     Constants.
!     ----------

      integer, parameter ::  &
     &  SIZSUB = 100  ! common block size constant


!     ----------------
!     Integer commons.
!     ----------------

!     ---------------------------------------------------------------------
!     Scalars that define the decomposition of a particular ACTM subdomain.
!     MUST BE COMMUNICATED TO SLAVE PROCESSORS.
!     ---------------------------------------------------------------------

      integer ::  &
     &  actm_package_number,  &
     &  actm_subdomain

      integer ::  &
     &  actm_east,  &
     &  actm_west,  &
     &  actm_north,  &
     &  actm_south

      integer ::  &
     &  eb,  &
     &  wb,  &
     &  nb,  &
     &  sb

      integer ::  &
     &  most_eastern_domain,  &
     &  most_western_domain

      integer :: gmi_sd_idum1
      integer :: gmi_sd_i1_array(SIZSUB)

!     ====================
      common  / gmisd_i1 /  &
!     ====================
     &  gmi_sd_idum1,  &
     &  actm_package_number, actm_subdomain,  &
     &  actm_east, actm_west, actm_north, actm_south,  &
     &  eb, wb, nb, sb,  &
     &  most_eastern_domain, most_western_domain

!     ===========
      equivalence (gmi_sd_i1_array, gmi_sd_idum1)
!     ===========


!     ----------------------------------------------------------------------
!     Scalars used for visiting groups of processors in an efficient manner.
!     MAINTAINED BY SLAVE PROCESSORS.
!     ----------------------------------------------------------------------

!     --------------
!     Communicators.
!     --------------

      integer ::  &
     &  commu_npole,  commu_spole,  &
     &  commu_slaves, commu_world

!     ====================
      common  / gmisd_i2 /  &
!     ====================
     &  commu_npole,  commu_spole,  &
     &  commu_slaves, commu_world

