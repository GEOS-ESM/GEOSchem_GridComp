
!=============================================================================
!
! $Id$
!
! CODE DEVELOPER
!   John Tannahill, LLNL
!   jrt@llnl.gov
!
! FILE
!   gmi_emiss_harvard.h
!
! DESCRIPTION
!   This include file contains declarations for the Harvard biogneic and
!   soil emissions.
!
! HISTORY
!   - August 12, 2005 * Jules Kouatchou
!     Since NLANDHAR is no longer a parameter, we declare the following
!     variables as allocatable pointer arrays: index_soil, soil_fert,
!     soil_precip, and soil_pulse.
!
!=============================================================================


!     ----------------
!     Integer commons.
!     ----------------

      integer ::  &
     &  ncon_soil (NVEGTYPE)     ! Olson -> soil type
                                 !    1 => water/desert/ice
                                 !    2 => tropical rain forest
                                 !    3 => conifers
                                 !    4 => dry deciduous
                                 !    5 => other deciduous
                                 !    6 => woodland
                                 !    7 => grassland
                                 !    8 => agriculture (other than rice)
                                 !    9 => rice paddies
                                 !   10 => wetland/tundra

      integer, pointer :: index_soil(:, :)
!      integer ::
!     &  index_soil(2, NLANDHAR)  ! i,j of the grid

!     ====================
      common  / gmieh_i1 /  &
!     ====================
     &  ncon_soil,  &
     &  index_soil


      integer, pointer ::  &
     &  ireg (:,:)         ! number of land types in a grid square

      integer, pointer ::  &
     &  iland(:,:,:),        & ! land type id in grid square for ireg land types
     &  iuse (:,:,:)       ! fraction of grid box area occupied by land type
                           ! (mil^-1?)

!     ====================
      common  / gmieh_i2 /  &
!     ====================
     &  ireg, iland, iuse


!     -------------
!     Real commons.
!     -------------

      real*8  ::  &
     &  coeff_isop   (NPOLY),       & ! coefficients used for polynomial fit
     &  convert_isop (NVEGTYPE),    & ! isoprene    emissions by landtype
                                  ! (atomsC/cm^2leaf/s)
     &  convert_monot(NVEGTYPE)   ! monoterpene emissions by landtype
                                  ! (atomsC/cm^2leaf/s)
!      real*8  ::
!     &  soil_fert    (NLANDHAR)   ! fertilizers (ng N/m^2/s)
!
!      real*8  ::
!     &  soil_precip (2, NLANDHAR),        ! two months of observed precip.
!                                          ! (mm/day/box)
!     &  soil_pulse  (NPULSE+1, NLANDHAR)  ! tracking of wet/dry & three types
!                                          ! of pulsing (Y&L, 94)
      real*8, pointer :: soil_fert   (:)
      real*8, pointer :: soil_precip (:,:)
      real*8, pointer :: soil_pulse  (:,:)

!     ====================
      common  / gmieh_r1 /  &
!     ====================
     &  coeff_isop,  &
     &  convert_isop, convert_monot,  &
     &  soil_fert,  &
     &  soil_precip, soil_pulse


!     ----------------------------------------------------------------
!     base_isop  : baseline emissions for isoprene     (kgC/box/step?)
!     base_monot : baseline emissions for monoterpenes (kgC/box/step?)
!
!     xlai       : leaf area index of land type for month #1
!     xlai2      : leaf area index of land type for month #2
!     ----------------------------------------------------------------

      real*8,  pointer ::  &
     &  base_isop (:,:,:),  &
     &  base_monot(:,:,:),  &
     &  xlai      (:,:,:),  &
     &  xlai2     (:,:,:)

!     ====================
      common  / gmieh_r2 /  &
!     ====================
     &  base_isop,  base_monot,  &
     &  xlai, xlai2

