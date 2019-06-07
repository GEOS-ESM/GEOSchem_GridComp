!------------------------------------------------------------------------------
!BOP
!
! !MODULE: TR_GravSetMod
!
      module TR_GravSetMod
!
! !USES:
!     use GmiSpcConcentrationMethod_mod, only : t_SpeciesConcentration
   USE ESMF
   USE MAPL_Mod

!
      implicit none
!
      private
!
! !PUBLIC MEMBER FUNCTIONS:
      public  :: TR_GMI_GravitationalSettling
!
! !AUTHOR:
!  Michael Manyin, SSAI/GSFC, michael.manyin@nasa.gov
!
!EOP
!------------------------------------------------------------------------------
      CONTAINS
!------------------------------------------------------------------------------
!BOP
!     
! !ROUTINE: TR_GMI_GravitationalSettling
!
! !INTERFACE:
!
      subroutine TR_GMI_GravitationalSettling ( species,  &
     &   AER_DENSITY, RADIUS_EFF, C1, C2, C3, C4,  tdt,  &
     &   grid_height, humidity, mass, press3e, kel,      &
     &   diffaer, s_radius, s_velocity,                  &
     &   i1, i2, j1, j2, km )

! A variant of code in GmiGravitationalSettling_mod.F90

      implicit none

#     include "gmi_phys_constants.h"

#if 0

#ifdef MICRO_AEROSOL
#     include "../GMIchem_GridComp/GMI_GridComp/GmiShared/GmiInclude/gmi_micro_aerosol.h"
#elif GOCARTaerosol
#     include "../GMIchem_GridComp/GMI_GridComp/GmiShared/GmiInclude/gocart_aerosol.h"
#else
#     include "../GMIchem_GridComp/GMI_GridComp/GmiShared/GmiInclude/gmi_aerosol.h"
#endif

#endif

!c?   Tight coupling to setkin?
! #     include "../GMIchem_GridComp/GMI_GridComp/GmiShared/GmiInclude/setkin_par.h"
! #     include "../GMIchem_GridComp/GMI_GridComp/GmiShared/GmiInclude/setkin_depos.h"

#if 0
#ifdef MICRO_AEROSOL
#     include "../GMIchem_GridComp/GMI_GridComp/GmiShared/GmiInclude/umaerosol.h"
#endif
#endif
!
! !INPUT PARAMETERS:
      real*4  , intent(in) :: AER_DENSITY, RADIUS_EFF, C1, C2, C3, C4
      real*4  , intent(in) :: tdt
      integer , intent(in) :: i1, i2, j1, j2, km

      ! NOTE  This was converted from the GMI vertical system (bottom-up)
      !       to the GEOS5 vertical system (top-down)

      real*4  , intent(in) :: grid_height(i1:i2, j1:j2, 1:km) ! height of each grid box (m)
      real*4  , intent(in) :: humidity   (i1:i2, j1:j2, 1:km) ! specific humidity
      real*4  , intent(in) :: mass       (i1:i2, j1:j2, 1:km) ! total mass of the atmosphere within each grid box (kg)
      real*4  , intent(in) :: press3e    (i1:i2, j1:j2, 0:km) ! atmospheric pressure at the edge of each grid box (mb)
      real*4  , intent(in) :: kel        (i1:i2, j1:j2, 1:km) ! temperature (degK)
!
! !OUTPUT PARAMETERS:
      real*8 , intent(out) :: diffaer    (i1:i2, j1:j2)       ! aerosol diffusivity at bottom layer
      real*8 , intent(out) :: s_radius   (i1:i2, j1:j2)       ! aerosol radius at bottom layer (m)
      real*8 , intent(out) :: s_velocity (i1:i2, j1:j2)       ! aerosol settling velocity at bottom layer (m/s)
!
! !INOUT/OUTPUT PARAMETERS:
      ! species concentration, known at zone centers (mixing ratio)
      REAL*4,    POINTER, DIMENSION(:,:,:), INTENT(INOUT)    :: species


! !DESCRIPTION:
!   This routine updates the gravitational settling of aerosols.
!
! !LOCAL VARIABLES:
      real*4  :: aer_den_swel(i1:i2, j1:j2, 1:km-1)   ! top down; values at interfaces
      real*4  :: c3t         (i1:i2, j1:j2, 1:km-1)
      real*4  :: del_grid_box(i1:i2, j1:j2, 1:km-1)
      real*4  :: kele        (i1:i2, j1:j2, 1:km-1)
      real*4  :: mfp         (i1:i2, j1:j2, 1:km-1)
      real*8  :: radius      (i1:i2, j1:j2, 1:km-1)
      real*4  :: relhume     (i1:i2, j1:j2, 1:km-1)
      real*8  :: slip        (i1:i2, j1:j2, 1:km-1)
      real*8  :: velocity    (i1:i2, j1:j2, 1:km-1)
      real*8  :: viscosity   (i1:i2, j1:j2, 1:km-1)

!EOP
!------------------------------------------------------------------------------
!BOC


! NOTE: lead had aero index 1 - sulfate
!       beryllium had aero index 15 - stratospheric sulfate (monodispersed)


!     ------------------------------------------------------------
!     Calculate all aerodynamic terms which are not dependent on 
!     physical characteristics of each aerosol type.
!     ------------------------------------------------------------

!     ----------------------------------------
!     Calculate temperature at grid interface.
!     ----------------------------------------

      kele(i1:i2,j1:j2,1:km-1) = (kel(i1:i2,j1:j2,1:km-1) +  &
     &                            kel(i1:i2,j1:j2,2:km  )) * 0.5e0

!     -----------------------------------------------------------------
!     First calculate relative humidity from Seinfeld (1986) p. 181.
!     The first  relhume is the temperature dependent parameter a.
!     The second relhume is the saturation vapor pressure of water.
!     The third  relhume is the actual relative humidity as a fraction.
!     Then make sure relhume is between 0 and 0.95 because swelling
!     parameterization goes too large at very high relative humidity.
!     -----------------------------------------------------------------

      relhume(:,:,:) = 1.0e0 - (373.15e0 / kele(:,:,:))

      relhume(:,:,:) =  &
     &  1013.25e0 * Exp (13.3185e0 * relhume(:,:,:)    -  &
     &                    1.9760e0 * relhume(:,:,:)**2 -  &
     &                    0.6445e0 * relhume(:,:,:)**3 -  &
     &                    0.1299e0 * relhume(:,:,:)**4)

      relhume(:,:,:) =  &
     &  (humidity(:,:,1:km-1) + humidity(:,:,2:km)) *  &
     &  0.5d0 * MWTAIR / 18.0d0 /  &
     &  GPKG * press3e(i1:i2,j1:j2,1:km-1) / relhume(:,:,:)

      relhume(:,:,:) = Max (Min (relhume(:,:,:), 0.95d0), 1.0d-30)

!     -----------------------------------------
!     Viscosity (kg/m/s) (fit to CRC handbook).
!     -----------------------------------------

      viscosity(:,:,:) = 1.8d-5 * (kele(:,:,:) / 298.0d0)**0.85d0

!     ----------------------------------------------------
!     Mean free path (m) (Seinfeld + Pandis 1998, p. 455).
!     ----------------------------------------------------

      mfp(:,:,:) =  &
     &  2.0d0 * viscosity(:,:,:) /  &
     &  (PASPMB * press3e(i1:i2,j1:j2,1:km-1) *  &
     &   Sqrt (8.0d0 * MWTAIR /  &
     &         (GMI_PI * GAS_CONST_J * 1.0d3 * kele(:,:,:))))


!       -----------------------------------------------------------
!       Swelling of aerosol due to humidity
!       (Gong, Barrie, and Blanchet 1997 JGR 3805-3818 equation 3).
!       -----------------------------------------------------------


        if (C1 == 0.0d0) then

          radius(:,:,:) = RADIUS_EFF

          aer_den_swel(:,:,:) = AER_DENSITY

        else

          c3t(:,:,:) =  &
     &      C3 *  &
     &      (1.0d0 + 0.004d0 * (298.0d0 - kele(:,:,:)))

          radius(:,:,:) =  &
     &      (C1 *  &
     &       (RADIUS_EFF*CMPM)**C2 /  &
     &      (c3t(:,:,:) *  &
     &       (RADIUS_EFF*CMPM)**C4 -  &
     &      Log10 (relhume(:,:,:))) +  &
     &      (RADIUS_EFF*CMPM)**3)**(1.0d0/3.0d0)


          aer_den_swel(:,:,:) = ((AER_DENSITY - 1.0d3) *  &
     &      (RADIUS_EFF*CMPM)**3 + 1.0d3 *  &
     &      radius(:,:,:)**3) / radius(:,:,:)**3

          radius(:,:,:) = radius(:,:,:) / CMPM

        end if

!       --------------------------------------------------------
!       Slip correction factor (Seinfeld + Pandis 1998, p. 464).
!       --------------------------------------------------------

        slip(:,:,:) =  &
     &    1.0d0 + (mfp(:,:,:) / radius(:,:,:)) *  &
     &    (1.257d0 +  &
     &     (0.4d0 * Exp (-(1.1d0 * radius(:,:,:) / mfp(:,:,:)))))


!       ---------------------------------------------------
!       Settling velocity (Seinfeld + Pandis 1998, p. 466).
!       ---------------------------------------------------

        velocity(:,:,:) =  &
     &    4.0d0 * radius(:,:,:)**2 *  &
     &    aer_den_swel(:,:,:) * GMI_G * slip(:,:,:) /  &
     &    (18.0d0 * viscosity(:,:,:))

!       ------------------------------------------------------
!       Aerosol diffusivity (Seinfeld + Pandis 1998, p. 474).
!       and velocity at k = km-1
!       ------------------------------------------------------

        diffaer(:,:) = BOLTZMN_J * kele(:,:,km-1) * slip(:,:,km-1) /  &
     &    (6.0d0 * GMI_PI * viscosity(:,:,km-1) * radius(:,:,km-1))

        s_radius(:,:) = radius(:,:,km-1)

        s_velocity(:,:) = velocity(:,:,km-1)

!       --------------------------------------------------------
!       Apply the settling velocity making sure that it does not
!       exceed the Courant limit.
!       --------------------------------------------------------

        del_grid_box(:,:,:) =  &
     &    velocity(:,:,:) * tdt / grid_height(:,:,1:km-1)

        where (del_grid_box(:,:,:) > 1.0d0)  &
     &    del_grid_box(:,:,:) = 1.0d0

        del_grid_box(:,:,:) =  &
     &    del_grid_box(:,:,:) * species(:,:,1:km-1)

        species(:,:,2:km) =  &
     &  species(:,:,2:km) +  &
     &    (del_grid_box(:,:,1:km-1) *  &
     &     mass(:,:,1:km-1) / mass(:,:,2:km))

        species(:,:,1:km-1) =  &
     &  species(:,:,1:km-1) -  &
     &    del_grid_box(:,:,1:km-1)

      return

      end subroutine TR_GMI_GravitationalSettling
!EOC
!------------------------------------------------------------------------------
      end module TR_GravSetMod
