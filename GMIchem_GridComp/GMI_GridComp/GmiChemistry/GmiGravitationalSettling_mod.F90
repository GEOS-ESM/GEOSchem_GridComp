!------------------------------------------------------------------------------
!BOP
!
! !MODULE: GmiGravitationalSettling_mod
!
      module GmiGravitationalSettling_mod
!
! !USES:
      use GmiArrayBundlePointer_mod    , only : t_GmiArrayBundle
      use GmiSpcConcentrationMethod_mod, only : t_SpeciesConcentration
      use GmiSpcConcentrationMethod_mod, only : Get_concentration
      use GmiSpcConcentrationMethod_mod, only : Set_concentration
      use GmiSpcConcentrationMethod_mod, only : isFixedConcentration
!
      implicit none
!
      private
!
! !PUBLIC MEMBER FUNCTIONS:
      public  :: updateGravitationalSettling
!
! !AUTHOR:
!  Dan Bergmann, LLNL, dbergmann@llnl.gov
!  Jules Kouatchou, GSFC, Jules.Kouatchou-1@nasa.gov
!
!EOP
!------------------------------------------------------------------------------
      CONTAINS
!------------------------------------------------------------------------------
!BOP
!     
! !ROUTINE: updateGravitationalSettling
!
! !INTERFACE:
!
      subroutine updateGravitationalSettling (SpeciesConcentration,            &
     &                 grid_height, humidity, mass, press3e, kel, diffaer,     &
     &                 s_radius, s_velocity, pr_diag, loc_proc, chem_opt, i1,  &
     &                 i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi, num_species,  &
     &                 tdt4)

      implicit none

#     include "gmi_phys_constants.h"

#ifdef MICRO_AEROSOL
#     include "gmi_micro_aerosol.h"
#elif GOCARTaerosol
#     include "gocart_aerosol.h"
#else
#     include "gmi_aerosol.h"
#endif

!c?   Tight coupling to setkin?
#     include "setkin_par.h"
#     include "setkin_depos.h"

#ifdef MICRO_AEROSOL
#     include "umaerosol.h"
#endif
!
! !INPUT PARAMETERS:
      logical          , intent(in) :: pr_diag
      integer          , intent(in) :: loc_proc
      integer          , intent(in) :: i1, i2, ju1, j2, k1, k2
      integer          , intent(in) :: ilo, ihi, julo, jhi
      integer          , intent(in) :: num_species
      integer          , intent(in) :: chem_opt
      ! height of each grid box (m)
      real*8           , intent(in) :: grid_height(i1:i2, ju1:j2, k1:k2)
      ! specific humidity
      real*8           , intent(in) :: humidity   (i1:i2, ju1:j2, k1:k2)
      ! total mass of the atmosphere within each grid box (kg)
      real*8           , intent(in) :: mass       (i1:i2, ju1:j2, k1:k2)
      ! atmospheric pressure at the edge of each grid box (mb)
      real*8           , intent(in) :: press3e    (ilo:ihi, julo:jhi, k1-1:k2)
      ! temperature (degK)
      real*8           , intent(in) :: kel        (ilo:ihi, julo:jhi, k1:k2)
      real             , intent(in) :: tdt4   ! timestep
!
! !OUTPUT PARAMETERS:
      ! aerosol diffusivity at bottom layer
      real*8 , intent(out) :: diffaer    (i1:i2, ju1:j2, num_species)
      ! aerosol radius at bottom layer (m)
      real*8 , intent(out) :: s_radius   (i1:i2, ju1:j2, num_species)
      ! aerosol settling velocity at bottom layer (m/s)
      real*8 , intent(out) :: s_velocity (i1:i2, ju1:j2, num_species)
!
! !INOUT/OUTPUT PARAMETERS:
      ! species concentration, known at zone centers (mixing ratio)
      type(t_SpeciesConcentration), intent(inOut) :: SpeciesConcentration

! !DESCRIPTION:
!   This routine updates the gravitational settling of aerosols.
!
!   These different aerosol types having settling velocities:
!  \begin{enumerate}
!      1) sulfate
!      2) organic carbon from biomass burning
!      3) black   carbon from biomass burning
!      4) organic carbon from fossil fuel
!      5) black   carbon from fossil fuel
!      6) organic carbon from natural sources
!      7) sea salt 1    (sub-micron)
!      8) sea salt 2    (super-micron)
!      9) winter dust 1 (sub-micron)
!     10) winter dust 2 (super-micron)
!     11) summer dust 1 (sub-micron)
!     12) summer dust 2 (super-micron)
!     13) desert dust 1 (sub-micron)
!     14) desert dust 2 (super-micron)
!     15) stratospheric sulfate (monodispersed)
!  \end{enumerate}
!
! !LOCAL VARIABLES:
      type (t_GmiArrayBundle), pointer :: concentration(:)
      integer :: ic

#ifdef MICRO_AEROSOL
      integer :: il, ij, ik
#endif

      real*8  :: aer_den_swel(i1:i2, ju1:j2, k1+1:k2)
      real*8  :: c3t         (i1:i2, ju1:j2, k1+1:k2)
      real*8  :: del_grid_box(i1:i2, ju1:j2, k1+1:k2)
      real*8  :: kele        (i1:i2, ju1:j2, k1+1:k2)
      real*8  :: mfp         (i1:i2, ju1:j2, k1+1:k2)
      real*8  :: radius      (i1:i2, ju1:j2, k1+1:k2)
      real*8  :: relhume     (i1:i2, ju1:j2, k1+1:k2)
      real*8  :: slip        (i1:i2, ju1:j2, k1+1:k2)
      real*8  :: velocity    (i1:i2, ju1:j2, k1+1:k2)
      real*8  :: viscosity   (i1:i2, ju1:j2, k1+1:k2)

#ifdef MICRO_AEROSOL
!-micro_aerosol--------begin---------------------------------------------
!c    for adding umaerosol

      real*8, parameter :: CMPM3 = CMPM * CMPM * CMPM

      real*8  :: h2osat  (i1:i2, ju1:j2, k1+1:k2)
      real*8  :: h2ogas  (i1:i2, ju1:j2, k1+1:k2)
      real*8  :: so4mfrac(i1:i2, ju1:j2, k1+1:k2)
      real*8  :: so4dens (i1:i2, ju1:j2, k1+1:k2)
      real*8  :: wetmas  (i1:i2, ju1:j2, k1+1:k2)
      real*8  :: so4aer  (i1:i2, ju1:j2, k1+1:k2, naer)
      real*8  :: so4radv (i1:i2, ju1:j2, k1+1:k2, nso4)
#ifndef NONON
      real*8  :: so4non  (i1:i2, ju1:j2, k1+1:k2, nnon)
      real*8  :: aernon  (i1:i2, ju1:j2, k1+1:k2, nnon)
      real*8  :: xnonum  (i1:i2, ju1:j2, k1+1:k2, nnon)
      real*8  :: so4vnon (i1:i2, ju1:j2, k1+1:k2, nnon)
#endif
!-micro_aerosol--------end-----------------------------------------------
#endif

!EOP
!------------------------------------------------------------------------------
!BOC

!     if (pr_diag) then
!       Write (6,*) 'Update_Grav_Settling called by ', loc_proc
!     end if

      call Get_concentration(SpeciesConcentration, concentration)

!     -------------------------------------------------------
!     Make sure aerosols are set correctly for radon/lead (1)
!     and beryllium (6) chemistry options.
!     -------------------------------------------------------

      if (chem_opt == 1) then

        aerosol(1) = 0
        aerosol(2) = 1

      else if (chem_opt == 6) then

        aerosol(1) = 15
        aerosol(2) = 15

      end if

!     ------------------------------------------------------------
!     Outside of the species loop, calculate all aerodynamic terms
!     which are not dependent on physical characteristics of each
!     aerosol type.
!     ------------------------------------------------------------

!     ----------------------------------------
!     Calculate temperature at grid interface.
!     ----------------------------------------

      kele(:,:,:) = (kel(i1:i2,ju1:j2,k1:k2-1) +  &
     &               kel(i1:i2,ju1:j2,k1+1:k2)) * 0.5d0

!     -----------------------------------------------------------------
!     First calculate relative humidity from Seinfeld (1986) p. 181.
!     The first  relhume is the temperature dependent parameter a.
!     The second relhume is the saturation vapor pressure of water.
!     The third  relhume is the actual relative humidity as a fraction.
!     Then make sure relhume is between 0 and 0.95 because swelling
!     parameterization goes too large at very high relative humidity.
!     -----------------------------------------------------------------

      relhume(:,:,:) = 1.0d0 - (373.15d0 / kele(:,:,:))

      relhume(:,:,:) =  &
     &  1013.25d0 * Exp (13.3185d0 * relhume(:,:,:)    -  &
     &                    1.9760d0 * relhume(:,:,:)**2 -  &
     &                    0.6445d0 * relhume(:,:,:)**3 -  &
     &                    0.1299d0 * relhume(:,:,:)**4)

      relhume(:,:,:) =  &
     &  (humidity(:,:,1:k2-1) + humidity(:,:,2:k2)) *  &
     &  0.5d0 * MWTAIR / 18.0d0 /  &
     &  GPKG * press3e(i1:i2,ju1:j2,k1:k2-1) / relhume(:,:,:)

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
     &  (PASPMB * press3e(i1:i2,ju1:j2,k1:k2-1) *  &
     &   Sqrt (8.0d0 * MWTAIR /  &
     &         (GMI_PI * GAS_CONST_J * 1.0d3 * kele(:,:,:))))

#ifdef MICRO_AEROSOL
!-micro_aerosol--------begin---------------------------------------------
!c    --------------------
!c    for adding umaerosol
!c    --------------------

!c    compute wet radius of so4 aerosol

      so4aer(:,:,:,1) = (concentration(ISO4M1)%pArray3D(:,:,k1:k2-1)+  &
     &                   concentration(ISO4M1)%pArray3D(:,:,k1+1:k2)) * 0.5d0
      so4aer(:,:,:,2) = (concentration(ISO4N1)%pArray3D(:,:,k1:k2-1)+  &
     &                   concentration(ISO4N1)%pArray3D(:,:,k1+1:k2)) * 0.5d0
      so4aer(:,:,:,3) = (concentration(ISO4M2)%pArray3D(:,:,k1:k2-1)+  &
     &                   concentration(ISO4M2)%pArray3D(:,:,k1+1:k2)) * 0.5d0
      so4aer(:,:,:,4) = (concentration(ISO4N2)%pArray3D(:,:,k1:k2-1)+  &
     &                   concentration(ISO4N2)%pArray3D(:,:,k1+1:k2)) * 0.5d0

      do il = i1, i2
        do ij = ju1, j2
          do ik = k1+1, k2

            h2osat(il,ij,ik) = h2osat_f(kele(il,ij,ik))
            h2ogas(il,ij,ik) = relhume(il,ij,ik)*h2osat(il,ij,ik)

            so4mfrac(il,ij,ik) = so4mfrac_f(kele(il,ij,ik),  &
     &                                      h2osat(il,ij,ik),  &
     &                                      h2ogas(il,ij,ik))
            so4dens(il,ij,ik) = so4dens_f(kele(il,ij,ik),  &
     &                                    so4mfrac(il,ij,ik))

          end do
        end do
      end do

      do iso4 = 1, nso4

        iso4n = iso4 * nmomso4
        iso4m = iso4n - 1

        wetmas(:,:,:) = max(r2so4min, so4aer(:,:,:,iso4m)  &
     &                / max(epsilo,   so4aer(:,:,:,iso4n)))  &
     &                / so4mfrac(:,:,:)

        so4radv(:,:,:,iso4) = (r3q * wetmas(:,:,:)  &
     &                      / (GMI_PI * so4dens(:,:,:)))**r1td
        so4radv(:,:,:,iso4) = max(so4radvmin,  &
     &                        min(so4radv(:,:,:,iso4),so4radvmax))

      end do

#ifndef NONON
      do ic = 1, nnon
         ! from ISO4NOC to ISO4S4
         so4non(:,:,:,ic) =  &
     &      (concentration(ISO4NOC+ic-1)%pArray3D(:,:,k1:k2-1) +  &
     &       concentration(ISO4NOC+ic-1)%pArray3D(:,:,k1+1:k2)) * 0.5d0

         ! from INOC to ISSLT4
         aernon(:,:,:,ic) =  &
     &      (concentration(INOC+ic-1)%pArray3D(:,:,k1:k2-1) +  &
     &       concentration(INOC+ic-1)%pArray3D(:,:,k1+1:k2)) * 0.5d0
      end do

!c    non-so4 aerosol radius (m) and mass (kg/particle)
      do inon = 1, nnon

        radvolm(:) = radgnon(:,inon)  &
     &             * exp(r3h*log(siggnon(:,inon))**2)
        radvnon(inon) = sum(fracnon(:,inon)*radvolm(:)**3)
        pmsnon(inon) = r4td * GMI_PI * rhonon(inon)  &
     &               * radvnon(inon)
        radvnon(inon) = radvnon(inon)**r1td

      end do

!c    non-so4 aerosol number concentration (#particles/kg air)
      do inon = 1, nnon
        xnonum(:,:,:,inon) = aernon(:,:,:,inon) / pmsnon(inon)
      end do

!c    so4 volume on each non-so4 aerosol surface (m3/particle)
      do inon = 1, nnon

        so4vnon(:,:,:,inon) =  &
     &    max(1.0d-30,so4non(:,:,:,inon))  &
     &    / max(r1,xnonum(:,:,:,inon))  &
     &    / (r4td*GMI_PI*so4mfrac(:,:,:)*so4dens(:,:,:))

      end do
#endif
!-micro_aerosol--------end-----------------------------------------------
#endif

!     ==============================
      ICLOOP: do ic = 1, num_species
!     ==============================

!       -------------------------------------------------
!       Loop over each species and do settling when not a
!       fixed species and when it is an aerosol.
!       -------------------------------------------------

        if (isFixedConcentration(ic) .or.  &
     &     (aerosol(ic) <= 0)  .or.  &
     &     (aerosol(ic) >  NUM_AEROSOL)) then
!         ============
          cycle ICLOOP
!         ============
        end if

!       -----------------------------------------------------------
!       Swelling of aerosol due to humidity
!       (Gong, Barrie, and Blanchet 1997 JGR 3805-3818 equation 3).
!       -----------------------------------------------------------
#ifndef MICRO_AEROSOL

        if (C1(aerosol(ic)) == 0.0d0) then

          radius(:,:,:) = RADIUS_EFF(aerosol(ic))

          aer_den_swel(:,:,:) = AER_DENSITY(aerosol(ic))

        else

          c3t(:,:,:) =  &
     &      C3(aerosol(ic)) *  &
     &      (1.0d0 + 0.004d0 * (298.0d0 - kele(:,:,:)))

          radius(:,:,:) =  &
     &      (C1(aerosol(ic)) *  &
     &       (RADIUS_EFF(aerosol(ic))*CMPM)**C2(aerosol(ic)) /  &
     &      (c3t(:,:,:) *  &
     &       (RADIUS_EFF(aerosol(ic))*CMPM)**C4(aerosol(ic)) -  &
     &      Log10 (relhume(:,:,:))) +  &
     &      (RADIUS_EFF(aerosol(ic))*CMPM)**3)**(1.0d0/3.0d0)

          aer_den_swel(:,:,:) = ((AER_DENSITY(aerosol(ic)) - 1.0d3) *  &
     &      (RADIUS_EFF(aerosol(ic))*CMPM)**3 + 1.0d3 *  &
     &      radius(:,:,:)**3) / radius(:,:,:)**3

          radius(:,:,:) = radius(:,:,:) / CMPM

        end if
#else
!-micro_aerosol--------begin---------------------------------------------
!       -------------------------------------------------
!       Calculate radius of sulfate for the gravitational
!       settlement of mass and number
!       -------------------------------------------------

        if (aerosol(ic) == 1 .or. aerosol(ic) == 2) then

          xlnsg = log(sigmod(aerosol(ic)))
          xlnsg2 = xlnsg * xlnsg

          if (ic == ISO4M1 .or. ic == ISO4M2) then

            radius(:,:,:) = so4radv(:,:,:,aerosol(ic))  &
     &                    * exp(r5h * xlnsg2)

          elseif (ic == ISO4N1 .or. ic == ISO4N2) then

            radius(:,:,:) = so4radv(:,:,:,aerosol(ic))  &
     &                    * exp(-r1h * xlnsg2)

          end if

          aer_den_swel(:,:,:) = so4dens(:,:,:)

        else

!       -----------------------------------------------------------
!       Swelling of aerosol due to humidity
!       (Gong, Barrie, and Blanchet 1997 JGR 3805-3818 equation 3).
!       -----------------------------------------------------------

#ifdef NONON
          if (C1(aerosol(ic)) == 0.0d0) then

            radius(:,:,:) = RADIUS_EFF(aerosol(ic))

            aer_den_swel(:,:,:) = AER_DENSITY(aerosol(ic))

          else
            c3t(:,:,:) =  &
     &        C3(aerosol(ic)) *  &
     &        (1.0d0 + 0.004d0 * (298.0d0 - kele(:,:,:)))

            radius(:,:,:) =  &
     &        (C1(aerosol(ic)) *  &
     &         (RADIUS_EFF(aerosol(ic))*CMPM)**C2(aerosol(ic)) /  &
     &        (c3t(:,:,:) *  &
     &         (RADIUS_EFF(aerosol(ic))*CMPM)**C4(aerosol(ic)) -  &
     &        Log10 (relhume(:,:,:))) +  &
     &        (RADIUS_EFF(aerosol(ic))*CMPM)**3)**(1.0d0/3.0d0)

            aer_den_swel(:,:,:) = ((AER_DENSITY(aerosol(ic)) - 1.0d3) *  &
     &        (RADIUS_EFF(aerosol(ic))*CMPM)**3 + 1.0d3 *  &
     &        radius(:,:,:)**3) / radius(:,:,:)**3

            radius(:,:,:) = radius(:,:,:) / CMPM

          end if

#else
          inon = aerosol(ic) - nso4

          if (C1(aerosol(ic)) == 0.0d0) then

            radius(:,:,:) =  &
     &        (RADIUS_EFF(aerosol(ic))**3 +  &
     &         so4vnon(:,:,:,inon))**(1.0d0/3.0d0)

            aer_den_swel(:,:,:) = (AER_DENSITY(aerosol(ic)) *  &
     &        RADIUS_EFF(aerosol(ic))**3 +  &
     &        so4dens(:,:,:) * so4vnon(:,:,:,inon)) /  &
     &        radius(:,:,:)**3

          else

            c3t(:,:,:) =  &
     &        C3(aerosol(ic)) *  &
     &        (1.0d0 + 0.004d0 * (298.0d0 - kele(:,:,:)))

            radius(:,:,:) =  &
     &        (C1(aerosol(ic)) *  &
     &         (RADIUS_EFF(aerosol(ic))*CMPM)**C2(aerosol(ic)) /  &
     &        (c3t(:,:,:) *  &
     &         (RADIUS_EFF(aerosol(ic))*CMPM)**C4(aerosol(ic)) -  &
     &        Log10 (relhume(:,:,:))) +  &
     &        (RADIUS_EFF(aerosol(ic))*CMPM)**3 +  &
     &        so4vnon(:,:,:,inon)*CMPM3)**(1.0d0/3.0d0)

            aer_den_swel(:,:,:) = ((AER_DENSITY(aerosol(ic)) - 1.0d3) *  &
     &        (RADIUS_EFF(aerosol(ic))*CMPM)**3 + 1.0d3 *  &
     &         radius(:,:,:)**3 +  &
     &        (so4dens(:,:,:) - 1.0d3) * so4vnon(:,:,:,inon)*CMPM3) /  &
     &        radius(:,:,:)**3

            radius(:,:,:) = radius(:,:,:) / CMPM

          endif
!-micro_aerosol--------end-----------------------------------------------
        end if
#endif
#endif

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
!       and velocity at k = 2
!       ------------------------------------------------------

        diffaer(:,:,ic) = BOLTZMN_J * kele(:,:,2) * slip(:,:,2) /  &
     &    (6.0d0 * GMI_PI * viscosity(:,:,2) * radius(:,:,2))

        s_radius(:,:,ic) = radius(:,:,2)

        s_velocity(:,:,ic) = velocity(:,:,2)

!       --------------------------------------------------------
!       Apply the settling velocity making sure that it does not
!       exceed the Courant limit.
!       --------------------------------------------------------

        del_grid_box(:,:,:) =  &
     &    velocity(:,:,:) * tdt4 / grid_height(:,:,k1+1:k2)

        where (del_grid_box(:,:,:) > 1.0d0)  &
     &    del_grid_box(:,:,:) = 1.0d0

        del_grid_box(:,:,:) =  &
     &    del_grid_box(:,:,:) * concentration(ic)%pArray3D(:,:,k1+1:k2)

        concentration(ic)%pArray3D(:,:,k1:k2-1) =  &
     &    concentration(ic)%pArray3D(:,:,k1:k2-1) +  &
     &    (del_grid_box(:,:,k1+1:k2) *  &
     &     mass(:,:,k1+1:k2) / mass(:,:,k1:k2-1))

        concentration(ic)%pArray3D(:,:,k1+1:k2) =  &
     &    concentration(ic)%pArray3D(:,:,k1+1:k2) -  &
     &    del_grid_box(:,:,k1+1:k2)

!     =============
      end do ICLOOP
!     =============

      call Set_concentration(SpeciesConcentration, concentration)

      return

      end subroutine updateGravitationalSettling
!EOC
!------------------------------------------------------------------------------
      end module GmiGravitationalSettling_mod
