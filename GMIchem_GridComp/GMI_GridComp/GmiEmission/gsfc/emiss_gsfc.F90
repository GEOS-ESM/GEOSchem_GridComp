!=============================================================================
!
! $Id$
!
! CODE DEVELOPER
!   Stephen Steenrod ; GSFC
!   steenrod@code916.gsfc.nasa.gov
!
! FILE
!   emiss_gsfc.F
!
! ROUTINES
!   Add_Emiss_Gsfc
!
!=============================================================================


!-----------------------------------------------------------------------------
!
! ROUTINE
!   Add_Emiss_Gsfc
!
! DESCRIPTION
!   This routine adds emissions to const.
!
!     1) Take emissions of moecules/s and multiply by the time step to get total
!        molecules of emission over time step.
!
!     2) Divide by mass of the zone to obtain emissions in terms of mixing
!        ratio; also multiply by the ratio of molecular weight of air to
!        molecular weight of the chemical emission to get volume mixing
!        ratio from mass mixing ratio.
!
!     3) Add emitted mixing ratio amount to existing mixing ratio of const.
!
! ARGUMENTS
!   mass           : total mass of the atmosphere within each grid box (kg)
!   const          : species concentration, known at zone centers
!                    (mixing ratio)
!   nymd           : model date
!   tdt            : model time step
!-----------------------------------------------------------------------------

      subroutine Add_Emiss_Gsfc  &
     &  (mass, concentration, nymd, tdt, &
     &   IAN, IMGAS, INO, &
     &   pr_diag, loc_proc, i1, i2, ju1, j2, k1, k2, num_species, &
     &   latdeg, press3c, kel )

      USE ESMF
      USE MAPL

      use GmiArrayBundlePointer_mod, only : t_GmiArrayBundle
      use gcr_mod                  , only : USE_GCR_DATA
      use gcr_mod                  , only : SET_GCR_EMISS

      implicit none

#include "gmi_phys_constants.h"      
!     ----------------------
!     Argument declarations.
!     ----------------------

      integer, intent(in   ) :: IAN, IMGAS, INO
      logical, intent(in   ) :: pr_diag
      integer, intent(in   ) :: loc_proc
      integer, intent(in   ) :: i1, i2, ju1, j2, k1, k2, num_species
      integer, intent(in   ) :: nymd
      real*8 , intent(in   ) :: tdt
      real*8 , intent(in   ) :: mass(i1:i2, ju1:j2, k1:k2)
      type (t_GmiArrayBundle), intent(inout) :: concentration(num_species)
      real*8 , intent(in   ) :: latdeg   (i1:i2, ju1:j2)
!     real   , intent(out  ) ::  gcrnox(i1:i2, ju1:j2, k1:k2)
      real*8 , intent(in   ) :: press3c(i1:i2, ju1:j2, k1:k2) ! hPa
      real*8 , intent(in   ) ::     kel(i1:i2, ju1:j2, k1:k2)


!     ----------------------
!     Variable declarations.
!     ----------------------

      integer :: year
      integer :: i, j, k

!     USE_GCR_DATA() will return values in these:
      REAL*8              :: sunspot
      REAL*8, ALLOCATABLE :: slope(:,:,:)
      REAL*8, ALLOCATABLE :: aintcp(:,:,:)

      REAL*8, ALLOCATABLE :: gcr_nox_em_rate(:,:,:)
      REAL*8, ALLOCATABLE :: air_dens(:,:,:)


!     ----------------
!     Begin execution.
!     ----------------

      if (pr_diag) then
        Write (6,*) 'Add_Emiss_Gsfc called by ', loc_proc
      endif

!     -------------------------------------------
!     Galactic Cosmic Ray emission of N and NO
!     -------------------------------------------
      year = int(nymd/10000)

      ALLOCATE(             slope(i1:i2, ju1:j2, k1:k2) )
      ALLOCATE(            aintcp(i1:i2, ju1:j2, k1:k2) )
      ALLOCATE(   gcr_nox_em_rate(i1:i2, ju1:j2, k1:k2) )

       slope(:,:,:)=0.0
      aintcp(:,:,:)=0.0

      CALL USE_GCR_DATA( year, press3c, latdeg, &
                         SUNSPOT, SLOPE, AINTCP, &
                         i1, i2, ju1, j2, k1, k2 )
!... This should be (molecules of NOx)/cm3/sec
      gcr_nox_em_rate = sunspot * slope + aintcp
!...   save it for later
      CALL  SET_GCR_EMISS ( gcr_nox_em_rate, &
                             i1, i2, ju1, j2, k1, k2 )
!... NOTE: IMGAS is in molec/cm3

!... Douglass/Oman/Manyin 9/11/13
!...   After seeing significant artifacts in N as a result of adding 55%
!...   (chemistry did not converge in nighttime areas)
!...   decided to shift 100% to NO

!... 55% of emission goes to N
!           if (IAN.gt.0) then
!             concentration(IAN)%pArray3D(i,j,k) = &
!             concentration(IAN)%pArray3D(i,j,k) + &
!               (0.55 * gcr_nox_em_rate(i,j,k) * tdt / concentration(IMGAS)%pArray3D(i,j,k) )
!           end if

!... 45% of emission goes to NO
            if (INO.gt.0) then
              concentration(INO)%pArray3D(i1:i2,ju1:j2,k1:k2) = &
              concentration(INO)%pArray3D(i1:i2,ju1:j2,k1:k2) + &
                ( gcr_nox_em_rate(i1:i2,ju1:j2,k1:k2) * tdt / concentration(IMGAS)%pArray3D(i1:i2,ju1:j2,k1:k2) )
            end if

      DEALLOCATE(    slope        )
      DEALLOCATE(   aintcp        )
      DEALLOCATE( gcr_nox_em_rate )

      return
      end
