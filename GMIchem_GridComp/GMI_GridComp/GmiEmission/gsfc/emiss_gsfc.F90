
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
     &   pr_diag, loc_proc, i1, i2, ju1, j2, k1, k2, num_species)

      use GmiArrayBundlePointer_mod, only : t_GmiArrayBundle
      use gcr_mod                    , only : NUM_GCR_YEARS
      use gcr_mod                    , only : gcr_date
      use gcr_mod                    , only : gcr_slope
      use gcr_mod                    , only : gcr_aintcp
      use gcr_mod                    , only : gcr_sunspot

      implicit none

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

!     ----------------------
!     Variable declarations.
!     ----------------------

      integer :: idxyr, inyr
      integer :: il, ij, ik, iij, iik

      real*8  :: temp


!     ----------------
!     Begin execution.
!     ----------------

      if (pr_diag) then
        Write (6,*) 'Add_Emiss_Gsfc called by ', loc_proc
      endif

!     -------------------------------------------
!     Galactic Cosmic Ray emission of N and NO
!     -------------------------------------------
      inyr = int(nymd/10000)
      if(inyr > 69) then
        inyr = inyr+1900
       else if(inyr < 70 .and. inyr < 100) then
        inyr = inyr+2000
       endif
      idxyr = inyr-gcr_date(1)+1
!... make sure indice is in range
      do while (idxyr <= 0)
        idxyr = idxyr+11
      enddo
      do while (idxyr >= NUM_GCR_YEARS)
        idxyr = idxyr-11
      enddo

      do ik = k1, k2

        do ij = ju1, j2
!... calc total emission of NOx (normalize to indices starting at 1)
          temp = gcr_sunspot(idxyr) * gcr_slope(ij,ik)  &
     &           + gcr_aintcp(ij,ik)

          do il = i1, i2
!... 55% of emission goes to N
            if(IAN.gt.0)  &
     &        concentration(IAN)%pArray3D(il,ij,ik) = concentration(IAN)%pArray3D(il,ij,ik) +  &
     &          (0.55 * temp * tdt / concentration(IMGAS)%pArray3D(il,ij,ik) )

!... 45% of emission goes to NO
            if(INO.gt.0)  &
     &        concentration(INO)%pArray3D(il,ij,ik) = concentration(INO)%pArray3D(il,ij,ik) +  &
     &          (0.45 *temp *tdt / concentration(IMGAS)%pArray3D(il,ij,ik) )

          enddo
        enddo
      enddo

      return
      end
