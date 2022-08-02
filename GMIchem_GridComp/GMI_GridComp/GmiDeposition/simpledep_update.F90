
!=============================================================================
!
! $Id$
!
! CODE DEVELOPER
!   John Tannahill, LLNL
!   jrt@llnl.gov
!
! FILE
!   simpledep_update.F
!
! ROUTINES
!   Update_Simpledep
!
!=============================================================================


!-----------------------------------------------------------------------------
!
! ROUTINE
!   Update_Simpledep
!
! DESCRIPTION
!   This routine does the parameterized deposition.
!
! ARGUMENTS
!   ibrono2_num   : const array index for BRONO2
!   ih2o2_num     : const array index for H2O2
!   ihcl_num      : const array index for HCl
!   ihno3_num     : const array index for HNO3
!   imgas_num     : const array index for air density
!   initrogen_num : const array index for nitrogen
!   ioxygen_num   : const array index for oxygen
!   num_ks_sdep   : number of vertical layers to apply 2 day loss factor to
!   tdt     : model time step (s)
!   press3c : atmospheric pressure at the center of each grid box (mb)
!   const   : species concentration, known at zone centers (mixing ratio)
!
!-----------------------------------------------------------------------------

      subroutine Update_Simpledep  &
     &  (ibrono2_num, ih2o2_num, ihcl_num, ihno3_num, imgas_num,  &
     &   initrogen_num, ioxygen_num, num_ks_sdep, tdt, press3c, concentration, &
     &   pr_diag, loc_proc, &
     &   i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi, num_species)

      use GmiArrayBundlePointer_mod, only : t_GmiArrayBundle
      use GmiSpcConcentrationMethod_mod, only : isFixedConcentration

      implicit none

#     include "gmi_time_constants.h"

!     ----------------------
!     Argument declarations.
!     ----------------------

      logical, intent(in   ) :: pr_diag
      integer, intent(in   ) :: loc_proc
      integer, intent(in   ) :: i1, i2, ju1, j2, k1, k2
      integer, intent(in   ) :: ilo, ihi, julo, jhi
      integer, intent(in   ) :: num_species
      integer, intent(in   ) :: ibrono2_num
      integer, intent(in   ) :: ih2o2_num
      integer, intent(in   ) :: ihcl_num
      integer, intent(in   ) :: ihno3_num
      integer, intent(in   ) :: imgas_num
      integer, intent(in   ) :: initrogen_num
      integer, intent(in   ) :: ioxygen_num
      integer, intent(in   ) :: num_ks_sdep
      real*8 , intent(in   ) :: tdt
      real*8 , intent(in   ) :: press3c(ilo:ihi, julo:jhi, k1:k2)
      type (t_GmiArrayBundle), intent(inout) :: concentration(num_species)

!     ----------------------
!     Variable declarations.
!     ----------------------

      logical, save :: first   = .true.

      integer :: ic

      real*8  :: rsecpdy

      real*8,  save :: loss_1day, loss_2day


!     ----------------
!     Begin execution.
!     ----------------

      if (pr_diag) then
        Write (6,*) 'Update_Simpledep called by ', loc_proc
      end if


      rsecpdy = SECPDY


!     ==========
      if (first) then
!     ==========

        first = .false.

        loss_1day = Exp (-tdt / rsecpdy)
        loss_2day = Exp (-tdt / (2.0d0 * rsecpdy))

      end if


      icloop: do ic = 1, num_species

!       ====================================
        if (isFixedConcentration(ic)) cycle icloop
!       ====================================

        if ((ic == ibrono2_num) .or.  &
     &      (ic == ih2o2_num)   .or.  &
     &      (ic == ihcl_num)    .or.  &
     &      (ic == ihno3_num)) then

!         ----------------------------------------------------------
!         Between 300 mb and 800 mb use a lifetime which is a third
!         order polynomial fit to the data P=800,700,600,500,400,300
!         and L=1,5,5,10,18,38 which results in the equation
!         L = 199.374 - .878334*P + .0013414*P^2 - 6.91682e-7*P^3.
!         ----------------------------------------------------------

          where ((press3c(i1:i2,ju1:j2,:) > 300.0d0) .and.  &
     &           (press3c(i1:i2,ju1:j2,:) < 800.0d0))

            concentration(ic)%pArray3D(:,:,:) =  &
     &        concentration(ic)%pArray3D(:,:,:) *  &
     &        Exp (-tdt /  &
     &             (rsecpdy *  &
     &              (199.374d0  -  &
     &               0.878334d0  * press3c(i1:i2,ju1:j2,:)    +  &
     &               0.0013414d0 * press3c(i1:i2,ju1:j2,:)**2 -  &
     &               6.91682d-7  * press3c(i1:i2,ju1:j2,:)**3)))

          end where

!         -----------------------------------------------
!         At pressures >= 800 mb use a lifetime of 1 day.
!         -----------------------------------------------

          where (press3c(i1:i2,ju1:j2,:) >= 800.d0)

            concentration(ic)%pArray3D(:,:,:) = concentration(ic)%pArray3D(:,:,:) * loss_1day

          end where

        else

          if ((ic /= initrogen_num) .and.  &
     &        (ic /= ioxygen_num)   .and.  &
     &        (ic /= imgas_num)) then

            concentration(ic)%pArray3D(:,:,1:num_ks_sdep) =  &
     &        concentration(ic)%pArray3D(:,:,1:num_ks_sdep) * loss_2day

          end if

        end if

      end do icloop


      return

      end subroutine Update_Simpledep

