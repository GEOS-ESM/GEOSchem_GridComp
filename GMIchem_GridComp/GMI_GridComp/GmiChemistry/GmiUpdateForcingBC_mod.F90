!-------------------------------------------------------------------------
! NASA GSFC - SIVO Code 610.3
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: GmiUpdateForcingBC_mod
!    
      module GmiUpdateForcingBC_mod
!    
! !USES:
      use GmiArrayBundlePointer_mod, only : t_GmiArrayBundle
      use GmiTimeControl_mod, only : GmiSplitDateTime
!     
      implicit none
! 
      private
!
! !PUBLIC MEMBER FUNCTIONS:
      public  :: updateForcingBC
!     
#     include "GmiParameters.h"
#     include "gmi_phys_constants.h"
#     include "gmi_time_constants.h"
!
! !AUTHOR:
!  John Tannahill, LLNL, jrt@llnl.gov
!  Jules Kouatchou, NASA/GSFC, Jules.Kouatchou-1@nasa.gov
!
! !REVISION HISTORY:
!   - October 31, 2007: Eric Nielsen
!     Modified determination of latitude nearest available boundary forcing
!     in Update_Forc_Bc to accomodate GEOS-5's ju1 = 1 always.
!   - April 9, 2010: Eric Nielsen
!     Allow possibility forc_bc_map(ic) <= 0
!   - December 15, 2011: Eric Nielsen
!     Generalized to 2D latitude input
!     
!EOP  
!-------------------------------------------------------------------------
      contains 
!-------------------------------------------------------------------------
!BOP  
!     
! !IROUTINE: updateForcingBC
!
! !INTERFACE:
!
      subroutine updateForcingBC (forc_bc_data, concentration, &
                       jlatmd, last_year, nymd,       &
     &                 gmi_sec,  fbc_j1, fbc_j2, forc_bc_num, forc_bc_kmax,    &
     &                 forc_bc_kmin, forc_bc_opt, forc_bc_map, forc_bc_incrpyr,&
     &                 forc_bc_start_num, forc_bc_years, &
     &                 i1, i2, ju1, j2, k1, k2, num_species)
!
      implicit none
!
#     include "gmi_forc_bc.h"
!
! !INPUT PARAMETERS:
      integer, intent(in) :: i1, i2, ju1, j2, k1, k2
      integer, intent(in) :: num_species
      integer, intent(in) :: forc_bc_kmax, forc_bc_kmin
      integer, intent(in) :: forc_bc_num
      integer, intent(in) :: forc_bc_opt
      integer, intent(in) :: forc_bc_years
      integer, intent(in) :: fbc_j1, fbc_j2
      integer, intent(in) :: forc_bc_map (forc_bc_num)
      real*8 , intent(in) :: forc_bc_incrpyr
      integer, intent(in) :: nymd
      real*8 , intent(in) :: gmi_sec
                             ! latitude of zone center in latitude direction
                             ! (rad)
      integer, intent(in) :: jlatmd(i1:i2,ju1:j2)
!
! !INPUT/OUTPUT PARAMETERS:
      integer, intent(inOut) :: last_year
      real*8                 forc_bc_data(:,:,:,:)
      integer, intent(inOut) :: forc_bc_start_num
                               ! species concentration, known at zone centers 
                               ! (mixing ratio)
      type (t_GmiArrayBundle), intent(inout) :: concentration(num_species)
!
! !DESCRIPTION:
!  This routine updates the forcing boundary conditions.
!                             
! !LOCAL VARIABLES:
      integer :: idumday
      integer :: ic, ii, ij
      integer :: jj, jx
      integer :: month
      integer :: year
      real*8  :: fbc_emiss_fac
      real*8  :: fbc_emission
      real*8  :: rjx
!
! !REVISION HISTORY:
!  Original code
!  Eric Nielsen 25 Jun 2008: Revision for option 1, which will now use the  values from year 
!                            forc_bc_start_num (of years 1 through forc_bc_years).
!
!  Eric Nielsen  2 Jan 2015: Revision for option 2: forc_bc_start_num is the year number of the 
!                            first year on the dataset, e.g. 1970.  Eliminated is the need to 
!                            update the RC file during multiple-year integrations with checkpoints.
!
!EOP
!------------------------------------------------------------------------------
!BOC

      call GmiSplitDateTime (nymd, year, month, idumday)

!     =====================
      if (forc_bc_opt == 1) then
!     =====================

        do ic = 1, forc_bc_num
          IF(forc_bc_map(ic) > 0) THEN
          do ij = ju1, j2
           do ii = i1, i2

            concentration(forc_bc_map(ic))%pArray3D(ii,ij,forc_bc_kmin:forc_bc_kmax) =  &
     &        forc_bc_data(jlatmd(ii,ij),month,forc_bc_start_num,ic)

           end do
          end do
          END IF
        end do


!     ==========================
      else if (forc_bc_opt == 2) then
!     ==========================

! ------------------------------- N O T I C E -------------------------------
! In this case forc_bc_start_num is the first year on the dataset, e.g. 1970.
! ---------------------------------------------------------------------------

        jx = year - forc_bc_start_num + 1
        IF(jx <             1) jx = 1
        IF(jx > forc_bc_years) jx = forc_bc_years

        do ic = 1, forc_bc_num
          IF(forc_bc_map(ic) > 0) THEN
          do ij = ju1, j2
           do ii = i1, i2

            concentration(forc_bc_map(ic))%pArray3D(ii,ij,forc_bc_kmin:forc_bc_kmax) =  &
     &        forc_bc_data(jlatmd(ii,ij),month,jx,ic)

           end do
          end do
          END IF
        end do


!     ==========================
      else if (forc_bc_opt == 3) then
!     ==========================

        fbc_emiss_fac = (gmi_sec / SECPYR) + (forc_bc_start_num - 1)
        fbc_emission  = fbc_emiss_fac * forc_bc_incrpyr * PPT_FAC

        do ic = 1, forc_bc_num
          IF(forc_bc_map(ic) > 0) THEN
          do ij = ju1, j2

            if ((ij >= fbc_j1) .and. (ij <= fbc_j2)) then

              concentration(forc_bc_map(ic))%pArray3D(:,ij,forc_bc_kmin:forc_bc_kmax)  &
     &          = fbc_emission

            end if

          end do
          END IF
        end do

      end if


      return
      end subroutine updateForcingBC
!EOC
!------------------------------------------------------------------------------

      end module GmiUpdateForcingBC_mod
