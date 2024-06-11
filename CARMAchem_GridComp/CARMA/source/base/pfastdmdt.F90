! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!!  This routine normalizes dmdt as to not make gas concentration go negative.
!!
!! @author Parker Case
!! @version Mar-2024
subroutine pfastdmdt(carma, cstate, iz, rc)

  ! types
  use carma_precision_mod
  use carma_enums_mod
  use carma_constants_mod
  use carma_types_mod
  use carmastate_mod
  use carma_mod

  implicit none

  type(carma_type), intent(in)         :: carma       !! the carma object
  type(carmastate_type), intent(inout) :: cstate      !! the carma state object
  integer, intent(in)                  :: iz          !! z index
  integer, intent(inout)               :: rc          !! return code, negative indicates failure

  ! Local declarations
  real(kind=f)  :: rvap
  real(kind=f)  :: gc_cgs
  real(kind=f)  :: gc_target
  real(kind=f)  :: scalefactor
  integer       :: igroup
  integer       :: ielem
  integer       :: igas
  integer       :: ibin

  do igroup = 1,NGROUP

    ielem = ienconc(igroup)     ! element of particle number concentration

    igas = igrowgas(ielem)      ! condensing gas

    if ((itype(ielem) == I_VOLATILE) .and. (igas /= 0)) then

      ! Calculate vapor pressures.
      call vaporp(carma, cstate, iz, igas, rc)

      ! Define gas constant for this gas
      rvap = RGAS / gwtmol(igas)

      ! Current gas concentration
      gc_cgs = gc(iz,igas) / (zmet(iz)*xmet(iz)*ymet(iz))

      ! Gas concentration at equilibrium
      gc_target = pvapl(iz,igas) / (rvap * t(iz))

      ! Determine total mass added to bins in implicit timestep
      scalefactor = sum((pc(iz,:,ielem) - pcl(iz,:,ielem))*rmass(:,igroup))/(gc_cgs - gc_target)

      ! Loop through bins and apply correction
      do ibin = 1,NBIN
        pc(iz,ibin,ielem) = pcl(iz,ibin,ielem) + (pc(iz,ibin,ielem) - pcl(iz,ibin,ielem)) / scalefactor
      end do
    end if

  end do

  return
end
