#include "MAPL_Generic.h"

module abortutils

  use cam_logfile, only: iulog
  use MAPL_Mod

  implicit none
  private

  integer, public, parameter :: CAM_CRITICAL_ERROR = 2014

  public :: endrun

contains
  
  subroutine endrun(msg,rc)

      character(len=*), optional, intent(in) :: msg
      integer, optional, intent(out)    :: rc

      __Iam__('CAM::endrun()')
      
      character(len=256) :: message

      if (present(msg)) then
          message = trim(msg)
      else
          message = 'Encountered critical error in the CAM aerosol microphysics code.'
      endif

      if (present(rc)) then
          rc = CAM_CRITICAL_ERROR
      end if

      __raise__(CAM_CRITICAL_ERROR, trim(message))
  end subroutine

end module abortutils
