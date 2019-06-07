!-------------------------------------------------------------------------
!BOP
!
! !MODULE: ReadSolarCycle_mod
!
! !INTERFACE:
!
    module ReadSolarCycle_mod
!
    implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
!
  private
  public  :: readSolarCycleData
!
! !DESCRIPTION:
! Routine to set or read solar cycle scaling factor data.
!
! !REVISION HISTORY:
! Luke Oman 4 Aug 2016 Added Solar Cycle Read file
!
!EOP
!-------------------------------------------------------------------------

  CONTAINS

!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: readSolarCycleData
!
! !INTERFACE:
!
      subroutine readSolarCycleData ( s_cycle_dates,s_cycle, sc_infile_name)
!
! !USES:
!
      use GmiASCIIoperations_mod,   only : AsciiOpenRead
      use GmiFastJX_ParametersMod,  only : NW
!
      implicit none
!
#     include "GmiParameters.h"
!
! !INPUT PARAMETERS:
      real, intent(inOut), dimension(2628)    :: s_cycle_dates        ! 2628 months : 1882 - 2100
      character (len=MAX_LENGTH_FILE_NAME), intent(in) :: sc_infile_name
      real, intent(inOut), dimension(NW,2628) :: s_cycle              ! 2628 months : 1882 - 2100

! !DESCRIPTION:
! This routine sets/reads the solar cycle scaling values.
!
! !LOCAL VARIABLES:
      integer :: ic, im
      integer :: lun
      integer :: fjx_bin = NW
!
!-----------------------------------------------------------------------------
!BOC

         call AsciiOpenRead (lun, sc_infile_name)

         do ic = 1, 2628
                  Read (lun, 900) s_cycle_dates(ic), (s_cycle(im,ic), im = 1, fjx_bin)
         end do

 900     format (19f10.5)

         Close (lun)

      return

      end subroutine readSolarCycleData
!EOC
!---------------------------------------------------------------------------

    end module ReadSolarCycle_mod
