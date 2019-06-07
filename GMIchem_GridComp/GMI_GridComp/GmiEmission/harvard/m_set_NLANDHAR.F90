!-------------------------------------------------------------------------
!
! !MODULE: m_set_NLANDHAR
!
! !INTERFACE:
!
      module m_set_NLANDHAR
!
! !USES:
!
      implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
!
      private
      public NLANDHAR_expected
!
! !PUBLIC MEMBER FUNCTIONS:
      public NLANDHAR

      integer :: NLANDHAR
!
! !DESCRIPTION:
!  This module contains a function to set the value
!  of NLANDHAR as function of the horizontal resolution.
!
! !AUTHOR:
!  Jules Kouatchou
!
! !REVISION HISTORY:
!  19August2005 Initial Code
!EOP
!-------------------------------------------------------------------------

      CONTAINS

!-------------------------------------------------------------------------
!BOP
!
! !IFUNCTION: NLANDHAR_expected
!
! !INTERFACE:
!
      integer function NLANDHAR_expected(horz_resolution)
!
      implicit none
!
! !INPUT PARAMETERS:
      integer, intent(in) :: horz_resolution
!
! !DEFINED PARAMETERS:
      integer, parameter :: NLANDHAR_1DEG = 17822 ! 1x1.25 horizontal resolution
      integer, parameter :: NLANDHAR_2DEG =  3920 ! 2x2.5  horizontal resolution
      integer, parameter :: NLANDHAR_4DEG =  1118 ! 4x5    horizontal resolution
!
! !DESCRIPTION:
!  This function get the global number of grid points (depends on the
!  horizontal resolution) in the longitude direction and returns the desired
!  value of NLANDHAR.
!
! !AUTHOR:
!  Jules Kouatchou
!
! !REVISION HISTORY:
!  19August2005  Original code.
!
!EOP
!-------------------------------------------------------------------------
!BOC

      Select Case (horz_resolution) ! horizontal resolution
      Case (288)
         NLANDHAR_expected = NLANDHAR_1DEG
      Case (144)
         NLANDHAR_expected = NLANDHAR_2DEG
      Case (72)
         NLANDHAR_expected = NLANDHAR_4DEG
      Case Default
         write (6,*) 'i2_gl = ', horz_resolution,' has an illegal value.'
         stop "Code stop in the function NLANDHAR_expected"
      End Select

      return
      end function NLANDHAR_expected
!EOC
!-------------------------------------------------------------------------

      end module m_set_NLANDHAR
