!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  GmiCheckRange_mod
!
! !INTERFACE:
!
   module GmiCheckRange_mod
!
! !USES:
!
   implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
!
   private
   public  CheckRange1d
   public  CheckRange2d
   public  CheckRange3d
   public  CheckRange4d
   public  CheckRange5d
!
! !DESCRIPTION:
!  Checks the range of arrays and prints an error message an array
!  is out of bound.
!
! !AUTHOR: 
!  Jules Kouatchou
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------

CONTAINS

!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: CheckRange1d
!
! !INTERFACE:
!
      subroutine CheckRange1d (var_name, loc_proc, d11, d12, var_1d, loval, hival)
!
! !USES:
!
      use GmiPrintError_mod
!
      implicit none
!
! !INPUT PARAMETERS:
!!  var_name : the string name of var_1d
!!  loc_proc : local processor #
!!  d11, d12 : dimensions for var_1d
!!  var_1d   : the specified 1D array
!!  loval    : the lowest  allowed value
!!  hival    : the highest allowed value
      character (len=*), intent(in)    :: var_name
      integer          , intent(in)    :: loc_proc
      integer          , intent(in)    :: d11, d12
      real*8           , intent(in)    :: var_1d(d11:d12)
      real*8           , intent(in)    :: loval, hival
!
! !DESCRIPTION:
!  Does a gross range check on the data values in a 1D array,
!  and exits if it finds one out of range.
!
! !LOCAL VARIABLES:
      character (len=75) :: err_msg

      integer :: max1(1), min1(1)
!
! !AUTHOR: 
!  John Tannahill (LLNL) and Jules Kouatchou
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC
      min1(:) = Minloc (var_1d(d11:d12))
      max1(:) = Maxloc (var_1d(d11:d12))

      min1(1) = min1(1) + d11 - 1
      max1(1) = max1(1) + d11 - 1
 
      if ((var_1d(min1(1)) < loval) .or. (var_1d(max1(1)) > hival)) then
        Write (6,*) min1(:), loc_proc
        Write (6,*) max1(:), loc_proc
 
        err_msg = 'Value out of range in CheckRange1d, ' // var_name
 
        call GmiPrintError (err_msg, .true., 0, 0, 0, 2, var_1d(min1(1)), var_1d(max1(1)))
      end if

      return

      end subroutine CheckRange1d
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: CheckRange2d
!
! !INTERFACE:
!
      subroutine CheckRange2d   &
        (var_name, loc_proc, d11, d12, d21, d22, var_2d, loval, hival)
!
! !USES:
!
      use GmiPrintError_mod
!
      implicit none
!
! !INPUT PARAMETERS:
!!    var_name : the string name of var_2d
!!    loc_proc : local processor #
!!    d11, d12, d21, d22 : dimensions for var_2d
!!    var_2d   : the specified 2D array
!!    loval    : the lowest  allowed value
!!    hival    : the highest allowed value
      character (len=*), intent(in)   :: var_name
      integer          , intent(in)   :: loc_proc
      integer          , intent(in)   :: d11, d12, d21, d22
      real*8           , intent(in)   :: var_2d(d11:d12, d21:d22)
      real*8           , intent(in)   :: loval, hival
!
! !DESCRIPTION:
!  Does a gross range check on the data values in a 2D array,
!  and exits if it finds one out of range.
!
! !LOCAL VARIABLES:
      character (len=75) :: err_msg

      integer :: max2(2), min2(2)
!
! !AUTHOR: 
!  John Tannahill (LLNL) and Jules Kouatchou
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC
      min2(:) = Minloc (var_2d(d11:d12,d21:d22))
      max2(:) = Maxloc (var_2d(d11:d12,d21:d22))

      min2(1) = min2(1) + d11 - 1
      min2(2) = min2(2) + d21 - 1

      max2(1) = max2(1) + d11 - 1
      max2(2) = max2(2) + d21 - 1
 
      if ((var_2d(min2(1),min2(2)) < loval) .or. (var_2d(max2(1),max2(2)) > hival)) then
        Write (6,*) min2(:), loc_proc
        Write (6,*) max2(:), loc_proc
 
        err_msg = 'Value out of range in CheckRange2d, ' // var_name
 
        call GmiPrintError (err_msg, .true., 0, 0, 0, 2,   &
           var_2d(min2(1),min2(2)), var_2d(max2(1),max2(2)))
      end if
 
      return

      end subroutine CheckRange2d
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: CheckRange3d
!
! !INTERFACE:
!

      subroutine CheckRange3d (var_name, loc_proc, d11, d12, d21, d22, &
                                 d31, d32, var_3d, loval, hival)
!
! !USES:
!
      use GmiPrintError_mod
!
      implicit none
!
! !INPUT PARAMETERS:
!!    var_name : the string name of var_3d
!!    loc_proc : local processor #
!!    d11, d12, d21, d22, d31, d32 : dimensions for var_3d
!!    var_3d   : the specified 3D array
!!    loval    : the lowest  allowed value
!!    hival    : the highest allowed value
      character (len=*) :: var_name
      integer :: loc_proc
      integer :: d11, d12, d21, d22, d31, d32
      real*8  :: var_3d(d11:d12, d21:d22, d31:d32)
      real*8  :: loval, hival
!
! !DESCRIPTION:
!   This routine does a gross range check on the data values in a 3D array,
!   and exits if it finds one out of range.
!
! !LOCAL VARIABLES:
      character (len=75) :: err_msg

      integer :: max3(3), min3(3)
!
! !AUTHOR: 
!  John Tannahill (LLNL) and Jules Kouatchou
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC
      min3(:) = Minloc (var_3d(d11:d12,d21:d22,d31:d32))
      max3(:) = Maxloc (var_3d(d11:d12,d21:d22,d31:d32))

      min3(1) = min3(1) + d11 - 1
      min3(2) = min3(2) + d21 - 1
      min3(3) = min3(3) + d31 - 1

      max3(1) = max3(1) + d11 - 1
      max3(2) = max3(2) + d21 - 1
      max3(3) = max3(3) + d31 - 1

 
      if ((var_3d(min3(1),min3(2),min3(3)) < loval) .or.   &
          (var_3d(max3(1),max3(2),max3(3)) > hival))     then
 
        Write (6,*) min3(:), loc_proc
        Write (6,*) max3(:), loc_proc
 
        err_msg = 'Value out of range in CheckRange3d, ' // var_name
 
        call GmiPrintError (err_msg, .true., 0, 0, 0, 2,         &
                         var_3d(min3(1),min3(2),min3(3)),     &
                         var_3d(max3(1),max3(2),max3(3)))
 
      end if

      return

      end subroutine CheckRange3d
!
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: CheckRange4d
!
! !INTERFACE:
!
      subroutine CheckRange4d (var_name, loc_proc, d11, d12, d21, &
                                 d22, d31, d32, d41, d42, var_4d, loval, hival)
!
! !USES:
!
      use GmiPrintError_mod
!
      implicit none
!
! !INPUT PARAMETERS:
!!  var_name : the string name of var_4d
!!  loc_proc : local processor #
!!  d11, d12, d21, d22, d31, d32, d41, d42 : dimensions for var_4d
!!  var_4d   : the specified 4D array
!!  loval    : the lowest  allowed value
!!  hival    : the highest allowed value
      character (len=*), intent(in)   :: var_name
      integer          , intent(in)    :: loc_proc
      integer          , intent(in)    :: d11, d12, d21, d22, d31, d32, d41, d42
      real*8           , intent(in)    :: var_4d(d11:d12, d21:d22, d31:d32, d41:d42)
      real*8           , intent(in)    :: loval, hival
!
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION
!   This routine does a gross range check on the data values in a 4D array,
!   and exits if it finds one out of range.
!
! !LOCAL VARIABLES:
      character (len=75) :: err_msg

      integer :: max4(4), min4(4)
!
! !AUTHOR: 
!  John Tannahill (LLNL) and Jules Kouatchou
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC 
      min4(:) = Minloc (var_4d(d11:d12,d21:d22,d31:d32,d41:d42))
      max4(:) = Maxloc (var_4d(d11:d12,d21:d22,d31:d32,d41:d42))

      min4(1) = min4(1) + d11 - 1
      min4(2) = min4(2) + d21 - 1
      min4(3) = min4(3) + d31 - 1
      min4(4) = min4(4) + d41 - 1

      max4(1) = max4(1) + d11 - 1
      max4(2) = max4(2) + d21 - 1
      max4(3) = max4(3) + d31 - 1
      max4(4) = max4(4) + d41 - 1

 
      if ((var_4d(min4(1),min4(2),min4(3),min4(4)) < loval) .or. &
          (var_4d(max4(1),max4(2),max4(3),max4(4)) > hival))    then
 
        Write (6,*) min4(:), loc_proc
        Write (6,*) max4(:), loc_proc
 
        err_msg = 'Value out of range in CheckRange4d, ' // var_name
 
        call GmiPrintError (err_msg, .true., 0, 0, 0, 2,              &
                         var_4d(min4(1),min4(2),min4(3),min4(4)),  &
                         var_4d(max4(1),max4(2),max4(3),max4(4)))
 
      end if

      return

      end subroutine CheckRange4d
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: CheckRange5d
!
! !INTERFACE:
!
      subroutine CheckRange5d (var_name, loc_proc, d11, d12, d21, d22,  &
                                 d31, d32, d41, d42, d51, d52, var_5d,    &
                                 loval, hival)
!
! !USES:
!
      use GmiPrintError_mod
!
      implicit none
!
! !INPUT PARAMETERS:
!!  loc_proc : local processor #
!!  d11, d12, d21, d22, d31, d32, d41, d42, d51, d52 : dimensions for var_5d
!!  var_5d   : the specified 5D array
!!  loval    : the lowest  allowed value
!!  hival    : the highest allowed value
      character (len=*) :: var_name
      integer :: loc_proc
      integer :: d11, d12, d21, d22, d31, d32, d41, d42, d51, d52
      real*8  :: var_5d(d11:d12, d21:d22, d31:d32, d41:d42, d51:d52)
      real*8  :: loval, hival
!
! !DESCRIPTION
!   This routine does a gross range check on the data values in a 5D array,
!   and exits if it finds one out of range.
!   var_name : the string name of var_5d
!
! !LOCAL VARIABLES:
      character (len=75) :: err_msg

      integer :: max5(5), min5(5)
!
! !AUTHOR: 
!  John Tannahill (LLNL) and Jules Kouatchou
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC
      min5(:) = Minloc (var_5d(d11:d12,d21:d22,d31:d32,d41:d42,d51:d52))
      max5(:) = Maxloc (var_5d(d11:d12,d21:d22,d31:d32,d41:d42,d51:d52))

      min5(1) = min5(1) + d11 - 1
      min5(2) = min5(2) + d21 - 1
      min5(3) = min5(3) + d31 - 1
      min5(4) = min5(4) + d41 - 1
      min5(5) = min5(5) + d51 - 1

      max5(1) = max5(1) + d11 - 1
      max5(2) = max5(2) + d21 - 1
      max5(3) = max5(3) + d31 - 1
      max5(4) = max5(4) + d41 - 1
      max5(5) = max5(5) + d51 - 1
 
      if ((var_5d(min5(1),min5(2),min5(3),min5(4),min5(5)) < loval) .or.  &
          (var_5d(max5(1),max5(2),max5(3),max5(4),max5(5)) > hival))      &
        then
 
        Write (6,*) min5(:), loc_proc
        Write (6,*) max5(:), loc_proc
 
        err_msg = 'Value out of range in CheckRange5d, ' // var_name
 
        call GmiPrintError                                        &
          (err_msg, .true., 0, 0, 0, 2,                        &
           var_5d(min5(1),min5(2),min5(3),min5(4),min5(5)),    &
           var_5d(max5(1),max5(2),max5(3),max5(4),max5(5)))
 
      end if
 
      return

      end subroutine CheckRange5d
!EOC
!-------------------------------------------------------------------------
end module GmiCheckRange_mod
