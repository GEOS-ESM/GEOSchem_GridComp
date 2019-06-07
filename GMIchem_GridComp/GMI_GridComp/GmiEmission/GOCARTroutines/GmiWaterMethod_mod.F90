!-------------------------------------------------------------------------
!  NASA/GFSC, SIVO, Code 610.3
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: GmiWaterMethod_mod
!
! !INTERFACE:
!
module GmiWaterMethod_mod
!
  implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
!
  private
  public  :: water
  public  :: Allocate_water
  public  :: Update_water

  REAL*8 , ALLOCATABLE :: water(:,:)
!
! !DESCRIPTION:
!  Provides routines to manipulate the variable "water".
!
! !AUTHOR: 
!  Jules Kouatchou
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------------
CONTAINS
!-------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Allocate_water
!
! !INTERFACE:
!
  subroutine Allocate_water(i1, i2, ju1, j2)
!
  implicit none
!
! !INPUT PARAMETERS:
  integer, intent(in) :: i1, i2, ju1, j2
!
! !DESCRIPTION:
! Allocate the variable "water".
!
! !AUTHOR: 
!  Jules Kouatchou
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC
  allocate (water (i1:i2, ju1:j2))
  return
  end subroutine Allocate_water
!EOC
!-------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Update_water
!
! !INTERFACE:
!
  subroutine Update_water (mcor, i1, i2, ju1, j2)
!
  implicit none
!
! !INPUT PARAMETERS:
  integer, intent(in) :: i1, i2, ju1, j2
  real*8 , intent(in) :: mcor(i1:i2, ju1:j2) ! area of grid box [m^2]
!
! !DESCRIPTION:
!  Ajusts the variable "water"
!  water: area of water surface per grid box (from 10'x10' Navy database)
!  water: fraction of water per grid box (from Paul's original file)
!
! !AUTHOR: 
!  Jules Kouatchou
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC

     ! convert water from area to fraction
!     DO j = ju1, j2
!        DO i = i1, i2
!           water(i,j) = water(i,j)/mcor(i,j)
!        END DO
!     END DO

  return
  
  end subroutine Update_water
!EOC
!-------------------------------------------------------------------------------
end module GmiWaterMethod_mod
