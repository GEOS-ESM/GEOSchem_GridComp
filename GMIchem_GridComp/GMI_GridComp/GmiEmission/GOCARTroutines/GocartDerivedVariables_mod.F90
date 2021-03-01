!-------------------------------------------------------------------------
!  NASA/GFSC, SIVO, Code 610.3
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: GocartDerivedVariables_mod
!
! !INTERFACE:
!
module GocartDerivedVariables_mod
!
  implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
!
  private
  public  :: w10m
  public  :: airden
  public  :: airmas
  public  :: gwet_1
  public  :: lwi_flags_1
  public  :: AllocateGocartDerivedVars
  public  :: SetGocartDerivedVars

  REAL*8 , pointer :: w10m  (:,:)      => null() ! vertical wind at 10m (m/s)
  REAL*8 , pointer :: gwet_1(:,:)      => null() ! ground wetness (fraction)
  REAL*8 , pointer :: airden(:,:,:)    => null() ! air density (kg/m3)
  REAL*8 , pointer :: airmas(:,:,:)    => null() ! air mass    (kg)
  INTEGER, pointer :: lwi_flags_1(:,:) => null() ! land, water, ice indicator
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
! !IROUTINE: AllocateGocartDerivedVars
!
! !INTERFACE:
!
  subroutine AllocateGocartDerivedVars(i1, i2, ju1, j2, k1, k2)
!
  implicit none
!
! !INPUT PARAMETERS:
  integer, intent(in) :: i1, i2, ju1, j2, k1, k2
!
! !DESCRIPTION:
! Allocate the variables 
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
  allocate (w10m   (i1:i2, ju1:j2))
  allocate (gwet_1 (i1:i2, ju1:j2))
  allocate (lwi_flags_1 (i1:i2, ju1:j2))
  allocate (airden (i1:i2, ju1:j2, k1:k2))
  allocate (airmas (i1:i2, ju1:j2, k1:k2))
  return
  end subroutine AllocateGocartDerivedVars
!EOC
!-------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: SetGocartDerivedVars
!
! !INTERFACE:
!
  subroutine SetGocartDerivedVars (u10m, v10m, gwet, press3c, &
                            kel, mass, mcor, grid_height, lwi_flags, &
                            i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi)
!
  implicit none

# include "gmi_phys_constants.h"

!
! !INPUT PARAMETERS:
  integer, intent(in) :: i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi
  real*8 , intent(in) :: u10m       (i1:i2, ju1:j2)
  real*8 , intent(in) :: v10m       (i1:i2, ju1:j2)
  real*8 , intent(in) :: gwet       (i1:i2, ju1:j2)
  integer, intent(in) :: lwi_flags  (i1:i2, ju1:j2)
  real*8 , intent(in) :: mcor       (i1:i2, ju1:j2) ! area of grid box [m^2]
  real*8 , intent(in) :: mass       (i1:i2, ju1:j2, k1:k2)
  real*8 , intent(in) :: grid_height(i1:i2, ju1:j2, k1:k2)
  real*8 , intent(in) :: kel        (ilo:ihi, julo:jhi, k1:k2)
  real*8 , intent(in) :: press3c    (ilo:ihi, julo:jhi, k1:k2)
!
! !DESCRIPTION:
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

  integer :: ivert

  ivert = k2-k1+1

  w10m(:,:)        = sqrt(u10m(:,:)**2 + v10m(:,:)**2)
  gwet_1(:,:)      = gwet(:,:)
  lwi_flags_1(:,:) = lwi_flags(:,:)
  airmas(:,:,:)    = mass(:,:,:)
!  airmas(:,:,:)    = press3c(i1:i2,ju1:j2,:) * MB2CGS / &
!                      (kel (i1:i2,ju1:j2,:) * BOLTZMN_E)
  airden(:,:,:)    = mass(:,:,:) / (spread(mcor(:,:),3,ivert) * &
                          grid_height(:,:,:))

  return
  
  end subroutine SetGocartDerivedVars
!EOC
!-------------------------------------------------------------------------------
end module GocartDerivedVariables_mod
