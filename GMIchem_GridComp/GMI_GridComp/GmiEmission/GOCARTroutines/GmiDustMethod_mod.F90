!-------------------------------------------------------------------------
!  NASA/GFSC, SIVO, Code 610.3
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: GmiDustMethod_mod
!
! !INTERFACE:
!
module GmiDustMethod_mod
!
  implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
!
  private
  public  :: Allocate_srcEmissDust
  public  :: Allocate_erod
  public  :: Allocate_erod_mod
  public  :: Update_erod
  public  :: InitializationDust
  public  :: SourceDust
!
! !PUBLIC DATA MEMBERS:
!
  public  :: ndcls, nDustBin, ndsrc
  public  :: srcEmissDust
  public  :: erod
  public  :: erod_mod
!
  INTEGER, PARAMETER :: ndcls=3
  INTEGER, PARAMETER :: nDustBin=5       ! number of dust bins
  INTEGER, PARAMETER :: ndsrc=1
  INTEGER, PARAMETER :: msub=8

  INTEGER            :: ipoint  (nDustBin)
  REAL*8             :: dustden (nDustBin)    ! Dust density            (kg/m3)
  REAL*8             :: dustreff(nDustBin)    ! Main effective radius   (m)
  REAL*8             :: frac_s  (nDustBin)
  REAL*8             :: vdif    (nDustBin)
  REAL*8             :: ch_dust (nDustBin,12) ! Constant to fudge the total emission of dust (s2/m2)
  REAL*8 , POINTER   :: erod    (:,:,:,:) => null() ! Fraction of erodible grid cell
                                                    ! for 1: Sand, 2: Silt, 3: Clay
  REAL*8 , POINTER   :: erod_mod(:,:)     => null()
  REAL*8 , POINTER   :: srcEmissDust(:,:,:) => null() ! dust emission computed from source
!
! !DESCRIPTION:
!  Provides routines to manipulate Dust related variables.
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
! !IROUTINE: Allocate_srcEmissDust
!
! !INTERFACE:
!
  subroutine Allocate_srcEmissDust(i1, i2, ju1, j2)
!
  IMPLICIT NONE
!
! !INPUT PARAMETERS:
  integer, intent(in) :: i1, i2, ju1, j2
!
! !DESCRIPTION:
! Allocate the variable "srcEmissDust".
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
  allocate (srcEmissDust (i1:i2, ju1:j2, nDustBin))
  srcEmissDust(:,:,:) = 0.0d0
  return
  end subroutine Allocate_srcEmissDust
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Allocate_erod
!
! !INTERFACE:
!
  subroutine Allocate_erod(i1, i2, ju1, j2)
!
  IMPLICIT NONE
!
! !INPUT PARAMETERS:
  integer, intent(in) :: i1, i2, ju1, j2
!
! !DESCRIPTION:
! Allocate the variable "erod".
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
  allocate (erod (i1:i2, ju1:j2, ndcls,ndsrc))
  return
  end subroutine Allocate_erod
!EOC
!-------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Allocate_erod_mod
!
! !INTERFACE:
!
  subroutine Allocate_erod_mod(i1, i2, ju1, j2)
!
  IMPLICIT NONE
!
! !INPUT PARAMETERS:
  integer, intent(in) :: i1, i2, ju1, j2
!
! !DESCRIPTION:
! Allocate the variable "erod_mod".
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
  allocate (erod_mod (i1:i2, ju1:j2))
  return
  end subroutine Allocate_erod_mod
!EOC
!-------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: InitializationDust
!
! !INTERFACE:
!
  subroutine InitializationDust( )
!
  IMPLICIT NONE
!
! !DESCRIPTION:
! Initialize dust related variables ("dustden", "dustreff",
! "ipoint", "frac_s", and "ch_dust"), depending of the number
! of dust bins (nDustBin).
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

  IF (nDustBin == 5) THEN
     dustden(1) = 2500.0  ! Clay
     dustden(2) = 2650.0  ! Small silt
     dustden(3) = 2650.0  ! Small silt
     dustden(4) = 2650.0  ! Small silt
     dustden(5) = 2650.0  ! Small silt
  END IF
  !       dustden(6) = 2650.0  ! Large silt
  !       dustden(7) = 2650.0  ! Sand
  
  ! Main effective radius (m)
  IF (nDustBin == 5) THEN
     dustreff(1) = 0.73E-6 ! Clay
     dustreff(2) = 1.4E-6  ! Small silt
     dustreff(3) = 2.4E-6  ! Small silt
     dustreff(4) = 4.5E-6  ! Small silt
     dustreff(5) = 8.0E-6  ! Small silt
  END IF
  !       dustreff(6) = 18.0E-6  ! Large silt
  !       dustreff(7) = 38.0E-6  ! Sand     
  
  IF (nDustBin == 5) THEN
     ipoint(1) = 3
     ipoint(2) = 2
     ipoint(3) = 2
     ipoint(4) = 2
     ipoint(5) = 2
     frac_s(1) = 0.1
     frac_s(2) = 0.25
     frac_s(3) = 0.25
     frac_s(4) = 0.25
     frac_s(5) = 0.25
  ELSE IF (nDustBin ==4) THEN
     ipoint(1) = 3
     ipoint(2) = 2
     ipoint(3) = 2
     ipoint(4) = 2
     frac_s(1) = 0.1
     frac_s(2) = 0.25
     frac_s(3) = 0.25
     frac_s(4) = 0.25
  END IF

!  IF (ndsrc == 1) THEN
!     ch_dust(1) = 1.0E-9      ! Transfert coeff for type natural source 
!                              ! (kg*s2/m5)
!     ch_dust(1) = 0.75E-9 ! Pete's modified value 
!  END IF

  ! set ch_dust; 
  ! ch_dust is used in source_du to fudge the total emission of dust;
  ! ch_dust_aerocom: total monthly emissions from aerocom run
  ! ch_dust_geos4_mod: total monthly emissions from geos4 run with 
  ! ch_dust=0.75E-9 and gwet=0.5
  ! this takes care of the different magnitude of w10m in geos3 and geos4;
  ! the threshold value for gwet still has to be modified in source_du in 
  ! order to get the right amount of emission per grid box (not just 
  ! globally)
  ch_dust(:,:) = 1.0E-9
!  ch_dust(:,:) = 0.75E-9
!       0.75E-9 * ( &
!       RESHAPE(SOURCE=ch_dust_aerocom_1(:),  SHAPE=(/5,12/)) / &
!       RESHAPE(SOURCE=ch_dust_geos4_mod_1(:),SHAPE=(/5,12/)) )
     
  return
  end subroutine InitializationDust
!EOC
!-------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: SourceDust
!
! !INTERFACE:
!
  SUBROUTINE SourceDust (mcor, tdt, nymd, i1, i2, ju1, j2, k1, k2)
!
  use GocartDerivedVariables_mod, only : w10m, airden, airmas
  use GocartDerivedVariables_mod, only : gwet => gwet_1
  use GmiWaterMethod_mod        , only : water
  use GmiTimeControl_mod        , only : GmiSplitDateTime
!
  IMPLICIT NONE

# include "gmi_phys_constants.h"
!
! !INPUT PARAMETERS:
  INTEGER, INTENT(IN)    :: i1, i2, ju1, j2, k1, k2
  INTEGER, INTENT(IN)    :: nymd
  REAL*8,  INTENT(IN)    :: tdt
  REAL*8,  INTENT(IN)    :: mcor  (i1:i2, ju1:j2)         ! area of grid box (m2)
!
! !OUTPUT PARAMETERS:
!  REAL*8,  INTENT(OUT)   :: bems      (i1:i2, ju1:j2, nDustBin)
  REAL*8 bems      (i1:i2, ju1:j2, nDustBin)
                                                          ! dust tracer mixing ratio at the
                                                          ! first model layer
!
! !DESCRIPTION:
!  Evaluate the source of each dust particles size classes  (kg/m3) 
!  by soil emission.
!
! !LOCAL VARIABLES:
!  REAL*8  :: gwet  (i1:i2, ju1:j2)         ! ground wetness (fraction)
  REAL*8  :: den(nDustBin), diam(nDustBin)
  REAL*8  :: tsrc, rhoa, u_ts0, cw, u_ts, dsrc, g, srce
  INTEGER :: i, j, n, m, k
  integer :: idumday, idumyear
  integer :: month_gmi
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

  srcEmissDust(:,:,:) = 0.0d0

  call GmiSplitDateTime (nymd, idumyear, month_gmi, idumday)

  DO n = 1, nDustBin
     ! Threshold velocity as a function of the dust density and the diameter
     ! from Bagnold (1941)
     den(n) = dustden(n)*1.0E-3
     diam(n) = 2.0*dustreff(n)*1.0E2
     g = GMI_G*1.0E2
     ! Pointer to the 3 classes considered in the source data files
     m = ipoint(n)
     tsrc = 0.0
     DO k = 1, ndsrc
        ! No flux if wet soil
        DO i = i1,i2
           DO j = ju1,j2
              rhoa = airden(i,j,1)*1.0E-3
              u_ts0 = 0.13*1.0E-2*SQRT(den(n)*g*diam(n)/rhoa)* &
                   SQRT(1.0+0.006/den(n)/g/(diam(n))**2.5)/ &
                   SQRT(1.928*(1331.0*(diam(n))**1.56+0.38)**0.092-1.0)
              ! Fraction of emerged surfaces (subtract lakes, coastal ocean,..)
              cw = 1.0 - water(i,j)
              ! Case of surface dry enough to erode
              IF (gwet(i,j) < 0.2) THEN
!              IF (gwet(i,j) < 0.5) THEN  !  Pete's modified value
                 u_ts = MAX(0.0,u_ts0*(1.2+0.2*LOG10(MAX(1.0E-3, gwet(i,j)))))
              ELSE
                 ! Case of wet surface, no erosion
                 u_ts = 100.0
              END IF
              srce = frac_s(n)*erod(i,j,m,k)*mcor(i,j)  ! (m2)
              dsrc = cw*ch_dust(n,month_gmi)*srce*w10m(i,j)**2 &
                   * (w10m(i,j) - u_ts)  ! (kg)
!                   * (w10m(i,j) - u_ts)*tdt  ! (kg)
!              if (erod(i,j,m,k) .ne. 0.0d0) then
!              end if
              IF (dsrc < 0.0) dsrc = 0.0

              ! Update dust mixing ratio at first model level.
              srcEmissDust(i,j,n) = srcEmissDust(i,j,n) + dsrc 
!              srcEmissDust(i,j,n) = srcEmissDust(i,j,n) + dsrc / airmas(i,j,1)

              
              bems(i,j,n) = dsrc
           END DO
        END DO
     END DO
  END DO

  return

  END SUBROUTINE SourceDust
!EOC
!-------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Update_erod
!
! !INTERFACE:
!
  subroutine Update_erod (latdeg, londeg, mcor, &
                          i1, i2, ju1, j2, i1_gl, i2_gl, ju1_gl, j2_gl)
!
  IMPLICIT none
!
! !INPUT PARAMETERS:
  integer, intent(in) :: i1, i2, ju1, j2, i1_gl, i2_gl, ju1_gl, j2_gl
  real*8 , intent(in) :: londeg(i1_gl:i2_gl)
  real*8 , intent(in) :: latdeg(ju1_gl:j2_gl)
  real*8 , intent(in) :: mcor  (i1:i2, ju1:j2)
!
! !DESCRIPTION:
! Adjust the variable erod and convert it from area to fraction.
! erod: 2D soil area per grid cell of sand, silt and clay which can erode
! erod_mod:   replace Asia source with Mian's modification
!
! !LOCAL VARIABLES:
  integer             :: i, j, k, m, n
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
!
  DO j = ju1, j2
!     IF (30.0 <= latdeg(j) .AND. latdeg(j) <= 50.0) THEN
        DO i = i1, i2
           !        IF ((30.0 <= latdeg(j) .AND. latdeg(j) <= 50.0) .AND. &
           !            (77.5 <= londeg(i) .AND. londeg(i) <= 130.0)) THEN
!        IF (i_erod_mod(i) == 1) THEN
           IF (erod_mod(i,j) > 0.0) THEN
              DO k = 1,ndcls
                 IF (k == 1) then
                    !! erod(104:125,j,k,1) = erod_mod(1:22)
                     erod(i,j,k,1) = erod_mod(i,j)
                 ELSE
                    !! erod(104:125,j,k,1) = erod_mod(1:22)/2.0 *1.5
                     erod(i,j,k,1) = erod_mod(i,j)/2.0 *1.5
                 END IF
              END DO
           END IF
        END DO
!     END IF
  END DO


  ! Reduce the Africa/Middle East sources
  DO j = ju1, j2
     DO i = i1, i2
        IF ((  0.0 <= latdeg(j) .AND. latdeg(j) <= 50.0) .AND. &
            (( 0.0 <= londeg(i) .AND. londeg(i) <= 60.0) .OR.  &
             (340.0 <= londeg(i) .AND. londeg(i) <= 360.0))) THEN
!             (-20.0 <= londeg(i) .AND. londeg(i) <= 60.0)) THEN
           erod(i,j,:,:) = erod(i,j,:,:)*0.8
        END IF
     END DO
  END DO

  ! convert erod from area to fraction
  DO n = 1,ndsrc
     DO m = 1,ndcls
        DO j = ju1, j2
           DO i = i1, i2
              erod(i,j,m,n) = erod(i,j,m,n)/mcor(i,j)
           END DO
        END DO
     END DO
  END DO

  return
  
  end subroutine Update_erod
!EOC
!-------------------------------------------------------------------------------
end module GmiDustMethod_mod
