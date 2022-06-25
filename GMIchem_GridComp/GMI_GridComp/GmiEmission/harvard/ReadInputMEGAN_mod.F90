!-------------------------------------------------------------------------
!  NASA/GFSC, SIVO, Code 610.3
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: ReadInputMEGAN_mod
!
! !INTERFACE:
!
  MODULE ReadInputMEGAN_mod

! !USES:
  USE GmiTimeControl_mod, ONLY : GmiSplitDateTime

  IMPLICIT NONE

! !PUBLIC MEMBER FUNCTIONS:

  PRIVATE
  PUBLIC  :: setMEGANisoLAI

#      include "gmi_time_constants.h"

! !DESCRIPTION:

! !AUTHOR:
!  Bob Yantosca
!  Jules Kouatchou

! !REVISION HISTORY:
!  Initial code.
!  20 Nov 2008 Eric Nielsen For use in GEOS-5.
!
!EOP
!------------------------------------------------------------------------------
     contains
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: setMEGANisoLAI
!
! !INTERFACE:
!
 SUBROUTINE setMEGANisoLAI(isoLai, isoLaiCurr, isoLaiPrev, isoLaiNext, &
                            days_btw_m, nymd, i1, i2, ju1, j2, i1_gl, ju1_gl)

  IMPLICIT NONE

! !INPUT PARAMETERS:
  INTEGER, INTENT(IN) :: i1, i2, ju1, j2, i1_gl, ju1_gl
  INTEGER, INTENT(IN) :: nymd
!
! !OUTPUT PARAMETERS:

  INTEGER, INTENT(OUT) :: days_btw_m ! days between midmonths

! !INPUT/OUTPUT PARAMETERS:

  REAL*8 , INTENT(INOUT) ::  isoLai(i1:i2,ju1:j2) ! Daily LAI data (interpolated)
  REAL*8 , INTENT(IN) :: isoLaiCurr(i1:i2,ju1:j2) !current  month
  REAL*8 , INTENT(IN) :: isoLaiPrev(i1:i2,ju1:j2) !previous month
  REAL*8 , INTENT(IN) :: isoLaiNext(i1:i2,ju1:j2) !next     month

! !DESCRIPTION: 
!  Sets isoLai daily.  The stored monthly LAI are used for the middle day in the
!  month and LAIs are interpolated for other days. (dsa, tmf, bmy, 10/20/05)

! !LOCAL VARIABLES:
  INTEGER :: thisMonth, thisYear, thisDayOfMonth, daysThisMonth
  INTEGER :: previousMonth, daysInFeb
  REAL :: fraction

! !DEFINED PARAMETERS:
  INTEGER, PARAMETER :: monthLength(MONTHS_PER_YEAR) = &
                           (/  31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)
!
! !AUTHOR:
!  Eric Nielsen
!
! !REVISION HISTORY:
!  20 Nov 2008 Eric Nielsen For use in GEOS-5.
!
!EOP
!------------------------------------------------------------------------------

! Extract year, month, and day
! ----------------------------
  CALL GmiSplitDateTime(nymd, thisYear, thisMonth, thisDayOfMonth)

! How many days in the current month?
! -----------------------------------
  IF(thisMonth == 2) THEN
   daysThisMonth = 28
   IF(MOD(thisYear,  4) == 0) daysThisMonth = 29
   IF(MOD(thisYear,100) == 0) daysThisMonth = 28
   IF(MOD(thisYear,400) == 0) daysThisMonth = 29
  ELSE
   daysThisMonth = monthLength(thisMonth)
  END IF

! Linearly interpolate using past month's LAI through
! the 15th and next month's LAI until month's end.
! ---------------------------------------------------
  IF(thisDayOfMonth < 16) THEN
   fraction = 0.50*(1.00+(thisDayOfMonth-1.00)/15.00)
   isoLai(:,:) = (1.00-fraction)*isoLaiPrev(:,:) + fraction*isoLaiCurr(:,:)
  ELSE
   fraction = 0.50*(1.00+(daysThisMonth+1.00-thisDayOfMonth)/(daysThisMonth-15.00))
   isoLai(:,:) = fraction*isoLaiCurr(:,:) + (1.00-fraction)*isoLaiNext(:,:)
  END IF

! Number of days between current and previous midmonth
! ----------------------------------------------------
  previousMonth = thisMonth-1
  IF(previousMonth == 0) previousMonth = 12

  IF(previousMonth == 2) THEN
   daysInFeb = 28
   IF(MOD(thisYear,  4) == 0) daysInFeb = 29
   IF(MOD(thisYear,100) == 0) daysInFeb = 28
   IF(MOD(thisYear,400) == 0) daysInFeb = 29
   days_btw_m = (daysInFeb+daysThisMonth)*0.50
  ELSE
   days_btw_m = (monthLength(previousMonth)+daysThisMonth)*0.50
  END IF

  RETURN

 END SUBROUTINE setMEGANisoLAI

END MODULE ReadInputMEGAN_mod
