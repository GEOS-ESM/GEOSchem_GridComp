MODULE GmiFileUnit_mod

 IMPLICIT NONE

 PRIVATE
 PUBLIC  :: GetFileUnitNumber

 INTEGER, PARAMETER :: MIN_UNIT_NUMBER =  7
 INTEGER, PARAMETER :: MAX_UNIT_NUMBER = 99

 CONTAINS

! -----------------------------------------------------------------

  SUBROUTINE GetFileUnitNumber(lun, error)

   IMPLICIT NONE

   INTEGER, INTENT(OUT) :: lun
   INTEGER, INTENT(OUT) :: error

   LOGICAL :: I_Am_OPENed, I_Am_Not_OPENed

   lun = MIN_UNIT_NUMBER
   error = 0

! Try to find an open logical unit number.
! ----------------------------------------
   Searching: DO

    INQUIRE(UNIT=lun, OPENED=I_Am_OPENed)
    
    I_Am_Not_OPENed = .NOT. I_Am_OPENed

! If an unopened unit number is found, we are done!
! -------------------------------------------------
    IF(I_Am_Not_OPENed) EXIT

! Otherwise increment and try again
! ---------------------------------
    lun = lun+1

! Did we exhaust the available unit numbers?
! ------------------------------------------
    IF(lun > MAX_UNIT_NUMBER) EXIT

   END DO Searching

! Return unopened unit number or print error message.
! ---------------------------------------------------
   IF(lun > MAX_UNIT_NUMBER) THEN
    error = -1
    PRINT *, "GetFileUnitNumber: Unable to find an available logical unit."
   END IF

   RETURN
  END SUBROUTINE GetFileUnitNumber 

! -----------------------------------------------------------------

END MODULE GmiFileUnit_mod
