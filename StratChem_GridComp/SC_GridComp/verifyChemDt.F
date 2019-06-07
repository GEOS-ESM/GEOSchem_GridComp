      SUBROUTINE verifyChemDt(examineDt,requiredPasses,chemDt,dt,requiredDt,startOver)
! --------------------------------------------------------------------------
! Purpose: SOLVEREST and/or SOLVERD2 return the length of the time step 
! that is required to keep the routines stable.  If that dt is less than
! the original dt, then this routine establishes a new chemistry dt, 
! determines the number of iterations, and sets a switch that tells chemdrv
! to restart the chemistry at the current grid box.
! --------------------------------------------------------------------------
      IMPLICIT NONE

      LOGICAL, INTENT(INOUT) :: startOver
      INTEGER, INTENT(INOUT) :: examineDt,requiredPasses
      REAL, INTENT(IN) :: requiredDt,chemDt
      REAL, INTENT(INOUT) :: dt

      REAL, PARAMETER :: tol = 1.00E-07

      IF(requiredDt < dt) THEN
       requiredPasses = INT(chemDt/(requiredDt+tol))+1
       dt = chemDt/FLOAT(requiredPasses)
       startOver = .TRUE.
      END IF

      examineDt = examineDt+1

      RETURN
      END SUBROUTINE verifyChemDt
