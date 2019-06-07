
!=============================================================================
!
! $Id$
!
! CODE DEVELOPER
!   John Tannahill, LLNL (Original code from Keith Grant, LLNL)
!   jrt@llnl.gov
!
! FILE
!   solrpos.F
!
! ROUTINES
!   Solrpos
!
!=============================================================================


!-----------------------------------------------------------------------------
!
! ROUTINE
!   Solrpos
!
! DESCRIPTION
!   Given the Julian day, this routine calculates the solar declination and
!   the square of the ratio of the mean Earth-Sun distance to the current
!   Earth-Sun distance.   Refer to:  Paltridge, G.W., and C.M.R. Platt, 1976:
!   "Radiative Processes in Meteorology and Climatology", Elsevier, pp. 57-63.
!
! ARGUMENTS
!   julday  : Julian day counting ranging from 0 on Jan. 1st to 364 on
!             Dec. 31st
!   decl    : solar declination (deg)
!   rdistsq : the square of the ratio of the mean Earth-Sun distance
!             to the current Earth-Sun distance
!
!-----------------------------------------------------------------------------

      subroutine Solrpos  &
     &  (julday, decl, rdistsq)

      implicit none

#     include "gmi_phys_constants.h"


!     ----------------------
!     Argument declarations.
!     ----------------------

      real*8  :: julday
      real*8  :: decl
      real*8  :: rdistsq


!     -----------------------
!     Parameter declarations.
!     -----------------------

      real*8, parameter :: TWO_PI = 2.0d0 * GMI_PI


!     ----------------------
!     Variable declarations.
!     ----------------------

      real*8  :: theta           ! year angle on specified day
      real*8  :: theta2, theta3


!     ----------------
!     Begin execution.
!     ----------------

!c?   Leap year?
      theta  = (TWO_PI * julday) / 365.0d0

      theta2 = 2.0d0 * theta
      theta3 = 3.0d0 * theta


      decl =  &
     &  (360.0d0 / TWO_PI) *  &
     &  (0.006918d0 -  &
     &   0.399912d0 * Cos (theta)  + 0.070257d0 * Sin (theta)  -  &
     &   0.006758d0 * Cos (theta2) + 0.000907d0 * Sin (theta2) -  &
     &   0.002697d0 * Cos (theta3) + 0.001480d0 * Sin (theta3))


      rdistsq =  &
     &  1.000110d0  +  &
     &  0.034221d0  * Cos (theta)  + 0.001280d0 * Sin (theta)  +  &
     &  0.000719d0  * Cos (theta2) + 0.000077d0 * Sin (theta2)


      return

      end

