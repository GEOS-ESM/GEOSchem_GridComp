
!=============================================================================
!
! $Id$
!
! CODE DEVELOPER
!   John Tannahill, LLNL (Original code from Keith Grant, LLNL)
!   jrt@llnl.gov
!
! FILE
!   solrza.F
!
! ROUTINES
!   Solrza
!
!=============================================================================


!-----------------------------------------------------------------------------
!
! ROUTINE
!   Solrza
!
! DESCRIPTION
!   Given a Greenwich time, a solar declination and reference latitude and
!   longitudes,  this routine returns a corresponding list of cosines of solar
!   zenith angles.  Refer to:  Paltridge, G.W., and C.M.R. Platt, 1976:
!   "Radiative Processes in Meteorology and Climatology", Elsevier, p. 62.
!
! ARGUMENTS
!   time   : Greenwich time since Jan 1, counting zero from midnight (days)
!   decl   : solar declination (deg)
!   lat    : latitude   (deg)
!   lon    : longitudes (deg)
!   nn     : number of longitudes for which to calculate cosines of the
!            solar zenith angle
!   cossza : cosines of the solar zenith angle (output)
!
!-----------------------------------------------------------------------------

      subroutine Solrza  &
     &  (time, decl, lat, lon, nn, cossza)

      implicit none

#     include "gmi_phys_constants.h"


!     ----------------------
!     Argument declarations.
!     ----------------------

      integer :: nn
      real*8  :: time
      real*8  :: decl
      real*8  :: lat
      real*8  :: lon   (nn)
      real*8  :: cossza(nn)


!     -----------------------
!     Parameter declarations.
!     -----------------------

      real*8,  parameter :: TWO_PI = 2.0d0 * GMI_PI
      real*8,  parameter :: PI_180 = TWO_PI / 360.0d0


!     ----------------------
!     Variable declarations.
!     ----------------------

      integer :: ii

      real*8  :: cosha


!     ----------------
!     Begin execution.
!     ----------------

      do ii = 1, nn

!       -------------------------------------------------------------
!       Calculate the cosine of the hour angle which is referenced to
!       noon.
!       -------------------------------------------------------------

        cosha =  &
     &    Cos (TWO_PI *  &
     &         (Mod (time + (lon(ii) / 360.d0), 1.0d0) - 0.5d0))


!       -----------------------------------------------
!       Calculate the cosine of the solar zenith angle.
!       -----------------------------------------------

        cossza(ii) =  &
     &    Sin (PI_180 * lat) * Sin (PI_180 * decl) +  &
     &    Cos (PI_180 * lat) * Cos (PI_180 * decl) * cosha

      end do


      return

      end

