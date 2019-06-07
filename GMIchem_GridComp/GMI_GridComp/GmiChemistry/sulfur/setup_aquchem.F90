!=============================================================================
!
! $Id$
!
! CODE DEVELOPER
!   originated from GRANTOUR programmed by John Walton
!   modified by Xiaohong Liu
!
! FILE
!   setup_aquchem.F
!
! ROUTINES
!   Do_Setup_Aquchem
!
!=============================================================================

!-----------------------------------------------------------------------------
!
! ROUTINE
!   Do_Setup_Aquchem
!
! DESCRIPTION
!   This is the routine set up for aqueous chemistry by
!     creating a subset of the grid points with cfp > 0.01 and cwp > 1.e-9 g/g
!
! ARGUMENTS
!   itloop     : # of zones (ilong * ilat * ivert)
!   numcl      : number of grid points in cloud
!   idxcl      : array of the numbers of the numcl points in cloud
!   cfp        : total cloud fraction [0 - 1]
!   cwp        : in-cloud cloud water mixing ratio (g/g air)
!
! ---------------------------------------------------------------------------

      subroutine Do_Setup_Aquchem  &
     &  (itloop, numcl, idxcl, cfp, cwp)

      implicit none

#     include "sulfchem.h"

!     ----------------------
!     GMIMOD-related values.
!     ----------------------

      integer :: itloop, numcl
      integer :: idxcl(itloop)

      real*8  :: cfp(itloop)
      real*8  :: cwp(itloop)


      numcl = 0

      do 100 ijk = 1, itloop

        if (cfp(ijk) .gt. 0.01 .and. cwp(ijk) .gt. 1.e-9) then
          idxcl(numcl+1) = ijk
          numcl = numcl + 1
        else
          cfp(ijk) = 0.
          cwp(ijk) = 0.
        endif

  100 continue


      return

      end

