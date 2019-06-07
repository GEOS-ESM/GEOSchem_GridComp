!=============================================================================
!
! $Id$
!
! CODE DEVELOPER
!   originated from GRANTOUR programmed by John Walton
!   modified by Xiaohong Liu
!
! FILE
!   const_sulf_update.F
!
! ROUTINES
!   Do_Const_Sulf_Update
!
!=============================================================================

!-----------------------------------------------------------------------------
!
! ROUTINE
!   Do_Const_Sulf_Update (chem40)
!
! DESCRIPTION
!
!     *****************************************************************
!     * advances the species concentrations and does mixing of clear
!     *   and cloudy regions
!     *****************************************************************
!
!     ********** input:
!
!     *     cp     : concentration (# cm-3 for gas-species, g/g for aerosols)
!     *     dcp    : storage for conc. change
!
!     *     fso4a  : aqueous fossil H2SO4 conc.
!     *     fso4n  : clear sky fossil H2SO4 conc.
!     *     fso2   : fossil SO2 conc.
!
!     *     nso4a  : aqueous natural H2SO4 conc.
!     *     nso4n  : clear sky natural H2SO4 conc.
!     *     nso2   : natural SO2 conc.
!     *     ndms   : natural DMS conc.
!
!     *     h2o2   : H2O2 conc.
!
!     *     dfso4a : aqueous fossil H2SO4 conc.
!     *     dfso4n : clear sky fossil H2SO4 conc.
!     *     dfso2  : fossil SO2 conc.
!
!     *     dnso4a : aqueous natural H2SO4 conc.
!     *     dnso4n : clear sky natural H2SO4 conc.
!     *     dnso2  : natural SO2 conc.
!     *     dndms  : natural DMS conc.
!
!     *     dh2o2  : H2O2 conc.
!
!
!     * For H2SO4 the aqueous and clear sky parts remain seperate.
!       However, the H2O2, DMS and SO2 are a combination of the aqueous
!       and clear sky parts, weighted by the cloud fraction cfp.
!
!     * The clear sky changes were computed and stored in chem11.
!       In chem30 the aqueous changes were computed and,
!       except for H2SO4, added to the clear sky.  e.g., with cfp
!       the cloud fraction at a grid point
!
!       dcp(ns) = (1-cfp) * d[clear sky] + cfp * d[aqueous]
!
!----------------------------------------------------------------------

      subroutine Do_Const_Sulf_Update  &
     &  (itloop, cp, dcp)

      implicit none

#     include "setkin_par.h"
#     include "sulfchem.h"

!     ----------------------
!     GMIMOD-related values.
!     ----------------------

      integer :: itloop

      real*8  :: cp (itloop, NSP)
      real*8  :: dcp(itloop, NSP)


      cp(:, IH2O2)   = cp(:, IH2O2) + dcp(:,IH2O2)
      cp(:, INDMS)   = cp(:, INDMS) + dcp(:,INDMS)
      cp(:, INSO2)   = cp(:, INSO2) + dcp(:,INSO2)
      cp(:, IFSO2)   = cp(:, IFSO2) + dcp(:,IFSO2)

      cp(:, IFSO4A)  = cp(:, IFSO4A ) + dcp(:,IFSO4A)
      cp(:, IFSO4N1) = cp(:, IFSO4N1) + dcp(:,IFSO4N1)
      cp(:, IFSO4N2) = cp(:, IFSO4N2) + dcp(:,IFSO4N2)
      cp(:, IFSO4N3) = cp(:, IFSO4N3) + dcp(:,IFSO4N3)

      cp(:, INSO4A)  = cp(:, INSO4A ) + dcp(:,INSO4A)
      cp(:, INSO4N1) = cp(:, INSO4N1) + dcp(:,INSO4N1)
      cp(:, INSO4N2) = cp(:, INSO4N2) + dcp(:,INSO4N2)
      cp(:, INSO4N3) = cp(:, INSO4N3) + dcp(:,INSO4N3)

      where(cp(:, IH2O2)    < 1.0d-30) cp(:, IH2O2) = 1.0d-30
      where(cp(:, INDMS)    < 1.0d-30) cp(:, INDMS) = 1.0d-30
      where(cp(:, INSO2)    < 1.0d-30) cp(:, INSO2) = 1.0d-30
      where(cp(:, IFSO2)    < 1.0d-30) cp(:, IFSO2) = 1.0d-30

      where(cp(:, IFSO4A)   < 1.0d-30) cp(:, IFSO4A ) = 1.0d-30
      where(cp(:, IFSO4N1)  < 1.0d-30) cp(:, IFSO4N1) = 1.0d-30
      where(cp(:, IFSO4N2)  < 1.0d-30) cp(:, IFSO4N2) = 1.0d-30
      where(cp(:, IFSO4N3)  < 1.0d-30) cp(:, IFSO4N3) = 1.0d-30

      where(cp(:, INSO4A)   < 1.0d-30) cp(:, INSO4A ) = 1.0d-30
      where(cp(:, INSO4N1)  < 1.0d-30) cp(:, INSO4N1) = 1.0d-30
      where(cp(:, INSO4N2)  < 1.0d-30) cp(:, INSO4N2) = 1.0d-30
      where(cp(:, INSO4N3)  < 1.0d-30) cp(:, INSO4N3) = 1.0d-30


      return

      end

