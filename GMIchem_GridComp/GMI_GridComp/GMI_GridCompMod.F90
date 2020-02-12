#include "MAPL_Generic.h"
!-------------------------------------------------------------------------
!NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1, GEOS/DAS!
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  GMI_GridCompMod --- GMI Grid Component Class
!
! Aggregated Grid Component class for the GMI combined stratopshere/troposphere 
! chemistry
!
! !INTERFACE:
!
      MODULE  GMI_GridCompMod

! !USES:

      USE ESMF
      USE MAPL
      USE Chem_Mod 	     ! Chemistry Base Class
      USE Chem_UtilMod

      USE GmiSAD_GCCMod
      USE GmiChem_GCCMod
      USE GmiDepos_GCCMod
      USE GmiEmiss_GCCMod
      USE GmiThermalRC_GCCMod
      USE GmiForcingBC_GCCMod
      USE GmiPhotolysis_GCCMod

      IMPLICIT NONE
!
! !DEFINED PARAMETERS:
      INTEGER, PARAMETER :: DBL = KIND(0.00D+00)

      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:

      PUBLIC  GMI_GridComp       ! The GMI object 
      PUBLIC  GMI_GridCompInitialize
      PUBLIC  GMI_GridCompRun1
      PUBLIC  GMI_GridCompRun2
      PUBLIC  GMI_GridCompRunOrig
      PUBLIC  GMI_GridCompFinalize
!
! !PUBLIC TYPES:
      TYPE GMI_GridComp
         CHARACTER(LEN=255) :: name = "GMI Stratospheric/Tropospheric Chemistry"

         ! Component derived type declarations
         ! -----------------------------------
         TYPE(GmiDepos_GridComp)      :: gcDepos 
         TYPE(GmiEmiss_GridComp)      :: gcEmiss
         TYPE(GmiSAD_GridComp)        :: gcSAD
         TYPE(GmiThermalRC_GridComp)  :: gcThermalRC
         TYPE(GmiChemistry_GridComp)  :: gcChem
         TYPE(GmiForcingBC_GridComp)  :: gcFBC
         TYPE(GmiPhotolysis_GridComp) :: gcPhot
      END TYPE GMI_GridComp
!
! !DESCRIPTION:
!  This module implements the GMI combined stratopshere/troposphere
!  chemistry. The methods (Initialize, Run and Finalize) of the
!  the following grid component classes are called here:
!  \begin{enumerate}
!  \item Deposition ({\em GmiDepos\_GridCompInitialize}, GmiDepos\_GridCompRun, GmiDepos\_GridCompFinalize)
!  \item Emission ({\em GmiEmiss\_GridCompInitialize}, GmiEmiss\_GridCompRun, GmiEmiss\_GridCompFinalize)
!  \item Surface Area Densities for Aerosols ({\em GmiSAD\_GridCompInitialize}, GmiSAD\_GridCompRun, GmiSAD\_GridCompFinalize)
!  \item Photolysis ({\em GmiPhotolysis\_GridCompInitialize}, GmiPhotolysis\_GridCompRun, GmiPhotolysis\_GridCompFinalize)
!  \item Forcing Boundary Conditions ({\em GmiForcingBC\_GridCompInitialize}, GmiForcingBC\_GridCompRun, GmiForcingBC\_GridCompFinalize)
!  \item Thermal Rate Constants ({\em GmiThermalRC\_GridCompInitialize}, GmiThermalRC\_GridCompRun, GmiThermalRC\_GridCompFinalize)
!  \item Chemistry Solver ({\em GmiChemistry\_GridCompInitialize}, GmiChemistry\_GridCompRun, GmiChemistry\_GridCompFinalize)
!  \end{enumerate}
!
! !REVISION HISTORY:
!
!  16Sep2003 da Silva  First crack.
!  24Jan2O05 Nielsen   Implementation of Code 916 chemistry
!  19Dec2005 d Silva   Minor portability mods.
!  30Oct2007 Nielsen   GMI Combo set up
!  09Sep2010 Kouatchou Added all the individuals compent class methods.
!
!EOP
!-------------------------------------------------------------------------
      CONTAINS
!-------------------------------------------------------------------------
!NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1, GEOS/DAS!
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  GMI_GridCompInitialize --- Initialize GMI_GridComp
!
! !INTERFACE:
!
   SUBROUTINE GMI_GridCompInitialize(gcGMI, w_c, impChem, expChem, nymd, nhms, cdt, gc, clock, rc)
!
   IMPLICIT none

! !INPUT PARAMETERS:

   TYPE(Chem_Bundle), INTENT(in) :: w_c         ! Chemical tracer fields, delp, +
   INTEGER,	      INTENT(in) :: nymd, nhms  ! Time from AGCM
   REAL, 	      INTENT(in) :: cdt         ! Chemistry time step (secs)

! !INPUT/OUTPUT PARAMETERS:

   TYPE(GMI_GridComp), INTENT(inOut)  :: gcGMI	    ! Grid Component
   TYPE(ESMF_State),   INTENT(inOut)  :: impChem    ! Import State
   TYPE(ESMF_State),   INTENT(inOut)  :: expChem    ! Export State
   type(ESMF_GridComp), intent(inout) :: gc      ! Grid Component
   type(ESMF_Clock),    intent(inout) :: clock   ! The clock
!
! !OUTPUT PARAMETERS:
   INTEGER, INTENT(out) ::  rc        ! Error return code:
                                      !  0 - all is well
                                      !  1 - 

! !DESCRIPTION: Initializes the GMI Grid Component. It primarily sets
!               the import state.
!
! !DEFINED PARAMETERS:
   INTEGER, PARAMETER :: DBL = KIND(0.00D+00)
   CHARACTER(LEN=*), PARAMETER :: IAm    = 'GMI_GridCompInitialize'
!
! !LOCAL VARIABLES:
   INTEGER :: ios, m, n, STATUS, procID
!
! !REVISION HISTORY:
!
!  18Sep2003 da Silva  First crack.
!  30Jun2007 Nielsen   GMI Combo set up
!   7Apr2016 Nielsen   Initialize rc
!
!EOP
!-------------------------------------------------------------------------
!BOC
      rc = 0

      CALL GmiEmiss_GridCompInitialize     (gcGMI%gcEmiss,     w_c, impChem, expChem, nymd, nhms, cdt, gc, clock, __RC__)

      CALL GmiDepos_GridCompInitialize     (gcGMI%gcDepos,     w_c, impChem, expChem, nymd, nhms, cdt,            __RC__)

      CALL GmiSAD_GridCompInitialize       (gcGMI%gcSAD,       w_c, impChem, expChem, nymd, nhms, cdt,            __RC__)

      CALL GmiPhotolysis_GridCompInitialize(gcGMI%gcPhot,      w_c, impChem, expChem, nymd, nhms, cdt,            __RC__)

      CALL GmiForcingBC_GridCompInitialize (gcGMI%gcFBC,       w_c, impChem, expChem, nymd, nhms, cdt,            __RC__)

      CALL GmiThermalRC_GridCompInitialize (gcGMI%gcThermalRC, w_c, impChem, expChem, nymd, nhms, cdt,            __RC__)

      CALL GmiChemistry_GridCompInitialize (gcGMI%gcChem,      w_c, impChem, expChem, nymd, nhms, cdt,            __RC__)

      CALL GmiEmiss_initSurfEmissBundle    (gcGMI%gcEmiss,     w_c,          expChem,                             __RC__)

      RETURN

      END SUBROUTINE GMI_GridCompInitialize
!EOC
!-------------------------------------------------------------------------
!NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1, GEOS/DAS!
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  GMI_GridCompRun1 --- The GMI Run 1 Driver
!
! !INTERFACE:
!
   SUBROUTINE GMI_GridCompRun1(gcGMI, w_c, impChem, expChem, nymd, nhms, tdt, clock, rc)

! !USES:

   IMPLICIT none

! !INPUT/OUTPUT PARAMETERS:

   TYPE(GMI_GridComp), INTENT(inOut) :: gcGMI   ! Grid Component
   TYPE(Chem_Bundle),  INTENT(inOut) :: w_c     ! Chemical tracer fields   
   TYPE(ESMF_State),   INTENT(inOut) :: impChem ! Import State
   TYPE(ESMF_State),   INTENT(inOut) :: expChem ! Export State
   TYPE(ESMF_Clock),   INTENT(inOut) :: clock   ! The clock

! !INPUT PARAMETERS:

   INTEGER, INTENT(in) :: nymd, nhms	      ! time
   REAL,    INTENT(in) :: tdt		      ! timestep (secs)

! !OUTPUT PARAMETERS:
   INTEGER, INTENT(out) ::  rc                ! Error return code:
                                              !  0 - all is well
                                              !  1 -

! !DESCRIPTION: This routine implements the GMI Strat/Trop Driver. That 
!               is, adds chemical tendencies to each of the constituents
!               Phase 1:  emissions  (done at the heartbeat)
!
! !IMPLEMENTATION NOTES:
!
!  No pointer is reservered in the export state for deposition of water.
!
! !DEFINED PARAMETERS:
      CHARACTER(LEN=*), PARAMETER :: IAm    = 'GMI_GridCompRun1'
!
! !LOCAL VARIABLES:
      INTEGER :: STATUS
      LOGICAL :: mixPBL  ! whether to explicitly distribute
                         ! aerosol emissions within the PBL
!
! !REVISION HISTORY:
!
!  18Sep2003 da Silva  First crack.
!  30Jun2007 Nielsen   GMI Combo set up
!   7Apr2016 Nielsen   Initialize rc
!
!EOP
!-------------------------------------------------------------------------
!BOC
      rc = 0

      mixPBL = .FALSE.

      CALL GmiEmiss_GridCompRun     (gcGMI%gcEmiss, w_c, impChem, expChem, nymd, nhms, tdt, clock, mixPBL, __RC__)

      CALL GmiForcingBC_GridCompRun (gcGMI%gcFBC,   w_c, impChem, expChem, nymd, nhms, tdt,                __RC__)

      RETURN

      END SUBROUTINE GMI_GridCompRun1

!EOC
!-------------------------------------------------------------------------
!NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1, GEOS/DAS!
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  GMI_GridCompRun2 --- The GMI Run 2 Driver
!
! !INTERFACE:
!
   SUBROUTINE GMI_GridCompRun2(gcGMI, w_c, impChem, expChem, nymd, nhms, tdt, cdt, doChem, rc)
! !USES:

   IMPLICIT none

! !INPUT/OUTPUT PARAMETERS:

   TYPE(GMI_GridComp), INTENT(inOut) :: gcGMI   ! Grid Component
   TYPE(Chem_Bundle),  INTENT(inOut) :: w_c     ! Chemical tracer fields   
   TYPE(ESMF_State),   INTENT(inOut) :: impChem ! Import State
   TYPE(ESMF_State),   INTENT(inOut) :: expChem ! Export State

! !INPUT PARAMETERS:

   INTEGER, INTENT(in) :: nymd, nhms          ! time
   REAL,    INTENT(in) :: tdt                 ! heartbeat (secs)
   REAL,    INTENT(in) :: cdt                 ! chemical timestep (secs)
   LOGICAL, INTENT(in) :: doChem              ! whether the alarm is ringing

! !OUTPUT PARAMETERS:
   INTEGER, INTENT(out) ::  rc                ! Error return code:
                                              !  0 - all is well
                                              !  1 -

! !DESCRIPTION: This routine implements the GMI Strat/Trop Driver. That 
!               is, adds chemical tendencies to each of the constituents
!               Phase 2:  deposition         (done at the heartbeat)
!                         rest of chemistry  (done at GMI timestep)
!
! !IMPLEMENTATION NOTES:
!
!  No pointer is reservered in the export state for deposition of water.
!
! !DEFINED PARAMETERS:
      CHARACTER(LEN=*), PARAMETER :: IAm    = 'GMI_GridCompRun2'
!
! !LOCAL VARIABLES:
      INTEGER :: STATUS
!
! !REVISION HISTORY:
!
!  18Sep2003 da Silva  First crack.
!  30Jun2007 Nielsen   GMI Combo set up
!   7Apr2016 Nielsen   Initialize rc
!
!EOP
!-------------------------------------------------------------------------
!BOC
      rc = 0

      CALL GmiDepos_GridCompRun       (gcGMI%gcDepos,     w_c, impChem, expChem, nymd, nhms, tdt, __RC__)  ! tdt

      IF ( doChem ) THEN

        CALL GmiSAD_GridCompRun       (gcGMI%gcSAD,       w_c, impChem, expChem, nymd, nhms, cdt, __RC__)

        CALL GmiPhotolysis_GridCompRun(gcGMI%gcPhot,      w_c, impChem, expChem, nymd, nhms, cdt, __RC__)

        CALL GmiThermalRC_GridCompRun (gcGMI%gcThermalRC, w_c, impChem, expChem, nymd, nhms, cdt, __RC__)

        CALL GmiChemistry_GridCompRun (gcGMI%gcChem,      w_c, impChem, expChem, nymd, nhms, cdt, __RC__)

      ENDIF

      RETURN

      END SUBROUTINE GMI_GridCompRun2
!EOC
!-------------------------------------------------------------------------
!NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1, GEOS/DAS!
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  GMI_GridCompRunOrig --- The ORIGINAL GMI Run Driver
!
! !INTERFACE:
!
   SUBROUTINE GMI_GridCompRunOrig(gcGMI, w_c, impChem, expChem, nymd, nhms, cdt, clock, rc)

! !USES:

   IMPLICIT none

! !INPUT/OUTPUT PARAMETERS:

   TYPE(GMI_GridComp), INTENT(inOut) :: gcGMI   ! Grid Component
   TYPE(Chem_Bundle),  INTENT(inOut) :: w_c     ! Chemical tracer fields   
   TYPE(ESMF_State),   INTENT(inOut) :: impChem ! Import State
   TYPE(ESMF_State),   INTENT(inOut) :: expChem ! Export State
   TYPE(ESMF_Clock),   INTENT(inOut) :: clock   ! The clock

! !INPUT PARAMETERS:

   INTEGER, INTENT(in) :: nymd, nhms	      ! time
   REAL,    INTENT(in) :: cdt		      ! chemical timestep (secs)

! !OUTPUT PARAMETERS:
   INTEGER, INTENT(out) ::  rc                ! Error return code:
                                              !  0 - all is well
                                              !  1 -

! !DESCRIPTION: This routine implements the GMI Strat/Trop Driver. That 
!               is, adds chemical tendencies to each of the constituents
!               This routine includes all GMI processes, and reflects the
!               order in which they were called *before* the Run1/Run2
!               approach was implemented.
!
! !IMPLEMENTATION NOTES:
!
!  No pointer is reservered in the export state for deposition of water.
!
! !DEFINED PARAMETERS:
      CHARACTER(LEN=*), PARAMETER :: IAm    = 'GMI_GridCompRunOrig'
!
! !LOCAL VARIABLES:
      INTEGER :: STATUS
      LOGICAL :: mixPBL  ! whether to explicitly distribute
                         ! aerosol emissions within the PBL
!
! !REVISION HISTORY:
!
!  18Sep2003 da Silva  First crack.
!  30Jun2007 Nielsen   GMI Combo set up
!   7Apr2016 Nielsen   Initialize rc
!
!EOP
!-------------------------------------------------------------------------
!BOC
      rc = 0

      mixPBL = .TRUE.

      CALL GmiDepos_GridCompRun     (gcGMI%gcDepos,     w_c, impChem, expChem, nymd, nhms, cdt,               __RC__)

      CALL GmiEmiss_GridCompRun     (gcGMI%gcEmiss,     w_c, impChem, expChem, nymd, nhms, cdt, clock, mixPBL, __RC__)

      CALL GmiSAD_GridCompRun       (gcGMI%gcSAD,       w_c, impChem, expChem, nymd, nhms, cdt,               __RC__)

      CALL GmiPhotolysis_GridCompRun(gcGMI%gcPhot,      w_c, impChem, expChem, nymd, nhms, cdt,               __RC__)

      CALL GmiForcingBC_GridCompRun (gcGMI%gcFBC,       w_c, impChem, expChem, nymd, nhms, cdt,               __RC__)

      CALL GmiThermalRC_GridCompRun (gcGMI%gcThermalRC, w_c, impChem, expChem, nymd, nhms, cdt,               __RC__)

      CALL GmiChemistry_GridCompRun (gcGMI%gcChem,      w_c, impChem, expChem, nymd, nhms, cdt,               __RC__)

      RETURN

      END SUBROUTINE GMI_GridCompRunOrig
!EOC
!-------------------------------------------------------------------------
!NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1, GEOS/DAS!
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  GMI_GridCompFinalize
!
! !INTERFACE:
!

   SUBROUTINE GMI_GridCompFinalize(gcGMI, w_c, impChem, expChem, nymd, nhms, cdt, rc)

   IMPLICIT none

! !INPUT/OUTPUT PARAMETERS:

   TYPE(GMI_GridComp), INTENT(inOut) :: gcGMI	! Grid Component
   TYPE(ESMF_State),   INTENT(inOut) :: impChem ! Import State
   TYPE(ESMF_State),   INTENT(inOut) :: expChem ! Import State

! !INPUT PARAMETERS:

   TYPE(Chem_Bundle), INTENT(in)  :: w_c      ! Chemical tracer fields   
   INTEGER, INTENT(in) :: nymd, nhms	      ! time
   REAL,    INTENT(in) :: cdt  	              ! chemical timestep (secs)


! !OUTPUT PARAMETERS:

   INTEGER, INTENT(out) ::  rc                  ! Error return code:
                                                !  0 - all is well
                                                !  1 -
! 
! !DESCRIPTION: 
!  This routine finalizes this Grid Component.
!
! !DEFINED PARAMETERS:
      CHARACTER(LEN=*), PARAMETER :: IAm    = 'GMI_GridCompFinalize'
!
! !LOCAL VARIABLES:
      INTEGER :: STATUS
!
! !REVISION HISTORY:
!
!  18Sep2003 da Silva  First crack.
!  30Jun2007 Nielsen   GMI Combo set up
!   7Apr2016 Nielsen   Initialize rc
!
!EOP
!-------------------------------------------------------------------------
!BOC
      rc = 0

      CALL GmiDepos_GridCompFinalize     (gcGMI%gcDepos,     w_c, impChem, expChem, nymd, nhms, cdt, __RC__)

      CALL GmiEmiss_GridCompFinalize     (gcGMI%gcEmiss,     w_c, impChem, expChem, nymd, nhms, cdt, __RC__)

      CALL GmiSAD_GridCompFinalize       (gcGMI%gcSAD,       w_c, impChem, expChem, nymd, nhms, cdt, __RC__)

      CALL GmiPhotolysis_GridCompFinalize(gcGMI%gcPhot,      w_c, impChem, expChem, nymd, nhms, cdt, __RC__)

      CALL GmiForcingBC_GridCompFinalize (gcGMI%gcFBC,       w_c, impChem, expChem, nymd, nhms, cdt, __RC__)

      CALL GmiThermalRC_GridCompFinalize (gcGMI%gcThermalRC, w_c, impChem, expChem, nymd, nhms, cdt, __RC__)

      CALL GmiChemistry_GridCompFinalize (gcGMI%gcChem,      w_c, impChem, expChem, nymd, nhms, cdt, __RC__)

      RETURN

      END SUBROUTINE GMI_GridCompFinalize
!EOC
!------------------------------------------------------------------------------
  
      END MODULE GMI_GridCompMod
