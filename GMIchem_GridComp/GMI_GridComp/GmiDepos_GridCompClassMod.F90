#include "MAPL_Generic.h"
!-------------------------------------------------------------------------
!NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1, GEOS/DAS!
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  GmiDepos_GCCMod --- GMI Deposition Grid Component Class
!
! Grid Component class for the GMI Deposition
!
! !INTERFACE:
!
   MODULE  GmiDepos_GCCMod
!
! !USES:
!
   USE ESMF
   USE MAPL
   USE Chem_Mod 	     ! Chemistry Base Class
   USE Chem_UtilMod

   USE Species_BundleMod

   USE GmiEmissionMethod_mod,         ONLY : t_Emission
   USE GmiDepositionMethod_mod,       ONLY : t_Deposition
   USE GmiSpcConcentrationMethod_mod, ONLY : t_SpeciesConcentration
   USE GmiGrid_mod,                   ONLY : t_gmiGrid
   USE GmiTimeControl_mod,            ONLY : t_GmiClock
   USE GmiESMFrcFileReading_mod,      ONLY : rcEsmfReadTable
   USE GmiSwapSpeciesBundlesMod,      ONLY : SwapSpeciesBundles, speciesReg_for_CCM
   USE VegLaiMod,                     ONLY : Decode_Land_Types, Decode_XLAI

   IMPLICIT NONE
   INTEGER, PARAMETER :: DBL = KIND(0.00D+00)

#include "setkin_par.h"
#include "GmiParameters.h"
#include "gmi_phys_constants.h"
#include "gmi_emiss_constants.h"
#include "setkin_mw.h"
#include "setkin_lchem.h"

! !TYPES:

   PRIVATE
   PUBLIC  GmiDepos_GridComp       ! The GMI object 

! !PUBLIC MEMBER FUNCTIONS:

   PUBLIC  GmiDepos_GridCompInitialize
   PUBLIC  GmiDepos_GridCompRun
   PUBLIC  GmiDepos_GridCompFinalize

!
! !DESCRIPTION:
!
!  This module implements the GMI Deposition.
!
! !REVISION HISTORY:
!
!  16Sep2003 da Silva  First crack.
!  24Jan2O05 Nielsen   Implementation of Code 916 chemistry
!  19Dec2005 d Silva   Minor portability mods.
!  30Oct2007 Nielsen   GMI Combo set up
!  01Jul2010 Kouatchou Created this Deposition class
!  01Jun2015 Nielsen   ExtData replaces Chem_UtilMPread
!
!EOP
!-------------------------------------------------------------------------

  TYPE GmiDepos_GridComp
   CHARACTER(LEN=255) :: name = "GMI Deposition"

! Does the gmichem_import_rst file exist?  If not, some import states will be
! unfilled, and the species will need to "freewheel" through the first time step.
! -------------------------------------------------------------------------------
   LOGICAL :: gotImportRst

! Set BCRealTime = .TRUE. when boundary conditions 
! must be for exact year of current calendar date.
! -------------------------------------------------
   LOGICAL :: BCRealTime

! Dimensions
! ----------
   INTEGER :: i1, i2, im, j1, j2, jm, km

! Various switches
! ----------------
   LOGICAL :: pr_dry_depos
   LOGICAL :: pr_wet_depos   
   LOGICAL :: pr_scav
   LOGICAL :: pr_diag
   LOGICAL :: do_drydep
   LOGICAL :: do_wetdep
   LOGICAL :: do_synoz
   LOGICAL :: do_grav_set

   INTEGER :: chem_opt
   character (len=max_LENGTH_MET_NAME) :: metdata_name_org
   character (len=max_LENGTH_MET_NAME) :: metdata_name_model

! Surface area of grid cells
! --------------------------
   REAL(KIND=DBL), POINTER :: cellArea(:,:)

! Longitudes and latitudes (radians)
! ----------------------------------
   REAL, POINTER :: lonRad(:,:)
   REAL, POINTER :: latRad(:,:)

! Extra diagnostics
! -----------------
   LOGICAL :: verbose

! Map GMI species indices to CCM indices
! --------------------------------------
   INTEGER, POINTER :: mapSpecies(:)

! Component derived type declarations
! -----------------------------------
   TYPE(t_Emission  )		:: Emission
   TYPE(t_Deposition)		:: Deposition
   TYPE(t_gmiGrid   )		:: gmiGrid
   TYPE(t_GmiClock  )           :: gmiClock
   TYPE(t_SpeciesConcentration) :: SpeciesConcentration

! Veg Fraction does not change, only read once
! --------------------------------------------
   LOGICAL :: veg_fraction_done
 
  END TYPE GmiDepos_GridComp
!EOP
!-------------------------------------------------------------------------
CONTAINS
!-------------------------------------------------------------------------
!NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1, GEOS/DAS!
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  GmiDepos_GridCompInitialize --- Initialize GmiDepos_GridComp
!
! !INTERFACE:
!

   SUBROUTINE GmiDepos_GridCompInitialize( self, bgg, bxx, impChem, expChem, nymd, nhms, &
                                      tdt, rc )

   USE GmiSpcConcentrationMethod_mod, ONLY : InitializeSpcConcentration
   USE GmiDepositionMethod_mod,       ONLY : InitializeDeposition
   USE GmiGrid_mod,		      ONLY : InitializeGmiGrid
   USE GmiTimeControl_mod,            ONLY : Set_curGmiDate, Set_curGmiTime
   USE GmiTimeControl_mod,            ONLY : Set_begGmiDate, Set_begGmiTime
   USE GmiTimeControl_mod,            ONLY : Set_numTimeSteps

   IMPLICIT none
   INTEGER, PARAMETER :: DBL = KIND(0.00D+00)


! !INPUT PARAMETERS:

   TYPE(Species_Bundle), INTENT(in) :: bgg             ! Transported GMI species
   TYPE(Species_Bundle), INTENT(in) :: bxx             ! Non-transported GMI species
   INTEGER, INTENT(IN) :: nymd, nhms		       ! Time from AGCM
   REAL,    INTENT(IN) :: tdt			       ! Chemistry time step (secs)

! !OUTPUT PARAMETERS:

   TYPE(GmiDepos_GridComp), INTENT(INOUT)  :: self      ! Grid Component
   TYPE(ESMF_State),   INTENT(INOUT)  :: impChem    ! Import State
   TYPE(ESMF_State),   INTENT(INOUT)  :: expChem    ! Export State

   INTEGER, INTENT(out) ::  rc                  ! Error return code:
                                                !  0 - all is well
                                                !  1 - 

! !DESCRIPTION: Initializes the GMI Grid Component. It primarily sets
!               the import state.
!
! !REVISION HISTORY:
!
!  18Sep2003 da Silva  First crack.
!  30Jun2007 Nielsen   GMI Combo set up
!
!EOP
!-------------------------------------------------------------------------

   CHARACTER(LEN=*), PARAMETER :: IAm = 'GmiDepos_GridCompInitialize'
   CHARACTER(LEN=255) :: rcfilen = 'GMI_GridComp.rc'
   CHARACTER(LEN=255) :: importRestartFile
   CHARACTER(LEN=255) :: string
   
   type (ESMF_Config) :: gmiConfigFile

   INTEGER :: ios, m, n, STATUS
   INTEGER :: i, i1, i2, ic, im, j, j1, j2, jm, k, km, kReverse

   INTEGER :: i1_gl, i2_gl, ju1_gl, j2_gl 
   INTEGER :: ju1, jv1, jv1_gl, j1p, j2p
   INTEGER :: k1, k2, k1_gl, k2_gl
   INTEGER :: ilong, ilat, ivert, itloop
   INTEGER :: NPIJ, NPI, NPJ
   INTEGER :: ilo, ihi, julo, jvlo, jhi
   INTEGER :: ilo_gl, ihi_gl, julo_gl, jvlo_gl, jhi_gl
   INTEGER :: gmi_nborder
   INTEGER :: NMR      ! number of species from the GMI_Mech_Registry.rc

   INTEGER :: loc_proc
   LOGICAL :: one_proc, rootProc
   LOGICAL :: exists,open,found
   
   REAL :: qmin, qmax, tokgCPerBox

   REAL, ALLOCATABLE :: veg2D(:,:)
   character (len=3) :: vegID
   character (len=15) :: vegName

! Grid cell area can be set by initialize
! ---------------------------------------
   REAL, POINTER, DIMENSION(:,:) :: cellArea

   self%name = 'GMI Deposition'

!  Initialize local variables
!  --------------------------
   rc = 0
   i1 = self%i1
   i2 = self%i2
   im = self%im
   
   j1 = self%j1
   j2 = self%j2
   jm = self%jm
   
   km = self%km

     !-------------------
     ! Load resource file
     !-------------------

      IF ( MAPL_AM_I_ROOT() ) THEN
         PRINT *," "
         PRINT *,TRIM(IAm)//":"
         PRINT *,"Starting Reading the GMI Deposition Resource File"
      ENDIF

      gmiConfigFile = ESMF_ConfigCreate(rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigLoadFile(gmiConfigFile, TRIM(rcfilen), rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(gmiConfigFile, importRestartFile, &
     &                label   = "importRestartFile:", &
     &                default = ' ', rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(gmiConfigFile, value=self%do_synoz, &
     &           label="do_synoz:", default=.false., rc=STATUS)
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(gmiConfigFile, self%chem_opt, &
     &                label   = "chem_opt:", &
     &                default = 2, rc=STATUS )
      VERIFY_(STATUS)

      !------------------------------
      ! Deposition related variables
      !------------------------------
      
      call ESMF_ConfigGetAttribute(gmiConfigFile, value=self%do_drydep, &
     &           label="do_drydep:", default=.false., rc=STATUS)
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(gmiConfigFile, value=self%do_wetdep, &
     &           label="do_wetdep:", default=.false., rc=STATUS)
      VERIFY_(STATUS)

      !------------------------------
      ! Diagnostics related variables
      !------------------------------

      call ESMF_ConfigGetAttribute(gmiConfigFile, value=self%do_grav_set, &
     &           label="do_grav_set:", default=.false., rc=STATUS)
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(gmiConfigFile, value=self%pr_diag, &
     &           label="pr_diag:", default=.false., rc=STATUS)
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(gmiConfigFile, value=self%verbose, &
     &           label="verbose:", default=.false., rc=STATUS)
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(gmiConfigFile, value=self%pr_scav, &
     &           label="pr_scav:", default=.false., rc=STATUS)
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(gmiConfigFile, value=self%pr_dry_depos, &
     &           label="pr_dry_depos:", default=.false., rc=STATUS)

      call ESMF_ConfigGetAttribute(gmiConfigFile, value=self%pr_wet_depos, &
     &           label="pr_wet_depos:", default=.false., rc=STATUS)
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(gmiConfigFile, self%metdata_name_org, &
     &                label   = "metdata_name_org:", &
     &                default = 'GMAO', rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(gmiConfigFile, self%metdata_name_model, &
     &                label   = "metdata_name_model:", &
     &                default = 'GEOS-5', rc=STATUS )
      VERIFY_(STATUS)


! Does the GMICHEM import restart file exist?  If not,
! the species must "freewheel" through the first time step.
! ---------------------------------------------------------
   INQUIRE(FILE=TRIM(importRestartFile),EXIST=self%gotImportRst)
   IF( MAPL_AM_I_ROOT() ) THEN
    PRINT *," ",TRIM(importRestartFile)," exists: ",self%gotImportRst
    PRINT *," "
   END IF

!  GMI grid specification
!  ----------------------
   gmi_nborder = 0
   i1_gl = 1
   i2_gl = i2
   ju1_gl = 1
   jv1_gl = 1
   j2_gl = j2
   ju1 = j1
   jv1 =j1
   k1 = 1
   k2 = km
   k1_gl = 1
   k2_gl = km
   NPIJ = 16     ! These three integers are irrelevant
   NPI = 4
   NPJ = 4
   ilo = i1 - gmi_nborder
   ihi = i2 + gmi_nborder 
   julo = ju1 - gmi_nborder
   jvlo = jv1 - gmi_nborder
   jhi = j2 + gmi_nborder
   ilo_gl = i1_gl  - gmi_nborder
   ihi_gl = i2_gl  + gmi_nborder
   julo_gl = ju1_gl - gmi_nborder
   jvlo_gl = jv1_gl - gmi_nborder
   jhi_gl = j2_gl  + gmi_nborder
   j1p = 0
   j2p = j2_gl - j1p + 1
   ilong = i2 - i1 + 1
   ilat = j2 - ju1 + 1
   ivert = k2 - k1 + 1
   itloop = ilat * ilong * ivert

   one_proc = .FALSE.
   loc_proc = -99
   
   rootProc=.FALSE.
   IF( MAPL_AM_I_ROOT() ) THEN
    rootProc=.TRUE.
   END IF

! Set GMI's clock
! ---------------
   CALL Set_begGmiDate(self%gmiClock, nymd)
   CALL Set_begGmiTime(self%gmiClock, nhms)
   CALL Set_curGmiDate(self%gmiClock, nymd)
   CALL Set_curGmiTime(self%gmiClock, nhms)
   CALL Set_numTimeSteps(self%gmiClock, 0)

! Grid box surface area, m^{2}
! ----------------------------
   CALL MAPL_GetPointer(impChem, cellArea, 'AREA', rc=STATUS)
   VERIFY_(STATUS)
   ALLOCATE(self%cellArea(i1:i2,j1:j2), STAT=STATUS)
   VERIFY_(STATUS)
   self%cellArea(i1:i2,j1:j2)=cellArea(i1:i2,j1:j2)

! Discretization
! --------------
   CALL InitializeGmiGrid(self%gmiGrid, NPIJ, NPI, NPJ, &
                          gmi_nborder, i1, i2, ju1, jv1, j2, k1, k2, &
                          i1_gl, i2_gl, ju1_gl, jv1_gl, j2_gl, k1_gl, k2_gl, &
                          ilo, ihi, julo, jvlo, jhi, &
                          ilo_gl, ihi_gl, julo_gl, jvlo_gl, jhi_gl, &
                          ilong, ilat, ivert, itloop, j1p, j2p)  

! Perform a consistency check with setkin_par.h.
!   NSP is the number of species in setkins
!   NMR is the number of species in the Mech Registry
!
!   H2O   is in setkins,           but not in GMI_Mech_Registry
!   AOA   is in GMI_Mech_Registry, but not in setkins
!   T2M15 is in GMI_Mech_Registry, but not in setkins
!
!   The number of species common to both is therefore
!     NMR - 2     ! number in GMI_Mech_Registry - 2
!     NSP - 1     ! number in setkins - 1
! --------------------------------------------------------------------------------
   NMR = bgg%nq + bxx%nq
   IF( NMR-2 /= NSP-1 ) THEN
    PRINT *,TRIM(IAm),': Number of species from GMI_Mech_Registry.rc does not match number in setkin_par.h'
    STATUS = 1
    VERIFY_(STATUS)
   END IF

! Allocate space, etc., but the initialization of the
! species from the internal state is left to the run method.
! ----------------------------------------------------------

      CALL InitializeSpcConcentration(self%SpeciesConcentration,              &
                     self%gmiGrid, gmiConfigFile, NSP, NMF, NCHEM,            &
                     loc_proc)

      CALL InitializeDeposition(self%Deposition, self%gmiGrid,               &
                     gmiConfigFile, NSP, loc_proc, self%pr_diag,      &
                     self%pr_dry_depos, self%pr_wet_depos, self%pr_scav)

      !#################################################################
      ! This section of the code was included to be able to read some 
      ! Emission related variables.
      ! Later, we plan to use the export state of the Emission component
      !#################################################################

      allocate(self%Emission%ireg (i1:i2, ju1:j2))
      allocate(self%Emission%iland(i1:i2, ju1:j2, NTYPE))
      allocate(self%Emission%iuse (i1:i2, ju1:j2, NTYPE))
      allocate(self%Emission%xlai (i1:i2, ju1:j2, NTYPE))
      allocate(self%Emission%xlai2(i1:i2, ju1:j2, NTYPE))

      self%Emission%ireg  = 0
      self%Emission%iland = 0
      self%Emission%iuse  = 0
      self%Emission%xlai  = 0.0d0
      self%Emission%xlai2 = 0.0d0

      call ESMF_ConfigGetAttribute(gmiConfigFile, self%Emission%emiss_opt, &
     &                label   = "emiss_opt:", &
     &                default = 0, rc=STATUS )
      VERIFY_(STATUS)

    !---------------------------------------------------------------
    ! Create and populate the array that maps GMI species indices to
    ! GEOS-5 species indices
    !---------------------------------------------------------------

    allocate(self%mapSpecies(NSP))
    self%mapSpecies(:) = speciesReg_for_CCM(lchemvar, NSP, bgg%reg%vname, bxx%reg%vname )


    self%veg_fraction_done = .FALSE.

  RETURN

  END SUBROUTINE GmiDepos_GridCompInitialize

!-------------------------------------------------------------------------
!NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1, GEOS/DAS!
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  GmiDepos_GridCompRun --- The GMI Driver
!
! !INTERFACE:
!

   SUBROUTINE GmiDepos_GridCompRun ( self, bgg, bxx, impChem, expChem, nymd, nhms, &
                                tdt, rc )

! !USES:

   USE GmiTimeControl_mod,            ONLY : Set_curGmiDate, Set_curGmiTime
   USE GmiTimeControl_mod,            ONLY : Set_numTimeSteps, Get_numTimeSteps
   USE GmiTimeControl_mod,            ONLY : Set_gmiSeconds, GetSecondsFromJanuary1
   USE GmiSpcConcentrationMethod_mod, ONLY : resetFixedConcentration
   USE GmiSolar_mod,                  ONLY : CalcCosSolarZenithAngle
   USE GmiDepositionMethod_mod,       ONLY : RunDryDeposition, RunWetDeposition
   USE GmiDepositionMethod_mod,       ONLY : Set_dry_depos, Set_wet_depos, Set_scav3d
   USE GmiDepositionMethod_mod,       ONLY : Get_dry_depos, Get_wet_depos, Get_scav3d
   USE GmiGravitationalSettling_mod,  ONLY : updateGravitationalSettling

   IMPLICIT none

! !INPUT/OUTPUT PARAMETERS:

   TYPE(GmiDepos_GridComp), INTENT(INOUT) :: self ! Grid Component
   TYPE(Species_Bundle), INTENT(INOUT) :: bgg             ! Transported GMI species
   TYPE(Species_Bundle), INTENT(INOUT) :: bxx             ! Non-transported GMI species

! !INPUT PARAMETERS:

   TYPE(ESMF_State), INTENT(INOUT) :: impChem ! Import State
   INTEGER, INTENT(IN) :: nymd, nhms	      ! time
   REAL,    INTENT(IN) :: tdt		      ! chemical timestep (secs)

! !OUTPUT PARAMETERS:

   TYPE(ESMF_State), INTENT(INOUT) :: expChem   ! Export State
   INTEGER, INTENT(OUT) ::  rc                  ! Error return code:
                                                !  0 - all is well
                                                !  1 -

! !DESCRIPTION: This routine implements the GMI Strat/Trop Driver. That 
!               is, adds chemical tendencies to each of the constituents
!
! !IMPLEMENTATION NOTES:
!
!  No pointer is reservered in the export state for deposition of water.
!
! !REVISION HISTORY:
!
!  18Sep2003 da Silva  First crack.
!  24Jan2005 Nielsen   Implementation of Code 916 chemistry
!  30Oct2007 Nielsen   Implementation of GMI cmbined 
!                       stratosphere/troposphere chemistry
!
!EOP
!-------------------------------------------------------------------------

   CHARACTER(LEN=*), PARAMETER :: IAm = 'GmiDepos_GridCompRun'

!  Imports from GEOS-5
!  -------------------
   REAL, POINTER, DIMENSION(:,:) :: cn_prcp, tprec, lwi, TsAir
   REAL, POINTER, DIMENSION(:,:) :: ustar, z0h, swndsrf, T2m
   REAL, POINTER, DIMENSION(:,:) :: cldtt, area, frlandice, snowdp

   REAL, POINTER, DIMENSION(:,:,:) :: airdens, ple, Q, T, zle
   REAL, POINTER, DIMENSION(:,:,:) :: dqdt,pfl_lsan,pfl_cn

!  Exports
!  -------
#include "Deposition_DeclarePointer___.h"

!  Local
!  -----
   INTEGER :: cymd, dymd, emiss_opt, hms
   INTEGER :: i, i1, i2, ic, im, iXj
   INTEGER :: j, j1, j2, jm
   INTEGER :: k, km, kReverse
   INTEGER :: loc_proc
   INTEGER :: n, STATUS

   INTEGER, PARAMETER :: ToGMI = 1
   INTEGER, PARAMETER :: FromGMI = -1

   REAL :: pi,degToRad,radToDeg,OneOverDt

   REAL, PARAMETER :: mwtAir = 28.9
   REAL, PARAMETER :: rStar = 8.314E+03
   REAL, PARAMETER :: Pa2hPa = 0.01
   REAL, PARAMETER :: ToGrPerKg = 1000.00
   REAL, PARAMETER :: secPerDay = 86400.00
   REAL, PARAMETER :: err = 1.00E-04

   REAL(KIND=DBL) :: chemDt, dayOfYear

   CHARACTER(LEN=255) :: speciesName
   CHARACTER(LEN=255) :: importName

   LOGICAL :: found, rootProc
   LOGICAL, PARAMETER :: doThis = .FALSE.
   
! Allocatables.  Use KIND=DBL where GMI expects REAL*8.
! -----------------------------------------------------
   INTEGER, ALLOCATABLE :: lwis_flags(:,:)

   REAL, ALLOCATABLE :: var3d(:,:,:)
   REAL, ALLOCATABLE :: pl(:,:,:)

   REAL(KIND=DBL), ALLOCATABLE :: lonDeg(:,:)
   REAL(KIND=DBL), ALLOCATABLE :: latDeg(:,:)

   REAL(KIND=DBL), ALLOCATABLE :: TwoMeter_air_temp(:,:)
   REAL(KIND=DBL), ALLOCATABLE :: fracCloudCover(:,:)
   REAL(KIND=DBL), ALLOCATABLE :: surf_rough(:,:)
   REAL(KIND=DBL), ALLOCATABLE :: cosSolarZenithAngle(:,:)
   REAL(KIND=DBL), ALLOCATABLE :: radswg(:,:)
   REAL(KIND=DBL), ALLOCATABLE :: frictionVelocity(:,:)
   REAL(KIND=DBL), ALLOCATABLE :: con_precip(:,:)
   REAL(KIND=DBL), ALLOCATABLE :: tot_precip(:,:)

   REAL(KIND=DBL), ALLOCATABLE :: var3dDBL(:,:,:)
   REAL(KIND=DBL), ALLOCATABLE :: var4dDBL(:,:,:,:)

   REAL(KIND=DBL), ALLOCATABLE :: mass(:,:,:)
   REAL(KIND=DBL), ALLOCATABLE :: press3c(:,:,:)
   REAL(KIND=DBL), ALLOCATABLE :: press3e(:,:,:)
   REAL(KIND=DBL), ALLOCATABLE :: gridBoxThickness(:,:,:)
   REAL(KIND=DBL), ALLOCATABLE :: kel(:,:,:)
   REAL(KIND=DBL), ALLOCATABLE :: humidity(:,:,:)
   REAL(KIND=DBL), ALLOCATABLE :: moistq(:,:,:)
   REAL(KIND=DBL), ALLOCATABLE :: rain3Dcn(:,:,:)
   REAL(KIND=DBL), ALLOCATABLE :: rain3Dls(:,:,:)

   REAL(KIND=DBL), ALLOCATABLE :: diffaer(:,:,:)
   REAL(KIND=DBL), ALLOCATABLE :: s_radius(:,:,:)
   REAL(KIND=DBL), ALLOCATABLE :: s_velocity(:,:,:)

   loc_proc = -99

!  Grid specs from Chem_Bundle%grid
!  --------------------------------
   rc = 0
   i1 = self%i1
   i2 = self%i2
   im = self%im
   
   j1 = self%j1
   j2 = self%j2
   jm = self%jm
   
   km = self%km
   
   iXj = (i2-i1+1)*(j2-j1+1)

!  Some real constants
!  -------------------
   pi = 4.00*ATAN(1.00)
   degToRad = pi/180.00
   radToDeg = 180.00/pi
   chemDt = tdt
   OneOverDt = 1.00/tdt

   rootProc=.FALSE.
   IF( MAPL_AM_I_ROOT() ) THEN
    rootProc=.TRUE.
   END IF

!  We need lots of pointers!
!  -------------------------
   CALL FindPointers(rc)

!  Reserve some local work space
!  -----------------------------
   ALLOCATE(lonDeg(i1:i2,j1:j2),STAT=STATUS)
   VERIFY_(STATUS)
   ALLOCATE(latDeg(i1:i2,j1:j2),STAT=STATUS)
   VERIFY_(STATUS)

   ALLOCATE(         lwis_flags(i1:i2,j1:j2),STAT=STATUS)
   VERIFY_(STATUS)
   ALLOCATE(     fracCloudCover(i1:i2,j1:j2),STAT=STATUS)
   VERIFY_(STATUS)
   ALLOCATE(         surf_rough(i1:i2,j1:j2),STAT=STATUS)
   VERIFY_(STATUS)
   ALLOCATE(cosSolarZenithAngle(i1:i2,j1:j2),STAT=STATUS)
   VERIFY_(STATUS)
   ALLOCATE(             radswg(i1:i2,j1:j2),STAT=STATUS)
   VERIFY_(STATUS)
   ALLOCATE(   frictionVelocity(i1:i2,j1:j2),STAT=STATUS)
   VERIFY_(STATUS)
   ALLOCATE(         con_precip(i1:i2,j1:j2),STAT=STATUS)
   VERIFY_(STATUS)
   ALLOCATE(         tot_precip(i1:i2,j1:j2),STAT=STATUS)
   VERIFY_(STATUS)
   ALLOCATE(  TwoMeter_air_temp(i1:i2,j1:j2),STAT=STATUS)
   VERIFY_(STATUS)

   ALLOCATE(                pl(i1:i2,j1:j2,1:km),STAT=STATUS)
   VERIFY_(STATUS)
   ALLOCATE(             var3d(i1:i2,j1:j2,1:km),STAT=STATUS)
   VERIFY_(STATUS)
   ALLOCATE(              mass(i1:i2,j1:j2,1:km),STAT=STATUS)
   VERIFY_(STATUS)
   ALLOCATE(           press3c(i1:i2,j1:j2,1:km),STAT=STATUS)
   VERIFY_(STATUS)
   ALLOCATE(           press3e(i1:i2,j1:j2,0:km),STAT=STATUS)
   VERIFY_(STATUS)
   ALLOCATE(  gridBoxThickness(i1:i2,j1:j2,1:km),STAT=STATUS)
   VERIFY_(STATUS)
   ALLOCATE(               kel(i1:i2,j1:j2,1:km),STAT=STATUS)
   VERIFY_(STATUS)
   ALLOCATE(          humidity(i1:i2,j1:j2,1:km),STAT=STATUS)
   VERIFY_(STATUS)
   ALLOCATE(            moistq(i1:i2,j1:j2,1:km),STAT=STATUS)
   VERIFY_(STATUS)
   ALLOCATE(          rain3Dcn(i1:i2,j1:j2,1:km),STAT=STATUS)
   VERIFY_(STATUS)
   ALLOCATE(          rain3Dls(i1:i2,j1:j2,1:km),STAT=STATUS)
   VERIFY_(STATUS)

   ALLOCATE(   diffaer(i1:i2,j1:j2,NSP),STAT=STATUS)
   VERIFY_(STATUS)
   ALLOCATE(  s_radius(i1:i2,j1:j2,NSP),STAT=STATUS)
   VERIFY_(STATUS)
   ALLOCATE(s_velocity(i1:i2,j1:j2,NSP),STAT=STATUS)
   VERIFY_(STATUS)

! Geolocation
! -----------
   lonDeg(i1:i2,j1:j2)=self%lonRad(i1:i2,j1:j2)*radToDeg
   latDeg(i1:i2,j1:j2)=self%latRad(i1:i2,j1:j2)*radToDeg

!  Layer mean pressures. NOTE: ple(:,:,0:km)
!  -----------------------------------------
   DO k=1,km
    pl(i1:i2,j1:j2,k)=(ple(i1:i2,j1:j2,k-1)+ple(i1:i2,j1:j2,k))*0.50
   END DO
   
! Set GMI's clock
! ---------------
   CALL Set_curGmiDate(self%gmiClock, nymd)
   CALL Set_curGmiTime(self%gmiClock, nhms)
   CALL Get_numTimeSteps(self%gmiClock, ic)
   CALL Set_numTimeSteps(self%gmiClock, ic+1)
   CALL Set_gmiSeconds(self%gmiClock, (ic+1)*chemDt)

! Update the following time-dependent boundary conditions:
!  Fixed concentration species
!  Stratospheric sulfate surface area
!  Emissions
! --------------------------------------------------------
   CALL Acquire_Clims(STATUS)
   VERIFY_(STATUS)


! Grab imports and do units conversions
! -------------------------------------
   CALL SatisfyImports(STATUS)
   VERIFY_(STATUS)

! Hand the species concentrations to GMI's bundle
! -----------------------------------------------
   IF (self%gotImportRst) THEN
      CALL SwapSpeciesBundles(ToGMI, self%SpeciesConcentration%concentration, &
               bgg%qa, bxx%qa, Q, self%mapSpecies, lchemvar, self%do_synoz, NSP, &
               STATUS)
      VERIFY_(STATUS)
   END IF

   DEALLOCATE(var3d, STAT=STATUS)
   VERIFY_(STATUS)

! Impose fixed concentrations
! ---------------------------
   IF(self%SpeciesConcentration%num_fixed_const > 0) THEN
    CALL resetFixedConcentration(self%SpeciesConcentration, self%gmiClock, self%gmiGrid, NSP)
   END IF

! ------------------------------------------------------------------------
! Gravitational settling
! ------------------------------------------------------------------------
   IF(self%do_grav_set .AND. self%gotImportRst) THEN
    CALL updateGravitationalSettling(self%SpeciesConcentration,                &
                              gridBoxThickness, humidity, mass, press3e, kel,  &
                              diffaer, s_radius, s_velocity,                   &
                              self%pr_diag, loc_proc, self%chem_opt, i1, i2,   &
                              j1, j2, 1, km, i1, i2, j1, j2, NSP, tdt)
   END IF

! ------------------------------------------------------------------------
! Deposition
!
! NOTE: In GEOS-5, the "instantaneous" deposition is held in dry_depos and
!       wet_depos.  The accumulated deposition is obtained by instantiating
!       time-averaged exports via the HISTORY.rc in the run script.
! ------------------------------------------------------------------------

   IF(self%Deposition%do_drydep .AND. self%gotImportRst) THEN

    IF(self%pr_dry_depos) THEN
     ALLOCATE(var3dDBL(i1:i2,j1:j2,1:NSP),STAT=STATUS)
     VERIFY_(STATUS)
     var3dDBL(i1:i2,j1:j2,1:NSP) = 0.00D+00
     Call Set_dry_depos(self%Deposition, var3dDBL)
     DEALLOCATE(var3dDBL, STAT=STATUS)
     VERIFY_(STATUS)
    END IF
    CALL RunDryDeposition(self%Deposition, self%Emission, self%SpeciesConcentration, &
                          self%gmiGrid, lwis_flags, self%cellArea,                   &
                          cosSolarZenithAngle, fracCloudCover, radswg,               &
                          TwoMeter_air_temp, surf_rough, frictionVelocity, mass,     &
                          diffaer, s_radius, s_velocity, gridBoxThickness(:,:,1),    &
                          loc_proc, mw_data, NSP, self%chem_opt, self%pr_dry_depos,  &
                          self%pr_diag, tdt)
   END IF

   IF(self%Deposition%do_wetdep .AND. self%gotImportRst) THEN

    IF(self%pr_wet_depos) THEN
     ALLOCATE(var3dDBL(i1:i2,j1:j2,1:NSP),STAT=STATUS)
     VERIFY_(STATUS)
     var3dDBL(i1:i2,j1:j2,1:NSP) = 0.00D+00
     Call Set_wet_depos(self%Deposition, var3dDBL)
     DEALLOCATE(var3dDBL, STAT=STATUS)
     VERIFY_(STATUS)
    END IF
    IF(self%pr_scav) THEN
     ALLOCATE(var4dDBL(i1:i2,j1:j2,1:km,1:NSP),STAT=STATUS)
     VERIFY_(STATUS)
     var4dDBL(i1:i2,j1:j2,1:km,1:NSP) = 0.00D+00
     Call Set_scav3d(self%Deposition, var4dDBL)
     DEALLOCATE(var4dDBL, STAT=STATUS)
     VERIFY_(STATUS)
    END IF

    CALL RunWetDeposition(self%Deposition, self%SpeciesConcentration,                  &
                          self%gmiGrid, IH2O2, IHNO3, mw_data, con_precip,             &
                          tot_precip, self%cellArea, gridBoxThickness, mass, moistq,   &
                          rain3Dcn, rain3Dls, kel, press3c, press3e, loc_proc, NSP,    &
                          self%chem_opt, self%pr_wet_depos, self%pr_scav, self%pr_diag,&
                          tdt)
   END IF

! Return species concentrations to the chemistry bundle
! -----------------------------------------------------
   IF (self%gotImportRst) then
      CALL SwapSpeciesBundles(FromGMI, self%SpeciesConcentration%concentration, &
               bgg%qa, bxx%qa, Q, self%mapSpecies, lchemvar, self%do_synoz, NSP,  &
               STATUS)
      VERIFY_(STATUS)
   END IF

! Export states
! -------------
   CALL FillExports(STATUS)
   VERIFY_(STATUS)

! Scratch local work space
! ------------------------
   DEALLOCATE(lonDeg, latDeg, STAT=STATUS)
   VERIFY_(STATUS)

   DEALLOCATE(lwis_flags, TwoMeter_air_temp, &
              fracCloudCover, surf_rough, cosSolarZenithAngle, &
              radswg, frictionVelocity, con_precip, tot_precip, STAT=STATUS)
   VERIFY_(STATUS)

   DEALLOCATE(pl, mass, press3c, press3e, gridBoxThickness, kel, humidity, &
              moistq, rain3Dcn, rain3Dls, STAT=STATUS)
   VERIFY_(STATUS)

   DEALLOCATE(diffaer, s_radius, s_velocity, STAT=STATUS)
   VERIFY_(STATUS)

! IMPORTANT: Reset this switch to .TRUE. after first pass.
! --------------------------------------------------------
   self%gotImportRst = .TRUE.

   RETURN

CONTAINS

!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  Acquire_Clims
!
! !INTERFACE:

  SUBROUTINE Acquire_Clims(rc)
  
  IMPLICIT NONE

  INTEGER, INTENT(OUT) :: rc
!
! !DESCRIPTION:
!
!  Obtain climatologies from files specified in ExtData.rc.
!
!  For details on decoding land use and leaf area indices, see Chem_Shared/VegLaiMod.F90

!EOP
!---------------------------------------------------------------------------

  CHARACTER(LEN=255) :: IAm
  CHARACTER(LEN=  3) :: vegID
  INTEGER :: ic, m
  INTEGER, ALLOCATABLE :: landNum(:)
  INTEGER, ALLOCATABLE :: milFrac(:)
  
  REAL, POINTER, DIMENSION(:,:,:) :: PTR3D

  rc = 0
  IAm = "Acquire_Clims"

! Update fixed concentration species
! ----------------------------------
  IF(self%SpeciesConcentration%num_fixed_const > 0) THEN
   DO ic = 1,self%SpeciesConcentration%num_fixed_const
    i = self%SpeciesConcentration%fixed_const_map(ic)

    speciesName = TRIM(lchemvar(i))
    importName = TRIM(speciesName)//'_FIXED'
    CALL MAPL_GetPointer(impChem, PTR3D, TRIM(importName), RC=STATUS)
    VERIFY_(STATUS)
    self%SpeciesConcentration%fixed_const(i1:i2,j1:j2,1:km,ic) = PTR3D(i1:i2,j1:j2,km:1:-1)
    NULLIFY(PTR3D)

   END DO
  END IF

! Update the vegetation & LAI info
! --------------------------------
  IF(BTEST(self%Emission%emiss_opt,1) .OR. self%do_drydep) THEN

   IF ( .NOT. self%veg_fraction_done ) THEN

     CALL MAPL_GetPointer(impChem, PTR3D, 'VEG_FRAC', __RC__ )

     CALL Decode_Land_Types(PTR3D, NTYPE,              &
                            i1, i2, j1, j2, NVEGTYPE-2,&
                            self%Emission%ireg,        &
                            self%Emission%iuse,        &
                            self%Emission%iland,       &
                            RC=STATUS)
     VERIFY_(STATUS)

     NULLIFY(PTR3D)

     self%veg_fraction_done = .TRUE.

   ENDIF

   CALL MAPL_GetPointer(impChem, PTR3D, 'LAI_FRAC', __RC__ )

   CALL Decode_XLAI(PTR3D, NTYPE,              &
                    i1, i2, j1, j2, NVEGTYPE-2,&
                    self%Emission%ireg,        &
                    self%Emission%iuse,        &
                    self%Emission%iland,       &
                    self%Emission%xlai,        &
                    RC=STATUS)
   VERIFY_(STATUS)

   NULLIFY(PTR3D)

  END IF

  RETURN
 END SUBROUTINE Acquire_Clims

!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  FillExports
!
! !INTERFACE:

  SUBROUTINE FillExports(rc)
  
  IMPLICIT NONE

  INTEGER, INTENT(OUT) :: rc
!
! !DESCRIPTION:
!
!  Fill the export state
!
!EOP
!---------------------------------------------------------------------------

  CHARACTER(LEN=*), PARAMETER :: IAm= 'GmiDepos_FillExports'

  rc=0

! ---------------------------------------------------------------------
! To conform with GEOS-5's other chemistry components, GMICHEM's ozone
! states are as follows:
!
!    NAME   STATE    Units             Comments
!  ------ --------   ----------------- --------------------------------
!      OX  Internal  mol/mol           Reqired name for ANALYSIS bundle
!      O3  Export    kg/kg             OX(vmr)*48/28.97
!  03PPMV  Export    ppmv              OX(vmr)*1.00E+06
!
! ---------------------------------------------------------------------

! Convert to kg m^{-2} s^{-1} for export.
! WARNING: NOT able to export deposition of water.
! ------------------------------------------------
   IF(self%do_drydep .OR. self%do_wetdep) THEN
    ALLOCATE(var3dDBL(i1:i2,j1:j2,1:NSP),STAT=STATUS)
    VERIFY_(STATUS)
    ALLOCATE(var3d(i1:i2,j1:j2,1:NSP),STAT=STATUS)
    VERIFY_(STATUS)
   END IF

   IF(self%do_drydep .AND. self%pr_dry_depos) THEN
    CALL Get_dry_depos(self%Deposition, var3dDBL)
    var3d(i1:i2,j1:j2,1:NSP)=var3dDBL(i1:i2,j1:j2,1:NSP)/tdt
#include "DD_FillExports___.h"
   END IF

   IF(self%do_wetdep .AND. self%pr_wet_depos) THEN
    CALL Get_wet_depos(self%Deposition, var3dDBL)
    var3d(i1:i2,j1:j2,1:NSP)=var3dDBL(i1:i2,j1:j2,1:NSP)/tdt
#include "WD_FillExports___.h"
   END IF

   IF(self%do_drydep .OR. self%do_wetdep) THEN
    DEALLOCATE(var3d, var3dDBL, STAT=STATUS)
    VERIFY_(STATUS)
   ENDIF

   IF(self%do_wetdep .AND. self%pr_scav) THEN
    ALLOCATE(var4dDBL(i1:i2,j1:j2,1:km,1:NSP),STAT=STATUS)
    VERIFY_(STATUS)
    CALL Get_scav3d(self%Deposition, var4dDBL)
    DO k=1,km
     kReverse = km-k+1
#include "SCAV_FillExports___.h"
    END DO
    DEALLOCATE(var4dDBL,STAT=STATUS)
    VERIFY_(STATUS)
   END IF

  RETURN
 END SUBROUTINE FillExports

!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  FindPointers
!
! !INTERFACE:

  SUBROUTINE FindPointers(rc)
  
  IMPLICIT NONE

  INTEGER, INTENT(OUT) :: rc
!
! !DESCRIPTION:
!
!  Find pointers to import and export states
!
!EOP
!---------------------------------------------------------------------------

  CHARACTER(LEN=*), PARAMETER :: IAm = 'GmiDepos_FindPointers'

  REAL :: qmin,qmax
  INTEGER :: STATUS
  
  rc=0

!  Pointers to imports
!  -------------------
   CALL MAPL_GetPointer(impChem,   cn_prcp,   'CN_PRCP', rc=STATUS)
   VERIFY_(STATUS)
   CALL MAPL_GetPointer(impChem,     tprec,     'TPREC', rc=STATUS)
   VERIFY_(STATUS)
   CALL MAPL_GetPointer(impChem,       lwi,       'LWI', rc=STATUS)
   VERIFY_(STATUS)
   CALL MAPL_GetPointer(impChem, frlandice, 'FRLANDICE', rc=STATUS)
   VERIFY_(STATUS)
   CALL MAPL_GetPointer(impChem,    snowdp,    'SNOWDP', rc=STATUS)
   VERIFY_(STATUS)
   CALL MAPL_GetPointer(impChem,     TsAir,        'TA', rc=STATUS)
   VERIFY_(STATUS)
   CALL MAPL_GetPointer(impChem,       T2m,       'T2M', rc=STATUS)
   VERIFY_(STATUS)
   CALL MAPL_GetPointer(impChem,     ustar,     'USTAR', rc=STATUS)
   VERIFY_(STATUS)
   CALL MAPL_GetPointer(impChem,       z0h,       'Z0H', rc=STATUS)
   VERIFY_(STATUS)
   CALL MAPL_GetPointer(impChem,   swndsrf,   'SWNDSRF', rc=STATUS)
   VERIFY_(STATUS)
   CALL MAPL_GetPointer(impChem,     cldtt,     'CLDTT', rc=STATUS)
   VERIFY_(STATUS)
  
   CALL MAPL_GetPointer(impChem,   airdens, 'AIRDENS', rc=STATUS)
   VERIFY_(STATUS)
   CALL MAPL_GetPointer(impChem,       ple,	'PLE', rc=STATUS)
   VERIFY_(STATUS)
   CALL MAPL_GetPointer(impChem,         Q,       'Q', rc=STATUS)
   VERIFY_(STATUS)
   CALL MAPL_GetPointer(impChem,	 T,	  'T', rc=STATUS)
   VERIFY_(STATUS)
   CALL MAPL_GetPointer(impChem,       zle,	'ZLE', rc=STATUS)
   VERIFY_(STATUS)
   CALL MAPL_GetPointer(impChem,      dqdt,    'DQDT', rc=STATUS)
   VERIFY_(STATUS)
   CALL MAPL_GetPointer(impChem,  pfl_lsan,'PFL_LSAN', rc=STATUS)
   VERIFY_(STATUS)
   CALL MAPL_GetPointer(impChem,    pfl_cn,  'PFL_CN', rc=STATUS)
   VERIFY_(STATUS)

!  Export state pointers
!  ---------------------
#include "Deposition_GetPointer___.h"

!  Validation
!  ----------
   Validate: IF(self%verbose) THEN
    IF(MAPL_AM_I_ROOT( )) PRINT *,TRIM(IAm),": Input ..."
    CALL pmaxmin('CN_PRCP:', cn_prcp, qmin, qmax, iXj, 1, 1. )
    CALL pmaxmin('TPREC:', tprec, qmin, qmax, iXj, 1, 1. )
    CALL pmaxmin('LWI:', lwi, qmin, qmax, iXj, 1, 1. )
    CALL pmaxmin('FRLANDICE:', frlandice, qmin, qmax, iXj, 1, 1. )
    CALL pmaxmin('SNOWDP:', snowdp, qmin, qmax, iXj, 1, 1. )
    CALL pmaxmin('TSAIR:', TsAir, qmin, qmax, iXj, 1, 1. )
    CALL pmaxmin('T2M:', T2m, qmin, qmax, iXj, 1, 1. )
    CALL pmaxmin('Q:', Q, qmin, qmax, iXj, km, 1. )
    CALL pmaxmin('T:', T, qmin, qmax, iXj, km, 1. )
    CALL pmaxmin('ZLE:', zle, qmin, qmax, iXj, km+1, 1. )
    CALL pmaxmin('PLE (hPa):', ple, qmin, qmax, iXj, km+1, 0.01 )
    CALL pmaxmin('AIRDENS:', airdens, qmin, qmax, iXj, km, 1. )
    CALL pmaxmin('DQ/DT:', dqdt, qmin, qmax, iXj, km, 1. )
    CALL pmaxmin('PFL_LSAN:', pfl_lsan, qmin, qmax, iXj, km+1, 1. )
    CALL pmaxmin('PFL_CN:', pfl_cn, qmin, qmax, iXj, km+1, 1. )

    CALL pmaxmin('USTAR:', ustar, qmin, qmax, iXj, 1, 1. )
    CALL pmaxmin('SURF_ROUGH:', z0h, qmin, qmax, iXj, 1, 1. )
    CALL pmaxmin('SWNDSRF:', swndsrf, qmin, qmax, iXj, 1, 1. )
    CALL pmaxmin('CLDTT:', cldtt, qmin, qmax, iXj, 1, 1. )

   END IF Validate

  RETURN
 END SUBROUTINE FindPointers

!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  SatisfyImports
!
! !INTERFACE:

  SUBROUTINE SatisfyImports(rc)
  
  IMPLICIT NONE

  INTEGER, INTENT(OUT) :: rc
!
! !DESCRIPTION:
!
!  Grab imports, change units if necessary, and convert to
!  REAL(KIND=DBL) when required.
!
!EOP
!---------------------------------------------------------------------------

  CHARACTER(LEN=255) :: IAm
  
  INTEGER, ALLOCATABLE :: on(:,:)
  
  rc=0
  IAm="SatisfyImports"

! ------------------------------------------------------------------------ 
! Imports to REAL*8, apply units conversion, and reverse vertical stacking
! ------------------------------------------------------------------------

! Singly-layered                                                            GEOS-5 Units       GMI Units
! --------------                                                            ------------       -------------
  TwoMeter_air_temp(i1:i2,j1:j2) = T2m(i1:i2,j1:j2)                         ! K
  fracCloudCover(i1:i2,j1:j2) = cldtt(i1:i2,j1:j2)                          ! fraction
  surf_rough(i1:i2,j1:j2) = z0h(i1:i2,j1:j2)                                ! m
  radswg(i1:i2,j1:j2) = swndsrf(i1:i2,j1:j2)                                ! w m^{-2}
  frictionVelocity(i1:i2,j1:j2) = ustar(i1:i2,j1:j2)                        ! m s^{-1}
  con_precip(i1:i2,j1:j2) = cn_prcp(i1:i2,j1:j2)*secPerDay                  ! kg m^{-2}s^{-1}   mm d^{-1}
  tot_precip(i1:i2,j1:j2) = tprec(i1:i2,j1:j2)*secPerDay                    ! kg m^{-2}s^{-1}   mm d^{-1}

! Layer means                                                               GEOS-5 Units       GMI Units
! -----------                                                               ------------       -------------
  DO k=1,km
   kReverse = km-k+1                                                        ! Lid-to-surf      Surf-to-lid
   press3c(i1:i2,j1:j2,kReverse) = pl(i1:i2,j1:j2,k)*Pa2hPa                 ! Pa               hPa
   kel(i1:i2,j1:j2,kReverse) = T(i1:i2,j1:j2,k)                             ! K
   humidity(i1:i2,j1:j2,kReverse) = Q(i1:i2,j1:j2,k)*ToGrPerKg              ! kg kg^{-1}       g kg^{-1}
   moistq(i1:i2,j1:j2,kReverse) = dqdt(i1:i2,j1:j2,k)*ToGrPerKg*secPerDay   ! kg kg^{-1}s^{-1} g kg^{-1}d^{-1}
  END DO

! Layer edges                                                               GEOS-5 Units       GMI Units
! -----------                                                               ------------       -------------
  DO k=0,km
   kReverse = km-k
   press3e(i1:i2,j1:j2,kReverse) = ple(i1:i2,j1:j2,k)*Pa2hPa                ! Pa               hPa
  END DO

  DO k=0,km-1
   kReverse = km-k
   rain3Dls(i1:i2,j1:j2,kReverse) = pfl_lsan(i1:i2,j1:j2,k)*secPerDay       ! kg m^{-2}s^{-1}  mm d^{-1}
   rain3Dcn(i1:i2,j1:j2,kReverse) = pfl_cn(i1:i2,j1:j2,k)*secPerDay         ! kg m^{-2}s^{-1}  mm d^{-1}
  END DO

! Retain land-water-ice flag numbering from GEOS
! Make land-water-ice-snow flag for Dry deposition
! specification: 0=water 1=land 2=ice 3=snow/glaciated
! ----------------------------------------------------
   lwis_flags(i1:i2,j1:j2)=FLOOR(lwi(i1:i2,j1:j2)+0.1)
!  Note - SNOWDP can be undefined (BIG)
   WHERE( (frlandice(i1:i2,j1:j2) > 0.5 .OR. (snowdp(i1:i2,j1:j2) >= 0.35 .AND.     &
                                              snowdp(i1:i2,j1:j2) <  1000.0)    ) ) &
          lwis_flags(i1:i2,j1:j2) = 3

! Cell mass and thickness                                                   GEOS-5 Units       GMI Units
! -----------------------                                                   ------------       -------------
   DO k=1,km
    kReverse = km-k+1
    mass(:,:,kReverse)=airdens(:,:,k)*self%cellArea(:,:)* &                 ! kg
                       (zle(:,:,k-1)-zle(:,:,k))
    gridBoxThickness(:,:,kReverse) = zle(:,:,k-1)-zle(:,:,k)                ! m
   END DO

! Obtain instantaneous apparent sun
! ---------------------------------
   CALL GetSecondsFromJanuary1(ic, nymd, nhms)
   dayOfYear = (1.00*ic)/secPerDay
   CALL CalcCosSolarZenithAngle(dayOfYear, latDeg, lonDeg, cosSolarZenithAngle, &
                                i1, i2, j1, j2)

  RETURN
 END SUBROUTINE SatisfyImports

 END SUBROUTINE GmiDepos_GridCompRun

!-------------------------------------------------------------------------
!NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1, GEOS/DAS!
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  GmiDepos_GridCompFinalize
!
! !INTERFACE:
!

   SUBROUTINE GmiDepos_GridCompFinalize ( self, impChem, expChem, &
                                     nymd, nhms, cdt, rc )

  IMPLICIT none

! !INPUT/OUTPUT PARAMETERS:

   TYPE(GmiDepos_GridComp), INTENT(inout) :: self ! Grid Component

! !INPUT PARAMETERS:

   INTEGER, INTENT(in) :: nymd, nhms	      ! time
   REAL,    INTENT(in) :: cdt  	              ! chemical timestep (secs)


! !OUTPUT PARAMETERS:

   TYPE(ESMF_State), INTENT(inout) :: impChem	! Import State
   TYPE(ESMF_State), INTENT(inout) :: expChem	! Import State
   INTEGER, INTENT(out) ::  rc                  ! Error return code:
                                                !  0 - all is well
                                                !  1 -
 
! !DESCRIPTION: This routine finalizes this Grid Component.
!
! !REVISION HISTORY:
!
!  18Sep2003 da Silva  First crack.
!
!EOP
!-------------------------------------------------------------------------

   CHARACTER(LEN=*), PARAMETER :: IAm = 'GmiDepos_GridCompFinalize'
   INTEGER :: STATUS
   rc=0
   DEALLOCATE(self%cellArea, self%lonRad, self%latRad, STAT=STATUS)
   VERIFY_(STATUS)
   RETURN

 END SUBROUTINE GmiDepos_GridCompFinalize
  
 END MODULE GmiDepos_GCCMod

