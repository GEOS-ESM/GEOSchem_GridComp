#include "MAPL_Generic.h"
!-------------------------------------------------------------------------
!NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1, GEOS/DAS!
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  GmiThermalRC_GCCMod --- GMI Thermal Rate Constants Grid Component Class
!
! Grid Component class for the GMI Thermal Rate Constants
!
! !INTERFACE:
!

   MODULE  GmiThermalRC_GCCMod

! !USES:

   USE ESMF
   USE MAPL
   USE Chem_Mod 	     ! Chemistry Base Class
   USE Chem_UtilMod

   USE Species_BundleMod

   USE GmiChemistryMethod_mod,        ONLY : t_Chemistry
   USE GmiEmissionMethod_mod,         ONLY : t_Emission
   USE GmiDepositionMethod_mod,       ONLY : t_Deposition
   USE GmiSpcConcentrationMethod_mod, ONLY : t_SpeciesConcentration
   USE GmiGrid_mod,                   ONLY : t_gmiGrid
   USE GmiTimeControl_mod,            ONLY : t_GmiClock
   USE GmiESMFrcFileReading_mod,      ONLY : rcEsmfReadTable
   use GmiArrayBundlePointer_mod, only : t_GmiArrayBundle, CleanArrayPointer
   use GmiFieldBundleESMF_mod,        ONLY : obtainTracerFromBundle
   use GmiFieldBundleESMF_mod,        ONLY : updateTracerToBundle
   use GmiFieldBundleESMF_mod,        ONLY : addTracerToBundle
   use GmiSwapSpeciesBundlesMod,      ONLY : SwapSpeciesBundles, speciesReg_for_CCM

   IMPLICIT NONE
   INTEGER, PARAMETER :: DBL = KIND(0.00D+00)

#include "setkin_par.h"
#include "GmiParameters.h"
#include "gmi_phys_constants.h"
#include "setkin_mw.h"
#include "setkin_lchem.h"

! !TYPES:

   PRIVATE
   PUBLIC  GmiThermalRC_GridComp       ! The GMI object

! !PUBLIC MEMBER FUNCTIONS:

   PUBLIC  GmiThermalRC_GridCompInitialize
   PUBLIC  GmiThermalRC_GridCompRun
   PUBLIC  GmiThermalRC_GridCompFinalize

!
! !DESCRIPTION:
!
!  This module implements the GMI Thermal Rate Constants calculations.
!
! !REVISION HISTORY:
!
!  16Sep2003 da Silva  First crack.
!  24Jan2O05 Nielsen   Implementation of Code 916 chemistry
!  19Dec2005 d Silva   Minor portability mods.
!  30Oct2007 Nielsen   GMI Combo set up
!  01Jul2010 Kouatchou Creation of Chemistry Grid Component class
!  01Jun2015 Nielsen   ExtData replaces Chem_UtilMPread
!
!EOP
!-------------------------------------------------------------------------

  TYPE GmiThermalRC_GridComp
   CHARACTER(LEN=255) :: name = "GMI Thermal Rate Constants"

! Does the gmichem_import_rst file exist?  If not, some import states will be
! unfilled, and the species will need to "freewheel" through the first time step.
! -------------------------------------------------------------------------------
   LOGICAL :: gotImportRst

! Various switches
! ----------------
   LOGICAL :: pr_diag
   LOGICAL :: do_synoz
   LOGICAL :: do_qqjk_inchem
   LOGICAL :: do_qqjk_reset
   LOGICAL :: pr_qqjk
   LOGICAL :: do_wetchem
   LOGICAL :: do_AerDust_Calc
   INTEGER :: phot_opt

! Dimensions
! ----------
   INTEGER :: i1, i2, im, j1, j2, jm, km

! Useful character strings
! ------------------------
   CHARACTER(LEN=255) :: chem_mecha
   CHARACTER(LEN=MAX_LENGTH_MET_NAME) :: metdata_name_org
   CHARACTER(LEN=MAX_LENGTH_MET_NAME) :: metdata_name_model

   REAL*8, POINTER     :: rxnr_adjust(:,:,:,:,:)   => null()
   INTEGER	       :: num_rxnr_adjust
   INTEGER	       :: rxnr_adjust_timpyr
   INTEGER, POINTER    :: rxnr_adjust_map(:) => null()
   LOGICAL	       :: do_rxnr_adjust
   CHARACTER (LEN=MAX_LENGTH_FILE_NAME) :: rxnr_adjust_infile_name
   CHARACTER (LEN=MAX_LENGTH_VAR_NAME)  :: rxnr_adjust_var_name

   TYPE (t_GmiArrayBundle), POINTER :: qkgmi(:) => null()

! Extra diagnostics
! -----------------
   LOGICAL :: verbose

! Map GMI species indices to CCM indices
! --------------------------------------
   INTEGER, POINTER :: mapSpecies(:)

! Component derived type declarations
! -----------------------------------
   TYPE(t_gmiGrid   )		:: gmiGrid
   TYPE(t_GmiClock  )           :: gmiClock
   TYPE(t_SpeciesConcentration) :: SpeciesConcentration
 
  END TYPE GmiThermalRC_GridComp

CONTAINS

!-------------------------------------------------------------------------
!NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1, GEOS/DAS!
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GmiThermalRC_GridCompInitialize --- Initialize GmiThermalRC_GridComp
!
! !INTERFACE:
!

   SUBROUTINE GmiThermalRC_GridCompInitialize( self, bgg, bxx, impChem, expChem, nymd, nhms, &
                                      tdt, rc )

   USE GmiSpcConcentrationMethod_mod, ONLY : InitializeSpcConcentration
   USE GmiGrid_mod,		      ONLY : InitializeGmiGrid
   USE GmiTimeControl_mod,            ONLY : Set_curGmiDate, Set_curGmiTime
   USE GmiTimeControl_mod,            ONLY : Set_begGmiDate, Set_begGmiTime
   USE GmiTimeControl_mod,            ONLY : Set_numTimeSteps

   IMPLICIT none
   INTEGER, PARAMETER :: DBL = KIND(0.00D+00)

! !INPUT PARAMETERS:

   TYPE(Species_Bundle), INTENT(in) :: bgg                ! Chemical tracer fields, delp, +
   TYPE(Species_Bundle), INTENT(in) :: bxx                ! Chemical tracer fields, delp, +
   INTEGER, INTENT(IN) :: nymd, nhms		       ! Time from AGCM
   REAL,    INTENT(IN) :: tdt			       ! Chemistry time step (secs)

! !OUTPUT PARAMETERS:

   TYPE(GmiThermalRC_GridComp), INTENT(INOUT)  :: self      ! Grid Component
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

   CHARACTER(LEN=*), PARAMETER :: IAm = 'GmiThermalRC_GridCompInitialize'
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
   INTEGER :: LogicalUnitNum

   INTEGER :: loc_proc, locGlobProc, commu_slaves
   LOGICAL :: one_proc, rootProc
   LOGICAL :: exists,open,found
   
   REAL :: qmin, qmax, tokgCPerBox
   REAL(KIND=DBL) :: tempLook

   real(rPrec), pointer :: var(:,:,:)
   type(ESMF_FieldBundle)      :: qkBundle
   integer                     :: numVars, ib
   character (len=4) :: binName
   character(len=ESMF_MAXSTR) :: varName

   self%name = 'GMI Thermal Rate Constants'

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
         PRINT *,"Starting Reading the GMI Resource File for ThermalRC"
      ENDIF

      gmiConfigFile = ESMF_ConfigCreate(rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigLoadFile(gmiConfigFile, TRIM(rcfilen), rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(gmiConfigFile, importRestartFile, &
     &                label   = "importRestartFile:", &
     &                default = ' ', rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(gmiConfigFile, self%chem_mecha, &
     &                label   = "chem_mecha:", &
     &                default = 'strat_trop', rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(gmiConfigFile, value=self%do_synoz, &
     &           label="do_synoz:", default=.false., rc=STATUS)
      VERIFY_(STATUS)

      !------------------------------
      ! Diagnostics related variables
      !------------------------------

      call ESMF_ConfigGetAttribute(gmiConfigFile, value=self%pr_diag, &
     &           label="pr_diag:", default=.false., rc=STATUS)
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(gmiConfigFile, value=self%verbose, &
     &           label="verbose:", default=.false., rc=STATUS)
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(gmiConfigFile, value=self%do_wetchem, &
     &           label="do_wetchem:", default=.false., rc=STATUS)
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(gmiConfigFile, self%metdata_name_org, &
     &                label   = "metdata_name_org:", &
     &                default = 'GMAO', rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(gmiConfigFile, self%metdata_name_model, &
     &                label   = "metdata_name_model:", &
     &                default = 'GEOS-5', rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(gmiConfigFile, self%phot_opt, &
     &                label   = "phot_opt:", &
     &                default = 1, rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(gmiConfigFile, value=self%do_AerDust_Calc, &
     &           label="do_AerDust_Calc:", default=.false., rc=STATUS)
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(gmiConfigFile, value=self%pr_qqjk, &
     &           label="pr_qqjk:", default=.false., rc=STATUS)
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(gmiConfigFile, value=self%do_qqjk_reset, &
     &           label="do_qqjk_reset:", default=.true., rc=STATUS)
      VERIFY_(STATUS)

      !----------------------------
      ! Chemistry Related Variables
      !----------------------------

      call ESMF_ConfigGetAttribute(gmiConfigFile, value=self%do_qqjk_inchem, &
     &           label="do_qqjk_inchem:", default=.false., rc=STATUS)
      VERIFY_(STATUS)

!     -------------------------
!     Reaction rate adjustment:
!     -------------------------

      call ESMF_ConfigGetAttribute(gmiConfigFile, value=self%do_rxnr_adjust, &
     &           label="do_rxnr_adjust:", default=.false., rc=STATUS)

      call ESMF_ConfigGetAttribute(gmiConfigFile, self%rxnr_adjust_infile_name, &
     &                label   = "rxnr_adjust_infile_name:", &
     &                default = '', rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(gmiConfigFile, self%rxnr_adjust_var_name, &
     &                label   = "rxnr_adjust_var_name:", &
     &                default = 'reac_rate_adj', rc=STATUS )
      VERIFY_(STATUS)

      if (self%do_rxnr_adjust) then
         IF(rootProc) THEN
           WRITE(6,*) 'Code not ready for self%do_rxnr_adjust=',self%do_rxnr_adjust
           WRITE(6,*) ' '
         END IF
         STOP
      else
         self%num_rxnr_adjust = 0
      end if


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
   locGlobProc = -99
   commu_slaves = -99
   
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

! Photolysis reaction list.  Read from kinetics
! text file after finding an available logical unit.
! --------------------------------------------------
   found = .FALSE.
   i = 11
   DO WHILE (.NOT. found .AND. i <= 99)
    INQUIRE(UNIT=i,EXIST=exists,OPENED=open)
    IF(exists .AND. .NOT. open) THEN
     found = .TRUE.
     LogicalUnitNum = i
    END IF
    i = i+1
   END DO

! Allocate space, etc., but the initialization of the
! species from the internal state is left to the run method.
! ----------------------------------------------------------

      CALL InitializeSpcConcentration(self%SpeciesConcentration,              &
                     self%gmiGrid, gmiConfigFile, NSP, NMF, NCHEM,            &
                     loc_proc)

      Allocate(self%qkgmi(NUM_K))
      do ic = 1, NUM_K
         Allocate(self%qkgmi(ic)%pArray3D(i1:i2, ju1:j2, k1:k2))
         self%qkgmi(ic)%pArray3D(i1:i2, ju1:j2, k1:k2) = 0.0d0
      end do

   !========================
   ! Initialize the qkBundle
   !========================

   call ESMF_StateGet(expChem, 'gmiQK' , qkBundle,   RC=STATUS)
   VERIFY_(STATUS)

   ! Add tracer to the bundle
   do ib = 1, NUM_K
      allocate( var(i1:i2, j1:j2, 1:km), STAT=STATUS)
      VERIFY_(STATUS)
      var(:,:,:)  = 0.0d0

      write (binName ,'(i4.4)') ib
      varName = 'qk'//binName

      call addTracerToBundle (qkBundle, var, bgg%grid_esmf, varName)
   end do

   ! Sanity check

   call ESMF_FieldBundleGet(qkBundle, fieldCount=numVars , rc=STATUS)
   VERIFY_(STATUS)
   _ASSERT(NUM_K == numVars,'needs informative message')

    !---------------------------------------------------------------
    ! Create and populate the array that maps GMI species indices to
    ! GEOS-5 species indices
    !---------------------------------------------------------------

    allocate(self%mapSpecies(NSP))
    self%mapSpecies(:) = speciesReg_for_CCM(lchemvar, NSP, bgg%reg%vname, bxx%reg%vname )

  RETURN

  END SUBROUTINE GmiThermalRC_GridCompInitialize

!-------------------------------------------------------------------------
!NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1, GEOS/DAS!
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  GmiThermalRC_GridCompRun --- The GMI Driver
!
! !INTERFACE:
!

   SUBROUTINE GmiThermalRC_GridCompRun ( self, bgg, bxx, impChem, expChem, nymd, nhms, &
                                tdt, rc )

! !USES:

   USE GmiTimeControl_mod,            ONLY : Set_curGmiDate, Set_curGmiTime
   USE GmiTimeControl_mod,            ONLY : Set_numTimeSteps, Get_numTimeSteps
   USE GmiTimeControl_mod,            ONLY : Set_gmiSeconds, GetSecondsFromJanuary1
   USE GmiSpcConcentrationMethod_mod, ONLY : resetFixedConcentration
   USE GmiSolar_mod,                  ONLY : CalcCosSolarZenithAngle
   use GmiThermalRateConstants_mod  , only : calcThermalRateConstants

   IMPLICIT none

! !INPUT/OUTPUT PARAMETERS:

   TYPE(GmiThermalRC_GridComp), INTENT(INOUT) :: self ! Grid Component
   TYPE(Species_Bundle), INTENT(INOUT) :: bgg    ! Chemical tracer fields   
   TYPE(Species_Bundle), INTENT(INOUT) :: bxx    ! Chemical tracer fields   

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

   CHARACTER(LEN=*), PARAMETER :: IAm = 'GmiThermalRC_GridCompRun'

!  Input fields from GEOS-5
!  ------------------------
   REAL, POINTER, DIMENSION(:,:) :: zpbl
   REAL, POINTER, DIMENSION(:,:) :: frland, frlandice, asnow

   REAL, POINTER, DIMENSION(:,:,:) :: ple, Q, T, zle
   REAL, POINTER, DIMENSION(:,:,:) :: ql,cnv_mfc
   REAL, POINTER, DIMENSION(:,:,:) :: rh2

!  Exports not part of internal state
!  ----------------------------------
   REAL, POINTER, DIMENSION(:,:,:) :: HO2PBLFlag

!  Local
!  -----
   INTEGER :: cymd, dymd, emiss_opt, hms
   INTEGER :: i, i1, i2, ic, im, iXj, iTile(1)
   INTEGER :: j, j1, j2, jm, jTile(1)
   INTEGER :: k, km, kReverse
   INTEGER :: i1_gl, i2_gl, ju1_gl, j2_gl, ilo, gmi_nborder
   INTEGER :: ihi, julo, jhi, ju1,  k1, k2, ilong, ilat, ivert
   INTEGER :: loc_proc
   INTEGER :: num_time_steps
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

   REAL(KIND=DBL) :: dehydmin = 0.00
   REAL(KIND=DBL) :: chemDt, dayOfYear

   CHARACTER(LEN=255) :: speciesName
   CHARACTER(LEN=255) :: importName

   LOGICAL :: found, rootProc
   LOGICAL, PARAMETER :: doThis = .FALSE.
   
! Allocatables.  Use KIND=DBL where GMI expects REAL*8.
! -----------------------------------------------------
   INTEGER, ALLOCATABLE :: conPBLFlag(:,:,:)

   REAL, ALLOCATABLE :: var3d(:,:,:)
   REAL, ALLOCATABLE :: pl(:,:,:)

   REAL(KIND=DBL), ALLOCATABLE :: tropopausePress(:,:)

   REAL(KIND=DBL), ALLOCATABLE :: press3c(:,:,:)
   REAL(KIND=DBL), ALLOCATABLE :: kel(:,:,:)
   REAL(KIND=DBL), ALLOCATABLE :: clwc(:,:,:)
   REAL(KIND=DBL), ALLOCATABLE :: cmf(:,:,:)
   REAL(KIND=DBL), ALLOCATABLE :: relativeHumidity(:,:,:)

   REAL(KIND=DBL), ALLOCATABLE :: tArea  (:,:,:,:)
   REAL(KIND=DBL), ALLOCATABLE :: eRadius(:,:,:,:)
   type (t_GmiArrayBundle), pointer :: gmiSAD(:) => null()

   loc_proc = -99

!  Grid specs
!  ----------
   rc = 0
   i1 = self%i1
   i2 = self%i2
   im = self%im
   
   j1 = self%j1
   j2 = self%j2
   jm = self%jm
   
   km = self%km
   
   iXj = (i2-i1+1)*(j2-j1+1)

   gmi_nborder = 0
   i1_gl  = 1
   i2_gl  = i2
   ju1_gl = 1
   j2_gl  = j2
   ju1    = j1
   k1     = 1
   k2     = km
   ilo    = i1  - gmi_nborder
   ihi    = i2  + gmi_nborder
   julo   = ju1 - gmi_nborder
   jhi    = j2  + gmi_nborder
   ilong  = i2 - i1  + 1
   ilat   = j2 - ju1 + 1
   ivert  = k2 - k1 + 1

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
   ALLOCATE(    tropopausePress(i1:i2,j1:j2),STAT=STATUS)
   VERIFY_(STATUS)

   ALLOCATE(                pl(i1:i2,j1:j2,1:km),STAT=STATUS)
   VERIFY_(STATUS)
   ALLOCATE(             var3d(i1:i2,j1:j2,1:km),STAT=STATUS)
   VERIFY_(STATUS)
   ALLOCATE(           press3c(i1:i2,j1:j2,1:km),STAT=STATUS)
   VERIFY_(STATUS)
   ALLOCATE(               kel(i1:i2,j1:j2,1:km),STAT=STATUS)
   VERIFY_(STATUS)
   ALLOCATE(              clwc(i1:i2,j1:j2,1:km),STAT=STATUS)
   VERIFY_(STATUS)
   ALLOCATE(               cmf(i1:i2,j1:j2,1:km),STAT=STATUS)
   VERIFY_(STATUS)
   ALLOCATE(  relativeHumidity(i1:i2,j1:j2,1:km),STAT=STATUS)
   VERIFY_(STATUS)
   ALLOCATE(        conPBLflag(i1:i2,j1:j2,1:km),STAT=STATUS)
   VERIFY_(STATUS)

   !---------------------------------
   ! Obtain data from the ESMF Bundle
   !---------------------------------

   call obtainBundles(expChem)

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
! --------------------------------------------------------

   CALL Acquire_Clims(STATUS)
   VERIFY_(STATUS)

! Diagnostics capabilities enabled?
! ---------------------------------
   IF(self%pr_qqjk .AND. self%do_qqjk_inchem) THEN
    IF( MAPL_AM_I_ROOT() ) THEN
     PRINT *,TRIM(IAm),': Code not ready for pr_qqjk=',self%pr_qqjk,' and do_qqjk_inchem=',self%do_qqjk_inchem
    END IF
    rc = 61
    RETURN
   END IF

! Grab imports and do units conversions
! -------------------------------------

   CALL SatisfyImports(STATUS)
   VERIFY_(STATUS)

! Hand the species concentrations to GMI's bundle
! -----------------------------------------------
   IF (self%gotImportRst) then
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
! GMI Thermal Rate Constants
! ------------------------------------------------------------------------

      IF (self%gotImportRst) THEN

         !--------------------------------------------------------
         ! Calculate the air density at the center of each grid box
         ! (molecules/cm^3).
         ! NOTE: BOLTZMN_E units = erg/K/molec
         !--------------------------------------------------------

          self%SpeciesConcentration%concentration(IMGAS)%pArray3D(:,:,:) =  &
     &                press3c(i1:i2,ju1:j2,:) * MB2CGS /  &
     &                (kel (i1:i2,ju1:j2,:) * BOLTZMN_E)

          self%SpeciesConcentration%concentration(IOXYGEN)%pArray3D(:,:,:) =  &
     &        self%SpeciesConcentration%concentration(IMGAS)%pArray3D(:,:,:) * MXRO2

          self%SpeciesConcentration%concentration(INITROGEN)%pArray3D(:,:,:) =  &
     &        self%SpeciesConcentration%concentration(IMGAS)%pArray3D(:,:,:) * MXRN2

          !===========================================
          ! Fixes for H2 - Provided by David Considine
          !===========================================
          self%SpeciesConcentration%concentration(IH2)%pArray3D(:,:,:) = MXRH2
          !=================
          ! end fixes for H2
          !=================

          CALL Get_numTimeSteps(self%gmiClock, num_time_steps)

          call calcThermalRateConstants (self%do_wetchem, self%chem_mecha,      &
     &             rootProc, num_time_steps, IH2O, IMGAS, nymd,                 &
     &             self%rxnr_adjust_map,     &
     &             press3c, tropopausePress, kel, clwc, cmf, gmiSAD,    &
     &             self%qkgmi, self%SpeciesConcentration%concentration,         &
     &             self%rxnr_adjust, eRadius, tArea, relativeHumidity,          &
     &             conPBLFlag, self%do_AerDust_Calc, self%phot_opt,             &
     &             self%pr_diag, loc_proc, self%num_rxnr_adjust,                &
     &             self%rxnr_adjust_timpyr, ivert, NSAD, NUM_K, NMF, NSP,       &
     &             ilo, ihi, julo, jhi, i1, i2, ju1, j2, k1, k2)
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

   CALL FillExports(rc)

   CALL populateBundleQK()

! Scratch local work space
! ------------------------
   DEALLOCATE(tropopausePress, STAT=STATUS)
   VERIFY_(STATUS)

   DEALLOCATE(pl, press3c, kel, clwc, cmf, relativeHumidity, &
              conPBLFlag, STAT=STATUS)
   VERIFY_(STATUS)

   DEALLOCATE(tArea, eRadius, STAT=STATUS)
   VERIFY_(STATUS)

   CALL CleanArrayPointer(gmiSAD, STATUS)
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
!  Obtain the following climatologies from files supplied by
!  the GMI project:
!   Fixed concentration species
!   Stratospheric sulfate surface area
!   Emissions
!
!EOP
!---------------------------------------------------------------------------

  CHARACTER(LEN=255) :: IAm
  INTEGER :: lightning_opt
  INTEGER, ALLOCATABLE :: flag(:)
  
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

  CHARACTER(LEN=255) :: IAm
  
  rc=0
  IAm="FillExports"

! This is the continental PBL flag (0 or 1) that is used to compute the
! loss of HO2 (Perhydroxyl radical) in aerosols.  See setkin_kcalc.F90.
! ---------------------------------------------------------------------
   IF(ASSOCIATED(HO2PBLFlag)) THEN
    DO k=1,km
     kReverse = km-k+1
     HO2PBLFlag(i1:i2,j1:j2,kReverse) = conPBLFlag(i1:i2,j1:j2,k)
    END DO
   END IF

  RETURN
 END SUBROUTINE FillExports

!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  populateBundleQK
!
! !INTERFACE:

      subroutine populateBundleQK()
!
      implicit none
!
! !DESCRIPTION:
!
! !LOCAL VARIABLES:
      integer :: STATUS, numVars, ib, rc
      real(rPrec), pointer, dimension(:,:,:)   :: ptr3D
      type(ESMF_FieldBundle)                ::      qkBundle
      character(len=ESMF_MAXSTR), parameter :: IAm = "populateBundleQK"
!
!EOP
!--------------------------------------------------------------------------------

      allocate(ptr3D(i1:i2, j1:j2, 1:km))

      !==================
      ! For the QK Bundle
      !==================

      call ESMF_StateGet(expChem, "gmiQK", qkBundle, rc=STATUS)
      VERIFY_(STATUS)

      call ESMF_FieldBundleGet(qkBundle, fieldCount=numVars, rc=STATUS)
      VERIFY_(STATUS)
      _ASSERT(numVars == NUM_K,'needs informative message')

      do ib = 1, numVars
         ptr3D(:,:,:) = self%qkgmi(ib)%pArray3D(:,:,km:1:-1)
         call updateTracerToBundle(qkBundle, ptr3D, ib)
      end do

      deallocate (ptr3D)

      return

      end subroutine populateBundleQK

!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  obtainBundles
!
! !INTERFACE:

  SUBROUTINE obtainBundles(state)
!
   USE ESMF

  implicit none
!
! !INPUT PARAMETERS:
  type(ESMF_State) , intent(in) :: state

!
! !DESCRIPTION:
! Extract arrays from the export state.
!
! !LOCAL VARIABLES:

      integer :: ib, numVars, rc
      type(ESMF_FieldBundle)      ::     sadBundle
      type(ESMF_FieldBundle)      ::   tAreaBundle
      type(ESMF_FieldBundle)      :: eRadiusBundle
      real(rPrec), pointer        :: ptr3D(:,:,:)

      integer :: STATUS
      character(len=ESMF_MAXSTR), parameter :: IAm = "obtainBundless"
!
!EOP
!-------------------------------------------------------------------------------

      !===============
      ! the SAD Bundle
      !===============

      call ESMF_StateGet (state, "gmiSAD", sadBundle, RC=STATUS )
      VERIFY_(STATUS)

      call ESMF_FieldBundleGet(sadBundle, fieldCount=numVars , rc=STATUS)
      VERIFY_(STATUS)

      ALLOCATE(gmiSAD(numVars),STAT=STATUS)

      do ib = 1, numVars
         CALL obtainTracerFromBundle(sadBundle, ptr3D, ib)

         Allocate(gmiSAD(ib)%pArray3D(i1:i2, j1:j2, 1:km))
         gmiSAD(ib)%pArray3D(:,:,km:1:-1) = ptr3D(:,:,:)
      end do

      !=================
      ! the tArea Bundle
      !=================

      call ESMF_StateGet (state, "gmiTAREA", tAreaBundle, RC=STATUS )
      VERIFY_(STATUS)

      call ESMF_FieldBundleGet(tAreaBundle, fieldCount=numVars , RC=STATUS)
      VERIFY_(STATUS)

      ALLOCATE(tArea(i1:i2, j1:j2, 1:km, numVars),STAT=STATUS)
      VERIFY_(STATUS)

      do ib = 1, numVars
         CALL obtainTracerFromBundle(tAreaBundle, ptr3D, ib)
         tArea(:,:,km:1:-1,ib) = ptr3D(:,:,:)
      end do

      !===================
      ! the eRadius Bundle
      !===================

      call ESMF_StateGet (state, "gmiERADIUS", eRadiusBundle, RC=STATUS )
      VERIFY_(STATUS)

      call ESMF_FieldBundleGet(eRadiusBundle, fieldCount=numVars , RC=STATUS)
      VERIFY_(STATUS)

      ALLOCATE(eRadius(i1:i2, j1:j2, 1:km, numVars),STAT=STATUS)
      VERIFY_(STATUS)

      do ib = 1, numVars
         CALL obtainTracerFromBundle(eRadiusBundle, ptr3D, ib)
         eRadius(:,:,km:1:-1,ib) = ptr3D(:,:,:)
      end do

      RETURN

  END SUBROUTINE obtainBundles

!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  FindPointers
!
! !INTERFACE:

  SUBROUTINE FindPointers(rc)
  
   USE ESMF

  IMPLICIT NONE

  INTEGER, INTENT(OUT) :: rc
!
! !DESCRIPTION:
!
!  Find pointers to import and export states
!
!EOP
!---------------------------------------------------------------------------

  CHARACTER(LEN=255) :: IAm
  REAL :: qmin,qmax
  
  rc=0
  IAm="FindPointers"

!  Pointers to imports
!  -------------------
   CALL MAPL_GetPointer(impChem,      zpbl,      'ZPBL', RC=STATUS)
   VERIFY_(STATUS)
   CALL MAPL_GetPointer(impChem,    frland,    'FRLAND', RC=STATUS)
   VERIFY_(STATUS)
   CALL MAPL_GetPointer(impChem, frlandice, 'FRLANDICE', RC=STATUS)
   VERIFY_(STATUS)
   CALL MAPL_GetPointer(impChem,     asnow,     'ASNOW', RC=STATUS)
   VERIFY_(STATUS)
  

   CALL MAPL_GetPointer(impChem,       ple,	'PLE', RC=STATUS)
   VERIFY_(STATUS)
   CALL MAPL_GetPointer(impChem,         Q,       'Q', RC=STATUS)
   VERIFY_(STATUS)
   CALL MAPL_GetPointer(impChem,	 T,	  'T', RC=STATUS)
   VERIFY_(STATUS)
   CALL MAPL_GetPointer(impChem,       zle,	'ZLE', RC=STATUS)
   VERIFY_(STATUS)
   CALL MAPL_GetPointer(impChem,	ql,	 'QL', RC=STATUS)
   VERIFY_(STATUS)
   CALL MAPL_GetPointer(impChem,   cnv_mfc, 'CNV_MFC', RC=STATUS)
   VERIFY_(STATUS)
   CALL MAPL_GetPointer(impChem,       rh2,	'RH2', RC=STATUS)
   VERIFY_(STATUS)

!  Export state pointers
!  ---------------------
   CALL MAPL_GetPointer(expChem,   HO2PBLFlag,  'HO2PBLFLAG', RC=STATUS)
   
!  Validation
!  ----------
   Validate: IF(self%verbose) THEN
    IF(MAPL_AM_I_ROOT( ))PRINT *,TRIM(IAm),": Input ..."
    CALL pmaxmin('ZPBL:', zpbl, qmin, qmax, iXj, 1, 1. )
    CALL pmaxmin('Q:', Q, qmin, qmax, iXj, km, 1. )
    CALL pmaxmin('T:', T, qmin, qmax, iXj, km, 1. )
    CALL pmaxmin('ZLE:', zle, qmin, qmax, iXj, km+1, 1. )
    CALL pmaxmin('PLE (hPa):', ple, qmin, qmax, iXj, km+1, 0.01 )
    CALL pmaxmin('QL:', ql, qmin, qmax, iXj, km, 1. )
    CALL pmaxmin('CNV_MFC:', cnv_mfc, qmin, qmax, iXj, km+1, 1. )
    CALL pmaxmin('RH2:', rh2, qmin, qmax, iXj, km, 1. )

    i = bxx%reg%nq
    CALL pmaxmin('TROPP:', bxx%qa(i)%data3d(:,:,km), qmin, qmax, iXj, 1, 0.01 )

    CALL pmaxmin('FRLAND:', frland, qmin, qmax, iXj, 1, 1. )
    CALL pmaxmin('FRLANDICE:', frlandice, qmin, qmax, iXj, 1, 1. )
    CALL pmaxmin('ASNOW:', asnow, qmin, qmax, iXj, 1, 1. )

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
! The most recent valid tropopause pressures are stored in T2M15D(:,:,km)
! -----------------------------------------------------------------------   ------------       -------------
  i = bxx%reg%nq
  tropopausePress(i1:i2,j1:j2) = bxx%qa(i)%data3d(i1:i2,j1:j2,km)*Pa2hPa    ! Pa               hPa

! Layer means                                                               GEOS-5 Units       GMI Units
! -----------                                                               ------------       -------------
  DO k=1,km
   kReverse = km-k+1                                                        ! Lid-to-surf      Surf-to-lid
   press3c(i1:i2,j1:j2,kReverse) = pl(i1:i2,j1:j2,k)*Pa2hPa                 ! Pa               hPa
   kel(i1:i2,j1:j2,kReverse) = T(i1:i2,j1:j2,k)                             ! K
   clwc(i1:i2,j1:j2,kReverse) = ql(i1:i2,j1:j2,k)*ToGrPerKg                 ! kg kg^{-1}       g kg^{-1}
   relativeHumidity(i1:i2,j1:j2,kReverse) = rh2(i1:i2,j1:j2,k)              ! fraction

  END DO

! These bounds are in Jules' RH code
! ----------------------------------
  WHERE(relativeHumidity < 0.00D+00) relativeHumidity = 0.00D+00
  WHERE(relativeHumidity > 0.95D+00) relativeHumidity = 0.95D+00

! Layer edges                                                               GEOS-5 Units       GMI Units
! -----------                                                               ------------       -------------
  DO k=0,km-1
   kReverse = km-k
   cmf(i1:i2,j1:j2,kReverse) = cnv_mfc(i1:i2,j1:j2,k)                       ! kg m^{-2}s^{-1}
  END DO

! Continental boundary layer flag for H02 net uptake reaction in setkin_kcalc.F90
! -------------------------------------------------------------------------------
   ALLOCATE(on(i1:i2,j1:j2))
   conPBLFlag(i1:i2,j1:j2,1:km) = 0

! Work up from the ground. Note: ZLE is referenced to MSL, not the ground!
! ------------------------------------------------------------------------
   DO k=km,1,-1
    kReverse = km-k+1

    on(i1:i2,j1:j2) = 0
    WHERE(zpbl(i1:i2,j1:j2) > (zle(i1:i2,j1:j2,k)-zle(i1:i2,j1:j2,km)) ) on(i1:i2,j1:j2) = 1

! Done when no more cells are in the PBL
! --------------------------------------
    n = SUM(on)
    IF(n == 0) EXIT

! Raise the flag when in the PBL and when more than half of the surface is  
! land, of which less than half of cell land surface covered with snow or ice.
! ----------------------------------------------------------------------------
    WHERE( on(i1:i2,j1:j2) == 1 .AND. frland(i1:i2,j1:j2) >= 0.50 .AND.   &
           (frlandice(i1:i2,j1:j2) < 0.5 .OR. asnow(i1:i2,j1:j2) < 0.5)  ) &
          conPBLFlag(i1:i2,j1:j2,kReverse) = 1
   END DO

   DEALLOCATE(on)

  RETURN
 END SUBROUTINE SatisfyImports

 END SUBROUTINE GmiThermalRC_GridCompRun

!-------------------------------------------------------------------------
!NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1, GEOS/DAS!
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  GmiThermalRC_GridCompFinalize
!
! !INTERFACE:
!

   SUBROUTINE GmiThermalRC_GridCompFinalize ( self, impChem, expChem, &
                                     nymd, nhms, cdt, rc )

  IMPLICIT none

! !INPUT/OUTPUT PARAMETERS:

   TYPE(GmiThermalRC_GridComp), INTENT(inout) :: self ! Grid Component

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

   CHARACTER(LEN=*), PARAMETER :: IAm = 'GmiThermalRC_GridCompFinalize'
   
   rc=0

   RETURN

 END SUBROUTINE GmiThermalRC_GridCompFinalize
  
 END MODULE GmiThermalRC_GCCMod

