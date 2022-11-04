#include "MAPL_Generic.h"
!-------------------------------------------------------------------------
!NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1, GEOS/DAS!
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  GmiEmiss_GCCMod --- GMI Grid Component Class
!
! Grid Component class for the GMI combined stratopshere/troposphere 
! chemistry
!
! !INTERFACE:
!

   MODULE  GmiEmiss_GCCMod

! !USES:

   USE ESMF
   USE MAPL
   USE Chem_Mod
   USE Chem_UtilMod

   USE Species_BundleMod

   USE GmiChemistryMethod_mod,        ONLY : t_Chemistry
   USE GmiEmissionMethod_mod,         ONLY : t_Emission
   USE GmiDepositionMethod_mod,       ONLY : t_Deposition
   USE GmiSpcConcentrationMethod_mod, ONLY : t_SpeciesConcentration
   USE GmiGrid_mod,                   ONLY : t_gmiGrid
   USE GmiTimeControl_mod,            ONLY : t_GmiClock
   USE GmiESMFrcFileReading_mod,      ONLY : rcEsmfReadTable
   USE GmiArrayBundlePointer_mod,     ONLY : t_GmiArrayBundle, CleanArrayPointer
   USE GmiFieldBundleESMF_mod,        ONLY : updateTracerToBundle
   USE GmiFieldBundleESMF_mod,        ONLY : addTracerToBundle
   USE GmiFieldBundleESMF_mod,        ONLY : obtainTracerFromBundle
   USE GmiShipEmission_mod,           ONLY : calcShipEmission
   USE GmiStateFieldESMF_mod,	      ONLY : getDataFromStateField
   USE GmiSurfaceEmissionInChemistry_mod, ONLY : updateSurfEmissionInChemistry
   USE GmiSwapSpeciesBundlesMod,      ONLY : SwapSpeciesBundles, speciesReg_for_CCM
   USE VegLaiMod,                     ONLY : Decode_Land_Types, Decode_XLAI
   USE OVP,                           ONLY : OVP_init, OVP_end_of_timestep_hms, OVP_mask, OVP_apply_mask

   IMPLICIT NONE
   INTEGER, PARAMETER :: DBL = KIND(0.00D+00)

#include "setkin_par.h"
#include "GmiParameters.h"
#include "gmi_phys_constants.h"
#include "gmi_emiss_constants.h"
#include "setkin_mw.h"
#include "setkin_lchem.h"
#include "gmi_sad_constants.h"

! !TYPES:

   PRIVATE

   INTEGER, SAVE, ALLOCATABLE :: MASK_10AM(:,:)
   INTEGER, SAVE, ALLOCATABLE :: MASK_2PM(:,:)
   INTEGER, SAVE              :: OVP_FIRST_HMS
   INTEGER, SAVE              :: OVP_RUN_DT
   INTEGER, SAVE              :: OVP_GC_DT
   INTEGER, SAVE              :: OVP_MASK_DT
   CHARACTER(LEN=ESMF_MAXSTR), save :: ptfile_save(MAX_NUM_CONST) = '/dev/null'

   PUBLIC  GmiEmiss_GridComp       ! The GMI object 
   PUBLIC  t_GmiPointEmiss         ! GMI point emission data type 

! !PUBLIC MEMBER FUNCTIONS:

   PUBLIC  GmiEmiss_GridCompInitialize, GmiEmiss_initSurfEmissBundle
   PUBLIC  GmiEmiss_GridCompRun
   PUBLIC  GmiEmiss_GridCompFinalize

!
! !DESCRIPTION:
!
!  This module implements the GMI combined stratopshere/troposphere
!  chemistry
!
! !REVISION HISTORY:
!
!  16Sep2003 da Silva  First crack.
!  24Jan2O05 Nielsen   Implementation of Code 916 chemistry
!  19Dec2005 d Silva   Minor portability mods.
!  30Oct2007 Nielsen   GMI Combo set up
!  01Jul2010 Kouatchou Created this Emission class
!  01Jun2015 Nielsen   ExtData replaces Chem_UtilMPread
!
!EOP
!-------------------------------------------------------------------------
!
!   private :: t_GmiPointEmiss
   TYPE t_GmiPointEmiss
     integer :: nPts
     integer, allocatable :: vStart(:)
     integer, allocatable :: vEnd(:)
     integer, allocatable :: iPoint(:)
     integer, allocatable :: jPoint(:)
     real, allocatable :: vLat(:)
     real, allocatable :: vLon(:)
     real, allocatable :: vBase(:)
     real, allocatable :: vTop(:)
     real, allocatable :: vEmis(:)
   END TYPE t_GmiPointEmiss

   TYPE GmiEmiss_GridComp
     CHARACTER(LEN=255) :: name = "GMI Stratospheric/Tropospheric Chemistry"

! Heartbeat (seconds)
! -------------------
     INTEGER :: heartBeat

! Is the GCM executing the predictor step of REPLAY at this moment?
! -----------------------------------------------------------------
     LOGICAL :: doingPredictorNow

! Does the gmichem_import_rst file exist?  If not, some import states will be
! unfilled, and the species will need to "freewheel" through the first time step.
! -------------------------------------------------------------------------------
     LOGICAL :: gotImportRst

! Set BCRealTime = .TRUE. when boundary conditions 
! must be for exact year of current calendar date.
! -------------------------------------------------
     LOGICAL :: BCRealTime

! Daily and monthly emissions
! ---------------------------
     INTEGER :: num_diurnal_emiss

! Various switches
! ----------------
     LOGICAL :: pr_dry_depos
     LOGICAL :: pr_wet_depos   
     LOGICAL :: pr_scav
     LOGICAL :: pr_diag
     LOGICAL :: do_drydep
     LOGICAL :: do_wetdep
     LOGICAL :: do_emission             
     LOGICAL :: do_aerocom             
     LOGICAL :: pr_const
     LOGICAL :: do_synoz
     LOGICAL :: do_gcr
     LOGICAL :: pr_surf_emiss
     LOGICAL :: pr_emiss_3d
     LOGICAL :: do_qqjk_inchem
     LOGICAL :: do_qqjk_reset
     LOGICAL :: pr_qqjk

     integer :: met_opt
     integer :: chem_opt
     integer :: trans_opt
     integer :: ic_NO_lgt

     character (len=MAX_LENGTH_MET_NAME) :: metdata_name_org
     character (len=MAX_LENGTH_MET_NAME) :: metdata_name_model

! Dimensions
! ----------
     INTEGER :: i1, i2, im, j1, j2, jm, km

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

! Information for the export states of emissions [EM_NO, for example].
! --------------------------------------------------------------------
     INTEGER :: numEM_Exports
     CHARACTER(LEN=ESMF_MAXSTR), ALLOCATABLE :: EM_ExportNames(:)
     CHARACTER(LEN=ESMF_MAXSTR), ALLOCATABLE :: EM_ExportUnits(:)

! Component derived type declarations
! -----------------------------------
     TYPE(t_Emission  )           :: Emission
     TYPE(t_Deposition)           :: Deposition
     TYPE(t_gmiGrid   )           :: gmiGrid
     TYPE(t_GmiClock  )           :: gmiClock
     TYPE(t_SpeciesConcentration) :: SpeciesConcentration

     integer :: ship_thisRecord
     integer :: ship_curRecord

! Veg Fraction does not change, only read once
! --------------------------------------------
     LOGICAL :: veg_fraction_done
!
! Point mission data
! --------------------------------------------
     TYPE(t_GmiPointEmiss), pointer :: GmiPointEmiss(:) => null()
 
   END TYPE GmiEmiss_GridComp

  CONTAINS

!-------------------------------------------------------------------------
!NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1, GEOS/DAS!
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  GmiEmiss_GridCompInitialize --- Initialize GmiEmiss_GridComp
!
! !INTERFACE:
!

   SUBROUTINE GmiEmiss_GridCompInitialize( self, bgg, bxx, impChem, expChem, nymd, nhms, &
                                      tdt, gc, clock, rc )

   USE GmiSpcConcentrationMethod_mod, ONLY : InitializeSpcConcentration
   USE GmiEmissionMethod_mod,         ONLY : InitializeEmission, initReadEmission
   USE GmiEmissionMethod_mod,         ONLY : Get_num_emiss, Get_lightning_opt
   USE GmiGrid_mod,                   ONLY : InitializeGmiGrid
   USE GmiTimeControl_mod,            ONLY : Set_curGmiDate, Set_curGmiTime
   USE GmiTimeControl_mod,            ONLY : Set_begGmiDate, Set_begGmiTime
   USE GmiTimeControl_mod,            ONLY : Set_numTimeSteps
   USE gcr_mod,                       ONLY : INIT_GCR_DIAG

   IMPLICIT none
   INTEGER, PARAMETER :: DBL = KIND(0.00D+00)

! !INPUT PARAMETERS:

   TYPE(Species_Bundle), INTENT(IN) :: bgg             ! GMI Species - transported
   TYPE(Species_Bundle), INTENT(IN) :: bxx             ! GMI Species - not transported
   INTEGER,              INTENT(IN) :: nymd, nhms      ! Time from AGCM
   REAL,                 INTENT(IN) :: tdt             ! Chemistry time step (secs)

! !OUTPUT PARAMETERS:

   TYPE(GmiEmiss_GridComp), INTENT(INOUT)  :: self      ! Grid Component
   TYPE(ESMF_State),   INTENT(INOUT)  :: impChem    ! Import State
   TYPE(ESMF_State),   INTENT(INOUT)  :: expChem    ! Export State
   type(ESMF_GridComp), intent(inout)  :: gc      ! Grid Component
   type(ESMF_Clock),  intent(inout) :: clock      ! The clock

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

   CHARACTER(LEN=*), PARAMETER :: IAm = 'GmiEmiss_GridCompInitialize'
   CHARACTER(LEN=255) :: rcfilen = 'GMI_GridComp.rc'
   CHARACTER(LEN=255) :: namelistFile
   CHARACTER(LEN=255) :: importRestartFile
   CHARACTER(LEN=255) :: string, fieldName
   
   CHARACTER(LEN=ESMF_MAXSTR), POINTER   :: itemNames(:)
   TYPE(ESMF_StateItem_Flag), POINTER :: itemTypes(:)
   TYPE(ESMF_Config)  :: gmiConfigFile
   TYPE(ESMF_Field)   :: FIELD

   INTEGER :: i, i1, i2, ic, im, j, j1, j2, jm, k, km, kReverse
   INTEGER :: m, n, scanNumber, STATUS

   INTEGER :: i1_gl, i2_gl, ju1_gl, j2_gl 
   INTEGER :: ju1, jv1, jv1_gl, j1p, j2p
   INTEGER :: k1, k2, k1_gl, k2_gl
   INTEGER :: ilong, ilat, ivert, itloop
   INTEGER :: NPIJ, NPI, NPJ
   INTEGER :: ilo, ihi, julo, jvlo, jhi
   INTEGER :: ilo_gl, ihi_gl, julo_gl, jvlo_gl, jhi_gl
   INTEGER :: gmi_nborder
   INTEGER :: NMR      ! number of species from the GMI_Mech_Registry.rc
   INTEGER :: lightning_opt, LogicalUnitNum
   INTEGER :: num_emiss

   INTEGER :: loc_proc, locGlobProc, commu_slaves
   LOGICAL :: one_proc, rootProc
   LOGICAL :: exists,open,found
   
! Grid cell area can be set by initialize
! ---------------------------------------
   REAL, POINTER, DIMENSION(:,:) :: cellArea

! Work space
! ----------
   REAL, ALLOCATABLE :: var2D(:,:)
   REAL, ALLOCATABLE :: var3D(:,:,:)

   real(rPrec), pointer :: var(:,:)
   type(ESMF_FieldBundle)      :: surfEmissBundle
   integer                     :: numVars, ib
   integer                     :: nTimes, begTime, incSecs
   character (len=4) :: binName
   character (len=15) :: vegName
   character(len=ESMF_MAXSTR) :: varName, speciesName, fileName

   REAL, POINTER, DIMENSION(:,:)   :: LONS

   self%name = 'GMI combined Stratosphere/Troposphere Chemistry'

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
     PRINT *,"Starting Reading the GMI Emission Resource File"
   ENDIF

   gmiConfigFile = ESMF_ConfigCreate(rc=STATUS )
   VERIFY_(STATUS)

   call ESMF_ConfigLoadFile(gmiConfigFile, TRIM(rcfilen), rc=STATUS )
   VERIFY_(STATUS)

   call ESMF_ConfigGetAttribute(gmiConfigFile, importRestartFile, &
                   label="importRestartFile:", default = ' ', rc=STATUS )
   VERIFY_(STATUS)

!------------------------------
! Deposition related variables
!------------------------------
    
    call ESMF_ConfigGetAttribute(gmiConfigFile, value=self%do_drydep, &
              label="do_drydep:", default=.false., rc=STATUS)
    VERIFY_(STATUS)

    call ESMF_ConfigGetAttribute(gmiConfigFile, value=self%do_wetdep, &
              label="do_wetdep:", default=.false., rc=STATUS)
    VERIFY_(STATUS)

!------------------------------
! Emission related variables
!------------------------------

    call ESMF_ConfigGetAttribute(gmiConfigFile, value=self%do_emission, &
              label="do_emission:", default=.false., rc=STATUS)
    VERIFY_(STATUS)

    call ESMF_ConfigGetAttribute(gmiConfigFile, value=self%do_synoz, &
              label="do_synoz:", default=.false., rc=STATUS)
    VERIFY_(STATUS)

    CALL ESMF_ConfigGetAttribute(gmiConfigFile, self%num_diurnal_emiss, &
              LABEL="Diurnal_Emission_Species:", DEFAULT=0, RC=STATUS)
    VERIFY_(STATUS)

    CALL ESMF_ConfigGetAttribute(gmiConfigFile, self%do_gcr, &
                  LABEL="do_gcr:", DEFAULT=.FALSE., RC=STATUS)
    VERIFY_(STATUS)

!------------------------------
! Diagnostics related variables
!------------------------------

    call ESMF_ConfigGetAttribute(gmiConfigFile, value=self%pr_diag, &
              label="pr_diag:", default=.false., rc=STATUS)
    VERIFY_(STATUS)

    call ESMF_ConfigGetAttribute(gmiConfigFile, value=self%verbose, &
              label="verbose:", default=.false., rc=STATUS)
    VERIFY_(STATUS)

    call ESMF_ConfigGetAttribute(gmiConfigFile, value=self%pr_surf_emiss, &
              label="pr_surf_emiss:", default=.false., rc=STATUS)
    VERIFY_(STATUS)

    call ESMF_ConfigGetAttribute(gmiConfigFile, value=self%pr_emiss_3d, &
              label="pr_emiss_3d:", default=.false., rc=STATUS)
    VERIFY_(STATUS)

    call ESMF_ConfigGetAttribute(gmiConfigFile, value=self%pr_const, &
              label="pr_const:", default=.false., rc=STATUS )
    VERIFY_(STATUS)

    call ESMF_ConfigGetAttribute(gmiConfigFile, value=self%do_aerocom, &
              label="do_aerocom:", default=.false., rc=STATUS )
    VERIFY_(STATUS)

    call ESMF_ConfigGetAttribute(gmiConfigFile, self%metdata_name_org, &
               label="metdata_name_org:", default = 'GMAO', rc=STATUS )
    VERIFY_(STATUS)

    call ESMF_ConfigGetAttribute(gmiConfigFile, self%metdata_name_model, &
              label="metdata_name_model:", default = 'GEOS-5', rc=STATUS )
    VERIFY_(STATUS)

    call ESMF_ConfigGetAttribute(gmiConfigFile, self%chem_opt, &
              label="chem_opt:", default = 2, rc=STATUS )
    VERIFY_(STATUS)

    call ESMF_ConfigGetAttribute(gmiConfigFile, self%trans_opt, &
              label="trans_opt:", default = 0, rc=STATUS )
    VERIFY_(STATUS)

    call ESMF_ConfigGetAttribute(gmiConfigFile, self%met_opt, &
              label="met_opt:", default = 3, rc=STATUS )
    VERIFY_(STATUS)


!-------------------------------------------
! Should BC files have current date and time?
! Useful for mission support and replays.
!--------------------------------------------
      
      call ESMF_ConfigGetAttribute(gmiConfigFile, value=self%BCRealTime, &
     &           label="BCRealTime:", default=.false., rc=STATUS)
      VERIFY_(STATUS)
     
      if (self%do_gcr) call INIT_GCR_DIAG(i1,i2,j1,j2,1,km)
      
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

      CALL InitializeEmission(self%Emission, self%SpeciesConcentration,       &
                     self%gmiGrid, gmiConfigFile, self%cellArea, IHNO3, IO3,  &
                     NSP, loc_proc, rootProc, self%chem_opt, self%trans_opt,  &
                     self%pr_diag,                                            &
                     self%pr_const, self%pr_surf_emiss, self%pr_emiss_3d, tdt)

      IF (BTEST(self%Emission%emiss_opt,1)) THEN

         allocate(self%Emission%soil_fert(i1:i2, ju1:j2))
         allocate(self%Emission%soil_precip(i1:i2, ju1:j2))

      END IF


    if (self%Emission%num_point_emiss .gt. 0) then 
      allocate(self%GmiPointEmiss(self%Emission%num_point_emiss))
    endif

      CALL initReadEmission(self%Emission, self%gmiClock, self%gmiGrid,     &
     &               self%cellArea, loc_proc, self%pr_diag, RC=STATUS)
      VERIFY_(STATUS)

!     !-----------------------------------
!     ! Initialize Surface Emission Bundle
!     !-----------------------------------
!     IF (self%Emission%do_semiss_inchem) THEN
!         call ESMF_StateGet(expChem, 'surfEmissForChem' , surfEmissBundle,   RC=STATUS)
!         VERIFY_(STATUS)
!
!         do ib = 1, NSP
!            allocate( var(i1:i2, j1:j2), STAT=STATUS)
!            VERIFY_(STATUS)
!            var(:,:)  = 0.0d0
!
!            write (binName ,'(i4.4)') ib
!            varName = 'surfEmiss'//binName
!
!            call addTracerToBundle (surfEmissBundle, var, w_c%grid_esmf, varName)
!         end do
!
!         call ESMF_FieldBundleGet(surfEmissBundle, fieldCount=numVars , rc=STATUS)
!         VERIFY_(STATUS)
!         _ASSERT(NSP == numVars,'needs informative message')
!     END IF

!---------------------------------------------------------------
! Create and populate the array that maps GMI species indices to
! GEOS-5 species indices
!---------------------------------------------------------------

    allocate(self%mapSpecies(NSP))
    self%mapSpecies(:) = speciesReg_for_CCM(lchemvar, NSP, bgg%reg%vname, bxx%reg%vname )

! Which (used) emission field is for lightning NO?
! ------------------------------------------------
    ic = 0
    self%ic_NO_lgt = -1
    CALL Get_num_emiss (self%Emission, num_emiss)
    DO i = 1,num_emiss
      speciesName = TRIM(self%Emission%emissionSpeciesNames(i))
      IF(self%Emission%emiss_map(i) >= 1) ic = ic+1
      IF(TRIM(speciesName) == 'NO_lgt') self%ic_NO_lgt = ic
    END DO

! Expose ambiguity in lightning emissions specification
! -----------------------------------------------------
    CALL Get_lightning_opt(self%Emission, lightning_opt)
    IF(lightning_opt == 0 .AND. self%ic_NO_lgt < 0) THEN
      IF( MAPL_AM_I_ROOT() ) THEN
        PRINT *,TRIM(Iam)//": Could not find emission specie name NO_lgt with lightning_opt = 0."
      END IF
      STATUS = 1
      VERIFY_(STATUS)
    END IF

! Grab the units for the export states of the EM_ emission species.  The
! first scan determines size of allocatable arrays.  Second scan fills them.
! --------------------------------------------------------------------------
    CALL ESMF_StateGet(expChem, ITEMCOUNT=n, RC=STATUS)
    VERIFY_(STATUS)
    ALLOCATE(itemNames(n), STAT=STATUS)
    VERIFY_(STATUS)
    ALLOCATE(itemTypes(n), STAT=STATUS)
    VERIFY_(STATUS)
    CALL ESMF_StateGet(expChem, ITEMNAMELIST=itemNames, ITEMTYPELIST=itemTypes, RC=STATUS)
    VERIFY_(STATUS)
   
    self%numEM_Exports = 0

    Scan: DO scanNumber = 1,2
      ic = 0

      Searching: DO m = 1,n

        TypeIsField: IF(itemTypes(m) == ESMF_StateItem_Field) THEN
          i = INDEX(TRIM(itemNames(m)), "EM_")

          Match: IF(i == 1) THEN
            CALL ESMF_StateGet(expChem, itemNames(m), FIELD, RC=STATUS)
            VERIFY_(STATUS)
            ic = ic+1
       
            IF(scanNumber == 2) THEN
              CALL ESMF_AttributeGet(FIELD, NAME='UNITS', VALUE=string, RC=status)
              VERIFY_(STATUS)
              self%EM_ExportNames(ic) = TRIM(itemNames(m))
              self%EM_ExportUnits(ic) = TRIM(string)
            END IF

          END IF Match

        END IF TypeIsField

      END DO Searching
    
      IF(scanNumber == 1) THEN
        self%numEM_Exports = ic
        ALLOCATE(self%EM_ExportNames(ic), STAT=STATUS)
        VERIFY_(STATUS)
        ALLOCATE(self%EM_ExportUnits(ic), STAT=STATUS)
        VERIFY_(STATUS)
      ELSE
        DEALLOCATE(itemNames, STAT=STATUS)
        VERIFY_(STATUS)
        DEALLOCATE(itemTypes, STAT=STATUS)
        VERIFY_(STATUS)
      END IF

    END DO Scan

! Use this space to hold the daily biomass burning emissions between
! Chem_UtilMPRead updates. The emissions are assumed to be singly layered.
! -----------------------------------------------------------------------

    self%veg_fraction_done = .FALSE.

!   Set up Overpass Masks
!   --------------------
    CALL OVP_init ( GC, "GMICHEM_DT:", LONS, OVP_RUN_DT, OVP_GC_DT, __RC__ ) !  Get LONS, timesteps

    IF(MAPL_AM_I_ROOT()) PRINT*,'GMI Emiss OVP_RUN_DT  OVP_GC_DT  :  ',OVP_RUN_DT, OVP_GC_DT

    ! In this case we update the Exports only after each GMI timestep:
    OVP_MASK_DT = OVP_GC_DT

    OVP_FIRST_HMS = OVP_end_of_timestep_hms( CLOCK, OVP_MASK_DT )
    IF(MAPL_AM_I_ROOT()) PRINT*,'GMI Emiss FIRST_HMS =',OVP_FIRST_HMS

    CALL OVP_mask ( LONS=LONS, DELTA_TIME=OVP_MASK_DT, OVERPASS_HOUR=10, MASK=MASK_10AM )
    CALL OVP_mask ( LONS=LONS, DELTA_TIME=OVP_MASK_DT, OVERPASS_HOUR=14, MASK=MASK_2PM  )

    RETURN
   
    END SUBROUTINE GmiEmiss_GridCompInitialize

!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GmiEmiss_initSurfEmissBundle
!
! !INTERFACE:
!
            SUBROUTINE GmiEmiss_initSurfEmissBundle (self, bgg, expChem, rc)
!
! !INPUT PARAMETERS:
   TYPE(Species_Bundle), INTENT(in) :: bgg      ! GMI Species - transported

! !OUTPUT PARAMETERS:
      INTEGER, INTENT(out) ::  rc               ! Error return code:
                                                !  0 - all is well
                                                !  1 - 
!
! !INPUT/OUTPUT PARAMETERS:
      TYPE(GmiEmiss_GridComp), INTENT(INOUT)  :: self      ! Grid Component
      TYPE(ESMF_State),        INTENT(INOUT)  :: expChem    ! Export State
!
! !DESCRIPTION: Initializes the surface emission bundle that is used
!      inside Chemistry. This is only done if the user wants surface
!      emission calculations to be done inside Chemistry.
!
! !REVISION HISTORY:
!    11May2017 Kouatchou   First crack.
!
! !LOCAL VARIABLES:
      INTEGER                    :: numVars, ib
      INTEGER                    :: STATUS
      character(len=4)           :: binName
      character(len=ESMF_MAXSTR) :: varName
      real(rPrec), pointer       :: var(:,:)
      type(ESMF_FieldBundle)     :: surfEmissBundle
      CHARACTER(LEN=*), PARAMETER :: IAm = 'GmiEmiss_initSurfEmissBundle'
!
!EOP
!-------------------------------------------------------------------------
!BOC
      RC = 0
      IF (self%Emission%do_semiss_inchem) THEN
         call ESMF_StateGet(expChem, 'surfEmissForChem' , surfEmissBundle, RC=STATUS)
         VERIFY_(STATUS)

         do ib = 1, NSP
            allocate( var(self%i1:self%i2, self%j1:self%j2), STAT=STATUS)
            VERIFY_(STATUS)
            var(:,:)  = 0.0d0

            write (binName ,'(i4.4)') ib
            varName = 'surfEmiss'//binName

            call addTracerToBundle (surfEmissBundle, var, bgg%grid_esmf, varName)
         end do

         call ESMF_FieldBundleGet(surfEmissBundle, fieldCount=numVars , RC=STATUS)
         VERIFY_(STATUS)
         ASSERT_(NSP == numVars)

         RC = STATUS
      END IF
      
      END SUBROUTINE GmiEmiss_initSurfEmissBundle
!EOC
!-------------------------------------------------------------------------

!-------------------------------------------------------------------------
!NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1, GEOS/DAS!
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  GmiEmiss_GridCompRun --- The GMI Driver
!
! !INTERFACE:
!

   SUBROUTINE GmiEmiss_GridCompRun ( self, bgg, bxx, impChem, expChem, nymd, nhms, &
                                     tdt, gc, clock, mixPBL, rc )

! !USES:

   USE GmiTimeControl_mod,            ONLY : Set_curGmiDate, Set_curGmiTime
   USE GmiTimeControl_mod,            ONLY : Set_numTimeSteps, Get_numTimeSteps
   USE GmiTimeControl_mod,            ONLY : Set_gmiSeconds, GetSecondsFromJanuary1
   USE GmiSpcConcentrationMethod_mod, ONLY : resetFixedConcentration
   USE GmiSolar_mod,                  ONLY : CalcCosSolarZenithAngle
   USE GmiEmissionMethod_mod,         ONLY : RunEmission
   USE GmiEmissionMethod_mod,         ONLY : Get_lightning_opt
!.sds
   USE Gmi_SEmissMethod_Mod
   use GmiSpeciesRegistry_mod, only : getSpeciesIndex
!... for point sources
   use m_StrTemplate,                 only : StrTemplate
   use GOCART2G_Process,              only : ReadPointEmissions
   use MAPL_BaseMod,                  only : MAPL_GetHorzIJIndex
!.sds.end
   IMPLICIT none

! !INPUT/OUTPUT PARAMETERS:

   TYPE(GmiEmiss_GridComp), INTENT(INOUT) :: self        ! Grid Component
   TYPE(Species_Bundle),    INTENT(INOUT) :: bgg         ! GMI Species - transported
   TYPE(Species_Bundle),    INTENT(INOUT) :: bxx         ! GMI Species - not transported
   TYPE(ESMF_State),        INTENT(INOUT) :: impChem     ! Import State
   TYPE(ESMF_State),        INTENT(INOUT) :: expChem     ! Export State
   TYPE(ESMF_Clock),        INTENT(INOUT) :: clock       ! The clock

!.sds
   type(ESMF_GridComp), intent(inout) :: gc     ! Gridded component 
   character (len=ESMF_MAXSTR)        :: COMP_NAME
   type(ESMF_Grid)                    :: grid
!.sds
! !INPUT PARAMETERS:

   INTEGER,                 INTENT(IN)    :: nymd, nhms  ! time
   REAL,                    INTENT(IN)    :: tdt         ! chemical timestep (secs)
   LOGICAL,                 INTENT(IN)    :: mixPBL      ! whether to explicitly distribute
                                                         ! aerosol emissions within the PBL
! !OUTPUT PARAMETERS:

   INTEGER,                 INTENT(OUT)   :: rc          ! Error return code:
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

   CHARACTER(LEN=*), PARAMETER :: IAm = 'GmiEmiss_GridCompRun'

!  Input fields from GEOS-5
!  ------------------------
   REAL, POINTER, DIMENSION(:,:) :: cn_prcp, tprec, lwi, zpbl, frlandice, snowdp
   REAL, POINTER, DIMENSION(:,:) :: T2m, u10m, v10m, ustar, z0h, swndsrf
   REAL, POINTER, DIMENSION(:,:) :: cldtt, wet1, dfpar, drpar

   REAL, POINTER, DIMENSION(:,:,:) :: airdens, ple, Q, T, zle
   REAL, POINTER, DIMENSION(:,:,:) :: cnv_mfc
!.sds
!... for DMS emissions
   real, pointer, dimension(:,:,:) :: delp
   real, pointer, dimension(:,:)   :: tskin
   real, pointer, dimension(:,:)   :: fr_ocean
   real, pointer, dimension(:,:)   :: DMS_ocean
!.sds.end

!  Exports not part of internal state
!  ----------------------------------
   REAL, POINTER, DIMENSION(:,:) :: emIsopSfc, emNOx, emMonot
   REAL, POINTER, DIMENSION(:,:) :: emBioCOMeth, emBioCOMonot, emBioPropene
   REAL, POINTER, DIMENSION(:,:) :: emSoilNOx, emShipHNO3, emShipO3

   REAL, POINTER, DIMENSION(:,:,:) :: EM_pointer

   REAL, POINTER, DIMENSION(:,:) :: ship_no
   REAL, POINTER, DIMENSION(:,:) :: jNO2val_phot

!  Local
!  -----
   INTEGER :: cymd, dymd, hms
   INTEGER :: i, i1, i2, ic, im
   INTEGER :: j, j1, j2, jm
   INTEGER :: k, km, kReverse
   INTEGER :: lightning_opt, loc_proc
   INTEGER :: m, n, STATUS
!  INTEGER :: num_emiss

   INTEGER, PARAMETER :: ToGMI = 1
   INTEGER, PARAMETER :: FromGMI = -1

   REAL :: mw, pi, degToRad, radToDeg, OneOverDt

   REAL, PARAMETER :: mwtAir = 28.9
   REAL, PARAMETER :: rStar = 8.314E+03
   REAL, PARAMETER :: Pa2hPa = 0.01
   REAL, PARAMETER :: ToGrPerKg = 1000.00
   REAL, PARAMETER :: secPerDay = 86400.00
   REAL, PARAMETER :: err = 1.00E-04

   REAL(KIND=DBL) :: chemDt, dayOfYear

   CHARACTER(LEN=ESMF_MAXSTR) :: speciesName
   CHARACTER(LEN=ESMF_MAXSTR) :: importName
   CHARACTER(LEN=ESMF_MAXSTR) :: fieldName
   CHARACTER(LEN=ESMF_MAXSTR) :: unitsName

   LOGICAL :: found, rootProc
!  LOGICAL, PARAMETER :: do_qqjk_reset = .TRUE.
   LOGICAL, PARAMETER :: doThis = .FALSE.
   
! Allocatables.  Use KIND=DBL where GMI expects REAL*8.
! -----------------------------------------------------
   INTEGER, ALLOCATABLE :: lwis_flags(:,:)

   REAL, ALLOCATABLE :: var2d(:,:)
   REAL, ALLOCATABLE :: var3d(:,:,:)
   REAL, ALLOCATABLE :: pl(:,:,:)

   REAL(KIND=DBL), ALLOCATABLE :: lonDeg(:,:)
   REAL(KIND=DBL), ALLOCATABLE :: latDeg(:,:)

   REAL(KIND=DBL), ALLOCATABLE :: TwoMeter_air_temp(:,:)
   REAL(KIND=DBL), ALLOCATABLE :: pctm2(:,:)
   REAL(KIND=DBL), ALLOCATABLE :: fracCloudCover(:,:)
   REAL(KIND=DBL), ALLOCATABLE :: surf_rough(:,:)
   REAL(KIND=DBL), ALLOCATABLE :: cosSolarZenithAngle(:,:)
   REAL(KIND=DBL), ALLOCATABLE :: radswg(:,:)
   REAL(KIND=DBL), ALLOCATABLE :: tenMeterU(:,:)
   REAL(KIND=DBL), ALLOCATABLE :: tenMeterV(:,:)
   REAL(KIND=DBL), ALLOCATABLE :: frictionVelocity(:,:)
   REAL(KIND=DBL), ALLOCATABLE :: con_precip(:,:)
   REAL(KIND=DBL), ALLOCATABLE :: tot_precip(:,:)
   REAL(KIND=DBL), ALLOCATABLE :: pbl(:,:)
   REAL(KIND=DBL), ALLOCATABLE :: soilWetness(:,:)
   REAL(KIND=DBL), ALLOCATABLE :: diffusePAR(:,:)
   REAL(KIND=DBL), ALLOCATABLE :: directPAR(:,:)
   REAL(KIND=DBL), ALLOCATABLE :: T_15_AVG(:,:)

   REAL(KIND=DBL), ALLOCATABLE :: var2dDBL(:,:)
   REAL(KIND=DBL), ALLOCATABLE :: var3dDBL(:,:,:)

   REAL(KIND=DBL), ALLOCATABLE :: mass(:,:,:)
   REAL(KIND=DBL), ALLOCATABLE :: press3c(:,:,:)
   REAL(KIND=DBL), ALLOCATABLE :: press3e(:,:,:)
   REAL(KIND=DBL), ALLOCATABLE :: height3e(:,:,:)
   REAL(KIND=DBL), ALLOCATABLE :: gridBoxThickness(:,:,:)
   REAL(KIND=DBL), ALLOCATABLE :: kel(:,:,:)
   REAL(KIND=DBL), ALLOCATABLE :: cmf(:,:,:)

!.sds
!... Workspace for point emissions
!   character(len=255)     :: point_emissions_srcfilen   ! filename for pointwise emissions
   integer :: nPts, l, irc
   real    :: ebot, etop, demissdz, z0, z1, dz, dPE    
   logical :: fileExists
!.sds.end
   REAL(KIND=DBL), ALLOCATABLE :: surfEmissForChem(:,:,:)
   INTEGER :: curRecord

   REAL*4, POINTER, DIMENSION(:,:,:) :: light_NO_prod => NULL()

   REAL, POINTER, DIMENSION(:,:,:) :: DATA_FOR_OVP_3D => NULL()
   REAL, POINTER, DIMENSION(:,:,:) :: OVP10_OUTPUT_3D => NULL()
   REAL, POINTER, DIMENSION(:,:,:) :: OVP14_OUTPUT_3D => NULL()
   INTEGER                         :: CURRENT_HMS  !  for the end of the timestep

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

!  Some real constants
!  -------------------
!   pi = 4.00*ATAN(1.00)
   pi = MAPL_PI
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
   if (self%Emission%do_ShipEmission) then
     ALLOCATE(jNO2val_phot(i1:i2,j1:j2),STAT=STATUS)
     VERIFY_(STATUS)
   endif

   CALL FindPointers(rc)

!  Reserve some local work space
!  -----------------------------
   ALLOCATE(lonDeg(i1:i2,j1:j2),STAT=STATUS)
   VERIFY_(STATUS)
   ALLOCATE(latDeg(i1:i2,j1:j2),STAT=STATUS)
   VERIFY_(STATUS)

   ALLOCATE(         lwis_flags(i1:i2,j1:j2),STAT=STATUS)
   VERIFY_(STATUS)
   ALLOCATE(              var2d(i1:i2,j1:j2),STAT=STATUS)
   VERIFY_(STATUS)
   ALLOCATE(              pctm2(i1:i2,j1:j2),STAT=STATUS)
   VERIFY_(STATUS)
   ALLOCATE(     fracCloudCover(i1:i2,j1:j2),STAT=STATUS)
   VERIFY_(STATUS)
   ALLOCATE(         surf_rough(i1:i2,j1:j2),STAT=STATUS)
   VERIFY_(STATUS)
   ALLOCATE(cosSolarZenithAngle(i1:i2,j1:j2),STAT=STATUS)
   VERIFY_(STATUS)
   ALLOCATE(             radswg(i1:i2,j1:j2),STAT=STATUS)
   VERIFY_(STATUS)
   ALLOCATE(          tenMeterU(i1:i2,j1:j2),STAT=STATUS)
   VERIFY_(STATUS)
   ALLOCATE(          tenMeterV(i1:i2,j1:j2),STAT=STATUS)
   VERIFY_(STATUS)
   ALLOCATE(   frictionVelocity(i1:i2,j1:j2),STAT=STATUS)
   VERIFY_(STATUS)
   ALLOCATE(         con_precip(i1:i2,j1:j2),STAT=STATUS)
   VERIFY_(STATUS)
   ALLOCATE(         tot_precip(i1:i2,j1:j2),STAT=STATUS)
   VERIFY_(STATUS)
   ALLOCATE(                pbl(i1:i2,j1:j2),STAT=STATUS)
   VERIFY_(STATUS)
   ALLOCATE(        soilWetness(i1:i2,j1:j2),STAT=STATUS)
   VERIFY_(STATUS)
   ALLOCATE(  TwoMeter_air_temp(i1:i2,j1:j2),STAT=STATUS)
   VERIFY_(STATUS)
   ALLOCATE(         diffusePAR(i1:i2,j1:j2),STAT=STATUS)
   VERIFY_(STATUS)
   ALLOCATE(          directPAR(i1:i2,j1:j2),STAT=STATUS)
   VERIFY_(STATUS)
   ALLOCATE(           T_15_AVG(i1:i2,j1:j2),STAT=STATUS)
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
   ALLOCATE(          height3e(i1:i2,j1:j2,0:km),STAT=STATUS)
   VERIFY_(STATUS)
   ALLOCATE(  gridBoxThickness(i1:i2,j1:j2,1:km),STAT=STATUS)
   VERIFY_(STATUS)
   ALLOCATE(               kel(i1:i2,j1:j2,1:km),STAT=STATUS)
   VERIFY_(STATUS)
   ALLOCATE(               cmf(i1:i2,j1:j2,1:km),STAT=STATUS)
   VERIFY_(STATUS)

   IF (self%Emission%do_semiss_inchem) THEN
     ALLOCATE(surfEmissForChem(i1:i2,j1:j2,1:NSP),STAT=STATUS)
     VERIFY_(STATUS)
   END IF

! Geolocation
! -----------
   lonDeg(i1:i2,j1:j2)=self%lonRad(i1:i2,j1:j2)*radToDeg
   latDeg(i1:i2,j1:j2)=self%latRad(i1:i2,j1:j2)*radToDeg

!  Layer mean pressures. NOTE: ple(:,:,0:km)
!  -----------------------------------------
   DO k=1,km
     pl(i1:i2,j1:j2,k) = (ple(i1:i2,j1:j2,k-1)+ple(i1:i2,j1:j2,k))*0.50
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
!  Vegetation Fraction
!  MEGAN emissions
! --------------------------------------------------------
   CALL Acquire_Clims(STATUS)
   VERIFY_(STATUS)

! Grab imports and do units conversions
! -------------------------------------
   CALL SatisfyImports(STATUS)
   VERIFY_(STATUS)

! Daily or monthly emissions inventories
! -----------------------------------------
   CALL Refresh_Daily(STATUS)
   VERIFY_(STATUS)

! Keep running average of T2M for the previous 15 days
! ----------------------------------------------------
   IF(.NOT. self%doingPredictorNow) THEN
     CALL MonitorT2M(STATUS)
     VERIFY_(STATUS)
   END IF
   T_15_AVG(i1:i2,j1:j2) = bxx%qa(bxx%reg%nq)%data3d(i1:i2,j1:j2,1)

! Hand the species concentrations to GMI's bundle
! -----------------------------------------------
   IF (self%gotImportRst) THEN
      CALL SwapSpeciesBundles(ToGMI, self%SpeciesConcentration%concentration,          &
               bgg%qa, bxx%qa, Q, self%mapSpecies, lchemvar, self%do_synoz, NSP, STATUS)
      VERIFY_(STATUS)
   END IF

   DEALLOCATE(var3D, STAT=STATUS)
   VERIFY_(STATUS)

! Impose fixed concentrations
! ---------------------------
   IF(self%SpeciesConcentration%num_fixed_const > 0) THEN
     CALL resetFixedConcentration(self%SpeciesConcentration, self%gmiClock, self%gmiGrid, NSP)
   END IF

! ------------------------------------------------------------------------
! Emission
!
! NOTE: In GEOS-5, the "instantaneous" emissions are held in surf_emiss_out2 and
!       and emiss_out_3d.  The accumulated emission is obtained by instantiating
!       time-averaged exports via the HISTORY.rc in the run script.
! ------------------------------------------------------------------------
   IF(self%do_emission .AND. self%gotImportRst) THEN

     IF(self%pr_emiss_3d)   self%Emission%emiss_3d_out    = 0.00D+00

     IF(self%pr_surf_emiss) self%Emission%surf_emiss_out2 = 0.00D+00

! Grab lightning option, and pass flash rates to Emission
!  0: Prescribed emission, NO_lgt
!  1: Parameterized emission, GmiEmissionLightning_mod.F90
!  2: None
! --------------------------------------------------------
     CALL Get_lightning_opt(self%Emission, lightning_opt)

     IF (self%Emission%do_ShipEmission) THEN

       CALL MAPL_GetPointer(impChem, ship_no, 'SHIP_NO', __RC__)

       call calcShipEmission (self%Emission%emiss_o3, &
                              self%Emission%emiss_hno3, &
                latDeg, jNO2val_phot, &
                self%Emission%emissionArray, ship_no,     &
                self%Emission%ship_o3_index, &
                self%Emission%ship_hno3_index, &
                self%cellArea, &
                self%gmiGrid%i1, self%gmiGrid%i2, self%gmiGrid%ju1, &
                self%gmiGrid%j2, self%Emission%num_emiss)
       DEALLOCATE(jNO2val_phot)

     ENDIF

     CALL MAPL_GetPointer(impChem, light_NO_prod, 'LIGHT_NO_PROD', __RC__)


!--------------------------------------------------------
! Calculate the air density at the center of each grid box
! (molecules/cm^3).
!--------------------------------------------------------
    self%SpeciesConcentration%concentration(IMGAS)%pArray3D(:,:,:) =  &
                     press3c(i1:i2,j1:j2,:) * MB2CGS /  &
                     (kel (i1:i2,j1:j2,:) * BOLTZMN_E)

    CALL RunEmission(self%Emission, self%SpeciesConcentration, self%gmiClock,  &
                     self%gmiGrid,                                             &
                     loc_proc, cosSolarZenithAngle, latdeg, self%cellArea,     &
                     mass, lwis_flags, radswg, TwoMeter_air_temp, surf_rough,  &
                     con_precip, tot_precip, frictionVelocity, fracCloudCover, &
                     kel, pbl, cmf, press3c, press3e, pctm2,                   &
                     tenMeterU, tenMeterV, soilWetness, gridBoxThickness,      &
                     mw_data, IBOC, IBBC, INOC, IFOC, IFBC, ISSLT1, ISSLT2,    &
                     ISSLT3, ISSLT4, IFSO2, INSO2, INDMS, IAN, IMGAS, INO,     &
                     IC5H8, INO, ICO, IC3H6, IHNO3, IO3, NSP, diffusePAR,      &
                     directPAR, T_15_AVG, self%met_opt, self%chem_opt,         &
                     self%trans_opt, self%do_aerocom, self%do_drydep,          &
                     self%pr_diag, self%pr_const, self%pr_surf_emiss,          &
                     self%pr_emiss_3d, self%metdata_name_org,                  &
                     self%metdata_name_model, tdt, mixPBL, light_NO_prod)

!... inchem emissions aren't working
     if (self%Emission%do_semiss_inchem) then
       call updateSurfEmissionInChemistry (self%pr_surf_emiss, self%pr_emiss_3d, &
                   self%Emission%emiss_isop, self%Emission%emiss_monot,           &
                   self%Emission%emiss_nox, self%Emission%do_ShipEmission,        &
                   self%Emission%emiss_hno3, self%Emission%emiss_o3, IHNO3, IO3,  &
                   self%cellArea, self%Emission%surf_emiss_out,                   &
                   self%Emission%surf_emiss_out2, self%Emission%emiss_3d_out,     &
                   self%Emission%emissionArray, surfEmissForChem,                 &
                   gridBoxThickness, self%Emission%emiss_timpyr,                  &
                   self%Emission%num_emiss, self%Emission%emiss_opt,              &
                   self%Emission%emiss_map, chemDt, nymd, ICO, INO, IC3H6, IC5H8, &
                   mw_data, self%pr_diag, loc_proc, self%gmiGrid%i1,              &
                   self%gmiGrid%i2, self%gmiGrid%ju1, self%gmiGrid%j2,            &
                   self%gmiGrid%k1, self%gmiGrid%k2, self%gmiGrid%ilo,            &
                   self%gmiGrid%ihi, self%gmiGrid%julo, self%gmiGrid%jhi, NSP)

        call populateBundle()
     endif
   ENDIF
 
! Return species concentrations to the chemistry bundle
! -----------------------------------------------------
   IF (self%gotImportRst) THEN
      CALL SwapSpeciesBundles(FromGMI, self%SpeciesConcentration%concentration,   &
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

   DEALLOCATE(lwis_flags, var2d, TwoMeter_air_temp, &
              pctm2, fracCloudCover, surf_rough, cosSolarZenithAngle, &
	      radswg, tenMeterU, tenMeterV, frictionVelocity, con_precip, &
	      tot_precip, pbl, soilWetness, diffusePAR, directPAR, T_15_AVG, STAT=STATUS)
   VERIFY_(STATUS)

   DEALLOCATE(pl, mass, press3c, press3e, height3e, gridBoxThickness, kel, cmf, &
              STAT=STATUS)
   VERIFY_(STATUS)

   IF (self%Emission%do_semiss_inchem) THEN
     DEALLOCATE(surfEmissForChem, STAT=STATUS)
     VERIFY_(STATUS)
   END IF

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
!  Obtain the following climatologies from files specified in ExtData.rc:
!   Fixed concentration species
!   Stratospheric sulfate surface area
!   Emissions

!  For details on decoding land use and leaf area indices, see Chem_Shared/VegLaiMod.F90

!EOP
!---------------------------------------------------------------------------

  CHARACTER(LEN=255) :: IAm
  CHARACTER(LEN=  3) :: numID

  INTEGER :: ic, p
  INTEGER :: lightning_opt
  INTEGER, ALLOCATABLE :: landNum(:)
  INTEGER, ALLOCATABLE :: milFrac(:)

  REAL :: tokgCPerBox

  REAL, POINTER, DIMENSION(:,:) :: PTR2D
  REAL, POINTER, DIMENSION(:,:,:) :: PTR3D

  REAL, PARAMETER :: VEG_MIN_VALID = 1.0

  rc = 0
  IAm = "Acquire_Clims"

! Grab lightning option
! ---------------------
  CALL Get_lightning_opt(self%Emission,lightning_opt)

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

! Update the vegetation boundary condition
! ----------------------------------------
  IF(BTEST(self%Emission%emiss_opt,1)) THEN

    IF ( .NOT. self%veg_fraction_done ) THEN

      importName = 'VEG_FRAC'
      CALL MAPL_GetPointer(impChem, PTR3D, TRIM(importName), RC=STATUS)
      VERIFY_(STATUS)

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

    importName = 'LAI_FRAC'
    CALL MAPL_GetPointer(impChem, PTR3D, importName, RC=STATUS)
    VERIFY_(STATUS)

    CALL Decode_XLAI(PTR3D, NTYPE,              &
                    i1, i2, j1, j2, NVEGTYPE-2,&
                    self%Emission%ireg,        &
                    self%Emission%iuse,        &
                    self%Emission%iland,       &
                    self%Emission%xlai,        &
                    RC=STATUS)
    VERIFY_(STATUS)

    NULLIFY(PTR3D)

     !------------------------------
     ! Read the Fertilizer File
     !------------------------------

    importName = 'SOILFERT'
    CALL MAPL_GetPointer(impChem, PTR2D, TRIM(importName), RC=STATUS)
    VERIFY_(STATUS)
    self%Emission%soil_fert(:,:) = PTR2D(:,:)
    NULLIFY(PTR2D)

     !------------------------------
     ! Read the Precipitation File
     !------------------------------

    importName = 'SOILPRECIP'
    CALL MAPL_GetPointer(impChem, PTR2D, TRIM(importName), RC=STATUS)
    VERIFY_(STATUS)
    self%Emission%soil_precip(:,:) = PTR2D(:,:)
    NULLIFY(PTR2D)

  END IF

! MEGAN annual emission factors and leaf-area-index
! -------------------------------------------------
  MEGAN: IF(self%Emission%doMEGANEmission) THEN

! Read emission factors and convert from [micrograms C m^{-2} hr^{-1}] to [kg C box^{-1} dt^{-1}]
! -----------------------------------------------------------------------------------------------
    tokgCPerBox = tdt*1.00E-09/3600.00

    IF(self%Emission%doMEGANviaHEMCO) THEN 
   ! get MEGAN emissions pointers from HEMCO  (already kgC/m2/s)
       CALL MAPL_GetPointer(impChem,   PTR2D, 'GMI_ISOPRENE', RC=STATUS)
       VERIFY_(STATUS)
!     CALL MAPL_MaxMin('emiss_isop ptr in gmiEmiss:', PTR2D)
       self%Emission%emiss_isop(:,:) = PTR2D(:,:)
       NULLIFY(PTR2D)
    END IF

    importName = 'MEGAN_ISOP'
    CALL MAPL_GetPointer(impChem, PTR2D, TRIM(importName), RC=STATUS)
    VERIFY_(STATUS)
    self%Emission%aefIsop(:,:) = PTR2D(:,:)*self%cellArea(:,:)*tokgCPerBox
    NULLIFY(PTR2D)

    importName = 'MEGAN_MBO'
    CALL MAPL_GetPointer(impChem, PTR2D, TRIM(importName), RC=STATUS)
    VERIFY_(STATUS)
    self%Emission%aefMbo(:,:) = PTR2D(:,:)*self%cellArea(:,:)*tokgCPerBox
    NULLIFY(PTR2D)

    importName = 'MEGAN_MPE'
    CALL MAPL_GetPointer(impChem, PTR2D, TRIM(importName), RC=STATUS)
    VERIFY_(STATUS)
    self%Emission%aefMonot(:,:) = PTR2D(:,:)*self%cellArea(:,:)*tokgCPerBox
    NULLIFY(PTR2D)

    importName = 'MEGAN_OVOC'
    CALL MAPL_GetPointer(impChem, PTR2D, TRIM(importName), RC=STATUS)
    VERIFY_(STATUS)
    self%Emission%aefOvoc(:,:) = PTR2D(:,:)*self%cellArea(:,:)*tokgCPerBox
    NULLIFY(PTR2D)

! Leaf area index: The MEGAN code requires the previous, current, and next months'
! LAI values without temporal interpolation. Given the correct refresh template in
! ExtData.rc ("-"), the following construct assures that ExtData initializes all twelve
! monthly values of isoLaiYear only on the first pass through the GmiEmission run method.
! ---------------------------------------------------------------------------------------
    DO i = 1,12
      WRITE(numID,'(I3.3)') i
      importName = 'MEGAN_LAI_'//numID
      CALL MAPL_GetPointer(impChem, PTR2D, TRIM(importName), RC=STATUS)
      VERIFY_(STATUS)
      self%Emission%isoLaiYear(:,:,i) = PTR2D(:,:)
      NULLIFY(PTR2D)
    END DO

    m = MOD(nymd,10000)/100

    IF(m == 1) THEN
      p = 12
    ELSE
      p = m-1
    END IF

    IF(m == 12) THEN
      n = 1
    ELSE
      n = m+1
    END IF

    self%Emission%isoLaiPrev(:,:) = self%Emission%isoLaiYear(:,:,p)
    self%Emission%isoLaiCurr(:,:) = self%Emission%isoLaiYear(:,:,m)
    self%Emission%isoLaiNext(:,:) = self%Emission%isoLaiYear(:,:,n)

  END IF MEGAN

  RETURN
  END SUBROUTINE Acquire_Clims

!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  Refresh_Daily
!
! !INTERFACE:

  SUBROUTINE Refresh_Daily(rc)
  
  IMPLICIT NONE

  INTEGER, INTENT(OUT) :: rc
!
! !DESCRIPTION:
!
!  Update emissions with daily inventories from A. Darmenov supplemented with monthly 
!  EDGAR or EDGAR/TRANSCOM CO, NO, and CH4 from L. Ott.
!
!  File names and update intervals are now specified by ExtData.rc
!
! !RESTRICTIONS:
!  Assumes that any non_zero entries of emiss_map() occur at the front of the vector.
!  Assumes diurnal emissions are singly-layered.
!  Propene is PRPE in GMI but is c3h6 on QFED2 emission files.
!  The field with short name "biomass" is the only one injested.
!  Each inventory has its own file.  On the ExtData.rc file it is
!   designated, for example:
!
!     ExtData/Y%y4/M%m2/qfed2.emis_xxx.005.%y4%m2%d2.nc4
!
!   For each daily inventory, the xxx is replaced with the lower case specie name.
!
!  Also update the fossil fuel inventories from monthly averages here instead
!  of in Acquire_Clims.  This is one-layer data, too.
!
!  WARNING: It is NOT required that daily and monthly emissions are sourced
!           from the same year. 
!
!EOP
!---------------------------------------------------------------------------

  CHARACTER(LEN=ESMF_MAXSTR) :: IAm
  INTEGER :: STATUS
  INTEGER :: i, k, kReverse
  INTEGER :: lightning_opt
!  INTEGER :: num_emiss
  REAL, ALLOCATABLE :: cellWeighting(:,:)
  REAL, ALLOCATABLE :: weightedField2D(:,:)

  REAL, POINTER, DIMENSION(:,:) :: PTR2D
  REAL, POINTER, DIMENSION(:,:,:) :: PTR3D
!.sds.. DMS emissions
  integer :: indDMS
  real, allocatable :: tbot(:,:)
  real, allocatable :: flux_DMS(:,:)
  real, allocatable :: DMS_atm(:,:,:)
  logical :: do_achem_dms_emiss
  character (len=ESMF_MAXSTR)  :: fname ! file name for point source emissions
!.sds.. point emissions
  integer :: num_point_emiss
!... local point source vrbls
  integer           :: nCols, itmp, sslen, selen
  type(ESMF_Config) :: cf
  CHARACTER(LEN=ESMF_MAXSTR) :: tmpstr
  real :: scale
!.sds.end point emissions

  rc = 0
  IAm = "Refresh_Daily"

! Note: emissionArray is allocated and set to zero in intializeEmission.

! ------------------------------- I M P O R T A N T --------------------------------
! The emissionArray units are kg s^{-1}, which is a legacy attribute. However, for
! proper mapping from regularly-spaced latitudes and longitudes to the cubed sphere
! the fluxes on the emissions file MUST be per unit area. When clim_emiss_by_area is 
! .TRUE., the units on the emission files are assumed kg m^{-2} s^{-1}, and 
! multiplication by cellArea follows data acquisition by MAPL_ExtData.
! ----------------------------------------------------------------------------------
  ALLOCATE(cellWeighting(i1:i2,j1:j2),weightedField2D(i1:i2,j1:j2),__STAT__)

  IF(self%Emission%clim_emiss_by_area) THEN
    cellWeighting(i1:i2,j1:j2) = self%cellArea(i1:i2,j1:j2)
    IF(MAPL_AM_I_ROOT()) PRINT *, TRIM(IAm)//": Emissions are per unit area."
  ELSE
    cellWeighting(i1:i2,j1:j2) = 1.00
    IF(MAPL_AM_I_ROOT()) PRINT *, TRIM(IAm)//": Emissions are NOT per unit area."
  END IF

  CALL Get_lightning_opt(self%Emission,lightning_opt)


! For each available inventory ...
! --------------------------------
!  num_emiss = count( self%Emission%emiss_map(:) > 0 )
  DO i = 1,self%Emission%num_point_start-1

    speciesName = TRIM(self%Emission%emissionSpeciesNames(i))

!... Special cases: 
! Parameterized ship emissions are handled elsewhere (calcShipEmission)
    IF ( TRIM(speciesName) == '*shipO3*'.OR. TRIM(speciesName) == '*shipHNO3*' ) CYCLE
!.sds.. DMS emission is special case, emissions handled elsewhere
    if(TRIM(speciesName) == 'DMS' ) cycle

    IF ( self%Emission%emiss_map(i) < 1 ) THEN
      PRINT *, ' '
      PRINT *,'GmiEmiss::'//TRIM(IAm)//' Species '//TRIM(speciesName)//' does not map to GMI species.'
      STATUS = 1
      VERIFY_(STATUS)
    END IF

! Select single- or multi-layered based on emissionSpeciesLayers
! --------------------------------------------------------------
    IF(self%Emission%emissionSpeciesLayers(i) == 1) THEN

      CALL MAPL_GetPointer(impChem, PTR2D, TRIM(speciesName), __RC__)

      IF(i <= self%num_diurnal_emiss) THEN
        weightedField2D(:,:) = PTR2D(:,:)*cellWeighting(:,:)
        CALL Chem_BiomassDiurnal(var2D,                      &
                               weightedField2D(:,:),         &
                               self%lonRad(:,:)*radToDeg,    &
                               self%latRad(:,:)*radToDeg,    &
                               nhms, tdt)
        self%Emission%emissionArray(i)%pArray3D(:,:,1) = var2D(:,:)
      ELSE
        self%Emission%emissionArray(i)%pArray3D(:,:,1) = PTR2D(:,:)*cellWeighting(:,:)
      END IF

      self%Emission%emissionArray(i)%pArray3D(:,:,2:km) = 0.00D+00

      NULLIFY(PTR2D)

    ELSE
  
      CALL MAPL_GetPointer(impChem, PTR3D, TRIM(speciesName), __RC__)

      IF(i <= self%num_diurnal_emiss) THEN
        PRINT *, ' '
        PRINT *,'GmiEmiss::'//TRIM(IAm)//':  Species '//TRIM(speciesName)//' cannot be diurnal and 3D.'
        STATUS = 1
        VERIFY_(STATUS)
      ELSE
        DO k = 1,km
          kReverse = km-k+1
          self%Emission%emissionArray(i)%pArray3D(:,:,kReverse) = PTR3D(:,:,k)*cellWeighting(:,:)
        END DO
      END IF

      NULLIFY(PTR3D)

! Use prescribed NO from lightning only if lightning_opt is zero
! lightning_opt = 0  -> read from file
! lightning_opt = 1  -> parameterized
! lightning_opt = 2  -> none
! --------------------------------------------------------------
      IF (TRIM(speciesName) == 'NO_lgt' .AND. lightning_opt >= 1) THEN
        self%Emission%emissionArray(i)%pArray3D(:,:,:) = 0.00
      END IF

    END IF

  END DO
!.sds
!... do point emissions
  if(self%Emission%num_point_emiss .gt. 0) then
!
    DO k=1,self%Emission%num_point_emiss
!... put point emissions into correct GMI Emission%emissionArray(ic)%pArray3D() slot
      ic = self%Emission%num_point_start + k-1
!... zero out initial pass
      self%Emission%emissionArray(ic)%pArray3D(:,:,:) = 0.0d0
!
      tmpstr = TRIM(self%Emission%emissionPointFilenames(k))
      call StrTemplate(fname, TRIM(tmpstr), xid='unknown', nymd=nymd, nhms=120000 )
      inquire(file=fname, exist=fileExists)
!
!... is it a new filename this time step?
!... read point emission data
      if ( fileExists .and. trim(ptfile_save(k)) .ne. trim(fname) ) then
        ptfile_save(k) = fname
!... get no. of emission records and allocate
        tmpstr = TRIM(self%Emission%num_point_type(k))
        cf = ESMF_ConfigCreate()
        call ESMF_ConfigLoadFile(cf, fileName=trim(fname), rc=STATUS )
        call ESMF_ConfigGetDim(cf, nPts, nCols, LABEL=TRIM(tmpstr)//'::', rc=STATUS)
        self%GmiPointEmiss(k)%nPts = nPts

        if(allocated(self%GmiPointEmiss(k)%iPoint)) then
          deallocate(self%GmiPointEmiss(k)%iPoint,__STAT__)
          deallocate(self%GmiPointEmiss(k)%jPoint,__STAT__)
          deallocate(self%GmiPointEmiss(k)%vLat,  __STAT__)
          deallocate(self%GmiPointEmiss(k)%vLon,  __STAT__)
          deallocate(self%GmiPointEmiss(k)%vBase, __STAT__)
          deallocate(self%GmiPointEmiss(k)%vTop,  __STAT__)
          deallocate(self%GmiPointEmiss(k)%vEmis, __STAT__)
          deallocate(self%GmiPointEmiss(k)%vStart,__STAT__)
          deallocate(self%GmiPointEmiss(k)%vEnd,  __STAT__)
        endif
!
        allocate(self%GmiPointEmiss(k)%iPoint(nPts),__STAT__)
        allocate(self%GmiPointEmiss(k)%jPoint(nPts),__STAT__)
        allocate(self%GmiPointEmiss(k)%vLat(nPts),  __STAT__)
        allocate(self%GmiPointEmiss(k)%vLon(nPts),  __STAT__)
        allocate(self%GmiPointEmiss(k)%vBase(nPts), __STAT__)
        allocate(self%GmiPointEmiss(k)%vTop(nPts),  __STAT__)
        allocate(self%GmiPointEmiss(k)%vEmis(nPts), __STAT__)
        allocate(self%GmiPointEmiss(k)%vStart(nPts),__STAT__)
        allocate(self%GmiPointEmiss(k)%vEnd(nPts),  __STAT__)
        self%GmiPointEmiss(k)%vStart(:) = 0
        self%GmiPointEmiss(k)%vEnd(:)   = 240000
!... read point emissions
        call ReadPointEmissions (nymd, fname, nPts, self%GmiPointEmiss(k)%vLat &
          , self%GmiPointEmiss(k)%vLon, self%GmiPointEmiss(k)%vBase, self%GmiPointEmiss(k)%vTop &
          , self%GmiPointEmiss(k)%vEmis, self%GmiPointEmiss(k)%vStart, self%GmiPointEmiss(k)%vEnd &
          , label=tmpstr, __RC__)
!... specified variable not in .rc point emission file
        if(nPts .le. 0) then
          if(mapl_am_i_root()) print *, 'Point Emission Variable ', trim(tmpstr) &
            ,' not found in ', trim(fname)
          nPts = -1
!... data read in .rc point emission file
        else
          self%GmiPointEmiss(k)%vLon(:) = self%GmiPointEmiss(k)%vLon(:)*real(MAPL_DEGREES_TO_RADIANS)
          self%GmiPointEmiss(k)%vLat(:) = self%GmiPointEmiss(k)%vLat(:)*real(MAPL_DEGREES_TO_RADIANS)
          if(mapl_am_i_root()) then
            selen = len_trim(fname)
            sslen = max(1,selen-70)
            print '(i8,'' emission points in: ...'',a72)', nPts, trim(fname(sslen:selen))
            print '('' '')'
          endif
        endif
!... no file...
      else if ( .not.fileExists .and. trim(ptfile_save(k)) .ne. trim(fname) ) then
        if(mapl_am_i_root()) print *,'GMI point emiss fn: ', trim(fname), ' not found; proceeding...'
        nPts = -1 ! set this back to -1 so the "if (self%nPts > 0)" conditional is not exercised.
!... retrieve number of point source points...
      else
        nPts = self%GmiPointEmiss(k)%nPts
      endif
!
!... Get indices for point emission locations
      if(nPts .gt. 0) then
        self%GmiPointEmiss(k)%iPoint(:) = -1
        self%GmiPointEmiss(k)%jPoint(:) = -1
        call ESMF_GridCompGet (gc, grid=grid, NAME=COMP_NAME, __RC__)
        call MAPL_GetHorzIJIndex(nPts, self%GmiPointEmiss(k)%iPoint &
          , self%GmiPointEmiss(k)%jPoint, grid=grid, lon=self%GmiPointEmiss(k)%vLon &
          , lat=self%GmiPointEmiss(k)%vLat, rc=status)
        if(status /= 0) then
          if(mapl_am_i_root()) print*, trim(Iam), ' - cannot get indices for point emissions'
          VERIFY_(status)
        endif
!... loop over all points in file
        do n = 1,nPts 
          i = self%GmiPointEmiss(k)%iPoint(n)
          j = self%GmiPointEmiss(k)%jPoint(n)
!... emission in sub-domain?
          if (i.lt.1 .or. j.lt.1) cycle
!... during emission period?
!    print '(''sds01: '',6i8,2f12.2)',k , n, i, j, vStart(n), vEnd(n), vBase(n), vTop(n)
          if(nhms .ge. self%GmiPointEmiss(k)%vStart(n) &
            .and. nhms .lt. self%GmiPointEmiss(k)%vEnd(n)) then
!... distribute in the vertical (m)
            ebot = self%GmiPointEmiss(k)%vBase(n)
            etop = self%GmiPointEmiss(k)%vTop(n)
!
!.sds.. can this be handled better from file input?
!
!... scale needed for emissions, volcano files are in kg(S) emitted into SO2 so scale .ne. 1.0
            scale = 1.0
!... volcanic point emissions are only in the top 1/3 of emission column, move ebot
            if (TRIM(tmpstr) .eq. 'volcano' .and. etop .ne. ebot) then
              ebot = etop - (etop-ebot)/3.
!... volcanic emissions are in kg(S), need kg(SO2)
              CALL getMW(TRIM('EM_SO2'), itmp, mw, rc)
              scale = mw/32.06
            endif
!
!... distribute into column (z1 is gridbox top, z0 is gridbox bottom, arrays are bottom-up)
            do l = 1, km
              z0 = height3e(i,j,l-1)
              z1 = height3e(i,j,l)
              dz = z1 - z0
              dPE = 0.0
!... emission (kg/s)/m if (ebot .ne. etop) 
              if (etop .ne. ebot) demissdz = self%GmiPointEmiss(k)%vEmis(n)/(etop-ebot)
!... emission is above this level "cycle"
              if (ebot .gt. z1) cycle
!... emission is below this level (except if below model surface) "cycle"
              if (z0 .gt. etop .and. l .ne. km) cycle
!... emiss all below model bottom, putting all in bottom level
              if ( l .eq. 1 .and. z1 .ge. etop ) then
                dPE = self%GmiPointEmiss(k)%vEmis(n)
!... emiss all in this single level
              elseif ( etop .le. z1 .and. ebot .ge. z0 ) then   ! completely in single layer
                dPE = self%GmiPointEmiss(k)%vEmis(n)
!... bottom of emiss cloud only in this level
              elseif ( (ebot .lt. z1 .and. ebot .ge. z0) .and. etop .gt. z1) then
                dPE = (z1-ebot) * demissdz
!... layer is filled with cloud which goes above and below
              elseif ( etop .gt. z1 .and. ebot .lt. z0) then
                dPE = (z1-z0) * demissdz
!... top of emiss cloud only in this level
              elseif ( (etop .le. z1 .and. etop .ge. z0) .and. ebot .lt. z0) then
                dPE = (etop-z0) * demissdz
              endif
!... add to GMI emission array (kg/gridbox/s)*timestep 
              self%Emission%emissionArray(ic)%pArray3D(i,j,l) = &
                self%Emission%emissionArray(ic)%pArray3D(i,j,l) + dPE * scale
            enddo
          endif
        enddo
!
      endif
    enddo
  endif
!
!.sds
!... do DMS emission if DMS is a carried species
!... get   constituent number
!  indDMS = getSpeciesIndex('DMS',.true.)
  if(self%Emission%GmiDMSEmissIndex .gt. 0) then
!
    allocate(tbot(i1:i2,j1:j2),  __STAT__)
    allocate(flux_DMS(i1:i2,j1:j2),  __STAT__)
    allocate(DMS_atm(i1:i2,j1:j2,1:km),  __STAT__)
!
!... moved to Refresh_Daily
    call MAPL_GetPointer(impChem, tskin,     'TS',        __RC__)
    call MAPL_GetPointer(impChem, fr_ocean,  'FROCEAN',   __RC__)
    call MAPL_GetPointer(impChem, DMS_ocean, 'DMS_OCEAN', __RC__)
!... sds end
!
!... do ACHEM like DMS flux calculation?
    do_achem_dms_emiss = .false.
!
    if(do_achem_dms_emiss) then
      tbot(:,:) = tskin(:,:)
    else
      tbot(:,:) = kel(:,:,1)
    endif
!... copy in DMS to temp array and flip in vertical
    DMS_atm(:,:,:) = self%SpeciesConcentration%concentration(self%Emission%GmiDMSEmissIndex)%pArray3D(:,:,km:1:-1)
!
!... calculate DMS flux emissions from the ocean
    call GmiDMS_emissions(do_achem_dms_emiss &
                     ,tbot        &
                     ,u10m        &
                     ,v10m        &
                     ,fr_ocean    &
                     ,DMS_ocean   &
                     ,DMS_atm     &
                     ,flux_DMS    &
                     ,irc)
!
!... flux_DMS in is kg m-2 s-1, emissionArray is kg box-1 s-1
!... Emissions diagnostic and Emissions Array are in units kg box-1 s-1, convert and assign
    self%Emission%emissionArray(self%Emission%GmiDMSEmissIndex)%pArray3D(i1:i2,j1:j2,1) &
      = flux_DMS(:,:) * self%cellArea(:,:)
!
    deallocate(flux_DMS)
    deallocate(DMS_atm)
    deallocate(tbot)
  endif
!
!.sds.end
  DEALLOCATE(cellWeighting, weightedField2D, __STAT__)

! LLNL emissions units conversion
! -------------------------------
  IF(self%Emission%emiss_conv_flag /= 0) THEN
    IF( MAPL_AM_I_ROOT() ) THEN
      PRINT *,TRIM(IAm),': Code not ready for emiss_conv_flag =',self%Emission%emiss_conv_flag
    END IF
    STATUS = 1
    VERIFY_(STATUS)
  END IF

  RETURN
  END SUBROUTINE Refresh_Daily

!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  FillExports
!
! !INTERFACE:

  SUBROUTINE FillExports(rc)

  USE gcr_mod,                       ONLY : GET_GCR_EMISS
    
  IMPLICIT NONE

   INTEGER, INTENT(OUT) :: rc
!
! !DESCRIPTION:
!
!  Fill the export state
!
!EOP
!---------------------------------------------------------------------------
   type(ESMF_FieldBundle)      :: sadBundle
   real(rPrec), pointer        :: ptr3D(:,:,:)
   real*4                      :: tmpisop(i1:i2,j1:j2) ! for pmaxmin

   CHARACTER(LEN=255) :: IAm
   rc = 0
   IAm = "FillExports"

! ------------------------------------------------------------------------
! Emissions
! ------------------------------------------------------------------------

! Surface emissions.  GmiEmission says that
! these are kg s^{-1}. Convert to kg m^{-2} s^{-1}.
! -------------------------------------------------
   IF(self%do_emission .AND. self%Emission%emiss_opt == 2) THEN
     IF(ASSOCIATED(emNOx))         emNOx(:,:) = self%Emission%emiss_nox(  :,:)/self%cellArea(:,:)
     IF(ASSOCIATED(emMonot))     emMonot(:,:) = self%Emission%emiss_monot(:,:)/self%cellArea(:,:)
     IF(ASSOCIATED(emIsopSfc)) emIsopSfc(:,:) = self%Emission%emiss_isop( :,:)/self%cellArea(:,:)

     IF(ASSOCIATED(emIsopSfc)) THEN
       tmpisop = emIsopSfc(i1:i2,j1:j2) 
       CALL MAPL_MaxMin('emIsopSfc in gmiEmiss:',tmpisop(:,:))
     END IF
   END IF

! Biogenic CO, soil NOx, and surface ship emissions.  GmiEmission says 
! that these are kg m^{-2} per time step.  Convert to kg m^{-2} s^{-1}.
! ---------------------------------------------------------------------
   IF(self%do_emission .AND. self%pr_surf_emiss) THEN
     IF(ASSOCIATED( emBioCOMeth))   emBioCOMeth(:,:) = self%Emission%surf_emiss_out2(:,:,1)*OneOverDt
     IF(ASSOCIATED(emBioCOMonot))  emBioCOMonot(:,:) = self%Emission%surf_emiss_out2(:,:,2)*OneOverDt
     IF(ASSOCIATED(emBioPropene))  emBioPropene(:,:) = self%Emission%surf_emiss_out2(:,:,3)*OneOverDt
     IF(ASSOCIATED(   emSoilNOx))     emSoilNOx(:,:) = self%Emission%surf_emiss_out2(:,:,4)*OneOverDt
     IF(self%Emission%do_ShipEmission) THEN
       IF(ASSOCIATED(  emShipHNO3))   emShipHNO3(:,:) = self%Emission%surf_emiss_out2(:,:,5)*OneOverDt
       IF(ASSOCIATED(    emShipO3))     emShipO3(:,:) = self%Emission%surf_emiss_out2(:,:,6)*OneOverDt
     END IF
   END IF

! EM_ emissions. Reverse in vertical and convert units if necessary.
! When lightning_opt = 0, the prescibed emissions are included in EM_NO.
! ----------------------------------------------------------------------
   IF(self%do_emission .AND. self%pr_emiss_3d .AND. self%numEM_Exports > 0) THEN

     ALLOCATE(var2dDBL(i1:i2,j1:j2), STAT=STATUS)
     VERIFY_(STATUS)

! For each EM_ export ...
! -----------------------
     Scan: DO n = 1, self%numEM_Exports

! Grab the molecular weight and find the pointer
! ----------------------------------------------
     fieldName = TRIM(self%EM_ExportNames(n))
     unitsName = TRIM(self%EM_ExportUnits(n))
     CALL getMW(TRIM(fieldName), ic, mw, rc)
     VERIFY_(rc)
     CALL ESMFL_StateGetPointerToData(expChem, EM_pointer, TRIM(fieldName), RC=STATUS)
     VERIFY_(STATUS)
     
     IsAssociated: IF(ASSOCIATED(EM_pointer)) THEN
! GCR emissions
! ------------
      IF(TRIM(fieldName) == "EM_GCR_NO") THEN
       ALLOCATE(var3dDBL(i1:i2,j1:j2,1:km),STAT=STATUS) 
       VERIFY_(STATUS)
       if (self%do_gcr) then
        CALL GET_GCR_EMISS ( var3dDBL, i1, i2, j1, j2, 1, km )  ! (molec/cm3/sec) (bottom up)

        SELECT CASE (unitsName)
        CASE ("molec cm-3 s-1")
         EM_pointer(:,:,km:1:-1) = var3dDBL(:,:,1:km)        
        CASE ("kg m-3 s-1")
         EM_pointer(:,:,km:1:-1) = var3dDBL(:,:,1:km) / MAPL_AVOGAD * mw * 1.e+6
        CASE ("kg m-2 s-1")
         EM_pointer(:,:,km:1:-1) = var3dDBL(:,:,1:km) / MAPL_AVOGAD * mw * 1.e+6  *  gridBoxThickness(:,:,1:km) 
        CASE DEFAULT
         PRINT *,TRIM(Iam)//": Modifications needed to export  "//TRIM(unitsName)//" for "//TRIM(fieldName)
         STATUS = -1
         VERIFY_(STATUS)
        END SELECT

       else        
        EM_pointer(:,:,km:1:-1) = 0.0 
       end if

       DEALLOCATE(var3dDBL)

! Lightning NO section
! --------------------
      ELSE IF(TRIM(fieldName) == "EM_LGTNO") THEN

       CALL Get_lightning_opt(self%Emission,lightning_opt)

       DO k = 1,km 
        kReverse = km-k+1

! Force kg m^{-3} s^{-1} for both prescribed and parameterized lightning NO
! -------------------------------------------------------------------------
             SELECT CASE (lightning_opt)
 	      CASE (0)
  	       EM_pointer(i1:i2,j1:j2,kReverse) = self%Emission%emissionArray(self%ic_NO_lgt)%pArray3D(i1:i2,j1:j2,k)/ &
	                                          (gridBoxThickness(:,:,k)*self%cellArea(:,:))
              CASE (1)
                EM_pointer(i1:i2,j1:j2,kReverse) = self%Emission%lightning_NO(i1:i2,j1:j2,k)
              CASE DEFAULT
	        EM_pointer(i1:i2,j1:j2,k) = 0.00
             END SELECT

! Convert units, if necessary
! ---------------------------
             SELECT CASE (unitsName)
              CASE ("kg m-2 s-1")
               EM_pointer(i1:i2,j1:j2,kReverse) = EM_pointer(i1:i2,j1:j2,kReverse)*gridBoxThickness(:,:,k)
              CASE ("kg m-3 s-1")
	        STATUS = 0
              CASE ("mol mol-1 s-1")
                var2dDBL(i1:i2,j1:j2) = MAPL_AIRMW/(AIRDENS(:,:,kReverse)*mw)
                EM_pointer(i1:i2,j1:j2,kReverse) = EM_pointer(i1:i2,j1:j2,kReverse)*var2dDBL(i1:i2,j1:j2)
              CASE DEFAULT
                PRINT *,TRIM(Iam)//": Modifications needed to export  "//TRIM(unitsName)//" for "//TRIM(fieldName)
                STATUS = -1 
                VERIFY_(STATUS)
             END SELECT

           END DO
      
!      Record the Overpass values   (note: AddExport done in GMIchem_GridCompMod.F90)
!      ------------------------------------------------------------------------------

           CURRENT_HMS = OVP_end_of_timestep_hms( CLOCK, OVP_RUN_DT )
!           IF(MAPL_AM_I_ROOT()) PRINT*,'GMI Emiss CURRENT_HMS =',CURRENT_HMS

!           EM_LGTNO overpass

           DATA_FOR_OVP_3D => EM_pointer

           CALL MAPL_GetPointer(expChem, OVP10_OUTPUT_3D, 'OVP10_EM_LGTNO', __RC__)
           CALL MAPL_GetPointer(expChem, OVP14_OUTPUT_3D, 'OVP14_EM_LGTNO', __RC__)

           CALL OVP_apply_mask( DATA_FOR_OVP_3D, OVP10_OUTPUT_3D, MASK_10AM, OVP_FIRST_HMS, CURRENT_HMS, K_EDGES=.FALSE., __RC__ )
           CALL OVP_apply_mask( DATA_FOR_OVP_3D, OVP14_OUTPUT_3D, MASK_2PM,  OVP_FIRST_HMS, CURRENT_HMS, K_EDGES=.FALSE., __RC__ )


! Emission exports other than lightning NO
! ----------------------------------------
         ELSE
      
           DO k = 1,km
             kReverse = km-k+1

! Convert units from kg m^{-2} per time step to specification
! -----------------------------------------------------------
        SELECT CASE (unitsName)
        CASE ("kg m-2 s-1")
         EM_pointer(:,:,kReverse) = self%Emission%emiss_3d_out(:,:,k,ic)*OneOverDt
        CASE ("kg m-3 s-1")
         EM_pointer(:,:,kReverse) = self%Emission%emiss_3d_out(:,:,k,ic)*OneOverDt/gridBoxThickness(:,:,k)
        CASE ("mol mol-1 s-1")
         var2dDBL(:,:) = OneOverDt*MAPL_AIRMW/(AIRDENS(:,:,kReverse)*mw*gridBoxThickness(:,:,k))
         EM_pointer(:,:,kReverse) = self%Emission%emiss_3d_out(:,:,k,ic)*var2dDBL(:,:)
        CASE DEFAULT
         PRINT *,TRIM(Iam)//": Modifications needed to export  "//TRIM(unitsName)//" for "//TRIM(fieldName)
         STATUS = -1
         VERIFY_(STATUS)
        END SELECT

       END DO

      END IF 

     END IF IsAssociated
    END DO Scan

!!!!!
! For CCMI

     CALL MAPL_GetPointer(expChem, EM_pointer, "BOX_HEIGHT", __RC__)
     IF(ASSOCIATED(EM_pointer)) EM_pointer(i1:i2,j1:j2,km:1:-1) = gridBoxThickness(i1:i2,j1:j2,1:km)

    ! To do units conversion for CCMI:  EM_species to EM_species_V
    ! Use this term to convert mol/mol/s to kg/m2/s in HISTORY:  EM_FIELD * AIRMASS * MW_species/MAPL_AIRMW
!!  ! previous:
!!  ! Use this term to convert mol/mol/s to kg/m2/s in HISTORY:  EM_FIELD * MOL_MOL_REVERSE * MW_species
!!  CALL MAPL_GetPointer(expChem, EM_pointer, "MOL_MOL_REVERSE", __RC__)
!!  IF(ASSOCIATED(EM_pointer)) EM_pointer(i1:i2,j1:j2,km:1:-1) = &
!!                               (AIRDENS(i1:i2,j1:j2,km:1:-1) * gridBoxThickness(i1:i2,j1:j2,1:km)) / MAPL_AIRMW

!! NATSAD and ICESAD
     call ESMF_StateGet (expChem, "gmiSAD", sadBundle, RC=STATUS )
     VERIFY_(STATUS)

     CALL MAPL_GetPointer(expChem, EM_pointer, "NATSAD", __RC__)
     IF(ASSOCIATED(EM_pointer)) THEN
       CALL obtainTracerFromBundle(sadBundle, ptr3D, INATSAD)
       EM_pointer(i1:i2,j1:j2,1:km) = ptr3D(i1:i2,j1:j2,1:km)
     END IF

     CALL MAPL_GetPointer(expChem, EM_pointer, "ICESAD", __RC__)
     IF(ASSOCIATED(EM_pointer)) THEN
       CALL obtainTracerFromBundle(sadBundle, ptr3D, IICESAD)
       EM_pointer(i1:i2,j1:j2,1:km) = ptr3D(i1:i2,j1:j2,1:km)
     END IF

     CALL MAPL_GetPointer(expChem, EM_pointer, "LBSSAD", __RC__)
     IF(ASSOCIATED(EM_pointer)) THEN
       CALL obtainTracerFromBundle(sadBundle, ptr3D, ILBSSAD)
       EM_pointer(i1:i2,j1:j2,1:km) = ptr3D(i1:i2,j1:j2,1:km)
     END IF

!! AIRMASS
     CALL MAPL_GetPointer(expChem, EM_pointer, "AIRMASS", __RC__)
     IF(ASSOCIATED(EM_pointer))   EM_pointer(i1:i2, j1:j2, 1:km) = AIRDENS(i1:i2, j1:j2, 1:km) * gridBoxThickness(i1:i2, j1:j2, km:1:-1)

!!!! >>>>>>>>>>>>>>>>  OVP
     IF(ASSOCIATED(EM_pointer)) THEN
       DATA_FOR_OVP_3D => EM_pointer

       CALL MAPL_GetPointer(expChem, OVP10_OUTPUT_3D, 'OVP10_AIRMASS', __RC__)
       CALL MAPL_GetPointer(expChem, OVP14_OUTPUT_3D, 'OVP14_AIRMASS', __RC__)

       CALL OVP_apply_mask( DATA_FOR_OVP_3D, OVP10_OUTPUT_3D, MASK_10AM, OVP_FIRST_HMS, CURRENT_HMS, K_EDGES=.FALSE., __RC__ )
       CALL OVP_apply_mask( DATA_FOR_OVP_3D, OVP14_OUTPUT_3D, MASK_2PM,  OVP_FIRST_HMS, CURRENT_HMS, K_EDGES=.FALSE., __RC__ )
     END IF
!!!! <<<<<<<<<<<<<<<<  OVP
!!!!!


! Clean up
! --------
     DEALLOCATE(var2dDBL, STAT=STATUS)
     VERIFY_(STATUS)

   END IF

   RETURN
   END SUBROUTINE FillExports

  SUBROUTINE getMW(name, i, mw, rc)
  CHARACTER(LEN=ESMF_MAXSTR) :: Iam = "getMW"
  CHARACTER(LEN=*), INTENT(IN) :: name
  INTEGER, INTENT(OUT) :: i, rc
  REAL, INTENT(OUT) :: mw
  INTEGER :: iii
  rc = 0
  SELECT CASE (TRIM(name))
   CASE ("EM_NO")
    mw = mw_data(  INO)
    i =   INO
   CASE ("EM_CO")
    mw = mw_data(  ICO)
    i =   ICO
   CASE ("EM_MEK")
    mw = mw_data( IMEK)
    i =  IMEK
   CASE ("EM_PRPE")
    mw = mw_data(IC3H6)
    i = IC3H6
   CASE ("EM_C2H6")
    mw = mw_data(IC2H6)
    i = IC2H6
   CASE ("EM_C3H8")
    mw = mw_data(IC3H8)
    i = IC3H8
   CASE ("EM_ALK4")
    mw = mw_data(IALK4)
    i = IALK4
   CASE ("EM_ALD2")
    mw = mw_data(IALD2)
    i = IALD2
   CASE ("EM_CH2O")
    mw = mw_data(ICH2O)
    i = ICH2O
   CASE ("EM_ACET")
    mw = mw_data(IACET)
    i = IACET
   CASE ("EM_CH4")
    mw = mw_data( ICH4)
    i =  ICH4
   CASE ("EM_LGTNO")
    mw = mw_data(  INO)
    i =  -1
   CASE ("EM_GCR_NO")
    mw = mw_data(  INO)
    i =  -1
   CASE DEFAULT
!... Find species in setkin_lchem.h list
    do iii=1,NSP
      if ( ESMF_UtilStringUpperCase(TRIM(name(4:))).eq.ESMF_UtilStringUpperCase(TRIM(lchemvar(iii))) ) then
        i = iii
        mw = mw_data(iii)
        RETURN
      endif
    enddo
    PRINT *,TRIM(Iam)//": Add "//TRIM(name)//" to molecular weight search list"
    STATUS = -1 
    VERIFY_(STATUS)
  END SELECT
  RETURN
  END SUBROUTINE getMW

!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  populateBundle
!
! !INTERFACE:

      subroutine populateBundle()
!
      implicit none
!
! !DESCRIPTION:
!
! !LOCAL VARIABLES:
      integer :: STATUS, numVars, ib, rc
      real(rPrec), pointer, dimension(:,:)  :: ptr2D
      type(ESMF_FieldBundle)                :: surfEmissBundle
      character(len=ESMF_MAXSTR), parameter :: IAm = "populateBundle"
!
!EOP
!--------------------------------------------------------------------------------

      allocate(ptr2D(i1:i2, j1:j2))

      !================================
      ! For the Surface Emission Bundle
      !================================

      call ESMF_StateGet(expChem, "surfEmissForChem", surfEmissBundle, rc=STATUS)
      VERIFY_(STATUS)

      call ESMF_FieldBundleGet(surfEmissBundle, fieldCount=numVars, rc=STATUS)
      VERIFY_(STATUS)
      _ASSERT(numVars == NSP,'needs informative message')

      do ib = 1, numVars
         ptr2D(:,:) = surfEmissForChem(:,:,ib)
         call updateTracerToBundle(surfEmissBundle, ptr2D, ib)
      end do

      deallocate (ptr2D)

      return

      end subroutine populateBundle

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

   CHARACTER(LEN=255) :: IAm
  
   rc=0
   IAm="FindPointers"

!  Pointers to imports
!  -------------------
   CALL MAPL_GetPointer(impChem,   cn_prcp,   'CN_PRCP', RC=STATUS)
   VERIFY_(STATUS)
   CALL MAPL_GetPointer(impChem,     tprec,     'TPREC', RC=STATUS)
   VERIFY_(STATUS)
   CALL MAPL_GetPointer(impChem,       lwi,       'LWI', RC=STATUS)
   VERIFY_(STATUS)
   CALL MAPL_GetPointer(impChem, frlandice, 'FRLANDICE', RC=STATUS)
   VERIFY_(STATUS)
   CALL MAPL_GetPointer(impChem,    snowdp,    'SNOWDP', RC=STATUS)
   VERIFY_(STATUS)
   CALL MAPL_GetPointer(impChem,       T2m,       'T2M', RC=STATUS)
   VERIFY_(STATUS)
   CALL MAPL_GetPointer(impChem,      zpbl,      'ZPBL', RC=STATUS)
   VERIFY_(STATUS)
   CALL MAPL_GetPointer(impChem,      u10m,      'U10M', RC=STATUS)
   VERIFY_(STATUS)
   CALL MAPL_GetPointer(impChem,      v10m,      'V10M', RC=STATUS)
   VERIFY_(STATUS)
   CALL MAPL_GetPointer(impChem,     ustar,     'USTAR', RC=STATUS)
   VERIFY_(STATUS)
   CALL MAPL_GetPointer(impChem,       z0h,       'Z0H', RC=STATUS)
   VERIFY_(STATUS)
   CALL MAPL_GetPointer(impChem,   swndsrf,   'SWNDSRF', RC=STATUS)
   VERIFY_(STATUS)
   CALL MAPL_GetPointer(impChem,     cldtt,     'CLDTT', RC=STATUS)
   VERIFY_(STATUS)
   CALL MAPL_GetPointer(impChem,      wet1,      'WET1', RC=STATUS)
   VERIFY_(STATUS)
   CALL MAPL_GetPointer(impChem,     dfpar,     'DFPAR', RC=STATUS)
   VERIFY_(STATUS)
   CALL MAPL_GetPointer(impChem,     drpar,     'DRPAR', RC=STATUS)
   VERIFY_(STATUS)
  
   CALL MAPL_GetPointer(impChem,   airdens, 'AIRDENS', RC=STATUS)
   VERIFY_(STATUS)
   CALL MAPL_GetPointer(impChem,       ple,	'PLE', RC=STATUS)
   VERIFY_(STATUS)
   CALL MAPL_GetPointer(impChem,         Q,       'Q', RC=STATUS)
   VERIFY_(STATUS)
   CALL MAPL_GetPointer(impChem,	 T,	  'T', RC=STATUS)
   VERIFY_(STATUS)
   CALL MAPL_GetPointer(impChem,       zle,	'ZLE', RC=STATUS)
   VERIFY_(STATUS)
   CALL MAPL_GetPointer(impChem,   cnv_mfc, 'CNV_MFC', RC=STATUS)
   VERIFY_(STATUS)

!  Export state pointers
!  ---------------------
   CALL MAPL_GetPointer(expChem,    emIsopSfc,   'EMISOPSFC', RC=STATUS)
   VERIFY_(STATUS)
   CALL MAPL_GetPointer(expChem,        emNOx,       'EMNOX', RC=STATUS)
   VERIFY_(STATUS)
   CALL MAPL_GetPointer(expChem,      emMonot,     'EMMONOT', RC=STATUS)
   VERIFY_(STATUS)
   CALL MAPL_GetPointer(expChem,  emBioCOMeth, 'EMBIOCOMETH', RC=STATUS)
   VERIFY_(STATUS)
   CALL MAPL_GetPointer(expChem, emBioCOMonot,'EMBIOCOMONOT', RC=STATUS)
   VERIFY_(STATUS)
   CALL MAPL_GetPointer(expChem, emBioPropene,'EMBIOPROPENE', RC=STATUS)
   VERIFY_(STATUS)
   CALL MAPL_GetPointer(expChem,    emSoilNOx,   'EMSOILNOX', RC=STATUS)
   VERIFY_(STATUS)
   CALL MAPL_GetPointer(expChem,   emShipHNO3,  'EMSHIPHNO3', RC=STATUS)
   VERIFY_(STATUS)
   CALL MAPL_GetPointer(expChem,     emShipO3,    'EMSHIPO3', RC=STATUS)
   VERIFY_(STATUS)

!  Validation
!  ----------
   Validate: IF(self%verbose) THEN
    IF(MAPL_AM_I_ROOT( ))PRINT *,TRIM(IAm),": Input ..."
    CALL MAPL_MaxMin('CN_PRCP:', cn_prcp)
    CALL MAPL_MaxMin('TPREC:', tprec)
    CALL MAPL_MaxMin('LWI:', lwi)
    CALL MAPL_MaxMin('FRLANDICE:', frlandice)
    CALL MAPL_MaxMin('SNOWDP:', snowdp)
    CALL MAPL_MaxMin('T2M:', T2m)
    CALL MAPL_MaxMin('ZPBL:', zpbl)
    CALL MAPL_MaxMin('Q:', Q)
    CALL MAPL_MaxMin('T:', T)
    CALL MAPL_MaxMin('ZLE:', zle)
    CALL MAPL_MaxMin('PLE (hPa):', ple)
    CALL MAPL_MaxMin('AIRDENS:', airdens)
    CALL MAPL_MaxMin('CNV_MFC:', cnv_mfc)

    CALL MAPL_MaxMin('U10M:', u10m)
    CALL MAPL_MaxMin('V10M:', v10m)
    CALL MAPL_MaxMin('USTAR:', ustar)
    CALL MAPL_MaxMin('SURF_ROUGH:', z0h)
    CALL MAPL_MaxMin('SWNDSRF:', swndsrf)
    CALL MAPL_MaxMin('CLDTT:', cldtt)
    CALL MAPL_MaxMin('WET1:', wet1)
    CALL MAPL_MaxMin('DFPAR:', dfpar)
    CALL MAPL_MaxMin('DRPAR:', drpar)

   END IF Validate

! Ship Emisssions
! ---------------
   if (self%Emission%do_ShipEmission) then
      jNO2val_phot(i1:i2,j1:j2) = bxx%qa(bxx%reg%nq)%data3d(i1:i2,j1:j2,km-1)
   end if

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
  pctm2(i1:i2,j1:j2) = ple(i1:i2,j1:j2,km)*Pa2hPa                           ! Pa               hPa
  fracCloudCover(i1:i2,j1:j2) = cldtt(i1:i2,j1:j2)                          ! fraction
  surf_rough(i1:i2,j1:j2) = z0h(i1:i2,j1:j2)                                ! m
  radswg(i1:i2,j1:j2) = swndsrf(i1:i2,j1:j2)                                ! w m^{-2}
  tenMeterU(i1:i2,j1:j2) = u10m(i1:i2,j1:j2)                                ! m s^{-1}
  tenMeterV(i1:i2,j1:j2) = v10m(i1:i2,j1:j2)                                ! m s^{-1}
  frictionVelocity(i1:i2,j1:j2) = ustar(i1:i2,j1:j2)                        ! m s^{-1}
  con_precip(i1:i2,j1:j2) = cn_prcp(i1:i2,j1:j2)*secPerDay                  ! kg m^{-2}s^{-1}   mm d^{-1}
  tot_precip(i1:i2,j1:j2) = tprec(i1:i2,j1:j2)*secPerDay                    ! kg m^{-2}s^{-1}   mm d^{-1}
  pbl(i1:i2,j1:j2) = zpbl(i1:i2,j1:j2)                                      ! m
  soilWetness(i1:i2,j1:j2) = wet1(i1:i2,j1:j2)                              ! fraction
  diffusePar(i1:i2,j1:j2) = dfpar(i1:i2,j1:j2)                              ! w m^{-2}
  directPar(i1:i2,j1:j2) = drpar(i1:i2,j1:j2)                               ! w m^{-2}

! Layer means                                                               GEOS-5 Units       GMI Units
! -----------                                                               ------------       -------------
  DO k=1,km
    kReverse = km-k+1                                                        ! Lid-to-surf      Surf-to-lid
    press3c(i1:i2,j1:j2,kReverse) = pl(i1:i2,j1:j2,k)*Pa2hPa                 ! Pa               hPa
    kel(i1:i2,j1:j2,kReverse) = T(i1:i2,j1:j2,k)                             ! K
  END DO

! Layer edges                                                               GEOS-5 Units       GMI Units
! -----------                                                               ------------       -------------
  DO k=0,km
    kReverse = km-k
    press3e(i1:i2,j1:j2,kReverse) = ple(i1:i2,j1:j2,k)*Pa2hPa                ! Pa               hPa
    height3e(i1:i2,j1:j2,kReverse) = zle(i1:i2,j1:j2,k)                      ! m                m
  END DO

  DO k=0,km-1
    kReverse = km-k
    cmf(i1:i2,j1:j2,kReverse) = cnv_mfc(i1:i2,j1:j2,k)                       ! kg m^{-2}s^{-1}
  END DO

! Incoming land-water-ice flag in GEOS original format:
!   0=water 1=land 2=ice
! Add snow flag (3) for Dry deposition
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

!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  MonitorT2M
!
! !INTERFACE:

  SUBROUTINE MonitorT2M(rc)
  
  IMPLICIT NONE

  INTEGER, INTENT(OUT) :: rc
!
! !DESCRIPTION:
!
!  Obtain the current and previous 15 days average T2M.
!
!  T2M15d is a 3D array that is not entirely filled and is concatenated to the 
!  XX (non-transported) species.  Each "layer" of T2M15d contains:

!     Index   Contents
!   --------  --------
!          1  Average T2M over the previous 15 whole days
!       2-15  Respective days' average T2M, with 15 being the most recent
!             and 2 being the least recent.
!         16  Today's average, which is being computed
!   17->km-2  Unused, but available for up to km-2 days
!       km-1  Place holder jNO2val
!         km  Place holder for most recent valid tropopause pressures (Pa)

!  When the day number changes, the values in position 16 must be 
!  moved to 15, 15 to 14, 14 to 13, and so on.  After this age-off is
!  complete, position 1 is replaced with the average of positions 1-15.
!  Position 1 is referenced as the 15-day average T2M.

!  USAGE NOTES:

!  Only whole calendar days are considered in the running averages, and 
!   age-off occurs on the first pass of the new calendar day, which occurs
!   at heartBeat seconds after midnight UTC.  Note that for 0:00 UTC restarts,
!   age-off occurs after restart, not before the end of the previous segment.

!  When the restart is not at 0:00 UTC, the running average for the current
!   calendar day continues to be computed using the value from the internal
!   restart file.

!  WARNING: Do not overwrite "layer" km, as it contains the most recent 
!   valid tropopause pressures (Pa). This is a kluge, but enables (1) the
!   reproducibility of GMIchem across various layouts, including the cubed
!   sphere, (2) eliminates the occasional MAPL_UNDEFs in TROPP, and (3) 
!   allows the information to be included the w_c bundle.

!EOP
!---------------------------------------------------------------------------

  CHARACTER(LEN=255) :: IAm
  REAL :: avgFactor
  INTEGER :: hms, hr, mi, r, sc

  rc = 0
  IAm = "MonitorT2M"
  avgFactor = tdt/86400.00

! Express tdt in hms format
! -------------------------
  hr = tdt/3600
  r = tdt-hr*3600
  mi = r/60
  sc = r-mi*60
  hms = 10000*hr+100*mi+sc

! Index of T2M15d in the non-transported species bundle
! -----------------------------------------------------
  i = bxx%reg%nq

! Sanity check
! ------------
  IF(km < 18) THEN
   ! Unable to do 15 day avg and store the other items we need
    IF( MAPL_AM_I_ROOT() ) THEN
      PRINT *,"GMICHEM::"//TRIM(IAm)//": Unable to perform 15 days average of T2M."
      PRINT *," "
    END IF
    rc = 51
  END IF

! Functions for startup and at change of day
! ------------------------------------------
  IF(nhms < hms) THEN

    IF( MAPL_AM_I_ROOT() ) THEN
      PRINT *," "
      PRINT *,"GMICHEM::"//TRIM(IAm)//": Doing T2M15d age-off and recalculation"  
      PRINT *," "
    END IF

! Age-off
! -------
   DO k = 1,15
    bxx%qa(i)%data3d(i1:i2,j1:j2,k) = bxx%qa(i)%data3d(i1:i2,j1:j2,k+1)
   END DO

! Calculate average T2M for previous 15 whole days
! ------------------------------------------------
   DO k = 2,15
    bxx%qa(i)%data3d(i1:i2,j1:j2,1) = bxx%qa(i)%data3d(i1:i2,j1:j2,1) + bxx%qa(i)%data3d(i1:i2,j1:j2,k)
   END DO
   bxx%qa(i)%data3d(i1:i2,j1:j2,1) = bxx%qa(i)%data3d(i1:i2,j1:j2,1)/15.00

! Initialize to zero for day that is just starting
! ------------------------------------------------
   bxx%qa(i)%data3d(:,:,16) = 0.00

  END IF

! Otherwise, keep running average for current day
! -----------------------------------------------
  bxx%qa(i)%data3d(i1:i2,j1:j2,16) = bxx%qa(i)%data3d(i1:i2,j1:j2,16) + T2m(i1:i2,j1:j2)*avgFactor
    
  RETURN
 END SUBROUTINE MonitorT2M

 END SUBROUTINE GmiEmiss_GridCompRun

!-------------------------------------------------------------------------
!NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1, GEOS/DAS!
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  GmiEmiss_GridCompFinalize
!
! !INTERFACE:
!

   SUBROUTINE GmiEmiss_GridCompFinalize ( self, impChem, expChem, &
                                     nymd, nhms, cdt, rc )

  USE gcr_mod,                       ONLY : Finalize_GCR
  IMPLICIT none

! !INPUT/OUTPUT PARAMETERS:

   TYPE(GmiEmiss_GridComp), INTENT(inout) :: self ! Grid Component

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

   CHARACTER(LEN=*), PARAMETER :: IAm = 'GmiEmiss_GridCompFinalize'
   INTEGER :: STATUS
   rc=0
   DEALLOCATE(self%cellArea, self%lonRad, self%latRad, STAT=STATUS)
   VERIFY_(STATUS)
   DEALLOCATE(self%EM_ExportNames, self%EM_ExportUnits, STAT=STATUS)
   VERIFY_(STATUS)

!  Free the masks
!  --------------------
   deallocate ( MASK_10AM, MASK_2PM, stat = STATUS )
   VERIFY_(STATUS)

   if (self%do_gcr) call Finalize_GCR()

   RETURN

 END SUBROUTINE GmiEmiss_GridCompFinalize
  
 END module GmiEmiss_GCCMod
