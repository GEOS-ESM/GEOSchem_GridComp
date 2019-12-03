#include "MAPL_Generic.h"
!-------------------------------------------------------------------------
!NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1, GEOS/DAS!
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  GmiPhotolysis_GCCMod --- GMI Photolysis Grid Component Class
!
! Grid Component class for the GMI Photolysis.
!
! !INTERFACE:
!
   MODULE  GmiPhotolysis_GCCMod

! !USES:

   USE ESMF
   USE MAPL
   USE Chem_Mod         ! Chemistry Base Class
   USE Chem_UtilMod

   USE GmiSpcConcentrationMethod_mod, ONLY : t_SpeciesConcentration
   USE GmiGrid_mod,                   ONLY : t_gmiGrid
   USE GmiTimeControl_mod,            ONLY : t_GmiClock
   USE GmiESMFrcFileReading_mod,      ONLY : rcEsmfReadLogical
   USE GmiESMFrcFileReading_mod,      ONLY : rcEsmfReadTable
   use GmiArrayBundlePointer_mod,     ONLY : t_GmiArrayBundle, CleanArrayPointer
   use GmiFieldBundleESMF_mod,        ONLY : updateTracerToBundle
   use GmiFieldBundleESMF_mod,        ONLY : addTracerToBundle
   use GmiStateFieldESMF_mod,         ONLY : setDataToStateField
   use GmiStateFieldESMF_mod,         ONLY : initDataInStateField
   use GmiSwapSpeciesBundlesMod,      ONLY : SwapSpeciesBundles, speciesReg_for_CCM
   USE GmiFastJX_includeMod,          ONLY : t_fastJXbundle
   USE GmiFastJX_ParametersMod,       ONLY : NW

   IMPLICIT NONE
   INTEGER, PARAMETER :: DBL = KIND(0.00D+00)

#include "setkin_par.h"
#include "GmiParameters.h"
#include "gmi_phys_constants.h"
#include "setkin_mw.h"
#include "setkin_lchem.h"
#include "gmi_AerDust_const.h"
#include "gmi_time_constants.h"
! !TYPES:

   PRIVATE
   PUBLIC  GmiPhotolysis_GridComp       ! The GMI object 

! !PUBLIC MEMBER FUNCTIONS:

   PUBLIC  GmiPhotolysis_GridCompInitialize
   PUBLIC  GmiPhotolysis_GridCompRun
   PUBLIC  GmiPhotolysis_GridCompFinalize

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
!  01Jul2010 Kouatchou Creation of Chemistry Grid Component class
!  01Jun2015 Nielsen   ExtData replaces Chem_UtilMPread
!  26May2016 Oman      Modified FastJX 6.5 initialization;
!                      also scale tau values by cloud fraction
!
!EOP
!-------------------------------------------------------------------------

  TYPE GmiPhotolysis_GridComp
   CHARACTER(LEN=255) :: name = "GMI Photolysis"

! Does the gmichem_import_rst file exist?  If not, some import states will be
! unfilled, and the species will need to "freewheel" through the first time step.
! -------------------------------------------------------------------------------
   LOGICAL :: gotImportRst

! Set BCRealTime = .TRUE. when boundary conditions 
! must be for exact year of current calendar date.
! -------------------------------------------------
   LOGICAL :: BCRealTime

! Switches for coupling with GOCART dust and aerosols
! ---------------------------------------------------
   LOGICAL :: usingGOCART_BC
   LOGICAL :: usingGOCART_DU
   LOGICAL :: usingGOCART_OC
   LOGICAL :: usingGOCART_SS
   LOGICAL :: usingGOCART_SU

! Perhaps GMICHEM is the AERO_PROVIDER
! ------------------------------------
   LOGICAL :: AM_I_AERO_PROVIDER
   CHARACTER(LEN=255) :: aeroProviderName

! Various switches
! ----------------
   LOGICAL :: pr_diag
   LOGICAL :: do_synoz

! Dimensions
! ----------
   INTEGER :: i1, i2, im, j1, j2, jm, km

! Useful character strings
! ------------------------
   CHARACTER(LEN=255) :: chem_mecha
   character(len=ESMF_MAXSTR) :: qj_labels(MAX_NUM_QJ) ! qj (photolysis) string labels

! Variables for Ship Emissions
! ----------------------------
   INTEGER :: jno2num
   LOGICAL :: do_ShipEmission
   LOGICAL :: do_semiss_inchem

! Surface area of grid cells
! --------------------------
   REAL(KIND=DBL), POINTER :: cellArea(:,:)

! Longitudes and latitudes (radians)
! ----------------------------------
   REAL, ALLOCATABLE, DIMENSION(:,:) :: lonRad
   REAL, ALLOCATABLE, DIMENSION(:,:) :: latRad

! Extra diagnostics
! -----------------
   LOGICAL :: verbose

! Map GMI species indices to CCM indices
! --------------------------------------
   INTEGER, POINTER :: mapSpecies(:)

! This is to store a generic name for each qjgmi entry
! ----------------------------------------------------
    character(len=MAX_STRING_LENGTH), pointer :: qjName(:) => null()

    character (len=MAX_LENGTH_VAR_NAME)  :: qj_var_name
    integer             :: phot_opt
    integer             :: fastj_opt
    logical             :: do_clear_sky
    real*8              :: fastj_offset_sec

    real*8              :: synoz_threshold
    integer             :: chem_mask_klo
    integer             :: chem_mask_khi

    real*8              :: qj_init_val
    integer             :: qj_timpyr
    character (len=MAX_LENGTH_FILE_NAME) :: qj_infile_name
    integer             :: num_qjs
    integer             :: num_qjo
    integer             :: sfalbedo_opt
    INTEGER             :: jNOindex
    REAL                :: jNOamp
    real*8              :: saldif_init_val
    real*8              :: saldir_init_val
    real*8              :: sasdif_init_val
    real*8              :: sasdir_init_val
            ! surface albedo data for diffuse nearIR (fraction 0-1)
    real*8, pointer     :: saldif_data(:,:,:) => null()
            ! surface albedo data for direct  nearIR (fraction 0-1)
    real*8, pointer     :: saldir_data(:,:,:) => null()
            ! surface albedo data for diffuse uv/vis (fraction 0-1)
    real*8, pointer     :: sasdif_data(:,:,:) => null()
            ! surface albedo data for direct  uv/vis (fraction 0-1)
    real*8, pointer     :: sasdir_data(:,:,:) => null()
    character (len=MAX_LENGTH_FILE_NAME) :: sfalbedo_infile_name
    integer             :: uvalbedo_opt
    real*8              :: uvalbedo_init_val
    real*8, pointer     :: uvalbedo_data(:,:,:) => null()
    character (len=MAX_LENGTH_FILE_NAME) :: uvalbedo_infile_name
    character (len=MAX_LENGTH_FILE_NAME) :: cross_section_file
    character (len=MAX_LENGTH_FILE_NAME) :: aerosolOpticalData_file
    character (len=MAX_LENGTH_FILE_NAME) :: rate_file
    character (len=MAX_LENGTH_FILE_NAME) :: T_O3_climatology_file
    character (len=MAX_LENGTH_FILE_NAME) :: scattering_data_file
    logical             :: do_solar_cycle
    character (len=MAX_LENGTH_FILE_NAME) :: sc_infile_name
    logical             :: do_ozone_inFastJX

    logical             :: pr_qj_o3_o1d
    logical             :: pr_qj_opt_depth

    integer             :: AerDust_Effect_opt
    logical             :: do_AerDust_Calc
    character (len=MAX_LENGTH_FILE_NAME) :: AerDust_infile_name
    character (len=MAX_LENGTH_FILE_NAME) :: Dust_infile_name
    character (len=MAX_LENGTH_FILE_NAME) :: Aerosol_infile_name

    type (t_GmiArrayBundle), pointer :: qjgmi(:) => null()
    real*8, pointer     :: overheadO3col(:,:,:)   => null()
    real*8, pointer     :: tArea    (:,:,:,:) => null()
    real*8, pointer     :: odAer    (:,:,:,:) => null()
    real*8, pointer     :: odMdust  (:,:,:,:) => null()
    real*8, pointer     :: eRadius  (:,:,:,:) => null()
    real*8, pointer     :: optDepth (:,:,:,:) => null()
    real*8, pointer     :: dust     (:,:,:,:) => null()
    real*8, pointer     :: wAersl   (:,:,:,:) => null()
    real*8, pointer     :: dAersl   (:,:,:,:) => null()


! Component derived type declarations
! -----------------------------------
   TYPE(t_gmiGrid   )		:: gmiGrid
   TYPE(t_GmiClock  )           :: gmiClock
   TYPE(t_SpeciesConcentration) :: SpeciesConcentration
 
   TYPE(t_fastJXbundle)   :: JXbundle

  END TYPE GmiPhotolysis_GridComp

CONTAINS

!-------------------------------------------------------------------------
!NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1, GEOS/DAS!
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  GmiPhotolysis_GridCompInitialize --- Initialize GmiPhotolysis_GridComp
!
! !INTERFACE:
!

   SUBROUTINE GmiPhotolysis_GridCompInitialize( self, w_c, impChem, expChem, nymd, nhms, &
                                      tdt, rc )

   USE GmiSpcConcentrationMethod_mod, ONLY : InitializeSpcConcentration
   USE GmiGrid_mod,		      ONLY : InitializeGmiGrid
   USE GmiTimeControl_mod,            ONLY : Set_curGmiDate, Set_curGmiTime
   USE GmiTimeControl_mod,            ONLY : Set_begGmiDate, Set_begGmiTime
   USE GmiTimeControl_mod,            ONLY : Set_numTimeSteps
   use GmiCheckNamelistFile_mod, only : CheckNamelistOptionRange
   use fastj           , only : InitializeFastj
   use fast_JX         , only : InitializeFastJX
   use Fast_JX53b      , only : InitializeFastJX53b
   use Fast_JX53c      , only : InitializeFastJX53c
   use fastJX65_mod,     only : initializeFastJX65
   USE ReadSolarCycle_mod,            ONLY : readSolarCycleData

   IMPLICIT none
   INTEGER, PARAMETER :: DBL = KIND(0.00D+00)

#     include "phot_lookup.h"
#     include "phot_monthly.h"
!
! !INPUT PARAMETERS:
!
   TYPE(Chem_Bundle), INTENT(in) :: w_c                ! Chemical tracer fields, delp, +
   INTEGER, INTENT(IN) :: nymd, nhms		       ! Time from AGCM
   REAL,    INTENT(IN) :: tdt			       ! Chemistry time step (secs)
!
! !OUTPUT PARAMETERS:
!
   TYPE(GmiPhotolysis_GridComp), INTENT(INOUT)  :: self      ! Grid Component
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

   CHARACTER(LEN=*), PARAMETER :: IAm = 'GmiPhotolysis_GridCompInitialize'
   CHARACTER(LEN=255) :: rcfilen = 'GMI_GridComp.rc'
   CHARACTER(LEN=255) :: kineticsTextFile
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
   INTEGER :: numSpecies
   INTEGER :: LogicalUnitNum
   INTEGER :: inyr,inmon,iscyr
   REAL, dimension(2628)     :: s_cycle_dates    ! 2628 months : 1882 - 2100
   REAL, dimension(NW,2628)  :: s_cycle          ! 2628 months : 1882 - 2100
   REAL(r8), POINTER         :: fjx_solar_cycle_param(:)


   INTEGER :: loc_proc, locGlobProc, commu_slaves
   LOGICAL :: one_proc, rootProc
   LOGICAL :: exists,open,found
   
   REAL :: qmin, qmax, tokgCPerBox

   real(rPrec), pointer :: var(:,:,:)
   type(ESMF_FieldBundle)      :: qjBundle
   type(ESMF_FieldBundle)      :: tAreaBundle
   type(ESMF_FieldBundle)      :: eRadiusBundle 
   integer                     :: numVars, ib
   character (len=4) :: binName
   character(len=ESMF_MAXSTR) :: varName
   real(r8)  :: hugeReal
   integer   :: badIndex = -9999
   integer   :: IXj
   REAL, POINTER, DIMENSION(:,:) :: jNO2val_phot

! Grid cell area can be set by initialize
! ---------------------------------------
   REAL, POINTER, DIMENSION(:,:) :: cellArea

! Work space
! ----------
   REAL, ALLOCATABLE :: var2D(:,:)
   REAL, ALLOCATABLE :: var3D(:,:,:)

   self%name = 'GMI Photolysis'

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
         PRINT *,"Starting Reading the GMI Resource File for Photolysis"
      ENDIF

      gmiConfigFile = ESMF_ConfigCreate(rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigLoadFile(gmiConfigFile, TRIM(rcfilen), rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(gmiConfigFile, kineticsTextFile, &
     &                label   = "kineticsTextFile:", &
     &                default = ' ', rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(gmiConfigFile, importRestartFile, &
     &                label   = "importRestartFile:", &
     &                default = ' ', rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(gmiConfigFile, self%chem_mecha, &
     &                label   = "chem_mecha:", &
     &                default = 'strat_trop', rc=STATUS )
      VERIFY_(STATUS)

      !------------------------------
      ! Emission related variables
      !------------------------------

      call rcEsmfReadLogical(gmiConfigFile, self%do_synoz, &
     &           "do_synoz:", default=.false., rc=STATUS)
      VERIFY_(STATUS)

      call rcEsmfReadLogical(gmiConfigFile, self%do_semiss_inchem, &
     &           "do_semiss_inchem:", default=.false., rc=STATUS)
      VERIFY_(STATUS)

      call rcEsmfReadLogical(gmiConfigFile, self%do_ShipEmission, &
     &           "do_ShipEmission:", default=.false., rc=STATUS)
      VERIFY_(STATUS)

      !------------------------------
      ! Diagnostics related variables
      !------------------------------

      call rcEsmfReadLogical(gmiConfigFile, self%pr_diag, &
     &           "pr_diag:", default=.false., rc=STATUS)
      VERIFY_(STATUS)

      call rcEsmfReadLogical(gmiConfigFile, self%verbose, &
     &           "verbose:", default=.false., rc=STATUS)
      VERIFY_(STATUS)

      !-------------------------------------------
      ! Should BC files have current date and time?
      ! Useful for mission support and replays.
      !--------------------------------------------
      
      call rcEsmfReadLogical(gmiConfigFile, self%BCRealTime, &
     &           "BCRealTime:", default=.false., rc=STATUS)
      VERIFY_(STATUS)
      
      !-----------------------------
      ! Photolysis Related Variables
      !-----------------------------

!     -----------------------------------------------------
!     phot_opt
!       0:  no photolysis
!       1:  set all qj values to qj_init_val
!       2:  read in qj values
!       3:  use fastj routine (for fastJ, fastJx, fastJx53b)
!           This option should be combined with fastj_opt.
!       4:  lookup table for qj (Kawa style)
!       5:  lookup table for qj (Kawa style) +
!           use ozone climatology for column ozone calc.
!       6:  calculate from table and Gmimod data (Quadchem)
!       7:  read in qj values (2-D, 12 months)
!       8:  use fast-JX routine (troposphere/stratosphere)
!     -----------------------------------------------------

      call ESMF_ConfigGetAttribute(gmiConfigFile, self%phot_opt, &
     &                label   = "phot_opt:", &
     &                default = 1, rc=STATUS )
      VERIFY_(STATUS)

!     -----------------------------------------------------
!     fastj_opt: set when phot_opt=3
!     0: for fastJ
!     1: for fastJx
!     2: for fastJx53b
!     2: for fastJx53c
!    -----------------------------------------------------

      call ESMF_ConfigGetAttribute(gmiConfigFile, self%fastj_opt, &
     &                label   = "fastj_opt:", &
     &                default = 1, rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(gmiConfigFile, self%do_clear_sky, "do_clear_sky:", &
     &                       default=.true., rc=STATUS)

      call ESMF_ConfigGetAttribute(gmiConfigFile, self%fastj_offset_sec, &
     &                label   = "fastj_offset_sec:", &
     &                default = 0.0d0, rc=STATUS )
      VERIFY_(STATUS)

      CALL ESMF_ConfigGetAttribute(gmiConfigFile, self%jNOamp, &
     &          LABEL="Prather_jNO_factor:", DEFAULT=1.00, RC=STATUS )
      VERIFY_(STATUS)

      found = .FALSE.
      i = 1
      DO WHILE (.NOT. found .AND. i <= NUM_J)
        IF ( TRIM(lqjchem(i)) == 'NO + hv = N + O' ) THEN
          self%jNOindex = i
          found = .TRUE.
        END IF
       i = i+1
      END DO
      IF ( .NOT. found ) THEN
        STATUS = 123
        VERIFY_(STATUS)
      END IF

!... do solar cycle in incoming solar flux?

      CALL rcEsmfReadLogical(gmiConfigFile, self%do_solar_cycle, "do_solar_cycle:", &
     &                       default=.false., rc=STATUS)
      VERIFY_(STATUS)

      CALL ESMF_ConfigGetAttribute(gmiConfigFile, self%sc_infile_name, &
     &                label   = "sc_infile_name:", &
     &                default = '', rc=STATUS )
      VERIFY_(STATUS)

!     ---------
!     qj / qqj:
!     ---------

      call ESMF_ConfigGetAttribute(gmiConfigFile, self%qj_init_val, &
     &                label   = "qj_init_val:", &
     &                default = 1.0d-30, rc=STATUS )
      VERIFY_(STATUS)

      ! sets of photolysis per year (1 => yearly, 12 => monthly)

      call ESMF_ConfigGetAttribute(gmiConfigFile, self%qj_timpyr, &
     &                label   = "qj_timpyr:", &
     &                default = 1, rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(gmiConfigFile, self%qj_infile_name, &
     &                label   = "qj_infile_name:", &
     &                default = '', rc=STATUS )
      VERIFY_(STATUS)

!     -------
!     albedo:
!     -------

!     ----------------------------------------------------------------------
!     sfalbedo_opt
!       0:  no sfalbedo
!       1:  set each type of sfalbedo to an intial value
!       2:  read in monthly sfalbedo values from a NetCDF file
!       3:  read in values of four types of surface albedo from the met data
!     ----------------------------------------------------------------------

      call ESMF_ConfigGetAttribute(gmiConfigFile, self%sfalbedo_opt, &
     &                label   = "sfalbedo_opt:", &
     &                default = 0, rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(gmiConfigFile, self%saldif_init_val, &
     &                label   = "saldif_init_val:", &
     &                default = 0.1d0, rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(gmiConfigFile, self%saldir_init_val, &
     &                label   = "saldir_init_val:", &
     &                default = 0.1d0, rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(gmiConfigFile, self%sasdif_init_val, &
     &                label   = "sasdif_init_val:", &
     &                default = 0.1d0, rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(gmiConfigFile, self%sasdir_init_val, &
     &                label   = "sasdir_init_val:", &
     &                default = 0.1d0, rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(gmiConfigFile, self%sfalbedo_infile_name, &
     &                label   = "sfalbedo_infile_name:", &
     &                default = '', rc=STATUS )
      VERIFY_(STATUS)

!     --------------------------------------------------------
!     uvalbedo_opt
!       0:  no uvalbedo
!       1:  set all uvalbedo values to uvalbedo_init_val
!       2:  read in monthly uvalbedo values from an ASCII file
!       3:  read in surface albedo values from the met data
!     --------------------------------------------------------

      call ESMF_ConfigGetAttribute(gmiConfigFile, self%uvalbedo_opt, &
     &                label   = "uvalbedo_opt:", &
     &                default = 0, rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(gmiConfigFile, self%uvalbedo_init_val, &
     &                label   = "uvalbedo_init_val:", &
     &                default = 0.1d0, rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(gmiConfigFile, self%uvalbedo_infile_name, &
     &                label   = "uvalbedo_infile_name:", &
     &                default = '', rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(gmiConfigFile, self%cross_section_file, &
     &                label   = "cross_section_file:", &
     &                default = '', rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(gmiConfigFile, self%aerosolOpticalData_file, &
     &                label   = "aerosolOpticalData_file:", &
     &                default = '', rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(gmiConfigFile, self%rate_file, &
     &                label   = "rate_file:", &
     &                default = '', rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(gmiConfigFile, self%T_O3_climatology_file, &
     &                label   = "T_O3_climatology_file:", &
     &                default = '', rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(gmiConfigFile, self%scattering_data_file, &
     &                label   = "scattering_data_file:", &
     &                default = '', rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(gmiConfigFile, self%do_ozone_inFastJX, &
     &              "do_ozone_inFastJX:", default=.false., rc=STATUS)
      VERIFY_(STATUS)

      !=================================================================
      ! do_AerDust_Calc is used to detrmine if aerosol/dust calculations
      !                 are done in the code. If set to FALSE, the code
      !                 will not read global aerosol/dust concentrations
      !                 and not do any aerosol/dust calculations.
      !=================================================================

      call rcEsmfReadLogical(gmiConfigFile, self%do_AerDust_Calc, "do_AerDust_Calc:", &
     &                       default=.false., rc=STATUS)

      !=================================================================
      ! AerDust_Effect_opt is used to select if the radiative effects
      !                    or/and heterogeneous chemistry on different
      !                    aerosols/dust are turned on/off.
      !     0: radiative effects on  and heterogeneous chemistry on
      !     1: radiative effects off and heterogeneous chemistry on
      !     2: radiative effects on  and heterogeneous chemistry off
      !     3: radiative effects off and heterogeneous chemistry off
      !=================================================================

      call ESMF_ConfigGetAttribute(gmiConfigFile, self%AerDust_Effect_opt, &
     &                label   = "AerDust_Effect_opt:", &
     &                default = 0, rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(gmiConfigFile, self%AerDust_infile_name, &
     &                label   = "AerDust_infile_name:", &
     &                default = '', rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(gmiConfigFile, self%Aerosol_infile_name, &
     &                label   = "Aerosol_infile_name:", &
     &                default = '', rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(gmiConfigFile, self%Dust_infile_name, &
     &                label   = "Dust_infile_name:", &
     &                default = '', rc=STATUS )
      VERIFY_(STATUS)

      call CheckNamelistOptionRange ('phot_opt', self%phot_opt, 0, 7)
      call CheckNamelistOptionRange ('fastj_opt', self%fastj_opt, 0, 4)
      call CheckNamelistOptionRange ('AerDust_Effect_opt', self%AerDust_Effect_opt, 0, 3)

      self%num_qjs       = NUM_J
      self%num_qjo       = NUM_J

      !!if (self%num_qjs > 0) self%qj_labels(1:self%num_qjs) = lqjchem(1:self%num_qjs)

      call rcEsmfReadLogical(gmiConfigFile, self%pr_qj_o3_o1d, "pr_qj_o3_o1d:", &
     &                       default=.false., rc=STATUS)

      call rcEsmfReadLogical(gmiConfigFile, self%pr_qj_opt_depth, "pr_qj_opt_depth:", &
     &                       default=.false., rc=STATUS)

!      if (self%pr_qj_o3_o1d) then
!         self%num_qjo = self%num_qjo + 1
!         self%qj_labels(self%num_qjo) = 'O3 + hv = O1D + O2'
!      end if
!
!      if (self%pr_qj_opt_depth) then
!         self%num_qjo = self%num_qjo + 1
!         self%qj_labels(self%num_qjo)   = 'optical depth'
!      end if

!     -----------------------------------------------------
!     chem_mask_klo, chem_mask_khi:
!       chemistry turned off where k is outside of range of
!       [chem_mask_klo, chem_mask_khi]
!     -----------------------------------------------------

      call ESMF_ConfigGetAttribute(gmiConfigFile, self%chem_mask_klo, &
     &                label   = "chem_mask_klo:", &
     &                default = 1, rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(gmiConfigFile, self%chem_mask_khi, &
     &                label   = "chem_mask_khi:", &
     &                default = km, rc=STATUS )
      VERIFY_(STATUS)

!     -------------------------------------------------------------------
!     synoz_threshold:  chemistry turned off where synoz > this threshold
!     -------------------------------------------------------------------

      hugeReal = Huge (hugeReal)
      
      call ESMF_ConfigGetAttribute(gmiConfigFile, self%synoz_threshold, &
     &                label   = "synoz_threshold:", &
     &                default = hugeReal, rc=STATUS )
      VERIFY_(STATUS)

      ! ----------------------------------------
      ! Do we want to couple to GOCART aerosols?
      ! ----------------------------------------
      
      call rcEsmfReadLogical(gmiConfigFile, self%usingGOCART_BC, &
     &           "usingGOCART_BC:", default=.false., rc=STATUS)
      VERIFY_(STATUS)
      
      call rcEsmfReadLogical(gmiConfigFile, self%usingGOCART_DU, &
     &           "usingGOCART_DU:", default=.false., rc=STATUS)
      VERIFY_(STATUS)
      
      call rcEsmfReadLogical(gmiConfigFile, self%usingGOCART_OC, &
     &           "usingGOCART_OC:", default=.false., rc=STATUS)
      VERIFY_(STATUS)
      
      call rcEsmfReadLogical(gmiConfigFile, self%usingGOCART_SS, &
     &           "usingGOCART_SS:", default=.false., rc=STATUS)
      VERIFY_(STATUS)
      
      call rcEsmfReadLogical(gmiConfigFile, self%usingGOCART_SU, &
     &           "usingGOCART_SU:", default=.false., rc=STATUS)
      VERIFY_(STATUS)


   IF( MAPL_AM_I_ROOT() ) THEN
    PRINT *," "
    PRINT *,TRIM(IAm)//":"
    PRINT *," Using GOCART   Black Carbon: ",self%usingGOCART_BC
    PRINT *," Using GOCART           Dust: ",self%usingGOCART_DU
    PRINT *," Using GOCART Organic Carbon: ",self%usingGOCART_OC
    PRINT *," Using GOCART       Sea Salt: ",self%usingGOCART_SS
    PRINT *," Using GOCART        Sulfate: ",self%usingGOCART_SU
    PRINT *," "
   END IF

! Perform consistency checks for aerosols
! ---------------------------------------
   IF(self%AM_I_AERO_PROVIDER) THEN
    IF(self%usingGOCART_BC .OR. self%usingGOCART_DU .OR. self%usingGOCART_OC .OR. &
       self%usingGOCART_SS .OR. self%usingGOCART_SU) THEN
     STATUS = 1
     IF( MAPL_AM_I_ROOT() ) THEN
      PRINT *," "
      PRINT *,TRIM(IAm)//":"
      PRINT *," AM_I_AERO_PROVIDER: ",self%AM_I_AERO_PROVIDER
      PRINT *," Cannot couple GOCART aerosols to GMICHEM when GMICHEM is the AERO_PROVIDER"
     END IF
     VERIFY_(STATUS)
    END IF
   END IF

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

!... solar cycle parameter
       if(self%do_solar_cycle) then

   CALL readSolarCycleData( s_cycle_dates, s_cycle, self%sc_infile_name )

!... figure out index for solar cycle array from year and month
          inyr = int(nymd/10000)
          inmon = int(nymd/100)-100*inyr
          iscyr = nint(((inyr+inmon/12.0)-s_cycle_dates(1))*12.)+1

    IF( MAPL_AM_I_ROOT() ) THEN
      PRINT *,"Solar cycle: ",s_cycle_dates(iscyr),s_cycle(:,iscyr)
    ENDIF
         self%JXbundle%fjx_solar_cycle_param(:) = s_cycle(:,iscyr)
       else
         self%JXbundle%fjx_solar_cycle_param(:) = 1.000
       endif

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

! Number of species and perform a consistency check with setkin_par.h.
! NOTES:
!  1. H2O is specie number 10 in the strat-trop mechanism, but will not be
!     found in w_c%reg%vname. H2O will be initialized from specific humidity, Q.
!  2. The GEOS-5 bundle has an Age-Of-Air tracer, which is not carried by GMI.
!  3. At the end of the XX (non-transported) species is a place holder for T2M15d.
! So w_c%reg%j_XX-w_c%reg%i_GMI must equal the parameter NSP = NCONST + NDYN.
! --------------------------------------------------------------------------------
   numSpecies = w_c%reg%j_XX-w_c%reg%i_GMI
   IF(numSpecies /= NSP) THEN
    PRINT *,TRIM(IAm),': Number of species from Chem_Registry.rc does not match number in setkin_par.h'
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
   IF(.NOT. found) THEN
    IF( MAPL_AM_I_ROOT() ) THEN
     PRINT *,TRIM(IAm),": Unable to find an OPEN logical unit to read qj_labels."
    END IF
    rc = 62
    RETURN
   END IF

   OPEN(UNIT=LogicalUnitNum,FILE=TRIM(kineticsTextFile),STATUS='old', &
        ACTION='read',FORM='formatted')
   DO
    READ(LogicalUnitNum,FMT="(A)") string
    IF(TRIM(string) == 'Photolyses:') EXIT
   END DO
   DO ic=1,NUM_J
    READ(LogicalUnitNum,FMT="(5X,A)") string
    self%qj_labels(ic)=TRIM(string)
   END DO

! Allocate space, etc., but the initialization of the
! species from the internal state is left to the run method.
! ----------------------------------------------------------

      CALL InitializeSpcConcentration(self%SpeciesConcentration,              &
     &               self%gmiGrid, gmiConfigFile, numSpecies, NMF, NCHEM,     &
     &               loc_proc)

      if (self%phot_opt /= 0) then
         Allocate(self%qjgmi(self%num_qjo))
         do ic = 1, self%num_qjo
            Allocate(self%qjgmi(ic)%pArray3D(i1:i2, ju1:j2, k1:k2))
            self%qjgmi(ic)%pArray3D(i1:i2, ju1:j2, k1:k2) = 0.0d0
         end do

         Allocate(self%overheadO3col(i1:i2, ju1:j2, k1:k2))
         self%overheadO3col = 0.0d0
      endif

      if (self%phot_opt == 2) then
         !Allocate (qjmon (i1:i2, ju1:j2, k1:k2, NUM_J, self%qj_timpyr))
         !qjmon    = 0.0d0
      elseif (self%phot_opt == 3) then
         !===========================
         select case (self%fastj_opt)
         !===========================

            !=======
            case (0)
            !=======
                call InitializeFastJ (self%cross_section_file, self%rate_file, &
     &                         self%T_O3_climatology_file, n_qj_O3_2OH,   &
     &                         NUM_J, self%chem_mask_khi, k2, k1)
            !=======
            case (1)
            !=======
               call InitializeFastJX (self%JXbundle, self%cross_section_file, &
                                      self%rate_file, &
     &                        self%T_O3_climatology_file, rootProc,            &
     &                        n_qj_O3_2OH, NUM_J, self%chem_mask_khi, k2, k1)
            !=======
            case (2)
            !=======
                call InitializeFastJX53b (k1, k2, self%chem_mask_khi,          &
     &                         NUM_J, self%cross_section_file,                 &
     &                         self%scattering_data_file, self%rate_file,      &
     &                         self%T_O3_climatology_file)
            !=======
            case (3)
            !=======
                call InitializeFastJX53c (k1, k2, self%chem_mask_khi,          &
     &                         NUM_J, self%cross_section_file,                 &
     &                         self%scattering_data_file, self%rate_file,      &
     &                         self%T_O3_climatology_file)
            !=======
            case (4)
            !=======
                call initializeFastJX65 (k1, k2, self%chem_mask_khi, NUM_J,    &
                                  self%cross_section_file,                     &
                               self%rate_file, self%T_O3_climatology_file, rootproc)
!                    call initializeFastJX65 (k1, k2, self%chem_mask_khi, NUM_J,    &
!                               self%aerosolOpticalData_file,                   &
!                               self%cross_section_file,                        &
!                               self%scattering_data_file,                      &
!                               self%rate_file, self%T_O3_climatology_file)
       !=========
         end select
         !=========
      end if

      if ((TRIM(self%chem_mecha) ==        'troposphere') .OR. &
     &    (TRIM(self%chem_mecha) ==         'strat_trop') .OR. &
     &    (TRIM(self%chem_mecha) == 'strat_trop_aerosol')) THEN
         if ((self%phot_opt == 3) .and. self%do_AerDust_Calc) then
            Allocate(self%dust(i1:i2, ju1:j2, k1:k2, nSADdust))
            self%dust = 0.0d0

            Allocate(self%dAersl(i1:i2, ju1:j2, k1:k2, 2))
            self%dAersl = 0.0d0

            Allocate(self%wAersl(i1:i2, ju1:j2, k1:k2, nSADaer))
            self%wAersl = 0.0d0

            Allocate(self%odAer(i1:i2, ju1:j2, k1:k2, nSADaer*NRH_b))
            self%odAer = 0.0d0

            Allocate(self%odMdust(i1:i2, ju1:j2, k1:k2, nSADdust))
            self%odMdust = 0.0d0

            Allocate(self%eRadius(i1:i2, ju1:j2, k1:k2, nSADdust+nSADaer))
            self%eRadius = 0.0d0

            Allocate(self%tArea(i1:i2, ju1:j2, k1:k2, nSADdust+nSADaer))
            self%tArea = 0.0d0

            Allocate(self%optDepth(i1:i2, ju1:j2, k1:k2, num_AerDust))
            self%optDepth = 0.0d0

            if(rootProc) PRINT *,"Chemistry initialize: allocated self%optDepth"
         end if
      end if

      if ((self%sfalbedo_opt == 1) .or. (self%sfalbedo_opt == 2)) then
         allocate(self%sasdir_data(i1:i2, ju1:j2, MONTHS_PER_YEAR))
         self%sasdir_data = 0.0d0

         allocate(self%sasdif_data(i1:i2, ju1:j2, MONTHS_PER_YEAR))
         self%sasdif_data = 0.0d0

         allocate(self%saldir_data(i1:i2, ju1:j2, MONTHS_PER_YEAR))
         self%saldir_data = 0.0d0

         allocate(self%saldif_data(i1:i2, ju1:j2, MONTHS_PER_YEAR))
         self%saldif_data = 0.0d0

         if (self%sfalbedo_opt  == 1) then
            self%sasdir_data(:,:,:) = self%sasdir_init_val
            self%sasdif_data(:,:,:) = self%sasdif_init_val
            self%saldir_data(:,:,:) = self%saldir_init_val
            self%saldif_data(:,:,:) = self%saldif_init_val
         end if
      end if


      !=======================
      ! Initialize the bundles
      !=======================

   ! Get the declared bundle from the state

   call ESMF_StateGet(expChem, 'gmiQJ' , qjBundle,   RC=STATUS)
   VERIFY_(STATUS)

   ! Add tracer to the bundle
   do ib = 1, self%num_qjo
      allocate( var(i1:i2, j1:j2, 1:km), STAT=STATUS)
      VERIFY_(STATUS)
      var(:,:,:)  = 0.0d0

      call addTracerToBundle (qjBundle, var, w_c%grid_esmf, self%qj_labels(ib))
   end do

   ! Sanity check

   call ESMF_FieldBundleGet(qjBundle, fieldCount=numVars , rc=STATUS)
   VERIFY_(STATUS)
   _ASSERT(self%num_qjo == numVars,'needs informative message')

   ! eRadius Bundle

   call ESMF_StateGet(expChem, 'gmiERADIUS' , eRadiusBundle,   RC=STATUS)
   VERIFY_(STATUS)

   ! Add tracer to the bundle
   do ib = 1, nSADdust+nSADaer
      allocate( var(i1:i2, j1:j2, 1:km), STAT=STATUS)
      VERIFY_(STATUS)
      var(:,:,:)  = 0.0d0

      write (binName ,'(i4.4)') ib
      varName = 'eRadius'//binName

      call addTracerToBundle (eRadiusBundle, var, w_c%grid_esmf, varName)
   end do

   ! Sanity check

   call ESMF_FieldBundleGet(eRadiusBundle, fieldCount=numVars , rc=STATUS)
   VERIFY_(STATUS)
   _ASSERT(nSADdust+nSADaer == numVars,'needs informative message')

   ! for tArea Bundle

   call ESMF_StateGet(expChem, 'gmiTAREA' , tAreaBundle,   RC=STATUS)
   VERIFY_(STATUS)

   ! Add tracer to the bundle
   do ib = 1, nSADdust+nSADaer
      allocate( var(i1:i2, j1:j2, 1:km), STAT=STATUS)
      VERIFY_(STATUS)
      var(:,:,:)  = 0.0d0

      write (binName ,'(i4.4)') ib
      varName = 'tArea'//binName

      call addTracerToBundle (tAreaBundle, var, w_c%grid_esmf, varName)
   end do

   ! Sanity check

   call ESMF_FieldBundleGet(tAreaBundle, fieldCount=numVars , rc=STATUS)
   VERIFY_(STATUS)
   _ASSERT(nSADdust+nSADaer == numVars,'needs informative message')

      !------------------------------------
      ! Initial Settings for Ship Emissions
      !------------------------------------

      self%jno2num = badIndex

      if (self%do_ShipEmission) then
         do ic=1, NUM_J
            if (TRIM(self%chem_mecha) ==         'strat_trop' .OR. &
	        TRIM(self%chem_mecha) == 'strat_trop_aerosol') THEN
               if (Trim(self%qj_labels(ic)) =='NO2 + hv = NO + O') self%jno2num = ic
            else if (TRIM(self%chem_mecha) == 'troposphere') then
               if (Trim(self%qj_labels(ic)) =='NO2 + hv = NO + O3') self%jno2num = ic
            endif
         end do
         if (self%jno2num .eq. badIndex) then
            print*,'jno2num not found in GmiPhotolysis_GridCompInitialize'
            stop
         endif

         ALLOCATE (jNO2val_phot(i1:i2,j1:j2))
         jNO2val_phot(i1:i2,j1:j2) = 0.0
         CALL initDataInStateField(expChem, w_c%grid_esmf, jNO2val_phot,  'jNO2val')
      endif

    !---------------------------------------------------------------
    ! Create and populate the array that maps GMI species indices to
    ! GEOS-5 species indices
    !---------------------------------------------------------------

    allocate(self%mapSpecies(NSP))
    self%mapSpecies(:) = speciesReg_for_CCM(lchemvar, w_c%reg%vname, NSP, w_c%reg%i_GMI, w_c%reg%j_XX)

  RETURN
   
  END SUBROUTINE GmiPhotolysis_GridCompInitialize

!-------------------------------------------------------------------------
!NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1, GEOS/DAS!
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  GmiPhotolysis_GridCompRun --- The GMI Driver
!
! !INTERFACE:
!

   SUBROUTINE GmiPhotolysis_GridCompRun ( self, w_c, impChem, expChem, nymd, nhms, &
                                tdt, rc )

! !USES:

   USE GmiTimeControl_mod,            ONLY : Set_curGmiDate, Set_curGmiTime
   USE GmiTimeControl_mod,            ONLY : Set_numTimeSteps, Get_numTimeSteps
   USE GmiTimeControl_mod,            ONLY : Set_gmiSeconds, GetSecondsFromJanuary1
   USE GmiTimeControl_mod,            ONLY : GetDaysFromJanuary1, ConvertTimeToSeconds
   USE GmiSpcConcentrationMethod_mod, ONLY : resetFixedConcentration
   USE GmiSolar_mod,		      ONLY : computeSolarZenithAngle_Photolysis
   USE GmiPhotRateConst_mod,	      ONLY : calcPhotolysisRateConstants

   IMPLICIT none

! !INPUT/OUTPUT PARAMETERS:

   TYPE(GmiPhotolysis_GridComp), INTENT(INOUT) :: self ! Grid Component
   TYPE(Chem_Bundle), INTENT(INOUT) :: w_c    ! Chemical tracer fields   

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

   CHARACTER(LEN=*), PARAMETER :: IAm = 'GmiPhotolysis_GridCompRun'

!  Input fields from GEOS-5
!  ------------------------
   REAL, POINTER, DIMENSION(:,:) :: cldtt, albvf

   REAL, POINTER, DIMENSION(:,:,:) :: airdens, ple, Q, T, zle
   REAL, POINTER, DIMENSION(:,:,:) :: fcld,taucli,tauclw,ql,cnv_mfc
   REAL, POINTER, DIMENSION(:,:,:) :: rh2,dqdt

!  Dust and aerosols.  May serve as imports from GOCART
!  or as exports to fill the AERO_BUNDLE, but not both.
!  ----------------------------------------------------
   REAL, POINTER, DIMENSION(:,:,:) :: BCphobic,BCphilic
   REAL, POINTER, DIMENSION(:,:,:) :: DU001,DU002,DU003,DU004,DU005
   REAL, POINTER, DIMENSION(:,:,:) :: OCphobic,OCphilic
   REAL, POINTER, DIMENSION(:,:,:) :: SS001,SS002,SS003,SS004,SS005
   REAL, POINTER, DIMENSION(:,:,:) :: DMS,SO2,SO4,MSA

!  Exports not part of internal state
!  ----------------------------------
   REAL, POINTER, DIMENSION(:,:) :: SZAPHOT
   REAL, POINTER, DIMENSION(:,:,:) :: DUSTOD, DUSTSA
   REAL, POINTER, DIMENSION(:,:,:) :: SO4OD, SO4HYGRO, SO4SA
   REAL, POINTER, DIMENSION(:,:,:) :: BCOD, BCHYGRO, BCSA
   REAL, POINTER, DIMENSION(:,:,:) :: OCOD, OCHYGRO, OCSA
   REAL, POINTER, DIMENSION(:,:,:) :: SSAOD, SSAHYGRO, SSASA, SSCOD, SSCHYGRO, SSCSA

!  Export for Ship Emissions
   REAL, POINTER, DIMENSION(:,:) :: jNO2val_phot

!  Local
!  -----
   INTEGER :: cymd, dymd, emiss_opt, hms
   INTEGER :: i, i1, i2, ic, idehyd_num, im, iXj, iTile(1)
   INTEGER :: j, j1, j2, jm, jTile(1)
   INTEGER :: k, km, kReverse
   INTEGER :: ilo, gmi_nborder
   INTEGER :: ihi, julo, jhi, ju1,  k1, k2, ilong, ilat
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

   REAL(KIND=DBL) :: dehydmin = 0.00
   REAL(KIND=DBL) :: chemDt, dayOfYear

   CHARACTER(LEN=255) :: speciesName
   CHARACTER(LEN=255) :: importName

   LOGICAL :: found, rootProc
   LOGICAL, PARAMETER :: doThis = .FALSE.
   
   integer :: nsec_jan1, jday
   real*8  :: rsec_jan1, time_sec

! Allocatables.  Use KIND=DBL where GMI expects REAL*8.
! -----------------------------------------------------
   REAL, ALLOCATABLE :: var2d(:,:)
   REAL, ALLOCATABLE :: var3d(:,:,:)
   REAL, ALLOCATABLE :: pl(:,:,:)

   REAL(KIND=DBL), ALLOCATABLE :: lonDeg(:,:)
   REAL(KIND=DBL), ALLOCATABLE :: latDeg(:,:)

   REAL(KIND=DBL), ALLOCATABLE :: averagePressEdge(:)

   REAL(KIND=DBL), ALLOCATABLE :: tropopausePress(:,:)
   REAL(KIND=DBL), ALLOCATABLE :: pctm2(:,:)
   REAL(KIND=DBL), ALLOCATABLE :: fracCloudCover(:,:)
   REAL(KIND=DBL), ALLOCATABLE :: surf_alb(:,:)

   REAL(KIND=DBL), ALLOCATABLE :: var4dDBL(:,:,:,:)

   REAL(KIND=DBL), ALLOCATABLE :: mass(:,:,:)
   REAL(KIND=DBL), ALLOCATABLE :: press3c(:,:,:)
   REAL(KIND=DBL), ALLOCATABLE :: press3e(:,:,:)
   REAL(KIND=DBL), ALLOCATABLE :: gridBoxThickness(:,:,:)
   REAL(KIND=DBL), ALLOCATABLE :: kel(:,:,:)
   REAL(KIND=DBL), ALLOCATABLE :: humidity(:,:,:)
   REAL(KIND=DBL), ALLOCATABLE :: totalCloudFraction(:,:,:)
   REAL(KIND=DBL), ALLOCATABLE :: tau_cloud(:,:,:)
   REAL(KIND=DBL), ALLOCATABLE :: tau_clw(:,:,:)
   REAL(KIND=DBL), ALLOCATABLE :: tau_cli(:,:,:)
   REAL(KIND=DBL), ALLOCATABLE :: clwc(:,:,:)
   REAL(KIND=DBL), ALLOCATABLE :: cmf(:,:,:)
   REAL(KIND=DBL), ALLOCATABLE :: relativeHumidity(:,:,:)
   REAL(KIND=DBL), ALLOCATABLE :: moistq(:,:,:)

   REAL(KIND=DBL), ALLOCATABLE :: solarZenithAngle(:,:)

   REAL*8  :: tdt8

   tdt8 = tdt

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
   ju1    = j1
   k1     = 1
   k2     = km
   ilo    = i1  - gmi_nborder
   ihi    = i2  + gmi_nborder
   julo   = ju1 - gmi_nborder
   jhi    = j2  + gmi_nborder
   ilong  = i2 - i1  + 1
   ilat   = j2 - ju1 + 1

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
   ALLOCATE(averagePressEdge(0:km),STAT=STATUS)
   VERIFY_(STATUS)

   ALLOCATE(              var2d(i1:i2,j1:j2),STAT=STATUS)
   VERIFY_(STATUS)
   ALLOCATE(    tropopausePress(i1:i2,j1:j2),STAT=STATUS)
   VERIFY_(STATUS)
   ALLOCATE(              pctm2(i1:i2,j1:j2),STAT=STATUS)
   VERIFY_(STATUS)
   ALLOCATE(     fracCloudCover(i1:i2,j1:j2),STAT=STATUS)
   VERIFY_(STATUS)
   ALLOCATE(           surf_alb(i1:i2,j1:j2),STAT=STATUS)
   VERIFY_(STATUS)
   ALLOCATE(   solarZenithAngle(i1:i2,j1:j2),STAT=STATUS)
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
   ALLOCATE(totalCloudFraction(i1:i2,j1:j2,1:km),STAT=STATUS)
   VERIFY_(STATUS)
   ALLOCATE(         tau_cloud(i1:i2,j1:j2,1:km),STAT=STATUS)
   VERIFY_(STATUS)
   ALLOCATE(              clwc(i1:i2,j1:j2,1:km),STAT=STATUS)
   VERIFY_(STATUS)
   ALLOCATE(               cmf(i1:i2,j1:j2,1:km),STAT=STATUS)
   VERIFY_(STATUS)
   ALLOCATE(  relativeHumidity(i1:i2,j1:j2,1:km),STAT=STATUS)
   VERIFY_(STATUS)
   ALLOCATE(            moistq(i1:i2,j1:j2,1:km),STAT=STATUS)
   VERIFY_(STATUS)
   ALLOCATE(           tau_clw(i1:i2,j1:j2,1:km),STAT=STATUS)
   VERIFY_(STATUS)
   ALLOCATE(           tau_cli(i1:i2,j1:j2,1:km),STAT=STATUS)
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

! Obtain dust and aerosols from either GOCART or
! from climatologies supplied by the GMI project
! ----------------------------------------------
   CALL Acquire_BC(STATUS)
   VERIFY_(STATUS)
   CALL Acquire_DU(STATUS)
   VERIFY_(STATUS)
   CALL Acquire_OC(STATUS)
   VERIFY_(STATUS)
   CALL Acquire_SS(STATUS)
   VERIFY_(STATUS)
   CALL Acquire_SU(STATUS)
   VERIFY_(STATUS)

! Grab imports and do units conversions
! -------------------------------------
   CALL SatisfyImports(STATUS)
   VERIFY_(STATUS)

! Hand the species concentrations to GMI's bundle
! -----------------------------------------------
   IF (self%gotImportRst) then
      CALL SwapSpeciesBundles(ToGMI, self%SpeciesConcentration%concentration, &
               w_c%qa, Q, self%mapSpecies, lchemvar, self%do_synoz, NSP, &
               STATUS)
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
! GMI Photolysis
! ------------------------------------------------------------------------

      IF (self%gotImportRst) THEN
         if ((self%phot_opt == 3) .or. (self%phot_opt == 4) .or.  &
             (self%phot_opt == 5)) then 

            call GetSecondsFromJanuary1 (nsec_jan1, nymd, nhms)
            rsec_jan1 = nsec_jan1

            call GetDaysFromJanuary1 (jday, nymd)
            time_sec = ConvertTimeToSeconds (nhms)

            solarZenithAngle(i1:i2,j1:j2) = &
                 computeSolarZenithAngle_Photolysis (jday, time_sec, &
                          self%fastj_offset_sec, latDeg, lonDeg, i1, i2, j1, j2)

            call calcPhotolysisRateConstants (self%JXbundle,                   &
                     self%chem_mecha, tropopausePress,                         &
     &               self%pr_qj_o3_o1d, self%pr_qj_opt_depth, tdt8,            &
     &               rsec_jan1, pctm2, mass, press3e, press3c, kel,            &
     &               self%SpeciesConcentration%concentration, solarZenithAngle,&
     &               self%cellArea, surf_alb, fracCloudCover,                  &
     &               tau_cloud, tau_clw, tau_cli,                              &
     &               self%overheadO3col, self%qjgmi, gridBoxThickness,         &
     &               self%optDepth, self%eRadius, self%tArea, self%odAer,      &
     &               relativeHumidity, self%odMdust, self%dust, self%wAersl,   &
     &               self%dAersl, humidity, num_AerDust, self%phot_opt,        &
     &               self%fastj_opt, self%fastj_offset_sec, self%do_clear_sky, &
     &               self%do_AerDust_Calc, self%do_ozone_inFastJX,             &
     &               self%do_synoz, self%qj_timpyr, IO3, IH2O, ISYNOZ,         &
     &               self%chem_mask_khi, nymd, nhms, self%pr_diag, loc_proc,   &
     &               self%synoz_threshold, self%AerDust_Effect_opt, NSP,       &
     &               NUM_J, self%num_qjo, ilo, ihi, julo, jhi, &
     &               i1, i2, ju1, j2, k1, k2, self%jNOindex, self%jNOamp)
         endif
      END IF

! Return species concentrations to the chemistry bundle
! -----------------------------------------------------
   IF (self%gotImportRst) then
      CALL SwapSpeciesBundles(FromGMI, self%SpeciesConcentration%concentration, &
               w_c%qa, Q, self%mapSpecies, lchemvar, self%do_synoz, NSP,  &
               STATUS)
      VERIFY_(STATUS)
   END IF

! Export states
! -------------

   CALL FillExports(STATUS)
   VERIFY_(STATUS)

   CALL populateBundleQJ ( )

! Scratch local work space
! ------------------------
   DEALLOCATE(lonDeg, latDeg, averagePressEdge, STAT=STATUS)
   VERIFY_(STATUS)

   DEALLOCATE(var2d, tropopausePress, pctm2, fracCloudCover, surf_alb, &
              solarZenithAngle, STAT=STATUS)
   VERIFY_(STATUS)

   DEALLOCATE(pl, mass, press3c, press3e, gridBoxThickness, kel, humidity, &
              totalCloudFraction, tau_cloud, clwc, cmf, relativeHumidity, &
              moistq, STAT=STATUS)
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
! !ROUTINE:  Acquire_BC
!
! !INTERFACE:

  SUBROUTINE Acquire_BC(rc)
  
  IMPLICIT NONE

  INTEGER, INTENT(OUT) :: rc
!
! !DESCRIPTION:
!
!  Obtain black carbon either by import from GOCART or from an external
!  file specified in ExtData.rc along with a refresh template.
!
!  Hydrophobic black carbon: self%dAersl(:,:,:,1)
!  Hydrophilic black carbon: self%wAersl(:,:,:,2)
!
!EOP
!---------------------------------------------------------------------------

  CHARACTER(LEN=255) :: IAm
  REAL, POINTER, DIMENSION(:,:,:) :: PTR3D
  REAL :: qmin,qmax
  
  TYPE(ESMF_State)           :: aero
  TYPE(ESMF_FieldBundle)     :: aerosols

  rc = 0
  IAm = "Acquire_BC"

  SELECT CASE (TRIM(self%aeroProviderName))

   CASE("GOCART.data")

    CALL ESMF_StateGet(impChem, 'AERO', aero, RC=STATUS)
    VERIFY_(STATUS)
    CALL ESMF_StateGet(aero, 'AEROSOLS', aerosols, RC=STATUS)
    VERIFY_(STATUS)

    CALL ESMFL_BundleGetPointertoData(aerosols, 'BCphobic', PTR3D, RC=STATUS)
    VERIFY_(STATUS)
    self%dAersl(:,:,km:1:-1,1) = PTR3D(:,:,1:km)*airdens(:,:,1:km)
    IF(self%verbose) CALL pmaxmin('BCphobic:', PTR3D, qmin, qmax, iXj, km, 1. )
    NULLIFY(PTR3D)

    CALL ESMFL_BundleGetPointertoData(aerosols, 'BCphilic', PTR3D, RC=STATUS)
    VERIFY_(STATUS)
    self%wAersl(:,:,km:1:-1,2) = PTR3D(:,:,1:km)*airdens(:,:,1:km)
    IF(self%verbose) CALL pmaxmin('BCphilic:', PTR3D, qmin, qmax, iXj, km, 1. )
    NULLIFY(PTR3D)

   CASE("GOCART")

    IF(self%usingGOCART_BC) THEN

     CALL MAPL_GetPointer(impChem, BCphobic, 'GOCART::BCphobic', RC=STATUS)
     VERIFY_(STATUS)
     CALL MAPL_GetPointer(impChem, BCphilic, 'GOCART::BCphilic', RC=STATUS)
     VERIFY_(STATUS)

     self%dAersl(:,:,km:1:-1,1) = BCphobic(:,:,1:km)*airdens(:,:,1:km)
     self%wAersl(:,:,km:1:-1,2) = BCphilic(:,:,1:km)*airdens(:,:,1:km)

     IF(self%verbose) THEN
      CALL pmaxmin('BCphobic:', BCphobic, qmin, qmax, iXj, km, 1. )
      CALL pmaxmin('BCphilic:', BCphilic, qmin, qmax, iXj, km, 1. )
     END IF

    END IF

   CASE("GMICHEM")

    CALL MAPL_GetPointer(impChem, PTR3D, 'BC1', RC=STATUS)
    VERIFY_(STATUS)
    self%dAersl(:,:,1:km,1) = PTR3D(:,:,km:1:-1)
    IF(self%verbose) CALL pmaxmin('BCphobic:', PTR3D, qmin, qmax, iXj, km, 1. )
    IF(ASSOCIATED(BCphobic) .AND. self%AM_I_AERO_PROVIDER) BCphobic(:,:,:) = PTR3D(:,:,:)/airdens(:,:,:)
    NULLIFY(PTR3D)

    CALL MAPL_GetPointer(impChem, PTR3D, 'BC2', RC=STATUS)
    VERIFY_(STATUS)
    self%wAersl(:,:,1:km,2) = PTR3D(:,:,km:1:-1)
    IF(self%verbose) CALL pmaxmin('BCphilic:', PTR3D, qmin, qmax, iXj, km, 1. )
    IF(ASSOCIATED(BCphilic) .AND. self%AM_I_AERO_PROVIDER) BCphilic(:,:,:) = PTR3D(:,:,:)/airdens(:,:,:)
    NULLIFY(PTR3D)

   CASE DEFAULT

    STATUS = 1
    VERIFY_(STATUS)

  END SELECT

  RETURN
 END SUBROUTINE Acquire_BC

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
!
!EOP
!---------------------------------------------------------------------------

  CHARACTER(LEN=255) :: IAm
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
! !ROUTINE:  Acquire_DU
!
! !INTERFACE:

  SUBROUTINE Acquire_DU(rc)
  
  IMPLICIT NONE

  INTEGER, INTENT(OUT) :: rc
!
! !DESCRIPTION:
!
!  Obtain dust either by import from GOCART or from an external
!  file specified in ExtData.rc along with a refresh template.
!  This version is expecting seven dust bins.  In the GOCART case,
!  the partitioning follows the recommendation of Pete Colarco.
!
!EOP
!---------------------------------------------------------------------------

  CHARACTER(LEN=255) :: IAm
  REAL, POINTER, DIMENSION(:,:,:) :: PTR3D
  REAL :: qmin,qmax
  
  TYPE(ESMF_State)           :: aero
  TYPE(ESMF_FieldBundle)     :: aerosols

  rc = 0
  IAm = "Acquire_DU"

  SELECT CASE (TRIM(self%aeroProviderName))

   CASE("GOCART.data")

    CALL ESMF_StateGet(impChem, 'AERO', aero, RC=STATUS)
    VERIFY_(STATUS)
    CALL ESMF_StateGet(aero, 'AEROSOLS', aerosols, RC=STATUS)
    VERIFY_(STATUS)

    CALL ESMFL_BundleGetPointertoData(aerosols, 'du001', PTR3D, RC=STATUS)
    VERIFY_(STATUS)
    self%dust(:,:,km:1:-1,1) = PTR3D(:,:,1:km)*airdens(:,:,1:km)*0.009
    self%dust(:,:,km:1:-1,2) = PTR3D(:,:,1:km)*airdens(:,:,1:km)*0.081
    self%dust(:,:,km:1:-1,3) = PTR3D(:,:,1:km)*airdens(:,:,1:km)*0.234
    self%dust(:,:,km:1:-1,4) = PTR3D(:,:,1:km)*airdens(:,:,1:km)*0.676
    IF(self%verbose) CALL pmaxmin('du001:', PTR3D, qmin, qmax, iXj, km, 1. )
    NULLIFY(PTR3D)

    CALL ESMFL_BundleGetPointertoData(aerosols, 'du002', PTR3D, RC=STATUS)
    VERIFY_(STATUS)
    self%dust(:,:,km:1:-1,5) = PTR3D(:,:,1:km)*airdens(:,:,1:km)
    IF(self%verbose) CALL pmaxmin('du002:', PTR3D, qmin, qmax, iXj, km, 1. )
    NULLIFY(PTR3D)

    CALL ESMFL_BundleGetPointertoData(aerosols, 'du003', PTR3D, RC=STATUS)
    VERIFY_(STATUS)
    self%dust(:,:,km:1:-1,6) = PTR3D(:,:,1:km)*airdens(:,:,1:km)
    IF(self%verbose) CALL pmaxmin('du003:', PTR3D, qmin, qmax, iXj, km, 1. )
    NULLIFY(PTR3D)

    CALL ESMFL_BundleGetPointertoData(aerosols, 'du004', PTR3D, RC=STATUS)
    VERIFY_(STATUS)
    self%dust(:,:,km:1:-1,7) = PTR3D(:,:,1:km)*airdens(:,:,1:km)
    IF(self%verbose) CALL pmaxmin('du004:', PTR3D, qmin, qmax, iXj, km, 1. )
    NULLIFY(PTR3D)


   CASE("GOCART")

    IF(self%usingGOCART_DU) THEN

     CALL MAPL_GetPointer(impChem, DU001, 'GOCART::du001', RC=STATUS)
     VERIFY_(STATUS)
     CALL MAPL_GetPointer(impChem, DU002, 'GOCART::du002', RC=STATUS)
     VERIFY_(STATUS)
     CALL MAPL_GetPointer(impChem, DU003, 'GOCART::du003', RC=STATUS)
     VERIFY_(STATUS)
     CALL MAPL_GetPointer(impChem, DU004, 'GOCART::du004', RC=STATUS)
     VERIFY_(STATUS)

     IF(self%verbose) THEN
      CALL pmaxmin('DU001:', DU001, qmin, qmax, iXj, km, 1. )
      CALL pmaxmin('DU002:', DU002, qmin, qmax, iXj, km, 1. )
      CALL pmaxmin('DU003:', DU003, qmin, qmax, iXj, km, 1. )
      CALL pmaxmin('DU004:', DU004, qmin, qmax, iXj, km, 1. )
     END IF

     self%dust(:,:,km:1:-1,1) = DU001(:,:,1:km)*airdens(:,:,1:km)*0.009
     self%dust(:,:,km:1:-1,2) = DU001(:,:,1:km)*airdens(:,:,1:km)*0.081
     self%dust(:,:,km:1:-1,3) = DU001(:,:,1:km)*airdens(:,:,1:km)*0.234
     self%dust(:,:,km:1:-1,4) = DU001(:,:,1:km)*airdens(:,:,1:km)*0.676
     self%dust(:,:,km:1:-1,5) = DU002(:,:,1:km)*airdens(:,:,1:km)
     self%dust(:,:,km:1:-1,6) = DU003(:,:,1:km)*airdens(:,:,1:km)
     self%dust(:,:,km:1:-1,7) = DU004(:,:,1:km)*airdens(:,:,1:km)

    END IF

   CASE("GMICHEM")

    DO i = 1,nSADdust

     WRITE(importName, FMT="('MDUST',I1)") i
     CALL MAPL_GetPointer(impChem, PTR3D, TRIM(importName), RC=STATUS)
     VERIFY_(STATUS)
     self%dust(:,:,1:km,i) = PTR3D(:,:,km:1:-1)
     IF(self%verbose) CALL pmaxmin(TRIM(importName)//":", PTR3D, qmin, qmax, iXj, km, 1. )

     SELECT CASE (i)
      CASE(5)
       IF(ASSOCIATED(DU002) .AND. self%AM_I_AERO_PROVIDER) DU002(:,:,:) = PTR3D(:,:,:)/airdens(:,:,:)
      CASE(6)
       IF(ASSOCIATED(DU003) .AND. self%AM_I_AERO_PROVIDER) DU003(:,:,:) = PTR3D(:,:,:)/airdens(:,:,:)
      CASE(7)
       IF(ASSOCIATED(DU004) .AND. self%AM_I_AERO_PROVIDER) DU004(:,:,:) = PTR3D(:,:,:)/airdens(:,:,:)
      CASE DEFAULT
       IF(ASSOCIATED(DU001) .AND. self%AM_I_AERO_PROVIDER) THEN
        IF(i == 1) DU001(:,:,:) = 0.00
        DU001(:,:,:) = DU001(:,:,:) + PTR3D(:,:,:)/airdens(:,:,:)
       END IF
     END SELECT

     NULLIFY(PTR3D)

    END DO

   CASE DEFAULT

    STATUS = 1
    VERIFY_(STATUS)

  END SELECT

  RETURN
 END SUBROUTINE Acquire_DU

!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  Acquire_OC
!
! !INTERFACE:

  SUBROUTINE Acquire_OC(rc)
  
  IMPLICIT NONE

  INTEGER, INTENT(OUT) :: rc
!
! !DESCRIPTION:
!
!  Obtain organic carbon either by import from GOCART or from an external
!  file specified in ExtData.rc along with a refresh template.
!
!  Hydrophobic organic carbon: OC001, self%dAersl(:,:,:,2)
!  Hydrophilic organic carbon: OC002, self%wAersl(:,:,:,3)
!
!EOP
!---------------------------------------------------------------------------

  CHARACTER(LEN=255) :: IAm
  REAL, POINTER, DIMENSION(:,:,:) :: PTR3D
  REAL :: qmin,qmax
  
  TYPE(ESMF_State)           :: aero
  TYPE(ESMF_FieldBundle)     :: aerosols

  rc = 0
  IAm = "Acquire_OC"

  SELECT CASE (TRIM(self%aeroProviderName))

   CASE("GOCART.data")

    CALL ESMF_StateGet(impChem, 'AERO', aero, RC=STATUS)
    VERIFY_(STATUS)
    CALL ESMF_StateGet(aero, 'AEROSOLS', aerosols, RC=STATUS)
    VERIFY_(STATUS)

    CALL ESMFL_BundleGetPointertoData(aerosols, 'OCphobic', PTR3D, RC=STATUS)
    VERIFY_(STATUS)
    self%dAersl(:,:,km:1:-1,2) = PTR3D(:,:,1:km)*airdens(:,:,1:km)
    IF(self%verbose) CALL pmaxmin('OCphobic:', PTR3D, qmin, qmax, iXj, km, 1. )
    NULLIFY(PTR3D)

    CALL ESMFL_BundleGetPointertoData(aerosols, 'OCphilic', PTR3D, RC=STATUS)
    VERIFY_(STATUS)
    self%wAersl(:,:,km:1:-1,3) = PTR3D(:,:,1:km)*airdens(:,:,1:km)
    IF(self%verbose) CALL pmaxmin('OCphilic:', PTR3D, qmin, qmax, iXj, km, 1. )
    NULLIFY(PTR3D)

   CASE("GOCART")

    IF(self%usingGOCART_OC) THEN

     CALL MAPL_GetPointer(impChem, OCphobic, 'GOCART::OCphobic', RC=STATUS)
     VERIFY_(STATUS)
     CALL MAPL_GetPointer(impChem, OCphilic, 'GOCART::OCphilic', RC=STATUS)
     VERIFY_(STATUS)

     self%dAersl(:,:,km:1:-1,2) = OCphobic(:,:,1:km)*airdens(:,:,1:km)
     self%wAersl(:,:,km:1:-1,3) = OCphilic(:,:,1:km)*airdens(:,:,1:km)

     IF(self%verbose) THEN
      CALL pmaxmin('OCphobic:', OCphobic, qmin, qmax, iXj, km, 1. )
      CALL pmaxmin('OCphilic:', OCphilic, qmin, qmax, iXj, km, 1. )
     END IF

    END IF

   CASE("GMICHEM")

    CALL MAPL_GetPointer(impChem, PTR3D, 'OC1', RC=STATUS)
    VERIFY_(STATUS)
    self%dAersl(:,:,1:km,2) = PTR3D(:,:,km:1:-1)
    IF(self%verbose) CALL pmaxmin('OCphobic:', PTR3D, qmin, qmax, iXj, km, 1. )
    IF(ASSOCIATED(OCphobic) .AND. self%AM_I_AERO_PROVIDER) OCphobic(:,:,:) = PTR3D(:,:,:)/airdens(:,:,:)
    NULLIFY(PTR3D)

    CALL MAPL_GetPointer(impChem, PTR3D, 'OC2', RC=STATUS)
    VERIFY_(STATUS)
    self%wAersl(:,:,1:km,3) = PTR3D(:,:,km:1:-1)
    IF(self%verbose) CALL pmaxmin('OCphilic:', PTR3D, qmin, qmax, iXj, km, 1. )
    IF(ASSOCIATED(OCphilic) .AND. self%AM_I_AERO_PROVIDER) OCphilic(:,:,:) = PTR3D(:,:,:)/airdens(:,:,:)
    NULLIFY(PTR3D)

   CASE DEFAULT

    STATUS = 1
    VERIFY_(STATUS)

  END SELECT

  RETURN
 END SUBROUTINE Acquire_OC

!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  Acquire_SS
!
! !INTERFACE:

  SUBROUTINE Acquire_SS(rc)
  
  IMPLICIT NONE

  INTEGER, INTENT(OUT) :: rc
!
! !DESCRIPTION:
!
!  Obtain sea salts either by import from GOCART or from an external
!  file specified in ExtData.rc along with a refresh template.
!
!  Sea salts: self%wAersl(:,:,:,4:5)
!
!  The following note is from Huisheng Bian describing the SS bin mapping
!  for both the monthly mean data and for GOCART imports:
!
!  Mian's GOCART uses 4 bins for SS simulation, while Pete's GOCART uses 5 
!  bins for SS simulation. Since Bryan obtained the data from Thomas, it 
!  should have 4 SS bins. In this case, the fine mode SS is only the bin 1 
!  and the rest bins belong to coarse mode. For the coupled model, the fine 
!  mode SS is the sum of bin 1 and 2.
!
!EOP
!---------------------------------------------------------------------------

  CHARACTER(LEN=255) :: IAm
  REAL, POINTER, DIMENSION(:,:,:) :: PTR3D
  REAL :: qmin,qmax
  
  TYPE(ESMF_State)           :: aero
  TYPE(ESMF_FieldBundle)     :: aerosols

  rc = 0
  IAm = "Acquire_SS"

  SELECT CASE (TRIM(self%aeroProviderName))

   CASE("GOCART.data")

    CALL ESMF_StateGet(impChem, 'AERO', aero, RC=STATUS)
    VERIFY_(STATUS)
    CALL ESMF_StateGet(aero, 'AEROSOLS', aerosols, RC=STATUS)
    VERIFY_(STATUS)

    CALL ESMFL_BundleGetPointertoData(aerosols, 'ss001', PTR3D, RC=STATUS)
    VERIFY_(STATUS)
    self%wAersl(:,:,km:1:-1,4) = PTR3D(:,:,1:km)*airdens(:,:,1:km)
    IF(self%verbose) CALL pmaxmin('ss001:', PTR3D, qmin, qmax, iXj, km, 1. )
    NULLIFY(PTR3D)

    CALL ESMFL_BundleGetPointertoData(aerosols, 'ss002', PTR3D, RC=STATUS)
    VERIFY_(STATUS)
    self%wAersl(:,:,km:1:-1,4) = self%wAersl(:,:,km:1:-1,4)+PTR3D(:,:,1:km)*airdens(:,:,1:km)
    IF(self%verbose) CALL pmaxmin('ss002:', PTR3D, qmin, qmax, iXj, km, 1. )
    NULLIFY(PTR3D)

    CALL ESMFL_BundleGetPointertoData(aerosols, 'ss003', PTR3D, RC=STATUS)
    VERIFY_(STATUS)
    self%wAersl(:,:,km:1:-1,5) = PTR3D(:,:,1:km)*airdens(:,:,1:km)
    IF(self%verbose) CALL pmaxmin('ss003:', PTR3D, qmin, qmax, iXj, km, 1. )
    NULLIFY(PTR3D)

    CALL ESMFL_BundleGetPointertoData(aerosols, 'ss004', PTR3D, RC=STATUS)
    VERIFY_(STATUS)
    self%wAersl(:,:,km:1:-1,5) = self%wAersl(:,:,km:1:-1,5)+PTR3D(:,:,1:km)*airdens(:,:,1:km)
    IF(self%verbose) CALL pmaxmin('ss004:', PTR3D, qmin, qmax, iXj, km, 1. )
    NULLIFY(PTR3D)

    CALL ESMFL_BundleGetPointertoData(aerosols, 'ss005', PTR3D, RC=STATUS)
    VERIFY_(STATUS)
    self%wAersl(:,:,km:1:-1,5) = self%wAersl(:,:,km:1:-1,5)+PTR3D(:,:,1:km)*airdens(:,:,1:km)
    IF(self%verbose) CALL pmaxmin('ss005:', PTR3D, qmin, qmax, iXj, km, 1. )
    NULLIFY(PTR3D)

   CASE("GOCART")

    IF(self%usingGOCART_SS) THEN

     CALL MAPL_GetPointer(impChem, SS001, 'GOCART::ss001', RC=STATUS)
     VERIFY_(STATUS)
     CALL MAPL_GetPointer(impChem, SS002, 'GOCART::ss002', RC=STATUS)
     VERIFY_(STATUS)
     CALL MAPL_GetPointer(impChem, SS003, 'GOCART::ss003', RC=STATUS)
     VERIFY_(STATUS)
     CALL MAPL_GetPointer(impChem, SS004, 'GOCART::ss004', RC=STATUS)
     VERIFY_(STATUS)
     CALL MAPL_GetPointer(impChem, SS005, 'GOCART::ss005', RC=STATUS)
     VERIFY_(STATUS)

     IF(self%verbose) THEN
      CALL pmaxmin('SS001:', SS001, qmin, qmax, iXj, km, 1. )
      CALL pmaxmin('SS002:', SS002, qmin, qmax, iXj, km, 1. )
      CALL pmaxmin('SS003:', SS003, qmin, qmax, iXj, km, 1. )
      CALL pmaxmin('SS004:', SS004, qmin, qmax, iXj, km, 1. )
      CALL pmaxmin('SS005:', SS005, qmin, qmax, iXj, km, 1. )
     END IF

! Accumulated
! -----------
     self%wAersl(:,:,km:1:-1,4) = (SS001(:,:,1:km)+SS002(:,:,1:km))*airdens(:,:,1:km)

! Coarse
! ------
     self%wAersl(:,:,km:1:-1,5) = (SS003(:,:,1:km)+SS004(:,:,1:km)+SS005(:,:,1:km))*airdens(:,:,1:km)

    END IF

   CASE("GMICHEM")

! Accumulated
! -----------
    CALL MAPL_GetPointer(impChem, PTR3D, 'SS1', RC=STATUS)
    VERIFY_(STATUS)
    self%wAersl(:,:,1:km,4) = PTR3D(:,:,km:1:-1)
    IF(self%verbose) CALL pmaxmin('SS1:', PTR3D, qmin, qmax, iXj, km, 1. )
    IF(ASSOCIATED(SS001) .AND. self%AM_I_AERO_PROVIDER) SS001(:,:,:) = PTR3D(:,:,:)/airdens(:,:,:)
    NULLIFY(PTR3D)

! Coarse
! ------
    self%wAersl(:,:,:,5) = 0.00
    DO i = 2,4

     WRITE(importName,FMT="('SS',I1)") i
     CALL MAPL_GetPointer(impChem, PTR3D, TRIM(importName), RC=STATUS)
     VERIFY_(STATUS)
     self%wAersl(:,:,1:km,5) = self%wAersl(:,:,1:km,5)+PTR3D(:,:,km:1:-1)
     IF(self%verbose) CALL pmaxmin(TRIM(importName)//':', PTR3D, qmin, qmax, iXj, km, 1. )

     SELECT CASE (i)
      CASE(2)
       IF(ASSOCIATED(SS003) .AND. self%AM_I_AERO_PROVIDER) SS003(:,:,:) = PTR3D(:,:,:)/airdens(:,:,:)
      CASE(3)
       IF(ASSOCIATED(SS004) .AND. self%AM_I_AERO_PROVIDER) SS004(:,:,:) = PTR3D(:,:,:)/airdens(:,:,:)
      CASE(4)
       IF(ASSOCIATED(SS005) .AND. self%AM_I_AERO_PROVIDER) SS005(:,:,:) = PTR3D(:,:,:)/airdens(:,:,:)
     END SELECT

     NULLIFY(PTR3D)

    END DO

   CASE DEFAULT

    STATUS = 1
    VERIFY_(STATUS)

  END SELECT

  RETURN
 END SUBROUTINE Acquire_SS

!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  Acquire_SU
!
! !INTERFACE:

  SUBROUTINE Acquire_SU(rc)

  IMPLICIT NONE

  INTEGER, INTENT(OUT) :: rc
!
! !DESCRIPTION:
!
!  Obtain sulfate either by import from GOCART or from an external
!  file specified in ExtData.rc along with a refresh template.
!
!  SO4: self%wAersl(:,:,:,1)
!
!EOP
!---------------------------------------------------------------------------

  CHARACTER(LEN=255) :: IAm
  REAL, POINTER, DIMENSION(:,:,:) :: PTR3D
  REAL :: qmin,qmax

  TYPE(ESMF_State)           :: aero
  TYPE(ESMF_FieldBundle)     :: aerosols

  rc = 0
  IAm = "Acquire_SU"

  SELECT CASE (TRIM(self%aeroProviderName))

   CASE("GOCART.data")

    CALL ESMF_StateGet(impChem, 'AERO', aero, RC=STATUS)
    VERIFY_(STATUS)
    CALL ESMF_StateGet(aero, 'AEROSOLS', aerosols, RC=STATUS)
    VERIFY_(STATUS)

    CALL ESMFL_BundleGetPointertoData(aerosols, 'SO4', PTR3D, RC=STATUS)
    VERIFY_(STATUS)
    self%wAersl(:,:,km:1:-1,1) = PTR3D(:,:,1:km)*airdens(:,:,1:km)
    IF(self%verbose) CALL pmaxmin('SO4:', PTR3D, qmin, qmax, iXj, km, 1. )
    NULLIFY(PTR3D)


   CASE("GOCART")

    IF(self%usingGOCART_SU) THEN

     CALL MAPL_GetPointer(impChem, SO4, 'GOCART::SO4', RC=STATUS)
     VERIFY_(STATUS)
     self%wAersl(:,:,km:1:-1,1) = SO4(:,:,1:km)*airdens(:,:,1:km)

     IF(self%verbose) THEN
      CALL pmaxmin('SO4:', SO4, qmin, qmax, iXj, km, 1. )
     END IF

     ! If volcanic SU exists, use it too:
     NULLIFY(SO4)
     CALL MAPL_GetPointer(impChem, SO4, 'GOCART::SO4v')
     IF( ASSOCIATED(SO4) ) THEN

       self%wAersl(:,:,km:1:-1,1) = &
       self%wAersl(:,:,km:1:-1,1) + SO4(:,:,1:km)*airdens(:,:,1:km)

       IF(self%verbose) THEN
        CALL pmaxmin('SO4v:', SO4, qmin, qmax, iXj, km, 1. )
       END IF

       NULLIFY(SO4)

     END IF

    END IF

   CASE("GMICHEM")

    CALL MAPL_GetPointer(impChem, PTR3D, 'SO4', RC=STATUS)
    VERIFY_(STATUS)
    self%wAersl(:,:,km:1:-1,1) = PTR3D(:,:,1:km)
    IF(self%verbose) CALL pmaxmin('SO4:', PTR3D, qmin, qmax, iXj, km, 1. )
    IF(ASSOCIATED(SO4) .AND. self%AM_I_AERO_PROVIDER) SO4(:,:,:) = PTR3D(:,:,:)/airdens(:,:,:)
    NULLIFY(PTR3D)

   CASE DEFAULT

    STATUS = 1
    VERIFY_(STATUS)

  END SELECT

  RETURN
 END SUBROUTINE Acquire_SU

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

! Solar zenith angle used in photolysis
! -------------------------------------
  IF(ASSOCIATED(SZAPHOT)) THEN
   IF( (self%phot_opt == 3 .OR. self%phot_opt == 4 .OR. self%phot_opt == 5) .AND. &
        self%gotImportRst) SZAPHOT(i1:i2,j1:j2) = solarZenithAngle(i1:i2,j1:j2)
  END IF

! Ship Emisssions
! ---------------
  if (self%do_ShipEmission)  then

     w_c%qa(w_c%reg%j_XX)%data3d(i1:i2,j1:j2,km-1) = self%qjgmi(self%jno2num)%pArray3D(:,:,1)

     ALLOCATE (jNO2val_phot(i1:i2,j1:j2),STAT=STATUS)
     VERIFY_(STATUS)

     jNO2val_phot(:,:) = self%qjgmi(self%jno2num)%pArray3D(:,:,1)
     CALL setDataToStateField(expChem,  jNO2val_phot,  'jNO2val')
     
     DEALLOCATE (jNO2val_phot,STAT=STATUS)
     VERIFY_(STATUS)
  end if

! ------------------------------------------------------------------------
! Generalized diagnostic for dust and aerosol. See GmiAerDustODSA_mod.F90.
! ------------------------------------------------------------------------
   IF(self%gotImportRst .AND. self%phot_opt == 3 .AND. self%do_AerDust_calc) THEN
     ALLOCATE(var4dDBL(i1:i2,j1:j2,1:km,1:num_AerDust),STAT=STATUS)
     VERIFY_(STATUS)
     var4dDBL(:,:,:,:) = self%optDepth(:,:,:,:)
     DO k=1,km
      kReverse = km-k+1
      IF(ASSOCIATED(  DUSTOD))   DUSTOD(i1:i2,j1:j2,kReverse) = var4dDBL(i1:i2,j1:j2,k, 4)
      IF(ASSOCIATED(  DUSTSA))   DUSTSA(i1:i2,j1:j2,kReverse) = var4dDBL(i1:i2,j1:j2,k, 5)
      IF(ASSOCIATED(   SO4OD))    SO4OD(i1:i2,j1:j2,kReverse) = var4dDBL(i1:i2,j1:j2,k, 6)
      IF(ASSOCIATED(SO4HYGRO)) SO4HYGRO(i1:i2,j1:j2,kReverse) = var4dDBL(i1:i2,j1:j2,k, 7)
      IF(ASSOCIATED(   SO4SA))    SO4SA(i1:i2,j1:j2,kReverse) = var4dDBL(i1:i2,j1:j2,k, 8)
      IF(ASSOCIATED(    BCOD))     BCOD(i1:i2,j1:j2,kReverse) = var4dDBL(i1:i2,j1:j2,k, 9)
      IF(ASSOCIATED( BCHYGRO))  BCHYGRO(i1:i2,j1:j2,kReverse) = var4dDBL(i1:i2,j1:j2,k,10)
      IF(ASSOCIATED(    BCSA))     BCSA(i1:i2,j1:j2,kReverse) = var4dDBL(i1:i2,j1:j2,k,11)
      IF(ASSOCIATED(    OCOD))     OCOD(i1:i2,j1:j2,kReverse) = var4dDBL(i1:i2,j1:j2,k,12)
      IF(ASSOCIATED( OCHYGRO))  OCHYGRO(i1:i2,j1:j2,kReverse) = var4dDBL(i1:i2,j1:j2,k,13)
      IF(ASSOCIATED(    OCSA))     OCSA(i1:i2,j1:j2,kReverse) = var4dDBL(i1:i2,j1:j2,k,14)
      IF(ASSOCIATED(   SSAOD))    SSAOD(i1:i2,j1:j2,kReverse) = var4dDBL(i1:i2,j1:j2,k,15)
      IF(ASSOCIATED(SSAHYGRO)) SSAHYGRO(i1:i2,j1:j2,kReverse) = var4dDBL(i1:i2,j1:j2,k,16)
      IF(ASSOCIATED(   SSASA))    SSASA(i1:i2,j1:j2,kReverse) = var4dDBL(i1:i2,j1:j2,k,17)
      IF(ASSOCIATED(   SSCOD))    SSCOD(i1:i2,j1:j2,kReverse) = var4dDBL(i1:i2,j1:j2,k,18)
      IF(ASSOCIATED(SSCHYGRO)) SSCHYGRO(i1:i2,j1:j2,kReverse) = var4dDBL(i1:i2,j1:j2,k,19)
      IF(ASSOCIATED(   SSCSA))    SSCSA(i1:i2,j1:j2,kReverse) = var4dDBL(i1:i2,j1:j2,k,20)
     END DO
     DEALLOCATE(var4dDBL, STAT=STATUS)
     VERIFY_(STATUS)
   END IF

  RETURN
 END SUBROUTINE FillExports

!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  populateBundleQJ
!
! !INTERFACE:

      subroutine populateBundleQJ()
!
      implicit none
!
! !DESCRIPTION:
!
! !LOCAL VARIABLES:
      integer :: STATUS, numVars, ib, rc
      real(rPrec), pointer, dimension(:,:,:)   :: ptr3D
      type(ESMF_FieldBundle)                ::      qjBundle
      type(ESMF_FieldBundle)                ::   tAreaBundle
      type(ESMF_FieldBundle)                :: eRadiusBundle
      character(len=ESMF_MAXSTR), parameter :: IAm = "populateBundleQJ"
!
!EOP
!--------------------------------------------------------------------------------

      allocate(ptr3D(i1:i2, j1:j2, 1:km))

      !==================
      ! For the QJ Bundle
      !==================

      call ESMF_StateGet(expChem, "gmiQJ", qjBundle, rc=STATUS)
      VERIFY_(STATUS)

      call ESMF_FieldBundleGet(qjBundle, fieldCount=numVars, rc=STATUS)
      VERIFY_(STATUS)
      _ASSERT(numVars == self%num_qjo,'needs informative message')

      do ib = 1, numVars
         ptr3D(:,:,:) = self%qjgmi(ib)%pArray3D(:,:,km:1:-1)
         call updateTracerToBundle(qjBundle, ptr3D, ib)
      end do

      !=====================
      ! For the tArea Bundle
      !=====================

      call ESMF_StateGet(expChem, "gmiTAREA", tAreaBundle, rc=STATUS)
      VERIFY_(STATUS)

      call ESMF_FieldBundleGet(tAreaBundle, fieldCount=numVars, rc=STATUS)
      VERIFY_(STATUS)
      _ASSERT(numVars == nSADdust+nSADaer,'needs informative message')

      do ib = 1, numVars
         ptr3D(:,:,:) = self%tArea(:,:,km:1:-1,ib)
         call updateTracerToBundle(tAreaBundle, ptr3D, ib)
      end do

      !=======================
      ! For the eRadius Bundle
      !=======================

      call ESMF_StateGet(expChem, "gmiERADIUS", eRadiusBundle, rc=STATUS)
      VERIFY_(STATUS)

      call ESMF_FieldBundleGet(eRadiusBundle, fieldCount=numVars, rc=STATUS)
      VERIFY_(STATUS)
      _ASSERT(numVars == nSADdust+nSADaer,'needs informative message')

      do ib = 1, numVars
         ptr3D(:,:,:) = self%eRadius(:,:,km:1:-1,ib)
         call updateTracerToBundle(eRadiusBundle, ptr3D, ib)
      end do

      deallocate (ptr3D)

      return

      end subroutine populateBundleQJ

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
  INTEGER :: STATUS
  
  rc=0
  IAm="FindPointers"

!  Pointers to imports
!  -------------------
   CALL MAPL_GetPointer(impChem,    cldtt,    'CLDTT', RC=STATUS)
   VERIFY_(STATUS)
   CALL MAPL_GetPointer(impChem,    albvf,    'ALBVF', RC=STATUS)
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
   CALL MAPL_GetPointer(impChem,      fcld,    'FCLD', RC=STATUS)
   VERIFY_(STATUS)
   CALL MAPL_GetPointer(impChem,    taucli,  'TAUCLI', RC=STATUS)
   VERIFY_(STATUS)
   CALL MAPL_GetPointer(impChem,    tauclw,  'TAUCLW', RC=STATUS)
   VERIFY_(STATUS)
   CALL MAPL_GetPointer(impChem,	ql,	 'QL', RC=STATUS)
   VERIFY_(STATUS)
   CALL MAPL_GetPointer(impChem,   cnv_mfc, 'CNV_MFC', RC=STATUS)
   VERIFY_(STATUS)
   CALL MAPL_GetPointer(impChem,       rh2,	'RH2', RC=STATUS)
   VERIFY_(STATUS)
   CALL MAPL_GetPointer(impChem,      dqdt,    'DQDT', RC=STATUS)
   VERIFY_(STATUS)

!  Export state pointers
!  ---------------------
   CALL MAPL_GetPointer(expChem,  SZAPHOT,  'SZAPHOT', RC=STATUS)
   VERIFY_(STATUS)
   CALL MAPL_GetPointer(expChem,   DUSTOD,   'DUSTOD', RC=STATUS)
   VERIFY_(STATUS)
   CALL MAPL_GetPointer(expChem,   DUSTSA,   'DUSTSA', RC=STATUS)
   VERIFY_(STATUS)
   CALL MAPL_GetPointer(expChem,    SO4OD,    'SO4OD', RC=STATUS)
   VERIFY_(STATUS)
   CALL MAPL_GetPointer(expChem, SO4HYGRO, 'SO4HYGRO', RC=STATUS)
   VERIFY_(STATUS)
   CALL MAPL_GetPointer(expChem,    SO4SA,    'SO4SA', RC=STATUS)
   VERIFY_(STATUS)
   CALL MAPL_GetPointer(expChem,     BCOD,     'BCOD', RC=STATUS)
   VERIFY_(STATUS)
   CALL MAPL_GetPointer(expChem,  BCHYGRO,  'BCHYGRO', RC=STATUS)
   VERIFY_(STATUS)
   CALL MAPL_GetPointer(expChem,     BCSA,     'BCSA', RC=STATUS)
   VERIFY_(STATUS)
   CALL MAPL_GetPointer(expChem,     OCOD,     'OCOD', RC=STATUS)
   VERIFY_(STATUS)
   CALL MAPL_GetPointer(expChem,  OCHYGRO,  'OCHYGRO', RC=STATUS)
   VERIFY_(STATUS)
   CALL MAPL_GetPointer(expChem,     OCSA,     'OCSA', RC=STATUS)
   VERIFY_(STATUS)
   CALL MAPL_GetPointer(expChem,    SSAOD,    'SSAOD', RC=STATUS)
   VERIFY_(STATUS)
   CALL MAPL_GetPointer(expChem, SSAHYGRO, 'SSAHYGRO', RC=STATUS)
   VERIFY_(STATUS)
   CALL MAPL_GetPointer(expChem,    SSASA,    'SSASA', RC=STATUS)
   VERIFY_(STATUS)
   CALL MAPL_GetPointer(expChem,    SSCOD,    'SSCOD', RC=STATUS)
   VERIFY_(STATUS)
   CALL MAPL_GetPointer(expChem, SSCHYGRO, 'SSCHYGRO', RC=STATUS)
   VERIFY_(STATUS)
   CALL MAPL_GetPointer(expChem,    SSCSA,    'SSCSA', RC=STATUS)
   VERIFY_(STATUS)
   
!  Validation
!  ----------
   Validate: IF(self%verbose) THEN
    IF(MAPL_AM_I_ROOT( ))PRINT *,TRIM(IAm),": Input ..."
    CALL pmaxmin('Q:', Q, qmin, qmax, iXj, km, 1. )
    CALL pmaxmin('T:', T, qmin, qmax, iXj, km, 1. )
    CALL pmaxmin('ZLE:', zle, qmin, qmax, iXj, km+1, 1. )
    CALL pmaxmin('PLE (hPa):', ple, qmin, qmax, iXj, km+1, 0.01 )
    CALL pmaxmin('AIRDENS:', airdens, qmin, qmax, iXj, km, 1. )
    CALL pmaxmin('FCLD:', fcld, qmin, qmax, iXj, km, 1. )
    CALL pmaxmin('TAUCLI:', taucli, qmin, qmax, iXj, km, 1. )
    CALL pmaxmin('TAUCLW:', tauclw, qmin, qmax, iXj, km, 1. )
    CALL pmaxmin('QL:', ql, qmin, qmax, iXj, km, 1. )
    CALL pmaxmin('CNV_MFC:', cnv_mfc, qmin, qmax, iXj, km+1, 1. )
    CALL pmaxmin('RH2:', rh2, qmin, qmax, iXj, km, 1. )
    CALL pmaxmin('DQ/DT:', dqdt, qmin, qmax, iXj, km, 1. )

    i = w_c%reg%j_XX
    CALL pmaxmin('TROPP:', w_c%qa(i)%data3d(:,:,km), qmin, qmax, iXj, 1, 0.01 )
    CALL pmaxmin('CLDTT:', cldtt, qmin, qmax, iXj, 1, 1. )
    CALL pmaxmin('ALBEDO:', albvf, qmin, qmax, iXj, 1, 1. )

   END IF Validate

!  Dust and aerosols are (instead!) part of the export state when GMICHEM is the AERO provider.
!  --------------------------------------------------------------------------------------------
   IF(self%AM_I_AERO_PROVIDER) THEN
    CALL MAPL_GetPointer(expChem, BCphobic, 'GMICHEM::BCphobic', ALLOC=.TRUE., RC=STATUS)
    VERIFY_(STATUS)
    CALL MAPL_GetPointer(expChem, BCphilic, 'GMICHEM::BCphilic', ALLOC=.TRUE., RC=STATUS)
    VERIFY_(STATUS)
    CALL MAPL_GetPointer(expChem,    DU001,    'GMICHEM::du001', ALLOC=.TRUE., RC=STATUS)
    VERIFY_(STATUS)
    CALL MAPL_GetPointer(expChem,    DU002,    'GMICHEM::du002', ALLOC=.TRUE., RC=STATUS)
    VERIFY_(STATUS)
    CALL MAPL_GetPointer(expChem,    DU003,    'GMICHEM::du003', ALLOC=.TRUE., RC=STATUS)
    VERIFY_(STATUS)
    CALL MAPL_GetPointer(expChem,    DU004,    'GMICHEM::du004', ALLOC=.TRUE., RC=STATUS)
    VERIFY_(STATUS)
    CALL MAPL_GetPointer(expChem, OCphobic, 'GMICHEM::OCphobic', ALLOC=.TRUE., RC=STATUS)
    VERIFY_(STATUS)
    CALL MAPL_GetPointer(expChem, OCphilic, 'GMICHEM::OCphilic', ALLOC=.TRUE., RC=STATUS)
    VERIFY_(STATUS)
    CALL MAPL_GetPointer(expChem,    SS001,    'GMICHEM::ss001', ALLOC=.TRUE., RC=STATUS)
    VERIFY_(STATUS)
    CALL MAPL_GetPointer(expChem,    SS003,    'GMICHEM::ss003', ALLOC=.TRUE., RC=STATUS)
    VERIFY_(STATUS)
    CALL MAPL_GetPointer(expChem,    SS004,    'GMICHEM::ss004', ALLOC=.TRUE., RC=STATUS)
    VERIFY_(STATUS)
    CALL MAPL_GetPointer(expChem,    SS005,    'GMICHEM::ss005', ALLOC=.TRUE., RC=STATUS)
    VERIFY_(STATUS)
    CALL MAPL_GetPointer(expChem,      SO4,      'GMICHEM::SO4', ALLOC=.TRUE., RC=STATUS)
    VERIFY_(STATUS)
   END IF

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
  i = w_c%reg%j_XX
  tropopausePress(i1:i2,j1:j2) = w_c%qa(i)%data3d(i1:i2,j1:j2,km)*Pa2hPa    ! Pa               hPa
  pctm2(i1:i2,j1:j2) = ple(i1:i2,j1:j2,km)*Pa2hPa                           ! Pa               hPa
  fracCloudCover(i1:i2,j1:j2) = cldtt(i1:i2,j1:j2)                          ! fraction

! Layer means                                                               GEOS-5 Units       GMI Units
! -----------                                                               ------------       -------------
  DO k=1,km
   kReverse = km-k+1                                                        ! Lid-to-surf      Surf-to-lid
   press3c(i1:i2,j1:j2,kReverse) = pl(i1:i2,j1:j2,k)*Pa2hPa                 ! Pa               hPa
   kel(i1:i2,j1:j2,kReverse) = T(i1:i2,j1:j2,k)                             ! K
   humidity(i1:i2,j1:j2,kReverse) = Q(i1:i2,j1:j2,k)*ToGrPerKg              ! kg kg^{-1}       g kg^{-1}
   totalCloudFraction(i1:i2,j1:j2,kReverse) = fcld(i1:i2,j1:j2,k)           ! fraction
   clwc(i1:i2,j1:j2,kReverse) = ql(i1:i2,j1:j2,k)*ToGrPerKg                 ! kg kg^{-1}       g kg^{-1}
   relativeHumidity(i1:i2,j1:j2,kReverse) = rh2(i1:i2,j1:j2,k)              ! fraction
   moistq(i1:i2,j1:j2,kReverse) = dqdt(i1:i2,j1:j2,k)*ToGrPerKg*secPerDay   ! kg kg^{-1}s^{-1} g kg^{-1}d^{-1}
   tau_cloud(i1:i2,j1:j2,kReverse) = (tauclw(i1:i2,j1:j2,k)+  &             ! 1
                                      taucli(i1:i2,j1:j2,k))*fcld(i1:i2,j1:j2,k)

!   tau_clw(i1:i2,j1:j2,kReverse) = tauclw(i1:i2,j1:j2,k)
!   tau_cli(i1:i2,j1:j2,kReverse) = taucli(i1:i2,j1:j2,k)
   tau_clw(i1:i2,j1:j2,kReverse) = tauclw(i1:i2,j1:j2,k)*fcld(i1:i2,j1:j2,k)
   tau_cli(i1:i2,j1:j2,kReverse) = taucli(i1:i2,j1:j2,k)*fcld(i1:i2,j1:j2,k)
  END DO

! These bounds are in Jules' RH code
! ----------------------------------
  WHERE(relativeHumidity < 0.00D+00) relativeHumidity = 0.00D+00
  WHERE(relativeHumidity > 0.95D+00) relativeHumidity = 0.95D+00

! Layer edges                                                               GEOS-5 Units       GMI Units
! -----------                                                               ------------       -------------
  DO k=0,km
   kReverse = km-k
   press3e(i1:i2,j1:j2,kReverse) = ple(i1:i2,j1:j2,k)*Pa2hPa                ! Pa               hPa
  END DO

  DO k=0,km-1
   kReverse = km-k
   cmf(i1:i2,j1:j2,kReverse) = cnv_mfc(i1:i2,j1:j2,k)                       ! kg m^{-2}s^{-1}
  END DO


! Surface Vis/UV albedoes are MAPL_UNDEF at night.  They
! are NOT weighted by cosine solar zenith angle during daytime.
! -------------------------------------------------------------
  surf_alb(i1:i2,j1:j2) = albvf(i1:i2,j1:j2)
  WHERE(surf_alb > 2.00D+00) surf_alb = 0.00D+00
  WHERE(surf_alb > 1.00D+00) surf_alb = 1.00D+00
  WHERE(surf_alb < 0.00D+00) surf_alb = 0.00D+00

! Cell mass and thickness                                                   GEOS-5 Units       GMI Units
! -----------------------                                                   ------------       -------------
   DO k=1,km
    kReverse = km-k+1
    mass(:,:,kReverse)=airdens(:,:,k)*self%cellArea(:,:)* &                 ! kg
                       (zle(:,:,k-1)-zle(:,:,k))
    gridBoxThickness(:,:,kReverse) = zle(:,:,k-1)-zle(:,:,k)                ! m
   END DO

! Set the dust and aerosol diagnostic to zero
! -------------------------------------------
   IF(TRIM(self%chem_mecha) ==        'troposphere' .OR. &
      TRIM(self%chem_mecha) ==         'strat_trop' .OR. &
      TRIM(self%chem_mecha) == 'strat_trop_aerosol') THEN
    IF(self%phot_opt == 3 .AND. self%do_AerDust_Calc) THEN
     self%optDepth(:,:,:,:)=0.00D+00
    END IF
   END IF

  RETURN
 END SUBROUTINE SatisfyImports

 END SUBROUTINE GmiPhotolysis_GridCompRun

!-------------------------------------------------------------------------
!NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1, GEOS/DAS!
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  GmiPhotolysis_GridCompFinalize
!
! !INTERFACE:
!

   SUBROUTINE GmiPhotolysis_GridCompFinalize ( self, w_c, impChem, expChem, &
                                     nymd, nhms, cdt, rc )

  IMPLICIT none

! !INPUT/OUTPUT PARAMETERS:

   TYPE(GmiPhotolysis_GridComp), INTENT(inout) :: self ! Grid Component

! !INPUT PARAMETERS:

   TYPE(Chem_Bundle), INTENT(in)  :: w_c      ! Chemical tracer fields   
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

   CHARACTER(LEN=*), PARAMETER :: IAm = 'GmiPhotolysis_GridCompFinalize'
   INTEGER :: STATUS
   rc = 0
   DEALLOCATE(self%cellArea, self%latRad, self%lonRad, STAT=STATUS)
   VERIFY_(STATUS)

   RETURN

 END SUBROUTINE GmiPhotolysis_GridCompFinalize
  
 END MODULE GmiPhotolysis_GCCMod

