#include "MAPL_Exceptions.h"
#include "MAPL_Generic.h"

!=============================================================================
!BOP

! !MODULE: GEOS_ChemGridCompMod -- Parent Aerosol/Chemistry Component

! !INTERFACE:

module GEOS_ChemGridCompMod

! !USES:

  use ESMF
  use MAPL
  use Chem_Mod
  use Chem_UtilMod
  use Bundle_IncrementMod

  use  GEOS_ChemEnvGridCompMod,  only :   ChemEnv_SetServices => SetServices
  use       GOCART_GridCompMod,  only :    GOCART_SetServices => SetServices
  use     GOCART2G_GridCompMod,  only :  GOCART2G_SetServices => SetServices !GOCART REFACTOR
  use    QuickChem_GridCompMod,  only : QUICKCHEM_SetServices => SetServices
  use    StratChem_GridCompMod,  only : StratChem_SetServices => SetServices
  use      GMIchem_GridCompMod,  only :       GMI_SetServices => SetServices
  use    CARMAchem_GridCompMod,  only :     CARMA_SetServices => SetServices
  use GEOSCHEMchem_GridCompMod,  only :    GCChem_SetServices => SetServices
  use   MATRIXchem_GridCompMod,  only :    MATRIX_SetServices => SetServices
  use      MAMchem_GridCompMod,  only :       MAM_SetServices => SetServices
  use    GEOS_PChemGridCompMod,  only :     PChem_SetServices => SetServices
  use    GEOS_AChemGridCompMod,  only :     AChem_SetServices => SetServices
  use         GAAS_GridCompMod,  only :      GAAS_SetServices => SetServices
  use          H2O_GridCompMod,  only :       H2O_SetServices => SetServices
  use           TR_GridCompMod,  only :        TR_SetServices => SetServices
  use          DNA_GridCompMod,  only :       DNA_SetServices => SetServices
  use        HEMCO_GridCompMod,  only :     HEMCO_SetServices => SetServices

  use    QuickChem_GridCompMod,  only : IS_QC_INSTANCE_RUNNING
  use     GOCART2G_GridCompMod,  only : IS_G2G_INSTANCE_RUNNING

  implicit none
  private

! !PUBLIC MEMBER FUNCTIONS:

  public SetServices

!                                             -----------
 
! Private state
! -------------
  TYPE GEOS_ChemGridComp
     PRIVATE
     type(ESMF_Config), pointer :: CF  ! Private Config
     LOGICAL :: enable_PCHEM
     LOGICAL :: enable_ACHEM
     LOGICAL :: enable_GOCART
     LOGICAL :: enable_GOCARTdata
     LOGICAL :: enable_GOCART2G        ! GOCART REFACTOR
     LOGICAL :: enable_QUICKCHEM
     LOGICAL :: enable_GAAS
     LOGICAL :: enable_H2O
     LOGICAL :: enable_STRATCHEM
     LOGICAL :: enable_GMICHEM
     LOGICAL :: enable_CARMA
     LOGICAL :: enable_GEOSCHEM
     LOGICAL :: enable_MATRIX
     LOGICAL :: enable_MAM
     LOGICAL :: enable_MAMdata
     LOGICAL :: enable_TR
     LOGICAL :: enable_DNA
     LOGICAL :: enable_HEMCO
     LOGICAL :: strict_child_timing    ! Add a barrier after each child is called
     INTEGER :: AERO_PROVIDER
     INTEGER :: RATS_PROVIDER          ! WARNING: May be multiple RATS_PROVIDERs 

     LOGICAL :: running_BC
     LOGICAL :: running_OC
     LOGICAL :: running_BR
     LOGICAL :: running_DU
     LOGICAL :: running_SU
     LOGICAL :: running_SS
     LOGICAL :: running_NI
     LOGICAL :: running_CO
     LOGICAL :: running_CH4
     LOGICAL :: running_OH
  END TYPE GEOS_ChemGridComp

! Hook for the ESMF
! -----------------
  TYPE GEOS_ChemGridComp_Wrap
     TYPE (GEOS_ChemGridComp), pointer :: PTR => null()
  END TYPE GEOS_ChemGridComp_Wrap

!=============================================================================

! !DESCRIPTION: This gridded component (GC) combines 
 
!EOP

! IMPORTANT: If adding a new component, make sure to update private function GetProvider_()
! -----------------------------------------------------------------------------------------
  integer ::      CHEMENV = -1
  integer ::        HEMCO = -1
  integer ::        PCHEM = -1
  integer ::        ACHEM = -1
  integer ::       GOCART = -1
  integer ::     GOCART2G = -1
  integer ::   GOCARTdata = -1
  integer ::    QUICKCHEM = -1
  integer ::         GAAS = -1
  integer ::          H2O = -1
  integer ::    STRATCHEM = -1
  integer ::      GMICHEM = -1
  integer ::        CARMA = -1
  integer ::     GEOSCHEM = -1
  integer ::       MATRIX = -1
  integer ::          MAM = -1
  integer ::      MAMdata = -1
  integer ::           TR = -1
  integer ::          DNA = -1

  INTEGER, PARAMETER :: DBL = KIND(0.00D+00)


contains

!BOP

! !IROUTINE: SetServices -- Sets ESMF services for this component

! !INTERFACE:

    subroutine SetServices ( GC, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
    integer,             intent(  OUT) :: RC  ! return code

! !DESCRIPTION:  The SetServices for the Chemistry GC needs to register its
!   Initialize and Run.  It uses the MAPL\_Generic construct for defining 
!   state specs and couplings among its children.  In addition, it creates the   
!   children GCs and runs their respective SetServices.

!EOP

!=============================================================================
!
! ErrLog Variables

    __Iam__('SetServices')      ! NOTE: this macro declares STATUS
                                ! ALSO: Never set Iam = TRIM(Iam) // suffix
                                !       because Iam is a SAVED varaible
    character(len=ESMF_MAXSTR) :: COMP_NAME

! Locals

   type (ESMF_GridComp), pointer :: GCS(:)

    integer                    :: I, RATS_PROVIDER, AERO_PROVIDER
    type (ESMF_Config), target :: CF, myCF

    integer                    :: n, id
    
    INTEGER, PARAMETER         :: numRATs = 8
    INTEGER                    :: RATsProviderNumber(numRATs)
    CHARACTER(LEN=ESMF_MAXSTR) :: RATsProviderName(numRATs)
    CHARACTER(LEN=ESMF_MAXSTR) :: speciesName(numRATs) = (/ "OX    ", "O3    ", "O3PPMV", "CH4   ", &
                                                            "N2O   ", "CFC11 ", "CFC12 ", "HCFC22"/)
    CHARACTER(LEN=ESMF_MAXSTR) :: providerName
    CHARACTER(LEN=ESMF_MAXSTR) :: shortName
    CHARACTER(LEN=ESMF_MAXSTR) :: str

    !GMI MEGAN isoprene related 
    CHARACTER(LEN=255) :: gmi_rcfilen = 'GMI_GridComp.rc'
    TYPE(ESMF_Config)  :: gmi_config
    LOGICAL  doMEGANemission,  doMEGANviaHEMCO

!   Private state
!   -------------
    type (GEOS_ChemGridComp), pointer :: myState   ! private, that is
    type (GEOS_ChemGridComp_Wrap)     :: wrap
    type(Chem_Registry)               :: chemReg
    character(len=ESMF_MAXSTR) :: chem_gridcomp_rc_file

!   Related to HEMCO
!   ----------------
    type (ESMF_Config)          :: HemcoCF
    character(len=ESMF_MAXSTR)  :: Label
    character(len=ESMF_MAXSTR)  :: ConfigFile
    integer                     :: nnInst
    logical                     :: GOCART_instance_of_HEMCO   ! TRUE if HEMCO is running a GOCART instance
    logical                     ::    GMI_instance_of_HEMCO   ! TRUE if HEMCO is running a GMI instance

!=============================================================================

! Begin...

! Get my name and set-up traceback handle
! ---------------------------------------
    call ESMF_GridCompGet( GC, NAME=COMP_NAME, CONFIG=CF, __RC__ )
    Iam = TRIM(COMP_NAME) // '::' //TRIM(Iam)

!   Wrap internal state for storing in GC; rename legacyState
!   ---------------------------------------------------------
    allocate ( myState, __STAT__ )
    wrap%ptr => myState

!   Load the Chemistry Registry
!   ---------------------------
    chemReg = Chem_RegistryCreate ( STATUS )
    VERIFY_(STATUS)

! Register services for this component
! ------------------------------------
    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_INITIALIZE, Init, __RC__ )
    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN, Run1, __RC__ )
    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN, Run2, __RC__ )

!   Store private state in GC
!   -------------------------
    call ESMF_UserCompSetInternalState ( GC, 'GEOSchem_GridComp_State', wrap, STATUS )
    VERIFY_(STATUS)

    call ESMF_ConfigGetAttribute(cf, chem_gridcomp_rc_file, label = "GEOS_ChemGridComp_RC_File:", &
         default = "GEOS_ChemGridComp.rc", __RC__)
  
! Identify which children to run
! ------------------------------
    myCF = ESMF_ConfigCreate(__RC__)

    call ESMF_ConfigLoadFile ( myCF, chem_gridcomp_rc_file, __RC__ )

    call ESMF_ConfigGetAttribute(myCF, myState%enable_PCHEM,        Default=.FALSE., Label="ENABLE_PCHEM:",        __RC__ )
    call ESMF_ConfigGetAttribute(myCF, myState%enable_ACHEM,        Default=.FALSE., Label="ENABLE_ACHEM:",        __RC__ )
    call ESMF_ConfigGetAttribute(myCF, myState%enable_GOCART,       Default=.FALSE., Label="ENABLE_GOCART:",       __RC__ )
    call ESMF_ConfigGetAttribute(myCF, myState%enable_GOCARTdata,   Default=.FALSE., Label="ENABLE_GOCART_DATA:",  __RC__ )
    call ESMF_ConfigGetAttribute(myCF, myState%enable_GOCART2G,     Default=.FALSE., Label="ENABLE_GOCART2G:",     __RC__ )
    call ESMF_ConfigGetAttribute(myCF, myState%enable_QUICKCHEM,    Default=.FALSE., Label="ENABLE_QUICKCHEM:",    __RC__ )
    call ESMF_ConfigGetAttribute(myCF, myState%enable_GAAS,         Default=.FALSE., Label="ENABLE_GAAS:",         __RC__ )
    call ESMF_ConfigGetAttribute(myCF, myState%enable_H2O,          Default=.FALSE., Label="ENABLE_H2O:",          __RC__ )
    call ESMF_ConfigGetAttribute(myCF, myState%enable_STRATCHEM,    Default=.FALSE., Label="ENABLE_STRATCHEM:",    __RC__ )
    call ESMF_ConfigGetAttribute(myCF, myState%enable_GMICHEM,      Default=.FALSE., Label="ENABLE_GMICHEM:",      __RC__ )
    call ESMF_ConfigGetAttribute(myCF, myState%enable_CARMA,        Default=.FALSE., Label="ENABLE_CARMA:",        __RC__ )
    call ESMF_ConfigGetAttribute(myCF, myState%enable_GEOSCHEM,     Default=.FALSE., Label="ENABLE_GEOSCHEM:",     __RC__ )
    call ESMF_ConfigGetAttribute(myCF, myState%enable_MATRIX,       Default=.FALSE., Label="ENABLE_MATRIX:",       __RC__ )
    call ESMF_ConfigGetAttribute(myCF, myState%enable_MAM,          Default=.FALSE., Label="ENABLE_MAM:",          __RC__ )
    call ESMF_ConfigGetAttribute(myCF, myState%enable_MAMdata,      Default=.FALSE., Label="ENABLE_MAM_DATA:",     __RC__ )
    call ESMF_ConfigGetAttribute(myCF, myState%enable_TR,           Default=.FALSE., Label="ENABLE_TR:",           __RC__ )
    call ESMF_ConfigGetAttribute(myCF, myState%enable_DNA,          Default=.FALSE., Label="ENABLE_DNA:",          __RC__ )
    call ESMF_ConfigGetAttribute(myCF, myState%enable_HEMCO,        Default=.FALSE., Label="ENABLE_HEMCO:",        __RC__ )

    call ESMF_ConfigGetAttribute(myCF, myState%strict_child_timing, Default=.FALSE., Label="strict_child_timing:", __RC__ )

!ALT: valgrind flagged a memory leak.    myState%CF => myCF ! save for later
    call ESMF_ConfigDestroy(myCF, __RC__)

! GOCART instances:
! -----------------
    myState%running_CO  = myState%enable_GOCART   .AND. chemReg%doing_CO    ! Make sure to overwrite the AMIP copy of the registry
    myState%running_CH4 = myState%enable_GOCART   .AND. chemReg%doing_CH4   ! Make sure to overwrite the AMIP copy of the registry

! GOCART2G instances:
! -------------------
    if ( myState%enable_GOCART2G ) then
      call IS_G2G_INSTANCE_RUNNING ( aerosol_name="CA", instance_name="CA.bc", running=myState%running_BC, __RC__ )
      call IS_G2G_INSTANCE_RUNNING ( aerosol_name="CA", instance_name="CA.oc", running=myState%running_OC, __RC__ )
      call IS_G2G_INSTANCE_RUNNING ( aerosol_name="CA", instance_name="CA.br", running=myState%running_BR, __RC__ )
      call IS_G2G_INSTANCE_RUNNING ( aerosol_name="DU", instance_name="DU"   , running=myState%running_DU, __RC__ )
      call IS_G2G_INSTANCE_RUNNING ( aerosol_name="SU", instance_name="SU"   , running=myState%running_SU, __RC__ )
      call IS_G2G_INSTANCE_RUNNING ( aerosol_name="SS", instance_name="SS"   , running=myState%running_SS, __RC__ )
      call IS_G2G_INSTANCE_RUNNING ( aerosol_name="NI", instance_name="NI"   , running=myState%running_NI, __RC__ )
    else
      myState%running_BC = .FALSE.
      myState%running_OC = .FALSE.
      myState%running_BR = .FALSE.
      myState%running_DU = .FALSE.
      myState%running_SU = .FALSE.
      myState%running_SS = .FALSE.
      myState%running_NI = .FALSE.
    end if

! QUICKCHEM instances:
! --------------------
    if ( myState%enable_QUICKCHEM ) then
      call IS_QC_INSTANCE_RUNNING ( species_name="OH", instance_name="OH", running=myState%running_OH, __RC__ )
    else
      myState%running_OH = .FALSE.
    end if

! Sanity checks:
! --------------
    if (myState%enable_GAAS) then
       _ASSERT(myState%enable_GOCART2G,'when GAAS is enabled GOCART-2G must be enabled as well.')
    end if

! Create children's gridded components and invoke their SetServices
! -----------------------------------------------------------------
    CHEMENV = MAPL_AddChild(GC, NAME='CHEMENV', SS=ChemEnv_SetServices, __RC__)

    if (myState%enable_HEMCO     )       HEMCO = MAPL_AddChild(GC, NAME='HEMCO',        SS=HEMCO_SetServices,     __RC__)
    if (myState%enable_PCHEM     )       PCHEM = MAPL_AddChild(GC, NAME='PCHEM',        SS=PChem_SetServices,     __RC__)
    if (myState%enable_ACHEM     )       ACHEM = MAPL_AddChild(GC, NAME='ACHEM',        SS=AChem_SetServices,     __RC__)
    if (myState%enable_GOCART    )      GOCART = MAPL_AddChild(GC, NAME='GOCART',       SS=GOCART_SetServices,    __RC__)
    if (myState%enable_GOCARTdata)  GOCARTdata = MAPL_AddChild(GC, NAME='GOCART.data',  SS=GOCART_SetServices,    __RC__)
    if (myState%enable_GOCART2G  )    GOCART2G = MAPL_AddChild(GC, NAME='GOCART2G',     SS=GOCART2G_SetServices,  __RC__)
    if (myState%enable_QUICKCHEM )   QUICKCHEM = MAPL_AddChild(GC, NAME='QUICKCHEM',    SS=QUICKCHEM_SetServices, __RC__)
    if (myState%enable_GAAS      )        GAAS = MAPL_AddChild(GC, NAME='GAAS',         SS=GAAS_SetServices,      __RC__)
    if (myState%enable_H2O       )         H2O = MAPL_AddChild(GC, NAME='H2O',          SS=H2O_SetServices,       __RC__)
    if (myState%enable_STRATCHEM )   STRATCHEM = MAPL_AddChild(GC, NAME='STRATCHEM',    SS=StratChem_SetServices, __RC__)
    if (myState%enable_GMICHEM   )     GMICHEM = MAPL_AddChild(GC, NAME='GMICHEM',      SS=GMI_SetServices,       __RC__)
    if (myState%enable_CARMA     )       CARMA = MAPL_AddChild(GC, NAME='CARMA',        SS=CARMA_SetServices,     __RC__)
    if (myState%enable_GEOSCHEM  )    GEOSCHEM = MAPL_AddChild(GC, NAME='GEOSCHEMCHEM', SS=GCChem_SetServices,    __RC__)
    if (myState%enable_MATRIX    )      MATRIX = MAPL_AddChild(GC, NAME='MATRIX',       SS=MATRIX_SetServices,    __RC__)
    if (myState%enable_MAM       )         MAM = MAPL_AddChild(GC, NAME='MAM',          SS=MAM_SetServices,       __RC__)
    if (myState%enable_MAMdata   )     MAMdata = MAPL_AddChild(GC, NAME='MAM.data',     SS=MAM_SetServices,       __RC__)
    if (myState%enable_TR        )          TR = MAPL_AddChild(GC, NAME='TR',           SS=TR_SetServices,        __RC__)
    if (myState%enable_DNA       )         DNA = MAPL_AddChild(GC, NAME='DNA',          SS=DNA_SetServices,       __RC__)


! A container for the friendly tracers
! ------------------------------------
    call MAPL_AddInternalSpec(GC,                                  &
         SHORT_NAME = 'CHEM_TRACERS',                              &
         LONG_NAME  = 'chemistry_friendly_tracers',                &
         UNITS      = 'X',                                         &
         DATATYPE   = MAPL_BundleItem,                             &
         __RC__  )

! A container for total chemistry tracer increments
! -------------------------------------------
    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME = 'CHEMTRI',                                   &
         LONG_NAME  = 'chemistry_friendly_tracers_increments',     &
         UNITS      = 'X',                                         &
         DATATYPE   = MAPL_BundleItem,                             &
         __RC__  )


! A container for chemistry tracer increments for Run1.
! ---------------------------------------------------------
    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME = 'CHEMTRIr1',                                 &
         LONG_NAME  = 'chemistry_friendly_tracers_increments_emissions',    &
         UNITS      = 'X',                                         &
         DATATYPE   = MAPL_BundleItem,                             &
         __RC__  )


! Radiatively Active Tracers (RATs).  If a RATS_PROVIDER is not
! specified in the AGCM.rc, then the provider defaults to PCHEM.
! --------------------------------------------------------------
  call GetProvider_(CF, Label='RATS_PROVIDER:', ID=RATS_PROVIDER, Name=providerName, Default='PCHEM', __RC__)

  RATsProviderNumber(1:numRATs) = RATS_PROVIDER
  RATsProviderName(1:numRATs)   = trim(providerName)

! Assume the H2O Tendency is available from the above provider
! ------------------------------------------------------------
  CALL MAPL_AddExportSpec ( GC, SHORT_NAME = 'H2O_TEND', &
                            CHILD_ID = RATS_PROVIDER, __RC__ )

! Priority for first three RATs, OX, O3 and O3PPMV, goes to the ANALYSIS_OX_PROVIDER.
! -----------------------------------------------------------------------------------
  call GetProvider_(CF, Label='ANALYSIS_OX_PROVIDER:', ID=i, Name=providerName, Default='PCHEM', __RC__)

  RATsProviderNumber(1:3) = i
  RATsProviderName(1:3)   = trim(providerName)

! Add export specs for the RATs ...
! ---------------------------------
  DO i = 1, numRATs
   CALL MAPL_AddExportSpec( GC, SHORT_NAME = TRIM(speciesName(i)), &
                            CHILD_ID = RATsProviderNumber(i), __RC__ )
  END DO

  IF(MAPL_AM_I_ROOT()) THEN
   PRINT *," "
   PRINT *, TRIM(Iam)//": RATs Provider List" 
    DO i = 1, numRATs
     PRINT *,"  "//TRIM(speciesName(i))//": "//TRIM(RATsProviderName(i))
    END DO
  END IF

! Choose the ozone mixing ratio tendency from same child as OX.
! -------------------------------------------------------------
  CALL MAPL_AddExportSpec ( GC, SHORT_NAME = 'OX_TEND', &
                            CHILD_ID = RATsProviderNumber(1), __RC__ )

! Aerosol for radiation.  If an AERO_PROVIDER is not specified
! in the AGCM.tmpl, then the provider defaults to GOCART2G.
! -----------------------------------------------------------

  call ESMF_ConfigGetAttribute(CF, providerName, Default='GOCART2G', &
                               Label="AERO_PROVIDER:", __RC__ )

  str = trim(providerName)
  str = ESMF_UtilStringLowerCase(str, __RC__)

  if (trim(str) == 'none') then
      AERO_PROVIDER = -1

      call MAPL_AddExportSpec(GC,                                        &
                              SHORT_NAME = 'AERO',                       &
                              LONG_NAME  = 'aerosol_mass_mixing_ratios', &
                              UNITS      = 'kg kg-1',                    &
                              DIMS       = MAPL_DimsHorzVert,            &
                              VLOCATION  = MAPL_VLocationCenter,         &
                              DATATYPE   = MAPL_StateItem, __RC__ )
#ifdef ENABLE_AERO_ACI
      ! IMPORTANT: This feature is currently disabled in Physics above,
      !            waiting for patches in the ESMF.
      call MAPL_AddExportSpec(GC,                                &
                              SHORT_NAME = 'AERO_ACI',                   &
                              LONG_NAME  = 'aerosol_mass_mixing_ratios', &
                              UNITS      = 'kg kg-1',                    &
                              DIMS       = MAPL_DimsHorzVert,            &
                              VLOCATION  = MAPL_VLocationCenter,         &
                              DATATYPE   = MAPL_StateItem, __RC__ )

#endif
      
      call MAPL_AddExportSpec(GC,                                &
                              SHORT_NAME = 'AERO_DP',            &
                              LONG_NAME  = 'aerosol_deposition', &
                              UNITS      = 'kg m-2 s-1',         &
                              DIMS       = MAPL_DimsHorzOnly,    &
                              DATATYPE   = MAPL_BundleItem, __RC__) 
  else

      ! Determine Id of the aerosol provider
      ! ------------------------------------
      call GetProvider_(CF, Label='AERO_PROVIDER:', ID=AERO_PROVIDER,  &
                            Name=providerName, Default='GOCART2G', __RC__)

!     Add export specs for aerosols and aerosol deposition
!     ----------------------------------------------------
      call MAPL_AddExportSpec ( GC, SHORT_NAME = 'AERO',    &
                                CHILD_ID = AERO_PROVIDER, __RC__  )

#ifdef ENABLE_AERO_ACI
      ! IMPORTANT: This feature is currently disabled in Physics above,
      !            waiting for patches in the ESMF. Exporting state
      !            of child with diffefrent name is niot suypported.

      ! GOCART-2G uses a single AERO state for both Radiation and Moist,
      ! so we special handle it here. This approach should be adopted
      ! by the other components.
      ! ----------------------------------------------------------------
      if ( AERO_PROVIDER == GOCART2G ) then
         call MAPL_AddExportSpec ( GC,    TO_NAME = 'AERO_ACI', & 
                                       SHORT_NAME = 'AERO',     &
                                   CHILD_ID = AERO_PROVIDER, __RC__  )
      else
      call MAPL_AddExportSpec ( GC, SHORT_NAME = 'AERO_ACI',&
                                CHILD_ID = AERO_PROVIDER, __RC__  )
      endif
#endif

      call MAPL_AddExportSpec ( GC, SHORT_NAME = 'AERO_DP', &
                                CHILD_ID = AERO_PROVIDER, __RC__  )
  end if

  if (myState%enable_GOCART .and. chemReg%doing_CO2) then
      call MAPL_AddExportSpec ( GC, SHORT_NAME = 'CO2SC001', &
                                CHILD_ID = GOCART, __RC__ )
  end if

! Save this information in private state for later use. 
! WARNING: Dangerous if there is more than one RATS_PROVIDER
! ----------------------------------------------------------
  myState%RATS_PROVIDER = RATS_PROVIDER
  myState%AERO_PROVIDER = AERO_PROVIDER


  IF(MAPL_AM_I_ROOT()) THEN
   PRINT *," "
   PRINT *,TRIM(Iam)//":"
   PRINT *," AERO Provider is ", TRIM(providerName)
   PRINT *," "
  END IF

! Connectivities between Children
! -------------------------------
  IF(myState%enable_GOCART) then
     CALL MAPL_AddConnectivity ( GC, &
          SHORT_NAME  = (/'AIRDENS     ','AIRDENS_DRYP', 'DELP        ', 'CN_PRCP     ', 'NCN_PRCP    '/), &
          DST_ID = GOCART, SRC_ID = CHEMENV, __RC__  )
  ENDIF

  IF(myState%enable_GOCARTdata) then
     CALL MAPL_AddConnectivity ( GC, &
          SHORT_NAME  = (/'DELP       '/), &
          DST_ID = GOCARTdata, SRC_ID = CHEMENV, __RC__  )
  ENDIF

  IF(myState%enable_GOCART2G) then
     CALL MAPL_AddConnectivity ( GC, &
          SHORT_NAME  = (/'DELP    ', 'AIRDENS ', 'NCN_PRCP' /), &
          DST_ID = GOCART2G, SRC_ID = CHEMENV, __RC__  )
  ENDIF

  IF(myState%enable_QUICKCHEM) then
     CALL MAPL_AddConnectivity ( GC, &
          SHORT_NAME  = (/'DELP    ', 'AIRDENS ' /), &
          DST_ID = QUICKCHEM, SRC_ID = CHEMENV, __RC__  )
  ENDIF

  IF(myState%enable_GOCART .AND. myState%enable_QUICKCHEM) then

     IF(myState%running_CO ) THEN
       CALL MAPL_AddConnectivity(GC, SRC_NAME=(/'GOCART::CO'/),  DST_NAME=(/'CO'/),        SRC_ID=GOCART, DST_ID=QUICKCHEM, __RC__)
       CALL MAPL_AddConnectivity(GC, SRC_NAME=(/'GOCART::CO'/),  DST_NAME=(/'CO_avg24'/),  SRC_ID=GOCART, DST_ID=QUICKCHEM, __RC__)
     ENDIF
     IF(myState%running_CH4) THEN
       CALL MAPL_AddConnectivity(GC, SRC_NAME=(/'GOCART::CH4'/), DST_NAME=(/'CH4'/),       SRC_ID=GOCART, DST_ID=QUICKCHEM, __RC__)
       CALL MAPL_AddConnectivity(GC, SRC_NAME=(/'GOCART::CH4'/), DST_NAME=(/'CH4_avg24'/), SRC_ID=GOCART, DST_ID=QUICKCHEM, __RC__)
     ENDIF

     IF(myState%running_OH ) THEN
       CALL MAPL_AddConnectivity(GC, SRC_NAME=(/'OH'/),          DST_NAME=(/'CO_OH'/),     SRC_ID=QUICKCHEM, DST_ID=GOCART, __RC__)
       CALL MAPL_AddConnectivity(GC, SRC_NAME=(/'OH'/),          DST_NAME=(/'CH4_oh'/),    SRC_ID=QUICKCHEM, DST_ID=GOCART, __RC__)
     ENDIF
  ENDIF

  IF(myState%enable_GAAS) then
     CALL MAPL_AddConnectivity ( GC, &
          SHORT_NAME  = (/ 'AIRDENS ', 'DELP    ' /), &
          DST_ID = GAAS, SRC_ID = CHEMENV, __RC__  )
          IF(myState%enable_GOCART2G) then
             call MAPL_AddConnectivity ( GC,        &
                SHORT_NAME  = (/'AERO'/),       &
                DST_ID      =  GAAS,                &
                SRC_ID      =  GOCART2G, __RC__ )
          ELSE
              __raise__(MAPL_RC_ERROR,"Cannot have GAAS enabled without GOCART")
          ENDIF
  ENDIF

  IF(myState%enable_CARMA) then
      CALL MAPL_AddConnectivity ( GC, &
           SHORT_NAME  = (/'AIRDENS ', 'CN_PRCP ', 'NCN_PRCP'/), &
           DST_ID = CARMA, SRC_ID = CHEMENV, __RC__  )
           
      if(myState%enable_GOCART) then
       if(chemReg%doing_SU) then
         CALL MAPL_AddConnectivity ( GC, &
              SHORT_NAME  = (/'PSO4TOT'/), &
              DST_ID=CARMA, SRC_ID=GOCART, __RC__)
       endif
      endif 
  ENDIF

  IF(myState%enable_STRATCHEM) then
     CALL MAPL_AddConnectivity ( GC, &
          SHORT_NAME  = (/ 'AIRDENS     ', 'AIRDENS_DRYP', 'DELP        ',    &
                           'CN_PRCP     ', 'NCN_PRCP    ', 'LFR         ' /), &
          DST_ID = STRATCHEM, SRC_ID = CHEMENV, __RC__  )
  ENDIF

  IF(myState%enable_GMICHEM) then
     CALL MAPL_AddConnectivity ( GC, &
          SHORT_NAME  = (/'AIRDENS      ', 'DELP         ', 'LFR          ', &
                          'TPREC        ', 'LIGHT_NO_PROD'/),                &
          DST_ID = GMICHEM, SRC_ID = CHEMENV, __RC__  )
  ENDIF

  IF(myState%enable_MATRIX) then
     CALL MAPL_AddConnectivity ( GC, &
          SHORT_NAME  = (/'AIRDENS', 'DELP   '/), &
          DST_ID = MATRIX, SRC_ID = CHEMENV, __RC__  )
  ENDIF

  IF(myState%enable_MAM) then
     CALL MAPL_AddConnectivity ( GC, &
          SHORT_NAME  = (/'AIRDENS ', 'DELP    ', 'CN_PRCP ', 'NCN_PRCP'/), &
          DST_ID = MAM, SRC_ID = CHEMENV, __RC__  )
  ENDIF

  IF(myState%enable_ACHEM) then
     CALL MAPL_AddConnectivity ( GC, &
          SHORT_NAME  = (/'AIRDENS ', 'DELP    ', 'CN_PRCP ', 'NCN_PRCP'/), &
          DST_ID = ACHEM, SRC_ID = CHEMENV, __RC__  )
  ENDIF

  IF(myState%enable_TR) then
     CALL MAPL_AddConnectivity ( GC, &
          SHORT_NAME  = (/'AIRDENS     ', 'AIRDENS_DRYP', 'DELP        ', &
                          'CN_PRCP     ', 'NCN_PRCP    '/), &
          DST_ID = TR, SRC_ID = CHEMENV, __RC__  )
  ENDIF

  IF(myState%enable_TR .AND. myState%enable_GMICHEM) then
     ! First test - add O3 and the species needed to compute O3 loss
     ! Later, parse the TR .rc files to determine the fields we need
     CALL MAPL_AddConnectivity ( GC, &
          SHORT_NAME  = (/'OX    ', 'QQK007', 'QQK027', 'QQK028', 'DD_OX ', 'QQK005', &
                          'QQK235', 'QQK170', 'QQK216', 'QQK179', 'QQK150'/), &
          DST_ID = TR, SRC_ID = GMICHEM, __RC__  )
  ENDIF

  IF(myState%enable_GEOSCHEM) then
     CALL MAPL_AddConnectivity ( GC, &
          SHORT_NAME  = (/'AIRDENS', 'DELP   ', 'LFR    ', 'BYNCY  '/), &
          DST_ID = GEOSCHEM, SRC_ID = CHEMENV, __RC__  )
  ENDIF

! Ozone mole fraction needed by GOCART for
! CFC-12 photolysis.  For GMICHEM case, see below.
! ------------------------------------------------
  IF(myState%enable_GOCART .AND. myState%enable_PCHEM) then
   CALL MAPL_AddConnectivity ( GC, &
        SRC_NAME  = (/"OX"/), &
        DST_NAME  = (/"O3"/), &
        DST_ID = GOCART, SRC_ID = PCHEM, __RC__  )
  ENDIF

  IF(myState%enable_GOCART .AND. myState%enable_STRATCHEM) then
   CALL MAPL_AddConnectivity ( GC, &
        SRC_NAME  = (/"O3CHEM"/), &
        DST_NAME  = (/"O3"/), &
        DST_ID = GOCART, SRC_ID = STRATCHEM, __RC__  )
  ENDIF

  IF(myState%enable_ACHEM .AND. myState%enable_PCHEM) then
   CALL MAPL_AddConnectivity ( GC, &
        SRC_NAME  = (/"O3"/), &
        DST_NAME  = (/"O3"/), &
        DST_ID = ACHEM, SRC_ID = PCHEM, __RC__  )
  ENDIF

  IF(myState%enable_MAM .AND. myState%enable_ACHEM) then
   CALL MAPL_AddConnectivity ( GC, &
        SRC_NAME  = (/'ACHEM::SO2     ', 'ACHEM::H2SO4   ', 'ACHEM::NH3     ', 'ACHEM::SOAG    ', 'pSO4_aq        ', 'pNH4_aq        ',   &
                      'DDT_DMS_gas    ', 'DDT_MSA_gas    ', 'DDT_SO2_gas    ', 'DDT_H2SO4_gas  ', 'DDT_NH3_gas    ', 'DDT_SOAG_gas   ',   &
                      'DDT_DMS_aq     ', 'DDT_MSA_aq     ', 'DDT_SO2_aq     ', 'DDT_H2SO4_aq   ', 'DDT_NH3_aq     ', 'DDT_SOAG_aq    ',   &
                      '_DMS_gas       ', '_MSA_gas       ', '_SO2_gas       ', '_H2SO4_gas     ', '_NH3_gas       ', '_SOAG_gas      ',   &
                      '_DMS_aq        ', '_MSA_aq        ', '_SO2_aq        ', '_H2SO4_aq      ', '_NH3_aq        ', '_SOAG_aq       '/), &
        DST_NAME  = (/'SO2            ', 'H2SO4          ', 'NH3            ', 'SOA_GAS        ', 'pSO4_aq        ', 'pNH4_aq        ',   &
                      'DDT_DMS_gas    ', 'DDT_MSA_gas    ', 'DDT_SO2_gas    ', 'DDT_H2SO4_gas  ', 'DDT_NH3_gas    ', 'DDT_SOA_GAS_gas',   &
                      'DDT_DMS_aq     ', 'DDT_MSA_aq     ', 'DDT_SO2_aq     ', 'DDT_H2SO4_aq   ', 'DDT_NH3_aq     ', 'DDT_SOA_GAS_aq ',   &
                      '_DMS_gas       ', '_MSA_gas       ', '_SO2_gas       ', '_H2SO4_gas     ', '_NH3_gas       ', '_SOA_GAS_gas   ',   &
                      '_DMS_aq        ', '_MSA_aq        ', '_SO2_aq        ', '_H2SO4_aq      ', '_NH3_aq        ', '_SOA_GAS_aq    '/), &
        DST_ID = MAM, SRC_ID = ACHEM, __RC__  )
  ENDIF

  IF(myState%enable_MATRIX .AND. myState%enable_ACHEM) then
   CALL MAPL_AddConnectivity ( GC, &
        SRC_NAME  = (/'SO2    ', 'H2SO4  ', 'NH3    ', 'SOAG   ', 'pSO4_aq', 'pNH4_aq'/), &
        DST_NAME  = (/'SO2    ', 'H2SO4  ', 'NH3    ', 'SOA_GAS', 'pSO4_aq', 'pNH4_aq'/), &
        DST_ID = MATRIX, SRC_ID = ACHEM, __RC__  )
  ENDIF

  IF(myState%enable_STRATCHEM .AND. myState%enable_ACHEM) then
   CALL MAPL_AddConnectivity ( GC, &
        SRC_NAME  = (/'O3P      ', 'OH       ', 'O3CHEM   ', 'OCS_JRATE'/), &
        DST_NAME  = (/'O3P      ', 'OHSTRAT  ', 'O3       ', 'OCS_JRATE'/), &
        DST_ID = ACHEM, SRC_ID = STRATCHEM, __RC__  )
  ENDIF
 

! GOCART2G <=> ACHEM (OCS CHEMISTRY)
! ---------------------------------
  IF(myState%enable_GOCART2G .AND. myState%enable_ACHEM) then
   IF(chemReg%doing_OCS) THEN
    CALL MAPL_AddConnectivity ( GC, &
         SHORT_NAME  = (/'pSO2_OCS'/), &
         DST_ID = GOCART, SRC_ID = ACHEM, __RC__  )
   ENDIF
   CALL MAPL_AddConnectivity ( GC, &
        SHORT_NAME  = (/'pSOA_ANTHRO_VOC', 'pSOA_BIOB_VOC  '/), &
        DST_ID = GOCART2G, SRC_ID = ACHEM, __RC__  )
  ENDIF

! GOCART2G <=> QUICKCHEM
! ----------------------
  IF (myState%enable_GOCART2G .AND. myState%enable_QUICKCHEM) THEN

     IF (myState%running_BC ) THEN
       CALL MAPL_AddConnectivity ( GC,   SRC_NAME=(/'CASCACOEFCA.bc'/), DST_NAME=(/'BCSCACOEF'/),       SRC_ID=GOCART2G, DST_ID=QUICKCHEM, __RC__ )
       CALL MAPL_AddConnectivity ( GC,   SRC_NAME=(/'CASCACOEFCA.bc'/), DST_NAME=(/'BCSCACOEF_avg24'/), SRC_ID=GOCART2G, DST_ID=QUICKCHEM, __RC__ )
       ! MAY NEED TO CHANGE SYNTAX in NEXT VERSION OF G2G:  SRC_NAME=(/'CA.bcSCACOEF'/)
     ENDIF

     IF (myState%running_OC ) THEN
       CALL MAPL_AddConnectivity ( GC,   SRC_NAME=(/'CASCACOEFCA.oc'/), DST_NAME=(/'OCSCACOEF'/),       SRC_ID=GOCART2G, DST_ID=QUICKCHEM, __RC__ )
       CALL MAPL_AddConnectivity ( GC,   SRC_NAME=(/'CASCACOEFCA.oc'/), DST_NAME=(/'OCSCACOEF_avg24'/), SRC_ID=GOCART2G, DST_ID=QUICKCHEM, __RC__ )
       ! MAY NEED TO CHANGE SYNTAX in NEXT VERSION OF G2G:  SRC_NAME=(/'CA.ocSCACOEF'/)
     ENDIF

     IF (myState%running_BR ) THEN
       CALL MAPL_AddConnectivity ( GC,   SRC_NAME=(/'CASCACOEFCA.br'/), DST_NAME=(/'BRSCACOEF'/),       SRC_ID=GOCART2G, DST_ID=QUICKCHEM, __RC__ )
       CALL MAPL_AddConnectivity ( GC,   SRC_NAME=(/'CASCACOEFCA.br'/), DST_NAME=(/'BRSCACOEF_avg24'/), SRC_ID=GOCART2G, DST_ID=QUICKCHEM, __RC__ )
       ! MAY NEED TO CHANGE SYNTAX in NEXT VERSION OF G2G:  SRC_NAME=(/'CA.brSCACOEF'/)
     ENDIF

     IF (myState%running_DU ) THEN
       CALL MAPL_AddConnectivity ( GC, SHORT_NAME=(/'DUSCACOEF'/),                                      SRC_ID=GOCART2G, DST_ID=QUICKCHEM, __RC__ )
       CALL MAPL_AddConnectivity ( GC,   SRC_NAME=(/'DUSCACOEF'/),      DST_NAME=(/'DUSCACOEF_avg24'/), SRC_ID=GOCART2G, DST_ID=QUICKCHEM, __RC__ )
     ENDIF

     IF (myState%running_SU ) THEN
       CALL MAPL_AddConnectivity ( GC, SHORT_NAME=(/'SUSCACOEF'/),                                      SRC_ID=GOCART2G, DST_ID=QUICKCHEM, __RC__ )
       CALL MAPL_AddConnectivity ( GC,   SRC_NAME=(/'SUSCACOEF'/),      DST_NAME=(/'SUSCACOEF_avg24'/), SRC_ID=GOCART2G, DST_ID=QUICKCHEM, __RC__ )
     ENDIF

     IF (myState%running_SS ) THEN
       CALL MAPL_AddConnectivity ( GC, SHORT_NAME=(/'SSSCACOEF'/),                                      SRC_ID=GOCART2G, DST_ID=QUICKCHEM, __RC__ )
       CALL MAPL_AddConnectivity ( GC,   SRC_NAME=(/'SSSCACOEF'/),      DST_NAME=(/'SSSCACOEF_avg24'/), SRC_ID=GOCART2G, DST_ID=QUICKCHEM, __RC__ )
     ENDIF

     IF (myState%running_NI ) THEN
       CALL MAPL_AddConnectivity ( GC, SHORT_NAME=(/'NISCACOEF'/),                                      SRC_ID=GOCART2G, DST_ID=QUICKCHEM, __RC__ )
       CALL MAPL_AddConnectivity ( GC,   SRC_NAME=(/'NISCACOEF'/),      DST_NAME=(/'NISCACOEF_avg24'/), SRC_ID=GOCART2G, DST_ID=QUICKCHEM, __RC__ )
     ENDIF

  ENDIF

! GOCART <=> StratChem coupling ...
! ---------------------------------
  IF(myState%enable_GOCART .AND. myState%enable_STRATCHEM) then

   IF(chemReg%doing_SU) THEN
     do n = chemReg%i_SU, chemReg%j_SU

        if(trim(chemReg%vname(n)) .eq. 'SO4') then
           CALL MAPL_AddConnectivity ( GC, &
            SHORT_NAME = (/ "GOCART::DMS", "GOCART::SO2", "GOCART::SO4", "GOCART::MSA" /), &
            DST_ID = STRATCHEM, SRC_ID = GOCART, __RC__)
           CALL MAPL_AddConnectivity ( GC, &
            SHORT_NAME  = (/"SO4SAREA"/), &
            DST_ID = STRATCHEM, SRC_ID = GOCART, __RC__)
        endif

        if(trim(chemReg%vname(n)) .eq. 'SO4v') then
           CALL MAPL_AddConnectivity ( GC, &
            SHORT_NAME = (/"GOCART::DMSv", "GOCART::SO2v", "GOCART::SO4v", "GOCART::MSAv" /), &
            DST_ID = STRATCHEM, SRC_ID = GOCART, __RC__)
           CALL MAPL_AddConnectivity ( GC, &
            SHORT_NAME  = (/"SO4SAREAvolc"/), &
            DST_ID = STRATCHEM, SRC_ID = GOCART, __RC__)
        endif

     enddo

   END IF

  END IF

! CARMA <=> StratChem coupling ...
! ---------------------------------
  IF(myState%enable_STRATCHEM .AND. TRIM(providerName) == "CARMA") then
   CALL MAPL_AddConnectivity ( GC, &
                   SRC_NAME  = (/ 'CARMA_SUSAREA ', 'CARMA_SUSAREAv' /), &
                   DST_NAME  = (/ 'SO4SAREA      ', 'SO4SAREAvolc  ' /), &
            DST_ID = STRATCHEM, SRC_ID = CARMA, __RC__)

  END IF

! GOCART.data <=> GMICHEM coupling ...
! ------------------------------------
  IF(myState%enable_GMICHEM .AND. TRIM(providerName) == "GOCART.data") THEN
   CALL MAPL_AddConnectivity ( GC, SHORT_NAME=(/ "AERO" /), DST_ID=GMICHEM, SRC_ID=GOCARTdata, __RC__)
  END IF

! GOCART <=> GMICHEM coupling ...
! -------------------------------
  IF(myState%enable_GMICHEM .AND. TRIM(providerName) == "GOCART") THEN

! GOCART connections to GMICHEM
! -----------------------------
   IF(chemReg%doing_BC) &
   CALL MAPL_AddConnectivity ( GC, &
     SHORT_NAME  = (/ "GOCART::BCphobic", "GOCART::BCphilic" /), &
     DST_ID=GMICHEM, SRC_ID=GOCART, __RC__)

   IF(chemReg%doing_DU) &
   CALL MAPL_AddConnectivity ( GC, &
     SHORT_NAME  = (/ "GOCART::du001", "GOCART::du002", "GOCART::du003", &
                      "GOCART::du004", "GOCART::du005" /), &
     DST_ID=GMICHEM, SRC_ID=GOCART, __RC__)

   IF(chemReg%doing_OC) &
   CALL MAPL_AddConnectivity ( GC, &
     SHORT_NAME  = (/ "GOCART::OCphobic", "GOCART::OCphilic" /), &
     DST_ID=GMICHEM, SRC_ID=GOCART, __RC__)

   IF(chemReg%doing_SS) &
   CALL MAPL_AddConnectivity ( GC, &
     SHORT_NAME  = (/ "GOCART::ss001", "GOCART::ss002", "GOCART::ss003", &
                      "GOCART::ss004", "GOCART::ss005" /), &
     DST_ID=GMICHEM, SRC_ID=GOCART, __RC__)

   IF(chemReg%doing_SU) &
   CALL MAPL_AddConnectivity ( GC, &
     SHORT_NAME  = (/ "GOCART::DMS", "GOCART::SO2", &
                      "GOCART::SO4", "GOCART::MSA" /), &
     DST_ID=GMICHEM, SRC_ID=GOCART, __RC__)

     Search_SO4v: DO n = chemReg%i_SU, chemReg%j_SU

      IF( TRIM(chemReg%vname(n)) == "SO4v" )  THEN
        CALL MAPL_AddConnectivity ( GC, &
          SHORT_NAME  = (/ "GOCART::SO4v" /), &
          DST_ID=GMICHEM, SRC_ID=GOCART, __RC__)
      END IF

     END DO Search_SO4v


! GMICHEM connections to GOCART
! -----------------------------
  IF(myState%enable_GOCART) THEN

! ... For GOCART::SU,NI
! ---------------------
   CALL MAPL_AddConnectivity ( GC, &
    SHORT_NAME  = (/ "OH  ", "CH4 ", "H2O2", "NO3 ", "HNO3" /), &
    DST_ID=GOCART, SRC_ID=GMICHEM, __RC__)

! ... For GOCART::CFC
! -------------------
   CALL MAPL_AddConnectivity ( GC, &
    SRC_NAME  = (/"OX"/), &
    DST_NAME  = (/"O3"/), &
    DST_ID=GOCART, SRC_ID=GMICHEM, __RC__)
  
  END IF

! GMICHEM connections to MAM
! -----------------------------
   IF(myState%enable_MAM) &
   CALL MAPL_AddConnectivity ( GC, &
    SHORT_NAME  = (/ "ALK4", "CH4 ", "C2H6", "C3H6", "C3H8", "ISOP" /), &
    DST_ID=MAM, SRC_ID=GMICHEM, __RC__)

! GMICHEM connections to ACHEM
! --------------------------------
   IF(myState%enable_ACHEM) THEN

   CALL MAPL_AddConnectivity ( GC, &
    SHORT_NAME  = (/ "OH  ", "H2O2", "NO3 " /), &
    DST_ID=ACHEM, SRC_ID=GMICHEM, __RC__)

   CALL MAPL_AddConnectivity ( GC, &
    SRC_NAME  = (/"OX"/), &
    DST_NAME  = (/"O3"/), &
    DST_ID=ACHEM, SRC_ID=GMICHEM, __RC__)

   END IF

  END IF


! PCHEM connections to GEOSCHEM
! -----------------------------
  IF(myState%enable_PCHEM .AND. myState%enable_GEOSCHEM) THEN

   CALL MAPL_AddConnectivity ( GC, &
    SHORT_NAME  = (/ "TO3" /), &
    DST_ID=GEOSCHEM, SRC_ID=PCHEM, __RC__)

   ! added by ckeller, 10/31/2018
   CALL MAPL_AddConnectivity ( GC, &
    SRC_NAME  = (/"O3"/), &
    DST_NAME  = (/"PCHEM_O3"/), &
    DST_ID = GEOSCHEM, SRC_ID = PCHEM, __RC__  )

  END IF

  ! GOCART needs ozone for CFC12 photolysis.
  ! Only provide it from GEOS-Chem if PCHEM is not running.
  IF(myState%enable_GOCART .AND. myState%enable_GEOSCHEM) then
   IF ( .NOT. myState%enable_PCHEM ) THEN
    CALL MAPL_AddConnectivity ( GC, &
       SRC_NAME  = (/"TRC_O3"/), &
       DST_NAME  = (/"O3"/), &
       DST_ID=GOCART, SRC_ID=GEOSCHEM, __RC__  )
   ENDIF
  ENDIF

! HEMCO connections to CHEMENV
! -----------------------------
  ! Default values:
       GMI_instance_of_HEMCO = .FALSE.
    GOCART_instance_of_HEMCO = .FALSE.

  IF( myState%enable_HEMCO ) THEN

    CALL MAPL_AddConnectivity ( GC, SHORT_NAME  = (/ 'AIRDENS', 'BYNCY  ' /), DST_ID=HEMCO, SRC_ID=CHEMENV, __RC__)

    !!!!!!!!!!!!!!!!!!
    ! Determine if GOCART or GMI expect data from HEMCO:  
    ! (adapted from HEMCO_GridCompMod.F90)
    !!!!!!!!!!!!!!!!!!

    ! Define ESMF config for HEMCO
    HemcoCF = ESMF_ConfigCreate(__RC__)
    CALL ESMF_ConfigLoadFile( HemcoCF, 'HEMCO_GridComp.rc', __RC__ )

    ! Get number of instances
    call ESMF_ConfigGetAttribute(HemcoCF, nnInst, Label="HEMCO_Instances:" , DEFAULT=1, __RC__)

    ! Verbose
    IF ( MAPL_Am_I_Root() ) WRITE(*,*) 'CHEMsetup - number of HEMCO instances: ', nnInst


    ! Identify HEMCO instances
    DO N = 1, nnInst

       ! Get HEMCO configuration file names
       WRITE(Label,'(a14,i3.3,a1)') 'HEMCO_CONFIG--',N,':'
       call ESMF_ConfigGetAttribute( HemcoCF, ConfigFile, Label=TRIM(Label), &
                                     DEFAULT="HEMCOsa_Config.rc", __RC__)

       IF ( TRIM(ConfigFile) ==    'HEMCOgmi_Config.rc' )    GMI_instance_of_HEMCO = .TRUE.
       IF ( TRIM(ConfigFile) == 'HEMCOgocart2g_Config.rc' )  GOCART_instance_of_HEMCO = .TRUE.

       ! Verbose
       IF ( MAPL_Am_I_Root() ) WRITE(*,'(a19,i3.3,a2,a)') '--> HEMCO instance ', N, ': ', TRIM(ConfigFile)

    ENDDO

    IF ( MAPL_Am_I_Root() ) WRITE(*,*) '-->    GMI_instance_of_HEMCO: ',    GMI_instance_of_HEMCO
    IF ( MAPL_Am_I_Root() ) WRITE(*,*) '--> GOCART_instance_of_HEMCO: ', GOCART_instance_of_HEMCO
  END IF

! HEMCO -> GOCART
! ---------------
  IF( myState%enable_HEMCO .AND. myState%enable_GOCART2G ) THEN
   CALL MAPL_AddConnectivity ( GC, &
    SHORT_NAME  = (/ 'OC_ISOPRENE', 'OC_MTPA    ', 'OC_MTPO    ', 'OC_LIMO    '/), &
    SRC_ID=HEMCO, DST_ID=GOCART2G, __RC__)
  END IF

! HEMCO -> GMI 
! ------------
  IF( myState%enable_HEMCO .AND. myState%enable_GMICHEM ) THEN

     ! MEM - NOTE: We have more flags than we need;
     ! if we are running a GMI instance of HEMCO then it stands to reason that
     ! doMEGANviaHEMCO and doMEGANemission are both TRUE.
     ! Enforce consistency!
     ! (We could just AddConnectivity for GMI_ISOPRENE without checking
     ! whether it will be used, but if settings are inconsistent
     ! the model may not be doing what the user expects.)

     ! read GMI configuration
     gmi_config = ESMF_ConfigCreate(__RC__ )
     call ESMF_ConfigLoadFile(gmi_config, TRIM(gmi_rcfilen), __RC__ )

     call ESMF_ConfigGetAttribute(gmi_config, value=doMEGANemission, &
                                             label="doMEGANemission:", default=.FALSE., __RC__ )

     call ESMF_ConfigGetAttribute(gmi_config, value=doMEGANviaHEMCO, &
                                             label="doMEGANviaHEMCO:", default=.FALSE., __RC__ )
       
     ! make sure we don't have inconsistent HEMCO flags
     IF ( GMI_instance_of_HEMCO .neqv. doMEGANviaHEMCO ) THEN
        PRINT*,'Inconsistency --- HEMCO GMI instance  = ', GMI_instance_of_HEMCO
        PRINT*,'              --- GMI:doMEGANviaHEMCO = ', doMEGANviaHEMCO
        STATUS=98
        VERIFY_(STATUS)
     END IF

     ! make sure we don't have inconsistent MEGAN flags
     IF ( doMEGANviaHEMCO .eqv. .TRUE.    .AND.  &
          doMEGANemission .eqv. .FALSE. ) THEN
        PRINT*,'Inconsistent GMI flags: doMEGANviaHEMCO==T, doMEGANemission==F'
        STATUS=99
        VERIFY_(STATUS)
     END IF

     ! connect HEMCO isoprene to GMI if doMEGANviaHEMCO is true
     IF (doMEGANviaHEMCO) THEN
        CALL MAPL_AddConnectivity ( GC,          &
             SHORT_NAME  = (/'GMI_ISOPRENE'/),   &
             SRC_ID=HEMCO, DST_ID=GMICHEM, __RC__)
     ENDIF

  END IF

! Finally, set the services
! -------------------------
  call MAPL_GenericSetServices ( GC, __RC__ )

  RETURN_(ESMF_SUCCESS)
  
  end subroutine SetServices


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!BOP
! !IROUTINE: Initialize -- Initialized method for composite Aero-Chemistry

! !INTERFACE:

  subroutine Init ( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:

  type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
  type(ESMF_State),    intent(inout) :: IMPORT ! Import state
  type(ESMF_State),    intent(inout) :: EXPORT ! Export state
  type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
  integer, optional,   intent(  out) :: RC     ! Error code

! !DESCRIPTION: The Initialize method of the Chemistry Composite Gridded 
!  Component. It acts as a driver for the initializtion of the children.

!EOP

! ErrLog Variables

  __Iam__('Init')
  character(len=ESMF_MAXSTR)           :: COMP_NAME

! Local derived type aliases

   type (MAPL_MetaComp),       pointer :: MAPL
   type (ESMF_GridComp),       pointer :: GCS(:)
   type (ESMF_State)                   :: INTERNAL
   type (ESMF_State),          pointer :: GEX(:)
   type (ESMF_FieldBundle)             :: fBUNDLE
   type (ESMF_State)                   :: AERO
   type (ESMF_Config)                  :: CF, myCF

!   Private state
!   -------------
    type (GEOS_ChemGridComp), pointer  :: myState   ! private, that is
    type (GEOS_ChemGridComp_Wrap)      :: wrap

!=============================================================================
 
! Begin... 

!   Get the target components name and set-up traceback handle.
!   -----------------------------------------------------------
    call ESMF_GridCompGet ( GC, name=COMP_NAME, CONFIG=CF, __RC__ )
    Iam = trim(COMP_NAME) // "::" // trim(Iam)

!   Get my internal state
!   ---------------------
    call ESMF_UserCompGetInternalState(gc, 'GEOSchem_GridComp_State', WRAP, STATUS)
    VERIFY_(STATUS)
    myState => wrap%ptr

!   Call GenericInitialize for every Child
!   --------------------------------------
    call MAPL_GenericInitialize ( GC, IMPORT, EXPORT, CLOCK, __RC__ )

!   Get my MAPL_Generic state
!   --------------------------
    call MAPL_GetObjectFromGC ( GC, MAPL, __RC__ )

!   Get children and their im/ex states from my generic state.
!   ----------------------------------------------------------
    call MAPL_Get ( MAPL, GCS=GCS, GEX=GEX,              &
                    INTERNAL_ESMF_STATE=INTERNAL, __RC__ )

!   Fill in INTERNAL friendly bundle: CHEM_TRACERS
!   VERY IMPORTANT: Only the RATS provider can make OX friendly
!                   to ANALYSIS.
!   -----------------------------------------------------------
    call ESMF_StateGet (INTERNAL, 'CHEM_TRACERS', fBUNDLE, __RC__)
    call MAPL_GridCompGetFriendlies(GCS,                       &
                                         (/ "ANALYSIS    ",    &
                                            "TURBULENCE  ",    &
                                            "DYNAMICS    ",    &
                                            "MOIST       ",    &
                                            "STRATCHEM   ",    &
                                            "GMICHEM     ",    &
                                            "GEOSCHEMCHEM",    &
                                            "TR          " /), &
                                           fBUNDLE, AddGCPrefix=.true., __RC__ )




!   Get tracer increment species count from AGCM.rc
!   -----------------------------------------------
    call Initialize_IncBundle_init(GC, INTERNAL, EXPORT, CHMinc, __RC__)


!   AERO State for AERO_PROVIDER set to NONE
!   NOTE: This is an architecture violation. The parent should NEVER override
!   the child attributes. Whether these properties are provided or not need
!   need to be determined inside the child component.
!   --------------------------------------------------------------------------
    if (myState%AERO_PROVIDER < 0) then
        ! Radiation will not call the aerosol optics method 
        ! unless this attribute is explicitly set to true.
        call ESMF_StateGet(EXPORT, 'AERO', AERO, __RC__)
        call ESMF_AttributeSet(AERO, name='implements_aerosol_optics_method', value=.false., __RC__)

        ! Moist will not call the aerosol activation method 
        ! unless this attribute is explicitly set to true.
        call ESMF_AttributeSet(AERO, name='implements_aerosol_activation_properties_method', value=.false., __RC__)
    end if

#ifdef PRINT_STATES

!   Print what my states are
!   ------------------------
    if ( MAPL_am_I_root() ) then

       print *,  trim(Iam)//": IMPORT State" 
                                             call ESMF_StatePrint ( IMPORT )
       print *,  trim(Iam)//": INTERNAL State" 
                                             call ESMF_StatePrint ( INTERNAL )
       print *,  trim(Iam)//": EXPORT State" 
                                             call ESMF_StatePrint ( EXPORT )

       print *,  trim(Iam)//": AERO State (EXPORT)"
                                             call ESMF_StateGet   ( EXPORT, 'AERO', AERO, __RC__ ) 
                                             call ESMF_StatePrint ( AERO, nestedFlag=.false., __RC__ )
       print *,  trim(Iam)//": AERO State (PROVIDER)",  myState%AERO_PROVIDER
                                         if (myState%AERO_PROVIDER > 0) then
                                             call ESMF_StateGet   ( GEX(myState%AERO_PROVIDER), 'AERO', AERO, __RC__ )
                                             call ESMF_StatePrint ( AERO, nestedFlag=.false., __RC__ )
                                         end if

       print *,  trim(Iam)//": Friendly Tracers (INTERNAL)"
                                             call ESMF_FieldBundlePrint ( fBUNDLE )
    end if

#endif

!   All Done
!   --------
    RETURN_(ESMF_SUCCESS)

  end subroutine Init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !IROUTINE: Run1 -- phase 1 run method for the composite Physics Gridded Component

! !INTERFACE:

  subroutine Run1 ( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:

  type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
  type(ESMF_State),    intent(inout) :: IMPORT ! Import state
  type(ESMF_State),    intent(inout) :: EXPORT ! Export state
  type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
  integer, optional,   intent(  out) :: RC     ! Error code

! !DESCRIPTION: The Run 1 method of the Chemistry Composite Gridded Component.
!  It acts as a driver for the run phase 1 of the children. If a child has only
!  one run phase, it is skipped (only called in Run2).
!  (ckeller, 2014/09/10)
!EOP

! ErrLog Variables

  __Iam__('Run1')
  character(len=ESMF_MAXSTR)           :: COMP_NAME

! Local derived type aliases

  type (MAPL_MetaComp),     pointer  :: MAPL
  type (ESMF_Alarm)                  :: ALARM

  type (ESMF_GridComp),     pointer     :: GCS(:)
  type (ESMF_State),        pointer     :: GIM(:)
  type (ESMF_State),        pointer     :: GEX(:)
  type (ESMF_State)                     :: INTERNAL
  type (MAPL_MetaComp),     pointer     :: CHLD
  integer,                  allocatable :: CHLDPHASES(:)
  integer                               :: I, NCHLD
  integer                               :: NPHASE, IPHASE
  integer                               :: userRC
  character(len=ESMF_MAXSTR)            :: CHILD_NAME
  real, pointer                         :: th(:,:,:) => NULL()

! Private state
  type (GEOS_ChemGridComp), pointer  :: myState   ! private, that is
  type (GEOS_ChemGridComp_Wrap)      :: wrap

! For a 'gather' operation:
  REAL(KIND=DBL)    ::  local_tally
  REAL(KIND=DBL)    ::  global_tally
  TYPE(ESMF_VM)     ::  vm


!=============================================================================

! Begin... 

! Get parameters from generic state. The RUNALARM is used to control
!  the calling of the full chemistry
!-------------------------------------------------------------------

   call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_Get(MAPL, RUNALARM = ALARM, RC=STATUS )
   VERIFY_(STATUS)

!  Get my internal state
!  ---------------------
   call ESMF_UserCompGetInternalState(gc, 'GEOSchem_GridComp_State', WRAP, STATUS)
   VERIFY_(STATUS)
   myState => wrap%ptr

!  Start timers
!  ------------
   call MAPL_TimerOn( MAPL, "TOTAL")

   if ( ESMF_AlarmIsRinging (ALARM, RC=STATUS) ) then

      ! Don't turn alarm off in phase 1, otherwise phase 2 will not be 
      ! executed!
!      call ESMF_AlarmRingerOff(ALARM, RC=STATUS)
!      VERIFY_(STATUS)

      ! Get the target components name and set-up traceback handle.
      ! -----------------------------------------------------------
      call ESMF_GridCompGet ( GC, name=COMP_NAME, __RC__ )
      Iam = trim(COMP_NAME) // "::Run1"

!     Error trap: make sure import TH is filled
!     -----------------------------------------
      call MAPL_GetPointer ( IMPORT, th, 'TH', __RC__ )
      if( sum(th) == 0.0 ) then
         if ( MAPL_am_I_Root() ) then
            write(*,*) '***********************************************************************'
            write(*,*) TRIM(Iam), ': Chemistry import TH is all zero - skip run phase 1'
            write(*,*) '***********************************************************************'
         endif
         CALL MAPL_TimerOff( MAPL, "TOTAL" )
         RETURN_(ESMF_SUCCESS)
      endif

      ! Call Run phase 1 for every child with two phases
      ! (ckeller, 2014/09/10)
      ! ------------------------------------------------
      ! Get the children's states
      call MAPL_Get ( MAPL, GCS=GCS, GIM=GIM, GEX=GEX,   &
                    INTERNAL_ESMF_STATE=INTERNAL, __RC__ )

!    !Initialize chemistry increment bundle
!    !---------------------------------------
      call Initialize_IncBundle_run(INTERNAL, EXPORT, CHMinc, __RC__)


      if(associated(GCS)) then
        NCHLD  = SIZE(GCS)  ! # of children
        IPHASE = 1          ! phase to be called
        ! do for every child: get child state, determine number of 
        ! of run phases and call phase one if more than one phase
        ! exists. Also updated MAPL_Get to accept the output 
        ! argument NumRunPhases (ckeller, 09/10/2014)
        ! --------------------------------------------------------
        do I=1,NCHLD
          call ESMF_GridCompGet( GCS(I), NAME=CHILD_NAME, __RC__ )
          call MAPL_GetObjectFromGC(GCS(I), CHLD, __RC__ )
          call MAPL_Get(CHLD, NumRunPhases=NPHASE, __RC__ )
          if ( NPHASE > 1 ) then
            call MAPL_TimerOn(MAPL,trim(CHILD_NAME))
            call ESMF_GridCompRun (GCS(I), &
                     importState = GIM(I), &
                     exportState = GEX(I), &
                           clock = CLOCK,  &
                           phase = IPHASE, &
                          userRC = userRC, &
                                     __RC__ )
            _ASSERT(userRC==ESMF_SUCCESS,'needs informative message')

            IF ( myState%strict_child_timing ) THEN
              IF ( MAPL_Am_I_Root() ) PRINT*,'CHEM Run1: Gather after CHILD ',i
              local_tally = 123.0
              call ESMF_VmGetCurrent(vm, __RC__)
              call MAPL_CommsAllReduceSum(vm, sendbuf=local_tally, recvbuf=global_tally, cnt=1, __RC__)
            END IF

            call MAPL_TimerOff(MAPL,trim(CHILD_NAME))
          endif
        enddo !I
      endif

!   ! Compute Run1 chemistry increment
    !----------------------------------
      call Compute_IncBundle(INTERNAL, EXPORT, CHMinc, MAPL, __RC__)

   endif ! alarm is ringing

!  Stop timers
!  ------------
   call MAPL_TimerOff( MAPL, "TOTAL")

!   All Done
!   --------
    RETURN_(ESMF_SUCCESS)

 end subroutine Run1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !IROUTINE: Run2 -- Phase 2 run method for the composite Physics Gridded Component

! !INTERFACE:

  subroutine Run2 ( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:

  type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
  type(ESMF_State),    intent(inout) :: IMPORT ! Import state
  type(ESMF_State),    intent(inout) :: EXPORT ! Export state
  type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
  integer, optional,   intent(  out) :: RC     ! Error code

! !DESCRIPTION: The Run 2 method of the Chemistry Composite Gridded Component.
!  It acts as a driver for the run phase 2 of the children. This will execute
!  the phase 2 run routine of all children with two (or more) run phases, and
!  the phase 1 run routine of all children with only one run phase.
!  (ckeller, 2014/09/10)
!EOP

! ErrLog Variables

  __Iam__('Run2')
  character(len=ESMF_MAXSTR)           :: COMP_NAME

! Local derived type aliases

  type (MAPL_MetaComp),     pointer  :: MAPL
  type (ESMF_Alarm)                  :: ALARM
  type (ESMF_State)                  :: INTERNAL
!=============================================================================

    INTEGER, PARAMETER :: numRATs = 8
    CHARACTER(LEN=ESMF_MAXSTR) :: speciesName(numRats) = (/"OX    ", "O3    ", "CH4   ", "N2O   ", &
                                                           "CFC11 ", "CFC12 ", "HCFC22", "O3PPMV"/)
    integer :: i, im, jm, km, ijm
    real, pointer :: rat(:,:,:)
    real :: qmax, qmin

    type (ESMF_GridComp),     pointer     :: GCS(:)
    type (ESMF_State),        pointer     :: GIM(:)
    type (ESMF_State),        pointer     :: GEX(:)
    type (MAPL_MetaComp),     pointer     :: CHLD
    integer,                  allocatable :: CHLDPHASES(:)
    integer                               :: NCHLD
    integer                               :: NPHASE, IPHASE
    integer                               :: userRC
    character(len=ESMF_MAXSTR)            :: CHILD_NAME

!   Private state
    type (GEOS_ChemGridComp), pointer  :: myState   ! private, that is
    type (GEOS_ChemGridComp_Wrap)      :: wrap

!   For a 'gather' operation:
    REAL(KIND=DBL)    ::  local_tally
    REAL(KIND=DBL)    ::  global_tally
    TYPE(ESMF_VM)     ::  vm

!-------------------------------------------------------------------
! Begin... 

! Get parameters from generic state. The RUNALARM is used to control
!  the calling of the full chemistry
!-------------------------------------------------------------------

   call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_Get(MAPL, RUNALARM = ALARM, __RC__ ) 

!  Get my internal state
!  ---------------------
   call ESMF_UserCompGetInternalState(gc, 'GEOSchem_GridComp_State', WRAP, STATUS)
   VERIFY_(STATUS)
   myState => wrap%ptr

!  Start timers
!  ------------
   call MAPL_TimerOn( MAPL, "TOTAL")

   ! Turn off alarm
   ! --------------
   if ( ESMF_AlarmIsRinging   (ALARM, RC=STATUS) ) then
      call ESMF_AlarmRingerOff(ALARM, RC=STATUS)
      VERIFY_(STATUS)

      ! Get the target components name and set-up traceback handle.
      ! -----------------------------------------------------------
      call ESMF_GridCompGet ( GC, name=COMP_NAME, __RC__ )
      Iam = trim(COMP_NAME) // "::Run2" 

      ! Call Run for every Child. This is either phase 1 (for components
      ! with only one phase) or phase 2 (for components with two phases).
      ! (ckeller, 2014/09/10)
      ! -----------------------------------------------------------------
      ! Get the children's states

      call MAPL_Get ( MAPL, GCS=GCS, GIM=GIM, GEX=GEX,   &
                    INTERNAL_ESMF_STATE=INTERNAL, __RC__ )


      ! Initialize Chemistry increment bundle
      !----------------------------------------------------
      call Initialize_IncBundle_run(INTERNAL, EXPORT, CHMincR2, __RC__)


      if(associated(GCS)) then
        NCHLD  = SIZE(GCS)  ! # of children
        ! do for every child: get child state, determine number of 
        ! run phases and phase to call, execute.
        ! --------------------------------------------------------
        do I=1,NCHLD
          call MAPL_GetObjectFromGC(GCS(I), CHLD, __RC__ )
          call ESMF_GridCompGet( GCS(I), NAME=CHILD_NAME, __RC__ )
          call MAPL_Get(CHLD, NumRunPhases=NPHASE, __RC__ )
          if ( NPHASE > 1 ) then
            IPHASE = 2
          else
            IPHASE = 1
          endif

          call MAPL_TimerOn(MAPL,trim(CHILD_NAME))
          call ESMF_GridCompRun (GCS(I), &
                   importState = GIM(I), &
                   exportState = GEX(I), &
                         clock = CLOCK,  &
                         phase = IPHASE, &
                        userRC = userRC, &
                                   __RC__ )
          _ASSERT(userRC==ESMF_SUCCESS,'needs informative message')

          call MAPL_GenericRunCouplers( MAPL, i, CLOCK, __RC__ )

          IF ( myState%strict_child_timing ) THEN
            IF ( MAPL_Am_I_Root() ) PRINT*,'CHEM Run2: Gather after CHILD ',i
            local_tally = 123.0
            call ESMF_VmGetCurrent(vm, __RC__)
            call MAPL_CommsAllReduceSum(vm, sendbuf=local_tally, recvbuf=global_tally, cnt=1, __RC__)
          END IF

          call MAPL_TimerOff(MAPL,trim(CHILD_NAME))
        enddo !I
      endif

      !Compute chemistry tracer increments 
      !---------------------------------------
      call Compute_IncBundle(INTERNAL, EXPORT, CHMincR2, MAPL, __RC__)


!   Check contents of RATS
!   ----------------------
#ifdef DEBUG
    if ( MAPL_AM_I_ROOT() ) then
       print *, '--------------------------------------------------------'
    end if
    do i = 1, numRATS
       rat => null()
       call MAPL_GetPointer ( EXPORT, rat, trim(speciesName(i)), __RC__)
       if ( .not. associated(rat) ) then
          print *, '******** RAT with null pointer: ', trim(speciesName(i))
       else
          im = ubound(rat,1)
          jm = ubound(rat,2)
          km = ubound(rat,3)
          ijm = im*jm
          call pmaxmin(speciesName(i), rat(:,:,:), qmin, qmax, ijm, km, 1. )
       end if
    end do
#endif

   endif ! alarm is ringing

!  Stop timers
!  ------------
   call MAPL_TimerOff( MAPL, "TOTAL")

!   All Done
!   --------
    RETURN_(ESMF_SUCCESS)

 end subroutine Run2

!-----------------------------------------------------------------------

     subroutine GetProvider_ ( CF, Label, ID, Name, Default, RC )
!
!    Returns provider name as per resource file.
!
     type (ESMF_Config), intent(inout) :: CF
     character(len=*),   intent(in)    :: Label
     character(len=*),   intent(in)    :: Default
     character(len=*),   intent(out)   :: Name
     integer, intent(out)              :: ID
     integer, intent(out)              :: RC

! ErrLog Variables

     __Iam__('GetProvider_')

     character(len=ESMF_MAXSTR)        :: message
     character(len=ESMF_MAXSTR)        :: providerName


     call ESMF_ConfigGetAttribute(CF, providerName, Default=trim(default), &
                                  Label=Label, __RC__)

     ID   = -1 
     name = trim(providerName)

     select case ( trim(name) )

           case ('PCHEM')
                                    ID = PCHEM
           case ('GOCART.data')
                                    ID = GOCARTdata
           case ('GOCART')
                                    ID = GOCART
           case ('GOCART2G')
                                    ID = GOCART2G
           case ('QUICKCHEM')
                                    ID = QUICKCHEM
           case ('GAAS')
                                    ID = GAAS
           case ('H2O')
                                    ID = H2O
           case ('STRATCHEM')
                                    ID = STRATCHEM
           case ('GMICHEM')
                                    ID = GMICHEM
           case ('CARMA')
                                    ID = CARMA
           case ('GEOSCHEMCHEM')
                                    ID = GEOSCHEM
           case ('MATRIX')
                                    ID = MATRIX
           case ('MAM')
                                    ID = MAM
           case ('MAM.data')
                                    ID = MAMdata
           case ('TR')
                                    ID = TR
           case ('DNA')
                                    ID = DNA
           case DEFAULT

              message = "unknown provider "//trim(name)
              __raise__(MAPL_RC_ERROR, message)

     end select
     
     if ( ID < 0 ) then
        message = "Component "//trim(name)//" is NOT enabled; cannot specify it for "//trim(label)
        __raise__(MAPL_RC_ERROR, message)
     end if

     RC = ESMF_SUCCESS

   end subroutine GetProvider_

end module GEOS_ChemGridCompMod


