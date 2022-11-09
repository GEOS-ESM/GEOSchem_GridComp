#include "MAPL_Generic.h"
!-------------------------------------------------------------------------
!NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1, GEOS/DAS!
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  SC_GridCompMod --- SC Grid Component Class
!
! Grid Component class for NASA GSFC Code 614 stratospheric chemistry
!
! !INTERFACE:
!

   MODULE  SC_GridCompMod

! !USES:

   USE ESMF
   USE MAPL
   USE Runtime_RegistryMod
   USE Species_BundleMod
   USE Chem_Mod 	     ! Chemistry Base Class
   USE Chem_StateMod	     ! Chemistry State
   USE Chem_ConstMod, only: undefval => undef         ! Constants !
   USE Chem_UtilMod
   USE m_inpak90	     ! Resource file management
   USE WetRemovalMod
   USE ConvectionMod         ! Offline convective mixing/scavenging


   IMPLICIT none
   INTEGER, PARAMETER :: DBL = KIND(0.00D+00)


! !PUBLIC TYPES:
!
   PRIVATE
   PUBLIC  SC_GridComp       ! The SC object 

!
! !PUBLIC MEMBER FUNCTIONS:
!

   PUBLIC  SC_GridCompInitialize
   PUBLIC  SC_GridCompRun
   PUBLIC  SC_GridCompFinalize

   include "netcdf.inc"

!
! !DESCRIPTION:
!
!  This module implements the Stratospheric Chemistry Grid Component. 
!
! !REVISION HISTORY:
!
!  16Sep2003 da Silva  First crack.
!  24Jan2O05 Nielsen   Implementation of Code 916 chemistry
!  19Dec2005 da Silva  Minor portability mods.
!  11May2012 Nielsen   Capability for FV cubed [Ganymed-1_0_UNSTABLE]
!   3Jun2015 Liang     Added of halons, HCFCs, and 5 VSLSs. NUMPHOTO is now 52.
!  21Sep2016 Nielsen   Reduced equation set implemented. See notes below.
!
!  Usage notes for the REDUCED equation set:
!
!  Qing Liang's halon, HCFC, and VSLS updates were first implemented in tag Heracles-5_3_JEN_SC-v2, which
!  increased by 19 the number of transported species with respect to Heracles-5_3_JEN_SC and previous tags.
!  At the time, the GMAO was contemplating the use of SC instead of PCHEM in the DAS. But because the 
!  short-term impacts of the update on ozone were predicted to be small, the increase in wallclock, not less
!  than 15 to 20%, was considered an unacceptable burden on the DAS production cycle.
!
!  To resolve the issue for the DAS, while still retaining Qing's updates for chemistry-climate and other
!  applications, an optional directive -DREDUCED can be added to the USER_FDEFS in the StratChem GNUmakefile.
!  Upon compilation, -DREDUCED disables certain sections of the code added by Qing that are now wrapped
!  in #ifndef REDUCED / #endif blocks, which "reduces" the chemistry to that of tags that preceded 
!  Heracles-5_3_JEN_SC-v2 with two exceptions: (1) Aggregated HCFC is obsolete and (2) H12_24 is replaced
!  with H1211.
!
!  Disabling transport is accomplished by running the stratchem_setup script in src/Applications/GEOSgcm_App.
!  It reconfigures the Chem_Registry.rc in the experiment RC directory by assigning the 19 species to the 
!  non-transported (XX) variable list when it detects the presence of -DREDUCED in the GNUmakefile. The 
!  above is necessary since the chemistry GC relies on the registries to order the species in the bsc bundle.
!  
!  It is important to note that the import and export states, ExtData, and SC_GridComp.rc are not altered
!  by -DREDUCED. However, the content of most export states related to the 19 species, tendencies and reaction
!  rates for example, as well as that in their internal states, though still available, will not be valid.
!  Note also that by intent, -DREDUCED is considered the exceptional case and, therefore, is NOT the default.
!  It must be manually added to the GNUmakefile before compiling.
!
!EOP
!-------------------------------------------------------------------------

  TYPE SC_GridComp
    CHARACTER(LEN=ESMF_MAXSTR) :: name = "NASA_GSFC_Stratopheric_Chemistry"

! Dimensions
! ----------
   INTEGER :: i1, i2, im, j1, j2, jm, km

! Surface area [m^{2}] of grid cells
! ----------------------------------
   REAL(KIND=DBL), POINTER :: cellArea(:,:)

! Longitudes and latitudes (radians)
! ----------------------------------
   REAL, POINTER :: lonRad(:,:)
   REAL, POINTER :: latRad(:,:)

! Switches set in SC_GridComp.rc
! ------------------------------
    LOGICAL :: doPSCs
    LOGICAL :: doSediment
    LOGICAL :: importSulfSA   ! bring in GOCART, etc., based sulfate surface area
    LOGICAL :: addSAclim      ! if importSulfSA, possibly add climatological SA to it
    LOGICAL :: verbose
    LOGICAL :: wordy
    LOGICAL :: spinup
    LOGICAL :: useSolCyc
    LOGICAL :: doFlux

! Largest solar zenith angle (degrees) allowed as daytime.
! --------------------------------------------------------
    REAL :: szaCutoff

!Eccentricity of the orbiting body.
! --------------------------------------------------------
    REAL :: ecc_SC

! Layer number (1=surface)below which the special
! treatment of the NO photorate will NOT be applied.
! --------------------------------------------------
    REAL :: kNOspec

! Number of years to add or subtract (must be greater than -1000) from the 
! current year for the surface mixing ratios of greenhouse gases and ozone depleting  
! depleting substances. To lock-in a specified year, enter the negative of the year number.
! -----------------------------------------------------------------------------------------
    INTEGER :: GHGYrAdj
    INTEGER :: ODSYrAdj

! Option to prescribe the year for sulface surface area
! -----------------------------------------------------
    INTEGER :: SO4saYr

! Two calendar elements
! ---------------------
    LOGICAL :: isLeapYear
    REAL :: dayOfYear

! Uppermost layer to apply chemistry (1=surface).  Note: This is a *BAD*
! construct attributable to the legacy code that should be changed.
! ----------------------------------------------------------------------
    INTEGER :: levels_cal

! Thermal reaction rates and table lookup parameters.
! Dynamic allocation is now done in rdkrate.F.
! --------------------------------------------------
    INTEGER :: numreacs
    INTEGER :: niters = 3
    INTEGER :: num1d
    INTEGER :: num2d
    INTEGER :: num3d
    INTEGER, POINTER :: indxs1(:), indxs2(:), indxs3(:)
    INTEGER, POINTER :: KRxnOrder(:)
    REAL, POINTER :: cnsttab(:,:)
    LOGICAL, POINTER :: standardKRxn(:)
    CHARACTER(LEN=ESMF_MAXSTR), POINTER :: KRxnName(:)
    LOGICAL :: didKRxnsAlloc = .FALSE.

! Photolysis tables
! -----------------
    INTEGER :: numphoto = 55
    INTEGER :: nxdo
    INTEGER :: nlam
    INTEGER :: nsza
    INTEGER :: numo3
    INTEGER :: nts
    INTEGER :: aqsize

! ODS emission years
    INTEGER :: nyemiss

    REAL, POINTER :: sdat(:,:,:,:)
    REAL, POINTER :: o2jdat(:,:,:)
    REAL, POINTER :: sza_tab(:)
    REAL, POINTER :: o3_tab(:,:)
    REAL, POINTER :: xtab(:,:,:)
    REAL, POINTER :: CH2O_aq(:)
    REAL, POINTER :: rlam(:)

! Degrees-to-radians
! ------------------
    REAL :: dtr

! Mixing ratios of gases.  Allow changes to CO2 in FixSSG.
! --------------------------------------------------------
    REAL :: o2 = 0.20946
    REAL :: n2 = 0.7808
    REAL :: NOxTRW
    REAL :: parts = 10.00

! Polar stratospheric cloud parameters
! ------------------------------------
!    PSCtmax - Maximum temperature for calculation of sads
!    PSCpmax - Maximum pressure for calculation of sads. If pmax is 
!              negative, then pmax will be the tropopause pressure.
!    PSCpmin - Minimum pressure for calculation of sads
    REAL :: PSCtmax = 240.
    REAL :: PSCpmax = 150.  !SC_GridComp.rc value overrides in Initialize
    REAL :: PSCpmin = 10.
    REAL :: PSClatlim = 45. !SC_GridComp.rc value overrides in Initialize
                            !Changed to radians in Initialize

! Enforce a maximum allowable condensed HNO3 (ppbv).
! ---------------------------------------------------------------
   INTEGER :: HNO3Ice_MAX

! calcsts - flag to set type 1 aerosols to STS (true) or NAT (false)
! ------------------------------------------------------------------
    LOGICAL :: calcsts = .FALSE.

! ------------------------------------------------------------------
!   Parameters controlling ice PSCs:
!
!    nice - number of ice particles/cm^3
!    densice - max density of ice in grams/cm^3 condensed volume
!    sigice - width of lognormal particle size distribution
!    rice - median radius of ice particle size distribution (microns)
!           used if constantnice is false
!    satratice - saturation ratio required before formation of ice pscs
!    constantnice - if true, set n and sigma, let r float; if false,
!                   set r and sigma, let n float.
! ------------------------------------------------------------------
    LOGICAL :: constantnice=.TRUE.
    REAL :: nice = 1.e-2
    REAL :: sigice = 1.6
    REAL :: rice = 10.
    REAL :: densice = 1.
    REAL :: satratice = 1.

! Parameters controlling nat aerosols:
! Definitions for nat are analogous to those for ice
! --------------------------------------------------
    LOGICAL :: constantnnat=.TRUE.
    REAL :: nnat = 0.1
    REAL :: signat = 1.6
    REAL :: rnat = 0.4
    REAL :: densnat = 1.6
    REAL :: satratnat=1.0

! Parameters controlling sts aerosols:
! Definitions for sts are analogous to those for ice
! No need for saturation ratio because sts are liquid.
! ----------------------------------------------------
    LOGICAL :: constantnsts=.TRUE.
    REAL :: nsts = 10.0
    REAL :: sigsts = 1.6
    REAL :: rsts = 0.4

! Parameters defining assumptions made about lbs aerosols.  The
! subroutine sadsts calls the subroutine calcsulf, which uses these
! parameters and the liquid binary sulfate surface area density to
! estimate the h2so4 mixing ratio needed to calculate the sts surface
! area density
! -------------------------------------------------------------------
    LOGICAL :: constantnlbs=.TRUE.
    REAL :: nlbs = 10.
    REAL :: siglbs = 1.6
    REAL :: rlbs = 0.1

! Species placement in bsc%qa(ic)%data3d 
! --------------------------------------
    INTEGER ::    iAoA,      iCO2,    iSF6
    INTEGER :: 	   iOx,      iNOx,    iHNO3,     iN2O5
    INTEGER :: iHO2NO2,   iClONO2,     iClx,      iHCl
    INTEGER ::   iHOCl,     iH2O2,     iBrx,      iN2O
    INTEGER :: 	  iCl2,     iOClO,    iBrCl,      iHBr
    INTEGER :: iBrONO2,      iCH4,    iHOBr,   iCH3OOH
    INTEGER :: 	   iCO,    iHNO3c,      iO3,      iO3p   
    INTEGER ::    iF11,      iF12,    iF113,     iF114
    INTEGER ::   iF115 
    INTEGER ::   iCCl4,  iCH3CCl3,   iCH3Cl,    iCH3Br
    INTEGER ::  iH1301,    iH1211,   iH1202,    iH2402
    INTEGER :: iHCFC22, iHCFC141b, iHCFC142b,   iCHBr3
    INTEGER :: iCH2Br2,  iCH2BrCl,  iCHBr2Cl, iCHBrCl2
    INTEGER ::  iHFC23,    iHFC32,  iHFC125
    INTEGER ::iHFC134a,  iHFC143a,  iHFC152a
    INTEGER :: 	  iO1d,        iN,      iNO,      iNO2
    INTEGER :: 	  iNO3,        iH,      iOH,      iHO2
    INTEGER :: 	   iCl,      iClO,     iBrO,       iBr
    INTEGER ::  iCl2O2,     iCH2O,   iCH3O2,    irO3Ox

! Climatologies
! -------------
   INTEGER :: BCnymd
   CHARACTER(LEN=ESMF_MAXSTR) :: climFileName

   INTEGER :: numClimLats
   REAL, POINTER :: climLatRad(:)

   INTEGER :: numGHGs = 3
   INTEGER :: numSSGs = 25
   INTEGER :: numMBCs = 3 
!  CHARACTER(LEN=ESMF_MAXSTR) :: SSGName(23) = (/    "CO2     ", &
!    "CH4     ", "N2O     ", "CFC11   ", "CFC12   ", "CFC113  ", &
!    "CCL4    ", "CH3CCL3 ", "HCFC22  ", "HCFC141B", "HCFC142B", &
!    "H1211   ", "H1301   ", "H1202   ", "H2402   ", "CH3BR   ", &
!    "CH3CL   ", "HFC23   ", "HFC32   ", "HFC125  ", "HFC134A ", &
!    "HFC143A ", "HFC152A " /)

! CO2 is reset by SSGSfcFlx
! -------------------------
   REAL :: CO2 = 358.00E-06

! 2-D Climatological fields
! -------------------------
   REAL, POINTER :: SSG(:)
   REAL, POINTER :: climRain(:,:)
   REAL, POINTER :: climO3(:,:)
   REAL, POINTER :: climOH(:,:)
   REAL, POINTER :: climH2O(:,:)
   REAL, POINTER :: climSF6LOSS(:,:)
   REAL, POINTER :: SO4sa(:,:)
   REAL, POINTER :: CO2MBC(:)
   REAL, POINTER :: SF6MBC(:)

! 3-D interplolated climatological fields
! ---------------------------------------
   REAL, POINTER :: climRain3D(:,:,:)
   REAL, POINTER :: climO33D(:,:,:)
   REAL, POINTER :: climOH3D(:,:,:)
   REAL, POINTER :: climH2O3D(:,:,:)
   REAL, POINTER :: climSF6LOSS3D(:,:,:)
   REAL, POINTER :: SO4sa3D(:,:,:)

! Surface CO2/CH4/SF6 boundary conditions
   REAL, POINTER :: CO2MBC2D(:,:)
   REAL, POINTER :: SF6MBC2D(:,:)

! Surface ODS emissioss
   REAL, POINTER ::  CFC11EMISS(:,:)
   REAL, POINTER ::  CFC12EMISS(:,:)
   REAL, POINTER :: CFC113EMISS(:,:)
   REAL, POINTER ::    MCFEMISS(:,:)
   REAL, POINTER :: HCFC22EMISS(:,:)
   REAL, POINTER ::   CCL4EMISS(:,:)
   REAL, POINTER ::  CHBR3EMISS(:,:)
   REAL, POINTER :: CH2BR2EMISS(:,:)

   REAL :: emfCFC11
   REAL :: emfCFC12
   REAL :: emfCFC113
   REAL :: emfMCF
   REAL :: emfHCFC22
   REAL :: emfCCL4
   
  END TYPE SC_GridComp

CONTAINS

!-------------------------------------------------------------------------
!NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1, GEOS/DAS!
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  SC_GridCompInitialize --- Initialize SC_GridComp
!
! !INTERFACE:
!

   SUBROUTINE SC_GridCompInitialize ( gcSC, scReg, bsc, impChem, expChem, &
                                      nymd, nhms, tdt, rc )

! !USES:

  IMPLICIT none

! !INPUT PARAMETERS:

   TYPE(Runtime_Registry), POINTER :: scReg ! Names  of StratChem Species
   TYPE(Species_Bundle),   POINTER :: bsc   ! Bundle of StratChem Species
   INTEGER, INTENT(IN) :: nymd, nhms	    ! time
   REAL,    INTENT(IN) :: tdt		    ! chemistry time step (secs)


! !OUTPUT PARAMETERS:

   TYPE(SC_GridComp), INTENT(INOUT) :: gcSC	! Grid Component
   TYPE(ESMF_State), INTENT(INOUT)  :: impChem  ! Import State
   TYPE(ESMF_State), INTENT(INOUT)  :: expChem  ! Export State
   INTEGER, INTENT(OUT) ::  rc                  ! Error return code:
                                                !  0 - all is well
                                                !  1 - 

! !DESCRIPTION: Initializes the SC Grid Component. It primarily sets
!               the import state for each active constituent package.
!
! !REVISION HISTORY:
!
!  18Sep2003 da Silva  First crack.
!  24Jan2O05 Nielsen   Implementation of Code 916 chemistry
!  11May2012 Nielsen   Capability for FV cubed [Ganymed-1_0_UNSTABLE]
!
!EOP
!-------------------------------------------------------------------------
#include "mpif.h"

   CHARACTER(LEN=*), PARAMETER :: Iam = 'SC_GridCompInitialize'
   CHARACTER(LEN=ESMF_MAXSTR) :: rcfilen = 'SC_GridComp.rc'

   TYPE(ESMF_VM) :: vm

   INTEGER :: i, i1, i2, im, j, j1, j2, jm, k, km
   INTEGER :: n, status

   CHARACTER(LEN=ESMF_MAXSTR) :: fileName
   CHARACTER(LEN=ESMF_MAXSTR) :: fnKRates
   CHARACTER(LEN=ESMF_MAXSTR) :: fnPhoto

   gcSC%name = 'Stratospheric_Chemistry'
   IF(MAPL_AM_I_ROOT()) THEN
    PRINT *,' '
    PRINT *,TRIM(Iam)//':'
   END IF

! Initialize local variables
! NOTE: gcSC%lonRad(i,j) and gcSC%latRad(i,j) are in radians
! ----------------------------------------------------------
   rc = 0
   i1 = gcSC%i1
   i2 = gcSC%i2
   im = gcSC%im
   
   j1 = gcSC%j1
   j2 = gcSC%j2
   jm = gcSC%jm
   
   km = gcSC%km
   gcSC%levels_cal = km
   gcSC%BCnymd = -1

! Grab the virtual machine
! ------------------------
   CALL ESMF_VMGetCurrent(vm, RC=status)
   VERIFY_(status)

! Load resource file
! ------------------
   CALL I90_loadf ( TRIM(rcfilen), status )
   VERIFY_(status)

! Parse resource file
! -------------------
   CALL I90_label ( 'solar_ZA_cutoff:', status )
   VERIFY_(status)
   gcSC%szaCutoff = I90_gfloat( status )
   VERIFY_(status)

   CALL I90_label ( 'NO_cutoff:', status )
   VERIFY_(status)
   gcSC%kNOspec = I90_gfloat( status )
   VERIFY_(status)

   CALL I90_label ( 'TRW_NOx:', status )
   VERIFY_(status)
   gcSC%NOxTRW = I90_gfloat( status )
   VERIFY_(status)

   CALL I90_label ( 'Eccentricity_SC:', status )
   VERIFY_(status)
   gcSC%ecc_SC = I90_gfloat( status )
   VERIFY_(status)

   CALL I90_label ( 'PSCpmax:', status )
   VERIFY_(status)
   gcSC%PSCpmax = I90_gfloat( status )
   VERIFY_(status)

   CALL I90_label ( 'PSClatlim:', status )
   VERIFY_(status)
   gcSC%PSClatlim = I90_gfloat( status )
   VERIFY_(status)

   CALL I90_label ( 'HNO3Ice_MAX:', status )
   VERIFY_(status)
   gcSC%HNO3Ice_MAX = I90_gfloat( status )
   VERIFY_(status)

   CALL I90_label ( 'doPSCs:', status )
   VERIFY_(status)
   i = I90_gint( status )
   VERIFY_(status)
   IF(i == 0) THEN
    gcSC%doPSCs = .FALSE.
   ELSE
    gcSC%doPSCs = .TRUE.
   END IF

   CALL I90_label ( 'importSulfateSA:', status )
   VERIFY_(status)
   i = I90_gint( status )
   VERIFY_(status)
   IF(i == 0) THEN
    gcSC%importSulfSA = .FALSE.
   ELSE
    gcSC%importSulfSA = .TRUE.
   END IF

   CALL I90_label ( 'addSAclim:', status )
   VERIFY_(status)
   i = I90_gint( status )
   VERIFY_(status)
   IF(i == 0) THEN
    gcSC%addSAclim = .FALSE.
   ELSE
    gcSC%addSAclim = .TRUE.
   END IF

   CALL I90_label ( 'spinup:', status )
   VERIFY_(status)
   i = I90_gint( status )
   VERIFY_(status)
   IF(i == 0) THEN
    gcSC%spinup = .FALSE.
   ELSE
    gcSC%spinup = .TRUE.
   END IF

   CALL I90_label ( 'doSediment:', status )
   VERIFY_(status)
   i = I90_gint( status )
   VERIFY_(status)
   IF(i == 0) THEN
    gcSC%doSediment = .FALSE.
   ELSE
    gcSC%doSediment = .TRUE.
   END IF
   
   CALL I90_label ( 'verbose:', status )
   VERIFY_(status)
   i = I90_gint( status )
   VERIFY_(status)
   IF(i == 0) THEN
    gcSC%verbose = .FALSE.
   ELSE
    gcSC%verbose = .TRUE.
   END IF
   
   CALL I90_label ( 'useSolCyc:', status )
   VERIFY_(status)
   i = I90_gint( status )
   VERIFY_(status)
   IF(i == 0) THEN
    gcSC%useSolCyc = .FALSE.
   ELSE
    gcSC%useSolCyc = .TRUE.
   END IF

   CALL I90_label ( 'doFlux:', status )
   VERIFY_(status)
   i = I90_gint( status )
   VERIFY_(status)
   IF(i == 0) THEN
    gcSC%doFlux = .FALSE.
   ELSE
    gcSC%doFlux = .TRUE.
   END IF

   CALL I90_Label ( 'GHGYrAdj:', status )
   VERIFY_(status)
   gcSC%GHGYrAdj = I90_Gint( status )
   VERIFY_(status)

   CALL I90_Label ( 'ODSYrAdj:', status )
   VERIFY_(status)
   gcSC%ODSYrAdj = I90_Gint( status )
   VERIFY_(status)

   CALL I90_Label ( 'SO4saYr:', status )
   VERIFY_(status)
   gcSC%SO4saYr = I90_Gint( status )
   VERIFY_(status)

   CALL I90_label ( 'KRateTables:', status )
   VERIFY_(status)
   CALL I90_Gtoken( fnKRates, status )
   VERIFY_(status)

   CALL I90_label ( 'photolysisFile:', status )
   VERIFY_(status)
   CALL I90_Gtoken( fnPhoto, status )
   VERIFY_(status)

   CALL I90_label ( 'climatologiesFile:', status )
   VERIFY_(status)
   CALL I90_Gtoken( gcSC%climFileName, status )
   VERIFY_(status)

! Release resource file
! ---------------------
   CALL I90_release()

! Obtain species order in bsc%qa(ic)%data3d from the registry
! -----------------------------------------------------------
   DO i = 1,scReg%nq
    IF(TRIM(scReg%vname(i)) ==      "OX") gcSC%iOx = i
    IF(TRIM(scReg%vname(i)) ==     "NOX") gcSC%iNOx = i
    IF(TRIM(scReg%vname(i)) ==    "HNO3") gcSC%iHNO3 = i
    IF(TRIM(scReg%vname(i)) ==    "N2O5") gcSC%iN2O5 = i
    IF(TRIM(scReg%vname(i)) ==  "HO2NO2") gcSC%iHO2NO2 = i
    IF(TRIM(scReg%vname(i)) ==  "CLONO2") gcSC%iClONO2 = i
    IF(TRIM(scReg%vname(i)) ==     "CLX") gcSC%iClx = i
    IF(TRIM(scReg%vname(i)) ==     "HCL") gcSC%iHCl = i
    IF(TRIM(scReg%vname(i)) ==    "HOCL") gcSC%iHOCl = i
    IF(TRIM(scReg%vname(i)) ==    "H2O2") gcSC%iH2O2 = i
    IF(TRIM(scReg%vname(i)) ==     "BRX") gcSC%iBrx = i
    IF(TRIM(scReg%vname(i)) ==     "N2O") gcSC%iN2O = i
    IF(TRIM(scReg%vname(i)) ==     "CL2") gcSC%iCl2 = i
    IF(TRIM(scReg%vname(i)) ==    "OCLO") gcSC%iOClO = i
    IF(TRIM(scReg%vname(i)) ==    "BRCL") gcSC%iBrCl = i
    IF(TRIM(scReg%vname(i)) ==     "HBR") gcSC%iHBr = i
    IF(TRIM(scReg%vname(i)) ==  "BRONO2") gcSC%iBrONO2 = i
    IF(TRIM(scReg%vname(i)) ==     "CH4") gcSC%iCH4 = i
    IF(TRIM(scReg%vname(i)) ==    "HOBR") gcSC%iHOBr = i
    IF(TRIM(scReg%vname(i)) ==  "CH3OOH") gcSC%iCH3OOH = i
    IF(TRIM(scReg%vname(i)) ==      "CO") gcSC%iCO = i
    IF(TRIM(scReg%vname(i)) =="HNO3COND") gcSC%iHNO3c = i
    IF(TRIM(scReg%vname(i)) ==   "CFC11") gcSC%iF11 = i
    IF(TRIM(scReg%vname(i)) ==   "CFC12") gcSC%iF12 = i
    IF(TRIM(scReg%vname(i)) ==  "CFC113") gcSC%iF113 = i
    IF(TRIM(scReg%vname(i)) ==  "CFC114") gcSC%iF114 = i
    IF(TRIM(scReg%vname(i)) ==  "CFC115") gcSC%iF115 = i
    IF(TRIM(scReg%vname(i)) ==  "HCFC22") gcSC%iHCFC22 = i
    IF(TRIM(scReg%vname(i)) =="HCFC141B") gcSC%iHCFC141b = i
    IF(TRIM(scReg%vname(i)) =="HCFC142B") gcSC%iHCFC142b = i
    IF(TRIM(scReg%vname(i)) ==    "CCL4") gcSC%iCCl4 = i
    IF(TRIM(scReg%vname(i)) == "CH3CCL3") gcSC%iCH3CCl3 = i
    IF(TRIM(scReg%vname(i)) ==   "CH3CL") gcSC%iCH3Cl = i
    IF(TRIM(scReg%vname(i)) ==   "CH3BR") gcSC%iCH3Br = i
    IF(TRIM(scReg%vname(i)) ==   "H1301") gcSC%iH1301 = i
    IF(TRIM(scReg%vname(i)) ==   "H1211") gcSC%iH1211 = i
    IF(TRIM(scReg%vname(i)) ==   "H1202") gcSC%iH1202 = i
    IF(TRIM(scReg%vname(i)) ==   "H2402") gcSC%iH2402 = i
    IF(TRIM(scReg%vname(i)) ==   "CHBR3") gcSC%iCHBr3 = i
    IF(TRIM(scReg%vname(i)) ==  "CH2BR2") gcSC%iCH2Br2 = i
    IF(TRIM(scReg%vname(i)) == "CH2BRCL") gcSC%iCH2BrCl = i
    IF(TRIM(scReg%vname(i)) == "CHBR2CL") gcSC%iCHBr2Cl = i
    IF(TRIM(scReg%vname(i)) == "CHBRCL2") gcSC%iCHBrCl2 = i
    IF(TRIM(scReg%vname(i)) ==   "HFC23") gcSC%iHFC23 = i
    IF(TRIM(scReg%vname(i)) ==   "HFC32") gcSC%iHFC32 = i
    IF(TRIM(scReg%vname(i)) ==  "HFC125") gcSC%iHFC125 = i
    IF(TRIM(scReg%vname(i)) == "HFC134A") gcSC%iHFC134a = i
    IF(TRIM(scReg%vname(i)) == "HFC143A") gcSC%iHFC143a = i
    IF(TRIM(scReg%vname(i)) == "HFC152A") gcSC%iHFC152a = i
    IF(TRIM(scReg%vname(i)) ==     "CO2") gcSC%iCO2 = i
    IF(TRIM(scReg%vname(i)) ==     "SF6") gcSC%iSF6 = i
    IF(TRIM(scReg%vname(i)) == "AOADAYS") gcSC%iAoA = i
    IF(TRIM(scReg%vname(i)) ==  "O3CHEM") gcSC%iO3 = i
    IF(TRIM(scReg%vname(i)) ==     "O3P") gcSC%iO3p = i
    IF(TRIM(scReg%vname(i)) ==     "O1D") gcSC%iO1d = i
    IF(TRIM(scReg%vname(i)) ==       "N") gcSC%iN = i
    IF(TRIM(scReg%vname(i)) ==      "NO") gcSC%iNO = i
    IF(TRIM(scReg%vname(i)) ==     "NO2") gcSC%iNO2 = i
    IF(TRIM(scReg%vname(i)) ==     "NO3") gcSC%iNO3 = i
    IF(TRIM(scReg%vname(i)) == "HATOMIC") gcSC%iH = i
    IF(TRIM(scReg%vname(i)) ==      "OH") gcSC%iOH = i
    IF(TRIM(scReg%vname(i)) ==     "HO2") gcSC%iHO2 = i
    IF(TRIM(scReg%vname(i)) ==      "CL") gcSC%iCl = i
    IF(TRIM(scReg%vname(i)) ==     "CLO") gcSC%iClO = i
    IF(TRIM(scReg%vname(i)) ==     "BRO") gcSC%iBrO = i
    IF(TRIM(scReg%vname(i)) ==      "BR") gcSC%iBr = i
    IF(TRIM(scReg%vname(i)) ==   "CL2O2") gcSC%iCl2O2 = i
    IF(TRIM(scReg%vname(i)) ==    "CH2O") gcSC%iCH2O = i
    IF(TRIM(scReg%vname(i)) ==   "CH3O2") gcSC%iCH3O2 = i
    IF(TRIM(scReg%vname(i)) ==   "RO3OX") gcSC%irO3Ox = i
   END DO

! The expanded ODS and GHG reaction set is the default
! ----------------------------------------------------
   IF(MAPL_AM_I_ROOT()) THEN
#ifdef REDUCED
    PRINT *," "
    PRINT *,"W A R N I N G: REDUCED ODS and GHG reaction set is enabled"
#endif
   END IF

! PSC configuration summary
! -------------------------
   EchoInfo: IF(MAPL_AM_I_ROOT()) THEN

    PRINT *," "
    IF(gcSC%SO4saYr == 0 ) THEN
     n = nymd/10000
    ELSE
     n = gcSC%SO4saYr
    END IF

    IF(gcSC%doPSCs) THEN
     PRINT *,'Reactions for PSCs are enabled'
     IF(gcSC%calcsts) THEN
      PRINT *,' Type 1 aerosols are STS'
      PRINT *,' Sulfate surface area chosen for ',n
     ELSE
      PRINT *,' Type 1 aerosols are NAT'
      PRINT *,' Sulfate surface area ignored'
     END IF
     PRINT *," High-pressure range limit for PSCs: ",gcSC%PSCpmax," hPa"
     PRINT *," Low-latitude range limit for PSCs: ",gcSC%PSClatlim," degrees"
    ELSE
     PRINT *,TRIM(Iam)//': Reactions for PSCs disabled'
    END IF

    PRINT *," "
    PRINT *,"Daytime solar zenith angle cutoff: ",gcSC%szaCutoff," degrees"
    PRINT *,"Lightning NOx source strength: ",gcSC%NOxTRW," cm^{-3} s^{-1}"
    PRINT *," "
   END IF EchoInfo

! Useful constants
! ----------------
   gcSC%dtr = MAPL_PI/180.
   gcSC%PSClatlim = gcSC%PSClatlim*gcSC%dtr

! Thermal reaction rates: Initialize from ASCII file and print summary
! --------------------------------------------------------------------
   fileName = TRIM(fnKRates)
   CALL rdkrate(gcSC, fileName, MAPL_AM_I_ROOT(), status)
   VERIFY_(status)

! Photolysis tables: Initialize from NetCDF file
! ----------------------------------------------
   fileName = TRIM(fnPhoto)
   CALL readPhotTables(fileName, status)
   VERIFY_(status)

! Grab climatological array sizes for subsequent allocation
! ---------------------------------------------------------
   CALL initClimVars(status)
   VERIFY_(status)

! Allocations for climatological fields
! -------------------------------------
   ALLOCATE(gcSC%SSG(gcSC%numSSGs), STAT=status)
   VERIFY_(status)
   ALLOCATE(gcSC%climRain(gcSC%numClimLats,km), STAT=status)
   VERIFY_(status)
   ALLOCATE(gcSC%climO3(gcSC%numClimLats,km), STAT=status)
   VERIFY_(status)
   ALLOCATE(gcSC%climOH(gcSC%numClimLats,km), STAT=status)
   VERIFY_(status)
   ALLOCATE(gcSC%climH2O(gcSC%numClimLats,km), STAT=status)
   VERIFY_(status)
   ALLOCATE(gcSC%climSF6LOSS(gcSC%numClimLats,km), STAT=status)
   VERIFY_(status)
   ALLOCATE(gcSC%CO2MBC(gcSC%numClimLats), STAT=status)
   VERIFY_(status)
   ALLOCATE(gcSC%SF6MBC(gcSC%numClimLats), STAT=status)
   VERIFY_(status)

   ALLOCATE(gcSC%climRain3D(i1:i2,j1:j2,km), STAT=status)
   VERIFY_(status)
   ALLOCATE(gcSC%climO33D(i1:i2,j1:j2,km), STAT=status)
   VERIFY_(status)
   ALLOCATE(gcSC%climOH3D(i1:i2,j1:j2,km), STAT=status)
   VERIFY_(status)
   ALLOCATE(gcSC%climH2O3D(i1:i2,j1:j2,km), STAT=status)
   VERIFY_(status)
   ALLOCATE(gcSC%climSF6LOSS3D(i1:i2,j1:j2,km), STAT=status)
   VERIFY_(status)

   ALLOCATE(gcSC%SO4sa(gcSC%numClimLats,km), STAT=status)
   VERIFY_(status)
   ALLOCATE(gcSC%SO4sa3D(i1:i2,j1:j2,km), STAT=status)
   VERIFY_(status)

   ALLOCATE(gcSC%CO2MBC2D(i1:i2,j1:j2), STAT=status)
   VERIFY_(status)
   ALLOCATE(gcSC%SF6MBC2D(i1:i2,j1:j2), STAT=status)
   VERIFY_(status)

! Allocate space
! --------------
   ALLOCATE(gcSC%CFC11EMISS(i1:i2,j1:j2), STAT=status )
   VERIFY_(status)
   ALLOCATE(gcSC%CFC12EMISS(i1:i2,j1:j2), STAT=status )
   VERIFY_(status)
   ALLOCATE(gcSC%CFC113EMISS(i1:i2,j1:j2), STAT=status )
   VERIFY_(status)
   ALLOCATE(gcSC%MCFEMISS(i1:i2,j1:j2), STAT=status )
   VERIFY_(status)
   ALLOCATE(gcSC%HCFC22EMISS(i1:i2,j1:j2), STAT=status )
   VERIFY_(status)
   ALLOCATE(gcSC%CCL4EMISS(i1:i2,j1:j2), STAT=status )
   VERIFY_(status)
   ALLOCATE(gcSC%CHBR3EMISS(i1:i2,j1:j2), STAT=status )
   VERIFY_(status)
   ALLOCATE(gcSC%CH2BR2EMISS(i1:i2,j1:j2), STAT=status )
   VERIFY_(status)

!  As a safety check, where value is undefined set to 0
   gcSC%CFC11EMISS(i1:i2,j1:j2) = 0.0
   gcSC%CFC12EMISS(i1:i2,j1:j2) = 0.0
   gcSC%CFC113EMISS(i1:i2,j1:j2) = 0.0
   gcSC%MCFEMISS(i1:i2,j1:j2) = 0.0
   gcSC%HCFC22EMISS(i1:i2,j1:j2) = 0.0
   gcSC%CCL4EMISS(i1:i2,j1:j2) = 0.0
   gcSC%CHBR3EMISS(i1:i2,j1:j2)  = 0.0
   gcSC%CH2BR2EMISS(i1:i2,j1:j2) = 0.0

  RETURN

  CONTAINS

!-------------------------------------------------------------------------
!NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1, GEOS/DAS!
!-------------------------------------------------------------------------
!BOP
!
! !INTERFACE:

 SUBROUTINE readPhotTables(fileName, rc)

! !USES:

  IMPLICIT NONE

! !DESCRIPTION:
!
! Read tables for photolysis in StratChem ... from a NetCDF file
!
! Input parameters:
!
  CHARACTER(LEN=*), INTENT(IN) :: fileName
!
! Output parameters:
!
  INTEGER, INTENT(OUT) :: rc
!
! Restrictions:
!  ASSERT that the number of pressure layers in the dataset equals km.
!
! !REVISION HISTORY:
!  Nielsen     11 May 2012: First crack.
!
!-----------------------------------------------------------------------

  CHARACTER(LEN=ESMF_MAXSTR) :: Iam = "SC::readPhotTables"

  INTEGER :: comm, info, unit, status
  INTEGER :: dimid, i, n

  INTEGER :: length

  INTEGER, PARAMETER :: nD = 7
  CHARACTER(LEN=ESMF_MAXSTR) :: dimName(nD)= (/"nsza  ", "numO3 ", "layers", &
                                               "nlam  ", "nts   ", "nxdo  ", "aqsize" /)

  INTEGER, PARAMETER :: nV = 7
  CHARACTER(LEN=ESMF_MAXSTR) :: varName(nV)= (/"sza    ", &
                        "lambda ", "O3TAB  ",  "SDAT   ", &
                        "O2JDAT ", "XTAB   ",  "CH2O_AQ" /)
  rc = 0

  CALL ESMF_VMGet(vm, MPICOMMUNICATOR=comm, rc=status)
  VERIFY_(status)

#ifdef H5_HAVE_PARALLEL

  CALL MPI_Info_create(info, status)
  VERIFY_(status)
  CALL MPI_Info_set(info, "romio_cb_read", "automatic", status)
  VERIFY_(status)

#ifdef NETCDF_NEED_NF_MPIIO
  status = NF_OPEN_PAR(TRIM(fileName), IOR(NF_NOWRITE,NF_MPIIO), comm, info, unit)
#else
  status = NF_OPEN_PAR(TRIM(fileName), NF_NOWRITE, comm, info, unit)
#endif

#else

  IF(MAPL_AM_I_ROOT(vm)) THEN 
   status = NF_OPEN(TRIM(fileName), NF_NOWRITE, unit)

#endif

   IF(status /= NF_NOERR) THEN
    PRINT *,'Error opening file ',TRIM(fileName), status
    PRINT *, NF_STRERROR(status)
    VERIFY_(status)
   END IF

   DO i = 1,nD

    status = NF_INQ_DIMID(unit, TRIM(dimName(i)), dimid)
    IF(status /= NF_NOERR) THEN
     PRINT *,"Error inquiring dimension ID for ", TRIM(dimName(i)), status
     PRINT *, NF_STRERROR(status)
     VERIFY_(status)
    END IF

    status = NF_INQ_DIMLEN(unit, dimid, n)
    IF(status /= NF_NOERR) THEN
     PRINT *,"Error inquiring  dimension length for ", TRIM(dimName(i)), status
     PRINT *, NF_STRERROR(status)
    END IF

    SELECT CASE (i)
     CASE (1)
      gcSC%nsza = n
     CASE (2)
      gcSC%numO3 = n
     CASE (3)
      _ASSERT(n == km,'needs informative message')
     CASE (4)
      gcSC%nlam = n
     CASE (5)
      gcSC%nts = n
     CASE (6)
      gcSC%nxdo = n
     CASE (7)
      gcSC%aqsize = n
     CASE DEFAULT
    END SELECT

   END DO

#ifndef H5_HAVE_PARALLEL

  END IF ! MAPL_AM_I_ROOT

  CALL MAPL_CommsBcast(vm, gcSC%nsza, 1, 0, RC=status)
  VERIFY_(status)
  CALL MAPL_CommsBcast(vm, gcSC%numO3, 1, 0, RC=status)
  VERIFY_(status)
  CALL MAPL_CommsBcast(vm, gcSC%nlam, 1, 0, RC=status)
  VERIFY_(status)
  CALL MAPL_CommsBcast(vm, gcSC%nts, 1, 0, RC=status)
  VERIFY_(status)
  CALL MAPL_CommsBcast(vm, gcSC%nxdo, 1, 0, RC=status)
  VERIFY_(status)
  CALL MAPL_CommsBcast(vm, gcSC%aqSize, 1, 0, RC=status)
  VERIFY_(status)

#endif

  ALLOCATE(gcSC%sdat(gcSC%nsza,gcSC%numo3,km,gcSC%nlam), STAT=status)
  VERIFY_(status)
  ALLOCATE(gcSC%o2jdat(gcSC%nsza,gcSC%numo3,km), STAT=status)
  VERIFY_(status)
  ALLOCATE(gcSC%o3_tab(gcSC%numo3,km), STAT=status)
  VERIFY_(status)
  ALLOCATE(gcSC%xtab(gcSC%nlam,gcSC%nxdo,gcSC%nts), STAT=status)
  VERIFY_(status)
  ALLOCATE(gcSC%sza_tab(gcSC%nsza), STAT=status)
  VERIFY_(status)
  ALLOCATE(gcSC%CH2O_aq(gcSC%aqSize), STAT=status)
  VERIFY_(status)
  ALLOCATE(gcSC%rlam(gcSC%nlam), STAT=status)
  VERIFY_(status)

#ifndef H5_HAVE_PARALLEL

  IF(MAPL_AM_I_ROOT()) THEN

#endif

   DO i = 1,nV

    status = NF_INQ_VARID(unit, TRIM(varName(i)), n)
    IF(status /= NF_NOERR) THEN
     PRINT *,"Error getting varid for ", TRIM(varName(i)), status
     PRINT *, NF_STRERROR(status)
     VERIFY_(status)
    END IF

    SELECT CASE (i)
     CASE (1)
      status = NF_GET_VAR_REAL(unit, n, gcSC%sza_tab)
     CASE (2)
      status = NF_GET_VAR_REAL(unit, n, gcSC%rlam)
     CASE (3)
      status = NF_GET_VAR_REAL(unit, n, gcSC%o3_tab)
     CASE (4)
      status = NF_GET_VAR_REAL(unit, n, gcSC%sdat)
     CASE (5)
      status = NF_GET_VAR_REAL(unit, n, gcSC%o2jdat)
     CASE (6)
      status = NF_GET_VAR_REAL(unit, n, gcSC%xtab)
     CASE (7)
      status = NF_GET_VAR_REAL(unit, n, gcSC%CH2O_aq)
     CASE DEFAULT
    END SELECT

    IF(status /= NF_NOERR) THEN
     PRINT *,"Error getting values for ", TRIM(varName(i)), status
     PRINT *, NF_STRERROR(status)
     VERIFY_(status)
    END IF

   END DO

  status = NF_CLOSE(unit)
  VERIFY_(status)

#ifdef H5_HAVE_PARALLEL

   CALL MPI_Info_free(info, status)
   VERIFY_(status)

#else

  END IF ! MAPL_AM_I_ROOT

  length = SIZE(gcSC%sza_tab)
  CALL MPI_Bcast(gcSC%sza_tab, length, MPI_REAL, 0, comm, status)
  VERIFY_(status)

  length = SIZE(gcSC%rlam)
  CALL MPI_Bcast(gcSC%rlam, length, MPI_REAL, 0, comm, status)
  VERIFY_(status)

  length = SIZE(gcSC%o3_tab)
  CALL MPI_Bcast(gcSC%o3_tab, length, MPI_REAL, 0, comm, status)
  VERIFY_(status)

  length = SIZE(gcSC%sdat)
  CALL MPI_Bcast(gcSC%sdat, length, MPI_REAL, 0, comm, status)
  VERIFY_(status)

  length = SIZE(gcSC%o2jdat)
  CALL MPI_Bcast(gcSC%o2jdat, length, MPI_REAL, 0, comm, status)
  VERIFY_(status)

  length = SIZE(gcSC%xtab)
  CALL MPI_Bcast(gcSC%xtab, length, MPI_REAL, 0, comm, status)
  VERIFY_(status)

  CALL MAPL_CommsBcast(vm, gcSC%CH2O_aq, gcSC%aqsize, 0, RC=status)
  VERIFY_(status)

#endif

  RETURN
 END SUBROUTINE readPhotTables

!-------------------------------------------------------------------------
!NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1, GEOS/DAS!
!-------------------------------------------------------------------------
!BOP
!
! !INTERFACE:

 SUBROUTINE initClimVars(rc)

! !USES:

  IMPLICIT NONE

! !DESCRIPTION:
!
! Read and broadcast the latitudes (radians)
! on which the climatological fields are stored. 
!
! Input parameters:
!
! NONE
!
! Output parameters:
!
  INTEGER, INTENT(OUT) :: rc
!
! Restrictions:
!  ASSERT that the number of pressure layers in the dataset equals km.
!
! !REVISION HISTORY:
!  Nielsen     11 May 2012: First crack.
!
!-----------------------------------------------------------------------
 CHARACTER(LEN=ESMF_MAXSTR) :: Iam = "SC::readClimVars"

  INTEGER :: comm, dimid, info, n, status, unit
  CHARACTER(LEN=ESMF_MAXSTR) :: fileName

  rc = 0

  CALL ESMF_VMGet(vm, MPICOMMUNICATOR=comm, rc=status)
  VERIFY_(status)
  fileName = TRIM(gcSC%climFileName)

#ifdef H5_HAVE_PARALLEL

  CALL MPI_Info_create(info, status)
  VERIFY_(status)
  CALL MPI_Info_set(info, "romio_cb_read", "automatic", status)
  VERIFY_(status)

#ifdef NETCDF_NEED_NF_MPIIO
  status = NF_OPEN_PAR(TRIM(fileName), IOR(NF_NOWRITE,NF_MPIIO), comm, info, unit)
#else
  status = NF_OPEN_PAR(TRIM(fileName), NF_NOWRITE, comm, info, unit)
#endif

#else

  IF(MAPL_AM_I_ROOT(vm)) THEN 
   status = NF_OPEN(TRIM(fileName), NF_NOWRITE, unit)

#endif

   IF(status /= NF_NOERR) THEN
    PRINT *,"Error opening file ",TRIM(fileName), status
    PRINT *, NF_STRERROR(status)
    VERIFY_(status)
   END IF

   status = NF_INQ_DIMID(unit, 'lat', dimid)
   IF(status /= NF_NOERR) THEN
    PRINT *,"Error inquiring dimid for lat", status
    PRINT *, NF_STRERROR(status)
    VERIFY_(status)
   END IF
   status = NF_INQ_DIMLEN(unit, dimid, gcSC%numClimLats)
   IF(status /= NF_NOERR) THEN
    PRINT *,"Error inquiring dimlen of lat", status
    PRINT *, NF_STRERROR(status)
    VERIFY_(status)
   END IF

   status = NF_INQ_DIMID(unit, 'lev', dimid)
   IF(status /= NF_NOERR) THEN
    PRINT *,"Error inquiring dimid for lev", status
    PRINT *, NF_STRERROR(status)
    VERIFY_(status)
   END IF
   status = NF_INQ_DIMLEN(unit, dimid, n)
   IF(status /= NF_NOERR) THEN
    PRINT *,"Error inquiring dimlen of lev", status
    PRINT *, NF_STRERROR(status)
    VERIFY_(status)
   END IF

   _ASSERT(n == gcSC%km,'needs informative message')

#ifndef H5_HAVE_PARALLEL

  END IF ! MAPL_AM_I_ROOT

  CALL MAPL_CommsBcast(vm, gcSC%numClimLats, 1, 0, RC=status)
  VERIFY_(status)

#endif

  ALLOCATE(gcSC%climLatRad(gcSC%numClimLats), STAT=status)
  VERIFY_(status)

#ifndef H5_HAVE_PARALLEL

  IF(MAPL_AM_I_ROOT()) THEN

#endif

   status = NF_INQ_VARID(unit, 'lat', n)
   IF(status /= NF_NOERR) THEN
    PRINT *,"Error getting varid for lat", status
    PRINT *, NF_STRERROR(status)
    VERIFY_(status)
   END IF

   status = NF_GET_VAR_REAL(unit, n, gcSC%climLatRad)
   IF(status /= NF_NOERR) THEN
    PRINT *,"Error getting values for climLatRad", status
    PRINT *, NF_STRERROR(status)
    VERIFY_(status)
   END IF

  status = NF_CLOSE(unit)
  VERIFY_(status)

#ifdef H5_HAVE_PARALLEL

   CALL MPI_Info_free(info, status)
   VERIFY_(status)

#else

  END IF ! MAPL_AM_I_ROOT

  CALL MAPL_CommsBcast(vm, gcSC%climLatRad, SIZE(gcSC%climLatRad), 0, RC=status)
  VERIFY_(status)

#endif

  RETURN
 END SUBROUTINE initClimVars

  END SUBROUTINE SC_GridCompInitialize

!-------------------------------------------------------------------------
!NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1, GEOS/DAS!
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  SC_GridCompRun --- Stratospheric chemistry mechanism
!
! !INTERFACE:
!

   SUBROUTINE SC_GridCompRun ( gcSC, scReg, bsc, impChem, expChem, nymd, nhms, &
                               tdt, rc )

! !USES:

  IMPLICIT none

! !INPUT/OUTPUT PARAMETERS:

   TYPE(SC_GridComp), INTENT(inout) :: gcSC   ! Grid Component
   TYPE(Runtime_Registry), POINTER :: scReg ! Names  of StratChem Species
   TYPE(Species_Bundle),   POINTER :: bsc   ! Bundle of StratChem Species

! !INPUT PARAMETERS:

   TYPE(ESMF_State), INTENT(inout) :: impChem ! Import State
   INTEGER, INTENT(IN) :: nymd, nhms	      ! time
   REAL,    INTENT(IN) :: tdt		      ! chemical timestep (secs)

! !OUTPUT PARAMETERS:

   TYPE(ESMF_State), INTENT(inout) :: expChem   ! Export State
   INTEGER, INTENT(out) ::  rc                  ! Error return code:
                                                !  0 - all is well
                                                !  1 -
 
! !DESCRIPTION: This routine implements the so-called SC Driver. That 
!               is, adds chemical tendencies to each of the constituents
!
! !REVISION HISTORY:
!
!  18Sep2003 da Silva  First crack.
!  24Jan2O05 Nielsen   Implementation of Code 916 chemistry
!  11May2012 Nielsen   Capability for FV cubed [Ganymed-1_0_UNSTABLE]
!
!EOP
!-------------------------------------------------------------------------
#include "mpif.h"

   CHARACTER(LEN=*), PARAMETER :: Iam = 'SC_GridCompRun'
   TYPE(ESMF_VM) :: vm

! Input fields from fvGCM
! -----------------------
   REAL, POINTER, DIMENSION(:,:,:) ::  specHum, T, qctot, ple, dqcond
   REAL, POINTER, DIMENSION(:,:,:) ::  rhoa, rhoadry, hghte
   REAL, POINTER, DIMENSION(:,:)   ::  cellArea
   REAL, POINTER, DIMENSION(:,:)   ::  oro   ! Surface ocean-land-ice mask
   REAL, POINTER, DIMENSION(:,:)   ::  precc, precl, pblh 
   REAL, POINTER, DIMENSION(:,:,:) ::  pfllsan, pfilsan 
   REAL, POINTER, DIMENSION(:,:)   ::  frlake, frocean, frseaice

!  Additional needs for GOCART convective diagnostic
   real, pointer, dimension(:,:,:)       ::  cmfmc, qccu, dtrain
   real*8, allocatable, dimension(:,:,:) ::  cmfmc_, qccu_, dtrain_, &
                                             airmass_, airmol_, vud_, &
                                             delz_, delp_, h2o2_, tmpu_, ple_
   real*8, allocatable                   ::  tc_(:,:,:,:), bcnv_(:,:,:)
   real*8, allocatable                   ::  area_(:,:), frlake_(:,:), &
                                             frocean_(:,:), frseaice_(:,:)
   integer*4                             ::  icdt
   
! Input fields from ExtData
   REAL, POINTER, DIMENSION(:,:)   ::  f11ems, f12ems, f113ems
   REAL, POINTER, DIMENSION(:,:)   ::  mcfems, ccl4ems, hcfc22ems
   REAL, POINTER, DIMENSION(:,:)   ::  chbr3ems, ch2br2ems
   REAL, POINTER, DIMENSION(:,:,:) ::  oh3dClim

! Exports not part of internal state
! ----------------------------------
   REAL, POINTER, DIMENSION(:,:,:) :: natsad,   icesad,   o3ppmv,    ozone
   REAL, POINTER, DIMENSION(:,:,:) :: dOxdt,    dQdt
   REAL, POINTER, DIMENSION(:,:)   :: n2oflx,   ch4flx
   REAL, POINTER, DIMENSION(:,:)   :: f12flx,   f11flx,  f113flx
   REAL, POINTER, DIMENSION(:,:)   :: ccl4flx,  mcfflx, ch3clflx, ch3brflx
   REAL, POINTER, DIMENSION(:,:)   :: h1301flx, h1211flx, h1202flx, h2402flx
   REAL, POINTER, DIMENSION(:,:)   :: hcfc22flx,  hcfc141bflx,  hcfc142bflx 
   REAL, POINTER, DIMENSION(:,:)   :: szarad,  scbasek
!  REAL, POINTER, DIMENSION(:,:,:) :: SO4SAimp  ! imported sulfate surface area density m2 m-3
   REAL, POINTER, DIMENSION(:,:,:) :: SO4SAv,SO4SAt  ! imported sulfate surface area density m2 m-3
   REAL, POINTER, DIMENSION(:,:)   :: mcfocnloss, ctcocnloss, ctclndloss
   REAL, POINTER, DIMENSION(:,:,:) :: jocs    ! j-rates of OCS
   REAL, POINTER, DIMENSION(:,:,:) :: rdair   ! ratio of dry air vs. moist air
   REAL, POINTER, DIMENSION(:,:,:) :: qqj001, qqj002, qqj003, qqj004, qqj005
   REAL, POINTER, DIMENSION(:,:,:) :: qqj006, qqj007, qqj008, qqj009, qqj010
   REAL, POINTER, DIMENSION(:,:,:) :: qqj011, qqj012, qqj013, qqj014, qqj015
   REAL, POINTER, DIMENSION(:,:,:) :: qqj016, qqj017, qqj018, qqj019, qqj020
   REAL, POINTER, DIMENSION(:,:,:) :: qqk001, qqk002, qqk003, qqk004, qqk005
   REAL, POINTER, DIMENSION(:,:,:) :: qqk006, qqk007, qqk008, qqk009, qqk010
   REAL, POINTER, DIMENSION(:,:,:) :: qqk011, qqk012, qqk013, qqk014, qqk015
   REAL, POINTER, DIMENSION(:,:,:) :: qqk016, qqk017, qqk018, qqk019, qqk020
   REAL, POINTER, DIMENSION(:,:,:) :: qqk021, qqk022, qqk023, qqk024, qqk025
   REAL, POINTER, DIMENSION(:,:,:) :: qqk026, qqk027, qqk028, qqk029, qqk030
   REAL, POINTER, DIMENSION(:,:,:) :: qqk031, qqk032, qqk033, qqk034, qqk035
    
   
   type(Chem_Array), pointer :: fluxout

! Local
! -----
   INTEGER, PARAMETER :: DBL = KIND(0.00D+00)
   INTEGER, PARAMETER :: FROM_BUNDLE = 1, TO_BUNDLE = -1
   REAL, PARAMETER :: MIN_VALUE = 1.00E-25

   TYPE(ESMF_Clock) :: SCClock
   TYPE(ESMF_Time)  :: timeNow
   TYPE(ESMF_TimeInterval) :: interval
   INTEGER :: ccyy, mm, dd, hh, m, s

   INTEGER :: i, i1, i2, ic, iXj, im, ijkl, ijk1l
   INTEGER :: j, j1, j2, jm, k, km
   INTEGER :: status

   INTEGER :: examineDt
   INTEGER :: iterCount
   INTEGER :: k1Strat, kRev
   INTEGER :: passNumber, requiredPasses

   REAL :: dt, qmax, qmin, requiredDt, szan, tropp, tx, rdistsq

   LOGICAL :: dayTime
   LOGICAL :: startOver
   LOGICAL :: tropo(gcSC%km)
   LOGICAL :: KIN

! Profiles
! --------
   REAL :: cellDepth(gcSC%km), cellVolume(gcSC%km)
   REAL :: climRainProfile(gcSC%km), climO3Profile(gcSC%km)
   REAL :: climH2OProfile(gcSC%km), climOHProfile(gcSC%km)
   REAL :: climSF6LOSSProfile(gcSC%km)
   REAL :: denssts(gcSC%km), edgePress(0:gcSC%km)
   REAL :: gcrProfile(gcSC%km), H2OcProfile(gcSC%km), HNO3cProfile(gcSC%km)
   REAL :: kel(gcSC%km), midPress(gcSC%km)
   REAL :: NOColumn(gcSC%km), NOProfile(gcSC%km)
   REAL :: numDens(gcSC%km), numDensDry(gcSC%km)
   REAL :: O2Column(gcSC%km), O3Column(gcSC%km), O3Profile(gcSC%km)
   REAL :: rmedice(gcSC%km), rmednat(gcSC%km), rmedsts(gcSC%km)
   REAL :: SO4saProfile(gcSC%km), vfall(gcSC%km)
   REAL :: mcfprofile(gcSC%km),ctcprofile(gcSC%km)

! Reaction rates
! --------------
   REAL(KIND=DBL) :: aj(gcSC%numphoto), ak(gcSC%km,gcSC%numreacs)

! Output loss rates
! --------------
   REAL(KIND=DBL) :: qqj(20), qqk(35)

! Initial state holding arrays
! ----------------------------
   REAL(KIND=DBL), ALLOCATABLE :: speciesBase(:)
   REAL(KIND=DBL), ALLOCATABLE :: ratesBase(:)

! Species, etc.
! -------------
   REAL(KIND=DBL) ::     h2o,     h2oc,       ox,     nox,    hno3
   REAL(KIND=DBL) ::    n2o5,   ho2no2,   clono2,     clx,     hcl
   REAL(KIND=DBL) ::    hocl,     h2o2,      brx,     n2o,     cl2
   REAL(KIND=DBL) ::    oclo,     brcl,      hbr,  brono2,     ch4
   REAL(KIND=DBL) ::    hobr,   ch3ooh,       co,   hno3c,    ccl4
   REAL(KIND=DBL) ::     f11,      f12,     f113,    f114,    f115
   REAL(KIND=DBL) ::  hcfc22, hcfc141b, hcfc142b,     co2,    sf6
   REAL(KIND=DBL) ::   h1301,    h1211,    h1202,   h2402
   REAL(KIND=DBL) ::   chbr3,   ch2br2,  ch2brcl, chbr2cl, chbrcl2
   REAL(KIND=DBL) :: ch3ccl3,    ch3cl,    ch3br,   hfc23,   hfc32
   REAL(KIND=DBL) ::  hfc125,  hfc134a,  hfc143a, hfc152a
   REAL(KIND=DBL) ::      o3,      o3p,      o1d,       n,      no
   REAL(KIND=DBL) ::     no2,      no3,        h,      oh,     ho2
   REAL(KIND=DBL) ::      cl,      clo,      bro,      br,   cl2o2
   REAL(KIND=DBL) ::    ch2o,    ch3o2,      ch3,    ch3o,     cho
   REAL(KIND=DBL) ::   saice,    sanat,  box_ro3ox

! Estimated species
! -----------------
   REAL(KIND=DBL) ::    o3pe,     o3e,     oxe,      ne,     noe,    no2e
   REAL(KIND=DBL) ::    no3e,    noxe,   n2o5e,   h2o2e,   hno3e, ho2no2e
   REAL(KIND=DBL) ::    hcle,   hocle,     cle,    cloe,   ocloe,  cl2o2e
   REAL(KIND=DBL) :: clono2e,    clxe,     bre,    broe, brono2e,    brxe
   REAL(KIND=DBL) ::    hbre,   hobre

! Ratios
! ------
   REAL(KIND=DBL) ::    ro3po3,    ro1do3,     ro3ox,      rnno
   REAL(KIND=DBL) ::    rnono2,   rno3no2,   rno2nox,    rclclo
   REAL(KIND=DBL) ::   rcloclx,  rocloclo, rcl2o2clo,   rbrobrx
   REAL(KIND=DBL) ::   rbrnbrx

! Production and loss
! -------------------
   REAL(KIND=DBL) ::  	    po3,        lo3,      ph2o2,      lh2o2,      phno3
   REAL(KIND=DBL) ::      lhno3,   phno3het,      pn2o5,      ln2o5,   ln2o5het
   REAL(KIND=DBL) ::    pho2no2,    lho2no2,       pnox,       lnox,     lnoxsq
   REAL(KIND=DBL) ::     pnoxcl,     lnoxcl,       pclx,       lclx,    pclono2
   REAL(KIND=DBL) ::    lclono2, lclono2het,       phcl,       lhcl,      phocl
   REAL(KIND=DBL) ::      lhocl,    pbrono2, lbrono2het,       pbrx,       lbrx
   REAL(KIND=DBL) ::       phbr,       lhbr,      phobr,      lhobr,       pno2
   REAL(KIND=DBL) ::       lno2,       pno3,       lno3,      poclo,      pnoxa
   REAL(KIND=DBL) ::      lnoxa,       pcl2,      pbrcl,     lo3hox,     lo3nox
   REAL(KIND=DBL) ::    lo3oxsq,      lo3cl,     lo3cly,      lo3oh,     lo3brx
   REAL(KIND=DBL) ::     lo3chx

! Total column density
   REAL(KIND=DBL) ::   mcfColumn, ctcColumn

! For calls to the old GOCART wet removal routine:
   TYPE(Chem_Array), pointer     :: qa_single(:)
   INTEGER :: species_index




! Grid specs from Chem_Bundle%grid
! --------------------------------
  rc = 0
  i1 = gcSC%i1
  i2 = gcSC%i2
  im = gcSC%im

  j1 = gcSC%j1
  j2 = gcSC%j2
  jm = gcSC%jm

  km = gcSC%km

  iXj   = (i2-i1+1)*(j2-j1+1)
  ijkl  = iXj * km
  ijk1l = iXj * (km+1) 

! Grab the virtual machine
! ------------------------
  CALL ESMF_VMGetCurrent(vm, RC=status)
  VERIFY_(status)

! Grab some calender elements
! ---------------------------
  ccyy = nymd/10000
  mm = (nymd-ccyy*10000)/100
  dd = nymd-ccyy*10000-mm*100
  hh = nhms/10000
  m = (nhms-hh*10000)/100
  s = nhms-hh*10000-m*100



! Create and set the clock
! ------------------------
  CALL ESMF_TimeSet(timeNow, YY=ccyy, MM=mm, DD=dd, H=hh, M=m, S=s, RC=status)
  VERIFY_(status)
  CALL ESMF_TimeIntervalSet(interval, S=s, RC=status)
  VERIFY_(status)
  SCClock = ESMF_ClockCreate(NAME="SCClock", TIMESTEP=interval, STARTTIME=timeNow, RC=status)
  VERIFY_(status)
  CALL ESMF_ClockSet(SCClock, CURRTIME=timeNow, RC=status)
  VERIFY_(status)

! Import state pointers
! ---------------------
  CALL MAPL_GetPointer(impChem,  specHum,	'Q', RC=status)
  VERIFY_(status)
  CALL MAPL_GetPointer(impChem,    qctot,   'QCTOT', RC=status)
  VERIFY_(status)
  CALL MAPL_GetPointer(impChem,     rhoa, 'AIRDENS', RC=status)
  VERIFY_(status)
  CALL MAPL_GetPointer(impChem,  rhoadry, 'AIRDENS_DRYP', RC=status)
  VERIFY_(status)

!  Get 2D Imports
!  --------------
   call MAPL_GetPointer ( impChem, pblh,     'ZPBL',    rc=status )
   VERIFY_(status)
   call MAPL_GetPointer ( impChem, oro,      'LWI',     rc=status )
   VERIFY_(status)
   call MAPL_GetPointer ( impChem, precc,    'CN_PRCP', rc=status )
   VERIFY_(status)
   call MAPL_GetPointer ( impChem, precl,    'NCN_PRCP',rc=status )
   VERIFY_(status)
   call MAPL_GetPointer ( impChem, cellArea, 'AREA',    rc=status )
   VERIFY_(status)
   call MAPL_GetPointer ( impChem, frlake,   'FRLAKE',  rc=status )
   VERIFY_(status)
   call MAPL_GetPointer ( impChem, frocean,  'FROCEAN', rc=status )
   VERIFY_(status)
   call MAPL_GetPointer ( impChem, frseaice, 'FRACI',   rc=status )
   VERIFY_(status)

!  Get 3D Imports
!  --------------
   call MAPL_GetPointer ( impChem,     T,   'T',       rc=status )
   VERIFY_(status)
   call MAPL_GetPointer ( impChem, hghte,  'ZLE',     rc=status )
   VERIFY_(status)
   call MAPL_GetPointer ( impChem, ple,    'PLE',     rc=status )
   VERIFY_(status)
   call MAPL_GetPointer ( impChem, qccu,   'CNV_QC',  rc=status )
   VERIFY_(status)
   call MAPL_GetPointer ( impChem, cmfmc,  'CNV_MFC', rc=status )
   VERIFY_(status)
   call MAPL_GetPointer ( impChem, dtrain, 'CNV_MFD', rc=status )
   VERIFY_(status)
   call MAPL_GetPointer ( impChem, pfllsan,'PFL_LSAN',rc=status )
   VERIFY_(status)
   call MAPL_GetPointer ( impChem, pfilsan,'PFI_LSAN',rc=status )
   VERIFY_(status)

  IF(gcSC%verbose) THEN
   ic = gcSC%irO3Ox
   CALL pmaxmin('SC:    TROPP ', bsc%qa(ic)%data3d(:,:,km), qmin, qmax, iXj, 1, 0.01 )
   CALL pmaxmin('SC:  	    T ', T,	   qmin, qmax, iXj, km,   1. )
   CALL pmaxmin('SC: 	    Q ', specHum,  qmin, qmax, iXj, km,   1. )
   CALL pmaxmin('SC:    QCTOT ', qctot,	   qmin, qmax, iXj, km,   1. )
   CALL pmaxmin('SC:      RHO ', rhoa,	   qmin, qmax, iXj, km,   1. )
   CALL pmaxmin('SC:     RHOD ', rhoadry,  qmin, qmax, iXj, km,   1. )
   CALL pmaxmin('SC:      PLE ', ple,	   qmin, qmax, iXj, km+1, 1. )
   CALL pmaxmin('SC:     AREA ', cellArea, qmin, qmax, iXj, 1,    1. )
   CALL pmaxmin('SC:      LWI ', oro,      qmin, qmax, iXj, 1,    1. )
   CALL pmaxmin('SC:   frlake ', frlake  , qmin, qmax, iXj, 1,    1. )
   CALL pmaxmin('SC:  frocean ', frocean , qmin, qmax, iXj, 1,    1. )
   CALL pmaxmin('SC: frseaice ', frseaice, qmin, qmax, iXj, 1,    1. )
   CALL pmaxmin('SC:    precl ', precl   , qmin, qmax, iXj, 1,    1. )
   CALL pmaxmin('SC:    precc ', precc   , qmin, qmax, iXj, 1,    1. )
   call pmaxmin('SC:     pblh ', pblh    , qmin, qmax, iXj, 1,    1. )
   call pmaxmin('SC:    hghte ', hghte   , qmin, qmax, iXj, km,   1. )
   call pmaxmin('SC:  pfllsan ', pfllsan , qmin, qmax, iXj, km+1, 1. )
   call pmaxmin('SC:  pfilsan ', pfilsan , qmin, qmax, iXj, km+1, 1. )
  END IF

  IF(gcSC%doFlux) THEN
   CALL MAPL_GetPointer(impChem,    f11ems,    'F11FLUX', RC=status)
   VERIFY_(status)
   CALL MAPL_GetPointer(impChem,    f12ems,    'F12FLUX', RC=status)
   VERIFY_(status)
   CALL MAPL_GetPointer(impChem,   f113ems,   'F113FLUX', RC=status)
   VERIFY_(status)
   CALL MAPL_GetPointer(impChem,    mcfems,    'MCFFLUX', RC=status)
   VERIFY_(status)
   CALL MAPL_GetPointer(impChem,   ccl4ems,   'CCL4FLUX', RC=status)
   VERIFY_(status)
   CALL MAPL_GetPointer(impChem, hcfc22ems, 'HCFC22FLUX', RC=status)
   VERIFY_(status)

#ifndef REDUCED
   CALL MAPL_GetPointer(impChem,  chbr3ems,  'CHBR3FLUX', RC=status)
   VERIFY_(status)
   CALL MAPL_GetPointer(impChem, ch2br2ems, 'CH2BR2FLUX', RC=status)
   VERIFY_(status)
#endif

   CALL readEMSFactor(status)
   VERIFY_(status)

   gcSC%CFC11EMISS(i1:i2,j1:j2)  =    f11ems(i1:i2,j1:j2)*gcSC%emfCFC11
   gcSC%CFC12EMISS(i1:i2,j1:j2)  =    f12ems(i1:i2,j1:j2)*gcSC%emfCFC12
   gcSC%CFC113EMISS(i1:i2,j1:j2) =   f113ems(i1:i2,j1:j2)*gcSC%emfCFC113
   gcSC%MCFEMISS(i1:i2,j1:j2)    =    mcfems(i1:i2,j1:j2)*gcSC%emfMCF
   gcSC%CCL4EMISS(i1:i2,j1:j2)   =   ccl4ems(i1:i2,j1:j2)*gcSC%emfCCL4
   gcSC%HCFC22EMISS(i1:i2,j1:j2) = hcfc22ems(i1:i2,j1:j2)*gcSC%emfHCFC22
#ifndef REDUCED
   gcSC%CHBR3EMISS(i1:i2,j1:j2)  =  chbr3ems(i1:i2,j1:j2)
   gcSC%CH2BR2EMISS(i1:i2,j1:j2) = ch2br2ems(i1:i2,j1:j2)
#endif

  ENDIF

! Read 3D OH climatology from Spivakovsky et al. 2000
  CALL MAPL_GetPointer(impChem, oh3dClim, 'OH3DClim', RC=status)
  VERIFY_(status)

! If importing SO4 surface area from GOCART::SU then query the import state. Units: m^{2} m^{-3}.
! -----------------------------------------------------------------------------------------------
  IF(gcSC%importSulfSA) THEN

   CALL MAPL_GetPointer(impChem,  SO4SAv,  'SO4SAREAvolc', RC=status)
   VERIFY_(status)
   IF(gcSC%verbose) THEN
    CALL pmaxmin('SC: SO4SAvolc ', SO4SAv, qmin, qmax, iXj, km,    1. )
   END IF

   CALL MAPL_GetPointer(impChem,  SO4SAt,  'SO4SAREA', RC=status)
   VERIFY_(status)
   IF(gcSC%verbose) THEN
   CALL pmaxmin('SC: SO4SA ', SO4SAt, qmin, qmax, iXj, km,    1. )
   END IF

  END IF

! Export state pointers
! ---------------------
  CALL MAPL_GetPointer(expChem,      szarad,      'SZARAD', RC=status)
  VERIFY_(status)
  CALL MAPL_GetPointer(expChem,     scbasek,     'SCBASEK', RC=status)
  VERIFY_(status)
  CALL MAPL_GetPointer(expChem,      n2oflx,      'N2OFLX', ALLOC=.TRUE., RC=status)
  VERIFY_(status)
  CALL MAPL_GetPointer(expChem,      ch4flx,      'CH4FLX', ALLOC=.TRUE., RC=status)
  VERIFY_(status)
  CALL MAPL_GetPointer(expChem,      f11flx,      'F11FLX', ALLOC=.TRUE., RC=status)
  VERIFY_(status)
  CALL MAPL_GetPointer(expChem,      f12flx,      'F12FLX', ALLOC=.TRUE., RC=status)
  VERIFY_(status)
  CALL MAPL_GetPointer(expChem,     f113flx,     'F113FLX', ALLOC=.TRUE., RC=status)
  VERIFY_(status)
  CALL MAPL_GetPointer(expChem,   hcfc22flx,   'HCFC22FLX', ALLOC=.TRUE., RC=status)
  VERIFY_(status)

#ifndef REDUCED
  CALL MAPL_GetPointer(expChem, hcfc141bflx, 'HCFC141BFLX', ALLOC=.TRUE., RC=status)
  VERIFY_(status)
  CALL MAPL_GetPointer(expChem, hcfc142bflx, 'HCFC142BFLX', ALLOC=.TRUE., RC=status)
  VERIFY_(status)
#endif

  CALL MAPL_GetPointer(expChem,     ccl4flx,     'CCL4FLX', ALLOC=.TRUE., RC=status)
  VERIFY_(status)
  CALL MAPL_GetPointer(expChem,      mcfflx,      'MCFFLX', ALLOC=.TRUE., RC=status)
  VERIFY_(status)
  CALL MAPL_GetPointer(expChem,    ch3clflx,    'CH3CLFLX', ALLOC=.TRUE., RC=status)
  VERIFY_(status)
  CALL MAPL_GetPointer(expChem,    ch3brflx,    'CH3BRFLX', ALLOC=.TRUE., RC=status)
  VERIFY_(status)
  CALL MAPL_GetPointer(expChem,    h1301flx,    'H1301FLX', ALLOC=.TRUE., RC=status)
  VERIFY_(status)
  CALL MAPL_GetPointer(expChem,    h1211flx,    'H1211FLX', ALLOC=.TRUE., RC=status)
  VERIFY_(status)

#ifndef REDUCED
  CALL MAPL_GetPointer(expChem,    h1202flx,    'H1202FLX', ALLOC=.TRUE., RC=status)
  VERIFY_(status)
  CALL MAPL_GetPointer(expChem,    h2402flx,    'H2402FLX', ALLOC=.TRUE., RC=status)
#endif

  IF(gcSC%doFlux) THEN
   CALL MAPL_GetPointer(expChem,    mcfocnloss,  'MCFOCNLOSS', ALLOC=.TRUE., RC=status)
   VERIFY_(status)
   CALL MAPL_GetPointer(expChem,    ctcocnloss,  'CTCOCNLOSS', ALLOC=.TRUE., RC=status)
   VERIFY_(status)
   CALL MAPL_GetPointer(expChem,    ctclndloss,  'CTCLNDLOSS', ALLOC=.TRUE., RC=status)
   VERIFY_(status)
  ENDIF

  CALL MAPL_GetPointer(expChem,    rdair,  'rDryAir', ALLOC=.TRUE., RC=status)
  VERIFY_(status)

  CALL MAPL_GetPointer(expChem,    qqj001,  'QQJ001', ALLOC=.TRUE., RC=status)
  VERIFY_(status)
  CALL MAPL_GetPointer(expChem,    qqj002,  'QQJ002', ALLOC=.TRUE., RC=status)
  VERIFY_(status)
  CALL MAPL_GetPointer(expChem,    qqj003,  'QQJ003', ALLOC=.TRUE., RC=status)
  VERIFY_(status)
  CALL MAPL_GetPointer(expChem,    qqj004,  'QQJ004', ALLOC=.TRUE., RC=status)
  VERIFY_(status)
  CALL MAPL_GetPointer(expChem,    qqj005,  'QQJ005', ALLOC=.TRUE., RC=status)
  VERIFY_(status)
  CALL MAPL_GetPointer(expChem,    qqj006,  'QQJ006', ALLOC=.TRUE., RC=status)
  VERIFY_(status)
  CALL MAPL_GetPointer(expChem,    qqj007,  'QQJ007', ALLOC=.TRUE., RC=status)
  VERIFY_(status)
  CALL MAPL_GetPointer(expChem,    qqj008,  'QQJ008', ALLOC=.TRUE., RC=status)
  VERIFY_(status)
  CALL MAPL_GetPointer(expChem,    qqj009,  'QQJ009', ALLOC=.TRUE., RC=status)
  VERIFY_(status)
  CALL MAPL_GetPointer(expChem,    qqj010,  'QQJ010', ALLOC=.TRUE., RC=status)
  VERIFY_(status)
  CALL MAPL_GetPointer(expChem,    qqj011,  'QQJ011', ALLOC=.TRUE., RC=status)
  VERIFY_(status)
  CALL MAPL_GetPointer(expChem,    qqj012,  'QQJ012', ALLOC=.TRUE., RC=status)
  VERIFY_(status)
  CALL MAPL_GetPointer(expChem,    qqj013,  'QQJ013', ALLOC=.TRUE., RC=status)
  VERIFY_(status)
  CALL MAPL_GetPointer(expChem,    qqj014,  'QQJ014', ALLOC=.TRUE., RC=status)
  VERIFY_(status)
  CALL MAPL_GetPointer(expChem,    qqj015,  'QQJ015', ALLOC=.TRUE., RC=status)
  VERIFY_(status)
  CALL MAPL_GetPointer(expChem,    qqj016,  'QQJ016', ALLOC=.TRUE., RC=status)
  VERIFY_(status)
  CALL MAPL_GetPointer(expChem,    qqj017,  'QQJ017', ALLOC=.TRUE., RC=status)
  VERIFY_(status)
  CALL MAPL_GetPointer(expChem,    qqj018,  'QQJ018', ALLOC=.TRUE., RC=status)
  VERIFY_(status)
  CALL MAPL_GetPointer(expChem,    qqj019,  'QQJ019', ALLOC=.TRUE., RC=status)
  VERIFY_(status)
  CALL MAPL_GetPointer(expChem,    qqj020,  'QQJ020', ALLOC=.TRUE., RC=status)
  VERIFY_(status)

  CALL MAPL_GetPointer(expChem,    qqk001,  'QQK001', ALLOC=.TRUE., RC=status)
  VERIFY_(status)
  CALL MAPL_GetPointer(expChem,    qqk002,  'QQK002', ALLOC=.TRUE., RC=status)
  VERIFY_(status)
  CALL MAPL_GetPointer(expChem,    qqk003,  'QQK003', ALLOC=.TRUE., RC=status)
  VERIFY_(status)
  CALL MAPL_GetPointer(expChem,    qqk004,  'QQK004', ALLOC=.TRUE., RC=status)
  VERIFY_(status)
  CALL MAPL_GetPointer(expChem,    qqk005,  'QQK005', ALLOC=.TRUE., RC=status)
  VERIFY_(status)
  CALL MAPL_GetPointer(expChem,    qqk006,  'QQK006', ALLOC=.TRUE., RC=status)
  VERIFY_(status)
  CALL MAPL_GetPointer(expChem,    qqk007,  'QQK007', ALLOC=.TRUE., RC=status)
  VERIFY_(status)
  CALL MAPL_GetPointer(expChem,    qqk008,  'QQK008', ALLOC=.TRUE., RC=status)
  VERIFY_(status)
  CALL MAPL_GetPointer(expChem,    qqk009,  'QQK009', ALLOC=.TRUE., RC=status)
  VERIFY_(status)
  CALL MAPL_GetPointer(expChem,    qqk010,  'QQK010', ALLOC=.TRUE., RC=status)
  VERIFY_(status)
  CALL MAPL_GetPointer(expChem,    qqk011,  'QQK011', ALLOC=.TRUE., RC=status)
  VERIFY_(status)
  CALL MAPL_GetPointer(expChem,    qqk012,  'QQK012', ALLOC=.TRUE., RC=status)
  VERIFY_(status)
  CALL MAPL_GetPointer(expChem,    qqk013,  'QQK013', ALLOC=.TRUE., RC=status)
  VERIFY_(status)
  CALL MAPL_GetPointer(expChem,    qqk014,  'QQK014', ALLOC=.TRUE., RC=status)
  VERIFY_(status)
  CALL MAPL_GetPointer(expChem,    qqk015,  'QQK015', ALLOC=.TRUE., RC=status)
  VERIFY_(status)
  CALL MAPL_GetPointer(expChem,    qqk016,  'QQK016', ALLOC=.TRUE., RC=status)
  VERIFY_(status)
  CALL MAPL_GetPointer(expChem,    qqk017,  'QQK017', ALLOC=.TRUE., RC=status)
  VERIFY_(status)
  CALL MAPL_GetPointer(expChem,    qqk018,  'QQK018', ALLOC=.TRUE., RC=status)
  VERIFY_(status)
  CALL MAPL_GetPointer(expChem,    qqk019,  'QQK019', ALLOC=.TRUE., RC=status)
  VERIFY_(status)
  CALL MAPL_GetPointer(expChem,    qqk020,  'QQK020', ALLOC=.TRUE., RC=status)
  VERIFY_(status)
  CALL MAPL_GetPointer(expChem,    qqk021,  'QQK021', ALLOC=.TRUE., RC=status)
  VERIFY_(status)
  CALL MAPL_GetPointer(expChem,    qqk022,  'QQK022', ALLOC=.TRUE., RC=status)
  VERIFY_(status)
  CALL MAPL_GetPointer(expChem,    qqk023,  'QQK023', ALLOC=.TRUE., RC=status)
  VERIFY_(status)
  CALL MAPL_GetPointer(expChem,    qqk024,  'QQK024', ALLOC=.TRUE., RC=status)
  VERIFY_(status)
  CALL MAPL_GetPointer(expChem,    qqk025,  'QQK025', ALLOC=.TRUE., RC=status)
  VERIFY_(status)
  CALL MAPL_GetPointer(expChem,    qqk026,  'QQK026', ALLOC=.TRUE., RC=status)
  VERIFY_(status)
  CALL MAPL_GetPointer(expChem,    qqk027,  'QQK027', ALLOC=.TRUE., RC=status)
  VERIFY_(status)
  CALL MAPL_GetPointer(expChem,    qqk028,  'QQK028', ALLOC=.TRUE., RC=status)
  VERIFY_(status)
  CALL MAPL_GetPointer(expChem,    qqk029,  'QQK029', ALLOC=.TRUE., RC=status)
  VERIFY_(status)
  CALL MAPL_GetPointer(expChem,    qqk030,  'QQK030', ALLOC=.TRUE., RC=status)
  VERIFY_(status)
  CALL MAPL_GetPointer(expChem,    qqk031,  'QQK031', ALLOC=.TRUE., RC=status)
  VERIFY_(status)
  CALL MAPL_GetPointer(expChem,    qqk032,  'QQK032', ALLOC=.TRUE., RC=status)
  VERIFY_(status)
  CALL MAPL_GetPointer(expChem,    qqk033,  'QQK033', ALLOC=.TRUE., RC=status)
  VERIFY_(status)
  CALL MAPL_GetPointer(expChem,    qqk034,  'QQK034', ALLOC=.TRUE., RC=status)
  VERIFY_(status)
  CALL MAPL_GetPointer(expChem,    qqk035,  'QQK035', ALLOC=.TRUE., RC=status)

  CALL MAPL_GetPointer(expChem,  o3ppmv,    'O3PPMV', RC=status)
  VERIFY_(status)
  CALL MAPL_GetPointer(expChem,   ozone,        'O3', RC=status)
  VERIFY_(status)
  CALL MAPL_GetPointer(expChem,   dOxdt,   'OX_TEND', RC=status)
  VERIFY_(status)
  CALL MAPL_GetPointer(expChem,    dQdt,  'H2O_TEND', RC=status)
  VERIFY_(status)
  CALL MAPL_GetPointer(expChem,    jocs, 'OCS_JRATE', ALLOC=.TRUE., RC=status)
  VERIFY_(status)

! PSC preliminaries
! -----------------
  IF(gcSC%doPSCs) THEN
   CALL MAPL_GetPointer(expChem, natsad, 'NATSAD', ALLOC=.TRUE., RC=status)
   VERIFY_(status)
   natsad(:,:,:) = 0.00
   CALL MAPL_GetPointer(expChem, icesad, 'ICESAD', ALLOC=.TRUE., RC=status)
   VERIFY_(status)
   icesad(:,:,:) = 0.00
  END IF

! Initialize tendency diagnostics
! -------------------------------
  IF(ASSOCIATED(dQdt)) dQdt(i1:i2,j1:j2,1:km) = specHum(i1:i2,j1:j2,1:km)
  IF(ASSOCIATED(dOxdt)) dOxdt(i1:i2,j1:j2,1:km) = bsc%qa(gcSC%iOx)%data3d(i1:i2,j1:j2,1:km)

! Periodically update the climatological boundary conditions
! ----------------------------------------------------------
  BCUpdates: IF(gcSC%BCnymd /= nymd) THEN

   IF(MAPL_AM_I_ROOT()) PRINT *,TRIM(Iam)//": Updating climatologies"

   CALL readClimVars(status)
   VERIFY_(status)

#ifndef REDUCED
   CALL mapToSurf(gcSC%CO2MBC, gcSC%CO2MBC2D, "CO2MBC", status)
   VERIFY_(status)
   CALL mapToSurf(gcSC%SF6MBC, gcSC%SF6MBC2D, "SF6MBC", status)
   VERIFY_(status)
#endif

   CALL mapToGrid(gcSC%climRain, gcSC%climRain3D, "Rain", status)
   VERIFY_(status)
   CALL mapToGrid(gcSC%climO3, gcSC%climO33D, "O3", status)
   VERIFY_(status)
   CALL mapToGrid(gcSC%climOH, gcSC%climOH3D, "OH", status)
   VERIFY_(status)
   CALL mapToGrid(gcSC%climH2O, gcSC%climH2O3D, "H2O", status)
   VERIFY_(status)

#ifndef REDUCED
   CALL mapToGrid(gcSC%climSF6LOSS, gcSC%climSF6LOSS3D, "SF6LOSS", status)
   VERIFY_(status)
#endif

   CALL mapToGrid(gcSC%SO4sa, gcSC%SO4sa3D, "SO4sa", status)
   VERIFY_(status)

   gcSC%BCnymd = nymd
  END IF BCUpdates

! Except for the ozone-to-odd oxygen ratio, constituents should be non-negative
! -----------------------------------------------------------------------------
  DO i = 1,scReg%nq-1
   WHERE(bsc%qa(i)%data3d(i1:i2,j1:j2,1:km) <= MIN_VALUE ) &
         bsc%qa(i)%data3d(i1:i2,j1:j2,1:km)  = MIN_VALUE
  END DO
  WHERE(specHum(i1:i2,j1:j2,1:km) <= 1.00E-25) specHum(i1:i2,j1:j2,1:km) = MIN_VALUE

! Cycle though the grid cells on this processor
! ---------------------------------------------
  Latitude: DO j = j1,j2
   Longitude: DO i = i1,i2

! Solar zenith angle (radians) and daytime indicator.
! ---------------------------------------------------
    CALL sza(gcSC%dayOfYear,gcSC%lonRad(i,j),gcSC%latRad(i,j),szan,gcSC%ecc_SC,rdistsq)
    dayTime = (szan <= gcSC%szaCutoff*MAPL_PI/180.00)
    IF(ASSOCIATED(szarad)) szarad(i,j) = szan

! Profiles. Temperature [K] and pressure [hPa]
! --------------------------------------------
    kel(1:km) = T(i,j,km:1:-1)
    edgePress(0:km) = ple(i,j,km:0:-1)*0.01
    midPress(1:km) = 0.50*(edgePress(0:km-1)+edgePress(1:km))

! Number density [cm^-3] in moist and dry air
! -------------------------------------------
!   numDens(1:km) = rhoa(i,j,km:1:-1)*MAPL_AVOGAD*1.00E-06/MAPL_AIRMW
    numDens(1:km) = (     rhoadry(i,j,km:1:-1)/MAPL_AIRMW + &
        specHum(i,j,km:1:-1)*rhoa(i,j,km:1:-1)/MAPL_H2OMW ) * MAPL_AVOGAD * 1.00E-06
    numDensDry(1:km) = rhoadry(i,j,km:1:-1)*MAPL_AVOGAD*1.00E-06/MAPL_AIRMW

    rdair(i,j,km:1:-1) = numDensDry(1:km)/numDens(1:km) 

! Cell depth [cm] and volume [m^3]
! --------------------------------
    cellDepth(1:km) = 100.00*MAPL_RGAS*kel(1:km)*LOG(edgePress(0:km-1)/edgePress(1:km))/MAPL_GRAV
    cellVolume(1:km) = 0.01*cellDepth(1:km)*cellArea(i,j)

! A logical flag to indicate presence in the troposphere.
! -------------------------------------------------------
    tropo(:) = .FALSE.
    ic = gcSC%irO3Ox
    tropp = bsc%qa(ic)%data3d(i,j,km)*0.01
    WHERE(edgePress(1:km) >= tropp) tropo(1:km) = .TRUE.

! Find layer number at base of stratosphere, and verify tropo=.TRUE. at each layer below
! --------------------------------------------------------------------------------------
    k1Strat = km
    DO k = km,1,-1
     IF(tropo(k)) EXIT
     k1Strat = k
    END DO
    tropo(1:k1Strat-1) = .TRUE.
    IF(ASSOCIATED(scbasek)) scbasek(i,j) = km-k1Strat+1

! Ozone and NO number densities
! -----------------------------
    O3profile(1:km) = bsc%qa(gcSC%iO3)%data3d(i,j,km:1:-1)*numDens(1:km)
    NOprofile(1:km) = bsc%qa(gcSC%iNO)%data3d(i,j,km:1:-1)*numDens(1:km)

! Overhead ozone and molecular oxygen
! -----------------------------------
    CALL O3Colum(km,edgePress,midPress,numDens,O3Profile,NOProfile,O3Column,O2Column,NOColumn)

! Condensed water [mole fraction] for SUBROUTINE sediment
! -------------------------------------------------------
    H2OcProfile(1:km) = qctot(i,j,km:1:-1)*MAPL_AIRMW/MAPL_H2OMW

! Profiles of the climatological fields
! -------------------------------------
    climRainProfile(1:km) = gcSC%climRain3D(i,j,km:1:-1)
    climO3Profile(1:km) = gcSC%climO33D(i,j,km:1:-1)

!    climOHProfile(1:km) = gcSC%climOH3D(i,j,km:1:-1)
!    climOHProfile(1:km) = OH3DClim(i,j,km:1:-1)
    climOHProfile(1:km) = oh3dClim(i,j,km:1:-1)*numDens(1:km)/numDensDry(1:km)

    climH2OProfile(1:km) = gcSC%climH2O3D(i,j,km:1:-1)

#ifndef REDUCED
    climSF6LOSSProfile(1:km) = gcSC%climSF6LOSS3D(i,j,km:1:-1)
#endif

    SO4saProfile(1:km) = gcSC%SO4sa3D(i,j,km:1:-1)

! Utilize Sulfate SA imported from GOCART::SU.
! Convert units from  m^{2} m^{-3} to cm^{2} cm^{-3}
! --------------------------------------------------
    IF(gcSC%importSulfSA) THEN
     IF(gcSC%addSAclim) THEN
!      SO4saProfile(1:km) = SO4saProfile(1:km) + 0.01 * SO4SAimp(i,j,km:1:-1)
      SO4saProfile(1:km) = SO4saProfile(1:km) + 0.01 * (SO4SAt(i,j,km:1:-1)+SO4SAv(i,j,km:1:-1))
     ELSE
!      SO4saProfile(1:km) =                      0.01 * SO4SAimp(i,j,km:1:-1)
      SO4saProfile(1:km) =                      0.01 * (SO4SAt(i,j,km:1:-1)+SO4SAv(i,j,km:1:-1))
     END IF
    END IF

    IF(gcSC%useSolCyc) THEN
     gcrProfile(1:km) = 0.00
!JEN gcrProfile(1:km) = gcSC%climGCR3D(i,j,km:1:-1) --- Awaits implementation updates JEN!
    ELSE
     gcrProfile(1:km) = 0.00
    END IF

! Calculate thermal reaction rates. SUBROUTINE khet3d
! recalculates reactions numbers 110-115 when PSCS are enabled.
! -------------------------------------------------------------
    CALL krates(gcSC,km,numDens,kel,ak)

! If use flux-based boundary conditions, then turn on surface losses for MCF and CCl4

    IF(gcSC%doFlux) THEN

! Calculate column density - unit: molec cm-2
! cellDepth is in unit of cm
     mcfColumn=0.0
     ctcColumn=0.0
     mcfprofile(1:km) = bsc%qa(gcSC%iCH3CCl3)%data3d(i,j,km:1:-1)*numDens(1:km)*cellDepth(1:km)
     ctcprofile(1:km) = bsc%qa(gcSC%iCCl4)%data3d(i,j,km:1:-1)*numDens(1:km)*cellDepth(1:km)
     DO k=1,km
      mcfColumn=mcfColumn+mcfprofile(k)
      ctcColumn=ctcColumn+ctcprofile(k)
     END DO

    END IF


! At each layer ...
! -----------------
    Layer: DO k=1,km
     kRev = km-k+1

! Place the species in reverse vertical order and convert to number density.
! --------------------------------------------------------------------------
     h2o = specHum(i,j,kRev)*numDens(k)*MAPL_AIRMW/MAPL_H2OMW
     h2oc =  qctot(i,j,kRev)*numDens(k)*MAPL_AIRMW/MAPL_H2OMW
     CALL swapSpecies(FROM_BUNDLE,status)
     VERIFY_(status)

! Initializations to complete a gcSC%spinup (rarely used):
!   1. Constrain HCl to be at least 1.e-11
!   2. Ox and O3 are initialized identically.
!   3. ClO and NO2 are used to determine ak(32).
! --------------------------------------------------------
     IF(gcSC%spinup) THEN
      hcl = MAX(hcl,1.00E-11*numDens(k))
       o3 = MAX(ox,1.00*numDens(k))
      clo = 0.00
      no2 = 0.00
     END IF

! Fix the surface source gases and infer the surface flux [kg s^{-1}]
! -------------------------------------------------------------------
     IF(k == 1) THEN
      CALL SSGSfcFlx(cellDepth(1),cellVolume(1),numDens(1),numDensDry(1),n2o,ch4, &
                     f11,f12,f113,f114,f115,ccl4,ch3ccl3,ch3cl,ch3br,h1301,h1211, &
                       h1202,h2402,hcfc22,hcfc141b,hcfc142b,chbr3,ch2br2,ch2brcl, &
                      chbrcl2,chbr2cl,co2,sf6,hfc23,hfc32,hfc125,hfc134a,hfc143a, &
                                     hfc152a,mcfcolumn,ctccolumn,oro(i,j),status)

      VERIFY_(status)
     END IF

! Update photolysis rates.
! ------------------------
     CALL jcalc4(k,szan,O3Column(k),midPress(k),kel(k),aj,gcSC,rdistsq)
     CALL meso_phot(k,km,0,NOColumn,O3Column,O2Column,szan,aj,numDens,gcSC,rdistsq)

     qqj(:) = 0.00
     qqk(:) = 0.00

! Lightning source of NOx between 20 S and 20 N in troposphere.
! Note: The rate can be adjusted at run time in SC_GridComp.rc.
! -------------------------------------------------------------
     IF(gcSC%latRad(i,j) >= -20.00*gcSC%dtr .AND. gcSC%latRad(i,j) <= 20.00*gcSC%dtr) THEN
      IF(k >= 5 .AND. tropo(k)) nox = nox+gcSC%NOxTRW*tdt
     END IF

! Store initial values in case time step is changed in the main chemistry loop.
! -----------------------------------------------------------------------------
     ALLOCATE(ratesBase(gcSC%numphoto+gcSC%numreacs),STAT=status)
     VERIFY_(status)
     ratesBase(1:gcSC%numphoto) = aj(1:gcSC%numphoto)
     ratesBase(gcSC%numphoto+1:gcSC%numphoto+gcSC%numreacs) = ak(k,1:gcSC%numreacs)
     ALLOCATE(speciesBase(scReg%nq+1),STAT=status)
     VERIFY_(status)
     CALL storeBaseChem(ox,nox,hno3,n2o5,ho2no2,clono2,clx,hcl,hocl,h2o2,brx,n2o,cl2, &
                 oclo,brcl,hbr,brono2,ch4,hobr,ch3ooh,co,h2o,hno3c,h2oc,f11,f12,f113, &
                f114,f115,ccl4,hcfc22,hcfc141b,hcfc142b,chbr3,ch2br2,ch2brcl,chbrcl2, &
                          chbr2cl,hfc23,hfc32,hfc125,hfc134a,hfc143a,hfc152a,co2,sf6, &
               ch3ccl3,ch3cl,ch3br,h1301,h1211,h1202,h2402,o3,o3p,o1d,n,no,no2,no3,h, & 
              oh,ho2,cl,clo,bro,br,cl2o2,ch2o,ch3o2,box_ro3ox,speciesBase,scReg%nq)

! Variables needed for chemistry time step length analysis.  These are changed
! within the while-do loop if the given time step is too long to ensure stability.
! --------------------------------------------------------------------------------
     dt = tdt
     requiredDt = tdt
     requiredPasses = 1
     passNumber = 1
     examineDt = 1
     startOver = .FALSE.

! For the current cell ...
! ------------------------
     ThisCell: DO WHILE(passNumber <= requiredPasses)
      
! Run the mechanism
! -----------------
      Chemistry: DO

       IF(startOver) THEN
        aj(1:gcSC%numphoto) = ratesBase(1:gcSC%numphoto)
        ak(k,1:gcSC%numreacs) = ratesBase(gcSC%numphoto+1:gcSC%numphoto+gcSC%numreacs)
        CALL getBaseChem(ox,nox,hno3,n2o5,ho2no2,clono2,clx,hcl,hocl,h2o2,brx,n2o,cl2, &
                  oclo,brcl,hbr,brono2,ch4,hobr,ch3ooh,co,h2o,hno3c,h2oc,f11,f12,f113, &
                 f114,f115,ccl4,hcfc22,hcfc141b,hcfc142b,chbr3,ch2br2,ch2brcl,chbrcl2, &
                          chbr2cl,hfc23,hfc32,hfc125,hfc134a,hfc143a,hfc152a,co2,sf6,  &
                ch3ccl3,ch3cl,ch3br,h1301,h1211,h1202,h2402,o3,o3p,o1d,n,no,no2,no3,h, &           
               oh,ho2,cl,clo,bro,br,cl2o2,ch2o,ch3o2,box_ro3ox,speciesBase,scReg%nq)
        passNumber = 1
       END IF

       startOver = .FALSE.
       CALL chemStage1(status)
       VERIFY_(status)
       IF(startOver) CYCLE

       CALL chemStage2(status)
       VERIFY_(status)
       IF(startOver) CYCLE
       EXIT

      END DO Chemistry

! Calculate SF6 photochemical loss using climatological loss rates from GSFC2D
! ----------------------------------------------------------------------------
#ifndef REDUCED
      sf6 = sf6 - sf6*climSF6LOSSProfile(k)*tdt
!   Should increase the loss when addressing "lifetime" concerns
!     sf6 = sf6 - sf6*10*climSF6LOSSProfile(k)*tdt 
#endif

! Parameterized scavenging of species in the troposphere
! ------------------------------------------------------
      IF(k < k1Strat) THEN
           tx = 1.00-climRainProfile(k)*tdt
          hcl =    hcl*tx
       clono2 = clono2*tx
         hno3 =   hno3*tx
         h2o2 =   h2o2*tx
       ho2no2 = ho2no2*tx
!         hbr =    hbr*tx
      END IF

! Relax tropospheric ozone to climatology with a time constant of 5 days.
! -----------------------------------------------------------------------
      IF(k < k1Strat) ox = ox+(climO3Profile(k)*numDens(k)-ox)*tdt/(5.00*86400.00)

! Archive thermal and photolysis loss rates for ODSs
! -----------------------------------------------------------------------
   qqj001(i,j,kRev) = qqj( 1)/(numDens(k)*tdt)
   qqj002(i,j,kRev) = qqj( 2)/(numDens(k)*tdt)
   qqj003(i,j,kRev) = qqj( 3)/(numDens(k)*tdt)
   qqj004(i,j,kRev) = qqj( 4)/(numDens(k)*tdt)
   qqj005(i,j,kRev) = qqj( 5)/(numDens(k)*tdt)
   qqj006(i,j,kRev) = qqj( 6)/(numDens(k)*tdt)
   qqj007(i,j,kRev) = qqj( 7)/(numDens(k)*tdt)
   qqj008(i,j,kRev) = qqj( 8)/(numDens(k)*tdt)
   qqj009(i,j,kRev) = qqj( 9)/(numDens(k)*tdt)
   qqj010(i,j,kRev) = qqj(10)/(numDens(k)*tdt)
   qqj011(i,j,kRev) = qqj(11)/(numDens(k)*tdt)
   qqj012(i,j,kRev) = qqj(12)/(numDens(k)*tdt)
   qqj013(i,j,kRev) = qqj(13)/(numDens(k)*tdt)
   qqj014(i,j,kRev) = qqj(14)/(numDens(k)*tdt)
   qqj015(i,j,kRev) = qqj(15)/(numDens(k)*tdt)
   qqj016(i,j,kRev) = qqj(16)/(numDens(k)*tdt)
   qqj017(i,j,kRev) = qqj(17)/(numDens(k)*tdt)
   qqj018(i,j,kRev) = qqj(18)/(numDens(k)*tdt)
   qqj019(i,j,kRev) = qqj(19)/(numDens(k)*tdt)
   qqj020(i,j,kRev) = qqj(20)/(numDens(k)*tdt)

   qqk001(i,j,kRev) = qqk( 1)/(numDens(k)*tdt)
   qqk002(i,j,kRev) = qqk( 2)/(numDens(k)*tdt)
   qqk003(i,j,kRev) = qqk( 3)/(numDens(k)*tdt)
   qqk004(i,j,kRev) = qqk( 4)/(numDens(k)*tdt)
   qqk005(i,j,kRev) = qqk( 5)/(numDens(k)*tdt)
   qqk006(i,j,kRev) = qqk( 6)/(numDens(k)*tdt)
   qqk007(i,j,kRev) = qqk( 7)/(numDens(k)*tdt)
   qqk008(i,j,kRev) = qqk( 8)/(numDens(k)*tdt)
   qqk009(i,j,kRev) = qqk( 9)/(numDens(k)*tdt)
   qqk010(i,j,kRev) = qqk(10)/(numDens(k)*tdt)
   qqk011(i,j,kRev) = qqk(11)/(numDens(k)*tdt)
   qqk012(i,j,kRev) = qqk(12)/(numDens(k)*tdt)
   qqk013(i,j,kRev) = qqk(13)/(numDens(k)*tdt)
   qqk014(i,j,kRev) = qqk(14)/(numDens(k)*tdt)
   qqk015(i,j,kRev) = qqk(15)/(numDens(k)*tdt)
   qqk016(i,j,kRev) = qqk(16)/(numDens(k)*tdt)
   qqk017(i,j,kRev) = qqk(17)/(numDens(k)*tdt)
   qqk018(i,j,kRev) = qqk(18)/(numDens(k)*tdt)
   qqk019(i,j,kRev) = qqk(19)/(numDens(k)*tdt)
   qqk020(i,j,kRev) = qqk(20)/(numDens(k)*tdt)
   qqk021(i,j,kRev) = qqk(21)/(numDens(k)*tdt)
   qqk022(i,j,kRev) = qqk(22)/(numDens(k)*tdt)
   qqk023(i,j,kRev) = qqk(23)/(numDens(k)*tdt)
   qqk024(i,j,kRev) = qqk(24)/(numDens(k)*tdt)
   qqk025(i,j,kRev) = qqk(25)/(numDens(k)*tdt)
   qqk026(i,j,kRev) = qqk(26)/(numDens(k)*tdt)
   qqk027(i,j,kRev) = qqk(27)/(numDens(k)*tdt)
   qqk028(i,j,kRev) = qqk(28)/(numDens(k)*tdt)
   qqk029(i,j,kRev) = qqk(29)/(numDens(k)*tdt)
   qqk030(i,j,kRev) = qqk(30)/(numDens(k)*tdt)
   qqk031(i,j,kRev) = qqk(31)/(numDens(k)*tdt)
   qqk032(i,j,kRev) = qqk(32)/(numDens(k)*tdt)
   qqk033(i,j,kRev) = qqk(33)/(numDens(k)*tdt)
   qqk034(i,j,kRev) = qqk(34)/(numDens(k)*tdt)
   qqk035(i,j,kRev) = qqk(35)/(numDens(k)*tdt)


! End of loop to update mixing ratios in the current cell
! -------------------------------------------------------
      passNumber = passNumber+1

     END DO ThisCell

! Housekeeping
! ------------
     DEALLOCATE(ratesBase,STAT=status)
     VERIFY_(status)
     DEALLOCATE(speciesBase,STAT=status)
     VERIFY_(status)

! Update the internal state bundle
! --------------------------------
     CALL swapSpecies(TO_BUNDLE,status)
     VERIFY_(status)

! Update specific humidity
! ------------------------
     specHum(i,j,kRev) = h2o*MAPL_H2OMW/(numDens(k)*MAPL_AIRMW)

!  OCS photolysis rates (for ACHEM)
!  --------------------------------
   IF(ASSOCIATED(jocs)) jocs(i,j,kRev) = aj(55) 

! End of chemistry at layer k of the current (i,j).
! -------------------------------------------------
    END DO Layer

! Sedimentation of HNO3c PSC particles. MOIST handles H2O sedimentation.
! ----------------------------------------------------------------------
    IF(gcSC%doPSCs .AND. gcSC%doSediment) THEN
     HNO3cProfile(1:km) = bsc%qa(gcSC%iHNO3c)%data3d(i,j,km:1:-1)
     CALL sediment(km,tdt,cellDepth,numDens,kel,HNO3CProfile,H2OcProfile, &
                   rmedsts,rmednat,rmedice,denssts,vfall,gcSC)
     bsc%qa(gcSC%iHNO3c)%data3d(i,j,km:1:-1) = HNO3cProfile(1:km)
    END IF

! Move to next cell on this processor
! -----------------------------------
   END DO Longitude
  END DO Latitude
 
!  Large-scale Wet Removal
!  --------------------------

   allocate( fluxout, stat=status )
   VERIFY_(status)
!  allocate( fluxout%data2d(i1:i2,j1:j2),stat=STATUS)  ALLOC IF WE PLAN TO USE IT
!  VERIFY_(STATUS)
   NULLIFY(  fluxout%data2d )

!  It would be nice if WetRemovalGOCART could take a Species_Array instead of Chem_Array
   allocate(qa_single(1),stat=status)
   VERIFY_(status)
   allocate(qa_single(1)%data3d(i1:i2,j1:j2,km),stat=status )
   VERIFY_(status)

   qa_single(1)%fwet = 0.0

!  HBr
   species_index = gcSC%iHBr

   qa_single(1)%data3d(:,:,:) = bsc%qa(species_index)%data3d(:,:,:)

   KIN = .TRUE.
   call WetRemovalGOCART(i1, i2, j1, j2, km, 1, 1, tdt, 'bromine', KIN,  &
              qa_single, ple, t, rhoa, pfllsan, pfilsan, &
              precc, precl, fluxout, rc )

   bsc%qa(species_index)%data3d(:,:,:) = qa_single(1)%data3d(:,:,:)

!  HOBr
   species_index = gcSC%iHOBr
   qa_single(1)%data3d(:,:,:) = bsc%qa(species_index)%data3d(:,:,:)
 
   KIN = .TRUE.
   call WetRemovalGOCART(i1, i2, j1, j2, km, 1, 1, tdt, 'bromine', KIN,  &
              qa_single, ple, t, rhoa, pfllsan, pfilsan, &
              precc, precl, fluxout, rc )

   bsc%qa(species_index)%data3d(:,:,:) = qa_single(1)%data3d(:,:,:)

!  BrONO2
   species_index = gcSC%iBrONO2

   qa_single(1)%data3d(:,:,:) = bsc%qa(species_index)%data3d(:,:,:)

   KIN = .TRUE.
   call WetRemovalGOCART(i1, i2, j1, j2, km, 1, 1, tdt, 'bromine', KIN,  &
              qa_single, ple, t, rhoa, pfllsan, pfilsan, &
              precc, precl, fluxout, rc )

   bsc%qa(species_index)%data3d(:,:,:) = qa_single(1)%data3d(:,:,:)

   deallocate( fluxout, stat=status )
   VERIFY_(status)

!   call pmaxmin('HBr: q_wet', bsc%qa(iHBr)%data3d(i1:i2,j1:j2,1:km), qmin, qmax, &
!                    iXj, km, 1. )

   deallocate(qa_single(1)%data3d,stat=status )
   VERIFY_(status)
   deallocate(qa_single,stat=status)
   VERIFY_(status)

! Convective-scale Wet Removal
   KIN = .TRUE.
   icdt = tdt
   allocate(cmfmc_(i1:i2,j1:j2,km+1), qccu_(i1:i2,j1:j2,km), &
            dtrain_(i1:i2,j1:j2,km), airmass_(i1:i2,j1:j2,km), &
            delz_(i1:i2,j1:j2,km), vud_(i1:i2,j1:j2,km), &
            tc_(i1:i2,j1:j2,km,3), delp_(i1:i2,j1:j2,km), &
            airmol_(i1:i2,j1:j2,km), tmpu_(i1:i2,j1:j2,km), &
            bcnv_(i1:i2,j1:j2,3), ple_(i1:i2,j1:j2,km+1), &
            area_(i1:i2,j1:j2), frlake_(i1:i2,j1:j2), &
            frocean_(i1:i2,j1:j2), frseaice_(i1:i2,j1:j2),&
            h2o2_(i1:i2,j1:j2,km), __STAT__ )

   area_     = cellArea
   frlake_   = frlake
   frocean_  = frocean
   frseaice_ = frseaice   

   do k = 1, km+1
    cmfmc_(:,:,k)   = cmfmc(:,:,km-k+1)
    ple_(:,:,k)     = ple(:,:,km-k+1)
   end do
   do k = 1, km
    dtrain_(:,:,k)  = dtrain(:,:,km-k+1)
    qccu_(:,:,k)    = qccu(:,:,km-k+1)
    delp_(:,:,k)    = bsc%delp(:,:,km-k+1)/100.
    airmass_(:,:,k) = bsc%delp(:,:,km-k+1)/MAPL_GRAV*area_
    airmol_(:,:,k)  = airmass_(:,:,k)*1000./28.966
    delz_(:,:,k)    = bsc%delp(:,:,km-k+1)/MAPL_GRAV/rhoa(:,:,km-k+1)
    tmpu_(:,:,k)    = T(:,:,km-k+1)
   enddo

   do k = 1, km
     tc_(:,:,k,1)   = bsc%qa(gcSC%iHBr)%data3d(:,:,km-k+1)
     tc_(:,:,k,2)   = bsc%qa(gcSC%iHOBr)%data3d(:,:,km-k+1)
     tc_(:,:,k,3)   = bsc%qa(gcSC%iBrONO2)%data3d(:,:,km-k+1)  
   enddo

   call set_vud(i1, i2, j1, j2, km, frlake_, frocean_, frseaice_, cmfmc_, qccu_, &
                airmass_, delz_, area_, vud_)
   call convection(i1, i2, j1, j2, km, 1, 3, icdt, 'bromine', kin, &
                   tc_, cmfmc_, dtrain_, area_, delz_, delp_, vud_, &
                   airmass_, airmol_, tmpu_, ple_, &
                   bcnv_)

! Return adjusted tracer to mixing ratio.
! ---------------------------------------
   do k = 1, km
     bsc%qa(gcSC%iHBr)%data3d(:,:,km-k+1) = tc_(:,:,k,1)
     bsc%qa(gcSC%iHOBr)%data3d(:,:,km-k+1) = tc_(:,:,k,2)
     bsc%qa(gcSC%iBrONO2)%data3d(:,:,km-k+1) = tc_(:,:,k,3)
   enddo


! This switch (rarely used) is set to false after one pass.
! ---------------------------------------------------------
  gcSC%spinup = .FALSE.

!  ------------------
!  Fill export states
!  ------------------

! ---------------------------------------------------------------------
!
!     NAME   STATE    Units             Comments
! -------- --------   ----------------- --------------------------------
!       OX  Internal  mol/mol           Reqired name for ANALYSIS bundle
!       O3  Export    kg/kg             O3CHEM(vmr)*48/28.97
!   03PPMV  Export    ppmv              O3CHEM(vmr)*1.00E+06
!  OX_TEND  Export    mol/(mol s)       Odd-oxygen tendency
! H2O_TEND  Export    kg/(kg s)         Specific humidity tendency
!
! ---------------------------------------------------------------------

!  Ozone mass mixing ratio
!  -----------------------
   IF(ASSOCIATED(ozone)) &
    ozone(i1:i2,j1:j2,1:km) = bsc%qa(gcSC%iO3)%data3d(i1:i2,j1:j2,1:km) * MAPL_O3MW/MAPL_AIRMW 

!  Ozone mole fraction in ppm
!  --------------------------
   IF(ASSOCIATED(o3ppmv)) &
    o3ppmv(i1:i2,j1:j2,1:km) = bsc%qa(gcSC%iO3)%data3d(i1:i2,j1:j2,1:km) * 1.00E+06

!  Odd-oxygen tendency, volume mixing ratio s^{-1}
!  -----------------------------------------------
   IF(ASSOCIATED(dOxdt)) dOxdt(i1:i2,j1:j2,1:km) = &
    (bsc%qa(gcSC%iOx)%data3d(i1:i2,j1:j2,1:km) - dOxdt(i1:i2,j1:j2,1:km))/tdt

! Specific humidity tendency due to chemistry
! -------------------------------------------
   IF(ASSOCIATED(dQdt)) &
    dQdt(i1:i2,j1:j2,1:km) = (specHum(i1:i2,j1:j2,1:km)-dQdt(i1:i2,j1:j2,1:km))/tdt


! Kill the clock
! --------------
   CALL ESMF_ClockDestroy(SCClock, RC=status)
   VERIFY_(status)

   RETURN

 CONTAINS

!-------------------------------------------------------------------------
!NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1, GEOS/DAS!
!-------------------------------------------------------------------------
!BOP
!
! !INTERFACE:

 SUBROUTINE readClimVars(rc)

! !USES:

  IMPLICIT NONE

! !DESCRIPTION:
!
! Read, interpolate in time, and broadcast for time-dependent boundary conditions
!
! Input parameters:
!
! NONE
!
! Output parameters:
!
  INTEGER, INTENT(OUT) :: rc
!
! Restrictions:
!  ASSERT that the number of pressure layers in the dataset equals km.
!
! !REVISION HISTORY:
!  Nielsen     11 May 2012: First crack.
!  Nielsen     15 Jan 2015: Separate year-number shifting for GHG and ODS
!
! Notes
!  MAPL_ClimInterpFac returns m1=12 and m2=1 during the first half of January and
!  the second half of December.
!
!-----------------------------------------------------------------------

  TYPE(ESMF_Time) :: timeNow
  CHARACTER(LEN=ESMF_MAXSTR) :: Iam = "SC::readClimVars"

  INTEGER :: beginYear, ccyy, ccyyAdjusted
  INTEGER :: comm, dimid, info, unit, status
  INTEGER :: start(3), cnt(3), startco2(2), cntco2(2)
  INTEGER :: i, j, k, kLast
  INTEGER :: m1, m2, mo, n, numSSGmonths, numSO4months
  INTEGER :: YrAdj
  
  REAL :: fac, q, r, s

  REAL, ALLOCATABLE :: x(:,:,:)
  REAL, ALLOCATABLE :: zm1(:,:),zm2(:,:)
  REAL, ALLOCATABLE :: sm1(:),sm2(:)

  CHARACTER(LEN=ESMF_MAXSTR) :: fileName, varName

  CHARACTER(LEN=ESMF_MAXSTR) :: SSGName(25) = (/    "CO2     ", &
    "CH4     ", "N2O     ", "CFC11   ", "CFC12   ", "CFC113  ", & 
    "CFC114  ", "CFC115  ", "CCL4    ", "CH3CCL3 ", "HCFC22  ", &
    "HCFC141B", "HCFC142B", "H1211   ", "H1301   ", "H1202   ", &
    "H2402   ", "CH3BR   ", "CH3CL   ", "HFC23   ", "HFC32   ", &
    "HFC125  ", "HFC134A ", "HFC143A ", "HFC152A " /)

  CHARACTER(LEN=ESMF_MAXSTR) :: MBCName(2) = (/ "CO2_NOAA", &
                                                "SF6_NOAA" /)

  CHARACTER(LEN=ESMF_MAXSTR) :: zmClimName(6) = (/ "RAIN   ", "OH     ", &
                                                   "H2O    ", "O3     ", &
                                                   "SF6LOSS", "SO4SA  "/)

  rc = 0

  CALL ESMF_VMGet(vm, MPICOMMUNICATOR=comm, rc=status)
  VERIFY_(status)

  CALL ESMF_ClockGet(SCClock, CURRTIME=timeNow, RC=STATUS)
  VERIFY_(STATUS)
  CALL ESMF_TimeGet(timeNow, YY=ccyy, MM=mo, RC=STATUS)
  VERIFY_(STATUS)
  CALL MAPL_ClimInterpFac(SCClock, m1, m2, fac, RC=status)
  VERIFY_(status) 

  IF(m1 > m2 .AND. mo ==  1) m1 =  0
  IF(m1 > m2 .AND. mo == 12) m2 = 13

  fileName = TRIM(gcSC%climFileName)

#ifdef H5_HAVE_PARALLEL

  CALL MPI_Info_create(info, status)
  VERIFY_(status)
  CALL MPI_Info_set(info, "romio_cb_read", "automatic", status)
  VERIFY_(status)

#ifdef NETCDF_NEED_NF_MPIIO
  status = NF_OPEN_PAR(TRIM(fileName), IOR(NF_NOWRITE,NF_MPIIO), comm, info, unit)
#else
  status = NF_OPEN_PAR(TRIM(fileName), NF_NOWRITE, comm, info, unit)
#endif

#else

  IF(MAPL_AM_I_ROOT()) THEN 
   status = NF_OPEN(TRIM(fileName), NF_NOWRITE, unit)

#endif

   IF(status /= NF_NOERR) THEN
    PRINT *,"Error opening file ",TRIM(fileName), status
    PRINT *, NF_STRERROR(status)
    VERIFY_(status)
   END IF

#ifndef H5_HAVE_PARALLEL

  END IF ! MAPL_AM_I_ROOT

  CALL MAPL_CommsBcast(vm, gcSC%numClimLats, 1, 0, RC=status)
  VERIFY_(status)

  IF(MAPL_AM_I_ROOT()) THEN

#endif

   status = NF_INQ_DIMID(unit, 'time_SSG', dimid)
   IF(status /= NF_NOERR) THEN
    PRINT *,"Error inquiring dimid for time_SSG", status
    PRINT *, NF_STRERROR(status)
    VERIFY_(status)
   END IF

   status = NF_INQ_DIMLEN(unit, dimid, numSSGmonths)
   IF(status /= NF_NOERR) THEN
    PRINT *,"Error inquiring dimlen of time_SSG", status
    PRINT *, NF_STRERROR(status)
    VERIFY_(status)
   END IF

   status = NF_INQ_VARID(unit, 'time_SSG', n)
   IF(status /= NF_NOERR) THEN
    PRINT *,"Error getting varid for time_SSG", status
    PRINT *, NF_STRERROR(status)
    VERIFY_(status)
   END IF

   status = NF_GET_ATT_INT(unit, n, "begin_date", i)
   IF(status /= NF_NOERR) THEN
    PRINT *,"Error getting begin_date for time_SSG", status
    PRINT *, NF_STRERROR(status)
    VERIFY_(status)
   END IF

   beginYear = i/10000

   DO k = 1,gcSC%numSSGs

    IF(k <= gcSC%numGHGs) THEN
     YrAdj = gcSC%GHGYrAdj
    ELSE
     YrAdj = gcSC%ODSYrAdj
    ENDIF

    IF(YrAdj < -1000) THEN
     ccyyAdjusted = ABS(YrAdj)
    ELSE
     ccyyAdjusted = ccyy+YrAdj
    END IF

    i = m1+(ccyyAdjusted-beginYear)*12
    IF(i <            1) i = 1
    IF(i > numSSGmonths) i = numSSGmonths
    j = m2+(ccyyAdjusted-beginYear)*12
    IF(j <            1) j = 1
    IF(j > numSSGmonths) j = numSSGmonths
   
    IF(MAPL_AM_I_ROOT() .AND. k == gcSC%numGHGs ) THEN
     PRINT *,TRIM(Iam)//":"
     PRINT *," GHGs updated to ",ccyyAdjusted,".  Indicies: ",i,j
    END IF
    IF(MAPL_AM_I_ROOT() .AND. k == gcSC%numSSGs ) THEN
     PRINT *," ODSs updated to ",ccyyAdjusted,".  Indicies: ",i,j
    END IF

    varName = TRIM(SSGName(k))

    status = NF_INQ_VARID(unit, TRIM(varName), n)
    IF(status /= NF_NOERR) THEN
     PRINT *,"Error getting varid for ",TRIM(varName), status
     PRINT *, NF_STRERROR(status)
     VERIFY_(status)
    END IF

    status = NF_GET_VARA_REAL(unit, n, i, 1, r)
    IF(status /= NF_NOERR) THEN
     PRINT *,"Error getting first month for ",TRIM(varName), status
     PRINT *, NF_STRERROR(status)
     VERIFY_(status)
    END IF
    
    status = NF_GET_VARA_REAL(unit, n, j, 1, s)
    IF(status /= NF_NOERR) THEN
     PRINT *,"Error getting second month for ",TRIM(varName), status
     PRINT *, NF_STRERROR(status)
     VERIFY_(status)
    END IF

    gcSC%SSG(k) = r*fac+s*(1-fac)

   END DO

   status = NF_INQ_VARID(unit, 'lat', n)
   IF(status /= NF_NOERR) THEN
    PRINT *,"Error getting varid for lat", status
    PRINT *, NF_STRERROR(status)
    VERIFY_(status)
   END IF

   status = NF_GET_VAR_REAL(unit, n, gcSC%climLatRad)
   IF(status /= NF_NOERR) THEN
    PRINT *,"Error getting values for climLatRad", status
    PRINT *, NF_STRERROR(status)
    VERIFY_(status)
   END IF

   startco2(1) = 1
   cntco2(1) = gcSC%numClimLats
   cntco2(2) = 1

   ALLOCATE(sm1(gcSC%numClimLats), STAT=status)
   VERIFY_(status)
   ALLOCATE(sm2(gcSC%numClimLats), STAT=status)
   VERIFY_(status)

   YrAdj = gcSC%GHGYrAdj
   IF(YrAdj < -1000) THEN
    ccyyAdjusted = ABS(YrAdj)
   ELSE
    ccyyAdjusted = ccyy+YrAdj
   END IF
   
   i = m1+(ccyyAdjusted-beginYear)*12
   IF(i <            1) i = 1
   IF(i > numSSGmonths) i = numSSGmonths
   j = m2+(ccyyAdjusted-beginYear)*12
   IF(j <            1) j = 1
   IF(j > numSSGmonths) j = numSSGmonths

   DO k=1,2 

    varName = TRIM(MBCName(k))
    status = NF_INQ_VARID(unit, TRIM(varName), n)
    IF(status /= NF_NOERR) THEN
     PRINT *,"Error getting varid for ",TRIM(varName), status
     PRINT *, NF_STRERROR(status)
     VERIFY_(status)
    END IF
   
    startco2(2)=i
    status = NF_GET_VARA_REAL(unit, n, startco2, cntco2, sm1)
    IF(status /= NF_NOERR) THEN
     PRINT *,"Error getting values for first month of ",TRIM(varName), status
     PRINT *, NF_STRERROR(status)
     VERIFY_(status)
    END IF
   
    startco2(2)=j
    status = NF_GET_VARA_REAL(unit, n, startco2, cntco2, sm2)
    IF(status /= NF_NOERR) THEN
     PRINT *,"Error getting values for second month of ",TRIM(varName), status
     PRINT *, NF_STRERROR(status)
     VERIFY_(status)
    END IF

    IF(k==1) gcSC%CO2MBC(:)  = sm1(:)*fac+sm2(:)*(1.00-fac)
    IF(k==2) gcSC%SF6MBC(:)  = sm1(:)*fac+sm2(:)*(1.00-fac)

   END DO

   DEALLOCATE(sm1, STAT=status)
   VERIFY_(status)
   DEALLOCATE(sm2, STAT=status)
   VERIFY_(status)

   start(1) = 1
   start(2) = 1
   cnt(1) = gcSC%numClimLats
   cnt(2) = gcSC%km
   cnt(3) = 1

   ALLOCATE(zm1(gcSC%numClimLats,gcSC%km), STAT=status)
   VERIFY_(status)
   ALLOCATE(zm2(gcSC%numClimLats,gcSC%km), STAT=status)
   VERIFY_(status)

   status = NF_INQ_DIMID(unit, 'time_SO4', dimid)
   IF(status /= NF_NOERR) THEN
    PRINT *,"Error inquiring dimid for time_SO4", status
    PRINT *, NF_STRERROR(status)
    VERIFY_(status)
   END IF

   status = NF_INQ_DIMLEN(unit, dimid, numSO4months)
   IF(status /= NF_NOERR) THEN
    PRINT *,"Error inquiring dimlen of time_SO4", status
    PRINT *, NF_STRERROR(status)
    VERIFY_(status)
   END IF

   status = NF_INQ_VARID(unit, 'time_SO4', n)
   IF(status /= NF_NOERR) THEN
    PRINT *,"Error getting varid for time_SO4", status
    PRINT *, NF_STRERROR(status)
    VERIFY_(status)
   END IF

   status = NF_GET_ATT_INT(unit, n, "begin_date", i)
   IF(status /= NF_NOERR) THEN
    PRINT *,"Error getting begin_date for time_SO4", status
    PRINT *, NF_STRERROR(status)
    VERIFY_(status)
   END IF

   beginYear = i/10000
   
   IF(gcSC%SO4saYr == 0) THEN
    ccyyAdjusted = ccyy
   ELSE
    ccyyAdjusted = gcSC%SO4saYr
   END IF

   i = m1+(ccyyAdjusted-beginYear)*12
   IF(i <	     1) i = 1
   IF(i > numSO4months) i = numSO4months
   j = m2+(ccyyAdjusted-beginYear)*12
   IF(j <	     1) j = 1
   IF(j > numSO4months) j = numSO4months

   IF(MAPL_AM_I_ROOT()) THEN
    PRINT *," SO4  updated to ",ccyyAdjusted,".  Indicies: ",i,j
   END IF

   IF(m1 ==  0) m1 = 12
   IF(m2 == 13) m2 =  1

   DO k = 1,6

    varName = TRIM(zmClimName(k))

    status = NF_INQ_VARID(unit, TRIM(varName), n)
    IF(status /= NF_NOERR) THEN
     PRINT *,"Error getting varid for ",TRIM(varName), status
     PRINT *, NF_STRERROR(status)
     VERIFY_(status)
    END IF

    IF(k <= 5) THEN
     start(3) = m1
    ELSE
     start(3) = i
    END IF

    status = NF_GET_VARA_REAL(unit, n, start, cnt, zm1)
    IF(status /= NF_NOERR) THEN
     PRINT *,"Error getting values for first month of ",TRIM(varName), status
     PRINT *, NF_STRERROR(status)
     VERIFY_(status)
    END IF

    IF(k <= 5) THEN
     start(3) = m2
    ELSE
     start(3) = j
    END IF

    status = NF_GET_VARA_REAL(unit, n, start, cnt, zm2)
    IF(status /= NF_NOERR) THEN
     PRINT *,"Error getting values for second month of ",TRIM(varName), status
     PRINT *, NF_STRERROR(status)
     VERIFY_(status)
    END IF

    SELECT CASE (k)
     CASE (1)
      gcSC%climRain(:,:) = zm1(:,:)*fac+zm2(:,:)*(1.00-fac)
     CASE (2)
      gcSC%climOH(:,:) = zm1(:,:)*fac+zm2(:,:)*(1.00-fac)
     CASE (3)
      gcSC%climH2O(:,:) = zm1(:,:)*fac+zm2(:,:)*(1.00-fac)
     CASE (4)
      gcSC%climO3(:,:) = zm1(:,:)*fac+zm2(:,:)*(1.00-fac)
     CASE (5)
      gcSC%climSF6LOSS(:,:) = zm1(:,:)*fac+zm2(:,:)*(1.00-fac)
     CASE (6)
      gcSC%SO4sa(:,:) = zm1(:,:)*fac+zm2(:,:)*(1.00-fac)
     CASE DEFAULT
    END SELECT

   END DO

   DEALLOCATE(zm1, STAT=status)
   VERIFY_(status)
   DEALLOCATE(zm2, STAT=status)
   VERIFY_(status)

   status = NF_CLOSE(unit)
   VERIFY_(status)

#ifdef H5_HAVE_PARALLEL

   CALL MPI_Info_free(info, status)
   VERIFY_(status)

#else

  END IF ! MAPL_AM_I_ROOT

  CALL MAPL_CommsBcast(vm, gcSC%SSG, gcSC%numSSGs, 0, RC=status)
  VERIFY_(status)
  CALL MPI_Bcast(gcSC%CO2MBC, SIZE(gcSC%CO2MBC), MPI_REAL, 0, comm, status)
  VERIFY_(status)
  CALL MPI_Bcast(gcSC%SF6MBC, SIZE(gcSC%SF6MBC), MPI_REAL, 0, comm, status)
  VERIFY_(status)
  CALL MPI_Bcast(gcSC%climOH, SIZE(gcSC%climOH), MPI_REAL, 0, comm, status)
  VERIFY_(status)
  CALL MPI_Bcast(gcSC%climH2O, SIZE(gcSC%climH2O), MPI_REAL, 0, comm, status)
  VERIFY_(status)
  CALL MPI_Bcast(gcSC%climO3, SIZE(gcSC%climO3), MPI_REAL, 0, comm, status)
  VERIFY_(status)
  CALL MPI_Bcast(gcSC%climSF6LOSS, SIZE(gcSC%climSF6LOSS), MPI_REAL, 0, comm, status)
  VERIFY_(status)
  CALL MPI_Bcast(gcSC%SO4sa, SIZE(gcSC%SO4sa), MPI_REAL, 0, comm, status)
  VERIFY_(status)

#endif

 RETURN
 END SUBROUTINE readClimVars

 SUBROUTINE readEMSFactor(rc)

! !USES:

  IMPLICIT NONE

! !DESCRIPTION:
!
! Read annual emission factor (with respect to 1951) for the flux-based ODSs
!
! Input parameters:
!
! NONE
!
! Output parameters:
!
  INTEGER, INTENT(OUT) :: rc
!
! !REVISION HISTORY:
!  Liang     11 Sep 2015: First crack.
!
!-----------------------------------------------------------------------
  
  TYPE(ESMF_Time) :: timeNow
  CHARACTER(LEN=ESMF_MAXSTR) :: Iam = "SC::readEMSFactor"

  INTEGER :: ccyy, ccyyAdjusted
  INTEGER :: comm, dimid, info, i, k, n, status, unit
  REAL    :: r
 
  CHARACTER(LEN=ESMF_MAXSTR) :: fileName, varName
  
  CHARACTER(LEN=ESMF_MAXSTR) :: EMSName(6) = (/"CFC11_emission_factor  ", &
                                               "CFC12_emission_factor  ", &
                                               "CFC113_emission_factor ", &
                                               "HCFC22_emission_factor ", &
                                               "CH3CCL3_emission_factor", &
                                               "CCL4_emission_factor   "  /)

  rc = 0

  CALL ESMF_VMGet(vm, MPICOMMUNICATOR=comm, rc=status)
  VERIFY_(status)
  fileName = TRIM(gcSC%climFileName)

   CALL ESMF_ClockGet(SCClock, CURRTIME=timeNow, RC=STATUS)
   VERIFY_(STATUS)
   CALL ESMF_TimeGet(timeNow, YY=ccyy, RC=STATUS)
   VERIFY_(STATUS)

#ifdef H5_HAVE_PARALLEL

  CALL MPI_Info_create(info, status)
  VERIFY_(status)
  CALL MPI_Info_set(info, "romio_cb_read", "automatic", status)
  VERIFY_(status)

#ifdef NETCDF_NEED_NF_MPIIO
  status = NF_OPEN_PAR(TRIM(fileName), IOR(NF_NOWRITE,NF_MPIIO), comm, info, unit)
#else
  status = NF_OPEN_PAR(TRIM(fileName), NF_NOWRITE, comm, info, unit)
#endif

#else

  IF(MAPL_AM_I_ROOT(vm)) THEN
   status = NF_OPEN(TRIM(fileName), NF_NOWRITE, unit)

#endif

   IF(status /= NF_NOERR) THEN
    PRINT *,"Error opening file ",TRIM(fileName), status
    PRINT *, NF_STRERROR(status)
    VERIFY_(status)
   END IF

   status = NF_INQ_DIMID(unit, 'emiss_year', dimid)
   IF(status /= NF_NOERR) THEN
    PRINT *,"Error inquiring dimid for year", status
    PRINT *, NF_STRERROR(status)
    VERIFY_(status)
   END IF
   status = NF_INQ_DIMLEN(unit, dimid, n)
   IF(status /= NF_NOERR) THEN
    PRINT *,"Error inquiring dimlen of year", status
    PRINT *, NF_STRERROR(status)
    VERIFY_(status)
   END IF

   gcSC%nyemiss = n

#ifndef H5_HAVE_PARALLEL

  END IF ! MAPL_AM_I_ROOT

  CALL MAPL_CommsBcast(vm, gcSC%nyemiss, 1, 0, RC=status)
  VERIFY_(status)

  IF(MAPL_AM_I_ROOT()) THEN

#endif

   i = ccyy-1950
 
   DO k = 1,6
    
    varName = TRIM(EMSName(k))

    status = NF_INQ_VARID(unit, TRIM(varName), n)
    IF(status /= NF_NOERR) THEN
     PRINT *,"Error getting varid for ",TRIM(varName), status
     PRINT *, NF_STRERROR(status)
     VERIFY_(status)
    END IF

    status = NF_GET_VARA_REAL(unit, n, i, 1, r)
    IF(status /= NF_NOERR) THEN
     PRINT *,"Error getting first month for ",TRIM(varName), status
     PRINT *, NF_STRERROR(status)
     VERIFY_(status)
    END IF

    IF(k == 1) THEN
     gcSC%emfCFC11 = r
    END IF
    IF(k == 2) THEN
     gcSC%emfCFC12 = r
    END IF
    IF(k == 3) THEN
     gcSC%emfCFC113 = r
    END IF
    IF(k == 4) THEN
     gcSC%emfHCFC22 = r
    END IF
    IF(k == 5) THEN
     gcSC%emfMCF = r
    END IF
    IF(k == 6) THEN
     gcSC%emfCCL4 = r
    END IF

   END DO

   status = NF_CLOSE(unit)
   VERIFY_(status)

#ifdef H5_HAVE_PARALLEL

   CALL MPI_Info_free(info, status)
   VERIFY_(status)

#else

  END IF ! MAPL_AM_I_ROOT

  CALL MAPL_CommsBcast(vm, gcSC%emfCFC11, 1, 0, RC=status)
  VERIFY_(status)
  CALL MAPL_CommsBcast(vm, gcSC%emfCFC12, 1, 0, RC=status)
  VERIFY_(status)
  CALL MAPL_CommsBcast(vm, gcSC%emfCFC113, 1, 0, RC=status)
  VERIFY_(status)
  CALL MAPL_CommsBcast(vm, gcSC%emfHCFC22, 1, 0, RC=status)
  VERIFY_(status)
  CALL MAPL_CommsBcast(vm, gcSC%emfMCF, 1, 0, RC=status)
  VERIFY_(status)
  CALL MAPL_CommsBcast(vm, gcSC%emfCCL4, 1, 0, RC=status)
  VERIFY_(status)

#endif

 RETURN
 END SUBROUTINE readEMSFactor
 
 SUBROUTINE mapToGrid(field2D, field3D, name, rc)
  IMPLICIT NONE
  REAL, INTENT(IN) :: field2D(gcSC%numClimLats,gcSC%km)
  REAL, INTENT(OUT) :: field3D(gcSC%i1:gcSC%i2,gcSC%j1:gcSC%j2,gcSC%km)
  INTEGER, INTENT(OUT) :: rc
  CHARACTER(LEN=*), INTENT(IN) :: name
  CHARACTER(LEN=ESMF_MAXSTR), PARAMETER :: Iam = 'SC::mapToGrid'

  rc = 0

  DO k = 1,gcSC%km
   DO j = gcSC%j1,gcSC%j2
    DO i = gcSC%i1,gcSC%i2
     CALL MAPL_INTERP(field3D(i,j,k), gcSC%latrad(i,j), field2D(:,k), gcSC%climLatRad)
    END DO
   END DO
  END DO

  RETURN
 END SUBROUTINE mapToGrid
  
 SUBROUTINE mapToSurf(field1D, field2D, name, rc)
  IMPLICIT NONE
  REAL, INTENT(IN) :: field1D(gcSC%numClimLats)
  REAL, INTENT(OUT) :: field2D(gcSC%i1:gcSC%i2,gcSC%j1:gcSC%j2)
  INTEGER, INTENT(OUT) :: rc
  CHARACTER(LEN=*), INTENT(IN) :: name
  CHARACTER(LEN=ESMF_MAXSTR), PARAMETER :: Iam = 'SC::mapToSurf'

  rc = 0

  DO j = gcSC%j1,gcSC%j2
   DO i = gcSC%i1,gcSC%i2
    CALL MAPL_INTERP(field2D(i,j), gcSC%latrad(i,j), field1D(:), gcSC%climLatRad)
   END DO
  END DO

  RETURN
 END SUBROUTINE mapToSurf

 SUBROUTINE SSGSfcFlx(dZ,dV,m,mdry,n2o,ch4,f11,f12,f113,f114,f115,ccl4,ch3ccl3, &
                  ch3cl,ch3br,h1301,h1211,h1202,h2402,hcfc22,hcfc141b,hcfc142b, &
                                 chbr3,ch2br2,ch2brcl,chbrcl2,chbr2cl,co2,sf6,  &
                                    hfc23,hfc32,hfc125,hfc134a,hfc143a,hfc152a, &
                                                          mcfcol,ctccol,lwi,rc)

! ----------------------------------------------------------------------
! Update the surface source gases mixing ratios, and find the surface 
! mass flux per unit time [ kg s^-1 ] that is required to account for 
! the change in mixing ratio.
!
! The species in the argument list are in units of number density.
!
! m is number density in cm^{-3}, dV is the volume of the cell in 
! m^{3}, and tdt is the time step length (s).
! ----------------------------------------------------------------------

  IMPLICIT NONE
  REAL, INTENT(IN) :: dZ    ! Cell depth [cm]
  REAL, INTENT(IN) :: dV    ! Cell volume [m^{3}]
  REAL, INTENT(IN) :: m     ! Number density  [cm^{-3}]
  REAL, INTENT(IN) :: mdry  ! Dry Air number density  [cm^{-3}]
  REAL, INTENT(IN) :: lwi   ! Surface type index, 0 - Ocean; 1 - Land
  REAL(KIND=DBL), INTENT(IN) :: mcfcol, ctccol  ! MCF and CTC total column density
  INTEGER, INTENT(OUT) :: rc

  INTEGER, PARAMETER :: DBL=KIND(0.00D+00)
  REAL(KIND=DBL), INTENT(INOUT) :: n2o,ch4,f11,f12,f113,f114,f115,ccl4,ch3ccl3
  REAL(KIND=DBL), INTENT(INOUT) :: ch3cl,ch3br,h1301,h1211,h1202,h2402
  REAL(KIND=DBL), INTENT(INOUT) :: hcfc22,hcfc141b,hcfc142b,co2,sf6
  REAL(KIND=DBL), INTENT(INOUT) :: chbr3,ch2br2,ch2brcl,chbrcl2,chbr2cl
  REAL(KIND=DBL), INTENT(INOUT) :: hfc23,hfc32,hfc125,hfc134a,hfc143a,hfc152a

  REAL :: cx,cf

! Molecular weights of 12 species [kg kmol^-1]
! --------------------------------------------
  REAL, PARAMETER ::      mwtN2O =  44.008
  REAL, PARAMETER ::      mwtCH4 =  16.043
  REAL, PARAMETER ::      mwtF11 = 137.36
  REAL, PARAMETER ::      mwtF12 = 120.91
  REAL, PARAMETER ::     mwtF113 = 187.37
  REAL, PARAMETER ::     mwtF114 = 170.92
  REAL, PARAMETER ::     mwtF115 = 154.47
  REAL, PARAMETER ::   mwtHCFC22 =  86.47
  REAL, PARAMETER :: mwtHCFC141b = 116.95
  REAL, PARAMETER :: mwtHCFC142b = 100.50
  REAL, PARAMETER ::     mwtCCl4 = 153.823
  REAL, PARAMETER ::  mwtCH3CCl3 = 133.405
  REAL, PARAMETER ::    mwtCH3Cl =  50.488
  REAL, PARAMETER ::    mwtCH3Br =  94.944
  REAL, PARAMETER ::    mwtH1301 = 148.92
  REAL, PARAMETER ::    mwtH1211 = 164.36
  REAL, PARAMETER ::    mwtH1202 = 209.81 
  REAL, PARAMETER ::    mwtH2402 = 259.81
  REAL, PARAMETER ::    mwtCHBr3 = 252.73
  REAL, PARAMETER ::   mwtCH2Br2 = 173.84
  REAL, PARAMETER ::  mwtCH2BrCl = 129.38
  REAL, PARAMETER ::  mwtCHBrCl2 = 163.82
  REAL, PARAMETER ::  mwtCHBr2Cl = 208.28

  rc = 0

! This conversion constant assumes that
! the species number densities are in cm^-3.
! ------------------------------------------
  cx = dV*1.00E+06/(tdt*MAPL_AVOGAD)

! CO2 (ppmv) is used only in solverd2.F, and
! multiplication by number density takes place there.
! ---------------------------------------------------
  gcSC%co2     = gcSC%SSG(1) * 1.00E-06

! N20 and CH4 ppbv
! ----------------
  ch4 = gcSC%SSG( 2) * mdry * 1.00E-09
  n2o = gcSC%SSG( 3) * mdry * 1.00E-09

! Others pptv
! -----------

#ifndef REDUCED
   f114     = gcSC%SSG( 7) * mdry * 1.00E-12
   f115     = gcSC%SSG( 8) * mdry * 1.00E-12
#endif

  ch3cl     = gcSC%SSG(19) * mdry * 1.00E-12

! Add 5 pptv to the methyl bromide boundary
! condition when running the REDUCED mechanism.
! ---------------------------------------------
#ifdef REDUCED
  ch3br     = (gcSC%SSG(18)+5.00) * mdry * 1.00E-12
#else
  ch3br     = gcSC%SSG(18) * mdry * 1.00E-12
#endif

  h1301     = gcSC%SSG(15) * mdry * 1.00E-12
  h1211     = gcSC%SSG(14) * mdry * 1.00E-12

#ifndef REDUCED
  h1202     = gcSC%SSG(16) * mdry * 1.00E-12
  h2402     = gcSC%SSG(17) * mdry * 1.00E-12
  hcfc141b  = gcSC%SSG(12) * mdry * 1.00E-12
  hcfc142b  = gcSC%SSG(13) * mdry * 1.00E-12

  hfc23     = gcSC%SSG(20) * mdry * 1.00E-12
  hfc32     = gcSC%SSG(21) * mdry * 1.00E-12
  hfc125    = gcSC%SSG(22) * mdry * 1.00E-12
  hfc134a   = gcSC%SSG(23) * mdry * 1.00E-12
  hfc143a   = gcSC%SSG(24) * mdry * 1.00E-12
  hfc152a   = gcSC%SSG(25) * mdry * 1.00E-12

  ch2brcl   = 0.14 * mdry * 1.00E-12
  chbrcl2   = 0.10 * mdry * 1.00E-12
  chbr2cl   = 0.20 * mdry * 1.00E-12
#endif

  IF(gcSC%doFlux) THEN

! Implement surface loss at the ocean surface layer for MCF for
! ocean lifetime of 94 yr = 2.96438e+9 sec
! This needs to divided by a factor of 0.6648 to account for the ocean surface fraction

   IF (lwi == 0) THEN
    cf = mcfcol/dZ/(2.96438e+09*0.6648)
    ch3ccl3 = ch3ccl3 - cf*tdt
    mcfocnloss(i,j) = (ch3ccl3-bsc%qa(gcSC%iCH3CCl3)%data3d(i,j,km)*m)*cx*mwtCH3CCl3
   END IF

! Implement surface loss at the ocean surface layer for CCl4 for
! ocean lifetime of 81 yr = 2.55442e+9 sec
! This needs to divided by a factor of 0.6648 to account for the ocean surface fraction
! Qing Liang - 2013/10/24
! Set loss rate to 5.23498e+09 sec (equavalent to lifetime of 165.9 yrs)

   IF (lwi == 0) THEN
    cf = ctccol/dZ/(5.23498e+09*0.6648)
    ccl4 = ccl4 - cf*tdt
    ctcocnloss(i,j) = (ccl4-bsc%qa(gcSC%iCCl4)%data3d(i,j,km)*m)*cx*mwtCCl4
   END IF

! Implement surface loss at the land surface layer for CCl4 for
! land lifetime of 195 yr = 6.1495e+9 sec
! loss rate = 2.52287e+10 sec (equavalent to lifetime of 800 yrs)
! This needs to divided by a factor of 0.2938 to account for the land surface fraction

   IF (lwi == 1) THEN
    cf = ctccol/dZ/(2.52287e+10*0.2938)
    ccl4 = ccl4 - cf*tdt
    ctclndloss(i,j) = (ccl4-bsc%qa(gcSC%iCCl4)%data3d(i,j,km)*m)*cx*mwtCCl4
   END IF

! Add surface emissions

   cf = MAPL_AVOGAD*1E-4*tdt/dZ

   f11     = f11     + gcSC%CFC11EMISS(i,j)  *cf/mwtF11
   f12     = f12     + gcSC%CFC12EMISS(i,j)  *cf/mwtF12
   f113    = f113    + gcSC%CFC113EMISS(i,j) *cf/mwtF113
   ccl4    = ccl4    + gcSC%CCL4EMISS(i,j)   *cf/mwtCCl4
   ch3ccl3 = ch3ccl3 + gcSC%MCFEMISS(i,j)    *cf/mwtCH3CCl3
   hcfc22  = hcfc22  + gcSC%HCFC22EMISS(i,j) *cf/mwtHCFC22

#ifndef REDUCED
   chbr3   = chbr3   + gcSC%CHBR3EMISS(i,j)  *cf/mwtCHBr3
   ch2br2  = ch2br2  + gcSC%CH2BR2EMISS(i,j) *cf/mwtCH2Br2
#endif

  ELSE

   f11       = gcSC%SSG( 4) * mdry * 1.00E-12
   f12       = gcSC%SSG( 5) * mdry * 1.00E-12
   f113      = gcSC%SSG( 6) * mdry * 1.00E-12
   ccl4      = gcSC%SSG( 9) * mdry * 1.00E-12
   ch3ccl3   = gcSC%SSG(10) * mdry * 1.00E-12
   hcfc22    = gcSC%SSG(11) * mdry * 1.00E-12

#ifndef REDUCED
   chbr3     = 0.61 * mdry * 1.00E-12
   ch2br2    = 0.92 * mdry * 1.00E-12
#endif

  END IF

#ifndef REDUCED
  co2   = gcSC%CO2MBC2D(i,j)  * mdry * 1.00E-06
  sf6   = gcSC%SF6MBC2D(i,j)  * mdry * 1.00E-12    
#endif
 
! Find the flux per unit time
! ---------------------------
      n2oflx(i,j) = (      n2o-bsc%qa(     gcSC%iN2O)%data3d(i,j,km)*m)*cx*mwtN2O
      ch4flx(i,j) = (      ch4-bsc%qa(     gcSC%iCH4)%data3d(i,j,km)*m)*cx*mwtCH4
      f11flx(i,j) = (      f11-bsc%qa(     gcSC%iF11)%data3d(i,j,km)*m)*cx*mwtF11
      f12flx(i,j) = (      f12-bsc%qa(     gcSC%iF12)%data3d(i,j,km)*m)*cx*mwtF12

     f113flx(i,j) = (     f113-bsc%qa(    gcSC%iF113)%data3d(i,j,km)*m)*cx*mwtF113
   hcfc22flx(i,j) = (   hcfc22-bsc%qa(  gcSC%iHCFC22)%data3d(i,j,km)*m)*cx*mwtHCFC22

#ifndef REDUCED
 hcfc141bflx(i,j) = ( hcfc141b-bsc%qa(gcSC%iHCFC141b)%data3d(i,j,km)*m)*cx*mwtHCFC141b
 hcfc142bflx(i,j) = ( hcfc142b-bsc%qa(gcSC%iHCFC142b)%data3d(i,j,km)*m)*cx*mwtHCFC142b
#endif

     ccl4flx(i,j) = (     ccl4-bsc%qa(    gcSC%iCCl4)%data3d(i,j,km)*m)*cx*mwtCCl4
      mcfflx(i,j) = (  ch3ccl3-bsc%qa( gcSC%iCH3CCl3)%data3d(i,j,km)*m)*cx*mwtCH3CCl3

    ch3clflx(i,j) = (    ch3cl-bsc%qa(   gcSC%iCH3Cl)%data3d(i,j,km)*m)*cx*mwtCH3Cl
    ch3brflx(i,j) = (    ch3br-bsc%qa(   gcSC%iCH3Br)%data3d(i,j,km)*m)*cx*mwtCH3Br
    h1301flx(i,j) = (    h1301-bsc%qa(   gcSC%iH1301)%data3d(i,j,km)*m)*cx*mwtH1301
    h1211flx(i,j) = (    h1211-bsc%qa(   gcSC%iH1211)%data3d(i,j,km)*m)*cx*mwtH1211

#ifndef REDUCED
    h1202flx(i,j) = (    h1202-bsc%qa(   gcSC%iH1202)%data3d(i,j,km)*m)*cx*mwtH1202
    h2402flx(i,j) = (    h2402-bsc%qa(   gcSC%iH2402)%data3d(i,j,km)*m)*cx*mwtH2402
#endif
  
  RETURN
 END SUBROUTINE SSGSfcFlx

 SUBROUTINE swapSpecies(direction,rc)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: direction
  INTEGER, INTENT(OUT) :: rc
  CHARACTER(LEN=ESMF_MAXSTR) :: Iam = "SC::swapSpecies"
  REAL :: m, rm

  m = numDens(k)
  rm = 1.00/m
  rc = 0

  SELECT CASE (direction)
  
   CASE(1)
 
           ox = bsc%qa(      gcSC%iOx)%data3d(i,j,kRev)*m
          nox = bsc%qa(     gcSC%iNOx)%data3d(i,j,kRev)*m
         hno3 = bsc%qa(    gcSC%iHNO3)%data3d(i,j,kRev)*m
         n2o5 = bsc%qa(    gcSC%iN2O5)%data3d(i,j,kRev)*m
       ho2no2 = bsc%qa(  gcSC%iHO2NO2)%data3d(i,j,kRev)*m
       clono2 = bsc%qa(  gcSC%iClONO2)%data3d(i,j,kRev)*m
          clx = bsc%qa(     gcSC%iClx)%data3d(i,j,kRev)*m
          hcl = bsc%qa(     gcSC%iHCl)%data3d(i,j,kRev)*m
         hocl = bsc%qa(    gcSC%iHOCl)%data3d(i,j,kRev)*m
         h2o2 = bsc%qa(    gcSC%iH2O2)%data3d(i,j,kRev)*m
          brx = bsc%qa(     gcSC%iBrx)%data3d(i,j,kRev)*m
          n2o = bsc%qa(     gcSC%iN2O)%data3d(i,j,kRev)*m
          cl2 = bsc%qa(     gcSC%iCl2)%data3d(i,j,kRev)*m
         oclo = bsc%qa(    gcSC%iOClO)%data3d(i,j,kRev)*m
         brcl = bsc%qa(    gcSC%iBrCl)%data3d(i,j,kRev)*m
          hbr = bsc%qa(     gcSC%iHBr)%data3d(i,j,kRev)*m
       brono2 = bsc%qa(  gcSC%iBrONO2)%data3d(i,j,kRev)*m
          ch4 = bsc%qa(     gcSC%iCH4)%data3d(i,j,kRev)*m
         hobr = bsc%qa(    gcSC%iHOBr)%data3d(i,j,kRev)*m
       ch3ooh = bsc%qa(  gcSC%iCH3OOH)%data3d(i,j,kRev)*m
           co = bsc%qa(      gcSC%iCO)%data3d(i,j,kRev)*m
        hno3c = bsc%qa(   gcSC%iHNO3c)%data3d(i,j,kRev)*m
          f11 = bsc%qa(     gcSC%iF11)%data3d(i,j,kRev)*m
          f12 = bsc%qa(     gcSC%iF12)%data3d(i,j,kRev)*m
         f113 = bsc%qa(    gcSC%iF113)%data3d(i,j,kRev)*m
         f114 = bsc%qa(    gcSC%iF114)%data3d(i,j,kRev)*m
         f115 = bsc%qa(    gcSC%iF115)%data3d(i,j,kRev)*m
       hcfc22 = bsc%qa(  gcSC%iHCFC22)%data3d(i,j,kRev)*m
     hcfc141b = bsc%qa(gcSC%iHCFC141b)%data3d(i,j,kRev)*m
     hcfc142b = bsc%qa(gcSC%iHCFC142b)%data3d(i,j,kRev)*m
         ccl4 = bsc%qa(    gcSC%iCCl4)%data3d(i,j,kRev)*m
      ch3ccl3 = bsc%qa( gcSC%iCH3CCl3)%data3d(i,j,kRev)*m
        ch3cl = bsc%qa(   gcSC%iCH3Cl)%data3d(i,j,kRev)*m
        ch3br = bsc%qa(   gcSC%iCH3Br)%data3d(i,j,kRev)*m
        h1301 = bsc%qa(   gcSC%iH1301)%data3d(i,j,kRev)*m
        h1211 = bsc%qa(   gcSC%iH1211)%data3d(i,j,kRev)*m
        h1202 = bsc%qa(   gcSC%iH1202)%data3d(i,j,kRev)*m
        h2402 = bsc%qa(   gcSC%iH2402)%data3d(i,j,kRev)*m
        chbr3 = bsc%qa(   gcSC%iCHBr3)%data3d(i,j,kRev)*m
       ch2br2 = bsc%qa(  gcSC%iCH2Br2)%data3d(i,j,kRev)*m
      ch2brcl = bsc%qa( gcSC%iCH2BrCl)%data3d(i,j,kRev)*m
      chbrcl2 = bsc%qa( gcSC%iCHBrCl2)%data3d(i,j,kRev)*m
      chbr2cl = bsc%qa( gcSC%iCHBr2Cl)%data3d(i,j,kRev)*m
        hfc23 = bsc%qa(   gcSC%iHFC23)%data3d(i,j,kRev)*m
        hfc32 = bsc%qa(   gcSC%iHFC32)%data3d(i,j,kRev)*m
       hfc125 = bsc%qa(  gcSC%iHFC125)%data3d(i,j,kRev)*m
      hfc134a = bsc%qa( gcSC%iHFC134a)%data3d(i,j,kRev)*m
      hfc143a = bsc%qa( gcSC%iHFC143a)%data3d(i,j,kRev)*m
      hfc152a = bsc%qa( gcSC%iHFC152a)%data3d(i,j,kRev)*m
          co2 = bsc%qa(     gcSC%iCO2)%data3d(i,j,kRev)*m
          sf6 = bsc%qa(     gcSC%iSF6)%data3d(i,j,kRev)*m
           o3 = bsc%qa(      gcSC%iO3)%data3d(i,j,kRev)*m
          o3p = bsc%qa(     gcSC%iO3p)%data3d(i,j,kRev)*m
          o1d = bsc%qa(     gcSC%iO1d)%data3d(i,j,kRev)*m
            n = bsc%qa(       gcSC%iN)%data3d(i,j,kRev)*m
           no = bsc%qa(      gcSC%iNO)%data3d(i,j,kRev)*m
          no2 = bsc%qa(     gcSC%iNO2)%data3d(i,j,kRev)*m
          no3 = bsc%qa(     gcSC%iNO3)%data3d(i,j,kRev)*m
            h = bsc%qa(       gcSC%iH)%data3d(i,j,kRev)*m
           oh = bsc%qa(      gcSC%iOH)%data3d(i,j,kRev)*m
          ho2 = bsc%qa(     gcSC%iHO2)%data3d(i,j,kRev)*m
           cl = bsc%qa(      gcSC%iCl)%data3d(i,j,kRev)*m
          clo = bsc%qa(     gcSC%iClO)%data3d(i,j,kRev)*m
          bro = bsc%qa(     gcSC%iBrO)%data3d(i,j,kRev)*m
           br = bsc%qa(      gcSC%iBr)%data3d(i,j,kRev)*m
        cl2o2 = bsc%qa(   gcSC%iCl2O2)%data3d(i,j,kRev)*m
         ch2o = bsc%qa(    gcSC%iCH2O)%data3d(i,j,kRev)*m
        ch3o2 = bsc%qa(   gcSC%iCH3O2)%data3d(i,j,kRev)*m
    box_ro3ox = bsc%qa(   gcSC%irO3Ox)%data3d(i,j,kRev)

   CASE(-1)

     bsc%qa(      gcSC%iOx)%data3d(i,j,kRev) =       ox*rm
     bsc%qa(     gcSC%iNOx)%data3d(i,j,kRev) =      nox*rm
     bsc%qa(    gcSC%iHNO3)%data3d(i,j,kRev) =     hno3*rm
     bsc%qa(    gcSC%iN2O5)%data3d(i,j,kRev) =     n2o5*rm
     bsc%qa(  gcSC%iHO2NO2)%data3d(i,j,kRev) =   ho2no2*rm
     bsc%qa(  gcSC%iClONO2)%data3d(i,j,kRev) =   clono2*rm
     bsc%qa(     gcSC%iClx)%data3d(i,j,kRev) =      clx*rm
     bsc%qa(     gcSC%iHCl)%data3d(i,j,kRev) =      hcl*rm
     bsc%qa(    gcSC%iHOCl)%data3d(i,j,kRev) =     hocl*rm
     bsc%qa(    gcSC%iH2O2)%data3d(i,j,kRev) =     h2o2*rm
     bsc%qa(     gcSC%iBrx)%data3d(i,j,kRev) =      brx*rm
     bsc%qa(     gcSC%iN2O)%data3d(i,j,kRev) =      n2o*rm
     bsc%qa(     gcSC%iCl2)%data3d(i,j,kRev) =      cl2*rm
     bsc%qa(    gcSC%iOClO)%data3d(i,j,kRev) =     oclo*rm
     bsc%qa(    gcSC%iBrCl)%data3d(i,j,kRev) =     brcl*rm
     bsc%qa(     gcSC%iHBr)%data3d(i,j,kRev) =      hbr*rm
     bsc%qa(  gcSC%iBrONO2)%data3d(i,j,kRev) =   brono2*rm
     bsc%qa(     gcSC%iCH4)%data3d(i,j,kRev) =      ch4*rm
     bsc%qa(    gcSC%iHOBr)%data3d(i,j,kRev) =     hobr*rm
     bsc%qa(  gcSC%iCH3OOH)%data3d(i,j,kRev) =   ch3ooh*rm
     bsc%qa(      gcSC%iCO)%data3d(i,j,kRev) =       co*rm
     bsc%qa(   gcSC%iHNO3c)%data3d(i,j,kRev) =    hno3c*rm
     bsc%qa(     gcSC%iF11)%data3d(i,j,kRev) =      f11*rm
     bsc%qa(     gcSC%iF12)%data3d(i,j,kRev) =      f12*rm
     bsc%qa(    gcSC%iF113)%data3d(i,j,kRev) =     f113*rm
     bsc%qa(    gcSC%iF114)%data3d(i,j,kRev) =     f114*rm
     bsc%qa(    gcSC%iF115)%data3d(i,j,kRev) =     f115*rm
     bsc%qa(  gcSC%iHCFC22)%data3d(i,j,kRev) =   hcfc22*rm
     bsc%qa(gcSC%iHCFC141b)%data3d(i,j,kRev) = hcfc141b*rm
     bsc%qa(gcSC%iHCFC142b)%data3d(i,j,kRev) = hcfc142b*rm
     bsc%qa(    gcSC%iCCl4)%data3d(i,j,kRev) =     ccl4*rm
     bsc%qa( gcSC%iCH3CCl3)%data3d(i,j,kRev) =  ch3ccl3*rm
     bsc%qa(   gcSC%iCH3Cl)%data3d(i,j,kRev) =    ch3cl*rm
     bsc%qa(   gcSC%iCH3Br)%data3d(i,j,kRev) =    ch3br*rm
     bsc%qa(   gcSC%iH1301)%data3d(i,j,kRev) =    h1301*rm
     bsc%qa(   gcSC%iH1211)%data3d(i,j,kRev) =    h1211*rm
     bsc%qa(   gcSC%iH1202)%data3d(i,j,kRev) =    h1202*rm
     bsc%qa(   gcSC%iH2402)%data3d(i,j,kRev) =    h2402*rm
     bsc%qa(   gcSC%iCHBr3)%data3d(i,j,kRev) =    chbr3*rm
     bsc%qa(  gcSC%iCH2Br2)%data3d(i,j,kRev) =   ch2br2*rm
     bsc%qa( gcSC%iCH2BrCl)%data3d(i,j,kRev) =  ch2brcl*rm
     bsc%qa( gcSC%iCHBrCl2)%data3d(i,j,kRev) =  chbrcl2*rm
     bsc%qa( gcSC%iCHBr2Cl)%data3d(i,j,kRev) =  chbr2cl*rm
     bsc%qa(   gcSC%iHFC23)%data3d(i,j,kRev) =    hfc23*rm
     bsc%qa(   gcSC%iHFC32)%data3d(i,j,kRev) =    hfc32*rm
     bsc%qa(  gcSC%iHFC125)%data3d(i,j,kRev) =   hfc125*rm
     bsc%qa( gcSC%iHFC134a)%data3d(i,j,kRev) =  hfc134a*rm
     bsc%qa( gcSC%iHFC143a)%data3d(i,j,kRev) =  hfc143a*rm
     bsc%qa( gcSC%iHFC152a)%data3d(i,j,kRev) =  hfc152a*rm
     bsc%qa(     gcSC%iCO2)%data3d(i,j,kRev) =      co2*rm
     bsc%qa(     gcSC%iSF6)%data3d(i,j,kRev) =      sf6*rm
     bsc%qa(      gcSC%iO3)%data3d(i,j,kRev) =       o3*rm
     bsc%qa(     gcSC%iO3p)%data3d(i,j,kRev) =      o3p*rm
     bsc%qa(     gcSC%iO1d)%data3d(i,j,kRev) =      o1d*rm
     bsc%qa(       gcSC%iN)%data3d(i,j,kRev) =        n*rm
     bsc%qa(      gcSC%iNO)%data3d(i,j,kRev) =       no*rm
     bsc%qa(     gcSC%iNO2)%data3d(i,j,kRev) =      no2*rm
     bsc%qa(     gcSC%iNO3)%data3d(i,j,kRev) =      no3*rm
     bsc%qa(       gcSC%iH)%data3d(i,j,kRev) =        h*rm
     bsc%qa(      gcSC%iOH)%data3d(i,j,kRev) =       oh*rm
     bsc%qa(     gcSC%iHO2)%data3d(i,j,kRev) =      ho2*rm
     bsc%qa(      gcSC%iCl)%data3d(i,j,kRev) =       cl*rm
     bsc%qa(     gcSC%iClO)%data3d(i,j,kRev) =      clo*rm
     bsc%qa(     gcSC%iBrO)%data3d(i,j,kRev) =      bro*rm
     bsc%qa(      gcSC%iBr)%data3d(i,j,kRev) =       br*rm
     bsc%qa(   gcSC%iCl2O2)%data3d(i,j,kRev) =    cl2o2*rm
     bsc%qa(    gcSC%iCH2O)%data3d(i,j,kRev) =     ch2o*rm
     bsc%qa(   gcSC%iCH3O2)%data3d(i,j,kRev) =    ch3o2*rm
     bsc%qa(   gcSC%irO3Ox)%data3d(i,j,kRev) =   box_ro3ox

   CASE DEFAULT

  END SELECT

 RETURN
 END SUBROUTINE swapSpecies

 SUBROUTINE chemStage1(rc)
  INTEGER, INTENT(OUT) :: rc
  CHARACTER(LEN=ESMF_MAXSTR) :: Iam = "SC::chemStage1"
  rc = 0

! Initialize ratios, production and loss rates and estimated species to zero.
! ---------------------------------------------------------------------------
  CALL localZeroSet(ro3po3,ro1do3,ro3ox,rnno,rnono2,rno3no2,rno2nox,rclclo, &
      rcloclx,rocloclo,rcl2o2clo,rbrobrx,rbrnbrx,po3,lo3,ph2o2,lh2o2,phno3, &
      lhno3,phno3het,pn2o5,ln2o5,ln2o5het,pho2no2,lho2no2,pnox,lnox,lnoxsq, &
        pnoxcl,lnoxcl,pclx,lclx,pclono2,lclono2,lclono2het,phcl,lhcl,phocl, &
        lhocl,pbrono2,lbrono2het,pbrx,lbrx,phbr,lhbr,phobr,lhobr,pno2,lno2, &
        pno3,lno3,poclo,pnoxa,lnoxa,pcl2,pbrcl,lo3hox,lo3nox,lo3oxsq,lo3cl, &
       lo3cly,lo3oh,lo3brx,lo3chx,o3pe,o3e,oxe,ne,noe,no2e,no3e,noxe,n2o5e, &
         h2o2e,hno3e,ho2no2e,hcle,hocle,cle,cloe,ocloe,cl2o2e,clono2e,clxe, &
                                          bre,broe,brono2e,brxe,hbre,hobre, &
                      rmedsts(k),rmednat(k),rmedice(k),denssts(k),vfall(k))

! Condensed water number density. QCTOT is mass fraction.
! -------------------------------------------------------
  h2oc = qctot(i,j,kRev)*numDens(k)*MAPL_AIRMW/MAPL_H2OMW

! Polar stratospheric clouds
! --------------------------
  DoingPSCs: IF(gcSC%doPSCs) THEN

   ic = gcSC%irO3Ox
   tropp = bsc%qa(ic)%data3d(i,j,km)*0.01

   CALL pscs(k,km,gcSC%latRad(i,j),tropp,midPress,edgePress, &
  	          kel,climH2OProfile,numDens,SO4saProfile,sanat,saice, &
  	     rmedsts,rmednat,rmedice,denssts,hno3,hno3c,h2o,h2oc,gcSC)

! Check for HNO3COND above chosen limit
! -------------------------------------
   qmax = 1.00E+09*hno3c/numDens(k)
   IF(qmax > gcSC%HNO3Ice_MAX) THEN
    PRINT *,TRIM(Iam)//": Found HNO3COND above limit: ",qmax," ppbv"
    status = 1
    VERIFY_(status)
   END IF

! Export states for NAT and ice surface area density [m^{-1}]
! -----------------------------------------------------------
   natsad(i,j,kRev) = 1.00E+02*sanat
   icesad(i,j,kRev) = 1.00E+02*saice

! Khet3d resets heterogeneous reaction rates 110-115.
! ---------------------------------------------------
   IF(k >= k1Strat) THEN
    CALL khet3d(k,km,gcSC%numreacs,gcSC%parts,SO4saProfile,sanat, &
		saice,numDens,kel,h2o,hcl,clono2,ak,status)
    VERIFY_(status)
   END IF

! Ad-hoc reaction rate adjustments.
! ---------------------------------
   IF(ak(k,111)* clono2 > 5.00E-04) ak(k,111) = 5.00E-04/clono2
   IF(ak(k,114)*   hocl > 5.00E-04) ak(k,114) = 5.00E-04/hocl
   IF(hcl/numDens(k) < 1.00E-12) THEN
    ak(k,111) = 0.00
    ak(k,114) = 0.00
   END IF

  END IF DoingPSCs

  IF(ak(k,31)*numDens(k)     > 5.00E-04) ak(k,31) = 5.00E-04/ numDens(k)
  IF(ak(k,32)*numDens(k)*clo > 5.00E-04) ak(k,32) = 5.00E-04/(numDens(k)*clo)
  IF(ak(k,32)*numDens(k)*no2 > 5.00E-04) ak(k,32) = 5.00E-04/(numDens(k)*no2)

! Constrain tropospheric OH
! -------------------------
  IF(k < k1Strat) oh = climOHProfile(k)

  CALL calc_qqjk(k,k1Strat,km,dt,daytime,aj,ak,qqj,qqk,oh,o1d,cl,ch4,n2o,f11,f12, &
            f113,f114,f115,ch3ccl3,ccl4,ch3cl,ch3br,chbr3,ch2br2,ch2brcl,chbrcl2, &
     chbr2cl,h1301,h1211,h1202,h2402,hcfc22,hcfc141b,hcfc142b,hfc23,hfc32,hfc125, & 
               hfc134a,hfc143a,hfc152a,gcSC%levels_cal,gcSC%numreacs,gcSC%numphoto) 

  IF(gcSC%spinup .AND. passNumber == 1) CALL sproic(k,km,midPress,numDens,daytime, &
            aj,ak,br,bro,brono2,brx,cl,cl2o2,clo,clx,h,h2o,hno3,ho2,n,n2o5,no,no2, &
                no3,nox,o1d,o3,o3p,oclo,oh,ox,gcSC%numphoto,gcSC%numreacs,gcSC%o2)

  CALL part(k,k1Strat,km,daytime,numDens,aj,ak,bro,brx,ch2o,ch3o2,ch4,clo,clono2, &
             h,hbr,hno3,ho2,hobr,n,n2o5,no,no2,no3,nox,o3,o3p,oh,rbrnbrx,rbrobrx, &
                   rcl2o2clo,rclclo,rcloclx,rnno,rno2nox,rno3no2,rnono2,rocloclo, &
                             gcSC%levels_cal,gcSC%numreacs,gcSC%numphoto,gcSC%o2) 

  CALL chxpart(k,k1Strat,km,daytime,numDens,aj,ak,brx,ch2o,ch3,ch3br,ch3o,ch3o2, &
           ch3ooh,ch4,cho,clx,ho2,nox,o1d,o3p,oh,rbrnbrx,rbrobrx,rclclo,rcloclx, &
             rno2nox,rnono2,gcSC%levels_cal,gcSC%numreacs,gcSC%numphoto,gcSC%o2)

  CALL hoxpart(k,k1Strat,km,daytime,numDens,tropo,climOHProfile,ak,aj,brx,ch2o,ch3o, &
              ch3o2,ch3ooh,ch4,cho,clono2,clx,co,h,hbr,hcl,h2o,h2o2,hno3,ho2,ho2no2, &
            hobr,hocl,nox,o1d,o3,o3p,oh,rbrnbrx,rbrobrx,rclclo,rcloclx,rnno,rno2nox, &
             rnono2,gcSC%levels_cal,gcSC%niters,gcSC%numreacs,gcSC%numphoto,gcSC%o2)

  CALL oxpandl(k,k1Strat,km,tropo,daytime,numDens,aj,ak,brx,ch2o,ch3,ch3ooh,ch4, &
      cl2o2,clono2,clx,h,h2o,h2o2,hbr,hcl,hno3,ho2,ho2no2,hocl,lh2o2,lo3,lo3brx, &
      lo3chx,lo3cl,lo3cly,lo3hox,lo3nox,lo3oh,lo3oxsq,n2o,n2o5,nox,o1d,o3,oh,ox, &
   ph2o2,po3,climRainProfile,rbrobrx,rclclo,rcloclx,rnno,rno2nox,rno3no2,rnono2, &
        ro1do3,ro3ox,ro3po3,gcSC%levels_cal,gcSC%numreacs,gcSC%numphoto,gcSC%o2)

  CALL brpandl(dt,k,k1Strat,km,daytime,ak,aj,brono2,brx,ch2o,ch3br,h1301,h1211, &
                 h1202,h2402,chbr3,ch2br2,ch2brcl,chbrcl2,chbr2cl,hbr,ho2,hobr, &
             lbrono2het,lbrx,lhbr,lhobr,o1d,oh,pbrx,phbr,phobr,climRainProfile, &
           rbrnbrx,rbrobrx,gcSC%levels_cal,gcSC%numreacs,gcSC%numphoto,gcSC%o2)

  CALL ncpandl(dt,k,k1Strat,km,daytime,numDens,ak,aj,ccl4,ch3ccl3,ch3cl,ch4,clono2, &
                clx,f11,f113,f12,f114,f115,h2o,h2o2,hcfc22,hcfc141b,hcfc142b,h1211, &
               ch2brcl,chbrcl2,chbr2cl,hcl,hno3,ho2,ho2no2,hocl,lbrono2het,lclono2, &
        lclono2het,lclx,lhcl,lhno3,lho2no2,lhocl,ln2o5,ln2o5het,lnox,lnoxcl,lnoxsq, &
         n,n2o,n2o5,no,no2,nox,o1d,o3p,oh,pbrono2,pclono2,pclx,phcl,phno3,phno3het, &
            pho2no2,phocl,pn2o5,pnox,pnoxcl,climRainProfile,rclclo,rcloclx,rno2nox, &
            rno3no2,gcSC%levels_cal,gcSC%numreacs,gcSC%numphoto,gcSC%o2,gcrProfile)

  CALL solverest(dt,requiredDt,examineDt,k,k1Strat,km,daytime,tropo,numDens,aj,ak,br, &
       brcl,bre,bro,broe,brono2,brono2e,brx,brxe,ccl4,ch3ccl3,ch3cl,cl,cle,cl2,cl2o2, &
           cl2o2e,clo,cloe,clono2,clono2e,clx,clxe,f11,f113,f12,f114,f115,h2o2,h2o2e, &
                     hbr,hbre,hcfc22,hcfc141b,hcfc142b,h1211,ch2brcl,chbr2cl,chbrcl2, &
             hcl,hcle,hno3,hno3e,ho2,ho2no2,ho2no2e,hobr,hobre,hocl,hocle,lbrono2het, &
        lbrx,lclono2,lclono2het,lclx,lh2o2,lhbr,lhcl,lhno3,lho2no2,lhobr,lhocl,ln2o5, &
      ln2o5het,lno2,lno3,lnox,lnoxa,lnoxcl,lnoxsq,lo3,n2o,n2o5,n2o5e,ne,no2,no2e,no3, &
    no3e,no,noe,nox,noxe,o1d,o3,o3e,o3p,o3pe,oclo,ocloe,oh,ox,oxe,pbrcl,pbrono2,pbrx, &
          pcl2,pclono2,pclx,ph2o2,phbr,phcl,pho2no2,phocl,phno3,phno3het,phobr,pn2o5, &
     pno2,pno3,pnox,pnoxa,pnoxcl,po3,poclo,climRainProfile,rbrnbrx,rbrobrx,rcl2o2clo, &
          rclclo,rcloclx,rnno,rnono2,rno2nox,rno3no2,ro3ox,ro3po3,rocloclo,box_ro3ox, &
                      gcSC%levels_cal,gcSC%numreacs,gcSC%numphoto,gcrProfile,gcSC%o2)

  IF(examineDt <= 3) CALL verifyChemDt(examineDt,requiredPasses,tdt,dt,requiredDt,startOver)

  RETURN
 END SUBROUTINE chemStage1

 SUBROUTINE chemStage2(rc)
  INTEGER, INTENT(OUT) :: rc
  CHARACTER(LEN=ESMF_MAXSTR) :: Iam = "SC::chemStage2"
  rc = 0

  CALL partest(k,k1Strat,km,daytime,numDens,aj,ak,broe,brxe,ch2o,ch3o2,ch4,cloe,clono2e, &
         h,hbre,hno3e,ho2,hobre,ne,n2o5e,noe,no2e,no3e,noxe,o3e,o3pe,oh,rbrnbrx,rbrobrx, &
     		          rcl2o2clo,rclclo,rcloclx,rnno,rno2nox,rno3no2,rnono2,rocloclo, &
                                    gcSC%levels_cal,gcSC%numreacs,gcSC%numphoto,gcSC%o2)

  CALL hoxpartest(k,k1Strat,km,daytime,numDens,tropo,climOHProfile,ak,aj,brxe,ch2o,ch3o, &
           ch3o2,ch3ooh,ch4,cho,clono2e,clxe,co,h,hbre,hcle,h2o,h2o2e,hno3e,ho2,ho2no2e, &
             hobre,hocle,noxe,o1d,o3,o3p,oh,rbrnbrx,rbrobrx,rclclo,rcloclx,rnno,rno2nox, &
                 rnono2,gcSC%levels_cal,gcSC%niters,gcSC%numreacs,gcSC%numphoto,gcSC%o2)

  CALL oxpandlest(k,k1Strat,km,tropo,daytime,numDens,aj,ak,brxe,ch2o,ch3,ch3ooh,ch4, &
 cl2o2e,clono2e,clxe,h,h2o,h2o2e,hbre,hcle,hno3e,ho2,ho2no2e,hocle,lh2o2,lo3,lo3brx, &
           lo3chx,lo3cl,lo3cly,lo3hox,lo3nox,lo3oh,lo3oxsq,n2o,n2o5e,noxe,oh,ox,oxe, &
       ph2o2,po3,climRainProfile,rbrobrx,rclclo,rcloclx,rnno,rno2nox,rno3no2,rnono2, &
            ro1do3,ro3ox,ro3po3,gcSC%levels_cal,gcSC%numreacs,gcSC%numphoto,gcSC%o2)

  CALL brpandlest(dt,k,k1Strat,km,daytime,ak,aj,brono2,brx,brxe,ch2o,ch3br,h1301,h1211, &
                  h1202,h2402,chbr3,ch2br2,ch2brcl,chbrcl2,chbr2cl,hbre,ho2,hobr,hobre, &
                     lbrono2het,lbrx,lhbr,lhobr,o3e,oh,pbrx,phbr,phobr,climRainProfile, &
            rbrnbrx,rbrobrx,ro1do3,gcSC%levels_cal,gcSC%numreacs,gcSC%numphoto,gcSC%o2)

  CALL ncpandlest(dt,k,k1Strat,km,daytime,numDens,ak,aj,brxe,ccl4,ch3ccl3,ch3cl,ch4,clono2e, &
                       clxe,f11,f113,f12,f114,f115,h2o,h2o2e,hcfc22,hcfc141b,hcfc142b,h1211, &
                    ch2brcl,chbrcl2,chbr2cl,hcl,hcle,hno3e,ho2,ho2no2,ho2no2e,hocle,lclono2, &
                 lclono2het,lclx,lhcl,lhno3,lho2no2,lhocl,ln2o5,ln2o5het,lnox,lnoxcl,lnoxsq, &
               n2o,n2o5,n2o5e,ne,noe,no2e,noxe,o3e,o3pe,oh,pclono2,pclx,phcl,phno3,phno3het, &
             pho2no2,phocl,pn2o5,pnox,pnoxcl,climRainProfile,rbrnbrx,rclclo,rcloclx,rno2nox, &
              rno3no2,ro1do3,gcSC%levels_cal,gcSC%numreacs,gcSC%numphoto,gcSC%o2,gcrProfile)

       CALL solverd2(dt,requiredDt,examineDt,k,k1Strat,km,daytime,tropo,numDens,aj,ak,br, &
           brcl,bro,brono2,brx,ccl4,ch2o,ch3br,ch3ccl3,ch3cl,ch3o2,ch3ooh,ch4,cho,cl,cl2, &
                   chbr3,ch2br2,ch2brcl,chbrcl2,chbr2cl,co2,sf6,cl2o2,clo,clono2,clono2e, &
               clx,clxe,co,f11,f113,f12,f114,f115,h,h1211,h1301,h1202,h2402,h2o,h2o2,hbr, &
            hcfc22,hcfc141b,hcfc142b,hfc23,hfc32,hfc125,hfc134a,hfc143a,hfc152a,hcl,hcle, &
      hno3,ho2,ho2no2,hobr,hocl,hocle,lbrono2het,lbrx,lclono2,lclono2het,lclx,lh2o2,lhbr, &
               lhcl,lhno3,lho2no2,lhobr,lhocl,ln2o5,ln2o5het,lno2,lno3,lnox,lnoxa,lnoxcl, &
          lnoxsq,lo3,n,n2o,n2o5,n2o5e,no2,no3,no,nox,o1d,o3,o3p,oclo,oh,ox,pbrcl,pbrono2, &
         pbrx,pcl2,pclono2,pclx,ph2o2,phbr,phcl,pho2no2,phocl,phno3,phno3het,phobr,pn2o5, &
         pno2,pno3,pnox,pnoxa,pnoxcl,po3,poclo,climRainProfile,rbrnbrx,rbrobrx,rcl2o2clo, &
       rclclo,rcloclx,rnno,rnono2,rno2nox,rno3no2,ro1do3,ro3ox,ro3po3,rocloclo,box_ro3ox, &
                     gcSC%levels_cal,gcSC%numreacs,gcSC%numphoto,gcrProfile,gcSC%o2,szan)


  IF(examineDt <= 3) CALL verifyChemDt(examineDt,requiredPasses,tdt,dt,requiredDt,startOver)

  RETURN
 END SUBROUTINE chemStage2

 END SUBROUTINE SC_GridCompRun

!-------------------------------------------------------------------------
!NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1, GEOS/DAS!
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  SC_GridCompFinalize --- The Chem Driver 
!
! !INTERFACE:
!

   SUBROUTINE SC_GridCompFinalize ( gcSC, bsc, impChem, expChem, &
                                    nymd, nhms, cdt, rc )

! !USES:

  IMPLICIT NONE

! !INPUT/OUTPUT PARAMETERS:

   TYPE(SC_GridComp), INTENT(inout) :: gcSC   ! Grid Component

! !INPUT PARAMETERS:

   TYPE(Species_Bundle),   POINTER :: bsc     ! Bundle of StratChem Species
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
!  11May2012 Nielsen   Capability for FV cubed [Ganymed-1_0_UNSTABLE]
!
!EOP
!-------------------------------------------------------------------------

   CHARACTER(LEN=*), PARAMETER :: Iam = 'SC_GridCompFinalize'
   INTEGER :: status
   
   rc = 0

   DEALLOCATE(gcSC%climLatRad, STAT=status)
   VERIFY_(status)

   DEALLOCATE(gcSC%sdat, STAT=status)
   VERIFY_(status)
   DEALLOCATE(gcSC%o2jdat, STAT=status)
   VERIFY_(status)
   DEALLOCATE(gcSC%o3_tab, STAT=status)
   VERIFY_(status)
   DEALLOCATE(gcSC%xtab, STAT=status)
   VERIFY_(status)
   DEALLOCATE(gcSC%sza_tab, STAT=status)
   VERIFY_(status)
   DEALLOCATE(gcSC%CH2O_aq, STAT=status)
   VERIFY_(status)
   DEALLOCATE(gcSC%rlam, STAT=status)
   VERIFY_(status)

   DEALLOCATE(gcSC%SSG, STAT=status)
   VERIFY_(status)
   DEALLOCATE(gcSC%climRain, STAT=status)
   VERIFY_(status)
   DEALLOCATE(gcSC%climO3, STAT=status)
   VERIFY_(status)
   DEALLOCATE(gcSC%climOH, STAT=status)
   VERIFY_(status)
   DEALLOCATE(gcSC%climH2O, STAT=status)
   VERIFY_(status)

   DEALLOCATE(gcSC%climRain3D, STAT=status)
   VERIFY_(status)
   DEALLOCATE(gcSC%climO33D, STAT=status)
   VERIFY_(status)
   DEALLOCATE(gcSC%climOH3D, STAT=status)
   VERIFY_(status)
   DEALLOCATE(gcSC%climH2O3D, STAT=status)
   VERIFY_(status)

   DEALLOCATE(gcSC%SO4sa, STAT=status)
   VERIFY_(status)
   DEALLOCATE(gcSC%SO4sa3D, STAT=status)
   VERIFY_(status)

   DEALLOCATE(gcSC%CO2MBC, STAT=status)
   VERIFY_(status)
   DEALLOCATE(gcSC%CO2MBC2D, STAT=status)
   VERIFY_(status)
   DEALLOCATE(gcSC%SF6MBC, STAT=status)
   VERIFY_(status)
   DEALLOCATE(gcSC%SF6MBC2D, STAT=status)
   VERIFY_(status)

   DEALLOCATE(gcSC%KRxnOrder, STAT=status)
   VERIFY_(status)
   DEALLOCATE(gcSC%standardKRxn, STAT=status)
   VERIFY_(status)
   DEALLOCATE(gcSC%KRxnName, STAT=status)
   VERIFY_(status)
   DEALLOCATE(gcSC%cnsttab, STAT=status)
   VERIFY_(status)
   DEALLOCATE(gcSC%indxs1, STAT=status)
   VERIFY_(status)
   DEALLOCATE(gcSC%indxs2, STAT=status)
   VERIFY_(status)
   DEALLOCATE(gcSC%indxs3, STAT=status)
   VERIFY_(status)

   DEALLOCATE(gcSC%CFC11EMISS, STAT=status)
   VERIFY_(status)
   DEALLOCATE(gcSC%CFC12EMISS, STAT=status)
   VERIFY_(status)
   DEALLOCATE(gcSC%CFC113EMISS, STAT=status)
   VERIFY_(status)
   DEALLOCATE(gcSC%MCFEMISS, STAT=status)
   VERIFY_(status)
   DEALLOCATE(gcSC%CCL4EMISS, STAT=status)
   VERIFY_(status)
   DEALLOCATE(gcSC%HCFC22EMISS, STAT=status)
   VERIFY_(status)
   DEALLOCATE(gcSC%CHBR3EMISS, STAT=status)
   VERIFY_(status)
   DEALLOCATE(gcSC%CH2BR2EMISS, STAT=status)
   VERIFY_(status)

   RETURN

 END SUBROUTINE SC_GridCompFinalize

 END MODULE SC_GridCompMod

