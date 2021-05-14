#include "MAPL_Generic.h"

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 910.1     !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: GMIchem_GridCompMod - The GMI COMBO Model Grid Component
!
! !INTERFACE:
!
   MODULE GMIchem_GridCompMod
!
! !USES:
!
   USE ESMF
   USE MAPL
   USE Chem_Mod 	                        ! Chemistry Base Class
   USE GMI_GridCompMod                          ! ESMF parent component
   USE Chem_UtilMod, ONLY : Chem_UtilNegFiller  ! Eliminates negative vmr
   USE Chem_GroupMod                            ! For Family Transport
   USE OVP,     ONLY:  OVP_init, OVP_end_of_timestep_hms, OVP_mask, OVP_apply_mask

   IMPLICIT NONE
   PRIVATE

   TYPE(Chem_Mie), DIMENSION(2), SAVE :: gocartMieTable

   INTEGER, SAVE, ALLOCATABLE :: MASK_10AM(:,:)
   INTEGER, SAVE, ALLOCATABLE :: MASK_2PM(:,:)
   INTEGER, SAVE              :: OVP_FIRST_HMS
   INTEGER, SAVE              :: OVP_RUN_DT
   INTEGER, SAVE              :: OVP_GC_DT
   INTEGER, SAVE              :: OVP_MASK_DT

!
! !PUBLIC MEMBER FUNCTIONS:
!
   PUBLIC SetServices

!
! !PRIVATE MEMBER FUNCTIONS:
!
   PRIVATE                    :: Initialize_    ! Init method
   PRIVATE                    :: Run_           ! Run method
   PRIVATE                    :: Finalize_      ! Finalize method

   PRIVATE                    :: Run1           ! Run wrapper phase 1
   PRIVATE                    :: Run2           ! Run wrapper phase 2

   PRIVATE                    :: extract_           ! Get values from ESMF
   PRIVATE                    :: aerosol_optics     ! Get params for aero optics
   PRIVATE                    :: secure_species_ptr ! Find species in w_c%qa


!
! !DESCRIPTION: 
!
!  {\tt GMIchem\_GridComp} is a ESMF gridded component for the Global Modeling
!  Initiative combined troposphere/stratospheric chemistry package.

! !REVISION HISTORY:

!  31Jul2006  da Silva  Created the GMI stub.
!  11Dec2007  Nielsen   Real code for Eros-beta7p17.
!  25Nov2011  Nielsen   Trying cubed sphere.
!  10Sep2013  Nielsen   Added run alarm, but allow for updating age-of-air and for 
!                       returning zero tendencies, etc., when alarm is not ringing.
!
!EOP
!-------------------------------------------------------------------------

  TYPE GMIchem_State
     PRIVATE
     TYPE(Chem_Registry), POINTER :: chemReg  => null()
     TYPE(GMI_GridComp),  POINTER :: gcGMI    => null()
     TYPE(Chem_Bundle),   POINTER :: w_c      => null()
  END TYPE GMIchem_State

  TYPE GMIchem_WRAP
     TYPE (GMIchem_State), pointer :: PTR => null()
  END TYPE GMIchem_WRAP

  ! Number of run phases. Can be set in the resource file (default is 2).
  INTEGER                          :: PHASE_COUNT

CONTAINS


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 910.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: SetServices --- Sets IRF services for GMIchem Grid Component
!
! !INTERFACE:

   SUBROUTINE SetServices ( GC, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
    integer, optional                  :: RC  ! return code

! !DESCRIPTION: Sets Initialize, Run and Finalize services. 
!
! !REVISION HISTORY:
!
!  31Jul2006  da Silva  First crack.
!
!EOP
!-------------------------------------------------------------------------

!   ErrLog Variables
!   ----------------
    character(len=ESMF_MAXSTR)      :: IAm
    integer                         :: STATUS
    character(len=ESMF_MAXSTR)      :: COMP_NAME

!   Local derived type aliases
!   --------------------------
    type (ESMF_Config)              :: CF
    type (GMIchem_State), pointer   :: state   ! internal, that is
    type (GMIchem_wrap)             :: wrap

    integer                         :: m, n, i_XX, j_XX
    CHARACTER(LEN=ESMF_MAXSTR)      :: FRIENDLIES
    CHARACTER(LEN=ESMF_MAXSTR)      :: providerName
    CHARACTER(LEN=ESMF_MAXSTR)      :: aeroProviderName

    LOGICAL :: searchForImports
    INTEGER, PARAMETER :: numAeros = 5
    CHARACTER(LEN=2) :: aeroID(numAeros) = (/ "BC", "DU", "OC", "SS", "SU" /)
    CHARACTER(LEN=2) :: leadChars
    CHARACTER(LEN=ESMF_MAXSTR) :: name

    LOGICAL :: do_ShipEmission
    INTEGER :: fastj_opt
    TYPE(ESMF_Config)  :: gmiConfig

    ! HEMCO isoprene related -sas
    CHARACTER(LEN=255) :: gmi_rcfilen = 'GMI_GridComp.rc'
    TYPE (ESMF_Config) :: gmi_config
    LOGICAL :: doMEGANviaHEMCO

!   Get my name and set-up traceback handle
!   ---------------------------------------
    call ESMF_GridCompGet( GC, NAME=COMP_NAME, CONFIG=CF, RC=STATUS )
    VERIFY_(STATUS)
    Iam = TRIM(COMP_NAME)//"::SetServices"

!   Wrap internal state for storing in GC; rename legacyState
!   -------------------------------------
    allocate ( state, stat=STATUS )
    VERIFY_(STATUS)
    wrap%ptr => state

!   Start by loading the Chem Registry
!   ----------------------------------
    allocate ( state%chemReg )
    state%chemReg = Chem_RegistryCreate ( STATUS )
    VERIFY_(STATUS)

!                       ------------------------
!                       ESMF Functional Services
!                       ------------------------


    IF(MAPL_AM_I_ROOT()) THEN
     PRINT *, TRIM(Iam)//': ACTIVE'
     CALL Chem_RegistryPrint ( state%chemReg )
    END IF


!   Set the Initialize, Run, Finalize entry points
!   ----------------------------------------------
    CALL MAPL_GridCompSetEntryPoint(GC, ESMF_METHOD_INITIALIZE, Initialize_, __RC__ )
    CALL MAPL_GridCompSetEntryPoint(GC, ESMF_METHOD_RUN,        Run1,        __RC__ )
    CALL MAPL_GridCompSetEntryPoint(GC, ESMF_METHOD_RUN,        Run2,        __RC__ )
    CALL MAPL_GridCompSetEntryPoint(GC, ESMF_METHOD_FINALIZE,   Finalize_,   __RC__ )

!   Store internal state in GC
!   --------------------------
    CALL ESMF_UserCompSetInternalState(GC, 'GMIchem_state', wrap, STATUS)
    VERIFY_(STATUS)

! ========================= IMPORT STATE =========================

    CALL ESMF_ConfigGetAttribute(CF, aeroProviderName, DEFAULT='none', LABEL="AERO_PROVIDER:", __RC__ )

    SELECT CASE (TRIM(aeroProviderName))

     CASE ("GOCART")

!   GOCART aerosols and dust.  With the exception of
!   SO4, select upon dust/aerosol two-letter identifier 
!   ---------------------------------------------------
      DO m=1,numAeros-1

       searchForImports = .FALSE.

       SELECT CASE (aeroID(m))
        CASE("BC")
         searchForImports = state%chemReg%doing_BC
        CASE("DU")
         searchForImports = state%chemReg%doing_DU
        CASE("OC")
         searchForImports = state%chemReg%doing_OC
        CASE("SS")
         searchForImports = state%chemReg%doing_SS
        CASE DEFAULT
         searchForImports = .FALSE.
       END SELECT

       Doing_Search: IF(searchForImports) THEN

        IF(MAPL_AM_I_ROOT() .AND. m == 1) PRINT *,"Adding the following from GOCART to GMICHEM import state:"

        DO n = state%chemReg%i_GOCART, state%chemReg%j_GOCART

         name = TRIM(state%chemReg%vname(n))
         leadChars = ESMF_UtilStringUpperCase(name(1:2))

         Match: IF(leadChars == aeroID(m)) THEN

          CALL MAPL_AddImportSpec(GC,                                  &
               SHORT_NAME  = "GOCART::"//TRIM(state%chemReg%vname(n)), &
               LONG_NAME   = state%chemReg%vtitle(n),		     &
               UNITS       = state%chemReg%vunits(n),		     &
               DIMS        = MAPL_DimsHorzVert,		             &
               VLOCATION   = MAPL_VLocationCenter,	  RC=STATUS  )
          VERIFY_(STATUS)

          IF(MAPL_AM_I_ROOT()) PRINT *,"  ",TRIM(state%chemReg%vname(n))

         END IF Match

        END DO 

       END IF Doing_Search

      END DO

!   This is the special case for SO4, which 
!   does not have "SU" as its leading two characters
!   ------------------------------------------------
      IF(state%chemReg%doing_SU) THEN

       Doing_SO4: DO n = state%chemReg%i_SU, state%chemReg%j_SU

        IF( (TRIM(state%chemReg%vname(n)) == "SO4" ) .OR.  &
            (TRIM(state%chemReg%vname(n)) == "SO4v")       )  THEN
         CALL MAPL_AddImportSpec(GC,                                  &
     	      SHORT_NAME  = "GOCART::"//TRIM(state%chemReg%vname(n)), &
     	      LONG_NAME   = state%chemReg%vtitle(n),                  &
     	      UNITS       = state%chemReg%vunits(n),                  &
     	      DIMS        = MAPL_DimsHorzVert,                        &
     	      VLOCATION   = MAPL_VLocationCenter,           RC=STATUS )
         VERIFY_(STATUS)
         IF(MAPL_AM_I_ROOT()) PRINT *,"  ",TRIM(state%chemReg%vname(n))
        END IF

       END DO Doing_SO4

      END IF

     CASE("GOCART.data")

      CALL MAPL_AddImportSpec(GC,				    &
	  SHORT_NAME	     = 'AERO',  			    &
	  LONG_NAME	     = 'aerosol_mass_mixing_ratios',	    &
	  UNITS 	     = 'kg kg-1',			    &
	  DIMS  	     = MAPL_DimsHorzVert,		    &
	  VLOCATION	     = MAPL_VLocationCenter,		    &
	  DATATYPE	     = MAPL_StateItem,  		    &
          RESTART            = MAPL_RestartSkip,                    &
							RC=STATUS  )
      VERIFY_(STATUS)


     CASE("GMICHEM")

      STATUS = 0

     CASE DEFAULT

      PRINT *, TRIM(Iam)//": Invalid AERO_PROVIDER when running GMIChem."
      STATUS = 1
      VERIFY_(STATUS)
    
    END SELECT

! Import NO from Ships, only if using the parameterization
! Import RI and RL, only if using Cloud-J

    gmiConfig = ESMF_ConfigCreate(__RC__)

    call ESMF_ConfigLoadFile(gmiConfig, 'GMI_GridComp.rc', __RC__)

!   This duplicates the call in the Emissions code; really should only be done once!
!   call rcEsmfReadLogical(gmiConfig, do_ShipEmission, "do_ShipEmission:", default=.false., __RC__)
    CALL ESMF_ConfigGetAttribute(gmiConfig, value= do_ShipEmission, Default=.false., &
                                            Label="do_ShipEmission:", __RC__)

    IF ( do_ShipEmission ) THEN
       call MAPL_AddImportSpec(GC,                         & 
          SHORT_NAME = 'SHIP_NO',                          &
          LONG_NAME  = 'NO from Ships',                    &
          UNITS      = 'kg NO m^(-2) s^(-1)',              &
          DIMS       = MAPL_DimsHorzOnly,                  &
          VLOCATION  = MAPL_VLocationNone,   __RC__) 
    END IF

    CALL ESMF_ConfigGetAttribute(gmiConfig, value= fastj_opt, Default=4, &
                                            Label="fastj_opt:", __RC__)

    ! We need RI and RL for Cloud-J
    ! The fields may not be available in CTM, so we import them conditionally
    IF ( fastj_opt == 5 ) THEN
       call MAPL_AddImportSpec(GC,                                           &
          SHORT_NAME         = 'RI',                                         &
          LONG_NAME          = 'ice_phase_cloud_particle_effective_radius',  &
          UNITS              = 'm',                                          &
          DIMS               = MAPL_DimsHorzVert,                            &
          VLOCATION          = MAPL_VLocationCenter,    __RC__)

       call MAPL_AddImportSpec(GC,                                           &
          SHORT_NAME         = 'RL',                                         &
          LONG_NAME          = 'liquid_cloud_particle_effective_radius',     &
          UNITS              = 'm',                                          &
          DIMS               = MAPL_DimsHorzVert,                            &
          VLOCATION          = MAPL_VLocationCenter,    __RC__)
    END IF

    call ESMF_ConfigDestroy(gmiConfig, __RC__)


! Future option - import OCS from ACHEM -  if (state%chemReg%doing_OCS) then import ACHEM::OCS

     call MAPL_AddImportSpec(GC,                           & 
        SHORT_NAME = 'OCS_CLIMO',                          &
        LONG_NAME  = 'Carbonyl Sulfide (OCS gas)',         &
        UNITS      = 'mol mol-1',                          &
        DIMS       = MAPL_DimsHorzVert,                    &
        VLOCATION  = MAPL_VLocationCenter,   __RC__) 

     call MAPL_AddImportSpec(GC,                           & 
        SHORT_NAME = 'CNV_FRC',                            &
        LONG_NAME  = 'convective_fraction',                &
        UNITS      = '',                                   &
        DIMS       = MAPL_DimsHorzOnly,                    &
        VLOCATION  = MAPL_VLocationNone,     __RC__) 

     ! add MEGAN emission imports -sas 
     ! note: might change this later to use GMICHEM_ImportSpec___.h
     gmi_config = ESMF_ConfigCreate(__RC__ )
     call ESMF_ConfigLoadFile(gmi_config, TRIM(gmi_rcfilen), __RC__ )
     call ESMF_ConfigGetAttribute(gmi_config, value= doMEGANviaHEMCO, &
                                              label="doMEGANviaHEMCO:", default=.false., __RC__ )

     IF ( doMEGANviaHEMCO ) THEN
        call MAPL_AddImportSpec(GC, &
             SHORT_NAME = 'GMI_ISOPRENE', &
             LONG_NAME  = 'isoprene emissions'  , &
             UNITS      = 'kgC/m2/s',                &
             DIMS       = MAPL_DimsHorzOnly,  &
             VLOCATION  = MAPL_VLocationNone, &
             RESTART    = MAPL_RestartSkip,   &
             RC         = STATUS)
        VERIFY_(STATUS)
     END IF

#include "GMICHEM_ImportSpec___.h"

! ======================== INTERNAL STATE =========================

! Is GMICHEM providing ozone (mole fraction) to the ANALYSIS bundle?
! ------------------------------------------------------------------
     CALL ESMF_ConfigGetAttribute(CF, providerName, Default="PCHEM", &
                                  Label="ANALYSIS_OX_PROVIDER:", RC=STATUS )
     VERIFY_(STATUS)


!   Species to be transported:
!   --------------------------
    DO n = state%chemReg%i_GMI, state%chemReg%j_GMI

	  IF(TRIM(state%chemReg%vname(n)) == "OX" .AND. TRIM(providerName) == "GMICHEM") THEN
           FRIENDLIES="ANALYSIS:DYNAMICS:TURBULENCE:MOIST"
	  ELSE
           FRIENDLIES="DYNAMICS:TURBULENCE:MOIST"
	  END IF

          FRIENDLIES = FRIENDLIES//':'//TRIM(COMP_NAME)
	 
          CALL MAPL_AddInternalSpec(GC,                                  &
               SHORT_NAME         = TRIM(state%chemReg%vname(n)),        &
               LONG_NAME          = TRIM(state%chemReg%vtitle(n)),       &
               UNITS              = TRIM(state%chemReg%vunits(n)),       &     
               FRIENDLYTO         = TRIM(FRIENDLIES),                    &
               DIMS               = MAPL_DimsHorzVert,                   &
               VLOCATION          = MAPL_VLocationCenter,     RC=STATUS  )
          VERIFY_(STATUS)

    END DO

!   Non-transported species
!   -----------------------
    DO n = state%chemReg%i_XX, state%chemReg%j_XX 

          CALL MAPL_AddInternalSpec(GC,                               &
               SHORT_NAME      = TRIM(state%chemReg%vname(n)),        &
               LONG_NAME       = TRIM(state%chemReg%vtitle(n)),       &
               UNITS           = TRIM(state%chemReg%vunits(n)),       &     
               FRIENDLYTO      = TRIM(COMP_NAME),                     &
               ADD2EXPORT      = .TRUE.,                              &
               DIMS            = MAPL_DimsHorzVert,                   &
               VLOCATION       = MAPL_VLocationCenter,     RC=STATUS  )
          VERIFY_(STATUS)

    END DO

!   Families to be transported:
!   --------------------------
          FRIENDLIES="DYNAMICS"
          CALL MAPL_AddInternalSpec(GC,                                  &
               SHORT_NAME         = "Bry",                               &
               LONG_NAME          = "Bry for family transport",          &
               UNITS              = "mol mol-1",                         &
               FRIENDLYTO         = TRIM(FRIENDLIES),                    &
               RESTART            = MAPL_RestartSkip,                    &
               DIMS               = MAPL_DimsHorzVert,                   &
               VLOCATION          = MAPL_VLocationCenter,     RC=STATUS  )
          VERIFY_(STATUS)
          CALL MAPL_AddInternalSpec(GC,                                  &
               SHORT_NAME         = "Cly",                               &
               LONG_NAME          = "Cly for family transport",          &
               UNITS              = "mol mol-1",                         &
               FRIENDLYTO         = TRIM(FRIENDLIES),                    &
               RESTART            = MAPL_RestartSkip,                    &
               DIMS               = MAPL_DimsHorzVert,                   &
               VLOCATION          = MAPL_VLocationCenter,     RC=STATUS  )
          VERIFY_(STATUS)


! ========================== EXPORT STATE =========================

    IF(TRIM(aeroProviderName) == "GMICHEM") THEN

!   This state is needed by radiation, and contains aerosols and aerosol optics
!   ---------------------------------------------------------------------------
     CALL MAPL_AddExportSpec(GC,                                   &
         SHORT_NAME         = 'AERO',                              &
         LONG_NAME          = 'aerosol_mass_mixing_ratios',        &
         UNITS              = 'kg kg-1',                           &
         DIMS               = MAPL_DimsHorzVert,                   &
         VLOCATION          = MAPL_VLocationCenter,                &
         DATATYPE           = MAPL_StateItem,                      &
                                                        RC=STATUS  )
     VERIFY_(STATUS)

!   This state is needed by MOIST - It should contain aerosol info.
!   For GMI aerosols, the functionality is turned off.
!   ---------------------------------------------------------------
     CALL MAPL_AddExportSpec(GC,                     &
         SHORT_NAME = 'AERO_ACI',                    &
         LONG_NAME  = 'aerosol_cloud_interaction',   &
         UNITS      = 'kg kg-1',                     &
         DIMS       = MAPL_DimsHorzVert,             &
         VLOCATION  = MAPL_VLocationCenter,          &
         DATATYPE   = MAPL_StateItem, __RC__)

! This bundle is reserved for SURFACE to update snow albedo due to
! aerosol settling and deposition. But since GMI's aerosols are prescribed,
! they are considered "data-driven," and the bundle is not filled by GMICHEM.
! ---------------------------------------------------------------------------
     CALL MAPL_AddExportSpec(GC,                                  &
         SHORT_NAME         = 'AERO_DP',                          &
         LONG_NAME          = 'aerosol_deposition',               &
         UNITS              = 'kg m-2 s-1',                       &
         DIMS               = MAPL_DimsHorzOnly,                  &
         VLOCATION          = MAPL_VLocationNone,                 &
         DATATYPE           = MAPL_BundleItem,                    &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

    END IF

! SAD related variables
! ---------------------

    CALL MAPL_AddExportSpec(GC,                                  &
       SHORT_NAME         = 'HNO3CONDsad'   ,                    &
       LONG_NAME          = 'condensed_phase_hno3',              &
       UNITS              = 'mixing_ratio',                      &
       DIMS               = MAPL_DimsHorzVert,                   &
       VLOCATION          = MAPL_VLocationCenter,                &
                                                      RC=STATUS  )
    VERIFY_(STATUS)

    CALL MAPL_AddExportSpec(GC,                                  &
       SHORT_NAME         = 'HNO3GASsad'   ,                     &
       LONG_NAME          = 'gas_phase-hno3',                    &
       UNITS              = 'mixing_ratio',                      &
       DIMS               = MAPL_DimsHorzVert,                   &
       VLOCATION          = MAPL_VLocationCenter,                &
                                                      RC=STATUS  )
    VERIFY_(STATUS)

! Ship Emissions
! --------------

    CALL MAPL_AddExportSpec(GC,                                  &
       SHORT_NAME         = 'jNO2val',                           &
       LONG_NAME          = 'photolysis_rate_constants_for_NO',  &
       UNITS              = 's^-1',                              &
       DIMS               = MAPL_DimsHorzOnly,                   &
       VLOCATION          = MAPL_VLocationNone,                  &
                                                      RC=STATUS  )
    VERIFY_(STATUS)

    CALL MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME         = 'surfEmissForChem',                 &
        LONG_NAME          = 'surface_emission_for_chemistry',   &
        UNITS              = 'kg m-2 s-1',                       &
        DIMS               = MAPL_DimsHorzOnly,                  &
        VLOCATION          = MAPL_VLocationNone,                 &
        DATATYPE           = MAPL_BundleItem,                    &
                                                      RC=STATUS  )
    VERIFY_(STATUS)

! Aerosol Surface Area Densities (SAD) bundle.
! --------------------------------------------

    CALL MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME         = 'gmiSAD',                           &
        LONG_NAME          = 'surface_area_densities',           &
        UNITS              = 'cm^2/cm^3',                        &
        DIMS               = MAPL_DimsHorzVert,                  &
        VLOCATION          = MAPL_VLocationCenter,               &
        DATATYPE           = MAPL_BundleItem,                    &
                                                      RC=STATUS  )
    VERIFY_(STATUS)

! Photolysis Rate Constants bundle.
! --------------------------------

    CALL MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME         = 'gmiQJ',                            &
        LONG_NAME          = 'photolysis_rate_constants',        &
        UNITS              = 'cm3 s-1',                          &
        DIMS               = MAPL_DimsHorzVert,                  &
        VLOCATION          = MAPL_VLocationCenter,               &
        DATATYPE           = MAPL_BundleItem,                    &
                                                      RC=STATUS  )
    VERIFY_(STATUS)

    CALL MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME         = 'gmiQQJ',                           &
        LONG_NAME          = 'photolysis_reaction_rates',        &
        UNITS              = 'cm-3 s-1',                         &
        DIMS               = MAPL_DimsHorzVert,                  &
        VLOCATION          = MAPL_VLocationCenter,               &
        DATATYPE           = MAPL_BundleItem,                    &
                                                      RC=STATUS  )
    VERIFY_(STATUS)

! Thermal Rate Constants bundle.
! --------------------------------

    CALL MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME         = 'gmiQK',                            &
        LONG_NAME          = 'thermal_rate_constants',           &
        UNITS              = '2-3body_varies',                   &
        DIMS               = MAPL_DimsHorzVert,                  &
        VLOCATION          = MAPL_VLocationCenter,               &
        DATATYPE           = MAPL_BundleItem,                    &
                                                      RC=STATUS  )
    VERIFY_(STATUS)

    CALL MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME         = 'gmiQQK',                           &
        LONG_NAME          = 'thermal_reaction_rates',           &
        UNITS              = 'cm-3 s-1',                         &
        DIMS               = MAPL_DimsHorzVert,                  &
        VLOCATION          = MAPL_VLocationCenter,               &
        DATATYPE           = MAPL_BundleItem,                    &
                                                      RC=STATUS  )
    VERIFY_(STATUS)

! Aerosol or Dust Radii bundle.
! -----------------------------

    CALL MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME         = 'gmiERADIUS',                       &
        LONG_NAME          = 'Aerosol_Dust_Radii',               &
        UNITS              = 'cm',                               &
        DIMS               = MAPL_DimsHorzVert,                  &
        VLOCATION          = MAPL_VLocationCenter,               &
        DATATYPE           = MAPL_BundleItem,                    &
                                                      RC=STATUS  )

! Surface Area of Aerosol or Dust bundle.
! ---------------------------------------

    CALL MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME         = 'gmiTAREA',                         &
        LONG_NAME          = 'surface_area_aerosol_dust',        &
        UNITS              = 'cm^2/cm^3',                        &
        DIMS               = MAPL_DimsHorzVert,                  &
        VLOCATION          = MAPL_VLocationCenter,               &
        DATATYPE           = MAPL_BundleItem,                    &
                                                      RC=STATUS  )
    VERIFY_(STATUS)

    CALL MAPL_AddExportSpec(GC,                                   &
        SHORT_NAME         = 'OVP10_NO2',                         &
        LONG_NAME          = 'Nitrogen_dioxide_10am_local',       &
        UNITS              = 'mol mol-1',                         &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                           __RC__ )

    CALL MAPL_AddExportSpec(GC,                                   &
        SHORT_NAME         = 'OVP14_NO2',                         &
        LONG_NAME          = 'Nitrogen_dioxide_2pm_local',        &
        UNITS              = 'mol mol-1',                         &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                           __RC__ )

    CALL MAPL_AddExportSpec(GC,                                   &
        SHORT_NAME         = 'OVP10_CH2O',                        &
        LONG_NAME          = 'Formaldehyde_10am_local',           &
        UNITS              = 'mol mol-1',                         &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                           __RC__ )

    CALL MAPL_AddExportSpec(GC,                                   &
        SHORT_NAME         = 'OVP14_CH2O',                        &
        LONG_NAME          = 'Formaldehyde_2pm_local',            &
        UNITS              = 'mol mol-1',                         &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                           __RC__ )

    CALL MAPL_AddExportSpec(GC,                                   &
        SHORT_NAME         = 'OVP10_CO',                          &
        LONG_NAME          = 'Carbon_monoxide_10am_local',        &
        UNITS              = 'mol mol-1',                         &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                           __RC__ )

    CALL MAPL_AddExportSpec(GC,                                   &
        SHORT_NAME         = 'OVP14_CO',                          &
        LONG_NAME          = 'Carbon_monoxide_2pm_local',         &
        UNITS              = 'mol mol-1',                         &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                           __RC__ )

    CALL MAPL_AddExportSpec(GC,                                   &
        SHORT_NAME         = 'OVP10_OX',                          &
        LONG_NAME          = 'Ozone_10am_local',                  &
        UNITS              = 'mol mol-1',                         &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                           __RC__ )

    CALL MAPL_AddExportSpec(GC,                                   &
        SHORT_NAME         = 'OVP14_OX',                          &
        LONG_NAME          = 'Ozone_2pm_local',                   &
        UNITS              = 'mol mol-1',                         &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                           __RC__ )

    CALL MAPL_AddExportSpec(GC,                                   &
        SHORT_NAME         = 'OVP10_CH4',                         &
        LONG_NAME          = 'Methane_10am_local',                &
        UNITS              = 'mol mol-1',                         &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                           __RC__ )

    CALL MAPL_AddExportSpec(GC,                                   &
        SHORT_NAME         = 'OVP14_CH4',                         &
        LONG_NAME          = 'Methane_2pm_local',                 &
        UNITS              = 'mol mol-1',                         &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                           __RC__ )

    CALL MAPL_AddExportSpec(GC,                                   &
        SHORT_NAME         = 'OVP14_BrO',                         &
        LONG_NAME          = 'Bromine_monoxide_2pm_local',        &
        UNITS              = 'mol mol-1',                         &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                           __RC__ )

    CALL MAPL_AddExportSpec(GC,                                   &
        SHORT_NAME         = 'OVP14_PAN',                         &
        LONG_NAME          = 'Peroxyacetyl_nitrate_2pm_local',    &
        UNITS              = 'mol mol-1',                         &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                           __RC__ )

    CALL MAPL_AddExportSpec(GC,                                   &
        SHORT_NAME         = 'OVP14_ISOP',                        &
        LONG_NAME          = 'Isoprene_2pm_local',        &
        UNITS              = 'mol mol-1',                         &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                           __RC__ )

    CALL MAPL_AddExportSpec(GC,                                   &
        SHORT_NAME         = 'OVP14_MOH',                         &
        LONG_NAME          = 'Methanol_2pm_local',                &
        UNITS              = 'mol mol-1',                         &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                           __RC__ )

    CALL MAPL_AddExportSpec(GC,                                   &
        SHORT_NAME         = 'OVP14_HCOOH',                       &
        LONG_NAME          = 'Formic_Acid_2pm_local',             &
        UNITS              = 'mol mol-1',                         &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                           __RC__ )

    CALL MAPL_AddExportSpec(GC,                                   &
        SHORT_NAME         = 'OVP14_N2O',                         &
        LONG_NAME          = 'Nitrous_oxide_2pm_local',           &
        UNITS              = 'mol mol-1',                         &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                           __RC__ )

    CALL MAPL_AddExportSpec(GC,                                   &
        SHORT_NAME         = 'OVP14_HCl',                         &
        LONG_NAME          = 'Hydrochloric_acid_2pm_local',       &
        UNITS              = 'mol mol-1',                         &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                           __RC__ )

    CALL MAPL_AddExportSpec(GC,                                   &
        SHORT_NAME         = 'OVP14_HNO3',                        &
        LONG_NAME          = 'Nitric_acid_2pm_local',             &
        UNITS              = 'mol mol-1',                         &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                           __RC__ )

    CALL MAPL_AddExportSpec(GC,                                   &
        SHORT_NAME         = 'OVP14_ClO',                         &
        LONG_NAME          = 'Chlorine_monoxide_radical_2pm_local', &
        UNITS              = 'mol mol-1',                         &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                           __RC__ )

    CALL MAPL_AddExportSpec(GC,                                   &
        SHORT_NAME         = 'OVP14_CH3Cl',                       &
        LONG_NAME          = 'Methyl_chloride_2pm_local',         &
        UNITS              = 'mol mol-1',                         &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                           __RC__ )

    CALL MAPL_AddExportSpec(GC,                                   &
        SHORT_NAME         = 'OVP10_NO',                          &
        LONG_NAME          = 'Nitric_oxide_10am_local',           &
        UNITS              = 'mol mol-1',                         &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                           __RC__ )

    CALL MAPL_AddExportSpec(GC,                                   &
        SHORT_NAME         = 'OVP14_NO',                          &
        LONG_NAME          = 'Nitric_oxide_2pm_local',            &
        UNITS              = 'mol mol-1',                         &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                           __RC__ )

    CALL MAPL_AddExportSpec(GC,                                   &
        SHORT_NAME         = 'OVP14_HOCl',                        &
        LONG_NAME          = 'Hypochlorous_acid_2pm_local',       &
        UNITS              = 'mol mol-1',                         &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                           __RC__ )

    CALL MAPL_AddExportSpec(GC,                                   &
        SHORT_NAME         = 'OVP14_ClONO2',                      &
        LONG_NAME          = 'Chlorine_nitrate_2pm_local',        &
        UNITS              = 'mol mol-1',                         &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                           __RC__ )

    CALL MAPL_AddExportSpec(GC,                                   &
        SHORT_NAME         = 'OVP14_Cl2O2',                       &
        LONG_NAME          = 'Chlorine_peroxide_2pm_local',       &
        UNITS              = 'mol mol-1',                         &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                           __RC__ )

    CALL MAPL_AddExportSpec(GC,                                   &
        SHORT_NAME         = 'OVP14_Cl',                          &
        LONG_NAME          = 'Ground_state_atomic_chlorine_2pm_local', &
        UNITS              = 'mol mol-1',                         &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                           __RC__ )

    CALL MAPL_AddExportSpec(GC,                                   &
        SHORT_NAME         = 'OVP14_Cl2',                         &
        LONG_NAME          = 'Molecular_chlorine_2pm_local',      &
        UNITS              = 'mol mol-1',                         &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                           __RC__ )

    CALL MAPL_AddExportSpec(GC,                                   &
        SHORT_NAME         = 'OVP14_OClO',                        &
        LONG_NAME          = 'Symmetrical_chlorine_dioxide_2pm_local', &
        UNITS              = 'mol mol-1',                         &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                           __RC__ )

    CALL MAPL_AddExportSpec(GC,                                   &
        SHORT_NAME         = 'OVP14_BrCl',                        &
        LONG_NAME          = 'Bromine_chloride_2pm_local',        &
        UNITS              = 'mol mol-1',                         &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                           __RC__ )

    CALL MAPL_AddExportSpec(GC,                                   &
        SHORT_NAME         = 'OVP14_Br',                          &
        LONG_NAME          = 'Ground_state_atomic_bromine_2pm_local', &
        UNITS              = 'mol mol-1',                         &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                           __RC__ )

    CALL MAPL_AddExportSpec(GC,                                   &
        SHORT_NAME         = 'OVP14_BrONO2',                      &
        LONG_NAME          = 'Bromine_nitrate_2pm_local',         &
        UNITS              = 'mol mol-1',                         &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                           __RC__ )

    CALL MAPL_AddExportSpec(GC,                                   &
        SHORT_NAME         = 'OVP14_HBr',                         &
        LONG_NAME          = 'Hydrogen_bromide_2pm_local',        &
        UNITS              = 'mol mol-1',                         &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                           __RC__ )

    CALL MAPL_AddExportSpec(GC,                                   &
        SHORT_NAME         = 'OVP14_HOBr',                        &
        LONG_NAME          = 'Hypobromous_acid_2pm_local',        &
        UNITS              = 'mol mol-1',                         &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                           __RC__ )

    CALL MAPL_AddExportSpec(GC,                                   &
        SHORT_NAME         = 'OVP14_HNO3COND',                    &
        LONG_NAME          = 'Condensed_nitric_acid_2pm_local',   &
        UNITS              = 'mol mol-1',                         &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                           __RC__ )

    CALL MAPL_AddExportSpec(GC,                                   &
        SHORT_NAME         = 'OVP14_OH',                          &
        LONG_NAME          = 'Hydroxyl_radical_2pm_local',        &
        UNITS              = 'mol mol-1',                         &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                           __RC__ )

    CALL MAPL_AddExportSpec(GC,                                   &
        SHORT_NAME         = 'OVP14_HO2',                         &
        LONG_NAME          = 'Perhydroxyl_radical_2pm_local',     &
        UNITS              = 'mol mol-1',                         &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                           __RC__ )
    CALL MAPL_AddExportSpec(GC,                                   &
        SHORT_NAME         = 'OVP14_R4N2',                         &
        LONG_NAME          = 'C4-C5_alkylnitrates_(C4H9O3N)_2pm_local',     &
        UNITS              = 'mol mol-1',                         &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                           __RC__ )

! Overpass fields for EM_LGTNO  (set in GmiEmiss_GridCompClassMod.F90)
! ----------------------------

    CALL MAPL_AddExportSpec(GC,                                         &
        SHORT_NAME         = 'OVP10_EM_LGTNO',                          &
        LONG_NAME          = 'NO_emissions_from_lightning_10am_local',  &
        UNITS              = 'mol mol-1 s-1',                           &
        DIMS               = MAPL_DimsHorzVert,                         &
        VLOCATION          = MAPL_VLocationCenter,                      &
                                                       RC=STATUS  )
    VERIFY_(STATUS)

    CALL MAPL_AddExportSpec(GC,                                         &
        SHORT_NAME         = 'OVP14_EM_LGTNO',                          &
        LONG_NAME          = 'NO_emissions_from_lightning_2pm_local',   &
        UNITS              = 'mol mol-1 s-1',                           &
        DIMS               = MAPL_DimsHorzVert,                         &
        VLOCATION          = MAPL_VLocationCenter,                      &
                                                       RC=STATUS  )
    VERIFY_(STATUS)

! Overpass fields for AIRMASS  (set in GmiEmiss_GridCompClassMod.F90)
! ----------------------------

    CALL MAPL_AddExportSpec(GC,                                         &
        SHORT_NAME         = 'OVP10_AIRMASS',                           &
        LONG_NAME          = 'mass_of_air_in_layer_10am_local',         &
        UNITS              = 'kg m-2',                                  &
        DIMS               = MAPL_DimsHorzVert,                         &
        VLOCATION          = MAPL_VLocationCenter,                      &
                                                       RC=STATUS  )
    VERIFY_(STATUS)

    CALL MAPL_AddExportSpec(GC,                                         &
        SHORT_NAME         = 'OVP14_AIRMASS',                           &
        LONG_NAME          = 'mass_of_air_in_layer_2pm_local',          &
        UNITS              = 'kg m-2',                                  &
        DIMS               = MAPL_DimsHorzVert,                         &
        VLOCATION          = MAPL_VLocationCenter,                      &
                                                       RC=STATUS  )
    VERIFY_(STATUS)

#include "GMICHEM_ExportSpec___.h"
#include "Deposition_ExportSpec___.h"
#include "Reactions_ExportSpec___.h"
#include "Tendency_ExportSpec___.h"

! =================================================================

!   Set the Profiling timers
!   ------------------------
    CALL MAPL_TimerAdd(GC, NAME="INITIALIZE", RC=STATUS)
    VERIFY_(STATUS)
    CALL MAPL_TimerAdd(GC, NAME="RUN", RC=STATUS)
    VERIFY_(STATUS)
    CALL MAPL_TimerAdd(GC, NAME="FINALIZE", RC=STATUS)
    VERIFY_(STATUS)

!   Generic Set Services
!   --------------------
    call MAPL_GenericSetServices ( GC, RC=STATUS )
    VERIFY_(STATUS)

    RETURN_(ESMF_SUCCESS)
  
  END SUBROUTINE SetServices


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Initialize_ --- Initialize GMIchem
!
! !INTERFACE:
!

   SUBROUTINE Initialize_ ( gc, impChem, expChem, clock, rc )

! !USES:

   implicit NONE

! !INPUT PARAMETERS:

   type(ESMF_Clock),  intent(inout) :: clock      ! The clock

! !OUTPUT PARAMETERS:

   type(ESMF_GridComp), intent(inout)  :: gc      ! Grid Component
   type(ESMF_State), intent(inout) :: impChem     ! Import State
   type(ESMF_State), intent(inout) :: expChem     ! Export State
   integer, intent(out) ::  rc                    ! Error return code:
                                                  !  0 - all is well
                                                  !  1 - 

! !DESCRIPTION: This is a simple ESMF wrapper.
!
! !REVISION HISTORY:
!
!  27Feb2005 da Silva  First crack.
!
!EOP
!-------------------------------------------------------------------------


!  ErrLog Variables
!  ----------------
   character(len=ESMF_MAXSTR)      :: IAm
   integer                         :: STATUS
   character(len=ESMF_MAXSTR)      :: COMP_NAME

   type(Chem_Registry), pointer    :: chemReg
   type(GMI_GridComp), pointer     :: gcGMI       ! Grid Component
   type(Chem_Bundle), pointer      :: w_c         ! Chemical tracer fields
   integer                         :: nymd, nhms  ! time of day
   real                            :: gmiDt       ! chemistry timestep (secs)
   real                            :: runDt       ! heartbeat (secs)

   type(ESMF_Config)               :: CF
   type(ESMF_Grid)                 :: grid
 
   integer                         :: i1=1, i2, ig=0, im  ! dist grid indices
   integer                         :: j1=1, j2, jg=0, jm  ! dist grid indices
   integer                         :: km                  ! dist grid indices
   integer                         :: dims(3), k, l, n

   type(Chem_Array), pointer       :: q(:)	   ! array of pointers
   type(MAPL_MetaComp), pointer    :: ggState	   ! GEOS Generic State
   type(ESMF_State)                :: internal
   type(MAPL_VarSpec), pointer     :: InternalSpec(:)

   REAL, POINTER, DIMENSION(:,:)   :: LATS
   REAL, POINTER, DIMENSION(:,:)   :: LONS

   CHARACTER(LEN=ESMF_MAXSTR)	   :: short_name
   CHARACTER(LEN=ESMF_MAXSTR)	   :: diurnal_bb
   CHARACTER(LEN=ESMF_MAXSTR)      :: providerName
   CHARACTER(LEN=ESMF_MAXSTR), POINTER, DIMENSION(:) :: fieldNames

   INTEGER, PARAMETER :: numAeroes = 13
   CHARACTER(LEN=ESMF_MAXSTR) :: aeroName(numAeroes) = (/"BCphobic","BCphilic", &
             "du001   ","du002   ","du003   ","du004   ","OCphobic","OCphilic", &
             "ss001   ","ss003   ","ss004   ","ss005   ","SO4     "/)

   rc = 0

!  Get my name and set-up traceback handle
!  ---------------------------------------
   call ESMF_GridCompGet( GC, NAME=COMP_NAME, CONFIG=CF, RC=STATUS )
   VERIFY_(STATUS)
   Iam = TRIM(COMP_NAME)//"::Initialize_"

!  Get my internal MAPL_Generic state
!  -----------------------------------
   call MAPL_GetObjectFromGC ( GC, ggState, RC=STATUS)
   VERIFY_(STATUS)

!  Start timers
!  ------------
   CALL MAPL_TimerOn(ggSTATE, "TOTAL")
   CALL MAPL_TimerOn(ggSTATE, "INITIALIZE")
   
!  Initialize GEOS Generic
!  ------------------------
   call MAPL_GenericInitialize ( gc, impChem, expChem, clock,  RC=STATUS )
   VERIFY_(STATUS)

!  Get parameters from gc and clock
!  --------------------------------
   call extract_ ( gc, clock, chemReg, gcGMI, w_c, nymd, nhms, gmiDt, runDt, STATUS )
   VERIFY_(STATUS)
   IF(MAPL_AM_I_ROOT()) THEN
    PRINT *," "
    PRINT *, TRIM(Iam)//": GMICHEM time step length: ",gmiDt," seconds"
   END IF

!  Aerosol for radiation
!  ---------------------
   CALL ESMF_ConfigGetAttribute(CF, providerName, Default="PCHEM", &
                                Label="AERO_PROVIDER:", __RC__ )
   gcGMI%gcPhot%aeroProviderName = TRIM(providerName)
   IF(TRIM(providerName) == "GMICHEM") THEN
    gcGMI%gcPhot%AM_I_AERO_PROVIDER = .TRUE.
   ELSE
    gcGMI%gcPhot%AM_I_AERO_PROVIDER = .FALSE.
   END IF

!  Create Chem Bundle
!  ------------------
   call ESMF_GridCompGet ( GC, GRID=grid, rc=STATUS)
   VERIFY_(STATUS)

   call MAPL_GridGet ( grid, globalCellCountPerDim=DIMS, RC=STATUS)
   VERIFY_(STATUS)

   im = dims(1)
   jm = dims(2)

   call ESMF_GridGet(GRID, localDE=0,staggerloc=ESMF_STAGGERLOC_CENTER, &
   		     computationalCount=DIMS, RC=STATUS)
   VERIFY_(STATUS)

!  Associate the Internal State fields with our legacy state 
!  ---------------------------------------------------------
   call MAPL_Get ( ggSTATE, INTERNALSPEC=InternalSpec, &
                   INTERNAL_ESMF_STATE=internal, &
                   LONS=LONS, LATS=LATS, RC=STATUS  )
   VERIFY_(STATUS)

!  Local sizes of three dimensions
!  --------------------------------
   i2 = dims(1)
   j2 = dims(2)
   km = dims(3)

!  Broadcast necessary information to individual GCs
!  -------------------------------------------------
   CALL sendToGCs(STATUS)
   VERIFY_(STATUS)

!  Initalize the legacy state but do not allocate memory for arrays
!  ----------------------------------------------------------------
   call Chem_BundleCreate_ ( chemReg, i1, i2, ig, im, j1, j2, jg, jm, km,  &
                             w_c, lon=lons, lat=lats, &
                             skipAlloc=.true., rc=STATUS )
   VERIFY_(STATUS)

   w_c%grid_esmf = grid
   ALLOCATE(w_c%delp(i1:i2,j1:j2,km),w_c%rh(i1:i2,j1:j2,km),__STAT__)

!  Allow user to specify 1 or 2 RUN phases   (set in AGCM.rc)
!  --------------------------------------------------------------------------
   CALL ESMF_ConfigGetAttribute(CF, PHASE_COUNT, LABEL="GMI_RUN_PHASES:", DEFAULT=2, __RC__ )
   ASSERT_(PHASE_COUNT==1.OR.PHASE_COUNT==2)

!  Activate or de-activate diurnal cycle for biomass burning. Default is OFF.
!  --------------------------------------------------------------------------
   CALL ESMF_ConfigGetAttribute(CF, diurnal_bb, LABEL="DIURNAL_BIOMASS_BURNING:", &
                                DEFAULT="NO", RC=STATUS )
   VERIFY_(STATUS)
   IF(diurnal_bb(1:3) == "yes" .OR. diurnal_bb(1:3) == "YES" .OR. diurnal_bb(1:3) == "Yes") THEN	
    short_name = "will be"
    w_c%diurnal_bb = .TRUE.
   ELSE
    short_name = "will not be"
    w_c%diurnal_bb = .FALSE.
   END IF
   IF(MAPL_AM_I_ROOT()) PRINT *, TRIM(Iam)//': Diurnal cycle '//TRIM(short_name)//" applied to biomass burning."

!  Consistency Checks
!  ------------------
   ASSERT_ ( chemReg%i_XX == (chemReg%j_GMI+1) )
   ASSERT_ ( size(InternalSpec) == (chemReg%n_GMI+chemReg%n_XX) + 2)  ! Bry and Cly families are extra internal species

   do L = 1, size(InternalSpec) - 2

      call MAPL_VarSpecGet ( InternalSpec(L),          &
                             SHORT_NAME = short_name,  &
                             RC=STATUS )
      VERIFY_(STATUS)

!     IF(MAPL_AM_I_ROOT()) print*,'GMI species SHORT NAME '//TRIM(short_name)

      N = chemReg%i_GMI + L - 1 ! Assumption: XX species immediately follow GMI species
      CALL MAPL_GetPointer ( internal, NAME=short_name, ptr=w_c%qa(N)%data3d, &
                             rc = STATUS )
      VERIFY_(STATUS)


      IF ( TRIM(short_name) == 'RCOOH' ) THEN
        IF ( MAXVAL(w_c%qa(N)%data3d) > 1.e-9 ) THEN
          PRINT*,'RCOOH values are too high (GT 1e-9), likely from an old RESTART'
          PRINT*,'Remove RCOOH from gmichem_internal_rst.'
          STATUS = 1
          VERIFY_(STATUS)
        ENDIF
      ENDIF

   end do

!  Call initialize
!  ---------------
   call GMI_GridCompInitialize(gcGMI, w_c, impChem, expChem, nymd, nhms, gmiDt, GC, clock, STATUS)
   VERIFY_(STATUS)

!  Set up Overpass Masks
!  --------------------
   CALL OVP_init ( GC, "GMICHEM_DT:", LONS, OVP_RUN_DT, OVP_GC_DT, __RC__ ) !  Get LONS, timesteps

   ! In this case we update the Exports at every timestep:
   OVP_MASK_DT = OVP_RUN_DT

   OVP_FIRST_HMS = OVP_end_of_timestep_hms( CLOCK, OVP_MASK_DT )
   IF(MAPL_AM_I_ROOT()) PRINT*,'GMICHEM FIRST_HMS =',OVP_FIRST_HMS

   CALL OVP_mask ( LONS=LONS, DELTA_TIME=OVP_MASK_DT, OVERPASS_HOUR=10, MASK=MASK_10AM )
   CALL OVP_mask ( LONS=LONS, DELTA_TIME=OVP_MASK_DT, OVERPASS_HOUR=14, MASK=MASK_2PM  )

!  Initialize the AERO state if GMIChem is the AERO_PROVIDER
!  ---------------------------------------------------------
   Building_AERO: IF(gcGMI%gcPhot%AM_I_AERO_PROVIDER) THEN

    CALL Aero_StateInitialize(STATUS)
    VERIFY_(STATUS)

   END IF Building_AERO

!  Init groups
!  -----------
   CALL Init_GMI_Chem_Groups()

!  Stop timers
!  -----------
   CALL MAPL_TimerOff(ggSTATE, "INITIALIZE")
   CALL MAPL_TimerOff(ggSTATE, "TOTAL")

   RETURN_(ESMF_SUCCESS)

  CONTAINS

   SUBROUTINE sendToGCs(rc)
   IMPLICIT NONE
   INTEGER, INTENT(OUT) :: rc
   rc = 0
   
   gcGMI%gcChem%i1 = i1
   gcGMI%gcChem%i2 = i2
   gcGMI%gcChem%im = im
   gcGMI%gcChem%j1 = j1
   gcGMI%gcChem%j2 = j2
   gcGMI%gcChem%jm = jm
   gcGMI%gcChem%km = km
   
   gcGMI%gcDepos%i1 = i1
   gcGMI%gcDepos%i2 = i2
   gcGMI%gcDepos%im = im
   gcGMI%gcDepos%j1 = j1
   gcGMI%gcDepos%j2 = j2
   gcGMI%gcDepos%jm = jm
   gcGMI%gcDepos%km = km
   
   gcGMI%gcEmiss%i1 = i1
   gcGMI%gcEmiss%i2 = i2
   gcGMI%gcEmiss%im = im
   gcGMI%gcEmiss%j1 = j1
   gcGMI%gcEmiss%j2 = j2
   gcGMI%gcEmiss%jm = jm
   gcGMI%gcEmiss%km = km
   gcGMI%gcEmiss%heartBeat = runDt
   
   gcGMI%gcFBC%i1 = i1
   gcGMI%gcFBC%i2 = i2
   gcGMI%gcFBC%im = im
   gcGMI%gcFBC%j1 = j1
   gcGMI%gcFBC%j2 = j2
   gcGMI%gcFBC%jm = jm
   gcGMI%gcFBC%km = km
   
   gcGMI%gcPhot%i1 = i1
   gcGMI%gcPhot%i2 = i2
   gcGMI%gcPhot%im = im
   gcGMI%gcPhot%j1 = j1
   gcGMI%gcPhot%j2 = j2
   gcGMI%gcPhot%jm = jm
   gcGMI%gcPhot%km = km
   
   gcGMI%gcSAD%i1 = i1
   gcGMI%gcSAD%i2 = i2
   gcGMI%gcSAD%im = im
   gcGMI%gcSAD%j1 = j1
   gcGMI%gcSAD%j2 = j2
   gcGMI%gcSAD%jm = jm
   gcGMI%gcSAD%km = km
   
   gcGMI%gcThermalRC%i1 = i1
   gcGMI%gcThermalRC%i2 = i2
   gcGMI%gcThermalRC%im = im
   gcGMI%gcThermalRC%j1 = j1
   gcGMI%gcThermalRC%j2 = j2
   gcGMI%gcThermalRC%jm = jm
   gcGMI%gcThermalRC%km = km

! Latitudes and longitudes (radians) where needed. Note: 
! Deallocations are done by each respective GridCompFinalize.
! -----------------------------------------------------------
   ALLOCATE(gcGMI%gcDepos%lonRad(1:i2,1:j2),STAT=STATUS)
   VERIFY_(STATUS)
   gcGMI%gcDepos%lonRad = LONS

   ALLOCATE(gcGMI%gcDepos%latRad(1:i2,1:j2),STAT=STATUS)
   VERIFY_(STATUS)
   gcGMI%gcDepos%latRad = LATS

    ALLOCATE(gcGMI%gcEmiss%lonRad(1:i2,1:j2),STAT=STATUS)
   VERIFY_(STATUS)
   gcGMI%gcEmiss%lonRad = LONS

   ALLOCATE(gcGMI%gcEmiss%latRad(1:i2,1:j2),STAT=STATUS)
   VERIFY_(STATUS)
   gcGMI%gcEmiss%latRad = LATS

  ALLOCATE(gcGMI%gcFBC%latRad(1:i2,1:j2),STAT=STATUS)
   VERIFY_(STATUS)
   gcGMI%gcFBC%latRad = LATS

   ALLOCATE(gcGMI%gcPhot%lonRad(1:i2,1:j2),STAT=STATUS)
   VERIFY_(STATUS)
   gcGMI%gcPhot%lonRad = LONS

   ALLOCATE(gcGMI%gcPhot%latRad(1:i2,1:j2),STAT=STATUS)
   VERIFY_(STATUS)
   gcGMI%gcPhot%latRad = LATS

   ALLOCATE(gcGMI%gcSAD%lonRad(1:i2,1:j2),STAT=STATUS)
   VERIFY_(STATUS)
   gcGMI%gcSAD%lonRad = LONS

   ALLOCATE(gcGMI%gcSAD%latRad(1:i2,1:j2),STAT=STATUS)
   VERIFY_(STATUS)
   gcGMI%gcSAD%latRad = LATS

   RETURN_(ESMF_SUCCESS)
   END SUBROUTINE sendToGCs

   SUBROUTINE Aero_StateInitialize(rc)
!  ---------------------------------------------------------------------------------
!   Substantially plagarized from GOCART to imitate, as closely as possible, 
!   building the AERO state, which is a container for aerosols and an optics method. 
!  ---------------------------------------------------------------------------------
   IMPLICIT NONE
   INTEGER, INTENT(OUT) :: rc
   CHARACTER(LEN=ESMF_MAXSTR) :: fieldName
   REAL, POINTER, DIMENSION(:,:,:) :: PTR3D

   TYPE(ESMF_State) :: aero
   type(ESMF_State) :: aero_aci
   TYPE(ESMF_Field) :: field
   TYPE(ESMF_Field) :: renamedField
   TYPE(ESMF_FieldBundle) :: aeroBundle
   
   INTEGER, PARAMETER :: numAttributes = 5
   CHARACTER(LEN=ESMF_MAXSTR) :: attName(numAttributes) = &
         (/'air_pressure_for_aerosol_optics            ', &
           'relative_humidity_for_aerosol_optics       ', &
           'extinction_in_air_due_to_ambient_aerosol   ', &
           'single_scattering_albedo_of_ambient_aerosol', &
           'asymmetry_parameter_of_ambient_aerosol     '/)
   CHARACTER(LEN=ESMF_MAXSTR) :: attValu(numAttributes) = (/'PLE','RH ','EXT','SSA','ASY'/)
   INTEGER :: loc(numAttributes) = (/ MAPL_VLocationEdge, MAPL_VLocationCenter, MAPL_VLocationCenter, &
                                                          MAPL_VLocationCenter, MAPL_VLocationCenter /)

!  GMI aerosols are prescribed, so consider our instantiation "data-driven"
!  ------------------------------------------------------------------------
    INTEGER, PARAMETER :: instance = 2

    rc = 0
    NULLIFY(PTR3D)

    CALL ESMF_StateGet(expChem, 'AERO', aero, __RC__)
    CALL ESMF_AttributeSet(aero, NAME='implements_aerosol_optics_method', VALUE=.TRUE., __RC__)
    aeroBundle = ESMF_FieldBundleCreate(NAME='AEROSOLS', __RC__)
    CALL MAPL_StateAdd(aero, aeroBundle, __RC__)

    DO n = 1,numAeroes
     CALL MAPL_GetPointer(expChem, PTR3D, "GMICHEM::"//TRIM(aeroName(n)), ALLOC=.TRUE., RC=STATUS)
     CALL ESMF_StateGet(expChem, "GMICHEM::"//TRIM(aeroName(n)), field, __RC__)
     renamedField = MAPL_FieldCreate(field, NAME=TRIM(aeroName(n)), __RC__)
     CALL MAPL_FieldBundleAdd(aeroBundle, renamedField, __RC__)
     NULLIFY(PTR3D)
    END DO

    gocartMieTable(instance) = Chem_MieCreate(CF, __RC__)

    CALL ESMF_AttributeSet(aero, NAME='mie_table_instance', VALUE=instance, __RC__)
    CALL ESMF_AttributeSet(aero, NAME='cloud_area_fraction_for_aerosol_optics', VALUE='', __RC__)
    CALL ESMF_AttributeSet(aero, NAME='band_for_aerosol_optics', VALUE=0, __RC__)

    DO n = 1,numAttributes
     CALL ESMF_AttributeSet(aero, NAME=TRIM(attName(n)), VALUE=TRIM(attValu(n)), __RC__)
    END DO

    DO n = 1,numAttributes
     CALL ESMF_AttributeGet(aero, NAME=TRIM(attName(n)), VALUE=fieldName, __RC__)
     IF(fieldName /= '') THEN
      field = MAPL_FieldCreateEmpty(TRIM(fieldName), w_c%grid_esmf, __RC__)
      CALL MAPL_FieldAllocCommit(field, DIMS=MAPL_DimsHorzVert, LOCATION=loc(n), TYPEKIND=MAPL_R4, HW=0, __RC__)
      CALL MAPL_StateAdd(aero, field, __RC__)
     END IF
    END DO

    CALL ESMF_MethodAdd(aero, LABEL='aerosol_optics', USERROUTINE=aerosol_optics, __RC__)

    IF(MAPL_AM_I_ROOT()) THEN
     PRINT *," "
     PRINT *, TRIM(Iam)//": AEROSOLS Bundle Members:" 
     CALL ESMF_FieldBundleGet(aeroBundle, FieldCount=n, RC=STATUS)
     VERIFY_(STATUS)
     ALLOCATE(fieldNames(n), STAT=STATUS)
     VERIFY_(STATUS)
     CALL ESMF_FieldBundleGet(aeroBundle, FieldNameList=fieldNames, RC=STATUS)
     VERIFY_(STATUS)
     WRITE(*,FMT="('  Number  Field name')")
     WRITE(*,FMT="('  ------  ------------------')")
     DO k = 1,n
      WRITE(*,FMT="(I7,3X,A)") k,TRIM(fieldNames(k))
     END DO
     PRINT *," "
     DEALLOCATE(fieldNames, STAT=STATUS)
     VERIFY_(STATUS)
    END IF

!   Turn off the aerosol-cloud interaction 
!   ---------------------------------------------------------------------
    call ESMF_StateGet(expChem, 'AERO_ACI', aero_aci, __RC__)

    ! This attribute indicates if the aerosol optics method is implemented or not. 
    ! Radiation will not call the aerosol optics method unless this attribute is 
    ! explicitly set to true.
    call ESMF_AttributeSet(aero_aci, name='implements_aerosol_activation_properties_method', value=.FALSE., __RC__)


   RETURN_(ESMF_SUCCESS)
   END SUBROUTINE Aero_StateInitialize

  END SUBROUTINE Initialize_


!-------------------------------------------------------------------------
!     NASA/GSFC, Atmospheric Chemistry and Dynamics Lab,  Code 614       !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Run1 --- PHASE 1 of GMI Run_
!
! !DESCRIPTION: Runs Deposition, Emission and Boundary Condition
!               components of GMI
!
! !INTERFACE:
!

   SUBROUTINE Run1 ( gc, impChem, expChem, clock, rc )

! !USES:

  implicit NONE

! !INPUT/OUTPUT PARAMETERS:

   type(ESMF_GridComp), intent(inout) :: gc       ! Grid Component
   type(ESMF_State),    intent(inout) :: impChem  ! Import State
   type(ESMF_State),    intent(inout) :: expChem  ! Export State
   type(ESMF_Clock),    intent(inout) :: clock    ! The clock

! !OUTPUT PARAMETERS:

   integer, intent(out) ::  rc                    ! Error return code:
                                                  !  ESMF_SUCCESS - all is well
                                                  !  any other value - failure
! !REVISION HISTORY:
!
!  21Feb2019 - Oman & Manyin - First crack
!
!EOP
!-------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                     :: phase

    CHARACTER(LEN=ESMF_MAXSTR)  :: Iam
    INTEGER                     :: STATUS

    !=======================================================================
    ! Run1 starts here 
    !=======================================================================

    ! Identify this routine to MAPL
    Iam = 'GMI::Run1'

    ! Call run routine phase 1
    ! Skip this step if the user has specified that only 1 phase is to be used.
    ! In this case, we do all chemistry related processes in Run2.
    IF ( PHASE_COUNT == 2 ) THEN
       phase = 1
       CALL Run_ ( gc, impChem, expChem, clock, phase, __RC__ )
    ENDIF

    ! Return w/ success
    RETURN_(ESMF_SUCCESS)

  END SUBROUTINE Run1
!EOC

!-------------------------------------------------------------------------
!     NASA/GSFC, Atmospheric Chemistry and Dynamics Lab,  Code 614       !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Run2 --- PHASE 2 of GMI Run_
!
! !DESCRIPTION: Runs SAD, Photolysis, ThermalRC and Chemistry components of GMI
!
! !INTERFACE:
!

   SUBROUTINE Run2 ( gc, impChem, expChem, clock, rc )

! !USES:

  implicit NONE

! !INPUT/OUTPUT PARAMETERS:

   type(ESMF_GridComp), intent(inout) :: gc       ! Grid Component
   type(ESMF_State),    intent(inout) :: impChem  ! Import State
   type(ESMF_State),    intent(inout) :: expChem  ! Export State
   type(ESMF_Clock),    intent(inout) :: clock    ! The clock

! !OUTPUT PARAMETERS:

   integer, intent(out) ::  rc                    ! Error return code:
                                                  !  ESMF_SUCCESS - all is well
                                                  !  any other value - failure
! !REVISION HISTORY:
!
!  21Feb2019 - Oman & Manyin - First crack
!
!EOP
!-------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                     :: phase

    CHARACTER(LEN=ESMF_MAXSTR)  :: Iam
    INTEGER                     :: STATUS

    !=======================================================================
    ! Run2 starts here 
    !=======================================================================

    ! Identify this routine to MAPL
    Iam = 'GMI::Run2'

    ! Determine phase number: use 2 for multi-phase runs, 99 otherwise.
    ! If set to 99, all processes are done (emissions, chemistry, etc.)

    IF ( PHASE_COUNT == 2 ) THEN
       phase = 2
    ELSE
       phase = 99
    ENDIF

    CALL Run_ ( gc, impChem, expChem, clock, phase, __RC__ )

    ! Return w/ success
    RETURN_(ESMF_SUCCESS)

  END SUBROUTINE Run2
!EOC

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Run_ --- Runs GMIchem
!
! !INTERFACE:
!

   SUBROUTINE Run_ ( gc, impChem, expChem, clock, phase, rc )

! !USES:

  implicit NONE

! !INPUT PARAMETERS:

   integer,             intent(in   )  :: phase   ! Run phase:
                                                  !  1 = do emission
                                                  !  2 = do chemistry
                                                  ! 99 = do both (original approach)
! !INPUT/OUTPUT PARAMETERS:

   type(ESMF_GridComp), intent(inout)  :: gc      ! Grid Component
   type(ESMF_State),    intent(inout)  :: impChem ! Import State
   type(ESMF_State),    intent(inout)  :: expChem ! Export State
   type(ESMF_Clock),    intent(inout)  :: clock   ! The clock

! !OUTPUT PARAMETERS:

   integer, intent(out) ::  rc                    ! Error return code:
                                                  !  ESMF_SUCCESS - all is well
                                                  !  any other value - failure

! !DESCRIPTION: This is a simple ESMF wrapper.
!
! !REVISION HISTORY:
!
!  27Feb2005 da Silva  First crack.
!  10Sep2013 Nielsen   Added run alarm, but allow for updating age-of-air and for 
!                      returning zero tendencies, etc., when alarm is not ringing.
!
!EOP
!-------------------------------------------------------------------------

!  ErrLog Variables
!  ----------------
   character(len=ESMF_MAXSTR)      :: IAm
   integer                         :: STATUS
   character(len=ESMF_MAXSTR)      :: COMP_NAME

   type(Chem_Registry), pointer    :: chemReg
   type(GMI_GridComp), pointer     :: gcGMI       ! Grid Component
   type(Chem_Bundle), pointer      :: w_c         ! Chemical tracer fields     
   integer                         :: nymd, nhms  ! time
   real                            :: gmiDt       ! chemistry timestep (secs)
   real                            :: runDt       ! heartbeat (secs)
   integer                         :: i, i2, iOX, iT2M, iOCS, j2, k, km, m, n
   LOGICAL                         :: RunGMINow

   type(ESMF_Config)               :: CF
   type(ESMF_Grid)                 :: grid        
   type(ESMF_Time)                 :: TIME
   TYPE(ESMF_Alarm)                :: ALARM

   type(MAPL_MetaComp), pointer    :: ggState      ! GEOS Generic State

   real, pointer, dimension(:,:,:) :: rh2
   real, pointer, dimension(:,:)   :: LATS
   real, pointer, dimension(:,:)   :: LONS

   REAL :: dtInverse
   REAL, POINTER, DIMENSION(:,:,:) :: WATER
   REAL, POINTER, DIMENSION(:,:,:) :: SPECHUM
   REAL, POINTER, DIMENSION(:,:,:) :: Q_TEND
   REAL, POINTER, DIMENSION(:,:,:) :: OX_TEND
   REAL, POINTER, DIMENSION(:,:,:) :: PLE
   REAL, POINTER, DIMENSION(:,:,:) :: DELP
   REAL, POINTER, DIMENSION(:,:,:) :: OCS_import
   REAL, POINTER, DIMENSION(:,:)   :: TROPP
   REAL, POINTER, DIMENSION(:,:)   :: AGCMTROPP
   REAL, POINTER, DIMENSION(:,:)   :: GMITROPP
   REAL, POINTER, DIMENSION(:,:)   :: TO3
   REAL, POINTER, DIMENSION(:,:)   :: TTO3
   REAL, ALLOCATABLE               :: wrk(:,:)
   REAL, ALLOCATABLE               :: wgt(:,:)

! Overpass Bundle
! ---------------
   REAL, POINTER, DIMENSION(:,:,:) :: DATA_FOR_OVP_3D => NULL()
   REAL, POINTER, DIMENSION(:,:,:) :: OVP10_OUTPUT_3D => NULL()
   REAL, POINTER, DIMENSION(:,:,:) :: OVP14_OUTPUT_3D => NULL()
   INTEGER                         :: CURRENT_HMS  ! for the end of the timestep

! Tendency Bundle
! ---------------
   TYPE(ESMF_Field)                  :: FIELD
   REAL, POINTER, DIMENSION(:,:,:,:) :: sInitial
   REAL, POINTER, DIMENSION(:,:,:)   :: sIncrement
   LOGICAL                           :: doingTendencies
   INTEGER                           :: nAlloc
   LOGICAL, ALLOCATABLE              :: doMyTendency(:)
   CHARACTER(LEN=ESMF_MAXSTR)        :: fieldName, incFieldName

! Replay mode detection
! ---------------------
   CHARACTER(LEN=ESMF_MAXSTR)        :: ReplayMode
   TYPE(ESMF_Alarm)                  :: PredictorIsActive
   LOGICAL                           :: doingPredictorNow

   LOGICAL                           :: RAS_NO_NEG  ! Whether RAS will guard against negatives

!  Get my name and set-up traceback handle
!  ---------------------------------------
   CALL ESMF_GridCompGet(GC, NAME=COMP_NAME, CONFIG=CF, GRID=grid, __RC__)
   Iam = TRIM(COMP_NAME)//"::Run_"

!  Get my internal MAPL_Generic state
!  -----------------------------------
   CALL MAPL_GetObjectFromGC(GC, ggState, __RC__)

!  Start a comprehensive timer
!  ---------------------------
   CALL MAPL_TimerOn(ggState, "TOTAL")

!  Get parameters from generic state
!  ---------------------------------
   CALL MAPL_Get(ggState, IM=i2, JM=j2, LM=km, LONS=LONS, LATS=LATS,  &
                 RUNALARM=ALARM, __RC__)

!  Query the alarm
!  ---------------
   RunGMINow = ESMF_AlarmIsRinging(ALARM, __RC__)

!  Get ESMF parameters from gc and clock
!  -------------------------------------
   CALL extract_(GC, clock, chemReg, gcGMI, w_c, nymd, nhms, gmiDt, runDt, STATUS)
   VERIFY_(STATUS)

  dtInverse = 1.00/runDt

!  Layer interface pressures
!  -------------------------
   CALL MAPL_GetPointer(impChem, PLE, 'PLE', __RC__)

!  Layer pressure thickness
!  ------------------------
   CALL MAPL_GetPointer(impChem, DELP, 'DELP', __RC__)
   w_c%delp = DELP

!  Fill in RH.  Note: Not converted to %
!  -------------------------------------
   CALL MAPL_GetPointer(impChem, rh2, 'RH2', __RC__)
   w_c%rh = rh2

!  Assure non-negative volumetric mixing ratios [mole fractions]
!  -------------------------------------------------------------
   CALL ESMF_ConfigGetAttribute( CF, RAS_NO_NEG, Label='RAS_NO_NEG:', default=.FALSE. , __RC__)
   IF ( .NOT. RAS_NO_NEG ) THEN
     DO n = ChemReg%i_GMI, ChemReg%j_GMI
! original approach - make all values no less than a tiny positive number:
       CALL Chem_UtilNegFiller(w_c%qa(n)%data3d, DELP, i2, j2, QMIN=TINY(1.0))
! debug print:
!      IF(  ANY(w_c%qa(n)%data3d < 1.0e-30) ) THEN
!        m = COUNT( w_c%qa(n)%data3d < 1.0e-30 )
!        PRINT*,'GMI SPECIES TOO SMALL (species,count):', n-ChemReg%i_GMI+1, m
!      ENDIF
! first take:
!     WHERE(w_c%qa(n)%data3d < 1.0e-30) w_c%qa(n)%data3d=1.0e-30
     END DO
   END IF

!  Occasionally, MAPL_UNDEFs appear in the imported tropopause pressures,
!  TROPP. To avoid encountering them, save the most recent valid tropopause 
!  pressures in an unused "layer" of T2M15D, i.e. w_c%qa(iT2M)%data3d(:,:,km). 
!  By passing the information through the internal state via w_c%qa (at least
!  for now), reproducibility across various layouts is assured, which is a 
!  requirement that Chem_UtilTroppFixer, used before Ganymed, cannot satisfy.
!  --------------------------------------------------------------------------
   CALL MAPL_GetPointer(impChem, TROPP, 'TROPP', __RC__)

   m = ChemReg%i_XX
   n = ChemReg%j_XX
   iT2M = -1
    
   DO i = m,n
    IF(TRIM(chemReg%vname(i)) == "T2M15d") iT2M = i
    IF(iT2M > 0) EXIT
   END DO

   IF(iT2M < 1) THEN
    PRINT *,TRIM(Iam)//": Invalid index for T2M15d (",iT2M,")"
    STATUS = 1
    VERIFY_(STATUS)
   END IF

   WHERE(tropp /= MAPL_UNDEF) w_c%qa(iT2M)%data3d(:,:,km) = TROPP

   IF( ANY(w_c%qa(iT2M)%data3d(:,:,km) == MAPL_UNDEF) ) THEN
    PRINT *,TRIM(Iam)//": At least one invalid tropopause pressure."
    STATUS = 1
    VERIFY_(STATUS)
   END IF

!  For comparison purposes, export both the "no MAPL_UNDEFs" TROPP and the
!  imported TROPP. Note that AGCM updates TROPP before HISTORY is written.
!  -----------------------------------------------------------------------
   CALL MAPL_GetPointer(expChem,  GMITROPP, 'GMITROPP',  __RC__)
   IF(ASSOCIATED( GMITROPP))  GMITROPP = w_c%qa(iT2M)%data3d(:,:,km)

   CALL MAPL_GetPointer(expChem, AGCMTROPP, 'AGCMTROPP', __RC__)
   IF(ASSOCIATED(AGCMTROPP)) AGCMTROPP = TROPP

! Are species tendencies requested?
! ---------------------------------
    m = ChemReg%i_GMI
    n = ChemReg%j_GMI
    ALLOCATE(doMyTendency(m:n), __STAT__)
    doMyTendency(:) = .FALSE.

    nAlloc = 0
    DO i = m,n

     fieldName = TRIM(chemReg%vname(i))
     incFieldName = TRIM(fieldName)//"_GMITEND"

     CALL MAPL_GetPointer(expChem, sIncrement, TRIM(incFieldName), __RC__)

     IF(ASSOCIATED(sIncrement)) THEN
      NULLIFY(sIncrement)
      nAlloc = nAlloc+1
      doMyTendency(i) = .TRUE.
     END IF

    END DO
    
    IF(nAlloc > 0) THEN
      doingTendencies = .TRUE.
    ELSE
      doingTendencies = .FALSE.
    ENDIF

!  Save current species configurations so chemical tendencies can be calculated.
!  NOTE: Restricted to transported species.
!  NOTE on PHASES:
!    phase  1: compute tendency and store in the export field
!    phase  2: compute tendency and add to   the export field
!    phase 99: compute tendency and store in the export field
!  ----------------------------------------------------------------------------
   StoreIC: IF(doingTendencies) THEN

    ALLOCATE(sInitial(1:i2,1:j2,km,nAlloc), __STAT__)

    k = 1
    m = ChemReg%i_GMI
    n = ChemReg%j_GMI
  
    DO i = m,n
     IF(doMyTendency(i)) THEN
      sInitial(:,:,:,k) = w_c%qa(i)%data3d
      k = k+1
     END IF
    END DO
   
   END IF StoreIC

! ... and for specific humidity
! -----------------------------
   IF ( phase == 2 .OR. phase == 99 ) THEN
     CALL MAPL_GetPointer(impChem, SPECHUM,       'Q', __RC__)
     CALL MAPL_GetPointer(expChem, Q_TEND, 'H2O_TEND', __RC__)
     IF(ASSOCIATED(Q_TEND)) Q_TEND = SPECHUM
   END IF

! ... and for OCS from ExtData or ACHEM
! -------------------------------------
   OCS: IF ( phase == 2 .OR. phase == 99 ) THEN

     CALL MAPL_GetPointer(impChem, OCS_import,  'OCS_CLIMO',  __RC__)
!    CALL MAPL_GetPointer(impChem, OCS_import,  'ACHEM::OCS', __RC__)

     m = ChemReg%i_XX
     n = ChemReg%j_XX
     iOCS = -1
    
     DO i = m,n
       IF(TRIM(chemReg%vname(i)) == "OCSg") THEN
         iOCS = i
         EXIT
       END IF
     END DO

     IF(iOCS < 1) THEN
      PRINT *,TRIM(Iam)//": Cannot find species OCSg in XX"
      STATUS = 1
      VERIFY_(STATUS)
     END IF

     w_c%qa(iOCS)%data3d(:,:,:) = OCS_import(:,:,:)

   END IF OCS

! Is replay running?
! ------------------
   doingPredictorNow = .FALSE.
   CALL ESMF_ConfigGetAttribute(CF, ReplayMode, Label='REPLAY_MODE:', DEFAULT="NoReplay", RC=STATUS)

   IF(ReplayMode == "Regular") THEN
    CALL ESMF_ClockGetAlarm(CLOCK, "PredictorActive", PredictorIsActive, RC=STATUS)

    IF(STATUS == 0) THEN
     doingPredictorNow = ESMF_AlarmIsRinging(PredictorIsActive, RC=STATUS)
     VERIFY_(STATUS)
     IF(MAPL_AM_I_ROOT()) PRINT *,TRIM(Iam)//": Replay predictor step detection: ",doingPredictorNow
    END IF

   END IF

! Pass to the underlying GCs, where needed
! ----------------------------------------
   gcGMI%gcEmiss%doingPredictorNow = doingPredictorNow



! At the Heartbeat do Run 1
! -------------------------
   Phase1: IF ( phase == 1 ) THEN

    IF(MAPL_AM_I_ROOT()) THEN
     PRINT *, " "
     PRINT *, TRIM(Iam)//": Running GMI Phase 1 (emissions) ..."
    END IF

    CALL MAPL_TimerOn( ggState, "RUN")

    CALL GMI_GridCompRun1(gcGMI, w_c, impChem, expChem, nymd, nhms, runDt, clock, STATUS)
    VERIFY_(STATUS)

    CALL MAPL_TimerOff(ggState, "RUN")

    IF(MAPL_AM_I_ROOT()) PRINT *, " "

    ! Also compute AOA (age of air) - see below

   END IF Phase1

! At the Heartbeat do Deposition, and at GMI timestep do the rest of chemistry
! ----------------------------------------------------------------------------
   Phase2: IF ( phase == 2 ) THEN

     IF(MAPL_AM_I_ROOT()) THEN
        PRINT *, " "
        IF(RunGMINow) THEN
          PRINT *, TRIM(Iam)//": Running GMI Phase 2 (deposition & chemistry) ..."
        ELSE
          PRINT *, TRIM(Iam)//": Running GMI Phase 2 (deposition only) ..."
        END IF
     END IF

     CALL MAPL_TimerOn(ggState, "RUN")

     CALL GMI_GridCompRun2(gcGMI, w_c, impChem, expChem, nymd, nhms, runDt, gmiDt, RunGMINow, STATUS)
     VERIFY_(STATUS)

     CALL MAPL_TimerOff(ggState, "RUN")

     IF(RunGMINow) CALL ESMF_AlarmRingerOff(ALARM, __RC__)

     IF(MAPL_AM_I_ROOT()) PRINT *, " "

     ! Also supply exports - see below
     ! Also compute OVP fields - see below

   END IF Phase2

! At GMI timestep do emissions and chemistry (old approach)
! ---------------------------------------------------------
   Phase99: IF ( phase == 99 ) THEN

     RunningGMI: IF(RunGMINow) THEN

       IF(MAPL_AM_I_ROOT()) THEN
          PRINT *, " "
          PRINT *, TRIM(Iam)//": Running GMI Phase 99 (emissions & chemistry) ..."
       END IF

       CALL MAPL_TimerOn(ggState, "RUN")

       CALL GMI_GridCompRunOrig(gcGMI, w_c, impChem, expChem, nymd, nhms, gmiDt, clock, STATUS)
       VERIFY_(STATUS)

       CALL MAPL_TimerOff(ggState, "RUN")

       CALL ESMF_AlarmRingerOff(ALARM, __RC__)

       IF(MAPL_AM_I_ROOT()) PRINT *, " "

     END IF RunningGMI

     ! Also compute AOA (age of air) - see below
     ! Also supply exports - see below
     ! Also compute OVP fields - see below

   END IF Phase99

!  Update age-of-air.
!  This transported species is at w_c%qa(ChemReg%i_GMI)%data3d.
!  This process is done at the Heartbeat (runDt)
!  --------------------------------------------------------------------------------
   IF ( phase == 1 .OR. phase == 99 ) THEN
     n = ChemReg%i_GMI
     w_c%qa(n)%data3d(:,:,:) = w_c%qa(n)%data3d(:,:,:)+runDt/86400.00
     w_c%qa(n)%data3d(:,:,km) = 0.00
   END IF

!  Gas-phase water in mole fraction.  Purpose: Allow plotting of mole 
!  fraction when using quickplot.  This avoids potential conflicts 
!  with existing directives for plotting specific humidity from MOIST. 
!  -------------------------------------------------------------------
   IF ( phase == 2 .OR. phase == 99 ) THEN
     CALL MAPL_GetPointer(expChem,   WATER, 'GMIH2O', __RC__)
     IF(ASSOCIATED(WATER)) WATER(:,:,:) = SPECHUM(:,:,:)*MAPL_AIRMW/MAPL_H2OMW
   END IF

!  Total ozone: In each layer
!   molecules m^{-2} =  O3(vmr) * Avogadros number * dp / ( mwt air * g )
!  The most recent valid tropopause pressures are stored in T2M15D(:,:,km)
!  -----------------------------------------------------------------------
   IF ( phase == 2 .OR. phase == 99 ) THEN

     CALL MAPL_GetPointer(expChem,  TO3,  'GMITO3', __RC__)
     IF(ASSOCIATED( TO3))  TO3 = 0.00

     CALL MAPL_GetPointer(expChem, TTO3, 'GMITTO3', __RC__)
     IF(ASSOCIATED(TTO3)) TTO3 = 0.00

     DoingTotalOzone: IF(ASSOCIATED(TTO3) .OR. ASSOCIATED(TO3)) THEN

      m = ChemReg%i_GMI
      n = ChemReg%j_GMI
      iOX = -1
    
      DO i = m,n
       IF(TRIM(chemReg%vname(i)) == "OX") iOX = i
       IF(iOx > 0) EXIT
      END DO

      IF(iOx < 1) THEN
       PRINT *,TRIM(Iam)//": Invalid index for Ox (",iOx,")"
       STATUS = 1
       VERIFY_(STATUS)
      END IF

      ALLOCATE(wrk(SIZE(LATS,1), SIZE(LATS,2)), __STAT__)
      ALLOCATE(wgt(SIZE(LATS,1), SIZE(LATS,2)), __STAT__)

      DO k = 1,km
       wrk = w_c%qa(iOX)%data3d(:,:,k)*(PLE(:,:,k)-PLE(:,:,k-1))*(MAPL_AVOGAD/2.69E+20)/(MAPL_AIRMW*MAPL_GRAV)
       IF(ASSOCIATED( TO3)) TO3 = TO3+wrk

       IF(ASSOCIATED(TTO3)) THEN
        wgt  = MAX(0.0,MIN(1.0,(PLE(:,:,k)-w_c%qa(iT2M)%data3d(:,:,km))/(PLE(:,:,k)-PLE(:,:,k-1))))
        TTO3 = TTO3+wrk*wgt
       END IF

      END DO

      DEALLOCATE(wrk)
      DEALLOCATE(wgt)

     END IF DoingTotalOzone
   END IF

!  Record the Overpass values
!  -------------------------------------------------------------------
   IF ( phase == 2 .OR. phase == 99 ) THEN

     CURRENT_HMS = OVP_end_of_timestep_hms( CLOCK, OVP_RUN_DT )
!    IF(MAPL_AM_I_ROOT()) PRINT*,'GMI CURRENT_HMS =',CURRENT_HMS

! NO2 overpass

     CALL secure_species_ptr( ChemReg, w_c, 'NO2', DATA_FOR_OVP_3D )

     CALL MAPL_GetPointer(expChem, OVP10_OUTPUT_3D, 'OVP10_NO2', __RC__)
     CALL MAPL_GetPointer(expChem, OVP14_OUTPUT_3D, 'OVP14_NO2', __RC__)

     CALL OVP_apply_mask( DATA_FOR_OVP_3D, OVP10_OUTPUT_3D, MASK_10AM, OVP_FIRST_HMS, CURRENT_HMS, K_EDGES=.FALSE., __RC__ )
     CALL OVP_apply_mask( DATA_FOR_OVP_3D, OVP14_OUTPUT_3D, MASK_2PM,  OVP_FIRST_HMS, CURRENT_HMS, K_EDGES=.FALSE., __RC__ )

! CH2O overpass

     CALL secure_species_ptr( ChemReg, w_c, 'CH2O', DATA_FOR_OVP_3D )

     CALL MAPL_GetPointer(expChem, OVP10_OUTPUT_3D, 'OVP10_CH2O', __RC__)
     CALL MAPL_GetPointer(expChem, OVP14_OUTPUT_3D, 'OVP14_CH2O', __RC__)

     CALL OVP_apply_mask( DATA_FOR_OVP_3D, OVP10_OUTPUT_3D, MASK_10AM, OVP_FIRST_HMS, CURRENT_HMS, K_EDGES=.FALSE., __RC__ )
     CALL OVP_apply_mask( DATA_FOR_OVP_3D, OVP14_OUTPUT_3D, MASK_2PM,  OVP_FIRST_HMS, CURRENT_HMS, K_EDGES=.FALSE., __RC__ )

! CO overpass

     CALL secure_species_ptr( ChemReg, w_c, 'CO', DATA_FOR_OVP_3D )

     CALL MAPL_GetPointer(expChem, OVP10_OUTPUT_3D, 'OVP10_CO', __RC__)
     CALL MAPL_GetPointer(expChem, OVP14_OUTPUT_3D, 'OVP14_CO', __RC__)

     CALL OVP_apply_mask( DATA_FOR_OVP_3D, OVP10_OUTPUT_3D, MASK_10AM, OVP_FIRST_HMS, CURRENT_HMS, K_EDGES=.FALSE., __RC__ )
     CALL OVP_apply_mask( DATA_FOR_OVP_3D, OVP14_OUTPUT_3D, MASK_2PM,  OVP_FIRST_HMS, CURRENT_HMS, K_EDGES=.FALSE., __RC__ )

! OX overpass

     CALL secure_species_ptr( ChemReg, w_c, 'OX', DATA_FOR_OVP_3D )

     CALL MAPL_GetPointer(expChem, OVP10_OUTPUT_3D, 'OVP10_OX', __RC__)
     CALL MAPL_GetPointer(expChem, OVP14_OUTPUT_3D, 'OVP14_OX', __RC__)

     CALL OVP_apply_mask( DATA_FOR_OVP_3D, OVP10_OUTPUT_3D, MASK_10AM, OVP_FIRST_HMS, CURRENT_HMS, K_EDGES=.FALSE., __RC__ )
     CALL OVP_apply_mask( DATA_FOR_OVP_3D, OVP14_OUTPUT_3D, MASK_2PM,  OVP_FIRST_HMS, CURRENT_HMS, K_EDGES=.FALSE., __RC__ )

! CH4 overpass

     CALL secure_species_ptr( ChemReg, w_c, 'CH4', DATA_FOR_OVP_3D )

     CALL MAPL_GetPointer(expChem, OVP10_OUTPUT_3D, 'OVP10_CH4', __RC__)
     CALL MAPL_GetPointer(expChem, OVP14_OUTPUT_3D, 'OVP14_CH4', __RC__)

     CALL OVP_apply_mask( DATA_FOR_OVP_3D, OVP10_OUTPUT_3D, MASK_10AM, OVP_FIRST_HMS, CURRENT_HMS, K_EDGES=.FALSE., __RC__ )
     CALL OVP_apply_mask( DATA_FOR_OVP_3D, OVP14_OUTPUT_3D, MASK_2PM,  OVP_FIRST_HMS, CURRENT_HMS, K_EDGES=.FALSE., __RC__ )

! BrO overpass

     CALL secure_species_ptr( ChemReg, w_c, 'BrO', DATA_FOR_OVP_3D )

     CALL MAPL_GetPointer(expChem, OVP14_OUTPUT_3D, 'OVP14_BrO', __RC__)
     CALL OVP_apply_mask( DATA_FOR_OVP_3D, OVP14_OUTPUT_3D, MASK_2PM,  OVP_FIRST_HMS, CURRENT_HMS, K_EDGES=.FALSE., __RC__ )

! PAN overpass

     CALL secure_species_ptr( ChemReg, w_c, 'PAN', DATA_FOR_OVP_3D )

     CALL MAPL_GetPointer(expChem, OVP14_OUTPUT_3D, 'OVP14_PAN', __RC__)
     CALL OVP_apply_mask( DATA_FOR_OVP_3D, OVP14_OUTPUT_3D, MASK_2PM,  OVP_FIRST_HMS, CURRENT_HMS, K_EDGES=.FALSE., __RC__ )

! ISOP overpass

     CALL secure_species_ptr( ChemReg, w_c, 'ISOP', DATA_FOR_OVP_3D )

     CALL MAPL_GetPointer(expChem, OVP14_OUTPUT_3D, 'OVP14_ISOP', __RC__)
     CALL OVP_apply_mask( DATA_FOR_OVP_3D, OVP14_OUTPUT_3D, MASK_2PM,  OVP_FIRST_HMS, CURRENT_HMS, K_EDGES=.FALSE., __RC__ )

! MOH overpass

     CALL secure_species_ptr( ChemReg, w_c, 'MOH', DATA_FOR_OVP_3D )

     CALL MAPL_GetPointer(expChem, OVP14_OUTPUT_3D, 'OVP14_MOH', __RC__)
     CALL OVP_apply_mask( DATA_FOR_OVP_3D, OVP14_OUTPUT_3D, MASK_2PM,  OVP_FIRST_HMS, CURRENT_HMS, K_EDGES=.FALSE., __RC__ )

! HCOOH overpass

     CALL secure_species_ptr( ChemReg, w_c, 'HCOOH', DATA_FOR_OVP_3D )

     CALL MAPL_GetPointer(expChem, OVP14_OUTPUT_3D, 'OVP14_HCOOH', __RC__)
     CALL OVP_apply_mask( DATA_FOR_OVP_3D, OVP14_OUTPUT_3D, MASK_2PM,  OVP_FIRST_HMS, CURRENT_HMS, K_EDGES=.FALSE., __RC__ )

! N2O overpass

     CALL secure_species_ptr( ChemReg, w_c, 'N2O', DATA_FOR_OVP_3D )

     CALL MAPL_GetPointer(expChem, OVP14_OUTPUT_3D, 'OVP14_N2O', __RC__)
     CALL OVP_apply_mask( DATA_FOR_OVP_3D, OVP14_OUTPUT_3D, MASK_2PM,  OVP_FIRST_HMS, CURRENT_HMS, K_EDGES=.FALSE., __RC__ )

! HCl overpass

     CALL secure_species_ptr( ChemReg, w_c, 'HCl', DATA_FOR_OVP_3D )

     CALL MAPL_GetPointer(expChem, OVP14_OUTPUT_3D, 'OVP14_HCl', __RC__)
     CALL OVP_apply_mask( DATA_FOR_OVP_3D, OVP14_OUTPUT_3D, MASK_2PM,  OVP_FIRST_HMS, CURRENT_HMS, K_EDGES=.FALSE., __RC__ )

! HNO3 overpass

     CALL secure_species_ptr( ChemReg, w_c, 'HNO3', DATA_FOR_OVP_3D )

     CALL MAPL_GetPointer(expChem, OVP14_OUTPUT_3D, 'OVP14_HNO3', __RC__)
     CALL OVP_apply_mask( DATA_FOR_OVP_3D, OVP14_OUTPUT_3D, MASK_2PM,  OVP_FIRST_HMS, CURRENT_HMS, K_EDGES=.FALSE., __RC__ )

! ClO overpass

     CALL secure_species_ptr( ChemReg, w_c, 'ClO', DATA_FOR_OVP_3D )

     CALL MAPL_GetPointer(expChem, OVP14_OUTPUT_3D, 'OVP14_ClO', __RC__)
     CALL OVP_apply_mask( DATA_FOR_OVP_3D, OVP14_OUTPUT_3D, MASK_2PM,  OVP_FIRST_HMS, CURRENT_HMS, K_EDGES=.FALSE., __RC__ )

! CH3Cl overpass

     CALL secure_species_ptr( ChemReg, w_c, 'CH3Cl', DATA_FOR_OVP_3D )

     CALL MAPL_GetPointer(expChem, OVP14_OUTPUT_3D, 'OVP14_CH3Cl', __RC__)
     CALL OVP_apply_mask( DATA_FOR_OVP_3D, OVP14_OUTPUT_3D, MASK_2PM,  OVP_FIRST_HMS, CURRENT_HMS, K_EDGES=.FALSE., __RC__ )

! NO overpass

     CALL secure_species_ptr( ChemReg, w_c, 'NO', DATA_FOR_OVP_3D )

     CALL MAPL_GetPointer(expChem, OVP10_OUTPUT_3D, 'OVP10_NO', __RC__)
     CALL MAPL_GetPointer(expChem, OVP14_OUTPUT_3D, 'OVP14_NO', __RC__)

     CALL OVP_apply_mask( DATA_FOR_OVP_3D, OVP10_OUTPUT_3D, MASK_10AM, OVP_FIRST_HMS, CURRENT_HMS, K_EDGES=.FALSE., __RC__ )
     CALL OVP_apply_mask( DATA_FOR_OVP_3D, OVP14_OUTPUT_3D, MASK_2PM,  OVP_FIRST_HMS, CURRENT_HMS, K_EDGES=.FALSE., __RC__ )

! HOCl overpass

     CALL secure_species_ptr( ChemReg, w_c, 'HOCl', DATA_FOR_OVP_3D )

     CALL MAPL_GetPointer(expChem, OVP14_OUTPUT_3D, 'OVP14_HOCl', __RC__)
     CALL OVP_apply_mask( DATA_FOR_OVP_3D, OVP14_OUTPUT_3D, MASK_2PM,  OVP_FIRST_HMS, CURRENT_HMS, K_EDGES=.FALSE., __RC__ )

! ClONO2 overpass

     CALL secure_species_ptr( ChemReg, w_c, 'ClONO2', DATA_FOR_OVP_3D )

     CALL MAPL_GetPointer(expChem, OVP14_OUTPUT_3D, 'OVP14_ClONO2', __RC__)
     CALL OVP_apply_mask( DATA_FOR_OVP_3D, OVP14_OUTPUT_3D, MASK_2PM,  OVP_FIRST_HMS, CURRENT_HMS, K_EDGES=.FALSE., __RC__ )

! Cl2O2 overpass

     CALL secure_species_ptr( ChemReg, w_c, 'Cl2O2', DATA_FOR_OVP_3D )

     CALL MAPL_GetPointer(expChem, OVP14_OUTPUT_3D, 'OVP14_Cl2O2', __RC__)
     CALL OVP_apply_mask( DATA_FOR_OVP_3D, OVP14_OUTPUT_3D, MASK_2PM,  OVP_FIRST_HMS, CURRENT_HMS, K_EDGES=.FALSE., __RC__ )

! Cl overpass

     CALL secure_species_ptr( ChemReg, w_c, 'Cl', DATA_FOR_OVP_3D )

     CALL MAPL_GetPointer(expChem, OVP14_OUTPUT_3D, 'OVP14_Cl', __RC__)
     CALL OVP_apply_mask( DATA_FOR_OVP_3D, OVP14_OUTPUT_3D, MASK_2PM,  OVP_FIRST_HMS, CURRENT_HMS, K_EDGES=.FALSE., __RC__ )

! Cl2 overpass

     CALL secure_species_ptr( ChemReg, w_c, 'Cl2', DATA_FOR_OVP_3D )

     CALL MAPL_GetPointer(expChem, OVP14_OUTPUT_3D, 'OVP14_Cl2', __RC__)
     CALL OVP_apply_mask( DATA_FOR_OVP_3D, OVP14_OUTPUT_3D, MASK_2PM,  OVP_FIRST_HMS, CURRENT_HMS, K_EDGES=.FALSE., __RC__ )

! OClO overpass

     CALL secure_species_ptr( ChemReg, w_c, 'OClO', DATA_FOR_OVP_3D )

     CALL MAPL_GetPointer(expChem, OVP14_OUTPUT_3D, 'OVP14_OClO', __RC__)
     CALL OVP_apply_mask( DATA_FOR_OVP_3D, OVP14_OUTPUT_3D, MASK_2PM,  OVP_FIRST_HMS, CURRENT_HMS, K_EDGES=.FALSE., __RC__ )

! BrCl overpass

     CALL secure_species_ptr( ChemReg, w_c, 'BrCl', DATA_FOR_OVP_3D )

     CALL MAPL_GetPointer(expChem, OVP14_OUTPUT_3D, 'OVP14_BrCl', __RC__)
     CALL OVP_apply_mask( DATA_FOR_OVP_3D, OVP14_OUTPUT_3D, MASK_2PM,  OVP_FIRST_HMS, CURRENT_HMS, K_EDGES=.FALSE., __RC__ )

! Br overpass

     CALL secure_species_ptr( ChemReg, w_c, 'Br', DATA_FOR_OVP_3D )

     CALL MAPL_GetPointer(expChem, OVP14_OUTPUT_3D, 'OVP14_Br', __RC__)
     CALL OVP_apply_mask( DATA_FOR_OVP_3D, OVP14_OUTPUT_3D, MASK_2PM,  OVP_FIRST_HMS, CURRENT_HMS, K_EDGES=.FALSE., __RC__ )

! BrONO2 overpass

     CALL secure_species_ptr( ChemReg, w_c, 'BrONO2', DATA_FOR_OVP_3D )

     CALL MAPL_GetPointer(expChem, OVP14_OUTPUT_3D, 'OVP14_BrONO2', __RC__)
     CALL OVP_apply_mask( DATA_FOR_OVP_3D, OVP14_OUTPUT_3D, MASK_2PM,  OVP_FIRST_HMS, CURRENT_HMS, K_EDGES=.FALSE., __RC__ )

! HBr overpass

     CALL secure_species_ptr( ChemReg, w_c, 'HBr', DATA_FOR_OVP_3D )

     CALL MAPL_GetPointer(expChem, OVP14_OUTPUT_3D, 'OVP14_HBr', __RC__)
     CALL OVP_apply_mask( DATA_FOR_OVP_3D, OVP14_OUTPUT_3D, MASK_2PM,  OVP_FIRST_HMS, CURRENT_HMS, K_EDGES=.FALSE., __RC__ )

! HOBr overpass

     CALL secure_species_ptr( ChemReg, w_c, 'HOBr', DATA_FOR_OVP_3D )

     CALL MAPL_GetPointer(expChem, OVP14_OUTPUT_3D, 'OVP14_HOBr', __RC__)
     CALL OVP_apply_mask( DATA_FOR_OVP_3D, OVP14_OUTPUT_3D, MASK_2PM,  OVP_FIRST_HMS, CURRENT_HMS, K_EDGES=.FALSE., __RC__ )

! OH overpass

     CALL secure_species_ptr( ChemReg, w_c, 'OH', DATA_FOR_OVP_3D )

     CALL MAPL_GetPointer(expChem, OVP14_OUTPUT_3D, 'OVP14_OH', __RC__)
     CALL OVP_apply_mask( DATA_FOR_OVP_3D, OVP14_OUTPUT_3D, MASK_2PM,  OVP_FIRST_HMS, CURRENT_HMS, K_EDGES=.FALSE., __RC__ )

! HO2 overpass

     CALL secure_species_ptr( ChemReg, w_c, 'HO2', DATA_FOR_OVP_3D )

     CALL MAPL_GetPointer(expChem, OVP14_OUTPUT_3D, 'OVP14_HO2', __RC__)
     CALL OVP_apply_mask( DATA_FOR_OVP_3D, OVP14_OUTPUT_3D, MASK_2PM,  OVP_FIRST_HMS, CURRENT_HMS, K_EDGES=.FALSE., __RC__ )

! R4N2 overpass

     CALL secure_species_ptr( ChemReg, w_c, 'R4N2', DATA_FOR_OVP_3D )

     CALL MAPL_GetPointer(expChem, OVP14_OUTPUT_3D, 'OVP14_R4N2', __RC__)
     CALL OVP_apply_mask( DATA_FOR_OVP_3D, OVP14_OUTPUT_3D, MASK_2PM,  OVP_FIRST_HMS, CURRENT_HMS, K_EDGES=.FALSE., __RC__ )

! HNO3COND overpass

     CALL secure_species_ptr( ChemReg, w_c, 'HNO3COND', DATA_FOR_OVP_3D )

     CALL MAPL_GetPointer(expChem, OVP14_OUTPUT_3D, 'OVP14_HNO3COND', __RC__)
     CALL OVP_apply_mask( DATA_FOR_OVP_3D, OVP14_OUTPUT_3D, MASK_2PM,  OVP_FIRST_HMS, CURRENT_HMS, K_EDGES=.FALSE., __RC__ )

   END IF



!  Obtain chemical tendencies and fill export states.  NOTE: Restricted to transported species.
!  NOTE on PHASES:
!    phase  1: compute tendency and store in the export field
!    phase  2: compute tendency and add to   the export field
!    phase 99: compute tendency and store in the export field
!  --------------------------------------------------------------------------------------------
   StoreTendencies: IF(doingTendencies) THEN

    k = 1
    m = ChemReg%i_GMI
    n = ChemReg%j_GMI

    DO i = m,n

     IF(doMyTendency(i)) THEN

      fieldName = TRIM(chemReg%vname(i))
      incFieldName = TRIM(fieldName)//"_GMITEND"

      CALL MAPL_GetPointer(expChem, sIncrement, TRIM(incFieldName), __RC__)

      IF(ASSOCIATED(sIncrement)) THEN
        IF ( phase == 1 .OR. phase == 99 ) THEN
          sIncrement =              (w_c%qa(i)%data3d - sInitial(:,:,:,k))*dtInverse
        END IF
        IF ( phase == 2 ) THEN
          sIncrement = sIncrement + (w_c%qa(i)%data3d - sInitial(:,:,:,k))*dtInverse
        END IF
      END IF

      NULLIFY(sIncrement)
      k = k+1

     END IF

    END DO

!  Clean up
!  --------
    DEALLOCATE(sInitial, STAT=STATUS)
    VERIFY_(STATUS)
   
   END IF StoreTendencies

   DEALLOCATE(doMyTendency, STAT=STATUS)
   VERIFY_(STATUS)

! ... and for specific humidity [kg kg^{-1} s^{-1}]
! -------------------------------------------------
   IF ( phase == 2 .OR. phase == 99 ) THEN
     IF(ASSOCIATED(Q_TEND)) Q_TEND = (SPECHUM-Q_TEND)*dtInverse
   END IF

!  Stop timer
!  ----------
   CALL MAPL_TimerOff(ggState, "TOTAL")

   RETURN_(ESMF_SUCCESS)

   END SUBROUTINE Run_

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Finalize_ --- Finalize GMIchem_GridComp (ESMF)
!
! !INTERFACE:
!

   SUBROUTINE Finalize_ ( gc, impChem, expChem, clock, rc )

! !USES:

  implicit NONE

! !INPUT PARAMETERS:

   type(ESMF_Clock),  intent(inout) :: clock      ! The clock

! !OUTPUT PARAMETERS:

   type(ESMF_GridComp), intent(inout)  :: gc      ! Grid Component
   type(ESMF_State), intent(inout) :: impChem     ! Import State
   type(ESMF_State), intent(inout) :: expChem     ! Export State
   integer, intent(out) ::  rc                    ! Error return code:
                                                  !  0 - all is well
                                                  !  1 - 

! !DESCRIPTION: This is a simple ESMF wrapper.
!
! !REVISION HISTORY:
!
!  27Feb2005 da Silva  First crack.
!
!EOP
!-------------------------------------------------------------------------


!  ErrLog Variables
!  ----------------
   character(len=ESMF_MAXSTR)      :: IAm
   integer                         :: STATUS
   character(len=ESMF_MAXSTR)      :: COMP_NAME

   type(Chem_Registry), pointer    :: chemReg
   type(GMI_GridComp), pointer     :: gcGMI       ! Grid Component
   type(Chem_Bundle), pointer      :: w_c         ! Chemical tracer fields     
   type(MAPL_MetaComp), pointer    :: ggState     ! GEOS Generic State
   integer                         :: nymd, nhms  ! time
   real                            :: gmiDt       ! chemistry timestep (secs)
   real                            :: runDt       ! heartbeat (secs)

    type(GMIchem_state), pointer   :: state

!  Get my name and set-up traceback handle
!  ---------------------------------------
   call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
   VERIFY_(STATUS)
   Iam = TRIM(COMP_NAME)//"Finalize_"

!  Get my internal MAPL_Generic state
!  -----------------------------------
   CALL MAPL_GetObjectFromGC(GC, ggState, RC=STATUS)

!  Start timers
!  ------------
   CALL MAPL_TimerOn(ggState, "TOTAL")
   CALL MAPL_TimerOn(ggState, "FINALIZE")

!  Get ESMF parameters from gc and clock
!  -------------------------------------
   call extract_(gc, clock, chemReg, gcGMI, w_c, nymd, nhms, gmiDt, runDt, STATUS, state = state)
   VERIFY_(STATUS)

!  Call ESMF version
!  -----------------
   call GMI_GridCompFinalize(gcGMI, w_c, impChem, expChem, nymd, nhms, gmiDt, STATUS)
   VERIFY_(STATUS)

!  Destroy Chem_Bundle
!  -------------------
   call Chem_BundleDestroy ( w_c, STATUS )
   VERIFY_(STATUS)

!  Destroy Chem_Registry
!  ---------------------
   call Chem_RegistryDestroy ( chemReg, STATUS ) 
   VERIFY_(STATUS)

!  Destroy Legacy state
!  --------------------
   deallocate ( state%chemReg, state%gcGMI, state%w_c, stat = STATUS )
   VERIFY_(STATUS)

!  Free the masks
!  --------------------
   deallocate ( MASK_10AM, MASK_2PM, stat = STATUS )
   VERIFY_(STATUS)

!  Stop timers
!  -----------
   CALL MAPL_TimerOff(ggState, "FINALIZE")
   CALL MAPL_TimerOff(ggState, "TOTAL")

!  Finalize MAPL Generic
!  ---------------------
   call MAPL_GenericFinalize ( gc, impChem, expChem, clock,  RC=STATUS )
   VERIFY_(STATUS)

   RETURN_(ESMF_SUCCESS)

   END SUBROUTINE Finalize_

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
    SUBROUTINE extract_(gc, clock, chemReg, gcGMI, w_c, nymd, nhms, gmiDt, runDt, rc, state)

    type(ESMF_GridComp), intent(inout) :: gc
    type(ESMF_Clock), intent(in)       :: clock
    type(Chem_Registry), pointer       :: chemReg
    type(GMI_GridComp), pointer        :: gcGMI
    type(Chem_Bundle), pointer         :: w_c
    integer, intent(out)               :: nymd, nhms
    real, intent(out)                  :: gmiDt
    real, intent(out)                  :: runDt
    integer, intent(out)               :: rc
    type(GMIchem_state), pointer, optional   :: state


    type(GMIchem_state), pointer    :: myState

!   ErrLog Variables
!   ----------------
    character(len=ESMF_MAXSTR)      :: IAm
    integer                         :: STATUS
    character(len=ESMF_MAXSTR)      :: COMP_NAME

    type(ESMF_Time)      :: TIME
    type(ESMF_Config)    :: CF
    type(GMIchem_Wrap)   :: wrap
    integer              :: IYR, IMM, IDD, IHR, IMN, ISC


!   Get my name and set-up traceback handle
!   ---------------------------------------
    call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // 'extract_'

    rc = 0

!   Get my internal state
!   ---------------------
    call ESMF_UserCompGetInternalState(gc, 'GMIchem_state', WRAP, STATUS)
    VERIFY_(STATUS)
    myState => wrap%ptr
    if ( present(state) ) then
         state => wrap%ptr
    end if

!   This is likely to be allocated during initialize only
!   -----------------------------------------------------
    if ( .not. associated(myState%chemReg) ) then
         allocate ( myState%chemReg, stat=STATUS )
         VERIFY_(STATUS)
    end if
    if ( .not. associated(myState%gcGMI) ) then
         allocate ( myState%gcGMI, stat=STATUS )
         VERIFY_(STATUS)
    end if
    if ( .not. associated(myState%w_c) ) then
         allocate ( myState%w_c, stat=STATUS )
         VERIFY_(STATUS)
    end if

    chemReg => myState%chemReg
    gcGMI  => myState%gcGMI
    w_c     => myState%w_c

!   Get the configuration
!   ---------------------
    call ESMF_GridCompGet ( GC, CONFIG = CF, RC=STATUS )
    VERIFY_(STATUS)

!   Get GEOS-5 time step
!   --------------------
    call ESMF_ConfigGetAttribute ( CF, runDt, LABEL="RUN_DT:", RC=STATUS )
    VERIFY_(STATUS)

!   Chemistry time step can be longer if GMICHEM_DT is set in AGCM.rc
!   -----------------------------------------------------------------
    CALL ESMF_ConfigGetAttribute ( CF, gmiDt, LABEL="GMICHEM_DT:", DEFAULT=runDt, RC=STATUS )
    VERIFY_(STATUS)
    
    IF(gmiDt < runDt) THEN
     IF(MAPL_AM_I_ROOT()) PRINT *,"GMICHEM_DT cannot be less than RUN_DT"
     STATUS = 1
     VERIFY_(STATUS)
    END IF

!   Extract nymd, nhms, day of year from clock
!   ------------------------------------------
    call ESMF_ClockGet(CLOCK,currTIME=TIME,rc=STATUS)
    VERIFY_(STATUS)

    call ESMF_TimeGet(TIME ,YY=IYR, MM=IMM, DD=IDD, H=IHR, M=IMN, S=ISC, rc=STATUS)
    VERIFY_(STATUS)

    call MAPL_PackTime(NYMD,IYR,IMM,IDD)
    call MAPL_PackTime(NHMS,IHR,IMN,ISC)

    RETURN_(ESMF_SUCCESS)

   END SUBROUTINE extract_

subroutine aerosol_optics(state, rc)

  implicit none

! Arguments
! ---------
  type(ESMF_State)     :: state
  integer, intent(out) :: rc


! Local
! ---------
  integer                                 :: n_aerosols
  character(len=ESMF_MAXSTR), allocatable :: aerosol_names(:)
  type(ESMF_FieldBundle)                  :: aerosols

  real, dimension(:,:,:), pointer         :: ple
  real, dimension(:,:,:), pointer         :: rh
  real, dimension(:,:,:), pointer         :: var
  real, dimension(:,:,:), pointer         :: q
  real, dimension(:,:,:,:), pointer       :: q_4d

  real, dimension(:,:,:), allocatable     :: dp, f_p

  character(len=ESMF_MAXSTR)              :: fld_name
  type(ESMF_Field)                        :: fld

  real, dimension(:,:,:,:), allocatable   :: ext, ssa, asy  ! (lon:,lat:,lev:,band:)

  integer                                 :: n
  integer                                 :: i1, j1, i2, j2, km

  integer                                 :: band, offset

  integer                                 :: instance

  integer                                 :: STATUS
  character(len=ESMF_MAXSTR)              :: Iam

  integer, parameter                      :: n_bands = 1

  real    :: x
  integer :: i, j, k

  Iam = 'GMI::aerosol_optics()'


! Mie Table instance/index
! ------------------------
  call ESMF_AttributeGet(state, name='mie_table_instance', value=instance, __RC__)

! Radiation band
! --------------
  band = 0
  call ESMF_AttributeGet(state, name='band_for_aerosol_optics', value=band, __RC__)
  offset = band - n_bands

! Pressure at layer edges 
! ------------------------
  call ESMF_AttributeGet(state, name='air_pressure_for_aerosol_optics', value=fld_name, __RC__)
  call MAPL_GetPointer(state, ple, trim(fld_name), __RC__)

  i1 = lbound(ple, 1); i2 = ubound(ple, 1)
  j1 = lbound(ple, 2); j2 = ubound(ple, 2)
                       km = ubound(ple, 3)

! Relative humidity
! -----------------
  call ESMF_AttributeGet(state, name='relative_humidity_for_aerosol_optics', value=fld_name, __RC__)
  call MAPL_GetPointer(state, rh, trim(fld_name), __RC__)

  i1 = lbound(rh, 1); i2 = ubound(rh, 1)
  j1 = lbound(rh, 2); j2 = ubound(rh, 2)
                      km = ubound(rh, 3)
  
  call ESMF_StateGet(state, 'AEROSOLS', aerosols, __RC__)
  call ESMF_FieldBundleGet(aerosols, fieldCount=n_aerosols, __RC__)

  allocate(aerosol_names(n_aerosols), __STAT__)
 
  call ESMF_FieldBundleGet(aerosols, FieldNameList=aerosol_names, __RC__)
 
  allocate(ext(i1:i2,j1:j2,km,n_bands), &
           ssa(i1:i2,j1:j2,km,n_bands), &
           asy(i1:i2,j1:j2,km,n_bands), __STAT__)

  allocate(q_4d(i1:i2,j1:j2,km,n_aerosols), __STAT__)

#if (0)
  allocate(dp(i1:i2,j1:j2,km), f_p(i1:i2,j1:j2,km), __STAT__)

  dp  = ple(:,:,1:km) - ple(:,:,0:km-1)
  f_p = dp / MAPL_GRAV

  do n = 1, n_aerosols
      call ESMF_FieldBundleGet(aerosols, trim(aerosol_names(n)), field=fld, __RC__)
      call ESMF_FieldGet(fld, farrayPtr=q, __RC__)

      q_4d(:,:,:,n) = f_p * q
  end do

  call ESMF_AttributeGet(state, name='mie_table_instance', value=instance, __RC__)
  call mie_(gocartMieTable(instance, aerosol_names, n_bands, offset, q_4d, rh, ext, ssa, asy, __RC__)

  deallocate(dp, f_p, __STAT__)
#else
  do n = 1, n_aerosols
      call ESMF_FieldBundleGet(aerosols, trim(aerosol_names(n)), field=fld, __RC__)
      call ESMF_FieldGet(fld, farrayPtr=q, __RC__)

      do k = 1, km
          do j = j1, j2
              do i = i1, i2
                  x = ((PLE(i,j,k) - PLE(i,j,k-1))*0.01)*(100./MAPL_GRAV)
                  q_4d(i,j,k,n) = x * q(i,j,k)
              end do
          end do
      end do
  end do

  call mie_(gocartMieTable(instance), aerosol_names, n_bands, offset, q_4d, rh, ext, ssa, asy, __RC__)
#endif
  
  call ESMF_AttributeGet(state, name='extinction_in_air_due_to_ambient_aerosol', value=fld_name, __RC__)
  if (fld_name /= '') then 
      call MAPL_GetPointer(state, var, trim(fld_name), __RC__)
      var = ext(:,:,:,1)
  end if

  call ESMF_AttributeGet(state, name='single_scattering_albedo_of_ambient_aerosol', value=fld_name, __RC__)
  if (fld_name /= '') then 
      call MAPL_GetPointer(state, var, trim(fld_name), __RC__)
      var = ssa(:,:,:,1)
  end if

  call ESMF_AttributeGet(state, name='asymmetry_parameter_of_ambient_aerosol', value=fld_name, __RC__)
  if (fld_name /= '') then 
      call MAPL_GetPointer(state, var, trim(fld_name), __RC__)
      var = asy(:,:,:,1)
  end if

  deallocate(aerosol_names, ext, ssa, asy, q_4d, __STAT__)

  RETURN_(ESMF_SUCCESS)

contains 

    subroutine mie_(mie_table, aerosol, nb, offset, q, rh, ext, ssa, asy, rc)
     
     implicit none

     type(Chem_Mie),    intent(inout):: mie_table    ! mie table
     character(len=*),  intent(in )  :: aerosol(:)   ! list of aerosols
     integer,           intent(in )  :: nb           ! number of bands
     integer,           intent(in )  :: offset       ! bands offset 
     real,              intent(in )  :: q(:,:,:,:)   ! aerosol mass mixing ratio, kg kg-1
     real,              intent(in )  :: rh(:,:,:)    ! relative humidity

     real,              intent(out)  :: ext(:,:,:,:) ! extinction
     real,              intent(out)  :: ssa(:,:,:,:) ! SSA
     real,              intent(out)  :: asy(:,:,:,:) ! asymmetry parameter

     integer,           intent(out)  :: rc

     ! local
     integer :: STATUS
     character(len=ESMF_MAXSTR) :: Iam='aerosol_optics::mie_' 

     integer :: l, idx, na

     real(kind=8) :: ext_(size(ext,1),size(ext,2),size(ext,3),size(ext,4))
     real(kind=8) :: ssa_(size(ext,1),size(ext,2),size(ext,3),size(ext,4))
     real(kind=8) :: asy_(size(ext,1),size(ext,2),size(ext,3),size(ext,4))

     na = size(aerosol)

     _ASSERT(na == size(q,4), 'needs informative message')

     ext_ = 0.0d0
     ssa_ = 0.0d0
     asy_ = 0.0d0

     do l = 1, na
        idx = Chem_MieQueryIdx(mie_table, trim(aerosol(l)), __RC__)

        call Chem_MieQueryAllBand4D(mie_table, idx, nb, offset, q(:,:,:,l), rh, ext, ssa, asy, __RC__)

        ext_ = ext_ +          ext     ! total extinction
        ssa_ = ssa_ +     (ssa*ext)    ! total scattering
        asy_ = asy_ + asy*(ssa*ext)    ! sum of (asy * sca)
     end do

     ext = ext_
     ssa = ssa_
     asy = asy_

     RETURN_(ESMF_SUCCESS)

    end subroutine mie_

end subroutine aerosol_optics

SUBROUTINE secure_species_ptr( ChemReg, w_c, SPECIES_NAME, DATA_PTR )

! !ARGUMENTS

     TYPE(Chem_Registry), POINTER,    INTENT(in)  :: chemReg
     TYPE(Chem_Bundle),   POINTER,    INTENT(in)  :: w_c
     CHARACTER(len=*),                INTENT(in)  :: SPECIES_NAME
     REAL, POINTER, DIMENSION(:,:,:), INTENT(out) :: DATA_PTR

! Locals

     INTEGER                         :: i, iSpecies

    iSpecies = -1

    DO i = ChemReg%i_GMI,ChemReg%j_GMI
     IF(TRIM(chemReg%vname(i)) == SPECIES_NAME) iSpecies = i
!    IF(iSpecies > 0) EXIT
    END DO

    IF(iSpecies < 1) THEN
      NULLIFY(DATA_PTR)
    ELSE
      DATA_PTR => w_c%qa(iSpecies)%data3d
    END IF

    RETURN

END SUBROUTINE secure_species_ptr

 END MODULE GMIchem_GridCompMod
