#include "MAPL_Generic.h"

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: GEOS_AChemGridCompMod -
!
! !INTERFACE:
!
   module GEOS_AChemGridCompMod
!
! !USES:
!
   use ESMF
   use MAPL

   use m_StrTemplate,         only: StrTemplate

   use DryDepositionMod,      only: DryDepositionGOCART

   use GACL_ConstantsMod,     only: pi, g_earth, N_avog, R_univ, &
                                    mw_air, mw_S, mw_SO2, mw_SO4, mw_H2SO4, &
                                    mw_DMS, mw_MSA, mw_OH, mw_NO3, mw_N, mw_NH3, mw_NH4, mw_SOAg

   use GACL_DryDepositionMod, only: DepositionVelocity

   use GACL_EmissionsMod,     only: NH3_Emissions,  &
                                    SO2_Emissions,  &
                                    DMS_Emissions,  &
                                    SOAG_Emissions, &
                                    VOC_Emissions


   use GACL_ReactionRatesMod, only: henry,                &
                                    H_SO2_298,  E_R_SO2,  &
                                    H_NH3_298,  E_R_NH3,  &
                                    H_H2O2_298, E_R_H2O2, &
                                    H_O3_298,   E_R_O3


   implicit none
   private
!
! !PUBLIC MEMBER FUNCTIONS:

   public SetServices
!
! !DESCRIPTION:
!
!  {\tt GEOS\_AChem} is an ESMF gridded component implementing gas and aqueous phase
!  chemistry in GEOS-5.
!
!  Developed for GEOS-5 release Fortuna 2.0 and later.
!
! !REVISION HISTORY:
!
!  08Aug2012  A. Darmenov  Cloned from MAM
!
!EOP
!-------------------------------------------------------------------------

! Legacy state
! ------------
  type AChem_State
     private
     type(ESMF_Config)   :: CF                        ! Private Config

     type(ESMF_Grid)     :: grid                      ! Grid

     logical             :: verbose                   ! turn on/off more verbose messages

     logical             :: mam_chem                  ! aerosol chemistry for MAM and alike
     logical             :: gas_phase_chem            ! enable/disable gas phase chemistry
     logical             :: aqu_phase_chem            ! enable/disable aqueous phase chemistry

     logical             :: ocs_chem                  ! enable/disable OCS chemistry mechanism
     real                :: ocs_surface_vmr = 0.0     ! OCS surface volume mixing ratio

     logical             :: voc_chem                  ! voc chemistry    ! turn on/off VOCs
     real                :: voc_BiomassBurnFactor = 0.0 ! conversion factor CO->VOC (BB)
     real                :: voc_AnthroFactor = 0.0    ! conversion factor CO->VOC (anthro)
     real                :: voc_MW = 0.0              ! molecular weight of VOC
     real                :: soa_MW = 0.0              ! molecular weight of SOA

     real                :: aqu_solver_max_dt         ! maximum time step used for integration in the aqueous-phase mechanism

     logical             :: apply_diurnal_cycle       ! flag that indicates if offline oxidant have to be temporally downscaled

     real, dimension(4)  :: aviation_layers           ! heights of the LTO, CDS and CRS layers

     integer             :: nymd_volcanic_emiss       ! nYMD of last volcanic emission update
     character(len=1024) :: volcanic_emiss_file       ! resource file with volcanic emissions data
     integer             :: n_volcanoes = 0           ! point wise location, amount, elevation, plume height, cell indexes
     real, pointer, dimension(:)    :: volc_lat   => null(), &
                                       volc_lon   => null(), &
                                       volc_SO2   => null(), &
                                       volc_elev  => null(), &
                                       volc_cloud => null()
     integer, pointer, dimension(:) :: volc_start => null(), &
                                       volc_end   => null(), &
                                       volc_i     => null(), &
                                       volc_j     => null()

     real, pointer, dimension(:,:,:) :: h2o2          ! buffer for H2O2 that is being replenished every 3 hours
                                                      ! if it is from climatology
  end type AChem_State

! Hook for the ESMF
! -----------------
  type AChem_Wrap
     type (AChem_State), pointer :: PTR => null()
  end type AChem_Wrap

contains


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 910.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: SetServices --- Sets IRF services for the AChem Grid Component
!
! !INTERFACE:

   subroutine SetServices(GC, RC)

! !ARGUMENTS:

    type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
    integer, optional                  :: RC  ! return code

! !DESCRIPTION: Sets Initialize, Run and Finalize services.
!
! !REVISION HISTORY:
!
!  1Dec2009  da Silva  First crack.
!
!EOP
!-------------------------------------------------------------------------

                            __Iam__('SetServices')

!   Local derived type aliases
!   --------------------------
    type (AChem_State), pointer  :: self   ! internal state
    type (AChem_Wrap)            :: wrap

    character(len=ESMF_MAXSTR) :: comp_name

!   Local variables
!   --------------------------
    integer :: n

!                              ------------

!   Get my name and set-up traceback handle
!   ---------------------------------------
    call ESMF_GridCompGet( GC, name=comp_name, __RC__ )
    Iam = TRIM(comp_name) // '::' // trim(Iam)

!   Wrap internal state for storing in GC; rename legacyState
!   -------------------------------------
    allocate(self, __STAT__)
    wrap%ptr => self

!   Load private Config Attributes
!   ------------------------------
    self%CF = ESMF_ConfigCreate(__RC__)
    call ESMF_ConfigLoadFile(self%CF, 'GEOS_AChemGridComp.rc', __RC__)

    call ESMF_ConfigGetAttribute(self%CF, self%verbose, label='verbose:', default=.false., __RC__)

    ! gas phase options
    call ESMF_ConfigGetAttribute(self%CF, self%gas_phase_chem, label='gas_chemistry:', default=.true., __RC__)

    ! aqueous phase options
    call ESMF_ConfigGetAttribute(self%CF, self%aqu_phase_chem, label='aqueous_chemistry:', default=.true., __RC__)
    call ESMF_ConfigGetAttribute(self%CF, self%aqu_solver_max_dt, label='aqueous_chemistry_solver_max_dt:', default=60.0, __RC__)
#if (0)
    ! combo(gas- and aqueous-phase) options
    call ESMF_ConfigGetAttribute(self%CF, self%combo_chem, label='combo_chemistry:', default=.true., __RC__)
    call ESMF_ConfigGetAttribute(self%CF, self%combo_solver_max_dt, label='combo_chemistry_max_dt:', default=10.0, __RC__)
#endif
    ! other options
    call ESMF_ConfigGetAttribute(self%CF, self%apply_diurnal_cycle, label='apply_diurnal_cycle:', default=.true., __RC__)

    ! volcanic emissions
    call ESMF_ConfigGetAttribute(self%CF, self%volcanic_emiss_file, Label='volcanoes:', default='/dev/null', __RC__)

    ! heights of aviation layers
    self%aviation_layers = 0.0
    call ESMF_ConfigFindLabel(self%CF, Label='aviation_vertical_layers:', __RC__)
    AVIATION_LAYERS: do n = 1, 4
        call ESMF_ConfigGetAttribute(self%CF, self%aviation_layers(n), __RC__)
    end do AVIATION_LAYERS

    ! OCS chemistry
    call ESMF_ConfigGetAttribute(self%CF, self%ocs_chem, Label='ocs_chemistry:', default=.false., __RC__)

    if (self%ocs_chem) then
        call ESMF_ConfigGetAttribute(self%CF, self%ocs_surface_vmr,  Label='ocs_surface_vmr:', __RC__)
    else
        self%ocs_surface_vmr = 0.0
    end if

    ! VOC chemistry
    call ESMF_ConfigGetAttribute(self%CF, self%voc_chem, Label='voc_chemistry:', default=.false., __RC__)

    if (self%voc_chem) then
        call ESMF_ConfigGetAttribute(self%CF, self%voc_BiomassBurnFactor, Label='voc_BiomassBurnFactor:', __RC__)
        call ESMF_ConfigGetAttribute(self%CF, self%voc_AnthroFactor, Label='voc_AnthroFactor:', __RC__)
        call ESMF_ConfigGetAttribute(self%CF, self%voc_MW, Label='voc_MW:', __RC__)
        call ESMF_ConfigGetAttribute(self%CF, self%soa_MW, Label='soa_MW:', __RC__)
    end if

    ! Minimalistic atmospheric mechanism for MAM and alike
    if (self%gas_phase_chem .or. self%aqu_phase_chem) then
        self%mam_chem = .true.
    else
        self%mam_chem = .false.
    end if

    if (MAPL_AM_I_ROOT()) then
        print *, trim(Iam)//': Configuration'
        print *, trim(Iam)//':     gas chemistry = ', self%gas_phase_chem
        print *, trim(Iam)//': aqueous chemistry = ', self%aqu_phase_chem
        print *, trim(Iam)//':     VOC chemistry = ', self%voc_chem
        print *, trim(Iam)//':     OCS chemistry = ', self%ocs_chem
        print *, ''
    end if


!                       ------------------------
!                       ESMF Functional Services
!                       ------------------------

!   Set the Initialize, Run, Finalize entry points
!   ----------------------------------------------
    call MAPL_GridCompSetEntryPoint(GC, ESMF_METHOD_INITIALIZE, Initialize_, __RC__)
    call MAPL_GridCompSetEntryPoint(GC, ESMF_METHOD_RUN,        Run_,        __RC__)
    call MAPL_GridCompSetEntryPoint(GC, ESMF_METHOD_FINALIZE,   Finalize_,   __RC__)

!   Store internal state in GC
!   --------------------------
    call ESMF_UserCompSetInternalState(GC, 'AChem_State', wrap, STATUS)
    VERIFY_(STATUS)

!                         ------------------
!                         MAPL Data Services
!                         ------------------

!BOS
!
! !IMPORT STATE:

    call MAPL_AddImportSpec(GC, &
                            SHORT_NAME = 'AIRDENS',  &
                            LONG_NAME  = 'Air density',  &
                            UNITS      = 'kg m-3', &
                            DIMS       = MAPL_DimsHorzVert,    &
                            VLOCATION  = MAPL_VLocationCenter,    &
                            __RC__)

    call MAPL_AddImportSpec(GC, &
                            SHORT_NAME = 'DELP',  &
                            LONG_NAME  = 'Pressure Thickness',  &
                            UNITS      = 'Pa', &
                            DIMS       = MAPL_DimsHorzVert,    &
                            VLOCATION  = MAPL_VLocationCenter,    &
                            __RC__)

    call MAPL_AddImportSpec(GC, &
                            SHORT_NAME = 'PLE',  &
                            LONG_NAME  = 'Edge pressure',  &
                            UNITS      = 'Pa', &
                            DIMS       = MAPL_DimsHorzVert,    &
                            VLOCATION  = MAPL_VLocationEdge,    &
                            __RC__)


    call MAPL_AddImportSpec(GC, &
                            SHORT_NAME = 'T',  &
                            LONG_NAME  = 'Air Temperature (from Dynamics)',  &
                            UNITS      = 'K', &
                            DIMS       = MAPL_DimsHorzVert,    &
                            VLOCATION  = MAPL_VLocationCenter,    &
                            __RC__)



    OPTIONAL_CHEM_IMPORT: if (self%mam_chem) then

    call MAPL_AddImportSpec(GC, &
                            SHORT_NAME = 'AREA',  &
                            LONG_NAME  = 'Cell area',  &
                            UNITS      = 'm2', &
                            DIMS       = MAPL_DimsHorzOnly,    &
                            VLOCATION  = MAPL_VLocationNone,    &
                            __RC__)

    call MAPL_AddImportSpec(GC, &
                            SHORT_NAME = 'ZLE',  &
                            LONG_NAME  = 'Edge heights',  &
                            UNITS      = 'm', &
                            DIMS       = MAPL_DimsHorzVert,    &
                            VLOCATION  = MAPL_VLocationEdge,    &
                            __RC__)

    call MAPL_AddImportSpec(GC, &
                            SHORT_NAME = 'QLTOT',  &
                            LONG_NAME  = 'Mass fraction of cloud liquid water',  &
                            UNITS      = 'kg kg-1', &
                            DIMS       = MAPL_DimsHorzVert,    &
                            VLOCATION  = MAPL_VLocationCenter,    &
                            __RC__)

    call MAPL_AddImportSpec(GC, &
                            SHORT_NAME = 'FCLD',  &
                            LONG_NAME  = 'Cloud fraction for radiation',  &
                            UNITS      = '1', &
                            DIMS       = MAPL_DimsHorzVert,    &
                            VLOCATION  = MAPL_VLocationCenter,    &
                            __RC__)

    call MAPL_AddImportSpec(GC, &
                            SHORT_NAME = 'U10N',  &
                            LONG_NAME  = 'Equivalent neutral 10 meter eastward wind',  &
                            UNITS      = 'm s-1', &
                            DIMS       = MAPL_DimsHorzOnly,    &
                            VLOCATION  = MAPL_VLocationNone,    &
                            __RC__)

    call MAPL_AddImportSpec(GC, &
                            SHORT_NAME = 'V10N',  &
                            LONG_NAME  = 'Equivalent neutral 10 meter northward wind',  &
                            UNITS      = 'm s-1', &
                            DIMS       = MAPL_DimsHorzOnly,    &
                            VLOCATION  = MAPL_VLocationNone,    &
                            __RC__)

    call MAPL_AddImportSpec(GC, &
                            SHORT_NAME = 'FROCEAN',  &
                            LONG_NAME  = 'Fraction of ocean',  &
                            UNITS      = '1', &
                            DIMS       = MAPL_DimsHorzOnly,    &
                            VLOCATION  = MAPL_VLocationNone,    &
                            __RC__)

    call MAPL_AddImportSpec(GC, &
                            SHORT_NAME = 'LWI',  &
                            LONG_NAME  = 'Land-water-ice flags',  &
                            UNITS      = '1', &
                            DIMS       = MAPL_DimsHorzOnly,    &
                            VLOCATION  = MAPL_VLocationNone,    &
                            __RC__)

    call MAPL_AddImportSpec(GC, &
                            SHORT_NAME = 'USTAR',  &
                            LONG_NAME  = 'Surface (friction) velocity scale',  &
                            UNITS      = 'm s-1', &
                            DIMS       = MAPL_DimsHorzOnly,    &
                            VLOCATION  = MAPL_VLocationNone,    &
                            __RC__)

    call MAPL_AddImportSpec(GC, &
                            SHORT_NAME = 'SH',  &
                            LONG_NAME  = 'Sensible heat flux',  &
                            UNITS      = 'W/m2', &
                            DIMS       = MAPL_DimsHorzOnly,    &
                            VLOCATION  = MAPL_VLocationNone,    &
                            __RC__)

    call MAPL_AddImportSpec(GC, &
                            SHORT_NAME = 'Z0H',  &
                            LONG_NAME  = 'Surface roughness for heat',  &
                            UNITS      = 'm', &
                            DIMS       = MAPL_DimsHorzOnly,    &
                            VLOCATION  = MAPL_VLocationNone,    &
                            __RC__)

    call MAPL_AddImportSpec(GC, &
                            SHORT_NAME = 'ZPBL',  &
                            LONG_NAME  = 'Height of PBL',  &
                            UNITS      = 'm', &
                            DIMS       = MAPL_DimsHorzOnly,    &
                            VLOCATION  = MAPL_VLocationNone,    &
                            __RC__)

    call MAPL_AddImportSpec(GC, &
                            SHORT_NAME = 'TS',  &
                            LONG_NAME  = 'Surface skin temperature',  &
                            UNITS      = 'K', &
                            DIMS       = MAPL_DimsHorzOnly,    &
                            VLOCATION  = MAPL_VLocationNone,    &
                            __RC__)

    call MAPL_AddImportSpec(GC, &
                            SHORT_NAME = 'H2O2',  &
                            LONG_NAME  = 'Hydrogen peroxide (H2O2)',  &
                            UNITS      = 'mol mol-1', &
                            DIMS       = MAPL_DimsHorzVert,    &
                            VLOCATION  = MAPL_VLocationCenter,    &
                            __RC__)

    call MAPL_AddImportSpec(GC, &
                            SHORT_NAME = 'OH',  &
                            LONG_NAME  = 'Hydroxyl radical (OH)',  &
                            UNITS      = 'mol mol-1', &
                            DIMS       = MAPL_DimsHorzVert,    &
                            VLOCATION  = MAPL_VLocationCenter,    &
                            __RC__)

    call MAPL_AddImportSpec(GC, &
                            SHORT_NAME = 'NO3',  &
                            LONG_NAME  = 'Nitrogen trixide (NO3)',  &
                            UNITS      = 'mol mol-1', &
                            DIMS       = MAPL_DimsHorzVert,    &
                            VLOCATION  = MAPL_VLocationCenter,    &
                            __RC__)

    call MAPL_AddImportSpec(GC, &
                            SHORT_NAME = 'O3',  &
                            LONG_NAME  = 'Ozone (mass mixing ratio)',  &
                            UNITS      = 'kg kg-1', &
                            DIMS       = MAPL_DimsHorzVert,    &
                            VLOCATION  = MAPL_VLocationCenter,    &
                            __RC__)

    call MAPL_AddImportSpec(GC, &
                            SHORT_NAME = 'DMS_CONC_OCEAN',  &
                            LONG_NAME  = 'Surface seawater concentration of DMS',  &
                            UNITS      = 'nmol L-1-1', &
                            DIMS       = MAPL_DimsHorzOnly,    &
                            VLOCATION  = MAPL_VLocationNone,    &
                            __RC__)

    call MAPL_AddImportSpec(GC, &
                            SHORT_NAME = 'SO2_EMIS_FIRES',  &
                            LONG_NAME  = 'SO2 emissions from biomass burning',  &
                            UNITS      = 'kg m-2 s-1', &
                            DIMS       = MAPL_DimsHorzOnly,    &
                            VLOCATION  = MAPL_VLocationNone,    &
                            __RC__)

    call MAPL_AddImportSpec(GC, &
                            SHORT_NAME = 'SO2_EMIS_NONENERGY',  &
                            LONG_NAME  = 'SO2 emissions from non-energy sectors',  &
                            UNITS      = 'kg m-2 s-1', &
                            DIMS       = MAPL_DimsHorzOnly,    &
                            VLOCATION  = MAPL_VLocationNone,    &
                            __RC__)

    call MAPL_AddImportSpec(GC, &
                            SHORT_NAME = 'SO2_EMIS_ENERGY',  &
                            LONG_NAME  = 'SO2 emissions from energy sector',  &
                            UNITS      = 'kg m-2 s-1', &
                            DIMS       = MAPL_DimsHorzOnly,    &
                            VLOCATION  = MAPL_VLocationNone,    &
                            __RC__)

    call MAPL_AddImportSpec(GC, &
                            SHORT_NAME = 'SO2_EMIS_SHIPPING',  &
                            LONG_NAME  = 'SO2 emissions from shipping sector',  &
                            UNITS      = 'kg m-2 s-1', &
                            DIMS       = MAPL_DimsHorzOnly,    &
                            VLOCATION  = MAPL_VLocationNone,    &
                            __RC__)

    call MAPL_AddImportSpec(GC, &
                            SHORT_NAME = 'SO2_EMIS_AIRCRAFT_LTO',  &
                            LONG_NAME  = 'SO2 emissions from aviation (LTO layer)',  &
                            UNITS      = 'kg m-2 s-1', &
                            DIMS       = MAPL_DimsHorzOnly,    &
                            VLOCATION  = MAPL_VLocationNone,    &
                            __RC__)

    call MAPL_AddImportSpec(GC, &
                            SHORT_NAME = 'SO2_EMIS_AIRCRAFT_CDS',  &
                            LONG_NAME  = 'SO2 emissions from aviation (CDS layer)',  &
                            UNITS      = 'kg m-2 s-1', &
                            DIMS       = MAPL_DimsHorzOnly,    &
                            VLOCATION  = MAPL_VLocationNone,    &
                            __RC__)

    call MAPL_AddImportSpec(GC, &
                            SHORT_NAME = 'SO2_EMIS_AIRCRAFT_CRS',  &
                            LONG_NAME  = 'SO2 emissions from aviation (CRS layer)',  &
                            UNITS      = 'kg m-2 s-1', &
                            DIMS       = MAPL_DimsHorzOnly,    &
                            VLOCATION  = MAPL_VLocationNone,    &
                            __RC__)

    call MAPL_AddImportSpec(GC, &
                            SHORT_NAME = 'NH3_EMIS',  &
                            LONG_NAME  = 'NH3 emissions - all sectors excluding biomass burning',  &
                            UNITS      = 'kg m-2 s-1', &
                            DIMS       = MAPL_DimsHorzOnly,    &
                            VLOCATION  = MAPL_VLocationNone,    &
                            __RC__)

    call MAPL_AddImportSpec(GC, &
                            SHORT_NAME = 'NH3_EMIS_FIRE',  &
                            LONG_NAME  = 'NH3 emissions - biomass burning',  &
                            UNITS      = 'kg m-2 s-1', &
                            DIMS       = MAPL_DimsHorzOnly,    &
                            VLOCATION  = MAPL_VLocationNone,    &
                            __RC__)

    call MAPL_AddImportSpec(GC, &
                            SHORT_NAME = 'SOAG_EMIS',  &
                            LONG_NAME  = 'SOA(gas) surface emissions',  &
                            UNITS      = 'm-2 s-1', &
                            DIMS       = MAPL_DimsHorzOnly,    &
                            VLOCATION  = MAPL_VLocationNone,    &
                            __RC__)

    end if OPTIONAL_CHEM_IMPORT



    OPTIONAL_VOC_IMPORTS: if (self%voc_chem) then

     if (.not. self%mam_chem) then
        call MAPL_AddImportSpec(GC, &
                                SHORT_NAME = 'OH',  &
                                LONG_NAME  = 'Hydroxyl radical (OH)',  &
                                UNITS      = 'mol mol-1', &
                                DIMS       = MAPL_DimsHorzVert,    &
                                VLOCATION  = MAPL_VLocationCenter,    &
                                __RC__)
     end if

        call MAPL_AddImportSpec(GC, &
                                SHORT_NAME = 'CO_BIOMASS_VOC',  &
                                LONG_NAME  = 'CO Biomass Burning Emissions',  &
                                UNITS      = 'kg m-2 s-1', &
                                DIMS       = MAPL_DimsHorzOnly,    &
                                VLOCATION  = MAPL_VLocationNone,   &
                                __RC__)

         call MAPL_AddImportSpec(GC, &
                                SHORT_NAME = 'CO_BF_VOC',  &
                                LONG_NAME  = 'CO Biofuel Emissions',  &
                                UNITS      = 'kg m-2 s-1', &
                                DIMS       = MAPL_DimsHorzOnly,    &
                                VLOCATION  = MAPL_VLocationNone,    &
                                __RC__)

         call MAPL_AddImportSpec(GC, &
                                SHORT_NAME = 'CO_FS_VOC',  &
                                LONG_NAME  = 'CO Fossil Fuel Emissions',  &
                                UNITS      = 'kg m-2 s-1', &
                                DIMS       = MAPL_DimsHorzOnly,    &
                                VLOCATION  = MAPL_VLocationNone,    &
                                __RC__)

    end if OPTIONAL_VOC_IMPORTS


    OPTIONAL_OCS_IMPORTS: if (self%ocs_chem) then
        call MAPL_AddImportSpec(GC, &
                                SHORT_NAME = 'TROPP',  &
                                LONG_NAME  = 'Tropopause pressure based on blended estimate',  &
                                UNITS      = 'Pa', &
                                DIMS       = MAPL_DimsHorzOnly,    &
                                VLOCATION  = MAPL_VLocationNone,    &
                                __RC__)

        call MAPL_AddImportSpec(GC,                                &
                                SHORT_NAME = 'OHSTRAT',            &
                                LONG_NAME  = 'Hydroxyl radical',   &
                                UNITS      = 'mol mol-1',          &
                                DIMS       = MAPL_DimsHorzVert,    &
                                VLOCATION  = MAPL_VLocationCenter, &
                                __RC__)

        call MAPL_AddImportSpec(GC,                                &
                                SHORT_NAME = 'O3P',                &
                                LONG_NAME  = 'O triplet P',        &
                                UNITS      = 'mol mol-1',          &
                                DIMS       = MAPL_DimsHorzVert,    &
                                VLOCATION  = MAPL_VLocationCenter, &
                                __RC__)

        call MAPL_AddImportSpec(GC,                                  &
                                SHORT_NAME = 'OCS_JRATE',            &
                                LONG_NAME  = 'OCS photolysis rates', &
                                UNITS      = 's-1',                  &
                                DIMS       = MAPL_DimsHorzVert,      &
                                VLOCATION  = MAPL_VLocationCenter,   &
                                __RC__)

    end if OPTIONAL_OCS_IMPORTS


! !INTERNAL STATE:

    OPTIONAL_CHEM_INTERNAL: if (self%mam_chem) then

        call MAPL_AddInternalSpec(GC,                                                  &
                                  SHORT_NAME = trim(comp_name)//'::'//'DMS',           &
                                  LONG_NAME  = 'Dimethyl sulfide (DMS)',               &
                                  UNITS      = 'mol mol-1',                            &
                                  DIMS       = MAPL_DimsHorzVert,                      &
                                  VLOCATION  = MAPL_VLocationCenter,                   &
                                  FRIENDLYTO = 'DYNAMICS:TURBULENCE:MOIST',            &
                                  ADD2EXPORT = .true.,                                 &
                                  __RC__)

        call MAPL_AddInternalSpec(GC,                                                  &
                                  SHORT_NAME = trim(comp_name)//'::'//'MSA',           &
                                  LONG_NAME  = 'Methanesulfonic acid (MSA)',           &
                                  UNITS      = 'mol mol-1',                            &
                                  DIMS       = MAPL_DimsHorzVert,                      &
                                  VLOCATION  = MAPL_VLocationCenter,                   &
                                  FRIENDLYTO = 'DYNAMICS:TURBULENCE:MOIST',            &
                                  ADD2EXPORT = .true.,                                 &
                                  __RC__)

        call MAPL_AddInternalSpec(GC,                                                  &
                                  SHORT_NAME = trim(comp_name)//'::'//'SO2',           &
                                  LONG_NAME  = 'Sulfur dioxide (SO2)',                 &
                                  UNITS      = 'mol mol-1',                            &
                                  DIMS       = MAPL_DimsHorzVert,                      &
                                  VLOCATION  = MAPL_VLocationCenter,                   &
                                  FRIENDLYTO = 'DYNAMICS:TURBULENCE:MOIST:MAM',        &
                                  ADD2EXPORT = .true.,                                 &
                                  __RC__)

        call MAPL_AddInternalSpec(GC,                                                  &
                                  SHORT_NAME = trim(comp_name)//'::'//'H2SO4',         &
                                  LONG_NAME  = 'Sulfuric acid (H2SO4 gas)',            &
                                  UNITS      = 'mol mol-1',                            &
                                  DIMS       = MAPL_DimsHorzVert,                      &
                                  VLOCATION  = MAPL_VLocationCenter,                   &
                                  FRIENDLYTO = 'DYNAMICS:TURBULENCE:MOIST:MAM',        &
                                  ADD2EXPORT = .true.,                                 &
                                  __RC__)

        call MAPL_AddInternalSpec(GC,                                                  &
                                  SHORT_NAME = trim(comp_name)//'::'//'NH3',           &
                                  LONG_NAME  = 'Ammonia (NH3)',                        &
                                  UNITS      = 'mol mol-1',                            &
                                  DIMS       = MAPL_DimsHorzVert,                      &
                                  VLOCATION  = MAPL_VLocationCenter,                   &
                                  FRIENDLYTO = 'DYNAMICS:TURBULENCE:MOIST:MAM',        &
                                  ADD2EXPORT = .true.,                                 &
                                  __RC__)

        call MAPL_AddInternalSpec(GC,                                                  &
                                  SHORT_NAME = trim(comp_name)//'::'//'SOAG',          &
                                  LONG_NAME  = 'Secondary Organic Aerosols (SOA gas)', &
                                  UNITS      = 'mol mol-1',                            &
                                  DIMS       = MAPL_DimsHorzVert,                      &
                                  VLOCATION  = MAPL_VLocationCenter,                   &
                                  FRIENDLYTO = 'DYNAMICS:TURBULENCE:MOIST:MAM',        &
                                  ADD2EXPORT = .true.,                                 &
                                  __RC__)
    end if OPTIONAL_CHEM_INTERNAL


    OPTIONAL_OCS_INTERNAL: if (self%ocs_chem) then

        call MAPL_AddInternalSpec(GC,                                                  &
                                  SHORT_NAME = trim(comp_name)//'::'//'OCS',           &
                                  LONG_NAME  = 'Carbonyl Sulfide (OCS gas)',           &
                                  UNITS      = 'mol mol-1',                            &
                                  DIMS       = MAPL_DimsHorzVert,                      &
                                  VLOCATION  = MAPL_VLocationCenter,                   &
                                  FRIENDLYTO = 'DYNAMICS:TURBULENCE:MOIST:MAM',        &
                                  ADD2EXPORT = .true.,                                 &
                                  __RC__)

    end if OPTIONAL_OCS_INTERNAL


    OPTIONAL_VOC_INTERNAL: if (self%voc_chem) then

        call MAPL_AddInternalSpec(GC,                                                  &
                                  SHORT_NAME = trim(comp_name)//'::'//'VOC',              &
                                  LONG_NAME  = 'Volatile Organic Compound (VOC) --anthropogenic sources', &
                                  UNITS      = 'mol mol-1',                            &
                                  DIMS       = MAPL_DimsHorzVert,                      &
                                  VLOCATION  = MAPL_VLocationCenter,                   &
                                  FRIENDLYTO = 'DYNAMICS:TURBULENCE:MOIST',            &
                                  ADD2EXPORT = .true.,                                 &
                                  __RC__)

        call MAPL_AddInternalSpec(GC,                                                  &
                                  SHORT_NAME = trim(comp_name)//'::'//'VOCbiob',       &
                                  LONG_NAME  = 'Volatile Organic Compound (VOC) -- biomass burning sources',&
                                  UNITS      = 'mol mol-1',                            &
                                  DIMS       = MAPL_DimsHorzVert,                      &
                                  VLOCATION  = MAPL_VLocationCenter,                   &
                                  FRIENDLYTO = 'DYNAMICS:TURBULENCE:MOIST',            &
                                  ADD2EXPORT = .true.,                                 &
                                  __RC__)

    end if OPTIONAL_VOC_INTERNAL

! !EXTERNAL STATE:
    OPTIONAL_CHEM_EXPORT: if (self%mam_chem) then
#include "GEOS_AChem_ExportSpec___.h"
    end if OPTIONAL_CHEM_EXPORT


    OPTIONAL_VOC_EXPORT: if (self%voc_chem) then

     if (.not. self%mam_chem) then
     call MAPL_AddExportSpec(GC,  &
        SHORT_NAME         = 'OH',  &
        LONG_NAME          = 'OH with imposed diurnal cycle',  &
        UNITS              = 'mol mol-1', &
        DIMS               = MAPL_DimsHorzVert,    &
        VLOCATION          = MAPL_VLocationCenter,    &
                                                       RC=STATUS  )
     VERIFY_(STATUS)
     end if


     call MAPL_AddExportSpec(GC,  &
        SHORT_NAME         = 'pSOA_ANTHRO_VOC',  &
        LONG_NAME          = 'Production of SOA from Anthropogenic + Biofuel Burning VOC',  &
        UNITS              = 'kg m-3 s-1', &
        DIMS               = MAPL_DimsHorzVert,    &
        VLOCATION          = MAPL_VLocationCenter,    &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,  &
        SHORT_NAME         = 'pSOA_ANTHRO_VOC_MMRday',  &
        LONG_NAME          = 'Production of SOA from Anthropogenic + Biofuel Burning VOC',  &
        UNITS              = 'kg m-3 d-1', &
        DIMS               = MAPL_DimsHorzVert,    &
        VLOCATION          = MAPL_VLocationCenter,    &
                                                       RC=STATUS  )
     VERIFY_(STATUS)


     call MAPL_AddExportSpec(GC,  &
        SHORT_NAME         = 'pSOA_BIOB_VOC',  &
        LONG_NAME          = 'Production of SOA from Biomass Burning VOC',  &
        UNITS              = 'kg m-3 s-1', &
        DIMS               = MAPL_DimsHorzVert,    &
        VLOCATION          = MAPL_VLocationCenter,    &
                                                       RC=STATUS  )
     VERIFY_(STATUS)


     call MAPL_AddExportSpec(GC,  &
        SHORT_NAME         = 'pSOA_BIOB_VOC_MMRday',  &
        LONG_NAME          = 'Production of SOA from Biomass Burning VOC',  &
        UNITS              = 'kg m-3 d-1', &
        DIMS               = MAPL_DimsHorzVert,    &
        VLOCATION          = MAPL_VLocationCenter,    &
                                                       RC=STATUS  )
     VERIFY_(STATUS)


    end if OPTIONAL_VOC_EXPORT


    OPTIONAL_OCS_EXPORT: if (self%ocs_chem) then
     call MAPL_AddExportSpec(GC,  &
        SHORT_NAME         = 'pSO2_OCS',  &
        LONG_NAME          = 'Production of SO2 from OCS',  &
        UNITS              = 'kg kg-1 s-1', &
        DIMS               = MAPL_DimsHorzVert,    &
        VLOCATION          = MAPL_VLocationCenter,    &
                                                       RC=STATUS  )
     VERIFY_(STATUS)


     call MAPL_AddExportSpec(GC,  &
        SHORT_NAME         = 'pSO2_OCS_OH',  &
        LONG_NAME          = 'Production of SO2 from OCS+OH',  &
        UNITS              = 'kg kg-1 s-1', &
        DIMS               = MAPL_DimsHorzVert,    &
        VLOCATION          = MAPL_VLocationCenter,    &
                                                       RC=STATUS  )
     VERIFY_(STATUS)


     call MAPL_AddExportSpec(GC,  &
        SHORT_NAME         = 'pSO2_OCS_O3p',  &
        LONG_NAME          = 'Production of SO2 from OCS+O3p',  &
        UNITS              = 'kg kg-1 s-1', &
        DIMS               = MAPL_DimsHorzVert,    &
        VLOCATION          = MAPL_VLocationCenter,    &
                                                       RC=STATUS  )
     VERIFY_(STATUS)


     call MAPL_AddExportSpec(GC,  &
        SHORT_NAME         = 'pSO2_OCS_jOCS',  &
        LONG_NAME          = 'Production of SO2 from OCS photolysis',  &
        UNITS              = 'kg kg-1 s-1', &
        DIMS               = MAPL_DimsHorzVert,    &
        VLOCATION          = MAPL_VLocationCenter,    &
                                                       RC=STATUS  )
     VERIFY_(STATUS)


     call MAPL_AddExportSpec(GC,  &
        SHORT_NAME         = 'lOCS',  &
        LONG_NAME          = 'Loss rate of OCS (molec cm-3 s-1)',  &
        UNITS              = 'cm-3 s-1', &
        DIMS               = MAPL_DimsHorzVert,    &
        VLOCATION          = MAPL_VLocationCenter,    &
                                                       RC=STATUS  )
     VERIFY_(STATUS)


     call MAPL_AddExportSpec(GC,  &
        SHORT_NAME         = 'lOCS_OH',  &
        LONG_NAME          = 'Loss rate of OCS from OCS+OH(molec cm-3 s-1)',  &
        UNITS              = 'cm-3 s-1', &
        DIMS               = MAPL_DimsHorzVert,    &
        VLOCATION          = MAPL_VLocationCenter,    &
                                                       RC=STATUS  )
     VERIFY_(STATUS)


     call MAPL_AddExportSpec(GC,  &
        SHORT_NAME         = 'lOCS_O3p',  &
        LONG_NAME          = 'Loss rate of OCS from OCS+O3p(molec cm-3 s-1)',  &
        UNITS              = 'cm-3 s-1', &
        DIMS               = MAPL_DimsHorzVert,    &
        VLOCATION          = MAPL_VLocationCenter,    &
                                                       RC=STATUS  )
     VERIFY_(STATUS)


     call MAPL_AddExportSpec(GC,  &
        SHORT_NAME         = 'lOCS_jOCS',  &
        LONG_NAME          = 'Loss rate of OCS from photolysis (molec cm-3 s-1)',  &
        UNITS              = 'cm-3 s-1', &
        DIMS               = MAPL_DimsHorzVert,    &
        VLOCATION          = MAPL_VLocationCenter,    &
                                                       RC=STATUS  )
     VERIFY_(STATUS)


     call MAPL_AddExportSpec(GC,  &
        SHORT_NAME         = 'pScl_OCS',  &
        LONG_NAME          = 'Production of SO2 from OCS (column integrated)',  &
        UNITS              = 'kg m-2', &
        DIMS               = MAPL_DimsHorzOnly,    &
        VLOCATION          = MAPL_VLocationNone,    &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

    end if OPTIONAL_OCS_EXPORT

!EOS

!   Set the Profiling timers
!   ------------------------
    call MAPL_TimerAdd(GC, name = 'TOTAL',                __RC__)
    call MAPL_TimerAdd(GC, name = 'RUN',                  __RC__)
    call MAPL_TimerAdd(GC, name = '-EMISSIONS',           __RC__)
    call MAPL_TimerAdd(GC, name = '-CHEMISTRY',           __RC__)
    call MAPL_TimerAdd(GC, name = '--CHEMISTRY_GAS',      __RC__)
    call MAPL_TimerAdd(GC, name = '--CHEMISTRY_AQUEOUS',  __RC__)
    call MAPL_TimerAdd(GC, name = '--CHEMISTRY_VOC',      __RC__)
    call MAPL_TimerAdd(GC, name = '--CHEMISTRY_OCS',      __RC__)
    call MAPL_TimerAdd(GC, name = 'INITIALIZE',           __RC__)


!   Generic Set Services
!   --------------------
    call MAPL_GenericSetServices(GC, __RC__)

!   All done
!   --------

    RETURN_(ESMF_SUCCESS)

  end subroutine SetServices


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Initialize_ --- Initialize AChem
!
! !INTERFACE:
!

   subroutine Initialize_(GC, IMPORT, EXPORT, CLOCK, rc)

! !USES:

   implicit none

! !INPUT PARAMETERS:

   type(ESMF_Clock),  intent(inout)   :: CLOCK   ! The clock

! !OUTPUT PARAMETERS:

   type(ESMF_GridComp), intent(inout) :: GC      ! Grid Component
   type(MAPL_MetaComp), pointer       :: mgState ! MAPL generic state
   type(ESMF_State), intent(inout)    :: IMPORT  ! Import State
   type(ESMF_State), intent(inout)    :: EXPORT  ! Export State
   integer, intent(out)               :: rc      ! Error return code:
                                                 !  0 - all is well
                                                 !  1 -

! !DESCRIPTION: This is a simple ESMF wrapper.
!
! !REVISION HISTORY:
!
!  01Dec2009 da Silva  First crack.
!
!EOP
!-------------------------------------------------------------------------

                              __Iam__('Initialize_')

    type(AChem_State), pointer      :: self        ! Legacy state
    type(ESMF_Grid)                 :: GRID        ! Grid
    type(ESMF_Config)               :: CF          ! Universal Config

    integer                         :: i1, i2, im  ! 3D Dimensions
    integer                         :: j1, j2, jm  !
    integer                         :: km          !

    integer                         :: nymd, nhms  ! date, time
    real                            :: cdt         ! time step in secs

    character(len=ESMF_MAXSTR)      :: comp_name   ! component's name

    real, pointer, dimension(:,:,:) :: q_H2O2      ! H2O2
    logical, parameter :: using_GMI_H2O2 = .false. ! coupling with GMI is not implemented


!  Declare pointers to IMPORT/EXPORT/INTERNAL states
!  -------------------------------------------------
#if(0)
#include "GEOS_AChem_DeclarePointer___.h"
#endif

!  Get my name and set-up traceback handle
!  ---------------------------------------
   call ESMF_GridCompGet(GC, name=comp_name, __RC__)
   Iam = trim(comp_name) // '::' // trim(Iam)

!                               --------
   if (MAPL_AM_I_ROOT()) then
      print *, trim(Iam)//': Starting...'
      print *, ''
   end if


!  Get my internal MAPL_Generic state
!  -----------------------------------
   call MAPL_GetObjectFromGC(GC, mgState, __RC__)

   call MAPL_TimerOn(mgState, 'TOTAL',      __RC__)
   call MAPL_TimerOn(mgState, 'INITIALIZE', __RC__)


!  Initialize MAPL Generic
!  -----------------------
   call MAPL_GenericInitialize(GC, IMPORT, EXPORT, clock, __RC__)

!  Get pointers to IMPORT/EXPORT/INTERNAL states
!  ---------------------------------------------
#if(0)
#include "GEOS_AChem_GetPointer___.h"
#endif

!  Extract relevant runtime information
!  ------------------------------------
   call extract_(GC, CLOCK, self, GRID, CF, i1, i2, im, j1, j2, jm, km, nymd, nhms, cdt, __RC__)


!  Set the grid
!  -------------------------------------------------
   self%grid = GRID

!  Initialize volcanic emissions timestamp
!  ---------------------------------------
   self%nymd_volcanic_emiss = -1

   nullify(self%volc_lat)
   nullify(self%volc_lon)
   nullify(self%volc_SO2)
   nullify(self%volc_elev)
   nullify(self%volc_cloud)
   nullify(self%volc_start)
   nullify(self%volc_end)
   nullify(self%volc_i)
   nullify(self%volc_j)

   INIT_H2O2: if (self%aqu_phase_chem) then
!  Initialize the internal copy of H2O2
!  ------------------------------------
   allocate(self%h2o2(i1:i2,j1:j2,km), __STAT__)

   if (using_GMI_H2O2) then
       self%h2o2 = 0.0 ! initial value is not important if H2O2 is from GMI
   else
       call MAPL_GetPointer(import, q_H2O2, 'H2O2', __RC__)
       self%h2o2 = q_H2O2
   end if
   end if INIT_H2O2

!  All done
!  --------
   call MAPL_TimerOff(mgState, 'INITIALIZE', __RC__)
   call MAPL_TimerOff(mgState, 'TOTAL',      __RC__)

   RETURN_(ESMF_SUCCESS)

   end subroutine Initialize_


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Run_ --- Runs AChem
!
! !INTERFACE:
!

   subroutine Run_(GC, IMPORT, EXPORT, CLOCK, rc)

! !USES:

   implicit none

! !INPUT PARAMETERS:

   type(ESMF_Clock),  intent(inout)   :: CLOCK  ! The clock

! !OUTPUT PARAMETERS:

   type(ESMF_GridComp), intent(inout) :: GC     ! Grid Component
   type(ESMF_State), intent(inout)    :: IMPORT ! Import State
   type(ESMF_State), intent(inout)    :: EXPORT ! Export State
   integer, intent(out)               :: rc     ! Error return code:
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

                              __Iam__('Run_')

   type(AChem_State), pointer     :: self        ! Legacy state
   type(ESMF_Grid)                :: GRID        ! Grid
   type(ESMF_Config)              :: CF          ! Universal Config

   type(MAPL_MetaComp), pointer   :: mgState     ! MAPL generic state
   type(ESMF_Alarm)               :: run_alarm
   logical                        :: run_alarm_ringing

   integer                        :: i1, i2, im  ! 3D Dimensions
   integer                        :: j1, j2, jm  !
   integer                        :: km          !

   real(ESMF_KIND_R4), pointer    :: lons(:,:)   ! Longitudes
   real(ESMF_KIND_R4), pointer    :: lats(:,:)   ! Latitudes

   integer                        :: nymd, nhms  ! date, time
   real                           :: cdt         ! time step in secs

   character(len=ESMF_MAXSTR)     :: comp_name

   integer                        :: k, k1

!  Input fields
!  ------------
   real, pointer, dimension(:,:)   :: cell_area => null()

   real, pointer, dimension(:,:,:) :: density_air => null()
   real, pointer, dimension(:,:,:) :: temperature => null()
   real, pointer, dimension(:,:,:) :: lwc => null()
   real, pointer, dimension(:,:,:) :: fcld => null()
   real, pointer, dimension(:,:,:) :: ple => null()
   real, pointer, dimension(:,:,:) :: delp => null()
   real, pointer, dimension(:,:,:) :: zle => null()
   real, pointer, dimension(:,:)   :: tropp => null()

   real, pointer, dimension(:,:)   :: u10n => null()
   real, pointer, dimension(:,:)   :: v10n => null()

   real, pointer, dimension(:,:)   :: tskin => null()

   real, pointer, dimension(:,:)   :: fr_ocean => null()

   real, pointer, dimension(:,:)   :: oro => null()
   real, pointer, dimension(:,:)   :: ustar => null()
   real, pointer, dimension(:,:)   :: shflux => null()
   real, pointer, dimension(:,:)   :: pblh => null()
   real, pointer, dimension(:,:)   :: z0h => null()

   real, pointer, dimension(:,:,:) :: q_OH => null()
   real, pointer, dimension(:,:,:) :: q_NO3 => null()
   real, pointer, dimension(:,:,:) :: q_H2O2 => null()
   real, pointer, dimension(:,:,:) :: q_O3 => null()

   real, pointer, dimension(:,:,:) :: q_OH_STRATCHEM => null()
   real, pointer, dimension(:,:,:) :: q_O3P_STRATCHEM => null()
   real, pointer, dimension(:,:,:) :: j_ocs => null()

   real, pointer, dimension(:,:)   :: DMS_ocean => null()

   real, pointer, dimension(:,:)   :: SO2_emiss_bb => null()
   real, pointer, dimension(:,:)   :: SO2_emiss_nonenergy => null()
   real, pointer, dimension(:,:)   :: SO2_emiss_energy => null()
   real, pointer, dimension(:,:)   :: SO2_emiss_shipping => null()
   real, pointer, dimension(:,:)   :: SO2_emiss_aviation_lto => null()
   real, pointer, dimension(:,:)   :: SO2_emiss_aviation_cds => null()
   real, pointer, dimension(:,:)   :: SO2_emiss_aviation_crs => null()

   real, pointer, dimension(:,:)   :: NH3_emiss => null()
   real, pointer, dimension(:,:)   :: NH3_emiss_bb => null()

   real, pointer, dimension(:,:)   :: SOAG_emiss => null()

   real,pointer,dimension(:,:)     :: co_biomass_voc => null()
   real,pointer,dimension(:,:)     :: co_bf_voc => null()
   real,pointer,dimension(:,:)     :: co_fs_voc => null()

!  Export fields
!  -------------
   type(ESMF_State)                :: internal

   real, pointer, dimension(:,:,:) :: ptr3d => null()
   real, pointer, dimension(:,:)   :: ptr2d => null()

   real, pointer, dimension(:,:,:) :: q_DMS => null()
   real, pointer, dimension(:,:,:) :: q_MSA => null()
   real, pointer, dimension(:,:,:) :: q_SO2 => null()
   real, pointer, dimension(:,:,:) :: q_H2SO4 => null()
   real, pointer, dimension(:,:,:) :: q_NH3 => null()
   real, pointer, dimension(:,:,:) :: q_SOAG => null()
   real, pointer, dimension(:,:,:) :: q_OCS => null()
   real, pointer, dimension(:,:,:) :: q_VOCanth => null()
   real, pointer, dimension(:,:,:) :: q_VOCbiob => null()

   real, allocatable, dimension(:,:,:) :: q_OAanth
   real, allocatable, dimension(:,:,:) :: q_OAanthmmrd
   real, allocatable, dimension(:,:,:) :: q_OAbiob
   real, allocatable, dimension(:,:,:) :: q_OAbiobmmrd

   real, allocatable, dimension(:,:,:) :: pSO2_OCS               ! production of S from OCS in the stratosphere
   real, allocatable, dimension(:,:,:) :: pSO2_OCS_OH            ! production of S from OCS+OH
   real, allocatable, dimension(:,:,:) :: pSO2_OCS_O3p           ! production of S from OCS+O3p
   real, allocatable, dimension(:,:,:) :: pSO2_OCS_jOCS          ! production of S from OCS photolysis
   real, allocatable, dimension(:,:,:) :: lOCS                   ! loss OCS, 'molecules cm-3 s-1'
   real, allocatable, dimension(:,:,:) :: lOCS_OH                ! loss rate of OCS from OCS+OH, 'molec cm-3 s-1'
   real, allocatable, dimension(:,:,:) :: lOCS_O3p               ! loss rate of OCS from OCS+O3p, 'molec cm-3 s-1'
   real, allocatable, dimension(:,:,:) :: lOCS_jOCS              ! loss rate of OCS from photolysis, 'molec cm-3 s-1'

   real, pointer, dimension(:,:)   :: dry_dep_DMS => null()  ! dry deposition fluxes
   real, pointer, dimension(:,:)   :: dry_dep_MSA => null()
   real, pointer, dimension(:,:)   :: dry_dep_SO2 => null()
   real, pointer, dimension(:,:)   :: dry_dep_H2SO4 => null()
   real, pointer, dimension(:,:)   :: dry_dep_NH3 => null()
   real, pointer, dimension(:,:)   :: dry_dep_SOAG => null()

   real, pointer, dimension(:,:,:) :: ddt_DMS_gas => null()  ! tendencies due to gas phase chemistry
   real, pointer, dimension(:,:,:) :: ddt_MSA_gas => null()
   real, pointer, dimension(:,:,:) :: ddt_SO2_gas => null()
   real, pointer, dimension(:,:,:) :: ddt_H2SO4_gas => null()
   real, pointer, dimension(:,:,:) :: ddt_NH3_gas => null()
   real, pointer, dimension(:,:,:) :: ddt_SOAG_gas => null()

   real, pointer, dimension(:,:,:) :: ddt_DMS_aq => null()   ! tendencies due to aqueous phase chemistry
   real, pointer, dimension(:,:,:) :: ddt_MSA_aq => null()
   real, pointer, dimension(:,:,:) :: ddt_SO2_aq => null()
   real, pointer, dimension(:,:,:) :: ddt_H2SO4_aq => null()
   real, pointer, dimension(:,:,:) :: ddt_NH3_aq => null()
   real, pointer, dimension(:,:,:) :: ddt_SOAG_aq => null()

   real, pointer, dimension(:,:,:) :: DMS_g_ => null()       ! tendencies due to gas phase chemistry
   real, pointer, dimension(:,:,:) :: MSA_g_ => null()
   real, pointer, dimension(:,:,:) :: SO2_g_ => null()
   real, pointer, dimension(:,:,:) :: H2SO4_g_ => null()
   real, pointer, dimension(:,:,:) :: NH3_g_ => null()
   real, pointer, dimension(:,:,:) :: SOAG_g_ => null()

   real, pointer, dimension(:,:,:) :: DMS_a_ => null()       ! tendencies due to aqueous phase chemistry
   real, pointer, dimension(:,:,:) :: MSA_a_ => null()
   real, pointer, dimension(:,:,:) :: SO2_a_ => null()
   real, pointer, dimension(:,:,:) :: H2SO4_a_ => null()
   real, pointer, dimension(:,:,:) :: NH3_a_ => null()
   real, pointer, dimension(:,:,:) :: SOAG_a_ => null()



!  Dry deposition frequency
!  ------------------------
   real, allocatable, dimension(:,:) :: dry_dep_frequency
   real, allocatable, dimension(:,:) :: dq


!  DMS flux
!  ---------
   real, allocatable, dimension(:,:) :: flux_DMS

!  Sulfur diagnostics
!  ------------------
   real, allocatable, dimension(:,:) :: SO2_emiss_total

!  VOC local arrays
!  ----------------
   real, allocatable, dimension(:,:,:) :: dVOC
   real, allocatable, dimension(:,:,:) :: dOAanth, dOAbiob
   real, allocatable, dimension(:,:,:) :: rk_OA_OH
   real, allocatable, dimension(:,:,:) :: fanth



!  Work buffers of oxidant fields
!  ------------------------------
   real, allocatable, dimension(:,:,:) :: q_OH_
   real, allocatable, dimension(:,:,:) :: q_NO3_


!  local
!  -----
   real, parameter :: ORO_OCEAN   = 0.0
   real, parameter :: ORO_LAND    = 1.0
   real, parameter :: ORO_SEA_ICE = 2.0

   integer :: doy                         ! day of year
   real    :: f_hour, x_hour              ! UTC hour

   logical, parameter :: using_GMI_H2O2 = .false.    ! coupling with GMI is not implemented
   logical, parameter :: using_GMI_OH   = .false.
   logical, parameter :: using_GMI_NO3  = .false.

   integer :: n, n_steps

   real, allocatable, dimension(:,:) :: SO2_emiss_volc_expl, SO2_emiss_volc_nonexpl

   real, allocatable, dimension(:,:) :: day_time, night_time     ! day time and night time durations, s
   real, allocatable, dimension(:,:) :: f_day_time, f_night_time ! day time and night time factors

   real, allocatable, dimension(:,:) :: sza
   real, allocatable, dimension(:,:) :: cos_sza
   real, allocatable, dimension(:,:) :: sum_cos_sza

   real, allocatable, dimension(:,:) :: cmd_S                    ! column mass density (i.e., column integrated mass loading) of S from gas species
   real, allocatable, dimension(:,:) :: cmd_DMS
   real, allocatable, dimension(:,:) :: cmd_MSA
   real, allocatable, dimension(:,:) :: cmd_SO2
   real, allocatable, dimension(:,:) :: cmd_H2SO4


   real, allocatable, dimension(:,:) :: cmd_NH3
   real, allocatable, dimension(:,:) :: cmd_N

   real, allocatable, dimension(:,:) :: cmd_SOAG

   real, allocatable, dimension(:,:,:) :: pSO4_aq                ! production rates from aquesous chemistry
   real, allocatable, dimension(:,:,:) :: pNH4_aq
   real, allocatable, dimension(:,:,:) :: pSO4_aq_SO2
   real, allocatable, dimension(:,:,:) :: pSO4_aq_H2SO4
   real, allocatable, dimension(:,:,:) :: pNH4_aq_NH3

   real, allocatable, dimension(:,:) :: cpl_NH3
   real, allocatable, dimension(:,:) :: cpl_DMS
   real, allocatable, dimension(:,:) :: cpl_MSA
   real, allocatable, dimension(:,:) :: cpl_SO2
   real, allocatable, dimension(:,:) :: cpl_H2SO4


!  Declare pointers to IMPORT/EXPORT/INTERNAL states
!  -------------------------------------------------
#if(0)
#include "GEOS_AChem_DeclarePointer___.h"
#endif

!  Get my name and set-up traceback handle
!  ---------------------------------------
   call ESMF_GridCompGet(GC, name=comp_name, __RC__)
   Iam = trim(comp_name) // '::' // trim(Iam)

!  Get pointers to IMPORT/EXPORT/INTERNAL states
!  ---------------------------------------------
#if(0)
   #include "GEOS_AChem_GetPointer___.h"
#else

#endif

!  Get my internal MAPL_Generic state
!  -----------------------------------
   call MAPL_GetObjectFromGC(GC, mgState, __RC__)

   call MAPL_TimerOn(mgState, 'TOTAL', __RC__)
   call MAPL_TimerOn(mgState, 'RUN',   __RC__)

!  Get parameters from generic state
!  ----------------------------------
   call MAPL_Get(mgState, LONS=lons, LATS=lats, RunAlarm=run_alarm, __RC__)


!  If it is time, update AChem state
!  ---------------------------------
   run_alarm_ringing = ESMF_AlarmIsRinging(run_alarm, __RC__)

   if (run_alarm_ringing) then
       call ESMF_AlarmRingerOff(run_alarm, __RC__)
   else
       RETURN_(ESMF_SUCCESS)
   endif


!  Extract relevant runtime information
!  ------------------------------------
   call extract_(GC, CLOCK, self, GRID, CF, i1, i2, im, j1, j2, jm, km, nymd, nhms, cdt, __RC__)


!  Get Imports
!  --------------
   call MAPL_GetPointer(import, density_air, 'AIRDENS', __RC__)
   call MAPL_GetPointer(import, delp,        'DELP',    __RC__)
   call MAPL_GetPointer(import, ple,         'PLE',     __RC__)
   call MAPL_GetPointer(import, temperature, 'T',       __RC__)


   if (self%mam_chem) then
       call MAPL_GetPointer(import, cell_area,   'AREA',    __RC__)
       call MAPL_GetPointer(import, zle,         'ZLE',     __RC__)
       call MAPL_GetPointer(import, lwc,         'QLTOT',   __RC__)
       call MAPL_GetPointer(import, fcld,        'FCLD',    __RC__)

       call MAPL_GetPointer(import, u10n,        'U10N',    __RC__)
       call MAPL_GetPointer(import, v10n,        'V10N',    __RC__)

       call MAPL_GetPointer(import, tskin,       'TS',      __RC__)

       call MAPL_GetPointer(import, fr_ocean,    'FROCEAN', __RC__)

       call MAPL_GetPointer(import, oro,         'LWI',     __RC__)
       call MAPL_GetPointer(import, shflux,      'SH',      __RC__)
       call MAPL_GetPointer(import, ustar,       'USTAR',   __RC__)
       call MAPL_GetPointer(import, z0h,         'Z0H',     __RC__)
       call MAPL_GetPointer(import, pblh,        'ZPBL',    __RC__)

       call MAPL_GetPointer(import, q_OH,        'OH',      __RC__)
       call MAPL_GetPointer(import, q_NO3,       'NO3',     __RC__)
       call MAPL_GetPointer(import, q_H2O2,      'H2O2',    __RC__)
       call MAPL_GetPointer(import, q_O3,        'O3',      __RC__)

       call MAPL_GetPointer(import, DMS_ocean, 'DMS_CONC_OCEAN', __RC__)

       call MAPL_GetPointer(import, SO2_emiss_bb,           'SO2_EMIS_FIRES',        __RC__)
       call MAPL_GetPointer(import, SO2_emiss_nonenergy,    'SO2_EMIS_NONENERGY',    __RC__)
       call MAPL_GetPointer(import, SO2_emiss_energy,       'SO2_EMIS_ENERGY',       __RC__)
       call MAPL_GetPointer(import, SO2_emiss_shipping,     'SO2_EMIS_SHIPPING',     __RC__)
       call MAPL_GetPointer(import, SO2_emiss_aviation_lto, 'SO2_EMIS_AIRCRAFT_LTO', __RC__)
       call MAPL_GetPointer(import, SO2_emiss_aviation_cds, 'SO2_EMIS_AIRCRAFT_CDS', __RC__)
       call MAPL_GetPointer(import, SO2_emiss_aviation_crs, 'SO2_EMIS_AIRCRAFT_CRS', __RC__)

       call MAPL_GetPointer(import, NH3_emiss,    'NH3_EMIS',      __RC__)
       call MAPL_GetPointer(import, NH3_emiss_bb, 'NH3_EMIS_FIRE', __RC__)

       call MAPL_GetPointer(import, SOAG_emiss,   'SOAG_EMIS',     __RC__)
   end if

   if (self%ocs_chem) then
       call MAPL_GetPointer(import, tropp,           'TROPP',      __RC__)
       call MAPL_GetPointer(import, q_O3p_STRATCHEM, 'O3P',        __RC__)
       call MAPL_GetPointer(import, q_OH_STRATCHEM,  'OHSTRAT',    __RC__)
       call MAPL_GetPointer(import, j_ocs,           'OCS_JRATE',  __RC__)
   end if

   if (self%voc_chem) then
       if (.not. associated(q_OH)) then
           call MAPL_GetPointer(import, q_OH,       'OH',         __RC__)
       end if

       call MAPL_GetPointer(import, co_biomass_voc, 'CO_BIOMASS_VOC', __RC__)
       call MAPL_GetPointer(import, co_bf_voc,      'CO_BF_VOC',      __RC__)
       call MAPL_GetPointer(import, co_fs_voc,      'CO_FS_VOC',      __RC__)
   end if

!  Get Exports
!  -------------
   if (self%mam_chem) then
       call MAPL_GetPointer(export, dry_dep_DMS,   'DRY_DEP_DMS',   __RC__)
       call MAPL_GetPointer(export, dry_dep_MSA,   'DRY_DEP_MSA',   __RC__)
       call MAPL_GetPointer(export, dry_dep_SO2,   'DRY_DEP_SO2',   __RC__)
       call MAPL_GetPointer(export, dry_dep_H2SO4, 'DRY_DEP_H2SO4', __RC__)
       call MAPL_GetPointer(export, dry_dep_NH3,   'DRY_DEP_NH3',   __RC__)
       call MAPL_GetPointer(export, dry_dep_SOAG,  'DRY_DEP_SOAG',  __RC__)

       call MAPL_GetPointer(export, ddt_DMS_gas,   'DDT_DMS_gas',   __RC__)
       call MAPL_GetPointer(export, ddt_MSA_gas,   'DDT_MSA_gas',   __RC__)
       call MAPL_GetPointer(export, ddt_SO2_gas,   'DDT_SO2_gas',   __RC__)
       call MAPL_GetPointer(export, ddt_H2SO4_gas, 'DDT_H2SO4_gas', __RC__)
       call MAPL_GetPointer(export, ddt_NH3_gas,   'DDT_NH3_gas',   __RC__)
       call MAPL_GetPointer(export, ddt_SOAG_gas,  'DDT_SOAG_gas',  __RC__)

       call MAPL_GetPointer(export, ddt_DMS_aq,    'DDT_DMS_aq',    __RC__)
       call MAPL_GetPointer(export, ddt_MSA_aq,    'DDT_MSA_aq',    __RC__)
       call MAPL_GetPointer(export, ddt_SO2_aq,    'DDT_SO2_aq',    __RC__)
       call MAPL_GetPointer(export, ddt_H2SO4_aq,  'DDT_H2SO4_aq',  __RC__)
       call MAPL_GetPointer(export, ddt_NH3_aq,    'DDT_NH3_aq',    __RC__)
       call MAPL_GetPointer(export, ddt_SOAG_aq,   'DDT_SOAG_aq',   __RC__)

       call MAPL_GetPointer(export, DMS_g_,        '_DMS_gas',      __RC__)
       call MAPL_GetPointer(export, MSA_g_,        '_MSA_gas',      __RC__)
       call MAPL_GetPointer(export, SO2_g_,        '_SO2_gas',      __RC__)
       call MAPL_GetPointer(export, H2SO4_g_,      '_H2SO4_gas',    __RC__)
       call MAPL_GetPointer(export, NH3_g_,        '_NH3_gas',      __RC__)
       call MAPL_GetPointer(export, SOAG_g_,       '_SOAG_gas',     __RC__)

       call MAPL_GetPointer(export, DMS_a_,        '_DMS_aq',       __RC__)
       call MAPL_GetPointer(export, MSA_a_,        '_MSA_aq',       __RC__)
       call MAPL_GetPointer(export, SO2_a_,        '_SO2_aq',       __RC__)
       call MAPL_GetPointer(export, H2SO4_a_,      '_H2SO4_aq',     __RC__)
       call MAPL_GetPointer(export, NH3_a_,        '_NH3_aq',       __RC__)
       call MAPL_GetPointer(export, SOAG_a_,       '_SOAG_aq',      __RC__)
   end if

!  Get Internals
!  -------------
   call MAPL_GetObjectFromGC(GC, mgState, __RC__)
   call MAPL_Get(mgState, INTERNAL_ESMF_STATE=internal, __RC__)

   if (self%mam_chem) then
       call MAPL_GetPointer(internal, q_DMS,    trim(comp_name)//'::'//'DMS',     __RC__)
       call MAPL_GetPointer(internal, q_MSA,    trim(comp_name)//'::'//'MSA',     __RC__)
       call MAPL_GetPointer(internal, q_SO2,    trim(comp_name)//'::'//'SO2',     __RC__)
       call MAPL_GetPointer(internal, q_H2SO4,  trim(comp_name)//'::'//'H2SO4',   __RC__)
       call MAPL_GetPointer(internal, q_NH3,    trim(comp_name)//'::'//'NH3',     __RC__)
       call MAPL_GetPointer(internal, q_SOAG,   trim(comp_name)//'::'//'SOAG',    __RC__)
   end if

   if (self%ocs_chem) then
       call MAPL_GetPointer(internal, q_OCS, trim(comp_name)//'::'//'OCS', __RC__)
   end if

   if (self%voc_chem) then
       call MAPL_GetPointer(internal, q_VOCanth, trim(comp_name)//'::'//'VOC', __RC__)
       call MAPL_GetPointer(internal, q_VOCbiob, trim(comp_name)//'::'//'VOCbiob', __RC__)
   end if


   call MAPL_TimerOn(mgState, '-EMISSIONS', __RC__)

   UPDATE_VOLCANIC_EMISSIONS: if (self%mam_chem) then
!  Update volcanic emissions if necessary (daily)
!  ----------------------------------------------
   if(self%nymd_volcanic_emiss .ne. nymd) then
       self%nymd_volcanic_emiss = nymd

       call GetVolcDailyTables(self%nymd_volcanic_emiss,       &
                               trim(self%volcanic_emiss_file), &
                               self%n_volcanoes,               &
                               self%volc_lat,                  &
                               self%volc_lon,                  &
                               self%volc_elev,                 &
                               self%volc_cloud,                &
                               self%volc_SO2,                  &
                               self%volc_start,                &
                               self%volc_end,                  &
                               __RC__)

       if (self%n_volcanoes > 0) then
           if (associated(self%volc_i)) deallocate(self%volc_i, __STAT__)
           allocate(self%volc_i(self%n_volcanoes), __STAT__)

           if (associated(self%volc_j)) deallocate(self%volc_j, __STAT__)
           allocate(self%volc_j(self%n_volcanoes), __STAT__)

           ! get indices for volcanic emissions
           call MAPL_GetHorzIJIndex(self%n_volcanoes,                      &
                                    self%volc_i, self%volc_j,              &
                                    grid = self%grid,                      &
                                    lon = self%volc_lon * (MAPL_PI/180.0), &
                                    lat = self%volc_lat * (MAPL_PI/180.0), &
                                    __RC__)
       end if
   end if
   end if UPDATE_VOLCANIC_EMISSIONS


!  Impose diurnal cycle to the offline OH and NO3 monthly mean fields
!  ------------------------------------------------------------------
   if (self%gas_phase_chem .or. self%voc_chem) then
       allocate(q_OH_(i1:i2,j1:j2,km),  __STAT__)
       q_OH_  = q_OH
   end if

   if (self%gas_phase_chem) then
       allocate(q_NO3_(i1:i2,j1:j2,km), __STAT__)
       q_NO3_ = q_NO3
   end if

   if (self%apply_diurnal_cycle .and. (self%mam_chem .or. self%voc_chem)) then
       ! find cos(SZA)
       doy = day_of_year(nymd)
       f_hour = ( real(nhms / 10000) * 3600 +         &
                  real(mod(nhms, 10000) / 100) * 60 + &
                  real(mod(nhms, 100)) ) / 3600


       ! want to find the sum of the cos(sza) for use in scaling OH diurnal variation
       allocate(sza(i1:i2,j1:j2),          __STAT__)
       allocate(cos_sza(i1:i2,j1:j2),      __STAT__)
       allocate(sum_cos_sza(i1:i2,j1:j2),  __STAT__)
       allocate(day_time(i1:i2,j1:j2),     __STAT__)
       allocate(night_time(i1:i2,j1:j2),   __STAT__)
       allocate(f_day_time(i1:i2,j1:j2),   __STAT__)
       allocate(f_night_time(i1:i2,j1:j2), __STAT__)


       x_hour = f_hour
       n_steps = 86400.0 / cdt
       sum_cos_sza(:,:) = 0.0
       day_time(:,:) = 0.0
       do n = 1, n_steps
           call solar_zenith_angle(doy, x_hour, (180.0/pi)*lons, (180.0/pi)*lats, sza, cos_sza)

           sum_cos_sza = sum_cos_sza + cos_sza

           x_hour = x_hour + cdt/3600.0
           if (x_hour > 24) x_hour = x_hour - 24

           ! find the daylight portion of the day
           where (cos_sza > 0.0)
               day_time = day_time + cdt
           end where
       end do

       night_time = 86400.0 - day_time

       call solar_zenith_angle(doy, f_hour, (180.0/pi)*lons, (180.0/pi)*lats, sza, cos_sza)

       where(sum_cos_sza > 0)
           f_day_time = (86400.0/cdt)*cos_sza / sum_cos_sza
       elsewhere
           f_day_time = 0.0
       end where


       ! scale OH
       do k = 1, km
           q_OH_(:,:,k) = q_OH_(:,:,k) * f_day_time(:,:)
       end do

       where(q_OH_ < 0.0) q_OH_ = 0.0


       ! set NO3 to 0 in sun lighten grid cells - average is
       ! distributed only over the night time portion
       where (cos_sza > 0 .or. night_time < tiny(0.0))
           f_night_time = 0.0
       elsewhere
           f_night_time = 86400.0 / night_time
       end where

       if (self%gas_phase_chem) then
           do k = 1, km
                q_NO3_(:,:,k) = q_NO3_(:,:,k) * f_night_time(:,:)
           end do

           where(q_NO3_ < 0.0) q_NO3_ = 0.0
       end if

       deallocate(sza,         __STAT__)
       deallocate(cos_sza,     __STAT__)
       deallocate(sum_cos_sza, __STAT__)
       deallocate(day_time,    __STAT__)
       deallocate(night_time,  __STAT__)
       deallocate(f_day_time,   __STAT__)
       deallocate(f_night_time, __STAT__)
   end if ! diurnal cycle of oxidants

   if (self%mam_chem .or. self%voc_chem) then
   call MAPL_GetPointer(export, ptr3d, 'OH', __RC__)
   if (associated(ptr3d)) then
       ptr3d = q_OH_
   end if
   end if

   if (self%mam_chem) then
   call MAPL_GetPointer(export, ptr3d, 'NO3', __RC__)
   if (associated(ptr3d)) then
       ptr3d = q_NO3_
   end if
   end if


!  If the H2O2 is from climatology, replenish it every 3 hours
!  -----------------------------------------------------------
   if (.not. using_GMI_H2O2 .and. self%aqu_phase_chem) then
       if (mod(nhms/10000, 3) == 0 .and. (nhms/10000*100 == nhms/100)) then
           self%h2o2 = q_H2O2
       end if
   end if


   UPDATE_CHEM_EMISSIONS: if (self%mam_chem) then
!
!  Ammonia emissions
!  -----------------
   call NH3_emissions(delp,         &
                      NH3_emiss,    &
                      NH3_emiss_bb, &
                      q_NH3,        &
                      cdt,          &
                      rc)


!  Sulfur emissions
!  ----------------
   allocate(flux_DMS(i1:i2,j1:j2),  __STAT__)

   call DMS_emissions(delp,        &
                      tskin,       &
                      u10n,        &
                      v10n,        &
                      fr_ocean,    &
                      DMS_ocean,   &
                      q_DMS,       &
                      flux_DMS,    &
                      cdt,         &
                      rc)

   call MAPL_GetPointer(export, ptr2d, 'EMIS_DMS', __RC__)
   if (associated(ptr2d)) then
       ! convert from 'mol-DMS m-2 s-1' to 'kg-DMS m-2 s-1'
       ptr2d = mw_DMS * flux_DMS
   end if

   call MAPL_GetPointer(export, ptr2d, 'EMIS_S_DMS', __RC__)
   if (associated(ptr2d)) then
       ! convert from 'mol-DMS m-2 s-1' to 'kg-S m-2 s-1'
       ptr2d = mw_S * flux_DMS
   end if

   deallocate(flux_DMS)


   allocate(SO2_emiss_total(i1:i2,j1:j2),        __STAT__)
   allocate(SO2_emiss_volc_expl(i1:i2,j1:j2),    __STAT__)
   allocate(SO2_emiss_volc_nonexpl(i1:i2,j1:j2), __STAT__)

   SO2_emiss_total        = 0.0
   SO2_emiss_volc_expl    = 0.0
   SO2_emiss_volc_nonexpl = 0.0

   ! 3D aircraft emissions
   ! TODO: move vertical distribution of 2D-layered emissions here
   ! calc_aviation_emissions(...)

   ! 3D volcanic emissions
   ! TODO: move calculation of area mean 3D volcanic emissions here
   ! calc_volcanic_emissions(...)

   call SO2_emissions(delp,                   &
                      zle,                    &
                      density_air,            &
                      SO2_emiss_bb,           &
                      SO2_emiss_nonenergy,    &
                      SO2_emiss_energy,       &
                      SO2_emiss_shipping,     &
                      SO2_emiss_aviation_lto, &
                      SO2_emiss_aviation_cds, &
                      SO2_emiss_aviation_crs, &
                      self%aviation_layers,   &
                      self%n_volcanoes,       &
                      self%volc_elev,         &
                      self%volc_cloud,        &
                      self%volc_SO2,          &
                      self%volc_start, self%volc_end,   &
                      self%volc_i, self%volc_j, &
                      SO2_emiss_volc_expl,    &
                      SO2_emiss_volc_nonexpl, &
                      SO2_emiss_total,        &
                      q_SO2,                  &
                      cell_area,              &
                      cdt,                    &
                      nymd,                   &
                      nhms,                   &
                      rc)


   call MAPL_GetPointer(export, ptr2d, 'EMIS_SO2', __RC__)
   if (associated(ptr2d)) then
       ptr2d = SO2_emiss_total
   end if

   call MAPL_GetPointer(export, ptr2d, 'EMIS_S_SO2', __RC__)
   if (associated(ptr2d)) then
       ptr2d = (mw_S/mw_SO2) * SO2_emiss_total
   end if

   call MAPL_GetPointer(export, ptr2d, 'EMIS_SO2_EXV', __RC__)
   if (associated(ptr2d)) then
       ptr2d = SO2_emiss_volc_expl
   end if

   call MAPL_GetPointer(export, ptr2d, 'EMIS_SO2_NXV', __RC__)
   if (associated(ptr2d)) then
       ptr2d = SO2_emiss_volc_nonexpl
   end if

   call MAPL_GetPointer(export, ptr2d, 'EMIS_SO2_VOLC', __RC__)
   if (associated(ptr2d)) then
       ptr2d = SO2_emiss_volc_expl + SO2_emiss_volc_nonexpl
   end if


   deallocate(SO2_emiss_total,        __STAT__)
   deallocate(SO2_emiss_volc_expl,    __STAT__)
   deallocate(SO2_emiss_volc_nonexpl, __STAT__)


!
!  SOAG emissions
!  --------------
   call SOAG_emissions(delp,       &
                       SOAG_emiss, &
                       q_SOAG,     &
                       cdt,        &
                       rc)

   end if UPDATE_CHEM_EMISSIONS


   UPDATE_VOC_EMISSIONS: if (self%voc_chem) then

      call VOC_Emissions(delp,                       &
                         self%voc_BiomassBurnFactor, &
                         self%voc_AnthroFactor,      &
                         co_biomass_voc,             &
                         co_bf_voc,                  &
                         co_fs_voc,                  &
                         self%voc_MW,                &
                         q_VOCanth, q_VOCbiob,       &
                         cdt,                        &
                         rc)

   end if UPDATE_VOC_EMISSIONS


   call MAPL_TimerOff(mgState, '-EMISSIONS', __RC__)


!  Dry deposition
!  --------------
   UPDATE_CHEM_DRY_DEP: if (self%mam_chem) then

   allocate(dry_dep_frequency(i1:i2,j1:j2), __STAT__)
   allocate(dq(i1:i2,j1:j2),   __STAT__)

   call DryDepositionGOCART(i1, i2, j1, j2, km, &
                            temperature, density_air, zle, oro, ustar, &
                            pblh, shflux, z0h, dry_dep_frequency, rc)

   ! DMS - no dry dep
   if (associated(dry_dep_DMS)) &
       dry_dep_DMS = 0.0

   ! MSA
   dq = -q_MSA(i1:i2,j1:j2,km) * (1.0 - exp(-dry_dep_frequency * cdt))

   q_MSA(i1:i2,j1:j2,km) = q_MSA(i1:i2,j1:j2,km) + dq

   if (associated(dry_dep_MSA)) &
       dry_dep_MSA = -dq / mw_air / cdt

   ! SO2
   where (abs(oro - ORO_OCEAN) < 0.5)  ! oro has descrete values 0, 1 or 2
       dq = -q_SO2(i1:i2,j1:j2,km) * (1.0 - exp(-(10.0 * dry_dep_frequency) * cdt))
   elsewhere
       dq = -q_SO2(i1:i2,j1:j2,km) * (1.0 - exp(-( 3.0 * dry_dep_frequency) * cdt))
   end where

   q_SO2(i1:i2,j1:j2,km) = q_SO2(i1:i2,j1:j2,km) + dq

   if (associated(dry_dep_SO2)) &
       dry_dep_SO2 = -dq / mw_air / cdt

   ! H2SO4
   where (abs(oro - ORO_OCEAN) < 0.5)  ! oro has descrete values 0, 1 or 2
       dq = -q_H2SO4(i1:i2,j1:j2,km) * (1.0 - exp(-(10.0 * dry_dep_frequency) * cdt))
   elsewhere
       dq = -q_H2SO4(i1:i2,j1:j2,km) * (1.0 - exp(-( 3.0 * dry_dep_frequency) * cdt))
   end where

   q_H2SO4(i1:i2,j1:j2,km) = q_H2SO4(i1:i2,j1:j2,km) + dq

   if (associated(dry_dep_H2SO4)) &
       dry_dep_H2SO4 = -dq / mw_air / cdt

   ! NH3
   dq = -q_NH3(i1:i2,j1:j2,km) * (1.0 - exp(-dry_dep_frequency * cdt))

   q_NH3(i1:i2,j1:j2,km) = q_NH3(i1:i2,j1:j2,km) + dq

   if (associated(dry_dep_NH3)) &
       dry_dep_NH3 = -dq / mw_air / cdt

   ! lumped precursor VOC (SOA gas)
   dq = -q_SOAG(i1:i2,j1:j2,km) * (1.0 - exp(-dry_dep_frequency * cdt))

   q_SOAG(i1:i2,j1:j2,km) = q_SOAG(i1:i2,j1:j2,km) + dq

   if (associated(dry_dep_SOAG)) &
       dry_dep_SOAG = -dq / mw_air / cdt

   deallocate(dry_dep_frequency, __STAT__)
   deallocate(dq, __STAT__)

   end if UPDATE_CHEM_DRY_DEP


!  Ensure positive values only
!  ---------------------------
   if (self%mam_chem) then
       where (q_NH3   < tiny(0.0))    q_NH3   = tiny(0.0)
       where (q_DMS   < tiny(0.0))    q_DMS   = tiny(0.0)
       where (q_MSA   < tiny(0.0))    q_MSA   = tiny(0.0)
       where (q_SO2   < tiny(0.0))    q_SO2   = tiny(0.0)
       where (q_H2SO4 < tiny(0.0))    q_H2SO4 = tiny(0.0)
       where (q_NO3_  < tiny(0.0))    q_NO3_  = tiny(0.0)
   end if

   if (self%mam_chem .or. self%voc_chem) then
       where (q_OH_   < tiny(0.0))    q_OH_   = tiny(0.0)
   end if


   call MAPL_TimerOn(mgState, '-CHEMISTRY', __RC__)

!  Gas-phase chemistry
!  -------------------
   call MAPL_TimerOn(mgState, '--CHEMISTRY_GAS', __RC__)

   UPDATE_CHEM_GAS_PHASE: if (self%mam_chem) then

   if (associated(DMS_g_))   DMS_g_   = q_DMS
   if (associated(MSA_g_))   MSA_g_   = q_MSA
   if (associated(SO2_g_))   SO2_g_   = q_SO2
   if (associated(H2SO4_g_)) H2SO4_g_ = q_H2SO4
   if (associated(NH3_g_))   NH3_g_   = q_NH3
   if (associated(SOAG_g_))  SOAG_g_  = q_SOAG

   allocate(cpl_NH3(i1:i2,j1:j2),   __STAT__)
   allocate(cpl_DMS(i1:i2,j1:j2),   __STAT__)
   allocate(cpl_MSA(i1:i2,j1:j2),   __STAT__)
   allocate(cpl_SO2(i1:i2,j1:j2),   __STAT__)
   allocate(cpl_H2SO4(i1:i2,j1:j2), __STAT__)

   cpl_NH3   = 0.0
   cpl_DMS   = 0.0
   cpl_MSA   = 0.0
   cpl_SO2   = 0.0
   cpl_H2SO4 = 0.0

   if (self%gas_phase_chem) then
       ! tendencies due to gas phase chemistry
       if (associated(ddt_DMS_gas))   ddt_DMS_gas   = q_DMS
       if (associated(ddt_MSA_gas))   ddt_MSA_gas   = q_MSA
       if (associated(ddt_SO2_gas))   ddt_SO2_gas   = q_SO2
       if (associated(ddt_H2SO4_gas)) ddt_H2SO4_gas = q_H2SO4
       if (associated(ddt_NH3_gas))   ddt_NH3_gas   = q_NH3
       if (associated(ddt_SOAG_gas))  ddt_SOAG_gas  = q_SOAG

       call gas_chemistry(ple,           &
                          temperature,   &
                          density_air,   &
                          q_NH3,         &
                          q_DMS,         &
                          q_MSA,         &
                          q_SO2,         &
                          q_H2SO4,       &
                          q_OH_,         &
                          q_NO3_,        &
                          cpl_NH3,       &
                          cpl_DMS,       &
                          cpl_MSA,       &
                          cpl_SO2,       &
                          cpl_H2SO4,     &
                          cdt,           &
                          rc)

       ! tendencies due to gas phase chemistry
       if (associated(ddt_DMS_gas))   ddt_DMS_gas   = (q_DMS   - ddt_DMS_gas)   / cdt
       if (associated(ddt_MSA_gas))   ddt_MSA_gas   = (q_MSA   - ddt_MSA_gas)   / cdt
       if (associated(ddt_SO2_gas))   ddt_SO2_gas   = (q_SO2   - ddt_SO2_gas)   / cdt
       if (associated(ddt_H2SO4_gas)) ddt_H2SO4_gas = (q_H2SO4 - ddt_H2SO4_gas) / cdt
       if (associated(ddt_NH3_gas))   ddt_NH3_gas   = (q_NH3   - ddt_NH3_gas)   / cdt
       if (associated(ddt_SOAG_gas))  ddt_SOAG_gas  = (q_SOAG  - ddt_SOAG_gas)  / cdt
   else
       ! set to zero the tendencies due to gas phase chemistry
       if (associated(ddt_DMS_gas))   ddt_DMS_gas   = 0.0
       if (associated(ddt_MSA_gas))   ddt_MSA_gas   = 0.0
       if (associated(ddt_SO2_gas))   ddt_SO2_gas   = 0.0
       if (associated(ddt_H2SO4_gas)) ddt_H2SO4_gas = 0.0
       if (associated(ddt_NH3_gas))   ddt_NH3_gas   = 0.0
       if (associated(ddt_SOAG_gas))  ddt_SOAG_gas  = 0.0
   end if

   ! column integrated P+L tendencies due to gas phase chemistry
   call MAPL_GetPointer(export, ptr2d, 'CPL_DMS_gas', __RC__)
   if (associated(ptr2d)) then
       ptr2d = cpl_DMS * mw_DMS      ! 'kg-DMS m-2 s-1'
   end if

   call MAPL_GetPointer(export, ptr2d, 'CPL_MSA_gas', __RC__)
   if (associated(ptr2d)) then
       ptr2d = cpl_MSA * mw_MSA      ! 'kg-MSA m-2 s-1'
   end if

   call MAPL_GetPointer(export, ptr2d, 'CPL_SO2_gas', __RC__)
   if (associated(ptr2d)) then
       ptr2d = cpl_SO2 * mw_SO2      ! 'kg-SO2 m-2 s-1'
   end if

   call MAPL_GetPointer(export, ptr2d, 'CPL_H2SO4_gas', __RC__)
   if (associated(ptr2d)) then
       ptr2d = cpl_H2SO4 * mw_H2SO4  ! 'kg-H2SO4 m-2 s-1'
   end if

   ! column integrated P+L tendencies due to gas phase chemistry
   call MAPL_GetPointer(export, ptr2d, 'CPL_S_DMS_gas', __RC__)
   if (associated(ptr2d)) then
       ptr2d = cpl_DMS * mw_S        ! 'kg-S m-2 s-1'
   end if

   call MAPL_GetPointer(export, ptr2d, 'CPL_S_MSA_gas', __RC__)
   if (associated(ptr2d)) then
       ptr2d = cpl_MSA * mw_S        ! 'kg-S m-2 s-1'
   end if

   call MAPL_GetPointer(export, ptr2d, 'CPL_S_SO2_gas', __RC__)
   if (associated(ptr2d)) then
       ptr2d = cpl_SO2 * mw_S        ! 'kg-S m-2 s-1'
   end if

   call MAPL_GetPointer(export, ptr2d, 'CPL_S_H2SO4_gas', __RC__)
   if (associated(ptr2d)) then
       ptr2d = cpl_H2SO4 * mw_S      ! 'kg-S m-2 s-1'
   end if

   call MAPL_GetPointer(export, ptr2d, 'CPL_S_gas', __RC__)
   if (associated(ptr2d)) then
       ptr2d = (cpl_DMS + cpl_MSA + cpl_SO2 + cpl_H2SO4) * mw_S      ! 'kg-S m-2 s-1'
   end if


   call MAPL_GetPointer(export, ptr2d, 'CPL_NH3_gas', __RC__)
   if (associated(ptr2d)) then
       ptr2d = cpl_NH3 * mw_NH3      ! 'kg-NH3 m-2 s-1'
   end if

   ! column integrated P+L tendencies due to gas phase chemistry
   call MAPL_GetPointer(export, ptr2d, 'CPL_N_NH3_gas', __RC__)
   if (associated(ptr2d)) then
       ptr2d = cpl_NH3 * mw_N        ! 'kg-N m-2 s-1'
   end if

   call MAPL_GetPointer(export, ptr2d, 'CPL_N_gas', __RC__)
   if (associated(ptr2d)) then
       ptr2d = (cpl_NH3) * mw_N      ! 'kg-N m-2 s-1'
   end if

   deallocate(cpl_NH3,   __STAT__)
   deallocate(cpl_DMS,   __STAT__)
   deallocate(cpl_MSA,   __STAT__)
   deallocate(cpl_SO2,   __STAT__)
   deallocate(cpl_H2SO4, __STAT__)

   end if UPDATE_CHEM_GAS_PHASE

   call MAPL_TimerOff(mgState, '--CHEMISTRY_GAS', __RC__)


!  Aqueous-phase chemistry
!  -----------------------
   call MAPL_TimerOn(mgState, '--CHEMISTRY_AQUEOUS', __RC__)

   UPDATE_CHEM_AQU_PHASE: if (self%mam_chem) then

   allocate(pSO4_aq(i1:i2,j1:j2,km),       __STAT__)
   allocate(pNH4_aq(i1:i2,j1:j2,km),       __STAT__)
   allocate(pSO4_aq_SO2(i1:i2,j1:j2,km),   __STAT__)
   allocate(pSO4_aq_H2SO4(i1:i2,j1:j2,km), __STAT__)
   allocate(pNH4_aq_NH3(i1:i2,j1:j2,km),   __STAT__)

   pSO4_aq       = 0.0
   pNH4_aq       = 0.0
   pSO4_aq_SO2   = 0.0
   pSO4_aq_H2SO4 = 0.0
   pNH4_aq_NH3   = 0.0

   if (associated(DMS_a_))   DMS_a_   = q_DMS
   if (associated(MSA_a_))   MSA_a_   = q_MSA
   if (associated(SO2_a_))   SO2_a_   = q_SO2
   if (associated(H2SO4_a_)) H2SO4_a_ = q_H2SO4
   if (associated(NH3_a_))   NH3_a_   = q_NH3
   if (associated(SOAG_a_))  SOAG_a_  = q_SOAG

   if (self%aqu_phase_chem) then
       ! tendencies due to aqueous chemistry
       if (associated(ddt_DMS_aq))   ddt_DMS_aq   = q_DMS
       if (associated(ddt_MSA_aq))   ddt_MSA_aq   = q_MSA
       if (associated(ddt_SO2_aq))   ddt_SO2_aq   = q_SO2
       if (associated(ddt_H2SO4_aq)) ddt_H2SO4_aq = q_H2SO4
       if (associated(ddt_NH3_aq))   ddt_NH3_aq   = q_NH3
       if (associated(ddt_SOAG_aq))  ddt_SOAG_aq  = q_SOAG

       call aqu_chemistry_fast(ple,           &
                               temperature,   &
                               density_air,   &
                               lwc,           &
                               0.5e-7,        &
                               fcld,          &
                               q_NH3,         &
                               q_SO2,         &
                               q_H2SO4,       &
                               self%h2o2,     &
                               pSO4_aq,       &
                               pNH4_aq,       &
                               pSO4_aq_SO2,   &
                               pSO4_aq_H2SO4, &
                               pNH4_aq_NH3,   &
                               cdt,           &
                               rc)

       ! tendencies due to aqueous chemistry
       if (associated(ddt_DMS_aq))   ddt_DMS_aq   = (q_DMS   - ddt_DMS_aq)   / cdt
       if (associated(ddt_MSA_aq))   ddt_MSA_aq   = (q_MSA   - ddt_MSA_aq)   / cdt
       if (associated(ddt_SO2_aq))   ddt_SO2_aq   = (q_SO2   - ddt_SO2_aq)   / cdt
       if (associated(ddt_H2SO4_aq)) ddt_H2SO4_aq = (q_H2SO4 - ddt_H2SO4_aq) / cdt
       if (associated(ddt_NH3_aq))   ddt_NH3_aq   = (q_NH3   - ddt_NH3_aq)   / cdt
       if (associated(ddt_SOAG_aq))  ddt_SOAG_aq  = (q_SOAG  - ddt_SOAG_aq)  / cdt
   else
       ! set to zero the tendencies due to aqueous chemistry
       if (associated(ddt_DMS_aq))   ddt_DMS_aq   = 0.0
       if (associated(ddt_MSA_aq))   ddt_MSA_aq   = 0.0
       if (associated(ddt_SO2_aq))   ddt_SO2_aq   = 0.0
       if (associated(ddt_H2SO4_aq)) ddt_H2SO4_aq = 0.0
       if (associated(ddt_NH3_aq))   ddt_NH3_aq   = 0.0
       if (associated(ddt_SOAG_aq))  ddt_SOAG_aq  = 0.0
   end if


   ! total production in aqueous phase
   call MAPL_GetPointer(export, ptr3d,  'pSO4_aq',  __RC__)
   if (associated(ptr3d)) ptr3d = pSO4_aq

   call MAPL_GetPointer(export, ptr3d,  'pNH4_aq',  __RC__)
   if (associated(ptr3d)) ptr3d = pNH4_aq

   ! contributions of production pathways in aqueous phase
   call MAPL_GetPointer(export, ptr3d,  'pSO4_aq_SO2',  __RC__)
   if (associated(ptr3d)) ptr3d = pSO4_aq_SO2

   call MAPL_GetPointer(export, ptr3d,  'pSO4_aq_H2SO4',  __RC__)
   if (associated(ptr3d)) ptr3d = pSO4_aq_H2SO4

   call MAPL_GetPointer(export, ptr3d,  'pNH4_aq_NH3',  __RC__)
   if (associated(ptr3d)) ptr3d = pNH4_aq_NH3

   deallocate(pSO4_aq,       __STAT__)
   deallocate(pNH4_aq,       __STAT__)
   deallocate(pSO4_aq_SO2,   __STAT__)
   deallocate(pSO4_aq_H2SO4, __STAT__)
   deallocate(pNH4_aq_NH3,   __STAT__)


!!     call aqu_chemistry(ple,                   &
!!                       temperature,            &
!!                       density_air,            &
!!                       lwc,                    &
!!                       fcld,                   &
!!                       q_NH3,                  &
!!                       q_SO2,                  &
!!                       q_SO4_aq,               &
!!                       q_NH4_aq,               &
!!                       q_H2O2,                 &
!!                       q_O3,                   &
!!                       cdt,                    &
!!                       self%aqu_solver_max_dt, &
!!                       0.5e-7,                 &
!!                       rc)

   end if UPDATE_CHEM_AQU_PHASE

   call MAPL_TimerOff(mgState, '--CHEMISTRY_AQUEOUS', __RC__)


   call MAPL_TimerOn(mgState, '--CHEMISTRY_VOC', __RC__)

   ! If doing VOC chemistry by OH to create OA
   ! -----------------------------------------
   ! Note that we do not update the OH concentration
   ! Could (should?) integrate this with gas_chemistry below
   ! Right now the VOC calculation is optional, specified in
   ! in GEOS_AchemGridComp.rc, but the export is coupled to
   ! GOCART if both GOCART and ACHEM are running, so fill
   ! in with zero in case it is requested.

   OPTIONAL_VOC_CHEMISTRY: if (self%voc_chem) then

       allocate(q_OAanth(i1:i2,j1:j2,1:km), &
                q_OAbiob(i1:i2,j1:j2,1:km), &
                q_OAanthmmrd(i1:i2,j1:j2,1:km), &
                q_OAbiobmmrd(i1:i2,j1:j2,1:km), __STAT__)

       allocate(dVOC(i1:i2,j1:j2,1:km), &
                dOAanth(i1:i2,j1:j2,1:km), &
                dOAbiob(i1:i2,j1:j2,1:km), &
                fanth(i1:i2,j1:j2,1:km), &
                rk_OA_OH(i1:i2,j1:j2,1:km),__STAT__)

       q_OAanth     = 0.0
       q_OAanthmmrd = 0.0
       q_OAbiob     = 0.0
       q_OAbiobmmrd = 0.0

       where (q_VOCanth < tiny(0.0)) q_VOCanth = tiny(0.0)
       where (q_VOCbiob < tiny(0.0)) q_VOCbiob = tiny(0.0)

       ! rate coefficient from Kim et al 2015
       rk_OA_OH = 1.25d-11*N_avog*q_OH_*density_air/mw_air*(1.0e-6)*cdt
       dVOC     = (q_VOCanth + q_VOCbiob)*(1.0-exp(-rk_OA_OH))     ! Loss of VOC (mol/mol air)
       dOAanth  = 0.0
       dOAbiob  = 0.0

       where (dVOC > 1.e-32)
          fanth     = q_VOCanth / (q_VOCanth + q_VOCbiob)          ! Anthropogenic fraction of total VOC
          dOAanth   = dVOC *fanth                                  ! Production of OA (mol/mol air)
          q_VOCanth = q_VOCanth - dVOC * fanth                     ! Update VOC (mol/mol air)
          dOAbiob   = dVOC * (1.0 -fanth)                          ! Production of OA (mol/mol air)
          q_VOCbiob = q_VOCbiob- dVOC * (1.0 -fanth)               ! Update VOC (mol/mol air)
       endwhere

       where (q_VOCanth<tiny(0.0)) q_VOCanth = tiny(0.0)
       where (q_VOCbiob<tiny(0.0)) q_VOCbiob = tiny(0.0)

       dOAanth = dOAanth*self%soa_MW/mw_air/cdt  ! kg kg-1 s-1
       dOAbiob = dOAbiob*self%soa_MW/mw_air/cdt  ! kg kg-1 s-1

       q_OAanthmmrd = dOAanth*86400.0            ! kg kg-1 day-1
       q_OAbiobmmrd = dOAbiob*86400.0            ! kg kg-1 day-1

       q_OAanth = dOAanth*density_air            ! kg m-3 s-1
       q_OAbiob = dOAbiob*density_air            ! kg m-3 s-1

       call MAPL_GetPointer(export, ptr3d, 'pSOA_ANTHRO_VOC', __RC__)
       if (associated(ptr3d)) ptr3d = q_OAanth

       call MAPL_GetPointer(export, ptr3d,  'pSOA_BIOB_VOC', __RC__)
       if (associated(ptr3d)) ptr3d = q_OAbiob

       call MAPL_GetPointer(export, ptr3d, 'pSOA_ANTHRO_VOC_MMRday', __RC__)
       if (associated(ptr3d)) ptr3d = q_OAanthmmrd

       call MAPL_GetPointer(export, ptr3d, 'pSOA_BIOB_VOC_MMRday', __RC__)
       if (associated(ptr3d)) ptr3d = q_OAbiobmmrd

       deallocate(q_OAanth, q_OAbiob, q_OAanthmmrd, q_OAbiobmmrd, __STAT__)
       deallocate(dVOC, dOAanth, dOAbiob, fanth, rk_OA_OH, __STAT__)
   end if OPTIONAL_VOC_CHEMISTRY

   call MAPL_TimerOff(mgState, '--CHEMISTRY_VOC', __RC__)


   call MAPL_TimerOn(mgState, '--CHEMISTRY_OCS', __RC__)

   OCS_CHEMISTRY_STRATOSPHERE: if (self%ocs_chem) then

       allocate(pSO2_OCS(i1:i2,j1:j2,km),      &
                pSO2_OCS_OH(i1:i2,j1:j2,km),   &
                pSO2_OCS_O3P(i1:i2,j1:j2,km),  &
                pSO2_OCS_jOCS(i1:i2,j1:j2,km), &
                lOCS(i1:i2,j1:j2,km),          &
                lOCS_OH(i1:i2,j1:j2,km),       &
                lOCS_O3p(i1:i2,j1:j2,km),      &
                lOCS_jOCS(i1:i2,j1:j2,km),     &
                __STAT__)

       call ocs_chemistry(ple,             &
                          temperature,     &
                          density_air,     &
                          self%ocs_surface_vmr, &
                          tropp,           &
                          q_OCS,           &
                          q_OH_STRATCHEM,  &
                          q_O3p_STRATCHEM, &
                          j_ocs,           & !photolysis rates
                          pSO2_OCS,        &
                          pSO2_OCS_OH,     &
                          pSO2_OCS_O3p,    &
                          pSO2_OCS_jOCS,   &
                          lOCS,            &
                          lOCS_OH,         &
                          lOCS_O3p,        &
                          lOCS_jOCS,       &
                          cdt,             &
                          rc)

      call MAPL_GetPointer(export, ptr3d, 'pSO2_OCS', __RC__)
       if (associated(ptr3d)) ptr3d = pSO2_OCS
      call MAPL_GetPointer(export, ptr3d, 'pSO2_OCS_OH', __RC__)
       if (associated(ptr3d)) ptr3d = pSO2_OCS_OH
      call MAPL_GetPointer(export, ptr3d, 'pSO2_OCS_O3p', __RC__)
       if (associated(ptr3d)) ptr3d = pSO2_OCS_O3p
      call MAPL_GetPointer(export, ptr3d, 'pSO2_OCS_jOCS', __RC__)
       if (associated(ptr3d)) ptr3d = pSO2_OCS_jOCS
      call MAPL_GetPointer(export, ptr3d, 'lOCS', __RC__)
       if (associated(ptr3d)) ptr3d = lOCS
      call MAPL_GetPointer(export, ptr3d, 'lOCS_OH', __RC__)
       if (associated(ptr3d)) ptr3d = lOCS_OH
      call MAPL_GetPointer(export, ptr3d, 'lOCS_O3p', __RC__)
       if (associated(ptr3d)) ptr3d = lOCS_O3p
      call MAPL_GetPointer(export, ptr3d, 'lOCS_jOCS', __RC__)
       if (associated(ptr3d)) ptr3d = lOCS_jOCS

      deallocate(pSO2_OCS, pSO2_OCS_OH, pSO2_OCS_O3p, pSO2_OCS_jOCS, &
                 lOCS, lOCS_OH, lOCS_O3p, lOCS_jOCS, __STAT__)

   end if OCS_CHEMISTRY_STRATOSPHERE

   call MAPL_TimerOff(mgState, '--CHEMISTRY_OCS', __RC__)


   ! could be used by the full and/or VOC chemistry mechanisms
   if (allocated(q_OH_ )) deallocate(q_OH_ , __STAT__)
   if (allocated(q_NO3_)) deallocate(q_NO3_, __STAT__)

   call MAPL_TimerOff(mgState, '-CHEMISTRY', __RC__)


!  Diagnostics
!  -----------
   UPDATE_CHEM_DIAGNOSTICS: if (self%mam_chem) then

!  3D concentrations in units '# molecules-constituent m-3'
!  -----------------------------------------------
   call MAPL_GetPointer(export, ptr3d, 'CONC_DMS', __RC__)
   if (associated(ptr3d)) then
       ptr3d = (N_avog / mw_air) * density_air(:,:,:) * q_DMS(:,:,:)
   end if

   call MAPL_GetPointer(export, ptr3d, 'CONC_MSA', __RC__)
   if (associated(ptr3d)) then
       ptr3d = (N_avog / mw_air) * density_air(:,:,:) * q_MSA(:,:,:)
   end if

   call MAPL_GetPointer(export, ptr3d, 'CONC_SO2', __RC__)
   if (associated(ptr3d)) then
       ptr3d = (N_avog / mw_air) * density_air(:,:,:) * q_SO2(:,:,:)
   end if

   call MAPL_GetPointer(export, ptr3d, 'CONC_H2SO4', __RC__)
   if (associated(ptr3d)) then
       ptr3d = (N_avog / mw_air) * density_air(:,:,:) * q_H2SO4(:,:,:)
   end if

   call MAPL_GetPointer(export, ptr3d, 'CONC_NH3', __RC__)
   if (associated(ptr3d)) then
       ptr3d = (N_avog / mw_air) * density_air(:,:,:) * q_NH3(:,:,:)
   end if

   call MAPL_GetPointer(export, ptr3d, 'CONC_SOAG', __RC__)
   if (associated(ptr3d)) then
       ptr3d = (N_avog / mw_air) * density_air(:,:,:) * q_SOAG(:,:,:)
   end if


!  Near-surface concentrations in units '# molecules-constituent m-3'
!  --------------------------------------------------------
   call MAPL_GetPointer(export, ptr2d, 'SFC_CONC_DMS', __RC__)
   if (associated(ptr2d)) then
       ptr2d = (N_avog / mw_air) * density_air(:,:,km) * q_DMS(:,:,km)
   end if

   call MAPL_GetPointer(export, ptr2d, 'SFC_CONC_MSA', __RC__)
   if (associated(ptr2d)) then
       ptr2d = (N_avog / mw_air) * density_air(:,:,km) * q_MSA(:,:,km)
   end if

   call MAPL_GetPointer(export, ptr2d, 'SFC_CONC_SO2', __RC__)
   if (associated(ptr2d)) then
       ptr2d = (N_avog / mw_air) * density_air(:,:,km) * q_SO2(:,:,km)
   end if

   call MAPL_GetPointer(export, ptr2d, 'SFC_CONC_H2SO4', __RC__)
   if (associated(ptr2d)) then
       ptr2d = (N_avog / mw_air) * density_air(:,:,km) * q_H2SO4(:,:,km)
   end if

   call MAPL_GetPointer(export, ptr2d, 'SFC_CONC_NH3', __RC__)
   if (associated(ptr2d)) then
       ptr2d = (N_avog / mw_air) * density_air(:,:,km) * q_NH3(:,:,km)
   end if

   call MAPL_GetPointer(export, ptr2d, 'SFC_CONC_SOAG', __RC__)
   if (associated(ptr2d)) then
       ptr2d = (N_avog / mw_air) * density_air(:,:,km) * q_SOAG(:,:,km)
   end if


!  Column mass densities
!  ---------------------
   allocate(cmd_S(i1:i2,j1:j2), __STAT__)
   cmd_S = 0.0

   ! column mass density of DMS in 'kg-DMS m-2' and 'kg-S m-2'
   allocate(cmd_DMS(i1:i2,j1:j2), __STAT__)
   cmd_DMS = ((mw_DMS/mw_air) / g_earth) * sum(q_DMS*delp, dim=3)
   cmd_S = cmd_S + (mw_S/mw_DMS) * cmd_DMS

   call MAPL_GetPointer(export, ptr2d, 'CMD_DMS', __RC__)
   if (associated(ptr2d)) ptr2d = cmd_DMS

   call MAPL_GetPointer(export, ptr2d, 'CMD_S_DMS', __RC__)
   if (associated(ptr2d)) ptr2d = (mw_S/mw_DMS) * cmd_DMS

   deallocate(cmd_DMS, __STAT__)

   ! column mass density of MSA in 'kg-MSA m-2' and 'kg-S m-2'
   allocate(cmd_MSA(i1:i2,j1:j2), __STAT__)
   cmd_MSA = ((mw_MSA/mw_air) / g_earth) * sum(q_MSA*delp, dim=3)
   cmd_S = cmd_S + (mw_S/mw_MSA) * cmd_MSA

   call MAPL_GetPointer(export, ptr2d, 'CMD_MSA', __RC__)
   if (associated(ptr2d)) ptr2d = cmd_MSA

   call MAPL_GetPointer(export, ptr2d, 'CMD_S_MSA', __RC__)
   if (associated(ptr2d)) ptr2d = (mw_S/mw_MSA) * cmd_MSA

   deallocate(cmd_MSA, __STAT__)


   ! column mass density of SO2 in 'kg-SO2 m-2' and 'kg-S m-2'
   allocate(cmd_SO2(i1:i2,j1:j2), __STAT__)
   cmd_SO2 = ((mw_SO2/mw_air) / g_earth) * sum(q_SO2*delp, dim=3)
   cmd_S = cmd_S + (mw_S/mw_SO2) * cmd_SO2

   call MAPL_GetPointer(export, ptr2d, 'CMD_SO2', __RC__)
   if (associated(ptr2d)) ptr2d = cmd_SO2

   call MAPL_GetPointer(export, ptr2d, 'CMD_S_SO2', __RC__)
   if (associated(ptr2d)) ptr2d = (mw_S/mw_SO2) * cmd_so2

   deallocate(cmd_SO2, __STAT__)

   ! column mass density of H2SO4 in 'kg-H2SO4 m-2' and 'kg-S m-2'
   allocate(cmd_H2SO4(i1:i2,j1:j2), __STAT__)
   cmd_H2SO4 = ((mw_H2SO4/mw_air) / g_earth) * sum(q_H2SO4*delp, dim=3)
   cmd_S = cmd_S + (mw_S/mw_H2SO4) * cmd_H2SO4

   call MAPL_GetPointer(export, ptr2d, 'CMD_H2SO4', __RC__)
   if (associated(ptr2d)) ptr2d = cmd_H2SO4

   call MAPL_GetPointer(export, ptr2d, 'CMD_S_H2SO4', __RC__)
   if (associated(ptr2d)) ptr2d = (mw_S/mw_H2SO4) * cmd_H2SO4

   deallocate(cmd_H2SO4, __STAT__)

   ! column mass density of sulfur (S) in trace gases (DMS, MSA, SO2, H2SO4)
   call MAPL_GetPointer(export, ptr2d, 'CMD_S_GAS', __RC__)
   if (associated(ptr2d)) ptr2d = cmd_S

   deallocate(cmd_S, __STAT__)


   allocate(cmd_N(i1:i2,j1:j2), __STAT__)
   cmd_N = 0.0

   allocate(cmd_NH3(i1:i2,j1:j2), __STAT__)
   cmd_NH3 = ((mw_NH3/mw_air) / g_earth) * sum(q_NH3*delp, dim=3)
   cmd_N = cmd_N + (mw_N/mw_NH3) * cmd_NH3

   call MAPL_GetPointer(export, ptr2d, 'CMD_NH3', __RC__)
   if (associated(ptr2d)) ptr2d = cmd_NH3

   call MAPL_GetPointer(export, ptr2d, 'CMD_N_NH3', __RC__)
   if (associated(ptr2d)) ptr2d = (mw_N/mw_NH3) * cmd_NH3

   deallocate(cmd_NH3, __STAT__)

   call MAPL_GetPointer(export, ptr2d, 'CMD_N_GAS', __RC__)
   if (associated(ptr2d)) ptr2d = cmd_N

   deallocate(cmd_N,   __STAT__)


   call MAPL_GetPointer(export, ptr2d, 'CMD_SOAG', __RC__)
   if (associated(ptr2d)) ptr2d = ((mw_SOAg/mw_air) / g_earth) * sum(q_SOAG*delp, dim=3)

   end if UPDATE_CHEM_DIAGNOSTICS

!  All done
!  --------
   call MAPL_TimerOff(mgState, 'RUN',   __RC__)
   call MAPL_TimerOff(mgState, 'TOTAL', __RC__)

   RETURN_(ESMF_SUCCESS)

   end subroutine Run_

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Finalize_ --- Finalize AChem
!
! !INTERFACE:
!

   subroutine Finalize_ ( GC, IMPORT, EXPORT, CLOCK, rc )

! !USES:

   implicit none

! !INPUT PARAMETERS:

   type(ESMF_Clock),  intent(inout)   :: CLOCK  ! The clock

! !OUTPUT PARAMETERS:

   type(ESMF_GridComp), intent(inout) :: gc     ! Grid Component
   type(ESMF_State), intent(inout)    :: IMPORT ! Import State
   type(ESMF_State), intent(inout)    :: EXPORT ! Export State
   integer, intent(out)               :: rc     ! Error return code:
                                                !  0 - all is well
                                                !  1 -

! !DESCRIPTION: This is a simple ESMF wrapper.
!
! !REVISION HISTORY:
!
!  08Aug2012  Darmenov  First crack.
!
!EOP
!-------------------------------------------------------------------------

                              __Iam__('Finalize_')

    type(AChem_State), pointer     :: self        ! Legacy state
    type(ESMF_Grid)                :: GRID        ! Grid
    type(ESMF_Config)              :: CF          ! Universal Config

    integer                        :: i1, i2, im  ! 3D Dimensions
    integer                        :: j1, j2, jm  !
    integer                        :: km          !

    integer                        :: nymd, nhms  ! date, time
    real                           :: cdt         ! time step in secs

    character(len=ESMF_MAXSTR)     :: COMP_NAME

!  Get my name and set-up traceback handle
!  ---------------------------------------
   call ESMF_GridCompGet( GC, name=comp_name, __RC__ )
   Iam = trim(comp_name) // '::' // trim(Iam)

!  Finalize MAPL Generic
!  ---------------------
    call MAPL_GenericFinalize ( GC, IMPORT, EXPORT, CLOCK,  __RC__ )

!  Extract relevant runtime information
!  ------------------------------------
   call extract_(GC, CLOCK, self, GRID, CF, i1, i2, im, j1, j2, jm, km, nymd, nhms, cdt, __RC__)


!  Free memory used to hold point volcanic emissions
!  -------------------------------------------------
   if (associated(self%volc_lat))    deallocate(self%volc_lat,   __STAT__)
   if (associated(self%volc_lon))    deallocate(self%volc_lon,   __STAT__)
   if (associated(self%volc_SO2))    deallocate(self%volc_SO2,   __STAT__)
   if (associated(self%volc_elev))   deallocate(self%volc_elev,  __STAT__)
   if (associated(self%volc_cloud))  deallocate(self%volc_cloud, __STAT__)
   if (associated(self%volc_i))      deallocate(self%volc_i,     __STAT__)
   if (associated(self%volc_j))      deallocate(self%volc_j,     __STAT__)


!  Free memory used to hold the internal copy of H2O2
!  --------------------------------------------------
   if (associated(self%h2o2))        deallocate(self%h2o2,       __STAT__)


!  Free memory for the internal the private state
!  ----------------------------------------------
   deallocate(self, __STAT__)

!  All done
!  --------
   RETURN_(ESMF_SUCCESS)

 end SUBROUTINE Finalize_

!.......................................................................

 subroutine extract_(GC, CLOCK, &
                         myState, GRID, CF, &
                         i1, i2, im,        &
                         j1, j2, jm,        &
                         km,                &
                         nymd, nhms,        &
                         cdt, rc)

    type(ESMF_GridComp), intent(INout)  :: GC           ! Grid Comp object
    type(ESMF_Clock), intent(in)        :: CLOCK        ! Clock

    type(AChem_State), pointer          :: myState      ! Legacy state
    type(ESMF_Grid),     intent(out)    :: GRID         ! Grid
    type(ESMF_Config),   intent(out)    :: CF           ! Universal Config

    integer, intent(out)                :: i1, i2, im   ! Dist grid indices
    integer, intent(out)                :: j1, j2, jm   !
    integer, intent(out)                :: km           !

    integer, intent(out)                :: nymd, nhms   ! date, time
    real, intent(out)                   :: cdt          ! time step in secs
    integer, intent(out), optional      :: rc

!                            ---

    __Iam__('extract_')

    character(len=ESMF_MAXSTR)          :: comp_name

    type(MAPL_MetaComp), pointer        :: mgState      ! MAPL generic state
    type(AChem_Wrap)                    :: wrap

    integer, dimension(3)               :: dims

    type(ESMF_Alarm)                    :: run_alarm
    type(ESMF_TimeInterval)             :: ring_interval
    real(ESMF_KIND_R8)                  :: time_step

    type(ESMF_Time)                     :: time
    integer                             :: iyr, imm, idd, ihr, imn, isc

    integer                             :: lm


!   Get my name and set-up traceback handle
!   ---------------------------------------
    call ESMF_GridCompGet(GC, name=comp_name, __RC__)
    Iam = trim(comp_name) // '::' // trim(Iam)

    rc = 0

!   Get my internal MAPL_Generic state
!   -----------------------------------
    call MAPL_GetObjectFromGC(GC, mgState, __RC__)

!   Get my internal state
!   ---------------------
    call ESMF_UserCompGetInternalState(GC, 'AChem_State', wrap, STATUS)
    VERIFY_(STATUS)
    myState => wrap%ptr

!   Get the configuration
!   ---------------------
    call ESMF_GridCompGet(GC, config=CF, __RC__)

!   Get time step
!   -------------
    call MAPL_Get(mgState, RunAlarm=run_alarm, __RC__)
    call ESMF_AlarmGet(run_alarm, ringInterval=ring_interval, __RC__)

    call ESMF_TimeIntervalGet(ring_interval, s_r8=time_step, __RC__)
    cdt = real(time_step)

!   Extract time as simple integers from clock
!   ------------------------------------------
    call ESMF_ClockGet(CLOCK, currTime=time, __RC__)
    call ESMF_TimeGet(TIME, yy=iyr, mm=imm, dd=idd, h=ihr,   m=imn,  s=isc, __RC__)

    call MAPL_PackTime(nymd, iyr, imm, idd)
    call MAPL_PackTime(nhms, ihr, imn, isc)

!   Extract the ESMF Grid
!   ---------------------
    call ESMF_GridCompGet(GC, grid=GRID, __RC__)

!   Local dimensions
!   ----------------
    call MAPL_GridGet ( grid, globalCellCountPerDim=DIMS, __RC__)

    im = dims(1)
    jm = dims(2)
    lm = dims(3)

    call ESMF_GridGet(GRID, localDE=0, &
                            staggerloc=ESMF_STAGGERLOC_CENTER, &
                            computationalCount=dims, __RC__)
    i1 = 1
    j1 = 1
    i2 = dims(1)
    j2 = dims(2)
    km = dims(3)


    RETURN_(ESMF_SUCCESS)

 end subroutine extract_


!-------------------------------------------------------------------------
!     NASA/GSFC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  GetVolcDailyTables - Get pointwise SO2 and altitude of volcanoes
!                                  from a daily file database

!
! !INTERFACE:
!

   subroutine GetVolcDailyTables(nymd, volc_emiss_file,   &
                                       nVolcPts,          &
                                       vLat, vLon,        &
                                       vElev,             &
                                       vCloud,            &
                                       vSO2,              &
                                       vStart, vEnd,      &
                                       rc)

! !USES:

   implicit none


!  Data for volcanic emissions comes from the daily inventory of all
!  volcanos (as represented by the text tables).  We return all the
!  volcanic emissions (as points, per volcano).

   integer, intent(in)            :: nymd
   character(len=*), intent(in)   :: volc_emiss_file
   integer, intent(out)           :: nVolcPts
   real, pointer, dimension(:)    :: vLat, vLon, vElev, vCloud, vSO2
   integer, pointer, dimension(:) :: vStart, vEnd

   integer, intent(out)           :: rc

   ! local
   integer                        :: i, j
   integer                        :: nLines, nCols
   integer                        :: nymd_volc, nhms_volc
   character(len=1024)            :: fname
   type(ESMF_Config)              :: cf
   real, pointer, dimension(:)    :: vData

                     __Iam__('GetVolcDailyTables')

   STATUS = ESMF_SUCCESS

!  If previous instance of volcano point data tables exist, deallocate it
!  to get the correct number of elements
   if (associated(vLat))   deallocate(vLat,   __STAT__)
   if (associated(vLon))   deallocate(vLon,   __STAT__)
   if (associated(vSO2))   deallocate(vSO2,   __STAT__)
   if (associated(vElev))  deallocate(vElev,  __STAT__)
   if (associated(vCloud)) deallocate(vCloud, __STAT__)
   if (associated(vStart)) deallocate(vStart, __STAT__)
   if (associated(vEnd))   deallocate(vEnd,   __STAT__)

   nVolcPts = 0

!  Daily files (e.g., from AEROCOM)
!  --------------------------------
!  Note: Volcanic emissions in these files are in mass of sulfur
!        Returned volcanic emissions (vSO2) are in mass of sulfur dioxide

   nymd_volc = nymd
   nhms_volc = 120000

   call StrTemplate(fname, trim(volc_emiss_file), xid='unknown', nymd=nymd_volc, nhms=nhms_volc)

   cf = ESMF_ConfigCreate()
   call ESMF_ConfigLoadFile(cf, fileName=trim(fname), __RC__)
   call ESMF_ConfigGetDim(cf, nLines, nCols, LABEL='volcano::', __RC__)

   nVolcPts = nLines

   PARSE_DATA: if (nVolcPts > 0) then
       call ESMF_ConfigFindLabel(cf, 'volcano::', __RC__)

       allocate(vData(nCols),   __STAT__)

       allocate(vLat(nLines),   __STAT__)
       allocate(vLon(nLines),   __STAT__)
       allocate(vSO2(nLines),   __STAT__)
       allocate(vElev(nLines),  __STAT__)
       allocate(vStart(nLines), __STAT__)
       allocate(vEnd(nLines),   __STAT__)
       allocate(vCloud(nLines), __STAT__)

       vStart = -1
       vEnd   = -1

       do i = 1, nLines
           call ESMF_ConfigNextLine(cf, __RC__)

           do j = 1, nCols
               call ESMF_ConfigGetAttribute(cf, vData(j), default=-1.0, __RC__)
           end do

           vLat(i)    = vData(1)
           vLon(i)    = vData(2)
           vSO2(i)    = vData(3) * mw_SO2 / mw_S
           vElev(i)   = vData(4)
           vCloud(i)  = vData(5)

           if(nCols >= 6) vStart(i)  = vData(6)
           if(nCols >= 7) vEnd(i)    = vData(7)
       end do

       where(vStart < 0) vStart = 000000
       where(vEnd   < 0) vEnd   = 240000
   end if PARSE_DATA

   call ESMF_ConfigDestroy(cf, __RC__)

   deallocate(vData, __STAT__)

   RETURN_(ESMF_SUCCESS)

 end subroutine GetVolcDailyTables


 subroutine gas_chemistry(delp,          &
                          temperature,   &
                          density_air,   &
                          q_NH3,         &
                          q_DMS,         &
                          q_MSA,         &
                          q_SO2,         &
                          q_H2SO4,       &
                          q_OH,          &
                          q_NO3,         &
                          cpl_NH3,       &
                          cpl_DMS,       &
                          cpl_MSA,       &
                          cpl_SO2,       &
                          cpl_H2SO4,     &
                          dt,            &
                          rc)

! !USES:

   use kpp_achem_gas_Precision,  only: kpp_r8 => dp

   use kpp_achem_gas_Global,     only: kpp_uf_conc       => CFACTOR, &
                                       kpp_conc          => C,       &
                                       kpp_reaction_rate => RCONST,  &
                                       kpp_sun           => SUN,     &
                                       kpp_temperature   => TEMP,    &
                                       kpp_conc_air      => c_air,   &
                                       kpp_conc_O2       => c_O2,    &
                                       kpp_time          => TIME,    &
                                       kpp_time_start    => TSTART,  &
                                       kpp_time_end      => TEND,    &
                                       kpp_dt            => DT,      &
                                       kpp_step_min      => STEPMIN, &
                                       kpp_step_max      => STEPMAX, &
                                       kpp_rtol          => RTOL,    &
                                       kpp_atol          => ATOL

   use kpp_achem_gas_Parameters, only: kpp_iNH3   => ind_NH3,        &
                                       kpp_iDMS   => ind_DMS,        &
                                       kpp_iSO2   => ind_SO2,        &
                                       kpp_iH2SO4 => ind_H2SO4,      &
                                       kpp_iMSA   => ind_MSA,        &
                                       kpp_iNO3   => ind_NO3,        &
                                       kpp_iOH    => ind_OH

   use kpp_achem_gas_Integrator, only: kpp_integrate => INTEGRATE

   use kpp_achem_gas_Rates,      only: kpp_update_sun => Update_SUN, &
                                       kpp_update_reaction_rates => Update_RCONST

   implicit none

! !INPUT PARAMETERS:

   real, dimension(:,:,:), intent(in)    :: delp
   real, dimension(:,:,:), intent(in)    :: temperature
   real, dimension(:,:,:), intent(in)    :: density_air

   real, dimension(:,:,:), intent(in)    :: q_OH
   real, dimension(:,:,:), intent(in)    :: q_NO3

   real, intent(in)                      :: dt

   integer, intent(out)                  :: rc

! !OUTPUT PARAMETERS:

   real, dimension(:,:,:), intent(inout) :: q_NH3
   real, dimension(:,:,:), intent(inout) :: q_DMS
   real, dimension(:,:,:), intent(inout) :: q_MSA
   real, dimension(:,:,:), intent(inout) :: q_SO2
   real, dimension(:,:,:), intent(inout) :: q_H2SO4

   real, dimension(:,:),   intent(inout) :: cpl_NH3
   real, dimension(:,:),   intent(inout) :: cpl_DMS
   real, dimension(:,:),   intent(inout) :: cpl_MSA
   real, dimension(:,:),   intent(inout) :: cpl_SO2
   real, dimension(:,:),   intent(inout) :: cpl_H2SO4



! !DESCRIPTION: Wrap the KPP generated code.
!
! !REVISION HISTORY:
!
!  13Aug2012 A. Darmenov   First crack.
!
!EOP
!-------------------------------------------------------------------------

                              __Iam__('gas_chemistry')

   ! parameters

   integer, parameter       :: r8 = kpp_r8

   real(kind=r8), parameter :: zero_concentration = 1d-6   ! very small concentration, # cm-3

   real(kind=r8), parameter :: rel_tolerance = 1.0d-2
   real(kind=r8), parameter :: abs_tolerance = 1.0d-2

   integer,       parameter :: ICNTRL_U(20) = (/1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/)
   real(kind=r8), parameter :: RCNTRL_U(20) = (/0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, &
                                                0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, &
                                                0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, &
                                                0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)

   ! local
   integer       :: ierr
   real(kind=r8) :: RSTATE(20)
   real(kind=r8) :: conc_air_inv
   real(kind=r8) :: time

   integer :: i, i1, i2
   integer :: j, j1, j2
   integer :: k, k1, km
   integer :: ijl, ijkl
   integer :: n

   real    :: ff


#ifdef DEBUG
   real, dimension(:,:,:), allocatable :: conc_air
#endif

   rc = 0

   i1 = lbound(density_air, dim=1); i2 = ubound(density_air, dim=1)
   j1 = lbound(density_air, dim=2); j2 = ubound(density_air, dim=2)
   k1 = lbound(density_air, dim=3); km = ubound(density_air, dim=3)

   ijl  = (i2 - i1 + 1) * (j2 - j1 + 1)
   ijkl = ijl * (km - k1 + 1)

#ifdef DEBUG
   ! air molecules concentrations, #molecules/cm-3
   allocate(conc_air(i1:i2,j1:j2,k1:km), __STAT__)
   conc_air = 1.0e-6 * (N_avog/mw_air) * density_air


   call write_parallel('[ DEBUG ]   ' // trim(Iam) // ': inputs to gas phase chemistry')

   call write_parallel('[ DEBUG ]   ' // 'volume mixing ratios:')
   call MAPL_MaxMin('OH    : ', q_OH(:,:,km))
   call MAPL_MaxMin('NO3   : ', q_NO3(:,:,km))
   call MAPL_MaxMin('DMS   : ', q_DMS(:,:,km))
   call MAPL_MaxMin('MSA   : ', q_MSA(:,:,km))
   call MAPL_MaxMin('SO2   : ', q_SO2(:,:,km))
   call MAPL_MaxMin('H2SO4 : ', q_H2SO4(:,:,km))
   call MAPL_MaxMin('NH3   : ', q_NH3(:,:,km))

   call write_parallel('[ DEBUG ]   ' // 'concentrations:')
   call MAPL_MaxMin('air   : ', conc_air(:,:,km))
   call MAPL_MaxMin('OH    : ', q_OH(:,:,km)    * conc_air(:,:,km))
   call MAPL_MaxMin('NO3   : ', q_NO3(:,:,km)   * conc_air(:,:,km))
   call MAPL_MaxMin('DMS   : ', q_DMS(:,:,km)   * conc_air(:,:,km))
   call MAPL_MaxMin('MSA   : ', q_MSA(:,:,km)   * conc_air(:,:,km))
   call MAPL_MaxMin('SO2   : ', q_SO2(:,:,km)   * conc_air(:,:,km))
   call MAPL_MaxMin('H2SO4 : ', q_H2SO4(:,:,km) * conc_air(:,:,km))
   call MAPL_MaxMin('NH3   : ', q_NH3(:,:,km)   * conc_air(:,:,km))
#endif


   cpl_NH3   = sum(q_NH3  * delp, dim=3)
   cpl_DMS   = sum(q_DMS  * delp, dim=3)
   cpl_MSA   = sum(q_MSA  * delp, dim=3)
   cpl_SO2   = sum(q_SO2  * delp, dim=3)
   cpl_H2SO4 = sum(q_H2SO4* delp, dim=3)

   do k = k1, km
       do j = j1, j2
           do i = i1, i2
               ! set the onversion factor for concentration units
               kpp_uf_conc = 1.0_r8

               ! set sunlight intensity for KPP
               kpp_sun = 1.0_r8

               ! initialize the KPP integration parameteres
               kpp_step_min = 0.0  ! seconds
               kpp_step_max = 0.0

               ! set the KPP realtive (RTOL) and absolute (ATOL) tolerances
               kpp_rtol(:) = rel_tolerance
               kpp_atol(:) = abs_tolerance

               ! set KPP times
               kpp_dt         = dt
               kpp_time_start = 0.0_r8
               kpp_time_end   = kpp_time_start + dt

               ! set temperature for KPP
               kpp_temperature      = temperature(i,j,k)

               ! set concentrations for KPP, #molecules/cm-3
               kpp_conc_air         = 1.0e-6 * (N_avog/mw_air) * density_air(i,j,k)
               kpp_conc_O2          = 0.20946 * kpp_conc_air

               kpp_conc(kpp_iNO3)   = kpp_conc_air * q_NO3(i,j,k)
               kpp_conc(kpp_iOH)    = kpp_conc_air * q_OH(i,j,k)

               kpp_conc(kpp_iDMS)   = kpp_conc_air * q_DMS(i,j,k)
               kpp_conc(kpp_iSO2)   = kpp_conc_air * q_SO2(i,j,k)
               kpp_conc(kpp_iH2SO4) = kpp_conc_air * q_H2SO4(i,j,k)
               kpp_conc(kpp_iMSA)   = kpp_conc_air * q_MSA(i,j,k)
               kpp_conc(kpp_iNH3)   = kpp_conc_air * q_NH3(i,j,k)

               time = kpp_time_start
               KPP_TIME_INTEGRATE: do while (time < kpp_time_end)

                   kpp_time = time

                   ! chemistry solver
                   call kpp_update_reaction_rates()

                   ! set sunlight intensity for KPP
                   kpp_sun = 1.0_r8

                   kpp_conc(kpp_iNO3)   = kpp_conc_air * q_NO3(i,j,k)
                   kpp_conc(kpp_iOH)    = kpp_conc_air * q_OH(i,j,k)


                   call kpp_integrate(TIN       = time,           &
                                      TOUT      = time+kpp_dt,    &
                                      RSTATUS_U = RSTATE,         &
                                      ICNTRL_U  = ICNTRL_U,       &
                                      RCNTRL_U  = RCNTRL_U,       &
                                      IERR_U    = ierr)

                   if (ierr < 0) then
                       STATUS = -1
                   else
                       STATUS = ESMF_SUCCESS
                   end if

                   VERIFY_(STATUS)

                   time = RSTATE(1)
               end do KPP_TIME_INTEGRATE



               ! update the model concentrations
               if (kpp_conc_air > zero_concentration) then
                   conc_air_inv = 1.0_r8 / kpp_conc_air
               else
                   conc_air_inv = 0.0_r8
               end if

               q_DMS(i,j,k)   = kpp_conc(kpp_iDMS)   * conc_air_inv
               q_SO2(i,j,k)   = kpp_conc(kpp_iSO2)   * conc_air_inv
               q_H2SO4(i,j,k) = kpp_conc(kpp_iH2SO4) * conc_air_inv
               q_MSA(i,j,k)   = kpp_conc(kpp_iMSA)   * conc_air_inv
               q_NH3(i,j,k)   = kpp_conc(kpp_iNH3)   * conc_air_inv
           end do
       end do
   end do

   ff = 1.0 / (mw_air * g_earth * dt)

   cpl_NH3   = (sum(q_NH3   * delp, dim=3) - cpl_NH3  ) * ff
   cpl_DMS   = (sum(q_DMS   * delp, dim=3) - cpl_DMS  ) * ff
   cpl_MSA   = (sum(q_MSA   * delp, dim=3) - cpl_MSA  ) * ff
   cpl_SO2   = (sum(q_SO2   * delp, dim=3) - cpl_SO2  ) * ff
   cpl_H2SO4 = (sum(q_H2SO4 * delp, dim=3) - cpl_H2SO4) * ff


#ifdef DEBUG
   call write_parallel('[ DEBUG ]   ' // trim(Iam) // ': fields after gas phase chemistry')

   call write_parallel('[ DEBUG ]   ' // 'volume mixing ratios:')
   call MAPL_MaxMin('OH    : ', q_OH(:,:,km))
   call MAPL_MaxMin('NO3   : ', q_NO3(:,:,km))
   call MAPL_MaxMin('DMS   : ', q_DMS(:,:,km))
   call MAPL_MaxMin('MSA   : ', q_MSA(:,:,km))
   call MAPL_MaxMin('SO2   : ', q_SO2(:,:,km))
   call MAPL_MaxMin('H2SO4 : ', q_H2SO4(:,:,km))
   call MAPL_MaxMin('NH3   : ', q_NH3(:,:,km))

   call write_parallel('[ DEBUG ]   ' // 'concentrations:')
   call MAPL_MaxMin('air   : ', conc_air(:,:,km))
   call MAPL_MaxMin('OH    : ', q_OH(:,:,km)    * conc_air(:,:,km))
   call MAPL_MaxMin('NO3   : ', q_NO3(:,:,km)   * conc_air(:,:,km))
   call MAPL_MaxMin('DMS   : ', q_DMS(:,:,km)   * conc_air(:,:,km))
   call MAPL_MaxMin('MSA   : ', q_MSA(:,:,km)   * conc_air(:,:,km))
   call MAPL_MaxMin('SO2   : ', q_SO2(:,:,km)   * conc_air(:,:,km))
   call MAPL_MaxMin('H2SO4 : ', q_H2SO4(:,:,km) * conc_air(:,:,km))
   call MAPL_MaxMin('NH3   : ', q_NH3(:,:,km)   * conc_air(:,:,km))

   deallocate(conc_air, __STAT__)
#endif

   where (q_DMS   < 0.0) q_DMS   = tiny(0.0)
   where (q_MSA   < 0.0) q_MSA   = tiny(0.0)
   where (q_SO2   < 0.0) q_SO2   = tiny(0.0)
   where (q_H2SO4 < 0.0) q_H2SO4 = tiny(0.0)
   where (q_NH3   < 0.0) q_NH3   = tiny(0.0)

   RETURN_(ESMF_SUCCESS)

 end subroutine gas_chemistry


 subroutine aqu_chemistry_fast(ple,           &
                               temperature,   &
                               density_air,   &
                               lwc,           &
                               lwc_min,       &
                               fcld,          &
                               q_NH3,         &
                               q_SO2,         &
                               q_H2SO4,       &
                               q_H2O2,        &
                               pSO4_aq,       &
                               pNH4_aq,       &
                               pSO4_aq_SO2,   &
                               pSO4_aq_H2SO4, &
                               pNH4_aq_NH3,   &
                               dt,            &
                               rc)

! !USES:

   implicit none

! !INPUT PARAMETERS:

   real, dimension(:,:,:), intent(in)    :: ple
   real, dimension(:,:,:), intent(in)    :: temperature
   real, dimension(:,:,:), intent(in)    :: density_air
   real, dimension(:,:,:), intent(in)    :: lwc
   real, intent(in)                      :: lwc_min
   real, dimension(:,:,:), intent(in)    :: fcld

   real, intent(in)                      :: dt

   integer, intent(out)                  :: rc

! !OUTPUT PARAMETERS:
   real, dimension(:,:,:), intent(inout) :: q_H2O2

   real, dimension(:,:,:), intent(inout) :: q_NH3
   real, dimension(:,:,:), intent(inout) :: q_SO2
   real, dimension(:,:,:), intent(inout) :: q_H2SO4

   real, dimension(:,:,:), intent(out)   :: pSO4_aq
   real, dimension(:,:,:), intent(out)   :: pNH4_aq

   real, dimension(:,:,:), intent(out)   :: pSO4_aq_SO2
   real, dimension(:,:,:), intent(out)   :: pSO4_aq_H2SO4
   real, dimension(:,:,:), intent(out)   :: pNH4_aq_NH3



! !DESCRIPTION: Super fast implementation of aqueous chemistry.
!
! !REVISION HISTORY:
!
!  09Nov2012 A. Darmenov   First crack.
!
!EOP
!-------------------------------------------------------------------------

                              __Iam__('aqu_chemistry_fast')

   ! parameters
   real, parameter :: T_freeze = 258.0            ! freezing point of supercooled cloud water, K


   ! local
   integer :: i, i1, i2
   integer :: j, j1, j2
   integer :: k, k1, km
   integer :: ijl, ijkl

   real    :: SO2, H2SO4, H2O2, NH3
   real    :: l_SO2, l_H2SO4, l_H2O2, l_NH3
   real    :: f


   rc = ESMF_SUCCESS

   i1 = lbound(q_SO2, dim=1); i2 = ubound(q_SO2, dim=1)
   j1 = lbound(q_SO2, dim=2); j2 = ubound(q_SO2, dim=2)
   k1 = lbound(q_SO2, dim=3); km = ubound(q_SO2, dim=3)

   ijl  = (i2 - i1 + 1) * (j2 - j1 + 1)
   ijkl = ijl * (km - k1 + 1)

   ! total production in aqueous phase
   pSO4_aq = 0.0
   pNH4_aq = 0.0

   ! contributions of production pathways in aqueous phase
   pSO4_aq_SO2   = 0.0
   pSO4_aq_H2SO4 = 0.0
   pNH4_aq_NH3   = 0.0

   do k = k1, km
       do j = j1, j2
           do i = i1, i2

               SO2   = q_SO2(i,j,k)
               H2SO4 = q_H2SO4(i,j,k)
               H2O2  = q_H2O2(i,j,k)
               NH3   = q_NH3(i,j,k)

               if ((lwc(i,j,k) > lwc_min) .and. (temperature(i,j,k) > T_freeze)) then

                   f = fcld(i,j,k)

                   if (SO2 > H2O2) then
                       f = f * (H2O2 / SO2)
                       H2O2 = H2O2 * (1 - fcld(i,j,k))
                   else
                       H2O2 = H2O2 * (1 - fcld(i,j,k)*(SO2/H2O2))
                   endif

                   ! aqueous loss
                   l_SO2   = f * SO2
                   l_H2SO4 = fcld(i,j,k) * H2SO4
                   l_NH3   = fcld(i,j,k) * NH3     ! all NH3 dissociates to NH4 i.e., [NH4] = [NH3] for pH < 5

                   ! update TMR
                   SO2   = SO2   * (1 - f)
                   H2SO4 = H2SO4 * (1 - fcld(i,j,k))
                   NH3   = NH3   * (1 - fcld(i,j,k))
               else
                   l_SO2   = 0.0
                   l_H2SO4 = 0.0
                   l_NH3   = 0.0
               endif

               ! H2O2 mixing ratio should be updated at this point
               ! and then reset it periodically
               SO2   = max(SO2,   tiny(SO2))
               H2SO4 = max(H2SO4, tiny(H2SO4))
               H2O2  = max(H2O2,  tiny(H2O2))
               NH3   = max(NH3,   tiny(NH3))

               q_SO2(i,j,k)   = SO2
               q_H2SO4(i,j,k) = H2SO4
               q_H2O2(i,j,k)  = H2O2
               q_NH3(i,j,k)   = NH3

               ! units are 'kg-SO4/kg-air/s'
               pSO4_aq_SO2(i,j,k)   = (mw_SO4/mw_air) * l_SO2   / dt
               pSO4_aq_H2SO4(i,j,k) = (mw_SO4/mw_air) * l_H2SO4 / dt

               pSO4_aq(i,j,k) = pSO4_aq(i,j,k) + pSO4_aq_SO2(i,j,k)
               pSO4_aq(i,j,k) = pSO4_aq(i,j,k) + pSO4_aq_H2SO4(i,j,k)

               ! units are 'kg-NH4/kg-air/s'
               pNH4_aq_NH3(i,j,k) = (mw_NH4/mw_air) * l_NH3 / dt
               pNH4_aq(i,j,k)     = pNH4_aq(i,j,k) + pNH4_aq_NH3(i,j,k)
           end do
       end do
   end do

   RETURN_(ESMF_SUCCESS)

 end subroutine aqu_chemistry_fast


 subroutine aqu_chemistry(ple,           &
                          temperature,   &
                          density_air,   &
                          lwc,           &
                          fcld,          &
                          q_NH3,         &
                          q_SO2,         &
                          q_SVI_aq,      &
                          q_NH4_aq,      &
                          q_H2O2,        &
                          q_O3,          &
                          dt,            &
                          solver_max_dt, &
                          lwc_min,       &
                          rc)

! !USES:

   implicit none

! !INPUT PARAMETERS:

   real, dimension(:,:,:), intent(in)    :: ple
   real, dimension(:,:,:), intent(in)    :: temperature
   real, dimension(:,:,:), intent(in)    :: density_air
   real, dimension(:,:,:), intent(in)    :: lwc
   real, dimension(:,:,:), intent(in)    :: fcld

   real, dimension(:,:,:), intent(in)    :: q_H2O2
   real, dimension(:,:,:), intent(in)    :: q_O3       ! has to be converted to VMR

   real, intent(in)                      :: dt
   real, intent(in)                      :: solver_max_dt
   real, intent(in)                      :: lwc_min

   integer, intent(out)                  :: rc

! !OUTPUT PARAMETERS:

   real, dimension(:,:,:), intent(inout) :: q_NH3
   real, dimension(:,:,:), intent(inout) :: q_SO2

   real, dimension(:,:,:), intent(inout) :: q_SVI_aq
   real, dimension(:,:,:), intent(inout) :: q_NH4_aq



! !DESCRIPTION: Fast implementation of aqueous chemistry.
!
! !REVISION HISTORY:
!
!  09Nov2012 A. Darmenov   First crack.
!
!EOP
!-------------------------------------------------------------------------

                              __Iam__('aqu_chemistry')

   ! parameters
   real, parameter :: atm = 1.01325e5             ! one atmosphere (pressure), Pa

   real, parameter :: zero_concentration = 1.0e-6 ! very small concentration (one molecule per cm3), # cm-3

   real, parameter :: Hp_0     = 10**-5.6
   real, parameter :: Hp_cloud = 10**-4.5
   real, parameter :: f_Hp     = 0.1


   ! local
   real    :: conc_air_inv

   integer :: i, i1, i2
   integer :: j, j1, j2
   integer :: k, k1, km
   integer :: ijl, ijkl

   integer :: n, n_steps                          ! number of time splitting sub-steps
   real    :: dt_step                             ! length of a time sub-step

   real    :: T                                   ! temperature
   real    :: Hp                                  ! hydrogen ion [H^+] concentration
   real    :: c_air, c_air_inv                    ! concentration of air and its reciprocal
   real    :: g_NH3, g_SO2, g_H2O2, g_O3          ! gas phase concentrations
   real    :: a_NH3, a_NH4p, a_SO2,             & ! aqueous phase concentrations
              a_OHm, a_H2O2, a_O3,              & ! ...
              a_HSO3m, a_SO3mm, a_HSO4m, a_SVI    ! ...
   real    :: p_NH3, p_SO2, p_H2O2, p_O3          ! partial pressure
   real    :: H_NH3, H_SO2, H_H2O2, H_O3          ! Henry's law constant

   real    :: f_pp                                ! partial pressure factor

   real    :: Kw, Ka1, Ks1, Ks2, Kso4             ! equilibrium constants for dissociation reactions

   real    :: k_SIV_O3, k_SIV_H2O2                ! S(IV)-O3 and S(IV)-H2O2 aqueous phase oxidation rates
   real    :: delta_a_SIV                         ! loss/production of S(IV)/S(VI)

   rc = 0

   i1 = lbound(density_air, dim=1); i2 = ubound(density_air, dim=1)
   j1 = lbound(density_air, dim=2); j2 = ubound(density_air, dim=2)
   k1 = lbound(density_air, dim=3); km = ubound(density_air, dim=3)

   ijl  = (i2 - i1 + 1) * (j2 - j1 + 1)
   ijkl = ijl * (km - k1 + 1)

   n_steps = max(1, nint(dt / solver_max_dt))
   dt_step = dt / n_steps

   do k = k1, km
       do j = j1, j2
           do i = i1, i2
               if (lwc(i,j,k) > lwc_min) then

                   T = temperature(i,j,k)

                   ! air molecules concentrations, #molecules/cm-3
                   c_air = 1.0e-6 * (N_avog/mw_air) * density_air(i,j,k)

                   if (c_air > zero_concentration) then
                       c_air_inv = 1.0 / c_air
                   else
                       c_air_inv = 0.0
                   end if

                   ! compute Henry's constants
                   H_NH3  = henry(H_NH3_298,  298.0, E_R_NH3,  T)
                   H_SO2  = henry(H_SO2_298,  298.0, E_R_SO2,  T)
                   H_H2O2 = henry(H_H2O2_298, 298.0, E_R_H2O2, T)
                   H_O3   = henry(H_O3_298,   298.0, E_R_O3,   T)

                   ! gas-phase concentrations and partial pressures
                   f_pp   = (R_univ / mw_air) * (density_air(i,j,k) * T) / atm

                   g_O3   = c_air * q_O3(i,j,k) * (MAPL_AIRMW / MAPL_O3MW)
                   p_O3   = f_pp  * q_O3(i,j,k) * (MAPL_AIRMW / MAPL_O3MW)

                   g_H2O2 = c_air * q_H2O2(i,j,k)
                   p_H2O2 = f_pp  * q_H2O2(i,j,k)

                   g_SO2  = c_air * q_SO2(i,j,k)
                   p_SO2  = f_pp  * q_SO2(i,j,k)

                   g_NH3  = c_air * q_NH3(i,j,k)
                   p_NH3  = f_pp  * q_NH3(i,j,k)

                   a_SVI  = 1e3 * q_SVI_aq(i,j,k) / (mw_air * lwc(i,j,k)) ! convert from mol/mol-air to mol L-1

                   ! dissociation  rates
                   Kw   = K_w(T)
                   Ka1  = K_a1(T)
                   Ks1  = K_s1(T)
                   Ks2  = K_s2(T)
                   Kso4 = K_1696(T)


                   ! time sub-splitting
                   do n = 1, n_steps
                       ! compute [H+]
#if (0)
                       Hp = Hp_0 + f_Hp * (2*c_SO4mm + c_HSO3m - c_NH4m) ! ~~ (c_SO4mm + c_HSO3m) ~~ (c_SO4mm + SO2_aq)
#else
                       Hp = Hp_cloud ! fixed pH
#endif
                       ! compute equilibrium aqueous phase concentrations, mol L-1
                       a_OHm   = Kw / Hp

                       a_NH4p  = H_NH3 * p_NH3 * (Ka1 / Kw) * Hp

                       a_HSO3m = H_SO2 * p_SO2 * Ks1 / Hp
                       a_SO3mm = H_SO2 * p_SO2 * Ks1 * Ks2 / Hp**2

                       a_HSO4m = Hp * a_SVI / (Hp + Kso4)

                       ! chemical reaction rates
                       k_SIV_H2O2 = 0.0
                       k_SIV_O3   = 0.0
                       k_SIV_H2O2 = 0.0

                       ! integrate
                       delta_a_SIV = 0.0 * dt_step

                       delta_a_SIV = 0.0 * dt_step


                       ! update model state
                       !q_SO2(i,j,k) = c_SO2 * conc_air_inv
                   end do

                   ! update model state
                   q_NH4_aq(i,j,k) = q_NH4_aq(i,j,k) + (lwc(i,j,k) * mw_air * a_NH4p)
               end if

           end do
       end do
   end do


   where (q_SO2   < 0.0) q_SO2   = tiny(0.0)
   where (q_NH3   < 0.0) q_NH3   = tiny(0.0)

   RETURN_(ESMF_SUCCESS)

   contains

   elemental real function K_w(T)
       implicit none
       real, intent(in) :: T

       real, parameter :: K_298 = 1.0e-14  ! M atm-1
       real, parameter :: dH_R  = 6710.0   ! K

       K_w = K_298 * exp(-dH_R * (1/T - 1/298.0))
   end function K_w

   elemental real function K_s1(T)
       implicit none
       real, intent(in) :: T

       real, parameter :: K_298 = 1.3e-2   ! M atm-1
       real, parameter :: dH_R  = -1960.0  ! K

       K_s1 = K_298 * exp(-dH_R * (1/T - 1/298.0))
   end function K_s1

   elemental real function K_s2(T)
       implicit none
       real, intent(in) :: T

       real, parameter :: K_298 = 6.6e-8   ! M atm-1
       real, parameter :: dH_R  = -1500.0  ! K

       K_s2 = K_298 * exp(-dH_R * (1/T - 1/298.0))
   end function K_s2

   elemental real function K_1696(T)
       implicit none
       real, intent(in) :: T

       real, parameter :: K_298 = 1.02e-2  ! M atm-1
       real, parameter :: dH_R  = -2720.0  ! K

       K_1696 = K_298 * exp(-dH_R * (1/T - 1/298.0))
   end function K_1696

   elemental real function K_a1(T)
       implicit none
       real, intent(in) :: T

       real, parameter :: K_298 = 1.7e-5   ! M atm-1
       real, parameter :: dH_R  = 450.0  ! K

       K_a1 = K_298 * exp(-dH_R * (1/T - 1/298.0))
   end function K_a1


 end subroutine aqu_chemistry


 subroutine ocs_chemistry(ple,             &
                          temperature,     &
                          density_air,     &
                          ocs_surface_vmr, &
                          tropp,           &
                          q_OCS,           &
                          q_OH,            &
                          q_O3p,           &
                          j_ocs,           &
                          pSO2_OCS,        &
                          pSO2_OCS_OH,         &
                          pSO2_OCS_O3p,        &
                          pSO2_OCS_jOCS,       &
                          lOCS,        &
                          lOCS_OH,         &
                          lOCS_O3p,        &
                          lOCS_jOCS,       &
                          dt,              &
                          rc)

! !USES:

   use kpp_achem_gas_Precision,  only: kpp_r8 => dp

   implicit none

! !INPUT PARAMETERS:

   real, dimension(:,:,:), intent(in) :: ple
   real, dimension(:,:,:), intent(in) :: temperature
   real, dimension(:,:,:), intent(in) :: density_air
   real                  , intent(in) :: ocs_surface_vmr

   real, intent(in)                   :: dt

   integer, intent(out)               :: rc

   real, dimension(:,:)  , intent(in) :: tropp

   real, dimension(:,:,:), intent(in) :: q_OH
   real, dimension(:,:,:), intent(in) :: q_O3p
   real, dimension(:,:,:), intent(in) :: j_ocs

! !OUTPUT PARAMETERS:
   real, dimension(:,:,:), intent(inout) :: q_OCS

   real, dimension(:,:,:), intent(inout) :: pSO2_OCS      ! production of S from OCS in the stratosphere
   real, dimension(:,:,:), intent(inout) :: pSO2_OCS_OH   ! production of S from OCS+OH
   real, dimension(:,:,:), intent(inout) :: pSO2_OCS_O3p  ! production of S from OCS+O3p
   real, dimension(:,:,:), intent(inout) :: pSO2_OCS_jOCS ! production of S from OCS photolysis
   real, dimension(:,:,:), intent(inout) :: lOCS          ! loss ocs, 'molecules cm-3 s-1'
   real, dimension(:,:,:), intent(inout) :: lOCS_OH       ! loss rate of OCS from OCS+OH, 'molec cm-3 s-1'
   real, dimension(:,:,:), intent(inout) :: lOCS_O3p      ! loss rate of OCS from OCS+O3p, 'molec cm-3 s-1'
   real, dimension(:,:,:), intent(inout) :: lOCS_jOCS     ! loss rate of OCS from photolysis, 'molec cm-3 s-1'


! !DESCRIPTION: Calculate the OCS chemistry
!
! !REVISION HISTORY:
!
!  16 April 2014 V. Aquila   First crack.
!
!EOP
!-------------------------------------------------------------------------

                              __Iam__('ocs_chemistry')

   ! parameters

   integer, parameter       :: r8 = kpp_r8
   real                     :: conc_air, c_air_inv         ! air concentration [molecules/cm3]
   real                     :: conc_OCS                    ! OCS concentration [molecules/cm3]
   real                     :: conc_OH                     ! OH concentration [molecules/cm3]
   real                     :: conc_O3p                    ! O3p concentration [molecules/cm3]
   real                     :: ltot_ocs                    ! OCS loss [molecules/cm3]

   real, parameter          :: mw_SO2   = 64.066           ! molar mass of sulfur dioxide, g mol-1
   real, parameter          :: mw_OH    = 17.01            ! molar mass of hydroxide ion, g mol -1
   real, parameter          :: mw_air   = 28.97            ! molar mass of dry air, g mol-1

   real, parameter          :: N_A      = 6.02214129e23    ! mol-1

   real(kind=r8), parameter :: zero_concentration = 1e-6   ! very small concentration, # cm-3

   real                     :: k_oh, k_o3p, kk              ! reaction rates

   !reactions parameters from Sander, S. P. et al. (2010), Chemical Kinetics and Photochemical
   !Data for Use in Atmospheric Studies (No. 17) NASA JPL.
   !Rate constant ak = A*exp[-ER/Temperature]
   real, parameter          :: A_oh      =  1.1e-13        ! Arrhenius A-factor OCS+OH
   real, parameter          :: ER_oh     =  1200           ! Temperature dependence OCS+OH
   real, parameter          :: A_o3p     =  2.1e-11        ! Arrhenius A-factor OCS+O3p
   real, parameter          :: ER_o3p    =  2200           ! Temperature dependence OCS+O3p

   real, allocatable, dimension(:, :,:) :: prod_SO2        ! production of SO2 from OCS
   real, allocatable, dimension(:, :,:) :: prod_SO2_OH     ! production of SO2 from OCS+OH
   real, allocatable, dimension(:, :,:) :: prod_SO2_O3p    ! production of SO2 from OCS+O3p
   real, allocatable, dimension(:, :,:) :: prod_SO2_jOCS   ! production of SO2 from OCS photolysis
   real, allocatable, dimension(:, :,:) :: loss_ocs           ! loss OCS. molecules cm-3
   real, allocatable, dimension(:, :,:) :: loss_ocs_OH        ! loss OCS. molecules cm-3
   real, allocatable, dimension(:, :,:) :: loss_ocs_JOCS      ! loss OCS. molecules cm-3
   real, allocatable, dimension(:, :,:) :: loss_ocs_O3p       ! loss OCS. molecules cm-3

   integer :: i, i1, i2
   integer :: j, j1, j2
   integer :: k, k1, km

   i1 = lbound(density_air, dim=1); i2 = ubound(density_air, dim=1)
   j1 = lbound(density_air, dim=2); j2 = ubound(density_air, dim=2)
   k1 = lbound(density_air, dim=3); km = ubound(density_air, dim=3)

   allocate(prod_SO2(i1:i2,j1:j2,k1:km), __STAT__)
   allocate(prod_SO2_OH(i1:i2,j1:j2,k1:km), __STAT__)
   allocate(prod_SO2_O3p(i1:i2,j1:j2,k1:km), __STAT__)
   allocate(prod_SO2_jOCS(i1:i2,j1:j2,k1:km), __STAT__)

   allocate(loss_ocs(i1:i2,j1:j2,k1:km), __STAT__)
   allocate(loss_ocs_OH(i1:i2,j1:j2,k1:km), __STAT__)
   allocate(loss_ocs_O3p(i1:i2,j1:j2,k1:km), __STAT__)
   allocate(loss_ocs_jOCS(i1:i2,j1:j2,k1:km), __STAT__)

   prod_SO2      = 0.0
   prod_SO2_OH   = 0.0
   prod_SO2_O3p  = 0.0
   prod_SO2_jOCS = 0.0

   loss_ocs      = 0.0
   loss_ocs_OH   = 0.0
   loss_ocs_O3p  = 0.0
   loss_ocs_jOCS = 0.0

   q_OCS(i1:i2,j1:j2,km) = ocs_surface_vmr

   do k = k1, km
       do j = j1, j2
           do i = i1, i2
               STRATOSPHERE: if (ple(i,j,k) <= tropp(i,j)) then
                   !transform from mol/mol to molecules/cm3
                   conc_air = dble(1e-3 * (N_A/mw_air) * density_air(i,j,k))

                   conc_OCS   = conc_air * dble(q_OCS(i,j,k))
                   conc_OH    = conc_air * dble(q_OH(i,j,k))
                   conc_O3p   = conc_air * dble(q_O3p(i,j,k))

                   k_oh  = A_oh  * exp(- ER_oh  / temperature(i,j,k))
                   k_o3p = A_o3p * exp(- ER_o3p / temperature(i,j,k))
                   kk = (k_oh * conc_OH + k_o3p * conc_O3p + j_ocs(i,j,k))

                   ltot_ocs = conc_OCS * (1 - exp(-kk * dt))

                   if (conc_air > zero_concentration) then
                       c_air_inv = 1.0 / conc_air
                   else
                       c_air_inv = 0.0
                   end if

                   q_OCS(i,j,k) = q_OCS(i,j,k) - ltot_ocs * c_air_inv
                   q_OCS(i,j,k) = max(q_OCS(i,j,k), tiny(0.0))

                   prod_SO2(i,j,k)      = ltot_ocs * c_air_inv * mw_SO2/mw_air / dt

                   prod_SO2_OH(i,j,k)   = (k_oh * conc_OH)   / kk * ltot_ocs * c_air_inv * mw_SO2/mw_air / dt
                   prod_SO2_O3p(i,j,k)  = (k_o3p * conc_O3p) / kk * ltot_ocs * c_air_inv * mw_SO2/mw_air / dt
                   prod_SO2_jOCS(i,j,k) = (j_OCS(i,j,k))     / kk * ltot_ocs * c_air_inv * mw_SO2/mw_air / dt

                   loss_ocs(i,j,k)      =  ltot_ocs / dt

                   loss_ocs_OH(i,j,k)   = (k_oh * conc_OH)   / kk * ltot_ocs / dt
                   loss_ocs_O3p(i,j,k)  = (k_o3p * conc_O3p) / kk * ltot_ocs / dt
                   loss_ocs_jOCS(i,j,k) = (j_OCS(i,j,k))     / kk * ltot_ocs / dt
               endif STRATOSPHERE
           enddo
        enddo
     enddo

     pSO2_OCS      = prod_SO2
     pSO2_OCS_OH   = prod_SO2_OH
     pSO2_OCS_O3p  = prod_SO2_O3p
     pSO2_OCS_jOCS = prod_SO2_jOCS

     lOCS          = loss_ocs
     lOCS_OH       = loss_ocs_OH
     lOCS_O3p      = loss_ocs_O3p
     lOCS_jOCS     = loss_ocs_jOCS


     deallocate(prod_SO2,      __STAT__)
     deallocate(prod_SO2_OH,   __STAT__)
     deallocate(prod_SO2_O3p,  __STAT__)
     deallocate(prod_SO2_jOCS, __STAT__)

     deallocate(loss_ocs,      __STAT__)
     deallocate(loss_ocs_OH,   __STAT__)
     deallocate(loss_ocs_O3p,  __STAT__)
     deallocate(loss_ocs_jOCS, __STAT__)

 end subroutine ocs_chemistry


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  solar_zenith_angle --- Given day of the year, UTC time and
!             geographical location computes solar zenith angle.
!
! !INTERFACE:
!

 subroutine solar_zenith_angle(doy, utc_hour, lon, lat, sza, cos_sza)

! !USES:

   implicit none

! !INPUT PARAMETERS:

   integer, intent(in)                 :: doy      ! day since the begining of the year
   real,    intent(in)                 :: utc_hour ! hour

   real, dimension(:,:), intent(in)    :: lon      ! longitudes, degrees
   real, dimension(:,:), intent(in)    :: lat      ! latitudes,  degrees

   real, dimension(:,:), intent(inout) :: sza      ! solar zenith angle, degrees
   real, dimension(:,:), intent(inout) :: cos_sza  ! cos(solar zenith angle)


! !OUTPUT PARAMETERS:


! !DESCRIPTION: Computes solar zenith angle
!
! !REVISION HISTORY:
!
!  28Sep2012 A. Darmenov   Addopted existing code from SulfateChemDriverMod.
!
!EOP
!-------------------------------------------------------------------------

   ! parameters
   real, parameter :: pi        = 3.1415926
   real, parameter :: f_rad2deg = 180.0/pi
   real, parameter :: f_deg2rad = pi/180.0

   real, parameter :: a0 = 0.006918
   real, parameter :: a1 = 0.399912
   real, parameter :: a2 = 0.006758
   real, parameter :: a3 = 0.002697
   real, parameter :: b1 = 0.070257
   real, parameter :: b2 = 0.000907
   real, parameter :: b3 = 0.000148

   ! local
   real    :: r
   real    :: solar_declination, sin_sd, cos_sd
   real    :: local_time
   real    :: hour_angle, cos_ha
   real    :: lat_
   integer :: i, i1, i2
   integer :: j, j1, j2


   r  = (2 * pi * (doy - 1))/365

   ! solar declination in radians
   solar_declination = a0 - a1*cos(  r) + b1*sin(  r) &
                          - a2*cos(2*r) + b2*sin(2*r) &
                          - a3*cos(3*r) + b3*sin(3*r)

   sin_sd = sin(solar_declination)
   cos_sd = cos(solar_declination)

   i1 = lbound(lon, dim=1); i2 = ubound(lon, dim=1)
   j1 = lbound(lon, dim=2); j2 = ubound(lon, dim=2)

   do j = j1, j2
       do i = i1, i2
           local_time = utc_hour + lon(i,j)/15

           if(local_time <  0) local_time = local_time + 24
           if(local_time > 24) local_time = local_time - 24

           hour_angle = (abs(local_time - 12) * 15) * f_deg2rad

           cos_ha = cos(hour_angle)

           lat_ = lat(i,j) * f_deg2rad
           cos_sza(i,j) = sin(lat_)*sin_sd + cos(lat_)*cos_sd*cos_ha
       end do
   end do

   sza = acos(cos_sza) * f_rad2deg
   where (cos_sza < 0) cos_sza = 0.0

 end subroutine solar_zenith_angle

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  day_of_year -- given nymd compute the day number of the yea
!
! !INTERFACE:
!

 pure function day_of_year(nymd) result (doy)

! !USES:

   implicit none

   integer :: doy

! !INPUT PARAMETERS:

   integer, intent(in) :: nymd

! !OUTPUT PARAMETERS:

! !DESCRIPTION: Computes day of year
!
! !REVISION HISTORY:
!
!  11Nov2012 A. Darmenov   Addopted existing code from SulfateChemDriverMod.
!
!EOP
!-------------------------------------------------------------------------
   ! parameters
   integer, parameter :: days(12) = (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)

   ! local
   integer :: yyyy
   integer :: mm
   integer :: dd
   integer :: month
   logical :: leap

   yyyy = nymd / 10000
   mm   = mod(nymd, 10000) / 100
   dd   = mod(nymd,   100)

   ! is it a leap year?
   leap = .false.

   if (mod(yyyy, 4) == 0) then
       leap = .true.
       if (mod(yyyy, 100) == 0) then
           leap = .false.
           if (mod(yyyy, 400) == 0) then
               leap = .true.
           endif
       endif
   endif

   ! calculate day of year
   doy = 0

   if (mm == 1) then
       doy = dd
   else
       do month = 1, mm - 1
           if ( (month == 2) .and. leap ) then
               doy = doy + 29
           else
               doy = doy + days(month)
           endif
       enddo

       doy = doy + dd
   endif

 end function day_of_year


 end module GEOS_AChemGridCompMod
