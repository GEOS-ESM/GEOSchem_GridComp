! $Id$

#include "MAPL_Generic.h"

!=============================================================================

module GEOS_PChemGridCompMod

!BOP

! !MODULE: GEOS_PChemGridCompMod

! !DESCRIPTION: 
!   {\tt GEOS\_PChem} is a proxy component for Aerochem that implements the 
!    specification or simple parameterization of the Aerochem species. It works
!    on three types of species: chemical species (oxygen, nitrous oxide, CFC-11,
!    CFC-12, CFC-22, methane, water vapor), diagnostic species (age-of-air),
!    and aerosols (arbitrary).
!
!    Each of the chemical species can be treated in one
!    of two ways: parameterized prediction from tabled zonally-symmetric production and 
!    loss (P-L) data, or specifycation from zonally-symmetric values (see Resources section
!    for how to control thios behavior).  A single flat file containing 
!    both the P-L and climatology data {\it must} be provided (see Resources section).
!    Aerosols are set to 3-dimensional climatological values.
!    The ``age-of-air'' is predicted by setting the surface values of this tracer
!    to the to zero and advancing other levels by dt.
!    All of these quantities except water vapor are INTERNAL state variables of
!    {\tt GEOS\_PChem}. Water vapor is assumed to be a Friendly Import and {\tt GEOS\_PChem}
!    leaves it unmodified below the tropopause, or the 200 hPa level if the tropopause
!    is below this level. 
!
!
!    For chemical species, the production rate is tabled directly.
!    For the loss, a rate coefficient is tabled.  Using Odd-oxygen $O_x$ as an example,
!    the species are updated as follows:
!    $$
!    \frac{\partial O_x}{\partial t} = 
!        \dot{Q}_{o} - \kappa_{o} O_x
!    $$
!    where $O_x$ is the specific mass of odd oxygen, $ \dot{Q}_{o}$ is the
!    odd oxygen production rate, $\kappa_{o}$ is the tabled rate 
!    coefficient for odd oxygen loss. This is finite differenced in time as:
!    $$
!      O_x^{n+1} =  \frac{O_x^{n} + \Delta t \dot{Q}_{o} }{1 + \Delta t \kappa_{o}}
!    $$
!
!    The component reads the monthly tables of the zonally averaged climatology of concentration
!    and production rates and loss frequencies in intialize. These are saved in the private 
!    internal state, which is static. The climatologies are interpolated to the natural 
!    locations and updated in the run method, and are kept in INTERNAL, an ESMF state attached to the
!    GEOS GENERIC object in the component. If no restart is specified for the
!    INTERNAL, the species are initialized to zero.
!
!    We have added a generalization that allows for a multiple-year climatology. Add
!     pchem\_clim\_years: nnn
!     pchem\_clim: dsn
!    to the AGCM.tmpl, where nnn is the number of years in the climatology and dsn is
!    the data set name.  We do not allow use of the production and loss parameterization when 
!    pchem\_clim\_years is greater than one.
!
!    Ozone is diagnosed from $O_x$ by assuming that it accounts for all
!    $O_x$ at pressures greater than 100 Pa (1 hPa) during the day and at all 
!    pressures at night. For those daylit cells where pressures are less than 1 hPa, 
!    we assume that the ozone fraction in $O_x$ decreases exponentially with decreasing
!    pressure.
!    
!    Aerosols are read from 3-dimensional data files that have to be on model levels
!    but may be on any regular lat-lon grid. These are hdf files and horizontal
!    interpolation is done through CFIO. The aerosols are read into a bundle that
!    is exported, and the number and names of the aerosols in the bundle are set
!    from the CFIO file. 
!\newline
!
!
!    \textbf{PC: ProTeX can parse the list of resources directly from the source code, if the
!    MAPL\_GetResouce calls are enclosed within the !BOR/!EOR. See, for example,} 
!    `\textit{RESOURCES:}' \textbf{below.}\\\\
!
!
!
!  !RESOURCES:
!     >> RUN\_DT:  none    real    seconds
!        Heartbeat. {\tt GEOS\_PChem} is called all the time.
!     >> pchem\_clim:                'pchem\_clim.dat'  string    none 
!           Zonally-symmetric chemistry data flat file.  
!     >> {\it{name}}\_FIXED\_VALUE:    use\_file          real      pppv 
!           Constant value at which to fix chemical species {\it{name}}.  If not specified, 
!         the file data will be used. If specified, {\it{name}}\_RELAXTIME: is ignored. 
!          {\it{name}} can be any of OX, CH4, N2O, CFC11, CFC12, HCFC22.  
!     >> {\it{name}}\_RELAXTIME:       0.0                real       seconds 
!         Timescale of relaxation to climatology on file for chemical species {\it{name}}. 
!         For values $<= 0$, the P-L parameterization will be used. To hold at the file's 
!         zonally-symmetric climatology, use a small positive number. 
!         {\it{name}} can be any of OX, CH4, N2O, CFC11, CFC12, HCFC22, H2O.  
!     >> {\it{name}}\_PCRIT:          1.e+16             real      Pa 
!           Pressure of level above which the relaxation to climatology is done. This is 
!         ignored if {\it{name}}\_RELAXTIME: is ignored or if {\it{name}}\_RELAXTIME: is $<= 0$. 
!         {\it{name}} can be any of OX, CH4, N2O, CFC11, CFC12, HCFC22.  
!     >> {\it{name}}\_DELP:           1.e-16             real      Pa 
!           Pressure interval over which the relaxation to climatology is ramped-in. This is 
!         ignored if {\it{name}}\_RELAXTIME: is ignored or if {\it{name}}\_RELAXTIME: is $<= 0$ 
!         {\it{name}} can be any of OX, CH4, N2O, CFC11, CFC12, HCFC22.  
!     >> {\it{name}}\_FRIENDLIES:     self              string    none 
!           String of colon separated component names to which this species is Friendly. 
!         {\it{name}} can be any of OX, CH4, N2O, CFC11, CFC12, HCFC22  
!     >> AOA\_FRIENDLIES:           'DYNAMICS:TURBULENCE'  string     none 
!           String of colon separated component names to which Age-of-Air is Friendly.  
!
! !USES:

  use ESMF
  use MAPL_Mod
  use Chem_Mod
  use ESMF_CFIOFileMOD
  use MAPL_CFIOMOD
  
  implicit none
  private
#include "mpif.h"

! !PUBLIC MEMBER FUNCTIONS:

  public SetServices

!EOP

!=============================================================================

  type T_Pchem_STATE
     private
     integer                             :: NLATS, NLEVS
     real, pointer, dimension(:)         :: LATS => null()
     real, pointer, dimension(:)         :: LEVS => null()
     real, pointer, dimension(:,:,:,:,:) :: MNPL => null() ! Production rates and loss frequencies
     real, pointer, dimension(:,:,:,:)   :: MNCV => null() ! Concentration (mole fraction)
     integer                             :: OX       = 1
     integer                             :: N2O      = 2
     integer                             :: CFC11    = 3
     integer                             :: CFC12    = 4
     integer                             :: CH4      = 5
     integer                             :: HCFC22   = 6
     integer                             :: H2O      = 7
     integer                             :: NSPECIES = 7
     character(len=ESMF_MAXSTR), dimension(7) :: ITEMNAMES= (/'OX    ','N2O   ','CFC11 ','CFC12 ','CH4   ','HCFC22','H2O   '/)

     INTEGER                             :: climYears
     INTEGER                             :: begClimYear
     INTEGER                             :: endClimYear

     INTEGER                             :: dayOfMonth = -1
     type(ESMF_Time)                     :: lastTimeHere
  end type T_Pchem_STATE

  type Pchem_WRAP
     type (T_Pchem_STATE), pointer :: PTR
  end type Pchem_WRAP

contains

!=============================================================================

!BOP

! !IROUTINE: SetServices

! !DESCRIPTION: Sets Initialize and Run services. 
! \newline
!

! !INTERFACE:

  subroutine SetServices ( GC, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
    integer, optional                  :: RC  ! return code

!EOP

!=============================================================================
!
! ErrLog Variables


    character(len=ESMF_MAXSTR)              :: IAm
    integer                                 :: STATUS
    character(len=ESMF_MAXSTR)              :: COMP_NAME

! Local derived type aliases

    type (ESMF_Config          )            :: CF

    type (T_PCHEM_STATE), pointer           :: PCHEM_state 
    type (T_PCHEM_STATE)                    :: DUMMY
    type (PCHEM_wrap)                       :: wrap
    type(Chem_Registry)                     :: chemReg   
    character(len=ESMF_MAXSTR)              :: OXFRIENDLY
    character(len=ESMF_MAXSTR)              :: N2OFRIENDLY
    character(len=ESMF_MAXSTR)              :: CFC11FRIENDLY
    character(len=ESMF_MAXSTR)              :: CFC12FRIENDLY
    character(len=ESMF_MAXSTR)              :: HCFC22FRIENDLY
    character(len=ESMF_MAXSTR)              :: CH4FRIENDLY
    character(len=ESMF_MAXSTR)              :: AOAFRIENDLY

    character(len=ESMF_MAXSTR)              :: FRIENDLIES
    character(len=ESMF_MAXSTR)              :: providerName

    LOGICAL :: Doing_RATs
    integer IO3AINC

!=============================================================================

! Begin...

! Get my name and set-up traceback handle
! ---------------------------------------

    call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME)//'::SetServices'

!   Start by loading the Chem Registry
!   ----------------------------------
    chemReg = Chem_RegistryCreate ( STATUS )
    VERIFY_(STATUS)

!   If not doing PChem, use GEOS Generic stubs from this point on
!   -------------------------------------------------------------
    if ( .NOT. chemReg%doing_PC ) then
       call MAPL_GenericSetServices ( GC, RC=STATUS )
       VERIFY_(STATUS)
       call Chem_RegistryDestroy ( chemReg, RC=STATUS )
       VERIFY_(STATUS)
       if (MAPL_AM_I_ROOT()) & 
           print *, trim(Iam)//': not ACTIVE, defaulting to GG stubs...'
       RETURN_(ESMF_SUCCESS)
    end if       

! Get the configuration
! ---------------------

    call ESMF_GridCompGet( GC, CONFIG = CF, RC=STATUS )
    VERIFY_(STATUS)

    call ESMF_ConfigGetAttribute( CF, IO3AINC, Label='ALLOW_PCHEM_AINC_UPDATE:', default=0,        RC=STATUS)
    VERIFY_(STATUS)

! Set the Initialize and Run entry point
! --------------------------------------

    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_INITIALIZE, Initialize, RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_GridCompSetEntryPoint ( GC,  ESMF_METHOD_RUN, Run,        RC=STATUS)
    VERIFY_(STATUS)

    if ( IO3AINC/=0 ) then
       call MAPL_GridCompSetEntryPoint ( GC,  ESMF_METHOD_RUN, AINC_UPDATE,    RC=STATUS)
       VERIFY_(STATUS)
    endif

! Allocate this instance of the internal state and put it in wrapper.
! -------------------------------------------------------------------

    allocate( PCHEM_state, stat=STATUS )
    VERIFY_(STATUS)

    WRAP%PTR => PCHEM_STATE
    PCHEM_STATE = DUMMY
 
! Save pointer to the wrapped internal state in the GC
! ----------------------------------------------------

    call ESMF_UserCompSetInternalState ( GC, 'Pchem_state', WRAP, STATUS )
    VERIFY_(STATUS)

! Set the state variable specs.
! -----------------------------

     FRIENDLIES = trim(COMP_NAME)

     call ESMF_ConfigGetAttribute(CF, OXFRIENDLY, Label='OX_FRIENDLIES:'      ,&
                                  default=FRIENDLIES                          ,&
                                                                     RC=STATUS )
     VERIFY_(STATUS)

     call ESMF_ConfigGetAttribute(CF, N2OFRIENDLY, Label='N2O_FRIENDLIES:'    ,&
                                  default=FRIENDLIES                          ,&
                                                                     RC=STATUS )
     VERIFY_(STATUS)

     call ESMF_ConfigGetAttribute(CF, CFC11FRIENDLY, Label='CFC11_FRIENDLIES:',&
                                  default=FRIENDLIES                          ,&
                                                                     RC=STATUS )
     VERIFY_(STATUS)

     call ESMF_ConfigGetAttribute(CF, CFC12FRIENDLY, Label='CFC12_FRIENDLIES:',&
                                  default=FRIENDLIES                          ,&
                                                                     RC=STATUS )
     VERIFY_(STATUS)

     call ESMF_ConfigGetAttribute(CF,HCFC22FRIENDLY,Label='HCFC22_FRIENDLIES:',&
                                  default=FRIENDLIES                          ,&
                                                                     RC=STATUS )
     VERIFY_(STATUS)

     call ESMF_ConfigGetAttribute(CF, CH4FRIENDLY, Label='CH4_FRIENDLIES:'    ,&
                                  default=FRIENDLIES                          ,&
                                                                     RC=STATUS )
     VERIFY_(STATUS)

     call ESMF_ConfigGetAttribute(CF, AOAFRIENDLY, Label='AOA_FRIENDLIES:'    ,&
                                  default=FRIENDLIES                          ,&
                                                                     RC=STATUS )
     VERIFY_(STATUS)


!BOS

! !IMPORT STATE:

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME = 'PLE',                                       &
        LONG_NAME  = 'air_pressure',                              &
        UNITS      = 'Pa',                                        &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationEdge,                          &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME ='Q',                                          &
        LONG_NAME  ='specific_humidity',                          &
        UNITS      ='kg kg-1',                                    &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,                        &  
                                                       RC=STATUS  )
     VERIFY_(STATUS)                                                                          

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'TROPP',                             &
        LONG_NAME          = 'tropopause_pressure',               &
        UNITS              = 'Pa',                                &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

    if ( IO3AINC /=0 ) then
       ! The following import is only for offline purposes
       ! NOTE: This is only used in offline applications so when adding new 
       !       fields to IMPORT state the suggestion is to add them BEFORE
       !       this state - unlike the usual procedure of always appending
       !       to the end of the state.
       call MAPL_AddImportSpec(GC,                                    &
            SHORT_NAME = 'O3AINC',                                    &
            LONG_NAME  = 'ozone_analysis_increment',                  &
            UNITS      = 'kg kg-1',                                   &
            default    = 0.0,                                         &
            DIMS       = MAPL_DimsHorzVert,                           &
            VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
       VERIFY_(STATUS)
    endif

! !INTERNAL STATE:

! For odd-oxygen only:  Ox is the second member of the analysis bundle, and 
! if the ANALYSIS_OX_PROVIDER is PCHEM, Ox must be friendly to "ANALYSIS".
! If Ox is not in TRANA, AGCM will fail after it counts the number of members 
! in TRANA [ASSERT(NumFriendly == 2)].
! --------------------------------------------------------------------------

     CALL ESMF_ConfigGetAttribute(CF, providerName, Default="PCHEM", &
                                  Label="ANALYSIS_OX_PROVIDER:", RC=STATUS )
     VERIFY_(STATUS)
     
     IF( providerName == "PCHEM" .AND. (INDEX(OXFRIENDLY,"ANALYSIS") == 0) ) THEN
      IF(MAPL_AM_I_ROOT()) THEN
       PRINT *," "
       PRINT *,TRIM(Iam)//": OX_FRIENDLIES in AGCM.tmpl invalid."
       PRINT *,"      You must at least specificy:  OX_FRIENDLIES:ANALYSIS"
       PRINT *,"      You may also add other component if desired. Example: DYNAMICS:ANALYSIS."
       PRINT *," "
      END IF
      STATUS = 1
      VERIFY_(STATUS)
     END IF

! Add species to the internal state only if PCHEM is the RATS provider
! --------------------------------------------------------------------
     CALL ESMF_ConfigGetAttribute(CF, providerName, Default="PCHEM", &
                                  Label="RATS_PROVIDER:", RC=STATUS )

     IF(providerName == "PCHEM") THEN
      Doing_RATs = .TRUE.
     ELSE
      Doing_RATs = .FALSE.
     END IF

     AddingRATS: IF(Doing_RATs) THEN

     call MAPL_AddInternalSpec(GC,                           &
        SHORT_NAME         = 'OX',                                &
        LONG_NAME          = 'odd_oxygen_volume_mixing_ratio',    &
        UNITS              = 'mol mol-1',                         &
        FRIENDLYTO         = OXFRIENDLY,                          &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddInternalSpec(GC,                           &
        SHORT_NAME         = 'N2O',                               &
        LONG_NAME          = 'nitrous_oxide_volume_mixing_ratio', &
        UNITS              = 'mol mol-1',                         &
        FRIENDLYTO         = N2OFRIENDLY,                         &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddInternalSpec(GC,                           &
        SHORT_NAME         = 'CFC11',                             &
        LONG_NAME          = 'CFC11_(CCl3F)_volume_mixing_ratio', &
        UNITS              = 'mol mol-1',                         &
        FRIENDLYTO         = CFC11FRIENDLY,                       &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddInternalSpec(GC,                           &
        SHORT_NAME         = 'CFC12',                             &
        LONG_NAME          = 'CFC12_(CCl2F2)_volume_mixing_ratio',&
        UNITS              = 'mol mol-1',                         &
        FRIENDLYTO         = CFC12FRIENDLY,                       &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddInternalSpec(GC,                           &
        SHORT_NAME         = 'HCFC22',                            &
        LONG_NAME          = 'HCFC22_(CHClF2)_volume_mixing_ratio', &
        UNITS              = 'mol mol-1',                         &
        FRIENDLYTO         = HCFC22FRIENDLY,                      &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddInternalSpec(GC,                           &
        SHORT_NAME         = 'CH4',                               &
        LONG_NAME          = 'methane_volume_mixing_ratio',       &
        UNITS              = 'mol mol-1',                         &
        FRIENDLYTO         = CH4FRIENDLY,                         &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     END IF AddingRATS

     call MAPL_AddInternalSpec(GC,                           &
        SHORT_NAME         = 'AOA',                               &
        LONG_NAME          = 'age_of_air',                        &
        UNITS              = 'days',                              &
        FRIENDLYTO         = AOAFRIENDLY,                         &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                       RC=STATUS  )
     VERIFY_(STATUS)


! !EXPORT STATE:

! Add species exports only if PCHEM is the RATS provider
! ------------------------------------------------------
     AddingRATsExports: IF(Doing_RATs) THEN

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'OX_TEND',                           &
        LONG_NAME          = 'tendency_of_odd_oxygen_mixing_ratio_due_to_chemistry', &
        UNITS              = 'mol mol-1 s-1',                     &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'H2O_TEND',                          &
        LONG_NAME          = 'tendency_of_water_vapor_mixing_ratio_due_to_chemistry', &
        UNITS              = 'kg kg-1 s-1',                       &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'OX_PROD',                           &
        LONG_NAME          = 'tendency_of_odd_oxygen_volume_mixing_ratio_due_to_production', &
        UNITS              = 'mol mol-1 s-1',                     &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'OX_LOSS',                           &
        LONG_NAME          = 'tendency_of_odd_oxygen_volume_mixing_ratio_due_to_loss', &
        UNITS              = 'mol mol-1 s-1',                     &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'N2O_PROD',                          &
        LONG_NAME          = 'tendency_of_nitrous_oxide_volume_mixing_ratio_due_to_production', &
        UNITS              = 'mol mol-1 s-1',                     &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'N2O_LOSS',                          &
        LONG_NAME          = 'tendency_of_nitrous_oxide_volume_mixing_ratio_due_to_loss', &
        UNITS              = 'mol mol-1 s-1',                     &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'CFC11_PROD',                        &
        LONG_NAME          = 'tendency_of_CFC11_volume_mixing_ratio_due_to_production', &
        UNITS              = 'mol mol-1 s-1',                     &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'CFC11_LOSS',                        &
        LONG_NAME          = 'tendency_of_CFC11_volume_mixing_ratio_due_to_loss', &
        UNITS              = 'mol mol-1 s-1',                     &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'CFC12_PROD',                        &
        LONG_NAME          = 'tendency_of_CFC12_volume_mixing_ratio_due_to_production', &
        UNITS              = 'mol mol-1 s-1',                     &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'CFC12_LOSS',                        &
        LONG_NAME          = 'tendency_of_CFC12_volume_mixing_ratio_due_to_loss', &
        UNITS              = 'mol mol-1 s-1',                     &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'HCFC22_PROD',                       &
        LONG_NAME          = 'tendency_of_HCFC22_volume_mixing_ratio_due_to_production', &
        UNITS              = 'mol mol-1 s-1',                     &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'HCFC22_LOSS',                       &
        LONG_NAME          = 'tendency_of_HCFC22_volume_mixing_ratio_due_to_loss', &
        UNITS              = 'mol mol-1 s-1',                     &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'CH4_PROD',                          &
        LONG_NAME          = 'tendency_of_methane_volume_mixing_ratio_due_to_production', &
        UNITS              = 'mol mol-1 s-1',                     &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'CH4_LOSS',                          &
        LONG_NAME          = 'tendency_of_methane_volume_mixing_ratio_due_to_loss', &
        UNITS              = 'mol mol-1 s-1',                     &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'H2O_PROD',                          &
        LONG_NAME          = 'tendency_of_specific_humidity_due_to_production', &
        UNITS              = 's-1',                               &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'H2O_LOSS',                          &
        LONG_NAME          = 'tendency_of_specific_humidity_due_to_loss', &
        UNITS              = 's-1',                               &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'O3',                                &
        LONG_NAME          = 'ozone_mass_mixing_ratio',           &
        UNITS              = 'kg kg-1',                           &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'O3PPMV',                            &
        LONG_NAME          = 'ozone_volume_mixing_ratio',         &
        UNITS              = 'ppmv',                              &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'TO3',                               &
        LONG_NAME          = 'total_column_ozone',                &
        UNITS              = 'Dobsons',                           &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'TTO3',                              &
        LONG_NAME          = 'tropospheric_column_ozone',         &
        UNITS              = 'Dobsons',                           &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     END IF AddingRATsExports



!EOS

! Set the Profiling timers
! ------------------------
    call MAPL_TimerAdd ( GC, name = "RUN",        RC=STATUS )
    VERIFY_(STATUS)
    call MAPL_TimerAdd ( GC, name = "-Read Species", RC=STATUS )
    VERIFY_(STATUS)
    call MAPL_TimerAdd ( GC, name = "INITIALIZE", RC=STATUS )
    VERIFY_(STATUS)
    call MAPL_TimerAdd ( GC, name = "-Read Header", RC=STATUS )
    VERIFY_(STATUS)

! Generic Set Services
! --------------------

    call MAPL_GenericSetServices ( GC,RC=STATUS )
    VERIFY_(STATUS)

    RETURN_(ESMF_SUCCESS)
  
  end subroutine SetServices


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!BOP

! !IROUTINE: INITIALIZE

! !DESCRIPTION: The Initialize method of Pchem gridded component.
!   It reads the production-loss file, which by default is 
!   {\tt pchem\_clim.dat}, but can be overridden from the configuration.
!   This version reads zonal-mean monthly climatologies and interpolates them to
!   the latitudes of the component's natural grid, which is the inherited grid.
!   \newline
!
  
! !INTERFACE:

  subroutine Initialize ( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
    type(ESMF_State),    intent(inout) :: IMPORT ! Import state
    type(ESMF_State),    intent(inout) :: EXPORT ! Export state
    type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
    integer, optional,   intent(  out) :: RC     ! Error code

!EOP

! ErrLog Variables

    character(len=ESMF_MAXSTR)              :: IAm 
    integer                                 :: STATUS
    character(len=ESMF_MAXSTR)              :: COMP_NAME
    
! Local derived type aliases

    type (MAPL_MetaComp), pointer           :: MAPL
    type (T_Pchem_STATE    ), pointer       :: Pchem_STATE 
    type (Pchem_wrap)                       :: WRAP
    type (ESMF_DELayout)                    :: layout
    type(ESMF_Time)                         :: PRVMONTH
    type(ESMF_State)                        :: INTERNAL
    type(ESMF_Alarm)                        :: PCHEM_ALARM
    type(ESMF_VM)                           :: VM

    character(len=ESMF_MAXSTR)              :: PCHEMFILE
    character(len=ESMF_MAXSTR)              :: providerName
    real, pointer                           :: LATS(:,:)
    integer                                 :: UNIT
    integer                                 :: NSPECIES
    integer                                 :: IM, JM, LM
    integer                                 :: dimid, varid, climYears, comm, info
    logical                                 :: Doing_RATs

!=============================================================================

! Begin... 

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

    call ESMF_GridCompGet ( GC, name=COMP_NAME, VM=VM, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME)//'::Initialize'

! Get my internal MAPL_Generic state
!-----------------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
    VERIFY_(STATUS)

! Start timer
!------------      

    call MAPL_TimerOn (MAPL,"INITIALIZE"  )

! Get my private state from the component
!----------------------------------------

    call ESMF_UserCompGetInternalState(gc, 'Pchem_state', WRAP, STATUS)
    VERIFY_(STATUS)

    Pchem_STATE => WRAP%PTR

! Call Generic Initialize
!------------------------

    call MAPL_GenericInitialize ( GC, IMPORT, EXPORT, CLOCK,  RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_TimerOn (MAPL,"TOTAL"  )

! Get latitudes from generic state.
!----------------------------------

    call MAPL_Get(MAPL,             &
         LAYOUT    = LAYOUT,                     &
         LATS      = LATS,                       &
         IM=IM, JM=JM, LM=LM,                    &
         INTERNAL_ESMF_STATE=INTERNAL,           &
                                       RC=STATUS )
    VERIFY_(STATUS)


! Is PCHEM the RATs provider?
! ---------------------------
    CALL MAPL_GetResource(MAPL, providerName, LABEL="RATS_PROVIDER:", DEFAULT="PCHEM", RC=STATUS )

    IF(providerName == "PCHEM") THEN
     Doing_RATs = .TRUE.
    ELSE
     Doing_RATs = .FALSE.
    END IF

    NeedRATsFile: IF(Doing_RATs) THEN

! Get file for monthly climatology of production rates and loss frequencies.
!---------------------------------------------------------------------------

    call MAPL_GetResource(MAPL, PCHEMFILE,'pchem_clim:' ,DEFAULT='pchem_clim.dat', RC=STATUS )
    VERIFY_(STATUS)

! Let us know if the above file contains more than one year of data
!------------------------------------------------------------------

    call MAPL_GetResource(MAPL, PCHEM_STATE%climYears, 'pchem_clim_years:' ,DEFAULT=1, RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_TimerOn (MAPL,"-Read Header"  )

#ifdef H5_HAVE_PARALLEL
    call MPI_Info_create(info, STATUS)
    call MPI_Info_set(info, "romio_cb_read", "automatic", STATUS)
    VERIFY_(STATUS)
    call ESMF_VMGet(vm, mpiCommunicator=comm, rc=STATUS)
    VERIFY_(STATUS)

#ifdef NETCDF_NEED_NF_MPIIO
    STATUS = NF_OPEN_PAR(trim(PCHEMFILE),IOR(NF_NOWRITE,NF_MPIIO),comm,info,UNIT)
#else
    STATUS = NF_OPEN_PAR(trim(PCHEMFILE),NF_NOWRITE,comm,info,UNIT)
#endif

#else
    if ( MAPL_am_I_root() ) then
       STATUS = NF_OPEN(trim(PCHEMFILE),NF_NOWRITE,UNIT)
#endif
    if(status /= nf_noerr) then
       print*,'Error opening file ',trim(PCHEMFILE), status
       print*, NF_STRERROR(status)
       stop
    endif

! Read various dimension and global information including
! number of levels, and the number of species in the PCHEMFILE.
!--------------------------------------------------------------

    STATUS = NF_INQ_DIMID(UNIT, 'lat', dimid)
    if(status /= nf_noerr) then
       print*,'Error getting dimid for lat', status
       print*, NF_STRERROR(status)
       stop
    endif
    STATUS = NF_INQ_DIMLEN(UNIT, dimid, PCHEM_STATE%NLATS)
    if(status /= nf_noerr) then
       print*,'Error getting dimlen for lat', status
       print*, NF_STRERROR(status)
       stop
    endif
    STATUS = NF_INQ_DIMID(UNIT, 'lev', dimid)
    if(status /= nf_noerr) then
       print*,'Error getting dimid for lev', status
       print*, NF_STRERROR(status)
       stop
    endif
    STATUS = NF_INQ_DIMLEN(UNIT, dimid, PCHEM_STATE%NLEVS)
    if(status /= nf_noerr) then
       print*,'Error getting dimlen for lev', status
       print*, NF_STRERROR(status)
       stop
    endif
    STATUS = NF_GET_ATT_INT(UNIT, NF_GLOBAL, 'NSPECIES', NSPECIES)
    if(status /= nf_noerr) then
       print*,'Error getting NSPECIES', status
       print*, NF_STRERROR(status)
       stop
    endif
    ASSERT_(PCHEM_STATE%NSPECIES==NSPECIES)

! If multiple climYears, the second record tells us starting and ending year and number of years
!-----------------------------------------------------------------------------------------------

    IF(PCHEM_STATE%climYears == 1) THEN
       PCHEM_STATE%begClimYear = 1
       PCHEM_STATE%endClimYear = 1
    ELSE
       STATUS = NF_GET_ATT_INT(UNIT, NF_GLOBAL, 'begClimYear', PCHEM_STATE%begClimYear)
       if(status /= nf_noerr) then
          print*,'Error getting begClimYear', status
          print*, NF_STRERROR(status)
          stop
       endif
       STATUS = NF_GET_ATT_INT(UNIT, NF_GLOBAL, 'endClimYear', PCHEM_STATE%endClimYear)
       if(status /= nf_noerr) then
          print*,'Error getting endClimYear', status
          print*, NF_STRERROR(status)
          stop
       endif
       STATUS = NF_GET_ATT_INT(UNIT, NF_GLOBAL, 'climYears'  , climYears)
       if(status /= nf_noerr) then
          print*,'Error getting climYears', status
          print*, NF_STRERROR(status)
          stop
       endif

       IF(climYears /= PCHEM_STATE%climYears) THEN
          PRINT *," "
          PRINT *,TRIM(Iam)//": Problem with "//TRIM(PCHEMFILE)
          PRINT *,"Expecting ",PCHEM_STATE%climYears," years but there are ",climYears
          STATUS=1
          VERIFY_(STATUS)
       END IF
    END IF

! Allocate and read PCHEMFILE's latitudes (radians) and pressures (Pa).
!----------------------------------------------------------------------
#ifndef H5_HAVE_PARALLEL
    endif ! MAPL_am_I_root

    call MAPL_CommsBcast (vm, PCHEM_STATE%NLATS   ,1, 0, rc=status)
    VERIFY_(STATUS)
    call MAPL_CommsBcast (vm, PCHEM_STATE%NLEVS   ,1, 0, rc=status)
    VERIFY_(STATUS)
    call MAPL_CommsBcast (vm, PCHEM_STATE%begClimYear   ,1, 0, rc=status)
    VERIFY_(STATUS)
    call MAPL_CommsBcast (vm, PCHEM_STATE%endClimYear   ,1, 0, rc=status)
    VERIFY_(STATUS)

#endif

    allocate ( PCHEM_STATE%LATS (PCHEM_STATE%NLATS), STAT=STATUS )
    VERIFY_(STATUS)
    allocate ( PCHEM_STATE%LEVS (PCHEM_STATE%NLEVS), STAT=STATUS )
    VERIFY_(STATUS)

#ifndef H5_HAVE_PARALLEL
    if ( MAPL_am_I_root() ) then
#endif

    STATUS = NF_INQ_VARID(UNIT, 'lat', varid)
    if(status /= nf_noerr) then
       print*,'Error getting varid for lat', status
       print*, NF_STRERROR(status)
       stop
    endif
    STATUS = NF_GET_VAR_REAL(UNIT, varid, PCHEM_STATE%LATS)
    if(status /= nf_noerr) then
       print*,'Error getting values for lat', status
       print*, NF_STRERROR(status)
       stop
    endif
    STATUS = NF_INQ_VARID(UNIT, 'lev', varid)
    if(status /= nf_noerr) then
       print*,'Error getting varid for lev', status
       print*, NF_STRERROR(status)
       stop
    endif
    STATUS = NF_GET_VAR_REAL(UNIT, varid, PCHEM_STATE%LEVS)
    if(status /= nf_noerr) then
       print*,'Error getting values for lev', status
       print*, NF_STRERROR(status)
       stop
    endif
#ifdef H5_HAVE_PARALLEL
    call MPI_Info_free(info, status)
    VERIFY_(STATUS)
#else
    endif ! MAPL_am_I_root

    call MAPL_CommsBcast (vm, PCHEM_STATE%LATS,size(PCHEM_STATE%LATS), 0, rc=status)
    VERIFY_(STATUS)
    call MAPL_CommsBcast (vm, PCHEM_STATE%LEVS,size(PCHEM_STATE%LEVS), 0, rc=status)
    VERIFY_(STATUS)

#endif
    STATUS = NF_CLOSE(UNIT)
    call MAPL_TimerOff (MAPL,"-Read Header"  )

! Allocate concentration and production rates and loss frequencies.
! For now, in the case of multiple climYears, we will not allow P and L.
!-----------------------------------------------------------------------

    ALLOCATE(PCHEM_STATE%MNCV(PCHEM_STATE%NLATS, PCHEM_STATE%NLEVS, PCHEM_STATE%NSPECIES, 2), stat=STATUS )
    VERIFY_(STATUS)
    PCHEM_STATE%MNCV = Z'7FA00000'

    IF(PCHEM_STATE%climYears == 1) THEN
       ALLOCATE(PCHEM_STATE%MNPL(PCHEM_STATE%NLATS, PCHEM_STATE%NLEVS, PCHEM_STATE%NSPECIES, 2, 2), stat=STATUS )
       VERIFY_(STATUS)
       PCHEM_STATE%MNPL = Z'7FA00000'
    ENDIF

! Setting the alarm to ringing will reinitialize all data during first run
!-------------------------------------------------------------------------

    call ESMF_TimeSet(PRVMONTH,YY=1869,MM=1,DD=1, H=0, M=0, S=0, RC=STATUS)
    VERIFY_(STATUS)
    PCHEM_ALARM = ESMF_AlarmCreate(name='REFRESH_PCHEM_SPECIES', clock=CLOCK,      &
                                ringTime=PRVMONTH, sticky=.false.,     RC=STATUS)
    VERIFY_(STATUS)
    call ESMF_AlarmRingerOn(PCHEM_ALARM, rc=status)
    VERIFY_(STATUS)

    END IF NeedRATsFile



! Time
!-----
    CALL ESMF_ClockGet(CLOCK, currTime=PCHEM_STATE%lastTimeHere, RC=STATUS)
    VERIFY_(STATUS)


#ifdef PRINT_STATES

!   Print what my states are
!   ------------------------
    if ( MAPL_am_I_root() ) then

       print *,  trim(Iam)//": IMPORT State" 
                                             call ESMF_StatePrint ( IMPORT)
       print *,  trim(Iam)//": INTERNAL State" 
                                             call ESMF_StatePrint ( INTERNAL )
       print *,  trim(Iam)//": EXPORT State" 
                                             call ESMF_StatePrint ( EXPORT )

    end if

#endif


    call MAPL_TimerOff (MAPL,"INITIALIZE" )
    call MAPL_TimerOff (MAPL,"TOTAL"      )

! All Done
!---------

    RETURN_(ESMF_SUCCESS)
  end subroutine Initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!BOP

! !IROUTINE: RUN

! !DESCRIPTION: Updates mixing ratios of Ox, CH4, N2O, CFC11, CFC12, HCFC22.
!    For each species, it either updates a state variable based on tabled
!    production rates and loss times or based on tabled climatological values.
!    The latter is performed when a non-zero relaxation time for the species
!    is found in the configuration. The relaxation times (in seconds) can be specified 
!    with the labels XX\_RELAXTIME:, for example CFC12\_RELAXTIME:. The default
!    is zero, which results in doing a production-loss calculation. To fix
!    the species to climatology, simply specify a very short relaxation time.
!    Since the udpate is done implicitly, this can be done safely.
! \newline
!

! !INTERFACE:

subroutine RUN ( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:

  type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
  type(ESMF_State),    intent(inout) :: IMPORT ! Import state
  type(ESMF_State),    intent(inout) :: EXPORT ! Export state
  type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
  integer, optional,   intent(  out) :: RC     ! Error code:

!EOP


! ErrLog Variables

  character(len=ESMF_MAXSTR)        :: IAm
  integer                           :: STATUS
  character(len=ESMF_MAXSTR)        :: COMP_NAME

! Local derived types

  type (MAPL_MetaComp),     pointer :: MAPL
  type (ESMF_State)                 :: INTERNAL
  type (T_Pchem_STATE),     pointer :: Pchem_STATE 
  type (Pchem_wrap)                 :: WRAP
  type (MAPL_SunOrbit)              :: ORBIT
  type (ESMF_Time)                  :: CurrTime, dummyTIME
  type(ESMF_Alarm)                  :: PCHEM_ALARM
  type(ESMF_Alarm)                  :: RUN_ALARM
  type(ESMF_TimeInterval)           :: RingInterval
  type(ESMF_VM)                     :: VM

! Local scalars

  integer                           :: IM, JM, LM, L
  integer                           :: INDX1, INDX2, NLEVS, NLATS
  integer                           :: N
  integer                           :: YY, MM, DD
  integer                           :: CCYY
  integer                           :: start(3), cnt(3), UNIT, K, varid, comm, info
  real                              :: FAC
  real                              :: DT

  character(len=ESMF_MAXSTR)        :: PCHEMFILE
  character(len=ESMF_MAXSTR)        :: FieldName
  
! pointers to import

  real, pointer, dimension(:,:,:)   :: PLE
  real, pointer, dimension(:,:,:)   :: H2O
  real, pointer, dimension(  :,:)   :: TROPP

! pointers to internal

  real, pointer, dimension(:,:,:)   :: OX
  real, pointer, dimension(:,:,:)   :: CFC11
  real, pointer, dimension(:,:,:)   :: CFC12
  real, pointer, dimension(:,:,:)   :: HCFC22
  real, pointer, dimension(:,:,:)   :: CH4
  real, pointer, dimension(:,:,:)   :: N2O
  real, pointer, dimension(:,:,:)   :: AOA

! pointers to export

  real, pointer, dimension(:,:,:)   :: O3
  real, pointer, dimension(:,:,:)   :: O3PPMV
  real, pointer, dimension(:,:)     :: TO3
  real, pointer, dimension(:,:)     :: TTO3

! scratch arrays
  
  real, allocatable :: PL      (:,:,:)
  real, allocatable :: PROD_INT(:,:,:)
  real, allocatable :: LOSS_INT(:,:,:)
  real, allocatable :: PROD    (:,:  )
  real, allocatable :: LOSS    (:,:  )
  real, allocatable :: PROD1   (:,:  )
  real, allocatable :: LOSS1   (:,:  )
  real, allocatable :: ZTH     (:,:  )
  real, allocatable :: O3VMR   (:,:  )
  real, allocatable :: WRK     (:,:  )
  real, allocatable :: WGT     (:,:  )

! scratch pointers

  real, pointer     :: LATS    (:,:  )
  real, pointer     :: LONS    (:,:  )

! There are 2.69E+20 molecules per Dobson
! ----------------------------------------

  real, parameter   :: DOBSONS_PER_MOLE=MAPL_AVOGAD/2.69E+20

  type(ESMF_Time        ) :: midMonth
  type(ESMF_TimeInterval) :: oneMonth

  character(len=ESMF_MAXSTR) :: providerName
  logical                    :: Doing_RATs
  real(ESMF_KIND_R8)         :: dt_r8

!=============================================================================

! Begin... 

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

    Iam = "Run"
    call ESMF_GridCompGet( GC, name=COMP_NAME, VM=VM, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME)//"::Run"

! Retrieve the pointer to the generic state
!------------------------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
    VERIFY_(STATUS)

! Start timer
!------------

    call MAPL_TimerOn (MAPL,"TOTAL")
    call MAPL_TimerOn (MAPL,"RUN"  )
 
! Get RUN_ALARM from MAPL
!------------------------

    call MAPL_Get( MAPL, RUNALARM = RUN_ALARM, RC=STATUS )
    VERIFY_(STATUS)

! Get the time step from the RUN_ALARM
! ------------------------------------

    call ESMF_AlarmGet ( RUN_ALARM, ringInterval=RingInterval,RC=STATUS )
    VERIFY_(STATUS)
    call ESMF_TimeIntervalGet( RingInterval, s_r8=dt_r8, RC=STATUS )
    VERIFY_(STATUS)

    DT = real(dt_r8)

! Retrieve the pointer to the private internal state
!---------------------------------------------------

    call ESMF_UserCompGetInternalState(GC, 'Pchem_state', WRAP, STATUS)
    VERIFY_(STATUS)

    Pchem_STATE => WRAP%PTR

! Get parameters from generic state.
!-----------------------------------

    call MAPL_Get(MAPL,                          &
         IM=IM, JM=JM, LM=LM,                    &
         LONS      = LONS,                       &
         LATS      = LATS,                       &
         ORBIT     = ORBIT,                      &
         INTERNAL_ESMF_STATE=INTERNAL,           &
                                       RC=STATUS )
    VERIFY_(STATUS)

! Time
!-----
    CALL ESMF_ClockGet(CLOCK, currTime=CurrTime, RC=STATUS)
    VERIFY_(STATUS)

! Is PCHEM the RATs provider?
! ---------------------------
    CALL MAPL_GetResource(MAPL, providerName, LABEL="RATS_PROVIDER:", DEFAULT="PCHEM", RC=STATUS )

    IF(providerName == "PCHEM") THEN
     Doing_RATs = .TRUE.
    ELSE
     Doing_RATs = .FALSE.
    END IF

    ProvidingRATs: IF(Doing_RATs) THEN

! Allocate space for interpolation
!---------------------------------

    NLEVS = PCHEM_STATE%NLEVS
    NLATS = PCHEM_STATE%NLATS

    allocate(PROD_INT(IM,JM,LM),stat=STATUS)
    VERIFY_(STATUS)
    allocate(LOSS_INT(IM,JM,LM),stat=STATUS)
    VERIFY_(STATUS)
    allocate(      PL(IM,JM,LM),stat=STATUS)
    VERIFY_(STATUS)
    allocate(    PROD(IM,NLEVS),stat=STATUS)
    VERIFY_(STATUS)
    allocate(    LOSS(IM,NLEVS),stat=STATUS)
    VERIFY_(STATUS)
    allocate(   PROD1(NLATS,NLEVS),stat=STATUS)
    VERIFY_(STATUS)
    allocate(   LOSS1(NLATS,NLEVS),stat=STATUS)
    VERIFY_(STATUS)

! Time interpolation parameters
! -----------------------------
    
    call MAPL_ClimInterpFac(CLOCK, INDX1, INDX2, FAC, RC=STATUS)
    VERIFY_(STATUS)

! Read bracketing months, make sure INDX1 and INDX2 are in range. Annual
! cycle is preserved for years that precede and succeed the climatology.
! ----------------------------------------------------------------------

       N = 12*PCHEM_STATE%climYears

       CALL ESMF_TimeGet(CurrTime, YY=YY, MM=MM, DD=DD, RC=STATUS)
       VERIFY_(STATUS)

       CCYY = YY
       IF(CCYY < PCHEM_STATE%begClimYear) CCYY = PCHEM_STATE%begClimYear
       IF(CCYY > PCHEM_STATE%endClimYear) CCYY = PCHEM_STATE%endClimYear

       INDX1 = INDX1+(CCYY-PCHEM_STATE%begClimYear)*12
       INDX2 = INDX2+(CCYY-PCHEM_STATE%begClimYear)*12

       IF(MM ==  1 .AND. INDX1 > INDX2) INDX1 = INDX1-12
       IF(MM == 12 .AND. INDX2 < INDX1) INDX2 = INDX2+12

       IF(INDX1 == 0) INDX1 = 12
       IF(YY < PCHEM_STATE%begClimYear .AND. INDX2 > 12) INDX2 = 1

       IF(INDX2 == N+1) INDX2 = N-11
       IF(YY > PCHEM_STATE%endClimYear .AND. INDX1 == N-12) INDX1 = N

       call ESMF_ClockGetAlarm(CLOCK,'REFRESH_PCHEM_SPECIES', PCHEM_ALARM,RC=STATUS)
       VERIFY_(STATUS)

       if (currTime < PCHEM_STATE%lastTimeHere) then
          ! this should have not happen, unless we are doing replay and rewind clock
          call ESMF_AlarmRingerOn(PCHEM_ALARM, RC=STATUS)
          VERIFY_(STATUS)
       end if

       if ( ESMF_AlarmIsRinging( PCHEM_ALARM ) ) then

          call ESMF_AlarmRingerOff(PCHEM_ALARM, RC=STATUS)
          VERIFY_(STATUS)

          call MAPL_TimerOff(MAPL,"RUN"  )
          call MAPL_TimerOn (MAPL,"-Read Species"  )
          call MAPL_GetResource(MAPL, PCHEMFILE,'pchem_clim:' ,DEFAULT='pchem_clim.dat', RC=STATUS )
          VERIFY_(STATUS)

          call ESMF_VMGet(vm, mpiCommunicator=comm, rc=STATUS)
          VERIFY_(STATUS)

#ifdef H5_HAVE_PARALLEL
          call MPI_Info_create(info, STATUS)
          VERIFY_(STATUS)
          call MPI_Info_set(info, "romio_cb_read", "automatic", STATUS)
          VERIFY_(STATUS)

#ifdef NETCDF_NEED_NF_MPIIO
          STATUS = NF_OPEN_PAR(trim(PCHEMFILE),IOR(NF_NOWRITE,NF_MPIIO),comm,info,UNIT)
#else
          STATUS = NF_OPEN_PAR(trim(PCHEMFILE),NF_NOWRITE,comm,info,UNIT)
#endif

#else
          if ( MAPL_am_I_root() ) then
             STATUS = NF_OPEN(trim(PCHEMFILE),NF_NOWRITE,UNIT)
#endif
          if(status /= nf_noerr) then
             print*,'Error opening file ',trim(PCHEMFILE), status
             print*, NF_STRERROR(status)
             stop
          endif

          start(1) = 1
          start(2) = 1
          cnt(1) = PCHEM_STATE%NLATS
          cnt(2) = PCHEM_STATE%NLEVS
          cnt(3) = 1

          DO K = 1,PCHEM_STATE%NSPECIES
             FieldName = PCHEM_STATE%ITEMNAMES(K)
             STATUS = NF_INQ_VARID(UNIT, trim(FieldName), varid)
             if(status /= nf_noerr) then
                print*,'Error getting varid for variable ',trim(FieldName), status
                print*, NF_STRERROR(status)
                stop
             endif
! Need two separate reads because INDX2 isn't always sequentially after INDX1, otherwise
! we could combine the reads into one
             start(3) = INDX1
             STATUS = NF_GET_VARA_REAL(UNIT, varid, start, cnt, PCHEM_STATE%MNCV(:,:,K,1))
             if(status /= nf_noerr) then
                print*,'Error reading lower bracket month ',status
                print*, NF_STRERROR(status)
                stop
             endif
             start(3) = INDX2
             STATUS = NF_GET_VARA_REAL(UNIT, varid, start, cnt, PCHEM_STATE%MNCV(:,:,K,2))
             if(status /= nf_noerr) then
                print*,'Error reading upper bracket month ',status
                print*, NF_STRERROR(status)
                stop
             endif

! Convert H2O to mass fraction.
!------------------------------
             IF(K == PCHEM_STATE%H2O) then
                PCHEM_STATE%MNCV(:,:,K,1) = PCHEM_STATE%MNCV(:,:,K,1)*(MAPL_H2OMW/MAPL_AIRMW)
                PCHEM_STATE%MNCV(:,:,K,2) = PCHEM_STATE%MNCV(:,:,K,2)*(MAPL_H2OMW/MAPL_AIRMW)
             endif

! Production rates and loss frequencies. If multiple climYears, simply set to zero.
! ---------------------------------------------------------------------------------
             IF(PCHEM_STATE%climYears == 1) THEN

                STATUS = NF_INQ_VARID(UNIT, trim(FieldName)//'_PROD', varid)
                if(status /= nf_noerr) then
                   print*,'Error getting varid for variable ',trim(FieldName)//'_PROD', status
                   print*, NF_STRERROR(status)
                   stop
                endif
                start(3) = INDX1
                STATUS = NF_GET_VARA_REAL(UNIT, varid, start, cnt, PCHEM_STATE%MNPL(:,:,K,1,1))
                if(status /= nf_noerr) then
                   print*,'Error reading lower bracket month for production ',status
                   print*, NF_STRERROR(status)
                   stop
                endif
                start(3) = INDX2
                STATUS = NF_GET_VARA_REAL(UNIT, varid, start, cnt, PCHEM_STATE%MNPL(:,:,K,1,2))
                if(status /= nf_noerr) then
                   print*,'Error reading upper bracket month for production ',status
                   print*, NF_STRERROR(status)
                   stop
                endif
                IF(K == PCHEM_STATE%H2O) PCHEM_STATE%MNPL(:,:,K,1,1) = PCHEM_STATE%MNPL(:,:,K,1,1)*(MAPL_H2OMW/MAPL_AIRMW)
                IF(K == PCHEM_STATE%H2O) PCHEM_STATE%MNPL(:,:,K,1,2) = PCHEM_STATE%MNPL(:,:,K,1,2)*(MAPL_H2OMW/MAPL_AIRMW)
! Loss
! ----
                STATUS = NF_INQ_VARID(UNIT, trim(FieldName)//'_LOSS', varid)
                if(status /= nf_noerr) then
                   print*,'Error getting varid for variable ',trim(FieldName)//'_LOSS', status
                   print*, NF_STRERROR(status)
                   stop
                endif
                start(3) = INDX1
                STATUS = NF_GET_VARA_REAL(UNIT, varid, start, cnt, PCHEM_STATE%MNPL(:,:,K,2,1))
                if(status /= nf_noerr) then
                   print*,'Error reading lower bracket month for loss ',status
                   print*, NF_STRERROR(status)
                   stop
                endif
                start(3) = INDX2
                STATUS = NF_GET_VARA_REAL(UNIT, varid, start, cnt, PCHEM_STATE%MNPL(:,:,K,2,2))
                if(status /= nf_noerr) then
                   print*,'Error reading upper bracket month for loss ',status
                   print*, NF_STRERROR(status)
                   stop
                endif
             ENDIF

          ENDDO

          STATUS = NF_CLOSE(UNIT)
          VERIFY_(STATUS)

#ifdef H5_HAVE_PARALLEL
          call MPI_Info_free(info, status)
          VERIFY_(STATUS)
#else
          endif ! MAPL_am_I_root
          call MPI_Bcast (PCHEM_STATE%MNCV, size(PCHEM_STATE%MNCV), MPI_REAL, 0, comm, STATUS)
          VERIFY_(STATUS)
          IF(PCHEM_STATE%climYears == 1) THEN
             call MPI_Bcast (PCHEM_STATE%MNPL, size(PCHEM_STATE%MNPL), MPI_REAL, 0, comm, STATUS)
             VERIFY_(STATUS)
          ENDIF
#endif

          call MAPL_TimerOff (MAPL,"-Read Species"  )
          call MAPL_TimerOn  (MAPL,"RUN"  )

          call ESMF_TimeIntervalSet(oneMonth, MM = 1, RC=STATUS )
          VERIFY_(STATUS)
          call ESMF_TimeGet(currTime, midMonth=midMonth, RC=STATUS )
          VERIFY_(STATUS)

          if( currTime < midMonth ) then
             dummyTIME = CurrTime
          else
             dummyTIME = CurrTime + OneMonth
          endif
          call ESMF_TimeGet (dummyTIME, midMonth=midMonth,    RC=STATUS)
          VERIFY_(STATUS)
          call ESMF_AlarmSet(PCHEM_ALARM, ringtime=midMonth, RC=STATUS)
          VERIFY_(STATUS)

#ifdef DEBUG
          if(MAPL_AM_I_ROOT()) then
             print*,'Next ring time for SPECIES Alarm is'
             call ESMF_TimePrint(midMonth, "string", rc)
          endif
#endif

       endif

! Verify INDX1 and INDX2 once a day
! ---------------------------------
       IF(PCHEM_STATE%dayOfMonth /= DD) THEN
         IF(MAPL_AM_I_ROOT()) THEN
           PRINT *," "
           PRINT *,TRIM(Iam)//": Indices selected from climatology are ",INDX1,INDX2
           PRINT *," "
         END IF
         PCHEM_STATE%dayOFMonth = DD
       END IF

! Obtain GEOS5 mid-layer pressures from layer edge pressures (Pa).
!-----------------------------------------------------------------
    
    call MAPL_GetPointer( IMPORT, TROPP,  'TROPP', RC=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetPointer( IMPORT,   PLE,    'PLE', RC=STATUS )
    VERIFY_(STATUS)

    PL = 0.5*(PLE(:,:,0:LM-1)+PLE(:,:,1:LM))

! Update
!-------

    call UPDATE(PCHEM_STATE%CH4,    'CH4',    CH4   ) ! CH4=   1.75e-6
    call UPDATE(PCHEM_STATE%N2O,    'N2O',    N2O   ) ! N2O=   0.28e-6
    call UPDATE(PCHEM_STATE%CFC11,  'CFC11',  CFC11 ) ! CFC11= 0.30e-9 
    call UPDATE(PCHEM_STATE%CFC12,  'CFC12',  CFC12 ) ! CFC12= 0.50e-9 
    call UPDATE(PCHEM_STATE%HCFC22, 'HCFC22', HCFC22) ! HCFC22= 0.20e-9 
    call UPDATE(PCHEM_STATE%OX,     'OX',     OX    )

! Water
!------

    if(MAPL_VerifyFriendly(IMPORT,'Q','CHEMISTRY')) then
       call UPDATE(PCHEM_STATE%H2O,  'H2O'  ,H2O  )
    endif

! Ozone
!------

    call MAPL_GetPointer ( EXPORT,     O3,     'O3', RC=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetPointer ( EXPORT, O3PPMV, 'O3PPMV', RC=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetPointer ( EXPORT,    TO3,    'TO3', RC=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetPointer ( EXPORT,   TTO3,   'TTO3', RC=STATUS )
    VERIFY_(STATUS)

    if(associated(O3) .or. associated(O3PPMV)) then

       allocate(  ZTH(IM,JM),stat=STATUS)
       VERIFY_(STATUS)
       allocate(O3VMR(IM,JM),stat=STATUS)
       VERIFY_(STATUS)
       allocate(  WRK(IM,JM),stat=STATUS)
       VERIFY_(STATUS)
       allocate(  WGT(IM,JM),stat=STATUS)
       VERIFY_(STATUS)

       call MAPL_SunGetInsolation(LONS, LATS,  &
            ORBIT, ZTH, O3VMR, CLOCK=CLOCK,    & ! I dont need SLR, but it is not optional.
            RC=STATUS  )
       VERIFY_(STATUS)

       if(associated( TO3)) TO3  = 0.00
       if(associated(TTO3)) TTO3 = 0.00

       do L=1,LM

! Eliminate mesospheric ozone when sun is out
!--------------------------------------------

          where(PL(:,:,L) < 100.0 .and. ZTH > 0.0)
             O3VMR = OX(:,:,L)*exp(-1.5*(log10(PL(:,:,L))-2.0)**2)
          elsewhere
             O3VMR = OX(:,:,L)
          end where

          if(associated(O3    )) O3    (:,:,L) = O3VMR * (MAPL_O3MW / MAPL_AIRMW)
          if(associated(O3PPMV)) O3PPMV(:,:,L) = O3VMR * 1.e6

! Total ozone.  Sum on layer the following:  
!   O3(vmr) * Avogadro's number * dp / ( mwt air * g ),
!   which yields the number of molecules m^{-2}.
! ---------------------------------------------------
          
          if(associated(TTO3).or.associated( TO3)) then
             WRK = O3VMR*(PLE(:,:,L)-PLE(:,:,L-1))*(DOBSONS_PER_MOLE/(MAPL_AIRMW*MAPL_GRAV))

             if(associated( TO3))  TO3 =  TO3+WRK

             if(associated(TTO3)) then
                WGT  = max(0.0,min(1.0,(PLE(:,:,L)-TROPP)/(PLE(:,:,L)-PLE(:,:,L-1))))
                TTO3 = TTO3+WRK*WGT
             end if
          end if

       end do

       if(associated(TTO3)) then
          where(TROPP == MAPL_Undef) TTO3 = MAPL_Undef
       endif

       deallocate(ZTH)
       deallocate(O3VMR)
       deallocate(WRK)
       deallocate(WGT)

    end if

    END IF ProvidingRATs

! Age of air
!-----------

    call MAPL_GetPointer ( INTERNAL, AOA, 'AOA', RC=STATUS )
    VERIFY_(STATUS)

    AOA         = AOA +  (DT/86400.0) 
    AOA(:,:,LM) = 0.0

    PCHEM_STATE%lastTimeHere = currTime


! Clean-up
!---------

    IF(Doing_RATs) THEN
     deallocate(PROD_INT)
     deallocate(LOSS_INT)
     deallocate(      PL)
     deallocate(    PROD)
     deallocate(    LOSS)
     deallocate(   PROD1)
     deallocate(   LOSS1)
    END IF

! Stop timer
!-----------

    call MAPL_TimerOff(MAPL,"RUN"  )
    call MAPL_TimerOff(MAPL,"TOTAL")

!  All done
!-----------

   RETURN_(ESMF_SUCCESS)

contains

  subroutine UPDATE(NN,NAME,XX)

    integer,          intent(IN) :: NN
    character(len=*), intent(IN) :: NAME
    real, pointer                :: XX(:,:,:)


    real, pointer, dimension(:,:,:)   :: XX_PROD
    real, pointer, dimension(:,:,:)   :: XX_LOSS
    real, pointer, dimension(:,:,:)   :: OX_TEND
    real, pointer, dimension(:,:,:)   :: H2O_TEND
    real                              :: TAU
    real                              :: VALUE
    real                              :: PCRIT
    real                              :: DELP
    integer                           :: I,J,L

    if (trim(NAME) == "H2O") then
       call MAPL_GetPointer ( IMPORT,   XX,  'Q', RC=STATUS )
       VERIFY_(STATUS)
       ASSERT_(associated(XX))
    else
       call MAPL_GetPointer ( INTERNAL, XX, NAME, RC=STATUS )
       VERIFY_(STATUS)
       ASSERT_(associated(XX))
       call MAPL_GetResource(MAPL, VALUE,LABEL=trim(NAME)//"_FIXED_VALUE:", default=-1., RC=STATUS)
       VERIFY_(STATUS)
       if(VALUE>=0.0) then
          XX = VALUE
          return
       end if
    endif
 
    call MAPL_GetResource(MAPL,   TAU,LABEL=trim(NAME)//"_RELAXTIME:", DEFAULT=0.0 ,RC=STATUS)
    VERIFY_(STATUS)

! If there are multiple climYears, we are not allowing production and loss.
! -------------------------------------------------------------------------
    IF(PCHEM_STATE%climYears > 1 .AND. TAU <= 0.0) THEN
     IF(MAPL_AM_I_ROOT()) THEN
      PRINT *,TRIM(Iam)//": Cannot run PCHEM in P & L mode with PCHEM_STATE%climYears > 1."
      PRINT *,"            "//TRIM(NAME)//"_RELAXTIME has value ",TAU
     END IF
     STATUS = 1
     VERIFY_(STATUS)
    END IF

    call MAPL_GetPointer ( EXPORT, XX_PROD, trim(NAME)//'_PROD', RC=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetPointer ( EXPORT, XX_LOSS, trim(NAME)//'_LOSS', RC=STATUS )
    VERIFY_(STATUS)

    if (TAU<=0.0) then  ! By convention this is the prod(index 1) and loss(index 2) case

       ASSERT_(trim(NAME)/="H2O")

       PROD1 = PCHEM_STATE%MNPL(:,:,NN,1,1)*FAC + PCHEM_STATE%MNPL(:,:,NN,1,2)*(1.-FAC)
       LOSS1 = PCHEM_STATE%MNPL(:,:,NN,2,1)*FAC + PCHEM_STATE%MNPL(:,:,NN,2,2)*(1.-FAC)

       do j=1,jm
          do l=1,nlevs
             call MAPL_INTERP( PROD(:,L), LATS(:,J), Prod1(:,L), PCHEM_STATE%LATS)
             call MAPL_INTERP( LOSS(:,L), LATS(:,J), Loss1(:,L), PCHEM_STATE%LATS)
          enddo
          do i=1,im
             call MAPL_INTERP( PROD_INT(i,j,:), PL(i,j,:), PROD(i,:), PCHEM_STATE%LEVS)
             call MAPL_INTERP( LOSS_INT(i,j,:), PL(i,j,:), LOSS(i,:), PCHEM_STATE%LEVS)
          enddo
       end do

       XX = (XX + DT*PROD_INT) / (1.0 + DT*LOSS_INT)

    else ! If the relaxation time is positive, relax to climatology.

       PROD1 = PCHEM_STATE%MNCV(:,:,NN,1)*FAC + PCHEM_STATE%MNCV(:,:,NN,2)*(1.-FAC)

       do j=1,jm
          do l=1,nlevs
             call MAPL_INTERP( PROD(:,L), LATS(:,J), Prod1(:,L), PCHEM_STATE%LATS)
          enddo
          do i=1,im
             call MAPL_INTERP( PROD_INT(i,j,:), PL(i,j,:), PROD(i,:), PCHEM_STATE%LEVS)
          enddo
       end do

       call MAPL_GetResource(MAPL, DELP,  LABEL=trim(NAME)//"_DELP:" , DEFAULT=5000. ,RC=STATUS)
       VERIFY_(STATUS)

       DELP = max(DELP, 1.e-16) ! avoid division by zero

       if(trim(NAME)=="H2O") then
          call MAPL_GetResource(MAPL, PCRIT, LABEL=trim(NAME)//"_PCRIT:", DEFAULT=20000. ,RC=STATUS)
          VERIFY_(STATUS)
          allocate(WRK(IM,JM),stat=STATUS)
          VERIFY_(STATUS)
          where (TROPP==MAPL_UNDEF)
             WRK = PCRIT
          elsewhere
             WRK = TROPP
          end where
          WRK = min(WRK, PCRIT)
          do L=1,LM
             LOSS_INT(:,:,L) = (1./TAU) * max( min( (WRK-PL(:,:,L))/DELP, 1.0), 0.0)
          end do
          deallocate(WRK)
       else
          call MAPL_GetResource(MAPL, PCRIT, LABEL=trim(NAME)//"_PCRIT:", DEFAULT=1.e+16 ,RC=STATUS)
          VERIFY_(STATUS)
          LOSS_INT = (1./TAU) * max( min( (PCRIT   -PL)/DELP, 1.0), 0.0)
       endif

       PROD_INT = LOSS_INT*PROD_INT

       XX = (XX + DT*PROD_INT) / (1.0 + DT*LOSS_INT)

    end if


    if(associated(XX_PROD)) XX_PROD =  PROD_INT
    if(associated(XX_LOSS)) XX_LOSS = -LOSS_INT*XX

    if(trim(NAME)=='OX') then
       call MAPL_GetPointer ( EXPORT, OX_TEND, 'OX_TEND', RC=STATUS )
       VERIFY_(STATUS)
       if(associated(OX_TEND)) OX_TEND = (PROD_INT - LOSS_INT*XX)
    end if

    if(trim(NAME)=='H2O') then
       call MAPL_GetPointer ( EXPORT, H2O_TEND, 'H2O_TEND', RC=STATUS )
       VERIFY_(STATUS)
       if(associated(H2O_TEND)) H2O_TEND = (PROD_INT - LOSS_INT*XX)
    end if

    return
  end subroutine UPDATE

end subroutine RUN

!BOP

! !IROUTINE: AINC_UPDATE -- Update OX with analysis increment

! !INTERFACE:

  subroutine AINC_UPDATE ( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:

  type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
  type(ESMF_State),    intent(inout) :: IMPORT ! Import state
  type(ESMF_State),    intent(inout) :: EXPORT ! Export state
  type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
  integer, optional,   intent(  out) :: RC     ! Error code:

! ! DESCRIPTION: This RUN method simply updates OX with the analysis increment.

!EOP

! ErrLog Variables

  character(len=ESMF_MAXSTR)          :: IAm
  integer                             :: STATUS
  character(len=ESMF_MAXSTR)          :: COMP_NAME

! Local derived type aliases

  type (MAPL_MetaComp),     pointer   :: MAPL
  type (ESMF_State)                   :: INTERNAL
  type (MAPL_SunOrbit)                :: ORBIT
  type (ESMF_TimeInterval)            :: DELT
  real, pointer, dimension(:,:)       :: LONS
  real, pointer, dimension(:,:)       :: LATS
  real, allocatable, dimension(:,:)   :: ZTH, SLR, ZTHN

  real, parameter :: czaLimit = 1.0e-5
  integer                             :: IM, JM, LM
  integer                             :: L, SUNFLAG
  real, allocatable, dimension(:,:,:) :: pl
  real, allocatable, dimension(:,:,:) :: ro3ox

  real, pointer, dimension(:,:,:)     :: ple
  real, pointer, dimension(:,:,:)     :: do3
  real, pointer, dimension(:,:,:)     :: ox

  type(ESMF_Grid)                     :: grid

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

   Iam = "AINC_UPDATE"
   call ESMF_GridCompGet( GC, name=COMP_NAME, RC=STATUS )
   VERIFY_(STATUS)
   Iam = trim(COMP_NAME) // Iam

   if ( MAPL_AM_I_ROOT() ) then
       print *, 'Now running ',trim(Iam)
   endif

! Retrieve the pointer to the state
!----------------------------------

   call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
   VERIFY_(STATUS)

! Local aliases to the state, grid, and configuration
! ---------------------------------------------------

!  call MAPL_TimerOn(MAPL,"TOTAL")

! Get parameters from generic state.
!-----------------------------------

    call MAPL_Get( MAPL, IM=IM, JM=JM, LM=LM,    &
                   LONS=LONS, LATS=LATS,         &
                   ORBIT               = ORBIT,  &
                   INTERNAL_ESMF_STATE=INTERNAL, &
                                       RC=STATUS )
    VERIFY_(STATUS)

    call ESMF_GridCompGet(GC, grid=grid, rc=status)
    VERIFY_(STATUS)

    call ESMF_ClockGet(CLOCK, TIMESTEP=DELT, RC=STATUS)
    VERIFY_(STATUS)

! **********************************************************************
! ****               Get Pointers to BKG Import Data                ****
! **********************************************************************
#if 0
    if ( MAPL_AM_I_ROOT() ) then
       call ESMF_StatePrint(IMPORT)
    end if
#endif

    call MAPL_GetResource( MAPL, SUNFLAG, 'SUN_FLAG:', DEFAULT=0 , RC=STATUS)
    VERIFY_(STATUS)

!   Get pointers to import variables
!   --------------------------------
    call MAPL_GetPointer(import,   do3, 'O3AINC',  RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(import,   ple,    'PLE',  RC=STATUS)
    VERIFY_(STATUS)

!   Get pointers to internal variables
!   ----------------------------------
    call MAPL_GetPointer(internal,   ox, 'OX',  RC=STATUS)
    VERIFY_(STATUS)

    allocate ( ro3ox(IM,JM,LM) )
    allocate (    pl(IM,JM,LM) )

    allocate (   zth(IM,JM) )
    allocate (  zthn(IM,JM) )
    allocate (   slr(IM,JM) )

    call MAPL_SunGetInsolation(LONS, LATS,      &
                               ORBIT, ZTH, SLR, &
!                              INTV  = DELT,    &
                               CLOCK = CLOCK,   &
!                              TIME  = SUNFLAG, &
!                              ZTHN = ZTHN,     &
                               RC=STATUS )
    VERIFY_(STATUS)

    ZTH = max(ZTH,0.0)
    ro3ox = 1.0
    pl = 0.5 * ( ple(:,:,1:LM) + ple(:,:,2:LM+1) )
    do L = 1, LM
       where (pl(:,:,L) < 100.0 .and. ZTH > czaLimit)
         ro3ox(:,:,L) = exp(-1.5*(log10(0.01*pl(:,:,L)))**2)
       end where
    enddo
    deallocate (   slr )
    deallocate (  zthn )
    deallocate (   zth )
    deallocate (    pl )

    ox = max(0.0,ox*ro3ox*1.0e6 + do3)    ! update o3
    ox =         ox      *1.0e-6/ro3ox    ! convert updated o3 to ox

    deallocate ( ro3ox )

  end subroutine AINC_UPDATE


end module GEOS_PChemGridCompMod
