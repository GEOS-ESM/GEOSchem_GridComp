#include "MAPL_Generic.h"

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 910.1     !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: MAMchem_GridCompMod - Implements MAM Chemistry
!
! !INTERFACE:
!
   module MAMchem_GridCompMod
!
! !USES:
!
   use ESMF
   use MAPL_Mod

   use MAPL_SimpleBundleMod

   use Chem_UtilMod,        only: Chem_UtilResVal

   USE modal_aero_amicphys, only: modal_aero_amicphys_intr

   use modal_aero_data, only: numptr_amode, lmassptr_amode
   USE modal_aero_calcsize, only: modal_aero_calcsize_sub

   use MAM3_DataMod
   use MAM7_DataMod

   use MAM_BaseMod
   use MAM_ComponentsDataMod

   use MAM_SizeMod,          only: MAM_DrySize, MAM_WetSize

   use MAM_SeasaltMod,       only: MAM_SS_Emission, MAM_SS_Diagnostics
   use MAM_DustMod,          only: MAM_DU_Emission, MAM_DU_Diagnostics
   use MAM_BlackCarbonMod,   only: MAM_BC_Emission
   use MAM_OrganicCarbonMod, only: MAM_OC_Emission
   use MAM_SulfateMod,       only: MAM_SO4_Emission

   use MAM_DryRemovalMod,    only: MAM_DryRemoval
   use MAM_WetRemovalMod,    only: MAM_WetRemoval

   use MAML_OpticsTableMod,  only: MAML_OpticsTable,        &
                                   MAML_OpticsTableCreate,  &
                                   MAML_OpticsTableDestroy, &
                                   MAML_OpticsTableRead
   use MAML_OpticsMod,       only: MAML_OpticsInterpolate


   implicit none
   private


   type(MAML_OpticsTable), save :: MAM7_MieTable(7)

   integer, parameter :: instanceComputational = 1
   integer, parameter :: instanceData          = 2


!
! !PUBLIC MEMBER FUNCTIONS:

   public SetServices
!
! !DESCRIPTION: 
!
!  {\tt MAMchem\_GridComp} is an ESMF gridded component implementing
!  the MAM aerosol microphysical processes.
!
!  Developed for GEOS-5 release Fortuna 2.0 and later.
!
! !REVISION HISTORY:
!
!  06Dec2009  da Silva     Created the MATRIX skeleton.
!  15Aug2011  A. Darmenov  Initial version of MAM
!
!EOP
!-------------------------------------------------------------------------

!  Legacy state
!  ------------
   type MAM_State
      private

      logical                     :: data_driven = .false.

      type(ESMF_Config)           :: CF           ! Private Config

      type(ESMF_Grid)             :: grid         ! Grid

      integer                     :: im_world     ! Horizontal dimensions - lon
      integer                     :: jm_world     ! Horizontal dimensions - lat

      type(MAPL_SimpleBundle)     :: qa           ! Interstitial aerosol species and absorbed water
      type(MAPL_SimpleBundle)     :: qc           ! Cloud-borne  aerosol species
      type(MAPL_SimpleBundle)     :: qg           ! Gas species 

      type(MAPL_SimpleBundle)     :: Da           ! Dry and 'wet' geometric mean diameter of interstitial aerosol number size distribution

      real                        :: dt           ! Model time step

      integer                     :: scheme_id    ! MAM7 or MAM3
      type(MAM_Scheme)            :: scheme       ! MAM scheme/configuration

      real                        :: femisSS      ! Seasalt emission tuning parameter
      real                        :: femisDU      ! Dust emission tuning parameter

      real                        :: pom_oc_ratio ! ratio of POM emissions to primary OC emissions

      logical                     :: dry_removal  ! turn on/off dry removal processes
      logical                     :: wet_removal  ! turn on/off wet removal processes
      logical                     :: nucleation   ! turn on/off nucleation process
      logical                     :: condensation ! turn on/off condensation process
      logical                     :: coagulation  ! turn on/off coagulation process
      logical                     :: rename       ! turn on/off rename manager

      logical                     :: microphysics ! turn on/off aerosol microphysics
      logical                     :: mode_merging ! turn on/off explicit mode merging

      real                        :: f_conv_scav_ait
      real                        :: f_conv_scav_acc
      real                        :: f_conv_scav_pcm
      real                        :: f_conv_scav_fss
      real                        :: f_conv_scav_css
      real                        :: f_conv_scav_fdu
      real                        :: f_conv_scav_cdu

      real                        :: f_wet_ait
      real                        :: f_wet_acc
      real                        :: f_wet_pcm
      real                        :: f_wet_fss
      real                        :: f_wet_css
      real                        :: f_wet_fdu
      real                        :: f_wet_cdu

      type(MAML_OpticsTable)      :: mie_ait
      type(MAML_OpticsTable)      :: mie_acc
      type(MAML_OpticsTable)      :: mie_pcm
      type(MAML_OpticsTable)      :: mie_fss
      type(MAML_OpticsTable)      :: mie_css
      type(MAML_OpticsTable)      :: mie_fdu
      type(MAML_OpticsTable)      :: mie_cdu

      logical                     :: verbose      ! verbosity flag
   end type MAM_State

!  Hook for the ESMF
!  -----------------
   type MAM_Wrap
      type (MAM_State), pointer :: PTR => null()
   end type MAM_Wrap


contains


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 910.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: SetServices --- Sets IRF services for the MAMchem Grid Component
!
! !INTERFACE:

   subroutine SetServices ( GC, RC )

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
    type (MAM_State), pointer  :: self   ! internal private state
    type (MAM_wrap)            :: wrap

    character(len=ESMF_MAXSTR) :: comp_name

!   Local variables
!   --------------------------
    character(len=ESMF_MAXSTR) :: scheme            ! name of MAM scheme/configuration

    character(len=MAM_MAXSTR)  :: field_name        ! field name
    character(len=MAM_MAXSTR)  :: field_long_name   ! field name
    character(len=MAM_MAXSTR)  :: mode_name         ! aerosol mode name
    character(len=MAM_MAXSTR)  :: mode_long_name    ! aerosol mode name
    character(len=MAM_MAXSTR)  :: species_name      ! aerosol species name
    character(len=MAM_MAXSTR)  :: attachment_state  ! attachment state of aerosols
    integer                    :: n_species         ! number of aerosol species
    integer                    :: m, s              ! mode and species indexes

!   local
!   -----
    type(ESMF_Config)          :: CF                ! Universal Config 
    character(len=1024)        :: mie_optics_file   ! MAM Mie optics table file

    character(len=ESMF_MAXSTR), parameter:: microphysics_process(4) = (/'GAEX', 'RNAM', 'NUCL', 'COND'/)
    character(len=ESMF_MAXSTR) :: process           ! abbreviation of the process
    integer                    :: i                 ! counter 

!   Get my name and set-up traceback handle
!   ---------------------------------------
    call ESMF_GridCompGet(GC, name=comp_name, __RC__)
    Iam = trim(comp_name) // '::' // trim(Iam)

!   Wrap internal state for storing in GC; rename legacyState
!   -------------------------------------
    allocate ( self, stat=STATUS )
    _VERIFY(STATUS)
    wrap%ptr => self


!   Is the component data driven
!   ----------------------------
    self%data_driven = isDataDrivenGC_(GC, __RC__)


!   Load private Config Attributes
!   ------------------------------
    self%CF = ESMF_ConfigCreate(__RC__)
    call ESMF_ConfigLoadFile ( self%CF, 'MAMchem_GridComp.rc', __RC__ )

    call ESMF_ConfigGetAttribute ( self%CF, self%verbose,      label='verbose:',       default=.false.,  __RC__ )

    call ESMF_ConfigGetAttribute ( self%CF, scheme,            label='scheme:',        default='MAM7' ,  __RC__ )

    call ESMF_ConfigGetAttribute ( self%CF, self%dry_removal,  label='dry_removal:',   default=.true.,   __RC__ )
    call ESMF_ConfigGetAttribute ( self%CF, self%wet_removal,  label='wet_removal:',   default=.true.,   __RC__ )
    call ESMF_ConfigGetAttribute ( self%CF, self%nucleation,   label='nucleation:',    default=.true.,   __RC__ )
    call ESMF_ConfigGetAttribute ( self%CF, self%condensation, label='condensation:',  default=.true.,   __RC__ )
    call ESMF_ConfigGetAttribute ( self%CF, self%coagulation,  label='coagulation:',   default=.true.,   __RC__ )
    call ESMF_ConfigGetAttribute ( self%CF, self%rename,       label='rename:',        default=.true.,   __RC__ )

    call ESMF_ConfigGetAttribute ( self%CF, self%microphysics, label='microphysics:',  default=.true.,   __RC__ )
    call ESMF_ConfigGetAttribute ( self%CF, self%mode_merging, label='mode_merging:',  default=.true.,   __RC__ )

    
!   call ESMF_ConfigGetAttribute ( self%CF, self%femisSS,      label='seasalt_femis:', default=1.0,      __RC__ )
!   call ESMF_ConfigGetAttribute ( self%CF, self%femisDU,      label='dust_femis:',    default=1.0,      __RC__ )

    call ESMF_ConfigGetAttribute ( self%CF, self%pom_oc_ratio, label='pom_oc_ratio:',  default=1.4,      __RC__ )

    call ESMF_ConfigGetAttribute ( self%CF, self%f_conv_scav_ait, label='f_conv_scav_ait:', default=0.0, __RC__ )
    call ESMF_ConfigGetAttribute ( self%CF, self%f_conv_scav_acc, label='f_conv_scav_acc:', default=0.0, __RC__ )
    call ESMF_ConfigGetAttribute ( self%CF, self%f_conv_scav_pcm, label='f_conv_scav_pcm:', default=0.0, __RC__ )
    call ESMF_ConfigGetAttribute ( self%CF, self%f_conv_scav_fss, label='f_conv_scav_fss:', default=0.0, __RC__ )
    call ESMF_ConfigGetAttribute ( self%CF, self%f_conv_scav_css, label='f_conv_scav_css:', default=0.0, __RC__ )
    call ESMF_ConfigGetAttribute ( self%CF, self%f_conv_scav_fdu, label='f_conv_scav_fdu:', default=0.0, __RC__ )
    call ESMF_ConfigGetAttribute ( self%CF, self%f_conv_scav_cdu, label='f_conv_scav_cdu:', default=0.0, __RC__ )

    call ESMF_ConfigGetAttribute ( self%CF, self%f_wet_ait, label='f_wet_ait:', default=0.0, __RC__ )
    call ESMF_ConfigGetAttribute ( self%CF, self%f_wet_acc, label='f_wet_acc:', default=0.0, __RC__ )
    call ESMF_ConfigGetAttribute ( self%CF, self%f_wet_pcm, label='f_wet_pcm:', default=0.0, __RC__ )
    call ESMF_ConfigGetAttribute ( self%CF, self%f_wet_fss, label='f_wet_fss:', default=0.0, __RC__ )
    call ESMF_ConfigGetAttribute ( self%CF, self%f_wet_css, label='f_wet_css:', default=0.0, __RC__ )
    call ESMF_ConfigGetAttribute ( self%CF, self%f_wet_fdu, label='f_wet_fdu:', default=0.0, __RC__ )
    call ESMF_ConfigGetAttribute ( self%CF, self%f_wet_cdu, label='f_wet_cdu:', default=0.0, __RC__ )


!   Set the MAM model scheme
!   -------------------------
    scheme = ESMF_UtilStringUpperCase(scheme, __RC__)

    select case (scheme)
        case ('MAM7')
        self%scheme_id = MAM7_SCHEME

        case default
        __raise__ (MAM_UNKNOWN_SCHEME_ERROR, "Unsupported MAM scheme: " // trim(scheme))
    end select



!   Set the profiling timers
!   ------------------------
    call MAPL_TimerAdd(GC, name='TOTAL',                       __RC__)
    call MAPL_TimerAdd(GC, name='INITIALIZE',                  __RC__)
    call MAPL_TimerAdd(GC, name='RUN',                         __RC__)
    call MAPL_TimerAdd(GC, name='-EMISSIONS',                  __RC__)
    call MAPL_TimerAdd(GC, name='-MICROPHYSICS',               __RC__)
    call MAPL_TimerAdd(GC, name='--MICROPHYSICS_POSITIVE',     __RC__)
    call MAPL_TimerAdd(GC, name='-AQUEOUS_CHEM',               __RC__)
    call MAPL_TimerAdd(GC, name='-SIZE',                       __RC__)
    call MAPL_TimerAdd(GC, name='--SIZE_DRY',                  __RC__)
    call MAPL_TimerAdd(GC, name='--SIZE_WET',                  __RC__)
    call MAPL_TimerAdd(GC, name='-MODE_MERGING',               __RC__)
    call MAPL_TimerAdd(GC, name='-REMOVAL',                    __RC__)
    call MAPL_TimerAdd(GC, name='--REMOVAL_DRY',               __RC__)
    call MAPL_TimerAdd(GC, name='---REMOVAL_DRY_SETTLING',     __RC__)
    call MAPL_TimerAdd(GC, name='---REMOVAL_DRY_DEPOSITION',   __RC__)
    call MAPL_TimerAdd(GC, name='---REMOVAL_DRY_SOLVER',       __RC__)
    call MAPL_TimerAdd(GC, name='--REMOVAL_WET',               __RC__)
    call MAPL_TimerAdd(GC, name='-HYGROSCOPIC_GROWTH',         __RC__)
    call MAPL_TimerAdd(GC, name='-DIAGNOSTICS',                __RC__)
    call MAPL_TimerAdd(GC, name='--DIAGNOSTICS_SEASALT',       __RC__)
    call MAPL_TimerAdd(GC, name='--DIAGNOSTICS_DUST',          __RC__)
    call MAPL_TimerAdd(GC, name='--DIAGNOSTICS_CIM',           __RC__)
    call MAPL_TimerAdd(GC, name='--DIAGNOSTICS_SFC',           __RC__)
    call MAPL_TimerAdd(GC, name='--DIAGNOSTICS_AOT',           __RC__)


!                       ------------------------
!                       ESMF Functional Services
!                       ------------------------

!   Set the Initialize, Run, Finalize entry points
!   ----------------------------------------------
    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_INITIALIZE, Initialize_, __RC__ )
    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,        Run_,        __RC__ )
    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_FINALIZE,   Finalize_,   __RC__ )
        
!   Store internal state in GC
!   --------------------------
    call ESMF_UserCompSetInternalState ( GC, 'MAM_state', wrap, STATUS )
    _VERIFY(STATUS)
  
!                         ------------------
!                         MAPL Data Services
!                         ------------------

!BOS
!
! !IMPORT STATE:

#include "MAMchem_ImportSpec___.h"

! !INTERNAL STATE:


! !EXTERNAL STATE:

#include "MAMchem_ExportSpec___.h"

!EOS

!   Set MAM infrastructure
!   -----------------------------
    call MAM_SchemeInit(self%scheme, self%scheme_id, __RC__)


!   Add interstitial aerosols to the internal state
!   -----------------------------------------------
    attachment_state = 'interstitial'

    do m = 1, self%scheme%n_modes
        call MAM_AerosolModeGet(self%scheme%mode(m), name=mode_name, long_name=mode_long_name, n_species=n_species)

        ! number mixing ratio
        field_name = 'NUM_A_' // trim(mode_name)
        field_long_name = 'number of ' // trim(attachment_state) // ' aerosol particles in ' // trim(mode_long_name) // ' mode'

        call MAPL_AddInternalSpec(GC, SHORT_NAME = trim(field_name),      &
                                      LONG_NAME  = trim(field_long_name), &
                                      UNITS      = '#/kg',                &
                                      DIMS       = MAPL_DimsHorzVert,     &
                                      VLOCATION  = MAPL_VLocationCenter,  &
                                      RESTART    = MAPL_RestartOptional,  &
                                      FRIENDLYTO = 'DYNAMICS:TURBULENCE:MOIST', __RC__)

        ! mass mixing ratios
        do s = 1, n_species
            species_name = self%scheme%mode(m)%species(s)%name

            field_name = trim(species_name) // '_A_' // trim(mode_name)
            field_long_name = 'mass mixing ratio of ' // trim(attachment_state) // ' aerosol particles in ' // trim(mode_long_name) // ' mode'

            call MAPL_AddInternalSpec(GC, SHORT_NAME = trim(field_name),      &
                                          LONG_NAME  = trim(field_long_name), &
                                          UNITS      = 'kg kg-1',             &
                                          DIMS       = MAPL_DimsHorzVert,     &
                                          VLOCATION  = MAPL_VLocationCenter,  &
                                          RESTART    = MAPL_RestartOptional,  &
                                          FRIENDLYTO = 'DYNAMICS:TURBULENCE:MOIST', __RC__)
        end do

        ! absorbed water
        field_name = 'WTR_A_' // trim(mode_name)
        field_long_name = 'mass mixing ratio of absorbed water by ' // trim(attachment_state) // ' aerosol particles in ' // trim(mode_long_name) // ' mode'

        call MAPL_AddInternalSpec(GC, SHORT_NAME = trim(field_name),      &
                                      LONG_NAME  = trim(field_long_name), &
                                      UNITS      = 'kg kg-1',             &
                                      DIMS       = MAPL_DimsHorzVert,     &
                                      VLOCATION  = MAPL_VLocationCenter,  &
                                      ADD2EXPORT = .true.,                &
                                      RESTART    = MAPL_RestartSkip, __RC__)

        ! dry size
        field_name = 'DGN_DRY_' // trim(mode_name)
        field_long_name = 'dry diameter of ' // trim(attachment_state) // ' aerosol particles in ' // trim(mode_long_name) // ' mode'

        call MAPL_AddInternalSpec(GC, SHORT_NAME = trim(field_name),      &
                                      LONG_NAME  = trim(field_long_name), &
                                      UNITS      = 'm',                   &
                                      DIMS       = MAPL_DimsHorzVert,     &
                                      VLOCATION  = MAPL_VLocationCenter,  &
                                      ADD2EXPORT = .true.,                &
                                      RESTART    = MAPL_RestartSkip, __RC__)

        ! wet size
        field_name = 'DGN_WET_' // trim(mode_name)
        field_long_name = 'wet diameter of ' // trim(attachment_state) // ' aerosol particles in ' // trim(mode_long_name) // ' mode'

        call MAPL_AddInternalSpec(GC, SHORT_NAME = trim(field_name),      &
                                      LONG_NAME  = trim(field_long_name), &
                                      UNITS      = 'm',                   &
                                      DIMS       = MAPL_DimsHorzVert,     &
                                      VLOCATION  = MAPL_VLocationCenter,  &
                                      ADD2EXPORT = .true.,                &
                                      RESTART    = MAPL_RestartSkip, __RC__)
    end do


!   Add cloud-borne aerosols to the internal state
!   ----------------------------------------------
    attachment_state = 'cloud-borne'

    do m = 1, self%scheme%n_modes
        call MAM_AerosolModeGet(self%scheme%mode(m), name=mode_name, long_name=mode_long_name, n_species=n_species)

        ! number mixing ratio
        field_name = 'NUM_C_' // trim(mode_name)
        field_long_name = 'number of ' // trim(attachment_state) // ' aerosol particles in ' // trim(mode_long_name) // ' mode'

        call MAPL_AddInternalSpec(GC, SHORT_NAME = trim(field_name),      &
                                      LONG_NAME  = trim(field_long_name), &
                                      UNITS      = '#/kg',                &
                                      DIMS       = MAPL_DimsHorzVert,     &
                                      VLOCATION  = MAPL_VLocationCenter,  &
                                      RESTART    = MAPL_RestartOptional,  &
                                      FRIENDLYTO = 'MOIST', __RC__)

        ! mass mixing ratios
        do s = 1, n_species
            species_name = self%scheme%mode(m)%species(s)%name

            field_name  = trim(species_name) // '_C_' // trim(mode_name)
            field_long_name = 'mass mixing ratio of ' // trim(attachment_state) // ' aerosol particles in ' // trim(mode_long_name) // ' mode'

            call MAPL_AddInternalSpec(GC, SHORT_NAME = trim(field_name),      &
                                          LONG_NAME  = trim(field_long_name), &
                                          UNITS      = 'kg kg-1',             &
                                          DIMS       = MAPL_DimsHorzVert,     &
                                          VLOCATION  = MAPL_VLocationCenter,  &
                                          RESTART    = MAPL_RestartOptional,  &
                                          FRIENDLYTO = 'MOIST', __RC__)
        end do
    end do

#if (0)
!   Add dry and wet aerosol sizes to the internal state
!   ---------------------------------------------------
    attachment_state = 'interstitial'

    do m = 1, self%scheme%n_modes
        call MAM_AerosolModeGet(self%scheme%mode(m), name=mode_name, long_name=mode_long_name, n_species=n_species)
      
        field_name = 'DGN_DRY_' // trim(mode_name)
        field_long_name = 'dry geometric mean diameter of ' // trim(attachment_state) // ' aerosol particles in ' // trim(mode_long_name) // ' mode'

        call MAPL_AddInternalSpec(GC, SHORT_NAME = trim(field_name),      &
                                      LONG_NAME  = trim(field_long_name), &
                                      UNITS      = 'm',                   &
                                      DIMS       = MAPL_DimsHorzVert,     &
                                      VLOCATION  = MAPL_VLocationCenter,  &
                                      RESTART    = MAPL_RestartSkip,      &
                                      FRIENDLYTO = trim(COMP_NAME), __RC__)

        
        field_name = 'DGN_WET_' // trim(mode_name)
        field_long_name = 'wet size of ' // trim(attachment_state) // ' aerosol particles in ' // trim(mode_long_name) // ' mode'

        call MAPL_AddInternalSpec(GC, SHORT_NAME = trim(field_name),      &
                                      LONG_NAME  = trim(field_long_name), &
                                      UNITS      = 'm',                   &
                                      DIMS       = MAPL_DimsHorzVert,     &
                                      VLOCATION  = MAPL_VLocationCenter,  &
                                      RESTART    = MAPL_RestartSkip,      &
                                      FRIENDLYTO = trim(COMP_NAME), __RC__)
    end do
#endif

!   This state is needed by radiation - It will contain 
!   aerosol number and mass mixing ratios, and aerosol optics
!   ---------------------------------------------------------------
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'AERO',                      &
                                LONG_NAME  = 'aerosol_mixing_ratios',     &
                                UNITS      = 'kg kg-1',                   &
                                DIMS       = MAPL_DimsHorzVert,           &
                                VLOCATION  = MAPL_VLocationCenter,        &
                                DATATYPE   = MAPL_StateItem, __RC__)


!   This state is needed by MOIST - It will contain
!   aerosol number concentrations and aerosol activation properties
!   ---------------------------------------------------------------
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'AERO_ACI',                  &
                                LONG_NAME  = 'aerosol_cloud_interaction', &
                                UNITS      = 'kg kg-1',                   &
                                DIMS       = MAPL_DimsHorzVert,           &
                                VLOCATION  = MAPL_VLocationCenter,        &
                                DATATYPE   = MAPL_StateItem, __RC__)


!
!   Diagnostics: Column-integrated tendencies due to 
!   gas-aerosol-exchange/condensation, rename, nucleation and 
!   coagulation
!   ---------------------------------------------------------------
    MICROPHYSICS_PROCESSES: do i = 1, size(microphysics_process)

    process = trim(microphysics_process(i))
  
    attachment_state = 'interstitial'

    do m = 1, self%scheme%n_modes
        call MAM_AerosolModeGet(self%scheme%mode(m), name=mode_name, long_name=mode_long_name, n_species=n_species)

        ! number mixing ratio
        field_name = 'DDT_NUM_A_' // trim(mode_name) // '_' // trim(process)
        field_long_name = 'column-integrated_tendency_due_to_' // trim(process)

        call MAPL_AddExportSpec(GC, SHORT_NAME = trim(field_name),      &
                                    LONG_NAME  = trim(field_long_name), &
                                    UNITS      = '# m-2 s-1',           &
                                    DIMS       = MAPL_DimsHorzOnly,     &
                                    VLOCATION  = MAPL_VLocationNone,    &
                                    __RC__)

        ! mass mixing ratios
        do s = 1, n_species
            species_name = self%scheme%mode(m)%species(s)%name

            field_name = 'DDT_' // trim(species_name) // '_A_' // trim(mode_name) // '_' // trim(process)
            field_long_name = 'column-integrated_tendency_due_to_' // trim(process)

            call MAPL_AddExportSpec(GC, SHORT_NAME = trim(field_name),      &
                                        LONG_NAME  = trim(field_long_name), &
                                        UNITS      = 'kg m-2 s-1',          &
                                        DIMS       = MAPL_DimsHorzOnly,     &
                                        VLOCATION  = MAPL_VLocationNone,    &
                                        __RC__)
        end do
    end do


    CLOUD_BOURNE_RENAME_DIAGNOSTICS: if (process == 'RNAM') then

    attachment_state = 'cloud-borne'

    do m = 1, self%scheme%n_modes
        call MAM_AerosolModeGet(self%scheme%mode(m), name=mode_name, long_name=mode_long_name, n_species=n_species)

        ! number mixing ratio
        field_name = 'DDT_NUM_C_' // trim(mode_name) // '_' // trim(process)
        field_long_name = 'column-integrated_tendency_due_to_' // trim(process)

        call MAPL_AddExportSpec(GC, SHORT_NAME = trim(field_name),      &
                                    LONG_NAME  = trim(field_long_name), &
                                    UNITS      = '# m-2 s-1',           &
                                    DIMS       = MAPL_DimsHorzOnly,     &
                                    VLOCATION  = MAPL_VLocationNone,    &
                                    __RC__)

        ! mass mixing ratios
        do s = 1, n_species
            species_name = self%scheme%mode(m)%species(s)%name

            field_name = 'DDT_' // trim(species_name) // '_C_' // trim(mode_name) // '_' // trim(process)
            field_long_name = 'column-integrated_tendency_due_to_' // trim(process)

            call MAPL_AddExportSpec(GC, SHORT_NAME = trim(field_name),      &
                                        LONG_NAME  = trim(field_long_name), &
                                        UNITS      = 'kg m-2 s-1',          &
                                        DIMS       = MAPL_DimsHorzOnly,     &
                                        VLOCATION  = MAPL_VLocationNone,    &
                                        __RC__)
        end do
    end do
    end if CLOUD_BOURNE_RENAME_DIAGNOSTICS

    end do MICROPHYSICS_PROCESSES
    

    
!   Create Mie tables for coupling with radiation
!   ---------------------------------------------
    call ESMF_GridCompGet(GC, config=CF, __RC__) ! read paths to Mie tables from global config

    call ESMF_ConfigGetAttribute(CF, mie_optics_file, label='MAM7_AIT_OPTICS:', __RC__)
    MAM7_MieTable(1) = MAML_OpticsTableCreate(trim(mie_optics_file), __RC__)
    call MAML_OpticsTableRead(MAM7_MieTable(1), __RC__)

    call ESMF_ConfigGetAttribute(CF, mie_optics_file, label='MAM7_ACC_OPTICS:', __RC__)
    MAM7_MieTable(2) = MAML_OpticsTableCreate(trim(mie_optics_file), __RC__)
    call MAML_OpticsTableRead(MAM7_MieTable(2), __RC__)

    call ESMF_ConfigGetAttribute(CF, mie_optics_file, label='MAM7_PCM_OPTICS:', __RC__)
    MAM7_MieTable(3) = MAML_OpticsTableCreate(trim(mie_optics_file), __RC__)
    call MAML_OpticsTableRead(MAM7_MieTable(3), __RC__)

    call ESMF_ConfigGetAttribute(CF, mie_optics_file, label='MAM7_FSS_OPTICS:', __RC__)
    MAM7_MieTable(4) = MAML_OpticsTableCreate(trim(mie_optics_file), __RC__)
    call MAML_OpticsTableRead(MAM7_MieTable(4), __RC__)

    call ESMF_ConfigGetAttribute(CF, mie_optics_file, label='MAM7_CSS_OPTICS:', __RC__)
    MAM7_MieTable(5) = MAML_OpticsTableCreate(trim(mie_optics_file), __RC__)
    call MAML_OpticsTableRead(MAM7_MieTable(5), __RC__)

    call ESMF_ConfigGetAttribute(CF, mie_optics_file, label='MAM7_FDU_OPTICS:', __RC__)
    MAM7_MieTable(6) = MAML_OpticsTableCreate(trim(mie_optics_file), __RC__)
    call MAML_OpticsTableRead(MAM7_MieTable(6), __RC__)

    call ESMF_ConfigGetAttribute(CF, mie_optics_file, label='MAM7_CDU_OPTICS:', __RC__)
    MAM7_MieTable(7) = MAML_OpticsTableCreate(trim(mie_optics_file), __RC__)
    call MAML_OpticsTableRead(MAM7_MieTable(7), __RC__)


!   This bundle is not filled in by MAM, just a place holder for now
!   ----------------------------------------------------------------
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'AERO_DP',              &
                                LONG_NAME  = 'aerosol_deposition',   &
                                UNITS      = 'kg m-2 s-1',           &
                                DIMS       = MAPL_DimsHorzOnly,      &
                                DATATYPE   = MAPL_BundleItem,  __RC__)

!   Generic Set Services
!   --------------------
    call MAPL_GenericSetServices ( GC, __RC__ )

!   All done
!   --------

    _RETURN(ESMF_SUCCESS)

  end subroutine SetServices


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Initialize_ --- Initialize MAMchem
!
! !INTERFACE:
!

   subroutine Initialize_ ( GC, IMPORT, EXPORT, CLOCK, rc )

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
!  01Dec2009 da Silva  First crack.
!
!EOP
!-------------------------------------------------------------------------

                              __Iam__('Initialize_')

    type(MAM_state), pointer      :: self               ! Legacy state
    type(ESMF_Grid)               :: GRID               ! Grid
    type(ESMF_Config)             :: CF                 ! Universal Config 

    type(ESMF_State)              :: aero               ! 
    type(ESMF_FieldBundle)        :: aero_state_aerosols!
    logical                       :: implements_aerosol_optics
    type(ESMF_Field)              :: field              ! field
    character(len=MAM_MAXSTR)     :: field_name         ! field name

    type(ESMF_State)              :: aero_aci           !
    logical                       :: implements_aap_method
    character(len=ESMF_MAXSTR), allocatable, dimension(:) :: aero_aci_modes

    integer                       :: im_World, jm_World ! Global 2D Dimensions
    integer                       :: im, jm, lm         ! 3D Dimensions

    integer, parameter            :: n_hres = 6         ! number of horizontal resolutions (a, b, c, d, e)
    real, dimension(n_hres)       :: f_hres             ! buffer for the resolution dependent factors
    integer                       :: n


    integer                       :: nymd, nhms         ! date, time
    real                          :: cdt                ! time step in secs

    character(len=1024)           :: var_names          ! comma separated names

    character(len=ESMF_MAXSTR)    :: comp_name          ! component's name

    character(len=MAM_MAXSTR)     :: mode_name          ! aerosol mode name
    character(len=MAM_MAXSTR)     :: species_name       ! aerosol species name
    integer                       :: n_species          ! number of aerosol species
    integer                       :: m, s               ! mode and species indexes

    real                          :: hygroscopicity     ! hygroscopicity of the aerosol species 
    real                          :: solubility         ! solubility of the aerosol specie

    real                          :: f_scav             ! globally uniform convective scavenging coefficient, km-1

    character(len=1024)           :: mie_optics_file

    real, parameter               :: f_scav_none = 0.0
    real, parameter               :: f_wet_none  = 0.0

    character(len=2), parameter   :: name_delimiter = ', '


!  Declare pointers to IMPORT/EXPORT/INTERNAL states 
!  -------------------------------------------------
    type(MAPL_MetaComp), pointer  :: mgState
    type(ESMF_State)              :: INTERNAL
  
!   Get my name and set-up traceback handle
!   ---------------------------------------
    call ESMF_GridCompGet( GC, name=COMP_NAME, __RC__ )
    Iam = trim(comp_name) // '::' // trim(Iam)

!                               --------
    if (MAPL_AM_I_ROOT()) then
       write (*,*) trim(Iam)//': Starting...'
       write (*,*)
    end if

    call MAPL_GetObjectFromGC ( GC, mgState, __RC__)

    call MAPL_TimerOn(mgState, 'TOTAL',      __RC__)
    call MAPL_TimerOn(mgState, 'INITIALIZE', __RC__)

!   Initialize MAPL Generic
!   -----------------------
    call MAPL_GenericInitialize ( GC, IMPORT, EXPORT, clock,  __RC__ )

!   Get pointers to IMPORT/EXPORT/INTERNAL states 
!   ---------------------------------------------
    call MAPL_Get (mgState, INTERNAL_ESMF_STATE=INTERNAL, __RC__)

!   Extract relevant runtime information
!   ------------------------------------
    call extract_ ( GC, CLOCK, self, GRID, CF, &
                    im_World, jm_World,        &
                    im, jm, lm,                &
                    nymd, nhms, cdt, __RC__ )

!   Aerosol Processes
!   -----------------
    if (self%verbose .and. MAPL_AM_I_ROOT()) then
       write (*,*) 'Aerosol processes:'
       call PrintProcessFlag('Dry removal ', self%dry_removal)
       call PrintProcessFlag('Wet removal ', self%wet_removal)

       if (.not. self%microphysics) then
           call PrintProcessFlag('Microphysics', self%microphysics)
       else
           call PrintProcessFlag('Condensation', self%condensation)
           call PrintProcessFlag('Nucleation  ', self%nucleation)
           call PrintProcessFlag('Coagulation ', self%coagulation)
           call PrintProcessFlag('Rename      ', self%rename)
       end if

       call PrintProcessFlag('Mode merging', self%mode_merging)

       write (*,*)
    end if


!   Set the grid and dimensions
!   ---------------------------
    self%grid = GRID

    self%im_world = im_World
    self%jm_world = jm_World


!   Set the time step
!   -----------------
    self%dt = cdt


!   Set resolution dependent parameters
!   -----------------------------------
    call ESMF_ConfigFindLabel(self%CF, 'seasalt_femis:', __RC__)
    do n = 1, size(f_hres)
        call ESMF_ConfigGetAttribute(self%CF, f_hres(n), __RC__)
    end do
    self%femisSS = Chem_UtilResVal(self%im_world, self%jm_world, f_hres(:), STATUS)
    _VERIFY(STATUS)

   
    call ESMF_ConfigFindLabel(self%CF, 'dust_femis:', __RC__)
    do n = 1, size(f_hres)
        call ESMF_ConfigGetAttribute(self%CF, f_hres(n), __RC__)
    end do
    self%femisDU = Chem_UtilResVal(self%im_world, self%jm_world, f_hres(:), STATUS)
    _VERIFY(STATUS)


!   Box model aerosol microphysics
!   ------------------------------
    call microphysics_initialize(imozart=1, verbose=.false., __RC__)


!   Bundle the interstitial and cloud-borne tracers
!   -------------------------------------------------
    var_names = ''
    do m = 1, self%scheme%n_modes
        call MAM_AerosolModeGet(self%scheme%mode(m), name=mode_name, n_species=n_species)

        ! number mixing ratio
        field_name = 'NUM_A_' // trim(mode_name)
        var_names  = trim(var_names) // trim(field_name) // trim(name_delimiter)

        ! absorbed water
        field_name = 'WTR_A_' // trim(mode_name)
        var_names  = trim(var_names) // trim(field_name) // trim(name_delimiter)

        ! mass mixing ratios
        do s = 1, n_species
            species_name = self%scheme%mode(m)%species(s)%name

            field_name  = trim(species_name) // '_A_' // trim(mode_name)
            var_names   = trim(var_names) // trim(field_name) // trim(name_delimiter)
        end do
    end do
    var_names = trim(var_names(1:len_trim(var_names)-1))

    self%qa = MAPL_SimpleBundleCreate(INTERNAL, name='MAM_INTERSTITIAL_AEROSOLS', &
                                                only_vars=var_names, __RC__)
    call MAPL_SimpleBundlePrint(self%qa)


    var_names = ''
    do m = 1, self%scheme%n_modes
        call MAM_AerosolModeGet(self%scheme%mode(m), name=mode_name, n_species=n_species)

        ! number mixing ratio
        field_name = 'NUM_C_' // trim(mode_name)
        var_names  = trim(var_names) // trim(field_name) // trim(name_delimiter)

        ! mass mixing ratios
        do s = 1, n_species
            species_name = self%scheme%mode(m)%species(s)%name

            field_name  = trim(species_name) // '_C_' // trim(mode_name)
            var_names   = trim(var_names) // trim(field_name) // trim(name_delimiter)
        end do
    end do
    var_names = trim(var_names(1:len_trim(var_names)-1))

    self%qc = MAPL_SimpleBundleCreate(INTERNAL, name='MAM_CLOUDBORNE_AEROSOLS', &
                                                only_vars=var_names, __RC__)
    call MAPL_SimpleBundlePrint(self%qc)


!   Bundle the dry and wet modal geometric diameters
!   -------------------------------------------------
    var_names = ''
    do m = 1, self%scheme%n_modes
        call MAM_AerosolModeGet(self%scheme%mode(m), name=mode_name, n_species=n_species)

        field_name = 'DGN_DRY_' // trim(mode_name)
        var_names  = trim(var_names) // trim(field_name) // trim(name_delimiter)

        field_name = 'DGN_WET_' // trim(mode_name)
        var_names  = trim(var_names) // trim(field_name) // trim(name_delimiter)
    end do
    var_names = trim(var_names(1:len_trim(var_names)-1))

    self%Da   = MAPL_SimpleBundleCreate(INTERNAL, name='MAM_INTERSTITIAL_AEROSOLS_DIAMETERS', &
                                                 only_vars=var_names, __RC__)
    call MAPL_SimpleBundlePrint(self%Da)

!   Bundle the gas tracers
!   -------------------------------------------------
    var_names = 'H2SO4, SO2, NH3, SOA_GAS'

    self%qg = MAPL_SimpleBundleCreate(IMPORT, name='MAM_GAS_SPECIES', &
                                              only_vars=var_names, __RC__)
    call MAPL_SimpleBundlePrint(self%qg)


!   Fill the AERO State with the aerosol mixing ratios
!   ---------------------------------------------------
    call ESMF_StateGet(EXPORT, 'AERO', aero, __RC__)

    ! This attribute indicates if the aerosol optics method is implemented or not. 
    ! Radiation will not call the aerosol optics method unless this attribute is 
    ! explicitly set to true.

    implements_aerosol_optics = .true. 

    call ESMF_AttributeSet(aero, name  = 'implements_aerosol_optics_method', & 
                                 value = implements_aerosol_optics, __RC__)
  
    COUPLING_TO_RADIATION: if (implements_aerosol_optics) then

        aero_state_aerosols = ESMF_FieldBundleCreate(name="AEROSOLS", __RC__)
        call MAPL_StateAdd(aero, aero_state_aerosols, __RC__)

        do m = 1, self%scheme%n_modes
            call MAM_AerosolModeGet(self%scheme%mode(m), name=mode_name, n_species=n_species)

            ! interstitial aerosol tracers
            field_name = 'NUM_A_' // trim(mode_name)
            call ESMF_StateGet(INTERNAL, trim(field_name), field, __RC__)
            call ESMF_FieldBundleAdd(aero_state_aerosols, (/field/), __RC__)

            field_name = 'WTR_A_' // trim(mode_name)
            call ESMF_StateGet(INTERNAL, trim(field_name), field, __RC__)
            call ESMF_FieldBundleAdd(aero_state_aerosols, (/field/), __RC__)

            field_name = 'DGN_WET_' // trim(mode_name)
            call ESMF_StateGet(INTERNAL, trim(field_name), field, __RC__)
            call ESMF_FieldBundleAdd(aero_state_aerosols, (/field/), __RC__)

            do s = 1, n_species
                species_name = self%scheme%mode(m)%species(s)%name
                field_name  = trim(species_name) // '_A_' // trim(mode_name)
           
                call ESMF_StateGet(INTERNAL, trim(field_name), field, __RC__)
                call ESMF_FieldBundleAdd(aero_state_aerosols, (/field/), __RC__)
            end do
        end do 


        ! TODO
        ! MAM_MieTable = MAM_MieCreate(CF, __RC__)

        ! state of the atmosphere
        call ESMF_AttributeSet(aero, name='air_pressure_for_aerosol_optics',             value='PLE', __RC__)
        call ESMF_AttributeSet(aero, name='relative_humidity_for_aerosol_optics',        value='RH',  __RC__)
        call ESMF_AttributeSet(aero, name='cloud_area_fraction_for_aerosol_optics',      value='',    __RC__) ! 'cloud_area_fraction_in_atmosphere_layer_for_aerosol_optics'

        ! aerosol optics
        call ESMF_AttributeSet(aero, name='band_for_aerosol_optics',                     value=0,     __RC__)
        call ESMF_AttributeSet(aero, name='extinction_in_air_due_to_ambient_aerosol',    value='EXT', __RC__)
        call ESMF_AttributeSet(aero, name='single_scattering_albedo_of_ambient_aerosol', value='SSA', __RC__)
        call ESMF_AttributeSet(aero, name='asymmetry_parameter_of_ambient_aerosol',      value='ASY', __RC__)

        ! add PLE to Aero state
        call ESMF_AttributeGet(aero, name='air_pressure_for_aerosol_optics', value=field_name, __RC__)
        if (field_name /= '') then
            field = MAPL_FieldCreateEmpty(trim(field_name), self%grid, __RC__)

            call MAPL_FieldAllocCommit(field, dims=MAPL_DimsHorzVert, location=MAPL_VLocationEdge, typekind=MAPL_R4, hw=0, __RC__)
            call MAPL_StateAdd(aero, field, __RC__)
        end if

        ! add RH to Aero state
        call ESMF_AttributeGet(aero, name='relative_humidity_for_aerosol_optics', value=field_name, __RC__)
        if (field_name /= '') then
            field = MAPL_FieldCreateEmpty(trim(field_name), self%grid, __RC__)

            call MAPL_FieldAllocCommit(field, dims=MAPL_DimsHorzVert, location=MAPL_VLocationCenter, typekind=MAPL_R4, hw=0, __RC__)
            call MAPL_StateAdd(aero, field, __RC__)
        end if

        ! add EXT to Aero state
        call ESMF_AttributeGet(aero, name='extinction_in_air_due_to_ambient_aerosol', value=field_name, __RC__)
        if (field_name /= '') then 
            field = MAPL_FieldCreateEmpty(trim(field_name), self%grid, __RC__)

            call MAPL_FieldAllocCommit(field, dims=MAPL_DimsHorzVert, location=MAPL_VLocationCenter, typekind=MAPL_R4, hw=0, __RC__)
            call MAPL_StateAdd(aero, field, __RC__)
        end if

        ! add SSA to aero state
        call ESMF_AttributeGet(aero, name='single_scattering_albedo_of_ambient_aerosol', value=field_name, __RC__)
        if (field_name /= '') then
            field = MAPL_FieldCreateEmpty(trim(field_name), self%grid, __RC__)

            call MAPL_FieldAllocCommit(field, dims=MAPL_DimsHorzVert, location=MAPL_VLocationCenter, typekind=MAPL_R4, hw=0, __RC__)
            call MAPL_StateAdd(aero, field, __RC__)
        end if

        ! add ASY to aero state
        call ESMF_AttributeGet(aero, name='asymmetry_parameter_of_ambient_aerosol', value=field_name, RC=STATUS)
        if (field_name /= '') then
            field = MAPL_FieldCreateEmpty(trim(field_name), self%grid, __RC__)

            call MAPL_FieldAllocCommit(field, dims=MAPL_DimsHorzVert, location=MAPL_VLocationCenter, typekind=MAPL_R4, hw=0, __RC__)
            call MAPL_StateAdd(aero, field, __RC__)
        end if
       
        ! attach the aerosol optics method
        call ESMF_MethodAdd(aero, label='aerosol_optics', userRoutine=aerosol_optics, __RC__)

    end if COUPLING_TO_RADIATION


!   Fill the AERO State with the aerosol mixing ratios
!   ---------------------------------------------------
    call ESMF_StateGet(EXPORT, 'AERO_ACI', aero_aci, __RC__)

    ! This attribute indicates if the aerosol optics method is implemented or not. 
    ! Radiation will not call the aerosol optics method unless this attribute is 
    ! explicitly set to true.

    implements_aap_method = .true. 

    call ESMF_AttributeSet(aero_aci, name  = 'implements_aerosol_activation_properties_method', & 
                                     value = implements_aap_method, __RC__)
  
    COUPLING_TO_CLOUD_MICROPHYSICS: if (implements_aap_method) then

        _ASSERT(self%scheme%n_modes > 0,'needs informative message')

        allocate(aero_aci_modes(self%scheme%n_modes), __STAT__)

        aero_state_aerosols = ESMF_FieldBundleCreate(name="AEROSOLS", __RC__)
        call MAPL_StateAdd(aero_aci, aero_state_aerosols, __RC__)

        do m = 1, self%scheme%n_modes
            call MAM_AerosolModeGet(self%scheme%mode(m), name=mode_name, n_species=n_species)

            ! add the mode name to the list of aerosol modes
            aero_aci_modes(m) = trim(mode_name)

            ! interstitial aerosol tracers
            field_name = 'NUM_A_' // trim(mode_name)
            call ESMF_StateGet(INTERNAL, trim(field_name), field, __RC__)
            call ESMF_FieldBundleAdd(aero_state_aerosols, (/field/), __RC__)

            do s = 1, n_species
                species_name = self%scheme%mode(m)%species(s)%name
                field_name  = trim(species_name) // '_A_' // trim(mode_name)
           
                call ESMF_StateGet(INTERNAL, trim(field_name), field, __RC__)
                call ESMF_FieldBundleAdd(aero_state_aerosols, (/field/), __RC__)
            end do
        end do

        ! Following the aerosol-cloud-interaction state protocol, next steps are:  
        !  - attach a list with the aerosol modes
        !  - attach required met fields
        !  - attach method that computes the aerosol activation properties
    
        call ESMF_AttributeSet(aero_aci, name='number_of_aerosol_modes', value=self%scheme%n_modes, __RC__)
        call ESMF_AttributeSet(aero_aci, name='aerosol_modes', itemcount=self%scheme%n_modes, valuelist=aero_aci_modes, __RC__)

        deallocate(aero_aci_modes, __STAT__)


        ! met fields and land fraction
        call ESMF_AttributeSet(aero_aci, name='air_pressure',                 value='',         __RC__)
        call ESMF_AttributeSet(aero_aci, name='air_temperature',              value='',         __RC__)
        call ESMF_AttributeSet(aero_aci, name='fraction_of_land_type',        value='',         __RC__)

        ! aerosol activation properties
        call ESMF_AttributeSet(aero_aci, name='width_of_aerosol_mode',        value='SIGMA',    __RC__)
        call ESMF_AttributeSet(aero_aci, name='aerosol_number_concentration', value='NUM',      __RC__)
        call ESMF_AttributeSet(aero_aci, name='aerosol_dry_size',             value='DGN',      __RC__)
        call ESMF_AttributeSet(aero_aci, name='aerosol_density',              value='density',  __RC__)
        call ESMF_AttributeSet(aero_aci, name='aerosol_hygroscopicity',       value='KAPPA',    __RC__)
        call ESMF_AttributeSet(aero_aci, name='fraction_of_dust_aerosol',     value='FDUST',    __RC__)
        call ESMF_AttributeSet(aero_aci, name='fraction_of_soot_aerosol',     value='FSOOT',    __RC__)
        call ESMF_AttributeSet(aero_aci, name='fraction_of_organic_aerosol',  value='FORGANIC', __RC__)


        ! add PLE to ACI state
        call ESMF_AttributeGet(aero_aci, name='air_pressure', value=field_name, __RC__)
        if (field_name /= '') then
            field = MAPL_FieldCreateEmpty(trim(field_name), self%grid, __RC__)

            call MAPL_FieldAllocCommit(field, dims=MAPL_DimsHorzVert, location=MAPL_VLocationEdge, typekind=MAPL_R4, hw=0, __RC__)
            call MAPL_StateAdd(aero_aci, field, __RC__)
        end if

        ! add T to ACI state
        call ESMF_AttributeGet(aero_aci, name='air_temperature', value=field_name, __RC__)
        if (field_name /= '') then
            field = MAPL_FieldCreateEmpty(trim(field_name), self%grid, __RC__)

            call MAPL_FieldAllocCommit(field, dims=MAPL_DimsHorzVert, location=MAPL_VLocationCenter, typekind=MAPL_R4, hw=0, __RC__)
            call MAPL_StateAdd(aero_aci, field, __RC__)
        end if

        ! add FRLAND to ACI state
        call ESMF_AttributeGet(aero_aci, name='fraction_of_land_type', value=field_name, __RC__)
        if (field_name /= '') then
            field = MAPL_FieldCreateEmpty(trim(field_name), self%grid, __RC__)

            call MAPL_FieldAllocCommit(field, dims=MAPL_DimsHorzOnly, location=MAPL_VLocationCenter, typekind=MAPL_R4, hw=0, __RC__)
            call MAPL_StateAdd(aero_aci, field, __RC__)
        end if

        ! add aerosol activation properties to ACI state
        call ESMF_AttributeGet(aero_aci, name='width_of_aerosol_mode', value=field_name, __RC__)
        if (field_name /= '') then
            field = MAPL_FieldCreateEmpty(trim(field_name), self%grid, __RC__)

            call MAPL_FieldAllocCommit(field, dims=MAPL_DimsHorzVert, location=MAPL_VLocationCenter, typekind=MAPL_R4, hw=0, __RC__)
            call MAPL_StateAdd(aero_aci, field, __RC__)
        end if

        call ESMF_AttributeGet(aero_aci, name='aerosol_number_concentration', value=field_name, __RC__)
        if (field_name /= '') then
            field = MAPL_FieldCreateEmpty(trim(field_name), self%grid, __RC__)

            call MAPL_FieldAllocCommit(field, dims=MAPL_DimsHorzVert, location=MAPL_VLocationCenter, typekind=MAPL_R4, hw=0, __RC__)
            call MAPL_StateAdd(aero_aci, field, __RC__)
        end if

        call ESMF_AttributeGet(aero_aci, name='aerosol_dry_size', value=field_name, __RC__)
        if (field_name /= '') then
            field = MAPL_FieldCreateEmpty(trim(field_name), self%grid, __RC__)

            call MAPL_FieldAllocCommit(field, dims=MAPL_DimsHorzVert, location=MAPL_VLocationCenter, typekind=MAPL_R4, hw=0, __RC__)
            call MAPL_StateAdd(aero_aci, field, __RC__)
        end if

        call ESMF_AttributeGet(aero_aci, name='aerosol_density', value=field_name, __RC__)
        if (field_name /= '') then
            field = MAPL_FieldCreateEmpty(trim(field_name), self%grid, __RC__)

            call MAPL_FieldAllocCommit(field, dims=MAPL_DimsHorzVert, location=MAPL_VLocationCenter, typekind=MAPL_R4, hw=0, __RC__)
            call MAPL_StateAdd(aero_aci, field, __RC__)
        end if

        call ESMF_AttributeGet(aero_aci, name='aerosol_hygroscopicity', value=field_name, __RC__)
        if (field_name /= '') then
            field = MAPL_FieldCreateEmpty(trim(field_name), self%grid, __RC__)

            call MAPL_FieldAllocCommit(field, dims=MAPL_DimsHorzVert, location=MAPL_VLocationCenter, typekind=MAPL_R4, hw=0, __RC__)
            call MAPL_StateAdd(aero_aci, field, __RC__)
        end if

        call ESMF_AttributeGet(aero_aci, name='fraction_of_dust_aerosol', value=field_name, __RC__)
        if (field_name /= '') then
            field = MAPL_FieldCreateEmpty(trim(field_name), self%grid, __RC__)

            call MAPL_FieldAllocCommit(field, dims=MAPL_DimsHorzVert, location=MAPL_VLocationCenter, typekind=MAPL_R4, hw=0, __RC__)
            call MAPL_StateAdd(aero_aci, field, __RC__)
        end if

        call ESMF_AttributeGet(aero_aci, name='fraction_of_soot_aerosol', value=field_name, __RC__)
        if (field_name /= '') then
            field = MAPL_FieldCreateEmpty(trim(field_name), self%grid, __RC__)

            call MAPL_FieldAllocCommit(field, dims=MAPL_DimsHorzVert, location=MAPL_VLocationCenter, typekind=MAPL_R4, hw=0, __RC__)
            call MAPL_StateAdd(aero_aci, field, __RC__)
        end if

        call ESMF_AttributeGet(aero_aci, name='fraction_of_organic_aerosol', value=field_name, __RC__)
        if (field_name /= '') then
            field = MAPL_FieldCreateEmpty(trim(field_name), self%grid, __RC__)

            call MAPL_FieldAllocCommit(field, dims=MAPL_DimsHorzVert, location=MAPL_VLocationCenter, typekind=MAPL_R4, hw=0, __RC__)
            call MAPL_StateAdd(aero_aci, field, __RC__)
        end if

        ! attach the aerosol optics method
        call ESMF_MethodAdd(aero_aci, label='aerosol_activation_properties', userRoutine=aerosol_activation_properties, __RC__)
     
    end if COUPLING_TO_CLOUD_MICROPHYSICS



!   Fill the scavenging coefficients
!   --------------------------------
    do m = 1, self%scheme%n_modes
        call MAM_AerosolModeGet(self%scheme%mode(m), name=mode_name, n_species=n_species)

        select case (trim(mode_name))
            case (trim(MAM7_AITKEN_MODE_NAME))
                call MAM_AerosolModeSet(self%scheme%mode(m), f_conv_scavenging = self%f_conv_scav_ait, &
                                                             f_wet             = self%f_wet_ait)

            case (trim(MAM7_ACCUMULATION_MODE_NAME))
                call MAM_AerosolModeSet(self%scheme%mode(m), f_conv_scavenging = self%f_conv_scav_acc, &
                                                             f_wet             = self%f_wet_acc)

            case (trim(MAM7_PRIMARY_CARBON_MODE_NAME))
                call MAM_AerosolModeSet(self%scheme%mode(m), f_conv_scavenging = self%f_conv_scav_pcm, &
                                                             f_wet             = self%f_wet_pcm)

            case (trim(MAM7_FINE_SEASALT_MODE_NAME))
                call MAM_AerosolModeSet(self%scheme%mode(m), f_conv_scavenging = self%f_conv_scav_fss, &
                                                             f_wet             = self%f_wet_fss)

            case (trim(MAM7_FINE_DUST_MODE_NAME))
                call MAM_AerosolModeSet(self%scheme%mode(m), f_conv_scavenging = self%f_conv_scav_fdu, &
                                                             f_wet             = self%f_wet_fdu)

            case (trim(MAM7_COARSE_SEASALT_MODE_NAME))
                call MAM_AerosolModeSet(self%scheme%mode(m), f_conv_scavenging = self%f_conv_scav_css, &
                                                             f_wet             = self%f_wet_css)

            case (trim(MAM7_COARSE_DUST_MODE_NAME))
                call MAM_AerosolModeSet(self%scheme%mode(m), f_conv_scavenging = self%f_conv_scav_cdu, &
                                                             f_wet             = self%f_wet_cdu)

            case default
                call MAM_AerosolModeSet(self%scheme%mode(m), f_conv_scavenging = f_scav_none, f_wet = f_wet_none)
        end select
    end do

    do m = 1, self%scheme%n_modes
        call MAM_AerosolModeGet(self%scheme%mode(m), name=mode_name, n_species=n_species, f_conv_scavenging=f_scav)

        if (self%verbose .and. MAPL_AM_I_ROOT()) then
            print *, trim(Iam)//': Convective Scavenging Parameter of '//trim(mode_name)//' : ', f_scav
        end if
 
        ! interstitial aerosol tracers
        field_name = 'NUM_A_' // trim(mode_name)
        call ESMF_StateGet(INTERNAL, trim(field_name), field, __RC__)

        call ESMF_AttributeSet(field, NAME='ScavengingFractionPerKm', VALUE=f_scav, __RC__)

        do s = 1, n_species
            species_name = self%scheme%mode(m)%species(s)%name
            field_name  = trim(species_name) // '_A_' // trim(mode_name)
            call ESMF_StateGet(INTERNAL, trim(field_name), field, __RC__)

            call ESMF_AttributeSet(field, NAME='ScavengingFractionPerKm', VALUE=f_scav, __RC__)
        end do

        ! cloud-borne aerosol tracers
        field_name = 'NUM_C_' // trim(mode_name)
        call ESMF_StateGet(INTERNAL, trim(field_name), field, __RC__)
        call ESMF_AttributeSet(field, NAME='ScavengingFractionPerKm', VALUE=f_scav_none, __RC__)

        do s = 1, n_species
            species_name = self%scheme%mode(m)%species(s)%name
            field_name  = trim(species_name) // '_C_' // trim(mode_name)
            call ESMF_StateGet(INTERNAL, trim(field_name), field, __RC__)
            call ESMF_AttributeSet(field, NAME='ScavengingFractionPerKm', VALUE=f_scav_none, __RC__)
        end do
    end do

#ifdef PRINT_STATES
    if (MAPL_AM_I_ROOT()) then
        print *, trim(Iam)//': AERO State during Initialize():'
        call ESMF_StatePrint(aero, nestedFlag=.true., __RC__)
    end if
#endif

!   Create narrow-band Mie tables
!   -----------------------------
    call ESMF_ConfigGetAttribute ( self%CF, mie_optics_file, Label='narrowband_optics_ait:', default='', __RC__ )
    self%mie_ait = MAML_OpticsTableCreate(trim(mie_optics_file), __RC__)
    call MAML_OpticsTableRead(self%mie_ait, __RC__)

    call ESMF_ConfigGetAttribute ( self%CF, mie_optics_file, Label='narrowband_optics_acc:', default='', __RC__ )
    self%mie_acc = MAML_OpticsTableCreate(trim(mie_optics_file), __RC__)
    call MAML_OpticsTableRead(self%mie_acc, __RC__)

    call ESMF_ConfigGetAttribute ( self%CF, mie_optics_file, Label='narrowband_optics_pcm:', default='', __RC__ )
    self%mie_pcm = MAML_OpticsTableCreate(trim(mie_optics_file), __RC__)
    call MAML_OpticsTableRead(self%mie_pcm, __RC__)

    call ESMF_ConfigGetAttribute ( self%CF, mie_optics_file, Label='narrowband_optics_fss:', default='', __RC__ )
    self%mie_fss = MAML_OpticsTableCreate(trim(mie_optics_file), __RC__)
    call MAML_OpticsTableRead(self%mie_fss, __RC__)

    call ESMF_ConfigGetAttribute ( self%CF, mie_optics_file, Label='narrowband_optics_css:', default='', __RC__ )
    self%mie_css = MAML_OpticsTableCreate(trim(mie_optics_file), __RC__)
    call MAML_OpticsTableRead(self%mie_css, __RC__)

    call ESMF_ConfigGetAttribute ( self%CF, mie_optics_file, Label='narrowband_optics_fdu:', default='', __RC__ )
    self%mie_fdu = MAML_OpticsTableCreate(trim(mie_optics_file), __RC__)
    call MAML_OpticsTableRead(self%mie_fdu, __RC__)

    call ESMF_ConfigGetAttribute ( self%CF, mie_optics_file, Label='narrowband_optics_cdu:', default='', __RC__ )
    self%mie_cdu = MAML_OpticsTableCreate(trim(mie_optics_file), __RC__)
    call MAML_OpticsTableRead(self%mie_cdu, __RC__)


!   All done
!   --------
    call MAPL_TimerOff(mgState, 'INITIALIZE', __RC__)
    call MAPL_TimerOff(mgState, 'TOTAL',      __RC__)

    _RETURN(ESMF_SUCCESS)

    contains
   
    subroutine PrintProcessFlag(name, state)
        implicit none        
        character(len=*), intent(in) :: name
        logical, intent(in)          :: state
       
        character(len=*), parameter  :: fmt_flag  = "(4X, A12, X, '-', X, A3)"

        if (state) then
            write (*, fmt_flag) trim(name), 'ON'
        else
            write (*, fmt_flag) trim(name), 'OFF'
        end if
    end subroutine PrintProcessFlag


    subroutine microphysics_initialize(imozart, verbose, rc)
        use MAPL_ConstantsMod, only: r8 => MAPL_R8

        use constituents,      only: pcnst, cnst_name, cnst_longname
        use chem_mods,         only: gas_pcnst, adv_mass

        use modal_aero_data,   only: ntot_amode, ntot_aspectype,       nspec_amode, numptr_amode, lmassptr_amode, numptrcw_amode, lmassptrcw_amode
        use modal_aero_initialize_data, only: modal_aero_register, modal_aero_initialize

        use MAM_ComponentsDataMod

        implicit none

        integer, intent(in)  :: imozart
        logical, intent(in)  :: verbose
        integer, intent(out) :: rc

        ! local
        real(r8) :: sigma(ntot_amode)
        real(r8) :: dgn(ntot_amode)
        real(r8) :: dgn_low(ntot_amode)
        real(r8) :: dgn_hi(ntot_amode)
        real(r8) :: rh_crystal(ntot_amode)
        real(r8) :: rh_deliques(ntot_amode)
        real(r8) :: spec_dens(ntot_aspectype)
        real(r8) :: spec_hygro(ntot_aspectype)

#ifdef DEBUG
        integer  :: m, l
#endif

        __Iam__('MAM::microphysics_initialize')


        ! this is probably done somewhere in CAM/chem or CAM/phys
        cnst_name(:)     = '__NONE__'
        cnst_longname(:) = '__NONE__'
        adv_mass(:)      = 0.0_r8

#if ( defined MODAL_AERO_7MODE )
        !
        ! based on files in pp_trop_mam7 from CESM-1.2.1
        !
        cnst_name(1:pcnst) = (/'H2O2            ','H2SO4           ','SO2             ','DMS             ','NH3             ', &
                               'SOAG            ','so4_a1          ','nh4_a1          ','pom_a1          ','soa_a1          ', &
                               'bc_a1           ','ncl_a1          ','num_a1          ','so4_a2          ','nh4_a2          ', &
                               'soa_a2          ','ncl_a2          ','num_a2          ','pom_a3          ','bc_a3           ', &
                               'num_a3          ','ncl_a4          ','so4_a4          ','nh4_a4          ','num_a4          ', &
                               'dst_a5          ','so4_a5          ','nh4_a5          ','num_a5          ','ncl_a6          ', &
                               'so4_a6          ','nh4_a6          ','num_a6          ','dst_a7          ','so4_a7          ', &
                               'nh4_a7          ','num_a7          ' /) 
    
        adv_mass(1:gas_pcnst) = (/ 34.013600_r8,    98.078400_r8,    64.064800_r8,    62.132400_r8,    17.028940_r8, &
                                   12.011000_r8,    96.063600_r8,    18.036340_r8,    12.011000_r8,    12.011000_r8, &
                                   12.011000_r8,    58.442468_r8,     1.007400_r8,    96.063600_r8,    18.036340_r8, &
                                   12.011000_r8,    58.442468_r8,     1.007400_r8,    12.011000_r8,    12.011000_r8, &
                                    1.007400_r8,    58.442468_r8,    96.063600_r8,    18.036340_r8,     1.007400_r8, &
                                  135.064039_r8,    96.063600_r8,    18.036340_r8,     1.007400_r8,    58.442468_r8, &
                                   96.063600_r8,    18.036340_r8,     1.007400_r8,   135.064039_r8,    96.063600_r8, &
                                   18.036340_r8,     1.007400_r8 /)


       sigma(:ntot_amode)          = (/ MAM7_ACCUMULATION_MODE_SIGMA, &
                                        MAM7_AITKEN_MODE_SIGMA, &
                                        MAM7_PRIMARY_CARBON_MODE_SIGMA, &
                                        MAM7_FINE_SEASALT_MODE_SIGMA, &
                                        MAM7_FINE_DUST_MODE_SIGMA, &
                                        MAM7_COARSE_SEASALT_MODE_SIGMA, &
                                        MAM7_COARSE_DUST_MODE_SIGMA /)

       dgn(:ntot_amode)            = (/ MAM7_ACCUMULATION_MODE_SIZE,&
                                        MAM7_AITKEN_MODE_SIZE, &
                                        MAM7_PRIMARY_CARBON_MODE_SIZE, &
                                        MAM7_FINE_SEASALT_MODE_SIZE, &
                                        MAM7_FINE_DUST_MODE_SIZE, &
                                        MAM7_COARSE_SEASALT_MODE_SIZE, &
                                        MAM7_COARSE_DUST_MODE_SIZE /)


       dgn_low(:ntot_amode)        = (/ MAM7_ACCUMULATION_MODE_SIZE_MIN, &
                                        MAM7_AITKEN_MODE_SIZE_MIN, &
                                        MAM7_PRIMARY_CARBON_MODE_SIZE_MIN, &
                                        MAM7_FINE_SEASALT_MODE_SIZE_MIN, &
                                        MAM7_FINE_DUST_MODE_SIZE_MIN, &
                                        MAM7_COARSE_SEASALT_MODE_SIZE_MIN, &
                                        MAM7_COARSE_DUST_MODE_SIZE_MIN /)

       dgn_hi(:ntot_amode)         = (/ MAM7_ACCUMULATION_MODE_SIZE_MAX, &
                                        MAM7_AITKEN_MODE_SIZE_MAX, &
                                        MAM7_PRIMARY_CARBON_MODE_SIZE_MAX, &
                                        MAM7_FINE_SEASALT_MODE_SIZE_MAX, &
                                        MAM7_FINE_DUST_MODE_SIZE_MAX, &
                                        MAM7_COARSE_SEASALT_MODE_SIZE_MAX, &
                                        MAM7_COARSE_DUST_MODE_SIZE_MAX /)


       rh_crystal(:ntot_amode)     = (/ MAM7_ACCUMULATION_MODE_RH_CRYSTALLIZATION, &
                                        MAM7_AITKEN_MODE_RH_CRYSTALLIZATION, &
                                        MAM7_PRIMARY_CARBON_MODE_RH_CRYSTALLIZATION, &
                                        MAM7_FINE_SEASALT_MODE_RH_CRYSTALLIZATION, &
                                        MAM7_FINE_DUST_MODE_RH_CRYSTALLIZATION, &
                                        MAM7_COARSE_SEASALT_MODE_RH_CRYSTALLIZATION, &
                                        MAM7_COARSE_DUST_MODE_RH_CRYSTALLIZATION /)

       rh_deliques(:ntot_amode)    = (/ MAM7_ACCUMULATION_MODE_RH_DELIQUESCENCE, &
                                        MAM7_AITKEN_MODE_RH_DELIQUESCENCE, &
                                        MAM7_PRIMARY_CARBON_MODE_RH_DELIQUESCENCE, &
                                        MAM7_FINE_SEASALT_MODE_RH_DELIQUESCENCE, &
                                        MAM7_FINE_DUST_MODE_RH_DELIQUESCENCE, &
                                        MAM7_COARSE_SEASALT_MODE_RH_DELIQUESCENCE, &
                                        MAM7_COARSE_DUST_MODE_RH_DELIQUESCENCE /)

#endif

       ! following the indexes in specname_amode(:ntot_aspectype)
       ! specname_amode(ntot_aspectype) = (/ 'sulfate   ', 'ammonium  ', 'nitrate   ', 'p-organic ', 's-organic ', 'black-c   ', 'seasalt   ', 'dust      '/)
       spec_dens(:ntot_aspectype)  = (/ MAM_SULFATE_COMPONENT_DENSITY, &
                                        MAM_AMMONIUM_COMPONENT_DENSITY, &
                                        MAM_SULFATE_COMPONENT_DENSITY, &  ! assigned value for nitrate: verify what is used in CAM
                                        MAM_POM_COMPONENT_DENSITY, &
                                        MAM_SOA_COMPONENT_DENSITY, &
                                        MAM_BLACK_CARBON_COMPONENT_DENSITY, &
                                        MAM_SEASALT_COMPONENT_DENSITY, &
                                        MAM_DUST_COMPONENT_DENSITY /)

       spec_hygro(:ntot_aspectype) = (/ MAM_SULFATE_COMPONENT_HYGROSCOPICITY, &
                                        MAM_AMMONIUM_COMPONENT_HYGROSCOPICITY, &
                                        MAM_SULFATE_COMPONENT_HYGROSCOPICITY, &  ! assigned value for nitrate: verify what is used in CAM
                                        MAM_POM_COMPONENT_HYGROSCOPICITY, &
                                        MAM_SOA_COMPONENT_HYGROSCOPICITY, &
                                        MAM_BLACK_CARBON_COMPONENT_HYGROSCOPICITY, &
                                        MAM_SEASALT_COMPONENT_HYGROSCOPICITY, &
                                        MAM_DUST_COMPONENT_HYGROSCOPICITY /)


        call modal_aero_register(verbose)

        call modal_aero_initialize(imozart, sigma, dgn, dgn_low, dgn_hi, rh_crystal, rh_deliques, spec_dens, spec_hygro, verbose)

#ifdef DEBUG
        if (MAPL_AM_I_ROOT()) then
       
            print *, 'Interstitial aerosols:'

            do m = 1, ntot_amode
                print *, 'mode:  numptr_amode(m) = ', m, numptr_amode(m)
            end do
            print *

            do m = 1, ntot_amode
                do l = 1, nspec_amode(m)
                    print *, 'mode, species: lmassptr_amode(l,m) = ', l, m,  lmassptr_amode(l,m)
                end do
                print *
            end do

            print *

            print *, 'Cloud-borne  aerosols:'

            do m = 1, ntot_amode
                print *, 'mode:  numptrcw_amode(m) = ', m, numptrcw_amode(m)
            end do
            print *

            do m = 1, ntot_amode
                do l = 1, nspec_amode(m)
                    print *, 'mode, species: lmassptrcw_amode(l,m) = ', l, m,  lmassptrcw_amode(l,m)
                end do
                print *
            end do

        end if
#endif


        _RETURN(ESMF_SUCCESS)

    end subroutine microphysics_initialize

   end subroutine Initialize_


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Run_ --- Runs MAMchem
!
! !INTERFACE:
!

   subroutine Run_ ( GC, IMPORT, EXPORT, CLOCK, rc )

! !USES:
    use MAPL_ConstantsMod,   only: r8 => MAPL_R8
    use cam_logfile,         only: iulog    
    use chem_mods,           only: gas_pcnst, adv_mass

    use modal_aero_data,     only: ntot_amode
    use modal_aero_amicphys, only: pcols, pver, pcnstxx, &
                                   nqtendaa, iqtend_cond, iqtend_rnam, iqtend_nnuc, iqtend_coag, &
                                   nqqcwtendaa, iqqcwtend_rnam, &
                                   modal_aero_amicphys_intr

    use MAM7_DataMod,        only: MAM7_AITKEN_MODE_SIZE, &
                                   MAM7_ACCUMULATION_MODE_SIZE, &
                                   MAM7_PRIMARY_CARBON_MODE_SIZE, &
                                   MAM7_FINE_SEASALT_MODE_SIZE, &
                                   MAM7_FINE_DUST_MODE_SIZE, &
                                   MAM7_COARSE_SEASALT_MODE_SIZE, &
                                   MAM7_COARSE_DUST_MODE_SIZE

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
   
    type(MAM_state), pointer      :: self               ! Legacy state
    type(ESMF_Grid)               :: GRID               ! Grid
    type(ESMF_Config)             :: CF                 ! Universal Config 

    type(MAPL_MetaComp), pointer  :: mgState            ! MAPL generic state
    type(ESMF_Alarm)              :: run_alarm
    logical                       :: run_alarm_ringing

    integer                       :: im_World, jm_World ! Global 2D Dimensions
    integer                       :: im, jm, lm         ! 3D Dimensions
    real(ESMF_KIND_R4), pointer   :: lons(:,:)          ! Longitudes
    real(ESMF_KIND_R4), pointer   :: lats(:,:)          ! Latitudes

    integer                       :: nymd, nhms         ! date, time
    real                          :: cdt                ! time step in secs

    character(len=ESMF_MAXSTR)    :: comp_name

    
    ! inputs to the aerosol microphysics core
    integer                         :: i, j
    integer                         :: iq

    real, pointer, dimension(:,:)   :: zpbl
    real, pointer, dimension(:,:,:) :: fcld
    real, pointer, dimension(:,:,:) :: Q
    real, pointer, dimension(:,:,:) :: T
    real, pointer, dimension(:,:,:) :: RH 
    real, pointer, dimension(:,:,:) :: delp
    real, pointer, dimension(:,:,:) :: ple 
    real, pointer, dimension(:,:,:) :: zle

!   real, pointer, dimension(:,:,:) :: h2o2
!   real, pointer, dimension(:,:,:) :: dms
!   real, pointer, dimension(:,:,:) :: msa
    real, pointer, dimension(:,:,:) :: so2
    real, pointer, dimension(:,:,:) :: h2so4
    real, pointer, dimension(:,:,:) :: nh3
    real, pointer, dimension(:,:,:) :: soa_g

    real, pointer, dimension(:,:,:) :: ddt_dms_gas
    real, pointer, dimension(:,:,:) :: ddt_msa_gas
    real, pointer, dimension(:,:,:) :: ddt_so2_gas
    real, pointer, dimension(:,:,:) :: ddt_h2so4_gas
    real, pointer, dimension(:,:,:) :: ddt_nh3_gas
    real, pointer, dimension(:,:,:) :: ddt_soa_g_gas

    real, pointer, dimension(:,:,:) :: ddt_dms_aq
    real, pointer, dimension(:,:,:) :: ddt_msa_aq
    real, pointer, dimension(:,:,:) :: ddt_so2_aq
    real, pointer, dimension(:,:,:) :: ddt_h2so4_aq
    real, pointer, dimension(:,:,:) :: ddt_nh3_aq
    real, pointer, dimension(:,:,:) :: ddt_soa_g_aq

    real, pointer, dimension(:,:,:) :: dms_g_
    real, pointer, dimension(:,:,:) :: msa_g_
    real, pointer, dimension(:,:,:) :: so2_g_
    real, pointer, dimension(:,:,:) :: h2so4_g_
    real, pointer, dimension(:,:,:) :: nh3_g_
    real, pointer, dimension(:,:,:) :: soa_g_g_

    real, pointer, dimension(:,:,:) :: dms_a_
    real, pointer, dimension(:,:,:) :: msa_a_
    real, pointer, dimension(:,:,:) :: so2_a_
    real, pointer, dimension(:,:,:) :: h2so4_a_
    real, pointer, dimension(:,:,:) :: nh3_a_
    real, pointer, dimension(:,:,:) :: soa_g_a_

    real, pointer, dimension(:,:,:) :: ait_a_num
    real, pointer, dimension(:,:,:) :: ait_a_so4
    real, pointer, dimension(:,:,:) :: ait_a_nh4
    real, pointer, dimension(:,:,:) :: ait_a_soa
    real, pointer, dimension(:,:,:) :: ait_a_ncl
    
    real, pointer, dimension(:,:,:) :: ait_c_num
    real, pointer, dimension(:,:,:) :: ait_c_so4
    real, pointer, dimension(:,:,:) :: ait_c_nh4
    real, pointer, dimension(:,:,:) :: ait_c_soa
    real, pointer, dimension(:,:,:) :: ait_c_ncl

    
    real, pointer, dimension(:,:,:) :: acc_a_num
    real, pointer, dimension(:,:,:) :: acc_a_so4
    real, pointer, dimension(:,:,:) :: acc_a_nh4
    real, pointer, dimension(:,:,:) :: acc_a_soa
    real, pointer, dimension(:,:,:) :: acc_a_pom
    real, pointer, dimension(:,:,:) :: acc_a_bc
    real, pointer, dimension(:,:,:) :: acc_a_ncl

    real, pointer, dimension(:,:,:) :: acc_c_num
    real, pointer, dimension(:,:,:) :: acc_c_so4
    real, pointer, dimension(:,:,:) :: acc_c_nh4
    real, pointer, dimension(:,:,:) :: acc_c_soa
    real, pointer, dimension(:,:,:) :: acc_c_pom
    real, pointer, dimension(:,:,:) :: acc_c_bc
    real, pointer, dimension(:,:,:) :: acc_c_ncl


    real, pointer, dimension(:,:,:) :: pcm_a_num
    real, pointer, dimension(:,:,:) :: pcm_a_pom
    real, pointer, dimension(:,:,:) :: pcm_a_bc

    real, pointer, dimension(:,:,:) :: pcm_c_num
    real, pointer, dimension(:,:,:) :: pcm_c_pom
    real, pointer, dimension(:,:,:) :: pcm_c_bc

 
    real, pointer, dimension(:,:,:) :: fdu_a_num
    real, pointer, dimension(:,:,:) :: fdu_a_dst
    real, pointer, dimension(:,:,:) :: fdu_a_so4
    real, pointer, dimension(:,:,:) :: fdu_a_nh4

    real, pointer, dimension(:,:,:) :: fdu_c_num
    real, pointer, dimension(:,:,:) :: fdu_c_dst
    real, pointer, dimension(:,:,:) :: fdu_c_so4
    real, pointer, dimension(:,:,:) :: fdu_c_nh4


    real, pointer, dimension(:,:,:) :: cdu_a_num
    real, pointer, dimension(:,:,:) :: cdu_a_dst
    real, pointer, dimension(:,:,:) :: cdu_a_so4
    real, pointer, dimension(:,:,:) :: cdu_a_nh4

    real, pointer, dimension(:,:,:) :: cdu_c_num
    real, pointer, dimension(:,:,:) :: cdu_c_dst
    real, pointer, dimension(:,:,:) :: cdu_c_so4
    real, pointer, dimension(:,:,:) :: cdu_c_nh4


    real, pointer, dimension(:,:,:) :: fss_a_num
    real, pointer, dimension(:,:,:) :: fss_a_ncl
    real, pointer, dimension(:,:,:) :: fss_a_so4
    real, pointer, dimension(:,:,:) :: fss_a_nh4

    real, pointer, dimension(:,:,:) :: fss_c_num
    real, pointer, dimension(:,:,:) :: fss_c_ncl
    real, pointer, dimension(:,:,:) :: fss_c_so4
    real, pointer, dimension(:,:,:) :: fss_c_nh4


    real, pointer, dimension(:,:,:) :: css_a_num
    real, pointer, dimension(:,:,:) :: css_a_ncl
    real, pointer, dimension(:,:,:) :: css_a_so4
    real, pointer, dimension(:,:,:) :: css_a_nh4

    real, pointer, dimension(:,:,:) :: css_c_num
    real, pointer, dimension(:,:,:) :: css_c_ncl
    real, pointer, dimension(:,:,:) :: css_c_so4
    real, pointer, dimension(:,:,:) :: css_c_nh4


    real, pointer, dimension(:,:,:) :: ait_a_wtr
    real, pointer, dimension(:,:,:) :: ait_dgn_dry
    real, pointer, dimension(:,:,:) :: ait_dgn_wet

    real, pointer, dimension(:,:,:) :: acc_a_wtr
    real, pointer, dimension(:,:,:) :: acc_dgn_dry
    real, pointer, dimension(:,:,:) :: acc_dgn_wet

    real, pointer, dimension(:,:,:) :: pcm_a_wtr
    real, pointer, dimension(:,:,:) :: pcm_dgn_dry
    real, pointer, dimension(:,:,:) :: pcm_dgn_wet

    real, pointer, dimension(:,:,:) :: fdu_a_wtr
    real, pointer, dimension(:,:,:) :: fdu_dgn_dry
    real, pointer, dimension(:,:,:) :: fdu_dgn_wet

    real, pointer, dimension(:,:,:) :: cdu_a_wtr
    real, pointer, dimension(:,:,:) :: cdu_dgn_dry
    real, pointer, dimension(:,:,:) :: cdu_dgn_wet

    real, pointer, dimension(:,:,:) :: fss_a_wtr 
    real, pointer, dimension(:,:,:) :: fss_dgn_dry
    real, pointer, dimension(:,:,:) :: fss_dgn_wet

    real, pointer, dimension(:,:,:) :: css_a_wtr
    real, pointer, dimension(:,:,:) :: css_dgn_dry
    real, pointer, dimension(:,:,:) :: css_dgn_wet

    ! pre-aqueous chemistry SO4 and NH4
    real, allocatable, dimension(:,:,:) :: ait_a_so4_
    real, allocatable, dimension(:,:,:) :: ait_a_nh4_
    
    real, allocatable, dimension(:,:,:) :: acc_a_so4_
    real, allocatable, dimension(:,:,:) :: acc_a_nh4_

    real, allocatable, dimension(:,:,:) :: fdu_a_so4_
    real, allocatable, dimension(:,:,:) :: fdu_a_nh4_

    real, allocatable, dimension(:,:,:) :: cdu_a_so4_
    real, allocatable, dimension(:,:,:) :: cdu_a_nh4_

    real, allocatable, dimension(:,:,:) :: fss_a_so4_
    real, allocatable, dimension(:,:,:) :: fss_a_nh4_

    real, allocatable, dimension(:,:,:) :: css_a_so4_
    real, allocatable, dimension(:,:,:) :: css_a_nh4_

    real, allocatable, dimension(:,:,:) :: css_c_so4_
    real, allocatable, dimension(:,:,:) :: css_c_nh4_

    real, allocatable, dimension(:,:,:) :: q_coltend_cond_
    real, allocatable, dimension(:,:,:) :: q_coltend_rename_
    real, allocatable, dimension(:,:,:) :: q_coltend_coag_
    real, allocatable, dimension(:,:,:) :: q_coltend_nucl_
    real, allocatable, dimension(:,:,:) :: qqcw_coltend_rename_
   
    real, pointer,     dimension(:,:)   :: q_coltend

    ! wrap the aerosol microphisics core 
    integer, parameter  :: ncol = 1                     ! number of atmospheric columns

    integer  :: amc_do_gasaerexch
    integer  :: amc_do_rename
    integer  :: amc_do_newnuc
    integer  :: amc_do_coag

    integer  :: amc_lchnk                               ! chunk identifier
    integer  :: amc_nstep                               ! model time-step number
    integer  :: amc_loffset                             ! offset applied to modal aero "ptrs"
    integer  :: amc_latndx(1), amc_lonndx(1)

    real(r8) :: amc_deltat                              ! time step (s)

    real(r8) :: amc_dqdt(ncol,pver,pcnstxx)             ! tendency of q
    logical  :: amc_dotend(pcnstxx)                     ! flag

    real(r8) :: amc_q(ncol,pver,pcnstxx)                ! current tracer mixing ratios (TMRs)
                                                        ! these values are updated (so out /= in)
                                                        ! *** MUST BE  #/kmol-air for number
                                                        ! *** MUST BE mol/mol-air for mass
                                                        ! *** NOTE ncol dimension

    real(r8) :: amc_qqcw(ncol,pver,pcnstxx)             ! like q but for cloud-borner tracers
                                                        ! these values are updated

    real(r8) :: amc_q_pregaschem(ncol,pver,pcnstxx)     ! q TMRs    before gas-phase chemistry
    real(r8) :: amc_q_precldchem(ncol,pver,pcnstxx)     ! q TMRs    before cloud chemistry
    real(r8) :: amc_qqcw_precldchem(ncol,pver,pcnstxx)  ! qqcw TMRs before cloud chemistry

    real(r8) :: amc_t(pcols,pver)                       ! temperature at model levels (K)
    real(r8) :: amc_pmid(pcols,pver)                    ! pressure at model level centers (Pa)
    real(r8) :: amc_pdel(pcols,pver)                    ! pressure thickness of levels (Pa)
    real(r8) :: amc_zm(pcols,pver)                      ! altitude (above ground) at level centers (m)
    real(r8) :: amc_pblh(pcols)                         ! planetary boundary layer depth (m)
    real(r8) :: amc_qv(pcols,pver)                      ! specific humidity (kg/kg)
    real(r8) :: amc_rh(pcols,pver)                      ! relative humidity (0, 1)
    real(r8) :: amc_cld(ncol,pver)                      ! cloud fraction (-) *** NOTE ncol dimension
    real(r8) :: amc_dgn_a_dry(pcols,pver,ntot_amode)
    real(r8) :: amc_dgn_a_wet(pcols,pver,ntot_amode)    ! dry & wet geo. mean dia. (m) of number distrib.
    real(r8) :: amc_wetdens_host(pcols,pver,ntot_amode) ! interstitial aerosol wet density (kg/m3)
    real(r8) :: amc_q_coltendaa(pcols,pcnstxx,nqtendaa) ! column-integrated tendencies for condensation, renaming, coagulation, and nucleation 
    real(r8) :: amc_qqcw_coltendaa(pcols,pcnstxx,nqqcwtendaa) ! --dito-- but for cloud-borne aerosols
    real(r8) :: amc_qaerwat(pcols,pver,ntot_amode)      ! optional, aerosol water mixing ratio (kg/kg)

    real(r8) :: tmp_min, tmp_max

    integer  :: i_H2O2, i_H2SO4, i_SO2, i_DMS, i_NH3, i_SOAG

    integer  :: i_so4_a1, i_nh4_a1, i_pom_a1, i_soa_a1, i_bc_a1, i_ncl_a1, i_num_a1
    integer  :: i_so4_a2, i_nh4_a2, i_soa_a2, i_ncl_a2, i_num_a2
    integer  :: i_pom_a3, i_bc_a3 , i_num_a3
    integer  :: i_ncl_a4, i_so4_a4, i_nh4_a4, i_num_a4
    integer  :: i_dst_a5, i_so4_a5, i_nh4_a5, i_num_a5
    integer  :: i_ncl_a6, i_so4_a6, i_nh4_a6, i_num_a6
    integer  :: i_dst_a7, i_so4_a7, i_nh4_a7, i_num_a7

    real, parameter :: mw_air = 28.97                   ! molar mass of dry air, g mol-1

!   Declare pointers to IMPORT/EXPORT/INTERNAL states 
!   -------------------------------------------------
    type(ESMF_State)              :: INTERNAL
  
!   Get my name and set-up traceback handle
!   ---------------------------------------
    call ESMF_GridCompGet( GC, name=comp_name, __RC__ )
    Iam = trim(comp_name) // '::' // trim(Iam)

!   Get my internal MAPL_Generic state
!   -----------------------------------
    call MAPL_GetObjectFromGC(GC, mgState, __RC__)

    call MAPL_TimerOn(mgState, 'TOTAL', __RC__)
    call MAPL_TimerOn(mgState, 'RUN',   __RC__)

!   Get pointers to IMPORT/EXPORT/INTERNAL states 
!   ---------------------------------------------
    call MAPL_Get (mgState, INTERNAL_ESMF_STATE=INTERNAL, __RC__)



!   Get parameters from generic state
!   ----------------------------------
    call MAPL_Get(mgState, LONS=lons, &
                           LATS=lats, &
                           RunAlarm=run_alarm, __RC__)


!   If it is time, update MAM state
!   -------------------------------
    run_alarm_ringing = ESMF_AlarmIsRinging(run_alarm, __RC__)

    if (run_alarm_ringing) then
        call ESMF_AlarmRingerOff(run_alarm, __RC__)
    else
        _RETURN(ESMF_SUCCESS)
    endif


!   Extract relevant runtime information
!   ------------------------------------
    call extract_(GC, CLOCK, self, GRID, CF, &
                  im_World, jm_World,        &
                  im, jm, lm,                &
                  nymd, nhms, cdt, __RC__)


!   Force non-negative mixing ratios
!   --------------------------------
    do iq = 1, self%qa%n3d
        where(self%qa%r3(iq)%q < 0) self%qa%r3(iq)%q = tiny(0.0)
    end do

    do iq = 1, self%qc%n3d
        where(self%qc%r3(iq)%q < 0) self%qc%r3(iq)%q = tiny(0.0)
    end do

    do iq = 1, self%qg%n3d
        where(self%qg%r3(iq)%q < 0) self%qg%r3(iq)%q = tiny(0.0)
    end do


!   Partition aqueous phase production of SO4 and NH4
!   -------------------------------------------------
    call MAPL_GetPointer(internal, ait_a_so4, 'SU_A_AIT' , __RC__)
    call MAPL_GetPointer(internal, ait_a_nh4, 'AMM_A_AIT', __RC__)
    call MAPL_GetPointer(internal, acc_a_so4, 'SU_A_ACC' , __RC__)
    call MAPL_GetPointer(internal, acc_a_nh4, 'AMM_A_ACC', __RC__)
    call MAPL_GetPointer(internal, fdu_a_so4, 'SU_A_FDU' , __RC__)
    call MAPL_GetPointer(internal, fdu_a_nh4, 'AMM_A_FDU', __RC__)
    call MAPL_GetPointer(internal, cdu_a_so4, 'SU_A_CDU' , __RC__)
    call MAPL_GetPointer(internal, cdu_a_nh4, 'AMM_A_CDU', __RC__)
    call MAPL_GetPointer(internal, fss_a_so4, 'SU_A_FSS' , __RC__)
    call MAPL_GetPointer(internal, fss_a_nh4, 'AMM_A_FSS', __RC__)
    call MAPL_GetPointer(internal, css_a_so4, 'SU_A_CSS' , __RC__)
    call MAPL_GetPointer(internal, css_a_nh4, 'AMM_A_CSS', __RC__)

    ! save pre-aqueous phase SO4 and NH4
    allocate(ait_a_so4_(im,jm,lm), __STAT__)
    allocate(ait_a_nh4_(im,jm,lm), __STAT__)    
    allocate(acc_a_so4_(im,jm,lm), __STAT__)
    allocate(acc_a_nh4_(im,jm,lm), __STAT__)
    allocate(fdu_a_so4_(im,jm,lm), __STAT__)
    allocate(fdu_a_nh4_(im,jm,lm), __STAT__)
    allocate(cdu_a_so4_(im,jm,lm), __STAT__)
    allocate(cdu_a_nh4_(im,jm,lm), __STAT__)
    allocate(fss_a_so4_(im,jm,lm), __STAT__)
    allocate(fss_a_nh4_(im,jm,lm), __STAT__)
    allocate(css_a_so4_(im,jm,lm), __STAT__)
    allocate(css_a_nh4_(im,jm,lm), __STAT__)

    ait_a_so4_ = ait_a_so4
    ait_a_nh4_ = ait_a_nh4
    acc_a_so4_ = acc_a_so4
    acc_a_nh4_ = acc_a_nh4
    fdu_a_so4_ = fdu_a_so4
    fdu_a_nh4_ = fdu_a_nh4
    cdu_a_so4_ = cdu_a_so4
    cdu_a_nh4_ = cdu_a_nh4
    fss_a_so4_ = fss_a_so4
    fss_a_nh4_ = fss_a_nh4
    css_a_so4_ = css_a_so4
    css_a_nh4_ = css_a_nh4

    call MAPL_TimerOn(mgState, '-AQUEOUS_CHEM', __RC__)
    call AqueousChemistry(self%scheme, import, export, self%qa, cdt, rc)
    call MAPL_TimerOff(mgState, '-AQUEOUS_CHEM', __RC__)


    ! colmn-integrated diagnostics: tendencies due to condensation, 
    ! rename, coagulation and nucleation
    allocate(q_coltend_cond_(im,jm,pcnstxx),      __STAT__)
    allocate(q_coltend_rename_(im,jm,pcnstxx),    __STAT__)
    allocate(q_coltend_coag_(im,jm,pcnstxx),      __STAT__)
    allocate(q_coltend_nucl_(im,jm,pcnstxx),      __STAT__)
    allocate(qqcw_coltend_rename_(im,jm,pcnstxx), __STAT__)

    q_coltend_cond_   = 0.0
    q_coltend_rename_ = 0.0
    q_coltend_coag_   = 0.0
    q_coltend_nucl_   = 0.0
    qqcw_coltend_rename_ = 0.0


    call MAPL_TimerOn(mgState, '-MICROPHYSICS', __RC__)

!   Calculate size of the particles --- CAM interface
!   ------------------------------------------------------
!   call CAM_CalculateSize(self%scheme, self%qa, self%qc, self%Dg_dry, self%Dg_wet, delp, self%dt, rc)



!   Aerosol microphysics core
!   -------------------------

    i_H2O2   = constituent_index_('H2O2'  , __RC__)
    i_H2SO4  = constituent_index_('H2SO4' , __RC__)
    i_SO2    = constituent_index_('SO2'   , __RC__)
    i_DMS    = constituent_index_('DMS'   , __RC__)
    i_NH3    = constituent_index_('NH3'   , __RC__)
    i_SOAG   = constituent_index_('SOAG'  , __RC__)
    i_so4_a1 = constituent_index_('so4_a1', __RC__)
    i_nh4_a1 = constituent_index_('nh4_a1', __RC__)
    i_pom_a1 = constituent_index_('pom_a1', __RC__)
    i_soa_a1 = constituent_index_('soa_a1', __RC__)
    i_bc_a1  = constituent_index_('bc_a1' , __RC__)
    i_ncl_a1 = constituent_index_('ncl_a1', __RC__)
    i_num_a1 = constituent_index_('num_a1', __RC__)
    i_so4_a2 = constituent_index_('so4_a2', __RC__)
    i_nh4_a2 = constituent_index_('nh4_a2', __RC__)
    i_soa_a2 = constituent_index_('soa_a2', __RC__)
    i_ncl_a2 = constituent_index_('ncl_a2', __RC__)
    i_num_a2 = constituent_index_('num_a2', __RC__)
    i_pom_a3 = constituent_index_('pom_a3', __RC__)
    i_bc_a3  = constituent_index_('bc_a3' , __RC__)
    i_num_a3 = constituent_index_('num_a3', __RC__)
    i_ncl_a4 = constituent_index_('ncl_a4', __RC__)
    i_so4_a4 = constituent_index_('so4_a4', __RC__)
    i_nh4_a4 = constituent_index_('nh4_a4', __RC__)
    i_num_a4 = constituent_index_('num_a4', __RC__)
    i_dst_a5 = constituent_index_('dst_a5', __RC__)
    i_so4_a5 = constituent_index_('so4_a5', __RC__)
    i_nh4_a5 = constituent_index_('nh4_a5', __RC__)
    i_num_a5 = constituent_index_('num_a5', __RC__)
    i_ncl_a6 = constituent_index_('ncl_a6', __RC__)
    i_so4_a6 = constituent_index_('so4_a6', __RC__)
    i_nh4_a6 = constituent_index_('nh4_a6', __RC__)
    i_num_a6 = constituent_index_('num_a6', __RC__)
    i_dst_a7 = constituent_index_('dst_a7', __RC__)
    i_so4_a7 = constituent_index_('so4_a7', __RC__)
    i_nh4_a7 = constituent_index_('nh4_a7', __RC__)
    i_num_a7 = constituent_index_('num_a7', __RC__)


    ! set initial values
    amc_do_gasaerexch = 0
    amc_do_rename     = 0
    amc_do_newnuc     = 0
    amc_do_coag       = 0

    amc_nstep   = 99                            ! model time-step number
    amc_lchnk   = 0                             ! chunk identifier
    amc_latndx  = 1 
    amc_lonndx  = 1                             ! lat and lon indices
    amc_loffset = 0


    if (self%condensation) amc_do_gasaerexch = 1
    if (self%rename)       amc_do_rename     = 1
    if (self%nucleation)   amc_do_newnuc     = 1
    if (self%coagulation)  amc_do_coag       = 1


    amc_deltat = self%dt                       ! time step (s)

    call MAPL_GetPointer(import, ple,   'PLE',   __RC__)
    call MAPL_GetPointer(import, delp,  'DELP',  __RC__)
    call MAPL_GetPointer(import, fcld,  'FCLD',  __RC__)
    call MAPL_GetPointer(import, Q,     'Q',     __RC__)
    call MAPL_GetPointer(import, RH,    'RH2',   __RC__)
    call MAPL_GetPointer(import, T,     'T',     __RC__)
    call MAPL_GetPointer(import, zle,   'ZLE',   __RC__)
    call MAPL_GetPointer(import, zpbl,  'ZPBL',  __RC__)

!   call MAPL_GetPointer(import, h2o2,  'H2O2',     __RC__)
    call MAPL_GetPointer(import, h2so4, 'H2SO4',    __RC__)
    call MAPL_GetPointer(import, so2,   'SO2',      __RC__)
!   call MAPL_GetPointer(import, dms,   'DMS',      __RC__)
    call MAPL_GetPointer(import, nh3,   'NH3',      __RC__)
    call MAPL_GetPointer(import, soa_g, 'SOA_GAS',  __RC__)

    call MAPL_GetPointer(import, ddt_dms_gas,   'DDT_DMS_gas',     __RC__)
    call MAPL_GetPointer(import, ddt_msa_gas,   'DDT_MSA_gas',     __RC__)
    call MAPL_GetPointer(import, ddt_so2_gas,   'DDT_SO2_gas',     __RC__)
    call MAPL_GetPointer(import, ddt_h2so4_gas, 'DDT_H2SO4_gas',   __RC__)
    call MAPL_GetPointer(import, ddt_nh3_gas,   'DDT_NH3_gas',     __RC__)
    call MAPL_GetPointer(import, ddt_soa_g_gas, 'DDT_SOA_GAS_gas', __RC__)

    call MAPL_GetPointer(import, ddt_dms_aq,    'DDT_DMS_aq',      __RC__)
    call MAPL_GetPointer(import, ddt_msa_aq,    'DDT_MSA_aq',      __RC__)
    call MAPL_GetPointer(import, ddt_so2_aq,    'DDT_SO2_aq',      __RC__)
    call MAPL_GetPointer(import, ddt_h2so4_aq,  'DDT_H2SO4_aq',    __RC__)
    call MAPL_GetPointer(import, ddt_nh3_aq,    'DDT_NH3_aq',      __RC__)
    call MAPL_GetPointer(import, ddt_soa_g_aq,  'DDT_SOA_GAS_aq',  __RC__)


    call MAPL_GetPointer(import, msa_g_,   '_MSA_gas',     __RC__)
    call MAPL_GetPointer(import, dms_g_,   '_DMS_gas',     __RC__)
    call MAPL_GetPointer(import, so2_g_,   '_SO2_gas',     __RC__)
    call MAPL_GetPointer(import, h2so4_g_, '_H2SO4_gas',   __RC__)
    call MAPL_GetPointer(import, nh3_g_,   '_NH3_gas',     __RC__)
    call MAPL_GetPointer(import, soa_g_g_, '_SOA_GAS_gas', __RC__)

    call MAPL_GetPointer(import, so2_a_,    '_SO2_aq',      __RC__)
    call MAPL_GetPointer(import, h2so4_a_,  '_H2SO4_aq',    __RC__)
    call MAPL_GetPointer(import, nh3_a_,    '_NH3_aq',      __RC__)
    call MAPL_GetPointer(import, soa_g_a_,  '_SOA_GAS_aq',  __RC__)


    ! aitken
    call MAPL_GetPointer(internal, ait_a_num, 'NUM_A_AIT', __RC__)
    call MAPL_GetPointer(internal, ait_a_so4, 'SU_A_AIT' , __RC__)
    call MAPL_GetPointer(internal, ait_a_nh4, 'AMM_A_AIT', __RC__)
    call MAPL_GetPointer(internal, ait_a_soa, 'SOA_A_AIT', __RC__)
    call MAPL_GetPointer(internal, ait_a_ncl, 'SS_A_AIT' , __RC__)

    call MAPL_GetPointer(internal, ait_c_num, 'NUM_C_AIT', __RC__)
    call MAPL_GetPointer(internal, ait_c_so4, 'SU_C_AIT' , __RC__)
    call MAPL_GetPointer(internal, ait_c_nh4, 'AMM_C_AIT', __RC__) 
    call MAPL_GetPointer(internal, ait_c_soa, 'SOA_C_AIT', __RC__)
    call MAPL_GetPointer(internal, ait_c_ncl, 'SS_C_AIT' , __RC__)
    
    ! accumulation
    call MAPL_GetPointer(internal, acc_a_num, 'NUM_A_ACC', __RC__)
    call MAPL_GetPointer(internal, acc_a_so4, 'SU_A_ACC' , __RC__)
    call MAPL_GetPointer(internal, acc_a_nh4, 'AMM_A_ACC', __RC__)
    call MAPL_GetPointer(internal, acc_a_soa, 'SOA_A_ACC', __RC__)
    call MAPL_GetPointer(internal, acc_a_pom, 'POM_A_ACC', __RC__)
    call MAPL_GetPointer(internal, acc_a_bc,  'BC_A_ACC' , __RC__)
    call MAPL_GetPointer(internal, acc_a_ncl, 'SS_A_ACC' , __RC__)
    
    call MAPL_GetPointer(internal, acc_c_num, 'NUM_C_ACC', __RC__)
    call MAPL_GetPointer(internal, acc_c_so4, 'SU_C_ACC' , __RC__)
    call MAPL_GetPointer(internal, acc_c_nh4, 'AMM_C_ACC', __RC__)
    call MAPL_GetPointer(internal, acc_c_soa, 'SOA_C_ACC', __RC__)
    call MAPL_GetPointer(internal, acc_c_pom, 'POM_C_ACC', __RC__)
    call MAPL_GetPointer(internal, acc_c_bc,  'BC_C_ACC' , __RC__)
    call MAPL_GetPointer(internal, acc_c_ncl, 'SS_C_ACC' , __RC__)

    ! primary carbon mode
    call MAPL_GetPointer(internal, pcm_a_num, 'NUM_A_PCM', __RC__)
    call MAPL_GetPointer(internal, pcm_a_pom, 'POM_A_PCM', __RC__)
    call MAPL_GetPointer(internal, pcm_a_bc,  'BC_A_PCM',  __RC__)

    call MAPL_GetPointer(internal, pcm_c_num, 'NUM_C_PCM', __RC__)
    call MAPL_GetPointer(internal, pcm_c_pom, 'POM_C_PCM', __RC__)
    call MAPL_GetPointer(internal, pcm_c_bc,  'BC_C_PCM',  __RC__)

    ! fine dust 
    call MAPL_GetPointer(internal, fdu_a_num, 'NUM_A_FDU', __RC__)
    call MAPL_GetPointer(internal, fdu_a_dst, 'DU_A_FDU' , __RC__)
    call MAPL_GetPointer(internal, fdu_a_so4, 'SU_A_FDU' , __RC__)
    call MAPL_GetPointer(internal, fdu_a_nh4, 'AMM_A_FDU', __RC__)

    call MAPL_GetPointer(internal, fdu_c_num, 'NUM_C_FDU', __RC__)
    call MAPL_GetPointer(internal, fdu_c_dst, 'DU_C_FDU' , __RC__)
    call MAPL_GetPointer(internal, fdu_c_so4, 'SU_C_FDU' , __RC__)
    call MAPL_GetPointer(internal, fdu_c_nh4, 'AMM_C_FDU', __RC__)

    ! caorse dust 
    call MAPL_GetPointer(internal, cdu_a_num, 'NUM_A_CDU', __RC__)
    call MAPL_GetPointer(internal, cdu_a_dst, 'DU_A_CDU' , __RC__)
    call MAPL_GetPointer(internal, cdu_a_so4, 'SU_A_CDU' , __RC__)
    call MAPL_GetPointer(internal, cdu_a_nh4, 'AMM_A_CDU', __RC__)

    call MAPL_GetPointer(internal, cdu_c_num, 'NUM_C_CDU', __RC__)
    call MAPL_GetPointer(internal, cdu_c_dst, 'DU_C_CDU' , __RC__)
    call MAPL_GetPointer(internal, cdu_c_so4, 'SU_C_CDU' , __RC__)
    call MAPL_GetPointer(internal, cdu_c_nh4, 'AMM_C_CDU', __RC__)

    ! fine seasalt
    call MAPL_GetPointer(internal, fss_a_num, 'NUM_A_FSS', __RC__)
    call MAPL_GetPointer(internal, fss_a_ncl, 'SS_A_FSS' , __RC__)
    call MAPL_GetPointer(internal, fss_a_so4, 'SU_A_FSS' , __RC__)
    call MAPL_GetPointer(internal, fss_a_nh4, 'AMM_A_FSS', __RC__)

    call MAPL_GetPointer(internal, fss_c_num, 'NUM_C_FSS', __RC__)
    call MAPL_GetPointer(internal, fss_c_ncl, 'SS_C_FSS' , __RC__)
    call MAPL_GetPointer(internal, fss_c_so4, 'SU_C_FSS' , __RC__)
    call MAPL_GetPointer(internal, fss_c_nh4, 'AMM_C_FSS', __RC__)

    ! caorse seasalt 
    call MAPL_GetPointer(internal, css_a_num, 'NUM_A_CSS', __RC__)
    call MAPL_GetPointer(internal, css_a_ncl, 'SS_A_CSS' , __RC__)
    call MAPL_GetPointer(internal, css_a_so4, 'SU_A_CSS' , __RC__)
    call MAPL_GetPointer(internal, css_a_nh4, 'AMM_A_CSS', __RC__)

    call MAPL_GetPointer(internal, css_c_num, 'NUM_C_CSS', __RC__)
    call MAPL_GetPointer(internal, css_c_ncl, 'SS_C_CSS' , __RC__)
    call MAPL_GetPointer(internal, css_c_so4, 'SU_C_CSS' , __RC__)
    call MAPL_GetPointer(internal, css_c_nh4, 'AMM_C_CSS', __RC__)

    ! size
    call MAPL_GetPointer(internal, ait_dgn_dry, 'DGN_DRY_AIT', __RC__)
    call MAPL_GetPointer(internal, ait_dgn_wet, 'DGN_WET_AIT', __RC__)
    call MAPL_GetPointer(internal, acc_dgn_dry, 'DGN_DRY_ACC', __RC__)  
    call MAPL_GetPointer(internal, acc_dgn_wet, 'DGN_WET_ACC', __RC__)
    call MAPL_GetPointer(internal, pcm_dgn_dry, 'DGN_DRY_PCM', __RC__)
    call MAPL_GetPointer(internal, pcm_dgn_wet, 'DGN_WET_PCM', __RC__)
    call MAPL_GetPointer(internal, fss_dgn_dry, 'DGN_DRY_FSS', __RC__)
    call MAPL_GetPointer(internal, fss_dgn_wet, 'DGN_WET_FSS', __RC__)
    call MAPL_GetPointer(internal, fdu_dgn_dry, 'DGN_DRY_FDU', __RC__)
    call MAPL_GetPointer(internal, fdu_dgn_wet, 'DGN_WET_FDU', __RC__)
    call MAPL_GetPointer(internal, css_dgn_dry, 'DGN_DRY_CSS', __RC__)
    call MAPL_GetPointer(internal, css_dgn_wet, 'DGN_WET_CSS', __RC__)
    call MAPL_GetPointer(internal, cdu_dgn_dry, 'DGN_DRY_CDU', __RC__)
    call MAPL_GetPointer(internal, cdu_dgn_wet, 'DGN_WET_CDU', __RC__)

    ! aerosol water
    call MAPL_GetPointer(internal, ait_a_wtr, 'WTR_A_AIT', __RC__)
    call MAPL_GetPointer(internal, acc_a_wtr, 'WTR_A_ACC', __RC__)
    call MAPL_GetPointer(internal, pcm_a_wtr, 'WTR_A_PCM', __RC__)
    call MAPL_GetPointer(internal, fdu_a_wtr, 'WTR_A_FDU', __RC__)
    call MAPL_GetPointer(internal, cdu_a_wtr, 'WTR_A_CDU', __RC__)
    call MAPL_GetPointer(internal, fss_a_wtr, 'WTR_A_FSS', __RC__)
    call MAPL_GetPointer(internal, css_a_wtr, 'WTR_A_CSS', __RC__)


    call MAPL_TimerOn(mgState, '--MICROPHYSICS_POSITIVE',   __RC__)

    where (ait_dgn_dry < tiny(0.0))
        ait_dgn_dry = MAM7_AITKEN_MODE_SIZE
    end where

    where (acc_dgn_dry < tiny(0.0))
        acc_dgn_dry = MAM7_ACCUMULATION_MODE_SIZE
    end where

    where (pcm_dgn_dry < tiny(0.0))
        pcm_dgn_dry = MAM7_PRIMARY_CARBON_MODE_SIZE
    end where

    where (fss_dgn_dry < tiny(0.0))
        fss_dgn_dry = MAM7_FINE_SEASALT_MODE_SIZE
    end where

    where (fdu_dgn_dry < tiny(0.0))
        fdu_dgn_dry = MAM7_FINE_SEASALT_MODE_SIZE
    end where

    where (css_dgn_dry < tiny(0.0))
        css_dgn_dry = MAM7_COARSE_SEASALT_MODE_SIZE
    end where

    where (cdu_dgn_dry < tiny(0.0))
        cdu_dgn_dry = MAM7_COARSE_DUST_MODE_SIZE
    end where


    where (ait_dgn_wet < tiny(0.0))
        ait_dgn_wet = MAM7_AITKEN_MODE_SIZE
    end where

    where (acc_dgn_wet < tiny(0.0))
        acc_dgn_wet = MAM7_ACCUMULATION_MODE_SIZE
    end where

    where (pcm_dgn_wet < tiny(0.0))
        pcm_dgn_wet = MAM7_PRIMARY_CARBON_MODE_SIZE
    end where

    where (fss_dgn_wet < tiny(0.0))
        fss_dgn_wet = MAM7_FINE_SEASALT_MODE_SIZE
    end where

    where (fdu_dgn_wet < tiny(0.0))
        fdu_dgn_wet = MAM7_FINE_SEASALT_MODE_SIZE
    end where

    where (css_dgn_wet < tiny(0.0))
        css_dgn_wet = MAM7_COARSE_SEASALT_MODE_SIZE
    end where

    where (cdu_dgn_wet < tiny(0.0))
        cdu_dgn_wet = MAM7_COARSE_DUST_MODE_SIZE
    end where

    call MAPL_TimerOff(mgState, '--MICROPHYSICS_POSITIVE',   __RC__)

#if (1)
    call MAPL_TimerOn(mgState, '-SIZE',   __RC__)

    call MAPL_TimerOn(mgState, '--SIZE_DRY',   __RC__)
    call MAM_DrySize(self%scheme, import, export, self%qa, self%Da, __RC__)
    call MAPL_TimerOff(mgState, '--SIZE_DRY',  __RC__)

    call MAPL_TimerOn(mgState, '--SIZE_WET',   __RC__)
    call MAM_WetSize(self%scheme, import, export, self%qa, self%Da, __RC__)
    call MAPL_TimerOff(mgState, '--SIZE_WET',   __RC__)

    call MAPL_TimerOff(mgState, '-SIZE',   __RC__)
#endif


    AEROSOL_MICROPHYSICS: if (self%microphysics) then

    do j = 1, jm
        do i = 1, im

            amc_t(ncol, :)    = T(i, j, :)              ! temperature at model levels (K)
            amc_pmid(ncol, :) = 0.5*(ple(i,j,0:lm-1)+ple(i,j,1:lm)) ! pressure at layer center (Pa)
            amc_pdel(ncol, :) = delp(i,j,:)             ! pressure thickness of layer (Pa)
            amc_zm(ncol, :)   = zle(i,j,1:lm)           ! altitude (above ground) at layer center (m)
            amc_pblh(ncol)    = zpbl(i,j)               ! planetary boundary layer depth (m)

            amc_qv(ncol, :)   = Q(i,j,:)                ! specific humidity (kg/kg)
            amc_cld(ncol, :)  = fcld(i,j,:)             ! cloud fraction
            amc_rh(ncol, :)   = RH(i,j,:)               ! relative humidity

            ! current tracer mixing ratios (TMRs)
            
            amc_qqcw(:ncol,:pver,:pcnstxx)            = tiny(0.0) 
            amc_qqcw_precldchem(:ncol,:pver,:pcnstxx) = tiny(0.0)  ! qqcw TMRs before cloud chemistry


            ! units mixing ratios should be 'mol/mol-air' and '#/kmol-air'
            amc_q(ncol,:,i_h2o2)   = tiny(0.0)               ! h2o2
            amc_q(ncol,:,i_h2so4)  = h2so4(i,j,:)
            amc_q(ncol,:,i_so2)    = so2(i,j,:)
            amc_q(ncol,:,i_dms)    = tiny(0.0)               ! dms
            amc_q(ncol,:,i_nh3)    = nh3(i,j,:)
            amc_q(ncol,:,i_soag)   = soa_g(i,j,:)
            ! accumulation mode 
            amc_q(ncol,:,i_so4_a1) = acc_a_so4(i,j,:) * (mw_air / adv_mass(i_so4_a1))
            amc_q(ncol,:,i_nh4_a1) = acc_a_nh4(i,j,:) * (mw_air / adv_mass(i_nh4_a1))
            amc_q(ncol,:,i_pom_a1) = acc_a_pom(i,j,:) * (mw_air / adv_mass(i_pom_a1))
            amc_q(ncol,:,i_soa_a1) = acc_a_soa(i,j,:) * (mw_air / adv_mass(i_soa_a1))
            amc_q(ncol,:,i_bc_a1)  = acc_a_bc(i,j,:)  * (mw_air / adv_mass(i_bc_a1))
            amc_q(ncol,:,i_ncl_a1) = acc_a_ncl(i,j,:) * (mw_air / adv_mass(i_ncl_a1))
            amc_q(ncol,:,i_num_a1) = acc_a_num(i,j,:) *  mw_air
            ! aitken mode
            amc_q(ncol,:,i_so4_a2) = ait_a_so4(i,j,:) * (mw_air / adv_mass(i_so4_a2))
            amc_q(ncol,:,i_nh4_a2) = ait_a_nh4(i,j,:) * (mw_air / adv_mass(i_nh4_a2))
            amc_q(ncol,:,i_soa_a2) = ait_a_soa(i,j,:) * (mw_air / adv_mass(i_soa_a2))
            amc_q(ncol,:,i_ncl_a2) = ait_a_ncl(i,j,:) * (mw_air / adv_mass(i_ncl_a2))
            amc_q(ncol,:,i_num_a2) = ait_a_num(i,j,:) *  mw_air
            ! primary carbon mode 
            amc_q(ncol,:,i_pom_a3) = pcm_a_pom(i,j,:) * (mw_air / adv_mass(i_pom_a3))
            amc_q(ncol,:,i_bc_a3)  = pcm_a_bc(i,j,:)  * (mw_air / adv_mass(i_bc_a3))
            amc_q(ncol,:,i_num_a3) = pcm_a_num(i,j,:) *  mw_air
            ! fine seasalt mode
            amc_q(ncol,:,i_ncl_a4) = fss_a_ncl(i,j,:) * (mw_air / adv_mass(i_ncl_a4))
            amc_q(ncol,:,i_so4_a4) = fss_a_so4(i,j,:) * (mw_air / adv_mass(i_so4_a4))
            amc_q(ncol,:,i_nh4_a4) = fss_a_nh4(i,j,:) * (mw_air / adv_mass(i_nh4_a4))
            amc_q(ncol,:,i_num_a4) = fss_a_num(i,j,:) *  mw_air
            ! fine dust mode
            amc_q(ncol,:,i_dst_a5) = fdu_a_dst(i,j,:) * (mw_air / adv_mass(i_dst_a5))
            amc_q(ncol,:,i_so4_a5) = fdu_a_so4(i,j,:) * (mw_air / adv_mass(i_so4_a5))
            amc_q(ncol,:,i_nh4_a5) = fdu_a_nh4(i,j,:) * (mw_air / adv_mass(i_nh4_a5))
            amc_q(ncol,:,i_num_a5) = fdu_a_num(i,j,:) *  mw_air
            ! coarse seasalt mode
            amc_q(ncol,:,i_ncl_a6) = css_a_ncl(i,j,:) * (mw_air / adv_mass(i_ncl_a6))
            amc_q(ncol,:,i_so4_a6) = css_a_so4(i,j,:) * (mw_air / adv_mass(i_so4_a6))
            amc_q(ncol,:,i_nh4_a6) = css_a_nh4(i,j,:) * (mw_air / adv_mass(i_nh4_a6))
            amc_q(ncol,:,i_num_a6) = css_a_num(i,j,:) *  mw_air
            ! coarse dust mode
            amc_q(ncol,:,i_dst_a7) = cdu_a_dst(i,j,:) * (mw_air / adv_mass(i_dst_a7))
            amc_q(ncol,:,i_so4_a7) = cdu_a_so4(i,j,:) * (mw_air / adv_mass(i_so4_a7))
            amc_q(ncol,:,i_nh4_a7) = cdu_a_nh4(i,j,:) * (mw_air / adv_mass(i_nh4_a7))
            amc_q(ncol,:,i_num_a7) = cdu_a_num(i,j,:) *  mw_air


            amc_q_pregaschem(:ncol,:pver,:pcnstxx) = amc_q      ! q TMRs    before gas-phase chemistry
#if (0)
            ! compute pregaschem using tendencies
            amc_q_pregaschem(ncol,:,i_h2so4) = h2so4(i,j,:) - (ddt_h2so4_gas(i,j,:) + ddt_h2so4_aq(i,j,:))*self%dt
            amc_q_pregaschem(ncol,:,i_so2)   = so2(i,j,:)   - (ddt_so2_gas(i,j,:)   + ddt_so2_aq(i,j,:)  )*self%dt
            amc_q_pregaschem(ncol,:,i_nh3)   = nh3(i,j,:)   - (ddt_nh3_gas(i,j,:)   + ddt_nh3_aq(i,j,:)  )*self%dt
#else
            ! ...or use the pregas exports
            amc_q_pregaschem(ncol,:,i_h2so4) = h2so4_g_(i,j,:)
            amc_q_pregaschem(ncol,:,i_so2)   = so2_g_(i,j,:)
            amc_q_pregaschem(ncol,:,i_nh3)   = nh3_g_(i,j,:)
#endif


            amc_q_precldchem(:ncol,:pver,:pcnstxx) = amc_q      ! q TMRs    before cloud chemistry
#if (0)
            ! compute preaqchem using tendencies
            amc_q_precldchem(ncol,:,i_h2so4)  = h2so4(i,j,:) - (ddt_h2so4_aq(i,j,:))*self%dt
            amc_q_precldchem(ncol,:,i_so2)    = so2(i,j,:)   - (ddt_so2_aq(i,j,:)  )*self%dt
            amc_q_precldchem(ncol,:,i_nh3)    = nh3(i,j,:)   - (ddt_nh3_aq(i,j,:)  )*self%dt
#else
            ! ...or use the preaq exports
            amc_q_precldchem(ncol,:,i_h2so4)  = h2so4_a_(i,j,:)
            amc_q_precldchem(ncol,:,i_so2)    = so2_a_(i,j,:)
            amc_q_precldchem(ncol,:,i_nh3)    = nh3_a_(i,j,:)
#endif
            amc_q_precldchem(ncol,:,i_so4_a1) = acc_a_so4_(i,j,:) * (mw_air / adv_mass(i_so4_a1))
            amc_q_precldchem(ncol,:,i_nh4_a1) = acc_a_nh4_(i,j,:) * (mw_air / adv_mass(i_nh4_a1))
            amc_q_precldchem(ncol,:,i_so4_a2) = ait_a_so4_(i,j,:) * (mw_air / adv_mass(i_so4_a2))
            amc_q_precldchem(ncol,:,i_nh4_a2) = ait_a_nh4_(i,j,:) * (mw_air / adv_mass(i_nh4_a2))
            amc_q_precldchem(ncol,:,i_so4_a4) = fss_a_so4_(i,j,:) * (mw_air / adv_mass(i_so4_a4))
            amc_q_precldchem(ncol,:,i_nh4_a4) = fss_a_nh4_(i,j,:) * (mw_air / adv_mass(i_nh4_a4))
            amc_q_precldchem(ncol,:,i_so4_a5) = fdu_a_so4_(i,j,:) * (mw_air / adv_mass(i_so4_a5))
            amc_q_precldchem(ncol,:,i_nh4_a5) = fdu_a_nh4_(i,j,:) * (mw_air / adv_mass(i_nh4_a5))
            amc_q_precldchem(ncol,:,i_so4_a6) = css_a_so4_(i,j,:) * (mw_air / adv_mass(i_so4_a6))
            amc_q_precldchem(ncol,:,i_nh4_a6) = css_a_nh4_(i,j,:) * (mw_air / adv_mass(i_nh4_a6))
            amc_q_precldchem(ncol,:,i_so4_a7) = cdu_a_so4_(i,j,:) * (mw_air / adv_mass(i_so4_a7))
            amc_q_precldchem(ncol,:,i_nh4_a7) = cdu_a_nh4_(i,j,:) * (mw_air / adv_mass(i_nh4_a7))

            
            amc_dgn_a_dry(pcols,:,1) = acc_dgn_dry(i,j,:)          ! dry geo. mean dia. (m) of number PSD
            amc_dgn_a_dry(pcols,:,2) = ait_dgn_dry(i,j,:)
            amc_dgn_a_dry(pcols,:,3) = pcm_dgn_dry(i,j,:)
            amc_dgn_a_dry(pcols,:,4) = fss_dgn_dry(i,j,:)
            amc_dgn_a_dry(pcols,:,5) = fdu_dgn_dry(i,j,:)
            amc_dgn_a_dry(pcols,:,6) = css_dgn_dry(i,j,:)
            amc_dgn_a_dry(pcols,:,7) = cdu_dgn_dry(i,j,:)

            amc_dgn_a_wet(pcols,:,1) = acc_dgn_wet(i,j,:)          ! wet geo. mean dia. (m) of number PSD
            amc_dgn_a_wet(pcols,:,2) = ait_dgn_wet(i,j,:)
            amc_dgn_a_wet(pcols,:,3) = pcm_dgn_wet(i,j,:)
            amc_dgn_a_wet(pcols,:,4) = fss_dgn_wet(i,j,:)
            amc_dgn_a_wet(pcols,:,5) = fdu_dgn_wet(i,j,:)
            amc_dgn_a_wet(pcols,:,6) = css_dgn_wet(i,j,:)
            amc_dgn_a_wet(pcols,:,7) = cdu_dgn_wet(i,j,:)


            amc_wetdens_host(:pcols,:pver,:ntot_amode) = 1.0e3     ! interstitial aerosol wet density (kg/m3)
            amc_qaerwat(:pcols,:pver,:ntot_amode)      = 0.0       ! optional, aerosol water mixing ratio (kg/kg)
            amc_qaerwat(pcols,:,1) = acc_a_wtr(i,j,:)              ! aerosol water
            amc_qaerwat(pcols,:,2) = ait_a_wtr(i,j,:)
            amc_qaerwat(pcols,:,3) = pcm_a_wtr(i,j,:)
            amc_qaerwat(pcols,:,4) = fss_a_wtr(i,j,:)
            amc_qaerwat(pcols,:,5) = fdu_a_wtr(i,j,:)
            amc_qaerwat(pcols,:,6) = css_a_wtr(i,j,:)
            amc_qaerwat(pcols,:,7) = cdu_a_wtr(i,j,:)

            amc_q_coltendaa        = 0.0d0                         ! column integrated tendencies diagnostics
            amc_qqcw_coltendaa     = 0.0d0                         ! --dito-- but for qqcw


            ! the modal_aero_amicphys_intr() subroutine does in the order listed below:
            !
            ! - in clear grid cells
            !    1. condensation / gas-aerosol-exchange of H2SO4, NH3, H2O
            !    2. renaming after "continuous growth"
            !    3. nucleation (new particle formation)
            !    4. coagulation
            !    5. primary carbon aging
            !
            ! - in cloudy grid cells
            !    1. condensation / gas-aerosol-exchange
            !    2. renaming after "continuous growth"
            !    3. primary carbon aging

            call modal_aero_amicphys_intr(amc_do_gasaerexch,   & 
                                          amc_do_rename,       &
                                          amc_do_newnuc,       &
                                          amc_do_coag,         &
                                          amc_lchnk,           &
                                          ncol,                &
                                          amc_nstep,           &
                                          amc_loffset,         &
                                          amc_deltat,          &
                                          amc_latndx,          &
                                          amc_lonndx,          &
                                          amc_t,               &
                                          amc_pmid,            &
                                          amc_pdel,            &
                                          amc_zm,              &
                                          amc_pblh,            &
                                          amc_qv,              &
                                          amc_cld,             &
                                          amc_rh,              &
                                          amc_q,               &
                                          amc_qqcw,            &
                                          amc_q_pregaschem,    &
                                          amc_q_precldchem,    &
                                          amc_qqcw_precldchem, &
                                          amc_dgn_a_dry,       &
                                          amc_dgn_a_wet,       &
                                          amc_wetdens_host,    &
                                          amc_q_coltendaa,     &
                                          amc_qqcw_coltendaa)! & amc_qaerwat -- optional)


            ! current tracer mixing ratios (TMRs)
!           h2o2             = amc_q(ncol,:,i_h2o2)
            h2so4(i,j,:)     = amc_q(ncol,:,i_h2so4)
            so2(i,j,:)       = amc_q(ncol,:,i_so2)
!           dms              = amc_q(ncol,:,i_dms)
            nh3(i,j,:)       = amc_q(ncol,:,i_nh3)
            soa_g(i,j,:)     = amc_q(ncol,:,i_soag)
            ! accumulation mode
            acc_a_so4(i,j,:) = amc_q(ncol,:,i_so4_a1) * (adv_mass(i_so4_a1) / mw_air)
            acc_a_nh4(i,j,:) = amc_q(ncol,:,i_nh4_a1) * (adv_mass(i_nh4_a1) / mw_air)
            acc_a_pom(i,j,:) = amc_q(ncol,:,i_pom_a1) * (adv_mass(i_pom_a1) / mw_air)
            acc_a_soa(i,j,:) = amc_q(ncol,:,i_soa_a1) * (adv_mass(i_soa_a1) / mw_air)
            acc_a_bc(i,j,:)  = amc_q(ncol,:,i_bc_a1)  * (adv_mass(i_bc_a1)  / mw_air)
            acc_a_ncl(i,j,:) = amc_q(ncol,:,i_ncl_a1) * (adv_mass(i_ncl_a1) / mw_air)
            acc_a_num(i,j,:) = amc_q(ncol,:,i_num_a1) / mw_air
            ! aitken mode
            ait_a_so4(i,j,:) = amc_q(ncol,:,i_so4_a2) * (adv_mass(i_so4_a2) / mw_air)
            ait_a_nh4(i,j,:) = amc_q(ncol,:,i_nh4_a2) * (adv_mass(i_nh4_a2) / mw_air)
            ait_a_soa(i,j,:) = amc_q(ncol,:,i_soa_a2) * (adv_mass(i_soa_a2) / mw_air)
            ait_a_ncl(i,j,:) = amc_q(ncol,:,i_ncl_a2) * (adv_mass(i_ncl_a2) / mw_air)
            ait_a_num(i,j,:) = amc_q(ncol,:,i_num_a2) / mw_air
            ! primary carbon mode 
            pcm_a_pom(i,j,:) = amc_q(ncol,:,i_pom_a3) * (adv_mass(i_pom_a3) / mw_air)
            pcm_a_bc(i,j,:)  = amc_q(ncol,:,i_bc_a3)  * (adv_mass(i_bc_a3)  / mw_air)
            pcm_a_num(i,j,:) = amc_q(ncol,:,i_num_a3) / mw_air
            ! fine seasalt mode
            fss_a_ncl(i,j,:) = amc_q(ncol,:,i_ncl_a4) * (adv_mass(i_ncl_a4) / mw_air)
            fss_a_so4(i,j,:) = amc_q(ncol,:,i_so4_a4) * (adv_mass(i_so4_a4) / mw_air)
            fss_a_nh4(i,j,:) = amc_q(ncol,:,i_nh4_a4) * (adv_mass(i_nh4_a4) / mw_air)
            fss_a_num(i,j,:) = amc_q(ncol,:,i_num_a4) / mw_air
            ! fine dust mode
            fdu_a_dst(i,j,:) = amc_q(ncol,:,i_dst_a5) * (adv_mass(i_dst_a5) / mw_air)
            fdu_a_so4(i,j,:) = amc_q(ncol,:,i_so4_a5) * (adv_mass(i_so4_a5) / mw_air)
            fdu_a_nh4(i,j,:) = amc_q(ncol,:,i_nh4_a5) * (adv_mass(i_nh4_a5) / mw_air)
            fdu_a_num(i,j,:) = amc_q(ncol,:,i_num_a5) / mw_air
            ! coarse seasalt mode
            css_a_ncl(i,j,:) = amc_q(ncol,:,i_ncl_a6) * (adv_mass(i_ncl_a6) / mw_air)
            css_a_so4(i,j,:) = amc_q(ncol,:,i_so4_a6) * (adv_mass(i_so4_a6) / mw_air)
            css_a_nh4(i,j,:) = amc_q(ncol,:,i_nh4_a6) * (adv_mass(i_nh4_a6) / mw_air)
            css_a_num(i,j,:) = amc_q(ncol,:,i_num_a6) / mw_air
            ! coarse dust mode
            cdu_a_dst(i,j,:) = amc_q(ncol,:,i_dst_a7) * (adv_mass(i_dst_a7) / mw_air)
            cdu_a_so4(i,j,:) = amc_q(ncol,:,i_so4_a7) * (adv_mass(i_so4_a7) / mw_air)
            cdu_a_nh4(i,j,:) = amc_q(ncol,:,i_nh4_a7) * (adv_mass(i_nh4_a7) / mw_air)
            cdu_a_num(i,j,:) = amc_q(ncol,:,i_num_a7) / mw_air

            ! save the colmn-integrated diagnostics
            q_coltend_cond_(i,j,:)   = amc_q_coltendaa(ncol,:,iqtend_cond)
            q_coltend_rename_(i,j,:) = amc_q_coltendaa(ncol,:,iqtend_rnam)
            q_coltend_nucl_(i,j,:)   = amc_q_coltendaa(ncol,:,iqtend_nnuc)
            q_coltend_coag_(i,j,:)   = amc_q_coltendaa(ncol,:,iqtend_coag)
            qqcw_coltend_rename_(i,j,:) = amc_q_coltendaa(ncol,:,iqqcwtend_rnam)
        end do
    end do

    end if AEROSOL_MICROPHYSICS

    call MAPL_TimerOff(mgState, '-MICROPHYSICS', __RC__)

    
    nullify(q_coltend); call MAPL_GetPointer(export, q_coltend, 'DDT_NUM_A_ACC_COND', __RC__)
    if (associated(q_coltend)) q_coltend = q_coltend_cond_(:,:,i_num_a1)
    
    nullify(q_coltend); call MAPL_GetPointer(export, q_coltend, 'DDT_NUM_A_AIT_COND', __RC__)
    if (associated(q_coltend)) q_coltend = q_coltend_cond_(:,:,i_num_a2)

    nullify(q_coltend); call MAPL_GetPointer(export, q_coltend, 'DDT_NUM_A_PCM_COND', __RC__)
    if (associated(q_coltend)) q_coltend = q_coltend_cond_(:,:,i_num_a3)

    nullify(q_coltend); call MAPL_GetPointer(export, q_coltend, 'DDT_NUM_A_FSS_COND', __RC__)
    if (associated(q_coltend)) q_coltend = q_coltend_cond_(:,:,i_num_a4)

    nullify(q_coltend); call MAPL_GetPointer(export, q_coltend, 'DDT_NUM_A_FDU_COND', __RC__)
    if (associated(q_coltend)) q_coltend = q_coltend_cond_(:,:,i_num_a5)

    nullify(q_coltend); call MAPL_GetPointer(export, q_coltend, 'DDT_NUM_A_CSS_COND', __RC__)
    if (associated(q_coltend)) q_coltend = q_coltend_cond_(:,:,i_num_a6)

    nullify(q_coltend); call MAPL_GetPointer(export, q_coltend, 'DDT_NUM_A_CDU_COND', __RC__)
    if (associated(q_coltend)) q_coltend = q_coltend_cond_(:,:,i_num_a7)
    

    ! free the memory used to hold the pre-aqueous SO4 and NH4
    deallocate(ait_a_so4_, __STAT__)
    deallocate(ait_a_nh4_, __STAT__)    
    deallocate(acc_a_so4_, __STAT__)
    deallocate(acc_a_nh4_, __STAT__)
    deallocate(fdu_a_so4_, __STAT__)
    deallocate(fdu_a_nh4_, __STAT__)
    deallocate(cdu_a_so4_, __STAT__)
    deallocate(cdu_a_nh4_, __STAT__)
    deallocate(fss_a_so4_, __STAT__)
    deallocate(fss_a_nh4_, __STAT__)
    deallocate(css_a_so4_, __STAT__)
    deallocate(css_a_nh4_, __STAT__)

    ! free the memory used to hold the column-integrated tendency diagnostics
    deallocate(q_coltend_cond_,   __STAT__)
    deallocate(q_coltend_rename_, __STAT__)
    deallocate(q_coltend_coag_,   __STAT__)
    deallocate(q_coltend_nucl_,   __STAT__)
    deallocate(qqcw_coltend_rename_, __STAT__)


!   Emissions:  note that emissions are done after the aerosol microphysics 
!   in order to obtain and apply the pre-gas and pre-aqueous phase mixing ratios
!   required by the later
!   ---------------------

    call MAPL_TimerOn(mgState, '-EMISSIONS', __RC__)

!   Seasalt emissions
!   -----------------
    call MAM_SS_Emission(self%scheme, import, export, self%qa, self%femisSS, self%dt, __RC__)

!   Dust emissions
!   -----------------
    call MAM_DU_Emission(self%scheme, import, export, self%qa, self%femisDU, self%dt, __RC__)

!   Black Carbon emissions
!   ----------------------
    call MAM_BC_Emission(self%scheme, import, export, self%qa, self%dt, __RC__)

!   Organic Carbon emissions
!   ----------------------
    call MAM_OC_Emission(self%scheme, import, export, self%qa, self%pom_oc_ratio, self%dt, __RC__)

!   Sulfate (SO4) emissions
!   -----------------------
    call MAM_SO4_Emission(self%scheme, import, export, self%qa, self%dt, __RC__)

    call MAPL_TimerOff(mgState, '-EMISSIONS', __RC__)


    call MAPL_TimerOn(mgState, '-MODE_MERGING', __RC__)

    AEROSOL_MODE_MERGING: if (self%mode_merging) then

    do j = 1, jm
        do i = 1, im

            ! current tracer mixing ratios (TMRs)
            amc_q(ncol,:,i_h2o2)   = tiny(0.0)
            amc_q(ncol,:,i_h2so4)  = h2so4(i,j,:)
            amc_q(ncol,:,i_so2)    = so2(i,j,:)
            amc_q(ncol,:,i_dms)    = tiny(0.0)               ! dms
            amc_q(ncol,:,i_nh3)    = nh3(i,j,:)
            amc_q(ncol,:,i_soag)   = soa_g(i,j,:)
            ! accumulation mode
            amc_q(ncol,:,i_so4_a1) = acc_a_so4(i,j,:) 
            amc_q(ncol,:,i_nh4_a1) = acc_a_nh4(i,j,:)
            amc_q(ncol,:,i_pom_a1) = acc_a_pom(i,j,:)
            amc_q(ncol,:,i_soa_a1) = acc_a_soa(i,j,:)
            amc_q(ncol,:,i_bc_a1)  = acc_a_bc(i,j,:)
            amc_q(ncol,:,i_ncl_a1) = acc_a_ncl(i,j,:)
            amc_q(ncol,:,i_num_a1) = acc_a_num(i,j,:)
            ! aitken mode
            amc_q(ncol,:,i_so4_a2) = ait_a_so4(i,j,:)
            amc_q(ncol,:,i_nh4_a2) = ait_a_nh4(i,j,:)
            amc_q(ncol,:,i_soa_a2) = ait_a_soa(i,j,:)
            amc_q(ncol,:,i_ncl_a2) = ait_a_ncl(i,j,:)
            amc_q(ncol,:,i_num_a2) = ait_a_num(i,j,:)
            ! primary carbon mode 
            amc_q(ncol,:,i_pom_a3) = pcm_a_pom(i,j,:)
            amc_q(ncol,:,i_bc_a3)  = pcm_a_bc(i,j,:)
            amc_q(ncol,:,i_num_a3) = pcm_a_num(i,j,:)
            ! fine seasalt mode
            amc_q(ncol,:,i_ncl_a4) = fss_a_ncl(i,j,:)
            amc_q(ncol,:,i_so4_a4) = fss_a_so4(i,j,:)
            amc_q(ncol,:,i_nh4_a4) = fss_a_nh4(i,j,:)
            amc_q(ncol,:,i_num_a4) = fss_a_num(i,j,:)
            ! fine dust mode
            amc_q(ncol,:,i_dst_a5) = fdu_a_dst(i,j,:)
            amc_q(ncol,:,i_so4_a5) = fdu_a_so4(i,j,:)
            amc_q(ncol,:,i_nh4_a5) = fdu_a_nh4(i,j,:)
            amc_q(ncol,:,i_num_a5) = fdu_a_num(i,j,:)
            ! coarse seasalt mode
            amc_q(ncol,:,i_ncl_a6) = css_a_ncl(i,j,:)
            amc_q(ncol,:,i_so4_a6) = css_a_so4(i,j,:)
            amc_q(ncol,:,i_nh4_a6) = css_a_nh4(i,j,:)
            amc_q(ncol,:,i_num_a6) = css_a_num(i,j,:)
            ! coarse dust mode
            amc_q(ncol,:,i_dst_a7) = cdu_a_dst(i,j,:)
            amc_q(ncol,:,i_so4_a7) = cdu_a_so4(i,j,:)
            amc_q(ncol,:,i_nh4_a7) = cdu_a_nh4(i,j,:) 
            amc_q(ncol,:,i_num_a7) = cdu_a_num(i,j,:)

            
            amc_qqcw(:ncol,:pver,:pcnstxx) = tiny(0.0) 

            
            amc_dgn_a_dry(pcols,:,1) = acc_dgn_dry(i,j,:)          ! dry geo. mean dia. (m) of number PSD
            amc_dgn_a_dry(pcols,:,2) = ait_dgn_dry(i,j,:)
            amc_dgn_a_dry(pcols,:,3) = pcm_dgn_dry(i,j,:)
            amc_dgn_a_dry(pcols,:,4) = fss_dgn_dry(i,j,:)
            amc_dgn_a_dry(pcols,:,5) = fdu_dgn_dry(i,j,:)
            amc_dgn_a_dry(pcols,:,6) = css_dgn_dry(i,j,:)
            amc_dgn_a_dry(pcols,:,7) = cdu_dgn_dry(i,j,:)

            amc_dgn_a_wet(pcols,:,1) = acc_dgn_wet(i,j,:)          ! wet geo. mean dia. (m) of number PSD
            amc_dgn_a_wet(pcols,:,2) = ait_dgn_wet(i,j,:)
            amc_dgn_a_wet(pcols,:,3) = pcm_dgn_wet(i,j,:)
            amc_dgn_a_wet(pcols,:,4) = fss_dgn_wet(i,j,:)
            amc_dgn_a_wet(pcols,:,5) = fdu_dgn_wet(i,j,:)
            amc_dgn_a_wet(pcols,:,6) = css_dgn_wet(i,j,:)
            amc_dgn_a_wet(pcols,:,7) = cdu_dgn_wet(i,j,:)

            amc_dotend = .false.
            amc_dqdt   = 0.0d0
            call modal_aero_calcsize_sub(1,             &
                                         1, 1, 72, 72,  &
                                         amc_pdel,      &
                                         amc_q,         &
                                         amc_qqcw,      &
                                         amc_dqdt,      &
                                         amc_dgn_a_dry, &
                                         amc_deltat,    &
                                         amc_dotend,    & 
                                         self%verbose)

            ! apply the tendencies due to the mode manager mechanism
            amc_q(ncol,:,7:pcnstxx) = amc_q(ncol,:,7:pcnstxx) + amc_dqdt(ncol,:,7:pcnstxx)*amc_deltat


                        ! current tracer mixing ratios (TMRs)
!           h2o2             = amc_q(ncol,:,i_h2o2)
            h2so4(i,j,:)     = amc_q(ncol,:,i_h2so4)
            so2(i,j,:)       = amc_q(ncol,:,i_so2)
!           dms              = amc_q(ncol,:,i_dms)
            nh3(i,j,:)       = amc_q(ncol,:,i_nh3)
            soa_g(i,j,:)     = amc_q(ncol,:,i_soag)

            ! accumulation mode
            acc_a_so4(i,j,:) = amc_q(ncol,:,i_so4_a1) 
            acc_a_nh4(i,j,:) = amc_q(ncol,:,i_nh4_a1)
            acc_a_pom(i,j,:) = amc_q(ncol,:,i_pom_a1)
            acc_a_soa(i,j,:) = amc_q(ncol,:,i_soa_a1)
            acc_a_bc(i,j,:)  = amc_q(ncol,:,i_bc_a1)
            acc_a_ncl(i,j,:) = amc_q(ncol,:,i_ncl_a1)
            acc_a_num(i,j,:) = amc_q(ncol,:,i_num_a1)
            ! aitken mode
            ait_a_so4(i,j,:) = amc_q(ncol,:,i_so4_a2)
            ait_a_nh4(i,j,:) = amc_q(ncol,:,i_nh4_a2)
            ait_a_soa(i,j,:) = amc_q(ncol,:,i_soa_a2)
            ait_a_ncl(i,j,:) = amc_q(ncol,:,i_ncl_a2)
            ait_a_num(i,j,:) = amc_q(ncol,:,i_num_a2)
            ! primary carbon mode 
            pcm_a_pom(i,j,:) = amc_q(ncol,:,i_pom_a3)
            pcm_a_bc(i,j,:)  = amc_q(ncol,:,i_bc_a3)
            pcm_a_num(i,j,:) = amc_q(ncol,:,i_num_a3)
            ! fine seasalt mode
            fss_a_ncl(i,j,:) = amc_q(ncol,:,i_ncl_a4) 
            fss_a_so4(i,j,:) = amc_q(ncol,:,i_so4_a4) 
            fss_a_nh4(i,j,:) = amc_q(ncol,:,i_nh4_a4)
            fss_a_num(i,j,:) = amc_q(ncol,:,i_num_a4)
            ! fine dust mode
            fdu_a_dst(i,j,:) = amc_q(ncol,:,i_dst_a5)
            fdu_a_so4(i,j,:) = amc_q(ncol,:,i_so4_a5)
            fdu_a_nh4(i,j,:) = amc_q(ncol,:,i_nh4_a5)
            fdu_a_num(i,j,:) = amc_q(ncol,:,i_num_a5)
            ! coarse seasalt mode
            css_a_ncl(i,j,:) = amc_q(ncol,:,i_ncl_a6)
            css_a_so4(i,j,:) = amc_q(ncol,:,i_so4_a6)
            css_a_nh4(i,j,:) = amc_q(ncol,:,i_nh4_a6)
            css_a_num(i,j,:) = amc_q(ncol,:,i_num_a6)
            ! coarse dust mode
            cdu_a_dst(i,j,:) = amc_q(ncol,:,i_dst_a7) 
            cdu_a_so4(i,j,:) = amc_q(ncol,:,i_so4_a7) 
            cdu_a_nh4(i,j,:) = amc_q(ncol,:,i_nh4_a7)
            cdu_a_num(i,j,:) = amc_q(ncol,:,i_num_a7)


#if (0)
            acc_a_wtr(i,j,:) = amc_qaerwat(pcols,:,1)     ! aerosol water
            ait_a_wtr(i,j,:) = amc_qaerwat(pcols,:,2)
            pcm_a_wtr(i,j,:) = amc_qaerwat(pcols,:,3)
            fss_a_wtr(i,j,:) = amc_qaerwat(pcols,:,4)
            fdu_a_wtr(i,j,:) = amc_qaerwat(pcols,:,5)
            css_a_wtr(i,j,:) = amc_qaerwat(pcols,:,6)
            cdu_a_wtr(i,j,:) = amc_qaerwat(pcols,:,7)
#endif

#if (0)
            acc_dgn_dry(i,j,:) = amc_dgn_a_dry(pcols,:,1)  ! dry geo. mean dia. (m) of number PSD
            ait_dgn_dry(i,j,:) = amc_dgn_a_dry(pcols,:,2)
            pcm_dgn_dry(i,j,:) = amc_dgn_a_dry(pcols,:,3)
            fss_dgn_dry(i,j,:) = amc_dgn_a_dry(pcols,:,4)
            fdu_dgn_dry(i,j,:) = amc_dgn_a_dry(pcols,:,5)
            css_dgn_dry(i,j,:) = amc_dgn_a_dry(pcols,:,6)
            cdu_dgn_dry(i,j,:) = amc_dgn_a_dry(pcols,:,7)

            acc_dgn_wet(i,j,:) = amc_dgn_a_wet(pcols,:,1)
            ait_dgn_wet(i,j,:) = amc_dgn_a_wet(pcols,:,2)
            pcm_dgn_wet(i,j,:) = amc_dgn_a_wet(pcols,:,3)
            fss_dgn_wet(i,j,:) = amc_dgn_a_wet(pcols,:,4)
            fdu_dgn_wet(i,j,:) = amc_dgn_a_wet(pcols,:,5)
            css_dgn_wet(i,j,:) = amc_dgn_a_wet(pcols,:,6)
            cdu_dgn_wet(i,j,:) = amc_dgn_a_wet(pcols,:,7)
#endif
        end do
    end do

    end if AEROSOL_MODE_MERGING

    call MAPL_TimerOff(mgState, '-MODE_MERGING', __RC__)


!
!   Dry removal - gravitational settling and dry deposition
!   -------------------------------------------------------

    if (self%dry_removal) then
#if (1)
        call MAPL_TimerOn(mgState, '-SIZE',   __RC__)

        call MAPL_TimerOn(mgState, '--SIZE_DRY', __RC__)
        call MAM_DrySize(self%scheme, import, export, self%qa, self%Da, __RC__)
        call MAPL_TimerOff(mgState, '--SIZE_DRY', __RC__)

        call MAPL_TimerOn(mgState, '--SIZE_WET', __RC__)
        call MAM_WetSize(self%scheme, import, export, self%qa, self%Da, __RC__)
        call MAPL_TimerOff(mgState, '--SIZE_WET', __RC__)

        call MAPL_TimerOff(mgState, '-SIZE',   __RC__)
#endif
        call MAPL_TimerOn(mgState, '-REMOVAL', __RC__)
        call MAPL_TimerOn(mgState, '--REMOVAL_DRY', __RC__)
        call MAM_DryRemoval(self%scheme, import, export, self%qa, self%Da, self%dt, __RC__)
        call MAPL_TimerOff(mgState, '--REMOVAL_DRY', __RC__)
        call MAPL_TimerOff(mgState, '-REMOVAL', __RC__)
    end if


!
!   Wet removal - large scale wet scavenging
!   ----------------------------------------
    call MAPL_TimerOn(mgState, '-REMOVAL',      __RC__) 
    call MAPL_TimerOn(mgState, '--REMOVAL_WET', __RC__)

    if (self%wet_removal) then
        call MAM_WetRemoval(self%scheme, import, export, self%qa, self%dt, __RC__)
    end if

    call MAPL_TimerOff(mgState, '--REMOVAL_WET', __RC__)
    call MAPL_TimerOff(mgState, '-REMOVAL',      __RC__)

!  
!   Update the aerosol size and absorbed water
!   ------------------------------------------
    call MAPL_TimerOn(mgState, '-HYGROSCOPIC_GROWTH', __RC__)
#if (1)
    call MAPL_TimerOn(mgState, '-SIZE',   __RC__)

    call MAPL_TimerOn(mgState, '--SIZE_DRY', __RC__)
    call MAM_DrySize(self%scheme, import, export, self%qa, self%Da, __RC__)
    call MAPL_TimerOff(mgState, '--SIZE_DRY', __RC__)

    call MAPL_TimerOn(mgState, '--SIZE_WET', __RC__)
    call MAM_WetSize(self%scheme, import, export, self%qa, self%Da, __RC__)
    call MAPL_TimerOff(mgState, '--SIZE_WET', __RC__)

    call MAPL_TimerOff(mgState, '-SIZE',   __RC__)
#endif
    call MAPL_TimerOff(mgState, '-HYGROSCOPIC_GROWTH', __RC__)


!   Diagnostics 
!   -----------------
!   NOTE : The order of which the processes are done will have 
!          some impact on the dignostic fields
!   -----------------------------------------------------------
    call MAPL_TimerOn(mgState, '-DIAGNOSTICS', __RC__)
    
    call MAPL_TimerOn(mgState,  '--DIAGNOSTICS_SEASALT', __RC__)
    call MAM_SS_Diagnostics(self%scheme, import, export, self%qa, self%dt, __RC__)
    call MAPL_TimerOff(mgState, '--DIAGNOSTICS_SEASALT', __RC__)

    call MAPL_TimerOn(mgState,  '--DIAGNOSTICS_DUST', __RC__)
    call MAM_DU_Diagnostics(self%scheme, import, export, self%qa, self%dt, __RC__)
    call MAPL_TimerOff(mgState, '--DIAGNOSTICS_DUST', __RC__)

    call MAPL_TimerOn(mgState,  '--DIAGNOSTICS_CIM', __RC__)
    call CIM_Diagnostics(self%scheme, import, export, self%qa, self%dt, __RC__)
    call MAPL_TimerOff(mgState, '--DIAGNOSTICS_CIM', __RC__)

    call MAPL_TimerOn(mgState,  '--DIAGNOSTICS_SFC', __RC__)
    call SFC_Diagnostics(self%scheme, import, export, self%qa, self%Da, self%dt, __RC__)
    call MAPL_TimerOff(mgState, '--DIAGNOSTICS_SFC', __RC__)


    call MAPL_TimerOn(mgState,  '--DIAGNOSTICS_AOT', __RC__)
    call AOT_Diagnostics(self%scheme, import, export, self%qa, self%Da, self%mie_ait, &
                                                                        self%mie_acc, &
                                                                        self%mie_pcm, &
                                                                        self%mie_fss, &
                                                                        self%mie_css, &
                                                                        self%mie_fdu, &
                                                                        self%mie_cdu, &
                                                                        self%dt, __RC__)
    call MAPL_TimerOff(mgState, '--DIAGNOSTICS_AOT', __RC__)

    call MAPL_TimerOff(mgState, '-DIAGNOSTICS', __RC__)

    call MAPL_TimerOff(mgState, 'RUN', __RC__)

    call MAPL_TimerOff(mgState, 'TOTAL', __RC__)

!   All done
!   --------
    _RETURN(ESMF_SUCCESS)

   end subroutine Run_


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Finalize_ --- Finalize MAMchem
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
!  01Dec2009  da Silva  First crack.
!
!EOP
!-------------------------------------------------------------------------

                              __Iam__('Finalize_')
   
    type(MAM_state), pointer      :: self               ! Legacy state
    type(ESMF_Grid)               :: GRID               ! Grid
    type(ESMF_Config)             :: CF                 ! Universal Config 

    integer                       :: im_World, jm_World ! Global 2D dimensions
    integer                       :: im, jm, lm         ! 3D Dimensions

    integer                       :: nymd, nhms         ! date, time
    real                          :: cdt                ! time step in secs

    character(len=ESMF_MAXSTR)    :: COMP_NAME

!  Get my name and set-up traceback handle
!  ---------------------------------------
   call ESMF_GridCompGet( GC, name=comp_name, __RC__ )
   Iam = trim(comp_name) // '::' // trim(Iam)

!  Finalize MAPL Generic
!  ---------------------
   call MAPL_GenericFinalize ( GC, IMPORT, EXPORT, CLOCK,  __RC__ )

!  Extract relevant runtime information
!  ------------------------------------
   call extract_ ( GC, CLOCK, self, GRID, CF, &
                   im_World, jm_World,        &
                   im, jm, lm,                &
                   nymd, nhms, cdt, __RC__ )


!  Delete the tracers bundles
!  ------------------------------------
   call MAPL_SimpleBundleDestroy(self%qa)
   call MAPL_SimpleBundleDestroy(self%qc)
   call MAPL_SimpleBundleDestroy(self%qg)

   call MAPL_SimpleBundleDestroy(self%Da)

!  Delete the broad-band optical tables
!  ------------------------------------
   call MAML_OpticsTableDestroy(MAM7_MieTable(1), __RC__)
   call MAML_OpticsTableDestroy(MAM7_MieTable(2), __RC__)
   call MAML_OpticsTableDestroy(MAM7_MieTable(3), __RC__)
   call MAML_OpticsTableDestroy(MAM7_MieTable(4), __RC__)
   call MAML_OpticsTableDestroy(MAM7_MieTable(5), __RC__)
   call MAML_OpticsTableDestroy(MAM7_MieTable(6), __RC__)
   call MAML_OpticsTableDestroy(MAM7_MieTable(7), __RC__)

!  Delete the narrow-band optical tables
!  ------------------------------------
   call MAML_OpticsTableDestroy(self%mie_ait, __RC__)
   call MAML_OpticsTableDestroy(self%mie_acc, __RC__)
   call MAML_OpticsTableDestroy(self%mie_pcm, __RC__)
   call MAML_OpticsTableDestroy(self%mie_fss, __RC__)
   call MAML_OpticsTableDestroy(self%mie_css, __RC__)
   call MAML_OpticsTableDestroy(self%mie_fdu, __RC__)
   call MAML_OpticsTableDestroy(self%mie_cdu, __RC__)

!  Delete the internal private state
!  ---------------------------------
   deallocate(self, __STAT__)


!  All done
!  --------
   _RETURN(ESMF_SUCCESS)

 end subroutine Finalize_

!.......................................................................

 subroutine extract_ (GC, CLOCK,             &
                          myState, GRID, CF, &
                          im_World, jm_World,&
                          im, jm, lm,        &
                          nymd, nhms,        &
                          cdt, rc)

    type(ESMF_GridComp), intent(inout) :: GC                 ! Grid Comp object
    type(ESMF_Clock), intent(in)       :: CLOCK              ! Clock

    type(MAM_state), pointer           :: myState            ! Legacy state
    type(ESMF_Grid),     intent(out)   :: GRID               ! Grid
    type(ESMF_Config),   intent(out)   :: CF                 ! Universal Config 

    integer, intent(out)               :: im_World, jm_World ! Global 2D Dimensions
    integer, intent(out)               :: im, jm, lm         ! 3D Dimensions

    integer, intent(out)               :: nymd, nhms         ! date, time
    real, intent(out)                  :: cdt                ! time step in secs
    integer, intent(out), optional     :: rc

!                            ---

    __Iam__('extract_')

    character(len=ESMF_MAXSTR)         :: comp_name

    type(MAPL_MetaComp), pointer       :: mgState      ! MAPL generic state
    type(MAM_Wrap)                     :: wrap

    integer, dimension(3)              :: dims

    type(ESMF_Alarm)                   :: run_alarm
    type(ESMF_TimeInterval)            :: ring_interval
    real(ESMF_KIND_R8)                 :: time_step

    type(ESMF_Time)                    :: time
    integer                            :: iyr, imm, idd, ihr, imn, isc


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
    call ESMF_UserCompGetInternalState(GC, 'MAM_state', wrap, STATUS)
    _VERIFY(STATUS)
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

!   Global dimensions
!   -----------------
    call MAPL_GridGet(GRID, globalCellCountPerDim=dims, __RC__)
    im_World = dims(1)
    jm_World = dims(2)

!   Local dimensions
!   ----------------
    call ESMF_GridGet(GRID, localDE=0, &
                            staggerloc=ESMF_STAGGERLOC_CENTER, &
                            computationalCount=dims, __RC__)
    im = dims(1)
    jm = dims(2)
    lm = dims(3)


    _RETURN(ESMF_SUCCESS)

 end subroutine extract_


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  MAM_GetTracerName --- returns short and long names of a tracer
!
! !INTERFACE:
!

   subroutine MAM_GetFieldName (mode_short_name, mode_long_name, species, attachment_state, type, short_name, long_name, rc)

! !USES:

    implicit none

! !INPUT/OUTPUT PARAMETERS:

! !INPUT PARAMETERS:
    character(len=MAM_MAXSTR), intent(in)  :: mode_short_name     ! mode short name
    character(len=MAM_MAXSTR), intent(in)  :: mode_long_name      ! mode long name
    character(len=MAM_MAXSTR), intent(in)  :: species             ! species name/alias
    character(len=MAM_MAXSTR), intent(in)  :: attachment_state    ! attachment state = {'interstitial' | 'cloud-borne'}
    character(len=MAM_MAXSTR), intent(in)  :: type                ! tracer type = {'number' | 'mass'}


! !OUTPUT PARAMETERS:
    character(len=MAM_MAXSTR), intent(out) :: short_name          ! short name of tracer
    character(len=MAM_MAXSTR), intent(out) :: long_name           ! long name of tracer

    integer, intent(out)                   :: rc                  ! error return code:
                                                                  !    0 - all is well
                                                                  !    1 -
 
! !DESCRIPTION: This routine constructs short and long names of aerosol tracer 
!               in one of the MAM modes depending on its type and 
!               attachment state.
!
! !REVISION HISTORY:
!
!  6 May 2014    A. Darmenov    Initial implementation.
!
!EOP
!-------------------------------------------------------------------------

    character(len=*), parameter :: Iam = 'MAM_GetTracerName'

    integer :: STATUS

    character(len=3)          :: state
    character(len=MAM_MAXSTR) :: name, buff


!   Initialize local variables
!   --------------------------
    rc = ESMF_SUCCESS

    _ASSERT(attachment_state == 'interstitial' .or. attachment_state == 'cloud-borne','needs informative message')
    _ASSERT(type == 'number' .or. type == 'mass','needs informative message')
    

    if (type == 'number') then
        name = 'NUM'
        buff = 'number of '
    else
        name = trim(species)
        buff = 'mass mixing ratio of '
    end if

    if (attachment_state == 'interstitial') then
        state = '_A_'
    else
        state = '_C_'
    end if
            
   
    short_name = trim(name) // state // trim(mode_short_name)
    long_name  = buff // trim(attachment_state) // ' aerosol particles in ' // trim(mode_long_name) // ' mode'


    _RETURN(ESMF_SUCCESS)

 end subroutine MAM_GetFieldName


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  AqueousChemistry --- Simplified treatment of aqueos chemistry
!
! !INTERFACE:
!

   subroutine AqueousChemistry (self, import, export, qa, cdt, rc)

! !USES:

    use Chem_ConstMod, only: grav

    implicit none

! !INPUT/OUTPUT PARAMETERS:
    type(MAPL_SimpleBundle), intent(inout) :: qa         ! interstitial aerosol tracer fields
    type(ESMF_State), intent(inout)        :: export     ! export fields

! !INPUT PARAMETERS:
    type(MAM_Scheme), intent(in)           :: self       ! MAM scheme/configuration

    type(ESMF_State), intent(inout)        :: import     ! import fields
    real,    intent(in)                    :: cdt        ! chemical timestep (secs)

! !OUTPUT PARAMETERS:
    integer, intent(out)                   :: rc         ! error return code:
                                                         !    0 - all is well
                                                         !    1 -
 
! !DESCRIPTION: This routine updates mass fields due to aqueos chemistry processes.
!
! !REVISION HISTORY:
!
!  2 June 2014    A. Darmenov    Initial implementation.
!
!EOP
!-------------------------------------------------------------------------

    character(len=*), parameter :: Iam = 'AqueousChemistry'

    integer :: STATUS
    integer :: i1, i2, j1, j2, k1, km, i, k, m, s
    integer :: ijl, ijkl

!   Input fields from fvGCM
!   -----------------------
    real, pointer, dimension(:,:,:) :: rhoa, ple, delp
    real, pointer, dimension(:,:,:) :: pSO4_aq, pNH4_aq

!   Exports - diagnostic fields
!   ---------------------------

!   Local
!   -----
    real, allocatable, dimension(:,:,:) :: num_aq       ! total number mixing ratio for modes affected by aqueos chemistry
    real, allocatable, dimension(:,:,:) :: f            ! fraction of aqueos production
    integer :: in_acc, in_fss, in_css, in_fdu, in_cdu   ! index of number mixing raio
    integer :: iq_acc, iq_fss, iq_css, iq_fdu, iq_cdu   ! index of mass mixing raio

!   Parameters 
!   ----------
    real, parameter :: mw_air = 28.965           ! molar mass of dry air, g mol-1
    real, parameter :: mw_SO4 = 96.07            ! molar mass of sulfate, g mol-1
    real, parameter :: mw_NH4 = 18.0385          ! molar mass of dry air, g mol-1


!   Initialize local variables
!   --------------------------
    rc = ESMF_SUCCESS

    _ASSERT(self%id == MAM7_SCHEME,'needs informative message') 
 

!   Get Imports
!   --------------
    call MAPL_GetPointer(import, rhoa,    'AIRDENS',   __RC__)
    call MAPL_GetPointer(import, ple,     'PLE',       __RC__)
    call MAPL_GetPointer(import, delp,    'DELP',      __RC__)

    call MAPL_GetPointer(import, pSO4_aq, 'pSO4_aq',   __RC__)
    call MAPL_GetPointer(import, pNH4_aq, 'pNH4_aq',   __RC__)

    
    if ((.not. associated(pSO4_aq)) .or. (.not. associated(pNH4_aq))) then
        print *, 'Skipping MAM::AqueousChemistry()'
        _RETURN(ESMF_SUCCESS)
    end if


!   Get Exports
!   -----------
    

!   Local dimensions
!   ----------------
    i1 = lbound(rhoa, 1)
    i2 = ubound(rhoa, 1)
    j1 = lbound(rhoa, 2)
    j2 = ubound(rhoa, 2)
    k1 = lbound(rhoa, 3)
    km = ubound(rhoa, 3)

    ijl  = (i2 - i1 + 1) * (j2 - j1 + 1)
    ijkl = ijl * km

#ifdef DEBUG
    call write_parallel(trim(Iam) // '::DEBUG::Indexes:')
    call write_parallel((/i1, i2/), format='(("i1, i2 = ", (X2I3)))')
    call write_parallel((/j1, j2/), format='(("j1, j2 = ", (XI3)))')
    call write_parallel((/k1, km/), format='(("k1, k2 = ", (XI3)))')

    call write_parallel(trim(Iam) // '::DEBUG::Inputs:')
    call write_parallel(self%id,    format='(("model = ", (I5)))')
#endif


    in_acc  = MAPL_SimpleBundleGetIndex(qa, 'NUM_A_ACC', 3, __RC__)
    in_fss  = MAPL_SimpleBundleGetIndex(qa, 'NUM_A_FSS', 3, __RC__)
    in_css  = MAPL_SimpleBundleGetIndex(qa, 'NUM_A_CSS', 3, __RC__)
    in_fdu  = MAPL_SimpleBundleGetIndex(qa, 'NUM_A_FDU', 3, __RC__)
    in_cdu  = MAPL_SimpleBundleGetIndex(qa, 'NUM_A_CDU', 3, __RC__)

    allocate(num_aq(i1:i2,j1:j2,km), __STAT__)
    allocate(f(i1:i2,j1:j2,km), __STAT__)

    num_aq = 0.0
    num_aq = ( qa%r3(in_acc)%q + &
               qa%r3(in_fss)%q + &
               qa%r3(in_css)%q + &
               qa%r3(in_fdu)%q + &
               qa%r3(in_cdu)%q )


    ! partition SO4
    iq_acc  = MAPL_SimpleBundleGetIndex(qa, 'SU_A_ACC',  3, __RC__)
    iq_fss  = MAPL_SimpleBundleGetIndex(qa, 'SU_A_FSS',  3, __RC__)
    iq_css  = MAPL_SimpleBundleGetIndex(qa, 'SU_A_CSS',  3, __RC__)
    iq_fdu  = MAPL_SimpleBundleGetIndex(qa, 'SU_A_FDU',  3, __RC__)
    iq_cdu  = MAPL_SimpleBundleGetIndex(qa, 'SU_A_CDU',  3, __RC__)

    ! this is not mass conservative if num_aq = 0, and needs to be revisited!!!
    ! i.e., if num_aq is small but the production from aq. processes is > 0,
    ! this mass will be lost. a better way is to perhaps create a new particle or
    ! and add the mass.

    f = 0.0
    where (num_aq * rhoa > 1.0e-3)  ! aerosol particle concentration larger than 1e-3 #/m-3
        f = (pSO4_aq * cdt) / num_aq
    end where

    qa%r3(iq_acc)%q = qa%r3(iq_acc)%q + (f * qa%r3(in_acc)%q)
    qa%r3(iq_fss)%q = qa%r3(iq_fss)%q + (f * qa%r3(in_fss)%q)
    qa%r3(iq_css)%q = qa%r3(iq_css)%q + (f * qa%r3(in_css)%q)
    qa%r3(iq_fdu)%q = qa%r3(iq_fdu)%q + (f * qa%r3(in_fdu)%q)
    qa%r3(iq_cdu)%q = qa%r3(iq_cdu)%q + (f * qa%r3(in_cdu)%q)

 
    ! partition NH4
    iq_acc = MAPL_SimpleBundleGetIndex(qa, 'AMM_A_ACC', 3, __RC__)
    iq_fss = MAPL_SimpleBundleGetIndex(qa, 'AMM_A_FSS', 3, __RC__)
    iq_css = MAPL_SimpleBundleGetIndex(qa, 'AMM_A_CSS', 3, __RC__)
    iq_fdu = MAPL_SimpleBundleGetIndex(qa, 'AMM_A_FDU', 3, __RC__)
    iq_cdu = MAPL_SimpleBundleGetIndex(qa, 'AMM_A_CDU', 3, __RC__)

    f = 0.0
    where (num_aq * rhoa > 1.0e-3)
        f = (pNH4_aq * cdt) / num_aq
    end where

    qa%r3(iq_acc)%q = qa%r3(iq_acc)%q + (f * qa%r3(in_acc)%q)
    qa%r3(iq_fss)%q = qa%r3(iq_fss)%q + (f * qa%r3(in_fss)%q)
    qa%r3(iq_css)%q = qa%r3(iq_css)%q + (f * qa%r3(in_css)%q)
    qa%r3(iq_fdu)%q = qa%r3(iq_fdu)%q + (f * qa%r3(in_fdu)%q)
    qa%r3(iq_cdu)%q = qa%r3(iq_cdu)%q + (f * qa%r3(in_cdu)%q)


    deallocate(num_aq, __STAT__)
    deallocate(f,      __STAT__)

    _RETURN(ESMF_SUCCESS)

 end subroutine AqueousChemistry


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  CIM_Diagnostics --- Column integrated mass diagnostics
!
! !INTERFACE:
!

   subroutine CIM_Diagnostics (self, import, export, qa, cdt, rc)

! !USES:

    use Chem_ConstMod, only: grav

    implicit none

! !INPUT/OUTPUT PARAMETERS:
    type(MAPL_SimpleBundle), intent(in)    :: qa         ! interstitial aerosol tracer fields
    type(ESMF_State), intent(inout)        :: export     ! export fields

! !INPUT PARAMETERS:
    type(MAM_Scheme), intent(in)           :: self       ! MAM scheme/configuration

    type(ESMF_State), intent(inout)        :: import     ! import fields
    real,    intent(in)                    :: cdt        ! chemical timestep (secs)

! !OUTPUT PARAMETERS:
    integer, intent(out)                   :: rc         ! error return code:
                                                         !    0 - all is well
                                                         !    1 -
 
! !DESCRIPTION: This routine computes column integrated (dry) mass fields.
!
! !REVISION HISTORY:
!
!  8 Mar 2013    A. Darmenov    Initial implementation.
!
!EOP
!-------------------------------------------------------------------------

    character(len=*), parameter :: Iam = 'CIM_Diagnostics'

    integer :: STATUS
    integer :: i1, i2, j1, j2, k1, km, i, k, m, s
    integer :: ijl, ijkl

!   Mode parameters
!   ------------------------
    character(len=MAM_MAXSTR)     :: mode_name      ! aerosol mode name
    character(len=MAM_MAXSTR)     :: species_name   ! aerosol species name
    character(len=MAM_MAXSTR)     :: field_name     ! field name
    integer                       :: n_species      ! number of aerosol species


!   Input fields from fvGCM
!   -----------------------
    real, pointer, dimension(:,:,:) :: rhoa, ple, delp

!   Exports - diagnostic fields
!   ---------------------------
    real, pointer, dimension(:,:)   :: colmass      ! column integrated mass density, kg m-2


!   Initialize local variables
!   --------------------------
    rc = ESMF_SUCCESS

    _ASSERT(self%id == MAM7_SCHEME .or. self%id == MAM3_SCHEME,'needs informative message') 
 

!   Get Imports
!   --------------
    call MAPL_GetPointer(import, rhoa,      'AIRDENS',   __RC__)
    call MAPL_GetPointer(import, ple,       'PLE',       __RC__)
    call MAPL_GetPointer(import, delp,      'DELP',      __RC__)

!   Get Exports
!   --------------


!   Local dimensions
!   ----------------
    i1 = lbound(rhoa, 1)
    i2 = ubound(rhoa, 1)
    j1 = lbound(rhoa, 2)
    j2 = ubound(rhoa, 2)
    k1 = lbound(rhoa, 3)
    km = ubound(rhoa, 3)

    ijl  = (i2 - i1 + 1) * (j2 - j1 + 1)
    ijkl = ijl * km

#ifdef DEBUG
    call write_parallel(trim(Iam) // '::DEBUG::Indexes:')
    call write_parallel((/i1, i2/), format='(("i1, i2 = ", (X2I3)))')
    call write_parallel((/j1, j2/), format='(("j1, j2 = ", (XI3)))')
    call write_parallel((/k1, km/), format='(("k1, k2 = ", (XI3)))')

    call write_parallel(trim(Iam) // '::DEBUG::Inputs:')
    call write_parallel(self%id,    format='(("model = ", (I5)))')
#endif


    do m = 1, self%n_modes
        call MAM_AerosolModeGet(self%mode(m), name=mode_name, n_species=n_species)

        call MAPL_GetPointer(export, colmass,  trim(mode_name) // 'CMASS',  __RC__)
        if (associated(colmass)) colmass = 0.0

        ! mass mixing ratios
        do s = 1, n_species
            species_name = self%mode(m)%species(s)%name

            field_name  = trim(species_name) // '_A_' // trim(mode_name)

            i = MAPL_SimpleBundleGetIndex(qa, trim(field_name), 3, __RC__)
           
            if (associated(colmass)) then
                do k = 1, km
                    colmass(:,:) = colmass(:,:) + qa%r3(i)%q(:,:,k) * delp(:,:,k)/grav
                end do
            end if
        end do
    end do

    _RETURN(ESMF_SUCCESS)

 end subroutine CIM_Diagnostics

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  SFC_Diagnostics --- Near-surface diagnostics
!
! !INTERFACE:
!

   subroutine SFC_Diagnostics (self, import, export, qa, Da, cdt, rc)

! !USES:

    use Chem_ConstMod, only: grav

    implicit none

! !INPUT/OUTPUT PARAMETERS:
    type(MAPL_SimpleBundle), intent(in)    :: qa         ! interstitial aerosol tracer fields
    type(MAPL_SimpleBundle), intent(in)    :: Da         ! interstitial aerosol size fields
    type(ESMF_State), intent(inout)        :: export     ! export fields

! !INPUT PARAMETERS:
    type(MAM_Scheme), intent(in)           :: self       ! MAM scheme/configuration

    type(ESMF_State), intent(inout)        :: import     ! import fields
    real,    intent(in)                    :: cdt        ! chemical timestep (secs)

! !OUTPUT PARAMETERS:
    integer, intent(out)                   :: rc         ! error return code:
                                                         !    0 - all is well
                                                         !    1 -
 
! !DESCRIPTION: This routine computes column integrated (dry) mass fields.
!
! !REVISION HISTORY:
!
!  12 Jun 2015    A. Darmenov    Initial implementation.
!
!EOP
!-------------------------------------------------------------------------

    character(len=*), parameter :: Iam = 'SFC_Diagnostics'

    integer :: STATUS
    integer :: i1, i2, j1, j2, k1, km, i, m, s
    integer :: ijl, ijkl

!   Mode parameters
!   ------------------------
    character(len=MAM_MAXSTR)     :: mode_name      ! aerosol mode name
    character(len=MAM_MAXSTR)     :: species_name   ! aerosol species name
    character(len=MAM_MAXSTR)     :: field_name     ! field name
    integer                       :: n_species      ! number of aerosol species


!   Input fields from fvGCM
!   -----------------------
    real, pointer, dimension(:,:,:) :: rhoa, delp

!   Exports - diagnostic fields
!   ---------------------------
    real, pointer, dimension(:,:)   :: ptr_2d


!   Initialize local variables
!   --------------------------
    rc = ESMF_SUCCESS

    _ASSERT(self%id == MAM7_SCHEME .or. self%id == MAM3_SCHEME,'needs informative message') 
 

!   Get Imports
!   --------------
    call MAPL_GetPointer(import, rhoa,      'AIRDENS',   __RC__)
    call MAPL_GetPointer(import, delp,      'DELP',      __RC__)

!   Get Exports
!   --------------


!   Local dimensions
!   ----------------
    i1 = lbound(rhoa, 1)
    i2 = ubound(rhoa, 1)
    j1 = lbound(rhoa, 2)
    j2 = ubound(rhoa, 2)
    k1 = lbound(rhoa, 3)
    km = ubound(rhoa, 3)

    ijl  = (i2 - i1 + 1) * (j2 - j1 + 1)
    ijkl = ijl * km

#ifdef DEBUG
    call write_parallel(trim(Iam) // '::DEBUG::Indexes:')
    call write_parallel((/i1, i2/), format='(("i1, i2 = ", (X2I3)))')
    call write_parallel((/j1, j2/), format='(("j1, j2 = ", (XI3)))')
    call write_parallel((/k1, km/), format='(("k1, k2 = ", (XI3)))')

    call write_parallel(trim(Iam) // '::DEBUG::Inputs:')
    call write_parallel(self%id,    format='(("model = ", (I5)))')
#endif


    ! surface number concentrations in '# m-3'
    do m = 1, self%n_modes
        call MAM_AerosolModeGet(self%mode(m), name=mode_name, n_species=n_species)

        ! surface number concentrations in '# m-3'
        call MAPL_GetPointer(export, ptr_2d,  'SFC_NUM_' // trim(mode_name),  __RC__)
        if (associated(ptr_2d)) then
            field_name  = 'NUM_A_' // trim(mode_name)
            i = MAPL_SimpleBundleGetIndex(qa, trim(field_name), 3, __RC__)
           
            ptr_2d(:,:) = qa%r3(i)%q(:,:,km) * rhoa(:,:,km)
        end if
    end do

    ! number concentrations in 'kg m-3'
    do m = 1, self%n_modes
        call MAM_AerosolModeGet(self%mode(m), name=mode_name, n_species=n_species)

        call MAPL_GetPointer(export, ptr_2d,  'SFC_WTR_' // trim(mode_name),  __RC__)
        if (associated(ptr_2d)) then
            field_name  = 'WTR_A_' // trim(mode_name)
            i = MAPL_SimpleBundleGetIndex(qa, trim(field_name), 3, __RC__)
           
            ptr_2d(:,:) = qa%r3(i)%q(:,:,km) * rhoa(:,:,km)
        end if
    end do


    ! dry size in 'm'
    do m = 1, self%n_modes
        call MAM_AerosolModeGet(self%mode(m), name=mode_name, n_species=n_species)

        call MAPL_GetPointer(export, ptr_2d,  'SFC_DGN_DRY_' // trim(mode_name),  __RC__)
        if (associated(ptr_2d)) then
            field_name  = 'DGN_DRY_' // trim(mode_name)
            i = MAPL_SimpleBundleGetIndex(Da, trim(field_name), 3, __RC__)
           
            ptr_2d(:,:) = Da%r3(i)%q(:,:,km)
        end if
    end do

    ! wet size in 'm'
    do m = 1, self%n_modes
        call MAM_AerosolModeGet(self%mode(m), name=mode_name, n_species=n_species)

        call MAPL_GetPointer(export, ptr_2d,  'SFC_DGN_WET_' // trim(mode_name),  __RC__)
        if (associated(ptr_2d)) then
            field_name  = 'DGN_WET_' // trim(mode_name)
            i = MAPL_SimpleBundleGetIndex(Da, trim(field_name), 3, __RC__)
           
            ptr_2d(:,:) = Da%r3(i)%q(:,:,km)
        end if
    end do


    _RETURN(ESMF_SUCCESS)

 end subroutine SFC_Diagnostics


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  AOT_Diagnostics --- Aerosol Optical Thickness
!
! !INTERFACE:
!
   subroutine AOT_Diagnostics (self, import, export, qa, Da, mie_ait, mie_acc, mie_pcm, mie_fss, mie_css, mie_fdu, mie_cdu, cdt, rc)

! !USES:

    use Chem_ConstMod, only: grav
    use MAM_ComponentsDataMod, only : MAM_WATER_COMPONENT_DENSITY,        &
                                      MAM_SULFATE_COMPONENT_DENSITY,      &
                                      MAM_AMMONIUM_COMPONENT_DENSITY,     &
                                      MAM_BLACK_CARBON_COMPONENT_DENSITY, &
                                      MAM_DUST_COMPONENT_DENSITY,         &
                                      MAM_SEASALT_COMPONENT_DENSITY,      &
                                      MAM_SOA_COMPONENT_DENSITY,          &
                                      MAM_POM_COMPONENT_DENSITY

    implicit none

! !INPUT/OUTPUT PARAMETERS:
    type(MAPL_SimpleBundle), intent(inout) :: qa         ! interstitial aerosol tracer fields
    type(MAPL_SimpleBundle), intent(inout) :: Da         ! dry and wet sizes of interstital aerosols
    type(ESMF_State), intent(inout)        :: export     ! export fields

! !INPUT PARAMETERS:
    type(MAM_Scheme), intent(in)           :: self       ! MAM scheme/configuration

    type(ESMF_State), intent(inout)        :: import     ! import fields

    type(MAML_OpticsTable)                 :: mie_ait
    type(MAML_OpticsTable)                 :: mie_acc
    type(MAML_OpticsTable)                 :: mie_pcm
    type(MAML_OpticsTable)                 :: mie_fss
    type(MAML_OpticsTable)                 :: mie_css
    type(MAML_OpticsTable)                 :: mie_fdu
    type(MAML_OpticsTable)                 :: mie_cdu

    real,    intent(in)                    :: cdt        ! chemical timestep (secs)

! !OUTPUT PARAMETERS:
    integer, intent(out)                   :: rc         ! error return code:
                                                         !    0 - all is well
                                                         !    1 -
 
! !DESCRIPTION: This routine computes aerosol optical thickness.
!
! !REVISION HISTORY:
!
!  23 May 2014    A. Darmenov    Initial implementation.
!
!EOP
!-------------------------------------------------------------------------

!   Input fields from fvGCM
!   -----------------------
    real, pointer, dimension(:,:,:)         :: delp

!   Exports - diagnostic fields
!   ---------------------------
    real, pointer, dimension(:,:)           :: tot_ext, tot_sca         ! aerosol optical thickness, '1'
    real, pointer, dimension(:,:)           :: ait_ext, ait_sca
    real, pointer, dimension(:,:)           :: acc_ext, acc_sca
    real, pointer, dimension(:,:)           :: pcm_ext, pcm_sca
    real, pointer, dimension(:,:)           :: fss_ext, fss_sca
    real, pointer, dimension(:,:)           :: css_ext, css_sca
    real, pointer, dimension(:,:)           :: fdu_ext, fdu_sca
    real, pointer, dimension(:,:)           :: cdu_ext, cdu_sca

    real, dimension(:,:,:), pointer         :: rh

    real, dimension(:,:,:,:), allocatable   :: qa_

    integer                                 :: i_mmr
    integer                                 :: i_dwet

    real(kind=8), dimension(:,:,:), allocatable :: ext, sca, asy, ssa   ! total, (lon:,lat:,lev:,band:)

    real, dimension(:,:,:), allocatable     :: ext_, sca_, asy_         ! mode,  (lon:,lat:,lev:,band:)
    real, dimension(:),     allocatable     :: density
    integer                                 :: nc

    integer                                 :: i1, j1, i2, j2, k1, km
    integer                                 :: ijl, ijkl

    integer                                 :: band

  
    integer                                 :: STATUS
    character(len=ESMF_MAXSTR)              :: Iam



                                   Iam = 'AOD_Diagnostics'

!   Get Imports
!   --------------
    call MAPL_GetPointer(import, delp,    'DELP',      __RC__)

!   Get Exports
!   --------------
    call MAPL_GetPointer(export, tot_ext, 'TOTEXTTAU', __RC__)
    call MAPL_GetPointer(export, tot_sca, 'TOTSCATAU', __RC__)

    call MAPL_GetPointer(export, ait_ext, 'AITEXTTAU', __RC__)
    call MAPL_GetPointer(export, ait_sca, 'AITSCATAU', __RC__)

    call MAPL_GetPointer(export, acc_ext, 'ACCEXTTAU', __RC__)
    call MAPL_GetPointer(export, acc_sca, 'ACCSCATAU', __RC__)

    call MAPL_GetPointer(export, pcm_ext, 'PCMEXTTAU', __RC__)
    call MAPL_GetPointer(export, pcm_sca, 'PCMSCATAU', __RC__)

    call MAPL_GetPointer(export, fss_ext, 'FSSEXTTAU', __RC__)
    call MAPL_GetPointer(export, fss_sca, 'FSSSCATAU', __RC__)

    call MAPL_GetPointer(export, css_ext, 'CSSEXTTAU', __RC__)
    call MAPL_GetPointer(export, css_sca, 'CSSSCATAU', __RC__)

    call MAPL_GetPointer(export, fdu_ext, 'FDUEXTTAU', __RC__)
    call MAPL_GetPointer(export, fdu_sca, 'FDUSCATAU', __RC__)

    call MAPL_GetPointer(export, cdu_ext, 'CDUEXTTAU', __RC__)
    call MAPL_GetPointer(export, cdu_sca, 'CDUSCATAU', __RC__)

!   Local dimensions
!   ----------------
    i1 = lbound(delp, 1)
    i2 = ubound(delp, 1)
    j1 = lbound(delp, 2)
    j2 = ubound(delp, 2)
    k1 = lbound(delp, 3)
    km = ubound(delp, 3)

    ijl  = (i2 - i1 + 1) * (j2 - j1 + 1)
    ijkl = ijl * km


    ! Radiation band
    ! --------------
    band = 7   ! 550nm 

    ! Pressure at layer edges 
    ! ------------------------



  allocate(ext(i1:i2,j1:j2,km), &
           sca(i1:i2,j1:j2,km), &
           ssa(i1:i2,j1:j2,km), &
           asy(i1:i2,j1:j2,km), __STAT__)

  allocate(ext_(i1:i2,j1:j2,km), &
           sca_(i1:i2,j1:j2,km), &
           asy_(i1:i2,j1:j2,km), __STAT__)
  

  ext = 0.0
  sca = 0.0
  ssa = 0.0
  asy = 0.0


  ! compute ext, sca and ssa from aerosols in aitken mode; su, amm, soa, ss
  nc = 4 + 1
  allocate(qa_(nc,i1:i2,j1:j2,km), density(nc), __STAT__)
 
  qa_ = 0.0 
  density = 0.0

  i_mmr = MAPL_SimpleBundleGetIndex(qa, 'WTR_A_AIT', 3, __RC__)
  qa_(1,:,:,:) = qa%r3(i_mmr)%q
  density(1)  = MAM_WATER_COMPONENT_DENSITY

  i_mmr = MAPL_SimpleBundleGetIndex(qa, 'SU_A_AIT',  3, __RC__)
  qa_(2,:,:,:) = qa%r3(i_mmr)%q
  density(2) = MAM_SULFATE_COMPONENT_DENSITY

  i_mmr = MAPL_SimpleBundleGetIndex(qa, 'AMM_A_AIT', 3, __RC__)
  qa_(3,:,:,:) = qa%r3(i_mmr)%q
  density(3)  = MAM_AMMONIUM_COMPONENT_DENSITY

  i_mmr = MAPL_SimpleBundleGetIndex(qa, 'SOA_A_AIT', 3, __RC__)
  qa_(4,:,:,:) = qa%r3(i_mmr)%q
  density(4)  = MAM_SOA_COMPONENT_DENSITY

  i_mmr = MAPL_SimpleBundleGetIndex(qa, 'SS_A_AIT',  3, __RC__)
  qa_(5,:,:,:) = qa%r3(i_mmr)%q
  density(5)  = MAM_SEASALT_COMPONENT_DENSITY

  i_dwet = MAPL_SimpleBundleGetIndex(Da, 'DGN_WET_AIT', 3, __RC__)

  ext_ = 0.0 
  sca_ = 0.0 
  asy_ = 0.0

  call MAML_OpticsInterpolate(mie_ait, band, qa_, density, Da%r3(i_dwet)%q, delp, ext_, sca_, asy_, nc, i1, i2, j1, j2, 1, km, rc)

  if (associated(ait_ext)) ait_ext = sum(ext_, dim=3)
  if (associated(ait_sca)) ait_sca = sum(sca_, dim=3)

  ext = ext + ext_
  sca = sca + sca_
  asy = asy + sca_*asy_

  deallocate(qa_, density, __STAT__)


  ! compute ext, sca and ssa from aerosols in accumulation mode: su, amm, soa, pom, bc, ss
  nc = 6 + 1 
  allocate(qa_(nc,i1:i2,j1:j2,km), density(nc), __STAT__)
 
  qa_ = 0.0 
  density = 0.0

  i_mmr = MAPL_SimpleBundleGetIndex(qa, 'WTR_A_ACC', 3, __RC__)
  qa_(1,:,:,:) = qa%r3(i_mmr)%q
  density(1)  = MAM_WATER_COMPONENT_DENSITY

  i_mmr = MAPL_SimpleBundleGetIndex(qa, 'SU_A_ACC',  3, __RC__)
  qa_(2,:,:,:) = qa%r3(i_mmr)%q
  density(2) = MAM_SULFATE_COMPONENT_DENSITY

  i_mmr = MAPL_SimpleBundleGetIndex(qa, 'AMM_A_ACC', 3, __RC__)
  qa_(3,:,:,:) = qa%r3(i_mmr)%q
  density(3)  = MAM_AMMONIUM_COMPONENT_DENSITY

  i_mmr = MAPL_SimpleBundleGetIndex(qa, 'SOA_A_ACC', 3, __RC__)
  qa_(4,:,:,:) = qa%r3(i_mmr)%q
  density(4)  = MAM_SOA_COMPONENT_DENSITY

  i_mmr = MAPL_SimpleBundleGetIndex(qa, 'POM_A_ACC',  3, __RC__)
  qa_(5,:,:,:) = qa%r3(i_mmr)%q
  density(5)  = MAM_POM_COMPONENT_DENSITY

  i_mmr = MAPL_SimpleBundleGetIndex(qa, 'BC_A_ACC',  3, __RC__)
  qa_(6,:,:,:) = qa%r3(i_mmr)%q
  density(6)  = MAM_BLACK_CARBON_COMPONENT_DENSITY

  i_mmr = MAPL_SimpleBundleGetIndex(qa, 'SS_A_ACC',  3, __RC__)
  qa_(7,:,:,:) = qa%r3(i_mmr)%q
  density(7)  = MAM_SEASALT_COMPONENT_DENSITY

  i_dwet = MAPL_SimpleBundleGetIndex(Da, 'DGN_WET_ACC', 3, __RC__)

  ext_ = 0.0 
  sca_ = 0.0 
  asy_ = 0.0

  call MAML_OpticsInterpolate(mie_acc, band, qa_, density, Da%r3(i_dwet)%q, delp, ext_, sca_, asy_, nc, i1, i2, j1, j2, 1, km, rc)

  if (associated(acc_ext)) acc_ext = sum(ext_, dim=3)
  if (associated(acc_sca)) acc_sca = sum(sca_, dim=3)

  ext = ext + ext_
  sca = sca + sca_
  asy = asy + sca_*asy_

  deallocate(qa_, density, __STAT__)


  ! compute ext, sca and ssa from aerosols in primary carbon mode: pom, bc
  nc = 2 + 1
  allocate(qa_(nc,i1:i2,j1:j2,km), density(nc), __STAT__)
 
  qa_ = 0.0 
  density = 0.0

  i_mmr = MAPL_SimpleBundleGetIndex(qa, 'WTR_A_PCM', 3, __RC__)
  qa_(1,:,:,:) = qa%r3(i_mmr)%q
  density(1)  = MAM_WATER_COMPONENT_DENSITY

  i_mmr = MAPL_SimpleBundleGetIndex(qa, 'POM_A_PCM',  3, __RC__)
  qa_(2,:,:,:) = qa%r3(i_mmr)%q
  density(2) = MAM_POM_COMPONENT_DENSITY

  i_mmr = MAPL_SimpleBundleGetIndex(qa, 'BC_A_PCM', 3, __RC__)
  qa_(3,:,:,:) = qa%r3(i_mmr)%q
  density(3)  = MAM_BLACK_CARBON_COMPONENT_DENSITY

  i_dwet = MAPL_SimpleBundleGetIndex(Da, 'DGN_WET_PCM', 3, __RC__)

  ext_ = 0.0 
  sca_ = 0.0 
  asy_ = 0.0

  call MAML_OpticsInterpolate(mie_pcm, band, qa_, density, Da%r3(i_dwet)%q, delp, ext_, sca_, asy_, nc, i1, i2, j1, j2, 1, km, rc)

  if (associated(pcm_ext)) pcm_ext = sum(ext_, dim=3)
  if (associated(pcm_sca)) pcm_sca = sum(sca_, dim=3)

  ext = ext + ext_
  sca = sca + sca_
  asy = asy + sca_*asy_

  deallocate(qa_, density, __STAT__)


  ! compute ext, sca and ssa from aerosols in fine seasalt mode: su, amm, ss
  nc = 3 + 1
  allocate(qa_(nc,i1:i2,j1:j2,km), density(nc), __STAT__)
 
  qa_ = 0.0 
  density = 0.0

  i_mmr = MAPL_SimpleBundleGetIndex(qa, 'WTR_A_FSS', 3, __RC__)
  qa_(1,:,:,:) = qa%r3(i_mmr)%q
  density(1)  = MAM_WATER_COMPONENT_DENSITY

  i_mmr = MAPL_SimpleBundleGetIndex(qa, 'SU_A_FSS',  3, __RC__)
  qa_(2,:,:,:) = qa%r3(i_mmr)%q
  density(2) = MAM_SULFATE_COMPONENT_DENSITY

  i_mmr = MAPL_SimpleBundleGetIndex(qa, 'AMM_A_FSS', 3, __RC__)
  qa_(3,:,:,:) = qa%r3(i_mmr)%q
  density(3)  = MAM_AMMONIUM_COMPONENT_DENSITY

  i_mmr = MAPL_SimpleBundleGetIndex(qa, 'SS_A_FSS', 3, __RC__)
  qa_(4,:,:,:) = qa%r3(i_mmr)%q
  density(4)  = MAM_SEASALT_COMPONENT_DENSITY

  i_dwet = MAPL_SimpleBundleGetIndex(Da, 'DGN_WET_FSS', 3, __RC__)

  ext_ = 0.0 
  sca_ = 0.0 
  asy_ = 0.0

  call MAML_OpticsInterpolate(mie_fss, band, qa_, density, Da%r3(i_dwet)%q, delp, ext_, sca_, asy_, nc, i1, i2, j1, j2, 1, km, rc)

  if (associated(fss_ext)) fss_ext = sum(ext_, dim=3)
  if (associated(fss_sca)) fss_sca = sum(sca_, dim=3)

  ext = ext + ext_
  sca = sca + sca_
  asy = asy + sca_*asy_

  deallocate(qa_, density, __STAT__)


  ! compute ext, sca and ssa from aerosols in coarse seasalt mode: su, amm, ss; lut = ('water', 'su', 'amm', 'ss')
  nc = 3 + 1
  allocate(qa_(nc,i1:i2,j1:j2,km), density(nc), __STAT__)
 
  qa_ = 0.0 
  density = 0.0

  i_mmr = MAPL_SimpleBundleGetIndex(qa, 'WTR_A_CSS', 3, __RC__)
  qa_(1,:,:,:) = qa%r3(i_mmr)%q
  density(1)  = MAM_WATER_COMPONENT_DENSITY

  i_mmr = MAPL_SimpleBundleGetIndex(qa, 'SU_A_CSS',  3, __RC__)
  qa_(2,:,:,:) = qa%r3(i_mmr)%q
  density(2) = MAM_SULFATE_COMPONENT_DENSITY

  i_mmr = MAPL_SimpleBundleGetIndex(qa, 'AMM_A_CSS', 3, __RC__)
  qa_(3,:,:,:) = qa%r3(i_mmr)%q
  density(3)  = MAM_AMMONIUM_COMPONENT_DENSITY

  i_mmr = MAPL_SimpleBundleGetIndex(qa, 'SS_A_CSS', 3, __RC__)
  qa_(4,:,:,:) = qa%r3(i_mmr)%q
  density(4)  = MAM_SEASALT_COMPONENT_DENSITY

  i_dwet = MAPL_SimpleBundleGetIndex(Da, 'DGN_WET_CSS', 3, __RC__)

  ext_ = 0.0 
  sca_ = 0.0 
  asy_ = 0.0

  call MAML_OpticsInterpolate(mie_css, band, qa_, density, Da%r3(i_dwet)%q, delp, ext_, sca_, asy_, nc, i1, i2, j1, j2, 1, km, rc)

  if (associated(css_ext)) css_ext = sum(ext_, dim=3)
  if (associated(css_sca)) css_sca = sum(sca_, dim=3)

  ext = ext + ext_
  sca = sca + sca_
  asy = asy + sca_*asy_

  deallocate(qa_, density, __STAT__)


  ! compute ext, sca and ssa from aerosols in fine dust mode: su, amm, du; lut = ('water', 'su', 'amm', 'du')
  nc = 3 + 1
  allocate(qa_(nc,i1:i2,j1:j2,km), density(nc), __STAT__)
 
  qa_ = 0.0 
  density = 0.0

  i_mmr = MAPL_SimpleBundleGetIndex(qa, 'WTR_A_FDU', 3, __RC__)
  qa_(1,:,:,:) = qa%r3(i_mmr)%q
  density(1)  = MAM_WATER_COMPONENT_DENSITY

  i_mmr = MAPL_SimpleBundleGetIndex(qa, 'SU_A_FDU',  3, __RC__)
  qa_(2,:,:,:) = qa%r3(i_mmr)%q
  density(2) = MAM_SULFATE_COMPONENT_DENSITY

  i_mmr = MAPL_SimpleBundleGetIndex(qa, 'AMM_A_FDU', 3, __RC__)
  qa_(3,:,:,:) = qa%r3(i_mmr)%q
  density(3)  = MAM_AMMONIUM_COMPONENT_DENSITY

  i_mmr = MAPL_SimpleBundleGetIndex(qa, 'DU_A_FDU', 3, __RC__)
  qa_(4,:,:,:) = qa%r3(i_mmr)%q
  density(4)  = MAM_DUST_COMPONENT_DENSITY

  i_dwet = MAPL_SimpleBundleGetIndex(Da, 'DGN_WET_FDU', 3, __RC__)

  ext_ = 0.0 
  sca_ = 0.0 
  asy_ = 0.0

  call MAML_OpticsInterpolate(mie_fdu, band, qa_, density, Da%r3(i_dwet)%q, delp, ext_, sca_, asy_, nc, i1, i2, j1, j2, 1, km, rc)

  if (associated(fdu_ext)) fdu_ext = sum(ext_, dim=3)
  if (associated(fdu_sca)) fdu_sca = sum(sca_, dim=3)

  ext = ext + ext_
  sca = sca + sca_
  asy = asy + sca_*asy_

  deallocate(qa_, density, __STAT__)


  ! compute ext, sca and ssa from aerosols in coarse dust mode: su, amm, du; lut = ('water', 'su', 'amm', 'du')
  nc = 3 + 1
  allocate(qa_(nc,i1:i2,j1:j2,km), density(nc), __STAT__)
 
  qa_ = 0.0 
  density = 0.0

  i_mmr = MAPL_SimpleBundleGetIndex(qa, 'WTR_A_CDU', 3, __RC__)
  qa_(1,:,:,:) = qa%r3(i_mmr)%q
  density(1)  = MAM_WATER_COMPONENT_DENSITY

  i_mmr = MAPL_SimpleBundleGetIndex(qa, 'SU_A_CDU',  3, __RC__)
  qa_(2,:,:,:) = qa%r3(i_mmr)%q
  density(2) = MAM_SULFATE_COMPONENT_DENSITY

  i_mmr = MAPL_SimpleBundleGetIndex(qa, 'AMM_A_CDU', 3, __RC__)
  qa_(3,:,:,:) = qa%r3(i_mmr)%q
  density(3)  = MAM_AMMONIUM_COMPONENT_DENSITY

  i_mmr = MAPL_SimpleBundleGetIndex(qa, 'DU_A_CDU', 3, __RC__)
  qa_(4,:,:,:) = qa%r3(i_mmr)%q
  density(4)  = MAM_DUST_COMPONENT_DENSITY

  i_dwet = MAPL_SimpleBundleGetIndex(Da, 'DGN_WET_CDU', 3, __RC__)

  ext_ = 0.0 
  sca_ = 0.0 
  asy_ = 0.0

  call MAML_OpticsInterpolate(mie_cdu, band, qa_, density, Da%r3(i_dwet)%q, delp, ext_, sca_, asy_, nc, i1, i2, j1, j2, 1, km, rc)

  if (associated(cdu_ext)) cdu_ext = sum(ext_, dim=3)
  if (associated(cdu_sca)) cdu_sca = sum(sca_, dim=3)

  ext = ext + ext_
  sca = sca + sca_
  asy = asy + sca_*asy_

  deallocate(qa_, density, __STAT__)


  ! total AOT
  if (associated(tot_ext)) tot_ext = sum(ext, dim=3)
  if (associated(tot_sca)) tot_sca = sum(sca, dim=3)


  deallocate(ext,  sca,  ssa,  asy, __STAT__)
  deallocate(ext_, sca_, asy_,      __STAT__)

  _RETURN(ESMF_SUCCESS)

 end subroutine AOT_Diagnostics


!! --


logical function isDataDrivenGC_(GC, rc)
   type(ESMF_GridComp), intent(inout) :: GC
   integer, intent(out)               :: rc
 
!  local 
   character(len=ESMF_MAXSTR)         :: IAm
   integer                            :: STATUS

   integer                            :: i
   character(len=ESMF_MAXSTR)         :: comp_name
   character(len=*), parameter        :: modifier = '.data'

   rc = ESMF_SUCCESS

   call ESMF_GridCompGet(GC, name=comp_name, __RC__)   
   i = index(trim(comp_name), trim(modifier), back=.true.)
 
   if (i > 0) then
       ! lets be strict
       if (comp_name(i:) == modifier) then
           isDataDrivenGC_ = .true.
       else
           isDataDrivenGC_ = .false.
       end if
   else
       isDataDrivenGC_ = .false.
   end if

   _RETURN(ESMF_SUCCESS)

end function isDataDrivenGC_



 subroutine aerosol_optics(state, rc)

  use MAM_ComponentsDataMod, only : MAM_WATER_COMPONENT_DENSITY,        &
                                    MAM_SULFATE_COMPONENT_DENSITY,      &
                                    MAM_AMMONIUM_COMPONENT_DENSITY,     &
                                    MAM_BLACK_CARBON_COMPONENT_DENSITY, &
                                    MAM_DUST_COMPONENT_DENSITY,         &
                                    MAM_SEASALT_COMPONENT_DENSITY,      &
                                    MAM_SOA_COMPONENT_DENSITY,          &
                                    MAM_POM_COMPONENT_DENSITY

  implicit none 



  

! Arguments
! ---------
  type(ESMF_State):: state
  integer, intent(out):: rc


! Local
! ---------
  type(ESMF_FieldBundle)                  :: aerosols

  real, dimension(:,:,:), pointer         :: ple
  real, dimension(:,:,:), pointer         :: rh
  real, dimension(:,:,:), pointer         :: var
  real, dimension(:,:,:), pointer         :: q
  real, dimension(:,:,:), pointer         :: dgn_wet

  real, dimension(:,:,:,:), allocatable   :: qa


  real, dimension(:,:,:), allocatable     :: dp

  character(len=ESMF_MAXSTR)              :: field_name

  real(kind=8), dimension(:,:,:), allocatable :: ext, sca, asy, ssa   ! total, (lon:,lat:,lev:,band:)

  real, dimension(:,:,:), allocatable     :: ext_, sca_, asy_     ! mode,  (lon:,lat:,lev:,band:)
  real, dimension(:),     allocatable     :: density
  integer                                 :: nc

  integer                                 :: i1, j1, i2, j2, km

  integer                                 :: band

  
  integer                                 :: STATUS
  character(len=ESMF_MAXSTR)              :: Iam



                                   Iam = 'MAM::aerosol_optics()'


! Radiation band
! --------------
  call ESMF_AttributeGet(state, name='band_for_aerosol_optics', value=band, __RC__)

! Pressure at layer edges 
! ------------------------
  call ESMF_AttributeGet(state, name='air_pressure_for_aerosol_optics', value=field_name, __RC__)
  call MAPL_GetPointer(state, ple, trim(field_name), __RC__)

  i1 = lbound(ple, 1); i2 = ubound(ple, 1)
  j1 = lbound(ple, 2); j2 = ubound(ple, 2)
                       km = ubound(ple, 3)

! Relative humidity
! -----------------
  call ESMF_AttributeGet(state, name='relative_humidity_for_aerosol_optics', value=field_name, __RC__)
  call MAPL_GetPointer(state, rh, trim(field_name), __RC__)

  i1 = lbound(rh, 1); i2 = ubound(rh, 1)
  j1 = lbound(rh, 2); j2 = ubound(rh, 2)
                      km = ubound(rh, 3)
  
  call ESMF_StateGet(state, 'AEROSOLS', aerosols, __RC__)


  allocate(dp(i1:i2,j1:j2,km), __STAT__)
  dp = (ple(:,:,1:km) - ple(:,:,0:km-1))


  allocate(ext(i1:i2,j1:j2,km), &
           sca(i1:i2,j1:j2,km), &
           ssa(i1:i2,j1:j2,km), &
           asy(i1:i2,j1:j2,km), __STAT__)

  allocate(ext_(i1:i2,j1:j2,km), &
           sca_(i1:i2,j1:j2,km), &
           asy_(i1:i2,j1:j2,km), __STAT__)
  

  ext = 0.0
  sca = 0.0
  ssa = 0.0
  asy = 0.0


  ! compute ext, sca and ssa from aerosols in aitken mode; su, amm, soa, ss
  nc = 4 + 1
  allocate(qa(nc,i1:i2,j1:j2,km), density(nc), __STAT__)
 
  qa = 0.0 
  density = 0.0

  call ESMFL_BundleGetPointerToData(aerosols, 'WTR_A_AIT', q, __RC__)
  qa(1,:,:,:) = q
  density(1)  = MAM_WATER_COMPONENT_DENSITY

  call ESMFL_BundleGetPointerToData(aerosols, 'SU_A_AIT',  q, __RC__)
  qa(2,:,:,:) = q
  density(2)  = MAM_SULFATE_COMPONENT_DENSITY

  call ESMFL_BundleGetPointerToData(aerosols, 'AMM_A_AIT', q, __RC__)
  qa(3,:,:,:) = q
  density(3)  = MAM_AMMONIUM_COMPONENT_DENSITY

  call ESMFL_BundleGetPointerToData(aerosols, 'SOA_A_AIT', q, __RC__)
  qa(4,:,:,:) = q
  density(4)  = MAM_SOA_COMPONENT_DENSITY

  call ESMFL_BundleGetPointerToData(aerosols, 'SS_A_AIT',  q, __RC__)
  qa(5,:,:,:) = q
  density(5)  = MAM_SEASALT_COMPONENT_DENSITY

  call ESMFL_BundleGetPointerToData(aerosols, 'DGN_WET_AIT', dgn_wet, __RC__)

  ext_ = 0.0 
  sca_ = 0.0 
  asy_ = 0.0

  call MAML_OpticsInterpolate(MAM7_MieTable(1), band, qa, density, dgn_wet, dp, ext_, sca_, asy_, nc, i1, i2, j1, j2, 1, km, rc)

  ext = ext + ext_
  sca = sca + sca_
  asy = asy + sca_*asy_
  deallocate(qa, density, __STAT__)


  ! compute ext, sca and ssa from aerosols in accumulation mode: su, amm, soa, pom, bc, ss
  nc = 6 + 1
  allocate(qa(nc,i1:i2,j1:j2,km), density(nc), __STAT__)

  qa      = 0.0
  density = 0.0

  call ESMFL_BundleGetPointerToData(aerosols, 'WTR_A_ACC', q, __RC__)
  qa(1,:,:,:) = q
  density(1)  = MAM_WATER_COMPONENT_DENSITY

  call ESMFL_BundleGetPointerToData(aerosols, 'SU_A_ACC',  q, __RC__)
  qa(2,:,:,:) = q
  density(2)  = MAM_SULFATE_COMPONENT_DENSITY

  call ESMFL_BundleGetPointerToData(aerosols, 'AMM_A_ACC', q, __RC__)
  qa(3,:,:,:) = q
  density(3)  = MAM_AMMONIUM_COMPONENT_DENSITY

  call ESMFL_BundleGetPointerToData(aerosols, 'SOA_A_ACC', q, __RC__)      
  qa(4,:,:,:) = q
  density(4)  = MAM_SOA_COMPONENT_DENSITY

  call ESMFL_BundleGetPointerToData(aerosols, 'POM_A_ACC', q, __RC__)
  qa(5,:,:,:) = q
  density(5)  = MAM_POM_COMPONENT_DENSITY

  call ESMFL_BundleGetPointerToData(aerosols, 'BC_A_ACC',  q, __RC__)
  qa(6,:,:,:) = q
  density(6)  = MAM_BLACK_CARBON_COMPONENT_DENSITY

  call ESMFL_BundleGetPointerToData(aerosols, 'SS_A_ACC',  q, __RC__)
  qa(7,:,:,:) = q
  density(7)  = MAM_SEASALT_COMPONENT_DENSITY

  call ESMFL_BundleGetPointerToData(aerosols, 'DGN_WET_ACC', dgn_wet, __RC__)

  ext_ = 0.0 
  sca_ = 0.0 
  asy_ = 0.0

  call MAML_OpticsInterpolate(MAM7_MieTable(2), band, qa, density, dgn_wet, dp, ext_, sca_, asy_, nc, i1, i2, j1, j2, 1, km, rc)

  ext = ext + ext_
  sca = sca + sca_
  asy = asy + sca_*asy_
  deallocate(qa, density, __STAT__)


  ! compute ext, sca and ssa from aerosols in primary carbon mode: pom, bc
  nc = 2 + 1
  allocate(qa(nc,i1:i2,j1:j2,km), density(nc), __STAT__)

  call ESMFL_BundleGetPointerToData(aerosols, 'WTR_A_PCM', q, __RC__)
  qa(1,:,:,:) = q
  density(1)  = MAM_WATER_COMPONENT_DENSITY

  call ESMFL_BundleGetPointerToData(aerosols, 'POM_A_PCM', q, __RC__)
  qa(2,:,:,:) = q
  density(2)  = MAM_POM_COMPONENT_DENSITY

  call ESMFL_BundleGetPointerToData(aerosols, 'BC_A_PCM',  q, __RC__)
  qa(3,:,:,:) = q
  density(3)  = MAM_BLACK_CARBON_COMPONENT_DENSITY

  call ESMFL_BundleGetPointerToData(aerosols, 'DGN_WET_PCM', dgn_wet, __RC__)

  ext_ = 0.0 
  sca_ = 0.0 
  asy_ = 0.0

  call MAML_OpticsInterpolate(MAM7_MieTable(3), band, qa, density, dgn_wet, dp, ext_, sca_, asy_, nc, i1, i2, j1, j2, 1, km, rc)

  ext = ext + ext_
  sca = sca + sca_
  asy = asy + sca_*asy_
  deallocate(qa, density, __STAT__)


  ! compute ext, sca and ssa from aerosols in fine seasalt mode: su, amm, ss
  nc = 3 + 1
  allocate(qa(nc,i1:i2,j1:j2,km), density(nc), __STAT__)

  call ESMFL_BundleGetPointerToData(aerosols, 'WTR_A_FSS', q, __RC__)
  qa(1,:,:,:) = q
  density(1)  = MAM_WATER_COMPONENT_DENSITY

  call ESMFL_BundleGetPointerToData(aerosols, 'SU_A_FSS',  q, __RC__)
  qa(2,:,:,:) = q
  density(2)  = MAM_SULFATE_COMPONENT_DENSITY

  call ESMFL_BundleGetPointerToData(aerosols, 'AMM_A_FSS', q, __RC__)
  qa(3,:,:,:) = q
  density(3)  = MAM_AMMONIUM_COMPONENT_DENSITY

  call ESMFL_BundleGetPointerToData(aerosols, 'SS_A_FSS',  q, __RC__)
  qa(4,:,:,:) = q
  density(4)  = MAM_SEASALT_COMPONENT_DENSITY

  call ESMFL_BundleGetPointerToData(aerosols, 'DGN_WET_FSS', dgn_wet, __RC__)

  ext_ = 0.0 
  sca_ = 0.0 
  asy_ = 0.0

  call MAML_OpticsInterpolate(MAM7_MieTable(4), band, qa, density, dgn_wet, dp, ext_, sca_, asy_, nc, i1, i2, j1, j2, 1, km, rc)

  ext = ext + ext_
  sca = sca + sca_
  asy = asy + sca_*asy_
  deallocate(qa, density, __STAT__)


  ! compute ext, sca and ssa from aerosols in coarse seasalt mode: su, amm, ss; lut = ('water', 'su', 'amm', 'ss')
  nc = 3 + 1
  allocate(qa(nc,i1:i2,j1:j2,km), density(nc), __STAT__)

  call ESMFL_BundleGetPointerToData(aerosols, 'WTR_A_CSS', q, __RC__)
  qa(1,:,:,:) = q
  density(1)  = MAM_WATER_COMPONENT_DENSITY

  call ESMFL_BundleGetPointerToData(aerosols, 'SU_A_CSS',  q, __RC__)
  qa(2,:,:,:) = q
  density(2)  = MAM_SULFATE_COMPONENT_DENSITY

  call ESMFL_BundleGetPointerToData(aerosols, 'AMM_A_CSS', q, __RC__)
  qa(3,:,:,:) = q
  density(3)  = MAM_AMMONIUM_COMPONENT_DENSITY

  call ESMFL_BundleGetPointerToData(aerosols, 'SS_A_CSS',  q, __RC__)
  qa(4,:,:,:) = q
  density(4)  = MAM_SEASALT_COMPONENT_DENSITY

  call ESMFL_BundleGetPointerToData(aerosols, 'DGN_WET_CSS', dgn_wet, __RC__)

  ext_ = 0.0 
  sca_ = 0.0 
  asy_ = 0.0

  call MAML_OpticsInterpolate(MAM7_MieTable(5), band, qa, density, dgn_wet, dp, ext_, sca_, asy_, nc, i1, i2, j1, j2, 1, km, rc)

  ext = ext + ext_
  sca = sca + sca_
  asy = asy + sca_*asy_
  deallocate(qa, density, __STAT__)


  ! compute ext, sca and ssa from aerosols in fine dust mode: su, amm, du; lut = ('water', 'su', 'amm', 'du')
  nc = 3 + 1
  allocate(qa(nc,i1:i2,j1:j2,km), density(nc), __STAT__)
 
  call ESMFL_BundleGetPointerToData(aerosols, 'WTR_A_FDU', q, __RC__)
  qa(1,:,:,:) = q
  density(1)  = MAM_WATER_COMPONENT_DENSITY

  call ESMFL_BundleGetPointerToData(aerosols, 'SU_A_FDU',  q, __RC__)
  qa(2,:,:,:) = q
  density(2)  = MAM_SULFATE_COMPONENT_DENSITY

  call ESMFL_BundleGetPointerToData(aerosols, 'AMM_A_FDU', q, __RC__)
  qa(3,:,:,:) = q
  density(3)  = MAM_AMMONIUM_COMPONENT_DENSITY

  call ESMFL_BundleGetPointerToData(aerosols, 'DU_A_FDU',  q, __RC__)
  qa(4,:,:,:) = q
  density(4)  = MAM_DUST_COMPONENT_DENSITY

  call ESMFL_BundleGetPointerToData(aerosols, 'DGN_WET_FDU', dgn_wet, __RC__)

  ext_ = 0.0 
  sca_ = 0.0 
  asy_ = 0.0
 
  call MAML_OpticsInterpolate(MAM7_MieTable(6), band, qa, density, dgn_wet, dp, ext_, sca_, asy_, nc, i1, i2, j1, j2, 1, km, rc)

  ext = ext + ext_
  sca = sca + sca_
  asy = asy + sca_*asy_
  deallocate(qa, density, __STAT__)


  ! compute ext, sca and ssa from aerosols in coarse dust mode: su, amm, du; lut = ('water', 'su', 'amm', 'du')
  nc = 3 + 1
  allocate(qa(nc,i1:i2,j1:j2,km), density(nc), __STAT__)

  call ESMFL_BundleGetPointerToData(aerosols, 'WTR_A_CDU', q, __RC__)
  qa(1,:,:,:) = q
  density(1)  = MAM_WATER_COMPONENT_DENSITY

  call ESMFL_BundleGetPointerToData(aerosols, 'SU_A_CDU',  q, __RC__)
  qa(2,:,:,:) = q
  density(2)  = MAM_SULFATE_COMPONENT_DENSITY

  call ESMFL_BundleGetPointerToData(aerosols, 'AMM_A_CDU', q, __RC__)
  qa(3,:,:,:) = q
  density(3)  = MAM_AMMONIUM_COMPONENT_DENSITY

  call ESMFL_BundleGetPointerToData(aerosols, 'DU_A_CDU',  q, __RC__)
  qa(4,:,:,:) = q
  density(4)  = MAM_DUST_COMPONENT_DENSITY

  call ESMFL_BundleGetPointerToData(aerosols, 'DGN_WET_CDU', dgn_wet, __RC__)

  ext_ = 0.0 
  sca_ = 0.0 
  asy_ = 0.0

  call MAML_OpticsInterpolate(MAM7_MieTable(7), band, qa, density, dgn_wet, dp, ext_, sca_, asy_, nc, i1, i2, j1, j2, 1, km, rc)

  ext = ext + ext_
  sca = sca + sca_
  asy = asy + sca_*asy_
  deallocate(qa, density, __STAT__)

 
  ! inputs for radiation:  ext, sca and sca*asy
  ! in the callback     : 'asy' = product of asy and sca
  !                       'ssa' = sca
  ssa = sca


  call ESMF_AttributeGet(state, name='extinction_in_air_due_to_ambient_aerosol', value=field_name, __RC__)
  if (field_name /= '') then 
      call MAPL_GetPointer(state, var, trim(field_name), __RC__)
      var = ext(:,:,:)
  end if

  call ESMF_AttributeGet(state, name='single_scattering_albedo_of_ambient_aerosol', value=field_name, __RC__)
  if (field_name /= '') then 
      call MAPL_GetPointer(state, var, trim(field_name), __RC__)
      var = ssa(:,:,:)
  end if

  call ESMF_AttributeGet(state, name='asymmetry_parameter_of_ambient_aerosol', value=field_name, __RC__)
  if (field_name /= '') then 
      call MAPL_GetPointer(state, var, trim(field_name), __RC__)
      var = asy(:,:,:)
  end if


  deallocate(dp,                    __STAT__) 
  deallocate(ext,  sca,  ssa,  asy, __STAT__)
  deallocate(ext_, sca_, asy_,      __STAT__)

  _RETURN(ESMF_SUCCESS)

 end subroutine aerosol_optics


subroutine aerosol_activation_properties(state, rc)

  implicit none

! Arguments
! ---------
  type(ESMF_State)     :: state
  integer, intent(out) :: rc


! Local
! ---------
  character(len=ESMF_MAXSTR)      :: mode              ! mode 
  type(ESMF_FieldBundle)          :: aerosols          ! field bundle containing the aerosol mass mixing ratios

  real, dimension(:,:,:), pointer :: q                 ! aerosol number or mass mixing ratio

  real, dimension(:,:,:), pointer :: num               ! number concentration of aerosol particles 
  real, dimension(:,:,:), pointer :: diameter          ! dry size of aerosol
  real, dimension(:,:,:), pointer :: sigma             ! width of aerosol mode
  real, dimension(:,:,:), pointer :: density           ! density of aerosol
  real, dimension(:,:,:), pointer :: hygroscopicity    ! hygroscopicity of aerosol 
  real, dimension(:,:,:), pointer :: f_dust            ! fraction of dust aerosol
  real, dimension(:,:,:), pointer :: f_soot            ! fraction of soot aerosol 
  real, dimension(:,:,:), pointer :: f_organic         ! fraction of organic aerosol

  real, allocatable, dimension(:,:,:,:) :: qa          ! temporary buffers
  real, allocatable, dimension(:) :: qa_density
  real, allocatable, dimension(:) :: qa_hygroscopicity
  real, allocatable, dimension(:) :: qa_f_dust
  real, allocatable, dimension(:) :: qa_f_soot
  real, allocatable, dimension(:) :: qa_f_organic
  

  character(len=ESMF_MAXSTR)      :: fld_name

  integer                         :: i1, j1, i2, j2, km
  integer                         :: i, j, k

  real                            :: sigma_
  integer                         :: ns

  integer                         :: STATUS
  character(len=ESMF_MAXSTR)      :: Iam


  Iam = 'MAM::aerosol_activation_properties()'


! Aerosol mode
! ------------
  call ESMF_AttributeGet(state, name='aerosol_mode', value=mode, __RC__)


! Aerosol mass mixing ratio and activation properties
! -------------------------
  call ESMF_StateGet(state, 'AEROSOLS', aerosols, __RC__)

! Activation activation properties
  call ESMF_AttributeGet(state, name='aerosol_number_concentration', value=fld_name, __RC__)
  call MAPL_GetPointer(state, num, trim(fld_name), __RC__)

  call ESMF_AttributeGet(state, name='aerosol_dry_size', value=fld_name, __RC__)
  call MAPL_GetPointer(state, diameter, trim(fld_name), __RC__)

  call ESMF_AttributeGet(state, name='width_of_aerosol_mode', value=fld_name, __RC__)
  call MAPL_GetPointer(state, sigma, trim(fld_name), __RC__) 
  
  call ESMF_AttributeGet(state, name='aerosol_density', value=fld_name, __RC__)
  call MAPL_GetPointer(state, density, trim(fld_name), __RC__)

  call ESMF_AttributeGet(state, name='aerosol_hygroscopicity', value=fld_name, __RC__)
  call MAPL_GetPointer(state, hygroscopicity, trim(fld_name), __RC__)

  call ESMF_AttributeGet(state, name='fraction_of_dust_aerosol', value=fld_name, __RC__)
  call MAPL_GetPointer(state, f_dust, trim(fld_name), __RC__)

  call ESMF_AttributeGet(state, name='fraction_of_soot_aerosol', value=fld_name, __RC__)
  call MAPL_GetPointer(state, f_soot, trim(fld_name), __RC__)

  call ESMF_AttributeGet(state, name='fraction_of_organic_aerosol', value=fld_name, __RC__)
  call MAPL_GetPointer(state, f_organic, trim(fld_name), __RC__)

  


! Obtain aerosol activation properties of this aerosol mode
! ---------------------------------------------------------
  i1 = lbound(num, 1); i2 = ubound(num, 1)
  j1 = lbound(num, 2); j2 = ubound(num, 2)
                       km = ubound(num, 3)

  call ESMF_StateGet(state, 'AEROSOLS', aerosols, __RC__)

  call ESMFL_BundleGetPointerToData(aerosols, 'NUM_A_'//trim(mode), q, __RC__)
  num = q

  select case(mode)
  case (MAM7_AITKEN_MODE_NAME)
      sigma_ = MAM7_AITKEN_MODE_SIGMA

      ns = 4

      allocate(qa(ns,i1:i2,j1:j2,km), __STAT__)
      allocate(qa_density(ns), qa_hygroscopicity(ns), __STAT__)
      allocate(qa_f_dust(ns), qa_f_soot(ns), qa_f_organic(ns), __STAT__)
 
      qa = 0.0
      qa_density = 0.0
      qa_hygroscopicity = 0.0
      qa_f_dust    = 0.0 
      qa_f_soot    = 0.0
      qa_f_organic = 0.0

      call ESMFL_BundleGetPointerToData(aerosols, 'SU_A_AIT',  q, __RC__)
      qa(1,:,:,:) = q(:,:,:)
      qa_density(1) = MAM_SULFATE_COMPONENT_DENSITY
      qa_hygroscopicity(1) = MAM_SULFATE_COMPONENT_HYGROSCOPICITY

      call ESMFL_BundleGetPointerToData(aerosols, 'AMM_A_AIT', q, __RC__)
      qa(2,:,:,:) = q(:,:,:)
      qa_density(2) = MAM_AMMONIUM_COMPONENT_DENSITY
      qa_hygroscopicity(2) = MAM_AMMONIUM_COMPONENT_HYGROSCOPICITY

      call ESMFL_BundleGetPointerToData(aerosols, 'SOA_A_AIT', q, __RC__)
      qa(3,:,:,:) = q(:,:,:)
      qa_density(3) = MAM_SOA_COMPONENT_DENSITY
      qa_hygroscopicity(3) = MAM_SOA_COMPONENT_HYGROSCOPICITY
      qa_f_organic(3) = 1.0


      call ESMFL_BundleGetPointerToData(aerosols, 'SS_A_AIT',  q, __RC__)
      qa(4,:,:,:) = q(:,:,:)
      qa_density(4) = MAM_SEASALT_COMPONENT_DENSITY
      qa_hygroscopicity(4) = MAM_SEASALT_COMPONENT_HYGROSCOPICITY

      call aap_(diameter, density, hygroscopicity, f_dust, f_soot, f_organic, &
                num, sigma_, &
                qa, qa_density, qa_hygroscopicity, qa_f_dust, qa_f_soot, qa_f_organic, ns, &
                i1, i2, j1, j2, km, rc)

      sigma = log(sigma_)

      deallocate(qa, __STAT__)
      deallocate(qa_density, qa_hygroscopicity, __STAT__)
      deallocate(qa_f_dust, qa_f_soot, qa_f_organic, __STAT__)
      
  case (MAM7_ACCUMULATION_MODE_NAME)
      sigma_ = MAM7_ACCUMULATION_MODE_SIGMA

      ns = 6

      allocate(qa(ns,i1:i2,j1:j2,km), __STAT__)
      allocate(qa_density(ns), qa_hygroscopicity(ns), __STAT__)
      allocate(qa_f_dust(ns), qa_f_soot(ns), qa_f_organic(ns), __STAT__)
 
      qa = 0.0
      qa_density = 0.0
      qa_hygroscopicity = 0.0
      qa_f_dust    = 0.0 
      qa_f_soot    = 0.0
      qa_f_organic = 0.0

      call ESMFL_BundleGetPointerToData(aerosols, 'SU_A_ACC',  q, __RC__)
      qa(1,:,:,:) = q(:,:,:)
      qa_density(1) = MAM_SULFATE_COMPONENT_DENSITY
      qa_hygroscopicity(1) = MAM_SULFATE_COMPONENT_HYGROSCOPICITY

      call ESMFL_BundleGetPointerToData(aerosols, 'AMM_A_ACC', q, __RC__)
      qa(2,:,:,:) = q(:,:,:)
      qa_density(2) = MAM_AMMONIUM_COMPONENT_DENSITY
      qa_hygroscopicity(2) = MAM_AMMONIUM_COMPONENT_HYGROSCOPICITY

      call ESMFL_BundleGetPointerToData(aerosols, 'SOA_A_AIT', q, __RC__)
      qa(3,:,:,:) = q(:,:,:)
      qa_density(3) = MAM_SOA_COMPONENT_DENSITY
      qa_hygroscopicity(3) = MAM_SOA_COMPONENT_HYGROSCOPICITY
      qa_f_organic(3) = 1.0

      call ESMFL_BundleGetPointerToData(aerosols, 'POM_A_ACC', q, __RC__)
      qa(4,:,:,:) = q(:,:,:)
      qa_density(4) = MAM_POM_COMPONENT_DENSITY
      qa_hygroscopicity(4) = MAM_POM_COMPONENT_HYGROSCOPICITY
      qa_f_organic(4) = 1.0

      call ESMFL_BundleGetPointerToData(aerosols, 'BC_A_ACC',  q, __RC__)
      qa(5,:,:,:) = q(:,:,:)
      qa_density(5) = MAM_BLACK_CARBON_COMPONENT_DENSITY
      qa_hygroscopicity(5) = MAM_BLACK_CARBON_COMPONENT_HYGROSCOPICITY
      qa_f_soot(5) = 1.0

      call ESMFL_BundleGetPointerToData(aerosols, 'SS_A_AIT',  q, __RC__)
      qa(6,:,:,:) = q(:,:,:)
      qa_density(6) = MAM_SEASALT_COMPONENT_DENSITY
      qa_hygroscopicity(6) = MAM_SEASALT_COMPONENT_HYGROSCOPICITY

      call aap_(diameter, density, hygroscopicity, f_dust, f_soot, f_organic, &
                num, sigma_, &
                qa, qa_density, qa_hygroscopicity, qa_f_dust, qa_f_soot, qa_f_organic, ns, &
                i1, i2, j1, j2, km, rc)

      sigma = log(sigma_)

      deallocate(qa, __STAT__)
      deallocate(qa_density, qa_hygroscopicity, __STAT__)
      deallocate(qa_f_dust, qa_f_soot, qa_f_organic, __STAT__)


  case (MAM7_PRIMARY_CARBON_MODE_NAME)
      sigma_ = MAM7_PRIMARY_CARBON_MODE_SIGMA

      ns = 2

      allocate(qa(ns,i1:i2,j1:j2,km), __STAT__)
      allocate(qa_density(ns), qa_hygroscopicity(ns), __STAT__)
      allocate(qa_f_dust(ns), qa_f_soot(ns), qa_f_organic(ns), __STAT__)
 
      qa = 0.0
      qa_density = 0.0
      qa_hygroscopicity = 0.0
      qa_f_dust    = 0.0 
      qa_f_soot    = 0.0
      qa_f_organic = 0.0

      call ESMFL_BundleGetPointerToData(aerosols, 'POM_A_PCM', q, __RC__)
      qa(1,:,:,:) = q(:,:,:)
      qa_density(1) = MAM_POM_COMPONENT_DENSITY
      qa_hygroscopicity(1) = MAM_POM_COMPONENT_HYGROSCOPICITY
      qa_f_organic(1) = 1.0


      call ESMFL_BundleGetPointerToData(aerosols, 'BC_A_PCM',  q, __RC__)
      qa(2,:,:,:) = q(:,:,:)
      qa_density(2) = MAM_BLACK_CARBON_COMPONENT_DENSITY
      qa_hygroscopicity(2) = MAM_BLACK_CARBON_COMPONENT_HYGROSCOPICITY
      qa_f_soot(2) = 1.0

      call aap_(diameter, density, hygroscopicity, f_dust, f_soot, f_organic, &
                num, sigma_, &
                qa, qa_density, qa_hygroscopicity, qa_f_dust, qa_f_soot, qa_f_organic, ns, &
                i1, i2, j1, j2, km, rc)

      sigma = log(sigma_)

      deallocate(qa, __STAT__)
      deallocate(qa_density, qa_hygroscopicity, __STAT__)
      deallocate(qa_f_dust, qa_f_soot, qa_f_organic, __STAT__)


  case (MAM7_FINE_SEASALT_MODE_NAME)
      sigma_ = MAM7_FINE_SEASALT_MODE_SIGMA

      ns = 3

      allocate(qa(ns,i1:i2,j1:j2,km), __STAT__)
      allocate(qa_density(ns), qa_hygroscopicity(ns), __STAT__)
      allocate(qa_f_dust(ns), qa_f_soot(ns), qa_f_organic(ns), __STAT__)
 
      qa = 0.0
      qa_density = 0.0
      qa_hygroscopicity = 0.0
      qa_f_dust    = 0.0 
      qa_f_soot    = 0.0
      qa_f_organic = 0.0

      call ESMFL_BundleGetPointerToData(aerosols, 'SU_A_FSS',  q, __RC__)
      qa(1,:,:,:) = q(:,:,:)
      qa_density(1) = MAM_SULFATE_COMPONENT_DENSITY
      qa_hygroscopicity(1) = MAM_SULFATE_COMPONENT_HYGROSCOPICITY

      call ESMFL_BundleGetPointerToData(aerosols, 'AMM_A_FSS', q, __RC__)
      qa(2,:,:,:) = q(:,:,:)
      qa_density(2) = MAM_AMMONIUM_COMPONENT_DENSITY
      qa_hygroscopicity(2) = MAM_AMMONIUM_COMPONENT_HYGROSCOPICITY

      call ESMFL_BundleGetPointerToData(aerosols, 'SS_A_FSS',  q, __RC__)
      qa(3,:,:,:) = q(:,:,:)
      qa_density(3) = MAM_SEASALT_COMPONENT_DENSITY
      qa_hygroscopicity(3) = MAM_SEASALT_COMPONENT_HYGROSCOPICITY

      call aap_(diameter, density, hygroscopicity, f_dust, f_soot, f_organic, &
                num, sigma_, &
                qa, qa_density, qa_hygroscopicity, qa_f_dust, qa_f_soot, qa_f_organic, ns, &
                i1, i2, j1, j2, km, rc)

      sigma = log(sigma_)

      deallocate(qa, __STAT__)
      deallocate(qa_density, qa_hygroscopicity, __STAT__)
      deallocate(qa_f_dust, qa_f_soot, qa_f_organic, __STAT__)


  case (MAM7_FINE_DUST_MODE_NAME)
      sigma_ = MAM7_FINE_DUST_MODE_SIGMA

      ns = 3

      allocate(qa(ns,i1:i2,j1:j2,km), __STAT__)
      allocate(qa_density(ns), qa_hygroscopicity(ns), __STAT__)
      allocate(qa_f_dust(ns), qa_f_soot(ns), qa_f_organic(ns), __STAT__)
 
      qa = 0.0
      qa_density = 0.0
      qa_hygroscopicity = 0.0
      qa_f_dust    = 0.0 
      qa_f_soot    = 0.0
      qa_f_organic = 0.0

      call ESMFL_BundleGetPointerToData(aerosols, 'SU_A_FDU',  q, __RC__)
      qa(1,:,:,:) = q(:,:,:)
      qa_density(1) = MAM_SULFATE_COMPONENT_DENSITY
      qa_hygroscopicity(1) = MAM_SULFATE_COMPONENT_HYGROSCOPICITY

      call ESMFL_BundleGetPointerToData(aerosols, 'AMM_A_FDU', q, __RC__)
      qa(2,:,:,:) = q(:,:,:)
      qa_density(2) = MAM_AMMONIUM_COMPONENT_DENSITY
      qa_hygroscopicity(2) = MAM_AMMONIUM_COMPONENT_HYGROSCOPICITY

      call ESMFL_BundleGetPointerToData(aerosols, 'DU_A_FDU',  q, __RC__)
      qa(3,:,:,:) = q(:,:,:)
      qa_density(3) = MAM_DUST_COMPONENT_DENSITY
      qa_hygroscopicity(3) = MAM_DUST_COMPONENT_HYGROSCOPICITY
      qa_f_dust(3) = 1.0

      call aap_(diameter, density, hygroscopicity, f_dust, f_soot, f_organic, &
                num, sigma_, &
                qa, qa_density, qa_hygroscopicity, qa_f_dust, qa_f_soot, qa_f_organic, ns, &
                i1, i2, j1, j2, km, rc)

      sigma = log(sigma_)

      deallocate(qa, __STAT__)
      deallocate(qa_density, qa_hygroscopicity, __STAT__)
      deallocate(qa_f_dust, qa_f_soot, qa_f_organic, __STAT__)


  case (MAM7_COARSE_SEASALT_MODE_NAME)
      sigma_ = MAM7_COARSE_SEASALT_MODE_SIGMA

      ns = 3

      allocate(qa(ns,i1:i2,j1:j2,km), __STAT__)
      allocate(qa_density(ns), qa_hygroscopicity(ns), __STAT__)
      allocate(qa_f_dust(ns), qa_f_soot(ns), qa_f_organic(ns), __STAT__)
 
      qa = 0.0
      qa_density = 0.0
      qa_hygroscopicity = 0.0
      qa_f_dust    = 0.0 
      qa_f_soot    = 0.0
      qa_f_organic = 0.0

      call ESMFL_BundleGetPointerToData(aerosols, 'SU_A_CSS',  q, __RC__)
      qa(1,:,:,:) = q(:,:,:)
      qa_density(1) = MAM_SULFATE_COMPONENT_DENSITY
      qa_hygroscopicity(1) = MAM_SULFATE_COMPONENT_HYGROSCOPICITY

      call ESMFL_BundleGetPointerToData(aerosols, 'AMM_A_CSS', q, __RC__)
      qa(2,:,:,:) = q(:,:,:)
      qa_density(2) = MAM_AMMONIUM_COMPONENT_DENSITY
      qa_hygroscopicity(2) = MAM_AMMONIUM_COMPONENT_HYGROSCOPICITY

      call ESMFL_BundleGetPointerToData(aerosols, 'SS_A_CSS',  q, __RC__)
      qa(3,:,:,:) = q(:,:,:)
      qa_density(3) = MAM_SEASALT_COMPONENT_DENSITY
      qa_hygroscopicity(3) = MAM_SEASALT_COMPONENT_HYGROSCOPICITY

      call aap_(diameter, density, hygroscopicity, f_dust, f_soot, f_organic, &
                num, sigma_, &
                qa, qa_density, qa_hygroscopicity, qa_f_dust, qa_f_soot, qa_f_organic, ns, &
                i1, i2, j1, j2, km, rc)

      sigma = log(sigma_)

      deallocate(qa, __STAT__)
      deallocate(qa_density, qa_hygroscopicity, __STAT__)
      deallocate(qa_f_dust, qa_f_soot, qa_f_organic, __STAT__)


  case (MAM7_COARSE_DUST_MODE_NAME)
      sigma_ = MAM7_COARSE_DUST_MODE_SIGMA

      ns = 3

      allocate(qa(ns,i1:i2,j1:j2,km), __STAT__)
      allocate(qa_density(ns), qa_hygroscopicity(ns), __STAT__)
      allocate(qa_f_dust(ns), qa_f_soot(ns), qa_f_organic(ns), __STAT__)
 
      qa = 0.0
      qa_density = 0.0
      qa_hygroscopicity = 0.0
      qa_f_dust    = 0.0 
      qa_f_soot    = 0.0
      qa_f_organic = 0.0

      call ESMFL_BundleGetPointerToData(aerosols, 'SU_A_CDU',  q, __RC__)
      qa(1,:,:,:) = q(:,:,:)
      qa_density(1) = MAM_SULFATE_COMPONENT_DENSITY
      qa_hygroscopicity(1) = MAM_SULFATE_COMPONENT_HYGROSCOPICITY

      call ESMFL_BundleGetPointerToData(aerosols, 'AMM_A_CDU', q, __RC__)
      qa(2,:,:,:) = q(:,:,:)
      qa_density(2) = MAM_AMMONIUM_COMPONENT_DENSITY
      qa_hygroscopicity(2) = MAM_AMMONIUM_COMPONENT_HYGROSCOPICITY

      call ESMFL_BundleGetPointerToData(aerosols, 'DU_A_CDU',  q, __RC__)
      qa(3,:,:,:) = q(:,:,:)
      qa_density(3) = MAM_DUST_COMPONENT_DENSITY
      qa_hygroscopicity(3) = MAM_DUST_COMPONENT_HYGROSCOPICITY
      qa_f_dust(3) = 1.0

      call aap_(diameter, density, hygroscopicity, f_dust, f_soot, f_organic, &
                num, sigma_, &
                qa, qa_density, qa_hygroscopicity, qa_f_dust, qa_f_soot, qa_f_organic, ns, &
                i1, i2, j1, j2, km, rc)

      sigma = log(sigma_)

      deallocate(qa, __STAT__)
      deallocate(qa_density, qa_hygroscopicity, __STAT__)
      deallocate(qa_f_dust, qa_f_soot, qa_f_organic, __STAT__)
  
  case default
      __raise__(MAM_UNKNOWN_AEROSOL_MODE_ERROR,"Unknown aerosol mode used in the MAM aerosol activation properties method: "//trim(mode))

  end select
  

  _RETURN(ESMF_SUCCESS)

contains

    subroutine aap_(diameter, density, hygroscopicity, &
                    f_dust, f_soot, f_organic, &
                    q_num, sigma, &
                    q, q_density, q_hygroscopicity, &
                    q_f_dust, q_f_soot, q_f_organic, & 
                    ns, &
                    i1, i2, j1, j2, km, &
                    rc)
     
     implicit none

     integer, intent(in) :: i1, i2                                    ! dimension bounds
     integer, intent(in) :: j1, j2                                    ! ... // ..
     integer, intent(in) :: km                                        ! ... // ..

     integer, intent(in) :: ns                                        ! number of species

     real, intent(in ), dimension(i1:i2,j1:j2,km)   :: q_num          ! number mixing ratio, #-particles kg-1
     real, intent(in )                              :: sigma          ! width of the mode

     real, intent(in ), dimension(ns,i1:i2,j1:j2,km):: q              ! aerosol mass mixing ratio, kg kg-1
     real, intent(in ), dimension(ns)               :: q_density      ! density of species 
     real, intent(in ), dimension(ns)               :: q_hygroscopicity
     real, intent(in ), dimension(ns)               :: q_f_dust       ! 
     real, intent(in ), dimension(ns)               :: q_f_soot       !
     real, intent(in ), dimension(ns)               :: q_f_organic    ! 

     real, intent(out), dimension(i1:i2,j1:j2,km)   :: diameter       ! dry size of aerosol
     real, intent(out), dimension(i1:i2,j1:j2,km)   :: density        ! density of aerosol
     real, intent(out), dimension(i1:i2,j1:j2,km)   :: hygroscopicity ! hygroscopicity of aerosol 
     real, intent(out), dimension(i1:i2,j1:j2,km)   :: f_dust         ! fraction of dust aerosol
     real, intent(out), dimension(i1:i2,j1:j2,km)   :: f_soot         ! fraction of soot aerosol 
     real, intent(out), dimension(i1:i2,j1:j2,km)   :: f_organic      ! fraction of organic aerosol

     integer, intent(out) :: rc                                       ! return code

     ! local
     real, dimension(ns) :: q_
     real, dimension(ns) :: v_

     real :: mass, vol
     real :: f

     integer :: i, j, k
     
     integer :: STATUS
     character(len=ESMF_MAXSTR) :: Iam = 'MAM::aerosol_activation_properties::aap_()'

     ! vol = number * (MAPL_PI/6) * Dgn**3 * exp(4.5 * log(sigma)**2)
     f = 1.0 / ((MAPL_PI/6) * exp(4.5 * log(sigma)*log(sigma)))
      
     do k = 1, km
         do j = j1, j2
             do i = i1, i2

                 q_ = q(:,i,j,k)
                 v_ = q_ / q_density

                 mass = sum(q_)
                 vol  = sum(v_)

                 if ((vol > 0) .and. (num(i,j,k) > 0)) then
                     f_dust(i,j,k)    = sum(v_ * q_f_dust)    / vol
                     f_soot(i,j,k)    = sum(v_ * q_f_soot)    / vol
                     f_organic(i,j,k) = sum(v_ * q_f_organic) / vol

                     hygroscopicity(i,j,k) = sum(v_ * q_hygroscopicity) / vol

                     density(i,j,k)   = mass / vol

                     ! num = (q * air_density) / ((MAPL_PI/6) * density * Dgn**3 * exp(4.5 * log(sigma)*log(sigma)))
                     diameter(i,j,k)  = ((vol / num(i,j,k)) * f)**(1.0/3.0)
                 else
!!!                  print *, 'DEBUG::MAM::aap_() ', num(i,j,k), vol

                     f_dust(i,j,k)    = 0.0
                     f_soot(i,j,k)    = 0.0
                     f_organic(i,j,k) = 0.0

                     density(i,j,k)   = q_density(1)
                     hygroscopicity(i,j,k) = q_hygroscopicity(1)
                     diameter(i,j,k)  = 0.0
                 end if

             end do
         end do
     end do

     _RETURN(ESMF_SUCCESS)

    end subroutine aap_
    
 end subroutine aerosol_activation_properties



 function constituent_index_(constituent_name, rc) result (i)
  
  use constituents,      only: pcnst, cnst_name
  
  character(len=*),  intent(in)  :: constituent_name
  integer, optional, intent(out) :: rc

  ! local
  integer                    :: n, i
  character(len=ESMF_MAXSTR) :: Iam
  integer                    :: status


  Iam = 'MAM::constituent_index_()'

  i = 0

  do n = 1, pcnst
      if (constituent_name == cnst_name(n)) then
          i = n
          exit
      end if
  end do

  if (i == 0) then
      __raise__ (MAM_UNKNOWN_AEROSOL_CONSTITUENT_ERROR, "MAM::Unknown constituent: " // trim(constituent_name))
  end if

  _RETURN(ESMF_SUCCESS)
 end function constituent_index_


 end module MAMchem_GridCompMod
