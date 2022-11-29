#include "MAPL_Generic.h"

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 910.1     !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: GAAS_GridCompMod - Implements GEOS-5 Aerosol Assimation 
!
! !INTERFACE:
!
   MODULE GAAS_GridCompMod
!
! !USES:
!
   Use ESMF
   Use MAPL
   Use m_StrTemplate  

   Use  LDE_Mod
   Use  Chem_SimpleBundleMod
   Use  Chem_RegistryMod
   Use  Chem_MieMod
   Use  Chem_AodMod
   Use  MAPL_GridManagerMod
   Use  MAPL_LatLonGridFactoryMod
   Use  MAPL_CubedSphereGridFactoryMod, only: CubedSphereGridFactory
   use mpi

   IMPLICIT NONE
   PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:

   PUBLIC SetServices
!
! !DESCRIPTION: 
!
!  {\tt GAAS\_GridComp} is an ESMF gridded component implementing
!  the GEOS-5 Aerosol Analysis System (GAAS). 
!
!  Developed for GEOS-5 release Fortuna 2.3 and later.
!
! !REVISION HISTORY:
!
!  30Nov2010 da Silva  Initial version.
!  26Mar2021 E.Sherman Revised to use ESMF/MAPL
!
!EOP
!-------------------------------------------------------------------------

! Legacy state
! ------------
  TYPE GAAS_State
     PRIVATE

     type (ESMF_Config)         :: CF           ! Private Config
     type (ESMF_Grid)           :: aodGrid      ! AOD Grid (Vertical is "channels")
     type (Chem_Mie)            :: Mie          ! Mie Tables, etc

     type (Chem_Registry)       :: aerReg       ! Registry with aerosol tracers
     type (Chem_Registry)       :: aodReg       ! Registry with single AOD tracer

     type(LDE)                  :: E            ! LDE object

     type (MAPL_SimpleBundle)   :: q_f          ! Concentration background 
     type (MAPL_SimpleBundle)   :: q_a          ! Concentration analysis

     type (MAPL_SimpleBundle)   :: y_f          ! On-line  AOD background
     type (MAPL_SimpleBundle)   :: z_f          ! Off-line AOD background
     type (MAPL_SimpleBundle)   :: z_a          ! off-line AOD analysis
     type (MAPL_SimpleBundle)   :: z_k          ! Averaging kernel approximation

     type (MAPL_SimpleBundle)   :: y_a          ! Background adjusted AOD analysis
     type (MAPL_SimpleBundle)   :: y_d          ! AOD Analysis increments

     logical                    :: verbose=.FALSE.

     real                       :: eps=0.01     ! parameter for Log(AOD+eps) transform
     real                       :: channel=550. ! Single channel to analyze

     logical                    :: no_fuss=.FALSE. ! do not fuss if analysis file is missing

  END TYPE GAAS_State

! Hook for the ESMF
! -----------------
  TYPE GAAS_Wrap
     TYPE (GAAS_State), pointer :: PTR => null()
  END TYPE GAAS_WRAP

CONTAINS

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 910.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: SetServices --- Sets IRF services for the GAAS Grid Component
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
!  30Nov2010 da Silva  Initial version.
!  26Mar20201 E.Sherman Added AERO state to IMPORT
!EOP
!-------------------------------------------------------------------------

!   Local derived type aliases
!   --------------------------
    type (GAAS_State), pointer  :: self   ! internal, that is
    type (GAAS_wrap)            :: wrap

    character(len=ESMF_MAXSTR) :: comp_name

                            __Iam__('SetServices')

!                              ------------

!   Get my name and set-up traceback handle
!   ---------------------------------------
    call ESMF_GridCompGet( GC, name=comp_name, __RC__ )
    Iam = TRIM(comp_name) // '::' // TRIM(Iam)

!   Greetings
!   ---------
    IF(MAPL_AM_I_ROOT()) THEN
         PRINT *, TRIM(Iam)//': ACTIVE'
         PRINT *,' '
    END IF

!   Wrap internal state for storing in GC; rename legacyState
!   -------------------------------------
    allocate ( self, stat=STATUS )
    VERIFY_(STATUS)
    wrap%ptr => self
 
!   Load private Config Attributes
!   ------------------------------
    self%CF = ESMF_ConfigCreate(__RC__)
    call ESMF_ConfigLoadFile ( self%CF,'GAAS_GridComp.rc',__RC__)
    call ESMF_ConfigGetAttribute(self%CF, self%verbose, Label='verbose:',  __RC__ )
    call ESMF_ConfigGetAttribute(self%CF, self%eps, Label='eps_for_log_transform_aod:', &
                                          default=0.01, __RC__)
    call ESMF_ConfigGetAttribute(self%CF, self%channel, Label='single_channel:',  &
                                          default = 550., __RC__ )

!                       ------------------------
!                       ESMF Functional Services
!                       ------------------------

!   Set the Initialize, Run, Finalize entry points
!   ----------------------------------------------
    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_INITIALIZE,  Initialize_, __RC__ )
    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,         Run_,        __RC__ )
    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_FINALIZE,    Finalize_,   __RC__ )
        
!   Store internal state in GC
!   --------------------------
    call ESMF_UserCompSetInternalState ( GC, 'GAAS_state', wrap, STATUS )
    VERIFY_(STATUS)
  
!                         ------------------
!                         MAPL Data Services
!                         ------------------

!!BOP
!
! !IMPORT STATE:

# include "GAAS_ImportSpec___.h"

    call MAPL_AddImportSpec(GC,                                   &
       LONG_NAME  = 'aerosols',                                   &
       UNITS      = 'kg kg-1',                                    &
       SHORT_NAME = 'AERO',                                       &
       DIMS       = MAPL_DimsHorzVert,                            &
       VLOCATION  = MAPL_VLocationCenter,                         &
       DATATYPE   = MAPL_StateItem,                               &
       RESTART    = MAPL_RestartSkip, __RC__)

! !INTERNAL STATE: (none for now)

! !EXTERNAL STATE:

#   include "GAAS_ExportSpec___.h"

!   Generic Set Services
!   --------------------
    call MAPL_GenericSetServices ( GC, __RC__ )

!   All done
!   --------
    RETURN_(ESMF_SUCCESS)

  END SUBROUTINE SetServices


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Initialize_ --- Initialize GAAS
!
! !INTERFACE:
!

   SUBROUTINE Initialize_ ( GC, IMPORT, EXPORT, CLOCK, rc )

! !USES:

   implicit NONE

! !INPUT PARAMETERS:

   type(ESMF_Clock),  intent(inout) :: CLOCK     ! The clock

! !OUTPUT PARAMETERS:

   type(ESMF_GridComp), intent(inout)  :: GC     ! Grid Component
   type(ESMF_State), intent(inout) :: IMPORT     ! Import State
   type(ESMF_State), intent(inout) :: EXPORT     ! Export State
   integer, intent(out)            ::  rc        ! Error return code:
                                                 !  0 - all is well
                                                 !  1 - 

! !DESCRIPTION: This is a simple ESMF wrapper.
!
! !REVISION HISTORY:
!
!  30Nov2010 da Silva   Initial version.
!  12Feb2021 E. Sherman Removed GEOS-4 legacy constructs. Uses ESMF/MAPL constructs.
!EOP
!-------------------------------------------------------------------------

    type(GAAS_state), pointer     :: self     ! Legacy state
    type(ESMF_Grid)               :: Grid     ! Grid
    type(ESMF_Config)             :: CF       ! Universal Config 
    type(ESMF_Time)               :: Time     ! Current time

    type(ESMF_State)              :: aero ! Aersol state
    type(ESMF_FieldBundle)        :: aeroSerialBundle ! serialized aerosol bundle

    integer                       :: nymd, nhms  ! date, time

    integer                       :: dims(3), IM_World, JM_World, CM_World, Nx, Ny
    logical                       :: isCubed
    type(LatLonGridFactory)       :: ll_factory
    type(CubedSphereGridFactory)  :: cs_factory

    character(len=ESMF_MAXSTR)    :: comp_name

    __Iam__('Initialize_')

!  Get my name and set-up traceback handle
!  ---------------------------------------
   call ESMF_GridCompGet( GC, name=COMP_NAME, __RC__ )
   Iam = trim(comp_name) // ':: Initialize'

!                               --------
  
   if (MAPL_AM_I_ROOT()) then
      PRINT *, TRIM(Iam)//': Starting...'
      PRINT *,' '
   end if

!  Initialize MAPL Generic
!  -----------------------
   call MAPL_GenericInitialize ( gc, IMPORT, EXPORT, clock,  __RC__ )

!  Extract relevant runtime information
!  ------------------------------------
   call extract_ ( GC, CLOCK, self, GRID, CF, time, nymd, nhms, __RC__)

!  Registries
!  ----------
   self%aerReg = Chem_RegistryCreate ( rcfile='GAAS_AerRegistry.rc', __RC__ )  ! REMOVE
   self%aodReg = Chem_RegistryCreate ( rcfile='GAAS_AodRegistry.rc', __RC__ )  ! REMOVE

!  Mie tables, etc
!  ---------------
   self%Mie = Chem_MieCreate(self%CF, chemReg=self%aerReg, __RC__)   !REMOVE

!  Grid size, etc
!  --------------
   call MAPL_GridGet ( GRID, globalCellCountPerDim=DIMS, RC=STATUS)
   IM_World = dims(1)
   JM_World = dims(2)
   CM_World = self%Mie%nch ! Mie tables determine number of channels
   call ESMF_ConfigGetAttribute(CF, Nx, Label='NX:', __RC__)
   call ESMF_ConfigGetAttribute(CF, Ny, Label='NY:', __RC__)

!  Lat lon or cubed sphere?
!  ------------------------
   if ( JM_World == 6*IM_World ) then
       isCubed = .True.
   else
       isCubed = .False.
   end if

!  Create AOD grid
!  ---------------
   if ( isCubed ) then
       cs_factory = CubedSphereGridFactory(im_world=im_world,lm=cm_world,nx=nx,ny=ny/6,__RC__)
       self%aodGrid = grid_manager%make_grid(cs_factory,__RC__)

   else
       ll_factory = LatLonGridFactory(grid_name='aodGrid',nx=nx,ny=ny, & 
                    im_world=im_world,jm_world=jm_world,lm=cm_world, &
                    pole='PC',dateline='DC',__RC__)
       self%aodGrid = grid_manager%make_grid(ll_factory,__RC__)
   end if
   call ESMF_GridValidate(self%aodGrid,__RC__)

!  Prepare aerosol SimpleBundle
!  -----------------------------
!  Execute AERO's serialize_bundle method
   call ESMF_StateGet(IMPORT, 'AERO', aero, __RC__)
   call ESMF_MethodExecute(aero, label="serialize_bundle", __RC__)
   call ESMF_StateGet(aero, 'serialized_aerosolBundle', aeroSerialBundle, __RC__)

!  Create SimpleBundle from aeroBundle
!  Associate SimpleBundle with concentration analysis/background
   self%q_f = MAPL_SimpleBundleCreate(aeroSerialBundle, name='q_f', __RC__)
   self%q_a = MAPL_SimpleBundleCreate(aeroSerialBundle, name='q_a', __RC__)

!  Create AOD Simple Bundles
!  -------------------------
   self%y_f = MAPL_SimpleBundleCreate (self%aodGrid, rc=status, name='y_f')
   self%y_f%n2d = 1
   allocate(self%y_f%r2(1), __STAT__)
   self%y_f%r2(1)%name='aod_bkg'

   self%y_a = MAPL_SimpleBundleCreate (self%aodGrid, rc=status, name='y_a')
   self%y_a%n2d = 1
   allocate(self%y_a%r2(1), __STAT__)
   self%y_a%r2(1)%name='aod_ana'

   self%y_d = MAPL_SimpleBundleCreate (self%aodGrid, rc=status, name='y_d')
   self%y_d%n2d = 1
   allocate(self%y_d%r2(1), __STAT__)
   self%y_d%r2(1)%name='aod_inc'

!  Create AOD Simple Bundles from off-line analysis
!  ------------------------------------------------
   self%z_f = MAPL_SimpleBundleCreate (self%aodGrid, rc=status, name='z_f')
   self%z_f%n2d = 1
   allocate(self%z_f%r2(1), __STAT__)
   self%z_f%r2(1)%name='z_f'

   self%z_a = MAPL_SimpleBundleCreate (self%aodGrid, rc=status, name='z_a')
   self%z_a%n2d = 1
   allocate(self%z_a%r2(1), __STAT__)
   self%z_a%r2(1)%name='z_a'

   self%z_k = MAPL_SimpleBundleCreate (self%aodGrid, rc=status, name='z_k')
   self%z_k%n2d = 1
   allocate(self%z_k%r2(1), __STAT__)
   self%z_k%r2(1)%name='z_k'

!  Create LDE object
!  -----------------
   call LDE_Create ( self%E, self%CF, self%aodGrid, __RC__ )

#ifdef PRINT_STATES

   if (MAPL_AM_I_ROOT()) then
       print *, trim(Iam)//': IMPORT   State during Initialize():'
       call ESMF_StatePrint ( IMPORT )
       print *, trim(Iam)//': EXPORT   State during Initialize():' 
       call ESMF_StatePrint ( EXPORT )
    end if

#endif


!  All done
!  --------
   RETURN_(ESMF_SUCCESS)

   END SUBROUTINE Initialize_


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Run_ --- Runs GAAS
!
! !INTERFACE:
!

   SUBROUTINE Run_ ( gc, IMPORT, EXPORT, CLOCK, rc )

! !USES:

  implicit NONE

! !INPUT PARAMETERS:

   type(ESMF_Clock),  intent(inout) :: CLOCK     ! The clock

! !OUTPUT PARAMETERS:

   type(ESMF_GridComp), intent(inout)  :: GC     ! Grid Component
   type(ESMF_State), intent(inout) :: IMPORT     ! Import State
   type(ESMF_State), intent(inout) :: EXPORT     ! Export State
   integer, intent(out) ::  rc                   ! Error return code:
                                                 !  0 - all is well
                                                 !  1 - 

! !DESCRIPTION: This is a simple ESMF wrapper.
!
! !REVISION HISTORY:
!
!  30Nov2010 da Silva  Initial version.
!  26Mar2021 E.Sherman Revised to use ESMF/MAPL

!EOP
!-------------------------------------------------------------------------

    type(GAAS_state), pointer     :: self     ! Legacy state

    type (MAPL_MetaComp), pointer :: MAPL
    type(ESMF_Grid)               :: Grid        ! Grid
    type(ESMF_Config)             :: CF          ! Universal Config 
    type(ESMF_Time)               :: Time     ! Current time
    type(ESMF_Alarm)              :: Alarm
    type(ESMF_Alarm)              :: Predictor_Alarm
    type(ESMF_Alarm)              :: ReplayShutOff_Alarm

    integer                       :: nymd, nhms, i550nm, izAOD, iyAOD
    logical                       :: analysis_time, fexists 
    logical                       :: PREDICTOR_STEP
    logical                       :: ReplayShutOff

    character(len=ESMF_MAXSTR)    :: comp_name

    !(stassi,14feb2012)--character(len=ESMF_MAXSTR)    :: TEMPLATE, filename, expid
    character(len=256)            :: TEMPLATE, filename, expid

    type(ESMF_State)              :: aero ! Aersol state
    character(len=ESMF_MAXSTR)    :: fieldName
    real, pointer, dimension(:,:) :: ptr2d
    real, pointer, dimension(:,:,:) :: ptr3d
    real, dimension(:,:), allocatable, target :: aodInt, aod_a_, aod_k_, aod_f_, &
                                                 y_a_, y_d_
    real, pointer, dimension(:,:,:)  :: DUsum, SSsum, SUsum, NIsum, CAOCsum, CABCsum, CABRsum
    type(ESMF_Field)              :: aod_field
    logical :: skip_analysis 
    integer :: hour,minute,second,year,month,day
    type(ESMF_Time) :: current_time
    type(ESMF_Field) :: aod_a_field

#   include "GAAS_DeclarePointer___.h"

                                __Iam__('Run_')

!  Get pointer for IMPORT/EXPORT/INTERNAL states 
!  ---------------------------------------------
#  include "GAAS_GetPointer___.h"
      
!  Set these exports to UNDEF
!  --------------------------
   if ( associated(aodana) ) aodana(:,:) = MAPL_UNDEF
   if ( associated(aodinc) ) aodinc(:,:) = MAPL_UNDEF

!  Get my name and set-up traceback handle
!  ---------------------------------------
   call ESMF_GridCompGet( GC, name=comp_name, __RC__ )
   Iam = trim(comp_name) // '::Run'

   call MAPL_GetObjectFromGC (GC, MAPL, __RC__)

!  Start timers
!  ------------
   call MAPL_TimerOn( MAPL, "TOTAL")

!  Extract relevant runtime information
!  ------------------------------------
   call extract_ ( GC, CLOCK, self, GRID, CF, time, nymd, nhms, __RC__)
!
   call ESMF_StateGet(import,"aod_a",aod_a_field,_RC)
!  Is it time for analysis?
!  ------------------------
   skip_analysis = ESMFL_field_is_undefined(aod_a_field,_RC)
   analysis_time = .not.skip_analysis

!  Is it time for analysis?
!  ------------------------
   call ESMF_ClockGetAlarm(Clock,'PredictorActive',Predictor_Alarm,__RC__)
   PREDICTOR_STEP = ESMF_AlarmIsRinging( Predictor_Alarm,__RC__)

   call ESMF_ClockGetAlarm(Clock,'ReplayShutOff',ReplayShutOff_Alarm,__RC__)
   ReplayShutOff  = ESMF_AlarmIsRinging( ReplayShutOff_Alarm,__RC__)

!  Stop here if it is NOT analysis time
!  -------------------------------------
   if ( PREDICTOR_STEP .or. ReplayShutOff .or. (.not. analysis_time) ) then
      RETURN_(ESMF_SUCCESS)
   end if

!  OK, let's assimilate the AOD analysis
!  -------------------------------------
   if (MAPL_AM_I_ROOT()) then
      PRINT *, TRIM(Iam)//': Starting Aerosol Assimilation at ', nymd, nhms
      PRINT *,' '
   end if

!  Retrieve AOD from AERO state and fill SimpleBundle
!  ------------------------------------------------------
   call ESMF_StateGet(import, 'AERO', aero, __RC__)

!  Set RH for aerosol optics
   call ESMF_AttributeGet(aero, name='relative_humidity_for_aerosol_optics', value=fieldName, __RC__)
   if (fieldName /= '') then
      call MAPL_GetPointer(aero, ptr3d, trim(fieldName), __RC__)
      ptr3d = rh2
   end if

   call ESMF_AttributeGet(aero, name='air_pressure_for_aerosol_optics', value=fieldName, __RC__)
   if (fieldName /= '') then
      call MAPL_GetPointer(aero, ptr3d, trim(fieldName), __RC__)
      ptr3d = ple
   end if

   ! Set wavelength at 550nm (550 should be a parameter called "monochromatic_wavelength_nm" defined in GAAS)
   call ESMF_AttributeSet(aero, name='wavelength_for_aerosol_optics', value=self%channel*1.0e-9, __RC__)

   ! execute the aero provider's optics method
   call ESMF_MethodExecute(aero, label="get_monochromatic_aop", __RC__)
 
   ! Retrieve vertically summed AOD from AERO
   allocate(aodInt(ubound(rh2,1), ubound(rh2,2)), __STAT__)
   aodInt = 0.0
   call ESMF_AttributeGet(aero, name='monochromatic_extinction_in_air_due_to_ambient_aerosol', value=fieldName, __RC__)
   if (fieldName /= '') then
      call MAPL_GetPointer(aero, ptr2d, trim(fieldName), __RC__)
      aodInt = ptr2d
   end if

!  Set AOD value in y_f
!   self%y_f%r2(1)%name='aod'
   self%y_f%r2(1)%qr4 => aodInt  ! vertically summed AOD
   self%y_f%r2(1)%q => self%y_f%r2(1)%qr4

   if ( self%verbose ) then
       call MAPL_SimpleBundlePrint(self%y_f)
       call MAPL_SimpleBundlePrint(self%q_f)
   end if


!  Read off-line AOD analysis, background and averaging kernel
!  -----------------------------------------------------------    
   self%z_a%r2(1)%qr4 => aod_a    !Move these pointer assignments to Initialize method? -ES
   self%z_a%r2(1)%q => self%z_a%r2(1)%qr4

   self%z_f%r2(1)%qr4 => aod_f
   self%z_f%r2(1)%q => self%z_f%r2(1)%qr4

   self%z_k%r2(1)%qr4 => aod_k
   self%z_k%r2(1)%q => self%z_k%r2(1)%qr4

!  Print summary of input
!  ----------------------
   if ( self%verbose ) then
      call MAPL_SimpleBundlePrint(self%z_f)
      call MAPL_SimpleBundlePrint(self%z_a)
      call MAPL_SimpleBundlePrint(self%z_k)
   end if

!  Convert AOD to Log(AOD+eps) for A.K. Adjustment
!  -----------------------------------------------
   self%z_a%r2(1)%q = Log(self%z_a%r2(1)%q + self%eps)
   self%z_f%r2(1)%q = Log(self%z_f%r2(1)%q + self%eps)
   self%y_f%r2(1)%q = Log(self%y_f%r2(1)%q + self%eps)

!  Background adjustment using averaging kernel
!   This must be done in the Log(AOD+eps) variable
!  -----------------------------------------------
   allocate(y_a_(ubound(rh2,1), ubound(rh2,2)), __STAT__)
   y_a_ = self%z_a%r2(1)%q &
        + (1.-self%z_k%r2(1)%q) &
        * ( self%y_f%r2(1)%q - self%z_f%r2(1)%q )

   self%y_a%r2(1)%qr4 => y_a_
   self%y_a%r2(1)%q => self%y_a%r2(1)%qr4

!  Convert from Log(AOD+eps) back to AOD
!  -------------------------------------
   self%y_a%r2(1)%q = Exp(self%y_a%r2(1)%q) - self%eps
   self%y_f%r2(1)%q = Exp(self%y_f%r2(1)%q) - self%eps

   allocate(y_d_(ubound(rh2,1), ubound(rh2,2)), __STAT__)
   y_d_ = self%y_a%r2(1)%q - self%y_f%r2(1)%q
   self%y_d%r2(1)%qr4 => y_d_
   self%y_d%r2(1)%q => self%y_d%r2(1)%qr4


   if ( self%verbose ) then
       call MAPL_SimpleBundlePrint(self%y_d)
       call MAPL_SimpleBundlePrint(self%y_a)
       call MAPL_SimpleBundlePrint(self%y_f)
   end if


!  Get sum of aerosol mixing ratios
   call get_aerosolSum (aero, 'dust', DUsum, __RC__)
   call get_aerosolSum (aero, 'seasalt', SSsum, __RC__)
   call get_aerosolSum (aero, 'sulfate', SUsum, __RC__)
   call get_aerosolSum (aero, 'nitrate', NIsum, __RC__)
   call get_aerosolSum (aero, 'organicCarbon', CAOCsum, __RC__)
   call get_aerosolSum (aero, 'blackCarbon', CABCsum, __RC__)
   call get_aerosolSum (aero, 'brownCarbon', CABRsum, __RC__)

!  Handle 3D exports (save bkg for increments)
!  -------------------------------------------
   if ( associated(duinc) ) duinc = DUsum
   if ( associated(ssinc) ) ssinc = SSsum
   if ( associated(niinc) ) niinc = NIsum
   if ( associated(bcinc) ) bcinc = CABCsum
   if ( associated(ocinc) ) ocinc = CAOCsum
   if ( associated(brinc) ) brinc = CABRsum
   if ( associated(suinc) ) suinc = SUsum

!  Create concetration analysis from AOD increments
!   Here we pass in the y_f and y_d in terms of AOD,
!   *not* Log(AOD+eps)
!  ------------------------------------------------
   call LDE_Projector1c ( self%E, self%q_a, self%q_f, self%y_f, self%y_d, self%verbose, __RC__ )

!  Get updated sum of aerosol mixing ratios
   call get_aerosolSum (aero, 'dust', DUsum, __RC__)
   call get_aerosolSum (aero, 'seasalt', SSsum, __RC__)
   call get_aerosolSum (aero, 'sulfate', SUsum, __RC__)
   call get_aerosolSum (aero, 'nitrate', NIsum, __RC__)
   call get_aerosolSum (aero, 'organicCarbon', CAOCsum, __RC__)
   call get_aerosolSum (aero, 'blackCarbon', CABCsum, __RC__)
   call get_aerosolSum (aero, 'brownCarbon', CABRsum, __RC__)

!  Handle 2D exports
!  -----------------
   if ( associated(aodana) ) aodana(:,:) = self%y_a%r2(1)%q(:,:)
   if ( associated(aodinc) ) aodinc(:,:) = self%y_d%r2(1)%q(:,:)

!  Handle 3D exports
!  -----------------
   if ( associated(duana) ) duana = DUsum
   if ( associated(ssana) ) ssana = SSsum
   if ( associated(niana) ) niana = NIsum
   if ( associated(bcana) ) bcana = CABCsum
   if ( associated(ocana) ) ocana = CAOCsum
   if ( associated(brana) ) brana = CABRsum
   if ( associated(suana) ) suana = SUsum

!  Compute increments
   if ( associated(duinc) ) duinc = DUsum - duinc
   if ( associated(ssinc) ) ssinc = SSsum - ssinc
   if ( associated(niinc) ) niinc = NIsum - niinc
   if ( associated(bcinc) ) bcinc = CABCsum - bcinc
   if ( associated(ocinc) ) ocinc = CAOCsum - ocinc
   if ( associated(brinc) ) brinc = CABRsum - brinc
   if ( associated(suinc) ) suinc = SUsum - suinc


#ifdef PRINT_STATES
   if (MAPL_AM_I_ROOT()) then
       print *, trim(Iam)//': IMPORT   State during Run():'
       call ESMF_StatePrint ( IMPORT )
       print *, trim(Iam)//': EXPORT   State during Run():' 
       call ESMF_StatePrint ( EXPORT )
    end if
#endif

!  Stop timers
!  ------------
   call MAPL_TimerOff( MAPL, "TOTAL")

if (mapl_am_i_root()) print*,'GAAS finished!'

!  All done
!  --------
   RETURN_(ESMF_SUCCESS)

Contains

   subroutine get_aerosolSum (state, aeroName, aeroSum, rc)

     implicit none
   
     !ARGUMENTS:
     type (ESMF_State),               intent(inout)    :: state
     character (len=*),               intent(in)       :: aeroName
     real, dimension(:,:,:), pointer, intent(out)      :: aeroSum
     integer, optional,               intent(out)      :: rc

     !LOCALS:
     integer                                          :: status
     character (len=ESMF_MAXSTR)                      :: fld_name
     character (len=ESMF_MAXSTR)                      :: aeroToken

     !Begin...

     select case (aeroName)
        case ('dust')
           aeroToken = 'DU'
        case ('seasalt')
           aeroToken = 'SS'
        case ('sulfate')
           aeroToken = 'SU'
        case ('nitrate')
           aeroToken = 'NI'
        case ('organicCarbon')
           aeroToken = 'CA.oc'
        case ('blackCarbon')
           aeroToken = 'CA.bc'
        case ('brownCarbon')
           aeroToken = 'CA.br'
     end select

     ! Set aerosol to retrieve sum for
     call ESMF_AttributeSet(state, name='aerosolName', value=trim(aeroName), __RC__)

     ! execute the aero provider's optics method
     call ESMF_MethodExecute(state, label="get_mixRatioSum", __RC__)

     call ESMF_AttributeGet(state, name='sum_of_internalState_aerosol_'//trim(aeroToken), value=fieldName, __RC__)
     if (fieldName /= '') then
        call MAPL_GetPointer(state, aeroSum, trim(fieldName), __RC__)
     end if

   end subroutine get_aerosolSum

   END SUBROUTINE Run_

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Finalize_ --- Finalize GAAS
!
! !INTERFACE:
!

   SUBROUTINE Finalize_ ( GC, IMPORT, EXPORT, CLOCK, rc )

! !USES:

  implicit NONE

! !INPUT PARAMETERS:

   type(ESMF_Clock),  intent(inout) :: CLOCK      ! The clock

! !OUTPUT PARAMETERS:

   type(ESMF_GridComp), intent(inout)  :: gc     ! Grid Component
   type(ESMF_State), intent(inout) :: IMPORT     ! Import State
   type(ESMF_State), intent(inout) :: EXPORT     ! Export State
   integer, intent(out) ::  rc                   ! Error return code:
                                                 !  0 - all is well
                                                 !  1 - 

! !DESCRIPTION: This is a simple ESMF wrapper.
!
! !REVISION HISTORY:
!
!  30Nov2010 da Silva  Initial version.
!
!EOP
!-------------------------------------------------------------------------

    type(GAAS_state), pointer     :: self     ! Legacy state
    type(ESMF_Grid)               :: Grid     ! Grid
    type(ESMF_Config)             :: CF       ! Universal Config 
    type(ESMF_Time)               :: Time     ! Current time

    integer                       :: nymd, nhms  ! date, time

    character(len=ESMF_MAXSTR)    :: COMP_NAME

#   include "GAAS_DeclarePointer___.h"

                                   __Iam__('Finalize_')
   
!  Declare and get pointer for IMPORT/EXPORT/INTERNAL states 
!  ---------------------------------------------------------
#  include "GAAS_GetPointer___.h"

!  Get my name and set-up traceback handle
!  ---------------------------------------
   call ESMF_GridCompGet( GC, name=comp_name, __RC__ )
   Iam = trim(comp_name) // '::Finalize'

!  Initialize MAPL Generic
!  -----------------------
   call MAPL_GenericFinalize ( GC, IMPORT, EXPORT, CLOCK,  __RC__ )

!  Extract relevant runtime information
!  ------------------------------------
   call extract_ ( GC, CLOCK, self, GRID, CF, time, nymd, nhms, __RC__)

!  Destroy LDE object
!  ------------------
   call LDE_Destroy ( self%E, __RC__ )

!  Destroy all simple bundles
!  --------------------------

   call MAPL_SimpleBundleDestroy(self%y_f, __RC__)
   call MAPL_SimpleBundleDestroy(self%y_a, __RC__)
   call MAPL_SimpleBundleDestroy(self%y_d, __RC__)

!ALT   call MAPL_SimpleBundleDestroy(self%q_f, __RC__)
!ALT   call MAPL_SimpleBundleDestroy(self%q_a, __RC__)
 
!  All done
!  --------
   RETURN_(ESMF_SUCCESS)

 end SUBROUTINE Finalize_

!.......................................................................

    subroutine extract_ ( GC, CLOCK, self, GRID, CF, time, nymd, nhms, rc)

    type(ESMF_GridComp), intent(INout)  :: GC           ! Grid Comp object
    type(ESMF_Clock), intent(in)        :: CLOCK        ! Clock

    type(GAAS_state), pointer           :: self         ! Legacy state
    type(ESMF_Grid),     intent(out)    :: GRID         ! Grid
    type(ESMF_Config),   intent(out)    :: CF           ! Universal Config 
    type(ESMF_TIME), intent(out)        :: Time         ! Time
    type(ESMF_TimeInterval)             :: TimeStep     ! used to define a clock
    integer, intent(out)                :: nymd, nhms   ! date, time
    integer, intent(out), optional      :: rc

!                                      ---

    character(len=ESMF_MAXSTR) :: comp_name
    
                                 __Iam__('extract_')

    type(MAPL_MetaComp), pointer  :: MC
    type(GAAS_Wrap)               :: wrap
    integer                       :: iyr, imm, idd, ihr, imn, isc

!   Get my name and set-up traceback handle
!   ---------------------------------------
    call ESMF_GridCompGet( GC, NAME=comp_name, __RC__ )
    Iam = trim(COMP_NAME) // '::extract_'

    rc = 0

!   Get my internal state
!   ---------------------
    call ESMF_UserCompGetInternalState(gc, 'GAAS_state', WRAP, STATUS)
    VERIFY_(STATUS)
    self => wrap%ptr

!   Get the configuration
!   ---------------------
    call ESMF_GridCompGet ( GC, config=CF, __RC__ )

!   Extract time as simple integers from clock
!   ------------------------------------------
    call ESMF_ClockGet(CLOCK, timeStep=TimeStep, currTIME=TIME,__RC__)

!   NOTE: we shift the time back one time step because the clock ticks
!   prior to writing HISTORY, so now we will write at 0Z what is intended
!   to be the analysis at 0Z.
    TIME = TIME + TimeStep

    call ESMF_TimeGet(TIME ,yy=iyr, mm=imm, dd=idd, h=ihr, m=imn, s=isc, __RC__)
    call MAPL_PackTime(nymd,iyr,imm,idd)
    call MAPL_PackTime(nhms,ihr,imn,isc)

!   Extract the ESMF Grid
!   ---------------------
    call ESMF_GridCompGet ( GC, grid=GRID, __RC__)

    RETURN_(ESMF_SUCCESS)

   end subroutine extract_

 END MODULE GAAS_GridCompMod
