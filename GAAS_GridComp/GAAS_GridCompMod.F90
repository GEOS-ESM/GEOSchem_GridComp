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
!
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
    call ESMF_ConfigGetAttribute(self%CF, self%no_fuss, Label='no_fuss_if_ana_missing:',__RC__ )
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
!  30Nov2010 da Silva  Initial version.
!
!EOP
!-------------------------------------------------------------------------

    type(GAAS_state), pointer     :: self     ! Legacy state
    type(ESMF_Grid)               :: Grid     ! Grid
    type(ESMF_Config)             :: CF       ! Universal Config 
    type(ESMF_Time)               :: Time     ! Current time

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
   self%aerReg = Chem_RegistryCreate ( rcfile='GAAS_AerRegistry.rc', __RC__ )
   self%aodReg = Chem_RegistryCreate ( rcfile='GAAS_AodRegistry.rc', __RC__ )

!  Mie tables, etc
!  ---------------
   self%Mie = Chem_MieCreate(self%CF, chemReg=self%aerReg, __RC__)

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

!  Create AOD Simple Bundles
!  -------------------------
   self%y_f = Chem_SimpleBundleCreate ('aod_bkg', self%aodReg, self%aodGrid, &
                                        Levs=1.e9*self%Mie%channels,         &
                                        LevUnits="nm", delp=null(), __RC__) 
   self%y_a = Chem_SimpleBundleCreate ('aod_ana', self%aodReg, self%aodGrid, &
                                        Levs=1.e9*self%Mie%channels,         &
                                        LevUnits="nm", delp=null(), __RC__) 
   self%y_d = Chem_SimpleBundleCreate ('aod_inc', self%aodReg, self%aodGrid, &
                                        Levs=1.e9*self%Mie%channels,         &
                                        LevUnits="nm", delp=null(), __RC__) 

!  Associate Import state with concentration analysis/background
!  -------------------------------------------------------------
   self%q_f = MAPL_SimpleBundleCreate(IMPORT,__RC__)
   self%q_a = MAPL_SimpleBundleCreate(IMPORT,__RC__) ! Check for Friendlies!!!

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
!
!EOP
!-------------------------------------------------------------------------

    type(GAAS_state), pointer     :: self     ! Legacy state

    type (MAPL_MetaComp), pointer :: MAPL
    type(ESMF_Grid)               :: Grid        ! Grid
    type(ESMF_Config)             :: CF          ! Universal Config 
    type(ESMF_Time)               :: Time     ! Current time
    type(ESMF_Alarm)              :: Alarm
    type(ESMF_Alarm)              :: Predictor_Alarm

    integer                       :: nymd, nhms, i550nm, izAOD, iyAOD
    logical                       :: analysis_time, fexists 
    logical                       :: PREDICTOR_STEP

    character(len=ESMF_MAXSTR)    :: comp_name

    !(stassi,14feb2012)--character(len=ESMF_MAXSTR)    :: TEMPLATE, filename, expid
    character(len=256)            :: TEMPLATE, filename, expid

#   include "GAAS_DeclarePointer___.h"

                                __Iam__('Run_')


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

!  Is it time for analysis?
!  ------------------------
   call MAPL_Get (MAPL, RUNALARM=Alarm, __RC__)
#if 0
   analysis_time = ESMF_AlarmIsRinging(Alarm,__RC__)
#endif

   call ESMF_ClockGetAlarm(Clock,'PredictorActive',Predictor_Alarm,__RC__)
   PREDICTOR_STEP = ESMF_AlarmIsRinging( Predictor_Alarm,__RC__)

!  For some reason the alarm above is not working.
!  For now, hardwire this...
!  -----------------------------------------------
   analysis_time =  nhms ==      0 .OR. &
                    nhms ==  30000 .OR. &
                    nhms ==  60000 .OR. &
                    nhms ==  90000 .OR. &
                    nhms == 120000 .OR. &
                    nhms == 150000 .OR. &
                    nhms == 180000 .OR. &
                    nhms == 210000

!  Stop here if it is NOT analysis time
!  -------------------------------------
   if ( PREDICTOR_STEP .or. (.not. analysis_time) ) then
      RETURN_(ESMF_SUCCESS)
   end if

!  If desired, just bail out if analysis file is not there
!  -------------------------------------------------------
   if ( self%no_fuss ) then
        call ESMF_ConfigGetAttribute(self%CF, expid, Label='EXPID:', default='unknown', __RC__)
        call ESMF_ConfigGetAttribute(self%CF, template, Label='aod_ana_filename:', __RC__)
        call StrTemplate ( filename, template, xid=expid, nymd=nymd, nhms=nhms )        
        inquire ( file=trim(filename), exist=fexists )
        if ( .not. fexists ) then
           if (MAPL_AM_I_ROOT()) then
              PRINT *, TRIM(Iam)//': Ignoring Aerosol Assimilation at ', nymd, nhms
              PRINT *, TRIM(Iam)//': Cannot find AOD analysis file ', trim(filename)
              PRINT *,' '
           end if
           RETURN_(ESMF_SUCCESS)
        end if
     end if

!  OK, let's assimilate the AOD analysis
!  -------------------------------------
   if (MAPL_AM_I_ROOT()) then
      PRINT *, TRIM(Iam)//': Starting Aerosol Assimilation at ', nymd, nhms
      PRINT *,' '
   end if

!  Run MAPL Generic
!  ----------------
!ALT   call MAPL_GenericRunChildren ( gc, IMPORT, EXPORT, clock,  __RC__ )

!  Get pointer for IMPORT/EXPORT/INTERNAL states 
!  ---------------------------------------------
#  include "GAAS_GetPointer___.h"

!  Calculate on-line AOD
!  ---------------------
   call Chem_AodCalculator (self%y_f, self%q_f, self%Mie, self%verbose, __RC__)

   if ( self%verbose ) then
       call MAPL_SimpleBundlePrint(self%y_f)
       call MAPL_SimpleBundlePrint(self%q_f)
   end if

!  Read off-line AOD analysis, background and averaging kernel
!  ----------------------------------------------------------- 
   self%z_f = Chem_SimpleBundleRead (self%CF, 'aod_bkg_filename', self%aodGrid, & 
                                     time=Time, only_vars='AOD', __RC__ )
   self%z_a = Chem_SimpleBundleRead (self%CF, 'aod_ana_filename', self%aodGrid, &
                                     time=Time, only_vars='AOD',  __RC__ )
   self%z_k = Chem_SimpleBundleRead (self%CF, 'aod_avk_filename', self%aodGrid, &
                                     time=Time, only_vars='AOD', __RC__ )

!  Print summary of input
!  ----------------------
   if ( self%verbose ) then
       call MAPL_SimpleBundlePrint(self%z_f)
       call MAPL_SimpleBundlePrint(self%z_a)
       call MAPL_SimpleBundlePrint(self%z_k)
   end if

!  Because the off-line analysis may have other fields in the Bundle,
!  we explicitly look for the AOD index to be safe
!  -----------------------------------------------------------------
   izAOD = MAPL_SimpleBundleGetIndex(self%z_f,'AOD',3,__RC__)
   iyAOD = MAPL_SimpleBundleGetIndex(self%y_f,'AOD',3,__RC__)
   _ASSERT(iyAOD==1,'needs informative message') ! what we have created must have only AOD

!  Convert AOD to Log(AOD+eps) for A.K. Adjustment
!  -----------------------------------------------
   self%z_a%r3(izAOD)%q = Log(self%z_a%r3(izAOD)%q + self%eps)
   self%z_f%r3(izAOD)%q = Log(self%z_f%r3(izAOD)%q + self%eps)
   self%y_f%r3(iyAOD)%q = Log(self%y_f%r3(iyAOD)%q + self%eps)

!  Background adjustment using averaging kernel
!   This must be done in the Log(AOD+eps) variable
!  -----------------------------------------------
   self%y_a%r3(iyAOD)%q = self%z_a%r3(izAOD)%q &
                        + (1.-self%z_k%r3(izAOD)%q) &
                        * ( self%y_f%r3(iyAOD)%q - self%z_f%r3(izAOD)%q )

!  Convert from Log(AOD+eps) back to AOD
!  -------------------------------------
   self%y_a%r3(iyAOD)%q = Exp(self%y_a%r3(iyAOD)%q) - self%eps
   self%y_f%r3(iyAOD)%q = Exp(self%y_f%r3(iyAOD)%q) - self%eps
   self%y_d%r3(iyAOD)%q = self%y_a%r3(iyAOD)%q - self%y_f%r3(iyAOD)%q 

   if ( self%verbose ) then
       call MAPL_SimpleBundlePrint(self%y_d)
   end if

!  Handle 3D exports (save bkg for increments)
!  -------------------------------------------
   if ( associated(duinc) ) duinc = du001+du002+du003+du004+du005
   if ( associated(ssinc) ) ssinc = ss001+ss002+ss003+ss004+ss005
   if ( associated(niinc) ) niinc = no3an1+no3an2+no3an3
   if ( associated(bcinc) ) bcinc = bcphobic + bcphilic
   if ( associated(ocinc) ) ocinc = ocphobic + ocphilic
   if ( associated(suinc) ) suinc = so4

!  Create concetration analysis from AOD increments
!   Here we pass in the y_f and y_d in terms of AOD,
!   *not* Log(AOD+eps)
!  ------------------------------------------------
   call LDE_Projector1c ( self%E, self%q_a, self%q_f, self%y_f, self%y_d, self%verbose, __RC__ )

!  Handle 2D exports
!  -----------------
   i550nm = getChannel_(self%y_a%coords%levs,__RC__)
   if ( associated(aodana) ) aodana(:,:) = self%y_a%r3(iyAOD)%q(:,:,i550nm)
   if ( associated(aodinc) ) aodinc(:,:) = self%y_d%r3(iyAOD)%q(:,:,i550nm)

!  Handle 3D exports
!  -----------------
   if ( associated(duana) ) duana = du001+du002+du003+du004+du005
   if ( associated(ssana) ) ssana = ss001+ss002+ss003+ss004+ss005
   if ( associated(niana) ) niana = no3an1+no3an2+no3an3
   if ( associated(bcana) ) bcana = bcphobic + bcphilic
   if ( associated(ocana) ) ocana = ocphobic + ocphilic
   if ( associated(suana) ) suana = so4

   if ( associated(duinc) ) duinc = du001+du002+du003+du004+du005 - duinc
   if ( associated(ssinc) ) ssinc = ss001+ss002+ss003+ss004+ss005 - ssinc
   if ( associated(niinc) ) niinc = no3an1+no3an2+no3an3 - niinc
   if ( associated(bcinc) ) bcinc = bcphobic + bcphilic - bcinc
   if ( associated(ocinc) ) ocinc = ocphobic + ocphilic - ocinc
   if ( associated(suinc) ) suinc = so4 - suinc

#ifdef PRINT_STATES

   if (MAPL_AM_I_ROOT()) then
       print *, trim(Iam)//': IMPORT   State during Run():'
       call ESMF_StatePrint ( IMPORT )
       print *, trim(Iam)//': EXPORT   State during Run():' 
       call ESMF_StatePrint ( EXPORT )
    end if

#endif

! Clean-up
! --------
    call MAPL_SimpleBundleDestroy(self%z_f, __RC__)
    call MAPL_SimpleBundleDestroy(self%z_a, __RC__)
    call MAPL_SimpleBundleDestroy(self%z_k, __RC__)

!  Stop timers
!  ------------
   call MAPL_TimerOff( MAPL, "TOTAL")

!  All done
!  --------
   RETURN_(ESMF_SUCCESS)

Contains

   function getChannel_(levs,rc) result (i)
     integer :: i, j, rc
     real :: levs(:)
     i = -1
     do j = 1, size(levs)
        if ( abs(levs(j)-self%channel) < 0.1 ) i = j
     end do
     if ( i<1 ) then
          rc = 1
     else       
          rc = 0
     end if
     return
        
   end function getChannel_

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
