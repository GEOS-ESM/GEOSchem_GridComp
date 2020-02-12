#include "MAPL_Generic.h"

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: DNA_GridCompMod - 
!
! !INTERFACE:
!
   module DNA_GridCompMod
!
! !USES:
!
   use ESMF
   use MAPL


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
  type DNA_State
     private
     type(ESMF_Config)   :: CF                        ! Private Config

     type(ESMF_Grid)     :: grid                      ! Grid

     logical             :: verbose                   ! turn on/off more verbose messages
  end type DNA_State

! Hook for the ESMF
! -----------------
  type DNA_Wrap
     type (DNA_State), pointer :: PTR => null()
  end type DNA_Wrap

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
    type (DNA_State), pointer  :: self   ! internal state
    type (DNA_Wrap)            :: wrap

    character(len=ESMF_MAXSTR) :: comp_name

!   Local variables
!   --------------------------
    character(len=ESMF_MAXSTR) :: field_name
    character(len=ESMF_MAXSTR) :: field_long_name
    integer                    :: band

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
    call ESMF_ConfigLoadFile(self%CF, 'DNA_GridComp.rc', __RC__)

    call ESMF_ConfigGetAttribute(self%CF, self%verbose, Label='verbose:', default=.false., __RC__)


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
    call ESMF_UserCompSetInternalState(GC, 'DNA_State', wrap, STATUS)
    VERIFY_(STATUS)
  
!                         ------------------
!                         MAPL Data Services
!                         ------------------

!BOS
!
! !IMPORT STATE:

   do band = 1, 18
       write (field_name, "('AEROSOL_EXT_BAND', I2.2)") band
       write (field_long_name, "('aerosol_extinction_for_band_', I2.2)") band
       call MAPL_AddImportSpec(GC,                         &
           SHORT_NAME = trim(field_name),                  &
           LONG_NAME  = trim(field_long_name),             &
           UNITS      = '',                                &
           DIMS       = MAPL_DimsHorzVert,                 &
           VLOCATION  = MAPL_VLocationCenter,              &
           __RC__)

       write (field_name, "('AEROSOL_SSA_BAND', I2.2)") band
       write (field_long_name, "('aerosol_single_scattering_albedo_for_band_', I2.2)") band
       call MAPL_AddImportSpec(GC,                         &
           SHORT_NAME = trim(field_name),                  &
           LONG_NAME  = trim(field_long_name),             &
           UNITS      = '',                                &
           DIMS       = MAPL_DimsHorzVert,                 &
           VLOCATION  = MAPL_VLocationCenter,              &
           __RC__)

       write (field_name, "('AEROSOL_ASY_BAND', I2.2)") band
       write (field_long_name, "('aerosol_asymmetry_parameter_for_band_', I2.2)") band
       call MAPL_AddImportSpec(GC,                         &
           SHORT_NAME = trim(field_name),                  &
           LONG_NAME  = trim(field_long_name),             &
           UNITS      = '',                                &
           DIMS       = MAPL_DimsHorzVert,                 &
           VLOCATION  = MAPL_VLocationCenter,              &
           __RC__)
   end do


! !INTERNAL STATE:
   do band = 1, 18
       write (field_name, "('DNA::AEROSOL_EXT_BAND', I2.2)") band
       write (field_long_name, "('aerosol_extinction_for_band_', I2.2)") band
       call MAPL_AddInternalSpec(GC,                       &
           SHORT_NAME = trim(field_name),                  &
           LONG_NAME  = trim(field_long_name),             &
           UNITS      = '',                                &
           DIMS       = MAPL_DimsHorzVert,                 &
           VLOCATION  = MAPL_VLocationCenter,              &
           ADD2EXPORT = .true.,                            &
           RESTART    = MAPL_RestartSkip,                  &
           __RC__)

       write (field_name, "('DNA::AEROSOL_SSA_BAND', I2.2)") band
       write (field_long_name, "('aerosol_single_scattering_albedo_for_band_', I2.2)") band
       call MAPL_AddInternalSpec(GC,                       &
           SHORT_NAME = trim(field_name),                  &
           LONG_NAME  = trim(field_long_name),             &
           UNITS      = '',                                &
           DIMS       = MAPL_DimsHorzVert,                 &
           VLOCATION  = MAPL_VLocationCenter,              &
           ADD2EXPORT = .true.,                            &
           RESTART    = MAPL_RestartSkip,                  &
           __RC__)

       write (field_name, "('DNA::AEROSOL_ASY_BAND', I2.2)") band
       write (field_long_name, "('aerosol_asymmetry_parameter_for_band_', I2.2)") band
       call MAPL_AddInternalSpec(GC,                       &
           SHORT_NAME = trim(field_name),                  &
           LONG_NAME  = trim(field_long_name),             &
           UNITS      = '',                                &
           DIMS       = MAPL_DimsHorzVert,                 &
           VLOCATION  = MAPL_VLocationCenter,              &
           ADD2EXPORT = .true.,                            &
           RESTART    = MAPL_RestartSkip,                  &
           __RC__)
   end do



! !EXTERNAL STATE:

!   This state is needed by radiation - It will contain 
!   aerosol number and mass mixing ratios, and aerosol optics
!   ---------------------------------------------------------
    call MAPL_AddExportSpec(GC,                            &
        SHORT_NAME = 'AERO',                               &
        LONG_NAME  = 'aerosol_mixing_ratios',              &
        UNITS      = 'kg kg-1',                            &
        DIMS       = MAPL_DimsHorzVert,                    &
        VLOCATION  = MAPL_VLocationCenter,                 &
        DATATYPE   = MAPL_StateItem,                       &
        __RC__)

!   This bundle is not filled in by DNA, just a place holder for now
!   ----------------------------------------------------------------
    call MAPL_AddExportSpec(GC,                            &
        SHORT_NAME = 'AERO_DP',                            &
        LONG_NAME  = 'aerosol_deposition',                 &
        UNITS      = 'kg m-2 s-1',                         &
        DIMS       = MAPL_DimsHorzOnly,                    &
        DATATYPE   = MAPL_BundleItem, __RC__)

!EOS

!   Set the Profiling timers
!   ------------------------
    call MAPL_TimerAdd(GC, name = "RUN",        __RC__)
    call MAPL_TimerAdd(GC, name = "INITIALIZE", __RC__)
    call MAPL_TimerAdd(GC, name = "FINALIZE",   __RC__)


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

    type(DNA_State), pointer        :: self        ! Legacy state
    type(ESMF_Grid)                 :: GRID        ! Grid
    type(ESMF_Config)               :: CF          ! Universal Config 

    type(MAPL_MetaComp), pointer    :: MetaComp
    type(ESMF_State)                :: INTERNAL

    integer                         :: i1, i2, im  ! 3D Dimensions
    integer                         :: j1, j2, jm  !
    integer                         :: km          !

    integer                         :: nymd, nhms  ! date, time
    real                            :: cdt         ! time step in secs

    character(len=ESMF_MAXSTR)      :: comp_name   ! component's name

    type(ESMF_State)                :: aero
    logical                         :: implements_aerosol_optics
    type(ESMF_FieldBundle)          :: radiation

    type(ESMF_Field)                :: field
    character(len=ESMF_MAXSTR)      :: field_name

    integer                         :: band


!  Get my name and set-up traceback handle
!  ---------------------------------------
   call ESMF_GridCompGet(GC, name=COMP_NAME, __RC__)
   Iam = trim(comp_name) // '::' // trim(Iam)

!                               --------
   if (MAPL_AM_I_ROOT()) then
      print *, trim(Iam)//': Starting...'
      print *, ''
   end if


!  Initialize MAPL Generic
!  -----------------------
   call MAPL_GenericInitialize(GC, IMPORT, EXPORT, clock, __RC__)


!  Get a pointer to the INTERNAL state
!  -----------------------------------
   call MAPL_GetObjectFromGC( GC, MetaComp, __RC__)
   call MAPL_Get(MetaComp, INTERNAL_ESMF_STATE=INTERNAL, __RC__)


!  Extract relevant runtime information
!  ------------------------------------
   call extract_(GC, CLOCK, self, GRID, CF, i1, i2, im, j1, j2, jm, km, nymd, nhms, cdt, __RC__)


!  Set the grid
!  -------------------------------------------------
   self%grid = GRID


!  Fill the AERO State
!  ------------------------------------------------
   call ESMF_StateGet(EXPORT, 'AERO', aero, __RC__)

   ! This attribute indicates if the aerosol optics method is implemented or not. 
   ! Radiation will not call the aerosol optics method unless this attribute is 
   ! explicitly set to true.

   implements_aerosol_optics = .true. 

   call ESMF_AttributeSet(aero, name  = 'implements_aerosol_optics_method', & 
                                value = implements_aerosol_optics, __RC__)
  
   COUPLING_TO_RADIATION: if (implements_aerosol_optics) then
       
       radiation = ESMF_FieldBundleCreate(name='RADIATION_PROPERTIES', __RC__)
       call MAPL_StateAdd(aero, radiation, __RC__)

       do band = 1, 18
           write (field_name, "('DNA::AEROSOL_EXT_BAND', I2.2)") band
           call ESMF_StateGet(internal, trim(field_name), field, __RC__)
           call ESMF_FieldBundleAdd(radiation, (/field/), __RC__)

           write (field_name, "('DNA::AEROSOL_SSA_BAND', I2.2)") band
           call ESMF_StateGet(internal, trim(field_name), field, __RC__)
           call ESMF_FieldBundleAdd(radiation, (/field/), __RC__)

           write (field_name, "('DNA::AEROSOL_ASY_BAND', I2.2)") band
           call ESMF_StateGet(internal, trim(field_name), field, __RC__)
           call ESMF_FieldBundleAdd(radiation, (/field/), __RC__)
       end do


       ! state of the atmosphere
       call ESMF_AttributeSet(aero, name='air_pressure_for_aerosol_optics',             value='',    __RC__)
       call ESMF_AttributeSet(aero, name='relative_humidity_for_aerosol_optics',        value='',    __RC__)
       call ESMF_AttributeSet(aero, name='cloud_area_fraction_for_aerosol_optics',      value='',    __RC__) ! 'cloud_area_fraction_in_atmosphere_layer_for_aerosol_optics'

       ! aerosol optics
       call ESMF_AttributeSet(aero, name='band_for_aerosol_optics',                     value=0,     __RC__)
       call ESMF_AttributeSet(aero, name='extinction_in_air_due_to_ambient_aerosol',    value='EXT', __RC__)
       call ESMF_AttributeSet(aero, name='single_scattering_albedo_of_ambient_aerosol', value='SSA', __RC__)
       call ESMF_AttributeSet(aero, name='asymmetry_parameter_of_ambient_aerosol',      value='ASY', __RC__)

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


!  All done
!  --------
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
   
   type(DNA_State), pointer       :: self        ! Legacy state
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

!  Input fields
!  ------------
   type(ESMF_State)                :: aero

   real, pointer, dimension(:,:,:) :: ext_
   real, pointer, dimension(:,:,:) :: ssa_
   real, pointer, dimension(:,:,:) :: asy_

!  Export fields
!  -------------
   type(ESMF_State)                :: internal

   real, pointer, dimension(:,:,:) :: ext
   real, pointer, dimension(:,:,:) :: ssa
   real, pointer, dimension(:,:,:) :: asy


!  local
!  -----
   real, parameter :: ORO_OCEAN   = 0.0
   real, parameter :: ORO_LAND    = 1.0
   real, parameter :: ORO_SEA_ICE = 2.0

   integer :: band

   character(len=ESMF_MAXSTR) :: field_name

!  Get my name and set-up traceback handle
!  ---------------------------------------
   call ESMF_GridCompGet(GC, name=comp_name, __RC__)
   Iam = trim(comp_name) // '::' // trim(Iam)


!  Get my internal MAPL_Generic state
!  -----------------------------------
   call MAPL_GetObjectFromGC(GC, mgState, __RC__)

!  Get a pointer to the INTERNAL state
!  -----------------------------------
   call MAPL_Get(mgState, INTERNAL_ESMF_STATE=INTERNAL, __RC__)

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


!  Update the 'AERO' state
!  -----------------------
   call ESMF_StateGet(export, 'AERO', aero, __RC__) 
   
   do band = 1, 18

       ! update aerosol extinction for this band
       write (field_name, "('AEROSOL_EXT_BAND', I2.2)") band
       call MAPL_GetPointer(import, ext_, trim(field_name), __RC__)
       write (field_name, "('DNA::AEROSOL_EXT_BAND', I2.2)") band
       call MAPL_GetPointer(internal, ext, trim(field_name), __RC__)

       ext = ext_

       ! update aerosol single scattering albedo for this band
       write (field_name, "('AEROSOL_SSA_BAND', I2.2)") band
       call MAPL_GetPointer(import, ssa_, trim(field_name), __RC__)
       write (field_name, "('DNA::AEROSOL_SSA_BAND', I2.2)") band
       call MAPL_GetPointer(internal, ssa, trim(field_name), __RC__)

       ssa = ssa_

       ! update aerosol asymmetry parameter for this band
       write (field_name, "('AEROSOL_ASY_BAND', I2.2)") band
       call MAPL_GetPointer(import, asy_, trim(field_name), __RC__)
       write (field_name, "('DNA::AEROSOL_ASY_BAND', I2.2)") band
       call MAPL_GetPointer(internal, asy, trim(field_name), __RC__)

       asy = asy_
   end do



!  All done
!  --------
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
   
    type(DNA_State), pointer       :: self        ! Legacy state
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


!  Free memory for the internal the private state
!  ----------------------------------------------
   deallocate(self, __STAT__)

!  All done
!  --------
   RETURN_(ESMF_SUCCESS)

 end SUBROUTINE Finalize_

!.......................................................................

 subroutine aerosol_optics(state, rc)

  implicit none 


! Arguments
! ---------
  type(ESMF_State):: state
  integer, intent(out):: rc


! Local
! ---------
  type(ESMF_FieldBundle)          :: radiation
  character(len=ESMF_MAXSTR)      :: field_name
  real, pointer, dimension(:,:,:) :: ext, ssa, asy   ! 3D total aerosol optical properties
  real, pointer, dimension(:,:,:) :: ptr3d

  integer :: band
  
  integer :: STATUS
  character(len=ESMF_MAXSTR) :: Iam



                                   Iam = 'GEOS_DNA::aerosol_optics()'


! Radiation band
! --------------
  call ESMF_AttributeGet(state, name='band_for_aerosol_optics', value=band, __RC__)

! Field bundle with radiation properties over all bands 
! -----------------------------------------------------
  call ESMF_StateGet(state, 'RADIATION_PROPERTIES', radiation, __RC__)


! Get the band specific properties
! --------------------------------
  write (field_name, "('DNA::AEROSOL_EXT_BAND', I2.2)") band
  call ESMFL_BundleGetPointerToData(radiation, trim(field_name), ext, __RC__)

  write (field_name, "('DNA::AEROSOL_SSA_BAND', I2.2)") band
  call ESMFL_BundleGetPointerToData(radiation, trim(field_name), ssa, __RC__)

  write (field_name, "('DNA::AEROSOL_ASY_BAND', I2.2)") band
  call ESMFL_BundleGetPointerToData(radiation, trim(field_name), asy, __RC__)


  ! inputs for radiation:  ext, sca and sca*asy
  ! in the callback     : 'asy' = product of asy and sca
  !                       'ssa' = sca


  call ESMF_AttributeGet(state, name='extinction_in_air_due_to_ambient_aerosol', value=field_name, __RC__)
  if (field_name /= '') then 
      call MAPL_GetPointer(state, ptr3d, trim(field_name), __RC__)
      ptr3d = ext
  end if

  call ESMF_AttributeGet(state, name='single_scattering_albedo_of_ambient_aerosol', value=field_name, __RC__)
  if (field_name /= '') then 
      call MAPL_GetPointer(state, ptr3d, trim(field_name), __RC__)
      ptr3d = ssa * ext
  end if

  call ESMF_AttributeGet(state, name='asymmetry_parameter_of_ambient_aerosol', value=field_name, __RC__)
  if (field_name /= '') then 
      call MAPL_GetPointer(state, asy, trim(field_name), __RC__)
      ptr3d = ssa * ext * asy
  end if


  RETURN_(ESMF_SUCCESS)

 end subroutine aerosol_optics



 subroutine extract_(GC, CLOCK, &
                         myState, GRID, CF, &
                         i1, i2, im,        &
                         j1, j2, jm,        &
                         km,                &
                         nymd, nhms,        &
                         cdt, rc)

    type(ESMF_GridComp), intent(INout)  :: GC           ! Grid Comp object
    type(ESMF_Clock), intent(in)        :: CLOCK        ! Clock

    type(DNA_State), pointer            :: myState      ! Legacy state
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
    type(DNA_Wrap)                      :: wrap

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
    call ESMF_UserCompGetInternalState(GC, 'DNA_State', wrap, STATUS)
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


 end module DNA_GridCompMod
