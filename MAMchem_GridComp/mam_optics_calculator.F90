#include "MAPL_Generic.h"

program mam_optics_calculator

!-----------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1   !
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: mam_optics_calculator --- Extinction calculator
!
! !INTERFACE:
!
!      Usage:  mam_optics_calculator.xx
!
! !USES:
!

  use ESMF

  use MAPL_Mod
  use MAPL_ShmemMod
  use MAPL_SimpleBundleMod

  use MAM_BaseMod
  use MAM3_DataMod, only: MAM3_MODE_NAME, MAM3_MODES 
  use MAM7_DataMod, only: MAM7_MODE_NAME, MAM7_MODES

  use MAML_OpticsTableMod
  use MAML_OpticsMod
 
  implicit none


! !DESCRIPTION: 2D/3D Aerosol Optics Calculator.
!
! !REVISION HISTORY:
!
!  27Mar2013  A. Darmenov  Initial Implementation
!
!EOP
!-----------------------------------------------------------------------

  character(len=*), parameter :: DEFAULT_CONFIG_FILE = 'mam_optics_calculator.rc'
  integer,          parameter :: MAX_STRFILE = 1024
 

  type MAM_OpticsCalculatorSetup

      type(ESMF_Config) :: config      ! private config

      logical           :: verbose     ! verbosity flag

      ! World
      type(ESMF_Grid)   :: grid        ! grid
      type(ESMF_Time)   :: time        ! time

      integer           :: im_world    ! global grid dimensions - lon
      integer           :: jm_world    ! global grid dimensions - lat
      integer           :: lm_world    ! global grid dimensions - vertical

      integer           :: Nx          ! layout 
      integer           :: Ny          !  ...

      ! MAM
      integer                                           :: scheme_id    ! MAM7 or MAM3
      character(len=ESMF_MAXSTR), pointer, dimension(:) :: mode         ! mode name

      ! Optics
      real,                       pointer, dimension(:) :: wavelength   ! wavelengths of channels/bands
      character(len=MAX_STRFILE), pointer, dimension(:) :: optics_lut   ! files with optics lookup tables

      ! Files
      character(len=MAX_STRFILE)                        :: aerosol_file ! input file with 3D aerosol number and mass mixing ratios
      character(len=MAX_STRFILE)                        :: optics_file  ! output file

  end type MAM_OpticsCalculatorSetup


                       __Iam__('mam_optics_calculator.xx')


  call main(DEFAULT_CONFIG_FILE, rc=status)  
  if (MAPL_VRFY(status, Iam, __LINE__)) call MAPL_Abort()

  call exit(status)

contains


subroutine main(config_file, rc)

  implicit none

  character(len=*),  intent( in)  :: config_file  ! config file
  integer, optional, intent(out)  :: rc           ! return code

! Local variables
! ---------------
  type(ESMF_VM)                   :: vm           ! global VM
  type(MAM_OpticsCalculatorSetup) :: setup        ! setup

  type(MAM_Scheme)                :: mam          ! MAM scheme/configuration

  type(ESMF_State)                :: aero_state   ! aerosol state

  type(MAPL_SimpleBundle)         :: q            ! aerosol mixing ratio
  type(MAPL_SimpleBundle)         :: o            ! 

  type(ESMF_Field)                :: field
  type(ESMF_FieldBundle)          :: optics       ! optics (extinction, etc.) parameters
  type(ESMF_FieldBundle)          :: bundle       ! field bundle


  type(MAML_OpticsTable)          :: lut          ! aerosol optics lookup table

  character(len=MAM_MAXSTR)       :: field_name   ! field name
  character(len=MAM_MAXSTR)       :: mode_name    ! aerosol mode name
  character(len=MAM_MAXSTR)       :: species_name ! aerosol species name
  integer                         :: n_species    ! number of aerosol species
  integer                         :: m, s         ! mode and species indexes
  integer                         :: iq           ! field index

  character(len=1024)             :: field_list   ! list of comma separated field names

  integer                         :: im, jm, km   ! local dim sizes

  integer                         :: i, j, k, n   ! loop counters 

#ifdef __PGI
  interface
     subroutine optics_compute(aero_state, rc)
     use ESMF_StateMod
     implicit none
     type(ESMF_State)     :: aero_state
     integer, intent(out) :: rc
     end subroutine
  end interface
#endif


                        __Iam__('mam_optics_calculator::main')

! Initialize the ESMF 
! -------------------
  call ESMF_Initialize(vm=vm, logKindFlag=ESMF_LOGKIND_NONE, __RC__)
  call ESMF_CalendarSetDefault(calkindflag=ESMF_CALKIND_GREGORIAN, __RC__)

! Show text banner
! ----------------
  if (MAPL_am_I_root()) then
      call text_banner()
  end if

! Initialize the setup
! --------------------
  call setup_initialize(setup, config_file, __RC__)


! Set MAM scheme
! --------------
  call MAM_SchemeInit(mam, setup%scheme_id, __RC__)


! Create a state
  aero_state = ESMF_StateCreate(name='AERO_STATE', __RC__)

! Create optics bundle and fill it with empty fields
! --------------------------------------------------
  optics = ESMF_FieldBundleCreate(name='MAM::OPTICS', __RC__)
  call optics_bundle_initialize(optics, setup%grid, mam, __RC__)

! Add the optics bundle to the aero state
! ---------------------------------------
  call ESMF_StateAdd(aero_state, (/optics/), __RC__)



  if (setup%verbose .and. .false.) then
      call ESMF_StateGet(aero_state, itemCount=n, __RC__)

      _ASSERT(n > 0,'needs informative message')

      call ESMF_StateGet(aero_state, 'MAM::OPTICS', bundle, __RC__)
      o = MAPL_SimpleBundleCreate(bundle, __RC__)

      call MAPL_SimpleBundlePrint(o)
      call MAPL_SimpleBundleDestroy(o, __RC__)
  end if


! Exercise callback mechanism
! ---------------------------
  call ESMF_MethodAdd(aero_state, label='OPTICS_COMPUTE', userRoutine=optics_compute, __RC__)

  call ESMF_MethodExecute(aero_state, label='OPTICS_COMPUTE', userRc=rc, __RC__)

  if (setup%verbose) then
      call ESMF_StateGet(aero_state, itemCount=n, __RC__)

      _ASSERT(n > 0,'needs informative message')

      call ESMF_StateGet(aero_state, 'MAM::OPTICS', bundle, __RC__)
      o = MAPL_SimpleBundleCreate(bundle, __RC__)

      call MAPL_SimpleBundlePrint(o)
      call MAPL_SimpleBundleDestroy(o, __RC__)
  end if



  loop_modes: do m = 1, mam%n_modes
  
      ! ------------------------------------------------------------------
      ! Read the 3D aerosol fields
      ! ------------------------------------------------------------------

      ! construct a list with 3D aerosol fields to read
      field_list = ''

      call MAM_AerosolModeGet(mam%mode(m), name=mode_name, n_species=n_species, __RC__)

      ! interstitial aerosol tracers 
      field_name = 'NUM_A_' // trim(mode_name)
      field_list = trim(field_list) // trim(field_name) // ','

      ! interstitial aerosol tracers
      field_name = 'WTR_A_' // trim(mode_name)
      field_list = trim(field_list) // trim(field_name) // ','

      do s = 1, n_species
          species_name = mam%mode(m)%species(s)%name
          field_name  = trim(species_name) // '_A_' // trim(mode_name)

          ! append to the list of fields to read
          field_list  = trim(field_list) // trim(field_name) // ','
      end do

      ! strip the trailing comma from the list
      field_list = trim(field_list(1:len_trim(field_list)-1))

      ! read the aerosol fields and bundle them
      q = MAPL_SimpleBundleRead(setup%aerosol_file,    &
                                setup%aerosol_file,    &
                                setup%grid,            &
                                setup%time,            &
                                only_vars=field_list,  &
                                verbose=setup%verbose, &
                                __RC__)

      if (setup%verbose) call MAPL_SimpleBundlePrint(q) 

      ! get the size of local dims
      im = size(q%coords%lons, 1)
      jm = size(q%coords%lons, 2)
      km = size(q%coords%levs)
 

      ! ------------------------------------------------------------------
      ! Read the optics lookup table
      ! ------------------------------------------------------------------
      _ASSERT(associated(setup%optics_lut),'needs informative message')
      _ASSERT(associated(setup%mode),'needs informative message')
      _ASSERT(size(setup%mode) == size(setup%optics_lut),'needs informative message')

      n = 0
      do n = 1, size(setup%mode)
          if (setup%mode(n) == mode_name) then
              exit
          end if
      end do

      _ASSERT(n > 0,'needs informative message')

      lut = MAML_OpticsTableCreate(setup%optics_lut(n), __RC__)

      call MAML_OpticsTableRead(lut, __RC__)


      ! ------------------------------------------------------------------ 
      ! Compute the aerosol optical quantities
      ! ------------------------------------------------------------------

!     call MAML_OpticsCalculator(o, q, lut, __RC__)

      call MAML_OpticsTableDestroy(lut, __RC__)

      !TODO - free all associated memory, e.g., ESMF fields, bundle, etc.
      call MAPL_SimpleBundleDestroy(q, __RC__)
  end do loop_modes

  ! ------------------------------------------------------------------
  ! Write aerosol optics quantities into a file
  ! ------------------------------------------------------------------
! call MAPL_SimpleBundleWrite(o, setup%optics_file, setup%time, __RC__)

  call MAPL_SimpleBundleDestroy(o, __RC__)


! Finalize self
! --------------
  call setup_finalize(setup, __RC__)

 
  call ESMF_StateDestroy(aero_state, __RC__)

! Finalize framework
! ------------------
  call ESMF_Finalize(__RC__)
  VERIFY_(status)
   

  RETURN_(ESMF_SUCCESS)

end subroutine main




subroutine setup_initialize(self, config_file, rc)
  implicit none

  type(MAM_OpticsCalculatorSetup), intent(inout) :: self
  character(len=*), intent(in)                   :: config_file
  integer, optional, intent(out)                 :: rc
 
                        __Iam__('mam_optics_calculator::setup_initialize')

! Local variables
! ---------------
  integer :: useShmem

! Set private config
! ------------------
  call setup_set_config_(self, config_file, __RC__)

! Set verbose flag
! ---------------- 
  call setup_set_verbosity_(self, __RC__)

! Set global grid and time
! ------------------------
  call setup_set_grid_(self, __RC__)
  call setup_set_time_(self, __RC__)

! Set MAM scheme
! --------------
  call setup_set_mam_scheme_id_(self, __RC__)

! Set modes and optics LUTs: these are tied together, i.e., 
! a hash table with keys modes and values the optics LUTs
! ---------------------------------------------------------
  call setup_set_mam_modes_(self, __RC__)
  call setup_set_mam_optics_lut_(self, __RC__)

! Set wavelengths of channel/bands
! --------------------------------
  call setup_set_wavelengths_(self, __RC__)

! Set I/O files
! -------------
  call setup_set_io_files_(self, __RC__)

! Check if user wants to use node shared memory (default is no)
! -------------------------------------------------------------
  call ESMF_ConfigGetAttribute(self%config, useShmem, label='USE_SHMEM:', default=0, __RC__)

  if (useShmem /= 0) then
      call MAPL_InitializeShmem(__RC__)
  end if

  RETURN_(ESMF_SUCCESS)
end subroutine setup_initialize


subroutine setup_finalize(self, rc)
  implicit none

  type(MAM_OpticsCalculatorSetup), intent(inout) :: self
  integer, optional, intent(out)                 :: rc
 
                        __Iam__('mam_optics_calculator::setup_finalize')


  if (associated(self%mode))       deallocate(self%mode)
  if (associated(self%optics_lut)) deallocate(self%optics_lut)
  if (associated(self%wavelength)) deallocate(self%wavelength)
   
  call ESMF_ConfigDestroy(self%config, __RC__)

  call MAPL_FinalizeShmem(__RC__)

  RETURN_(ESMF_SUCCESS)
end subroutine setup_finalize


subroutine optics_bundle_initialize(optics, grid, mam, rc)
  implicit none

  type(ESMF_FieldBundle)         :: optics    ! optics (extinction, etc.) parameters
  type(ESMF_Grid)                :: grid      ! grid
  type(MAM_Scheme)               :: mam       ! MAM scheme/configuration
  integer, optional, intent(out) :: rc        ! return code


                __Iam__('mam_optics_calculator::optics_bundle_initialize')

  ! local
  type(ESMF_Field)               :: field
  character(len=MAM_MAXSTR)      :: field_name   ! field name
  character(len=MAM_MAXSTR)      :: mode_name    ! aerosol mode name
  integer                        :: m            ! mode index


  ! Fill in optics bundle with empty fields
  ! ---------------------------------------

  ! field for total quantity, e.g., extinction
  field_name = 'ext'
  field_name = ESMF_UtilStringLowerCase(field_name, __RC__)
  field = MAPL_FieldCreateEmpty(trim(field_name), grid, __RC__)
  call MAPL_FieldAllocCommit(field, dims=MAPL_DimsHorzVert, location=MAPL_VLocationCenter, typekind=MAPL_R4, hw=0, __RC__)
  call MAPL_FieldBundleAdd(optics, field, __RC__)

  ! fields for mode quantities
  do m = 1, mam%n_modes
      call MAM_AerosolModeGet(mam%mode(m), name=mode_name, __RC__)

      ! create an empty field
      field_name = 'ext_' // trim(mode_name)
      field_name = ESMF_UtilStringLowerCase(field_name, __RC__)
      field = MAPL_FieldCreateEmpty(trim(field_name), grid, __RC__)
      call MAPL_FieldAllocCommit(field, dims=MAPL_DimsHorzVert, location=MAPL_VLocationCenter, typekind=MAPL_R4, hw=0, __RC__)
      call MAPL_FieldBundleAdd(optics, field, __RC__)
  end do

  RETURN_(ESMF_SUCCESS)

end subroutine optics_bundle_initialize


subroutine optics_compute(aero_state, rc)
  use ESMF_StateMod

  implicit none

  type(ESMF_State)     :: aero_state ! aerosol state
  integer, intent(out) :: rc         ! return code

          __Iam__('mam_optics_calculator::optics_bundle_initialize')

  
  !local
  type(ESMF_FieldBundle)          :: optics
  type(ESMF_Field)                :: field

  character(len=ESMF_MAXSTR)      :: name

  real, pointer, dimension(:,:,:) :: q

  
  call ESMF_StateGet(aero_state, name=name, __RC__)

  call ESMF_StateGet(aero_state, 'MAM::OPTICS', optics, __RC__)

  call ESMFL_BundleGetPointerToData(optics, 'ext', q, __RC__)
  
  q = q + 5.0 

  RETURN_(ESMF_SUCCESS) 

end subroutine optics_compute



subroutine setup_set_config_(self, config_file, rc)
  implicit none
 
  type(MAM_OpticsCalculatorSetup), intent(inout) :: self
  character(len=*), intent(in)                   :: config_file
  integer, optional, intent(out)                 :: rc
 
                        __Iam__('mam_optics_calculator::setup_set_config_')


  self%config = ESMF_ConfigCreate(__RC__)
  call ESMF_ConfigLoadFile(self%config, fileName=trim(config_file), __RC__)
 
  RETURN_(ESMF_SUCCESS)
end subroutine setup_set_config_


subroutine setup_set_verbosity_(self, rc)
  implicit none
 
  type(MAM_OpticsCalculatorSetup), intent(inout) :: self
  integer, optional, intent(out)                 :: rc
 
                        __Iam__('mam_optics_calculator::setup_set_verbosity_')

  call ESMF_ConfigGetAttribute(self%config, self%verbose, label='verbose:', __RC__)
 
  RETURN_(ESMF_SUCCESS)
end subroutine setup_set_verbosity_


subroutine setup_set_io_files_(self, rc)
  implicit none
 
  type(MAM_OpticsCalculatorSetup), intent(inout) :: self
  integer, optional, intent(out)                 :: rc
 
                        __Iam__('mam_optics_calculator::setup_set_io_files_')

  call ESMF_ConfigGetAttribute(self%config, self%aerosol_file, label='aerosol_file:', __RC__)
  call ESMF_ConfigGetAttribute(self%config, self%optics_file,  label='optics_file:',  __RC__)
 
  RETURN_(ESMF_SUCCESS)
end subroutine setup_set_io_files_


subroutine setup_set_grid_(self, rc)
  implicit none
 
  type(MAM_OpticsCalculatorSetup), intent(inout) :: self
  integer, optional, intent(out)                 :: rc
 
                        __Iam__('mam_optics_calculator::setup_set_grid_')


! World grid dimensions and layout
! --------------------------------
  call ESMF_ConfigGetAttribute(self%config, self%im_world, label='WORLD_IM:', __RC__)
  call ESMF_ConfigGetAttribute(self%config, self%jm_world, label='WORLD_JM:', __RC__)
  call ESMF_ConfigGetAttribute(self%config, self%lm_world, label='WORLD_LM:', __RC__)
  call ESMF_ConfigGetAttribute(self%config, self%Nx,       label='NX:',       __RC__)
  call ESMF_ConfigGetAttribute(self%config, self%Ny,       label='NY:',       __RC__)

! Create global grid
! ------------------
  self%grid = MAPL_LatLonGridCreate(name     = 'etaGrid',     &
                                    Nx       = self%Nx,       &
                                    Ny       = self%Ny,       &
                                    IM_World = self%im_world, &
                                    JM_World = self%jm_world, &
                                    LM_World = self%lm_world, &
                                    __RC__)

! Validate grid
! -------------
  call ESMF_GridValidate(self%grid, __RC__)


  RETURN_(ESMF_SUCCESS)
end subroutine setup_set_grid_


subroutine setup_set_time_(self, rc)
  implicit none
 
  type(MAM_OpticsCalculatorSetup), intent(inout) :: self
  integer, optional, intent(out)                 :: rc
 
                        __Iam__('mam_optics_calculator::setup_set_time_')

! Local variables
! ---------------
  integer :: year, month, day
  integer :: hours, minutes, seconds
  integer :: nymd, nhms

! Initialize the date/time
! ------------------------
  nymd = 0
  nhms = 0

! Get date/time from config
! -------------------------
  call ESMF_ConfigGetAttribute(self%config, nymd, label='nymd:', __RC__)
  call ESMF_ConfigGetAttribute(self%config, nhms, label='nhms:', __RC__)

! Set ESMF Time
! -------------
  year  = nymd / 10000; month   = (nymd - 10000 *  year) / 100; day     = nymd - (10000 * year   + 100 * month)
  hours = nhms / 10000; minutes = (nhms - 10000 * hours) / 100; seconds = nhms - (10000 * hours  + 100 * minutes)

  call ESMF_TimeSet(self%time, yy=year, mm=month, dd=day, h=hours, m=minutes, s=seconds, __RC__)

  RETURN_(ESMF_SUCCESS)
end subroutine setup_set_time_


subroutine setup_set_mam_scheme_id_(self, rc)
  implicit none

  type(MAM_OpticsCalculatorSetup), intent(inout) :: self
  integer, optional, intent(out)                 :: rc

                        __Iam__('mam_optics_calculator::setup_set_mam_scheme_id_')

! Local variables
! ---------------
  character(len=ESMF_MAXSTR) :: scheme      ! name of MAM scheme/configuration


  call ESMF_ConfigGetAttribute(self%config, scheme, label='scheme:', default='MAM7',  __RC__)
  
  scheme = ESMF_UtilStringUpperCase(scheme, __RC__)

  select case (scheme)
      case ('MAM7')
      self%scheme_id = MAM7_SCHEME

      case default
      __raise__(MAM_UNKNOWN_SCHEME_ERROR, "Unsupported MAM scheme: " // trim(scheme))
  end select

  RETURN_(ESMF_SUCCESS)
end subroutine setup_set_mam_scheme_id_


subroutine setup_set_mam_modes_(self, rc)
  implicit none

  type(MAM_OpticsCalculatorSetup), intent(inout) :: self
  integer, optional, intent(out)                 :: rc

                        __Iam__('mam_optics_calculator::setup_set_mam_modes_')

! Local variables
! ---------------
  character(len=ESMF_MAXSTR) :: mode_
  logical :: flag
  integer :: i, n


  call ESMF_ConfigFindLabel(self%config, 'modes:', __RC__)

  n = 0
  do while (.true.)
      call ESMF_ConfigGetAttribute(self%config, mode_, rc=status)
      if (status == ESMF_SUCCESS) then
          n = n + 1
      else
          exit
      end if
  end do

  _ASSERT(n > 0,'needs informative message')
  _ASSERT(n < (MAM_MAX_NUMBER_MODES + 1),'needs informative message')
  
  allocate(self%mode(n), __STAT__)
  
  call ESMF_ConfigFindLabel(self%config, 'modes:', __RC__)
  do i = 1, n
      call ESMF_ConfigGetAttribute(self%config, self%mode(i), __RC__)
  end do

  select case (self%scheme_id)
      case (MAM7_SCHEME)
      do i = 1, n 
          flag = any(MAM7_MODE_NAME(:) .eq. self%mode(i))
          if (flag .eqv. .false.) then
              mode_ = self%mode(i)
              deallocate(self%mode, __STAT__)
              __raise__(MAM_UNKNOWN_AEROSOL_MODE_ERROR, "Unknown MAM mode " // trim(mode_))
          end if   
      end do    

      case (MAM3_SCHEME)
      do i = 1, n 
          flag = any(MAM3_MODE_NAME(:) .eq. self%mode(i))
          if (flag .eqv. .false.) then
              mode_ = self%mode(i)
              deallocate(self%mode, __STAT__)
              __raise__(MAM_UNKNOWN_AEROSOL_MODE_ERROR, "Unknown MAM mode " // trim(mode_))
          end if    
      end do

      case default
      __raise__(MAM_UNKNOWN_SCHEME_ERROR, "Unsupported MAM scheme")
  end select

  RETURN_(ESMF_SUCCESS)
end subroutine setup_set_mam_modes_


subroutine setup_set_mam_optics_lut_(self, rc)
  implicit none

  type(MAM_OpticsCalculatorSetup), intent(inout) :: self
  integer, optional, intent(out)                 :: rc

                        __Iam__('mam_optics_calculator::setup_set_mam_optics_lut_')

! Local variables
! ---------------
  character(len=MAX_STRFILE) :: optics_lut_
  character(len=ESMF_MAXSTR) :: optics_lut_label
  integer :: i, n

  if (associated(self%mode)) then
      n = size(self%mode)
  else
      n = 0
  end if

  _ASSERT(n > 0,'needs informative message')
  _ASSERT(n < (MAM_MAX_NUMBER_MODES + 1),'needs informative message')

  allocate(self%optics_lut(n), __STAT__)

  do i = 1, n
      optics_lut_label = 'MAM_' // trim(self%mode(i)) // '_OPTICS:'
      call ESMF_ConfigGetAttribute(self%config, self%optics_lut(i), label=trim(optics_lut_label), __RC__)
  end do

  RETURN_(ESMF_SUCCESS)
end subroutine setup_set_mam_optics_lut_


subroutine setup_set_wavelengths_(self, rc)
  implicit none

  type(MAM_OpticsCalculatorSetup), intent(inout) :: self
  integer, optional, intent(out)                 :: rc

                        __Iam__('mam_optics_calculator::setup_set_wavelengths_')

! Local variables
! ---------------
  real    :: wavelength_
  integer :: i, n

  call ESMF_ConfigFindLabel(self%config, 'wavelength:', __RC__)

  n = 0
  do while (.true.)
      call ESMF_ConfigGetAttribute(self%config, wavelength_, rc=status)
      if (status == ESMF_SUCCESS) then
          n = n + 1
      else
          exit
      end if
  end do


  _ASSERT(n > 0,'needs informative message')

  allocate(self%wavelength(n), __STAT__)
  
  call ESMF_ConfigFindLabel(self%config, 'wavelength:', __RC__)
  do i = 1, n
      call ESMF_ConfigGetAttribute(self%config, self%wavelength(i), __RC__)
  end do

  RETURN_(ESMF_SUCCESS)
end subroutine setup_set_wavelengths_


subroutine text_banner()
  implicit none

  print *
  print *, '     --------------------------------------'
  print *, '         MAM - 3D Extinction Calculator    '
  print *, '     --------------------------------------'
  print *

end subroutine text_banner



end program mam_optics_calculator

