#include "MAPL_Generic.h"

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 910.1     !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: MATRIXchem_GridCompMod - Implements MATRIX Chemistry
!
! !INTERFACE:
!
module MATRIXchem_GridCompMod
!
! !USES:
!
   use ESMF
   use MAPL_Mod

   use Chem_UtilMod,       only: Chem_UtilResVal

   use SeasaltEmissionMod, only: SeasaltEmission

   use DustEmissionMod,    only: MAM_DustEmissionGOCART, MAM_DustEmission

   use aero_config,        only: MATRIX_CONFIGURATION       => MECH,          &
                                 MATRIX_N_AEROSOL_MODES     => NMODES,        &
                                 MATRIX_N_AEROSOLS          => NAEROBOX

   use aero_param,         only: MATRIX_N_GASES             => NGASES,        &
                                 MATRIX_N_EMIS_SPECIES      => NEMIS_SPCS,    &
                                 MATRIX_N_AEROSOL_DIAG      => NDIAG_AERO,    &
                                 MATRIX_N_MASS_SPECIES      => NMASS_SPCS,    & 
                                 MATRIX_I                   => IXXX,          &
                                 MATRIX_J                   => IYYY,          &
                                 MATRIX_L                   => ILAY,          &
                                 MATRIX_MW_H2SO4            => MW_H2SO4,      &
                                 MATRIX_MW_NH3              => MW_NH3,        &
                                 MATRIX_MASS_NO3            => MASS_NO3,      &
                                 MATRIX_MASS_NH4            => MASS_NH4,      &
                                 MATRIX_MASS_H2O            => MASS_H2O,      &  
                                 MATRIX_NUMB_AKK            => NUMB_AKK_1,    & 
                                 MATRIX_MASS_AKK_SU         => MASS_AKK_SULF, &  
                                 MATRIX_NUMB_ACC            => NUMB_ACC_1,    &
                                 MATRIX_MASS_ACC_SU         => MASS_ACC_SULF, &
                                 MATRIX_NUMB_DD1            => NUMB_DD1_1,    &
                                 MATRIX_MASS_DD1_SU         => MASS_DD1_SULF, & 
                                 MATRIX_MASS_DD1_DU         => MASS_DD1_DUST, & 
                                 MATRIX_NUMB_DS1            => NUMB_DS1_1,    &
                                 MATRIX_MASS_DS1_SU         => MASS_DS1_SULF, & 
                                 MATRIX_MASS_DS1_DU         => MASS_DS1_DUST, & 
                                 MATRIX_NUMB_DD2            => NUMB_DD2_1,    & 
                                 MATRIX_MASS_DD2_SU         => MASS_DD2_SULF, & 
                                 MATRIX_MASS_DD2_DU         => MASS_DD2_DUST, & 
                                 MATRIX_NUMB_DS2            => NUMB_DS2_1,    & 
                                 MATRIX_MASS_DS2_SU         => MASS_DS2_SULF, & 
                                 MATRIX_MASS_DS2_DU         => MASS_DS2_DUST, & 
                                 MATRIX_NUMB_SSA            => NUMB_SSA_1,    & 
                                 MATRIX_MASS_SSA_SU         => MASS_SSA_SULF, & 
                                 MATRIX_MASS_SSA_SS         => MASS_SSA_SEAS, & 
                                 MATRIX_NUMB_SSC            => NUMB_SSC_1,    &
                                 MATRIX_MASS_SSC_SU         => MASS_SSC_SULF, & 
                                 MATRIX_MASS_SSC_SS         => MASS_SSC_SEAS, &
                                 MATRIX_NUMB_OCC            => NUMB_OCC_1,    &
                                 MATRIX_MASS_OCC_SU         => MASS_OCC_SULF, & 
                                 MATRIX_MASS_OCC_OC         => MASS_OCC_OCAR, &
                                 MATRIX_NUMB_BC1            => NUMB_BC1_1,    &
                                 MATRIX_MASS_BC1_SU         => MASS_BC1_SULF, & 
                                 MATRIX_MASS_BC1_BC         => MASS_BC1_BCAR, &
                                 MATRIX_NUMB_BC2            => NUMB_BC2_1,    &
                                 MATRIX_MASS_BC2_SU         => MASS_BC2_SULF, & 
                                 MATRIX_MASS_BC2_BC         => MASS_BC2_BCAR, &
                                 MATRIX_NUMB_BC3            => NUMB_BC3_1,    &
                                 MATRIX_MASS_BC3_SU         => MASS_BC3_SULF, & 
                                 MATRIX_MASS_BC3_BC         => MASS_BC3_BCAR, &
                                 MATRIX_NUMB_DBC            => NUMB_DBC_1,    &
                                 MATRIX_MASS_DBC_SU         => MASS_DBC_SULF, & 
                                 MATRIX_MASS_DBC_BC         => MASS_DBC_BCAR, & 
                                 MATRIX_MASS_DBC_DU         => MASS_DBC_DUST, &
                                 MATRIX_NUMB_BOC            => NUMB_BOC_1,    &
                                 MATRIX_MASS_BOC_SU         => MASS_BOC_SULF, & 
                                 MATRIX_MASS_BOC_BC         => MASS_BOC_BCAR, & 
                                 MATRIX_MASS_BOC_OC         => MASS_BOC_OCAR, &
                                 MATRIX_NUMB_BCS            => NUMB_BCS_1,    &
                                 MATRIX_MASS_BCS_SU         => MASS_BCS_SULF, &
                                 MATRIX_MASS_BCS_BC         => MASS_BCS_BCAR, &
                                 MATRIX_NUMB_MXX            => NUMB_MXX_1,    &
                                 MATRIX_MASS_MXX_SU         => MASS_MXX_SULF, & 
                                 MATRIX_MASS_MXX_BC         => MASS_MXX_BCAR, & 
                                 MATRIX_MASS_MXX_OC         => MASS_MXX_OCAR, & 
                                 MATRIX_MASS_MXX_DU         => MASS_MXX_DUST, & 
                                 MATRIX_MASS_MXX_SS         => MASS_MXX_SEAS

   use aero_setup,         only: MATRIX_SU_MAP              => SULF_MAP,     & 
                                 MATRIX_DU_MAP              => DUST_MAP,     &
                                 MATRIX_SS_MAP              => SEAS_MAP,     &
                                 MATRIX_OC_MAP              => OCAR_MAP,     &
                                 MATRIX_BC_MAP              => BCAR_MAP,     &
                                 MATRIX_setup_config        => SETUP_CONFIG, &
                                 MATRIX_setup_species_maps  => SETUP_SPECIES_MAPS, &
                                 MATRIX_setup_dp0           => SETUP_DP0,    &
                                 MATRIX_setup_aero_mass_map => SETUP_AERO_MASS_MAP, &
                                 MATRIX_setup_coag_tensors  => SETUP_COAG_TENSORS, &
                                 MATRIX_setup_emis          => SETUP_EMIS,   &
                                 MATRIX_setup_kci           => SETUP_KCI

   use aero_coag,          only: MATRIX_setup_kij           => SETUP_KIJ

   use aero_npf,           only: MATRIX_setup_npfmass       => SETUP_NPFMASS

   use amp_aerosol,        only: MATRIX_VDDEP_AERO          => VDDEP_AERO, &
                                 MATRIX_DIAMETER            => DIAM,       &
                                 MATRIX_NACTV               => NACTV,      &
                                 MATRIX_CCNSS               => CCNSS


   use constant,           only: MATRIX_MW_AIR              => mair

   use aero_diam,          only: MATRIX_setup_diam          => SETUP_DIAM

   implicit none
   private
!
! !PUBLIC MEMBER FUNCTIONS:

   public SetServices
!
! !DESCRIPTION: 
!
!  {\tt MATRIXchem\_GridComp} is an ESMF gridded component implementing
!  the MATRIX aerosol microphysical processes.
!
!  Developed for GEOS-5 release Fortuna 2.0 and later.
!
! !REVISION HISTORY:
!
!  13 Mar 2015    A.Darmenov  Coupled the MATRIX module with GEOS-5.
!  06 Dec 2009    da Silva    Created the MATRIX skeleton.
!
!EOP
!-------------------------------------------------------------------------

!  Legacy state
!  ------------
   type MATRIX_State
       private

       integer                     :: configuration = 1 ! MATRIX configuration

       type(ESMF_Config)           :: CF                ! Private Config
       type(ESMF_Grid)             :: grid              ! Grid

       integer                     :: im_world          ! Horizontal dimensions - lon
       integer                     :: jm_world          ! Horizontal dimensions - lat

       real                        :: dt                ! Model time step

       real                        :: f_emiss_seasalt   ! Global seasalt emissions tuning parameter
       real                        :: f_emiss_dust      ! Global dust    emissions tuning parameter

       logical                     :: verbose           ! verbosity flag
   end type MATRIX_State

!  Hook for the ESMF
!  -----------------
   type MATRIX_Wrap
       type (MATRIX_State), pointer :: PTR => null()
   end type MATRIX_Wrap

contains


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 910.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: SetServices --- Sets IRF services for the MATRIXchem Grid Component
!
! !INTERFACE:

   subroutine SetServices(GC, rc)

! !ARGUMENTS:

    type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
    integer, optional                  :: rc  ! return code

! !DESCRIPTION: Sets Initialize, Run and Finalize services. 
!
! !REVISION HISTORY:
!
!  13 Mar 2015    A.Darmenov  Coupled the MATRIX module with GEOS-5.
!  01 Dec 2009    da Silva    First crack.
!
!EOP
!-------------------------------------------------------------------------

                            __Iam__('SetServices')

!   Local derived type aliases
!   --------------------------
    type (MATRIX_State), pointer  :: self   ! internal private, that is
    type (MATRIX_wrap)            :: wrap

    character(len=ESMF_MAXSTR)    :: comp_name

!                              ------------

!   Get my name and set-up traceback handle
!   ---------------------------------------
    call ESMF_GridCompGet(GC, name=comp_name, __RC__)
    Iam = trim(comp_name) // '::' // trim(Iam)


!   Wrap internal state for storing in GC; rename legacyState
!   -------------------------------------
    allocate(self, __STAT__)
    wrap%ptr => self


!   Load private Config Attributes
!   ------------------------------
    self%CF = ESMF_ConfigCreate(__RC__)
    call ESMF_ConfigLoadFile(self%CF, 'MATRIXchem_GridComp.rc', __RC__)

    call ESMF_ConfigGetAttribute(self%CF, self%verbose,       Label='verbose:', default=.false.,  __RC__)

    call ESMF_ConfigGetAttribute(self%CF, self%configuration, Label='matrix:',  default=1,        __RC__)
    ASSERT_(self%configuration == MATRIX_CONFIGURATION)



!   Set the profiling timers
!   ------------------------
    call MAPL_TimerAdd(GC, name='RUN',                         __RC__)
#ifdef __MATRIX_TODO__
    call MAPL_TimerAdd(GC, name='-EMISSIONS',                  __RC__)
    call MAPL_TimerAdd(GC, name='-MICROPHYSICS',               __RC__)
    call MAPL_TimerAdd(GC, name='-REMOVAL',                    __RC__)
    call MAPL_TimerAdd(GC, name='--REMOVAL_DRY',               __RC__)
    call MAPL_TimerAdd(GC, name='--REMOVAL_WET',               __RC__)
    call MAPL_TimerAdd(GC, name='-DIAGNOSTICS',                __RC__)
    call MAPL_TimerAdd(GC, name='--DIAGNOSTICS_SEASALT',       __RC__)
    call MAPL_TimerAdd(GC, name='--DIAGNOSTICS_DUST',          __RC__)
    call MAPL_TimerAdd(GC, name='--DIAGNOSTICS_CIM',           __RC__)
    call MAPL_TimerAdd(GC, name='--DIAGNOSTICS_AOT',           __RC__)
#endif

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
    call ESMF_UserCompSetInternalState(GC, 'MATRIX_state', wrap, STATUS)
    VERIFY_(STATUS)
  
!                         ------------------
!                         MAPL Data Services
!                         ------------------

!BOP
!
! !IMPORT STATE:

#include "MATRIXchem_ImportSpec___.h"

! !INTERNAL STATE:

#include "MATRIXchem_InternalSpec___.h"

! !EXTERNAL STATE:

#include "MATRIXchem_ExportSpec___.h"

!   Generic Set Services
!   --------------------
    call MAPL_GenericSetServices(GC, __RC__)

!   Anounce that MATRIX is active
!   -----------------------------
    if (MAPL_AM_I_ROOT()) then
        write (*,*) trim(Iam)//': ACTIVE'
        write (*,*)
    end if

!   All done
!   --------

    RETURN_(ESMF_SUCCESS)

  end subroutine SetServices


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Initialize_ --- Initialize MATRIXchem
!
! !INTERFACE:
!

   subroutine Initialize_(GC, IMPORT, EXPORT, clock, rc)

! !USES:

    implicit none

! !INPUT PARAMETERS:

    type(ESMF_Clock),    intent(inout)  :: clock  ! The clock

! !OUTPUT PARAMETERS:

    type(ESMF_GridComp), intent(inout)  :: GC     ! Grid Component
    type(ESMF_State),    intent(inout)  :: IMPORT ! Import State
    type(ESMF_State),    intent(inout)  :: EXPORT ! Export State
    integer,             intent(out)    :: rc     ! Error return code:
                                                  !  0 - all is well
                                                  !  1 - 

! !DESCRIPTION: This is a simple ESMF wrapper.
!
! !REVISION HISTORY:
!
!  13 Mar 2015    A.Darmenov  Coupled the MATRIX module with GEOS-5.
!  01 Dec 2009    da Silva    First crack.
!
!EOP
!-------------------------------------------------------------------------

                              __Iam__('Initialize_')

    type(MATRIX_state), pointer   :: self               ! Legacy state
    type(ESMF_Grid)               :: GRID               ! Grid
    type(ESMF_Config)             :: CF                 ! Universal Config 

    integer                       :: im_World, jm_World ! Global 2D Dimensions
    integer                       :: im, jm, lm         ! 3D Dimensions
    real(ESMF_KIND_R4), pointer   :: lons(:,:)          ! Longitudes
    real(ESMF_KIND_R4), pointer   :: lats(:,:)          ! Latitudes

    integer                       :: nymd, nhms         ! date, time
    real                          :: cdt                ! time step in secs

    character(len=ESMF_MAXSTR)    :: comp_name          ! name of the component

    integer, parameter            :: n_res = 6          ! number of horizontal resolutions (a, b, c, d, e)
    real, dimension(n_res)        :: f_res              ! buffer for the resolution dependent factors
    integer                       :: n                  ! counter

    integer, allocatable          :: matrix_aerosol_indexes(:)
   

!   Declare pointers to IMPORT/EXPORT/INTERNAL states 
!   -------------------------------------------------
#include "MATRIXchem_DeclarePointer___.h"
  
!   Get my name and set-up traceback handle
!   ---------------------------------------
    call ESMF_GridCompGet(GC, name=comp_name, __RC__)
    Iam = trim(comp_name) // trim(Iam)

!                               --------
    if (MAPL_AM_I_ROOT()) then
        write (*,*) trim(Iam)//': Starting...'
        write (*,*)
    end if
  

!   Initialize MAPL Generic
!   -----------------------
    call MAPL_GenericInitialize (GC, IMPORT, EXPORT, clock,  __RC__)

!   Get pointers to IMPORT/EXPORT/INTERNAL states 
!   ---------------------------------------------
#include "MATRIXchem_GetPointer___.h"

!   Extract relevant runtime information
!   ------------------------------------
    call extract_(GC, clock, self, GRID, CF, &
                  im_World, jm_World,        &
                  im, jm, lm, lons, lats,    &
                  nymd, nhms, cdt, __RC__)

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
    call ESMF_ConfigFindLabel(self%CF, 'f_emissions_seasalt:', __RC__)
    do n = 1, n_res
        call ESMF_ConfigGetAttribute(self%CF, f_res(n), __RC__)
    end do
    self%f_emiss_seasalt = Chem_UtilResVal(self%im_world, self%jm_world, f_res(:), __RC__)

    call ESMF_ConfigFindLabel(self%CF, 'f_emissions_dust:', __RC__)
    do n = 1, n_res
        call ESMF_ConfigGetAttribute(self%CF, f_res(n), __RC__)
    end do
    self%f_emiss_dust = Chem_UtilResVal(self%im_world, self%jm_world, f_res(:), __RC__)


!   MATRIX core
!   -----------
    call MATRIX_setup_config
    call MATRIX_setup_species_maps
    call MATRIX_setup_dp0
    call MATRIX_setup_aero_mass_map
    call MATRIX_setup_coag_tensors 
    call MATRIX_setup_dp0
    call MATRIX_setup_kij
    call MATRIX_setup_emis
    call MATRIX_setup_kci
    call MATRIX_setup_npfmass
!   call MATRIX_setup_diam

#ifdef  __MATRIX_TODO__
!     CALL SETUP_RAD
#endif


#ifdef  __MATRIX_TODO__
    1. 'FIRSTIME' code goes here and gets dissabled in the MATRIX core.
    2. ...might need to add fields to the private internal state.
#endif

    if (MAPL_AM_I_ROOT() .and. self%verbose) then
        write (*,*) 'MODES    : ', MATRIX_N_AEROSOL_MODES
        write (*,*) 'AEROSOLS : ', MATRIX_N_AEROSOLS
        write (*,*) 'MAP SU   : ', MATRIX_SU_MAP
        write (*,*) 'MAP OC   : ', MATRIX_OC_MAP
        write (*,*) 'MAP BC   : ', MATRIX_BC_MAP
        write (*,*) 'MAP SS   : ', MATRIX_SS_MAP
        write (*,*) 'MAP DU   : ', MATRIX_DU_MAP
        write (*,*)
        write (*,*) 'MATRIX_MASS_NO3    : ', MATRIX_MASS_NO3
        write (*,*) 'MATRIX_MASS_NH4    : ', MATRIX_MASS_NH4
        write (*,*) 'MATRIX_MASS_H2O    : ', MATRIX_MASS_H2O
        write (*,*) 'MATRIX_NUMB_AKK    : ', MATRIX_NUMB_AKK
        write (*,*) 'MATRIX_MASS_AKK_SU : ', MATRIX_MASS_AKK_SU
        write (*,*) 'MATRIX_NUMB_ACC    : ', MATRIX_NUMB_ACC
        write (*,*) 'MATRIX_MASS_ACC_SU : ', MATRIX_MASS_ACC_SU
        write (*,*) 'MATRIX_NUMB_DD1    : ', MATRIX_NUMB_DD1
        write (*,*) 'MATRIX_MASS_DD1_SU : ', MATRIX_MASS_DD1_SU
        write (*,*) 'MATRIX_MASS_DD1_DU : ', MATRIX_MASS_DD1_DU
        write (*,*) 'MATRIX_NUMB_DS1    : ', MATRIX_NUMB_DS1
        write (*,*) 'MATRIX_MASS_DS1_SU : ', MATRIX_MASS_DS1_SU 
        write (*,*) 'MATRIX_MASS_DS1_DU : ', MATRIX_MASS_DS1_DU
        write (*,*) 'MATRIX_NUMB_DD2    : ', MATRIX_NUMB_DD2
        write (*,*) 'MATRIX_MASS_DD2_SU : ', MATRIX_MASS_DD2_SU
        write (*,*) 'MATRIX_MASS_DD2_DU : ', MATRIX_MASS_DD2_DU 
        write (*,*) 'MATRIX_NUMB_DS2    : ', MATRIX_NUMB_DS2
        write (*,*) 'MATRIX_MASS_DS2_SU : ', MATRIX_MASS_DS2_SU
        write (*,*) 'MATRIX_MASS_DS2_DU : ', MATRIX_MASS_DS2_DU
        write (*,*) 'MATRIX_NUMB_SSA    : ', MATRIX_NUMB_SSA
        write (*,*) 'MATRIX_MASS_SSA_SU : ', MATRIX_MASS_SSA_SU
        write (*,*) 'MATRIX_MASS_SSA_SS : ', MATRIX_MASS_SSA_SS
        write (*,*) 'MATRIX_NUMB_SSC    : ', MATRIX_NUMB_SSC
        write (*,*) 'MATRIX_MASS_SSC_SU : ', MATRIX_MASS_SSC_SU
        write (*,*) 'MATRIX_MASS_SSC_SS : ', MATRIX_MASS_SSC_SS
        write (*,*) 'MATRIX_NUMB_OCC    : ', MATRIX_NUMB_OCC
        write (*,*) 'MATRIX_MASS_OCC_SU : ', MATRIX_MASS_OCC_SU
        write (*,*) 'MATRIX_MASS_OCC_OC : ', MATRIX_MASS_OCC_OC
        write (*,*) 'MATRIX_NUMB_BC1    : ', MATRIX_NUMB_BC1
        write (*,*) 'MATRIX_MASS_BC1_SU : ', MATRIX_MASS_BC1_SU
        write (*,*) 'MATRIX_MASS_BC1_BC : ', MATRIX_MASS_BC1_BC
        write (*,*) 'MATRIX_NUMB_BC2    : ', MATRIX_NUMB_BC2
        write (*,*) 'MATRIX_MASS_BC2_SU : ', MATRIX_MASS_BC2_SU
        write (*,*) 'MATRIX_MASS_BC2_BC : ', MATRIX_MASS_BC2_BC
        write (*,*) 'MATRIX_NUMB_BC3    : ', MATRIX_NUMB_BC3
        write (*,*) 'MATRIX_MASS_BC3_SU : ', MATRIX_MASS_BC3_SU
        write (*,*) 'MATRIX_MASS_BC3_BC : ', MATRIX_MASS_BC3_BC
        write (*,*) 'MATRIX_NUMB_DBC    : ', MATRIX_NUMB_DBC
        write (*,*) 'MATRIX_MASS_DBC_SU : ', MATRIX_MASS_DBC_SU
        write (*,*) 'MATRIX_MASS_DBC_BC : ', MATRIX_MASS_DBC_BC
        write (*,*) 'MATRIX_MASS_DBC_DU : ', MATRIX_MASS_DBC_DU
        write (*,*) 'MATRIX_NUMB_BOC    : ', MATRIX_NUMB_BOC
        write (*,*) 'MATRIX_MASS_BOC_SU : ', MATRIX_MASS_BOC_SU
        write (*,*) 'MATRIX_MASS_BOC_BC : ', MATRIX_MASS_BOC_BC
        write (*,*) 'MATRIX_MASS_BOC_OC : ', MATRIX_MASS_BOC_OC
        write (*,*) 'MATRIX_NUMB_BCS    : ', MATRIX_NUMB_BCS
        write (*,*) 'MATRIX_MASS_BCS_SU : ', MATRIX_MASS_BCS_SU
        write (*,*) 'MATRIX_MASS_BCS_BC : ', MATRIX_MASS_BCS_BC
        write (*,*) 'MATRIX_NUMB_MXX    : ', MATRIX_NUMB_MXX
        write (*,*) 'MATRIX_MASS_MXX_SU : ', MATRIX_MASS_MXX_SU
        write (*,*) 'MATRIX_MASS_MXX_BC : ', MATRIX_MASS_MXX_BC
        write (*,*) 'MATRIX_MASS_MXX_OC : ', MATRIX_MASS_MXX_OC
        write (*,*) 'MATRIX_MASS_MXX_DU : ', MATRIX_MASS_MXX_DU
        write (*,*) 'MATRIX_MASS_MXX_SS : ', MATRIX_MASS_MXX_SS
    end if

    allocate(matrix_aerosol_indexes(MATRIX_N_AEROSOLS), __STAT__)
    matrix_aerosol_indexes = 0

    ! verify MATRIX::AEROSOL array indexes
    matrix_aerosol_indexes = (/MATRIX_MASS_NO3, MATRIX_MASS_NH4,    MATRIX_MASS_H2O,                        &   
                               MATRIX_NUMB_AKK, MATRIX_MASS_AKK_SU,                                         &
                               MATRIX_NUMB_ACC, MATRIX_MASS_ACC_SU,                                         &
                               MATRIX_NUMB_DD1, MATRIX_MASS_DD1_SU, MATRIX_MASS_DD1_DU,                     &
                               MATRIX_NUMB_DS1, MATRIX_MASS_DS1_SU, MATRIX_MASS_DS1_DU,                     &
                               MATRIX_NUMB_DD2, MATRIX_MASS_DD2_SU, MATRIX_MASS_DD2_DU,                     &
                               MATRIX_NUMB_DS2, MATRIX_MASS_DS2_SU, MATRIX_MASS_DS2_DU,                     &
                               MATRIX_NUMB_SSA, MATRIX_MASS_SSA_SU, MATRIX_MASS_SSA_SS,                     &
                               MATRIX_NUMB_SSC, MATRIX_MASS_SSC_SU, MATRIX_MASS_SSC_SS,                     &
                               MATRIX_NUMB_OCC, MATRIX_MASS_OCC_SU, MATRIX_MASS_OCC_OC,                     &
                               MATRIX_NUMB_BC1, MATRIX_MASS_BC1_SU, MATRIX_MASS_BC1_BC,                     &
                               MATRIX_NUMB_BC2, MATRIX_MASS_BC2_SU, MATRIX_MASS_BC2_BC,                     &
                               MATRIX_NUMB_BC3, MATRIX_MASS_BC3_SU, MATRIX_MASS_BC3_BC,                     &
                               MATRIX_NUMB_DBC, MATRIX_MASS_DBC_SU, MATRIX_MASS_DBC_BC, MATRIX_MASS_DBC_DU, &
                               MATRIX_NUMB_BOC, MATRIX_MASS_BOC_SU, MATRIX_MASS_BOC_BC, MATRIX_MASS_BOC_OC, &
                               MATRIX_NUMB_BCS, MATRIX_MASS_BCS_SU, MATRIX_MASS_BCS_BC,                     &
                               MATRIX_NUMB_MXX, MATRIX_MASS_MXX_SU, MATRIX_MASS_MXX_BC, MATRIX_MASS_MXX_OC, &
                                                MATRIX_MASS_MXX_DU, MATRIX_MASS_MXX_SS/)

    ASSERT_(any(matrix_aerosol_indexes > 0))
    deallocate(matrix_aerosol_indexes, __STAT__)
    

!   All done
!   --------
    RETURN_(ESMF_SUCCESS)

   end subroutine Initialize_


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Run_ --- Runs MATRIXchem
!
! !INTERFACE:
!

   subroutine Run_(GC, IMPORT, EXPORT, clock, rc)

! !USES:

    implicit none

! !INPUT PARAMETERS:

   type(ESMF_Clock),    intent(inout) :: clock  ! The clock

! !OUTPUT PARAMETERS:

   type(ESMF_GridComp), intent(inout) :: GC     ! Grid Component
   type(ESMF_State),    intent(inout) :: IMPORT ! Import State
   type(ESMF_State),    intent(inout) :: EXPORT ! Export State
   integer,             intent(out)   :: rc     ! Error return code:
                                                !  0 - all is well
                                                !  1 - 

! !DESCRIPTION: This is a simple ESMF wrapper.
!
! !REVISION HISTORY:
!
!  13 Mar 2015    A.Darmenov  Coupled the MATRIX module with GEOS-5.
!  27 Feb 2005    da Silva    First crack.
!
!EOP
!-------------------------------------------------------------------------

                              __Iam__('Run_')
   
    type(MATRIX_state), pointer   :: self               ! Legacy state
    type(ESMF_Grid)               :: GRID               ! Grid
    type(ESMF_Config)             :: CF                 ! Universal Config 

    integer                       :: im_World, jm_World ! Global 2D Dimensions
    integer                       :: im, jm, lm         ! 3D Dimensions
    real(ESMF_KIND_R4), pointer   :: lons(:,:)          ! Longitudes
    real(ESMF_KIND_R4), pointer   :: lats(:,:)          ! Latitudes

    integer                       :: nymd, nhms         ! date, time
    real                          :: cdt                ! time step in secs

    character(len=ESMF_MAXSTR)    :: comp_name          ! name of the component

    integer                       :: i, j, l            ! current 3D indexes

!   --- MATRIX ---
    real(8) :: T_            ! absolute temperature [K]
    real(8) :: RH_           ! relative humidity [0-1]
    real(8) :: P_            ! ambient pressure [Pa]  
    real(8) :: so4_aq_rate_  ! in-cloud SO4 production rate [ug/m^3/s]
    real(8) :: w_updraft_    ! cloud updraft velocity [m/s]

    real(8) :: aerosol_(MATRIX_N_AEROSOLS)                             ! aerosol conc. [ug/m^3] or [#/m^3]
    real(8) :: gas_(MATRIX_N_GASES)                                    ! gas-phase conc. [ug/m^3]
    real(8) :: emissions_(MATRIX_N_EMIS_SPECIES)                       ! mass emission rates [ug/m^3/s]
    real(8) :: diagnostics_(MATRIX_N_AEROSOL_DIAG, MATRIX_N_AEROSOLS)  ! budget or tendency diagnostics [ug/m^3/s] or [#/m^3/s]
    real(8) :: species_mass_(MATRIX_N_MASS_SPECIES+2)                  ! total mass concentration of each model species (SU, BC, OC, DU, SS, NO3, NH4; but not water)

    real(8) :: f_emiss                                                 ! units factor

    integer, parameter              :: ss_emiss_method = 1
    real, parameter,   dimension(2) :: ssa_size_range = (/0.05, 0.5/)  ! lower and upper size|diameter range in 'um'
    real, parameter,   dimension(2) :: ssc_size_range = (/ 0.5, 8.0/)  ! lower and upper size|diameter range in 'um'
    real, pointer, dimension(:,:)   :: ssa_emiss_mass, ssa_emiss_num   ! fine
    real, pointer, dimension(:,:)   :: ssc_emiss_mass, ssc_emiss_num   ! coarse
    real, pointer, dimension(:,:)   :: f_grid_efficiency
    real, pointer, dimension(:,:)   :: w10m

    real, parameter,   dimension(2) :: duf_size_range = (/0.1, 2.5/)  ! lower and upper size|diameter range in 'um'
    real, parameter,   dimension(2) :: duc_size_range = (/2.0,10.0/)  ! lower and upper size|diameter range in 'um'
    real, pointer, dimension(:,:)   :: duf_emiss_mass, duf_emiss_num  ! fine
    real, pointer, dimension(:,:)   :: duc_emiss_mass, duc_emiss_num  ! coarse 
    real, pointer, dimension(:,:)   :: dust_emiss_tot                 ! total is the sum of the emissions in the 5 GOCART bins




!   Declare pointers to IMPORT/EXPORT/INTERNAL states 
!   -------------------------------------------------
#include "MATRIXchem_DeclarePointer___.h"
  
!   Get my name and set-up traceback handle
!   ---------------------------------------
    call ESMF_GridCompGet(GC, name=comp_name, __RC__)
    Iam = trim(comp_name) // trim(Iam)

!   Get pointers to IMPORT/EXPORT/INTERNAL states 
!   ---------------------------------------------
#include "MATRIXchem_GetPointer___.h"


!   Start the main timer 
!   --------------------
    call MAPL_TimerOn(MetaComp, 'RUN', __RC__)


!   Extract relevant runtime information
!   ------------------------------------
    call extract_(GC, CLOCK, self, GRID, CF, &
                  im_World, jm_World,        &
                  im, jm, lm, lons, lats,    &
                  nymd, nhms, cdt, __RC__)



!   Aerosol emissions
!   -----------------
    call MAPL_TimerOn(MetaComp, '-EMISSIONS', __RC__)


!   Seasalt emissions
!   -----------------
    allocate(ssa_emiss_mass(im,jm),    __STAT__)
    allocate(ssa_emiss_num (im,jm),    __STAT__)
    
    allocate(ssc_emiss_mass(im,jm),    __STAT__)
    allocate(ssc_emiss_num (im,jm),    __STAT__)

    allocate(f_grid_efficiency(im,jm), __STAT__) 
    allocate(w10m(im,jm),              __STAT__)

    !  define 10-m wind speed
    w10m = sqrt(u10m*u10m + v10m*v10m)

    !  fefine grid emission efficiency
    f_grid_efficiency = 0.0
    where(LWI < 1) f_grid_efficiency = 1.0           ! water points (should use fractional water cover)

    ! seasalt emissions: accumulation mode
    ssa_emiss_num  = 0.0
    ssa_emiss_mass = 0.0

    call SeasaltEmission(0.5*ssa_size_range(1), 0.5*ssa_size_range(2),  ss_emiss_method, &
                         w10m, ustar, ssa_emiss_mass, ssa_emiss_num, __RC__)

    ssa_emiss_mass = self%f_emiss_seasalt * f_grid_efficiency * ssa_emiss_mass
    ssa_emiss_num  = self%f_emiss_seasalt * f_grid_efficiency * ssa_emiss_num

    ! seasalt emissions: coarse mode
    ssc_emiss_num  = 0.0
    ssc_emiss_mass = 0.0

    call SeasaltEmission(0.5*ssc_size_range(1), 0.5*ssc_size_range(2), ss_emiss_method, &
                         w10m, ustar, ssc_emiss_mass, ssc_emiss_num, __RC__)

    ssc_emiss_mass = self%f_emiss_seasalt * f_grid_efficiency * ssc_emiss_mass
    ssc_emiss_num  = self%f_emiss_seasalt * f_grid_efficiency * ssc_emiss_num
 
    if (associated(emiss_ssa)) emiss_ssa = ssa_emiss_mass
    if (associated(emiss_ssc)) emiss_ssc = ssc_emiss_mass
    if (associated(emiss_ss )) emiss_ss  = ssa_emiss_mass + ssc_emiss_mass


!   Dust emissions
!   -----------------
    allocate(dust_emiss_tot(im,jm), __STAT__)
    allocate(duf_emiss_mass(im,jm), __STAT__)
    allocate(duf_emiss_num (im,jm), __STAT__)
    allocate(duc_emiss_mass(im,jm), __STAT__)
    allocate(duc_emiss_num (im,jm), __STAT__)


    dust_emiss_tot = 0.0
    call MAM_DustEmissionGOCART(1, im, 1, jm, lm, frlake, wet1, lwi, u10m, v10m, dust_emiss_tot, __RC__)

    ! apply the dust emission tuning coefficient [kg s2 m-5] and Ginoux dust source function
    dust_emiss_tot = (self%f_emiss_dust * 1e-9) * GINOUX_DU * dust_emiss_tot

    call MAM_DustEmission(1, im, 1, jm, lm, 1e-6*duf_size_range(1), 1e-6*duf_size_range(2), dust_emiss_tot, duf_emiss_mass, duf_emiss_num, __RC__)
    call MAM_DustEmission(1, im, 1, jm, lm, 1e-6*duc_size_range(1), 1e-6*duc_size_range(2), dust_emiss_tot, duc_emiss_mass, duc_emiss_num, __RC__)

    if (associated(emiss_duf)) emiss_duf = duf_emiss_mass
    if (associated(emiss_duc)) emiss_duc = duc_emiss_mass
    if (associated(emiss_du )) emiss_du  = duf_emiss_mass + duc_emiss_mass



#ifdef __MATRIX_TODO__
!   Dust emissions
!   -----------------
    call MX_DU_Emissions(self%scheme, import, export, self%qa, self%femisDU, self%dt, __RC__)

!   Black Carbon emissions
!   ----------------------
    call MX_BC_Emissions(self%scheme, import, export, self%qa, self%dt, __RC__)

!   Organic Carbon emissions
!   ----------------------
    call MX_OC_Emissions(self%scheme, import, export, self%qa, self%pom_oc_ratio, self%dt, __RC__)
#endif
    call MAPL_TimerOff(MetaComp, '-EMISSIONS', __RC__)



!   Aerosol microphysics
!   --------------------
    call MAPL_TimerOn(MetaComp, '-MICROPHYSICS', __RC__)

    allocate(MATRIX_VDDEP_AERO(im,jm,MATRIX_N_AEROSOL_MODES,2), __STAT__)  ! [m/s]
    allocate(MATRIX_DIAMETER(im,jm,lm,MATRIX_N_AEROSOL_MODES),  __STAT__)  ! [m  ]
    allocate(MATRIX_NACTV(im,jm,lm,MATRIX_N_AEROSOL_MODES),     __STAT__)  ! [#/m^3]
    allocate(MATRIX_CCNSS(im,jm,lm,MATRIX_N_AEROSOL_MODES,3),   __STAT__)  ! [#/m^3]


    if (associated(total_su)) total_su = 0.0
    if (associated(total_du)) total_du = 0.0
    if (associated(total_ss)) total_ss = 0.0
    if (associated(total_oc)) total_oc = 0.0
    if (associated(total_bc)) total_bc = 0.0

    do l = 1, lm
        do j = 1, jm
            do i = 1, im
                
                ! set indexes for MATRIX
                MATRIX_I = i
                MATRIX_J = j
                MATRIX_L = lm - l + 1

                ! set atm state
                T_         = T(i,j,l)
                RH_        = min(RH2(i,j,l), 0.95)
                P_         = 0.5*(PLE(i,j,l) + PLE(i,j,l-1))
                w_updraft_ = 0.05d0

                ! aerosol emissions and microphysics
                aerosol_   = 0.0d0
                gas_       = 0.0d0
                emissions_ = 0.0d0

                
                ! in-cloud SO4 production rate [ug/m^3/s] 
                so4_aq_rate_ = tiny(0.0d0)  !FCLD(i,j,l) * 1e3*4.43D-11      !0.1 * FCLD(i,j,l) * SO2(i,j,l) * 1.0d9 * (64.066 / MATRIX_MW_AIR) * airdens(i,j,l) / self%dt

                ! set concetrations of gases, 'ug(constituent) m-3'
                gas_(1) = tiny(0.0d0) ! h2so4(i,j,l) * 1.0d9 * (MATRIX_MW_H2SO4 / MATRIX_MW_AIR) * airdens(i,j,l)    
                gas_(2) = tiny(0.0d0) ! 0.01 * gas_(1)!1.0d9 * (MATRIX_MW_HNO3  / MATRIX_MW_AIR) * airdens(i,j,l)    ! ug(HNO3)  m-3
                gas_(3) = tiny(0.0d0) ! NH3(i,j,l)   * 1.0d9 * (MATRIX_MW_NH3   / MATRIX_MW_AIR) * airdens(i,j,l)

                where (gas_ < tiny(0.0d0))
                    gas_ = 1e3*tiny(0.0d0)
                end where


                ! set concentrations of aerosols: mass is in 'ug m-3', number is in '# m-3'
                ! total
                aerosol_(MATRIX_MASS_NO3)    = M_NO3(i,j,l)
                aerosol_(MATRIX_MASS_NH4)    = M_NH4(i,j,l)
                aerosol_(MATRIX_MASS_H2O)    = M_H2O(i,j,l)
                ! AKK
                aerosol_(MATRIX_MASS_AKK_SU) = M_AKK_SU(i,j,l)
                aerosol_(MATRIX_NUMB_AKK)    = N_AKK(i,j,l)
                ! ACC 
                aerosol_(MATRIX_MASS_ACC_SU) = M_ACC_SU(i,j,l)
                aerosol_(MATRIX_NUMB_ACC)    = N_ACC(i,j,l)
                ! DD1
                aerosol_(MATRIX_MASS_DD1_SU) = M_DD1_SU(i,j,l)
                aerosol_(MATRIX_MASS_DD1_DU) = M_DD1_DU(i,j,l)
                aerosol_(MATRIX_NUMB_DD1)    = N_DD1(i,j,l)
                ! DS1
                aerosol_(MATRIX_MASS_DS1_SU) = M_DS1_SU(i,j,l)
                aerosol_(MATRIX_MASS_DS1_DU) = M_DS1_DU(i,j,l)
                aerosol_(MATRIX_NUMB_DS1)    = N_DS1(i,j,l)
                ! DD2
                aerosol_(MATRIX_MASS_DD2_SU) = M_DD2_SU(i,j,l)
                aerosol_(MATRIX_MASS_DD2_DU) = M_DD2_DU(i,j,l)
                aerosol_(MATRIX_NUMB_DD2)    = N_DD2(i,j,l)
                ! DS2
                aerosol_(MATRIX_MASS_DS2_SU) = M_DS2_SU(i,j,l)
                aerosol_(MATRIX_MASS_DS2_DU) = M_DS2_DU(i,j,l)
                aerosol_(MATRIX_NUMB_DS2)    = N_DS2(i,j,l)

                ! SSA
                aerosol_(MATRIX_MASS_SSA_SU) = M_SSA_SU(i,j,l)
                aerosol_(MATRIX_MASS_SSA_SS) = M_SSA_SS(i,j,l)
                aerosol_(MATRIX_NUMB_SSA)    = N_SSA(i,j,l)
                 !SSC
                aerosol_(MATRIX_MASS_SSC_SU) = M_SSC_SU(i,j,l)
                aerosol_(MATRIX_MASS_SSC_SS) = M_SSC_SS(i,j,l)
                aerosol_(MATRIX_NUMB_SSC)    = N_SSC(i,j,l)
                ! OCC
                aerosol_(MATRIX_MASS_OCC_SU) = M_OCC_SU(i,j,l)
                aerosol_(MATRIX_MASS_OCC_OC) = M_OCC_OC(i,j,l)
                aerosol_(MATRIX_NUMB_OCC)    = N_OCC(i,j,l)
                ! BC1
                aerosol_(MATRIX_MASS_BC1_SU) = M_BC1_SU(i,j,l)
                aerosol_(MATRIX_MASS_BC1_BC) = M_BC1_BC(i,j,l)
                aerosol_(MATRIX_NUMB_BC1)    = N_BC1(i,j,l)
                ! BC2
                aerosol_(MATRIX_MASS_BC2_SU) = M_BC2_SU(i,j,l)
                aerosol_(MATRIX_MASS_BC2_BC) = M_BC2_BC(i,j,l)
                aerosol_(MATRIX_NUMB_BC2)    = N_BC2(i,j,l)
                ! BC3
                aerosol_(MATRIX_MASS_BC3_SU) = M_BC3_SU(i,j,l)
                aerosol_(MATRIX_MASS_BC3_BC) = M_BC3_BC(i,j,l)
                aerosol_(MATRIX_NUMB_BC3)    = N_BC3(i,j,l)
                ! DBC
                aerosol_(MATRIX_MASS_DBC_SU) = M_DBC_SU(i,j,l)
                aerosol_(MATRIX_MASS_DBC_BC) = M_DBC_BC(i,j,l)
                aerosol_(MATRIX_MASS_DBC_DU) = M_DBC_DU(i,j,l)
                aerosol_(MATRIX_NUMB_DBC)    = N_DBC(i,j,l)
                ! BOC
                aerosol_(MATRIX_MASS_BOC_SU) = M_BOC_SU(i,j,l)
                aerosol_(MATRIX_MASS_BOC_BC) = M_BOC_BC(i,j,l)
                aerosol_(MATRIX_MASS_BOC_OC) = M_BOC_OC(i,j,l)
                aerosol_(MATRIX_NUMB_BOC)    = N_BOC(i,j,l)
                ! BCS
                aerosol_(MATRIX_MASS_BCS_SU) = M_BCS_SU(i,j,l)
                aerosol_(MATRIX_MASS_BCS_BC) = M_BCS_BC(i,j,l)
                aerosol_(MATRIX_NUMB_BCS)    = N_BCS(i,j,l)
                ! MXX
                aerosol_(MATRIX_MASS_MXX_SU) = M_MXX_SU(i,j,l)
                aerosol_(MATRIX_MASS_MXX_BC) = M_MXX_BC(i,j,l)
                aerosol_(MATRIX_MASS_MXX_OC) = M_MXX_OC(i,j,l)
                aerosol_(MATRIX_MASS_MXX_DU) = M_MXX_DU(i,j,l)
                aerosol_(MATRIX_MASS_MXX_SS) = M_MXX_SS(i,j,l)
                aerosol_(MATRIX_NUMB_MXX)    = N_MXX(i,j,l)

                if (l == lm) then
                    aerosol_ = aerosol_ * exp(-0.5d-2 * self%dt * (MAPL_GRAV * airdens(i,j,l) / delp(i,j,l)))
                endif

                where (aerosol_ < tiny(0.0d0))
                    aerosol_ = 1e3*tiny(0.0d0)
                end where

             


                ! set emissions
                if (l == lm) then

                ! did I get the emissions indexes right?
                !    REAL, DIMENSION(NEMIS_SPCS) :: EMIS_DENS = (/  EMIS_DENS_SULF,
                !                                                   EMIS_DENS_SULF, 
                !                                                   EMIS_DENS_BCAR, 
                !                                                   EMIS_DENS_OCAR,
                !                                                   EMIS_DENS_DUST, 
                !                                                   EMIS_DENS_SEAS, 
                !                                                   EMIS_DENS_SEAS,
                !                                                   EMIS_DENS_BOCC, 
                !                                                   EMIS_DENS_BOCC, 
                !                                                   EMIS_DENS_DUST /)
                ! 
                !    CHEM_SPC_NAME(NMASS_SPCS) = (/'SULF','BCAR','OCAR','DUST','SEAS'/)
                !    EMIS_SPCS_MAP             = (/1,1,2,3,4,5,5,2,3,4/)

                    f_emiss  = 1.0d9 * MAPL_GRAV * airdens(i,j,l) / delp(i,j,l)

                    emissions_(1) = f_emiss * 0.025*(SO2_EMIS_FIRES(i,j)     + &
                                                     SO2_EMIS_NONENERGY(i,j) + &
                                                     SO2_EMIS_ENERGY(i,j)    + &
                                                     SO2_EMIS_SHIPPING(i,j)) ! AKK SU: volc + biomass
                    emissions_(2) = f_emiss * SO4_EMIS_SHIP(i,j)             ! ACC SU: volc + biomass
                    emissions_(3) = f_emiss * BC_EMIS_FIRE(i,j)              ! BC1 BC: 
                    emissions_(4) = f_emiss * OC_EMIS_FIRE(i,j)              ! OC 
                    emissions_(5) = f_emiss * duf_emiss_mass(i,j)            ! DU 
                    emissions_(6) = f_emiss * ssa_emiss_mass(i,j)            ! Accumulation mode: 0.01 -- 0.5 microns
                    emissions_(7) = f_emiss * ssc_emiss_mass(i,j)            ! Coarse mode:       0.50 -- 8.0 microns
                    emissions_(8) = f_emiss * (BC_EMIS_BIOFUEL(i,j) + BC_EMIS_FOSSILFUEL(i,j) + BC_EMIS_SHIP(i,j))
                    emissions_(9) = f_emiss * (OC_EMIS_BIOFUEL(i,j) + OC_EMIS_FOSSILFUEL(i,j) + OC_EMIS_SHIP(i,j))
                    emissions_(10)= f_emiss * duc_emiss_mass(i,j)
                end if


                where (aerosol_ < tiny(0.0d0))
                    aerosol_ = tiny(0.0d0)
                end where


                MATRIX_DIAMETER(MATRIX_I,MATRIX_J,MATRIX_L, 1) = DGN_AKK(i,j,l)
                MATRIX_DIAMETER(MATRIX_I,MATRIX_J,MATRIX_L, 2) = DGN_ACC(i,j,l)
                MATRIX_DIAMETER(MATRIX_I,MATRIX_J,MATRIX_L, 3) = DGN_DD1(i,j,l)
                MATRIX_DIAMETER(MATRIX_I,MATRIX_J,MATRIX_L, 4) = DGN_DS1(i,j,l)
                MATRIX_DIAMETER(MATRIX_I,MATRIX_J,MATRIX_L, 5) = DGN_DD2(i,j,l)
                MATRIX_DIAMETER(MATRIX_I,MATRIX_J,MATRIX_L, 6) = DGN_DS2(i,j,l)
                MATRIX_DIAMETER(MATRIX_I,MATRIX_J,MATRIX_L, 7) = DGN_SSA(i,j,l)
                MATRIX_DIAMETER(MATRIX_I,MATRIX_J,MATRIX_L, 8) = DGN_SSC(i,j,l)
                MATRIX_DIAMETER(MATRIX_I,MATRIX_J,MATRIX_L, 9) = DGN_OCC(i,j,l)
                MATRIX_DIAMETER(MATRIX_I,MATRIX_J,MATRIX_L,10) = DGN_BC1(i,j,l)
                MATRIX_DIAMETER(MATRIX_I,MATRIX_J,MATRIX_L,11) = DGN_BC2(i,j,l)
                MATRIX_DIAMETER(MATRIX_I,MATRIX_J,MATRIX_L,12) = DGN_BC3(i,j,l)
                MATRIX_DIAMETER(MATRIX_I,MATRIX_J,MATRIX_L,13) = DGN_DBC(i,j,l)
                MATRIX_DIAMETER(MATRIX_I,MATRIX_J,MATRIX_L,14) = DGN_BOC(i,j,l)
                MATRIX_DIAMETER(MATRIX_I,MATRIX_J,MATRIX_L,15) = DGN_BCS(i,j,l)
                MATRIX_DIAMETER(MATRIX_I,MATRIX_J,MATRIX_L,16) = DGN_MXX(i,j,l)


                diagnostics_ = tiny(0.0)

                call SPCMASSES(aerosol_, gas_, species_mass_)
                call MATRIX(aerosol_, gas_, emissions_, self%dt, T_, RH_, P_, so4_aq_rate_, w_updraft_, diagnostics_)

                ! diagnostics
                if (associated(total_su)) total_su(i,j,l) = sum(aerosol_(MATRIX_SU_MAP(:)))
                if (associated(total_du)) total_du(i,j,l) = sum(aerosol_(MATRIX_DU_MAP(:)))
                if (associated(total_ss)) total_ss(i,j,l) = sum(aerosol_(MATRIX_SS_MAP(:)))
                if (associated(total_oc)) total_oc(i,j,l) = sum(aerosol_(MATRIX_OC_MAP(:)))
                if (associated(total_bc)) total_bc(i,j,l) = sum(aerosol_(MATRIX_BC_MAP(:)))

                ! total
                M_NO3(i,j,l)    = aerosol_(MATRIX_MASS_NO3)
                M_NH4(i,j,l)    = aerosol_(MATRIX_MASS_NH4)
                M_H2O(i,j,l)    = aerosol_(MATRIX_MASS_H2O)
                ! AKK
                M_AKK_SU(i,j,l) = aerosol_(MATRIX_MASS_AKK_SU)
                N_AKK(i,j,l)    = aerosol_(MATRIX_NUMB_AKK)
                ! ACC 
                M_ACC_SU(i,j,l) = aerosol_(MATRIX_MASS_ACC_SU)
                N_ACC(i,j,l)    = aerosol_(MATRIX_NUMB_ACC)
                ! DD1
                M_DD1_SU(i,j,l) = aerosol_(MATRIX_MASS_DD1_SU)
                M_DD1_DU(i,j,l) = aerosol_(MATRIX_MASS_DD1_DU)
                N_DD1(i,j,l)    = aerosol_(MATRIX_NUMB_DD1)
                ! DS1
                M_DS1_SU(i,j,l) = aerosol_(MATRIX_MASS_DS1_SU)
                M_DS1_DU(i,j,l) = aerosol_(MATRIX_MASS_DS1_DU)
                N_DS1(i,j,l)    = aerosol_(MATRIX_NUMB_DS1)
                ! DD2
                M_DD2_SU(i,j,l) = aerosol_(MATRIX_MASS_DD2_SU)
                M_DD2_DU(i,j,l) = aerosol_(MATRIX_MASS_DD2_DU)
                N_DD2(i,j,l)    = aerosol_(MATRIX_NUMB_DD2)
                ! DS2
                M_DS2_SU(i,j,l) = aerosol_(MATRIX_MASS_DS2_SU)
                M_DS2_DU(i,j,l) = aerosol_(MATRIX_MASS_DS2_DU)
                N_DS2(i,j,l)    = aerosol_(MATRIX_NUMB_DS2)
                ! SSA
                M_SSA_SU(i,j,l) = aerosol_(MATRIX_MASS_SSA_SU)
                M_SSA_SS(i,j,l) = aerosol_(MATRIX_MASS_SSA_SS)
                N_SSA           = aerosol_(MATRIX_NUMB_SSA)
                 !SSC
                M_SSC_SU(i,j,l) = aerosol_(MATRIX_MASS_SSC_SU)
                M_SSC_SS(i,j,l) = aerosol_(MATRIX_MASS_SSC_SS)
                N_SSC (i,j,l)   = aerosol_(MATRIX_NUMB_SSC)
                ! OCC
                M_OCC_SU(i,j,l) = aerosol_(MATRIX_MASS_OCC_SU)
                M_OCC_OC(i,j,l) = aerosol_(MATRIX_MASS_OCC_OC)
                N_OCC(i,j,l)    = aerosol_(MATRIX_NUMB_OCC)
                ! BC1
                M_BC1_SU(i,j,l) = aerosol_(MATRIX_MASS_BC1_SU)
                M_BC1_BC(i,j,l) = aerosol_(MATRIX_MASS_BC1_BC)
                N_BC1(i,j,l)    = aerosol_(MATRIX_NUMB_BC1)
                ! BC2
                M_BC2_SU(i,j,l) = aerosol_(MATRIX_MASS_BC2_SU)
                M_BC2_BC(i,j,l) = aerosol_(MATRIX_MASS_BC2_BC)
                N_BC2(i,j,l)    = aerosol_(MATRIX_NUMB_BC2)
                ! BC3
                M_BC3_SU(i,j,l) = aerosol_(MATRIX_MASS_BC3_SU)
                M_BC3_BC(i,j,l) = aerosol_(MATRIX_MASS_BC3_BC)
                N_BC3(i,j,l)    = aerosol_(MATRIX_NUMB_BC3)
                ! DBC
                M_DBC_SU(i,j,l) = aerosol_(MATRIX_MASS_DBC_SU)
                M_DBC_BC(i,j,l) = aerosol_(MATRIX_MASS_DBC_BC)
                M_DBC_DU(i,j,l) = aerosol_(MATRIX_MASS_DBC_DU)
                N_DBC(i,j,l)    = aerosol_(MATRIX_NUMB_DBC)
                ! BOC
                M_BOC_SU(i,j,l) = aerosol_(MATRIX_MASS_BOC_SU)
                M_BOC_BC(i,j,l) = aerosol_(MATRIX_MASS_BOC_BC)
                M_BOC_OC(i,j,l) = aerosol_(MATRIX_MASS_BOC_OC)
                N_BOC(i,j,l)    = aerosol_(MATRIX_NUMB_BOC)
                ! BCS
                M_BCS_SU(i,j,l) = aerosol_(MATRIX_MASS_BCS_SU)
                M_BCS_BC(i,j,l) = aerosol_(MATRIX_MASS_BCS_BC)
                N_BCS(i,j,l)    = aerosol_(MATRIX_NUMB_BCS)
                ! MXX
                M_MXX_SU(i,j,l) = aerosol_(MATRIX_MASS_MXX_SU)
                M_MXX_BC(i,j,l) = aerosol_(MATRIX_MASS_MXX_BC)
                M_MXX_OC(i,j,l) = aerosol_(MATRIX_MASS_MXX_OC)
                M_MXX_DU(i,j,l) = aerosol_(MATRIX_MASS_MXX_DU)
                M_MXX_SS(i,j,l) = aerosol_(MATRIX_MASS_MXX_SS)
                N_MXX(i,j,l)    = aerosol_(MATRIX_NUMB_MXX)

                ! update sizes
                DGN_AKK(i,j,l)  = MATRIX_DIAMETER(MATRIX_I,MATRIX_J,MATRIX_L, 1) 
                DGN_ACC(i,j,l)  = MATRIX_DIAMETER(MATRIX_I,MATRIX_J,MATRIX_L, 2)
                DGN_DD1(i,j,l)  = MATRIX_DIAMETER(MATRIX_I,MATRIX_J,MATRIX_L, 3)
                DGN_DS1(i,j,l)  = MATRIX_DIAMETER(MATRIX_I,MATRIX_J,MATRIX_L, 4)
                DGN_DD2(i,j,l)  = MATRIX_DIAMETER(MATRIX_I,MATRIX_J,MATRIX_L, 5)
                DGN_DS2(i,j,l)  = MATRIX_DIAMETER(MATRIX_I,MATRIX_J,MATRIX_L, 6)
                DGN_SSA(i,j,l)  = MATRIX_DIAMETER(MATRIX_I,MATRIX_J,MATRIX_L, 7)
                DGN_SSC(i,j,l)  = MATRIX_DIAMETER(MATRIX_I,MATRIX_J,MATRIX_L, 8)
                DGN_OCC(i,j,l)  = MATRIX_DIAMETER(MATRIX_I,MATRIX_J,MATRIX_L, 9)
                DGN_BC1(i,j,l)  = MATRIX_DIAMETER(MATRIX_I,MATRIX_J,MATRIX_L,10)
                DGN_BC2(i,j,l)  = MATRIX_DIAMETER(MATRIX_I,MATRIX_J,MATRIX_L,11)
                DGN_BC3(i,j,l)  = MATRIX_DIAMETER(MATRIX_I,MATRIX_J,MATRIX_L,12)
                DGN_DBC(i,j,l)  = MATRIX_DIAMETER(MATRIX_I,MATRIX_J,MATRIX_L,13)
                DGN_BOC(i,j,l)  = MATRIX_DIAMETER(MATRIX_I,MATRIX_J,MATRIX_L,14)
                DGN_BCS(i,j,l)  = MATRIX_DIAMETER(MATRIX_I,MATRIX_J,MATRIX_L,15)
                DGN_MXX(i,j,l)  = MATRIX_DIAMETER(MATRIX_I,MATRIX_J,MATRIX_L,16)
             end do 
        end do
    end do

    deallocate(MATRIX_NACTV,      __STAT__)
    deallocate(MATRIX_CCNSS,      __STAT__)
    deallocate(MATRIX_DIAMETER,   __STAT__)
    deallocate(MATRIX_VDDEP_AERO, __STAT__)
    

    call MAPL_TimerOff(MetaComp, '-MICROPHYSICS', __RC__)

    deallocate(ssa_emiss_mass,    __STAT__)
    deallocate(ssa_emiss_num,     __STAT__)
    deallocate(ssc_emiss_mass,    __STAT__)
    deallocate(ssc_emiss_num,     __STAT__)
    deallocate(f_grid_efficiency, __STAT__) 
    deallocate(w10m,              __STAT__)

    deallocate(dust_emiss_tot, __STAT__)
    deallocate(duf_emiss_mass, __STAT__)
    deallocate(duf_emiss_num,  __STAT__)
    deallocate(duc_emiss_mass, __STAT__)
    deallocate(duc_emiss_num,  __STAT__)


!   Stop the main timer 
!   -------------------
    call MAPL_TimerOff(MetaComp, 'RUN', __RC__)

!   All done
!   --------
    RETURN_(ESMF_SUCCESS)

   end subroutine Run_

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Finalize_ --- Finalize MATRIXchem
!
! !INTERFACE:
!

   subroutine Finalize_(GC, IMPORT, EXPORT, clock, rc)

! !USES:

    implicit none

! !INPUT PARAMETERS:

    type(ESMF_Clock),    intent(inout) :: clock  ! The clock

! !OUTPUT PARAMETERS:

    type(ESMF_GridComp), intent(inout) :: GC     ! Grid Component
    type(ESMF_State),    intent(inout) :: IMPORT ! Import State
    type(ESMF_State),    intent(inout) :: EXPORT ! Export State
    integer,             intent(out)   :: rc     ! Error return code:
                                                 !  0 - all is well
                                                 !  1 - 

! !DESCRIPTION: This is a simple ESMF wrapper.
!
! !REVISION HISTORY:
!
!  01Dec2009    da Silva  First crack.
!
!EOP
!-------------------------------------------------------------------------

                              __Iam__('Finalize_')
   
    type(MATRIX_state), pointer   :: self               ! Legacy state
    type(ESMF_Grid)               :: GRID               ! Grid
    type(ESMF_Config)             :: CF                 ! Universal Config 

    integer                       :: im_World, jm_World ! Global 2D Dimensions
    integer                       :: im, jm, lm         ! 3D Dimensions
    real(ESMF_KIND_R4), pointer   :: lons(:,:)          ! Longitudes
    real(ESMF_KIND_R4), pointer   :: lats(:,:)          ! Latitudes

    integer                       :: nymd, nhms         ! date, time
    real                          :: cdt                ! time step in secs

    character(len=ESMF_MAXSTR)    :: comp_name          ! name of the component

!   Get my name and set-up traceback handle
!   ---------------------------------------
    call ESMF_GridCompGet(GC, name=comp_name, __RC__)
    Iam = trim(comp_name) // trim(Iam)

!   Finalize MAPL Generic
!   ---------------------
    call MAPL_GenericFinalize(GC, IMPORT, EXPORT, clock,  __RC__)

!   Extract relevant runtime information
!   ------------------------------------
    call extract_(GC, CLOCK, self, GRID, CF, &
                  im_World, jm_World,        &
                  im, jm, lm, lons, lats,    &
                  nymd, nhms, cdt, __RC__)

!   All done
!   --------
    RETURN_(ESMF_SUCCESS)

   end subroutine Finalize_

!.......................................................................

   subroutine extract_(GC, clock,              &
                           myState, GRID, CF,  &
                           im_World, jm_World, &
                           im, jm, lm,         &
                           lons, lats,         &
                           nymd, nhms,         &
                           cdt, rc)

    type(ESMF_GridComp),         intent(inout) :: GC                 ! Grid Comp object
    type(ESMF_Clock),            intent(in)    :: clock              ! Clock

    type(MATRIX_state), pointer, intent(out)   :: myState            ! Legacy state
    type(ESMF_Grid),             intent(out)   :: GRID               ! Grid
    type(ESMF_Config),           intent(out)   :: CF                 ! Universal Config 

    integer,                     intent(out)   :: im_World, jm_World ! Global 2D Dimensions
    integer,                     intent(out)   :: im, jm, lm         ! 3D Dimensions

    real(ESMF_KIND_R4), pointer                :: lons(:,:)          ! Longitudes
    real(ESMF_KIND_R4), pointer                :: lats(:,:)          ! Latitudes



    integer,                     intent(out)   :: nymd, nhms         ! date, time
    real,                        intent(out)   :: cdt                ! time step in secs
    integer, optional,           intent(out)   :: rc

!                            ---

                      __Iam__('extract_')

    character(len=ESMF_MAXSTR)   :: comp_name

    type(MAPL_MetaComp), pointer :: mgState      ! MAPL generic state
    type(MATRIX_Wrap)            :: wrap

    integer, dimension(3)        :: dims

    type(ESMF_Alarm)             :: run_alarm
    type(ESMF_TimeInterval)      :: ring_interval
    real(ESMF_KIND_R8)           :: time_step

    type(ESMF_Time)              :: time
    integer                      :: iyr, imm, idd, ihr, imn, isc


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
    call ESMF_UserCompGetInternalState(GC, 'MATRIX_state', wrap, STATUS)
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
    call ESMF_TimeGet(TIME, yy=iyr, mm=imm, dd=idd, h=ihr, m=imn, s=isc, __RC__)

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
    call ESMF_GridGet(GRID, localDE=0, staggerloc=ESMF_STAGGERLOC_CENTER, &
                                       computationalCount=dims, __RC__)
    im = dims(1)
    jm = dims(2)
    lm = dims(3)

!  Get horizontal coordinate variables
!  -----------------------------------
   call MAPL_Get(mgState, lons=lons, lats=lats,  __RC__)



    RETURN_(ESMF_SUCCESS)

   end subroutine extract_

end module MATRIXchem_GridCompMod
