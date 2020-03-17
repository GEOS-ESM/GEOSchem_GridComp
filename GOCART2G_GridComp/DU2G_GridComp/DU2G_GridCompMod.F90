#include "MAPL_Generic.h"

!=============================================================================
!BOP

! !MODULE: DU2G_GridCompMod - GOCART refactoring of the DU gridded component 

! !INTERFACE:
module DU2G_GridCompMod

!  !USES:
   use ESMF
   use MAPL
   use Chem_MieMod2G
   use Chem_AeroGeneric
   use iso_c_binding, only: c_loc, c_f_pointer, c_ptr

   use Chem_UtilMod          
   use GriddedEmission       ! Emissions
   
   implicit none
   private

   integer, parameter :: instanceComputational = 1
   integer, parameter :: instanceData          = 2

! !PUBLIC MEMBER FUNCTIONS:
   public  SetServices


! !DESCRIPTION: This module implements GOCART2G's Dust (DU) Gridded Component.

! !REVISION HISTORY:
! 16Oct2019  Sherman, da Silva, Darmenov, Clune -  First attempt at refactoring.

!EOP
!===========================================================================

   integer, parameter         :: NHRES = 6  ! DEV NOTE!!! should this be allocatable, and not a parameter?
   real, parameter :: radToDeg = 57.2957795

!  !Dust state
   type DU2G_GridComp
       type(Chem_Mie), dimension(2)    :: rad_MieTable, diag_MieTable
       real, allocatable      :: radius(:)      ! particle effective radius [um]
       real, allocatable      :: rlow(:)        ! particle effective radius lower bound [um]
       real, allocatable      :: rup(:)         ! particle effective radius upper bound [um]
       real, allocatable      :: sfrac(:)       ! fraction of total source
       real, allocatable      :: rhop(:)        ! soil class density [kg m-3]
       real                   :: Ch_DU(NHRES)   ! dust emission tuning coefficient [kg s2 m-5].
       real, allocatable      :: fscav(:)       ! scavenging efficiency
       real, allocatable      :: molwght(:)     ! molecular weight
       real, allocatable      :: fnum(:)        ! number of particles per kg mass
       real                   :: maringFlag     ! maring settling velocity correction

!      !Workspae for point emissions
       logical                :: doing_point_emissions = .FALSE.
       character(len=255)     :: point_emissions_srcfilen   ! filename for pointwise emissions
       integer                         :: nPts = -1
       integer, pointer, dimension(:)  :: pstart => null(), pend => null()
       real, pointer, dimension(:)     :: pLat  => null(), &
                                          pLon  => null(), &
                                          pBase => null(), &
                                          pTop  => null(), &
                                          pEmis => null()
   end type DU2G_GridComp

   type wrap_
      type (DU2G_GridComp), pointer     :: PTR => null()
   end type wrap_

contains

!============================================================================
!BOP

! !IROUTINE: SetServices 

! !INTERFACE:
  subroutine SetServices (GC, RC)

!   !ARGUMENTS:
    type (ESMF_GridComp), intent(INOUT)   :: GC  ! gridded component
    integer,              intent(  OUT)   :: RC  ! return code


!   !DESCRIPTION: This version uses MAPL_GenericSetServices, which sets
!     the Initialize and Finalize services to generic versions. It also
!     allocates our instance of a generic state and puts it in the 
!     gridded component (GC). Here we only set the two-stage run method
!     and declare the data services.

!   !REVISION HISTORY: 
!   16oct2019   E.Sherman  First attempt at refactoring

!EOP
!============================================================================

!
!   !Locals
    character (len=ESMF_MAXSTR)                 :: COMP_NAME
    type (ESMF_Config)                          :: cfg
    type (wrap_)                                :: wrap
    type (DU2G_GridComp), pointer               :: self

    character (len=ESMF_MAXSTR)                 :: field_name
    character (len=ESMF_MAXSTR), allocatable    :: aerosol_names(:)

    integer                                     :: n, i, nCols, n_bins
    real                                        :: DEFVAL
    logical                                     :: data_driven = .true.
    integer, parameter                          :: bins = 5  ! This should equal the number of bins. Is this how we want to handle this?

    __Iam__('SetServices')

!****************************************************************************
!   Begin...

!   Get my name and set-up traceback handle
!   ---------------------------------------
    call ESMF_GridCompGet (GC, NAME=COMP_NAME, __RC__)
    Iam = trim(COMP_NAME) //'::'// Iam

if (mapl_am_I_root()) print*,' test DU2G SetServices COMP_NAME = ',trim(COMP_NAME)


!   Wrap internal state for storing in GC
!   -------------------------------------
    allocate (self, __STAT__)
    wrap%ptr => self

!   Load resource file  
!   -------------------
    cfg = ESMF_ConfigCreate (__RC__)
    call ESMF_ConfigLoadFile (cfg, 'DU2G_GridComp_'//trim(COMP_NAME)//'.rc', rc=status)
    if (status /= 0) then
        if (mapl_am_i_root()) print*,'DU2G_GridComp_'//trim(COMP_NAME)//'.rc does not exist! Loading DU2G_GridComp_DU.rc instead'
        call ESMF_ConfigLoadFile (cfg, 'DU2G_GridComp_DU.rc', __RC__)
    end if


!   Get names of aerosols and write to aerosol_names to add to AERO State
!   ----------------------------------------------------------------------
    call ESMF_ConfigGetDim (cfg, n_bins, nCols, label=('variable_table::'), __RC__ )
    allocate (aerosol_names(n_bins), __STAT__)
    call ESMF_ConfigFindLabel (cfg, 'variable_table::', __RC__ )

    do i = 1, n_bins
        call ESMF_ConfigNextLine( cfg, __RC__ )
        call ESMF_ConfigGetAttribute( cfg, aerosol_names(i), __RC__ )
    end do

!   Parse config file into private internal state
!   ----------------------------------------------
    allocate(self%radius(n_bins), self%rlow(n_bins), self%rup(n_bins), self%sfrac(n_bins), &
             self%rhop(n_bins), self%fscav(n_bins), self%molwght(n_bins), self%fnum(n_bins), __STAT__)

    call ESMF_ConfigGetAttribute (cfg, self%radius,     label='particle_radius_microns:', __RC__)
    call ESMF_ConfigGetAttribute (cfg, self%rlow,       label='radius_lower:', __RC__)
    call ESMF_ConfigGetAttribute (cfg, self%rup,        label='radius_upper:', __RC__)
    call ESMF_ConfigGetAttribute (cfg, self%sfrac,      label='source_fraction:', __RC__)
    call ESMF_ConfigGetAttribute (cfg, self%rhop,       label='soil_density:', __RC__)
    call ESMF_ConfigGetAttribute (cfg, self%Ch_DU,      label='Ch_DU:', __RC__)
    call ESMF_ConfigGetAttribute (cfg, self%fscav,      label='fscav:', __RC__)
    call ESMF_ConfigGetAttribute (cfg, self%molwght,    label='molecular_weight:', __RC__)
    call ESMF_ConfigGetAttribute (cfg, self%fnum,       label='fnum:', __RC__)
    call ESMF_ConfigGetAttribute (cfg, self%maringFlag, label='maringFlag:', __RC__)
    call ESMF_ConfigGetAttribute (cfg, self%point_emissions_srcfilen, &
                                  label='point_emissions_srcfilen:', rc=status)
    if (status /= 0) then
        self%doing_point_emissions = .false.
    else
        if ( (index(self%point_emissions_srcfilen,'/dev/null')>0) ) then
            self%doing_point_emissions = .FALSE. ! disable it if no file specified
        else
            self%doing_point_emissions = .TRUE.  ! we are good to go
        end if
    end if

if(mapl_am_i_root()) print*,'DU2G self%Ch_DU',self%Ch_DU


!   Set entry points
!   ------------------------
    call MAPL_GridCompSetEntryPoint (GC, ESMF_Method_Initialize,  Initialize, __RC__)
    call MAPL_GridCompSetEntryPoint (GC, ESMF_Method_Run, Run, __RC__)
    call MAPL_GridCompSetEntryPoint (GC, ESMF_Method_Run, Run2, __RC__)

!   Is DU data driven?
!   ------------------
    call data_driven_ (COMP_NAME, data_driven, __RC__)


!   INTERNAL STATE
!   ---------------
!   Default internal state values
!   -----------------------------
    DEFVAL = 0.0

!   Aerosol Tracers to be transported
!   ---------------------------------
    do i = 1, n_bins
        call MAPL_AddInternalSpec(GC,                                    &
          short_name = trim(comp_name)//'::'//trim(aerosol_names(i)),    &
          long_name  = 'Dust Mixing Ratio (bin '//trim(aerosol_names(i))//')', &
          units      = 'kg kg-1',                                        &
          restart    = MAPL_RestartOptional,                             &
          default    = DEFVAL,                                           &
          friendlyto = trim(comp_name),                                  &
          dims       = MAPL_DimsHorzVert,                                &
          vlocation  = MAPL_VLocationCenter, __RC__)

if(mapl_am_i_root()) print*,'DU2G setservices shortname = ',trim(COMP_NAME)//'::'//trim(aerosol_names(i))
    end do


!   IMPORT STATE
!   -------------

!     Pressure thickness
!     ------------------
       call MAPL_AddImportSpec(GC,                            &
          SHORT_NAME = 'DELP',                                &
          LONG_NAME  = 'pressure_thickness',                  &
          UNITS      = 'Pa',                                  &
          DIMS       = MAPL_DimsHorzVert,                     &
          VLOCATION  = MAPL_VLocationCenter,                  &
          RESTART    = MAPL_RestartSkip,     __RC__)


    if (data_driven) then

!      Pressure at layer edges
!      -----------------------
       call MAPL_AddImportSpec(GC,                            &
          SHORT_NAME = 'PLE',                                 &
          LONG_NAME  = 'air_pressure',                        &
          UNITS      = 'Pa',                                  &
          DIMS       = MAPL_DimsHorzVert,                     &
          VLOCATION  = MAPL_VLocationEdge,                    &
          RESTART    = MAPL_RestartSkip,     __RC__)

!      RH: is between 0 and 1
!      ----------------------
       call MAPL_AddImportSpec(GC,                            &
          SHORT_NAME = 'RH2',                                 &
          LONG_NAME  = 'Rel_Hum_after_moist',                 &
          UNITS      = '1',                                   &
          DIMS       = MAPL_DimsHorzVert,                     &
          VLOCATION  = MAPL_VLocationCenter,                  &
          RESTART    = MAPL_RestartSkip,     __RC__)

       do i = 1, bins 
           write (field_name, '(A, I0.3)') '', i
            call MAPL_AddImportSpec(GC,                                     &
              short_name = 'climdu'//trim(field_name),                        &
              long_name  = 'Dust Mixing Ratio (bin '//trim(field_name)//')',  &
              units      = 'kg kg-1',                                         &
              restart    = MAPL_RestartSkip,                                  &
              dims       = MAPL_DimsHorzVert,                                 &
              vlocation  = MAPL_VLocationCenter, __RC__)

!           ! dry deposition
            call MAPL_AddImportSpec(GC,                                       &
              short_name = 'climDUDP'//trim(field_name),                      &
              long_name  = 'Dust Mixing Ratio (bin '//trim(field_name)//')',  &
              units      = 'kg kg-1',                                         &
              dims       = MAPL_DimsHorzOnly,                                 &
              vlocation  = MAPL_VLocationCenter,                              &
              restart    = MAPL_RestartSkip, __RC__)

!           ! wet deposition    
            call MAPL_AddImportSpec(GC,                                       &
               short_name = 'climDUWT'//trim(field_name),                     &
               long_name  = 'Dust Mixing Ratio (bin '//trim(field_name)//')', &
               units      = 'kg kg-1',                                        &
               dims       = MAPL_DimsHorzOnly,                                &
               vlocation  = MAPL_VLocationCenter,                             &
               restart    = MAPL_RestartSkip, __RC__)

!           ! gravitational settling
            call MAPL_AddImportSpec(GC,                                       &
               short_name = 'climDUSD'//trim(field_name),                     &
               long_name  = 'Dust Mixing Ratio (bin '//trim(field_name)//')', &
               units      = 'kg kg-1',                                        &
               dims       = MAPL_DimsHorzOnly,                                &
               vlocation  = MAPL_VLocationCenter,                             &
               restart    = MAPL_RestartSkip, __RC__)

!           convective scavenging
            call MAPL_AddImportSpec(GC,                                       &
               short_name = 'climDUSV'//trim(field_name),                     &
               long_name  = 'Dust Mixing Ratio (bin '//trim(field_name)//')', &
               units      = 'kg kg-1',                                        &
               dims       = MAPL_DimsHorzOnly,                                &
               vlocation  = MAPL_VLocationCenter,                             &
               restart    = MAPL_RestartSkip, __RC__)
        end do
    else

       call MAPL_AddImportSpec(GC,                                       &
          short_name = 'DU_SRC',                        &
          long_name  = 'erod'  ,                                         &
          units      = '1',                                              &
          dims       = MAPL_DimsHorzOnly,                                &
          vlocation  = MAPL_VLocationNone,                               &
          restart    = MAPL_RestartSkip, __RC__)

!      FRACLAKE
!      --------
       call MAPL_AddImportSpec(GC,                           &
          SHORT_NAME = 'FRLAKE',                             &
          LONG_NAME  = 'fraction_of_lake',                   &
          UNITS      = '1',                                  &
          DIMS       = MAPL_DimsHorzOnly,                    &
          VLOCATION  = MAPL_VLocationNone,                   &
          RESTART    = MAPL_RestartSkip,    __RC__)

!      GWETTOP
!      -------
       call MAPL_AddImportSpec(GC,                           &
          SHORT_NAME = 'WET1',                               &
          LONG_NAME  = 'surface_soil_wetness',               &
          UNITS      = '1',                                  &
          DIMS       = MAPL_DimsHorzOnly,                    &
          VLOCATION  = MAPL_VLocationNone, __RC__)

!      LWI
!      ---
       call MAPL_AddImportSpec(GC,                           &
          SHORT_NAME = 'LWI',                                &
          LONG_NAME  = 'land-ocean-ice_mask',                &
          UNITS      = '1',                                  &
          DIMS       = MAPL_DimsHorzOnly,                    &
          VLOCATION  = MAPL_VLocationNone,                   &
                                          __RC__)
!      U10M
!      ----
       call MAPL_AddImportSpec(GC,                           &
          SHORT_NAME = 'U10M',                               &
          LONG_NAME  = '10-meter_eastward_wind',             &
          UNITS      = 'm s-1',                              &
          DIMS       = MAPL_DimsHorzOnly,                    &
          VLOCATION  = MAPL_VLocationNone,                   &
          RESTART    = MAPL_RestartSkip,   __RC__)

!      V10M
!      ----
       call MAPL_AddImportSpec(GC,                           &
          SHORT_NAME = 'V10M',                               &
          LONG_NAME  = '10-meter_northward_wind',            &
          UNITS      = 'm s-1',                              &
          DIMS       = MAPL_DimsHorzOnly,                    &
          VLOCATION  = MAPL_VLocationNone,                   &
          RESTART    = MAPL_RestartSkip,   __RC__)

!      AIRDENS: moist air density
!      --------------------------
       call MAPL_AddImportSpec(GC,                           &
          SHORT_NAME = 'AIRDENS',                            &
          LONG_NAME  = 'moist_air_density',                  &
          UNITS      = 'kg/m^3',                             &
          DIMS       = MAPL_DimsHorzVert,                    &
          VLOCATION  = MAPL_VLocationCenter,                 &
          RESTART    = MAPL_RestartSkip,    __RC__)

!      Cell area
!      ---------
       call MAPL_AddImportSpec(GC,                            &
           SHORT_NAME = 'AREA',                               &
           LONG_NAME  = 'agrid_cell_area',                    &
           UNITS      = 'm^2',                                &
           DIMS       = MAPL_DimsHorzOnly,                    &
           VLOCATION  = MAPL_VLocationNone,                   &
           RESTART    = MAPL_RestartSkip,   __RC__)

    end if ! (data_driven)


!   EXPORT STATE
!   -------------
    if (.not. data_driven) then
#       include "DU2G_ExportSpec___.h"
    end if

!   This state holds fields needed by radiation
!   ---------------------------------------------
    call MAPL_AddExportSpec (GC,                             &
       short_name = trim(COMP_NAME)//'_AERO',                &
       long_name  = 'aerosols_from_'//trim(COMP_NAME),       &
       units      = 'kg kg-1',                               &
       dims       = MAPL_DimsHorzVert,                       &
       vlocation  = MAPL_VLocationCenter,                    &
       datatype   = MAPL_StateItem, __RC__)

!   This state is needed by MOIST - It will contain aerosols
!   ----------------------------------------------------------
    call MAPL_AddExportSpec (GC,                                                  &
       short_name = trim(COMP_NAME)//'_AERO_ACI',                                 &
       long_name  = 'aerosol_cloud_interaction_aerosols_from_'//trim(COMP_NAME),  &
       units      = 'kg kg-1',                                                    &
       dims       = MAPL_DimsHorzVert,                                            &
       vlocation  = MAPL_VLocationCenter,                                         &
       datatype   = MAPL_StateItem, __RC__)

!   This bundle is needed by surface for snow albedo modification
!   by aerosol settling and deposition
!   ~~~DEVELOPERS NOTE~~~ Change to StateItem when possible
!                         This will require refactoring Radiation
!   ---------------------------------------------------------------
    call MAPL_AddExportSpec (GC,                                  &
       short_name = trim(COMP_NAME)//'_AERO_DP',                                    &
       long_name  = 'aerosol_deposition_from_'//trim(COMP_NAME),  &
       units      = 'kg m-2 s-1',                                 &
       dims       = MAPL_DimsHorzOnly,                            &
       datatype   = MAPL_BundleItem, __RC__)



!  Store internal state in GC
!  --------------------------
   call ESMF_UserCompSetInternalState ( GC, 'DU2G_GridComp', wrap, STATUS )
   VERIFY_(STATUS)

!   Set generic services
!   ----------------------------------
    call MAPL_GenericSetServices (GC, __RC__)


    RETURN_(ESMF_SUCCESS)

  end subroutine SetServices

!============================================================================
!BOP

! !IROUTINE: Initialize 

! !INTERFACE:
  subroutine Initialize (GC, import, export, clock, RC)

!   !ARGUMENTS:
    type (ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
    type (ESMF_State),    intent(inout) :: import ! Import state
    type (ESMF_State),    intent(inout) :: export ! Export state
    type (ESMF_Clock),    intent(inout) :: clock  ! The clock
    integer, optional,    intent(  out) :: RC     ! Error code

! !DESCRIPTION: This initializes DU's Grid Component. It primaryily fills 
!               GOCART's AERO states with its dust fields. 

! !REVISION HISTORY: 
! 16oct2019   E.Sherman  First attempt at refactoring

!EOP
!============================================================================
!   !Locals
    character (len=ESMF_MAXSTR)          :: COMP_NAME
    type (MAPL_MetaComp),       pointer  :: MAPL
    type (ESMF_Grid)                     :: grid
    type (ESMF_State)                    :: internal
    type (ESMF_State)                    :: aero, aero_aci
    type (ESMF_State)                    :: providerState
    type (ESMF_Config)                   :: cfg
    type (ESMF_FieldBundle)              :: Bundle_DP
    type (wrap_)                         :: wrap
    type (DU2G_GridComp), pointer        :: self

    integer, allocatable                 :: mieTable_pointer(:)
    integer                              :: i, n_bins, nCols, dims(3)
    integer                              :: instance
    type (ESMF_Field)                    :: field, fld
    character (len=ESMF_MAXSTR)          :: field_name, prefix
    character (len=ESMF_MAXSTR), allocatable           :: aerosol_names(:)

    logical                              :: data_driven
    integer                              :: NUM_BANDS

    __Iam__('Initialize')

!****************************************************************************

!   Begin... 

!   Get the target components name and set-up traceback handle.
!   -----------------------------------------------------------
    call ESMF_GridCompGet (GC, grid=grid, name=COMP_NAME, __RC__)
    Iam = trim(COMP_NAME) // '::' //trim(Iam)


if (mapl_am_I_root()) print*,trim(comp_name),' INIT BEGIN'


!   Get my internal MAPL_Generic state
!   -----------------------------------
    call MAPL_GetObjectFromGC (GC, MAPL, __RC__)

!   Get my internal state
!   ---------------------
    call ESMF_UserCompGetInternalState(GC, 'DU2G_GridComp', wrap, STATUS)
    VERIFY_(STATUS)
    self => wrap%ptr

if (mapl_am_I_root()) print*,'DU2G Ch_DU before = ', self%Ch_DU

    call MAPL_GridGet ( grid, globalCellCountPerDim=dims, __RC__ )

!  Dust emission tuning coefficient [kg s2 m-5]. NOT bin specific.
!  ---------------------------------------------------------------
   self%Ch_DU = Chem_UtilResVal(dims(1), dims(2), self%Ch_DU(:), __RC__)
   self%Ch_DU = self%Ch_DU * 1.00E-09
if (mapl_am_I_root()) print*,'DU2G dims(1) = ', dims(1)
if (mapl_am_I_root()) print*,'DU2G dims(2) = ', dims(2)
if (mapl_am_I_root()) print*,'DU2G Ch_DU after = ', self%Ch_DU


!   Load resource file  
!   -------------------
    cfg = ESMF_ConfigCreate (__RC__)
    call ESMF_ConfigLoadFile (cfg, 'DU2G_GridComp_'//trim(COMP_NAME)//'.rc', RC=STATUS)
    if (STATUS /= 0) then
        if (mapl_am_i_root()) print*,'DU2G_GridComp_'//trim(COMP_NAME)//'.rc does not exist! &
                                      loading DU2G_GridComp_DU.rc instead'
        call ESMF_ConfigLoadFile (cfg, 'DU2G_GridComp_DU.rc', __RC__)
    end if

!   Get names of aerosols and write to list to add to AERO State
!   -------------------------------------------------------------
    call ESMF_ConfigGetDim (cfg, n_bins, nCols, label=('variable_table::'), __RC__ )
    allocate (aerosol_names(n_bins), __STAT__)
    call ESMF_ConfigFindLabel (cfg, 'variable_table::', __RC__ )

    do i = 1, n_bins
        call ESMF_ConfigNextLine (cfg, __RC__ )
        call ESMF_ConfigGetAttribute (cfg, aerosol_names(i), __RC__ )
    end do


!   Call Generic Initialize 
!   ------------------------
    call MAPL_GenericInitialize (GC, import, export, clock, __RC__)

!   Get parameters from generic state.
!   -----------------------------------
    call MAPL_Get (MAPL, INTERNAL_ESMF_STATE=internal, __RC__)

!   Is DU data driven?
!   ------------------
    call data_driven_ (COMP_NAME, data_driven, __RC__)


!   If this is a data component, the data is provided in the import
!   state via ExtData instead of the actual GOCART children
!   ---------------------------------------------------------------- 
    if ( data_driven ) then
        providerState = import
        prefix = 'clim'
    else
        providerState = export
        prefix = ''
    end if

if(mapl_am_i_root()) print*,'DU2G INIT data_driven = ',data_driven


!   Fill AERO States with dust fields
!   ------------------------------------
    call ESMF_StateGet (export, trim(COMP_NAME)//'_AERO'    , aero    , __RC__)
    call ESMF_StateGet (export, trim(COMP_NAME)//'_AERO_ACI', aero_aci, __RC__)
    call ESMF_StateGet (export, trim(COMP_NAME)//'_AERO_DP' , Bundle_DP, __RC__)

!   Add list of aerosol names to AERO state 
    call ESMF_AttributeSet(aero, name='aerosol_names', valueList=aerosol_names, itemCount=size(aerosol_names), __RC__)


!   Fill AERO States and Bundle
!   ------------------------------
    do i = 1, size(aerosol_names)
        write (field_name, '(A, I0.3)') '', i
        call ESMF_StateGet (internal, trim(COMP_NAME)//'::'//trim(aerosol_names(i)), field, __RC__)
        fld = MAPL_FieldCreate (field, name='du'//trim(field_name), __RC__)
        call MAPL_StateAdd (aero    , fld, __RC__)
        call MAPL_StateAdd (aero_aci, fld, __RC__)

!       Dry deposition
!       ---------------
        call ESMF_StateGet (providerState, trim(prefix)//'DUDP'//trim(field_name), field, __RC__)
        call MAPL_AllocateCoupling (field, __RC__)
        call MAPL_FieldBundleAdd (Bundle_DP, field, __RC__)

!        Wet deposition (Convective scavenging)
!        --------------------------------------
        call ESMF_StateGet (providerState, trim(prefix)//'DUSV'//trim(field_name), field, __RC__)
        call MAPL_AllocateCoupling (field, __RC__)
        call MAPL_FieldBundleAdd (Bundle_DP, field, __RC__)

!        Wet deposition
!        ---------------
        call ESMF_StateGet (providerState, trim(prefix)//'DUWT'//trim(field_name), field, __RC__)
        call MAPL_AllocateCoupling (field, __RC__)
        call MAPL_FieldBundleAdd (Bundle_DP, field, __RC__)

!        Gravitational Settling
!        ----------------------
        call ESMF_StateGet (providerState, trim(prefix)//'DUSD'//trim(field_name), field, __RC__)
        call MAPL_AllocateCoupling (field, __RC__)
        call MAPL_FieldBundleAdd (Bundle_DP, field, __RC__)
    end do

!   Set AERO States' attributes
!   ----------------------------
    if (data_driven) then
        instance = instanceData
    else
        instance = instanceComputational
    end if


    call MAPL_GetResource (MAPL, NUM_BANDS, 'NUM_BANDS:', __RC__)
    self%rad_MieTable(instance) = Chem_MieCreateng (cfg, NUM_BANDS,__RC__)

    ! Mie Table instance/index
    call ESMF_AttributeSet (aero, name='mie_table_instance', value=instance, __RC__)

    ! Add variables to DU instance's AERO state. This is used in aerosol optics calculations
    call addAero (aero, label='air_pressure_for_aerosol_optics',             label2='PLE', grid=grid, typekind=MAPL_R4, __RC__)
    call addAero (aero, label='relative_humidity_for_aerosol_optics',        label2='RH',  grid=grid, typekind=MAPL_R4, __RC__)
!   call ESMF_StateGet (import, 'PLE', field, __RC__)
!   call MAPL_StateAdd (aero, field, __RC__)

!   call ESMF_StateGet (import, 'RH2', field, __RC__)
!   call MAPL_StateAdd (aero, field, __RC__)


    call addAero (aero, label='extinction_in_air_due_to_ambient_aerosol',    label2='EXT', grid=grid, typekind=MAPL_R8, __RC__)
    call addAero (aero, label='single_scattering_albedo_of_ambient_aerosol', label2='SSA', grid=grid, typekind=MAPL_R8, __RC__)
    call addAero (aero, label='asymmetry_parameter_of_ambient_aerosol',      label2='ASY', grid=grid, typekind=MAPL_R8, __RC__)

    call ESMF_AttributeSet(aero, name='band_for_aerosol_optics',             value=0,     __RC__)

!    mieTable_pointer = transfer(c_loc(DU2G_GridComp), [1])
    mieTable_pointer = transfer(c_loc(self), [1])

    call ESMF_AttributeSet(aero, name='mieTable_pointer', valueList=mieTable_pointer, itemCount=size(mieTable_pointer), __RC__)


    ! attach the aerosol optics method
    call ESMF_MethodAdd (aero, label='aerosol_optics', userRoutine=aerosol_optics, __RC__)

if (mapl_am_I_root()) print*,trim(comp_name),' INIT END'

    RETURN_(ESMF_SUCCESS)

  end subroutine Initialize


!============================================================================
!BOP
! !IROUTINE: Run 

! !INTERFACE:
  subroutine Run (GC, import, export, clock, RC)

!   !ARGUMENTS:
    type (ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
    type (ESMF_State),    intent(inout) :: import ! Import state
    type (ESMF_State),    intent(inout) :: export ! Export state
    type (ESMF_Clock),    intent(inout) :: clock  ! The clock
    integer, optional,    intent(  out) :: RC     ! Error code:

!   !DESCRIPTION: Run method for the Dust Grid Component. Determines whether to
!                 run data or computational run method.

!EOP
!============================================================================
!   !Locals
    character (len=ESMF_MAXSTR)       :: COMP_NAME
    type (MAPL_MetaComp), pointer     :: MAPL
    type (ESMF_State)                 :: internal

    logical                           :: data_driven

    __Iam__('Run')

!*****************************************************************************
!   Begin... 

!   Get my name and set-up traceback handle
!   ---------------------------------------
    call ESMF_GridCompGet (GC, NAME=COMP_NAME, __RC__)
    Iam = trim(COMP_NAME) //'::'// Iam

if (mapl_am_I_root()) print*,trim(comp_name),' Run BEGIN'


!   Get my internal MAPL_Generic state
!   -----------------------------------
    call MAPL_GetObjectFromGC (GC, MAPL, __RC__)

!   Get parameters from generic state.
!   -----------------------------------
    call MAPL_Get (MAPL, INTERNAL_ESMF_STATE=internal, __RC__)

!   Is DU data driven?
!   ------------------
    call data_driven_ (COMP_NAME, data_driven, __RC__)

!   Update INTERNAL state variables with ExtData
!   ---------------------------------------------
    if (data_driven) then
        call Run_data (GC, import, export, internal, __RC__)
    else
        call Run1 (GC, import, export, clock, __RC__)
    end if

if (mapl_am_I_root()) print*,trim(comp_name),' Run END'

    RETURN_(ESMF_SUCCESS)

  end subroutine Run




!============================================================================
!BOP
! !IROUTINE: Run1 

! !INTERFACE:
  subroutine Run1 (GC, import, export, clock, RC)

!   !ARGUMENTS:
    type (ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
    type (ESMF_State),    intent(inout) :: import ! Import state
    type (ESMF_State),    intent(inout) :: export ! Export state
    type (ESMF_Clock),    intent(inout) :: clock  ! The clock
    integer, optional,    intent(  out) :: RC     ! Error code:

! !DESCRIPTION:  Computes emissions/sources for Dust

!EOP
!============================================================================
!   !Locals
    character (len=ESMF_MAXSTR)       :: COMP_NAME
    type (MAPL_MetaComp), pointer     :: mapl
    type (ESMF_State)                 :: internal
    type (ESMF_Grid)                  :: grid
    type (wrap_)                      :: wrap
    type (DU2G_GridComp), pointer     :: self
    type(ESMF_Time)                   :: time

    logical                           :: data_driven
    integer                           :: n_bins, dims(3), import_shape(2)
    integer                           :: i1=1, i2   ! dist grid indices
    integer                           :: j1=1, j2   ! dist grid indices
    integer                           :: km, nq              ! dist grid indices
    real                              :: CDT         ! chemistry timestep (secs)
    integer                           :: HDT         ! model     timestep (secs)

    real, pointer, dimension(:,:,:)   :: ptr3d_int, delp, rhoa
    real(ESMF_KIND_R4), pointer, dimension(:,:) :: cell_area
    real, pointer, dimension(:,:)     :: fraclake, gwettop, oro, u10m, v10m, du_src
    real, pointer                     :: emissions(:,:,:), dqa(:,:,:)
    real, allocatable                 :: radius(:)

    integer                           :: nymd, nhms, iyr, imm, idd, ihr, imn, isc

!   !Indices for point emissions
    integer, pointer, dimension(:)    :: iPoint, jPoint
    real, dimension(:), allocatable   :: point_column_emissions

    integer :: n, i, j, ii

    !REPLACE undef and GRAV WITH MAPL
    real, parameter ::  UNDEF  = 1.e15
    real, parameter ::  GRAV   = 9.80616 ! USE MAPL_GRAV????


   type(MAPL_VarSpec), pointer     :: InternalSpec(:)
   character(len=ESMF_MAXSTR)      :: short_name

   __Iam__('Run1')

!*****************************************************************************
!   Begin... 


!#include "DU_GetPointer___.h"


!   Get my name and set-up traceback handle
!   ---------------------------------------
    call ESMF_GridCompGet (GC, grid=grid, NAME=COMP_NAME, __RC__)
    Iam = trim(COMP_NAME) //'::'// Iam

if (mapl_am_I_root()) print*,trim(comp_name),' Run1 BEGIN'

!   Get my internal MAPL_Generic state
!   -----------------------------------
    call MAPL_GetObjectFromGC (GC, mapl, __RC__)

!   Get parameters from generic state.
!   -----------------------------------
    call MAPL_Get (mapl, INTERNAL_ESMF_STATE=internal, INTERNALSPEC = InternalSpec, __RC__)

!   Get my private internal state
!   ------------------------------
    call ESMF_UserCompGetInternalState(GC, 'DU2G_GridComp', wrap, STATUS)
    VERIFY_(STATUS)
    self => wrap%ptr

!   Get DTs
!   -------
    call MAPL_GetResource(mapl, HDT, Label='RUN_DT:', __RC__)
    call MAPL_GetResource(mapl, CDT, Label='GOCART_DT:', default=real(HDT), __RC__)

!   Get 2D Imports
!   --------------
    call MAPL_GetPointer (import, fraclake, 'FRLAKE', __RC__)
    call MAPL_GetPointer (import, gwettop,  'WET1',   __RC__)
    call MAPL_GetPointer (import, oro,      'LWI',    __RC__)
    call MAPL_GetPointer (import, u10m,     'U10M',   __RC__)
    call MAPL_GetPointer (import, v10m,     'V10M',   __RC__)
    call MAPL_GetPointer (import, delp,     'DELP',   __RC__)
    call MAPL_GetPointer (import, cell_area,'AREA',   __RC__)
    call MAPL_GetPointer (import, du_src,   'DU_SRC', __RC__)

!   Set du_src to 0 where undefined
!   --------------------------------
    where (1.01*du_src > UNDEF) du_src = 0.

!   Get 3D Imports
!   --------------
    call MAPL_GetPointer (import, rhoa, 'AIRDENS',  __RC__ )

!   Get dimensions
!   ---------------
    call MAPL_GridGet (grid, globalCellCountPerDim=dims, __RC__)
    km = dims(3)

    import_shape = shape(gwettop)
    i2 = import_shape(1)
    j2 = import_shape(2)

    n_bins = size(self%radius)

!   Extract nymd(yyyymmdd) from clock
!   ---------------------------------
    call ESMF_ClockGet (clock, currTime=time, __RC__)
    call ESMF_TimeGet (time ,YY=iyr, MM=imm, DD=idd, H=ihr, M=imn, S=isc, __RC__)
    call MAPL_PackTime (nymd, iyr, imm , idd)
    call MAPL_PackTime (nhms, ihr, imn, isc)

!   Dust particle radius [m] and density [kg m-3]
!   ---------------------------------------------
    allocate(radius(n_bins), __STAT__ )
    radius = 1.e-6 * self%radius
!   DU_rhop   = self%rhop
    allocate(emissions(i1:i2,j1:j2,n_bins), dqa(i1:i2,j1:j2,n_bins), __STAT__)

!   Implement gridded emission dust source
!   --------------------------------------
    emissions = 0.0
    dqa = 0.0

    call DustEmissionGOCART2G(radius,fraclake, gwettop, oro, u10m, v10m, &
                              self%Ch_DU(1), self%sfrac, du_src, GRAV, &
                              emissions, rc )

    do i = 1, n_bins
        dqa(:,:,i) = emissions(:,:,i) * CDT * GRAV / delp(:,:,km)

        call MAPL_VarSpecGet(InternalSpec(i), SHORT_NAME=short_name, __RC__)
        call MAPL_GetPointer(internal, NAME=short_name, ptr=ptr3d_int, __RC__)
!       ! update internal pointer with emission
        ptr3d_int(:,:,km) = ptr3d_int(:,:,km) + dqa(:,:,i)
    end do

!   Read pointwise emissions, if requested
!   ---------------------------------------
    if(self%doing_point_emissions) then
        call Chem_UtilPointEmissions( nymd, self%point_emissions_srcfilen, &
                                      self%nPts, self%pLat, self%pLon, &
                                      self%pBase, self%pTop, self%pEmis, &
                                      self%pStart, self%pEnd )

!   In case pStart or pEnd were not specified in the file set to defaults
        where(self%pStart < 0) self%pStart = 000000
        where(self%pEnd < 0)   self%pEnd   = 240000
    endif

!   Distribute pointwise sources if requested
!   -----------------------------------------
    POINTWISE_SOURCES: if( self%doing_point_emissions .and. self%nPts > 0) then

!   Get indices for point emissions
!   -------------------------------
    allocate(iPoint(self%nPts), jPoint(self%nPts), point_column_emissions(km), __STAT__)

! DEV NOTE - radToDeg is a defined parameter. Is there a MAPL equivalent?
    call MAPL_GetHorzIJIndex(self%nPts, iPoint, jPoint, &
                             grid = grid,               &
                             lon  = self%pLon/radToDeg, &
                             lat  = self%pLat/radToDeg, &
                             rc   = rc)
    if ( rc /= 0 ) then
        if (mapl_am_i_root()) print*, trim(Iam), ' - cannot get indices for point emissions'
        VERIFY_(rc)
    end if

    do ii = 1, self%nPts
        i = iPoint(ii)
        j = jPoint(ii)

!       point emission not in this sub-domain
        if( i<1 .or. j<1 ) cycle
!       Emissions not occurring in current time step
        if(nhms < self%pStart(ii) .or. nhms >= self%pEnd(ii)) cycle

           call DistributePointEmission(km, delp(i,j,:), rhoa(i,j,:), self%pBase(ii), & 
                                        self%pTop(ii), GRAV, self%pEmis(ii), &
                                        point_column_emissions, rc)

           do n = 1, n_bins
               call MAPL_VarSpecGet(InternalSpec(n), SHORT_NAME=short_name, __RC__)
               call MAPL_GetPointer(internal, NAME=short_name, ptr=ptr3d_int, __RC__)
!              ! update internal pointer with emission
               ptr3d_int(i,j,:) = ptr3d_int(i,j,:) + cdt * grav / delp(i,j,:) * &
               self%sfrac(n) * point_column_emissions / cell_area(i,j) 
           enddo ! do n
    enddo ! do ii

    deallocate(iPoint, jPoint, __STAT__)

    endif POINTWISE_SOURCES

    RETURN_(ESMF_SUCCESS)

  end subroutine Run1

!============================================================================
!BOP
! !IROUTINE: Run2 

! !INTERFACE:

  subroutine Run2 (GC, import, export, clock, RC)

    ! !ARGUMENTS:

    type (ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
    type (ESMF_State),    intent(inout) :: import ! Import state
    type (ESMF_State),    intent(inout) :: export ! Export state
    type (ESMF_Clock),    intent(inout) :: clock  ! The clock
    integer, optional,    intent(  out) :: RC     ! Error code:

! !DESCRIPTION: Run2 method for the Dust Grid Component.

!EOP
!============================================================================
! Locals
    character (len=ESMF_MAXSTR)       :: COMP_NAME
    type (MAPL_MetaComp), pointer     :: MAPL
    type (ESMF_State)                 :: internal
    logical                           :: data_driven

    __Iam__('Run2')

!*****************************************************************************
!   Begin... 

!   Get my name and set-up traceback handle
!   ---------------------------------------
    call ESMF_GridCompGet (GC, NAME=COMP_NAME, __RC__)
    Iam = trim(COMP_NAME) // '::' // Iam


!   Get my internal MAPL_Generic state
!   -----------------------------------
    call MAPL_GetObjectFromGC (GC, MAPL, __RC__)

!   Get parameters from generic state.
!   -----------------------------------
!    call MAPL_Get (MAPL, INTERNAL_ESMF_STATE=internal, __RC__)



    RETURN_(ESMF_SUCCESS)

  end subroutine Run2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!BOP
! !IROUTINE: Run_data -- ExtData Dust Grid Component

! !INTERFACE:

  subroutine Run_data (GC, import, export, internal, RC)

    ! !ARGUMENTS:

    type (ESMF_GridComp), intent(inout) :: GC       ! Gridded component 
    type (ESMF_State),    intent(inout) :: import   ! Import state
    type (ESMF_State),    intent(inout) :: export   ! Export state
    type (ESMF_State),    intent(inout) :: internal ! Interal state
    integer, optional,    intent(  out) :: RC       ! Error code:

! !DESCRIPTION: Updates pointers in Internal state with fields from ExtData. 

!EOP

!****************************************************************************
! Locals
    character (len=ESMF_MAXSTR)        :: COMP_NAME
    type (MAPL_MetaComp), pointer      :: MAPL
    type (ESMF_State)                  :: aero 

    integer                            :: i, n
    character (len=ESMF_MAXSTR)        :: field_name

    real, pointer, dimension(:,:,:)    :: ptr3d_int
    real, pointer, dimension(:,:,:)    :: ptr3d_imp

    __Iam__('Run_data')


!*****************************************************************************
!   Begin... 

!   Get my name and set-up traceback handle
!   ---------------------------------------
    call ESMF_GridCompGet (GC, NAME=COMP_NAME, __RC__)

if (mapl_am_I_root()) print*,trim(comp_name),' Run_data BEGIN'


    call ESMF_StateGet (export, trim(COMP_NAME)//'_AERO', aero, __RC__)
    call ESMF_AttributeGet (aero, name='aerosol_names', itemCount=n, __RC__)

!   Update interal data pointers with ExtData
!   -----------------------------------------
    do i = 1, n
    write (field_name, '(A, I0.3)') 'du', i
        call MAPL_GetPointer (internal, NAME=trim(COMP_NAME)//'::'//trim(field_name), ptr=ptr3d_int, __RC__)
        call MAPL_GetPointer (import,  NAME='clim'//trim(field_name), ptr=ptr3d_imp, __RC__)

        ptr3d_int = ptr3d_imp
    end do


if (mapl_am_I_root()) print*,trim(comp_name),' Run_data END'


    RETURN_(ESMF_SUCCESS)

  end subroutine Run_data


!-------------------------------------------------------------------------------------
  subroutine aerosol_optics(state, rc)

    implicit none

!   !ARGUMENTS:
    type (ESMF_State)                                :: state
    integer,            intent(out)                  :: rc

!   !Local
    integer, parameter                               :: DP=kind(1.0d0)
    real, dimension(:,:,:), pointer                  :: ple, rh
    real(kind=DP), dimension(:,:,:), pointer         :: var
    real, dimension(:,:,:), pointer                  :: q
    real, dimension(:,:,:,:), pointer                :: q_4d
    integer, allocatable                             :: opaque_self(:)
    type(C_PTR)                                      :: address
    type(DU2G_GridComp), pointer                     :: self

    character (len=ESMF_MAXSTR)                      :: fld_name
    type(ESMF_Field)                                 :: fld
    character (len=ESMF_MAXSTR)                      :: COMP_NAME
    character (len=ESMF_MAXSTR), allocatable         :: aerosol_names(:)

    real(kind=DP), dimension(:,:,:), allocatable     :: ext_s, ssa_s, asy_s  ! (lon:,lat:,lev:)

    integer                                          :: instance
    integer                                          :: n
    integer                                          :: i1, j1, i2, j2, km
    integer                                          :: band, offset
    integer, parameter                               :: n_bands = 1

    real    :: x
    integer :: i, j, k

    __Iam__('DU2G::aerosol_optics')


!   Begin... 

!   Mie Table instance/index
!   ------------------------
    call ESMF_AttributeGet (state, name='mie_table_instance', value=instance, __RC__)

!   Get aerosol names in this state
!   -------------------------------
    call ESMF_AttributeGet (state, name='aerosol_names', itemCount=n, __RC__)
    allocate (aerosol_names(n), __STAT__)
    call ESMF_AttributeGet (state, name='aerosol_names', valueList=aerosol_names, __RC__)


!   Radiation band
!   --------------
    band = 0
    call ESMF_AttributeGet (state, name='band_for_aerosol_optics', value=band, __RC__)
    offset = band - n_bands


!   Pressure at layer edges 
!   ------------------------
    call ESMF_AttributeGet (state, name='air_pressure_for_aerosol_optics', value=fld_name, __RC__)
    call MAPL_GetPointer (state, ple, trim(fld_name), __RC__)

!    call MAPL_GetPointer (state, ple, 'PLE', __RC__)

    i1 = lbound(ple, 1); i2 = ubound(ple, 1)
    j1 = lbound(ple, 2); j2 = ubound(ple, 2)
                         km = ubound(ple, 3)


!   Relative humidity
!   -----------------
    call ESMF_AttributeGet (state, name='relative_humidity_for_aerosol_optics', value=fld_name, __RC__)
    call MAPL_GetPointer (state, rh, trim(fld_name), __RC__)

!    call MAPL_GetPointer (state, rh, 'RH2', __RC__)


    allocate(ext_s(i1:i2, j1:j2, km), &
             ssa_s(i1:i2, j1:j2, km), &
             asy_s(i1:i2, j1:j2, km), __STAT__)

    allocate(q_4d(i1:i2, j1:j2, km, size(aerosol_names)), __STAT__)

    do n = 1, size(aerosol_names)
        call ESMF_StateGet (state, trim(aerosol_names(n)), field=fld, __RC__)
        call ESMF_FieldGet (fld, farrayPtr=q, __RC__)

        do k = 1, km
            do j = j1, j2
                do i = i1, i2
                    x = ((PLE(i,j,k) - PLE(i,j,k-1))*0.01)*(100./MAPL_GRAV)
                    q_4d(i,j,k,n) = x * q(i,j,k)
                end do
            end do
        end do
    end do

!if(mapl_am_i_root()) print*,'DU2G CALL BACK TEST'

    call ESMF_AttributeGet(state, name='mieTable_pointer', itemCount=n, __RC__)
    allocate (opaque_self(n), __STAT__)
    call ESMF_AttributeGet(state, name='mieTable_pointer', valueList=opaque_self, __RC__)

    address = transfer(opaque_self, address)
    call c_f_pointer(address, self)


!    call mie_ (DU2G_GridComp%rad_MieTable(instance), aerosol_names, n_bands, offset, q_4d, rh, ext_s, ssa_s, asy_s, __RC__)
    call mie_ (self%rad_MieTable(instance), aerosol_names, n_bands, offset, q_4d, rh, ext_s, ssa_s, asy_s, __RC__)


    call ESMF_AttributeGet (state, name='extinction_in_air_due_to_ambient_aerosol', value=fld_name, __RC__)
    if (fld_name /= '') then
        call MAPL_GetPointer (state, var, trim(fld_name), __RC__)
        var = ext_s(:,:,:)
    end if

    call ESMF_AttributeGet (state, name='single_scattering_albedo_of_ambient_aerosol', value=fld_name, __RC__)
    if (fld_name /= '') then
        call MAPL_GetPointer (state, var, trim(fld_name), __RC__)
        var = ssa_s(:,:,:)
    end if

    call ESMF_AttributeGet (state, name='asymmetry_parameter_of_ambient_aerosol', value=fld_name, __RC__)
    if (fld_name /= '') then
        call MAPL_GetPointer (state, var, trim(fld_name), __RC__)
        var = asy_s(:,:,:)
    end if


    deallocate(ext_s, ssa_s, asy_s, __STAT__)


    RETURN_(ESMF_SUCCESS)

  contains

!    subroutine mie_(mie_table, aerosol_names, nb, offset, q, rh, bext_, bssa_, gasym_, rc)
    subroutine mie_(mie_table, aerosol_names, nb, offset, q, rh, bext_s, bssa_s, basym_s, rc)

    implicit none

    type(Chem_Mie),                intent(inout) :: mie_table        ! mie table
    character(len=ESMF_MAXSTR),    intent(in )   :: aerosol_names(:) ! names of aerosols
    integer,                       intent(in )   :: nb               ! number of bands
    integer,                       intent(in )   :: offset           ! bands offset 
    real,                          intent(in )   :: q(:,:,:,:)       ! aerosol mass mixing ratio, kg kg-1
    real,                          intent(in )   :: rh(:,:,:)        ! relative humidity
    real(kind=DP), intent(  out) :: bext_s (size(ext_s,1),size(ext_s,2),size(ext_s,3))
    real(kind=DP), intent(  out) :: bssa_s (size(ext_s,1),size(ext_s,2),size(ext_s,3))
    real(kind=DP), intent(  out) :: basym_s(size(ext_s,1),size(ext_s,2),size(ext_s,3))
    integer,                       intent(  out) :: rc

    ! local
    integer                           :: l, n_bins, idx
    real                              :: bext (size(ext_s,1),size(ext_s,2),size(ext_s,3))  ! extinction
    real                              :: bssa (size(ext_s,1),size(ext_s,2),size(ext_s,3))  ! SSA
    real                              :: gasym(size(ext_s,1),size(ext_s,2),size(ext_s,3))  ! asymmetry parameter

    __Iam__('DU2G::aerosol_optics::mie_')


     n_bins = size(aerosol_names)

     ASSERT_ (n_bins == size(q,4))

     bext_s  = 0.0d0
     bssa_s  = 0.0d0
     basym_s = 0.0d0

     do l = 1, n_bins
         idx = Chem_MieQueryIdx(mie_table, trim(aerosol_names(l)), __RC__)
         call Chem_MieQuery(mie_table, idx, real(offset+1.), q(:,:,:,l), rh, bext, gasym=gasym, ssa=bssa)

         bext_s  = bext_s  +             bext     ! extinction
         bssa_s  = bssa_s  +       (bssa*bext)    ! scattering extinction
         basym_s = basym_s + gasym*(bssa*bext)    ! asymetry parameter multiplied by scatering extiction 
     end do


     RETURN_(ESMF_SUCCESS)

    end subroutine mie_

  end subroutine aerosol_optics



!---------------------------------------------------------------------------------------
  subroutine data_driven_(COMP_NAME, data_driven, RC)

    !ARGUMENTS:
    integer, optional,               intent(  out)   :: RC          ! Error code:
    character (len=ESMF_MAXSTR),     intent(in   )   :: COMP_NAME
    logical,                         intent(  out)   :: data_driven

    !Local
    integer                                          :: i

!   Description: Determines whether Dust is data driven or not.

     __Iam__('data_driven_')

!   Begin... 

!   Is DU data driven?
!   ------------------
    data_driven = .false.

    i = index(COMP_NAME, 'data')
    if (i > 0) then
      data_driven = .true.
    end if

    RETURN_(ESMF_SUCCESS)

  end subroutine data_driven_

!-------------------------------------------------------------------------------------


end module DU2G_GridCompMod

