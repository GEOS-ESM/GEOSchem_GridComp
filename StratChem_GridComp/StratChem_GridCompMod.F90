#include "MAPL_Generic.h"

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 910.1     !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: StratChem_GridCompMod - The StratChem Grid Component
!
! !INTERFACE:
!
   MODULE StratChem_GridCompMod
!
! !USES:
!
   USE ESMF
   USE MAPL
   USE Chem_Mod 	                        ! Chemistry Base Class
   USE SC_GridCompMod                           ! ESMF parent component


   IMPLICIT NONE
   PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:

   PUBLIC SetServices
!
! !DESCRIPTION: 
!
!  {\tt StratChem} is a ESMF gridded component implemented Code 613.3
!  Stratospheric Chemistry package.
!  This code derives from the pre-ESMF SC component from GEOS-4. This
!  GEOS-4 "component" used ESMF like constructs (Chem component class, 
!  import/export states, etc) but no ESMF specific data types because of 
!  an odd incompatibility with the fvGCM code (the so-called 
!  {\tt oldworld} library. Unlike GEOS-4, here the Stratospheric Chemistry
!  component is treated separately.
!
! !REVISION HISTORY:
!
!  25feb2005  da Silva  First crack.
!  19jul2006  da Silva  First separate StratChem component.
!
!EOP
!-------------------------------------------------------------------------

  type StratChem_State
     private
     type(Chem_Registry), pointer :: chemReg => null()
     type(SC_GridComp), pointer :: gcChem => null()
     type(Chem_Bundle), pointer   :: w_c    => null()
   end type StratChem_State

  type StratChem_WRAP
     type (StratChem_State), pointer :: PTR => null()
  end type StratChem_WRAP

CONTAINS


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 910.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: SetServices --- Sets IRF services for StratChem Grid Component
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
!  25feb2005  da Silva  First crack.
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
    type (ESMF_Config)               :: CF
    type (StratChem_State), pointer  :: state   ! internal, that is
    type (StratChem_wrap)            :: wrap

    integer                          :: n, nq, i_XX, j_XX
    CHARACTER(LEN=ESMF_MAXSTR)       :: FRIENDLIES, providerName

!                              ------------


!   Get my name and set-up traceback handle
!   ---------------------------------------
    call ESMF_GridCompGet( GC, NAME=COMP_NAME, CONFIG=CF, RC=STATUS )
    VERIFY_(STATUS)
    Iam = TRIM(COMP_NAME)//"::SetServices"

!   Wrap internal state for storing in GC. Rename legacy state
!   ----------------------------------------------------------
    allocate ( state, stat=STATUS )
    VERIFY_(STATUS)
    wrap%ptr => state

!   Start by loading the Chem Registry
!   ----------------------------------
    allocate ( state%chemReg )
    state%chemReg = Chem_RegistryCreate ( STATUS )
    VERIFY_(STATUS)

    nq = state%chemReg%nq     ! total number of chemical tracers
    


!                       ------------------------
!                       ESMF Functional Services
!                       ------------------------


 if ( state%chemReg%doing_SC ) then

    IF(MAPL_AM_I_ROOT()) THEN
     PRINT *, TRIM(Iam)//': ACTIVE'
     CALL Chem_RegistryPrint ( state%chemReg )
    END IF


!   Set the Initialize, Run, Finalize entry points
!   ----------------------------------------------

    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_INITIALIZE,  Initialize_, &
                                      RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_GridCompSetEntryPoint ( GC,  ESMF_METHOD_RUN,  Run_,        &
                                      RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_GridCompSetEntryPoint ( GC,  ESMF_METHOD_FINALIZE,  Finalize_,  &
                                      RC=STATUS)
    VERIFY_(STATUS)


!   Store internal state in GC
!   --------------------------
    call ESMF_UserCompSetInternalState ( GC, 'StratChem_state', wrap, STATUS )
    VERIFY_(STATUS)


 else

        if (MAPL_AM_I_ROOT()) &
        print *, trim(Iam)//': NOT ACTIVE.  Defaulting to generic no-op stubs.'

 endif

! ======================== IMPORT STATE ===========================

#include "StratChem_ImportSpec___.h"

!  Imported aerosol surface area for heterogeneous chemistry
!  provided if using an interactive aerosol code
!  So: in StratChem_ExtData.rc these are provided as /dev/null
!      but if GOCART running and appropriate tracers running 
!      these are provided in GEOS_ChemGridComp (alternatively,
!      if CARMA is aero provider it will provide these)
   CALL MAPL_AddImportSpec(GC,  &
         SHORT_NAME         = 'SO4SAREAvolc',  &
         LONG_NAME          = 'SO4 aerosol surface area (Volcanic)',  &
         UNITS              = 'm2 m-3', &
         DIMS               = MAPL_DimsHorzVert,    &
         VLOCATION          = MAPL_VLocationCenter,    &
                                                        RC=STATUS  )
   VERIFY_(STATUS)
   CALL MAPL_AddImportSpec(GC,  &
         SHORT_NAME         = 'SO4SAREA',  &
         LONG_NAME          = 'SO4 aerosol surface area (non-Volcanic)',  &
         UNITS              = 'm2 m-3', &
         DIMS               = MAPL_DimsHorzVert,    &
         VLOCATION          = MAPL_VLocationCenter,    &
                                                        RC=STATUS  )
   VERIFY_(STATUS)

! ======================== INTERNAL STATE =========================


! Is STRATCHEM providing O3, in Ox's stead, to the ANALYSIS bundle?
! -----------------------------------------------------------------
     CALL ESMF_ConfigGetAttribute(CF, providerName, Default="PCHEM", &
                                  Label="ANALYSIS_OX_PROVIDER:", RC=STATUS )
     VERIFY_(STATUS)

!   Species to be transported
!   -------------------------
    do n = state%chemReg%i_SC, state%chemReg%j_SC 

	  IF(TRIM(state%chemReg%vname(n)) == "OX" .AND. TRIM(providerName) == "STRATCHEM") THEN
           FRIENDLIES="ANALYSIS:DYNAMICS:TURBULENCE:MOIST"
          ELSE IF(TRIM(ESMF_UtilStringUpperCase(state%chemReg%vname(n))) == "HBR"  .OR. &
                  TRIM(ESMF_UtilStringUpperCase(state%chemReg%vname(n))) == "HOBR" .OR. &
                  TRIM(ESMF_UtilStringUpperCase(state%chemReg%vname(n))) == "BRONO2" ) THEN
           FRIENDLIES="DYNAMICS:TURBULENCE"
          ELSE
           FRIENDLIES="DYNAMICS:TURBULENCE:MOIST"
	  END IF

          FRIENDLIES = FRIENDLIES//':'//TRIM(COMP_NAME)
	 
          call MAPL_AddInternalSpec(GC,                          &
               SHORT_NAME = TRIM(state%chemReg%vname(n)),        &
               LONG_NAME  = TRIM(state%chemReg%vtitle(n)),       &
               UNITS      = TRIM(state%chemReg%vunits(n)),       &
               FRIENDLYTO = FRIENDLIES,                          &
               DIMS       = MAPL_DimsHorzVert,                   &
               VLOCATION  = MAPL_VLocationCenter,     RC=STATUS  )
          VERIFY_(STATUS)

    end do

!   Non-transported species
!   -----------------------
    do n = state%chemReg%i_XX, state%chemReg%j_XX 

          call MAPL_AddInternalSpec(GC,                          &
               SHORT_NAME = TRIM(state%chemReg%vname(n)),	 &
               LONG_NAME  = TRIM(state%chemReg%vtitle(n)),	 &
               UNITS      = TRIM(state%chemReg%vunits(n)),	 &
               FRIENDLYTO = TRIM(COMP_NAME),                     &
               ADD2EXPORT = .TRUE.,				 &
               DIMS       = MAPL_DimsHorzVert,  		 &
               VLOCATION  = MAPL_VLocationCenter,     RC=STATUS  )
          VERIFY_(STATUS)

    end do

! ======================== EXPORT STATE ===========================

    CALL ESMF_ConfigGetAttribute(CF, providerName, Default="GOCART.data", Label="AERO_PROVIDER:", RC=STATUS )
    VERIFY_(STATUS)

    IF(TRIM(providerName) == "STRATCHEM") THEN

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

! ----------------------------------------------------------------
! To conform with GEOS-5's other chemistry components, STRATCHEM's
! ozone exports are as follows:
!
!    NAME   STATE  Units             Comments
!  ------ -------- ----------------- ----------------------
!      OX Internal mol/mol           Formerly named OXSTRAT
!      O3 Export   kg/kg             O3CHEM(vmr)*48/28.97
!  03PPMV Export   ppmv              O3CHEM(vmr)*1.00E+06
!
! ----------------------------------------------------------------
#include "StratChem_ExportSpec___.h"

! =================================================================

!   Set the profiling timers
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
! !IROUTINE:  Initialize_ --- Initialize StratChem
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
   type(SC_GridComp), pointer      :: gcChem      ! Grid Component
   type(Chem_Bundle), pointer      :: w_c         ! Chemical tracer fields     
   integer                         :: nymd, nhms  ! time
   real                            :: cdt         ! chemistry timestep (secs)
   real                            :: rdt         ! run (heartbeat) timestep (secs)

   type(ESMF_Grid)                 :: grid
 
   integer                         :: i1=1, i2, ig=0, im  ! dist grid indices
   integer                         :: j1=1, j2, jg=0, jm  ! dist grid indices
   integer                         :: km, nq              ! dist grid indices
   integer                         :: k, n, dims(3), l, split

   type(Chem_Array), pointer       :: q(:)         ! array of pointers
   type(MAPL_MetaComp), pointer    :: ggState	   ! GEOS Generic State
   type (ESMF_State)               :: internal
   type(MAPL_VarSpec), pointer     :: InternalSpec(:)

   real, pointer, dimension(:,:)   :: LATS
   real, pointer, dimension(:,:)   :: LONS
   real, pointer, dimension(:,:)   :: AREA
   REAL, POINTER, DIMENSION(:,:)   :: depo

   TYPE(ESMF_Field)                :: field
   TYPE(ESMF_Config)               :: CF
   TYPE(ESMF_FieldBundle)          :: bundle
   CHARACTER(LEN=ESMF_MAXSTR)      :: short_name
   CHARACTER(LEN=ESMF_MAXSTR)      :: providerName
   CHARACTER(LEN=ESMF_MAXSTR), POINTER, DIMENSION(:) :: fieldNames

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
   CALL MAPL_TimerOn(ggState, "TOTAL")
   CALL MAPL_TimerOn(ggState, "INITIALIZE")
   
!  Initialize GEOS Generic
!  ------------------------
   call MAPL_GenericInitialize ( gc, impChem, expChem, clock,  RC=STATUS )
   VERIFY_(STATUS)

!  Get pre-ESMF parameters from gc and clock
!  -----------------------------------------
   call extract_ ( gc, clock, chemReg, gcChem, w_c, nymd, nhms, rdt, split, STATUS )
   VERIFY_(STATUS)

!  Chemistry time step length (seconds)
!  ------------------------------------
   STATUS = 0
   IF(split < 1) THEN
    PRINT *,TRIM(Iam)//": SC_SPLIT cannot be less than 1"
    VERIFY_(STATUS)
   END IF
   cdt = rdt/split
   IF(MAPL_AM_I_ROOT()) PRINT *,TRIM(Iam)//": StratChem time step length: ",cdt," seconds"

!  Create Chem Bundle
!  ------------------
   call ESMF_GridCompGet ( GC, grid=grid, rc=STATUS)
   VERIFY_(STATUS)

   call MAPL_GridGet ( grid, globalCellCountPerDim=DIMS, RC=STATUS)
   VERIFY_(STATUS)

   im = dims(1)
   jm = dims(2)
   nq = chemReg%nq

   call ESMF_GridGet(GRID, localDE=0, &
        staggerloc=ESMF_STAGGERLOC_CENTER, &
        computationalCount=DIMS, RC=STATUS)
   VERIFY_(STATUS)

!  Associate the Internal State fields with our legacy state 
!  ---------------------------------------------------------
   call MAPL_Get ( ggSTATE, INTERNALSPEC=InternalSpec, &
                                        INTERNAL_ESMF_STATE=internal, &
                                        LONS=LONS, &
                                        LATS=LATS, &
                                        RC=STATUS  )
   VERIFY_(STATUS)

! Local sizes of three dimensions
!--------------------------------
   i2 = dims(1)
   j2 = dims(2)
   km = dims(3)

! Broadcast necessary information the SC GridComp
! -----------------------------------------------
   gcChem%i1 = i1
   gcChem%i2 = i2
   gcChem%im = im
   gcChem%j1 = j1
   gcChem%j2 = j2
   gcChem%jm = jm
   gcChem%km = km

! Latitudes and longitudes [radians]
! ----------------------------------
   ALLOCATE(gcChem%lonRad(1:i2,1:j2),STAT=STATUS)
   VERIFY_(STATUS)
   gcChem%lonRad = LONS

   ALLOCATE(gcChem%latRad(1:i2,1:j2),STAT=STATUS)
   VERIFY_(STATUS)
   gcChem%latRad = LATS

! Cell area [m^{2}] from import state
! -----------------------------------
   CALL MAPL_GetPointer(impChem, AREA, 'AREA', RC=STATUS)
   VERIFY_(STATUS)
   ALLOCATE(gcChem%cellArea(1:i2,1:j2),STAT=STATUS)
   VERIFY_(STATUS)
   gcChem%cellArea = AREA

!  Initalize the legacy state but do not allocate memory for arrays
!  ----------------------------------------------------------------
   call Chem_BundleCreate_ ( chemReg, i1, i2, ig, im, j1, j2, jg, jm, km,  &
                             w_c, lon=lons(:,:), lat=lats(:,:), &
                             skipAlloc=.true., rc=STATUS )
   VERIFY_(STATUS)

   w_c%grid_esmf = grid  ! Will need this for I/O later

   allocate(w_c%delp(1:i2,1:j2,km),stat=STATUS)
   VERIFY_(STATUS)

!  Consistency Checks
!  ------------------
   _ASSERT( chemReg%i_XX == (chemReg%j_SC+1), 'needs informative message' )
   _ASSERT( size(InternalSpec) == (chemReg%n_SC+chemReg%n_XX), 'needs informative message' )

   do L = 1, size(InternalSpec)

      call MAPL_VarSpecGet ( InternalSpec(L),          &
                             SHORT_NAME = short_name,  &
                             RC=STATUS )
      VERIFY_(STATUS)

      N = chemReg%i_SC + L - 1 !  Assumption: SC and XX bins are contiguous

      call MAPL_GetPointer ( internal, NAME=short_name, ptr=w_c%qa(N)%data3d, &
                             rc = STATUS )
      VERIFY_(STATUS)

   end do

! STRATCHEM cannot be the AERO_PROVIDER
! -------------------------------------
     CALL ESMF_ConfigGetAttribute(CF, providerName, Default="GOCART.data", Label="AERO_PROVIDER:", RC=STATUS)
     VERIFY_(STATUS)
     IF(TRIM(providerName) == "STRATCHEM") THEN
      PRINT *,TRIM(Iam)//": STRATCHEM cannot be the AERO_PROVIDER"
      VERIFY_(ESMF_FAILURE)
     END IF

!  Call Legacy Initialize
!  ----------------------
   CALL SC_GridCompInitialize(gcChem, w_c, impChem, expChem, nymd, nhms, cdt, STATUS)
   VERIFY_(STATUS)

!  Stop timers
!  -----------
   CALL MAPL_TimerOff(ggState, "INITIALIZE")
   CALL MAPL_TimerOff(ggState, "TOTAL")

   RETURN_(ESMF_SUCCESS)

  END SUBROUTINE Initialize_


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Run_ --- Runs StratChem
!
! !INTERFACE:
!

   SUBROUTINE Run_ ( gc, impChem, expChem, clock, rc )

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
   type(SC_GridComp), pointer      :: gcChem      ! Grid Component
   type(Chem_Bundle), pointer      :: w_c         ! Chemical tracer fields     
   integer                         :: nymd, nhms  ! time
   real                            :: cdt         ! chemistry timestep (secs)
   real                            :: rdt         ! run (heatbeat) timestep (secs)
   real                            :: dtInverse   ! inverse chemistry timestep (1/s)
   integer                         :: i, iOX, irO3Ox, k, m, n, split

   type(ESMF_Config)               :: CF
   type(ESMF_Grid)                 :: grid
   type(ESMF_Time)                 :: TIME

   type(MAPL_MetaComp), pointer    :: ggState     ! GEOS Generic State

   real, pointer, dimension(:,:)   :: LATS
   real, pointer, dimension(:,:)   :: LONS
   type (MAPL_SunOrbit)            :: ORBIT
   real, allocatable, target       :: ZTH(:,:)
   real, allocatable               :: SLR(:,:)

   REAL, POINTER, DIMENSION(:,:,:) :: PLE
   REAL, POINTER, DIMENSION(:,:,:) :: DELP
   REAL, POINTER, DIMENSION(:,:)   :: TROPP
   REAL, POINTER, DIMENSION(:,:)   :: AGCMTROPP
   REAL, POINTER, DIMENSION(:,:)   :: SCTROPP
   REAL, POINTER, DIMENSION(:,:)   :: TO3
   REAL, POINTER, DIMENSION(:,:)   :: TTO3
   REAL, ALLOCATABLE               :: wrk(:,:)
   REAL, ALLOCATABLE               :: wgt(:,:)

   INTEGER :: i1, i2, j1, j2, km
   LOGICAL :: isLeapYear
   REAL :: dayOfYear, SCBaseP
   REAL(ESMF_KIND_R8) :: dayOfYear_r8

   TYPE(ESMF_Field)                :: FIELD

   REAL, POINTER, DIMENSION(:,:,:,:) :: sInitial
   REAL, POINTER, DIMENSION(:,:,:)   :: sIncrement
   LOGICAL                           :: doingTendencies
   INTEGER                           :: nAlloc
   LOGICAL, ALLOCATABLE              :: doMyTendency(:)
   CHARACTER(LEN=ESMF_MAXSTR)        :: fieldName, incFieldName

!  Get my name and set-up traceback handle
!  ---------------------------------------
   call ESMF_GridCompGet( GC, NAME=COMP_NAME, CONFIG=CF, GRID=grid, RC=STATUS )
   VERIFY_(STATUS)
   Iam = TRIM(COMP_NAME)//"::Run_"

!  Get my internal MAPL_Generic state
!  -----------------------------------
   call MAPL_GetObjectFromGC ( GC, ggState, RC=STATUS)
   VERIFY_(STATUS)

!  Start timers
!  ------------
   CALL MAPL_TimerOn(ggState, "TOTAL")
   CALL MAPL_TimerOn(ggState, "RUN"  )

! Get parameters from generic state.
!-----------------------------------
   call MAPL_Get(ggState, LONS=LONS, LATS=LATS, ORBIT=ORBIT, RC=STATUS)
   VERIFY_(STATUS)

   allocate(ZTH(SIZE(LATS,1), SIZE(LATS,2)), STAT=STATUS)
   VERIFY_(STATUS)
   allocate(SLR(SIZE(LATS,1), SIZE(LATS,2)), STAT=STATUS)
   VERIFY_(STATUS)

!  Update solar zenith angle
!  --------------------------
   call MAPL_SunGetInsolation(LONS, LATS,  &
        ORBIT, ZTH, SLR, CLOCK=CLOCK,      &
        RC=STATUS  )
   VERIFY_(STATUS)

!  Get pre-ESMF parameters from gc and clock
!  -----------------------------------------
   call extract_ ( gc, clock, chemReg, gcChem, w_c, nymd, nhms, rdt, split, STATUS )
   VERIFY_(STATUS)

!  Obtain day of year. ESMF returns 1.00-366.x, which is what we want.
!  -------------------------------------------------------------------
   CALL ESMF_ClockGet(CLOCK, currTIME=TIME, rc=STATUS)
   VERIFY_(STATUS)
   CALL ESMF_TimeGet(TIME, dayOfYear_r8=dayOfYear_r8, rc=STATUS)
   VERIFY_(STATUS)
   dayOfYear = dayOfYear_r8

!  Is this a leap year?
!  --------------------
   isLeapYear = ESMF_TimeIsLeapYear(TIME,rc=STATUS)
   VERIFY_(STATUS)

   gcChem%dayOfYear  = dayOfYear
   gcChem%isLeapYear = isLeapYear

!  Set pointers for cosine zenith angle
!  ------------------------------------
   w_c%cosz => zth

!  Layer interface pressures
!  -------------------------
   CALL MAPL_GetPointer(impChem, PLE, 'PLE', RC=STATUS)
   VERIFY_(STATUS)

!  Layer pressure thickness
!  ------------------------
   CALL MAPL_GetPointer(impChem, DELP, 'DELP', RC=STATUS)
   VERIFY_(STATUS)
   w_c%delp = DELP

   i1 = gcChem%i1
   i2 = gcChem%i2
   j1 = gcChem%j1
   j2 = gcChem%j2
   km = gcChem%km

!  The objective is to operate stratospheric chemistry "down" to the tropopause. However,
!  MAPL_UNDEFs occasionally appear in the imported tropopause pressures, TROPP.  Where
!  they are present, replace them with the most recent valid tropopause pressures, and
!  save the updated values in the surface layer of RO3OX, i.e. w_c%qa(irO3Ox)%data3d(:,:,km). 
!  ------------------------------------------------------------------------------------------
   irO3Ox = gcChem%irO3Ox
   CALL MAPL_GetPointer(impChem, TROPP, 'TROPP', RC=STATUS)
   VERIFY_(STATUS)
   WHERE(TROPP /= MAPL_UNDEF) w_c%qa(irO3Ox)%data3d(:,:,km) = TROPP
   
!  It appears that other anomalously high TROPP can be generated that SC cannot accommodate.
!  So impose an upper limit (Pa) to the pressure range through which SC is to be applied.
!  -----------------------------------------------------------------------------------------
   CALL ESMF_ConfigGetAttribute(CF, SCBaseP, LABEL="SC_BASE_PRESSURE:", DEFAULT=40000.00, RC=STATUS)
   VERIFY_(STATUS)
   WHERE(TROPP > SCBaseP) w_c%qa(irO3Ox)%data3d(:,:,km) = SCBaseP

   IF( ANY(w_c%qa(irO3Ox)%data3d(:,:,km) == MAPL_UNDEF) ) THEN
    PRINT *,TRIM(Iam)//": At least one invalid tropopause pressure."
    STATUS = 1
    VERIFY_(STATUS)
   END IF

!  For comparison purposes, export the SC base pressures and the
!  imported TROPP. Note that AGCM updates TROPP before HISTORY is written.
!  -----------------------------------------------------------------------
   CALL MAPL_GetPointer(expChem, SCTROPP, 'SCTROPP', RC=STATUS)
   VERIFY_(STATUS)
   IF(ASSOCIATED(SCTROPP)) SCTROPP = w_c%qa(irO3Ox)%data3d(:,:,km)
   CALL MAPL_GetPointer(expChem, AGCMTROPP, 'AGCMTROPP', RC=STATUS)
   VERIFY_(STATUS)
   IF(ASSOCIATED(AGCMTROPP)) AGCMTROPP = TROPP

! Are species tendencies requested?
! ---------------------------------
   m = ChemReg%i_SC
   n = ChemReg%j_SC
   ALLOCATE(doMyTendency(m:n), STAT=STATUS)
   VERIFY_(STATUS)
   doMyTendency(:) = .FALSE.

   nAlloc = 0
   DO i = m,n

    fieldName = TRIM(chemReg%vname(i))
    incFieldName = TRIM(fieldName)//"_SCTEND"

    CALL MAPL_GetPointer(expChem, sIncrement, TRIM(incFieldName), RC=STATUS)
    IF( STATUS /= ESMF_SUCCESS ) THEN
      print*,'SCTEND problem Getting Pointer for '//TRIM(incFieldName)
    ENDIF
    VERIFY_(STATUS)

    IF(ASSOCIATED(sIncrement)) THEN
     NULLIFY(sIncrement)
     nAlloc = nAlloc+1
     doMyTendency(i) = .TRUE.
    END IF

   END DO
    
   IF(nAlloc > 0) doingTendencies = .TRUE.

!  Save current species configurations so chemical tendencies can be calculated.
!  NOTE: Restricted to transported species.
!  ----------------------------------------------------------------------------
   StoreIC: IF(doingTendencies) THEN

    ALLOCATE(sInitial(1:i2,1:j2,km,nAlloc), STAT=STATUS)
    VERIFY_(STATUS)

    k = 1
    m = ChemReg%i_SC
    n = ChemReg%j_SC
  
    DO i = m,n
     IF(doMyTendency(i)) THEN
      sInitial(:,:,:,k) = w_c%qa(i)%data3d
      k = k+1
     END IF
    END DO
   
   END IF StoreIC

!  Chemistry time step length (seconds)
!  ------------------------------------
   cdt = rdt/split

!  Run the chemistry
!  -----------------
   DO i = 1,split
    CALL SC_GridCompRun(gcChem, w_c, impChem, expChem, nymd, nhms, cdt, STATUS)
    VERIFY_(STATUS)
   END DO

!  Update age-of-air (days)
!  ------------------------
   n = gcChem%iAOA
   w_c%qa(n)%data3d(:,:,:) = w_c%qa(n)%data3d(:,:,:)+rdt/86400.00
   w_c%qa(n)%data3d(:,:,km) = 0.00

   DEALLOCATE(SLR)
   DEALLOCATE(ZTH)

!  Total ozone: In each layer
!   molecules m^{-2} =  O3(vmr) * Avogadro's number * dp / ( mwt air * g )
!  The most recent valid tropopause pressures are stored in RO3OX(:,:,km)
!  -----------------------------------------------------------------------
   CALL MAPL_GetPointer(impChem, PLE, 'PLE', RC=STATUS)
   VERIFY_(STATUS)
   CALL MAPL_GetPointer(expChem,  TO3,  'SCTO3', RC=STATUS)
   VERIFY_(STATUS)
   IF(ASSOCIATED( TO3))  TO3 = 0.00
   CALL MAPL_GetPointer(expChem, TTO3, 'SCTTO3', RC=STATUS)
   VERIFY_(STATUS)
   IF(ASSOCIATED(TTO3)) TTO3 = 0.00

   DoingTotalOzone: IF(ASSOCIATED(TTO3) .OR. ASSOCIATED(TO3)) THEN

    iOX = gcChem%iOx
    irO3Ox = gcChem%irO3Ox
    
    ALLOCATE(wrk(SIZE(LATS,1), SIZE(LATS,2)), STAT=STATUS)
    VERIFY_(STATUS)
    ALLOCATE(wgt(SIZE(LATS,1), SIZE(LATS,2)), STAT=STATUS)
    VERIFY_(STATUS)

    DO k = 1,km
     wrk = w_c%qa(iOX)%data3d(:,:,k)*(PLE(:,:,k)-PLE(:,:,k-1))*(MAPL_AVOGAD/2.69E+20)/(MAPL_AIRMW*MAPL_GRAV)
     IF(ASSOCIATED( TO3)) TO3 = TO3+wrk

     IF(ASSOCIATED(TTO3)) THEN
      wgt  = MAX(0.0,MIN(1.0,(PLE(:,:,k)-w_c%qa(irO3Ox)%data3d(:,:,km))/(PLE(:,:,k)-PLE(:,:,k-1))))
      TTO3 = TTO3+wrk*wgt  
     END IF

    END DO

    DEALLOCATE(wrk)
    DEALLOCATE(wgt)

   END IF DoingTotalOzone

!  Obtain chemical tendencies and fill export states.  NOTE: Restricted to transported species.
!  --------------------------------------------------------------------------------------------
   StoreTendencies: IF(doingTendencies) THEN

    k = 1
    m = ChemReg%i_SC
    n = ChemReg%j_SC
    dtInverse = 1.00/rdt

    DO i = m,n

     IF(doMyTendency(i)) THEN

      fieldName = TRIM(chemReg%vname(i))
      incFieldName = TRIM(fieldName)//"_SCTEND"

      CALL MAPL_GetPointer(expChem, sIncrement, TRIM(incFieldName), RC=STATUS)
      VERIFY_(STATUS)
      IF(ASSOCIATED(sIncrement)) sIncrement = (w_c%qa(i)%data3d - sInitial(:,:,:,k))*dtInverse

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

!  Stop timers
!  -----------
   CALL MAPL_TimerOff(ggState, "RUN"  )
   CALL MAPL_TimerOff(ggState, "TOTAL")

   RETURN_(ESMF_SUCCESS)

   END SUBROUTINE Run_

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Finalize_ --- Finalize SC_GridComp (ESMF)
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
   type(SC_GridComp), pointer      :: gcChem      ! Grid Component
   type(Chem_Bundle), pointer      :: w_c         ! Chemical tracer fields     
   type(MAPL_MetaComp), pointer    :: ggState     ! GEOS Generic State
   integer                         :: nymd, nhms  ! time
   integer                         :: split
   real                            :: rdt         ! chemistry timestep (secs)

    type(StratChem_state), pointer  :: state

!  Get my name and set-up traceback handle
!  ---------------------------------------
   call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
   VERIFY_(STATUS)
   Iam = TRIM(COMP_NAME)//"::Finalize_"

!  Get my internal MAPL_Generic state
!  -----------------------------------
   CALL MAPL_GetObjectFromGC(GC, ggState, RC=STATUS)

!  Start timers
!  ------------
   CALL MAPL_TimerOn(ggState, "TOTAL")
   CALL MAPL_TimerOn(ggState, "FINALIZE")

!  Get parameters from gc and clock
!  --------------------------------
   call extract_ ( gc, clock, chemReg, gcChem, w_c, nymd, nhms, rdt, split, STATUS, &
                   state = state )
   VERIFY_(STATUS)

!  Finalize
!  --------
   call SC_GridCompFinalize ( gcChem, w_c, impChem, expChem, nymd, nhms, rdt, STATUS )
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
   deallocate ( state%chemReg, state%gcChem, state%w_c, stat = STATUS )
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

    SUBROUTINE extract_ ( gc, clock, chemReg, gcChem, w_c, nymd, nhms, rdt, split, rc, state )

    type(ESMF_GridComp), intent(INOUT) :: gc
    type(ESMF_Clock), intent(in)     :: clock
    type(Chem_Registry), pointer     :: chemReg
    type(SC_GridComp), pointer       :: gcChem
    type(Chem_Bundle), pointer       :: w_c
    integer, intent(out)             :: nymd, nhms
    real, intent(out)                :: rdt
    integer, intent(out)             :: split, rc
    type(StratChem_state), pointer, optional   :: state


    type(StratChem_state), pointer    :: myState

!   ErrLog Variables
!   ----------------
    character(len=ESMF_MAXSTR)      :: IAm
    integer                         :: STATUS
    character(len=ESMF_MAXSTR)      :: COMP_NAME

    type(ESMF_Time)      :: TIME
    type(ESMF_Config)    :: CF
    type(StratChem_Wrap)  :: wrap
    integer              :: IYR, IMM, IDD, IHR, IMN, ISC


!   Get my name and set-up traceback handle
!   ---------------------------------------
    call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // 'extract_'

    rc = 0

!   Get my internal state
!   ---------------------
    call ESMF_UserCompGetInternalState(gc, 'StratChem_state', WRAP, STATUS)
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
    if ( .not. associated(myState%gcChem) ) then
         allocate ( myState%gcChem, stat=STATUS )
         VERIFY_(STATUS)
    end if
    if ( .not. associated(myState%w_c) ) then
         allocate ( myState%w_c, stat=STATUS )
         VERIFY_(STATUS)
    end if

    chemReg => myState%chemReg
    gcChem  => myState%gcChem
    w_c     => myState%w_c

!   Get the configuration
!   ---------------------
    call ESMF_GridCompGet ( GC, CONFIG = CF, RC=STATUS )
    VERIFY_(STATUS)

!   Get time step
!   -------------
    call ESMF_ConfigGetAttribute ( CF, rdt, LABEL="RUN_DT:", RC=STATUS )
    VERIFY_(STATUS)
    call ESMF_ConfigGetAttribute ( CF, split, LABEL="SC_SPLIT:", DEFAULT=1, RC=STATUS )
    VERIFY_(STATUS)

!   Get time stamp
!   --------------
    call ESMF_ClockGet(CLOCK,currTIME=TIME,rc=STATUS)
    VERIFY_(STATUS)

    call ESMF_TimeGet(TIME ,YY=IYR, MM=IMM, DD=IDD, H=IHR, M=IMN, S=ISC, rc=STATUS)
    VERIFY_(STATUS)

    call MAPL_PackTime(NYMD,IYR,IMM,IDD)
    call MAPL_PackTime(NHMS,IHR,IMN,ISC)

    RETURN_(ESMF_SUCCESS)

   END SUBROUTINE extract_

 END MODULE StratChem_GridCompMod
