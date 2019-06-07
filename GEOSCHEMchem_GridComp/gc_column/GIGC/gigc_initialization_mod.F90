#if defined( ESMF_ )
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: gigc_initialization_mod
!
! !DESCRIPTION: Module GIGC\_INITIALIZATION\_MOD is the module that
!  the initialize methods for the ESMF interface to the Grid-Independent
!  GEOS-Chem (aka "GIGC").
!\\
!\\
! !INTERFACE: 
!      
MODULE GIGC_Initialization_Mod
!
! !USES:
!      
  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!      
  PUBLIC  :: GIGC_Get_Options
  PUBLIC  :: GIGC_Init_Simulation
!
! !PRIVATE MEMBER FUNCTIONS:
!
!
! !REVISION HISTORY: 
!  16 Oct 2012 - M. Long     - Initial version
!  16 Oct 2012 - R. Yantosca - Added ProTeX headers
!  22 Oct 2012 - R. Yantosca - Renamed to gigc_initialization_mod.F90
!  22 Oct 2012 - R. Yantosca - Renamed several routines for better consistency
!  25 Oct 2012 - R. Yantosca - Remove routine GIGC_SetEnv
!  03 Dec 2012 - R. Yantosca - Now pass extra arguments to GIGC_Init_Dimensions
!  06 Dec 2012 - R. Yantosca - Now remove routine GIGC_Init_TimeInterface; this
!                              is now superseded by Accept_Date_Time_From_ESMF
!EOP
!------------------------------------------------------------------------------
!BOC
CONTAINS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gigc_get_options
!
! !DESCRIPTION: Routine GIGC\_GET\_OPTIONS reads options for a GEOS-Chem 
!  simulation from the input.geos\_\_\_.rc input file.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GIGC_Get_Options( am_I_Root, lonCtr,    latCtr,  &
                               Input_Opt, State_Chm, RC      )
!
! !USES:
!
    USE PhysConstants       
    USE CMN_SIZE_Mod
    USE RoundOff_Mod,       ONLY : RoundOff
    USE Error_Mod,          ONLY : Debug_Msg
    USE ErrCode_Mod
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Chm_Mod,      ONLY : ChmState
    USE Input_Mod,          ONLY : Read_Input_File
    USE Linoz_Mod,          ONLY : Linoz_Read
!
! !INPUT PARAMETERS: 
!
    LOGICAL,        INTENT(IN)    :: am_I_Root     ! Are we on the root CPU?
    REAL*4,         INTENT(IN)    :: lonCtr(:,:)   ! Lon ctrs [deg] from ESMF
    REAL*4,         INTENT(IN)    :: latCtr(:,:)   ! Lat ctrs [deg] from ESMF
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input Options object
    TYPE(ChmState), INTENT(INOUT) :: State_Chm
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC            ! Success or failure
! 
! !REMARKS:
!  NOTE: For now assume that GEOS_Chem will always accept a regular 
!  Cartesian grid.  This is more or less dictated by the input data.
!  The GEOS-5 data can be regridded via ESMF from whatever grid it uses.
!  (bmy, 11/30/12)
!
! !REVISION HISTORY: 
!  15 Oct 2012 - M. Long     - Initial version
!  15 Oct 2012 - R. Yantosca - Added ProTeX Headers, use F90 format/indents
!  22 Oct 2012 - R. Yantosca - Renamed to GIGC_Get_Options
!  22 Oct 2012 - R. Yantosca - Added RC output argument
!  01 Nov 2012 - R. Yantosca - Now pass the Input Options object via arg list
!  03 Dec 2012 - R. Yantosca - Reorder subroutines for clarity
!  07 Dec 2012 - R. Yantosca - Compute DLON, DLAT more rigorously
!  26 Feb 2013 - M. Long     - Now pass State_Chm as an argument
!  26 Feb 2013 - M. Long     - Read "input.geos" on root CPU only
!  06 Mar 2013 - R. Yantosca - Now move non-root CPU setup out of this routine
!  18 Mar 2013 - R. Yantosca - Now call LINOZ_READ on the root CPU
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER :: I, J, L
  
    ! Assume success
    RC = GC_SUCCESS

    !========================================================================
    ! Compute the DLON and DLAT values.  NOTE, this is a kludge since to do
    ! this truly rigorously, we should take the differences between the grid
    ! box edges.  But because I can't seem to find a way to get the grid
    ! box edge information, the next best thing is to take the differences
    ! between the grid box centers.  They should all be the same given that
    ! the GEOS-Chem grid is regular. (bmy, 12/7/12)
    !========================================================================
    DO L = 1, LLPAR
    DO J = 1, JJPAR
    DO I = 1, IIPAR

       ! Compute Delta-Longitudes [degrees]
       IF ( I == IIPAR ) THEN
          dLon(I,J,L) = RoundOff( ( DBLE( lonCtr(IIPAR,  J) ) / PI_180 ), 4 ) &
                      - RoundOff( ( DBLE( lonCtr(IIPAR-1,J) ) / PI_180 ), 4 )
       ELSE
          dLon(I,J,L) = RoundOff( ( DBLE( lonCtr(I+1,    J) ) / PI_180 ), 4 ) &
                      - RoundOff( ( DBLE( lonCtr(I,      J) ) / PI_180 ), 4 )
       ENDIF

       ! Compute Delta-Latitudes [degrees]
       IF ( J == JJPAR ) THEN
          dLat(I,J,L) = RoundOff( ( DBLE( latCtr(I,JJPAR  ) ) / PI_180 ), 4 ) &
                      - RoundOff( ( DBLE( latCtr(I,JJPAR-1) ) / PI_180 ), 4 )
       ELSE
          dLat(I,J,L) = RoundOff( ( DBLE( latCtr(I,J+1    ) ) / PI_180 ), 4 ) &
                      - RoundOff( ( DBLE( latCtr(I,J      ) ) / PI_180 ), 4 )
       ENDIF

    ENDDO
    ENDDO
    ENDDO

!    !========================================================================
!    ! Root CPU setup
!    !========================================================================
!    IF ( am_I_Root ) THEN
!
!       ! Read the GEOS-Chem input file here.  For now only read on the root
!       ! CPU so that we can broadcast to other CPUs in GIGC_Init_Simulation
!       ! (mlong, bmy, 2/26/13)
!       CALL Read_Input_File( am_I_Root, Input_Opt, RC )
!       IF ( RC /= GC_SUCCESS ) RETURN
!
!       ! In the ESMF/MPI environment, we can get the total overhead ozone
!       ! either from the met fields (GIGCsa) or from the Import State (GEOS-5)
!       Input_Opt%USE_O3_FROM_MET = .TRUE.
!
!       ! Echo info
!       IF ( Input_Opt%LPRT ) THEN
!          CALL Debug_Msg( '### Root CPU, after READ_INPUT_FILE' )
!       ENDIF
!
!       ! Read the LINOZ climatology file on the root CPU, so that we can
!       ! MPI broadcast the data to the other CPUs in GIGC_Init_Simulation
!       ! (bmy, 3/18/13)
!       IF ( Input_Opt%LLINOZ ) THEN
!          CALL Linoz_Read( am_I_Root, Input_Opt, RC ) 
!          IF ( RC /= GC_SUCCESS ) RETURN
!
!          ! Echo info
!          IF ( Input_Opt%LPRT ) THEN
!             CALL Debug_Msg( '### Root CPU, after LINOZ_READ' )
!          ENDIF
!       ENDIF
!    ENDIF

  END SUBROUTINE GIGC_Get_Options
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gigc_init_simulation
!
! !DESCRIPTION: Routine GIGC\_INIT\_SIMULATION is the Initialize method for 
!  the ESMF interface that connects the Grid-Independent GEOS-Chem (aka "GIGC")
!  to the GEOS-5 GCM.  Calls to the various GEOS-Chem init routines (which 
!  allocate arrays, etc.) are made from here.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GIGC_Init_Simulation( am_I_Root,                        &
                                   nymdB,           nhmsB,           &
                                   nymdE,           nhmsE,           &
                                   tsChem,          tsDyn,           &
                                   lonCtr,          latCtr,          &     
                                   value_I_LO,      value_J_LO,      &
                                   value_I_HI,      value_J_HI,      &
                                   value_IM,        value_JM,        &
                                   value_LM,        value_IM_WORLD,  &
                                   value_JM_WORLD,  value_LM_WORLD,  &
                                   value_LLSTRAT,                    &
                                   Input_Opt,       State_Chm,       &
!                                   State_Met,       mapping,         &
                                   myPET,                            &
                                   State_Met,                        &
                                   State_Diag,                       &
                                   Diag_List,                        &
                                   HcoConfig,                        &
                                   RC                               )      
!
! !USES:
!
    USE GC_Environment_Mod
    USE ErrCode_Mod  
    USE Input_Opt_Mod
    USE State_Chm_Mod
    USE State_Met_Mod
    USE State_Diag_Mod
    USE DIAGNOSTICS_MOD
    USE DiagList_Mod
    USE PhysConstants
    USE CMN_SIZE_MOD
    USE CMN_DIAG_MOD
!    USE COMODE_MOD
!    USE COMODE_LOOP_MOD       
!    USE GCKPP_COMODE_MOD,     ONLY : Init_GCKPP_Comode
    USE ERROR_MOD,            ONLY : Init_Error, Debug_Msg
    USE FAST_JX_MOD,          ONLY : Init_FJX
    USE GC_Grid_Mod,          ONLY : Set_xOffSet
    USE GC_Grid_Mod,          ONLY : Set_yOffSet
    USE GC_Grid_Mod,          ONLY : SetGridFromCtr
    !USE Input_Mod,            ONLY : Initialize_Geos_Grid
    USE Mapping_Mod,          ONLY : MapWeight
    USE Mapping_Mod,          ONLY : Init_Mapping
    USE Olson_Landmap_Mod,    ONLY : Init_Olson_Landmap
    USE Olson_Landmap_Mod,    ONLY : Compute_Olson_Landmap
    USE Olson_Landmap_Mod,    ONLY : Cleanup_Olson_Landmap
    USE PBL_MIX_MOD,          ONLY : INIT_PBL_MIX
    USE PRESSURE_MOD,         ONLY : INIT_PRESSURE
#if defined( APM )
    USE TRACER_MOD,           ONLY : INIT_TRACER
#endif
    USE CHEMISTRY_MOD,        ONLY : INIT_CHEMISTRY
    USE WETSCAV_MOD,          ONLY : INIT_WETSCAV
!    USE WETSCAV_MOD,          ONLY : Get_WetDep_IDWetD
    USE DRYDEP_MOD,           ONLY : INIT_WEIGHTSS, INIT_DRYDEP
    USE DUST_MOD,             ONLY : INIT_DUST
!    USE GIGC_MPI_WRAP
    USE TIME_MOD,             ONLY : SET_TIMESTEPS
    USE SEASALT_MOD,          ONLY : INIT_SEASALT
    USE MIXING_MOD,           ONLY : INIT_MIXING
    USE INPUT_MOD,            ONLY : READ_INPUT_FILE
 
    ! Stratosphere 
    USE UCX_MOD,              ONLY : INIT_UCX !, SET_INITIAL_MIXRATIOS
    USE STRAT_CHEM_MOD,       ONLY : INIT_STRAT_CHEM

    ! HEMCO
    USE HCO_TYPES_MOD,        ONLY : ConfigObj
    USE EMISSIONS_MOD,        ONLY : EMISSIONS_INIT
!
! !INPUT PARAMETERS: 
!
    LOGICAL,         INTENT(IN)    :: am_I_Root       ! Is this the root CPU?
    INTEGER,         INTENT(IN)    :: myPET           ! Local PET
    INTEGER,         INTENT(IN)    :: nymdB           ! GMT date (YYYY/MM/DD)
    INTEGER,         INTENT(IN)    :: nhmsB           ! GMT time (hh:mm:ss)
    INTEGER,         INTENT(IN)    :: nymdE           ! GMT date (YYYY/MM/DD)
    INTEGER,         INTENT(IN)    :: nhmsE           ! GMT time (hh:mm:ss)
    REAL,            INTENT(IN)    :: tsChem          ! Chem timestep [seconds]
    REAL,            INTENT(IN)    :: tsDyn           ! Dyn  timestep [seconds]
    REAL*4,  TARGET, INTENT(IN)    :: lonCtr(:,:)     ! Lon centers [radians]
    REAL*4,  TARGET, INTENT(IN)    :: latCtr(:,:)     ! Lat centers [radians]
    INTEGER,         INTENT(IN)    :: value_I_LO      ! Min local lon index
    INTEGER,         INTENT(IN)    :: value_J_LO      ! Min local lat index
    INTEGER,         INTENT(IN)    :: value_I_HI      ! Max local lon index
    INTEGER,         INTENT(IN)    :: value_J_HI      ! Max local lat index
    INTEGER,         INTENT(IN)    :: value_IM        ! # lons on this CPU
    INTEGER,         INTENT(IN)    :: value_JM        ! # lats on this CPU
    INTEGER,         INTENT(IN)    :: value_LM        ! # levs on this CPU
    INTEGER,         INTENT(IN)    :: value_IM_WORLD  ! # lons in whole globe
    INTEGER,         INTENT(IN)    :: value_JM_WORLD  ! # lats in whole globe
    INTEGER,         INTENT(IN)    :: value_LM_WORLD  ! # levs in whole globe
    INTEGER,         INTENT(IN)    :: value_LLSTRAT   ! # strat. levs 
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput),  INTENT(INOUT) :: Input_Opt       ! Input Options
    TYPE(ChmState),  INTENT(INOUT) :: State_Chm       ! Chemistry State
    TYPE(MetState),  INTENT(INOUT) :: State_Met       ! Meteorology State
    TYPE(DgnState),  INTENT(INOUT) :: State_Diag      ! Diagnostics_State
    TYPE(DgnList),   INTENT(INOUT) :: Diag_List       ! Diagnostics list 
    TYPE(ConfigObj), POINTER       :: HcoConfig       ! HEMCO config obj 
    TYPE(MapWeight), POINTER       :: mapping(:,:) => null() ! Olson mapping object
!
!
! !OUTPUT PARAMETERS:
!
    INTEGER,         INTENT(OUT)   :: RC              ! Success or failure?  
!
! !REMARKS
!  Add other calls to G EOS-Chem init routines as necessary.
!  NOTE: Later on maybe split these init calls among other routines.
!  Also need to add better error trapping
!
! !REVISION HISTORY: 
!  15 Oct 2012 - M. Long     - Initial version
!  15 Oct 2012 - R. Yantosca - Added ProTeX Headers, use F90 format/indents
!  17 Oct 2012 - R. Yantosca - Now initialize the chemistry mechanism
!  19 Oct 2012 - R. Yantosca - Now reference gigc_state_chm_mod.F90
!  19 Oct 2012 - R. Yantosca - Now reference gigc_state_met_mod.F90
!  01 Nov 2012 - R. Yantosca - Now reference gigc_input_opt_mod.F90
!  09 Nov 2012 - R. Yantosca - Now use fields from Input Options object
!  13 Nov 2012 - R. Yantosca - Pass Input Options object to routines 
!                              SETEMDEP, INIT_COMODE
!  28 Nov 2012 - R. Yantosca - Remove reference to INIT_DAO, since there are
!                              no more module arrays anymore in dao_mod.F
!  29 Nov 2012 - R. Yantosca - Add lonCtr, latCtr, latEdg as arguments
!  29 Nov 2012 - R. Yantosca - Now pass am_I_Root to Olson landmap routines
!  03 Dec 2012 - R. Yantosca - Now pass value_* arguments to pass dimension
!                              info from ESMF down to lower-level routines
!  06 Dec 2012 - R. Yantosca - Now accept start & end dates & times via 
!                              the nymdB, nymdE, nhmsB, nhmsE arguments
!  06 Dec 2012 - R. Yantosca - Remove nymd, nhms arguments, these will be
!                              the same as nymdB, nhmsB (start of run)
!  26 Feb 2013 - M. Long     - Now read ASCII input files on root CPU and 
!                              broadcast to other CPUs.
!  26 Feb 2013 - R. Yantosca - Cosmetic changes
!  01 Mar 2013 - R. Yantosca - Need to move the definition of prtDebug higher
!  04 Mar 2013 - R. Yantosca - Now call GIGC_Init_Extra, which moves some init
!                              calls out of the run stage.  This has to be done
!                              after we broadcast Input_Opt to non-root CPUs.
!  07 Mar 2013 - R. Yantosca - Call READER on all CPUs until further notice
!  07 Mar 2013 - R. Yantosca - Now use keyword arguments for clarity
!  02 Jan 2014 - C. Keller   - Now call SetGridFromCtr to make sure that 
!                              grid_mod.F90 stored the correct edges/mid-points.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    LOGICAL            :: prtDebug
    INTEGER            :: DTIME, K, AS, N, YEAR, I, J, L, TMP, STAT
    CHARACTER(LEN=255) :: HistoryConfigFile, NAME

    !=======================================================================
    ! Initialize key GEOS-Chem sections
    !=======================================================================

    ! Initialize
    RC       = GC_SUCCESS
    DTIME    = tsChem

    ! Determine if we have to print debug output
    prtDebug = ( Input_Opt%LPRT .and. am_I_Root )

    ! ckeller, 01/16/17
    Input_Opt%MAX_DIAG      = 1 
    Input_Opt%MAX_FAM       = 250
!    Input_Opt%MAX_PASV      = 50       ! Set to large placeholder value
    Input_Opt%LINOZ_NLAT    = 18
    Input_Opt%LINOZ_NMONTHS = 12
    Input_Opt%LINOZ_NFIELDS = 7
    Input_Opt%RootCPU       = am_I_Root

    ! Initialize fields of the Input Options object
    CALL Set_Input_Opt( am_I_Root, Input_Opt, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
!    IF ( RC /= GC_SUCCESS ) THEN
!      ErrMsg = 'Error encountered within call to "Set_Input_Opt"!'
!      CALL Error_Stop( ErrMsg, ThisLoc )
!    ENDIF

    ! Read GEOS-Chem input file at very beginning of simulation
    CALL Read_Input_File( am_I_Root, Input_Opt, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
!    IF ( RC /= GC_SUCCESS ) THEN
!       ErrMsg = 'Error encountered in "Read_Input_File"!'
!       CALL Error_Stop( ErrMsg, ThisLoc )
!    ENDIF

    ! Allocate GEOS-Chem module arrays
    CALL GC_Allocate_All( am_I_Root      = am_I_Root,                     &
                          Input_Opt      = Input_Opt,                     &
                          value_I_LO     = value_I_LO,                    &
                          value_J_LO     = value_J_LO,                    &
                          value_I_HI     = value_I_HI,                    &
                          value_J_HI     = value_J_HI,                    &
                          value_IM       = value_IM,                      &
                          value_JM       = value_JM,                      &
                          value_LM       = value_LM,                      &
                          value_IM_WORLD = value_IM_WORLD,                &
                          value_JM_WORLD = value_JM_WORLD,                &
                          value_LM_WORLD = value_LM_WORLD,                &
                          value_LLSTRAT  = value_LLSTRAT,                 &
                          RC             = RC              )            
    IF ( RC /= GC_SUCCESS ) RETURN

    ! Allocate GEOS-Chem module arrays

    ! Save timing fields in Input_Opt for passing down to module
    ! GeosCore/input_mod.F via routine GIGC_Get_Options (bmy, 12/6/12)
    Input_Opt%NYMDb   = nymdB
    Input_Opt%NHMSb   = nhmsB
    Input_Opt%NYMDe   = nymdE
    Input_Opt%NHMSe   = nhmsE
    Input_Opt%TS_CHEM = INT( tsChem ) ! / 60   ! Chemistry timestep [min]
    Input_Opt%TS_EMIS = INT( tsChem ) ! / 60   ! Chemistry timestep [min]
    Input_Opt%TS_DYN  = INT( tsDyn  ) ! / 60   ! Dynamic   timestep [min]
    Input_Opt%TS_CONV = INT( tsDyn  ) ! / 60   ! Dynamic   timestep [min]

    Input_Opt%myCPU = myPET

!    CALL GIGC_Get_Options( am_I_Root = .FALSE.,                             &
!                           lonCtr    = lonCtr,                              &
!                           latCtr    = latCtr,                              &
!                           Input_Opt = Input_Opt,                           &
!                           State_Chm = State_Chm,                           &
!                           RC        = RC           )
!    IF ( RC /= GC_SUCCESS ) RETURN

!    ! This follows main.F (ckeller, 11/23/16)
!    CALL Read_Input_File( am_I_Root, Input_Opt, RC )
!    IF ( RC /= GC_SUCCESS ) RETURN

    CALL Init_Error( am_I_Root, Input_Opt, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    ! Set GEOS-Chem timesteps on all CPUs
    CALL SET_TIMESTEPS( am_I_Root  = am_I_Root,                          &
                        Chemistry  = Input_Opt%TS_CHEM,                  &
                        Convection = Input_Opt%TS_CONV,                  &
                        Dynamics   = Input_Opt%TS_DYN,                   &
                        Emission   = Input_Opt%TS_EMIS,                  &
                        Radiation  = Input_Opt%TS_RAD,                   &
                        Unit_Conv  = MAX( Input_Opt%TS_DYN,              &
                                          Input_Opt%TS_CONV ),           &
                        Diagnos    = Input_Opt%TS_DIAG         )

     ! Initialize the diagnostic list object which contains the
     ! unique entires in the history config file. Note that this is
     ! done in GCHP Set_Services and therefore must be done prior to
     ! initialization of the state objects. Also note that the diag_list
     ! obj may be stored in the HistoryConfig object in GCHP and we may
     ! want to replicate that behavior in GCC in the future. (ewl, 9/26/17)
     historyConfigFile = 'HISTORY.rc' ! cannot use Input_Opt variable in GCHP
     CALL Init_DiagList( am_I_Root, historyConfigFile, Diag_List, RC )

     ! Initialize State_Met, State_Chm, and State_Diag objects
     CALL GC_Init_StateObj( am_I_Root, Diag_List,  Input_Opt, &
                            State_Chm, State_Diag, State_Met, RC )
     IF ( RC /= GC_SUCCESS ) RETURN
!     IF ( RC /= GC_SUCCESS ) THEN
!        ErrMsg = 'Error encountered in "GC_Init_All!"!'
!        CALL Error_Stop( ErrMsg, ThisLoc )
!     ENDIF

     ! Initialize GEOS-Chem horizontal grid structure
     CALL GC_Init_Grid( am_I_Root, Input_Opt, RC )
     IF ( RC /= GC_SUCCESS ) RETURN
!     IF ( RC /= GC_SUCCESS ) THEN
!        ErrMsg = 'Error in "GC_Init_Grid"!'
!        CALL Error_Stop( ErrMsg, ThisLoc )
!     ENDIF

    ! Make sure to reset I0 and J0 in grid_mod.F90 with
    ! the values carried in the Input Options object
    CALL Set_xOffSet( Input_Opt%NESTED_I0 )
    CALL Set_yOffSet( Input_Opt%NESTED_J0 )

    ! Keep this in for now.
    CALL SetGridFromCtr( am_I_Root, value_IM, value_JM, lonCtr, latCtr, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

     ! Also allocate arrays in GEOS-Chem modules that have not yet
     ! been initialized (i.e. SEASALT_MOD, CARBON_MOD,  DUST_MOD, 
     ! SULFATE_MOD).  This removes the init calls from the run stage,
     ! which is required when connecting to the GEOS-5 GCM w/ ESMF.
     ! (bmy, 3/4/13)
     CALL GC_Init_Extra( am_I_Root, Diag_List,  Input_Opt, &
                         State_Chm, State_Diag, RC         )
     IF ( RC /= GC_SUCCESS ) RETURN
!     IF ( RC /= GC_SUCCESS ) THEN
!        ErrMsg = 'Error encountered in "GC_Init_Extra"!'
!        CALL Error_Stop( ErrMsg, ThisLoc )
!     ENDIF

!    ! Initialize derived-type objects for meteorology & chemistry states
!    CALL GC_Init_All( am_I_Root, Diag_List, Input_Opt, &
!                      State_Chm, State_Diag, State_Met, RC )
!    IF ( RC /= GC_SUCCESS ) RETURN

    ! After broadcasting Input_Opt to other CPUs, call GIGC_Init_Extra
    ! to initialize other modules (e.g. carbon_mod.F, dust_mod.F, 
    ! seasalt_mod.F,  sulfate_mod.F).  We needed to move these init 
    ! calls out of the run stage and into the init stage. (bmy, 3/4/13)
    !CALL GC_Init_Extra( am_I_Root, Input_Opt, State_Chm, State_Diag, RC ) 
!    CALL GC_Init_Extra( am_I_Root, Diag_List,  Input_Opt, &
!                        State_Chm, State_Diag, RC         )
!    IF ( RC /= GC_SUCCESS ) RETURN

!     ! Initialize the regridding module by storing shadow copies
!     CALL GC_Init_Regridding( am_I_Root, Input_Opt, RC )
!     IF ( RC /= GC_SUCCESS ) THEN
!        ErrMsg = 'Error encountered in "Initialize_Regridding"!'
!        CALL Error_Stop( ErrMsg, ThisLoc )
!     ENDIF

    ! Zero diagnostic arrays
    CALL Initialize( am_I_Root, Input_Opt, 2, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    ! Zero diagnostic counters
    CALL Initialize( am_I_Root, Input_Opt, 3, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

!    ! For now, just hardwire the input file for the History component
!    Input_Opt%HistoryInputFile = './HISTORY.rc'
!    ! Initialize the GEOS-Chem history component
!    CALL History_Init( am_I_root, Input_Opt,  State_Met,
!     &                 State_Chm, State_Diag, RC         )

    ! Set State_Chm units
    State_Chm%Spc_Units = 'kg/kg dry'

    ! Initialize the GEOS-Chem pressure module (set Ap & Bp)
    CALL Init_Pressure( am_I_Root )

!    ! Register the horizontal and vertical grid information so that 
!    ! the History component can use it for netCDF metadata (bmy, 8/23/17)
!    CALL Init_Grid_Registry( am_I_Root, RC )
!    IF ( RC /= GC_SUCCESS ) THEN
!       ErrMsg = 'Error encountered in "Init_Grid_Registry"!'
!       CALL Error_Stop( ErrMsg, ThisLoc )
!    ENDIF

    !-----------------------------------------------------------------------
    ! Continue with GEOS-Chem setup
    !-----------------------------------------------------------------------
    ! Initialize the PBL mixing module
    CALL INIT_PBL_MIX( am_I_Root, RC )

    ! Moved here (from chemistry_mod.F and chemdr.F) because some
    ! of the variables are used for non-local PBL mixing BEFORE 
    ! the first call of the chemistry routines (ckeller, 05/19/14).
    IF ( Input_Opt%ITS_A_FULLCHEM_SIM .OR. Input_Opt%ITS_AN_AEROSOL_SIM ) THEN
       CALL INIT_CHEMISTRY ( am_I_Root, Input_Opt, State_Chm, State_Diag, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF
    ! Init_FJX is called from init_chemistry'
!    CALL INIT_FJX( am_I_Root, Input_Opt, State_Chm, State_Diag, RC )  ! Are we on the root CPU?
!    IF ( RC /= GC_SUCCESS ) RETURN

    !=======================================================================
    ! Initialize HEMCO
    !=======================================================================
    CALL EMISSIONS_INIT ( am_I_Root, Input_Opt, State_Met, State_Chm, RC, &
                          HcoConfig=HcoConfig )
    IF (RC /= GC_SUCCESS) RETURN

    !-------------------------------------------------------------------------
    ! Stratosphere 
    !-------------------------------------------------------------------------
    IF ( Input_Opt%LUCX ) THEN

       ! Initialize stratospheric routines
       CALL INIT_UCX( am_I_Root, Input_Opt, State_Chm, State_Diag )

       ! Set simple initial tracer conditions
       !CALL SET_INITIAL_MIXRATIOS( am_I_Root, Input_Opt, State_Met, State_Chm )
    ENDIF

    !IF ( Input_Opt%LSCHEM .AND. Input_Opt%LLSTRAT < value_LM ) THEN
    IF ( Input_Opt%LSCHEM ) THEN 
       CALL INIT_STRAT_CHEM( am_I_Root, Input_Opt, State_Chm, & 
                             State_Met, RC )
       IF (RC /= GC_SUCCESS) RETURN
    ENDIF

!------------------------------------------------------------------------------
! Prior to 3/7/13:
! NOTE: for now, just call INPHOT on all CPUs.  Try to figure out how
! to MPI broadcast later.  This could be very difficult. (bmy, mlong, 3/7/13)
!    ENDIF
!
!    Broadcast FAST-J inputs to other CPUs
!    CALL GIGC_Inphot_Bcast(  Input_Opt, RC )
!------------------------------------------------------------------------------

    ! Return w/ success
    RC = GC_Success

  END SUBROUTINE GIGC_Init_Simulation
!EOC
END MODULE GIGC_Initialization_Mod
#endif
