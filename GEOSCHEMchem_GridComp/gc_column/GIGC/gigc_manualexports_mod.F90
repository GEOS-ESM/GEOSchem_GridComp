!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: gigc_geos5exports_mod.F90
!
! !DESCRIPTION: Module GIGC\_GEOS5Exports\_Mod serves as the interface
! between the HISTORY configuration file, GEOS-Chem diagnostics defined for 
! GEOS-5 only, and the ESMF Export State. 
!\\
!\\
! !INTERFACE: 
!
MODULE GIGC_GEOS5Exports_Mod
!
! !USES:
!
#include "MAPL_Generic.h"
  USE DiagList_Mod
  USE ErrCode_Mod
  USE Precision_Mod
  USE MAPL_Mod

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: GEOS5Exports_SetServices
  PUBLIC :: GEOS5Exports_SetDataPointers
  PUBLIC :: CopyData2Exports
  PUBLIC :: Destroy_GEOS5diagConfig 
!
! !PRIVATE:
!
  PRIVATE :: Init_Geos5diagconfig
  PRIVATE :: Init_Geos5Export
  PRIVATE :: Init_GEOS5ExportsList
  PRIVATE :: Append_GEOS5ExportsList
  PRIVATE :: Check_GEOS5ExportsList
  PRIVATE :: Print_GEOS5ExportsList
!
! !PUBLIC TYPES
!
  ! History Configuration Object 
  TYPE, PUBLIC :: Geos5diagconfigObj

     CHARACTER(LEN=255)                   :: ROOT ! TODO: needed?
     CHARACTER(LEN=255)                   :: ConfigFileName
     LOGICAL                              :: ConfigFileRead
     TYPE(GEOS5ExportsListObj), POINTER   :: GEOS5ExportsList
     TYPE(DgnList)                        :: DiagList
 
 END TYPE Geos5diagconfigObj
!
! !PRIVATE TYPES
!
  ! History Exports Linked List
  TYPE :: GEOS5ExportsListObj

     TYPE(HistoryExportObj), POINTER :: head
     INTEGER                         :: numExports

  END TYPE GEOS5ExportsListObj

  ! History Export Object
  TYPE :: HistoryExportObj

     CHARACTER(LEN=255)              :: name 
     CHARACTER(LEN=255)              :: metadataID
     CHARACTER(LEN=255)              :: registryID
     CHARACTER(LEN=255)              :: long_name  
     CHARACTER(LEN=255)              :: units       
     INTEGER                         :: vloc
     INTEGER                         :: rank        
     INTEGER                         :: type
     LOGICAL                         :: isMet
     LOGICAL                         :: isChem
     LOGICAL                         :: isDiag
     TYPE(HistoryExportObj), POINTER :: next

     ! Pointers to ESMF Export and GEOS-Chem State
     ! TODO: for now, include all possible data types in the registry. 
     REAL,     POINTER :: ExportData2d(:,:)
     REAL,     POINTER :: ExportData3d(:,:,:)
     REAL(fp), POINTER :: GCStateData0d
     REAL(fp), POINTER :: GCStateData1d(:)
     REAL(fp), POINTER :: GCStateData2d(:,:)
     REAL(fp), POINTER :: GCStateData3d(:,:,:)
     REAL(f4), POINTER :: GCStateData0d_4
     REAL(f4), POINTER :: GCStateData1d_4(:)
     REAL(f4), POINTER :: GCStateData2d_4(:,:)
     REAL(f4), POINTER :: GCStateData3d_4(:,:,:) 
     REAL(f8), POINTER :: GCStateData0d_8
     REAL(f8), POINTER :: GCStateData1d_8(:)
     REAL(f8), POINTER :: GCStateData2d_8(:,:)
     REAL(f8), POINTER :: GCStateData3d_8(:,:,:)
     INTEGER,  POINTER :: GCStateData0d_I
     INTEGER,  POINTER :: GCStateData1d_I(:)
     INTEGER,  POINTER :: GCStateData2d_I(:,:)
     INTEGER,  POINTER :: GCStateData3d_I(:,:,:)

  END TYPE HistoryExportObj
!
! !PUBLIC PARAMETERS
!
  ! Prefix of the species names in the internal state and HISTORY.rc
  CHARACTER(LEN=4), PUBLIC, PARAMETER  :: SPFX = 'TRC_'
!
! !REVISION HISTORY:
!  01 Sep 2017 - E. Lundgren - Initial version
!  08 Mar 2018 - E. Lundgren - Define the internal state prefix expected in 
!                              HISTORY.rc within this module
!EOP
!------------------------------------------------------------------------------
!BOC

CONTAINS
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Init_Geos5diagconfig 
!
! !DESCRIPTION:
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_Geos5diagconfig ( am_I_Root, Geos5diagconfig, configFile, RC )
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
    LOGICAL,             INTENT(IN) :: am_I_Root
    CHARACTER(LEN=*),    INTENT(IN) :: configFile
!
! !OUTPUT PARAMETERS:
!
    TYPE(Geos5diagconfigObj), POINTER :: Geos5diagconfig
    INTEGER, INTENT(OUT)            :: RC 
!
! !REVISION HISTORY:
!  01 Sep 2017 - E. Lundgren - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
    __Iam__('Init_Geos5diagconfig (gigc_geos5exports_mod.F90)')
    ALLOCATE(Geos5diagconfig)
    Geos5diagconfig%ROOT               =  ''
    Geos5diagconfig%ConfigFileName     =  TRIM(configFile)
    Geos5diagconfig%ConfigFileRead     =  .FALSE.
    CALL Init_DiagList( am_I_Root, configFile, Geos5diagconfig%DiagList, RC )
    IF ( RC == GC_FAILURE ) THEN
       ASSERT_(.FALSE.)
       RETURN
    ENDIF
    CALL Print_DiagList( am_I_Root, Geos5diagconfig%DiagList, RC )

    CALL Init_GEOS5ExportsList( am_I_Root, Geos5diagconfig, RC )
    IF ( RC == GC_FAILURE ) THEN
       ASSERT_(.FALSE.)
       RETURN
    ENDIF
    CALL Print_GEOS5ExportsList( am_I_Root, Geos5diagconfig, RC )

  END SUBROUTINE Init_Geos5diagconfig
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Init_GEOS5ExportsList 
!
! !DESCRIPTION: 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_GEOS5ExportsList ( am_I_Root, Geos5diagconfig, RC )
!
! !USES:
!
    USE State_Chm_Mod,    ONLY: Get_Metadata_State_Chm
    USE State_Diag_Mod,   ONLY: Get_Metadata_State_Diag 
    USE State_Met_Mod,    ONLY: Get_Metadata_State_Met
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN)    :: am_I_Root
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(Geos5diagconfigObj), POINTER :: Geos5diagconfig
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT)   :: RC
!
! !REVISION HISTORY:
!  01 Sep 2017 - E. Lundgren - initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER               :: N, rank, vloc, type
    CHARACTER(LEN=255)    :: ErrMsg, desc, units, tag
    LOGICAL               :: isMet, isChem, isDiag, found
    TYPE(HistoryExportObj),  POINTER :: NewHistExp
    TYPE(DgnItem),           POINTER :: current

    ! ================================================================
    ! Init_GEOS5ExportsList begins here
    ! ================================================================
    __Iam__('Init_GEOS5ExportsList (gigc_geos5exports_mod.F90)')

    ! Init
    NewHistExp => NULL()

    ! Create GEOS5ExportsList object
    ALLOCATE(Geos5diagconfig%GEOS5ExportsList)
    Geos5diagconfig%GEOS5ExportsList%numExports = 0
    Geos5diagconfig%GEOS5ExportsList%head => NULL()

    ! Loop over entries in DiagList
    current => Geos5diagconfig%DiagList%head
    DO WHILE ( ASSOCIATED( current ) )

       ! Skip State_Chm%Species entries since in internal state
       ! TODO: In GCHP this would appear with prefix stored in SPFX
       !       Need to make GCC and GCHP more consistent in future
       IF ( INDEX( current%name,  TRIM(SPFX) ) > 0 ) THEN
          current => current%next
          CYCLE
       ENDIF

       ! Skip emissions diagnostics since handled by HEMCO
       ! Will need to revisit this since name may change
       IF ( INDEX( current%name,  'EMIS' ) == 1 ) THEN
          current => current%next
          CYCLE
       ENDIF

       ! Check history exports list to see if already added (unless wildcard)
       ! TODO: consider making the call a function that returns a logical
       IF ( .NOT. current%isWildcard ) THEN
          CALL Check_GEOS5ExportsList( am_I_Root, current%name,           &
                                         Geos5diagconfig%GEOS5ExportsList,  &
                                         found, RC                         )
          IF ( found ) THEN
             current => current%next
             CYCLE
          ENDIF
       ENDIF

       ! Get metadata using metadataID and state
       ! If isTagged, then append to description
       ! If isWildcard, shouldn't get here
       ! The name of the export is simply name
       Found = .TRUE.
       isMet  = .FALSE.
       isChem = .FALSE.
       isDiag = .FALSE.
       IF ( TRIM(current%state) == 'MET' ) THEN
          isMet = .TRUE.
          CALL Get_Metadata_State_Met( am_I_Root, current%metadataID,     &
                                       Found, RC, desc=desc, units=units, &
                                       rank=rank, type=type, vloc=vloc )
          ! TODO: need to add found to outputs of get_metadata_state_met
       ELSEIF ( TRIM(current%state) == 'CHEM' ) THEN
          isCHEM = .TRUE.
          CALL Get_Metadata_State_Chm( am_I_Root, current%metadataID,     &
                                       Found, RC, desc=desc, units=units, &
                                       rank=rank, type=type, vloc=vloc )
       ELSEIF ( TRIM(current%state) == 'DIAG' ) THEN
          isDIAG = .TRUE.
          CALL Get_Metadata_State_Diag( am_I_Root, current%metadataID,     &
                                        Found, Rc, desc=desc, units=units, &
                                        rank=rank, type=type, vloc=vloc )
       ELSEIF ( TRIM(current%state) == 'GEOS5' ) THEN
          ! Skip it
          current => current%next
          CYCLE
       ELSE
          ErrMsg = "Unknown state of item " // TRIM(current%name) // &
                   " in DiagList: " // TRIM(current%state)
          EXIT
       ENDIF
       IF ( Found .eqv. .FALSE. ) THEN
          ErrMsg = "Metadata not found for " // TRIM(current%name)
          EXIT
       ENDIF

       ! If wildcard is present
       IF ( current%isWildcard ) THEN
          ! Do nothing. This should never happen at this point since
          ! Init_DiagList will exit with an error if wildcard is
          ! encountered in HISTORY.rc while compiling with ESMF_.

          ! When it comes time to implement, create exports in a loop,
          ! either for all species or for advected species only. Include 
          ! a check that the export was not already created. Loop over 
          ! AdvNames if wildcard is ADV. Loop over SpecNames for all other 
          ! cases, passing not found = OK so that not all are necessarily 
          ! output. Later on, after species database is initialized, exports 
          ! for only species in the specific wildcard will be associated 
          ! with data and thus included in the output file.

          ! If the meantime, skip wildcards if it gets here.
          current => current%next
          CYCLE
       ENDIF

       ! If this item is for a specific tag, append description. 
       ! This will need revisiting since there may be tag-dependent
       ! strings to append to long names
       IF ( current%isTagged ) THEN
          desc = TRIM(desc) // " for " // TRIM(current%tag)
       ENDIF

       ! Create a new HistoryExportObj object
       CALL Init_HistoryExport( am_I_Root, NewHistExp,         &
                                name=current%name,             &
                                metadataID=current%metadataID, &
                                registryID=current%registryID, &
                                long_name=desc,                &
                                units=units,                   &
                                vloc=vloc,                     &
                                rank=rank,                     &
                                type=type,                     &
                                isMet=isMet,                   &
                                isChem=isChem,                 &
                                isDiag=isDiag,                 &
                                RC=RC )
       IF ( RC == GC_FAILURE ) THEN
          ErrMsg = "History export init fail for " // TRIM(current%name)
          EXIT
       ENDIF
       
       ! Add new HistoryExportObj to linked list
       CALL Append_GEOS5ExportsList( am_I_Root,     NewHistExp, &
                                       Geos5diagconfig, RC       )
       IF ( RC == GC_FAILURE ) THEN
          ErrMsg = "History export append fail for " // TRIM(current%name)
          EXIT
       ENDIF

       ! Set up for next item in DiagList
       current => current%next

    ENDDO
    current => NULL()

    IF ( RC == GC_SUCCESS ) THEN
       Geos5diagconfig%ConfigFileRead = .TRUE.
    ELSE
       CALL GC_ERROR( ErrMsg, RC, Iam )
       ASSERT_(.FALSE.)
       RETURN
    ENDIF

  END SUBROUTINE Init_GEOS5ExportsList
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Init_HistoryExport 
!
! !DESCRIPTION: 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_HistoryExport ( am_I_Root,  NewHistExp, name,         &
                                  metadataID, registryID, long_name,    &
                                  units,      vloc,       rank,         &
                                  type,       isMet,      isChem,       &
                                  isDiag,     RC  )
!
! !INPUT PARAMETERS:
!
    LOGICAL,             INTENT(IN)    :: am_I_Root
!
! !OUTPUT PARAMETERS:
!
    TYPE(HistoryExportObj), POINTER :: NewHistExp
    CHARACTER(LEN=*), OPTIONAL      :: name
    CHARACTER(LEN=*), OPTIONAL      :: metadataID
    CHARACTER(LEN=*), OPTIONAL      :: registryID
    CHARACTER(LEN=*), OPTIONAL      :: long_name
    CHARACTER(LEN=*), OPTIONAL      :: units
    INTEGER,          OPTIONAL      :: vloc 
    INTEGER,          OPTIONAL      :: rank
    INTEGER,          OPTIONAL      :: type
    LOGICAL,          OPTIONAL      :: isMet
    LOGICAL,          OPTIONAL      :: isChem
    LOGICAL,          OPTIONAL      :: isDiag
    INTEGER,          OPTIONAL      :: RC
!
! !REVISION HISTORY:
!  01 Sep 2017 - E. Lundgren - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
    __Iam__('Init_HistoryExport (gigc_geos5exports_mod.F90)')
    ALLOCATE(NewHistExp)
    NewHistExp%name        = TRIM(name)
    NewHistExp%metadataID  = TRIM(metadataID)
    NewHistExp%registryID  = TRIM(registryID)
    NewHistExp%long_name   = TRIM(long_name)
    NewHistExp%units       = TRIM(units)
    NewHistExp%vloc        = vloc
    NewHistExp%rank        = rank
    NewHistExp%type        = type
    NewHistExp%isMet       = isMet
    NewHistExp%isChem      = isChem
    NewHistExp%isDiag      = isDiag
    NewHistExp%next          => NULL()
    NewHistExp%ExportData2d  => NULL()
    NewHistExp%ExportData3d  => NULL()
    NewHistExp%GCStateData0d => NULL()
    NewHistExp%GCStateData1d => NULL()
    NewHistExp%GCStateData2d => NULL()
    NewHistExp%GCStateData3d => NULL()
    NewHistExp%GCStateData0d_4 => NULL()
    NewHistExp%GCStateData1d_4 => NULL()
    NewHistExp%GCStateData2d_4 => NULL()
    NewHistExp%GCStateData3d_4 => NULL()
    NewHistExp%GCStateData0d_8 => NULL()
    NewHistExp%GCStateData1d_8 => NULL()
    NewHistExp%GCStateData2d_8 => NULL()
    NewHistExp%GCStateData3d_8 => NULL()
    NewHistExp%GCStateData0d_I => NULL()
    NewHistExp%GCStateData1d_I => NULL()
    NewHistExp%GCStateData2d_I => NULL()
    NewHistExp%GCStateData3d_I => NULL()

  END SUBROUTINE Init_HistoryExport
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Append_GEOS5ExportsList 
!
! !DESCRIPTION: 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Append_GEOS5ExportsList ( am_I_Root,     HistoryExport, &
                                         Geos5diagconfig, RC        )
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN)    :: am_I_Root
    TYPE(HistoryExportObj), POINTER :: HistoryExport
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(Geos5diagconfigObj), POINTER :: Geos5diagconfig
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT)   :: RC
!
! !REVISION HISTORY:
!  01 Sep 2017 - E. Lundgren - initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    TYPE(HistoryExportObj),  POINTER :: NewHistExp

    ! ================================================================
    ! Append_GEOS5ExportsList begins here
    ! ================================================================
    __Iam__('Append_GEOS5ExportsList (gigc_geos5exports_mod.F90)')
    RC = GC_SUCCESS

    ! Add new object to the beginning of the linked list
    HistoryExport%next => Geos5diagconfig%GEOS5ExportsList%head
    Geos5diagconfig%GEOS5ExportsList%head => HistoryExport

    ! Update # of list items
    Geos5diagconfig%GEOS5ExportsList%numExports = &
         Geos5diagconfig%GEOS5ExportsList%numExports + 1

  END SUBROUTINE Append_GEOS5ExportsList
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Check_GEOS5ExportsList
!
! !DESCRIPTION:
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Check_GEOS5ExportsList ( am_I_Root, name,  &
                                        ExportsList, found, RC )
!
! !INPUT PARAMETERS:
!
    LOGICAL,           INTENT(IN)        :: am_I_Root
    CHARACTER(LEN=*),  INTENT(IN)        :: name
    TYPE(GEOS5ExportsListObj), POINTER :: ExportsList
!
! !OUTPUT PARAMETERS:
!
    LOGICAL, INTENT(OUT)               :: found
    INTEGER, INTENT(OUT)               :: RC 
!
! !REVISION HISTORY:
!  12 Sep 2017 - E. Lundgren - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    TYPE(HistoryExportObj), POINTER :: current

    __Iam__('Check_GEOS5ExportsList (gigc_geos5exports_mod.F90)')
    RC = GC_SUCCESS

    ! Assume not found
    found = .False.

    current => ExportsList%head
    DO WHILE ( ASSOCIATED( current ) )
       IF ( current%name == name ) THEN
          found = .TRUE.
          RETURN
       ENDIF
       current => current%next    
    ENDDO
    current => NULL()

  END SUBROUTINE Check_GEOS5ExportsList
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !!
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GEOS5Exports_SetServices
!
  ! !DESCRIPTION:
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GEOS5Exports_SetServices( am_I_Root, config_file, GC,  &
                                       Geos5diagconfig,  RC )
!
! !USES:
!
    USE ESMF, ONLY : ESMF_GridComp
    USE Registry_Params_Mod
!
! !INPUT PARAMETERS:
!
    LOGICAL,             INTENT(IN)    :: am_I_Root
    CHARACTER(LEN=*),    INTENT(IN)    :: config_file
!
! !INPUT AND OUTPUT PARAMETERS:
!
    TYPE(ESMF_GridComp), INTENT(INOUT) :: GC             ! Gridded Component
!
! !OUTPUT PARAMETERS:
!
    TYPE(Geos5diagconfigObj), POINTER    :: Geos5diagconfig  ! History config object
    INTEGER,             INTENT(OUT)   :: RC
!
! !REMARKS:
!  !
! !REVISION HISTORY: 
!  01 Sep 2017 - E. Lundgren - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=255)                   :: ErrMsg
    TYPE(HistoryExportObj),      POINTER :: current

    ! ================================================================
    ! GEOS5Exports_SetServices begins here
    ! ================================================================

    ! For MAPL/ESMF error handling (defines Iam and STATUS)
    __Iam__('GEOS5Exports_SetServices (gigc_geos5exports_mod.F90)')

    ! Create a config object if it does not already exist
    IF ( .NOT. ASSOCIATED(Geos5diagconfig) ) THEN
       CALL Init_Geos5diagconfig( am_I_Root, Geos5diagconfig, config_file, RC )
       IF ( RC == GC_FAILURE ) THEN
          ASSERT_(.FALSE.)
          RETURN
       ENDIF
    ENDIF

    ! Loop over the History Exports list to add one export per item
    IF ( am_I_Root ) THEN
       WRITE(6,*) " "
       WRITE(6,*) "Adding history variables to GCHP Export State:"
    ENDIF
    current => Geos5diagconfig%GEOS5ExportsList%head
    DO WHILE ( ASSOCIATED( current ) )
       ! Create an export for this item
       IF ( current%rank == 3 ) THEN
          IF ( current%vloc == VLocationCenter ) THEN
             IF ( am_I_Root ) PRINT *, "adding export: ", TRIM(current%name)
             CALL MAPL_AddExportSpec(GC,                                     &
                                     SHORT_NAME = TRIM(current%name),        &
                                     LONG_NAME  = TRIM(current%long_name),   &
                                     UNITS      = TRIM(current%units),       &
                                     DIMS       = MAPL_DimsHorzVert,         &
                                     VLOCATION  = MAPL_VLocationCenter,      &
                                     RC         = RC                         )
          IF ( RC == GC_FAILURE ) THEN
             ErrMsg =  "Problem adding 3D export for " // TRIM(current%name)
             EXIT
          ENDIF
         ELSEIF ( current%vloc == VLocationEdge ) THEN
            CALL MAPL_AddExportSpec(GC,                                     &
                                    SHORT_NAME = TRIM(current%name), &
                                    LONG_NAME  = TRIM(current%long_name),   &
                                    UNITS      = TRIM(current%units),       &
                                    DIMS       = MAPL_DimsHorzVert,         &
                                    VLOCATION  = MAPL_VLocationEdge,        &
                                    RC         = STATUS                    )
         ELSE
            IF ( am_I_Root ) THEN
               PRINT *, "Unknown vertical location for ", &
                        TRIM(current%name)
            ENDIF
         ENDIF
       ELSEIF ( current%rank == 2 ) THEN
          CALL MAPL_AddExportSpec(GC,                                     &
                                  SHORT_NAME = TRIM(current%name), &
                                  LONG_NAME  = TRIM(current%long_name),   &
                                  UNITS      = TRIM(current%units),       &
                                  DIMS       = MAPL_DimsHorzOnly,         &
                                  RC         = RC                        )
          IF ( RC == GC_FAILURE ) THEN
             ErrMsg =  "Problem adding 2D export for " // TRIM(current%name)
             EXIT
          ENDIF
       ELSE
          ErrMsg = "Problem adding export for " // TRIM(current%name) // &
                   ". Rank is only implemented for 2 or 3!"
          RC = GC_FAILURE
          EXIT
       ENDIF

       current => current%next    
    ENDDO
    current => NULL()

    IF ( RC == GC_FAILURE ) THEN
       CALL GC_ERROR( ErrMsg, RC, Iam )
       ASSERT_(.FALSE.)
       RETURN
    ENDIF
    
  END SUBROUTINE GEOS5Exports_SetServices
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: CopyData2Exports
!
! !DESCRIPTION:
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CopyData2Exports( am_I_Root, Input_Opt, Geos5diagconfig, RC )
!
! !USES:
!
  USE HCOI_GC_Main_Mod, ONLY : HCOI_GC_WriteDiagn
  USE Input_Opt_Mod,    ONLY : OptInput
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)    :: am_I_Root
    TYPE(OptInput), INTENT(IN)    :: Input_Opt
!
! !INPUT AND OUTPUT PARAMETERS:
!
    TYPE(Geos5diagconfigObj), POINTER :: Geos5diagconfig  ! History config object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,            INTENT(OUT) :: RC
!
! !REMARKS:
!  !
! !REVISION HISTORY: 
!  01 Sep 2017 - E. Lundgren - Initial version
!  02 Nov 2017 - E. Lundgren - Copy HEMCO data to emissions exports
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=255)              :: ErrMsg
    TYPE(HistoryExportObj), POINTER :: current

    ! ================================================================
    ! CopyData2Exports begins here
    ! ================================================================
    __Iam__('CopyData2Exports (gigc_geos5exports_mod.F90)')

    ! Loop over the History Exports list
    current => Geos5diagconfig%GEOS5ExportsList%head
    DO WHILE ( ASSOCIATED( current ) )

       IF ( current%rank == 2 ) THEN
          IF ( ASSOCIATED ( current%GCStateData2d ) ) THEN
             current%ExportData2d = current%GCStateData2d
          ELSEIF ( ASSOCIATED ( current%GCStateData2d_4 ) ) THEN
             current%ExportData2d = current%GCStateData2d_4
          ELSEIF ( ASSOCIATED ( current%GCStateData2d_8 ) ) THEN
             current%ExportData2d = current%GCStateData2d_8
          ELSEIF ( ASSOCIATED ( current%GCStateData2d_I ) ) THEN
             ! Convert integer to float (integers not allowed in MAPL exports)
             current%ExportData2d = FLOAT(current%GCStateData2d_I)
          ELSE
             RC = GC_FAILURE
             ErrMsg = "No GC 2D pointer found for " // TRIM(current%name)
             EXIT
          ENDIF
       ELSEIF ( current%rank == 3 ) THEN
          IF ( ASSOCIATED ( current%GCStateData3d ) ) THEN
             current%ExportData3d = current%GCStateData3d
          ELSEIF ( ASSOCIATED ( current%GCStateData3d_4 ) ) THEN
             current%ExportData3d = current%GCStateData3d_4
          ELSEIF ( ASSOCIATED ( current%GCStateData3d_8 ) ) THEN
             current%ExportData3d = current%GCStateData3d_8
          ELSEIF ( ASSOCIATED ( current%GCStateData3d_I ) ) THEN
             current%ExportData3d = FLOAT(current%GCStateData3d_I)
          ELSE
             RC = GC_FAILURE
             ErrMsg = "No GC 3D pointer found for " // TRIM(current%name)
             EXIT
          ENDIF
       ENDIF

       current => current%next    
    ENDDO
    current => NULL()

    ! Error handling
    IF ( RC == GC_FAILURE ) THEN
       CALL GC_ERROR( ErrMsg, RC, Iam )
       ASSERT_(.FALSE.)
       RETURN
    ENDIF

    ! Copy emissions data to MAPL exports via HEMCO
    CALL HCOI_GC_WriteDiagn( am_I_Root, Input_Opt, .FALSE., RC )
    IF ( RC == GC_FAILURE ) THEN
       ErrMsg = "Error copying emissions data to MAPL via HEMCO"
       CALL GC_ERROR( ErrMsg, RC, Iam )
       ASSERT_(.FALSE.)
       RETURN
    ENDIF
    
  END SUBROUTINE CopyData2Exports
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Print_GEOS5ExportsList
!
! !DESCRIPTION:
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Print_GEOS5ExportsList( am_I_Root, Geos5diagconfig, RC )
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
    LOGICAL,             INTENT(IN) :: am_I_Root
!
! !INPUT AND OUTPUT PARAMETERS:
!
    TYPE(Geos5diagconfigObj), POINTER :: Geos5diagconfig  ! History config object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,            INTENT(OUT) :: RC
!
! !REMARKS:
!  !
! !REVISION HISTORY: 
!  01 Sep 2017 - E. Lundgren - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    TYPE(HistoryExportObj), POINTER :: current

    ! ================================================================
    ! Print_GEOS5ExportsList begins here
    ! ================================================================
    __Iam__('Print_GEOS5ExportsList (gigc_geos5exports_mod.F90)')

    RC = GC_SUCCESS

    ! Loop over the History Exports list
    current => Geos5diagconfig%GEOS5ExportsList%head
    IF ( am_I_Root ) PRINT *, '==========================='
    IF ( am_I_Root ) PRINT *, 'History Exports List:'
    IF ( am_I_Root ) PRINT *, ' '
    DO WHILE ( ASSOCIATED( current ) )
       IF ( am_I_Root ) THEN
          PRINT *, "Name:        ",   TRIM(current%name) 
          PRINT *, " MetadataID: ",   TRIM(current%metadataID)
          PRINT *, " RegistryID: ",   TRIM(current%registryID)
          PRINT *, " Long name:  ",   TRIM(current%long_name)  
          PRINT *, " Units:      ",   TRIM(current%units)       
          PRINT *, " Vert loc:   ",   current%vloc
          PRINT *, " Rank:       ",   current%rank        
          PRINT *, " Type:       ",   current%type
          PRINT *, " isMet:      ",   current%isMet
          PRINT *, " isChem:     ",   current%isChem
          PRINT *, " isDiag:     ",   current%isDiag
          PRINT *, " "
       ENDIF
       current => current%next    
    ENDDO
    IF ( am_I_Root ) PRINT *, '==========================='
    current => NULL()
    
  END SUBROUTINE Print_GEOS5ExportsList
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GEOS5Exports_SetDataPointers
!
! !DESCRIPTION:
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GEOS5Exports_SetDataPointers( am_I_Root,     EXPORT,    &
                                             Geos5diagconfig, State_Chm, &
                                             State_Diag,    State_Met, &
                                             RC                       )
!
! !USES:
!
    USE ESMF,           ONLY : ESMF_State
    USE Registry_Mod,   ONLY : Registry_Lookup
    USE State_Chm_Mod,  ONLY : ChmState
    USE State_Diag_Mod, ONLY : DgnState
    USE State_Met_Mod,  ONLY : MetState
!
! !INPUT PARAMETERS:
!
    LOGICAL,         INTENT(IN)     :: am_I_Root
!
! !INPUT AND OUTPUT PARAMETERS:
!
    TYPE(ESMF_State), INTENT(INOUT), TARGET :: EXPORT ! Export state obj
    TYPE(Geos5diagconfigObj), POINTER :: Geos5diagconfig  ! History config obj
    TYPE(ChmState),   INTENT(INOUT) :: State_Chm      ! Chemistry State obj
    TYPE(MetState),   INTENT(INOUT) :: State_Met      ! Meteorology State obj
    TYPE(DgnState),   INTENT(INOUT) :: State_Diag     ! Diagnostics State obj
!
! !OUTPUT PARAMETERS:
!
    INTEGER,             INTENT(OUT)   :: RC
!
! !REMARKS:
!  !
! !REVISION HISTORY: 
!  01 Sep 2017 - E. Lundgren - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=255)              :: ErrMsg
    TYPE(HistoryExportObj), POINTER :: current

    ! ================================================================
    ! GEOS5Exports_SetDataPointers begins here
    ! ================================================================
    __Iam__('GEOS5Exports_SetDataPointers')

    IF ( am_I_Root ) THEN
       WRITE(6,*) " "
       WRITE(6,*) "Setting history variable pointers to GC and Export States:"
    ENDIF

    ! Loop over the History Exports list
    current => Geos5diagconfig%GEOS5ExportsList%head
    DO WHILE ( ASSOCIATED( current ) )

       ! Get pointer to GC state data
       !IF ( am_I_Root ) WRITE(6,*) current%name
       IF ( current%isMET ) THEN
          CALL Registry_Lookup( am_I_Root = am_I_Root,               &
                                Registry  = State_Met%Registry,      &
                                State     = State_Met%State,         &
                                Variable  = current%registryID,      &
                                Ptr2d     = current%GCStateData2d,   &
                                Ptr2d_4   = current%GCStateData2d_4, &
                                Ptr2d_8   = current%GCStateData2d_8, &
                                Ptr2d_I   = current%GCStateData2d_I, &
                                Ptr3d     = current%GCStateData3d,   &
                                Ptr3d_4   = current%GCStateData3d_4, &
                                Ptr3d_8   = current%GCStateData3d_8, &
                                Ptr3d_I   = current%GCStateData3d_I, &
                                RC        = RC                      )
       ELSEIF ( current%isChem ) THEN
          CALL Registry_Lookup( am_I_Root = am_I_Root,               &
                                Registry  = State_Chm%Registry,      &
                                State     = State_Chm%State,         &
                                Variable  = current%registryID,      &
                                Ptr2d     = current%GCStateData2d,   &
                                Ptr2d_4   = current%GCStateData2d_4, &
                                Ptr2d_8   = current%GCStateData2d,   &
                                Ptr2d_I   = current%GCStateData2d_I, &
                                Ptr3d     = current%GCStateData3d,   &
                                Ptr3d_4   = current%GCStateData3d_4, &
                                Ptr3d_8   = current%GCStateData3d,   &
                                Ptr3d_I   = current%GCStateData3d_I, &
                                RC        = RC                      )
       ELSEIF ( current%isDiag ) THEN
          CALL Registry_Lookup( am_I_Root = am_I_Root,               &
                                Registry  = State_Diag%Registry,     &
                                State     = State_Diag%State,        &
                                Variable  = current%registryID,      &
                                Ptr2d     = current%GCStateData2d,   &
                                Ptr2d_4   = current%GCStateData2d_4, &
                                Ptr2d_8   = current%GCStateData2d,   &
                                Ptr2d_I   = current%GCStateData2d_I, &
                                Ptr3d     = current%GCStateData3d,   &
                                Ptr3d_4   = current%GCStateData3d_4, &
                                Ptr3d_8   = current%GCStateData3d,   &
                                Ptr3d_I   = current%GCStateData3d_I, &
                                RC        = RC                      )
       ENDIF
       IF ( RC == GC_FAILURE ) THEN
          ErrMsg = "Registry pointer not found for " // TRIM(current%name) // &
                   ". Check that the tag (e.g. species) is valid "         // &
                   "for this diagnostic."
          EXIT
       ENDIF

       ! For MAPL export, need to pass a pointer of the right dimension
       IF ( current%rank == 2 ) THEN
          CALL MAPL_GetPointer ( EXPORT, current%ExportData2d, &
                                 current%name, __RC__ )
       ELSEIF ( current%rank == 3 ) THEN
          CALL MAPL_GetPointer ( EXPORT, current%ExportData3d, &
                                 current%name, __RC__ )
       ENDIF

       IF ( Am_I_Root) THEN
          WRITE(6,*) TRIM(current%name)
       ENDIF

       current => current%next    
    ENDDO
    current => NULL()

    IF ( RC == GC_FAILURE ) THEN
       CALL GC_ERROR( ErrMsg, RC, Iam )
       ASSERT_(.FALSE.)
       RETURN
    ENDIF
    
  END SUBROUTINE GEOS5Exports_SetDataPointers
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Destroy_Geos5diagconfig 
!
! !DESCRIPTION: Subroutine Destroy_Geos5diagconfig deallocates a Geos5diagconfig
!  object and all of its member objects including the linked list of
!  HistoryExport objects.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Destroy_Geos5diagconfig ( am_I_Root, Geos5diagconfig, RC )
!
! !INPUT PARAMETERS:
!
    LOGICAL,            INTENT(IN)    :: am_I_Root     ! root CPU?
    TYPE(Geos5diagconfigObj), POINTER   :: Geos5diagconfig
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,            INTENT(INOUT) :: RC            ! Success?
!
! !REVISION HISTORY:
!  01 Sep 2017 - E. Lundgren - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    TYPE(HistoryExportObj), POINTER :: current
    TYPE(HistoryExportObj), POINTER :: next

    ! ================================================================
    ! Destroy_Geos5diagconfig begins here
    ! ================================================================
    __Iam__('Destroy_Geos5diagconfig (gigc_geos5exports_mod.F90)')

    current => NULL()
    next => NULL()

    ! Destroy each item in the linked list of HistoryExport objects
    current => Geos5diagconfig%GEOS5ExportsList%head
    IF ( ASSOCIATED( current ) ) next => current%next
    DO WHILE ( ASSOCIATED( current ) )
       DEALLOCATE( current, STAT=RC )
       ASSERT_( RC == GC_SUCCESS )
       IF ( .NOT. ASSOCIATED ( next ) ) EXIT
       current => next
       next => current%next
    ENDDO

    ! Deallocate the GEOS5ExportsList object
    DEALLOCATE( Geos5diagconfig%GEOS5ExportsList, STAT=RC )
    ASSERT_( RC == GC_SUCCESS )

    ! Deallocate the Geos5diagconfig object
    DEALLOCATE( Geos5diagconfig, STAT=RC )
    ASSERT_( RC == GC_SUCCESS )

    ! Final cleanup
    current => NULL()
    next    => NULL()

  END SUBROUTINE Destroy_Geos5diagconfig
!EOC
!!------------------------------------------------------------------------------
!!          Harvard University Atmospheric Chemistry Modeling Group            !
!!------------------------------------------------------------------------------
!!BOP
!!
!! !IROUTINE: StateDiag2Export 
!!
!! !DESCRIPTION: 
!!\\
!!\\
!! !INTERFACE:
!!
!  SUBROUTINE StateDiag2Export ( am_I_Root,  Input_Opt, State_Met, State_Chm, &
!                                State_Diag, Diag_List, EXPORT,    IsChemTime, &
!                                DoDryDep,   Phase, IM, JM, LM, RC )
!!
!! !USES:
!!
!    USE ErrCode_Mod
!    USE Precision_Mod
!    USE Input_Opt_Mod,   ONLY : OptInput
!    USE State_Met_Mod,   ONLY : MetState
!    USE State_Chm_Mod,   ONLY : ChmState, Ind_
!    USE State_Diag_Mod,  ONLY : DgnState
!    USE DiagList_Mod,    ONLY : DgnList 
!    USE DIAG_MOD,        ONLY : AD21
!    USE DIAG_OH_MOD,     ONLY : OH_MASS, AIR_MASS 
!!
!! !INPUT PARAMETERS:
!!
!    LOGICAL,          INTENT(IN)    :: am_I_Root     ! Are we on the root CPU?
!    TYPE(OptInput),   INTENT(IN)    :: Input_Opt     ! Input Options object
!    LOGICAL,          INTENT(IN)    :: IsChemTime    ! Chemistry time? 
!    LOGICAL,          INTENT(IN)    :: DoDryDep      ! Do dry deposition? 
!    INTEGER,          INTENT(IN)    :: Phase         ! Run phase 
!    INTEGER,          INTENT(IN)    :: IM, JM, LM    ! Grid dimensions 
!!
!! !INPUT/OUTPUT PARAMETERS:
!!
!    TYPE(MetState),   INTENT(INOUT) :: State_Met     ! Meteorology State object
!    TYPE(ChmState),   INTENT(INOUT) :: State_Chm     ! Chem state object
!    TYPE(DgnState),   INTENT(INOUT) :: State_Diag    ! Diagnostics state object
!    TYPE(DgnList),    INTENT(INOUT) :: Diag_List     ! Diagnostics List 
!    TYPE(ESMF_State), INTENT(INOUT) :: Export
!!
!! !OUTPUT PARAMETERS:
!!
!    INTEGER,          INTENT(OUT)   :: RC            ! Success or failure
!!
!! !REMARKS:
! This was Christoph's implementation for using diaglist. It was called
! in gigc_chunk_run after the call to HCOI_GC_WriteDiagn. 
!    ! Pass arrays from State_Diag to MAPL Export
!    CALL StateDiag2Export ( am_I_Root, Input_Opt, State_Met, State_Chm, &
!                            State_Diag, Diag_List, HcoState%EXPORT, &
!                            DoChem, DoDryDep, Phase, IM, JM, LM, RC )
!    ASSERT_(RC==GC_SUCCESS)
!! 
!! !REVISION HISTORY: 
!!  22 Nov 2017 - C. Keller   - Initial version
!!EOP
!!------------------------------------------------------------------------------
!!BOC
!!
!! !LOCAL VARIABLES:
!!
!    INTEGER                    :: I, L, N
!    CHARACTER(LEN=  3)         :: III
!    CHARACTER(LEN=255)         :: DiagName
!    REAL, POINTER              :: Ptr2D(:,:)  
!    REAL, POINTER              :: Ptr3D(:,:,:)  
!    REAL(f4), POINTER          :: PM25ptr(:,:,:)
!    REAL, ALLOCATABLE, TARGET  :: Diag2D(:,:)
!
!    INTEGER                    :: STATUS
!    CHARACTER(LEN=ESMF_MAXSTR) :: Iam
!
!    !=======================================================================
!    ! StateDiag2Export begins here 
!    !=======================================================================
!
!    ! Error trap
!    Iam = 'StateDiag2Export (gigc_chunk_mod.F90)'
! 
!    ! Assume success
!    RC = GC_SUCCESS
!
!    ! DryDepVel
!    IF ( ASSOCIATED(State_Diag%DryDepVel) .AND. DoDryDep ) THEN
!       DO I = 1, State_Chm%nDryDep
!          N        = State_Chm%Map_DryDep(I)
!          DiagName = 'DryDepVel_'//TRIM(State_Chm%SpcData(N)%Info%Name)
!          CALL MAPL_GetPointer( Export, Ptr2D, TRIM(DiagName), NotFoundOk=.TRUE., __RC__ )
!          IF ( ASSOCIATED(Ptr2D) ) THEN
!             Ptr2D(:,:) = State_Diag%DryDepVel(:,:,I)
!             Ptr2D => NULL()
!          END IF
!       ENDDO
!    ENDIF
!
!!    ! Aerodynamic resistance --> now computed in GCC_GridCompMod.F90
!!    IF ( ASSOCIATED(State_Diag%DryDepRa2m) .AND. DoDryDep ) THEN
!!       DiagName = 'DryDepRa2m_GCC'
!!       CALL MAPL_GetPointer( Export, Ptr2D, TRIM(DiagName), NotFoundOk=.TRUE., __RC__ )
!!       IF ( ASSOCIATED(Ptr2D) ) THEN
!!          Ptr2D(:,:) = State_Diag%DryDepRa2m(:,:)
!!          Ptr2D => NULL()
!!       END IF
!!    ENDIF
!!    IF ( ASSOCIATED(State_Diag%DryDepRa10m) .AND. DoDryDep ) THEN
!!       DiagName = 'DryDepRa10m_GCC'
!!       CALL MAPL_GetPointer( Export, Ptr2D, TRIM(DiagName), NotFoundOk=.TRUE., __RC__ )
!!       IF ( ASSOCIATED(Ptr2D) ) THEN
!!          Ptr2D(:,:) = State_Diag%DryDepRa10m(:,:)
!!          Ptr2D => NULL()
!!       END IF
!!    ENDIF
!
!    ! DryDepFlux_Mix
!    IF ( ASSOCIATED(State_Diag%DryDepMix) .AND. DoDryDep ) THEN
!       DO I = 1, State_Chm%nDryDep
!          N        = State_Chm%Map_DryDep(I)
!          DiagName = 'DryDep_'//TRIM(State_Chm%SpcData(N)%Info%Name)
!          CALL MAPL_GetPointer( Export, Ptr2D, TRIM(DiagName), NotFoundOk=.TRUE., __RC__ )
!          IF ( ASSOCIATED(Ptr2D) ) THEN
!             Ptr2D(:,:) = State_Diag%DryDep(:,:,I)
!             Ptr2D => NULL()
!          END IF
!       ENDDO
!    ENDIF
!
!    ! Monin Obukhov length
!    IF ( ASSOCIATED(State_Diag%MoninObukhov) ) THEN
!       DiagName = 'MoninObukhov'
!       CALL MAPL_GetPointer( Export, Ptr2D, TRIM(DiagName), NotFoundOk=.TRUE., __RC__ )
!       IF ( ASSOCIATED(Ptr2D) ) THEN
!          Ptr2D(:,:) = State_Diag%MoninObukhov
!          Ptr2D => NULL()
!       ENDIF
!    ENDIF
!
!    ! Wet deposition in convection (collapse to 2D)
!    IF ( Phase /= 2 ) THEN
!       IF ( ASSOCIATED(State_Diag%WetLossConv) ) THEN
!          DO I = 1, State_Chm%nWetDep
!             N = State_Chm%Map_WetDep(I)
!             DiagName = 'WetLossConv_'//TRIM(State_Chm%SpcData(N)%Info%Name)
!             CALL MAPL_GetPointer( Export, Ptr2D, TRIM(DiagName), NotFoundOk=.TRUE., __RC__ )
!             IF ( ASSOCIATED(Ptr2D) ) THEN
!                Ptr2D(:,:) = SUM(State_Diag%WetLossConv(:,:,:,I),DIM=3)
!                Ptr2D(:,:) = Ptr2D(:,:) / State_Met%AREA_M2(:,:,1)
!                Ptr2D => NULL()
!             END IF
!             DiagName = 'WetLossConv3D_'//TRIM(State_Chm%SpcData(N)%Info%Name)
!             CALL MAPL_GetPointer( Export, Ptr3D, TRIM(DiagName), NotFoundOk=.TRUE., __RC__ )
!             IF ( ASSOCIATED(Ptr3D) ) THEN
!                DO L = 1,LM
!                   Ptr3D(:,:,LM-L+1) = State_Diag%WetLossConv(:,:,L,I) / State_Met%AREA_M2(:,:,1)
!                ENDDO
!                Ptr3D => NULL()
!             END IF
!          ENDDO
!       ENDIF
!       ! Reset diagnostics
!       !!!State_Diag%WetLossConv(:,:,:,:) = 0.0_f4
!    ENDIF
!
!    ! Wet deposition in large scale precip/washout(collapse to 2D)
!    IF ( Phase /= 1 ) THEN
!       IF ( ASSOCIATED(State_Diag%WetLossLS) ) THEN
!          DO I = 1, State_Chm%nWetDep
!             N = State_Chm%Map_WetDep(I)
!             DiagName = 'WetLossLS_'//TRIM(State_Chm%SpcData(N)%Info%Name)
!             CALL MAPL_GetPointer( Export, Ptr2D, TRIM(DiagName), NotFoundOk=.TRUE., __RC__ )
!             IF ( ASSOCIATED(Ptr2D) ) THEN
!                Ptr2D(:,:) = SUM(State_Diag%WetLossLS(:,:,:,I),DIM=3)
!                Ptr2D(:,:) = Ptr2D(:,:) / State_Met%AREA_M2(:,:,1)
!                Ptr2D => NULL()
!             END IF
!             DiagName = 'WetLossLS3D_'//TRIM(State_Chm%SpcData(N)%Info%Name)
!             CALL MAPL_GetPointer( Export, Ptr3D, TRIM(DiagName), NotFoundOk=.TRUE., __RC__ )
!             IF ( ASSOCIATED(Ptr3D) ) THEN
!                DO L = 1,LM
!                   Ptr3D(:,:,LM-L+1) = State_Diag%WetLossLS(:,:,L,I) / State_Met%AREA_M2(:,:,1)
!                ENDDO
!                Ptr3D => NULL()
!             END IF
!          ENDDO
!       ENDIF
!       ! Reset diagnostics
!       !!!State_Diag%WetLossLS(:,:,:,:) = 0.0_f4
!    ENDIF
!
!    ! Strat. aerosol surface densities 
!    IF ( ASSOCIATED(State_Diag%AerSurfAreaDust) ) THEN
!       DiagName = 'AerSurfAreaDust'
!       CALL MAPL_GetPointer( Export, Ptr3D, TRIM(DiagName), NotFoundOk=.TRUE., __RC__ )
!       IF ( ASSOCIATED(Ptr3D) ) THEN
!          Ptr3D(:,:,LM:1:-1) = State_Diag%AerSurfAreaDust
!          Ptr3D => NULL()
!       ENDIF
!    ENDIF
!    IF ( ASSOCIATED(State_Diag%AerSurfAreaSLA) ) THEN
!       DiagName = 'AerSurfAreaStratLiquid'
!       CALL MAPL_GetPointer( Export, Ptr3D, TRIM(DiagName), NotFoundOk=.TRUE., __RC__ )
!       ! testing only
!       IF ( ASSOCIATED(Ptr3D) ) THEN
!          Ptr3D(:,:,LM:1:-1) = State_Diag%AerSurfAreaSLA
!          Ptr3D => NULL()
!       ENDIF
!    ENDIF
!    IF ( ASSOCIATED(State_Diag%AerSurfAreaPSC) ) THEN
!       DiagName = 'AerSurfAreaPolarStratCloud'
!       CALL MAPL_GetPointer( Export, Ptr3D, TRIM(DiagName), NotFoundOk=.TRUE., __RC__ )
!       IF ( ASSOCIATED(Ptr3D) ) THEN
!          Ptr3D(:,:,LM:1:-1) = State_Diag%AerSurfAreaPSC
!          Ptr3D => NULL()
!       ENDIF
!    ENDIF
!    IF ( ASSOCIATED(State_Diag%AerNumDenSLA) ) THEN
!       DiagName = 'AerNumDensityStratLiquid'
!       CALL MAPL_GetPointer( Export, Ptr3D, TRIM(DiagName), NotFoundOk=.TRUE., __RC__ )
!       IF ( ASSOCIATED(Ptr3D) ) THEN
!          Ptr3D(:,:,LM:1:-1) = State_Diag%AerNumDenSLA
!          Ptr3D => NULL()
!       ENDIF
!    ENDIF
!    IF ( ASSOCIATED(State_Diag%AerNumDenPSC) ) THEN
!       DiagName = 'AerNumDensityStratParticulate'
!       CALL MAPL_GetPointer( Export, Ptr3D, TRIM(DiagName), NotFoundOk=.TRUE., __RC__ )
!       IF ( ASSOCIATED(Ptr3D) ) THEN
!          Ptr3D(:,:,LM:1:-1) = State_Diag%AerNumDenPSC
!          Ptr3D => NULL()
!       ENDIF
!    ENDIF
!
!    ! PM25
!    DO I = 1,8
!       SELECT CASE ( I )
!          CASE ( 1 )
!             DiagName =  'PM25'
!             PM25ptr  => State_Diag%PM25
!          CASE ( 2 )
!             DiagName =  'PM25ni'
!             PM25ptr  => State_Diag%PM25ni
!          CASE ( 3 )
!             DiagName =  'PM25su'
!             PM25ptr  => State_Diag%PM25su
!          CASE ( 4 )
!             DiagName =  'PM25oc'
!             PM25ptr  => State_Diag%PM25oc
!          CASE ( 5 )
!             DiagName =  'PM25bc'
!             PM25ptr  => State_Diag%PM25bc
!          CASE ( 6 )
!             DiagName =  'PM25ss'
!             PM25ptr  => State_Diag%PM25ss
!          CASE ( 7 )
!             DiagName =  'PM25du'
!             PM25ptr  => State_Diag%PM25du
!          CASE ( 8 )
!             DiagName =  'PM25soa'
!             PM25ptr  => State_Diag%PM25soa
!          CASE DEFAULT 
!             DiagName =  'dummy'
!             PM25ptr  => NULL() 
!       END SELECT
!       IF ( ASSOCIATED ( PM25ptr ) ) THEN   
!          CALL MAPL_GetPointer( Export, Ptr3D, TRIM(DiagName), NotFoundOk=.TRUE., __RC__ )
!          IF ( ASSOCIATED(Ptr3D) ) THEN
!             Ptr3D(:,:,LM:1:-1) = PM25ptr(:,:,:)
!             Ptr3D   => NULL()
!             PM25ptr => NULL()
!          ENDIF
!       ENDIF
!    ENDDO
!
!    ! OH concentration after chemistry
!    IF ( Input_Opt%LCHEM .AND. IsChemTime ) THEN
!       IF ( ASSOCIATED(State_Diag%OHconcAfterChem) ) THEN
!          DiagName = 'OHconcAfterChem'
!          CALL MAPL_GetPointer( Export, Ptr3D, TRIM(DiagName), NotFoundOk=.TRUE., __RC__ )
!          IF ( ASSOCIATED(Ptr3D) ) THEN
!             Ptr3D = State_Diag%OHconcAfterChem(:,:,LM:1:-1)
!          ENDIF
!       ENDIF
!       IF ( ASSOCIATED(State_Diag%O3concAfterChem) ) THEN
!          DiagName = 'O3concAfterChem'
!          CALL MAPL_GetPointer( Export, Ptr3D, TRIM(DiagName), NotFoundOk=.TRUE., __RC__ )
!          IF ( ASSOCIATED(Ptr3D) ) THEN
!             Ptr3D = State_Diag%O3concAfterChem(:,:,LM:1:-1)
!          ENDIF
!       ENDIF
!       IF ( ASSOCIATED(State_Diag%RO2concAfterChem) ) THEN
!          DiagName = 'RO2concAfterChem'
!          CALL MAPL_GetPointer( Export, Ptr3D, TRIM(DiagName), NotFoundOk=.TRUE., __RC__ )
!          IF ( ASSOCIATED(Ptr3D) ) THEN
!             Ptr3D = State_Diag%RO2concAfterChem(:,:,LM:1:-1)
!          ENDIF
!       ENDIF
!    ENDIF
!
!    ! CH4 pseudo flux to balance chemistry
!    DiagName = 'CH4pseudoFlux'
!    CALL MAPL_GetPointer( Export, Ptr2D, TRIM(DiagName), NotFoundOk=.TRUE., __RC__ )
!    IF ( ASSOCIATED(Ptr2D) ) THEN
!       Ptr2D = State_Diag%CH4pseudoFlux
!    ENDIF
!
!    ! AOD 
!    !IF ( ND21 > 0 .AND. Input_Opt%LCHEM .AND. IsChemTime ) THEN
!    IF ( Input_Opt%ND21 > 0 ) THEN
!       ALLOCATE( Diag2D(IM,JM) )
!       DO I = 1,14
!          SELECT CASE ( I )
!             CASE ( 1 )
!                DiagName    = 'AOD550_CLOUD'
!                Diag2D(:,:) = SUM(State_Met%OPTD,DIM=3)
!             CASE ( 2 )
!                DiagName    = 'AOD550_DUST'
!                Diag2D(:,:) = SUM(AD21(:,:,:,4),DIM=3)
!             CASE ( 3 )
!                DiagName    = 'AOD550_SULFATE'
!                Diag2D(:,:) = SUM(AD21(:,:,:,6),DIM=3)
!             CASE ( 4 )
!                DiagName    = 'AOD550_BC'
!                Diag2D(:,:) = SUM(AD21(:,:,:,9),DIM=3)
!             CASE ( 5 )
!                DiagName    = 'AOD550_OC'
!                Diag2D(:,:) = SUM(AD21(:,:,:,12),DIM=3)
!             CASE ( 6 )
!                DiagName    = 'AOD550_SALA'
!                Diag2D(:,:) = SUM(AD21(:,:,:,15),DIM=3)
!             CASE ( 7 )
!                DiagName    = 'AOD550_SALC'
!                Diag2D(:,:) = SUM(AD21(:,:,:,18),DIM=3)
!             CASE ( 8 )
!                DiagName    = 'AOD550_DST1'
!                Diag2D(:,:) = SUM(AD21(:,:,:,21),DIM=3)
!             CASE ( 9 )
!                DiagName    = 'AOD550_DST2'
!                Diag2D(:,:) = SUM(AD21(:,:,:,22),DIM=3)
!             CASE ( 10 )
!                DiagName    = 'AOD550_DST3'
!                Diag2D(:,:) = SUM(AD21(:,:,:,23),DIM=3)
!             CASE ( 11 )
!                DiagName    = 'AOD550_DST4'
!                Diag2D(:,:) = SUM(AD21(:,:,:,24),DIM=3)
!             CASE ( 12 )
!                DiagName    = 'AOD550_DST5'
!                Diag2D(:,:) = SUM(AD21(:,:,:,25),DIM=3)
!             CASE ( 13 )
!                DiagName    = 'AOD550_DST6'
!                Diag2D(:,:) = SUM(AD21(:,:,:,26),DIM=3)
!             CASE ( 14 )
!                DiagName    = 'AOD550_DST7'
!                Diag2D(:,:) = SUM(AD21(:,:,:,27),DIM=3)
!             CASE DEFAULT
!                DiagName    = 'AOD550_DUMMY'
!          END SELECT
!          CALL MAPL_GetPointer( Export, Ptr2D, TRIM(DiagName), NotFoundOk=.TRUE., __RC__ )
!          IF ( ASSOCIATED(Ptr2D) ) THEN
!             Ptr2D(:,:) = Diag2D
!             Ptr2D => NULL()
!          END IF
!       ENDDO
!       DEALLOCATE( Diag2D )
!    ENDIF
!
!    ! Some met variables 
!    DiagName = 'GCC_AIRNUMDEN'
!    CALL MAPL_GetPointer( Export, Ptr3D, TRIM(DiagName), NotFoundOk=.TRUE., __RC__ )
!    IF ( ASSOCIATED(Ptr3D) ) THEN
!       ! AIRNUMDEN is in molec cm-3, convert to mol cm-3
!       Ptr3D = State_Met%AIRNUMDEN(:,:,LM:1:-1)/6.022e23
!    ENDIF
!
!    DiagName = 'GCC_AIRVOL'
!    CALL MAPL_GetPointer( Export, Ptr3D, TRIM(DiagName), NotFoundOk=.TRUE., __RC__ )
!    IF ( ASSOCIATED(Ptr3D) ) THEN
!       ! AIRVOL is [cm3], convert to km3
!       Ptr3D = State_Met%AIRVOL(:,:,LM:1:-1)/1.0e15
!    ENDIF
!
!    ! Reaction rates
!    IF ( Input_Opt%LCHEM .AND. IsChemTime ) THEN
!    IF ( Input_Opt%NN_RxnRates > 0 ) THEN
!       DO I = 1, Input_Opt%NN_RxnRates
!          WRITE(III,'(i3.3)') Input_Opt%RxnRates_IDs(I)
!          DiagName = 'GCC_RR_'//TRIM(III)
!          CALL MAPL_GetPointer( Export, Ptr3D, TRIM(DiagName), NotFoundOK=.TRUE., __RC__ )
!          IF ( ASSOCIATED(Ptr3D) ) THEN
!             Ptr3D = State_Diag%RxnRates(:,:,LM:1:-1,I)
!          ENDIF 
!       ENDDO
!    ENDIF
!    ENDIF
!
!    ! Reaction rate constants
!    IF ( Input_Opt%LCHEM .AND. IsChemTime ) THEN
!    IF ( Input_Opt%NN_RxnRconst > 0 ) THEN
!       DO I = 1, Input_Opt%NN_RxnRconst
!          WRITE(III,'(i3.3)') Input_Opt%RxnRconst_IDs(I)
!          DiagName = 'GCC_RC_'//TRIM(III)
!          CALL MAPL_GetPointer( Export, Ptr3D, TRIM(DiagName), NotFoundOK=.TRUE., __RC__ )
!          IF ( ASSOCIATED(Ptr3D) ) THEN
!             Ptr3D = State_Diag%RxnRconst(:,:,LM:1:-1,I)
!     
!          ENDIF 
!       ENDDO
!    ENDIF
!    ENDIF
!
!    ! J-values 
!    IF ( Input_Opt%LCHEM .AND. IsChemTime ) THEN
!    IF ( Input_Opt%NN_Jvals > 0 ) THEN
!       DO I = 1, Input_Opt%NN_Jvals
!          WRITE(III,'(i3.3)') Input_Opt%Jval_IDs(I)
!          DiagName = 'GCC_JVAL_'//TRIM(III)
!          CALL MAPL_GetPointer( Export, Ptr3D, TRIM(DiagName), NotFoundOk=.TRUE., __RC__ )
!          IF ( ASSOCIATED(Ptr3D) ) THEN
!             Ptr3D = State_Diag%JvalIndiv(:,:,LM:1:-1,I)
!          ENDIF 
!       ENDDO
!    ENDIF
!    ENDIF
!
!  END SUBROUTINE StateDiag2Export
!!EOC
END MODULE GIGC_GEOS5Exports_Mod
