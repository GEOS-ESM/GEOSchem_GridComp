#include "MAPL_Generic.h"
!-------------------------------------------------------------------------
!NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1, GEOS/DAS!
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  GmiChem_GCCMod --- GMI Chemistry Grid Component Class
!
! Grid Component class for the GMI combined stratopshere/troposphere 
! chemistry
!
! !INTERFACE:
!

   MODULE  GmiChem_GCCMod

! !USES:

   USE ESMF
   USE MAPL
   USE Chem_Mod 	     ! Chemistry Base Class
   USE Chem_UtilMod

   USE Species_BundleMod

   USE GmiChemistryMethod_mod,        ONLY : t_Chemistry
   USE GmiDepositionMethod_mod,       ONLY : t_Deposition
   USE GmiSpcConcentrationMethod_mod, ONLY : t_SpeciesConcentration
   USE GmiGrid_mod,                   ONLY : t_gmiGrid
   USE GmiTimeControl_mod,            ONLY : t_GmiClock
   USE GmiESMFrcFileReading_mod,      ONLY : rcEsmfReadTable
   use GmiArrayBundlePointer_mod,     ONLY : t_GmiArrayBundle, CleanArrayPointer
   use GmiFieldBundleESMF_mod,        ONLY : obtainTracerFromBundle
   use GmiFieldBundleESMF_mod,        ONLY : addTracerToBundle
   use GmiStateFieldESMF_mod,         ONLY : getDataFromStateField
   use GmiSwapSpeciesBundlesMod,      ONLY : SwapSpeciesBundles, speciesReg_for_CCM

   IMPLICIT NONE
   INTEGER, PARAMETER :: DBL = KIND(0.00D+00)

#include "setkin_par.h"
#include "GmiParameters.h"
#include "gmi_phys_constants.h"
#include "setkin_mw.h"
#include "setkin_lchem.h"

! !TYPES:

   PRIVATE
   PUBLIC  GmiChemistry_GridComp       ! The GMI object 

! !PUBLIC MEMBER FUNCTIONS:

   PUBLIC  GmiChemistry_GridCompInitialize
   PUBLIC  GmiChemistry_GridCompRun
   PUBLIC  GmiChemistry_GridCompFinalize

!
! !DESCRIPTION:
!
!  This module implements the GMI combined stratopshere/troposphere
!  chemistry
!
! !REVISION HISTORY:
!
!  16Sep2003 da Silva  First crack.
!  24Jan2O05 Nielsen   Implementation of Code 916 chemistry
!  19Dec2005 d Silva   Minor portability mods.
!  30Oct2007 Nielsen   GMI Combo set up
!  01Jul2010 Kouatchou Creation of Chemistry Grid Component class
!  01Jun2015 Nielsen   ExtData replaces Chem_UtilMPread
!
!EOP
!-------------------------------------------------------------------------

  INTEGER, PARAMETER :: RXN_NAME_LENGTH     = 16
  INTEGER, PARAMETER :: RXN_LONGNAME_LENGTH = 256
  INTEGER, PARAMETER :: TOKEN_LENGTH  = 256

  TYPE GmiChemistry_GridComp
   CHARACTER(LEN=255) :: name = "GMI Stratospheric/Tropospheric Chemistry"

! Does the gmichem_import_rst file exist?  If not, some import states will be
! unfilled, and the species will need to "freewheel" through the first time step.
! -------------------------------------------------------------------------------
   LOGICAL :: gotImportRst

! Various switches
! ----------------
   LOGICAL :: pr_diag
   LOGICAL :: do_synoz
   LOGICAL :: do_semiss_inchem
   LOGICAL :: pr_surf_emiss
   LOGICAL :: pr_emiss_3d
   LOGICAL :: do_qqjk_inchem
   LOGICAL :: do_qqjk_reset
   LOGICAL :: pr_qqjk
   LOGICAL :: do_grav_set
   REAL*8  :: pr_nc_period
   LOGICAL :: rd_restart
   LOGICAL :: do_ftiming
   LOGICAL :: pr_smv2
   LOGICAL :: pr_const
   LOGICAL :: pr_qj_o3_o1d
   LOGICAL :: pr_qj_opt_depth
   CHARACTER (LEN=MAX_LENGTH_MET_NAME) :: metdata_name_org
   CHARACTER (LEN=MAX_LENGTH_MET_NAME) :: metdata_name_model

! Dimensions
! ----------
   INTEGER :: i1, i2, im, j1, j2, jm, km

! Surface area of grid cells
! --------------------------
   REAL(KIND=DBL), POINTER :: cellArea(:,:)

! for computing tropospheric OX loss
! ----------------------------------
   INTEGER :: stOX_rxn_count
   CHARACTER(LEN=RXN_NAME_LENGTH),     pointer :: rname(:)  ! vector of reaction short names
   REAL(KIND=DBL), allocatable                 :: rmult(:)  ! vector of multipliers
   CHARACTER(LEN=RXN_LONGNAME_LENGTH), pointer :: rdesc(:)  ! vector of reaction long  names

! Extra diagnostics
! -----------------
   LOGICAL :: verbose

! Map GMI species indices to CCM indices
! --------------------------------------
   INTEGER, POINTER :: mapSpecies(:)

! Component derived type declarations
! -----------------------------------
   TYPE(t_Chemistry )		:: Chemistry
   TYPE(t_Deposition)		:: Deposition
   TYPE(t_gmiGrid   )		:: gmiGrid
   TYPE(t_GmiClock  )           :: gmiClock
   TYPE(t_SpeciesConcentration) :: SpeciesConcentration
 
  END TYPE GmiChemistry_GridComp

CONTAINS


!-------------------------------------------------------------------------
!NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1, GEOS/DAS!
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  GmiChemistry_GridCompInitialize --- Initialize GmiChemistry_GridComp
!
! !INTERFACE:
!

   SUBROUTINE GmiChemistry_GridCompInitialize( self, bgg, bxx, impChem, expChem, nymd, nhms, &
                                      tdt, rc )

   USE GmiSpcConcentrationMethod_mod, ONLY : InitializeSpcConcentration
   USE GmiDepositionMethod_mod,       ONLY : InitializeDeposition
   USE GmiChemistryMethod_mod,        ONLY : InitializeChemistry, initReadChemistry
   USE GmiGrid_mod,		      ONLY : InitializeGmiGrid
   USE GmiTimeControl_mod,            ONLY : Set_curGmiDate, Set_curGmiTime
   USE GmiTimeControl_mod,            ONLY : Set_begGmiDate, Set_begGmiTime
   USE GmiTimeControl_mod,            ONLY : Set_numTimeSteps

   IMPLICIT none
   INTEGER, PARAMETER :: DBL = KIND(0.00D+00)

! !INPUT PARAMETERS:

   TYPE(Species_Bundle), INTENT(in) :: bgg                ! GMI Species - transported
   TYPE(Species_Bundle), INTENT(in) :: bxx                ! GMI Species - not transported
   INTEGER,              INTENT(IN) :: nymd, nhms         ! Time from AGCM
   REAL,                 INTENT(IN) :: tdt                ! Chemistry time step (secs)

! !OUTPUT PARAMETERS:

   TYPE(GmiChemistry_GridComp), INTENT(INOUT)  :: self      ! Grid Component
   TYPE(ESMF_State),   INTENT(INOUT)  :: impChem    ! Import State
   TYPE(ESMF_State),   INTENT(INOUT)  :: expChem    ! Export State

   INTEGER, INTENT(out) ::  rc                  ! Error return code:
                                                !  0 - all is well
                                                !  1 - 

! !DESCRIPTION: Initializes the GMI Grid Component. It primarily sets
!               the import state.
!
! !REVISION HISTORY:
!
!  18Sep2003 da Silva  First crack.
!  30Jun2007 Nielsen   GMI Combo set up
!
!EOP
!-------------------------------------------------------------------------

   INTEGER, PARAMETER :: RC_DATA_LINE    = 1
   INTEGER, PARAMETER :: RC_END_OF_TABLE = 2
   INTEGER, PARAMETER :: RC_END_OF_FILE  = 3

   CHARACTER(LEN=*), PARAMETER :: IAm    = 'GmiChem_GridCompClassInitialize'
   CHARACTER(LEN=255) :: rcfile = 'GMI_GridComp.rc'
   CHARACTER(LEN=255) :: namelistFile
   CHARACTER(LEN=255) :: importRestartFile
   CHARACTER(LEN=255) :: string
   
   type (ESMF_Config) :: gmiConfig

   INTEGER :: ios, m, n, STATUS
   INTEGER :: i, i1, i2, ic, im, j, j1, j2, jm, k, km, kReverse

   INTEGER :: i1_gl, i2_gl, ju1_gl, j2_gl 
   INTEGER :: ju1, jv1, jv1_gl, j1p, j2p
   INTEGER :: k1, k2, k1_gl, k2_gl
   INTEGER :: ilong, ilat, ivert, itloop
   INTEGER :: NPIJ, NPI, NPJ
   INTEGER :: ilo, ihi, julo, jvlo, jhi
   INTEGER :: ilo_gl, ihi_gl, julo_gl, jvlo_gl, jhi_gl
   INTEGER :: gmi_nborder
   INTEGER :: NMR      ! number of species from the GMI_Mech_Registry.rc
   INTEGER :: LogicalUnitNum

   INTEGER :: loc_proc, locGlobProc, commu_slaves
   LOGICAL :: one_proc, rootProc
   LOGICAL :: exists,open,found
   
   REAL :: qmin, qmax, tokgCPerBox

   TYPE(ESMF_FieldBundle) :: qqkBundle
   TYPE(ESMF_FieldBundle) :: qqjBundle
   REAL(rPrec), POINTER   :: var(:,:,:)
   INTEGER                :: numVars, ib
   CHARACTER(LEN=4)       :: binName
   CHARACTER(LEN=ESMF_MAXSTR) :: varName

   CHARACTER(LEN=TOKEN_LENGTH) ::  str_arr(3)   ! strings for rxn_name, multiplier and rxn_longname
   INTEGER                     :: nrxn, nx, item_count, retcode
   CHARACTER(LEN=32)           :: table_name

! Grid cell area can be set by initialize
! ---------------------------------------
   REAL, POINTER, DIMENSION(:,:) :: cellArea

   self%name = 'GMI combined Stratosphere/Troposphere Chemistry'

!  Initialize local variables
!  --------------------------
   rc = 0
   i1 = self%i1
   i2 = self%i2
   im = self%im
   
   j1 = self%j1
   j2 = self%j2
   jm = self%jm
   
   km = self%km

     !-------------------
     ! Load resource file
     !-------------------

      IF ( MAPL_AM_I_ROOT() ) THEN
         PRINT *," "
         PRINT *,TRIM(IAm)//":"
         PRINT *,"Starting Reading the GMI Resource File for Chemistry"
      ENDIF

      gmiConfig = ESMF_ConfigCreate( __RC__ )

      call ESMF_ConfigLoadFile(gmiConfig, TRIM(rcfile), __RC__ )

!     stOX Reactions
!     --------------
      table_name = 'stOX_loss_reactions::'

      call ESMF_ConfigGetDim(gmiConfig, nrxn, nx, label=TRIM(table_name), rc=rc)
      _ASSERT(rc==0, TRIM(Iam)//': Cannot get dims for table '//TRIM(table_name)//' in '//TRIM(rcfile))

      call ESMF_ConfigFindLabel(gmiConfig, TRIM(table_name), rc=rc)
      _ASSERT(rc==0, TRIM(Iam)//': Cannot find '//TRIM(table_name)//' in file '//TRIM(rcfile))

!     Allocate memory
!     ---------------
      self%stOX_rxn_count = nrxn
      allocate ( self%rname(nrxn), self%rmult(nrxn), self%rdesc(nrxn), __STAT__ )

!     Read the reactions
!     ------------------
      do i=1,nrxn
         call get_line ( gmiConfig, 3, str_arr, item_count, retcode )
         select case( retcode )
           case( RC_END_OF_FILE  )
             _FAIL(TRIM(Iam)//': early EOF in file '//TRIM(rcfile))
           case( RC_END_OF_TABLE )
             _FAIL(TRIM(Iam)//': table too short '//TRIM(table_name)//' in file '//TRIM(rcfile))
           case( RC_DATA_LINE    )
             _ASSERT(item_count==3, TRIM(Iam)//': fewer than 3 entries in '//TRIM(table_name)//' in file '//TRIM(rcfile))
             self%rname(i) = str_arr(1)
             READ(           str_arr(2),*,IOSTAT=rc) self%rmult(i)
             _ASSERT(rc==0, TRIM(Iam)//': Bad scaling factor for rxn '//TRIM(str_arr(1))//' in file '//TRIM(rcfile)//': '//TRIM(str_arr(2)))
             self%rdesc(i) = str_arr(3)
         end select
      end do

      IF( MAPL_AM_I_ROOT() ) THEN
        PRINT*,'stOX reactions>>>'
        do i=1,self%stOX_rxn_count
          PRINT*,i,TRIM(self%rname(i)), self%rmult(i), TRIM(self%rdesc(i))
        end do
        PRINT*,'stOX reactions<<<'
      END IF



      call ESMF_ConfigGetAttribute(gmiConfig, importRestartFile, &
     &                label   = "importRestartFile:",            &
     &                default = ' ', rc=STATUS )
      VERIFY_(STATUS)

      !------------------------------
      ! Advection related variables
      !------------------------------

      call ESMF_ConfigGetAttribute(gmiConfig, value=self%do_grav_set, &
     &           label="do_grav_set:", default=.false., rc=STATUS)
      VERIFY_(STATUS)

      !------------------------------
      ! Emission related variables
      !------------------------------

      call ESMF_ConfigGetAttribute(gmiConfig, value=self%do_synoz, &
     &           label="do_synoz:", default=.false., rc=STATUS)
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(gmiConfig, value=self%do_semiss_inchem, &
     &           label="do_semiss_inchem:", default=.false., rc=STATUS)
      VERIFY_(STATUS)

      !------------------------------
      ! Diagnostics related variables
      !------------------------------

      call ESMF_ConfigGetAttribute(gmiConfig, value=self%pr_diag, &
     &           label="pr_diag:", default=.false., rc=STATUS)
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(gmiConfig, value=self%verbose, &
     &           label="verbose:", default=.false., rc=STATUS)
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(gmiConfig, value=self%pr_surf_emiss, &
     &           label="pr_surf_emiss:", default=.false., rc=STATUS)
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(gmiConfig, value=self%pr_emiss_3d, &
     &           label="pr_emiss_3d:", default=.false., rc=STATUS)
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(gmiConfig, value=self%pr_qqjk, &
     &           label="pr_qqjk:", default=.false., rc=STATUS)
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(gmiConfig, value=self%do_qqjk_reset, &
     &           label="do_qqjk_reset:", default=.true., rc=STATUS)
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(gmiConfig, self%pr_nc_period, &
     &                label   = "pr_nc_period:", &
     &                default = -1.0d0, rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(gmiConfig, value=self%do_ftiming, &
     &               label="do_ftiming:", default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(gmiConfig, value=self%rd_restart, &
     &               label="rd_restart:", default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(gmiConfig, value=self%pr_qj_o3_o1d, &
     &               label="pr_qj_o3_o1d:", default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(gmiConfig, value=self%pr_qj_opt_depth, &
     &               label="pr_qj_opt_depth:", default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(gmiConfig, value=self%pr_smv2, &
     &               label="pr_smv2:", default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(gmiConfig, value=self%pr_const, &
     &               label="pr_const:", default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(gmiConfig, self%metdata_name_org, &
     &                label   = "metdata_name_org:", &
     &                default = 'GMAO', rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(gmiConfig, self%metdata_name_model, &
     &                label   = "metdata_name_model:", &
     &                default = 'GEOS-5', rc=STATUS )


      !----------------------------
      ! Chemistry Related Variables
      !----------------------------

      call ESMF_ConfigGetAttribute(gmiConfig, value=self%do_qqjk_inchem, &
     &           label="do_qqjk_inchem:", default=.false., rc=STATUS)
      VERIFY_(STATUS)

! Does the GMICHEM import restart file exist?  If not,
! the species must "freewheel" through the first time step.
! ---------------------------------------------------------
   INQUIRE(FILE=TRIM(importRestartFile),EXIST=self%gotImportRst)
   IF( MAPL_AM_I_ROOT() ) THEN
    PRINT *," ",TRIM(importRestartFile)," exists: ",self%gotImportRst
    PRINT *," "
   END IF

!  GMI grid specification
!  ----------------------
   gmi_nborder = 0
   i1_gl = 1
   i2_gl = i2
   ju1_gl = 1
   jv1_gl = 1
   j2_gl = j2
   ju1 = j1
   jv1 =j1
   k1 = 1
   k2 = km
   k1_gl = 1
   k2_gl = km
   NPIJ = 16     ! These three integers are irrelevant
   NPI = 4
   NPJ = 4
   ilo = i1 - gmi_nborder
   ihi = i2 + gmi_nborder 
   julo = ju1 - gmi_nborder
   jvlo = jv1 - gmi_nborder
   jhi = j2 + gmi_nborder
   ilo_gl = i1_gl  - gmi_nborder
   ihi_gl = i2_gl  + gmi_nborder
   julo_gl = ju1_gl - gmi_nborder
   jvlo_gl = jv1_gl - gmi_nborder
   jhi_gl = j2_gl  + gmi_nborder
   j1p = 0
   j2p = j2_gl - j1p + 1
   ilong = i2 - i1 + 1
   ilat = j2 - ju1 + 1
   ivert = k2 - k1 + 1
   itloop = ilat * ilong * ivert

   one_proc = .FALSE.
   loc_proc = -99
   locGlobProc = -99
   commu_slaves = -99
   
   rootProc=.FALSE.
   IF( MAPL_AM_I_ROOT() ) THEN
    rootProc=.TRUE.
   END IF

! Set GMI's clock
! ---------------
   CALL Set_begGmiDate(self%gmiClock, nymd)
   CALL Set_begGmiTime(self%gmiClock, nhms)
   CALL Set_curGmiDate(self%gmiClock, nymd)
   CALL Set_curGmiTime(self%gmiClock, nhms)
   CALL Set_numTimeSteps(self%gmiClock, 0)

! Grid box surface area, m^{2}
! ----------------------------
   CALL MAPL_GetPointer(impChem, cellArea, 'AREA', rc=STATUS)
   VERIFY_(STATUS)
   ALLOCATE(self%cellArea(i1:i2,j1:j2), STAT=STATUS)
   VERIFY_(STATUS)
   self%cellArea(i1:i2,j1:j2)=cellArea(i1:i2,j1:j2)

! Discretization
! --------------
   CALL InitializeGmiGrid(self%gmiGrid, NPIJ, NPI, NPJ, &
                          gmi_nborder, i1, i2, ju1, jv1, j2, k1, k2, &
                          i1_gl, i2_gl, ju1_gl, jv1_gl, j2_gl, k1_gl, k2_gl, &
                          ilo, ihi, julo, jvlo, jhi, &
                          ilo_gl, ihi_gl, julo_gl, jvlo_gl, jhi_gl, &
                          ilong, ilat, ivert, itloop, j1p, j2p)  

! Perform a consistency check with setkin_par.h.
!   NSP is the number of species in setkins
!   NMR is the number of species in the Mech Registry
!
!   H2O   is in setkins,           but not in GMI_Mech_Registry
!   AOA   is in GMI_Mech_Registry, but not in setkins
!   T2M15 is in GMI_Mech_Registry, but not in setkins
!
!   The number of species common to both is therefore
!     NMR - 2     ! number in GMI_Mech_Registry - 2
!     NSP - 1     ! number in setkins - 1
! --------------------------------------------------------------------------------
   NMR = bgg%nq + bxx%nq
   IF( NMR-2 /= NSP-1 ) THEN
    PRINT *,TRIM(IAm),': Number of species from GMI_Mech_Registry.rc does not match number in setkin_par.h'
    STATUS = 1
    VERIFY_(STATUS)
   END IF

! Photolysis reaction list.  Read from kinetics
! text file after finding an available logical unit.
! --------------------------------------------------
   found = .FALSE.
   i = 11
   DO WHILE (.NOT. found .AND. i <= 99)
    INQUIRE(UNIT=i,EXIST=exists,OPENED=open)
    IF(exists .AND. .NOT. open) THEN
     found = .TRUE.
     LogicalUnitNum = i
    END IF
    i = i+1
   END DO

! Allocate space, etc., but the initialization of the
! species from the internal state is left to the run method.
! ----------------------------------------------------------

      CALL InitializeSpcConcentration(self%SpeciesConcentration,              &
                     self%gmiGrid, gmiConfig, NSP, NMF, NCHEM,                &
                     loc_proc)

      CALL InitializeChemistry(self%Chemistry, self%gmiGrid,                  &
                     gmiConfig, loc_proc, NSP, self%pr_diag,                  &
                     self%pr_qqjk, self%do_qqjk_inchem, self%pr_smv2,         &
                     rootProc, tdt)

  CALL ESMF_ConfigDestroy(gmiConfig, __RC__ )
  IF(self%pr_qqjk .AND. .NOT. self%do_qqjk_inchem) THEN
   IF(MAPL_AM_I_ROOT()) THEN
    PRINT *,TRIM(IAm)//": Initializing reaction rate bundles"
    PRINT *," "
   END IF

! Initialize the qqjBundle
! ------------------------
   CALL ESMF_StateGet(expChem, 'gmiQQJ' , qqjBundle,   RC=STATUS)
   VERIFY_(STATUS)

   ! Allocate space to the bundle
   DO ib = 1, NUM_J
      ALLOCATE( var(i1:i2, j1:j2, 1:km), STAT=STATUS)
      VERIFY_(STATUS)
      var(:,:,:)  = 0.0d0

      WRITE (binName ,'(i4.4)') ib
      varName = 'qqj'//binName

      CALL addTracerToBundle (qqjBundle, var, bgg%grid_esmf, varName)
   END DO

   ! Sanity check

   CALL ESMF_FieldBundleGet(qqjBundle, fieldCount=numVars , rc=STATUS)
   VERIFY_(STATUS)
   _ASSERT(NUM_J == numVars,'needs informative message')

! Initialize the qqkBundle
! ------------------------
   CALL ESMF_StateGet(expChem, 'gmiQQK' , qqkBundle,   RC=STATUS)
   VERIFY_(STATUS)

   ! Allocate space to the bundle
   DO ib = 1, NUM_K
      ALLOCATE( var(i1:i2, j1:j2, 1:km), STAT=STATUS)
      VERIFY_(STATUS)
      var(:,:,:)  = 0.0d0

      WRITE (binName ,'(i4.4)') ib
      varName = 'qqk'//binName

      CALL addTracerToBundle (qqkBundle, var, bgg%grid_esmf, varName)
   END DO

   ! Sanity check

   CALL ESMF_FieldBundleGet(qqkBundle, fieldCount=numVars , rc=STATUS)
   VERIFY_(STATUS)
   _ASSERT(NUM_K == numVars,'needs informative message')

  END IF ! pr_qqjk is .TRUE.

    !---------------------------------------------------------------
    ! Create and populate the array that maps GMI species indices to
    ! GEOS-5 species indices
    !---------------------------------------------------------------

    allocate(self%mapSpecies(NSP))
    self%mapSpecies(:) = speciesReg_for_CCM(lchemvar, NSP, bgg%reg%vname, bxx%reg%vname )

  RETURN

   CONTAINS

!     -------------------
!     GET_LINE
!     Advance one line and then try to read <expected_entries> items
!     Retcode will be set to one of these values:
!       RC_END_OF_FILE    - cannot advance a line
!       RC_END_OF_TABLE   - first item is '::'
!       RC_DATA_LINE      - at least one entry has been put into str_arr
!     Note that ESMF automatically skips over blank lines and comment lines
!     -------------------
      subroutine get_line ( cf, expected_entries, str_arr, item_count, retcode )
!     -------------------
      type(ESMF_Config),           intent(inout)  :: cf
      integer,                     intent(in)     :: expected_entries  ! read this many items
      character(len=TOKEN_LENGTH), intent(inout)  :: str_arr(*)   ! space for one or more items
      integer,                     intent(out)    :: item_count   ! how many were successfully read
      integer,                     intent(out)    :: retcode      ! see possible values above

      integer :: i

      call ESMF_ConfigNextLine(cf, rc=rc)
      if ( rc/=0 ) then
        retcode = RC_END_OF_FILE
        return
      end if

      ! Because ESMF skips over blank lines, rc should always be 0:
      call ESMF_ConfigGetAttribute(cf, str_arr(1), rc=rc)
      if ( rc/=0 ) then
        retcode = RC_END_OF_FILE
        return
      end if
      if ( INDEX(str_arr(1), '::' ) == 1 ) then
        retcode = RC_END_OF_TABLE
        return
      end if

      retcode = RC_DATA_LINE
      item_count = 1
      do i=2,expected_entries
        call ESMF_ConfigGetAttribute(cf, str_arr(i), rc=rc)
        if (rc==0) item_count = item_count + 1
      end do

      return

      end subroutine get_line


  END SUBROUTINE GmiChemistry_GridCompInitialize

!-------------------------------------------------------------------------
!NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1, GEOS/DAS!
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  GmiChemistry_GridCompRun --- The GMI Driver
!
! !INTERFACE:
!

   SUBROUTINE GmiChemistry_GridCompRun ( self, bgg, bxx, impChem, expChem, nymd, nhms, &
                                tdt, rc )

! !USES:

   USE GmiTimeControl_mod,            ONLY : Set_curGmiDate, Set_curGmiTime
   USE GmiTimeControl_mod,            ONLY : Set_numTimeSteps, Get_numTimeSteps
   USE GmiTimeControl_mod,            ONLY : Set_gmiSeconds, GetSecondsFromJanuary1
   USE GmiSpcConcentrationMethod_mod, ONLY : resetFixedConcentration
   USE GmiSolar_mod,                  ONLY : CalcCosSolarZenithAngle
   USE GmiChemistryMethod_mod,        ONLY : RunChemistry

   IMPLICIT none

! !INPUT/OUTPUT PARAMETERS:

   TYPE(GmiChemistry_GridComp), INTENT(INOUT) :: self   ! Grid Component
   TYPE(Species_Bundle),        INTENT(INOUT) :: bgg    ! GMI Species - transported
   TYPE(Species_Bundle),        INTENT(INOUT) :: bxx    ! GMI Species - not transported

! !INPUT PARAMETERS:

   TYPE(ESMF_State), INTENT(INOUT) :: impChem ! Import State
   INTEGER, INTENT(IN) :: nymd, nhms          ! time
   REAL,    INTENT(IN) :: tdt                 ! chemical timestep (secs)

! !OUTPUT PARAMETERS:

   TYPE(ESMF_State), INTENT(INOUT) :: expChem   ! Export State
   INTEGER, INTENT(OUT) ::  rc                  ! Error return code:
                                                !  0 - all is well
                                                !  1 -

! !DESCRIPTION: This routine implements the GMI Strat/Trop Driver. That 
!               is, adds chemical tendencies to each of the constituents
!
! !IMPLEMENTATION NOTES:
!
!  No pointer is reservered in the export state for deposition of water.
!
! !REVISION HISTORY:
!
!  18Sep2003 da Silva  First crack.
!  24Jan2005 Nielsen   Implementation of Code 916 chemistry
!  30Oct2007 Nielsen   Implementation of GMI cmbined 
!                       stratosphere/troposphere chemistry
!
!EOP
!-------------------------------------------------------------------------

   CHARACTER(LEN=*), PARAMETER :: IAm    = 'GmiChemistry_GridCompRun'

!  Imports from GEOS-5
!  -------------------
   REAL, POINTER, DIMENSION(:,:) ::  zpbl

   REAL, POINTER, DIMENSION(:,:,:) :: airdens, ple, Q, T, zle

!  Exports not part of internal state
!  ----------------------------------
   REAL, POINTER, DIMENSION(:,:,:) :: O3ppmv, O3, stOX_loss

!  Exports for reactions diagnostics
!  ---------------------------------
#include "Reactions_DeclarePointer___.h"

! SAD related variables coming from the SAD module
! ------------------------------------------------
   REAL(rPrec), POINTER, DIMENSION(:,:,:) :: HNO3GASsad

!  Local
!  -----
   INTEGER :: cymd, dymd, emiss_opt, hms
   INTEGER :: i, i1, i2, ic, im, iXj
   INTEGER :: j, j1, j2, jm
   INTEGER :: k, km, kReverse
   INTEGER :: loc_proc
   INTEGER :: n, STATUS

   INTEGER, PARAMETER :: ToGMI = 1
   INTEGER, PARAMETER :: FromGMI = -1

   REAL :: pi,degToRad,radToDeg,OneOverDt

   REAL, PARAMETER :: mwtAir = 28.9
   REAL, PARAMETER :: rStar = 8.314E+03
   REAL, PARAMETER :: Pa2hPa = 0.01
   REAL, PARAMETER :: ToGrPerKg = 1000.00
   REAL, PARAMETER :: secPerDay = 86400.00
   REAL, PARAMETER :: err = 1.00E-04

   REAL(KIND=DBL) :: chemDt, dayOfYear

   CHARACTER(LEN=255) :: speciesName
   CHARACTER(LEN=255) :: importName

   LOGICAL :: found, rootProc
   LOGICAL, PARAMETER :: doThis = .FALSE.
   
! Allocatables.  Use KIND=DBL where GMI expects REAL*8.
! -----------------------------------------------------
   REAL, ALLOCATABLE :: var3d(:,:,:)
   REAL, ALLOCATABLE :: pl(:,:,:)

   REAL(KIND=DBL), ALLOCATABLE :: tropopausePress(:,:)
   REAL(KIND=DBL), ALLOCATABLE :: pctm2(:,:)

   REAL(KIND=DBL), ALLOCATABLE :: mass(:,:,:)
   REAL(KIND=DBL), ALLOCATABLE :: press3c(:,:,:)
   REAL(KIND=DBL), ALLOCATABLE :: press3e(:,:,:)
   REAL(KIND=DBL), ALLOCATABLE :: gridBoxThickness(:,:,:)
   REAL(KIND=DBL), ALLOCATABLE :: kel(:,:,:)
   REAL(KIND=DBL), ALLOCATABLE :: humidity(:,:,:)

   REAL(KIND=DBL), ALLOCATABLE :: HNO3GAS(:,:,:)
   REAL(KIND=DBL), ALLOCATABLE :: surfEmissForChem(:,:,:)

   TYPE (t_GmiArrayBundle), POINTER :: gmiQJ(:) => null()
   TYPE (t_GmiArrayBundle), POINTER :: gmiQK(:) => null()

   TYPE (t_GmiArrayBundle), POINTER :: gmiQQJ(:) => null()
   TYPE (t_GmiArrayBundle), POINTER :: gmiQQK(:) => null()

   loc_proc = -99

!  Grid specs from Chem_Bundle%grid
!  --------------------------------
   rc = 0
   i1 = bgg%grid%i1
   i2 = bgg%grid%i2
   im = bgg%grid%im
   
   j1 = bgg%grid%j1
   j2 = bgg%grid%j2
   jm = bgg%grid%jm
   
   km = bgg%grid%km
   
   iXj = (i2-i1+1)*(j2-j1+1)

!  Some real constants
!  -------------------
   pi = 4.00*ATAN(1.00)
   degToRad = pi/180.00
   radToDeg = 180.00/pi
   chemDt = tdt
   OneOverDt = 1.00/tdt

   rootProc=.FALSE.
   IF( MAPL_AM_I_ROOT() ) THEN
    rootProc=.TRUE.
   END IF

!  We need lots of pointers!
!  -------------------------
   CALL FindPointers(STATUS)
   VERIFY_(STATUS)

!  Reserve some local work space
!  -----------------------------
   ALLOCATE(    tropopausePress(i1:i2,j1:j2),STAT=STATUS)
   ALLOCATE(              pctm2(i1:i2,j1:j2),STAT=STATUS)

   ALLOCATE(                pl(i1:i2,j1:j2,1:km),STAT=STATUS)
   VERIFY_(STATUS)
   ALLOCATE(             var3d(i1:i2,j1:j2,1:km),STAT=STATUS)
   VERIFY_(STATUS)
   ALLOCATE(              mass(i1:i2,j1:j2,1:km),STAT=STATUS)
   VERIFY_(STATUS)
   ALLOCATE(           press3c(i1:i2,j1:j2,1:km),STAT=STATUS)
   VERIFY_(STATUS)
   ALLOCATE(           press3e(i1:i2,j1:j2,0:km),STAT=STATUS)
   VERIFY_(STATUS)
   ALLOCATE(  gridBoxThickness(i1:i2,j1:j2,1:km),STAT=STATUS)
   VERIFY_(STATUS)
   ALLOCATE(               kel(i1:i2,j1:j2,1:km),STAT=STATUS)
   VERIFY_(STATUS)
   ALLOCATE(          humidity(i1:i2,j1:j2,1:km),STAT=STATUS)
   VERIFY_(STATUS)

   ALLOCATE(        HNO3GAS(i1:i2,j1:j2,1:km),STAT=STATUS)
   VERIFY_(STATUS)

   !---------------------------------
   ! Obtain data from the ESMF Bundle
   !---------------------------------

   call obtainBundles(expChem)

!  Layer mean pressures. NOTE: ple(:,:,0:km)
!  -----------------------------------------
   DO k=1,km
    pl(i1:i2,j1:j2,k)=(ple(i1:i2,j1:j2,k-1)+ple(i1:i2,j1:j2,k))*0.50
   END DO
   
! Set GMI's clock
! ---------------
   CALL Set_curGmiDate(self%gmiClock, nymd)
   CALL Set_curGmiTime(self%gmiClock, nhms)
   CALL Get_numTimeSteps(self%gmiClock, ic)
   CALL Set_numTimeSteps(self%gmiClock, ic+1)
   CALL Set_gmiSeconds(self%gmiClock, (ic+1)*chemDt)

! Update the following time-dependent boundary conditions:
!  Fixed concentration species
!  Stratospheric sulfate surface area
!  Emissions
! --------------------------------------------------------

   CALL Acquire_Clims(STATUS)
   VERIFY_(STATUS)


! -----------------------------------------
! Pass the following tests before proceding 
! -----------------------------------------

! Chemistry option
! ----------------
   IF(self%Chemistry%chem_opt /= 2) THEN
    IF( MAPL_AM_I_ROOT() ) THEN
     PRINT *,TRIM(IAm),': Code not ready for chem_opt =',self%Chemistry%chem_opt
    END IF
    rc = 61
    RETURN
   END IF


! Diagnostics capabilities enabled?
! ---------------------------------
   IF(self%pr_qqjk .AND. self%do_qqjk_inchem) THEN
    IF( MAPL_AM_I_ROOT() ) THEN
     PRINT *,TRIM(IAm),': Code not ready for pr_qqjk=',self%pr_qqjk,' and do_qqjk_inchem=',self%do_qqjk_inchem
    END IF
    rc = 61
    RETURN
   END IF

! Grab imports and do units conversions
! -------------------------------------

   CALL SatisfyImports(STATUS)
   VERIFY_(STATUS)

! Hand the species concentrations to GMI's bundle
! -----------------------------------------------
   IF (self%gotImportRst) THEN
      CALL SwapSpeciesBundles(ToGMI, self%SpeciesConcentration%concentration, &
               bgg%qa, bxx%qa, Q, self%mapSpecies, lchemvar, self%do_synoz, NSP, &
               STATUS)
      VERIFY_(STATUS)
   END IF

! Impose fixed concentrations
! ---------------------------

   IF(self%SpeciesConcentration%num_fixed_const > 0) THEN
    CALL resetFixedConcentration(self%SpeciesConcentration, self%gmiClock, self%gmiGrid, NSP)
   END IF

! ------------------------------------------------------------------------
! Chemistry
! ------------------------------------------------------------------------

      IF (self%gotImportRst) &
         CALL RunChemistry(self%Chemistry, self%SpeciesConcentration,          &
                 self%gmiClock, self%gmiGrid, press3c, press3e,                &
                 gridBoxThickness, self%cellArea, mass, kel, humidity, pctm2,  &
                 loc_proc, NSP, self%do_qqjk_reset, HNO3GAS, gmiQK,            &
                 gmiQQK, gmiQJ, gmiQQJ, surfEmissForChem, self%pr_diag,        &
                 self%do_ftiming, self%do_qqjk_inchem, self%pr_qqjk,           &
                 self%do_semiss_inchem, self%pr_smv2, self%pr_nc_period,       &
                 rootProc, self%metdata_name_org,             &
                 self%metdata_name_model, tdt)

! Return species concentrations to the chemistry bundle
! -----------------------------------------------------
   IF (self%gotImportRst) then
      CALL SwapSpeciesBundles(FromGMI, self%SpeciesConcentration%concentration, &
               bgg%qa, bxx%qa, Q, self%mapSpecies, lchemvar, self%do_synoz, NSP,  &
               STATUS)
      VERIFY_(STATUS)
   END IF

! Export states
! -------------
   CALL FillExports(STATUS)
   VERIFY_(STATUS)

! Scratch local work space
! ------------------------
   DEALLOCATE(tropopausePress, pctm2, STAT=STATUS)
   VERIFY_(STATUS)
   DEALLOCATE(pl, mass, press3c, press3e, gridBoxThickness, kel, humidity, &
              var3d, HNO3GAS, STAT=STATUS)
   VERIFY_(STATUS)

   call CleanArrayPointer(gmiQJ, STATUS)
   VERIFY_(STATUS)
   call CleanArrayPointer(gmiQK, STATUS)
   VERIFY_(STATUS)
   IF (self%pr_qqjk) THEN
      call CleanArrayPointer(gmiQQJ, STATUS)
      VERIFY_(STATUS)
      call CleanArrayPointer(gmiQQK, STATUS)
      VERIFY_(STATUS)
   END IF

! IMPORTANT: Reset this switch to .TRUE. after first pass.
! --------------------------------------------------------
   self%gotImportRst = .TRUE.

   RETURN

CONTAINS

!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  Acquire_Clims
!
! !INTERFACE:

  SUBROUTINE Acquire_Clims(rc)
  
  IMPLICIT NONE

  INTEGER, INTENT(OUT) :: rc
!
! !DESCRIPTION:
!
!  Obtain the following climatologies from files supplied by
!  the GMI project:
!   Fixed concentration species
!   Stratospheric sulfate surface area
!   Emissions
!
!EOP
!---------------------------------------------------------------------------

  CHARACTER(LEN=255) :: IAm
  INTEGER :: lightning_opt
  INTEGER, ALLOCATABLE :: flag(:)

  REAL, POINTER, DIMENSION(:,:,:) :: PTR3D
  
  rc = 0
  IAm = "Acquire_Clims"

! Update fixed concentration species
! ----------------------------------
  IF(self%SpeciesConcentration%num_fixed_const > 0) THEN
   DO ic = 1,self%SpeciesConcentration%num_fixed_const
    i = self%SpeciesConcentration%fixed_const_map(ic)

    speciesName = TRIM(lchemvar(i))
    importName = TRIM(speciesName)//'_FIXED'
    CALL MAPL_GetPointer(impChem, PTR3D, TRIM(importName), RC=STATUS)
    VERIFY_(STATUS)
    self%SpeciesConcentration%fixed_const(i1:i2,j1:j2,1:km,ic) = PTR3D(i1:i2,j1:j2,km:1:-1)
    NULLIFY(PTR3D)

   END DO
  END IF

  RETURN
 END SUBROUTINE Acquire_Clims

!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  FillExports
!
! !INTERFACE:

  SUBROUTINE FillExports(rc)
  
  IMPLICIT NONE

  INTEGER, INTENT(OUT) :: rc
!
! !DESCRIPTION:
!
!  Fill the export state
!
!EOP
!---------------------------------------------------------------------------

  INTEGER :: rxn_index
  CHARACTER(LEN=RXN_NAME_LENGTH) :: one_name
  REAL(KIND=DBL), allocatable :: stOX_loss_dbl(:,:,:)

  CHARACTER(LEN=255) :: IAm
  
  rc=0
  IAm = "FillExports"

! ---------------------------------------------------------------------
! To conform with GEOS-5's other chemistry components, GMICHEM's ozone
! states are as follows:
!
!    NAME   STATE    Units             Comments
!  ------ --------   ----------------- --------------------------------
!      OX  Internal  mol/mol           Reqired name for ANALYSIS bundle
!      O3  Export    kg/kg             OX(vmr)*48/28.97
!  03PPMV  Export    ppmv              OX(vmr)*1.00E+06
!
! ---------------------------------------------------------------------

! Note: SwapSpeciesBundles(ToGMI,rc) must already be done.
! --------------------------------------------------------
   DO i=1,bgg%reg%nq
    IF(TRIM(bgg%reg%vname(i)) == "OX") ic = i
   END DO

   IF(ASSOCIATED(O3ppmv)) &
    O3ppmv(i1:i2,j1:j2,1:km) = bgg%qa(ic)%data3d(i1:i2,j1:j2,1:km)*1.00E+06

   IF(ASSOCIATED(O3)) &
        O3(i1:i2,j1:j2,1:km) = bgg%qa(ic)%data3d(i1:i2,j1:j2,1:km)*(MAPL_O3MW/MAPL_AIRMW)

! --------------------------------------------------------------------
! Reaction rate constants (q) and rates (qq)
! --------------------------------------------------------------------

#include "QJ_FillExports___.h"
#include "QK_FillExports___.h"

   IF(ASSOCIATED(stOX_loss)) THEN
     _ASSERT(self%pr_qqjk, TRIM(Iam)//': to compute stOX_loss, set GMI pr_qqjk = TRUE')
   END IF

   IF(self%pr_qqjk) THEN

#include "QQK_FillExports___.h"
#include "QQJ_FillExports___.h"

     IF(ASSOCIATED(stOX_loss)) THEN
       ALLOCATE(stOX_loss_dbl(i1:i2,j1:j2,1:km),__STAT__)

       stOX_loss_dbl(:,:,:) = 0.0
       do i=1,self%stOX_rxn_count

         one_name = self%rname(i)

         ! Assume each name is QQKnnn or QQJnnn
         READ( one_name(4:6),*,IOSTAT=rc) rxn_index
         _ASSERT(rc==0, TRIM(Iam)//': trouble extracting index from '//TRIM(self%rname(i)))

!        IF(MAPL_AM_I_ROOT( ))PRINT *,TRIM(IAm),': Add to stOX_loss from '//TRIM(one_name)//' index ',rxn_index

         SELECT CASE (one_name(1:3))

           CASE("QQK")
             stOX_loss_dbl(i1:i2,j1:j2,1:km) = &
             stOX_loss_dbl(i1:i2,j1:j2,1:km) + self%rmult(i) * gmiQQK(rxn_index)%pArray3D(i1:i2,j1:j2,km:1:-1)
!                                                                 ^^^
           CASE("QQJ")
             stOX_loss_dbl(i1:i2,j1:j2,1:km) = &
             stOX_loss_dbl(i1:i2,j1:j2,1:km) + self%rmult(i) * gmiQQJ(rxn_index)%pArray3D(i1:i2,j1:j2,km:1:-1)
!                                                                 ^^^
           CASE DEFAULT
             _ASSERT(.FALSE., TRIM(Iam)//': reaction must be QQK or QQJ : '//TRIM(one_name))

         END SELECT

       end do

       stOX_loss = stOX_loss_dbl
       DEALLOCATE( stOX_loss_dbl, __STAT__ )

     END IF  ! stOX_loss is needed

   END IF  ! QQJ and QQK

  RETURN
 END SUBROUTINE FillExports

!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  obtainBundles
!
! !INTERFACE:

  SUBROUTINE obtainBundles(state)
!
   USE ESMF

  implicit none
!
! !INPUT PARAMETERS:
  type(ESMF_State) , intent(in) :: state

!
! !DESCRIPTION:
! Extract arrays from the export state.
!
! !LOCAL VARIABLES:

      integer :: ib, numVars, rc
      type(ESMF_FieldBundle)      :: qjBundle
      type(ESMF_FieldBundle)      :: qkBundle
      type(ESMF_FieldBundle)      :: qqjBundle
      type(ESMF_FieldBundle)      :: qqkBundle
      type(ESMF_FieldBundle)      :: surfEmissBundle
      real(rPrec), pointer        :: ptr3D(:,:,:)
      real(rPrec), pointer        :: ptr2D(:,:)

      integer :: STATUS
      character(len=ESMF_MAXSTR), parameter :: IAm = "obtainBundless"
!
!EOP
!-------------------------------------------------------------------------------
      !==============
      ! the QJ Bundle
      !==============

      call ESMF_StateGet (state, "gmiQJ", qjBundle, RC=STATUS )
      VERIFY_(STATUS)

      call ESMF_FieldBundleGet(qjBundle, fieldCount=numVars , RC=STATUS)
      VERIFY_(STATUS)

      ALLOCATE(gmiQJ(numVars),STAT=STATUS)
      VERIFY_(STATUS)

      do ib = 1, numVars
         CALL obtainTracerFromBundle(qjBundle, ptr3D, ib)
         ALLOCATE(gmiQJ(ib)%pArray3D(i1:i2, j1:j2, 1:km))
         gmiQJ(ib)%pArray3D(:,:,km:1:-1) = ptr3D(:,:,1:km)
      end do

      !==============
      ! the QK Bundle
      !==============

      call ESMF_StateGet (state, "gmiQK", qkBundle, RC=STATUS )
      VERIFY_(STATUS)

      call ESMF_FieldBundleGet(qkBundle, fieldCount=numVars , RC=STATUS)
      VERIFY_(STATUS)

      ALLOCATE(gmiQK(numVars),STAT=STATUS)
      VERIFY_(STATUS)

      do ib = 1, numVars
         CALL obtainTracerFromBundle(qkBundle, ptr3D, ib)
         ALLOCATE(gmiQK(ib)%pArray3D(i1:i2, j1:j2, 1:km))
         gmiQK(ib)%pArray3D(:,:,km:1:-1) = ptr3D(:,:,1:km)
      end do

      IF (self%pr_qqjk) THEN
         !===============
         ! the QQJ Bundle
         !===============

         call ESMF_StateGet (state, "gmiQQJ", qqjBundle, RC=STATUS )
         VERIFY_(STATUS)

         call ESMF_FieldBundleGet(qqjBundle, fieldCount=numVars , RC=STATUS)
         VERIFY_(STATUS)

         ALLOCATE(gmiQQJ(numVars),STAT=STATUS)
         VERIFY_(STATUS)

         do ib = 1, numVars
            CALL obtainTracerFromBundle(qqjBundle, ptr3D, ib)
            ALLOCATE(gmiQQJ(ib)%pArray3D(i1:i2, j1:j2, 1:km))
            gmiQQJ(ib)%pArray3D(:,:,km:1:-1) = ptr3D(:,:,1:km)
         end do

         !===============
         ! the QQK Bundle
         !===============

         call ESMF_StateGet (state, "gmiQQK", qqkBundle, RC=STATUS )
         VERIFY_(STATUS)

         call ESMF_FieldBundleGet(qqkBundle, fieldCount=numVars , RC=STATUS)
         VERIFY_(STATUS)

         ALLOCATE(gmiQQK(numVars),STAT=STATUS)
         VERIFY_(STATUS)

         do ib = 1, numVars
            CALL obtainTracerFromBundle(qqkBundle, ptr3D, ib)
            ALLOCATE(gmiQQK(ib)%pArray3D(i1:i2, j1:j2, 1:km))
            gmiQQK(ib)%pArray3D(:,:,km:1:-1) = ptr3D(:,:,1:km)
         end do
      END IF

      !============================
      ! the Surface Emission Bundle
      !============================

      if (self%do_semiss_inchem) then
         call ESMF_StateGet (state, "surfEmissForChem", surfEmissBundle, RC=STATUS )
         VERIFY_(STATUS)

         call ESMF_FieldBundleGet(surfEmissBundle, fieldCount=numVars , rc=STATUS)
         VERIFY_(STATUS)

         Allocate(surfEmissForChem(i1:i2, j1:j2, 1:NSP))

         do ib = 1, numVars
            CALL obtainTracerFromBundle(surfEmissBundle, ptr2D, ib)

            surfEmissForChem(:,:,ib) = ptr2D(:,:)
         end do
      end if


      RETURN

  END SUBROUTINE obtainBundles

!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  FindPointers
!
! !INTERFACE:

  SUBROUTINE FindPointers(rc)
  
   USE ESMF

  IMPLICIT NONE

  INTEGER, INTENT(OUT) :: rc
!
! !DESCRIPTION:
!
!  Find pointers to import and export states
!
!EOP
!---------------------------------------------------------------------------

  CHARACTER(LEN=255) :: IAm
  REAL :: qmin,qmax
  
  rc=0
  IAm = "FindPointers"

!  Pointers to imports and exports
!  -------------------------------
   CALL MAPL_GetPointer(impChem,      zpbl,       'ZPBL', RC=STATUS)
   VERIFY_(STATUS)

   CALL MAPL_GetPointer(impChem,   airdens,    'AIRDENS', RC=STATUS)
   VERIFY_(STATUS)
   CALL MAPL_GetPointer(impChem,       ple,	   'PLE', RC=STATUS)
   VERIFY_(STATUS)
   CALL MAPL_GetPointer(impChem,         Q,          'Q', RC=STATUS)
   VERIFY_(STATUS)
   CALL MAPL_GetPointer(impChem,	 T,	     'T', RC=STATUS)
   VERIFY_(STATUS)
   CALL MAPL_GetPointer(impChem,       zle,	   'ZLE', RC=STATUS)
   VERIFY_(STATUS)

   CALL MAPL_GetPointer(expChem,    O3ppmv,     'O3PPMV', RC=STATUS)
   VERIFY_(STATUS)
   CALL MAPL_GetPointer(expChem,        O3,         'O3', RC=STATUS)
   VERIFY_(STATUS)
   CALL MAPL_GetPointer(expChem, stOX_loss,  'stOX_loss', RC=STATUS)
   VERIFY_(STATUS)

#include "Reactions_GetPointer___.h"

!  Validation
!  ----------
   Validate: IF(self%verbose) THEN
    IF(MAPL_AM_I_ROOT( ))PRINT *,TRIM(IAm),": Input ..."
    i = bxx%reg%nq
    CALL pmaxmin('TROPP:', bxx%qa(i)%data3d(:,:,km), qmin, qmax, iXj, 1, 0.01 )
    CALL pmaxmin('ZPBL:', zpbl, qmin, qmax, iXj, 1, 1. )
    CALL pmaxmin('Q:', Q, qmin, qmax, iXj, km, 1. )
    CALL pmaxmin('T:', T, qmin, qmax, iXj, km, 1. )
    CALL pmaxmin('ZLE:', zle, qmin, qmax, iXj, km+1, 1. )
    CALL pmaxmin('PLE (hPa):', ple, qmin, qmax, iXj, km+1, 0.01 )
    CALL pmaxmin('AIRDENS:', airdens, qmin, qmax, iXj, km, 1. )
   END IF Validate

!  Grab these data from the export state
!  -------------------------------------
   CALL getDataFromStateField(expChem,  HNO3GASsad,  'HNO3GASsad')

  RETURN
 END SUBROUTINE FindPointers

!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  SatisfyImports
!
! !INTERFACE:

  SUBROUTINE SatisfyImports(rc)
  
  IMPLICIT NONE

  INTEGER, INTENT(OUT) :: rc
!
! !DESCRIPTION:
!
!  Grab imports, change units if necessary, and convert to
!  REAL(KIND=DBL) when required.
!
!EOP
!---------------------------------------------------------------------------

  CHARACTER(LEN=255) :: IAm
  
  INTEGER, ALLOCATABLE :: on(:,:)
  
  rc=0
  IAm="SatisfyImports"

! ------------------------------------------------------------------------ 
! Imports to REAL*8, apply units conversion, and reverse vertical stacking
! ------------------------------------------------------------------------

! Singly-layered                                                            GEOS-5 Units       GMI Units
! The most recent valid tropopause pressures are stored in T2M15D(:,:,km)
! -----------------------------------------------------------------------   ------------       -------------
  i = bxx%reg%nq
  tropopausePress(i1:i2,j1:j2) = bxx%qa(i)%data3d(i1:i2,j1:j2,km)*Pa2hPa    ! Pa               hPa
  pctm2(i1:i2,j1:j2) = ple(i1:i2,j1:j2,km)*Pa2hPa                           ! Pa               hPa

! Layer means                                                               GEOS-5 Units       GMI Units
! -----------                                                               ------------       -------------
  DO k=1,km
   kReverse = km-k+1                                                        ! Lid-to-surf      Surf-to-lid
   press3c(i1:i2,j1:j2,kReverse) = pl(i1:i2,j1:j2,k)*Pa2hPa                 ! Pa               hPa
   kel(i1:i2,j1:j2,kReverse) = T(i1:i2,j1:j2,k)                             ! K
   humidity(i1:i2,j1:j2,kReverse) = Q(i1:i2,j1:j2,k)*ToGrPerKg              ! kg kg^{-1}       g kg^{-1}

   HNO3GAS(i1:i2,j1:j2,kReverse)  = HNO3GASsad(i1:i2,j1:j2,k)
  END DO

! Layer edges                                                               GEOS-5 Units       GMI Units
! -----------                                                               ------------       -------------
  DO k=0,km
   kReverse = km-k
   press3e(i1:i2,j1:j2,kReverse) = ple(i1:i2,j1:j2,k)*Pa2hPa                ! Pa               hPa
  END DO

! Cell mass and thickness                                                   GEOS-5 Units       GMI Units
! -----------------------                                                   ------------       -------------
   DO k=1,km
    kReverse = km-k+1
    mass(:,:,kReverse)=airdens(:,:,k)*self%cellArea(:,:)* &                 ! kg
                       (zle(:,:,k-1)-zle(:,:,k))
    gridBoxThickness(:,:,kReverse) = zle(:,:,k-1)-zle(:,:,k)                ! m
   END DO

  RETURN
 END SUBROUTINE SatisfyImports

 END SUBROUTINE GmiChemistry_GridCompRun

!-------------------------------------------------------------------------
!NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1, GEOS/DAS!
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  GmiChemistry_GridCompFinalize
!
! !INTERFACE:
!

   SUBROUTINE GmiChemistry_GridCompFinalize ( self, impChem, expChem, &
                                     nymd, nhms, cdt, rc )

  IMPLICIT none

! !INPUT/OUTPUT PARAMETERS:

   TYPE(GmiChemistry_GridComp), INTENT(inout) :: self ! Grid Component

! !INPUT PARAMETERS:

   INTEGER, INTENT(in) :: nymd, nhms	      ! time
   REAL,    INTENT(in) :: cdt  	              ! chemical timestep (secs)


! !OUTPUT PARAMETERS:

   TYPE(ESMF_State), INTENT(inout) :: impChem	! Import State
   TYPE(ESMF_State), INTENT(inout) :: expChem	! Import State
   INTEGER, INTENT(out) ::  rc                  ! Error return code:
                                                !  0 - all is well
                                                !  1 -
 
! !DESCRIPTION: This routine finalizes this Grid Component.
!
! !REVISION HISTORY:
!
!  18Sep2003 da Silva  First crack.
!
!EOP
!-------------------------------------------------------------------------

   CHARACTER(LEN=*), PARAMETER :: IAm = 'GmiChemistry_GridCompFinalize'
   INTEGER :: STATUS
   rc=0
   DEALLOCATE(self%cellArea, STAT=STATUS)
   VERIFY_(STATUS)
   RETURN

 END SUBROUTINE GmiChemistry_GridCompFinalize
  
 END MODULE GmiChem_GCCMod

