#include "MAPL_Generic.h"
!-------------------------------------------------------------------------
!NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1, GEOS/DAS!
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  GmiForcingBC_GCCMod --- GMI Forcing BC Grid Component Class
!
! Grid Component class for the GMI Forcing Boundary Conditions.
!
! !INTERFACE:
!

   MODULE  GmiForcingBC_GCCMod

! !USES:

   USE ESMF
   USE MAPL_Mod
   USE Chem_Mod 	     ! Chemistry Base Class
   USE Chem_UtilMod

   USE GmiSpcConcentrationMethod_mod, ONLY : t_SpeciesConcentration
   USE GmiGrid_mod,                   ONLY : t_gmiGrid
   USE GmiTimeControl_mod,            ONLY : t_GmiClock
   USE GmiESMFrcFileReading_mod,      ONLY : rcEsmfReadLogical
   USE GmiESMFrcFileReading_mod,      ONLY : rcEsmfReadTable
   use GmiArrayBundlePointer_mod,     only : t_GmiArrayBundle, CleanArrayPointer
   use GmiSpeciesRegistry_mod,        only : getSpeciesIndex, UNKNOWN_SPECIES
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
   PUBLIC  GmiForcingBC_GridComp       ! The GMI object 

! !PUBLIC MEMBER FUNCTIONS:

   PUBLIC  GmiForcingBC_GridCompInitialize
   PUBLIC  GmiForcingBC_GridCompRun
   PUBLIC  GmiForcingBC_GridCompFinalize

!
! !DESCRIPTION:
!
!  This module implements the Forcing Boundary Conditions.
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

  TYPE GmiForcingBC_GridComp
   CHARACTER(LEN=255) :: name = "GMI Forcing Boundary Conditions"

! Does the gmichem_import_rst file exist?  If not, some import states will be
! unfilled, and the species will need to "freewheel" through the first time step.
! -------------------------------------------------------------------------------
   LOGICAL :: gotImportRst

! Set BCRealTime = .TRUE. when boundary conditions 
! must be for exact year of current calendar date.
! -------------------------------------------------
   LOGICAL :: BCRealTime

    integer             :: forc_bc_opt
    integer             :: fbc_j1
    integer             :: fbc_j2
    integer             :: forc_bc_years
    integer             :: forc_bc_start_num
    integer             :: forc_bc_kmin
    integer             :: forc_bc_kmax
    integer             :: forc_bc_num
    integer             :: forc_bc_map (MAX_NUM_CONST)
    real*8              :: forc_bc_init_val
    real*8              :: forc_bc_incrpyr
    real*8              :: forc_bc_lz_val
    ! forcing boundary condition data (ppmv)
    real*8, pointer     :: forc_bc_data(:,:,:,:)

    integer             :: last_year
    integer, pointer    :: jlatmd(:,:)
    character (len=MAX_LENGTH_FILE_NAME) :: forc_bc_infile_name

! Various switches
! ----------------
   LOGICAL :: pr_diag
   LOGICAL :: do_synoz

! Dimensions
! ----------
   INTEGER :: i1, i2, im, j1, j2, jm, km

! Useful character strings
! ------------------------
   CHARACTER(LEN=255) :: chem_mecha

! Latitudes (radians)
! -------------------
   REAL, POINTER :: latRad(:,:)

! Extra diagnostics
! -----------------
   LOGICAL :: verbose

! Map GMI species indices to CCM indices
! --------------------------------------
   INTEGER, POINTER :: mapSpecies(:)

! Component derived type declarations
! -----------------------------------
   TYPE(t_gmiGrid   )		:: gmiGrid
   TYPE(t_GmiClock  )           :: gmiClock
   TYPE(t_SpeciesConcentration) :: SpeciesConcentration
 
  END TYPE GmiForcingBC_GridComp

CONTAINS

!-------------------------------------------------------------------------
!NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1, GEOS/DAS!
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  GmiForcingBC_GridCompInitialize --- Initialize GmiForcingBC_GridComp
!
! !INTERFACE:
!

   SUBROUTINE GmiForcingBC_GridCompInitialize( self, w_c, impChem, expChem, nymd, nhms, &
                                      tdt, rc )

   USE GmiSpcConcentrationMethod_mod, ONLY : InitializeSpcConcentration
   USE GmiGrid_mod,		      ONLY : InitializeGmiGrid
   USE GmiTimeControl_mod
   use GmiStringManipulation_mod, only : constructListNames
   use ReadForcedBC_mod, only : readForcedBcData

   IMPLICIT none

#     include "gmi_forc_bc.h"

   INTEGER, PARAMETER :: DBL = KIND(0.00D+00)

! !INPUT PARAMETERS:

   TYPE(Chem_Bundle), INTENT(in) :: w_c                ! Chemical tracer fields, delp, +
   INTEGER, INTENT(IN) :: nymd, nhms		       ! Time from AGCM
   REAL,    INTENT(IN) :: tdt			       ! Chemistry time step (secs)

! !OUTPUT PARAMETERS:

   TYPE(GmiForcingBC_GridComp), INTENT(INOUT)  :: self      ! Grid Component
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

   CHARACTER(LEN=*), PARAMETER :: IAm = 'GmiForcingBC_GridCompInitialize'
   CHARACTER(LEN=255) :: rcfilen = 'GMI_GridComp.rc'
   CHARACTER(LEN=255) :: namelistFile
   CHARACTER(LEN=255) :: kineticsTextFile
   CHARACTER(LEN=255) :: importRestartFile
   CHARACTER(LEN=255) :: string
      character (len=MAX_LENGTH_SPECIES_NAME), pointer :: tempListNames(:)
      character (len=MAX_STRING_LENGTH      ) :: forcedBcSpeciesNames
   
   type (ESMF_Config) :: gmiConfigFile

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
   INTEGER :: numSpecies
   INTEGER :: LogicalUnitNum

   INTEGER :: loc_proc, locGlobProc, commu_slaves
   LOGICAL :: one_proc, rootProc
   LOGICAL :: exists,open,found
   
   REAL :: qmin, qmax, tokgCPerBox

   integer :: start_ymd, idumday, month, year, ii, ij, jj
   INTEGER :: fbc_lat_num(1)
   REAL(KIND=DBL) :: fbc_lat(FBC_LATDIM)
   REAL(KIND=DBL) :: delta(FBC_LATDIM)
   REAL(KIND=DBL) :: dfbc_lat
   REAL(KIND=DBL), ALLOCATABLE :: latRad(:,:)

! Grid cell area can be set by initialize
! ---------------------------------------
   REAL, POINTER, DIMENSION(:,:) :: cellArea

   integer           :: ib
   character (len=4) :: binName

   self%name = 'GMI Forcing Boundary Conditions'

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

   rootProc=.FALSE.
   IF( MAPL_AM_I_ROOT() ) THEN
    rootProc=.TRUE.
   END IF 

     !-------------------
     ! Load resource file
     !-------------------

      IF ( MAPL_AM_I_ROOT() ) THEN
         PRINT *," "
         PRINT *,TRIM(IAm)//":"
         PRINT *,"Starting Reading the GMI Resource File for ForcingBC"
      ENDIF

      gmiConfigFile = ESMF_ConfigCreate(rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigLoadFile(gmiConfigFile, TRIM(rcfilen), rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(gmiConfigFile, importRestartFile, &
     &                label   = "importRestartFile:", &
     &                default = ' ', rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(gmiConfigFile, self%chem_mecha, &
     &                label   = "chem_mecha:", &
     &                default = 'strat_trop', rc=STATUS )
      VERIFY_(STATUS)

      !------------------------------
      ! Diagnostics related variables
      !------------------------------

      call rcEsmfReadLogical(gmiConfigFile, self%pr_diag, &
     &           "pr_diag:", default=.false., rc=STATUS)
      VERIFY_(STATUS)

      call rcEsmfReadLogical(gmiConfigFile, self%verbose, &
     &           "verbose:", default=.false., rc=STATUS)
      VERIFY_(STATUS)

      call rcEsmfReadLogical(gmiConfigFile, self%do_synoz, &
     &           "do_synoz:", default=.false., rc=STATUS)
      VERIFY_(STATUS)

!     ---------------------------
!     Forcing boundary condition:
!     ---------------------------

!     --------------------------------------------------
!     forc_bc_opt
!	1:  Time-slice     Annually repeating, corresponding to forc_bc_start_num.
!	2:  Time-dependent Start at forc_bc_start_num, subsequent year at forc_bc_start_num+1, etc.
!                          WARNING: Limited by range of years in forc_bc_infile_name.
!	3:  Calculate	   No horizontal variation.
!     --------------------------------------------------

      call ESMF_ConfigGetAttribute(gmiConfigFile, self%forc_bc_opt, &
     &                label   = "forc_bc_opt:", &
     &                default = 1, rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(gmiConfigFile, self%fbc_j1, &
     &                label   = "fbc_j1:", &
     &                default = ju1_gl, rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(gmiConfigFile, self%fbc_j2, &
     &                label   = "fbc_j2:", &
     &                default = j2_gl, rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(gmiConfigFile, self%forc_bc_years, &
     &                label   = "forc_bc_years:", &
     &                default = 1, rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(gmiConfigFile, self%forc_bc_start_num, &
     &                label   = "forc_bc_start_num:", &
     &                default = 1, rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(gmiConfigFile, self%forc_bc_kmin, &
     &                label   = "forc_bc_kmin:", &
     &                default = 1, rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(gmiConfigFile, self%forc_bc_kmax, &
     &                label   = "forc_bc_kmax:", &
     &                default = 1, rc=STATUS )
      VERIFY_(STATUS)

      self%forc_bc_map(:)      = 0

      call rcEsmfReadTable(gmiConfigFile, forcedBcSpeciesNames, &
     &                     "forcedBcSpeciesNames::", rc=STATUS)

      call ESMF_ConfigGetAttribute(gmiConfigFile, self%forc_bc_init_val, &
     &                label   = "forc_bc_init_val:", &
     &                default = 0.0d0, rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(gmiConfigFile, self%forc_bc_incrpyr, &
     &                label   = "forc_bc_incrpyr:", &
     &                default = 0.3d0, rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(gmiConfigFile, self%forc_bc_lz_val, &
     &                label   = "forc_bc_lz_val:", &
     &                default = 0.0d0, rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(gmiConfigFile, self%forc_bc_infile_name, &
     &                label   = "forc_bc_infile_name:", &
     &                default = 'forc_bc_co2.asc', rc=STATUS )
      VERIFY_(STATUS)

      ! Set the initial value of the list
      allocate(tempListNames(NSP))
      tempListNames(:) = ''

      ! Construct the list of names using the long string
      call constructListNames(tempListNames, forcedBcSpeciesNames)

      self%forc_bc_num = Count (tempListNames(:) /= '')
      if (self%forc_bc_num > 0) then
         do ic = 1, self%forc_bc_num
            self%forc_bc_map(ic) = getSpeciesIndex(tempListNames(ic))
         end do

         IF (rootProc) THEN
            PRINT *," "
            PRINT *,"Number of GHG and ODP gases: ",self%forc_bc_num
            PRINT *," "
            PRINT *,"forc_bc_map:"
            PRINT *, self%forc_bc_map
            PRINT *," "
         END IF

      end if

      deallocate(tempListNames)

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
   
! Set GMI's clock
! ---------------
   CALL Set_begGmiDate(self%gmiClock, nymd)
   CALL Set_begGmiTime(self%gmiClock, nhms)
   CALL Set_curGmiDate(self%gmiClock, nymd)
   CALL Set_curGmiTime(self%gmiClock, nhms)
   CALL Set_numTimeSteps(self%gmiClock, 0)

! Discretization
! --------------
   CALL InitializeGmiGrid(self%gmiGrid, NPIJ, NPI, NPJ, &
                          gmi_nborder, i1, i2, ju1, jv1, j2, k1, k2, &
                          i1_gl, i2_gl, ju1_gl, jv1_gl, j2_gl, k1_gl, k2_gl, &
                          ilo, ihi, julo, jvlo, jhi, &
                          ilo_gl, ihi_gl, julo_gl, jvlo_gl, jhi_gl, &
                          ilong, ilat, ivert, itloop, j1p, j2p)  

! Number of species and perform a consistency check with setkin_par.h.
! NOTES:
!  1. H2O is specie number 10 in the strat-trop mechanism, but will not be
!     found in w_c%reg%vname. H2O will be initialized from specific humidity, Q.
!  2. The GEOS-5 bundle has an Age-Of-Air tracer, which is not carried by GMI.
!  3. At the end of the XX (non-transported) species is a place holder for T2M15d.
! So w_c%reg%j_XX-w_c%reg%i_GMI must equal the parameter NSP = NCONST + NDYN.
! --------------------------------------------------------------------------------
   numSpecies = w_c%reg%j_XX-w_c%reg%i_GMI
   IF(numSpecies /= NSP) THEN
    PRINT *,TRIM(IAm),': Number of species from Chem_Registry.rc does not match number in setkin_par.h'
    STATUS = 1
    VERIFY_(STATUS)
    RETURN
   END IF

! Allocate space, etc., but the initialization of the
! species from the internal state is left to the run method.
! ----------------------------------------------------------

      CALL InitializeSpcConcentration(self%SpeciesConcentration,              &
     &               self%gmiGrid, gmiConfigFile, numSpecies, NMF, NCHEM,     &
     &               loc_proc)

      ! Boundary forcing Data
      if (self%forc_bc_num > 0) then
         Allocate (self%forc_bc_data(FBC_LATDIM, FBC_MONDIM,  &
         &            self%forc_bc_years, self%forc_bc_num))
         self%forc_bc_data = 0.0d0 
         call readForcedBcData (self%pr_diag, loc_proc, self%forc_bc_opt, &
                             self%forc_bc_years, self%forc_bc_num, &
                             self%forc_bc_data, &
                             self%forc_bc_init_val, self%forc_bc_infile_name)
      end if

      if (self%forc_bc_opt <= 2) then

          Allocate (self%jlatmd(i1:i2,ju1:j2))
          self%jlatmd(i1:i2,ju1:j2) = 0

          Allocate (latRad(i1:i2,ju1:j2))
          latRad(:,:) = self%latRad(:,:)

          call Get_begGmiDate  (self%gmiClock, start_ymd)
          call GmiSplitDateTime (start_ymd, year, month, idumday)

          self%last_year = year

          ! -------------------------------------------------------------
          ! The forc_bc values are known only at certain locations.  The
          ! following code determines which of these locations is
          ! closest to the simulation latitude zone.  No interpolation is
          ! done, just pick closest.
          ! -------------------------------------------------------------

          ! -------------------------------------------------------------
          ! Find latitudes (degrees) of forc_bc values, assuming index 1
          ! is at the South Pole and constant increments to the North 
          ! Pole.
          ! -------------------------------------------------------------
          dfbc_lat = 180.00D+00/(FBC_LATDIM-1.00D+00)
          DO ij = 1,FBC_LATDIM
             fbc_lat(ij) = -90.00D+00 + dfbc_lat*(ij-1)
          END DO

          ! -------------------------------------------------------------
          ! For each latitude circle on this processor, (1) find its
          ! displacement from EACH latitude of forc_bc values, and 
          ! then (2) find which index holds the minimum displacement. 
          ! -------------------------------------------------------------
          DO jj = ju1,j2
             DO ii = i1,i2
                delta(1:FBC_LATDIM) = ABS(fbc_lat(1:FBC_LATDIM)-latRad(ii,jj)*DEGPRAD)
                fbc_lat_num = MINLOC(delta)
                self%jlatmd(ii,jj) = fbc_lat_num(1)
             END DO
          END DO
          deallocate(latRad)

      end if

    !---------------------------------------------------------------
    ! Create and populate the array that maps GMI species indices to 
    ! GEOS-5 species indices
    !---------------------------------------------------------------

    allocate(self%mapSpecies(NSP))
    self%mapSpecies(:) = speciesReg_for_CCM(lchemvar, w_c%reg%vname, NSP, w_c%reg%i_GMI, w_c%reg%j_XX)

  RETURN

  END SUBROUTINE GmiForcingBC_GridCompInitialize

!-------------------------------------------------------------------------
!NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1, GEOS/DAS!
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  GmiForcingBC_GridCompRun --- The GMI Driver
!
! !INTERFACE:
!

   SUBROUTINE GmiForcingBC_GridCompRun ( self, w_c, impChem, expChem, nymd, nhms, &
                                tdt, rc )

! !USES:

   USE GmiTimeControl_mod,            ONLY : Set_curGmiDate, Set_curGmiTime
   USE GmiTimeControl_mod,            ONLY : Set_numTimeSteps, Get_numTimeSteps
   USE GmiTimeControl_mod,            ONLY : Set_gmiSeconds, GetSecondsFromJanuary1
   USE GmiTimeControl_mod,            ONLY : Get_gmiSeconds, Get_begGmiDate
   USE GmiSpcConcentrationMethod_mod, ONLY : resetFixedConcentration
   USE GmiSolar_mod,                  ONLY : CalcCosSolarZenithAngle
   use GmiUpdateForcingBC_mod       , only : updateForcingBC

   IMPLICIT none

! !INPUT/OUTPUT PARAMETERS:

   TYPE(GmiForcingBC_GridComp), INTENT(INOUT) :: self ! Grid Component
   TYPE(Chem_Bundle), INTENT(INOUT) :: w_c    ! Chemical tracer fields   

! !INPUT PARAMETERS:

   TYPE(ESMF_State), INTENT(INOUT) :: impChem ! Import State
   INTEGER, INTENT(IN) :: nymd, nhms	      ! time
   REAL,    INTENT(IN) :: tdt		      ! chemical timestep (secs)

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

   CHARACTER(LEN=*), PARAMETER :: IAm = 'GmiForcingBC_GridCompRun'

!  Input fields from GEOS-5
!  ------------------------
   REAL, POINTER, DIMENSION(:,:,:) :: Q

!  Local
!  -----
   INTEGER :: cymd, dymd, hms
   INTEGER :: i, i1, i2, ic, im, iXj, iTile(1)
   INTEGER :: j, j1, j2, jm, jTile(1)
   INTEGER :: k, km
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
   REAL(KIND=DBL) :: gmi_sec
   integer        :: start_ymd

   CHARACTER(LEN=255) :: speciesName
   CHARACTER(LEN=255) :: importName

   LOGICAL :: found, rootProc
   LOGICAL, PARAMETER :: doThis = .FALSE.
   
! Allocatables.  Use KIND=DBL where GMI expects REAL*8.
! -----------------------------------------------------
   REAL, ALLOCATABLE :: var3d(:,:,:)
   REAL(KIND=DBL), ALLOCATABLE :: latRad(:,:)

   loc_proc = -99

!  Grid specs
!  ----------
   rc = 0
   i1 = self%i1
   i2 = self%i2
   im = self%im

   j1 = self%j1
   j2 = self%j2
   jm = self%jm

   km = self%km

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
   ALLOCATE(latRad(i1:i2,j1:j2),STAT=STATUS)
   VERIFY_(STATUS)
   ALLOCATE(var3d(i1:i2,j1:j2,1:km),STAT=STATUS)
   VERIFY_(STATUS)

! Latitudes (radians)
! -------------------
   latRad(i1:i2,j1:j2)=self%latRad(i1:i2,j1:j2)

! Set GMI's clock
! ---------------
   CALL Set_curGmiDate(self%gmiClock, nymd)
   CALL Set_curGmiTime(self%gmiClock, nhms)
   CALL Get_numTimeSteps(self%gmiClock, ic)
   CALL Set_numTimeSteps(self%gmiClock, ic+1)
   CALL Set_gmiSeconds(self%gmiClock, (ic+1)*chemDt)

   call Get_begGmiDate  (self%gmiClock, start_ymd)
   call Get_gmiSeconds  (self%gmiClock, gmi_sec)

! Update the following time-dependent boundary conditions:
!  Fixed concentration species
! --------------------------------------------------------

   CALL Acquire_Clims(STATUS)
   VERIFY_(STATUS)

! Hand the species concentrations to GMI's bundle
! -----------------------------------------------
   IF (self%gotImportRst) then
      CALL SwapSpeciesBundles(ToGMI, self%SpeciesConcentration%concentration, &
               w_c%qa, Q, self%mapSpecies, lchemvar, self%do_synoz, NSP, &
               STATUS)
      VERIFY_(STATUS)
   END IF

   DEALLOCATE(var3D, STAT=STATUS)
   VERIFY_(STATUS)

! Impose fixed concentrations
! ---------------------------

   IF(self%SpeciesConcentration%num_fixed_const > 0) THEN
    CALL resetFixedConcentration(self%SpeciesConcentration, self%gmiClock, self%gmiGrid, NSP)
   END IF

! ------------------------------------------------------------------------
! Update Forcing Boundary Conditions
! ------------------------------------------------------------------------

      IF (self%gotImportRst) then
         if (self%forc_bc_num > 0) then
            call updateForcingBC (self%forc_bc_data, &
                       self%SpeciesConcentration%concentration, self%jlatmd, &
     &                 self%last_year, nymd, gmi_sec, self%fbc_j1, self%fbc_j2, &
     &                 self%forc_bc_num, self%forc_bc_kmax, self%forc_bc_kmin, &
     &                 self%forc_bc_opt, self%forc_bc_map, self%forc_bc_incrpyr, &
     &                 self%forc_bc_start_num, self%forc_bc_years, self%gmiGrid%i1, &
     &                 self%gmiGrid%i2, self%gmiGrid%ju1, self%gmiGrid%j2,       &
     &                 self%gmiGrid%k1, self%gmiGrid%k2, NSP)
         end if
      END IF

! Return species concentrations to the chemistry bundle
! -----------------------------------------------------
   IF (self%gotImportRst) then
      CALL SwapSpeciesBundles(FromGMI, self%SpeciesConcentration%concentration, &
               w_c%qa, Q, self%mapSpecies, lchemvar, self%do_synoz, NSP,  &
               STATUS)
      VERIFY_(STATUS)
   END IF

! Export states
! -------------

   CALL FillExports(STATUS)
   VERIFY_(STATUS)

! Scratch local work space
! ------------------------
   DEALLOCATE(latRad, STAT=STATUS)
   VERIFY_(STATUS)

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
!
!EOP
!---------------------------------------------------------------------------

  CHARACTER(LEN=255) :: IAm
  INTEGER :: lightning_opt
  INTEGER, ALLOCATABLE :: flag(:)
  
  REAL, POINTER, DIMENSION(:,:,:) :: PTR3D
  
  rc = 0
  IAm = "Acquire_Clims"

! Fixed concentration species
! ---------------------------
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

  CHARACTER(LEN=255) :: IAm
  
  rc=0
  IAm="FillExports"


  RETURN
 END SUBROUTINE FillExports

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
  IAm="FindPointers"

!  Pointers to imports
!  -------------------
   CALL MAPL_GetPointer(impChem, Q, 'Q', RC=STATUS)
   VERIFY_(STATUS)

!  Validation
!  ----------
   Validate: IF(self%verbose) THEN
    IF(MAPL_AM_I_ROOT( ))PRINT *,TRIM(IAm),": Input ..."
    CALL pmaxmin('Q:', Q, qmin, qmax, iXj, km, 1. )
   END IF Validate


  RETURN
 END SUBROUTINE FindPointers

 END SUBROUTINE GmiForcingBC_GridCompRun

!-------------------------------------------------------------------------
!NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1, GEOS/DAS!
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  GmiForcingBC_GridCompFinalize
!
! !INTERFACE:
!

   SUBROUTINE GmiForcingBC_GridCompFinalize ( self, w_c, impChem, expChem, &
                                     nymd, nhms, cdt, rc )

  IMPLICIT none

! !INPUT/OUTPUT PARAMETERS:

   TYPE(GmiForcingBC_GridComp), INTENT(inout) :: self ! Grid Component

! !INPUT PARAMETERS:

   TYPE(Chem_Bundle), INTENT(in)  :: w_c      ! Chemical tracer fields   
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

   CHARACTER(LEN=*), PARAMETER :: IAm = 'GmiForcingBC_GridCompFinalize'
   INTEGER :: STATUS
   rc=0
   DEALLOCATE(self%latRad, STAT=STATUS)
   VERIFY_(STATUS)
   RETURN

 END SUBROUTINE GmiForcingBC_GridCompFinalize
  
 END MODULE GmiForcingBC_GCCMod

