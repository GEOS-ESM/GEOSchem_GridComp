#include "MAPL_Generic.h"
!-------------------------------------------------------------------------
!NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1, GEOS/DAS!
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  GmiSAD_GCCMod --- GMI SAD Grid Component Class
!
! Grid Component class for Surface Area Densities (SAD) module.
!
! !INTERFACE:
!

   MODULE  GmiSAD_GCCMod

! !USES:

   USE ESMF
   USE MAPL
   USE Chem_Mod            ! SAD Base Class
   USE Chem_UtilMod

   USE GmiSpcConcentrationMethod_mod, ONLY : t_SpeciesConcentration
   USE GmiGrid_mod,                   ONLY : t_gmiGrid
   USE GmiTimeControl_mod,            ONLY : t_GmiClock
   USE GmiESMFrcFileReading_mod,      ONLY : rcEsmfReadLogical
   USE GmiESMFrcFileReading_mod,      ONLY : rcEsmfReadTable
   use GmiArrayBundlePointer_mod, only : t_GmiArrayBundle
   use GmiPrintError_mod,        only : GmiPrintError
   use GmiCheckNamelistFile_mod, only : CheckNamelistOptionRange
   use GmiFieldBundleESMF_mod,        ONLY : addTracerToBundle
   use GmiFieldBundleESMF_mod,        ONLY : updateTracerToBundle
   use GmiStateFieldESMF_mod,         ONLY : setDataToStateField
   use GmiStateFieldESMF_mod,         ONLY : initDataInStateField
   use GmiSwapSpeciesBundlesMod,      ONLY : SwapSpeciesBundles, speciesReg_for_CCM

   IMPLICIT NONE
   INTEGER, PARAMETER :: DBL = KIND(0.00D+00)

#include "setkin_par.h"
#include "GmiParameters.h"
#include "gmi_phys_constants.h"
#include "setkin_mw.h"
#include "setkin_lchem.h"
#include "gmi_time_constants.h"

! !TYPES:

   PRIVATE
   PUBLIC  GmiSAD_GridComp       ! The GMI object 

! !PUBLIC MEMBER FUNCTIONS:

   PUBLIC  GmiSAD_GridCompInitialize
   PUBLIC  GmiSAD_GridCompRun
   PUBLIC  GmiSAD_GridCompFinalize

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

  TYPE GmiSAD_GridComp
   CHARACTER(LEN=255) :: name = "GMI Surface Area Densities"

! Does the gmichem_import_rst file exist?  If not, some import states will be
! unfilled, and the species will need to "freewheel" through the first time step.
! -------------------------------------------------------------------------------
   LOGICAL :: gotImportRst

! PSC chemistry is not run in latitude (deg) range +/- NoPSCZone.  Default: Zero.
! -------------------------------------------------------------------------------
   INTEGER :: NoPSCZone

! Upper pressure limit for PSC chemistry (hPa).  If <= 0, tropopause pressure is
! used.  Otherwise, the MIN of tropopause pressure and PSCMaxP is chosen.  Defaut: Zero.
! --------------------------------------------------------------------------------------
   INTEGER :: PSCMaxP

! Enforce a maximum allowable condensed HNO3 (ppbv).  Default: 25
! ---------------------------------------------------------------
   INTEGER :: HNO3Ice_MAX

! Enforce a maximum allowable HCl (ppbv).  Default: 5
! ---------------------------------------------------------------
   REAL :: HCl_MAX

! Various switches
! ----------------
   LOGICAL :: pr_sad
   LOGICAL :: pr_diag
   LOGICAL :: do_synoz
   LOGICAL :: do_qqjk_inchem
   LOGICAL :: do_qqjk_reset
   LOGICAL :: pr_qqjk

! Dimensions
! ----------
   INTEGER :: i1, i2, im, j1, j2, jm, km

! Useful character strings
! ------------------------
   CHARACTER(LEN=255) :: chem_mecha
   character (len=MAX_LENGTH_MET_NAME) :: metdata_name_model

! Longitudes and latitudes (radians)
! ----------------------------------
   REAL, POINTER :: lonRad(:,:)
   REAL, POINTER :: latRad(:,:)

! Extra diagnostics
! -----------------
   LOGICAL :: verbose

! Map GMI species indices to CCM indices
! --------------------------------------
   INTEGER, POINTER :: mapSpecies(:)

    real(KIND=DBL)          :: dehydmin ! minimum dehyd value (mixing ratio)
    integer                 :: ihno3cond_num  
    integer                 :: idehyd_num
    integer                 :: ih2oaircr_num
    integer                 :: ich4_num  
    integer                 :: ihno3_num  
    integer                 :: ih2o_num  
    integer                 :: num_sad

    type (t_GmiArrayBundle), pointer :: sadgmi(:) => null()
    real(KIND=DBL), pointer :: lbssad   (:,:,:,:) => null()
    real(KIND=DBL), pointer :: h2oclim  (:,:,:,:) => null()
    real(KIND=DBL), pointer :: ch4clim  (:,:,:,:) => null()
    real(KIND=DBL), pointer :: loss_freq(:,:,:,:) => null()
    real(KIND=DBL), pointer :: hno3cond (:,:,:)   => null()
    real(KIND=DBL), pointer :: hno3gas  (:,:,:)   => null()
    real(KIND=DBL), pointer :: h2oback  (:,:,:)   => null()
    real(KIND=DBL), pointer :: h2ocond  (:,:,:)   => null()
    real(KIND=DBL), pointer :: reffice  (:,:,:)   => null()
    real(KIND=DBL), pointer :: reffsts  (:,:,:)   => null()
    real(KIND=DBL), pointer :: vfall    (:,:,:)   => null()

    integer             :: dehyd_opt
    integer             :: sad_opt, chem_opt
    integer             :: h2oclim_opt
    integer             :: h2oclim_timpyr
    real*8              :: ch4clim_init_val
    real*8              :: h2oclim_init_val
    integer             :: lbssad_opt
    integer             :: lbssad_timpyr
    real*8              :: lbssad_init_val

    ! This is to store a generic name for each SAD entry
    character(len=ESMF_MAXSTR), pointer :: sadName(:) => null()

! Component derived type declarations
! -----------------------------------
   TYPE(t_gmiGrid   )		:: gmiGrid
   TYPE(t_GmiClock  )           :: gmiClock
   TYPE(t_SpeciesConcentration) :: SpeciesConcentration
 
  END TYPE GmiSAD_GridComp

CONTAINS

!-------------------------------------------------------------------------
!NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1, GEOS/DAS!
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  GmiSAD_GridCompInitialize --- Initialize GmiSAD_GridComp
!
! !INTERFACE:
!

   SUBROUTINE GmiSAD_GridCompInitialize( self, w_c, impChem, expChem, nymd, nhms, &
                                      tdt, rc )

   USE GmiSpcConcentrationMethod_mod, ONLY : InitializeSpcConcentration
   USE GmiGrid_mod,		      ONLY : InitializeGmiGrid
   USE GmiTimeControl_mod,            ONLY : Set_curGmiDate, Set_curGmiTime
   USE GmiTimeControl_mod,            ONLY : Set_begGmiDate, Set_begGmiTime
   USE GmiTimeControl_mod,            ONLY : Set_numTimeSteps

   IMPLICIT none
   INTEGER, PARAMETER :: DBL = KIND(0.00D+00)

! !INPUT PARAMETERS:

   TYPE(Chem_Bundle), INTENT(in) :: w_c                ! Chemical tracer fields, delp, +
   INTEGER, INTENT(IN) :: nymd, nhms		       ! Time from AGCM
   REAL,    INTENT(IN) :: tdt			       ! time step (secs)

! !OUTPUT PARAMETERS:

   TYPE(GmiSAD_GridComp), INTENT(INOUT)  :: self      ! Grid Component
   TYPE(ESMF_State),   INTENT(INOUT)  :: impChem    ! Import State
   TYPE(ESMF_State),   INTENT(INOUT)  :: expChem    ! Export State

   INTEGER, INTENT(out) ::  rc                  ! Error return code:
                                                !  0 - all is well
                                                !  1 - 

! !DESCRIPTION: Initializes the GMI SAD Grid Component. It primarily sets
!               the import state.
!
! !REVISION HISTORY:
!
!  18Sep2003 da Silva  First crack.
!  30Jun2007 Nielsen   GMI Combo set up
!
!EOP
!-------------------------------------------------------------------------

   CHARACTER(LEN=*), PARAMETER :: IAm = 'GmiSAD_GridCompInitialize'
   CHARACTER(LEN=255) :: rcfilen = 'GMI_GridComp.rc'
   CHARACTER(LEN=255) :: importRestartFile
   CHARACTER(LEN=255) :: string, err_msg
   
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
   INTEGER :: numSpecies, ib
   INTEGER :: LogicalUnitNum

   character (len=4) :: binName

   INTEGER :: loc_proc, locGlobProc, commu_slaves
   LOGICAL :: one_proc, rootProc
   LOGICAL :: exists,open,found
   
   REAL :: qmin, qmax, tokgCPerBox
   integer  :: dims1(3), dims2(3)

   REAL(rPrec), POINTER, DIMENSION(:,:,:) :: HNO3CONDsad
   REAL(rPrec), POINTER, DIMENSION(:,:,:) :: HNO3GASsad
   real(rPrec), pointer :: var(:,:,:)

   type(ESMF_FieldBundle)      :: sadBun
   integer                     :: numVars

   self%name = 'GMI Surface Area Densities'

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
         PRINT *,"Starting Reading the GMI Resource File for SAD"
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

      call ESMF_ConfigGetAttribute(gmiConfigFile, self%metdata_name_model, &
     &                label   = "metdata_name_model:", &
     &                default = 'GEOS-5', rc=STATUS )

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

      call ESMF_ConfigGetAttribute(gmiConfigFile, self%chem_opt, &
     &                label   = "chem_opt:", &
     &                default = 2, rc=STATUS )
      VERIFY_(STATUS)

      !---------------------------------------------------
      ! PSC exclusion zone, +/- NoPSCZone degrees latitude
      !---------------------------------------------------
      CALL ESMF_ConfigGetAttribute(gmiConfigFile, self%NoPSCZone, &
     &                LABEL   = "NoPSCZone:", &
     &                DEFAULT = 0, RC=STATUS )
      VERIFY_(STATUS)

      !--------------------------------------------------
      ! Highest pressure at which PSCs are computed. If <= 0, 
      ! tropopause pressure is used.  Otherwise, the MIN of 
      ! PSCMaxP and the tropopause pressure is chosen.
      !--------------------------------------------------
      CALL ESMF_ConfigGetAttribute(gmiConfigFile, self%PSCMaxP, &
     &                LABEL   = "PSC_Max_P_hPa:", &
     &                DEFAULT = 0, RC=STATUS )
      VERIFY_(STATUS)

      !--------------------------------------------------
      ! Enforce a maximum allowable condensed HNO3 (ppbv)
      !--------------------------------------------------
      CALL ESMF_ConfigGetAttribute(gmiConfigFile, self%HNO3Ice_MAX, &
     &                LABEL   = "Condensed_HNO3_limit:", &
     &                DEFAULT = 25, RC=STATUS )
      VERIFY_(STATUS)

      !--------------------------------------------------
      ! Enforce a maximum allowable HCl (ppbv)
      !--------------------------------------------------
      CALL ESMF_ConfigGetAttribute(gmiConfigFile, self%HCl_MAX, &
     &                LABEL   = "HCl_limit:", &
     &                DEFAULT = 5.0, RC=STATUS )
      VERIFY_(STATUS)

!     ----------------------------------------------------
!     sad_opt
!       0:  do not allocate or process SAD array
!       1:  allocate, but zero out SAD array
!       2:  call Considine code (i.e., Condense)
!       3:  read SAD array from a file of monthly averages
!     ----------------------------------------------------

      call ESMF_ConfigGetAttribute(gmiConfigFile, self%sad_opt, &
     &                label   = "sad_opt:", &
     &                default = 0, rc=STATUS )
      VERIFY_(STATUS)

!     ------------------------------------------------
!     h2oclim_opt
!       1:  set all h2oclim values to h2oclim_init_val
!       2:  read in h2oclim
!       3:  h2oclim, ch4clim not used.  Instead, transported
!           H2O and CH4 are provided by the host AGCM.
!     ------------------------------------------------

      call ESMF_ConfigGetAttribute(gmiConfigFile, self%h2oclim_opt, &
     &                label   = "h2oclim_opt:", &
     &                default = 2, rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(gmiConfigFile, self%h2oclim_timpyr, &
     &                label   = "h2oclim_timpyr:", &
     &                default = MONTHS_PER_YEAR, rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(gmiConfigFile, self%ch4clim_init_val, &
     &                label   = "ch4clim_init_val:", &
     &                default = 0.0d0, rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(gmiConfigFile, self%h2oclim_init_val, &
     &                label   = "h2oclim_init_val:", &
     &                default = 0.0d0, rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(gmiConfigFile, self%pr_sad, &
     &           "pr_sad:", default=.false., rc=STATUS)
      VERIFY_(STATUS)

!     ------------------------------------------------
!     dehyd_opt
!       0:  dehyd is not available and is set to zero.
!           In general this is the AGCM case where H20 and
!           condensed H2O are transported species.
!       1:  dehyd is a transported specie and is available.
!           In general, this is the CTM case.
!     ------------------------------------------------

      call ESMF_ConfigGetAttribute(gmiConfigFile, self%dehyd_opt, &
     &                label   = "dehyd_opt:", &
     &                default = 1, rc=STATUS )
      VERIFY_(STATUS)

!     ----------------------------------------------
!     lbssad_opt
!       1:  set all lbssad values to lbssad_init_val
!       2:  read in lbssad 3d fields
!       3:  read in lbssad zonal average fields
!       4:  lbssad provided by AGCM, current month only.
!     ----------------------------------------------

      call ESMF_ConfigGetAttribute(gmiConfigFile, self%lbssad_opt, &
     &                label   = "lbssad_opt:", &
     &                default = 2, rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(gmiConfigFile, self%lbssad_timpyr, &
     &                label   = "lbssad_timpyr:", &
     &                default = MONTHS_PER_YEAR, rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(gmiConfigFile, self%lbssad_init_val, &
     &                label   = "lbssad_init_val:", &
     &                default = 0.0d0, rc=STATUS )
      VERIFY_(STATUS)

      call CheckNamelistOptionRange ('h2oclim_opt', self%h2oclim_opt, 1, 3)
      call CheckNamelistOptionRange ('lbssad_opt', self%lbssad_opt, 1, 4)
      call CheckNamelistOptionRange ('sad_opt', self%sad_opt, 0, 3)
      call CheckNamelistOptionRange ('dehyd_opt', self%dehyd_opt, 0, 1)

      IF (self%lbssad_opt == 4) self%lbssad_timpyr = 1

      if (self%sad_opt /= 0) then
         self%num_sad = NSAD

         allocate (self%sadName(NSAD))
         do ib = 1, NSAD
            write (binName ,'(i4.4)') ib 
            self%sadName(ib) = 'gmiSAD'//binName
         end do

      end if

      self%ich4_num      = ICH4
      self%ih2o_num      = IH2O
      self%idehyd_num    = 0                               ! JK 12June2010
      if (self%dehyd_opt == 0) self%idehyd_num    = IDEHYD ! JK 12June2010
      self%ihno3_num     = IHNO3
      self%ihno3cond_num = IHNO3COND    ! JK 12June2010

      !===============
      ! Error Checking
      !===============

      if (self%pr_sad .and. ((self%sad_opt == 0) .or. (self%num_sad == 0))) then
         err_msg = 'pr_sad/sad_opt/num_sad problem.'
         call GmiPrintError  &
     &    (err_msg, .true., 2, self%sad_opt, self%num_sad, 0, 0.0d0, 0.0d0)
      end if

      if ((self%sad_opt == 1) .or. (self%sad_opt == 2)) then

         if ((self%ich4_num == 0) .or. (self%ich4_num /= ICH4)) then
            err_msg = 'chem_opt/ich4_num problem in the rc File.'
            call GmiPrintError  &
     &        (err_msg, .true., 2, self%chem_opt, self%ich4_num, 0, 0.0d0, 0.0d0)
         end if

!.sds.. unnecessary test
!.sds..         if ((self%idehyd_num == 0) .or. (self%idehyd_num /= IDEHYD)) then
!.sds..            err_msg = 'chem_opt/idehyd_num problem in the rc File.'
!.sds..            call GmiPrintError  &
!.sds..     &       (err_msg, .true., 2, self%chem_opt, self%idehyd_num, 0, 0.0d0, 0.0d0)
!.sds..         end if

         if ((self%ih2o_num == 0) .or. (self%ih2o_num /= IH2O)) then
            err_msg = 'chem_opt/ih2o_num problem in the rc File.'
            call GmiPrintError  &
     &        (err_msg, .true., 2, self%chem_opt, self%ih2o_num, 0, 0.0d0, 0.0d0)
         end if

      end if

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
   i1_gl = i1
   i2_gl = i2
   ju1_gl = i1
   jv1_gl = j1
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

      !===================
      ! Allocate Variables
      !===================

      if (self%sad_opt /= 0) then
         Allocate(self%sadgmi(self%num_sad))
         do ic = 1, self%num_sad
            Allocate(self%sadgmi(ic)%pArray3D(i1:i2, ju1:j2, k1:k2))
            self%sadgmi(ic)%pArray3D(i1:i2, ju1:j2, k1:k2) = 0.0d0
         end do

         if ((self%sad_opt == 1) .or. (self%sad_opt == 2)) then

            Allocate(self%hno3gas(i1:i2, ju1:j2, k1:k2))
            self%hno3gas = 0.0d0

            Allocate(self%hno3cond(i1:i2, ju1:j2, k1:k2))
            self%hno3cond = 0.0d0

            if (self%sad_opt == 2) then
               Allocate(self%h2ocond(i1:i2, ju1:j2, k1:k2))
               self%h2ocond = 0.0d0

               Allocate(self%h2oback(i1:i2, ju1:j2, k1:k2))
               self%h2oback = 0.0d0

               Allocate(self%vfall(i1:i2, ju1:j2, k1:k2))
               self%vfall = 0.0d0

               Allocate(self%reffice(i1:i2, ju1:j2, k1:k2))
               self%reffice = 0.0d0

               Allocate(self%reffsts(i1:i2, ju1:j2, k1:k2))
               self%reffsts = 0.0d0

               Allocate(self%h2oclim(i1:i2, ju1:j2, k1:k2, self%h2oclim_timpyr))
               self%h2oclim = 0.0d0

               Allocate(self%ch4clim(i1:i2, ju1:j2, k1:k2, self%h2oclim_timpyr))
               self%ch4clim = 0.0d0

               Allocate(self%lbssad(i1:i2, ju1:j2, k1:k2, self%lbssad_timpyr))
               self%lbssad = 0.0d0
            end if
          elseif (self%sad_opt == 3) then
            Allocate(self%lbssad(i1:i2, ju1:j2, k1:k2, self%lbssad_timpyr))
            self%lbssad = 0.0d0
          end if
      end if

! Allocate space, etc., but the initialization of the
! species from the internal state is left to the run method.
! ----------------------------------------------------------

      CALL InitializeSpcConcentration(self%SpeciesConcentration,              &
     &               self%gmiGrid, gmiConfigFile, numSpecies, NMF, NCHEM,     &
     &               loc_proc)

      ALLOCATE (HNO3CONDsad(i1:i2,j1:j2,1:km))
      HNO3CONDsad(i1:i2,j1:j2,1:km) = 0.0
      CALL initDataInStateField(expChem, w_c%grid_esmf, HNO3CONDsad, 'HNO3CONDsad')

      ALLOCATE (HNO3GASsad(i1:i2,j1:j2,1:km))
      HNO3GASsad(i1:i2,j1:j2,1:km) = 0.0
      CALL initDataInStateField(expChem, w_c%grid_esmf, HNO3GASsad,  'HNO3GASsad')

   !------------------------------
   ! Add Fields to the ESMF Bundle
   !------------------------------

   ! Get the declared bundle from the state

   call ESMF_StateGet(expChem, 'gmiSAD' , sadBun,   RC=STATUS)
   VERIFY_(STATUS)

   ! Add tracer to the bundle
   do ib = 1, NSAD
      allocate( var(i1:i2, j1:j2, 1:km), STAT=STATUS)
      VERIFY_(STATUS)

      var(:,:,:)  = 0.0d0
      
      call addTracerToBundle (sadBun, var, w_c%grid_esmf, self%sadName(ib))
   end do

   ! Sanity check

   call ESMF_FieldBundleGet(sadBun, fieldCount=numVars , rc=STATUS)
   VERIFY_(STATUS)
   IF (MAPL_AM_I_ROOT()) THEN
      PRINT*,"  Number of fields in the bundle gmiSAD: ", numVars
   END IF
   _ASSERT(NSAD == numVars,'needs informative message')

    !---------------------------------------------------------------
    ! Create and populate the array that maps GMI species indices to
    ! GEOS-5 species indices
    !---------------------------------------------------------------

    allocate(self%mapSpecies(NSP))
    self%mapSpecies(:) = speciesReg_for_CCM(lchemvar, w_c%reg%vname, NSP, w_c%reg%i_GMI, w_c%reg%j_XX)

  RETURN

  END SUBROUTINE GmiSAD_GridCompInitialize

!-------------------------------------------------------------------------
!NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1, GEOS/DAS!
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  GmiSAD_GridCompRun --- The GMI Driver
!
! !INTERFACE:
!

   SUBROUTINE GmiSAD_GridCompRun ( self, w_c, impChem, expChem, nymd, nhms, &
                                tdt, rc )

! !USES:

   USE GmiTimeControl_mod,            ONLY : Set_curGmiDate, Set_curGmiTime
   USE GmiTimeControl_mod,            ONLY : Set_numTimeSteps, Get_numTimeSteps
   USE GmiTimeControl_mod,            ONLY : Set_gmiSeconds, GetSecondsFromJanuary1
   USE GmiSolar_mod,                  ONLY : CalcCosSolarZenithAngle
   USE GmiSpcConcentrationMethod_mod, ONLY : resetFixedConcentration
   USE GmiUpdateSAD_mod,       ONLY : updateSurfaceAreaDensities

   IMPLICIT none

! !INPUT/OUTPUT PARAMETERS:

   TYPE(GmiSAD_GridComp), INTENT(INOUT) :: self ! Grid Component
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

   CHARACTER(LEN=*), PARAMETER :: IAm = 'GmiSAD_GridCompRun'

!  Input fields from GEOS-5
!  ------------------------
   REAL, POINTER, DIMENSION(:,:,:) :: ple, T, Q, qctot
   REAL, POINTER, DIMENSION(:,:,:) :: reffice, reffsts, fallVel

!  Exports not part of internal state
!  ----------------------------------
   REAL(rPrec), POINTER, DIMENSION(:,:,:) :: HNO3CONDsad
   REAL(rPrec), POINTER, DIMENSION(:,:,:) :: HNO3GASsad

!  Local
!  -----
   INTEGER :: cymd, dymd, emiss_opt, hms
   INTEGER :: i, i1, i2, ic, idehyd_num, im, iXj, iTile(1)
   INTEGER :: j, j1, j2, jm, jTile(1)
   INTEGER :: k, km, kReverse
   INTEGER :: i1_gl, i2_gl, ju1_gl, j2_gl, ilo, gmi_nborder
   INTEGER :: ihi, julo, jhi, ju1,  k1, k2, ilong, ilat
   INTEGER :: loc_proc
   INTEGER :: n
   INTEGER :: STATUS

   INTEGER, PARAMETER :: ToGMI = 1
   INTEGER, PARAMETER :: FromGMI = -1

   REAL :: pi,degToRad,radToDeg,OneOverDt,r

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
!  LOGICAL, PARAMETER :: do_qqjk_reset = .TRUE.
   LOGICAL, PARAMETER :: doThis = .FALSE.
   
! Allocatables.  Use KIND=DBL where GMI expects REAL*8.
! -----------------------------------------------------
   REAL, ALLOCATABLE :: var3d(:,:,:)
   REAL, ALLOCATABLE :: pl(:,:,:)

   REAL(KIND=DBL), ALLOCATABLE :: lonDeg(:,:)
   REAL(KIND=DBL), ALLOCATABLE :: latDeg(:,:)

   REAL(KIND=DBL), ALLOCATABLE :: tropopausePress(:,:)

   REAL(KIND=DBL), ALLOCATABLE :: press3c(:,:,:)
   REAL(KIND=DBL), ALLOCATABLE :: press3e(:,:,:)
   REAL(KIND=DBL), ALLOCATABLE :: kel(:,:,:)

   REAL*8                      :: tdt8

   tdt8 = tdt

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

   gmi_nborder = 0
   i1_gl  = i1
   i2_gl  = i2
   ju1_gl = j1
   j2_gl  = j2
   ju1    = j1
   k1     = 1
   k2     = km
   ilo    = i1  - gmi_nborder
   ihi    = i2  + gmi_nborder
   julo   = ju1 - gmi_nborder
   jhi    = j2  + gmi_nborder
   ilong  = i2 - i1  + 1
   ilat   = j2 - ju1 + 1

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
   ALLOCATE(lonDeg(i1:i2,j1:j2),STAT=STATUS)
   VERIFY_(STATUS)
   ALLOCATE(latDeg(i1:i2,j1:j2),STAT=STATUS)
   VERIFY_(STATUS)

   ALLOCATE(        tropopausePress(i1:i2,j1:j2),STAT=STATUS)
   VERIFY_(STATUS)

   ALLOCATE(                pl(i1:i2,j1:j2,1:km),STAT=STATUS)
   VERIFY_(STATUS)
   ALLOCATE(             var3d(i1:i2,j1:j2,1:km),STAT=STATUS)
   VERIFY_(STATUS)
   ALLOCATE(           press3c(i1:i2,j1:j2,1:km),STAT=STATUS)
   VERIFY_(STATUS)
   ALLOCATE(           press3e(i1:i2,j1:j2,0:km),STAT=STATUS)
   VERIFY_(STATUS)
   ALLOCATE(               kel(i1:i2,j1:j2,1:km),STAT=STATUS)
   VERIFY_(STATUS)

! Geolocation
! -----------
   lonDeg(i1:i2,j1:j2)=self%lonRad(i1:i2,j1:j2)*radToDeg
   latDeg(i1:i2,j1:j2)=self%latRad(i1:i2,j1:j2)*radToDeg

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

! Grab imports and do units conversions
! -------------------------------------
   CALL SatisfyImports(STATUS)
   VERIFY_(STATUS)

! Hand the species concentrations to GMI's bundle
! -----------------------------------------------
   IF (self%gotImportRst) then
      CALL SwapSpeciesBundles(ToGMI, self%SpeciesConcentration%concentration, &
               w_c%qa, Q, self%mapSpecies, lchemvar, self%do_synoz, NSP, &
               STATUS)
      VERIFY_(STATUS)
   END IF

   DEALLOCATE(var3d, STAT=STATUS)
   VERIFY_(STATUS)

   IF (self%gotImportRst) THEN

! Impose fixed concentrations
! ---------------------------
    IF(self%SpeciesConcentration%num_fixed_const > 0) THEN
     CALL resetFixedConcentration(self%SpeciesConcentration, self%gmiClock, self%gmiGrid, NSP)
    END IF

! Update
! ---------------------------
    self%dehydmin = 0.00

    CALL updateSurfaceAreaDensities (tdt8, tropopausePress, press3c,      &
    	       press3e, kel, self%SpeciesConcentration%concentration,	  &
    	       self%ch4clim, self%h2oclim, self%hno3cond, self%hno3gas,   &
    	       self%lbssad, self%sadgmi, self%h2oback, self%h2ocond,	  &
    	       self%reffice, self%reffsts, self%vfall, self%dehydmin,	  &
    	       self%dehyd_opt, self%h2oclim_opt, self%lbssad_opt,	  &
    	       self%sad_opt, IHNO3COND, self%idehyd_num, ICH4,  	  &
    	       IHNO3, IH2O, nymd, self%pr_diag, loc_proc, NSP,  	  &
    	       NSAD, self%lbssad_timpyr, self%h2oclim_timpyr,		  &
    	       ilo, ihi, julo, jhi, i1, i2, ju1, j2, k1, k2, londeg,	  &
    	       latdeg, self%NoPSCZone, self%PSCMaxP, self%chem_mecha)

! Return species concentrations to the chemistry bundle
! -----------------------------------------------------
    CALL SwapSpeciesBundles(FromGMI, self%SpeciesConcentration%concentration, &
    	     w_c%qa, Q, self%mapSpecies, lchemvar, self%do_synoz, NSP,  &
    	     STATUS)
    VERIFY_(STATUS)

! Check for HNO3COND above chosen limit
! -------------------------------------
    DO ic = 1,NSP
     IF(TRIM(lchemvar(ic)) == "HNO3COND") EXIT
    END DO
    r = MAXVAL(self%SpeciesConcentration%concentration(ic)%pArray3D(:,:,:))
    IF(r > self%HNO3Ice_MAX*1.00E-09) THEN
     PRINT *,TRIM(Iam)//": Found HNO3COND above limit: ",r*1.00E+09," ppbv"
     status = 1
     VERIFY_(status)
    END IF

! Check for HCl above chosen limit
! -------------------------------------
    DO ic = 1,NSP
     IF(TRIM(lchemvar(ic)) == "HCl") EXIT
    END DO
    r = MAXVAL(self%SpeciesConcentration%concentration(ic)%pArray3D(:,:,:))
    IF(r > self%HCl_MAX*1.00E-09) THEN
     PRINT *,TRIM(Iam)//": Found HCl above limit: ",r*1.00E+09," ppbv"
     status = 1
     VERIFY_(status)
    END IF

   END IF

! Export states
! -------------
   CALL FillExports(STATUS)
   VERIFY_(STATUS)

   CALL populateBundleSAD (expChem)

! Scratch local work space
! ------------------------
   DEALLOCATE(tropopausePress, STAT=STATUS)
   VERIFY_(STATUS)

   DEALLOCATE(pl, press3c, press3e, kel, STAT=STATUS)
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
!   Stratospheric sulfate surface area
!   Emissions
!
!EOP
!---------------------------------------------------------------------------

  CHARACTER(LEN=255) :: IAm
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

! Sulfate surface area. 12 month climatology, zonal average.
! ----------------------------------------------------------
  self%lbssad(:,:,:,1)=self%lbssad_init_val
  IF(self%lbssad_opt >= 2) THEN
    importName = 'SAD'
    CALL MAPL_GetPointer(impChem, PTR3D, TRIM(importName), RC=STATUS)
    VERIFY_(STATUS)
    self%lbssad(:,:,1:km,1) = PTR3D(:,:,km:1:-1)
    NULLIFY(PTR3D)
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

      ALLOCATE (HNO3CONDsad(i1:i2,j1:j2,1:km))
      HNO3CONDsad(i1:i2,j1:j2,km:1:-1) = self%hno3cond(i1:i2,j1:j2,1:km)
      CALL setDataToStateField(expChem, HNO3CONDsad, 'HNO3CONDsad')

      ALLOCATE (HNO3GASsad(i1:i2,j1:j2,1:km))
      HNO3GASsad(i1:i2,j1:j2,km:1:-1) = self%hno3gas(i1:i2,j1:j2,1:km)
      CALL setDataToStateField(expChem,  HNO3GASsad,  'HNO3GASsad')

      DEALLOCATE (HNO3CONDsad, HNO3GASsad)

! Ice effective radius
! --------------------
   IF(self%gotImportRst .AND. ASSOCIATED(reffice)) &
    reffice(i1:i2,j1:j2,km:1:-1) = self%reffice(i1:i2,j1:j2,1:km)

! STS effective radius
! --------------------
   IF(self%gotImportRst .AND. ASSOCIATED(reffsts)) &
    reffsts(i1:i2,j1:j2,km:1:-1) = self%reffsts(i1:i2,j1:j2,1:km)

! Effective aerosol fall velocity
! -------------------------------
   IF(self%gotImportRst .AND. ASSOCIATED(fallVel)) &
    fallVel(i1:i2,j1:j2,km:1:-1) = self%vfall(i1:i2,j1:j2,1:km)

  RETURN
 END SUBROUTINE FillExports

!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  populateBundleSAD(state)
!
! !INTERFACE:

      subroutine populateBundleSAD(state)
!
      implicit none
!
! !INPUT/OUTPUT PARAMETERS:
      type(ESMF_State), intent(inOut) :: state
!
! !DESCRIPTION:
!
! !LOCAL VARIABLES:
      integer :: STATUS, numVars, ib, rc
      real(rPrec), pointer, dimension(:,:,:) :: ptr3Dreal
      type(ESMF_FieldBundle)                :: sadBun
      character(len=ESMF_MAXSTR), parameter :: IAm = "populateBundleSAD"
!
!EOP
!--------------------------------------------------------------------------------

      ! Get the bundle from the state

      call ESMF_StateGet(state, "gmiSAD", sadBun, rc=STATUS)
      VERIFY_(STATUS)

      call ESMF_FieldBundleGet(sadBun, fieldCount=numVars, rc=STATUS)
      VERIFY_(STATUS)

      ! Verify that the number of fields in the bundle is equal to the number
      ! of SAD variables.

      _ASSERT(numVars == NSAD,'needs informative message')

      allocate(ptr3Dreal(i1:i2, j1:j2, 1:km))

      do ib = 1, numVars
         ptr3Dreal(:,:,:) = self%sadgmi(ib)%pArray3D(:,:,km:1:-1)

         call updateTracerToBundle(sadBun, ptr3Dreal, ib)
      end do

      return

      end subroutine populateBundleSAD

!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  FindPointers
!
! !INTERFACE:

  SUBROUTINE FindPointers(rc)
  
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
   CALL MAPL_GetPointer(impChem,       ple,	'PLE', RC=STATUS)
   VERIFY_(STATUS)
   CALL MAPL_GetPointer(impChem,         Q,       'Q', RC=STATUS)
   VERIFY_(STATUS)
   CALL MAPL_GetPointer(impChem,	 T,	  'T', RC=STATUS)
   VERIFY_(STATUS)
   CALL MAPL_GetPointer(impChem,     qctot,   'QCTOT', RC=STATUS)
   VERIFY_(STATUS)

   CALL MAPL_GetPointer(expChem,   reffice, 'REFFICE', RC=STATUS)
   VERIFY_(STATUS)
   CALL MAPL_GetPointer(expChem,   reffsts, 'REFFSTS', RC=STATUS)
   VERIFY_(STATUS)
   CALL MAPL_GetPointer(expChem,   fallVel,   'VFALL', RC=STATUS)
   VERIFY_(STATUS)

!  Validation
!  ----------
   Validate: IF(self%verbose) THEN
    IF(MAPL_AM_I_ROOT( ))PRINT *,TRIM(IAm),": Input ..."
    i = w_c%reg%j_XX
    CALL pmaxmin('TROPP:', w_c%qa(i)%data3d(:,:,km), qmin, qmax, iXj, 1, 0.01 )
    CALL pmaxmin('Q:', Q, qmin, qmax, iXj, km, 1. )
    CALL pmaxmin('T:', T, qmin, qmax, iXj, km, 1. )
    CALL pmaxmin('PLE (hPa):', ple, qmin, qmax, iXj, km+1, 0.01 )
    CALL pmaxmin('QCTOT:', qctot, qmin, qmax, iXj, km, 1. )
   END IF Validate

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
  
  rc=0
  IAm="SatisfyImports"

! ------------------------------------------------------------------------ 
! Imports to REAL*8, apply units conversion, and reverse vertical stacking
! ------------------------------------------------------------------------

! Singly-layered                                                            GEOS-5 Units       GMI Units
! The most recent valid tropopause pressures are stored in T2M15D(:,:,km)
! -----------------------------------------------------------------------   ------------       -------------
  i = w_c%reg%j_XX
  tropopausePress(i1:i2,j1:j2) = w_c%qa(i)%data3d(i1:i2,j1:j2,km)*Pa2hPa    ! Pa               hPa

! Layer means                                                               GEOS-5 Units       GMI Units
! -----------                                                               ------------       -------------
   press3c(i1:i2,j1:j2,km:1:-1) = pl(i1:i2,j1:j2,1:km)*Pa2hPa               ! Pa  hPa
   kel(i1:i2,j1:j2,km:1:-1) = T(i1:i2,j1:j2,1:km)                           ! K

! Set the condensed water to QCTOT, the total cloud water                    GEOS-5 Units       GMI Units
! -------------------------------------------------------                    ------------       -------------
   self%h2ocond(i1:i2,j1:j2,km:1:-1) = qctot(i1:i2,j1:j2,1:km)*(MWTAIR/MWTH2O) ! kg kg^{-1}    mol mol^{-1}


! Layer edges                                                               GEOS-5 Units       GMI Units
! -----------                                                               ------------       -------------
   press3e(i1:i2,j1:j2,km:0:-1) = ple(i1:i2,j1:j2,0:km)*Pa2hPa               ! Pa hPa

  RETURN
 END SUBROUTINE SatisfyImports

 END SUBROUTINE GmiSAD_GridCompRun

!-------------------------------------------------------------------------
!NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1, GEOS/DAS!
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  GmiSAD_GridCompFinalize
!
! !INTERFACE:
!

   SUBROUTINE GmiSAD_GridCompFinalize ( self, w_c, impChem, expChem, &
                                     nymd, nhms, cdt, rc )

  IMPLICIT none

! !INPUT/OUTPUT PARAMETERS:

   TYPE(GmiSAD_GridComp), INTENT(inout) :: self ! Grid Component

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

   CHARACTER(LEN=*), PARAMETER :: IAm = 'GmiSAD_GridCompFinalize'
   INTEGER :: status
   rc=0

   DEALLOCATE(self%lonRad, self%latRad, STAT=status)
   VERIFY_(status)

   RETURN

 END SUBROUTINE GmiSAD_GridCompFinalize
  
 END MODULE GmiSAD_GCCMod

