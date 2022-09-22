!-------------------------------------------------------------------------
! NASA GSFC - SIVO Code 610.3
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: GmiChemistryMethod_mod
!
! !INTERFACE:
!
#include "MAPL_Generic.h"
!
      module GmiChemistryMethod_mod
!
! !USES:
      use ESMF
      use MAPL
      use GmiGrid_mod             , only : t_gmiGrid, Get_i1, Get_i2
      use GmiGrid_mod             , only : Get_ju1, Get_j2, Get_k1, Get_k2
      use GmiGrid_mod             , only : Get_ilo, Get_ihi, Get_julo, Get_jhi
      use GmiGrid_mod             , only : Get_ilong, Get_ilat, Get_ivert, Get_itloop
      use GmiTimeControl_mod,       only : t_GmiClock, Get_curGmiDate, GmiSplitDateTime
      use GmiSpeciesRegistry_mod,   only : UNKNOWN_SPECIES
      use GmiPrintError_mod,        only : GmiPrintError
      use GmiESMFrcFileReading_mod, only : rcEsmfReadTable
      use GmiArrayBundlePointer_mod, only : t_GmiArrayBundle
      use GmiSavedVariables_mod,     only : t_ChemistrySaved

!
  implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
      private
      public  :: InitializeChemistry, RunChemistry     , FinalizeChemistry
      public  :: initReadChemistry
      public  :: Get_oz_eq_synoz_opt
      public  :: Get_do_chem_grp 
      public  :: Get_t_cloud_ice  , Get_synoz_threshold
      public  :: Get_qqjda        , Get_qqkda 
      public  :: Get_qqjgmi       , Get_qqkgmi
      public  :: Get_yda    
      public  :: Set_qqjda        , Set_qqkda 
      public  :: Set_qqjgmi       , Set_qqkgmi
      public  :: Set_yda      
!
! !PUBLIC DATA MEMBERS:

      public  :: t_Chemistry

# include "GmiParameters.h"
# include "gmi_AerDust_const.h"
# include "setkin_par.h"
# include "setkin_lchem.h"
# include "setkin_mw.h"
# include "setkin_surf_bc.h"

 
  type t_Chemistry
!    private
     type(t_ChemistrySaved) :: savedVars

    real*8              :: dehydmin ! minimum dehyd value (mixing ratio)
    integer             :: ibrono2_num, ich4_num, in2o_num, idehyd_num, ih2_num
    integer             :: ih2o_num, ih2o2_num, ihcl_num, ico_num
    integer             :: ihno3_num, imgas_num, initrogen_num, ioxygen_num
    integer             :: ihno3cond_num, io3_num, isynoz_num
    integer             :: ino_num, iacetone_num, iisoprene_num, ipropene_num
    integer             :: num_nox ! number of Nodoz NOx species
    integer             :: num_noy ! number of Nodoz NOy species
    integer             :: nox_map(MAX_NUM_SMARRAY) ! mapping of NOx spec. # to const spc. #
    integer             :: noy_map(MAX_NUM_SMARRAY) ! mapping of NOy spec. # to const spc. #
    logical             :: do_synoz  ! do Synoz?
    logical             :: do_nodoz  ! do Nodoz?
    integer             :: cloudDroplet_opt
    real*8, pointer     :: cloud_param(:,:,:,:) => null()
    real*8, pointer     :: cloud_tau  (:,:,:) => null()
    real*8, pointer     :: flux_gt    (:,:,:) => null()

    real*8, pointer     :: qqjgmi  (:,:,:,:) => null()
    real*8, pointer     :: qqkgmi  (:,:,:,:) => null()
    real*8, pointer     :: qqjda   (:,:,:,:) => null()
    real*8, pointer     :: qqkda   (:,:,:,:) => null()
    real*8, pointer     :: yda     (:,:,:,:) => null()

    integer             :: chem_opt
    real*8              :: chem_cycle
    integer             :: chem_mask_klo
    integer             :: chem_mask_khi
    real*8              :: synoz_threshold
    real*8              :: t_cloud_ice
!
    integer             :: oz_eq_synoz_opt
!
    logical             :: do_chem_grp
    logical             :: do_smv_reord
    logical             :: do_wetchem
!
!
    integer             :: sad_opt
    integer             :: h2oclim_opt
!
    character (len=MAX_LENGTH_VAR_NAME)  :: qj_var_name
    integer             :: phot_opt
    logical             :: do_clear_sky
!
    logical             :: do_full_chem
    real*8 , pointer    :: mw (:) => null()
    integer             :: num_sad
    integer             :: num_qjs
    integer             :: num_qjo
    integer             :: num_qks
    integer             :: num_chem
    integer             :: num_active
    integer             :: num_ks_sbc
    integer             :: num_spc_sbc
    integer             :: num_molefrac
    integer, pointer    :: surf_bc_map(:) => null()
    character(len=MAX_LENGTH_SPECIES_NAME), pointer :: const_labels(:) => null() ! constituent string labels
!
  end type t_Chemistry

  INTEGER :: STATUS, RC
!
! !DESCRIPTION:
! Initialization and run methods for the Chemistry component.
!
! !AUTHOR:
! Jules Kouatchou, NASA/GSFC, Jules.Kouatchou-1@nasa.gov
!
! !REVISION HISTORY:
!
!EOP
!-------------------------------------------------------------------------
  CONTAINS
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: readChemistryResourceFile
!
! !INTERFACE:
!
      subroutine readChemistryResourceFile (self, config, loc_proc,   &
     &               rootProc, numSpecies, i1, i2, ju1, j2, k1, k2, pr_diag)
!
! !USES:
      use GmiCheckNamelistFile_mod, only : CheckNamelistOptionRange
      use GmiStringManipulation_mod, only : constructListNames
!
! !INPUT PARAMETERS:
      logical,         intent(in) :: rootProc, pr_diag
      integer,         intent(in) :: loc_proc, numSpecies
      integer,         intent(in) :: i1, i2, ju1, j2, k1, k2
!
! !INPUT/OUTPUT PARAMETERS:
      type (t_Chemistry), intent(inOut) :: self
      type (ESMF_Config), intent(inOut) :: config
!
! !DESCRIPTION:
! Reads in Chemistry related variables from the resource file.
!
! !LOCAL VARIABLES:
      integer :: ic
      real    :: tempR4
      real*8  :: hugeReal
      character (len=MAX_LENGTH_SPECIES_NAME), pointer :: tempListNames(:)
      character (len=MAX_STRING_LENGTH      ) :: forcedBcSpeciesNames
      character(len=ESMF_MAXSTR) :: IAm, err_msg
!
!EOP
!------------------------------------------------------------------------------
!BOC
      IAm = "readChemistryResourceFile"

      if (pr_diag) Write(6,*) IAm, 'called by ', loc_proc

      allocate(tempListNames(numSpecies), STAT=STATUS)
      VERIFY_(STATUS)

      !################################
      ! Begin reading the resource file
      !################################

!     -------------------------------------
!     chem_opt
!       0:  no chemistry
!       1:  call Radon/Lead chemistry
!       2:  call Smvgear2
!       3:  call simple loss
!       4:  call forcing boundary condition
!       5:  call Synoz tracer
!       6:  call Beryllium chemistry
!       7:  call Quadchem
!       8:  call Sulfur chemistry
!     -------------------------------------

      call ESMF_ConfigGetAttribute(config, self%chem_opt, &
     &                label   = "chem_opt:", &
     &                default = 2, rc=STATUS )
      VERIFY_(STATUS)

!     -------------------------------------------------------------
!     chem_cycle:  number of time steps to cycle chemistry calls on
!       < 1.0:  chemistry will subcycle
!         1.0:  chemistry called each time step
!     -------------------------------------------------------------

      call ESMF_ConfigGetAttribute(config, self%chem_cycle, &
     &                label   = "chem_cycle:", &
     &                default = 1.0d0, rc=STATUS )
      VERIFY_(STATUS)

!     -----------------------------------------------------
!     chem_mask_klo, chem_mask_khi:
!       chemistry turned off where k is outside of range of
!       [chem_mask_klo, chem_mask_khi]
!     -----------------------------------------------------

      call ESMF_ConfigGetAttribute(config, self%chem_mask_klo, &
     &                label   = "chem_mask_klo:", &
     &                default = k1, rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, self%chem_mask_khi, &
     &                label   = "chem_mask_khi:", &
     &                default = k2, rc=STATUS )
      VERIFY_(STATUS)

!     -------------------------------------------------------------------
!     synoz_threshold:  chemistry turned off where synoz > this threshold
!     -------------------------------------------------------------------

      hugeReal = Huge (hugeReal)

      call ESMF_ConfigGetAttribute(config, self%synoz_threshold, &
     &                label   = "synoz_threshold:", &
     &                default = hugeReal, rc=STATUS )
      VERIFY_(STATUS)

!    cloudDroplet_opt = 1: Boucher and LohMan    Correlation (default)
!                     = 2: Nenes and Seinfeld    Parameterization
!                     = 3: Abdul-Razzak and Ghan Parameterization
!                     = 4: Segal amd Khain       Correllation
!   The variable is only use when the GT cloud module is employed.

      call ESMF_ConfigGetAttribute(config, self%cloudDroplet_opt, &
     &                label   = "cloudDroplet_opt:", &
     &                default = 1, rc=STATUS )
      VERIFY_(STATUS)
      
!     -----------------------------------------------------
!     oz_eq_synoz_opt 
!       0:  do not use stratospheric loss in gmi_step.F
!       1:  use stratospheric loss in gmi_step.F
!     -----------------------------------------------------

      call ESMF_ConfigGetAttribute(config, self%oz_eq_synoz_opt, &
     &                label   = "oz_eq_synoz_opt:", &
     &                default = 1, rc=STATUS )
      VERIFY_(STATUS)
      
!     ------------------------------------------------------------
!     t_cloud_ice:  temperature for ice formation in clouds (degK)
!     ------------------------------------------------------------

      call ESMF_ConfigGetAttribute(config, self%t_cloud_ice, &
     &                label   = "t_cloud_ice:", &
     &                default = 263.0d0, rc=STATUS )
      VERIFY_(STATUS) 

      call ESMF_ConfigGetAttribute(config, value=self%do_chem_grp,  label="do_chem_grp:",  &
     &                       default=.false., __RC__ )
     
      call ESMF_ConfigGetAttribute(config, value=self%do_smv_reord, label="do_smv_reord:", &
     &                       default=.false., __RC__ )

      call ESMF_ConfigGetAttribute(config, value=self%do_wetchem,   label="do_wetchem:",   &
     &                       default=.false., __RC__ )

!     ---------------------------
!     Surface Area Density (SAD):
!     ---------------------------

!     ----------------------------------------------------
!     sad_opt
!       0:  do not allocate or process SAD array
!       1:  allocate, but zero out SAD array
!       2:  call Considine code (i.e., Condense)
!       3:  read SAD array from a file of monthly averages
!     ----------------------------------------------------

      call ESMF_ConfigGetAttribute(config, self%sad_opt, &
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

      call ESMF_ConfigGetAttribute(config, self%h2oclim_opt, &
     &                label   = "h2oclim_opt:", &
     &                default = 2, rc=STATUS )
      VERIFY_(STATUS)


!     =========       
!     nlGmiPhotolysis 
!     =========

!!     -----------------------------------------------------
!!     phot_opt
!!       0:  no photolysis
!!       1:  set all qj values to qj_init_val 
!!       2:  read in qj values 
!!       3:  use fastj routine (for fastJX65 or CloudJ)
!!           This option should be combined with fastj_opt.
!!       4:  lookup table for qj (Kawa style)
!!       5:  lookup table for qj (Kawa style) +
!!           use ozone climatology for column ozone calc.
!!       6:  calculate from table and Gmimod data (Quadchem)
!!       7:  read in qj values (2-D, 12 months)
!!       8:  use fast-JX routine (troposphere/stratosphere)
!!     -----------------------------------------------------
!
      call ESMF_ConfigGetAttribute(config, self%phot_opt, &
     &                label   = "phot_opt:", &
     &                default = 1, rc=STATUS )
      VERIFY_(STATUS)
!
!     ---------
!     qj / qqj:
!     ---------

      allocate(self%mw(numSpecies), STAT=STATUS)
      VERIFY_(STATUS)
      self%mw(1:numSpecies) = mw_data(1:numSpecies)

      call rcEsmfReadTable(config, self%mw, "mw::", rc=STATUS)

      allocate(self%const_labels(numSpecies), STAT=STATUS)
      VERIFY_(STATUS)
      
      self%const_labels(1:numSpecies) = lchemvar(1:numSpecies)
     
      call rcEsmfReadTable(config, self%const_labels, &
     &                     "const_labels::", rc=STATUS)
      
    ! ---------------------------------------------------------------
    ! Check option ranges.  Note that as new options are added, these
    ! range checks will have to be modified.
    ! ---------------------------------------------------------------
      
      call CheckNamelistOptionRange ('chem_opt', self%chem_opt, 0, 8)
      call CheckNamelistOptionRange ('h2oclim_opt', self%h2oclim_opt, 1, 3)
      call CheckNamelistOptionRange ('phot_opt', self%phot_opt, 0, 3)
      call CheckNamelistOptionRange ('sad_opt', self%sad_opt, 0, 3)
      call CheckNamelistOptionRange ('oz_eq_synoz_opt', self%oz_eq_synoz_opt, 0, 1)
      call CheckNamelistOptionRange ('cloudDroplet_opt', self%cloudDroplet_opt, 1, 4)
!     
      if ((self%chem_opt == 2) .or. (self%chem_opt == 7) .or. &
     &    (self%chem_opt == 8)) then
          self%do_full_chem = .true.
      
          self%mw          (1:numSpecies) = mw_data (1:numSpecies)
          self%const_labels(1:numSpecies) = lchemvar(1:numSpecies)
      else 
          self%do_full_chem = .false.

         if ((self%mw(1) == 0.0d0) .and. (numSpecies == NSP)) then
            self%mw(1:numSpecies) = mw_data(1:numSpecies)
         end if            
      end if               
      
      if (self%do_full_chem) then
         self%num_active    = NACT
         self%num_chem      = NCHEM
         self%num_molefrac  = NMF
         self%num_qjs       = NUM_J
         self%num_qjo       = NUM_J
         self%num_qks       = NUM_K

         self%io3_num       = IO3
         self%ibrono2_num   = IBRONO2
         self%ich4_num      = ICH4
         self%in2o_num      = IN2O
         self%ih2_num       = IH2
         self%ih2o_num      = IH2O
         self%ih2o2_num     = IH2O2
         self%ihcl_num      = IHCL
         self%ihno3_num     = IHNO3
         self%imgas_num     = IMGAS
         self%initrogen_num = INITROGEN
         self%ioxygen_num   = IOXYGEN
         self%ico_num       = ICO
         self%ino_num       = INO
         self%ihno3cond_num = IHNO3COND    ! JK 12June2010

         self%isynoz_num    = ISYNOZ       ! JK 12June2010

         self%iacetone_num  = IACET
         self%ico_num       = ICO
         self%iisoprene_num = IC5H8
         self%ipropene_num  = IC3H6

         !if (pr_qj_o3_o1d)    self%num_qjo = self%num_qjo + 1
         !if (pr_qj_opt_depth) self%num_qjo = self%num_qjo + 1

         self%num_ks_sbc    = K_SBC

         if (NUM_SBC > 0) then
            allocate(self%surf_bc_map(1:NUM_SBC), STAT=STATUS)
            VERIFY_(STATUS)
            self%surf_bc_map(:)         = 0
            self%surf_bc_map(1:NUM_SBC) = sbc_map(1:NUM_SBC)

            if (self%surf_bc_map(1) /= 0) then
              self%num_spc_sbc = NUM_SBC
            else
              self%num_spc_sbc = 0
            end if
         end if

        !if (self%do_chem_grp) then
        !   call setupChemicalGroup (i1, i2, ju1, j2, k1, k2)
        !end if

!     ====
      else
!     ====

         self%num_chem      = numSpecies
         self%num_molefrac  = numSpecies
         self%num_qjs       = 0
         self%num_qks       = 0

         self%num_ks_sbc    = 0
         self%num_spc_sbc   = 0

         self%ih2o2_num     = 0
         self%ih2_num       = 0
         self%ih2o_num      = 0
         self%ihno3_num     = 0
         self%idehyd_num    = 0
         self%imgas_num     = 0
         self%ich4_num      = 0
         self%ihno3cond_num = 0
         self%ico_num       = 0
         self%ino_num       = 0
         self%initrogen_num = 0
         self%ioxygen_num   = 0
         self%ihcl_num      = 0
         self%ibrono2_num   = 0
         self%in2o_num      = 0

         call ESMF_ConfigGetAttribute(config, self%io3_num, &
     &                label   = "io3_num:", &
     &                default = 0, rc=STATUS )
         VERIFY_(STATUS)

         call ESMF_ConfigGetAttribute(config, self%iacetone_num, &
     &                label   = "iacetone_num:", &
     &                default = 0, rc=STATUS )
         VERIFY_(STATUS)

         call ESMF_ConfigGetAttribute(config, self%ipropene_num, &
     &                label   = "ipropene_num:", &
     &                default = 0, rc=STATUS )
         VERIFY_(STATUS)

         call ESMF_ConfigGetAttribute(config, self%iisoprene_num, &
     &                label   = "iisoprene_num:", &
     &                default = 0, rc=STATUS )
         VERIFY_(STATUS)
      end if

      self%do_synoz   = .false.
      self%do_nodoz   = .false.

      self%isynoz_num = 0

      self%num_nox    = 0
      self%num_noy    = 0

      self%nox_map(:) = 0
      self%noy_map(:) = 0

      !if (self%do_full_chem .or. (self%chem_opt == 5)) then
      !   call Setup_Syn_Spc (self, numSpecies)
      !end if
      
      if (self%sad_opt /= 0) then
         self%num_sad = NSAD
      else
         self%num_sad = 0
      end if
         
      return
         
      end subroutine readChemistryResourceFile
!EOC     
!-------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: InitializeChemistry
!
! !INTERFACE:
!
      subroutine InitializeChemistry (self, gmiGrid, config, loc_proc,&
     &                     numSpecies, pr_diag, pr_qqjk, do_qqjk_inchem,&
                           pr_smv2, rootProc, tdt4)
!
! !USES:
      use GmiGrid_mod     , only : t_gmiGrid
      use GmiGrid_mod     , only : Get_i1, Get_i2, Get_ju1, Get_j2, Get_k1, Get_k2
      use ReadForcedBC_mod, only : readForcedBcData
!
!
! !INPUT PARAMETERS:
      type (t_gmiGrid),  intent(in) :: gmiGrid  
      integer          , intent(in) :: numSpecies
      integer          , intent(in) :: loc_proc
      logical          , intent(in) :: pr_diag, do_qqjk_inchem, pr_qqjk
      logical          , intent(in) :: pr_smv2
      logical          , intent(in) :: rootProc
      real             , intent(in) :: tdt4
!
! !INPUT/OUTPUT PARAMETERS:
      type (ESMF_Config), intent(inOut) :: config
      type (t_Chemistry), intent(inOut) :: self

!
! !DESCRIPTION:
!  Initializes the Chemistry operator.
!
! !LOCAL VARIABLES:
      integer :: ncid_rra, cnt1d(1), strt1d(1)
      integer  :: i1, i2, ju1, j2, k1, k2, ic, itloop
      character (len=MAX_LENGTH_FILE_NAME) :: smv_filnam
      real*8 :: chemintv
!
!EOP
!------------------------------------------------------------------------------
!BOC
      call Get_i1    (gmiGrid, i1)
      call Get_i2    (gmiGrid, i2)
      call Get_ju1   (gmiGrid, ju1)
      call Get_j2    (gmiGrid, j2)
      call Get_k1    (gmiGrid, k1)
      call Get_k2    (gmiGrid, k2)
      call Get_itloop(gmiGrid, itloop)

      !=======================
      ! Read the resource file
      !=======================

      call readChemistryResourceFile (self, config, loc_proc, rootProc, &
     &         numSpecies, i1, i2, ju1, j2, k1, k2, pr_diag)


      if (self%chem_opt /= 0) then
         if (self%do_full_chem) then

            if (pr_qqjk) then
               if (do_qqjk_inchem) then
                  call Allocate_qqjda(self, i1,i2, ju1,j2, k1, k2)
                  call Allocate_qqkda(self, i1,i2, ju1,j2, k1, k2)
                  call Allocate_yda  (self, i1,i2, ju1,j2, k1, k2)
               end if
            end if

         end if
      end if

      if (pr_smv2) then
         smv_filnam = 'smv2chem.asc'
      end if

      chemintv = tdt4 * self%chem_cycle

      call Do_Smv2_Init (self%savedVars, smv_filnam, self%do_smv_reord, &
                  pr_diag, pr_smv2, loc_proc, itloop, chemintv)

      return

      end subroutine InitializeChemistry
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: initReadChemistry
!
! !INTERFACE:
!
      subroutine initReadChemistry(self, gmiGrid, loc_proc)
!

!
! !INPUT PARAMETERS:
      type  (t_gmiGrid  ), intent(in) :: gmiGrid
      integer            , intent(in) :: loc_proc
!
! !INPUT/OUTPUT PARAMETERS:
      type  (t_Chemistry), intent(inOut) :: self
!
! !DESCRIPTION:
! This routine reads in chemistry related files that contain geographical
! related data. The read is done once during initialization procedures
! for the Chemistry component.
!
! !LOCAL VARIABLES:
      integer :: ncid_rra, cnt1d(1), strt1d(1)
      integer  :: i1, i2, ju1, j2, k1, k2
!
!EOP
!-------------------------------------------------------------------------
!BOC
!
      call Get_i1    (gmiGrid, i1)
      call Get_i2    (gmiGrid, i2)
      call Get_ju1   (gmiGrid, ju1)
      call Get_j2    (gmiGrid, j2)
      call Get_k1    (gmiGrid, k1)
      call Get_k2    (gmiGrid, k2)

      return

      end subroutine initReadChemistry
!
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: runChemistry
!
! !INTERFACE:
!
      subroutine runChemistry (self, SpeciesConcentration, gmiClock, gmiGrid,   &
     &              press3c, press3e, gridBoxHeight, mcor, mass, kel, humidity, &
     &              pctm2, loc_proc, num_species, do_qqjk_reset, HNO3CONDsad,   &
     &              HNO3GASsad, gmiQK, gmiQQK, gmiQJ, gmiQQJ, surfEmissForChem, &
     &              pr_diag, do_ftiming, do_qqjk_inchem, pr_qqjk,               &
     &              do_semiss_inchem, pr_smv2, pr_nc_period,         &
     &              rootProc, metdata_name_org, metdata_name_model, tdt4)

! !USES:
      use GmiUpdateChemistry_mod, only : updateChemistry
      use GmiSpcConcentrationMethod_mod, only : t_SpeciesConcentration
      use GmiSpcConcentrationMethod_mod, only : Get_concentration, Set_concentration
      use GmiTimeControl_mod   , only : GmiSplitDateTime    , GetSecondsFromJanuary1
      use GmiTimeControl_mod   , only : t_GmiClock      , Get_begGmiDate
      use GmiTimeControl_mod   , only : Get_curGmiDate  , Get_curGmiTime
      use GmiTimeControl_mod   , only : Get_gmiSeconds  , Get_numTimeSteps
!
#     include "gmi_sad_constants.h"
!
! !INPUT PARAMETERS:
      real*8 ,          intent(in) :: pr_nc_period
      logical,          intent(in) :: pr_diag, rootProc, do_semiss_inchem
      logical,          intent(in) :: do_ftiming, do_qqjk_inchem, pr_qqjk
      logical,          intent(in) :: pr_smv2
      integer,          intent(in) :: loc_proc
      integer,          intent(in) :: num_species
      real*8 ,          intent(in) :: humidity(:,:,:), pctm2(:,:)
      real*8 ,          intent(in) :: kel(:,:,:), mcor(:,:), mass(:,:,:)
      real*8 ,          intent(in) :: press3c(:,:,:), press3e(:,:,:)
      real*8 ,          intent(in) :: gridBoxHeight(:,:,:)
      type(t_gmiGrid ), intent(in) :: gmiGrid 
      type(t_GmiClock), intent(in) :: gmiClock
      real*8 ,          intent(in) :: surfEmissForChem(:,:,:)
      real*8 ,          intent(in) :: HNO3CONDsad(:,:,:)
      real*8 ,          intent(in) :: HNO3GASsad(:,:,:)
      type (t_GmiArrayBundle), intent(in) :: gmiQJ(:)
      type (t_GmiArrayBundle), intent(in) :: gmiQK(:)
      type (t_GmiArrayBundle), intent(inout) :: gmiQQJ(:)
      type (t_GmiArrayBundle), intent(inout) :: gmiQQK(:)
      character (len=MAX_LENGTH_MET_NAME), intent(in) :: metdata_name_org
      character (len=MAX_LENGTH_MET_NAME), intent(in) :: metdata_name_model
      real,             intent(in) :: tdt4
!
! !INPUT/OUTPUT PARAMETERS:
      logical,                      intent(inOut) :: do_qqjk_reset
      type(t_Chemistry),            intent(inOut) :: self
      type(t_SpeciesConcentration), intent(inOut) :: SpeciesConcentration
!
! !DESCRIPTION:
! Run method for the Chemistry component.
!
! !LOCAL VARIABLES:
!
      integer :: ic

      integer       :: ydummy, thisDay, thisMonth, thisDate, ddummy
      integer       :: nymd, nhms, num_time_steps, start_ymd
      real*8        :: tdt8, gmi_sec
      integer       :: ilo, ihi, julo, jhi
      integer       :: ilong, ilat, ivert, itloop
      integer       :: i1, i2, ju1, j2, k1, k2
      type (t_GmiArrayBundle), pointer :: concentration(:)

      LOGICAL, PARAMETER :: doThis=.FALSE.
!EOP
!-------------------------------------------------------------------------
!EOC
      tdt8 = tdt4

      self%dehydmin = 0.00

      ! Get the GMI grid information
      call Get_i1    (gmiGrid, i1   )
      call Get_i2    (gmiGrid, i2   )
      call Get_ju1   (gmiGrid, ju1  )
      call Get_j2    (gmiGrid, j2   )
      call Get_k1    (gmiGrid, k1   )
      call Get_k2    (gmiGrid, k2   )
      call Get_ilo   (gmiGrid, ilo  )
      call Get_ihi   (gmiGrid, ihi  )
      call Get_julo  (gmiGrid, julo )
      call Get_jhi   (gmiGrid, jhi  )
      call Get_ilong (gmiGrid, ilong)
      call Get_ilat  (gmiGrid, ilat )
      call Get_ivert (gmiGrid, ivert)
      call Get_itloop(gmiGrid, itloop)

      ! Obtain model time information

      call Get_begGmiDate  (gmiClock, start_ymd     )
      call Get_curGmiDate  (gmiClock, nymd	    )
      call Get_curGmiTime  (gmiClock, nhms	    )
      call Get_numTimeSteps(gmiClock, num_time_steps)
      call Get_gmiSeconds  (gmiClock, gmi_sec	    )

       call Get_concentration(SpeciesConcentration, concentration)


      ! Call the Chemistry control routine
      call updateChemistry (self%savedVars, rootProc, do_ftiming, &
     &           TRIM(metdata_name_org), TRIM(metdata_name_model),             &
     &           do_qqjk_inchem, do_qqjk_reset, pr_qqjk,  surfEmissForChem,    &
     &           press3c, press3e, pr_smv2, pr_nc_period, mass, concentration, &
     &           gmiQJ, gmiQK, kel, humidity, pctm2, gmiQQJ, gmiQQK,           &
     &           self%yda, self%qqkda, self%qqjda,                             &
     &           HNO3GASsad, HNO3CONDsad, self%h2oclim_opt, self%chem_opt,     &
     &           self%sad_opt, self%phot_opt, self%do_smv_reord, self%do_synoz,&
     &           do_semiss_inchem, self%do_wetchem, nymd, nhms,                &
     &           gmi_sec, tdt8, pr_diag, loc_proc, self%synoz_threshold,       &
     &           self%chem_cycle, self%chem_mask_klo, self%chem_mask_khi,      &
     &           self%ih2_num, self%ih2o_num, self%ihno3_num, self%ich4_num,   &
     &           self%imgas_num, self%initrogen_num, self%ioxygen_num,         &
     &           self%isynoz_num, num_species, self%num_qks, self%num_qjs,     &
     &           self%num_qjo, self%num_sad, self%num_molefrac, self%num_chem, &
     &           self%num_active, ilong, ilat, ivert, itloop, ilo, ihi, julo,  &
     &           jhi, i1, i2, ju1, j2, k1, k2)

      call Set_concentration(SpeciesConcentration, concentration)

      return

      end subroutine RunChemistry
!EOC
!-------------------------------------------------------------------------

  subroutine FinalizeChemistry (self)

  type (t_Chemistry), intent(inout) :: self

  PRINT*,'  Finalize Chemistry'
  call finalizeJsparse(self%savedVars)
  call finalizeKsparse(self%savedVars)
  call finalizeSMV2solver(self%savedVars)
  call finalizePhysproc(self%savedVars)
!!$  call finalizeSMV2diag(self%savedVars)
  return

  end subroutine FinalizeChemistry
!EOC
!-------------------------------------------------------------------------
      subroutine setupChemicalGroup (i1, i2, ju1, j2, k1, k2)
      
#     include "setkin_group.h"

      integer :: num_cgrp
      integer :: num_cgrp_elem  (MAX_NUM_CGRP)
      integer :: cgrp_elem_map  (MAX_NUM_SMARRAY, MAX_NUM_CGRP)
      real*8  :: cgrp_fac  (MAX_NUM_SMARRAY, MAX_NUM_CGRP)
      
      integer, intent(in) :: i1, i2, ju1, j2, k1, k2
      
!     ----------------------
!     Variable declarations.
!     ----------------------

      integer :: ii, jj
      

!     ----------------
!     Begin Execution.
!     ----------------

      num_cgrp = NUMGRP
      
      num_cgrp_elem(:) = 0
      

      do jj = 1, num_cgrp
      
        IILOOP1: do ii = 1, MAXGRP_ELEM
        
          if (sgrp_elem_map(ii,jj) /= 0) then
          
            num_cgrp_elem(jj)    = num_cgrp_elem(jj) + 1
            
            cgrp_elem_map(ii,jj) = sgrp_elem_map(ii,jj)

            cgrp_fac(ii,jj) = sgrp_fac(ii,jj)

          else

            exit IILOOP1

          end if

        end do IILOOP1

      end do

      return

      end subroutine setupChemicalGroup

!-----------------------------------------------------------------------------
!   This routine sets up the synthetic species info.
!-----------------------------------------------------------------------------

      subroutine Setup_Syn_Spc (self, numSpecies)

#     include "setkin_synspc.h"

      integer, intent(in) :: numSpecies
      type (t_Chemistry), intent(inOut) :: self

      integer :: ii

!      self%do_synoz   = .false.
!      self%do_nodoz   = .false.
!
!      self%isynoz_num = 0
!
!      self%num_nox    = 0
!      self%num_noy    = 0
!
!      self%nox_map(:) = 0
!      self%noy_map(:) = 0


!     ==================
      if (self%chem_opt == 2) then
!     ==================

        if (USE_SYNOZ) then

          self%do_synoz   = .true.

          self%isynoz_num = ISYNOZ

          if (USE_NODOZ) then

            self%do_nodoz = .true.

            NOXLOOP: do ii = 1, MAXNODOZ_ELEM

              if (NOX_ELEM_MAP(ii) /= 0) then

                self%num_nox = self%num_nox + 1

                self%nox_map(ii) = NOX_ELEM_MAP(ii)

              else

                exit NOXLOOP

              end if

            end do NOXLOOP

            NOYLOOP: do ii = 1, MAXNODOZ_ELEM

              if (NOY_ELEM_MAP(ii) /= 0) then

                self%num_noy = self%num_noy + 1

                self%noy_map(ii) = NOY_ELEM_MAP(ii)

              else

                exit NOYLOOP

              end if

            end do NOYLOOP

          end if

        end if

!     =======================
      else if (self%chem_opt == 5) then
!     =======================

        self%do_synoz   = .true.

        self%isynoz_num = 1

        if (numSpecies == 2) then

          self%do_nodoz = .true.

          self%ihno3_num = 2

        end if

!     ======
      end if
!     ======

      return

      end subroutine Setup_Syn_Spc



!-----------------------------------------------------------------------------
! Does error checking to make sure that resource file variables are properly
! set.
!-----------------------------------------------------------------------------

      subroutine rcCheckChemistrySetting(self, pr_sad, pr_qj, pr_qk, pr_qqjk, &
     &       do_qqjk_inchem, do_mean, do_semiss_inchem, numSpecies,&
     &        met_opt)
!
#     include "smv2chem_par.h"
!
      logical, intent(in) :: pr_sad, pr_qj, pr_qk, pr_qqjk
      logical, intent(in) :: do_qqjk_inchem, do_mean
      logical          , intent(in) :: do_semiss_inchem
      integer          , intent(in) :: met_opt
      integer          , intent(in) :: numSpecies
      type(t_Chemistry), intent(in) :: self
!
      character(len=90) :: err_msg
      integer           :: rc, icycles, ic
      real*8            :: rcycles
!
      !-------------------------------------------------------------------
      ! molecular weights must always be supplied for all selected species
      !-------------------------------------------------------------------

      do ic = 1, numSpecies
         if (self%mw(ic) == 0.0d0) then
            err_msg = 'mw problem in rc File.'
            call GmiPrintError (err_msg, .true., 2, 1, ic, 1, self%mw(ic), 0.0d0)
         end if
      end do

!
      if (pr_sad .and. ((self%sad_opt == 0) .or. (self%num_sad == 0))) then
        err_msg = 'pr_sad/sad_opt/num_sad problem.'
        call GmiPrintError  &
     &    (err_msg, .true., 2, self%sad_opt, self%num_sad, 0, 0.0d0, 0.0d0)
      end if

      if (pr_qj .and. ((self%phot_opt == 0) .or. (self%num_qjs == 0))) then
        err_msg = 'pr_qj/phot_opt/num_qj problem.'
        call GmiPrintError  &
     &    (err_msg, .true., 2, self%phot_opt, self%num_qjs, 0, 0.0d0, 0.0d0)
      end if

      if (pr_qk .and. (self%num_qks == 0)) then
        err_msg = 'pr_qk/num_qk problem.'
        call GmiPrintError  &
     &    (err_msg, .true., 1, self%num_qks, 0, 0, 0.0d0, 0.0d0)
      end if

      if (pr_qqjk) then
        if ((self%num_qjs == 0) .or. (self%num_qks == 0)) then
          err_msg = 'pr_qqjk/num_qj/k problem.'
          call GmiPrintError  &
     &      (err_msg, .true., 2, self%num_qjs, self%num_qks, 0, 0.0d0, 0.0d0)
        end if

        if (self%phot_opt == 0) then
          err_msg = 'pr_qqjk/phot_opt problem.'
          call GmiPrintError  &
     &      (err_msg, .true., 1, self%phot_opt, 0, 0, 0.0d0, 0.0d0)
        end if

        if (.not. self%do_full_chem) then
          err_msg = 'pr_qqjk/chem_opt problem.'
          call GmiPrintError  &
     &      (err_msg, .true., 1, self%chem_opt, 0, 0, 0.0d0, 0.0d0)
        end if

        if (do_qqjk_inchem .and. do_mean) then
          err_msg = 'do_qqjk_inchem/do_mean problem.'
          call GmiPrintError (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
        end if
      end if

!     ==================
      if (self%chem_opt == 1) then
!     ==================
         if (numSpecies /= 2) then
            err_msg = 'The Radon/Lead chemistry requires two species.'
            call GmiPrintError (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
         end if
      end if
            
!     ==================
      if (self%chem_opt == 2) then
!     ==================
      
        if (numSpecies < IGAS) then
          err_msg = 'num_species/IGAS problem.'
          call GmiPrintError (err_msg, .true., 2, numSpecies, IGAS, 0, 0.0d0, 0.0d0)
        end if
      
        if (self%num_qjs /= IPHOT) then
          err_msg = 'num_qjs/IPHOT problem.'
          call GmiPrintError  &
     &      (err_msg, .true., 2, self%num_qjs, IPHOT, 0, 0.0d0, 0.0d0)
        end if
      
        if (self%num_qks /= ITHERM) then
          err_msg = 'num_qks/ITHERM problem in.' 
          call GmiPrintError  &
     &      (err_msg, .true., 2, self%num_qks, ITHERM, 0, 0.0d0, 0.0d0)
        end if
      
      end if
      
!     ==================
      if (self%chem_opt == 6) then
!     ==================
     
        if (numSpecies /= 2) then
          err_msg = 'The Beryllium chemistry requires two species.'
          call GmiPrintError  &
     &      (err_msg, .true., 1, numSpecies, 0, 0, 0.0d0, 0.0d0)
        end if 
        
      end if

!     =================
      if (self%do_full_chem) then
!     =================

        if (self%phot_opt == 0) then
          err_msg = 'chem_opt/phot_opt problem in the rc File.'
          call GmiPrintError  &
     &      (err_msg, .true., 2, self%chem_opt, self%phot_opt, 0, 0.0d0, 0.0d0)
        end if

        if ((self%sad_opt == 0) .and. (met_opt /= 3)) then
          err_msg = 'sad_opt/met_opt problem in the rc File.'
          call GmiPrintError  &
     &      (err_msg, .true., 2, self%sad_opt, met_opt, 0, 0.0d0, 0.0d0)
        end if

        if ((self%sad_opt == 1) .or. (self%sad_opt == 2)) then

          if ((self%ich4_num == 0) .or. (self%ich4_num /= ICH4)) then
            err_msg = 'chem_opt/ich4_num problem in the rc File.'
            call GmiPrintError  &
     &        (err_msg, .true., 2, self%chem_opt, self%ich4_num, 0, 0.0d0, 0.0d0)
          end if

          if ((self%idehyd_num == 0) .or. (self%idehyd_num /= IDEHYD)) then
            err_msg = 'chem_opt/idehyd_num problem in the rc File.'
            call GmiPrintError  &
     &       (err_msg, .true., 2, self%chem_opt, self%idehyd_num, 0, 0.0d0, 0.0d0)
          end if

          if ((self%ih2o_num == 0) .or. (self%ih2o_num /= IH2O)) then
            err_msg = 'chem_opt/ih2o_num problem in the rc File.'
            call GmiPrintError  &
     &        (err_msg, .true., 2, self%chem_opt, self%ih2o_num, 0, 0.0d0, 0.0d0)
          end if

        end if

        if ((self%imgas_num == 0) .or. (self%imgas_num /= IMGAS)) then
          err_msg = 'chem_opt/imgas_num problem in the rc File.'
          call GmiPrintError  &
     &      (err_msg, .true., 2, self%chem_opt, self%imgas_num, 0, 0.0d0, 0.0d0)
        end if

        if (self%do_synoz) then
          if (self%isynoz_num == 0) then
            err_msg = 'do_synoz/isynoz_num problem in the rc File.'
            call GmiPrintError  &
     &        (err_msg, .true., 1, self%isynoz_num, 0, 0, 0.0d0, 0.0d0)
          end if
        end if

        if (self%do_nodoz) then
          if (.not. self%do_synoz) then
             err_msg = 'do_nodoz/do_synoz problem in the rc File.'
          call GmiPrintError  &
     &      (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
          end if

          if ((self%num_nox == 0) .or. (self%num_noy == 0)) then
            err_msg = 'num_nox/num_noy problem in the rc File.'
            call GmiPrintError  &
     &        (err_msg, .true., 2, self%num_nox, self%num_noy, 0, 0.0d0, 0.0d0)
          end if

          if (self%ihno3_num == 0) then
            err_msg = 'do_nodoz/ihno3_num problem in the rc File.'
            call GmiPrintError  &
     &        (err_msg, .true., 1, self%ihno3_num, 0, 0, 0.0d0, 0.0d0)
          end if
        end if
        
        if (((self%io3_num  == 0) .or. (self%io3_num  /= IO3)) .and.  &
     &      (self%phot_opt == 4)) then
          err_msg = 'io3_num/IO3 problem in the rc File.'
          call GmiPrintError  &
     &      (err_msg, .true., 2, self%io3_num, IO3, 0, 0.0d0, 0.0d0)
        end if

        if ((self%num_chem == 0) .or. (self%num_chem /= NCHEM)) then
          err_msg = 'chem_opt/num_chem problem in the rc File.'
          call GmiPrintError  & 
     &      (err_msg, .true., 2, self%chem_opt, self%num_chem, 0, 0.0d0, 0.0d0)
        end if

        if ((self%num_molefrac == 0) .or. (self%num_molefrac /= NMF)) then
          err_msg = 'chem_opt/num_molefrac problem in the rc File.'
          call GmiPrintError  & 
     &      (err_msg, .true., 2, self%chem_opt, self%num_molefrac,  &
     &       0, 0.0d0, 0.0d0)
        end if
          
        if ((self%num_qjs == 0) .or. (self%num_qjs /= NUM_J)) then
          err_msg = 'chem_opt/num_qjs problem in the rc File.'
          call GmiPrintError  & 
     &      (err_msg, .true., 2, self%chem_opt, self%num_qjs, 0, 0.0d0, 0.0d0)
        end if
        
        if ((self%num_qks == 0) .or. (self%num_qks /= NUM_K)) then
          err_msg = 'chem_opt/num_qks problem in the rc File.'
          call GmiPrintError  &
     &      (err_msg, .true., 2, self%chem_opt, self%num_qks, 0, 0.0d0, 0.0d0)
        end if 
     
        if (numSpecies /= NSP) then
          err_msg = 'chem_opt/num_species problem in the rc File.'
          call GmiPrintError  &
     &      (err_msg, .true., 2, self%chem_opt, numSpecies, 0, 0.0d0, 0.0d0)
        end if
            
        if (met_opt /= 3) then
          if ((self%num_sad == 0) .or. (self%num_sad /= NSAD)) then
            err_msg = 'met_opt/num_sad problem in the rc File.'
            call GmiPrintError  &
     &        (err_msg, .true., 2, met_opt, self%num_sad, 0, 0.0d0, 0.0d0)
          end if
        end if
          
        if ((self%chem_cycle <= 0.0d0) .or. (self%chem_cycle > 1.0d0)) then
          
          err_msg = 'chem_cycle must be > 0 .and. <= 1.'
          call GmiPrintError  &
     &      (err_msg, .true., 0, 0, 0, 1, self%chem_cycle, 0.0d0)
            
        else if (self%chem_cycle < 1.0d0) then  ! subcycle chemistry

          rcycles = 1.0d0 / self%chem_cycle
          icycles = Nint (rcycles)

          if ((Abs (rcycles - icycles)) > 0.000001d0) then
            err_msg =  &
     &    'if chem_cycle < 1, 1/chem_cycle must be a whole real number.'
            call GmiPrintError  &
     &        (err_msg, .true., 0, 0, 0, 1, self%chem_cycle, 0.0d0)
          end if

        end if

      end if

      return

      end subroutine rcCheckChemistrySetting
!-----------------------------------------------------------------------------
!-------------------------------------------------------------------------
  subroutine Allocate_qqjda (self, i1, i2, ju1, j2, k1, k2)
    integer           , intent(in   )  :: i1, i2, ju1, j2, k1, k2
    type (t_Chemistry), intent(inout)  :: self
    Allocate(self%qqjda(i1:i2, ju1:j2, k1:k2, self%num_qjs))
    self%qqjda = 0.0d0
    return
  end subroutine Allocate_qqjda
!-------------------------------------------------------------------------
  subroutine Allocate_qqkda (self, i1, i2, ju1, j2, k1, k2)
    integer           , intent(in   )  :: i1, i2, ju1, j2, k1, k2
    type (t_Chemistry), intent(inout)  :: self
    Allocate(self%qqkda(i1:i2, ju1:j2, k1:k2, self%num_qks))
    self%qqkda = 0.0d0
    return
  end subroutine Allocate_qqkda
!-------------------------------------------------------------------------
  subroutine Allocate_qqjgmi (self, i1, i2, ju1, j2, k1, k2)
    integer           , intent(in   )  :: i1, i2, ju1, j2, k1, k2
    type (t_Chemistry), intent(inout)  :: self
    Allocate(self%qqjgmi(i1:i2, ju1:j2, k1:k2, self%num_qjo))
    self%qqjgmi = 0.0d0
    return
  end subroutine Allocate_qqjgmi
!-------------------------------------------------------------------------
  subroutine Allocate_qqkgmi (self, i1, i2, ju1, j2, k1, k2)
    integer           , intent(in   )  :: i1, i2, ju1, j2, k1, k2
    type (t_Chemistry), intent(inout)  :: self
    Allocate(self%qqkgmi(i1:i2, ju1:j2, k1:k2, self%num_qks))
    self%qqkgmi = 0.0d0
    return
  end subroutine Allocate_qqkgmi
!-------------------------------------------------------------------------
  subroutine Allocate_yda (self, i1, i2, ju1, j2, k1, k2)
    integer           , intent(in   )  :: i1, i2, ju1, j2, k1, k2
    type (t_Chemistry), intent(inout)  :: self
    Allocate(self%yda(i1:i2, ju1:j2, k1:k2, self%num_active))
    self%yda = 0.0d0
    return
  end subroutine Allocate_yda
!-------------------------------------------------------------------------
  subroutine Set_qqjda (self, qqjda)
    real*8          , intent(in)  :: qqjda(:,:,:,:)
    type (t_Chemistry), intent(inout)   :: self
    self%qqjda(:,:,:,:) = qqjda(:,:,:,:)
    return
  end subroutine Set_qqjda
!-----------------------------------------------------------------------
  subroutine Set_qqkda (self, qqkda)
    real*8          , intent(in)  :: qqkda(:,:,:,:)
    type (t_Chemistry), intent(inout)   :: self
    self%qqkda(:,:,:,:) = qqkda(:,:,:,:)
    return
  end subroutine Set_qqkda
!-----------------------------------------------------------------------
  subroutine Set_yda (self, yda)
    real*8          , intent(in)  :: yda(:,:,:,:)
    type (t_Chemistry), intent(inout)   :: self
    self%yda(:,:,:,:) = yda(:,:,:,:)
    return
  end subroutine Set_yda
!-----------------------------------------------------------------------
  subroutine Set_qqjgmi (self, qqjgmi)
    real*8          , intent(in)  :: qqjgmi(:,:,:,:)
    type (t_Chemistry), intent(inout)   :: self
    self%qqjgmi(:,:,:,:) = qqjgmi(:,:,:,:)
    return
  end subroutine Set_qqjgmi
!-----------------------------------------------------------------------
  subroutine Set_qqkgmi (self, qqkgmi)
    real*8          , intent(in)  :: qqkgmi(:,:,:,:)
    type (t_Chemistry), intent(inout)   :: self
    self%qqkgmi(:,:,:,:) = qqkgmi(:,:,:,:)
    return
  end subroutine Set_qqkgmi
!-----------------------------------------------------------------------
  subroutine Get_qqjda (self, qqjda)
    real*8          , intent(out)  :: qqjda(:,:,:,:)
    type (t_Chemistry), intent(in)   :: self
    qqjda(:,:,:,:) = self%qqjda(:,:,:,:)
    return
  end subroutine Get_qqjda
!-----------------------------------------------------------------------
  subroutine Get_qqkda (self, qqkda)
    real*8          , intent(out)  :: qqkda(:,:,:,:)
    type (t_Chemistry), intent(in)   :: self
    qqkda(:,:,:,:) = self%qqkda(:,:,:,:)
    return
  end subroutine Get_qqkda
!-----------------------------------------------------------------------
  subroutine Get_yda (self, yda)
    real*8          , intent(out)  :: yda(:,:,:,:)
    type (t_Chemistry), intent(in)   :: self
    yda(:,:,:,:) = self%yda(:,:,:,:)
    return
  end subroutine Get_yda
!-----------------------------------------------------------------------
  subroutine Get_qqjgmi (self, qqjgmi)
    real*8          , intent(out)  :: qqjgmi(:,:,:,:)
    type (t_Chemistry), intent(in)   :: self
    qqjgmi(:,:,:,:) = self%qqjgmi(:,:,:,:)
    return
  end subroutine Get_qqjgmi
!-----------------------------------------------------------------------
  subroutine Get_qqkgmi (self, qqkgmi)
    real*8          , intent(out)  :: qqkgmi(:,:,:,:)
    type (t_Chemistry), intent(in)   :: self
    qqkgmi(:,:,:,:) = self%qqkgmi(:,:,:,:)
    return
  end subroutine Get_qqkgmi
!-----------------------------------------------------------------------
  subroutine Get_oz_eq_synoz_opt (self, oz_eq_synoz_opt)
    integer         , intent(out)  :: oz_eq_synoz_opt
    type (t_Chemistry), intent(in)   :: self
    oz_eq_synoz_opt = self%oz_eq_synoz_opt
    return
  end subroutine Get_oz_eq_synoz_opt
!-------------------------------------------------------------------------
  subroutine Get_sad_opt (self, sad_opt)
    integer         , intent(out)  :: sad_opt
    type (t_Chemistry), intent(in)   :: self
    sad_opt = self%sad_opt
    return
  end subroutine Get_sad_opt
!-------------------------------------------------------------------------
  subroutine Get_phot_opt (self, phot_opt)
    integer         , intent(out)  :: phot_opt
    type (t_Chemistry), intent(in)   :: self
    phot_opt = self%phot_opt
    return
  end subroutine Get_phot_opt
!-------------------------------------------------------------------------
  subroutine Get_do_chem_grp (self, do_chem_grp)
    logical          , intent(out)  :: do_chem_grp
    type (t_Chemistry), intent(in )  :: self
    do_chem_grp = self%do_chem_grp
    return
  end subroutine Get_do_chem_grp
!-------------------------------------------------------------------------
  subroutine Get_synoz_threshold (self, synoz_threshold)
    real*8           , intent(out)  :: synoz_threshold
    type (t_Chemistry), intent(in )  :: self
    synoz_threshold = self%synoz_threshold
    return
  end subroutine Get_synoz_threshold
!-------------------------------------------------------------------------
  subroutine Get_t_cloud_ice (self, t_cloud_ice)
    real*8           , intent(out)  :: t_cloud_ice
    type (t_Chemistry), intent(in )  :: self
    t_cloud_ice = self%t_cloud_ice
    return
  end subroutine Get_t_cloud_ice
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
  end module GmiChemistryMethod_mod
