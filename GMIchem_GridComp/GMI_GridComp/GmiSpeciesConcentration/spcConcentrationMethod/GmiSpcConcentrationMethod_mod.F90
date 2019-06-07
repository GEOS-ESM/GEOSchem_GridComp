!-------------------------------------------------------------------------
!BOP
!
! !MODULE: GmiSpcConcentrationMethod_mod
!
! !INTERFACE:
!
#include "MAPL_Generic.h"
!
  module GmiSpcConcentrationMethod_mod
!
! !USES:
      use ESMF
      use MAPL_Mod
      use GmiESMFrcFileReading_mod, only : rcEsmfReadTable, rcEsmfReadLogical
      use GmiArrayBundlePointer_mod, only : t_GmiArrayBundle, CleanArrayPointer
      use GmiArrayBundlePointer_mod, only : setArrayPointer
      use GmiGrid_mod              , only : t_gmiGrid
      use GmiGrid_mod              , only : Get_i1, Get_i2, Get_ju1, Get_j2, Get_k1, Get_k2
      use GmiGrid_mod              , only : Get_ilong, Get_ilat, Get_ivert, Get_i1_gl, Get_ju1_gl
      use GmiPrintError_Mod,         only : GmiPrintError
      use GmiSpeciesRegistry_mod   , only : set_labelsSpecies, set_numSpecies
!
  implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
!
!  private
  public  :: InitializeSpcConcentration, FinalizeSpcConcentration
  public  :: ResetFixedConcentration   , constructConcentrationBundle
  public  :: isfixedConcentration      , Set_const_label
  public  :: Get_concentration      , Set_concentration
  public  :: Get_const_opt          !, Get_mw
  public  :: Get_const_init_val     , Get_const_infile_name
  public  :: Get_const_var_name     , Get_const_labels
  public  :: Get_num_fixed_const
  public  :: Get_fixed_const_map    , Get_fixed_const_infile_name
  public  :: Get_fixed_const        , Get_io3_num           
  public  :: Get_tracer_opt         , Get_efol_time           
  public  :: Get_tr_source_land     , Get_tr_source_ocean        
!
! !PUBLIC MEMBER DATA:
!
  public  :: t_SpeciesConcentration
!
! !PRIVATE DATA MEMBERS:
  integer, pointer :: fixedConstMap(:)
  integer          :: numFixedConst
  integer          :: numMoleFrac
  integer          :: numChem

#     include "GmiParameters.h"
#     include "setkin_par.h"
#     include "setkin_lchem.h"

  type t_SpeciesConcentration
!    private
    ! For tracer experiments
                                           ! tracer run option
    integer                                :: tracer_opt
                                           ! e-folding time of the tracer (in days)
    real*8                                 :: efol_time
                                           ! land  source of the tracer
    real*8                                 :: tr_source_land 
                                           ! ocean source of the tracer
    real*8                                 :: tr_source_ocean
                                           ! spec. conc. option
    integer                                :: const_opt        
                                           ! number of species to input
    integer                                :: num_const_inrecs
                                           ! number of fixed constituents
    integer                                :: num_fixed_const
                                           !  mapping of fixed spc. conc. number to spc. conc. number
    integer, pointer                       :: fixed_const_map(:)
                                           ! array of values to initialize each spc. conc. to mixing ratio.
    real*8                , pointer        :: const_init_val(:)
                                           ! Spc. conc. string labels
    character(len=MAX_LENGTH_VAR_NAME)     :: const_var_name    
                                           ! spc. conc. input file name
    character(len=MAX_LENGTH_FILE_NAME)    :: const_infile_name 
                                           ! fixed spc. conc. input file name
    character(len=MAX_LENGTH_FILE_NAME)    :: fixed_const_infile_name
                                           ! constituent string labels
    character(len=MAX_LENGTH_SPECIES_NAME), pointer :: const_labels(:) => null()
                                           ! index of ozone
    integer                                :: io3_num  
                          ! value of production of a species
    real*8                , pointer        :: prod(:,:,:) => null()
                          ! value of loss       of a species 
    real*8                , pointer        :: loss(:,:,:) => null()
    real*8                , pointer        :: fixed_const(:,:,:,:) => null()
                                           ! Species concentration of the surface
    real*8                , pointer        :: concentrationSurf    (:,:,:) => null()
                                           ! Column species concentration for troposphere
    real*8                , pointer        :: concentrationColTrop (:,:,:) => null()
                                           ! Column species concentration for strat/trop
    real*8                , pointer        :: concentrationColCombo(:,:,:) => null()
    type(t_GmiArrayBundle), pointer        :: concentration(:) => null()
    type(t_GmiArrayBundle), pointer        :: net_cum_tend(:,:) => null()
  end type t_SpeciesConcentration

  logical, save :: pr_diag

! !DESCRIPTION:
!
! !AUTHOR:
!
! !REVISION HISTORY:
!
!EOP
!-------------------------------------------------------------------------
  CONTAINS
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  ReadSpcConcentrationResourceFile
!
! !INTERFACE:
!
      subroutine ReadSpcConcentrationResourceFile (self, config, loc_proc,     &
     &               rootProc, numSpecies)
!
! !USES:
      use GmiCheckNamelistFile_mod, only : CheckNamelistOptionRange
      use GmiSpeciesRegistry_mod, only : getSpeciesIndex, UNKNOWN_SPECIES
      use GmiStringManipulation_mod, only : constructListNames
!
      implicit none
!
! !INPUT PARAMETERS:
      logical,          intent(in) :: rootProc
      integer,          intent(in) :: loc_proc, numSpecies
!
! !INPUT/OUTPUT PARAMETERS:
      type (ESMF_Config)           , intent(inOut) :: config
      type (t_SpeciesConcentration), intent(inOut) :: self
!
! !DESCRIPTION:
! Reads in Species concentration related variables from the resource file.
!
! !LOCAL VARIABLES:
      integer :: STATUS, RC, ic, neg_marker
      character (len=MAX_LENGTH_SPECIES_NAME) :: tempWord
      character (len=MAX_LENGTH_SPECIES_NAME) :: tempListNames(NSP)
      character (len=MAX_STRING_LENGTH      ) :: fixedConcentrationSpeciesNames
      character(len=ESMF_MAXSTR) :: IAm, err_msg
!EOP
!-------------------------------------------------------------------------
!BOC
      IAm = "ReadSpcConcentrationResourceFile"

      call rcEsmfReadLogical(config, pr_diag, &
     &               "pr_diag:", default = .false., rc=STATUS )
      VERIFY_(STATUS)

      if (pr_diag) Write(6,*) IAm, 'called by ', loc_proc

      !################################
      ! Begin reading the resource file
      !################################

      IF(rootProc) THEN
       WRITE (6,*) 'Reading section GmiSpeciesConcentration from rc file'
      END IF

      ! --------------------------------------------------------------
      ! const_opt
      !   0:  No fixed concentration species
      !   1:  set const values to const_init_val
      !   2:  read in const values [default]
      !   3:  solid body rotation
      !   4:  dummy test pattern with linear slope in each dimension
      !   5:  exponential in vertical (decays with height)
      !   6:  sin in latitude (largest at equator)
      !   7:  linear vertical gradient
      !   8:  sin in latitude (largest at equator) + vertical gradient
      ! --------------------------------------------------------------
      
      call ESMF_ConfigGetAttribute(config, self%const_opt, &
     &               label   = "const_opt:",&
     &               default = 2, rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, self%const_infile_name, &
     &               label   = "const_infile_name:",&
     &               default = ' ', rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, self%const_var_name, &
     &               label   = "const_var_name:",& 
     &               default = 'const', rc=STATUS )
      VERIFY_(STATUS)

      allocate(self%const_init_val(numSpecies)) 
      self%const_init_val(:)  = 1.0d-30

      call rcEsmfReadTable(config, self%const_init_val, &
     &                     "const_init_val::", rc=STATUS)
      
      call rcEsmfReadTable(config, fixedConcentrationSpeciesNames, &
     &               "fixedConcentrationSpeciesNames::", rc=STATUS)
      
      call ESMF_ConfigGetAttribute(config, self%fixed_const_infile_name, &
     &                label   = "fixed_const_infile_name:", &
     &                default = '', rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, self%io3_num, &
     &                label   = "io3_num:", &
     &                default = 0, rc=STATUS )
      VERIFY_(STATUS)

      allocate(self%const_labels(numSpecies))

      self%const_labels(1:numSpecies) = lchemvar(1:numSpecies)

      call rcEsmfReadTable(config, self%const_labels, &
     &                     "const_labels::", rc=STATUS)

      ! Copy const_labels to labelsSpecies
      ! ----------------------------------
        CALL set_numSpecies(numSpecies)
        CALL set_labelsSpecies(self%const_labels)

!     -----------------------------------------------------
!     tracer_opt
!       0:  no tracer run
!       1:  tracer run with decay of efol_time
!       2:  Age of Air tracer run
!       3:  HTAP Tagged CO tracer run
!       4:  YY Tagged CH4 tracer run
!       5:  BND Tagged CO tracer run
!       6:  BND Tagged CO tracer with additional 5, 10, & 15 day decays also run
!     -----------------------------------------------------

      call ESMF_ConfigGetAttribute(config, self%tracer_opt, &
     &                label   = "tracer_opt:", &
     &                default = 0, rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, self%efol_time, &
     &                label   = "efol_time:", &
     &                default = 0.0d0, rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, self%tr_source_land, &
     &                label   = "tr_source_land:", &
     &                default = 0.0d0, rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, self%tr_source_ocean, &
     &                label   = "tr_source_ocean:", &
     &                default = 0.0d0, rc=STATUS )
      VERIFY_(STATUS)

      ! ---------------------------------------------------------------
      ! Check option ranges.  Note that as new options are added, these
      ! range checks will have to be modified.
      ! ---------------------------------------------------------------

      call CheckNamelistOptionRange ('const_opt', self%const_opt, 0, 8)
      call CheckNamelistOptionRange ('tracer_opt', self%tracer_opt, 0, 6)

      if (self%const_opt == 1) then
!       ---------------------------------------------------------
!       Check for any negative const_init_val() marker set in the
!       resource input file.  If one is found, set all the
!       const_init_val's from the negative value on to the value
!       preceding the negative value.
!       ---------------------------------------------------------
        neg_marker = INT_DUM_VALUE
        do ic = 2, numSpecies
          if ((neg_marker == INT_DUM_VALUE) .and. (self%const_init_val(ic) < 0)) then
            neg_marker = ic
          end if
          if (neg_marker /= INT_DUM_VALUE) then
            self%const_init_val(ic) = self%const_init_val(neg_marker-1)
          end if
        end do 
      end if

      self%num_fixed_const = 0

      DoingSomeFixed: IF(self%const_opt > 0) THEN
      
      ! Set the initial value of the list
       tempListNames(:) = ''

      ! Construct the list of names using the long string
       call constructListNames(tempListNames, fixedConcentrationSpeciesNames)
     
       self%num_fixed_const = count(tempListNames /= '')

       IF (rootProc .AND. self%num_fixed_const > 0) THEN
         PRINT *," "
         PRINT *,"Found ",self%num_fixed_const," fixed-concentration specie(s)"
         PRINT *," Specie names:"
         DO ic = 1,self%num_fixed_const
            PRINT *,"  ",ic,"  ",TRIM(tempListNames(ic))
         END DO
       END IF

       if (self%num_fixed_const > 0) then
         allocate(self%fixed_const_map(self%num_fixed_const))
         self%fixed_const_map(:) = 0
     
         do ic = 1, self%num_fixed_const
            self%fixed_const_map(ic) = getSpeciesIndex(tempListNames(ic))
         end do
       end if
     
       IF (rootProc) THEN
         PRINT *,"Fixed_const_map:"
         PRINT *,self%fixed_const_map
         PRINT *," "
       END IF

      ELSE

       IF (rootProc) THEN
         PRINT *," "
         PRINT *," No fixed-concentration species"
         PRINT *," "
       END IF

      END IF DoingSomeFixed

      return

      end subroutine ReadSpcConcentrationResourceFile
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: InitializeSpcConcentration
!
! !INTERFACE:
!
      subroutine InitializeSpcConcentration (self, gmiGrid, config, numSpecies,&
     &                     num_molefrac, num_chem, loc_proc )
!
  implicit none
!
! !INPUT PARAMETERS:
       integer          , intent(in) :: numSpecies, loc_proc
       integer          , intent(in) :: num_molefrac, num_chem
       type (t_gmiGrid ), intent(in) :: gmiGrid 
!
! !INPUT/OUTPUT PARAMETERS:
      type (ESMF_Config)           , intent(inOut) :: config
      type (t_SpeciesConcentration), intent(inOut) :: self
!
! !LOCAL VARIABLES:
      integer :: ix, ic, id
      integer :: i1, i2, ju1, j2, k1, k2, ilong, ilat, ivert, i1_gl, ju1_gl
      logical :: rootProc
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
      call Get_ilong (gmiGrid, ilong)
      call Get_ilat  (gmiGrid, ilat )
      call Get_ivert (gmiGrid, ivert)
      call Get_i1_gl (gmiGrid, i1_gl)
      call Get_ju1_gl(gmiGrid, ju1_gl)

      rootProc = MAPL_AM_I_ROOT()

      call ReadSpcConcentrationResourceFile (self, config, loc_proc, rootProc, &
     &            numSpecies)

      call Allocate_concentration(self, i1, i2, ju1, j2, k1, k2, numSpecies)

      IF(self%num_fixed_const > 0) THEN
         IF(rootProc) THEN
          PRINT *,"Allocating fixed_const for "
          PRINT *," Number of species: ",self%num_fixed_const
          PRINT *," "
         END IF
         CALL Allocate_fixed_const(self, i1, i2, ju1, j2, k1, k2)
      END IF

      ! Initialization needed for the function isFixedConcentration

      allocate(fixedConstMap(1:self%num_fixed_const))
      fixedConstMap(:) = self%fixed_const_map(1:self%num_fixed_const)
      numFixedConst    = self%num_fixed_const
      numMoleFrac      = num_molefrac
      numChem          = num_chem

      return

      end subroutine InitializeSpcConcentration
!EOC
!-------------------------------------------------------------------------
  subroutine FinalizeSpcConcentration (self)

  implicit none

  type (t_SpeciesConcentration)   , intent(inOut) :: self

  PRINT*,'  Finalize SpcConcentration'

  return

  end subroutine FinalizeSpcConcentration
!-------------------------------------------------------------------------
!BOP
  subroutine resetFixedConcentration (self, gmiClock, gmiGrid, numSpecies)
!
! !USES:
!
  use GmiTimeControl_mod, only : GmiSplitDateTime, Get_curGmiDate, t_gmiClock
      implicit none
! !INPUT PARAMETERS:
  integer            , intent(in) :: numSpecies
  type (t_gmiGrid )  , intent(in) :: gmiGrid 
  type (t_gmiClock)  , intent(in) :: gmiClock 
!
! !INPUT/OUTPUT PARAMETERS:
  type (t_SpeciesConcentration), intent(inOut) :: self
!
! !DESCRIPTION:
! This routines fixes values of certain species concentrations.
!
! !LOCAL VARIABLES:
  integer :: i1, i2, ju1, j2, k1, k2, nymd
  integer :: ic, im, icx, idumday, idumyear
!
! !AUTHOR:
! !HISTORY:
!BOP
!-------------------------------------------------------------------------
!BOC
  call Get_i1 (gmiGrid, i1 )
  call Get_i2 (gmiGrid, i2 )
  call Get_ju1(gmiGrid, ju1)
  call Get_j2 (gmiGrid, j2 )
  call Get_k1 (gmiGrid, k1 )
  call Get_k2 (gmiGrid, k2 )

  call Get_curGmiDate(gmiClock, nymd)

  if (self%num_fixed_const /= 0) then
     do ic = 1, self%num_fixed_const
        icx = self%fixed_const_map(ic)
        self%concentration(icx)%pArray3D(:,:,:) = self%fixed_const(:,:,:,ic)
     end do
  end if

  return

  end subroutine resetFixedConcentration

!EOC
!-------------------------------------------------------------------------
!BOP
  subroutine constructConcentrationBundle &
          (self, concArray, i1, i2, ju1, j2, k1, k2, numSpecies)
  implicit none
! !INPUT PARAMETERS:
  integer                      , intent(in) :: i1, i2, ju1, j2, k1, k2
  integer                      , intent(in) :: numSpecies
  real*8                       , intent(in) :: concArray(i1:i2,ju1:j2,k1:k2,numSpecies)
!
! !INPUT/OUTPUT PARAMETERS:
  type (t_SpeciesConcentration), intent(inOut) :: self
!
! !DESCRIPTION:
!
! !LOCAL VARIABLES:
  integer :: i
!
!EOP
!-------------------------------------------------------------------------
!BOC
  do i = 1, numSpecies
     call setArrayPointer(self%concentration(i), concArray(i1,ju1,k1,i), &
                 i1 ,i2, ju1, j2, k1, k2)
  end do  
  return
  end subroutine constructConcentrationBundle
!EOC
!-------------------------------------------------------------------------
!BOP
  function isFixedConcentration (ic)
  implicit none
! !INPUT PARAMETERS:
  integer :: ic  ! specie concentration array index
  logical :: isFixedConcentration
! !DESCRIPTION:
! This function checks to see if a particular specie concentration
! has fixed values.
!
! !LOCAL VARIABLES:
  integer :: ix
!
!EOP
!-------------------------------------------------------------------------
!BOC
  isFixedConcentration = .false.
  ixloop: do ix = 1, numFixedConst
     if (fixedConstMap(ix) == ic) then
        isFixedConcentration = .true.
        exit ixloop
     end if
  end do ixloop
  if ((ic > numMoleFrac) .and. (ic <= numChem)) then
!------------------------------------------------------------------
!Take care of the special case of fixed species such as O2, N2, ad.
!------------------------------------------------------------------
     isFixedConcentration = .true.
   end if
   return
  end function isFixedConcentration
!EOC
!-------------------------------------------------------------------------
  subroutine Allocate_concentration (self, i1, i2, ju1, j2, k1, k2, numSpecies)
    implicit none
    integer          , intent(in   ) :: i1, i2, ju1, j2, k1, k2, numSpecies
    type (t_SpeciesConcentration), intent(inOut) :: self
    integer :: ic
    Allocate(self%concentration(numSpecies))
    do ic = 1, numSpecies
       Allocate(self%concentration(ic)%pArray3D(i1:i2, ju1:j2, k1:k2))
!       self%concentration(ic)%pArray3D(i1:i2, ju1:j2, k1:k2) = 0.0d0
    end do
    return
  end subroutine Allocate_concentration
!-------------------------------------------------------------------------
  subroutine Get_concentration (self, concentration)
    implicit none
    type (t_GmiArrayBundle), pointer :: concentration (:)
    type (t_SpeciesConcentration), intent(in)   :: self
    concentration => self%concentration
    return
  end subroutine Get_concentration
!-------------------------------------------------------------------------
  subroutine Set_concentration (self, concentration)
    implicit none
    type (t_GmiArrayBundle), pointer :: concentration (:)
    type (t_SpeciesConcentration), intent(inOut) :: self
    self%concentration => concentration
    return
  end subroutine Set_concentration
!-------------------------------------------------------------------------
  subroutine Allocate_fixed_const (self, i1, i2, ju1, j2, k1, k2)
    implicit none
    integer        , intent(in )  :: i1, i2, ju1, j2, k1, k2
    type (t_SpeciesConcentration), intent(inOut)   :: self
    Allocate(self%fixed_const(i1:i2, ju1:j2, k1:k2, self%num_fixed_const))
    self%fixed_const = 0.0d0
    return
  end subroutine Allocate_fixed_const
!-------------------------------------------------------------------------
  subroutine Get_fixed_const (self, fixed_const)
    implicit none
    real*8         , intent(out)  :: fixed_const (:,:,:,:)
    type (t_SpeciesConcentration), intent(in)   :: self
    fixed_const(:,:,:,:) = self%fixed_const(:,:,:,:)
    return
  end subroutine Get_fixed_const
!-------------------------------------------------------------------------
  subroutine Allocate_fixed_const_map (self, numSpecies)
    implicit none
    integer        , intent(in )  :: numSpecies
    type (t_SpeciesConcentration), intent(inOut)   :: self
    Allocate(self%fixed_const_map(1:numSpecies))
    self%fixed_const_map = 0
    return
  end subroutine Allocate_fixed_const_map
!-------------------------------------------------------------------------
  subroutine Get_const_init_val (self, const_init_val, numSpecies)
    implicit none
    integer        , intent(in )  :: numSpecies
    real*8         , intent(out)  :: const_init_val (:)
    type (t_SpeciesConcentration), intent(in)   :: self
    const_init_val(1:numSpecies) = self%const_init_val(1:numSpecies)
    return
  end subroutine Get_const_init_val
!-------------------------------------------------------------------------
  subroutine Get_num_fixed_const (self, num_fixed_const)
    implicit none
    integer        , intent(out)  :: num_fixed_const 
    type (t_SpeciesConcentration), intent(in)   :: self
    num_fixed_const = self%num_fixed_const
    return
  end subroutine Get_num_fixed_const
!-------------------------------------------------------------------------
  subroutine Get_fixed_const_map (self, fixed_const_map, numSpecies)
    implicit none
    integer        , intent(in )  :: numSpecies
    integer        , intent(out)  :: fixed_const_map (:)
    type (t_SpeciesConcentration), intent(in)   :: self
    fixed_const_map(1:numSpecies) = self%fixed_const_map(1:numSpecies)
    return
  end subroutine Get_fixed_const_map
!-------------------------------------------------------------------------
  subroutine Get_const_opt (self, const_opt)
    implicit none
    integer        , intent(out)  :: const_opt 
    type (t_SpeciesConcentration), intent(in)   :: self
    const_opt = self%const_opt
    return
  end subroutine Get_const_opt
!-------------------------------------------------------------------------
  subroutine Get_io3_num (self, io3_num)
    implicit none
    integer        , intent(out)  :: io3_num
    type (t_SpeciesConcentration), intent(in)   :: self
    io3_num = self%io3_num
    return
  end subroutine Get_io3_num
!-------------------------------------------------------------------------
  subroutine Get_fixed_const_infile_name (self, fixed_const_infile_name)
    implicit none
    character (len=*), intent(out)  :: fixed_const_infile_name
    type (t_SpeciesConcentration)    , intent(in)   :: self
    fixed_const_infile_name = self%fixed_const_infile_name
    return
  end subroutine Get_fixed_const_infile_name
!-------------------------------------------------------------------------
  subroutine Get_const_infile_name (self, const_infile_name)
    implicit none
    character (len=*), intent(out)  :: const_infile_name
    type (t_SpeciesConcentration)    , intent(in)   :: self
    const_infile_name = self%const_infile_name
    return
  end subroutine Get_const_infile_name
!-------------------------------------------------------------------------
  subroutine Get_const_labels (self, const_labels, numSpecies)
    implicit none
    integer          , intent(in )  :: numSpecies
    character (len=*), intent(out)  :: const_labels(numSpecies)
    type (t_SpeciesConcentration)    , intent(in)   :: self
    const_labels(1:numSpecies) = self%const_labels(1:numSpecies)
    return
  end subroutine Get_const_labels
!-------------------------------------------------------------------------
  subroutine Set_const_label (self, ic, label)
    implicit none
    integer          , intent(in)  :: ic
    character (len=*), intent(in)  :: label
    type (t_SpeciesConcentration)    , intent(inout)   :: self
    self%const_labels(ic) = TRIM(label)
    return
  end subroutine Set_const_label
!-------------------------------------------------------------------------
  subroutine Get_const_var_name (self, const_var_name)
    implicit none
    character (len=*), intent(out)  :: const_var_name
    type (t_SpeciesConcentration)    , intent(in)   :: self
    const_var_name = self%const_var_name
    return
  end subroutine Get_const_var_name
!-------------------------------------------------------------------------
  subroutine Get_tracer_opt (self, tracer_opt)
    implicit none
    integer        , intent(out)  :: tracer_opt
    type (t_SpeciesConcentration), intent(in)   :: self
    tracer_opt = self%tracer_opt
    return
  end subroutine Get_tracer_opt
!-------------------------------------------------------------------------
  subroutine Get_efol_time (self, efol_time)
    implicit none
    real*8         , intent(out)  :: efol_time
    type (t_SpeciesConcentration), intent(in)   :: self
    efol_time = self%efol_time
    return
  end subroutine Get_efol_time
!-------------------------------------------------------------------------
  subroutine Get_tr_source_land (self, tr_source_land)
    implicit none
    real*8         , intent(out)  :: tr_source_land
    type (t_SpeciesConcentration), intent(in)   :: self
    tr_source_land = self%tr_source_land
    return
  end subroutine Get_tr_source_land
!-------------------------------------------------------------------------
  subroutine Get_tr_source_ocean (self, tr_source_ocean)
    implicit none
    real*8         , intent(out)  :: tr_source_ocean
    type (t_SpeciesConcentration), intent(in)   :: self
    tr_source_ocean = self%tr_source_ocean
    return
  end subroutine Get_tr_source_ocean
!-------------------------------------------------------------------------
  end module GmiSpcConcentrationMethod_mod
