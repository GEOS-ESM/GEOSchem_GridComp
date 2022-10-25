#include "MAPL_Generic.h"
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: GmiEmissionMethod_mod
!
! !INTERFACE:
!
  module GmiEmissionMethod_mod
!
! !USES:

      use ESMF, only : ESMF_Config, ESMF_MAXSTR, ESMF_ConfigGetAttribute
      use MAPL
      use GmiTimeControl_mod  , only : t_GmiClock, Get_curGmiDate
      use GmiGrid_mod, only : t_gmiGrid
      use GmiGrid_mod, only : Get_i1, Get_i2, Get_ju1, Get_j2, Get_k1, Get_k2
      use GmiGrid_mod, only : Get_i1_gl, Get_i2_gl, Get_ju1_gl, Get_j2_gl
      use GmiGrid_mod, only : Get_ilo, Get_ihi, Get_julo, Get_jhi
      use GmiGrid_mod, only : Get_ilong, Get_ilat, Get_ivert
      use GmiPrintError_mod, only : GmiPrintError
      use GmiCheckRange_mod, only : CheckRange3d
      use GmiESMFrcFileReading_mod, only : rcEsmfReadTable
      use GmiArrayBundlePointer_mod, only : t_GmiArrayBundle, CleanArrayPointer
      use GmiSpcConcentrationMethod_mod, only : t_SpeciesConcentration
      use GmiSpcConcentrationMethod_mod, only : Get_concentration, Set_concentration
      use GmiSpcConcentrationMethod_mod, only : Get_tracer_opt, Get_tr_source_land
      use GmiSpcConcentrationMethod_mod, only : Get_tr_source_ocean
!
      use, intrinsic :: iso_fortran_env, only: REAL64
!
      implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
!
      private
      public  :: InitializeEmission, initReadEmission
      public  :: RunEmission
      public  :: FinalizeEmission
    
!     public  :: Get_emissionArray
      public  :: Get_lightning_opt       , Get_emiss_aero_opt
      public  :: Get_emiss_dust_opt      , Get_num_emiss
      public  :: Get_ndust               , Get_naero
      public  :: Get_do_gcr
      public  :: Get_do_semiss_inchem
!
! !PUBLIC MEMBER DATA:
      public  :: t_Emission

# include "GmiParameters.h"
# include "gmi_emiss_constants.h"
# include "gmi_phys_constants.h"
# include "gmi_time_constants.h"

  type t_Emission
    integer             :: idaySoilType
    logical             :: firstBiogenicBase
                           ! current record number for the
    integer             :: curEmissionFileRecord  
                           ! main emission file.
    integer             :: semiss_inchem_flag
                           ! for fossil fuel scaling
    logical             :: doScaleNOffEmiss    
                           ! for biomass burning scaling
    logical             :: doScaleNObbEmiss 
    real*8 , pointer    :: scFacNOff(:,:,:) => null()
    real*8 , pointer    :: scFacNObb(:,:,:) => null()
                           ! scale factor NO fossil fuel
    character (len=MAX_LENGTH_FILE_NAME) :: scFactorNOff_infile_name
                           ! scale factor NO biomass burning
    character (len=MAX_LENGTH_FILE_NAME) :: scFactorNObb_infile_name
                           ! Flag to indicate climatological emissions are per unit area (default T)
    LOGICAL             :: clim_emiss_by_area
                           ! Names of the emitted species
    character (len=10 ) :: emittedSpeciesNames (MAX_NUM_CONST)
                           ! emission conversion flag
    integer             :: emiss_conv_flag 
                           ! emission conversion factor (s^-1)
    real*8              :: emiss_conv_fac 
                           ! value to initialize emiss array to (kg/s)
    real*8              :: emiss_init_val
                           ! sulfur aerosol emiss option
    integer             :: emiss_aero_opt   
                           ! number of aerosols
    integer             :: naero           
                           ! mapping of aerosol emiss number to const species #
    integer             :: emiss_map_aero(MAX_NUM_CONST) 
                           ! aerosol (sulf. code) input file name
    character (len=MAX_LENGTH_FILE_NAME) :: emiss_aero_infile_name 
                           ! dust (sulf. code) input file name
    character (len=MAX_LENGTH_FILE_NAME) :: emiss_dust_infile_name 
                           ! do Galactic Cosmic Ray
    logical             :: do_gcr 
                           ! ! Galactic Cosmic Ray  input file name
    character (len=MAX_LENGTH_FILE_NAME) :: gcr_infile_name 
    character (len=MAX_LENGTH_FILE_NAME) :: soil_infile_name              ! soil type            input file name
    character (len=MAX_LENGTH_FILE_NAME) :: isopconv_infile_name          ! isoprene convert     input file name
    character (len=MAX_LENGTH_FILE_NAME) :: monotconv_infile_name         ! monoterpene convert  input file name
    real*8              :: isop_scale (12)               ! array of monthly isoprene scaling coefficients

    CHARACTER(LEN=64):: emissionSpeciesNames(MAX_NUM_CONST)  ! Names of emssions on netCDF/hdf file #

    integer          :: lightning_opt		         ! lightning option
    real*8 , pointer :: lightning_NO   (:,:,:) => null() ! NO produced by parameterized lightning

    real*8 , pointer    :: emissAero   (:,:,:)   => null() ! used in sulfur chemistry
    real*8 , pointer    :: emissDust_t (:,:,:,:) => null() ! used in sulfur chemistry
    real*8 , pointer    :: emissAero_t (:,:,:,:) => null() ! used in sulfur chemistry

    integer          :: emiss_timpyr             ! emission times per year
                                                 ! (1 => yearly, 12 => monthly)
    integer          :: emiss_opt                ! emission option
    integer          :: emiss_in_opt             ! emission input option
    integer          :: emiss_map(MAX_NUM_CONST) ! mapping of emission number to const species #
    integer          :: emiss_map_dust(MAX_NUM_CONST) ! mapping of dust emiss number to const species #
    integer          :: emissionSpeciesLayers(MAX_NUMBER_SPECIES)
    integer          :: ndust                     ! number of dust
    integer          :: nst_dust                  ! starting index for dust
    integer          :: nt_dust                   ! ending   index for dust
    logical          :: do_ShipEmission           ! do ship emissions?
    logical          :: doMEGANemission           ! do MEGAN emissions?
    logical          :: doMEGANviaHEMCO           ! use HEMCO for MEGAN? 
!.sds
    integer          :: GmiDMSEmissIndex          ! GMI DMS index in emissionArray
    CHARACTER(LEN=MAX_LENGTH_FILE_NAME):: emissionPointFilenames(MAX_NUM_CONST) ! Names of files of point emissions
    integer          :: num_point_emiss           ! number of point emissions species
    integer          :: num_point_start           ! number in emissions list that point sources start
    CHARACTER(LEN=32):: num_point_type(MAX_NUM_CONST)! number in emissions list that point sources start
!.sds.end
    character (len=MAX_LENGTH_FILE_NAME) ::      laiMEGAN_InfileName ! Inpuf file name for AVHRR leaf-area-indices
    character (len=MAX_LENGTH_FILE_NAME) ::   aefMboMEGAN_InfileName ! Annual emission factor for methyl butenol input file name
    character (len=MAX_LENGTH_FILE_NAME) ::  aefIsopMEGAN_InfileName ! Annual emission factor for isoprene input file name
    character (len=MAX_LENGTH_FILE_NAME) :: aefMonotMEGAN_InfileName ! Annual emission factor for monoterpenes input file name
    character (len=MAX_LENGTH_FILE_NAME) ::  aefOvocMEGAN_InfileName ! Annual emission factor for other biogenic VOCs input file name
    integer          :: days_btw_m     ! days between midmonths in the LAI data
    real*8 , pointer :: isoLai    (:,:) => null()  ! AVHRR LAI data for the current day
    real*8 , pointer :: isoLaiPrev(:,:) => null()  ! AVHRR LAI data for the previous month
    real*8 , pointer :: isoLaiCurr(:,:) => null()  ! AVHRR LAI data for the current  month
    real*8 , pointer :: isoLaiNext(:,:) => null()  ! AVHRR LAI data for the next     month
    real   , pointer :: isoLaiYear(:,:,:) => null()! AVHRR LAI data storage (accommodates use of ExtData.rc)
    real*8 , pointer :: aefMbo   (:,:) => null()   ! Annual emission factor for methyl butenol
    real*8 , pointer :: aefIsop  (:,:) => null()   ! Annual emission factor for isoprene
    real*8 , pointer :: aefOvoc  (:,:) => null()   ! Annual emission factor for other biogenic VOCs
    real*8 , pointer :: aefMonot (:,:) => null()   ! Annual emission factor for monoterpenes
    logical          :: do_semiss_inchem           ! do surface emissions inside chemistry solver?
    integer          :: emiss_dust_opt                ! sulfur dust emiss option
    integer          :: num_emiss                     ! number of emissions
    integer          :: ship_o3_index                 ! index of ship O3   in emissionArray
    integer          :: ship_hno3_index               ! index of ship HNO3 in emissionArray
    real*8 , pointer :: emiss          (:,:,:,:)   => null() ! array of emissions (kg/s)
    real*8 , pointer :: emissDust     (:,:,:)     => null() ! used in sulfur chemistry
    real*8 , pointer :: emiss_isop     (:,:)       => null() ! isoprene    emissions (kg/s)
    real*8 , pointer :: emiss_monot    (:,:)       => null() ! monoterpene emissions (kg/s)
    real*8 , pointer :: emiss_nox      (:,:)       => null() ! NOx         emissions (kg/s)
    real*8 , pointer :: emiss_o3       (:,:)       => null() ! ozone       emissions
    real*8 , pointer :: emiss_hno3     (:,:)       => null() ! hno3        emissions
    integer          :: ncon_soil (NVEGTYPE)           ! Olson -> soil type
                                                       !    1 => water/desert/ice
                                                       !    2 => tropical rain forest
                                                       !    3 => conifers
                                                       !    4 => dry deciduous
                                                       !    5 => other deciduous
                                                       !    6 => woodland
                                                       !    7 => grassland
                                                       !    8 => agriculture (other than rice)
                                                       !    9 => rice paddies
                                                       !   10 => wetland/tundra
    real*8 :: exp_fac(NPULSE)    ! pulsing decay per time step (day^-1)
    integer, pointer :: index_soil  (:,:,:) => null()
    integer, pointer :: ireg        (:,:)   => null()  ! number of land types in a grid square
    integer, pointer :: iuse        (:,:,:) => null()  ! fraction of grid box area occupied by land type
    integer, pointer :: iland       (:,:,:) => null()  ! land type id in grid square for ireg land types
    real*8           :: coeff_isop   (NPOLY)           ! coefficients used for polynomial fit
    real*8           :: convert_monot(NVEGTYPE)        ! monoterpene emissions by landtype (atomsC/cm^2leaf/s)
    real*8           :: convert_isop(NVEGTYPE)        ! isoprene emissions by landtype (atomsC/cm^2leaf/s)
    real*8 , pointer :: soil_fert   (:,:)   => null()  ! fertilizers (ng N/m^2/s)
    real*8 , pointer :: soil_precip (:,:)   => null()  ! two months of observed precip (mm/day/box)
    real*8 , pointer :: soil_pulse  (:,:,:) => null()  ! tracking of wet/dry & three types of pulsing (Y&L, 94)
    real*8 , pointer :: base_isop   (:,:,:) => null()  ! baseline emissions for isoprene     (kgC/box/step?)
    real*8 , pointer :: base_monot  (:,:,:) => null()  ! baseline emissions for monoterpenes (kgC/box/step?)
    real*8 , pointer :: xlai        (:,:,:) => null()  ! leaf area index of land type for month #1
    real*8 , pointer :: xlai2       (:,:,:) => null()  ! leaf area index of land type for month #2
    type(t_GmiArrayBundle), pointer :: emissionArray(:)

    real*8 , pointer    :: surf_emiss_out     (:,:,:)     => null()
    real*8 , pointer    :: surf_emiss_out2    (:,:,:)     => null()
    real*8 , pointer    :: emiss_3d_out       (:,:,:,:)   => null()
    real*8 , pointer    :: aerosolEmiss3D     (:,:,:,:)   => null()
    real*8 , pointer    :: aerosolSurfEmiss   (:,:,:)     => null()
                   ! Array storing the surface emission diagnostics for
                   ! aerosols. FSO2, NSO2, DMS, NSO4A, FSO4A, aerosol
                   ! (OC, BC, sea salt) and dust.
    integer, pointer    :: aerosolSurfEmissMap(:)     => null()
                   ! mapping of SurfEmiss number for aerosols to const species #
    character (len=MAX_LENGTH_FILE_NAME) :: MEGAN_infile_name       ! File for MEGAN AEF and LAI
  end type t_Emission

!
! !DESCRIPTION:
! Initialization and run methods for the Emission component.
!
! !AUTHOR:
!
! !REVISION HISTORY:
!   - January 15, 2008 * Eric Nielsen
!     Two-meter air temperature replaces surface air temperature
!
!EOP
!-------------------------------------------------------------------------
  CONTAINS
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: readEmissionResourceFile
!
! !INTERFACE:
!
      subroutine readEmissionResourceFile(self, config, procID, numSpecies, k2, &
                     pr_diag, chem_opt)
!
! !USES:
      use GmiCheckNamelistFile_mod, only : CheckNamelistOptionRange
      use GmiSpeciesRegistry_mod, only : getSpeciesIndex, UNKNOWN_SPECIES
      use GmiStringManipulation_mod, only : constructListNames
!
! !INPUT PARAMETERS:
      logical, intent(in) :: pr_diag
      integer, intent(in) :: numSpecies, procID, chem_opt, k2
!
! !INPUT/OUTPUT PARAMETERS:
      type (t_Emission), intent(inOut) :: self
      type (ESMF_Config), intent(inOut) :: config
!
! !DESCRIPTION:
! Reads in Emission related variables from the resource file.
!
! !LOCAL VARIABLES:
      integer ::  STATUS, RC, ic, ios, n
      character (len=MAX_LENGTH_SPECIES_NAME) :: tempListNames(MAX_NUM_CONST)
      character (len=MAX_STRING_LENGTH      ) :: emissionSpeciesNames
      character (len=MAX_STRING_LENGTH      ) :: emissionDustSpeciesNames
      character (len=MAX_STRING_LENGTH      ) :: emissionAeroSpeciesNames
      character(len=ESMF_MAXSTR) :: IAm, err_msg, sp_name
!.sds.. point emissions variables
      character (len=MAX_STRING_LENGTH      ) :: emissionPointTemp(MAX_NUM_CONST)
      character (len=MAX_STRING_LENGTH      ) :: emissionPointFilenames
!.sds.end
!EOP
!--------------------------------------------------------------------
!BOC
      IAm = "readEmissionResourceFile"

      if (pr_diag) Write(6,*) IAm, 'called by ', procID

!################################
! Begin reading the resource file
!################################

! ---------------------------------     
! emiss_opt 
!   (set up with bit switches, can turn on more than one option by adding
!    the numbers togeether, ie 3 will turn on llnl and harvard emissions)
!   0:  no emissions
!   1:  do LLNL emissions
!   2:  do Harvard emissions
!   3:  do LLNL and Harvard emissions
!   4:  do GSFC emissions (Galactic Cosmic Rays only right now)
!   5:  do LLNL and GSFC emissions
!   6:  do Harvard and GSFC emissions
!   7:  do LLNL and Harvard and GSFC emissions
! ---------------------------------

      call ESMF_ConfigGetAttribute(config, self%emiss_opt, &
                     label   = "emiss_opt:", &
                     default = 0, rc=STATUS )
      VERIFY_(STATUS)
      
! --------------------------------------------
! emiss_in_opt
!   0:  no emissions data
!   1:  set all emiss values to emiss_init_val
!   2:  read in emiss values
!   3:  set in emiss values using constant source terms
! --------------------------------------------

      call ESMF_ConfigGetAttribute(config, self%emiss_in_opt, &
                     label   = "emiss_in_opt:", &
                     default = 0, rc=STATUS )
      VERIFY_(STATUS)

! ------------------------------------
! emiss_conv_flag
!   0:  no conversion performed
!   1:  use emiss_conv_fac
!   2:  convert from kg/km2-hr to kg/s
! ------------------------------------

      call ESMF_ConfigGetAttribute(config, self%emiss_conv_flag, &
                     label   = "emiss_conv_flag:", &
                     default = 0, rc=STATUS )
      VERIFY_(STATUS)

! ------------------------------------------------------------------
! semiss_inchem_flag
!   <0:  if emissions are on, surface emissions will be done in
!        Smvgear chemistry if it is on; outside of chemistry if
!        Smvgear chemistry is off (i.e., "auto set")
!    0:  if emissions are on, surface emissions will be done outside
!        of chemistry
!   >0:  if emissions are on, surface emissions will be done in
!        Smvgear chemistry
! ------------------------------------------------------------------
  
      call ESMF_ConfigGetAttribute(config, self%semiss_inchem_flag, &
                     label   = "semiss_inchem_flag:", &
                     default = -1, rc=STATUS )
      VERIFY_(STATUS)
      
! sets of emissons per year (1 => yearly, 12 => monthly)
    
      call ESMF_ConfigGetAttribute(config, self%emiss_timpyr, &
                     label   = "emiss_timpyr:", &
                     default = 1, rc=STATUS )
      VERIFY_(STATUS)
      
      self%emiss_map(:) =  0  
     
      call rcEsmfReadTable(config, emissionSpeciesNames, &
                          "emissionSpeciesNames::", rc=STATUS)

      call ESMF_ConfigGetAttribute(config, self%emiss_conv_fac, &
                     label   = "emiss_conv_fac:", &
                     default = 1.0d0, rc=STATUS )
      VERIFY_(STATUS)
    
      call ESMF_ConfigGetAttribute(config, self%emiss_init_val, &
                     label   = "emiss_init_val:", &
                     default = 1.0d0, rc=STATUS )
      VERIFY_(STATUS)

      CALL ESMF_ConfigGetAttribute(config, value=self%clim_emiss_by_area, &
                label="clim_emiss_by_area:", DEFAULT=.TRUE., RC=STATUS)
      VERIFY_(STATUS)

! Save the number of emitting layers for each emissionSpeciesName
! ---------------------------------------------------------------

      self%emissionSpeciesLayers(:) = 0

      call rcEsmfReadTable(config, self%emissionSpeciesLayers, &
                          "emissionSpeciesLayers::", rc=STATUS)

! --------------------------------------
! Aerosols and Sulfur from Penner et al.
! --------------------------------------

! emiss_aero_opt   0: for no     aerosol emissions
!                  1: for GMI    aerosol emissions
!                  2: for GOCART aerosol emissions
  
      call ESMF_ConfigGetAttribute(config, self%emiss_aero_opt, &
                     label   = "emiss_aero_opt:", &
                     default = 0, rc=STATUS )
      VERIFY_(STATUS)

! number of aerosol emissions
     
      call ESMF_ConfigGetAttribute(config, self%naero, &
                     label   = "naero:", &
                     default = 1, rc=STATUS )
      VERIFY_(STATUS) 
     
      self%emiss_map_aero(:) = 0    ! map emissions to species #

      call rcEsmfReadTable(config, emissionAeroSpeciesNames, &
                          "emissionAeroSpeciesNames::", rc=STATUS)

! aerosol emission input file

      call ESMF_ConfigGetAttribute(config, self%emiss_aero_infile_name, &
                     label   = "emiss_aero_infile_name:", &
                     default = ' ', rc=STATUS )
      VERIFY_(STATUS)

! emiss_dust_opt   0: for no     dust emissions
!                  1: for GMI    dust emissions
!                  2: for GOCART dust emissions

      call ESMF_ConfigGetAttribute(config, self%emiss_dust_opt, &
                     label   = "emiss_dust_opt:", &
                     default = 0, rc=STATUS )
      VERIFY_(STATUS)

! number of dust bins

      call ESMF_ConfigGetAttribute(config, self%ndust, &
                     label   = "ndust:", &
                     default = 1, rc=STATUS )
      VERIFY_(STATUS)

! number of starting point for new dust emiss.
      
      call ESMF_ConfigGetAttribute(config, self%nst_dust, &
                     label   = "nst_dust:", &
                     default = 1, rc=STATUS )
      VERIFY_(STATUS)
      
! number of dust emissions per emiss. file
     
      call ESMF_ConfigGetAttribute(config, self%nt_dust, &
                     label   = "nt_dust:", &
                     default = 1, rc=STATUS )
      VERIFY_(STATUS)   
     
      self%emiss_map_dust(:) = 0    ! map emissions to species #
      
      emissionDustSpeciesNames = ""
     
      call rcEsmfReadTable(config, emissionDustSpeciesNames, &
                          "emissionDustSpeciesNames::", rc=STATUS)
      
! dust emission input file
      
      call ESMF_ConfigGetAttribute(config, self%emiss_dust_infile_name, &
                     label   = "emiss_dust_infile_name:", &
                     default = ' ', rc=STATUS )
      VERIFY_(STATUS)
     
! -----------------------------------------------------
! Information on scale factors for various NO emissions
! -----------------------------------------------------

! Fossil fuel
      call ESMF_ConfigGetAttribute(config, value=self%doScaleNOffEmiss, &
                label="doScaleNOffEmiss:", default=.false., rc=STATUS)

      call ESMF_ConfigGetAttribute(config, self%scFactorNOff_infile_name, &
                     label   = "scFactorNOff_infile_name:", &
                     default = ' ', rc=STATUS )
      VERIFY_(STATUS)

! Biomass burning
      call ESMF_ConfigGetAttribute(config, value=self%doScaleNObbEmiss, &
                label="doScaleNObbEmiss:", default=.false., rc=STATUS)
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, self%scFactorNObb_infile_name, &
                     label   = "scFactorNObb_infile_name:", &
                     default = ' ', rc=STATUS )
      VERIFY_(STATUS)

! ---------------------------------------------------
! Harvard emissions:  acetone, isoprene, propene, NO.
! ---------------------------------------------------

      self%isop_scale(:) = 1.0d0
     
      call rcEsmfReadTable(config, self%isop_scale, &
                          "isop_scale::", rc=STATUS)
     
!... do Galactic Cosmic Rays source of NOx?
    
!... turn on Galactic Cosmic ray emmission of N and NO

      call ESMF_ConfigGetAttribute(config, value=self%do_gcr, &
                label="do_gcr:", default=.false., __RC__ )

      call ESMF_ConfigGetAttribute(config, self%gcr_infile_name, &
                     label   = "gcr_infile_name:", &
                     default = ' ', __RC__ )

!... do Ship Emission Calculations

      call ESMF_ConfigGetAttribute(config, value=self%do_ShipEmission, &
                label="do_ShipEmission:", default=.false., __RC__ )

!.sds
!... Allow Point Emission Calculations

!      call ESMF_ConfigGetAttribute(config, value=self%do_PointEmission, &
!     &           label="do_PointEmission:", default=.false., __RC__ )
     
!... setup point emissions
      call rcEsmfReadTable(config, emissionPointFilenames &
                          ,"emissionPointFilenames::", __RC__)

!.sds.end
!
!... do MEGAN emission calculations

      call ESMF_ConfigGetAttribute(config, value=self%doMEGANemission, &
                label="doMEGANemission:", default=.false., __RC__ )

      call ESMF_ConfigGetAttribute(config, value=self%doMEGANviaHEMCO, &  
                label="doMEGANviaHEMCO:", default=.false., __RC__ )

      call ESMF_ConfigGetAttribute(config, self%MEGAN_infile_name, &
                     label   = "MEGAN_infile_name:", &
                     default = ' ', __RC__ )

      call ESMF_ConfigGetAttribute(config, self%laiMEGAN_InfileName, &
                     label   = "laiMEGAN_InfileName:", &
                     default = ' ', __RC__ )

      call ESMF_ConfigGetAttribute(config, self%aefMboMEGAN_InfileName, &
                     label   = "aefMboMEGAN_InfileName:", &
                     default = ' ', __RC__ )
     
      call ESMF_ConfigGetAttribute(config, self%aefIsopMEGAN_InfileName, &
                     label   = "aefIsopMEGAN_InfileName:", &
                     default = ' ', __RC__ )
      
      call ESMF_ConfigGetAttribute(config, self%aefOvocMEGAN_InfileName, &
                     label   = "aefOvocMEGAN_InfileName:", &
                     default = ' ', __RC__ )
     
      call ESMF_ConfigGetAttribute(config, self%aefMonotMEGAN_InfileName, &
                     label   = "aefMonotMEGAN_InfileName:", &
                     default = ' ', __RC__ )
     
      call ESMF_ConfigGetAttribute(config, self%soil_infile_name, &
                     label   = "soil_infile_name:", &
                     default = 'soiltype.asc', __RC__ )

      call ESMF_ConfigGetAttribute(config, self%isopconv_infile_name, &
                     label   = "isopconv_infile_name:", &
                     default = 'isopconvtable.asc', __RC__ )

      call ESMF_ConfigGetAttribute(config, self%monotconv_infile_name, &
                     label   = "monotconv_infile_name:", &
                     default = 'monotconvtable.asc', __RC__ )

    ! ---------------------------------------------------------------------------------------
    ! lightning_opt = 0 --> default lightning
    !               = 1 --> variable_lightning (Dale Allen's algorithm)
    !               = 2 --> no lightning
    !               = in options 1 and 2, NO_lgt (in gmi_emiss.F) is made zero
    ! ---------------------------------------------------------------------------------------

      call ESMF_ConfigGetAttribute(config, self%lightning_opt, &
                     label   = "lightning_opt:", &
                     default = 0, rc=STATUS )
      VERIFY_(STATUS)
     
      ! ---------------------------------------------------------------
      ! Check option ranges.  Note that as new options are added, these
      ! range checks will have to be modified.
      ! ---------------------------------------------------------------
    
      call CheckNamelistOptionRange ('emiss_opt'     , self%emiss_opt     , 0, 2)
      call CheckNamelistOptionRange ('emiss_aero_opt', self%emiss_aero_opt, 0, 2)
      call CheckNamelistOptionRange ('emiss_dust_opt', self%emiss_dust_opt, 0, 2)
      call CheckNamelistOptionRange ('emiss_in_opt'  , self%emiss_in_opt  , 0, 3)
      call CheckNamelistOptionRange ('lightning_opt' , self%lightning_opt , 0, 2)

      if (self%emiss_opt == 0) then
         self%emiss_in_opt    = 0
         self%emiss_conv_flag = 0
         self%emiss_map(:)    = 0
         self%num_emiss = 0
!.sds
         self%num_point_emiss = 0
         self%num_point_start = 1
!.sds.end
      else

         ! Set the initial value of the list
         tempListNames(:) = ''

         ! Construct the list of names using the long string
         call constructListNames(tempListNames, emissionSpeciesNames)

         self%num_emiss = Count (tempListNames(:) /= '')
!
!.sds.. Is DMS in mechanism? If so, need to turn on emissions
         ic = getSpeciesIndex('DMS',NOSTOP=.true.)
         if(ic .gt. 0) then
           self%GMIDMSEmissIndex = self%num_emiss + 1
           self%num_emiss = self%num_emiss + 1
           self%emiss_map(self%num_emiss) = ic
           self%emissionSpeciesLayers(self%num_emiss) = 1
!... Add another entry to tempListNames
           tempListNames(self%num_emiss) = 'DMS'
         else
           self%GMIDMSEmissIndex = 0
         endif
!.sds.end
! MEM
         IF ( self%do_ShipEmission ) THEN
           ! Add 2 special entries
           tempListNames( self%num_emiss + 1 ) = '*shipO3*'
           tempListNames( self%num_emiss + 2 ) = '*shipHNO3*'
           self%num_emiss = self%num_emiss + 2
         END IF

         ! The species name precedes the underscore. In GEOS-5 we need to keep
         ! the entire emissionSpeciesName for use in the initialize method.
         ! -------------------------------------------------------------------
         IF (self%num_emiss > 0) THEN
            DO ic = 1, self%num_emiss
               sp_name = tempListNames(ic)
               self%emissionSpeciesNames(ic) = TRIM(sp_name)
               ios = INDEX(sp_name,'_')
               IF(ios > 0) sp_name = sp_name(1:ios-1)

               IF      ( TRIM(sp_name) ==        '*shipO3*'   ) THEN
                 self%emiss_map(ic) = getSpeciesIndex('O3')
                ELSE IF ( TRIM(sp_name) ==       '*shipHNO3*' ) THEN
                 self%emiss_map(ic) = getSpeciesIndex('HNO3')
                ELSE
                 self%emiss_map(ic) = getSpeciesIndex(TRIM(sp_name))
                ENDIF

            END DO
         END IF
!.sds
!... point emissions filenames and species
         emissionPointTemp(:) = ''
         call constructListNames(emissionPointTemp, emissionPointFilenames)
         self%num_point_emiss = Count (emissionPointTemp(:) .ne. '')
!... get point emission prelim stuff set up
         self%num_point_start = self%num_emiss + 1
         if(self%num_point_emiss .gt. 0) then
           do n = 1, self%num_point_emiss
             sp_name = emissionPointTemp(n)
             ios = INDEX(sp_name,':',.true.)
             if(ios .gt. 0) then
               self%emissionPointFilenames(n) = sp_name(ios+1:)
               sp_name = sp_name(1:ios-1)
               ios = INDEX(sp_name,':',.true.)
               self%num_point_type(n) = sp_name(ios+1:)
               sp_name = sp_name(1:ios-1)
             else
               err_msg = 'emissionPointFilenames problem in rc File.'
               call GmiPrintError (err_msg, .true., ios, self%num_point_emiss, &
                 self%num_emiss, 0, 0.0d0, 0.0d0)
             endif
!... get index number of species
             ic = getSpeciesIndex(sp_name,NOSTOP=.true.)
!... add to emissionSpeciesNames, etc
             self%num_emiss = self%num_emiss + 1
             self%emissionSpeciesNames(self%num_emiss) = TRIM(sp_name)
             self%emiss_map(self%num_emiss) = ic
             self%emissionSpeciesLayers(self%num_emiss) = k2
           enddo
         endif
!.sds.end
         self%num_emiss = count( self%emiss_map(:) > 0 )
!
         IF(MAPL_AM_I_ROOT()) THEN
           print '('' '')'
           print '(''Num_emiss='',i5)', self%num_emiss
           print '(''Num_point_emiss='',i5)', self%num_point_emiss
           print '('' '')'
           print '(''EmissNo  InputName                  SpeciesNo   NumOfLevels'')'
            do ic = 1, self%num_point_start-1
             print '(i7,2x,a24,2x,i10,2x,i12)', ic, self%emissionSpeciesNames(ic) &
              , self%emiss_map(ic), self%emissionSpeciesLayers(ic)
           enddo
            do ic = self%num_point_start,self%num_emiss
             print '(i7,2x,a10,1x,a18,1x,i6,2x,i12)', ic, self%num_point_type(ic-self%num_point_start+1) &
              , self%emissionSpeciesNames(ic), self%emiss_map(ic), self%emissionSpeciesLayers(ic)
           enddo
         endif
!
         if (self%emiss_dust_opt > 0) then
            ! Set the initial value of the list
            tempListNames(:) = ''

            ! Construct the list of names using the long string
            call constructListNames(tempListNames, emissionDustSpeciesNames)
         
            do ic = 1, self%ndust
               self%emiss_map_dust(ic) = getSpeciesIndex(tempListNames(ic))
            end do
         end if
         
         if (self%emiss_aero_opt > 0) then
            ! Set the initial value of the list
            tempListNames(:) = ''
         
            ! Construct the list of names using the long string
            call constructListNames(tempListNames, emissionAeroSpeciesNames)

            do ic = 1, self%naero
               self%emiss_map_aero(ic) = getSpeciesIndex(tempListNames(ic))
            end do
         end if 
         
      endif

      if (self%semiss_inchem_flag < 0) then  ! "auto set"
         if ((self%emiss_opt /= 0) .and. (chem_opt == 2)) then
            self%do_semiss_inchem = .true.
         else
            self%do_semiss_inchem = .false.
         end if
      else if (self%semiss_inchem_flag == 0) then
         self%do_semiss_inchem   = .false.
      else if (self%semiss_inchem_flag  > 0) then
         self%do_semiss_inchem   = .true.
      end if

      !###############
      ! Error Checking
      !###############

      if ((numSpecies > MAX_NUM_CONST) .or.  &
          (self%num_emiss   > MAX_NUM_CONST)) then
         err_msg = 'MAX_NUM_CONST problem in rc File.'
         call GmiPrintError (err_msg, .true., 2, numSpecies, &
                self%num_emiss, 0, 0.0d0, 0.0d0)
      end if

      if ((self%emiss_opt /= 0) .and. (self%num_emiss == 0)) then
        err_msg = 'emiss_opt/num_emiss problem.'
        call GmiPrintError  &
         (err_msg, .true., 2, self%emiss_opt, self%num_emiss, 0, 0.0d0, 0.0d0)
      end if

      if ((self%emiss_in_opt == 0) .and. (self%emiss_opt /= 0)) then
        err_msg = 'emiss_in_opt/emiss_opt problem.'
        call GmiPrintError  &
         (err_msg, .true., 2, self%emiss_in_opt, self%emiss_opt, 0, 0.0d0, 0.0d0)
      end if
         
      if ((self%emiss_timpyr /= 1) .and.  &
         (self%emiss_timpyr /= MONTHS_PER_YEAR)) then
        err_msg = 'emiss_timpyr range problem.'
        call GmiPrintError  &
         (err_msg, .true., 1, self%emiss_timpyr, 0, 0, 0.0d0, 0.0d0)
      end if
      
      if ((self%emiss_in_opt /= 2) .and.  &
         (self%emiss_timpyr == MONTHS_PER_YEAR)) then
        err_msg = 'emiss_in_opt/emiss_timpyr problem.'
        call GmiPrintError (err_msg, .true., 2, self%emiss_in_opt,             &
                           self%emiss_timpyr, 0, 0.0d0, 0.0d0)
      end if     
      
  return
      
  end subroutine readEmissionResourceFile
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: initReadEmission
!
! !INTERFACE:
!
  subroutine initReadEmission (self, gmiClock, gmiGrid, mcor, loc_proc, &
                               pr_diag, rc )
!
! !USES:
  use gcr_mod             , only : READ_GCR_FILE

! !INPUT PARAMETERS:
  logical,           intent(in) :: pr_diag
  integer,           intent(in) :: loc_proc
  real*8,            intent(in) :: mcor(:,:)
  type (t_gmiGrid),  intent(in) :: gmiGrid
  type (t_GmiClock), intent(in) :: gmiClock
!
! !INPUT/OUTPUT PARAMETERS:
  type (t_Emission), intent(inOut) :: self
  integer,           intent(out)   :: rc

!
! !DESCRIPTION:
! This routines reads in (daily or monthly) emission related files.
!
! !LOCAL VARIABLES:
  character (len=75) :: err_msg
  integer       :: il, ij, num_emiss, ic
  real*8        :: units_fac
  integer       :: nymd, ydummy, thisDay, thisMonth, thisDate, ddummy
  integer       :: i1, i2, ju1, j2, k1, k2
  integer       :: i1_gl, i2_gl, ju1_gl, j2_gl, ilong, ilat, ivert
  logical       :: rootProc
!
!EOP
!----------------------------------------------------------------------
!BOC

  call Get_i1    (gmiGrid, i1)
  call Get_i2    (gmiGrid, i2)
  call Get_ju1   (gmiGrid, ju1)
  call Get_j2    (gmiGrid, j2)
  call Get_k1    (gmiGrid, k1)
  call Get_k2    (gmiGrid, k2)
  call Get_i1_gl (gmiGrid, i1_gl)
  call Get_ju1_gl(gmiGrid, ju1_gl)
  call Get_ilong (gmiGrid, ilong)
  call Get_ilat  (gmiGrid, ilat)

  rootProc = MAPL_AM_I_ROOT()
  rc = 0 

  if (self%emiss_dust_opt == 1 ) then  ! GMI sulfurous dust emissions
    call InitEmissDust (self, mcor, loc_proc, rootProc, &
                        i1, i2, ju1, j2, i1_gl, ju1_gl, ilong, ilat)
  end if

  if (self%emiss_aero_opt /= 0 ) then
    call InitEmissAero (self, mcor, loc_proc, rootProc, &
                        i1, i2, ju1, j2, i1_gl, ju1_gl, ilong, ilat)
  end if

  call Get_curGmiDate(gmiClock, nymd)

  if (self%do_gcr) then
    CALL READ_GCR_FILE( self%gcr_infile_name, rc )   ! Store data in the GCR module
  end if
  
  return

  end subroutine initReadEmission
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: InitializeEmission
!
! !INTERFACE:
!
  subroutine InitializeEmission (self, SpeciesConcentration, gmiGrid, config,              &
                       mcor, ihno3_num, io3_num, numSpecies, loc_proc, rootProc,           &
                       chem_opt, trans_opt, pr_diag, pr_const, pr_surf_emiss, pr_emiss_3d, &
                       tdt)
!
! !USES:
  use ReadOtherEmissionData_mod, only : readLightData
  use ReadOtherEmissionData_mod, only : readIsopreneConvertData, readSoilData
  use ReadOtherEmissionData_mod, only : readMonoterpeneConvertData
  use GocartDerivedVariables_mod, only : AllocateGocartDerivedVars
  use GmiSeaSaltMethod_mod      , only : InitializationSeaSalt
  use GmiDustMethod_mod         , only : InitializationDust

!
! !INPUT PARAMETERS:
  logical,           intent(in) :: pr_diag, pr_const, pr_surf_emiss, pr_emiss_3d
  logical,           intent(in) :: rootProc
  integer          , intent(in) :: numSpecies, io3_num, ihno3_num
  integer          , intent(in) :: loc_proc, chem_opt, trans_opt
  real*8           , intent(in) :: mcor(:, :)
  type (t_gmiGrid ), intent(in) :: gmiGrid 
  type (t_SpeciesConcentration), intent(in) :: SpeciesConcentration 
  real             , intent(in) :: tdt
!
! !INPUT/OUTPUT VARIABLES:
  type (ESMF_Config), intent(inOut) :: config
  type (t_Emission) , intent(inOut) :: self
!
! !DESCRIPTION:
! Initialize the Emission component.
!
! !LOCAL VARIABLES:
  integer :: ix, ic, id
  integer :: i1, i2, ju1, j2, k1, k2
  integer :: i1_gl, i2_gl, ju1_gl, j2_gl
  integer :: ilong, ilat, ivert
  integer :: tracer_opt
  integer :: badIndex = -9999
  real*8  :: tr_source_land, tr_source_ocean
  real*8  :: tdt_days, rsecpdy
!EOP
!------------------------------------------------------------------------------
!BOC
  call Get_i1    (gmiGrid, i1)
  call Get_i2    (gmiGrid, i2)
  call Get_ju1   (gmiGrid, ju1)
  call Get_j2    (gmiGrid, j2)
  call Get_k1    (gmiGrid, k1)
  call Get_k2    (gmiGrid, k2)
  call Get_i1_gl (gmiGrid, i1_gl)
  call Get_i2_gl (gmiGrid, i2_gl)
  call Get_ju1_gl(gmiGrid, ju1_gl)
  call Get_j2_gl (gmiGrid, j2_gl )
  call Get_ilong (gmiGrid, ilong )
  call Get_ilat  (gmiGrid, ilat  )
  call Get_ivert (gmiGrid, ivert )

  call Get_tracer_opt     (SpeciesConcentration, tracer_opt   )
  call Get_tr_source_land (SpeciesConcentration, tr_source_land )
  call Get_tr_source_ocean(SpeciesConcentration, tr_source_ocean)

  call readEmissionResourceFile(self, config, loc_proc, numSpecies, k2, pr_diag, chem_opt)
  self%idaySoilType = badIndex
  self%firstBiogenicBase = .TRUE.

  if ((tracer_opt == 0) .and. btest(self%emiss_opt,1)) then
    allocate (self%index_soil(2, i1:i2, ju1:j2))
    allocate (self%soil_pulse(NPULSE+1, i1:i2, ju1:j2))
    self%soil_pulse(1:NPULSE+1, i1:i2, ju1:j2) = 0.00D+00

    rsecpdy  = SECPDY
    tdt_days = tdt / rsecpdy

    self%exp_fac(1) = Exp (-PULSE_DECAY(1) * tdt_days)
    self%exp_fac(2) = Exp (-PULSE_DECAY(2) * tdt_days)
    self%exp_fac(3) = Exp (-PULSE_DECAY(3) * tdt_days)

    call Allocate_ireg       (self, i1, i2, ju1, j2)
    call Allocate_iuse       (self, i1, i2, ju1, j2)
    call Allocate_iland      (self, i1, i2, ju1, j2)
    call Allocate_xlai       (self, i1, i2, ju1, j2)
    call Allocate_xlai2      (self, i1, i2, ju1, j2)

    IF(self%doMEGANemission) THEN
      CALL Allocate_isoLaiCurr  (self, i1, i2, ju1, j2)
      CALL Allocate_isoLaiPrev  (self, i1, i2, ju1, j2)
      CALL Allocate_isoLaiNext  (self, i1, i2, ju1, j2)
      CALL Allocate_isoLaiYear  (self, i1, i2, ju1, j2)
      CALL Allocate_aefMbo      (self, i1, i2, ju1, j2)
      CALL Allocate_aefIsop     (self, i1, i2, ju1, j2)
      CALL Allocate_aefOvoc     (self, i1, i2, ju1, j2)
      CALL Allocate_aefMonot    (self, i1, i2, ju1, j2)
    ELSE
      CALL Allocate_base_isop  (self, i1, i2, ju1, j2)
      CALL Allocate_base_monot (self, i1, i2, ju1, j2)

      CALL readIsopreneConvertData &
     	      (self%isopconv_infile_name, self%convert_isop, pr_diag, rootProc)

      CALL readMonoterpeneConvertData &
              (self%monotconv_infile_name, self%convert_monot, pr_diag, rootProc)
    END IF

    call readSoilData &
             (self%soil_infile_name, self%ncon_soil, pr_diag, rootProc)
  end if

  IF(rootProc) THEN
    WRITE(6,*) ' '
  END IF

  if (self%emiss_in_opt /= 0) then
    call Allocate_emissionArray (self, i1, i2, ju1, j2, k1, k2)
  endif

! Diagnostic variables: surface and 3D emissions

  if (pr_const) then
    if (pr_surf_emiss) then
      call Allocate_surf_emiss_out (self, i1, i2, ju1, j2, numSpecies)
      call Allocate_surf_emiss_out2(self, i1, i2, ju1, j2)
   
!        if ((self%emiss_aero_opt > 0) .or. (self%emiss_dust_opt > 0)) then
!           call Allocate_aerosolSurfEmissMap(self)
!           call setAerosolSurfEmissMap(self%aerosolSurfEmissMap, &
!                                       self%emiss_map_aero, &
!                                       self%emiss_map_dust, &
!                                       self%naero, &
!                                       self%ndust)
!           call Allocate_aerosolSurfEmiss(self, i1, i2, ju1, j2)
!        endif
    endif
   
    if (pr_emiss_3d) then
      call Allocate_emiss_3d_out (self, i1, i2, ju1, j2, k1, k2, numSpecies)
      if ((self%emiss_aero_opt > 0) .or. (self%emiss_dust_opt > 0)) then
        call Allocate_aerosolEmiss3D(self, i1, i2, ju1, j2, k1, k2)
      end if
    end if
  end if

  if (self%emiss_opt == 2) then
    call Allocate_emiss_monot (self, i1, i2, ju1, j2)
    call Allocate_emiss_isop  (self, i1, i2, ju1, j2)
    call Allocate_emiss_nox   (self, i1, i2, ju1, j2)
  endif

  if (self%do_ShipEmission) then
    call Allocate_emiss_o3(self, i1, i2, ju1, j2)
    call Allocate_emiss_hno3 (self, i1, i2, ju1, j2)
     
    self%ship_o3_index   = badIndex
    self%ship_hno3_index = badIndex
     
    do ix = 1, numSpecies
      if ( self%emissionSpeciesNames(ix) == '*shipO3*'   )  self%ship_o3_index   = ix
      if ( self%emissionSpeciesNames(ix) == '*shipHNO3*' )  self%ship_hno3_index = ix
    end do

    if (self%ship_o3_index   == badIndex) stop "shipO3   is not in emissionArray"
    if (self%ship_hno3_index == badIndex) stop "shipHNO3 is not in emissionArray"

  endif

  if (self%emiss_aero_opt /= 0) then
    call Allocate_emissAero   (self, i1, i2, ju1, j2)
    call Allocate_emissAero_t (self, i1, i2, ju1, j2)
  endif

  if (self%emiss_dust_opt /= 0) then
    call Allocate_emissDust   (self, i1, i2, ju1, j2)
    call Allocate_emissDust_t (self, i1, i2, ju1, j2)
  endif

  if (self%lightning_opt == 1) then
     call Allocate_lightning_NO (self, i1, i2, ju1, j2, k1, k2)
  endif

  if ((self%emiss_dust_opt == 2) .or. (self%emiss_aero_opt == 2))  then
    call AllocateGocartDerivedVars(i1, i2, ju1, j2, k1, k2)
    if (self%emiss_aero_opt == 2) call InitializationSeaSalt()
    if (self%emiss_dust_opt == 2) call InitializationDust   ()
  end if

! Set the intial emission array

  call InitEmiss (self, mcor, tr_source_ocean, tr_source_land, &
                        loc_proc, rootProc, pr_diag, trans_opt, &
                        i1, i2, ju1, j2, k1, k2, i1_gl, i2_gl, ju1_gl, j2_gl)

  return

  end subroutine InitializeEmission
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINES: RunEmission
!
! !INTERFACE:
!
      subroutine RunEmission(self, SpeciesConcentration, gmiClock, gmiGrid,    &
                   loc_proc, cosSolarZenithAngle, latdeg, mcor, mass,         &
                   lwis_flags, radswg, TwoMeter_air_temp, surf_rough,         &
                   con_precip, tot_precip, ustar, fracCloudCover, kel, pbl,   &
                   cmf, press3c, press3e, pctm2, u10m, v10m,                  &
                   gwet, gridBoxHeight, mw, IBOC, IBBC, INOC, IFOC, IFBC,     &
                   ISSLT1, ISSLT2, ISSLT3, ISSLT4, IFSO2, INSO2, INDMS, IAN,  &
                   IMGAS, INO, iisoprene_num, ino_num, ico_num, ipropene_num, &
                   ihno3_num, io3_num,  numSpecies, pardif, pardir, T_15_AVG, &
                   met_opt, chem_opt, trans_opt, do_aerocom, do_drydep,       &
                   pr_diag, pr_const, pr_surf_emiss, pr_emiss_3d,             &
                   metdata_name_org, metdata_name_model, tdt4, mixPBL,        &
                   light_NO_prod )
!
! !USES:
      use GmiTimeControl_mod, only : Get_numTimeSteps
      use GocartDerivedVariables_mod, only : SetGocartDerivedVars
!
! !INPUT PARAMETERS:
      type(t_GmiClock), intent(in) :: gmiClock
      type (t_gmiGrid), intent(in) :: gmiGrid

      integer, intent(in) :: met_opt, chem_opt, trans_opt
      logical, intent(in) :: do_aerocom, do_drydep
      logical, intent(in) :: pr_diag, pr_const, pr_surf_emiss, pr_emiss_3d
      integer, intent(in) :: IBOC, IBBC, INOC, IFOC, IFBC, ISSLT1, ISSLT2, ISSLT3, ISSLT4
      integer, intent(in) :: IFSO2, INSO2, INDMS, IAN, IMGAS, INO, iisoprene_num
      integer, intent(in) :: ino_num, ico_num, ipropene_num, ihno3_num, io3_num
      integer, intent(in) :: loc_proc
      integer, intent(in) :: numSpecies
      real*8 , intent(in) :: mw(numSpecies)
      real*8 , intent(in) :: cosSolarZenithAngle(:, :), latdeg(:,:), mass(:,:,:)
      real*8 , intent(in) :: mcor(:, :)

      real*8 , intent(in) :: pctm2  (:, :), press3c(:, :, :), press3e(:, :, :)

      integer, intent(in) :: lwis_flags(:, :)            ! 0=water; 1=land; 2=ice; 3=snow

      real*8 , intent(in) :: radswg(:, :), TwoMeter_air_temp(:, :), surf_rough(:, :)
      real*8 , intent(in) :: pardif(:, :), pardir(:, :), T_15_AVG(:, :)
      real*8 , intent(in) :: con_precip(:, :), tot_precip(:, :), ustar(:, :)
      real*8 , intent(in) :: kel(:, :, :), pbl(:, :), u10m(:,:), v10m(:,:), gwet(:,:)
      real*8 , intent(in) :: fracCloudCover (:,:) , cmf(:,:,:)
      real*8 , intent(in) :: gridBoxHeight (:,:,:)

      character (len=MAX_LENGTH_MET_NAME), intent(in) :: metdata_name_org
      character (len=MAX_LENGTH_MET_NAME), intent(in) :: metdata_name_model
      real   , intent(in) :: tdt4
      logical, intent(in) :: mixPBL   ! whether to explicitly distribute
                                      ! aerosol emissions within the PBL
      real*4 , intent(in) :: light_NO_prod(:, :, :)   ! (m-3 s-1) TOP-DOWN !!
!
! !INPUT/OUTPUT VARIABLES:
      type (t_Emission),            intent(inOut) :: self
      type(t_SpeciesConcentration), intent(inOut) :: SpeciesConcentration
!
! !DESCRIPTION:
! Run method for the Emission component.
!
! !LOCAL VARIABLES:
      real*4 , allocatable :: productionNO(:,:,:)
      type(t_GmiArrayBundle), pointer :: concentration(:)

      character (len=128) err_msg
      integer       :: ydummy, thisDay, thisMonth, thisDate, ddummy, ic, curRecord
      logical       :: rootProc
      integer       :: nymd, num_time_steps, ndt
      real*8        :: tdt8
      integer       :: i1, i2, ju1, j2, k, k1, k2, ilo, ihi, julo, jhi
      integer       :: i1_gl, i2_gl, ju1_gl, j2_gl, ilong, ilat, ivert
      real*8        :: tr_source_ocean, tr_source_land
      integer       :: rc, tracer_opt

      real, parameter :: MWT_NO       = 30.0064     ! molecular weight of NO

!
!EOP
!------------------------------------------------------------------------------
!BOC

!... Get the GMI grid information
      call Get_i1    (gmiGrid, i1   )
      call Get_i2    (gmiGrid, i2   )
      call Get_ju1   (gmiGrid, ju1  )
      call Get_j2    (gmiGrid, j2   )
      call Get_k1    (gmiGrid, k1   )
      call Get_k2    (gmiGrid, k2   )
      call Get_i1_gl (gmiGrid, i1_gl)
      call Get_i2_gl (gmiGrid, i2_gl)
      call Get_ju1_gl(gmiGrid, ju1_gl)
      call Get_j2_gl (gmiGrid, j2_gl )
      call Get_ilo   (gmiGrid, ilo  )
      call Get_ihi   (gmiGrid, ihi  )
      call Get_julo  (gmiGrid, julo )
      call Get_jhi   (gmiGrid, jhi  )
      call Get_ilong (gmiGrid, ilong)
      call Get_ilat  (gmiGrid, ilat )
      call Get_ivert (gmiGrid, ivert)

!... Obtain model time information
      call Get_curGmiDate  (gmiClock, nymd          )
      call Get_numTimeSteps(gmiClock, num_time_steps)

      ndt = Nint (tdt4)
      tdt8 = tdt4

      call Get_tracer_opt     (SpeciesConcentration, tracer_opt     )
      call Get_tr_source_land (SpeciesConcentration, tr_source_land )
      call Get_tr_source_ocean(SpeciesConcentration, tr_source_ocean)
      call Get_concentration  (SpeciesConcentration, concentration  )

     if (self%emiss_opt == 2) then
        self%emiss_monot = 0.0d0
        self%emiss_nox   = 0.0d0
        IF ( .not. self%doMEGANviaHEMCO ) THEN 
           self%emiss_isop  = 0.0d0 
        END IF
     end if

!... For GOCART emission
     if ((self%emiss_dust_opt == 2) .or. (self%emiss_aero_opt == 2))  then
        call SetGocartDerivedVars (u10m, v10m, gwet, press3c, &
                         kel, mass, mcor, gridBoxHeight,     &
                         i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi)
     end if

!... Emission control routines
      rootProc = MAPL_AM_I_ROOT()
      call Update_Emiss (gmiGrid, self%idaySoilType, self%firstBiogenicBase,   &
                 lwis_flags, cosSolarZenithAngle, latdeg,                      &
                 mcor, self%emiss_isop, self%emiss_monot, self%emiss_nox,      &
                 self%do_ShipEmission, self%emiss_hno3, self%emiss_o3,         &
                 ihno3_num, io3_num, radswg, TwoMeter_air_temp, surf_rough,    &
                 con_precip, tot_precip, ustar, mass, fracCloudCover, kel,     &
                 self%surf_emiss_out, self%surf_emiss_out2, self%emiss_3d_out, &
                 self%aerosolEmiss3D, self%aerosolSurfEmiss,                   &
                 self%aerosolSurfEmissMap, concentration, self%emissionArray,  &
                 self%emissDust_t, self%emissDust, self%emissAero_t,           &
                 self%emissAero, pbl, gridBoxHeight, self%index_soil,          &
                 self%ncon_soil, self%soil_fert, self%soil_precip,             &
                 self%soil_pulse, self%ireg, self%iland, self%iuse,            &
                 self%convert_isop, self%convert_monot, self%coeff_isop,       &
                 self%base_isop, self%base_monot, self%xlai, IBOC, IBBC, INOC, &
                 IFOC, IFBC, ISSLT1, ISSLT2, ISSLT3, ISSLT4, IFSO2, INSO2,     &
                 INDMS, IAN, IMGAS, INO, iisoprene_num, ino_num, ico_num,      &
                 ipropene_num, pr_surf_emiss, pr_emiss_3d, pr_diag, loc_proc,  &
                 rootProc, met_opt, self%emiss_opt, chem_opt, trans_opt,       &
                 self%emiss_aero_opt, self%emiss_dust_opt, do_aerocom,         &
                 self%do_semiss_inchem, self%do_gcr, do_drydep, self%emiss_map,&
                 self%emiss_map_dust, self%emiss_map_aero, self%ndust,         &
                 self%nst_dust, self%nt_dust, self%naero, nymd, num_time_steps,&
                 mw, tdt8, ndt, self%emiss_timpyr, self%num_emiss,             &
                 self%isop_scale, i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi,&
                 i1_gl, i2_gl, ju1_gl, j2_gl, ilong, numSpecies,               &
                 self%doMEGANemission, self%doMEGANviaHEMCO, self%aefIsop,     &
                 self%aefMbo, self%aefMonot, self%isoLaiPrev, self%isoLaiCurr, &
                 self%isoLaiNext, pardif, pardir, T_15_AVG,                    &
                 self%emissionSpeciesLayers, self%exp_fac, mixPBL, press3c )
! NO production from lightning, with flashRate imported from MOIST
! ----------------------------------------------------------------
      IF(self%lightning_opt == 1) THEN

        ALLOCATE(productionNO(i1:i2, ju1:j2, k1:k2))

! Flip from Top-Down to Bottom-Up:

        productionNO(:,:,k1:k2) = light_NO_prod(:,:,k2:k1:-1)

!       NOTE: MAPL_AVOGAD is [molec/kmol]
        self%lightning_NO(:,:,:) = DBLE(productionNO(:,:,:)*MWT_NO/MAPL_AVOGAD)

! Convert production rate, currently m^{-3} s^{-1}, to volume
! mixing ratio (mole fraction) tendency. Note that mass has units of kg.
! ----------------------------------------------------------------------
        DO k = k1,k2
!         NOTE: MAPL_AVOGAD is [molec/kmol]
          productionNO(:,:,k) = productionNO(:,:,k)*MWTAIR*mcor(:,:)* &
                               gridBoxHeight(:,:,k)/(MAPL_AVOGAD*mass(:,:,k))
        END DO

! Update NO mole fraction and save the nitrogen density tendency for export
! -------------------------------------------------------------------------
     	concentration(ino_num)%pArray3D(:,:,:) = &
     	 concentration(ino_num)%pArray3D(:,:,:) + productionNO(:,:,:) * tdt8

     	DEALLOCATE(productionNO)

     END IF

     CALL Set_concentration(SpeciesConcentration, concentration)

     RETURN

     END SUBROUTINE RunEmission
!EOC
!-------------------------------------------------------------------------

  subroutine FinalizeEmission (self)

  type (t_Emission)   , intent(inOut) :: self

  PRINT*,'  Finalize Emission'

  return

  end subroutine FinalizeEmission

!-------------------------------------------------------------------------

  subroutine InitEmiss (self, mcor, tr_source_ocean, tr_source_land, &
                        loc_proc, rootProc, pr_diag, trans_opt, &
                        i1, i2, ju1, j2, k1, k2, i1_gl, i2_gl, ju1_gl, j2_gl)


  type (t_Emission)  , intent(inOut) :: self
  integer, intent(in) :: i1, i2, ju1, j2, k1, k2
  integer, intent(in) :: i1_gl, i2_gl, ju1_gl, j2_gl
  logical, intent(in) :: rootProc, pr_diag
  integer, intent(in) :: loc_proc, trans_opt
  real*8,  intent(in) :: mcor(i1:i2, ju1:j2)
  real*8 , intent(in) :: tr_source_ocean, tr_source_land

! ----------------------
! Variable declarations.
! ----------------------

  character (len=75) :: err_msg
  integer :: il, ij, num_emiss, ic
  real*8  :: tr_source
  real*8  :: units_fac

! ----------------
! Begin execution.
! ----------------

  if (rootProc .AND. pr_diag) then
    Write (6,*) 'InitEmiss called by ', loc_proc
  end if

  num_emiss = self%num_emiss

  if (self%emiss_in_opt == 1) then

    do ic = 1, num_emiss
      self%emissionArray(ic)%pArray3D(:,:,:) = self%emiss_init_val
    end do

  end if

! ----------------------------
! Convert emiss, if necessary.
! ----------------------------

  if (self%emiss_conv_flag /= 0) then

    if (self%emiss_conv_flag == 1) then

      do ic = 1, num_emiss
        self%emissionArray(ic)%pArray3D(:,:,:) = &
          self%emissionArray(ic)%pArray3D(:,:,:) * self%emiss_conv_fac
      end do

    else if (self%emiss_conv_flag == 2) then

      if (trans_opt == 2) then
        IF(rootProc) THEN
          WRITE(6,*) "GMI initEmis: emiss_conv_flag=2 and trans_opt=2 not allowed."
        END IF
        STOP
      end if

      units_fac = (1.0d0 / SECPHR) * (KMPM * KMPM)

      do ic = 1, num_emiss
        do ij = ju1, j2
          do il = i1, i2
            self%emissionArray(ic)%pArray3D(il,ij,:) =  &
              self%emissionArray(ic)%pArray3D(il,ij,:) * mcor(il,ij) * units_fac
          end do
        end do
      end do

    end if

  end if

  do ic = 1, num_emiss
    call CheckRange3d  &
      ('emissionArray', loc_proc, i1, i2, ju1, j2, k1, k2, &
        self%emissionArray(ic)%pArray3D(:,:,:), -1.0d20, 1.0d20)
  end do
  
  return

  end subroutine InitEmiss

!-------------------------------------------------------------------------
! This routine sets the aerosol emissions.
!-------------------------------------------------------------------------

      subroutine InitEmissAero (self, mcor, loc_proc, rootProc, &
                                 i1, i2, ju1, j2, i1_gl, ju1_gl, ilong, ilat)

!     ----------------------
!     Argument declarations.
!     ----------------------

      type (t_Emission), intent(inOut) :: self
      integer          , intent(in   ) :: i1, i2, ju1, j2, i1_gl, ju1_gl, ilong, ilat
      integer          , intent(in   ) :: loc_proc
      logical          , intent(in   ) :: rootProc
      real*8           , intent(in   ) :: mcor(i1:i2, ju1:j2)

!     ----------------------
!     Variable declarations.
!     ----------------------

      integer :: il, ij

!     ----------------
!     Begin execution.
!     ----------------

      if (self%emiss_aero_opt /= 0) then
      
       IF(rootProc) THEN
        WRITE(6,*) 'GMI InitEmissAero: Must use Chem_UtilMPread from run method'
       END IF
       STOP

!       -----------------------------------------------------------------
!       Other aero (carbon & sslt) emissions in kg/m^2/s, change to kg/s.
!       -----------------------------------------------------------------

        do ij = ju1, j2
          do il = i1, i2
            self%emissAero_t(il,ij,:,:) = self%emissAero_t(il,ij,:,:) * mcor(il,ij)
          end do
        end do

      end if

      return

      end subroutine InitEmissAero

!-------------------------------------------------------------------------
! This routine sets the dust emissions.
!-------------------------------------------------------------------------

      subroutine InitEmissDust  (self, mcor, loc_proc, rootProc, &
                                i1, i2, ju1, j2, i1_gl, ju1_gl, ilong, ilat)

!     ----------------------
!     Argument declarations.
!     ----------------------

      type (t_Emission)  , intent(inOut) :: self
      integer, intent(in) :: i1, i2, ju1, j2
      integer, intent(in) :: i1_gl, ju1_gl, ilong, ilat
      integer, intent(in) :: loc_proc
      logical, intent(in) :: rootProc
      real*8 , intent(in) :: mcor(i1:i2, ju1:j2)

!     ----------------------
!     Variable declarations.
!     ----------------------

      integer :: il, ij

!     ----------------
!     Begin execution.
!     ----------------

      if (self%emiss_dust_opt == 1) then  ! GMI dust emissions

       IF(rootProc) THEN
        WRITE(6,*) 'GMI InitEmissDust: Must use Chem_UtilMPread from run method'
       END IF
       STOP


!       -------------------------------------------
!       Dust emissions in kg/m^2/s, change to kg/s.
!       -------------------------------------------

        do ij = ju1, j2
          do il = i1, i2
            self%emissDust_t(il,ij,:,:) = self%emissDust_t(il,ij,:,:) * mcor(il,ij)
          end do
        end do

      end if

      return

      end subroutine InitEmissDust

!-------------------------------------------------------------------------
  subroutine Allocate_emissionArray (self, i1, i2, ju1, j2, k1, k2)
    integer          , intent(in   ) :: i1, i2, ju1, j2, k1, k2
    type (t_Emission), intent(inOut) :: self
    integer                          :: i, num_emiss
    num_emiss = self%num_emiss
    Allocate(self%emissionArray(num_emiss))
    do i = 1, num_emiss
       Allocate(self%emissionArray(i)%pArray3D(i1:i2, ju1:j2, k1:k2))
       self%emissionArray(i)%pArray3D(i1:i2, ju1:j2, k1:k2) = 0.0d0
    end do
    return
  end subroutine Allocate_emissionArray
!-------------------------------------------------------------------------
!!   subroutine Get_emissionArray (self, emissionArray)
!! !    type (t_GmiArrayBundle), pointer, intent(out) :: emissionArray (:)
!!     type (t_GmiArrayBundle), pointer :: emissionArray (:)
!!     type (t_Emission), intent(in)   :: self
!! !    emissionArray(:) = self%emissionArray(:)
!!     emissionArray => self%emissionArray
!!     return
!!   end subroutine Get_emissionArray
!-------------------------------------------------------------------------
  subroutine Allocate_base_isop (self, i1, i2, ju1, j2)
    integer          , intent(in   ) :: i1, i2, ju1, j2
    type (t_Emission), intent(inOut) :: self
    Allocate(self%base_isop(i1:i2, ju1:j2, NTYPE))
    self%base_isop = 0.0d0
    return
  end subroutine Allocate_base_isop
!-------------------------------------------------------------------------
  subroutine Allocate_base_monot (self, i1, i2, ju1, j2)
    integer          , intent(in   ) :: i1, i2, ju1, j2
    type (t_Emission), intent(inOut) :: self
    Allocate(self%base_monot(i1:i2, ju1:j2, NTYPE))
    self%base_monot = 0.0d0
    return
  end subroutine Allocate_base_monot
!-------------------------------------------------------------------------
  subroutine Allocate_xlai (self, i1, i2, ju1, j2)
    integer          , intent(in   ) :: i1, i2, ju1, j2
    type (t_Emission), intent(inOut) :: self
    Allocate(self%xlai(i1:i2, ju1:j2, NTYPE))
    self%xlai = 0.0d0
    return
  end subroutine Allocate_xlai
!-------------------------------------------------------------------------
  subroutine Allocate_xlai2 (self, i1, i2, ju1, j2)
    integer          , intent(in   ) :: i1, i2, ju1, j2
    type (t_Emission), intent(inOut) :: self
    Allocate(self%xlai2(i1:i2, ju1:j2, NTYPE))
    self%xlai2 = 0.0d0
    return
  end subroutine Allocate_xlai2
!-------------------------------------------------------------------------
  subroutine Allocate_iland (self, i1, i2, ju1, j2)
    integer          , intent(in   ) :: i1, i2, ju1, j2
    type (t_Emission), intent(inOut) :: self
    Allocate(self%iland(i1:i2, ju1:j2, NTYPE))
    self%iland = 0
    return
  end subroutine Allocate_iland
!-------------------------------------------------------------------------
  subroutine Allocate_iuse (self, i1, i2, ju1, j2)
    integer          , intent(in   ) :: i1, i2, ju1, j2
    type (t_Emission), intent(inOut) :: self
    Allocate(self%iuse(i1:i2, ju1:j2, NTYPE))
    self%iuse = 0
    return
  end subroutine Allocate_iuse
!-------------------------------------------------------------------------
  subroutine Allocate_ireg (self, i1, i2, ju1, j2)
    integer          , intent(in   ) :: i1, i2, ju1, j2
    type (t_Emission), intent(inOut) :: self
    Allocate(self%ireg(i1:i2, ju1:j2))
    self%ireg = 0
    return
  end subroutine Allocate_ireg
!-------------------------------------------------------------------------
  subroutine Allocate_emissAero (self, i1, i2, ju1, j2)
    integer          , intent(in   ) :: i1, i2, ju1, j2
    type (t_Emission)       , intent(inOut) :: self
    Allocate(self%emissAero(i1:i2, ju1:j2, self%naero))
    self%emissAero = 0.0d0
    return
  end subroutine Allocate_emissAero
!-------------------------------------------------------------------------
  subroutine Allocate_emissDust (self, i1, i2, ju1, j2)
    integer          , intent(in   ) :: i1, i2, ju1, j2
    type (t_Emission)       , intent(inOut) :: self
    Allocate(self%emissDust(i1:i2, ju1:j2, 1:self%ndust))
    self%emissDust = 0.0d0
    return
  end subroutine Allocate_emissDust
!-------------------------------------------------------------------------
  subroutine Allocate_lightning_NO (self, i1, i2, ju1, j2, k1, k2)
    integer          , intent(in   ) :: i1, i2, ju1, j2, k1, k2
    type (t_Emission)      , intent(inOut) :: self
    Allocate(self%lightning_NO(i1:i2, ju1:j2, k1:k2))
    self%lightning_NO = 0.0d0
    return
  end subroutine Allocate_lightning_NO
!-------------------------------------------------------------------------
  subroutine Allocate_emiss_isop (self, i1, i2, ju1, j2)
    integer          , intent(in   ) :: i1, i2, ju1, j2
    type (t_Emission)      , intent(inOut) :: self
    Allocate(self%emiss_isop(i1:i2, ju1:j2))
    self%emiss_isop = 0.0d0
    return
  end subroutine Allocate_emiss_isop
!-------------------------------------------------------------------------
  subroutine Allocate_emiss_monot (self, i1, i2, ju1, j2)
    integer          , intent(in   ) :: i1, i2, ju1, j2
    type (t_Emission)      , intent(inOut) :: self
    Allocate(self%emiss_monot(i1:i2, ju1:j2))
    self%emiss_monot = 0.0d0
    return
  end subroutine Allocate_emiss_monot
!-------------------------------------------------------------------------
  subroutine Allocate_emiss_nox (self, i1, i2, ju1, j2)
    integer          , intent(in   ) :: i1, i2, ju1, j2
    type (t_Emission)      , intent(inOut) :: self
    Allocate(self%emiss_nox(i1:i2, ju1:j2))
    self%emiss_nox = 0.0d0
    return
  end subroutine Allocate_emiss_nox
!-------------------------------------------------------------------------
  subroutine Allocate_emiss_hno3 (self, i1, i2, ju1, j2)
    integer          , intent(in   ) :: i1, i2, ju1, j2
    type (t_Emission)      , intent(inOut) :: self
    Allocate(self%emiss_hno3(i1:i2, ju1:j2))
    self%emiss_hno3 = 0.0d0
    return
  end subroutine Allocate_emiss_hno3
!-------------------------------------------------------------------------
  subroutine Allocate_emiss_o3 (self, i1, i2, ju1, j2)
    integer          , intent(in   ) :: i1, i2, ju1, j2
    type (t_Emission)      , intent(inOut) :: self
    Allocate(self%emiss_o3(i1:i2, ju1:j2))
    self%emiss_o3 = 0.0d0
    return
  end subroutine Allocate_emiss_o3
!-------------------------------------------------------------------------
  subroutine Allocate_surf_emiss_out (self, i1, i2, ju1, j2, numSpecies)
    integer          , intent(in   ) :: i1, i2, ju1, j2, numSpecies
    type (t_Emission)      , intent(inOut) :: self
    Allocate(self%surf_emiss_out(i1:i2, ju1:j2, numSpecies))
    self%surf_emiss_out = 0.0d0
    return
  end subroutine Allocate_surf_emiss_out
!-------------------------------------------------------------------------
  subroutine Allocate_surf_emiss_out2 (self, i1, i2, ju1, j2)
    integer          , intent(in   ) :: i1, i2, ju1, j2
    type (t_Emission)      , intent(inOut) :: self
    Allocate(self%surf_emiss_out2(i1:i2, ju1:j2, 6))
    self%surf_emiss_out2 = 0.0d0
    return
  end subroutine Allocate_surf_emiss_out2
!-------------------------------------------------------------------------
  subroutine Allocate_emiss_3d_out (self, i1, i2, ju1, j2, k1, k2, numSpecies)
    integer          , intent(in   ) :: i1, i2, ju1, j2, k1, k2, numSpecies
    type (t_Emission), intent(inOut) :: self
    Allocate(self%emiss_3d_out(i1:i2, ju1:j2, k1:k2, numSpecies))
    self%emiss_3d_out = 0.0d0
    return
  end subroutine Allocate_emiss_3d_out
!-------------------------------------------------------------------------
  subroutine Allocate_aerosolEmiss3D (self, i1, i2, ju1, j2, k1, k2)
    integer          , intent(in   ) :: i1, i2, ju1, j2, k1, k2
    type (t_Emission)      , intent(inOut) :: self
    Allocate(self%aerosolEmiss3D(i1:i2, ju1:j2, k1:k2, 5))
    self%aerosolEmiss3D = 0.0d0
    return
  end subroutine Allocate_aerosolEmiss3D
!-------------------------------------------------------------------------
  subroutine Allocate_aerosolSurfEmiss (self, i1, i2, ju1, j2)
    integer          , intent(in   ) :: i1, i2, ju1, j2
    type (t_Emission)      , intent(inOut) :: self
    Allocate(self%aerosolSurfEmiss(i1:i2, ju1:j2, &
        self%ndust + self%naero + 5))
    self%aerosolSurfEmiss = 0.0d0
    return
  end subroutine Allocate_aerosolSurfEmiss
!-------------------------------------------------------------------------
  subroutine Allocate_aerosolSurfEmissMap (self)
    type (t_Emission)      , intent(inOut) :: self
    Allocate(self%aerosolSurfEmissMap &
        (self%ndust + self%naero + 5))
    return
  end subroutine Allocate_aerosolSurfEmissMap
!-------------------------------------------------------------------------
  subroutine Allocate_emissDust_t (self, i1, i2, ju1, j2)
    integer          , intent(in   ) :: i1, i2, ju1, j2
    type (t_Emission), intent(inOut) :: self
    integer                          :: ndust, nt_dust, nst_dust
    ndust    = self%ndust
    nst_dust = self%nst_dust
    nt_dust  = self%nt_dust
    Allocate(self%emissDust_t(i1:i2, ju1:j2, 1:ndust, nst_dust:nst_dust+nt_dust-1))    
    self%emissDust_t = 0.0d0
    return
  end subroutine Allocate_emissDust_t
!-------------------------------------------------------------------------
  subroutine Allocate_emissAero_t (self, i1, i2, ju1, j2)
    integer          , intent(in   ) :: i1, i2, ju1, j2
    type (t_Emission), intent(inOut) :: self
    integer                          :: naero, emiss_timpyr
    naero        = self%naero
    emiss_timpyr = self%emiss_timpyr
    Allocate(self%emissAero_t(i1:i2, ju1:j2, 1:naero, emiss_timpyr))    
    self%emissAero_t = 0.0d0
    return
  end subroutine Allocate_emissAero_t
!-------------------------------------------------------------------------
  subroutine Get_do_gcr (self, do_gcr)
    logical        , intent(out)  :: do_gcr
    type (t_Emission), intent(in)   :: self
    do_gcr = self%do_gcr
    return
  end subroutine Get_do_gcr
!-------------------------------------------------------------------------
  subroutine Get_do_semiss_inchem (self, do_semiss_inchem)
    logical        , intent(out)  :: do_semiss_inchem
    type (t_Emission), intent(in)   :: self
    do_semiss_inchem = self%do_semiss_inchem
    return
  end subroutine Get_do_semiss_inchem
!-------------------------------------------------------------------------
  subroutine Allocate_isoLaiNext(self, i1, i2, ju1, j2)
    integer, intent(in) :: i1, i2, ju1, j2
    type (t_Emission), intent(inOut)   :: self
    allocate(self%isoLaiNext(i1:i2,ju1:j2))
    self%isoLaiNext = 0.0d0
    return
  end subroutine Allocate_isoLaiNext
!-------------------------------------------------------------------------
  subroutine Allocate_isoLaiCurr(self, i1, i2, ju1, j2)
    integer, intent(in) :: i1, i2, ju1, j2
    type (t_Emission), intent(inOut)   :: self
    allocate(self%isoLaiCurr(i1:i2,ju1:j2))
    self%isoLaiCurr = 0.0d0
    return
  end subroutine Allocate_isoLaiCurr
!-------------------------------------------------------------------------
  subroutine Allocate_isoLaiPrev(self, i1, i2, ju1, j2)
    integer, intent(in) :: i1, i2, ju1, j2
    type (t_Emission), intent(inOut)   :: self
    allocate(self%isoLaiPrev(i1:i2,ju1:j2))
    self%isoLaiPrev = 0.0d0
    return
  end subroutine Allocate_isoLaiPrev
!-------------------------------------------------------------------------
  subroutine Allocate_isoLaiYear(self, i1, i2, ju1, j2)
    integer, intent(in) :: i1, i2, ju1, j2
    type (t_Emission), intent(inOut)   :: self
    allocate(self%isoLaiYear(i1:i2,ju1:j2,12))
    self%isoLaiYear = 0.0d0
    return
  end subroutine Allocate_isoLaiYear
!-------------------------------------------------------------------------
  subroutine Allocate_aefMbo(self, i1, i2, ju1, j2)
    integer, intent(in) :: i1, i2, ju1, j2
    type (t_Emission), intent(inOut)   :: self
    allocate(self%aefMbo(i1:i2,ju1:j2))
    self%aefMbo = 0.0d0
    return
  end subroutine Allocate_aefMbo
!-------------------------------------------------------------------------
  subroutine Allocate_aefIsop(self, i1, i2, ju1, j2)
    integer, intent(in) :: i1, i2, ju1, j2
    type (t_Emission), intent(inOut)   :: self
    allocate(self%aefIsop(i1:i2,ju1:j2))
    self%aefIsop = 0.0d0
    return
  end subroutine Allocate_aefIsop
!-------------------------------------------------------------------------
  subroutine Allocate_aefOvoc(self, i1, i2, ju1, j2)
    integer, intent(in) :: i1, i2, ju1, j2
    type (t_Emission), intent(inOut)   :: self
    allocate(self%aefOvoc(i1:i2,ju1:j2))
    self%aefOvoc = 0.0d0
    return
  end subroutine Allocate_aefOvoc
!-------------------------------------------------------------------------
  subroutine Allocate_aefMonot(self, i1, i2, ju1, j2)
    integer, intent(in) :: i1, i2, ju1, j2
    type (t_Emission), intent(inOut)   :: self
    allocate(self%aefMonot(i1:i2,ju1:j2))
    self%aefMonot = 0.0d0
    return
  end subroutine Allocate_aefMonot
!-------------------------------------------------------------------------
  subroutine Get_lightning_opt (self, lightning_opt)
    integer        , intent(out)  :: lightning_opt
    type (t_Emission), intent(in)   :: self
    lightning_opt = self%lightning_opt
    return
  end subroutine Get_lightning_opt
!-------------------------------------------------------------------------
  subroutine Get_emiss_aero_opt (self, emiss_aero_opt)
    integer        , intent(out)  :: emiss_aero_opt
    type (t_Emission), intent(in)   :: self
    emiss_aero_opt = self%emiss_aero_opt
    return
  end subroutine Get_emiss_aero_opt
!-------------------------------------------------------------------------
  subroutine Get_emiss_dust_opt (self, emiss_dust_opt)
    integer        , intent(out)  :: emiss_dust_opt
    type (t_Emission), intent(in)   :: self
    emiss_dust_opt = self%emiss_dust_opt
    return
  end subroutine Get_emiss_dust_opt
!-------------------------------------------------------------------------
  subroutine Get_num_emiss (self, num_emiss)
    integer        , intent(out)  :: num_emiss
    type (t_Emission), intent(in)   :: self
    num_emiss = self%num_emiss
    return
  end subroutine Get_num_emiss
!-------------------------------------------------------------------------
  subroutine Get_ndust (self, ndust)
    integer        , intent(out)  :: ndust
    type (t_Emission), intent(in)   :: self
    ndust = self%ndust
    return
  end subroutine Get_ndust
!-------------------------------------------------------------------------
  subroutine Get_naero (self, naero)
    integer        , intent(out)  :: naero
    type (t_Emission), intent(in)   :: self
    naero = self%naero
    return
  end subroutine Get_naero
!-------------------------------------------------------------------------
  end module GmiEmissionMethod_mod
