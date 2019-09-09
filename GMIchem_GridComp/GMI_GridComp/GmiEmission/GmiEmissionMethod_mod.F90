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
      use MAPL_Mod
      use GmiTimeControl_mod  , only : t_GmiClock, Get_curGmiDate, GmiSplitDateTime
      use GmiGrid_mod, only : t_gmiGrid
      use GmiGrid_mod, only : Get_i1, Get_i2, Get_ju1, Get_j2, Get_k1, Get_k2
      use GmiGrid_mod, only : Get_i1_gl, Get_i2_gl, Get_ju1_gl, Get_j2_gl
      use GmiGrid_mod, only : Get_ilo, Get_ihi, Get_julo, Get_jhi
      use GmiGrid_mod, only : Get_ilong, Get_ilat, Get_ivert
      use GmiPrintError_mod, only : GmiPrintError
      use GmiCheckRange_mod, only : CheckRange3d
      use GmiESMFrcFileReading_mod, only : rcEsmfReadTable, rcEsmfReadLogical
      use GmiArrayBundlePointer_mod, only : t_GmiArrayBundle, CleanArrayPointer
      use GmiSpcConcentrationMethod_mod, only : t_SpeciesConcentration
      use GmiSpcConcentrationMethod_mod, only : Get_concentration, Set_concentration
      use GmiSpcConcentrationMethod_mod, only : Get_tracer_opt, Get_tr_source_land
      use GmiSpcConcentrationMethod_mod, only : Get_tr_source_ocean
!
      implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
!
      private
      public  :: InitializeEmission, initReadEmission
      public  :: RunEmission
      public  :: FinalizeEmission
    
      public  :: Set_emissionArray
      public  :: Set_emiss_isop  , Set_emiss_monot, Set_emiss_nox
      public  :: Set_emiss_o3    , Set_emiss_hno3
      public  :: Set_lightning_NO, Set_flashrate
      public  :: Set_emissDust   , Set_emissAero
      public  :: Set_emissDust_t , Set_emissAero_t
      public  :: Set_isop_scale  , Set_aerosolSurfEmiss, Set_aerosolEmiss3D
      public  :: Set_emiss_3d_out, Set_surf_emiss_out  , Set_surf_emiss_out2

      public  :: Get_emissionArray
      public  :: Get_emiss_isop          , Get_emiss_monot , Get_emiss_nox
      public  :: Get_emiss_o3            , Get_emiss_hno3
      public  :: Get_emiss_opt           , Get_emiss_in_opt
      public  :: Get_emiss_map           , Get_emiss_timpyr
      public  :: Get_lightning_NO        , Get_flashrate
      public  :: Get_emissDust           , Get_emissAero
      public  :: Get_emissDust_t         , Get_emissAero_t
      public  :: Get_fertscal_infile_name
      public  :: Get_lai_infile_name     , Get_precip_infile_name
      public  :: Get_soil_infile_name    , Get_veg_infile_name
      public  :: Get_isopconv_infile_name, Get_monotconv_infile_name
      public  :: Get_lightning_opt       , Get_emiss_aero_opt
      public  :: Get_emiss_dust_opt      , Get_num_emiss
      public  :: Get_ndust               , Get_naero
      public  :: Get_nst_dust            , Get_nt_dust 
      public  :: Get_emiss_map_dust      , Get_emiss_map_aero
      public  :: Get_isop_scale          , Get_do_gcr
      public  :: Get_do_semiss_inchem    , Get_do_ShipEmission, Get_doMEGANemission
      public  :: Get_o3_index            , Get_hno3_index
      public  :: Get_doReadDailyEmiss
      public  :: Get_lightNOampFactor    , Get_numberNOPerFlash, Get_minDeepCloudTop
      public  :: Get_begDailyEmissRec    , Get_endDailyEmissRec
      public  :: Get_aerosolSurfEmiss    , Get_aerosolSurfEmissMap, Get_aerosolEmiss3D
      public  :: Get_emiss_3d_out        , Get_surf_emiss_out     , Get_surf_emiss_out2
      public  :: Get_ireg        , Get_iuse
      public  :: Get_iland       , Get_ncon_soil
      public  :: Get_xlai        , Get_xlai2
      public  :: Get_base_isop   , Get_base_monot
      public  :: Get_coeff_isop  , Get_convert_isop, Get_convert_monot
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
    character (len=MAX_LENGTH_FILE_NAME) :: fertscal_infile_name          ! fertilizer scale     input file name
    character (len=MAX_LENGTH_FILE_NAME) :: lai_infile_name               ! leaf area index      input file name
    character (len=MAX_LENGTH_FILE_NAME) :: precip_infile_name            ! precipitation        input file name
    character (len=MAX_LENGTH_FILE_NAME) :: soil_infile_name              ! soil type            input file name
    character (len=MAX_LENGTH_FILE_NAME) :: veg_infile_name               ! vegetation type      input file name
    character (len=MAX_LENGTH_FILE_NAME) :: isopconv_infile_name          ! isoprene convert     input file name
    character (len=MAX_LENGTH_FILE_NAME) :: monotconv_infile_name         ! monoterpene convert  input file name
    character (len=MAX_LENGTH_FILE_NAME) :: GOCARTerod_infile_name
    character (len=MAX_LENGTH_FILE_NAME) :: GOCARTocean_infile_name
    character (len=MAX_LENGTH_FILE_NAME) :: GOCARTerod_mod_infile_name
    real*8              :: isop_scale (12)               ! array of monthly isoprene scaling coefficients

    CHARACTER(LEN=64):: emissionSpeciesNames(MAX_NUM_CONST)  ! Names of emssions on netCDF/hdf file #

    integer          :: lightning_opt		         ! lightning option
    real             :: lightNOampFactor	         ! Lightning NO production amplification/suppression factor 
    real             :: numberNOPerFlash	         ! NO molecules generated by each flash
    real             :: minDeepCloudTop	                 ! Minimum cloud top [km] for selecting deep convection profiles
    real*8 , pointer :: flashrate      (:,:) => null()   ! Imported 2D flash rates
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
    character (len=MAX_LENGTH_FILE_NAME) :: laiMEGAN_InfileName    ! Inpuf file name for AVHRR
                                                  ! leaf-area-indices
    character (len=MAX_LENGTH_FILE_NAME) :: aefMboMEGAN_InfileName ! Annual emission factor for
                                                  ! methyl butenol input file name
    character (len=MAX_LENGTH_FILE_NAME) :: aefIsopMEGAN_InfileName ! Annual emission factor for
                                                   ! isoprene input file name
    character (len=MAX_LENGTH_FILE_NAME) :: aefMonotMEGAN_InfileName ! Annual emission factor for
                                                    ! monoterpenes input file name
    character (len=MAX_LENGTH_FILE_NAME) :: aefOvocMEGAN_InfileName ! Annual emission factor for other
                                                   ! biogenic VOCs input file name
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
    integer          :: o3_index   ! index of O3   in emiss_map
    integer          :: hno3_index ! index of HNO3 in emiss_map
    real*8 , pointer :: emiss          (:,:,:,:)   => null() ! array of emissions (kg/s)
    real*8 , pointer :: emissDust     (:,:,:)     => null() ! used in sulfur chemistry
    real*8 , pointer :: emiss_isop     (:,:)       => null() ! isoprene    emissions (kg/s)
    real*8 , pointer :: emiss_monot    (:,:)       => null() ! monoterpene emissions (kg/s)
    real*8 , pointer :: emiss_nox      (:,:)       => null() ! NOx         emissions (kg/s)
    real*8 , pointer :: emiss_o3    (:,:)       => null() ! ozone       emissions
    real*8 , pointer :: emiss_ozone (:,:)       => null() ! ozone ship emissions
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

    logical       :: doReadDailyEmiss
    integer       :: begDailyEmissRec, endDailyEmissRec

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
      subroutine readEmissionResourceFile(self, config, procID, numSpecies, &
                     pr_diag, chem_opt)
!
! !USES:
      use GmiCheckNamelistFile_mod, only : CheckNamelistOptionRange
      use GmiSpeciesRegistry_mod, only : getSpeciesIndex, UNKNOWN_SPECIES
      use GmiStringManipulation_mod, only : constructListNames
!
! !INPUT PARAMETERS:
      logical, intent(in) :: pr_diag
      integer, intent(in) :: numSpecies, procID, chem_opt
!
! !INPUT/OUTPUT PARAMETERS:
      type (t_Emission), intent(inOut) :: self
      type (ESMF_Config), intent(inOut) :: config
!
! !DESCRIPTION:
! Reads in Emission related variables from the resource file.
!
! !LOCAL VARIABLES:
      integer ::  STATUS, RC, ic, ios
      character (len=MAX_LENGTH_SPECIES_NAME) :: tempListNames(MAX_NUM_CONST)
      character (len=MAX_STRING_LENGTH      ) :: emissionSpeciesNames
      character (len=MAX_STRING_LENGTH      ) :: emissionDustSpeciesNames
      character (len=MAX_STRING_LENGTH      ) :: emissionAeroSpeciesNames
      character(len=ESMF_MAXSTR) :: IAm, err_msg

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
     &                label   = "emiss_opt:", &
     &                default = 0, rc=STATUS )
      _VERIFY(STATUS)
      
    ! --------------------------------------------
    ! emiss_in_opt
    !   0:  no emissions data
    !   1:  set all emiss values to emiss_init_val
    !   2:  read in emiss values
    !   3:  set in emiss values using constant source terms
    ! --------------------------------------------

      call ESMF_ConfigGetAttribute(config, self%emiss_in_opt, &
     &                label   = "emiss_in_opt:", &
     &                default = 0, rc=STATUS )
      _VERIFY(STATUS)

    ! ------------------------------------
    ! emiss_conv_flag
    !   0:  no conversion performed
    !   1:  use emiss_conv_fac
    !   2:  convert from kg/km2-hr to kg/s
    ! ------------------------------------

      call ESMF_ConfigGetAttribute(config, self%emiss_conv_flag, &
     &                label   = "emiss_conv_flag:", &
     &                default = 0, rc=STATUS )
      _VERIFY(STATUS)

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
     &                label   = "semiss_inchem_flag:", &
     &                default = -1, rc=STATUS )
      _VERIFY(STATUS)
      
      ! sets of emissons per year (1 => yearly, 12 => monthly)
    
      call ESMF_ConfigGetAttribute(config, self%emiss_timpyr, &
     &                label   = "emiss_timpyr:", &
     &                default = 1, rc=STATUS )
      _VERIFY(STATUS)
      
      self%emiss_map(:) =  0  
     
      call rcEsmfReadTable(config, emissionSpeciesNames, &
     &                     "emissionSpeciesNames::", rc=STATUS)
    
      call ESMF_ConfigGetAttribute(config, self%emiss_conv_fac, &
     &                label   = "emiss_conv_fac:", &
     &                default = 1.0d0, rc=STATUS )
      _VERIFY(STATUS)
    
      call ESMF_ConfigGetAttribute(config, self%emiss_init_val, &
     &                label   = "emiss_init_val:", &
     &                default = 1.0d0, rc=STATUS )
      _VERIFY(STATUS)

      CALL rcEsmfReadLogical(config, self%clim_emiss_by_area, &
     &           "clim_emiss_by_area:", DEFAULT=.TRUE., RC=STATUS)
      _VERIFY(STATUS)

     ! Save the number of emitting layers for each emissionSpeciesName
     ! ---------------------------------------------------------------

      self%emissionSpeciesLayers(:) = 0

      call rcEsmfReadTable(config, self%emissionSpeciesLayers, &
     &                     "emissionSpeciesLayers::", rc=STATUS)

    !---------------------------------------
    ! Reading of daily emission file options
    !---------------------------------------

      call rcEsmfReadLogical(config, self%doReadDailyEmiss, &
     &           "doReadDailyEmiss:", default=.false., rc=STATUS)
      _VERIFY(STATUS)

      call ESMF_ConfigGetAttribute(config, self%begDailyEmissRec, &
     &                label   = "begDailyEmissRec:", &
     &                default = 1, rc=STATUS )
      _VERIFY(STATUS)

      call ESMF_ConfigGetAttribute(config, self%endDailyEmissRec, &
     &                label   = "endDailyEmissRec:", &
     &                default = 366, rc=STATUS )
      _VERIFY(STATUS)
      
    ! --------------------------------------
    ! Aerosols and Sulfur from Penner et al.
    ! --------------------------------------

     ! emiss_aero_opt   0: for no     aerosol emissions
     !                  1: for GMI    aerosol emissions
     !                  2: for GOCART aerosol emissions
      
      call ESMF_ConfigGetAttribute(config, self%emiss_aero_opt, &
     &                label   = "emiss_aero_opt:", &
     &                default = 0, rc=STATUS )
      _VERIFY(STATUS)

      ! number of aerosol emissions
     
      call ESMF_ConfigGetAttribute(config, self%naero, &
     &                label   = "naero:", &
     &                default = 1, rc=STATUS )
      _VERIFY(STATUS) 
     
      self%emiss_map_aero(:) = 0    ! map emissions to species #

      call rcEsmfReadTable(config, emissionAeroSpeciesNames, &
     &                     "emissionAeroSpeciesNames::", rc=STATUS)

      ! aerosol emission input file

      call ESMF_ConfigGetAttribute(config, self%emiss_aero_infile_name, &
     &                label   = "emiss_aero_infile_name:", &
     &                default = ' ', rc=STATUS )
      _VERIFY(STATUS)

     ! emiss_dust_opt   0: for no     dust emissions
     !                  1: for GMI    dust emissions
     !                  2: for GOCART dust emissions

      call ESMF_ConfigGetAttribute(config, self%emiss_dust_opt, &
     &                label   = "emiss_dust_opt:", &
     &                default = 0, rc=STATUS )
      _VERIFY(STATUS)

      ! number of dust bins

      call ESMF_ConfigGetAttribute(config, self%ndust, &
     &                label   = "ndust:", &
     &                default = 1, rc=STATUS )
      _VERIFY(STATUS)

     ! number of starting point for new dust emiss.
      
      call ESMF_ConfigGetAttribute(config, self%nst_dust, &
     &                label   = "nst_dust:", &
     &                default = 1, rc=STATUS )
      _VERIFY(STATUS)
      
      ! number of dust emissions per emiss. file
     
      call ESMF_ConfigGetAttribute(config, self%nt_dust, &
     &                label   = "nt_dust:", &
     &                default = 1, rc=STATUS )
      _VERIFY(STATUS)   
     
      self%emiss_map_dust(:) = 0    ! map emissions to species #
      
      emissionDustSpeciesNames = ""
     
      call rcEsmfReadTable(config, emissionDustSpeciesNames, &
     &                     "emissionDustSpeciesNames::", rc=STATUS)
      
      ! dust emission input file
      
      call ESMF_ConfigGetAttribute(config, self%emiss_dust_infile_name, &
     &                label   = "emiss_dust_infile_name:", &
     &                default = ' ', rc=STATUS )
      _VERIFY(STATUS)
     
     ! -----------------------------------------------------
     ! Information on scale factors for various NO emissions
     ! -----------------------------------------------------

     ! Fossil fuel
      call rcEsmfReadLogical(config, self%doScaleNOffEmiss, &
     &           "doScaleNOffEmiss:", default=.false., rc=STATUS)

      call ESMF_ConfigGetAttribute(config, self%scFactorNOff_infile_name, &
     &                label   = "scFactorNOff_infile_name:", &
     &                default = ' ', rc=STATUS )
      _VERIFY(STATUS)

     ! Biomass burning
      call rcEsmfReadLogical(config, self%doScaleNObbEmiss, &
     &           "doScaleNObbEmiss:", default=.false., rc=STATUS)
      _VERIFY(STATUS)

      call ESMF_ConfigGetAttribute(config, self%scFactorNObb_infile_name, &
     &                label   = "scFactorNObb_infile_name:", &
     &                default = ' ', rc=STATUS )
      _VERIFY(STATUS)

    ! ---------------------------------------------------
    ! Harvard emissions:  acetone, isoprene, propene, NO.
    ! ---------------------------------------------------

      self%isop_scale(:) = 1.0d0
     
      call rcEsmfReadTable(config, self%isop_scale, &
     &                     "isop_scale::", rc=STATUS)
     
    !     --------------------------------------
    !     Aerosols and Sulfur from GOCART
    !     --------------------------------------
      
      call ESMF_ConfigGetAttribute(config, self%GOCARTerod_infile_name, &
     &                label   = "GOCARTerod_infile_name:", &
     &                default = ' ', rc=STATUS )
      _VERIFY(STATUS)
     
      call ESMF_ConfigGetAttribute(config, self%GOCARTocean_infile_name, &
     &                label   = "GOCARTocean_infile_name:", &
     &                default = ' ', rc=STATUS )
      _VERIFY(STATUS)
      
      call ESMF_ConfigGetAttribute(config, self%GOCARTerod_mod_infile_name, &
     &                label   = "GOCARTerod_mod_infile_name:", &
     &                default = ' ', rc=STATUS )
      _VERIFY(STATUS)
    
    !... do Galactic Cosmic Rays source of NOx?
    
      ! turn on Galactic Cosmic ray emmission of N and NO

      call rcEsmfReadLogical(config, self%do_gcr, &
     &           "do_gcr:", default=.false., rc=STATUS)

      call ESMF_ConfigGetAttribute(config, self%gcr_infile_name, &
     &                label   = "gcr_infile_name:", &
     &                default = ' ', rc=STATUS )
      _VERIFY(STATUS)

    !... do Ship Emission Calculations

      call rcEsmfReadLogical(config, self%do_ShipEmission, &
     &           "do_ShipEmission:", default=.false., rc=STATUS)

    !... do MEGAN emission calculations

      call rcEsmfReadLogical(config, self%doMEGANemission, &
     &           "doMEGANemission:", default=.false., rc=STATUS)
      _VERIFY(STATUS)

      call ESMF_ConfigGetAttribute(config, self%MEGAN_infile_name, &
     &                label   = "MEGAN_infile_name:", &
     &                default = ' ', rc=STATUS )
      _VERIFY(STATUS)

      call ESMF_ConfigGetAttribute(config, self%laiMEGAN_InfileName, &
     &                label   = "laiMEGAN_InfileName:", &
     &                default = ' ', rc=STATUS )
      _VERIFY(STATUS)

      call ESMF_ConfigGetAttribute(config, self%aefMboMEGAN_InfileName, &
     &                label   = "aefMboMEGAN_InfileName:", &
     &                default = ' ', rc=STATUS )
      _VERIFY(STATUS)
     
      call ESMF_ConfigGetAttribute(config, self%aefIsopMEGAN_InfileName, &
     &                label   = "aefIsopMEGAN_InfileName:", &
     &                default = ' ', rc=STATUS )
      _VERIFY(STATUS) 
      
      call ESMF_ConfigGetAttribute(config, self%aefOvocMEGAN_InfileName, &
     &                label   = "aefOvocMEGAN_InfileName:", &
     &                default = ' ', rc=STATUS )
      _VERIFY(STATUS)
     
      call ESMF_ConfigGetAttribute(config, self%aefMonotMEGAN_InfileName, &
     &                label   = "aefMonotMEGAN_InfileName:", &
     &                default = ' ', rc=STATUS )
      _VERIFY(STATUS)
     
      call ESMF_ConfigGetAttribute(config, self%fertscal_infile_name, &
     &                label   = "fertscal_infile_name:", &
     &                default = 'fertscale_4x5_dao.asc', rc=STATUS )
      _VERIFY(STATUS)

      call ESMF_ConfigGetAttribute(config, self%lai_infile_name, &
     &                label   = "lai_infile_name:", &
     &                default = 'lai_4x5_dao.asc', rc=STATUS )
      _VERIFY(STATUS)

      call ESMF_ConfigGetAttribute(config, self%precip_infile_name, &
     &                label   = "precip_infile_name:", &
     &                default = 'precip_4x5_dao.asc', rc=STATUS )
      _VERIFY(STATUS)

      call ESMF_ConfigGetAttribute(config, self%soil_infile_name, &
     &                label   = "soil_infile_name:", &
     &                default = 'soiltype.asc', rc=STATUS )
      _VERIFY(STATUS)

      call ESMF_ConfigGetAttribute(config, self%veg_infile_name, &
     &                label   = "veg_infile_name:", &
     &                default = 'vegtype_4x5_dao.asc', rc=STATUS )
      _VERIFY(STATUS)

      call ESMF_ConfigGetAttribute(config, self%isopconv_infile_name, &
     &                label   = "isopconv_infile_name:", &
     &                default = 'isopconvtable.asc', rc=STATUS )
      _VERIFY(STATUS)

      call ESMF_ConfigGetAttribute(config, self%monotconv_infile_name, &
     &                label   = "monotconv_infile_name:", &
     &                default = 'monotconvtable.asc', rc=STATUS )
      _VERIFY(STATUS)

    ! ---------------------------------------------------------------------------------------
    ! lightning_opt = 0 --> default lightning
    !               = 1 --> variable_lightning (Dale Allen's algorithm)
    !               = 2 --> no lightning
    !               = in options 1 and 2, NO_lgt (in gmi_emiss.F) is made zero
    !
    ! lightNOampFactor = 1.0 --> NO production amplification/suppression factor, > 0
    ! numberNOperFlash = 1.50E+26 --> NO molecules generated by each flash
    ! minDeepCloudTop = 7.0 --> Minimum cloud top [km] for selecting deep convection profiles
    ! ---------------------------------------------------------------------------------------

      call ESMF_ConfigGetAttribute(config, self%lightning_opt, &
     &                label   = "lightning_opt:", &
     &                default = 0, rc=STATUS )
      _VERIFY(STATUS)
     
      call ESMF_ConfigGetAttribute(config, self%lightNOampFactor, &
     &                label   = "lightNOampFactor:", &
     &                default = 1.0, rc=STATUS )
      _VERIFY(STATUS)
     
      call ESMF_ConfigGetAttribute(config, self%numberNOperFlash, &
     &                label   = "numberNOperFlash:", &
     &                default = 1.50E+26, rc=STATUS )
      _VERIFY(STATUS)
    
      call ESMF_ConfigGetAttribute(config, self%minDeepCloudTop, &
     &                label   = "minDeepCloudTop:", &
     &                default = 7.0, rc=STATUS )
      _VERIFY(STATUS)
    
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
      else

         ! Set the initial value of the list
         tempListNames(:) = ''

         ! Construct the list of names using the long string
         call constructListNames(tempListNames, emissionSpeciesNames)

         self%num_emiss = Count (tempListNames(:) /= '')

         ! The specie name precedes the underscore. In GEOS-5 we need to keep
         ! the entire emissionSpeciesName for use in the initialize method.
         ! -------------------------------------------------------------------
         IF (self%num_emiss > 0) THEN
            DO ic = 1, self%num_emiss
               err_msg = tempListNames(ic)
               self%emissionSpeciesNames(ic) = TRIM(err_msg)
               ios = INDEX(err_msg,'_')
               IF(ios > 0) err_msg = err_msg(1:ios-1)
               self%emiss_map(ic) = getSpeciesIndex(TRIM(err_msg))
            END DO
         END IF

         self%num_emiss = count( self%emiss_map(:) > 0 )

         IF(MAPL_AM_I_ROOT()) THEN
          PRINT *," "
          PRINT *,"Num_emiss=",self%num_emiss
          PRINT *," "
          PRINT *,"Emiss_map:"
          PRINT *, self%emiss_map
          PRINT *," "
          PRINT *,"Emission species layers:"
          PRINT *, self%emissionSpeciesLayers
          PRINT *," "
         END IF

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
     &    (self%num_emiss   > MAX_NUM_CONST)) then
         err_msg = 'MAX_NUM_CONST problem in rc File.'
         call GmiPrintError (err_msg, .true., 2, numSpecies, &
     &           self%num_emiss, 0, 0.0d0, 0.0d0)
      end if

      if ((self%emiss_opt /= 0) .and. (self%num_emiss == 0)) then
        err_msg = 'emiss_opt/num_emiss problem.'
        call GmiPrintError  &
     &    (err_msg, .true., 2, self%emiss_opt, self%num_emiss, 0, 0.0d0, 0.0d0)
      end if

      if ((self%emiss_in_opt == 0) .and. (self%emiss_opt /= 0)) then
        err_msg = 'emiss_in_opt/emiss_opt problem.'
        call GmiPrintError  &
     &    (err_msg, .true., 2, self%emiss_in_opt, self%emiss_opt, 0, 0.0d0, 0.0d0)
      end if
         
      if ((self%emiss_timpyr /= 1) .and.  &
     &    (self%emiss_timpyr /= MONTHS_PER_YEAR)) then
        err_msg = 'emiss_timpyr range problem.'
        call GmiPrintError  &
     &    (err_msg, .true., 1, self%emiss_timpyr, 0, 0, 0.0d0, 0.0d0)
      end if
      
      if ((self%emiss_in_opt /= 2) .and.  &
     &    (self%emiss_timpyr == MONTHS_PER_YEAR)) then
        err_msg = 'emiss_in_opt/emiss_timpyr problem.'
        call GmiPrintError (err_msg, .true., 2, self%emiss_in_opt,             &
     &                      self%emiss_timpyr, 0, 0.0d0, 0.0d0)
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
                     pr_diag)
!
! !USES:
      use ReadOtherEmissionData_mod, only : readFertilizerData, readPrecipitationData

! !INPUT PARAMETERS:
      logical,           intent(in) :: pr_diag
      integer,           intent(in) :: loc_proc
      real*8,            intent(in) :: mcor(:,:)
      type (t_gmiGrid),  intent(in) :: gmiGrid
      type (t_GmiClock), intent(in) :: gmiClock
!
! !INPUT/OUTPUT PARAMETERS:
      type (t_Emission)  , intent(inOut) :: self
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

  if (self%emiss_dust_opt == 1 ) then  ! GMI sulfurous dust emissions
     call InitEmissDust (self, mcor, loc_proc, rootProc, &
                         i1, i2, ju1, j2, i1_gl, ju1_gl, ilong, ilat)
  end if

  if (self%emiss_aero_opt /= 0 ) then
     call InitEmissAero (self, mcor, loc_proc, rootProc, &
                        i1, i2, ju1, j2, i1_gl, ju1_gl, ilong, ilat)
  end if

  call Get_curGmiDate(gmiClock, nymd)

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

  call readEmissionResourceFile(self, config, loc_proc, numSpecies, pr_diag, chem_opt)
  self%idaySoilType = -999
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
     CALL Allocate_aefMbo  (self, i1, i2, ju1, j2)
     CALL Allocate_aefIsop (self, i1, i2, ju1, j2)
     CALL Allocate_aefOvoc (self, i1, i2, ju1, j2)
     CALL Allocate_aefMonot(self, i1, i2, ju1, j2)
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
     Allocate(self%emiss_ozone(i1:i2, ju1:j2))
     
     id                       = 0
     self%o3_index   = 999
     self%hno3_index = 999
     
     if (io3_num   == 0) stop "the index of   O3 should be non zero"
     if (ihno3_num == 0) stop "the index of HNO3 should be non zero"

     do ix = 1, numSpecies
        ic = self%emiss_map(ix)
        if (ic > 0) id = id + 1
        if (ic == io3_num  ) self%o3_index   = id
        if (ic == ihno3_num) self%hno3_index = id
     end do

     if (self%o3_index   == 999) stop "O3   is not defined in emiss_map"
     if (self%hno3_index == 999) stop "HNO3 is not defined in emiss_map"
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
     call Allocate_flashrate    (self, i1, i2, ju1, j2)
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
     &              loc_proc, cosSolarZenithAngle, latdeg, mcor, mass,         &
     &              lwi_flags, radswg, TwoMeter_air_temp, surf_rough,          &
     &              con_precip, tot_precip, ustar, fracCloudCover, kel, pbl,   &
     &              cmf, press3c, press3e, dtrn, pctm2, u10m, v10m, &
     &              gwet, gridBoxHeight, mw, IBOC, IBBC, INOC, IFOC, IFBC,     &
     &              ISSLT1, ISSLT2, ISSLT3, ISSLT4, IFSO2, INSO2, INDMS, IAN,  &
     &              IMGAS, INO, iisoprene_num, ino_num, ico_num, ipropene_num, &
     &              ihno3_num, io3_num,  numSpecies, pardif, pardir, T_15_AVG, &
     &              met_opt, chem_opt, trans_opt, do_aerocom, do_drydep,       &
     &              pr_diag, pr_const, pr_surf_emiss, pr_emiss_3d,             &
     &              metdata_name_org, metdata_name_model, tdt4, mixPBL)
!
! !USES:
      use GmiTimeControl_mod, only : GmiSplitDateTime, GetDaysFromJanuary1
      use GmiTimeControl_mod, only : t_GmiClock
      use GmiTimeControl_mod, only : Get_curGmiDate  , Get_numTimeSteps
      use GmiEmissionLightning_mod, only  : emiss_lightning
      use GmiArrayBundlePointer_mod, only : t_GmiArrayBundle
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
      integer, intent(in) :: lwi_flags(:, :)
      real*8 , intent(in) :: radswg(:, :), TwoMeter_air_temp(:, :), surf_rough(:, :)
      real*8 , intent(in) :: pardif(:, :), pardir(:, :), T_15_AVG(:, :)
      real*8 , intent(in) :: con_precip(:, :), tot_precip(:, :), ustar(:, :)
      real*8 , intent(in) :: kel(:, :, :), pbl(:, :), u10m(:,:), v10m(:,:), gwet(:,:)
      real*8 , intent(in) :: fracCloudCover (:,:) , dtrn(:,:,:), cmf(:,:,:)
      real*8 , intent(in) :: gridBoxHeight (:,:,:) 
      character (len=MAX_LENGTH_MET_NAME), intent(in) :: metdata_name_org
      character (len=MAX_LENGTH_MET_NAME), intent(in) :: metdata_name_model
      real   , intent(in) :: tdt4
      logical, intent(in) :: mixPBL   ! whether to explicitly distribute
                                      ! aerosol emissions within the PBL
!
! !INPUT/OUTPUT VARIABLES:
      type (t_Emission),            intent(inOut) :: self
      type(t_SpeciesConcentration), intent(inOut) :: SpeciesConcentration
!
! !DESCRIPTION:
! Run method for the Emission component.
!
! !LOCAL VARIABLES:
      real*8 , allocatable :: productionNO(:,:,:)
      type(t_GmiArrayBundle), pointer :: concentration(:)

      character (len=128) err_msg
      integer       :: ydummy, thisDay, thisMonth, thisDate, ddummy, ic, curRecord
      logical       :: doReadDailyEmiss, rootProc
      integer       :: begDailyEmissRec, endDailyEmissRec
      integer       :: hno3_index
      integer       :: nymd, num_time_steps, ndt
      real*8        :: tdt8
      integer       :: i1, i2, ju1, j2, k, k1, k2, ilo, ihi, julo, jhi
      integer       :: i1_gl, i2_gl, ju1_gl, j2_gl, ilong, ilat, ivert
      real*8        :: tr_source_ocean, tr_source_land
      integer       :: rc, tracer_opt
!
!EOP
!------------------------------------------------------------------------------
!BOC
      ! Get the GMI grid information
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

      ! Obtain model time information

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
        self%emiss_isop  = 0.0d0
        self%emiss_nox   = 0.0d0
     end if

     if (self%do_ShipEmission) then
        self%emiss_hno3(:,:) = self%emissionArray(hno3_index)%pArray3D(:,:,1)
     end if

! For GOCART emission

     if ((self%emiss_dust_opt == 2) .or. (self%emiss_aero_opt == 2))  then
        call SetGocartDerivedVars (u10m, v10m, gwet, press3c, &
     &                    kel, mass, mcor, gridBoxHeight, lwi_flags,  &
     &                    i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi)
     end if

! Emission control routines

      rootProc = MAPL_AM_I_ROOT()

      call Update_Emiss (gmiGrid, self%idaySoilType, self%firstBiogenicBase,   &
                 lwi_flags, cosSolarZenithAngle, latdeg,      &
     &           mcor, self%emiss_isop, self%emiss_monot, self%emiss_nox,      &
     &           self%do_ShipEmission, self%emiss_hno3, self%emiss_o3,         &
     &           ihno3_num, io3_num, radswg, TwoMeter_air_temp, surf_rough,    &
     &           con_precip, tot_precip, ustar, mass, fracCloudCover, kel,     &
     &           self%surf_emiss_out, self%surf_emiss_out2, self%emiss_3d_out, &
     &           self%aerosolEmiss3D, self%aerosolSurfEmiss,                   &
     &           self%aerosolSurfEmissMap, concentration, self%emissionArray,  &
     &           self%emissDust_t, self%emissDust, self%emissAero_t,           &
     &           self%emissAero, pbl, gridBoxHeight, self%index_soil,          &
     &           self%ncon_soil, self%soil_fert, self%soil_precip,             &
     &           self%soil_pulse, self%ireg, self%iland, self%iuse,            &
     &           self%convert_isop, self%convert_monot, self%coeff_isop,       &
     &           self%base_isop, self%base_monot, self%xlai, IBOC, IBBC, INOC, &
     &           IFOC, IFBC, ISSLT1, ISSLT2, ISSLT3, ISSLT4, IFSO2, INSO2,     &
     &           INDMS, IAN, IMGAS, INO, iisoprene_num, ino_num, ico_num,      &
     &           ipropene_num, pr_surf_emiss, pr_emiss_3d, pr_diag, loc_proc,  &
     &           rootProc, met_opt, self%emiss_opt, chem_opt, trans_opt,       &
     &           self%emiss_aero_opt, self%emiss_dust_opt, do_aerocom,         &
     &           self%do_semiss_inchem, self%do_gcr, do_drydep, self%emiss_map,&
     &           self%emiss_map_dust, self%emiss_map_aero, self%ndust,         &
     &           self%nst_dust, self%nt_dust, self%naero, nymd, num_time_steps,&
     &           mw, tdt8, ndt, self%emiss_timpyr, self%num_emiss,             &
     &           self%isop_scale, i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi,&
     &           i1_gl, i2_gl, ju1_gl, j2_gl, ilong, numSpecies,               &
     &           self%doMEGANemission, self%aefIsop, self%aefMbo,              &
     &           self%aefMonot, self%isoLaiPrev, self%isoLaiCurr,              &
     &           self%isoLaiNext, pardif, pardir, T_15_AVG,                    &
     &           self%emissionSpeciesLayers, self%exp_fac, mixPBL)

! NO production from lightning, with flashRate imported from MOIST
! ----------------------------------------------------------------
      IF(self%lightning_opt == 1) THEN

     	ALLOCATE(productionNO(i1:i2, ju1:j2, k1:k2))
	productionNO(:,:,:) = 0.00

     	CALL emiss_lightning(i1, i2, ju1, j2, k1, k2, self%minDeepCloudTop, self%lightNOampFactor, &
        		     self%numberNOperFlash, lwi_flags, self%flashRate, gridBoxHeight, dtrn, &
     			     productionNO, self%lightning_NO, rc)

! Convert production rate, currently m^{-3} s^{-1}, to volume
! mixing ratio (mole fraction) tendency. Note that mass has units of kg.
! ----------------------------------------------------------------------
        DO k = k1,k2
	 productionNO(i1:i2,ju1:j2,k) = productionNO(i1:i2,ju1:j2,k)*MWTAIR*mcor(i1:i2,ju1:j2)* &
	                                gridBoxHeight(i1:i2,ju1:j2,k)/(1000.00*AVOGAD*mass(i1:i2,ju1:j2,k))
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
  integer :: lwi_flags(i1:i2, ju1:j2)

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
     &                self%emissionArray(ic)%pArray3D(il,ij,:) * mcor(il,ij) * units_fac
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
  subroutine Get_emissionArray (self, emissionArray)
!    type (t_GmiArrayBundle), pointer, intent(out) :: emissionArray (:)
    type (t_GmiArrayBundle), pointer :: emissionArray (:)
    type (t_Emission), intent(in)   :: self
!    emissionArray(:) = self%emissionArray(:)
    emissionArray => self%emissionArray
    return
  end subroutine Get_emissionArray
!-------------------------------------------------------------------------
  subroutine Set_emissionArray (self, emissionArray)
    type (t_GmiArrayBundle), pointer :: emissionArray (:)
    type (t_Emission), intent(inOut) :: self
    self%emissionArray => emissionArray
    return
  end subroutine Set_emissionArray
!-------------------------------------------------------------------------
  subroutine Allocate_base_isop (self, i1, i2, ju1, j2)
    integer          , intent(in   ) :: i1, i2, ju1, j2
    type (t_Emission), intent(inOut) :: self
    Allocate(self%base_isop(i1:i2, ju1:j2, NTYPE))
    self%base_isop = 0.0d0
    return
  end subroutine Allocate_base_isop
!-------------------------------------------------------------------------
  subroutine Get_base_isop (self, base_isop)
    real*8          , intent(out)  :: base_isop (:,:,:)
    type (t_Emission), intent(in)   :: self
    base_isop(:,:,:) = self%base_isop(:,:,:)
    return
  end subroutine Get_base_isop
!-------------------------------------------------------------------------
  subroutine Set_base_isop (self, base_isop)
    real*8          , intent(in)  :: base_isop (:,:,:)
    type (t_Emission), intent(inOut) :: self
    self%base_isop(:,:,:) = base_isop(:,:,:)
    return
  end subroutine Set_base_isop
!-------------------------------------------------------------------------
  subroutine Allocate_base_monot (self, i1, i2, ju1, j2)
    integer          , intent(in   ) :: i1, i2, ju1, j2
    type (t_Emission), intent(inOut) :: self
    Allocate(self%base_monot(i1:i2, ju1:j2, NTYPE))
    self%base_monot = 0.0d0
    return
  end subroutine Allocate_base_monot
!-------------------------------------------------------------------------
  subroutine Get_base_monot (self, base_monot)
    real*8          , intent(out)  :: base_monot (:,:,:)
    type (t_Emission), intent(in)   :: self
    base_monot(:,:,:) = self%base_monot(:,:,:)
    return
  end subroutine Get_base_monot
!-------------------------------------------------------------------------
  subroutine Set_base_monot (self, base_monot)
    real*8          , intent(in)  :: base_monot (:,:,:)
    type (t_Emission), intent(inOut) :: self
    self%base_monot(:,:,:) = base_monot(:,:,:)
    return
  end subroutine Set_base_monot
!-------------------------------------------------------------------------
  subroutine Allocate_xlai (self, i1, i2, ju1, j2)
    integer          , intent(in   ) :: i1, i2, ju1, j2
    type (t_Emission), intent(inOut) :: self
    Allocate(self%xlai(i1:i2, ju1:j2, NTYPE))
    self%xlai = 0.0d0
    return
  end subroutine Allocate_xlai
!-------------------------------------------------------------------------
  subroutine Get_xlai (self, xlai)
    real*8          , intent(out)  :: xlai (:,:,:)
    type (t_Emission), intent(in)   :: self
    xlai(:,:,:) = self%xlai(:,:,:)
    return
  end subroutine Get_xlai
!-------------------------------------------------------------------------
  subroutine Set_xlai (self, xlai)
    real*8          , intent(in)  :: xlai (:,:,:)
    type (t_Emission), intent(inOut) :: self
    self%xlai(:,:,:) = xlai(:,:,:)
    return
  end subroutine Set_xlai
!-------------------------------------------------------------------------
  subroutine Allocate_xlai2 (self, i1, i2, ju1, j2)
    integer          , intent(in   ) :: i1, i2, ju1, j2
    type (t_Emission), intent(inOut) :: self
    Allocate(self%xlai2(i1:i2, ju1:j2, NTYPE))
    self%xlai2 = 0.0d0
    return
  end subroutine Allocate_xlai2
!-------------------------------------------------------------------------
  subroutine Get_xlai2 (self, xlai2)
    real*8          , intent(out)  :: xlai2 (:,:,:)
    type (t_Emission), intent(in)   :: self
    xlai2(:,:,:) = self%xlai2(:,:,:)
    return
  end subroutine Get_xlai2
!-------------------------------------------------------------------------
  subroutine Set_xlai2 (self, xlai2)
    real*8          , intent(in)  :: xlai2 (:,:,:)
    type (t_Emission), intent(inOut) :: self
    self%xlai2(:,:,:) = xlai2(:,:,:)
    return
  end subroutine Set_xlai2
!-------------------------------------------------------------------------
  subroutine Get_coeff_isop (self, coeff_isop)
    real*8          , intent(out)  :: coeff_isop (:)
    type (t_Emission), intent(in)   :: self
    coeff_isop(:) = self%coeff_isop(:)
    return
  end subroutine Get_coeff_isop
!-------------------------------------------------------------------------
  subroutine Set_coeff_isop (self, coeff_isop)
    real*8          , intent(in)  :: coeff_isop (:)
    type (t_Emission), intent(inOut) :: self
    self%coeff_isop(:) = coeff_isop(:)
    return
  end subroutine Set_coeff_isop
!-------------------------------------------------------------------------
  subroutine Get_convert_isop (self, convert_isop)
    real*8          , intent(out)  :: convert_isop (:)
    type (t_Emission), intent(in)   :: self
    convert_isop(:) = self%convert_isop(:)
    return
  end subroutine Get_convert_isop
!-------------------------------------------------------------------------
  subroutine Set_convert_isop (self, convert_isop)
    real*8           , intent(in)  :: convert_isop (:)
    type (t_Emission), intent(inOut) :: self
    self%convert_isop(:) = convert_isop(:)
    return
  end subroutine Set_convert_isop
!-------------------------------------------------------------------------
  subroutine Get_convert_monot (self, convert_monot)
    real*8           , intent(out)  :: convert_monot (:)
    type (t_Emission), intent(in)   :: self
    convert_monot(:) = self%convert_monot(:)
    return
  end subroutine Get_convert_monot
!-------------------------------------------------------------------------
  subroutine Set_convert_monot (self, convert_monot)
    real*8           , intent(in)  :: convert_monot (:)
    type (t_Emission), intent(inOut) :: self
    self%convert_monot(:) = convert_monot(:)
    return
  end subroutine Set_convert_monot
!-------------------------------------------------------------------------
  subroutine Allocate_iland (self, i1, i2, ju1, j2)
    integer          , intent(in   ) :: i1, i2, ju1, j2
    type (t_Emission), intent(inOut) :: self
    Allocate(self%iland(i1:i2, ju1:j2, NTYPE))
    self%iland = 0
    return
  end subroutine Allocate_iland
!-------------------------------------------------------------------------
  subroutine Get_iland (self, iland)
    integer          , intent(out)  :: iland (:,:,:)
    type (t_Emission), intent(in)   :: self
    iland(:,:,:) = self%iland(:,:,:)
    return
  end subroutine Get_iland
!-------------------------------------------------------------------------
  subroutine Set_iland (self, iland)
    integer          , intent(in)  :: iland (:,:,:)
    type (t_Emission), intent(inOut) :: self
    self%iland(:,:,:) = iland(:,:,:)
    return
  end subroutine Set_iland
!-------------------------------------------------------------------------
  subroutine Allocate_iuse (self, i1, i2, ju1, j2)
    integer          , intent(in   ) :: i1, i2, ju1, j2
    type (t_Emission), intent(inOut) :: self
    Allocate(self%iuse(i1:i2, ju1:j2, NTYPE))
    self%iuse = 0
    return
  end subroutine Allocate_iuse
!-------------------------------------------------------------------------
  subroutine Get_iuse (self, iuse)
    integer          , intent(out)  :: iuse (:,:,:)
    type (t_Emission), intent(in)   :: self
    iuse(:,:,:) = self%iuse(:,:,:)
    return
  end subroutine Get_iuse
!-------------------------------------------------------------------------
  subroutine Set_iuse (self, iuse)
    integer          , intent(in)  :: iuse (:,:,:)
    type (t_Emission), intent(inOut) :: self
    self%iuse(:,:,:) = iuse(:,:,:)
    return
  end subroutine Set_iuse
!-------------------------------------------------------------------------
  subroutine Allocate_ireg (self, i1, i2, ju1, j2)
    integer          , intent(in   ) :: i1, i2, ju1, j2
    type (t_Emission), intent(inOut) :: self
    Allocate(self%ireg(i1:i2, ju1:j2))
    self%ireg = 0
    return
  end subroutine Allocate_ireg
!-------------------------------------------------------------------------
  subroutine Get_ireg (self, ireg)
    integer          , intent(out)  :: ireg (:,:)
    type (t_Emission), intent(in)   :: self
    ireg(:,:) = self%ireg(:,:)
    return
  end subroutine Get_ireg
!-------------------------------------------------------------------------
  subroutine Set_ireg (self, ireg)
    integer          , intent(in)  :: ireg (:,:)
    type (t_Emission), intent(inOut) :: self
    self%ireg(:,:) = ireg(:,:)
    return
  end subroutine Set_ireg
!-------------------------------------------------------------------------
  subroutine Get_ncon_soil (self, ncon_soil)
    integer          , intent(out)  :: ncon_soil (:)
    type (t_Emission), intent(in)   :: self
    ncon_soil(:) = self%ncon_soil(:)
    return
  end subroutine Get_ncon_soil
!-------------------------------------------------------------------------
  subroutine Set_ncon_soil (self, ncon_soil)
    integer          , intent(in)  :: ncon_soil (:)
    type (t_Emission), intent(inOut) :: self
    self%ncon_soil(:) = ncon_soil(:)
    return
  end subroutine Set_ncon_soil
!-------------------------------------------------------------------------
  subroutine Allocate_emissAero (self, i1, i2, ju1, j2)
    integer          , intent(in   ) :: i1, i2, ju1, j2
    type (t_Emission)       , intent(inOut) :: self
    Allocate(self%emissAero(i1:i2, ju1:j2, self%naero))
    self%emissAero = 0.0d0
    return
  end subroutine Allocate_emissAero
!-------------------------------------------------------------------------
  subroutine Get_emissAero (self, emiss_aero)
    real*8          , intent(out)  :: emiss_aero (:,:,:)
    type (t_Emission), intent(in)   :: self
    emiss_aero(:,:,:) = self%emissAero(:,:,:)
    return
  end subroutine Get_emissAero
!-------------------------------------------------------------------------
  subroutine Set_emissAero (self, emiss_aero)
    real*8          , intent(in)  :: emiss_aero (:,:,:)
    type (t_Emission), intent(inOut) :: self
    self%emissAero(:,:,:) = emiss_aero(:,:,:) 
    return
  end subroutine Set_emissAero
!-------------------------------------------------------------------------
  subroutine Allocate_emissDust (self, i1, i2, ju1, j2)
    integer          , intent(in   ) :: i1, i2, ju1, j2
    type (t_Emission)       , intent(inOut) :: self
    Allocate(self%emissDust(i1:i2, ju1:j2, 1:self%ndust))
    self%emissDust = 0.0d0
    return
  end subroutine Allocate_emissDust
!-------------------------------------------------------------------------
  subroutine Get_emissDust (self, emiss_dust)
    real*8           , intent(out)  :: emiss_dust (:,:,:)
    type (t_Emission), intent(in )  :: self
    emiss_dust(:,:,:) = self%emissDust(:,:,:)
    return
  end subroutine Get_emissDust
!-------------------------------------------------------------------------
  subroutine Set_emissDust (self, emiss_dust)
    real*8           , intent(in   ) :: emiss_dust (:,:,:)
    type (t_Emission), intent(inOut) :: self
    self%emissDust(:,:,:) = emiss_dust(:,:,:) 
    return
  end subroutine Set_emissDust
!-------------------------------------------------------------------------
  subroutine Allocate_lightning_NO (self, i1, i2, ju1, j2, k1, k2)
    integer          , intent(in   ) :: i1, i2, ju1, j2, k1, k2
    type (t_Emission)      , intent(inOut) :: self
    Allocate(self%lightning_NO(i1:i2, ju1:j2, k1:k2))
    self%lightning_NO = 0.0d0
    return
  end subroutine Allocate_lightning_NO
!-------------------------------------------------------------------------
  subroutine Get_lightning_NO (self, lightning_NO)
    real*8          , intent(out)  :: lightning_NO (:,:,:)
    type (t_Emission), intent(in)   :: self
    lightning_NO(:,:,:) = self%lightning_NO(:,:,:)
    return
  end subroutine Get_lightning_NO
!-------------------------------------------------------------------------
  subroutine Set_lightning_NO (self, lightning_NO)
    real*8          , intent(in)  :: lightning_NO (:,:,:)
    type (t_Emission), intent(inOut) :: self
    self%lightning_NO(:,:,:) = lightning_NO(:,:,:) 
    return
  end subroutine Set_lightning_NO
!-------------------------------------------------------------------------
  subroutine Allocate_emiss_isop (self, i1, i2, ju1, j2)
    integer          , intent(in   ) :: i1, i2, ju1, j2
    type (t_Emission)      , intent(inOut) :: self
    Allocate(self%emiss_isop(i1:i2, ju1:j2))
    self%emiss_isop = 0.0d0
    return
  end subroutine Allocate_emiss_isop
!-------------------------------------------------------------------------
  subroutine Get_emiss_isop (self, emiss_isop)
    real*8          , intent(out)  :: emiss_isop (:,:)
    type (t_Emission), intent(in)   :: self
    emiss_isop(:,:) = self%emiss_isop(:,:)
    return
  end subroutine Get_emiss_isop
!-------------------------------------------------------------------------
  subroutine Set_emiss_isop (self, emiss_isop)
    real*8          , intent(in)  :: emiss_isop (:,:)
    type (t_Emission), intent(inOut) :: self
    self%emiss_isop(:,:) = emiss_isop(:,:)
    return
  end subroutine Set_emiss_isop
!-------------------------------------------------------------------------
  subroutine Allocate_emiss_monot (self, i1, i2, ju1, j2)
    integer          , intent(in   ) :: i1, i2, ju1, j2
    type (t_Emission)      , intent(inOut) :: self
    Allocate(self%emiss_monot(i1:i2, ju1:j2))
    self%emiss_monot = 0.0d0
    return
  end subroutine Allocate_emiss_monot
!-------------------------------------------------------------------------
  subroutine Get_emiss_monot (self, emiss_monot)
    real*8          , intent(out)  :: emiss_monot (:,:)
    type (t_Emission), intent(in)   :: self
    emiss_monot(:,:) = self%emiss_monot(:,:)
    return
  end subroutine Get_emiss_monot
!-------------------------------------------------------------------------
  subroutine Set_emiss_monot (self, emiss_monot)
    real*8          , intent(in)  :: emiss_monot (:,:)
    type (t_Emission), intent(inOut) :: self
    self%emiss_monot(:,:) = emiss_monot(:,:)
    return
  end subroutine Set_emiss_monot
!-------------------------------------------------------------------------
  subroutine Allocate_emiss_nox (self, i1, i2, ju1, j2)
    integer          , intent(in   ) :: i1, i2, ju1, j2
    type (t_Emission)      , intent(inOut) :: self
    Allocate(self%emiss_nox(i1:i2, ju1:j2))
    self%emiss_nox = 0.0d0
    return
  end subroutine Allocate_emiss_nox
!-------------------------------------------------------------------------
  subroutine Get_emiss_nox (self, emiss_nox)
    real*8          , intent(out)  :: emiss_nox (:,:)
    type (t_Emission), intent(in)   :: self
    emiss_nox(:,:) = self%emiss_nox(:,:)
    return
  end subroutine Get_emiss_nox
!-------------------------------------------------------------------------
  subroutine Set_emiss_nox (self, emiss_nox)
    real*8          , intent(in)  :: emiss_nox (:,:)
    type (t_Emission), intent(inOut) :: self
    self%emiss_nox(:,:) = emiss_nox(:,:)
    return
  end subroutine Set_emiss_nox
!-------------------------------------------------------------------------
  subroutine Allocate_emiss_hno3 (self, i1, i2, ju1, j2)
    integer          , intent(in   ) :: i1, i2, ju1, j2
    type (t_Emission)      , intent(inOut) :: self
    Allocate(self%emiss_hno3(i1:i2, ju1:j2))
    self%emiss_hno3 = 0.0d0
    return
  end subroutine Allocate_emiss_hno3
!-------------------------------------------------------------------------
  subroutine Get_emiss_hno3 (self, emiss_hno3)
    real*8          , intent(out)  :: emiss_hno3 (:,:)
    type (t_Emission), intent(in)   :: self
    emiss_hno3(:,:) = self%emiss_hno3(:,:)
    return
  end subroutine Get_emiss_hno3
!-------------------------------------------------------------------------
  subroutine Set_emiss_hno3 (self, emiss_hno3)
    real*8          , intent(in)  :: emiss_hno3 (:,:)
    type (t_Emission), intent(inOut) :: self
    self%emiss_hno3(:,:) = emiss_hno3(:,:)
    return
  end subroutine Set_emiss_hno3
!-------------------------------------------------------------------------
  subroutine Allocate_emiss_o3 (self, i1, i2, ju1, j2)
    integer          , intent(in   ) :: i1, i2, ju1, j2
    type (t_Emission)      , intent(inOut) :: self
    Allocate(self%emiss_o3(i1:i2, ju1:j2))
    self%emiss_o3 = 0.0d0
    return
  end subroutine Allocate_emiss_o3
!-------------------------------------------------------------------------
  subroutine Get_emiss_o3 (self, emiss_o3)
    real*8          , intent(out)  :: emiss_o3 (:,:)
    type (t_Emission), intent(in)   :: self
    emiss_o3(:,:) = self%emiss_o3(:,:)
    return
  end subroutine Get_emiss_o3
!-------------------------------------------------------------------------
  subroutine Set_emiss_o3 (self, emiss_o3)
    real*8          , intent(in)  :: emiss_o3 (:,:)
    type (t_Emission), intent(inOut) :: self
    self%emiss_o3(:,:) = emiss_o3(:,:)
    return
  end subroutine Set_emiss_o3
!-------------------------------------------------------------------------
  subroutine Allocate_surf_emiss_out (self, i1, i2, ju1, j2, numSpecies)
    integer          , intent(in   ) :: i1, i2, ju1, j2, numSpecies
    type (t_Emission)      , intent(inOut) :: self
    Allocate(self%surf_emiss_out(i1:i2, ju1:j2, numSpecies))
    self%surf_emiss_out = 0.0d0
    return
  end subroutine Allocate_surf_emiss_out
!-------------------------------------------------------------------------
  subroutine Get_surf_emiss_out (self, surf_emiss_out)
    real*8          , intent(out)  :: surf_emiss_out (:,:,:)
    type (t_Emission), intent(in)   :: self
    surf_emiss_out(:,:,:) = self%surf_emiss_out(:,:,:)
    return
  end subroutine Get_surf_emiss_out
!-------------------------------------------------------------------------
  subroutine Set_surf_emiss_out (self, surf_emiss_out)
    real*8          , intent(in)  :: surf_emiss_out (:,:,:)
    type (t_Emission), intent(inOut) :: self
    self%surf_emiss_out(:,:,:) = surf_emiss_out(:,:,:)
    return
  end subroutine Set_surf_emiss_out
!-------------------------------------------------------------------------
  subroutine Allocate_surf_emiss_out2 (self, i1, i2, ju1, j2)
    integer          , intent(in   ) :: i1, i2, ju1, j2
    type (t_Emission)      , intent(inOut) :: self
    Allocate(self%surf_emiss_out2(i1:i2, ju1:j2, 6))
    self%surf_emiss_out2 = 0.0d0
    return
  end subroutine Allocate_surf_emiss_out2
!-------------------------------------------------------------------------
  subroutine Get_surf_emiss_out2 (self, surf_emiss_out2)
    real*8          , intent(out)  :: surf_emiss_out2 (:,:,:)
    type (t_Emission), intent(in)   :: self
    surf_emiss_out2(:,:,:) = self%surf_emiss_out2(:,:,:)
    return
  end subroutine Get_surf_emiss_out2
!-------------------------------------------------------------------------
  subroutine Set_surf_emiss_out2 (self, surf_emiss_out2)
    real*8          , intent(in)  :: surf_emiss_out2 (:,:,:)
    type (t_Emission), intent(inOut) :: self
    self%surf_emiss_out2(:,:,:) = surf_emiss_out2(:,:,:)
    return
  end subroutine Set_surf_emiss_out2
!-------------------------------------------------------------------------
  subroutine Allocate_emiss_3d_out (self, i1, i2, ju1, j2, k1, k2, numSpecies)
    integer          , intent(in   ) :: i1, i2, ju1, j2, k1, k2, numSpecies
    type (t_Emission), intent(inOut) :: self
    Allocate(self%emiss_3d_out(i1:i2, ju1:j2, k1:k2, numSpecies))
    self%emiss_3d_out = 0.0d0
    return
  end subroutine Allocate_emiss_3d_out
!-------------------------------------------------------------------------
  subroutine Get_emiss_3d_out (self, emiss_3d_out)
    real*8          , intent(out)  :: emiss_3d_out (:,:,:,:)
    type (t_Emission), intent(in)   :: self
    emiss_3d_out(:,:,:,:) = self%emiss_3d_out(:,:,:,:)
    return
  end subroutine Get_emiss_3d_out
!-------------------------------------------------------------------------
  subroutine Set_emiss_3d_out (self, emiss_3d_out)
    real*8          , intent(in)  :: emiss_3d_out (:,:,:,:)
    type (t_Emission), intent(inOut) :: self
    self%emiss_3d_out(:,:,:,:) = emiss_3d_out(:,:,:,:)
    return
  end subroutine Set_emiss_3d_out
!-------------------------------------------------------------------------
  subroutine Allocate_aerosolEmiss3D (self, i1, i2, ju1, j2, k1, k2)
    integer          , intent(in   ) :: i1, i2, ju1, j2, k1, k2
    type (t_Emission)      , intent(inOut) :: self
    Allocate(self%aerosolEmiss3D(i1:i2, ju1:j2, k1:k2, 5))
    self%aerosolEmiss3D = 0.0d0
    return
  end subroutine Allocate_aerosolEmiss3D
!-------------------------------------------------------------------------
  subroutine Get_aerosolEmiss3D (self, aerosolEmiss3D)
    real*8          , intent(out)  :: aerosolEmiss3D (:,:,:,:)
    type (t_Emission), intent(in)   :: self
    aerosolEmiss3D(:,:,:,:) = self%aerosolEmiss3D(:,:,:,:)
    return
  end subroutine Get_aerosolEmiss3D
!-------------------------------------------------------------------------
  subroutine Set_aerosolEmiss3D (self, aerosolEmiss3D)
    real*8          , intent(in)  :: aerosolEmiss3D (:,:,:,:)
    type (t_Emission), intent(inOut) :: self
    self%aerosolEmiss3D(:,:,:,:) = aerosolEmiss3D(:,:,:,:)
    return
  end subroutine Set_aerosolEmiss3D
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
  subroutine Get_aerosolSurfEmiss (self, aerosolSurfEmiss)
    real*8          , intent(out)  :: aerosolSurfEmiss (:,:,:)
    type (t_Emission), intent(in)   :: self
    aerosolSurfEmiss(:,:,:) = self%aerosolSurfEmiss(:,:,:)
    return
  end subroutine Get_aerosolSurfEmiss
!-------------------------------------------------------------------------
  subroutine Set_aerosolSurfEmiss (self, aerosolSurfEmiss)
    real*8          , intent(in)  :: aerosolSurfEmiss (:,:,:)
    type (t_Emission), intent(inOut) :: self
    self%aerosolSurfEmiss(:,:,:) = aerosolSurfEmiss(:,:,:)
    return
  end subroutine Set_aerosolSurfEmiss
!-------------------------------------------------------------------------
  subroutine Allocate_aerosolSurfEmissMap (self)
    type (t_Emission)      , intent(inOut) :: self
    Allocate(self%aerosolSurfEmissMap &
        (self%ndust + self%naero + 5))
    return
  end subroutine Allocate_aerosolSurfEmissMap
!-------------------------------------------------------------------------
  subroutine Get_aerosolSurfEmissMap (self, aerosolSurfEmissMap)
    integer          , intent(out)  :: aerosolSurfEmissMap (:)
    type (t_Emission), intent(in)   :: self
    aerosolSurfEmissMap(:) = self%aerosolSurfEmissMap(:)
    return
  end subroutine Get_aerosolSurfEmissMap
!-------------------------------------------------------------------------
  subroutine Allocate_flashrate (self, i1, i2, ju1, j2)
    integer          , intent(in   ) :: i1, i2, ju1, j2
    type (t_Emission)      , intent(inOut) :: self
    Allocate(self%flashrate(i1:i2, ju1:j2))
    self%flashrate = 0.0d0
    return
  end subroutine Allocate_flashrate
!-------------------------------------------------------------------------
  subroutine Get_flashrate (self, flashrate)
    real*8          , intent(out)  :: flashrate (:,:)
    type (t_Emission), intent(in)   :: self
    flashrate(:,:) = self%flashrate(:,:)
    return
  end subroutine Get_flashrate
!-------------------------------------------------------------------------
  subroutine Set_flashrate (self, flashrate)
    real*8          , intent(in)  :: flashrate (:,:)
    type (t_Emission), intent(inOut) :: self
    self%flashrate(:,:) = flashrate(:,:) 
    return
  end subroutine Set_flashrate
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
  subroutine Get_emissDust_t (self, emiss_dust_t)
    real*8          , intent(out)  :: emiss_dust_t (:,:,:,:)
    type (t_Emission), intent(in)   :: self
    emiss_dust_t(:,:,:,:) = self%emissDust_t(:,:,:,:)
    return
  end subroutine Get_emissDust_t
!-------------------------------------------------------------------------
  subroutine Set_emissDust_t (self, emiss_dust_t)
    real*8          , intent(in)  :: emiss_dust_t (:,:,:,:)
    type (t_Emission), intent(inOut) :: self
    self%emissDust_t(:,:,:,:) = emiss_dust_t(:,:,:,:)
    return
  end subroutine Set_emissDust_t
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
  subroutine Get_emissAero_t (self, emiss_aero_t)
    real*8          , intent(out)  :: emiss_aero_t (:,:,:,:)
    type (t_Emission), intent(in)   :: self
    emiss_aero_t(:,:,:,:) = self%emissAero_t(:,:,:,:)
    return
  end subroutine Get_emissAero_t
!-------------------------------------------------------------------------
  subroutine Set_emissAero_t (self, emiss_aero_t)
    real*8          , intent(in)  :: emiss_aero_t (:,:,:,:)
    type (t_Emission), intent(inOut) :: self
    self%emissAero_t(:,:,:,:) = emiss_aero_t(:,:,:,:)
    return
  end subroutine Set_emissAero_t
!-------------------------------------------------------------------------
  subroutine Get_emiss_map (self, emiss_map, numSpecies)
    integer        , intent(in )  :: numSpecies
    integer        , intent(out)  :: emiss_map (:)
    type (t_Emission), intent(in)   :: self
    emiss_map(1:numSpecies) = self%emiss_map(1:numSpecies)
    return
  end subroutine Get_emiss_map
!-------------------------------------------------------------------------
  subroutine Get_emiss_timpyr (self, emiss_timpyr)
    integer        , intent(out)  :: emiss_timpyr 
    type (t_Emission), intent(in)   :: self
    emiss_timpyr = self%emiss_timpyr
    return
  end subroutine Get_emiss_timpyr
!-------------------------------------------------------------------------
  subroutine Set_emiss_timpyr (self, emiss_timpyr)
    integer        , intent(in)  :: emiss_timpyr 
    type (t_Emission), intent(inOut) :: self
    self%emiss_timpyr = emiss_timpyr
    return
  end subroutine Set_emiss_timpyr
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
  subroutine Get_do_ShipEmission (self, do_ShipEmission)
    logical        , intent(out)  :: do_ShipEmission
    type (t_Emission), intent(in)   :: self
    do_ShipEmission = self%do_ShipEmission
    return
  end subroutine Get_do_ShipEmission
!-------------------------------------------------------------------------
  subroutine Get_doMEGANemission (self, doMEGANemission)
    logical        , intent(out)  :: doMEGANemission
    type (t_Emission), intent(in)   :: self
    doMEGANemission = self%doMEGANemission
    return
  end subroutine Get_doMEGANemission
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
  subroutine Get_aefMbo(self, aefMbo)
    real*8 , intent(out) :: aefMbo(:,:)
    type (t_Emission), intent(in)   :: self
    aefMbo(:,:) = self%aefMbo(:,:)
    return
  end subroutine Get_aefMbo
!-------------------------------------------------------------------------
  subroutine Allocate_aefIsop(self, i1, i2, ju1, j2)
    integer, intent(in) :: i1, i2, ju1, j2
    type (t_Emission), intent(inOut)   :: self
    allocate(self%aefIsop(i1:i2,ju1:j2))
    self%aefIsop = 0.0d0
    return
  end subroutine Allocate_aefIsop
!-------------------------------------------------------------------------
  subroutine Get_aefIsop(self, aefIsop)
    real*8 , intent(out) :: aefIsop(:,:)
    type (t_Emission), intent(in)   :: self
    aefIsop(:,:) = self%aefIsop(:,:)
    return
  end subroutine Get_aefIsop
!-------------------------------------------------------------------------
  subroutine Allocate_aefOvoc(self, i1, i2, ju1, j2)
    integer, intent(in) :: i1, i2, ju1, j2
    type (t_Emission), intent(inOut)   :: self
    allocate(self%aefOvoc(i1:i2,ju1:j2))
    self%aefOvoc = 0.0d0
    return
  end subroutine Allocate_aefOvoc
!-------------------------------------------------------------------------
  subroutine Get_aefOvoc(self, aefOvoc)
    real*8 , intent(out) :: aefOvoc(:,:)
    type (t_Emission), intent(in)   :: self
    aefOvoc(:,:) = self%aefOvoc(:,:)
    return
  end subroutine Get_aefOvoc
!-------------------------------------------------------------------------
  subroutine Allocate_aefMonot(self, i1, i2, ju1, j2)
    integer, intent(in) :: i1, i2, ju1, j2
    type (t_Emission), intent(inOut)   :: self
    allocate(self%aefMonot(i1:i2,ju1:j2))
    self%aefMonot = 0.0d0
    return
  end subroutine Allocate_aefMonot
!-------------------------------------------------------------------------
  subroutine Get_aefMonot(self, aefMonot)
    real*8 , intent(out) :: aefMonot(:,:)
    type (t_Emission), intent(in)   :: self
    aefMonot(:,:) = self%aefMonot(:,:)
    return
  end subroutine Get_aefMonot
!-------------------------------------------------------------------------
  subroutine Get_emiss_in_opt (self, emiss_in_opt)
    integer        , intent(out)  :: emiss_in_opt 
    type (t_Emission), intent(in)   :: self
    emiss_in_opt = self%emiss_in_opt
    return
  end subroutine Get_emiss_in_opt
!-------------------------------------------------------------------------
  subroutine Get_emiss_opt (self, emiss_opt)
    integer        , intent(out)  :: emiss_opt 
    type (t_Emission), intent(in)   :: self
    emiss_opt = self%emiss_opt
    return
  end subroutine Get_emiss_opt
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
  subroutine Get_o3_index (self, o3_index)
    integer        , intent(out)  :: o3_index
    type (t_Emission), intent(in)   :: self
    o3_index = self%o3_index
    return
  end subroutine Get_o3_index
!-------------------------------------------------------------------------
  subroutine Get_hno3_index (self, hno3_index)
    integer        , intent(out)  :: hno3_index
    type (t_Emission), intent(in)   :: self
    hno3_index = self%hno3_index
    return
  end subroutine Get_hno3_index
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
  subroutine Get_nst_dust (self, nst_dust)
    integer        , intent(out)  :: nst_dust
    type (t_Emission), intent(in)   :: self
    nst_dust = self%nst_dust
    return
  end subroutine Get_nst_dust
!-------------------------------------------------------------------------
  subroutine Get_nt_dust (self, nt_dust)
    integer        , intent(out)  :: nt_dust
    type (t_Emission), intent(in)   :: self
    nt_dust = self%nt_dust
    return
  end subroutine Get_nt_dust
!-------------------------------------------------------------------------
  subroutine Get_emiss_map_dust (self, emiss_map_dust, ndust)
    integer        , intent(in )  :: ndust
    integer        , intent(out)  :: emiss_map_dust(:)
    type (t_Emission), intent(in)   :: self
    emiss_map_dust(1:ndust) = self%emiss_map_dust(1:ndust)
    return
  end subroutine Get_emiss_map_dust
!-------------------------------------------------------------------------
  subroutine Get_emiss_map_aero (self, emiss_map_aero, naero)
    integer        , intent(in )  :: naero
    integer        , intent(out)  :: emiss_map_aero(:)
    type (t_Emission), intent(in)   :: self
    emiss_map_aero(1:naero) = self%emiss_map_aero(1:naero)
    return
  end subroutine Get_emiss_map_aero
!-------------------------------------------------------------------------
  subroutine Get_fertscal_infile_name (self, fertscal_infile_name)
    character (len=MAX_LENGTH_FILE_NAME), intent(out)  :: fertscal_infile_name
    type (t_Emission)    , intent(in)   :: self
    fertscal_infile_name = self%fertscal_infile_name
    return
  end subroutine Get_fertscal_infile_name
!-------------------------------------------------------------------------
  subroutine Get_lai_infile_name (self, lai_infile_name)
    character (len=MAX_LENGTH_FILE_NAME), intent(out)  :: lai_infile_name
    type (t_Emission)    , intent(in)   :: self
    lai_infile_name = self%lai_infile_name
    return
  end subroutine Get_lai_infile_name
!-------------------------------------------------------------------------
  subroutine Get_precip_infile_name (self, precip_infile_name)
    character (len=MAX_LENGTH_FILE_NAME), intent(out)  :: precip_infile_name
    type (t_Emission)    , intent(in)   :: self
    precip_infile_name = self%precip_infile_name
    return
  end subroutine Get_precip_infile_name
!-------------------------------------------------------------------------
  subroutine Get_soil_infile_name (self, soil_infile_name)
    character (len=MAX_LENGTH_FILE_NAME), intent(out)  :: soil_infile_name
    type (t_Emission)    , intent(in)   :: self
    soil_infile_name = self%soil_infile_name
    return
  end subroutine Get_soil_infile_name
!-------------------------------------------------------------------------
  subroutine Get_veg_infile_name (self, veg_infile_name)
    character (len=MAX_LENGTH_FILE_NAME), intent(out)  :: veg_infile_name
    type (t_Emission)    , intent(in)   :: self
    veg_infile_name = self%veg_infile_name
    return
  end subroutine Get_veg_infile_name
!-------------------------------------------------------------------------
  subroutine Get_isopconv_infile_name (self, isopconv_infile_name)
    character (len=MAX_LENGTH_FILE_NAME), intent(out)  :: isopconv_infile_name
    type (t_Emission)    , intent(in)   :: self
    isopconv_infile_name = self%isopconv_infile_name
    return
  end subroutine Get_isopconv_infile_name
!-------------------------------------------------------------------------
  subroutine Get_monotconv_infile_name (self, monotconv_infile_name)
    character (len=MAX_LENGTH_FILE_NAME), intent(out)  :: monotconv_infile_name
    type (t_Emission)    , intent(in)   :: self
    monotconv_infile_name = self%monotconv_infile_name
    return
  end subroutine Get_monotconv_infile_name
!-------------------------------------------------------------------------
  subroutine Get_GOCARTerod_infile_name (self, GOCARTerod_infile_name)
    character (len=MAX_LENGTH_FILE_NAME), intent(out)  :: GOCARTerod_infile_name
    type (t_Emission)    , intent(in)   :: self
    GOCARTerod_infile_name = self%GOCARTerod_infile_name
    return
  end subroutine Get_GOCARTerod_infile_name
!-------------------------------------------------------------------------
  subroutine Get_GOCARTocean_infile_name (self, GOCARTocean_infile_name)
    character (len=MAX_LENGTH_FILE_NAME), intent(out)  :: GOCARTocean_infile_name
    type (t_Emission)    , intent(in)   :: self
    GOCARTocean_infile_name = self%GOCARTocean_infile_name
    return
  end subroutine Get_GOCARTocean_infile_name
!-------------------------------------------------------------------------
  subroutine Get_GOCARTerod_mod_infile_name (self, GOCARTerod_mod_infile_name)
    character (len=MAX_LENGTH_FILE_NAME), intent(out)  :: GOCARTerod_mod_infile_name
    type (t_Emission)    , intent(in)   :: self
    GOCARTerod_mod_infile_name = self%GOCARTerod_mod_infile_name
    return
  end subroutine Get_GOCARTerod_mod_infile_name
!-------------------------------------------------------------------------
  subroutine Get_lightNOampFactor (self, lightNOampFactor)
    implicit none
    real           , intent(out)  :: lightNOampFactor
    type (t_Emission), intent(in)   :: self
    lightNOampFactor = self%lightNOampFactor
    return
  end subroutine Get_lightNOampFactor
!-------------------------------------------------------------------------
  subroutine Get_numberNOperFlash (self, numberNOperFlash)
    implicit none
    real           , intent(out)  :: numberNOperFlash
    type (t_Emission), intent(in)   :: self
    numberNOperFlash = self%numberNOperFlash
    return
  end subroutine Get_numberNOperFlash
!-------------------------------------------------------------------------
  subroutine Get_minDeepCloudTop (self, minDeepCloudTop)
    implicit none
    real           , intent(out)  :: minDeepCloudTop
    type (t_Emission), intent(in)   :: self
    minDeepCloudTop = self%minDeepCloudTop
    return
  end subroutine Get_minDeepCloudTop
!-------------------------------------------------------------------------
  subroutine Get_isop_scale (self, isop_scale)
    real*8         , intent(out)  :: isop_scale(:)
    type (t_Emission), intent(in)   :: self
    isop_scale(:) = self%isop_scale(:)
    return
  end subroutine Get_isop_scale
!-------------------------------------------------------------------------
  subroutine Set_isop_scale (self, isop_scale)
    real*8         , intent(in   )  :: isop_scale(:)
    type (t_Emission), intent(inOut)  :: self
    self%isop_scale(:) = isop_scale(:)
    return
  end subroutine Set_isop_scale
!-------------------------------------------------------------------------
  subroutine Get_begDailyEmissRec (self, begDailyEmissRec)
    integer          , intent(out)  :: begDailyEmissRec 
    type (t_Emission), intent(in)   :: self
    begDailyEmissRec = self%begDailyEmissRec
    return
  end subroutine Get_begDailyEmissRec
!-------------------------------------------------------------------------
  subroutine Get_endDailyEmissRec (self, endDailyEmissRec)
    integer          , intent(out)  :: endDailyEmissRec 
    type (t_Emission), intent(in)   :: self
    endDailyEmissRec = self%endDailyEmissRec
    return
  end subroutine Get_endDailyEmissRec
!-------------------------------------------------------------------------
  subroutine Get_doReadDailyEmiss (self, doReadDailyEmiss)
    logical          , intent(out)  :: doReadDailyEmiss 
    type (t_Emission), intent(in)   :: self
    doReadDailyEmiss = self%doReadDailyEmiss
    return
  end subroutine Get_doReadDailyEmiss
!-------------------------------------------------------------------------
  end module GmiEmissionMethod_mod
