#include "MAPL_Generic.h"
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: GmiDepositionMethod_mod 
!
! !INTERFACE:
!
  module GmiDepositionMethod_mod
!
! !USES:
      use Esmf
      use MAPL
      use GmiGrid_mod                  , only : t_gmiGrid, Get_i1, Get_i2
      use GmiGrid_mod                  , only : Get_ju1, Get_j2, Get_k1, Get_k2
      use GmiGrid_mod                  , only : Get_ilo, Get_ihi
      use GmiGrid_mod                  , only : Get_julo, Get_jhi
      use GmiGrid_mod                  , only : Get_ilong, Get_ivert
      use GmiPrintError_mod            , only : GmiPrintError
      use GmiEmissionMethod_mod        , only : t_Emission
      use GmiESMFrcFileReading_mod     , only : rcEsmfReadTable
      use GmiArrayBundlePointer_mod    , only : t_GmiArrayBundle
      use GmiSpcConcentrationMethod_mod, only : t_SpeciesConcentration
      use GmiSpcConcentrationMethod_mod, only : Get_concentration
      use GmiSpcConcentrationMethod_mod, only : Set_concentration
!
  implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
!
      private 
      public  :: InitializeDeposition
      public  :: RunDryDeposition
      public  :: RunWetDeposition
!     public  :: RunSimpleDeposition
      public  :: FinalizeDeposition

      public  :: Get_num_ks_sdep, Get_do_drydep, Get_do_wetdep
      public  :: Get_do_simpledep, Get_wetdep_eff
      public  :: Get_scav3d, Get_dry_depos, Get_wet_depos
      public  :: Set_scav3d, Set_dry_depos, Set_wet_depos
!
! !PUBLIC DATA MEMBERS:
!
  public  :: t_Deposition

#     include "GmiParameters.h"
#     include "gmi_emiss_constants.h"

      type t_Deposition
        !private
        integer :: num_ks_sdep    ! number of vertical layers to apply 2 
                                  ! day loss factor to in simple deposition
        logical :: do_drydep      ! do dry    deposition?
        logical :: do_wetdep      ! do wet    deposition?
        logical :: do_simpledep   ! do simple deposition?
        real*8  :: wetdep_eff(MAX_NUM_CONST)
                                  ! array of wet deposition (scavenging) efficiencies
        real*8, pointer :: dry_depos(:,:,:) => null()
        real*8, pointer :: wet_depos(:,:,:) => null()
        real*8, pointer ::  scav3d(:,:,:,:) => null()
      end type t_Deposition
!
! !DESCRIPTION:
! Initialization and run methods for the Deposition component.
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
! !IROUTINE: readDepositionResourceFile
!
! !INTERFACE:
!
      subroutine readDepositionResourceFile (self, procID, config, pr_diag)
!
      implicit none
!
! !INPUT PARAMETERS:
      integer, intent(in) :: procID
      logical, intent(in) :: pr_diag
!
! !INPUT/OUTPUT PARAMETERS:
      type (ESMF_Config), intent(inOut) :: config
      type (t_Deposition), intent(inOut) :: self
!
! !DESCRIPTION:
! Reads the Deposition section of the resource file.
!
! !LOCAL VARIABLES:
      integer :: STATUS, RC
      character(len=ESMF_MAXSTR) :: err_msg, IAm
!EOP
!-----------------------------------------------------------------------------
!BOC
      Iam = "readDepositionResourceFile"

      if (pr_diag) Write(6,*) IAm, 'called by ', procID

     !-------------------------------------
     ! Reading Deposition Related Variables
     !-------------------------------------

      call ESMF_ConfigGetAttribute(config, self%num_ks_sdep, &
     &                label   = "num_ks_sdep:",&
     &                default = 1, rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, value=self%do_drydep,    &
     &               label="do_drydep:",    default=.false., __RC__ )

      call ESMF_ConfigGetAttribute(config, value=self%do_wetdep,    &
     &               label="do_wetdep:",    default=.false., __RC__ )

      call ESMF_ConfigGetAttribute(config, value=self%do_simpledep, &
     &               label="do_simpledep:", default=.false., __RC__ )

     ! ----------------------------------------------------------------
     ! wetdep_eff : wet deposition (scavenging) efficiencies; should be
     !              set to values between 0.0 and 1.0 for each species
     ! ----------------------------------------------------------------
      
      self%wetdep_eff(:) = 0.0d0
      
      call rcEsmfReadTable(config, self%wetdep_eff, "wetdep_eff::", rc=STATUS)
      !VERIFY_(STATUS)

      !##############################
      ! End reading the resource file
      !##############################
      
      if (self%do_simpledep .and. self%do_drydep) then
         err_msg = 'do_simpledep/do_drydep problem in Check_Nlvalue.'
         call GmiPrintError (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
      end if

      if (self%do_simpledep .and. self%do_wetdep) then
         err_msg = 'do_simpledep/do_wetdep problem in Check_Nlvalue.'
         call GmiPrintError (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
      end if
      
      return         
      
      end subroutine readDepositionResourceFile
!EOC 
!-------------------------------------------------------------------------------
!BOC
!
! !ROUTINE: InitializeDeposition
!
! !INTERFACE
!
      subroutine InitializeDeposition (self, gmiGrid, config, numSpecies, &
                           loc_proc, pr_diag, pr_dry_depos, pr_wet_depos, pr_scav)
!
      implicit none
!
! !INPUT PARAMETERS:
      logical, intent(in) :: pr_diag, pr_dry_depos, pr_wet_depos, pr_scav
  integer         , intent(in) :: numSpecies
  integer         , intent(in) :: loc_proc
  type (t_gmiGrid), intent(in) :: gmiGrid   
!
! !INPUT/OUTPUT PARAMETERS:
      type (ESMF_Config), intent(inOut) :: config
      type (t_Deposition), intent(inOut) :: self
!
! !DESCRIPTION:
!  Initialize the Deposition component.
!
! !LOCAL VARIABLES
      integer :: i1, i2, ju1, j2, k1, k2
!
!EOP
!------------------------------------------------------------------------------
!BOC
      call Get_i1 (gmiGrid, i1 )
      call Get_i2 (gmiGrid, i2 )
      call Get_ju1(gmiGrid, ju1)
      call Get_j2 (gmiGrid, j2 )
      call Get_k1 (gmiGrid, k1 )
      call Get_k2 (gmiGrid, k2 )

      call readDepositionResourceFile (self, loc_proc, config, pr_diag)

      if (pr_dry_depos) then
         call Allocate_dry_depos(self, i1, i2, ju1, j2, numSpecies)
      end if

      if (pr_wet_depos) then
         call Allocate_wet_depos(self, i1, i2, ju1, j2, numSpecies)
      end if

      if (pr_scav) then
         call Allocate_scav3d(self, i1, i2, ju1, j2, k1, k2, numSpecies)
      end if

      return

      end subroutine InitializeDeposition
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: RunDryDeposition
!
! !INTERFACE:
!
      subroutine RunDryDeposition (self, Emission, SpeciesConcentration,       &
     &              gmiGrid, lwi_flags, mcor, cosSolarZenithAngle,             &
     &              fracCloudCover, radswg, surf_air_temp, surf_rough, ustar,  &
     &              mass, diffaer, s_radius, s_velocity, BoxHeightEdge,        &
     &              loc_proc, mw, numSpecies, chem_opt, pr_dry_depos, pr_diag, &
     &              tdt4)
!
      implicit none
!
! !INPUT PARAMETERS:
      logical,              intent(in) :: pr_diag
      logical,              intent(in) :: pr_dry_depos
      integer,              intent(in) :: loc_proc
      integer,              intent(in) :: chem_opt
      integer,              intent(in) :: numSpecies
      real*8 ,              intent(in) :: mw(numSpecies)
      integer,              intent(in) :: lwi_flags (:,:)
      real*8 ,              intent(in) :: mcor(:,:)
      real*8 ,              intent(in) :: mass(:,:,:)
      real*8 ,              intent(in) :: cosSolarZenithAngle (:,:)
      real*8 ,              intent(in) :: fracCloudCover (:,:)
      real*8 ,              intent(in) :: radswg(:,:)
      real*8 ,              intent(in) :: surf_air_temp(:,:)
      real*8 ,              intent(in) :: surf_rough(:,:)
      real*8 ,              intent(in) :: ustar(:,:)
      real*8 ,              intent(in) :: BoxHeightEdge  (:,:)
      real*8 ,              intent(in) :: diffaer   (:,:,:)
      real*8 ,              intent(in) :: s_radius  (:,:,:)
      real*8 ,              intent(in) :: s_velocity(:,:,:)
      type(t_gmiGrid),      intent(in) :: gmiGrid
      real   ,              intent(in) :: tdt4
!
! !INPUT/OUTPUT PRAMETERS:
      type(t_Deposition),           intent(inOut) :: self
      type(t_Emission),             intent(inOut) :: Emission
      type(t_SpeciesConcentration), intent(inOut) :: SpeciesConcentration
!
! !DESCRIPTION:
! Run method for the Dry Deposition component
!
! !LOCAL VARIABLES:
      real*8 , allocatable :: dry_depos (:,:,:)
      real*8  :: tdt8
      integer :: i1, i2, ju1, j2, k1, k2, ilong
      integer :: ilo, ihi, julo, jhi
      type (t_GmiArrayBundle), pointer :: concentration(:)
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
      call Get_ilo   (gmiGrid, ilo  )
      call Get_ihi   (gmiGrid, ihi  )
      call Get_julo  (gmiGrid, julo )
      call Get_jhi   (gmiGrid, jhi  )
      call Get_ilong (gmiGrid, ilong)

      tdt8 = tdt4

      call Get_concentration(SpeciesConcentration, concentration)

      call Update_Drydep (pr_dry_depos, lwi_flags, mcor, cosSolarZenithAngle,  &
     &           fracCloudCover, radswg, surf_air_temp, surf_rough, ustar,     &
     &           mass, concentration, self%dry_depos, diffaer, s_radius,       &
     &           s_velocity, BoxHeightEdge, Emission%ireg, Emission%iland,     &
     &           Emission%iuse, Emission%xlai, pr_diag,                        &
     &           loc_proc, chem_opt, tdt8, i1, i2, ju1, j2, k1, k2, ilong, ilo,&
     &           ihi, julo, jhi, mw, numSpecies)

      if (pr_dry_depos) then
         where (self%dry_depos(:,:,:) < 1.0d-30)
              self%dry_depos(:,:,:) = 1.0d-30
         end where
      end if

      call Set_concentration(SpeciesConcentration, concentration)

      return

      end subroutine RunDryDeposition
!EOC
!-------------------------------------------------------------------------
!BOP  
!     
! !ROUTINE: RunWetDeposition
!     
! !INTERFACE:
!     
      subroutine RunWetDeposition (self, SpeciesConcentration,                 &
     &              gmiGrid,  ih2o2_num, ihno3_num, mw, con_precip,            &
     &              tot_precip, mcor, gridBoxHeight, mass, moistq, rain_cn,    &
     &              rain_ls, kel, press3c, press3e, loc_proc, numSpecies,      &
     &              chem_opt, pr_wet_depos, pr_scav, pr_diag, tdt4)
!
      implicit none
!
! !INPUT PARAMETERS:
      integer           , intent(in) :: loc_proc
      logical           , intent(in) :: pr_wet_depos, pr_scav, pr_diag
      integer           , intent(in) :: chem_opt
      integer           , intent(in) :: numSpecies
      integer           , intent(in) :: ih2o2_num
      integer           , intent(in) :: ihno3_num
      real*8            , intent(in) :: mw(numSpecies)
      real*8            , intent(in) :: con_precip (:,:)
      real*8            , intent(in) :: tot_precip(:,:)
      real*8            , intent(in) :: mcor(:,:)
      real*8            , intent(in) :: mass(:,:,:)
      real*8            , intent(in) :: gridBoxHeight(:,:,:)
      real*8            , intent(in) :: moistq(:,:,:)
      real*8            , intent(in) :: rain_cn(:,:,:)
      real*8            , intent(in) :: rain_ls(:,:,:)
      real*8            , intent(in) :: kel(:,:,:)
      real*8            , intent(in) :: press3c(:,:,:)
      real*8            , intent(in) :: press3e(:,:,:)
      type(t_gmiGrid   ), intent(in) :: gmiGrid   
      real              , intent(in) :: tdt4

!
! !INPUT/OUTPUT PARAMETERS:
      type(t_Deposition),           intent(inOut) :: self
      type(t_SpeciesConcentration), intent(inOut) :: SpeciesConcentration
!
! !DESCRIPTION:
! Run method for the Wet Deposition component.
! 
! !LOCAL VARIABLES:
      real*8  :: tdt8
      integer :: i1, i2, ju1, j2, k1, k2, ilong, ivert
      integer :: ilo, ihi, julo, jhi
      type (t_GmiArrayBundle), pointer :: concentration(:)
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
      call Get_ilo   (gmiGrid, ilo  )
      call Get_ihi   (gmiGrid, ihi  )
      call Get_julo  (gmiGrid, julo )
      call Get_jhi   (gmiGrid, jhi  )
      call Get_ilong (gmiGrid, ilong)
      call Get_ivert (gmiGrid, ivert)

      tdt8 = tdt4

      call Get_concentration(SpeciesConcentration, concentration)

      call Update_Wetdep (pr_wet_depos, pr_scav, chem_opt, ih2o2_num,          &
     &           ihno3_num, tdt8, mw, con_precip, tot_precip, mcor,            &
     &           gridBoxHeight, mass, moistq, rain_cn, rain_ls, kel, press3c,  &
     &           press3e, concentration, self%wet_depos, self%scav3d, pr_diag, &
     &           loc_proc, i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi,       &
     &           ilong, ivert, numSpecies)

      if (pr_wet_depos) then
         where (self%wet_depos(:,:,:) < 1.0d-30)
              self%wet_depos(:,:,:) = 1.0d-30
         end where
      end if

      call Set_concentration(SpeciesConcentration, concentration)

      return

      end subroutine RunWetDeposition
!EOC
!-------------------------------------------------------------------------
!!BOP
!!
!! !ROUTINE: RunSimpleDeposition
!!
!! !INTERFACE:
!!
!      subroutine RunSimpleDeposition (self, SpeciesConcentration, gmiClock,    &
!     &              gmiGrid, ibrono2_num, ih2o2_num, ihcl_num, ihno3_num,      &
!     &              imgas_num, initrogen_num, ioxygen_num, press3c, loc_proc,  &
!     &              numSpecies, pr_diag)
!!
!      implicit none
!!
!! !INPUT PARAMETERS:
!      logical,           intent(in) :: pr_diag
!      integer,           intent(in) :: loc_proc, numSpecies
!      integer,           intent(in) :: ibrono2_num, ih2o2_num
!      integer,           intent(in) :: ihcl_num, ihno3_num, imgas_num
!      integer,           intent(in) :: initrogen_num, ioxygen_num
!      real*8 ,           intent(in) :: press3c(:,:,:)
!      type (t_GmiClock), intent(in) :: gmiClock
!      type (t_gmiGrid ), intent(in) :: gmiGrid 
!!
!! !INPUT/OUTPUT VARIABLES:
!      type (t_Deposition),          intent(inOut) :: self
!      type(t_SpeciesConcentration), intent(inOut) :: SpeciesConcentration
!!
!! !DESCRIPTION:
!!  Run method for the Simple Deposition component.
!!
!! !LOCAL VARIABLES:
!      real*8  :: tdt
!      integer :: i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi
!      type (t_GmiArrayBundle), pointer :: concentration(:)
!!
!!EOP
!!------------------------------------------------------------------------------
!!BOC
!      ! Get the GMI grid information
!      call Get_i1    (gmiGrid, i1   )
!      call Get_i2    (gmiGrid, i2   )
!      call Get_ju1   (gmiGrid, ju1  )
!      call Get_j2    (gmiGrid, j2   )
!      call Get_k1    (gmiGrid, k1   )
!      call Get_k2    (gmiGrid, k2   )
!      call Get_ilo   (gmiGrid, ilo  )
!      call Get_ihi   (gmiGrid, ihi  )
!      call Get_julo  (gmiGrid, julo )
!      call Get_jhi   (gmiGrid, jhi  )
!
!      ! Obtain the model time step
!      call Get_gmiTimeStep(gmiClock, tdt)
!
!      call Get_concentration(SpeciesConcentration, concentration)
!
!      call Update_Simpledep (ibrono2_num, ih2o2_num, ihcl_num, ihno3_num,      &
!     &           imgas_num, initrogen_num, ioxygen_num, self%num_ks_sdep, tdt, &
!     &           press3c, concentration, pr_diag, loc_proc, i1, i2, ju1, j2,   &
!     &           k1, k2, ilo, ihi, julo, jhi, numSpecies)
!
!      call Set_concentration(SpeciesConcentration, concentration)
!
!      return
!
!      end subroutine RunSimpleDeposition
!!EOC
!-------------------------------------------------------------------------

  subroutine FinalizeDeposition (self)

  type (t_Deposition)   , intent(inout) :: self

  PRINT*,'  Finalize Deposition'

  return

  end subroutine FinalizeDeposition

!-------------------------------------------------------------------------
  subroutine Allocate_dry_depos (self, i1, i2, ju1, j2, numSpecies)
    integer            , intent(in)  :: i1, i2, ju1, j2, numSpecies
    type (t_Deposition), intent(inout)   :: self
    allocate(self%dry_depos(i1:i2, ju1:j2, 1:numSpecies))
    self%dry_depos = 0.0d0
    return
  end subroutine Allocate_dry_depos
!-------------------------------------------------------------------------
  subroutine Allocate_wet_depos (self, i1, i2, ju1, j2, numSpecies)
    integer            , intent(in)  :: i1, i2, ju1, j2, numSpecies
    type (t_Deposition), intent(inout)   :: self
    allocate(self%wet_depos(i1:i2, ju1:j2, 1:numSpecies))
    self%wet_depos = 0.0d0
    return
  end subroutine Allocate_wet_depos
!-------------------------------------------------------------------------
  subroutine Get_dry_depos (self, dry_depos)
    real*8             , intent(out)  :: dry_depos(:,:,:)
    type (t_Deposition), intent(in)   :: self
    dry_depos(:,:,:) = self%dry_depos(:,:,:)
    return
  end subroutine Get_dry_depos
!-------------------------------------------------------------------------
  subroutine Get_wet_depos (self, wet_depos)
    real*8             , intent(out)  :: wet_depos(:,:,:)
    type (t_Deposition), intent(in)   :: self
    wet_depos(:,:,:) = self%wet_depos(:,:,:)
    return
  end subroutine Get_wet_depos
!-------------------------------------------------------------------------
  subroutine Set_dry_depos (self, dry_depos)
    real*8             , intent(in)    :: dry_depos(:,:,:)
    type (t_Deposition), intent(inout) :: self
    self%dry_depos(:,:,:) = dry_depos(:,:,:)
    return
  end subroutine Set_dry_depos
!-------------------------------------------------------------------------
  subroutine Set_wet_depos (self, wet_depos)
    real*8             , intent(in)    :: wet_depos(:,:,:)
    type (t_Deposition), intent(inout) :: self
    self%wet_depos(:,:,:) = wet_depos(:,:,:)
    return
  end subroutine Set_wet_depos
!-------------------------------------------------------------------------
  subroutine Allocate_scav3d (self, i1, i2, ju1, j2, k1, k2, numSpecies)
    integer            , intent(in)  :: i1, i2, ju1, j2, k1, k2, numSpecies
    type (t_Deposition), intent(inout)   :: self
    allocate(self%scav3d(i1:i2, ju1:j2, k1:k2, 1:numSpecies))
    self%scav3d = 0.0d0
    return
  end subroutine Allocate_scav3d
!-------------------------------------------------------------------------
  subroutine Set_scav3d (self, scav3d)
    real*8             , intent(in)    :: scav3d(:,:,:,:)
    type (t_Deposition), intent(inout) :: self
    self%scav3d(:,:,:,:) = scav3d(:,:,:,:)
    return
  end subroutine Set_scav3d
!-------------------------------------------------------------------------
  subroutine Get_scav3d (self, scav3d)
    real*8             , intent(out)  :: scav3d(:,:,:,:)
    type (t_Deposition), intent(in)   :: self
    scav3d(:,:,:,:) = self%scav3d(:,:,:,:)
    return
  end subroutine Get_scav3d
!-------------------------------------------------------------------------
  subroutine Get_num_ks_sdep (self, num_ks_sdep)
    integer         , intent(out)  :: num_ks_sdep
    type (t_Deposition), intent(in)   :: self
    num_ks_sdep = self%num_ks_sdep
    return
  end subroutine Get_num_ks_sdep
!-------------------------------------------------------------------------
  subroutine Get_do_drydep (self, do_drydep)
    logical          , intent(out)  :: do_drydep
    type (t_Deposition), intent(in )  :: self
    do_drydep = self%do_drydep
    return
  end subroutine Get_do_drydep
!-------------------------------------------------------------------------
  subroutine Get_do_wetdep (self, do_wetdep)
    logical          , intent(out)  :: do_wetdep
    type (t_Deposition), intent(in )  :: self
    do_wetdep = self%do_wetdep
    return
  end subroutine Get_do_wetdep
!-------------------------------------------------------------------------
  subroutine Get_do_simpledep (self, do_simpledep)
    logical          , intent(out)  :: do_simpledep
    type (t_Deposition), intent(in )  :: self
    do_simpledep = self%do_simpledep
    return
  end subroutine Get_do_simpledep
!-------------------------------------------------------------------------
  subroutine Get_wetdep_eff (self, wetdep_eff)
    real*8           , intent(out)  :: wetdep_eff(:)
    type (t_Deposition), intent(in )  :: self
    wetdep_eff(:) = self%wetdep_eff(:) 
    return
  end subroutine Get_wetdep_eff
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
  end module GmiDepositionMethod_mod
