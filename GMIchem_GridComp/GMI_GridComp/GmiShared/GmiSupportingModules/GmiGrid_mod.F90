!-------------------------------------------------------------------------
!BOP
!
! !MODULE: GmiGrid_mod
!
! !INTERFACE:
!
  module GmiGrid_mod
!
! !USES:
!
  implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
!
  private
  public  :: InitializeGmiGrid
  public  :: Get_NPIJ_actm, Get_NPI_actm, Get_NPJ_actm
  public  :: Get_i1, Get_i2, Get_ju1, Get_jv1, Get_j2, Get_k1, Get_k2
  public  :: Get_i1_gl, Get_i2_gl, Get_ju1_gl, Get_jv1_gl, Get_j2_gl, Get_k1_gl, Get_k2_gl
  public  :: Get_ilo, Get_ihi, Get_julo, Get_jvlo, Get_jhi
  public  :: Get_ilo_gl, Get_ihi_gl, Get_julo_gl, Get_jvlo_gl, Get_jhi_gl
  public  :: Get_ilong, Get_ilat, Get_ivert, Get_itloop
  public  :: Get_j1p, Get_j2p

! !PUBLIC DATA MEMBERS:

  public  :: t_gmiGrid

  type t_gmiGrid
!    private
    integer :: NPIJ_actm   ! Total number of computing processors 
    integer :: NPI_actm    ! number of processors in the i (i.e., longitude) direction
    integer :: NPJ_actm    ! number of processors in the j (i.e., latitude)  direction

    integer :: i1_gl       ! index of first global longitude      (no ghost zones)
    integer :: i2_gl       ! index of last  global longitude      (no ghost zones)
    integer :: ju1_gl      ! index of first global "u"   latitude (no ghost zones)
    integer :: jv1_gl      ! index of first global "v"   latitude (no ghost zones)
    integer :: j2_gl       ! index of last  global "u&v" latitude (no ghost zones)
    integer :: k1_gl       ! index of first global altitude       (no ghost zones)
    integer :: k2_gl       ! index of last  global altitude       (no ghost zones)

    integer :: ilo_gl      ! i1_gl  - gmi_nborder (has ghost zones)
    integer :: ihi_gl      ! i2_gl  + gmi_nborder (has ghost zones)
    integer :: julo_gl     ! ju1_gl - gmi_nborder (has ghost zones)
    integer :: jvlo_gl     ! jv1_gl - gmi_nborder (has ghost zones)
    integer :: jhi_gl      ! j2_gl  + gmi_nborder (has ghost zones)

    integer :: i1          ! index of first local longitude      (no ghost zones)
    integer :: i2          ! index of last  local longitude      (no ghost zones)
    integer :: ju1         ! index of first local "u"   latitude (no ghost zones)
    integer :: jv1         ! index of first local "v"   latitude (no ghost zones)
    integer :: j2          ! index of last  local "u&v" latitude (no ghost zones)
    integer :: k1          ! index of first local altitude       (no ghost zones)
    integer :: k2          ! index of last  local altitude       (no ghost zones)

    integer :: ilo         ! i1  - gmi_nborder (has ghost zones)
    integer :: ihi         ! i2  + gmi_nborder (has ghost zones)
    integer :: julo        ! ju1 - gmi_nborder (has ghost zones)
    integer :: jvlo        ! jv1 - gmi_nborder (has ghost zones)
    integer :: jhi         ! j2  + gmi_nborder (has ghost zones)

    integer :: ilat        ! number of latitudes
    integer :: ilong       ! number of longitudes
    integer :: ivert       ! number of vertical layers
    integer :: itloop      ! number of zones (ilat * ilong * ivert)

    integer :: gmi_nborder ! number of longitude and latitude ghost zones

    integer :: j1p         ! determines size of the Polar cap
    integer :: j2p         ! j2_gl - j1p + 1
  end type t_gmiGrid

! !DESCRIPTION:
!  Domain decomposition information.
!
! !AUTHOR:
!  Jules Kouatchou
!
! !REVISION HISTORY:
! 5April2007 Initial implementation.
!
!EOP
!-------------------------------------------------------------------------

  CONTAINS

!-------------------------------------------------------------------------
  subroutine Get_i1 (self, i1)
    integer         , intent(out)  :: i1
    type (t_gmiGrid), intent(in)   :: self
    i1 = self%i1
    return
  end subroutine Get_i1
!------------------------------------------------------------------------
  subroutine Get_i2 (self, i2)
    integer         , intent(out)  :: i2
    type (t_gmiGrid), intent(in)   :: self
    i2 = self%i2
    return
  end subroutine Get_i2
!------------------------------------------------------------------------
  subroutine Get_ju1 (self, ju1)
    integer         , intent(out)  :: ju1
    type (t_gmiGrid), intent(in)   :: self
    ju1 = self%ju1
    return
  end subroutine Get_ju1
!------------------------------------------------------------------------
  subroutine Get_jv1 (self, jv1)
    integer         , intent(out)  :: jv1
    type (t_gmiGrid), intent(in)   :: self
    jv1 = self%jv1
    return
  end subroutine Get_jv1
!------------------------------------------------------------------------
  subroutine Get_j2 (self, j2)
    integer         , intent(out)  :: j2
    type (t_gmiGrid), intent(in)   :: self
    j2 = self%j2
    return
  end subroutine Get_j2
!------------------------------------------------------------------------
  subroutine Get_k1 (self, k1)
    integer         , intent(out)  :: k1
    type (t_gmiGrid), intent(in)   :: self
    k1 = self%k1
    return
  end subroutine Get_k1
!------------------------------------------------------------------------
  subroutine Get_k2 (self, k2)
    integer         , intent(out)  :: k2
    type (t_gmiGrid), intent(in)   :: self
    k2 = self%k2
    return
  end subroutine Get_k2
!------------------------------------------------------------------------
  subroutine Get_i1_gl (self, i1_gl)
    integer         , intent(out)  :: i1_gl
    type (t_gmiGrid), intent(in)   :: self
    i1_gl = self%i1_gl
    return
  end subroutine Get_i1_gl
!------------------------------------------------------------------------
  subroutine Get_i2_gl (self, i2_gl)
    integer         , intent(out)  :: i2_gl
    type (t_gmiGrid), intent(in)   :: self
    i2_gl = self%i2_gl
    return
  end subroutine Get_i2_gl
!------------------------------------------------------------------------
  subroutine Get_ju1_gl (self, ju1_gl)
    integer         , intent(out)  :: ju1_gl
    type (t_gmiGrid), intent(in)   :: self
    ju1_gl = self%ju1_gl
    return
  end subroutine Get_ju1_gl
!------------------------------------------------------------------------
  subroutine Get_jv1_gl (self, jv1_gl)
    integer         , intent(out)  :: jv1_gl
    type (t_gmiGrid), intent(in)   :: self
    jv1_gl = self%jv1_gl
    return
  end subroutine Get_jv1_gl
!------------------------------------------------------------------------
  subroutine Get_j2_gl (self, j2_gl)
    integer         , intent(out)  :: j2_gl
    type (t_gmiGrid), intent(in)   :: self
    j2_gl = self%j2_gl
    return
  end subroutine Get_j2_gl
!------------------------------------------------------------------------
  subroutine Get_k1_gl (self, k1_gl)
    integer         , intent(out)  :: k1_gl
    type (t_gmiGrid), intent(in)   :: self
    k1_gl = self%k1_gl
    return
  end subroutine Get_k1_gl
!------------------------------------------------------------------------
  subroutine Get_k2_gl (self, k2_gl)
    integer         , intent(out)  :: k2_gl
    type (t_gmiGrid), intent(in)   :: self
    k2_gl = self%k2_gl
    return
  end subroutine Get_k2_gl
!------------------------------------------------------------------------
  subroutine Get_ilo (self, ilo)
    integer         , intent(out)  :: ilo
    type (t_gmiGrid), intent(in)   :: self
    ilo = self%ilo
    return
  end subroutine Get_ilo
!------------------------------------------------------------------------
  subroutine Get_ihi (self, ihi)
    integer         , intent(out)  :: ihi
    type (t_gmiGrid), intent(in)   :: self
    ihi = self%ihi
    return
  end subroutine Get_ihi
!------------------------------------------------------------------------
  subroutine Get_julo (self, julo)
    integer         , intent(out)  :: julo
    type (t_gmiGrid), intent(in)   :: self
    julo = self%julo
    return
  end subroutine Get_julo
!------------------------------------------------------------------------
  subroutine Get_jvlo (self, jvlo)
    integer         , intent(out)  :: jvlo
    type (t_gmiGrid), intent(in)   :: self
    jvlo = self%jvlo
    return
  end subroutine Get_jvlo
!------------------------------------------------------------------------
  subroutine Get_jhi (self, jhi)
    integer         , intent(out)  :: jhi
    type (t_gmiGrid), intent(in)   :: self
    jhi = self%jhi
    return
  end subroutine Get_jhi
!------------------------------------------------------------------------
  subroutine Get_ilo_gl (self, ilo_gl)
    integer         , intent(out)  :: ilo_gl
    type (t_gmiGrid), intent(in)   :: self
    ilo_gl = self%ilo_gl
    return
  end subroutine Get_ilo_gl
!------------------------------------------------------------------------
  subroutine Get_ihi_gl (self, ihi_gl)
    integer         , intent(out)  :: ihi_gl
    type (t_gmiGrid), intent(in)   :: self
    ihi_gl = self%ihi_gl
    return
  end subroutine Get_ihi_gl
!------------------------------------------------------------------------
  subroutine Get_julo_gl (self, julo_gl)
    integer         , intent(out)  :: julo_gl
    type (t_gmiGrid), intent(in)   :: self
    julo_gl = self%julo_gl
    return
  end subroutine Get_julo_gl
!------------------------------------------------------------------------
  subroutine Get_jvlo_gl (self, jvlo_gl)
    integer         , intent(out)  :: jvlo_gl
    type (t_gmiGrid), intent(in)   :: self
    jvlo_gl = self%jvlo_gl
    return
  end subroutine Get_jvlo_gl
!------------------------------------------------------------------------
  subroutine Get_jhi_gl (self, jhi_gl)
    integer         , intent(out)  :: jhi_gl
    type (t_gmiGrid), intent(in)   :: self
    jhi_gl = self%jhi_gl
    return
  end subroutine Get_jhi_gl
!------------------------------------------------------------------------
  subroutine Get_ilong (self, ilong)
    integer         , intent(out)  :: ilong
    type (t_gmiGrid), intent(in)   :: self
    ilong = self%ilong
    return
  end subroutine Get_ilong
!------------------------------------------------------------------------
  subroutine Get_ilat (self, ilat)
    integer         , intent(out)  :: ilat
    type (t_gmiGrid), intent(in)   :: self
    ilat = self%ilat
    return
  end subroutine Get_ilat
!------------------------------------------------------------------------
  subroutine Get_ivert (self, ivert)
    integer         , intent(out)  :: ivert
    type (t_gmiGrid), intent(in)   :: self
    ivert = self%ivert
    return
  end subroutine Get_ivert
!------------------------------------------------------------------------
  subroutine Get_itloop (self, itloop)
    integer         , intent(out)  :: itloop
    type (t_gmiGrid), intent(in)   :: self
    itloop = self%itloop
    return
  end subroutine Get_itloop
!------------------------------------------------------------------------
  subroutine Get_j1p (self, j1p)
    integer         , intent(out)  :: j1p
    type (t_gmiGrid), intent(in)   :: self
    j1p = self%j1p
    return
  end subroutine Get_j1p
!------------------------------------------------------------------------
  subroutine Get_j2p (self, j2p)
    integer         , intent(out)  :: j2p
    type (t_gmiGrid), intent(in)   :: self
    j2p = self%j2p
    return
  end subroutine Get_j2p
!------------------------------------------------------------------------
  subroutine Get_NPIJ_actm (self, NPIJ_actm)
    integer         , intent(out)  :: NPIJ_actm
    type (t_gmiGrid), intent(in)   :: self
    NPIJ_actm = self%NPIJ_actm
    return
  end subroutine Get_NPIJ_actm
!------------------------------------------------------------------------
  subroutine Get_NPI_actm (self, NPI_actm)
    integer         , intent(out)  :: NPI_actm
    type (t_gmiGrid), intent(in)   :: self
    NPI_actm = self%NPI_actm
    return
  end subroutine Get_NPI_actm
!------------------------------------------------------------------------
  subroutine Get_NPJ_actm (self, NPJ_actm)
    integer         , intent(out)  :: NPJ_actm
    type (t_gmiGrid), intent(in)   :: self
    NPJ_actm = self%NPJ_actm
    return
  end subroutine Get_NPJ_actm
!------------------------------------------------------------------------
  subroutine InitializeGmiGrid (self, &
               NPIJ_actm, NPI_actm, NPJ_actm, gmi_nborder, &
               i1, i2, ju1, jv1, j2, k1, k2, &
               i1_gl, i2_gl, ju1_gl, jv1_gl, j2_gl, k1_gl, k2_gl, &
               ilo, ihi, julo, jvlo, jhi, &
               ilo_gl, ihi_gl, julo_gl, jvlo_gl, jhi_gl, &
               ilong, ilat, ivert, itloop, j1p, j2p)  

  implicit none

  integer          , intent(in   ) :: NPIJ_actm, NPI_actm, NPJ_actm, gmi_nborder
  integer          , intent(in   ) :: i1, i2, ju1, jv1, j2, k1, k2
  integer          , intent(in   ) :: i1_gl, i2_gl, ju1_gl, jv1_gl, j2_gl, k1_gl, k2_gl
  integer          , intent(in   ) :: ilo, ihi, julo, jvlo, jhi
  integer          , intent(in   ) :: ilo_gl, ihi_gl, julo_gl, jvlo_gl, jhi_gl
  integer          , intent(in   ) :: ilong, ilat, ivert, itloop, j1p, j2p

  type  (t_gmiGrid), intent(inout) :: self

  self%NPI_actm   = NPI_actm
  self%NPJ_actm   = NPJ_actm
  self%NPIJ_actm  = NPIJ_actm
  self%gmi_nborder= gmi_nborder
  self%i1         = i1
  self%i2         = i2
  self%ju1        = ju1
  self%jv1        = jv1
  self%j2         = j2
  self%k1         = k1
  self%k2         = k2
  self%i1_gl      = i1_gl
  self%i2_gl      = i2_gl
  self%ju1_gl     = ju1_gl
  self%jv1_gl     = jv1_gl
  self%j2_gl      = j2_gl
  self%k1_gl      = k1_gl
  self%k2_gl      = k2_gl
  self%ilo        = ilo
  self%ihi        = ihi
  self%julo       = julo
  self%jvlo       = jvlo
  self%jhi        = jhi
  self%ilo_gl     = ilo_gl
  self%ihi_gl     = ihi_gl
  self%julo_gl    = julo_gl
  self%jvlo_gl    = jvlo_gl
  self%jhi_gl     = jhi_gl
  self%ilong      = ilong
  self%ilat       = ilat
  self%ivert      = ivert
  self%itloop     = itloop
  self%j1p        = j1p
  self%j2p        = j2p

  return

  end subroutine InitializeGmiGrid

!-------------------------------------------------------------------------


  end module GmiGrid_mod
