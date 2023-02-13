#include "MAPL_Generic.h"
!-------------------------------------------------------------------------
!NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1, GEOS/DAS!
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  CARMA_GridCompMod --- CARMA Grid Component Class
!
! Grid Component class for the Community Aerosol and Radiation
! Model for Atmospheres aerosol/cloud microphysics package.
!
! !INTERFACE:
!

   MODULE  CARMA_GridCompMod

! !USES:

   USE ESMF
   USE MAPL
   USE Chem_Mod
   USE Chem_UtilMod
   USE m_inpak90	     ! Resource file management

!  Utility Modules
   use DryDepositionMod      ! Aerosol Dry Deposition
   use WetRemovalMod         ! Aerosol Wet Removal
   use DustEmissionMod, only: KokSizeDistribution


!  CARMA Specific Methods
   use carma_precision_mod 
   use carma_constants_mod 
   use carma_enums_mod 
   use carma_types_mod 
   use carmaelement_mod
   use carmagroup_mod
   use carmagas_mod
   use carmastate_mod
   use carma_mod

   IMPLICIT NONE
   INTEGER, PARAMETER :: DBL = KIND(0.00D+00)

! !TYPES:

   PRIVATE
   PUBLIC  CARMA_GridComp       ! The CARMA object 
   PUBLIC  CARMA_Registry

! !PUBLIC MEMBER FUNCTIONS:

   PUBLIC  CARMA_GridCompInitialize
   PUBLIC  CARMA_GridCompRun
   PUBLIC  CARMA_GridCompFinalize
   PRIVATE dumpGas
   PRIVATE dumpElement
   PRIVATE dumpGroup

!
! !DESCRIPTION:
!
!  This module implements the CARMA aerosol & cloud microphysics model
!
! !REVISION HISTORY:
!
!  16Sep2003 da Silva  First crack.
!  24Jan2O05 Nielsen   Implementation of Code 916 chemistry
!  19Dec2005 d Silva   Minor portability mods.
!  30Oct2007 Nielsen   GMI Combo set up
!  18May2009 Colarco   Follow GMI setup to implement CARMA
!
!EOP
!-------------------------------------------------------------------------

  TYPE CARMA_GridComp
   CHARACTER(LEN=255) :: name = "CARMA aerosol/cloud microphysics"
   type(CARMA_Registry), pointer :: CARMAreg => null()
   type(carma_type), pointer     :: carma
   type(Chem_Mie), pointer       :: CARMAmie     ! GOCART style Mie lookup tables
   integer :: i1 = 1, i2, im, j1 = 1, j2, jm, km
   type(ESMF_grid) :: grid
   real, pointer, dimension(:,:) :: LONS, LATS
   integer :: nymd_bc = 1

!  Pointers to species specific emissions

!  Dust
   real, pointer, dimension(:,:) :: dust_source => null()

!  Smoke

!  Sulfate
   real, pointer :: vLat(:)    => null(), &
                    vLon(:)    => null(), &
                    vSO2(:)    => null(), &
                    vElev(:)   => null(), &
                    vCloud(:)  => null()

! Component derived type declarations
! -----------------------------------
!   TYPE(t_Chemistry )		:: Chemistry
 
  END TYPE CARMA_GridComp

  TYPE CARMA_Registry
     logical                     :: doing_CARMA = .false.
     integer                     :: nq
     character(len=255), pointer :: vname(:)  ! variable name (groupname::elemname::XXX)
     CHARACTER(LEN=255)          :: rcfilen = 'CARMAchem_Registry.rc'
     integer                     :: NBIN, NGROUP, NELEM, NSOLUTE, NGAS, NWAVE
     REAL(kind=f), pointer       :: RMRAT(:)        =>null(), &
                                    RMIN(:)         =>null(), &
                                    RHOP(:)         =>null(), &
                                    ESHAPE(:)       =>null(), &
                                    FSCAV(:)        =>null()
     INTEGER, pointer            :: IGROUP(:)       =>null(), &
                                    IRHSWELL(:)     =>null(), &
                                    IRHSWCOMP(:)    =>null(), &
                                    ISHAPE(:)       =>null(), &
                                    ICOMPOSITION(:) =>null(), &
                                    ITYPE(:)        =>null()
     character(len=255), pointer :: GROUPNAME(:)    =>null(), &
                                    ELEMNAME(:)     =>null()

!    Gases
     character(len=255), pointer :: GASNAME(:)      => null()
     integer, pointer            :: IGCOMP(:)       => null(), &
                                    IGVAPREQ(:)     => null()

     logical :: do_cnst_rlh = .false.
     logical :: do_coag = .false.       !! do coagulation?
     logical :: do_detrain = .false.
     logical :: do_fixedinit = .false.
     logical :: do_grow = .false.       !! do nucleation, growth and evaporation?
     logical :: do_incloud = .false.
     logical :: do_explised = .false.
     logical :: do_print_init = .false.
     logical :: do_substep = .false.    !! do substepping
     logical :: do_thermo = .false.     !! do thermodynamics
     logical :: do_vdiff = .false.      !! do Brownin diffusion
     logical :: do_vtran  = .false.     !! do sedimentation
     real(kind=f)  :: vf_const = 0._f   !! if specified and non-zero, constant
                                        !! fall velocity for all particles [cm/s]
     integer :: minsubsteps = 1         !! minimum number of substeps, default = 1
     integer :: maxsubsteps = 32        !! maximum number of substeps, default = 32
     integer :: maxretries = 16         !! maximum number of substep retries, default = 16
     real(kind=f)  :: conmax  = 0.1_f   !! minimum relative concentration to 
                                        !! consider, default = 1e-1
!    Species specific information
     integer :: igrp_mixed    = -1      !! mixed group
     integer :: mixedcorecomp = -1      !! mixed core element (sulfate)
     integer :: igrp_sulfate  = -1      !! sulfate group
     integer :: ielm_sulfate  = -1      !! sulfate pc element
     integer :: igrp_dust  = -1         !! dust group
     integer :: ielm_dust  = -1         !! dust pc element
     integer :: igrp_seasalt  = -1      !! seasalt group
     integer :: ielm_seasalt  = -1      !! seasalt pc element
     integer :: igrp_smoke  = -1        !! smoke group
     integer :: ielm_smoke  = -1        !! smoke pc element
     integer :: igrp_black_carbon  = -1 !! black carbon group
     integer :: ielm_black_carbon  = -1 !! black carbon pc element
     integer :: igrp_ash  = -1          !! ash group
     integer :: ielm_ash  = -1          !! ash pc element
     integer :: ielm_mxpc     = -1      !! mixed group pc element
     integer :: ielm_mxsulfate= -1      !! mixed group sulfate core element
     integer :: ielm_mxdust   = -1      !! mixed group dust core element
     integer :: ielm_mxsmoke  = -1      !! mixed group smoke core element
     integer :: ielm_mxseasalt= -1      !! mixed group seasalt core element
     integer :: ielm_mxbc     = -1      !! mixed group black carbon core element
     integer :: ielm_mxash    = -1      !! mixed group ash core element
     integer :: igas_h2o      = -1      !! water vapor
     integer :: igas_h2so4    = -1      !! sulfuric acid gas
     integer :: igas_hno3     = -1      !! nitric acid gas
     integer :: ifallrtn      =  1      !! default fall velocity routine for particles

!    Dust
     real               :: dust_emissions_fudgefactor
     real, pointer      :: dmass_dust(:) => null()  !! dust emission size distribution

!    Sea Salt
     real               :: seasalt_emissions_fudgefactor

!    Smoke
     real               :: organic_matter_to_organic_carbon_ratio
     real               :: fraction_terpene_to_organic_carbon

!    Black Carbon

!    Sulfate
     logical :: gmi_chem_provider=.false. ! direct exchange of h2so4, hno3, ...

!    GOCART-style Mie Lookup Tables
     integer            :: nchannels, nmoments
     real, pointer      :: channels(:)
     character(len=255) :: du_optics_file
     character(len=255) :: ss_optics_file
     character(len=255) :: bc_optics_file
     character(len=255) :: sm_optics_file
     character(len=255) :: su_optics_file

!    Workspace for any requested point emissions
!    Sulfate
     logical :: doing_point_emissions_sulfate=.FALSE.         ! Providing pointwise emissions
     character(len=255) :: point_emissions_srcfilen_sulfate   ! filename for pointwise emissions
     integer                         :: nPts_sulfate = -1
     integer, pointer, dimension(:)  :: vstart_sulfate => null(), &
                                        vend_sulfate   => null()
     real, pointer, dimension(:)     :: vLat_sulfate   => null(), &
                                        vLon_sulfate   => null(), &
                                        vBase_sulfate  => null(), &
                                        vTop_sulfate   => null(), &
                                        vEmis_sulfate  => null()
!    Ash
     logical :: doing_point_emissions_ash=.FALSE.         ! Providing pointwise emissions
     character(len=255) :: point_emissions_srcfilen_ash   ! filename for pointwise emissions
     integer                         :: nPts_ash = -1
     integer, pointer, dimension(:)  :: vstart_ash => null(), &
                                        vend_ash   => null()
     real, pointer, dimension(:)     :: vLat_ash   => null(), &
                                        vLon_ash   => null(), &
                                        vBase_ash  => null(), &
                                        vTop_ash   => null(), &
                                        vEmis_ash  => null()
!    Dust
     logical :: doing_point_emissions_dust=.FALSE.         ! Providing pointwise emissions
     character(len=255) :: point_emissions_srcfilen_dust   ! filename for pointwise emissions
     integer                         :: nPts_dust = -1
     integer, pointer, dimension(:)  :: vstart_dust => null(), &
                                        vend_dust   => null()
     real, pointer, dimension(:)     :: vLat_dust   => null(), &
                                        vLon_dust   => null(), &
                                        vBase_dust  => null(), &
                                        vTop_dust   => null(), &
                                        vEmis_dust  => null()


  END TYPE CARMA_Registry

CONTAINS

!-------------------------------------------------------------------------
!NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1, GEOS/DAS!
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  CARMA_GridCompInitialize --- Initialize CARMA_GridComp
!
! !INTERFACE:
!

   SUBROUTINE CARMA_GridCompInitialize( gcCARMA, impChem, expChem, nymd, nhms, cdt, &
                                        rc )

   IMPLICIT none
   INTEGER, PARAMETER :: DBL = KIND(0.00D+00)

! !INPUT PARAMETERS:

   INTEGER, INTENT(IN) :: nymd, nhms		       ! Time from AGCM
   REAL,    INTENT(IN) :: cdt			       ! Chemistry time step (secs)

! !OUTPUT PARAMETERS:

   TYPE(CARMA_GridComp), INTENT(INOUT) :: gcCARMA    ! Grid Component
   TYPE(ESMF_State),   INTENT(INOUT)   :: impChem    ! Import State
   TYPE(ESMF_State),   INTENT(INOUT)   :: expChem    ! Export State

   INTEGER, INTENT(out) ::  rc                  ! Error return code:
                                                !  0 - all is well
                                                !  1 - 

! !DESCRIPTION: Initializes the CARMA Grid Component. It primarily sets
!               the import state.
!
! !REVISION HISTORY:
!
!  18Sep2003 da Silva  First crack.
!  30Jun2007 Nielsen   GMI Combo set up
!  18May2009 Colarco   Adapt to use for CARMA
!
!EOP
!-------------------------------------------------------------------------

   CHARACTER(LEN=*), PARAMETER :: IAm = 'CARMA_GridCompInitialize'
   CHARACTER(LEN=255) :: rcfilen = 'CARMAchem_Registry.rc'
   CHARACTER(LEN=255) :: string

   INTEGER :: ios, n
   INTEGER, ALLOCATABLE :: ier(:)
   INTEGER :: i, i1, i2, ic, im, j, j1, j2, jm, km
   INTEGER :: nbins, n1, n2
   INTEGER :: STATUS

   INTEGER :: NBIN, NGROUP, NELEM, NSOLUTE, NWAVE, NGAS
   REAL(kind=f), allocatable    :: rmrat(:), rmin(:), rhop(:), &
                                   eshape(:), ishape(:),       &
                                   radius_(:), rlow_(:), rup_(:), &
                                   rhod_(:), rhog_(:)
   real, allocatable            :: radius(:),  rlow(:),  rup(:)
   REAL(kind=f)                 :: gwtmol
   REAL                         :: rhod, rhog
   INTEGER, allocatable :: IGROUP(:)
   character(len=255)              :: groupname, elemname, gasname
   type(carmagroup_type)               :: cgroup
   type(carmaelement_type)             :: celement
   type(carmastate_type), allocatable  :: cstate(:)
   logical :: do_coag = .false.       !! do coagulation?
   logical :: do_grow = .false.       !! do nucleation, growth and evaporation?
   logical :: do_implised = .false.   !! do sedimentation with substepping
   logical :: do_substep = .false.    !! do substepping
   logical :: do_thermo = .false.     !! do thermodynamics
   logical :: do_vdiff = .false.      !! do Brownin diffusion
   logical :: do_vtran  = .false.     !! do sedimentation
   logical :: is_sulfate = .false.    !! special handling for sulfate aerosol case
   real(kind=f) :: vf_const = 0._f    !! if specified and non-zero, constant fall velocity for all particles [cm/s]
   integer :: minsubsteps = 1  !! minimum number of substeps, default = 1
   integer :: maxsubsteps = 1  !! maximum number of substeps, default = 1
   integer :: maxretries = 5   !! maximum number of substep retries, default = 5
   real(kind=f) :: conmax  = 0.1_f      !! minimum relative concentration to consider, default = 1e-1
   type(CARMA_Registry), pointer :: reg => null()
   type(carma_type), pointer     :: r => null()


! This is for initializing condensed H2O to zero, for now
! -------------------------------------------------------
   REAL, ALLOCATABLE :: h2ocond(:,:,:)
   REAL, ALLOCATABLE :: cellArea(:,:)

   reg => gcCARMA%CARMAreg
   gcCARMA%name = 'CARMA aerosol/cloud microphysics'

!  Initialize local variables
!  --------------------------
   rc = 0

   CALL init_()
   IF ( rc /= 0 ) RETURN
   ier(:)=0

!!  Check on the parameters and if they agree with the Chem_Registry
!!  ----------------------------------------------------------------
!!  n_CARMA = NELEM*NBIN + NGAS
!   if(nbins .ne. (nbin*nelem + ngas) ) then
!    call final_(25)
!    return
!   endif

!  Establish the CARMA structure
!  -----------------------------
   allocate(gcCARMA%carma, stat=ios)
   if(ios /= 0) then
    call final_(103)
    return
   endif
   if(MAPL_AM_I_ROOT()) then
    call CARMA_Create(gcCARMA%carma, reg%NBIN, reg%NELEM, reg%NGROUP, &
                                     reg%NSOLUTE, reg%NGAS, reg%NWAVE, rc, &
                                     LUNOPRT=6)
   else 
    call CARMA_Create(gcCARMA%carma, reg%NBIN, reg%NELEM, reg%NGROUP, &
                                     reg%NSOLUTE, reg%NGAS, reg%NWAVE, rc )
   endif

   if (rc /=0) then
    call final_(rc)
    return
   endif

   r   => gcCARMA%carma


!  Establish the groups
!  --------------------
!  NOTE: Hard coded optionals and parameters here
   do j = 1, reg%NGROUP
    is_sulfate = .false.
!   Assumes MIXEDP has sulfate as PC
    if(  ESMF_UtilStringUpperCase(trim(reg%groupname(j))) == 'SULFATE' .or. &
         ESMF_UtilStringUpperCase(trim(reg%groupname(j))) == 'MIXEDP' ) is_sulfate = .true.
    if(reg%ishape(j) .eq. 1) then
     reg%ifallrtn = I_FALLRTN_STD
    else
     reg%ifallrtn = I_FALLRTN_STD_SHAPE
    endif
    call CARMAGROUP_Create(r, j, reg%groupname(j), reg%rmin(j), reg%rmrat(j), &
                           reg%ishape(j), reg%eshape(j), .FALSE., rc, ifallrtn=reg%ifallrtn, &
                           irhswell=reg%irhswell(j), irhswcomp=reg%irhswcomp(j), is_sulfate=is_sulfate)
    if(rc /=0) then
     call final_(rc)
     return
    endif
   end do

!  Establish the elements
!  ----------------------
!  NOTE: Hard coded optionals and parameters here
   do i = 1, reg%NELEM
    call CARMAELEMENT_Create(r, i, reg%igroup(i), &
                             reg%elemname(i), reg%rhop(i), reg%itype(i), reg%icomposition(i), rc)
    if(rc /=0) then
     call final_(rc)
     return
    endif
   end do

!  Establish the gases
!  -------------------
   do i = 1, reg%NGAS
    select case (reg%igcomp(i))
     case (1) 
      gwtmol = WTMOL_H2O
     case (2)
      gwtmol = WTMOL_H2SO4
     case (3)
      gwtmol = WTMOL_SO2
     case (4)
      gwtmol = WTMOL_HNO3
     case default
      print *, 'Unknown gas IGCOMP from CARMAchem_Registry.rc for gas ',i
      call final_(-100)
      return
    end select
    call CARMAGAS_Create(r, i, reg%gasname(i), &
                         gwtmol, reg%igvapreq(i), reg%igcomp(i), rc, ds_threshold=-0.2_f)
    if(rc /=0) then
     call final_(rc)
     return
    endif
   end do

!  Check the group/element/gas names and assign some indices
!  ---------------------------------------------------------
!  Look for pure groups, mixed group, and sulfuric acid gas
   do i = 1, reg%NELEM
    j = reg%igroup(i)
    groupname = ESMF_UtilStringUpperCase(trim(reg%groupname(j)))
    elemname  = ESMF_UtilStringUpperCase(trim(reg%elemname(i)))
    if(groupname == 'SULFATE') then
     reg%igrp_sulfate = j
     if(elemname == 'PC') reg%ielm_sulfate = i
    endif
    if(groupname == 'DUST') then
     reg%igrp_dust = j
     if(elemname == 'PC') reg%ielm_dust = i
    endif
    if(groupname == 'SEASALT') then
     reg%igrp_seasalt = j
     if(elemname == 'PC') reg%ielm_seasalt = i
    endif
    if(groupname == 'SMOKE') then
     reg%igrp_smoke = j
     if(elemname == 'PC') reg%ielm_smoke = i
    endif
    if(groupname == 'BLACK_CARBON') then
     reg%igrp_black_carbon = j
     if(elemname == 'PC') reg%ielm_black_carbon = i
    endif
    if(groupname == 'ASH') then
     reg%igrp_ash = j
     if(elemname == 'PC') reg%ielm_ash = i
    endif
!   Mixed group may contain sulfate element
    if(groupname == 'MIXEDP') then
     reg%igrp_mixed = j
     if(elemname == 'PC')           reg%ielm_mxpc      = i
     if(elemname == 'SULFATE')      reg%ielm_mxsulfate = i
     if(elemname == 'DUST')         reg%ielm_mxdust    = i
     if(elemname == 'SEASALT')      reg%ielm_mxseasalt = i
     if(elemname == 'SMOKE')        reg%ielm_mxsmoke   = i
     if(elemname == 'ASH')          reg%ielm_mxash     = i
     if(elemname == 'BLACK_CARBON') reg%ielm_mxbc      = i
    endif
   end do

   do i = 1, reg%NGAS
    gasname = ESMF_UtilStringUpperCase(trim(reg%gasname(i)))
    if(gasname == 'H2SO4') reg%igas_h2so4 = i
    if(gasname == 'H2O'  ) reg%igas_h2o   = i
    if(gasname == 'HNO3'  ) reg%igas_hno3 = i
   end do

!  NEED:
!  Hooks to CARMA_Solute


!  Setup Growth/Nucleation
!  -----------------------
!  Check that growth is correctly implemented based on elements/gases
   if(reg%do_grow) then
    if(reg%igrp_sulfate < 0 .or. reg%ielm_sulfate < 0 .or. &
       reg%igas_h2so4 < 0 .or. reg%igas_h2o < 0 ) then
       reg%do_grow = .false.
       print *, 'Not set up correctly for growth; do_grow set false'
    endif

   endif

   if(reg%do_grow) then
!   Pure sulfate group
    call CARMA_AddGrowth(r, reg%ielm_sulfate, reg%igas_h2so4, rc)
    if(rc /=0) then
     call final_(rc)
     return
    endif

!   Add growth to the sulfate element of the mixed group (which is nominally
!   the particle concentration element)
    if(reg%ielm_mxsulfate > 0) then
     call CARMA_AddGrowth(r, reg%ielm_mxsulfate, reg%igas_h2so4, rc)
     if(rc /=0) then
      call final_(rc)
      return
     endif
    endif

!   Add nucleation
    call CARMA_AddNucleation(r, reg%ielm_sulfate, reg%ielm_sulfate, &
                             I_HOMNUC, 0._f, rc, igas=reg%igas_h2so4)
    if(reg%igrp_mixed > 0) then
     if(reg%ielm_dust > 0 .and. reg%ielm_mxdust > 0) &
      call CARMA_AddNucleation(r, reg%ielm_dust, reg%ielm_mxdust, &
                               I_HETNUCSULF, 0._f, rc, igas=reg%igas_h2so4, &
                               ievp2elem=reg%ielm_dust)
     if(reg%ielm_smoke > 0 .and. reg%ielm_mxsmoke > 0) &
      call CARMA_AddNucleation(r, reg%ielm_smoke, reg%ielm_mxsmoke, &
                               I_HETNUCSULF, 0._f, rc, igas=reg%igas_h2so4, &
                               ievp2elem=reg%ielm_smoke)
     if(reg%ielm_seasalt > 0 .and. reg%ielm_mxseasalt > 0) &
      call CARMA_AddNucleation(r, reg%ielm_seasalt, reg%ielm_mxseasalt, &
                               I_HETNUCSULF, 0._f, rc, igas=reg%igas_h2so4, &
                               ievp2elem=reg%ielm_seasalt)
    endif
    if(rc /=0) then
     call final_(rc)
     return
    endif
   endif


!  Setup Coagulation
!  --------------------
!  We set up self coagulation for pure SULFATE and SMOKE groups.
!  If there is MIXEDP we allow coagulation of SULFATE with MIXEDP
   if(reg%do_coag) then
    do i = 1, reg%NELEM
     j = reg%igroup(i)
     groupname = ESMF_UtilStringUpperCase(trim(reg%groupname(j)))
     elemname  = ESMF_UtilStringUpperCase(trim(reg%elemname(i)))
     if(groupname == 'MIXEDP' .and. elemname == 'SULFATE') reg%mixedcorecomp = reg%icomposition(i)

!    This block adds the self coagulation of the pure sulfate group and the
!    coagulation of the pure sulfate to the mixed group.
     if( groupname == 'SULFATE' ) then
      call CARMA_AddCoagulation(r, j, j, j, I_COLLEC_FUCHS, rc )
      if(reg%igrp_mixed > 0 .AND. reg%icomposition(i) == reg%mixedcorecomp) then
       call CARMA_AddCoagulation(r, j, reg%igrp_mixed, reg%igrp_mixed, I_COLLEC_FUCHS, rc )
      endif
     endif
     if( groupname == 'SMOKE' ) call CARMA_AddCoagulation(r, j, j, j, I_COLLEC_FUCHS, rc )
     if(rc /=0) then
      call final_(rc)
      return
     endif
    end do
   endif



!  Initialize CARMA
!  ----------------
   call CARMA_Initialize(r, rc, &
                         do_cnst_rlh=reg%do_cnst_rlh, do_coag=reg%do_coag, &
                         do_detrain=reg%do_detrain, do_fixedinit=reg%do_fixedinit, &
                         do_grow=reg%do_grow, do_incloud=reg%do_incloud, &
                         do_explised=reg%do_explised, do_print_init=reg%do_print_init, &
                         do_substep=reg%do_substep, &
                         do_thermo=reg%do_thermo, do_vdiff=reg%do_vdiff, &
                         do_vtran=reg%do_vtran, vf_const=reg%vf_const, conmax=reg%conmax, &
                         minsubsteps=reg%minsubsteps, maxsubsteps=reg%maxsubsteps, &
                         maxretries=reg%maxretries, dt_threshold=1._f )

!  Get the dust emissions size fraction
!  -----------------------
!   Look for dust aerosol group / element
    do i = 1, reg%NELEM
     j = reg%igroup(i)
     groupname  = ESMF_UtilStringUpperCase(trim(reg%groupname(j)))
     if(groupname == 'DUST' .OR. ESMF_UtilStringUpperCase(trim(reg%elemname(i))) == 'DUST') then
      allocate(radius_(reg%NBIN), rlow_(reg%NBIN), rup_(reg%NBIN), __STAT__)
      allocate(radius(reg%NBIN),  rlow(reg%NBIN),  rup(reg%NBIN), __STAT__)
      allocate(rhod_(reg%NBIN),  rhog_(reg%NBIN),  __STAT__)
      call CARMAGroup_Get(r, j, rc, r=radius_, rlow=rlow_, rup=rup_)
      radius = radius_ / 100.  ! go from CARMA cm -> m
      rlow   = rlow_   / 100.
      rup    = rup_    / 100.
      call CARMAElement_Get(r, i, rc, rho=rhod_)
      call CARMAElement_Get(r, r%f_group(j)%f_ienconc, rc, rho=rhog_)
      rhod   = rhod_(1) * 1000.   ! go from CARMA to MKS
      rhog   = rhog_(1) * 1000.   ! go from CARMA to MKS
      call KokSizeDistribution(radius, rlow, rup, reg%dmass_dust, rhod=rhod, rhog=rhog)
      deallocate(radius, rlow, rup, radius_, rlow_, rup_, rhod_, rhog_, __STAT__)
     endif
    enddo



!  Print information
!  -----------------
   IF( MAPL_AM_I_ROOT() ) THEN
    call dumpGroup(r, rc)
    if(rc /=0) then
     call final_(104)
     return
    endif
    call dumpElement(r, rc)
    if(rc /=0) then
     call final_(105)
     return
    endif
   END IF

!! Housekeeping
!! ------------
!   deallocate ( r, stat=ios )
!   if ( ios /= 0) then
!    call final_(200)
!    return
!   endif
!   ier(:)=0

  RETURN

CONTAINS

   SUBROUTINE init_()
   INTEGER :: ios, n
   n=128
   ios=0
   ALLOCATE ( ier(n), stat=ios )
   IF ( ios /= 0 ) rc = 100
   END SUBROUTINE init_

   SUBROUTINE final_(ierr)
   INTEGER :: ios, ierr
   DEALLOCATE ( r, ier, stat=ios )
   CALL I90_release()
   rc = ierr
   END SUBROUTINE final_
   
!-------------------------------------------------------------------------
!NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1, GEOS/DAS!
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  TFQuery: Find whether a token in T or F.
!
! !INTERFACE:
!
   LOGICAL FUNCTION TFQuery(string,fn)

! !USES:

    IMPLICIT NONE

! !INPUT PARAMETERS:

    CHARACTER(LEN=*), INTENT(IN) :: string, fn

! !OUTPUT PARAMETERS:

! !DESCRIPTION: Return the value (T or F) of a particular token (string)
!               in the file fn.

! !REVISION HISTORY:
!
!  15 Aug 2007  Nielsen    First version.
!EOP
!-------------------------------------------------------------------------
    CHARACTER(LEN=8) :: tOrF

    INTEGER :: rc
    rc = 0
    tOrF = ' '
    TFQuery = .FALSE.

    CALL I90_label ( TRIM(string), rc )
    IF(rc .NE. 0) THEN
     PRINT *,'Could not find ',TRIM(string),' in ',TRIM(fn)
     CALL final_(99)
    END IF

    CALL I90_Gtoken( tOrF, rc )
    IF(TRIM(tOrF) ==    'T' .OR. TRIM(tOrF) ==    't' .OR. &
       TRIM(tOrF) == 'TRUE' .OR. TRIM(tOrF) == 'true') TFQuery = .TRUE.

   END FUNCTION TFQuery

!-------------------------------------------------------------------------
!NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1, GEOS/DAS!
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  setmcor: Find area of tiles.
!
! !INTERFACE:
!
  SUBROUTINE setmcor(i1,i2,j1,j2,im,jm,lats,cellArea)

! !USES:

  IMPLICIT NONE

! !INPUT PARAMETERS:

  INTEGER, INTENT(IN) :: i1,i2,j1,j2,im,jm
  REAL, INTENT(IN) :: lats(i1:i2,j1:j2)   !radians

! !OUTPUT PARAMETERS:

  REAL, INTENT(OUT) :: cellArea(i1:i2,j1:j2)   !m^2

! !DESCRIPTION: Find the horizontal surface area (m^2) of each cell.  
!               In testing with 8-byte words, the total surface area 
!               was 4.0000508 PI R^2.

! !REVISION HISTORY:

!  15 Aug 2007  Nielsen    First version.
!EOP
!-------------------------------------------------------------------------
  REAL, PARAMETER :: ae=6.371E+06

  INTEGER :: i,j
  REAL :: scale,dlat,f,arg,err,pi

  err=1.00E-05
  pi=4.00*ATAN(1.00)
  dlat=pi/FLOAT(jm-1)
  scale = 2.00*pi*ae*ae/FLOAT(im)

  DO j=j1,j2
   DO i=i1,i2

! South pole

    IF( lats(i,j) < -0.50*pi+err ) THEN
     f=0.25
     arg=0.50*(lats(i,j)+lats(i,j+1))

! North pole

    ELSE IF( lats(i,j) > 0.50*pi-err ) THEN
     f=0.25
     arg=0.50*(lats(i,j-1)+lats(i,j))

! Interior

    ELSE
     f=1.00
     arg=lats(i,j)
    END IF

    cellArea(i,j)=scale*dlat*f*cos(arg)

   END DO
  END DO

  RETURN 
  END SUBROUTINE setmcor

   END SUBROUTINE CARMA_GridCompInitialize

!-------------------------------------------------------------------------
!NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1, GEOS/DAS!
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  CARMA_GridCompRun --- The CARMA Driver 
!
! !INTERFACE:
!

   SUBROUTINE CARMA_GridCompRun ( gcCARMA, qa, impChem, expChem, nymd, nhms, &
                                  cdt, rc )

   IMPLICIT none

! !INPUT/OUTPUT PARAMETERS:

   TYPE(CARMA_GridComp), INTENT(INOUT) :: gcCARMA ! Grid Component
   TYPE(Chem_Array), pointer           :: qa(:)   ! tracer array will go here

! !INPUT PARAMETERS:

   TYPE(ESMF_State), INTENT(INOUT) :: impChem ! Import State
   INTEGER, INTENT(IN) :: nymd, nhms	      ! time
   REAL,    INTENT(IN) :: cdt		      ! chemical timestep (secs)

! !OUTPUT PARAMETERS:

   TYPE(ESMF_State), INTENT(INOUT) :: expChem   ! Export State
   INTEGER, INTENT(OUT) ::  rc                  ! Error return code:
                                                !  0 - all is well
                                                !  1 -

! !DESCRIPTION: This routine implements the CARMA driver
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
!  12Aug2009 Colarco   First crack at CARMA run method
!
!EOP
!-------------------------------------------------------------------------

   CHARACTER(LEN=*), PARAMETER :: IAm = 'CARMA_GridCompRun'
   INTEGER                     :: STATUS

!  Input fields from GEOS-5
!  ------------------------
   REAL, POINTER, DIMENSION(:,:,:) :: p, ple, rhoa, tmpu, zc, zl, q, zle, &
                                      rh, u, v, su_nuc, &
                                      zsubsteps, sath2so4, su_sareav, &
                                      su_sarea, su_numd, su_reff, &
                                      du_sarea, du_numd, du_reff, &
                                      ash_sarea, ash_numd, ash_reff, &
                                      ss_sarea, ss_numd, ss_reff, &
                                      sm_sarea, sm_numd, sm_reff, &
                                      mx_sarea, mx_numd, mx_reff, &
                                      su_mass, hno3, h2so4
   REAL, POINTER, DIMENSION(:,:)   :: gwettop, fraclake, oro, u10m, v10m, &
                                      ustar, pblh, z0h, shflux, precc, precl, &
                                      substeps, retries
   real, pointer, dimension(:,:)   :: du_sed, su_sed, ss_sed, bc_sed, ash_sed, sm_sed, &
                                      mxdu_sed, mxsu_sed, mxss_sed, mxbc_sed, mxash_sed, mxsm_sed
   type(Chem_Array), pointer       :: suvf(:), mxvf(:)


!  Local
!  -----
   INTEGER :: i, i1, i2, ic, ier(512), im, ijl, ios
   INTEGER :: j, j1, j2, jm
   INTEGER :: ielem, ibin, igrp, igas
   INTEGER :: k, km, kReverse
   INTEGER :: n, n2, nbegin, nend, nCARMABegin, nCARMAEnd
   INTEGER :: nymd1, nhms1
   INTEGER :: substep_int, last_sub
   real(kind=f) :: retry_real, last_ret
   logical :: rootproc
   real(kind=f) :: dtime
   character(len=ESMF_MAXSTR)        :: binstr

   INTEGER, PARAMETER :: ToCARMA = 1
   INTEGER, PARAMETER :: FromCARMA = -1
!  We are using the CARMA constants here (CGS units) but need
!  MKS values to go back to GEOS-5
   REAL, PARAMETER    :: grav_mks = grav/100.
   integer :: igroup
   CHARACTER(LEN=255) :: groupname, elemname, gasname

   REAL :: qmax,qmin
   real(kind=f) :: lon, lat
   real(kind=f), allocatable :: xc(:), dx(:), yc(:), dy(:)
   real(kind=f), allocatable :: p_(:), ple_(:), tmpu_(:), zc_(:), zl_(:), &
                                q_(:), rh_(:), nuc_(:), sarea_(:), numd_(:), &
                                r_wet_(:), reff_num(:), reff_den(:), vf_(:), &
                                zsubsteps_(:)
   real(kind=f), allocatable :: satice_(:,:), satliq_(:,:), told_(:), gasold_(:,:)
   real(kind=f) :: dq_

   type(carmastate_type) :: cstate
   type(carma_type), pointer :: r => null()
   type(CARMA_Registry), pointer :: reg => null()
   character(len=255) :: string

!  For a reference atmosphere we'll choose some values
   real, parameter, dimension(73) :: pleRef = & 
         (/ 1, 2, 3, 4, 6, 8, 11, 15, 21, 27, 36, 47, 61, 79, 101, 130,       &
            165, 208, 262, 327, 407, 504, 621, 761, 929, 1127, 1364, 1645,  &
            1979, 2373, 2836, 3381, 4017, 4764, 5638, 6660, 7851, 9236,     &
            10866, 12783, 15039, 17693, 20792, 24398, 28606, 33388, 37003,  &
            40612, 44214, 47816, 51405, 54997, 58584, 62170, 65769, 68147,  &
            70540, 72931, 75313, 77711, 79623, 81046, 82485, 83906, 85344,  &
            86765, 88201, 89636, 91071, 92516, 93921, 95376, 100000 /) 
   real, parameter, dimension(72) :: tmpuRef = & 
         (/ 219, 221, 223, 228, 230, 230, 232, 238, 245, 253, 259, 263, &
            264, 262, 258, 253, 247, 239, 233, 229, 227, 227, 226, 223, &
            222, 221, 220, 219, 218, 217, 216, 215, 214, 213, 212, 212, &
            214, 214, 216, 219, 219, 210, 210, 218, 227, 234, 240, 245, &
            250, 254, 257, 260, 262, 263, 265, 266, 267, 268, 269, 270, &
            270, 270, 270, 270, 271, 271, 271, 270, 267, 265, 266, 266 /)
   real, parameter, dimension(72) :: rhRef = 1e-6 * &
         (/ 1, 2, 2, 2, 3, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 6, 18, 51,          &
            129, 267, 394, 502, 682, 1135, 1603, 2076, 2820, 3792, 5120,     &
            6806, 8912, 11597, 15397, 20386, 28168, 29755, 28748, 33875,     &
            34058, 28657, 43458, 401856, 947266, 932618, 902344, 657227,     &
            371583, 203370, 235108, 317872, 413086, 511719, 691407, 686524,  &
            601563, 456055, 475098, 626954, 590821, 483399, 380860, 297852,  &
            230958, 183594, 144288, 111084, 96558, 136963, 369629, 770508,   &
            793946, 799805 /) 

   ier(:) = 0

!  Short-hand to object
   r   => gcCARMA%carma
   reg => gcCARMA%CARMAreg

!  Grid specs from Chem_Bundle%grid
!  --------------------------------
   rc = 0
   i1 = gcCARMA%i1
   i2 = gcCARMA%i2
   im = gcCARMA%im
   
   j1 = gcCARMA%j1
   j2 = gcCARMA%j2
   jm = gcCARMA%jm
   
   km = gcCARMA%km
   
   ijl = (i2-i1+1)*(j2-j1+1)

   dtime = cdt

!  Location of species from Chem_Bundle%registry.
!  ----------------------------------------------
   nCARMABegin =  1
   nCARMAEnd   =  gcCARMA%CARMAreg%nq

   rootProc=.FALSE.
   IF( MAPL_AM_I_ROOT() ) THEN
    rootProc=.TRUE.
   END IF

!  Allocate
!  --------
   allocate(p(i1:i2,j1:j2,km), __STAT__ )
   allocate(xc(km), dx(km), yc(km), dy(km), &
            p_(km), ple_(km+1), tmpu_(km), zc_(km), zl_(km+1), &
            q_(km), rh_(km), nuc_(km), sarea_(km), numd_(km), &
            r_wet_(km), reff_num(km), reff_den(km), vf_(km+1), &
            zsubsteps_(km), __STAT__ )
   allocate(told_(km), gasold_(km,reg%NGAS), satice_(km,reg%NGAS), satliq_(km,reg%NGAS), __STAT__ )

!  Get Imports
!  -----------
   call MAPL_GetPointer ( impChem, rhoa, 'AIRDENS', __RC__)
   call MAPL_GetPointer ( impChem, ple, 'PLE', __RC__)
   call MAPL_GetPointer ( impChem, zle, 'ZLE', __RC__)
   call MAPL_GetPointer ( impChem, q, 'Q', __RC__)
   call MAPL_GetPointer ( impChem, rh, 'RH2', __RC__)
   call MAPL_GetPointer ( impChem, tmpu, 'T', __RC__)
   call MAPL_GetPointer ( impChem, ustar, 'USTAR', __RC__)
   call MAPL_GetPointer ( impChem, fraclake, 'FRLAKE', __RC__)
   call MAPL_GetPointer ( impChem, gwettop, 'WET1', __RC__)
   call MAPL_GetPointer ( impChem, u10m, 'U10M', __RC__)
   call MAPL_GetPointer ( impChem, v10m, 'V10M', __RC__)
   call MAPL_GetPointer ( impChem, pblh, 'ZPBL', __RC__)
   call MAPL_GetPointer ( impChem, z0h, 'Z0H', __RC__)
   call MAPL_GetPointer ( impChem, shflux, 'SH', __RC__)
   call MAPL_GetPointer ( impChem, precl, 'NCN_PRCP', __RC__)
   call MAPL_GetPointer ( impChem, precc, 'CN_PRCP', __RC__)
   call MAPL_GetPointer ( impChem, oro, 'LWI', __RC__)
   call MAPL_GetPointer ( impChem, u, 'U', __RC__)
   call MAPL_GetPointer ( impChem, v, 'V', __RC__)
   call MAPL_GetPointer ( impChem, hno3,  'CARMA_HNO3',    notFoundOK=.TRUE., __RC__)
   call MAPL_GetPointer ( impChem, h2so4, 'CARMA_H2SO4',   notFoundOK=.TRUE., __RC__)

!  Fill the internal state with direct gas species from GMI
!  Expectation is species are in VMR and needed in MMR for CARMA
!  -----------
   if(reg%gmi_chem_provider .and. reg%NGAS > 0) then
    do igas = 1, reg%NGAS
     n  = nCARMAbegin + reg%NBIN*reg%NELEM - 1 + igas
     gasname = ESMF_UtilStringUpperCase(reg%gasname(igas))
     if(gasname == 'H2SO4') qa(n)%data3d = h2so4*WTMOL_H2SO4/WTMOL_AIR
     if(gasname == 'HNO3' ) qa(n)%data3d = hno3 *WTMOL_HNO3 /WTMOL_AIR
    enddo
   endif

!  Get Exports
!  -----------
!  Mixed Group
   call MAPL_GetPointer(expChem, mx_sarea,  'CARMA_MXSAREA',  __RC__)
   call MAPL_GetPointer(expChem, mx_numd,   'CARMA_MXNUMD',   __RC__)
   call MAPL_GetPointer(expChem, mx_reff,   'CARMA_MXREFF',   __RC__)
!  Dust
   call MAPL_GetPointer(expChem, du_sed,    'CARMA_DUSD',   __RC__)
   call MAPL_GetPointer(expChem, du_sarea,  'CARMA_DUSAREA',  __RC__)
   call MAPL_GetPointer(expChem, du_numd,   'CARMA_DUNUMD',   __RC__)
   call MAPL_GetPointer(expChem, du_reff,   'CARMA_DUREFF',   __RC__)
   call MAPL_GetPointer(expChem, mxdu_sed,  'CARMA_MXDUSD',   __RC__)
!  Ash
   call MAPL_GetPointer(expChem, ash_sed,    'CARMA_ASHSD',   __RC__)
   call MAPL_GetPointer(expChem, ash_sarea,  'CARMA_ASHSAREA',  __RC__)
   call MAPL_GetPointer(expChem, ash_numd,   'CARMA_ASHNUMD',   __RC__)
   call MAPL_GetPointer(expChem, ash_reff,   'CARMA_ASHREFF',   __RC__)
   call MAPL_GetPointer(expChem, mxash_sed,  'CARMA_MXASHSD',   __RC__)
!  Sulfate
   call MAPL_GetPointer(expChem, su_sed,    'CARMA_SUSD',   __RC__)
   call MAPL_GetPointer(expChem, su_nuc,    'CARMA_SUNUC',  __RC__)
   call MAPL_GetPointer(expChem, su_sarea,  'CARMA_SUSAREA',  __RC__)
   call MAPL_GetPointer(expChem, su_mass,   'CARMA_SUMASS',  __RC__)
   call MAPL_GetPointer(expChem, su_numd,   'CARMA_SUNUMD',   __RC__)
   call MAPL_GetPointer(expChem, su_reff,   'CARMA_SUREFF',   __RC__)
   call MAPL_GetPointer(expChem, mxsu_sed,  'CARMA_MXSUSD',   __RC__)
!  Sea salt
   call MAPL_GetPointer(expChem, ss_sed,    'CARMA_SSSD',   __RC__)
   call MAPL_GetPointer(expChem, ss_sarea,  'CARMA_SSSAREA',  __RC__)
   call MAPL_GetPointer(expChem, ss_numd,   'CARMA_SSNUMD',   __RC__)
   call MAPL_GetPointer(expChem, ss_reff,   'CARMA_SSREFF',   __RC__)
   call MAPL_GetPointer(expChem, mxss_sed,  'CARMA_MXSSSD',   __RC__)
!  Smoke
   call MAPL_GetPointer(expChem, sm_sed,    'CARMA_SMSD',   __RC__)
   call MAPL_GetPointer(expChem, sm_sarea,  'CARMA_SMSAREA',  __RC__)
   call MAPL_GetPointer(expChem, sm_numd,   'CARMA_SMNUMD',   __RC__)
   call MAPL_GetPointer(expChem, sm_reff,   'CARMA_SMREFF',   __RC__)
   call MAPL_GetPointer(expChem, mxsm_sed,  'CARMA_MXSMSD',   __RC__)
!  Other
   call MAPL_GetPointer(expChem, bc_sed,    'CARMA_BCSD',   __RC__)
   call MAPL_GetPointer(expChem, mxbc_sed,  'CARMA_MXBCSD',   __RC__)
   call MAPL_GetPointer(expChem, substeps,  'CARMA_SUBSTEPS',   __RC__)
   call MAPL_GetPointer(expChem, retries,   'CARMA_RETRIES',    __RC__)
   call MAPL_GetPointer(expChem, zsubsteps, 'CARMA_ZSUBSTEPS',  __RC__)
   call MAPL_GetPointer(expChem, sath2so4,  'CARMA_SATH2SO4',   __RC__)
   call MAPL_GetPointer(expChem, su_sareav,  'CARMA_SUSAREAv',  __RC__)


!  Allocate space for fall velocity diagnostic and see if requested
!  ----------------------------------------------------------------
   allocate(suvf(reg%NBIN), mxvf(reg%NBIN), __STAT__)
   do ibin = 1, reg%NBIN
    write(binstr,'(i3)') ibin
    binstr = adjustl(binstr)
    if(ibin .lt. 10)  binstr = '0'//binstr
    if(ibin .lt. 100) binstr = '0'//binstr
!    if(MAPL_AM_I_ROOT()) print *, 'CARMA_SUVF'//trim(binstr), 'CARMA_MXVF'//trim(binstr)
!    call MAPL_GetPointer(expChem, suvf(ibin)%data3d, 'CARMA_SUVF'//trim(binstr), __RC__)
!    call MAPL_GetPointer(expChem, mxvf(ibin)%data3d, 'CARMA_MXVF'//trim(binstr), __RC__)
   enddo

!  Get the mid-point pressure
!  --------------------------
   DO k=1,km
    p(i1:i2,j1:j2,k)=exp((log(ple(i1:i2,j1:j2,k-1))+log(ple(i1:i2,j1:j2,k)) )*0.50)
   END DO

!   call pmaxmin('CARMA::U:       ', u(i1:i2,j1:j2,1:km), qmin, qmax, ijl, km, 1. )
!   call pmaxmin('CARMA::V:       ', v(i1:i2,j1:j2,1:km), qmin, qmax, ijl, km, 1. )
!   call pmaxmin('CARMA::AIRDENS: ', rhoa(i1:i2,j1:j2,1:km), qmin, qmax, ijl, km, 1. )
!   call pmaxmin('CARMA::Q:       ', q(i1:i2,j1:j2,1:km), qmin, qmax, ijl, km, 1. )
!   call pmaxmin('CARMA::RH:      ', rh(i1:i2,j1:j2,1:km), qmin, qmax, ijl, km, 1. )
!   call pmaxmin('CARMA::P:       ', p(i1:i2,j1:j2,1:km), qmin, qmax, ijl, km, 1. )
!   call pmaxmin('CARMA::ZLE:     ', zle(i1:i2,j1:j2,1:km), qmin, qmax, ijl, km, 1. )
!   call pmaxmin('CARMA::T:       ', tmpu(i1:i2,j1:j2,1:km), qmin, qmax, ijl, km, 1. )
!   call pmaxmin('CARMA::ORO:     ', oro(i1:i2,j1:j2), qmin, qmax, ijl, 1, 1. )
!   call pmaxmin('CARMA::USTAR:   ', ustar(i1:i2,j1:j2), qmin, qmax, ijl, 1, 1. )
!   call pmaxmin('CARMA::FRACLAKE:', fraclake(i1:i2,j1:j2), qmin, qmax, ijl, 1, 1. )
!   call pmaxmin('CARMA::U10M:    ', u10m(i1:i2,j1:j2), qmin, qmax, ijl, 1, 1. )
!   call pmaxmin('CARMA::V10M:    ', v10m(i1:i2,j1:j2), qmin, qmax, ijl, 1, 1. )
!   call pmaxmin('CARMA::PBLH:    ', pblh(i1:i2,j1:j2), qmin, qmax, ijl, 1, 1. )
!   call pmaxmin('CARMA::Z0H:     ', z0h(i1:i2,j1:j2), qmin, qmax, ijl, 1, 1. )
!   call pmaxmin('CARMA::SHFLUX:  ', shflux(i1:i2,j1:j2), qmin, qmax, ijl, 1, 1. )
!   call pmaxmin('CARMA::PRECL:   ', precl(i1:i2,j1:j2), qmin, qmax, ijl, 1, 1. )
!   call pmaxmin('CARMA::PRECC:   ', precc(i1:i2,j1:j2), qmin, qmax, ijl, 1, 1. )
!   call pmaxmin('CARMA::GWETTOP: ', gwettop(i1:i2,j1:j2), qmin, qmax, ijl, 1, 1. )

   IF( ANY(ier(:) /= 0) ) THEN
    PRINT *,Iam,': Failed on MAPL_GetPointer for imports in CARMA_GridCompRun.'
    rc = 11 
    RETURN
   END IF
   ier(:)=0

! For substepping you want to remember the old temperature.
! This is set in the internal_spec but possibly is bootstrapped.
! If bootstrapped set to current temperature.
  n = nCARMAbegin + reg%NBIN*reg%NELEM + reg%NGAS
  if(qa(n)%data3d(i1,j2,km) < 0.) qa(n)%data3d = tmpu
! And same for gases -- first, initialize water vapor to current
  do igas = 1, reg%NGAS
   n  = nCARMAbegin + reg%NBIN*reg%NELEM - 1 + igas
   if(trim(reg%gasname(igas)) == 'h2o' .or. trim(reg%gasname(igas)) == 'H2O') qa(n)%data3d = q
   n2 = nCARMAbegin + reg%NBIN*reg%NELEM + reg%NGAS + igas
   if(qa(n2)%data3d(i1,j2,km) < 0.) qa(n2)%data3d = qa(n)%data3d
  enddo

#ifdef DEBUG
if(reg%NGAS > 0) then
   n = reg%NBIN*reg%NELEM + 1
   call pmaxmin('CARMA::h2o_0:           ', qa(n)%data3d(i1:i2,j1:j2,1:km), qmin, qmax, ijl, km, 1. )
   n = reg%NBIN*reg%NELEM + reg%NGAS
   call pmaxmin('CARMA::h2so4_0:         ', qa(n)%data3d(i1:i2,j1:j2,1:km), qmin, qmax, ijl, km, 1. )
   call pmaxmin('CARMA::su001_0:         ', qa(1)%data3d(i1:i2,j1:j2,1:km), qmin, qmax, ijl, km, 1. )
   n = reg%NBIN*reg%NELEM + reg%NGAS + 1
   call pmaxmin('CARMA::told_0:          ', qa(n)%data3d(i1:i2,j1:j2,1:km), qmin, qmax, ijl, km, 1. )
   n = reg%NBIN*reg%NELEM + reg%NGAS + reg%NGAS
   call pmaxmin('CARMA::h2o_old_0:       ', qa(n)%data3d(i1:i2,j1:j2,1:km), qmin, qmax, ijl, km, 1. )
   n = reg%NBIN*reg%NELEM + reg%NGAS + reg%NGAS + 1
   call pmaxmin('CARMA::h2so4_old_0:     ', qa(n)%data3d(i1:i2,j1:j2,1:km), qmin, qmax, ijl, km, 1. )
   n = reg%NBIN*reg%NELEM + reg%NGAS + reg%NGAS + reg%NGAS + 1
   call pmaxmin('CARMA::satliq2_old_0:   ', qa(n)%data3d(i1:i2,j1:j2,1:km), qmin, qmax, ijl, km, 1. )
   n = reg%NBIN*reg%NELEM + reg%NGAS + reg%NGAS + reg%NGAS + reg%NGAS + 1
   call pmaxmin('CARMA::satice2_old_0:   ', qa(n)%data3d(i1:i2,j1:j2,1:km), qmin, qmax, ijl, km, 1. )
endif
#endif

!  ====================  CARMA Step ================================
!  Establish the CARMA state, do a step, and retain diagnostic
   if( associated(DU_sed))    DU_sed(:,:) = 0.
   if( associated(SU_sed))    SU_sed(:,:) = 0.
   if( associated(SS_sed))    SS_sed(:,:) = 0.
   if( associated(BC_sed))    BC_sed(:,:) = 0.
   if( associated(SM_sed))    SM_sed(:,:) = 0.
   if( associated(ASH_sed))   ASH_sed(:,:) = 0.
   if( associated(MXDU_sed))  MXDU_sed(:,:) = 0.
   if( associated(MXSU_sed))  MXSU_sed(:,:) = 0.
   if( associated(MXSM_sed))  MXSM_sed(:,:) = 0.
   if( associated(MXSS_sed))  MXSS_sed(:,:) = 0.
   if( associated(MXBC_sed))  MXBC_sed(:,:) = 0.
   if( associated(MXASH_sed)) MXASH_sed(:,:) = 0.
   if( associated(SU_sarea))  SU_sarea(:,:,:) = 0.
   if( associated(SU_numd))   SU_numd(:,:,:) = 0.
   if( associated(SU_reff))   SU_reff(:,:,:) = 0.
   if( associated(DU_sarea))  DU_sarea(:,:,:) = 0.
   if( associated(DU_numd))   DU_numd(:,:,:) = 0.
   if( associated(DU_reff))   DU_reff(:,:,:) = 0.
   if( associated(SS_sarea))  SS_sarea(:,:,:) = 0.
   if( associated(SS_numd))   SS_numd(:,:,:) = 0.
   if( associated(SS_reff))   SS_reff(:,:,:) = 0.
   if( associated(SM_sarea))  SM_sarea(:,:,:) = 0.
   if( associated(SM_numd))   SM_numd(:,:,:) = 0.
   if( associated(SM_reff))   SM_reff(:,:,:) = 0.
   if( associated(SU_nuc))    SU_nuc(:,:,:) = 0. 
   if( associated(substeps))  substeps(:,:)    = 0.
   if( associated(retries))   retries(:,:)     = 0.
   if( associated(zsubsteps)) zsubsteps(:,:,:) = 0.
   if( associated(sath2so4))  sath2so4(:,:,:)  = 0.
   if( associated(su_sareav)) su_sareav(:,:,:) = 0.
   if( associated(su_mass))   su_mass(:,:,:)  = 0.
   last_sub = 0
   last_ret = 0

!  Possibly create a CARMA reference state column for 1,1 column in tile
   if(reg%do_fixedinit) then
!    dx and dy are hack as if for a "b" resolution grid
     dx(:) = 2.5
     dy(:) = 2.
     xc(:) = 0.
     yc(:) = 0.
     p_(1:km) = exp((log(pleRef(1:km))+log(pleRef(2:km+1)) )*0.50)
     ple_(1:km+1) = pleRef
     zc_(:) = p_(:)/pleRef(km+1)
     zl_(1:km+1) = pleRef(1:km+1)/pleRef(km+1)
     tmpu_(:) = tmpuRef(:)
     rh_(:) = rhRef(:)
     call CARMASTATE_CreateFromReference(cstate, r, 1._f, dtime, km, &
                              I_HYBRID, I_CART, lat, lon,  &
                              xc, dx, yc, dy, &
                              zc_, zl_, p_, ple_, tmpu_, rc, &
                              relhum=rh_)
   endif

   do j = j1, j2
    do i = i1, i2
     
     lon = gcCARMA%lons(i,j)
     lat = gcCARMA%lats(i,j)
!    dx and dy are hack as if for a "b" resolution grid
     dx(:) = 2.5
     dy(:) = 2.
     xc(:) = lon
     yc(:) = lat
     p_(:) = p(i,j,:)
     ple_(1:km+1) = ple(i,j,0:km)
     zc_(:) = p(i,j,:)/ple(i,j,km)
     zl_(1:km+1) = ple(i,j,0:km)/ple(i,j,km)
     tmpu_(:) = tmpu(i,j,:)
     rh_(:) = rh(i,j,:)
!    prior time step values for sub-stepping
     n = reg%NBIN*reg%NELEM + reg%NGAS + 1
     told_(:) = qa(n)%data3d(i,j,:)
     do igas = 1, reg%NGAS
      gasold_(:,igas) = qa(n+igas)%data3d(i,j,:)
      satliq_(:,igas) = qa(n+reg%NGAS+igas)%data3d(i,j,:)
      satice_(:,igas) = qa(n+reg%NGAS+reg%NGAS+igas)%data3d(i,j,:)
     enddo

     call CARMASTATE_Create(cstate, r, 1._f, dtime, km, &
                            I_HYBRID, I_CART, lat, lon,  &
                            xc, dx, yc, dy, &
                            zc_, zl_, p_, ple_, tmpu_, rc, &
                            relhum=rh_, told=told_)

   ! Map the model MMR to CARMA
     do ielem = 1, reg%NELEM
      do ibin = 1, reg%NBIN
       n = nCARMAbegin + (ielem-1)*reg%NBIN + ibin - 1
       q_(:) = qa(n)%data3d(i,j,:)
       where(q_ < 1.e-32) q_ = 1.e-32
       call CARMASTATE_SetBin(cstate, ielem, ibin, q_, rc)
      end do
     end do

   ! Map the model gases to CARMA
     if(reg%NGAS > 0) then
     do igas = 1, reg%NGAS
      n = nCARMAbegin + reg%NELEM*reg%NBIN - 1 + igas
      q_(:) = qa(n)%data3d(i,j,:)
      where(q_ < 1.e-32) q_ = 1.e-32
!!     HACK: Want to put in 10 Tg S of H2SO4 per year total into band 30N - 30S
!!     all longitude between 20 - 25 km altitude.  So this is 2.e10 kg H2SO4
!!     per year over an area of 2.626e14 m2 over 32.65 hPa depth for levels 30 - 34.
!      if( (trim(reg%gasname(igas)) == 'h2so4' .or. trim(reg%gasname(igas)) == 'H2SO4') .and. &
!         lat >= -30. .and. lat <= 30.) then
!       do k = 30,34
!        q_(k) = q_(k) + 2.e10 / (365.*86400.)*cdt  / 2.626e14 / 3265. * 9.8
!       enddo
!      endif

      call CARMASTATE_SetGas(cstate, igas, q_(:), rc, &
                             mmr_old = gasold_(:,igas), satice_old=satice_(:,igas), &
                             satliq_old=satliq_(:,igas) )
     end do
     endif
			
   ! Execute the step
     call CARMASTATE_Step(cstate, ios)

   ! Map CARMA back to model MMR
     do ielem = 1, reg%NELEM
      do ibin = 1, reg%NBIN
       n = nCARMAbegin + (ielem-1)*reg%NBIN + ibin - 1
       call CARMASTATE_GetBin(cstate, ielem, ibin, &
                              q_, rc)
       where(q_ < 1.e-32) q_ = 1.e-32
       qa(n)%data3d(i,j,:) = q_(:)
      end do
     end do

   ! Map CARMA back to model gas
     if(reg%NGAS > 0) then
     do igas = 1, reg%NGAS
      n = nCARMAbegin + reg%NELEM*reg%NBIN - 1 + igas
      call CARMASTATE_GetGas(cstate, igas, q_, rc, &
                             satice=satice_(:,igas), satliq=satliq_(:,igas))
      where(q_ < 1.e-32) q_ = 1.e-32
      qa(n)%data3d(i,j,:) = q_(:)
!     Save current gas mixing ratio and saturations for "old" values of next step
      n = reg%NBIN*reg%NELEM + reg%NGAS + 1
      qa(n+igas)%data3d(i,j,:) = q_(:)
      qa(n+reg%NGAS+igas)%data3d(i,j,:) = satliq_(:,igas)
      qa(n+reg%NGAS+reg%NGAS+igas)%data3d(i,j,:) = satice_(:,igas)
!     Save h2so4 supersaturation if it's asked for
      if( ESMF_UtilStringUpperCase(trim(reg%gasname(igas))) == 'H2SO4' .and. &
          associated(sath2so4) )   sath2so4(i,j,:) = satliq_(:,igas)
     end do

!    Hack - for now we assume gas does not change temperature, save told
     n = nCARMAbegin + reg%NELEM*reg%NBIN + reg%NGAS
     qa(n)%data3d(i,j,:) = tmpu_
     endif

!    Get requested sedimentation flux diagnostics per element
     do ielem = 1, reg%NELEM
      igroup = gcCARMA%CARMAreg%igroup(ielem)
      groupname = ESMF_UtilStringUpperCase(trim(gcCARMA%CARMAreg%groupname(igroup)))
      elemname  = ESMF_UtilStringUpperCase(trim(gcCARMA%CARMAreg%elemname(ielem)))
      do ibin = 1, reg%NBIN
       n = nCARMAbegin + (ielem-1)*reg%NBIN + ibin - 1
       call CARMASTATE_GetBin(cstate, ielem, ibin, &
                              q_, rc, sedimentationflux=dq_)
       if(associated(DU_sed)  .and. igroup .eq. reg%igrp_dust)         DU_sed(i,j) = DU_sed(i,j) + dq_
       if(associated(SS_sed)  .and. igroup .eq. reg%igrp_seasalt)      SS_sed(i,j) = SS_sed(i,j) + dq_
       if(associated(SM_sed)  .and. igroup .eq. reg%igrp_smoke)        SM_sed(i,j) = SM_sed(i,j) + dq_
       if(associated(SU_sed)  .and. igroup .eq. reg%igrp_sulfate)      SU_sed(i,j) = SU_sed(i,j) + dq_
       if(associated(BC_sed)  .and. igroup .eq. reg%igrp_black_carbon) BC_sed(i,j) = BC_sed(i,j) + dq_
       if(associated(ASH_sed) .and. igroup .eq. reg%igrp_ash)          ASH_sed(i,j) = ASH_sed(i,j) + dq_
!      Mixed group -- assume "pc" element is sulfate and subtract cores
       if(igroup .eq. reg%igrp_mixed) then
        if(associated(MXDU_sed) .and. ielem .eq. reg%ielm_mxdust)      MXDU_sed(i,j) = MXDU_sed(i,j) + dq_
        if(associated(MXSS_sed) .and. ielem .eq. reg%ielm_mxseasalt)   MXSS_sed(i,j) = MXSS_sed(i,j) + dq_
        if(associated(MXSM_sed) .and. ielem .eq. reg%ielm_mxsmoke)     MXSM_sed(i,j) = MXSM_sed(i,j) + dq_
        if(associated(MXSU_sed) .and. ielem .eq. reg%ielm_mxsulfate)   MXSU_sed(i,j) = MXSU_sed(i,j) + dq_
        if(associated(MXBC_sed) .and. ielem .eq. reg%ielm_mxbc)        MXBC_sed(i,j) = MXBC_sed(i,j) + dq_
        if(associated(MXASH_sed) .and. ielem .eq. reg%ielm_mxash)      MXASH_sed(i,j) = MXASH_sed(i,j) + dq_
!       subtract cores
        if(associated(MXSU_sed) .and. ielem .eq. reg%ielm_mxdust)      MXSU_sed(i,j) = MXSU_sed(i,j) - dq_
        if(associated(MXSU_sed) .and. ielem .eq. reg%ielm_mxsmoke)     MXSU_sed(i,j) = MXSU_sed(i,j) - dq_
        if(associated(MXSU_sed) .and. ielem .eq. reg%ielm_mxseasalt)   MXSU_sed(i,j) = MXSU_sed(i,j) - dq_
        if(associated(MXSU_sed) .and. ielem .eq. reg%ielm_mxbc)        MXSU_sed(i,j) = MXSU_sed(i,j) - dq_
        if(associated(MXSU_sed) .and. ielem .eq. reg%ielm_mxash)       MXSU_sed(i,j) = MXSU_sed(i,j) - dq_
       endif
      end do
     end do

   ! If the fall velocity diagnostic is asked for, get it (note conversion
   ! of sign change to define positive)
     do ielem = 1, reg%NELEM
      igroup = gcCARMA%CARMAreg%igroup(ielem)
      groupname = ESMF_UtilStringUpperCase(trim(gcCARMA%CARMAreg%groupname(igroup)))
      elemname  = ESMF_UtilStringUpperCase(trim(gcCARMA%CARMAreg%elemname(ielem)))
      if(ielem /= r%f_group(igroup)%f_ienconc ) cycle
      do ibin = 1, reg%NBIN
       if(groupname == 'SULFATE' .and. associated(suvf(ibin)%data3d)) then
        call CARMASTATE_GetBin(cstate, ielem, ibin, &
                               q_, rc, vf=vf_)
        suvf(ibin)%data3d(i,j,:) = -1. * vf_
       endif
       if( (groupname == 'MIXEDP' .or. groupname == 'DUST') .and. &
           associated(mxvf(ibin)%data3d)) then
        call CARMASTATE_GetBin(cstate, ielem, ibin, &
                               q_, rc, vf=vf_)
        mxvf(ibin)%data3d(i,j,:) = -1. * vf_
       endif
      end do
     end do

!    Get the nucleation rate if it is asked for (m-3 s-1)
     do ielem = 1, reg%NELEM
      igroup = gcCARMA%CARMAreg%igroup(ielem)
      groupname = trim(gcCARMA%CARMAreg%groupname(igroup))
      if(groupname /= 'sulfate' .AND. groupname /= 'SULFATE') cycle
      if(ielem /= r%f_group(igroup)%f_ienconc ) cycle
      if(.not.associated(SU_nuc)) cycle
      do ibin = 1, reg%NBIN
       call CARMASTATE_GetBin(cstate, ielem, ibin, &
                              q_, rc, nucleationrate=nuc_)
       SU_nuc(i,j,:) = SU_nuc(i,j,:) + nuc_
      enddo
     enddo

!    Get the group effective wet radius (m), surface area, and number density
     do ielem = 1, reg%NELEM
      igroup = gcCARMA%CARMAreg%igroup(ielem)
      groupname = ESMF_UtilStringUpperCase(trim(gcCARMA%CARMAreg%groupname(igroup)))
      if(ielem /= r%f_group(igroup)%f_ienconc ) cycle
      reff_num = 0.
      reff_den = 0.
      do ibin = 1, reg%NBIN
        call CARMASTATE_GetBin(cstate, ielem, ibin, q_, rc, &
         r_wet=r_wet_, numberdensity=numd_, areadensity=sarea_)
        if(associated(MX_sarea) .and. igroup .eq. reg%igrp_mixed)   MX_sarea(i,j,:) = MX_sarea(i,j,:) + sarea_
        if(associated(MX_numd)  .and. igroup .eq. reg%igrp_mixed)   MX_numd(i,j,:)  = MX_numd(i,j,:)  + numd_
        if(associated(SU_sarea) .and. igroup .eq. reg%igrp_sulfate) SU_sarea(i,j,:) = SU_sarea(i,j,:) + sarea_
        if(associated(SU_numd)  .and. igroup .eq. reg%igrp_sulfate) SU_numd(i,j,:)  = SU_numd(i,j,:)  + numd_
        if(associated(SU_mass)  .and. igroup .eq. reg%igrp_sulfate) SU_mass(i,j,:)  = SU_mass(i,j,:)  + q_
        if(associated(DU_sarea) .and. igroup .eq. reg%igrp_dust)    DU_sarea(i,j,:) = DU_sarea(i,j,:) + sarea_
        if(associated(DU_numd)  .and. igroup .eq. reg%igrp_dust)    DU_numd(i,j,:)  = DU_numd(i,j,:)  + numd_
        if(associated(ASH_sarea) .and. igroup .eq. reg%igrp_ash)    ASH_sarea(i,j,:) = ASH_sarea(i,j,:) + sarea_
        if(associated(ASH_numd)  .and. igroup .eq. reg%igrp_ash)    ASH_numd(i,j,:)  = ASH_numd(i,j,:)  + numd_
        if(associated(SM_sarea) .and. igroup .eq. reg%igrp_smoke)   SM_sarea(i,j,:) = SM_sarea(i,j,:) + sarea_
        if(associated(SM_numd)  .and. igroup .eq. reg%igrp_smoke)   SM_numd(i,j,:)  = SM_numd(i,j,:)  + numd_
        if(associated(SS_sarea) .and. igroup .eq. reg%igrp_seasalt) SS_sarea(i,j,:) = SS_sarea(i,j,:) + sarea_
        if(associated(SS_numd)  .and. igroup .eq. reg%igrp_seasalt) SS_numd(i,j,:)  = SS_numd(i,j,:)  + numd_
        reff_num = reff_num + r_wet_**3.*numd_
        reff_den = reff_den + r_wet_**2.*numd_
      enddo
      if(associated(MX_reff) .and. igroup .eq. reg%igrp_mixed)      where(reff_den > 0) MX_reff(i,j,:) = reff_num / reff_den
      if(associated(SM_reff) .and. igroup .eq. reg%igrp_smoke)      where(reff_den > 0) SM_reff(i,j,:) = reff_num / reff_den
      if(associated(DU_reff) .and. igroup .eq. reg%igrp_dust)       where(reff_den > 0) DU_reff(i,j,:) = reff_num / reff_den
      if(associated(ASH_reff) .and. igroup .eq. reg%igrp_ash)       where(reff_den > 0) ASH_reff(i,j,:) = reff_num / reff_den
      if(associated(SU_reff) .and. igroup .eq. reg%igrp_sulfate)    where(reff_den > 0) SU_reff(i,j,:) = reff_num / reff_den
      if(associated(SS_reff) .and. igroup .eq. reg%igrp_seasalt)    where(reff_den > 0) SS_reff(i,j,:) = reff_num / reff_den
     enddo


!    Get the number of substeps, retries from CARMA state
     if(reg%do_grow) then
      call CARMASTATE_Get(cstate, rc, nsubstep=substep_int, &
                          nretry=retry_real, zsubsteps=zsubsteps_)
      if(associated(substeps))  substeps(i,j)    = REAL(substep_int-last_sub, kind=f)
      if(associated(retries))   retries(i,j)     = retry_real-last_ret
      if(associated(zsubsteps)) zsubsteps(i,j,:) = zsubsteps_
      last_sub = substep_int
      last_ret = retry_real
     endif

!  Hack -- for now don't change temperature
!   ! Get the updated temperature.
!     call CARMASTATE_GetState(cstate, rc, t=tmpu_)
!     tmpu(i,j,:) = tmpu_(:)

    end do
   end do

!  Return the updated gas species to GMI from the internal state
!  Expectation is species are in MMR and needed in VMR for GMI
!  -----------
   if(reg%gmi_chem_provider .and. reg%NGAS > 0) then
    do igas = 1, reg%NGAS
     n  = nCARMAbegin + reg%NBIN*reg%NELEM - 1 + igas
     gasname = ESMF_UtilStringUpperCase(reg%gasname(igas))
     if(gasname == 'H2SO4') h2so4 = qa(n)%data3d*WTMOL_AIR/WTMOL_H2SO4
!     Don't update HNO3
!     if(gasname == 'HNO3' ) hno3  = qa(n)%data3d*WTMOL_AIR /WTMOL_HNO3
    enddo
   endif



 ! Cleanup the carma state objects
   call CARMASTATE_Destroy(cstate, rc)


#ifdef DEBUG
if(reg%NGAS > 0) then
   n = reg%NBIN*reg%NELEM + reg%NGAS
   call pmaxmin('CARMA::h2o_1:           ', qa(n-1)%data3d(i1:i2,j1:j2,1:km), qmin, qmax, ijl, km, 1. )
   call pmaxmin('CARMA::h2so4_1:         ', qa(n)%data3d(i1:i2,j1:j2,1:km), qmin, qmax, ijl, km, 1. )
   call pmaxmin('CARMA::su001_1:         ', qa(1)%data3d(i1:i2,j1:j2,1:km), qmin, qmax, ijl, km, 1. )
   n = reg%NBIN*reg%NELEM + reg%NGAS + 1
   call pmaxmin('CARMA::told_1:          ', qa(n)%data3d(i1:i2,j1:j2,1:km), qmin, qmax, ijl, km, 1. )
   n = reg%NBIN*reg%NELEM + reg%NGAS + reg%NGAS
   call pmaxmin('CARMA::h2o_old_1:       ', qa(n)%data3d(i1:i2,j1:j2,1:km), qmin, qmax, ijl, km, 1. )
   n = reg%NBIN*reg%NELEM + reg%NGAS + reg%NGAS + 1
   call pmaxmin('CARMA::h2so4_old_1:     ', qa(n)%data3d(i1:i2,j1:j2,1:km), qmin, qmax, ijl, km, 1. )
   n = reg%NBIN*reg%NELEM + reg%NGAS + reg%NGAS + reg%NGAS + 1
   call pmaxmin('CARMA::satliq2_old_1:   ', qa(n)%data3d(i1:i2,j1:j2,1:km), qmin, qmax, ijl, km, 1. )
   n = reg%NBIN*reg%NELEM + reg%NGAS + reg%NGAS + reg%NGAS + reg%NGAS + 1
   call pmaxmin('CARMA::satice2_old_1:   ', qa(n)%data3d(i1:i2,j1:j2,1:km), qmin, qmax, ijl, km, 1. )
endif
#endif
!  =================  END CARMA Step ================================

!  Deallocate
!  --------
   deallocate(p, xc, dx, yc, dy, gasold_, told_, satice_, satliq_, &
              p_, ple_, tmpu_, zc_, zl_, q_, rh_, nuc_, sarea_, numd_, &
              r_wet_, reff_num, reff_den, vf_, zsubsteps_, suvf, mxvf, stat=ios)


! ------------------------------------------------------------------------
! Export states
! ------------------------------------------------------------------------



   RETURN

 END SUBROUTINE CARMA_GridCompRun

!-------------------------------------------------------------------------
!NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1, GEOS/DAS!
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  CARMA_GridCompFinalize
!
! !INTERFACE:
!

   SUBROUTINE CARMA_GridCompFinalize ( gcCARMA, impChem, expChem, &
                                       nymd, nhms, cdt, rc )

  IMPLICIT none

! !INPUT/OUTPUT PARAMETERS:

   TYPE(CARMA_GridComp), INTENT(inout) :: gcCARMA ! Grid Component

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

   CHARACTER(LEN=*), PARAMETER   :: IAm = 'CARMA_GridCompFinalize'
   integer :: STATUS

   rc=0

   deallocate ( gcCARMA%carma )

   RETURN

 END SUBROUTINE CARMA_GridCompFinalize
 

 
  subroutine dumpElement(carma, rc)

  type(carma_type), intent(in)     :: carma              !! the carma object
  integer, intent(inout)           :: rc                 !! return code, negative indicates failure
  
  ! Local Variables
  integer                          :: i
	
  write(*,*)  ""
  write(*,*)  "Element Information"
  
  do i = 1, carma%f_NELEM
    call CARMAELEMENT_Print(carma, i, rc)
    if (rc /=0) write(carma%f_LUNOPRT, *) "    *** FAILED ***, rc=", rc
    write(carma%f_LUNOPRT,*) ""  
  end do
 
  write(carma%f_LUNOPRT,*) ""
  return
  end subroutine dumpElement




  subroutine dumpGas(carma, rc)

  type(carma_type), pointer, intent(inout)  :: carma              !! the carma object
  integer, intent(inout)                    :: rc                 !! return code, negative indicates failure
  
  ! Local Variables
  integer                       :: i
  type(carmagas_type), pointer  :: cgas
  character(len=255)            :: gasname
  real(kind=f)                  :: gwtmol
	
  write(*,*)  ""
  write(*,*)  "Gas Information"
  
  do i = 1, carma%f_NGAS
!   call CARMA_GetGas(carma, i, cgas, rc)
   if (rc /=0) write(*, *) "    *** FAILED ***, rc=", rc
   
!   call CARMAGAS_Print(cgas, carma, rc)
   if (rc /=0) write(*, *) "    *** FAILED ***, rc=", rc
   
   write(*,*) ""  
  end do
 
  write(*,*) ""  
  end subroutine dumpGas




  subroutine dumpGroup(carma, rc)

  type(carma_type), intent(in)     :: carma              !! the carma object
  integer, intent(inout)           :: rc                 !! return code, negative indicates failure
  
  ! Local Variables
  integer                          :: i
	
  write(*,*)  ""
  write(*,*)  "Group Information"
  
  do i = 1, carma%f_NGROUP
    call CARMAGROUP_Print(carma, i, rc)
    if (rc /=0) write(carma%f_LUNOPRT, *) "    *** FAILED ***, rc=", rc
    
    write(carma%f_LUNOPRT,*) ""  
  end do
 
  write(carma%f_LUNOPRT,*) ""
  return
  end subroutine dumpGroup

 END MODULE CARMA_GridCompMod

