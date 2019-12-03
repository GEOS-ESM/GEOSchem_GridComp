!
! Implements Local Displacement Ensembles. It can handle the GEOS-5 lat/lon
! or cubed-sphere grids.
!
! Arlindo da Silva <arlindo.dasilva@nasa.gov>, April 2010
! Cubed sphere support added July 2012.
!----------------------------------------------------------------------------

#  include "MAPL_Generic.h"

module LDE_Mod

   use ESMF
   use MAPL

   use Chem_SimpleBundleMod
   use m_Random
   use m_MergeSorts

   implicit NONE

   private
   public LDE
   public LDE_Create
   public LDE_Destroy
   public LDE_Generate
   public LDE_Projector1c

   type LDE

      type(ESMF_Config), pointer :: CF

      type(ESMF_Grid), pointer :: Grid
      type(ESMF_VM)            :: VM

      integer          :: IM_World=-1   ! global number of lons
      integer          :: JM_World=-1   ! global number of lats
      integer          :: EM_World=-1   ! Max Size of ensemble dimension = Nx*Ny

      real             :: channel             ! single channel to analyze

      real             :: Delta=-1            ! weight parameter 

      integer          :: ks=1                ! top vertical layer

      integer          :: EM=-1         ! Number of esembles to keep
      integer          :: Nx=-1, Ny=-1  ! Stencil (lon,lat) sizes
      real             :: R=-1.0        ! stencil Radius

      logical          :: isCubed       ! Either Cubed Sphere or LatLon

!     Lat/Lon indices
!     ---------------
      integer, pointer :: Ie(:,:) => null()   ! (JM_World,EM_World)
      integer, pointer :: Je(:,:) => null()   ! (JM_World,EM_World)

!     Cubed Sphere Indices
!     --------------------
      integer, pointer :: Indx(:) => null()      ! (EM_WORLD)

   end type LDE

   interface LDE_Generate
      module procedure LDE_Generate2d
   end interface

   interface LDE_Projector
      module procedure LDE_Projector1c
   end interface

   integer, parameter :: OCEAN = 0, LAND = 1, SEA_ICE = 2

   include "mpif.h"
CONTAINS

   subroutine LDE_Create ( self, CF, Grid, rc )
!
!    Initialize ensemble parameters, including ensemble size and
!    ensemble indices (on root PE only). Notice that it is implicitly 
!    assumed a GEOS-5 lat/lon grid for now. 
!
     type(LDE)        , intent(inout)         :: self
     type(ESMF_Config), intent(inout), target :: CF 
     type(ESMF_Grid),   intent(inout), target :: Grid
     integer, intent(out)                     :: rc
!                       ---

                        __Iam__('LDE_Create')

     integer :: dims(3)

     self%Grid => Grid

!    World coordinates
!    -----------------
     call MAPL_GridGet ( grid, globalCellCountPerDim=DIMS, __RC__ )
     self%IM_WORLD = dims(1)
     self%JM_WORLD = dims(2)

!    Cubed sphere or lat/lon grid?
!    -----------------------------
     if ( self%JM_WORLD == self%IM_WORLD * 6 ) then ! Heuristic
          self%isCubed = .TRUE.
     else
          self%isCubed = .FALSE.
     end if
  
!    Stencil properties
!    ------------------
     call ESMF_ConfigGetAttribute(CF, self%R,  Label='stencil_radius_in_km:', __RC__ )
     self%R = 1000. * self%R  ! now in meters

!    Get ensemble size
!    -----------------
     if ( self%isCubed ) then
        call getEnsSizeCubed_ ( self%IM_World, self%JM_World, self%R, self%Nx, self%Ny, __RC__ )
        self%EM_World = self%Nx * self%Ny
     else
        call getEnsSizeLatLon_ ( self%IM_World, self%JM_World, self%R, self%Nx, self%Ny, __RC__ )
        self%EM_World = nint(0.75 * self%Nx * self%Ny) ! pi/4 = 0.7853 ~ 0.75
     end if

!    User may elect a smaller number of ensemble members
!    ---------------------------------------------------
     call ESMF_ConfigGetAttribute(CF, self%EM, default=self%EM_World, &
                                  Label='number_of_ensemble_members:', __RC__ )
     self%EM = min(self%EM,self%EM_World)

!    Get VM for later
!    ----------------
     call ESMF_VMGetGlobal(self%VM,__RC__)

!    Set Ensemble indices on root PE
!    -------------------------------
     if ( self%isCubed ) then

        allocate(self%Indx(self%EM_World), __STAT__ )
        
        if ( MAPL_am_I_root() ) then
             call getEnsIndicesCubed_ ( self%EM_World, self%Indx, __RC__ )
        endif

        call MAPL_CommsBcast (self%vm, self%Indx, size(self%Indx), 0, __RC__)

!    Lat/Lon Grid
!    ------------
     else

        allocate(self%Ie(self%JM_World,self%EM_World), &
                 self%Je(self%JM_World,self%EM_World), __STAT__ )
        
        if ( MAPL_am_I_root() ) then
             call getEnsIndicesLatLon_ ( self%IM_World, self%JM_World, self%EM_World, &
                                         self%R, self%Ie, self%Je, __RC__)
        endif

        call MAPL_CommsBcast (self%vm, self%Ie, size(self%Ie), 0, __RC__)
        call MAPL_CommsBcast (self%vm, self%Je, size(self%Je), 0, __RC__)

     end if

!    Land/ocean channels for single channel formulation
!    --------------------------------------------------
     call ESMF_ConfigGetAttribute(CF, self%channel, Label='single_channel:',  __RC__) 

!    Weight parameter (Delta<0 means do not weight ensemble members)
!    ---------------------------------------------------------------
     call ESMF_ConfigGetAttribute(CF, self%delta, Label='aod_weight_delta:', &
                                                  default=-1.0, __RC__)

!    Top vertical layer
!    ------------------
     call ESMF_ConfigGetAttribute(CF, self%ks, Label='top_vertical_layer:', &
                                               default=1, __RC__)

     if ( MAPL_AM_I_Root() ) then
        print *
        if ( self%isCubed ) then
           print *, 'Initialized LDE on Cubed Sphere with ', self%em, ' ensemble members'
        else
           print *, 'Initialized LDE on LatLon with ', self%em, ' ensemble members'
        end if
     end if

   end subroutine LDE_Create

!.........................................................................

   subroutine LDE_Destroy ( self, rc )
!
!    Initialize ensemble parameters, including ensemble size and
!    ensemble indices (on root PE only). Notice that it is implicitly 
!    assumed a GEOS-5 lat/lon grid for now. 
!
     type(LDE)        , intent(inout)      :: self
     integer, intent(out)                  :: rc
!                       ---

                   __Iam__('LDE_Destroy')

     rc = 0
     self%IM_WORLD = -1
     self%JM_WORLD = -1
     self%EM       = -1
     self%Nx       = -1
     self%Ny       = -1
     self%R        = -1.0

     if ( self%isCubed ) then
        deallocate(self%Indx, __STAT__ )
     else
        deallocate(self%Ie, self%Je, __STAT__ )
     end if

   end subroutine LDE_Destroy

!.........................................................................

   subroutine LDE_Generate2d ( self, e, a, rc )
!
!    Generates LDE based on a 2D horizontal (lon,lat) array.
!
     type(LDE), intent(in) :: self
     real, pointer         :: a(:,:)     ! input distributed horizontal array
     real, pointer         :: e(:,:,:)   ! ensemble of distributed horizontal array
     integer, intent(out)  :: rc

!                             ---

                            __Iam__('LDE_Generate2d')

   if ( self%isCubed ) then 
        call LDE_Generate2d_Cubed_ ( self, e, a, rc )
   else
        call LDE_Generate2d_LatLon_ ( self, e, a, rc )
   end if

   end subroutine LDE_Generate2d

!.........................................................................

  subroutine LDE_Qinc_Global ( self, x_d, x_f,  V, rc )
!
!    This is a single PE routine for computing LDEs and computing the
!    Q analysis increments.
!    All arrays are global.
!
     type(LDE), intent(in) :: self
!     real(kind=ESMF_KIND_R8), pointer :: x_d(:,:) ! Output: Increments
     real(kind=ESMF_KIND_R4), pointer :: x_d(:,:) ! Output: Increments
     real, pointer         :: x_f(:,:)   ! Input:  Bkg on a single level
     real, pointer         :: V(:,:,:)   ! Input:  Ensemble of V 
     integer, intent(out)  :: rc

!                             ---

                            __Iam__('LDE_Qinc_Global')

   if ( self%isCubed ) then 
!ALT        call LDE_Qinc_Global_Cubed_ ( self, x_d, x_f,  V, rc )
   else
        call LDE_Qinc_Global_LatLon_ ( self, x_d, x_f,  V, rc )
   end if

 end subroutine LDE_Qinc_Global

! .............................................................................

   subroutine LDE_Projector1c ( self, bQ_a, bQ_f, bY_f, bY_d, verbose, rc )
!
!    Uses Lagrangian Displacement Ensembles to produce aerosol mixing ratio
!    analysis given AOD (or log-transformed AOD) background and analysis 
!    increments, along with the concentrations background.
!
!    This is the SINGLE CHANNEL version, with same channel used over land 
!    and ocean.

     type(LDE),               intent(inout) :: self
     type(MAPL_SimpleBundle), intent(inout) :: bQ_f ! Aerosol concentration Background
     type(MAPL_SimpleBundle), intent(inout) :: bY_f ! AOD background (or log-AOD)
     type(MAPL_SimpleBundle), intent(inout) :: bY_d ! AOD increment  (or log-AOD)

     type(MAPL_SimpleBundle), intent(inout) :: bQ_a ! Aerosol concentration Analysis; may share
                                                    !  storage with bQ_f

     logical,   OPTIONAL,     intent(in)    :: verbose
     integer,                 intent(out)   :: rc  ! error code


   if ( self%isCubed ) then 
        call LDE_Projector1c_Cubed_ ( self, bQ_a, bQ_f, bY_f, bY_d, verbose, rc )
   else
        call LDE_Projector1c_Latlon_ ( self, bQ_a, bQ_f, bY_f, bY_d, verbose, rc )
   end if

   end subroutine LDE_Projector1c

! .............................................................................

!                       ---------------------
!                       Lat/Lon Grid Routines
!                       ---------------------

   subroutine getEnsSizeLatLon_ ( im, jm, R, Nx, Ny, rc )
     integer, intent(in) :: im
     integer, intent(in) :: jm 
     real,    intent(in)  :: R       ! stencil radius
     integer, intent(out) :: Nx      ! stencil size in lon, at the equator
     integer, intent(out) :: Ny      ! stencil size in lat, away from poles
     integer, intent(out) :: rc 
     
     real*8 :: dx, dy

                          __Iam__('getEnsSizeLatLon_')

     rc = 0     
     dx = 2. * MAPL_Radius * MAPL_PI / im
     dy = MAPL_Radius * MAPL_PI / ( jm - 1 )
     Nx = 2 * (nint(R/dx)-1) + 1 
     Ny = 2 * (nint(R/dy)-1) + 1 

     if ( Nx<3 .OR. Ny<3 ) then
          rc = 1
          return
     end if
     
     if ( Nx*Ny > im*jm ) then
          rc = 2
          return
     end if
     
   end subroutine getEnsSizeLatLon_

   subroutine getEnsIndicesLatLon_ ( im, jm, em, R, Ie, Je, rc )
     integer, intent(in)  :: im
     integer, intent(in)  :: jm 
     integer, intent(in)  :: em      ! ensemble size 
     real,    intent(in)  :: R       ! stencil radius
!                                   --- Ensemble Coordinates ----
     integer, intent(out) :: Ie(jm,em) ! symmetric in longitude
     integer, intent(out) :: Je(jm,em) ! symmetric in longitude
     integer, intent(out) :: rc 
!                      ----

     real*8 :: lat, lon, dLat, dLon ! in radians
     real*8 :: coslon(im), sinlon(im), coslat(jm), sinlat(jm)
     real*8 :: x, y, z, xs, ys, zs, D2_, d2
     integer :: dJ, i, j, is, js, Nx, Ny, ne
     integer :: Kx, Mx, j_deficit, n_deficit

     integer :: Ie_(im*jm), Je_(im*jm), indx(im*jm) ! Ensemble index before thinning
     real*8  :: rn(im*jm)               ! random numbers

                        __Iam__('getEnsIndicesLatLon_')

     rc = 0

!    Check consistency of ensemble size
!    ----------------------------------
     call getEnsSizeLatLon_ ( im, jm, R, Nx, Ny, __RC__ )
     if ( em > Nx * Ny ) then
        print *, trim(Iam)//': inconsistent em, Nx, Ny', em, Nx, Ny
        STATUS = 1
        VERIFY_(STATUS)
     end if

!    Trig - Assumes GEOS-5 lat/lon grid
!    ----------------------------------
     dLon = 2. * MAPL_PI / im
     dLat = MAPL_PI / ( jm - 1 )
     do i = 1, im
        lon = -MAPL_pi + i * dLon
        coslon(i) = cos(lon)
        sinlon(i) = sin(lon)
     end do
     do j = 1, jm
        lat = -MAPL_pi/2. + (j-1)*dLat
        coslat(j) = cos(lat)
        sinlat(j) = sin(lat)
     end do

!    Half threshold distance squared on units of radius
!    --------------------------------------------------
     D2_ = (R / MAPL_Radius)**2

!    Build patch around (is,js) with points that are within L
!    --------------------------------------------------------
     Ie = 0
     Je = 0
     dJ = Ny / 2  ! Ny is always odd

     call zufalli(0) ! initialize random number generator with default seed

     is = 1 ! symmetric in longitude 
     js_: do js = 1, jm

!       Initialize indices for this latitude
!       ------------------------------------
        iE_ = 0
        jE_ = 0
        ne  = 0  

!       Coordinates of reference point on unit sphere
!       ---------------------------------------------
        xs = coslat(js) * coslon(is)
        ys = coslat(js) * sinlon(is)
        zs = sinlat(js)
        
!       Find those grid points that are close enough
!       --------------------------------------------
        jj_: do j = js-dJ, js+dJ
           if ( j < 1 .OR. j > jm ) cycle jj_
           ii_: do i = 1, im                       ! it works at the poles as well
              if ( i==is .AND. j==js ) cycle ii_   ! eliminate central point
              x = coslat(j) * coslon(i)
              y = coslat(j) * sinlon(i)
              z = sinlat(j)
              d2 = 2. * abs(1.0 - (x*xs + y*ys + z*zs)) ! chordal distance
              if ( d2 <= D2_ ) then
                 ne = ne + 1
                 iE_(ne) = i  ! record this longitude index
                 jE_(ne) = j  ! record this latitude  index
              end if
           end do ii_ 
        end do jj_

!       Consistency check, should never happen
!       --------------------------------------
        if ( ne < em ) then ! recall that we skipped middle point (ZERO perturbation)
           print *, trim(Iam)//': not enough ensemble members: ', &
                js, ne, em, (em-ne) 
           STATUS = 3
           VERIFY_(STATUS)
        end if

!       Final shuffle so that we can select fewer members later,
!       say, the first 100 member will be uniformly distributed
!       in space. ZERO member will be added to the end
!       NOTE: Sampling is biased towards polar latitudes
!       -----------------------------------------------------
        call zufall ( ne, rn) 
        call IndexSet ( ne, indx )
        call IndexSort ( ne, indx, rn, descend=.false.)
        Ie(js,1:em) = Ie_ ( (/ (indx(i), i=1,em) /) )
        Je(js,1:em) = Je_ ( (/ (indx(i), i=1,em) /) )

     end do js_

     rc = 0

   end subroutine getEnsIndicesLatLon_

!.........................................................................

   subroutine LDE_Generate2d_LatLon_ ( self, e, a, rc )
!
!    Generates LDE based on a 2D horizontal (lon,lat) array.
!
     type(LDE), intent(in) :: self
     real, pointer         :: a(:,:)     ! input distributed horizontal array
     real, pointer         :: e(:,:,:)   ! ensemble of distributed horizontal array
     integer, intent(out)  :: rc

!                             ---

                            __Iam__('LDE_Generate2d_LatLon')

     type(ESMF_Grid), pointer :: Grid

!    Global version of arrays (root PE only) 
!    ---------------------------------------
     integer :: IM_World, JM_World, EM, i, j, n
     real, pointer :: a_world(:,:)    => null() 
     real, pointer :: e_world(:,:)    => null() 
     real, pointer :: ie(:)           => null() 
     real, pointer :: je(:)           => null() 

     Grid => self%Grid

!    Allocate buffers
!    ----------------
     IM_World = self%IM_World
     JM_World = self%JM_World
     EM = self%EM
     allocate(a_world(IM_World,JM_World), &
              e_world(IM_World,JM_World), &
              ie(JM_World), je(JM_World), &
              __STAT__)

!    Gather input array
!    ------------------
     call ArrayGather ( a, a_world, Grid, __RC__ )

!    Generate LDE
!    ------------
Ens: do n = 1, EM

!       Generate this ensemble member on root PE
!       ----------------------------------------
        if ( MAPL_AM_I_ROOT() ) then
           je = self%je(:,n)
zonal:     do i = 1, IM_World
              ie = self%ie(:,n) + i - 1 ! shift zonal indices
              where ( ie <        1 ) ie = IM_World + ie
              where ( ie > IM_World ) ie = ie - IM_World
merid:        do j = 1, JM_World
                 e_world(i,j) = a_world(ie(j),je(j)) 
              end do merid
           end do zonal
!!!           e_world = e_world - a_world ! displacement from central point
        end if

!       Scatter this member
!       -------------------
        call ArrayScatter ( e(:,:,n), e_world, Grid, __RC__ )

!       Compute displacement from central point
!       ---------------------------------------
        e(:,:,n) = e(:,:,n) - a(:,:)  

     end do Ens

!    Free memory
!    -----------
     deallocate(a_world,e_world,ie,je,__STAT__)

   end subroutine LDE_Generate2d_LatLon_

! .............................................................................

   subroutine LDE_Projector1c_LatLon_ ( self, bQ_a, bQ_f, bY_f, bY_d, verbose, rc )
!
!    Uses Lagrangian Displacement Ensembles to produce aerosol mixing ratio
!    analysis given AOD (or log-transformed AOD) background and analysis 
!    increments, along with the concentrations background.
!
!    This is the SINGLE CHANNEL version, with same channel used over land 
!    and ocean.
!
!    IMPORTANT: This routine also works for the cubed spehere, albeit not as
!               efficiently as LDE_Projector1c_Cubed_(). Once the cubed version is
!               made to handle Lat-Lon as well we should drop this rotuine.
!
     type(LDE),               intent(inout) :: self
     type(MAPL_SimpleBundle), intent(inout) :: bQ_f ! Aerosol concentration Background
     type(MAPL_SimpleBundle), intent(inout) :: bY_f ! AOD background (or log-AOD)
     type(MAPL_SimpleBundle), intent(inout) :: bY_d ! AOD increment  (or log-AOD)

     type(MAPL_SimpleBundle), intent(inout) :: bQ_a ! Aerosol concentration Analysis; may share
                                                    !  storage with bQ_f

     logical,   OPTIONAL,     intent(in)    :: verbose
     integer,                 intent(out)   :: rc  ! error code

     integer :: i, j, k, e, s, im, jm, km, em
     integer :: ifAOD, idAOD
     logical :: verbose_, missing_f, missing_d

!     real(kind=ESMF_KIND_R8), pointer :: vnorm(:,:), x_d(:,:,:)  ! accumulators
     real(kind=ESMF_KIND_R8), pointer :: vnorm(:,:)  ! accumulators
     real(kind=ESMF_KIND_R4), pointer :: x_d(:,:,:)  ! accumulators
     real, pointer                    :: q_f(:,:), y_f(:,:), y_d(:,:)  ! 2D single instances
     real, pointer                    :: X(:,:,:), V(:,:,:), W(:,:,:)  ! 2D ensemble variables

     real(kind=ESMF_KIND_R4), pointer :: x_d_World(:,:)
     real(kind=ESMF_KIND_R4), pointer :: x_d_World3d(:,:,:)
     real, pointer                    :: q_f_World(:,:) 
     real, pointer                    :: q_f_World3d(:,:,:) 
     real, pointer                    :: V_World(:,:,:)
     real, pointer                    :: a(:,:), a_World(:,:)
     integer, ALLOCATABLE             :: krank(:)
     integer                          :: mype, npes, nn, color, comm, lde_comm
     integer                          :: nnodes

                        __Iam__('LDE_Projector1c')
                        
    if ( present(verbose) ) then
       verbose_ = verbose
    else
       verbose_ = .FALSE.
    end if

    im = size(bQ_f%r3(1)%q,1)
    jm = size(bQ_f%r3(1)%q,2)
    km = ubound(bQ_f%r3(1)%q,3)
    em = self%em

!   Allocate workspace
!   ------------------
    allocate ( y_f(im,jm),   &
               y_d(im,jm),   & 
               vnorm(im,jm), &
               X(im,jm,em),  &
               V(im,jm,em),  &
               W(im,jm,em),  &
               __STAT__ )

     allocate(x_d(im,jm,self%ks:km),  __STAT__)

!    Determine convenience indices
!    -----------------------------
     ifAOD = MAPL_SimpleBundleGetIndex(bY_f,'AOD',3,__RC__)
     idAOD = MAPL_SimpleBundleGetIndex(bY_d,'AOD',3,__RC__)

!    Use single channel
!    ------------------
     _ASSERT(size(bY_f%coords%levs) == size(bY_d%coords%levs),'needs informative message')
     missing_f = .TRUE.
     missing_d = .TRUE.
     do k = 1, size(bY_f%coords%levs)
        if ( abs(bY_f%coords%levs(k)-self%channel) < 0.01 ) then
             y_f = bY_f%r3(ifAOD)%q(:,:,k)
             missing_f = .FALSE.
        end if
        if ( abs(bY_d%coords%levs(k)-self%channel) < 0.01 ) then
             y_d = bY_d%r3(idAOD)%q(:,:,k)
             missing_d = .FALSE.
        end if
     end do
     if ( missing_f ) then
        __raise__(MAPL_RC_ERROR,"could not find matching channel for <y_f>")
     end if
     if ( missing_d ) then
        __raise__(MAPL_RC_ERROR,"could not find matching channel for <y_d>")
     end if

#ifdef DEBUG
     if ( MAPL_AM_I_Root() .and. verbose_ ) print *
     call MAPL_MaxMin('y_f',y_f)
     call MAPL_MaxMin('y_d',y_d)
#endif

!    Generate ensembles of AOD backgrounds
!    -------------------------------------
     call LDE_Generate2d ( self, V, y_f, __RC__ ) 

#ifdef DEBUG
     if ( MAPL_AM_I_Root() .and. verbose_ ) print *
     call MAPL_MaxMin(' V ',V)
#endif

!    Create ensemble weights
!    -----------------------
     if ( self%Delta <= 0.0 ) then
        W = 1.0 ! ensemble members are equal-probable
     else
        do e = 1, em
           vnorm = ((V(:,:,e)-y_d(:,:))/self%Delta)**2
           where(vnorm<20.) ! underflow protection
              W(:,:,e) = exp(-vnorm)
           elsewhere
              W(:,:,e) = exp(-20.)
           end where
        end do
     end if

#ifdef DEBUG
     call MAPL_MaxMin(' W ',W)
#endif
        
!    Normalized AOD ensembles
!
!      v{e} = y_f{e} * y_d / <y_f,y_f>
!
!    for each ensemble member {e}
!    ---------------------------------
     vnorm = 0.0
     do e = 1, em
        vnorm = vnorm + W(:,:,e) * V(:,:,e)**2
     end do
     where ( vnorm==0.0 ) ! division by zero protection
        y_d = 0.0
     elsewhere
        y_d = y_d / vnorm
     end where
     do e = 1, em
        V(:,:,e) = W(:,:,e) * V(:,:,e) * y_d(:,:)
     end do

#ifdef DEBUG
     call MAPL_MaxMin(' V ',V)
#endif

!    Gather V to all processes that will participate in the analysis
!    First we make a sub-communicator containing those processes (lde_comm)
!    ----------------------------------------------------------------------
     call ESMF_VMGet       (self%VM, mpiCommunicator=COMM, localpet=MYPE, petcount=NPES,  __RC__)
     allocate(krank(self%ks:km))
     nNodes = size(MAPL_NodeRankList)
     call MAPL_RoundRobinPEList(krank, nNodes, __RC__)
     color = MPI_UNDEFINED
     do k = self%ks, km
        if( krank(k)==mype ) then
           color = 0
        endif
     enddo

     CALL MPI_COMM_SPLIT(COMM, color, mype, lde_comm, STATUS)
     VERIFY_(STATUS)

!    Allocate V_World on all processes that will participate in the analysis
!    -----------------------------------------------------------------------
     if(color /= MPI_UNDEFINED) then
        allocate(V_World(self%IM_World,self%JM_World,EM), __STAT__)
     endif

!    Gather the distributed V to form the global V_World
!    All process participate in this call
!    ArrayGather gathers to MAPL_Root by default. To be safe,
!    we specify the first rank from MAPL_RoundRobinPEList
!    --------------------------------------------------------
     do e = 1, em
        call ArrayGather(V(:,:,e),V_World(:,:,e), self%Grid, depe=krank(self%ks), __RC__)
     end do

!    Now broadcast from the gather rank to all others in the lde_comm communicator
!    -----------------------------------------------------------------------------
     if (color /= MPI_UNDEFINED) &
         call MPI_Bcast(V_World, size(V_World), MPI_REAL, krank(self%ks), lde_comm, STATUS)

!    Next compute aerosol concentration analysis for each level, species
!             q_a = q_f + <X,V>
!    -------------------------------------------------------------------
     if ( MAPL_AM_I_Root() .and. verbose_ ) then
        if ( self%isCubed ) then
           print *, 'Calculating LDE increments on Cubed Sphere with ', em, ' ensemble members'
        else
           print *, 'Calculating LDE increments on LatLon with ', em, ' ensemble members'
        end if
     end if

     do s = 1, bQ_f%n3d
        if ( .not. isAerosol_(trim(bQ_f%r3(s)%name)) ) cycle

        if ( MAPL_AM_I_Root() .and. verbose_ ) &
           print *, ' [ ] Working on <'//trim(bQ_f%r3(s)%name)//'>'

#ifdef DEBUG
        call MAPL_MaxMin('      q_f',bQ_f%r3(s)%q(:,:,self%ks:km))
#endif

!       Gather distributed levels onto a single processor
!       Level to processor assignment occurs inside MAPL_CollectiveGather3D
!       Our call to MAPL_RoundRobinPEList is assumed to mimic the assignment
!       achieved in MAPL_CollectiveGather3D
!       --------------------------------------------------------------------
        call MAPL_CollectiveGather3D(self%Grid, bQ_f%r3(s)%q(:,:,self%ks:km), &
                                     q_f_World3d, __RC__)

!       Allocate work space depending on level assignment
!       -------------------------------------------------
        allocate(x_d_World3d(self%IM_World,self%JM_World,size(q_f_World3d,3)),  __STAT__)

!       Each process does the analysis on its assigned level
!       ----------------------------------------------------
        nn=0
        do k = self%ks, km
           if( krank(k)==mype ) then
              nn=nn+1
              x_d_World => x_d_World3d(:,:,nn)
              q_f_World => q_f_World3d(:,:,nn)
              call LDE_Qinc_Global(self, x_d_World, q_f_World, V_World, __RC__)
           endif
        end do
        deallocate(q_f_World3d)

!       Scatter the analysis
!       --------------------
        call MAPL_CollectiveScatter3D(self%Grid,  x_d_World3d(:,:,:nn), x_d(:,:,self%ks:km), &
                                      __RC__)
        deallocate(x_d_World3d)

!       Add analysis increments to q
!       ----------------------------
        do k = self%ks, km
           bQ_a%r3(s)%q(:,:,k) = bQ_f%r3(s)%q(:,:,k) + x_d(:,:,k)
        end do

!       Zero increments above top analysis level
!       ----------------------------------------
        do k = 1,self%ks-1 
           bQ_a%r3(s)%q(:,:,k) = bQ_f%r3(s)%q(:,:,k)
        end do
        
#ifdef DEBUG
        call MAPL_MaxMin('      q_a',bQ_a%r3(s)%q(:,:,self%ks:km))
#endif

        bQ_a%r3(s)%q(:,:,self%ks:km) = max(0.0,bQ_a%r3(s)%q(:,:,self%ks:km)) ! fix q<0

     end do ! variable loop

!    Free memory
!    -----------
     deallocate ( y_f, y_d, vnorm, X, V, W, x_d, __STAT__ )

     if(color /= MPI_UNDEFINED) deallocate ( V_World, __STAT__)
                                deallocate ( krank  , __STAT__)

     rc = 0
     return

   end subroutine LDE_Projector1c_LatLon_

!.........................................................................

   subroutine LDE_Qinc_Global_LatLon_ ( self, x_d, x_f,  V, rc )
!
!    This is a single PE routine for computing LDEs and computing the
!    Q analysis increments.
!    All arrays are global.
!
     type(LDE), intent(in) :: self
!     real(kind=ESMF_KIND_R8), pointer :: x_d(:,:) ! Output: Increments
     real(kind=ESMF_KIND_R4), pointer :: x_d(:,:) ! Output: Increments
     real, pointer         :: x_f(:,:)   ! Input:  Bkg on a single level
     real, pointer         :: V(:,:,:)   ! Input:  Ensemble of V 
     integer, intent(out)  :: rc

!                             ---

                            __Iam__('LDE_Qinc_Global_LatLon')

     integer :: IM_World, JM_World, EM, i, j, n
     real, pointer :: X(:,:) => null()  ! LDE based on x_f
     real, pointer :: ie(:)  => null() 
     real, pointer :: je(:)  => null() 

!    Allocate buffers
!    ----------------
     IM_World = self%IM_World
     JM_World = self%JM_World
     EM = self%EM
     allocate(X(IM_World,JM_World), &
              ie(JM_World), je(JM_World), &
              __STAT__)

!    Generate LDE
!    ------------
     x_d = 0.0
Ens: do n = 1, EM

!       Generate this ensemble member on root PE
!       ----------------------------------------
        je = self%je(:,n)
zonal:  do i = 1, IM_World
           ie = self%ie(:,n) + i - 1 ! shift zonal indices
           where ( ie <        1 ) ie = IM_World + ie
           where ( ie > IM_World ) ie = ie - IM_World
merid:     do j = 1, JM_World
              X(i,j) = x_f(ie(j),je(j)) 
           end do merid
        end do zonal

        x_d(:,:) = x_d(:,:) + (X(:,:) - x_f(:,:)) * V(:,:,n)

     end do Ens

!    Free memory
!    -----------
     deallocate(X,ie,je,__STAT__)

   end subroutine LDE_Qinc_Global_LatLon_


!                       ---------------------
!                       Cubed-Sphere Routines
!                       ---------------------

   subroutine getEnsSizeCubed_ ( im, jm, R, Nx, Ny, rc )
     integer, intent(in) :: im
     integer, intent(in) :: jm 
     real,    intent(in)  :: R       ! stencil radius

     integer, intent(out) :: Nx      ! stencil size in X, always odd
     integer, intent(out) :: Ny      ! same as Nx for cubed-sphere
     integer, intent(out) :: rc
     
!                                 ---
     real*8 :: dArea, dx

                          __Iam__('getEnsSizeCubed_')

!  Compute average grid box sizes
!  ------------------------------
   dArea = 4. * MAPL_PI * MAPL_RADIUS**2 / ( im * jm ) ! mean area
   dx = sqrt(dArea) ! assumes square "mean" grid-boxes
     
!  Compute stencil size, making sure it is a odd number for symmetry
!  -----------------------------------------------------------------
   Nx = 2 * nint(R/dx) + 1
   Ny = Nx

   if ( Nx<3 ) then
      rc = 1
   end if

   end subroutine getEnsSizeCubed_

!............................................................
   subroutine getEnsIndicesCubed_ ( EM_World, Indx, rc )
     integer, intent(in)  :: EM_World        ! maximum ensemble size 
     integer, intent(out) :: Indx(EM_World)  ! Randomized indices
     integer, intent(out) :: rc 
!                      ----

     real*8  :: rn(EM_World)               ! random numbers

                        __Iam__('getEnsIndicesCubed_')

     rc = 0
     call zufalli(0) ! initialize random number generator with default seed
     call zufall(EM_World,rn) ! sample 
     call IndexSet ( EM_World, Indx )
     call IndexSort ( EM_WORLD, Indx, rn, descend=.false.)

   end subroutine getEnsIndicesCubed_

!............................................................

   subroutine LDE_HaloedFace_ ( im, nH, iFace, hA, A, rc )

     integer, intent(in)  :: im         ! x/y size for a (square) cube face 
     integer, intent(in)  :: nH         ! number of grid-points in halo
     integer, intent(in)  :: iFace      ! which face of the cube to halo 
     real,    intent(in)  :: A(im,im,6) ! global array on cobed-sphere 

     real, intent(out)    :: hA(-nH+1:im+nH,-nH+1:im+nH)  ! haloed array on face iFace
     integer, intent(out) :: rc
!
!    Given a global array on the cubed sphere, returns a haloed array on a single face.
!    For face 1 we will have:
!
!                    x | 3 | x
!                  ----|---|---    
!                    5 | 1 | 2 
!                  ----|---|---    
!                    x | 6 | x
!
!    where we have indicated the relevant faces. The nearby faces with "x" are the so-called
!    dead zones and values in these regions will be set to UNDEF because they are already
!    included in other tiles. Recall that axis of adjascent faces may or may not to be rotated
!    to align properly. I relied on a paper cut out of the cubed-sphere to have this coded up.
!
!    Example for im=5 and nH=3:
!
!                                    | 1 2 3 4 5 |
!            (LEFT)          -2 -1 0 | 1 2 3 4 5 | 6 7 8       (RIGHT)
!                                    | 1 2 3 4 5 |
!
!                             ----

                        __Iam__('getHaloedFace_')

     rc = 0

!    Start with all UNDEFs, dead zones will remain UNDEF
!    ---------------------------------------------------
     hA = MAPL_UNDEF

!    Fill in the core of the domain
!    ------------------------------
     hA(1:im,1:im) = A(1:im,1:im,iFace)

!    Special handle each face
!    ------------------------
     if ( iFace == 1 ) then
          call fill_ ( im, 6,   0, 'bottom', hA )  
          call fill_ ( im, 3, -90, 'top'   , hA )           
          call fill_ ( im, 5, +90, 'left'  , hA ) 
          call fill_ ( im, 2,   0, 'right' , hA ) 
     else if ( iFace == 2 ) then
          call fill_ ( im, 6, -90, 'bottom', hA )  
          call fill_ ( im, 3,   0, 'top'   , hA )           
          call fill_ ( im, 1,   0, 'left'  , hA ) 
          call fill_ ( im, 4, +90, 'right' , hA ) 
     else if ( iFace == 3 ) then
          call fill_ ( im, 2,   0, 'bottom', hA )  
          call fill_ ( im, 5, -90, 'top'   , hA )           
          call fill_ ( im, 1, +90, 'left'  , hA ) 
          call fill_ ( im, 4,   0, 'right' , hA ) 
     else if ( iFace == 4 ) then
          call fill_ ( im, 2, -90, 'bottom', hA )  
          call fill_ ( im, 5,   0, 'top'   , hA )           
          call fill_ ( im, 3,   0, 'left'  , hA ) 
          call fill_ ( im, 6, +90, 'right' , hA ) 
     else if ( iFace == 5 ) then
          call fill_ ( im, 4,   0, 'bottom', hA )  
          call fill_ ( im, 1, -90, 'top'   , hA )           
          call fill_ ( im, 3, +90, 'left'  , hA ) 
          call fill_ ( im, 6,   0, 'right' , hA ) 
     else if ( iFace == 6 ) then
          call fill_ ( im, 4, -90, 'bottom', hA )  
          call fill_ ( im, 1,   0, 'top'   , hA )           
          call fill_ ( im, 5,   0, 'left'  , hA ) 
          call fill_ ( im, 2, +90, 'right' , hA ) 
     end if

     contains

           subroutine fill_ ( im, jFace, iRot, location, hA ) 
           integer, intent(in) :: im
           integer, intent(in) :: jFace 
           integer, intent(in) :: iRot ! whether or not to rotate array
           character(len=*), intent(in) :: location
           real, intent(out)   :: hA(-nH+1:im+nH,-nH+1:im+nH)  ! haloed array on face iFace

!
!             Fills top, bottom, left or right halo
!
           integer i, j
           real :: x(im,im)

!          Rotate adjascent face as necessary
!          ----------------------------------           
           if ( iRot==0 ) then
              x = A(:,:,jFace)
           else if ( iRot == +90 ) then ! clockwise
              do j = 1, im
                 do i = 1, im
                    x(i,j) = A(j,im-i+1,jFace)
                 end do
              end do
           else if ( iRot == -90 ) then ! counter-clockwise
              do j = 1, im
                 do i = 1, im
                    x(i,j) = A(im-j+1,i,jFace)
                 end do
              end do
           end if

!          Fill in this halo segment
!          -------------------------
           if (      location == 'bottom' ) then 
                                                 hA(1:im,-nh+1:0)    = x(1:im,im-nh+1:im) 
           else if ( location == 'top' )    then
                                                 hA(1:im,im+1:im+nh) = x(1:im,1:nh)
           else if ( location == 'left' )   then
                                                 hA(-nh+1:0   ,1:im) = x(im-nh+1:im,1:im)
           else if ( location == 'right' )  then
                                                 hA(im+1:im+nh,1:im) = x(1:nH,1:im)
           end if

         end subroutine fill_

   end subroutine LDE_HaloedFace_

!........................................................

   subroutine LDE_Qinc_Distrib_Cubed_ ( x_d, a, V, indx, im, jm, em, IM_World, EM_World, nh, self, rc ) 
      integer, intent(in)  :: IM_World ! number of x,y gridpoints on face (global)
      integer, intent(in)  :: EM_World ! maximum ensemble size 

      integer, intent(in)  :: im, jm   ! local dimensions (distributed)
      integer, intent(in)  :: em       ! desired number of ensemble members
      integer, intent(in)  :: nH       ! halo size

      integer, intent(in)  :: indx(EM_World) ! Randomized indices

      real,    intent(in)  :: a(im,jm)  ! first guess (distributed)
      type(LDE), intent(in)      :: self


      real,    intent(out) :: x_d(im,jm) ! Analysis increments (distributed)
      integer, intent(out) :: rc

      real                 :: V(im,jm,em) ! distributed RHS
      real                 :: X(im,jm,em) ! Analysis increments (distributed)
      integer              :: k

      call LDE_Generate2d_Cubed_Core_ ( X, a, indx, im, jm, em, IM_World, EM_World, nh, self, rc ) 

      x_d = 0.0
      do k = 1, em
         x_d(:,:) = x_d(:,:) + X(:,:,k) * V(:,:,k)
      end do
   end subroutine LDE_Qinc_Distrib_Cubed_

   subroutine LDE_Generate2d_Cubed_Core_ ( X, a, indx, im, jm, em, IM_World, EM_World, nh, self, rc ) 

      integer, intent(in)  :: IM_World ! number of x,y gridpoints on face (global)
      integer, intent(in)  :: EM_World ! maximum ensemble size 

      integer, intent(in)  :: im, jm   ! local dimensions (distributed)
      integer, intent(in)  :: em       ! desired number of ensemble members
      integer, intent(in)  :: nH       ! halo size

      integer, intent(in)  :: indx(EM_World) ! Randomized indices

      real,    intent(in)  :: a(im,jm)  ! first guess (distributed)
      type(LDE), intent(in)         :: self

      real,    intent(out) :: X(im,jm,em) ! Analysis increments (distributed)
      integer, intent(out) :: rc
!
!      Returns LDEs given a 2D array. This is a single PE routine; all arrays are global.
!
!                        ----
     real    ::  hA(-nH+1:IM_World+nH,-nH+1:IM_World+nH)  ! haloed array on a single face
     real    ::  e_(EM_World)   ! Maximum possible ensembles for a given point
     integer :: iFace, is, js, i, j, k, nH2, j2, d2, ne
     integer :: ig, jg
     integer :: myFace, i_, j_



     real  :: A_World(IM_World,IM_World*6) ! Input 2D array on the cubed sphere
     real  :: V(im,jm,em) ! distributed RHS
     integer :: i1, in, j1, jn
     logical :: inside_domain
     integer :: status 
     character(len=ESMF_MAXSTR) :: Iam = 'LDE_Generate2d_Cubed_Core'

     rc = 0
     nH2 = nH**2

! gather a to A_World
!--------------------
     call ArrayGather ( a, a_world, self%Grid, __RC__ )
     call MAPL_CommsBcast (self%vm, a_world, size(a_world), 0, __RC__)

     X = 0.0 ! Just in case 

!    Consistency check
!    -----------------
     if ( EM_World /= (2*nH+1)**2 ) then
        print *, 'Very strange: Inconsistent nH, EM_World =', nH, EM_World 
        rc = 1
        return
     end if

!    Loop over faces of the cube
!    ---------------------------
!    Determine myface
! ---------------------
     ! get lower left and upper right corners of my domain
     call MAPL_GRID_INTERIOR(self%Grid,I1,IN,J1,JN)
     
     myFace = (J1-1)/IM_World + 1
     ! Saniny checking: make sure the upper right corner is on the same face
     _ASSERT(myFace == (jn -1)/IM_World+1,'needs informative message')

     face: do iFace = myface, myface

!       Haloed array on this face
!       -------------------------
        call LDE_HaloedFace_ ( IM_World, nH, iFace, hA, A_World, rc )

!       For grid point (is,js) on this face...
!       --------------------------------------
        js_: do js = 1, IM_World
           is_: do is = 1, IM_World

! convert (is,js,iFace) --> global (ig,jg)
              ig = is
              jg = (iFace-1)*IM_World + js
! convert (ig, jg) --> local (i_,j_)
              i_ = ig - i1 + 1
              j_ = jg - j1 + 1

! Check if (is,js,iFace) is in domain
              inside_domain = ig >= i1 .and. ig <=in .and.  jg>=j1 .and. jg<=jn
              if ( .not. inside_domain  ) cycle 

!             Look around for ensemble members
!             --------------------------------
              ne = 0
              e_ = MAPL_UNDEF
              jj_: do j = js-nH, js+nH
                 j2 = (j-js)**2
                 ii_: do i = is-nH, is+nH
                    d2 = j2 + (i-is)**2 ! distance squared from (is,js)
                    ne = ne + 1
                    if ( (i==is .AND. j==js) .OR. (d2>nH2) ) then 
                         e_(ne) = MAPL_UNDEF ! do not include central point or outside "circle"
                    else
                         e_(ne) = hA(i,j)
                    end if
                 end do ii_
              end do jj_

              e_ = e_(indx(:)) ! randomize ensembles, including dead zones

!             Select only "em" defined ensembles (exclude dead zones, central point
!              and points outsize the circle)
!             ---------------------------------------------------------------------
              ne = 0
              ens: do k = 1, EM_World
                      if ( e_(k) /= MAPL_UNDEF ) then ! skip undefined memebers
                           ne = ne + 1
                           if ( ne > em ) exit ens
                           X(i_,j_,ne) = e_(k) - a(i_,j_)
                      end if
              end do ens
              if ( ne < em ) then
                   print *, 'Very strange: not enough ensemble members found - ne, em =', ne, em
                   rc = 2  ! not enough members found, should not happen
                   return
              end if

           end do is_
        end do js_

     end do face

   end subroutine LDE_Generate2d_Cubed_Core_

!..............................................................................

     subroutine LDE_Generate2d_Cubed_ ( self, e, a, rc )
!
!    Generates LDE based on a 2D horizontal (lon,lat) array.
!
     type(LDE), intent(in) :: self
     real, pointer         :: a(:,:)     ! input distributed horizontal array
     real, pointer         :: e(:,:,:)   ! ensemble of distributed horizontal array
     integer, intent(out)  :: rc

!                             ---

                            __Iam__('LDE_Generate2d_Cubed')

     type(ESMF_Grid), pointer :: Grid

     integer :: IM, JM, EM
     integer :: IM_World, JM_World, EM_World, n, nH

!    Aliases
!    ----------------
     Grid => self%Grid
     nH = self%Nx / 2  ! hallo size

     im = size(a,1)
     jm = size(a,2)
     em = self%em

     IM_World = self%IM_World
     JM_World = self%JM_World
     EM_World = self%EM_World

     call LDE_Generate2d_Cubed_Core_ ( e, a, self%indx, im, jm, em, IM_World, EM_World, nh, self, rc ) 


    end subroutine LDE_Generate2d_Cubed_

!.........................................................................

   subroutine LDE_Projector1c_Cubed_ ( self, bQ_a, bQ_f, bY_f, bY_d, verbose, rc )
!
!    Uses Lagrangian Displacement Ensembles to produce aerosol mixing ratio
!    analysis given AOD (or log-transformed AOD) background and analysis 
!    increments, along with the concentrations background.
!
!    This is the SINGLE CHANNEL version, with same channel used over land 
!    and ocean.
!
!    IMPORTANT: This routine does not yet work for Lat-Lon; this could be accomplished
!               by implementing a version of LDE_Qinc_Distrib() that can handle the
!               cubed sphere.
!

     type(LDE),               intent(inout) :: self
     type(MAPL_SimpleBundle), intent(inout) :: bQ_f ! Aerosol concentration Background
     type(MAPL_SimpleBundle), intent(inout) :: bY_f ! AOD background (or log-AOD)
     type(MAPL_SimpleBundle), intent(inout) :: bY_d ! AOD increment  (or log-AOD)

     type(MAPL_SimpleBundle), intent(inout) :: bQ_a ! Aerosol concentration Analysis; may share
                                                    !  storage with bQ_f

     logical,   OPTIONAL,     intent(in)    :: verbose
     integer,                 intent(out)   :: rc  ! error code

!                             ----

     integer :: i, j, k, e, s, im, jm, km, em
     integer :: ifAOD, idAOD
     logical :: verbose_, missing_f, missing_d

!     real(kind=ESMF_KIND_R8), pointer :: vnorm(:,:), x_d(:,:,:)  ! accumulators
     real(kind=ESMF_KIND_R8), pointer :: vnorm(:,:)  ! accumulators
     real(kind=ESMF_KIND_R4), pointer :: x_d(:,:,:)  ! accumulators
     real, pointer                    :: q_f(:,:), y_f(:,:), y_d(:,:)  ! 2D single instances
     real, pointer                    :: X(:,:,:), V(:,:,:), W(:,:,:)  ! 2D ensemble variables
     real, pointer                    :: x_2d(:,:)

     real(kind=ESMF_KIND_R4), pointer :: x_d_World(:,:)
     real(kind=ESMF_KIND_R4), pointer :: x_d_World3d(:,:,:)
     real, pointer                    :: q_f_World(:,:) 
     real, pointer                    :: V_World(:,:,:)
     real, pointer                    :: a(:,:), a_World(:,:)
     integer, ALLOCATABLE             :: krank(:)
     integer                          :: mype, npes, nn, color, comm, lde_comm
     integer                          :: im_world, jm_world, em_world
     integer                          :: nH

                        __Iam__('LDE_Projector1c')
                        
    if ( present(verbose) ) then
       verbose_ = verbose
    else
       verbose_ = .FALSE.
    end if

    im = ubound(bQ_f%r3(1)%q,1)
    jm = ubound(bQ_f%r3(1)%q,2)
    km = ubound(bQ_f%r3(1)%q,3)
    em = self%em
    im_world = self%im_world
    jm_world = self%jm_world
    em_world = self%em_world

    nH = self%Nx/2

!   Allocate workspace
!   ------------------
    allocate ( y_f(im,jm),   &
               y_d(im,jm),   & 
               vnorm(im,jm), &
               X(im,jm,em),  &
               V(im,jm,em),  &
               W(im,jm,em),  &
               __STAT__ )

     allocate(x_d(im,jm,self%ks:km),  __STAT__)

     V = 0.0 ! ALT: Initialize just in case

!    Determine convenience indices
!    -----------------------------
     ifAOD = MAPL_SimpleBundleGetIndex(bY_f,'AOD',3,__RC__)
     idAOD = MAPL_SimpleBundleGetIndex(bY_d,'AOD',3,__RC__)

!    Use single channel
!    ------------------
     _ASSERT(size(bY_f%coords%levs) == size(bY_d%coords%levs),'needs informative message')
     missing_f = .TRUE.
     missing_d = .TRUE.
     do k = 1, size(bY_f%coords%levs)
        if ( abs(bY_f%coords%levs(k)-self%channel) < 0.01 ) then
             y_f = bY_f%r3(ifAOD)%q(:,:,k)
             missing_f = .FALSE.
        end if
        if ( abs(bY_d%coords%levs(k)-self%channel) < 0.01 ) then
             y_d = bY_d%r3(idAOD)%q(:,:,k)
             missing_d = .FALSE.
        end if
     end do
     if ( missing_f ) then
        __raise__(MAPL_RC_ERROR,"could not find matching channel for <y_f>")
     end if
     if ( missing_d ) then
        __raise__(MAPL_RC_ERROR,"could not find matching channel for <y_d>")
     end if

#ifdef DEBUG
     if ( MAPL_AM_I_Root() .and. verbose_ ) print *
     call MAPL_MaxMin('y_f',y_f)
     call MAPL_MaxMin('y_d',y_d)
#endif

!    Generate ensembles of AOD backgrounds
!    -------------------------------------
     call LDE_Generate2d ( self, V, y_f, __RC__ ) 

#ifdef DEBUG
     if ( MAPL_AM_I_Root() .and. verbose_ ) print *
     call MAPL_MaxMin(' V ',V)
#endif

!    Create ensemble weights
!    -----------------------
     if ( self%Delta <= 0.0 ) then
        W = 1.0 ! ensemble members are equal-probable
     else
        do e = 1, em
           vnorm = ((V(:,:,e)-y_d(:,:))/self%Delta)**2
           where(vnorm<20.) ! underflow protection
              W(:,:,e) = exp(-vnorm)
           elsewhere
              W(:,:,e) = exp(-20.)
           end where
        end do
     end if

#ifdef DEBUG
     call MAPL_MaxMin(' W ',W)
#endif
        
!    Normalized AOD ensembles
!
!      v{e} = y_f{e} * y_d / <y_f,y_f>
!
!    for each ensemble member {e}
!    ---------------------------------
     vnorm = 0.0
     do e = 1, em
        vnorm = vnorm + W(:,:,e) * V(:,:,e)**2
     end do
     where ( vnorm==0.0 ) ! division by zero protection
        y_d = 0.0
     elsewhere
        y_d = y_d / vnorm
     end where
     do e = 1, em
        V(:,:,e) = W(:,:,e) * V(:,:,e) * y_d(:,:)
     end do

#ifdef DEBUG
     call MAPL_MaxMin(' V ',V)
#endif

!    Gather V to all processes that will participate in the analysis
!    First we make a sub-communicator containing those processes (lde_comm)
!    ----------------------------------------------------------------------
     call ESMF_VMGet       (self%VM, mpiCommunicator=COMM, localpet=MYPE, petcount=NPES,  __RC__)

!    Next compute aerosol concentration analysis for each level, species
!             q_a = q_f + <X,V>
!    -------------------------------------------------------------------
     if ( MAPL_AM_I_Root() .and. verbose_ ) then
        if ( self%isCubed ) then
           print *, 'Calculating LDE increments on Cubed Sphere with ', em, ' ensemble members'
        else
           print *, 'Calculating LDE increments on LatLon with ', em, ' ensemble members'
        end if
     end if

     do s = 1, bQ_f%n3d
        if ( .not. isAerosol_(trim(bQ_f%r3(s)%name)) ) cycle

        if ( MAPL_AM_I_Root() .and. verbose_ ) &
           print *, ' [ ] Working on <'//trim(bQ_f%r3(s)%name)//'>'

#ifdef DEBUG
        call MAPL_MaxMin('      q_f',bQ_f%r3(s)%q(:,:,self%ks:km))
#endif

!       Each process does the analysis on its assigned level
!       ----------------------------------------------------
        do k = self%ks, km
           q_f => bQ_f%r3(s)%q(:,:,k)
           x_2d => x_d(:,:,k)
           call LDE_Qinc_Distrib_Cubed_(x_2d, q_f, V, self%indx, im, jm, em, IM_World, EM_World, nh, self, __RC__ ) 
        end do

!       Add analysis increments to q
!       ----------------------------
        do k = self%ks, km
           bQ_a%r3(s)%q(:,:,k) = bQ_f%r3(s)%q(:,:,k) + x_d(:,:,k)
        end do

!       Zero increments above top analysis level
!       ----------------------------------------
        do k = 1,self%ks-1 
           bQ_a%r3(s)%q(:,:,k) = bQ_f%r3(s)%q(:,:,k)
        end do
        
#ifdef DEBUG
        call MAPL_MaxMin('      q_a',bQ_a%r3(s)%q(:,:,self%ks:km))
#endif

        bQ_a%r3(s)%q(:,:,self%ks:km) = max(0.0,bQ_a%r3(s)%q(:,:,self%ks:km)) ! fix q<0

     end do ! variable loop

!    Free memory
!    -----------
     deallocate ( y_f, y_d, vnorm, X, V, W, x_d, __STAT__ )

     rc = 0
     return

   end subroutine LDE_Projector1c_Cubed_

!.........................................................................

   logical function isAerosol_ ( name )

     character(len=*), intent(in) :: name

                       __Iam__('isAerosol_')
         
     if ( ESMF_UtilStringUpperCase(name(1:2))=='DU'       .OR.  &
          ESMF_UtilStringUpperCase(name(1:2))=='SS'       .OR.  &
          ESMF_UtilStringUpperCase(name(1:2))=='NI'       .OR.  &
          ESMF_UtilStringUpperCase(name)     =='SO4'      .OR.  &
          ESMF_UtilStringUpperCase(name)     =='BCPHOBIC' .OR.  &
          ESMF_UtilStringUpperCase(name)     =='BCPHILIC' .OR.  &
          ESMF_UtilStringUpperCase(name)     =='OCPHOBIC' .OR.  &
          ESMF_UtilStringUpperCase(name)     =='OCPHILIC'       ) then
            
          isAerosol_ = .TRUE.

     else

          isAerosol_ = .FALSE.

     end if

   end function isAerosol_

 end module LDE_Mod
