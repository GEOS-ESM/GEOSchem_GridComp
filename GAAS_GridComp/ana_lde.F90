! Ana_LDE: Produces 3D Aerosol Anaysis Increments using Lagrangian Displacent
!          Ensembles and pre-computed AOD analysis increments.
!
! Arlindo da Silva <arlindo.dasilva@nasa.gov>, December 2009
!----------------------------------------------------------------------------

#  include "MAPL_Generic.h"

!............................................................................................

   Program Ana_LDE

   use ESMF
   use MAPL

   use LDE_Mod
   use Chem_SimpleBundleMod

   implicit NONE

!  Basic ESMF objects
!  ------------------
   type(ESMF_Config)       :: CF       ! Resource file
   type(ESMF_Grid)         :: etaGrid  ! Eta Grid (lon, lat, eta)
   type(ESMF_Grid)         :: aodGrid  ! AOD Grid (lon, lat, channel)
   type(ESMF_Time)         :: Time    

!  Main data objects
!  -----------------
   type(MAPL_SimpleBundle), target  :: q_f  ! 3D backround aerosol concentration
   type(MAPL_SimpleBundle), target  :: q_b  ! 3D backround aerosol bias
   type(MAPL_SimpleBundle), pointer :: q_a  ! 3D analyzed  aerosol concentration

   type(MAPL_SimpleBundle) :: y_f  ! 2D AOD background (log-transform)
   type(MAPL_SimpleBundle) :: y_d  ! 2D AOD analysis increment (log-transform)
   type(MAPL_SimpleBundle) :: y_b  ! 2D AOD bias (log-transform)

   type(LDE)               :: E    ! LDE object

!  Basic information about the parallel environment
!         PET = Persistent Execution Threads
!  In the current implementation, a PET is equivalent 
!  to an MPI process
!  ------------------------------------------------
   integer :: myPET   ! The local PET number
   integer :: nPET    ! The total number of PETs you are running on

   integer :: status, rc
   integer :: i, j, n, im, jm, km

   integer :: Nx, Ny                           ! Layout
   integer :: IM_World, JM_World, LM_WORLD     ! Global Grid dimensions
   integer :: CM_World                         ! Number of channels
   integer :: nymd, nhms
   integer :: yy, mm, dd, h, m, s, iAOD
   real    :: eps  ! eps for log-transform
   real    :: channel

   logical :: verbose, perform_lde_in_aod_variable, isCubed
 
   type (CubedSphereGridFactory) :: cs_factory
   type (LatlonGridFactory) :: ll_factory

!  Coordinate variables
!  --------------------
   character(len=ESMF_MAXSTR)    :: name
   character(len=ESMF_MAXSTR)    :: Iam = "ana_lde"

!                             -----

    call Main()

CONTAINS

    subroutine Main()

!   For now concentration analysis/background share storage
!   -------------------------------------------------------
    q_a => q_f  

!                                --------------
!                                Initialization
!                                --------------


!   Initialize the ESMF. For performance reasons, it is important
!    to turn OFF ESMF's automatic logging feature
!   -------------------------------------------------------------
    call ESMF_Initialize (LogKindFlag=ESMF_LOGKIND_NONE, __RC__)
    call ESMF_CalendarSetDefault ( ESMF_CALKIND_GREGORIAN, __RC__ )

    if ( MAPL_am_I_root() ) then
         print *
         print *, '        ----------------'
         print *, '        Starting ' // trim(Iam)
         print *, '        ----------------'
         print *
    end if

!   Load resources
!   --------------
    CF = ESMF_ConfigCreate(__RC__)
    call ESMF_ConfigLoadFile(CF, fileName='lde.rc', __RC__)
    call ESMF_ConfigGetAttribute(CF, verbose, Label='verbose:',  __RC__ )
    call ESMF_ConfigGetAttribute(CF, eps, Label='eps_for_log_transform_aod:', &
                                default=-0.01, __RC__)
    call ESMF_ConfigGetAttribute(CF, perform_lde_in_aod_variable, &
                                 Label='perform_lde_in_aod_variable:', __RC__)


!   World grid dimensions and layout
!   --------------------------------
    call ESMF_ConfigGetAttribute(CF, IM_World, Label='IM_World:',  __RC__ )
    call ESMF_ConfigGetAttribute(CF, JM_World, Label='JM_World:',  __RC__ )
    call ESMF_ConfigGetAttribute(CF, LM_World, Label='LM_World:',  __RC__ )
    call ESMF_ConfigGetAttribute(CF, CM_World, Label='CM_World:',  __RC__ )
    call ESMF_ConfigGetAttribute(CF, Nx,       Label='Layout_Nx:', __RC__ )
    call ESMF_ConfigGetAttribute(CF, Ny,       Label='Layout_Ny:', __RC__ )

!   Lat lon or cubed sphere?
!   ------------------------
    if ( JM_World == 6*IM_World ) then
       isCubed = .True.
    else
       isCubed = .False.
    end if

!   Create global grids:
!   -------------------
    if ( isCubed ) then ! Cubed-sphere

       cs_factory = CubedSphereGridFactory(im_world=im_world,lm=lm_world,nx=nx,ny=ny/6,__RC__)
       etaGrid = cs_factory%make_grid(__RC__)
       cs_factory = CubedSphereGridFactory(im_world=im_world,lm=cm_world,nx=nx,ny=ny/6,__RC__)
       aodGrid = cs_factory%make_grid(__RC__)

    else ! Lat Lon Grid

       ll_factory = LatLonGridFactory(grid_name='etaGrid', &
                        Nx = Nx, Ny = Ny,   &
                        IM_World = IM_World,  &
                        JM_World = JM_World,  &
                        LM = LM_WORLD, pole='PC', dateline='DC',  &
                               __RC__)
       etaGrid = ll_factory%make_grid(__RC__)
       ll_factory = LatLonGridFactory(grid_name='etaGrid', &
                        Nx = Nx, Ny = Ny,   &
                        IM_World = IM_World,  &
                        JM_World = JM_World,  &
                        LM = CM_WORLD, pole='PC', dateline='DC',  &
                               __RC__)
       aodGrid = ll_factory%make_grid(__RC__)

    end if

!   Validate grid
!   -------------
    call ESMF_GridValidate(etaGrid,__RC__)
    call ESMF_GridValidate(aodGrid,__RC__)

!   Get date/time from CF
!   ---------------------
    call ESMF_ConfigGetAttribute(CF, nymd, Label='nymd:', __RC__ )
    call ESMF_ConfigGetAttribute(CF, nhms, Label='nhms:', __RC__ )

     call ESMF_ConfigGetAttribute(CF, channel, Label='single_channel:',  __RC__) 

!   Create ESMF Time
!   ----------------
    yy = nymd/10000; mm = (nymd-yy*10000) / 100; dd = nymd - (10000*yy + mm*100)
    h  = nhms/10000;  m = (nhms - h*10000) / 100;  s = nhms - (10000*h  +  m*100)
    call ESMF_TimeSet(Time, yy=yy, mm=mm, dd=dd,  h=h,  m=m, s=s)

!                                --------------
!                                  Read Files
!                                --------------

    y_f = Chem_SimpleBundleRead (CF, 'aod_bkg_filename',         aodGrid, time=Time, __RC__ )
    y_d = Chem_SimpleBundleRead (CF, 'aod_inc_filename',         aodGrid, time=Time, __RC__ )
    q_f = Chem_SimpleBundleRead (CF, 'aer_bkg_filename',         etaGrid, time=Time, verbose=.TRUE., __RC__ )
#if 0
    q_b = Chem_SimpleBundleRead (CF, 'aerbias_internal_restart', etaGrid, __RC__ )
    y_b = Chem_SimpleBundleRead (CF, 'aodbias_internal_restart', aodGrid, __RC__ )
#endif

    y_f%coords%levs(1) = channel
    y_d%coords%levs(1) = channel

    if ( verbose ) then
       call MAPL_SimpleBundlePrint(y_f)
       call MAPL_SimpleBundlePrint(y_d)
       call MAPL_SimpleBundlePrint(q_f)
    endif

!   Here we assume that the analysis increments coming in are log(AOD+eps)
!    while y_a and y_f are always in terms of AOD. When doing the LDE the
!    *perform_lde_in_aod_variable* on the rc files specifies whether we 
!    generate the LDE ensembles based on AOD or Log(AOD+eps), regardless 
!    of how the analysis was performed.
!
!    TO DO: 
!        Figure out an automatic  way to check whether analysis in log-AOD
!   -----------------------------------------------------------------------
     if ( perform_lde_in_aod_variable ) then ! use AOD for LDE calculation
          iAOD = MAPL_SimpleBundleGetIndex(y_d,'AOD',3,__RC__)
                                                                           !Note: eta = log(tau+eps)
          y_d%r3(iAOD)%q(:,:,:) = Exp(                                   &
                                       Log(y_f%r3(iAOD)%q(:,:,:) + eps)  & ! eta_f
                                       + y_d%r3(iAOD)%q(:,:,:)           & ! + deta_a
                                     ) - eps                             & ! tau_a = exp(eta_a) - eps
                                  - y_f%r3(iAOD)%q(:,:,:)                  ! - tau_e

     else ! use log(AOD+eps) for LDE calculation
          iAOD = MAPL_SimpleBundleGetIndex(y_f,'AOD',3,__RC__)
          y_f%r3(iAOD)%q(:,:,:) = Log(y_f%r3(iAOD)%q(:,:,:) + eps)
     end if
          
!                                 -----------
!                                 Calculation
!                                 -----------

!   Create concetration analysis from AOD increments
!   ------------------------------------------------
    call LDE_Create ( E, CF, aodGrid, __RC__ )
    call LDE_Projector1c ( E, q_f, q_f, y_f, y_d, verbose, __RC__ )
    call LDE_Destroy ( E, __RC__ )

!                                --------------
!                                 Write Files
!                                --------------

    call Chem_SimpleBundleWrite (q_a, CF, 'aer_ana_filename',            Time, __RC__ )
#if 0
    call Chem_SimpleBundleWrite (q_b, CF, 'aerbias_internal_checkpoint', Time, __RC__ )
    call Chem_SimpleBundleWrite (y_b, CF, 'aodbias_internal_checkpoint', Time, __RC__ )
#endif

    if ( verbose ) then
       call MAPL_SimpleBundlePrint(q_a)
    endif



!   All done
!   --------
    call ESMF_Finalize(__RC__)

  end subroutine Main

  end Program Ana_LDE

