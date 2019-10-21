#include "MAPL_Exceptions.h"
#include "MAPL_Generic.h"


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: MAM_WetRemovalMod - Scavenging in convective updrafts and first-order
!                              rainout and washout by precipitation
!
! !INTERFACE:
!
   module MAM_WetRemovalMod
!
! !USES:
!
   use ESMF

   use MAPL_Mod

   use MAPL_SimpleBundleMod

   use MAML_WetRemovalMod
   use MAM_BaseMod

   implicit NONE
   private

!
! !PUBLIC MEMBER FUNCTIONS:

   public MAM_WetRemoval

!
! !PUBLIC PARAMETERS:

!
! !PRIVATE PARAMETERS:
   real, private, parameter :: pi = MAPL_PI
   real, private, parameter :: density_water = MAPL_RHOWTR    ! density of water,  'kg m-3'


!
! !DESCRIPTION: 
!
!  {\tt MAML\_WetRemovalMod} provides a collection of methods for 
!  modeling wet removal of aerosol particles.
!
!
! !REVISION HISTORY:
!
!  29Nov2012  A. Darmenov  Initial version
!
!EOP
!-------------------------------------------------------------------------


   contains

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: MAM_WetRemoval --- 
!
! !INTERFACE:

   subroutine MAM_WetRemoval(self, import, export, qa, cdt, rc)
! !USES:

   implicit NONE

! !INPUT/OUTPUT PARAMETERS:
   type(MAPL_SimpleBundle), intent(inout) :: qa         ! number/mass mixing ratio
   type(ESMF_State), intent(inout)        :: export     ! export state
   integer, optional, intent(inout)       :: rc         ! return code

! !INPUT PARAMETERS:
   type(MAM_Scheme), intent(in)           :: self       ! MAM scheme/configuration

   type(ESMF_State), intent(inout)        :: import     ! import state

   real,    intent(in)                    :: cdt        ! chemical timestep (secs)

! !OUTPUT PARAMETERS:


! !DESCRIPTION: Wet removal of aerosol particles.
!
! !REVISION HISTORY:
!
! 29Nov2011  A. Darmenov  First crack.
!
!EOP
!-------------------------------------------------------------------------

              __Iam__('MAM_WetRemoval')

   ! Mode parameters
   ! ---------------
   real                            :: f_wet
   character(len=MAM_MAXSTR)       :: aero_type
   logical                         :: kin
   character(len=MAM_MAXSTR)       :: mode_name

   integer                         :: n_species
   character(len=MAM_MAXSTR)       :: species_name

   character(len=MAM_MAXSTR)       :: field_name

   integer                         :: m, s
   integer                         :: i_qa


   !  Input fields from fvGCM
   !  -----------------------
   real, pointer, dimension(:,:,:) :: rhoa
   real, pointer, dimension(:,:,:) :: ple
   real, pointer, dimension(:,:,:) :: T
  
   real, pointer, dimension(:,:,:) :: pfllsan
   real, pointer, dimension(:,:,:) :: pfilsan
   real, pointer, dimension(:,:,:) :: qccu
   real, pointer, dimension(:,:,:) :: cmfmc
   real, pointer, dimension(:,:,:) :: dtrain

   real, pointer, dimension(:,:)   :: precc
   real, pointer, dimension(:,:)   :: precl

   !  Exports
   !  -----------------------
   real, pointer, dimension(:,:)   :: flux => null()
   

   !  Get Imports
   !  --------------
   call MAPL_GetPointer(import, rhoa,    'AIRDENS',  __RC__)
   call MAPL_GetPointer(import, ple,     'PLE',      __RC__)
   call MAPL_GetPointer(import, T,       'T',        __RC__)
   call MAPL_GetPointer(import, precc,   'CN_PRCP',  __RC__)
   call MAPL_GetPointer(import, precl,   'NCN_PRCP', __RC__)
   call MAPL_GetPointer(import, pfllsan, 'PFL_LSAN', __RC__)
   call MAPL_GetPointer(import, pfilsan, 'PFI_LSAN', __RC__)
   call MAPL_GetPointer(import, cmfmc,   'CNV_MFC',  __RC__)
   call MAPL_GetPointer(import, dtrain,  'CNV_MFD',  __RC__)

   ! large-scale wet removal
   aero_type = 'sulfate'
   kin       = .true.
   
   do m = 1, self%n_modes
       call MAM_AerosolModeGet(self%mode(m), name      = mode_name, &
                                             n_species = n_species, &
                                             f_wet     = f_wet)

       ! number mixing ratio
       field_name = 'NUM_A_' // trim(mode_name)
       i_qa = MAPL_SimpleBundleGetIndex(qa, field_name, 3, __RC__)

       call MAML_WetRemoval(qa%r3(i_qa)%q, f_wet, aero_type, kin, ple, T, rhoa, &
                            pfllsan, pfilsan, precc, precl, cdt, flux, __RC__) 

 
       do s = 1, n_species
           species_name = self%mode(m)%species(s)%name
           field_name  = trim(species_name) // '_A_' // trim(mode_name)
           i_qa      = MAPL_SimpleBundleGetIndex(qa, field_name, 3, __RC__)

           call MAML_WetRemoval(qa%r3(i_qa)%q, f_wet, aero_type, kin, ple, T, rhoa, &
                                pfllsan, pfilsan, precc, precl, cdt, flux, __RC__) 
       end do       
   end do


   RETURN_(ESMF_SUCCESS)

   end subroutine MAM_WetRemoval


end module MAM_WetRemovalMod
