#include "MAPL_Exceptions.h"
#include "MAPL_Generic.h"


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: MAM_NucleationMod - Nucleation (formation and growth) of new 
!                              aerosol particles
!
! !INTERFACE:
!
   module MAM_NucleationMod
!
! !USES:
!

   use ESMF

   use MAPL_Mod
   use MAPL_SimpleBundleMod

   use MAML_NucleationMod

   use MAM_ConstituentsDataMod
   use MAM_ComponentsDataMod

   use MAM3_DataMod 
   use MAM7_DataMod

   use MAM_BaseMod


   implicit NONE
   private

!
! !PUBLIC MEMBER FUNCTIONS:

   public MAM_Nucleation

!
! !PUBLIC PARAMETERS:

!
! !PRIVATE PARAMETERS:


!
! !DESCRIPTION: 
!
!  {\tt MAM\_NucleationMod} provides methods for computing the nucleation 
!  rates and changes of number and mass mixing ratios.
!
!
! !REVISION HISTORY:
!
!  27Jan2012  A. Darmenov  Initial version
!
!EOP
!-------------------------------------------------------------------------


   contains


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: MAM_Nucleation --- models the nucleation of Aitken mode 
!
!
! !INTERFACE:

   subroutine MAM_Nucleation(self, import, export, qa, qg, Da, cdt, rc)
! !USES:

   implicit NONE

! !INPUT/OUTPUT PARAMETERS:
   type(MAPL_SimpleBundle), intent(inout) :: qa         ! number/mass mixing ratio
   type(MAPL_SimpleBundle), intent(inout) :: qg         ! gas mixing ratio
   type(ESMF_State), intent(inout)        :: export     ! export state
   integer, optional, intent(inout)       :: rc         ! return code

! !INPUT PARAMETERS:
   type(MAPL_SimpleBundle), intent(in)    :: Da         ! dry(geometric mean) and wet diameter 
                                                        ! of number size distribution

   type(MAM_Scheme), intent(in)           :: self       ! MAM scheme/configuration

   type(ESMF_State), intent(inout)        :: import     ! import state

   real,    intent(in)                    :: cdt        ! chemical timestep (secs)

! !OUTPUT PARAMETERS:


! !DESCRIPTION: Nucleation: 
!               - MAM7 scheme models nucleation by treating ...
!               - MAM3 scheme ...
!
! !REVISION HISTORY:
!
!  09Jan2012  A. Darmenov  First crack.
!
!EOP
!-------------------------------------------------------------------------

              __Iam__('MAM_Nucleation')


   ! Mode parameters
   ! ---------------
   character(len=MAM_MAXSTR) :: mode_name

   integer                   :: n_species
   character(len=MAM_MAXSTR) :: species_name

   character(len=MAM_MAXSTR) :: field_name

   ! indexes 
   integer :: i, im, j, jm, k, km

   real    :: temperature                                               ! local temperature
   real    :: pressure                                                  ! local mid-level pressure
   real    :: density_air                                               ! local air density
   real    :: rel_humidity                                              ! local RH
   real    :: z                                                         ! local mid-level height
   real    :: f_cld                                                     ! local cloud fraction
   real    :: pbl_height                                                ! PBL height

   real    :: diameter_dry                                              ! local Aitken mode dry size 
   real    :: diameter_dry_min, diameter_dry_max                        ! Aitken mode dry size limits

   real    :: q_number                                                  ! local number mixing ratio
   real    :: q_amm                                                     ! local ammonium mass mixing ratio
   real    :: q_su                                                      ! local sulfate mass mixing ratio
   real    :: q_h2so4                                                   ! local H2SO4 volume mixing ratio
   real    :: q_nh3                                                     ! local NH3 volume mixing ratio

   real    :: dq_h2so4_gasprod                                          ! H2SO4 
   real    :: dq_h2so4_aeruptk                                          ! H2SO4

   integer :: ait_index  = 0                                            ! index of the Aitken mode
   integer :: iq_dgn_dry = 0                                            ! index of the Aitken mode dry diameter
   integer :: iq_number  = 0                                            ! index of the Aitken mode number mixing ratio
   integer :: iq_amm     = 0                                            ! index of the Aitken mode ammonium aerosol mass mixing ratio
   integer :: iq_su      = 0                                            ! index of the Aitken mode sulfate aerosol mass mixing ratio
   integer :: iq_h2so4   = 0                                            ! index of the H2SO4 volume mixing ratio
   integer :: iq_nh3     = 0                                            ! index of the NH3 volume mixing ratio

   logical :: do_nh3
 

   ! other derived variables
   real,    allocatable, dimension(:)   :: dz
   

   !  Input fields from fvGCM
   !  -----------------------
   real, pointer, dimension(:,:,:) :: rhoa
   real, pointer, dimension(:,:,:) :: delp
   real, pointer, dimension(:,:,:) :: ple
   real, pointer, dimension(:,:,:) :: zle
   real, pointer, dimension(:,:,:) :: T
   real, pointer, dimension(:,:,:) :: rh
   real, pointer, dimension(:,:,:) :: fcld
   real, pointer, dimension(:,:)   :: zpbl

   !  Exports
   !  -----------------------
   real, pointer, dimension(:,:)   :: flux

   !  Parameters
   !  ----------
   real, parameter :: mw_su  = MAM_SULFATE_COMPONENT_MOLECULAR_WEIGHT
   real, parameter :: mw_amm = MAM_AMMONIUM_COMPONENT_MOLECULAR_WEIGHT
   real, parameter :: density_su = MAM_SULFATE_COMPONENT_DENSITY



   !  Get Imports
   !  --------------
   call MAPL_GetPointer(import, rhoa, 'AIRDENS', __RC__)
   call MAPL_GetPointer(import, delp, 'DELP',    __RC__)
   call MAPL_GetPointer(import, ple,  'PLE',     __RC__)
   call MAPL_GetPointer(import, zle,  'ZLE',     __RC__)
   call MAPL_GetPointer(import, T,    'T',       __RC__)
   call MAPL_GetPointer(import, rh,   'RH2',     __RC__)
   call MAPL_GetPointer(import, fcld, 'FCLD',    __RC__)
   call MAPL_GetPointer(import, zpbl, 'ZPBL',    __RC__)

   !  Local dimensions
   !  ----------------
   im = size(rhoa, 1)
   jm = size(rhoa, 2)
   km = size(rhoa, 3)


#ifdef CAM
   !
   ! init ------------------------------- 
   ! 

   mait = modeptr_aitken

   if (mait > 0) then
        lnumait = numptr_amode(mait)
        lso4ait = lptr_so4_a_amode(mait)
        lnh4ait = lptr_nh4_a_amode(mait)
   end if

   if ((l_h2so4  <= 0) .or. (l_h2so4 > pcnst)) then
        write(*,'(/a/)')  '*** modal_aero_newnuc bypass -- l_h2so4 <= 0'
        return
   else if ((lso4ait <= 0) .or. (lso4ait > pcnst)) then
        write(*,'(/a/)')  '*** modal_aero_newnuc bypass -- lso4ait <= 0'
        return
   else if ((lnumait <= 0) .or. (lnumait > pcnst)) then
        write(*,'(/a/)')  '*** modal_aero_newnuc bypass -- lnumait <= 0'
        return
   else if ((mait <= 0) .or. (mait > ntot_amode)) then
        write(*,'(/a/)')  '*** modal_aero_newnuc bypass -- modeptr_aitken <= 0'
        return
   end if

   !
   ! run ------------------------------------------
   !

   !   skip if no aitken mode OR if no h2so4 species
   if ((l_h2so4 <= 0) .or. (lso4ait <= 0) .or. (lnumait <= 0)) then
       return
   end if

   lnh4ait = lptr_nh4_a_amode(mait) - loffset
   if ((l_nh3   > 0) .and. (l_nh3   <= pcnst) .and. &
        (lnh4ait > 0) .and. (lnh4ait <= pcnst)) then
        do_nh3 = .true.
        dotend(lnh4ait) = .true.
        dotend(l_nh3) = .true.
    else
        do_nh3 = .false.
    end if
#endif

   !  Find the nucleation mode
   !  ------------------------

   mode_name = ''

   if (self%id == MAM7_SCHEME) then
       mode_name = trim(MAM7_AITKEN_MODE_NAME)
   else if (self%id == MAM3_SCHEME) then
       mode_name = trim(MAM3_AITKEN_MODE_NAME)
   else
       return
   end if


   ait_index = 0
   ait_index = MAM_SchemeGetModeIndex(self, trim(mode_name), __RC__)

   call MAM_AerosolModeGet(self%mode(ait_index), size_min  = diameter_dry_min, &
                                                 size_max  = diameter_dry_max, &
                                                 n_species = n_species)


   !  find the indexes of the aerosol species 
   !  ---------------------------------------
   field_name = 'DGN_DRY_' // trim(mode_name)
   iq_dgn_dry = MAPL_SimpleBundleGetIndex(Da, field_name, 3, __RC__)

   field_name = trim(MAM_NUMBER_PARTICLES_NAME) // '_A_' // trim(mode_name)
   iq_number = MAPL_SimpleBundleGetIndex(qa, field_name, 3, __RC__)

   field_name = trim(MAM_AMMONIUM_CONSTITUENT_NAME) // '_A_' // trim(mode_name)
   iq_amm = MAPL_SimpleBundleGetIndex(qa, field_name, 3, __RC__)

   field_name = trim(MAM_SULFATE_CONSTITUENT_NAME) // '_A_' // trim(mode_name)
   iq_su = MAPL_SimpleBundleGetIndex(qa, field_name, 3, __RC__)

   
   ! find the indexes of the gas species
   ! -----------------------------------
   field_name = trim(MAM_H2SO4_CONSTITUENT_NAME)
   iq_h2so4 = MAPL_SimpleBundleGetIndex(qg, field_name, 3, __RC__)

   field_name = trim(MAM_NH3_CONSTITUENT_NAME)
   iq_nh3   = MAPL_SimpleBundleGetIndex(qg, field_name, 3, __RC__)



   ! nucleation calculation and column integrated number flux due to it
   allocate(dz(km), __STAT__)

   do j = 1, jm
       do i = 1, im

           do k = 1, km

               ! mid level pressure, temperature and air density
               pressure = 0.5 * (ple(i,j,k-1) + ple(i,j,k))
               z        = 0.5 * (zle(i,j,k-1) + zle(i,j,k))

               temperature  = T(i,j,k)
               density_air  = rhoa(i,j,k)
               rel_humidity = rh(i,j,k)
               pbl_height   = zpbl(i,j)

               ! aerosol dry size
               diameter_dry     = Da%r3(iq_dgn_dry)%q(i,j,k)

               ! aerosol species mass mixing ratios
               q_number = qa%r3(iq_number)%q(i,j,k)
               q_amm    = qa%r3(iq_amm)%q(i,j,k)
               q_su     = qa%r3(iq_su)%q(i,j,k)

               q_h2so4  = qg%r3(iq_h2so4)%q(i,j,k)
               q_nh3    = qg%r3(iq_nh3)%q(i,j,k)

               dq_h2so4_gasprod = 0.0
               dq_h2so4_aeruptk = 0.0

               do_nh3 = .true.

               ! nucleation
               call MAML_Nucleation(pressure,         &
                                    temperature,      &
                                    density_air,      &
                                    rel_humidity,     &
                                    f_cld,            &
                                    z,                &
                                    pbl_height,       &
                                    q_number,         &
                                    q_amm,            &
                                    q_su,             &
                                    q_h2so4,          &
                                    q_nh3,            &
                                    do_nh3,           &
                                    diameter_dry,     &
                                    diameter_dry_min, &
                                    diameter_dry_max, &
                                    density_su,       &
                                    mw_su,            &
                                    mw_amm,           &
                                    dq_h2so4_gasprod, &
                                    dq_h2so4_aeruptk, &
                                    cdt)


               ! update the number and mass mixing ratios
               qa%r3(iq_number)%q(i,j,k) = q_number

               qa%r3(iq_amm)%q(i,j,k)    = q_amm
               qa%r3(iq_su)%q(i,j,k)     = q_su

               qg%r3(iq_h2so4)%q(i,j,k)  = q_h2so4
               qg%r3(iq_nh3)%q(i,j,k)    = q_nh3

           end do ! k

           dz(:) = delp(i,j,:) / (MAPL_GRAV * rhoa(i,j,:))

       end do ! i
   end do ! j

   deallocate(dz, __STAT__)

   RETURN_(ESMF_SUCCESS)

   end subroutine MAM_Nucleation


end module MAM_NucleationMod
