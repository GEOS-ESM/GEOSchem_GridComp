#include "MAPL_Exceptions.h"
#include "MAPL_Generic.h"


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: MAM_gasAerosolExchange - Gas condensation on aerosol particles 
!
!
! !INTERFACE:
!
   module MAM_GasAerosolExchangeMod
!
! !USES:
!

   use ESMF

   use MAPL_Mod
   use MAPL_SimpleBundleMod

!  use MAML_GasAerosolExchangeMod

   use MAM_ConstituentsDataMod
   use MAM_ComponentsDataMod

   use MAM3_DataMod
   use MAM7_DataMod

   use MAM_BaseMod


   implicit NONE
   private

!
! !PUBLIC MEMBER FUNCTIONS:

   public MAM_GasAerosolExchange

!
! !PUBLIC PARAMETERS:

!
! !PRIVATE PARAMETERS:


!
! !DESCRIPTION: 
!
!  {\tt MAM\_GasAerosolEchangeMod} provides methods for computing the ... 
!  rates and changes of number and mass mixing ratios.
!
!
! !REVISION HISTORY:
!
!  29Jun2012  A. Darmenov  Initial version
!
!EOP
!-------------------------------------------------------------------------


   contains


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: MAM_GasAerosolExchange --- models ...
!
!
! !INTERFACE:

   subroutine MAM_GasAerosolExchange(self, qa, qc, qg, Da, import, export, cdt, rc)
! !USES:

   implicit NONE

! !INPUT/OUTPUT PARAMETERS:
   type(MAPL_SimpleBundle), intent(inout) :: qa         ! number/mass mixing ratio of interstitial aerosols
   type(MAPL_SimpleBundle), intent(inout) :: qc         ! number/mass mixing ratio of cloud-borne  aerosols
   type(MAPL_SimpleBundle), intent(inout) :: qg         ! mixing ratio of gas species
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
!               - MAM7 and MAM3 schemes model gas aerosol exchange ...
!
! !REVISION HISTORY:
!
!  29Jun2012  A. Darmenov  First crack.
!
!EOP
!-------------------------------------------------------------------------

              __Iam__('MAM_GasAerosolExchange')
#ifdef WORK_IN_PROGRESS

   ! Mode parameters
   ! ---------------
   character(len=MAM_MAXSTR) :: species_name

   character(len=MAM_MAXSTR) :: field_name

   ! indexes 
   integer :: i, im, j, jm, k, km
   integer :: s, m

   real    :: temperature                               ! local temperature
   real    :: pressure                                  ! local mid-level pressure
   real    :: density_air                               ! local air density
   real    :: rel_humidity                              ! local RH
   real    :: z                                         ! local mid-level height
   real    :: f_cld                                     ! local cloud fraction

   real    :: qg_h2so4                                  ! H2SO4    volume mixing ratio
   real    :: qg_nh3                                    ! NH3      volume mixing ratio
   real    :: qg_soag                                   ! SOA(gas) volume mixing ratio

   integer :: n_modes                                   ! number of aerosol modes
   integer :: n_species                                 ! number of species in a mode
   integer :: n_max_species                             ! largest number of species

   integer :: iq_h2so4                                  ! index of the H2SO4 volume mixing ratio
   integer :: iq_nh3                                    ! index of the NH3 volume mixing ratio
   integer :: iq_soag                                   ! index of the SOA(gas) volume mixing ratio


!! integer, dimension(:,:), allocatable :: id_species   ! IDs of species in aerosol modes
!! integer, dimension(:),   allocatable :: id_mode      ! IDs of modes
!!
!! real,    dimension(:),   allocatable :: qa_number    ! interstitial aerosols - number mixing ratio
!! real,    dimension(:,:), allocatable :: qa_mass      !                       - mass mixing ratio
!!
!! real,    dimension(:),   allocatable :: qc_number    ! cloud-borne  aerosols - number mixing ratio
!! real,    dimension(:,:), allocatable :: qc_mass      !                         mass mixing ratio
!!
!! real,    dimension(:),   allocatable :: D_dry        ! dry size (geometric mean diameter of number size distribution)
!! real,    dimension(:),   allocatable :: D_wet        ! wet size

   
#if (0)
   real    :: dq_h2so4_gasprod                          ! H2SO4 
   real    :: dq_h2so4_aeruptk                          ! H2SO4

   integer :: ait_index                                 ! index of the Aitken mode
   integer :: acc_index                                 ! index of the accumulation mode

   integer :: iq_dgn_dry = 0                            ! index of the Aitken mode dry diameter
   integer :: iq_number  = 0                            ! index of the Aitken mode number mixing ratio
   integer :: iq_amm     = 0                            ! index of the Aitken mode ammonium aerosol mass mixing ratio
   integer :: iq_su      = 0                            ! index of the Aitken mode sulfate aerosol mass mixing ratio
   integer :: iq_h2so4   = 0                            ! index of the H2SO4 volume mixing ratio
   integer :: iq_nh3     = 0                            ! index of the NH3 volume mixing ratio

   logical :: do_nh3                                    !

   ! other derived variables
   real,    allocatable, dimension(:)   :: dz
#endif



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
#if (0)   
   real, parameter :: mw_su    = MAM_SULFATE_COMPONENT_MOLECULAR_WEIGHT
   real, parameter :: mw_amm   = MAM_AMMONIUM_COMPONENT_MOLECULAR_WEIGHT
   real, parameter :: mw_soa   = MAM_SOA_COMPONENT_MOLECULAR_WEIGHT

   real, parameter :: mw_h2so4 = MAM_H2SO4_COMPONENT_MOLECULAR_WEIGHT
   real, parameter :: mw_nh3   = MAM_NH3_COMPONENT_MOLECULAR_WEIGHT
   real, parameter :: mw_soag  = MAM_SOA_GAS_COMPONENT_MOLECULAR_WEIGHT
#endif




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
   !

   ! define "from mode" and "to mode" for primary carbon aging
   !
   ! skip (turn off) aging if either is absent, 
   ! or if accum mode so4 is absent
   !
   modefrm_pcage = -999888777
   modetoo_pcage = -999888777

   if ((modeptr_pcarbon <= 0) .or. (modeptr_accum <= 0)) goto 15000

   l = lptr_so4_a_amode(modeptr_accum)
   if ((l < 1) .or. (l > pcnst)) goto 15000

   modefrm_pcage = modeptr_pcarbon
   modetoo_pcage = modeptr_accum

   !
   ! define species involved in each primary carbon aging pairing
   !	(include aerosol water)
   !   
   !
   mfrm = modefrm_pcage
   mtoo = modetoo_pcage

   nspec = 0
   aa_iqfrm: do iqfrm = -1, nspec_amode(mfrm)

       if (iqfrm == -1) then
           lsfrm = numptr_amode(mfrm)
           lstoo = numptr_amode(mtoo)
       else if (iqfrm == 0) then
           ! bypass transfer of aerosol water due to primary-carbon aging
           cycle aa_iqfrm
           ! lsfrm = lwaterptr_amode(mfrm)
           ! lstoo = lwaterptr_amode(mtoo)
       else
           lsfrm = lmassptr_amode(iqfrm,mfrm)
           lstoo = 0
       end if

       if ((lsfrm < 1) .or. (lsfrm > pcnst)) cycle aa_iqfrm

       if (lsfrm>0 .and. iqfrm>0 ) then
           ! find "too" species having same lspectype_amode as the "frm" species
           do iqtoo = 1, nspec_amode(mtoo)
               if ( lspectype_amode(iqtoo,mtoo) .eq. lspectype_amode(iqfrm,mfrm) ) then
                   lstoo = lmassptr_amode(iqtoo,mtoo)
                   exit
               end if
           end do
       end if

       if ((lstoo < 1) .or. (lstoo > pcnst)) lstoo = 0
       nspec = nspec + 1
       lspecfrm_pcage(nspec) = lsfrm
       lspectoo_pcage(nspec) = lstoo
   end do aa_iqfrm

   nspecfrm_pcage = nspec


15000 continue

   ! set gas species indices
   call cnst_get_ind( 'H2SO4', l_so4g, .false. )
   call cnst_get_ind( 'NH3',   l_nh4g, .false. )
   call cnst_get_ind( 'MSA',   l_msag, .false. )
   call cnst_get_ind( 'SOAG',  l_soag, .false. )

   if ((l_so4g <= 0) .or. (l_so4g > pcnst)) then
       write( *, '(/a/a,2i7)' )   &
            '*** modal_aero_gasaerexch_init -- cannot find H2SO4 species',   &
            '    l_so4g=', l_so4g
      call endrun( 'modal_aero_gasaerexch_init error' )
   end if

   do_nh4g = .false.
   do_msag = .false.
   do_soag = .false.

   if ((l_nh4g > 0) .and. (l_nh4g <= pcnst)) do_nh4g = .true.
   if ((l_msag > 0) .and. (l_msag <= pcnst)) do_msag = .true.
   if ((l_soag > 0) .and. (l_soag <= pcnst)) do_soag = .true.


   ! set tendency flags
   dotend(:) = .false.
   dotend(l_so4g) = .true.

   if ( do_nh4g ) dotend(l_nh4g) = .true.
   if ( do_msag ) dotend(l_msag) = .true.
   if ( do_soag ) dotend(l_soag) = .true.

   do n = 1, ntot_amode
       l = lptr_so4_a_amode(n)
       if ((l > 0) .and. (l <= pcnst)) then
           dotend(l) = .true.

           if ( do_nh4g ) then
               l = lptr_nh4_a_amode(n)

               if ((l > 0) .and. (l <= pcnst)) dotend(l) = .true.
           end if
       end if

       l = lptr_soa_a_amode(n)
       if ((l > 0) .and. (l <= pcnst)) then
           dotend(l) = .true.
       end if
   end do

   if (modefrm_pcage > 0) then
       do iq = 1, nspecfrm_pcage
           lsfrm = lspecfrm_pcage(iq)
           lstoo = lspectoo_pcage(iq)

           if ((lsfrm > 0) .and. (lsfrm <= pcnst)) then
               dotend(lsfrm) = .true.

               if ((lstoo > 0) .and. (lstoo <= pcnst)) then
                   dotend(lstoo) = .true.
               end if

           end if
       end do
   end if


   !  define history fields for aitken-->accum renaming
   dotend(:) = .false.
   dotendqqcw(:) = .false.

   do ipair = 1, npair_renamexf
       do iq = 1, nspecfrm_renamexf(ipair)

           lsfrm = lspecfrma_renamexf(iq,ipair)
           lstoo = lspectooa_renamexf(iq,ipair)

           if ((lsfrm > 0) .and. (lsfrm <= pcnst)) then
               dotend(lsfrm) = .true.

               if ((lstoo > 0) .and. (lstoo <= pcnst)) then
                   dotend(lstoo) = .true.
               end if
           end if

           lsfrm = lspecfrmc_renamexf(iq,ipair)
           lstoo = lspectooc_renamexf(iq,ipair)

           if ((lsfrm > 0) .and. (lsfrm <= pcnst)) then
               dotendqqcw(lsfrm) = .true.
               if ((lstoo > 0) .and. (lstoo <= pcnst)) then
                   dotendqqcw(lstoo) = .true.
               end if
           end if

       end do ! iq = ...
   end do ! ipair = ...


   ! calculate soa_equivso4_factor
   ! if do_soag == .false., then set it to zero as a safety measure
   soa_equivso4_factor = 0.0

   if ( do_soag ) then
       tmp1 = -1.0 ; tmp2 = -1.0

       do l = 1, ntot_aspectype
           if (specname_amode(l) == 's-organic') tmp1 = spechygro(l)
           if (specname_amode(l) == 'sulfate'  ) tmp2 = spechygro(l)
       end do

       if ((tmp1 > 0.0_r8) .and. (tmp2 > 0.0_r8)) then
           soa_equivso4_factor = tmp1/tmp2
       else
           write(*,'(a/a,1p,2e10.2)') '*** subr modal_aero_gasaerexch_init', &
                                      '    cannot find hygros - tmp1/2 =', tmp1, tmp2
           call endrun()
       end if
   end if


   !
   ! run ------------------------------------------
   !

   call cnst_get_ind( 'H2SO4', l_so4g, .false. )
   call cnst_get_ind( 'NH3',   l_nh4g, .false. )
   call cnst_get_ind( 'MSA',   l_msag, .false. )
   call cnst_get_ind( 'SOAG',  l_soag, .false. )

   l_so4g = l_so4g - loffset
   l_nh4g = l_nh4g - loffset
   l_msag = l_msag - loffset
   l_soag = l_soag - loffset

   if ((l_so4g <= 0) .or. (l_so4g > pcnstxx)) then
       write( *, '(/a/a,2i7)' )   '*** modal_aero_gasaerexch_sub -- cannot find H2SO4 species',   &
                                  '    l_so4g, loffset =', l_so4g, loffset
       call endrun( 'modal_aero_gasaerexch_sub error' )
   end if

   do_nh4g = .false.
   do_msag = .false.

   if ((l_nh4g > 0) .and. (l_nh4g <= pcnstxx)) do_nh4g = .true.
   if ((l_msag > 0) .and. (l_msag <= pcnstxx)) do_msag = .true.

   do_soag = .false.

   if ((method_soa == 1) .or. (method_soa == 2)) then
       if ((l_soag > 0) .and. (l_soag <= pcnstxx)) do_soag = .true.
   else if (method_soa /= 0) then
       write(*,'(/a,1x,i10)') '*** modal_aero_gasaerexch_sub - bad method_soa =', method_soa
       call endrun( 'modal_aero_gasaerexch_sub error' )
   end if

   ! set tendency flags
   dotend(:) = .false.
   dotendqqcw(:) = .false.
   ido_so4a(:) = 0
   ido_nh4a(:) = 0
   ido_soaa(:) = 0

   dotend(l_so4g) = .true.
   if ( do_nh4g ) dotend(l_nh4g) = .true.
   if ( do_msag ) dotend(l_msag) = .true.
   if ( do_soag ) dotend(l_soag) = .true.

   ntot_soamode = 0
   do n = 1, ntot_amode

       l = lptr_so4_a_amode(n)-loffset

       if ((l > 0) .and. (l <= pcnstxx)) then

           dotend(l) = .true.
           ido_so4a(n) = 1

           if ( do_nh4g ) then
               l = lptr_nh4_a_amode(n)-loffset
               if ((l > 0) .and. (l <= pcnstxx)) then
                   dotend(l) = .true.
                   ido_nh4a(n) = 1
               end if
           end if
       end if

       if ( do_soag ) then

           l = lptr_soa_a_amode(n)-loffset

           if ((l > 0) .and. (l <= pcnstxx)) then
               dotend(l) = .true.
               ido_soaa(n) = 1
               ntot_soamode = n
           end if
       end if
   end do

   if ( do_soag ) ntot_soamode = max( ntot_soamode, modefrm_pcage )

   if (modefrm_pcage > 0) then
       ido_so4a(modefrm_pcage) = 2

       if (ido_nh4a(modetoo_pcage) == 1) ido_nh4a(modefrm_pcage) = 2
       if (ido_soaa(modetoo_pcage) == 1) ido_soaa(modefrm_pcage) = 2

       do iq = 1, nspecfrm_pcage
           lsfrm = lspecfrm_pcage(iq)-loffset
           lstoo = lspectoo_pcage(iq)-loffset

           if ((lsfrm > 0) .and. (lsfrm <= pcnst)) then

               dotend(lsfrm) = .true.

               if ((lstoo > 0) .and. (lstoo <= pcnst)) then
                   dotend(lstoo) = .true.
               end if
           end if
       end do
   end if
#endif

   ! Find the number of modes 
   ! ------------------------
   n_modes = self%n_modes

   
   ! Find the largest number of aerosol species in a single mode
   ! -----------------------------------------------------------
   n_max_species = 0
   do m = 1, n_modes
       call MAM_AerosolModeGet(self%mode(m), n_species=n_species)
       n_max_species = max(n_max_species, n_species)
   end do

   
   ! Allocate memory for the buffers
!! allocate(qa_number(n_modes),  __STAT__)               ! number mixing ratio of interstitial aerosols
!! allocate(qc_number(n_modes),  __STAT__)               ! number mixing ratio of cloud-borne  aerosols   
!!
!! allocate(D_dry(n_modes), __STAT__)                    ! dry size - geometric mean diameter of number size distribution
!! allocate(D_wet(n_modes), __STAT__)                    ! wet size
!!
!! allocate(qa_mass(n_max_species,n_modes), __STAT__)    !
!! allocate(qc_mass(n_max_species,n_modes), __STAT__)    !
!!
!! allocate(id_species(n_max_species,n_modes), __STAT__) !
!! allocate(id_mode(n_modes), __STAT__)

#if (0)
   !  Find aitken and accumulation modes
   !  ----------------------------------

   mode_name_ait = ''
   mode_name_acc = ''

   if (self%id == MAM7_MODEL) then
       mode_name_ait = trim(MAM7_AITKEN_MODE_NAME)
       mode_name_acc = trim(MAM7_ACCUMULATION_MODE_NAME) 
   else if (self%id == MAM3_MODEL) then
       mode_name_ait = trim(MAM3_AITKEN_MODE_NAME)
       mode_name_acc = trim(MAM3_ACCUMULATION_MODE_NAME)
   else
       return
   end if
#endif
   


   ! Find the indexes of the gas species
   ! -----------------------------------
   field_name = trim(MAM_H2SO4_CONSTITUENT_NAME)
   iq_h2so4   = MAPL_SimpleBundleGetIndex(qg, field_name, 3, __RC__)

   field_name = trim(MAM_NH3_CONSTITUENT_NAME)
   iq_nh3     = MAPL_SimpleBundleGetIndex(qg, field_name, 3, __RC__)

   field_name = trim(MAM_NH3_CONSTITUENT_NAME)
   iq_soag    = MAPL_SimpleBundleGetIndex(qg, field_name, 3, __RC__)

   
   ! Gas-aerosol exchange
   ! --------------------
   do j = 1, jm
       do i = 1, im

           do k = 1, km

               ! mid level pressure, temperature and air density
               pressure = 0.5 * (ple(i,j,k-1) + ple(i,j,k))
               z        = 0.5 * (zle(i,j,k-1) + zle(i,j,k))

               temperature  = T(i,j,k)
               density_air  = rhoa(i,j,k)
               rel_humidity = rh(i,j,k)

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

#if(0)
               ! gas-aerosol exchange
               MAML_GasAerosolExchange(pressure,          &
                                       temperature,       &
                                       density_air,       &
                                       rel_humidity,      &
                                       f_cld,             &
                                       z,                 &
                                       pbl_height,        &
                                       q_number,          &
                                       q_nh4,             &
                                       q_so4,             &
                                       q_h2so4,           &
                                       q_nh3,             &
                                       do_nh3,            &
                                       do_nh4g,           &
                                       do_msag,           &
                                       do_soag,           &
                                       Dg,                &
                                       Dg_min,            &
                                       Dg_max,            &
                                       density_so4,       &
                                       mw_so4a,           &
                                       mw_nh4a,           &
                                       dq_h2so4_gasprod,  &
                                       dq_h2so4_aeruptk,  &
                                       cdt)

               ! update the number and mass mixing ratios
               qa%r3(iq_number)%q(i,j,k) = q_number

               qa%r3(iq_amm)%q(i,j,k)    = q_amm
               qa%r3(iq_su)%q(i,j,k)     = q_su

               qg%r3(iq_h2so4)%q(i,j,k)  = q_h2so4
               qg%r3(iq_nh3)%q(i,j,k)    = q_nh3

           end do ! k
#endif

       end do ! i
   end do ! j
#endif
   end subroutine MAM_GasAerosolExchange


end module MAM_GasAerosolExchangeMod
           
