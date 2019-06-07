#include "MAPL_Exceptions.h"
#include "MAPL_Generic.h"


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: GACL_EmissionsMod - Primary emissions of aerosol precursor gases
!
! !INTERFACE:
!
   module GACL_EmissionsMod
!
! !USES:
!

   use ESMF, only : ESMF_Grid, ESMF_Config,          &
                               ESMF_ConfigCreate,    &
                               ESMF_ConfigDestroy,   &
                               ESMF_ConfigLoadFile,  &
                               ESMF_ConfigGetDim,    &
                               ESMF_ConfigFindLabel, &
                               ESMF_ConfigNextLine,  & 
                               ESMF_ConfigGetAttribute

   use MAPL_Mod 

   use m_StrTemplate, only : StrTemplate

   use GACL_ConstantsMod, only : g_earth, T_ice, N_avog, &
                                 mw_air, mw_S, mw_SO2, mw_NH3, mw_DMS, mw_H2O


   implicit NONE
   private

!
! !PUBLIC MEMBER FUNCTIONS:

   public NH3_Emissions
   public SO2_Emissions
   public DMS_Emissions
   public SOAG_Emissions
   public VOC_Emissions
!
! !PUBLIC PARAMETERS:

!
! !PRIVATE PARAMETERS:

!
! !DESCRIPTION: 
!
!  {\tt GACL\_EmissionsMod} - emissions of sulfate and ammonia.
!
! !REVISION HISTORY:
!
!  29Sep2011  A. Darmenov  Initial version
!
!EOP
!-------------------------------------------------------------------------


   contains

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NH3_emissions --- Ammonia (NH3) emissions from natural and 
!            anthropogenic sources. Emissions are injected in the 
!            surface model layer.
!
! !INTERFACE:

   subroutine NH3_Emissions(delp,          &
                            emiss_lumped,  &
                            emiss_bb,      &
                            q,             &
                            cdt,           &
                            rc)
! !USES:

   implicit None

! !INPUT/OUTPUT PARAMETERS:

! !INPUT PARAMETERS:

   real, dimension(:,:,:), intent(in)    :: delp            ! pressure level thickness, Pa

   real, dimension(:,:),   intent(in)    :: emiss_lumped    ! emissions of NH3 (not including biomass burning)
   real, dimension(:,:),   intent(in)    :: emiss_bb        ! emissions of NH3 from biomass burning

   real, intent(in)                      :: cdt             ! time step


! !OUTPUT PARAMETERS:
   
   real, dimension(:,:,:), intent(inout) :: q               ! NH3 mixing ratio, mol mol-1

   integer, intent(out)                  :: rc              ! return code   

! !DESCRIPTION: Emissions of NH3. Emissions are injected in the first surface layer.
!
!
! !REVISION HISTORY:
!
!  29Sep2012  A. Darmenov
!
!EOP
!-------------------------------------------------------------------------

                   __Iam__('NH3_emissions')

   ! local 
   real    :: f
   integer :: k1, km
   
   rc = 0

   k1 = lbound(q, 3)
   km = ubound(q, 3)
   
   f = (mw_air / mw_NH3) * g_earth * cdt
   q(:,:,km) = q(:,:,km) + f * (emiss_lumped(:,:) + emiss_bb(:,:)) / delp(:,:,km)

   end subroutine NH3_Emissions



!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: SO2_emissions --- Sulfur dioxide (SO2) emissions from natural and 
!            anthropogenic sources. For now anthropogenic (except aviation sector) 
!            and fire emissions are injected in the surface model layer.
!
! !INTERFACE:

   subroutine SO2_Emissions(delp,                          &
                            zle,                           & 
                            rho_air,                       &
                            emiss_bb,                      &
                            emiss_nonenergy,               &
                            emiss_energy,                  &
                            emiss_ship,                    &
                            emiss_aircraft_lto,            &
                            emiss_aircraft_cds,            &
                            emiss_aircraft_crs,            &
                            aviation_layers,               &
                            n_volcanos,                    &
                            volc_elev,                     &
                            volc_cloud,                    &
                            volc_SO2,                      &
                            volc_start, volc_end,          &
                            volc_i, volc_j,                &
                            emiss_volcanic_expl,           &
                            emiss_volcanic_nexp,           &
                            emiss_tot,                     &
                            q,                             &
                            cell_area,                     &
                            cdt,                           &
                            nymd,                          &
                            nhms,                          &
                            rc)

! !USES:

   implicit None

! !INPUT/OUTPUT PARAMETERS:

! !INPUT PARAMETERS:

   real, dimension(:,:,:), intent(in)    :: delp                ! pressure level thickness, Pa
   real, dimension(:,:,:), intent(in)    :: rho_air             ! air density, kg m-3
   real, dimension(:,:,0:),intent(in)    :: zle                 ! edge heights, m

   real, dimension(:,:),   intent(in)    :: emiss_bb            ! emissions of SO2 from biomass burning
   real, dimension(:,:),   intent(in)    :: emiss_nonenergy     ! emissions of SO2 from non-energy sectors
   real, dimension(:,:),   intent(in)    :: emiss_energy        ! emissions of SO2 from energy sectors
   real, dimension(:,:),   intent(in)    :: emiss_ship          ! emissions of SO2 from ships

   real, dimension(:,:),   intent(in)    :: emiss_aircraft_lto  ! emissions of SO2 from aircraft - LTO layer
   real, dimension(:,:),   intent(in)    :: emiss_aircraft_cds  ! emissions of SO2 from aircraft - CDS layer
   real, dimension(:,:),   intent(in)    :: emiss_aircraft_crs  ! emissions of SO2 from aircraft - CRS layer
   real, dimension(4),     intent(in)    :: aviation_layers     ! extend of the LTO, CDS amd CRS layers

   integer                               :: n_volcanos
   real,    pointer, dimension(:)        :: volc_elev, volc_SO2, volc_cloud
   integer, pointer, dimension(:)        :: volc_start, volc_end
   integer, pointer, dimension(:)        :: volc_i, volc_j

   integer, intent(in)                   :: nymd                ! current date
   integer, intent(in)                   :: nhms                ! current time

   real, dimension(:,:),   intent(in)    :: cell_area           !


   real, intent(in)                      :: cdt                 ! time step

! !OUTPUT PARAMETERS:
   real, dimension(:,:), intent(inout)   :: emiss_volcanic_expl ! emissions from explosive volcanoes
   real, dimension(:,:), intent(inout)   :: emiss_volcanic_nexp ! emissions from non-explosive volcanoes

   real, dimension(:,:), intent(inout)   :: emiss_tot           ! diagnostics: total emissions of SO2, 'kg m-2 s-1'
   real, dimension(:,:,:), intent(inout) :: q                   ! SO2 mixing ratio, mol mol-1

   integer, intent(out)                  :: rc                  ! return code

! !DESCRIPTION: 
!
!
! !REVISION HISTORY:
!
!  29Sep2012  A. Darmenov
!
!EOP
!-------------------------------------------------------------------------

                   __Iam__('SO2_Emissions')

   ! local 
   real               :: f

   integer            :: i, i1, i2
   integer            :: j, j1, j2
   integer            :: k, k1, km

   integer            :: it

   real               :: z_lto_bot, z_lto_top
   real               :: z_cds_bot, z_cds_top
   real               :: z_crs_bot, z_crs_top

   real, allocatable, dimension(:,:,:) :: emiss_aviation_layer
   real, allocatable, dimension(:,:,:) :: emiss_aviation

   real                                :: so2volcano  

   real, allocatable, dimension(:,:)   :: z0

   real                                :: hup, hlow, dz_volc
   real                                :: dz, z1
   real                                :: deltaSO2v

   rc = 0

   i1 = lbound(q, 1); i2 = ubound(q, 1)
   j1 = lbound(q, 2); j2 = ubound(q, 2)
   k1 = lbound(q, 3); km = ubound(q, 3)
   
   f = (mw_air / mw_SO2) * g_earth * cdt

   ! for now inject emissions in the first model layer
   q(:,:,km) = q(:,:,km) + f * (emiss_bb(:,:)        + &
                                emiss_nonenergy(:,:) + &
                                emiss_energy(:,:)    + &
                                emiss_ship(:,:)) / delp(:,:,km)
   
   ! aircraft emissions: LTO, CDS and CRS layers
   z_lto_bot = max(1e-3, aviation_layers(1))
   z_lto_top = max(2e-3, aviation_layers(2))

   z_cds_bot = max(2e-3, aviation_layers(2))
   z_cds_top = max(3e-3, aviation_layers(3))

   z_crs_bot = max(3e-3, aviation_layers(3))
   z_crs_top = max(4e-3, aviation_layers(4))

   allocate(emiss_aviation_layer(i1:i2,j1:j2,km), __STAT__)
   allocate(emiss_aviation(i1:i2,j1:j2,km),       __STAT__)

   emiss_aviation_layer = 0.0
   emiss_aviation       = 0.0

   call distribute_aviation_emissions(delp, rho_air, z_lto_bot, z_lto_top, emiss_aircraft_lto, emiss_aviation_layer, i1, i2, j1, j2, km)
   emiss_aviation = emiss_aviation + emiss_aviation_layer

   call distribute_aviation_emissions(delp, rho_air, z_cds_bot, z_cds_top, emiss_aircraft_cds, emiss_aviation_layer, i1, i2, j1, j2, km)
   emiss_aviation = emiss_aviation + emiss_aviation_layer

   call distribute_aviation_emissions(delp, rho_air, z_crs_bot, z_crs_top, emiss_aircraft_crs, emiss_aviation_layer, i1, i2, j1, j2, km)
   emiss_aviation = emiss_aviation + emiss_aviation_layer

   deallocate(emiss_aviation_layer, __STAT__)

   q(:,:,:) = q(:,:,:)  + f * emiss_aviation(:,:,:) / delp(:,:,:)

   ! volcanic emissions
   ! Point source volcanos (loop over each volcano)
   allocate(z0(i1:i2,j1:j2), __STAT__)
   z0 = zle(:,:,km)

   emiss_volcanic_expl = 0.0
   emiss_volcanic_nexp = 0.0

   if (n_volcanos > 0) then

       VOLCANOES: do it = 1, n_volcanos

           i = volc_i(it)
           j = volc_j(it)
         
           ! skip this volcano?
           if ((i < 1) .or. (j < 1)) cycle ! volcano not in sub-domain

           ! check time against time range of eruption
           if (nhms < volc_start(it) .or. nhms >= volc_end(it)) cycle

           so2volcano = 0.0

           ! emissions per volcano
           if (cell_area(i,j) > 1.0) then                    ! omit volcanos in very small grid boxes
               so2volcano = volc_so2(it) / cell_area(i,j)    ! to 'kg(SO2) s-1 m-2'
               so2volcano = max(so2volcano, tiny(so2volcano))
           endif

           ! distribute in the vertical
           hup  = volc_cloud(it)
           hlow = volc_elev(it)

           if (hup .ne. hlow) then
               hlow = hup - (hup - hlow)/3.0
           endif

           ! diagnostic - sum of volcanos
           if (hup .eq. hlow) then
               emiss_volcanic_nexp(i,j) = emiss_volcanic_nexp(i,j) + so2volcano
           else
               emiss_volcanic_expl(i,j) = emiss_volcanic_expl(i,j) + so2volcano
           end if

           dz_volc = hup - hlow

           VERTICAL_LEVELS: do k = km, 1, -1
               z1 = zle(i,j,k-1)
               dz = z1 - z0(i,j)
               deltaSO2v = 0.0

               ! volcano is above this level
               if(z1 .lt. hlow) then
                   z0(i,j) = z1
                   cycle
               end if

               ! volcano is below this level
               if (z0(i,j) .gt. hup) then
                   z0(i,j) = z1
                   cycle
               end if

               ! volcano is in this level
               if ((k .eq. km .and. z0(i,j) .gt. hup) .or. &      ! below surface
                   (z0(i,j) .le. hlow .and. z1 .ge. hup))  then   ! in level
                   deltaSO2v = so2volcano

               ! volcano only partly in level                     ! cell:
               else if (z0(i,j) .lt. hlow .and. z1 .lt. hup) then ! has bottom of cloud
                   deltaSO2v = (z1 - hlow)/dz_volc * so2volcano
 
               else if (z0(i,j) .gt. hlow .and. z1 .gt. hup) then ! has top of cloud
                   deltaSO2v = (hup - z0(i,j))/dz_volc * so2volcano
 
               else                                               ! is filled with cloud
                   deltaSO2v = dz/dz_volc * so2volcano
               end if

               z0(i,j) = z1

               q(i,j,k) = q(i,j,k) + (mw_air / mw_SO2)*deltaSO2v*cdt*g_earth/delp(i,j,k)
           end do VERTICAL_LEVELS
       end do VOLCANOES

   endif
     
   ! diagnostics - total SO2 emissions
   emiss_tot = (emiss_bb            + &
                emiss_nonenergy     + &
                emiss_energy        + &
                emiss_ship          + &
                emiss_volcanic_expl + &
                emiss_volcanic_nexp + &
                sum(emiss_aviation, dim=3))


   deallocate(emiss_aviation, __STAT__)
   deallocate(z0,             __STAT__)

contains

    subroutine distribute_aviation_emissions(delp, rhoa, z_bot, z_top, emissions_layer, emissions, i1, i2, j1, j2, km)

    implicit none

    integer, intent(in) :: i1, i2, j1, j2, km

    real, dimension(:,:,:), intent(in) :: delp
    real, dimension(:,:,:), intent(in) :: rhoa
    real, dimension(:,:),   intent(in) :: emissions_layer
    real, intent(in)                   :: z_bot
    real, intent(in)                   :: z_top
    real, dimension(:,:,:), intent(out):: emissions
    
!   local
    integer :: i, j, k
    integer :: k_bot, k_top
    real    :: z_
    real, dimension(km) :: z, dz, w_
    
    do j = j1, j2
        do i = i1, i2
            ! find level height
            z = 0.0
            z_= 0.0 

            do k = km, 1, -1
                dz(k) = delp(i,j,k)/rhoa(i,j,k)/g_earth
                z_    = z_ + dz(k)
                z(k)  = z_
            end do

            ! find the bottom level
            do k = km, 1, -1
                if (z(k) >= z_bot) then
                    k_bot = k
                    exit
                end if
            end do
            
            ! find the top level
            do k = k_bot, 1, -1
                if (z(k) >= z_top) then
                    k_top = k
                    exit
                end if
            end do

            ! find the weights
            w_ = 0

!           if (k_top > k_bot) then
!               need to bail - something went wrong here
!           end if

            if (k_bot .eq. k_top) then
                w_(k_bot) = z_top - z_bot
            else
                do k = k_bot, k_top, -1
                    if ((k < k_bot) .and. (k > k_top)) then
                        w_(k) = dz(k)
                    else
                        if (k == k_bot) then
                            w_(k) = (z(k) - z_bot)
                        end if

                        if (k == k_top) then
                            w_(k) = z_top - (z(k)-dz(k))
                        end if
                    end if
                end do
            end if
           
            ! distribute emissions in the vertical 
            emissions(i,j,:) = (w_ / sum(w_)) * emissions_layer(i,j)
        end do 
    end do

    return 

    end subroutine distribute_aviation_emissions

   end subroutine SO2_Emissions



!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: DMS_emissions --- DMS emissions from ocean. Emissions are 
!            injected in the surface model layer.
!
! !INTERFACE:

   subroutine DMS_emissions(delp,        &
                            t_skin,      &
                            u10n,        &
                            v10n,        &
                            fr_ocean,    &
                            DMS_ocean,   &
                            q,           &
                            flux,        &
                            cdt,         &
                            rc)
! !USES:

   implicit None

! !INPUT/OUTPUT PARAMETERS:

! !INPUT PARAMETERS:

   real, dimension(:,:,:), intent(in)    :: delp          ! pressure level thickness, Pa

   real, dimension(:,:),   intent(in)    :: t_skin        ! skin teneperature, K

   real, dimension(:,:),   intent(in)    :: u10n          ! equivalient neutral wind speed at 10m 
   real, dimension(:,:),   intent(in)    :: v10n          ! equivalient neutral wind speed at 10m

   real, dimension(:,:),   intent(in)    :: fr_ocean      ! fraction of ocean 

   real, dimension(:,:),   intent(in)    :: DMS_ocean     ! sea surface concentrations of DMS

   real, intent(in)                      :: cdt           ! time step

! !OUTPUT PARAMETERS:
   
   real, dimension(:,:,:), intent(inout) :: q             ! DMS mixing ratio, mol mol-1
   real, dimension(:,:),   intent(inout) :: flux          ! DMS flux, mol m-2 s-1

   integer, intent(out)                  :: rc            ! return code

! !DESCRIPTION: 
!
!
! !REVISION HISTORY:
!
!  14Aug2012  A. Darmenov    Fisrst crack
!
!EOP
!-------------------------------------------------------------------------


                     __Iam__('DMS_Emissions')

   ! local
   real    :: f

   integer :: i, i1, i2
   integer :: j, j1, j2
   integer :: k1, km

   rc = 0

   i1 = lbound(q, 1); i2 = ubound(q, 1)
   j1 = lbound(q, 2); j2 = ubound(q, 2)
   k1 = lbound(q, 3); km = ubound(q, 3)

   f = mw_air * cdt * g_earth
   flux = 0.0
   forall (i = i1:i2, j = j1:j2, (fr_ocean(i,j) > 0.0) .and. (t_skin(i,j) > T_ice))
       flux(i,j) = DMS_flux(DMS_ocean(i,j), &
                            q(i,j,km),      &
                            u10n(i,j),      &
                            v10n(i,j),      &
                            t_skin(i,j))

       q(i,j,km) = q(i,j,km) + f * flux(i,j) / delp(i,j,km)
   end forall

   end subroutine DMS_Emissions
 


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: DMS_flux --- Computes sea-to-air DMS flux. 
!
! !INTERFACE:

   pure real function DMS_flux(DMS_ocean, DMS_atmosphere, u10n, v10n, SST)

! !USES:

   implicit None

! !INPUT/OUTPUT PARAMETERS:

! !INPUT PARAMETERS:
   
   real, intent(in) :: DMS_ocean            ! ocean DMS concentration, mol m-3
   real, intent(in) :: DMS_atmosphere       ! DMS concentration, mol mol-1

   real, intent(in) :: u10n, v10n           ! equivalient neutral wind speed at 10m, m s-1

   real, intent(in) :: SST                  ! sea surface temperature (SST), K

! !OUTPUT PARAMETERS:

! !DESCRIPTION: Computes sea-to-air DMS flux.
!
!
! !REVISION HISTORY:
!
!  14Aug2012  A. Darmenov    First crack
!
!EOP
!-------------------------------------------------------------------------

   !               __Iam__('DMS_flux')

   ! parameters
   real, parameter :: f_a = 659 * sqrt(mw_DMS / mw_H2O)

   ! local
   real :: k_Sc600       ! gas transfer coefficient for a Schmidt number of 600 
   real :: k_w           ! water side DMS gas transfer velocity, cm h-1
   real :: k_a           ! airside DMS gas transfer velosity, cm h-1

   real :: k             ! total gas transfer velocity, cm h-1

   real :: gamma_a       ! airside gradient fraction 
   real :: alpha         ! Ostwald solubility coefficient, alpha = H * (R * T * water_density), kH is Henry's law coefficient
   real :: Sc_DMS        ! Schmidt number for DMS
  
   real :: SST_C         ! SST, C
   real :: w10n          ! equivalient neutral wind speed at 10m, m s-1
   

   ! equivalent neutral wind speed at 10 meters
   w10n = sqrt(u10n*u10n + v10n*v10n)

   ! water side DMS gas transfer velocity is based on the 10 m wind‐speed‐based 
   ! parameterization of Nightingale et al. [2000]
   k_Sc600 = 0.222*w10n**2 + 0.333*w10n
   SST_c   = min(max(0.0, SST - 273.15), 35.0)
   Sc_DMS  = 2764.0 + SST_c*(-147.12 + SST_c*(3.726 + SST_c*(-0.038)))

   k_w = k_Sc600 * sqrt(Sc_DMS / 600.0)

   ! Ostwald solubility coefficient for DMS
   alpha = exp(3525.0/SST - 9.464)

   ! airside transfer velocity
   k_a = f_a * w10n             ! k_a = (659 * w10n) * sqrt(mw_DMS / mw_H2O)

   ! atmospheric gradient fraction
   gamma_a = 1.0/(1.0 + k_a / (alpha * k_w))

   ! total gas transfer velocity
   k = k_w * (1.0 - gamma_a)   ! cm h-1
   k = (1e-2/3600) * k         ! converted to m s-1

   ! DMS emission flux, mol m-2 s-1
#if(1)
   DMS_flux = k * (DMS_ocean - alpha * DMS_atmosphere)
#else
   DMS_flux = k * DMS_ocean
#endif

   DMS_flux = max(0.0, DMS_flux)

   end function DMS_flux



!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: SOAG_emissions --- SOA(gas) emissions from natural and 
!            anthropogenic sources. Emissions are injected in the 
!            surface model layer.
!
! !INTERFACE:

   subroutine SOAG_Emissions(delp,          &
                             emiss_lumped,  &
                             q,             &
                             cdt,           &
                             rc)
! !USES:

   implicit None

! !INPUT/OUTPUT PARAMETERS:

! !INPUT PARAMETERS:

   real, dimension(:,:,:), intent(in)    :: delp            ! pressure level thickness, Pa

   real, dimension(:,:),   intent(in)    :: emiss_lumped    ! lumped emissions of SOA(gas), 'molecules m-2 s-1'

   real, intent(in)                      :: cdt             ! time step


! !OUTPUT PARAMETERS:
   
   real, dimension(:,:,:), intent(inout) :: q               ! SOA(gas) mixing ratio, mol mol-1

   integer, intent(out)                  :: rc              ! return code   

! !DESCRIPTION: 
!
!
! !REVISION HISTORY:
!
!  29Sep2012  A. Darmenov
!
!EOP
!-------------------------------------------------------------------------

                   __Iam__('SOAG_emissions')

   ! local 
   real    :: f
   integer :: k1, km
   
   rc = 0

   k1 = lbound(q, 3)
   km = ubound(q, 3)
   
   f = (mw_air / N_avog) * g_earth * cdt
   q(:,:,km) = q(:,:,km) + f * emiss_lumped(:,:) / delp(:,:,km)

   end subroutine SOAG_Emissions


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: VOC_emissions ---  
!            VOC emissions from biomass burning and anthropogenic sources,
!            based on corresponding CO emissions. Emissions are injected in
!            the model surface layer.
!
! !INTERFACE:

   subroutine VOC_Emissions(delp,          &
                            voc_BiomassBurnFactor,  &
			    voc_AnthroFactor,  &
			    co_biomass_voc, &
			    co_bf_voc, &
			    co_fs_voc, &
			    voc_MW, &
                            q, qb,         &
                            cdt,           &
                            rc)
! !USES:

   implicit None

! !INPUT/OUTPUT PARAMETERS:

! !INPUT PARAMETERS:

   real, dimension(:,:,:), intent(in)    :: delp            ! pressure level thickness, Pa

   real, dimension(:,:),   intent(in)    :: co_biomass_voc  ! CO biomass burning emissions, kg m-2 s-1
   real, dimension(:,:),   intent(in)    :: co_bf_voc       ! CO biofuel emissions, kg m-2 s-1
   real, dimension(:,:),   intent(in)    :: co_fs_voc       ! CO fossil fuel emissions, kg m-2 s-1

   real, intent(in) :: voc_BiomassBurnFactor                ! 'g/g CO'
   real, intent(in) :: voc_AnthroFactor                     ! 'g/g CO'
   real, intent(in) :: voc_MW
   
   real, intent(in)                      :: cdt             ! time step


! !OUTPUT PARAMETERS:
   real, dimension(:,:,:), intent(inout) :: q               ! VOC mixing ratio, mol mol-1 (anthro)
   real, dimension(:,:,:), intent(inout) :: qb              ! VOC mixing ratio, mol mol-1 (biob)
   integer, intent(out)                  :: rc              ! return code   

! !DESCRIPTION: 
!
!
! !REVISION HISTORY:
!
!  07Oct2016  M.S. Johnson/P.R. Colarco
!
!-------------------------------------------------------------------------
                     __Iam__('VOC_emissions')

   ! local variables
   real                                  :: f
   real, allocatable, dimension(:,:)     :: dvoc
   integer                               :: i1, i2, j1, j2, k1, km
   
   rc = 0

   i1 = lbound(q, 1); i2 = ubound(q, 1)
   j1 = lbound(q, 2); j2 = ubound(q, 2)
   k1 = lbound(q, 3); km = ubound(q, 3)
   
   ! The scaling here results in a change in the volume mixing ratio of VOC
   f =  (mw_air / voc_MW) * g_earth * cdt
   allocate(dvoc(i1:i2,j1:j2), __STAT__)
   ! Anthropogenic + Biofuel
   dvoc = f * (co_bf_voc + co_fs_voc) * voc_AnthroFactor / delp(:,:,km)
   q(:,:,km) = q(:,:,km) + dvoc
   ! Biomass burning
   dvoc = f *  co_biomass_voc * voc_BiomassBurnFactor / delp(:,:,km)
   qb(:,:,km) = qb(:,:,km) + dvoc
   deallocate(dvoc, __STAT__)
   
   end subroutine VOC_Emissions
!
!EOP
!-------------------------------------------------------------------------


!-------------------------------------------------------------------------
!     NASA/GSFC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  distribute_aviation_emissions - Distributes 2D aviation emissions
!             in the vertical column  

!
! !INTERFACE:
!
    subroutine distribute_aviation_emissions(delp, rhoa, z_bot, z_top, emissions_layer, emissions, i1, i2, j1, j2, km)

    implicit none

    integer, intent(in) :: i1, i2, j1, j2, km

    real, dimension(:,:,:), intent(in) :: delp
    real, dimension(:,:,:), intent(in) :: rhoa
    real, dimension(:,:),   intent(in) :: emissions_layer
    real, intent(in)                   :: z_bot
    real, intent(in)                   :: z_top
    real, dimension(:,:,:), intent(out):: emissions
    
!   local
    integer :: i, j, k
    integer :: k_bot, k_top
    real    :: z_
    real, dimension(km) :: z, dz, w_
    
    do j = j1, j2
        do i = i1, i2
            ! find level height
            z = 0.0
            z_= 0.0 

            do k = km, 1, -1
                dz(k) = delp(i,j,k)/rhoa(i,j,k)/g_earth
                z_    = z_ + dz(k)
                z(k)  = z_
            end do

            ! find the bottom level
            do k = km, 1, -1
                if (z(k) >= z_bot) then
                    k_bot = k
                    exit
                end if
            end do
            
            ! find the top level
            do k = k_bot, 1, -1
                if (z(k) >= z_top) then
                    k_top = k
                    exit
                end if
            end do

            ! find the weights
            w_ = 0

!           if (k_top > k_bot) then
!               need to bail - something went wrong here
!           end if

            if (k_bot .eq. k_top) then
                w_(k_bot) = z_top - z_bot
            else
                do k = k_bot, k_top, -1
                    if ((k < k_bot) .and. (k > k_top)) then
                        w_(k) = dz(k)
                    else
                        if (k == k_bot) then
                            w_(k) = (z(k) - z_bot)
                        end if

                        if (k == k_top) then
                            w_(k) = z_top - (z(k)-dz(k))
                        end if
                    end if
                end do
            end if
           
            ! distribute emissions in the vertical 
            emissions(i,j,:) = (w_ / sum(w_)) * emissions_layer(i,j)
        end do 
    end do

    end subroutine distribute_aviation_emissions

  end module GACL_EmissionsMod 
