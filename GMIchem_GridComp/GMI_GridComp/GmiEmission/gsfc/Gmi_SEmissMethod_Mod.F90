#include "MAPL_Generic.h"
!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: Gmi_SEmissMethod_Mod - DMS emissions of sulfur precursor gases
!
! !INTERFACE:
!
 module Gmi_SEmissMethod_Mod
!
! !USES:
!
!   use GACL_ConstantsMod, only : g_earth, T_ice, mw_air, mw_H2O, mw_DMS
! From GACL_ConstantsMod USES:
   use MAPL_ConstantsMod, only : MAPL_GRAV, MAPL_AIRMW, MAPL_TICE, MAPL_H2OMW
!
   IMPLICIT NONE
   real, parameter, public :: g_earth  = MAPL_GRAV              ! standard gravity,                   'm s-2'
   real, parameter, public :: T_ice    = MAPL_TICE              ! Melting point of ice,               'K'
   real, parameter, public :: mw_air   = 1e-3 * MAPL_AIRMW      ! molar mass of dry air,              'kg mol-1'
   real, parameter, public :: mw_H2O   = 1e-3 * MAPL_H2OMW      ! molar mass of water,                'kg mol-1'
   real, parameter, public :: mw_DMS   = 62.13e-3               ! molar mass of dimethyl sulfide,     'kg mol-1'
!
   public GmiDMS_Emissions
!   public GMI_SO2Volc_Emissions
   private GmiDMS_flux
   private GOCARTDMS_flux
!
   contains
!
!
!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: DMS_emissions --- DMS emissions from ocean. Emissions are 
!            injected in the surface model layer.
!
! !INTERFACE:

   subroutine GmiDMS_emissions(do_achem_dms_emiss &
                            ,t_bottom    &
                            ,u10n        &
                            ,v10n        &
                            ,fr_ocean    &
                            ,DMS_ocean   &
                            ,DMS_atm     &
                            ,flux        &
                            ,rc)
!
   implicit None

! !INPUT/OUTPUT PARAMETERS:

! !INPUT PARAMETERS:
  logical,   intent(in)    :: do_achem_dms_emiss        ! which DMS flux calc to use
   real, dimension(:,:),   intent(in)  :: t_bottom      ! skin temperature, K
   real, dimension(:,:),   intent(in)  :: u10n          ! equivalent neutral wind speed at 10m 
   real, dimension(:,:),   intent(in)  :: v10n          ! equivalent neutral wind speed at 10m
   real, dimension(:,:),   intent(in)  :: fr_ocean      ! fraction of ocean 
   real, dimension(:,:),   intent(in)  :: DMS_ocean     ! sea surface concentrations of DMS (kg/m^3)
   real, dimension(:,:,:), intent(in)  :: DMS_atm       ! DMS mixing ratio, mol mol-1
!
! !OUTPUT PARAMETERS:  
   real, dimension(:,:),   intent(out) :: flux          ! DMS flux, mol m-2 s-1
   integer               , intent(out) :: rc            ! return code

! !DESCRIPTION: 
!
!
! !REVISION HISTORY:
!
!  10Apr2021  S. Steenrod
!
!EOP
!-------------------------------------------------------------------------
!
! local
   real    :: tbot
!
   integer :: i, i1, i2
   integer :: j, j1, j2
   integer :: k1
!
   rc = 0
!
   i1 = lbound(DMS_atm, 1)
   i2 = ubound(DMS_atm, 1)
   j1 = lbound(DMS_atm, 2)
   j2 = ubound(DMS_atm, 2)
   k1 = lbound(DMS_atm, 3)
!
   flux = 0.0
!   forall (i = i1:i2, j = j1:j2, (fr_ocean(i,j) > 0.0) .and. (t_bottom(i,j) > T_ice))
   do j=j1,j2
     do i=i1,i2
       tbot = t_bottom(i,j)
       if ((fr_ocean(i,j) > 0.0) .and. (tbot > T_ice)) then
!
         if(do_achem_dms_emiss) then
!... ACHEM DMS emission method - emissions way too high, unit issue?
           if(tbot > T_ice) &
            flux(i,j) = GmiDMS_flux(DMS_ocean(i,j), DMS_atm(i,j,k1), u10n(i,j), v10n(i,j), tbot)
         else
!... GOCART DMS emission method
           flux(i,j) = GOCARTDMS_flux (DMS_ocean(i,j), u10n(i,j), v10n(i,j), tbot)
         endif
!
       else
         flux(i,j) = 0.0
       endif
!
     enddo
   enddo
!
   end subroutine GmiDMS_Emissions
 

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: DMS_flux --- Computes sea-to-air DMS flux. 
!
! !INTERFACE:

   pure real function GmiDMS_flux(DMS_ocean, DMS_atmosphere, u10n, v10n, SST)

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

! Ostwald solubility coefficient for DMS (order 5-50)
   alpha = exp(3525.0/SST - 9.464)

! airside transfer velocity
   k_a = f_a * w10n             ! k_a = (659 * w10n) * sqrt(mw_DMS / mw_H2O)

! atmospheric gradient fraction
   gamma_a = 1.0/(1.0 + k_a / (alpha * k_w))

! total gas transfer velocity
   k = k_w * (1.0 - gamma_a)   ! cm h-1
   k = (1e-2/3600) * k         ! converted to m s-1
!
! DMS emission flux, mol m-2 s-1
#if(1)
   GmiDMS_flux = k * (DMS_ocean - alpha * DMS_atmosphere)
#else
   GmiDMS_flux = k * DMS_ocean
#endif
!
   GmiDMS_flux = max(0.0, GmiDMS_flux)
!
   end function GmiDMS_flux
!
!
!... GOCART DMS emission
!
!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!
! !INTERFACE:
   pure real function GOCARTDMS_flux(DMS_ocean, u10n, v10n, Tsfc)
!
! !USES:
   implicit None
!
! !INPUT/OUTPUT PARAMETERS:
!
! !INPUT PARAMETERS:
   real, intent(in) :: DMS_ocean            ! ocean DMS concentration, kg m-3
   real, intent(in) :: u10n, v10n           ! equivalient neutral wind speed at 10m, m s-1
   real, intent(in) :: Tsfc                 ! surface layer temperature, K
!
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   Computes sea-to-air DMS flux like GOCART
!
! !REVISION HISTORY:
!  2021May  S. Steenrod
!
!EOP
!-------------------------------------------------------------------------
!
! parameters
   real, parameter :: sCO2 = 600.

! local
   real :: sst           ! SST, C
   real :: w10n          ! equivalient neutral wind speed at 10m, m s-1
   real :: akw           ! total gas transfer velocity, cm h-1
   real :: schmidt       ! Schmidt number for DMS  
!
!... Tsfc from K to C and cap at 28C
   sst = Tsfc - 273.15                                         
   sst = min(sst,28.0)                                       
!
! equivalent neutral wind speed at 10 meters
   w10n = sqrt(u10n*u10n + v10n*v10n)
!... only valid for ocean and warm enough temperatures              
   if(sst .gt. -20.) then
     schmidt = 2764.0 - 147.12*sst + 3.726*(sst*sst) - 0.038*(sst*sst*sst)
!
     if(w10n .le. 3.6) then                                           
       akw = 0.17*w10n*((sCO2/schmidt)**0.667)
     else if (w10n .le. 13.) then
       akw = (2.85*w10n - 9.65)*sqrt(sCO2/schmidt)
     else
       akw = (5.90*w10n - 49.3)*sqrt(sCO2/schmidt)
     endif
!... This parameterization has put akw in units cm hr-1 -> goto m s-1
     akw = akw/100./3600.
   else
     akw = 0.0
   endif
!... DMS_ocean concentration is kg m-3
!... output units: kg m-2 s-1
   GOCARTDMS_flux = akw * DMS_ocean
!... floor
   GOCARTDMS_flux = max(0.0, GOCARTDMS_flux)
!
   end function GOCARTDMS_flux
!
 end module Gmi_SEmissMethod_Mod
