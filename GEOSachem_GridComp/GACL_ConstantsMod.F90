!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: GACL_ConstantsMod
!
! !INTERFACE:
!
   module GACL_ConstantsMod
!
! !USES:
!
#ifdef MAPL
   use MAPL_ConstantsMod, only : MAPL_PI, MAPL_GRAV, MAPL_AVOGAD, MAPL_RUNIV, &
                                 MAPL_AIRMW, MAPL_H2OMW, MAPL_O3MW, &
                                 MAPL_TICE
#endif

   implicit none
   private

!
! !PUBLIC MEMBER FUNCTIONS:

!
! !PUBLIC PARAMETERS:
#ifdef MAPL
   real, parameter, public :: pi       = MAPL_PI                ! pi 
   real, parameter, public :: N_avog   = 1e-3 * MAPL_AVOGAD     ! Avogadro's constant,                '1 mol-1'
   real, parameter, public :: R_univ   = 1e-3 * MAPL_RUNIV      ! Universal/ideal gas constant,       'J K-1 mol-1'
   real, parameter, public :: g_earth  = MAPL_GRAV              ! standard gravity,                   'm s-2'
   real, parameter, public :: T_ice    = MAPL_TICE              ! melting point of ice,               'K'
   real, parameter, public :: mw_air   = 1e-3 * MAPL_AIRMW      ! molar mass of dry air,              'kg mol-1'
   real, parameter, public :: mw_H2O   = 1e-3 * MAPL_H2OMW      ! molar mass of water,                'kg mol-1'
   real, parameter, public :: mw_O3    = 1e-3 * MAPL_O3MW       ! molar mass of ozone,                'kg mol-1'
#else
   real, parameter, public :: pi       = 3.141592653589793      ! pi 
   real, parameter, public :: N_avog   = 6.022e23               ! Avogadro's constant,                '1 mol-1'
   real, parameter, public :: R_univ   = 8.31447                ! Universal/ideal gas constant,       'J K-1 mol-1'
   real, parameter, public :: g_earth  = 9.80665                ! standard gravity,                   'm s-2'
   real, parameter, public :: T_ice    = 273.16                 ! Melting point of ice,               'K'
   real, parameter, public :: mw_air   = 28.965e-3              ! molar mass of dry air,              'kg mol-1'
   real, parameter, public :: mw_H2O   = 18.015e-3              ! molar mass of water,                'kg mol-1'
   real, parameter, public :: mw_O3    = 47.9982e-3             ! molar mass of ozone,                'kg mol-1'
#endif

   real, parameter, public :: mw_S     = 32.065e-3              ! atomic mass of sulfur,              'kg mol-1'
   real, parameter, public :: mw_SO2   = 64.066e-3              ! molar mass of sulfur dioxide,       'kg mol-1'
   real, parameter, public :: mw_SO4   = 96.07e-3               ! molar mass of sulfate,              'kg mol-1'
   real, parameter, public :: mw_H2SO4 = 98.079e-3              ! molar mass of sulfuric acid,        'kg mol-1'
   real, parameter, public :: mw_DMS   = 62.13e-3               ! molar mass of dimethyl sulfide,     'kg mol-1'
   real, parameter, public :: mw_MSA   = 96.11e-3               ! molar mass of methanesulfonic acid, 'kg mol-1'
   real, parameter, public :: mw_OH    = 17.01e-3               ! molar mass of hydroxyl radical,     'kg mol-1'
   real, parameter, public :: mw_HO2   = 33.01e-3               ! molar mass of hydroperoxyl radical, 'kg mol-1'
   real, parameter, public :: mw_H2O2  = 34.0147e-3             ! molar mass of hydrogen peroxide,    'kg mol-1'
   real, parameter, public :: mw_N     = 14.007e-3              ! atomic mass of nitrogen,            'kg mol-1'
   real, parameter, public :: mw_NO3   = 62.0049e-3             ! molar mass of nitrate ion,          'kg mol-1'
   real, parameter, public :: mw_NH3   = 17.031e-3              ! molar mass of ammonia,              'kg mol-1'
   real, parameter, public :: mw_NH4   = 18.0385e-3             ! molar mass of ammonium,             'kg mol-1'
   real, parameter, public :: mw_OCS   = 60.075e-3              ! molar mass of carbonyl sulfide,     'kg mol-1'
   real, parameter, public :: mw_SOAg  = 12.0e-3                ! molar mass of SOA(gas),             'kg mol-1'

!
! !PRIVATE PARAMETERS:

!
! !DESCRIPTION: 
!
!  {\tt GACL\_ConstantsMod} defines physics and chemistry constants.
!
!
! !REVISION HISTORY:
!
!  06June2014  A. Darmenov  Initial version
!
!EOP
!-------------------------------------------------------------------------

   end module GACL_ConstantsMod

