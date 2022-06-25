!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: MAM_ComponentsDataMod - MAM chem/aerosol components properties
!
! !INTERFACE:
!
   module MAM_ComponentsDataMod
!
! !USES:
!

   implicit NONE
   private
!
! !PUBLIC MEMBER FUNCTIONS:

!
! !DESCRIPTION: 
!
!  {\tt MAM\_ComponentsData} describes properties (density, hygroscopicity, etc.) 
!  of chem/aerosol components in MAM.
!
!
! !REVISION HISTORY:
!
!  21Nov2011  A. Darmenov  Initial version
!
!EOP
!-------------------------------------------------------------------------

   ! Component Names
   ! ---------------
   character(len=*), public, parameter :: MAM_WATER_COMPONENT_NAME         = 'water'
   character(len=*), public, parameter :: MAM_SULFATE_COMPONENT_NAME       = 'sulfate'
   character(len=*), public, parameter :: MAM_AMMONIUM_COMPONENT_NAME      = 'ammonium'
   character(len=*), public, parameter :: MAM_NITRATE_COMPONENT_NAME       = 'nitrate'
   character(len=*), public, parameter :: MAM_BLACK_CARBON_COMPONENT_NAME  = 'black carbon'
   character(len=*), public, parameter :: MAM_DUST_COMPONENT_NAME          = 'dust'
   character(len=*), public, parameter :: MAM_SEASALT_COMPONENT_NAME       = 'seasalt'
   character(len=*), public, parameter :: MAM_SOA_COMPONENT_NAME           = 'secondary organic'
   character(len=*), public, parameter :: MAM_POM_COMPONENT_NAME           = 'primary organic'

   character(len=*), public, parameter :: MAM_H2O2_COMPONENT_NAME          = 'hydrogen peroxide'
   character(len=*), public, parameter :: MAM_H2SO4_COMPONENT_NAME         = 'sulfuric acid'
   character(len=*), public, parameter :: MAM_SO2_COMPONENT_NAME           = 'sulfur dioxide'
   character(len=*), public, parameter :: MAM_DMS_COMPONENT_NAME           = 'dimethyl sulfide'
   character(len=*), public, parameter :: MAM_NH3_COMPONENT_NAME           = 'ammonia'
   character(len=*), public, parameter :: MAM_SOA_GAS_COMPONENT_NAME       = 'SOA gas'


   ! Component Bulk Density, 'kg m-3'
   ! --------------------------------
   real, public, parameter :: MAM_WATER_COMPONENT_DENSITY        = 1000.0
   real, public, parameter :: MAM_SULFATE_COMPONENT_DENSITY      = 1770.0
   real, public, parameter :: MAM_AMMONIUM_COMPONENT_DENSITY     = 1770.0
   real, public, parameter :: MAM_BLACK_CARBON_COMPONENT_DENSITY = 1700.0
   real, public, parameter :: MAM_DUST_COMPONENT_DENSITY         = 2600.0
   real, public, parameter :: MAM_SEASALT_COMPONENT_DENSITY      = 1900.0
   real, public, parameter :: MAM_SOA_COMPONENT_DENSITY          = 1000.0
   real, public, parameter :: MAM_POM_COMPONENT_DENSITY          = 1000.0


   ! Component Molecular Weight, 'kg Kmol-1'
   ! ---------------------------------------
   real, public, parameter :: MAM_WATER_COMPONENT_MOLECULAR_WEIGHT        = 18.0
   real, public, parameter :: MAM_SULFATE_COMPONENT_MOLECULAR_WEIGHT      = 96.0
   real, public, parameter :: MAM_AMMONIUM_COMPONENT_MOLECULAR_WEIGHT     = 18.0
   real, public, parameter :: MAM_NITRATE_COMPONENT_MOLECULAR_WEIGHT      = 62.0
   real, public, parameter :: MAM_BLACK_CARBON_COMPONENT_MOLECULAR_WEIGHT = 12.0
   real, public, parameter :: MAM_DUST_COMPONENT_MOLECULAR_WEIGHT         = 135.0
   real, public, parameter :: MAM_SEASALT_COMPONENT_MOLECULAR_WEIGHT      = 58.5
   real, public, parameter :: MAM_SOA_COMPONENT_MOLECULAR_WEIGHT          = 12.0
   real, public, parameter :: MAM_POM_COMPONENT_MOLECULAR_WEIGHT          = 12.0

   real, public, parameter :: MAM_H2O2_COMPONENT_MOLECULAR_WEIGHT         = 34.0147
   real, public, parameter :: MAM_H2SO4_COMPONENT_MOLECULAR_WEIGHT        = 98.07848
   real, public, parameter :: MAM_SO2_COMPONENT_MOLECULAR_WEIGHT          = 64.064
   real, public, parameter :: MAM_DMS_COMPONENT_MOLECULAR_WEIGHT          = 62.1324
   real, public, parameter :: MAM_NH3_COMPONENT_MOLECULAR_WEIGHT          = 17.03052
   real, public, parameter :: MAM_SOA_GAS_COMPONENT_MOLECULAR_WEIGHT      = 12.011



   ! Component Hygroscopicity
   ! ------------------------
   real, public, parameter :: MAM_SULFATE_COMPONENT_HYGROSCOPICITY      = 0.507
   real, public, parameter :: MAM_AMMONIUM_COMPONENT_HYGROSCOPICITY     = 0.507
   real, public, parameter :: MAM_BLACK_CARBON_COMPONENT_HYGROSCOPICITY = 1e-10
   real, public, parameter :: MAM_DUST_COMPONENT_HYGROSCOPICITY         = 0.068
   real, public, parameter :: MAM_SEASALT_COMPONENT_HYGROSCOPICITY      = 1.16
   real, public, parameter :: MAM_SOA_COMPONENT_HYGROSCOPICITY          = 0.14
   real, public, parameter :: MAM_POM_COMPONENT_HYGROSCOPICITY          = 0.10


   ! Component Solubility
   ! --------------------
   real, public, parameter :: MAM_SULFATE_COMPONENT_SOLUBILITY      = 1.00
   real, public, parameter :: MAM_AMMONIUM_COMPONENT_SOLUBILITY     = 1.00
   real, public, parameter :: MAM_BLACK_CARBON_COMPONENT_SOLUBILITY = 0.00
   real, public, parameter :: MAM_DUST_COMPONENT_SOLUBILITY         = 0.20
   real, public, parameter :: MAM_SEASALT_COMPONENT_SOLUBILITY      = 1.00
   real, public, parameter :: MAM_SOA_COMPONENT_SOLUBILITY          = 0.05
   real, public, parameter :: MAM_POM_COMPONENT_SOLUBILITY          = 0.05

   end module MAM_ComponentsDataMod
