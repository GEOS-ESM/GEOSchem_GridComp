!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: MAM_ConstituentsDataMod - MAM chem/aerosol constituents
!
! !INTERFACE:
!
   module MAM_ConstituentsDataMod
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
!  {\tt MAM\_ConstituentsData} describes the chem/aerosol constituents in MAM.
!
!
! !REVISION HISTORY:
!
!  21Nov2011  A. Darmenov  Initial version
!
!EOP
!-------------------------------------------------------------------------

   ! Number of Particles Name
   ! ------------------------
   character(len=*), public, parameter :: MAM_NUMBER_PARTICLES_NAME         = 'NUM'

   ! Absorbed Water Name
   ! -------------------
   character(len=*), public, parameter :: MAM_ABSORBED_WATER_NAME           = 'WTR'

   ! Constituent names
   ! -----------------
   character(len=*), public, parameter :: MAM_SULFATE_CONSTITUENT_NAME      = 'SU'
   character(len=*), public, parameter :: MAM_AMMONIUM_CONSTITUENT_NAME     = 'AMM'
   character(len=*), public, parameter :: MAM_BLACK_CARBON_CONSTITUENT_NAME = 'BC'
   character(len=*), public, parameter :: MAM_DUST_CONSTITUENT_NAME         = 'DU'
   character(len=*), public, parameter :: MAM_SEASALT_CONSTITUENT_NAME      = 'SS'
   character(len=*), public, parameter :: MAM_SOA_CONSTITUENT_NAME          = 'SOA'
   character(len=*), public, parameter :: MAM_POM_CONSTITUENT_NAME          = 'POM'

   character(len=*), public, parameter :: MAM_H2SO4_CONSTITUENT_NAME        = 'H2SO4'
   character(len=*), public, parameter :: MAM_SO2_CONSTITUENT_NAME          = 'SO2'
   character(len=*), public, parameter :: MAM_NH3_CONSTITUENT_NAME          = 'NH3'
   character(len=*), public, parameter :: MAM_SOA_GAS_CONSTITUENT_NAME      = 'SOA_GAS'
   

   character(len=*), public, parameter :: MAM_SU_CONSTITUENT_NAME  = MAM_SULFATE_CONSTITUENT_NAME 
   character(len=*), public, parameter :: MAM_AMM_CONSTITUENT_NAME = MAM_AMMONIUM_CONSTITUENT_NAME
   character(len=*), public, parameter :: MAM_BC_CONSTITUENT_NAME  = MAM_BLACK_CARBON_CONSTITUENT_NAME
   character(len=*), public, parameter :: MAM_DU_CONSTITUENT_NAME  = MAM_DUST_CONSTITUENT_NAME
   character(len=*), public, parameter :: MAM_SS_CONSTITUENT_NAME  = MAM_SEASALT_CONSTITUENT_NAME

   end module MAM_ConstituentsDataMod
