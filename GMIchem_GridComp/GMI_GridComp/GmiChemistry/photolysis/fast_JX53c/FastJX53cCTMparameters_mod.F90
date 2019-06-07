!-----------------------------------------------------------------------
!                     NASA/GSFC - SIVO ASTG Code 610.3                 !
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: FastJX53cCTMparameters_mod
!
! !INTERFACE:
      module FastJX53cCTMparameters_mod
!
      IMPLICIT NONE
!
! !PUBLIC MEMBER FUNCTIONS:
      private
      public  :: SetFastjCTMdims
!
! !PUBLIC DATA MEMBERS:
!
      public  :: I_     , J_    , L_   , LWE_
      public  :: JVL_   , JVN_  , W_   , X_     , A_
      public  :: L1_    , L2_  
      public  :: MX     , NTR_
      public  :: SZAMAX
      
      integer            :: I_    !  longitude dim of CTM grid
      integer            :: J_    !  latitude  dim of CTM grid
      integer            :: L_    !  altitude(levels) dim of CTM grid
      integer            :: JVL_  !  vertical(levels) dim for J-values
      integer            :: JVN_  !  no. of J-values
      integer            :: L1_
      integer            :: L2_   ! 2*L_+2 !  2*L1_ = 2*L_ + 2 = no. levels in the basic Fast-JX grid (mid-level)
      integer, parameter ::  LWE_   = 37     ! for EC T42L37 - altitude(level) dim for trop processes (clouds, rain)
      integer, parameter ::  A_     = 21     !  dim = no. of Aerosol/cloud Mie sets (input data)
      integer, parameter ::  W_     = 18     !  dim = no. of Wavelength bins
      integer, parameter ::  X_     = 64     !  dim = no. of X-section data sets (input data)
      integer, parameter ::  MX     = 4      !  no. of aerosol/cloud types supplied from CTM
      integer, parameter ::  NTR_   = 1      !  no. of CTM tracers
      real*8,  parameter ::  SZAMAX = 98.0d0 !  Solar zenith angle cut-off, above which to skip calculation
!      integer, parameter ::  I_     = 128    !  longitude dim of CTM grid
!      integer, parameter ::  J_     = 64     !  latitude  dim of CTM grid
!      integer, parameter ::  L_     = 37     !  altitude(levels) dim of CTM grid
!      integer, parameter ::  L1_    = L_+1
!      integer, parameter ::  L2_    = 2*L_+2 !  2*L1_ = 2*L_ + 2 = no. levels in the basic Fast-JX grid (mid-level)
!      integer, parameter ::  JVL_   = 37     !  vertical(levels) dim for J-values
!      integer, parameter ::  JVN_   = 62     !  no. of J-values
!
! !DESCRIPTION:
!  Defined CTM parameters.
!
! !AUTHOR:
! Michael Prather
!
! !HISTORY:
!
!EOP
!-----------------------------------------------------------------------
      CONTAINS
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: SetFastjCTMdims
!
! !INTERFACE:
!
      subroutine SetFastjCTMdims(iI_, jJ_, lL_, jJVL_, jJVN_)
!
! !INPUT PARAMETERS:
      integer, intent(in) :: iI_, jJ_, lL_, jJVL_, jJVN_
! !DESCRIPTION:
!  Sets the values of some integer parameters.
!
! !AUTHOR:
!  Jules Kouatchou
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC
      I_   = iI_
      J_   = jJ_
      L_   = lL_
      JVL_ = jJVL_  
      JVN_ = jJVN_  
      L1_  =   L_ + 1
      L2_  = 2*L_ + 2

      return

      end subroutine SetFastjCTMdims
!EOC
!-----------------------------------------------------------------------

      end module FastJX53cCTMparameters_mod

