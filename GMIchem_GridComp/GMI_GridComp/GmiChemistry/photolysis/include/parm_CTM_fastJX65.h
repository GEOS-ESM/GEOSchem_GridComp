!    include 'parm_CTM.h'  for fast-JX code v5.3+ (prather 6/05)
!
!     I_ = longitude dim of CTM grid
!     J_ = latitude  dim of CTM grid
!     L_ = altitude(levels) dim of CTM grid
!     LWE_ = altitude(level) dim for trop processes (clouds, rain)
!     JVL_ = vertical(levels) dim for J-values
!     L2_  = 2*L1_ = 2*L_ + 2 = no. levels in the basic Fast-JX grid (mi
!     JVN_ =  no. of J-values
!     W_   = dim = no. of Wavelength bins   (now in parm_MIE.f)
!     WX_  = dim = no. of wavelengths in input file
!     X_   = dim = no. of X-section data sets (input data)
!     A_   = dim = no. of Aerosol/cloud Mie sets (input data)
!     MX   = no. of aerosol/cloud types supplied from CTM
!     NTR_ = no. of CTM tracers
!     SZAMAX    Solar zenith angle cut-off, above which to skip calculat
!
!-----------------------------------------------------------------------
!Old Setting
!      integer, parameter :: I_ = 128, J_ = 64, L_ = 37, LWE_ = 37  !for EC T42L
!      integer, parameter :: JVL_ = 37, JVN_ = 62, X_ = 64, A_ = 40, WX_ = 18
!      integer, parameter :: L1_ = L_ + 1, L2_ = 2 * L_ + 2  
!      integer, parameter :: MX = 4, NTR_ = 1  
!
!New Setting
      real*8,  parameter :: SZAMAX=98.0d0
      integer, parameter :: MAX_JVL_ = 70
      integer, parameter :: MAX_JVN_ = 82
      integer, parameter :: MAX_L_   = 132
      integer, parameter :: MAX_L1_  = MAX_L_ + 1
      integer, parameter :: MAX_L2_  = 2*MAX_L_ + 2

      integer, parameter :: I_ = 1,  J_   = 1
!      integer, parameter :: X_ = 64, A_   = 40, WX_ = 18
      integer, parameter :: X_ = 63, A_   = 56, WX_ = 18
      integer, parameter :: MX = 45,  NTR_ = 1  
!integer, parameter :: MX = 4,  NTR_ = 1  
      integer            ::  L_, JVL_, JVN_, L1_, L2_, NCTA

      common /dimsFastJX65/ L_, JVL_, JVN_, L1_, L2_, NCTA
