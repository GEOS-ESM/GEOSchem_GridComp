!    include 'parm_CTM.f'  for fast-JX code v5.3 (prather 6/05)
!
!     I_ = longitude dim of CTM grid
!     J_ = latitude  dim of CTM grid
!     L_ = altitude(levels) dim of CTM grid
!     LWE_ = altitude(level) dim for trop processes (clouds, rain)
!     JVL_ = vertical(levels) dim for J-values
!     L2_  = 2*L1_ = 2*L_ + 2 = no. levels in the basic Fast-JX grid (mid-level_
!     JVN_ =  no. of J-values
!     W_   = dim = no. of Wavelength bins
!     X_   = dim = no. of X-section data sets (input data)
!     A_   = dim = no. of Aerosol/cloud Mie sets (input data)
!     MX   = no. of aerosol/cloud types supplied from CTM
!     NTR_ = no. of CTM tracers
!     SZAMAX    Solar zenith angle cut-off, above which to skip calculation
!
!-----------------------------------------------------------------------
!      integer, parameter ::  I_=128, J_=64, L_=37, LWE_=37  !for EC T42L37
      integer, parameter ::  I_=1, J_=1, MAX_L_=70  !for EC T42L37
      integer, parameter ::  MAX_JVL_=70, MAX_JVN_=82, W_=18, X_=64
!      integer, parameter ::  A_=21
      integer, parameter ::  A_=56
      integer, parameter ::  MAX_L1_=MAX_L_+1, MAX_L2_=2*MAX_L_+2
      integer, parameter ::  MX=4, NTR_=1
      real*8,  parameter ::  SZAMAX=98.0d0
      integer            ::  L_, JVL_, JVN_, L1_, L2_
      integer            ::  NCTA
      common /dims/ L_, JVL_, JVN_, L1_, L2_, NCTA
!-----------------------------------------------------------------------

