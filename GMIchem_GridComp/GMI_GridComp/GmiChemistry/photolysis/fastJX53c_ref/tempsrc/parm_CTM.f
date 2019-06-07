c    include 'parm_CTM.f'  for fast-JX code v5.3 (prather 6/05)
c
c     I_ = longitude dim of CTM grid
c     J_ = latitude  dim of CTM grid
c     L_ = altitude(levels) dim of CTM grid
c     LWE_ = altitude(level) dim for trop processes (clouds, rain)
c     JVL_ = vertical(levels) dim for J-values
c     L2_  = 2*L1_ = 2*L_ + 2 = no. levels in the basic Fast-JX grid (mid-level_
c     JVN_ =  no. of J-values
c     W_   = dim = no. of Wavelength bins
c     X_   = dim = no. of X-section data sets (input data)
c     A_   = dim = no. of Aerosol/cloud Mie sets (input data)
c     MX   = no. of aerosol/cloud types supplied from CTM
c     NTR_ = no. of CTM tracers
c     SZAMAX    Solar zenith angle cut-off, above which to skip calculation
c
c-----------------------------------------------------------------------
      integer, parameter ::  I_=128, J_=64, L_=37, LWE_=37  !for EC T42L37
      integer, parameter ::  JVL_=37, JVN_=62, W_=18, X_=64, A_=21
      integer, parameter ::  L1_=L_+1, L2_=2*L_+2 
      integer, parameter ::  MX=4, NTR_=1
      real*8,  parameter ::  SZAMAX=98.0d0
c-----------------------------------------------------------------------

