!    include 'parm_mie.f'  for fast-JX code v5.3+ (prather 6/05)
!
!     N_  = no. of levels in Mie scattering arrays
!         = 2*NC+1 = 4*LPAR + 5 + 2*sum(JADDLV)
!     M_  = no. of Gauss points used, must = 4 in fast_JX (no option)
!     M2_ = 2*M_ = 8, replaces MFIT
!     W_   = dim = no. of Wavelength bins:  =18 std, =12 trop only
!
!-----------------------------------------------------------------------
      integer, parameter::N_ = 751, M_ = 4, M2_ = 2 * M_, W_ = 18  
!      integer, parameter ::   N_=501, M_=4, M2_=2*M_, W_=12
!-----------------------------------------------------------------------
!    4 Gauss pts = 8-stream
      real * 8, dimension (M_), parameter :: EMU = [.06943184420297d0, &
             .33000947820757d0, .66999052179243d0, .93056815579703d0]
      real*8, dimension(M_), parameter  :: WT  = [.17392742256873d0, &
     &   .32607257743127d0, .32607257743127d0, .17392742256873d0]
