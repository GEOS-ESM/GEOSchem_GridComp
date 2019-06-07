!    include 'parm_mie.f'  for fast-JX code v5.3 (prather 6/05)
!
!     N_  = no. of levels in Mie scattering arrays
!         = 2*NC+1 = 4*LPAR + 5 + 2*sum(JADDLV)
!     M_  = no. of Gauss points used, must = 4 in fast_JX (no option)
!
!-----------------------------------------------------------------------
      integer, parameter ::   N_=501, M_=4
!-----------------------------------------------------------------------

