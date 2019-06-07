!----jv_mie.h-----COMMON BLOCKS for FAST-J code: 4x4x85 (prather 4/96)
!
!  Parameters
!  ----------
!
!     NL    Maximum number of levels after insertion of extra Mie levels
!     N__   Number of levels in Mie grid: 2*(2*lpar+2+jaddto(1))+3
!     M__   Number of Gauss points used
!
!-----------------------------------------------------------------------
      INTEGER    NL, N__, M__
      PARAMETER (NL=1000, N__=2*NL, M__=4)
      REAL*8 A,B,C1,H,AA,CC,S,W,UX1,VX1,WT,EMU,PM,PM0,POMEGA
      REAL*8 ZTAU,FZ,FJ,DD,RR,ZREFL,ZFLUX,RADIUS,ZU0
      INTEGER ND,N,M,MFIT
      COMMON/MIEBLK/ A(M__),B(M__,M__),C1(M__),H(M__),AA(M__,M__),  &
     &  CC(M__,M__),S(M__,M__),W(M__,M__),UX1(M__,M__),VX1(M__),WT(M__),  &
     &   EMU(M__),PM(M__,2*M__),PM0(2*M__),POMEGA(2*M__,N__),ZTAU(N__),  &
     &   FZ(N__),FJ(N__),DD(M__,M__,N__),RR(M__,N__),  &
     &   ZREFL,ZFLUX,RADIUS,ZU0
      COMMON/MINDEX/ ND,N,M,MFIT
!-----------------------------------------------------------------------
