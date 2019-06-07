!-------------------------------------------------------------------------
!  NASA/GFSC, SIVO, Code 610.3
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  FastJX53cCoreScattering_mod
!
! !INTERFACE:
!
      module FastJX53cCoreScattering_mod
!
      implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
!
      private
      public  :: MIESCT
!
! !DESCRIPTION:
!  Contains core scattering routines.
!
! !AUTHOR:
!  Michael Prather
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
      CONTAINS
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: MIESCT
!
! !INTERFACE:
!
      subroutine MIESCT(FJ,FJT,FJB,POMEGA,FZ,ZTAU,ZFLUX,ZREFL,ZU0,MFIT,ND)
!
! !USES:
      use FastJX53cMIEparameters_mod, only : M_, N_
!
      implicit none
!
! !INPUT PARAMETERS:
      integer, intent(in)  :: MFIT
      integer, intent(in)  :: ND
      real*8 , intent(in)  :: POMEGA(2*M_,N_)
      real*8 , intent(in)  :: FZ(N_)
      real*8 , intent(in)  :: ZTAU(N_)
      real*8 , intent(in)  :: ZREFL, ZU0, ZFLUX
!
! !OUTPUT PARAMETERS:
      real*8 , intent(out) :: FJ(N_), FJT,FJB
!
! !DESCRIPTION:
! \begin{verbatim}
!-----------------------------------------------------------------------
!   This is an adaption of the Prather radiative transfer code, (mjp, 10/95)
!     Prather, 1974, Astrophys. J. 192, 787-792.
!         Sol_n of inhomogeneous Rayleigh scattering atmosphere. 
!         (original Rayleigh w/ polarization)
!     Cochran and Trafton, 1978, Ap.J., 219, 756-762.
!         Raman scattering in the atmospheres of the major planets.
!         (first use of anisotropic code)
!     Jacob, Gottlieb and Prather, 1989, J.Geophys.Res., 94, 12975-13002.
!         Chemistry of a polluted cloudy boundary layer,
!         (documentation of extension to anisotropic scattering)
!
!    takes atmospheric structure and source terms from std J-code
!    ALSO limited to 4 Gauss points, only calculates mean field!
!
!   mean rad. field ONLY (M=1)
! \end{verbatim}
!
! !LOCAL VARIABLES:
      real*8  WT(M_),EMU(M_),PM(M_,2*M_),PM0(2*M_),CMEQ1
      integer I, ID, IM, M, N
!
! !AUTHOR:
! Michael Prather
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC
!-----------------------------------------------------------------------
!---fix scattering to 4 Gauss pts = 8-stream
      call GAUSSP (N,EMU,WT)

!---calc in OPMIE:  ZFLUX = (ZU0*FZ(ND)*ZREFL)/(1.0d0+ZREFL)
      M = 1
      do I = 1,N
        call LEGND0 (EMU(I),PM0,MFIT)
        do IM = M,MFIT
          PM(I,IM) = PM0(IM)
        enddo
      enddo
!
      CMEQ1 = 0.25d0
      call LEGND0 (-ZU0,PM0,MFIT)
      do IM=M,MFIT
        PM0(IM) = CMEQ1*PM0(IM)
      enddo
!
      call BLKSLV(FJ,POMEGA,FZ,ZTAU,ZFLUX,ZREFL,WT,EMU,PM,PM0 &
     &           ,FJT,FJB,M,N,MFIT,ND)
!
      do ID=1,ND,2
        FJ(ID) = 4.0d0*FJ(ID) + FZ(ID)
      enddo

      return
      end subroutine MIESCT
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: BLKSLV
!
! !INTERFACE:
!
      subroutine BLKSLV(FJ,POMEGA,FZ,ZTAU,ZFLUX,ZREFL,WT,EMU,PM,PM0 &
     &                 ,FJTOP,FJBOT, M,N,MFIT,ND)
!
! !USES:
      use FastJX53cMIEparameters_mod, only : M_, N_
!
      implicit none
!
! !INPUT PARAMETERS:
      integer, intent(in)  :: M, N, MFIT, ND
      real*8 , intent(in)  :: POMEGA(2*M_,N_)
      real*8 , intent(in)  :: FZ(N_)
      real*8 , intent(in)  :: ZTAU(N_) 
      real*8 , intent(in)  :: WT(M_)
      real*8 , intent(in)  :: EMU(M_)
      real*8 , intent(in)  :: PM(M_,2*M_)
      real*8 , intent(in)  :: PM0(2*M_)
      real*8 , intent(in)  :: ZFLUX
      real*8 , intent(in)  :: ZREFL
!
! !OUTPUT PARAMETERS:
      real*8 , intent(out) :: FJ(N_),FJTOP,FJBOT
!
! !DESCRIPTION:
!  Solves the block tri-diagonal system:
!  \begin{verbatim}
!               A(I)*X(I-1) + B(I)*X(I) + C(I)*X(I+1) = H(I)
!  \end{verbatim}
!
! !LOCAL VARIABLES:
      real*8, dimension(M_)    :: A, C1, H
      real*8, dimension(M_,M_) :: B, AA, CC
      real*8                      DD(M_,M_,N_), RR(M_,N_)
      real*8  SUMM
      integer I, J, K, ID
!
! !AUTHOR:
! Michael Prather
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC
!
!-----------UPPER BOUNDARY ID=1
      call GEN(POMEGA,FZ,ZTAU,ZFLUX,ZREFL,WT,EMU,PM,PM0 &
     &        ,B,CC,AA,A,H,C1,M,N,MFIT,ND,1)
      call MATIN4 (B)
      do I = 1,N
         RR(I,1) = 0.0d0
        do J = 1,N
          SUMM = 0.0d0
         do K = 1,N
          SUMM = SUMM - B(I,K)*CC(K,J)
         enddo
         DD(I,J,1) = SUMM
         RR(I,1) = RR(I,1) + B(I,J)*H(J)
        enddo
      enddo
!----------CONTINUE THROUGH ALL DEPTH POINTS ID=2 TO ID=ND-1
      do ID = 2,ND-1
        call GEN(POMEGA,FZ,ZTAU,ZFLUX,ZREFL,WT,EMU,PM,PM0 &
     &          ,B,CC,AA,A,H,C1,M,N,MFIT,ND,ID)
        do I = 1,N
          do J = 1,N
          B(I,J) = B(I,J) + A(I)*DD(I,J,ID-1)
          enddo
          H(I) = H(I) - A(I)*RR(I,ID-1)
        enddo
        call MATIN4 (B)
        do I = 1,N
          RR(I,ID) = 0.0d0
          do J = 1,N
          RR(I,ID) = RR(I,ID) + B(I,J)*H(J)
          DD(I,J,ID) = - B(I,J)*C1(J)
          enddo
        enddo
      enddo
!---------FINAL DEPTH POINT: ND
      call GEN(POMEGA,FZ,ZTAU,ZFLUX,ZREFL,WT,EMU,PM,PM0 &
     &        ,B,CC,AA,A,H,C1,M,N,MFIT,ND,ND)
      do I = 1,N
        do J = 1,N
          SUMM = 0.0d0
          do K = 1,N
          SUMM = SUMM + AA(I,K)*DD(K,J,ND-1)
          enddo
        B(I,J) = B(I,J) + SUMM
        H(I) = H(I) - AA(I,J)*RR(J,ND-1)
        enddo
      enddo
      call MATIN4 (B)
      do I = 1,N
        RR(I,ND) = 0.0d0
        do J = 1,N
        RR(I,ND) = RR(I,ND) + B(I,J)*H(J)
        enddo
      enddo
!-----------BACK SOLUTION
      do ID = ND-1,1,-1
       do I = 1,N
        do J = 1,N
         RR(I,ID) = RR(I,ID) + DD(I,J,ID)*RR(J,ID+1)
        enddo
       enddo
      enddo
!----------MEAN J & H
      do ID = 1,ND,2
        FJ(ID) = 0.0d0
       do I = 1,N
        FJ(ID) = FJ(ID) + RR(I,ID)*WT(I)
       enddo
      enddo
      do ID = 2,ND,2
        FJ(ID) = 0.0d0
       do I = 1,N
        FJ(ID) = FJ(ID) + RR(I,ID)*WT(I)*EMU(I)
       enddo
      enddo

        FJTOP = 0.0d0
        FJBOT = 0.0d0
       do I = 1,N
        FJTOP = FJTOP + RR(I,1)*WT(I)*EMU(I)
!       FJBOT = FJBOT + RR(I,1)*WT(I)*EMU(I) !this is not correct, includes up+down
       enddo
        FJTOP = 4.d0*FJTOP

      return
      end subroutine BLKSLV
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: GEN
!
! !INTERFACE:
!
      subroutine GEN(POMEGA,FZ,ZTAU,ZFLUX,ZREFL,WT,EMU,PM,PM0 &
     &              ,B,CC,AA,A,H,C1,M,N,MFIT,ND,ID)
!
! !USES:
      use FastJX53cMIEparameters_mod, only : M_, N_
!
      implicit none
!
! !INPUT PARAMETERS:
      integer, intent(in)  :: M
      integer, intent(in)  :: N
      integer, intent(in)  :: MFIT
      integer, intent(in)  :: ND
      integer, intent(in)  :: ID
      real*8 , intent(in)  :: POMEGA(2*M_,N_)
      real*8 , intent(in)  :: FZ(N_)
      real*8 , intent(in)  :: ZTAU(N_)
      real*8 , intent(in)  :: WT(M_)
      real*8 , intent(in)  :: EMU(M_)
      real*8 , intent(in)  :: PM(M_,2*M_)
      real*8 , intent(in)  :: PM0(2*M_)
      real*8 , intent(in)  :: ZFLUX
      real*8 , intent(in)  :: ZREFL
!
! !OUTPUT PARAMETERS:
      real*8 , intent(out) :: B (M_,M_)
      real*8 , intent(out) :: AA(M_,M_)
      real*8 , intent(out) :: CC(M_,M_)
      real*8 , intent(out) :: A (M_)
      real*8 , intent(out) :: C1(M_)
      real*8 , intent(out) :: H (M_)
!
! !DESCRIPTION:
!  Generates coefficient matrices for the block tri-diagonal system:
! \begin{verbatim}
!               A(I)*X(I-1) + B(I)*X(I) + C(I)*X(I+1) = H(I)
! \end{verbatim}
!
! !LOCAL VARIABLES:
      integer ID0, ID1, IM, I, J, K, MSTART
      real*8  SUM0, SUM1, SUM2, SUM3
      real*8  DELTAU, D1, D2, SURFAC
!
      real*8  S(M_,M_), W(M_,M_), U1(M_,M_), V1(M_)
!
! !AUTHOR:
! Michael Prather
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC
      if (ID.eq.1 .or. ID.eq.ND) then
!---------calculate generic 2nd-order terms for boundaries
       ID0 = ID
       ID1 = ID+1
       if (ID .ge. ND) ID1 = ID-1

       do I = 1,N
          SUM0 = 0.0d0
          SUM1 = 0.0d0
          SUM2 = 0.0d0
          SUM3 = 0.0d0
        do IM = M,MFIT,2
          SUM0 = SUM0 + POMEGA(IM,ID0)*PM(I,IM)*PM0(IM)
          SUM2 = SUM2 + POMEGA(IM,ID1)*PM(I,IM)*PM0(IM)
        enddo
        do IM = M+1,MFIT,2
          SUM1 = SUM1 + POMEGA(IM,ID0)*PM(I,IM)*PM0(IM)
          SUM3 = SUM3 + POMEGA(IM,ID1)*PM(I,IM)*PM0(IM)
        enddo
         H(I) = 0.5d0*(SUM0*FZ(ID0) + SUM2*FZ(ID1))
         A(I) = 0.5d0*(SUM1*FZ(ID0) + SUM3*FZ(ID1))
        do J = 1,I
          SUM0 = 0.0d0
          SUM1 = 0.0d0
          SUM2 = 0.0d0
          SUM3 = 0.0d0
         do IM = M,MFIT,2
          SUM0 = SUM0 + POMEGA(IM,ID0)*PM(I,IM)*PM(J,IM)
          SUM2 = SUM2 + POMEGA(IM,ID1)*PM(I,IM)*PM(J,IM)
         enddo
         do IM = M+1,MFIT,2
          SUM1 = SUM1 + POMEGA(IM,ID0)*PM(I,IM)*PM(J,IM)
          SUM3 = SUM3 + POMEGA(IM,ID1)*PM(I,IM)*PM(J,IM)
         enddo
         S(I,J) = - SUM2*WT(J)
         S(J,I) = - SUM2*WT(I)
         W(I,J) = - SUM1*WT(J)
         W(J,I) = - SUM1*WT(I)
         U1(I,J) = - SUM3*WT(J)
         U1(J,I) = - SUM3*WT(I)
          SUM0 = 0.5d0*(SUM0 + SUM2)
         B(I,J) = - SUM0*WT(J)
         B(J,I) = - SUM0*WT(I)
        enddo
         S(I,I) = S(I,I) + 1.0d0
         W(I,I) = W(I,I) + 1.0d0
         U1(I,I) = U1(I,I) + 1.0d0
         B(I,I) = B(I,I) + 1.0d0
       enddo

       do I = 1,N
         SUM0 = 0.0d0
        do J = 1,N
         SUM0 = SUM0 + S(I,J)*A(J)/EMU(J)
        enddo
        C1(I) = SUM0
       enddo
       do I = 1,N
        do J = 1,N
          SUM0 = 0.0d0
          SUM2 = 0.0d0
         do K = 1,N
          SUM0 = SUM0 + S(J,K)*W(K,I)/EMU(K)
          SUM2 = SUM2 + S(J,K)*U1(K,I)/EMU(K)
         enddo
         A(J) = SUM0
         V1(J) = SUM2
        enddo
        do J = 1,N
         W(J,I) = A(J)
         U1(J,I) = V1(J)
        enddo
       enddo
       if (ID .eq. 1) then
!-------------upper boundary, 2nd-order, C-matrix is full (CC)
        DELTAU = ZTAU(2) - ZTAU(1)
        D2 = 0.25d0*DELTAU
        do I = 1,N
          D1 = EMU(I)/DELTAU
          do J = 1,N
           B(I,J) = B(I,J) + D2*W(I,J)
           CC(I,J) = D2*U1(I,J)
          enddo
          B(I,I) = B(I,I) + D1
          CC(I,I) = CC(I,I) - D1
!         H(I) = H(I) + 2.0d0*D2*C1(I) + D1*SISOTP
          H(I) = H(I) + 2.0d0*D2*C1(I)
          A(I) = 0.0d0
        enddo
       else
!-------------lower boundary, 2nd-order, A-matrix is full (AA)
        DELTAU = ZTAU(ND) - ZTAU(ND-1)
        D2 = 0.25d0*DELTAU
        SURFAC = 4.0d0*ZREFL/(1.0d0 + ZREFL)
        do I = 1,N
          D1 = EMU(I)/DELTAU
          H(I) = H(I) - 2.0d0*D2*C1(I)
           SUM0 = 0.0d0
          do J = 1,N
           SUM0 = SUM0 + W(I,J)
          enddo
           SUM0 = D1 + D2*SUM0
           SUM1 = SURFAC*SUM0
          do J = 1,N
           B(I,J) = B(I,J) + D2*W(I,J) - SUM1*EMU(J)*WT(J)
          enddo
          B(I,I) = B(I,I) + D1
          H(I) = H(I) + SUM0*ZFLUX
          do J = 1,N
           AA(I,J) = - D2*U1(I,J)
          enddo
           AA(I,I) = AA(I,I) + D1
           C1(I) = 0.0d0
        enddo
       endif
!------------intermediate points:  can be even or odd, A & C diagonal
      else
        DELTAU = ZTAU(ID+1) - ZTAU(ID-1)
        MSTART = M + MOD(ID+1,2)
        do I = 1,N
          A(I) = EMU(I)/DELTAU
          C1(I) = -A(I)
           SUM0 = 0.0d0
          do IM = MSTART,MFIT,2
           SUM0 = SUM0 + POMEGA(IM,ID)*PM(I,IM)*PM0(IM)
          enddo
          H(I) = SUM0*FZ(ID)
          do J=1,I
            SUM0 = 0.0d0
           do IM = MSTART,MFIT,2
            SUM0 = SUM0 + POMEGA(IM,ID)*PM(I,IM)*PM(J,IM)
           enddo
            B(I,J) =  - SUM0*WT(J)
            B(J,I) =  - SUM0*WT(I)
          enddo
          B(I,I) = B(I,I) + 1.0d0
        enddo
      endif
      return
      end subroutine GEN
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: LEGND0
!
! !INTERFACE:
!
      subroutine LEGND0 (X,PL,N)
!
      implicit none
!
! !IMPUT PARAMETERS:
      integer, intent(in)  :: N
      real*8 , intent(in)  :: X
!
! !OUTPUT PARAMETERS:
      real*8 , intent(out) :: PL(N)

! !DESCRIPTION:
!  Calculates ORDINARY Legendre fns of X (real). 
!  \begin{verbatim}
!   from P[0] = PL(1) = 1,  P[1] = X, .... P[N-1] = PL(N)
!  \end{verbatim}
!
! !LOCAL VARIABLES:
      integer I
      real*8  DEN
!
! !AUTHOR:
! Michael Prather
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC
!---Always does PL(2) = P[1]
        PL(1) = 1.d0
        PL(2) = X
        do I = 3,N
         DEN = (I-1)
         PL(I) = PL(I-1)*X*(2.d0-1.0/DEN) - PL(I-2)*(1.d0-1.d0/DEN)
        enddo
      return
      end subroutine LEGND0
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: MATIN4
!
! !INTERFACE:
!
      subroutine MATIN4 (A)
!
      implicit none
!
! !INPUT/OUTPUT PARAMETERS:
     real*8, intent(INout)  ::  A(4,4)
!
! !DESCRIPTION:
!  Invert $4 \times 4$ matrix A(4,4) in place with L-U decomposition (mjp, old...)
!
! !AUTHOR:
! Michael Prather
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC
!---SETUP L AND U
      A(2,1) = A(2,1)/A(1,1)
      A(2,2) = A(2,2)-A(2,1)*A(1,2)
      A(2,3) = A(2,3)-A(2,1)*A(1,3)
      A(2,4) = A(2,4)-A(2,1)*A(1,4)
      A(3,1) = A(3,1)/A(1,1)
      A(3,2) = (A(3,2)-A(3,1)*A(1,2))/A(2,2)
      A(3,3) = A(3,3)-A(3,1)*A(1,3)-A(3,2)*A(2,3)
      A(3,4) = A(3,4)-A(3,1)*A(1,4)-A(3,2)*A(2,4)
      A(4,1) = A(4,1)/A(1,1)
      A(4,2) = (A(4,2)-A(4,1)*A(1,2))/A(2,2)
      A(4,3) = (A(4,3)-A(4,1)*A(1,3)-A(4,2)*A(2,3))/A(3,3)
      A(4,4) = A(4,4)-A(4,1)*A(1,4)-A(4,2)*A(2,4)-A(4,3)*A(3,4)
!---INVERT L
      A(4,3) = -A(4,3)
      A(4,2) = -A(4,2)-A(4,3)*A(3,2)
      A(4,1) = -A(4,1)-A(4,2)*A(2,1)-A(4,3)*A(3,1)
      A(3,2) = -A(3,2)
      A(3,1) = -A(3,1)-A(3,2)*A(2,1)
      A(2,1) = -A(2,1)
!---INVERT U
      A(4,4) = 1.d0/A(4,4)
      A(3,4) = -A(3,4)*A(4,4)/A(3,3)
      A(3,3) = 1.d0/A(3,3)
      A(2,4) = -(A(2,3)*A(3,4)+A(2,4)*A(4,4))/A(2,2)
      A(2,3) = -A(2,3)*A(3,3)/A(2,2)
      A(2,2) = 1.d0/A(2,2)
      A(1,4) = -(A(1,2)*A(2,4)+A(1,3)*A(3,4)+A(1,4)*A(4,4))/A(1,1)
      A(1,3) = -(A(1,2)*A(2,3)+A(1,3)*A(3,3))/A(1,1)
      A(1,2) = -A(1,2)*A(2,2)/A(1,1)
      A(1,1) = 1.d0/A(1,1)
!---MULTIPLY (U-INVERSE)*(L-INVERSE)
      A(1,1) = A(1,1)+A(1,2)*A(2,1)+A(1,3)*A(3,1)+A(1,4)*A(4,1)
      A(1,2) = A(1,2)+A(1,3)*A(3,2)+A(1,4)*A(4,2)
      A(1,3) = A(1,3)+A(1,4)*A(4,3)
      A(2,1) = A(2,2)*A(2,1)+A(2,3)*A(3,1)+A(2,4)*A(4,1)
      A(2,2) = A(2,2)+A(2,3)*A(3,2)+A(2,4)*A(4,2)
      A(2,3) = A(2,3)+A(2,4)*A(4,3)
      A(3,1) = A(3,3)*A(3,1)+A(3,4)*A(4,1)
      A(3,2) = A(3,3)*A(3,2)+A(3,4)*A(4,2)
      A(3,3) = A(3,3)+A(3,4)*A(4,3)
      A(4,1) = A(4,4)*A(4,1)
      A(4,2) = A(4,4)*A(4,2)
      A(4,3) = A(4,4)*A(4,3)
      return
      end subroutine MATIN4
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: GAUSSP
!
! !INTERFACE:
!
      subroutine GAUSSP (N,XPT,XWT)
!
      implicit none
!
! !OUTPUT PARAMETERS:
      integer, intent(out) :: N
      real*8 , intent(out) :: XPT(*)
      real*8 , intent(out) :: XWT(*)
!
! !DESCRIPTION:
!  Loads in pre-set Gauss points for 4 angles from 0 to +1 in cos(theta)=mu
!
! !LOCAL VARIABLES:
      real*8   GPT4(4),GWT4(4)
      integer  I
      data GPT4/.06943184420297d0,.33000947820757d0,.66999052179243d0, &
     &          .93056815579703d0/
      data GWT4/.17392742256873d0,.32607257743127d0,.32607257743127d0, &
     &          .17392742256873d0/
!
! !AUTHOR:
! Michael Prather
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC
      N = 4
      do I = 1,N
        XPT(I) = GPT4(I)
        XWT(I) = GWT4(I)
      enddo
      return
      end subroutine GAUSSP
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: EFOLD
!
! !INTERFACE:
!
      subroutine EFOLD (F0, F1, N, F)
!
      implicit none
!
! !INPUT PARAMETERS:
      integer, intent(in)  :: N
      real*8 , intent(in)  :: F0
      real*8 , intent(in)  :: F1
!
! !OUTPUT PARAMETERS:
      real*8 , intent(out) :: F(101)     !F(N+1)  
!
! !DESCRIPTION:
! \begin{verbatim}
!     ***not used in fast-J, part of original scattering code
!---Speciality subroutine for calculating consistent exp(-tau/mu0)
!---  values on the tau grid so that photons are conserved.
!---  ***only works for plane-parallel, NOT psuedo-spherical atmos
!
!---  calculate the e-fold between two boundaries, given the value
!---     at both boundaries F0(x=0) = top, F1(x=1) = bottom.
!---  presume that F(x) proportional to exp[-A*x] for x=0 to x=1
!---          d2F/dx2 = A*A*F  and thus expect F1 = F0 * exp[-A]
!---           alternatively, could define A = ln[F0/F1]
!---  let X = A*x, d2F/dX2 = F
!---  assume equal spacing (not necessary, but makes this easier)
!---      with N-1 intermediate points (and N layers of thickness dX = A/N)
!---
!---  2nd-order finite difference:  (F(i-1) - 2F(i) + F(i+1)) / dX*dX = F(i)
!---      let D = 1 / dX*dX:
!
!  1  |   1        0        0        0        0        0   |    | F0 |
!     |                                                    |    | 0  |
!  2  |  -D      2D+1      -D        0        0        0   |    | 0  |
!     |                                                    |    | 0  |
!  3  |   0       -D      2D+1      -D        0        0   |    | 0  |
!     |                                                    |    | 0  |
!     |   0        0       -D      2D+1      -D        0   |    | 0  |
!     |                                                    |    | 0  |
!  N  |   0        0        0       -D      2D+1      -D   |    | 0  |
!     |                                                    |    | 0  |
! N+1 |   0        0        0        0        0        1   |    | F1 |
!      
!-----------------------------------------------------------------------
!  Advantage of scheme over simple attenuation factor: conserves total
!  number of photons - very useful when using scheme for heating rates.
!  Disadvantage: although reproduces e-folds very well for small flux
!  differences, starts to drift off when many orders of magnitude are
!  involved.
!-----------------------------------------------------------------------
! \end{verbatim}
!
! !LOCAL VARIABLES:
      integer I
      real*8 A,DX,D,DSQ,DDP1, B(101),R(101)
!
! !AUTHOR:
! Michael Prather
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC
!
      if (F0 .eq. 0.d0) then
        do I = 1,N
          F(I)=0.d0
        enddo
        return
      elseif (F1.eq.0.d0) then
        A = log(F0/1.d-250)
      else
        A = log(F0/F1)
      endif
      DX = float(N)/A
      D = DX*DX
      DSQ = D*D
      DDP1 = D+D+1.d0
      B(2) = DDP1
      R(2) = +D*F0
      do I = 3,N
        B(I) = DDP1 - DSQ/B(I-1)
        R(I) = +D*R(I-1)/B(I-1)
      enddo
      F(N+1) = F1
      do I = N,2,-1
        F(I) = (R(I) + D*F(I+1))/B(I)
      enddo
      F(1) = F0
      return
      end subroutine EFOLD
!EOC
!-------------------------------------------------------------------------

      end module FastJX53cCoreScattering_mod

