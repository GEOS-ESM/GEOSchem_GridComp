      SUBROUTINE GAUSS(N,X,R,W,ZF)
!----------------------------------------------------------------------------------------------------------------------
!@sum    Perform an n-point quadrature from 2n moments to yield  n abscissas and n weights.
!@auth   Susanne Bauer/Doug Wright
!----------------------------------------------------------------------------------------------------------------------
      IMPLICIT NONE

      ! Arguments.

      INTEGER, INTENT(IN   ) :: N               ! number of quadrature points
      REAL(8), INTENT(INOUT) :: X(2*N)          ! moments
      REAL(8), INTENT(  OUT) :: R(N)            ! abscissas
      REAL(8), INTENT(  OUT) :: W(N)            ! weights
      REAL(8), INTENT(  OUT) :: ZF              ! =0.0 successful quadrature, =1.0 failed quadrature

      ! Local variables.

      INTEGER :: IFAILTQL
      REAL(8) :: A(N),B(N),ANU(2*N),AMU0

      ZF = 0.0D+00                              ! successful quadrature
      AMU0 = X(1)                               ! normalizing moment
      ANU(:) = X(:)/AMU0                        ! normalize the moments
      CALL ORTHOG(N,ANU,A,B)
      CALL GAUCOF(N,A,B,AMU0,R,W,IFAILTQL)
      IF(     MINVAL(R(:)) .LT. 0.0D+00 ) THEN  ! failed quadrature
        ZF = 1.0D+00
      ELSEIF( MINVAL(W(:)) .LT. 0.0D+00 ) THEN  ! failed quadrature
        ZF = 1.0D+00
      ELSEIF( IFAILTQL .GT. 0 ) THEN            ! failed quadrature
        ZF = 1.0D+00
      ENDIF
      RETURN
      END


      SUBROUTINE ORTHOG(N,ANU,A,B)

      IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------
!     SEE NUMERICAL RECIPES, W. PRESS ET AL., 2ND EDITION.
!----------------------------------------------------------------------------------------------------------------------
      INTEGER :: L,N,K,NMAX
      PARAMETER (NMAX=30)
      REAL(8) :: A(N),ANU(2*N),B(N),SIG(2*NMAX+1,2*NMAX+1)

      DO 11 L=3,2*N
        SIG(1,L)=0.D+00
11    CONTINUE
      DO 12 L=2,2*N+1
        SIG(2,L)=ANU(L-1)
12    CONTINUE
      A(1)=ANU(2)/ANU(1)
      B(1)=0.D+00
      DO 14 K=3,N+1
        DO 13 L=K,2*N-K+3
          SIG(K,L)=SIG(K-1,L+1)-A(K-2)*SIG(K-1,L)-B(K-2)*SIG(K-2,L)
          IF(SIG(K,K).LE.0.D+00) SIG(K,K) = 1.D-20
13      CONTINUE
        A(K-1)=SIG(K,K+1)/SIG(K,K)-SIG(K-1,K)/SIG(K-1,K-1)
        B(K-1)=SIG(K,K)/SIG(K-1,K-1)
14    CONTINUE
      RETURN
      END SUBROUTINE ORTHOG


      SUBROUTINE GAUCOF(N,A,B,AMU0,X,W,IFAILTQL)

      IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------
!     SEE NUMERICAL RECIPES, W. PRESS ET AL., 2ND EDITION.
!----------------------------------------------------------------------------------------------------------------------
      INTEGER :: I,J,N,NMAX, IFAILTQL
      PARAMETER (NMAX=30)
      REAL(8) :: A(N),B(N),W(N),X(N),Z(NMAX,NMAX),AMU0
      IFAILTQL = 0
      DO 12 I=1,N
        IF(I.NE.1)B(I)=SQRT(B(I))
        DO 11 J=1,N
          IF(I.EQ.J)THEN
            Z(I,J)=1.D+00
          ELSE
            Z(I,J)=0.D+00
          ENDIF
11      CONTINUE
12    CONTINUE
      CALL TQLI(A,B,N,NMAX,Z,IFAILTQL)
      IF(IFAILTQL.GT.0) RETURN
!----------------------------------------------------------------------------------------------------------------------
!     Ordering of the abscissas is usually not needed.
!----------------------------------------------------------------------------------------------------------------------
!     CALL EIGSRT(A,Z,N,NMAX)
!----------------------------------------------------------------------------------------------------------------------
      DO 13 I=1,N
        X(I)=A(I)
        W(I)=AMU0*Z(1,I)**2
        !--------------------------------------------------------------------------------------------------------------
        ! AVOID ZERO WEIGHTS.
        !--------------------------------------------------------------------------------------------------------------
        ! IF(W(I).EQ.0.D+00) W(I) = 1.D-30
        !--------------------------------------------------------------------------------------------------------------
13    CONTINUE
      RETURN
      END SUBROUTINE GAUCOF


      SUBROUTINE TQLI(D,E,N,NP,Z,IFAILTQL)
      IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------
!     SEE NUMERICAL RECIPES, W. PRESS ET AL., 2ND EDITION.
!----------------------------------------------------------------------------------------------------------------------
      INTEGER :: I,N,NP,M,L,K,ITER,IFAILTQL
      REAL(8) :: D(NP),E(NP),Z(NP,NP),DD,G,R,S,C,P,F,B,PYTHAG
      DO 11 I=2,N
        E(I-1)=E(I)
11    CONTINUE
      E(N)=0.D+00
      DO 15 L=1,N
        ITER=0
1       DO 12 M=L,N-1
          DD=ABS(D(M))+ABS(D(M+1))
          IF (ABS(E(M))+DD.EQ.DD) GOTO 2
12      CONTINUE
        M=N
2       IF(M.NE.L)THEN
          IF(ITER.EQ.300) THEN
            IFAILTQL = 1
            RETURN
          ENDIF
          ITER=ITER+1
          G=(D(L+1)-D(L))/(2.D+00*E(L))
          R=PYTHAG(G,1.0D0)
          G=D(M)-D(L)+E(L)/(G+SIGN(R,G))
          S=1.D+00
          C=1.D+00
          P=0.D+00
          DO 14 I=M-1,L,-1
            F=S*E(I)
            B=C*E(I)
            R=PYTHAG(F,G)
            E(I+1)=R
            IF(R.EQ.0.D0)THEN
              D(I+1)=D(I+1)-P
              E(M)=0.D+00
              GOTO 1
            ENDIF
            S=F/R
            C=G/R
            G=D(I+1)-P
            R=(D(I)-G)*S+2.D0*C*B
            P=S*R
            D(I+1)=G+P
            G=C*R-B
C     OMIT LINES FROM HERE ...
            DO 13 K=1,N
              F=Z(K,I+1)
              Z(K,I+1)=S*Z(K,I)+C*F
              Z(K,I)=C*Z(K,I)-S*F
13          CONTINUE
C     ... TO HERE WHEN FINDING ONLY EIGENVALUES.
14        CONTINUE
          D(L)=D(L)-P
          E(L)=G
          E(M)=0.D+00
          GOTO 1
        ENDIF
15    CONTINUE
      RETURN
      END SUBROUTINE TQLI


      DOUBLE PRECISION FUNCTION PYTHAG(A,B)

      IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------
!     SEE NUMERICAL RECIPES, W. PRESS ET AL., 2ND EDITION.
!----------------------------------------------------------------------------------------------------------------------
      REAL(8) :: ABSA,ABSB, A,B,PHYTAG
      ABSA=ABS(A)
      ABSB=ABS(B)
      IF(ABSA.GT.ABSB)THEN
        PYTHAG=ABSA*SQRT(1.D+00+(ABSB/ABSA)**2)
      ELSE
        IF(ABSB.EQ.0.D+00)THEN
          PYTHAG=0.D+00
        ELSE
          PYTHAG=ABSB*SQRT(1.D+00+(ABSA/ABSB)**2)
        ENDIF
      ENDIF
      RETURN
      END


      SUBROUTINE EIGSRT(D,V,N,NP)

      IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------
!     SEE NUMERICAL RECIPES, W. PRESS ET AL., 2ND EDITION.
!----------------------------------------------------------------------------------------------------------------------
      INTEGER :: N,NP,K,J,I
      REAL(8) :: D(NP),V(NP,NP),P
      DO 13 I=1,N-1
        K=I
        P=D(I)
        DO 11 J=I+1,N
          IF(D(J).GE.P)THEN
            K=J
            P=D(J)
          ENDIF
11      CONTINUE
        IF(K.NE.I)THEN
          D(K)=D(I)
          D(I)=P
          DO 12 J=1,N
            P=V(J,I)
            V(J,I)=V(J,K)
            V(J,K)=P
12        CONTINUE
        ENDIF
13    CONTINUE
      RETURN
      END SUBROUTINE EIGSRT


      SUBROUTINE GAUSSINV(N,X,W,U)
!----------------------------------------------------------------------------------------------------------------------
!     DLW: 091206: Computes the moments from the abscissas and weights.
!                  This routine is independent of the units used.
!----------------------------------------------------------------------------------------------------------------------
      IMPLICIT NONE

      ! Arguments.

      INTEGER :: N          ! number of quadrature points
      REAL(8) :: U(2*N)     ! moments
      REAL(8) :: X(N)       ! abscissas
      REAL(8) :: W(N)       ! weights
 
      ! Local variables.

      INTEGER :: I 

      DO I=1, 2*N
        U(I) = SUM( W(:) * ( X(:)**(I-1) ) )
      ENDDO

      RETURN
      END SUBROUTINE GAUSSINV


      SUBROUTINE TEST_QUAD
!----------------------------------------------------------------------------------------------------------------------
!     DLW, 091306: Check quadrature routines. 
!----------------------------------------------------------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER, PARAMETER :: N = 6               ! number of quadrature points
      REAL(8)            :: X(2*N)              ! moments
      REAL(8)            :: R(N)                ! abscissas
      REAL(8)            :: W(N)                ! weights
      REAL(8)            :: ZF                  ! =0.0 successful quadrature, =1.0 failed quadrature
      INTEGER            :: I,K
      REAL(8)            :: N0,DG,SIGMAG,SG 

      ZF = 0.0D+00
      N0 = 1.0D+03
      DG = 0.1D+00
      SIGMAG = 1.6D+00
      SG = EXP( 0.5D+00 * ( LOG(SIGMAG) )**2 )
      DO I=1, 2*N
        K = I-1
        X(I) = N0 * DG**K * SG**(K*K)
      ENDDO
      WRITE(*,90) X(:)
      CALL GAUSS(N,X,R,W,ZF)
      WRITE(*,90) R(:),W(:)
      CALL GAUSSINV(N,R,W,X)
      WRITE(*,90) X(:)

90    FORMAT(50D18.10)
      RETURN
      END SUBROUTINE TEST_QUAD


