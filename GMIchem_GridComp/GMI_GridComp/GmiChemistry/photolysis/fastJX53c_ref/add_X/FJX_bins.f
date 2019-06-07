      program FJX_bins
C
C   ******************************************************************
C   *                                                                *
C   *   GENERATING SOLAR FLUX AVERAGED CROSS SECTION IN GIVEN BINS   *
C   *   >>>>>>>from ref code 77-bin to 18-bin fast-JX v5.3 <<<<<<<   *
C   *                                                                *
C   ******************************************************************
C
C-----------------------------------------------------------------------
C
C   Combines the 77-bin std Xsects (including XX SR ODF bins) = 145 total
C      into the fast-J(2) bins (7 at trop wavels, 18 for strat+trop)
C      according to a supplied "FJX_bins.dat" file for the Fast-JX code
C   Output = FJX_05.dat
C
C   history          Standard model:        Xin Zhu
C                    Fast-J & Distribution: Oliver Wild   (9/98)
C                    Fast-J2:               Huisheng Bian (8/02)
C                    new Fast-J(2)          Prather (7/04)
C                    new Fast-JX v-5.3      Prather (6/05)
C
c      open (10, file='fort10.x',status='OLD')
c      open ( 2, file='FJX_bins.dat',status='OLD')
c      open ( 3, file='FJX_spec.dat',status='NEW')


C-----------------------------------------------------------------------
      include 'cmn_bins.f'

      character*8 TITJ(64)
      data TITJ/
     1'J-NO    ','J-O2    ','J-O3    ','J-O3(1D)','J-H2COa ','J-H2COb ',
     2'J-H2O2  ','J-CH3OOH','J-NO2   ','J-NO3   ','J-N2O5  ','J-HNO2  ',
     3'J-HNO3  ','J-HNO4  ','J-ClNO3a','J-ClNO3b','J-Cl2   ','J-HOCl  ',
     4'J-OClO  ','J-Cl2O2 ','J-ClO   ','J-BrO   ','J-BrNO3 ','J-HOBr  ',
     5'J-BrCl  ','J-N2O   ','J-CFCl3 ','J-CF2Cl2','J-F113  ','J-F114  ',
     6'J-F115  ','J-CCl4  ','J-CH3Cl ','J-MeCCL3','J-CH2Cl2','J-CHF2Cl',
     7'J-F123  ','J-F141b ','J-F142b ','J-CH3Br ','J-H1211 ','J-H1301 ',
     8'J-H2402 ','J-CH2Br2','J-CHBr3 ','J-CH3I  ','J-CF3I  ','J-OCS   ',
     9'J-PAN   ','J-CH3NO3','J-ActAld','J-MeVK  ','J-MeAcr ','J-GlyAld',
     A'J-MEKeto','J-EAld  ','J-MGlyxl','J-Glyxla','J-Glyxlb','J-XXXX  ',
     B'J-C3H6O ','J Q2-Ac ','J Q1A-Ac','J Q1B-Ac'/

C---------77-bin spectral data----------UNIT=10----------------------------
      open (10, file='fort10.x',status='OLD')

      do J=1,40
       do K=1,3
        TQQ(K,J) = 0.0
       enddo
      enddo

      read(10,'(A)') TITLE0
        WRITE(6,'(1X,A)') TITLE0
      read(10,'(10x,8i5)') NJVAL,NQQQ,NWWW,NW1,NW2,NWSRB,NSR,NODF

c----NJVAL =no. J values
C----NQQQ = no. additional X-sects:  # J-values = NQQQ+4+additional (p-dep)
C             e.g., 2 additional J's from pressure dep q-yld (MeVK, Acetone)

      read(10,'(a)')
      read(10,'(8e10.3)') (WBIN(IW),IW=1,NWWW+1)
        do IW=1,NWWW
          WL(IW) = 0.5*(WBIN(IW)+WBIN(IW+1))
        enddo
      read(10,'(a)') 
      read(10,'(8e10.3)') (FL(IW),IW=1,NWWW)
      read(10,'(a)') 
      do L = 1,NSR
        read(10,'(6F10.1,I3)')  (ODF(I,L),I=1,NODF), ISR(L)
      enddo
      read(10,'(a20)') TITLEJ(1,1)
      read(10,'(a)') 
          TITLEJ(2,1) = ' '
          TITLEJ(3,1) = ' '
      do L = 1,NSR
        read(10,'(7F10.1)') FNO(L),(QNO(I,L),I=1,NODF)
      enddo

C-------NO X-sections from oscillator strength: 8.85E-13=(pi*e^2)/(m*c^2) (cm)
      do L = 1,NSR
       K = L+NWSRB
       if(FNO(L).GT.0.0) then
         CNO = 8.85E-13*FNO(L) / (1.E7*(1./WBIN(K) - 1./WBIN(K+1)))
         do I = 1,NODF
           QNO(I,L) = CNO*QNO(I,L)/ODF(I,L)
         enddo
       endif
      enddo
C--------O2:  X-sections & S-R band opacity distribution functions
      do K = 1, 3
        read(10,'(a20,f5.0)') TITLEJ(K,2),TQQ(K,2)
        read(10,'(8e10.3)') (QO2(IW,K),IW=1,NWWW)
        read(10,'(a)') 
        do L = 1,NSR
          read(10,'(8e10.3)') (O2X(I,L,K),I=1,NODF)
        enddo
      enddo
      do K = 1, 3
        read(10,'(A20,F5.0)') TITLEJ(K,3),TQQ(K,3)
        read(10,'(8e10.3)') (QO3(IW,K),IW=1,NWWW)
      enddo
      do K = 1, 3
        read(10,'(A20,F5.0)') TITLEJ(K,4),TQQ(K,4)
        read(10,'(8e10.3)') (Q1D(IW,K),IW=1,NWWW)
      enddo
C-------rest of X-sections
      do JQ=5,NQQQ
        read(10,'(a20,f5.0)') TITLEJ(1,JQ),TQQ(1,JQ)
        read(10,'(8e10.3)') (QQQ(IW,1,JQ),IW=1,NWWW)
        read(10,'(a20,f5.0)') TITLEJ(2,JQ),TQQ(2,JQ)
        read(10,'(8e10.3)') (QQQ(IW,2,JQ),IW=1,NWWW)
      enddo

       write(6,'(I5,A10,2x,a20,2x,a20,2x,3F5.0)') 
     &  (JQ,TITJ(JQ),TITLEJ(1,JQ),TITLEJ(2,JQ),
     &   TQQ(1,JQ),TQQ(2,JQ),TQQ(3,JQ),  JQ=1,NJVAL)
       write(6,'(a3,I2,A10,2x,a20,2x,a20,2x,3F5.0)') 
     &  (' xx',JQ,' q-ylds',TITLEJ(1,JQ),TITLEJ(2,JQ),
     &   TQQ(1,JQ),TQQ(2,JQ),TQQ(3,JQ),  JQ=NJVAL+1,NQQQ)

      close (10)

        NSPEC = NQQQ
        NBIN = NWWW

C****************now combine 77-bins+ODFs (145 bins) to get fast-J*********

      call COMB(18)

      stop
      end


      subroutine COMB(NN)
C-----------------------------------------------------------------------
C  Combine x-sections from 'NBIN' to 'NN' bins  (77 to 18 for Fast-J2)
C-----------------------------------------------------------------------
      include 'cmn_bins.f'
      integer   NN,KT,IG,JG,NG(25,25)
      real*8    RAYLAY
c
C-----Read in wavelength bins to combine:
      open ( 2, file='FJX_bins.dat',status='OLD')
      do IG = 1,NN
        read(2,'(25I4)') (NG(JG,IG),JG=1,25)
      enddo
      close(2)
c
      do K=NW1,NW2
C----calculate true Rayleigh scattering
c        POW=-4.0+0.3228-(3.89e-04*WL(K))-(94.26/WL(K))
c        if(WL(K).gt.550.0) POW=-4.04
c        QRL(K) = 4.02e-28*(WL(K)*0.001)**POW

        QRL(K) = RAYLAY(WL(K))

      enddo
c
C----re-distribution of S-R sub-bins into one dimension
      N=0
      do I=1,NSR
        M=N
        do J=1,ISR(I)
          N=M+J
          FFL(N)=FL(I)*ODF(J,I)
C----weight the wavelengths by approx. Rayleigh scattering (1/w^4)
          WWL(N)=FL(I)*ODF(J,I)/WL(I)**4
          QQRL(N)=QRL(I)
          QQNO(N)=QNO(J,I)
          do K=1,3
            QQO2(N,K)=O2X(J,I,K)
            QQO3(N,K)=QO3(I,K)
            QQ1D(N,K)=Q1D(I,K)
          enddo
          do JQ=1,NQQQ
            QQQQ(N,1,JQ)=QQQ(I,1,JQ)
            QQQQ(N,2,JQ)=QQQ(I,2,JQ)
          enddo
        enddo
      enddo
c
C----distribute the rest of 77 bins into one dimension above
      NT = N + NBIN - NSR
      WRITE(6,*) ' TOTAL BINS',NT,NSR
      do I=N+1,NT
        FFL(I)=FL(I-N+NSR)
        WWL(I)=FL(I-N+NSR)/WL(I-N+NSR)**4
        QQRL(I)=QRL(I-N+NSR)
        QQNO(I)=0.0
        do K=1,3
          QQO2(I,K)=QO2((I-N+NSR),K)
          QQO3(I,K)=QO3((I-N+NSR),K)
          QQ1D(I,K)=Q1D((I-N+NSR),K)
        enddo
        do JQ=1,NQQQ
c  NQQQ instead
          QQQQ(I,1,JQ)=QQQ((I-N+NSR),1,JQ)
          QQQQ(I,2,JQ)=QQQ((I-N+NSR),2,JQ)
        enddo
      enddo
C
      do I=1,NN
        SFL(I)=0.D0
        SWL(I)=0.D0
        SQRL(I)=0.D0
        SQNO(I)=0.D0
        do K=1,3
          SQO2(I,K)=0.D0
          SQO3(I,K)=0.D0
          SQ1D(I,K)=0.D0
        enddo
        do JQ=1,NQQQ
          SQQQ(I,1,JQ)=0.D0
          SQQQ(I,2,JQ)=0.D0
        enddo
      enddo
C
C----replace all X-sect_s by X-sect * Flux
      do N=1,NT
        do IG=1,NN
          do JG=1,25
            if(N.EQ.NG(JG,IG)) then
              call SUMM(IG,N)
            endif
          enddo
        enddo
      enddo
C
      do I=1,NN
        FL(I)=SFL(I)
        WL(I)=1.0/SQRT(SQRT(SWL(I)/FL(I)))
        QRL(I)=SQRL(I)/FL(I)
        QQNO(I)=SQNO(I)/FL(I)
        do K=1,3
          QO2(I,K)=SQO2(I,K)/FL(I)
          Q1D(I,K)=SQ1D(I,K)/SQO3(I,K)
          QO3(I,K)=SQO3(I,K)/FL(I)
        enddo
        do JQ=1,NQQQ
          QQQ(I,1,JQ)=SQQQ(I,1,JQ)/FL(I)
          QQQ(I,2,JQ)=SQQQ(I,2,JQ)/FL(I)
        enddo
      enddo

C-----Open output files
      open ( 3, file='FJX_spec.dat',status='NEW')

c  output in Fast-JX 18-wavelength format
c  (JX_spec.dat) UCI fastJX-5.3: JPL02+irHNO4+IUPAC/NO2/VOC+Bltz (J-8.6, 6/05)
      write(3,'(A)') TITLE0
      write(3,'(a10,5I5,a)') 'NW-JValues',NJVAL,NQQQ,NN,1,NN,
     &    '  NJVAL,NQQQ, NWWW, NW1:NW2=01:18 or 12:18'
      write(3,1030) (WL(K),K=NW1,NN)
      write(3,1040) (FL(K),K=NW1,NN)
      write(3,1050) (QRL(K),K=NW1,NN) ! Rough values only
      do KT=1,3
        write(3,1100) TITLEJ(KT,2)(1:7), nint(TQQ(KT,2)),
     &                                    (QO2(K,KT),K=NW1,NN)
      enddo
      do KT=1,3
        write(3,1100) TITLEJ(KT,3)(1:7), nint(TQQ(KT,3)),
     &                                    (QO3(K,KT),K=NW1,NN)
      enddo
      do KT=1,3
        write(3,1100) TITLEJ(KT,4)(1:7), nint(TQQ(KT,4)),
     &                                    (Q1D(K,KT),K=NW1,NN)
      enddo
      write(3,1070) (QQNO(K),K=NW1,NN) 
      write(3,1070) (QQNO(K),K=NW1,NN) 
      do JQ=5,NQQQ
        do I=1,2
          if(TQQ(1,JQ).gt.TQQ(2,JQ)) then
            J=mod(I,2)+1
          else
            J=I
          endif
          write(3,1100)  TITLEJ(J,JQ)(1:7), nint(TQQ(J,JQ)),
     &                                    (QQQ(K,J,JQ),K=NW1,NN)
        enddo
      enddo
c
 1030 format('w-eff (nm)',6f10.0/(10x,6f10.0)/(10x,6f10.0))
 1040 format('SOL#/cm2/s',6(1pE10.3)/(10x,6(1pE10.3))/(10x,6(1pE10.3)))
 1050 format('Raylay cm2',6(1pE10.3)/(10x,6(1pE10.3))/(10x,6(1pE10.3)))
 1070 format('NO     300',6(1pE10.3)/(10x,6(1pE10.3))/(10x,6(1pE10.3)))
 1100 format(a7,i3,6(1pE10.3)/(10x,6(1pE10.3))/(10x,6(1pE10.3)))
c
      close(3)
   99 continue
      return
      end

      
      subroutine SUMM(I,N)
      include 'cmn_bins.f'
         SFL(I)=SFL(I)+FFL(N)
         SWL(I)=SWL(I)+WWL(N)
         SQRL(I)=SQRL(I)+QQRL(N)*FFL(N)
         SQNO(I)=SQNO(I)+QQNO(N)*FFL(N)
        do K=1,3
          SQO2(I,K)=SQO2(I,K)+QQO2(N,K)*FFL(N)
          SQO3(I,K)=SQO3(I,K)+QQO3(N,K)*FFL(N)
          SQ1D(I,K)=SQ1D(I,K)+QQ1D(N,K)*QQO3(N,K)*FFL(N)
        enddo
        do J=1,NQQQ
          SQQQ(I,1,J)=SQQQ(I,1,J)+QQQQ(N,1,J)*FFL(N)
          SQQQ(I,2,J)=SQQQ(I,2,J)+QQQQ(N,2,J)*FFL(N)
        enddo
      return
      end


      function RAYLAY(WAVE)
      real*8 WAVE, RAYLAY
C-----CALCULATE RAYLEIGH CROSS-SECTION AT WAVE (ANGSTROM)
C---RAYLEIGH+RAMAN CROSS-SECTION (INCLUDE FOR ALL WAVELENGTHS)
      if(WAVE.LT.170.) then
          RAYLAY = 1.E-24
      else
       WSQI = 1.E6/(WAVE*WAVE)
       REFRM1 = 1.0E-6*(64.328+29498.1/(146.-WSQI)+255.4/(41.-WSQI))
       RAYLAY = 5.40E-21*(REFRM1*WSQI)**2
      endif
      return
      end

