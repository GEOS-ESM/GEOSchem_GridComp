!-------------------------------------------------------------------------
!  NASA/GFSC, SIVO, Code 610.3
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: FastJX53cCoreFastj_mod
!
! !INTERFACE:
!
      module FastJX53cCoreFastj_mod
!
      implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
!
      private
      public  :: JRATET
      public  :: JP_ATM
      public  :: RD_XXX
      public  :: RD_MIE
      public  :: FLINT
      public  :: SOLARZ
      public  :: SPHERE
      public  :: EXTRAL
      public  :: OPMIE
!
! !DESCRIPTION:
! Contains core fast-J routines.
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
! !IROUTINE: JRATET
!
! !INTERFACE:
!
      subroutine JRATET(PPJ,TTJ,FFF, VALJL)
!
! !USES:
      use FastJX53cJvaluesVars_mod  , only : NJVAL, TQQ, QO2, QO3, Q1D, TITLEJ
      use FastJX53cJvaluesVars_mod  , only : NW1, NW2, QQQ
      use FastJX53cCTMparameters_mod, only : L1_, W_, JVL_, X_

      implicit none
!
! !INPUT PARAMETERS:
      real*8, intent(in)  :: PPJ(L1_+1)        ! pressure profile at edges
      real*8, intent(in)  :: TTJ(L1_)          ! temperatures at mid-level
      real*8, intent(in)  :: FFF(W_,JVL_)      ! mean actinic flux
!
! !OUTPUT PARAMETERS:
      real*8, intent(out) :: VALJL(JVL_,NJVAL) ! JVL_ = no of levels
!
! !DESCRIPTION:
!  Calculate J-value, called by PTOTOJ.
!
! !LOCAL VARIABLES:
      real*8  VALJ(X_)          ! temp for call J's at one L
      real*8  QO2TOT, QO3TOT, QO31DY, QO31D, QQQT, TFACT
      real*8  TT,PP,DD,TT200,TFACA,TFAC0,TFAC1,TFAC2,QQQA,QQ2,QQ1A,QQ1B
      integer J,K,L, IV
!
! !AUTHOR:
!  Michael Prather
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC
      do L = 1,JVL_    ! master loop over layer = L

!---need temperature and density (for some quantum yields):
!---in this case the Pressures PPJ are defined at the boundaries,
!---                Temperatures in the middle of each layer
        TT   = TTJ(L)
        PP  = (PPJ(L)+PPJ(L+1))*0.5d0
          if (L.eq.1) PP = PPJ(1)
        DD = 7.24e18*PP/TT

        do J = 1,NJVAL
          VALJ(J) = 0.d0
        enddo

        do K = NW1,NW2                    ! Using model 'T's here
           QO3TOT = FLINT(TT,TQQ(1,2),TQQ(2,2),TQQ(3,2) &
     &                       ,QO3(K,1),QO3(K,2),QO3(K,3))
           QO2TOT = FLINT(TT,TQQ(1,1),TQQ(2,1),TQQ(3,1) &
     &                       ,QO2(K,1),QO2(K,2),QO2(K,3))
           QO31DY = FLINT(TT,TQQ(1,3),TQQ(2,3),TQQ(3,3) &
     &                       ,Q1D(K,1),Q1D(K,2),Q1D(K,3))
           QO31D  = QO31DY*QO3TOT
          VALJ(1) = VALJ(1) + QO2TOT*FFF(K,L)
          VALJ(2) = VALJ(2) + QO3TOT*FFF(K,L)
          VALJ(3) = VALJ(3) + QO31D*FFF(K,L)
        enddo
        VALJ(2) = VALJ(2) - VALJ(3) ! Kouatchou To see if O3 -> O(3P)+O2 will be close to from previous fastj versions

        do J = 4,NJVAL

          if (TQQ(2,J) .gt. TQQ(1,J)) then
           TFACT = max(0.d0,min(1.d0,(TT-TQQ(1,J))/(TQQ(2,J)-TQQ(1,J))))
          else
           TFACT = 0.d0
          endif

          do K = NW1,NW2
            QQQT    = QQQ(K,1,J) + (QQQ(K,2,J) - QQQ(K,1,J))*TFACT
            VALJ(J) = VALJ(J) + QQQT*FFF(K,L)
          enddo

! #52 Methylvinyl ketone   'MeVK  '     q(M) = 1/(1 + 1.67e-19*[M])
          if (TITLEJ(J).eq.'MeVK  ') then
            VALJ(J) = VALJ(J)/(1.0 + 1.67e-19*DD)
          endif
! #55 Methylethyl ketone   MEKeto     q(M) = 1/(1 + 2.0*[M/2.5e19])
          if (TITLEJ(J).eq.'MEKeto') then
            VALJ(J) = VALJ(J)/(1.0 + 0.80E-19*DD)
          endif
! #57 Methyl glyoxal       MGlyxl     q(M) = 1/(1 + 4.15*[M/2.5E19])
          if (TITLEJ(J).eq.'MGlyxl') then
            VALJ(J) = VALJ(J)/(1.0 + 1.66e-19*DD)
          endif

        enddo

      if (TITLEJ(NJVAL-1).eq.'Acet-a') then
!--------------J-ref v8.3 includes Blitz ACETONE q-yields--------------
!---Acetone is a special case:   (as per Blitz et al GRL, 2004)
!---     61 = NJVAL-1 = J1(acetone-a) ==> CH3CO + CH3
!---     62 = NJVAL   = J2(acetone-b) ==> CH3 + CO + CH3
          VALJ(NJVAL-1) = 0.d0
          VALJ(NJVAL)   = 0.d0
!---IV=NJVAL-1 = Xsect (total abs) for Acetone - pre-calc Temp interp factors
        IV    = NJVAL-1
        TFACA = (TT-TQQ(1,IV))/(TQQ(2,IV)-TQQ(1,IV))
        TFACA = max(0.d0, min(1.d0, TFACA))
!---IV=NJVAL = Q2 for Acetone=>(2), specifically designed for quadratic interp.
!---      but force to Q2=0 by 210K
        IV    = NJVAL
        TFAC0 = ( (TT-TQQ(1,IV))/(TQQ(2,IV)-TQQ(1,IV)) )**2
        if (TT .lt. TQQ(1,IV)) then
          TFAC0 = (TT - 210.d0)/(TQQ(1,IV)-210.d0)
        endif
        TFAC0 = max(0.d0, min(1.d0, TFAC0))
!---IV=NJVAL+1 = Q1A for Acetone => (1), allow full range of T = 200K-300K
        IV    = NJVAL+1
        TT200 = min(300.d0, max(200.d0, TT))
        TFAC1 = (TT200-TQQ(1,IV))/(TQQ(2,IV)-TQQ(1,IV))
!---IV=NJVAL+2 = Q1B for Acetone => (1)
        IV    = NJVAL+2
        TFAC2 = (TT200-TQQ(1,IV))/(TQQ(2,IV)-TQQ(1,IV))

!---now integrate over wavelengths
        do K = NW1,NW2
!---NJVAL-1 = Xsect (total abs) for Acetone
          IV   = NJVAL-1
          QQQA = QQQ(K,1,IV) + (QQQ(K,2,IV)-QQQ(K,1,IV))*TFACA
!---NJVAL   = Q2 for Acetone=>(2), specifically designed for quadratic interp.
          IV   = NJVAL
          QQ2  = QQQ(K,1,IV) + (QQQ(K,2,IV)-QQQ(K,1,IV))*TFAC0
          if (TT .lt. TQQ(1,IV)) then
            QQ2 = QQQ(K,1,IV)*TFAC0
          endif
!---NJVAL+1 = Q1A for Acetone => (1), allow full range of T = 200K-300K
          IV   = NJVAL+1
          QQ1A = QQQ(K,1,IV) + (QQQ(K,2,IV)-QQQ(K,1,IV))*TFAC1
!---NJVAL+2 = Q1B for Acetone => (1)   ! scaled to [M]=2.5e19
          IV   = NJVAL+2
          QQ1B = QQQ(K,1,IV) + (QQQ(K,2,IV)-QQQ(K,1,IV))*TFAC2
          QQ1B = QQ1B*4.d-20
!---J(61)
          VALJ(NJVAL-1) = VALJ(NJVAL-1) &
     &         + FFF(K,L)*QQQA*(1.d0-QQ2)/(QQ1A + QQ1B*DD)
!---J(62)
          VALJ(NJVAL) = VALJ(NJVAL) + FFF(K,L)*QQQA*QQ2

        enddo    !K
!-----------end v-8.3 includes Blitz ACETONE q-yields--------------
      endif

!----Load array of J-values in native order, need to be indexed/scaled
!    by ASAD-related code later: ZPJ(L,JJ) = VALJL(L,JIND(JJ))*JFACTA(JJ)
        do J=1,NJVAL
          VALJL(L,J) = VALJ(J)
        enddo

      enddo    ! master loop over L=1,JVL_
      return
      end subroutine JRATET
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: JP_ATM
!
! !INTERFACE:
!
      subroutine JP_ATM(PPJ,TTJ,DDJ,ZZJ,ZHL,ZZHT,AER,ABX,ADX,JXTRA)
!
! !USES:
!
      use FastJX53cCTMparameters_mod, only : L1_, L2_
!
      implicit none
!
! !INPUT PARAMETERS:
!--------key amtospheric data needed to solve plane-parallel J---------
      real*8, dimension(5,L1_) :: AER
      integer,dimension(5,L1_) :: ADX
      real*8, dimension(L1_)   :: ABX, TTJ,DDJ,ZZJ,ZHL
      real*8, dimension(L1_+1) :: PPJ 
      integer,dimension(L2_+1) :: JXTRA
      real*8                   :: ZZHT
!
! !DESCRIPTION:
!  print out atmosphere used in J-value calc.
!
! !LOCAL VARIABLES:
      integer  I,J,K,L
      real*8   COL(4),COLO2,COLO3,ZKM,DELZ,ZTOP
!
! !AUTHOR:
!  Michael Prather
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC
        write(6,100)
  100 format('   L z(km)     p      T   ','    d(air)   d(O3)', &
     & '  col(O2)  col(O3)   ndx  col(aer/cld)')

          COLO2 = 0.d0
          COLO3 = 0.d0
         do I=1,4
          COL(I) = 0.d0
         enddo
          ZTOP = ZHL(L1_) + ZZHT

        do L = L1_,1,-1

          do I=1,4
          COL(I) = COL(I) + AER(I,L)
          enddo
          COLO2 = COLO2 + DDJ(L)*0.20948d0  
          COLO3 = COLO3 + ZZJ(L)
          DELZ = ZTOP-ZHL(L)
          ZTOP = ZHL(L)
          ZKM = ZHL(L)*1.d-5

        write(6,110) L,ZKM,PPJ(L),TTJ(L),DDJ(L)/DELZ,ZZJ(L)/DELZ, &
     &       COLO2,COLO3, &
     &      (ADX(I,L),COL(I), I=1,4), JXTRA(L+L),JXTRA(L+L-1)
  110 format(1x,i3,0p,f6.2,f10.3,f7.2,1p,4e9.2,4(i3,e9.2),2i3)

        enddo

      return
      end subroutine JP_ATM
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: RD_XXX
!
! !INTERFACE:
!
      subroutine RD_XXX(NJ1,NAMFIL)
!
! !USES:
      use FastJX53cCTMparameters_mod, only : X_
      use FastJX53cJvaluesVars_mod  , only : TQQ   , QQQ   , QRAYL
      use FastJX53cJvaluesVars_mod  , only : QO2   , QO3   , Q1D
      use FastJX53cJvaluesVars_mod  , only : WBIN  , WL    , FL 
      use FastJX53cJvaluesVars_mod  , only : NJVAL , NW1   , NW2
      use FastJX53cJvaluesVars_mod  , only : TITLE0, TITLEJ, TITLEJ2, TITLEJ3
!
      implicit none
!
! !INPUT PARAMETERS:
      integer     , intent(in) :: NJ1
      character(*), intent(in) ::  NAMFIL
!
! !DESCRIPTION:
!  Read in wavelength bins, solar fluxes, Rayleigh parameters, 
!      T-dependent X-sections. 
! \newline
!  {\bf Current code revised to JPL-02 ver 8.5 (5/05)}
! \newline
!
! \begin{verbatim}
!     NAMFIL   Name of spectral data file (j2_spec.dat) >> j2 for fast-J2
!     NJ1      Channel number for reading data file
!
!     NJVAL    Number of species to calculate J-values for
!     NWWW     Number of wavelength bins, from 1:NWWW
!     WBIN     Boundaries of wavelength bins
!     WL       Centres of wavelength bins - 'effective wavelength'
!     FL       Solar flux incident on top of atmosphere (cm-2.s-1)
!     QRAYL    Rayleigh parameters (effective cross-section) (cm2)
!     QO2      O2 cross-sections
!     QO3      O3 cross-sections
!     Q1D      O3 => O(1D) quantum yield
!     TQQ      Temperature for supplied cross sections
!     QQQ      Supplied cross sections in each wavelength bin (cm2)
! \end{verbatim}
!
! !LOCAL VARIABLES:
      integer  I, J, K, IW, NQQQ, NWWW
!
! !AUTHOR:
!  Michael Prather
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC

      TQQ(:,:) = 0.d0

!----------spectral data----set for new format data J-ver8.3------------------
!         note that NJVAL = # J-values, but NQQQ (>NJVAL) = # Xsects read in
!         for 2005a data, NJVAL = 62 (including a spare XXXX) and 
!              NQQQ = 64 so that 4 wavelength datasets read in for acetone
!         note NQQQ is not used outside this subroutine!

      open (NJ1,FILE=NAMFIL,status='old',form='formatted')
      read (NJ1,100) TITLE0
      read (NJ1,101) NJVAL,NQQQ, NWWW,NW1,NW2
      if (NJVAL.gt.X_ .or. NQQQ.gt.X_) then
        write(6,201) NJVAL,X_
        stop
      endif
!      write(6,'(1X,A)') TITLE0
!----J-values:  1=O2, 2=O3P,3=O3D 4=readin Xsects
      read (NJ1,102) (WL(IW),IW=1,NWWW)
      read (NJ1,102) (FL(IW),IW=1,NWWW)
      read (NJ1,102) (QRAYL(IW),IW=1,NWWW)

!---Read O2 X-sects, O3 X-sects, O3=>O(1D) quant yields (each at 3 temps)
      read (NJ1,103) TITLEJ(1),TQQ(1,1), (QO2(IW,1),IW=1,NWWW)
      read (NJ1,103) TITLEJ2,  TQQ(2,1), (QO2(IW,2),IW=1,NWWW)
      read (NJ1,103) TITLEJ3,  TQQ(3,1), (QO2(IW,3),IW=1,NWWW)

      read (NJ1,103) TITLEJ(2),TQQ(1,2), (QO3(IW,1),IW=1,NWWW)
      read (NJ1,103) TITLEJ2,  TQQ(2,2), (QO3(IW,2),IW=1,NWWW)
      read (NJ1,103) TITLEJ3,  TQQ(3,2), (QO3(IW,3),IW=1,NWWW)

      read (NJ1,103) TITLEJ(3),TQQ(1,3), (Q1D(IW,1),IW=1,NWWW)
      read (NJ1,103) TITLEJ2,  TQQ(2,3), (Q1D(IW,2),IW=1,NWWW)
      read (NJ1,103) TITLEJ3,  TQQ(3,3), (Q1D(IW,3),IW=1,NWWW)

!      do J = 1,3
!        write(6,200) TITLEJ(J),(TQQ(I,J),I=1,3)
!      enddo

!---Read remaining species:  X-sections at 2 T_s
      do J = 4,NQQQ
        read (NJ1,103) TITLEJ(J),TQQ(1,J),(QQQ(IW,1,J),IW=1,NWWW)
        read (NJ1,103) TITLEJ2,  TQQ(2,J),(QQQ(IW,2,J),IW=1,NWWW)
!          write(6,200) TITLEJ(J),(TQQ(I,J),I=1,2)
      enddo

!  Reset the titles for NJVAL-1 & NJVAL to be the two acetone J_s
!   61: C3H6O  = Acet-a     (CH3CO + CH3) 
!   62: Q2-Ac  = Acet-b     (CH3 + CO + CH3)

      TITLEJ(NJVAL-1) = 'Acet-a'
      TITLEJ(NJVAL)   = 'Acet-b'
      
      close(NJ1)
      
  100 format(a)
  101 format(10x,5i5)
  102 format(10x,    6e10.3/(10x,6e10.3)/(10x,6e10.3))
  103 format(a7,f3.0,6e10.3/(10x,6e10.3)/(10x,6e10.3))
  200 format(1x,' x-sect:',a10,3(3x,f6.2))
  201 format(' Number of x-sections supplied to Fast-J2: ',i3,/, &
     &       ' Maximum number allowed (X_) only set to: ',i3, &
     &       ' - increase in cmn_jv.f')

      return
      end subroutine RD_XXX
!EOC      
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: RD_MIE
!
! !INTERFACE:
!
      subroutine RD_MIE(NJ1,NAMFIL)
!
! !USES:
!
      use FastJX53cJvaluesVars_mod  , only : NAA, QAA, WAA, PAA, RAA, SSA
      use FastJX53cJvaluesVars_mod  , only : JTAUMX, ATAU, ATAU0
      use FastJX53cJvaluesVars_mod  , only : TITLE0, TITLEA

      implicit none
!
! !INPUT PARAMETERS:
      integer     , intent(in) :: NJ1
      character(*), intent(in) ::  NAMFIL
!
! !DESCRIPTION:
!  Aerosols/cloud scattering data set for fast-JX (ver 5.3).
!  \newline
!  {\bf Spectral data rev to J-ref ver8.5 (5/05)}
!
!  \begin{verbatim}
!     NAMFIL   Name of scattering data file (e.g., FJX_scat.dat)
!     NJ1      Channel number for reading data file
!     NAA      Number of categories for scattering phase functions
!     QAA      Aerosol scattering phase functions
!     NK       Number of wavelengths at which functions supplied (set as 4)
!     WAA      Wavelengths for the NK supplied phase functions
!     PAA      Phase function: first 8 terms of expansion
!     RAA      Effective radius associated with aerosol type
!     SSA      Single scattering albedo
!  \end{verbatim}
!
! !LOCAL VARIABLES:
      integer  I, J, K
!
! !AUTHOR:
!  Michael Prather
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC
      open (NJ1,FILE=NAMFIL,status='old',form='formatted')

      read (NJ1,'(i2,a78)') NAA,TITLE0
      read (NJ1,'(5x,i5,2f10.5)') JTAUMX,ATAU,ATAU0
!        write(6,'(a,2f9.5,i5)') ' ATAU/ATAU0/JMX',ATAU,ATAU0,JTAUMX
      read (NJ1,*)
      
      do J = 1,NAA
          read (NJ1,'(3x,a20)') TITLEA(J)
        do K = 1,4     ! Fix number of aerosol wavelengths at 4 
          read (NJ1,'(f5.0,f8.1,f7.3,f8.4,f7.3,7f6.3)')  &
     &      WAA(K,J),QAA(K,J),RAA(K,J),SSA(K,J),(PAA(I,K,J),I=1,8)
        enddo
      enddo

      close(NJ1)

!        write(6,*) 'Aerosol phase functions & wavelengths'
!        write(6,*) TITLE0
!      do J=1,NAA
!        write(6,'(1x,A8,I2,A,9F8.1)') &
!     &                   TITLEA(J),J,'  wavel=',(WAA(K,J),K=1,4)
!        write(6,'(9x,I2,A,9F8.4)') J,'  Qext =',(QAA(K,J),K=1,4)
!      enddo

      return
      end subroutine RD_MIE
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !IFUNCTION: FLINT
!
! !INTERFACE:
!
      real*8 FUNCTION FLINT (TINT,T1,T2,T3,F1,F2,F3)
!
      implicit none
!
! !INPUT PARAMETERS:
      real*8  TINT,T1,T2,T3,F1,F2,F3
!
! !DESCRIPTION:
!  Three-point linear interpolation function
!
! !AUTHOR:
!  Michael Prather
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC
      if (TINT .le. T2)  then
        if (TINT .le. T1)  then
          FLINT = F1
        else
          FLINT = F1 + (F2 - F1)*(TINT -T1)/(T2 -T1)
        endif
      else
        if (TINT .ge. T3)  then
          FLINT = F3
        else
          FLINT = F2 + (F3 - F2)*(TINT -T2)/(T3 -T2)
        endif
      endif
      return
      end FUNCTION FLINT
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: SOLARZ
!
! !INTERFACE:
!
      subroutine SOLARZ(GMTIME,NDAY,YGRDJ,XGRDI, SZA,COSSZA,SOLFX)
!
      implicit none
!
! !INPUT PARAMETERS:
      real*8, intent(in) ::   GMTIME,YGRDJ,XGRDI, SZA
      integer, intent(in) ::  NDAY
!
! !OUTPUT PARAMETERS:
      real*8, intent(out) ::  COSSZA,SOLFX
!
! !DESCRIPTION:
!  Calculate the Solar Zenith Angle and Solar Flux factor for given 
!  lat/lon/UT.
!  \begin{verbatim}
!     GMTIME = UT for when J-values are wanted 
!           (for implicit solver this is at the end of the time step)
!     NDAY   = integer day of the year (used for solar lat and declin)
!     YGRDJ  = laitude (radians) for grid (I,J)
!     XGDRI  = longitude (radians) for grid (I,J)
!
!     SZA = solar zenith angle in degrees
!     COSSZA = U0 = cos(SZA)
!  \end{verbatim}
!
! !LOCAL VARIABLES:
      real*8  PI, PI180, LOCT
      real*8  SINDEC, SOLDEK, COSDEC, SINLAT, SOLLAT, COSLAT, COSZ
!
! !AUTHOR:
!  Michael Prather
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC
!
      PI     = 3.141592653589793d0
      PI180  = PI/180.d0
      !SINDEC = 0.3978d0*sin(0.9863d0*(dble(NDAY)-80.d0)*PI180)
      !SOLDEK = asin(SINDEC)
      !COSDEC = cos(SOLDEK)
      !SINLAT = sin(YGRDJ)
      !SOLLAT = asin(SINLAT)
      !COSLAT = cos(SOLLAT)
!
      !LOCT   = (((GMTIME)*15.d0)-180.d0)*PI180 + XGRDI
      !COSSZA = COSDEC*COSLAT*cos(LOCT) + SINDEC*SINLAT
      !SZA    = acos(COSSZA)/PI180

      COSSZA = cos(SZA*PI180)

      SOLFX  = 1.d0-(0.034d0*cos(dble(NDAY-186)*2.d0*PI/365.d0))

      return
      end subroutine SOLARZ
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: SPHERE
!
! !INTERFACE:
!
      subroutine SPHERE(GMU,RAD,ZHL,ZZHT,AMF,L1_)
!
      implicit none
!
! !INPUT PARAMETERS:
      integer, intent(in)  :: L1_      ! dimension of CTM = levels +1
      real*8 , intent(in)  :: GMU      ! MU0 = cos(solar zenith angle)
      real*8 , intent(in)  :: RAD      ! radius of Earth mean sea level (cm)
      real*8 , intent(in)  :: ZHL(L1_) ! height (cm) of the bottome edge of CTM level L1_
      real*8 , intent(in)  :: ZZHT     ! scale height (cm) used above top of CTM (ZHL(L_+1)
!
! !OUTPUT PARAMETERS:
      real*8 , intent(out) :: AMF(L1_,L1_) 
                           ! air mass factor for CTM level I for sunlight reaching J
!
! !DESCRIPTION:
!  Calculation of spherical geometry; derive tangent heights, slant path
!  lengths and air mass factor for each layer. Not called when
!  $SZA \;\; > \;\; 98$ degrees.  Beyond 90 degrees, include treatment 
!  of emergent  beam (where tangent height is below altitude J-value 
!  desired at).
!
! !LOCAL VARIABLES:
      integer  I, J, K, II
      real*8   AIRMAS
      real*8   XMU1, XMU2
      real*8   DIFF, Ux
      real*8   H, ZBYR
      real*8   RZ(L1_)   ! Distance from centre of Earth to each point (cm)
      real*8   RQ(L1_)   ! Square of radius ratios
      real*8   TANHT     ! Tangent height for the current SZA
      real*8   XL        ! Slant path between points
!
! !AUTHOR:
!  Michael Prather
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC
!
!  Inlined air mass factor function for top of atmosphere
      AIRMAS(Ux,H) = (1.0d0+H)/SQRT(Ux*Ux+2.0d0*H*(1.0d0- &
     &         0.6817d0*EXP(-57.3d0*abs(Ux)/SQRT(1.0d0+5500.d0*H))/ &
     &                                             (1.0d0+0.625d0*H)))
!
      RZ(1) = RAD + ZHL(1)
      ZBYR  = ZZHT/RAD
      do II = 2,L1_
        RZ(II)   = RAD + ZHL(II)
        RQ(II-1) = (RZ(II-1)/RZ(II))**2
      enddo
      if (GMU .lt. 0.0d0)  then
        TANHT = RZ(1)/dsqrt(1.0d0-GMU**2)
      else
        TANHT = RZ(1)
      endif
!
!  Go up from the surface calculating the slant paths between each level
!  and the level above, and deriving the appropriate Air Mass Factor
      do 16 J = 1,L1_
        do I = 1,L1_
          AMF(I,J) = 0.d0
        enddo
!
!  Air Mass Factors all zero if below the tangent height
        if (RZ(J) .lt. TANHT) goto 16
!  Ascend from layer J calculating AMFs
        XMU1 = abs(GMU)
        do I = J,L1_-1
          XMU2     = dsqrt(1.0d0 - RQ(I)*(1.0d0-XMU1**2))
          XL       = RZ(I+1)*XMU2 - RZ(I)*XMU1
          AMF(I,J) = XL / (RZ(I+1)-RZ(I))
          XMU1     = XMU2
        enddo
!  Use function and scale height to provide AMF above top of model
        AMF(L1_,J)  = AIRMAS(XMU1,ZBYR)
!
!  Twilight case - Emergent Beam
        if (GMU .ge. 0.0d0) goto 16
        XMU1       = abs(GMU)
!  Descend from layer J
        do II = J-1,1,-1
          DIFF        = RZ(II+1)*sqrt(1.0d0-XMU1**2)-RZ(II)
          if (II.eq.1)  DIFF = max(DIFF,0.d0)   ! filter
!  Tangent height below current level - beam passes through twice
          if (DIFF .lt. 0.0d0)  then
            XMU2      = sqrt(1.0d0 - (1.0d0-XMU1**2)/RQ(II))
            XL        = abs(RZ(II+1)*XMU1-RZ(II)*XMU2)
            AMF(II,J) = 2.d0*XL/(RZ(II+1)-RZ(II))
            XMU1      = XMU2
!  Lowest level intersected by emergent beam
          else
            XL        = RZ(II+1)*XMU1*2.0d0
            AMF(II,J) = XL/(RZ(II+1)-RZ(II))
            goto 16
          endif
        enddo

   16 continue
      return
      end subroutine SPHERE
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: EXTRAL
!
! !INTERFACE:
!
      subroutine EXTRAL(AER,ADX,L1X,L2X,NX,JTAUMX,ATAU,ATAU0, JXTRA)
!
      implicit none
!
      integer, intent(in) ::  JTAUMX,L1X,L2X  !index of cloud/aerosol
      integer, intent(in) ::  NX              !Mie scattering array size
      real*8,  intent(in) ::  AER(5,L1X)      !cloud+3aerosol OD in each layer
      real*8,  intent(in) ::  ATAU,ATAU0
      integer, intent(in) ::  ADX(5,L1X)      !index of cloud/aerosol
      integer, intent(out) ::  JXTRA(L2X+1)   !number of sub-layers to be added
!
! !DESCRIPTION:
!    New version 5.3, add sub-layers (JXTRA) to thick cloud/aerosol layers.
!    This version sets up log-spaced sub-layers of increasing thickness ATAU.
!
!   \begin{verbatim}
!     AER(5,L=1:L1X) = Optical Depth in layer L (general visible OD)
!        This is not used in the calculation of J's but in calculating
!        the number in levels to insert in each layer L
!        Set for log-spacing of tau levels, increasing top-down.
!     AER(1:4,L) = aerosol+cloud OD (up to 4 types, index type = ADX(1:4,L)
!     AER(5,L) = Rayleigh-not used here
!   \end{verbatim}
!
!     N.B. the TTAU, etc caluclate here are not really used elsewhere.
!     \newline \newline
!
!   The log-spacing parameters have been tested for convergence and chosen
!   to be within $0.5\%$ for ranges OD=1-500, rflect=$0-100\%$, mu0=0.1-1.0
!   use of ATAU = 1.18 and min = 0.01, gives at most +135 pts for OD=100 
!
! !LOCAL VARIABLES:
      integer JTOTL,I,L,L2 
      real*8  DTAUX(L1X),TTAU(L2X+1),DTAUJ, ATAU1,ATAULN,ATAUM,ATAUN1
!
! !AUTHOR:
!  Michael Prather
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC
!
!---Reinitialize arrays
      TTAU(:)  = 0.d0
      JXTRA(:) = 0
!
!---Set up total optical depth over each CTM level (L=1:L1X)
!---   DTAUX(L) = bulk properties (OD) of each CTM layer

      do L = 1,L1X
!  Total optical depth from all elements I=1:4 = clouds  + 3 aerosol types
!    Assume here that AER(1:4, 1:L1X) is 'visible' Optical Depth = 600 nm
!    do NOT count if ADX = 0
          DTAUX(L)   = 0.d0
        do I = 1,4
          if (ADX(I,L) .gt. 0) DTAUX(L) = DTAUX(L) + AER(I,L)
        enddo
      enddo
!
!---combine these edge- and mid-layer points into grid of size:
!---              L2X+1 = 2*L1X+1 = 2*L_+3
!---calculate column optical depths above each level, TTAU(1:L2X+1)
!---      note that TTAU(L2X+1)=0 and TTAU(1)=total OD
!
!---Divide thick layers to achieve better accuracy in the scattering code
!---In the original fast-J, equal sub-layers were chosen, this is wasteful
!---and this new code (ver 5.3) uses log-scale:  
!---        Each succesive layer (down) increase thickness by ATAU > 1
!---        e.g., if ATAU = 2, a layer with OD = 15 could be divided into
!---        4 sub-layers with ODs = 1 - 2 - 4 - 8
!---The key parameters are:
!---        ATAU = factor increase from one layer to the next
!---        ATAUMN = the smallest OD layer desired
!---        JTAUMX = maximum number of divisions (i.e., may not get to ATAUMN)
!---These are hardwired below, can be changed, but have been tested/optimized

      ATAU1  = ATAU - 1.d0
      ATAULN = log(ATAU)
        TTAU(L2X+1)  = 0.0d0
      do L2 = L2X,1,-1
        L         = (L2+1)/2
        DTAUJ     = 0.5d0 * DTAUX(L)
        TTAU(L2)  = TTAU(L2+1) + DTAUJ
!---Now compute the number of log-spaced sub-layers to be added in
!---   the interval TTAU(L2) > TTAU(L2+1)
!---The objective is to have successive TAU-layers increasing by factor ATAU >1
!---the number of sub-layers + 1
        if (TTAU(L2) .lt. ATAU0) then
          JXTRA(L2) = 0
        else
          ATAUM    = max(ATAU0, TTAU(L2+1))
          ATAUN1 = log(TTAU(L2)/ATAUM) / ATAULN
          JXTRA(L2) = min(JTAUMX, max(0, int(ATAUN1 - 0.5d0)))
        endif
      enddo

!---check on overflow of arrays, cut off JXTRA at lower L if too many levels
      JTOTL    = L2X + 2
      do L2 = L2X,1,-1
        JTOTL  = JTOTL + JXTRA(L2)
        if (JTOTL .gt. NX/2)  then
          write(6,'(A,2I5,F9.2)') 'N_/L2_/L2-cutoff JXTRA:',NX,L2X,L2
          do L = L2,1,-1
            JXTRA(L) = 0
          enddo
          go to 10
        endif
      enddo
  10  continue

      return
      end subroutine EXTRAL
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: OPMIE
!
! !INTERFACE:
!
      subroutine OPMIE &
     &   (KW,KM,WAVEL,ABX,AER,ADX,U0,RFLECT,AMF,JXTRA,FMEAN, FJTOP,FJBOT)
!
! !USES:
!
      use FastJX53cJvaluesVars_mod   , only : QAA, SSA, PAA, ATAU0
      use FastJX53cCTMparameters_mod , only : L1_, L2_, L_
      use FastJX53cMIEparameters_mod , only : M_, N_
      use FastJX53cCoreScattering_mod, only : MIESCT
!
      implicit none
!
!
! !INPUT PARAMETERS:
      real*8 , intent(in)  ::  WAVEL
               ! wavelength of bin (in nm) - not now used
      real*8 , intent(in)  ::  AER(5,L1_)
               ! 5 vertical profiles of Optical Depth in each layer
      real*8 , intent(in)  ::  ABX(L1_)
               ! vertical profiles of ABSORPTION Optical Depth in each layer
               ! includes O2 and O3 for now (BC under aerosols)
      real*8 , intent(in)  ::  AMF(L1_,L1_)
      real*8 , intent(in)  ::  U0
               ! cos (SZA)
      real*8 , intent(in)  ::  RFLECT
               ! Lambertian albedo of surface
      integer, intent(in)  ::  KW           
               ! wavelength bin # (1:18)
      integer, intent(in)  ::  KM           
               ! wavelength index for Mie properties (1:4 = 300-400-600-999 nm)
      integer, intent(in)  ::  ADX(5,L1_)
               ! integer index of the scattering properties of each AER
               ! 1:4 are reserved for aerosols and clouds
               ! 5 is meant only for Rayleigh scattering!
      integer, intent(in)  ::  JXTRA(L2_+1)
               ! number 0:J = no. of additional levels to be inserted
!
! !OUTPUT PARAMETERS:
      real*8 , intent(out) ::  FMEAN(L_)
               ! mean actinic flux at standard CTM levels (mid-layer)
      real*8 , intent(out) ::  FJTOP,FJBOT
! !DESCRIPTION:
! \begin{verbatim}
!
!     DTAUX    Local optical depth of each CTM level
!     TTAU     Optical depth of air vertically above each point (to top of atm)
!     FTAU     Attenuation of solar beam
!     POMEGA   Scattering phase function
!
!---new Ver 5.3 code adds sub-layers (# = JXTRA(L2)) using ATAU as the 
!   factor increase from sub-layer to sub-layer 

!  fast-J Mie code for J_s, only uses 8-term expansion, 4-Gauss pts
!  Currently allow up to A_ aerosol phase functions (at all altitudes) to
!  be associated with optical depth AER(1:L2_) = aerosol opt.depth @ 1000 nm
!  
!  Pick Mie-wavelength with phase function and Qext: e.g., 
!  01 RAYLE = Rayleigh phase
!  02 ISOTR = isotropic
!  03 ABSRB = fully absorbing 'soot', wavelength indep.
!  04 X_Bkg = backgrnd stratospheric sulfate (n=1.46,log-norm:r=.09um/sigma=.6)
!  05 X_Vol = volcanic stratospheric sulfate (n=1.46,log-norm:r=.08um/sigma=.8)
!. . .
!  11 W_C13 = water cloud (C1/Deirm.) (n=1.335, gamma:  r-mode=13.3um /alpha=6)
!. . .
!  13 Ice-H = hexagonal ice cloud (Mishchenko)
!  14 Ice-I = irregular ice cloud (Mishchenko)
!
!  Choice of the 4 aerosol/cloud indices ADX is made earlier
!  Optical depths for the 4 (aerosol+clouds) = AER
!  \end{verbatim}
!
! !LOCAL VARIABLES:
      integer JNDLEV(L_),JADDLV(L2_+1),JADDTO(L2_+1),L2LEV(L2_+1) &
     &       ,JTOTL,I,II,J,K,L,IX,JK,   L2,L2L,L22,LZ,LZZ,NDZ
!
      real*8  QXMIE(5),XLAER(5),SSALB(5),DTAUX(L1_),PIAER(5,L1_) &
     &       ,POMEGAJ(2*M_,L2_+1),TTAU(L2_+1),FTAU(L1_) &
     &       ,FTAU2(L2_+1),DTTAU,DFTAU,dPOMEGA(2*M_) &
     &       ,XLO2,XLO3,XLRAY,XLTAU,ZK,TAUDN,TAUUP,ZK2 &
     &       ,DTAUJ,DTAUAD,TAUBTM,TAUTOP,FBTM,FTOP,POMEGAB(2*M_) &
     &       ,ATAUA,ATAUZ
!--- variables used in mie code-----------------------------------------
      real*8   FJ(N_),POMEGA(2*M_,N_),FZ(N_),ZTAU(N_),ZREFL,ZU0,ZFLUX
      real*8   FJT,FJB
      integer  MFIT, ND
!
! !AUTHOR:
!  Michael Prather
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC
!---Reinitialize arrays
      ZTAU(:)     = 0.d0
      FZ(:)       = 0.d0
      POMEGA(:,:) = 0.d0


!---Set up optical depth DTAUX(L), and scattering fraction PIAER(1:5,L)
!---    where L = 1:L1_ = bulk properties of each CTM layer.
      do L = 1,L1_

        do I = 1,4
          if (ADX(I,L) .eq. 0)  then
            QXMIE(I) = 0.d0
            SSALB(I) = 0.d0
          else
!---for Mie code scale extinction at 600 nm to wavelength WAVEL (QXMIE)
            QXMIE(I) = QAA(KM,ADX(I,L))/QAA(3,ADX(I,L))
            SSALB(I) = SSA(KM,ADX(I,L))
          endif
        enddo
!---special case for Rayleigh scattering
            QXMIE(5) = 1.d0
            SSALB(5) = 1.d0
        do I = 1,5
          XLAER(I) = AER(I,L)*QXMIE(I)
        enddo

          DTAUX(L) = ABX(L)
        do I = 1,5
          DTAUX(L) = DTAUX(L) + XLAER(I)
        enddo
!---fractional extinction for Rayleigh scattering and each aerosol type
        do I = 1,5
          PIAER(I,L) = SSALB(I)*XLAER(I)/DTAUX(L)
        enddo

      enddo

!---Calculate attenuated incident beam exp(-TTAU/U0 = DTAUX * AirMassFactor)
!---      at the edges of the CTM layers L=1:L1_
      do L = 1,L1_
        if (AMF(L,L) .gt. 0.0d0) then
            XLTAU = 0.0d0
          do I = 1,L1_
            XLTAU = XLTAU + DTAUX(I)*AMF(I,L)
          enddo
          if (XLTAU .gt. 76.d0) then   ! zero out flux at 1e-33
            FTAU(L) = 0.d0
          else
            FTAU(L) = exp(-XLTAU)
          endif
        else
          FTAU(L)   = 0.0d0
        endif
      enddo
!
!---Define the total scattering phase fn for each CTM layer L=1:L_+1
!---   from a DTAUX-wt_d mix of aerosols, cloud & Rayleigh
!---No. of quadrature pts fixed at 4(M_), expansion of phase fn @ 8
      MFIT = 2*M_
      do L = 1,L1_
        do I = 1,MFIT
           POMEGAJ(I,L) = 0.d0
          do K = 1,5
           POMEGAJ(I,L)=POMEGAJ(I,L) + PIAER(K,L)*PAA(I,KM,ADX(K,L))
          enddo
        enddo
      enddo

!--printout diagnositcs of different aerosol+cld contrib to OD
!      if(KW.eq.18) then
!      do L=1,L1_
!      write(6,'(3i5,1p,20e12.5)') L,JXTRA(L+L),JXTRA(L+L-1),
!     &   DTAUX(L),(AER(I,L),I=1,5)
!      enddo
!      endif

!
!------------------------------------------------------------------------
!  Take optical properties on CTM layers and convert to a photolysis
!  level grid corresponding to layer centres and boundaries. This is
!  required so that J-values can be calculated for the centre of CTM
!  layers; the index of these layers is kept in the JNDLEV array.
!------------------------------------------------------------------------
!
!---Now combine the CTM layer edges (1:L_+2) with the CTM mid-layer
!---    points (1:L_) plus 1 for the mid point of added top layer.

!---combine these edge- and mid-layer points into grid of size:
!---              L2_+1 = 2*L1_+1 = 2*L_+3
!---calculate column optical depths above each level, TTAU(1:L2_+1)
!---      note that TTAU(L2_+1)=0 and TTAU(1)=total OD
        TTAU(L2_+1) = 0.0d0
      do L2 = L2_,1,-1
        L          = (L2+1)/2
        DTAUJ      = 0.5d0 * DTAUX(L)
        TTAU(L2)   = TTAU(L2+1) + DTAUJ
      enddo
!
!---reflected flux from surface
      if (U0 .gt. 0.d0) then
        ZFLUX = U0*FTAU(1)*RFLECT/(1.d0+RFLECT)
      else
        ZFLUX = 0.d0
      endif
!
!---calculate attenuated beam FTAU2 on the new doubled-levels L2=1:L2_+1
!---       calculate FTAU2 at CTM mid-layers from sqrt
!---       L2_ = 2*L1_ = 2*L_+2

        FTAU2(L2_+1) = 1.0d0
        FTAU2(L2_)   = sqrt(1.0d0*FTAU(L1_))
        FTAU2(L2_-1) = FTAU(L1_)
      do L2 = L2_-3,1,-2
        L           = (L2+1)/2
        FTAU2(L2)   = FTAU(L)
        FTAU2(L2+1) = sqrt(FTAU(L+1)*FTAU(L))
      enddo
!
!  Calculate scattering properties, level centres then level boundaries
! ***be careful of order, we are shifting the 'POMEGAJ' upward in index***
      do L2 = L2_,2,-2
        L   = L2/2
        do I = 1,MFIT
          POMEGAJ(I,L2) = POMEGAJ(I,L)
        enddo
      enddo
!---lower boundary value is set (POMEGAJ(I,1), but set upper:
        do I = 1,MFIT
          POMEGAJ(I,L2_+1) = POMEGAJ(I,L2_)
        enddo
!---now have POMEGAJ filled at even points from L2=3:L2_-1
!---use inverse interpolation for correct tau-weighted values at edges
      do L2 = 3,L2_-1,2
        TAUDN = TTAU(L2-1)-TTAU(L2)
        TAUUP = TTAU(L2)-TTAU(L2+1)
        do I = 1,MFIT
          POMEGAJ(I,L2) = (POMEGAJ(I,L2-1)*TAUDN +  &
     &             POMEGAJ(I,L2+1)*TAUUP) / (TAUDN+TAUUP)
        enddo
      enddo
!---at this point FTAU2(1:L2_+1) and POMEAGJ(1:8, 1:L2_+1)
!---    where FTAU2(L2_+1) = 1.0 = top-of-atmos, FTAU2(1) = surface
!
!------------------------------------------------------------------------
!  Calculate cumulative total and define levels we want J-values at.
!  Sum upwards for levels, and then downwards for Mie code readjustments.
!
!     JXTRA(L2)  Number of new levels to add between (L2) and (L2+1)
!           ***JXTRA(1:L2_+1) is calculated based on the aerosol+cloud OD_s
!     JADDLV(L2)  Number of new levels actually added at each wavelength
!            where JADDLV = 0 when there is effectively no FTAU2 
!     JADDTO(L2)   Total number of new levels to add to and above level (L2)
!     JNDLEV(L) = L2 index that maps on CTM mid-layer L
!
!------------------------------------------------------------------------

!---JADDLV(L2=1:L2_) = number of levels to add between TTAU2(L2) and TTAU(L2+1)
!---    JADDLV is taken from JXTRA, which is based on visible OD.
!---    JADDTO(L2=1:L2_+1) is the cumulative number of levels to be added
!---note that JADDLV and JADDTO will change with wavelength and solar zenith

!--now try to insert additional levels for thick clouds, ONLY IF FTAU2 > 1.e-8
!-- this will cut off additional levels where the solar beam is negligible.

      do L2 = 1,L2_
        if (FTAU2(L2+1) .lt. 1.d-30) then
          JADDLV(L2) = 0
        else
          JADDLV(L2) = JXTRA(L2)
        endif
      enddo
        JADDTO(L2_+1) = 0
      do L2 = L2_,1,-1
        JADDTO(L2) = JADDTO(L2+1) + JADDLV(L2)
      enddo

!---expanded grid now included CTM edge and mid layers plus expanded 
!---    grid to allow for finer delta-tau at tops of clouds.
!---    DIM of new grid = L2_ + JADDTO(1) + 1

!---L2LEV(L2) = L2-index for old level L2 in expanded J-grid (w/JADDLV)
!     in absence of JADDLV, L2LEV(L2) = L2
        L2LEV(1)  = 1
      do L2 = 2,L2_+1
        L2LEV(L2) = L2LEV(L2-1) + 1 + JADDLV(L2-1)
      enddo

!---JNDLEV(L=1:L_) = L2-index in expanded grid for CTM mid-layer L 
      do L = 1,L_
        JNDLEV(L) = L2LEV(2*L)
      enddo

!---------------------SET UP FOR MIE CODE-------------------------------
!
!  Transpose the ascending TTAU grid to a descending ZTAU grid.
!  Double the resolution - TTAU points become the odd points on the
!  ZTAU grid, even points needed for asymm phase fn soln, contain 'h'.
!  Odd point added at top of grid for unattenuated beam   (Z='inf')
!
!  The following mapping holds for JADDLV=0
!        Surface:   TTAU(1)    ==> ZTAU(2*L2_+1)
!        Top:       TTAU(L2_)  ==> ZTAU(3)
!        Infinity:     0.0     ==> ZTAU(1)
!        index: 2*(L2_+1-L2)+1 ==> LZ
!
!  Mie scattering code only used from surface to level L2_
!------------------------------------------------------------------------
!
!------------------------------------------------------------------------
!  Insert new levels, working downwards from the top of the atmosphere
!  to the surface (down in 'LZ', up in 'L2'). This allows ztau and pomega
!  to be incremented linearly (in a +ve sense), and the flux fz to be
!  attenuated top-down (avoiding problems where lower level fluxes are
!  zero).
!
!    zk        fractional increment in level
!    dTTAU     change in ttau per increment    (linear, positive)
!    dPOMEGA   change in pomega per increment  (linear)
!    ftaulog   change in ftau per increment    (exponential, normally < 1)
!
!------------------------------------------------------------------------
!
!  Ascend through atmosphere transposing grid and adding extra points
!  remember L2=1 is surface of CTM, but last layer (LZ) in scattering code.
!  there are twice the number of layers in the LZ arrays (2*L2_ + 2*JADDTO + 1)
!    because we need to insert the intermediate layers (even LZ) for the 
!    asymmetric scattering code.


!  Transfer the L2=1:L2_+1 values (TTAU,FTAU2,POMEGAJ) onto the reverse
!    order, expanded, doubled-level scatter grid. 
!    Note that we need to deal with the expansion by JADD levels (L2L).
!      These JADDLV levels are skipped and need to be interpolated later.
!    Note that only odd LZ levels are filled, 

      NDZ = 2*L2_ + 2*JADDTO(1) + 1

!   Note that the successive sub-layers have the ratio in OD of ATAU
!      ATAUA = (ATAU - 1.d0)/ATAU     ! this is the limit for L22=>inf

      do L2 = 1,L2_+1          ! L2 = index of CTM edge- and mid-layers
        L2L = L2LEV(L2)        ! L2L = index for L2 in expanded scale(JADD)
        LZ  = NDZ + 2 - 2*L2L  ! LZ = index for L2 in scatt arrays
          ZTAU(LZ) = TTAU(L2)
          FZ(LZ)   = FTAU2(L2)
        do I=1,MFIT
          POMEGA(I,LZ) = POMEGAJ(I,L2)
        enddo
      enddo

!   Now go thru the pairs of L2 levels to see if we need JADD levels
      do L2 = 1,L2_             ! L2 = index of CTM edge- and mid-layers
        L2L = L2LEV(L2)         ! L2L = index for L2 in expanded scale(JADD)
        LZ  = NDZ + 2 - 2*L2L   ! LZ = index for L2 in scatt arrays
        L22 = L2LEV(L2+1) - L2LEV(L2) - 1   ! L22 = 0 if no added levels
       if (L22 .gt. 0) then
          TAUBTM = TTAU(L2)
          TAUTOP = TTAU(L2+1)
          FBTM   = FTAU2(L2)
          FTOP   = FTAU2(L2+1)
         do I = 1,MFIT
          POMEGAB(I) = POMEGAJ(I,L2)
         enddo

!---to fit L22 new layers between TAUBOT > TAUTOP, calculate new 1/ATAU factor
!---  such that TAU(just above TAU-btm) = ATUAZ * TAUBTM < TAUBTM

          ATAUZ = exp(-log(TAUBTM/max(TAUTOP,ATAU0))/float(L22+1))

        do L = 1,L22           ! add odd levels between L2LEV(L2) & L2LEV(L2+1)
          LZZ = LZ - 2*L       ! LZZ = index(odd) of added level in scatt arrays
          ZTAU(LZZ) = TAUBTM * ATAUZ
          ATAUA=(TAUBTM-ZTAU(LZZ))/(TAUBTM-TAUTOP) !fraction from TAUBTM=>TAUTOP
          if (FBTM .lt. 1.d-32) then
            FZ(LZZ) = 0.d0
          else    
            FZ(LZZ) = FBTM * (FTOP/FBTM)**ATAUA
          endif
          do I = 1,MFIT
            POMEGA(I,LZZ) = POMEGAB(I) +  &
     &               ATAUA*(POMEGAJ(I,L2+1)-POMEGAB(I))
          enddo
            TAUBTM       = ZTAU(LZZ)
            FBTM         = FZ(LZZ)
          do I = 1,MFIT
            POMEGAB(I) = POMEGA(I,LZZ)
          enddo
        enddo
       endif
      enddo

!   Now fill in the even points with simple interpolation in scatter arrays:
      do LZ = 2,NDZ-1,2
        ZTAU(LZ) = 0.5d0*(ZTAU(LZ-1)+ZTAU(LZ+1))
        FZ(LZ)   = sqrt(FZ(LZ-1)*FZ(LZ+1))
       do I=1,MFIT
        POMEGA(I,LZ) = 0.5d0*(POMEGA(I,LZ-1)+POMEGA(I,LZ+1))
       enddo
      enddo
      
      ND = NDZ
      ZU0 = U0
      ZREFL = RFLECT

!---PRINT diagnostics
!      write(6,'(A,2I6)') 'Total levels of photolysis in OPMIE',ND
      if(ND .gt. N_) then
        write(6,'(a,2i9)') ' overflow of scatter arrays:',ND,N_
        stop
      endif

!-----------------------------------------------------------------------
      call MIESCT(FJ,FJT,FJB,POMEGA,FZ,ZTAU,ZFLUX,ZREFL,ZU0,MFIT,ND)
!-----------------------------------------------------------------------

!---Move mean intensity from scatter array FJ(LZ=1:ND) 
!---              to CTM mid-level array FMEAN(L=1:L_)
      do L = 1,L_
        L2L = JNDLEV(L)
        LZ  = ND+2 - 2*L2L
        FMEAN(L) = FJ(LZ)
      enddo

!---fluxes reflected at top, incident at bottom (1/2 layer flux)
        FJTOP = FJT
        FJBOT = FJB

      return
      end subroutine OPMIE
!EOC
!-------------------------------------------------------------------------

      end module FastJX53cCoreFastj_mod
