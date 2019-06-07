c  >>>>>>>>>>>>>>>>current code revised to JX ver 5.3c (4/06)<<<<<<<<<<<<
c version 5.3c adds
c      calculates reflected UV-vis solar energy (relative to 1.0)
c      new solar spectrum (J-O2 increases in strat by 10%, J-NO by 15+%)
c
c version 5.3b changes include:
c      new data files for specral Xsection and mie-scattering.
c      add sub-layers (JXTRA) to thick cloud/aerosol layers,
c           sets up log-spaced sub-layers of increasing thickness ATAU
c      correction 'b' does massive clean up of the linking code,
c           now the only subroutine that has access to CTM arrays is PHOTOJ
c           Also, the access to the cmn_JVdat.f is 'read-only' after init.
c           This should enable safe openMP/MPI coding.
c
c common files and what they mean:
c   parm_CTM.f  dimensions & params for code (CTM and fast-JX)
c   parm_MIE.f  dimensions for mie code variables.
c   cmn_metdat.f  CTM 3-D arrays, time of day, grid,  etc.
c   cmn_JVdat.f   Xsects, Mie, etc., (initialized and then read-only)
c
c<<<<<<<<<<<<<<<<<<<<<begin CTM-specific subroutines<<<<<<<<<<<<<<<<<<<<
c subroutines:
c
c     SET_ATM(GMTAU)
c           set ups atmosphere (p,T,O3,airmass, etc) for time GMTAU
c              COMMON BLOCKS: cmn_metdat.f
c
c     SET_AER(GMTAU)
c              set ups aerosols for time GMTAU = DUMMY
c                          
c     SET_CLD(GMTAU)
c           set ups clouds for time GMTAU = DUMMY
c
c     INPHOT:  Init. photolysis rate, called once by CHMSET
c              COMMON BLOCKS: cmn_metdat.f, cmn_JVdat.f
c              Input files: ECT42_grid.dat
c
c     RD_JS(NJ1,NAMFIL):  Read labels of photo. rates, called once by INPHOT.
c              COMMON BLOCKS: cmn_metdat.f, cmn_JVdat.f
c              Input files: ratj.dat
c
c     RD_PROF(NJ2,NAMFIL):  Read T & O3 climatology, called once by INPHOT.
c              COMMON BLOCKS: cmn_metdat.f
c              Input files: atmos_std.dat
c
c     SET_CLD0(TINIT,ODINIT,ODINDX,ALBEDO)
c              Initialize cloud and surface properties, called by MAIN.
c              COMMON BLOCKS: cmn_metdat.f
c
c     SET_AER0:  Iniitalize (climatology) aerosol OD and types (3 arrays)
c                     called by MAIN, CHMSET
c              COMMON BLOCKS: cmn_metdat.f
c
c     SET_ATM0:  Initialize climatologies for T & O3, set up atmospheric profiles
c              COMMON BLOCKS: cmn_metdat.f
c
c<<<<<<<<<<<<<<<<<<<<<begin CTM-fastJX linking subroutines<<<<<<<<<<<<<<
c  
c     PHOTOJ(UTIME,IDAY,ILNG,JLAT, SZA,ZPJ)
c              Gateway to fast-JX, Update the photolysis rates
c              COMMON BLOCKS: cmn_metdat.f, cmn_JVdat.f
c
c<<<<<<<<<<<<<<<<<<<<<begin core fast-J subroutines<<<<<<<<<<<<<<<<<<<<<
c  N.B. all these need access to cmn_JVdat.f, but do NOT write into it.
c           also have no need to access cmn_metdat.f
c
c     JRATET(PPJ,TTJ,FFF, VALJL):  Calculate J-value, called by PTOTOJ.
c              COMMON BLOCKS: cmn_JVdat.f
c
c     JP_ATM(PPJ,TTJ,DDJ,ZZJ,ZHL,ZZHT,AER,ABX,ADX,JXTRA)
c              print out atmosphere used in J-value calc.            
c              COMMON BLOCKS: cmn_JVdat.f
c
c
c     RD_XXX(NJ1,NAMFIL):  Read wavelength bins, solar fluxes, Rayleigh
c             parameters, TEM-dependent X-sections, called once by INPHOT.
c              COMMON BLOCKS: cmn_JVdat.f
c              Input files: FJX_spec.dat
c
c     RD_MIE(NJ1,NAMFIL):  Set aerosols/cloud scattering, called once by INPHOT
c              COMMON BLOCKS: cmn_JVdat.f
c              Input files: FJX_scat.dat
c
c     FUNCTION FLINT(TINT,T1,T2,T3,F1,F2,F3)
c
c     SOLARZ(GMTIME,NDAY,YGRDJ,XGRDI, SZA,COSSZA,SOLFX)
c              calc SZA and Solar Flux factor for given lat/lon/UT
c
c     SPHERE(U0,RAD,ZHL,ZZHT,AMF,L1_):  
c              calculate spherical geometry, air-mass factors
c
c     EXTRAL(AER,ADX,L1X,L2X,NX,JTAUMX,ATAU,ATAU0, JXTRA)
c              add sub-layers (JXTRA) to thick cloud/aerosol layers
c
c     OPMIE (KW,KM,WAVEL,ABX,AER,ADX,U0,RFLECT,AMF,JXTRA,FMEAN)
c              calculate mean actinic flux at desired levels
c              COMMON BLOCKS: cmn_JVdat.f
c
c<<<<<<<<<<<<<<<<<<<<<<<begin core scattering subroutines<<<<<<<<<<<<<<<
c
c      MIESCT (FJ,POMEGA,FZ,ZTAU,ZFLUX,ZREFL,ZU0,MFIT,ND)
c            include 'parm_MIE.f' = dimension parameters
c
c      BLKSLV (FJ,POMEGA,FZ,ZTAU,ZFLUX,ZREFL,WT,EMU,PM,PM0,M,N,MFIT,ND)
c              PARAMETER FILE: parm_MIE.f
c
c      GEN (POMEGA,FZ,ZTAU,ZFLUX,ZREFL,WT,EMU,PM,PM0,B,CC,AA,A,H,C1
c            ,M,N,MFIT,ND,ID)
c              PARAMETER FILE: parm_MIE.f
c
c      LEGND0 (X,PL,N)
c
c      MATIN4 (A)
c
c      GAUSSP (N,XPT,XWT)
c
c      EFOLD  (F0, F1, N, F)
c
c-----------------------------------------------------------------------
      program standalone
c-----------------------------------------------------------------------
c  Routine to test Fast-J code by simulating model calls
c
c  >>>>>overall program revised to run with FAST-JX  J-code ver 5.3 (6/05)
c  >>>>>  now designed to use log spacing in TAU to get <0.5% accuracy
c         much better accuracy and fewer added levels vs. old linear adds.
c
c-----------------------------------------------------------------------
c
c<<<<<<<<<<<<<<<<<<<<<begin CTM-specific subroutines<<<<<<<<<<<<<<<<<<<<
c
c  these are mesy and meant to link with or be replace by the CTM code
c  SOME of the important subroutines - eg, SET_CLD need to be updated.
      
      implicit none
      include 'parm_CTM.f'
      include 'cmn_metdat.f'
      include 'cmn_JVdat.f'

      real*8  ZPJ(JVL_,JVN_)    !2-D array of J's indexed to CTM chemistry!
      real*8  FREFL,U0             ! fraction of flux(energy-wtd) reflected
      real*8 GMTAU,DELTAU, TINIT(L_),ODINIT(L_),ALBEDO,SZA
      integer I,J,K,L,NSTEP   ,ILNG,JLAT,ODINDX, IDAY
c-----------------------------------------------------------------------

      write(6,*) ' fast-JX ver-5.3 standalone CTM code'

c  Read in input file
      open(1,file='FJX_test.dat',status='old')
      read(1,*)
      read(1,*)
      read(1,*)
      read(1,*) IDAY
          MONTH=int(dble(IDAY)*12.d0/365.d0)+1    !  Approximately
      read(1,*) NSTEP
      read(1,*) GMTAU
      read(1,*) DELTAU
      read(1,*) JLAT
      read(1,*) ILNG
      read(1,*) P(ILNG,JLAT)
      read(1,*) ALBEDO
      read(1,*) ODINDX
      read(1,*)
c--only for L_=37-layer EC
      do L=1,L_
        read(1,*) J,ETAA(J),ETAB(J),TINIT(J),ODINIT(J)
      enddo
        read(1,*) J,ETAA(J),ETAB(J)    !top layer
      close(1)
c
c  Initial call to Fast-J to set things up

C-----------------------------------------------------------------------
      call INPHOT
C-----------------------------------------------------------------------
      call SET_CLD0(TINIT,ODINIT,ODINDX,ALBEDO)
C-----------------------------------------------------------------------
      call SET_ATM0
          PMEAN(ILNG,JLAT) = P(ILNG,JLAT)
C-----------------------------------------------------------------------
      call SET_AER0
C-----------------------------------------------------------------------


c  Step through SZAs throughout the day:

         GMTAU = GMTAU - DELTAU

      do I=1,NSTEP
         
         GMTAU = GMTAU + DELTAU

         if(GMTAU.ge.24.0) then
           IDAY=mod(IDAY,365)+1
           GMTAU=mod(GMTAU,24.d0)
           MONTH=int(dble(IDAY)*12.d0/365.d0)+1   !do a better calendar 
         endif

c---reset the atmosphere, aerosols, clouds for the time step (now = dummy)
C-----------------------------------------------------------------------
         call SET_ATM(GMTAU)
         call SET_AER(GMTAU)
         call SET_CLD(GMTAU)
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
        call PHOTOJ(GMTAU,IDAY,ILNG,JLAT, SZA,U0,FREFL,ZPJ)
C-----------------------------------------------------------------------

c---printout J's:
         write(6,1001) I,IDAY,GMTAU,SZA, YDGRD(JLAT),XDGRD(ILNG)
         write(6,1002) (JLABEL(K), K=1,JVN_)
        do L=JVL_,1,-1
          write(6,1003) L,(ZPJ(L,K),K=1,JVN_)
        enddo
 1001 format('Step=',i3,'  Day=',i3,' UT(hr)=',f4.1,'  SZA=',f7.2,
     &       '  LAT x LONG=',2f7.2)
 1002 format(1x,'L=  ',64(a7,2x))
 1003 format(i3,1p, 64e9.2)

      enddo
      stop
      end


C-----------------------------------------------------------------------
      subroutine SET_ATM(GMTAU)
C-----------------------------------------------------------------------
      implicit none
      include 'parm_CTM.f'
      include 'cmn_metdat.f'

      real*8 MASFAC,SCALEH,GMTAU, PJ(L1_+1)
      integer I,J,L,N
c
c  Mass factor - delta-Pressure (mbars) to delta-Column (molecules.cm-2)
      MASFAC = 100.d0*6.022d+23/(28.97d0*9.8d0*10.d0)

      do J = 1,J_
      do I = 1,I_
        P(I,J) = PMEAN(I,J)  
        do L = 1,L1_
          PJ(L) = ETAA(L) + ETAB(L)*P(I,J)
        enddo
          PJ(L1_+1) = 0.d0
        do L = 1,L1_
          DM(I,J,L)  = (PJ(L)-PJ(L+1))*MASFAC
        enddo
        do L = 1,L_
          TJ(I,J,L)  = T(I,J,L)
        enddo
c-------calculate effective altitude of each CTM level edge
          ZH(I,J,1) = 16d5*log10(1000.d0/P(I,J))
        do L = 1,L_
          SCALEH      = 1.3806d-19*MASFAC*TJ(I,J,L)
          ZH(I,J,L+1) = ZH(I,J,L) -(log(PJ(L+1)/PJ(L))*SCALEH)
        enddo
      enddo
      enddo

c---load O3 from CTM is being calculated:
      do N = 1,NTR_
        if (TCNAME(N) .eq. 'O3') then
          do J = 1,J_
          do I = 1,I_
          do L = 1,L_         ! unit of DO3:  # molecules/cm^2
             DO3(I,J,L) = 6.023d26*STT(I,J,L,N)/48.d0
     &                    *1.d-4 /AREAXY(I,J)
            enddo
          enddo
          enddo
        endif
      enddo

      return 
      end

      subroutine SET_AER(GMTAU)
      return 
      end

      subroutine SET_CLD(GMTAU)
      return 
      end
      

c-----------------------------------------------------------------------
      subroutine INPHOT
C-----------------------------------------------------------------------
C  Routine to initialise photolysis rate data, called directly from the
C  cinit routine in ASAD. Currently use it to read the JPL spectral data
C  and standard O3 and T profiles and to set the appropriate reaction index.
C-----------------------------------------------------------------------
c
c     IPH       Channel number for reading all data files
c     RAD       Radius of Earth (cm)
c     ZZHT      Effective scale height above top of atmosphere (cm)
c     DATUMX    Maximum opt.depth above which sub layers should be inserted
c     SZAMAX    Solar zenith angle cut-off, above which to skip calculation
c
C-----------------------------------------------------------------------
c
      implicit none
      include 'parm_CTM.f'
      include 'cmn_metdat.f'
      include 'cmn_JVdat.f'

      integer  IPH, I, J, K
      real*8   PI
c
      PI   = 3.141592653589793d0
c Use channel 8 to read files at the moment
      IPH  = 8
c
c Defaults & constants
      RAD  = 6375.d5
      ZZHT = 5.d5
      STT(:,:,:,:) = 1.d6
      TCNAME(:)    = 'CO2'
c
c Read in annual mean pressure field
      open (IPH,file='ECT42_grid.dat',form='formatted',status='old')
      read (IPH,*)
      read (IPH,'(16F8.2)') (XDGRD(I),I=1,I_)
      read (IPH,*)
      read (IPH,'(16F8.2)') (YDGRD(J),J=1,J_)
      read (IPH,*)
      read (IPH,'(16F8.2)') (AREAXY(1,J),J=1,J_)
      read (IPH,*)
      read (IPH,'(16F8.1)') ((PMEAN(I,J),I=1,I_),J=1,J_)
      close(IPH)
      
      do J=1,J_
        AREAXY(1,J) = AREAXY(1,J)*1.d9 
        do I=2,I_
        AREAXY(I,J) = AREAXY(1,J)
        enddo
      enddo

      do I = 1,I_
        XGRD(I) = XDGRD(I) *PI/180.d0
      enddo
      do J = 1,J_
        YGRD(J) = YDGRD(J) *PI/180.d0
      enddo
c
c Read in fast-J X-sections (spectral data) <<<<<<<<<<<<<< new fast-JX
      call RD_XXX(IPH,'FJX_spec.dat')

c Read in aerosol/cloud scattering data <<<<<<<<<<<<<<<<<< new fast-JX
      call RD_MIE(IPH,'FJX_scat.dat')

c Read in labels of photolysis rates required   >>>>> keyed to chem code
      call RD_JS(IPH,'ratj.dat')
c
c Read in T & O3 climatology                    >>>> general backup clim.
      call RD_PROF(IPH,'atmos_std.dat')
c
      return
      end


      subroutine RD_JS(NJ1,NAMFIL)
C-----------------------------------------------------------------------
c  Reread the ratj.dat file to map photolysis rate to reaction
c  Read in quantum yield 'jfacta' and fastj2 label 'jlabel'
C-----------------------------------------------------------------------
c
c     jfacta    Quantum yield (or multiplication factor) for photolysis
c     jlabel    Reference label identifying appropriate J-value to use
c     ipr       Photolysis reaction counter - should total 'JVN_'
c
C-----------------------------------------------------------------------
      implicit none
      include 'parm_CTM.f'
      include 'cmn_metdat.f'
      include 'cmn_JVdat.f'
c
      integer, intent(in) ::  NJ1
      character(*), intent(in) ::  NAMFIL

      integer  IPR, I, J, K
      character*120 CLINE
c
c Reread the ratj.dat file to map photolysis rate to reaction
c                     Read in quantum yield jfacta and fastj2 label jlabel
      IPR = 0
      open (NJ1,file=NAMFIL,status='old',form='formatted')
 10   read (NJ1,'(A)',err=20)  CLINE
      if (IPR .eq. JVN_) goto 20

      if (CLINE(2:5).eq.'9999') then
        go to 20
      elseif (CLINE(1:1).eq.'#') then
        go to 10
      elseif (CLINE(5:5).eq.'$') then
        go to 10
      else
        IPR = IPR+1
        read (CLINE(79:83),'(F5.1)') JFACTA(IPR)
        read (CLINE(86:92),'(A7)')   JLABEL(IPR)
        JFACTA(IPR) = JFACTA(IPR)/100.d0
        go to 10
      endif
 20   close(NJ1)

      NRATJ = IPR

C-----------------------------------------------------------------------
c  compare Xsections titles with J-values listed in chem code (jratd.dat)
c  map the J-values needed for chemistry (ratj.dat) onto the fast-JX rates
c  >>>>>>>>>>>>>>>>current code revised to JPL-02 ver 8.5 (5/05)<<<<<<<<<
c          >>>this must now follow the read in of Xsects, etc<<<
C-----------------------------------------------------------------------

c---Zero / Set index arrays that map Jvalue(j) onto rates
      do J = 1,JVN_
        JIND(J) = 0
      enddo
      do J = 1,NJVAL
      do K = 1,NRATJ
        if (JLABEL(K) .eq. TITLEJ(J)) JIND(K)=J
      enddo
      enddo

      write(6,'(a,i4,a)') ' Photochemistry Scheme with ',IPR,' J-values'
      do K=1,NRATJ
        J = JIND(K)
        if (J.eq.0) then
         write(6,'(i5,a9,f6.2,a,i4,a9)') K,JLABEL(K),JFACTA(K),
     &         ' has no mapping onto onto fast-JX'
        else
         write(6,'(i5,a9,f6.2,a,i4,a9)') K,JLABEL(K),JFACTA(K),
     &         ' mapped onto fast-JX:',J,TITLEJ(J)
        endif
      enddo  

      return
      end


C-----------------------------------------------------------------------
      subroutine RD_PROF(NJ2,NAMFIL)
C-----------------------------------------------------------------------
c  Routine to input T and O3 reference profiles
C-----------------------------------------------------------------------
      implicit none
      include 'parm_CTM.f'
      include 'cmn_metdat.f'

      integer, intent(in) ::  NJ2
      character(*), intent(in) ::  NAMFIL
c
      integer IA, I, M, L, LAT, MON, NTLATS, NTMONS, N216
      real*8  OFAC, OFAK
      character*72 TITLE0
c
      open (NJ2,file=NAMFIL,status='old',form='formatted')
      read (NJ2,'(A)') TITLE0
      read (NJ2,'(2I5)') NTLATS,NTMONS
c      write(6,'(1X,A)') TITLE0
      write(6,1000) NTLATS,NTMONS
      N216  = min(216, NTLATS*NTMONS)
      do IA = 1,N216
        read (NJ2,'(1X,I3,3X,I2)') LAT, MON
        M = min(12, max(1, MON))
        L = min(18, max(1, (LAT+95)/10))
        read (NJ2,'(3X,11F7.1)') (TREF(I,L,M), I=1,41)
        read (NJ2,'(3X,11F7.4)') (OREF(I,L,M), I=1,31)
      enddo
      close (NJ2)

c  Extend climatology to 100 km
      OFAC = exp(-2.d5/5.d5)
      do I = 32,51
        OFAK = OFAC**(I-31)
        do M = 1,NTMONS
        do L = 1,NTLATS
          OREF(I,L,M) = OREF(31,L,M)*OFAK
        enddo
        enddo
      enddo
      do L = 1,NTLATS
      do M = 1,NTMONS
      do I = 42,51
        TREF(I,L,M) = TREF(41,L,M)
      enddo
      enddo
      enddo

      return
 1000 format(1x,'std atmos profiles: ',i3,' lat x ',i2,' mon')
      end

      

c-----------------------------------------------------------------------
      subroutine SET_CLD0(TINIT,ODINIT,ODINDX,ALBEDO)
c-----------------------------------------------------------------------
c  Routine to set cloud and surface properties: now loads the input prof's
c           
c  >>>>>> this subroutine will need to be customized
c  >>>>>> it is separate from aerosols since it comes from met fields
c-----------------------------------------------------------------------
      implicit none
      include 'parm_CTM.f'
      include 'cmn_metdat.f'
      include 'cmn_JVdat.f'

      integer, intent(in) ::  ODINDX
      real*8, intent(in) :: TINIT(L_),ODINIT(L_),ALBEDO
      integer  I, J, L
c
      do J = 1,J_
      do I = 1,I_
        do L = 1,L_
          ODCLD(I,J,L) = ODINIT(L)
          T(I,J,L)     = TINIT(L)
        enddo
          ODCLD(I,J,L1_)  = 0.d0
          SA(I,J) = ALBEDO
      enddo
      enddo
          NCLDX(:,:,:) = ODINDX
c
      return
      end


C-----------------------------------------------------------------------
      subroutine SET_AER0
C-----------------------------------------------------------------------
c  Set up aerosols >>>>customize for climatology or CTM as source.
C-----------------------------------------------------------------------
      implicit none
      include 'parm_CTM.f'
      include 'cmn_metdat.f'
      include 'cmn_JVdat.f'

      integer I, J, L, K

c  Add aerosol OD and indices
c    DAERn(L=1:L_) = aerosol/cloud Optical Depth (OD) within CTM layer L
c    DAERn(L1_=L_+1)= aerosol/cloud OD in layer above CTM (L1_)

          NAER1(:,:,:) = 0
          NAER2(:,:,:) = 0
          NAER3(:,:,:) = 0
      do J = 1,J_
      do I = 1,I_
        do L = 1,L_+1
          DAER1(I,J,L) = 0.d0
          DAER2(I,J,L) = 0.d0
          DAER3(I,J,L) = 0.d0
        enddo
      enddo
      enddo

      end


c-----------------------------------------------------------------------
      subroutine SET_ATM0
C-----------------------------------------------------------------------
c  Routine to set up atmospheric profiles required by Fast-J2 using a
c  doubled version of the level scheme used in the CTM. First pressure
c  and z* altitude are defined, then O3 and T are taken from the supplied
c  climatology and integrated to the CTM levels (may be overwritten with
c  values directly from the CTM, if desired).
c                                       Oliver (04/07/99) & MJP (7/05)
C-----------------------------------------------------------------------
c
c
c     PJ       Pressure at boundaries of model levels (hPa)
c     MASFAC   Conversion factor for pressure to column density
c
c     TJ       Temperature profile on model grid
c     DM       Air column for each model level (molecules.cm-2)
c     DO3      Ozone column for each model level (molecules.cm-2)
c     ZH       Altitude of boundaries of model levels (cm)
c     PSTD     Approximate pressures of levels for supplied climatology
C-----------------------------------------------------------------------
      implicit none
      include 'parm_CTM.f'
      include 'cmn_metdat.f'

      integer  I, J, K, L, M, N
      real*8   PSTD(52),OREF2(51),TREF2(51),PJ(L1_+1)
      real*8   DLOGP,F0,T0,PB,PC,XC,MASFAC,SCALEH
c
c  Select appropriate month
      M = max(1,min(12,MONTH))

c  Mass factor - delta-Pressure (mbars) to delta-Column (molecules.cm-2)
      MASFAC = 100.d0*6.022d+23/(28.97d0*9.8d0*10.d0)

C-----------------------------------------------------------------------
      do J = 1,J_

c  Select appropriate latitudinal profiles
          N = max(1, min(18, (int(YDGRD(J))+99)/10 ))
c  Temporary zonal arrays for climatology data
          do K = 1,51
            OREF2(K) = OREF(K,N,M)
            TREF2(K) = TREF(K,N,M)
          enddo

        do I = 1,I_

c  Apportion O3 and T on supplied climatology z* levels onto CTM levels
c  with mass (pressure) weighting, assuming constant mixing ratio and
c  temperature half a layer on either side of the point supplied.
c  L1_ = L_+1:  
c       PJ(L=1:L1_) = pressure at CTM layer edge, PJ(L1_+1)=0 (top-of-atmos)

           PJ(L1_+1) = 0.d0
         do K = 1,L1_
           PJ(K) = ETAA(K) + ETAB(K)*PMEAN(I,J)
         enddo

c  Set up pressure levels for O3/T climatology - assume that value
c  given for each 2 km z* level applies from 1 km below to 1 km above,
c  so select pressures at these boundaries. Surface level values at
c  1000 mb are assumed to extend down to the actual P(ILNG,JLAT).
c
           PSTD(1) = max(PJ(1),1000.d0)
           PSTD(2) = 1000.d0 * 10.d0**(-1.d0/16.d0)
           DLOGP   = 10.d0**(-2.d0/16.d0)
         do K = 3,51
           PSTD(K) = PSTD(K-1)*DLOGP
         enddo
           PSTD(52)  = 0.d0

         do L = 1,L1_
           F0 = 0.d0
           T0 = 0.d0
           do K = 1,51
             PC   = min(PJ(L),PSTD(K))
             PB   = max(PJ(L+1),PSTD(K+1))
             if (PC .gt. PB) then
               XC = (PC-PB)/(PJ(L)-PJ(L+1))
               F0 = F0 + OREF2(K)*XC
               T0 = T0 + TREF2(K)*XC
             endif
           enddo
           TJ(I,J,L)  = T0
           DM(I,J,L)  = (PJ(L)-PJ(L+1))*MASFAC
           DO3(I,J,L) = F0*1.d-6*DM(I,J,L)
         enddo

c  Calculate effective altitudes using scale height at each level
           ZH(I,J,1) = 0.d0
         do L = 1,L_
           SCALEH      = 1.3806d-19*MASFAC*TJ(I,J,L)
           ZH(I,J,L+1) = ZH(I,J,L) -( LOG(PJ(L+1)/PJ(L)) * SCALEH )
         enddo

        enddo
      enddo

      return
      end

c<<<<<<<<<<<<<<<<<<<<<<<<end CTM-specific subroutines<<<<<<<<<<<<<<<<<<<






c<<<<<<<<<<<<<<<<<<<<<begin CTM-fastJX linking subroutines<<<<<<<<<<<<<<
      
c-----------------------------------------------------------------------
      subroutine PHOTOJ(UTIME,IDAY,ILNG,JLAT, SZA,U0,FREFL,ZPJ)
c-----------------------------------------------------------------------
c
c  PHOTOJ is the gateway to fast-JX calculations:
c        only access to CTM 3-D GLOBAL arrays
c        sets up the 1-D column arrays for calculating J's

c
C-----------------------------------------------------------------------
c     AVGF   Attenuation of beam at each level for each wavelength
c     FFF    Actinic flux at each desired level
c     XQO2   Absorption cross-section of O2
c     XQO3   Absorption cross-section of O3
C-----------------------------------------------------------------------
      implicit none
      include 'parm_CTM.f'
      include 'parm_MIE.f'
      include 'cmn_metdat.f'
      include 'cmn_JVdat.f'

      real*8, intent(in)  ::  UTIME,SZA
      integer,intent(in)  ::  IDAY,ILNG,JLAT
      real*8, intent(out) ::  ZPJ(JVL_,JVN_)   !2-D array of J's indexed to CTM chemistry!
      real*8, intent(out) ::  FREFL,U0         !fraction of energy reflected
      
c-----------------------------------------------------------------------
c--------key amtospheric data needed to solve plane-parallel J---------
      real*8, dimension(5,L1_) :: AER
      integer,dimension(5,L1_) :: ADX
      real*8, dimension(L1_)   :: ABX, TTJ,DDJ,ZZJ,ZHL
      real*8, dimension(L1_+1) :: PPJ 
c     integer,dimension(L1_)   :: 
      integer,dimension(L2_+1) :: JXTRA
      
      real*8                   :: RFLECT,SOLF,FREFS,FTOP,FBOT
      real*8                   :: AMF(L1_,L1_)
c-----------------------------------------------------------------------
c------------key arrays AFTER solving for J's---------------------------
      real*8  FFF(W_,JVL_),VALJ(X_)
      real*8  VALJL(JVL_,NJVAL) !2-D array of J_s returned by JRATET

      integer  I,J,K,L,M,KM
      real*8   AVGF(L_),XQO3,XQO2    ,WAVE, FLINT, TTT
C-----------------------------------------------------------------------
      ZPJ(:,:)  = 0.d0
      FFF(:,:) = 0.d0
      FREFL = 0.d0
      FREFS = 0.d0
C-----------------------------------------------------------------------
      call SOLARZ(UTIME,IDAY,YGRD(JLAT),XGRD(ILNG), SZA,U0,SOLF)
C-----------------------------------------------------------------------
      SOLF = 1.d0   ! this needs to be dropped to include 6.7% annual cycle

c---check for dark conditions SZA > 98.0 deg => tan ht = 63 km
c                        or         99.                  80 km
      if (SZA .gt. SZAMAX) goto 99

c---load the amtospheric column data
      do L = 1,L1_
        PPJ(L) = ETAA(L) + ETAB(L)*P(ILNG,JLAT)
        TTJ(L) = TJ(ILNG,JLAT,L)
        DDJ(L) = DM(ILNG,JLAT,L)
        ZZJ(L) = DO3(ILNG,JLAT,L)
      enddo
        PPJ(L1_+1) = 0.d0

c---calculate spherical weighting functions (AMF: Air Mass Factor)
      do L = 1,L1_
        ZHL(L) = ZH(ILNG,JLAT,L)
      enddo

C-----------------------------------------------------------------------
      call SPHERE(U0,RAD,ZHL,ZZHT,AMF,L1_)
C-----------------------------------------------------------------------

c---load the profiles of aerosols & clouds: treated as the same from here
      do L = 1,L1_
        AER(1,L)  = DAER1(ILNG,JLAT,L)  ! Opt. Depth aerosol 1 in layer L
        AER(2,L)  = DAER2(ILNG,JLAT,L)  ! Opt. Depth aerosol 2 in layer L
        AER(3,L)  = DAER3(ILNG,JLAT,L)  ! Opt. Depth aerosol 3 in layer L
        AER(4,L)  = ODCLD(ILNG,JLAT,L)  ! cloud Opt. Depth in L
        AER(5,L)  = 0.d0                ! save space for Rayleigh
        ADX(1,L)  = NAER1(ILNG,JLAT,L)  ! index for aerosol 1 at layer L
        ADX(2,L)  = NAER2(ILNG,JLAT,L)  ! index for aerosol 2 at layer L
        ADX(3,L)  = NAER3(ILNG,JLAT,L)  ! index for aerosol 3 at layer L
        ADX(4,L)  = NCLDX(ILNG,JLAT,L)  ! index for cloud at layer L
        ADX(5,L)  = 1                   ! index for Rayleigh phase in L
      enddo

      do L = 1,L1_
       do I=1,4
        ADX(I,L) = min(NAA, max(0, ADX(I,L)))
       enddo
      enddo

c---Now given the aerosol+cloud OD/layer in visible (600 nm) can calculate
c        how to add additonal levels at top of clouds (now uses log spacing)
C-----------------------------------------------------------------------
      call EXTRAL(AER,ADX,L1_,L2_,N_,JTAUMX,ATAU,ATAU0, JXTRA)
C-----------------------------------------------------------------------

c---set surface reflectance
        RFLECT = SA(ILNG,JLAT)
        RFLECT = max(0.d0,min(1.d0,RFLECT))

C---Loop over all wavelength bins to calc mean actinic flux AVGF(L)
      do K = NW1,NW2
        WAVE = WL(K)
C---Pick nearest Mie wavelength, no interpolation--------------
                               KM=1  ! use 300 nm aerosol properties for <355 nm
        if( WAVE .gt. 355.d0 ) KM=2  ! use 400 nm prop for 355-500 nm
        if( WAVE .gt. 500.d0 ) KM=3
        if( WAVE .gt. 800.d0 ) KM=4 

c---Loop over CTM layers L=1:L1_ = 1:L_+1,
c     values at L1_=L_+1 are a pseudo layer above the top CTM layer (L_)
        do L = 1,L1_
         TTT     = TTJ(L)
         XQO3 = FLINT(TTT,TQQ(1,2),TQQ(2,2),TQQ(3,2)
     &                      ,QO3(K,1),QO3(K,2),QO3(K,3))
         XQO2 = FLINT(TTT,TQQ(1,1),TQQ(2,1),TQQ(3,1)
     &                      ,QO2(K,1),QO2(K,2),QO2(K,3))

         ABX(L) = XQO3*ZZJ(L) + XQO2*DDJ(L)*0.20948d0

         AER(5,L) = DDJ(L)*QRAYL(K)
        enddo

C-----------------------------------------------------------------------
      call OPMIE
     & (K,KM, WAVE, ABX,AER,ADX,U0,RFLECT,AMF,JXTRA,AVGF,FTOP,FBOT)
C-----------------------------------------------------------------------

        do L = 1,JVL_
          FFF(K,L) = FFF(K,L) + SOLF*FL(K)*AVGF(L)
        enddo
          FREFL = FREFL + SOLF*FL(K)*FTOP/WAVE
          FREFS = FREFS + SOLF*FL(K)/WAVE

      enddo       ! loop over wavelength K

          FREFL = FREFL/FREFS      !calculate reflected flux (energy weighted)

C-----------------------------------------------------------------------
      call JRATET(PPJ,TTJ,FFF, VALJL)
C-----------------------------------------------------------------------

c---map the J-values from fast-JX onto ASAD ones (use JIND & JFACTA)
      do L = 1,JVL_
        do J = 1,NRATJ
          if (JIND(J).gt.0) then 
            ZPJ(L,J) = VALJL(L,JIND(J))*JFACTA(J)
          else
            ZPJ(L,J) = 0.d0
          endif
        enddo
      enddo

c---diagnostics that are NOT returned to the CTM code

C-----------------------------------------------------------------------
      write(6,*) 'fast-JX(5.3b)--------------internal print------------'

      call JP_ATM(PPJ,TTJ,DDJ,ZZJ,ZHL,ZZHT,AER,ABX,ADX,JXTRA)

C---Print solar flux terms
      write(6,'(a,5f10.4)') ' fast-JX: SZA/u0/Frefl/albedo/SOL*fact/',
     &  SZA,U0,FREFL,RFLECT,SOLF
      write(6,'(a5,18i9)')   ' bin:',(K, K=NW1,NW2)
      write(6,'(a5,18f9.2)') ' wvl:',(WL(K), K=NW1,NW2)
      do L = JVL_,1,-1
        write(6,'(i3,2x,18f9.6)') L,(FFF(K,L)/FL(K),K=NW1,NW2)
      enddo
      write(6,*) 'fast-JX(5.3b)--------------internal print------------'
C-----------------------------------------------------------------------
      
   99 continue
      return
      end

c<<<<<<<<<<<<<<<<<<<<<<<end CTM-fastJX linking subroutines<<<<<<<<<<<<<<



c<<<<<<<<<<<<<<<<<<<<<begin core fast-J subroutines<<<<<<<<<<<<<<<<<<<<<


c-----------------------------------------------------------------------
      subroutine JRATET(PPJ,TTJ,FFF, VALJL)
c-----------------------------------------------------------------------
c in:
c        PPJ(L1_+1) = pressure profile at edges
c        TTJ(L1_) = = temperatures at mid-level
c        FFF(K=1:NW, L=1:JVL_) = mean actinic flux 
c out:
c        VALJ(JVL_,NJVAL)  JVL_ = no of levels
C-----------------------------------------------------------------------
      implicit none
      include 'parm_CTM.f'
      include 'cmn_JVdat.f'

      real*8, intent(in)  ::  PPJ(L1_+1),TTJ(L1_)
      real*8, intent(in)  ::  FFF(W_,JVL_)
      real*8, intent(out) ::  VALJL(JVL_,NJVAL)

      real*8  FLINT             ! external function for X-sections
      real*8  VALJ(X_)          ! temp for call J's at one L
      real*8  QO2TOT, QO3TOT, QO31DY, QO31D, QQQT, TFACT
      real*8  TT,PP,DD,TT200,TFACA,TFAC0,TFAC1,TFAC2,QQQA,QQ2,QQ1A,QQ1B
      integer J,K,L, IV
c
      do L = 1,JVL_    ! master loop over layer = L

c---need temperature and density (for some quantum yields):
c---in this case the Pressures PPJ are defined at the boundaries,
c---                Temperatures in the middle of each layer
        TT   = TTJ(L)
        PP  = (PPJ(L)+PPJ(L+1))*0.5d0
          if (L.eq.1) PP = PPJ(1)
        DD = 7.24e18*PP/TT

        do J = 1,NJVAL
          VALJ(J) = 0.d0
        enddo

        do K = NW1,NW2                    ! Using model 'T's here
           QO3TOT = FLINT(TT,TQQ(1,2),TQQ(2,2),TQQ(3,2)
     &                       ,QO3(K,1),QO3(K,2),QO3(K,3))
           QO2TOT = FLINT(TT,TQQ(1,1),TQQ(2,1),TQQ(3,1)
     &                       ,QO2(K,1),QO2(K,2),QO2(K,3))
           QO31DY = FLINT(TT,TQQ(1,3),TQQ(2,3),TQQ(3,3)
     &                       ,Q1D(K,1),Q1D(K,2),Q1D(K,3))
           QO31D  = QO31DY*QO3TOT
          VALJ(1) = VALJ(1) + QO2TOT*FFF(K,L)
          VALJ(2) = VALJ(2) + QO3TOT*FFF(K,L)
          VALJ(3) = VALJ(3) + QO31D*FFF(K,L)
        enddo

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

c #52 Methylvinyl ketone   'MeVK  '     q(M) = 1/(1 + 1.67e-19*[M])
          if (TITLEJ(J).eq.'MeVK  ') then
            VALJ(J) = VALJ(J)/(1.0 + 1.67e-19*DD)
          endif
c #55 Methylethyl ketone   MEKeto     q(M) = 1/(1 + 2.0*[M/2.5e19])
          if (TITLEJ(J).eq.'MEKeto') then
            VALJ(J) = VALJ(J)/(1.0 + 0.80E-19*DD)
          endif
c #57 Methyl glyoxal       MGlyxl     q(M) = 1/(1 + 4.15*[M/2.5E19])
          if (TITLEJ(J).eq.'MGlyxl') then
            VALJ(J) = VALJ(J)/(1.0 + 1.66e-19*DD)
          endif

        enddo

      if (TITLEJ(NJVAL-1).eq.'Acet-a') then
c--------------J-ref v8.3 includes Blitz ACETONE q-yields--------------
c---Acetone is a special case:   (as per Blitz et al GRL, 2004)
c---     61 = NJVAL-1 = J1(acetone-a) ==> CH3CO + CH3
c---     62 = NJVAL   = J2(acetone-b) ==> CH3 + CO + CH3
          VALJ(NJVAL-1) = 0.d0
          VALJ(NJVAL)   = 0.d0
c---IV=NJVAL-1 = Xsect (total abs) for Acetone - pre-calc Temp interp factors
        IV    = NJVAL-1
        TFACA = (TT-TQQ(1,IV))/(TQQ(2,IV)-TQQ(1,IV))
        TFACA = max(0.d0, min(1.d0, TFACA))
c---IV=NJVAL = Q2 for Acetone=>(2), specifically designed for quadratic interp.
c---      but force to Q2=0 by 210K
        IV    = NJVAL
        TFAC0 = ( (TT-TQQ(1,IV))/(TQQ(2,IV)-TQQ(1,IV)) )**2
        if (TT .lt. TQQ(1,IV)) then
          TFAC0 = (TT - 210.d0)/(TQQ(1,IV)-210.d0)
        endif
        TFAC0 = max(0.d0, min(1.d0, TFAC0))
c---IV=NJVAL+1 = Q1A for Acetone => (1), allow full range of T = 200K-300K
        IV    = NJVAL+1
        TT200 = min(300.d0, max(200.d0, TT))
        TFAC1 = (TT200-TQQ(1,IV))/(TQQ(2,IV)-TQQ(1,IV))
c---IV=NJVAL+2 = Q1B for Acetone => (1)
        IV    = NJVAL+2
        TFAC2 = (TT200-TQQ(1,IV))/(TQQ(2,IV)-TQQ(1,IV))

c---now integrate over wavelengths
        do K = NW1,NW2
c---NJVAL-1 = Xsect (total abs) for Acetone
          IV   = NJVAL-1
          QQQA = QQQ(K,1,IV) + (QQQ(K,2,IV)-QQQ(K,1,IV))*TFACA
c---NJVAL   = Q2 for Acetone=>(2), specifically designed for quadratic interp.
          IV   = NJVAL
          QQ2  = QQQ(K,1,IV) + (QQQ(K,2,IV)-QQQ(K,1,IV))*TFAC0
          if (TT .lt. TQQ(1,IV)) then
            QQ2 = QQQ(K,1,IV)*TFAC0
          endif
c---NJVAL+1 = Q1A for Acetone => (1), allow full range of T = 200K-300K
          IV   = NJVAL+1
          QQ1A = QQQ(K,1,IV) + (QQQ(K,2,IV)-QQQ(K,1,IV))*TFAC1
c---NJVAL+2 = Q1B for Acetone => (1)   ! scaled to [M]=2.5e19
          IV   = NJVAL+2
          QQ1B = QQQ(K,1,IV) + (QQQ(K,2,IV)-QQQ(K,1,IV))*TFAC2
          QQ1B = QQ1B*4.d-20
c---J(61)
          VALJ(NJVAL-1) = VALJ(NJVAL-1)
     &         + FFF(K,L)*QQQA*(1.d0-QQ2)/(QQ1A + QQ1B*DD)
c---J(62)
          VALJ(NJVAL) = VALJ(NJVAL) + FFF(K,L)*QQQA*QQ2

        enddo    !K
c-----------end v-8.3 includes Blitz ACETONE q-yields--------------
      endif

c----Load array of J-values in native order, need to be indexed/scaled
c    by ASAD-related code later: ZPJ(L,JJ) = VALJL(L,JIND(JJ))*JFACTA(JJ)
        do J=1,NJVAL
          VALJL(L,J) = VALJ(J)
        enddo

      enddo    ! master loop over L=1,JVL_
      return
      end


c-----------------------------------------------------------------------
      subroutine JP_ATM(PPJ,TTJ,DDJ,ZZJ,ZHL,ZZHT,AER,ABX,ADX,JXTRA)
c-----------------------------------------------------------------------
      implicit none
      include 'parm_CTM.f'
c-----------------------------------------------------------------------
c--------key amtospheric data needed to solve plane-parallel J---------
      real*8, dimension(5,L1_) :: AER
      integer,dimension(5,L1_) :: ADX
      real*8, dimension(L1_)   :: ABX, TTJ,DDJ,ZZJ,ZHL
      real*8, dimension(L1_+1) :: PPJ 
      integer,dimension(L2_+1) :: JXTRA
      real*8                   :: ZZHT
c-----------------------------------------------------------------------
      integer  I,J,K,L
      real*8   COL(4),COLO2,COLO3,ZKM,DELZ,ZTOP

        write(6,100)
  100 format('   L z(km)     p      T   ','    d(air)   d(O3)',
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

        write(6,110) L,ZKM,PPJ(L),TTJ(L),DDJ(L)/DELZ,ZZJ(L)/DELZ,
     &       COLO2,COLO3,
     &      (ADX(I,L),COL(I), I=1,4), JXTRA(L+L),JXTRA(L+L-1)
  110 format(1x,i3,0p,f6.2,f10.3,f7.2,1p,4e9.2,4(i3,e9.2),2i3)

        enddo

      return
      end


C-----------------------------------------------------------------------
      subroutine RD_XXX(NJ1,NAMFIL)
C-----------------------------------------------------------------------
c  Read in wavelength bins, solar fluxes, Rayleigh parameters, 
c      T-dependent X-sections. 

c  >>>current code revised to JPL-02 ver 8.5 (5/05)<<<<<

C-----------------------------------------------------------------------
c     NAMFIL   Name of spectral data file (j2_spec.dat) >> j2 for fast-J2
c     NJ1      Channel number for reading data file
c
c     NJVAL    Number of species to calculate J-values for
c     NWWW     Number of wavelength bins, from 1:NWWW
c     WBIN     Boundaries of wavelength bins
c     WL       Centres of wavelength bins - 'effective wavelength'
c     FL       Solar flux incident on top of atmosphere (cm-2.s-1)
c     QRAYL    Rayleigh parameters (effective cross-section) (cm2)
c     QO2      O2 cross-sections
c     QO3      O3 cross-sections
c     Q1D      O3 => O(1D) quantum yield
c     TQQ      Temperature for supplied cross sections
c     QQQ      Supplied cross sections in each wavelength bin (cm2)
C-----------------------------------------------------------------------
      implicit none
      include 'parm_CTM.f'
      include 'cmn_JVdat.f'

      integer, intent(in) :: NJ1
      character(*), intent(in) ::  NAMFIL

      integer  I, J, K, IW, NQQQ, NWWW

      TQQ(:,:) = 0.d0

C----------spectral data----set for new format data J-ver8.3------------------
c         note that NJVAL = # J-values, but NQQQ (>NJVAL) = # Xsects read in
c         for 2005a data, NJVAL = 62 (including a spare XXXX) and 
c              NQQQ = 64 so that 4 wavelength datasets read in for acetone
c         note NQQQ is not used outside this subroutine!

      open (NJ1,FILE=NAMFIL,status='old',form='formatted')
      read (NJ1,100) TITLE0
      read (NJ1,101) NJVAL,NQQQ, NWWW,NW1,NW2
      if (NJVAL.gt.X_ .or. NQQQ.gt.X_) then
        write(6,201) NJVAL,X_
        stop
      endif
      write(6,'(1X,A)') TITLE0
C----J-values:  1=O2, 2=O3P,3=O3D 4=readin Xsects
      read (NJ1,102) (WL(IW),IW=1,NWWW)
      read (NJ1,102) (FL(IW),IW=1,NWWW)
      read (NJ1,102) (QRAYL(IW),IW=1,NWWW)

C---Read O2 X-sects, O3 X-sects, O3=>O(1D) quant yields (each at 3 temps)
      read (NJ1,103) TITLEJ(1),TQQ(1,1), (QO2(IW,1),IW=1,NWWW)
      read (NJ1,103) TITLEJ2,  TQQ(2,1), (QO2(IW,2),IW=1,NWWW)
      read (NJ1,103) TITLEJ3,  TQQ(3,1), (QO2(IW,3),IW=1,NWWW)

      read (NJ1,103) TITLEJ(2),TQQ(1,2), (QO3(IW,1),IW=1,NWWW)
      read (NJ1,103) TITLEJ2,  TQQ(2,2), (QO3(IW,2),IW=1,NWWW)
      read (NJ1,103) TITLEJ3,  TQQ(3,2), (QO3(IW,3),IW=1,NWWW)

      read (NJ1,103) TITLEJ(3),TQQ(1,3), (Q1D(IW,1),IW=1,NWWW)
      read (NJ1,103) TITLEJ2,  TQQ(2,3), (Q1D(IW,2),IW=1,NWWW)
      read (NJ1,103) TITLEJ3,  TQQ(3,3), (Q1D(IW,3),IW=1,NWWW)

      do J = 1,3
        write(6,200) TITLEJ(J),(TQQ(I,J),I=1,3)
      enddo

C---Read remaining species:  X-sections at 2 T_s
      do J = 4,NQQQ
        read (NJ1,103) TITLEJ(J),TQQ(1,J),(QQQ(IW,1,J),IW=1,NWWW)
        read (NJ1,103) TITLEJ2,  TQQ(2,J),(QQQ(IW,2,J),IW=1,NWWW)
          write(6,200) TITLEJ(J),(TQQ(I,J),I=1,2)
      enddo

c  Reset the titles for NJVAL-1 & NJVAL to be the two acetone J_s
c   61: C3H6O  = Acet-a     (CH3CO + CH3) 
c   62: Q2-Ac  = Acet-b     (CH3 + CO + CH3)

      TITLEJ(NJVAL-1) = 'Acet-a'
      TITLEJ(NJVAL)   = 'Acet-b'
      
      close(NJ1)
      
  100 format(a)
  101 format(10x,5i5)
  102 format(10x,    6e10.3/(10x,6e10.3)/(10x,6e10.3))
  103 format(a7,f3.0,6e10.3/(10x,6e10.3)/(10x,6e10.3))
  200 format(1x,' x-sect:',a10,3(3x,f6.2))
  201 format(' Number of x-sections supplied to Fast-J2: ',i3,/,
     &       ' Maximum number allowed (X_) only set to: ',i3,
     &       ' - increase in cmn_jv.f')

      return
      end
      
      
C-----------------------------------------------------------------------
      subroutine RD_MIE(NJ1,NAMFIL)
C-----------------------------------------------------------------------
C-------aerosols/cloud scattering data set for fast-JX (ver 5.3)
c  >>>>>>>>>>>>>>>>spectral data rev to J-ref ver8.5 (5/05)<<<<<<<<<<<<
C-----------------------------------------------------------------------
c     NAMFIL   Name of scattering data file (e.g., FJX_scat.dat)
c     NJ1      Channel number for reading data file
c     NAA      Number of categories for scattering phase functions
c     QAA      Aerosol scattering phase functions
c     NK       Number of wavelengths at which functions supplied (set as 4)
c     WAA      Wavelengths for the NK supplied phase functions
c     PAA      Phase function: first 8 terms of expansion
c     RAA      Effective radius associated with aerosol type
c     SSA      Single scattering albedo
C-----------------------------------------------------------------------
      implicit none
      include 'parm_CTM.f'
      include 'cmn_JVdat.f'

      integer, intent(in) :: NJ1
      character(*), intent(in) ::  NAMFIL
      integer  I, J, K

      open (NJ1,FILE=NAMFIL,status='old',form='formatted')

      read (NJ1,'(i2,a78)') NAA,TITLE0
      read (NJ1,'(5x,i5,2f10.5)') JTAUMX,ATAU,ATAU0
        write(6,'(a,2f9.5,i5)') ' ATAU/ATAU0/JMX',ATAU,ATAU0,JTAUMX
      read (NJ1,*)
      
      do J = 1,NAA
          read (NJ1,'(3x,a20)') TITLEA(J)
        do K = 1,4     ! Fix number of aerosol wavelengths at 4 
          read (NJ1,'(f5.0,f8.1,f7.3,f8.4,f7.3,7f6.3)') 
     &      WAA(K,J),QAA(K,J),RAA(K,J),SSA(K,J),(PAA(I,K,J),I=1,8)
        enddo
      enddo

      close(NJ1)

        write(6,*) 'Aerosol phase functions & wavelengths'
        write(6,*) TITLE0
      do J=1,NAA
        write(6,'(1x,A8,I2,A,9F8.1)')
     &                   TITLEA(J),J,'  wavel=',(WAA(K,J),K=1,4)
        write(6,'(9x,I2,A,9F8.4)') J,'  Qext =',(QAA(K,J),K=1,4)
      enddo

      return
      end


c-----------------------------------------------------------------------
      real*8 FUNCTION FLINT (TINT,T1,T2,T3,F1,F2,F3)
C-----------------------------------------------------------------------
c  Three-point linear interpolation function
C-----------------------------------------------------------------------
      real*8  TINT,T1,T2,T3,F1,F2,F3
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
      end


c-----------------------------------------------------------------------
      subroutine SOLARZ(GMTIME,NDAY,YGRDJ,XGRDI, SZA,COSSZA,SOLFX)
C-----------------------------------------------------------------------
c     GMTIME = UT for when J-values are wanted 
c           (for implicit solver this is at the end of the time step)
c     NDAY   = integer day of the year (used for solar lat and declin)
c     YGRDJ  = laitude (radians) for grid (I,J)
c     XGDRI  = longitude (radians) for grid (I,J)
c
c     SZA = solar zenith angle in degrees
c     COSSZA = U0 = cos(SZA)
C-----------------------------------------------------------------------
      implicit none
      real*8, intent(in) ::   GMTIME,YGRDJ,XGRDI
      integer, intent(in) ::  NDAY
      real*8, intent(out) ::  SZA,COSSZA,SOLFX

      real*8  PI, PI180, LOCT
      real*8  SINDEC, SOLDEK, COSDEC, SINLAT, SOLLAT, COSLAT, COSZ
c
      PI     = 3.141592653589793d0
      PI180  = PI/180.d0
      SINDEC = 0.3978d0*sin(0.9863d0*(dble(NDAY)-80.d0)*PI180)
      SOLDEK = asin(SINDEC)
      COSDEC = cos(SOLDEK)
      SINLAT = sin(YGRDJ)
      SOLLAT = asin(SINLAT)
      COSLAT = cos(SOLLAT)
c
      LOCT   = (((GMTIME)*15.d0)-180.d0)*PI180 + XGRDI
      COSSZA = COSDEC*COSLAT*cos(LOCT) + SINDEC*SINLAT
      SZA    = acos(COSSZA)/PI180

      SOLFX  = 1.d0-(0.034d0*cos(dble(NDAY-186)*2.d0*PI/365.d0))

      return
      end


c-----------------------------------------------------------------------
      subroutine SPHERE(GMU,RAD,ZHL,ZZHT,AMF,L1_)
C-----------------------------------------------------------------------
c  Calculation of spherical geometry; derive tangent heights, slant path
c  lengths and air mass factor for each layer. Not called when
c  SZA > 98 degrees.  Beyond 90 degrees, include treatment of emergent
c  beam (where tangent height is below altitude J-value desired at).
C-----------------------------------------------------------------------
c in:
c     GMU     = MU0 = cos(solar zenith angle)
c     RAD     radius of Earth mean sea level (cm)
c     ZHL(L)  height (cm) of the bottome edge of CTM level L
c     ZZHT    scale height (cm) used above top of CTM (ZHL(L_+1)
c     L1_     dimension of CTM = levels +1
c out:
c     AMF(I,J) = air mass factor for CTM level I for sunlight reaching J
C-----------------------------------------------------------------------
      implicit none
      integer, intent(in) ::   L1_
      real*8, intent(in)  ::   GMU,RAD,ZHL(L1_),ZZHT
      real*8, intent(out) ::   AMF(L1_,L1_)

c     RZ      Distance from centre of Earth to each point (cm)
c     RQ      Square of radius ratios
c     TANHT   Tangent height for the current SZA
c     XL      Slant path between points

      integer  I, J, K, II
      real*8   AIRMAS, XMU1, XMU2, XL, DIFF, TANHT
     &        ,Ux,H,RZ(L1_),RQ(L1_),ZBYR
c
c  Inlined air mass factor function for top of atmosphere
      AIRMAS(Ux,H) = (1.0d0+H)/SQRT(Ux*Ux+2.0d0*H*(1.0d0-
     &         0.6817d0*EXP(-57.3d0*abs(Ux)/SQRT(1.0d0+5500.d0*H))/
     &                                             (1.0d0+0.625d0*H)))
c
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
c
c  Go up from the surface calculating the slant paths between each level
c  and the level above, and deriving the appropriate Air Mass Factor
      do 16 J = 1,L1_
        do I = 1,L1_
          AMF(I,J) = 0.d0
        enddo
c
c  Air Mass Factors all zero if below the tangent height
        if (RZ(J) .lt. TANHT) goto 16
c  Ascend from layer J calculating AMFs
        XMU1 = abs(GMU)
        do I = J,L1_-1
          XMU2     = dsqrt(1.0d0 - RQ(I)*(1.0d0-XMU1**2))
          XL       = RZ(I+1)*XMU2 - RZ(I)*XMU1
          AMF(I,J) = XL / (RZ(I+1)-RZ(I))
          XMU1     = XMU2
        enddo
c  Use function and scale height to provide AMF above top of model
        AMF(L1_,J)  = AIRMAS(XMU1,ZBYR)
c
c  Twilight case - Emergent Beam
        if (GMU .ge. 0.0d0) goto 16
        XMU1       = abs(GMU)
c  Descend from layer J
        do II = J-1,1,-1
          DIFF        = RZ(II+1)*sqrt(1.0d0-XMU1**2)-RZ(II)
          if (II.eq.1)  DIFF = max(DIFF,0.d0)   ! filter
c  Tangent height below current level - beam passes through twice
          if (DIFF .lt. 0.0d0)  then
            XMU2      = sqrt(1.0d0 - (1.0d0-XMU1**2)/RQ(II))
            XL        = abs(RZ(II+1)*XMU1-RZ(II)*XMU2)
            AMF(II,J) = 2.d0*XL/(RZ(II+1)-RZ(II))
            XMU1      = XMU2
c  Lowest level intersected by emergent beam
          else
            XL        = RZ(II+1)*XMU1*2.0d0
            AMF(II,J) = XL/(RZ(II+1)-RZ(II))
            goto 16
          endif
        enddo

   16 continue
      return
      end
      

C-----------------------------------------------------------------------
      subroutine EXTRAL(AER,ADX,L1X,L2X,NX,JTAUMX,ATAU,ATAU0, JXTRA)
C-----------------------------------------------------------------------
c
c    new version 5.3, add sub-layers (JXTRA) to thick cloud/aerosol layers
c    this version sets up log-spaced sub-layers of increasing thickness ATAU
c
c     AER(5,L=1:L1X) = Optical Depth in layer L (general visible OD)
c        This is not used in the calculation of J's but in calculating
c        the number in levels to insert in each layer L
c        Set for log-spacing of tau levels, increasing top-down.
c     AER(1:4,L) = aerosol+cloud OD (up to 4 types, index type = ADX(1:4,L)
c     AER(5,L) = Rayleigh-not used here
c
c     N.B. the TTAU, etc caluclate here are not really used elsewhere

c---The log-spacing parameters have been tested for convergence and chosen
c---  to be within 0.5% for ranges OD=1-500, rflect=0-100%, mu0=0.1-1.0
c---  use of ATAU = 1.18 and min = 0.01, gives at most +135 pts for OD=100 
C-----------------------------------------------------------------------
c
      implicit none
      integer, intent(in) ::  JTAUMX,L1X,L2X  !index of cloud/aerosol
      integer, intent(in) ::  NX              !Mie scattering array size
      real*8,  intent(in) ::  AER(5,L1X)      !cloud+3aerosol OD in each layer
      real*8,  intent(in) ::  ATAU,ATAU0
      integer, intent(in) ::  ADX(5,L1X)      !index of cloud/aerosol
      integer, intent(out) ::  JXTRA(L2X+1)   !number of sub-layers to be added
c
      integer JTOTL,I,L,L2 
      real*8  DTAUX(L1X),TTAU(L2X+1),DTAUJ, ATAU1,ATAULN,ATAUM,ATAUN1
c
C---Reinitialize arrays
      TTAU(:)  = 0.d0
      JXTRA(:) = 0
c
C---Set up total optical depth over each CTM level (L=1:L1X)
c---   DTAUX(L) = bulk properties (OD) of each CTM layer

      do L = 1,L1X
c  Total optical depth from all elements I=1:4 = clouds  + 3 aerosol types
c    Assume here that AER(1:4, 1:L1X) is 'visible' Optical Depth = 600 nm
c    do NOT count if ADX = 0
          DTAUX(L)   = 0.d0
        do I = 1,4
          if (ADX(I,L) .gt. 0) DTAUX(L) = DTAUX(L) + AER(I,L)
        enddo
      enddo
c
c---combine these edge- and mid-layer points into grid of size:
c---              L2X+1 = 2*L1X+1 = 2*L_+3
c---calculate column optical depths above each level, TTAU(1:L2X+1)
c---      note that TTAU(L2X+1)=0 and TTAU(1)=total OD
c
c---Divide thick layers to achieve better accuracy in the scattering code
c---In the original fast-J, equal sub-layers were chosen, this is wasteful
c---and this new code (ver 5.3) uses log-scale:  
c---        Each succesive layer (down) increase thickness by ATAU > 1
c---        e.g., if ATAU = 2, a layer with OD = 15 could be divided into
c---        4 sub-layers with ODs = 1 - 2 - 4 - 8
c---The key parameters are:
c---        ATAU = factor increase from one layer to the next
c---        ATAUMN = the smallest OD layer desired
c---        JTAUMX = maximum number of divisions (i.e., may not get to ATAUMN)
c---These are hardwired below, can be changed, but have been tested/optimized

      ATAU1  = ATAU - 1.d0
      ATAULN = log(ATAU)
        TTAU(L2X+1)  = 0.0d0
      do L2 = L2X,1,-1
        L         = (L2+1)/2
        DTAUJ     = 0.5d0 * DTAUX(L)
        TTAU(L2)  = TTAU(L2+1) + DTAUJ
c---Now compute the number of log-spaced sub-layers to be added in
c---   the interval TTAU(L2) > TTAU(L2+1)
c---The objective is to have successive TAU-layers increasing by factor ATAU >1
c---the number of sub-layers + 1
        if (TTAU(L2) .lt. ATAU0) then
          JXTRA(L2) = 0
        else
          ATAUM    = max(ATAU0, TTAU(L2+1))
          ATAUN1 = log(TTAU(L2)/ATAUM) / ATAULN
          JXTRA(L2) = min(JTAUMX, max(0, int(ATAUN1 - 0.5d0)))
        endif
      enddo

c---check on overflow of arrays, cut off JXTRA at lower L if too many levels
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
      end


C-----------------------------------------------------------------------
      subroutine OPMIE
     &  (KW,KM,WAVEL,ABX,AER,ADX,U0,RFLECT,AMF,JXTRA,FMEAN,FJTOP,FJBOT)
C-----------------------------------------------------------------------
c in:    
c     KW = wavelength bin # (1:18)
c     KM = wavelength index for Mie properties (1:4 = 300-400-600-999 nm)
c     WAVEL = wavelength of bin (in nm) - not now used
c     ABX(L1_) = vertical profiles of ABSORPTION Optical Depth in each layer
c                    includes O2 and O3 for now (BC under aerosols)
c     AER(1:5,1:L1_) = 5 vertical profiles of Optical Depth in each layer
c     ADX(1:5,1:L1_) = integer index of the scattering properties of each AER
c            1:4 are reserved for aerosols and clouds
c            5 is meant only for Rayleigh scattering!
c     JXTRA(1:L1_) = number 0:J = no. of additional levels to be inserted
c     U0  = cos (SZA)
c     RFLECT = Lambertian albedo of surface
c out:
c     FMEAN(1:L_) = mean actinic flux at standard CTM levels (mid-layer)
C-----------------------------------------------------------------------
c
c     DTAUX    Local optical depth of each CTM level
c     TTAU     Optical depth of air vertically above each point (to top of atm)
c     FTAU     Attenuation of solar beam
c     POMEGA   Scattering phase function
c
c---new Ver 5.3 code adds sub-layers (# = JXTRA(L2)) using ATAU as the 
c   factor increase from sub-layer to sub-layer 

C  fast-J Mie code for J_s, only uses 8-term expansion, 4-Gauss pts
C  Currently allow up to A_ aerosol phase functions (at all altitudes) to
C  be associated with optical depth AER(1:L2_) = aerosol opt.depth @ 1000 nm
C  
C  Pick Mie-wavelength with phase function and Qext: e.g., 
C  01 RAYLE = Rayleigh phase
C  02 ISOTR = isotropic
C  03 ABSRB = fully absorbing 'soot', wavelength indep.
C  04 X_Bkg = backgrnd stratospheric sulfate (n=1.46,log-norm:r=.09um/sigma=.6)
C  05 X_Vol = volcanic stratospheric sulfate (n=1.46,log-norm:r=.08um/sigma=.8)
C. . .
C  11 W_C13 = water cloud (C1/Deirm.) (n=1.335, gamma:  r-mode=13.3um /alpha=6)
C. . .
C  13 Ice-H = hexagonal ice cloud (Mishchenko)
C  14 Ice-I = irregular ice cloud (Mishchenko)
C
C  Choice of the 4 aerosol/cloud indices ADX is made earlier
c  Optical depths for the 4 (aerosol+clouds) = AER
c
C-----------------------------------------------------------------------
      implicit none
      include 'parm_CTM.f'
      include 'parm_MIE.f'
      include 'cmn_JVdat.f'

      real*8, intent(in) ::   WAVEL,AER(5,L1_),ABX(L1_),AMF(L1_,L1_)
      real*8, intent(in) ::   U0,RFLECT
      integer, intent(in) ::  KW,KM,ADX(5,L1_)
      integer, intent(in) ::  JXTRA(L2_+1)
      real*8, intent(out) ::  FMEAN(L_),FJTOP,FJBOT
c
      integer JNDLEV(L_),JADDLV(L2_+1),JADDTO(L2_+1),L2LEV(L2_+1)
     &       ,JTOTL,I,II,J,K,L,IX,JK,   L2,L2L,L22,LZ,LZZ,NDZ
c
      real*8  QXMIE(5),XLAER(5),SSALB(5),DTAUX(L1_),PIAER(5,L1_)
     &       ,POMEGAJ(2*M_,L2_+1),TTAU(L2_+1),FTAU(L1_)
     &       ,FTAU2(L2_+1),DTTAU,DFTAU,dPOMEGA(2*M_)
     &       ,XLO2,XLO3,XLRAY,XLTAU,ZK,TAUDN,TAUUP,ZK2
     &       ,DTAUJ,DTAUAD,TAUBTM,TAUTOP,FBTM,FTOP,POMEGAB(2*M_)
     &       ,ATAUA,ATAUZ
c--- variables used in mie code-----------------------------------------
      real*8   FJ(N_),POMEGA(2*M_,N_),FZ(N_),ZTAU(N_),ZREFL,ZU0,ZFLUX
      real*8   FJT,FJB
      integer  MFIT, ND

C---Reinitialize arrays
      ZTAU(:)     = 0.d0
      FZ(:)       = 0.d0
      POMEGA(:,:) = 0.d0


C---Set up optical depth DTAUX(L), and scattering fraction PIAER(1:5,L)
c---    where L = 1:L1_ = bulk properties of each CTM layer.
      do L = 1,L1_

        do I = 1,4
          if (ADX(I,L) .eq. 0)  then
            QXMIE(I) = 0.d0
            SSALB(I) = 0.d0
          else
C---for Mie code scale extinction at 600 nm to wavelength WAVEL (QXMIE)
            QXMIE(I) = QAA(KM,ADX(I,L))/QAA(3,ADX(I,L))
            SSALB(I) = SSA(KM,ADX(I,L))
          endif
        enddo
c---special case for Rayleigh scattering
            QXMIE(5) = 1.d0
            SSALB(5) = 1.d0
        do I = 1,5
          XLAER(I) = AER(I,L)*QXMIE(I)
        enddo

          DTAUX(L) = ABX(L)
        do I = 1,5
          DTAUX(L) = DTAUX(L) + XLAER(I)
        enddo
c---fractional extinction for Rayleigh scattering and each aerosol type
        do I = 1,5
          PIAER(I,L) = SSALB(I)*XLAER(I)/DTAUX(L)
        enddo

      enddo

C---Calculate attenuated incident beam exp(-TTAU/U0 = DTAUX * AirMassFactor)
c---      at the edges of the CTM layers L=1:L1_
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
c
C---Define the total scattering phase fn for each CTM layer L=1:L_+1
c---   from a DTAUX-wt_d mix of aerosols, cloud & Rayleigh
C---No. of quadrature pts fixed at 4(M_), expansion of phase fn @ 8
      MFIT = 2*M_
      do L = 1,L1_
        do I = 1,MFIT
           POMEGAJ(I,L) = 0.d0
          do K = 1,5
           POMEGAJ(I,L)=POMEGAJ(I,L) + PIAER(K,L)*PAA(I,KM,ADX(K,L))
          enddo
        enddo
      enddo

c--printout diagnositcs of different aerosol+cld contrib to OD
c      if(KW.eq.18) then
c      do L=1,L1_
c      write(6,'(3i5,1p,20e12.5)') L,JXTRA(L+L),JXTRA(L+L-1),
c     &   DTAUX(L),(AER(I,L),I=1,5)
c      enddo
c      endif

c
C------------------------------------------------------------------------
c  Take optical properties on CTM layers and convert to a photolysis
c  level grid corresponding to layer centres and boundaries. This is
c  required so that J-values can be calculated for the centre of CTM
c  layers; the index of these layers is kept in the JNDLEV array.
C------------------------------------------------------------------------
c
c---Now combine the CTM layer edges (1:L_+2) with the CTM mid-layer
c---    points (1:L_) plus 1 for the mid point of added top layer.

c---combine these edge- and mid-layer points into grid of size:
c---              L2_+1 = 2*L1_+1 = 2*L_+3
c---calculate column optical depths above each level, TTAU(1:L2_+1)
c---      note that TTAU(L2_+1)=0 and TTAU(1)=total OD
        TTAU(L2_+1) = 0.0d0
      do L2 = L2_,1,-1
        L          = (L2+1)/2
        DTAUJ      = 0.5d0 * DTAUX(L)
        TTAU(L2)   = TTAU(L2+1) + DTAUJ
      enddo
c
C---reflected flux from surface
      if (U0 .gt. 0.d0) then
        ZFLUX = U0*FTAU(1)*RFLECT/(1.d0+RFLECT)
      else
        ZFLUX = 0.d0
      endif
c
c---calculate attenuated beam FTAU2 on the new doubled-levels L2=1:L2_+1
c---       calculate FTAU2 at CTM mid-layers from sqrt
c---       L2_ = 2*L1_ = 2*L_+2

        FTAU2(L2_+1) = 1.0d0
        FTAU2(L2_)   = sqrt(1.0d0*FTAU(L1_))
        FTAU2(L2_-1) = FTAU(L1_)
      do L2 = L2_-3,1,-2
        L           = (L2+1)/2
        FTAU2(L2)   = FTAU(L)
        FTAU2(L2+1) = sqrt(FTAU(L+1)*FTAU(L))
      enddo
c
c  Calculate scattering properties, level centres then level boundaries
c ***be careful of order, we are shifting the 'POMEGAJ' upward in index***
      do L2 = L2_,2,-2
        L   = L2/2
        do I = 1,MFIT
          POMEGAJ(I,L2) = POMEGAJ(I,L)
        enddo
      enddo
c---lower boundary value is set (POMEGAJ(I,1), but set upper:
        do I = 1,MFIT
          POMEGAJ(I,L2_+1) = POMEGAJ(I,L2_)
        enddo
c---now have POMEGAJ filled at even points from L2=3:L2_-1
c---use inverse interpolation for correct tau-weighted values at edges
      do L2 = 3,L2_-1,2
        TAUDN = TTAU(L2-1)-TTAU(L2)
        TAUUP = TTAU(L2)-TTAU(L2+1)
        do I = 1,MFIT
          POMEGAJ(I,L2) = (POMEGAJ(I,L2-1)*TAUDN + 
     &             POMEGAJ(I,L2+1)*TAUUP) / (TAUDN+TAUUP)
        enddo
      enddo
C---at this point FTAU2(1:L2_+1) and POMEAGJ(1:8, 1:L2_+1)
c---    where FTAU2(L2_+1) = 1.0 = top-of-atmos, FTAU2(1) = surface
c
C------------------------------------------------------------------------
c  Calculate cumulative total and define levels we want J-values at.
c  Sum upwards for levels, and then downwards for Mie code readjustments.
c
c     JXTRA(L2)  Number of new levels to add between (L2) and (L2+1)
c           ***JXTRA(1:L2_+1) is calculated based on the aerosol+cloud OD_s
c     JADDLV(L2)  Number of new levels actually added at each wavelength
c            where JADDLV = 0 when there is effectively no FTAU2 
c     JADDTO(L2)   Total number of new levels to add to and above level (L2)
c     JNDLEV(L) = L2 index that maps on CTM mid-layer L
c
C------------------------------------------------------------------------

c---JADDLV(L2=1:L2_) = number of levels to add between TTAU2(L2) and TTAU(L2+1)
c---    JADDLV is taken from JXTRA, which is based on visible OD.
c---    JADDTO(L2=1:L2_+1) is the cumulative number of levels to be added
c---note that JADDLV and JADDTO will change with wavelength and solar zenith

c--now try to insert additional levels for thick clouds, ONLY IF FTAU2 > 1.e-8
c-- this will cut off additional levels where the solar beam is negligible.

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

c---expanded grid now included CTM edge and mid layers plus expanded 
c---    grid to allow for finer delta-tau at tops of clouds.
c---    DIM of new grid = L2_ + JADDTO(1) + 1

c---L2LEV(L2) = L2-index for old level L2 in expanded J-grid (w/JADDLV)
c     in absence of JADDLV, L2LEV(L2) = L2
        L2LEV(1)  = 1
      do L2 = 2,L2_+1
        L2LEV(L2) = L2LEV(L2-1) + 1 + JADDLV(L2-1)
      enddo

c---JNDLEV(L=1:L_) = L2-index in expanded grid for CTM mid-layer L 
      do L = 1,L_
        JNDLEV(L) = L2LEV(2*L)
      enddo

C---------------------SET UP FOR MIE CODE-------------------------------
c
c  Transpose the ascending TTAU grid to a descending ZTAU grid.
c  Double the resolution - TTAU points become the odd points on the
c  ZTAU grid, even points needed for asymm phase fn soln, contain 'h'.
c  Odd point added at top of grid for unattenuated beam   (Z='inf')
c
c  The following mapping holds for JADDLV=0
c        Surface:   TTAU(1)    ==> ZTAU(2*L2_+1)
c        Top:       TTAU(L2_)  ==> ZTAU(3)
c        Infinity:     0.0     ==> ZTAU(1)
c        index: 2*(L2_+1-L2)+1 ==> LZ
c
c  Mie scattering code only used from surface to level L2_
C------------------------------------------------------------------------
c
C------------------------------------------------------------------------
c  Insert new levels, working downwards from the top of the atmosphere
c  to the surface (down in 'LZ', up in 'L2'). This allows ztau and pomega
c  to be incremented linearly (in a +ve sense), and the flux fz to be
c  attenuated top-down (avoiding problems where lower level fluxes are
c  zero).
c
c    zk        fractional increment in level
c    dTTAU     change in ttau per increment    (linear, positive)
c    dPOMEGA   change in pomega per increment  (linear)
c    ftaulog   change in ftau per increment    (exponential, normally < 1)
c
C------------------------------------------------------------------------
c
c  Ascend through atmosphere transposing grid and adding extra points
c  remember L2=1 is surface of CTM, but last layer (LZ) in scattering code.
c  there are twice the number of layers in the LZ arrays (2*L2_ + 2*JADDTO + 1)
c    because we need to insert the intermediate layers (even LZ) for the 
c    asymmetric scattering code.


c  Transfer the L2=1:L2_+1 values (TTAU,FTAU2,POMEGAJ) onto the reverse
c    order, expanded, doubled-level scatter grid. 
c    Note that we need to deal with the expansion by JADD levels (L2L).
c      These JADDLV levels are skipped and need to be interpolated later.
c    Note that only odd LZ levels are filled, 

      NDZ = 2*L2_ + 2*JADDTO(1) + 1

c   Note that the successive sub-layers have the ratio in OD of ATAU
c      ATAUA = (ATAU - 1.d0)/ATAU     ! this is the limit for L22=>inf

      do L2 = 1,L2_+1          ! L2 = index of CTM edge- and mid-layers
        L2L = L2LEV(L2)        ! L2L = index for L2 in expanded scale(JADD)
        LZ  = NDZ + 2 - 2*L2L  ! LZ = index for L2 in scatt arrays
          ZTAU(LZ) = TTAU(L2)
          FZ(LZ)   = FTAU2(L2)
        do I=1,MFIT
          POMEGA(I,LZ) = POMEGAJ(I,L2)
        enddo
      enddo

c   Now go thru the pairs of L2 levels to see if we need JADD levels
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

c---to fit L22 new layers between TAUBOT > TAUTOP, calculate new 1/ATAU factor
c---  such that TAU(just above TAU-btm) = ATUAZ * TAUBTM < TAUBTM

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
            POMEGA(I,LZZ) = POMEGAB(I) + 
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

c   Now fill in the even points with simple interpolation in scatter arrays:
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

c---PRINT diagnostics
c      write(6,'(A,2I6)') 'Total levels of photolysis in OPMIE',ND
      if(ND .gt. N_) then
        write(6,'(a,2i9)') ' overflow of scatter arrays:',ND,N_
        stop
      endif

C-----------------------------------------------------------------------
      call MIESCT(FJ,FJT,FJB,POMEGA,FZ,ZTAU,ZFLUX,ZREFL,ZU0,MFIT,ND)
C-----------------------------------------------------------------------

c---Move mean intensity from scatter array FJ(LZ=1:ND) 
c---              to CTM mid-level array FMEAN(L=1:L_)
      do L = 1,L_
        L2L = JNDLEV(L)
        LZ  = ND+2 - 2*L2L
        FMEAN(L) = FJ(LZ)
      enddo
c---fluxes reflected at top, incident at bottom (1/2 layer flux)
        FJTOP = FJT
        FJBOT = FJB

      return
      end

c<<<<<<<<<<<<<<<<<<<<<<<end core fast-J subroutines<<<<<<<<<<<<<<<<<<<<<


c<<<<<<<<<<<<<<<<<<<<<<<begin core scattering subroutines<<<<<<<<<<<<<<<
C-----------------------------------------------------------------------
      subroutine MIESCT(FJ,FJT,FJB,
     &                    POMEGA,FZ,ZTAU,ZFLUX,ZREFL,ZU0,MFIT,ND)
C-----------------------------------------------------------------------
C   This is an adaption of the Prather radiative transfer code, (mjp, 10/95)
C     Prather, 1974, Astrophys. J. 192, 787-792.
C         Sol_n of inhomogeneous Rayleigh scattering atmosphere. 
C         (original Rayleigh w/ polarization)
C     Cochran and Trafton, 1978, Ap.J., 219, 756-762.
C         Raman scattering in the atmospheres of the major planets.
C         (first use of anisotropic code)
C     Jacob, Gottlieb and Prather, 1989, J.Geophys.Res., 94, 12975-13002.
C         Chemistry of a polluted cloudy boundary layer,
C         (documentation of extension to anisotropic scattering)
C
C    takes atmospheric structure and source terms from std J-code
C    ALSO limited to 4 Gauss points, only calculates mean field!
C
C   mean rad. field ONLY (M=1)
C-----------------------------------------------------------------------
      implicit none
      include 'parm_MIE.f'
c--- expect parameters M_, N_ in parm_MIE.f------------------------------
c
      integer, intent(in) ::  MFIT, ND
      real*8, intent(in)  ::  POMEGA(2*M_,N_),FZ(N_),ZTAU(N_)
     &                       ,ZREFL,ZU0,ZFLUX
      real*8, intent(out) ::  FJ(N_),FJT,FJB
c
      real*8  WT(M_),EMU(M_),PM(M_,2*M_),PM0(2*M_),CMEQ1
      integer I, ID, IM, M, N
C-----------------------------------------------------------------------
C---fix scattering to 4 Gauss pts = 8-stream
      call GAUSSP (N,EMU,WT)
c---calc in OPMIE:  ZFLUX = (ZU0*FZ(ND)*ZREFL)/(1.0d0+ZREFL)
      M = 1
      do I = 1,N
        call LEGND0 (EMU(I),PM0,MFIT)
        do IM = M,MFIT
          PM(I,IM) = PM0(IM)
        enddo
      enddo
C
      CMEQ1 = 0.25d0
      call LEGND0 (-ZU0,PM0,MFIT)
      do IM=M,MFIT
        PM0(IM) = CMEQ1*PM0(IM)
      enddo
C
      call BLKSLV(FJ,POMEGA,FZ,ZTAU,ZFLUX,ZREFL,WT,EMU,PM,PM0
     &            ,FJT,FJB,  M,N,MFIT,ND)
C
      do ID=1,ND,2
        FJ(ID) = 4.0d0*FJ(ID) + FZ(ID)
      enddo

      return
      end


      subroutine BLKSLV(FJ,POMEGA,FZ,ZTAU,ZFLUX,ZREFL,WT,EMU,PM,PM0
     +                 ,FJTOP,FJBOT,  M,N,MFIT,ND)
C-----------------------------------------------------------------------
C  Solves the block tri-diagonal system:
C               A(I)*X(I-1) + B(I)*X(I) + C(I)*X(I+1) = H(I)
C-----------------------------------------------------------------------
      implicit none
      include 'parm_MIE.f'
c--- expect parameters M_, N_ in parm_MIE.f------------------------------
c
      integer, intent(in) ::  M, N, MFIT, ND
      real*8, intent(in)  ::  POMEGA(2*M_,N_),FZ(N_),ZTAU(N_)
     &                       ,WT(M_),EMU(M_),PM(M_,2*M_),PM0(2*M_)
     &                       ,ZFLUX,ZREFL
      real*8, intent(out) ::  FJ(N_),FJTOP,FJBOT
c
      real*8, dimension(M_)    :: A, C1, H
      real*8, dimension(M_,M_) :: B, AA, CC
      real*8                      DD(M_,M_,N_), RR(M_,N_)
      real*8  SUMM
      integer I, J, K, ID
c
C-----------UPPER BOUNDARY ID=1
      call GEN(POMEGA,FZ,ZTAU,ZFLUX,ZREFL,WT,EMU,PM,PM0
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
C----------CONTINUE THROUGH ALL DEPTH POINTS ID=2 TO ID=ND-1
      do ID = 2,ND-1
        call GEN(POMEGA,FZ,ZTAU,ZFLUX,ZREFL,WT,EMU,PM,PM0
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
C---------FINAL DEPTH POINT: ND
      call GEN(POMEGA,FZ,ZTAU,ZFLUX,ZREFL,WT,EMU,PM,PM0
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
C-----------BACK SOLUTION
      do ID = ND-1,1,-1
       do I = 1,N
        do J = 1,N
         RR(I,ID) = RR(I,ID) + DD(I,J,ID)*RR(J,ID+1)
        enddo
       enddo
      enddo
C----------MEAN J & H
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
c       FJBOT = FJBOT + RR(I,1)*WT(I)*EMU(I) !this is not correct, includes up+down
       enddo
        FJTOP = 4.d0*FJTOP

      return
      end


      subroutine GEN(POMEGA,FZ,ZTAU,ZFLUX,ZREFL,WT,EMU,PM,PM0
     &              ,B,CC,AA,A,H,C1,M,N,MFIT,ND,ID)
C-----------------------------------------------------------------------
C  Generates coefficient matrices for the block tri-diagonal system:
C               A(I)*X(I-1) + B(I)*X(I) + C(I)*X(I+1) = H(I)
C-----------------------------------------------------------------------
      implicit none
      include 'parm_MIE.f'
c--- expect parameters M_, N_ in parm_MIE.f------------------------------
c
      integer, intent(in) ::  M, N, MFIT, ND, ID
      real*8, intent(in)  ::  POMEGA(2*M_,N_),FZ(N_),ZTAU(N_)
     &                       ,WT(M_),EMU(M_),PM(M_,2*M_),PM0(2*M_)
     &                       ,ZFLUX,ZREFL
      real*8, intent(out) ::  B(M_,M_),AA(M_,M_),CC(M_,M_),A(M_),C1(M_)
     &                       ,H(M_)
c
      integer ID0, ID1, IM, I, J, K, MSTART
      real*8  SUM0, SUM1, SUM2, SUM3
      real*8  DELTAU, D1, D2, SURFAC
c
      real*8  S(M_,M_), W(M_,M_), U1(M_,M_), V1(M_)
C---------------------------------------------
      if (ID.eq.1 .or. ID.eq.ND) then
C---------calculate generic 2nd-order terms for boundaries
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
C-------------upper boundary, 2nd-order, C-matrix is full (CC)
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
C         H(I) = H(I) + 2.0d0*D2*C1(I) + D1*SISOTP
          H(I) = H(I) + 2.0d0*D2*C1(I)
          A(I) = 0.0d0
        enddo
       else
C-------------lower boundary, 2nd-order, A-matrix is full (AA)
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
C------------intermediate points:  can be even or odd, A & C diagonal
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
      end


      subroutine LEGND0 (X,PL,N)
C---Calculates ORDINARY Legendre fns of X (real) 
C---   from P[0] = PL(1) = 1,  P[1] = X, .... P[N-1] = PL(N)
      implicit none
      integer N,I
      real*8 X,PL(N),DEN
C---Always does PL(2) = P[1]
        PL(1) = 1.d0
        PL(2) = X
        do I = 3,N
         DEN = (I-1)
         PL(I) = PL(I-1)*X*(2.d0-1.0/DEN) - PL(I-2)*(1.d0-1.d0/DEN)
        enddo
      return
      end


      subroutine MATIN4 (A)
C-----------------------------------------------------------------------
C  invert 4x4 matrix A(4,4) in place with L-U decomposition (mjp, old...)
C-----------------------------------------------------------------------
      implicit none
      real*8, intent(INout)  ::  A(4,4)
C---SETUP L AND U
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
C---INVERT L
      A(4,3) = -A(4,3)
      A(4,2) = -A(4,2)-A(4,3)*A(3,2)
      A(4,1) = -A(4,1)-A(4,2)*A(2,1)-A(4,3)*A(3,1)
      A(3,2) = -A(3,2)
      A(3,1) = -A(3,1)-A(3,2)*A(2,1)
      A(2,1) = -A(2,1)
C---INVERT U
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
C---MULTIPLY (U-INVERSE)*(L-INVERSE)
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
      end


      subroutine GAUSSP (N,XPT,XWT)
C-----------------------------------------------------------------------
C  Loads in pre-set Gauss points for 4 angles from 0 to +1 in cos(theta)=mu
C-----------------------------------------------------------------------
      implicit none
      integer, intent(out) :: N
      real*8, intent(out)  ::  XPT(*),XWT(*)
      real*8   GPT4(4),GWT4(4)
      integer  I
      data GPT4/.06943184420297d0,.33000947820757d0,.66999052179243d0,
     &          .93056815579703d0/
      data GWT4/.17392742256873d0,.32607257743127d0,.32607257743127d0,
     &          .17392742256873d0/
      N = 4
      do I = 1,N
        XPT(I) = GPT4(I)
        XWT(I) = GWT4(I)
      enddo
      return
      end


      subroutine EFOLD (F0, F1, N, F)
C-----------------------------------------------------------------------
C     ***not used in fast-J, part of original scattering code
C---Speciality subroutine for calculating consistent exp(-tau/mu0)
C---  values on the tau grid so that photons are conserved.
C---  ***only works for plane-parallel, NOT psuedo-spherical atmos
C
C---  calculate the e-fold between two boundaries, given the value
C---     at both boundaries F0(x=0) = top, F1(x=1) = bottom.
C---  presume that F(x) proportional to exp[-A*x] for x=0 to x=1
C---          d2F/dx2 = A*A*F  and thus expect F1 = F0 * exp[-A]
C---           alternatively, could define A = ln[F0/F1]
C---  let X = A*x, d2F/dX2 = F
C---  assume equal spacing (not necessary, but makes this easier)
C---      with N-1 intermediate points (and N layers of thickness dX = A/N)
C---
C---  2nd-order finite difference:  (F(i-1) - 2F(i) + F(i+1)) / dX*dX = F(i)
C---      let D = 1 / dX*dX:
C
C  1  |   1        0        0        0        0        0   |    | F0 |
C     |                                                    |    | 0  |
C  2  |  -D      2D+1      -D        0        0        0   |    | 0  |
C     |                                                    |    | 0  |
C  3  |   0       -D      2D+1      -D        0        0   |    | 0  |
C     |                                                    |    | 0  |
C     |   0        0       -D      2D+1      -D        0   |    | 0  |
C     |                                                    |    | 0  |
C  N  |   0        0        0       -D      2D+1      -D   |    | 0  |
C     |                                                    |    | 0  |
C N+1 |   0        0        0        0        0        1   |    | F1 |
C      
C-----------------------------------------------------------------------
C  Advantage of scheme over simple attenuation factor: conserves total
C  number of photons - very useful when using scheme for heating rates.
C  Disadvantage: although reproduces e-folds very well for small flux
C  differences, starts to drift off when many orders of magnitude are
C  involved.
C-----------------------------------------------------------------------
      implicit none
      integer, intent(in) :: N
      real*8, intent(in)  :: F0,F1
      real*8, intent(out) :: F(101)     !F(N+1)  
      integer I
      real*8 A,DX,D,DSQ,DDP1, B(101),R(101)
C
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
      end

c<<<<<<<<<<<<<<<<<<<<<<<<<end core scattering subroutines<<<<<<<<<<<<<<<
