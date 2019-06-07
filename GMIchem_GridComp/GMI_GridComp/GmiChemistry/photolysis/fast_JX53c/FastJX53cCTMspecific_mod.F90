!-------------------------------------------------------------------------
!  NASA/GFSC, SIVO, Code 610.3
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: FastJX53cCTMspecific_mod
!
! !INTERFACE:
!
       module FastJX53cCTMspecific_mod
!
       implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
!
       private
       public  :: SET_ATM
       public  :: SET_AER
       public  :: SET_CLD
       public  :: INPHOT
       public  :: SET_CLD0
       public  :: SET_AER0
       public  :: SET_ATM0
!
! !DESCRIPTION:
!  Contains CTM-specific routines.
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
! !IROUTINE: SET_ATM
!
! !INTERFACE:
!
      subroutine SET_ATM(GMTAU)
!
! !USES:
!
      use FastJX53cMetVars_mod      , only : T, P, PJ, PMEAN, TJ, DM, ETAA, ETAB
      use FastJX53cMetVars_mod      , only : ZH, STT, DO3, AREAXY, TCNAME
      use FastJX53cCTMparameters_mod, only : I_, J_, L_, L1_, NTR_
!
      implicit none
!
! !INPUT PARAMETERS:
      real*8,  intent(in) :: GMTAU  ! time
!
! !DESCRIPTION:
!  Sets up atmosphere (p,T,O3,airmass, etc) for time GMTAU.
!
! !LOCAL VARIABLES:
      real*8 MASFAC,SCALEH
!      real*8 MASFAC,SCALEH,PJ(L1_+1)
      integer I,J,L,N
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
!  Mass factor - delta-Pressure (mbars) to delta-Column (molecules.cm-2)
      MASFAC = 100.d0*6.022d+23/(28.97d0*9.8d0*10.d0)

      do J = 1,J_
      do I = 1,I_
        P(I,J) = PMEAN(I,J)  
!Begin Jules Kouatchou 08/24/2006
!        do L = 1,L1_
!          PJ(L) = ETAA(L) + ETAB(L)*P(I,J)
!        enddo
!          PJ(L1_+1) = 0.d0
!End Jules Kouatchou 08/24/2006
        do L = 1,L1_
          DM(I,J,L)  = (PJ(L)-PJ(L+1))*MASFAC
          DO3(I,J,L) = DO3(I,J,L)*DM(I,J,L)      ! Jules Kouatchou 11/02/06
        enddo
        do L = 1,L_
          TJ(I,J,L)  = T(I,J,L)
        enddo
!-------calculate effective altitude of each CTM level edge
          ZH(I,J,1) = 16d5*log10(1000.d0/P(I,J))
        do L = 1,L_
          SCALEH      = 1.3806d-19*MASFAC*TJ(I,J,L)
          ZH(I,J,L+1) = ZH(I,J,L) -(log(PJ(L+1)/PJ(L))*SCALEH)
        enddo
      enddo
      enddo

!Begin Jules Kouatchou 08/24/2006
!!---load O3 from CTM is being calculated:
!      do N = 1,NTR_
!        if (TCNAME(N) .eq. 'O3') then
!          do J = 1,J_
!          do I = 1,I_
!          do L = 1,L_         ! unit of DO3:  # molecules/cm^2
!             DO3(I,J,L) = 6.023d26*STT(I,J,L,N)/48.d0 &
!     &                    *1.d-4 /AREAXY(I,J)
!            enddo
!          enddo
!          enddo
!        endif
!      enddo
!End Jules Kouatchou 08/24/2006

      return 

      end subroutine SET_ATM
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: SET_AER
!
! !INTERFACE:
!
      subroutine SET_AER(GMTAU)
!
      implicit none
!
! !INPUT PARAMETERS:
     real*8 ,intent(in) :: GMTAU   ! dummy time
!
! !DESCRIPTION:
!  Sets up aerosols for time GMTAU.
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
      return 
      end subroutine SET_AER
!EOC      
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: SET_CLD
!
! !INTERFACE:
!
      subroutine SET_CLD(GMTAU)
!
      implicit none
!
! !INPUT PARAMETERS:
     real*8 ,intent(in) :: GMTAU   ! dummy time
!
! !DESCRIPTION:
!  Sets up cloud for time GMTAU.
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
      return
      end subroutine SET_CLD
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: INPHOT
!
! !INTERFACE:
!
      subroutine INPHOT(CrossSection_InfileName,     &
     &                  ScatteringData_InfileName,   &
     &                  PhotRateLabels_InfileName,   &
     &                  T_O3_Climatology_InfileName)
!
! !USES:
!
      use FastJX53cCoreFastj_mod, only : RD_XXX, RD_MIE
      use FastJX53cJvaluesVars_mod, only : RAD, ZZHT
      use FastJX53cMetVars_mod    , only : TCNAME, STT
!
      implicit none
!
! !INPUT PARAMETERS:
      character (len=128), intent (in) :: CrossSection_InfileName
                                        ! fast-J X-sections (spectral data) input file name
      character (len=128), intent (in) :: ScatteringData_InfileName
                                        ! Aerosol/cloud scattering data input file name
      character (len=128), intent (in) :: PhotRateLabels_InfileName
                                        ! Labels of photolysis rates required input file name
                                        ! keyed to chem code
      character (len=128), intent (in) :: T_O3_Climatology_InfileName
                                        ! Read in T & O3 climatology input file name
                                        ! general backup clim.
!
! !DESCRIPTION:
!  Routine to initialise photolysis rate data, called directly from the
!  cinit routine in ASAD. Currently use it to read the JPL spectral data
!  and standard O3 and T profiles and to set the appropriate reaction index.
!  \begin{verbatim}
!     IPH       Channel number for reading all data files
!     RAD       Radius of Earth (cm)
!     ZZHT      Effective scale height above top of atmosphere (cm)
!     DATUMX    Maximum opt.depth above which sub layers should be inserted
!     SZAMAX    Solar zenith angle cut-off, above which to skip calculation
!  \end{verbatim}
!
!
! !LOCAL VARIABLES:
      integer  IPH
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
! Defaults & constants
      RAD  = 6375.d5
      ZZHT = 5.d5
      STT(:,:,:,:) = 1.d6
      TCNAME(:)    = 'CO2'

! Use channel 8 to read files at the moment
      IPH  = 8

! Read in fast-J X-sections (spectral data) <<<<<<<<<<<<<< new fast-JX
      call RD_XXX(IPH, CrossSection_InfileName)

! Read in aerosol/cloud scattering data <<<<<<<<<<<<<<<<<< new fast-JX
      call RD_MIE(IPH, ScatteringData_InfileName)

! Read in labels of photolysis rates required   >>>>> keyed to chem code
      call RD_JS(IPH, PhotRateLabels_InfileName)
!
! Read in T & O3 climatology                    >>>> general backup clim.
      call RD_PROF(IPH, T_O3_Climatology_InfileName)
!
      return
      end subroutine INPHOT
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: RD_JS
!
! !INTERFACE:
!
      subroutine RD_JS(NJ1,NAMFIL)
!
! !USES:
!
      use FastJX53cJvaluesVars_mod  , only : JFACTA, JLABEL, TITLEJ
      use FastJX53cJvaluesVars_mod  , only : JIND, NRATJ, NJVAL
      use FastJX53cCTMparameters_mod, only : JVN_
!
      implicit none
!
!
! !INPUT PARAMETERS:
      integer     , intent(in) ::  NJ1
      character(*), intent(in) ::  NAMFIL
!
! !DESCRIPTION:
!  Reread the ratj.dat file to map photolysis rate to reaction
!  Read in quantum yield 'jfacta' and fastj2 label 'jlabel'
!  \begin{verbatim}
!     jfacta    Quantum yield (or multiplication factor) for photolysis
!     jlabel    Reference label identifying appropriate J-value to use
!     ipr       Photolysis reaction counter - should total 'JVN_'
!  \end{verbatim}
!
! !LOCAL VARIABLES:
      integer        :: IPR, I, J, K
      character*120  :: CLINE
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
! Reread the ratj.dat file to map photolysis rate to reaction
!                     Read in quantum yield jfacta and fastj2 label jlabel
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
!Beg Jules Kouatchou 08/29/06
!        if (JLABEL(IPR) == "C3H6O") JLABEL(IPR) = "Acet-a"
!End Jules Kouatchou 08/29/06
        go to 10
      endif
 20   close(NJ1)

      NRATJ = IPR

!-----------------------------------------------------------------------
!  compare Xsections titles with J-values listed in chem code (jratd.dat)
!  map the J-values needed for chemistry (ratj.dat) onto the fast-JX rates
!  >>>>>>>>>>>>>>>>current code revised to JPL-02 ver 8.5 (5/05)<<<<<<<<<
!          >>>this must now follow the read in of Xsects, etc<<<
!-----------------------------------------------------------------------

!---Zero / Set index arrays that map Jvalue(j) onto rates
      do J = 1,JVN_
        JIND(J) = 0
      enddo
      do J = 1,NJVAL
      do K = 1,NRATJ
        if (JLABEL(K) .eq. TITLEJ(J)) JIND(K)=J
      enddo
      enddo

!      write(6,'(a,i4,a)') ' Photochemistry Scheme with ',IPR,' J-values'
      do K=1,NRATJ
        J = JIND(K)
        if (J.eq.0) then
         write(6,'(i5,a9,f6.2,a,i4,a9)') K,JLABEL(K),JFACTA(K), &
     &         ' has no mapping onto onto fast-JX'
        else
!         write(6,'(i5,a9,f6.2,a,i4,a9)') K,JLABEL(K),JFACTA(K), &
!     &         ' mapped onto fast-JX:',J,TITLEJ(J)
        endif
      enddo  

      return

      end subroutine RD_JS
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:  RD_PROF
!
! !INTERFACE:
!
      subroutine RD_PROF(NJ2,NAMFIL)
!
! !USES:
!
      use FastJX53cMetVars_mod, only : OREF, TREF
!
      implicit none
!
! !INPUT PARAMETERS:
      integer     , intent(in) ::  NJ2
      character(*), intent(in) ::  NAMFIL
!
! !DESCRIPTION:
!  Routine to input T and O3 reference profiles
!
! !LOCAL VARIABLES:
      integer IA, I, M, L, LAT, MON, NTLATS, NTMONS, N216
      real*8  OFAC, OFAK
      character*72 TITLE0
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
      open (NJ2,file=NAMFIL,status='old',form='formatted')
      read (NJ2,'(A)') TITLE0
      read (NJ2,'(2I5)') NTLATS,NTMONS
!      write(6,'(1X,A)') TITLE0
!      write(6,1000) NTLATS,NTMONS
      N216  = min(216, NTLATS*NTMONS)
      do IA = 1,N216
        read (NJ2,'(1X,I3,3X,I2)') LAT, MON
        M = min(12, max(1, MON))
        L = min(18, max(1, (LAT+95)/10))
        read (NJ2,'(3X,11F7.1)') (TREF(I,L,M), I=1,41)
        read (NJ2,'(3X,11F7.4)') (OREF(I,L,M), I=1,31)
      enddo
      close (NJ2)

!  Extend climatology to 100 km
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
! 1000 format(1x,'std atmos profiles: ',i3,' lat x ',i2,' mon')

      end subroutine RD_PROF
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: SET_CLD0
!
! !INTERFACE:
!
      subroutine SET_CLD0(TINIT,ODINIT,ODINDX,ALBEDO)
!
! !USES:
!
      use FastJX53cMetVars_mod      , only : T, ODCLD, SA, NCLDX
      use FastJX53cCTMparameters_mod, only : I_, J_, L_, L1_
!
      implicit none
!
! !INPUT PARAMETERS:
      integer, intent(in) :: ODINDX
      real*8 , intent(in) :: TINIT(L_), ODINIT(L_), ALBEDO
!
! !DESCRIPTION:
!  Routine to set cloud and surface properties: now loads the 
!  input prof's.
!           
!  {\bf This subroutine will need to be customized.
!  It is separate from aerosols since it comes from met fields.}
!
! !LOCAL VARIABLES:
      integer  I, J, L
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
      do J = 1,J_
         do I = 1,I_
            do L = 1,L_
               ODCLD(I,J,L) = ODINIT(L)
               T    (I,J,L) = TINIT(L)
            enddo
            ODCLD(I,J,L1_)  = 0.d0
            SA   (I,J)     = ALBEDO
         enddo
      enddo
      NCLDX(:,:,:) = ODINDX

      return

      end subroutine SET_CLD0
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: SET_AER0
!
! !INTERFACE:
!
      subroutine SET_AER0
!
! !USES:
      use FastJX53cMetVars_mod      , only : NAER1, NAER2, NAER3
      use FastJX53cMetVars_mod      , only : DAER1, DAER2, DAER3
      use FastJX53cCTMparameters_mod, only : I_, J_, L_
!
      implicit none
!
! !DESCRIPTION:
!  Set up aerosols.
!
! {\bf Customize for climatology or CTM as source.}
!
! !LOCAL VARIABLES:
      integer I, J, L, K
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

!  Add aerosol OD and indices
!    DAERn(L=1:L_) = aerosol/cloud Optical Depth (OD) within CTM layer L
!    DAERn(L1_=L_+1)= aerosol/cloud OD in layer above CTM (L1_)

         NAER1(:,:,:) = 0
         NAER2(:,:,:) = 0
         NAER3(:,:,:) = 0
!      NAER1(:,:,:) =  3    !  Black carbon absorber
!      NAER2(:,:,:) = 10    !  Water Cloud (Deirmenjian 8 micron)
!      NAER3(:,:,:) = 14    !  Irregular Ice Cloud (Mishchenko)
      do J = 1,J_
      do I = 1,I_
        do L = 1,L_+1
          DAER1(I,J,L) = 0.d0
          DAER2(I,J,L) = 0.d0
          DAER3(I,J,L) = 0.d0
        enddo
      enddo
      enddo

      return

      end subroutine SET_AER0
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: SET_ATM0
!
! !INTERFACE:
!
      subroutine SET_ATM0
!
! !USES:
      use FastJX53cMetVars_mod      , only : OREF, TREF, DM, DO3, TJ, ZH, doModelClima
      use FastJX53cMetVars_mod      , only : ETAA, ETAB, PMEAN, MONTH, YDGRD
      use FastJX53cCTMparameters_mod, only : I_, J_, L_,  L1_
!
      implicit none
!
! !INPUT PARAMETERS:

!
! !DESCRIPTION:
!  Routine to set up atmospheric profiles required by Fast-J2 using a
!  doubled version of the level scheme used in the CTM. First pressure
!  and z* altitude are defined, then O3 and T are taken from the supplied
!  climatology and integrated to the CTM levels (may be overwritten with
!  values directly from the CTM, if desired). \newline
!                                       Oliver (04/07/99) \& MJP (7/05)
!
!  \begin{verbatim}
!     PJ       Pressure at boundaries of model levels (hPa)
!     MASFAC   Conversion factor for pressure to column density
!
!     TJ       Temperature profile on model grid
!     DM       Air column for each model level (molecules.cm-2)
!     DO3      Ozone column for each model level (molecules.cm-2)
!     ZH       Altitude of boundaries of model levels (cm)
!     PSTD     Approximate pressures of levels for supplied climatology
!  \end{verbatim}
!
! !LOCAL VARIABLES:
      integer  I, J, K, L, M, N
      real*8   PSTD(52),OREF2(51),TREF2(51),PJ(L1_+1)
      real*8   DLOGP,F0,T0,PB,PC,XC,MASFAC,SCALEH
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
!  Select appropriate month
      M = max(1,min(12,MONTH))

!  Mass factor - delta-Pressure (mbars) to delta-Column (molecules.cm-2)
      MASFAC = 100.d0*6.022d+23/(28.97d0*9.8d0*10.d0)

      do J = 1,J_

!  Select appropriate latitudinal profiles
          N = max(1, min(18, (int(YDGRD(J))+99)/10 ))
!  Temporary zonal arrays for climatology data
          do K = 1,51
            if (.not. doModelClima) then ! Jules Kouatchou 08/31/06
               OREF2(K) = OREF(K,N,M)
               TREF2(K) = TREF(K,N,M)
            end if                       ! Jules Kouatchou 08/31/06
          enddo

        do I = 1,I_

!  Apportion O3 and T on supplied climatology z* levels onto CTM levels
!  with mass (pressure) weighting, assuming constant mixing ratio and
!  temperature half a layer on either side of the point supplied.
!  L1_ = L_+1:  
!       PJ(L=1:L1_) = pressure at CTM layer edge, PJ(L1_+1)=0 (top-of-atmos)

!Begin Jules Kouatchou 08/24/2006
!           PJ(L1_+1) = 0.d0
!         do K = 1,L1_
!           PJ(K) = ETAA(K) + ETAB(K)*PMEAN(I,J)
!         enddo
!End Jules Kouatchou 08/24/2006

!  Set up pressure levels for O3/T climatology - assume that value
!  given for each 2 km z* level applies from 1 km below to 1 km above,
!  so select pressures at these boundaries. Surface level values at
!  1000 mb are assumed to extend down to the actual P(ILNG,JLAT).
!
           PSTD(1) = max(PJ(1),1000.d0)
           PSTD(2) = 1000.d0 * 10.d0**(-1.d0/16.d0)
           DLOGP   = 10.d0**(-2.d0/16.d0)
         do K = 3,51
           PSTD(K) = PSTD(K-1)*DLOGP
         enddo
           PSTD(52)  = 0.d0

         do L = 1,L1_
            if (.not. doModelClima) then ! Jules Kouatchou 08/31/06
               F0 = 0.d0
               T0 = 0.d0
            end if                       ! Jules Kouatchou 08/31/06
           do K = 1,51
             PC   = min(PJ(L),PSTD(K))
             PB   = max(PJ(L+1),PSTD(K+1))
             if (PC .gt. PB) then
               XC = (PC-PB)/(PJ(L)-PJ(L+1))
               if (.not. doModelClima) then ! Jules Kouatchou 08/31/06
                  F0 = F0 + OREF2(K)*XC
                  T0 = T0 + TREF2(K)*XC
               end if                       ! Jules Kouatchou 08/31/06
             endif
           enddo
           if (.not. doModelClima) then ! Jules Kouatchou 08/31/06
              TJ(I,J,L)  = T0
              DM(I,J,L)  = (PJ(L)-PJ(L+1))*MASFAC
              DO3(I,J,L) = F0*1.d-6*DM(I,J,L)
           end if                       ! Jules Kouatchou 08/31/06
         enddo

!  Calculate effective altitudes using scale height at each level
           ZH(I,J,1) = 0.d0
         do L = 1,L_
           SCALEH      = 1.3806d-19*MASFAC*TJ(I,J,L)
           ZH(I,J,L+1) = ZH(I,J,L) -( LOG(PJ(L+1)/PJ(L)) * SCALEH )
         enddo

        enddo
      enddo

      return
      end subroutine SET_ATM0
!EOC
!-------------------------------------------------------------------------

      end module FastJX53cCTMspecific_mod
