!=============================================================================
! CODE DEVELOPER
!    A. Oloso, GSFC, 9/2005
!   (Adapted from "program standalone" in fastJX.f of fast-JX53b distributed by
!    Michael Prather, UCI)
!
! FILE
!   fast_JX53b_control.F
!
! MODULES
!    Fast_JX53b
!
! ROUTINES
!    Control_Fast_JX53b
!
!=============================================================================
!-----------------------------------------------------------------------------
! MODULE
!    Fast_JX53b
!
! DESCRIPTION
!    This module contains one public routine Control_Fast_JX53b.
!    Everything else in the module is private.
!
! ROUTINE
!   Control_Fast_JX53b
!
! DESCRIPTION
!   This routine acts as the interface between the GMI model and Fast-JX53b.
!
! ARGUMENTS
!
!   INPUTS
!     k1            : first interface index [lpar = k2 - k1 + 1]
!     k2            : last  interface index [lpar = k2 - k1 + 1]
!     chem_mask_khi : number of chemistry levels     [JVL_]
!     num_qjs       : number of photolysis reactions [JVN_]
!     month_gmi     : number of month (1- 12)        [month]
!     jday          : day    of year  (1-365)        [iday]
!     time_sec      : time of day in model (s, GMT)  [GMTAU (hrs)]
!     londeg_i      : longitude (midpoint, deg) [XDGRD]
!     latdeg_j      : latitude  (midpoint, deg) [YDGRD]
!     pctm_ij       : surface pressure (mb) [PMEAN]
!     kel_ij        : temperature at box centers (degK) [T]
!     optdepth_ij   : optical depth in box   (unitless) [ODCLD]
!     surf_alb_ij   : surface albedo     (fraction 0-1) [SA]
!     ozone_ij      : mixing ratio for ozone [DO3, if use model ozone]
!     sza_ij        : solar zenith angle
!   OUTPUTS
!     qjgmi_ij      : jvalues at grid centers (s^-1?)   [zpj]
!
!-----------------------------------------------------------------------------
!  >>>>>>>>>>>>>>>>current code revised to JX ver 5.3b (7/05)<<<<<<<<<<<<
! version 5.3b changes include:
!      new data files for specral Xsection and mie-scattering.
!      add sub-layers (JXTRA) to thick cloud/aerosol layers,
!           sets up log-spaced sub-layers of increasing thickness ATAU
!      correction 'b' does massive clean up of the linking code,
!           now the only subroutine that has access to CTM arrays is PHOTOJ
!           Also, the access to the cmn_JVdat.h is 'read-only' after init.
!           This should enable safe openMP/MPI coding.
!
! common files and what they mean:
!   parm_CTM.h  dimensions & params for code (CTM and fast-JX)
!   parm_MIE.h  dimensions for mie code variables.
!   cmn_metdat.h  CTM 3-D arrays, time of day, grid,  etc.
!   cmn_JVdat.h   Xsects, Mie, etc., (initialized and then read-only)
!
!<<<<<<<<<<<<<<<<<<<<<begin CTM-specific subroutines<<<<<<<<<<<<<<<<<<<<
! subroutines:
!
!     SET_ATM(GMTAU)
!           set ups atmosphere (p,T,O3,airmass, etc) for time GMTAU
!              COMMON BLOCKS: cmn_metdat.h
!
!     SET_AER(GMTAU)
!              set ups aerosols for time GMTAU = DUMMY
!
!     SET_CLD(GMTAU)
!           set ups clouds for time GMTAU = DUMMY
!
!     INPHOT:  Init. photolysis rate, called once by CHMSET
!              COMMON BLOCKS: cmn_metdat.h, cmn_JVdat.h
!              Input files: ECT42_grid.dat
!
!     RD_JS(NJ1,NAMFIL):  Read labels of photo. rates, called once by INPHOT.
!              COMMON BLOCKS: cmn_metdat.h, cmn_JVdat.h
!              Input files: ratj.dat
!
!     RD_PROF(NJ2,NAMFIL):  Read T & O3 climatology, called once by INPHOT.
!              COMMON BLOCKS: cmn_metdat.h
!              Input files: atmos_std.dat
!
!     SET_CLD0(TINIT,ODINIT,ODINDX,ALBEDO)
!              Initialize cloud and surface properties, called by MAIN.
!              COMMON BLOCKS: cmn_metdat.h
!
!     SET_AER0:  Iniitalize (climatology) aerosol OD and types (3 arrays)
!                     called by MAIN, CHMSET
!              COMMON BLOCKS: cmn_metdat.h
!
!     SET_ATM0:  Initialize climatologies for T & O3, set up atmospheric profiles if model ozone is not used
!              COMMON BLOCKS: cmn_metdat.h
!
!<<<<<<<<<<<<<<<<<<<<<begin CTM-fastJX linking subroutines<<<<<<<<<<<<<<
!
!     PHOTOJ(UTIME,IDAY,ILNG,JLAT, SZA,ZPJ)
!              Gateway to fast-JX, Update the photolysis rates
!              COMMON BLOCKS: cmn_metdat.h, cmn_JVdat.h
!
!<<<<<<<<<<<<<<<<<<<<<begin core fast-J subroutines<<<<<<<<<<<<<<<<<<<<<
!  N.B. all these need access to cmn_JVdat.h, but do NOT write into it.
!           also have no need to access cmn_metdat.h
!
!     JRATET(PPJ,TTJ,FFF, VALJL):  Calculate J-value, called by PTOTOJ.
!              COMMON BLOCKS: cmn_JVdat.h
!
!     JP_ATM(PPJ,TTJ,DDJ,ZZJ,ZHL,ZZHT,AER,ABX,ADX,JXTRA)
!              print out atmosphere used in J-value calc.
!              COMMON BLOCKS: cmn_JVdat.h
!
!
!     RD_XXX(NJ1,NAMFIL):  Read wavelength bins, solar fluxes, Rayleigh
!             parameters, TEM-dependent X-sections, called once by INPHOT.
!              COMMON BLOCKS: cmn_JVdat.h
!              Input files: FJX_spec.dat
!
!     RD_MIE(NJ1,NAMFIL):  Set aerosols/cloud scattering, called once by INPHOT
!              COMMON BLOCKS: cmn_JVdat.h
!              Input files: FJX_scat.dat
!
!     FUNCTION FLINT(TINT,T1,T2,T3,F1,F2,F3)
!
!     SOLARZ(GMTIME,NDAY,YGRDJ,XGRDI, SZA,COSSZA,SOLFX)
!              calc SZA and Solar Flux factor for given lat/lon/UT
!
!     SPHERE(U0,RAD,ZHL,ZZHT,AMF,L1_):
!              calculate spherical geometry, air-mass factors
!
!     EXTRAL(AER,ADX,L1X,L2X,NX,JTAUMX,ATAU,ATAU0, JXTRA)
!              add sub-layers (JXTRA) to thick cloud/aerosol layers
!
!     OPMIE (KW,KM,WAVEL,ABX,AER,ADX,U0,RFLECT,AMF,JXTRA,FMEAN)
!              calculate mean actinic flux at desired levels
!              COMMON BLOCKS: cmn_JVdat.h
!
!<<<<<<<<<<<<<<<<<<<<<<<begin core scattering subroutines<<<<<<<<<<<<<<<
!
!      MIESCT (FJ,POMEGA,FZ,ZTAU,ZFLUX,ZREFL,ZU0,MFIT,ND)
!            include 'parm_MIE.h' = dimension parameters
!
!      BLKSLV (FJ,POMEGA,FZ,ZTAU,ZFLUX,ZREFL,WT,EMU,PM,PM0,M,N,MFIT,ND)
!              PARAMETER FILE: parm_MIE.h
!
!      GEN (POMEGA,FZ,ZTAU,ZFLUX,ZREFL,WT,EMU,PM,PM0,B,CC,AA,A,H,C1
!            ,M,N,MFIT,ND,ID)
!              PARAMETER FILE: parm_MIE.h
!
!      LEGND0 (X,PL,N)
!
!      MATIN4 (A)
!
!      GAUSSP (N,XPT,XWT)
!
!      EFOLD  (F0, F1, N, F)
!
!-----------------------------------------------------------------------
      module Fast_JX53b
      private
      public :: Control_Fast_JX53b
      public :: InitializeFastJX53b
      public :: GetQAA_RAAinFastJX53b
      contains
      subroutine Control_Fast_JX53b  &
     & (k1, k2, chem_mask_khi, num_qjs, month_gmi, jday, time_sec,  &
     &       sza_ij, press3c_ij, pctm_ij, kel_ij,  &
     &       optdepth_ij, surf_alb_ij, qjgmi_ij, overheadO3col_ij, ozone_ij)
!-----------------------------------------------------------------------
!  Routine to test Fast-J code by simulating model calls
!
!  >>>>>overall program revised to run with FAST-JX  J-code ver 5.3 (6/05)
!  >>>>>  now designed to use log spacing in TAU to get <0.5% accuracy
!         much better accuracy and fewer added levels vs. old linear adds.
!
!-----------------------------------------------------------------------
!
!<<<<<<<<<<<<<<<<<<<<<begin CTM-specific subroutines<<<<<<<<<<<<<<<<<<<<
!
!  these are mesy and meant to link with or be replace by the CTM code
!  SOME of the important subroutines - eg, SET_CLD need to be updated.

      implicit none
#     include "fastJX53b_parm_CTM.h"
#     include "fastJX53b_cmn_metdat.h"
#     include "fastJX53b_cmn_JVdat.h"

      integer, intent(in)  :: k1                 ! first interface index
      integer, intent(in)  :: k2                 ! last  interface index
      integer, intent(in)  :: chem_mask_khi      ! number of chemistry levels     [jpnl]
      integer, intent(in)  :: num_qjs            ! number of photolysis reactions [jppj]
      integer, intent(in)  :: month_gmi          ! number of month (1- 12)        [month]
      integer, intent(in)  :: jday
      real*8,  intent(in)  :: kel_ij(k1:k2)      ! temperature at box centers (degK) [t]
      real*8,  intent(in)  :: time_sec           ! time of day in model (s, GMT)  [gmtau (hrs)]
      real*8,  intent(in)  :: sza_ij           ! solar zenith angle
      !real*8,  intent(in)  :: londeg_i           ! longitude (midpoint, deg) [xdgrd]
      !real*8,  intent(in)  :: latdeg_j           ! latitude  (midpoint, deg) [ydgrd]
      real*8,  intent(in)  :: press3c_ij(k1:k2)  
      real*8,  intent(in)  :: pctm_ij            ! surface pressure (mb) [p]
      real*8,  intent(in)  :: optdepth_ij(k1:k2) ! optical depth in box   (unitless) [od]
      real*8,  intent(in)  :: surf_alb_ij        ! surface albedo     (fraction 0-1) [sa]
      real*8, optional, intent(in) :: ozone_ij(k1:k2)      ! mixing ratio for ozone
!-----------------------------------------------------------------------

      real*8,  intent(out) :: overheadO3col_ij(k1:k2)
      real*8,  intent(out) :: qjgmi_ij(k1:chem_mask_khi, num_qjs)

      real*8  ::  amm(k1-1:k2)
      real*8  ::  bmm(k1-1:k2)

      real*8  ZPJ(chem_mask_khi-k1+1, num_qjs)    !2-D array of J's indexed to CTM chemistry!
!      real*8 GMTAU,DELTAU, TINIT(L_),ODINIT(L_),ALBEDO,SZA
      real*8 GMTAU,DELTAU, ALBEDO,SZA
!      integer I,J,K,L,NSTEP   ,ILNG,JLAT,ODINDX, IDAY
      integer I,J,K,L,ILNG,JLAT,ODINDX, IDAY, il, ik
      real*8, parameter :: PI=3.141592653589793d0
      logical, save :: first = .true.

      ILNG                = I_
      JLAT                = J_

!      L_                  = k2 - k1 + 1
!      JVL_                = chem_mask_khi
!      JVN_                = num_qjs
!      L1_                 = L_ + 1
!      L2_                 = 2 * L1_


      MONTH               = month_gmi
      IDAY                = jday
      GMTAU               = time_sec / 3600.0d0

      !XDGRD(I_)           = londeg_i
      !YDGRD(J_)           = latdeg_j
!      amm(k1-1)           = 0.0d0
!      bmm(k1-1)           = 1.0d0
!      amm(k1:k2)          = am(:)
!      bmm(k1:k2)          = bm(:)

!      etaa(1:L1_)         = amm(k1-1:k2) * pt
!      etab(1:L1_)         = bmm(k1-1:k2)

      PJ(1)               = pctm_ij
      PJ(2:L1_)           = press3c_ij(k1:k2)
      PJ(L1_+1)           = 0.0d0

      PMEAN(I_, J_)       = pctm_ij
      T(I_, J_, 1:L_)     = kel_ij(k1:k2)
      ODCLD(I_, J_, 1:L_) = optdepth_ij(k1:k2)
      ODCLD(I_, J_, L1_)  = 0.0
      SA(I_, J_)          = surf_alb_ij
      NCLDX(I_, J_, 1:L1_) = 10

      XGRD(I_) = XDGRD(I_) *PI/180.d0
      YGRD(J_) = YDGRD(J_) *PI/180.d0

!  Initial call to Fast-J to set things up
!      if (first) then
!        first = .false.
!!-----------------------------------------------------------------------
!        call INPHOT
!!-----------------------------------------------------------------------
!!        call SET_CLD0(TINIT,ODINIT,ODINDX,ALBEDO)
!!-----------------------------------------------------------------------
!!          PMEAN(ILNG,JLAT) = P(ILNG,JLAT)
!!-----------------------------------------------------------------------
!        call SET_AER0
!!-----------------------------------------------------------------------
!      endif

      if (.not. present(ozone_ij)) then
        call SET_ATM0
      endif

!---reset the atmosphere, aerosols, clouds for the time step (now = dummy)
!-----------------------------------------------------------------------
      if (present(ozone_ij)) then
        DO3(I_, J_, 1:L_)   = ozone_ij(k1:k2)
        DO3(I_, J_, L1_)   = ozone_ij(k2)
        call SET_ATM(GMTAU)
      endif
      call SET_AER(GMTAU)
      call SET_CLD(GMTAU)
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
      call PHOTOJ(GMTAU,IDAY,ILNG,JLAT, sza_ij, ZPJ)
!-----------------------------------------------------------------------

      qjgmi_ij(:,:) = ZPJ(:,:)

!!!!
!!!! Overhead ozone column
!!!!
      overheadO3col_ij = 0.0d0
      overheadO3col_ij(k2) = DO3(I_,J_,k2)
      do il = k1, k2-1
         do ik = il+1, k2
            overheadO3col_ij(il) = overheadO3col_ij(il) + DO3(I_,J_,ik)
         end do
      end do

      return

      end subroutine Control_Fast_JX53b

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

      subroutine InitializeFastJX53b &
                 (k1, k2, chem_mask_khi, num_qjs, &
                   CrossSection_InfileName  , ScatteringData_InfileName,   &
                   PhotRateLabels_InfileName, T_O3_Climatology_InfileName)

      implicit none

      integer            , intent(in) :: k1, k2
      integer            , intent(in) :: chem_mask_khi ! number of chemistry levels [JVL_]
      integer            , intent(in) :: num_qjs       ! number of photolysis reactions [JVN_]
      character (len=128), intent(in) :: CrossSection_InfileName
!                             ! fast-J X-sections (spectral data) input file name
      character (len=128), intent(in) :: ScatteringData_InfileName
!                             ! Aerosol/cloud scattering data input file name
      character (len=128), intent(in) :: PhotRateLabels_InfileName
!                             ! Labels of photolysis rates required input file name
!                             ! keyed to chem code
      character (len=128), intent(in) :: T_O3_Climatology_InfileName
!                             ! Read in T & O3 climatology input file name
!                             ! general backup clim.
      
#     include "fastJX53b_parm_CTM.h"
#     include "fastJX53b_cmn_JVdat.h"

!      I_                  = 1
!      J_                  = 1
      L_                  = k2 - k1 + 1
      JVL_                = chem_mask_khi
      JVN_                = num_qjs
      L1_                 = L_ + 1
      L2_                 = 2 * L1_

      ! Read in the input files
      call INPHOT (CrossSection_InfileName,     &
     &             ScatteringData_InfileName,   &
     &             PhotRateLabels_InfileName,   &
     &             T_O3_Climatology_InfileName)

      ! Set up aerosol
      call SET_AER0

      return

      end subroutine InitializeFastJX53b

!-----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine GetQAA_RAAinFastJX53b (bRAA, bQAA)

      implicit none

#     include "fastJX53b_parm_CTM.h"
#     include "fastJX53b_cmn_JVdat.h"

      real*8 , intent(out) :: bRAA(:,:), bQAA(:,:)

      bRAA(:,:) = RAA(:,:)
      bQAA(:,:) = QAA(:,:)

      return

      end subroutine GetQAA_RAAinFastJX53b

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-----------------------------------------------------------------------

#     include "fastJX53b.code"

      end module Fast_JX53b
