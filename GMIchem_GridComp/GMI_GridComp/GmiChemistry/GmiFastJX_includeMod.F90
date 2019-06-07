
      module GmiFastJX_includeMod

      USE GmiFastJX_ParametersMod

      private
      PUBLIC  :: t_fastJXbundle

      !!fast_JX_cmn_t.h
      !PUBLIC  :: nslat, nslon
!
      !!fast_JX_jv_mie.h
      !PUBLIC  :: A,B,C1,H,AA,CC,S,W,U1,V1,WT,EMU,PM,PM0,POMEGA
      !PUBLIC  :: ZTAU,FZ,FJ,DD,RR,ZREFL,ZFLUX,RADIUS,ZU0
      !PUBLIC  :: ND,N,M,MFIT
!
      !!fast_JX_cmn_w.h
      !PUBLIC  :: P,T,OD,SA, optaer, optdust
!!
      !fast_JX_cmn_h.h
      !!PUBLIC  :: lpar, jpnl, jppj
      !PUBLIC  :: xgrd,ygrd,ydgrd,etaa,etab
      !PUBLIC  :: tau,month,iday
      !PUBLIC  :: do_model_clima
!
      !!fast_JX_jv_cmn.h
      !PUBLIC  :: ncat_acetone_trop,ncat_acetone_stra  
      !!PUBLIC  :: ncat_met_vinyl_ketone   
      !PUBLIC  :: ncat_met_ethyl_ketone                   
      !PUBLIC  :: ncat_methyl_glyoxal 
      !PUBLIC  :: dtaumax,szamax
      !PUBLIC  :: zj,jfacta,zpdep,npdep,jpdep,jind,jlabel
      !PUBLIC  :: jadsub,dtausub,dsubdiv
      !PUBLIC  :: TITLE0,TITLEJ,TITLEA
      !PUBLIC  :: WBIN, WL, FL, QO2, QO3, Q1D, QQQ, QRAYL, TQQ, FFF
      !PUBLIC  :: VALJ, WAA, QAA, PAA, RAA, SSA, QBC
      !PUBLIC  :: TJ, PJ, DM, DO3, DBC, Z, AER, AMF
      !PUBLIC  :: TREF, OREF, BREF
      !PUBLIC  :: NJVAL,NW1,NW2,MIEDX,NAA,NLBATM
      !PUBLIC  :: RAD,RFLECT,SZA,U0,TANHT,ZZHT
!
      !PUBLIC  :: LPAR_MAX, JPNL_MAX, JPPJ_MAX, IPAR, JPAR, LDEG45
      !!PUBLIC  :: NB, NC, NS, NW, NP, NH, MX
      !PUBLIC  :: M__, N__, NL

# include "gmi_AerDust_const.h"
# include "setkin_par.h"

TYPE t_fastJXbundle
!----------------------------------------------------------------------------
! fast_JX_jv_cmn.h
!     NB    Number of levels in CTM plus one for above model top
!     NC    Number of levels in the fundamental Fast-J grid
!-----------------------------------------------------------------------
      INTEGER      NB, NC
      integer :: ncat_acetone_trop,ncat_acetone_stra ! Added for GMI to {AOO, 8/04}
      integer :: ncat_met_vinyl_ketone               ! Added for GMI to {AOO, 8/04}
      integer :: ncat_met_ethyl_ketone               ! Added for GMI to {AOO, 8/04}
      integer :: ncat_methyl_glyoxal                 ! Added for GMI to {AOO, 8/04}
      CHARACTER*20 TITLEA(NP)
      CHARACTER*78 TITLE0
      CHARACTER*7  TITLEJ(3,NS), jlabel(jppj_max), hzlab(nh)
      INTEGER jind(jppj_max),jadsub(nc_max),nhz,hzind(nh)
      INTEGER NJVAL,NW1,NW2,MIEDX(MX),NAA,NLBATM,npdep,jpdep(NS)
      REAL*8 dtaumax,szamax,zj(jpnl_max,jppj_max),jfacta(jppj_max)
      REAL*8 dtausub,dsubdiv
      REAL*8  :: RAD,RFLECT,SZA,U0,TANHT,ZZHT
      REAL*8  :: TREF(51,18,12),OREF(51,18,12),BREF(51)
      REAL*8  :: TJ(NB_MAX),PJ(NB_MAX+1),DM(NB_MAX),DO3(NB_MAX),  &
     &              DBC(NB_MAX),Z(NB_MAX),AER(MX,NB_MAX),  &
     &              AMF(NB_MAX,NB_MAX)
      REAL*8 zpdep(NW,3)
      REAL*8  :: WBIN(NW+1),WL(NW),FL(NW),QO2(NW,3),QO3(NW,3),  &
                 Q1D(NW,3),QQQ(NW,2,NS-3),QRAYL(NW+1),TQQ(3,NS),  &
                 FFF(NW,jpnl_max),VALJ(NS),WAA(4,NP),QAA(4,NP),  &
                 PAA(8,4,NP),RAA(4,NP),SSA(4,NP),QBC(NW)

!----------------------------------------------------------------------------
!fast_JX_cmn_h.h

      integer  lpar
      integer  jpnl, jppj
      real*8  xgrd(ipar)         !  Longitude (midpoint, radians)
      real*8  ygrd(jpar)         !  Latitude  (midpoint, radians)
      real*8  ydgrd(jpar)        !  Latitude  (midpoint, degrees)
      real*8  etaa(lpar_max+1)   !  Eta(a) value for level boundaries
      real*8  etab(lpar_max+1)   !  Eta(b) value for level boundaries
      real*8  tau                !  Time of Day (hours, GMT)
      integer month              !  Number of month (1-12)
      integer iday               !  Day of year
      logical do_model_clima     ! determines if climatology data come
                                 ! from the model.

!----------------------------------------------------------------------------
! fast_JX_cmn_w.h

      real*8  P(ipar,jpar)                   ! Surface pressure
      real*8  T(ipar,jpar,lpar_max+1)        ! Temperature profile
!     
      real*8  OD(ipar,jpar,lpar_max)         ! Optical Depth profile
      real*8  SA(ipar,jpar)                  ! Surface Albedo
!
      real*8  optaer(lpar_max,NSADaer*nrh_b) ! Opt. Depth profile for aerosols
      real*8  optdust(lpar_max,NSADdust)     ! Opt. Depth profile for dust
!
!----------------------------------------------------------------------------
! fast_JX_jv_mie.h

      REAL*8  :: ZREFL,ZFLUX,RADIUS,ZU0
      INTEGER :: ND,N,M,MFIT

      REAL*8  :: A(M__), B(M__,M__), C1(M__), H(M__), AA(M__,M__)
      REAL*8  :: CC(M__,M__), S(M__,M__), W(M__,M__), U1(M__,M__)
      REAL*8  :: V1(M__), WT(M__), EMU(M__), PM(M__,2*M__), PM0(2*M__)
      REAL*8  :: POMEGA(2*M__,N__), ZTAU(N__), FZ(N__), FJ(N__)
      REAL*8  :: DD(M__,M__,N__), RR(M__,N__)

!----------------------------------------------------------------------------
!fast_JX_cmn_t.h

      integer nslat               !  Latitude of current profile point
      integer nslon               !  Longitude of current profile point 

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!solar cycle capability

       REAL*8 :: fjx_solar_cycle_param(NW)

!----------------------------------------------------------------------------

END TYPE t_fastJXbundle
      end module GmiFastJX_includeMod
