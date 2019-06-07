!----jv_cmn.h---COMMON BLOCKS for new FAST-J code (wild/prather 7/99)
!
!  Parameters
!  ----------
!
!     NB    Number of levels in CTM plus one for above model top
!     NC    Number of levels in the fundamental Fast-J grid
!     NS    Maximum number of species which require J-values calculating
!     NW    Maximum number of wavelength bins that can be used
!     NP    Maximum number of aerosol/cloud types that can be used
!     NH    Maximum number of Herzberg X-sections that can be used
!     MX    Number of aerosol/cloud types supplied from CTM
!
!                                        Note: THETA(NL) no longer used
!
! HISTORY
!   * February 24, 2004 - Jules Kouatchou
!     Modified the value of NP (from 21 to 56)
!     Modified the value of MX (from 3 to 35)
!-----------------------------------------------------------------------
      INTEGER      NB, NC, NS, NW, NP, NH, MX
!      PARAMETER   (NB=LPAR+1, NC=2*NB, NS=61, NW=18, NP=21, NH=7, MX=3)
!      PARAMETER   (NS=61, NW=18, NP=21, NH=7, MX=3)
      PARAMETER   (NS=61, NW=18, NP=56, NH=7, MX=35)
      INTEGER, parameter :: NB_MAX = LPAR_MAX+1   ! Added for GMI  {AOO, 8/04}
      INTEGER, parameter :: NC_MAX = 2*NB_MAX     ! Added for GMI  {AOO, 8/04}
      CHARACTER*20 TITLEA(NP)
      CHARACTER*78 TITLE0
      CHARACTER*7  TITLEJ(3,NS), jlabel(jppj_max), hzlab(nh)
      INTEGER jind(jppj_max),jadsub(nc_max),nhz,hzind(nh)
      INTEGER NJVAL,NW1,NW2,MIEDX,NAA,NLBATM,npdep,jpdep(NS)
      REAL*8 TJ,PJ,DM,DO3,Z,AER,AMF,RAD,RFLECT,SZA,U0,TANHT,ZZHT
      REAL*8 WBIN,WL,FL,QO2,QO3,Q1D,QQQ,QRAYL,TQQ,FFF,VALJ,WAA,QAA,PAA
      REAL*8 RAA,SSA,TREF,OREF,BREF,QBC,DBC,zpdep(NW,3)
      REAL*8 dtaumax,szamax,zj(jpnl_max,jppj_max),jfacta(jppj_max)
      REAL*8 dtausub,dsubdiv
      integer :: ncat_acetone_trop,ncat_acetone_stra  ! Added for GMI to {AOO, 8/04}
      integer :: ncat_met_vinyl_ketone                   ! Added for GMI to {AOO, 8/04}
      integer :: ncat_met_ethyl_ketone                   ! Added for GMI to {AOO, 8/04}
      integer :: ncat_methyl_glyoxal                   ! Added for GMI to {AOO, 8/04}
      COMMON /TITLS/TITLE0,TITLEJ,TITLEA
      COMMON /ATMOS/TJ(NB_MAX),PJ(NB_MAX+1),DM(NB_MAX),DO3(NB_MAX),  &
     &              DBC(NB_MAX),Z(NB_MAX),AER(MX,NB_MAX),  &
     &              AMF(NB_MAX,NB_MAX),RAD,RFLECT,SZA,U0,TANHT,ZZHT
      COMMON /CCWVL/WBIN(NW+1),WL(NW),FL(NW),QO2(NW,3),QO3(NW,3),  &
     &              Q1D(NW,3),QQQ(NW,2,NS-3),QRAYL(NW+1),TQQ(3,NS),  &
     &              FFF(NW,jpnl_max),VALJ(NS),WAA(4,NP),QAA(4,NP),  &
     &              PAA(8,4,NP),RAA(4,NP),SSA(4,NP),QBC(NW),  &
     &              NJVAL,NW1,NW2,MIEDX(MX),NAA,NLBATM
      COMMON /CLIM/ TREF(51,18,12),OREF(51,18,12),BREF(51)
      COMMON /JCNTR/dtaumax,szamax
      COMMON /JVALS/zj,jfacta,zpdep,npdep,jpdep,jind,jlabel
      COMMON /JVSUB/jadsub,dtausub,dsubdiv
      COMMON /dims2/ NB, NC            ! Added for GMI {AOO, 8/04}
      common /acetone_cat/ ncat_acetone_trop,ncat_acetone_stra  ! Added for GMI to {AOO, 8/04}
      common /met_vinyl_ketone/ ncat_met_vinyl_ketone   ! Added for GMI to {AOO, 8/04}
      common /met_ethyl_ketone/ ncat_met_ethyl_ketone                   ! Added for GMI to {AOO, 8/04}
      common /methyl_glyoxal/ ncat_methyl_glyoxal                   ! Added for GMI to {AOO, 8/04}
!-----------------------------------------------------------------------
