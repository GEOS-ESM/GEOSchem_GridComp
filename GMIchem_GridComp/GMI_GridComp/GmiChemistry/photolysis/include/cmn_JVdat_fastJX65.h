      !-----------------------------------------------------------------------
      !    include 'cmn_JVdat.h'  for fast-JX code v 6.0+ (prather 4/07)
      !
      ! NB - ALL of these common variables are set paramters,
      !    They are NOT to be used as variables for a local solution
      !    Thus this entire set is 'in' only after it is initialized
      !-----------------------------------------------------------------------
      real (8) :: RAD, ZZHT, ATAU, ATAU0  
      real (8) :: WBIN (WX_ + 1)  , QBC(WX_)
      real (8) :: WL (WX_), FL (WX_), QO2 (WX_, 3), QO3 (WX_, 3), &
                  Q1D (WX_, 3)
      real (8) :: QQQ (WX_, 2, X_), QRAYL (WX_ + 1), TQQ (3, X_)  
      real (8) :: WAA (5, A_), QAA (5, A_), PAA (8, 5, A_), SAA (5, A_)
      real (8) :: DAA (A_)
      real (8) :: RAA(5, A_)
      real (8) :: JFACTA (MAX_JVN_)  

      real (8) :: UMAER (3, 5, 21, 33)  
      integer :: ncat_acetone_trop, ncat_acetone_stra
      integer :: ncat_met_vinyl_ketone, ncat_met_ethyl_ketone
      integer :: ncat_methyl_glyoxal
      integer :: jppj
      integer :: JIND (MAX_JVN_), NRATJ, NJVAL, NW1, NW2, NAA, JTAUMX  
      character (len=20) :: TITLAA (A_), TITLUM (33)  
      character (len=78) :: TITLE0  
      character (len=7) :: TITLEJ (X_), TITLEJ2, TITLEJ3, JLABEL (MAX_JVN_)  
!
      common / jvchem / JFACTA, JIND, NRATJ, JLABEL, &
                        ncat_acetone_trop,ncat_acetone_stra, &
                        ncat_met_vinyl_ketone, ncat_met_ethyl_ketone, &
                        ncat_methyl_glyoxal, &
                        jppj
      common /jvdat/WBIN,WL,FL,QO2,QO3,Q1D,QQQ,QRAYL,TQQ,QBC, &
     &        ZZHT,ATAU,ATAU0, WAA,QAA,PAA,SAA,RAA,DAA,RAD, UMAER, &
     &        JTAUMX, NJVAL,NW1,NW2,NAA ,TITLE0,TITLEJ,TITLAA,TITLUM

