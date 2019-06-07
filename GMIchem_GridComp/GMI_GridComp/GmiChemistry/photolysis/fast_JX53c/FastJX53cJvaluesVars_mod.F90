!-----------------------------------------------------------------------
!                     NASA/GSFC - SIVO ASTG Code 610.3                 !
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: FastJX53cJvaluesVars_mod
!
! !INTERFACE:
      module FastJX53cJvaluesVars_mod
!
! !USES:
      use FastJX53cCTMparameters_mod, only : A_, W_, X_, JVN_
!
      IMPLICIT NONE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      private
      public  :: InitializationFastjJvaluesVars
!
! !PUBLIC DATA MEMBERS:
!
      public  :: RAD   , ZZHT   , ATAU   , ATAU0
      public  :: WBIN  , WL     , FL     , QO2   , QO3, Q1D
      public  :: QQQ   , QRAYL  , TQQ
      public  :: WAA   , QAA    , PAA    , RAA   , SSA
      public  :: JFACTA
      public  :: JIND  , NRATJ  , NJVAL
      public  :: NW1   , NW2    , NAA    , JTAUMX
      public  :: TITLEA
      public  :: TITLE0
      public  :: TITLEJ, TITLEJ2, TITLEJ3, JLABEL

      real*8       :: RAD, ZZHT, ATAU, ATAU0
      real*8       :: WBIN(W_+1)  , WL(W_)     , FL(W_)
      real*8       :: QO2(W_,3)   , QO3(W_,3)  , Q1D(W_,3)
      real*8       :: QQQ(W_,2,X_), QRAYL(W_+1), TQQ(3,X_)
      real*8       :: WAA(4,A_)   , QAA(4,A_)  , PAA(8,4,A_), RAA(4,A_), SSA(4,A_)
!      real*8       :: JFACTA(JVN_)
!      integer      :: JIND(JVN_)
      integer      :: NRATJ, NJVAL, NW1, NW2, NAA, JTAUMX
      character*20 :: TITLEA(A_)
      character*78 :: TITLE0
      character*7  :: TITLEJ(X_), TITLEJ2, TITLEJ3
!      character*7  :: JLABEL(JVN_)

      real*8     , pointer :: JFACTA(:) => null()
      integer    , pointer :: JIND  (:) => null()
      character*7, pointer :: JLABEL(:) => null()
!
! !DESCRIPTION:
!  Declare J-values related variables.
!
! !AUTHOR:    
! Michael Prather
!
! !HISTORY:
!
!EOP
!-----------------------------------------------------------------------
   CONTAINS
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:  InitializationFastjJvaluesVars
!
! !INTERFACE:
!
      subroutine InitializationFastjJvaluesVars()
!
      implicit none
!
! !DESCRIPTION:
!  Allocates JvaluesVars related variables.
!
! !AUTHOR:
!  Jules Kouatchou
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC
      allocate( JFACTA(JVN_))
      allocate( JIND  (JVN_))
      allocate( JLABEL(JVN_))

      JFACTA(:) = 0.0d0

      return

      end subroutine InitializationFastjJvaluesVars
!EOC
!-----------------------------------------------------------------------

 end module FastJX53cJvaluesVars_mod
