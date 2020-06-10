!------------------------------------------------------------------------------
!  GMI_fjx_mod module is to get old (FasyJX 6.5 and earlier aerosol parameters
!-----------------------------------------------------------------------
!
! !INTERFACE:
!
      MODULE GMI_CMN_MOD
!
      implicit none
      public
!
      integer, parameter :: NK_=4, NgmiA_=56, Nwaer_=5, NgmiRH_=7
!
      real*8 :: gmiWAA (NK_, NgmiA_)
      real*8 :: gmiQAA (NK_, NgmiA_)
      real*8 :: gmiRAA (NK_, NgmiA_) 
      real*8 :: gmiDAA (NK_, NgmiA_) 
      real*8 :: gmiSAA (NK_, NgmiA_)
      real*8 :: gmiPAA (8, NK_, NgmiA_)
!
      real*8 :: gmiRW (Nwaer_, NgmiRH_)
      real*8 :: gmiQW (Nwaer_, NgmiRH_)
!
      real*8, parameter :: tRH(NgmiRH_)=[0.0,0.50,0.70,0.80,0.90,0.95,0.99]
!
      integer*8, parameter :: iDRYwaer(Nwaer_)=[22,29,36,43,50]
!
      END MODULE GMI_CMN_MOD
!EOP
