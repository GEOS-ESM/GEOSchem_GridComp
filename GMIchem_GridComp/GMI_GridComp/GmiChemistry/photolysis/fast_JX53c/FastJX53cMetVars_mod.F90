!-----------------------------------------------------------------------
!                     NASA/GSFC - SIVO ASTG Code 610.3                 !
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: FastJX53cMetVars_mod
!
! !INTERFACE:
      module FastJX53cMetVars_mod
!
! !USES:
      use FastJX53cCTMparameters_mod, only : I_, J_, L_, LWE_, L1_, NTR_

      IMPLICIT NONE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      private
      public  :: InitializationFastjMetVars
!
! !PUBLIC DATA MEMBERS:
!
      public  :: doModelClima
      public  :: P     , T     , OD   , PJ
      public  :: TJ    , DM    , DO3  , ZH
      public  :: SA    , STT   , PMEAN
      public  :: TREF  , OREF
      public  :: XGRD  , XDGRD
      public  :: YGRD  , YDGRD
      public  :: ETAA  , ETAB
      public  :: MONTH , NSLAT , NSLON
      public  :: DAER1 , DAER2 , DAER3, ODCLD
      public  :: NAER1 , NAER2 , NAER3, NCLDX
      public  :: AREAXY, TCNAME

      logical          :: doModelClima
      integer          :: MONTH
      integer          :: NSLAT                   ! Latitude(J) index of current column
      integer          :: NSLON                   ! Longitude(I) index of current column
      real*8           :: TREF(51,18,12)
      real*8           :: OREF(51,18,12)
      character*10     :: TCNAME(NTR_)

      real*8 , pointer :: PJ(:)        => null()  !  Column pressure
      real*8 , pointer :: P(:,:)       => null()  !  Surface pressure
      real*8 , pointer :: T(:,:,:)     => null()  !  Temperature profile
      real*8 , pointer :: OD(:,:,:)    => null()  !  Optical Depth profile
      real*8 , pointer :: XGRD(:)      => null()  !  Longitude (midpoint, radians)
      real*8 , pointer :: XDGRD(:)     => null()
      real*8 , pointer :: YGRD(:)      => null()  !  Latitude  (midpoint, radians)
      real*8 , pointer :: YDGRD(:)     => null() 
      real*8 , pointer :: ETAA(:)      => null()  !  Eta(a) value for level boundaries
      real*8 , pointer :: ETAB(:)      => null()  !  Eta(b) value for level boundaries
      real*8 , pointer :: AREAXY(:,:)  => null()  !  area (m^2)
      real*8 , pointer :: TJ(:,:,:)    => null()
      real*8 , pointer :: DM(:,:,:)    => null()
      real*8 , pointer :: DO3(:,:,:)   => null()
      real*8 , pointer :: ZH(:,:,:)    => null()
      real*8 , pointer :: DAER1(:,:,:) => null()
      real*8 , pointer :: DAER2(:,:,:) => null()
      real*8 , pointer :: DAER3(:,:,:) => null()
      real*8 , pointer :: ODCLD(:,:,:) => null()
      integer, pointer :: NAER1(:,:,:) => null()
      integer, pointer :: NAER2(:,:,:) => null()
      integer, pointer :: NAER3(:,:,:) => null()
      integer, pointer :: NCLDX(:,:,:) => null()
      real*8 , pointer :: PMEAN(:,:)   => null()
      real*8 , pointer :: SA(:,:)      => null()
      real*8 , pointer :: STT(:,:,:,:) => null()

!
! !DESCRIPTION:
!  Delivers $p$ , $T$, Surf Albedo, and Optical Depth from CTM to fastJX.
!  This is for standalone fast-JX ver 5.3   (6.05)
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
! !IROUTINE:  InitializationFastjMetVars
! 
! !INTERFACE:
!
      subroutine InitializationFastjMetVars()
!
      implicit none
!
! !DESCRIPTION:
!  Allocates CTM related variables.
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
      doModelClima = .false.

      Allocate(PJ(L1_+1))
      PJ     = 0.0d0
      Allocate(P(I_,J_))
      P      = 0.0d0
      Allocate(T(I_,J_,L_))
      T      = 0.0d0
      Allocate(OD(I_,J_,LWE_))
      OD     = 0.0d0

      Allocate(XGRD (I_))
      XGRD   = 0.0d0
      Allocate(XDGRD(I_))
      XDGRD  = 0.0d0

      Allocate(YGRD (J_))
      YGRD   = 0.0d0
      Allocate(YDGRD(J_))
      YDGRD  = 0.0d0

      Allocate(ETAA(L_+1))
      ETAA   = 0.0d0
      Allocate(ETAB(L_+1))
      ETAB   = 0.0d0

      Allocate(AREAXY(I_,J_))
      AREAXY = 0.0d0

      Allocate(TJ   (I_,J_,L1_))
      TJ    = 0.0d0
      Allocate(DM   (I_,J_,L1_))
      DM    = 0.0d0
      Allocate(DO3  (I_,J_,L1_))
      DO3   = 0.0d0
      Allocate(ZH   (I_,J_,L1_))
      ZH    = 0.0d0

      Allocate(DAER1(I_,J_,L1_))
      DAER1 = 0.0d0
      Allocate(DAER2(I_,J_,L1_))
      DAER2 = 0.0d0
      Allocate(DAER3(I_,J_,L1_))
      DAER3 = 0.0d0
      Allocate(ODCLD(I_,J_,L1_))
      ODCLD = 0.0d0

      Allocate(NAER1(I_,J_,L1_))
      Allocate(NAER2(I_,J_,L1_))
      Allocate(NAER3(I_,J_,L1_))
      Allocate(NCLDX(I_,J_,L1_))

      Allocate(PMEAN(I_,J_))
      PMEAN = 0.0d0
      Allocate(SA   (I_,J_))
      SA    = 0.0d0

      Allocate(STT(I_,J_,L_,NTR_))
      STT   = 0.0d0

      return

      end subroutine InitializationFastjMetVars
!EOC
!-----------------------------------------------------------------------

 end module FastJX53cMetVars_mod
