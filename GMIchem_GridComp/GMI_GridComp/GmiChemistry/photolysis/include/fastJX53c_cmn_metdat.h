!-----------------------------------------------------------------------
!    include 'cmn_metdat.f'  for fast-JX code v5.3 (prather 6/05)
!
!        needs 'parm_ctm.f' for dimensions
!        delivers p, T, Surf Albedo, and Optical Depth from CTM to fastJX
!        >>>>this is for standalone fast-JX ver 5.3   (6.05)
!-----------------------------------------------------------------------
!
      real*8  P(I_,J_)         !  Surface pressure
      real*8  T(I_,J_,MAX_L_)      !  Temperature profile
!      real*8  OD(I_,J_,LWE_)   !  Optical Depth profile
!
      real*8  XGRD(I_)         !  Longitude (midpoint, radians)
      real*8  XDGRD(I_)
      real*8  YGRD(J_)         !  Latitude  (midpoint, radians)
      real*8  YDGRD(J_)
      real*8  PJ(MAX_L_+1)    
      real*8  ETAA(MAX_L_+1)       !  Eta(a) value for level boundaries
      real*8  ETAB(MAX_L_+1)       !  Eta(b) value for level boundaries
      real*8  AREAXY(I_,J_)    !  area (m^2)
      integer  MONTH
      integer  NSLAT         ! Latitude(J) index of current column
      integer  NSLON         ! Longitude(I) index of current column

      real*8, dimension(I_,J_,MAX_L1_) :: TJ, DM, DO3, ZH
      real*8, dimension(I_,J_,MAX_L1_) :: DAER1, DAER2, DAER3, ODCLD
      integer,dimension(I_,J_,MAX_L1_) :: NAER1, NAER2, NAER3, NCLDX
      real*8, dimension(I_,J_)     :: PMEAN, SA
!
      real*8   STT(I_,J_,MAX_L_,NTR_)
      real*8   TREF(51,18,12),OREF(51,18,12)
      character*10  TCNAME(NTR_)

      real*8,  pointer :: optDust(:,:)
      real*8,  pointer :: optAer (:,:)
      real*8,  pointer :: AER    (:,:)
      integer, pointer :: ADX    (:,:)
      real*8, dimension(MAX_L1_) :: ODCOL

      common/metdatAerosol/optDust, optAer, AER, ADX

!      common/metdat/P,T,OD, STT,
      common/metdat_JX53c/P,T,STT,  &
     &  XGRD,YGRD,XDGRD,YDGRD,ETAA,ETAB,AREAXY,  &
     &  TREF,OREF,PMEAN,SA, MONTH,NSLAT,NSLON, TCNAME
!
      common /jvdatIJ_JX53c/TJ,DM,DO3,ZH, PJ, &
     &                  DAER1,DAER2,DAER3,ODCLD,  &
     &                  NAER1,NAER2,NAER3,NCLDX, &
     &                  ODCOL
!-----------------------------------------------------------------------

