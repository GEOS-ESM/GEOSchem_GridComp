      !-----------------------------------------------------------------------
      !    include 'cmn_metdat.f'  for fast-JX code v5.3+ (prather 6/05)
      !
      !        needs 'parm_ctm.f' for dimensions
      !        delivers p, T, Surf Albedo, and Optical Depth from CTM to fastJ
      !        >>>>this is for standalone fast-JX ver 5.3   (6.05)
      !-----------------------------------------------------------------------
      real (8) :: P (I_, J_)           !  Surface pressure
      real (8) :: T (I_, J_, MAX_L_)   !  Temperature profile
      real (8) :: XGRD (I_)            !  Longitude (midpoint, radians)
      real (8) :: XDGRD (I_)           !  Latitude  (midpoint, radians)
      real (8) :: YGRD (J_)  
      real (8) :: YDGRD (J_)  
      real (8) :: ETAA (MAX_L1_)       !  Eta(a) value for level boundaries
      real (8) :: ETAB (MAX_L1_)       !  Eta(b) value for level boundaries
      real (8) :: AREAXY (I_, J_)      !  area (m^2)
      integer :: MONTH  
      integer :: NSLAT                 ! Latitude(J) index of current column
      integer :: NSLON                 ! Longitude(I) index of current column

      !----------------------------------------------------------------------
      ! Note that the climatology for O3, T, clouds, etc include a layer abo
      ! the CTM (L1_=L_+1) for calculating the J's
      !----------------------------------------------------------------------

      real (8), dimension (MAX_L1_)         :: PJ 
      real (8), dimension (MAX_L1_)         :: RELH  !  Rel Hum (0.00 to 1.00)
      real (8), dimension (I_, J_, MAX_L1_) :: TJ, DM, DO3, ZH  
      real (8), dimension (I_, J_, MAX_L1_) :: CLDLWP, AER1P, AER2P  
      integer , dimension (I_, J_, MAX_L1_) :: CLDNDX, AER1N, AER2N  
      integer , dimension (MX) :: AERindex

      real (8), dimension (I_, J_)   :: PMEAN, SA  
      !
      real (8) :: STT (I_, J_, MAX_L_, NTR_)  
      real (8) :: TREF (51, 18, 12), OREF (51, 18, 12)  
      character (len=10) :: TCNAME (NTR_)  
!
      integer, pointer            :: CLDWNDX(:,:,:)
      integer, pointer            :: CLDINDX(:,:,:)
      real*8,  pointer            :: CLDW(:,:,:)
      real*8,  pointer            :: CLDI(:,:,:)
      real*8,  pointer            :: optDust(:,:)
      real*8,  pointer            :: optAer (:,:)
      real*8,  pointer            :: AER    (:,:)
      integer, pointer            :: ADX    (:,:)
      real*8,  pointer            :: AER1P_aero(:,:)
      real*8,  pointer            :: AER1P_dust(:,:)

      common /metdatAerosol/ optDust, optAer, AER, ADX, &
                             AER1P_aero, AER1P_dust, AERindex
!
      common / metdat / P, T, STT, XGRD, YGRD, XDGRD, YDGRD, ETAA, ETAB, &
               AREAXY, TREF, OREF, PMEAN, SA, MONTH, NSLAT, NSLON, TCNAME, &
               CLDW, CLDI, CLDWNDX, CLDINDX

      common /jvdatIJ/ TJ, DM, DO3, ZH, PJ, AER1P, AER2P, &
     &                  AER1N, AER2N

