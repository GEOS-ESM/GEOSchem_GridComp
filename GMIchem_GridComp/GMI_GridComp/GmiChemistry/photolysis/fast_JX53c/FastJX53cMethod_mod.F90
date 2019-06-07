!-------------------------------------------------------------------------
!  NASA/GFSC, SIVO, Code 610.3
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: FastJX53cMethod_mod
!
! !INTERFACE:
      module FastJX53cMethod_mod
!
! !USES:
      use FastJX53cMetVars_mod, only : MONTH, PMEAN
!
      implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
!
      private
      public  :: InitializeFastJX53c
      public  :: RunFastJX53c
      public  :: FinalizeFastJX53c
!
! !DESCRIPTION:
!  Routines to Initialize, Run and Finalize the package fastJ.
!
! !AUTHOR:
! Jules Kouatchou
!
! !HISTORY:
!
!EOP
!-----------------------------------------------------------------------
      CONTAINS
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: InitializeFastJX53c
!
! !INTERFACE:
!
      subroutine InitializeFastJX53c &
                 (k1, k2, chem_mask_khi, num_qjs, &
                   CrossSection_InfileName  , ScatteringData_InfileName,   &
                   PhotRateLabels_InfileName, T_O3_Climatology_InfileName)
!
! !USES:
      use FastJX53cCTMparameters_mod, only : SetFastjCTMdims
      use FastJX53cMetVars_mod      , only : InitializationFastjMetVars
      use FastJX53cCTMspecific_mod  , only : INPHOT, SET_CLD0, SET_AER0, SET_ATM0
      use FastJX53cJvaluesVars_mod  , only : InitializationFastjJvaluesVars

      implicit none
!
! !INPUT PARAMETERS:
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
!
! !DESCRIPTION:
!  Carry out the initialization of the package FastJX53c.
!
! !LOCAL VARIABLES:
      integer :: LongDim     ! Longitude dim of CTM grid
      integer :: LatDim      ! latidude  dim of CTM grid
      integer :: VertDim     ! altitude  dim of CTM grid
      integer :: NumJval     ! num. photolysis reactions 
      integer :: VertDimJval ! vertical(levels) dim for J-values
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

      longDim     = 1
      latDim      = 1
      vertDim     = k2-k1+1
      numJval     = num_qjs
      vertDimJval = chem_mask_khi

      ! Set values of some CTM related dimensions
      call SetFastjCTMdims (LongDim, Latdim, VertDim, VertDimJval, NumJval)

      ! Initialize MetFields related variables.
      call InitializationFastjMetVars ()
 
      ! Initialize Jvalues related variables.
      call InitializationFastjJvaluesVars()

      ! Read in the input files
      call INPHOT (CrossSection_InfileName,     &
     &             ScatteringData_InfileName,   &
     &             PhotRateLabels_InfileName,   &
     &             T_O3_Climatology_InfileName)

!      call SET_CLD0(Temp, CloudOptDepth, IndexCloud, SurfAlbedo)

      ! Set up aerosol
      call SET_AER0

      return

      end subroutine InitializeFastJX53c
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  RunFastJX53c
!
! !INTERFACE:
!
      subroutine RunFastJX53c(k1, k2, IDAY, timeSec, iMONTH, sza_ij, &
                          press3c_ij, pctm_ij, kel_ij, optdepth_ij, &
                          surf_alb_ij, PhotRate, ozone_ij)
!
! !USE:
      use FastJX53cMetVars_mod         , only : PJ, PMEAN, MONTH, SA, T, NCLDX, ODCLD
      use FastJX53cMetVars_mod         , only : XDGRD, YDGRD, XGRD, YGRD, DO3, doModelClima
      use FastJX53cCTM_FastjLinking_mod, only : PHOTOJ
      use FastJX53cCTMspecific_mod     , only : SET_ATM0, SET_ATM, SET_AER, SET_CLD
      use FastJX53cCTMparameters_mod   , only : I_, J_, L_, JVL_, JVN_, L1_

      implicit none
      
! !INPUT PARAMETERS:
      integer, intent(in)  :: k1, k2
      real*8 , intent(in)  :: timeSec            ! Time of day in model (s)
      integer, intent(in)  :: IDAY               ! Current day
      integer, intent(in)  :: imonth             ! Current Month
      real*8,  intent(in)  :: kel_ij(k1:k2)      ! temperature at box centers (degK) [t]
      !real*8,  intent(in)  :: londeg_i           ! longitude (midpoint, deg) [xdgrd]
      !real*8,  intent(in)  :: latdeg_j           ! latitude  (midpoint, deg) [ydgrd]
      real*8,  intent(in)  :: sza_ij             ! solar zenith angle
      real*8,  intent(in)  :: press3c_ij(k1:k2)
      real*8,  intent(in)  :: pctm_ij            ! surface pressure (mb) [p]
      real*8,  intent(in)  :: optdepth_ij(k1:k2) ! optical depth in box   (unitless) [od]
      real*8,  intent(in)  :: surf_alb_ij        ! surface albedo     (fraction 0-1) [sa]
      real*8, optional, intent(in) :: ozone_ij(k1:k2)      ! mixing ratio for ozone
!
! !OUPUT PARAMETERS:
      real*8 , intent(out) :: PhotRate(JVL_,JVN_)  ! array of J's indexed to CTM chemistry
!
! !LOCAL VARIABLES:
      integer :: ILNG, JLAT
      real*8  :: GMTAU       ! Time of day (hr)
      real*8  :: PI
      real*8  :: FREFL, U0            ! fraction of flux(energy-wtd) reflected
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
      
      PI                   = 4.0d0*atan(1.0d0)
      ILNG                 = I_
      JLAT                 = J_

      MONTH                = imonth                 ! current month
      GMTAU                = timeSec / 3600.0d0     ! time (hr)

      !XDGRD(I_)            = londeg_i
      !YDGRD(J_)            = latdeg_j
      !XGRD(I_)             = XDGRD(I_) *PI/180.d0
      !YGRD(J_)             = YDGRD(J_) *PI/180.d0

      PJ(1)                = pctm_ij
      PJ(2:L1_)            = press3c_ij(k1:k2)      ! pressure at the edge
      PJ(L1_+1)            = 0.0d0

      PMEAN(I_, J_)        = pctm_ij                ! surface pressure
      T(I_, J_, 1:L_)      = kel_ij(k1:k2)          ! temperature
      ODCLD(I_, J_, 1:L_)  = optdepth_ij(k1:k2)     ! cloud optical depth
      ODCLD(I_, J_, L1_)   = 0.0
      SA(I_, J_)           = surf_alb_ij            ! surface albedo
      NCLDX(I_, J_, 1:L1_) = 10

!---reset the atmosphere, aerosols, clouds for the time step (now = dummy)
!-----------------------------------------------------------------------

      ! Set up atmospheric variables at the current time if using
      ! fastJ climatology.
      if (.not. present(ozone_ij)) then
        call SET_ATM0
      endif

      ! Set up atmospheric variables at the current time if using CTM
      ! climatology
      if (present(ozone_ij)) then
         doModelClima = .true.
         DO3(I_, J_, 1:L_)  = ozone_ij(k1:k2)
         DO3(I_, J_, L1_ )  = ozone_ij(k2)
         call SET_ATM(GMTAU)
      endif

      ! Set up aerosols at the current time
      call SET_AER(GMTAU)

      ! Set up cloud at the current time
      call SET_CLD(GMTAU)

      ! Calculate the photolysis rate
      call PHOTOJ(GMTAU,IDAY,ILNG,JLAT, sza_ij, U0, FREFL, PhotRate)

      return

      end subroutine RunFastJX53c
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  FinalizeFastJX53c
!
! !INTERFACE:
!
      subroutine FinalizeFastJX53c()
!
      implicit none
!
! !DESCRIPTION:
!  Finalize the package FastJX53c by deallocating (if necessary) variables.
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

      return

      end subroutine FinalizeFastJX53c
!EOC
!-----------------------------------------------------------------------

 end module FastJX53cMethod_mod
