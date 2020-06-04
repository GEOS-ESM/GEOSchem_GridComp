!---------------------------------------------------------------------------
!BOC
      module CloudJ_mod
        USE FJX_CMN_MOD
        USE FJX_SUB_MOD
        USE FJX_INIT_MOD
        USE CLD_SUB_MOD, ONLY : CLOUD_JX
        USE OSA_SUB_MOD
!.sds  added old GMI aerosol parameters
        USE GMI_CMN_MOD
!.sds  added MAPL physical constants being used
        use MAPL_ConstantsMod, only: MAPL_PI, MAPL_PI_R8
!.sds  added LDRADIUS[3,4] to GEOS_UtilitiesMod.F90
!       use GEOS_UtilsMod, only: LDRADIUS3, LDRADIUS4
!
        public :: controlFastJX74
        public :: initializeFastJX74
        public :: RD_gmi
!
!---------------------------------------------------------------------------
      contains
!---------------------------------------------------------------------------
!
!>>>>>>>>  Cloud-J version 8.0   new lower albedo can be angle/wavelength dependent
!
      subroutine controlFastJX74 (k1, k2, chem_mask_khi, lat_ij, num_qjs, &
     &                  month_gmi, jday, time_sec, do_clear_sky, cldflag, gridBoxHeight_ij, &
     &                  SZA_ij, cloudfrac_ij, qi_ij, ql_ij, ri_ij, rl_ij, &
!    &                  cnv_frc_ij, frland_ij, &
     &                  press3e_ij, pctm_ij, kel_ij,                   &
     &                  surf_alb_ij, qjgmi_ij, relHumidity_ij,         &
     &                  overheadO3col_ij, ODAER_ij, ODMDUST_ij, ODcAER_ij, HYGRO_ij,    &
     &                  do_AerDust_Calc, AerDust_Effect_opt, cldOD_ij, eradius_ij, tArea_ij,      &
     &                  fjx_solar_cycle_param, CH4_ij, H2O_ij, ozone_ij)
!
! USES:
!
      implicit none
!
# include "gmi_AerDust_const.h"
# include "gmi_phys_constants.h"
# include "gmi_time_constants.h"
# include "setkin_par.h"
!---------------key params in/out of CLOUD_J-------------------------
      integer, intent(in) :: k1, k2, num_qjs, chem_mask_khi
      integer, intent(in) :: jday, month_gmi, cldflag
      integer, intent(in) :: AerDust_Effect_opt
      logical, intent(in) :: do_clear_sky, do_AerDust_Calc
      real*8, intent(in) :: lat_ij  !... grid point latitude (degrees)
      real*8, intent(in) :: time_sec, SZA_ij, pctm_ij, surf_alb_ij
!     real*8, intent(in) :: cnv_frc_ij, frland_ij
      real*8, intent(in), dimension(k1-1:k2) :: press3e_ij
      real*8, intent(in), dimension(k1:k2)   :: kel_ij, relHumidity_ij
      real*8, intent(in), dimension(k1:k2)   :: gridBoxHeight_ij
      real*8, intent(in), dimension(k1:k2)   :: cloudfrac_ij, qi_ij, ql_ij, ri_ij, rl_ij
      real*8, intent(in), dimension(W_)      :: fjx_solar_cycle_param
      real*8, intent(in), dimension(k1:k2)   :: CH4_ij, H2O_ij(k1:k2)  ! mixing ratio for CH4 & H2O
      real*8, intent(in), dimension(k1:k2), optional :: ozone_ij      ! mixing ratio for ozone
      real*8, intent(out), dimension(k1:k2)                    :: overheadO3col_ij
      real*8, intent(out), dimension(k1:chem_mask_khi,num_qjs) :: qjgmi_ij
      real*8, intent(inout), dimension(k1:k2,NSADaer*nrh_b)    :: ODAER_ij
      real*8, intent(inout), dimension(k1:k2,NSADaer)          :: HYGRO_ij
      real*8, intent(inout), dimension(k1:k2,2)                :: ODcAER_ij
      real*8, intent(inout), dimension(k1:k2,NSADdust)         :: ODMDUST_ij
      real*8, intent(inout), dimension(k1:k2,NSADdust+NSADaer) :: ERADIUS_ij, TAREA_ij
      real*8, intent(inout), dimension(k1:k2)                  :: cldOD_ij
!
      logical, parameter           :: LPRTJ=.false.
      integer                      :: TCLDFLAG, NRANDO, IRAN, LNRG
      integer                      :: NICA, JCOUNT
      real*8                       :: U0, SZA, CLDCOR
      integer, dimension(1:k2-k1+2)    :: CLDIW
      integer, dimension(1:k2-k1+2,AN_) :: IDXAER
      real*8, dimension(1:k2-k1+2,AN_) :: AERSP, aerOD_out
      real*8, dimension(W_)        :: SOLF
      real*8, dimension(0:k2-k1+2) :: PPP,ZZZ
      real*8, dimension(1:k2-k1+2) :: TTT, DDD, RRR, OOO
      real*8, dimension(1:k2-k1+2) :: LWP, IWP, REFFL, REFFI, cldOD_out
      real*8, dimension(1:k2-k1+1,num_qjs) :: VALJXX
      real*8, dimension(5,W_+W_r)  :: RFL
!
!-------------local use-----------------------
      integer :: MONTH
      integer :: I, J, K, L, LL, N
      integer :: LTOP  !, NJXX
      integer :: kall, num_aer, loc_naer, gmi_naer
      real :: ANGLES(5), OWAVEL, OWIND, OCHLR, OSA_dir(5), factor
      real :: PDEL, ZDEL, ICWC, F1, WLC4
!     real :: PMID, NN3, NN4, NN_LAND, NN_OCEAN, TTT4, WIC4
      real*8  :: GMTAU, CO2_IJ, baseSOLF, rho
      real*8, dimension(L_) :: WLC, WIC, CLF
      real, dimension(L_) :: CLDFRW, CLDIWCW, CLDLWCW, CLDREFFI, CLDREFFL
      logical :: oldgmi_aero=.true.
!
!. tropSO4, BC, OC, SeaSalt(accum), SeaSalt(coarse)
      integer :: IR, IRH, ioffd, kdry
      integer, parameter :: idxAtype(NSADaer)=(/-1, -14, -22, -4, -5/), NRH_GMI=(7)
!. old tropSO4, BC, OC, SeaSalt(accum), SeaSalt(coarse) densities
      real, dimension(NSADdust), parameter :: dMSDENS=(/2500.0, 2500.0, 2500.0, 2500.0, 2650.0, 2650.0, 2650.0/)
      real, dimension(NSADaer), parameter :: wMSDENS=(/1700.0, 1000.0, 1800.0, 2200.0, 2200.0/)
!      INTEGER, parameter :: IND(NSADaer) = (/22, 29, 36, 43, 50/)
      REAL*8 :: FRAC, REFF, scaleR, scaleQ
!
!
!
! Cldflag can be from 1 to 8.  The recommended one is 7.   Here are the references for the choice of cldflag =7, 
!
! M.J. Prather (2015) Photolysis rates in correlated overlapping cloud fields: Cloud-J 7.3c,
!  Geosci. Model Dev., 8, 2587-2595, doi:10.5194/gmd-8-2587-2015. 
!
! Neu, J.L., M.J. Prather, J.E. Penner (2007) Global atmospheric chemistry:  integrating over
!  fractional cloud cover, J. Geophys. Res., 112, D11306. 
!
! Cldflag = 1 
!  completely ignore cloud conditions and reset LWP =IWP=0, clear sky case 
!  (aerosols are not excluded if they exist) LWP=liquid water path, IWP=ice water path.

! Cldflag =2 (average cloud, or overcast cloud --take the input LWP and IWP as the in-cloud
!  condition, reset them by multiplying them with the cloud fraction, clf, and then reset clf to
!  be one. It's like spreading the cloud to the whole grid to appear overcast) 
!
! Cldflag=3 similar to cldflag=2 but multiply by clf^{3/2), instead of clf^{1}
!
! Cldflag =1 to 3  consider only one atmosphere per (i,j) grid point. One single atmosphere is
!  passed to do the radiative transfer calculation (call photo_jx) --  either clear sky (1) or
!  overcast cloud (2, 3) 

! Cldflag =4, 5, 6, 7, 8, however, consider multiple possible atmospheres per (i,j) grid point
!  by using the information of cloud fractions.  Lets imagine the horizontal grid is large (say 1
!  degree by 1 degree), and you look up at surface, the clouds are typically broken (not filled),
!  and some are aligned vertically and some aren't.  Since we don't have the exact information
!  how the clouds are aligned vertically,  we use the vertical distribution of cloud fractions to
!  come up with possible combinations of so called ICAs (independent column atmospheres), at each
!  grid point.    

! The way to decide how may ICAs are generated is decided by (a) first organizing them into
!  groups from all tropospheric layers of a column, and (b) between these sorted groups from (a),
!  decide how they are overlapped or not by all possible combinations.  

! Grouping (a) can be based on the assumption of the maximum overlap between adjacent levels
!  of clouds. The break between different groups are the existences of zero clouds. For example,
!  there might be cloud free layers in the middle of the column between low-level stratus and
!  high level cirrus. The most known scheme for this is called the maximum-random overlap.  
!  Michael came up with another algorithm (correlated cloud overlap) by using the observed
!  established correlations between vertically adjacent clouds as the basis for grouping (see his
!  2015 paper above).  
!
! Once you obtain N-ICAs from (a) and (b), it can be very expensive if you call the
!  radiation-transfer subroutine (call photo_jx)  N-times and then average the N results for the
!  final radiation budget of the grid point (although this is precisely the choice for
!  cldflag=8).  

! There are cheaper options computationally for dealing with that.  Basically N- ICAs can be
!  sorted out into four quadratures according to the total column cloud optical depth of each ICA
!  (tau = 0-0.5, 0.5-4, 4-30, >30). Within each quadrature, we can choose one representative ICA
!  (a) by averaging among them by the weight function for each ICA (b) picking the mid-point of
!  its CDF distribution.  For (a), cldflag =7, for (b) cldflag=6.  Please see Neu et al. (2007)
!  for details.  After obtaining the representative atmosphere of each quadrature, one only has
!  to call photo_jx four times and then average them for the radiation budgets.  
!
! You can also just pick a subset of ICAs randomly by drawing NRANDO, = 6  ICAs (for example)
!  from the cdf distributions of ICAs, and this option is cldflag=5. Although I'm not sure what
!  RAN4 is (a fortran intrinsic function?)  
!
!
! from Michael:
!
! According to what I have in the ATom J-paper, you are running:
!   Fast-J v6.5, liquid cloud C1 (6 um) and ice cloud hexagonal (50 um), Briegleb
!   averaging 
!
! Thus if you pick CLDFLAG=3, it should come very close to your old cloud averaging. You should
!  start with this one. 
! There may be differences in the cloud opacities, they are determined primarily by your Reff,
!  which says 6 microns above (that is on the low side = higher opacities). 
!
! As Juno says, CLDFLAG = 7 is recommended, 6 is OK, but has some larger noise, 8 is too
!  costly. 
!
! For the cloud overlap algorithm that generates the ICAs, you should probably use
!   LNRG = 6 (MAX-COR), & CLDCOR = 0.33 (de-correlation factor)
!   ATM0 runs from 0(flat) 1(spherical) 2(+refraction) 3(+geometric expansion)
!  I am not sure the version you have does 3, so I would stick with 1 = ATM0 
!
! For J-values:
!  I would stick with the 2 classic ones we used in the ATom J-paper (attached)- J-O1D and J-NO2
!  They cover a range between being sensitive to O3 to just being Rayleigh and cloud driven.
!  You ran these for me for the paper. And we have a history of comparing them (PHOTOCOMP, ..
!
! So first run with CLDFLAG=3 and ATM0=1, that should be closest to what the earlier Fst-J
!  codes might give. 
!
! Michael
!
!
!      call CLOUD_JX (U0, SZA, RFL, SOLF, LPRTJ, PPP, ZZZ, TTT,  &
!             DDD, RRR, OOO, LWP, IWP, REFFL, REFFI, CLF, CLDCOR, CLDIW,  &
!             AERSP, IDXAER, L_+1, AN_, VALJXX, num_qjs,  &
!             TCLDFLAG, NRANDO, IRAN, LNRG, NICA, JCOUNT)
!
! PPP     = gridbox edge pressures (hPa)
! TTT     = gridbox temperatures (K)
! RRR     = Relative humidity (%)
! ZZZ     = gridbox bottom height [edge] (cm)
! U0      = cos(SZA)
! SZA     = solar zenith angle (degrees) - PASS IN
! RFL     = Lambertian albedo of surface for angles 1:4 & U0 (#5)
! SOLF    = solar flux
! LPRTJ   = print diagnostic - disable
! PPP     = gridbox edge pressures (hPa)
! ZZZ     = gridbox edge heights (cm)
! TTT     = gridbox temperatures (K)
! DDD     = column density - delp*MASFAC (molecules/cm3?)
! RRR     = Relative humidity (%)
! OOO     = gridbox ozone (vmr, converted below to molecules/cm3)
! LWP     = liquid water path (in layer, in cloud) in g/m**2
! IWP     = ice water path (in layer, in cloud) in g/m**2
! REFFL   = effective radius of liquid cloud (microns)
! REFFI   = effective radius of ice cloud (microns)
! CLF     = Cloud fraction
! CLDCOR  = Cloud correction for CLDFLAG=5?
! CLDIW   = integer flag: 1 = water cld, 2 = ice cloud, 3 = both
! AERSP   = aerosol species array - path (g/m2) of aerosol
! IDXAER  = aerosol index for assigning optical characteristics
! L1_     = 
! AN_     = number of aerosol types
! VALJXX  = 
! JVN_    = 
! CLDFLAG = type of cloud handling
! NRANDO  = parameter for CLDFLAG = 5
! IRAN    = parameter for CLDFLAG = 5
! LNRG    = cloud overlap param
! NICA    = output - No. ICAs
! JCOUNT  = output
!
      MONTH = month_gmi
      GMTAU = time_sec / SECPHR
!
!... map in CTM date to arrays to input to CloudJ, add extra level on top!
!... gridbox edge pressures (hPa) - top to 0.0!!
      PPP(0:k2-k1+1) = press3e_ij(k1-1:k2)
      PPP(k2-k1+2) = 0.0d0
!-------calculate effective altitude of each CTM level edge (cm)
      ZZZ(0) = 0.0
      do L=1,k2-k1+1
        ZZZ(L) = ZZZ(L-1) + gridBoxHeight_ij(L) * 100.d0
      enddo
!.???.. what to do?
      ZZZ(k2+1) = ZZZ(k2) + ZZZ(k2)
!
!... gridbox temperatures (K)
      TTT(1:k2-k1+1) = kel_ij(k1:k2)
      TTT(k2-k1+2) = kel_ij(k2)   !/2.0d0
!... Relative humidity (%)
      RRR(1:k2-k1+1) = relHumidity_ij(k1:k2)
      RRR(k2-k1+2) = relHumidity_ij(k2)   !/2.0d0
!... solarZenithAngle
      SZA = SZA_ij
!
!... Set surface albedo, ignoring wavelength dependence
      RFL(:,:) = max (0.d0, min (1.d0, surf_alb_ij) )  
!
!---set T, D & Z
!... molecules/cm2
      do L=1,k2-k1+2
        DDD(L) = (PPP(L-1)-PPP(L)) * MASFAC
      enddo
! 
!---set O3
!... gridbox ozone (vmr?)
      if (present(ozone_ij)) then
        OOO(1:k2-k1+1) = ozone_ij(k1:k2) * DDD(k1:k2)
        OOO(k2-k1+2) = ozone_ij(k2) * DDD(k2-k1+2)
      else
!... read in Michael's T and O3 climotology
        call ACLIM_FJX (lat_ij, MONTH, PPP, TTT, OOO, L_+1)
      endif
!
!... set to 18 to turn off LW fluxes
      NSBIN = 18
!---sets longwave important process (CO2)
      if (LRRTMG) then
!... set to allow solar LW fluxes
        NSBIN = 27
        CO2_IJ = 350.e-6
        colo3(:)  = OOO(:)
        colco2(:) = CO2_IJ * DDD(:)               !CO2 
        colch4(1:k2-k1+1) = CH4_ij(k1:k2) * 1.0d9 !CH4 vmr in ppb
        colch4(k2+1) = CH4_ij(k2) * 1.0d9 
        colo2(:)  = MXRO2 * DDD(:)
        colh2o(1:k2-k1+1) = H2O_IJ(k1:k2)         !H2O in kg/kg
        colh2o(k2+1) = H2O_IJ(k2)
        coldry(:) = DDD(:) - colh2o(:)
      endif
!
!!! set up aerosols
!... make sure surface area is initialized
      ERADIUS_ij(:,:) = 0.0D0
      TAREA_ij(:,:) = 0.0D0
      HYGRO_ij(:,:) = 0.0D0
      cldOD_out(:) = 0.0D0
      aerOD_out(:,:) = 0.0D0
      num_aer = 0

!.sds... NEED TO SET UP AEROSOLS YET!!!!!!!!!!!!!1
!---U Michigan aerosol data sets, this generate fast-JX data formats.
!---Approximates the Legendre expansion(L) of the scattering phase fn as (2*L+1)*g**L
!---UMAER(I,J,K,L):
!   I=1:3 = [SSAbldeo, g, k-ext(m2/g)]
!   J=1:5 = [200, 300, 400, (550,) 600 , 1000 nm]
!   K=1:21= [0, 5, 10, 15, ..., 90, 95, 99 %RelHum]
!   L=1:33= UM aerosol types [SULF, SS-1,-2,-3,-4, DD-1,-2,-3,-4, FF00(0%BC),
!                      FF02, ...FF14(14%BC), BB00, BB02, ...BB30(30%BC)]
!... aerosols zeroed out now
! Old number:  4,5, 15,16,17,18,19,20,21
! New number:  1,2, 11,12,13,14,15,16,17
!... units: path (g/m2) of aerosol/cloud
      AERSP(:,:)  = 0.d0
      IDXAER(:,:) = 0
!... turn on radiative effects of tropo aerosols
      if( do_AerDust_Calc .and. AerDust_Effect_opt.lt.3 ) then
!
!... dust map to UCI tables and convert to g/m^2
!.... GOCART Bins    Reff                  UMa idx 
!     dustreff(1) = 0.73E-6 ! Clay         -7
!     dustreff(2) = 1.4E-6  ! Small silt   -8
!     dustreff(3) = 2.4E-6  ! Small silt   -8
!     dustreff(4) = 4.5E-6  ! Small silt   -9
!     dustreff(5) = 8.0E-6  ! Small silt   -9
  ! - du001 - du005 for the following five dust bins (see DU_GridComp.rc in
  !   GOCART):
  !   radius_lower: 0.1 1.0 1.8 3.0 6.0
  !   radius_upper: 1.0 1.8 3.0 6.0 10.0
!. scat-UMa tables - all wl=(200/300/400/550/600/1000nm) RH=(0%-5%-10%-...-90%-95%-99%)
!    (-6)  DD-1     dust_aop.r0.05to0.63.dat
!    (-7)  DD-2     dust_aop.r0.63to1.26.dat
!    (-8)  DD-3     dust_aop.r1.26to2.5.dat
!    (-9)  DD-4     dust_aop.r2.5to10.dat
!
! for subroutine OPTICA
!01 S-Bkg   LOGN:r=.090 s=.600 n=1.514/.../1.435  reff=0.221___G=.0523_rho=1.630
!02 S-Vol   LOGN:r=.080 s=.800 n=1.514/.../1.435  reff=0.386___G=.0721_rho=1.630
!03 UT-sulfate LOGN:r=0.05 s=.693 n=1.44          reff=0.166___G=.0205_rho=1.769
!04 UT-sulfate LOGN:r=0.05 s=.693 n=1.46          reff=0.166___G=.0205_rho=1.769
!05 UT-sulfatM LOGN:r=.050 s=.642 n=1.53          reff=0.140___G=.0179_rho=1.769
!06 UM-BC1     LOGN:r=.050 s=.642 n=1.80+0.50i    reff=0.140___G=.0179_rho=1.500
!07 UM-BC2     LOGN:r=.080 s=.501 n=1.80+0.50i    reff=0.150___G=.0332_rho=1.500
!08 UM-BB08 (%BC)LOGN:r=.080 s=.500 n=1.552+0.04i reff=0.149___G=.0331_rho=1.230
!09 UM-FF04 (%BC)LOGN:r=.050 s=.642 n=1.541+0.02i reff=0.140___G=.0179_rho=1.212
!10 UM-FF10 (%BC)LOGN:r=.050 s=.642 n=1.557+0.05i reff=0.140___G=.0179_rho=1.230
!11 MDust.15  (R.V. Martin generated phase fns)   reff=0.150___G=1.000_rho=2.600
!12 MDust.25  (R.V. Martin generated phase fns)   reff=0.250___G=1.000_rho=2.600
!13 MDust.40  (R.V. Martin generated phase fns)   reff=0.400___G=1.000_rho=2.600
!14 MDust.80  (R.V. Martin generated phase fns)   reff=0.800___G=1.000_rho=2.600
!15 MDust1.5  (R.V. Martin generated phase fns)   reff=1.500___G=1.000_rho=2.600
!16 MDust2.5  (R.V. Martin generated phase fns)   reff=2.500___G=1.000_rho=2.600
!17 MDust4.0  (R.V. Martin generated phase fns)   reff=4.000___G=1.000_rho=2.600
!
!... remapped 7 bins = [0.150(11) 0.250(12) 0.400(13)  0.800(14) 1.500(15) 2.500(16) 4.000(17)]
!
! new_type  dAersl 
!    6      1: BC hydrophobic 29    MSDENS(1) = 1000.0   raa=0.039
!    8?     2: OC hydrophobic 36    MSDENS(1) = 1800.0   raa=0.070
! wAersl
!    -1     1: SO4            22+   MSDENS(1) = 1700.0   raa=0.159 0.217 0.241 0.260 0.297 0.347 0.498
!    -14    2: BC             29+   MSDENS(2) = 1000.0   raa=0.039 0.039 0.039 0.047 0.055 0.059 0.075
!    -22    3: OC             36+   MSDENS(3) = 1800.0   raa=0.070 0.087 0.095 0.102 0.116 0.133 0.177
!    -4     4: SS (accum)     43+   MSDENS(4) = 2200.0   raa=0.732 1.177 1.324 1.457 1.740 2.119 3.484
!    -5     5: SS (coarse)    50+   MSDENS(5) = 2200.0   raa=5.674 9.024 10.107 10.879 12.372 14.057 18.159
!
!... tropospheric SO4
!. wAersl(1) = "SO4"+"SO4v"  == #22
!.old tables: 
!   22 S00(rvm) Trop sulfate at RH=00 (n@400=1.44 log-norm: r=.05um/sigma=2.0)
!     w(nm)    Q    r-eff  ss-alb  pi(0) pi(1) pi(2) pi(3) pi(4) pi(5) pi(6) pi(7)
!     300  2.1411  0.159  1.0000  1.000 2.107 2.570 2.325 2.050 1.662 1.353 1.095
!   23 S50(rvm) Trop sulfate at RH=50 (n@400=1.38 log-norm: r=.06um/sigma=2.0)
!   24 S70(rvm) Trop sulfate at RH=70 (n@400=1.36 log-norm: r=.07um/sigma=2.0)
!   25 S80(rvm) Trop sulfate at RH=80 (n@400=1.36 log-norm: r=.07um/sigma=2.0)
!   26 S90(rvm) Trop sulfate at RH=90 (n@400=1.35 log-norm: r=.08um/sigma=2.0)
!   27 S95(rvm) Trop sulfate at RH=95 (n@400=1.35 log-norm: r=.10um/sigma=2.0)
!   28 S99(rvm) Trop sulfate at RH=99 (n@400=1.34 log-norm: r=.14um/sigma=2.0)
!... 
!. New tables
!   (scat-aer)
!    (03)    UT-sulf1| 0.166 1.769| upper-trop sulf1 logN:r=0.05:s=.693  n=1.44
!        300 2.1375 1.0000 2.165 2.685 2.493 2.224 1.832 1.498 1.224
!    (04)    UT-sulf2| 0.166 1.769| upper-trop sulf2 logN:r=0.05:s=.693  n=1.46
!        300 2.2061 1.0000 2.129 2.621 2.399 2.131 1.741 1.424 1.15
!    (05)    UT-sulfM| 0.140 1.769| upper-trop sulf (UMich) logN:r=.050:s=.642 n=1.53
!        300 2.2755 1.0000 2.002 2.337 1.993 1.678 1.278 0.999 0.752
!.  scat-UMa tables - all wl=(200/300/400/550/600/1000nm) RH=(0%-5%-10%-...-90%-95%-99%)
!    (-1)  SULF     sulfate_all.dat
!       (300, 0%) 1.00000  0.67176  6.88417 !(SSA, g, k-ext)
! ODAER_ij
!
!.old tables:
!... BC
!. dAersl(1) = "BCphobic"  == #29
!. wAersl(2) = "BCphilic"  == #29
!   29 BC00(rvm) Black C, RH=00 (n@400=1.75-.46i log-norm: r=.01um/sigma=2.0)
!     300  1.0643  0.039  0.3128  1.000 1.361 1.201 0.710 0.410 0.231 0.135 0.080
!   30 BC50(rvm) Black C, RH=50 (n@400=1.75-.46i log-norm: r=.01um/sigma=2.0)
!   31 BC70(rvm) Black C, RH=70 (n@400=1.75-.46i log-norm: r=.01um/sigma=2.0)
!   32 BC80(rvm) Black C, RH=80 (n@400=1.57-.27i log-norm: r=.01um/sigma=2.0)
!   32 BC80(rvm) Black C, RH=80 (n@400=1.57-.27i log-norm: r=.01um/sigma=2.0)
!   33 BC90(rvm) Black C, RH=90 (n@400=1.48-.17i log-norm: r=.02um/sigma=2.0)
!   34 BC95(rvm) Black C, RH=95 (n@400=1.45-.14i log-norm: r=.02um/sigma=2.0)
!   35 BC99(rvm) Black C, RH=99 (n@400=1.39-.07i log-norm: r=.02um/sigma=2.0)
!... OC
!. dAersl(2) = "OCphobic"  == #36
!. wAersl(3) = "OCphilic"  == #36
!   36 OC00(rvm) Organic C, RH=00 (n@400=1.53-.005i log-norm: r=.02um/sigma=2.0)
!     300  1.0162  0.070  0.9530  1.000 1.871 1.938 1.467 1.056 0.706 0.480 0.316
!   37 OC50(rvm) Organic C, RH=50 (n@400=1.44-.003i log-norm: r=.03um/sigma=2.0)
!   38 OC70(rvm) Organic C, RH=70 (n@400=1.42-.002i log-norm: r=.03um/sigma=2.0)
!   39 OC80(rvm) Organic C, RH=80 (n@400=1.40-.002i log-norm: r=.03um/sigma=2.0)
!   40 OC90(rvm) Organic C, RH=90 (n@400=1.38-.001i log-norm: r=.03um/sigma=2.0)
!   41 OC95(rvm) Organic C, RH=95 (n@400=1.37-.001i log-norm: r=.04um/sigma=2.0)
!   42 OC99(rvm) Organic C, RH=99 (n@400=1.350.000i log-norm: r=.05um/sigma=2.0)
!
!. New tables
!   (scat-aer)
!    (06)    UM-BC1  | 0.140 1.500| UMich w/Mie code logN:r=.050:s=.642  n=1.80+0.50i
!    (07)    UM-BC2  | 0.140 1.500| UMich w/Mie code logN:r=.080:s=.501  n=1.80+0.50i
!    (08)    UM-BB08C| 0.149 1.230| UMich w/Mie code 8%BC 0%RH logN:r=.080:s=.500 n=1.552+0.04i
!    (09)    UM-FF04C| 0.140 1.212| UMich w/Mie code 4%BC 0%RH logN:r=.050:s=.642 n=1.541+0.02i
!    (10)    UM-FF10C| 0.140 1.230| UMich w/Mie code 10%BC 0%RH logN:r=.050 s=.642 n=1.557+0.05i
!.  scat-UMa tables - all wl=(200/300/400/550/600/1000nm) RH=(0%-5%-10%-...-90%-95%-99%)
!    (-10) FF00     ffi0c.RH99.r000.dat
!    (-11) FF02     ffi0c.RH99.r002.dat
!    (-12) FF04     ffi0c.RH99.r004.dat
!    (-13) FF06     ffi0c.RH99.r006.dat
!    (-14) FF08     ffi0c.RH99.r008.dat
!    (-15) FF10     ffi0c.RH99.r010.dat
!    (-16) FF12     ffi0c.RH99.r012.dat
!    (-17) FF14     ffi0c.RH99.r014.dat
!    (-18) BB00     bio0c.RH99.r000.dat
!    (-19) BB02     bio0c.RH99.r002.dat
!    (-20) BB04     bio0c.RH99.r004.dat
!    (-21) BB06     bio0c.RH99.r006.dat
!    (-22) BB08     bio0c.RH99.r008.dat
!    (-23) BB10     bio0c.RH99.r010.dat
!    (-24) BB12     bio0c.RH99.r012.dat
!    (-25) BB14     bio0c.RH99.r014.dat
!    (-26) BB16     bio0c.RH99.r016.dat
!    (-27) BB20     bio0c.RH99.r020.dat
!    (-28) BB22     bio0c.RH99.r022.dat
!    (-29) BB24     bio0c.RH99.r024.dat
!    (-30) BB26     bio0c.RH99.r026.dat
!    (-31) BB28     bio0c.RH99.r028.dat
!    (-32) BB30     bio0c.RH99.r030.dat
!
!... SS
!. current method:
! - ss001-ss005 for the following five sea salt aerosol bins (see SS_GridComp.rc in GOCART):
!    radius_lower: 0.03 0.1 0.5 1.5 5.0
!    radius_upper: 0.1  0.5 1.5 5.0 10.0
!.... GMI usage:
!...accum wAersl(4) = "SS001"+"SS002"  == #43
!    (43) SSa00(rvm) Sea Salt (accum), RH=00 (n@400=1.50 log-norm: r=.21um/sigma=2.0)
!       300  2.4998  0.732  0.9999  1.000 2.122 3.044 3.083 3.688 3.676 4.110 4.155
!    (44) SSa50(rvm) Sea Salt (accum), RH=50 (n@400=1.38 log-norm: r=.34um/sigma=2.0)
!    (45) SSa70(rvm) Sea Salt (accum), RH=70 (n@400=1.37 log-norm: r=.38um/sigma=2.0)
!    (46) SSa80(rvm) Sea Salt (accum), RH=80 (n@400=1.36 log-norm: r=.42um/sigma=2.0)
!    (47) SSa90(rvm) Sea Salt (accum), RH=90 (n@400=1.35 log-norm: r=.50um/sigma=2.0)
!    (48) SSa95(rvm) Sea Salt (accum), RH=95 (n@400=1.35 log-norm: r=.60um/sigma=2.0)
!    (49) SSa99(rvm) Sea Salt (accum), RH=99 (n@400=1.34 log-norm: r=.99um/sigma=2.0)
!...coarse wAersl(5) = "SS003"+"SS004"+"SS005"  == #50
!    (50) SSc00(rvm) Sea Salt (coarse), RH=00 (n@400=1.50 log-norm: r=1.8um/sigma=2.0)
!       300  2.1044  5.674  0.9995  1.000 2.411 3.610 3.964 4.970 5.369 6.477 7.058
!    (51) SSc50(rvm) Sea Salt (coarse), RH=50 (n@400=1.38 log-norm: r=2.8um/sigma=2.0)
!    (52) SSc70(rvm) Sea Salt (coarse), RH=70 (n@400=1.37 log-norm: r=3.2um/sigma=2.0)
!    (53) SSc80(rvm) Sea Salt (coarse), RH=80 (n@400=1.36 log-norm: r=3.5um/sigma=2.0)
!    (54) SSc90(rvm) Sea Salt (coarse), RH=90 (n@400=1.35 log-norm: r=4.2um/sigma=2.0)
!    (55) SSc95(rvm) Sea Salt (coarse), RH=95 (n@400=1.34 log-norm: r=5.1um/sigma=2.0)
!    (56) SSc99(rvm) Sea Salt (coarse), RH=99 (n@400=1.34 log-norm: r=8.6um/sigma=2.0)
!. New tables all wl=(200/300/400/550/600/1000nm) RH=(0%-5%-10%-...-90%-95%-99%)
!    (-2)  SS-1     ss_aop.r0.05to0.63.dat
!    (-3)  SS-2     ss_aop.r0.63to1.26.dat
!    (-4)  SS-3     ss_aop.r1.26to2.50.dat
!    (-5)  SS-4     ss_aop.r2.50to10.dat
!
        if(oldgmi_aero) then 
          ioffd = 2000   ! aertype > 1000 is use old GMI parameters, aertype > 2000 is hydrophobic
          gmi_naer = 14  ! add 14 for offset in orig GMI table to first dust type
        else
          ioffd = 10     ! use UCI aerosol types
          gmi_naer = 0   ! use UCI aerosol types, no additional offset
        endif
!... dust (hydrophobic)
        do N=1,NSADdust
          num_aer = num_aer+1
          if(oldgmi_aero) then 
            rho = gmiDAA(4,gmi_naer+N)*1000.0
          else
            rho = DAA(10+N)*1.0d+3
          endif
          do L = k1,k2
            LL = L-k1+1
            IDXAER(LL,num_aer) = ioffd+gmi_naer+N
!... capture tropospheric aerosol parameters for diagnostic output
            if(oldgmi_aero) then 
              ERADIUS_ij(l,num_aer) = gmiRAA(4,gmi_naer+N) * 1.0D-4
            else
              ERADIUS_ij(l,num_aer) = RAA(10+N) * 1.0D-4
            endif
!... het reactions turned on, otherwise leave 0.0
            if(AerDust_Effect_opt.eq.0.or.AerDust_Effect_opt.eq.2) then
              TAREA_ij(l,num_aer) = 3.D0 * ODMDUST_ij(L,N) / ( eradius_ij(L,num_aer) * rho )
            endif
!... convert to approriate units for FastJX
            AERSP(LL,num_aer) = ODMDUST_ij(L,N) * gridBoxHeight_ij(L) * 1.0d3
          enddo
!
        enddo
!
!... remapped NSADaer (hydrophyllic)
        do N=1,NSADaer
          num_aer = num_aer+1
          kdry = iDRYwaer(N)
          rho = gmiDAA(4,kdry)*1.0D+3
          do L = k1,k2
            LL = L-k1+1
!... capture tropospheric aerosol parameters for diagnostic output
!... old method for aerosol surface area
            IF ( relHumidity_ij(L) <= tRH(2) ) THEN
               IRH = 1
             ELSE IF ( relHumidity_ij(L) <= tRH(3) ) THEN
               IRH = 2
             ELSE IF ( relHumidity_ij(L) <= tRH(4) ) THEN
               IRH = 3
             ELSE IF ( relHumidity_ij(L) <= tRH(5) ) THEN
               IRH = 4
             ELSE IF ( relHumidity_ij(L) <= tRH(6) ) THEN
               IRH = 5
             ELSE IF ( relHumidity_ij(L) <= tRH(7) ) THEN
               IRH = 6
             ELSE
               IRH = 7
            ENDIF
!... if sulfate, distinguish between strat and trop parameters?
!... trop
!            if(PPP(LL).ge.tropp_ij) then
            if(oldgmi_aero) then 
              ioffd = 1000   ! aertype > 1000 is use old GMI parameters
              IDXAER(LL,num_aer) = ioffd+kdry+IRH-1
            else
              IDXAER(LL,num_aer) = idxAtype(n)
            endif
!... strat SO4
!            else
!              if(oldgmi_aero) then 
!                IDXAER(LL,num_aer) = ioffd+kdry+IRH-1
!              else
!                IDXAER(LL,num_aer) = 1
!              endif
!            endif
            if(IRH.ge.NgmiRH_) then
              SCALEQ = gmiQW(N,NgmiRH_) / gmiQW(N,1)  !gmiQW(N,1) is dry extinction eff. for wAersl "N"
              REFF   = gmiRW(N,NgmiRH_)
            else
              FRAC = min( ((relHumidity_ij(L)-tRH(IRH)) / (tRH(IRH+1)-tRH(IRH))) ,1.0d0)
              scaleQ = (FRAC*gmiQW(N,IRH+1) + (1.d0-FRAC)*gmiQW(N,IRH)) / gmiQW(N,1)
              REFF   =  FRAC*gmiRW(N,IRH+1) + (1.d0-FRAC)*gmiRW(N,IRH)
            endif
            scaleR = REFF / (gmiRW(N,1))
!... calc diagnostic: effective radius
            eradius_ij (L,num_aer) = 1.0D-4 * REFF
!... calc diagnostic: hygroscopic growth
            HYGRO_ij(L,N) = scaleQ * scaleR * scaleR
!... het reactions turned off
            if(AerDust_Effect_opt.eq.0.or.AerDust_Effect_opt.eq.2) then
              TAREA_ij(L,num_aer) = 3.D0 * ODAER_ij(L,N) * scaleR**3 &
                                    / ( eradius_ij(L,num_aer) * rho )
            endif
!
!... convert aerosols to approriate units in photo_jx
            AERSP(LL,num_aer) = ODAER_ij(L,N) * gridBoxHeight_ij(L) * 1.0d3
!... end old method
!
          enddo
        enddo
!
!... BCphobic Aerosols
        num_aer = num_aer+1
        loc_naer = 6    ! UM-BC1  | 0.140 1.500| UMich w/Mie code logN:r=.050:s=.642  n=1.80+0.50i
        gmi_naer = iDRYwaer(2)   !
!
        if(oldgmi_aero) then 
          ioffd = 2000   ! aertype > 1000 is use old GMI parameters, aertype > 2000 is hydrophobic
          IDXAER(:,num_aer) = ioffd+gmi_naer
        else
          IDXAER(:,num_aer) = loc_naer
        endif
!... capture tropospheric aerosol parameters for diagnostic output
        rho = gmiDAA(4,gmi_naer)*1.0D+3
        do l = k1,k2
          ll = l-k1+1
!          ERADIUS_ij(ll,NSADdust+2) = gmiRAA(4,gmi_naer) * 1.0d-4
          if(ODcAER_ij(l,1).gt.0.0.and.(AerDust_Effect_opt.eq.0.or.AerDust_Effect_opt.eq.2)) then
             TAREA_ij(ll,NSADdust+2) = TAREA_ij(ll,NSADdust+2) &
                + 3.D0*ODcAER_ij(l,1) / (ERADIUS_ij(ll,NSADdust+2)*rho)
          endif
          AERSP(ll,num_aer) = ODcAER_ij(l,1) * gridBoxHeight_ij(l) * 1.0d3
        enddo
!
!... OCphobic Aerosols
        num_aer = num_aer+1
        loc_naer = 8    ! UM-BB08C| 0.149 1.230| UMich w/Mie code 8%BC 0%RH logN:r=.080:s=.500 n=1.552+0.04i
        gmi_naer = iDRYwaer(3)   !
!
        if(oldgmi_aero) then 
          IDXAER(:,num_aer) = ioffd+gmi_naer
        else
          IDXAER(:,num_aer) = loc_naer
        endif
!... capture tropospheric aerosol parameters for diagnostic output
        rho = gmiDAA(4,gmi_naer)*1.0D+3
        do l = k1,k2
          ll = l-k1+1
!          ERADIUS_ij(ll,NSADdust+3) = gmiRAA(4,gmi_naer) * 1.0d-4
!... het reactions turned off
          if(ODcAER_ij(l,2).gt.0.0.and.(AerDust_Effect_opt.eq.0.or.AerDust_Effect_opt.eq.2)) then
             TAREA_ij(ll,NSADdust+3) = TAREA_ij(ll,NSADdust+3)  &
                + 3.D0*ODcAER_ij(l,2) / (ERADIUS_ij(ll,NSADdust+3)*rho)
          endif
          AERSP(ll,num_aer) = ODcAER_ij(l,2) * gridBoxHeight_ij(l) * 1.0d3
        enddo
!... end aerosol initialization
      endif 
!
!!! set up clouds
!
      LTOP = k2-k1+2
!
!
!!! begin call to Cloud_J
!--CLDFLAG=   different cloud schemes (5:8 require max-ran overlap algorithm)
      if(do_clear_sky) then
        TCLDFLAG = 1 ! clear sky
        LWP(:) = 0.0d0
        IWP(:) = 0.0d0 
        REFFL(:) = 0.0d0 
        REFFI(:) = 0.0d0 
        CLF(:) = 0.0d0 
        CLDIW(:) = 0
        CLDCOR = 0.0d0 
      else
!        CLDFLAG = 1  !  clear sky
!        CLDFLAG = 2  !  grid-box avg clouds cloud: fract*(in cloud ODs) (minamal overlap?)
!        CLDFLAG = 3  !  cloud-fract**3/2*(in cloud ODs) (random overlap?) 
!        CLDFLAG = 4  !  NOT ALLOWED
!        CLDFLAG = 5  !  Random select NRANDO ICA's (Independent Column Atmos.) from all
!        CLDFLAG = 6  !  Use all (up to 4) quadrature cloud cover QCAs (mid-pts of bin)
!        CLDFLAG = 7  !  Use all (up to 4) QCAs (average clouds within each Q-bin)
!        CLDFLAG = 8  !  Calculate Js for ALL ICAs (up to 20,000 per cell!)
        TCLDFLAG = cldflag
!... cloud info ,need mean path and effective radii for ice and liquid given 
!...  gridbox cloud mixing ratios and cloud frac
!... cloud fraction
        CLDFRW(:) = 0.0d0
        CLDFRW(1:k2-k1+1) = cloudfrac_ij(k1:k2)
!
!... liquid cloud water content (kg/kg)
        CLDLWCW(:) = 0.0d0
        CLDLWCW(1:k2-k1+1) = ql_ij(k1:k2)
!... ice cloud water content
        CLDIWCW(:) = 0.0d0
        CLDIWCW(1:k2-k1+1) = qi_ij(k1:k2)
!... ice phase cloud particle effective radius
        CLDREFFI(:) = 0.0d0
        CLDREFFI(1:k2-k1+1) = ri_ij(k1:k2)
!... liquid cloud particle effective radius
        CLDREFFL(:) = 0.0d0
        CLDREFFL(1:k2-k1+1) = rl_ij(k1:k2)
!
        CLDIW(:) = 0
        CLF(:) = 0.d0
        WLC(:) = 0.d0
        WIC(:) = 0.d0
!... Do I need to adjust Anvil fractions for warm clouds?? see "subroutine RADCOUPLE"
!
        do L = 1,LTOP
          if (CLDFRW(L) .gt. 0.005d0) then
            CLF(L) = CLDFRW(L)
            WLC(L) = CLDLWCW(L) / CLF(L)
            WIC(L) = CLDIWCW(L) / CLF(L)
! ... CLDIW is an integer flag: 1 = water cld, 2 = ice cloud, 3 = both
            if (WLC(L) .gt. 1.d-11) CLDIW(L) = 1
            if (WIC(L) .gt. 1.d-11) CLDIW(L) = CLDIW(L) + 2
          else
            CLF(L) = 0.d0
            WLC(L) = 0.d0
            WIC(L) = 0.d0
          endif
        enddo
!
!... for LDRADIUS3:
!       NN3 = 50.*1.0e6
!... for LDRADIUS4:
!       NN_LAND  =  30.0e6*cnv_frc_ij + 150.0e6*(1.0-cnv_frc_ij)
!       NN_OCEAN =  30.0e6
!... Number Concentration adjustment Assumptions Over Land v Over Ocean
!       NN4 = frland_ij*NN_LAND + (1.0-frland_ij)*NN_OCEAN
!
!---derive R-effective for clouds:  the current UCI algorithm - use your own
        do L = 1,LTOP
!... temp variables for LDRADIUSn
          PDEL = (PPP(L-1) - PPP(L))*100.0d0
!         PMID = 0.5d0*(PPP(L-1)+PPP(L))*100.0d0        ! Pa
!         TTT4 = TTT(L)
!         WIC4 = WIC(L)
          WLC4 = WLC(L)
!---ice clouds
          if (WIC(L) .gt. 1.d-12) then
            IWP(L) = (1.0d3*WIC(L)*PDEL/MAPL_GRAV)  !*G100) ! / CLF(L) g/m2
!... GEOS routine, with conversion to microns
!           REFFI(L) = 1.0d6 * LDRADIUS3(PMID,TTT4,WIC4,NN3)
!           REFFI(L) = 1.0d6 * LDRADIUS4(PDEL,TTT4,WIC4,NN4,2)
            REFFI(L) = 1.0d6 * CLDREFFI(L)
          else
            IWP(L) = 0.d0
            REFFI(L) = 0.d0
          endif
!---water clouds
          if (WLC4 .gt. 1.d-12) then
            LWP(L) = (1.0d3*WLC4*PDEL/MAPL_GRAV)  !*G100) ! / CLF(L)  g/m2
!GEOS: taucld2=(((dp(k)*1.0e3)/MAPL_GRAV)*hydromets(k,2))*(awb_uv(1)+awb_uv(2)/reff(k,2))
!... GEOS routine, with conversion to microns
!           REFFL(L) = 1.0d6 * LDRADIUS3(PMID,TTT4,WLC4,NN3)
!           REFFL(L) = 1.0d6 * LDRADIUS4(PDEL,TTT4,WLC4,NN4,1)
            REFFL(L) = 1.0d6 * CLDREFFL(L)
          else
            LWP(L) = 0.d0
            REFFL(L) = 0.d0
          endif
        enddo
!
      endif
!
!... no orbital effects
!      baseSOLF = 1.d0
!... from orig SOLAR_JX - orbital effects in incoming radiation
      baseSOLF = 1.d0 - (0.034d0 * cos (dble (jday - 186) * 2.d0 * MAPL_PI_R8 / 365.d0) )
!
      do k=1,W_
!... solar cycle
        SOLF(k) = baseSOLF * fjx_solar_cycle_param(k)
      enddo
!
!... need to provide accurately
!... cos(SZA)
      U0 = cos(SZA*CPI180)
!... set to 8 or 12 to do tropo only
      NWBIN = 0
!
!================================================
!
! Cloud Correl Factor decreases with gap in G6 groups in ICA_ALL
      CLDCOR = 0.33  ! Prather recommended
!... from standalone input:  06   05 0.33   LNRG/RANDO/CLDCOR 
      LNRG = 6  ! Prather recommended
!... if using CLDFLAG = 5, NRANDO suggested by Jun Oh - IRAN should be ?
      NRANDO = 6
      IRAN = 2
!
!=======================================================================
!
      call CLOUD_JX (U0, SZA, RFL, SOLF, LPRTJ, PPP, ZZZ, TTT,  &
             DDD, RRR, OOO, LWP, IWP, REFFL, REFFI, CLF, CLDCOR, CLDIW,  &
             AERSP, IDXAER, LTOP, num_aer, VALJXX, num_qjs,  &
             TCLDFLAG, NRANDO, IRAN, LNRG, NICA, JCOUNT, cldOD_out, aerOD_out)
!     
!... send CLOUD_JX calcd cloud optical depth of 400nm back for diagnostic output w no FJX top layer
     cldOD_ij(k1:k2) = cldOD_out(1:k2-k1+1)
!
!... send aerosol optical depth of 400nm back for diagnostic output
     num_aer = 0
     do N=1,NSADdust
       num_aer = num_aer+1
       ODMDUST_ij(k1:k2,N) = aerOD_out(1:k2-k1+1,num_aer)
     enddo
     do N=1,NSADaer
        num_aer = num_aer+1
        ODAER_ij(k1:k2,N)  = aerOD_out(1:k2-k1+1,num_aer)
     enddo
!... add in hydrophobic BC and OC
     do N=1,2
        num_aer = num_aer+1
        ODcAER_ij(k1:k2,N) = aerOD_out(1:k2-k1+1,num_aer)
     enddo
!... map FastJX's Jrates to our order
      kall = chem_mask_khi-k1+1
      do n=1,num_qjs
!
!.sds.. as per previous 6.5 code:
!... * original fast-J/2 code had VALJ(2) as J[O3 -> O(3P)+O2]
!...     but now J[O3] is total O3 rate [O(1D) and O(3P)], correct
        if(JVMAP(n).eq.'O3') then
          qjgmi_ij(k1:chem_mask_khi,n) = (JFACTA(n) &
            * (VALJXX(1:kall,JINDO3) - VALJXX(1:kall,JINDO1D)))
!.sds.. as per previous 6.5 code:
!...     need to add both branches of Acet for GMI
        elseif(JVMAP(n).eq.'Acet-a') then
          qjgmi_ij(k1:chem_mask_khi,n) = (JFACTA(n) &
            * (VALJXX(1:kall,JINDAceta) + VALJXX(1:kall,JINDAcetb)))
        else
!.sds.. all others
          qjgmi_ij(k1:chem_mask_khi,n) = JFACTA(n) * VALJXX(1:kall,JIND(n))
        endif
      enddo
!
!----------------------------------------------
! Compute the Overhead ozone column diagnostics
!----------------------------------------------
      overheadO3col_ij(:) = 0.0d0
!      overheadO3col_ij(k2) = OOO(k2-k1+2)
      do l = k1, k2
         do k = l+1, k2+1
            overheadO3col_ij(l) = overheadO3col_ij(l) + OOO(k-k1+1)
         end do
      end do
!
      return
!
      end subroutine controlFastJX74
!EOP
!---------------------------------------------------------------------------
!BOP
!
! !IROUTINE: initializeFastJX74
!
! !INTERFACE:
!
      subroutine initializeFastJX74 (k1, k2, chem_mask_khi, num_qjs &
       , cross_section_file, cloud_scat_file, ssa_scat_file &
       , aer_scat_file, UMaer_scat_file, GMI_scat_file, T_O3_climatology_file &
       , H2O_CH4_climatology_file, cldflag, rootProc)
!
      implicit none
!
#     include "gmi_AerDust_const.h"
#     include "setkin_par.h"
#     include "setkin_fastj.h"
#     include "setkin_lchem.h"  
!
! !INPUT PARAMETERS:
      logical            , intent(in) :: rootProc
      integer            , intent(in) :: k1, k2
      integer            , intent(in) :: chem_mask_khi ! number of chemistry levels [JVL_]
      integer            , intent(in) :: num_qjs       ! number of photolysis reactions [JVN_]
      integer            , intent(in) :: cldflag       ! type of cloud OD calc for J rates
!... fast-J X-sections (spectral data) input file name
      character (len=*), intent(in) :: cross_section_file
!... CloudJ input files for scattering info
      character (len=*), intent(in) :: cloud_scat_file, ssa_scat_file, aer_scat_file &
        , UMaer_scat_file, GMI_scat_file
! Read in T & O3 climatology input file name general backup clim.
      character (len=*), intent(in) :: T_O3_climatology_file
! Read in H2O & CH4 climatology input file name general backup clim.
      character (len=*), intent(in) :: H2O_CH4_climatology_file
!... Local
      integer j, jj, k, JXUNIT
!
!----------------------------------------------------------------------
!
!---read in & store all fast-JX data:
!-----------------------------------------------------------------------
      call INIT_FJX (cross_section_file, cloud_scat_file, ssa_scat_file &
       , aer_scat_file, UMaer_scat_file, GMI_scat_file, T_O3_climatology_file &
       , H2O_CH4_climatology_file, rootProc)
!-----------------------------------------------------------------------
!
      JINDO3 = 0
      JINDO1D = 0
      JINDAceta = 0
      JINDAcetb = 0
!... need to capture some rates so that we can use for our chem,
!...   1) our J(O3) is for just O3P, UCI is total, we need to difference
!...   2) must combine both acetone branches for us
      do J = 1,NJX
        if (TITLEJX(J) .eq. 'O3') JINDO3 = J
        if (TITLEJX(J) .eq. 'O3(1D)') JINDO1D = J
        if (TITLEJX(J) .eq. 'Acet-a') JINDAceta = J
        if (TITLEJX(J) .eq. 'Acet-b') JINDAcetb = J
      enddo
!... exit if we fail to find needed cross-sections
      if(JINDO3 .eq. 0 .or. JINDO1D .eq. 0 .or. JINDAceta .eq. 0 .or. JINDAceta .eq. 0) &
        call EXITC(' CLOUD_JX: missing X-sections')
!
      JVMAP(1:NUM_J) = fastj_lookup(1:NUM_J)
      JFACTA(1:NUM_J) = fastj_yield(1:NUM_J)
      JIND(:) = 0
      NRATJ = NUM_J
!
!---Zero / Set index arrays that map Jvalue(j) onto our order
      do jj = 1,NUM_J
        do J = 1,NJX
          if (JVMAP(jj) .eq. TITLEJX(J)) JIND(jj) = J
        enddo
      enddo
!
      IF (rootProc) THEN
        write(6,'(a,i4,a,i2)')' CloudJ Photochemistry Scheme with',NRATJ,' J-values using cldflag=',cldflag
        do K=1,NRATJ
          if (JVMAP(K) .ne. '------' ) then
            J = JIND(K)
            if (J.eq.0) then
              write(6,'(i5,a9,f7.3,'' no mapping onto CloudJ'')') K,JVMAP(K),JFACTA(K)
            else
              write(6,'(i5,a9,f6.2,a18,i4,1x,a90)') J,JVMAP(K), JFACTA(K) &
                    , ' mapped to CloudJ:', K, lqjchem(K)
            endif
          endif
        enddo
      END IF
!
! Read in GMI aerosol scattering data
      JXUNIT  = 8
      call RD_GMI (JXUNIT,GMI_scat_file)
!
      return
!
      end subroutine initializeFastJX74
!BOC
!-----------------------------------------------------------------------
      subroutine RD_GMI(NUN,NAMFIL)
!-----------------------------------------------------------------------
!-------UMich aerosol optical data for fast-JX (ver 6.1+)
!-----------------------------------------------------------------------
!     NAMFIL   Name of scattering data file (e.g., FJX_scat.dat)
!     NUN      Channel number for reading data file
!-----------------------------------------------------------------------
      implicit none
!
      integer, intent(in) :: NUN
      character(*), intent(in) :: NAMFIL
!
      integer :: i, J, K
      integer :: NAA
      real*8  :: FWET 
!
      character*10 TITLE0
!      character (len=20) :: TITLAA (NgmiA_)
!
!
!...start execution
      open (NUN,FILE=NAMFIL,status='old',form='formatted')
!
      read(NUN,'(A10)') TITLE0
      read(NUN,'(A10,I5,/)') TITLE0, NAA
!      read(NUN,'(A10)') TITLE0  !... '/' above skips this line
!
!... Read aerosol phase functions:
      do j=1,NAA
        read(NUN,'(a10)') TITLE0
!
!e.g.:   300  3.0039  0.150  2500.0  0.884 1.000 2.030 2.363 2.082 1.710 1.290  0.941  0.638
       do k=1,NK_
          read(NUN,*) gmiWAA(k,j), gmiQAA(k,j), gmiRAA(k,j), gmiDAA(k,j) &
            , gmiSAA(k,j), (gmiPAA(i,k,j),i=1,8)
        enddo
      enddo
      close(NUN)
!... convert units
      gmiDAA(:,:) = gmiDAA(:,:)*1e-3
!
!... calc derived parameters
      DO k=1,Nwaer_
        DO i=1,NgmiRH_
! Wet radius in "jv_spec.dat"
          gmiRW(k,i) = gmiRAA(4,iDRYwaer(k)+i-1)
! Wet frac of aerosol
          FWET    = (gmiRW(k,i)**3 - gmiRW(k,1)**3) / gmiRW(k,i)**3
! Extinction efficiency Q for each RH bin
          gmiQW(k,i) = gmiQAA(4,iDRYwaer(k)+i-1)*FWET + gmiQAA(4,iDRYwaer(k))*(1.d0-FWET)
        ENDDO
      ENDDO

      return
!
      END SUBROUTINE RD_GMI
!------------------------------------------------------------------------------
!EOP
      END MODULE CloudJ_mod
!EOC
