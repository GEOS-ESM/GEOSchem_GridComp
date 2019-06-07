      module fast_JX

      use GmiPrintError_mod, only : GmiPrintError
      use GmiASCIIoperations_mod, only : AsciiOpenRead
      USE GmiFastJX_ParametersMod
      USE GmiFastJX_includeMod, ONLY : t_fastJXbundle

      implicit none

      private
      public :: Control_Fast_JX
      public :: InitializeFastJX
      public :: GetQAA_RAAinFastJX
      contains
!=============================================================================
!
! $Id$
!
! CODE DEVELOPER
!   Original code by Philip Cameron-Smith, LLNL
!   (Adapted from standalone.f in Fast-J v3 distributed by Oliver Wild, UCI)
!   cameronsmith1@llnl.gov
!   Modifeid by A. Oloso, GSFC, 8/04, for Fast-JX v2 distributed by M. Prather, UCI, 7/04)
!
! FILE
!   fast_JX_control.F
!
! ROUTINES
!   Control_Fast_JX
!
! HISTORY
!   - March 9, 2005 * Jules Kouatchou
!     The variable "ozone_ij" is now an optional arguments. When it is available,
!     the code uses the model climatology and the variable "do3" is set.
!     The variables "ODAER_ij" and "ODMDUST_ij" were added as arguments.
!   - May 25, 2016 * Luke Oman
!     Now use Pressures on Edges instead of on Centers (fixed long-standing bug)
!
!=============================================================================


!-----------------------------------------------------------------------------
!
! ROUTINE
!   Control_Fast_JX
!
! DESCRIPTION
!   This routine acts as the interface between the GMI model and Fast-JX.
!
! ARGUMENTS
!
!   INPUTS
!     cross_section_file    : X-Section quantum yield
!     rate_file             : Master rate file
!     T_O3_climatology_file : T & O3 climatology
!     k1            : first interface index [lpar = k2 - k1 + 1]
!     k2            : last  interface index [lpar = k2 - k1 + 1]
!     chem_mask_khi : number of chemistry levels     [jpnl]
!     num_qjs       : number of photolysis reactions [jppj]
!     month_gmi     : number of month (1- 12)        [month]
!     jday          : day    of year  (1-365)        [iday]
!     time_sec      : time of day in model (s, GMT)  [tau (hrs)]
!     fastj_offset_sec : offset from tau at which to do photolysis (s)
!                        [timej (hrs)]
!     londeg_i      : longitude (midpoint, deg) [xgrd (rad)]
!     latdeg_j      : latitude  (midpoint, deg) [ygrd (rad), ydgrd (deg)]
!     am            : hybrid am levels (mid-point of level) [etaa = pt*am]
!     pt            : scale factor for am's [etaa = pt*am]
!     bm            : hybrid bm levels (mid-point of level)      [etab]
!     pctm_ij       : surface pressure (mb) [p]
!     kel_ij        : temperature at box centers (degK) [t]
!     optdepth_ij   : optical depth in box   (unitless) [od]
!     surf_alb_ij   : surface albedo     (fraction 0-1) [sa]
!     ODAER_ij      : optical depth for aerosol
!     ODMDUST_ij    : optical depth for mineral dust
!     ozone_ij      : mixing ratio for ozone
!       sza_ij      : solar zenith angle
!   OUTPUTS
!     qjgmi_ij      : jvalues at grid centers (s^-1?)   [zpj]
!     n_qj_O3_2OH : index of O3 -> O1d photolysis IF it is actually
!                        intended to be corrected to be O3 -> 2OH photolysis
!
!-----------------------------------------------------------------------------

      subroutine Control_Fast_JX  &
     &   (JXbundle, k1, k2, chem_mask_khi, num_qjs,  &
     &   month_gmi, jday, time_sec, fastj_offset_sec, sza_ij, &
     &   press3e_ij, pctm_ij, kel_ij, optdepth_ij,  &
     &   surf_alb_ij, qjgmi_ij, overheadO3col_ij,  &
     &   ODAER_ij, ODMDUST_ij, ozone_ij)


      implicit none

# include "gmi_AerDust_const.h"
# include "setkin_par.h"

!     ----------------------
!     Argument declarations.
!     ----------------------
      integer, intent(in)  ::  k1
      integer, intent(in)  ::  k2
      integer, intent(in)  ::  chem_mask_khi
      integer, intent(in)  ::  num_qjs
      integer, intent(in)  ::  month_gmi
      integer, intent(in)  ::  jday
      real*8,  intent(in)  ::  time_sec
      real*8,  intent(in)  ::  fastj_offset_sec
      real*8,  intent(in)  ::  sza_ij
      real*8,  intent(in)  ::  press3e_ij(k1:k2)
      real*8,  intent(in)  ::  pctm_ij
      real*8,  intent(in)  ::  kel_ij(k1:k2)
      real*8,  intent(in)  ::  optdepth_ij(k1:k2)
      real*8,  intent(in)  ::  surf_alb_ij
      real*8,  intent(out) ::  overheadO3col_ij(k1:k2)
      real*8,  intent(out) ::  qjgmi_ij(k1:chem_mask_khi, num_qjs)
      real*8,  intent(in)  ::  ODAER_ij(k1:k2,NSADaer*nrh_b)
      real*8,  intent(in)  ::  ODMDUST_ij(k1:k2,NSADdust)
      real*8, intent(in), optional ::  ozone_ij   (k1:k2)
      TYPE(t_fastJXbundle), intent(inOut) :: JXbundle
!
! !LOCAL VARIABLES:
      integer :: il, ik
      real*8  :: pi
      real*8  :: timej
      real*8  ::  amm(k1-1:k2)
      real*8  ::  bmm(k1-1:k2)
      real*8  :: zpj(chem_mask_khi-k1+1, num_qjs)
!--------------------------------------------------------------------------
!BOC
      pi = 3.141592653589793d0  ! FLAG {PJC}

!     -------------------------------------------------
!     Set up Fast-JX variables from GMI variables {PJC}.
!     -------------------------------------------------

      if ( present(ozone_ij) ) then
         JXbundle%do_model_clima = .true.
      else
         JXbundle%do_model_clima = .false.
      endif

      JXbundle%month = month_gmi
      JXbundle%iday  = jday
      JXbundle%tau   = time_sec         / 3600.0d0     ! convert to hours {PJC}
      timej = fastj_offset_sec / 3600.0d0     ! convert to hours {PJC}

      JXbundle%pj(1)                    = pctm_ij
      JXbundle%pj(2:JXbundle%NB)                 = press3e_ij(k1:k2)
      JXbundle%pj(JXbundle%NB+1)                 = 0.0d0

!      p (nslon, nslat)         = pctm_ij
      JXbundle%t (JXbundle%nslon, JXbundle%nslat, 1:JXbundle%lpar) = kel_ij     (k1:k2)
      JXbundle%od(JXbundle%nslon, JXbundle%nslat, 1:JXbundle%lpar) = optdepth_ij(k1:k2)
      if (JXbundle%do_model_clima) then
         JXbundle%do3(1:JXbundle%lpar)              = ozone_ij(k1:k2)
      endif
      JXbundle%sa(JXbundle%nslon, JXbundle%nslat)         = surf_alb_ij

      ! Aerosol OD profile [unitless]
      !optaer(nslon, nslat, 1:lpar,:)  = ODAER_ij(k1:k2,:)
       JXbundle%optaer(1:JXbundle%lpar,:)  = ODAER_ij(k1:k2,:)

      ! Mineral dust OD profile [unitless]
      ! optdust(nslon, nslat, 1:lpar,:) =ODMDUST_ij(k1:k2,:)
       JXbundle%optdust(1:JXbundle%lpar,:) =ODMDUST_ij(k1:k2,:)

      zpj(:,:) = -1.0d0

!     --------------------------------------------------------------------
!     call main Fast-JX routine, which actually calculates photolysis rates.
!     --------------------------------------------------------------------

      JXbundle%sza = sza_ij    ! Oloso added sza_ij 12/2010

!     ===========
      call Photo (JXbundle, zpj, timej)
!     ===========

      qjgmi_ij(:,:) = zpj(:,:)

!!!!
!!!! Overhead ozone column
!!!!
      overheadO3col_ij = 0.0d0
      overheadO3col_ij(k2) = JXbundle%DO3(k2)
      do il = k1, k2-1
         do ik = il+1, k2
            overheadO3col_ij(il) = overheadO3col_ij(il) + JXbundle%DO3(ik)
         end do
      end do

      return

      end subroutine Control_Fast_JX

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine InitializeFastJX (JXbundle, &
               cross_section_file,rate_file, T_O3_climatology_file, &
     &     rootProc,n_qj_O3_2OH, num_qjs, chem_mask_khi, k2, k1)

      implicit none

      integer  :: n_qj_O3_2OH, num_qjs, chem_mask_khi, k2, k1
      character (len=128), intent(in) :: cross_section_file
      character (len=128), intent(in) :: rate_file
      character (len=128), intent(in) :: T_O3_climatology_file
      logical,             intent(in) :: rootProc
      TYPE(t_fastJXbundle), intent(inOut) :: JXbundle

      JXbundle%nslat = 1
      JXbundle%nslon = 1

!     -------------------------
!     Set up Fast-JX dimensions.
!     -------------------------

      JXbundle%lpar = k2 - k1 + 1

      JXbundle%jpnl = chem_mask_khi
      JXbundle%jppj = num_qjs

      JXbundle%nb   = JXbundle%lpar + 1
      JXbundle%nc   = 2 * JXbundle%nb

!       ----------------------------------------
!       Initial call to Fast-JX to set things up.
!       ----------------------------------------
!       ===========
        call Inphot  &
!       ===========
     &    (JXbundle, rootProc,cross_section_file,rate_file,  &
     &     T_O3_climatology_file, n_qj_O3_2OH)

      return 

      end subroutine InitializeFastJX

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine GetQAA_RAAinFastJX (JXbundle, bRAA, bQAA)

!      USE GmiFastJX_includeMod, ONLY : RAA, QAA

      implicit none

      real*8 , intent(out) :: bRAA(:,:), bQAA(:,:)
      TYPE(t_fastJXbundle), intent(inOut) :: JXbundle

      bRAA(:,:) = JXbundle%RAA(:,:)
      bQAA(:,:) = JXbundle%QAA(:,:)

      return

      end subroutine GetQAA_RAAinFastJX

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# include "newcol.code"

      end module fast_JX
