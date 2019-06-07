! $Id$
      module fastj

      use GmiPrintError_mod, only : GmiPrintError
      use GmiASCIIoperations_mod, only : AsciiOpenRead

      private
      public :: Control_Fastj
      public :: InitializeFastj
      public :: GetQAA_RAAinFastj
      contains
!=============================================================================
!
!
! CODE DEVELOPER
!   Original code by Philip Cameron-Smith, LLNL
!   (Adapted from standalone.f in Fast-J v3 distributed by Oliver Wild, UCI)
!   cameronsmith1@llnl.gov
!
! FILE
!   fastj_control.F
!
! ROUTINES
!   Control_Fastj
!
! HISTORY
!   - January 28, 2005 * Jules Kouatchou
!     Added "ODAER_ij" and "ODMDUST_ij" as arguments of Control_Fastj.
!     Assigned ODAER_ij to "optaer" and "ODMDUST_ij" to "optdust".
!=============================================================================


!-----------------------------------------------------------------------------
!
! ROUTINE
!   Control_Fastj
!
! DESCRIPTION
!   This routine acts as the interface between the GMI model and Fast-J.
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
!     press3e_ij    : pressure at box edge (mb) [p]
!     kel_ij        : temperature at box centers (degK) [t]
!     optdepth_ij   : optical depth in box   (unitless) [od]
!     surf_alb_ij   : surface albedo     (fraction 0-1) [sa]
!     ODAER_ij      : optical depth for aerosol
!     ODMDUST_ij    : optical depth for mineral dust
!     sza_ij        : solar zenith angle
!   OUTPUTS
!     qjgmi_ij      : jvalues at grid centers (s^-1?)   [zpj]
!     n_qj_O3_2OH   : index of O3 -> O1d photolysis IF it is actually
!                        intended to be corrected to be O3 -> 2OH photolysis
!
!-----------------------------------------------------------------------------

      subroutine Control_Fastj  &
     &  (k1, k2, chem_mask_khi, num_qjs,  &
     &   month_gmi, jday, time_sec, fastj_offset_sec, sza_ij,  &
     &   press3e_ij, kel_ij, optdepth_ij,  &
     &   surf_alb_ij, qjgmi_ij, overheadO3col_ij, ODAER_ij, ODMDUST_ij)

      implicit none

#     include "gmi_AerDust_const.h"
#     include "setkin_par.h"
#     include "fastj_cmn_h.h"
#     include "fastj_cmn_t.h"
#     include "fastj_cmn_w.h"
#     include "fastj_jv_cmn.h"


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
      real*8,  intent(in)  ::  press3e_ij(k1-1:k2)
      real*8,  intent(in)  ::  kel_ij(k1:k2)
      real*8,  intent(in)  ::  optdepth_ij(k1:k2)
      real*8,  intent(in)  ::  surf_alb_ij
!
      real*8,  intent(in)  ::  ODAER_ij(k1:k2,NSADaer*nrh_b)
      real*8,  intent(in)  ::  ODMDUST_ij(k1:k2,NSADdust)
!
      real*8,  intent(out) ::  overheadO3col_ij(k1:k1)
      real*8,  intent(out) ::  qjgmi_ij(k1:chem_mask_khi, num_qjs)


!     ----------------------
!     Variable declarations.
!     ----------------------

      character (len=75) :: err_msg

      logical, save :: first = .true.

      integer :: il, ik

      real*8  :: pi
      real*8  :: timej

      real*8  :: zpj(chem_mask_khi-k1+1, num_qjs)


!     ----------------
!     Begin execution.
!     ----------------

!c    if (pr_diag) then
!c      Write (6,*) 'Control_Fastj called by ', loc_proc
!c    end if


!     -------------------------
!     Just do one column {PJC}.
!     -------------------------

!      nslat = 1
!      nslon = 1

      pi = 3.141592653589793d0  ! FLAG {PJC}


!     -----------------------------------------------------
!     Check that k1=1, otherwise Fast-J may not work {PJC}.
!     -----------------------------------------------------

      if (k1 /= 1) then
        err_msg =  &
     &    'Control_Fastj: Fast-J may not work correctly if k1 /= 1.'
        call GmiPrintError (err_msg, .true., 1, k1, 0, 0, 0.0d0, 0.0d0)
      end if

!     ----------------------------------------------------------------
!     Check dimensions against dimension maximums specified in include
!     files {PJC}.
!     ----------------------------------------------------------------

      if (lpar > LPAR_MAX) then
        err_msg = 'Control_Fastj: lpar exceeds LPAR_MAX. '//  &
     &            'Increase LPAR_MAX and recompile.'
        call GmiPrintError  &
     &    (err_msg, .true., 2, lpar, LPAR_MAX, 0, 0.0d0, 0.0d0)
      end if

      if (jpnl > JPNL_MAX) then
        err_msg = 'Control_Fastj: jpnl exceeds JPNL_MAX. '//  &
     &            'Increase JPNL_MAX and recompile.'
        call GmiPrintError  &
     &    (err_msg, .true., 2, jpnl, JPNL_MAX, 0, 0.0d0, 0.0d0)
      end if

      if (jppj > JPPJ_MAX) then
        err_msg = 'Control_Fastj: jppj exceeds JPPJ_MAX. '//  &
     &            'Increase JPPJ_MAX and recompile.'
        call GmiPrintError  &
     &    (err_msg, .true., 2, jppj, JPPJ_MAX, 0, 0.0d0, 0.0d0)
      end if

      if (nb > NB_MAX) then
        err_msg = 'Control_Fastj: nb exceeds NB_MAX. '//  &
     &            'Increase NB_MAX and recompile.'
        call GmiPrintError  &
     &    (err_msg, .true., 2, nb,   NB_MAX,   0, 0.0d0, 0.0d0)
      end if

      if (nc > NC_MAX) then
        err_msg = 'Control_Fastj: nc exceeds NC_MAX. '//  &
     &            'Increase NC_MAX and recompile.'
        call GmiPrintError  &
     &    (err_msg, .true., 2, nc,   NC_MAX,   0, 0.0d0, 0.0d0)
      end if


!     -------------------------------------------------
!     Set up Fast-J variables from GMI variables {PJC}.
!     -------------------------------------------------

      month = month_gmi
      iday  = jday
      tau   = time_sec         / 3600.0d0     ! convert to hours {PJC}
      timej = fastj_offset_sec / 3600.0d0     ! convert to hours {PJC}

      !xgrd (nslon) = londeg_i * pi / 180.0d0
      !ygrd (nslat) = latdeg_j * pi / 180.0d0
      !ydgrd(nslat) = latdeg_j

      pj(1:NB)                 = press3e_ij(k1-1:k2)
      pj(NB+1)                 = 0.0d0

      t (nslon, nslat, 1:lpar) = kel_ij     (k1:k2)
      od(nslon, nslat, 1:lpar) = optdepth_ij(k1:k2)
      sa(nslon, nslat)         = surf_alb_ij

      ! Aerosol OD profile [unitless]
      !optaer(nslon, nslat, 1:lpar,:)  = ODAER_ij(k1:k2,:)
       optaer(1:lpar,:)  = ODAER_ij(k1:k2,:)

      ! Mineral dust OD profile [unitless]
      ! optdust(nslon, nslat, 1:lpar,:) =ODMDUST_ij(k1:k2,:)
       optdust(1:lpar,:) =ODMDUST_ij(k1:k2,:)

      zpj(:,:) = -1.0d0

!     ----------------------------------
!     Call main Fast-J routine (Photoj).
!     ----------------------------------

!c    Open (2, file = 'output')


!c    Write (6,900) 1, iday, tau


!     --------------------------------------------------------------------
!     Call main Fast-J routine, which actually calulates photolysis rates.
!     --------------------------------------------------------------------

!c    Write (6,*) 'Calling Photoj.'

      sza = sza_ij    ! Oloso added sza_ij 12/2010

!     ===========
      call Photoj (zpj, timej)
!     ===========

      qjgmi_ij(:,:) = zpj(:,:)
!c    qjgmi_ij(:,:) = 0.0d0
!!!!
!!!! Overhead ozone column
!!!! 
      overheadO3col_ij = 0.0d0
      overheadO3col_ij(k2) = DO3(k2)
      do il = k1, k2-1
         do ik = il+1, k2
            overheadO3col_ij(il) = overheadO3col_ij(il) + DO3(ik)
         end do
      end do 


      return
      end subroutine Control_Fastj

! ---------------------------------------------------------------------
! ---------------------------------------------------------------------

      subroutine InitializeFastj &
     &    (cross_section_file, rate_file, T_O3_climatology_file, &
     &     n_qj_O3_2OH, num_qjs, chem_mask_khi, k2, k1)

      implicit none

      integer                         :: n_qj_O3_2OH, num_qjs, chem_mask_khi, k2, k1
      character (len=128), intent(in) :: cross_section_file
      character (len=128), intent(in) :: rate_file
      character (len=128), intent(in) :: T_O3_climatology_file

#     include "fastj_cmn_h.h"
#     include "fastj_cmn_t.h"
#     include "fastj_jv_cmn.h"

      nslat = 1
      nslon = 1

!       ----------------------------------------
!       Initial call to Fast-J to set things up.
!       ----------------------------------------

      lpar = k2 - k1 + 1

      jpnl = chem_mask_khi
      jppj = num_qjs

      nb   = lpar + 1
      nc   = 2 * nb

        n_qj_O3_2OH = 0

!       ===========
        call Inphot  &
!       ===========
     &    (cross_section_file, rate_file,  &
     &     T_O3_climatology_file, n_qj_O3_2OH)

      return

      end subroutine InitializeFastj

! ---------------------------------------------------------------------

      subroutine GetQAA_RAAinFastj (bRAA, bQAA)

      implicit none

#     include "fastj_cmn_h.h"
#     include "fastj_jv_cmn.h"

      real*8 , intent(out) :: bRAA(:,:), bQAA(:,:)

      bRAA(:,:) = RAA(:,:)
      bQAA(:,:) = QAA(:,:)

      return

      end subroutine GetQAA_RAAinFastj

! ---------------------------------------------------------------------
#     include "newcol.code"

      end module fastj

