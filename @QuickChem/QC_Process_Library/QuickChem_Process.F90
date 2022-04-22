#define __SUCCESS__ 0
#define __FAIL__ 1
#define __VERIFY__(x) if(x/=0) then; if(present(rc)) rc=x; return; endif
#define __VERIFY_NO_OPT__(x) if(x/=0) then; rc=x; return; endif
#define __RC__ rc=status); __VERIFY__(status
#define __RC_NO_OPT__ rc=status); __VERIFY_NO_OPT__(status
#define __STAT__ stat=status); __VERIFY__(status
#define __IOSTAT__ iostat=status); __VERIFY__(status
#define __RETURN__(x) if (present(rc)) rc=x; return
#define __ASSERT__(expr) if(.not. (expr)) then; if (present(rc)) rc=-1; return; endif
!-------------------------------------------------------------------------
!
! !MODULE: QuickChem_Process -- QuickChem process library
!
! !INTERFACE:
   module  QuickChem_Process

! !USES:
!  Only instrinsic fortran types and functions are allowed.
   use, intrinsic :: iso_fortran_env, only: IOSTAT_END

   implicit none
   private

!
! !PUBLIC MEMBER FUNCTIONS:
!

!! THESE ARE LEFT OVER FROM GOCART2G:
   public Chem_UtilResVal     !!  Useful for Resolution Vectors
!  public Chem_UtilIdow  !! integer day of week
!  public Chem_UtilCdow  !! char    day of week
!  public Chem_BiomassDiurnal  !! if diurnal is needed, see G2G
!  public ReadPointEmissions   !! may be useful
!  public EmissionReader    !! for point emissions (ASCII files)

!  type :: EmissionReader
!     private
!     integer, allocatable :: unit
!  contains
!     procedure :: open
!     procedure :: close
!     procedure :: rewind => rewind_reader
!     procedure :: is_end_marker
!     procedure :: read_table
!     procedure :: next_line
!     procedure :: count_words
!     procedure :: scan_to_label
!     procedure :: get_dims
!  end type EmissionReader

!  type KeywordEnforcer
!  end type KeywordEnforcer

!
! !DESCRIPTION:
!
!  This module contains and implements all necessary process calculations for GOCART.
!
! !REVISION HISTORY:
!
!  11Feb2020  E.Sherman, A.da Silva, T.Clune, A.Darmenov - Ported/consolidated/refactored GOCART
!                   physics and chemistry code into a single process library that only uses
!                   intrinsic Fortran functions.
!
!  01Apr2021  R.Montuoro/NOAA - Added FENGSHA dust scheme and related methods.
!
!
!EOP
!-------------------------------------------------------------------------
CONTAINS


!==================================================================================
!BOP
! !IROUTINE: idaynum

   integer function idaynum (nymd)

! !USES:
   implicit NONE

! !INPUT PARAMETERS:
  integer :: nymd

! !OUTPUT PARAMETERS:

! !DESCRIPTION: Given nymd compute the day number of the year.

!
! !REVISION HISTORY:
! 29July2004 P.Colarco - Legacy code
! 23July2020 E.Sherman - moved from SulfateChemDriverMod.F90 for use in process library.

! !Local Variables

   integer :: yyyy, mm, dd, imon, isleapyr
   integer :: ndays(12)

   data ndays /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/

   yyyy = nymd / 10000
   mm = mod(nymd,10000) / 100
   dd = mod(nymd,100)

!EOP
!-------------------------------------------------------------------------
!  Begin...

!  Is it a leap year?
   isleapyr = 0
   if(mod(yyyy,4) .eq. 0) then
    isleapyr = 1
    if(mod(yyyy,100) .eq. 0) then
     isleapyr = 0
     if(mod(yyyy,400) .eq. 0) then
      isleapyr = 1
     endif
    endif
   endif

!  What day number
   idaynum = 0
   if(mm .eq. 1) then
    idaynum = dd
   else
    do imon = 1, mm-1
     if(imon .eq. 2 .and. isleapyr .eq. 1) then
      idaynum = idaynum+29
     else
      idaynum = idaynum + ndays(imon)
     endif
    enddo
    idaynum = idaynum + dd
   endif

   return
   end function idaynum

!==================================================================================
!BOP
!
! !IROUTINE:  Chem_UtilResVal --- returns resolution dependent value
!
! !INTERFACE:
!
   function Chem_UtilResVal( im_World, jm_World, res_value, rc ) result (val)

! !USES:

   implicit NONE

   real :: val                                ! resolution dependent value

! !INPUT/OUTPUT PARAMETERS:
   integer, intent(in) :: im_World, jm_World  ! number of global grid cells
   real,    intent(in) :: res_value(:)        ! array with the resolution dependent values:
                                              ! the 'a', 'b', ..., 'e' resolution values have 
                                              ! indexes 1, 2, ..., 5.

! !OUTPUT PARAMETERS:
   integer, intent(inout) :: rc               ! return code


! !DESCRIPTION: 
!
! !REVISION HISTORY:
!
! 13 Feb2012   Anton Darmenov  First crack.
! 25 Oct2012   Anton Darmenov  Added support for FV3 resolutions.
! 19 Aug2020   E. Sherman - moved from Chem_UtilMod.F90 to process library
!
!EOP
!-------------------------------------------------------------------------
       character(len=*), parameter :: Iam = 'Chem_UtilResVal'

       integer            :: i_res
       integer, parameter :: res_a = 1  ! 'a' to 'e' resolution indexes
       integer, parameter :: res_b = 2  !
       integer, parameter :: res_c = 3  !
       integer, parameter :: res_d = 4  !
       integer, parameter :: res_e = 5  !
       integer, parameter :: res_f = 6  !

       i_res = 0

       if ((im_World < 1) .or. (jm_World < 1)) then
!           call die(Iam, 'incorrect model resolution')
           print*,'QuickChem_Process::Chem_UtilResVal - incorrect model resolution'
           return
       end if

       if (jm_World == 6*im_World) then
           if (im_World <= 24) then
               i_res = res_a
           else if (im_World <=  48) then
               i_res = res_b
           else if (im_World <=  90) then
               i_res = res_c
           else if (im_World <= 180) then
               i_res = res_d
           else if (im_World <= 360) then
               i_res = res_e
           else if (im_World <= 720) then
               i_res = res_f
           else
               i_res = res_f
           end if
       else
           if ((im_World <= 72) .and. (jm_World <= 46)) then
               i_res = res_a
           else if ((im_World <=  144) .and. (jm_World <=  91)) then
               i_res = res_b
           else if ((im_World <=  288) .and. (jm_World <= 181)) then
               i_res = res_c
           else if ((im_World <=  576) .and. (jm_World <= 361)) then
               i_res = res_d
           else if ((im_World <= 1152) .and. (jm_World <= 721)) then
               i_res = res_e
           else if ((im_World <= 2304) .and. (jm_World <=1441)) then
               i_res = res_f
           else
               i_res = res_f
           end if


       end if

       if ((i_res < 1) .or. (i_res > size(res_value))) then
           val = 0.0
           rc  = __FAIL__
       else
           val = res_value(i_res)
           rc  = __SUCCESS__
       end if

   end function Chem_UtilResVal

!!==================================================================================
!
!   function Chem_UtilIdow(nymd) result (idow)
!     implicit NONE
!     integer, intent(in) :: nymd
!     integer :: idow ! day of the week: Sun=1, Mon=2, etc.
!     integer :: y, m, d
!     integer, parameter :: t(0:11) = (/ 0, 3, 2, 5, 0, 3, 5, 1, 4, 6, 2, 4 /)
!     y = nymd / 10000
!     m = (nymd - y*10000)/100
!     d = nymd - (y*10000 + m*100)
!     if ( m<3 ) then
!        y = y - 1
!     end if
!     idow = 1+mod(y + y/4 - y/100 + y/400 + t(m-1) + d,7)
!     return
!   end function Chem_UtilIdow
!
!   function Chem_UtilCdow(nymd) result (cdow)
!     implicit NONE
!     integer, intent(in) :: nymd
!     character(len=3) :: cdow ! day of the week: Sun, Mon, etc.
!     character(len=3) :: cday(7) = (/ 'Sun','Mon', 'Tue', 'Wed', 'Thu', 'Fri', 'Sat' /)
!     cdow = cday(Chem_UtilIdow(nymd))
!     return
!   end function Chem_UtilCdow


!!==================================================================================
!
!!BOP
!!
!! !ROUTINE:  Chem_BiomassDiurnal - Applies diurnal cycle to biomass emissions.
!!
!! !INTERFACE:
!     subroutine Chem_BiomassDiurnal ( Eout, Ein, lons, lats, nhms, cdt)
!
!! !USES:
!
!  IMPLICIT NONE
!
!! !ARGUMENTS:
!
!       real, intent(out)   :: Eout(:,:) ! Emissions valid at NHMS
!       real, intent(in)    :: Ein(:,:)  ! Daily-mean emissions
!       real, intent(in)    :: lons(:,:) ! Latitudes in degrees
!       real, intent(in)    :: lats(:,:) ! Latitudes in degrees
!       integer, intent(in) :: nhms
!       real, intent(in)    :: cdt       ! time step in seconds
!
!! !DESCRIPTION: 
!!
!!      Applies diurnal cycle to biomass emissions.       
!!
!! !DESCRIPTION:
!!
!!  This module implements assorted odds & ends for fvChem.
!!
!! !REVISION HISTORY:
!!
!!  13nov2009  da Silva  First crack.
!!  19Aug2020  E. Sherman - moved from Chem_UtilMod.F90 to process library
!!
!!EOP
!!-------------------------------------------------------------------------
!
!!      Hardwired diurnal cycle (multiplied by 100)
!!      These numbers were derived from GOES-12
!!      fire counts for 2003-2007.
!!      -------------------------------------------
!       integer, parameter :: N = 240
!       real,    parameter :: DT = 86400. / N
!
!!      Apply flat diurnal cycle for boreal forests as a 
!!      temporary solution to prevent very high aerosol
!!      optical depth during the day
!       real,    parameter :: Boreal(N) = 1.0
!!      real,    parameter :: Boreal(N) = &
!!      (/ 0.0277, 0.0292, 0.0306, 0.0318, 0.0327, 0.0335, &
!!         0.0340, 0.0342, 0.0341, 0.0338, 0.0333, 0.0326, &
!!         0.0316, 0.0305, 0.0292, 0.0278, 0.0263, 0.0248, &
!!         0.0233, 0.0217, 0.0202, 0.0187, 0.0172, 0.0158, &
!!         0.0145, 0.0133, 0.0121, 0.0110, 0.0100, 0.0091, &
!!         0.0083, 0.0075, 0.0068, 0.0062, 0.0056, 0.0051, &
!!         0.0046, 0.0042, 0.0038, 0.0035, 0.0032, 0.0030, &
!!         0.0028, 0.0026, 0.0025, 0.0024, 0.0024, 0.0024, &
!!         0.0024, 0.0026, 0.0027, 0.0030, 0.0033, 0.0036, &
!!         0.0041, 0.0046, 0.0052, 0.0060, 0.0069, 0.0079, &
!!         0.0090, 0.0104, 0.0119, 0.0137, 0.0157, 0.0180, &
!!         0.0205, 0.0235, 0.0268, 0.0305, 0.0346, 0.0393, &
!!         0.0444, 0.0502, 0.0565, 0.0634, 0.0711, 0.0794, &
!!         0.0884, 0.0982, 0.1087, 0.1201, 0.1323, 0.1453, &
!!         0.1593, 0.1742, 0.1900, 0.2069, 0.2249, 0.2439, &
!!         0.2642, 0.2858, 0.3086, 0.3329, 0.3587, 0.3860, &
!!         0.4149, 0.4455, 0.4776, 0.5115, 0.5470, 0.5840, &
!!         0.6227, 0.6628, 0.7043, 0.7470, 0.7908, 0.8355, &
!!         0.8810, 0.9271, 0.9735, 1.0200, 1.0665, 1.1126, &
!!         1.1580, 1.2026, 1.2460, 1.2880, 1.3282, 1.3664, &
!!         1.4023, 1.4356, 1.4660, 1.4933, 1.5174, 1.5379, &
!!         1.5548, 1.5679, 1.5772, 1.5826, 1.5841, 1.5818, &
!!         1.5758, 1.5661, 1.5529, 1.5365, 1.5169, 1.4944, &
!!         1.4693, 1.4417, 1.4119, 1.3801, 1.3467, 1.3117, &
!!         1.2755, 1.2383, 1.2003, 1.1616, 1.1225, 1.0832, &
!!         1.0437, 1.0044, 0.9653, 0.9265, 0.8882, 0.8504, &
!!         0.8134, 0.7771, 0.7416, 0.7070, 0.6734, 0.6407, &
!!         0.6092, 0.5787, 0.5493, 0.5210, 0.4939, 0.4680, &
!!         0.4433, 0.4197, 0.3974, 0.3763, 0.3565, 0.3380, &
!!         0.3209, 0.3051, 0.2907, 0.2777, 0.2662, 0.2561, &
!!         0.2476, 0.2407, 0.2352, 0.2313, 0.2289, 0.2279, &
!!         0.2283, 0.2300, 0.2329, 0.2369, 0.2417, 0.2474, &
!!         0.2536, 0.2602, 0.2670, 0.2738, 0.2805, 0.2869, &
!!         0.2927, 0.2979, 0.3024, 0.3059, 0.3085, 0.3101, &
!!         0.3107, 0.3102, 0.3087, 0.3061, 0.3026, 0.2983, &
!!         0.2931, 0.2871, 0.2806, 0.2735, 0.2659, 0.2579, &
!!         0.2497, 0.2412, 0.2326, 0.2240, 0.2153, 0.2066, &
!!         0.1979, 0.1894, 0.1809, 0.1726, 0.1643, 0.1562, &
!!         0.1482, 0.1404, 0.1326, 0.1250, 0.1175, 0.1101, &
!!         0.1028, 0.0956, 0.0886, 0.0818, 0.0751, 0.0687 /)       
!       real,    parameter :: NonBoreal(N) = &
!       (/ 0.0121, 0.0150, 0.0172, 0.0185, 0.0189, 0.0184, &
!          0.0174, 0.0162, 0.0151, 0.0141, 0.0133, 0.0126, &
!          0.0121, 0.0117, 0.0115, 0.0114, 0.0114, 0.0116, &
!          0.0120, 0.0126, 0.0133, 0.0142, 0.0151, 0.0159, &
!          0.0167, 0.0174, 0.0180, 0.0184, 0.0187, 0.0189, &
!          0.0190, 0.0190, 0.0191, 0.0192, 0.0192, 0.0193, &
!          0.0194, 0.0194, 0.0193, 0.0192, 0.0190, 0.0187, &
!          0.0185, 0.0182, 0.0180, 0.0178, 0.0177, 0.0176, &
!          0.0174, 0.0172, 0.0169, 0.0166, 0.0162, 0.0158, &
!          0.0153, 0.0149, 0.0144, 0.0138, 0.0132, 0.0126, &
!          0.0118, 0.0109, 0.0101, 0.0092, 0.0085, 0.0081, &
!          0.0080, 0.0083, 0.0091, 0.0102, 0.0117, 0.0135, &
!          0.0157, 0.0182, 0.0210, 0.0240, 0.0273, 0.0308, &
!          0.0345, 0.0387, 0.0432, 0.0483, 0.0540, 0.0606, &
!          0.0683, 0.0775, 0.0886, 0.1022, 0.1188, 0.1388, &
!          0.1625, 0.1905, 0.2229, 0.2602, 0.3025, 0.3500, &
!          0.4031, 0.4623, 0.5283, 0.6016, 0.6824, 0.7705, &
!          0.8650, 0.9646, 1.0676, 1.1713, 1.2722, 1.3662, &
!          1.4491, 1.5174, 1.5685, 1.6014, 1.6173, 1.6200, &
!          1.6150, 1.6082, 1.6040, 1.6058, 1.6157, 1.6353, &
!          1.6651, 1.7045, 1.7513, 1.8024, 1.8541, 1.9022, &
!          1.9429, 1.9738, 1.9947, 2.0072, 2.0132, 2.0141, &
!          2.0096, 1.9994, 1.9829, 1.9604, 1.9321, 1.8977, &
!          1.8562, 1.8052, 1.7419, 1.6646, 1.5738, 1.4734, &
!          1.3693, 1.2676, 1.1724, 1.0851, 1.0052, 0.9317, &
!          0.8637, 0.8004, 0.7414, 0.6862, 0.6348, 0.5871, &
!          0.5434, 0.5037, 0.4682, 0.4368, 0.4097, 0.3864, &
!          0.3667, 0.3499, 0.3355, 0.3231, 0.3123, 0.3029, &
!          0.2944, 0.2862, 0.2773, 0.2670, 0.2547, 0.2402, &
!          0.2238, 0.2061, 0.1882, 0.1712, 0.1562, 0.1434, &
!          0.1332, 0.1251, 0.1189, 0.1141, 0.1103, 0.1071, &
!          0.1043, 0.1018, 0.0996, 0.0979, 0.0968, 0.0964, &
!          0.0966, 0.0970, 0.0973, 0.0970, 0.0959, 0.0938, &
!          0.0909, 0.0873, 0.0831, 0.0784, 0.0732, 0.0676, &
!          0.0618, 0.0565, 0.0521, 0.0491, 0.0475, 0.0473, &
!          0.0480, 0.0492, 0.0504, 0.0514, 0.0519, 0.0521, &
!          0.0520, 0.0517, 0.0513, 0.0510, 0.0507, 0.0507, &
!          0.0508, 0.0512, 0.0515, 0.0518, 0.0519, 0.0518, &
!          0.0513, 0.0506, 0.0496, 0.0482, 0.0465, 0.0443, &
!          0.0418, 0.0387, 0.0351, 0.0310, 0.0263, 0.0214 /)
!
!!      Fixed normalization factors; a more accurate normalization would take
!!      in consideration longitude and time step
!!      ---------------------------------------------------------------------
!       real*8, save :: fBoreal = -1., fNonBoreal = -1
!       real,   save :: fDT=-1
!
!       integer :: hh, mm, ss, ndt, i, j, k
!       integer :: NN
!       real :: secs, secs_local, aBoreal, aNonBoreal, alpha
!
!!                              -----
!
!!      Normalization factor depends on timestep
!!      ----------------------------------------
!       if ( fDT /= cdt ) then
!            fBoreal = 0.0
!            fNonBoreal = 0.0
!            NN = 0
!            ndt = max(1,nint(cdt/DT))
!
!            do k = 1, N, ndt
!               NN = NN + 1
!               fBoreal    = fBoreal    + Boreal(k)
!               fNonBoreal = fNonBoreal + NonBoreal(k)
!            end do
!
!            fBoreal    = fBoreal / NN
!            fnonBoreal = fnonBoreal / NN
!            fDT = cdt ! so it recalculates only if necessary
!       end if
!
!
!!      Find number of secs since begining of the day (GMT)
!!      ---------------------------------------------------
!       hh = nhms/10000
!       mm = (nhms - 10000*hh) / 100
!       ss = nhms - 10000*hh - 100*mm
!       secs = 3600.*hh + 60.*mm + ss
!
!!      Apply factors depending on latitude
!!      -----------------------------------
!       do j = lbound(Ein,2), ubound(Ein,2)
!         do i = lbound(Ein,1), ubound(Ein,1)
!
!!            Find corresponding index in hardwired diurnal cycle
!!            240 = 24 * 60 * 60 secs / 360 deg
!!            ---------------------------------------------------
!             secs_local = secs + 240. * lons(i,j)
!             k = 1 + mod(nint(secs_local/DT),N)
!             if ( k < 1 ) k = N + k
!
!!            Apply diurnal cycle
!!            -------------------
!             aBoreal = Boreal(k) / fBoreal
!             aNonBoreal = NonBoreal(k) / fNonBoreal
!
!                if ( lats(i,j) >= 50. ) then
!                   Eout(i,j) = aBoreal    * Ein(i,j)
!                else if ( lats(i,j) >= 30. ) then
!                   alpha = (lats(i,j) - 30. ) / 20.
!                   Eout(i,j) = (1-alpha) * aNonBoreal * Ein(i,j) + &
!                                  alpha  * aBoreal    * Ein(i,j)
!                else
!                   Eout(i,j) = aNonBoreal * Ein(i,j)
!                end if
!          end do
!       end do
!
!     end subroutine Chem_BiomassDiurnal
!!==================================================================================
!
!   subroutine ReadPointEmissions( nymd, filename, nPts, vLat, vLon, vBase, vTop, vEmis, vStart, vEnd, unusable, label, rc)
!      integer, intent(in)            :: nymd
!      character(*), intent(in) :: filename
!      integer, intent(out)           :: nPts
!      real, allocatable, dimension(:), intent(out)    :: vLat, vLon, vTop, vBase, vEmis
!      integer, allocatable, dimension(:), intent(out) :: vStart, vEnd
!
!      type(KeywordEnforcer), optional, intent(in) :: unusable
!      character(*), optional, intent(in) :: label
!      integer, optional, intent(out) :: rc
!
!      ! Local arguments
!      type(EmissionReader) :: reader
!      character(:), allocatable :: label_
!      real, allocatable :: table(:,:)
!      integer :: nCols
!      integer :: status
!
!      if (present(label)) then
!         label_ = trim(label)
!      else
!         label_ = 'source'
!      end if
!
!      reader = EmissionReader()
!      call reader%open(filename, __RC__)
!      table = reader%read_table(label=label_, __RC__)
!      call reader%close(__RC__)
!
!      nCols = size(table,1)
!      nPts = size(table,2)
!      vStart = spread(-1, 1, nPts)
!      vEnd = spread(-1, 1, nPts)
!
!      vLat  = table(1,:)
!      vLon  = table(2,:)
!      vEmis = table(3,:)
!      vBase = table(4,:)
!      vTop  = table(5,:)
!      if (nCols >= 6) vStart = table(6,:)
!      if (nCols >= 7) vEnd = table(7,:)
!
!      where(vStart < 0) vStart = 000000
!      where(vEnd < 0)   vEnd   = 240000
!      call reader%close()
!
!      __RETURN__(__SUCCESS__)
!   end subroutine ReadPointEmissions
!
!!==================================================================================
!   subroutine open(this, filename, rc)
!      class(EmissionReader), intent(inout) :: this
!      character(*), intent(in) :: filename
!      integer, optional, intent(out) :: rc
!
!      integer :: status
!
!      __ASSERT__(.not. allocated(this%unit))
!      allocate(this%unit)
!
!      open(newunit=this%unit, file=filename,  &
!           form='formatted', access = 'sequential', status='old', &
!           action='read', __IOSTAT__)
!
!      __RETURN__(__SUCCESS__)
!   end subroutine open
!
!
!   subroutine close(this, rc)
!      class(EmissionReader), intent(inout) :: this
!      integer, optional, intent(out) :: rc
!
!      integer :: status
!
!      __ASSERT__(allocated(this%unit))
!      close(this%unit, __IOSTAT__)
!      deallocate(this%unit)
!
!      __RETURN__(__SUCCESS__)
!   end subroutine close
!
!
!   subroutine rewind_reader(this, rc)
!      class(EmissionReader), intent(in) :: this
!      integer, optional, intent(out) :: rc
!
!      integer :: status
!
!      __ASSERT__(allocated(this%unit))
!      rewind(this%unit, __IOSTAT__)
!
!      __RETURN__(__SUCCESS__)
!   end subroutine rewind_reader
!
!   function get_dims(this, label, rc) result(dims)
!      integer :: dims(2)
!      class(EmissionReader), intent(in) :: this
!      character(*), intent(in) :: label
!      integer, optional, intent(out) :: rc
!
!      integer :: status
!      logical :: eof
!      character(:), allocatable :: line
!      integer :: n_words
!
!      call this%rewind(__RC__)
!      call this%scan_to_label(label, __RC__)
!!      print*,__FILE__,__LINE__, ' found label'
!
!      dims = 0
!      do
!         line = this%next_line(eof=eof, __RC__)
!         __ASSERT__(.not. eof)
!         if (this%is_end_marker(line)) exit
!
!         dims(2) = dims(2) + 1
!
!         n_words = this%count_words(line)
!         dims(1) = max(dims(1), n_words)
!      end do
!
!      __RETURN__(__SUCCESS__)
!   end function get_dims
!
!   integer function count_words(this, line) result(n_words)
!      class(EmissionReader), intent(in) :: this
!      character(*), intent(in) :: line
!
!      integer :: idx, i0
!
!      n_words = 0
!      i0 = 0
!      do
!         ! scan to start of next word
!         idx = verify(line(i0+1:), ' ')
!
!         n_words = n_words + 1
!         i0 = i0 + idx
!
!         ! scan to end of current word
!         idx = index(line(i0+1:), ' ')
!         i0 = i0 + idx
!         if (idx == 0) exit
!
!      end do
!
!      return
!   end function count_words
!
!   logical function is_end_marker(this, line)
!      class(EmissionReader), intent(in) :: this
!      character(*), intent(in) :: line
!
!      is_end_marker = (line == '::')
!
!   end function is_end_marker
!
!   function read_table(this, label, rc) result(table)
!      class(EmissionReader), intent(in) :: this
!      real, allocatable :: table(:,:)
!      character(*), intent(in) :: label
!      integer, optional, intent(out) :: rc
!
!      integer :: i, j
!      integer :: dims(2)
!      integer :: status
!      logical :: eof
!      character(:), allocatable :: line
!
!      dims = this%get_dims(label, __RC__)
!      call this%scan_to_label(label, __RC__)
!
!      associate (n_words => dims(1), n_lines => dims(2))
!        allocate(table(n_words, n_lines), __STAT__)
!
!        do j = 1, n_lines
!           line = this%next_line(eof=eof)
!           __ASSERT__(.not. eof)
!
!           read(line,*, iostat=status) (table(i,j),i=1,n_words)
!           __VERIFY__(status)
!        end do
!
!      end associate
!
!      __RETURN__(__SUCCESS__)
!   end function read_table
!
!   function next_line(this, eof, rc) result(line)
!      character(:), allocatable :: line
!      class(EmissionReader), intent(in) :: this
!      logical, intent(out) :: eof
!      integer, optional, intent(out) :: rc
!
!      integer, parameter :: MAX_LINE_LEN=1024
!      character(len=MAX_LINE_LEN) :: buffer
!      integer :: idx
!      integer :: status
!
!      eof = .false.
!      do
!
!         read(this%unit,'(a)', iostat=status) buffer
!         if (status == IOSTAT_END) then
!            eof = .true.
!            __RETURN__(__SUCCESS__)
!         end if
!         __VERIFY__(status)
!
!         idx = index(buffer, '#')
!         if (idx == 0) idx = len(buffer)
!
!         line = trim(buffer(:idx-1))
!         if (line /= '')  exit
!
!      end do
!
!      __RETURN__(__SUCCESS__)
!   end function next_line
!
!   subroutine scan_to_label(this, label, rc)
!      class(EmissionReader), intent(in) :: this
!      character(*), intent(in) :: label
!      integer, optional, intent(out) :: rc
!
!      integer :: status
!      logical :: eof
!      character(:), allocatable :: line
!
!      call this%rewind(__RC__)
!      do
!         line = this%next_line(eof=eof, __RC__)
!         if (line == label // '::') exit
!      end do
!
!      __RETURN__(__SUCCESS__)
!   end subroutine scan_to_label
 end module QuickChem_Process
