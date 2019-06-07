!
!  Main header file required by CTM - provide just enough stuff here
!  to allow model subroutines to be run as a stand-alone code.
!
!                                                 Oliver (9 July 99)
!-----------------------------------------------------------------------
!
      integer  ipar, jpar, lpar
      integer lpar_max, jpnl_max, jppj_max         !Added for GMI {PJC}
      integer  jpnl, jppj
      logical  ldeg45
      parameter(ipar=1,jpar=1)   !  Number of lat/lon points in CTM.
!      parameter(lpar=9,jpnl=8)   !  Number of levels in CTM and number.
!c      parameter(lpar=23,jpnl=14) !    of levels requiring chemistry.
      parameter(lpar_max = 70)  !  Max number of levels in CTM.
      parameter(jpnl_max = lpar_max) ! Max num. of levels doing chemistry.
!      parameter(jppj=3)          !  Number of photolytic reactions supplied.
      parameter(jppj_max = 55)  !  Max number of photolysis reactions.
      parameter(ldeg45=.false.)  !  Logical flag for degraded CTM resolution.

!
      real*8  xgrd(ipar)         !  Longitude (midpoint, radians)
      real*8  ygrd(jpar)         !  Latitude  (midpoint, radians)
      real*8  ydgrd(jpar)        !  Latitude  (midpoint, degrees)
!      real*8  etaa(lpar+1)       !  Eta(a) value for level boundaries
!      real*8  etab(lpar+1)       !  Eta(b) value for level boundaries
      real*8  etaa(lpar_max+1)   !  Eta(a) value for level boundaries
      real*8  etab(lpar_max+1)   !  Eta(b) value for level boundaries
!
      real*8  tau                !  Time of Day (hours, GMT)
      integer month              !  Number of month (1-12)
      integer iday               !  Day of year
!
      common/dims/lpar, jpnl, jppj     ! Added for GMI. {PJC}
      common/grid/xgrd,ygrd,ydgrd,etaa,etab
      common/fastj_time/tau,month,iday
!
!-----------------------------------------------------------------------
