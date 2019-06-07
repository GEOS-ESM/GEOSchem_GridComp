!
!  Main header file required by CTM - provide just enough stuff here
!  to allow model subroutines to be run as a stand-alone code.
!
!                                                 Oliver (9 July 99)
!-----------------------------------------------------------------------
!
      integer  ipar, jpar, lpar
      integer lpar_max, jpnl_max, jppj_max         !Added for GMI {AOO, 8/04}
      integer  jpnl, jppj
      logical  ldeg45
      parameter(ipar=1,jpar=1)   !  Number of lat/lon points in CTM

!cc Use lpar=9 for GISS 9-layer CTM
!cc      parameter(lpar=9,jpnl=9)   !  Number of levels in CTM and number
!cc      parameter(lpar=9,jpnl=8)   !  Number of levels in CTM and number

!cc Use lpar=23 for the GISS 23-layer CTM
!cc      parameter(lpar=23,jpnl=23) !    of levels requiring chemistry

!cc Use lpar=40 to do US std atmosphere
      parameter(lpar_max=72) !    of levels requiring chemistry
      parameter(jpnl_max = lpar_max)

!ccc      parameter(jppj=46)         !  Number of photolytic reactions supplied
      parameter(jppj_max=82)         !  Number of photolytic reactions supplied


      parameter(ldeg45=.false.)  !  Logical flag for degraded CTM resolution
!
      real*8  xgrd(ipar)         !  Longitude (midpoint, radians)
      real*8  ygrd(jpar)         !  Latitude  (midpoint, radians)
      real*8  ydgrd(jpar)        !  Latitude  (midpoint, degrees)
!      real*8  etaa(lpar+1)       !  Eta(a) value for level boundaries
!      real*8  etab(lpar+1)       !  Eta(b) value for level boundaries
      real*8  etaa(lpar_max+1)       !  Eta(a) value for level boundaries
      real*8  etab(lpar_max+1)       !  Eta(b) value for level boundaries
!
      real*8  tau                !  Time of Day (hours, GMT)
      integer month              !  Number of month (1-12)
      integer iday               !  Day of year

      logical do_model_clima     ! determines if climatology data come
                                 ! from the model.

!      real*8  DUSSTD,OUSSTD,TUSSTD,PUSSTD

!
      common/dims/lpar, jpnl, jppj     ! Added for GMI. {AOO, 8/04}
      common/grid/xgrd,ygrd,ydgrd,etaa,etab
      common/fast_JX_time/tau,month,iday
      common/fast_JX_clima/do_model_clima

!      common/usstd/ DUSSTD(52),OUSSTD(52),TUSSTD(52),PUSSTD(52)

!
!-----------------------------------------------------------------------
