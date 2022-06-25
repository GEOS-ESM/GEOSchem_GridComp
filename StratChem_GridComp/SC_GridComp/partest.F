      SUBROUTINE partest(k,k1Strat,levels,daytime,m,aj,ak,broe,brxe,ch2o,ch3o2,ch4,cloe,clono2e,
     &           h,hbre,hno3e,ho2,hobre,ne,n2o5e,noe,no2e,no3e,noxe,o3e,o3pe,oh,rbrnbrx,rbrobrx,
     &                            rcl2o2clo,rclclo,rcloclx,rnno,rno2nox,rno3no2,rnono2,rocloclo,
     &                                                          levels_cal,numreacs,numphoto,o2)

      IMPLICIT NONE
      INTEGER, PARAMETER :: DBL = KIND(0.00D+00)

      INTEGER, INTENT(IN) :: k,k1Strat,levels,levels_cal,numreacs,numphoto
      LOGICAL, INTENT(IN) :: daytime
      REAL, INTENT(IN) :: m(levels),o2

      REAL(KIND=DBL), INTENT(IN) :: aj(numphoto)
      REAL(KIND=DBL), INTENT(IN) :: broe,brxe,ch2o,ch3o2,ch4,cloe,clono2e
      REAL(KIND=DBL), INTENT(IN) :: hbre,hno3e,ho2,hobre,ne,noe,h,n2o5e
      REAL(KIND=DBL), INTENT(IN) :: noxe,o3e,o3pe,oh

      REAL(KIND=DBL), INTENT(INOUT) :: rbrnbrx,rbrobrx,rcl2o2clo,rclclo,rcloclx
      REAL(KIND=DBL), INTENT(INOUT) :: rnno,rno2nox,rno3no2,rnono2,rocloclo
      REAL(KIND=DBL), INTENT(INOUT) :: no2e,no3e,ak(levels,numreacs)

      INTEGER :: iter
      REAL :: o2m
      REAL(KIND=DBL) :: asdfgh,asdfgk,asdfgj,asdfgl
      REAL(KIND=DBL) :: brotmp,rbrbro,rbrnbro,temp

c     define rnno=[n]/[no], rnono2=[no]/[no2], rno3no2=[no3]/[no2], and
c     rno2nox = [no2]/[nox] assuming p.e. 
c     assume [o(3p)] at dawn equals [o(3p)] at dusk. 
c     assume [ho2]   at dawn equals [ho2]   at dusk. 
c     assume [n]     at dawn equals [n]     at dusk. 
c     note: [no]=[no2]*rnono2 
      
c     define rclclo=[cl]/[clo] and rcloclx=[clo]/[clx] 
      
! The following comments were in the GEOS-4 version of the code:
! 
! "aj(17)*no3/no2 has been removed from the numerator and ak(k,74)*no3 has 
! been removed from the denominator of rnono2 to avoid night/day transition 
! problems."
!
! However, these components were actually still included in the code!  When
! implemented in GEOS-5, problems arose at sunset in the upper troposphere
! as no2e and no3e sometimes went to zero.  So, we will set them to a 
! minimum value.  JEN, 25 Sep 2006.
      
      if(k .ge. k1Strat .and. k .le. levels_cal) then

               if(daytime) then
                  
                  o2m = o2*m(k) 
                  
c     nox partitioning

                  IF(no2e < 1.00E-25) no2e = 1.00E-25
                  IF(no3e < 1.00E-25) no3e = 1.00E-25
                  
                  rnono2 = (aj(7) + ak(k,43)*o3pe
     &                 +ak(k,125)*h
     &                 +aj(17)*no3e/no2e
     &                 +ak(k,20)*ne*o2m/no2e)
     &                 /(aj(16)+aj(39)+ak(k,9)*o3e
     &                 +ak(k,38)*ho2
     &                 +ak(k,28)*cloe
     &                 +ak(k,15)*ch3o2
     &                 +ak(k,66)*broe
     &                 +ak(k,74)*no3e
     &                 +ak(k,77)*o3pe*m(k)+ak(k,47)*ne)
                  
                  rnno = (aj(16)+aj(39))/
     &                 (ak(k,47)*rnono2*no2e+ak(k,81)*no2e
     &                 +ak(k,20)*o2m)
                  
                  rno3no2 = (ak(k,10)*o3e
     &                 +ak(k,44)*o3pe*m(k)
     &                 +(aj(15)+ak(k,13)*oh
     &                 +ak(k,35)*o3pe)*clono2e/no2e
     &                 +(ak(k,31)*m(k)+aj(9))
     &                 *n2o5e/no2e
     &                 +ak(k,37)*oh*hno3e/no2e)
     &                 /(aj(5)+aj(17)
     &                 +ak(k,74)*rnono2*no2e
     &                 +ak(k,46)*no2e*m(k))
                  
                  rno2nox = 1./
     &                 (1.+rnno+rnono2+rno3no2)
                  
c...  clo partitioning
                  
                  rclclo = (ak(k,25)*o3pe
     &                 +ak(k,28)*noe
     &                 +ak(k,60)*oh
     &                 +ak(k, 6)*ho2)
     &                 /(ak(k,24)*o3e
     &                 +ak(k,67)*ho2)
                  
                  rcl2o2clo = (ak(k,107)*cloe*m(k))
     &                 /(ak(k,108)*m(k)+(aj(43)))
                  
c...  corrected partitioning error 9/28/94.
c     rcloclx = 1./(rclclo+rcl2o2clo+1.) 
                  rocloclo = ak(k,102)*broe/aj(37)
                  
                  rcloclx = 1./
     &                 (rclclo+2.*rcl2o2clo+rocloclo+1.) 
                  
c     .sds fix to ak(26) so o3(small) doesn't cause negative hcl and source of cl
                  temp = ch4*rclclo*rcloclx
                  if(ak(k,26)*temp.gt.5.e-04) ak(k,26) = 5.e-04/temp
c     .sds end fix
                  
c...  brx partitioning, steady state ratio
                  rbrnbro = 
     &                 ak(k,96)*m(k)*rno2nox*noxe
     &                 /(aj(28)+ak(k,115))
                  rbrnbro = max(rbrnbro,0.)
                  
                  brotmp = broe
                  asdfgh = ak(k,69)*ch2o
     &                 + ak(k,91)*o3e
     &                 + ak(k,92)*ho2
                  
                  asdfgl = aj(27) + cloe
     &                 *(ak(k,102)+ak(k,103)+ak(k,104))
     &                 +ak(k,118)*o3pe
     &                 +ak(k,66)*noe
                  
                  do iter=1,3
                     asdfgj = 2.*ak(k,94)*brotmp + asdfgl
                     if (brotmp.gt.1.) then
                        asdfgk =  (aj(25)*hobre
     &                     +ak(k,95)*oh*hbre)/brotmp
                     else 
                        asdfgk = 0
                     endif
                     rbrbro = ( asdfgj + asdfgk ) / asdfgh
                     rbrbro = max(0.,rbrbro)
                     rbrobrx = 1./(1.+rbrbro+rbrnbro)
                     brotmp = rbrobrx*brxe
                  enddo
                  rbrnbrx = rbrobrx*rbrnbro
               endif

      endif

      RETURN
      END SUBROUTINE partest
