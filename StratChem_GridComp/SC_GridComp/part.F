      SUBROUTINE part(k,k1Strat,levels,daytime,m,aj,ak,bro,brx,ch2o,ch3o2,ch4,clo,clono2,
     &               h,hbr,hno3,ho2,hobr,n,n2o5,no,no2,no3,nox,o3,o3p,oh,rbrnbrx,rbrobrx,
     &                     rcl2o2clo,rclclo,rcloclx,rnno,rno2nox,rno3no2,rnono2,rocloclo,
     &                                                   levels_cal,numreacs,numphoto,o2)

      IMPLICIT NONE
      INTEGER, PARAMETER :: DBL = KIND(0.00D+00)

      INTEGER, INTENT(IN) :: k,k1Strat,levels,levels_cal,numreacs,numphoto
      LOGICAL, INTENT(IN) :: daytime
      REAL, INTENT(IN) :: m(levels),o2
      
      REAL(KIND=DBL), INTENT(IN) :: aj(numphoto)
      REAL(KIND=DBL), INTENT(IN) :: bro,brx,ch2o,ch3o2,ch4,clo,clono2
      REAL(KIND=DBL), INTENT(IN) :: hbr,hno3,ho2,hobr,n,no,h,n2o5
      REAL(KIND=DBL), INTENT(IN) :: no2,no3,nox,o3,o3p,oh

      REAL(KIND=DBL), INTENT(INOUT) :: rbrnbrx,rbrobrx,rcl2o2clo,rclclo,rcloclx
      REAL(KIND=DBL), INTENT(INOUT) :: rnno,rno2nox,rno3no2,rnono2,rocloclo
      REAL(KIND=DBL), INTENT(INOUT) :: ak(levels,numreacs)

crk 941108 added hobr,  941111 added brono2+h2o
crk 951106 updated to modified br scheme and completed terms in no/no2 
c        partitioning
 
      INTEGER :: iter
      REAL :: o2m
      REAL(KIND=DBL) :: asdfgh,asdfgk,asdfgj,asdfgl
      REAL(KIND=DBL) :: brotmp,rbrbro,rbrnbro,temp

c... define rnno=[n]/[no], rnono2=[no]/[no2], rno3no2=[no3]/[no2], and
c... rno2nox = [no2]/[nox] assuming p.e. 
c... rbrnbrx= [brono2]/brx
c... assume [o(3p)] at dawn equals [o(3p)] at dusk. 
c... assume [ho2]   at dawn equals [ho2]   at dusk. 
c... assume [n]     at dawn equals [n]     at dusk. 
c... note: [no]=[no2]*rnono2 

c... define rclclo=[cl]/[clo] and rcloclx=[clo]/[clx] 
c... define rbrobrx=bro/(bro+br), brcl(day) will be added to brx
c... and set to 0. in solver.

! The following comments were in the GEOS-4 version of the code:
! 
! "aj(17)*no3/no2 has been removed from the numerator and ak(k,74)*no3 has 
! been removed from the denominator of rnono2 to avoid night/day transition 
! problems."
!
! However, these components were actually still included in the code!  It
! works if no2 and no3 are set to some mimimum value before enering this
! routine.  See comments in partest.  JEN, 25 Sep 2006.

      if(k .ge. k1Strat .and. k .le. levels_cal) then
               if (daytime) then
                  o2m = o2*m(k)

c...  Nitrogen species
                  
                  rnono2 = (aj(7)
     &                 +ak(k,43)*o3p
     &                 +ak(k,125)*h
     &                 +aj(17)*no3/no2
     &                 +ak(k,20)*n*o2m/no2)
     &                 /(ak(k,9)*o3
     &                 +ak(k,38)*ho2
     &                 +aj(16)+aj(39)
     &                 +ak(k,28)*clo
     &                 +ak(k,15)*ch3o2
     &                 +ak(k,66)*bro
     &                 +ak(k,74)*no3
     &                 +ak(k,77)*o3p*m(k)+ak(k,47)*n)
                  
                  rnno = (aj(16)+aj(39))
     &                 /(ak(k,47)*rnono2*no2+ak(k,81)*no2
     &                 +ak(k,20)*o2m)
                  
                  rno3no2= (ak(k,10)*o3 
     &                 +ak(k,44)*o3p*m(k)
     &                 +(aj(15)+ak(k,13)*oh
     &                 +ak(k,35)*o3p)*clono2/no2
     &                 +(ak(k,31)*m(k)+aj(9))
     &                 *n2o5/no2
     &                 +ak(k,37)*oh*hno3/no2)
     &                 /(aj(5) + aj(17)
     &                 +ak(k,46)*no2*m(k)
     &                 +ak(k,74)*rnono2*no2)
                  
                  rno2nox = 1./
     &                 (1.+rnno+rnono2+rno3no2)
                  
c...  Chlorine species    
                  rclclo = (ak(k,25)*o3p
     &                 +ak(k,28)*no
     &                 +ak(k,60)*oh
     &                 +ak(k,6)*ho2)
     &                 /(ak(k,24)*o3
     &                 +ak(k,67)*ho2)
                  
                  rcl2o2clo = 
     &                 (ak(k,107)*clo*m(k))
     &                 /(ak(k,108)*m(k)+(aj(43)))
                  
c...  corrected partitioning error (9/28/94)
c     rcloclx = 1./(rclclo+rcl2o2clo+1.) 
                  rocloclo = 
     &                 ak(k,102)*bro/aj(37)
                  
                  rcloclx = 1./
     &                 (rclclo+2.*rcl2o2clo+rocloclo+1.) 
                  
c     .sds fix to ak(26) so o3(small) doesn't cause negative hcl and source of cl
                  temp = ch4*rclclo*rcloclx
                  if(ak(k,26)*temp.gt.5.e-4) ak(k,26) = 5.e-4/temp
c     .sds end fix
                  
c...  Bromine species partitioning (steady state ratio)
                  rbrnbro = ak(k,96)*m(k)*rno2nox
     &                 *nox/(aj(28)+ak(k,115))
                  rbrnbro = max(rbrnbro,0.)
                  
                  brotmp = bro
                  asdfgh = ak(k,69)*ch2o
     &                 + ak(k,91)*o3
     &                 + ak(k,92)*ho2
                  
                  asdfgl = aj(27) + clo
     &                 *(ak(k,102)+ak(k,103)+ak(k,104))
     &                 +ak(k,118)*o3p
     &                 +ak(k,118)*o3p
     &                 +ak(k,66)*no

                  do iter=1,3
                     asdfgj = 2.*ak(k,94)*brotmp + asdfgl
                     if (brotmp.gt.1.) then
                        asdfgk =  (aj(25)*hobr
     &                      +ak(k,95)*oh*hbr)/brotmp
                     else 
                        asdfgk = 0
                     endif
                     rbrbro = ( asdfgj + asdfgk ) / asdfgh
                     rbrbro = max(0.,rbrbro)
                     rbrobrx = 1./(1.+rbrbro+rbrnbro)    
                     brotmp = rbrobrx*brx
                  enddo
                  rbrnbrx = rbrobrx*rbrnbro
               endif
      endif
c     
      RETURN
      END SUBROUTINE part
