      SUBROUTINE oxpandl(k,k1Strat,levels,tropo,daytime,m,aj,ak,brx,ch2o,ch3,ch3ooh,ch4,
     &        cl2o2,clono2,clx,h,h2o,h2o2,hbr,hcl,hno3,ho2,ho2no2,hocl,lh2o2,lo3,lo3brx,
     &        lo3chx,lo3cl,lo3cly,lo3hox,lo3nox,lo3oh,lo3oxsq,n2o,n2o5,nox,o1d,o3,oh,ox,
     &                ph2o2,po3,rain,rbrobrx,rclclo,rcloclx,rnno,rno2nox,rno3no2,rnono2,
     &                              ro1do3,ro3ox,ro3po3,levels_cal,numreacs,numphoto,o2)

      IMPLICIT NONE
      INTEGER, PARAMETER :: DBL = KIND(0.00D+00)

      INTEGER, INTENT(IN) :: k,k1Strat,levels,levels_cal,numreacs,numphoto
      LOGICAL, INTENT(IN) :: tropo(levels),daytime
      REAL, INTENT(IN) :: m(levels),rain(levels),o2

      REAL(KIND=DBL), INTENT(IN) :: aj(numphoto),ak(levels,numreacs)
      REAL(KIND=DBL), INTENT(IN) :: brx,ch2o,ch3,ch3ooh,ch4,cl2o2,clono2
      REAL(KIND=DBL), INTENT(IN) :: clx,h,h2o,h2o2,hbr,hcl,hno3,ho2
      REAL(KIND=DBL), INTENT(IN) :: ho2no2,hocl,n2o,n2o5,nox,o1d,o3,oh,ox
      REAL(KIND=DBL), INTENT(IN) :: rbrobrx,rclclo,rcloclx,rnno
      REAL(KIND=DBL), INTENT(IN) :: rno2nox,rno3no2,rnono2

      REAL(KIND=DBL), INTENT(INOUT) :: lh2o2,lo3,lo3brx,lo3chx,lo3cl,lo3cly
      REAL(KIND=DBL), INTENT(INOUT) :: lo3hox,lo3nox,lo3oh,lo3oxsq,ph2o2,po3
      REAL(KIND=DBL), INTENT(INOUT) :: ro1do3,ro3ox,ro3po3

c...  Latest update 27 October 1994 
C     rk  951107 Updated to include CHx chem and new Br scheme
C     sds 980317 Updated to include tropospheric loss

      REAL :: o2m
      REAL(KIND=DBL) :: po3hox,po3phot,po3ox,po3cl,po3noy,po3het
      REAL(KIND=DBL) :: po3chx,tbro,tcl,tclo,tn,tno,tno2,tno3
      REAL(KIND=DBL) :: to1d,to3,to3p
      
c...  define ro3po3=[o]/[o3], ro3ox=[o3]/[ox], and ro1do3=[o(1d)]/[o3]   
      
c...  Daytime Production and Loss for Ox

      if(k .le. levels_cal) then
c...  in the troposphere, now relaxing to climatology
          if(tropo(k)) then
            lo3 = 0.00

             if(daytime .and. k .ge. k1Strat) then !Selecting upper trop
              o2m = o2*m(k)

              ro3po3 = (aj(2)+aj(3)
     &         +(ak(k,124)*o1d))
     &         /(ak(k,1)*o2m*m(k))

	      ro3ox = 1. / (ro3po3 + 1. )

              ro1do3 =  aj(2)/(ak(k,8)*m(k)+ak(k,124)*o3)

	      ph2o2 = ak(k,12)*m(k)*oh**2
     &         +(ak(k,19)+ak(k,64)*m(k))*ho2**2

	      to3p = ro3po3*ro3ox*ox
	      tcl = rclclo*rcloclx*clx

	      lh2o2 = ( ak(k,29)*oh
     &         +ak(k,7)*tcl
     &         +ak(k,33)*to3p
     &         +aj(8) +rain(k) )*h2o2
	     endif
c... in the stratosphere prod/loss
           else
            if(daytime .and. k .ge. k1Strat) then !Above i1chem and not trop

               o2m = o2*m(k)

               ro3po3 = (aj(2)+aj(3)
     &          +(ak(k,124)*o1d))
     &          /(ak(k,1)*o2m*m(k)) 

               ro3ox = 1. / (ro3po3 + 1. )

               ro1do3 =  aj(2)/(ak(k,8)*m(k)+ak(k,124)*o3)

               to3 = ro3ox*ox
               to3p = ro3po3*to3
               to1d = ro1do3*to3
               tno2 = rno2nox*nox
               tno3 = rno3no2*tno2
               tno = rnono2*tno2
               tn = rnno*tno
               tclo = rcloclx*clx
               tcl = rclclo*tclo
               tbro = rbrobrx*brx

               po3hox = ( 2.
     &          *( ak(k,19)+ak(k,64)*m(k) )*ho2
     &          +( ak(k,72)+2.*ak(k,73) )*h )
     &          *ho2

               po3phot =  2.0*(aj(1)+aj(41))*o2m
     &          +aj(4)*h2o
     &          +aj(14)*n2o 
     &          +(aj(38)+aj(40))*h2o

               po3ox = ( ak(k,39)*h2o
     &          +ak(k,45)*n2o )*to1d

               po3cl = 2.*ak(k,67)*ho2*tcl
      
               po3noy = 2.*
     &          (ak(k,20)*o2m*tn
     &          +ak(k,38)*ho2*tno) 

               po3het = ak(k,112)*clono2
     &          +ak(k,110)*n2o5

               po3chx = 2.*ak(k,50)*ch3*o2m*m(k)

               po3 = po3hox+po3phot+po3ox+po3cl+po3noy+po3het
     &          +po3chx

               lo3hox = 2.*(ak(k,4)*oh*to3
     &          +( ak(k,33)*h2o2
     &          +ak(k,41)*oh )*to3p)

               lo3nox = 2.*(ak(k,43)*to3p*tno2
     &          +aj(17)*tno3)+ak(k,81)*rnno*tno*tno2

               lo3oxsq = 2.*(ak(k,2)*to3*to3p
     &          +ak(k,90)*m(k)*to3p*to3p)
     &          +2.*ak(k,123)*to3*to1d

               lo3cl = 2.*(ak(k,7)*h2o2*tcl
     &          +ak(k,25)*to3p*tclo
     &          +(ak(k,60)+ak(k,61))*oh*tclo)
        
               lo3cly = 2.*aj(36)*cl2o2
     &          +(ak(k,27)*hcl
     &          +ak(k,62)*hocl )*oh

               lo3brx = 2.*((ak(k,103)+ak(k,104))*tclo
     &          +ak(k,94)*tbro+ak(k,118)*to3p)*tbro
     &          +ak(k,95)*hbr*oh 

               lo3oh = ( ak(k,14)*ch4
     &          +3.*ak(k,29)*h2o2
     &          +ak(k,37)*hno3
     &          +ak(k,40)*ho2
     &          +ak(k,55)*ho2no2
     &          +ak(k,59)*oh
     &          +ak(k,65)*ho2*m(k) )*oh
 
               lo3chx= ( ak(k,51)*ch2o
     &          +ak(k,58)*ch3ooh )*oh

               lo3 = lo3hox+lo3nox+
     &           lo3oxsq+lo3cl+lo3cly+
     &           lo3oh+lo3brx+lo3chx
c
c... Daytime Production and Loss terms for H2O2
c
                ph2o2 = ak(k,12)*m(k)*oh**2
     &         +(ak(k,19)+ak(k,64)*m(k))*ho2**2

                lh2o2 = ( ak(k,29)*oh
     &           +ak(k,7)*tcl
     &           +ak(k,33)*to3p
     &           +aj(8) +rain(k) )*h2o2
             else
                lo3 = 0.0
                po3 = 0.0
             endif
           endif

       endif

       RETURN
       END SUBROUTINE oxpandl
