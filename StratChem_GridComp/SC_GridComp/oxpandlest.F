      SUBROUTINE oxpandlest(k,k1Strat,levels,tropo,daytime,m,aj,ak,brxe,ch2o,ch3,ch3ooh,ch4,
     &   cl2o2e,clono2e,clxe,h,h2o,h2o2e,hbre,hcle,hno3e,ho2,ho2no2e,hocle,lh2o2,lo3,lo3brx,
     &             lo3chx,lo3cl,lo3cly,lo3hox,lo3nox,lo3oh,lo3oxsq,n2o,n2o5e,noxe,oh,ox,oxe,
     &                    ph2o2,po3,rain,rbrobrx,rclclo,rcloclx,rnno,rno2nox,rno3no2,rnono2,
     &                                  ro1do3,ro3ox,ro3po3,levels_cal,numreacs,numphoto,o2)

      IMPLICIT NONE
      INTEGER, PARAMETER :: DBL = KIND(0.00D+00)

      INTEGER, INTENT(IN) :: k,k1Strat,levels,levels_cal,numreacs,numphoto
      LOGICAL, INTENT(IN) :: tropo(levels),daytime
      REAL, INTENT(IN) :: m(levels),rain(levels),o2

      REAL(KIND=DBL), INTENT(IN) :: aj(numphoto),ak(levels,numreacs)
      REAL(KIND=DBL), INTENT(IN) :: brxe,ch2o,ch3,ch3ooh,ch4,cl2o2e,clono2e
      REAL(KIND=DBL), INTENT(IN) :: clxe,h,h2o,h2o2e,hbre,hcle,hno3e,ho2
      REAL(KIND=DBL), INTENT(IN) :: ho2no2e,hocle,n2o,n2o5e,noxe,oh,ox,oxe
      REAL(KIND=DBL), INTENT(IN) :: rbrobrx,rclclo,rcloclx,rnno,rno2nox,rno3no2
      REAL(KIND=DBL), INTENT(IN) :: rnono2,ro1do3,ro3ox,ro3po3

      REAL(KIND=DBL), INTENT(INOUT) :: lh2o2,lo3,lo3brx,lo3chx,lo3cl,lo3cly
      REAL(KIND=DBL), INTENT(INOUT) :: lo3hox,lo3nox,lo3oh,lo3oxsq,ph2o2,po3

      REAL :: o2m
      REAL(KIND=DBL) :: p1,p2,po3hox,po3phot,po3ox,po3cl,po3noy,po3het
      REAL(KIND=DBL) :: po3chx,tbroe,tcle,tcloe,tne,tnoe,tno2e,tno3e
      REAL(KIND=DBL) :: to1de,to3e,to3pe
      
c... define  =[o]/[o3], ro3ox=[o3]/[ox], 
c...  and ro1do3=[o(1d)]/[o3]   
c...  these don't need to be changed as they don't involve estimates!

      if(k .le. levels_cal) then

c... in the troposphere now relaxing to climatology
          if(tropo(k)) then
            lo3 = 0.00
            if(daytime .and. k .ge. k1Strat) then !Selecting upper trop

c... Solve for the H2O2 concentration 
c... Calculate average production and loss for predictor corrector
              ph2o2 = .5*( ph2o2
     &         +ak(k,12)*m(k)*oh**2
     &         +( ak(k,19)+ak(k,64)*m(k) )
     &         *ho2**2 )

              to3pe = ro3po3*ro3ox*oxe
              tcle = rclclo*rcloclx*clxe

              lh2o2 = .5*( lh2o2
     &         +( ak(k,29)*oh
     &         +ak(k,7)*tcle
     &         +ak(k,33)*to3pe
     &         +aj(8)+rain(k) )*h2o2e )
             endif
c... in the stratosphere prod/loss
           else
             if(daytime .and. k .ge. k1Strat) then !Strat where k>=i1chem
               o2m = o2*m(k)
               to3e = ro3ox*oxe
               to3pe = ro3po3*to3e
               to1de = ro1do3*to3e
               tcloe = rcloclx*clxe
               tcle = rclclo*tcloe
               tno2e = rno2nox*noxe
               tnoe = rnono2*tno2e
               tno3e = rno3no2*tno2e
               tbroe = rbrobrx*brxe

               po3hox = 2.*(ak(k,19)+ak(k,64)*m(k))
     &          *ho2**2 + (ak(k,72)   
     &          +2.*ak(k,73))*h*ho2

               po3phot =  2.0*(aj(1)+aj(41))*o2m
     &          +aj(4)*h2o
     &          +aj(14)*n2o 
     &          +(aj(38)+aj(40))*h2o

               po3ox = to1de* 
     &          (ak(k,39)*h2o+ak(k,45)*n2o) 

               po3cl = 2.*ak(k,67)*ho2*tcle

               po3noy = 2.*( ak(k,20)*o2m*rnno 
     &          +ak(k,38)*ho2 )*tnoe

c... average the production and loss calculated before and after the
c... estimates!

               po3het = ak(k,112)*clono2e
     &          +ak(k,110)*n2o5e

               po3chx = 2.*ak(k,50)*ch3*o2m*m(k)

               po3 = .5*( po3+po3hox+po3phot+po3ox+
     &          po3cl+po3noy+po3het+po3chx )

c... calc loss
               lo3hox = 2.*(ak(k,4)*to3e*oh
     &           +ak(k,33)*to3pe*h2o2e
     &           +ak(k,41)*to3pe*oh)

                lo3oh = (ak(k,14)*ch4
     &           +3.*ak(k,29)*h2o2e
     &           +ak(k,37)*hno3e
     &           +ak(k,40)*ho2
     &           +ak(k,55)*ho2no2e
     &           +ak(k,59)*oh
     &           +ak(k,65)*ho2*m(k) )*oh

                lo3nox = 2.*(ak(k,43)*to3pe*tno2e
     &           +aj(17)*tno3e)+ak(k,81)*rnno*tnoe*tno2e

                lo3oxsq = 2.*(ak(k,2)*to3e
     &           +ak(k,90)*m(k)*to3pe )*to3pe
     &           +2.*ak(k,123)*to3e*to1de


                lo3cl = 2.
     &           *( ak(k,7)*tcle*h2o2e
     &           +ak(k,25)*to3pe*tcloe
     &           +(ak(k,60)+ak(k,61))*oh*tcloe )

                lo3cly = ( ak(k,27)*hcle
     &           +ak(k,62)*hocle )*oh
     &           +2.*aj(36)*cl2o2e

                lo3brx = 2.*((ak(k,103)+ak(k,104))*tcloe
     &           +ak(k,94)*tbroe+ak(k,118)*to3pe)*tbroe
     &           +ak(k,95)*hbre*oh 

                lo3chx= ( ak(k,51)*ch2o
     &           +ak(k,58)*ch3ooh )*oh

                lo3 = .5*( lo3 +lo3hox+lo3oh+
     &            lo3nox+lo3oxsq+lo3cl+
     &            lo3cly+lo3brx+lo3chx )

c... Solve for the H2O2 concentration 
c... Calculate average production and loss for predictor corrector
                ph2o2 = .5*( ph2o2
     &           +ak(k,12)*m(k)*oh**2
     &           +( ak(k,19)+ak(k,64)*m(k) )
     &           *ho2**2 )

                lh2o2 = .5*( lh2o2
     &           +( ak(k,29)*oh
     &           +ak(k,7)*tcle
     &           +ak(k,33)*to3pe
     &           +aj(8)+rain(k) )*h2o2e )
             else
                lo3 = 0.0
                po3 = 0.0
              endif
           endif  
         
       endif

      RETURN
      END SUBROUTINE oxpandlest
