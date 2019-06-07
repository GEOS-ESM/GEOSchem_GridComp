      SUBROUTINE hoxpartest(k,k1Strat,levels,daytime,m,tropo,harvard_oh,ak,aj,brxe,ch2o,ch3o,
     &          ch3o2,ch3ooh,ch4,cho,clono2e,clxe,co,h,hbre,hcle,h2o,h2o2e,hno3e,ho2,ho2no2e,
     &            hobre,hocle,noxe,o1d,o3,o3p,oh,rbrnbrx,rbrobrx,rclclo,rcloclx,rnno,rno2nox,
     &                                         rnono2,levels_cal,niters,numreacs,numphoto,o2)

      IMPLICIT NONE
      INTEGER, PARAMETER :: DBL = KIND(0.00D+00)

      INTEGER, INTENT(IN) :: k,k1Strat,levels
      INTEGER, INTENT(IN) :: levels_cal,niters,numreacs,numphoto
      LOGICAL, INTENT(IN) :: daytime,tropo(levels)
      REAL, INTENT(IN) :: m(levels),harvard_oh(levels),o2

      REAL(KIND=DBL), INTENT(IN) :: aj(numphoto),ak(levels,numreacs)
      REAL(KIND=DBL), INTENT(IN) :: brxe,ch2o,ch3o,ch3o2,ch3ooh,ch4,cho
      REAL(KIND=DBL), INTENT(IN) :: clono2e,clxe,co,hbre,hcle,h2o,h2o2e
      REAL(KIND=DBL), INTENT(IN) :: hno3e,ho2no2e,hobre,hocle,noxe,o1d
      REAL(KIND=DBL), INTENT(IN) :: o3,o3p,rbrobrx,rclclo,rcloclx,rnno
      REAL(KIND=DBL), INTENT(IN) :: rno2nox,rnono2,rbrnbrx

      REAL(KIND=DBL), INTENT(INOUT) :: h,ho2,oh

c... Iteration is used to approximate the [H], [OH],
c... and [HO2].
 
      INTEGER :: iter

      REAL :: o2m
      REAL(KIND=DBL) :: a,b,c,gbrx,gch4,gclx,gh2o,ghc,gnox,goxhox,rad
      REAL(KIND=DBL) :: rho2oh,rhoh,tbre,tbroe,tbrono2e,tcle,tcloe,term
      REAL(KIND=DBL) :: tne,tnoe,tno2e

c... Niters iterations are performed to estimate the [H], [OH],
c... and [HO2].  At dawn the concentrations of H, OH, HO2, O3P,
c... O1D, and NO are equal to their evening concentrations. 
c... O2 mixing ratio is specified

      if(k .ge. k1Strat .and. k .le. levels_cal) then

          if(daytime) then  
            o2m = o2*m(k)
            tcloe = rcloclx*clxe
            tcle = rclclo*tcloe

            tbroe = rbrobrx*brxe
            tbrono2e = rbrnbrx*brxe
            tbre = brxe-tbroe-tbrono2e

            tno2e = rno2nox*noxe        
            tnoe = rnono2*tno2e
            tne = rnno*tnoe
            do iter=1,niters
c... Define rho2oh=[HO2]/[OH] and rhoh=[H]/[OH].             
              rho2oh = ho2/oh 
              rhoh = h/oh
      

c... Odd hydrogen(H+OH+HO2) equation is written in the form
c...  g7*oh**2 + g8*oh + g9 and solved for oh. 
              a = -4.*(ak(k,59)
     &         +ak(k,12)*m(k)
     &         +(ak(k,40) +ak(k,65)*m(k)
     &         +(ak(k,71) +ak(k,72))*rhoh
     &         +(ak(k,19) +ak(k,64)*m(k))*rho2oh)
     &         *rho2oh)

              goxhox = 0.

              gnox = ak(k,37)*hno3e
     &         +ak(k,55)*ho2no2e
     &         +ak(k,14)*ch4
     &         +( ak(k,18)+ak(k,34)*rho2oh )*tno2e*m(k)

Crk rho2oh appended to ak(53) and ak(6) terms below

              gclx = ( ak(k,6)*tcloe+ak(k,53)*tcle )*rho2oh
     &         +ak(k,61)*tcloe
     &         +ak(k,13)*clono2e
     &         +ak(k,27)*hcle
     &         +ak(k,62)*hocle  

              gbrx = ak(k,95)*hbre
     &         +( ak(k,109)*tbroe+ak(k,92)*tbre )*rho2oh
          
              ghc= ak(k,22)*ch3o2*rho2oh
     &         +0.7*ak(k,58)*ch3ooh 
     &         +ak(k,51)*ch2o

              b = -goxhox -gnox -gclx -gbrx -ghc
        
              goxhox = 2.*(
     &          (aj(4)+ak(k,39)*o1d)*h2o
     &         +(aj(8)+ak(k,33)*o3p)*h2o2e)

              gnox = aj(6)*hno3e
     &         +( aj(24)+ak(k,70)*m(k) )*ho2no2e
     &         +ak(k,125)*h*tno2e

Crk *o1d term appended to following
              gch4 = ak(k,49)*ch4*o1d 

              gclx = aj(18)*hcle
     &         +aj(23)*hocle 
     &         +ak(k,7)*tcle*h2o2e

              gbrx = aj(25)*hobre                           !Crk
  
              ghc = aj(13)*ch3ooh 
     &         +( ak(k,52)*cho
     &         +ak(k,17)*ch3o )*o2m
     &         +( aj(10)+ak(k,21)*o3p )
     &         *ch2o

              gh2o = aj(38)*h2o

              c = goxhox +gnox +gch4 +gclx +gbrx +ghc +gh2o
         
c... Average oh concentrations at time t and t+tdt to
c... speed convergence
c... took factor of 2*a to above eqn for a - 8/98(sds)
              rad = abs(b**2-2.*a*c)
              rad = sqrt(rad)  
              term = max(1.,-b/a +abs(rad/a)) 
c... use Harvard climo for OH in the troposphere
              if(tropo(k) .and. k .le. levels_cal) then
                oh = harvard_oh(k)
               else
                oh = 0.5*( oh+term )
               endif

c... Solve for H concentration assuming P.E.
              h  = (aj(18)*hcle 
     &         +( ak(k,68)*tne
     &         +ak(k,41)*o3p
     &         +ak(k,36)*co )*oh 
     &         +ak(k,119)*ch4*o1d
     &         +aj(42)*ch4
     &         +aj(4)*h2o 
     &         +aj(10)*ch2o
     &         +aj(38)*h2o )
     &         /(ak(k,3)*o2m*m(k)
     &         +ak(k,125)*tno2e
     &         +ak(k,11)*o3
     &         +( ak(k,71)+ak(k,72)+ak(k,73) )
     &          *oh*rho2oh )
 
c... HO2 equation is written in the form g4*ho2**2+g5*ho2+g6
c... and solved for HO2. 

              goxhox =
     &         ak(k,3)*h*o2m*m(k)
     &         +ak(k,4)*oh*o3
     &         +( ak(k,29)*oh
     &         +ak(k,33)*o3p )*h2o2e  

Crk         g6nox = ak(k,70)*ho2no2e*m(k) 
              gnox = ( ak(k,70)*m(k)+aj(24) )
     &         *ho2no2e

              gclx = ak(k,7)*tcle*h2o2e
     &         +ak(k,60)*tcloe*oh

              ghc = ( ak(k,52)*cho
     &         +ak(k,17)*ch3o )*o2m

              c = goxhox +gnox +gclx +ghc
   
              goxhox = ak(k,5)*o3
     &         +ak(k,42)*o3p
     &         +( ak(k,40)+ak(k,65)*m(k) )*oh
     &         +( ak(k,71)+ak(k,72)+ak(k,73) )*h

              gnox = ak(k,34)*tno2e*m(k)
     &         +ak(k,38)*tnoe

              gclx = ( ak(k,53)+ak(k,67) )*tcle
     &         +ak(k,6)*tcloe

              gbrx = ak(k,109)*tbroe +ak(k,92)*tbre

              ghc = ak(k,22)*ch3o2
c... added ghc in (was calcd but not used - sds 8/98)
              b = -goxhox -gnox -gclx -gbrx -ghc

              a = -4.*(ak(k,19)+ak(k,64)*m(k))
         
c... took factor of 2*a to above eqn for a - 8/98(sds)
              rad = abs(b**2-2.*a*c)
              rad = sqrt(rad)  
              ho2 = max(1.,-b/a+abs(rad/a)) 

             enddo
           endif

       endif

      RETURN
      END SUBROUTINE hoxpartest
