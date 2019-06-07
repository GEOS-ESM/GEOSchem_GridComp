      SUBROUTINE hoxpart(k,k1Strat,levels,daytime,m,tropo,harvard_oh,ak,aj,brx,ch2o,ch3o,
     &             ch3o2,ch3ooh,ch4,cho,clono2,clx,co,h,hbr,hcl,h2o,h2o2,hno3,ho2,ho2no2,
     &           hobr,hocl,nox,o1d,o3,o3p,oh,rbrnbrx,rbrobrx,rclclo,rcloclx,rnno,rno2nox,
     &                                     rnono2,levels_cal,niters,numreacs,numphoto,o2)

      IMPLICIT NONE
      INTEGER, PARAMETER :: DBL = KIND(0.00D+00)

      INTEGER, INTENT(IN) :: k,k1Strat,levels
      INTEGER, INTENT(IN) :: levels_cal,niters,numreacs,numphoto
      LOGICAL, INTENT(IN) :: daytime,tropo(levels)
      REAL, INTENT(IN) :: m(levels),harvard_oh(levels),o2

      REAL(KIND=DBL), INTENT(IN) :: aj(numphoto),ak(levels,numreacs)
      REAL(KIND=DBL), INTENT(IN) :: brx,ch2o,ch3o,ch3o2,ch3ooh,ch4,cho
      REAL(KIND=DBL), INTENT(IN) :: clono2,clx,co,hbr,hcl,h2o,h2o2,hno3
      REAL(KIND=DBL), INTENT(IN) :: ho2no2,hobr,hocl,nox,o1d,o3,o3p
      REAL(KIND=DBL), INTENT(IN) :: rbrobrx,rclclo,rcloclx,rnno,rno2nox
      REAL(KIND=DBL), INTENT(IN) :: rnono2,rbrnbrx

      REAL(KIND=DBL), INTENT(INOUT) :: h,ho2,oh

c     Latest revision 26 October 1994 
Crk   Revision 9 November 1994 to add HOBr
Crk   Revision 951107 to add HC chem and update to new Br scheme

c... Iterated to approximate H,OH,HO2.  

      INTEGER :: iter

      REAL :: o2m
      REAL(KIND=DBL) :: a,b,c,gbrx,gch4,gclx
      REAL(KIND=DBL) :: gh2o,ghc,gnox,goxhox,rad
      REAL(KIND=DBL) :: rho2oh,rhoh,tbr,tbro,tbrono2
      REAL(KIND=DBL) :: tcl,tclo,term,tn,tno,tno2

c... Niters iterations are performed to estimate the [H], [OH],
c... and [HO2].  At dawn the concentrations of H, OH, HO2, O3P,
c... O1D, and NO are equal to their evening concentrations.
c... O2 mixing ratio is specified

      if(k .ge. k1Strat .and. k .le. levels_cal) then

          if (daytime) then   

            o2m = o2*m(k)
            tclo = rcloclx*clx
            tcl = rclclo*tclo

            tbro = rbrobrx*brx
            tbrono2 = rbrnbrx*brx
            tbr = brx-tbro-tbrono2

            tno2 = rno2nox*nox    
            tno = rnono2*tno2
            tn = rnno*tno
c... Define rho2oh=[HO2]/[OH] and rhoh=[H]/[OH].      
            do iter=1,niters
              rho2oh = ho2 / oh 
              rhoh =   h / oh

c... Odd hydrogen(H+OH+HO2) equation is written in the form
c... a*oh**2 + g8*oh + c and solved for oh

              a = -4.*( ak(k,59)
     &         +ak(k,12)*m(k)
     &         +( ak(k,40)
     &         +ak(k,65)*m(k)
     &         +(ak(k,71)+ak(k,72))*rhoh
     &         +(ak(k,19)+ak(k,64)*m(k) )*rho2oh )
     &         *rho2oh )

c... oxhox term
              goxhox = 0.

c... nox term
              gnox = ak(k,37)*hno3
     &         +( ak(k,18) +ak(k,34)*rho2oh )*tno2*m(k)
     &         +ak(k,55)*ho2no2
     &         +ak(k,14)*ch4

c... clx term
C rk  rho2oh added to ak(53) and ak(6) terms below
              gclx = ( ak(k,6)*tclo+ak(k,53)*tcl )*rho2oh
     &         +ak(k,13)*clono2
     &         +ak(k,27)*hcl
     &         +ak(k,62)*hocl
     &         +ak(k,61)*tclo

c... brx term
              gbrx = ( ak(k,109)*tbro +ak(k,92)*tbr )*rho2oh
     &         +ak(k,95)*hbr

c... hc term
              ghc = ak(k,22)*ch3o2*rho2oh
     &         +0.7*ak(k,58)*ch3ooh 
     &         +ak(k,51)*ch2o
c... reversed sign of g8, except in rad calc below (sds 8/98)
              b = -goxhox -gnox -gclx -gbrx -ghc
       
              goxhox = 2.*(
     &         ( aj(4)+ak(k,39)*o1d )*h2o
     &         +( aj(8)+ak(k,33)*o3p )*h2o2)

              gnox = aj(6)*hno3
     &         +( aj(24)+ak(k,70)*m(k) )*ho2no2
     &         +ak(k,125)*h*tno2

C rk *o1d term appended to following

              gch4 = ak(k,49)*ch4*o1d

              gclx = ak(k,7)*tcl*h2o2
     &         +aj(18)*hcl
     &         +aj(23)*hocl 

              gbrx = aj(25)*hobr      !Crk

              ghc = aj(13)*ch3ooh
     &        +(ak(k,52)*cho+ak(k,17)*ch3o)
     &         *o2m
     &        +(aj(10)+ak(k,21)*o3p)*ch2o

              gh2o = aj(38)*h2o

              c = goxhox +gnox +gch4 +gclx +gbrx +ghc +gh2o

c        oh  = (oh +quad(a,g8,c)) * 0.5  
c... Average oh concentrations at time t and t+tdt to speed convergence
c
              rad = abs(b**2-2.*a*c)
              rad = sqrt(rad)  
              term = max(1.,-b/a+abs(rad/a)) 
c              oh = 0.5*( oh+term )
c... use Harvard climo for OH in the troposphere
              if(tropo(k) .and. k .le. levels_cal) then
                oh = harvard_oh(k)
               else
                oh = 0.5*( oh+term )
               endif

c... Solve for H concentration assuming P.E.
              h  = (aj(18)*hcl
     &         +( ak(k,68)*tn
     &         +ak(k,41)*o3p
     &         +ak(k,36)*co )*oh
     &         +aj(42)*ch4
     &         +ak(k,119)*ch4*o1d
     &         +aj(4)*h2o 
     &         +aj(10)*ch2o
     &         +aj(38)*h2o )
     &         /(ak(k,3)*o2m*m(k)
     &         +ak(k,11)*o3
     &         +ak(k,125)*tno2
     &         +( ak(k,71)+ak(k,72)+ak(k,73) )
     &         *oh*rho2oh)
     
c... HO2 equation is written in the form g4*ho2**2+g5*ho2+g6
c... and solved for HO2
              goxhox =
     &         ak(k,3)*h*o2m*m(k)
     &         +ak(k,4)*oh*o3
     &         +ak(k,29)*oh*h2o2
     &         +ak(k,33)*o3p*h2o2  

c...         gnox = ak(k,70)*ho2no2*m(k) 
              gnox = ( ak(k,70)*m(k)+aj(24) )
     &         *ho2no2

              gclx = ak(k,7)*tcl*h2o2
     &         +ak(k,60)*tclo*oh

              ghc = ( ak(k,52)*cho
     &         +ak(k,17)*ch3o )*o2m

              c = goxhox +gnox +gclx +ghc
    
              goxhox = ak(k,5)*o3
     &         +ak(k,42)*o3p
     &         +( ak(k,40)+ak(k,65)*m(k) )*oh
     &         +( ak(k,71)+ak(k,72)+ak(k,73) )*h

              gnox = ak(k,34)*m(k)*tno2 +ak(k,38)*tno

              gclx = ak(k,6)*tclo+( ak(k,53)+ak(k,67))*tcl

              gbrx = ak(k,109)*tbro +ak(k,92)*tbr

              ghc = ak(k,22)*ch3o2

c... reversed sign of b, except in rad calc below (sds 8/98)
              b = -goxhox -gnox -gclx -gbrx -ghc
 
              a = -4.*(ak(k,19)+ak(k,64)*m(k))
  
              rad = abs(b**2-2.*a*c)
              rad = sqrt(rad)
              ho2 = max(1.,-b/a+abs(rad/a))
             enddo
           endif  

       endif

      RETURN
      END SUBROUTINE hoxpart
