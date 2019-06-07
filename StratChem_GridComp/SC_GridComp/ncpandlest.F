      SUBROUTINE ncpandlest(dt,k,k1Strat,levels,daytime,m,ak,aj,brxe,ccl4,ch3ccl3,ch3cl,ch4,clono2e,
     &                         clxe,f11,f113,f12,f114,f115,h2o,h2o2e,hcfc22,hcfc141b,hcfc142b,h1211,
     &                      ch2brcl,chbrcl2,chbr2cl,hcl,hcle,hno3e,ho2,ho2no2,ho2no2e,hocle,lclono2,
     &                   lclono2het,lclx,lhcl,lhno3,lho2no2,lhocl,ln2o5,ln2o5het,lnox,lnoxcl,lnoxsq,
     &                 n2o,n2o5,n2o5e,ne,noe,no2e,noxe,o3e,o3pe,oh,pclono2,pclx,phcl,phno3,phno3het,
     &                          pho2no2,phocl,pn2o5,pnox,pnoxcl,rain,rbrnbrx,rclclo,rcloclx,rno2nox,
     &                                           rno3no2,ro1do3,levels_cal,numreacs,numphoto,o2,gcr)

c... Subroutine to calculate final production and loss terms for 
c... NOy, Cly species;  the final value is the AVERAGE of a production and
c... loss term using the estimated concentrations and the production and loss
c... terms calculated using the previous time step concentrations
c... Recall that Br P and L have already been finalized.

      IMPLICIT NONE
      INTEGER, PARAMETER :: DBL = KIND(0.00D+00)

      INTEGER, INTENT(IN) :: k,k1Strat,levels,levels_cal,numreacs,numphoto
      LOGICAL, INTENT(IN) :: daytime
      REAL, INTENT(IN) :: dt,m(levels),rain(levels),o2

      REAL(KIND=DBL), INTENT(IN) :: aj(numphoto),ak(levels,numreacs)
      REAL(KIND=DBL), INTENT(IN) :: brxe,ccl4,ch3ccl3,ch3cl,ch4,clono2e,clxe
      REAL(KIND=DBL), INTENT(IN) :: f11,f113,f12,f114,f115
      REAL(KIND=DBL), INTENT(IN) :: h2o,h2o2e,hcl,hcle,hno3e
      REAL(KIND=DBL), INTENT(IN) :: hcfc22,hcfc141b,hcfc142b,h1211
      REAL(KIND=DBL), INTENT(IN) :: ch2brcl,chbrcl2,chbr2cl
      REAL(KIND=DBL), INTENT(IN) :: ho2,ho2no2,ho2no2e,hocle,n2o,n2o5,n2o5e,ne
      REAL(KIND=DBL), INTENT(IN) :: noe,no2e,noxe,o3e,o3pe,oh,rclclo,rcloclx
      REAL(KIND=DBL), INTENT(IN) :: rno2nox,rno3no2,rbrnbrx,ro1do3

      REAL, INTENT(IN)  :: gcr(levels)

      REAL(KIND=DBL), INTENT(INOUT) :: pclono2,pnoxcl,pclx,phcl,phno3,phno3het
      REAL(KIND=DBL), INTENT(INOUT) :: pho2no2,phocl,pn2o5,pnox,lclono2
      REAL(KIND=DBL), INTENT(INOUT) :: lclono2het,lclx,lhcl,lhno3,lho2no2,lhocl
      REAL(KIND=DBL), INTENT(INOUT) :: ln2o5,ln2o5het,lnox,lnoxcl,lnoxsq

      REAL(KIND=DBL) :: o1de,tcle,tcloe,tno2e,tno3e
 
      if(k .ge. k1Strat .and. k .le. levels_cal) then
          if(daytime) then 
c... set up estimated o3p
            tno2e = rno2nox*noxe
            tno3e = rno3no2*tno2e
            tcloe = rcloclx*clxe
            tcle = rclclo*tcloe
            o1de = ro1do3*o3e
c... det prod and loss of hocl, hcl, clx, and clono2
c... AVERAGE for this solution scheme

            phocl = .5*( phocl
     &        +ak(k,6)*ho2*tcloe
     &        +( ak(k,13)*oh
     &        +ak(k,112) )*clono2e )

            lhocl = .5*( lhocl
     &        +( ak(k,62)*oh
     &        +aj(23)
     &        +ak(k,114)*hcle )*hocle )
      
             phcl= .5*( phcl
     &        +( ak(k,7)*h2o2e
     &        +ak(k,26)*ch4
     &        +ak(k,53)*ho2 )*tcle
     &        +ak(k,61)*oh*tcloe )
         
            lhcl = .5*( lhcl
     &       +( ak(k,27)*oh
     &       +aj(18)
     &       +ak(k,111)*clono2e
     &       +ak(k,113)*n2o5e
     &       +ak(k,114)*hocle
     &       +rain(k) )*hcle )

            pclono2 =.5*( pclono2
     &       +ak(k,32)*tcloe*tno2e*m(k) )

            lclono2 = .5*( lclono2
     &       +( ak(k,35)*o3pe
     &       +ak(k,13)*oh
     &       +aj(15)
     &       +rain(k) )*clono2e )

            lclono2het = .5*( lclono2het
     &       +(ak(k,111)*hcle
     &       +ak(k,112))*clono2e )

c... clx prod from hcl ,clono2, and hocl       
c... substitute following if rain
            pclx = .5*( pclx
     &       +( ak(k,27)*oh +aj(18)
     &       +2.*ak(k,111)*clono2e
     &       +ak(k,113)*n2o5e
     &       +2.*ak(k,114)*hocle )*hcle
     &       +(ak(k,35)*o3pe+aj(15))*clono2e
     &       +(ak(k,62)*oh+aj(23))*hocle
c... was clxsource from 2d model
     &       +ak(k,16)*ch3cl*oh 
     &       +3*ak(k,75)*ch3ccl3*oh
     &       +  ak(k,98)*hcfc22*oh

#ifndef REDUCED
     &       +2*ak(k,129)*hcfc141b*oh
     &       +  ak(k,130)*hcfc142b*oh
#endif

     &       +2*ak(k,80)*f12*o1de
     &       +3*ak(k,83)*f11*o1de
     &       +3*ak(k,99)*f113*o1de

#ifndef REDUCED
     &       +3*ak(k,148)*f114*o1de
     &       +3*ak(k,149)*f115*o1de
#endif

     &       +4*ak(k,121)*ccl4*o1de
     &       +  ak(k,126)*hcfc22*o1de

#ifndef REDUCED
     &       +2*ak(k,127)*hcfc141b*o1de
     &       +  ak(k,128)*hcfc142b*o1de
     &       +  ak(k,117)*h1211*o1de
#endif

     &       +3*aj(29)*f11
     &       +2*aj(30)*f12
     &       +4*aj(31)*ccl4
     &       +3*aj(32)*ch3ccl3
     &       +  aj(33)*hcfc22
     &       +3*aj(34)*f113
     &       +  aj(35)*ch3cl 
     &       +  aj(22)*h1211

#ifndef REDUCED
     &       +  (aj(46)+ak(k,133)*oh)*ch2brcl
     &       +2*(aj(47)+ak(k,134)*oh)*chbrcl2 
     &       +  (aj(48)+ak(k,135)*oh)*chbr2cl
#endif

     &       ) 
c... end clxsource

            lclx= .5*( lclx
     &       +( ak(k,6)*ho2
     &       +ak(k,32)*m(k)*tno2e
     &       +ak(k,61)*oh )*tcloe
     &       +( ak(k,7)*h2o2e
     &       +ak(k,26)*ch4
     &       +ak(k,53)*ho2 )*tcle )

c... NOx species parts
c... Define HNO3, N2O5, and HO2NO2 production and loss.  
c... AVERAGE production and loss 

            phno3 = .5*( phno3
     &       +ak(k,18)*oh*tno2e*m(k) 
     &       +2.*ak(k,79)*n2o5e*h2o)  

            phno3het=.5*( phno3het
     &       +2.*ak(k,110)*n2o5e
     &       +ak(k,113)*hcle*n2o5e
     &       +ak(k,111)*hcle*clono2e
     &       +ak(k,112)*clono2e
     &       +ak(k,115)*rbrnbrx*brxe)
	 
            lhno3 = .5*( lhno3
     &       +( ak(k,37)*oh
     &       +ak(k,54)*o3pe
     &       +aj(6) +rain(k) )*hno3e )
     
            pn2o5 = .5*( pn2o5
     &       +ak(k,46)*m(k)*tno2e*tno3e)

            ln2o5 = .5*( ln2o5
     &       +(ak(k,31)*m(k)
     &       +aj(9)
     &       +ak(k,79)*h2o )*n2o5e )

c... note:  Reaction 110 assumes rxn is n2o5 limited so ak(110)
c... is like a first order loss. (do not need to mult. by h2o) 	 

            ln2o5het = .5*( ln2o5het
     &       +(ak(k,110)+ak(k,113)*hcle )*n2o5e)

            pho2no2 = .5*( pho2no2
     &       +ak(k,34)*ho2*tno2e*m(k) )

            lho2no2 = .5*( lho2no2
     &       +( ak(k,55)*oh
     &       +ak(k,70)*m(k)
     &       +aj(24) +rain(k) )*ho2no2e )

c... constrain loss of ho2no2
            if(lho2no2*dt.ge.ho2no2) 
     &        lho2no2 = ho2no2/dt

            pnox = 2.*ln2o5+lhno3+lho2no2
     &       +0.5*ak(k,113)*( n2o5e*hcle
     &       +n2o5*hcl ) +gcr(k)
     &       +2.*ak(k,45)*n2o*o1de

c... production of nox from chlorine species
c... take care of brono2 term in solver

            pnoxcl = lclono2

c.. loss terms

            lnox = phno3+pho2no2
     &        +2.0*ak(k,47)*ne*noe+2.0*ak(k,81)*ne*tno2e
  
c... loss of nox to form chlorine species, do brono2 in solver
            lnoxcl = pclono2
         
            lnoxsq = 2.*pn2o5

           endif 

       endif

      RETURN
      END SUBROUTINE ncpandlest
