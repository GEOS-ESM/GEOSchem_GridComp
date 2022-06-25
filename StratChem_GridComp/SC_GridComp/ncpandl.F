      SUBROUTINE ncpandl(dt,k,k1Strat,levels,daytime,m,ak,aj,ccl4,ch3ccl3,ch3cl,ch4,clono2,
     &          clx,f11,f113,f12,f114,f115,h2o,h2o2,hcfc22,hcfc141b,hcfc142b,h1211,ch2brcl,
     &              chbrcl2,chbr2cl,hcl,hno3,ho2,ho2no2,hocl,lbrono2het,lclono2,lclono2het,
     &                     lclx,lhcl,lhno3,lho2no2,lhocl,ln2o5,ln2o5het,lnox,lnoxcl,lnoxsq,
     &           n,n2o,n2o5,no,no2,nox,o1d,o3p,oh,pbrono2,pclono2,pclx,phcl,phno3,phno3het,
     &                         pho2no2,phocl,pn2o5,pnox,pnoxcl,rain,rclclo,rcloclx,rno2nox,
     &                                         rno3no2,levels_cal,numreacs,numphoto,o2,gcr)

c  First subroutine to calculate the production and loss terms for
c  noy, cly species
Crk  951107  Updated to include new Br scheme.  Br p and l now done separately.

      IMPLICIT NONE
      INTEGER, PARAMETER :: DBL = KIND(0.00D+00)

      INTEGER, INTENT(IN) :: k,k1Strat,levels,levels_cal,numreacs,numphoto
      LOGICAL, INTENT(IN) :: daytime
      REAL, INTENT(IN) :: dt,m(levels),rain(levels),o2

      REAL(KIND=DBL), INTENT(IN) :: aj(numphoto),ak(levels,numreacs)
      REAL(KIND=DBL), INTENT(IN) :: ccl4,ch3ccl3,ch3cl,ch4,clono2,clx
      REAL(KIND=DBL), INTENT(IN) :: f11,f12,f113,f114,f115
      REAL(KIND=DBL), INTENT(IN) :: h2o,h2o2,hcl,hno3,ho2,ho2no2
      REAL(KIND=DBL), INTENT(IN) :: hcfc22,hcfc141b,hcfc142b,h1211
      REAL(KIND=DBL), INTENT(IN) :: ch2brcl,chbrcl2,chbr2cl
      REAL(KIND=DBL), INTENT(IN) :: hocl,lbrono2het,n,n2o,n2o5,no,no2,nox
      REAL(KIND=DBL), INTENT(IN) :: o1d,o3p,oh,pbrono2,rclclo,rcloclx
      REAL(KIND=DBL), INTENT(IN) :: rno2nox,rno3no2

      REAL, INTENT(IN) :: gcr(levels)

      REAL(KIND=DBL), INTENT(INOUT) :: lclono2,lclono2het,lclx,lhcl,lhno3,lho2no2
      REAL(KIND=DBL), INTENT(INOUT) :: lhocl,ln2o5,ln2o5het,lnox,lnoxcl,lnoxsq
      REAL(KIND=DBL), INTENT(INOUT) :: pclono2,pclx,phcl,phno3,phno3het
      REAL(KIND=DBL), INTENT(INOUT) :: pho2no2,phocl,pn2o5,pnox,pnoxcl

      REAL(KIND=DBL) :: tcl,tclo,tno2,tno3
      
      if(k .ge. k1Strat .and. k .le. levels_cal) then

          if (daytime) then 
            tno2 = rno2nox*nox
            tno3 = rno3no2*tno2
            tclo = rcloclx*clx
            tcl = rclclo*tclo

c det prod and loss of hocl, hcl, clx, and clono2

            phocl = ak(k,6)*ho2*tclo
     &       +( ak(k,13)*oh
     &       +ak(k,112) )*clono2

            lhocl = ( ak(k,62)*oh
     &       +aj(23)
     &       +ak(k,114)*hcl )*hocl
      
            phcl = ( ak(k,7)*h2o2
     &       +ak(k,26)*ch4
     &       +ak(k,53)*ho2 )*tcl
     &       +ak(k,61)*oh*tclo
     
            lhcl = ( ak(k,27)*oh
     &       +aj(18)
     &       +ak(k,111)*clono2
     &       +ak(k,113)*n2o5
     &       +ak(k,114)*hocl
     &       +rain(k) )*hcl
c
            pclono2 = ak(k,32)*m(k)*tclo*tno2

            lclono2 = ( ak(k,35)*o3p
     &        +ak(k,13)*oh
     &        +aj(15) +rain(k) )*clono2

            lclono2het = ( ak(k,111)*hcl
     &         +ak(k,112) )*clono2

c... clx prod from hcl ,clono2, and hocl 
c... NOTE!!! Cl2 produced from het rxns 111 and 114 rapidly
c... photolyzes and leads to a production of 2 Clx's.  The
c... ones produced by 114 are included in lhcl and lhocl. 
c... 2 is not included in front of ak111 since one
c... clx molecule is already included in lhcl. 

            pclx = lhcl+lhocl
c... was clxsource from 2d model
     &       +ak(k,16)*ch3cl*oh
     &       +3*ak(k,75)*ch3ccl3*oh
     &       +  ak(k,98)*hcfc22*oh
     &       +2*ak(k,129)*hcfc141b*oh
     &       +  ak(k,130)*hcfc142b*oh
     &       +2*ak(k,80)*f12*o1d
     &       +3*ak(k,83)*f11*o1d
     &       +3*ak(k,99)*f113*o1d

#ifndef REDUCED
     &       +2*ak(k,148)*f114*o1d
     &       +  ak(k,149)*f115*o1d
#endif

     &       +4*ak(k,121)*ccl4*o1d
     &       +  ak(k,126)*hcfc22*o1d

#ifndef REDUCED
     &       +2*ak(k,127)*hcfc141b*o1d
     &       +  ak(k,128)*hcfc142b*o1d
#endif

     &       +  ak(k,117)*h1211*o1d
     &       +3*aj(29)*f11
     &       +2*aj(30)*f12
     &       +4*aj(31)*ccl4
     &       +3*aj(32)*ch3ccl3
     &       +  aj(33)*hcfc22

#ifndef REDUCED
     &       +2*aj(51)*hcfc141b
     &       +  aj(52)*hcfc142b
#endif

     &       +3*aj(34)*f113

#ifndef REDUCED
     &       +2*aj(53)*f114
     &       +  aj(54)*f115
#endif

     &       +  aj(35)*ch3cl
     &       +  aj(22)*h1211

#ifndef REDUCED
     &       +  (aj(46)+ak(k,133)*oh)*ch2brcl
     &       +2*(aj(47)+ak(k,134)*oh)*chbrcl2
     &       +  (aj(48)+ak(k,135)*oh)*chbr2cl
#endif

c... end clxsource
     &       +( ak(k,35)*o3p +aj(15)
     &       +ak(k,111)*hcl )*clono2
     &       -rain(k)*hcl 

            lclx = ak(k,6)*ho2*tclo
     &       +pclono2+phcl
    
c... Define HNO3, N2O5, and HO2NO2 production and loss.  

            phno3 = 2.*ak(k,79)*n2o5*h2o
     &       +ak(k,18)*oh*tno2*m(k)

            phno3het = ( 2.*ak(k,110)
     &       +ak(k,113)*hcl )*n2o5
     &       +lclono2het+lbrono2het
            
            lhno3 = ( ak(k,37)*oh
     &       +ak(k,54)*o3p
     &       +aj(6) +rain(k) )*hno3

            pn2o5 = ak(k,46)*m(k)*tno2*tno3

            ln2o5 = ( ak(k,31)*m(k)
     &       +aj(9)
     &       +ak(k,79)*h2o )*n2o5

c... note:  Reaction 110 assumes rxn is n2o5 limited so ak(110)
c... is like a first order loss. (do not need to mult. by h2o) 

            ln2o5het = ( ak(k,110)
     &       +ak(k,113)*hcl )*n2o5

            pho2no2 = ak(k,34)*ho2
     &       *tno2*m(k) 

            lho2no2 = ( ak(k,55)*oh
     &       +ak(k,70)*m(k)
     &       +aj(24) +rain(k) )*ho2no2
c... constrain loss of ho2no2
            if(lho2no2*dt.ge.ho2no2) 
     &        lho2no2 = ho2no2/dt

c... Production of nox is defined. 
            pnox = 2.*ln2o5+lhno3+lho2no2
     &       -rain(k)*( hno3+ho2no2 )
     &       +ak(k,113)*n2o5*hcl 

            pnox = pnox +gcr(k)
     &       +2.*ak(k,45)*n2o*o1d 
c
c... production of NOx from clono2 
C... handle brono2 in solver
c
            pnoxcl = lclono2 -rain(k)*clono2
						
            lnox = phno3 +pho2no2
     &        +2.0*ak(k,47)*n*no+2.0*ak(k,81)*n*tno2

c... loss of NOx to chlorine species, not brono2
            lnoxcl = pclono2
c...            lnoxcl =pclono2+pbrono2
         
            lnoxsq=2.*pn2o5

           endif 

       endif

      RETURN
      END SUBROUTINE ncpandl
