      SUBROUTINE solverd2(dt,requiredDt,examineDt,k,k1Strat,levels,daytime,tropo,m,aj,ak,br,
     &        brcl,bro,brono2,brx,ccl4,ch2o,ch3br,ch3ccl3,ch3cl,ch3o2,ch3ooh,ch4,cho,cl,cl2, 
     &              chbr3,ch2br2,ch2brcl,chbrcl2,chbr2cl,co2b,sf6b,cl2o2,clo,clono2,clono2e, 
     &            clx,clxe,co,f11,f113,f12,f114,f115,h,h1211,h1301,h1202,h2402,h2o,h2o2,hbr, 
     &         hcfc22,hcfc141b,hcfc142b,hfc23,hfc32,hfc125,hfc134a,hfc143a,hfc152a,hcl,hcle, 
     &   hno3,ho2,ho2no2,hobr,hocl,hocle,lbrono2het,lbrx,lclono2,lclono2het,lclx,lh2o2,lhbr, 
     &            lhcl,lhno3,lho2no2,lhobr,lhocl,ln2o5,ln2o5het,lno2,lno3,lnox,lnoxa,lnoxcl, 
     &       lnoxsq,lo3,n,n2o,n2o5,n2o5e,no2,no3,no,nox,o1d,o3,o3p,oclo,oh,ox,pbrcl,pbrono2, 
     &      pbrx,pcl2,pclono2,pclx,ph2o2,phbr,phcl,pho2no2,phocl,phno3,phno3het,phobr,pn2o5, 
     &                 pno2,pno3,pnox,pnoxa,pnoxcl,po3,poclo,rain,rbrnbrx,rbrobrx,rcl2o2clo, 
     &    rclclo,rcloclx,rnno,rnono2,rno2nox,rno3no2,ro1do3,ro3ox,ro3po3,rocloclo,box_ro3ox, 
     &                                             levels_cal,numreacs,numphoto,gcr,o2,szan)

c... Currently no production or loss of ox at night (rxn 112)

      IMPLICIT NONE
      INTEGER, PARAMETER :: DBL = KIND(0.00D+00)
      REAL, PARAMETER :: PI = 3.14159265

      INTEGER, INTENT(IN) :: k,k1Strat,examineDt,levels
      INTEGER, INTENT(IN) :: levels_cal,numreacs,numphoto
      LOGICAL, INTENT(IN) :: daytime,tropo(levels)
      REAL, INTENT(IN) :: dt,o2,m(levels),rain(levels),gcr(levels),szan

      REAL(KIND=DBL), INTENT(IN) :: aj(numphoto),ak(levels,numreacs)
      REAL(KIND=DBL), INTENT(IN) :: ch2o,ch3o2,cho,clono2e,clxe,h,hcle,ho2,hocle
      REAL(KIND=DBL), INTENT(IN) :: lbrx,lhbr,lhobr,lo3,n2o5e,oh,pbrx,ph2o2
      REAL(KIND=DBL), INTENT(IN) :: phbr,phobr,pnoxcl,po3,rbrobrx,rclclo
      REAL(KIND=DBL), INTENT(IN) :: rcloclx,rnno,rnono2,rno2nox,rno3no2,ro3ox
      REAL(KIND=DBL), INTENT(IN) :: ro3po3,rocloclo,ro1do3

      REAL, INTENT(INOUT) :: requiredDt

      REAL(KIND=DBL), INTENT(INOUT) :: br,brcl,bro,brono2,brx,ccl4,ch3br,ch3ccl3
      REAL(KIND=DBL), INTENT(INOUT) :: ch3cl,ch3ooh,ch4,cl,cl2,cl2o2,clo,clono2
      REAL(KIND=DBL), INTENT(INOUT) :: clx,co,f11,f113,f12,f114,f115,hcfc22
      REAL(KIND=DBL), INTENT(INOUT) :: hcfc141b,hcfc142b,h1301,h1211,h1202,h2402
      REAL(KIND=DBL), INTENT(INOUT) :: hfc23,hfc32,hfc125,hfc134a
      REAL(KIND=DBL), INTENT(INOUT) :: hfc143a,hfc152a,co2b,sf6b
      REAL(KIND=DBL), INTENT(INOUT) :: chbr3,ch2br2,ch2brcl,chbrcl2,chbr2cl 
      REAL(KIND=DBL), INTENT(INOUT) :: h2o,h2o2,hbr,hcl,hno3,ho2no2
      REAL(KIND=DBL), INTENT(INOUT) :: hobr,hocl,lbrono2het,lclono2,lclono2het
      REAL(KIND=DBL), INTENT(INOUT) :: lclx,lh2o2,lhcl,lhno3,lho2no2,lhocl,ln2o5
      REAL(KIND=DBL), INTENT(INOUT) :: ln2o5het,lno2,lno3,lnox,lnoxa,lnoxcl
      REAL(KIND=DBL), INTENT(INOUT) :: lnoxsq,n,n2o,n2o5,no2,no3,no,nox,o1d
      REAL(KIND=DBL), INTENT(INOUT) :: o3,o3p,oclo,ox,pbrcl,pbrono2,pcl2,pclono2
      REAL(KIND=DBL), INTENT(INOUT) :: pclx,phcl,pho2no2,phocl,phno3,phno3het
      REAL(KIND=DBL), INTENT(INOUT) :: pn2o5,pno2,pno3,pnox,pnoxa
      REAL(KIND=DBL), INTENT(INOUT) :: poclo,rbrnbrx,rcl2o2clo,box_ro3ox

      LOGICAL :: runVariableDt

      INTEGER :: itercount

      REAL(KIND=DBL) :: alimit,brntmp,cnoxbr,dbrn,dbrn0,h2oprod,h2oloss
      REAL(KIND=DBL) :: rat,rbrbro,sav,testclx,testhcl,testnox,tloss
      REAL(KIND=DBL) :: tmplclx

      REAL :: h2,hdt,minAllowDt,shortDt(3)

      hdt = 0.5*dt
      h2 = 5.00E-07

c... O2 mixing ratio is specified
c... NOX fix for denitrified conditions.  Scale nox loss

      if(k .ge. k1Strat .and. k .le. levels_cal) then

c...  daytime chemistry
               if(daytime) then 
                  
                  brx = brx+(pbrx-lbrx)*dt

c Prevent daytime negative Brx when forcing term is within truncation
c error of input Brx. MG/JEN 23 Feb 2004
                  if(brx .lt. 1.00e-20) brx=1.00e-20
                  
                  dbrn = (rbrnbrx*brx-brono2)/dt
c...  gas phase
                  dbrn0 = dbrn+lbrono2het
c...  [nox] lost to BrONO2
                  cnoxbr = max(dbrn0,0.)
                  
                  pnoxa = (pnox+pnoxcl)*dt 
                  lnoxa = (lnox+lnoxcl
     &                 +lnoxsq)*dt+cnoxbr*dt
                  testnox = nox+pnoxa-lnoxa
                  testclx = clx+dt*(pclx-lclx)
                  testhcl = hcl+dt*(phcl-lhcl)
                  itercount = 0

! If the tests indicate a negative number density will be produced,
! then find the time step length that guarantees a positive result,
! and return to start over ...

! ... unless the minimum allowable dt has already been requested.  In
! that case we hope the traditional method of 5 iterations works.

      minAllowDt = 60.00

      IF(dt <= minAllowDt) THEN
       runVariableDt=.FALSE.
      ELSE
       runVariableDt=.TRUE.
      END IF
      
      IF(runVariableDt) THEN

       if((testclx .le. 0.00 .or. testhcl .le. 0.00 .or. 
     &     testnox .le. 0.00) .and. examineDt .le. 3) then
        shortDt(1:3)=1.00e+10
        if(testnox .le. 0.00) shortDt(1)=nox/abs(pnox+pnoxcl-
     &                                   (lnox+lnoxcl+lnoxsq+cnoxbr))
        if(testclx .le. 0.00) shortDt(2)=clx/abs(pclx-lclx)
        if(testhcl .le. 0.00) shortDt(3)=hcl/abs(phcl-lhcl)
!       requiredDt=minval(shortDt)
        requiredDt=max(minval(shortDt),minAllowDt)
        return
       endif

      ELSE

       do while((testclx.lt.0..or.testhcl.lt.0..or.testnox.lt.0)
     &      .and.itercount.le.5)

                     if(testnox.lt.0) then
                        rat = nox/lnoxa  
                        
                        lnox = rat*lnox       
                        lnoxcl = rat*lnoxcl
                        lnoxsq = rat*lnoxsq 
                        lclx = lclx-(1.-rat)*pclono2
                        phno3 = rat*phno3
                        pho2no2 = rat*pho2no2

                        pclono2 = rat*pclono2
                        pn2o5 = rat*pn2o5
c...  NOx is lost to BrONO2
                        if (dbrn0.gt.0.) then
                           brntmp = rbrnbrx
                           dbrn0 = rat*dbrn0
                           dbrn = dbrn0-lbrono2het
                           rbrnbrx = (dbrn*dt+brono2)/brx
                           lbrono2het = rbrnbrx/brntmp
                        endif
                        testclx = clx+dt*(pclx-lclx)
                     endif      
                     
                     if (testclx < 0) then
                        sav = lclx 
                        lclx = clx/dt*
     &                       (1.-exp(-lclx*dt/clx))
                        rat = lclx / sav 
                        pclono2 = rat*pclono2
                        lnoxcl = pclono2
c...  taken out at some point          +pbrono2
                        phocl = rat*
     &                       ak(k,6)*ho2*rcloclx*clxe
     &                       +(ak(k,13)*oh
     &                       +ak(k,112))*clono2e
                        phcl = rat*phcl
c...  phcl changed, redo
                        testhcl = hcl+dt*(phcl-lhcl)
                     endif
                     
                     if(testhcl < 0) then
                        rat = (lhcl-hcl/dt)/
     &                       ((ak(k,27)*oh
     &                       +ak(k,111)*clono2e
     &                       +ak(k,113)*n2o5e
     &                       +ak(k,114)*hocle
     &                       +rain(k)
     &                       +aj(18))*hcle)
                        
                        sav = lhcl
                        lhcl = hcl/dt*
     &                       (1.-EXP(-lhcl*dt/hcl))
                        rat = lhcl/sav
                        lclono2het = lclono2het -rat * 
     &                       ak(k,111)*clono2e*hcle
                        phno3het = phno3het - rat * 
     &                       (ak(k,111)*clono2e*hcle+
     &                       ak(k,113)*n2o5e*hcle)   
                        pnox = pnox - rat * 
     &                       ak(k,113)*n2o5e*hcle
                        ln2o5het = ln2o5het - rat * 
     &                       ak(k,113)*n2o5e*hcle 
                        lhocl = lhocl - rat * 
     &                       ak(k,114)*hcle*hocle 
                        pclx = pclx
     &                       -rat*(ak(k,27)*hcle*oh 
     &                       +2.*ak(k,111)*clono2e*hcle 
     &                       +ak(k,113)*n2o5e*hcle 
     &                       +2.*ak(k,114)*hocle*hcle 
     &                       +aj(18)*hcle)
                     endif 
                     
                     pnoxa = (pnox+pnoxcl )*dt 
                     lnoxa = (lnox+lnoxcl
     &                    +lnoxsq)*dt+cnoxbr*dt
                     testnox = nox+pnoxa-lnoxa
                     testclx = clx+dt*(pclx-lclx)
                     testhcl = hcl+dt*(phcl-lhcl)

! MEM Added If statement to only print at levels well above tropopause
! LDO Added If statement to just print diagnostics if szan < 90 degress
! MEM Changed from itercount > 3  to  itercount > 4
    
       IF(k .ge. k1Strat+3) THEN
       IF(szan < PI*0.50) THEN
        IF(itercount > 4) THEN
         PRINT 345,k,itercount,testnox,testclx,testhcl
!        PRINT 346,nox,pnoxa,lnoxa,clx,pclx,lclx,hcl,phcl,lhcl
        END IF
       END IF
       END IF

  345 FORMAT(1h ,'SOLVERD2: (',i4,') count=',i1,' testnox=',1pe10.3,
     &       '  testclx=',1pe10.3,'  testhcl=',1pe10.3)
  346 FORMAT(1h ,' NOx=',1pe10.3,' pnoxa=',1pe10.3,' lnoxa=',1pe10.3,/,
     &       '  Clx=',1pe10.3,'  pclx=',1pe10.3,'  lclx=',1pe10.3,/,
     &       '  HCl=',1pe10.3,'  phcl=',1pe10.3,'  lhcl=',1pe10.3)

       itercount = itercount+1
       END DO

      END IF

c...  update N2O5,HO2NO2, H2O2, and HNO3
                  
                  n2o5 = n2o5 +(pn2o5-
     &                 (ln2o5+ln2o5het))*dt
                  
c     .sds... check for n2o5 < 0 and set small
                  if(n2o5 < 0.) then
                     n2o5 = 1.0e-15*m(k)
                  endif
                  
                  hno3 = hno3 + (phno3
     &                 +phno3het- lhno3)*dt
                  
                  h2o2 = h2o2 + (ph2o2-
     &                 lh2o2)*dt 
                  
                  ho2no2 = ho2no2+(pho2no2
     &                 -lho2no2)*dt 
                  
c...  update hydrocarbon, chlorine  and bromine species
                  n2o = n2o - (aj(14) + (ak(k,45)+
     &                 ak(k,78))*o1d) * n2o*dt
     &                 +ak(k,81)*n*no2*dt 

                  co2b = co2b-aj(12)*co2b
     &                  +(ak(k,36)*oh+ak(k,122)*o3p)*co *dt 

                  ch4 = ch4 -
     &                 (ak(k,14)*oh+aj(42)
     &                 +ak(k,26)*cl
     &                 +(ak(k,49)+ak(k,57)+ak(k,119))*o1d )
     &                 *ch4 *dt 
                  
                  co = co + (aj(11)*ch2o
     &                 +ak(k,52)*cho*m(k)*o2
     &                 +aj(12)*co2b
     &                 -(ak(k,36)*oh
     &                 +ak(k,122)*o3p)*co )*dt
                  
                  ch3ooh = ch3ooh + (ak(k,22)*
     &                 ch3o2*ho2 - (ak(k,58)*
     &                 oh + aj(13))*ch3ooh) * dt 
                  
                  clx = clx +(pclx-lclx)*dt
                  
                  clono2=clono2+(pclono2
     &                 -(lclono2+lclono2het))*dt
                  
c     .fxc... check for clono2 < 0 and reset small
                  if(clono2 < 0.) then
                     clono2 = 1.0e-15*m(k)
                  endif
                  
                  hcl = hcl+(phcl-lhcl)*dt
c     .sds... check for hcl < 0 and set small
                  if(hcl < 0.) then
                     hcl = 1.0e-15*m(k)
                  endif
                  
                  hocl = hocl+(phocl
     &                 -lhocl)*dt 
                  
                  hbr = hbr+(phbr-lhbr)*dt 
                  
                  hobr = hobr+(phobr-lhobr)*dt 
                  
c...  Update Ox and NOx
                  if(tropo(k)) then
                   ox = ox - lo3*dt
                  else
                   ox = ox +(po3-lo3)*dt
                  endif
c     .sds... check for ox very small and constrain
                  if(ox < 1.0e-10*m(k)) then
                     ox = 1.0e-10*m(k)
                  endif
                  
                  nox = nox+(pnox-dbrn0+pnoxcl
     &                 -(lnox+lnoxcl+lnoxsq))*dt
c     .sds... check for nox < 0 and set small
                  if(nox < 0.) then
                     nox = 1.0e-15*m(k)
                  endif
 
c      Update water (sds-10/10/2001).
c      Latest fixes: Mohan Gupta, 4 Jan 2005

                  h2oprod = ak(k,14)*ch4*oh*dt+
     &                 ak(k,16)*ch3cl*oh*dt+
     &                 ak(k,51)*ch2o*oh*dt+
     &                 ak(k,58)*ch3ooh*oh*dt+
     &                 ak(k,59)*oh*oh*dt+
     &                 ak(k,37)*oh*hno3*dt+
     &                 ak(k,29)*oh*h2o2*dt+
     &                 ak(k,62)*hocl*oh*dt+
     &                 ak(k,27)*hcl*oh*dt+
     &                 ak(k,40)*oh*ho2*dt+
     &                (ak(k,57)*o1d*ch4+ak(k,71)*h*ho2+aj(11)*ch2o)*dt+
     &                 ak(k,95)*oh*hbr*dt+
     &                 ak(k,72)*h*ho2*dt+
     &                 ak(k,55)*oh*ho2no2*dt+
     &                 ak(k,30)*oh*h2*m(k)*dt

                  h2oloss = (aj(38)+aj(40)+ak(k,39)*o1d)*h2o*dt

                  h2o = h2o + (h2oprod-h2oloss)

c...  Update Cl and Br source gas species (sds-11/06/02)
                  f11 = f11 - (aj(29)
     &                 +ak(k,83)*o1d )*f11 *dt 
                  f12 = f12 - (aj(30)
     &                 +ak(k,80)*o1d )*f12 *dt 
                  f113 = f113 - (aj(34)
     &                 +ak(k,99)*o1d )*f113 *dt 
                  f114 = f114 - (aj(53)
     &                 +ak(k,148)*o1d )*f114*dt
                  f115 = f115 - (aj(54)
     &                 +ak(k,149)*o1d )*f115*dt
                  hcfc22 = hcfc22 - (aj(33)
     &                 +ak(k,126)*o1d+ak(k,98)*oh)
     &                 *hcfc22 *dt 

#ifndef REDUCED
                  hcfc141b = hcfc141b - (aj(51)
     &                 +ak(k,127)*o1d+ak(k,129)*oh)
     &                 *hcfc141b *dt
                  hcfc142b = hcfc142b - (aj(52)
     &                 +ak(k,128)*o1d+ak(k,130)*oh)
     &                 *hcfc142b *dt 
#endif
                  ccl4 = ccl4 - (aj(31)
     &                 +ak(k,121)*o1d )*ccl4 *dt 
                  ch3ccl3 = ch3ccl3 - (aj(32)
     &                 +ak(k,75)*oh )*ch3ccl3 *dt 
                  ch3cl = ch3cl - (aj(35)
     &                 +ak(k,16)*oh )*ch3cl *dt 
                  ch3br = ch3br - (aj(20)
     &                 +ak(k,97)*oh+ak(k,120)*o1d)
     &                 *ch3br *dt 
                  h1301 = h1301 - (aj(21)
     &                 +ak(k,116)*o1d )*h1301 *dt 
                  h1211 = h1211 - (aj(22)
     &                 +ak(k,117)*o1d )*h1211 *dt 

#ifndef REDUCED
                  h1202 = h1202 - aj(49)*h1202*dt
                  h2402 = h2402 - aj(50)*h2402*dt
                  chbr3 = chbr3 - (aj(44)
     &                  +ak(k,131)*oh)*chbr3 *dt
                  ch2br2 = ch2br2 -(aj(45)
     &                  +ak(k,132)*oh)*ch2br2 *dt
                  ch2brcl = ch2brcl - (aj(46)
     &                  +ak(k,133)*oh)*ch2brcl *dt
                  chbrcl2 = chbrcl2 - (aj(47)
     &                  +ak(k,134)*oh)*chbrcl2 *dt
                  chbr2cl = chbr2cl - (aj(48)
     &                  +ak(k,135)*oh)*chbr2cl *dt
                  hfc23 = hfc23 - (ak(k,136)*o1d
     &                  +ak(k,142)*oh)*hfc23 *dt
                  hfc32 = hfc32 - (ak(k,137)*o1d
     &                  +ak(k,143)*oh)*hfc32 *dt
                  hfc125 = hfc125 - (ak(k,138)*o1d
     &                  +ak(k,144)*oh)*hfc125 *dt
                  hfc134a = hfc134a - (ak(k,139)*o1d
     &                  +ak(k,145)*oh)*hfc134a *dt
                  hfc143a = hfc143a - (ak(k,140)*o1d
     &                  +ak(k,146)*oh)*hfc143a *dt
                  hfc152a = hfc152a - (ak(k,141)*o1d
     &                  +ak(k,147)*oh)*hfc152a *dt
#endif


c...  Update inferred species
c...  ox species
                  o3 = ro3ox*ox 
                  o3p = ro3po3*o3
                  o1d = ro1do3*o3
                  box_ro3ox = o3/ox
c...  nox species
                  no2 = rno2nox*nox    
                  no = rnono2*no2
                  no3 = rno3no2*no2
                  n = rnno*no
c...  Chlorine species
                  clo= rcloclx*clx 

                  cl = rclclo*clo
                  oclo = rocloclo*clo

                  cl2o2 = rcl2o2clo*clo

c...  Bromine species
                  bro = rbrobrx*brx
                  if(bro .lt. 1.00e-20) bro = 1e-20
                  
                  brono2 = rbrnbrx*brx
     
                  br = brx-bro-brono2

                  brcl = 0.0
                  
c...  BEGIN NIGHTTIME CHEMISTRY. 
c...  Remember that the first call to
c...  SOLVER takes care of the dusk repartitioning
               else
                  
c...  Calculate NO3 and NO2 concentrations by assuming N and NO
c...  concentrations are zero, using transported [NOx], and using
c...  NO3 and NO2 partitioning from previous step. 
                  
                  if (no3.gt.0.) then  
                     no3 = (1./(1.+no2/no3))
     &                    *nox 
                     no2 = nox-no3
                  else 
                     no3 = 0. 
                     no2 = nox 
                  endif
                  
c...  Ox, O3, O3P.  Latest revision: 14 Oct 2004 from Gupta

                  IF(box_ro3ox .NE. -1.00) THEN
                     o3  = box_ro3ox*ox
                     o3p = ox-o3
		     alimit = ak(k,1)*m(k)*o2*o3p*dt*m(k)
                     IF(alimit >= o3p) THEN
                        o3 = ox
                        o3p = 1.00E-30
                        box_ro3ox = 1.00
                     ELSE
                        o3  = o3 + alimit
                        o3p = ox-o3
                        box_ro3ox = o3/ox
                     END IF
                  ELSE
                     o3 = ox   
                  END IF

                  pno3 = ak(k,10)*no2*o3
     &                 +ak(k,31)*n2o5*m(k)
                  
                  pbrono2 = ak(k,96)*bro*m(k)
     &                 *no2
                  
                  lbrono2het = ak(k,115)*brono2      
                  
                  pbrcl = ak(k,104)*bro*clo
                  
c     fxc brx conservation fix
                  
                  if((pbrono2+pbrcl)*hdt .gt. brx-
     &                 brono2) then  
                     tloss = pbrono2+pbrcl
                     
                     pbrono2 = (brx-brono2)/(hdt)*
     &                    pbrono2/tloss
                     pbrcl = (brx-brono2)/(hdt)*
     &                    pbrcl/tloss
                  endif
                  
                  if(brono2.lt.
     &                 (lbrono2het-pbrono2)*hdt)
     &                 lbrono2het = brono2/(hdt)+pbrono2
                  
                  lno3 = (ak(k,46)*m(k))*no3*
     &                 no2
                  
                  pno2 = ak(k,31)*n2o5*m(k)+
     &                 gcr(k) 
                  
                  pclono2 = ak(k,32)*m(k)*clo*
     &                 no2 
                  
                  pn2o5 = (ak(k,46)*m(k))*no3*
     &                 no2
                  
                  lno2 = ak(k,10)*o3*no2+
     &                 pn2o5 + pclono2 + pbrono2
                  
                  ln2o5 = (ak(k,31)*m(k)+ak(k,110)
     &                 +ak(k,113)*hcl)*n2o5 
                  
                  lclono2= (ak(k,111)*hcl+rain(k))*
     &                 clono2
                  
                  lclono2het = ak(k,112)*clono2 
                  
c...  The addition of heterogeneous rxn 113:
c...  N2O5 + HCL(aer) --> HNO3 + ClONO is not straightforward since
c...  ClONO is not carried in the CTM.  ClONO accumulates at night
c...  and  quickly photolyzes into Cl and NO2 during the morning. 
c     
c...  Therefore this reaction leads to the production of HNO3, Clx,
c...  and  NOx during the daytime. 
c...  At night it is unclear what to do with the ClONO. The reaction
c...  itself is only marginally important. What to do?
c...  1) Set reaction to zero (both day and night)?
c...  2) Assume ClONO  is stored as ClONO2?
c...  3) Assume ClONO is stored as Cl2 and HNO3 i.e. assume 
c...  N2O5+HCl --> 2HNO3 + 1/2Cl2  (at night)
c...  We chose option 3. 
c     
c...  PROBLEM: Counting oxygen atoms indicates that 5+0 .ne. 6+0.
c...  In other words this reaction produces an odd oxygen atom that
c...  we will not count in the PO3 subs. 
                  
                  phno3 = 2.*ak(k,110)*n2o5
     &                 +2.*ak(k,113)*n2o5*hcl
     &                 +lclono2het+ak(k,111)*hcl*clono2
     &                 +lbrono2het !Crk
                  
                  lhno3 = rain(k)*hno3 
                  
                  poclo = ak(k,102)*clo*bro
                  
                  rcl2o2clo = ak(k,107)*clo
     &                 /ak(k,108)
                  
                  rbrbro = (2.*ak(k,94)*bro 
     &                 +(ak(k,102)+ak(k,103))*clo
     &                 +ak(k,118)*o3p )
     &                 /(ak(k,91)*o3)
                  
c...  note  count clono as .5*cl2
                  pcl2 = (ak(k,111)*clono2*hcl
     &                 +0.5*ak(k,113)*n2o5*hcl
     &                 +ak(k,114)*hocl*hcl)
                  
                  lclx = pclono2+pbrcl
                  
                  lhcl = hcl*(ak(k,111)*clono2
     &                 +ak(k,113)*n2o5+ak(k,114)*hocl
     &                 +rain(k)) 
                  
                  lh2o2 = rain(k)*h2o2 
                  
                  lho2no2 = rain(k)*ho2no2 
                  
                  phocl = ak(k,112)*clono2
                  
                  lhocl = ak(k,114)*hocl*hcl
                  
c     .fxc... Constrain clx loss rate
                  if(lclx*hdt.gt.clx) then
                     tmplclx = lclx 
                     lclx = (clx*0.99)/(hdt)
                     rat = lclx / tmplclx 
                     
                     pclono2 = pclono2 * rat 
                     pbrcl = pbrcl * rat
                     
                     lno2 = (ak(k,10)*o3
     &                    +(ak(k,46)*m(k))*no3)
     &                    *no2+pclono2+pbrono2
                     
                  endif
                  
c...  Negative NOx mixing ratios have become a problem for conditions of
c...  low  NOx.  
c...  NOx loss is constrained to keep NOx at least zero after time step.
c...  To conserve nitrogen and chlorine species.  The production and loss
c...  of other species involved in the loss of NO2 must also be adjusted. 
c...  rat gives the ratio of lno2 after constraint to lno2 before con-
c...  straint. 
                  
                  if((no2+hdt*(pno2-lno2))
     &                 .lt.0.) then 
                     
                     rat = lno2
                     lno2 = no2 / (hdt) 
                     rat = lno2 / rat 
                     
                     pno3 = ak(k,31)*n2o5*m(k)
     &                    +rat*ak(k,10)*no2*o3
                     
                     lno3 = rat * lno3
                     
                     pn2o5 = rat * pn2o5 
                     pclono2 = rat*pclono2
                     
                     lclx = pclono2+pbrcl 
                     pbrono2 = rat*pbrono2        
                     
                  endif 
                  
c...  To conserve nitrogen and chlorine species.  The production and loss
c...  of other species involved in the loss of NO3 must also be adjusted. 
c...  rat gives the ratio of lno3 after constraint to lno3 before con-
c...  straint. 
                  
                  pnoxa = no3
     &                 +hdt*(pno3-lno3)  
                  if(pnoxa.lt.0.) then 
                     
                     rat = lno3 
                     lno3 = no3 / (hdt)  
                     rat = lno3 / rat 
                     
                     pn2o5 = rat * pn2o5       
                     lno2 = (ak(k,10)*o3+rat*
     &                    ak(k,46)*m(k)*no3)*no2
     &                    +pclono2+pbrono2
                     
                  endif 
                  
c...  Negative mixing ratios of HCl were encountered in well-processed
c...  air.   The following fix will reduce nighttime processing when 
c...  updated mixing ratio of HCl is .lt. zero
                  if(hcl.lt.hdt*lhcl) then 
                     
c...  store old lhcl
                     rat = lhcl
                     
c...  Set new lhcl to 0 since phcl=0 at night. 
                     lhcl = dim(hcl,1.0e6)/(hdt)   
                     if(lhcl*hdt .gt. hcl) 
     &                    write(6,*) 'lhcl > hcl in solverd2: ',k
                     
c...  Calculate ratio of new to old lhcl
                     rat = lhcl / rat   
                     
                     lclono2 = ( rain(k)
     &                    +rat*ak(k,111)*hcl )*clono2
                     
                     lhocl = rat*hocl*ak(k,114)*
     &                    hcl
                     
                     phno3 = 2.*( ak(k,110)
     &                    +rat*ak(k,113)*hcl )*n2o5
     &                    +rat*ak(k,111)*hcl*clono2
     &                    +lclono2het +lbrono2het !Crk   
                     
                     ln2o5 = (ak(k,31)*m(k)+ak(k,110)
     &                    +rat*ak(k,113)*hcl)*n2o5 
                     
                     pcl2 = rat*pcl2 
                     
                  endif 
                  
c...  UPDATE NIGHTTIME CONCENTRATIONS. 
                  
c...  N2O5+H2O -> 2HNO3 actually produces 1 Ox but we will
c...  not include it for now.           
                  
c...  note  that for nighttime, production and loss terms contain the het
c...  processes;these are not kept separate as during the daytime
                  
                  no2 = no2
     &                 +(pno2-lno2)*hdt 
                  
                  no3 = no3
     &                 +(pno3-lno3)*hdt  
                  
                  hno3 = hno3
     &                 +(phno3-lhno3)*hdt 
                  
                  n2o5 = n2o5+(pn2o5-ln2o5)*hdt
                  
                  nox = no3 + no2 
                  
                  clx = clx-lclx*hdt
                  
                  oclo = oclo+poclo*hdt

                  clono2 = clono2
     &                 +(pclono2-lclono2-lclono2het)*hdt 
                  
                  hcl = hcl - lhcl*hdt 
                  
                  hocl = hocl+(phocl-lhocl)*hdt
                  
c...  note  0.5 for time step and 0.5 b/c 1 per 2 clx or hcl
                  cl2 = cl2+pcl2*hdt
                  
c...  note: Partitioning is used to determine cl2o2 
c...  clx   and oclo are calulated using P and L s.   clo is calc.
c...  by solving clx = clo+2*cl2o2+oclo for clo. 
                  
                  clo = (clx-oclo)
     &                 /(1.+2.*rcl2o2clo)
c
                  if (clo.lt.0.) then 
                     clo = 0.
                     oclo = 0. 
                  endif         
                  
                  cl2o2 = rcl2o2clo*clo
                  
                  brono2 = brono2
     &                 +(pbrono2-lbrono2het)*hdt       
                  
                  brcl = brcl+pbrcl*hdt
                  
                  hobr = hobr+lbrono2het*hdt
                  
                  hbr = hbr*(1.-rain(k)*hdt)
                  
                  brx = brx - (lbrono2het+pbrcl)*hdt
                  
                  bro = (brx-brono2)/(1.+rbrbro)

c Prevent negative Bro MG/JEN 23 Feb 2004
                  if(bro .le. 1.00e-20) bro = 1.00e-20
                  
                  br = bro*rbrbro
                  
                  h2o2 = h2o2-lh2o2*hdt   
                  
                  ho2no2 = ho2no2-lho2no2*hdt
                  
c...  Update Ox - nighttime loss in troposphere
                  if(tropo(k)) 
     &                 ox = ox - lo3*dt
c     .sds... check for ox very small and constrain
                  if(ox < 1.0e-10*m(k)) then
                     ox = 1.0e-10*m(k)
                  endif
                  
               endif
      endif
      
c...  do loss of surface source gases in lower layers
      
      if(k .lt. k1Strat) then
c...  Update stratospheric Ox below "chemistry"
               ox = ox - lo3*dt
               o3 = ox
c...  OH reaction with source gases.  CH4 added Tue 16 Mar 2004 JEN.
c...  Qing Liang -- 3 Jun 2015
c...  Multiply OH by 2 to account for 24-hour averaged OH climatology 
c...  Add HCFCs, VSLSs, extend Halon photolysis to troposphere 
               if(daytime) then
                  ch3ccl3 = ch3ccl3 - (aj(32)
     &                 +ak(k,75)*2*oh )*ch3ccl3 *dt 
                  ch3cl = ch3cl - (aj(35)
     &                 +ak(k,16)*2*oh )*ch3cl *dt 
                  ch3br = ch3br - (aj(20)
     &                 +ak(k,97)*2*oh)*ch3br *dt 
                  ch4 = ch4 - (ak(k,14)*2*oh*ch4*dt)
                  hcfc22 = hcfc22 - (aj(33)
     &                 +ak(k,126)*o1d+ak(k,98)*2*oh)
     &                 *hcfc22 *dt
                  hcfc141b = hcfc141b - (aj(51)
     &                 +ak(k,127)*o1d+ak(k,129)*2*oh)
     &                 *hcfc141b *dt
                  hcfc142b = hcfc142b - (aj(52)
     &                 +ak(k,128)*o1d+ak(k,130)*2*oh)
     &                 *hcfc142b *dt
                  h1211 = h1211 - (aj(22)
     &                 +ak(k,117)*o1d )*h1211 *dt
                  h1202 = h1202 - aj(49)*h1202*dt
                  h2402 = h2402 - aj(50)*h2402*dt
                  chbr3 = chbr3 - (aj(44)
     &                  +ak(k,131)*2*oh)*chbr3 *dt
                  ch2br2 = ch2br2 -(aj(45)
     &                  +ak(k,132)*2*oh)*ch2br2 *dt
                  ch2brcl = ch2brcl - (aj(46)
     &                  +ak(k,133)*2*oh)*ch2brcl *dt
                  chbrcl2 = chbrcl2 - (aj(47)
     &                  +ak(k,134)*2*oh)*chbrcl2 *dt
                  chbr2cl = chbr2cl - (aj(48)
     &                  +ak(k,135)*2*oh)*chbr2cl *dt
                  hfc23 = hfc23 - (ak(k,136)*o1d
     &                  +ak(k,142)*2*oh)*hfc23 *dt
                  hfc32 = hfc32 - (ak(k,137)*o1d
     &                  +ak(k,143)*2*oh)*hfc32 *dt
                  hfc125 = hfc125 - (ak(k,138)*o1d
     &                  +ak(k,144)*2*oh)*hfc125 *dt
                  hfc134a = hfc134a - (ak(k,139)*o1d
     &                  +ak(k,145)*2*oh)*hfc134a *dt
                  hfc143a = hfc143a - (ak(k,140)*o1d
     &                  +ak(k,146)*2*oh)*hfc143a *dt
                  hfc152a = hfc152a - (ak(k,141)*o1d
     &                  +ak(k,147)*2*oh)*hfc152a *dt

               endif

      endif

      RETURN
      END SUBROUTINE solverd2
