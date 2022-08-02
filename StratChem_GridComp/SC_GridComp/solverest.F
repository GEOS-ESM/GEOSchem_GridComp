      SUBROUTINE solverest(dt,requiredDt,examineDt,k,k1Strat,levels,daytime,tropo,m,aj,ak,br,
     &         brcl,bre,bro,broe,brono2,brono2e,brx,brxe,ccl4,ch3ccl3,ch3cl,cl,cle,cl2,cl2o2,
     &    cl2o2e,clo,cloe,clono2,clono2e,clx,clxe,f11,f113,f12,f114,f115,h2o2,h2o2e,hbr,hbre,
     &                                hcfc22,hcfc141b,hcfc142b,h1211,ch2brcl,chbr2cl,chbrcl2,
     &               hcl,hcle,hno3,hno3e,ho2,ho2no2,ho2no2e,hobr,hobre,hocl,hocle,lbrono2het,
     &          lbrx,lclono2,lclono2het,lclx,lh2o2,lhbr,lhcl,lhno3,lho2no2,lhobr,lhocl,ln2o5,
     &        ln2o5het,lno2,lno3,lnox,lnoxa,lnoxcl,lnoxsq,lo3,n2o,n2o5,n2o5e,ne,no2,no2e,no3,
     &      no3e,no,noe,nox,noxe,o1d,o3,o3e,o3p,o3pe,oclo,ocloe,oh,ox,oxe,pbrcl,pbrono2,pbrx,
     &            pcl2,pclono2,pclx,ph2o2,phbr,phcl,pho2no2,phocl,phno3,phno3het,phobr,pn2o5,
     &                  pno2,pno3,pnox,pnoxa,pnoxcl,po3,poclo,rain,rbrnbrx,rbrobrx,rcl2o2clo,
     &            rclclo,rcloclx,rnno,rnono2,rno2nox,rno3no2,ro3ox,ro3po3,rocloclo,box_ro3ox,
     &                                                   levels_cal,numreacs,numphoto,gcr,o2)

c     This is the FIRST solution.  Daytime, this calculation is a prediction;
c     these arrays are named for the constituent with an "e" added; these
c     are 2d arrays.  P&L used here for the estimate will be averaged with
c     values calculated using the estimates.  Nighttime, this is an explicit 
c     calculation, but with dtnight= .5*dtday.  
      
c     Latest modification  14 July 1994 
C     rk 951107 Updated to include CHx chemistry and new Br scheme
      
c...  Currently no production or loss of ox at night (rxn 19,112)

      IMPLICIT NONE
      INTEGER, PARAMETER :: DBL = KIND(0.00D+00)

      INTEGER, INTENT(IN) :: examineDt,k,k1Strat,levels
      INTEGER, INTENT(IN) :: levels_cal,numreacs,numphoto
      LOGICAL, INTENT(IN) :: daytime,tropo(levels)
      REAL, INTENT(IN) :: dt,m(levels),rain(levels),o2

      REAL, INTENT(IN) :: gcr(levels)

      REAL(KIND=DBL), INTENT(IN) :: aj(numphoto),ak(levels,numreacs)
      REAL(KIND=DBL), INTENT(IN) :: ccl4,ch3ccl3,ch3cl,f11,f113,f12,f114,f115
      REAL(KIND=DBL), INTENT(IN) :: h1211,hcfc22,hcfc141b,hcfc142b
      REAL(KIND=DBL), INTENT(IN) :: ch2brcl,chbr2cl,chbrcl2
      REAL(KIND=DBL), INTENT(IN) :: ho2,lbrx,lhbr,lhobr,lo3,n2o,o1d,o3p
      REAL(KIND=DBL), INTENT(IN) :: oh,ox,pbrx,ph2o2,phbr,phobr,pnoxcl,po3
      REAL(KIND=DBL), INTENT(IN) :: rbrobrx,rclclo,rcloclx,rnno,rnono2,rno2nox
      REAL(KIND=DBL), INTENT(IN) :: rno3no2,ro3ox,ro3po3,rocloclo

      REAL, INTENT(INOUT) :: requiredDt

      REAL(KIND=DBL), INTENT(INOUT) :: br,brcl,bre,bro,broe,brono2,brono2e,brx
      REAL(KIND=DBL), INTENT(INOUT) :: brxe,cl,cle,cl2,cl2o2,cl2o2e,clo
      REAL(KIND=DBL), INTENT(INOUT) :: cloe,clono2,clono2e,clx,clxe,h2o2,h2o2e
      REAL(KIND=DBL), INTENT(INOUT) :: hbr,hbre,hcl,hcle,hno3,hno3e,ho2no2
      REAL(KIND=DBL), INTENT(INOUT) :: ho2no2e,hobr,hobre,hocl,hocle
      REAL(KIND=DBL), INTENT(INOUT) :: lbrono2het,lclono2,lclono2het
      REAL(KIND=DBL), INTENT(INOUT) :: lclx,lh2o2,lhcl,lhno3,lho2no2,lhocl
      REAL(KIND=DBL), INTENT(INOUT) :: ln2o5,ln2o5het,lno2,lno3,lnox
      REAL(KIND=DBL), INTENT(INOUT) :: lnoxa,lnoxcl,lnoxsq,n2o5,n2o5e
      REAL(KIND=DBL), INTENT(INOUT) :: ne,no2,no2e,no3,no3e,no,noe,nox
      REAL(KIND=DBL), INTENT(INOUT) :: noxe,o3,o3e,o3pe,oclo,ocloe,oxe,pbrcl
      REAL(KIND=DBL), INTENT(INOUT) :: pbrono2,pcl2,pclono2,pclx,phcl,pho2no2
      REAL(KIND=DBL), INTENT(INOUT) :: phocl,phno3
      REAL(KIND=DBL), INTENT(INOUT) :: phno3het,pn2o5,pno2,pno3,pnox,pnoxa
      REAL(KIND=DBL), INTENT(INOUT) :: poclo,rbrnbrx,rcl2o2clo,box_ro3ox

      LOGICAL :: runVariableDt

      INTEGER :: itercount

      REAL(KIND=DBL) :: alimit,brntmp,cnoxbr,dbrn,dbrn0,rat,rbrbro
      REAL(KIND=DBL) :: sav,testclx,testhcl,testnox,tloss,tmplclx

      REAL :: hdt,minAllowDt,shortDt(3)
      
      hdt = 0.5*dt

      if(k .ge. k1Strat .and. k .le. levels_cal) then

               if(daytime) then
c     fxc... put in fix to keep transported brono2 <= brx
                  if(brono2.gt.brx) brono2 = brx
                  
c...  Add Cl2 and brcl back into Clx and Brx for the daytime calculation
c...  Since these are transported separately, this addition should
c...  come here, so anything transported out of the polar night will be
c...  added in correctly.
c...  O2 mixing ratio is specified.
      
                  clx = clx+2.*cl2+brcl
                  brx = brx+brcl
                  cl2 = 0. 
                  brcl = 0.     
                  
                  brxe = brx+(pbrx-lbrx)*dt

! Prevent daytime negative Brx when forcing term is within truncation
! error of input Brx. JEN 23 Feb 2004
                  if (brxe .lt. 1.00e-20) brxe=1.00e-20

                  dbrn = (rbrnbrx*brxe-brono2)/dt
c...  gas phase
                  dbrn0 = dbrn+lbrono2het
c...  [nox] lost to BrONO2
                  cnoxbr = max(dbrn0,0.)
                  
c...  NOX fix for denitrified conditions.  Scale losses to [NOx]/dt.
C...  Assume pnox is enough to produce BrONO2 w/o limiting 
                  pnoxa = (pnox+pnoxcl)*dt 
                  lnoxa = (lnox+lnoxcl+lnoxsq)*dt+cnoxbr*dt
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

       IF((testclx .LE. 0.00 .OR. testhcl .LE. 0.00 .OR. 
     &     testnox .LE. 0.00) .AND. examineDt .LE. 3) THEN
     
        shortDt(1:3)=1.00e+10
        IF(testnox .LE. 0.00) shortDt(1)=nox/ABS(pnox+pnoxcl-
     &    	        		 (lnox+lnoxcl+lnoxsq+cnoxbr))
        IF(testclx .LE. 0.00) shortDt(2)=clx/ABS(pclx-lclx)
        IF(testhcl .LE. 0.00) shortDt(3)=hcl/ABS(phcl-lhcl)
        requiredDt=MAX(MINVAL(shortDt),minAllowDt)
     	RETURN

       ENDIF

      ELSE
                  
       do while ((testclx < 0. .or. testhcl < 0. .or. testnox < 0)
     &      .and. itercount < 5)
                     
                     if(testnox < 0) then
                        
                        rat = nox / lnoxa  
                        lnox = rat*lnox    
                        lnoxcl = rat*lnoxcl
                        lnoxsq= rat*lnoxsq 
                        
                        lclx=lclx-(1.-rat)*pclono2
                        phno3 = rat*phno3
                        pho2no2 = rat*pho2no2
                        pclono2 = rat*pclono2
                        pn2o5 = rat*pn2o5
c...  NOx is lost to BrONO2
                        if(dbrn0.gt.0.) then
                           brntmp = rbrnbrx
                           dbrn0 = rat*dbrn0
                           dbrn = dbrn0-lbrono2het
                           rbrnbrx = (dbrn*dt+brono2)/brxe

       if(brntmp .le. 0.) print 345, k,brntmp,dbrn,dbrn0,
     &                    rbrnbrx,lbrono2het,brxe
  345 format(' ','SOLVEREST:(',i4,') brntmp=',1pe10.3,' dbrn=',
     &       1pe10.3,' dbrn0=',1pe10.3,' rbrnbrx=',1pe10.3,
     &       ' lbrono2het=',1pe10.3,' brxe=',1pe10.3)

                           lbrono2het = rbrnbrx/brntmp
                        endif
c...  lclx changed above, recalculate testclx
                        testclx = clx+dt*(pclx-lclx)
                     endif 
                     
                     if(testclx < 0) then
                        sav = lclx
c...  lclx = clx/dt  
                        lclx = clx/dt*
     &                       (1.-EXP(-lclx*dt/clx))
                        rat = lclx / sav  
                        pclono2 = rat*pclono2
                        lnoxcl = pclono2
c...  taken out at some point            +pbrono2
                        phocl = rat*ak(k,6)*ho2
     &                       *rcloclx*clx+(ak(k,13)*oh
     &                       +ak(k,112))*clono2  
                        phcl = rat*phcl 
c...  phcl changed, recalc testhcl
                        testhcl = hcl+dt*(phcl-lhcl)
                     endif
                     
                     if(testhcl < 0) then
c...  store old lhcl
                        sav = lhcl
                        lhcl = hcl/dt*
     &                       (1.-exp(-lhcl*dt/hcl))
                        rat = lhcl / sav  
                        lclono2het = (rat*ak(k,111)*hcl
     &                       +ak(k,112))*clono2
                        phno3het = 2.*ak(k,110)*n2o5
     &                       +lclono2het+lbrono2het
     &                       +rat*ak(k,113)*n2o5*hcl 
                        pnox = 2.*ln2o5+lhno3
     &                       +lho2no2
     &                       +rat*ak(k,113)*n2o5*hcl
     &                       +gcr(k)+2.*ak(k,45)*n2o*o1d 
                        ln2o5het = ak(k,110)*n2o5
     &                       +rat*ak(k,113)*n2o5*hcl 
                        lhocl = (ak(k,62)*oh
     &                       +aj(23)
     &                       +rat*ak(k,114)*hcl)*hocl
                        pclx=lhcl-rat*rain(k)*hcl
     &                       +lhocl
c...  was clxsource from 2d model
     &                       +ak(k,16)*ch3cl*oh
     &                       +3*ak(k,75)*ch3ccl3*oh
     &                       +  ak(k,98)*hcfc22*oh

#ifndef REDUCED
     &                       +2*ak(k,129)*hcfc141b*oh
     &                       +  ak(k,130)*hcfc142b*oh    
     &                       +  ak(k,133)*ch2brcl*oh
     &                       +2*ak(k,134)*chbrcl2*oh
     &                       +  ak(k,135)*chbr2cl*oh
#endif

     &                       +2*ak(k,80)*f12*o1d
     &                       +3*ak(k,83)*f11*o1d

#ifndef REDUCED
     &                       +2*ak(k,148)*f114*o1d
     &                       +  ak(k,149)*f115*o1d
#endif

     &                       +  ak(k,126)*hcfc22*o1d

#ifndef REDUCED
     &                       +2*ak(k,127)*hcfc141b*o1d
     &                       +  ak(k,128)*hcfc142b*o1d
#endif

     &                       +  ak(k,117)*h1211*o1d
     &                       +3*ak(k,99)*f113*o1d
     &                       +4*ak(k,121)*ccl4*o1d
     &                       +3*aj(29)*f11
     &                       +2*aj(30)*f12
     &                       +4*aj(31)*ccl4
     &                       +3*aj(32)*ch3ccl3
     &                       +  aj(33)*hcfc22

#ifndef REDUCED
     &                       +2*aj(51)*hcfc141b
     &                       +  aj(52)*hcfc142b
#endif

     &                       +3*aj(34)*f113

#ifndef REDUCED
     &                       +2*aj(53)*f114
     &                       +  aj(54)*f115
#endif

     &                       +  aj(35)*ch3cl
     &                       +  aj(22)*h1211

#ifndef REDUCED
     &                       +  aj(46)*ch2brcl
     &                       +2*aj(47)*chbrcl2 
     &                       +  aj(48)*chbr2cl 
#endif

c...  end clxsource
     &                       +(ak(k,35)*o3p+aj(15)
     &                       +rat*ak(k,111)*hcl)*clono2
                     endif 
                     
          pnoxa = (pnox+pnoxcl )*dt 
          lnoxa = (lnox+lnoxcl
     &         +lnoxsq )*dt +cnoxbr*dt
          testnox = nox+pnoxa-lnoxa
          testclx = clx+dt*(pclx-lclx)
          testhcl = hcl+dt*(phcl-lhcl)
          itercount = itercount+1

       enddo
      END IF

c...  Update daytime transported species estimates
                  
                  clxe = clx+(pclx-lclx)*dt
                  
                  hbre = hbr+(phbr-lhbr)*dt
                  
                  hobre = hobr+(phobr-lhobr)*dt
                  clono2e = clono2
     &                 +(pclono2-(lclono2+lclono2het))*dt
                  
c     fxc check for clono2 < 0
                  if(clono2e < 0.) then
                     clono2e = 1.0e-15*m(k)
                  endif
c     fxc end check for clono2 < 0
                  
                  hcle = hcl+(phcl-lhcl)*dt 
c     .sds... check for hcl < 0 and set small
                  if(hcle < 0.) then
                     hcle = 1.0e-15*m(k)
                  endif
c     .sds.
                  
                  hocle = hocl+(phocl-lhocl)*dt
                  
                  n2o5e = n2o5+(pn2o5
     &                 -(ln2o5+ln2o5het))*dt
                  
c     .sds... check for n2o5 < 0 and set small
                  if(n2o5e < 0.) then
                     n2o5e = 1.0e-15*m(k)
                  endif
c     .sds.
                  hno3e= hno3 
     &                 +(phno3+phno3het-lhno3)*dt
                  
                  h2o2e = h2o2+(ph2o2-lh2o2)*dt
                  
                  ho2no2e = ho2no2
     &                 +(pho2no2-lho2no2)*dt 
                  
c...  Update Ox and NOx
                  if(tropo(k)) then
                   oxe = ox - lo3*dt
                  else
                   oxe = ox + (po3-lo3)*dt
                  endif

c     .sds... check for ox very small and constrain
                  if(oxe < 1.0e-10*m(k)) then
                     oxe = 1.0e-10*m(k)
                  endif
                  
                  noxe = nox
     &                 +(pnox-dbrn0+pnoxcl
     &                 -(lnox+lnoxcl+lnoxsq))*dt 
c     .sds... check for nox < 0 and set small
                  if(noxe < 0.) then
                     noxe = 1.0e-15*m(k)
                  endif
                  
c...  Update inferred species
                  o3e = ro3ox * oxe
                  o3pe = ro3po3 * o3e
                  
                  no2e = rno2nox * noxe 
                  noe = rnono2 * no2e
                  no3e = rno3no2 * no2e
                  ne = rnno * noe
                  
                  cloe = rcloclx*clxe
c                  cle = rclclo*cloe
c                  ocloe = rocloclo*cloe
                  cl2o2e = rcl2o2clo*cloe 
                  
                  broe = rbrobrx*brxe

c                  brono2e = rbrnbrx*brxe
c                  bre = brxe-broe-brono2e

c...  Don't worry about photolyzing BrCl; this process is so rapid, just
c...  set BrCl to zero above.
                  
c...  BEGIN NIGHTIME CHEMISTRY
                  
               else
c...  NOx is transported; NO3 is approximated from its previous
c...  value at this grid point
                  if(nox.gt.no3) then 
                     no2 = nox-no3
                  else 
                     no2 = 0. 
                     no3 = nox 
                  endif 
                  
c...  Assume all Cl + ClO goes to ClO
                  clo = cl + clo
                  
                  cl = 0. 
                  no = 0. 
                  
c...  Some species are transported as individuals and in families.
c...  In the daytime, the transport is ignored and they are calculated
c...  from photochemical equilibrium.  In the nighttime, they are
c...  subtracted from the Whole Family and then the remainder of the
c...  family is calculated from chemical production and loss.
c...  These are OClO and BrONO2.
                  
c...  Tie bromine species to transported brx. Retain old ratios. 
c...  Transport BrCl separately    
                  
                  
c     fxc put in fix to keep transported brono2 <= brx
                  
                  if(brono2.gt.brx) brono2 = brx
                  
c...  Partitioning from previous step and NOX after transport
c...  is used to infer no3 and no2 concentrations.  
c...  no3 after transport and BEFORE chemical step
                  
c...  approximate NO2 and NO3 partitioning as for previous timestep
                  
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
                     o3e  = box_ro3ox*ox
                     o3pe = ox-o3e
                     alimit = ak(k,1)*m(k)*o2*o3pe*dt*m(k)
                     IF(alimit >= o3pe) THEN
                        o3e = ox
                        o3pe = 1.00E-30
                        box_ro3ox = 1.0
                     ELSE
                        o3e  = o3e + alimit
                        o3pe = ox-o3e
                        box_ro3ox = o3e/ox
                     END IF
                  ELSE
                     o3e = ox
                  END IF

                  lh2o2 = rain(k)*h2o2 
                  
                  lho2no2 = rain(k)*ho2no2 
                  
                  pno3 = ak(k,10)*no2*o3+
     &                 ak(k,31)*n2o5*m(k)
                  
                  lno3 = (ak(k,46)*m(k))*no2*
     &                 no3 
                  
                  pno2 = ak(k,31)*n2o5*m(k) + 
     &                 gcr(k) 
                  
                  pbrono2 = ak(k,96)*bro*m(k)
     &                 *no2

                  lbrono2het = ak(k,115)*brono2

                  pbrcl = ak(k,104)*bro*clo 
                  
c...  solve brx conservation without messing up noy conservation
c...  (limit BrO loss since BrONO2 production is very rapid if NO2 is high)
                  
                  if((pbrono2+pbrcl)*hdt .gt. 
     &                 brx-brono2) then
                     
                     tloss = pbrono2+pbrcl

                     pbrono2 = (brx-brono2)/(hdt)*
     &                    pbrono2/tloss
                     pbrcl = (brx-brono2)/(hdt)*pbrcl/tloss

                  endif
                  
                  if(brono2.lt.
     &                 (lbrono2het-pbrono2)*hdt)
     &                 lbrono2het = brono2/(hdt)+pbrono2
                  
                  pn2o5 = (ak(k,46)*m(k))*no2*
     &                 no3
                  
                  ln2o5 = (ak(k,31)*m(k)+ak(k,110)
     &                 +ak(k,113)*hcl )*n2o5  
                  
                  pclono2 = ak(k,32)*m(k)*clo*
     &                 no2 
                  
                  poclo = ak(k,102)*clo*bro 
                  
                  lclono2 = (ak(k,111)*hcl+rain(k))*
     &                 clono2
                  
                  lclono2het = ak(k,112)*clono2 
                  
                  lhcl = hcl*(ak(k,111)*clono2
     &                 +ak(k,113)*n2o5+ak(k,114)*hocl
     &                 +rain(k)) 
                  
                  phocl = ak(k,112)*clono2 
                  
                  lhocl = ak(k,114)*hocl*hcl
                  
c...  The addition of heterogeneous rxn 113:
c...  N2O5 + HCL(aer) --> HNO3 + ClONO is not straightforward since
c...  ClONO is not carried in the CTM.  ClONO accumulates at night
c...  and quickly photolyzes into Cl and NO2 during the morning. 
c     
c...  Therefore this reaction leads to the production of HNO3, Clx,
c...  and NOx during the daytime. 
c...  At night it is unclear what to do with the ClONO. The reaction
c...  itself is only marginally important. What to do?
c...  &) Set reaction to zero (both day and night)?
c...  2) Assume ClONO  is stored as ClONO2?
c...  3) Assume ClONO is stored as Cl2 and HNO3 i.e. assume 
c...  N2O5+HCl --> 2HNO3 + 1/2Cl2  (at night)
c...  We chose option 3. 
c     
c...  PROBLEM: Counting oxygen atoms indicates that 5+0 .ne. 6+0.
c...  This minor source of odd oxygen is not counted in the PO3 
c...  subs
                  
                  phno3 = 2.*( ak(k,110)
     &                 +ak(k,113)*hcl )*n2o5
     &                 +ak(k,111)*hcl*clono2
     &                 +lclono2het+lbrono2het
                  
                  lhno3 = rain(k)*hno3 
                  
c...  Assume Cl2O2 and ClO are equilibrated
                  
                  rcl2o2clo = ak(k,107)*clo/ak(k,108)
                  
                  lclx = pclono2+pbrcl
                  
c...  count clono as .5*cl2
                  
                  pcl2 = ( ak(k,111)*clono2
     &                 +0.5*ak(k,113)*n2o5
     &                 +ak(k,114)*hocl )*hcl
                  
                  lno2 = ( ak(k,10)*o3
     &                 +ak(k,46)*m(k)*no3 )*no2
     &                 +pclono2 +pbrono2
                  
c     fxc Constrain clx loss rate
                  
                  if(lclx*hdt .gt. clx) then
                     tmplclx = lclx  
c     lclx = clx/(hdt) 
                     lclx = (clx*0.99)/(hdt)
                     rat = lclx / tmplclx
                     
                     pclono2 = pclono2*rat
                     pbrcl = pbrcl*rat
                     
                     lno2 = ( ak(k,10)*o3
     &                    +ak(k,46)*m(k)*no3 )*no2
     &                    +pclono2 +pbrono2
                     
                  endif
                  
                  rbrbro = (2.*ak(k,94)*bro
     &                 +(ak(k,102)+ak(k,103))*clo
     &                 +ak(k,118)*o3p )
     &                 / (ak(k,91)*o3)
                  
c...  Negative NOx mixing ratios have become a problem for conditions of
c...  low NOx.  
c...  NOx loss is constrained to keep NOx at least zero after time step.
c...  To conserve nitrogen and chlorine species.  The production and loss
c...  of other species involved in the loss of NO2 must also be adjusted. 
c...  rat gives the ratio of lno2 after constraint to lno2 before con-
c...  straint. 
                  
                  if((no2+hdt*(pno2-lno2)).lt.0.)
     &                 then 
                     rat = lno2
                     lno2 = no2 / (hdt) 
                     rat = lno2 / rat 
                     
                     pno3 = ak(k,31)*n2o5*m(k)
     &                    +rat*ak(k,10)*no2*o3
                     lno3 = rat*lno3
                     pn2o5 = rat*pn2o5 
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
                  
                  if (pnoxa.lt.0.) then
                     rat = lno3
                     lno3 = no3 / (hdt) 
                     rat = lno3 / rat 
                     
                     pn2o5 = rat * pn2o5       
                     lno2 = ( ak(k,10)*o3
     &                    +rat*ak(k,46)*m(k)*no3 )
     &                    *no2
     &                    +pclono2 +pbrono2
                  endif 
                  
c...  Negative mixing ratios of HCl were encountered in well-processed
c...  air.  The following fix will reduce nighttime processing when 
c...  updated mixing ratio of HCl is .lt. zero. 
                  if(hcl-hdt*(lhcl).lt.0.) then 
                     
c...  store old lhcl
                     rat = lhcl
                     
c...  Set new lhcl to 0 since phcl=0 at night. 
                     lhcl = dim(hcl,1.0e6)/hdt 
                     if(lhcl*hdt > hcl) then
                        print *,'SOLVEREST: lhcl > hcl: ',k,hcl
                     endif
                     
c...  Calculate ratio of new to old lhcl
                     rat = lhcl / rat   
                     
                     lclono2 = rat*ak(k,111)*hcl*
     &                    clono2 + rain(k)*clono2 
                     
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
                  
                  no2 = no2+(pno2-lno2)*hdt
                  
                  no3 = no3+(pno3-lno3)*hdt
                  
                  hno3 = hno3+(phno3-lhno3)*hdt
                  
                  n2o5 = n2o5+(pn2o5-ln2o5)*hdt
                  
                  nox = no3+no2 
                  
c...  Update Ox - tropospheric loss only at night
                  if(tropo(k)) 
     &                 oxe = ox - lo3*hdt
                  
c...  0.25 is used because only 1 cl2 for every hcl and clono2
c...  from hcl+clono2. 
                  cl2 = cl2+hdt*pcl2 
                  
                  clx = clx-lclx*hdt
                  
                  oclo = oclo+poclo*hdt
                  
c     This trap was introduced to make the code more stable in 
c     the middle latitude upper trop/lower strat where Clx is very small
c     and NO2 is large enough that ClONO2 formation is expected at night.
                  
                  clono2 = clono2 
     &                 +(pclono2-lclono2-lclono2het)*hdt 
                  
                  hcl = hcl-lhcl*hdt 
                  
                  hocl = hocl+hdt*(phocl-lhocl)
     &                 
                  
                  clo = (clx-oclo)
     &                 /(1.+2.*rcl2o2clo)

                  if (clx-oclo .lt.0.) then 
                     clo = 0.
                     oclo = clx 
                  endif  
                  cl2o2 = rcl2o2clo*clo

                  brono2 = brono2
     &                 +(pbrono2-lbrono2het)*hdt !Crk
                  
                  hobr = hobr+lbrono2het*hdt
                  
                  hbr = hbr*(1.-rain(k)*hdt)
                  
                  brcl = brcl+pbrcl*hdt
                  
                  brx = brx -(lbrono2het+pbrcl)*hdt
                  bro = (brx-brono2)/(1.+rbrbro)

c Prevent small negative BrO MG/JEN 23 Feb 2004
                  if(bro.lt. 1.00e-20) bro = 1.00e-20
                  
                  br = bro*rbrbro

                  h2o2 = h2o2-lh2o2*hdt 
                  
                  ho2no2 = ho2no2-lho2no2*hdt 
                  
               endif 

      endif
      RETURN
      END SUBROUTINE solverest
