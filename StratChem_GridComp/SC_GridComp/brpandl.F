      SUBROUTINE brpandl(dt,k,k1Strat,levels,daytime,ak,aj,brono2,brx,ch2o,ch3br,h1301,h1211,
     &                         h1202,h2402,chbr3,ch2br2,ch2brcl,chbrcl2,chbr2cl,hbr,ho2,hobr,
     &                                lbrono2het,lbrx,lhbr,lhobr,o1d,oh,pbrx,phbr,phobr,rain,
     &                                       rbrnbrx,rbrobrx,levels_cal,numreacs,numphoto,o2)

      IMPLICIT NONE
      INTEGER, PARAMETER :: DBL = KIND(0.00D+00)

      INTEGER, INTENT(IN) :: k,k1Strat,levels,levels_cal,numreacs,numphoto
      LOGICAL, INTENT(IN) :: daytime
      REAL, INTENT(IN) :: dt,o2,rain(levels)

      REAL(KIND=DBL), INTENT(IN) :: aj(numphoto),ak(levels,numreacs)
      REAL(KIND=DBL), INTENT(IN) :: brono2,brx,ch2o,ch3br
      REAL(KIND=DBL), INTENT(IN) :: h1301,h1211,h1202,h2402
      REAL(KIND=DBL), INTENT(IN) :: chbr3,ch2br2,ch2brcl,chbrcl2,chbr2cl
      REAL(KIND=DBL), INTENT(IN) :: hbr,ho2,hobr,o1d,oh,rbrnbrx,rbrobrx

      REAL(KIND=DBL), INTENT(INOUT) :: lbrono2het,lbrx,lhbr
      REAL(KIND=DBL), INTENT(INOUT) :: lhobr,pbrx,phbr,phobr

      REAL(KIND=DBL) :: tbr,tbro,tbrono2,tmp,tmp2

c  First subroutine to calculate the production and loss terms for
c   bry species

Crk  951005 adapted from ncbpandl and part subroutines.
C In this version the Bry species: BrONO2(night), HOBr, HBr, BrCl(night), 
C and Brx (=BrO + Br + BrONO2(day)), are transported and 
C solved explicitly from P and L.  BrCl is 0 in daylight.
C BrO/Brx and BrONO2/Brx are set in part.
c...  O2 mixing ratio is specified.
Crk  951107  adapted to full 3-D
C  20180119  Took out the rain term from lhbr and pbrx  (MEM per QL)

      if(k .ge. k1Strat .and. k .le. levels_cal) then

               if(daytime) then
                  
                  tbrono2 = rbrnbrx*brx
                  tbro = rbrobrx*brx
                  tbr = brx-tbro-tbrono2
                  
c...  brono2
                  lbrono2het = ak(k,115)*tbrono2
c
                  if(brono2.lt.lbrono2het*dt)
     &                 lbrono2het = brono2/dt
                  
c...  hbr
                  phbr = ( ak(k,92)*ho2
     &                 +ak(k,69)*ch2o )*tbr
                  
c                 lhbr = ( ak(k,95)*oh
c    &                 +rain(k) )*hbr
                  lhbr = ak(k,95)*oh*hbr
                  
                  if(hbr.lt.(lhbr-phbr)*dt) 
     &                 lhbr = hbr/dt+phbr
                  
c...  hobr
                  phobr = ak(k,109)*ho2*tbro
     &                 +lbrono2het
                  
c     .sds... new formulation for lhobr
                  lhobr = (hobr-phobr/aj(25))
     &                 *(1.-exp(-aj(25)*dt))/dt + phobr
                  
c...  brx
c                  pbrx = lhobr+(lhbr
c    &                  -rain(k)*hbr)
                   pbrx = lhobr+lhbr
c...  rest was brysource
     &                 +ak(k,117)*h1211*o1d
     &                 +ak(k,116)*h1301*o1d
     &                 +(ak(k,97)*oh+ak(k,120)*o1d)*ch3br
     &                 +  aj(20)*ch3br
     &                 +  aj(21)*h1301
     &                 +  aj(22)*h1211

#ifndef REDUCED
     &                 +  aj(49)*h1202
     &                 +2*aj(50)*h2402
     &                 +3*(aj(44)+ak(k,131)*oh)*chbr3
     &                 +2*(aj(45)+ak(k,132)*oh)*ch2br2
     &                 +  (aj(46)+ak(k,133)*oh)*ch2brcl
     &                 +  (aj(47)+ak(k,134)*oh)*chbrcl2
     &                 +2*(aj(48)+ak(k,135)*oh)*chbr2cl 
#endif
                  
                  lbrx = phobr +phbr
                  if(brx.lt.(lbrx-pbrx)*dt) then
                     tmp = max(lbrx,1.0)
                     lbrx = brx/dt+pbrx
                     tmp2 = lbrx/tmp
                     phobr = phobr*tmp2
                     lbrono2het = lbrono2het*tmp2
                     phbr = phbr*tmp2
                  endif
               endif

      endif

      RETURN
      END SUBROUTINE brpandl
