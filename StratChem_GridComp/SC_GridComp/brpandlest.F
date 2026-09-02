      SUBROUTINE brpandlest(dt,k,k1Strat,levels,daytime,ak,aj,brono2,brx,brxe,ch2o,ch3br,h1301,h1211,
     &                          h1202,h2402,chbr3,ch2br2,ch2brcl,chbrcl2,chbr2cl,hbre,ho2,hobr,hobre,
     &                                        lbrono2het,lbrx,lhbr,lhobr,o3e,oh,pbrx,phbr,phobr,rain,
     &                                        rbrnbrx,rbrobrx,ro1do3,levels_cal,numreacs,numphoto,o2)

      IMPLICIT NONE
      INTEGER, PARAMETER :: DBL = KIND(0.00D+00)

      INTEGER, INTENT(IN) :: k,k1Strat,levels,levels_cal,numreacs,numphoto
      LOGICAL, INTENT(IN) :: daytime
      REAL, INTENT(IN) :: dt,o2,rain(levels)

      REAL(KIND=DBL), INTENT(IN) :: aj(numphoto),ak(levels,numreacs)
      REAL(KIND=DBL), INTENT(IN) :: brono2,brx,brxe,ch2o,ch3br
      REAL(KIND=DBL), INTENT(IN) :: h1301,h1211,h1202,h2402
      REAL(KIND=DBL), INTENT(IN) :: chbr3,ch2br2,ch2brcl,chbrcl2,chbr2cl  
      REAL(KIND=DBL), INTENT(IN) :: hbre,ho2,hobr,hobre,o3e,oh
      REAL(KIND=DBL), INTENT(IN) :: rbrnbrx,rbrobrx,ro1do3
 
      REAL(KIND=DBL), INTENT(INOUT) :: lbrx,lbrono2het,lhbr,lhobr
      REAL(KIND=DBL), INTENT(INOUT) :: pbrx,phbr,phobr

      REAL(KIND=DBL) :: adj,o1de,tbre,tbroe,tbrono2e,tmp,tmp2

c  Subroutine to calculate final production and loss terms for 
c  Bry species;  the final value is the AVERAGE of a production and
c  loss term using the estimated concentrations and the production and loss
c  terms calculated using the previous time step concentrations.
c  O2 mixing ratio is specified.

      if(k .ge. k1Strat .and. k .le. levels_cal) then

          if(daytime) then 

            tbrono2e = rbrnbrx*brxe
            tbroe = rbrobrx*brxe
            tbre = brxe-tbroe-tbrono2e
            o1de = ro1do3*o3e

c det prod and loss of hobr, hbr, brx, and brono2 het
c AVERAGE for this solution scheme

c  brono2 must be handled this way to prevent instability in high sun

            lbrono2het = .5*( lbrono2het
     &       +ak(k,115)*tbrono2e ) 
     
cfxc replace the following fix with one using brono2 instead
cfxc of brono2e to be consistent with the way hobr is treated.
c         if (brono2e.LT.lbrono2het*dt) lbrono2het
c    &      = brono2e/dt
            if (brono2.lt.lbrono2het*dt) 
     &         lbrono2het = brono2/dt
     
            phbr = .5*( phbr
     &       +( ak(k,92)*ho2
     &       +ak(k,69)*ch2o )*tbre )

            lhbr = .5*(lhbr
     &       + ak(k,95)*oh*hbre)        
c     &       +( ak(k,95)*oh
c     &       +rain(k) )*hbre)

            phobr = .5*( phobr
     &       +ak(k,109)*ho2*tbroe
     &       +ak(k,115)*tbrono2e )
     &       
            lhobr = .5*(lhobr+aj(25)*hobre)
            adj = 0.
c... adjustment to loss
            if (hobr.lt.(lhobr-phobr)*dt) then
              tmp = lhobr
c.sds changed when changed lhobr calc in brpandl on 1999-08-06
              lhobr = (hobr-phobr/aj(25))
     &          *(1.-exp(-aj(25)*dt))/dt + phobr
              adj = tmp-lhobr
            endif

            pbrx = .5*( pbrx
     &       +aj(25)*hobre
     &       +ak(k,95)*oh*hbre 
c... rest was brysource (except "adj")
     &       +ak(k,117)*h1211*o1de
     &       +ak(k,116)*h1301*o1de
     &       +(ak(k,97)*oh+ak(k,120)*o1de)*ch3br
     &       +  aj(20)*ch3br
     &       +  aj(21)*h1301
     &       +  aj(22)*h1211

#ifndef REDUCED
     &       +  aj(49)*h1202
     &       +2*aj(50)*h2402
     &       +3*(aj(44)+ak(k,131)*oh)*chbr3
     &       +2*(aj(45)+ak(k,132)*oh)*ch2br2
     &       +  (aj(46)+ak(k,133)*oh)*ch2brcl
     &       +  (aj(47)+ak(k,134)*oh)*chbrcl2
     &       +2*(aj(48)+ak(k,135)*oh)*chbr2cl
#endif

     &        )-adj

            lbrx = phobr + phbr        
            if (brx.lt.(lbrx-pbrx)*dt) then
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
      END SUBROUTINE brpandlest
