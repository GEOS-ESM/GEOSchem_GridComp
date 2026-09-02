      SUBROUTINE chxpart(k,k1Strat,levels,daytime,m,aj,ak,brx,ch2o,ch3,ch3br,ch3o,ch3o2,
     &             ch3ooh,ch4,cho,clx,ho2,nox,o1d,o3p,oh,rbrnbrx,rbrobrx,rclclo,rcloclx,
     &                                   rno2nox,rnono2,levels_cal,numreacs,numphoto,o2)

      IMPLICIT NONE
      INTEGER, PARAMETER :: DBL = KIND(0.00D+00)

      INTEGER, INTENT(IN) :: k,k1Strat,levels,levels_cal,numreacs,numphoto
      LOGICAL, INTENT(IN) :: daytime
      REAL, INTENT(IN) :: m(levels),o2

      REAL(KIND=DBL), INTENT(IN) :: aj(numphoto),ak(levels,numreacs)
      REAL(KIND=DBL), INTENT(IN) :: brx,ch3br,ch3ooh,ch4,clx,ho2,nox
      REAL(KIND=DBL), INTENT(IN) :: o1d,o3p,oh,rbrnbrx,rbrobrx,rclclo
      REAL(KIND=DBL), INTENT(IN) :: rcloclx,rno2nox,rnono2

      REAL(KIND=DBL), INTENT(INOUT) :: ch2o,ch3,ch3o,ch3o2,cho

      REAL :: o2m
      REAL(KIND=DBL) :: tbr,tcl,tno

c  partition ch4 oxidation products. p and l for ch4,co,ch3ooh are in solver.

      if(k .ge. k1Strat .and. k .le. levels_cal) then
          if(daytime) then
            o2m = o2*m(k)
            tno = max(1.,rnono2*rno2nox*nox)
            tcl = rclclo*rcloclx*clx
            tbr = (1.-rbrobrx-rbrnbrx)*brx
        
            ch3 = ( ak(k,14)*oh+aj(42)
     &        +ak(k,26)*tcl
     &        +ak(k,49)*o1d )*ch4
     &        /(ak(k,50)*m(k)*o2m)
     
            ch3o2 = 
     &        ( ak(k,50)*m(k)*o2m*ch3
     &        +0.7*ak(k,58)*oh*ch3ooh)
     &        /( ak(k,15)*tno+ak(k,22)*ho2 )

            ch3o = ( ak(k,15)*tno*ch3o2
     &        +aj(13)*ch3ooh
     &        +ak(k,119)*o1d*ch4
     &        +ak(k,120)*o1d*ch3br )
     &        /(ak(k,17)*o2m)
    
c... Going to a time dependent ch2o will fix evening and morning diffs
c...  from bench mark, but this has very small effect on hox.  
c... Keep ch2o in steady-state to avoid having to transport it.  

            ch2o = ( ak(k,17)*o2m*ch3o
     &         +ak(k,57)*ch4*o1d
     &         +0.3*ak(k,58)*oh*ch3ooh ) 
     &         /(ak(k,21)*o3p 
     &         +ak(k,51)*oh
     &         +ak(k,63)*tcl
     &         +aj(10)
     &         +aj(11)
     &         +ak(k,69)*tbr)

            cho = ( aj(10)
     &         +ak(k,21)*o3p
     &         +ak(k,51)*oh
     &         +ak(k,63)*tcl
     &         +ak(k,69)*tbr )
     &         *ch2o/(ak(k,52)*o2m)
     
           else
             ch3 = 0.0
             ch3o2 = 0.0
             ch3o = 0.0
             cho = 0.0
c... remove this line for time-dependent ch2o if used
c            ch2o = 0.0
           endif

       endif

       RETURN
       END SUBROUTINE chxpart
