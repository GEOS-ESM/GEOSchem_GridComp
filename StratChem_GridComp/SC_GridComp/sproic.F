      SUBROUTINE sproic(k,levels,press,m,daytime,
     &  aj,ak,br,bro,brono2,brx,cl,cl2o2,clo,clx,h,h2o,hno3,ho2,n,n2o5,no,no2,
     &                     no3,nox,o1d,o3,o3p,oclo,oh,ox,numphoto,numreacs,o2)

      IMPLICIT NONE
      INTEGER, PARAMETER :: DBL = KIND(0.00D+00)

      INTEGER, INTENT(IN) :: k,levels,numphoto,numreacs
      LOGICAL, INTENT(IN) :: daytime
      REAL, INTENT(IN) :: press(levels),m(levels),o2

      REAL(KIND=DBL), INTENT(IN) :: aj(numphoto),ak(levels,numreacs)
      REAL(KIND=DBL), INTENT(IN) :: brx,clx,nox,o3,hno3,h2o,n2o5

      REAL(KIND=DBL), INTENT(INOUT) :: br,bro,brono2,cl,cl2o2,clo,h,ho2
      REAL(KIND=DBL), INTENT(INOUT) :: n,no,no2,no3,o1d,o3p,oclo,oh,ox

c Initial conditions for inferred species are specified using
c initial conditions for transported species and O3. 
c Modified 7/96 by Randy Kawa and Frank Crum for winter 1995-96 simulation
c... redesigned coding by SDS May 1998

c... number density, reaction and photolysis rates

      REAL(KIND=DBL) :: rprime, a, b, c, rho2oh, term
      REAL(KIND=DBL) :: ratbrbro, ratbrono2bro, h0, oh0
      REAL(KIND=DBL) :: ho20
c   
      h0 = 100000./press(k)
      oh0 = 100000./(press(k)**.666)
      ho20 = 1.e5

c...  Specify O1d and O3p initial conditions. 
c...  daytime o1d ic is determined assuming p.e.  nightime o1d ic is
c...  zero since o3 photolysis is zero.  
      o1d = (aj(2)*o3
     &      +aj(40)*h2o)
     &      /(ak(k,8)*m(k)+ak(k,124)*o3) 	
      o3p = ((aj(2)+aj(3))*o3
     &      +aj(40)*h2o
     &      +2.0*ak(k,124)*o1d*o3)
     &      /(ak(k,1)*m(k)*o2*m(k))
               
c...  Adjust ox
               ox = o3+o1d+o3p 
               
c...  Cl species
               clo = clx
               cl = 0.
               cl2o2 = 0.
               oclo = 0.
               
c...  n, no, no2, and no3 ic's are reconstructed from the nox ic.
c...  for daytime conditions: 
c...  initially assume that [nox]=[no]+[no2].  with the additional
c...  assumption that [ho2] = 0 an expression for [n]/[no2] can be
c...  obtained which can be used with [nox] to estimate [no] and [no2]
c...  daytime [n] and [no3] are calculated by assuming pe. and using
c...  the estimated [no] and [no2]. 
               
c...  for nighttime conditions assume [n] and [no] are zero. assume
c...  [no2] = [nox] and est. [no3] by assuming p.e. and using the est.
c...  [no2]. 
               
c...  finally adjust to conserve [nox]. 
               
               if (daytime) then 
                  rprime = (aj(7)+ak(k,43)*o3p
     &                    +ak(k,125)*h)/
     &                    (aj(16)+aj(39)+ak(k,9)*o3)  
                  no2 = nox/(1.+rprime)  
                  no = rprime*no2  

                 n = (aj(16)+aj(39))*no
     &               /(ak(k,47)*no+ak(k,81)*no2+ak(k,20)*o2*m(k))
                 no3 = (ak(k,10)*no2*o3)
     &               /(aj(5)+aj(17)+ak(k,74)*no)	    
                  
c...  Adjust number densities of n,no,no2, and no3 to ensure that
c...  total nox is unchanged. 
                  rprime = nox / (n+no+
     &                 no2+no3) 
                  n   = n  *rprime
                  no  = no *rprime
                  no2 = no2*rprime
                  no3 = no3*rprime
           
                  rho2oh = ( ak(k,4)*o3
     &                 +ak(k,41)*o3p
     &                 -(ak(k,11)*o3+ak(k,41)*o3p
     &                 +ak(k,125)*no2)
     &                /(ak(k,11)*o3+ak(k,125)*no2
     &                  +ak(k,3)*o2*m(k)**2))
     &                 /(ak(k,5)*o3+ak(k,38)*no
     &                 +ak(k,42)*o3p)

c...  oh and ho2 ic's are defined.
                  a = 2.*ak(k,59)+2.*ak(k,40)*rho2oh
                  b = ak(k,18)*no2*m(k) 
     &              +ak(k,37) * hno3
                  c = -2.*ak(k,39)*h2o*o1d
     &                 -aj(6) * hno3
     &                 -aj(38)*h2o
                  
                  term = sqrt(abs(b*b-4.*a*c)) 
                  oh = max(oh0,-b/(2.*a)+abs(0.5*term/a))

                  if(oh.le.0.0) print *,'SPROIC: oh ',oh,a,c,b
                  
                  ho2 = rho2oh * oh 
                  
c...  h ic is defined
                  h = (ak(k,41)*oh*o3p)/
     &                 ( ak(k,11)*o3+ak(k,125)*no2
     &                  +ak(k,3)*o2*m(k)**2 )
                  
c...  Determine initial conditions for BrO, Br and BrONO2
                  ratbrbro = aj(27)/(ak(k,91)*o3)
                  ratbrono2bro = ak(k,96)*m(k)*no2
     &                 /(aj(28)+ak(k,115))
                  bro = brx/(ratbrbro+1.+ratbrono2bro)
                  if(bro.le.0.0) bro = 1e-31
                  br = ratbrbro*bro

                  brono2 = ratbrono2bro*bro
               else 
                  no = 0. 
                  no2 = nox 
                  
                  n = 0.   
                  
                  no3 = (ak(k,10)*no2*o3
     &                 +ak(k,31)*n2o5*m(k))
     &                 /(ak(k,46)*no2*m(k)) 
                  
c...  Adjust number densities of n,no,no2, and no3 to ensure that
c...  total nox is unchanged. 
                  rprime = nox / (n+no+
     &                 no2+no3) 
                  n   = n  *rprime
                  no  = no *rprime
                  no2 = no2*rprime
                  no3 = no3*rprime
                  oh = oh0  
                  ho2 = ho20  
c...  h ic is defined
                  h = h0  
c...  Determine initial conditions for BrO, Br and BrONO2
                  bro = 0.5*brx
                  if(bro.le.0.0) bro = 1e-31 
                  br =  0.  
               endif

      RETURN
      END SUBROUTINE sproic
