      SUBROUTINE khet3d(k,levels,numreacs,parts,sasulf,sanat,
     &                  saice,m,kel,h2o,hcl,clono2,ak,rc)

! Calculate heterogeneous rate coefficients and put them in rate array ak.  
! The arrays gamma and sa for appropriate surface are determined at calling 
! program.  

!  Inputs: ak (modified as below)-rate array,
!	  kel,
!	  parts - number of particles (cm-3)
!	  sasulf - array of sulfate surface areas,
!	  sanat - array of NAT surface areas,
!	  saice - array of ice surface areas,
!	   hcl - hcl number density
!	   h2o - water vapor concentration (molec cm-3)
!	   j - latitude index
!    clono2 - ClONO2 volume mixing ratio (for Shi/JPL02 gammas), now num. den.
               
!  DATA: gamma-array of sticking coefficients for 3 sfcs x 5 reactions
!	gamma(1,j)->ak(110)=N2O5+H2O
!	gamma(2,j)->ak(112)=ClONO2+H2O
!	gamma(3,j)->ak(111)=ClONO2+HCl
!	gamma(4,j)->ak(113)=N2O5+HCl
!	gamma(5,j)->ak(114)=HOCl+HCl
!	gamma(6,j)->ak(115)=BrONO2+H2O
!	       j: 1=sulfate    <->     itype=0
!		  2=NAT        <->     itype=1
!		  3=ice        <->     itype=2
!
! Mod:930409, 940112
!rk    941111 added BrONO2+H2O
!rk    941121 upgraded ClONO2+H2O to Hanson et al, JGR, 1994 formulation
!rk    951107 upgraded HOCl+HCl to Hanson et al, JGR, 1994 formulation
!             and ClONO2+H2O,HCl to Hanson & Ravi, JPC 1994
!rk 030408 upgraded to JPL 02 forms and rates - as implemented by SDS
!en 20031208 changed from const to named constituents, input with number density

! Rates with H2O are first order equivalents.

! Reaction coefficients for those involving HCl are limited by amount of HCl
! available to adsorb to the surface.  If gas phase HCl is not enough to 
! cover the surface, the rate is reduced by the ratio of [HCl] to the 
! amount of HCl required to cover the particle surface at equilibrium.  The 
! NAT and ice surface coverages are taken from Hanson and Ravi, J. Phys. 
! Chem., 96, 1992.  Sulfate coverage is assumed to be negligble (1.0).  The
! NAT surface coverage is assumed temperature/RH dependent according to the 
! T-dependence of gamma HCl-ClONO2 given by Hanson and Ravi, JGR, 1993. 

! T/RH-dependent gammas for HCl-ClONO2 and ClONO2-H2O on NAT are
! also taken from the latter reference.

! Gamma for BrONO2+H2O taken as 0.4 from Hansen and Ravi, GRL, 1995
!  and assumed constant with wt % and particle type for now. 
!  Later updated to JPL02.


      IMPLICIT NONE

      INTEGER, INTENT(IN) :: k,levels,numreacs
      INTEGER, INTENT(OUT) :: rc
      REAL, INTENT(IN) :: parts,sasulf(levels)
      REAL, INTENT(IN) :: kel(levels)
      REAL, INTENT(IN) :: m(levels)

      INTEGER, PARAMETER :: DBL = KIND(0.00D+00)
      REAL(KIND=DBL), INTENT(IN) :: hcl,h2o,clono2,sanat,saice
      REAL(KIND=DBL), INTENT(INOUT) :: ak(levels,numreacs)

      INTEGER :: i,ik,il,ir,itype

!     (cm2/cm3)

      REAL(KIND=DBL) :: a,ai,ak0,ak1,ak2,akcheck,akmax,av,b,bi
      REAL(KIND=DBL) :: cfreq,cpr,gamma1,gamma2,gamma3,gamma4,gamma5,gamma6
      REAL(KIND=DBL) :: g2,g3,g5,gr,h2omr,hcle,pr,r,rad,sa1,sparts,t1,w
      REAL(KIND=DBL) :: wtbrn,wtcln,wthocl,wtn5,z,hclmr,clono2mr

      LOGICAL :: gotHCle
      
!     sulfate,NAT,ice sfc HCl uptake (molec/cm2) -H&R,1992
!     Sulfate uptake changed from 1.e10 to 1.(951107) to 1e6 (Randy 6/96)
      REAL, SAVE :: hclsfc(3) = (/1.e6, 2.e14, 5.e14/)
      
!     coefficients for JPL00 N2O5+H2O wt pct dependence
      REAL, SAVE :: c0(0:3) = (/ -25.5265,-0.133188,0.00930846,-9.0194e-5 /)
      REAL, SAVE :: c1(0:3) = (/ 9283.76,115.345,-5.19258,0.0483464 /)
      REAL, SAVE :: c2(0:3) = (/ -851801.,-22191.2,766.916,-6.85427 /)
      REAL, SAVE :: gamma(6,3) ! Sulfate, NAT, ice

      REAL, PARAMETER :: PI = 3.14159265

      gamma = RESHAPE((/ 0.1   , 0.006, 0.0006, 0.0  , 0.0, 0.4,
     &                   0.0004, 0.004, 0.2   , 0.003, 0.1, 0.3,
     &                   0.02  , 0.3  , 0.3   , 0.03 , 0.2, 0.3 /),(/6,3/))

      rc = 0

!     Set ak(*,110) - ak(*,115) to zero.  They are set in this sub
      
      do ir=110,115
       ak(k,ir) = 0.0
      enddo

      wtn5 = 0.00925926         ! = 1./108.
      wtcln = 0.0102617         ! = 1./97.45
      wthocl =  0.0190658       ! = 1./52.45
      wtbrn = 0.00704722        ! = 1./141.9
!     Marti&Mauersberger,GRL,1993
      ai = -2663.5
      bi = 12.537
!     coeff's for T-dep of gamma3 on NAT from H&R, 1993 preprint
      a = 0.7022
      b = -0.518
!     (J/K-mole),(molec/mole)
      r = 8.3143
      av = 6.022e23
!     rate limit
      akmax = 5.e-4

      hcle=0.00
      gotHCle=.false.
      t1 = kel(k)
               
!     mixing ratio to pressure(Pa)
      cpr = m(k)*r*t1*1.e6/av

      if(saice .GT. 0.00) then
         itype = 2
         sa1 = saice
         gamma1 = gamma(1,3) 
         gamma2 = gamma(2,3) 
         gamma3 = gamma(3,3)
         gamma5 = gamma(5,3)
         gamma6 = gamma(6,3)
!     max [HCl] on surface    
         hcle = sa1*hclsfc(3)
         gotHCle=.true.
      endif

!     NAT T/RH-dep
      if(sanat .GT. 1.e-9 .and. .not. gotHCle) then
         itype = 1
         sa1 = sanat
                  
!     calculate relative humidity dependence of HCl+ClONO2 on NAT
!     set T-indpnt values
         gamma2 = gamma(2,2)
!     JPL 02
         gamma3 = gamma(3,2)
!     set T-indpnt values
         gamma1 = gamma(1,2)
         gamma5 = gamma(5,2)
         gamma6 = gamma(6,2)
!     T-dep of sfc coverage
         hcle = sa1*hclsfc(2)
!     hcle = sa1*hclsfc(2)*gamma3/gamma(3,2)
         gotHCle=.true.
      endif
               
      if(sasulf(k) .gt. 0.00 .and. .not. gotHCle) then
!     sulfate
        sa1 = sasulf(k)
        itype = 0
        sparts = parts
!     rk  Begin Hanson et al 1994 formulation for reactions in sulfate
!     ln pressure in mb
       z = LOG(cpr*h2o/(m(k)*100.))
!     wt% H2SO4(Hanson& 94)
        w = ((-14.458+0.62456*z)*t1+3565.)/
     &       (44.777+1.3204*z-0.19988*t1) 
!     (rk) Found w < 0 when H2O is high, T low, e.g., 10, 185; 100,200
        w = min(100.,max(w,0.01))
!     JPL00 form for N2O5+H2O(l)
        ak0 = c0(0)
        ak1 = c1(0)
        ak2 = c2(0)
        DO i=1,3
          ak0 = ak0+c0(i)*w**i
          ak1 = ak1+c1(i)*w**i
          ak2 = ak2+c2(i)*w**i
        ENDDO
        gamma1 = EXP(ak0+ak1/t1+ak2/(t1**2.))
!     JPL00 for BrONO2+H2O(l)
        gr = EXP(29.24-.396*w)+0.114
        gamma6 = 1./(1./.805+1./gr)
               
!     radius (um)
        rad = sa1/sparts/4./pi
        rad = sqrt(rad)
!     mbar
        pr = cpr/100.
        h2omr = h2o/m(k)
        hclmr = hcl/m(k)
        clono2mr = clono2/m(k)
!     JPL 02

        CALL gamma_shi(t1,pr,h2omr,hclmr,
     &                 clono2mr,rad,g2,g3,g5)

!     gamma ClONO2+HCl
        gamma3 = g3
!     gamma ClONO2+H2O
        gamma2 = g2
               
!     HOCl+HCl
        IF(hcl .GT. 0.00) THEN 
          gamma5 = min(g5,1.)
        ELSE
          gamma5 = 0.
        ENDIF
        hcle = sa1*hclsfc(1)
        gotHCle=.true.

      endif  !finish sasulf conditional

!     Do the work when an HCLe has been defined

      if(gotHCle) then
!     max hcl for reaction
        if (hcl .gt. hcle) hcle = hcl

        IF(hcle <= 0.00) THEN
	 rc = 1
	 PRINT *,"SC::khet3d found hcle <= 0.00"
	 RETURN
	END IF

!     N2O5 (s-1)
        cfreq = SQRT(3.*8314.*t1*wtn5)*100.*sa1*.25
        ak(k,110) = gamma1*cfreq			 
        ak(k,113) = gamma(4,itype+1)*cfreq/hcle
               
!     ClONO2
        cfreq = SQRT(3.*8314.*t1*wtcln)*100.*sa1*.25
!     (s-1)
        ak(k,112) = gamma2*cfreq

!     (cm3/molec-s)      
        ak(k,111) = gamma3*cfreq/hcle
               
!     HOCl
        cfreq = SQRT(3.*8314.*t1*wthocl)*100.*sa1*.25
        ak(k,114) = gamma5*cfreq/hcle
               
!     BrONO2
        cfreq = SQRT(3.*8314.*t1*wtbrn)*100.*sa1*.25
!     (s-1)
        ak(k,115) = gamma6*cfreq
               
!     Constrain reaction rates to max. of akmax
        ak(k,110) = min(ak(k,110),akmax)
        ak(k,112) = min(ak(k,112),akmax)
        ak(k,115) = min(ak(k,115),akmax)

        akcheck = akmax / hcle
        ak(k,111) = min(ak(k,111),akcheck)
        ak(k,113) = min(ak(k,113),akcheck)
        ak(k,114) = min(ak(k,114),akcheck)
        
        ak(k,110) = max(ak(k,110),0.0)
        ak(k,111) = max(ak(k,111),0.0)
        ak(k,112) = max(ak(k,112),0.0)
        ak(k,113) = max(ak(k,113),0.0)
        ak(k,114) = max(ak(k,114),0.0)
        ak(k,115) = max(ak(k,115),0.0)
        
      endif

      RETURN
      END SUBROUTINE khet3d
