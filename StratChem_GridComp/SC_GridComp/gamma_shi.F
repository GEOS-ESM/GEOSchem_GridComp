      SUBROUTINE gamma_shi(t,pr,h2o,hcl,cln,rad,g2,g3,g5)
!
! Stratospheric Reaction Coefficients of ClONO2 and HOCl
! Shi, Q., P. Davidovits, Boston College, Chestnut Hill, MA 02167
! D. R. Worsnop, T. J. Jayne, C. E. Kolb, Aerodyne Research, Inc., Billerica,
! MA 01821
! worsnop@aerodyne.com//shiq@aerodyne.com
! INSTRUCTIONS FOR USE:
!
! *INPUT*	
!	t -- The user must provide Temperature (K, scalar or vector)
! --scalars or vector of same length as t --
!	pr -- Pressure (mb)
!	h2o -- water vapor mixing ratio
!	hcl -- HCl mixing ratio
!	cln -- ClONO2 mixing ratio
!	rad -- Aerosol radius (cm), typical value is 1e-5
!	 
! *OUTPUT*
! Heterogeneous reaction probabilities for:
! ClONO2 reacting with HCl and H2O
! HOCl reacting with HCl.
! These are symbolized in the following code as:
! gClONO2_HCl, gClONO2_H2O, gHOCl
!	 -> output g3, g2, g5 respectively
!
! Convert to Water vapor partial press (mbar)
! HCl and ClONO2 partial pressure (atm).
! These are symbolized in the following code, respectively, as:
! T, r, PH2O, PHCl, PClONO2

      IMPLICIT NONE
      INTEGER, PARAMETER :: DBL = KIND(0.00D+00)

      REAL(KIND=DBL), INTENT(IN) :: t,pr,h2o,hcl,cln,rad
      REAL(KIND=DBL), INTENT(OUT) :: g2,g3,g5

      REAL(KIND=DBL) :: ph2o,phcl,pclono2,p0H2O,aw,a1,b1,c1,d1
      REAL(KIND=DBL) :: a2,b2,c2,d2,y1,y2,am,wt,z1,z2,z3,rho
      REAL(KIND=DBL) :: aMSO4,X,A,To,h,aH,c_ClONO2,SClONO2,HClONO2
      REAL(KIND=DBL) :: DClONO2,akH2O,akH,akhydr,GbH2O,HHCl,aMHCl
      REAL(KIND=DBL) :: akHCl,alClONO2,fClONO2,GClONO2rxn,GbHCl
      REAL(KIND=DBL) :: Gs,FHCl,Gsp,GbHClp,Gb,gClONO2,gClONO2_HCl
      REAL(KIND=DBL) :: gClONO2_H2O,gHOCl,c_HOCl,SHOCl,HHOCl
      REAL(KIND=DBL) :: DHOCl,akHOCl_HCl,GHOClrxn,alHOCl,fHOCl,eval

! Bound water

      aw=max(5e-11,h2o)
      aw=min(aw,7e-03)

      ph2o = pr*aw
      phcl = pr*hcl/1013.
      pclono2 = pr*cln/1013.

! //Table 1: H2SO4 wt% from T and pH2O

      p0H2O = EXP(18.452406985-3505.1578807/T-330918.55082/T**2.
     &  +12725068.262/T**3.)
      aw = pH2O/ p0H2O
      aw = min(aw,1.)			!limit fit to 1

      if (aw < 0.05) THEN
         a1 = 12.37208932
         b1 = -0.16125516114
         c1 = -30.490657554
         d1 = -2.1133114241
         a2 = 13.455394705
         b2 = -0.1921312255
         c2 = -34.285174607
         d2 = -1.7620073078
      else if (aw >= 0.85) THEN 
         a1 = -180.06541028
         b1 = -0.38601102592
         c1 = -93.317846778
         d1 = 273.88132245
         a2 = -176.95814097
         b2 = -0.3625704815410
         c2 = -90.469744201
         d2 = 267.45509988
!     if (aw >= 0.05) THEN
      else
         a1 = 11.820654354
         b1 = -0.20786404244
         c1 = -4.807306373
         d1 = -5.1727540348
         a2 = 12.891938068
         b2 = -0.23233847708
         c2 = -6.4261237757
         d2 = -4.9005471319
      endif

      y1 = a1*aw**b1+c1*aw+d1
      y2 = a2*aw**b2+c2*aw+d2
      am = y1+(T-190.)*(y2-y1)/70.
      wt = 9800.*am/(98.*am+1000.)
      wt = min(100.,max(wt,1.e-7))

! //Table 2: Parameters for H2SO4 Solution
      z1 = 0.12364-5.6e-7*T**2
      z2 = -0.02954+1.814e-7*T**2
      z3 = 2.343e-3-1.487e-6*T-1.324e-8*T**2
      rho = 1+z1*am+z2*am**1.5+z3*am**2
      aMSO4 = rho*wt/9.8
      X = wt/(wt+(100.-wt)*98./18.)
      A = 169.5+5.18*wt-0.0825*wt**2.+3.27e-3*wt**3.
      To = 144.11+0.166*wt-0.015*wt**2.+2.18e-4*wt**3.
      eval = max(T-To,1.0)
      h = A*T**(-1.43)*exp(448.00/eval)
      aH = exp(60.51-0.095*wt+0.0077*wt**2.-1.61e-5*wt**3.-
     &  (1.76+2.52e-4*wt**2)*sqrt(T)+(-805.89+253.05*wt**0.076)/sqrt(T))

! //Table 3. ClONO2 + H2O and ClONO2+HCl
      c_ClONO2 = 1474*sqrt(T)
      SClONO2 = 0.306+24.0/T
      HClONO2 = 1.6e-6*exp(4710./T)*exp(-SClONO2*aMSO4)
      DClONO2 = 5e-8*T/h
      akH2O = 1.95e10*exp(-2800./T)
      akH = 1.22e12*exp(-6200./T)
      akhydr = akH2O*aw + akH*aH*aw
      GbH2O = 4.*HClONO2*0.082*T*sqrt(DClONO2*akhydr)/c_ClONO2
      HHCl = (0.094-0.61*X+1.2*X**2.)*exp(-8.68+(8515-10718*X**0.7)/T)
      aMHCl = HHCl *PHCl
      akHCl = 7.9e11*aH*DClONO2*amhcl
      alClONO2 = sqrt(DClONO2/(akhydr+akHCl))
      fClONO2 = 1/tanh(rad/alClONO2)- alClONO2/rad
      GClONO2rxn = fClONO2*GbH2O *sqrt(1+akHCl/akhydr)
      GbHCl = GClONO2rxn* akHCl/( akHCl+ akhydr)
      Gs = 66.12*exp(-1374/T)*HClONO2*amhcl
      if(phcl.eq.0.0) then 
         FHCl = 1
       else
         FHCl = 1/(1+0.612*(Gs+GbHCl)* pClONO2/ pHCl)
       endif

      Gsp = FHCl*Gs
      GbHClp = FHCl*GbHCl
      Gb = GbHClp + GClONO2rxn* akhydr/( akHCl+ akhydr)
      if((Gsp + Gb).eq.0.0) then
        gClONO2 = 0.0
        gClONO2_HCl = 0.0
       else
        gClONO2 = 1./(1.+1./(Gsp + Gb))
        gClONO2_HCl = gClONO2 *(Gsp + GbHClp)/(Gsp + Gb)
       endif
      if(gClONO2_HCl.gt.gClONO2) then
        gClONO2_H2O = 0.0
       else
        gClONO2_H2O = gClONO2 - gClONO2_HCl
       endif

! //Table 4. HOCl + HCl
      if(phcl.eq.0.0) then 
         gHOCl = 0.0
       else
         c_HOCl = 2009*sqrt(T)
         SHOCl = 0.0776+59.18/T
         HHOCl = 5862.4/t
         HHOCl = 1.91e-6*exp(HHOCl)
         HHOCl = HHOCl*exp(-SHOCl*aMSO4)
         DHOCl = 6.4e-8*t/h
         akHOCl_HCl = 1.25e9*aH*DHOCl*amhcl
         GHOClrxn = 4*HHOCl*0.082*T
         GHOClrxn = GHOClrxn*sqrt(DHOCl*akHOCl_HCl)
         GHOClrxn = GHOClrxn/c_HOCl
         alHOCl = sqrt(DHOCl/akHOCl_HCl)
         fHOCl = 1/tanh(rad/alHOCl)- alHOCl/rad
         if((fHOCl*GHOClrxn*FHCl).le.0.0) then
            gHOCl = 0.0
         else
            gHOCl = 1/(1+1/ ( fHOCl*GHOClrxn*FHCl))
         endif
       endif
 
      g2 = gClONO2_H2O
      g3 = gClONO2_HCl
      g5 = gHOCl

      RETURN
      END SUBROUTINE gamma_shi
