      SUBROUTINE jcalc4(k,szan,o3column,press,kel,aj,gcSC,rdistsq)
! ---------------------------------------------------------------------------------
! NAME: jcalc4
! PURPOSE:
!   Calculate photolysis rates
! INPUT:
!   k         Current layer number
!   levels    Number of layers
!   szan      Solar zenith angle (radians)
!   o3column  Overhead O3 values
!   press     Mid-layer pressure (hPa)
!   kel       Mid-layer temperature (K)
!   rdistsq   Scale solar by eccentricity factor
! OUTPUT:
!   aj        Array of photolysis rates
! RESTRICTIONS:
!   Currently set up for 23-J set (see var gcSC%nxdo)
! REQUIRED ROUTINES:
!   interp_s
! MODIFICATION HISTORY: 
!   26 Aug 1993 Kawa    Created
!   23 Nov 1993 Kawa    Remade xtab to do multiplication by solar flux beforehand 
!                        and removed inputs.
!   25 Feb 1994         Add 3 additional Js, incl N2O
!   18 Sep 1995         Add 2 additional Js, up to 22, and do CH2O special
!   13 May 1996 Crum    Removed fossils, move toward Fortran 90
!   10 Jul 1996         Modified to handle J(O2) separately and use 28 levels
!    1 Apr 2009 Nielsen GEOS-5 form with standardized SC_GridComp interface.
!    1 Jun 2009 Nielsen Updated to JPL 2006
!   12 Dec 2010 Nielsen Updated to JPL 2010 following Luke Oman's testing.
!   11 May 2012 Nielsen Accomodation for GEOS-5 FV cubed release
!    2 Sep 2014 Oman    Added factor to account for eccentricity
!    3 Jun 2015 Liang   Updated to the new 50-slot table with addition of halons,
!                       HCFCs, and 5 VSLSs
!                       numphoto is now updated to 52
!
! WARNING: Photolysis reaction rate numbers 38-42 are calculated in MESO_PHOT.
! ---------------------------------------------------------------------------------

      USE SC_GridCompMod
      IMPLICIT NONE
      INTEGER, PARAMETER :: DBL = KIND(0.00D+00)

      TYPE(SC_GridComp), INTENT(INOUT) :: gcSC   ! Grid Component

      INTEGER, INTENT(IN) :: k
      REAL, INTENT(IN) :: szan, o3column, press, kel,rdistsq
      REAL(KIND=DBL), INTENT(OUT) :: aj(gcSC%numphoto)

      INTEGER :: ilam,indt,ix

      REAL :: alpha300, alphat, jo2, rjs(gcSC%nxdo), q1, q2, r1mq1
      REAL :: s(gcSC%nlam), sx(2,gcSC%nlam), tfac, wvl

! Start with a clean slate
! ------------------------
      aj(1:gcSC%numphoto) = 0.

! Interpolate radiative flux function values to model conditions
! --------------------------------------------------------------
      CALL interp_s(k,szan,o3column,s,jo2,gcSC)
      indt = kel-148.5
      indt = MAX(1,indt)
      indt = MIN(indt,200)

! Preliminaries for CH2O quantum yield dependence on m, T, wavelength
! -------------------------------------------------------------------
      tfac = (kel-80.0)/80.0

      DO ilam=1,gcSC%nlam

       ZeroS: IF(s(ilam) == 0.) THEN
        sx(1,ilam) = 0.00
	sx(2,ilam) = 0.00
       ELSE 

        wvl = gcSC%rlam(ilam)*0.10

        IF(wvl < 250.00) THEN
         q1 = 0.24
        ELSE IF(wvl >= 339.00) THEN
         q1 = 0.00
        ELSE
         q1 = gcSC%CH2O_aq(1) + gcSC%CH2O_aq(2)*wvl +
     &                          gcSC%CH2O_aq(3)*wvl*wvl +
     &                          gcSC%CH2O_aq(4)*wvl*wvl*wvl + 
     &                          gcSC%CH2O_aq(5)*wvl*wvl*wvl*wvl
        END IF

        r1mq1 = 1./(1.-q1)

        IF(wvl < 330.00) THEN
!        q2 = gcSC%xtab(ilam,22,indt)
         q2 = gcSC%xtab(ilam,22,indt)*rdistsq
        ELSE IF(wvl > 360.00) THEN
	 q2 = 0.00
	ELSE
!        alpha300 = 1.00E-03*(1./gcSC%xtab(ilam,22,1)-r1mq1)
         alpha300 = 1.00E-03*(1./(gcSC%xtab(ilam,22,1)*rdistsq)-r1mq1)
         alphat = alpha300*(1.+0.05*(wvl-329.)*((300.-kel)/80.))
         q2 = 1.00/(r1mq1+alphat*press)
        END IF

        IF(wvl .LT. 250.00) q2=0.5

!       sx(2,ilam) = s(ilam)*gcSC%xtab(ilam,21,indt)*q2
!       sx(1,ilam) = s(ilam)*gcSC%xtab(ilam,21,indt)*q1
        sx(2,ilam) = s(ilam)*gcSC%xtab(ilam,21,indt)*q2*rdistsq
        sx(1,ilam) = s(ilam)*gcSC%xtab(ilam,21,indt)*q1*rdistsq

       END IF ZeroS

      END DO

! J(BrONO2) through J(OCLO)
! -------------------------
      DO ix=1,14
         rjs(ix) = 0.

         DO ilam=1,gcSC%nlam
!           rjs(ix) = rjs(ix)+s(ilam)*gcSC%xtab(ilam,ix,indt)
            rjs(ix) = rjs(ix)+s(ilam)*gcSC%xtab(ilam,ix,indt)*rdistsq
         END DO

      END DO

! J(O2)
! -----
      rjs(15) = jo2

! J(O3_O1D) through J(N2O)
! ------------------------
      DO ix=16,20
       rjs(ix) = 0.

       DO ilam=1,gcSC%nlam
!       rjs(ix) = rjs(ix)+s(ilam)*gcSC%xtab(ilam,ix,indt)
        rjs(ix) = rjs(ix)+s(ilam)*gcSC%xtab(ilam,ix,indt)*rdistsq
       END DO

      END DO

! J(CH2O)
! -------
      rjs(21) = 0.
      rjs(22) = 0.
      DO ilam=1,gcSC%nlam
       rjs(21) = rjs(21)+sx(1,ilam)
       rjs(22) = rjs(22)+sx(2,ilam)
      END DO


! J(CO2 -> CO + O) through xH1211
! -------------------------------
      DO ix=23,gcSC%nxdo
       rjs(ix) = 0.

       DO ilam=1,gcSC%nlam
!       rjs(ix) = rjs(ix)+s(ilam)*gcSC%xtab(ilam,ix,indt)
        rjs(ix) = rjs(ix)+s(ilam)*gcSC%xtab(ilam,ix,indt)*rdistsq
       END DO

      END DO
               
! ---------------------------------------------------------------
! Order photolysis rates to match order in full chemistry model.  
! Sort rjs into CTM photolysis rate array, aj.  Order of rjs:
!
!  1-J(BrONO2)
!  2-J(BrO)
!  3-J(Cl2O2)
!  4-J(ClONO2)
!  5-J(H2O2)
!  6-J(HCl)
!  7-J(HNO3)
!  8-J(HO2NO2)
!  9-J(HOCl)
! 10-J(N2O5)
! 11-J(NO2)
! 12-J(NO3_NO)
! 13-J(NO3_NO2)
! 14-J(OClO)
! 15-J(O2)
! 16-J(O3_O1D)
! 17-J(O3_3P)
! 18-J(HOBr)
! 19-J(CH3OOH)
! 20-J(N2O)
! 21-J(CH2O_HCO)
! 22-J(CH2O_CO)
! 23-J(CO2 -> CO + O)
! 24-xCFC-11
! 25-xCFC-12
! 26-xCCl4
! 27-xCH3CCl3
! 28-xHCFC-22
! 29-xCFC-113
! 30-xCH3Cl
! 31-xCH3Br
! 32-xH1301
! 33-xH1211 
! 34-xH1202
! 35-xH2402
! 36-xCHBr3
! 37-xCH2Br2
! 38-xCH2ClBr
! 39-xCHClBr2
! 40-xCHCl2Br
! 41-xHCFC-141b
! 42-xHCFC-142b
! 43-xCFC-114 
! 44-xCFC-115
! 45-xOCS
! 46-
! 47-
! 48-
! 49-
! 50-
! ---------------------------------------------------------------
! Solar cycle goes here when ready  
!     aj( 1) = rjs(15)*gcSC%s_cycle(3,gcSC%iscyr)
! ----------------------------------------------------------------
      aj( 1) = rjs(15)
      aj( 2) = rjs(16)
      aj( 3) = rjs(17)
! H2O
! ---
      aj( 4) = 0.
      aj( 5) = rjs(13)
      aj( 6) = rjs(7)
      aj( 7) = rjs(11)
      aj( 8) = rjs(5)
      aj( 9) = rjs(10)
      aj(10) = rjs(21)
      aj(11) = rjs(22)
      aj(12) = rjs(23)
      aj(13) = rjs(19)
      aj(14) = rjs(20)
      aj(15) = rjs(4)
      aj(16) = 0.
      aj(17) = rjs(12)
      aj(18) = rjs(6)
      aj(19) = 0.

! CH3Br(20) H1301(21) H12_24(22)
! ------------------------------
      aj(20) = rjs(31)
      aj(21) = rjs(32)
      aj(22) = rjs(33)
      aj(23) = rjs(9)
      aj(24) = rjs(8)
      aj(25) = rjs(18)
      aj(26) = 0.
      aj(27) = rjs(2)
      aj(28) = rjs(1)

! F11(29) F12(30) CCl4(31) CHCCl3(32) HCFC(33) F113(34) CH3Cl(35)
! ---------------------------------------------------------------
      aj(29) = rjs(24)
      aj(30) = rjs(25)
      aj(31) = rjs(26)
      aj(32) = rjs(27)
      aj(33) = rjs(28)
      aj(34) = rjs(29)
      aj(35) = rjs(30)
      aj(36) = rjs(3)
      aj(37) = rjs(14)

! ------------------------------------------
! WARNING: Photolysis reaction rate
! numbers 38-42 are calculated in MESO_PHOT.
! ------------------------------------------
! Add aj(43) which is J(Cl2O2) for partitioning but not Ox loss 
! which is aj(36). In lookup table J(Cl2O2) is J*qy where qy is 0.8 
! so multiply by 1.25 to equal J and used in part.F and partest.F

      aj(43) = rjs(3)*1.25

! QingLiang -- 06/03/2015
! CHBr3(44) CH2Br2(45) CH2BrCl(46) CHBrCl2(47) CHBr2Cl(48)
      aj(44) = rjs(36)
      aj(45) = rjs(37)
      aj(46) = rjs(38)
      aj(47) = rjs(39)
      aj(48) = rjs(40)

! QingLiang -- 06/03/2015
! Add two new halons: H-1202 (49) H2402 (50) 
! and two new HCFCs: HCFC-141b (51) HCFC-142b (52) 
      aj(49) = rjs(34)
      aj(50) = rjs(35)
      aj(51) = rjs(41)
      aj(52) = rjs(42)

! QingLiang -- 02/05/2016
! Add CFC-114 and CFC-115
! Add OCS for GOCART module
      aj(53) = rjs(43)
      aj(54) = rjs(44)
      aj(55) = rjs(45)
!      aj(53) = rjs(34)
!      aj(54) = rjs(34)
!      aj(55) = rjs(34)

      RETURN
      END SUBROUTINE jcalc4
