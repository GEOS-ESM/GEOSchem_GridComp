!=======================================================================
!
! $Id$
!
! ROUTINE
!   kcalc - GMIMOD 3D model (setkin_kcalc.F90)
!   1 AUG 02 - PSC
!
! DESCRIPTION
!   Calculates and returns rate constants for thermal
!   reactions for the temperatures and pressures supplied
!
! ARGUMENTS
!  INPUT
!   pressure    : mb - profile of pressures
!   temperature : K - profile of temperatures
!   adcol       : molecules cm-3 - profile of total number density
!   specarr     : molecules cm-3 - profiles of species concentrations
!  OUTPUT
!   rcarr       : cm3 molecule-1 s-1 - rate constant values in units as
!                 appropriate
!
!  Chemistry input file:    08/2018
!  Reaction dictionary:     GMI_Combo_rxns_119species_SO2_JPL15_OCS.db
!  Setkin files generated:  Tue Sep 18 18:25:41 2018
!
!=======================================================================
      subroutine kcalc( npres0,sadcol,sadcol2,pressure,ptrop,cPBLcol, &
     &  temperature,lwc,adcol,specarr,rcarr,radA,FRH)

      implicit none

#     include "gmi_phys_constants.h"
#     include "setkin_par.h"
#     include "setkin_mw.h"

!.... Argument declarations

      INTEGER, INTENT (IN)  :: npres0
      INTEGER, INTENT (IN)  :: cPBLcol(npres0)

      REAL*8,  INTENT (IN)  :: ptrop

      REAL*8,  INTENT (IN)  :: adcol       (       npres0)
      REAL*8,  INTENT (IN)  :: lwc         (       npres0)
      REAL*8,  INTENT (IN)  :: pressure    (       npres0)
      REAL*8,  INTENT (IN)  :: temperature (       npres0)
      REAL*8,  INTENT (IN)  :: sadcol      (NSAD  ,npres0)
      REAL*8,  INTENT (IN)  :: sadcol2     (NSADaer+NSADdust,npres0)
      REAL*8,  INTENT (IN)  :: specarr     (NMF   ,npres0)
      REAL*8,  INTENT (IN)  :: FRH         (       npres0)
      REAL*8,  INTENT (IN)  :: radA        (NSADaer+NSADdust,npres0)

      REAL*8,  INTENT (OUT) :: rcarr       (NUM_K ,npres0)

!.... Local variable declarations

      INTEGER               :: naltmax

      real*8 &
     &  nitrogen (npres0) &
     & ,oxygen   (npres0) &
     & ,water    (npres0)

! The sad_* variables are not used in the Tropospheric Mechanism.
      real*8 &
     &  sad_ice  (npres0) &
     & ,sad_lbs  (npres0) &
     & ,sad_nat  (npres0) &
     & ,sad_soot (npres0) &
     & ,sad_sts  (npres0)

      real*8 mw(NSP)

      mw(:) = mw_data(:)

      naltmax     = npres0
!===================================================
! Variables needed for gas/heterogeneous chemistry.
! adcol = air density (molec/cm3)
! FRH = relative humidity fraction (0-1)
! radA = effective radius of aerosol (cm)
! sadcol = surface area of aerosol/volume of air (cm2/cm3)
!===================================================

!....          Define molecular nitrogen and oxygen number densities
!
      nitrogen(:) = adcol(:) * MXRN2
      oxygen(:)   = adcol(:) * MXRO2
      water(:)    = specarr(10 ,:)

      sad_lbs(:)  = sadcol(1, :)
      sad_sts(:)  = sadcol(2, :)
      sad_nat(:)  = sadcol(3, :)
      sad_ice(:)  = sadcol(4, :)
      sad_soot(:) = sadcol(5, :)
!
!....          Start thermal rate constants
!
!....           O + O2 = O3
!
      rcarr(1,:) = skterlp(  6.000D-34 ,2.40D+00 ,0.0D0 &
     &                     ,temperature ,adcol)
!
!....           O + O3 = 2 O2
!
      rcarr(2,:) = skarr(  8.000D-12 ,2060.0D+00 ,temperature)
!
!....           N2 + O1D = N2 + O
!
      rcarr(3,:) = skarr(  2.150D-11 ,-110.0D+00 ,temperature)
!
!....           O1D + O2 = O + O2
!
      rcarr(4,:) = skarr(  3.300D-11 ,-55.0D+00 ,temperature)
!
!....           O1D + O3 = 2 O2
!
      rcarr(5,:) = skarr(  1.200D-10 ,0.0D+00 ,temperature)
!
!....           O1D + O3 = 2 O + O2
!
      rcarr(6,:) = skarr(  1.200D-10 ,0.0D+00 ,temperature)
!
!....           H2O + O1D = 2 OH
!
      rcarr(7,:) = skarr(  1.630D-10 ,-60.0D+00 ,temperature)
!
!....           H2 + O1D = H + OH
!
      rcarr(8,:) = skarr(  1.200D-10 ,0.0D+00 ,temperature)
!
!....           N2O + O1D = N2 + O2
!
      rcarr(9,:) = skarr(  4.640D-11 ,-20.0D+00 ,temperature)
!
!....           N2O + O1D = 2 NO
!
      rcarr(10,:) = skarr(  7.260D-11 ,-20.0D+00 ,temperature)
!
!....           CH4 + O1D = MO2 + OH
!
      rcarr(11,:) = skarr(  1.310D-10 ,0.0D+00 ,temperature)
!
!....           CH4 + O1D = CH2O + H + HO2
!
      rcarr(12,:) = skarr(  3.500D-11 ,0.0D+00 ,temperature)
!
!....           CH4 + O1D = CH2O + H2
!
      rcarr(13,:) = skarr(  8.750D-12 ,0.0D+00 ,temperature)
!
!....           CF2Cl2 + O1D = 2 Cl
!
      rcarr(14,:) = skarr(  1.400D-10 ,-25.0D+00 ,temperature)
!
!....           CFC113 + O1D = 3 Cl
!
      rcarr(15,:) = skarr(  2.320D-10 ,0.0D+00 ,temperature)
!
!....           CFC114 + O1D = 2 Cl
!
      rcarr(16,:) = skarr(  1.300D-10 ,-25.0D+00 ,temperature)
!
!....           CFC115 + O1D = Cl
!
      rcarr(17,:) = skarr(  5.400D-11 ,-30.0D+00 ,temperature)
!
!....           HCFC22 + O1D = Cl
!
      rcarr(18,:) = skarr(  1.020D-10 ,0.0D+00 ,temperature)
!
!....           HCFC141b + O1D = 2 Cl
!
      rcarr(19,:) = skarr(  2.600D-10 ,0.0D+00 ,temperature)
!
!....           HCFC142b + O1D = Cl
!
      rcarr(20,:) = skarr(  2.000D-10 ,0.0D+00 ,temperature)
!
!....           H + O2 = HO2
!
      rcarr(21,:) = sktroe(  4.400D-32 ,1.30D0 & 
     &                     , 7.500D-11 ,-0.20D0 ,0.0D0 &
     &                     ,temperature ,adcol)
!
!....           H + O3 = O2 + OH
!
      rcarr(22,:) = skarr(  1.400D-10 ,470.0D+00 ,temperature)
!
!....           O + OH = H + O2
!
      rcarr(23,:) = skarr(  1.800D-11 ,-180.0D+00 ,temperature)
!
!....           HO2 + O = O2 + OH
!
      rcarr(24,:) = skarr(  3.000D-11 ,-200.0D+00 ,temperature)
!
!....           H + HO2 = 2 OH
!
      rcarr(25,:) = skarr(  7.200D-11 ,0.0D+00 ,temperature)
!
!....           NO + O3 = NO2 + O2
!
      rcarr(26,:) = skarr(  3.000D-12 ,1500.0D+00 ,temperature)
!
!....           O3 + OH = HO2 + O2
!
      rcarr(27,:) = skarr(  1.700D-12 ,940.0D+00 ,temperature)
!
!....           HO2 + O3 = 2 O2 + OH
!
      rcarr(28,:) = skarr(  1.000D-14 ,490.0D+00 ,temperature)
!
!....           NO2 + O3 = NO3 + O2
!
      rcarr(29,:) = skarr(  1.200D-13 ,2450.0D+00 ,temperature)
!
!....           OH + OH = H2O + O
!
      rcarr(30,:) = skarr(  1.800D-12 ,0.0D+00 ,temperature)
!
!....           OH + OH = H2O2
!
      rcarr(31,:) = sktroe(  6.900D-31 ,1.00D0 & 
     &                     , 2.600D-11 ,0.00D0 ,0.0D0 &
     &                     ,temperature ,adcol)
!
!....           HO2 + OH = H2O + O2
!
      rcarr(32,:) = skarr(  4.800D-11 ,-250.0D+00 ,temperature)
!
!....           H2O2 + OH = H2O + HO2
!
      rcarr(33,:) = skarr(  1.800D-12 ,0.0D+00 ,temperature)
!
!....           HO2 + NO = NO2 + OH
!
      rcarr(34,:) = skarr(  3.300D-12 ,-270.0D+00 ,temperature)
!
!....           HO2 + HO2 = H2O2 + O2
!
      rcarr(35,:) = skho2dis (temperature ,adcol)
!
!....           H2O + HO2 + HO2 = H2O + H2O2 + O2
!
      rcarr(36,:) = skho2h2o (temperature ,adcol)
!
!....           H2 + OH = H + H2O
!
      rcarr(37,:) = skarr(  2.800D-12 ,1800.0D+00 ,temperature)
!
!....           CO + OH = H
!
      rcarr(38,:) = skcooh (temperature ,adcol)
!
!....           CH4 + OH = H2O + MO2
!
      rcarr(39,:) = skohch4 (temperature)
!
!....           MO2 + NO = CH2O + HO2 + NO2
!
      rcarr(40,:) = skarr(  2.800D-12 ,-300.0D+00 ,temperature)
!
!....           ClO + MO2 = CH2O + Cl + HO2 + O2
!
      rcarr(41,:) = skarr(  1.800D-12 ,600.0D+00 ,temperature)
!
!....           HO2 + MO2 = MP + O2
!
      rcarr(42,:) = skarr(  4.100D-13 ,-750.0D+00 ,temperature)
!
!....           MO2 + MO2 = CH2O + MOH + O2
!
      rcarr(43,:) = skmo2dis_1 (temperature)
!
!....           MO2 + MO2 = 2 CH2O + 2 HO2
!
      rcarr(44,:) = skmo2dis_2 (temperature)
!
!....           MP + OH = H2O + MO2
!
      rcarr(45,:) = skarr(  2.660D-12 ,-200.0D+00 ,temperature)
!
!....           MP + OH = CH2O + H2O + OH
!
      rcarr(46,:) = skarr(  1.140D-12 ,-200.0D+00 ,temperature)
!
!....           CH2O + OH = CO + H2O + HO2
!
      rcarr(47,:) = skarr(  5.500D-12 ,-125.0D+00 ,temperature)
!
!....           N + O2 = NO + O
!
      rcarr(48,:) = skarr(  1.500D-11 ,3600.0D+00 ,temperature)
!
!....           N + NO = N2 + O
!
      rcarr(49,:) = skarr(  2.100D-11 ,-100.0D+00 ,temperature)
!
!....           NO2 + O = NO + O2
!
      rcarr(50,:) = skarr(  5.100D-12 ,-210.0D+00 ,temperature)
!
!....           NO3 + O = NO2 + O2
!
      rcarr(51,:) = skarr(  1.000D-11 ,0.0D+00 ,temperature)
!
!....           NO2 + OH = HNO3
!
      rcarr(52,:) = sktroe(  1.800D-30 ,3.00D0 & 
     &                     , 2.800D-11 ,0.00D0 ,0.0D0 &
     &                     ,temperature ,adcol)
!
!....           HNO3 + OH = H2O + NO3
!
      rcarr(53,:) = skohhno3 (temperature ,adcol)
!
!....           NO + OH = HNO2
!
      rcarr(54,:) = sktroe(  7.000D-31 ,2.60D0 & 
     &                     , 3.600D-11 ,0.10D0 ,0.0D0 &
     &                     ,temperature ,adcol)
!
!....           HNO2 + OH = H2O + NO2
!
      rcarr(55,:) = skarr(  1.800D-11 ,390.0D+00 ,temperature)
!
!....           HO2 + NO2 = HNO4
!
      rcarr(56,:) = sktroe(  1.900D-31 ,3.40D0 & 
     &                     , 4.000D-12 ,0.30D0 ,0.0D0 &
     &                     ,temperature ,adcol)
!
!....           HNO4 = HO2 + NO2
!
      rcarr(57,:) = sktroe(  9.050D-05 ,3.4D0 &
     &                     , 1.900D+15 ,0.30D0 ,10900.0D+00 &
     &                     ,temperature ,adcol)
!
!....           HNO4 + OH = H2O + NO2 + O2
!
      rcarr(58,:) = skarr(  1.300D-12 ,-380.0D+00 ,temperature)
!
!....           HO2 + NO3 = NO2 + O2 + OH
!
      rcarr(59,:) = skarr(  3.500D-12 ,0.0D+00 ,temperature)
!
!....           NO + NO3 = 2 NO2
!
      rcarr(60,:) = skarr(  1.500D-11 ,-170.0D+00 ,temperature)
!
!....           NO3 + OH = HO2 + NO2
!
      rcarr(61,:) = skarr(  2.200D-11 ,0.0D+00 ,temperature)
!
!....           NO2 + NO3 = N2O5
!
      rcarr(62,:) = sktroe(  2.400D-30 ,3.00D0 & 
     &                     , 1.600D-12 ,-0.10D0 ,0.0D0 &
     &                     ,temperature ,adcol)
!
!....           N2O5 = NO2 + NO3
!
      rcarr(63,:) = sktroe(  4.140D-04 ,3.0D0 &
     &                     , 2.760D+14 ,-0.10D0 ,10840.0D+00 &
     &                     ,temperature ,adcol)
!
!....           HCOOH + OH = H2O + HO2
!
      rcarr(64,:) = skarr(  4.000D-13 ,0.0D+00 ,temperature)
!
!....           MOH + OH = CH2O + HO2
!
      rcarr(65,:) = skarr(  2.900D-12 ,345.0D+00 ,temperature)
!
!....           NO2 + NO3 = NO + NO2 + O2
!
      rcarr(66,:) = skarr(  4.500D-14 ,1260.0D+00 ,temperature)
!
!....           CH2O + NO3 = CO + HNO3 + HO2
!
      rcarr(67,:) = skarr(  5.800D-16 ,0.0D+00 ,temperature)
!
!....           Cl + O3 = ClO + O2
!
      rcarr(68,:) = skarr(  2.300D-11 ,200.0D+00 ,temperature)
!
!....           Cl + H2 = H + HCl
!
      rcarr(69,:) = skarr(  3.050D-11 ,2270.0D+00 ,temperature)
!
!....           Cl + H2O2 = HCl + HO2
!
      rcarr(70,:) = skarr(  1.100D-11 ,980.0D+00 ,temperature)
!
!....           Cl + HO2 = HCl + O2
!
      rcarr(71,:) = skarr(  1.400D-11 ,-270.0D+00 ,temperature)
!
!....           Cl + HO2 = ClO + OH
!
      rcarr(72,:) = skarr(  3.600D-11 ,375.0D+00 ,temperature)
!
!....           ClO + O = Cl + O2
!
      rcarr(73,:) = skarr(  2.800D-11 ,-85.0D+00 ,temperature)
!
!....           ClO + OH = Cl + HO2
!
      rcarr(74,:) = skarr(  7.400D-12 ,-270.0D+00 ,temperature)
!
!....           ClO + OH = HCl + O2
!
      rcarr(75,:) = skarr(  6.000D-13 ,-230.0D+00 ,temperature)
!
!....           ClO + HO2 = HOCl + O2
!
      rcarr(76,:) = skarr(  2.600D-12 ,-290.0D+00 ,temperature)
!
!....           ClO + HO2 = HCl + O3
!
      rcarr(77,:) = skarr(  0.000D+00 ,0.0D+00 ,temperature)
!
!....           ClO + NO = Cl + NO2
!
      rcarr(78,:) = skarr(  6.400D-12 ,-290.0D+00 ,temperature)
!
!....           ClO + NO2 = ClONO2
!
      rcarr(79,:) = sktroe(  1.800D-31 ,3.40D0 & 
     &                     , 1.500D-11 ,1.90D0 ,0.0D0 &
     &                     ,temperature ,adcol)
!
!....           ClO + ClO = 2 Cl + O2
!
      rcarr(80,:) = skarr(  3.000D-11 ,2450.0D+00 ,temperature)
!
!....           ClO + ClO = Cl2 + O2
!
      rcarr(81,:) = skarr(  1.000D-12 ,1590.0D+00 ,temperature)
!
!....           ClO + ClO = Cl + OClO
!
      rcarr(82,:) = skarr(  3.500D-13 ,1370.0D+00 ,temperature)
!
!....           ClO + ClO = Cl2O2
!
      rcarr(83,:) = sktroe(  1.900D-32 ,3.60D0 & 
     &                     , 3.700D-12 ,1.60D0 ,0.0D0 &
     &                     ,temperature ,adcol)
!
!....           Cl2O2 = 2 ClO
!
      rcarr(84,:) = sktroe(  8.800D-06 ,3.6D0 &
     &                     , 1.710D+15 ,1.60D0 ,8537.0D+00 &
     &                     ,temperature ,adcol)
!
!....           HCl + OH = Cl + H2O
!
      rcarr(85,:) = skarr(  1.800D-12 ,250.0D+00 ,temperature)
!
!....           HOCl + OH = ClO + H2O
!
      rcarr(86,:) = skarr(  3.000D-12 ,500.0D+00 ,temperature)
!
!....           ClONO2 + O = ClO + NO3
!
      rcarr(87,:) = skarr(  3.600D-12 ,840.0D+00 ,temperature)
!
!....           ClONO2 + OH = HOCl + NO3
!
      rcarr(88,:) = skarr(  1.200D-12 ,330.0D+00 ,temperature)
!
!....           Cl + ClONO2 = Cl2 + NO3
!
      rcarr(89,:) = skarr(  6.500D-12 ,-135.0D+00 ,temperature)
!
!....           Br + O3 = BrO + O2
!
      rcarr(90,:) = skarr(  1.600D-11 ,780.0D+00 ,temperature)
!
!....           Br + HO2 = HBr + O2
!
      rcarr(91,:) = skarr(  4.800D-12 ,310.0D+00 ,temperature)
!
!....           Br + CH2O = CO + HBr + HO2
!
      rcarr(92,:) = skarr(  1.700D-11 ,800.0D+00 ,temperature)
!
!....           BrO + O = Br + O2
!
      rcarr(93,:) = skarr(  1.900D-11 ,-230.0D+00 ,temperature)
!
!....           BrO + HO2 = HOBr + O2
!
      rcarr(94,:) = skarr(  4.500D-12 ,-460.0D+00 ,temperature)
!
!....           BrO + NO = Br + NO2
!
      rcarr(95,:) = skarr(  8.800D-12 ,-260.0D+00 ,temperature)
!
!....           BrO + NO2 = BrONO2
!
      rcarr(96,:) = sktroe(  5.400D-31 ,3.10D0 & 
     &                     , 6.500D-12 ,2.90D0 ,0.0D0 &
     &                     ,temperature ,adcol)
!
!....           BrO + ClO = Br + OClO
!
      rcarr(97,:) = skarr(  9.500D-13 ,-550.0D+00 ,temperature)
!
!....           BrO + ClO = Br + Cl + O2
!
      rcarr(98,:) = skarr(  2.300D-12 ,-260.0D+00 ,temperature)
!
!....           BrO + ClO = BrCl + O2
!
      rcarr(99,:) = skarr(  4.100D-13 ,-290.0D+00 ,temperature)
!
!....           BrO + BrO = 2 Br + O2
!
      rcarr(100,:) = skbrodis (temperature)
!
!....           HBr + OH = Br + H2O
!
      rcarr(101,:) = skarr(  5.500D-12 ,-200.0D+00 ,temperature)
!
!....           CH2O + O = CO + HO2 + OH
!
      rcarr(102,:) = skarr(  3.400D-11 ,1600.0D+00 ,temperature)
!
!....           CH4 + Cl = HCl + MO2
!
      rcarr(103,:) = skarr(  7.100D-12 ,1270.0D+00 ,temperature)
!
!....           CH2O + Cl = CO + HCl + HO2
!
      rcarr(104,:) = skarr(  8.100D-11 ,30.0D+00 ,temperature)
!
!....           CH3Cl + OH = Cl + H2O + HO2
!
      rcarr(105,:) = skarr(  1.960D-12 ,1200.0D+00 ,temperature)
!
!....           CH3CCl3 + OH = 3 Cl + H2O
!
      rcarr(106,:) = skarr(  1.640D-12 ,1520.0D+00 ,temperature)
!
!....           HCFC22 + OH = Cl + H2O
!
      rcarr(107,:) = skarr(  9.200D-13 ,1560.0D+00 ,temperature)
!
!....           HCFC141b + OH = 2 Cl + H2O
!
      rcarr(108,:) = skarr(  1.250D-12 ,1600.0D+00 ,temperature)
!
!....           HCFC142b + OH = Cl + H2O
!
      rcarr(109,:) = skarr(  1.300D-12 ,1770.0D+00 ,temperature)
!
!....           CH3Cl + Cl = CO + 2 HCl + HO2
!
      rcarr(110,:) = skarr(  2.030D-11 ,1110.0D+00 ,temperature)
!
!....           CH3Br + OH = Br + H2O + HO2
!
      rcarr(111,:) = skarr(  1.420D-12 ,1150.0D+00 ,temperature)
!
!....           A3O2 + HO2 = RA3P
!
      rcarr(112,:) = ska3o2_ho2 (temperature)
!
!....           A3O2 + MO2 =  0.75 CH2O + HO2 +  0.25 MOH +  0.75 RCHO +  0.25 ROH
!
      rcarr(113,:) = skarr(  5.920D-13 ,0.0D+00 ,temperature)
!
!....           A3O2 + NO = HO2 + NO2 + RCHO
!
      rcarr(114,:) = skarr(  2.900D-12 ,-350.0D+00 ,temperature)
!
!....           ACET + OH = ATO2 + H2O
!
      rcarr(115,:) = skacetoh (temperature)
!
!....           ACTA + OH = H2O + MO2
!
      rcarr(116,:) = skarr(  3.150D-14 ,-920.0D+00 ,temperature)
!
!....           ALD2 + NO3 = HNO3 + MCO3
!
      rcarr(117,:) = skarr(  1.400D-12 ,1900.0D+00 ,temperature)
!
!....           ALD2 + OH =  0.05 CH2O +  0.05 CO + H2O +  0.05 HO2 +  0.95 MCO3
!
      rcarr(118,:) = skarr(  4.630D-12 ,-350.0D+00 ,temperature)
!
!....           ALK4 + NO3 = HNO3 + R4O2
!
      rcarr(119,:) = skarr(  2.800D-12 ,3280.0D+00 ,temperature)
!
!....           ALK4 + OH = R4O2
!
      rcarr(120,:) = skarr(  9.100D-12 ,405.0D+00 ,temperature)
!
!....           ATO2 + HO2 = MCO3 + MO2
!
      rcarr(121,:) = skarr(  8.600D-13 ,-700.0D+00 ,temperature)
!
!....           ATO2 + MCO3 = ACTA + MEK
!
      rcarr(122,:) = skarr(  1.870D-13 ,-500.0D+00 ,temperature)
!
!....           ATO2 + MCO3 =  0.20 CH2O +  0.80 HO2 +  0.20 MCO3 +  0.80 MGLY + MO2
!
      rcarr(123,:) = skarr(  1.680D-12 ,-500.0D+00 ,temperature)
!
!....           ATO2 + MO2 =  0.50 CH2O +  0.20 HAC +  0.30 HO2 +  0.30 MCO3 +  0.50 MGLY +  0.50 MOH
!
      rcarr(124,:) = skarr(  7.500D-13 ,-500.0D+00 ,temperature)
!
!....           ATO2 + NO =  0.96 CH2O +  0.96 MCO3 +  0.96 NO2 +  0.04 R4N2
!
      rcarr(125,:) = skarr(  2.800D-12 ,-300.0D+00 ,temperature)
!
!....           B3O2 + HO2 = RB3P
!
      rcarr(126,:) = skb3o2_ho2 (temperature)
!
!....           B3O2 + MCO3 = ACET + ACTA
!
      rcarr(127,:) = skarr(  1.870D-13 ,-500.0D+00 ,temperature)
!
!....           B3O2 + MCO3 = ACET + HO2 + MO2
!
      rcarr(128,:) = skarr(  1.680D-12 ,-500.0D+00 ,temperature)
!
!....           B3O2 + MO2 =  0.75 ACET +  0.75 CH2O + HO2 +  0.25 MOH +  0.25 ROH
!
      rcarr(129,:) = skarr(  8.370D-14 ,0.0D+00 ,temperature)
!
!....           B3O2 + NO = ACET + HO2 + NO2
!
      rcarr(130,:) = skarr(  2.700D-12 ,-350.0D+00 ,temperature)
!
!....           C2H6 + NO3 = ETO2 + HNO3
!
      rcarr(131,:) = skarr(  1.400D-18 ,0.0D+00 ,temperature)
!
!....           C2H6 + OH = ETO2 + H2O
!
      rcarr(132,:) = skarr(  7.660D-12 ,1020.0D+00 ,temperature)
!
!....           C2H6 + Cl = ETO2 + HCl
!
      rcarr(133,:) = skarr(  7.200D-11 ,70.0D+00 ,temperature)
!
!....           C3H8 + OH = A3O2
!
      rcarr(134,:) = skc3h8ox_2 (temperature)
!
!....           C3H8 + OH = B3O2
!
      rcarr(135,:) = skc3h8ox_1 (temperature)
!
!....           EOH + OH = ALD2 + HO2
!
      rcarr(136,:) = skarr(  3.350D-12 ,0.0D+00 ,temperature)
!
!....           ETO2 + ETO2 = 2 ALD2 + 2 HO2
!
      rcarr(137,:) = skarr(  4.100D-14 ,0.0D+00 ,temperature)
!
!....           ETO2 + ETO2 = ALD2 + EOH
!
      rcarr(138,:) = skarr(  2.700D-14 ,0.0D+00 ,temperature)
!
!....           ETO2 + NO = ALD2 + HO2 + NO2
!
      rcarr(139,:) = skarr(  2.600D-12 ,-365.0D+00 ,temperature)
!
!....           ETP + OH =  0.50 ALD2 +  0.50 ETO2 +  0.50 OH
!
      rcarr(140,:) = skarr(  5.180D-12 ,-200.0D+00 ,temperature)
!
!....           GLYC + OH =  0.73 CH2O +  0.50 CO +  0.13 GLYX +  0.13 HCOOH +  0.77 HO2 +  0.23 OH
!
      rcarr(141,:) = skarr(  8.000D-12 ,0.0D+00 ,temperature)
!
!....           GLYC + OH = CO + HCOOH + OH
!
      rcarr(142,:) = skarr(  8.000D-12 ,0.0D+00 ,temperature)
!
!....           GLYX + NO3 = 2 CO + HNO3 + HO2
!
      rcarr(143,:) = skno3glyx (temperature ,oxygen)
!
!....           GLYX + OH = 2 CO + HO2
!
      rcarr(144,:) = skohglyx (temperature ,oxygen)
!
!....           HAC + OH = HO2 + MGLY
!
      rcarr(145,:) = skarr(  2.150D-12 ,-305.0D+00 ,temperature)
!
!....           HAC + OH =  0.50 ACTA +  0.50 CO +  0.50 HCOOH +  0.50 MO2 + OH
!
      rcarr(146,:) = skarr(  2.150D-12 ,-305.0D+00 ,temperature)
!
!....           ETO2 + HO2 = ETP
!
      rcarr(147,:) = skarr(  2.910D-13 ,-1300.0D+00 ,temperature)
!
!....           HO2 + MCO3 = ACTA + O3
!
      rcarr(148,:) = skho2mco3_1 (temperature)
!
!....           HO2 + MCO3 = MAP
!
      rcarr(149,:) = skho2mco3_2 (temperature)
!
!....           IALD + O3 =  0.12 CH2O +  0.28 GLYC +  0.20 GLYX +  0.20 HAC +  0.20 HCOOH +  0.60 MGLY +  0.30 O3 +  0.10 OH
!
      rcarr(150,:) = skarr(  6.160D-15 ,1814.0D+00 ,temperature)
!
!....           IALD + OH =  0.15 HO2 +  0.44 IAO2 +  0.41 MAO3
!
      rcarr(151,:) = skarr(  3.700D-11 ,0.0D+00 ,temperature)
!
!....           HO2 + IAO2 = IAP
!
      rcarr(152,:) = skarr(  2.910D-13 ,-1300.0D+00 ,temperature)
!
!....           IAO2 + MCO3 = ACTA + MEK
!
      rcarr(153,:) = skarr(  1.870D-13 ,-500.0D+00 ,temperature)
!
!....           IAO2 + MCO3 =  0.40 CH2O +  0.29 CO +  0.26 GLYC +  0.18 GLYX +  0.36 HAC + HO2 +  0.58 MGLY + MO2
!
      rcarr(154,:) = skarr(  1.680D-12 ,-500.0D+00 ,temperature)
!
!....           IAO2 + MO2 =  0.95 CH2O +  0.15 CO +  0.13 GLYC +  0.09 GLYX +  0.18 HAC + HO2 +  0.25 MEK +  0.29 MGLY +  0.25 MOH +  0.25 ROH
!
      rcarr(155,:) = skarr(  1.300D-12 ,0.0D+00 ,temperature)
!
!....           IAO2 + NO =  0.35 CH2O +  0.27 CO +  0.24 GLYC +  0.17 GLYX +  0.33 HAC +  0.08 HNO3 +  0.92 HO2 +  0.53 MGLY +  0.92 NO2
!
      rcarr(156,:) = skarr(  2.700D-12 ,-350.0D+00 ,temperature)
!
!....           IAP + OH =  0.50 IAO2 +  0.50 OH +  0.50 RCHO
!
      rcarr(157,:) = skarr(  5.310D-12 ,-200.0D+00 ,temperature)
!
!....           HO2 + INO2 = INPN
!
      rcarr(158,:) = skino2_ho2 (temperature)
!
!....           INO2 + MCO3 =  0.15 CH2O +  0.85 HNO3 +  0.80 HO2 +  0.10 MACR + MO2 +  0.05 MVK +  0.15 NO2
!
      rcarr(159,:) = skarr(  1.680D-12 ,-500.0D+00 ,temperature)
!
!....           INO2 + MCO3 = ACTA + NO2 + RCHO
!
      rcarr(160,:) = skarr(  1.870D-13 ,-500.0D+00 ,temperature)
!
!....           INO2 + MO2 =  0.83 CH2O +  0.43 HNO3 +  0.90 HO2 +  0.05 MACR +  0.25 MOH +  0.03 MVK +  0.57 NO2 +  0.25 RCHO +  0.25 ROH
!
      rcarr(161,:) = skarr(  1.300D-12 ,0.0D+00 ,temperature)
!
!....           INO2 + NO =  0.15 CH2O +  0.85 HNO3 +  0.80 HO2 +  0.10 MACR +  0.05 MVK +  1.15 NO2
!
      rcarr(162,:) = skarr(  2.700D-12 ,-350.0D+00 ,temperature)
!
!....           INPN + OH = INO2
!
      rcarr(163,:) = skarr(  5.180D-12 ,-200.0D+00 ,temperature)
!
!....           HO2 + ISN1 = ISNP
!
      rcarr(164,:) = skarr(  2.910D-13 ,-1300.0D+00 ,temperature)
!
!....           ISN1 + MCO3 = GLYC + HAC + MO2 + NO2
!
      rcarr(165,:) = skarr(  1.680D-12 ,-500.0D+00 ,temperature)
!
!....           ISN1 + MCO3 = ACTA + NO2 + RCHO
!
      rcarr(166,:) = skarr(  1.870D-13 ,-500.0D+00 ,temperature)
!
!....           ISN1 + MO2 =  0.75 CH2O +  0.50 GLYC +  0.50 HAC +  0.50 HO2 +  0.25 MOH + NO2 +  0.25 RCHO +  0.25 ROH
!
      rcarr(167,:) = skarr(  8.370D-14 ,0.0D+00 ,temperature)
!
!....           ISNP + OH =  0.50 ISN1 +  0.50 NO2 +  0.50 OH +  0.50 RCHO
!
      rcarr(168,:) = skarr(  4.750D-12 ,-200.0D+00 ,temperature)
!
!....           ISOP + NO3 = INO2
!
      rcarr(169,:) = skarr(  3.300D-12 ,450.0D+00 ,temperature)
!
!....           ISOP + O3 =  0.90 CH2O +  0.05 CO +  0.06 HO2 +  0.39 MACR +  0.16 MVK +  0.10 O3 +  0.27 OH +  0.07 PRPE
!
      rcarr(170,:) = skarr(  1.000D-14 ,1970.0D+00 ,temperature)
!
!....           ISOP + OH = RIO2
!
      rcarr(171,:) = skarr(  3.100D-11 ,-350.0D+00 ,temperature)
!
!....           HO2 + KO2 = MGLY + MO2
!
      rcarr(172,:) = skko2_ho2 (temperature)
!
!....           KO2 + MCO3 = ACTA + MEK
!
      rcarr(173,:) = skarr(  1.870D-13 ,-500.0D+00 ,temperature)
!
!....           KO2 + MCO3 = ALD2 + MCO3 + MO2
!
      rcarr(174,:) = skarr(  1.680D-12 ,-500.0D+00 ,temperature)
!
!....           KO2 + MO2 =  0.50 ALD2 +  0.75 CH2O +  0.50 HO2 +  0.50 MCO3 +  0.25 MEK +  0.25 MOH +  0.25 ROH
!
      rcarr(175,:) = skarr(  8.370D-14 ,0.0D+00 ,temperature)
!
!....           KO2 + NO =  0.93 ALD2 +  0.93 MCO3 +  0.93 NO2 +  0.07 R4N2
!
      rcarr(176,:) = skarr(  2.700D-12 ,-350.0D+00 ,temperature)
!
!....           MACR + NO3 = MAN2
!
      rcarr(177,:) = skarr(  2.300D-15 ,0.0D+00 ,temperature)
!
!....           MACR + NO3 = HNO3 + MAO3
!
      rcarr(178,:) = skarr(  1.100D-15 ,0.0D+00 ,temperature)
!
!....           MACR + O3 =  0.70 CH2O +  0.20 CO +  0.28 HO2 +  0.80 MGLY +  0.20 O3 +  0.22 OH
!
      rcarr(179,:) = skarr(  1.400D-15 ,2100.0D+00 ,temperature)
!
!....           MACR + OH =  0.53 MAO3 +  0.47 MRO2
!
      rcarr(180,:) = skarr(  8.000D-12 ,-380.0D+00 ,temperature)
!
!....           HO2 + MAN2 = ISNP
!
      rcarr(181,:) = skko2_ho2 (temperature)
!
!....           MAN2 + MCO3 = CH2O + MGLY + MO2 + NO2
!
      rcarr(182,:) = skarr(  1.680D-12 ,-500.0D+00 ,temperature)
!
!....           MAN2 + MCO3 = ACTA + NO2 + RCHO
!
      rcarr(183,:) = skarr(  1.870D-13 ,-500.0D+00 ,temperature)
!
!....           MAN2 + MO2 =  1.25 CH2O +  0.50 HO2 +  0.50 MGLY +  0.25 MOH + NO2 +  0.25 RCHO +  0.25 ROH
!
      rcarr(184,:) = skarr(  8.370D-14 ,0.0D+00 ,temperature)
!
!....           MAN2 + NO = CH2O + MGLY + 2 NO2
!
      rcarr(185,:) = skarr(  2.700D-12 ,-350.0D+00 ,temperature)
!
!....           HO2 + MAO3 =  0.59 CH2O +  0.39 CO +  0.41 MAOP +  0.39 MO2 +  0.15 O3 +  0.44 OH
!
      rcarr(186,:) = skarr(  4.300D-13 ,-1040.0D+00 ,temperature)
!
!....           MAO3 + MCO3 = CH2O + MCO3 + MO2
!
      rcarr(187,:) = skarr(  2.500D-12 ,-500.0D+00 ,temperature)
!
!....           MAO3 + MO2 = 2 CH2O + HO2 + MCO3
!
      rcarr(188,:) = skarr(  1.680D-12 ,-500.0D+00 ,temperature)
!
!....           MAO3 + MO2 = CH2O + RCOOH
!
      rcarr(189,:) = skarr(  1.870D-13 ,-500.0D+00 ,temperature)
!
!....           MAO3 + NO2 = PMN
!
      rcarr(190,:) = sktroe(  9.000D-28 ,8.90D0 & 
     &                     , 7.700D-12 ,0.20D0 ,0.0D0 &
     &                     ,temperature ,adcol)
!
!....           MAO3 + NO =  0.50 CH2O +  0.50 CO +  0.50 MCO3 +  0.50 MO2 + NO2
!
      rcarr(191,:) = skarr(  6.700D-12 ,-340.0D+00 ,temperature)
!
!....           MAOP + OH = MAO3
!
      rcarr(192,:) = skarr(  3.600D-12 ,-380.0D+00 ,temperature)
!
!....           A3O2 + MCO3 = ACTA + RCHO
!
      rcarr(193,:) = skarr(  1.870D-12 ,-500.0D+00 ,temperature)
!
!....           A3O2 + MCO3 = HO2 + MO2 + RCHO
!
      rcarr(194,:) = skarr(  1.680D-12 ,-500.0D+00 ,temperature)
!
!....           ETO2 + MCO3 = ACTA + ALD2
!
      rcarr(195,:) = skarr(  1.870D-13 ,-500.0D+00 ,temperature)
!
!....           ETO2 + MCO3 = ALD2 + HO2 + MO2
!
      rcarr(196,:) = skarr(  1.680D-12 ,-500.0D+00 ,temperature)
!
!....           MCO3 + MCO3 = 2 MO2
!
      rcarr(197,:) = skarr(  2.500D-12 ,-500.0D+00 ,temperature)
!
!....           MCO3 + MO2 = ACTA + CH2O
!
      rcarr(198,:) = skarr(  2.000D-13 ,-500.0D+00 ,temperature)
!
!....           MCO3 + MO2 = CH2O + HO2 + MO2
!
      rcarr(199,:) = skarr(  1.800D-12 ,-500.0D+00 ,temperature)
!
!....           MCO3 + NO2 = PAN
!
      rcarr(200,:) = sktroe(  9.700D-29 ,5.60D0 & 
     &                     , 9.300D-12 ,1.50D0 ,0.0D0 &
     &                     ,temperature ,adcol)
!
!....           MCO3 + NO = MO2 + NO2
!
      rcarr(201,:) = skarr(  8.100D-12 ,-270.0D+00 ,temperature)
!
!....           MCO3 + PO2 = ACTA +  0.65 HAC +  0.35 RCHO
!
      rcarr(202,:) = skarr(  1.870D-12 ,-500.0D+00 ,temperature)
!
!....           MCO3 + PO2 = ALD2 + CH2O + HO2 + MO2
!
      rcarr(203,:) = skarr(  1.680D-12 ,-500.0D+00 ,temperature)
!
!....           MEK + NO3 = HNO3 + KO2
!
      rcarr(204,:) = skarr(  8.000D-16 ,0.0D+00 ,temperature)
!
!....           MEK + OH = H2O + KO2
!
      rcarr(205,:) = skohmek (temperature)
!
!....           MGLY + NO3 = CO + HNO3 + MCO3
!
      rcarr(206,:) = skarr(  3.360D-12 ,1860.0D+00 ,temperature)
!
!....           MGLY + OH = CO + MCO3
!
      rcarr(207,:) = skarr(  1.500D-11 ,0.0D+00 ,temperature)
!
!....           ETO2 + MO2 =  0.75 ALD2 +  0.75 CH2O +  0.25 EOH + HO2 +  0.25 MOH
!
      rcarr(208,:) = skarr(  3.000D-13 ,0.0D+00 ,temperature)
!
!....           HO2 + MRO2 = MRP
!
      rcarr(209,:) = skko2_ho2 (temperature)
!
!....           MCO3 + MRO2 = ACTA + MEK
!
      rcarr(210,:) = skarr(  1.870D-13 ,-500.0D+00 ,temperature)
!
!....           MCO3 + MRO2 =  0.17 CH2O +  0.83 CO +  0.83 HAC + HO2 +  0.17 MGLY + MO2
!
      rcarr(211,:) = skarr(  1.680D-12 ,-500.0D+00 ,temperature)
!
!....           MO2 + MRO2 = CH2O +  0.60 CO +  0.60 HAC + HO2 +  0.25 MGLY +  0.15 ROH
!
      rcarr(212,:) = skarr(  8.370D-14 ,0.0D+00 ,temperature)
!
!....           MRO2 + NO = HNO3
!
      rcarr(213,:) = skro2noadd_3 (temperature ,adcol)
!
!....           MRO2 + NO =  0.17 CH2O +  0.83 CO +  0.83 HAC + HO2 +  0.17 MGLY + NO2
!
      rcarr(214,:) = skro2noabs_3 (temperature ,adcol)
!
!....           MRP + OH = MRO2
!
      rcarr(215,:) = skarr(  1.840D-12 ,-200.0D+00 ,temperature)
!
!....           MVK + O3 =  0.04 ALD2 +  0.80 CH2O +  0.05 CO +  0.06 HO2 +  0.82 MGLY +  0.20 O3 +  0.08 OH
!
      rcarr(216,:) = skarr(  8.500D-16 ,1520.0D+00 ,temperature)
!
!....           MVK + OH = VRO2
!
      rcarr(217,:) = skarr(  2.600D-12 ,-610.0D+00 ,temperature)
!
!....           MAP + OH =  0.50 CH2O +  0.50 MCO3 +  0.50 OH
!
      rcarr(218,:) = skarr(  6.130D-13 ,-200.0D+00 ,temperature)
!
!....           OH + RCHO = H2O + RCO3
!
      rcarr(219,:) = skarr(  6.000D-12 ,-410.0D+00 ,temperature)
!
!....           PAN = MCO3 + NO2
!
      rcarr(220,:) = skpanan (temperature,adcol)
!
!....           PMN = MAO3 + NO2
!
      rcarr(221,:) = skppndecomp (temperature,adcol)
!
!....           O3 + PMN =  0.60 CH2O + HO2 + NO2
!
      rcarr(222,:) = skarr(  8.200D-18 ,0.0D+00 ,temperature)
!
!....           OH + PMN = CO + HAC + NO2
!
      rcarr(223,:) = skarr(  2.900D-11 ,0.0D+00 ,temperature)
!
!....           HO2 + PO2 = PP
!
      rcarr(224,:) = ska3o2_ho2 (temperature)
!
!....           MO2 + PO2 =  0.50 ALD2 + CH2O +  0.16 HAC + HO2 +  0.25 MOH +  0.09 RCHO +  0.25 ROH
!
      rcarr(225,:) = skarr(  5.920D-13 ,0.0D+00 ,temperature)
!
!....           NO + PO2 = ALD2 + CH2O + HO2 + NO2
!
      rcarr(226,:) = skarr(  2.700D-12 ,-350.0D+00 ,temperature)
!
!....           PPN = NO2 + RCO3
!
      rcarr(227,:) = skppndecomp (temperature,adcol)
!
!....           OH + PP =  0.79 HAC +  0.79 OH +  0.21 PO2
!
      rcarr(228,:) = skarr(  8.780D-12 ,-200.0D+00 ,temperature)
!
!....           HO2 + PRN1 = PRPN
!
      rcarr(229,:) = ska3o2_ho2 (temperature)
!
!....           MCO3 + PRN1 = ALD2 + CH2O + MO2 + NO2
!
      rcarr(230,:) = skarr(  1.680D-12 ,-500.0D+00 ,temperature)
!
!....           MCO3 + PRN1 = ACTA + NO2 + RCHO
!
      rcarr(231,:) = skarr(  1.870D-13 ,-500.0D+00 ,temperature)
!
!....           MO2 + PRN1 =  0.50 ALD2 +  1.25 CH2O +  0.50 HO2 +  0.25 MOH + NO2 +  0.25 RCHO +  0.25 ROH
!
      rcarr(232,:) = skarr(  8.370D-14 ,0.0D+00 ,temperature)
!
!....           NO + PRN1 = ALD2 + CH2O + 2 NO2
!
      rcarr(233,:) = skarr(  2.700D-12 ,-350.0D+00 ,temperature)
!
!....           NO3 + PRPE = PRN1
!
      rcarr(234,:) = skarr(  4.590D-13 ,1156.0D+00 ,temperature)
!
!....           O3 + PRPE =  0.50 ALD2 +  0.54 CH2O +  0.42 CO +  0.06 H2 +  0.30 HO2 +  0.31 MO2 +  0.14 OH
!
      rcarr(235,:) = skarr(  5.500D-15 ,1880.0D+00 ,temperature)
!
!....           OH + PRPE = PO2
!
      rcarr(236,:) = sktroe(  8.000D-27 ,3.50D0 & 
     &                     , 3.000D-11 ,1.00D0 ,0.0D0 &
     &                     ,temperature ,adcol)
!
!....           OH + PRPN = PRN1
!
      rcarr(237,:) = skarr(  8.780D-12 ,-200.0D+00 ,temperature)
!
!....           HO2 + R4N1 = R4N2
!
      rcarr(238,:) = skarr(  7.400D-13 ,-700.0D+00 ,temperature)
!
!....           MCO3 + R4N1 =  0.75 ALD2 +  0.39 CH2O + MO2 + NO2 +  0.30 R4O2 +  0.57 RCHO
!
      rcarr(239,:) = skarr(  1.680D-12 ,-500.0D+00 ,temperature)
!
!....           MCO3 + R4N1 = ACTA + NO2 + RCHO
!
      rcarr(240,:) = skarr(  1.870D-13 ,-500.0D+00 ,temperature)
!
!....           MO2 + R4N1 =  0.38 ALD2 +  0.95 CH2O +  0.50 HO2 +  0.25 MOH + NO2 +  0.15 R4O2 +  0.54 RCHO +  0.25 ROH
!
      rcarr(241,:) = skarr(  8.370D-14 ,0.0D+00 ,temperature)
!
!....           NO + R4N1 =  0.75 ALD2 +  0.39 CH2O + 2 NO2 +  0.30 R4O2 +  0.57 RCHO
!
      rcarr(242,:) = skarr(  2.700D-12 ,-350.0D+00 ,temperature)
!
!....           OH + R4N2 = H2O + R4N1
!
      rcarr(243,:) = skarr(  1.600D-12 ,0.0D+00 ,temperature)
!
!....           HO2 + R4O2 = R4P
!
      rcarr(244,:) = skarr(  7.400D-13 ,-700.0D+00 ,temperature)
!
!....           MCO3 + R4O2 = ACTA + MEK
!
      rcarr(245,:) = skarr(  1.870D-13 ,-500.0D+00 ,temperature)
!
!....           MCO3 + R4O2 =  0.05 A3O2 +  0.32 ACET +  0.32 ALD2 +  0.18 B3O2 +  0.32 ETO2 +  0.27 HO2 +  0.19 MEK +  1.18 MO2 +  0.13 RCHO
!
      rcarr(246,:) = skarr(  1.680D-12 ,-500.0D+00 ,temperature)
!
!....           MO2 + R4O2 =  0.03 A3O2 +  0.16 ACET +  0.16 ALD2 +  0.09 B3O2 +  0.75 CH2O +  0.16 ETO2 +  0.64 HO2 +  0.35 MEK +  0.09 MO2 +  0.25 MOH +  0.07 RCHO +  0.25 ROH
!
      rcarr(247,:) = skarr(  8.370D-14 ,0.0D+00 ,temperature)
!
!....           NO + R4O2 =  0.05 A3O2 +  0.32 ACET +  0.32 ALD2 +  0.18 B3O2 +  0.32 ETO2 +  0.27 HO2 +  0.19 MEK +  0.18 MO2 + NO2 +  0.13 RCHO
!
      rcarr(248,:) = skro2noabs_1 (temperature ,adcol)
!
!....           NO + R4O2 = R4N2
!
      rcarr(249,:) = skro2noadd_1 (temperature ,adcol)
!
!....           OH + R4P =  0.50 OH +  0.50 R4O2 +  0.50 RCHO
!
      rcarr(250,:) = skarr(  8.780D-12 ,-200.0D+00 ,temperature)
!
!....           OH + RA3P =  0.50 A3O2 +  0.50 OH +  0.50 RCHO
!
      rcarr(251,:) = skarr(  5.180D-12 ,-200.0D+00 ,temperature)
!
!....           OH + RB3P =  0.79 ACET +  0.21 B3O2 +  0.79 OH
!
      rcarr(252,:) = skarr(  8.780D-12 ,-200.0D+00 ,temperature)
!
!....           NO3 + RCHO = HNO3 + RCO3
!
      rcarr(253,:) = skarr(  6.500D-15 ,0.0D+00 ,temperature)
!
!....           OH + RCOOH = ETO2
!
      rcarr(254,:) = skarr(  1.200D-12 ,0.0D+00 ,temperature)
!
!....           HO2 + RCO3 =  0.30 O3 +  0.30 RCOOH +  0.70 RP
!
      rcarr(255,:) = skarr(  4.300D-13 ,-1040.0D+00 ,temperature)
!
!....           MCO3 + RCO3 = ETO2 + MO2
!
      rcarr(256,:) = skarr(  2.500D-12 ,-500.0D+00 ,temperature)
!
!....           MO2 + RCO3 = CH2O + ETO2 + HO2
!
      rcarr(257,:) = skarr(  1.680D-12 ,-500.0D+00 ,temperature)
!
!....           MO2 + RCO3 = CH2O + RCOOH
!
      rcarr(258,:) = skarr(  1.870D-13 ,-500.0D+00 ,temperature)
!
!....           NO2 + RCO3 = PPN
!
      rcarr(259,:) = sktroe(  9.000D-28 ,8.90D0 & 
     &                     , 7.700D-12 ,0.20D0 ,0.0D0 &
     &                     ,temperature ,adcol)
!
!....           NO + RCO3 = ETO2 + NO2
!
      rcarr(260,:) = skarr(  6.700D-12 ,-340.0D+00 ,temperature)
!
!....           HO2 + RIO1 = RIP
!
      rcarr(261,:) = skarr(  2.910D-13 ,-1300.0D+00 ,temperature)
!
!....           MCO3 + RIO1 = ACTA + MEK
!
      rcarr(262,:) = skarr(  1.870D-13 ,-500.0D+00 ,temperature)
!
!....           MCO3 + RIO1 =  0.75 CH2O + HO2 + IALD + MO2
!
      rcarr(263,:) = skarr(  1.680D-12 ,-500.0D+00 ,temperature)
!
!....           MO2 + RIO1 =  1.13 CH2O + HO2 +  0.50 IALD +  0.25 MEK +  0.25 MOH +  0.25 ROH
!
      rcarr(264,:) = skarr(  8.370D-14 ,0.0D+00 ,temperature)
!
!....           NO + RIO1 = HNO3
!
      rcarr(265,:) = skro2noadd_2 (temperature ,adcol)
!
!....           NO + RIO1 =  0.75 CH2O + HO2 + IALD + NO2
!
      rcarr(266,:) = skro2noabs_2 (temperature ,adcol)
!
!....           HO2 + RIO2 = RIP
!
      rcarr(267,:) = skino2_ho2 (temperature)
!
!....           MCO3 + RIO2 = ACTA + MEK
!
      rcarr(268,:) = skarr(  1.870D-13 ,-500.0D+00 ,temperature)
!
!....           MCO3 + RIO2 =  0.69 CH2O +  0.86 HO2 +  0.13 IALD +  0.29 MACR + MO2 +  0.40 MVK +  0.14 RIO1
!
      rcarr(269,:) = skarr(  1.680D-12 ,-500.0D+00 ,temperature)
!
!....           MO2 + RIO2 =  1.10 CH2O +  0.93 HO2 +  0.06 IALD +  0.14 MACR +  0.25 MEK +  0.25 MOH +  0.20 MVK +  0.07 RIO1 +  0.25 ROH
!
      rcarr(270,:) = skarr(  8.370D-14 ,0.0D+00 ,temperature)
!
!....           NO + RIO2 = HNO3
!
      rcarr(271,:) = skro2noadd_2 (temperature ,adcol)
!
!....           NO + RIO2 =  0.69 CH2O +  0.86 HO2 +  0.13 IALD +  0.29 MACR +  0.40 MVK + NO2 +  0.14 RIO1
!
      rcarr(272,:) = skro2noabs_2 (temperature ,adcol)
!
!....           OH + RIP =  0.50 IAO2 +  0.10 RIO1 +  0.40 RIO2
!
      rcarr(273,:) = skarr(  4.750D-12 ,-200.0D+00 ,temperature)
!
!....           OH + ROH = HO2 + RCHO
!
      rcarr(274,:) = skarr(  4.600D-12 ,-70.0D+00 ,temperature)
!
!....           OH + RP =  0.50 ALD2 +  0.50 OH +  0.50 RCO3
!
      rcarr(275,:) = skarr(  6.130D-13 ,-200.0D+00 ,temperature)
!
!....           HO2 + VRO2 =  0.10 CH2O +  0.58 GLYC +  0.10 HO2 +  0.58 MCO3 +  0.19 MEK +  0.10 MGLY +  0.68 OH +  0.03 RCHO +  0.10 VRP
!
      rcarr(276,:) = skko2_ho2 (temperature)
!
!....           MCO3 + VRO2 = ACTA + MEK
!
      rcarr(277,:) = skarr(  1.870D-13 ,-500.0D+00 ,temperature)
!
!....           MCO3 + VRO2 =  0.28 CH2O +  0.72 GLYC +  0.28 HO2 +  0.72 MCO3 +  0.28 MGLY + MO2
!
      rcarr(278,:) = skarr(  1.680D-12 ,-500.0D+00 ,temperature)
!
!....           MO2 + VRO2 =  0.89 CH2O +  0.36 GLYC +  0.64 HO2 +  0.36 MCO3 +  0.25 MEK +  0.14 MGLY +  0.25 MOH +  0.25 ROH
!
      rcarr(279,:) = skarr(  8.370D-14 ,0.0D+00 ,temperature)
!
!....           NO + VRO2 = HNO3
!
      rcarr(280,:) = skro2noadd_3 (temperature ,adcol)
!
!....           NO + VRO2 =  0.28 CH2O +  0.72 GLYC +  0.28 HO2 +  0.72 MCO3 +  0.28 MGLY + NO2
!
      rcarr(281,:) = skarr(  2.700D-12 ,-350.0D+00 ,temperature)
!
!....           OH + VRP =  0.50 OH +  0.50 RCHO +  0.50 VRO2
!
      rcarr(282,:) = skarr(  8.780D-12 ,-200.0D+00 ,temperature)
!
!....           N2O5 = 2 HNO3
!
      rcarr(283,:) = sksts_n2o5 (temperature  & 
     &           ,pressure ,sad_lbs ,ptrop)
!
!....           ClONO2 = HNO3 + HOCl
!
      rcarr(284,:) = sksts_clono2 (temperature  & 
     &           ,adcol ,pressure ,sad_lbs ,specarr(  37,:) ,water  & 
     &           ,ptrop)
!
!....           BrONO2 = HNO3 + HOBr
!
      rcarr(285,:) = sksts_brono2 (temperature  & 
     &           ,pressure ,sad_lbs ,ptrop)
!
!....           ClONO2 + HCl = Cl2 + HNO3
!
      rcarr(286,:) = sksts_clono2_hcl (temperature  & 
     &           ,adcol ,pressure ,sad_lbs ,specarr(    36,:)  & 
     &           ,specarr(  37,:) ,water ,ptrop)
!
!....           HCl + HOCl = Cl2 + H2O
!
      rcarr(287,:) = sksts_hocl_hcl (temperature  & 
     &           ,adcol ,pressure ,sad_lbs ,specarr(  38,:)  & 
     &           ,specarr( 37,:) ,water ,ptrop)
!
!....           HCl + HOBr = BrCl + H2O
!
      rcarr(288,:) = sksts_hobr_hcl (temperature  & 
     &           ,adcol ,pressure ,sad_lbs ,specarr(  31,:)  & 
     &           ,specarr( 37,:) ,water ,ptrop)
!
!....           N2O5 = 2 HNO3
!
      rcarr(289,:) = sksts_n2o5 (temperature  & 
     &           ,pressure ,sad_sts ,ptrop)
!
!....           ClONO2 = HNO3 + HOCl
!
      rcarr(290,:) = sksts_clono2 (temperature  & 
     &           ,adcol ,pressure ,sad_sts ,specarr(  37,:) ,water  & 
     &           ,ptrop)
!
!....           BrONO2 = HNO3 + HOBr
!
      rcarr(291,:) = sksts_brono2 (temperature  & 
     &           ,pressure ,sad_sts ,ptrop)
!
!....           ClONO2 + HCl = Cl2 + HNO3
!
      rcarr(292,:) = sksts_clono2_hcl (temperature  & 
     &           ,adcol ,pressure ,sad_sts ,specarr(    36,:)  & 
     &           ,specarr(  37,:) ,water ,ptrop)
!
!....           HCl + HOCl = Cl2 + H2O
!
      rcarr(293,:) = sksts_hocl_hcl (temperature  & 
     &           ,adcol ,pressure ,sad_sts ,specarr(  38,:)  & 
     &           ,specarr( 37,:) ,water ,ptrop)
!
!....           HCl + HOBr = BrCl + H2O
!
      rcarr(294,:) = sksts_hobr_hcl (temperature  & 
     &           ,adcol ,pressure ,sad_sts ,specarr(  31,:)  & 
     &           ,specarr( 37,:) ,water ,ptrop)
!
!....           ClONO2 = HNO3 + HOCl
!
      rcarr(295,:) = sknat_clono2 (temperature  & 
     &           ,pressure ,sad_nat ,ptrop)
!
!....           BrONO2 = HNO3 + HOBr
!
      rcarr(296,:) = sknat_brono2 (temperature  & 
     &           ,pressure ,sad_nat ,ptrop)
!
!....           ClONO2 + HCl = Cl2 + HNO3
!
      rcarr(297,:) = sknat_hcl_clono2 (temperature  & 
     &           ,pressure ,sad_nat ,specarr(  37,:) ,ptrop)
!
!....           HCl + HOCl = Cl2 + H2O
!
      rcarr(298,:) = sknat_hcl_hocl (temperature  & 
     &           ,pressure ,sad_nat ,specarr(  37,:) ,ptrop)
!
!....           BrONO2 + HCl = BrCl + HNO3
!
      rcarr(299,:) = sknat_hcl_brono2 (temperature  & 
     &           ,pressure ,sad_nat ,specarr(  37,:) ,ptrop)
!
!....           HCl + HOBr = BrCl + H2O
!
      rcarr(300,:) = sknat_hcl_hobr (temperature  & 
     &           ,pressure ,sad_nat ,specarr(  37,:) ,ptrop)
!
!....           ClONO2 = HNO3 + HOCl
!
      rcarr(301,:) = skice_clono2 (temperature  & 
     &           ,pressure ,sad_ice ,ptrop)
!
!....           BrONO2 = HNO3 + HOBr
!
      rcarr(302,:) = skice_brono2 (temperature  & 
     &           ,pressure ,sad_ice ,ptrop)
!
!....           ClONO2 + HCl = Cl2 + HNO3
!
      rcarr(303,:) = skice_hcl_clono2 (temperature  & 
     &           ,pressure ,sad_ice ,specarr(  37,:) ,ptrop)
!
!....           HCl + HOCl = Cl2 + H2O
!
      rcarr(304,:) = skice_hcl_hocl (temperature  & 
     &           ,pressure ,sad_ice ,specarr(  37,:) ,ptrop)
!
!....           BrONO2 + HCl = BrCl + HNO3
!
      rcarr(305,:) = skice_hcl_brono2 (temperature  & 
     &           ,pressure ,sad_ice ,specarr(  37,:) ,ptrop)
!
!....           HCl + HOBr = BrCl + H2O
!
      rcarr(306,:) = skice_hcl_hobr (temperature  & 
     &           ,pressure ,sad_ice ,specarr(  37,:) ,ptrop)
!
!....           HNO3 = NO2 + OH
!
      rcarr(307,:) = sksoot_hno3 (temperature ,sad_soot)
!
!....           NO3 + NO3 = 2 NO2 + O2
!
      rcarr(308,:) = skarr(  8.500D-13 ,2450.0D+00 ,temperature)
!
!....           HO2 =  0.50 H2O
!
      rcarr(309,:) = sktrs_ho2 (temperature, & 
     &           sadcol2,adcol,radA,NSADaer,NSADdust,cPBLcol,pressure)
!
!....           NO2 =  0.50 HNO2 +  0.50 HNO3
!
      rcarr(310,:) = sktrs_no2 (temperature, & 
     &            sadcol2, adcol, radA, NSADaer,NSADdust,ptrop & 
     &           , pressure)
!
!....           NO3 = HNO3
!
      rcarr(311,:) = sktrs_no3 (temperature, & 
     &            sadcol2, adcol, radA,NSADaer,NSADdust,ptrop & 
     &           , pressure)
!
!....           N2O5 = 2 HNO3
!
      rcarr(312,:) = sktrs_n2o5 (temperature, & 
     &           sadcol2,adcol,radA,FRH,NSADaer,NSADdust, ptrop & 
     &           , pressure)
!
!....          End thermal rate constants
!
      CONTAINS
        FUNCTION skarr (af,ae,tk)
          real*8 &
     &      af ,ae ,tk(:)
          real*8, dimension(size(tk)) :: skarr
          skarr(:) = af * exp(-ae / tk(:))
        END FUNCTION skarr
        FUNCTION sklp (af,npwr,tk,ad)
          real*8 &
     &      af ,npwr ,tk(:) ,ad(:)
          real*8, dimension(size(tk)) :: sklp
          sklp(:) = ad(:) * af * (300.0d0/tk(:))**npwr
        END FUNCTION sklp
        FUNCTION skhp (ai,mpwr,tk)
          real*8 &
     &      ai ,mpwr ,tk(:)
          real*8, dimension(size(tk)) :: skhp
          skhp(:) = ai * (300.0d0/tk(:))**mpwr
        END FUNCTION skhp
        FUNCTION skfo (af,npwr,ai,mpwr,tk,ad)
          real*8 &
     &      af ,npwr ,ai ,mpwr ,tk(:) ,ad(:)
          real*8, dimension(size(tk)) :: skfo
          skfo(:) = sklp(af,npwr,tk,ad) / skhp(ai,mpwr,tk)
        END FUNCTION skfo
        FUNCTION skterlp (af,npwr,ae,tk,ad)
          real*8 &
     &      af ,npwr ,ae ,tk(:) ,ad(:)
          real*8, dimension(size(tk)) :: skterlp
          skterlp(:) = sklp(af,npwr,tk,ad) * exp(-ae / tk(:))
        END FUNCTION skterlp
        FUNCTION sktroe (af,npwr,ai,mpwr,ae,tk,ad,fc)
          real*8 &
     &      af ,npwr ,ai ,mpwr ,ae ,tk(:) ,ad(:)
          real*8, OPTIONAL :: fc
          real*8, dimension(size(tk)) :: sktroe
          real*8 fsubc
          real*8 skfo_local(size(tk))
          if (present (fc)) then; fsubc=fc; else; fsubc=0.6d0; end if
          skfo_local(:) = skfo(af,npwr,ai,mpwr,tk,ad)
          sktroe(:) = skterlp(af,npwr,ae,tk,ad) * fsubc** &
     &                (1.0d0/(1.0d0+log10(skfo_local(:))**2)) / &
     &                (1.0d0+skfo(af,npwr,ai,mpwr,tk,ad))
        END FUNCTION sktroe
!
!.... skho2dis (temperature ,adcol)
!
!_1_
!
!.... JPL 10-6
!
        FUNCTION skho2dis (tk,ad)
!
!....       HO2 + HO2 = H2O2 + O2
!....
!.... This routine returns a bimolecular rate constant that accounts
!.... for the pressure, but not the H2O, dependence of the reaction.
!.... The H2O dependence is treated separately.
!
          real*8  tk(:) ,ad(:)
          real*8, DIMENSION (size(tk)) :: skho2dis
          skho2dis(:) = 3.0d-13 * exp(460.0d0 / tk(:)) +  &
     &                  2.1d-33 * ad(:) * exp(920.0d0 / tk(:))
        END FUNCTION skho2dis
!
!.... skho2h2o (temperature ,adcol)
!
!_2_
!
!.... JPL 10-6
!
        FUNCTION skho2h2o (tk,ad)
!
!....       HO2 + HO2 + H2O = H2O2 + O2 + H2O
! NOT USED : multiplication factor is given in B13
!            not right here as there is H2O dependence
!....
!
          real*8  tk(:) ,ad(:)
          real*8, DIMENSION (size(tk)) :: skho2h2o
          skho2h2o(:) = skho2dis(tk ,ad) * 1.4d-21 * exp(2200.d0 / tk(:))
        END FUNCTION skho2h2o
!
!.... skcooh (temperature ,adcol)
!
!_1_
!
!.... JPL 10-6 ; now CO + OH is composed of two separate reactions!
!.... Now it's density and temperature dependent.
!         M
! OH + CO -> HOCO , but HOCO + O2 -> HO2 + CO2 quickly ; termolecular
!         M
! OH + CO -> H + CO2, but H + O2 -> HO2 quickly ; chemical activation reaction
!
        FUNCTION skcooh (tk,ad)
!
!                  M
!....      1) OH + CO = H + CO2; assume H+O2->HO2 quick
!....      2) OH + CO = HO2 + CO2
!....
!.... Pressure in hPa
!
          real*8  tk(:), ad(:), af ,npwr ,ai ,mpwr ,ae
          real*8, DIMENSION (size(tk)) :: skcooh
          real*8, DIMENSION (size(tk)) :: skcoohlp
          real*8, DIMENSION (size(tk)) :: skcoohhp
          real*8, DIMENSION (size(tk)) :: skcooh1
          real*8, DIMENSION (size(tk)) :: skcooh2
          real*8 fsubc
!
! 1)
          fsubc = 0.6d0
          af = 1.5d-13
          npwr = -0.6
          ai = 2.1d9
          mpwr = -6.1
!
          skcoohlp(:) = af * (300.0d0/tk(:))**npwr
          skcoohhp(:) = skhp(ai,mpwr,tk)
          skcooh2(:)  =  skcoohlp(:)/( skcoohhp(:)/ad(:) )
          skcooh(:)   = ( skcoohlp(:)/(1.0d0+skcooh2(:)) ) * fsubc**  &
     &                  ( 1.0d0/( 1.0d0+(log10(skcooh2(:)))**2 ) )
!
          fsubc = 0.6d0
          af = 5.9d-33
          npwr = 1.4
          ai = 1.1d-12
          mpwr = -1.3
!
! 2)
          skcooh(:) = skcooh(:) + sktroe(af,npwr,ai,mpwr,0.0D0,tk(:),ad(:))
!
        END FUNCTION skcooh
!
!.... skohch4 (temperature)
! _1_
!
!.... Harvard/GMI JPL 10-6
!
        FUNCTION skohch4 (tk)
!
!....      OH + CH4 = MO2 + H2O
!
          real*8  tk(:)
          real*8, DIMENSION(size(tk)) :: skohch4
!
          skohch4(:) = 2.80D-14 * tk(:)**0.667d0 * exp(-1575.0d0 / tk(:))
!
        END FUNCTION skohch4
!
!.... skmo2dis_1 (temperature)
! _2_
!
!.... Harvard/GMI
!
        FUNCTION skmo2dis_1 (tk)
!
!....      MO2 + MO2 = MOH + CH2O + O2
!....
!======================================================================
!
          real*8  tk(:)
          real*8, DIMENSION(size(tk)) :: skmo2dis_1
!
          skmo2dis_1(:) = 9.50D-14 * exp(390.0d0 / tk(:)) * (1.0d0 /  &
     &                    (1.0d0 + 26.2d0 * exp(-1130.0d0 / tk(:))))
!
        END FUNCTION skmo2dis_1
!
!.... skmo2dis_2 (temperature)
! _3_
!
!.... Harvard/GMI
!
        FUNCTION skmo2dis_2 (tk)
!
!....      MO2 + MO2 = 2 CH2O + 2 HO2
!....
!======================================================================
!
          real*8  tk(:)
          real*8, DIMENSION(size(tk)) :: skmo2dis_2
!
          skmo2dis_2(:) = 9.50D-14 * exp(390.0d0 / tk(:)) * (1.0d0 /  &
     &                    (1.0d0 + 0.04d0 * exp(1130.0d0 / tk(:))))
!
        END FUNCTION skmo2dis_2
!
!.... skohhno3 (temperature ,adcol)
!
!_1_
!
!.... JPL 10-6
!
        FUNCTION skohhno3 (tk,ad)
!
!                    M
!....      OH + HNO3 = NO3 + H2O
!....
          real*8  tk(:) ,ad(:)
          real*8, DIMENSION (size(tk)) :: skohhno3
          real*8  xxx1(size(tk)) ,xxx2(size(tk))
!
          xxx1(:)     = ad(:) * 6.5d-34 * exp(1335.0d0 / tk(:))
          xxx2(:)     = 2.7d-17 * exp(2199.0d0 / tk(:))
          skohhno3(:) = 2.4d-14 * exp(460.0d0 / tk(:)) + xxx1(:) / (1.0d0 + xxx1(:) / xxx2(:))
!
        END FUNCTION skohhno3
!
!.... skbrodis (temperature)
!
!_1_
!
!.... JPL 10-6 (unchanged from JPL 00-003)
!
        FUNCTION skbrodis (tk)
!
!....       BrO + BrO = 2 Br + O2
!....
!.... PSC - 8/23/2002
!....
!.... Combined two product channels into one.
!.... Appears as special function to avoid confusion with
!.... actual elementary reaction in reaction database.
!
          real*8  tk(:)
          real*8, DIMENSION (size(tk)) :: skbrodis
          skbrodis(:) = 1.5D-12 * exp(230.0d0 / tk(:))
        END FUNCTION skbrodis
!
!.... ska3o2_ho2 (temperature)
! _27_
!
!.... Harvard/GMI
!
        FUNCTION ska3o2_ho2 (tk)
!
!.... A3O2 + HO2 = RA3P
!
          real*8  tk(:)
          real*8, DIMENSION(size(tk)) :: ska3o2_ho2
!
! A3O2 +  HO2 => RA3P : 
!A  486 2.91E-13  0.0E+00   1300 1 HR  0.00     0.     0.         
!       3.00E+00  0.0E+00      0 0     0.00     0.     0.         
!
       ska3o2_ho2(:) = 2.91D-13 * exp(1300.0d0 / tk(:)) / (1.0D0+3.00E+00)
!
        END FUNCTION ska3o2_ho2
!
!.... skacetoh (temperature)
! _24_
!
!.... Harvard/GMI JPL 10-6
!
        FUNCTION skacetoh (tk)
!
!....      ACET + OH = ATO2 + H2O
!....
!======================================================================
!
          real*8  tk(:)
          real*8, dimension(size(tk)) :: skacetoh
!
      skacetoh(:) = skarr(3.8200D-11 ,2000.0D+00 ,tk(:)) +1.330D-13
!
        END FUNCTION skacetoh
!
!.... skb3o2_ho2 (temperature)
! _26_
!
!.... Harvard/GMI
!
        FUNCTION skb3o2_ho2 (tk)
!
!.... B3O2 + HO2 = RB3P
!
          real*8  tk(:)
          real*8, DIMENSION(size(tk)) :: skb3o2_ho2
!
! B3O2 +  HO2 => RB3P : 
!A  472 2.91E-13  0.0E+00   1300 1 HR  0.00     0.     0.         
!       3.00E+00  0.0E+00      0 0     0.00     0.     0.         
!
       skb3o2_ho2(:) = 2.91D-13 * exp(1300.0d0 / tk(:)) / (1.0D0+3.00E+00)
!
        END FUNCTION skb3o2_ho2
!
!.... skc3h8ox_2 (temperature)
! _23_
!
!.... Harvard/GMI
!
        FUNCTION skc3h8ox_2 (tk)
!
! ***********************NOT USED********************************
!        
!....      C3H8 + OH = A3O2
!....
!======================================================================
!
          real*8  tk(:)
          real*8, dimension(size(tk)) :: skc3h8ox_2
!
          skc3h8ox_2(:) = 7.60D-12 * exp(-585.0d0 / tk(:)) * (1.0d0 /  &
     &                    (1.0d0 + 0.17d0 * (300.0d0/tk(:))**(-0.64d0) * exp(816.0d0 / tk(:))))
!
        END FUNCTION skc3h8ox_2
!
!.... skc3h8ox_1 (temperature)
! _22_
!
!.... Harvard/GMI
!
        FUNCTION skc3h8ox_1 (tk)
!
! ***********************NOT USED********************************
!
!....      C3H8 + OH = B3O2
!....
!======================================================================
!
          real*8  tk(:)
          real*8, dimension(size(tk)) :: skc3h8ox_1
!
          skc3h8ox_1(:) = 7.60D-12 * exp(-585.0d0 / tk(:)) * (1.0d0 /  &
     &                    (1.0d0 + 5.87d0 * (300.0d0/tk(:))**(0.64d0) *  exp(-816.0d0 / tk(:))))
!
        END FUNCTION skc3h8ox_1
!
!.... skno3glyx (temperature ,oxygen)
! _5_
!
!.... Harvard/GMI
!
        FUNCTION skno3glyx (tk,o2)
!
!                     O2
!....      NO3 + GLYX = HO2 + 2 CO
!....
!.... PSC - 8/8/2002
!....
!
          real*8  tk(:) ,o2(:)
          real*8, DIMENSION(size(tk)) :: skno3glyx
!
          skno3glyx(:) = 1.40D-12 * exp(-1860.0d0 / tk(:)) *  &
     &                   (o2(:) + 3.5D+18) / (2.0d0 * o2(:) + 3.5D+18)
!
        END FUNCTION skno3glyx
!
!.... skohglyx (temperature ,oxygen)
! _4_
!
!.... Harvard/GMI
!
        FUNCTION skohglyx (tk,o2)
!
!                    O2
!....      OH + GLYX = HO2 + 2 CO
!....
!.... PSC - 8/8/2002
!....
!
          real*8  tk(:) ,o2(:)
          real*8, DIMENSION(size(tk)) :: skohglyx
!
          skohglyx(:) = 1.10D-11 * (o2(:) + 3.5D+18) / (2.0d0 * o2(:) + 3.5D+18)
!
        END FUNCTION skohglyx
!
!.... skho2mco3_1 (temperature)
! _7_
!
!.... Harvard/GMI JPL 10-6
!
        FUNCTION skho2mco3_1 (tk)
!
!....      HO2 + MCO3 = ACTA + O3
!
          real*8  tk(:)
          real*8, DIMENSION(size(tk)) :: skho2mco3_1
!
          skho2mco3_1(:) = 4.30D-13 * exp(1040.0d0 / tk(:)) * (1.0d0 /  &
     &                     (1.0d0 + 37.0d0 * exp(-660.0d0 / tk(:))))
!
        END FUNCTION skho2mco3_1
!
!.... skho2mco3_2 (temperature)
! _8_
!
!.... Harvard/GMI JPS 10-6
!
        FUNCTION skho2mco3_2 (tk)
!
!....      HO2 + MCO3 = MAP
!
          real*8  tk(:)
          real*8, DIMENSION(size(tk)) :: skho2mco3_2
!
          skho2mco3_2(:) = 4.30D-13 * exp(1040.0d0 / tk(:)) * (1.0d0 /  &
     &                     (1.0d0 + 2.70D-02 * exp(660.0d0 / tk(:))))
!
        END FUNCTION skho2mco3_2
!
!.... skino2_ho2 (temperature)
! _28_
!
!.... Harvard/GMI
!
        FUNCTION skino2_ho2 (tk)
!
!.... INO2 + HO2 = INPN
!
          real*8  tk(:)
          real*8, DIMENSION(size(tk)) :: skino2_ho2
!
! INO2 +  HO2 => INPN : 
!A  472 2.91E-13  0.0E+00   1300 1 HR  0.00     0.     0.         
!       5.00E+00  0.0E+00      0 0     0.00     0.     0.         
!
       skino2_ho2(:) = 2.91D-13 * exp(1300.0d0 / tk(:)) / (1.0D0+5.00E+00)
!
        END FUNCTION skino2_ho2
!
!.... skko2_ho2 (temperature)
! _29_
!
!.... Harvard/GMI
!
        FUNCTION skko2_ho2 (tk)
!
!.... KO2 + HO2 = MO2 + MGLY
!.... MAN2 + HO2 = ISNP
!.... MRO2 + HO2 = MRP 
!.... O2 + HO2 = MO2 + MGLY
!
          real*8  tk(:)
          real*8, DIMENSION(size(tk)) :: skko2_ho2
!
! KO2 +  HO2 => MO2 + MGLY (ours)
! KO2 +  HO2 => 0.15 OH + 0.15 ALD2 + 0.15 MCO3 + 0.85 ATOOH  (GEOSCHEM)
!A  472 2.91E-13  0.0E+00   1300 1 HR  0.00     0.     0.         
!       4.00E+00  0.0E+00      0 0     0.00     0.     0.         
!
       skko2_ho2(:) = 2.91D-13 * exp(1300.0d0 / tk(:)) / (1.0D0+4.00E+00)
!
        END FUNCTION skko2_ho2
!
!.... skohmek (temperature)
! _6_
!
!.... Harvard/GMI
!
        FUNCTION skohmek (tk)
!
!....      OH + MEK = KO2 + H2O
!....
!.... PSC - 8/8/2002
!....
!
          real*8  tk(:)
          real*8, DIMENSION(size(tk)) :: skohmek
!
          skohmek(:) = 2.92D-12 * (300.0d0 / tk(:))**(-2.0d0) * exp(414.0d0 / tk(:))
!
        END FUNCTION skohmek
!
!.... skro2noadd_3 (temperature ,adcol)
! _18_
!
!.... Harvard/GMI
!
        FUNCTION skro2noadd_3 (tk,ad)
!
!.... VRO2 + NO = HNO3
!.... MRO2 + NO = HNO3
!....
!.... PSC - 8/8/2002
!.... SDS - 7/6/2016 updated rate
!....
!
          real*8  tk(:) ,ad(:)
          real*8, DIMENSION(size(tk)) :: skro2noadd_3
!
          skro2noadd_3(:) = 2.70D-12 * exp(350.0d0 / tk(:)) * fyrno3(4.0d0,tk,ad)
!
        END FUNCTION skro2noadd_3
!
!.... skro2noabs_3 (temperature ,adcol)
! _15_
!
!.... Harvard/GMI
!
        FUNCTION skro2noabs_3 (tk,ad)
!
!.... VRO2 + NO = NO2 + 0.280 HO2 + 0.280 CH2O +
!....                   0.720 MCO3 + 0.720 GLYC +
!....                   0.280 MGLY
!.... MRO2 + NO = NO2 + HO2 + 0.170 MGLY + 0.830 HAC +
!....                         0.830 CO + 0.170 CH2O
!....
!.... PSC - 8/8/2002
!.... SDS - 7/6/2016 updated rate
!....
!
          real*8  tk(:) ,ad(:)
          real*8, DIMENSION(size(tk)) :: skro2noabs_3
!
          skro2noabs_3(:) = 2.70D-12 * exp(350.0d0 / tk(:)) * (1.0d0 - fyrno3(4.0d0,tk,ad))
!
        END FUNCTION skro2noabs_3
!
!.... skpanan (temperature,adcol)
! _19_
!
!.... Harvard/GMI JPL 10-6
!
        FUNCTION skpanan (tk,ad)
!
!.... PAN = MCO3 + NO2
!....
          real*8  tk(:),ad(:)
          real*8, DIMENSION(size(tk)) :: skpanan
!
       skpanan(:) =  sktroe(9.700D-29,5.60D0,9.300D-12,1.50D0,0.0D0,tk(:),ad(:))/  &
     &  skarr(9.000D-29, -14000.0D0,tk(:))
!
        END FUNCTION skpanan
!
!.... skppndecomp (temperature,adcol)
! _25_
!
!.... Harvard/GMI JPL 10-6
!
        FUNCTION skppndecomp (tk,ad)
!
!....     PPN   = NO2 + RCO3
!....     GPAN  = NO2 + GCO3
!....     PMN   = NO2 + MAO3
!....
!======================================================================
!
          real*8  tk(:), ad(:)
          real*8, dimension(size(tk)) :: skppndecomp
!
      skppndecomp(:) = sktroe(9.000D-28, 8.9d0,7.700D-12, 0.2d0, 0.0d0,  &
     &      tk(:), ad(:))/skarr(9.00D-29 ,-14000.0D+00, tk(:))
!
        END FUNCTION skppndecomp
!
!.... skro2noabs_1 (temperature ,adcol)
! _13_
!
!.... Harvard/GMI
!
        FUNCTION fyrno3(xcarbn,tk,ad)
          real*8  xcarbn
          real*8  tk(:) ,ad(:)
          real*8, DIMENSION(size(tk)) :: fyrno3
          real*8  aaa(size(tk)) ,rarb(size(tk)) ,xxyn(size(tk)) ,yyyn(size(tk)) ,zzyn(size(tk))
!
          xxyn(:)   = 1.94D-22 * exp(0.97d0 * xcarbn) * ad(:) * (300.0d0 / tk(:))**0.0d0
          yyyn(:)   = 0.826d0 * (300.0d0 / tk(:))**8.1d0
          aaa(:)    = log10(xxyn(:) / yyyn(:))
          zzyn(:)   = 1.0d0 / (1.0d0 + aaa(:) * aaa(:))
          rarb(:)   = (xxyn(:) / (1.0d0 + (xxyn(:) / yyyn(:)))) * 0.411d0**zzyn(:)
          fyrno3(:) = rarb(:) / (1.0d0 + rarb(:))
!
        END FUNCTION fyrno3
        FUNCTION skro2noabs_1 (tk,ad)
!
!....      R4O2 + NO = NO2 + 0.320 ACET + 0.190 MEK +
!....                        0.180 MO2 + 0.270 HO2 +
!....                        0.320 ALD2 + 0.130 RCHO +
!....                        0.050 A3O2 + 0.180 B3O2 +
!....                        0.320 ETO2
!....
!.... PSC - 8/8/2002
!.... SDS - 7/6/2016 updated rate
!....
!
          real*8  tk(:) ,ad(:)
          real*8, DIMENSION(size(tk)) :: skro2noabs_1
!
          skro2noabs_1(:) = 2.70D-12 * exp(350.0d0 / tk(:)) * (1.0d0 - fyrno3(4.5d0,tk,ad))
!
        END FUNCTION skro2noabs_1
!
!.... skro2noadd_1 (temperature ,adcol)
! _16_
!
!.... Harvard/GMI
!
        FUNCTION skro2noadd_1 (tk,ad)
!
!.... R4O2 + NO = R4N2
!....
!.... PSC - 8/8/2002
!.... SDS - 7/6/2016 updated rate
!....
!
          real*8  tk(:) ,ad(:)
          real*8, DIMENSION(size(tk)) :: skro2noadd_1
!
          skro2noadd_1(:) = 2.70D-12 * exp(350.0d0 / tk(:)) * fyrno3(4.5d0,tk,ad)
!
        END FUNCTION skro2noadd_1
!
!.... skro2noadd_2 (temperature ,adcol)
! _17_
!
!.... Harvard/GMI
!
        FUNCTION skro2noadd_2 (tk,ad)
!
!.... RIO2 + NO = HNO3
!.... RIO1 + NO = HNO3
!....
!.... PSC - 8/8/2002
!.... SDS - 7/6/2016 updated rate
!....
!
          real*8  tk(:) ,ad(:)
          real*8, DIMENSION(size(tk)) :: skro2noadd_2
!
          skro2noadd_2(:) = 2.70D-12 * exp(350.0d0 / tk(:)) * fyrno3(5.0d0,tk,ad)
!
        END FUNCTION skro2noadd_2
!
!.... skro2noabs_2 (temperature ,adcol)
! _14_
!
!.... Harvard/GMI
!
        FUNCTION skro2noabs_2 (tk,ad)
!
!.... RIO2 + NO = NO2 + 0.864 HO2 + 0.690 CH2O +
!....                   0.402 MVK + 0.288 MACR +
!....                   0.136 RIO1 + 0.127 IALD
!.... RIO1 + NO = NO2 + IALD + HO2 + 0.750 CH2O
!....
!.... PSC - 8/8/2002
!.... SDS - 7/6/2016 updated rate
!....
!
          real*8  tk(:) ,ad(:)
          real*8, DIMENSION(size(tk)) :: skro2noabs_2
!
          skro2noabs_2(:) = 2.70D-12 * exp(350.0d0 / tk(:)) * (1.0d0 - fyrno3(5.0d0,tk,ad))
!
        END FUNCTION skro2noabs_2
!
!.... sksts_n2o5 (temperature ,pressure ,sad_lbs ,ptrop)
!
!_1_
!
!.... JPL 97-4
!
        FUNCTION sksts_n2o5 (tk,pr,sad,ptrop)
          real*8, OPTIONAL :: ptrop
          real*8  tk(:) ,pr(:) ,sad(:)
          real*8, DIMENSION (size(tk)) :: sksts_n2o5
          real*8  gamma ,pi
          real*8  avgvel(size(tk))
!
          pi = acos(-1.0d0)
!
!=======================================================================
!     N2O5 + stratospheric sulfate aerosol = 2 HNO3
!=======================================================================
!
!.... First order reaction rate constant
!.... PSC 3/30/99
!
          gamma    = 0.10d0
!
          avgvel(:)= 100.0d0 * (8.0d0 * 8.31448d0 * tk(:) * 1000.0d0 /  &
     &               (pi * mw(IN2O5)))**0.5d0
!
          where( sad > 0.0d0 )
            sksts_n2o5   = 0.25d0 * gamma * avgvel * sad
          elsewhere
            sksts_n2o5   = 0.0d0
          end where
!
          if ( present(ptrop) ) then
            where( pr > ptrop ) sksts_n2o5 = 0.0d0
          end if
!
        END FUNCTION sksts_n2o5
!
!.... sksts_clono2 (temperature ,adcol ,pressure ,sad_lbs ,specarr( HCl,:) ,water ,ptrop)
!
!_2_
!
!.... JPL 97-4
!
        FUNCTION sksts_clono2 (tk,ad,pr,sad,hcl,h2o,ptrop)
          real*8, OPTIONAL :: ptrop
          real*8  tk(:), ad(:), pr(:), sad(:), h2o(:), hcl(:)
          real*8, DIMENSION (size(tk)) :: sksts_clono2
          real*8  adrop ,alpha ,ksur ,minconc ,pi ,ro
          real*8  adivl(size(tk)) ,ah2o(size(tk)) ,avgvel(size(tk))  &
     &     ,fterm(size(tk))  &
     &     ,gamma(size(tk)) ,gam0(size(tk)) ,gcalc(size(tk))  &
     &     ,gprob_hcl(size(tk)) ,gprob_tot(size(tk)) ,gsurf(size(tk))  &
     &     ,hstar(size(tk))  &
     &     ,ph2o(size(tk)) ,phcl(size(tk)) ,prate(size(tk))  &
     &     ,tk_150(size(tk))
!
          pi = acos(-1.0d0)
!
!=======================================================================
!     ClONO2 + stratospheric sulfate aerosol = HOCl + HNO3
!=======================================================================
!
!.... First order reaction rate constant
!.... PSC 3/30/99
!
!....   NOTE: alpha is modified from 0.3 in Table 2, JPC,
!....         Hanson and Ravi, 98, 5734, 1994 to 1.0 based on
!....         Ravi and Hanson, 101, pg 3887, JGR, 1996.
!
          minconc = 1.0d0
          alpha   = 1.0d0
          ksur    = 576.0d0
          ro      = 2000.0d0
          adrop   = 1.0d-05
!
!....    NOTE: Partial pressure of HCl and H2O (in atmospheres)
!
          ph2o(:) = ((h2o(:) / ad(:)) * pr(:)) * (1.0d0 / 1013.25d0)
!
          where( hcl(:) <= minconc )
            phcl  = ((1.0d0 / ad(:)) * pr(:)) * (1.0d0 / 1013.25d0)
          elsewhere
            phcl  = ((hcl(:) / ad(:)) * pr(:)) * (1.0d0 / 1013.25d0)
          end where
!
!....    NOTE: Activity of H2O not allowed to exceed 1.1
!....          ah2o and Hstar taken from Table 2, JPC,
!....          Hanson and Ravi, 98, 5734, 1994
!
          ah2o(:) = 1013.25d0 * ph2o(:) / 10.0d0**  &
     &                       (9.217d0 - (2190.0d0 / (tk(:) - 12.7d0)))
!
          where( ah2o(:) > 1.1d0)
            ah2o = 1.1d0
          end where
!
          tk_150(:) = tk(:)
!
          where( tk(:) < 150.0d0 )
            tk_150 = 150.0d0
          endwhere
!
          hstar(:) = exp((6250.0d0 / tk_150(:)) - 10.414d0) * (ah2o(:)**3.49d0)
!
          gsurf(:) = ah2o(:) * ksur * hstar(:) * phcl(:)
!
          prate(:) = ro * hstar(:) * phcl(:) / ah2o(:)
!
          gam0(:)  = 1.18d-04 + (9.1d-03 * ah2o(:)) + (0.5d0 * ah2o(:)**2.0d0)
!
          gcalc(:) = gam0(:) * sqrt(1.0d0 + prate(:))
!
          adivl(:) = adrop / (1.4d-06 * sqrt(1.0d0 / ah2o(:)))
!
          fterm(:) = ((exp(adivl(:)) + exp(-adivl(:))) /  &
     &                (exp(adivl(:)) - exp(-adivl(:)))) -  &
     &                (1.0d0 / adivl(:))
!
!....   NOTE: gprob_tot is the overall uptake coeff for ClONO2
!
          gprob_tot(:) = 1.0d0 / (1.0d0 / (gsurf(:) + fterm(:) * gcalc(:)) + 1.0d0 / alpha)
!
          gprob_hcl(:) = gprob_tot(:) *  &
     &                   (gsurf(:) +  &
     &                    fterm(:) * gcalc(:) * prate(:) /  &
     &                   (1.0d0 + prate(:))) /  &
     &                    (gsurf(:) + fterm(:) * gcalc(:))
!
          gamma(:) = gprob_tot(:) - gprob_hcl(:)
!
          avgvel(:) = 100.0d0 * (8.0d0 * 8.31448d0 * tk(:) * 1000.0d0 / (pi * mw(ICLONO2)))**0.5d0
!
          sksts_clono2 = 0.25d0 * gamma * avgvel * sad
!
          where( sad < 0.0d0 ) sksts_clono2  = 0.0d0
!
          if ( present(ptrop) ) then
            where( pr > ptrop ) sksts_clono2 = 0.0d0
          end if
!
        END FUNCTION sksts_clono2
!
!.... sksts_brono2 (temperature ,pressure ,sad_lbs ,ptrop)
!
!_3_
!
!.... JPL 97-4
!
        FUNCTION sksts_brono2 (tk,pr,sad,ptrop)
          real*8, OPTIONAL :: ptrop
          real*8  tk(:) ,pr(:) ,sad(:)
          real*8, DIMENSION (size(tk)) :: sksts_brono2
          real*8  gamma ,pi
          real*8  avgvel(size(tk))
!
          pi = acos(-1.0d0)
!
!=======================================================================
!     BrONO2 + stratospheric sulfate aerosol = HOBr + HNO3
!=======================================================================
!
!.... First order reaction rate constant
!
!.... David Hanson, personal communication, May 13, 1997
!
          avgvel = 100.0d0 * (8.0d0 * 8.31448d0 * tk(:) * 1000.0d0 /  &
     &              (pi * mw(IBRONO2)))**0.5d0
!
          gamma  = 0.8d0
!
          sksts_brono2 = 0.25d0 * gamma * avgvel * sad
!
          where( sad < 0.0d0 ) sksts_brono2   = 0.0d0
!
          if ( present(ptrop) ) then
            where( pr > ptrop ) sksts_brono2 = 0.0d0
          end if
!
        END FUNCTION sksts_brono2
!
!.... sksts_clono2_hcl (temperature ,adcol ,pressure ,sad_lbs ,specarr(ClONO2,:) ,specarr( HCl,:) ,water ,ptrop)
!
!_4_
!
!.... JPL 97-4
!
        FUNCTION sksts_clono2_hcl (tk,ad,pr,sad,clono2,hcl,h2o,ptrop)
          real*8, OPTIONAL :: ptrop
          real*8  tk(:) ,ad(:) ,pr(:) ,sad(:) ,h2o(:) ,hcl(:) ,clono2(:)
          real*8, DIMENSION (size(tk)) :: sksts_clono2_hcl
          real*8  adrop ,alpha ,ksur ,minconc ,pi ,ro
          real*8  adivl(size(tk)) ,ah2o(size(tk)) ,avgvel(size(tk))  &
     &     ,fterm(size(tk))  &
     &     ,gam0(size(tk)) ,gcalc(size(tk))  &
     &     ,gprob_hcl(size(tk)) ,gprob_tot(size(tk)) ,gsurf(size(tk))  &
     &     ,hstar(size(tk))  &
     &     ,ph2o(size(tk)) ,phcl(size(tk)) ,prate(size(tk))  &
     &     ,tk_150(size(tk))
!
          pi = acos(-1.0d0)
!
!=======================================================================
!     ClONO2 + HCl on stratospheric sulfate aerosol = Cl2 + HNO3
!=======================================================================
!
!.... First order reaction rate constant
!
!.... Hanson and Ravi, JPC, 98, 5728, 1994
!.... DEK, 1/10/97
!
!....   NOTE: alpha is modified from 0.3 in Table 2, JPC,
!....         Hanson and Ravi, 98, 5734, 1994 to 1.0 based on
!....         Ravi and Hanson, 101, pg 3887, JGR, 1996.
!
          minconc = 1.0d0
          alpha   = 1.0d0
          ksur    = 576.0d0
          ro      = 2000.0d0
          adrop   = 1.0d-05
!
!....    NOTE: Partial pressure of HCl and H2O (in atmospheres)
!
          ph2o(:) = ((h2o(:)/ad(:)) * pr(:)) * (1.0d0/1013.25d0)
!
          where( hcl(:) <= minconc )
            phcl = ((1.0d0 / ad(:)) * pr(:)) * (1.0d0 / 1013.25d0)
          elsewhere
            phcl = ((hcl(:) / ad(:)) * pr(:)) * (1.0d0 / 1013.25d0)
          end where
!
!....    NOTE: Activity of H2O not allowed to exceed 1.1
!....          ah2o and Hstar taken from Table 2, JPC,
!....          Hanson and Ravi, 98, 5734, 1994
!
          ah2o(:) = 1013.25d0 * ph2o(:) / 10.0d0**(9.217d0 - (2190.0d0 / (tk(:) - 12.7d0)))
!
          where( ah2o(:) > 1.1d0)
            ah2o = 1.1d0
          end where
!
          tk_150(:) = tk(:)
!
          where( tk(:) < 150.0d0 )
            tk_150 = 150.0d0
          endwhere
!
          hstar(:) = exp((6250.0d0 / tk_150(:)) - 10.414d0) * (ah2o(:)**3.49d0)
!
          gsurf(:) = ah2o(:) * ksur * hstar(:) * phcl(:)
!
          prate(:) = ro * hstar(:) * phcl(:) / ah2o(:)
!
          gam0(:)  = 1.18d-04 + (9.1d-03 * ah2o(:)) + (0.5d0 * ah2o(:)**2.0d0)
!
          gcalc(:) = gam0(:) * sqrt(1.0d0 + prate(:))
!
          adivl(:) = adrop / (1.4d-06 * sqrt(1.0d0 / ah2o(:)))
!
          fterm(:) = ((exp(adivl(:)) + exp(-adivl(:))) /  &
     &                (exp(adivl(:)) - exp(-adivl(:)))) - (1.0d0 / adivl(:))
!
!....   NOTE: gprob_tot is the overall uptake coeff for ClONO2
!
!....   NOTE: alpha is modified from 0.3 in Table 2, JPC,
!....         Hanson and Ravi, 98, 5734, 1994 to 1.0 based on
!....         Ravi and Hanson, 101, pg 3887, JGR, 1996.
!
          gprob_tot(:) = 1.0d0 / (1.0d0 / (gsurf(:) + fterm(:) * gcalc(:)) +  &
     &                    1.0d0 / alpha)
!
          gprob_hcl(:) = gprob_tot(:) * (gsurf(:) + fterm(:) * gcalc(:) * prate(:) /  &
     &                   (1.0d0 + prate(:))) / (gsurf(:) + fterm(:) * gcalc(:))
!
          avgvel(:)    = 100.0d0 * (8.0d0 * 8.31448d0 * tk(:) * 1000.0d0 /  &
     &                   (pi * mw(ICLONO2)))**0.5d0
!
          where( hcl > minconc )
            sksts_clono2_hcl = 0.25d0 * gprob_hcl * avgvel * sad / hcl
          elsewhere
            sksts_clono2_hcl = 0.25d0 * gprob_hcl * avgvel * sad
          end where
!
          where( sad < 0.0d0 ) sksts_clono2_hcl = 0.0d0
!
          if ( present(ptrop) ) then
            where( pr > ptrop ) sksts_clono2_hcl = 0.0d0
          end if
!
        END FUNCTION sksts_clono2_hcl
!
!.... sksts_hocl_hcl (temperature ,adcol ,pressure ,sad_lbs ,specarr(HOCl,:) ,specarr(HCl,:) ,water ,ptrop)
!
!_5_
!
!.... JPL 97-4
!
        FUNCTION sksts_hocl_hcl (tk,ad,pr,sad,hocl,hcl,h2o,ptrop)
          real*8, OPTIONAL :: ptrop
          real*8  tk(:) ,ad(:) ,pr(:) ,sad(:) ,h2o(:) ,hcl(:) ,hocl(:)
          real*8, DIMENSION (size(tk)) :: sksts_hocl_hcl
          real*8  adrop ,alpha ,d1 ,minconc ,pi,minadivl
          real*8  adivl(size(tk)) ,ah2o(size(tk)) ,avgvel(size(tk))  &
     &     ,c1(size(tk)) ,c2(size(tk)) ,c3(size(tk)) ,conv(size(tk))  &
     &     ,fterm(size(tk))  &
     &     ,gcalc(size(tk))  &
     &     ,gprob_tot(size(tk))  &
     &     ,hhuth(size(tk)) ,hm(size(tk))  &
     &     ,hsqrtd(size(tk)) ,hstar(size(tk)) ,hstar_hocl(size(tk))  &
     &     ,k(size(tk)) ,kii(size(tk))  &
     &     ,mterm(size(tk))  &
     &     ,ph2o(size(tk)) ,phcl(size(tk))  &
     &     ,rho(size(tk))  &
     &     ,tk_150(size(tk))  &
     &     ,wtper(size(tk))  &
     &     ,z(size(tk))
!
          pi = acos(-1.0d0)
!
!=======================================================================
!     HOCl + HCl on stratospheric sulfate aerosol = Cl2 + H2O
!=======================================================================
!
!.... First order reaction rate constant
!
!.... Hanson and Ravi, JPC, 98, 5728, 1994
!.... DEK, 1/10/97
!
!....   NOTE: alpha is modified from 0.3 in Table 2, JPC,
!....         Hanson and Ravi, 98, 5734, 1994 to 1.0 based on
!....         Ravi and Hanson, 101, pg 3887, JGR, 1996.
!
          minconc  = 1.0d0
          adrop    = 1.0D-05
          alpha    = 1.0d0
          d1       = 9.0D-09
          minadivl = 1.00D-15
!
!....    NOTE: Partial pressure of HCl and H2O (in atmospheres)
!
          ph2o(:) = (h2o(:) / ad(:)) * pr(:)
!
          z(:) = log(ph2o(:))
!
          ph2o(:) = ph2o(:) / 1013.25d0
!
!
          where( hcl(:) <= minconc )
            phcl = ((1.0d0 / ad(:)) * pr(:)) * (1.0d0 / 1013.25d0)
          elsewhere
            phcl = ((hcl(:) / ad(:)) * pr(:)) * (1.0d0 / 1013.25d0)
          end where
!
!....    NOTE: Activity of H2O not allowed to exceed 1.1
!....          ah2o and Hstar taken from Table 2, JPC,
!....          Hanson and Ravi, 98, 5734, 1994
!
          ah2o(:) = 1013.25d0 * ph2o(:) / 10.0d0**(9.217d0 - (2190.0d0 / (tk(:) - 12.7d0)))
!
          where( ah2o(:) > 1.1d0)
            ah2o = 1.1d0
          end where
!
          wtper(:) = ((-14.0508d0 + 0.708928d0 * z(:)) * tk(:) + 3578.6d0) /  &
     &               (45.5374d0 + 1.55981d0 * z(:) - 0.197298d0 * tk(:))
!
          where( wtper < 40.0d0 )
            wtper = 40.0d0
          end where
!
          where( wtper > 80.0d0 )
            wtper = 80.0d0
          end where
!
          kii(:) = exp(2.303d0 * (6.08d0 - 1050.0d0 / tk(:) + 0.0747d0 * wtper(:)))
!
          tk_150(:) = tk(:)
!
          where( tk(:) < 150.0d0 )
            tk_150 = 150.0d0
          endwhere
!
          hstar(:) = exp((6250.0d0 / tk_150(:)) - 10.414d0) * (ah2o(:)**3.49d0)
!
          k(:)     = kii(:) * hstar(:) * phcl(:)
!
          adivl(:) = adrop / sqrt(d1 / k(:))
!
           where( adivl(:) < minadivl)
              adivl = minadivl
           endwhere
!
          fterm(:) = ((exp(adivl(:)) + exp(-adivl(:))) /  &
     &                (exp(adivl(:)) - exp(-adivl(:)))) -  &
     &               (1.0d0 / adivl(:))
!
          mterm(:) = 10.196d0 * wtper(:) / (100.0d0 - wtper(:))
!
          c1(:)    = 123.64d0 - 5.6D-04 * tk(:)**2.0d0
!
          c2(:)    = -29.54d0 + 1.814D-04 * tk(:)**2.0d0
!
          c3(:)    = 2.243d0 - 1.487D-03 * tk(:) + 1.324D-05 * tk(:)**2.0d0
!
          rho(:)   = 1000.0d0 + c1(:) * mterm(:) +  &
     &                          c2(:) * mterm(:)**1.5d0 +  &
     &                          c3(:) * mterm(:)**2.0d0
!
          conv(:)  = (rho(:) / 1000.0d0) / (1.0d0 + mterm(:) * 0.09808d0)
!
          hhuth(:) = exp(6.4946d0 - mterm(:) * (-0.04107d0 + 54.56d0 / tk(:)) -  &
     &                   5862.0d0 * (1.0d0 / 298.15d0 - 1.0d0 / tk(:)))
!
          hm(:)    = hhuth(:) * conv(:)
!
          hstar_hocl(:) = hm(:) * (1.0d0 + 1.052d0 *  &
     &                     exp(0.273d0 * (wtper(:) - 65.66d0)))
!
          hsqrtd(:) = hstar_hocl(:) * sqrt(d1)
!
          gcalc(:)  = 2.2548D-05 * hsqrtd(:) * sqrt(tk(:) * mw(IHOCL) * k(:))
!
!
!....   NOTE: gprob_tot is the overall uptake coeff for HOCl
!
          where( fterm > 0.0d0 )
            gprob_tot = 1.0d0 / (1.0d0 / (fterm(:) * gcalc(:)) + 1.0d0 / alpha)
          elsewhere
            gprob_tot = 0.0d0
          end where
!
          avgvel(:)   = 100.0d0 * (8.0d0 * 8.31448d0 * tk(:) * 1000.0d0 / (pi * mw(IHOCL)))**0.5d0
!
          where( hcl > minconc )
            sksts_hocl_hcl = 0.25d0 * gprob_tot * avgvel * sad / hcl
          elsewhere
            sksts_hocl_hcl = 0.25d0 * gprob_tot * avgvel * sad
          end where
!
          where( sad < 0.0d0 ) sksts_hocl_hcl   = 0.0d0
!
          if ( present(ptrop) ) then
            where( pr > ptrop ) sksts_hocl_hcl = 0.0d0
          end if
!
        END FUNCTION sksts_hocl_hcl
!
!.... sksts_hobr_hcl (temperature ,adcol ,pressure ,sad_lbs ,specarr(HOBr,:) ,specarr(HCl,:) ,water ,ptrop)
!
!_6_
!
!.... JPL 97-4
!
        FUNCTION sksts_hobr_hcl (tk,ad,pr,sad,hobr,hcl,h2o,ptrop)
          real*8, OPTIONAL :: ptrop
          real*8  tk(:) ,ad(:) ,pr(:) ,sad(:) ,h2o(:) ,hcl(:) ,hobr(:)
          real*8, DIMENSION (size(tk)) :: sksts_hobr_hcl
          real*8 adrop ,alpha ,d1 ,hsqrtd ,kii ,minconc ,pi,minadivl
          real*8 adivl(size(tk)) ,ah2o(size(tk)) ,avgvel(size(tk)) &
     &     ,fterm(size(tk)) &
     &     ,gcalc(size(tk)) &
     &     ,gprob_tot(size(tk)) &
     &     ,hstar(size(tk)) &
     &     ,k(size(tk)) &
     &     ,ph2o(size(tk)) ,phcl(size(tk)) &
     &     ,tk_150(size(tk))
!
          pi = acos(-1.0d0)
!
!=======================================================================
!     HOBr + HCl on stratospheric sulfate aerosol = BrCl + H2O
!=======================================================================
!
!.... First order reaction rate constant
!
!.... Hanson and Ravi, JPC, 98, 5728, 1994
!.... DEK, 1/10/97
!
!....   NOTE: alpha is modified from 0.3 in Table 2, JPC,
!....         Hanson and Ravi, 98, 5734, 1994 to 1.0 based on
!....         Ravi and Hanson, 101, pg 3887, JGR, 1996.
!
          minconc  = 1.0d0
          adrop    = 1.0D-05
          alpha    = 1.0d0
          d1       = 1.2D-08
          minadivl = 1.00D-15
!
!....    NOTE: Partial pressure of HCl and H2O (in atmospheres)
!
          ph2o(:) = (h2o(:) / ad(:)) * pr(:)
!
          ph2o(:) = ph2o(:) / 1013.25d0
!
          where( hcl(:) <= minconc )
            phcl   = ((1.0d0 / ad(:)) * pr(:)) * (1.0d0 / 1013.25d0)
          elsewhere
            phcl   = ((hcl(:) / ad(:)) * pr(:)) * (1.0d0 / 1013.25d0)
          end where
!
!....    NOTE: Activity of H2O not allowed to exceed 1.1
!....          ah2o and Hstar taken from Table 2, JPC,
!....          Hanson and Ravi, 98, 5734, 1994
!
          ah2o(:) = 1013.25d0 * ph2o(:) / 10.0d0**(9.217d0 - (2190.0d0 / (tk(:) - 12.7d0)))
!
          where( ah2o(:) > 1.1d0)
            ah2o    = 1.1d0
          end where
!
          tk_150(:) = tk(:)
!
          where( tk(:) < 150.0d0 )
            tk_150 = 150.0d0
          endwhere
!
          hstar(:) = exp((6250.0d0 / tk_150(:)) - 10.414d0) * (ah2o(:)**3.49d0)
!
          kii      = 1.0D+05
!
          k(:)     = kii * hstar(:) * phcl(:)
!
          hsqrtd   = 110.0d0
!
          gcalc(:) = 2.2548D-05 * hsqrtd * sqrt(tk(:) * mw(IHOBR) * k(:))
!
          adivl(:) = adrop / sqrt(d1 / k(:))
!
           where( adivl(:) < minadivl)
              adivl = minadivl
           endwhere
!
          fterm(:) = ((exp(adivl(:)) + exp(-adivl(:))) / &
     &                (exp(adivl(:)) - exp(-adivl(:)))) - &
     &               (1.0d0 / adivl(:))
!
!....   NOTE: gprob_tot is the overall uptake coeff for HOCl
!
          where( fterm > 0.0d0 )
            gprob_tot = 1.0d0 / (1.0d0 / (fterm(:) * gcalc(:)) + 1.0d0 / alpha)
          elsewhere
            gprob_tot = 0.0d0
          end where
!
          avgvel(:)   = 100.0d0 * (8.0d0 * 8.31448d0 * tk(:) * 1000.0d0 / &
     &                  (pi * mw(IHOBR)))**0.5d0
!
          where( hcl > minconc )
            sksts_hobr_hcl = 0.25d0 * gprob_tot * avgvel * sad / hcl
         elsewhere
            sksts_hobr_hcl = 0.25d0 * gprob_tot * avgvel * sad
         end where
!
          where( sad < 0.0d0 ) sksts_hobr_hcl = 0.0d0
!
          if ( present(ptrop) ) then
            where( pr > ptrop ) sksts_hobr_hcl = 0.0d0
          endif
!
!... JPL02 has lots of caveats and uncertainity - ignore for now
!!!!!!          sksts_hobr_hcl = 0.0d0
!
        END FUNCTION sksts_hobr_hcl
!
!.... sknat_clono2 (temperature ,pressure ,sad_nat ,ptrop)
!
!_1_
!
!.... (1) JPL 00-003
!
        FUNCTION sknat_clono2 (tk,pr,sad,ptrop)
          real*8, OPTIONAL :: ptrop
          real*8  tk(:) ,pr(:) ,sad(:)
          real*8, DIMENSION(size(tk)) :: sknat_clono2
          real*8  gprob ,pi
          real*8  avgvel(size(tk))
!
          pi = acos(-1.0d0)
!
!=======================================================================
!     ClONO2 + PSC Type I NAT particles = HOCl + HNO3
!=======================================================================
!
!.... First order reaction rate constant
!.... JPL 00-003
!.... PSC 1/16/2002
!
          gprob    = 0.004d0
          avgvel(:)= 100.0d0 * (8.0d0 * 8.31448d0 * tk(:) * 1000.0d0 /  &
     &               (pi * mw(ICLONO2)))**0.5d0
!
          sknat_clono2(:) = 0.25d0 * gprob * avgvel(:) * sad(:)
           if ( present(ptrop) ) then
            where( pr > ptrop ) sknat_clono2 = 0.0d0
          end if
        END FUNCTION sknat_clono2
!
!.... sknat_brono2 (temperature ,pressure ,sad_nat ,ptrop)
!
!_2_
!
!.... (2) JPL 00-003
!
        FUNCTION sknat_brono2 (tk,pr,sad,ptrop)
          real*8, OPTIONAL :: ptrop
          real*8  tk(:) ,pr(:) ,sad(:)
          real*8, DIMENSION(size(tk)) :: sknat_brono2
          real*8  gprob ,pi
          real*8  avgvel(size(tk))
!
          pi = acos(-1.0d0)
!
!=======================================================================
!     BrONO2 + PSC Type I NAT particles = HOBr + HNO3
!=======================================================================
!
!.... First order reaction rate constant
!.... JPL 00-003
!.... PSC 1/16/2002
!
          gprob     = 0.004d0
          avgvel(:) = 100.0d0 * (8.0d0 * 8.31448d0 * tk(:) * 1000.0d0 /  &
     &                (pi * mw(IBRONO2)))**0.5d0
!
          sknat_brono2(:) = 0.25d0 * gprob * avgvel(:) * sad(:)
          if ( present(ptrop) ) then
            where( pr > ptrop ) sknat_brono2 = 0.0d0
          end if
        END FUNCTION sknat_brono2
!
!.... sknat_hcl_clono2 (temperature ,pressure ,sad_nat ,specarr( HCl,:) ,ptrop)
!
!_3_
!
!.... (3) JPL 00-003
!
        FUNCTION sknat_hcl_clono2 (tk,pr,sad,hcl,ptrop)
          real*8, OPTIONAL :: ptrop
          real*8  tk(:) ,pr(:) ,sad(:) ,hcl(:)
          real*8, DIMENSION(size(tk)) :: sknat_hcl_clono2
          real*8  gprob ,minconc ,pi
          real*8  avgvel(size(tk))
!
          pi = acos(-1.0d0)
!
!=======================================================================
!     HCl + ClONO2 on PSC Type I NAT particles = Cl2 + HNO3
!=======================================================================
!
!.... Pseudo first order reaction rate constant
!.... JPL 00-003
!.... PSC 1/16/2002
!
          minconc  = 1.0d0
          gprob    = 0.20d0
          avgvel(:) = 100.0d0 *  &
     &               (8.0d0 * 8.31448d0 * tk(:) * 1000.0d0 /  &
     &               (pi * mw(ICLONO2)))**0.5d0
!
          sknat_hcl_clono2(:) = 0.25d0 * gprob * avgvel(:) * sad(:)
!
          where (hcl < minconc)
            sknat_hcl_clono2  = sknat_hcl_clono2 / minconc
          elsewhere
            sknat_hcl_clono2  = sknat_hcl_clono2 / hcl
          end where
          if ( present(ptrop) ) then
            where( pr > ptrop ) sknat_hcl_clono2 = 0.0d0
          end if
        END FUNCTION sknat_hcl_clono2
!
!.... sknat_hcl_hocl (temperature ,pressure ,sad_nat ,specarr( HCl,:) ,ptrop)
!
!_4_
!
!.... (4) JPL 00-003
!
        FUNCTION sknat_hcl_hocl (tk,pr,sad,hcl,ptrop)
          real*8, OPTIONAL :: ptrop
          real*8  tk(:) ,pr(:) ,sad(:) ,hcl(:)
          real*8, DIMENSION(size(tk)) :: sknat_hcl_hocl
          real*8  gprob ,minconc ,pi
          real*8  avgvel(size(tk))
!
          pi = acos(-1.0d0)
!
!=======================================================================
!     HCl + HOCl on PSC Type I NAT particles = Cl2 + H2O
!=======================================================================
!
!.... Pseudo first order reaction rate constant
!.... JPL 00-003
!.... PSC 1/16/2002
!
          minconc   = 1.0d0
          gprob     = 0.10d0
          avgvel(:) = 100.0d0 *  &
     &                (8.0d0 * 8.31448d0 * tk(:) * 1000.0d0 /  &
     &                (pi * mw(IHOCL)))**0.5d0
!
          sknat_hcl_hocl(:) = 0.25d0 * gprob * avgvel(:) * sad(:)
!
          where (hcl < minconc)
            sknat_hcl_hocl  = sknat_hcl_hocl / minconc
          elsewhere
            sknat_hcl_hocl  = sknat_hcl_hocl / hcl
          end where
          if ( present(ptrop) ) then
            where( pr > ptrop ) sknat_hcl_hocl = 0.0d0
          end if
        END FUNCTION sknat_hcl_hocl
!
!.... sknat_hcl_brono2 (temperature ,pressure ,sad_nat ,specarr( HCl,:) ,ptrop)
!
!_5_
!
!.... (5) JPL 00-003
!
        FUNCTION sknat_hcl_brono2 (tk,pr,sad,hcl,ptrop)
          real*8, OPTIONAL :: ptrop
          real*8  tk(:) ,pr(:) ,sad(:) ,hcl(:)
          real*8, DIMENSION(size(tk)) :: sknat_hcl_brono2
          real*8  gprob ,minconc ,pi
          real*8  avgvel(size(tk))
!
          pi = acos(-1.0d0)
!
!=======================================================================
!     HCl + BrONO2 on PSC Type I NAT particles = BrCl + HNO3
!=======================================================================
!
!.... Pseudo first order reaction rate constant
!.... JPL 00-003
!.... PSC 1/16/2002
!
          minconc   = 1.0d0
          gprob     = 0.20d0
          avgvel(:) = 100.0d0 *  &
     &                (8.0d0 * 8.31448d0 * tk(:) * 1000.0d0 /  &
     &                (pi * mw(IBRONO2)))**0.5d0
!
          sknat_hcl_brono2(:) = 0.25d0 * gprob * avgvel(:) * sad(:)
!
          where (hcl < minconc)
            sknat_hcl_brono2  = sknat_hcl_brono2 / minconc
          elsewhere
            sknat_hcl_brono2  = sknat_hcl_brono2 / hcl
          end where
          if ( present(ptrop) ) then
            where( pr > ptrop ) sknat_hcl_brono2 = 0.0d0
          end if
        END FUNCTION sknat_hcl_brono2
!
!.... sknat_hcl_hobr (temperature ,pressure ,sad_nat ,specarr( HCl,:) ,ptrop)
!
!_6_
!
!.... (6) JPL 00-003
!
        FUNCTION sknat_hcl_hobr (tk,pr,sad,hcl,ptrop)
          real*8, OPTIONAL :: ptrop
          real*8  tk(:) ,pr(:) ,sad(:) ,hcl(:)
          real*8, DIMENSION(size(tk)) :: sknat_hcl_hobr
          real*8  gprob ,minconc ,pi
          real*8  avgvel(size(tk))
!
          pi = acos(-1.0d0)
!
!=======================================================================
!     HCl + HOBr on PSC Type I NAT particles = BrCl + H2O
!=======================================================================
!
!.... Pseudo first order reaction rate constant
!.... JPL 00-003
!.... PSC 1/16/2002
!
          minconc   = 1.0d0
          gprob     = 0.10d0
          avgvel(:) = 100.0d0 *  &
     &                (8.0d0 * 8.31448d0 * tk(:) * 1000.0d0 /  &
     &                (pi * mw(IHOBR)))**0.5d0
!
          sknat_hcl_hobr(:) = 0.25d0 * gprob * avgvel(:) * sad(:)
!
          where (hcl < minconc)
            sknat_hcl_hobr  = sknat_hcl_hobr / minconc
          elsewhere
            sknat_hcl_hobr  = sknat_hcl_hobr / hcl
          end where
          if ( present(ptrop) ) then
            where( pr > ptrop ) sknat_hcl_hobr = 0.0d0
          end if
        END FUNCTION sknat_hcl_hobr
!
!.... skice_clono2 (temperature ,pressure ,sad_ice ,ptrop)
!
!_1_
!
!.... (1) JPL 00-003
!
        FUNCTION skice_clono2 (tk,pr,sad,ptrop)
          real*8, OPTIONAL :: ptrop
          real*8  tk(:) ,pr(:) ,sad(:)
          real*8, DIMENSION(size(tk)) :: skice_clono2
          real*8  gprob ,pi
          real*8  avgvel(size(tk))
!
          pi               = acos(-1.0d0)
!
!=======================================================================
!     ClONO2 + PSC Type II ice particles = HOCl + HNO3
!=======================================================================
!
!.... First order reaction rate constant
!.... JPL 00-003
!.... PSC 1/16/2002
!
          gprob            = 0.30d0
          avgvel(:)        = 100.0d0 *  &
     &                       (8.0d0 * 8.31448d0 * tk(:) * 1000.0d0 /  &
     &                       (pi * mw(ICLONO2)))**0.5d0
!
          skice_clono2(:)  = 0.25d0 * gprob * avgvel(:) * sad(:)
          where( sad < 0.0d0 ) skice_clono2   = 0.0d0
          if ( present(ptrop) ) then
            where( pr > ptrop ) skice_clono2 = 0.0d0
          end if
!
        END FUNCTION skice_clono2
!
!.... skice_brono2 (temperature ,pressure ,sad_ice ,ptrop)
!
!_2_
!
!.... (2) JPL 00-003
!
        FUNCTION skice_brono2 (tk,pr,sad,ptrop)
          real*8, OPTIONAL :: ptrop
          real*8  tk(:) ,pr(:) ,sad(:)
          real*8, DIMENSION(size(tk)) :: skice_brono2
          real*8  gprob ,pi
          real*8  avgvel(size(tk))
!
          pi = acos(-1.0d0)
!
!=======================================================================
!     BrONO2 + PSC Type II ice particles = HOBr + HNO3
!=======================================================================
!
!.... First order reaction rate constant
!.... JPL 00-003
!.... PSC 1/16/2002
!
          gprob            = 0.30d0
          avgvel(:)        = 100.0d0 * (8.0d0 * 8.31448d0 * tk(:) * 1000.0d0 /  &
     &                       (pi * mw(IBRONO2)))**0.5d0
!
          skice_brono2(:)  = 0.25d0 * gprob * avgvel(:) * sad(:)
          where( sad < 0.0d0 ) skice_brono2   = 0.0d0
          if ( present(ptrop) ) then
            where( pr > ptrop ) skice_brono2 = 0.0d0
          end if
!
        END FUNCTION skice_brono2
!
!.... skice_hcl_clono2 (temperature ,pressure ,sad_ice ,specarr( HCl,:) ,ptrop)
!
!_3_
!
!.... (3) JPL 00-003
!
        FUNCTION skice_hcl_clono2 (tk,pr,sad,hcl,ptrop)
          real*8, OPTIONAL :: ptrop
          real*8  tk(:) ,pr(:) ,sad(:) ,hcl(:)
          real*8, DIMENSION(size(tk)) :: skice_hcl_clono2
          real*8  gprob ,minconc ,pi
          real*8  avgvel(size(tk))
!
          pi = acos(-1.0d0)
!
!=======================================================================
!     HCl + ClONO2 on PSC Type II ice particles = Cl2 + HNO3
!=======================================================================
!
!.... Pseudo first order reaction rate constant
!.... JPL 00-003
!.... PSC 1/16/2002
!
          minconc              = 1.0d0
          gprob                = 0.30d0
          avgvel(:)            = 100.0d0 * (8.0d0 * 8.31448d0 * tk(:) *  &
     &                            1000.0d0 / (pi * mw(ICLONO2)))**0.5d0
!
          skice_hcl_clono2(:)  = 0.25d0 * gprob * avgvel(:) * sad(:)
!
          where (hcl < minconc)
            skice_hcl_clono2   = skice_hcl_clono2 / minconc
          elsewhere
            skice_hcl_clono2   = skice_hcl_clono2 / hcl
          end where
          where( sad < 0.0d0 ) skice_hcl_clono2   = 0.0d0
          if ( present(ptrop) ) then
            where( pr > ptrop ) skice_hcl_clono2 = 0.0d0
          end if
!
        END FUNCTION skice_hcl_clono2
!
!.... skice_hcl_hocl (temperature ,pressure ,sad_ice ,specarr( HCl,:) ,ptrop)
!
!_4_
!
!.... (4) JPL 00-003
!
        FUNCTION skice_hcl_hocl (tk,pr,sad,hcl,ptrop)
          real*8, OPTIONAL :: ptrop
          real*8  tk(:) ,pr(:) ,sad(:) ,hcl(:)
          real*8, DIMENSION(size(tk)) :: skice_hcl_hocl
          real*8  gprob ,minconc ,pi
          real*8  avgvel(size(tk))
!
          pi = acos(-1.0d0)
!
!=======================================================================
!     HCl + HOCl on PSC Type II ice particles = Cl2 + H2O
!=======================================================================
!
!.... Pseudo first order reaction rate constant
!.... JPL 00-003
!.... PSC 1/16/2002
!
          minconc            = 1.0d0
          gprob              = 0.20d0
          avgvel(:)          = 100.0d0 *  &
     &                         (8.0d0 * 8.31448d0 * tk(:) * 1000.0d0 /  &
     &                         (pi * mw(IHOCL)))**0.5d0
!
          skice_hcl_hocl(:)  = 0.25d0 * gprob * avgvel(:) * sad(:)
!
          where (hcl < minconc)
            skice_hcl_hocl   = skice_hcl_hocl / minconc
          elsewhere
            skice_hcl_hocl   = skice_hcl_hocl / hcl
          end where
          where( sad < 0.0d0 ) skice_hcl_hocl   = 0.0d0
          if ( present(ptrop) ) then
            where( pr > ptrop ) skice_hcl_hocl = 0.0d0
          end if
!
        END FUNCTION skice_hcl_hocl
!
!.... skice_hcl_brono2 (temperature ,pressure ,sad_ice ,specarr( HCl,:) ,ptrop)
!
!_5_
!
!.... (5) JPL 00-003
!
        FUNCTION skice_hcl_brono2 (tk,pr,sad,hcl,ptrop)
          real*8, OPTIONAL :: ptrop
          real*8  tk(:) ,pr(:) ,sad(:) ,hcl(:)
          real*8, DIMENSION(size(tk)) :: skice_hcl_brono2
          real*8  gprob ,minconc ,pi
          real*8  avgvel(size(tk))
!
          pi = acos(-1.0d0)
!
!=======================================================================
!     HCl + BrONO2 on PSC Type II ice particles = BrCl + HNO3
!=======================================================================
!
!.... Pseudo first order reaction rate constant
!.... JPL 00-003
!.... PSC 1/16/2002
!
          minconc              = 1.0d0
          gprob                = 0.30d0
          avgvel(:)            = 100.0d0 *  &
     &                           (8.0d0 * 8.31448d0 * tk(:) *  &
     &                            1000.0d0 /  &
     &                           (pi * mw(IBRONO2)))**0.5d0
!
          skice_hcl_brono2(:)  = 0.25d0 * gprob * avgvel(:) * sad(:)
!
          where (hcl < minconc)
            skice_hcl_brono2   = skice_hcl_brono2 / minconc
          elsewhere
            skice_hcl_brono2   = skice_hcl_brono2 / hcl
          end where
          where( sad < 0.0d0 ) skice_hcl_brono2   = 0.0d0
          if ( present(ptrop) ) then
            where( pr > ptrop ) skice_hcl_brono2 = 0.0d0
          end if
!
        END FUNCTION skice_hcl_brono2
!
!.... skice_hcl_hobr (temperature ,pressure ,sad_ice ,specarr( HCl,:) ,ptrop)
!
!_6_
!
!.... (6) JPL 00-003
!
        FUNCTION skice_hcl_hobr (tk,pr,sad,hcl,ptrop)
          real*8, OPTIONAL :: ptrop
          real*8  tk(:) ,pr(:) ,sad(:) ,hcl(:)
          real*8, DIMENSION(size(tk)) :: skice_hcl_hobr
          real*8  gprob ,minconc ,pi
          real*8  avgvel(size(tk))
!
          pi = acos(-1.0d0)
!
!=======================================================================
!     HCl + HOBr on PSC Type II ice particles = BrCl + H2O
!=======================================================================
!
!.... Pseudo first order reaction rate constant
!.... JPL 00-003
!.... PSC 1/16/2002
!
          minconc            = 1.0d0
          gprob              = 0.20d0
          avgvel(:)          = 100.0d0 * (8.0d0 * 8.31448d0 * tk(:) * 1000.0d0 /  &
     &                         (pi * mw(IHOBR)))**0.5d0
!
          skice_hcl_hobr(:)  = 0.25d0 * gprob * avgvel(:) * sad(:)
!
          where (hcl < minconc)
            skice_hcl_hobr   = skice_hcl_hobr / minconc
          elsewhere
            skice_hcl_hobr   = skice_hcl_hobr / hcl
          end where
          where( sad < 0.0d0 ) skice_hcl_hobr   = 0.0d0
          if ( present(ptrop) ) then
            where( pr > ptrop ) skice_hcl_hobr = 0.0d0
          end if
!
        END FUNCTION skice_hcl_hobr
!
!.... sksoot_hno3 (temperature ,sad_soot)
!
!.... (1)
!
        FUNCTION sksoot_hno3 (tk,sad)
          real*8  tk(:) ,sad(:)
          real*8, DIMENSION(size(tk)) :: sksoot_hno3
          real*8  gprob ,pi
          real*8  avgvel(size(tk))
!
          pi = acos(-1.0d0)
!
!=======================================================================
!     HNO3 + soot particles = OH + NO2
!=======================================================================
!
!.... First order reaction rate constant
!.... PSC 1/16/2002
!
          gprob           = 0.0d0
          avgvel(:)       = 100.0d0 * (8.0d0 * 8.31448d0 * tk(:) * 1000.0d0 /  &
     &                      (pi * mw(IHNO3)))**0.5d0
!
          sksoot_hno3(:) = 0.25d0 * gprob * avgvel(:) * sad(:)
!
        END FUNCTION sksoot_hno3
!
!.... sktrs_ho2 (temperature,sadcol2,adcol,radA,NSADaer,NSADdust,cPBLcol,pressure)
!
!_1_
!
!.... Harvard/GMI Tropospheric Chemistry
!
        FUNCTION sktrs_ho2 (tk, sad, ad, radA, NSADaer,NSADdust, continental_pbl, pr)
          real*8  tk(:) ,sad(:,:), ad(:), radA(:,:), pr(:)
          real*8, DIMENSION (size(tk)) :: sktrs_ho2, gamma
          real*8  pi, sqm
          real*8  avgvel(size(tk)), dfkg(size(tk))
          integer jj,k,ntotA,NSADaer,NSADdust
          integer, DIMENSION (size(tk)) :: continental_pbl
	
          ntotA=NSADaer+NSADdust
!
!=======================================================================
!     HO2 + tropospheric aerosol = 0.5 H2O2
!=======================================================================
! ntotA = # dust bins + # aerosol bins
! radA = radius of aerosol (cm)
! ad = molec/cm3 air
! tk = temperature (K)
! sad = surface area of aerosols/volume of air (cm2/cm3)
! conPBLFlag = 1 if in continental PBL, 0 if not.
!
! loss rate (k = 1/s) of species on aerosol surfaces
!
! k = sad * [ radA/Dg +4/(vL) ]^(-1)
!
! where
! Dg = gas phase diffusion coefficient (cm2/s)
! L = sticking coefficient (unitless)  = gamma
! v = mean molecular speed (cm/s) = [ 8RT / (pi*M) ]^1/2
!
! radA/Dg = uptake by gas-phase diffusion to the particle surface
! 4/(vL) = uptake by free molecular collisions of gas molecules with the surface
!=======================================================================
!
          sktrs_ho2(:)= 0.d0
          pi          = acos(-1.0d0)
          sqm         = SQRT(mw(IHO2))

! calculate gas phase diffusion coefficient (cm2/s)
          dfkg (:) = 9.45D17 / ad(:) * (tk(:))**0.5d0 * ( 3.472D-2  &
     &               + 1.D0/mw(IHO2) )**0.5d0
!
! calculate mean molecular speed (cm/s)
          avgvel(:) = 100.0d0 *  &
     &                (8.0d0 * 8.31448d0 * tk(:) * 1000.0d0 /  &
     &                (pi * mw(IHO2)))**0.5d0
!
! loop over dust and aerosols, summing up loss rate
      do jj = 1,ntotA

       gamma(:) = 0.00D+00
     
 ! obtain gamma when aerosol or dust is present

       DO k=1,size(tk)
	  IF(sad(jj,k) > 0.0d0) &
     &     gamma(k) = HO2(radA(jj,k), tk(k), ad(k), sqm,   &
     &  		  specarr(IHO2,k), jj, continental_pbl(k))
       END DO

! Bound gamma to avoid division by zero immediately below

       WHERE(gamma(:) < 1.00D-05) gamma(:) = 1.00D-05

       WHERE( sad(jj,:) > 0.0d0 )
         sktrs_ho2(:) = sktrs_ho2(:) +  &
     &         sad(jj,:) * ( 4.0d0 / ( gamma(:) * avgvel(:) )+  &
     &         radA(jj,:) / dfkg(:) )**(-1.0d0)
       END WHERE

      END DO
!
        END FUNCTION sktrs_ho2
!
!.... sktrs_no2 (temperature, sadcol2, adcol, radA, NSADaer,NSADdust,ptrop, pressure)
!
!_2_
!
!.... Harvard/GMI Tropospheric Chemistry
!
        FUNCTION sktrs_no2 (tk, sad, ad, radA, NSADaer,NSADdust, ptrop, pr)
          real*8, OPTIONAL :: ptrop
          real*8  tk(:) ,sad(:,:), ad(:), radA(:,:), pr(:)
          real*8, DIMENSION (size(tk)) :: sktrs_no2
          real*8  gamma ,pi
          real*8  avgvel(size(tk)), dfkg(size(tk))
          integer ntotA,jj,NSADaer,NSADdust
!
           ntotA=NSADaer+NSADdust
!
!=======================================================================
!     NO2 + tropospheric aerosol = 0.5 HNO3 + 0.5 HONO
!=======================================================================
! ntotA = # dust bins + # aerosol bins
! radA = radius of aerosol (cm)
! ad = molec/cm3 air
! tk = temperature (K)
! sad = surface area of aerosols/volume of air (cm2/cm3)
!
! loss rate (k = 1/s) of species on aerosol surfaces
!
! k = sad * [ radA/Dg +4/(vL) ]^(-1)
!
! where
! Dg = gas phase diffusion coefficient (cm2/s)
! L = sticking coefficient (unitless)  = gamma
! v = mean molecular speed (cm/s) = [ 8RT / (pi*M) ]^1/2
!
! radA/Dg = uptake by gas-phase diffusion to the particle surface
! 4/(vL) = uptake by free molecular collisions of gas molecules with the surface
!=======================================================================
!
          sktrs_no2(:)= 0.d0
!
          gamma       = 1.0d-04
          pi          = acos(-1.0d0)
!
! calculate gas phase diffusion coefficient (cm2/s)
          dfkg (:) = 9.45D17 / ad(:) * (tk(:))**0.5d0 * ( 3.472D-2  &
     &              + 1.D0/mw(INO2) )**0.5d0
!
! calculate mean molecular speed (cm/s)
          avgvel(:) = 100.0d0 *  &
     &               (8.0d0 * 8.31448d0 * tk(:) * 1000.0d0 /  &
     &               (pi * mw(INO2)))**0.5d0
!
! loop over dust and aerosols, summing up loss rate
      do jj = 1,ntotA
!
          where( sad(jj,:) > 0.0d0 )
            sktrs_no2(:) = sktrs_no2(:) +  &
     &            sad(jj,:) * ( 4.0d0 / ( gamma * avgvel(:) )+  &
     &            radA(jj,:) / dfkg(:) )**(-1.0d0)
                end where
!
       enddo
!... reaction only used in troposphere if ptrop is passed in
          if ( present(ptrop) ) then
            where( pr <= ptrop ) sktrs_no2 = 0.0d0
          end if
!
        END FUNCTION sktrs_no2
!
!.... sktrs_no3 (temperature, sadcol2, adcol, radA,NSADaer,NSADdust,ptrop, pressure)
!
!_3_
!
!.... Harvard/GMI Tropospheric Chemistry
!
        FUNCTION sktrs_no3 (tk, sad, ad, radA, NSADaer,NSADdust, ptrop, pr)
          real*8, OPTIONAL :: ptrop
          real*8  tk(:) ,sad(:,:), ad(:), radA(:,:), pr(:)
          real*8, DIMENSION (size(tk)) :: sktrs_no3
          real*8  gamma ,pi
          real*8  avgvel(size(tk)), dfkg(size(tk))
          integer NSADaer,NSADdust,ntotA,jj
!
          ntotA=NSADaer+NSADdust
!
!=======================================================================
!     NO3 + tropospheric aerosol = HNO3
!=======================================================================
! ntotA = # dust bins + # aerosol bins
! radA = radius of aerosol (cm)
! ad = molec/cm3 air
! tk = temperature (K)
! sad = surface area of aerosols/volume of air (cm2/cm3)
!
! loss rate (k = 1/s) of species on aerosol surfaces
!
! k = sad * [ radA/Dg +4/(vL) ]^(-1)
!
! where
! Dg = gas phase diffusion coefficient (cm2/s)
! L = sticking coefficient (unitless)  = gamma
! v = mean molecular speed (cm/s) = [ 8RT / (pi*M) ]^1/2
!
! radA/Dg = uptake by gas-phase diffusion to the particle surface
! 4/(vL) = uptake by free molecular collisions of gas molecules with the surface
!=======================================================================
!
          sktrs_no3(:)= 0.d0
!
          gamma       = 1.0d-03
          pi          = acos(-1.0d0)
!
! calculate gas phase diffusion coefficient (cm2/s)
          dfkg(:) = 9.45D17 / ad(:) * (tk(:))**0.5d0 * ( 3.472D-2  &
     &              + 1.D0/mw(INO3) )**0.5d0
!
! calculate mean molecular speed (cm/s)
          avgvel(:) = 100.0d0 *  &
     &               (8.0d0 * 8.31448d0 * tk(:) * 1000.0d0 /  &
     &               (pi * mw(INO3)))**0.5d0
!
! loop over dust and aerosols, summing up loss rate
      do jj = 1,ntotA
!
          where( sad(jj,:) > 0.0d0 )
            sktrs_no3(:) =  sktrs_no3(:) +  &
     &           sad(jj,:) * ( 4.0d0 / ( gamma * avgvel(:) )+  &
     &           radA(jj,:) / dfkg(:) )**(-1.0d0)
          endwhere
!
       enddo
!... reaction only used in troposphere if ptrop is passed in
          if ( present(ptrop) ) then
            where( pr <= ptrop ) sktrs_no3 = 0.0d0
          end if
!
        END FUNCTION sktrs_no3
!
!.... sktrs_n2o5 (temperature,sadcol2,adcol,radA,FRH,NSADaer,NSADdust, ptrop, pressure)
!
!_4_
!
!.... Harvard/GMI Tropospheric Chemistry
!
        FUNCTION sktrs_n2o5 (tk, sad, ad, radA,FRH, NSADaer,NSADdust,ptrop, pr)
          real*8, OPTIONAL :: ptrop
          real*8  tk(:) ,sad(:,:), ad(:), radA(:,:), FRH(:), pr(:)
          real*8, DIMENSION (size(tk)) :: sktrs_n2o5
          real*8  pi, FRH_P(size(tk)),ttk(size(tk)),fact(size(tk))
          real*8  avgvel(size(tk)),dfkg(size(tk)),gamma(size(tk))
          integer NSADaer,NSADdust,ntotA,jj
!
          ntotA=NSADaer+NSADdust
!=======================================================================
! N2O5 + tropospheric aerosol = 2 HNO3                  Branch 4
!=======================================================================
! FRH = relative humidity fraction (0-1)
! ntotA = # dust bins + # aerosol bins
! radA = radius of aerosol (cm)
! ad = molec/cm3 air
! tk = temperature (K)
! sad = surface area of aerosols/volume of air (cm2/cm3)
!
! loss rate (k = 1/s) of species on aerosol surfaces
!
! k = sad * [ radA/Dg +4/(vL) ]^(-1)
!
! where
! Dg = gas phase diffusion coefficient (cm2/s)
! L = sticking coefficient (unitless)  = gamma
! v = mean molecular speed (cm/s) = [ 8RT / (pi*M) ]^1/2
!
! radA/Dg = uptake by gas-phase diffusion to the particle surface
! 4/(vL) = uptake by free molecular collisions of gas molecules with the surface
!=======================================================================
!
          sktrs_n2o5(:)= 0.d0
!
          pi           = acos(-1.0d0)
!
! calculate gas phase diffusion coefficient (cm2/s)
          dfkg (:) =  9.45D17 / ad(:) * (tk(:))**0.5d0 * ( 3.472D-2  &
     &                + 1.D0/mw(IN2O5) )**0.5d0
!
! calculate mean molecular speed (cm/s)
          avgvel(:)    = 100.0d0 *  &
     &                   (8.0d0 * 8.31448d0 * tk(:) * 1000.0d0 /  &
     &                   (pi * mw(IN2O5)))**0.5d0
!
! loop over dust and aerosols, summing up loss rate
      do jj = 1,ntotA
!
!***********************************************************
! calculate gamma which is function of aerosol type, T, & RH
! following Evans and Jacob, "Impact of new laboratory studies of N2O5
! hydrolysis on global model budgets of tropospheric nitrogen oxides,
! ozone, and OH"
!***********************************************************
      ! Convert RH to % (max = 100%)
      FRH_P(:)  = FRH(:) * 100.d0
      where( FRH_P(:) > 100.d0) FRH_P(:) = 100.d0
!
        gamma(:) = 0d0
!
! DUST
!    Based on unpublished Crowley work
      if(jj.le.7) gamma(:) = 0.01d0
! SULFATE
      if(jj.eq.8) then
!===========================================================
! RH dependence from Kane et al., Heterogenous uptake of
! gaseous N2O5 by (NH4)2SO4, NH4HSO4 and H2SO4 aerosols
! J. Phys. Chem. A , 2001, 105, 6465-6470
!===========================================================
            gamma(:) = 2.79d-4 + FRH_P(:)*(  1.30d-4 +  &
     &                        FRH_P(:)*( -3.43d-6 +  &
     &                        FRH_P(:)*(  7.52d-8 ) ) )
!
!===========================================================
! Temperature dependence factor (Cox et al, Cambridge UK)
! is of the form:
!
!          10^( LOG10( G294 ) - 0.04 * ( TTEMP - 294 ) )
! FACT = -------------------------------------------------
!                     10^( LOG10( G294 ) )
!
! Where G294 = 1e-2 and TTEMP is MAX( TEMP, 282 ).
!
! For computational speed, replace LOG10( 1e-2 ) with -2
! and replace 10^( LOG10( G294 ) ) with G294
!===========================================================
            ttk(:) = tk(:)
            where( ttk(:) < 282d0) ttk(:) = 282d0
            fact(:) = 10d0**( -2d0 - 4d-2*(ttk(:) - 294.d0))/1d-2
!
            ! Apply temperature dependence
            gamma(:) = gamma(:) * fact(:)
       endif
! BLACK CARBON
!     from IUPAC
      if(jj.eq.9) gamma(:) = 0.005d0
! ORGANIC CARBON
      if(jj.eq.10) then
!===========================================================
! Based on Thornton, Braban and Abbatt, 2003
! N2O5 hydrolysis on sub-micron organic aerosol: the effect
! of relative humidity, particle phase and particle size
!===========================================================
            where ( FRH_P(:) >= 57d0 ) gamma(:) = 0.03d0
            where ( FRH_P(:) <  57d0 ) gamma(:) = FRH_P(:) * 5.2d-4
!
! Bryan & Jules 01/14/05
! Set gamma to very samll number when gamma=0, which occurs
! when relative humidity is 0 (e.g., UT/LS polar night).
!
            where ( gamma(:) == 0.0d0) gamma(:) = 0.03d-2
!
       endif
! SEA SALT
!     Based on IUPAC recomendation
      if(jj.ge.11) then
            where ( FRH_P(:) >= 62 ) gamma(:) = 0.03d0
            where ( FRH_P(:) <  62 ) gamma(:) = 0.005d0
      endif
!
! end calculation of gamma
!***********************************************************
! calculate loss rate
!***********************************************************
          where( sad(jj,:) > 0.0d0 )
            sktrs_n2o5(:) = sktrs_n2o5(:) +  &
     &            sad(jj,:) * ( 4.0d0 / ( gamma(:) * avgvel(:) )+  &
     &            radA(jj,:) / dfkg(:) )**(-1.0d0)
                end where
!
       enddo
!... reaction only used in troposphere if ptrop is passed in
          if ( present(ptrop) ) then
            where( pr <= ptrop ) sktrs_n2o5 = 0.0d0
          end if
!
        END FUNCTION sktrs_n2o5
  
      FUNCTION HO2( RADIUS, TEMP, DENAIR, SQM, HO2DENS, &
                    AEROTYPE, CONTINENTAL_PBL ) RESULT( GAMMA )

      !=================================================================
      ! Internal function HO2 computes the GAMMA reaction probability
      ! for HO2 loss in aerosols based on the recommendation of 
      ! Thornton, Jaegle, and McNeill, 
      ! "Assessing Known Pathways For HO2 Loss in Aqueous Atmospheric
      !  Aerosols: Regional and Global Impacts on Tropospheric Oxidants"
      !  J. Geophys. Res.,  doi:10.1029/2007JD009236, 2008  
      !
      ! gamma(HO2) is a function of aerosol type, radius, temperature
      !
      ! jaegle 01/22/2008
      ! 
      ! Arguments as Input:
      ! ----------------------------------------------------------------
      ! (1 ) RADIUS   (REAL*8 ) : Aerosol radius [cm]
      ! (2 ) TEMP     (REAL*8 ) : Temperature [K]
      ! (3 ) DENAIR   (REAL*8 ) : Air Density [molec/cm3]
      ! (4 ) HO2DENS  (REAL*8 ) : HO2 Number Density [molec/cm3]
      ! (5 ) SQM      (REAL*8 ) : Square root of molecular weight [g/mole]
      ! (6 ) AEROTYPE (INTEGER) : # denoting aerosol type (cf FAST-J)
      ! (7 ) CONTINENTAL_PBL (INTEGER)  : Flag set to 1 if the
      !         box is located in the continenal boundary layer,
      !         otherwise it is zero. Also check for ICE/SNOW (to
      !         disable this at high latitudes)
      !
      ! NOTES:
      !=================================================================
      

      ! Arguments
      REAL*8,  INTENT(IN) :: RADIUS, TEMP, DENAIR, HO2DENS, SQM
      INTEGER, INTENT(IN) :: AEROTYPE, CONTINENTAL_PBL

      ! Local variables
      REAL*8              :: ALPHA
      REAL*8              :: delG, Keq, w, H_eff
      REAL*8              :: A1, B1, k1, k2, A, B, C
      REAL*8              :: kaq, kmt, o2_ss, fluxrxn, DFKG
      REAL*8              :: TEST


      ! Avogadro's number
      REAL*8,  PARAMETER   :: Na = 6.022d23

      ! Ideal gas constant [atm cm3/mol/K], Raq
      REAL*8,  PARAMETER   :: Raq=82.d0

      ! Function return value
      REAL*8              :: GAMMA
      
      !=================================================================
      ! HO2 begins here!
      !=================================================================

      ! Default value
      GAMMA = 0.0d0

      ! Special handling for various aerosols
      SELECT CASE ( AEROTYPE )

         !----------------
         ! Dust 
         !----------------
         CASE ( 1, 2, 3, 4, 5, 6, 7 )      
                                
            ! Assume default gamma=0.1 on dust aerosols
            ! This is tentative as no lab measurements presently exist
            ! for gamma(HO2) on dust aerosols. We assume the rate to
            ! be fast on dust aerosols as transition metal ion induced
            ! chemistry is likely to occur in a thin aqueous surface layer.
            GAMMA = 0.1d0

         !----------------
         ! For Sulfate(8), Black Carbon (9), Organic Carbon (10),
         ! Sea-salt accum & coarse (11,12) calculate the 
         ! reaction probability due to self reaction 
         ! by using the algebraic expression in Thornton et al.  (2008)
         ! (equation 7) which is a function of temperature, aerosol radius,
         ! air density and HO2 concentration. 
         !
         ! Transition metal ions (such as copper and iron) in sea-salt and 
         ! carbonaceous aerosols are complexed to ligands and/or exist at 
         ! a concentration too low to catalyze HO2 loss efficiently, so we 
         ! apply the HO2 self reaction expression directly for these aerosols.
         ! 
         ! In the case of sulfate aerosol, the aerosols likely
         ! contain copper in the continental boundary layer and
         ! HO2 uptake proceeds rapidly. To account for the metal catalyzed
         ! uptake, we assume gamma(HO2)=0.07 (in the mid-range of the recommended
         ! 0.04-0.1 by Thornton et al, based on observed copper concentrations
         ! in the US boundary layer). Outside the continental boundary layer, we
         ! use the HO2-only algebraic expression.
         !
         !----------------
         CASE ( 8, 9, 10, 11, 12)  

            ! Mean molecular speed [cm/s]
            w = 14550.5d0 * sqrt(TEMP/(SQM*SQM))

            ! DFKG = Gas phase diffusion coeff [cm2/s]
            DFKG  = 9.45D17/DENAIR * SQRT(TEMP) * SQRT(3.472D-2 + 1.D0/(SQM*SQM))

            !calculate T-dependent solubility and aq. reaction rate constants
            ! hydronium ion concentration
            ! A1 = 1.+(Keq/hplus) 
            ! with Keq = 2.1d-5 [M], Equilibrium constant for 
            ! HO2aq = H+ + O2- (Jacob, 2000)
            !      hplus=10.d0^(-pH), with pH = 5
            ! B1 = Req * TEMP
            ! with Req = 1.987d-3 [kcal/K/mol], Ideal gas constant
            ! Note that we assume a constant pH of 5.
            A1 = 1.+ (2.1d-5 / (10.d0**(-5) ) )
            B1 = 1.987d-3 * TEMP

            ! Free energy change for solvation of HO2 (Hanson 1992, Golden 1991)
            ! in [kcal/mol]:
            ! delG = -4.9-(TEMP-298d0)*delS
            ! with delS=-0.023  [kcal/mol/K],  Entropy change for solvation of HO2
            delG  = -4.9d0 - (TEMP-298.d0) * (-0.023)
            H_eff = exp( -delG / B1 ) * A1

            ! Estimated temp dependent value for HO2 + O2- (k1) and 
            ! HO2+HO2 (see Jacob 1989)
            k1  =   1.58d10 * exp( -3. / B1 )
            k2  =   2.4d9   * exp( -4.7 / B1 )
            kaq = ( k1 * (A1 - 1.d0) + k2) / (A1**2)

            ! Calculate the mass transfer rate constant and s.s. conc. of 
            ! total HO2 in the aqueous phase:
            ! kmt = (RADIUS/DFKG + 4d0/w/alpha)^(-1)
            ! with alpha = mass accomodation coefficient, assumed 
            ! to be 1 (Thornton et al.)
            kmt = 1.d0/( RADIUS/DFKG + 4d0/w/1. )

            !use quadratic formula to obtain [O2-] in particle of radius RADIUS
            A = -2d0 * kaq
            B = -3d0 * kmt / RADIUS / (H_eff * 0.082 * TEMP)
            C =  3d0 * kmt * HO2DENS * 1000d0 / RADIUS / Na

            ! Error check that B^2-(4d0*A*C) is not negative
            TEST= B**2-(4d0*A*C)
            IF ( TEST < 0d0 ) THEN
                GAMMA = 0d0
            ELSE
                ! Calculate the concentration of O2- in the aerosol
                o2_ss= ( -B  -sqrt(B**2-(4d0*A*C)) )/(2d0*A)

                ! Calculate the reactive flux
                fluxrxn = kmt*HO2DENS - o2_ss*Na*kmt/H_eff/Raq/TEMP

                IF ( fluxrxn <= 0d0 ) THEN
                   GAMMA = 0d0
                ELSE
                   ! Gamma for HO2 at TEMP, ho2, and RADIUS given
                   GAMMA = 1./( ( ( HO2DENS/fluxrxn ) - ( RADIUS/DFKG ) ) * w / 4.d0 )
                ENDIF
            ENDIF
            ! For sulfate aerosols, check whether we are in
            ! the continental boundary layer, in which case
            ! copper catalyzed HO2 uptake likely dominates and
            ! speeds up the reaction: we assume gamma=0.07,
            ! which is in the middle of the 0.04-0.1 range recommended
            ! by Thornton et al. (2008)
            !
            IF ( AEROTYPE == 8 .and. CONTINENTAL_PBL == 1) THEN
                GAMMA = 0.07
            ENDIF 

         !----------------
         ! Default
         !----------------
         CASE DEFAULT
            WRITE (6,*) 'Not a suitable aerosol surface '
            WRITE (6,*) 'for HO2 uptake'
            WRITE (6,*) 'AEROSOL TYPE =',AEROTYPE
            STOP

      END SELECT
     
      ! If negative value is calculated, set it to zero
      IF ( GAMMA  <= 0d0 ) GAMMA = 0d0

      ! Return to CALCRATE
      END FUNCTION HO2
      END SUBROUTINE kcalc
