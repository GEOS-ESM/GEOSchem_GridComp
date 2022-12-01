!=======================================================================
!
! $Id: $
!
! FILE
!   setkin_lchem.h - character labels for species and reactions
!             (setkin_lchem.h)
!   12 JUN 02 - PSC
!
! DESCRIPTION
!   Include file that provides ascii strings identifying reactions
!   and species
!
!  Input mechanism:        GeosCCM_Combo_2015mechanism.txt
!  Reaction dictionary:    GMI_reactions_JPL15.db
!  Setkin files generated: Tue Mar 30 16:58:02 2021
!
!=======================================================================


      integer kmg_i

      character*16 lchemvar(NSP)
      character*180 lqkchem(NUM_K),  lqjchem(NUM_J)
!
!.... All species labels
!
      data lchemvar(1) /"A3O2"/
      data lchemvar(2) /"ACTA"/
      data lchemvar(3) /"ALD2"/
      data lchemvar(4) /"ALK4"/
      data lchemvar(5) /"ATO2"/
      data lchemvar(6) /"B3O2"/
      data lchemvar(7) /"Br"/
      data lchemvar(8) /"BrCl"/
      data lchemvar(9) /"BrO"/
      data lchemvar(10) /"BrONO2"/
      data lchemvar(11) /"C2H6"/
      data lchemvar(12) /"C3H8"/
      data lchemvar(13) /"CCl4"/
      data lchemvar(14) /"CF2Br2"/
      data lchemvar(15) /"CFC12"/
      data lchemvar(16) /"CF2ClBr"/
      data lchemvar(17) /"CF3Br"/
      data lchemvar(18) /"CFC113"/
      data lchemvar(19) /"CFC114"/
      data lchemvar(20) /"CFC115"/
      data lchemvar(21) /"CFC11"/
      data lchemvar(22) /"CH2Br2"/
      data lchemvar(23) /"CH2O"/
      data lchemvar(24) /"CH3Br"/
      data lchemvar(25) /"CH3CCl3"/
      data lchemvar(26) /"CH3Cl"/
      data lchemvar(27) /"CH4"/
      data lchemvar(28) /"CHBr3"/
      data lchemvar(29) /"Cl2"/
      data lchemvar(30) /"Cl2O2"/
      data lchemvar(31) /"Cl"/
      data lchemvar(32) /"ClO"/
      data lchemvar(33) /"ClONO2"/
      data lchemvar(34) /"CO"/
      data lchemvar(35) /"DMS"/
      data lchemvar(36) /"EOH"/
      data lchemvar(37) /"ETO2"/
      data lchemvar(38) /"ETP"/
      data lchemvar(39) /"GLYC"/
      data lchemvar(40) /"GLYX"/
      data lchemvar(41) /"H2402"/
      data lchemvar(42) /"H2"/
      data lchemvar(43) /"H2O2"/
      data lchemvar(44) /"H2O"/
      data lchemvar(45) /"HAC"/
      data lchemvar(46) /"HBr"/
      data lchemvar(47) /"HCFC141b"/
      data lchemvar(48) /"HCFC142b"/
      data lchemvar(49) /"HCFC22"/
      data lchemvar(50) /"HCl"/
      data lchemvar(51) /"HCOOH"/
      data lchemvar(52) /"HFC125"/
      data lchemvar(53) /"HFC134a"/
      data lchemvar(54) /"HFC143a"/
      data lchemvar(55) /"HFC152a"/
      data lchemvar(56) /"HFC23"/
      data lchemvar(57) /"HFC32"/
      data lchemvar(58) /"H"/
      data lchemvar(59) /"HNO2"/
      data lchemvar(60) /"HNO3"/
      data lchemvar(61) /"HNO4"/
      data lchemvar(62) /"HO2"/
      data lchemvar(63) /"HOBr"/
      data lchemvar(64) /"HOCl"/
      data lchemvar(65) /"IALD"/
      data lchemvar(66) /"IAO2"/
      data lchemvar(67) /"IAP"/
      data lchemvar(68) /"INO2"/
      data lchemvar(69) /"INPN"/
      data lchemvar(70) /"ISN1"/
      data lchemvar(71) /"ISNP"/
      data lchemvar(72) /"ISOP"/
      data lchemvar(73) /"KO2"/
      data lchemvar(74) /"MACR"/
      data lchemvar(75) /"MAN2"/
      data lchemvar(76) /"MAO3"/
      data lchemvar(77) /"MAOP"/
      data lchemvar(78) /"MAP"/
      data lchemvar(79) /"MCO3"/
      data lchemvar(80) /"MEK"/
      data lchemvar(81) /"MGLY"/
      data lchemvar(82) /"MO2"/
      data lchemvar(83) /"MOH"/
      data lchemvar(84) /"MP"/
      data lchemvar(85) /"MRO2"/
      data lchemvar(86) /"MRP"/
      data lchemvar(87) /"MVK"/
      data lchemvar(88) /"N2O5"/
      data lchemvar(89) /"N2O"/
      data lchemvar(90) /"N"/
      data lchemvar(91) /"NO2"/
      data lchemvar(92) /"NO3"/
      data lchemvar(93) /"NO"/
      data lchemvar(94) /"O1D"/
      data lchemvar(95) /"O3"/
      data lchemvar(96) /"OClO"/
      data lchemvar(97) /"OCSg"/
      data lchemvar(98) /"OH"/
      data lchemvar(99) /"O"/
      data lchemvar(100) /"PAN"/
      data lchemvar(101) /"PMN"/
      data lchemvar(102) /"PO2"/
      data lchemvar(103) /"PP"/
      data lchemvar(104) /"PPN"/
      data lchemvar(105) /"PRN1"/
      data lchemvar(106) /"PRPE"/
      data lchemvar(107) /"PRPN"/
      data lchemvar(108) /"R4N1"/
      data lchemvar(109) /"R4N2"/
      data lchemvar(110) /"R4O2"/
      data lchemvar(111) /"R4P"/
      data lchemvar(112) /"RA3P"/
      data lchemvar(113) /"RB3P"/
      data lchemvar(114) /"RCHO"/
      data lchemvar(115) /"RCO3"/
      data lchemvar(116) /"RCOOH"/
      data lchemvar(117) /"RIO1"/
      data lchemvar(118) /"RIO2"/
      data lchemvar(119) /"RIP"/
      data lchemvar(120) /"ROH"/
      data lchemvar(121) /"RP"/
      data lchemvar(122) /"VRO2"/
      data lchemvar(123) /"VRP"/
      data lchemvar(124) /"SO2"/
      data lchemvar(125) /"H2SO4"/
      data lchemvar(126) /"ACET"/
      data lchemvar(127) /"N2"/
      data lchemvar(128) /"O2"/
      data lchemvar(129) /"NUMDENS"/
      data lchemvar(130) /"HNO3COND"/
!
!.... Thermal reaction labels
!
      data (lqkchem(kmg_i), kmg_i=1,10) / &
     & 'O + O2 = O3', &
     & 'O + O3 = 2 O2', &
     & 'N2 + O1D = N2 + O', &
     & 'O1D + O2 = O + O2', &
     & 'O1D + O3 = 2 O2', &
     & 'O1D + O3 = 2 O + O2', &
     & 'H2O + O1D = 2 OH', &
     & 'H2 + O1D = H + OH', &
     & 'N2O + O1D = N2 + O2', &
     & 'N2O + O1D = 2 NO' /

      data (lqkchem(kmg_i), kmg_i=11,20) / &
     & 'CH4 + O1D = MO2 + OH', &
     & 'CH4 + O1D = CH2O + H + HO2', &
     & 'CH4 + O1D = CH2O + H2', &
     & 'CFC12 + O1D = 2 Cl', &
     & 'CFC113 + O1D = 3 Cl', &
     & 'CFC114 + O1D = 2 Cl', &
     & 'CFC115 + O1D = Cl', &
     & 'HCFC22 + O1D = Cl', &
     & 'HCFC141b + O1D = 2 Cl', &
     & 'HCFC142b + O1D = Cl' /

      data (lqkchem(kmg_i), kmg_i=21,30) / &
     & 'H + O2 = HO2', &
     & 'H + O3 = O2 + OH', &
     & 'O + OH = H + O2', &
     & 'HO2 + O = O2 + OH', &
     & 'H + HO2 = 2 OH', &
     & 'NO + O3 = NO2 + O2', &
     & 'O3 + OH = HO2 + O2', &
     & 'HO2 + O3 = 2 O2 + OH', &
     & 'NO2 + O3 = NO3 + O2', &
     & 'OH + OH = H2O + O' /

      data (lqkchem(kmg_i), kmg_i=31,40) / &
     & 'OH + OH = H2O2', &
     & 'HO2 + OH = H2O + O2', &
     & 'H2O2 + OH = H2O + HO2', &
     & 'HO2 + NO = NO2 + OH', &
     & 'HO2 + HO2 = H2O2 + O2', &
     & 'H2O + HO2 + HO2 = H2O + H2O2 + O2', &
     & 'H2 + OH = H + H2O', &
     & 'CO + OH = H', &
     & 'CH4 + OH = H2O + MO2', &
     & 'MO2 + NO = CH2O + HO2 + NO2' /

      data (lqkchem(kmg_i), kmg_i=41,50) / &
     & 'ClO + MO2 = CH2O + Cl + HO2 + O2', &
     & 'HO2 + MO2 = MP + O2', &
     & 'MO2 + MO2 = CH2O + MOH + O2', &
     & 'MO2 + MO2 = 2 CH2O + 2 HO2', &
     & 'MP + OH = H2O + MO2', &
     & 'MP + OH = CH2O + H2O + OH', &
     & 'CH2O + OH = CO + H2O + HO2', &
     & 'N + O2 = NO + O', &
     & 'N + NO = N2 + O', &
     & 'N + NO2 = N2O + O' /

      data (lqkchem(kmg_i), kmg_i=51,60) / &
     & 'NO2 + O = NO + O2', &
     & 'NO3 + O = NO2 + O2', &
     & 'NO2 + OH = HNO3', &
     & 'HNO3 + OH = H2O + NO3', &
     & 'NO + OH = HNO2', &
     & 'HNO2 + OH = H2O + NO2', &
     & 'HO2 + NO2 = HNO4', &
     & 'HNO4 = HO2 + NO2', &
     & 'HNO4 + OH = H2O + NO2 + O2', &
     & 'HO2 + NO3 = NO2 + O2 + OH' /

      data (lqkchem(kmg_i), kmg_i=61,70) / &
     & 'NO + NO3 = 2 NO2', &
     & 'NO3 + OH = HO2 + NO2', &
     & 'NO2 + NO3 = N2O5', &
     & 'N2O5 = NO2 + NO3', &
     & 'HCOOH + OH = H2O + HO2', &
     & 'MOH + OH = CH2O + HO2', &
     & 'NO2 + NO3 = NO + NO2 + O2', &
     & 'CH2O + NO3 = CO + HNO3 + HO2', &
     & 'Cl + O3 = ClO + O2', &
     & 'Cl + H2 = H + HCl' /

      data (lqkchem(kmg_i), kmg_i=71,80) / &
     & 'Cl + H2O2 = HCl + HO2', &
     & 'Cl + HO2 = HCl + O2', &
     & 'Cl + HO2 = ClO + OH', &
     & 'ClO + O = Cl + O2', &
     & 'ClO + OH = Cl + HO2', &
     & 'ClO + OH = HCl + O2', &
     & 'ClO + HO2 = HOCl + O2', &
     & 'ClO + NO = Cl + NO2', &
     & 'ClO + NO2 = ClONO2', &
     & 'ClO + ClO = 2 Cl + O2' /

      data (lqkchem(kmg_i), kmg_i=81,90) / &
     & 'ClO + ClO = Cl2 + O2', &
     & 'ClO + ClO = Cl + OClO', &
     & 'ClO + ClO = Cl2O2', &
     & 'Cl2O2 = 2 ClO', &
     & 'HCl + OH = Cl + H2O', &
     & 'HOCl + OH = ClO + H2O', &
     & 'ClONO2 + O = ClO + NO3', &
     & 'ClONO2 + OH = HOCl + NO3', &
     & 'Cl + ClONO2 = Cl2 + NO3', &
     & 'Br + O3 = BrO + O2' /

      data (lqkchem(kmg_i), kmg_i=91,100) / &
     & 'Br + HO2 = HBr + O2', &
     & 'Br + CH2O = CO + HBr + HO2', &
     & 'BrO + O = Br + O2', &
     & 'BrO + HO2 = HOBr + O2', &
     & 'BrO + NO = Br + NO2', &
     & 'BrO + NO2 = BrONO2', &
     & 'BrO + ClO = Br + OClO', &
     & 'BrO + ClO = Br + Cl + O2', &
     & 'BrO + ClO = BrCl + O2', &
     & 'BrO + BrO = 2 Br + O2' /

      data (lqkchem(kmg_i), kmg_i=101,110) / &
     & 'HBr + OH = Br + H2O', &
     & 'CHBr3 + OH = 3 Br', &
     & 'CH2Br2 + OH = 2 Br', &
     & 'CH2O + O = CO + HO2 + OH', &
     & 'CH4 + Cl = HCl + MO2', &
     & 'CH2O + Cl = CO + HCl + HO2', &
     & 'CH3Cl + OH = Cl + H2O + HO2', &
     & 'CH3CCl3 + OH = 3 Cl + H2O', &
     & 'HCFC22 + OH = Cl + H2O', &
     & 'HCFC141b + OH = 2 Cl + H2O' /

      data (lqkchem(kmg_i), kmg_i=111,120) / &
     & 'HCFC142b + OH = Cl + H2O', &
     & 'CH3Cl + Cl = CO + 2 HCl + HO2', &
     & 'CH3Br + OH = Br + H2O + HO2', &
     & 'HFC23 + O1D =  0.25 H2O +  0.75 O', &
     & 'HFC32 + O1D =  0.30 H2O +  0.70 O', &
     & 'HFC125 + O1D =  0.15 H2O +  0.25 O +  0.60 OH', &
     & 'HFC134a + O1D =  0.11 H2O +  0.65 O +  0.24 OH', &
     & 'HFC143a + O1D =  0.27 H2O +  0.35 O +  0.38 OH', &
     & 'HFC152a + O1D =  0.40 H2O +  0.45 O +  0.15 OH', &
     & 'HFC23 + OH = H2O' /

      data (lqkchem(kmg_i), kmg_i=121,130) / &
     & 'HFC32 + OH = H2O', &
     & 'HFC125 + OH = H2O', &
     & 'HFC134a + OH = H2O', &
     & 'HFC143a + OH = H2O', &
     & 'HFC152a + OH = H2O', &
     & 'A3O2 + HO2 = RA3P', &
     & 'A3O2 + MO2 =  0.75 CH2O + HO2 +  0.25 MOH +  0.75 RCHO +  0.25 ROH', &
     & 'A3O2 + NO = HO2 + NO2 + RCHO', &
     & 'ACET + OH = ATO2 + H2O', &
     & 'ACTA + OH = H2O + MO2' /

      data (lqkchem(kmg_i), kmg_i=131,140) / &
     & 'ALD2 + NO3 = HNO3 + MCO3', &
     & 'ALD2 + OH =  0.05 CH2O +  0.05 CO + H2O +  0.05 HO2 +  0.95 MCO3', &
     & 'ALK4 + NO3 = HNO3 + R4O2', &
     & 'ALK4 + OH = R4O2', &
     & 'ATO2 + HO2 = MCO3 + MO2', &
     & 'ATO2 + MCO3 = ACTA + MEK', &
     & 'ATO2 + MCO3 =  0.20 CH2O +  0.80 HO2 +  0.20 MCO3 +  0.80 MGLY + MO2', &
     & 'ATO2 + MO2 =  0.50 CH2O +  0.20 HAC +  0.30 HO2 +  0.30 MCO3 +  0.50 MGLY +  0.50 MOH', &
     & 'ATO2 + NO =  0.96 CH2O +  0.96 MCO3 +  0.96 NO2 +  0.04 R4N2', &
     & 'B3O2 + HO2 = RB3P' /

      data (lqkchem(kmg_i), kmg_i=141,150) / &
     & 'B3O2 + MCO3 = ACET + ACTA', &
     & 'B3O2 + MCO3 = ACET + HO2 + MO2', &
     & 'B3O2 + MO2 =  0.75 ACET +  0.75 CH2O + HO2 +  0.25 MOH +  0.25 ROH', &
     & 'B3O2 + NO = ACET + HO2 + NO2', &
     & 'C2H6 + NO3 = ETO2 + HNO3', &
     & 'C2H6 + OH = ETO2 + H2O', &
     & 'C2H6 + Cl = ETO2 + HCl', &
     & 'C3H8 + OH = A3O2', &
     & 'C3H8 + OH = B3O2', &
     & 'EOH + OH = ALD2 + HO2' /

      data (lqkchem(kmg_i), kmg_i=151,160) / &
     & 'ETO2 + ETO2 = 2 ALD2 + 2 HO2', &
     & 'ETO2 + ETO2 = ALD2 + EOH', &
     & 'ETO2 + NO = ALD2 + HO2 + NO2', &
     & 'ETP + OH =  0.50 ALD2 +  0.50 ETO2 +  0.50 OH', &
     & 'GLYC + OH =  0.73 CH2O +  0.50 CO +  0.13 GLYX +  0.13 HCOOH +  0.77 HO2 +  0.23 OH', &
     & 'GLYC + OH = CO + HCOOH + OH', &
     & 'GLYX + NO3 = 2 CO + HNO3 + HO2', &
     & 'GLYX + OH = 2 CO + HO2', &
     & 'HAC + OH = HO2 + MGLY', &
     & 'HAC + OH =  0.50 ACTA +  0.50 CO +  0.50 HCOOH +  0.50 MO2 + OH' /

      data (lqkchem(kmg_i), kmg_i=161,170) / &
     & 'ETO2 + HO2 = ETP', &
     & 'HO2 + MCO3 = ACTA + O3', &
     & 'HO2 + MCO3 = MAP', &
     & 'IALD + O3 =  0.12 CH2O +  0.28 GLYC +  0.20 GLYX +  0.20 HAC +  0.20 HCOOH +  0.60 MGLY +  0.30 O3 +  0.10 OH', &
     & 'IALD + OH =  0.15 HO2 +  0.44 IAO2 +  0.41 MAO3', &
     & 'HO2 + IAO2 = IAP', &
     & 'IAO2 + MCO3 = ACTA + MEK', &
     & 'IAO2 + MCO3 =  0.40 CH2O +  0.29 CO +  0.26 GLYC +  0.18 GLYX +  0.36 HAC + HO2 +  0.58 MGLY + MO2', &
     & 'IAO2 + MO2 =  0.95 CH2O +  0.15 CO +  0.13 GLYC +  0.09 GLYX +  0.18 HAC + HO2 +  0.25 MEK +  0.29 MGLY +  0.25 MOH +  0.25 ROH', &
     & 'IAO2 + NO =  0.35 CH2O +  0.27 CO +  0.24 GLYC +  0.17 GLYX +  0.33 HAC +  0.08 HNO3 +  0.92 HO2 +  0.53 MGLY +  0.92 NO2' /

      data (lqkchem(kmg_i), kmg_i=171,180) / &
     & 'IAP + OH =  0.50 IAO2 +  0.50 OH +  0.50 RCHO', &
     & 'HO2 + INO2 = INPN', &
     & 'INO2 + MCO3 =  0.15 CH2O +  0.85 HNO3 +  0.80 HO2 +  0.10 MACR + MO2 +  0.05 MVK +  0.15 NO2', &
     & 'INO2 + MCO3 = ACTA + NO2 + RCHO', &
     & 'INO2 + MO2 =  0.83 CH2O +  0.43 HNO3 +  0.90 HO2 +  0.05 MACR +  0.25 MOH +  0.03 MVK +  0.57 NO2 +  0.25 RCHO +  0.25 ROH', &
     & 'INO2 + NO =  0.15 CH2O +  0.85 HNO3 +  0.80 HO2 +  0.10 MACR +  0.05 MVK +  1.15 NO2', &
     & 'INPN + OH = INO2', &
     & 'HO2 + ISN1 = ISNP', &
     & 'ISN1 + MCO3 = GLYC + HAC + MO2 + NO2', &
     & 'ISN1 + MCO3 = ACTA + NO2 + RCHO' /

      data (lqkchem(kmg_i), kmg_i=181,190) / &
     & 'ISN1 + MO2 =  0.75 CH2O +  0.50 GLYC +  0.50 HAC +  0.50 HO2 +  0.25 MOH + NO2 +  0.25 RCHO +  0.25 ROH', &
     & 'ISNP + OH =  0.50 ISN1 +  0.50 NO2 +  0.50 OH +  0.50 RCHO', &
     & 'ISOP + NO3 = INO2', &
     & 'ISOP + O3 =  0.90 CH2O +  0.05 CO +  0.06 HO2 +  0.39 MACR +  0.16 MVK +  0.10 O3 +  0.27 OH +  0.07 PRPE', &
     & 'ISOP + OH = RIO2', &
     & 'HO2 + KO2 = MGLY + MO2', &
     & 'KO2 + MCO3 = ACTA + MEK', &
     & 'KO2 + MCO3 = ALD2 + MCO3 + MO2', &
     & 'KO2 + MO2 =  0.50 ALD2 +  0.75 CH2O +  0.50 HO2 +  0.50 MCO3 +  0.25 MEK +  0.25 MOH +  0.25 ROH', &
     & 'KO2 + NO =  0.93 ALD2 +  0.93 MCO3 +  0.93 NO2 +  0.07 R4N2' /

      data (lqkchem(kmg_i), kmg_i=191,200) / &
     & 'MACR + NO3 = MAN2', &
     & 'MACR + NO3 = HNO3 + MAO3', &
     & 'MACR + O3 =  0.70 CH2O +  0.20 CO +  0.28 HO2 +  0.80 MGLY +  0.20 O3 +  0.22 OH', &
     & 'MACR + OH =  0.53 MAO3 +  0.47 MRO2', &
     & 'HO2 + MAN2 = ISNP', &
     & 'MAN2 + MCO3 = CH2O + MGLY + MO2 + NO2', &
     & 'MAN2 + MCO3 = ACTA + NO2 + RCHO', &
     & 'MAN2 + MO2 =  1.25 CH2O +  0.50 HO2 +  0.50 MGLY +  0.25 MOH + NO2 +  0.25 RCHO +  0.25 ROH', &
     & 'MAN2 + NO = CH2O + MGLY + 2 NO2', &
     & 'HO2 + MAO3 =  0.59 CH2O +  0.39 CO +  0.41 MAOP +  0.39 MO2 +  0.15 O3 +  0.44 OH' /

      data (lqkchem(kmg_i), kmg_i=201,210) / &
     & 'MAO3 + MCO3 = CH2O + MCO3 + MO2', &
     & 'MAO3 + MO2 = 2 CH2O + HO2 + MCO3', &
     & 'MAO3 + MO2 = CH2O + RCOOH', &
     & 'MAO3 + NO2 = PMN', &
     & 'MAO3 + NO =  0.50 CH2O +  0.50 CO +  0.50 MCO3 +  0.50 MO2 + NO2', &
     & 'MAOP + OH = MAO3', &
     & 'A3O2 + MCO3 = ACTA + RCHO', &
     & 'A3O2 + MCO3 = HO2 + MO2 + RCHO', &
     & 'ETO2 + MCO3 = ACTA + ALD2', &
     & 'ETO2 + MCO3 = ALD2 + HO2 + MO2' /

      data (lqkchem(kmg_i), kmg_i=211,220) / &
     & 'MCO3 + MCO3 = 2 MO2', &
     & 'MCO3 + MO2 = ACTA + CH2O', &
     & 'MCO3 + MO2 = CH2O + HO2 + MO2', &
     & 'MCO3 + NO2 = PAN', &
     & 'MCO3 + NO = MO2 + NO2', &
     & 'MCO3 + PO2 = ACTA +  0.65 HAC +  0.35 RCHO', &
     & 'MCO3 + PO2 = ALD2 + CH2O + HO2 + MO2', &
     & 'MEK + NO3 = HNO3 + KO2', &
     & 'MEK + OH = H2O + KO2', &
     & 'MGLY + NO3 = CO + HNO3 + MCO3' /

      data (lqkchem(kmg_i), kmg_i=221,230) / &
     & 'MGLY + OH = CO + MCO3', &
     & 'ETO2 + MO2 =  0.75 ALD2 +  0.75 CH2O +  0.25 EOH + HO2 +  0.25 MOH', &
     & 'HO2 + MRO2 = MRP', &
     & 'MCO3 + MRO2 = ACTA + MEK', &
     & 'MCO3 + MRO2 =  0.17 CH2O +  0.83 CO +  0.83 HAC + HO2 +  0.17 MGLY + MO2', &
     & 'MO2 + MRO2 = CH2O +  0.60 CO +  0.60 HAC + HO2 +  0.25 MGLY +  0.15 ROH', &
     & 'MRO2 + NO = HNO3', &
     & 'MRO2 + NO =  0.17 CH2O +  0.83 CO +  0.83 HAC + HO2 +  0.17 MGLY + NO2', &
     & 'MRP + OH = MRO2', &
     & 'MVK + O3 =  0.04 ALD2 +  0.80 CH2O +  0.05 CO +  0.06 HO2 +  0.82 MGLY +  0.20 O3 +  0.08 OH' /

      data (lqkchem(kmg_i), kmg_i=231,240) / &
     & 'MVK + OH = VRO2', &
     & 'MAP + OH =  0.50 CH2O +  0.50 MCO3 +  0.50 OH', &
     & 'OH + RCHO = H2O + RCO3', &
     & 'OH + RCOOH = ETO2', &
     & 'PAN = MCO3 + NO2', &
     & 'PMN = MAO3 + NO2', &
     & 'O3 + PMN =  0.60 CH2O + HO2 + NO2', &
     & 'OH + PMN = CO + HAC + NO2', &
     & 'HO2 + PO2 = PP', &
     & 'MO2 + PO2 =  0.50 ALD2 + CH2O +  0.16 HAC + HO2 +  0.25 MOH +  0.09 RCHO +  0.25 ROH' /

      data (lqkchem(kmg_i), kmg_i=241,250) / &
     & 'NO + PO2 = ALD2 + CH2O + HO2 + NO2', &
     & 'PPN = NO2 + RCO3', &
     & 'OH + PP =  0.79 HAC +  0.79 OH +  0.21 PO2', &
     & 'HO2 + PRN1 = PRPN', &
     & 'MCO3 + PRN1 = ALD2 + CH2O + MO2 + NO2', &
     & 'MCO3 + PRN1 = ACTA + NO2 + RCHO', &
     & 'MO2 + PRN1 =  0.50 ALD2 +  1.25 CH2O +  0.50 HO2 +  0.25 MOH + NO2 +  0.25 RCHO +  0.25 ROH', &
     & 'NO + PRN1 = ALD2 + CH2O + 2 NO2', &
     & 'NO3 + PRPE = PRN1', &
     & 'O3 + PRPE =  0.50 ALD2 +  0.54 CH2O +  0.42 CO +  0.06 H2 +  0.30 HO2 +  0.31 MO2 +  0.14 OH' /

      data (lqkchem(kmg_i), kmg_i=251,260) / &
     & 'OH + PRPE = PO2', &
     & 'OH + PRPN = PRN1', &
     & 'HO2 + R4N1 = R4N2', &
     & 'MCO3 + R4N1 =  0.75 ALD2 +  0.39 CH2O + MO2 + NO2 +  0.30 R4O2 +  0.57 RCHO', &
     & 'MCO3 + R4N1 = ACTA + NO2 + RCHO', &
     & 'MO2 + R4N1 =  0.38 ALD2 +  0.95 CH2O +  0.50 HO2 +  0.25 MOH + NO2 +  0.15 R4O2 +  0.54 RCHO +  0.25 ROH', &
     & 'NO + R4N1 =  0.75 ALD2 +  0.39 CH2O + 2 NO2 +  0.30 R4O2 +  0.57 RCHO', &
     & 'OH + R4N2 = H2O + R4N1', &
     & 'HO2 + R4O2 = R4P', &
     & 'MCO3 + R4O2 = ACTA + MEK' /

      data (lqkchem(kmg_i), kmg_i=261,270) / &
     & 'MCO3 + R4O2 =  0.05 A3O2 +  0.32 ACET +  0.32 ALD2 +  0.18 B3O2 +  0.32 ETO2 +  0.27 HO2 +  0.19 MEK +  1.18 MO2 +  0.13 RCHO', &
     & 'MO2 + R4O2 =  0.03 A3O2 +  0.16 ACET +  0.16 ALD2 +  0.09 B3O2 +  0.75 CH2O +  0.16 ETO2 +  0.64 HO2 +  0.35 MEK +  0.09 MO2 +  0.25 MOH +  0.07 RCHO +  0.25 ROH', &
     & 'NO + R4O2 =  0.05 A3O2 +  0.32 ACET +  0.32 ALD2 +  0.18 B3O2 +  0.32 ETO2 +  0.27 HO2 +  0.19 MEK +  0.18 MO2 + NO2 +  0.13 RCHO', &
     & 'NO + R4O2 = R4N2', &
     & 'OH + R4P =  0.50 OH +  0.50 R4O2 +  0.50 RCHO', &
     & 'OH + RA3P =  0.50 A3O2 +  0.50 OH +  0.50 RCHO', &
     & 'OH + RB3P =  0.79 ACET +  0.21 B3O2 +  0.79 OH', &
     & 'NO3 + RCHO = HNO3 + RCO3', &
     & 'HO2 + RCO3 =  0.30 O3 +  0.30 RCOOH +  0.70 RP', &
     & 'MCO3 + RCO3 = ETO2 + MO2' /

      data (lqkchem(kmg_i), kmg_i=271,280) / &
     & 'MO2 + RCO3 = CH2O + ETO2 + HO2', &
     & 'MO2 + RCO3 = CH2O + RCOOH', &
     & 'NO2 + RCO3 = PPN', &
     & 'NO + RCO3 = ETO2 + NO2', &
     & 'HO2 + RIO1 = RIP', &
     & 'MCO3 + RIO1 = ACTA + MEK', &
     & 'MCO3 + RIO1 =  0.75 CH2O + HO2 + IALD + MO2', &
     & 'MO2 + RIO1 =  1.13 CH2O + HO2 +  0.50 IALD +  0.25 MEK +  0.25 MOH +  0.25 ROH', &
     & 'NO + RIO1 = HNO3', &
     & 'NO + RIO1 =  0.75 CH2O + HO2 + IALD + NO2' /

      data (lqkchem(kmg_i), kmg_i=281,290) / &
     & 'HO2 + RIO2 = RIP', &
     & 'MCO3 + RIO2 = ACTA + MEK', &
     & 'MCO3 + RIO2 =  0.69 CH2O +  0.86 HO2 +  0.13 IALD +  0.29 MACR + MO2 +  0.40 MVK +  0.14 RIO1', &
     & 'MO2 + RIO2 =  1.10 CH2O +  0.93 HO2 +  0.06 IALD +  0.14 MACR +  0.25 MEK +  0.25 MOH +  0.20 MVK +  0.07 RIO1 +  0.25 ROH', &
     & 'NO + RIO2 = HNO3', &
     & 'NO + RIO2 =  0.69 CH2O +  0.86 HO2 +  0.13 IALD +  0.29 MACR +  0.40 MVK + NO2 +  0.14 RIO1', &
     & 'OH + RIP =  0.50 IAO2 +  0.10 RIO1 +  0.40 RIO2', &
     & 'OH + ROH = HO2 + RCHO', &
     & 'OH + RP =  0.50 ALD2 +  0.50 OH +  0.50 RCO3', &
     & 'HO2 + VRO2 =  0.10 CH2O +  0.58 GLYC +  0.10 HO2 +  0.58 MCO3 +  0.19 MEK +  0.10 MGLY +  0.68 OH +  0.03 RCHO +  0.10 VRP' /

      data (lqkchem(kmg_i), kmg_i=291,300) / &
     & 'MCO3 + VRO2 = ACTA + MEK', &
     & 'MCO3 + VRO2 =  0.28 CH2O +  0.72 GLYC +  0.28 HO2 +  0.72 MCO3 +  0.28 MGLY + MO2', &
     & 'MO2 + VRO2 =  0.89 CH2O +  0.36 GLYC +  0.64 HO2 +  0.36 MCO3 +  0.25 MEK +  0.14 MGLY +  0.25 MOH +  0.25 ROH', &
     & 'NO + VRO2 = HNO3', &
     & 'NO + VRO2 =  0.28 CH2O +  0.72 GLYC +  0.28 HO2 +  0.72 MCO3 +  0.28 MGLY + NO2', &
     & 'OH + VRP =  0.50 OH +  0.50 RCHO +  0.50 VRO2', &
     & 'N2O5 = 2 HNO3', &
     & 'ClONO2 = HNO3 + HOCl', &
     & 'BrONO2 = HNO3 + HOBr', &
     & 'ClONO2 + HCl = Cl2 + HNO3' /

      data (lqkchem(kmg_i), kmg_i=301,310) / &
     & 'HCl + HOCl = Cl2 + H2O', &
     & 'HCl + HOBr = BrCl + H2O', &
     & 'N2O5 = 2 HNO3', &
     & 'ClONO2 = HNO3 + HOCl', &
     & 'BrONO2 = HNO3 + HOBr', &
     & 'ClONO2 + HCl = Cl2 + HNO3', &
     & 'HCl + HOCl = Cl2 + H2O', &
     & 'HCl + HOBr = BrCl + H2O', &
     & 'ClONO2 = HNO3 + HOCl', &
     & 'BrONO2 = HNO3 + HOBr' /

      data (lqkchem(kmg_i), kmg_i=311,320) / &
     & 'ClONO2 + HCl = Cl2 + HNO3', &
     & 'HCl + HOCl = Cl2 + H2O', &
     & 'BrONO2 + HCl = BrCl + HNO3', &
     & 'HCl + HOBr = BrCl + H2O', &
     & 'ClONO2 = HNO3 + HOCl', &
     & 'BrONO2 = HNO3 + HOBr', &
     & 'ClONO2 + HCl = Cl2 + HNO3', &
     & 'HCl + HOCl = Cl2 + H2O', &
     & 'BrONO2 + HCl = BrCl + HNO3', &
     & 'HCl + HOBr = BrCl + H2O' /

      data (lqkchem(kmg_i), kmg_i=321,330) / &
     & 'HNO3 = NO2 + OH', &
     & 'NO3 + NO3 = 2 NO2 + O2', &
     & 'HO2 =  0.50 H2O', &
     & 'NO2 =  0.50 HNO2 +  0.50 HNO3', &
     & 'NO3 = HNO3', &
     & 'N2O5 = 2 HNO3', &
     & 'DMS + OH = O2 + SO2', &
     & 'DMS + NO3 = HNO3 + SO2', &
     & 'O + SO2 = H2SO4', &
     & 'OH + SO2 = H2SO4' /

      data (lqkchem(kmg_i), kmg_i=331,333) / &
     & 'O3 + SO2 = H2SO4', &
     & 'O + OCSg = CO + SO2', &
     & 'OCSg + OH = SO2' /
!
!.... Photolytic reaction labels
!
      data (lqjchem(kmg_i), kmg_i=1,10) / &
     & 'O2 + hv = 2 O', &
     & 'O3 + hv = O + O2', &
     & 'O3 + hv = O1D + O2', &
     & 'NO + hv = N + O', &
     & 'N2O + hv = N2 + O1D', &
     & 'NO2 + hv = NO + O', &
     & 'H2O2 + hv = 2 OH', &
     & 'MP + hv = CH2O + HO2 + OH', &
     & 'CH2O + hv = CO + H + HO2', &
     & 'CH2O + hv = CO + H2' /

      data (lqjchem(kmg_i), kmg_i=11,20) / &
     & 'HNO3 + hv = NO2 + OH', &
     & 'HNO2 + hv = NO + OH', &
     & 'HNO4 + hv = NO3 + OH', &
     & 'HNO4 + hv = HO2 + NO2', &
     & 'NO3 + hv = NO2 + O3', &
     & 'NO3 + hv = NO + O2', &
     & 'N2O5 + hv = NO2 + NO3', &
     & 'Cl2 + hv = 2 Cl', &
     & 'OClO + hv = ClO + O', &
     & 'Cl2O2 + hv = 2 Cl + O2' /

      data (lqjchem(kmg_i), kmg_i=21,30) / &
     & 'HOCl + hv = Cl + OH', &
     & 'ClONO2 + hv = Cl + NO3', &
     & 'ClONO2 + hv = ClO + NO2', &
     & 'BrCl + hv = Br + Cl', &
     & 'BrO + hv = Br + O', &
     & 'HOBr + hv = Br + OH', &
     & 'BrONO2 + hv = Br + NO3', &
     & 'CHBr3 + hv = 3 Br', &
     & 'CH2Br2 + hv = 2 Br', &
     & 'CH3Cl + hv = Cl + MO2' /

      data (lqjchem(kmg_i), kmg_i=31,40) / &
     & 'CCl4 + hv = 4 Cl', &
     & 'CH3CCl3 + hv = 3 Cl', &
     & 'CFC11 + hv = 3 Cl', &
     & 'CFC12 + hv = 2 Cl', &
     & 'CFC113 + hv = 3 Cl', &
     & 'CFC114 + hv = 2 Cl', &
     & 'CFC115 + hv = Cl', &
     & 'HCFC141b + hv = 2 Cl', &
     & 'HCFC142b + hv = Cl', &
     & 'CH3Br + hv = Br + MO2' /

      data (lqjchem(kmg_i), kmg_i=41,50) / &
     & 'CF3Br + hv = Br', &
     & 'CF2Br2 + hv = 2 Br', &
     & 'H2402 + hv = 2 Br', &
     & 'CF2ClBr + hv = Br + Cl', &
     & 'ALD2 + hv = CO + HO2 + MO2', &
     & 'PAN + hv = MCO3 + NO2', &
     & 'RCHO + hv = CO + ETO2 + HO2', &
     & 'ACET + hv = MCO3 + MO2', &
     & 'MEK + hv = ETO2 + MCO3', &
     & 'GLYC + hv = CH2O + CO + 2 HO2' /

      data (lqjchem(kmg_i), kmg_i=51,60) / &
     & 'GLYX + hv = 2 CO + H2', &
     & 'GLYX + hv = 2 CO + 2 HO2', &
     & 'GLYX + hv = CH2O + CO', &
     & 'MGLY + hv = CO + HO2 + MCO3', &
     & 'MVK + hv = CO + PRPE', &
     & 'MVK + hv = MAO3 + MO2', &
     & 'MACR + hv =  0.20 CH2O + CO +  1.80 HO2 +  0.20 MCO3 +  0.80 MGLY', &
     & 'HAC + hv = CH2O + HO2 + MCO3', &
     & 'INPN + hv = HO2 + NO2 + OH + RCHO', &
     & 'PRPN + hv = HO2 + NO2 + OH + RCHO' /

      data (lqjchem(kmg_i), kmg_i=61,70) / &
     & 'ETP + hv = ALD2 + HO2 + OH', &
     & 'RA3P + hv = HO2 + OH + RCHO', &
     & 'RB3P + hv = HO2 + OH + RCHO', &
     & 'R4P + hv = HO2 + OH + RCHO', &
     & 'PP + hv = HO2 + OH + RCHO', &
     & 'RP + hv = ALD2 + HO2 + OH', &
     & 'RIP + hv =  0.69 CH2O +  0.86 HO2 +  0.13 IALD +  0.29 MACR +  0.40 MVK + OH +  0.14 RIO1', &
     & 'IAP + hv =  0.67 CO +  0.26 GLYC +  0.19 H2 +  0.36 HAC + HO2 +  0.58 MGLY + OH', &
     & 'ISNP + hv = HO2 + NO2 + OH + RCHO', &
     & 'VRP + hv =  0.28 CH2O +  0.72 GLYC +  0.28 HO2 +  0.72 MCO3 +  0.28 MGLY + OH' /

      data (lqjchem(kmg_i), kmg_i=71,76) / &
     & 'MRP + hv =  0.17 CH2O +  0.83 CO +  0.83 HAC + HO2 +  0.17 MGLY + OH', &
     & 'MAOP + hv = HO2 + OH + RCHO', &
     & 'R4N2 + hv =  0.05 A3O2 +  0.32 ACET +  0.32 ALD2 +  0.18 B3O2 +  0.32 ETO2 +  0.27 HO2 +  0.19 MEK +  0.18 MO2 + NO2 +  0.13 RCHO', &
     & 'MAP + hv = MO2 + OH', &
     & 'OCSg + hv = CO + SO2', &
     & 'H2SO4 + hv = H2O + SO2' /

!                                  --^--

