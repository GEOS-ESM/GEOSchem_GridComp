!=======================================================================
!
! $Id: setkin_fastj.h,v 1.1 2019/12/17 17:10:31 mmanyin Exp $
!
! FILE
!   setkin_fastj.h - lookup matrix for GMI vs FastJX cross-sections
!             (setkin_fastj.h)
!   Oct 2018 - SDS
!
! DESCRIPTION
!   Include file that provides ascii strings identifying
!   GMI vs FastJX photo reaction names
!
!  Input mechanism file:    GMI_Combo_2015mechanism.txt
!  Reaction dictionary:     GMI_reactions_JPL15.db
!  Setkin files generated:  Tue Nov 20 18:21:26 2018
!
!=======================================================================
!
!.... 
      integer idxph
!
!... Reference label for FastJX
      character*6 :: fastj_lookup(NUM_J)
      data (fastj_lookup(idxph), idxph=1,NUM_J) / &
     & "O2    ", "O3    ", "O3(1D)", "NO    ", "N2O   ",  &
     & "NO2   ", "H2O2  ", "CH3OOH", "H2COa ", "H2COb ",  &
     & "HNO3  ", "HNO2  ", "HNO4  ", "NO3   ", "NO3   ",  &
     & "N2O5  ", "HNO4  ", "Cl2   ", "OClO  ", "Cl2O2 ",  &
     & "HOCl  ", "ClNO3a", "ClNO3b", "BrCl  ", "BrO   ",  &
     & "HOBr  ", "BrNO3 ", "CH3Cl ", "CCl4  ", "MeCCl3",  &
     & "CFC11 ", "CFC12 ", "F113  ", "F114  ", "F115  ",  &
     & "F141b ", "F142b ", "CH3Br ", "H1301 ", "H1211 ",  &
     & "H2402 ", "H1211 ", "ActAld", "PAN   ", "PrAld ",  &
     & "Acet-a", "MEKeto", "GlyAld", "Glyxlb", "Glyxla",  &
     & "Glyxlb", "MGlyxl", "MeVK  ", "MeVK  ", "MeAcr ",  &
     & "Acet-a", "CH3OOH", "CH3OOH", "CH3OOH", "CH3OOH",  &
     & "CH3OOH", "CH3OOH", "CH3OOH", "CH3OOH", "CH3OOH",  &
     & "CH3OOH", "CH3OOH", "CH3OOH", "CH3OOH", "CH3OOH",  &
     & "CH3NO3", "CH3OOH", "OCS   " /
!... Quantum yield for photolysis
      real*8 :: fastj_yield(NUM_J)
      data (fastj_yield(idxph), idxph=1,NUM_J) / &
     &   1.000D+00,   1.000D+00,   1.000D+00,   1.000D+00,   1.000D+00,   &
     &   1.000D+00,   1.000D+00,   1.000D+00,   1.000D+00,   1.000D+00,   &
     &   1.000D+00,   1.000D+00,   3.330D-01,   8.860D-01,   1.140D-01,   &
     &   1.000D+00,   6.670D-01,   1.000D+00,   1.000D+00,   1.000D+00,   &
     &   1.000D+00,   1.000D+00,   1.000D+00,   1.000D+00,   1.000D+00,   &
     &   1.000D+00,   1.000D+00,   1.000D+00,   1.000D+00,   1.000D+00,   &
     &   1.000D+00,   1.000D+00,   1.000D+00,   1.000D+00,   1.000D+00,   &
     &   1.000D+00,   1.000D+00,   1.000D+00,   1.000D+00,   1.000D+00,   &
     &   1.000D+00,   1.000D+00,   1.000D+00,   1.000D+00,   1.000D+00,   &
     &   1.000D+00,   1.000D+00,   1.000D+00,   1.000D+00,   1.000D+00,   &
     &   1.000D+00,   1.000D+00,   6.000D-01,   4.000D-01,   1.000D+00,   &
     &   1.000D+00,   1.000D+00,   1.000D+00,   1.000D+00,   1.000D+00,   &
     &   1.000D+00,   1.000D+00,   1.000D+00,   1.000D+00,   1.000D+00,   &
     &   1.000D+00,   1.000D+00,   1.000D+00,   1.000D+00,   1.000D+00,   &
     &   1.000D+00,   1.000D+00,   1.000D+00  /
!
!... end setkin_fastj.h
