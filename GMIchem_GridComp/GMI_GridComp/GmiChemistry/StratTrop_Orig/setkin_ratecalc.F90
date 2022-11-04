!=======================================================================
!
! $Id$
!
! ROUTINE
!   Calc_rate_Setkin - GMI version (setkin_ratecalc.F)
!   13 JUN 02 - PSC
!
! DESCRIPTION
!   Calculates and returns rates of kinetic processes for block
!   of boxes
!
! ARGUMENTS
!   INPUT
!    nblock  : number of blocks
!    numj    : number of photolysis rate constants
!    numk    : number of thermal rate constants
!    numspc  ; number of species
!    qk     : thermal reaction rate constants (cm^3 molecule^-1 s^-1)
!    qj     : photolysis rate constants (s^-1)
!    y      : concentrations of transported species (molecules cm^-3)
!   OUTPUT
!    qqk    : rates of thermal processes (molecules cm^-3 s^-1)
!    qqj    : rates of photolytic processes (molecules cm^-3 s^-1)
!
!  Chemistry input file:    08/2018
!  Reaction dictionary:     GMI_Combo_rxns_119species_SO2_JPL15_OCS.db
!  Setkin files generated:  Tue Sep 18 18:25:41 2018
!
!=======================================================================
      subroutine Calc_rate_Setkin &
     &  (nblock, numj, numk, numspc, &
     &   qk, qj, y, qqk, qqj)

      implicit none

      integer &
     &  nblock, &
     &  numj, numk, numspc

      real*8 &
     &  qj    (nblock, numj), &
     &  qk    (nblock, numk), &
     &  qqj   (nblock, numj), &
     &  qqk   (nblock, numk), &
     &  y     (nblock, numspc)

      integer &
     &  kloop

      do kloop = 1, nblock
!
!                  start thermal reactions
!
!....         O + O2 = O3
!
      qqk(kloop,1)=qk(kloop,1)*y(kloop,22)*y(kloop,117)
!
!....         O + O3 = 2 O2
!
      qqk(kloop,2)=qk(kloop,2)*y(kloop,22)*y(kloop,24)
!
!....         N2 + O1D = N2 + O
!
      qqk(kloop,3)=qk(kloop,3)*y(kloop,116)*y(kloop,23)
!
!....         O1D + O2 = O + O2
!
      qqk(kloop,4)=qk(kloop,4)*y(kloop,23)*y(kloop,117)
!
!....         O1D + O3 = 2 O2
!
      qqk(kloop,5)=qk(kloop,5)*y(kloop,23)*y(kloop,24)
!
!....         O1D + O3 = 2 O + O2
!
      qqk(kloop,6)=qk(kloop,6)*y(kloop,23)*y(kloop,24)
!
!....         H2O + O1D = 2 OH
!
      qqk(kloop,7)=qk(kloop,7)*y(kloop,10)*y(kloop,23)
!
!....         H2 + O1D = H + OH
!
      qqk(kloop,8)=qk(kloop,8)*y(kloop,5)*y(kloop,23)
!
!....         N2O + O1D = N2 + O2
!
      qqk(kloop,9)=qk(kloop,9)*y(kloop,17)*y(kloop,23)
!
!....         N2O + O1D = 2 NO
!
      qqk(kloop,10)=qk(kloop,10)*y(kloop,17)*y(kloop,23)
!
!....         CH4 + O1D = MO2 + OH
!
      qqk(kloop,11)=qk(kloop,11)*y(kloop,2)*y(kloop,23)
!
!....         CH4 + O1D = CH2O + H + HO2
!
      qqk(kloop,12)=qk(kloop,12)*y(kloop,2)*y(kloop,23)
!
!....         CH4 + O1D = CH2O + H2
!
      qqk(kloop,13)=qk(kloop,13)*y(kloop,2)*y(kloop,23)
!
!....         CFC12 + O1D = 2 Cl
!
      qqk(kloop,14)=qk(kloop,14)*y(kloop,45)*y(kloop,23)
!
!....         CFC113 + O1D = 3 Cl
!
      qqk(kloop,15)=qk(kloop,15)*y(kloop,46)*y(kloop,23)
!
!....         CFC114 + O1D = 2 Cl
!
      qqk(kloop,16)=qk(kloop,16)*y(kloop,47)*y(kloop,23)
!
!....         CFC115 + O1D = Cl
!
      qqk(kloop,17)=qk(kloop,17)*y(kloop,48)*y(kloop,23)
!
!....         HCFC22 + O1D = Cl
!
      qqk(kloop,18)=qk(kloop,18)*y(kloop,49)*y(kloop,23)
!
!....         HCFC141b + O1D = 2 Cl
!
      qqk(kloop,19)=qk(kloop,19)*y(kloop,50)*y(kloop,23)
!
!....         HCFC142b + O1D = Cl
!
      qqk(kloop,20)=qk(kloop,20)*y(kloop,51)*y(kloop,23)
!
!....         H + O2 = HO2
!
      qqk(kloop,21)=qk(kloop,21)*y(kloop,4)*y(kloop,117)
!
!....         H + O3 = O2 + OH
!
      qqk(kloop,22)=qk(kloop,22)*y(kloop,4)*y(kloop,24)
!
!....         O + OH = H + O2
!
      qqk(kloop,23)=qk(kloop,23)*y(kloop,22)*y(kloop,25)
!
!....         HO2 + O = O2 + OH
!
      qqk(kloop,24)=qk(kloop,24)*y(kloop,11)*y(kloop,22)
!
!....         H + HO2 = 2 OH
!
      qqk(kloop,25)=qk(kloop,25)*y(kloop,4)*y(kloop,11)
!
!....         NO + O3 = NO2 + O2
!
      qqk(kloop,26)=qk(kloop,26)*y(kloop,18)*y(kloop,24)
!
!....         O3 + OH = HO2 + O2
!
      qqk(kloop,27)=qk(kloop,27)*y(kloop,24)*y(kloop,25)
!
!....         HO2 + O3 = 2 O2 + OH
!
      qqk(kloop,28)=qk(kloop,28)*y(kloop,11)*y(kloop,24)
!
!....         NO2 + O3 = NO3 + O2
!
      qqk(kloop,29)=qk(kloop,29)*y(kloop,19)*y(kloop,24)
!
!....         OH + OH = H2O + O
!
      qqk(kloop,30)=qk(kloop,30)*y(kloop,25)*y(kloop,25)
!
!....         OH + OH = H2O2
!
      qqk(kloop,31)=qk(kloop,31)*y(kloop,25)*y(kloop,25)
!
!....         HO2 + OH = H2O + O2
!
      qqk(kloop,32)=qk(kloop,32)*y(kloop,11)*y(kloop,25)
!
!....         H2O2 + OH = H2O + HO2
!
      qqk(kloop,33)=qk(kloop,33)*y(kloop,12)*y(kloop,25)
!
!....         HO2 + NO = NO2 + OH
!
      qqk(kloop,34)=qk(kloop,34)*y(kloop,11)*y(kloop,18)
!
!....         HO2 + HO2 = H2O2 + O2
!
      qqk(kloop,35)=qk(kloop,35)*y(kloop,11)*y(kloop,11)
!
!....         H2O + HO2 + HO2 = H2O + H2O2 + O2
!
      qqk(kloop,36)=qk(kloop,36)*y(kloop,10)*y(kloop,11)*y(kloop,11)
!
!....         H2 + OH = H + H2O
!
      qqk(kloop,37)=qk(kloop,37)*y(kloop,5)*y(kloop,25)
!
!....         CO + OH = H
!
      qqk(kloop,38)=qk(kloop,38)*y(kloop,3)*y(kloop,25)
!
!....         CH4 + OH = H2O + MO2
!
      qqk(kloop,39)=qk(kloop,39)*y(kloop,2)*y(kloop,25)
!
!....         MO2 + NO = CH2O + HO2 + NO2
!
      qqk(kloop,40)=qk(kloop,40)*y(kloop,13)*y(kloop,18)
!
!....         ClO + MO2 = CH2O + Cl + HO2 + O2
!
      qqk(kloop,41)=qk(kloop,41)*y(kloop,34)*y(kloop,13)
!
!....         HO2 + MO2 = MP + O2
!
      qqk(kloop,42)=qk(kloop,42)*y(kloop,11)*y(kloop,13)
!
!....         MO2 + MO2 = CH2O + MOH + O2
!
      qqk(kloop,43)=qk(kloop,43)*y(kloop,13)*y(kloop,13)
!
!....         MO2 + MO2 = 2 CH2O + 2 HO2
!
      qqk(kloop,44)=qk(kloop,44)*y(kloop,13)*y(kloop,13)
!
!....         MP + OH = H2O + MO2
!
      qqk(kloop,45)=qk(kloop,45)*y(kloop,15)*y(kloop,25)
!
!....         MP + OH = CH2O + H2O + OH
!
      qqk(kloop,46)=qk(kloop,46)*y(kloop,15)*y(kloop,25)
!
!....         CH2O + OH = CO + H2O + HO2
!
      qqk(kloop,47)=qk(kloop,47)*y(kloop,1)*y(kloop,25)
!
!....         N + O2 = NO + O
!
      qqk(kloop,48)=qk(kloop,48)*y(kloop,16)*y(kloop,117)
!
!....         N + NO = N2 + O
!
      qqk(kloop,49)=qk(kloop,49)*y(kloop,16)*y(kloop,18)
!
!....         NO2 + O = NO + O2
!
      qqk(kloop,50)=qk(kloop,50)*y(kloop,19)*y(kloop,22)
!
!....         NO3 + O = NO2 + O2
!
      qqk(kloop,51)=qk(kloop,51)*y(kloop,20)*y(kloop,22)
!
!....         NO2 + OH = HNO3
!
      qqk(kloop,52)=qk(kloop,52)*y(kloop,19)*y(kloop,25)
!
!....         HNO3 + OH = H2O + NO3
!
      qqk(kloop,53)=qk(kloop,53)*y(kloop,8)*y(kloop,25)
!
!....         NO + OH = HNO2
!
      qqk(kloop,54)=qk(kloop,54)*y(kloop,18)*y(kloop,25)
!
!....         HNO2 + OH = H2O + NO2
!
      qqk(kloop,55)=qk(kloop,55)*y(kloop,7)*y(kloop,25)
!
!....         HO2 + NO2 = HNO4
!
      qqk(kloop,56)=qk(kloop,56)*y(kloop,11)*y(kloop,19)
!
!....         HNO4 = HO2 + NO2
!
      qqk(kloop,57)=qk(kloop,57)*y(kloop,9)
!
!....         HNO4 + OH = H2O + NO2 + O2
!
      qqk(kloop,58)=qk(kloop,58)*y(kloop,9)*y(kloop,25)
!
!....         HO2 + NO3 = NO2 + O2 + OH
!
      qqk(kloop,59)=qk(kloop,59)*y(kloop,11)*y(kloop,20)
!
!....         NO + NO3 = 2 NO2
!
      qqk(kloop,60)=qk(kloop,60)*y(kloop,18)*y(kloop,20)
!
!....         NO3 + OH = HO2 + NO2
!
      qqk(kloop,61)=qk(kloop,61)*y(kloop,20)*y(kloop,25)
!
!....         NO2 + NO3 = N2O5
!
      qqk(kloop,62)=qk(kloop,62)*y(kloop,19)*y(kloop,20)
!
!....         N2O5 = NO2 + NO3
!
      qqk(kloop,63)=qk(kloop,63)*y(kloop,21)
!
!....         HCOOH + OH = H2O + HO2
!
      qqk(kloop,64)=qk(kloop,64)*y(kloop,6)*y(kloop,25)
!
!....         MOH + OH = CH2O + HO2
!
      qqk(kloop,65)=qk(kloop,65)*y(kloop,14)*y(kloop,25)
!
!....         NO2 + NO3 = NO + NO2 + O2
!
      qqk(kloop,66)=qk(kloop,66)*y(kloop,19)*y(kloop,20)
!
!....         CH2O + NO3 = CO + HNO3 + HO2
!
      qqk(kloop,67)=qk(kloop,67)*y(kloop,1)*y(kloop,20)
!
!....         Cl + O3 = ClO + O2
!
      qqk(kloop,68)=qk(kloop,68)*y(kloop,32)*y(kloop,24)
!
!....         Cl + H2 = H + HCl
!
      qqk(kloop,69)=qk(kloop,69)*y(kloop,32)*y(kloop,5)
!
!....         Cl + H2O2 = HCl + HO2
!
      qqk(kloop,70)=qk(kloop,70)*y(kloop,32)*y(kloop,12)
!
!....         Cl + HO2 = HCl + O2
!
      qqk(kloop,71)=qk(kloop,71)*y(kloop,32)*y(kloop,11)
!
!....         Cl + HO2 = ClO + OH
!
      qqk(kloop,72)=qk(kloop,72)*y(kloop,32)*y(kloop,11)
!
!....         ClO + O = Cl + O2
!
      qqk(kloop,73)=qk(kloop,73)*y(kloop,34)*y(kloop,22)
!
!....         ClO + OH = Cl + HO2
!
      qqk(kloop,74)=qk(kloop,74)*y(kloop,34)*y(kloop,25)
!
!....         ClO + OH = HCl + O2
!
      qqk(kloop,75)=qk(kloop,75)*y(kloop,34)*y(kloop,25)
!
!....         ClO + HO2 = HOCl + O2
!
      qqk(kloop,76)=qk(kloop,76)*y(kloop,34)*y(kloop,11)
!
!....         ClO + HO2 = HCl + O3
!
      qqk(kloop,77)=qk(kloop,77)*y(kloop,34)*y(kloop,11)
!
!....         ClO + NO = Cl + NO2
!
      qqk(kloop,78)=qk(kloop,78)*y(kloop,34)*y(kloop,18)
!
!....         ClO + NO2 = ClONO2
!
      qqk(kloop,79)=qk(kloop,79)*y(kloop,34)*y(kloop,19)
!
!....         ClO + ClO = 2 Cl + O2
!
      qqk(kloop,80)=qk(kloop,80)*y(kloop,34)*y(kloop,34)
!
!....         ClO + ClO = Cl2 + O2
!
      qqk(kloop,81)=qk(kloop,81)*y(kloop,34)*y(kloop,34)
!
!....         ClO + ClO = Cl + OClO
!
      qqk(kloop,82)=qk(kloop,82)*y(kloop,34)*y(kloop,34)
!
!....         ClO + ClO = Cl2O2
!
      qqk(kloop,83)=qk(kloop,83)*y(kloop,34)*y(kloop,34)
!
!....         Cl2O2 = 2 ClO
!
      qqk(kloop,84)=qk(kloop,84)*y(kloop,35)
!
!....         HCl + OH = Cl + H2O
!
      qqk(kloop,85)=qk(kloop,85)*y(kloop,37)*y(kloop,25)
!
!....         HOCl + OH = ClO + H2O
!
      qqk(kloop,86)=qk(kloop,86)*y(kloop,38)*y(kloop,25)
!
!....         ClONO2 + O = ClO + NO3
!
      qqk(kloop,87)=qk(kloop,87)*y(kloop,36)*y(kloop,22)
!
!....         ClONO2 + OH = HOCl + NO3
!
      qqk(kloop,88)=qk(kloop,88)*y(kloop,36)*y(kloop,25)
!
!....         Cl + ClONO2 = Cl2 + NO3
!
      qqk(kloop,89)=qk(kloop,89)*y(kloop,32)*y(kloop,36)
!
!....         Br + O3 = BrO + O2
!
      qqk(kloop,90)=qk(kloop,90)*y(kloop,26)*y(kloop,24)
!
!....         Br + HO2 = HBr + O2
!
      qqk(kloop,91)=qk(kloop,91)*y(kloop,26)*y(kloop,11)
!
!....         Br + CH2O = CO + HBr + HO2
!
      qqk(kloop,92)=qk(kloop,92)*y(kloop,26)*y(kloop,1)
!
!....         BrO + O = Br + O2
!
      qqk(kloop,93)=qk(kloop,93)*y(kloop,28)*y(kloop,22)
!
!....         BrO + HO2 = HOBr + O2
!
      qqk(kloop,94)=qk(kloop,94)*y(kloop,28)*y(kloop,11)
!
!....         BrO + NO = Br + NO2
!
      qqk(kloop,95)=qk(kloop,95)*y(kloop,28)*y(kloop,18)
!
!....         BrO + NO2 = BrONO2
!
      qqk(kloop,96)=qk(kloop,96)*y(kloop,28)*y(kloop,19)
!
!....         BrO + ClO = Br + OClO
!
      qqk(kloop,97)=qk(kloop,97)*y(kloop,28)*y(kloop,34)
!
!....         BrO + ClO = Br + Cl + O2
!
      qqk(kloop,98)=qk(kloop,98)*y(kloop,28)*y(kloop,34)
!
!....         BrO + ClO = BrCl + O2
!
      qqk(kloop,99)=qk(kloop,99)*y(kloop,28)*y(kloop,34)
!
!....         BrO + BrO = 2 Br + O2
!
      qqk(kloop,100)=qk(kloop,100)*y(kloop,28)*y(kloop,28)
!
!....         HBr + OH = Br + H2O
!
      qqk(kloop,101)=qk(kloop,101)*y(kloop,30)*y(kloop,25)
!
!....         CH2O + O = CO + HO2 + OH
!
      qqk(kloop,102)=qk(kloop,102)*y(kloop,1)*y(kloop,22)
!
!....         CH4 + Cl = HCl + MO2
!
      qqk(kloop,103)=qk(kloop,103)*y(kloop,2)*y(kloop,32)
!
!....         CH2O + Cl = CO + HCl + HO2
!
      qqk(kloop,104)=qk(kloop,104)*y(kloop,1)*y(kloop,32)
!
!....         CH3Cl + OH = Cl + H2O + HO2
!
      qqk(kloop,105)=qk(kloop,105)*y(kloop,41)*y(kloop,25)
!
!....         CH3CCl3 + OH = 3 Cl + H2O
!
      qqk(kloop,106)=qk(kloop,106)*y(kloop,42)*y(kloop,25)
!
!....         HCFC22 + OH = Cl + H2O
!
      qqk(kloop,107)=qk(kloop,107)*y(kloop,49)*y(kloop,25)
!
!....         HCFC141b + OH = 2 Cl + H2O
!
      qqk(kloop,108)=qk(kloop,108)*y(kloop,50)*y(kloop,25)
!
!....         HCFC142b + OH = Cl + H2O
!
      qqk(kloop,109)=qk(kloop,109)*y(kloop,51)*y(kloop,25)
!
!....         CH3Cl + Cl = CO + 2 HCl + HO2
!
      qqk(kloop,110)=qk(kloop,110)*y(kloop,41)*y(kloop,32)
!
!....         CH3Br + OH = Br + H2O + HO2
!
      qqk(kloop,111)=qk(kloop,111)*y(kloop,40)*y(kloop,25)
!
!....         A3O2 + HO2 = RA3P
!
      qqk(kloop,112)=qk(kloop,112)*y(kloop,56)*y(kloop,11)
!
!....         A3O2 + MO2 =  0.75 CH2O + HO2 +  0.25 MOH +  0.75 RCHO +  0.25 ROH
!
      qqk(kloop,113)=qk(kloop,113)*y(kloop,56)*y(kloop,13)
!
!....         A3O2 + NO = HO2 + NO2 + RCHO
!
      qqk(kloop,114)=qk(kloop,114)*y(kloop,56)*y(kloop,18)
!
!....         ACET + OH = ATO2 + H2O
!
      qqk(kloop,115)=qk(kloop,115)*y(kloop,115)*y(kloop,25)
!
!....         ACTA + OH = H2O + MO2
!
      qqk(kloop,116)=qk(kloop,116)*y(kloop,57)*y(kloop,25)
!
!....         ALD2 + NO3 = HNO3 + MCO3
!
      qqk(kloop,117)=qk(kloop,117)*y(kloop,58)*y(kloop,20)
!
!....         ALD2 + OH =  0.05 CH2O +  0.05 CO + H2O +  0.05 HO2 +  0.95 MCO3
!
      qqk(kloop,118)=qk(kloop,118)*y(kloop,58)*y(kloop,25)
!
!....         ALK4 + NO3 = HNO3 + R4O2
!
      qqk(kloop,119)=qk(kloop,119)*y(kloop,59)*y(kloop,20)
!
!....         ALK4 + OH = R4O2
!
      qqk(kloop,120)=qk(kloop,120)*y(kloop,59)*y(kloop,25)
!
!....         ATO2 + HO2 = MCO3 + MO2
!
      qqk(kloop,121)=qk(kloop,121)*y(kloop,60)*y(kloop,11)
!
!....         ATO2 + MCO3 = ACTA + MEK
!
      qqk(kloop,122)=qk(kloop,122)*y(kloop,60)*y(kloop,84)
!
!....         ATO2 + MCO3 =  0.20 CH2O +  0.80 HO2 +  0.20 MCO3 +  0.80 MGLY + MO2
!
      qqk(kloop,123)=qk(kloop,123)*y(kloop,60)*y(kloop,84)
!
!....         ATO2 + MO2 =  0.50 CH2O +  0.20 HAC +  0.30 HO2 +  0.30 MCO3 +  0.50 MGLY +  0.50 MOH
!
      qqk(kloop,124)=qk(kloop,124)*y(kloop,60)*y(kloop,13)
!
!....         ATO2 + NO =  0.96 CH2O +  0.96 MCO3 +  0.96 NO2 +  0.04 R4N2
!
      qqk(kloop,125)=qk(kloop,125)*y(kloop,60)*y(kloop,18)
!
!....         B3O2 + HO2 = RB3P
!
      qqk(kloop,126)=qk(kloop,126)*y(kloop,61)*y(kloop,11)
!
!....         B3O2 + MCO3 = ACET + ACTA
!
      qqk(kloop,127)=qk(kloop,127)*y(kloop,61)*y(kloop,84)
!
!....         B3O2 + MCO3 = ACET + HO2 + MO2
!
      qqk(kloop,128)=qk(kloop,128)*y(kloop,61)*y(kloop,84)
!
!....         B3O2 + MO2 =  0.75 ACET +  0.75 CH2O + HO2 +  0.25 MOH +  0.25 ROH
!
      qqk(kloop,129)=qk(kloop,129)*y(kloop,61)*y(kloop,13)
!
!....         B3O2 + NO = ACET + HO2 + NO2
!
      qqk(kloop,130)=qk(kloop,130)*y(kloop,61)*y(kloop,18)
!
!....         C2H6 + NO3 = ETO2 + HNO3
!
      qqk(kloop,131)=qk(kloop,131)*y(kloop,62)*y(kloop,20)
!
!....         C2H6 + OH = ETO2 + H2O
!
      qqk(kloop,132)=qk(kloop,132)*y(kloop,62)*y(kloop,25)
!
!....         C2H6 + Cl = ETO2 + HCl
!
      qqk(kloop,133)=qk(kloop,133)*y(kloop,62)*y(kloop,32)
!
!....         C3H8 + OH = A3O2
!
      qqk(kloop,134)=qk(kloop,134)*y(kloop,63)*y(kloop,25)
!
!....         C3H8 + OH = B3O2
!
      qqk(kloop,135)=qk(kloop,135)*y(kloop,63)*y(kloop,25)
!
!....         EOH + OH = ALD2 + HO2
!
      qqk(kloop,136)=qk(kloop,136)*y(kloop,64)*y(kloop,25)
!
!....         ETO2 + ETO2 = 2 ALD2 + 2 HO2
!
      qqk(kloop,137)=qk(kloop,137)*y(kloop,65)*y(kloop,65)
!
!....         ETO2 + ETO2 = ALD2 + EOH
!
      qqk(kloop,138)=qk(kloop,138)*y(kloop,65)*y(kloop,65)
!
!....         ETO2 + NO = ALD2 + HO2 + NO2
!
      qqk(kloop,139)=qk(kloop,139)*y(kloop,65)*y(kloop,18)
!
!....         ETP + OH =  0.50 ALD2 +  0.50 ETO2 +  0.50 OH
!
      qqk(kloop,140)=qk(kloop,140)*y(kloop,66)*y(kloop,25)
!
!....         GLYC + OH =  0.73 CH2O +  0.50 CO +  0.13 GLYX +  0.13 HCOOH +  0.77 HO2 +  0.23 OH
!
      qqk(kloop,141)=qk(kloop,141)*y(kloop,67)*y(kloop,25)
!
!....         GLYC + OH = CO + HCOOH + OH
!
      qqk(kloop,142)=qk(kloop,142)*y(kloop,67)*y(kloop,25)
!
!....         GLYX + NO3 = 2 CO + HNO3 + HO2
!
      qqk(kloop,143)=qk(kloop,143)*y(kloop,68)*y(kloop,20)
!
!....         GLYX + OH = 2 CO + HO2
!
      qqk(kloop,144)=qk(kloop,144)*y(kloop,68)*y(kloop,25)
!
!....         HAC + OH = HO2 + MGLY
!
      qqk(kloop,145)=qk(kloop,145)*y(kloop,69)*y(kloop,25)
!
!....         HAC + OH =  0.50 ACTA +  0.50 CO +  0.50 HCOOH +  0.50 MO2 + OH
!
      qqk(kloop,146)=qk(kloop,146)*y(kloop,69)*y(kloop,25)
!
!....         ETO2 + HO2 = ETP
!
      qqk(kloop,147)=qk(kloop,147)*y(kloop,65)*y(kloop,11)
!
!....         HO2 + MCO3 = ACTA + O3
!
      qqk(kloop,148)=qk(kloop,148)*y(kloop,11)*y(kloop,84)
!
!....         HO2 + MCO3 = MAP
!
      qqk(kloop,149)=qk(kloop,149)*y(kloop,11)*y(kloop,84)
!
!....         IALD + O3 =  0.12 CH2O +  0.28 GLYC +  0.20 GLYX +  0.20 HAC +  0.20 HCOOH +  0.60 MGLY +  0.30 O3 +  0.10 OH
!
      qqk(kloop,150)=qk(kloop,150)*y(kloop,70)*y(kloop,24)
!
!....         IALD + OH =  0.15 HO2 +  0.44 IAO2 +  0.41 MAO3
!
      qqk(kloop,151)=qk(kloop,151)*y(kloop,70)*y(kloop,25)
!
!....         HO2 + IAO2 = IAP
!
      qqk(kloop,152)=qk(kloop,152)*y(kloop,11)*y(kloop,71)
!
!....         IAO2 + MCO3 = ACTA + MEK
!
      qqk(kloop,153)=qk(kloop,153)*y(kloop,71)*y(kloop,84)
!
!....         IAO2 + MCO3 =  0.40 CH2O +  0.29 CO +  0.26 GLYC +  0.18 GLYX +  0.36 HAC + HO2 +  0.58 MGLY + MO2
!
      qqk(kloop,154)=qk(kloop,154)*y(kloop,71)*y(kloop,84)
!
!....         IAO2 + MO2 =  0.95 CH2O +  0.15 CO +  0.13 GLYC +  0.09 GLYX +  0.18 HAC + HO2 +  0.25 MEK +  0.29 MGLY +  0.25 MOH +  0.25 ROH
!
      qqk(kloop,155)=qk(kloop,155)*y(kloop,71)*y(kloop,13)
!
!....         IAO2 + NO =  0.35 CH2O +  0.27 CO +  0.24 GLYC +  0.17 GLYX +  0.33 HAC +  0.08 HNO3 +  0.92 HO2 +  0.53 MGLY +  0.92 NO2
!
      qqk(kloop,156)=qk(kloop,156)*y(kloop,71)*y(kloop,18)
!
!....         IAP + OH =  0.50 IAO2 +  0.50 OH +  0.50 RCHO
!
      qqk(kloop,157)=qk(kloop,157)*y(kloop,72)*y(kloop,25)
!
!....         HO2 + INO2 = INPN
!
      qqk(kloop,158)=qk(kloop,158)*y(kloop,11)*y(kloop,73)
!
!....         INO2 + MCO3 =  0.15 CH2O +  0.85 HNO3 +  0.80 HO2 +  0.10 MACR + MO2 +  0.05 MVK +  0.15 NO2
!
      qqk(kloop,159)=qk(kloop,159)*y(kloop,73)*y(kloop,84)
!
!....         INO2 + MCO3 = ACTA + NO2 + RCHO
!
      qqk(kloop,160)=qk(kloop,160)*y(kloop,73)*y(kloop,84)
!
!....         INO2 + MO2 =  0.83 CH2O +  0.43 HNO3 +  0.90 HO2 +  0.05 MACR +  0.25 MOH +  0.03 MVK +  0.57 NO2 +  0.25 RCHO +  0.25 ROH
!
      qqk(kloop,161)=qk(kloop,161)*y(kloop,73)*y(kloop,13)
!
!....         INO2 + NO =  0.15 CH2O +  0.85 HNO3 +  0.80 HO2 +  0.10 MACR +  0.05 MVK +  1.15 NO2
!
      qqk(kloop,162)=qk(kloop,162)*y(kloop,73)*y(kloop,18)
!
!....         INPN + OH = INO2
!
      qqk(kloop,163)=qk(kloop,163)*y(kloop,74)*y(kloop,25)
!
!....         HO2 + ISN1 = ISNP
!
      qqk(kloop,164)=qk(kloop,164)*y(kloop,11)*y(kloop,75)
!
!....         ISN1 + MCO3 = GLYC + HAC + MO2 + NO2
!
      qqk(kloop,165)=qk(kloop,165)*y(kloop,75)*y(kloop,84)
!
!....         ISN1 + MCO3 = ACTA + NO2 + RCHO
!
      qqk(kloop,166)=qk(kloop,166)*y(kloop,75)*y(kloop,84)
!
!....         ISN1 + MO2 =  0.75 CH2O +  0.50 GLYC +  0.50 HAC +  0.50 HO2 +  0.25 MOH + NO2 +  0.25 RCHO +  0.25 ROH
!
      qqk(kloop,167)=qk(kloop,167)*y(kloop,75)*y(kloop,13)
!
!....         ISNP + OH =  0.50 ISN1 +  0.50 NO2 +  0.50 OH +  0.50 RCHO
!
      qqk(kloop,168)=qk(kloop,168)*y(kloop,76)*y(kloop,25)
!
!....         ISOP + NO3 = INO2
!
      qqk(kloop,169)=qk(kloop,169)*y(kloop,77)*y(kloop,20)
!
!....         ISOP + O3 =  0.90 CH2O +  0.05 CO +  0.06 HO2 +  0.39 MACR +  0.16 MVK +  0.10 O3 +  0.27 OH +  0.07 PRPE
!
      qqk(kloop,170)=qk(kloop,170)*y(kloop,77)*y(kloop,24)
!
!....         ISOP + OH = RIO2
!
      qqk(kloop,171)=qk(kloop,171)*y(kloop,77)*y(kloop,25)
!
!....         HO2 + KO2 = MGLY + MO2
!
      qqk(kloop,172)=qk(kloop,172)*y(kloop,11)*y(kloop,78)
!
!....         KO2 + MCO3 = ACTA + MEK
!
      qqk(kloop,173)=qk(kloop,173)*y(kloop,78)*y(kloop,84)
!
!....         KO2 + MCO3 = ALD2 + MCO3 + MO2
!
      qqk(kloop,174)=qk(kloop,174)*y(kloop,78)*y(kloop,84)
!
!....         KO2 + MO2 =  0.50 ALD2 +  0.75 CH2O +  0.50 HO2 +  0.50 MCO3 +  0.25 MEK +  0.25 MOH +  0.25 ROH
!
      qqk(kloop,175)=qk(kloop,175)*y(kloop,78)*y(kloop,13)
!
!....         KO2 + NO =  0.93 ALD2 +  0.93 MCO3 +  0.93 NO2 +  0.07 R4N2
!
      qqk(kloop,176)=qk(kloop,176)*y(kloop,78)*y(kloop,18)
!
!....         MACR + NO3 = MAN2
!
      qqk(kloop,177)=qk(kloop,177)*y(kloop,79)*y(kloop,20)
!
!....         MACR + NO3 = HNO3 + MAO3
!
      qqk(kloop,178)=qk(kloop,178)*y(kloop,79)*y(kloop,20)
!
!....         MACR + O3 =  0.70 CH2O +  0.20 CO +  0.28 HO2 +  0.80 MGLY +  0.20 O3 +  0.22 OH
!
      qqk(kloop,179)=qk(kloop,179)*y(kloop,79)*y(kloop,24)
!
!....         MACR + OH =  0.53 MAO3 +  0.47 MRO2
!
      qqk(kloop,180)=qk(kloop,180)*y(kloop,79)*y(kloop,25)
!
!....         HO2 + MAN2 = ISNP
!
      qqk(kloop,181)=qk(kloop,181)*y(kloop,11)*y(kloop,80)
!
!....         MAN2 + MCO3 = CH2O + MGLY + MO2 + NO2
!
      qqk(kloop,182)=qk(kloop,182)*y(kloop,80)*y(kloop,84)
!
!....         MAN2 + MCO3 = ACTA + NO2 + RCHO
!
      qqk(kloop,183)=qk(kloop,183)*y(kloop,80)*y(kloop,84)
!
!....         MAN2 + MO2 =  1.25 CH2O +  0.50 HO2 +  0.50 MGLY +  0.25 MOH + NO2 +  0.25 RCHO +  0.25 ROH
!
      qqk(kloop,184)=qk(kloop,184)*y(kloop,80)*y(kloop,13)
!
!....         MAN2 + NO = CH2O + MGLY + 2 NO2
!
      qqk(kloop,185)=qk(kloop,185)*y(kloop,80)*y(kloop,18)
!
!....         HO2 + MAO3 =  0.59 CH2O +  0.39 CO +  0.41 MAOP +  0.39 MO2 +  0.15 O3 +  0.44 OH
!
      qqk(kloop,186)=qk(kloop,186)*y(kloop,11)*y(kloop,81)
!
!....         MAO3 + MCO3 = CH2O + MCO3 + MO2
!
      qqk(kloop,187)=qk(kloop,187)*y(kloop,81)*y(kloop,84)
!
!....         MAO3 + MO2 = 2 CH2O + HO2 + MCO3
!
      qqk(kloop,188)=qk(kloop,188)*y(kloop,81)*y(kloop,13)
!
!....         MAO3 + MO2 = CH2O + RCOOH
!
      qqk(kloop,189)=qk(kloop,189)*y(kloop,81)*y(kloop,13)
!
!....         MAO3 + NO2 = PMN
!
      qqk(kloop,190)=qk(kloop,190)*y(kloop,81)*y(kloop,19)
!
!....         MAO3 + NO =  0.50 CH2O +  0.50 CO +  0.50 MCO3 +  0.50 MO2 + NO2
!
      qqk(kloop,191)=qk(kloop,191)*y(kloop,81)*y(kloop,18)
!
!....         MAOP + OH = MAO3
!
      qqk(kloop,192)=qk(kloop,192)*y(kloop,82)*y(kloop,25)
!
!....         A3O2 + MCO3 = ACTA + RCHO
!
      qqk(kloop,193)=qk(kloop,193)*y(kloop,56)*y(kloop,84)
!
!....         A3O2 + MCO3 = HO2 + MO2 + RCHO
!
      qqk(kloop,194)=qk(kloop,194)*y(kloop,56)*y(kloop,84)
!
!....         ETO2 + MCO3 = ACTA + ALD2
!
      qqk(kloop,195)=qk(kloop,195)*y(kloop,65)*y(kloop,84)
!
!....         ETO2 + MCO3 = ALD2 + HO2 + MO2
!
      qqk(kloop,196)=qk(kloop,196)*y(kloop,65)*y(kloop,84)
!
!....         MCO3 + MCO3 = 2 MO2
!
      qqk(kloop,197)=qk(kloop,197)*y(kloop,84)*y(kloop,84)
!
!....         MCO3 + MO2 = ACTA + CH2O
!
      qqk(kloop,198)=qk(kloop,198)*y(kloop,84)*y(kloop,13)
!
!....         MCO3 + MO2 = CH2O + HO2 + MO2
!
      qqk(kloop,199)=qk(kloop,199)*y(kloop,84)*y(kloop,13)
!
!....         MCO3 + NO2 = PAN
!
      qqk(kloop,200)=qk(kloop,200)*y(kloop,84)*y(kloop,19)
!
!....         MCO3 + NO = MO2 + NO2
!
      qqk(kloop,201)=qk(kloop,201)*y(kloop,84)*y(kloop,18)
!
!....         MCO3 + PO2 = ACTA +  0.65 HAC +  0.35 RCHO
!
      qqk(kloop,202)=qk(kloop,202)*y(kloop,84)*y(kloop,92)
!
!....         MCO3 + PO2 = ALD2 + CH2O + HO2 + MO2
!
      qqk(kloop,203)=qk(kloop,203)*y(kloop,84)*y(kloop,92)
!
!....         MEK + NO3 = HNO3 + KO2
!
      qqk(kloop,204)=qk(kloop,204)*y(kloop,85)*y(kloop,20)
!
!....         MEK + OH = H2O + KO2
!
      qqk(kloop,205)=qk(kloop,205)*y(kloop,85)*y(kloop,25)
!
!....         MGLY + NO3 = CO + HNO3 + MCO3
!
      qqk(kloop,206)=qk(kloop,206)*y(kloop,86)*y(kloop,20)
!
!....         MGLY + OH = CO + MCO3
!
      qqk(kloop,207)=qk(kloop,207)*y(kloop,86)*y(kloop,25)
!
!....         ETO2 + MO2 =  0.75 ALD2 +  0.75 CH2O +  0.25 EOH + HO2 +  0.25 MOH
!
      qqk(kloop,208)=qk(kloop,208)*y(kloop,65)*y(kloop,13)
!
!....         HO2 + MRO2 = MRP
!
      qqk(kloop,209)=qk(kloop,209)*y(kloop,11)*y(kloop,87)
!
!....         MCO3 + MRO2 = ACTA + MEK
!
      qqk(kloop,210)=qk(kloop,210)*y(kloop,84)*y(kloop,87)
!
!....         MCO3 + MRO2 =  0.17 CH2O +  0.83 CO +  0.83 HAC + HO2 +  0.17 MGLY + MO2
!
      qqk(kloop,211)=qk(kloop,211)*y(kloop,84)*y(kloop,87)
!
!....         MO2 + MRO2 = CH2O +  0.60 CO +  0.60 HAC + HO2 +  0.25 MGLY +  0.15 ROH
!
      qqk(kloop,212)=qk(kloop,212)*y(kloop,13)*y(kloop,87)
!
!....         MRO2 + NO = HNO3
!
      qqk(kloop,213)=qk(kloop,213)*y(kloop,87)*y(kloop,18)
!
!....         MRO2 + NO =  0.17 CH2O +  0.83 CO +  0.83 HAC + HO2 +  0.17 MGLY + NO2
!
      qqk(kloop,214)=qk(kloop,214)*y(kloop,87)*y(kloop,18)
!
!....         MRP + OH = MRO2
!
      qqk(kloop,215)=qk(kloop,215)*y(kloop,88)*y(kloop,25)
!
!....         MVK + O3 =  0.04 ALD2 +  0.80 CH2O +  0.05 CO +  0.06 HO2 +  0.82 MGLY +  0.20 O3 +  0.08 OH
!
      qqk(kloop,216)=qk(kloop,216)*y(kloop,89)*y(kloop,24)
!
!....         MVK + OH = VRO2
!
      qqk(kloop,217)=qk(kloop,217)*y(kloop,89)*y(kloop,25)
!
!....         MAP + OH =  0.50 CH2O +  0.50 MCO3 +  0.50 OH
!
      qqk(kloop,218)=qk(kloop,218)*y(kloop,83)*y(kloop,25)
!
!....         OH + RCHO = H2O + RCO3
!
      qqk(kloop,219)=qk(kloop,219)*y(kloop,25)*y(kloop,104)
!
!....         PAN = MCO3 + NO2
!
      qqk(kloop,220)=qk(kloop,220)*y(kloop,90)
!
!....         PMN = MAO3 + NO2
!
      qqk(kloop,221)=qk(kloop,221)*y(kloop,91)
!
!....         O3 + PMN =  0.60 CH2O + HO2 + NO2
!
      qqk(kloop,222)=qk(kloop,222)*y(kloop,24)*y(kloop,91)
!
!....         OH + PMN = CO + HAC + NO2
!
      qqk(kloop,223)=qk(kloop,223)*y(kloop,25)*y(kloop,91)
!
!....         HO2 + PO2 = PP
!
      qqk(kloop,224)=qk(kloop,224)*y(kloop,11)*y(kloop,92)
!
!....         MO2 + PO2 =  0.50 ALD2 + CH2O +  0.16 HAC + HO2 +  0.25 MOH +  0.09 RCHO +  0.25 ROH
!
      qqk(kloop,225)=qk(kloop,225)*y(kloop,13)*y(kloop,92)
!
!....         NO + PO2 = ALD2 + CH2O + HO2 + NO2
!
      qqk(kloop,226)=qk(kloop,226)*y(kloop,18)*y(kloop,92)
!
!....         PPN = NO2 + RCO3
!
      qqk(kloop,227)=qk(kloop,227)*y(kloop,94)
!
!....         OH + PP =  0.79 HAC +  0.79 OH +  0.21 PO2
!
      qqk(kloop,228)=qk(kloop,228)*y(kloop,25)*y(kloop,93)
!
!....         HO2 + PRN1 = PRPN
!
      qqk(kloop,229)=qk(kloop,229)*y(kloop,11)*y(kloop,95)
!
!....         MCO3 + PRN1 = ALD2 + CH2O + MO2 + NO2
!
      qqk(kloop,230)=qk(kloop,230)*y(kloop,84)*y(kloop,95)
!
!....         MCO3 + PRN1 = ACTA + NO2 + RCHO
!
      qqk(kloop,231)=qk(kloop,231)*y(kloop,84)*y(kloop,95)
!
!....         MO2 + PRN1 =  0.50 ALD2 +  1.25 CH2O +  0.50 HO2 +  0.25 MOH + NO2 +  0.25 RCHO +  0.25 ROH
!
      qqk(kloop,232)=qk(kloop,232)*y(kloop,13)*y(kloop,95)
!
!....         NO + PRN1 = ALD2 + CH2O + 2 NO2
!
      qqk(kloop,233)=qk(kloop,233)*y(kloop,18)*y(kloop,95)
!
!....         NO3 + PRPE = PRN1
!
      qqk(kloop,234)=qk(kloop,234)*y(kloop,20)*y(kloop,96)
!
!....         O3 + PRPE =  0.50 ALD2 +  0.54 CH2O +  0.42 CO +  0.06 H2 +  0.30 HO2 +  0.31 MO2 +  0.14 OH
!
      qqk(kloop,235)=qk(kloop,235)*y(kloop,24)*y(kloop,96)
!
!....         OH + PRPE = PO2
!
      qqk(kloop,236)=qk(kloop,236)*y(kloop,25)*y(kloop,96)
!
!....         OH + PRPN = PRN1
!
      qqk(kloop,237)=qk(kloop,237)*y(kloop,25)*y(kloop,97)
!
!....         HO2 + R4N1 = R4N2
!
      qqk(kloop,238)=qk(kloop,238)*y(kloop,11)*y(kloop,98)
!
!....         MCO3 + R4N1 =  0.75 ALD2 +  0.39 CH2O + MO2 + NO2 +  0.30 R4O2 +  0.57 RCHO
!
      qqk(kloop,239)=qk(kloop,239)*y(kloop,84)*y(kloop,98)
!
!....         MCO3 + R4N1 = ACTA + NO2 + RCHO
!
      qqk(kloop,240)=qk(kloop,240)*y(kloop,84)*y(kloop,98)
!
!....         MO2 + R4N1 =  0.38 ALD2 +  0.95 CH2O +  0.50 HO2 +  0.25 MOH + NO2 +  0.15 R4O2 +  0.54 RCHO +  0.25 ROH
!
      qqk(kloop,241)=qk(kloop,241)*y(kloop,13)*y(kloop,98)
!
!....         NO + R4N1 =  0.75 ALD2 +  0.39 CH2O + 2 NO2 +  0.30 R4O2 +  0.57 RCHO
!
      qqk(kloop,242)=qk(kloop,242)*y(kloop,18)*y(kloop,98)
!
!....         OH + R4N2 = H2O + R4N1
!
      qqk(kloop,243)=qk(kloop,243)*y(kloop,25)*y(kloop,99)
!
!....         HO2 + R4O2 = R4P
!
      qqk(kloop,244)=qk(kloop,244)*y(kloop,11)*y(kloop,100)
!
!....         MCO3 + R4O2 = ACTA + MEK
!
      qqk(kloop,245)=qk(kloop,245)*y(kloop,84)*y(kloop,100)
!
!....         MCO3 + R4O2 =  0.05 A3O2 +  0.32 ACET +  0.32 ALD2 +  0.18 B3O2 +  0.32 ETO2 +  0.27 HO2 +  0.19 MEK +  1.18 MO2 +  0.13 RCHO
!
      qqk(kloop,246)=qk(kloop,246)*y(kloop,84)*y(kloop,100)
!
!....         MO2 + R4O2 =  0.03 A3O2 +  0.16 ACET +  0.16 ALD2 +  0.09 B3O2 +  0.75 CH2O +  0.16 ETO2 +  0.64 HO2 +  0.35 MEK +  0.09 MO2 +  0.25 MOH +  0.07 RCHO +  0.25 ROH
!
      qqk(kloop,247)=qk(kloop,247)*y(kloop,13)*y(kloop,100)
!
!....         NO + R4O2 =  0.05 A3O2 +  0.32 ACET +  0.32 ALD2 +  0.18 B3O2 +  0.32 ETO2 +  0.27 HO2 +  0.19 MEK +  0.18 MO2 + NO2 +  0.13 RCHO
!
      qqk(kloop,248)=qk(kloop,248)*y(kloop,18)*y(kloop,100)
!
!....         NO + R4O2 = R4N2
!
      qqk(kloop,249)=qk(kloop,249)*y(kloop,18)*y(kloop,100)
!
!....         OH + R4P =  0.50 OH +  0.50 R4O2 +  0.50 RCHO
!
      qqk(kloop,250)=qk(kloop,250)*y(kloop,25)*y(kloop,101)
!
!....         OH + RA3P =  0.50 A3O2 +  0.50 OH +  0.50 RCHO
!
      qqk(kloop,251)=qk(kloop,251)*y(kloop,25)*y(kloop,102)
!
!....         OH + RB3P =  0.79 ACET +  0.21 B3O2 +  0.79 OH
!
      qqk(kloop,252)=qk(kloop,252)*y(kloop,25)*y(kloop,103)
!
!....         NO3 + RCHO = HNO3 + RCO3
!
      qqk(kloop,253)=qk(kloop,253)*y(kloop,20)*y(kloop,104)
!
!....         OH + RCOOH = ETO2
!
      qqk(kloop,254)=qk(kloop,254)*y(kloop,25)*y(kloop,106)
!
!....         HO2 + RCO3 =  0.30 O3 +  0.30 RCOOH +  0.70 RP
!
      qqk(kloop,255)=qk(kloop,255)*y(kloop,11)*y(kloop,105)
!
!....         MCO3 + RCO3 = ETO2 + MO2
!
      qqk(kloop,256)=qk(kloop,256)*y(kloop,84)*y(kloop,105)
!
!....         MO2 + RCO3 = CH2O + ETO2 + HO2
!
      qqk(kloop,257)=qk(kloop,257)*y(kloop,13)*y(kloop,105)
!
!....         MO2 + RCO3 = CH2O + RCOOH
!
      qqk(kloop,258)=qk(kloop,258)*y(kloop,13)*y(kloop,105)
!
!....         NO2 + RCO3 = PPN
!
      qqk(kloop,259)=qk(kloop,259)*y(kloop,19)*y(kloop,105)
!
!....         NO + RCO3 = ETO2 + NO2
!
      qqk(kloop,260)=qk(kloop,260)*y(kloop,18)*y(kloop,105)
!
!....         HO2 + RIO1 = RIP
!
      qqk(kloop,261)=qk(kloop,261)*y(kloop,11)*y(kloop,107)
!
!....         MCO3 + RIO1 = ACTA + MEK
!
      qqk(kloop,262)=qk(kloop,262)*y(kloop,84)*y(kloop,107)
!
!....         MCO3 + RIO1 =  0.75 CH2O + HO2 + IALD + MO2
!
      qqk(kloop,263)=qk(kloop,263)*y(kloop,84)*y(kloop,107)
!
!....         MO2 + RIO1 =  1.13 CH2O + HO2 +  0.50 IALD +  0.25 MEK +  0.25 MOH +  0.25 ROH
!
      qqk(kloop,264)=qk(kloop,264)*y(kloop,13)*y(kloop,107)
!
!....         NO + RIO1 = HNO3
!
      qqk(kloop,265)=qk(kloop,265)*y(kloop,18)*y(kloop,107)
!
!....         NO + RIO1 =  0.75 CH2O + HO2 + IALD + NO2
!
      qqk(kloop,266)=qk(kloop,266)*y(kloop,18)*y(kloop,107)
!
!....         HO2 + RIO2 = RIP
!
      qqk(kloop,267)=qk(kloop,267)*y(kloop,11)*y(kloop,108)
!
!....         MCO3 + RIO2 = ACTA + MEK
!
      qqk(kloop,268)=qk(kloop,268)*y(kloop,84)*y(kloop,108)
!
!....         MCO3 + RIO2 =  0.69 CH2O +  0.86 HO2 +  0.13 IALD +  0.29 MACR + MO2 +  0.40 MVK +  0.14 RIO1
!
      qqk(kloop,269)=qk(kloop,269)*y(kloop,84)*y(kloop,108)
!
!....         MO2 + RIO2 =  1.10 CH2O +  0.93 HO2 +  0.06 IALD +  0.14 MACR +  0.25 MEK +  0.25 MOH +  0.20 MVK +  0.07 RIO1 +  0.25 ROH
!
      qqk(kloop,270)=qk(kloop,270)*y(kloop,13)*y(kloop,108)
!
!....         NO + RIO2 = HNO3
!
      qqk(kloop,271)=qk(kloop,271)*y(kloop,18)*y(kloop,108)
!
!....         NO + RIO2 =  0.69 CH2O +  0.86 HO2 +  0.13 IALD +  0.29 MACR +  0.40 MVK + NO2 +  0.14 RIO1
!
      qqk(kloop,272)=qk(kloop,272)*y(kloop,18)*y(kloop,108)
!
!....         OH + RIP =  0.50 IAO2 +  0.10 RIO1 +  0.40 RIO2
!
      qqk(kloop,273)=qk(kloop,273)*y(kloop,25)*y(kloop,109)
!
!....         OH + ROH = HO2 + RCHO
!
      qqk(kloop,274)=qk(kloop,274)*y(kloop,25)*y(kloop,110)
!
!....         OH + RP =  0.50 ALD2 +  0.50 OH +  0.50 RCO3
!
      qqk(kloop,275)=qk(kloop,275)*y(kloop,25)*y(kloop,111)
!
!....         HO2 + VRO2 =  0.10 CH2O +  0.58 GLYC +  0.10 HO2 +  0.58 MCO3 +  0.19 MEK +  0.10 MGLY +  0.68 OH +  0.03 RCHO +  0.10 VRP
!
      qqk(kloop,276)=qk(kloop,276)*y(kloop,11)*y(kloop,112)
!
!....         MCO3 + VRO2 = ACTA + MEK
!
      qqk(kloop,277)=qk(kloop,277)*y(kloop,84)*y(kloop,112)
!
!....         MCO3 + VRO2 =  0.28 CH2O +  0.72 GLYC +  0.28 HO2 +  0.72 MCO3 +  0.28 MGLY + MO2
!
      qqk(kloop,278)=qk(kloop,278)*y(kloop,84)*y(kloop,112)
!
!....         MO2 + VRO2 =  0.89 CH2O +  0.36 GLYC +  0.64 HO2 +  0.36 MCO3 +  0.25 MEK +  0.14 MGLY +  0.25 MOH +  0.25 ROH
!
      qqk(kloop,279)=qk(kloop,279)*y(kloop,13)*y(kloop,112)
!
!....         NO + VRO2 = HNO3
!
      qqk(kloop,280)=qk(kloop,280)*y(kloop,18)*y(kloop,112)
!
!....         NO + VRO2 =  0.28 CH2O +  0.72 GLYC +  0.28 HO2 +  0.72 MCO3 +  0.28 MGLY + NO2
!
      qqk(kloop,281)=qk(kloop,281)*y(kloop,18)*y(kloop,112)
!
!....         OH + VRP =  0.50 OH +  0.50 RCHO +  0.50 VRO2
!
      qqk(kloop,282)=qk(kloop,282)*y(kloop,25)*y(kloop,113)
!
!....         N2O5 = 2 HNO3
!
      qqk(kloop,283)=qk(kloop,283)*y(kloop,21)
!
!....         ClONO2 = HNO3 + HOCl
!
      qqk(kloop,284)=qk(kloop,284)*y(kloop,36)
!
!....         BrONO2 = HNO3 + HOBr
!
      qqk(kloop,285)=qk(kloop,285)*y(kloop,29)
!
!....         ClONO2 + HCl = Cl2 + HNO3
!
      qqk(kloop,286)=qk(kloop,286)*y(kloop,36)*y(kloop,37)
!
!....         HCl + HOCl = Cl2 + H2O
!
      qqk(kloop,287)=qk(kloop,287)*y(kloop,37)*y(kloop,38)
!
!....         HCl + HOBr = BrCl + H2O
!
      qqk(kloop,288)=qk(kloop,288)*y(kloop,37)*y(kloop,31)
!
!....         N2O5 = 2 HNO3
!
      qqk(kloop,289)=qk(kloop,289)*y(kloop,21)
!
!....         ClONO2 = HNO3 + HOCl
!
      qqk(kloop,290)=qk(kloop,290)*y(kloop,36)
!
!....         BrONO2 = HNO3 + HOBr
!
      qqk(kloop,291)=qk(kloop,291)*y(kloop,29)
!
!....         ClONO2 + HCl = Cl2 + HNO3
!
      qqk(kloop,292)=qk(kloop,292)*y(kloop,36)*y(kloop,37)
!
!....         HCl + HOCl = Cl2 + H2O
!
      qqk(kloop,293)=qk(kloop,293)*y(kloop,37)*y(kloop,38)
!
!....         HCl + HOBr = BrCl + H2O
!
      qqk(kloop,294)=qk(kloop,294)*y(kloop,37)*y(kloop,31)
!
!....         ClONO2 = HNO3 + HOCl
!
      qqk(kloop,295)=qk(kloop,295)*y(kloop,36)
!
!....         BrONO2 = HNO3 + HOBr
!
      qqk(kloop,296)=qk(kloop,296)*y(kloop,29)
!
!....         ClONO2 + HCl = Cl2 + HNO3
!
      qqk(kloop,297)=qk(kloop,297)*y(kloop,36)*y(kloop,37)
!
!....         HCl + HOCl = Cl2 + H2O
!
      qqk(kloop,298)=qk(kloop,298)*y(kloop,37)*y(kloop,38)
!
!....         BrONO2 + HCl = BrCl + HNO3
!
      qqk(kloop,299)=qk(kloop,299)*y(kloop,29)*y(kloop,37)
!
!....         HCl + HOBr = BrCl + H2O
!
      qqk(kloop,300)=qk(kloop,300)*y(kloop,37)*y(kloop,31)
!
!....         ClONO2 = HNO3 + HOCl
!
      qqk(kloop,301)=qk(kloop,301)*y(kloop,36)
!
!....         BrONO2 = HNO3 + HOBr
!
      qqk(kloop,302)=qk(kloop,302)*y(kloop,29)
!
!....         ClONO2 + HCl = Cl2 + HNO3
!
      qqk(kloop,303)=qk(kloop,303)*y(kloop,36)*y(kloop,37)
!
!....         HCl + HOCl = Cl2 + H2O
!
      qqk(kloop,304)=qk(kloop,304)*y(kloop,37)*y(kloop,38)
!
!....         BrONO2 + HCl = BrCl + HNO3
!
      qqk(kloop,305)=qk(kloop,305)*y(kloop,29)*y(kloop,37)
!
!....         HCl + HOBr = BrCl + H2O
!
      qqk(kloop,306)=qk(kloop,306)*y(kloop,37)*y(kloop,31)
!
!....         HNO3 = NO2 + OH
!
      qqk(kloop,307)=qk(kloop,307)*y(kloop,8)
!
!....         NO3 + NO3 = 2 NO2 + O2
!
      qqk(kloop,308)=qk(kloop,308)*y(kloop,20)*y(kloop,20)
!
!....         HO2 =  0.50 H2O
!
      qqk(kloop,309)=qk(kloop,309)*y(kloop,11)
!
!....         NO2 =  0.50 HNO2 +  0.50 HNO3
!
      qqk(kloop,310)=qk(kloop,310)*y(kloop,19)
!
!....         NO3 = HNO3
!
      qqk(kloop,311)=qk(kloop,311)*y(kloop,20)
!
!....         N2O5 = 2 HNO3
!
      qqk(kloop,312)=qk(kloop,312)*y(kloop,21)
!
!                  start photolytic reactions
!
!....  O2 + hv = 2 O
!
      qqj(kloop,1)=qj(kloop,1)*y(kloop,117)
!
!....  O3 + hv = O + O2
!
      qqj(kloop,2)=qj(kloop,2)*y(kloop,24)
!
!....  O3 + hv = O1D + O2
!
      qqj(kloop,3)=qj(kloop,3)*y(kloop,24)
!
!....  NO + hv = N + O
!
      qqj(kloop,4)=qj(kloop,4)*y(kloop,18)
!
!....  N2O + hv = N2 + O1D
!
      qqj(kloop,5)=qj(kloop,5)*y(kloop,17)
!
!....  NO2 + hv = NO + O
!
      qqj(kloop,6)=qj(kloop,6)*y(kloop,19)
!
!....  H2O2 + hv = 2 OH
!
      qqj(kloop,7)=qj(kloop,7)*y(kloop,12)
!
!....  MP + hv = CH2O + HO2 + OH
!
      qqj(kloop,8)=qj(kloop,8)*y(kloop,15)
!
!....  CH2O + hv = CO + H + HO2
!
      qqj(kloop,9)=qj(kloop,9)*y(kloop,1)
!
!....  CH2O + hv = CO + H2
!
      qqj(kloop,10)=qj(kloop,10)*y(kloop,1)
!
!....  HNO3 + hv = NO2 + OH
!
      qqj(kloop,11)=qj(kloop,11)*y(kloop,8)
!
!....  HNO2 + hv = NO + OH
!
      qqj(kloop,12)=qj(kloop,12)*y(kloop,7)
!
!....  HNO4 + hv = NO3 + OH
!
      qqj(kloop,13)=qj(kloop,13)*y(kloop,9)
!
!....  NO3 + hv = NO2 + O3
!
      qqj(kloop,14)=qj(kloop,14)*y(kloop,20)
!
!....  NO3 + hv = NO + O2
!
      qqj(kloop,15)=qj(kloop,15)*y(kloop,20)
!
!....  N2O5 + hv = NO2 + NO3
!
      qqj(kloop,16)=qj(kloop,16)*y(kloop,21)
!
!....  HNO4 + hv = HO2 + NO2
!
      qqj(kloop,17)=qj(kloop,17)*y(kloop,9)
!
!....  Cl2 + hv = 2 Cl
!
      qqj(kloop,18)=qj(kloop,18)*y(kloop,33)
!
!....  OClO + hv = ClO + O
!
      qqj(kloop,19)=qj(kloop,19)*y(kloop,39)
!
!....  Cl2O2 + hv = 2 Cl + O2
!
      qqj(kloop,20)=qj(kloop,20)*y(kloop,35)
!
!....  HOCl + hv = Cl + OH
!
      qqj(kloop,21)=qj(kloop,21)*y(kloop,38)
!
!....  ClONO2 + hv = Cl + NO3
!
      qqj(kloop,22)=qj(kloop,22)*y(kloop,36)
!
!....  ClONO2 + hv = ClO + NO2
!
      qqj(kloop,23)=qj(kloop,23)*y(kloop,36)
!
!....  BrCl + hv = Br + Cl
!
      qqj(kloop,24)=qj(kloop,24)*y(kloop,27)
!
!....  BrO + hv = Br + O
!
      qqj(kloop,25)=qj(kloop,25)*y(kloop,28)
!
!....  HOBr + hv = Br + OH
!
      qqj(kloop,26)=qj(kloop,26)*y(kloop,31)
!
!....  BrONO2 + hv = Br + NO3
!
      qqj(kloop,27)=qj(kloop,27)*y(kloop,29)
!
!....  CH3Cl + hv = Cl + MO2
!
      qqj(kloop,28)=qj(kloop,28)*y(kloop,41)
!
!....  CCl4 + hv = 4 Cl
!
      qqj(kloop,29)=qj(kloop,29)*y(kloop,43)
!
!....  CH3CCl3 + hv = 3 Cl
!
      qqj(kloop,30)=qj(kloop,30)*y(kloop,42)
!
!....  CFC11 + hv = 3 Cl
!
      qqj(kloop,31)=qj(kloop,31)*y(kloop,44)
!
!....  CFC12 + hv = 2 Cl
!
      qqj(kloop,32)=qj(kloop,32)*y(kloop,45)
!
!....  CFC113 + hv = 3 Cl
!
      qqj(kloop,33)=qj(kloop,33)*y(kloop,46)
!
!....  CFC114 + hv = 2 Cl
!
      qqj(kloop,34)=qj(kloop,34)*y(kloop,47)
!
!....  CFC115 + hv = Cl
!
      qqj(kloop,35)=qj(kloop,35)*y(kloop,48)
!
!....  HCFC141b + hv = 2 Cl
!
      qqj(kloop,36)=qj(kloop,36)*y(kloop,50)
!
!....  HCFC142b + hv = Cl
!
      qqj(kloop,37)=qj(kloop,37)*y(kloop,51)
!
!....  CH3Br + hv = Br + MO2
!
      qqj(kloop,38)=qj(kloop,38)*y(kloop,40)
!
!....  CF3Br + hv = Br
!
      qqj(kloop,39)=qj(kloop,39)*y(kloop,54)
!
!....  CF2Br2 + hv = 2 Br
!
      qqj(kloop,40)=qj(kloop,40)*y(kloop,52)
!
!....  H2402 + hv = 2 Br
!
      qqj(kloop,41)=qj(kloop,41)*y(kloop,55)
!
!....  CF2ClBr + hv = Br + Cl
!
      qqj(kloop,42)=qj(kloop,42)*y(kloop,53)
!
!....  ALD2 + hv = CO + HO2 + MO2
!
      qqj(kloop,43)=qj(kloop,43)*y(kloop,58)
!
!....  PAN + hv = MCO3 + NO2
!
      qqj(kloop,44)=qj(kloop,44)*y(kloop,90)
!
!....  RCHO + hv = CO + ETO2 + HO2
!
      qqj(kloop,45)=qj(kloop,45)*y(kloop,104)
!
!....  ACET + hv = MCO3 + MO2
!
      qqj(kloop,46)=qj(kloop,46)*y(kloop,115)
!
!....  MEK + hv = ETO2 + MCO3
!
      qqj(kloop,47)=qj(kloop,47)*y(kloop,85)
!
!....  GLYC + hv = CH2O + CO + 2 HO2
!
      qqj(kloop,48)=qj(kloop,48)*y(kloop,67)
!
!....  GLYX + hv = 2 CO + H2
!
      qqj(kloop,49)=qj(kloop,49)*y(kloop,68)
!
!....  GLYX + hv = 2 CO + 2 HO2
!
      qqj(kloop,50)=qj(kloop,50)*y(kloop,68)
!
!....  GLYX + hv = CH2O + CO
!
      qqj(kloop,51)=qj(kloop,51)*y(kloop,68)
!
!....  MGLY + hv = CO + HO2 + MCO3
!
      qqj(kloop,52)=qj(kloop,52)*y(kloop,86)
!
!....  MVK + hv = CO + PRPE
!
      qqj(kloop,53)=qj(kloop,53)*y(kloop,89)
!
!....  MVK + hv = MAO3 + MO2
!
      qqj(kloop,54)=qj(kloop,54)*y(kloop,89)
!
!....  MACR + hv =  0.20 CH2O + CO +  1.80 HO2 +  0.20 MCO3 +  0.80 MGLY
!
      qqj(kloop,55)=qj(kloop,55)*y(kloop,79)
!
!....  HAC + hv = CH2O + HO2 + MCO3
!
      qqj(kloop,56)=qj(kloop,56)*y(kloop,69)
!
!....  INPN + hv = HO2 + NO2 + OH + RCHO
!
      qqj(kloop,57)=qj(kloop,57)*y(kloop,74)
!
!....  PRPN + hv = HO2 + NO2 + OH + RCHO
!
      qqj(kloop,58)=qj(kloop,58)*y(kloop,97)
!
!....  ETP + hv = ALD2 + HO2 + OH
!
      qqj(kloop,59)=qj(kloop,59)*y(kloop,66)
!
!....  RA3P + hv = HO2 + OH + RCHO
!
      qqj(kloop,60)=qj(kloop,60)*y(kloop,102)
!
!....  RB3P + hv = HO2 + OH + RCHO
!
      qqj(kloop,61)=qj(kloop,61)*y(kloop,103)
!
!....  R4P + hv = HO2 + OH + RCHO
!
      qqj(kloop,62)=qj(kloop,62)*y(kloop,101)
!
!....  PP + hv = HO2 + OH + RCHO
!
      qqj(kloop,63)=qj(kloop,63)*y(kloop,93)
!
!....  RP + hv = ALD2 + HO2 + OH
!
      qqj(kloop,64)=qj(kloop,64)*y(kloop,111)
!
!....  RIP + hv =  0.69 CH2O +  0.86 HO2 +  0.13 IALD +  0.29 MACR +  0.40 MVK + OH +  0.14 RIO1
!
      qqj(kloop,65)=qj(kloop,65)*y(kloop,109)
!
!....  IAP + hv =  0.67 CO +  0.26 GLYC +  0.19 H2 +  0.36 HAC + HO2 +  0.58 MGLY + OH
!
      qqj(kloop,66)=qj(kloop,66)*y(kloop,72)
!
!....  ISNP + hv = HO2 + NO2 + OH + RCHO
!
      qqj(kloop,67)=qj(kloop,67)*y(kloop,76)
!
!....  VRP + hv =  0.28 CH2O +  0.72 GLYC +  0.28 HO2 +  0.72 MCO3 +  0.28 MGLY + OH
!
      qqj(kloop,68)=qj(kloop,68)*y(kloop,113)
!
!....  MRP + hv =  0.17 CH2O +  0.83 CO +  0.83 HAC + HO2 +  0.17 MGLY + OH
!
      qqj(kloop,69)=qj(kloop,69)*y(kloop,88)
!
!....  MAOP + hv = HO2 + OH + RCHO
!
      qqj(kloop,70)=qj(kloop,70)*y(kloop,82)
!
!....  R4N2 + hv =  0.05 A3O2 +  0.32 ACET +  0.32 ALD2 +  0.18 B3O2 +  0.32 ETO2 +  0.27 HO2 +  0.19 MEK +  0.18 MO2 + NO2 +  0.13 RCHO
!
      qqj(kloop,71)=qj(kloop,71)*y(kloop,99)
!
!....  MAP + hv = MO2 + OH
!
      qqj(kloop,72)=qj(kloop,72)*y(kloop,83)
!
!....  OCSg + hv = CO
!
      qqj(kloop,73)=qj(kloop,73)*y(kloop,114)
      enddo
!
!.... End of kloop
!
      return
      end
