!=======================================================================
!
! $Id$
!
! CODE DEVELOPER
!   Peter Connell, LLNL
!   connell2@llnl.gov
!
! FILE
!   setkin_ibcb.h
!
! DESCRIPTION
!   This include file contains information about treatment of surface
!   boundary conditions.
!
!  Chemistry input file:    08/2018
!  Reaction dictionary:     GMI_Combo_rxns_119species_SO2_JPL15_OCS.db
!  Setkin files generated:  Tue Sep 18 18:25:41 2018
!
!=======================================================================
!
!.... Set default boundary condition types
!
!.... Type 1 means fixed concentration
!.... surface boundary condition, Type 2 means
!.... fixed flux surface boundary condition
!
      ibcb(:)           = 2

      ibcb(NACT+1:IGAS) = 1
!
!.... Reset boundary condition type for special cases
!                CH4
!      ibcb(2) = 1
!                CO
!      ibcb(3) = 1
!                H2
!      ibcb(5) = 1
!                H2O
!      ibcb(10) = 1
!                N2O
!      ibcb(17) = 1
!                HCl
!      ibcb(37) = 1
!                CH3Br
!      ibcb(40) = 1
!                CH3Cl
!      ibcb(41) = 1
!                CH3CCl3
!      ibcb(42) = 1
!                CCl4
!      ibcb(43) = 1
!                CFC11
!      ibcb(44) = 1
!                CFC12
!      ibcb(45) = 1
!                CFC113
!      ibcb(46) = 1
!                CFC114
!      ibcb(47) = 1
!                CFC115
!      ibcb(48) = 1
!                HCFC22
!      ibcb(49) = 1
!                HCFC141b
!      ibcb(50) = 1
!                HCFC142b
!      ibcb(51) = 1
!                CF2Br2
!      ibcb(52) = 1
!                CF2ClBr
!      ibcb(53) = 1
!                CF3Br
!      ibcb(54) = 1
!                H2402
!      ibcb(55) = 1
!                C2H6
!      ibcb(62) = 1
!                C3H8
!      ibcb(63) = 1
!                ISOP
!      ibcb(77) = 1
!                                  --^--
