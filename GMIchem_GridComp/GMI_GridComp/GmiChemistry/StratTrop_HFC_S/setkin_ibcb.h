!=======================================================================
!
! $Id: setkin_ibcb.h $
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
!  Input mechanism:        GeosCCM_Combo_2015mechanism.txt
!  Reaction dictionary:    GMI_reactions_JPL15.db
!  Setkin files generated: Tue Mar 30 16:58:02 2021
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
!                C2H6
!      ibcb(11) = 1
!                C3H8
!      ibcb(12) = 1
!                CCl4
!      ibcb(13) = 1
!                CF2Br2
!      ibcb(14) = 1
!                CFC12
!      ibcb(15) = 1
!                CF2ClBr
!      ibcb(16) = 1
!                CF3Br
!      ibcb(17) = 1
!                CFC113
!      ibcb(18) = 1
!                CFC114
!      ibcb(19) = 1
!                CFC115
!      ibcb(20) = 1
!                CFC11
!      ibcb(21) = 1
!                CH3Br
!      ibcb(24) = 1
!                CH3CCl3
!      ibcb(25) = 1
!                CH3Cl
!      ibcb(26) = 1
!                CH4
!      ibcb(27) = 1
!                CO
!      ibcb(34) = 1
!                H2402
!      ibcb(41) = 1
!                H2
!      ibcb(42) = 1
!                H2O
!      ibcb(44) = 1
!                HCFC141b
!      ibcb(47) = 1
!                HCFC142b
!      ibcb(48) = 1
!                HCFC22
!      ibcb(49) = 1
!                HCl
!      ibcb(50) = 1
!                ISOP
!      ibcb(72) = 1
!                N2O
!      ibcb(89) = 1
!                                  --^--
