!=======================================================================
!
! $Id: $
!
! FILE
!   setkin_smv2par.h
!   12 JUN 02 - PSC
!
! DESCRIPTION
!   This include file sets symbolic constants for the photochemical
!   mechanism.
!
!  Chemistry input file:    GeosCCM_Combo_2015mechanism.txt
!  Reaction dictionary:     GMI_reactions_JPL15.db
!  Setkin files generated:  Tue Mar 30 16:58:02 2021
!
!========1=========2=========3=========4=========5=========6=========7==

      integer &
     &  SK_IGAS &
     & ,SK_IPHOT &
     & ,SK_ITHERM &
     & ,SK_NACT

      parameter (SK_IGAS   = 129)
      parameter (SK_IPHOT  =  76)
      parameter (SK_ITHERM = 333)
      parameter (SK_NACT   = 125)

!                                  --^--

