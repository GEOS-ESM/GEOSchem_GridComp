!=======================================================================
!
! $Id$
!
! FILE
!   setkin_par.h
!   12 JUN 02 - PSC
!
! DESCRIPTION
!   This include file sets symbolic constants for the photochemical
!   mechanism.
!
!  Chemistry input file:    08/2018
!  Reaction dictionary:     GMI_Combo_rxns_119species_SO2_JPL15_OCS.db
!  Setkin files generated:  Tue Sep 18 18:25:41 2018
!
!========1=========2=========3=========4=========5=========6=========7==


      integer &
     &  NSP &
     & ,NUM_J &
     & ,NUM_K

      parameter (NSP   = 119)
      parameter (NUM_J =  73)
      parameter (NUM_K = 312)

!.... Species affiliations and assignments

      integer &
     &  NACT &
     & ,NCHEM &
     & ,NCONST &
     & ,NDYN &
     & ,NFAM &
     & ,NMF &
     & ,NSS

      parameter (NACT   = 114) !Used with ibcb in setkin_ibcb.h
      parameter (NCHEM  = 118) !Apparently not used
      parameter (NCONST =   4) !Apparently not used
      parameter (NDYN   = 115) !Sizes ldynvar(NDYN) in setkin_lchem.h
      parameter (NFAM   =   2) !Probably replaced by NUMGRP in setkin_group.h
      parameter (NMF    = 115) !Sizes specarr in setkin_kcalc.F90
      parameter (NSS    =   0)

!.... Ancillary input tables

      integer &
     &  NSAD, NSADdust, NSADaer

      parameter (NSAD  =   5)
      parameter (NSADdust =   7)
      parameter (NSADaer =   5)


!.... Species individual identifications

! -----------------------------------------------------------------
! WARNING: Integers IBBC through ISSLT4 are erroneously set to one
! to enable bounds checking at the compile step. Several options, 
! among which are (but not limited to!) chem_opt, do_aerocom, and
! do_dust_emiss, control whether these integers are used in array
! references.
! -----------------------------------------------------------------

      integer IBBC
      parameter (IBBC      =   1)

      integer IBOC
      parameter (IBOC      =   1)

      integer IDUST1
      parameter (IDUST1    =   1)

      integer IDUST2
      parameter (IDUST2    =   1)

      integer IDUST3
      parameter (IDUST3    =   1)

      integer IDUST4
      parameter (IDUST4    =   1)

      integer IFBC
      parameter (IFBC      =   1)

      integer IFOC
      parameter (IFOC      =   1)

      integer IFSO2
      parameter (IFSO2     =   1)

      integer IFSO4A
      parameter (IFSO4A    =   1)

      integer IFSO4D1
      parameter (IFSO4D1   =   1)

      integer IFSO4D2
      parameter (IFSO4D2   =   1)

      integer IFSO4D3
      parameter (IFSO4D3   =   1)

      integer IFSO4D4
      parameter (IFSO4D4   =   1)

      integer IFSO4N1
      parameter (IFSO4N1   =   1)

      integer IFSO4N2
      parameter (IFSO4N2   =   1)

      integer IFSO4N3
      parameter (IFSO4N3   =   1)

      integer IFSO4N4
      parameter (IFSO4N4   =   1)

      integer IFSO4S1
      parameter (IFSO4S1   =   1)

      integer IFSO4S2
      parameter (IFSO4S2   =   1)

      integer IFSO4S3
      parameter (IFSO4S3   =   1)

      integer IFSO4S4
      parameter (IFSO4S4   =   1)

      integer INDMS
      parameter (INDMS     =   1)

      integer INOC
      parameter (INOC      =   1)

      integer INSO2
      parameter (INSO2     =   1)

      integer INSO4A
      parameter (INSO4A    =   1)

      integer INSO4D1
      parameter (INSO4D1   =   1)

      integer INSO4D2
      parameter (INSO4D2   =   1)

      integer INSO4D3
      parameter (INSO4D3   =   1)

      integer INSO4D4
      parameter (INSO4D4   =   1)

      integer INSO4N1
      parameter (INSO4N1   =   1)

      integer INSO4N2
      parameter (INSO4N2   =   1)

      integer INSO4N3
      parameter (INSO4N3   =   1)

      integer INSO4N4
      parameter (INSO4N4   =   1)

      integer INSO4S1
      parameter (INSO4S1   =   1)

      integer INSO4S2
      parameter (INSO4S2   =   1)

      integer INSO4S3
      parameter (INSO4S3   =   1)

      integer INSO4S4
      parameter (INSO4S4   =   1)

      integer ISSLT1
      parameter (ISSLT1    =   1)

      integer ISSLT2
      parameter (ISSLT2    =   1)

      integer ISSLT3
      parameter (ISSLT3    =   1)

      integer ISSLT4
      parameter (ISSLT4    =   1)

      integer ICH2O
      parameter (ICH2O     =   1)

      integer ICH4
      parameter (ICH4      =   2)

      integer ICO
      parameter (ICO       =   3)

      integer IH2
      parameter (IH2       =   5)

      integer IHNO2
      parameter (IHNO2     =   7)

      integer IHNO3
      parameter (IHNO3     =   8)

      integer IHNO4
      parameter (IHNO4     =   9)

      integer IH2O
      parameter (IH2O      =  10)

      integer IHO2
      parameter (IHO2      =  11)

      integer IH2O2
      parameter (IH2O2     =  12)

      integer IMP
      parameter (IMP       =  15)

      integer IAN
      parameter (IAN       =  16)

      integer IN2O
      parameter (IN2O      =  17)

      integer INO
      parameter (INO       =  18)

      integer INO2
      parameter (INO2      =  19)

      integer INO3
      parameter (INO3      =  20)

      integer IN2O5
      parameter (IN2O5     =  21)

      integer IO3
      parameter (IO3       =  24)

      integer IOH
      parameter (IOH       =  25)

      integer IBRCL
      parameter (IBRCL     =  27)

      integer IBRONO2
      parameter (IBRONO2   =  29)

      integer IHOBR
      parameter (IHOBR     =  31)

      integer ICL2
      parameter (ICL2      =  33)

      integer ICLO
      parameter (ICLO      =  34)

      integer ICL2O2
      parameter (ICL2O2    =  35)

      integer ICLONO2
      parameter (ICLONO2   =  36)

      integer IHCL
      parameter (IHCL      =  37)

      integer IHOCL
      parameter (IHOCL     =  38)

      integer IOCLO
      parameter (IOCLO     =  39)

      integer IALD2
      parameter (IALD2     =  58)

      integer IALK4
      parameter (IALK4     =  59)

      integer IC2H6
      parameter (IC2H6     =  62)

      integer IC3H8
      parameter (IC3H8     =  63)

      integer IC5H8
      parameter (IC5H8     =  77)

      integer IMACR
      parameter (IMACR     =  79)

      integer IMCO3
      parameter (IMCO3     =  84)

      integer IMEK
      parameter (IMEK      =  85)

      integer IMVK
      parameter (IMVK      =  89)

      integer IPMN
      parameter (IPMN      =  91)

      integer IC3H6
      parameter (IC3H6     =  96)

      integer IR4N2
      parameter (IR4N2     =  99)

      integer IRCHO
      parameter (IRCHO     = 104)

      integer IACET
      parameter (IACET     = 115)

      integer INITROGEN
      parameter (INITROGEN = 116)

      integer IOXYGEN
      parameter (IOXYGEN   = 117)

      integer IHNO3COND
      parameter (IHNO3COND = 119)

      integer IMGAS
      parameter (IMGAS     = 118)

!.... Additional model-required parameters

      integer IN2
      parameter (IN2      = 116)

      integer IO2
      parameter (IO2      = 117)

! -----------------------------------------------------------------
! WARNING: In GEOS-5 we have abandoned synthetic ozone as a
! discriminator of stratospheric air.  To enable bounds checking,
! ISYNOZ is erroneously set to one, and do_synoz must always be
! .FALSE. 
! -----------------------------------------------------------------

      integer ISYNOZ
      parameter (ISYNOZ   =   1)

! -----------------------------------------------------------------
! WARNING: In GEOS-5 we no longer transport dehyd and its only 
! array reference in sad_update.F90 has been commented out. 
! -----------------------------------------------------------------

      integer IDEHYD
      parameter (IDEHYD   =   -99)

! -----------------------------------------------------------------
! The following two parameters appear to never be referenced.
! -----------------------------------------------------------------

      integer JBRONO21
      parameter (JBRONO21 =   0)

      integer KHNO31
      parameter (KHNO31   =   0)


!                                  --^--

