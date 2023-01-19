!=======================================================================
!
! $Id: $
!
! FILE
!   setkin_par.h
!   12 JUN 02 - PSC
!
! DESCRIPTION
!   This include file sets symbolic constants for the photochemical
!   mechanism.
!
!  Input mechanism:        GeosCCM_Combo_2015mechanism.txt
!  Reaction dictionary:    GMI_reactions_JPL15.db
!  Setkin files generated: Tue Mar 30 16:58:02 2021
!
!========1=========2=========3=========4=========5=========6=========7==


      integer &
     &  NSP &
     & ,NUM_J &
     & ,NUM_K

      parameter (NSP   = 130)
      parameter (NUM_J =  76)
      parameter (NUM_K = 333)

!.... Species affiliations and assignments

      integer &
     &  NACT &
     & ,NCHEM &
     & ,NCONST &
     & ,NFAM &
     & ,NMF &
     & ,NSS

      parameter (NACT   = 125)
      parameter (NCHEM  = 129)
      parameter (NCONST =   4)
      parameter (NFAM   =   2)
      parameter (NMF    = 126)
      parameter (NSS    =   0)

!.... Ancillary input tables

      integer &
     &  nSAD, nSADdust, nSADaer

      parameter (nSAD  =   5)
      parameter (nSADdust =   7)
      parameter (nSADaer =   5)


!.... Species individual identifications

      integer IBBC
      parameter (IBBC      =   0)

      integer IBOC
      parameter (IBOC      =   0)

      integer IDUST1
      parameter (IDUST1    =   0)

      integer IDUST2
      parameter (IDUST2    =   0)

      integer IDUST3
      parameter (IDUST3    =   0)

      integer IDUST4
      parameter (IDUST4    =   0)

      integer IFBC
      parameter (IFBC      =   0)

      integer IFOC
      parameter (IFOC      =   0)

      integer IFSO2
      parameter (IFSO2     =   0)

      integer IFSO4A
      parameter (IFSO4A    =   0)

      integer IFSO4D1
      parameter (IFSO4D1   =   0)

      integer IFSO4D2
      parameter (IFSO4D2   =   0)

      integer IFSO4D3
      parameter (IFSO4D3   =   0)

      integer IFSO4D4
      parameter (IFSO4D4   =   0)

      integer IFSO4N1
      parameter (IFSO4N1   =   0)

      integer IFSO4N2
      parameter (IFSO4N2   =   0)

      integer IFSO4N3
      parameter (IFSO4N3   =   0)

      integer IFSO4N4
      parameter (IFSO4N4   =   0)

      integer IFSO4S1
      parameter (IFSO4S1   =   0)

      integer IFSO4S2
      parameter (IFSO4S2   =   0)

      integer IFSO4S3
      parameter (IFSO4S3   =   0)

      integer IFSO4S4
      parameter (IFSO4S4   =   0)

      integer INDMS
      parameter (INDMS     =   0)

      integer INOC
      parameter (INOC      =   0)

      integer INSO2
      parameter (INSO2     =   0)

      integer INSO4A
      parameter (INSO4A    =   0)

      integer INSO4D1
      parameter (INSO4D1   =   0)

      integer INSO4D2
      parameter (INSO4D2   =   0)

      integer INSO4D3
      parameter (INSO4D3   =   0)

      integer INSO4D4
      parameter (INSO4D4   =   0)

      integer INSO4N1
      parameter (INSO4N1   =   0)

      integer INSO4N2
      parameter (INSO4N2   =   0)

      integer INSO4N3
      parameter (INSO4N3   =   0)

      integer INSO4N4
      parameter (INSO4N4   =   0)

      integer INSO4S1
      parameter (INSO4S1   =   0)

      integer INSO4S2
      parameter (INSO4S2   =   0)

      integer INSO4S3
      parameter (INSO4S3   =   0)

      integer INSO4S4
      parameter (INSO4S4   =   0)

      integer ISSLT1
      parameter (ISSLT1    =   0)

      integer ISSLT2
      parameter (ISSLT2    =   0)

      integer ISSLT3
      parameter (ISSLT3    =   0)

      integer ISSLT4
      parameter (ISSLT4    =   0)

      integer IALD2
      parameter (IALD2     =   3)

      integer IALK4
      parameter (IALK4     =   4)

      integer IBRCL
      parameter (IBRCL     =   8)

      integer IBRONO2
      parameter (IBRONO2   =  10)

      integer IC2H6
      parameter (IC2H6     =  11)

      integer IC3H8
      parameter (IC3H8     =  12)

      integer ICH2O
      parameter (ICH2O     =  23)

      integer ICH4
      parameter (ICH4      =  27)

      integer ICL2
      parameter (ICL2      =  29)

      integer ICL2O2
      parameter (ICL2O2    =  30)

      integer ICLO
      parameter (ICLO      =  32)

      integer ICLONO2
      parameter (ICLONO2   =  33)

      integer ICO
      parameter (ICO       =  34)

      integer IH2
      parameter (IH2       =  42)

      integer IH2O2
      parameter (IH2O2     =  43)

      integer IH2O
      parameter (IH2O      =  44)

      integer IHCL
      parameter (IHCL      =  50)

      integer IHNO2
      parameter (IHNO2     =  59)

      integer IHNO3
      parameter (IHNO3     =  60)

      integer IHNO4
      parameter (IHNO4     =  61)

      integer IHO2
      parameter (IHO2      =  62)

      integer IHOBR
      parameter (IHOBR     =  63)

      integer IHOCL
      parameter (IHOCL     =  64)

      integer IC5H8
      parameter (IC5H8     =  72)

      integer IMACR
      parameter (IMACR     =  74)

      integer IMCO3
      parameter (IMCO3     =  79)

      integer IMEK
      parameter (IMEK      =  80)

      integer IMP
      parameter (IMP       =  84)

      integer IMVK
      parameter (IMVK      =  87)

      integer IN2O5
      parameter (IN2O5     =  88)

      integer IN2O
      parameter (IN2O      =  89)

      integer IAN
      parameter (IAN       =  90)

      integer INO2
      parameter (INO2      =  91)

      integer INO3
      parameter (INO3      =  92)

      integer INO
      parameter (INO       =  93)

      integer IO3
      parameter (IO3       =  95)

      integer IOCLO
      parameter (IOCLO     =  96)

      integer IOH
      parameter (IOH       =  98)

      integer IPMN
      parameter (IPMN      = 101)

      integer IC3H6
      parameter (IC3H6     = 106)

      integer IR4N2
      parameter (IR4N2     = 109)

      integer IRCHO
      parameter (IRCHO     = 114)

      integer IACET
      parameter (IACET     = 126)

      integer INITROGEN
      parameter (INITROGEN = 127)

      integer IOXYGEN
      parameter (IOXYGEN   = 128)

      integer IHNO3COND
      parameter (IHNO3COND = 130)

      integer IMGAS
      parameter (IMGAS     = 129)

!.... Additional model-required parameters

      integer IN2
      parameter (IN2      = 127)

      integer IO2
      parameter (IO2      = 128)

      integer IDEHYD
      parameter (IDEHYD   =   0)

      integer ISYNOZ
      parameter (ISYNOZ   =   0)

      integer JBRONO21
      parameter (JBRONO21 =   0)

      integer KHNO31
      parameter (KHNO31   =   0)


!                                  --^--

