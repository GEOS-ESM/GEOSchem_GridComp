
!=============================================================================
!
! $Id$
!
! FILE
!   sulfchem.h
!
! CODE DEVELOPER
!   Originated from GRANTOUR programmed by John Walton.
!   Modified by Xiaohong Liu.
!
! DESCRIPTION
!   Do not use implicit statement.
!   Default :
!     implicit integer(i-n)
!     implicit real*8 (a-h,o-z)
!
!=============================================================================

!     ---------------
!     Integer values.
!     ---------------

      integer ::  &
     &  i, ic, j, k, m, n, iw, nc, ijk,  &
     &  ndone, ndt, nr1, numwk, num1


!     -------------
!     Real values.
!     -------------

      real*8  ::  &
     &  tcl,                       & ! accumulated time for cloudy mixing
     &  hpp, rhpp,  &
     &  dmol0, temp0, pres0,  &
     &  dtcyc,  &
     &  conv, convdust, convsslt,  &
     &  ftd, fpi, ftdpi,  &
     &  rhop,  &
     &  alphtrm, dmoltrm,  &
     &  prh2o2, prfso2, prndms, prnso2,  &
     &  dsh2o2, dsfso2, dsndms, dsnso2,  &
     &  zh2o2,  zfso2,  znso2,  zndms,  &
     &  zsfso2, zsndms, zsnso2,  &
     &  za1, za2, za3, za4, za5, za6, za7, za8, za9


!     -----------
!     Parameters.
!     -----------

      integer, parameter ::  &
     &  NOCYC = 5                ! number of sub-cycles in chem60 dtchem


      real*8, parameter ::  &
     &  ALONUM  =   2.68719d19,    & ! Loschmidt's Number (molecules/cm^3)
     &  FRDUST  =   0.05,          & ! fraction of dust     molecules available
                                 !   to dissolve SO4
     &  FRSSLT  =   0.005,         & ! fraction of sea salt molecules available
                                 !   to dissolve SO4
     &  GC      =   8.317d7,       & ! gas constant (ergs/K/mole or cm3*dynes/cm^3)
     &  PI      =   3.141592653589793d0,  &
     &  WTMAIR  =  28.88,          & ! molecular weight of air              (g/mol)
     &  WTMDUST = 100.09,          & ! molecular weight of dust     (CaCO3) (g/mol)
     &  WTMSO2  =  64.06,          & ! molecular weight of SO2              (g/mol)
     &  WTMSO4  =  98.08,          & ! molecular weight of H2SO4            (g/mol)
     &  WTMSSLT =  58.45         ! molecular weight of sea salt (NaCl)  (g/mol)


!     ----------------
!     Data statements.
!     ----------------

      real*8 :: ppr(5)           ! pressure (bar) for 2-d qj H2O2 rates

      data ppr / 0.0, 0.1, 0.4, 0.7, 1.0 /


      real*8 :: ffra_so4(3), nfra_so4(3)  ! relative vol. fraction of 5 bin SO4

      data ffra_so4                & ! from the fossil  (contin.) SO4 size distr.
     &  / 0.125, 0.825, 0.050 /

      data nfra_so4                & ! from the natural (marine)  SO4 size distr.
     &  / 0.100, 0.888, 0.012 /


!     --------
!     Commons.
!     --------

      integer :: nrss

      common / seasalt1 / nrss(5)


      real*8  :: qjh2o2_2d

      common / qjrates / qjh2o2_2d(91,5,12)


      real*8  :: rdryss, rss, dndrtss

      common / seasalt2 / rdryss(2001,5), rss(2001,5), dndrtss(2001,5)

