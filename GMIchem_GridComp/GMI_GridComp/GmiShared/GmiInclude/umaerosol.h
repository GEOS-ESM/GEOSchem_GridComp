
!=============================================================================
!
! $Id$
!
! FILE
!   umaerosol.h
!
! CODE DEVELOPER
!   Xiaohong Liu, University of Michigan
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

      integer :: ng, nso4, nmomso4, nnon, nmomnon, naer
      integer :: iso4, iso4n, iso4m, inon
      integer :: irg, nrg

!     -----------
!     Parameters.
!     -----------

      parameter (ng=3,nso4=2,nmomso4=2,nnon=13,nmomnon=1)
      parameter (naer=nso4*nmomso4)
      parameter (nrg=21)

      real*8, parameter ::  &
     &        r1       = 1.0d0,  &
     &        r2       = 2.0d0,  &
     &        so4mol   = 98.08d-3,  &
     &        avno     = 6.02252d23,  &
     &        so4min   = so4mol/avno,  &
     &        rhoso4   = 1313.d0,  &
     &        r1td     = 1.d0/3.d0,  &
     &        r4td     = 4.d0/3.d0,  &
     &        r3q      = 0.75d0,  &
     &        r1h      = 0.5d0,  &
     &        r3h      = 1.5d0,  &
     &        r5h      = 2.5d0,  &
     &        r10      = 10.0d0,  &
     &        r3000    = 3000.0d0,  &
     &        epsilo   = 1.0d-20,  &
     &        r2so4min = r2 * so4min,   & ! min mass (kg) of one so4 particle
     &        rso4max  = 5.5d-12,       & ! max mass (kg) of one so4 particle (r=10um)
     &        so4radvmin = 1.0d-9,      & ! min volume mean so4 radius (0.001um)
     &        so4radvmax = 1.0d-5     ! max volume mean so4 radius (10 um)

      real*8, parameter :: fso4crit = 10.0d0

      real*8, parameter :: a11 = -1.8436d0,  &
     &                     a12 = -9.21d0

!     -------------
!     Real values.
!     -------------

      real*8 :: sigmod(nso4),radmerg(nso4),  &
     &          radgsrc(ng),siggsrc(ng),  &
     &          radgso4(ng),siggso4(ng),fracso4(ng),  &
     &          radgnon(ng,nnon),siggnon(ng,nnon),fracnon(ng,nnon),  &
     &          rhonon(nnon),xmolnon(nnon)

      integer :: nrh_flg(nnon)

      real*8 :: crh1(nnon), crh2(nnon), crh3(nnon), crh4(nnon)

      real*8 :: radvolm(ng), radvnon(nnon), pmsnon(nnon), pmssrc(ng)

      real*8 :: xlnsg, xlnsg2

      real*8 :: rg(nrg), k0_q0(nrg), k1_q3(nrg), k2_q3(nrg)

      real*8, external :: h2osat_f, so4mfrac_f, so4dens_f

!     ----------------
!     Data statements.
!     ----------------

      data sigmod  /1.514d0,1.77623d0/
      data radmerg /5.d-8,5.d-7/

      data radgsrc /0.005d-6,0.035d-6,0.035d-6/           ! so4 near source
      data siggsrc /1.6d0,2.0d0,2.0d0/

      data radgso4  /0.0210d-6,0.0650d-6,0.3515d-6/       ! so4
      data siggso4  /1.514d0 ,1.778d0 ,1.230d0 /
      data fracso4  /0.8994d0,0.1002d0,0.0004d0/

      data rhonon(1)  /1500.d0/                           ! natural oc
      data radgnon(:,1) /0.0774d-6,0.336d-6,0.9577d-6/
      data siggnon(:,1) /1.402d0,  1.383d0, 1.425d0  /
      data fracnon(:,1) /0.9987d0, 1.306d-3,2.830d-4 /

      data rhonon(2)  /1500.d0/                           ! fossil fuel oc
      data radgnon(:,2) /0.005d-6,0.080d-6,2.500d-6/
      data siggnon(:,2) /1.500d0 ,1.700d0 ,1.650d0 /
      data fracnon(:,2) /0.428571d0,0.571428d0,1.d-6/

      data rhonon(3)  /1500.d0/                           ! fossil fuel bc
      data radgnon(:,3) /0.005d-6,0.080d-6,2.500d-6/
      data siggnon(:,3) /1.500d0 ,1.700d0 ,1.650d0 /
      data fracnon(:,3) /0.428571d0,0.571428d0,1.d-6/

      data rhonon(4)  /1500.d0/                           ! biomass oc
      data radgnon(:,4) /0.0774d-6,0.336d-6,0.9577d-6/
      data siggnon(:,4) /1.402d0,  1.383d0, 1.425d0  /
      data fracnon(:,4) /0.9987d0, 1.306d-3,2.830d-4 /

      data rhonon(5)  /1500.d0/                           ! biomass bc
      data radgnon(:,5) /0.0774d-6,0.336d-6,0.9577d-6/
      data siggnon(:,5) /1.402d0,  1.383d0, 1.425d0  /
      data fracnon(:,5) /0.9987d0, 1.306d-3,2.830d-4 /

      data rhonon(6)  /2500.d0/                           ! dst1
      data radgnon(:,6) /0.201d-6, 2*0.d0/
      data siggnon(:,6) /1.0d0   , 2*1.d0/
      data fracnon(:,6) /1.0d0   , 2*0.d0/

      data rhonon(7)  /2500.d0/                           ! dst2
      data radgnon(:,7) /0.900d-6, 2*0.d0/
      data siggnon(:,7) /1.0d0   , 2*1.d0/
      data fracnon(:,7) /1.0d0   , 2*0.d0/

      data rhonon(8)  /2650.d0/                           ! dst3
      data radgnon(:,8) /1.744d-6, 2*0.d0/
      data siggnon(:,8) /1.0d0   , 2*1.d0/
      data fracnon(:,8) /1.0d0   , 2*0.d0/

      data rhonon(9)  /2650.d0/                           ! dst4
      data radgnon(:,9) /4.137d-6,2*0.d0/
      data siggnon(:,9) /1.0d0   ,2*1.d0/
      data fracnon(:,9) /1.0d0   ,2*0.d0/

      data rhonon(10) /2216.d0/                           ! slt1
      data radgnon(:,10) /0.181d-6,2*0.d0/
      data siggnon(:,10) /1.0d0   ,2*1.d0/
      data fracnon(:,10) /1.0d0   ,2*0.d0/

      data rhonon(11) /2216.d0/                           ! slt2
      data radgnon(:,11) /0.851d-6,2*0.d0/
      data siggnon(:,11) /1.0d0   ,2*1.d0/
      data fracnon(:,11) /1.0d0   ,2*0.d0/

      data rhonon(12) /2216.d0/                           ! slt3
      data radgnon(:,12) /1.568d-6,2*0.d0/
      data siggnon(:,12) /1.0d0   ,2*1.d0/
      data fracnon(:,12) /1.0d0   ,2*0.d0/

      data rhonon(13) /2216.d0/                           ! slt4
      data radgnon(:,13) /2.958d-6,2*0.d0/
      data siggnon(:,13) /1.0d0   ,2*1.d0/
      data fracnon(:,13) /1.0d0   ,2*0.d0/

      data xmolnon(1)  /16.800d-3/                        ! noc
      data xmolnon(2)  /16.800d-3/                        ! foc
      data xmolnon(3)  /12.000d-3/                        ! fbc
      data xmolnon(4)  /16.800d-3/                        ! boc
      data xmolnon(5)  /12.000d-3/                        ! bbc
      data xmolnon(6)  /60.000d-3/                        ! dst1
      data xmolnon(7)  /60.000d-3/                        ! dst2
      data xmolnon(8)  /60.000d-3/                        ! dst3
      data xmolnon(9)  /60.000d-3/                        ! dst4
      data xmolnon(10) /58.450d-3/                        ! slt1
      data xmolnon(11) /58.450d-3/                        ! slt2
      data xmolnon(12) /58.450d-3/                        ! slt3
      data xmolnon(13) /58.450d-3/                        ! slt4

      data nrh_flg(:) /5*0, 4*0, 4*1/

      data crh1(:) /  &
     &     0.4736d0, 0.6251d0, 0.6251d0, 0.4736d0, 0.4736d0,  &
     &     0.0000d0, 0.0000d0, 0.0000d0, 0.0000d0,  &
     &     1.1041d0, 1.1041d0, 1.1041d0, 1.1041d0 /

      data crh2(:) /  &
     &     3.115d0, 3.101d0, 3.101d0, 3.115d0, 3.115d0,  &
     &     0.000d0, 0.000d0, 0.000d0, 0.000d0,  &
     &     3.079d0, 3.079d0, 3.079d0, 3.079d0 /

      data crh3(:) /  &
     &     8.622d-14, 6.520d-14, 6.520d-14, 8.622d-14, 8.622d-14,  &
     &     0.000d-00, 0.000d-00, 0.000d-00, 0.000d-00,  &
     &     3.651d-14, 3.651d-14, 3.651d-14, 3.651d-14 /

      data crh4(:) /  &
     &     -1.399d0, -1.404d0, -1.404d0, -1.399d0, -1.399d0,  &
     &      0.000d0,  0.000d0,  0.000d0,  0.000d0,  &
     &     -1.424d0, -1.424d0, -1.424d0, -1.424d0 /


!c    give below-cloud first-order washout loss rate constant
!c    using the data Dana and Hales (1976) A.E., Fig.2 for q=3

!c    unit: um
      data rg /  &
     &    0.001d0,  0.002d0,  0.004d0,  0.006d0,  0.008d0,  &
     &     0.01d0,   0.02d0,   0.04d0,   0.06d0,   0.08d0,  &
     &      0.1d0,    0.2d0,    0.4d0,    0.6d0,    0.8d0,  &
     &      1.0d0,    2.0d0,    4.0d0,    6.0d0,    8.0d0,  &
     &     10.0d0 /

!c    unit: mm^-1
!c    sigma = 1.5
      data k1_q3 /  &
     &     0.04d0,  0.015d0, 0.0056d0, 0.0032d0, 0.0023d0,  &
     &   0.0018d0, 0.0012d0, 0.0013d0, 0.0016d0, 0.0019d0,  &
     &   0.0023d0, 0.0060d0,  0.025d0,  0.053d0,   0.13d0,  &
     &     0.26d0,   0.74d0,   1.25d0,   1.80d0,    2.0d0,  &
     &      2.2d0 /

!c    sigma = 1.75
      data k2_q3 /  &
     &    0.025d0,  0.009d0,  0.004d0, 0.0023d0, 0.0018d0,  &
     &   0.0013d0, 0.0011d0, 0.0018d0, 0.0023d0, 0.0040d0,  &
     &   0.0052d0,   0.02d0,   0.13d0,   0.32d0,   0.50d0,  &
     &     0.65d0,    1.1d0,    1.4d0,   1.80d0,    2.0d0,  &
     &      2.2d0 /

!c    use Fig.4 for q=0,
!c    when rg < 0.1 um, use sigma=2,n=0,
!c    when rg > 0.1 um, use in-between sigma=2, n=0, and sigma=1
      data k0_q0 /  &
     &     0.12d0,  0.042d0,  0.016d0, 0.0093d0, 0.0065d0,  &
     &   0.0048d0, 0.0023d0, 0.0014d0, 0.0016d0, 0.0018d0,  &
     &    0.002d0,  0.004d0,  0.013d0,  0.035d0,   0.06d0,  &
     &      0.1d0,    0.5d0,    0.8d0,    1.2d0,    1.4d0,  &
     &      1.6d0 /

!     --------
!     Commons.
!     --------

