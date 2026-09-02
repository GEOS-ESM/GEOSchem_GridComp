      SUBROUTINE meso_phot(k,levels,notify,nocolumn,o3column,o2column,
     &                     szan,aj,m,gcSC,rdistsq)

      USE SC_GridCompMod
      IMPLICIT none
      INTEGER, PARAMETER :: DBL = KIND(0.00D+00)

      TYPE(SC_GridComp), INTENT(IN) :: gcSC     ! Grid Component

      INTEGER, INTENT(IN) :: k,levels,notify
      REAL, INTENT(IN) :: szan,rdistsq
      REAL, INTENT(IN) :: o3column(levels),nocolumn(levels)
      REAL, INTENT(IN) :: o2column(levels),m(levels)

      REAL(KIND=DBL), INTENT(INOUT) :: aj(gcSC%numphoto)

! -------------------------------------------------------------------------
! Compute photolysis rates for:
!		     Number in     aj
!     Reaction       2d model	 Number
!  ----------------- ---------   ------
!  h2o+hv = h+oh	 4	   38
!  no +hv = n+o 	16	   39
!  h2o+hv = h2+o(1d)	25	   40
!  o2 +hv = o+o(1d)	46	   41
!  ch4+hv = 2h2o+co	59	   42

! All data included here are taken from GSFC 2-D model.  Assume branching 
! ratios for these photoreactions is 1.0

! Wavelength data is in Angstrom
! -------------------------------------------------------------------------

      INTEGER, PARAMETER :: nW = 19

      REAL, PARAMETER ::   wavel(nW) = (/1215.7, 1709.5, 1731.5, 1746.5, 1762.0, 1778.0, 1794.0, 1810.0,
     &           1826.5, 1843.5, 1860.5, 1878.0, 1896.0, 1914.0, 1932.5, 1951.5, 1970.5, 1990.0, 2010.0/)

      REAL, PARAMETER ::   sflux(nW) = (/4.006E+11, 1.764E+11, 1.017E+11, 1.302E+11, 1.722E+11, 2.200E+11,
     &                                   2.438E+11, 3.215E+11, 3.691E+11, 3.549E+11, 4.190E+11, 5.555E+11,
     &    	              6.464E+11, 7.342E+11, 7.709E+11, 1.055E+12, 1.189E+12, 1.333E+12, 1.610E+12/)

      REAL, PARAMETER ::   o2xs1(nW) = (/3.000E-21, 1.336E-19, 1.000E-28, 1.000E-28, 1.500E-19, 1.200E-20,
     &                                   6.500E-20, 5.000E-20, 3.000E-20, 1.400E-20, 7.600E-21, 3.500E-21,
     &    	              1.400E-21, 6.000E-22, 5.700E-22, 1.900E-22, 5.400E-23, 8.130E-24, 8.012E-24/)

      REAL, PARAMETER ::   o2xs2(nW) = (/7.000E-21, 3.851E-18, 4.580E-19, 2.740E-19, 0.000E+00, 0.000E+00,
     &    	                         0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     &    	              0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00/)

      REAL, PARAMETER ::   o3xs1(nW) = (/0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 7.297E-19, 7.189E-19, 
     &    	                         7.070E-19, 6.867E-19, 6.560E-19, 6.180E-19, 5.598E-19, 5.183E-19, 
     &    	              4.726E-19, 4.284E-19, 3.852E-19, 3.441E-19, 3.123E-19, 2.909E-19, 2.826E-19/)

      REAL, PARAMETER ::   o3xs2(nW) = (/0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 8.108E-20, 7.987E-20,
     &    	                         7.855E-20, 7.630E-20, 7.289E-20, 6.867E-20, 6.220E-20, 5.759E-20, 
     &    	              5.251E-20, 4.760E-20, 4.281E-20, 3.824E-20, 3.470E-20, 3.233E-20, 3.140E-20/)

      REAL, PARAMETER ::  h2o4xs(nW) = (/1.440E-17, 0.000E+00, 0.000E+00, 2.750E-18, 2.367E-18, 1.719E-18,
     &                                   1.020E-18, 5.041E-19, 2.190E-19, 8.417E-20, 3.166E-20, 1.438E-20,
     &                        6.262E-21, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00/)

      REAL, PARAMETER :: h2o25xs(nW) = (/1.600E-18, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     &    	                         0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     &    	              0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00/)

      REAL, PARAMETER ::    noxs(nW) = (/2.420E-18, 8.681E-19, 4.000E-19, 0.000E+00, 0.000E+00, 0.000E+00,
     &    	                         0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     &    	              0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00/)

      REAL, PARAMETER ::  o246xs(nW) = (/7.000E-21, 3.851E-18, 4.580E-19, 2.740E-19, 0.000E+00, 0.000E+00,
     &    	                         0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     &    	              0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00/)

      REAL, PARAMETER ::   ch4xs(nW) = (/2.000E-17, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     &                                   0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     &                        0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00/)

      INTEGER, PARAMETER :: nC = 6

      REAL, PARAMETER ::   smo2a(nC) = (/1.120E-23, 2.450E-23, 7.190E-23, 3.040E-22, 1.750E-21, 1.110E-20/)
      REAL, PARAMETER ::   smwa1(nC) = (/0.000E+00, 5.120E-02, 1.360E-01, 1.650E-01, 1.410E-01, 4.500E-02/)
      REAL, PARAMETER ::  smnoa1(nC) = (/0.000E+00, 1.320E-18, 6.350E-19, 7.090E-19, 2.180E-19, 4.670E-19/)
      REAL, PARAMETER ::   smwa2(nC) = (/0.000E+00, 5.680E-03, 1.520E-02, 1.830E-02, 1.570E-02, 5.000E-03/)
      REAL, PARAMETER ::  smnoa2(nC) = (/0.000E+00, 4.410E-17, 4.450E-17, 4.500E-17, 2.940E-17, 4.350E-17/)
      REAL, PARAMETER ::   smo2b(nC) = (/1.350E-22, 2.990E-22, 7.330E-22, 3.070E-21, 1.690E-20, 1.660E-19/)
      REAL, PARAMETER ::   smwb1(nC) = (/0.000E+00, 0.000E+00, 1.930E-03, 9.730E-02, 9.750E-02, 3.480E-02/)
      REAL, PARAMETER ::  smnob1(nC) = (/0.000E+00, 0.000E+00, 3.050E-21, 5.760E-19, 2.290E-18, 2.210E-18/)
      REAL, PARAMETER ::   smwb2(nC) = (/0.000E+00, 0.000E+00, 2.140E-04, 1.080E-02, 1.080E-02, 3.860E-03/)
      REAL, PARAMETER ::  smnob2(nC) = (/0.000E+00, 0.000E+00, 3.200E-21, 5.710E-17, 9.090E-17, 6.000E-17/)
      REAL, PARAMETER ::   smo2c(nC) = (/2.970E-22, 5.830E-22, 2.050E-21, 8.190E-21, 4.800E-20, 2.660E-19/)
      REAL, PARAMETER ::   smwc1(nC) = (/4.500E-02, 1.800E-01, 2.250E-01, 2.250E-01, 1.800E-01, 4.500E-02/)
      REAL, PARAMETER ::  smnoc1(nC) = (/1.800E-18, 1.500E-18, 5.010E-19, 7.200E-20, 6.720E-20, 1.490E-21/)
      REAL, PARAMETER ::   smwc2(nC) = (/5.000E-03, 2.000E-02, 2.500E-02, 2.500E-02, 2.000E-02, 5.000E-03/)
      REAL, PARAMETER ::  smnoc2(nC) = (/1.400E-16, 1.520E-16, 7.000E-17, 2.830E-17, 2.730E-17, 6.570E-18/)

! LDO Add Chapman function for high solar zenith angle grazing from 2D model
! --------------------------------------------------------------------------
      REAL :: cosbc,r,rtd,sinc,chida,zgrz,sfaca,s

      REAL, PARAMETER :: d1 = 1.060693
      REAL, PARAMETER :: d2 = 0.55643831
      REAL, PARAMETER :: d3 = 1.0619896
      REAL, PARAMETER :: d4 = 1.7245609
      REAL, PARAMETER :: d5 = 0.56498823
      REAL, PARAMETER :: d6 = 0.06651874

      REAL, PARAMETER :: hbar = 6.79
      REAL, PARAMETER :: zbar = 30.0
      REAL, PARAMETER :: r0   = 6.371E3
      REAL, PARAMETER :: zp   = 60.0

      REAL, PARAMETER :: PI = 3.14159265

! End of parameters for chapman function from 2D model
! ----------------------------------------------------

      REAL :: opt_dep, sfluxsc(nW), sfluxz(nW)
      REAL :: wavel1, afac, raycs
      REAL :: enop, spdis, emission, rkq, denn2, pdis, to3col, colno, colo2
      REAL :: terma, termb, termc, tmno, frno, fracwi
      INTEGER :: l, i, jl

      rtd  = 180.0/PI

! WHS: Scale solar flux for solar cycle.
! --------------------------------------
!     sfluxsc(1:nW) = sflux(1:nW)  !no solar cycle scaling
      sfluxsc(1:nW) = sflux(1:nW)*rdistsq  !no solar cycle scaling

! Calculate the height and wavelength dependent optical depth and solar 
! intensity for each cell. Rayleigh cross-section is in cm-2 (from 
! Brasseur-Solomon, 1986; pp. 107)
! ---------------------------------------------------------------------
      aj(38)=0.0
      aj(39)=0.0
      aj(40)=0.0
      aj(41)=0.0
      aj(42)=0.0

      IF(notify == 1) PRINT *,'MESO_PHOT: CORRECTED photolysis for the mesosphere.'

! LDO Allow solar zenith angles up to the designated value
! --------------------------------------------------------

      Daytime: IF(szan < gcSC%szaCutoff*PI/180.00) THEN
   
! LDO Chapman Function Calculation from 2D model
! ----------------------------------------------
       cosbc = COS(szan)
       cosbc = MAX(-1.00, cosbc)
       cosbc = MIN(cosbc,  1.00)

       r = ACOS(cosbc)
       sinc = SIN(r)
       chida = r*rtd

       IF(chida <= 90.0) THEN
        zgrz = 1000.0
       ELSE
        zgrz = sinc*(zp+r0)-r0
       ENDIF

       IF(chida < 70.0) THEN
        sfaca = 1.0/cosbc
       ELSE IF (zgrz > 0.0) THEN
        r = SQRT(0.50*r0/hbar)
        s = r*ABS(cosbc)

        IF(s <= 8.0) THEN
         s = (d1+d2*s)/(d3+d4*s+s**2)
        ELSE
         s = d5/(d6+s)
        END IF

        r = r*SQRT(PI)
        sfaca = r*s

        IF(chida > 90.0) THEN
         sfaca = 2.0*r*EXP((r0+zbar)*(1.0-sinc)/hbar) - sfaca
        END IF
       END IF

! LDO End of Chapman Function Calculation 
! ---------------------------------------

       DO l = 1, nW

        wavel1 = wavel(l)*1.0E-4
        afac = 3.916+0.074*wavel1+0.05/wavel1
        afac = wavel1**(afac)
        raycs = (4.0E-28)/afac

        opt_dep = o3column(k)*(o3xs1(l)+o3xs2(l))+
     &  	  o2column(k)*(o2xs1(l)+o2xs2(l))+
     &            o2column(k)*(1.0/0.22)*raycs

! LDO now calculating airmass factor using Chapman Function from 2D model
! -------------------------------------------------------------------------
        sfluxz(l) = sfluxsc(l)*EXP(-opt_dep*sfaca)

        aj(38) = aj(38)+sfluxz(l)* h2o4xs(l)
        aj(39) = aj(39)+sfluxz(l)*   noxs(l)
        aj(40) = aj(40)+sfluxz(l)*h2o25xs(l)
        aj(41) = aj(41)+sfluxz(l)* o246xs(l)
        aj(42) = aj(42)+sfluxz(l)*  ch4xs(l)

       END DO

! -------------------------------------------------------------------------
! Special treatment for calculation of aj(39) - NO photorate borrowed 
! from Eric Fleming and Charley Jackman.  Use Murray et al., J. Chem. Phys., 
! 101, 62-73, 1994 delta band oscillator strengths.
!
! Restrict to stratosphere and mesosphere
! -------------------------------------------------------------------------
       Jackman: IF(k >= gcSC%kNOspec) THEN

        enop = aj(39)

! SPDIS is the rate of spontaneous dissociation
! ---------------------------------------------
        spdis = 1.65E+09

! Emission from the C(2PIv=0) state
! ---------------------------------
        emission = 5.1E+07

! RKQ is the quenching rate constant from N2
! ------------------------------------------
        rkq = 1.5E-09

! DENN2 is the number density of N2
! ---------------------------------
        denn2 = 0.79054*m(k)

! PDIS is the rate of pontaneous dissociation
! -------------------------------------------
        pdis = spdis/(emission + spdis + (rkq*denn2))

! LDO now calculating airmass factor using Chapman Function from 2D model
! -----------------------------------------------------------------------
        to3col = o3column(k)*sfaca
        colno = nocolumn(k)*sfaca
        colo2 = o2column(k)*sfaca

        terma = 0.00
        DO i = 1,nC
         tmno = EXP(-smo2a(i)*colo2) *
     &           (smwa1(i)*smnoa1(i)*EXP(-smnoa1(i)*colno) +
     &            smwa2(i)*smnoa2(i)*EXP(-smnoa2(i)*colno))
         terma = terma+tmno
        END DO

        termb = 0.00
        DO i = 1,nC
         tmno = EXP(-smo2b(i)*colo2) *
     &           (smwb1(i)*smnob1(i)*EXP(-smnob1(i)*colno) +
     &            smwb2(i)*smnob2(i)*EXP(-smnob2(i)*colno))
         termb = termb+tmno
        END DO

        termc = 0.00
        DO i = 1,nC
         tmno = EXP(-smo2c(i)*colo2) *
     &           (smwc1(i)*smnoc1(i)*EXP(-smnoc1(i)*colno) +
     &            smwc2(i)*smnoc2(i)*EXP(-smnoc2(i)*colno))
         termc = termc+tmno
        END DO

        DO jl = 8,9
         frno = pdis * EXP(-(o3xs1(jl)+o3xs2(jl))*to3col)
         IF(jl == 8) fracwi = (181.8-181.6)/(181.8-180.2)
         IF(jl == 9) fracwi = (183.1-181.8)/(183.5-181.8)
         enop = enop + fracwi * sfluxsc(jl) * frno * termc
        END DO

        DO jl = 9,10
         frno=pdis * EXP(-(o3xs1(jl)+o3xs2(jl))*to3col)
         IF(jl ==  9) fracwi = (183.5-183.1)/(183.5-181.8)
         IF(jl == 10) fracwi = (184.6-183.5)/(185.2-183.5)
         enop = enop + fracwi * sfluxsc(jl) * frno * termb
        END DO

        DO jl = 13,15
         frno = pdis * EXP(-(o3xs1(jl)+o3xs2(jl))*to3col)
         IF(jl == 13) fracwi = (190.5-190.2)/(190.5-188.7)
         IF(jl == 14) fracwi = 1.0E0
         IF(jl == 15) fracwi = (192.5-192.3)/(194.2-192.3)

! Reduce photorate from 1.45 to 1.00 following the
! practice of Jackman in the 2-D model.  JEN 18-Sep-2006
! ------------------------------------------------------
         enop = enop + 1.00 * fracwi * sfluxsc(jl) * frno * terma

        END DO

        aj(39) = enop

       END IF Jackman

      END IF Daytime

      RETURN
      END SUBROUTINE meso_phot

