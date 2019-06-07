      SUBROUTINE sza(dayOfYear,rlong,rlat,szan,ecc_SC,rdistsq)
!--------------------------------------------------------------------------
! Compute the solar zenith angle from day-of-year, latitude, and longitude.
!
!  dayOfYear  INPUT   REAL   Day-of-year (1.0-365.99)
!  rlong      INPUT   REAL   Longitude in radians
!  rlat       INPUT   REAL   Latitude in radians
!  sza	      OUTPUT  REAL   Solar zenith angle (radians)
!  ecc_SC     INPUT   REAL   Eccentricity of the Earth's orbit
!  rdistsq    OUTPUT  REAL   Square of the ratio of the mean Earth-Sun
!                            distance to the current Earth-Sun distance
!--------------------------------------------------------------------------
      IMPLICIT NONE
      REAL, INTENT(IN) :: dayOfYear, rlong, rlat,ecc_SC
      REAL, INTENT(OUT) :: szan,rdistsq

      REAL, PARAMETER :: PI   = 3.14159265
      REAL, PARAMETER :: DEC1 = 4.8795
      REAL, PARAMETER :: DEC2 = 0.0172027
      REAL, PARAMETER :: DEC3 = 6.23025
      REAL, PARAMETER :: DEC4 = 0.0172019
      REAL, PARAMETER :: DEC5 = 0.03344
      REAL, PARAMETER :: DEC6 = 0.000349
      REAL, PARAMETER :: DEC7 = 0.39782

      INTEGER :: iday
      REAL :: ameanl, ameana, arg, cossza
      REAL :: cosd, decl, eclipl, sind, t
      REAL :: e_new, meana

! Truncate day-of-year
! --------------------
      iday = dayOfYear

! Calculate solar declination
! ---------------------------
      ameanl = DEC1 + DEC2*iday
      ameana = DEC3 + DEC4*iday
      eclipl = ameanl + DEC5*SIN(ameana) + DEC6*SIN(2*ameana)
      sind = DEC7*SIN(eclipl)    
      decl = ASIN(sind) 
      cosd = COS(decl)

! Calculate hour angle
! --------------------
      t = dayOfYear-iday
      arg = rlong + pi*(t+t-1.)

! Zenith angle (radians), with precision check
! --------------------------------------------
      cossza = COS(arg)*COS(rlat)*cosd +SIN(rlat)*sind
      
      cossza = MAX(-1.00, cossza)
      cossza = MIN(cossza,  1.00)

      szan = ACOS(cossza)

! LDO 9/2/14 Add in eccentricity impact on solar flux from 2-D model

      e_new=ecc_SC

      meana = dec3+dec4*iday
      rdistsq = 1.0/(1.0-e_new*COS(meana))**2

      RETURN
      END SUBROUTINE sza
