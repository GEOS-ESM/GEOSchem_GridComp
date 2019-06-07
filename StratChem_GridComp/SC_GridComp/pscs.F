       SUBROUTINE pscs(k,levels,latrad,tropPres,press,prslay,kel,
     &                 h2o_clim,m,sasulf,sanat,saice,rmedsts,
     &                 rmednat,rmedice,denssts,hno3,hno3c,h2o,h2oc,gcSC)
! --------------------------------------------------------------------------------
! This code predicts surface area densities of Type 1 and Type 2 PSCs and
! sediments the HNO3 and H2O condensed into PSCs.  The code was written
! for the Goddard Space Flight Center three dimensional chemistry and
! transport model (GSFC 3D CTM).  It was modified from a PSC
! parameterization written for the Global Modeling Initiative 3D CTM,
! originally named CONDENSE. This is FORTRAN 90 code (although the FORTRAN
! 77 format is used - 5 label columns, a continuation column, 7-72 code
! This version is modified from original to not use dehyd and background
! H2O, but to instead use transported gas phase and condensed H2O.
!
! Author:
!
!    David B. Considine 
!    Earth System Science Interdisciplinary Center
!    and Department of Meteorology
!    University of Maryland
!    College Park, MD 20742
!    email: dbc@welkin.gsfc.nasa.gov
!    phone: (301) 405-5282

! Subroutine Arguments:

!   hno3   - hno3 gas into pscs (mixing ratio units), now num.den.
!   h2o    - h2o gas into pscs (mixing ratio units), now num.den.
!   hno3c  - hno3 condensed into pscs (mixing ratio units), now num.den.
!   h2oc   - condensed phase h2o (mixing ratio), now num.den.
!   sasulf - liquid binary sulfate aerosol surface area density (input)
!   sanat  - nitric acid trihydrate sad (cm-1 units)
!   saice  - ice PSC sad (cm-1 units)
!   stssad - sts psc sad (cm-1 units)

! 12 Dec 2011: In GEOS-5, H2O microphysics is performed by the MOIST grid
!              component. Do not allow updates to H2O or H2O condensed.
! --------------------------------------------------------------------------------

      USE SC_GridCompMod
      IMPLICIT NONE

      TYPE(SC_GridComp), INTENT(INOUT) :: gcSC  ! SC Grid Component

      INTEGER, INTENT(IN) :: k,levels
      REAL, INTENT(IN) :: h2o_clim(levels),kel(levels)
      REAL, INTENT(IN) :: m(levels),latrad,tropPres,press(levels)
      REAL, INTENT(IN) :: prslay(0:levels)
      REAL, INTENT(INOUT) :: sasulf(levels)

      INTEGER, PARAMETER :: DBL = KIND(0.00D+00)
      REAL(KIND=DBL), INTENT(INOUT) :: hno3,hno3c,h2o,h2oc,saice,sanat

      REAL, INTENT(INOUT) :: rmedsts(levels),rmednat(levels)
      REAL, INTENT(INOUT) :: rmedice(levels),denssts(levels)

      REAL, PARAMETER :: PI = 3.14159265
      REAL, PARAMETER :: s = 75.00
      REAL :: stssad, h2so4gas,pscPmax,r,coLat
      REAL(KIND=DBL) :: hno3mr,hno3cmr,h2omr,h2ocmr
      INTEGER :: sign

! Zero out arrays which don't carry info into subroutine:
! rmedice - median radius of ICE aerosols (cm)
! rmednat - median radius of NAT aerosols (cm)
! rmedsts - median radius of STS aerosols (cm)
! -------------------------------------------------------
            rmedice(k) = 0.0
            rmednat(k) = 0.0
            rmedsts(k) = 0.0
            denssts(k) = 0.0
                 saice = 0.0
                 sanat = 0.0

! --------------------------
! Choose high pressure limit
! --------------------------

! Select tropopause pressure if gcSC%PSCpmax is less than zero.
! -------------------------------------------------------------
      IF(gcSC%PSCpmax < 0.00) THEN
       pscPmax = tropPres

      ELSE

! Select the lesser of the tropopause pressure and
! the quasi-tropopause-following function, r, below.
! --------------------------------------------------
       IF(latrad >= 0.00) THEN
        sign = 1
       ELSE
        sign = -1
       END IF

       IF(ABS(latrad) <= PI/6.00) THEN
        coLat = 0.00
       ELSE
        coLat = latrad+sign*2.00*(sign*latrad-PI/6.00)-sign*PI/6.00
       END IF

        r = s+(0.50*(gcSC%PSCpmax-s))*(1.00-COS(coLat))
        pscPmax = MIN(r,tropPres)

      END IF

! Range check
! -----------
      IF(kel(k) <= gcSC%PSCtmax .AND. ABS(latrad) >= gcSC%PSClatlim .AND.
     &  (press(k) <= pscPmax .AND. press(k) >= gcSC%PSCpmin) ) THEN

              h2omr=h2o/m(k)
              h2ocmr=h2oc/m(k)
              hno3mr=hno3/m(k)
              hno3cmr=hno3c/m(k)

              CALL sadice(kel(k),press(k),PI,h2omr,h2ocmr,saice,
     &             rmedice(k),gcSC%constantnice,gcSC%satratice,
     &             gcSC%nice,gcSC%densice,gcSC%sigice,gcSC%rice,m(k))

! If no ice forms (h2ocond = 0) then check to see if STS or
! NAT will form at this location. Note that sadsts is being
! passed ambient h2o ...
! WARNING: ... CALCSULF IS GIVEN H2O_CLIM, NOT AMBIENT H2O!
! ---------------------------------------------------------
               IF(gcSC%calcsts .AND. h2ocmr == 0.00) THEN

! Set fixed species h2o_clim to mixing ratio from h2o
! ---------------------------------------------------
                 stssad = 0.0

		 CALL sadsts(kel(k),press(k),PI,
     &  	      h2omr,hno3mr,hno3cmr,stssad,rmedsts(k),denssts(k),
     &  	      gcSC%constantnsts,gcSC%nsts,gcSC%sigsts,gcSC%rsts,sasulf(k),
     &  	      h2so4gas,gcSC%rlbs,gcSC%siglbs,gcSC%nlbs,gcSC%constantnlbs,
     &  	      h2o_clim(k),m(k))

               ELSE

! In the nat/ice calculation (gcSC%calcsts = false), this code is set up so that
! sadnat will calculate the amount of HNO3 removed from gas phase even if
! ice has formed at this location.  Subroutine sets NAT surface area density
! to zero if ice exists at gridpoint.  The basic assumption begin made here
! is that in the case of ice formation, a NAT phase cocondenses with the ice.
! --------------------------------------------------------------------------
                 CALL sadnat(kel(k),press(k),PI,
     &                h2omr,h2ocmr,hno3mr,hno3cmr,sanat,
     &                rmednat(k),gcSC%constantnnat,gcSC%satratnat,gcSC%nnat,
     &                gcSC%densnat,gcSC%signat,gcSC%rnat,m(k))

               END IF

! Update the HNO3 number densities only. In GEOS-5, the MOIST
! grid component is responsible for water ice microphysics.
! -----------------------------------------------------------
               hno3=hno3mr*m(k)
               hno3c=hno3cmr*m(k)

      ELSE

! Temperature and pressure ranges are not correct. In this case, we
! want to make sure that any values that on a previous step might
! have been assigned a value and that might be used elsewhere
! are now set to 0.
! -----------------------------------------------------------------
               hno3=hno3+hno3c
               hno3c=0.00

      END IF

      RETURN
      END SUBROUTINE pscs
