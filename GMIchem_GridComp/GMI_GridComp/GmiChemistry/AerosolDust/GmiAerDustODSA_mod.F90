       module GmiAerDustODSA_mod

       implicit none

       private
       public   :: Aero_OptDep_SurfArea
       public   :: Dust_OptDep_SurfArea

       CONTAINS

!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Aero_OptDep_SurfArea
!
! !INTERFACE:
!
      subroutine Aero_OptDep_SurfArea(gridBoxHeight, concentration, tropp, pres3c, &
     &    optDepth, eRadius, tArea, Odaer, relativeHumidity, dAersl, wAersl, &
     &    raa_b, qaa_b, do_synoz, isynoz_num, synoz_threshold, AerDust_Effect_opt, &
     &    i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi, num_species, num_AerDust)
!
! !USES:
      use GmiArrayBundlePointer_mod, only : t_GmiArrayBundle
!
      implicit none

#     include "gmi_phys_constants.h"
#     include "setkin_par.h"
#     include "gmi_AerDust_const.h"
!
! !INPUT PARAMETERS:
  integer                , intent(in) :: i1, i2, ju1, j2, k1, k2
  integer                , intent(in) :: ilo, ihi, julo, jhi
  integer                , intent(in) :: num_species
  integer                , intent(in) :: num_AerDust
  integer                , intent(in) :: AerDust_Effect_opt
  logical                , intent(in) :: do_synoz
  integer                , intent(in) :: isynoz_num
  real*8                 , intent(in) :: synoz_threshold
  real*8                 , intent(in) :: gridBoxHeight   (i1:i2,   ju1:j2,  k1:k2)
  real*8                 , intent(in) :: relativeHumidity(i1:i2, ju1:j2, k1:k2)
  REAL*8                 , intent(in) :: raa_b(4, NP_b)
  REAL*8                 , intent(in) :: qaa_b(4, NP_b)
  type (t_GmiArrayBundle), intent(in) :: concentration(num_species)

  REAL*8, INTENT(IN) :: tropp(i1:i2, ju1:j2)         !Tropopause pressure (hPa)
  REAL*8, INTENT(IN) :: pres3c(i1:i2, ju1:j2, k1:k2) !Layer mean pressure (hPa)

  REAL*8 :: optDepth(i1:i2, ju1:j2, k1:k2, num_AerDust)
  REAL*8 :: eRadius (i1:i2, ju1:j2, k1:k2, NSADdust+NSADaer)
  REAL*8 :: tArea   (i1:i2, ju1:j2, k1:k2, NSADdust+NSADaer)
  REAL*8 :: odAer   (i1:i2, ju1:j2, k1:k2, NSADaer*NRH_b)
  REAL*8 :: dAersl  (i1:i2, ju1:j2, k1:k2, 2)
  REAL*8 :: wAersl  (i1:i2, ju1:j2, k1:k2, NSADaer )

!
! !DESCRIPTION:
!   Calculate optical depth and surface area at each timestep
!   to account for the change in relative humidity.
!
!   For the optical depth calculation, this involves carrying the
!   optical depth at each RH as separate aerosols.
!
!   Scaling is sufficient for the surface area calculation.
!
!   This routine is a version of RDAER from GEOS-CHEM.
!   Here we do not read global aerosol concentrations.
!
!   We only compute the optical depth and surface area at tropospheric
!   levels. We use the conditions:
!
!                WHERE (const(:,:,:,isynoz_num) < synoz_threshold))
!
! !LOCAL VARIABLES:
!!    FWET:     Fraction of aerosol from H2O
!!    RW:       Effective radius at RH bins read in from "jv_spec.dat"
!!    REFF:     Effective radius at RH after interpolation
!!    QW:       Q at different RH bins read in from "jv_spec.dat"
!!    FRAC:     Used to interpolate between sizes
!!    SCALEQ:   Change in Q (extinction efficiency)
!!    SCALER:   Change in Radius with RH
!!    SCALEOD:  Change in Optical Depth vs RH
!!    SCALEVOL: Change in Vol vs RH
!!    RH:       Relative Humidity bins

      REAL*8,  SAVE  :: RH(NRH_b)   = (/0d0,0.5d0,0.7d0,0.8d0,0.9d0/)
      REAL*8         :: FWET, REFF, FRAC, SCALEQ, SCALER
      REAL*8         :: RW(NRH_b)
      REAL*8         :: QW(NRH_b)
      REAL*8         :: SCALEOD (i1:i2, ju1:j2, k1:k2, NRH_b     )
      REAL*8         :: SCALEVOL(i1:i2, ju1:j2, k1:k2            )
      REAL*8         :: MSDENS(NSADaer), XTAU
      REAL*8         :: dryArea

!   Index to aerosol types in jv_spec.dat
!   The following are ordered according to the mass densities below
      INTEGER, SAVE  :: IND(NSADaer) = (/22, 29, 36, 43, 50/)

      INTEGER        :: I, J, R, L, N, IRH, IRHN, IU_cross
      LOGICAL        :: do_calc
!
! !REVISION HISTORY:
!   February2005, Jules Kouatchou (Jules.Kouatchou.1@gsfc.nasa.gov)
!     Original code.
!   January 2008, Eric Nielsen 
!     Enable use of troposphere flag when do_synoz=.FALSE.
!
!EOP
!-------------------------------------------------------------------------
!BOC

      MSDENS(1) = 1700.0    !SO4
      MSDENS(2) = 1000.0    !BC
      MSDENS(3) = 1800.0    !OC
      MSDENS(4) = 2200.0    !SS (accum)
      MSDENS(5) = 2200.0    !SS (coarse)

      ! Loop over types of aerosol
      DO N = 1, NSADaer

         ! Zero array
         SCALEOD(:,:,:,:) = 0.0d0

         !==============================================================
         ! Determine aerosol growth rates from the relative
         ! humidity in each box
         !
         ! The optical depth scales with the radius and Q alone
         ! since SCALEDENS cancels as follows
         !
         !    SCALER    = RW / RDRY
         !    SCALEDENS = DENSWET / DENSDRY
         !    SCALEM    = SCALEDENS * SCALER**3
         !    SCALEOD   = (SCALEQ * SCALEM) / (SCALEDENS * SCALER)
         !              = SCALEQ * SCALER**2
         !
         ! Cap aerosol values at 90% relative humidity since
         ! aerosol growth at that point becomes highly nonlinear and
         ! relative humidities above this value essentially mean
         ! there is a cloud in that grid box
         !
         ! Q is the extinction efficiency
         !
         ! Each grid box (I,J,L) will fall into one of the RH bins,
         ! since each grid box will have a different RH value.  So,
         ! for SCALEOD(I,J,L,:), only one of the IRH bins will contain
         ! nonzero data, while the other IRH bins will all be zero.
         !==============================================================

         ! Loop over relative humidity bins
         DO R = 1, NRH_b

            ! Wet radius in "jv_spec.dat"
            RW(R) = raa_b(4,IND(N)+R-1)

            ! Wet frac of aerosol
            FWET  = (RW(R)**3 - RW(1)**3) / RW(R)**3

            ! Extinction efficiency Q for each RH bin
            QW(R) = qaa_b(4,IND(N)+R-1)*FWET + qaa_b(4,IND(N))*(1.d0-FWET)
         ENDDO

         ! Loop over SMVGEAR grid boxes
         DO L = k1, k2
            DO J = ju1, j2
               DO I = i1, i2

                  ! Sort into relative humidity bins
                  IF (      relativeHumidity(I,J,L) <= RH(2) ) THEN
                     IRH = 1
                  ELSE IF ( relativeHumidity(I,J,L) <= RH(3) ) THEN
                     IRH = 2
                  ELSE IF ( relativeHumidity(I,J,L) <= RH(4) ) THEN
                     IRH = 3
                  ELSE IF ( relativeHumidity(I,J,L) <= RH(5) ) THEN
                     IRH = 4
                  ELSE
                     IRH = 5
                  ENDIF

                  ! For the NRHth bin, we don't have to interpolate
                  ! For the other bins, we have to interpolate
                  IF ( IRH == NRH_b ) THEN
                     SCALEQ = QW(NRH_b) / QW(1)  !QW(1) is dry extinction eff.
                     REFF   = RW(NRH_b)

                  ELSE

                     ! Interpolate between different RH
                     FRAC = (relativeHumidity(I,J,L)-RH(IRH)) / (RH(IRH+1)-RH(IRH))
                     IF ( FRAC > 1.0d0 ) FRAC = 1.0d0

                     SCALEQ = (FRAC*QW(IRH+1) + (1.d0-FRAC)*QW(IRH)) / QW(1)
                     REFF   =  FRAC*RW(IRH+1) + (1.d0-FRAC)*RW(IRH)

                  ENDIF

                  SCALER                    = REFF / RW(1)
                  SCALEOD (I,J,L,IRH)       = SCALEQ * SCALER * SCALER
                  SCALEVOL(I,J,L)           = SCALER**3
                  eradius (I,J,L,NSADdust+N) = 1.0D-4 * REFF

                  !==============================================================
                  ! optDepth Diagnostic:
                  !
                  ! Computed here:
                  ! --------------
                  ! #7  Hygroscopic growth of SO4                [unitless]
                  ! #10 Hygroscopic growth of Black Carbon       [unitless]
                  ! #13 Hygroscopic growth of Organic Carbon     [unitless]
                  ! #16 Hygroscopic growth of Sea Salt (accum)   [unitless]
                  ! #19 Hygroscopic growth of Sea Salt (coarse)  [unitless]
                  !==============================================================
		  do_calc = .FALSE.

                  IF(do_synoz) THEN
		   IF(concentration(isynoz_num)%pArray3D(I,J,L) < synoz_threshold) do_calc = .TRUE.
		  ELSE
		   IF(pres3c(I,J,L) >= tropp(I,J)) do_calc = .TRUE.
                  END IF

		  IF(do_calc) optDepth(I,J,L,4+3*N) = optDepth(I,J,L,4+3*N) + SCALEOD(I,J,L,IRH)

               ENDDO
            ENDDO
         ENDDO

         !==============================================================
         ! Convert concentration [kg/m3] to optical depth [unitless].
         !
         ! odAer = ( 0.75 * gridBoxHeight * AERSL * qaa ) /
         !         ( MSDENS * raa * 1e-6 )
         ! (see Tegen and Lacis, JGR, 1996, 19237-19244, eq. 1)
         !
         ! Units ==> AERSL    [ kg/m3    ]
         !           MSDENS   [ kg/m3    ]
         !           raa      [ um       ]
         !           gridBoxHeight [ m        ]
         !           qaa      [ unitless ]
         !           odAer    [ unitless ]
         !
         ! NOTES:
         ! (1 ) Do the calculation at qaa(4,:) (i.e. 999 nm).
         ! (2 ) raa is the 'effective radius', Hansen and Travis, 1974
         ! (3 ) Report at the more relevant qaa(2,:) (i.e. 400 nm)
         !       Although SCALEOD would be slightly different at 400nm
         !       than at 1000nm as done here, FAST-J does currently
         !       allow one to provide different input optical depths at
         !       different wavelengths.  Therefore the reported value at
         !       determined with qaa(2,:) is as used in FAST-J.
         ! (4 ) Now use explicit indices in parallel DO-loops, since
         !       some compilers may not like array masks in parallel
         !       regions (bmy, 2/28/02)
         !==============================================================

         DO R = 1, NRH_b

            ! Bin for aerosol type and relative humidity
            IRHN = ( (N-1) * NRH_b ) + R

            ! Save aerosol optical depth for each combination
            ! of aerosol type and relative humidity into odAer,
            ! which will get passed to the FAST-J routines
            DO L = k1, k2
               DO J = ju1, j2
                  DO I = i1, i2
                     odAer(I,J,L,IRHN) = SCALEOD(I,J,L,R)  &
     &                           * 0.75d0 * gridBoxHeight(I,J,L)  &
     &                           * wAersl(I,J,L,N) * qaa_b(4,IND(N)) /  &
     &                           ( MSDENS(N) * raa_b(4,IND(N)) * 1.0D-6 )
                  ENDDO
               ENDDO
            ENDDO
         ENDDO

         
         !==============================================================
         !  Calculate Aerosol Surface Area
         !
         !  Units ==> AERSL    [ kg aerosol m^-3 air ]
         !            MSDENS   [ kg aerosol m^-3 aerosol ]
         !            eradius  [ cm      ]
         !            tArea    [ cm^2 dry aerosol/cm^3 air ]
         !
         !  Note: first find volume of aerosol (cm^3 arsl/cm^3 air), then
         !        multiply by 3/radius to convert to surface area in cm^2
         !
         !  Wet Volume = AERSL * SCALER**3 / MSDENS
         !  Wet Surface Area = 3 * (Wet Volume) / eradius
         !
         !  Effective radius for surface area and optical depths
         !  are identical.
         !==============================================================


         ! Store aerosol surface areas in tArea, and be sure
         ! to list them following the dust surface areas
         DO L = k1, k2
            DO J = ju1, j2
               DO I = i1, i2
                  tArea(I,J,L,N+NSADdust) = 3.D0*wAersl(I,J,L,N) *  &
     &                 SCALEVOL(I,J,L) /  &
     &                    ( eradius (I,J,L,NSADdust+N) * MSDENS(N) )
               ENDDO
            ENDDO
         ENDDO

      ENDDO  !Loop over NSADaer


      !==============================================================
      ! Account for hydrophobic aerosols (BC and OC), N=2 and N=3
      !==============================================================
      DO N = 2, 3

         ! Index for combination of aerosol type and RH
         IRHN = ( (N-1) * NRH_b ) + 1

         ! Aerosol optical depth
         DO L = k1, k2
            DO J = ju1, j2
               DO I = i1, i2
                  odAer(I,J,L,IRHN) = odAer(I,J,L,IRHN) +  &
     &                0.75d0 * gridBoxHeight(I,J,L) *  &
     &                dAersl(I,J,L,N-1) * qaa_b(4,IND(N))   /  &
     &                ( MSDENS(N) * raa_b(4,IND(N)) * 1.0D-6 )
               ENDDO
            ENDDO
         ENDDO

         ! Effective radius
         REFF = 1.0D-4 * raa_b(4,IND(N))

         DO L = k1, k2
            DO J = ju1, j2
               DO I = i1, i2

	          do_calc = .FALSE.

	          IF(do_synoz) THEN
                   IF (concentration(isynoz_num)%pArray3D(I,J,L) < synoz_threshold) do_calc = .TRUE.
                  ELSE
		   IF(pres3c(I,J,L) >= tropp(I,J)) do_calc = .TRUE.
                  END IF


                  IF(do_calc) THEN
                     ! Dry surface area
                     dryArea = 3.D0 * dAersl(I,J,L,N-1) / ( REFF * MSDENS(N) )

                     ! Add surface area to tArea array
                     tArea(I,J,L,NSADdust+N) = tArea(I,J,L,NSADdust+N) + dryArea

                     IF ((dryArea /= 0.0d0) .and. (tArea(I,J,L,NSADdust+N) /= 0.0d0)) THEN
                        ! Define a new effective radius that accounts
                        ! for the hydrophobic aerosol
                        eradius(I,J,L,NSADdust+N) = ( eradius(I,J,L,NSADdust+N) *  &
     &                                   tArea  (I,J,L,NSADdust+N) +  &
     &                                   REFF * dryArea)          /  &
     &                                 ( tArea  (I,J,L,NSADdust+N) + dryArea )
                     END IF

                  END IF

               END DO
            END DO
         END DO

      ENDDO

      !==============================================================
      ! optDepth Diagnostic: Aerosol OD's, Growth Rates, Surface Areas
      !
      ! Computed at the end of this routine
      ! -----------------------------------
      ! #1: Cloud optical depths (1000 nm)
      ! #2: Max Overlap Cld Frac
      ! #3: Random Overlap Cld Frac
      !
      ! Computed in the routine Dust_OptDep_SurfArea
      ! --------------------------------------------
      ! #4: Dust optical depths (400 nm)
      ! #5: Dust surface areas
      !
      ! Computed previously (above)
      ! ---------------------------------
      ! #7  Hygroscopic growth of SO4                [unitless]
      ! #10 Hygroscopic growth of Black Carbon       [unitless]
      ! #13 Hygroscopic growth of Organic Carbon     [unitless]
      ! #16 Hygroscopic growth of Sea Salt (accum)   [unitless]
      ! #19 Hygroscopic growth of Sea Salt (coarse)  [unitless]
      !
      ! Computed here:
      ! ---------------------------------
      ! #6  Sulfate Optical Depth (400 nm)           [unitless]
      ! #8  Sulfate Surface Area                     [cm2/cm3 ]
      ! #9  Black Carbon Optical Depth (400 nm)      [unitless]
      ! #11 Black Carbon Surface Area                [cm2/cm3 ]
      ! #12 Organic Carbon Optical Depth (400 nm)    [unitless]
      ! #14 Organic Carbon Surface Area              [cm2/cm3 ]
      ! #15 Sea Salt (accum) Opt Depth (400 nm)      [unitless]
      ! #17 Sea Salt (accum) Surface Area            [cm2/cm3 ]
      ! #18 Sea Salt (coarse) Opt Depth(400 nm)      [unitless]
      ! #20 Sea Salt (coarse) Surface Area           [cm2/cm3 ]
      !
      ! NOTE: The cloud optical depths are actually recorded at
      !       1000 nm, but vary little with wavelength.
      !==============================================================

      ! Loop over aerosol types
      DO N = 1, NSADaer

         !------------------------------------
         ! Optical Depths
         ! Scale of optical depths w/ RH
         !------------------------------------
         DO R = 1, NRH_b
            ! Index for type of aerosol and RH value
            IRHN = ( (N-1) * NRH_b ) + R

            ! Optical Depths
	   IF(do_synoz) THEN
            WHERE (concentration(isynoz_num)%pArray3D(:,:,:) < synoz_threshold)
               optDepth(:,:,:,3+3*N) = optDepth(:,:,:,3+3*N) +  &
     &                                    odAer(:,:,:,IRHN ) *  &
     &                                  qaa_b(2,IND(N)) /qaa_b(4,IND(N))
            END WHERE
	   ELSE
            DO L = k1, k2
             WHERE(pres3c(:,:,L) >= tropp(:,:))
              optDepth(:,:,L,3+3*N) = optDepth(:,:,L,3+3*N) +  &
     &                                    odAer(:,:,L,IRHN ) *  &
     &                                  qaa_b(2,IND(N)) /qaa_b(4,IND(N))
             END WHERE
            END DO
	   END IF
         ENDDO

         !------------------------------------
         ! Surface areas
         !------------------------------------
	 IF(do_synoz) THEN
          WHERE (concentration(isynoz_num)%pArray3D(:,:,:) < synoz_threshold)
             optDepth(:,:,:,5+3*N) = optDepth(:,:,:,5+3*N) +  &
     &                               tArea(:,:,:,N+NSADdust)
          END WHERE
         ELSE
          DO L = k1, k2
           WHERE(pres3c(:,:,L) >= tropp(:,:))
             optDepth(:,:,L,5+3*N) = optDepth(:,:,L,5+3*N) +  &
     &                               tArea(:,:,L,N+NSADdust)
           END WHERE
          END DO
         END IF

      ENDDO

      !=================================================================
      ! Choose if the radiative effects or/and heterogeneous chemistry
      ! on different aerosols/dust are turned on/off.
      !  AerDust_Effect_opt =
      !     0: radiative effects on  and heterogeneous chemistry on
      !     1: radiative effects off and heterogeneous chemistry on
      !     2: radiative effects on  and heterogeneous chemistry off
      !     3: radiative effects off and heterogeneous chemistry off
      !=================================================================

      ! Turn off the radiative effects of different aerososl
      IF ((AerDust_Effect_opt == 1) .or. (AerDust_Effect_opt == 3)) THEN
         DO R = 1,NRH_b
            odAer(:,:,:,R)         = 0.d0  !sulfate
            odAer(:,:,:,R+NRH_b)   = 0.d0  !BC
            odAer(:,:,:,R+2*NRH_b) = 0.d0  !OC
            odAer(:,:,:,R+3*NRH_b) = 0.d0  !SS(accum)
            odAer(:,:,:,R+4*NRH_b) = 0.d0  !SS(coarse)
         ENDDO
      END IF

      ! To turn off heterogeneous chemistry on different aerosols
      IF ((AerDust_Effect_opt == 2) .or. (AerDust_Effect_opt == 3)) THEN
         tArea(:,:,:,NSADdust+1) = 0.d0  !Sulfate
         tArea(:,:,:,NSADdust+2) = 0.d0  !BC
         tArea(:,:,:,NSADdust+3) = 0.d0  !OC
         tArea(:,:,:,NSADdust+4) = 0.d0  !SS (accum)
         tArea(:,:,:,NSADdust+5) = 0.d0  !SS (coarse)
      END IF

      return

      end subroutine Aero_OptDep_SurfArea
!EOC
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Dust_OptDep_SurfArea
!
! !INTERFACE:
!
      subroutine Dust_OptDep_SurfArea(gridBoxHeight, concentration, &
     &    tropp, pres3c, optDepth, eradius, tArea, Odmdust, Dust, &
     &    raa_b, qaa_b, do_synoz, isynoz_num, synoz_threshold, AerDust_Effect_opt, &
     &    i1, i2, ju1, j2, k1, k2, num_species, num_AerDust)
!
! !USES:
      use GmiArrayBundlePointer_mod, only : t_GmiArrayBundle
!
      implicit none

#     include "setkin_par.h"
#     include "gmi_AerDust_const.h"
!
! !INPUT PARAMETERS:
  integer, intent(in)    :: i1, i2, ju1, j2, k1, k2
  integer, intent(in)    :: num_species
  integer, intent(in)    :: num_AerDust
  integer, intent(in)    :: AerDust_Effect_opt
  logical, intent(in)    :: do_synoz
  integer, intent(in)    :: isynoz_num
  real*8 , intent(in)    :: synoz_threshold
  real*8 , intent(in)    :: gridBoxHeight(i1:i2,   ju1:j2,  k1:k2)
  type (t_GmiArrayBundle), intent(in) :: concentration(num_species)

  REAL*8, INTENT(IN) :: tropp(i1:i2, ju1:j2)         !Tropopause pressure (hPa)
  REAL*8, INTENT(IN) :: pres3c(i1:i2, ju1:j2, k1:k2) !Layer mean pressure (hPa)

  REAL*8 , intent(inOut) :: optDepth(i1:i2, ju1:j2, k1:k2, num_AerDust)
  REAL*8 , intent(inOut) :: eradius (i1:i2, ju1:j2, k1:k2, NSADdust+NSADaer)
  REAL*8 , intent(inOut) :: tArea   (i1:i2, ju1:j2, k1:k2, NSADdust+NSADaer)
  REAL*8 , intent(inOut) :: ODmdust (i1:i2, ju1:j2, k1:k2, NSADdust)
  REAL*8 :: Dust    (i1:i2, ju1:j2, k1:k2, NSADdust)
  REAL*8 :: raa_b   (4, NP_b)
  REAL*8 :: qaa_b   (4, NP_b)
!
! !DESCRIPTION: Use the mineral dust concentrations to calculate
!   the dust optical depth at each level.
!
!   This routine is a version of RDUST from GEOS-CHEM.
!   Here we do not read global mineral dust concentrations.
!
!   We only compute the optical depth and surface area at tropospheric
!   levels. We use the conditions:
!
!                WHERE (const(:,:,:,isynoz_num) < synoz_threshold)
!
! !LOCAL VARIABLES:
      REAL*8       :: MSDENS(NSADdust)

      integer      :: i, j , k, l, n, r
!
! !REVISION HISTORY:
!   February2005, Jules Kouatchou (Jules.Kouatchou.1@gsfc.nasa.gov)
!     Original code.
!   January 2008, Eric Nielsen 
!     Enable use of troposphere flag when do_synoz=.FALSE.
!EOP
!-------------------------------------------------------------------------
!BOC

      !==============================================================
      ! Convert concentration [kg/m3] to optical depth [unitless].
      !
      ! odmDust = ( 0.75 * BXHEIGHT * CONC * qaa ) /
      !           ( MSDENS * raa * 1e-6 )
      ! (see Tegen and Lacis, JGR, 1996, 19237-19244, eq. 1)
      !
      !  Units ==> dust     [ kg/m3    ]
      !            MSDENS   [ kg/m3    ]
      !            raa      [ um       ]
      !            BXHEIGHT [ m        ]
      !            qaa      [ unitless ]
      !            odmDust  [ unitless ]
      !
      ! NOTES:
      ! (1) Do the calculation at qaa(4,:) (i.e. 999 nm).
      !==============================================================
      MSDENS(1) = 2500.0
      MSDENS(2) = 2500.0
      MSDENS(3) = 2500.0
      MSDENS(4) = 2500.0
      MSDENS(5) = 2650.0
      MSDENS(6) = 2650.0
      MSDENS(7) = 2650.0

      DO N = 1, NSADdust
         DO L = k1, k2
            DO J = ju1, j2
               DO I = i1, i2
                  odmDust(I,J,L,N) = 0.75d0 * gridBoxHeight(I,J,L) *  &
     &                      dust(I,J,L,N) * qaa_b(4,14+N)     /  &
     &                     ( MSDENS(N) * raa_b(4,14+N) * 1.0D-6 )
               ENDDO
            ENDDO
         ENDDO
      ENDDO


      !==============================================================
      ! Calculate Dust Surface Area
      !
      ! Units ==> dust     [ kg dust/m^3 air    ]
      !           MSDENS   [ kg dust/m^3 dust   ]
      !           raa      [ um                 ]
      !           tArea    [ cm^2 dust/cm^3 air ]
      !           eradius  [ cm                 ]
      !
      ! NOTE: first find volume of dust (cm3 dust/cm3 air), then
      !       multiply by 3/radius to convert to surface area in cm2
      !
      ! tArea(:,1:NSADdust) and eradius(:,1:NSADdust) are for
      ! the NSADdustAST-J dust wavelength bins (read into dust)
      !==============================================================
      DO N = 1, NSADdust
         DO L = k1, k2
            DO J = ju1, j2
               DO I = i1, i2
                  eradius(I,J,L,N) = raa_b(4,14+N) * 1.0D-4

                  tArea(I,J,L,N)   = 3.D0 / eradius(I,J,L,N) *  &
     &                      dust(I,J,L,N) / MSDENS(N)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      !==============================================================
      ! optDepth Diagnostic:
      !
      ! Tracer #1: Cloud optical depths    (from Aero_OptDep_SurfArea)
      ! Tracer #2: Max Overlap Cld Frac    (from Aero_OptDep_SurfArea)
      ! Tracer #3: Random Overlap Cld Frac (from Aero_OptDep_SurfArea)
      ! Tracer #4: Dust optical depths at 400 nm (from all size bins)
      ! Tracer #5: Dust surface areas (from all size bins)
      !==============================================================

      DO N = 1, NSADdust
        IF(do_synoz) THEN
         WHERE (concentration(isynoz_num)%pArray3D(:,:,:) < synoz_threshold)
             !--------------------------------------
             ! optDepth tracer #4: Dust optical depths
             !--------------------------------------
             optDepth(:,:,:,4) = optDepth(:,:,:,4) +  &
     &             (odmDust(:,:,:,N) * qaa_b(2,14+N) / qaa_b(4,14+N))
             !--------------------------------------
             ! optDepth tracer #5: Dust surface areas
             !--------------------------------------
             optDepth(:,:,:,5) = optDepth(:,:,:,5) + tArea(:,:,:,N)
         END WHERE
	ELSE
         DO L = k1, k2
	   WHERE(pres3c(:,:,L) >= tropp(:,:))
             !--------------------------------------
             ! optDepth tracer #4: Dust optical depths
             !--------------------------------------
             optDepth(:,:,L,4) = optDepth(:,:,L,4) +  &
     &             (odmDust(:,:,L,N) * qaa_b(2,14+N) / qaa_b(4,14+N))
             !--------------------------------------
             ! optDepth tracer #5: Dust surface areas
             !--------------------------------------
             optDepth(:,:,L,5) = optDepth(:,:,L,5) + tArea(:,:,L,N)
	   END WHERE
         END DO
	END IF
      ENDDO

      !=================================================================
      ! Choose if the radiative effects or/and heterogeneous chemistry
      ! on different aerosols/dust are turned on/off.
      !  AerDust_Effect_opt =
      !     0: radiative effects on  and heterogeneous chemistry on
      !     1: radiative effects off and heterogeneous chemistry on
      !     2: radiative effects on  and heterogeneous chemistry off
      !     3: radiative effects off and heterogeneous chemistry off
      !=================================================================

      ! To turn off the radiative effects of dust
      IF ((AerDust_Effect_opt == 1) .or. (AerDust_Effect_opt == 3)) THEN
         DO R = 1,NSADdust
            odmDust(:,:,:,R)       = 0.d0  !dust
         ENDDO
      ENDIF

      ! To turn off heterogeneous chemistry on dust
      IF ((AerDust_Effect_opt == 2) .or. (AerDust_Effect_opt == 3)) THEN
         DO R = 1,NSADdust
           tArea(:,:,:,R) = 0.d0  !dust
         ENDDO
      END IF
!
      RETURN

      END subroutine Dust_OptDep_SurfArea
!EOC
!-------------------------------------------------------------------------
      end module GmiAerDustODSA_mod
