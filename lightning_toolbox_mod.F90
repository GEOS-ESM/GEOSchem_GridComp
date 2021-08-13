!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: lightning_toolbox_mod.F90
!
! !DESCRIPTION: Module lightning\_toolbox\_mod.F90 is a wrapper module 
!  that contains routines to compute lightning flash rates as used by
!  HEMCO. It has no dependencies and can thus also be used outside of 
!  HEMCO. 
!\\
!\\
! References:
! \begin{itemize}
! \item Murray, L. T., Jacob, D. J., Logan, J. A., Hudman, R. C., and
!       Koshak, W. J.: \emph{Optimized regional and interannual variability
!       of lightning in a global chemical transport model con- strained
!       by LIS/OTD satellite data}, \underline{J. Geophys. Res.},
!       Atmospheres, 117, 2012.
! \item Ott, L. E., K. E. Pickering, G. L. Stenchikov, D. J. Allen,
!       A. J. DeCaria, B. Ridley, R.-F. Lin, S. Lang, and W.-K. Tao,
!       \emph{Production of lightning NOx and its vertical distribution
!       calculated  from three-dimensional cloud-scale chemical transport
!       model simulations}, \underline{J. Geophys. Res.}, 115, D04301, 2010.
! \end{itemize}
!
! !INTERFACE:
!
MODULE Lightning_Toolbox_Mod
!
! !USES:
!
  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: CalcFlashRate 
  PUBLIC  :: LightDist 
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: FLASHES_CTH
  PRIVATE :: GET_IC_CG_RATIO
!
! !PUBLIC DATA MEMBERS:
!
! !REVISION HISTORY:
!  14 Apr 2004 - L. Murray, R. Hudman - Initial version
!  10 Jun 2019 - C. Keller   - Separated flash rate calculation from
!                              hcox_lightnox_mod.F90 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
  INTEGER, PARAMETER :: hp = KIND( REAL( 0.0, 8 ) )

  ! Parameter for LFR calculation
  REAL*8,  PARAMETER            :: RFLASH_MIDLAT = 3.011d26   ! 500 mol/flash
  REAL*8,  PARAMETER            :: RFLASH_TROPIC = 1.566d26   ! 260 mol/flash
  REAL*8,  PARAMETER            :: EAST_WEST_DIV = -30d0
  REAL*8,  PARAMETER            :: WEST_NS_DIV   =  35d0
  REAL*8,  PARAMETER            :: EAST_NS_DIV   =  35d0
  REAL*8,  PARAMETER            :: T_NEG_BOT     = 273.0d0    !   0 C 
  REAL*8,  PARAMETER            :: T_NEG_CTR     = 258.0d0    ! -15 C
  REAL*8,  PARAMETER            :: T_NEG_TOP     = 233.0d0    ! -40 C
  REAL*8,  PARAMETER            :: Rdg0          = 287.0e+0_hp / 9.80665e+0_hp ! Rd/g0 

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: CalcFlashRate 
!
! !DESCRIPTION: Subroutine CalcFlashRate computes the flash rate as described
! in Murray et al. 2012. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CalcFlashRate( am_I_Root, IM, JM, LM, AREA_M2, LWI2D, & 
                            OTDLIS2D, XMID2D, YMID2D, TK3D, PEDGE3D, &
                            BXHEIGHT3D, CNVMFC3D, BYNCY3D, CNVFRC2D,   &
                            LCNVFRC, OTD_LIS_SCALE, LNOX2D, LFR2D, IC2D, &
                            CG2D, LTOP2D, H02D, RC )
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
    LOGICAL,         INTENT(IN   )  :: am_I_Root             ! Root CPU?
    INTEGER,         INTENT(IN   )  :: IM, JM, LM            ! # of x,y,z levels
    REAL(hp),        INTENT(IN   )  :: AREA_M2(IM,JM)        ! Surface area in m2
    INTEGER,         INTENT(IN   )  :: LWI2D(IM,JM)          ! 0=land, 1=water, 2=ice 
    REAL(hp),        INTENT(IN   )  :: OTDLIS2D(IM,JM)       ! OTD-LIS redistribution factor. Ignored if negative. 
    REAL(hp),        INTENT(IN   )  :: XMID2D(IM,JM)         ! Grid box longitude (degrees) 
    REAL(hp),        INTENT(IN   )  :: YMID2D(IM,JM)         ! Grid box latitude (degrees)
    REAL(hp),        INTENT(IN   )  :: TK3D(IM,JM,LM)        ! Air temperature (Kelvin) 
    REAL(hp),        INTENT(IN   )  :: PEDGE3D(IM,JM,LM+1)   ! Pressure edges (Pa) 
    REAL(hp),        INTENT(IN   )  :: BXHEIGHT3D(IM,JM,LM)  ! Grid box heights (m) 
    REAL(hp),        INTENT(IN   )  :: CNVMFC3D(IM,JM,LM+1)  ! convective cloud mass flux (kg/m2/s) 
    REAL(hp),        INTENT(IN   )  :: BYNCY3D(IM,JM,LM)     ! buoyancy (ms-2)
    REAL(hp),        INTENT(IN   )  :: CNVFRC2D(IM,JM)       ! convective fraction (unitless) 
    LOGICAL ,        INTENT(IN   )  :: LCNVFRC               ! determine H0 from convective fraction and buoyancy? 
    REAL(hp),        INTENT(IN   )  :: OTD_LIS_SCALE         ! Global LFR scale factor 
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,         INTENT(INOUT)  :: RC              ! Return code (0=success, 1=fail)
!
! !OUTPUT PARAMETERS:
!
    REAL(hp),        INTENT(  OUT)  :: LNOX2D(IM,JM)    ! Total lightning NOx molecules per 6 hours 
    REAL(hp),        INTENT(  OUT)  :: LFR2D(IM,JM)    ! Flash rate (flashes/min/km2)
    REAL(hp),        INTENT(  OUT)  :: IC2D(IM,JM)     ! Intra-cloud flash rate (flashes/min/km2)
    REAL(hp),        INTENT(  OUT)  :: CG2D(IM,JM)     ! Cloud-to-ground flash rate (flashes/min/km2)
    INTEGER,         INTENT(  OUT)  :: LTOP2D(IM,JM)   ! Lightning top level 
    REAL(hp),        INTENT(  OUT)  :: H02D(IM,JM)     ! Convective cloud top height (m) 
! 
! !REVISION HISTORY: 
!  10 May 2006 - L. Murray - Initial version  
!  10 Jun 2019 - C. Keller - Separated flash rate calculation from
!                            subroutine LIGHTNOX in hcox_lightnox_mod.F90 
!  8 Aug 2019  - M. Damon -  Put flash rate routine into GMAO Shared.
!                            Wrapper routine from "Lightning_mod.F90"
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!   
    INTEGER           :: I,         J,           L,        LCHARGE
    INTEGER           :: LMAX,      LTOP,        LBOTTOM,  A_KM2
    INTEGER           :: LTOP1,     LTOP2
    REAL*8            :: XMID,      YMID,        CC,       DLNP
    REAL*8            :: DZ,        FLASHRATE,   H0,       HBOTTOM
    REAL*8            :: HCHARGE,   IC_CG_RATIO, Z_CG,     ZUP
    REAL*8            :: RATE_SAVE, REDIST,      T1,       T2
    REAL*8            :: TOTAL,     TOTAL_CG,    TOTAL_IC, X
    REAL*8            :: P1,        P2,          P3,       RATE
    REAL*8            :: Z_IC
    INTEGER           :: SFCTYPE

    !=================================================================
    ! CalcFlasRate begins here!
    !=================================================================

    ! Assume success
    RC = 0

!     print*, "Hola from HEMCO lightning: ", IM, JM, LM, am_I_Root
!     print*, "HEMCO surface area: ", minval(AREA_M2), maxval(AREA_M2)
!     print*, "HEMCO LWI: ", minval(LWI2D), maxval(LWI2D)
!     print*, "HEMCO TROPP: ", minval(TROPP2D), maxval(TROPP2D)
!     print*, "HEMCO OTDLIS2D: ", minval(OTDLIS2D), maxval(OTDLIS2D) 
!     print*, "HEMCO lon: ", minval(XMID2D), maxval(XMID2D)
!     print*, "HEMCO lat: ", minval(YMID2D), maxval(YMID2D)
!     print*, "HEMCO Temp: ", minval(TK3D), maxval(TK3D)
!     print*, "HEMCO Pedge: ", minval(PEDGE3D), maxval(PEDGE3D)
!     print*, "HEMCO box height: ", minval(BXHEIGHT3D), maxval(BXHEIGHT3D)
!     print*, "HEMCO CNVMFC: ", minval(CNVMFC3D), maxval(CNVMFC3D)
!     print*, "HEMCO buoyancy: ", minval(BYNCY3D), maxval(BYNCY3D)
!     print*, "HEMCO CNVFRC2D: ", minval(CNVFRC2D), maxval(CNVFRC2D)
!     print*, "HEMCO LCNVFRC: ", LCNVFRC
!     print*, "HEMCO OTD_LIS_SCALE: ", OTD_LIS_SCALE


    ! Reset arrays 
    LNOX2D = 0.0_hp
    LFR2D  = 0.0_hp
    IC2D   = 0.0_hp
    CG2D   = 0.0_hp
    H02D   = 0.0_hp
    LTOP2D = 0 


    ! LMAX: the highest L-level to look for lightnox (usually LLPAR-1)
    LMAX   = LM - 1

    !=================================================================
    ! Compute flashrate for each (I,J) column
    !=================================================================

!$OMP PARALLEL DO                                                   &
!$OMP DEFAULT( SHARED )                                             &
!$OMP PRIVATE( I,           J,        L,        A_KM2             ) &
!$OMP PRIVATE( LCHARGE,     P1,       P2,       XMID,   YMID      ) &
!$OMP PRIVATE( T1,          T2,       DLNP,     DZ,     P3        ) &
!$OMP PRIVATE( ZUP,         HCHARGE,  LTOP,     H0,     Z_CG      ) &
!$OMP PRIVATE( Z_IC,        LBOTTOM,  HBOTTOM,  CC,     FLASHRATE ) &
!$OMP PRIVATE( IC_CG_RATIO, RATE                                  ) &
!$OMP PRIVATE( X,           TOTAL_IC, TOTAL_CG, TOTAL,  REDIST    ) &
!$OMP PRIVATE( RATE_SAVE,   SFCTYPE                               ) &
!$OMP PRIVATE( LTOP1,       LTOP2                                 ) &
!$OMP SCHEDULE( DYNAMIC )

    ! Loop over surface boxes
    DO J = 1, JM 
    DO I = 1, IM 

       ! Grid box surface area in [km2]
       A_KM2    = AREA_M2( I, J ) / 1d6

       ! Grid box latitude and longitude [degrees]
       YMID     = YMID2D( I, J )
       XMID     = XMID2D( I, J )

       ! Make sure xmid is between -180 and +180
       IF ( XMID >= 180.0_hp ) XMID = XMID - 360.0_hp

       ! Get surface type. Note that these types are different than 
       ! the types used elsewhere: 0 = land, 1=water, 2=ice!
       SFCTYPE = LWI2D(I,J)

       ! Initialize
       LBOTTOM       = 0 
       LCHARGE       = 0
       CC            = 0d0
       HCHARGE       = 0d0
       HBOTTOM       = 0d0
       TOTAL         = 0d0
       TOTAL_IC      = 0d0
       TOTAL_CG      = 0d0
       RATE          = 0d0
       H0            = 0d0
       IC_CG_RATIO   = 0d0
       
       ! Get factors for OTD-LIS local redistribution or none.
       ! This constrains the seasonality and spatial distribution
       ! of the parameterized lightnox to match the HRMC v2.2
       ! product from LIS/OTD, while still allowing the model to
       ! place lightnox locally within deep convective events.
       ! (ltm, bmy, 1/31/07)
       REDIST = OTDLIS2D(I,J)
       IF ( REDIST < 0.0d0 ) THEN
          REDIST = 1.0d0
       ENDIF

       !===========================================================
       ! (1) FIND NEGATIVE CHARGE LAYER
       !
       ! LCHARGE is the L-value where the negative charge layer is
       ! found.  According to Williams (1985), the negative charge
       ! layer occurs where T is between 0 C and -40 C.  The 
       ! original model code set this at -10 C, but according to 
       ! Houze (1993), a good proxy for the negative charge layer
       ! maximum density is at -15 C.
       !
       ! Also of interest for later, will be the bottom of the
       ! negative charge layer (i.e., temp = 0 C) in calculating 
       ! the cold cloud depth.
       !
       ! If LCHARGE=1, then it is too cold to have water droplets
       ! in the column, so there will be no lightnox events,
       ! and we go to the next (I,J) box.
       !
       ! (ltm, bmy, 5/10/06, 12/11/06)
       !===========================================================

       ! Find negative charge layer
       DO L = 1, LMAX

          IF ( TK3D(I,J,L) <= T_NEG_CTR ) THEN
             LCHARGE = L
             EXIT
          ENDIF
       ENDDO

       ! Error check LCHARGE
       LCHARGE = MIN(LCHARGE,LMAX)


       IF ( LCHARGE > 1 ) THEN



          !-----------------------------------------------------------
          ! (1a) Define more quantities
          !-----------------------------------------------------------
              
          ! Pressure [Pa] at the centers of grid
          ! boxes (I,J,LCHARGE-1) and (I,J,LCHARGE)
          ! Now calculate from grid edges (ckeller, 10/06/2014)
          P1 = ( PEDGE3D(I,J,LCHARGE-1) &
               + PEDGE3D(I,J,LCHARGE  ) ) / 2.0_hp
          P2 = ( PEDGE3D(I,J,LCHARGE  ) &
               + PEDGE3D(I,J,LCHARGE+1) ) / 2.0_hp
   
          ! Temperatures [K] at the centers of grid
          ! boxes (I,J,LCHARGE-1) and (I,J,LCHARGE)
          T1   = TK3D(I,J,LCHARGE-1)
          T2   = TK3D(I,J,LCHARGE  )
    
          ! DZ is the height [m] from the center of box (I,J,LCHARGE-1)
          ! to the negative charge layer.  It may be found in either
          ! the (LCHARGE)th sigma layer or the (LCHARGE-1)th layer.
          ! We use the hypsometric eqn to find the distance between
          ! the center of (LCHARGE)th and (LCHARGE-1)th boxes, then
          ! assume a linear temp distribution to scale between the two.
          DLNP = LOG( P1 / P2 ) / ( T1 - T2 ) * ( T1 - T_NEG_CTR )
          DZ   = Rdg0 * ( ( T1 + T2 ) / 2d0 ) * DLNP
   
          ! Pressure [Pa] at the bottom edge of box (I,J,LCHARGE),
          ! or, equivalently, the top edge of box (I,J,LCHARGE-1).
          P3   = PEDGE3D( I, J, LCHARGE )
   
          ! Height [m] from the center of grid box (I,J,LCHARGE-1) 
          ! to the top edge of grid box (I,J,LCHARGE-1)
          ZUP  = Rdg0 * T1 * LOG( P1 / P3 )
   
          !-----------------------------------------------------------
          ! (1b) HCHARGE is the height of the negative charge layer 
          ! above the bottom edge of box (I,J,LCHARGE).  
          ! 
          ! If DZ < ZUP, then DZ is in grid box (I,J,LCHARGE-1);
          ! therefore subtract 1 from LCHARGE and compute HCHARGE 
          ! accordingly.  
          !
          ! In this case, please note that BXHEIGHT(I,J,LCHARGE)-ZUP 
          ! is the distance from the bottom edge of the grid box to 
          ! the center of the newly defined (LCHARGE)th layer.
          !-----------------------------------------------------------
          IF ( DZ >= ZUP ) THEN
             HCHARGE = DZ - ZUP
          ELSE
             LCHARGE = LCHARGE - 1
             HCHARGE = (BXHEIGHT3D(I,J,LCHARGE)-ZUP) + DZ
          ENDIF
    
          !===========================================================
          ! (2) COMPUTE CONVECTIVE CLOUD TOP HEIGHT
          !
          ! LTOP is the L-layer where the convective cloud top is 
          ! found.  The cloud top is located at the highest sigma 
          ! level for which the cloud mass flux is nonzero.  Since 
          ! GMAO cloud mass flux is defined at the top of each sigma 
          ! level, the convective cloud top is located at the top 
          ! edge of layer LTOP.
          !
          ! For lightnox to exist, the cloud must straddle the 
          ! negative charge layer (in other words, at the very 
          ! minimum, the cloud bottom must occur in the LCHARGEth 
          ! layer).  If LTOP < LCHARGE go to the next (I,J) location.
          ! 
          ! Additionally, because the negative charge layer extends 
          ! from 0 C to around -40 C (Williams 1985), any cloud type 
          ! heights that are not colder than -40 C will be considered 
          ! unable to create the necessary dipole.  Therefore, if 
          ! T(I,J,LTOP) >= -40 C, go to the next (I,J) location. 
          !
          ! To be easily translatable to an ESMF environment, we now 
          ! use the convective cloud mass flux to determine LTOP.
          ! Use the same definition as used in GEOS-Chem.
          !
          ! (ltm, bmy, 5/10/06, 12/11/06)
          !
          ! GEOS-FP turns off convection in grid boxes where vertical
          ! transport is explicitly resolved. The convective mass flux
          ! (which is computed within the convection code) is then zero
          ! in these grid boxes, even though convection did occur at 
          ! these places. 
          ! This may become increasingly relevant as GEOS-FP operates 
          ! at even higher resolutions. 
          ! If available, also determine cloud top height from 
          ! buoyancy and the convective fraction. Define it as the
          ! highest level with non-negative buoyancy and for columns 
          ! with non-zero convective fraction (ckeller, 3/04/16).
          !===========================================================
   
          ! 'Traditional definition of cloud top level
          LTOP1 = 1
          DO L = LM, 1, -1
             IF ( CNVMFC3D(I,J,L) > 0.0_hp ) THEN
                LTOP1 = L + 1
                EXIT
             ENDIF
          ENDDO 
   
          ! To determine cloud top height from buoyancy for all grid
          ! boxes with non-zero convective fraction (define cloud top
          ! as top level with positive buoyancy).
          LTOP2 = 0
          IF ( LCNVFRC ) THEN 
             IF ( CNVFRC2D(I,J) > 0.0_hp ) THEN
                DO L = LM, 1, -1
                   IF ( BYNCY3D(I,J,L) >= 0.0_hp ) THEN 
                      LTOP2= L + 1
                      EXIT
                   ENDIF
                ENDDO
             ENDIF 
          ENDIF 
   
          ! Take whichever value is higher
          LTOP = MAX(LTOP1,LTOP2)
   
          !----------------------------------------------------------------
          ! Error checks for LTOP 
          !----------------------------------------------------------------
   
          ! Error check LTOP
          !IF ( LTOP == 0 ) CYCLE
          IF ( LTOP > 0 ) THEN



             ! Error check LTOP as described above
             IF ( LTOP >  LMAX    ) LTOP = LMAX
             IF ( LTOP >= LCHARGE ) THEN
   
                ! Diagnose used LTOP
                LTOP2D(I,J) = LTOP
         
                ! H0 is the convective cloud top height [m].  This is the
                ! distance from the surface to the top edge of box (I,J,LTOP).
                H0   = SUM(BXHEIGHT3D(I,J,1:LTOP))
         
                ! Z_CG is the cloud-ground path (ground --> HCHARGE) [m]
                Z_CG = SUM(BXHEIGHT3D(I,J,1:LCHARGE-1)) + HCHARGE
         
                ! Z_IC is the intra-cloud path (HCHARGE --> cloud top) [m]
                Z_IC = SUM(BXHEIGHT3D(I,J,LCHARGE:LTOP)) - HCHARGE
         
                !===========================================================
                ! (3) COMPUTE COLD CLOUD THICKNESS
                ! 
                ! Find the cold cloud thickness (CC) -- the distance from 
                ! where the temperature is 0 C up to the top of the cloud.  
                ! This is necessary for calculating the f_CG/f_IC ratio as 
                ! per Price and Rind 1993.  
                !
                ! This is a clone of the method above to find height to 
                ! HCHARGE, and we can recycle many of the same variables 
                ! that aren't used again.
                ! 
                ! Grid box (I,J,LBOTTOM) is the model layer where the 
                ! temperature of the cloud is 0C.
                !
                ! NOTE: If no temperature in the column is above 0 C, it 
                ! moves on to the next (I,J) box as before with the -15 C.
                !
                ! (ltm, bmy, 5/10/06, 12/11/06)
                !===========================================================
         
                ! Find the level where T = 0 C
                DO L = 1, LMAX
                   IF ( TK3D(I,J,L) <= T_NEG_BOT ) THEN
                      LBOTTOM = L
                      EXIT
                   ENDIF
                ENDDO
         
                ! Error check LBOTTOM as described above
                LBOTTOM = MIN(LBOTTOM,LMAX)
                IF ( LBOTTOM > 1    ) THEN 



                   !-----------------------------------------------------------
                   ! (3a) Define more quantities
                   !-----------------------------------------------------------
            
                   ! Pressure [Pa] at the centers of grid
                   ! boxes (I,J,LBOTTOM-1) and (I,J,LBOTTOM)
                   ! Now calculate from grid edges (ckeller, 10/06/2014)
                   P1 = ( PEDGE3D(I,J,LBOTTOM-1) &
                        + PEDGE3D(I,J,LBOTTOM  ) ) / 2.0_hp
                   P2 = ( PEDGE3D(I,J,LBOTTOM  ) &
                        + PEDGE3D(I,J,LBOTTOM+1) ) / 2.0_hp
            
                   ! Temperature [K] at the centers of grid
                   ! boxes (I,J,LBOTTOM-1) and (I,J,LBOTTOM)
                   T1   = TK3D(I,J,LBOTTOM-1)
                   T2   = TK3D(I,J,LBOTTOM  )
                   
                   ! DZ is the height [m] from the center of box (I,J,LCHARGE-1)
                   ! to the negative charge layer.  It may be found in either
                   ! the (LCHARGE)th sigma layer or the (LCHARGE-1)th layer.
                   ! We use the hypsometric eqn to find the distance between
                   ! the center of (LCHARGE)th and (LCHARGE-1)th boxes, then
                   ! assume a linear temp distribution to scale between the two.
                   DLNP = LOG( P1 / P2 ) / ( T1 - T2 ) * ( T1 - T_NEG_BOT )
                   DZ   = Rdg0 * ( ( T1 + T2 ) / 2d0 ) * DLNP
            
                   ! Pressure [Pa] at the bottom edge of box (I,J,LBOTTOM),
                   ! or, equivalently, the top edge of box (I,J,BOTTOM-1).
                   P3   = PEDGE3D( I, J, LBOTTOM )
            
                   ! Height [m] from the center of grid box (I,J,LBOTTOM-1) 
                   ! to the top edge of grid box (I,J,LBOTTOM-1)
                   ZUP  = Rdg0 * T1 * LOG( P1 / P3 )
            
                   !-----------------------------------------------------------
                   ! (3b) HBOTTOM is the height of the 0 C layer above the 
                   ! bottom edge of box (I,J,LBOTTOM).  
                   ! 
                   ! If DZ < ZUP, then DZ is in grid box (I,J,LBOTTOM-1);
                   ! therefore subtract 1 from LBOTTOM and compute HBOTTOM
                   ! accordingly.
                   !
                   ! In this case, please note that BXHEIGHT(I,J,LBOTTOM)-ZUP 
                   ! is the distance from the bottom edge of the grid box to 
                   ! the center of the newly defined (LBOTTOM)th layer.
                   !-----------------------------------------------------------
                   IF ( DZ >= ZUP ) THEN
                      HBOTTOM = DZ - ZUP
                   ELSE
                      LBOTTOM = LBOTTOM - 1
                      HBOTTOM = (BXHEIGHT3D(I,J,LBOTTOM) - ZUP) + DZ
                   ENDIF
           
                   ! Cold cloud thickness is difference of cloud top 
                   ! height (H0) and the height to the bottom.
                   CC = H0 - SUM(BXHEIGHT3D(I,J,1:LBOTTOM-1) ) - HBOTTOM 
            
                   !===========================================================
                   ! (4) COMPUTE IC/CG FLASH_RATIO FROM COLD-CLOUD DEPTH
                   !
                   ! This is necessary as an input for the MFLUX and PRECON
                   ! parameterizations, as well as for determining the fraction 
                   ! of LNOX generated by either type of flash, and will
                   ! eventually be used for separate vertical distributions
                   ! when they become available.  (ltm, bmy, 12/11/06)
                   !===========================================================
            
                   ! Get Inter-Cloud/Cloud-Ground flash ratio [unitless]
                   IC_CG_RATIO = GET_IC_CG_RATIO( CC )
            
                   !===========================================================
                   ! (5) COMPUTE LIGHTNOX FLASH RATES
                   !
                   ! Now that we have computed the the ratio of intra-cloud
                   ! flashes to cloud-ground flashes, compute the lightnox
                   ! flash rate via one of these parameterizations:
                   !
                   ! (a) Cloud top height (CTH)
                   ! (b) Mass flux (MFLUX)
                   ! (c) Convective Precpitation (PRECON)
                   ! 
                   ! (ltm, bmy, 5/10/06, 12/11/06)
                   !===========================================================
            
                   !--------------------------------------------------------
                   ! (5a) CLOUD TOP HEIGHT PARAMETERIZATION (all met fields)
                   !
                   ! Based on Price & Rind [1992,1993,1994].
                   !--------------------------------------------------------
            
                   ! Get lightnox flash rate per minute and IC/CG ratio
                   CALL FLASHES_CTH( I, J, H0, FLASHRATE, SFCTYPE ) 

                   !===========================================================
                   ! (6) COMPUTE TOTAL LNOx AND PARTITION INTO VERTICAL LAYERS
                   ! 
                   ! (6a) We convert FLASHRATE (computed above) to units of
                   ! [flashes/6h] and store in the RATE variable.
                   !
                   ! We then multiply RATE by a scale factor based on 
                   ! OTD/LIS observations.  This is necessary in order to make 
                   ! sure that the lightnox flashes happen in the correct 
                   ! locations as diagnosed by OTD/LIS satellite observations.  
                   ! There are two redistribution options:
                   !
                   !   (1) Apply regional scale factors based on OTD/LIS
                   !        observations (method of L. Jourdain et al)
                   !
                   !   (2) Apply box-by-box scale scale factors based on
                   !        OTD/LIS observations (method of B. Sauvage)
                   !
                   ! NOTE: As of 3/27/07, only method (1) is implemented.
                   ! 
                   ! (6b) We then compute X, which is the ratio
                   !   [cloud-ground flashes / total flashes].
                   !
                   ! The amount of lightnox released will depend whether we
                   ! are in the tropics or in mid-latitudes.
                   !
                   !
                   ! (6c) LIGHTNOX NOx EMISSIONS IN THE TROPICS:
                   ! ----------------------------------------------------------
                   ! N. American / S. American     tropics: lat <= 23 N 
                   ! African / Oceanian / Eurasian tropics: lat <= 35 N
                   !
                   ! The lightnox NOx released in the inter-cloud (IC) and 
                   ! cloud-ground (CG) paths are given by:
                   !
                   !   TOTAL_IC = RFLASH_TROPIC * RATE * (1-X) * Z_IC
                   !   TOTAL_CG = RFLASH_TROPIC * RATE * (  X) * Z_CG
                   !
                   ! where:
                   !   RFLASH_TROPIC = # of NOx molecules released per flash 
                   !                    per meter (same as in previous code)
                   !   RATE          = lightnox flashes / 6h computed above
                   !   Z_IC          = IC pathway in meters (from the negative    
                   !                    cloud layer to the cloud top)
                   !   Z_CG          = CG pathway in meters (from the negative
                   !                    cloud layer to the ground surface)
                   !   
                   ! We also apply a top-down final global scaling factor, 
                   ! calculated by previously bringing total global LNOx to 
                   ! 6 Tg N/yr  (2x2.5: 0.3683, 4x5: 0.8996).  In 2004, the 
                   ! tropics-only contribution to LNOx was 4.5379 Tg N.
                   !
                   ! 
                   ! (6d) LIGHTING NOx EMISSIONS AT MIDLATITUDES:
                   ! ----------------------------------------------------------
                   ! N. American midlatitudes : lat > 23N
                   ! Eurasian    midlatitudes : lat > 35N
                   !
                   ! The lightnox NOx released at midlatitudes is independent
                   ! of path length.  Thus:
                   !
                   !   TOTAL_IC = RFLASH_MIDLAT * RATE * (1-X) * MID_LAT_SCALE
                   !   TOTAL_CG = RFLASH_MIDLAT * RATE *    X  * MID_LAT_SCALE
                   !
                   ! where 
                   !   RFLASH_MIDLAT = # of NOx molecules released per flash 
                   !                    per meter (based on 500 mol/flash)
                   !   RATE          = lightnox flashes / 6h computed above
                   !   Z_IC          = IC pathway in meters (from the negative    
                   !                    cloud layer to the cloud top)
                   !   Z_CG          = CG pathway in meters (from the negative
                   !                    cloud layer to the ground surface)
                   !
                   ! We now emit at the Northern Mid-latitudes using an RFLASH
                   ! value of 500 mol/flash.  This is independent of path 
                   ! length.  
                   !
                   ! NOTE: The OTD-LIS local redistribution method was expanded
                   ! upon from Sauvage et al, 2007, ACP.
                   ! http://www.atmos-chem-phys.net/7/815/2007/acp-7-815-2007.pdf
                   !
                   ! (6e) The total lightnox is the sum of IC+CG paths:
                   !   TOTAL = TOTAL_IC + TOTAL_CG
                   !
                   ! (6g) We then partition the NOx into each of the vertical 
                   ! grid boxes within the column with a Ken Pickering PDF
                   ! (see comments below).
                   !
                   ! (ltm, rch, bmy, 5/10/06, 3/27/07)
                   !===========================================================
            
                   !-----------------------------------------------------------
                   ! (6a) Compute flash rate and apply OTD/LIS redistribution
                   !-----------------------------------------------------------
            
                   ! Convert [flashes/min] to [flashes/6h]
                   RATE     = FLASHRATE * 360.0d0
                   
                   
                   ! Apply regional or local OTD-LIS redistribution so that the 
                   ! flashes occur in the right place. 
                   RATE = RATE * REDIST
            
                   ! Apply scaling factor to make sure annual average flash rate 
                   ! equals that of the climatology. (ltm, 09/24/07)
                   RATE = RATE * OTD_LIS_SCALE
            
                END IF ! LBOTTOM >  1
             END IF ! LTOP    >= LCHARGE 
          END IF ! LTOP    >  0 
       END IF ! LCHARGE >  1 
            


       ! Do not allow flash density to become unrealistically high
       ! Globally limit the flash rate to its highest observed value
       ! of 4.2e-3 flashes / km2 / s from the ENTLN global product
       ! (ltm, mps, 8/9/18)
       IF ( ( RATE / 21600. / A_KM2 ) > 0.004177159 ) THEN
          RATE = 0.004177159 * 21600. * A_KM2
       END IF
   
       ! Ratio of cloud-to-ground flashes to total # of flashes
       X    = 1d0 / ( 1d0 + IC_CG_RATIO )
   
       ! Compute LNOx emissions for tropics or midlats 
       IF ( XMID > EAST_WEST_DIV ) THEN
   
          !--------------------------------------------------------
          ! (6c,6d) We are in EURASIA
          !--------------------------------------------------------
          IF ( YMID > EAST_NS_DIV ) THEN


             ! 6d: Eurasian Mid-Latitudes
             TOTAL_IC = RFLASH_MIDLAT * RATE * ( 1d0 - X )
             TOTAL_CG = RFLASH_MIDLAT * RATE * X
          ELSE

             ! 6c: Eurasian Tropics
             TOTAL_IC = RFLASH_TROPIC * RATE * ( 1d0 - X )
             TOTAL_CG = RFLASH_TROPIC * RATE * X

          ENDIF

       ELSE

          !--------------------------------------------------------
          ! (6c,6d) We are in the AMERICAS
          !--------------------------------------------------------
          IF ( YMID > WEST_NS_DIV ) THEN

             ! 6d: American Mid-Latitudes
             TOTAL_IC = RFLASH_MIDLAT * RATE * ( 1d0 - X )
             TOTAL_CG = RFLASH_MIDLAT * RATE * X

          ELSE

             ! 6c: American Tropics
             TOTAL_IC = RFLASH_TROPIC * RATE * ( 1d0 - X )
             TOTAL_CG = RFLASH_TROPIC * RATE * X

          ENDIF
       ENDIF

       !-----------------------------------------------------------
       ! (6e) Compute total lightnox
       !-----------------------------------------------------------

       ! Sum of IC + CG [molec/6h]
       LNOX2D(I,J) = TOTAL_IC + TOTAL_CG

       ! Convert to flashes per minute per km2
       LFR2D(I,J) = RATE / A_KM2 / 360.0_hp

       ! IC flash rate
       IC2D(I,J)  = RATE * ( 1d0 - X ) / A_KM2 / 360.0_hp

       ! CG flash rates
       CG2D(I,J)  = RATE * X / A_KM2 / 360.0_hp

       ! Convective cloud top height in meters
       H02D(I,J) = H0

    ENDDO
    ENDDO
!$OMP END PARALLEL DO

    ! Return with success
    RC = 0

  END SUBROUTINE CalcFlashRate 
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Flashes_CTH
!
! !DESCRIPTION: Subroutine Flashes\_CTH determines the rate of lightnox 
!  flashes per minute based on the height of convective cloud tops, and the 
!  intra-cloud to cloud-ground strike ratio.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Flashes_CTH( I, J, HEIGHT, FLASHRATE, SFCTYPE ) 
!
! !INPUT PARAMETERS: 
!
    INTEGER, INTENT(IN)  :: I           ! Longitude index
    INTEGER, INTENT(IN)  :: J           ! Latitude index
    REAL*8,  INTENT(IN)  :: HEIGHT      ! Height of conv cloud top [m]
    INTEGER, INTENT(IN)  :: SFCTYPE     ! Surface type (0=land, 1=water, 2=ice)
!
! !OUTPUT PARAMETERS:
!
    REAL*8,  INTENT(OUT) :: FLASHRATE   ! LightNOX flash rate [flashes/min]
!
! !REVISION HISTORY: 
!  10 May 2006 - L. Murray - Initial version
!  (1  ) Subroutine renamed from FLASHES (ltm, bmy, 5/10/06)
!  (2  ) Remove CCTHICK, IC_CG_RATIO as arguments.  Remove computation of
!         IC_CG_RATIO and move that to GET_IC_CG_RATIO. (ltm, bmy, 12/11/06)
!  (3  ) Remove the near-land formulation (i.e. use function IS_LAND 
!         instead of IS_NEAR).(ltm, bmy, 9/24/07)
!  10 Nov 2010 - R. Yantosca - Added ProTeX headers
!  22 Oct 2013 - C. Keller   - Now a HEMCO extension.
!  10 Jun 2019 - C. Keller   - Moved to ligthning_toolbox_mod.F90
!EOP
!------------------------------------------------------------------------------
!BOC
    !================================================================
    ! FLASHES_CTH begins here!
    !
    ! COMPUTE LIGHTNOX FLASH RATE / MINUTE
    !
    ! Price & Rind (1992) give the following parameterizations for
    ! lightnox flash rates as a function of convective cloud top
    ! height [km]:
    !
    !    FLAND  = 3.44e-5 * ( CLDTOP HEIGHT [km] ^ 4.9  )
    !    FOCEAN = 6.4e-4  * ( CLDTOP HEIGHT [km] ^ 1.73 )
    !
    ! LightNOX will therefore occur much more often on land.  It 
    ! goes as approx. the 5th power of height, as opposed to approx. 
    ! the 2nd power of height over oceans.
    !
    ! We suppress lightnox where the surface is mostly ice.  
    !
    ! (ltm, bmy, 5/10/06, 12/11/06)
    !================================================================

    ! Test for land type
    IF ( SFCTYPE == 0 ) THEN

       ! Flashes/min over land boxes
       FLASHRATE   = 3.44d-5 * ( ( HEIGHT * 1d-3 )**4.9d0  )

    ELSE IF ( SFCTYPE == 1 ) THEN

       ! Flahes/min over water
       FLASHRATE   = 6.4d-4  * ( ( HEIGHT * 1d-3 )**1.73d0 )

    ELSE IF ( SFCTYPE == 2 ) THEN

       ! Suppress lightnox over snow/ice
       FLASHRATE   = 0d0

    ENDIF

  END SUBROUTINE Flashes_CTH
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_IC_CG_Ratio
!
! !DESCRIPTION: Function Get\_IC\_CG\_Ratio calculates the Intra-Cloud (IC) 
!  and Cloud-to-Ground (CG) lightnox flash ratio based on the method of 
!  Price and Rind 1993, which is calculated from the cold-cloud depth 
!  (CCTHICK).
!\\
!\\
! !INTERFACE:
!
  FUNCTION Get_IC_CG_ratio( CCTHICK ) RESULT( IC_CG_RATIO )
!
! !INPUT PARAMETERS: 
!
    REAL*8,  INTENT(IN) :: CCTHICK       ! Cold cloud thickness [m]
!
! !RETURN VALUE:
!
    REAL*8              :: IC_CG_RATIO   ! Intra-cloud/cloud-ground ratio
!
! !REVISION HISTORY: 
!  11 Dec 2006 - R. Yantosca - Initial version
!  (1 ) Split off from FLASHES_CTH, FLASHES_MFLUX, FLASHES_PRECON into this
!        separate function (ltm, bmy, 12/11/06)
!  (2 ) Bug fix for XLF compiler (morin, bmy, 7/8/09)
!  10 Nov 2010 - R. Yantosca - Added ProTeX headers
!  22 Oct 2013 - C. Keller   - Now a HEMCO extension.
!  10 Jun 2019 - C. Keller   - Moved to ligthning_toolbox_mod.F90
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL*8 :: CC, F_CG

    !=================================================================
    ! GET_IC_CG_RATIO begins here!
    !
    ! COMPUTE INTRA-CLOUD / CLOUD-GROUND FLASH RATIO 
    !
    ! Price & Rind (1993) compute the ratio of Cloud-Ground 
    ! to Total Flashes by the parameterization:
    !
    ! For 5.5 < dz < 14:
    !
    !     f_CG = 1 / (( A*dz^4 + B*dz^3 + C*dz^2 + D*dz + E ) + 1 )
    !
    ! For dz > 14:
    !
    !     f_CG = 0.02
    !                                        
    ! Where:
    !
    ! (1) dz is the depth [km] of the cloud above the freezing 
    !     level.  The cold-cloud thickness (dz) is the depth of 
    !     the layer between the cloud top and the center of the 
    !     highest layer for which the temperature exceeds 273 K. 
    !     The cold-cloud thickness is set to 5.5 km at grid points 
    !     where it is less than 5.5 km.
    !
    ! (2) The polynomial coefficients are:
    !        A=0.021,  B=-0.648,  C=7.493,  D=-36.54,  E=63.09
    !
    ! 
    ! Note: f_IC = 1 - f_CG
    ! 
    ! And hence,
    !
    !     IC_CG_RATIO = ( 1 - f_CG ) / f_CG
    !
    !
    ! IC_CG_RATIO is passed back to routine the LIGHTNOX_NL, where
    ! it is passed to FLASHES_MFLUX and FLASHES_PRECON.  In these
    ! routines, the fraction of total lightnox flashes that are 
    ! cloud-ground (CG) flashes is computed by:
    ! 
    !     F_CG        = 1d0 / ( 1d0 + IC_CG_RATIO )
    !
    ! and the fraction of the total lightnox flashes that are
    ! intra-cloud (IC) flashes is computed by:
    !
    !     F_IC        = 1d0 - 1d0 / ( 1d0 + IC_CG_RATIO )
    !=====================================================================

    ! Convert cold cloud thickness from [m] to [km] (min value: 5.5 km)
    CC = MAX( CCTHICK * 1d-3, 5.5d0 )

    ! Compute cloud-ground flash ratio as described above
    IF ( CC > 14d0 ) THEN

       ! Constant value above 14 km
       F_CG = 0.02d0

    ELSE

       ! First create the polynomial expression
       F_CG = 63.09d0 + CC * ( -36.54d0  + &
                        CC * (   7.493d0 + &
                        CC * (  -0.648d0 + &
                        CC * (   0.021d0 ) ) ) )

       ! Then put it in the denominator
       F_CG = 1d0 / ( F_CG + 1d0 )
                  
    ENDIF

    ! Intra-Cloud / Cloud-Ground flash ratio
    IC_CG_RATIO = ( 1d0 - F_CG ) / F_CG

  END FUNCTION Get_IC_CG_Ratio
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: LightDist
!
! !DESCRIPTION: Subroutine LightDist reads in the CDF used to partition the 
!  column lightnox NOx into the GEOS-Chem vertical layers. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE LightDist( NZ, NNLIGHT, NLTYPE, LTOP, H0, XLAT, TOTAL, &
                        VERTPROF, SFCTYPE, cMt, MTYPE, BXHEIGHT_M, PROFILE )
!
! !INPUT PARAMETERS: 
!
    INTEGER,         INTENT(IN)  :: NZ             ! Number of vertical levels 
    INTEGER,         INTENT(IN)  :: NNLIGHT        ! Dimension 1 of profile array
    INTEGER,         INTENT(IN)  :: NLTYPE         ! Dimension 2 of profile array 
    INTEGER,         INTENT(IN)  :: LTOP           ! Level of conv cloud top
    REAL*8,          INTENT(IN)  :: H0             ! Conv cloud top height [m]
    REAL*8,          INTENT(IN)  :: XLAT           ! Latitude value [degrees]
    REAL*8,          INTENT(IN)  :: TOTAL          ! Column Total # of LNOx molec 
    REAL(hp),        INTENT(IN)  :: BXHEIGHT_M(NZ) ! Boxheights in meter 
    INTEGER,         INTENT(IN)  :: SFCTYPE        ! Surface type 
    INTEGER,         INTENT(IN)  :: cMt            ! Current month 
    REAL*8,          INTENT(IN)  :: PROFILE(NNLIGHT, NLTYPE ) ! Profile tables
!
! !OUTPUT PARAMETERS:
!
    REAL*8,          INTENT(OUT) :: VERTPROF(NZ)   ! Vertical profile 
    INTEGER,         INTENT(OUT) :: MTYPE          ! lightning type 
!
! !REMARKS:
!  References:
!  ============================================================================
!  (1 ) Pickering et al., JGR 103, 31,203 - 31,316, 1998.
!  (2 ) Ott et al., JGR, 2010
!  (3 ) Allen et al., JGR, 2010
! 
! !REVISION HISTORY: 
!  18 Sep 2002 - M. Evans - Initial version (based on Yuhang Wang's code)
!  (1 ) Use functions IS_LAND and IS_WATER to determine if the given grid
!        box is over land or water.  These functions work for all DAO met
!        field data sets. (bmy, 4/2/02)
!  (2 ) Renamed M2 to LTOP and THEIGHT to H0 for consistency w/ variable names
!        w/in "lightnox.f".  Now read the "light_dist.dat.geos3" file for 
!        GEOS-3 directly from the DATA_DIR/lightnox_NOx_200203/ subdirectory.
!        Now read the "light_dist.dat" file for GEOS-1, GEOS-STRAT directly 
!        from the DATA_DIR/lightnox_NOx_200203/ subdirectory.  Added 
!        descriptive comment header.  Now trap I/O errors across all 
!        platforms with subroutine "ioerror.f".  Updated comments, cosmetic 
!        changes.  Redimension FRAC(NNLIGHT) to FRAC(LLPAR). (bmy, 4/2/02)
!  (3 ) Deleted obsolete code from April 2002.  Now reference IU_FILE and
!        IOERROR from "file_mod.f".  Now use IU_FILE instead of IUNIT as the
!        file unit number. (bmy, 6/27/02)
!  (4 ) Now reference BXHEIGHT from "dao_mod.f" (bmy, 9/18/02)
!  (5 ) Bug fix: add GEOS_4 to the #if block (bmy, 3/4/04)
!  (6 ) Now bundled into "lightnox_mod.f".  CDF's are now read w/in
!        routine INIT_LIGHTNOX to allow parallelization (bmy, 4/14/04)
!  (7 ) Now references DATA_DIR from "directory_mod.f" (bmy, 7/20/04)
!  (8 ) Now uses near-land formulation (ltm, bmy, 5/10/06)
!  (9 ) Added extra safety check for pathological boxes (bmy, 12/11/06)
!  (10) Remove the near-land formulation, except for PRECON (ltm, bmy, 9/24/07)
!  (11) Now use the Ott et al. [2010] profiles, and apply consistently with
!        GMI model [Allen et al., 2010] (ltm, bmy, 1/25/11).
!  10 Nov 2010 - R. Yantosca - Added ProTeX headers
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_CM2(I,J,L) from grid_mod.F90
!  15 Jun 2012 - Nielsen - INQUIRE finds free logical unit number for IU_FILE
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  22 Oct 2013 - C. Keller   - Now a HEMCO extension.
!  10 Jun 2019 - C. Keller   - Moved to lightning_toolbox_mod.F90.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: L
    REAL*8  :: ZHEIGHT, YMID
    REAL*8  :: FRAC(NZ)

      !=================================================================
      ! LIGHTDIST begins here!
      !=================================================================

    ! Initialize 
    MTYPE    = 0
    VERTPROF = 0d0

    !%%% NOTE: Use L=1 for GRID_MOD functions.  This is OK for the
    !%%% existing GEOS-Chem with a pure cartesian grid, but may be an
    !%%% issue when interfaced with a GCM with a non-regular grid
    !%%% (bmy, 3/1/12)
    YMID = XLAT 

    !=================================================================
    ! Test whether location (I,J) is continental, marine, or snow/ice
    !
    ! Depending on the combination of land/water and latitude, 
    ! assign a flag describing the type of lightnox:
    !
    !   MTYPE = 1: ocean lightnox
    !   MTYPE = 2: tropical continental lightnox
    !   MTYPE = 3: midlatitude continental lightnox 
    !   MTYPE = 4: subtropical lightnox
    !             
    ! (ltm, bmy, 1/25/11)
    !=================================================================

    ! Assign profile kind to grid box, following Allen et al. 
    ! [JGR, 2010] (ltm, 1/25,11)
!    MONTH = cMt

    SELECT CASE (cMt)

       ! Southern Hemisphere Summer
       CASE ( 1,2,3,12 )

           IF ( ABS(YMID) .le. 15 ) THEN
              IF ( SFCTYPE == 0 ) THEN
                 MTYPE = 2        ! Tropical continental
              ELSE
                 MTYPE = 1        ! Tropical marine
              ENDIF
           ELSE IF ( ( YMID .gt. 15. ) .and. ( YMID .le. 30. ) ) THEN
              MTYPE = 4           ! N. Subtropics
           ELSE IF ( ( YMID .ge. -40. ) .and. ( YMID .lt. -15. ) ) THEN
              MTYPE = 4           ! S. Subtropics
           ELSE
              MTYPE = 3           ! Midlatitude
           ENDIF

        ! Equinox months
        CASE ( 4,5,10,11 )

           IF ( ABS(YMID) .le. 15 ) THEN
              IF ( SFCTYPE == 0 ) THEN
                 MTYPE = 2        ! Tropical continental
              ELSE
                 MTYPE = 1        ! Tropical marine
              ENDIF
           ELSE IF ( ABS(YMID) .le. 30 ) THEN
              MTYPE = 4           ! Subtropics
           ELSE
              MTYPE = 3           ! Midlatitude
           ENDIF

        ! Northern Hemisphere Summer
        CASE ( 6,7,8,9 )

           IF ( ABS(YMID) .le. 15 ) THEN
              IF ( SFCTYPE == 0 ) THEN
                 MTYPE = 2        ! Tropical continental
              ELSE
                 MTYPE = 1        ! Tropical marine
              ENDIF
           ELSE IF ( ( YMID .gt. 15. ) .and. ( YMID .le. 40. ) ) THEN
              MTYPE = 4           ! N. Subtropics
           ELSE IF ( ( YMID .ge. -30. ) .and. ( YMID .lt. -15. ) ) THEN
              MTYPE = 4           ! S. Subtropics
           ELSE
              MTYPE = 3           ! Midlatitude
           ENDIF
         
    END SELECT

    ! Extra safety check for pathological grid boxes (bmy, 11/29/06)
    IF ( MTYPE == 0 ) RETURN

    !=================================================================
    ! Use the CDF for this type of lightnox to partition the total
    ! column lightnox into the layers
    !=================================================================
    ZHEIGHT = 0.0

    ! Compute the height [km] at the top of each vertical level.
    ! Look up the cumulative fraction of NOx for each vertical level
    DO L = 1, LTOP
       ZHEIGHT = ZHEIGHT + BXHEIGHT_M(L)
       FRAC(L) = PROFILE( NINT( ( ZHEIGHT/H0 )*3200. ), MTYPE ) *0.01
    ENDDO

    ! Convert from cumulative fraction to fraction for each level
    DO L = LTOP, 2, - 1
       FRAC(L) = FRAC(L) - FRAC(L-1)
    ENDDO
      
    ! Partition lightnox NOx by layer into VERTPROF
    DO L = 1, LTOP
       VERTPROF(L) = ( FRAC(L) * TOTAL )
    ENDDO

  END SUBROUTINE LightDist
!EOC
END MODULE Lightning_Toolbox_Mod 
