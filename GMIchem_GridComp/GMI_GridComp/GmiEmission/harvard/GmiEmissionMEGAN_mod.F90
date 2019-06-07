!-------------------------------------------------------------------------
!  NASA/GFSC, SIVO, Code 610.3
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: GmiEmissionMEGAN_mod
!
! !INTERFACE:
!
      module GmiEmissionMEGAN_mod
!
! !USES:

      implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
      private
      public  :: calcBiogenicMEGANemission
!
#    include "gmi_phys_constants.h"
!
! !DESCRIPTION:
!  Contains routines specifying the algorithms that control the MEGAN inventory 
!  of biogenic emissions.
!
! \begin{verbatim}
!  Variables needed by this module:
!
!     * surfTemp           : Current surf. temp. 
!     * T_15_AVG           : 24h average surf. temp. over past 15 days 
!     * aefIsop            : Annual emission factor for isoprene
!     * aefMonot           : Annual emission factor for monoterpenes
!     * aefMbo             : Annual emission factor for methyl butenol
!     * aefOvoc            : Annual emission factor for other biogenic VOCs
!     * pardif             : Diffuse photosynthetically active radiation (0.35-0.70 um)
!     * pardir             : Direct  photosynthetically active radiation (0.35-0.70 um)
!     * isoLai             : Daily value of LAI data (interpolation)
!     * isoLaiCurr         : AVHRR LAI data for the current  month
!     * isoLaiPrev         : AVHRR LAI data for the previous month
!     * DAYS_BTW_M         : days between midmonths in the LAI data
!     * cosSolarZenithAngle: cosine of the solar zenith angle
!
!  References:
!  ============================================================================
!  (1) Guenther, A., et al., A global model of natural volatile organic 
!      commpound emissions, J.Geophys. Res., 100, 8873-8892, 1995.
!  (2) Wang, Y., D. J. Jacob, and J. A. Logan, Global simulation of 
!      tropospheric O3-Nox-hydrocarbon chemistry: 1. Model formulation, J. 
!      Geophys. Res., 103, D9, 10713-10726, 1998.
!  (3) Guenther, A., B. Baugh, G. Brasseur, J. Greenberg, P. Harley, L. 
!      Klinger, D. Serca, and L. Vierling, Isoprene emission estimates and 
!      uncertanties for the Central African EXPRESSO study domain, J. 
!      Geophys. Res., 104, 30,625-30,639, 1999.
!  (4) Guenther, A. C., T. Pierce, B. Lamb, P. Harley, and R. Fall, Natural 
!      emissions of non-methane volatile organic compounds, carbon 
!      monoxide, and oxides of nitrogen from North America, Atmos. Environ.,
!      34, 2205-2230, 2000.
!  (5) Guenther, A., and C. Wiedinmyer, User's guide to Model of Emissions of 
!      Gases and Aerosols from Nature. http://cdp.ucar.edu. (Nov. 3, 2004) 
!  (6) Guenther, A., AEF for methyl butenol, personal commucation. (Nov, 2004)
! \end{verbatim}
!
! !DEFINED PARAMETERS:
      REAL*8,  PARAMETER  :: WM2_TO_UMOLM2S = 4.766d0
!
! !AUTHOR:
!  Bob Yantosca and Jules Kouatchou
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
      CONTAINS
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: calcBiogenicMEGANemission
!
! !INTERFACE:
!
      subroutine calcBiogenicMEGANemission  &
     &      (emiss_isop, emiss_monot, DAYS_BTW_M, &
     &       aefIsop, aefMbo, aefMonot, isoLai, isoLaiCurr, isoLaiPrev, &
     &       surfTemp, T_15_AVG, pardif, pardir, cosSolarZenithAngle, &
     &       pr_diag, procID, i1, i2, ju1, j2, PT_15isOK)
!
      implicit none

#     include "gmi_emiss_constants.h"
!
! !INPUT PARAMETERS:
      logical, intent(in) :: pr_diag
      integer, intent(in) :: procID
      integer, intent(in) :: i1, i2, ju1, j2
                             ! surface air temperature  (degK)
      real*8 , intent(in) :: surfTemp            (i1:i2, ju1:j2)
                             ! Daily surface temperature average over 15 days
      real*8 , intent(in) :: T_15_AVG            (i1:i2, ju1:j2)
                             ! Isoprene annual emission factors
      real*8 , intent(in) :: aefIsop             (i1:i2, ju1:j2)
                             ! Methyl butenol annual emission factors
      real*8 , intent(in) :: aefMbo              (i1:i2, ju1:j2)
                             ! Monoterpene annual emission factors
      real*8 , intent(in) :: aefMonot            (i1:i2, ju1:j2)
                             ! days between midmonths in the LAI data
      integer, intent(in) :: DAYS_BTW_M
                             ! Daily value of LAI data (interpolation)
      real*8 , intent(in) :: isoLai              (i1:i2, ju1:j2)
                             ! AVHRR LAI data for the current  month
      real*8 , intent(in) :: isoLaiCurr          (i1:i2, ju1:j2)
                             ! AVHRR LAI data for the previous month
      real*8 , intent(in) :: isoLaiPrev          (i1:i2, ju1:j2)
                             ! Diffuse photosynthetically active radiation (0.35-0.70 um)
      real*8 , intent(in) :: pardif              (i1:i2, ju1:j2)
                             ! Direct  photosynthetically active radiation (0.35-0.70 um)
      real*8 , intent(in) :: pardir              (i1:i2, ju1:j2)
                             ! Cosines of the solar zenith angle
      real*8 , intent(in) :: cosSolarZenithAngle (i1:i2, ju1:j2)
      LOGICAL, INTENT(IN) :: PT_15isOK ! If .TRUE. use Eq. 5b, 5c for E_OPT and T_OPT
!
! !OUTPUT PARAMETERS:
                              ! isoprene    emissions    (kg/s)
      real*8 , intent(out) :: emiss_isop (i1:i2, ju1:j2)
                              ! monoterpene emissions    (kg/s)
      real*8 , intent(out) :: emiss_monot(i1:i2, ju1:j2)
!
! !DESCRIPTION:
! Updates the MEGAN biogenic emissions (isoprene, monoterpene).
!
! !DEFINED PARAMETERS:
      ! Molecules C / kg C = { [ molec C /mole C ] / [ kg C / mole C ] }
      real*8, parameter :: XNUMOL_C = 6.022d+23 / 12d-3     
!
! !LOCAL VARIABLES:
      integer :: ix, iy
      real*8  :: Q_DIR, Q_DIFF, SUNCOS
      real*8  :: TS_ij, T_15_AVG_ij, aefIsop_ij, aefMbo_ij, aefMonot_ij
      real*8  :: isoLai_ij, isoLaiCurr_ij, isoLaiPrev_ij
!
! !AUTHOR:
! Jules Kouatchou, NASA GSFC, Jules.Kouatchou@nasa.gov
!
!
!EOP
!------------------------------------------------------------------------------
!BOC
      if (pr_diag) Write(6,*) ' calcBiogenicMEGANemission called by ', procID

      do iy = ju1, j2
         do ix = i1, i2

            Q_DIR         = pardir             (ix, iy)
            Q_DIFF        = pardif             (ix, iy)
            SUNCOS        = cosSolarZenithAngle(ix, iy)
            TS_ij         = surfTemp           (ix, iy)
            T_15_AVG_ij   = T_15_AVG           (ix, iy)
            aefIsop_ij    = aefIsop            (ix, iy)
            aefMbo_ij     = aefMbo             (ix, iy)
            aefMonot_ij   = aefMonot           (ix, iy)
            isoLai_ij     = isoLai             (ix, iy)
            isoLaiCurr_ij = isoLaiCurr         (ix, iy)
            isoLaiPrev_ij = isoLaiPrev         (ix, iy)

            ! Isoprene
            emiss_isop (ix,iy) =  GET_EMISOP_MEGAN &
     &                               (T_15_AVG_ij, SUNCOS, aefIsop_ij, &
     &                                isoLai_ij, isoLaiCurr_ij, isoLaiPrev_ij, &
     &                                DAYS_BTW_M, TS_ij, XNUMOL_C, Q_DIR, Q_DIFF, PT_15isOK)

            ! Monoterpenes
            emiss_monot(ix,iy) = GET_EMMONOT_MEGAN &
     &                              (T_15_AVG_ij, TS_ij, XNUMOL_C, aefMonot_ij, &
     &                               isoLai_ij, isoLaiCurr_ij, isoLaiPrev_ij, &
     &                               DAYS_BTW_M, PT_15isOK)

!            ! Methyl butenol
!            emiss_meth (ix,iy) = GET_EMMBO_MEGAN &
!     &                              (T_15_AVG_ij, SUNCOS, aefMbo_ij, &
!     &                               isoLai_ij, isoLaiCurr_ij, isoLaiPrev_ij, &
!     &                               DAYS_BTW_M, TS_ij, XNUMOL_C, Q_DIR, Q_DIFF, PT_15isOK )

         end do
      end do

      return

      end subroutine calcBiogenicMEGANemission
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GET_EMISOP_MEGAN
!
! !INTERFACE:
!
      FUNCTION GET_EMISOP_MEGAN(T_15_AVG_ij, SUNCOS, aefIsop_ij, &
     &                           isoLai_ij, isoLaiCurr_ij, isoLaiPrev_ij, &
     &                            DAYS_BTW_M, TS, XNUMOL, Q_DIR, Q_DIFF, PT_15isOK) &
     &         RESULT( EMISOP )

     implicit none
!
! !INPUT PARAMETERS:
      REAL*8,  INTENT(IN) :: SUNCOS      ! Cos( solar zenith angle )
      REAL*8,  INTENT(IN) :: TS          ! local surf air temperature [K]
      REAL*8,  INTENT(IN) :: XNUMOL      ! Number of atoms C / kg C
      REAL*8,  INTENT(IN) :: Q_DIR       ! Flux of direct PAR above canopy [W/m2]
      REAL*8,  INTENT(IN) :: Q_DIFF      ! Flux of diffuser PAR above canopy [W/m2]
      INTEGER, INTENT(IN) :: DAYS_BTW_M
      real*8,  intent(in) :: isoLai_ij
      real*8,  intent(in) :: isoLaiCurr_ij
      real*8,  intent(in) :: isoLaiPrev_ij
      real*8,  intent(in) :: aefIsop_ij
      real*8 , intent(in) :: T_15_AVG_ij
      LOGICAL, INTENT(IN) :: PT_15isOK ! If .TRUE. use Eq. 5b, 5c for E_OPT and T_OPT
!
! !RETURN VALUE:
      REAL*8              :: EMISOP
!
! !DESCRIPTION:
!  Computes ISOPRENE EMISSIONS in units of 
!  [atoms C/box] using the MEGAN inventory. (dsa, tmf, 9/03, 10/24/05)  
!
! \begin{verbatim}
!  References (see above for full citations):
!  ============================================================================
!  (1) Guenther et al, 1995
!  (2) Wang,    et al, 1998
!  (3) Guenther et al, 1999
!  (4) Guenther et al, 2000
!  (5) Guenther et al, 2004
!
!  Notes:
!  (1) Original code by Dorian Abbot (9/2003).  Updated to the latest 
!      algorithm and modified for the standard code by May Fu (11/20/04)
!  (2) All MEGAN biogenic emission are currently calculated using TS from DAO 
!      met field. TS is the surface air temperature, which should be 
!      carefully distinguished from TSKIN. (tmf, 11/20/04)
! \end{verbatim}
!
! !LOCAL VARIABLES:
      REAL*8              :: MEA,     DEA,     HEA_TL 
      REAL*8              :: D_BTW_M, Q_DIR_2, Q_DIFF_2
!
!
! !REVISION HISTORY:
!  Initial code.
!            
!EOP
!------------------------------------------------------------------------------
!BOC
      ! Initialize
      EMISOP   = 0.d0
      MEA      = 0.d0
      DEA      = 0.d0
      HEA_TL   = 0.d0

      ! Number of days between isoLaiCurr_ij and isoLaiPrev_ij
      D_BTW_M  = DBLE( DAYS_BTW_M )

      ! Convert Q_DIR and Q_DIFF from (W/m2) to (micromol/m2/s)
      Q_DIR_2  = Q_DIR  * WM2_TO_UMOLM2S
      Q_DIFF_2 = Q_DIFF * WM2_TO_UMOLM2S

      !---------------------------------------------------
      ! Do only during day and over continent
      ! Only interested in terrestrial biosphere (pip)
      !---------------------------------------------------
      IF ( SUNCOS > 0d0 ) THEN

         ! If (local LAI != 0 .AND. baseline emission !=0 ) 
         IF ( isoLai_ij * aefIsop_ij > 0d0 ) THEN

            ! Hourly exchange activitiy for temp & light
            HEA_TL = GET_HEA_TL( TS, T_15_AVG_ij, &
     &                           isoLai_ij, SUNCOS, Q_DIR_2, Q_DIFF_2, PT_15isOK)

            ! Daily exchange activity.  Alex Guenther advised us 
            ! to set DEA = 1.d0 for now (tmf, 10/24/05)
            DEA    = 1.d0

            ! Monthly exchange activity
            MEA    = GET_MEA( isoLaiCurr_ij, isoLaiPrev_ij, D_BTW_M )

         ELSE

            ! If it's night or over ocean, set activity factors to zero
            HEA_TL = 0.d0
            DEA    = 0.d0
            MEA    = 0.d0

         ENDIF
    
         ! Isoprene emission is the product of all these
         EMISOP    = aefIsop_ij * HEA_TL * DEA * MEA
!%
!% NOTE: you can convert to whichever unit the GMI model needs.
!% GEOS-Chem requires atoms C/box.
!%
         ! CONVERT from [kg/box] to [atoms C/box]
         EMISOP    = EMISOP * XNUMOL
         
      ENDIF

      ! return to calling program
      END FUNCTION GET_EMISOP_MEGAN

!------------------------------------------------------------------------------
!BOC
!
! !IROUTINE: GET_EMMONOT_MEGAN
!
! !INTERFACE:
!
      FUNCTION GET_EMMONOT_MEGAN(T_15_AVG_ij, TS, XNUMOL, aefMonot_ij, &
     &                           isoLai_ij, isoLaiCurr_ij, isoLaiPrev_ij, DAYS_BTW_M, PT_15isOK) &
     &                           RESULT( EMMONOT )

      implicit none
!
! !INPUT PARAMETERS:
      REAL*8,  INTENT(IN) :: TS          ! local surf air temperature [K]
      REAL*8,  INTENT(IN) :: XNUMOL      ! Number of atoms C / kg C
      INTEGER, INTENT(IN) :: DAYS_BTW_M
      real*8,  intent(in) :: aefMonot_ij
      real*8,  intent(in) :: isoLai_ij
      real*8,  intent(in) :: isoLaiCurr_ij
      real*8,  intent(in) :: isoLaiPrev_ij
      real*8,  intent(in) :: T_15_AVG_ij
      LOGICAL, INTENT(IN) :: PT_15isOK ! If .TRUE. use Eq. 5b, 5c for E_OPT and T_OPT
!
! !RETURN VALUE:
      REAL*8              :: EMMONOT
!
! !DESCRIPTION:
!  Subroutine GET_EMMONOT_MEGAN computes MONOTERPENE EMISSIONS in units of 
!  [atoms C/box] using the MEGAN inventory. (dsa, tmf, 9/03, 11/20/04)  
!
! \begin{verbatim}
!  References (see above for full citations):
!  ============================================================================
!  (1) Guenther et al, 1995
!  (2) Wang,    et al, 1998
!  (3) Guenther et al, 1999
!  (4) Guenther et al, 2000
!  (5) Guenther et al, 2004
!
!  Notes:
!  (1) Original code by Dorian Abbot (9/2003).  Updated to the latest 
!      algorithm and modified for the standard code by May Fu (11/20/04)
!  (2) All MEGAN biogenic emission are currently calculated using TS from DAO 
!      met field. TS is the surface air temperature, which should be 
!      carefully distinguished from TSKIN. (tmf, 11/20/04)
! \end{verbatim}
!
! !LOCAL VARIABLES:
      REAL*8              :: MEA, DEA, HEA_T, D_BTW_M
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC
      ! Initialize
      EMMONOT  = 0.d0
      MEA      = 0.d0
      DEA      = 0.d0
      HEA_T    = 0.d0

      ! Number of days between isoLaiCurr_ij and isoLaiPrev_ij
      D_BTW_M  = DBLE( DAYS_BTW_M )

      !-----------------------------------------------------
      ! Only interested in terrestrial biosphere (pip)
      ! If (local LAI != 0 .AND. baseline emission !=0 ) 
      !-----------------------------------------------------
      IF ( isoLai_ij * aefMonot_ij > 0d0 ) THEN

         ! Hourly exchange activity (temp only)
         HEA_T = GET_HEA_T( TS, T_15_AVG_ij, PT_15isOK ) 

         ! Daily exchange activity.  Alex Guenther advised us to 
         ! set DEA = 1.d0 for now (tmf, 10/24/05)
         DEA   = 1.d0

         ! Monthly exchange activity
         MEA   = GET_MEA( isoLaiCurr_ij, isoLaiPrev_ij, D_BTW_M )

      ELSE

         ! Otherwise set all activities to zero
         HEA_T = 0.d0
         DEA   = 0.d0
         MEA   = 0.d0

      ENDIF    

      ! Monoterpene emissions [kg/box]
      EMMONOT  = aefMonot_ij * HEA_T * DEA * MEA
!%
!% NOTE: you can convert to whichever unit the GMI model needs.
!% GEOS-Chem requires atoms C/box.
!%
      ! Convert from [kg/box] to [atoms C/box]
      EMMONOT  = EMMONOT * XNUMOL

      END FUNCTION GET_EMMONOT_MEGAN
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GET_EMMBO_MEGAN
!
! !INTERFACE:
!
      FUNCTION GET_EMMBO_MEGAN(T_15_AVG_ij, SUNCOS, aefMbo_ij, &
     &                          isoLai_ij, isoLaiCurr_ij, isoLaiPrev_ij, DAYS_BTW_M, & 
     &                          TS, XNUMOL, Q_DIR, Q_DIFF, PT_15isOK)  &
     &         RESULT( EMMBO )

      implicit none
!
! !INPUT PARAMETERS:
      REAL*8,  INTENT(IN) :: TS          ! local surf air temperature [K]
      REAL*8,  INTENT(IN) :: XNUMOL      ! Number of atoms C / kg C
      REAL*8,  INTENT(IN) :: SUNCOS      ! Cos( solar zenith angle )
      REAL*8,  INTENT(IN) :: Q_DIR       ! Flux of direct PAR above canopy [W/m2]
      REAL*8,  INTENT(IN) :: Q_DIFF      ! Flux of diffuser PAR above canopy [W/m2]
      INTEGER, INTENT(IN) :: DAYS_BTW_M
      real*8,  intent(in) :: aefMbo_ij
      real*8,  intent(in) :: isoLai_ij
      real*8,  intent(in) :: isoLaiCurr_ij
      real*8,  intent(in) :: isoLaiPrev_ij
      real*8,  intent(in) :: T_15_AVG_ij
      LOGICAL, INTENT(IN) :: PT_15isOK ! If .TRUE. use Eq. 5b, 5c for E_OPT and T_OPT
!
! !RETURN VALUE:
      REAL*8              :: EMMBO
!
! !DESCRIPTION:
!  Computes METHYLBUTENOL EMISSIONS in units of 
!  [atoms C/box] using the MEGAN inventory. (dsa, tmf, bmy, 9/03, 10/24/05)  
!
! \begin{verbatim}
!  References (see above for full citations):
!  ============================================================================
!  (1) Guenther et al, 1995
!  (2) Wang     et al, 1998
!  (3) Guenther et al, 1999
!  (4) Guenther et al, 2000
!  (5) Guenther et al, 2004
! \end{verbatim}
!
! !LOCAL VARIABLES:
      REAL*8              :: MEA,     DEA,     HEA_TL
      REAL*8              :: D_BTW_M, Q_DIR_2, Q_DIFF_2
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC
      ! Initialize
      EMMBO    = 0.d0
      MEA      = 0.d0
      DEA      = 0.d0
      HEA_TL   = 0.d0

      ! Number of days between isoLaiCurr_ij and isoLaiPrev
      D_BTW_M  = DBLE( DAYS_BTW_M )

      ! Convert Q_DIR and Q_DIFF from [W/m2] to [umol/m2/s]
      Q_DIR_2  = Q_DIR  * WM2_TO_UMOLM2S
      Q_DIFF_2 = Q_DIFF * WM2_TO_UMOLM2S

      !-------------------------------------------------
      ! Do only during day and over continent
      ! Only interested in terrestrial biosphere (pip)
      !-------------------------------------------------
      IF ( SUNCOS > 0d0 ) THEN

         ! If local LAI > 0 and baseline emission > 0 ...
         IF ( isoLai_ij * aefMbo_ij > 0d0 ) THEN

            ! Hourly exchange activity (based on temp & light)
            HEA_TL    = GET_HEA_TL( TS, T_15_AVG_ij, isoLai_ij, SUNCOS, &
     &                              Q_DIR_2, Q_DIFF_2, PT_15isOK )

            ! Daily exchange activity.  Alex Guenther advised us 
            ! to set DEA = 1.d0 for now (tmf, 11/20/05)
            DEA    = 1.d0

            ! Monthly exchange activity
            MEA    = GET_MEA( isoLaiCurr_ij, isoLaiPrev_ij, D_BTW_M )
       
         ELSE

            ! Otherwise set activities to zero
            HEA_TL = 0.d0
            DEA    = 0.d0
            MEA    = 0.d0

         ENDIF

         ! MBO emissions in [kg/box]
         EMMBO = aefMbo_ij * HEA_TL * DEA * MEA
!%
!% NOTE: you can convert to whichever unit the GMI model needs.
!% GEOS-Chem requires atoms C/box.
!%
         ! Convert from [atoms C/box] to [kg/box]
         EMMBO = EMMBO * XNUMOL
         
      ENDIF

      END FUNCTION GET_EMMBO_MEGAN
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GET_EMOVOC_MEGAN
!
! !INTERFACE:
!
      FUNCTION GET_EMOVOC_MEGAN(T_15_AVG_ij, TS, XNUMOL, aefOvoc_ij, &
     &                           isoLai_ij, isoLaiCurr_ij, isoLaiPrev_ij, DAYS_BTW_M, PT_15isOK) &
     &                          RESULT( EMOVOC )

      implicit none
!
! !INPUT PARAMETERS:
      REAL*8,  INTENT(IN) :: TS          ! local surf air temperature [K]
      REAL*8,  INTENT(IN) :: XNUMOL      ! Number of atoms C / kg C
      INTEGER, INTENT(IN) :: DAYS_BTW_M
      real*8,  intent(in) :: aefOvoc_ij
      real*8,  intent(in) :: isoLai_ij
      real*8,  intent(in) :: isoLaiCurr_ij
      real*8,  intent(in) :: isoLaiPrev_ij
      real*8,  intent(in) :: T_15_AVG_ij
      LOGICAL, INTENT(IN) :: PT_15isOK ! If .TRUE. use Eq. 5b, 5c mods for E_OPT and T_OPT
!
! !RETURN VALUE:
      REAL*8              :: EMOVOC
!
! !DESCRIPTION:
!  Computes other BVOC EMISSIONS in units of 
!  [atoms C/box] using the MEGAN inventory. (dsa, tmf, 9/03, 11/20/04)  
!
! \begin{verbatim}
!  References (see above for full citations):
!  ============================================================================
!  (1) Guenther et al, 1995
!  (2) Wang     et al, 1998
!  (3) Guenther et al, 1999
!  (4) Guenther et al, 2000
!  (5) Guenther et al, 2004
!  (6) Guenther pers.comm, 2004
! \end{verbatim}
!
! !LOCAL VARIABLES:
      REAL*8              :: MEA, DEA, HEA_T, D_BTW_M
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC
      ! Initialize
      EMOVOC  = 0.d0
      MEA     = 0.d0
      DEA     = 0.d0
      HEA_T   = 0.d0

      ! Number of days between isoLaiCurr_ij and isoLaiPrev
      D_BTW_M = DBLE( DAYS_BTW_M )

      !--------------------------------------------------
      ! Only interested in terrestrial biosphere (pip)
      ! If (local LAI != 0 .AND. baseline emission !=0 ) 
      !--------------------------------------------------
      IF ( isoLai_ij * aefOvoc_ij > 0d0 ) THEN

         ! Hourly exchange activity (temp only)
         HEA_T = GET_HEA_T( TS, T_15_AVG_ij, PT_15isOK ) 

         ! Daily exchange activity.  Alex Guenther advised us 
         ! to set DEA = 1.d0 for now (tmf, 10/24/05)
         DEA   = 1.d0

         ! Monthly exchange activity
         MEA   = GET_MEA( isoLaiCurr_ij, isoLaiPrev_ij, D_BTW_M )
       
      ELSE
         
         ! Otherwise set activities to zero
         HEA_T = 0.d0
         DEA   = 0.d0
         MEA   = 0.d0

      ENDIF

      ! OVOC emissions [kg/box]
      EMOVOC = aefOvoc_ij * HEA_T * DEA * MEA
!%
!% NOTE: you can convert to whichever unit the GMI model needs.
!% GEOS-Chem requires atoms C/box.
!%
      ! Convert from [kg/box] to [atoms C/box]
      EMOVOC = EMOVOC * XNUMOL

      ! return to calling program
      END FUNCTION GET_EMOVOC_MEGAN
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GET_HEA_TL
!
! !INTERFACE:
!
      FUNCTION GET_HEA_TL( T,       PT_15,   LAI, &
     &                     SUNCOS1, Q_DIR_2, Q_DIFF_2 ,PT_15isOK) RESULT( HEA_TL )

      implicit none 
!
! !INPUT PARAMETERS:
      REAL*8,  INTENT(IN) :: T         ! Current leaf temperature, the surface air temp
                                       ! field (TS) is assumed equivalent to the leaf
                                       ! temperature over forests.
      REAL*8,  INTENT(IN) :: PT_15     ! Average leaf tempearture over the past 15 days
      REAL*8,  INTENT(IN) :: LAI       ! cumulative leaf area index above leaf
      REAL*8,  INTENT(IN) :: SUNCOS1   ! Cosine of solar zenith angle
      REAL*8,  INTENT(IN) :: Q_DIR_2   ! flux of direct PAR above canopy [umol m-2 s-1]
      REAL*8,  INTENT(IN) :: Q_DIFF_2  ! flux of diffuser PAR above canopy [umol m-2 s-1]
      LOGICAL, INTENT(IN) :: PT_15isOK ! If .TRUE. use Eq. 5b, 5c mods for E_OPT and T_OPT
!
! !RETURN VALUE:
      REAL*8              :: HEA_TL
!
! !DESCRIPTION:
!  Computes the hourly exchange activity factor ({\em HEA\_TL}) with sensitivity to 
!  both TEMPERATURE and LIGHT for ISOP emission.  (tmf, 11/18/04, 10/24/05)
!
! $$ HEA\_TL = C\_T * C\_PPFD $$
! %
! \begin{description}   
!  \item[C\_T:] Effect of temperature on leaf isoprene emission, including effect 
!          of average temperature over previous 15 days, based on Eq 5abc 
!          from Guenther et al. (1999)
!
!  \item[C\_PPFD:] Effect of increasing PPFD up to a saturation point, where emission 
!          level off, based on Eq 4abc from Guenther et al. (1999)
!          In addition, a 5 layered canopy model based on Eqs 12-16 
!          from Guenther et al. (1995) is included to correct for light 
!          attenuation in the canopy.
! \end{description}   
!
! \begin{verbatim}
!  References (see above for full citations):
!  ============================================================================
!  (1) Guenther et al, 1995
!  (2) Wang     et al, 1998
!  (3) Guenther et al, 1999
!  (5) Guenther et al, 2004
! \end{verbatim}
!
! !LOCAL VARIABLES:
      !-----------------------------------------------------------------
      ! Based on Eq 5a, 5b, and 5c in Guenther et al. (1999)
      !-----------------------------------------------------------------

      ! C_T: Effect of temperature on leaf BVOC emission, including 
      ! effect of average temperature over previous 15 days, based on 
      ! Eq 5a, 5b, 5c from Guenther et al, 1999.
      REAL*8              :: C_T

      ! E_OPT: maximum normalized emission capacity
      REAL*8              :: E_OPT

      ! T_OPT: temperature at which E_OPT occurs
      REAL*8              :: T_OPT

      ! X: variable related to temp.
      REAL*8              :: X

      ! CT1, CT2: energy of activation and deactivation, respectively 
      REAL*8              :: CT1, CT2

      ! R - ideal gas constant (J mol-1 K-1)
      REAL*8, PARAMETER   :: R = 8.314d-3


      !-----------------------------------------------------------------
      ! Canopy model variables (Eqs 12-16 from Guenther et al, 1995)
      !-----------------------------------------------------------------

      ! C_PPFD: Effect of increasing PPFD up to a saturation point, where 
      ! emissions level off, based on Eq 4abc from Guenther et al. (1999)
      ! In addition, a 5 layered canopy model based on Eqs 12-16 
      ! from Guenther et al. (1995) is included to correct for light 
      ! attenuation in the canopy.
      REAL*8             :: C_PPFD

      ! LAYERS: number of layers in canopy model
      INTEGER, PARAMETER :: LAYERS = 5

      ! A: mean leaf-Sun angle, 60 degrees represents a 
      !    spherical leaf angle distribution
      REAL*8,  PARAMETER  :: A = 6.0d1 * GMI_PI / 1.8d2 

      ! WEIGHTGAUSS : weights for gaussian integration
      REAL*8,  PARAMETER :: WEIGHTGAUSS(LAYERS) = (/ 1.1846d-1,  &
     &                                                2.3931d-1,  &
     &                                                2.8444d-1,  &
     &                                                2.3931d-1,  &
     &                                                1.1846d-1 /)

      ! DISTGAUSS: points to evaluate fcn for gaussian integration
      REAL*8,  PARAMETER :: DISTGAUSS(LAYERS)   = (/ 4.6910d-2,  &
     &                                                2.3075d-1,  &
     &                                                5.0000d-1,  &
     &                                                7.6924d-1,  &
     &                                                9.5308d-1 /)

      ! SCAT: Scattering coefficient
      REAL*8,  PARAMETER :: SCAT = 2.0d-1

      ! REFLD: Reflection coefficient for diffuse light
      REAL*8,  PARAMETER :: REFLD = 5.7d-2

      ! CLUSTER: clustering coefficient (accounts for leaf clumping 
      ! influence on mean projected leaf area in the direction 
      ! of the sun's beam) use 0.85 for default
      REAL*8, PARAMETER   :: CLUSTER = 8.5d-1

      ! NORMAL_FACTOR : C_PPFD calculated with LAI = 5, and above canopy 
      ! total PPFD = 1500 umol/m2/s, and Q_DIFF = 0.2 * above canopy total 
      ! PPFD.  May Fu calculated this to be 1.8967d0
      REAL*8,  PARAMETER :: NORMAL_FACTOR = 1.8967d0

      REAL*8             :: F_SUN     ! Fraction leaves sunlit
      REAL*8             :: F_SHADE   ! Fraction leaves shaded
      REAL*8             :: LAI_DEPTH ! Cumulative LAI above current layer
      REAL*8             :: Q_SUN     ! Flux density of PAR on sunny leaves [umol/m2/s]
      REAL*8             :: Q_SHADE   ! Flux density of PAR on shaded leaves [umol/m2/s]
      REAL*8             :: Q_SHADE_1
      REAL*8             :: Q_SHADE_2
      REAL*8             :: KB        ! Extinction coefficient for black 
                                      ! leaves for direct (beam)
      REAL*8             :: KBP
      REAL*8             :: KD        ! Extinction coefficient for black 
                                      ! leaves for diffuse light
      REAL*8             :: KDP
      REAL*8             :: REFLB
      REAL*8             :: C_P_SUN   ! C_PPFD at layer I for sunlit leaves
      REAL*8             :: C_P_SHADE ! C_PPFD at layer I for shaded leaves

      REAL*8             :: C_PPFD_I  ! C_PPFD at layer I
      REAL*8             :: ALPHA, CL ! Empirical functions 
                                      ! (Eq 4a, 4b, 4c from Guenther et al, 1999)
      REAL*8             :: P
      INTEGER            :: I ! Index for layers
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC
      !-----------------------------------------------------------------
      ! Compute C_T
      !
      ! NOTES:
      ! (1) Algorithm is from Eq 5a, 5b, and 5c in Guenther et al, 1999
      ! (2) Alex Guenther suggests that we set DEA=1 for now
      ! (3) E_OPT and T_OPT depend on PT_15 according to Eq 5b and 5c
      !-----------------------------------------------------------------

      ! E_OPT: maximum normalized emission capacity

      IF(PT_15isOK) THEN
       E_OPT = 1.9d0 * EXP( 1.25d-1 * ( PT_15 - 3.01d2 ) )
      ELSE
       E_OPT = 1.9d0
      END IF

      ! T_OPT: temperature at which E_OPT occurs

      IF(PT_15isOK) THEN
       T_OPT = 3.16d2 + ( 5.0d-1 * ( PT_15 - 3.01d2 ) )
      ELSE
       T_OPT = 3.16d2
      END IF

      ! Energy of activation 
      CT1   = 76d0

      ! Energy of deactivation 
      CT2   = 160d0

      ! Variable related to temperature 
      X     = ( 1.d0/T_OPT - 1.d0/T ) / R

      ! C_T: Effect of temperature on leaf BVOC emission, including 
      ! effect of average temperature over previous 15 days, based on 
      ! Eq 5a, 5b, 5c from Guenther et al, 1999.
      C_T   = E_OPT * CT2 * EXP( CT1 * X ) /  &
     &        ( CT2 - CT1 * ( 1.d0 - EXP( CT2 * X ) ) )

      !-----------------------------------------------------------------
      ! Compute C_PPFD
      !-----------------------------------------------------------------  

      ! Initialize
      C_PPFD  = 0.d0

      ! 0.5 and 0.8 assume a spherical leaf-angle distribution
      KB      = 0.5d0 * CLUSTER / SUNCOS1
      KD      = 0.8d0 * CLUSTER
      P       = SQRT( 1.d0 - SCAT )
      REFLB   = 1.d0 -  &
     &          EXP(-2.d0*KB * ( 1.d0-P ) / ( 1.d0+P ) / ( 1.d0+KB ) )
      KBP     = KB * P
      KDP     = KD * P

      ! 5-layer Gaussian integration over canopy
      DO I = 1, LAYERS 

         ! Cumulative LAI above layer I
         LAI_DEPTH  = LAI * DISTGAUSS( I )

         ! Fraction sun and shade leaves at layer I
         F_SUN      = EXP( -KB * LAI_DEPTH ) 
         F_SHADE    = 1.d0 - F_SUN

         ! For PAR on shaded leaves
         Q_SHADE_1  = Q_DIFF_2 * KDP * ( 1.d0 - REFLD ) * &
     &                EXP( -KDP * LAI_DEPTH ) 

         ! For PAR on shaded leaves
         Q_SHADE_2  = Q_DIR_2 * ( KBP * ( 1.d0 - REFLB ) *  &
     &                EXP( -KBP * LAI_DEPTH ) - KB * ( 1.d0 - SCAT ) *  &
     &                EXP( -KB * LAI_DEPTH ) )

         ! PAR on shaded leaves
         Q_SHADE    = ( Q_SHADE_1 + Q_SHADE_2 ) / ( 1.d0 - SCAT )

         ! PAR on sunlit leaves
         Q_SUN      = Q_SHADE + KB * Q_DIR_2

         ! Update C_P_SUN and C_P_SHADE at layer I
         ! The following already accounts for canopy attenuation
         ! (Guenther et al, 1999)
         ALPHA      = 1.0d-3 + ( 8.5d-4 ) * LAI_DEPTH
         CL         = ( 1.42d0 ) * EXP( -( 3.0d-1 ) * LAI_DEPTH )

         C_P_SUN    = ALPHA * CL * Q_SUN / &
     &                SQRT( 1.d0 + ALPHA*ALPHA * Q_SUN*Q_SUN )

         C_P_SHADE  = ALPHA * CL * Q_SHADE / &
     &                SQRT( 1.d0 + ALPHA*ALPHA * Q_SHADE*Q_SHADE )

         ! Update C_PPFD_I at layer I
         C_PPFD_I   = F_SUN * C_P_SUN + F_SHADE * C_P_SHADE

         ! Add on to total C_PPFD
         C_PPFD     = C_PPFD + WEIGHTGAUSS( I ) * C_PPFD_I * LAI

      ENDDO

      ! Finally, HEA_TL.  Prevent negative values.
      HEA_TL = MAX( C_T * C_PPFD / NORMAL_FACTOR, 0d0 )

      ! return to calling program
      END FUNCTION GET_HEA_TL
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GET_HEA_T
!
! !INTERFACE:
!
      FUNCTION GET_HEA_T( T, PT_15, PT_15isOK) RESULT( HEA_T )

      implicit none
!
! !INPUT PARAMETERS:
      REAL*8,  INTENT(IN) :: T       ! Current leaf temperature, the surface air
                                     ! temperature field (TS) is assumed equivalent to
                                     ! the leaf temperature over forests.
      REAL*8,  INTENT(IN) :: PT_15   ! Average leaf temperature over the past 15 days
      LOGICAL, INTENT(IN) :: PT_15isOK ! If .TRUE. use Eq. 5b, 5c for E_OPT and T_OPT
!
! !RETURN VALUE:
      REAL*8              :: HEA_T
!
! !DESCRIPTION:
!  Computes the hourly exchange activity factor ({\em HEA\_T}) with sensitivity 
!  to only TEMPERATURE for all BVOC emission except ISOP. 
!  (tmf, bmy, 11/25/04, 10/24/05)
!
! \begin{verbatim}
!  References (see above for full citations):
!  ============================================================================
!  (1) Guenther et al, 1995
!  (2) Wang     et al, 1998
!  (3) Guenther et al, 1999
!  (5) Guenther et al, 2004
!
!  Notes:
!  (1) Original code by Dorian Abbot (9/2003).  Updated to the latest 
!      algorithm and modified for the standard code by May Fu (11/2004)
!  (2) Algorithm is based on Guenther (1999) instead of the latest User's 
!      Guide, since the User's Guide does not contain sensitivity to !
!      temperature in the past 15 days. (tmf, 11/19/04)
!  (3) Algorithm updated to Users' Guide as of Nov 3, 2004.  Sensitivity on 
!      temperature history is included in DEA from GET_DEA. (tmf, 11/19/04)
!  (4) All MEGAN biogenic emission are currently calculated using TS from DAO 
!      met field. TS is the surface air temperature, which should be 
!      carefully distinguished from TSKIN. (tmf, 11/20/04)
!  (5) Switch back to Guenther et al. (1999) algorithm and account for !
!      tempearture history in HEA.  This is accompanied by setting DEA = 1. 
!      (Alex Guenther, personal communication, 11/25/04)
!  (6) For consistency, C_T is calculated according to Guenther et al. (1999) 
!      algorithm.  Alex Guenther suggested that we use HEA TYPE 2 algorithm 
!      in Eq 12 from User's guide.  However, this makes incorporation of 
!      temperature history more difficult.  tmf checked the difference 
!      between Eq 9 and 12 from User's guide and they are essentially the 
!      same. Therefore we see no need to differentiate the C_T in HEA_T 
!      and HEA_TL.
! \end{verbatim}
!
! !LOCAL VARIABLES:
      REAL*8              :: C_T,   CT1,   CT2
      REAL*8              :: E_OPT, T_OPT, X

      ! Ideal gas constant [J/mol/K]
      REAL*8, PARAMETER   :: R   = 8.314d-3
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC
      !=================================================================
      ! NOTES:
      ! (1) Algorithm is from Eq 5a, 5b, and 5c in Guenther et al, 1999
      ! (2) Alex Guenther suggests that we set DEA=1 for now
      ! (3) E_OPT and T_OPT depend on PT_15 according to Eq 5b and 5c
      !================================================================= 

      ! E_OPT: maximum normalized emission capacity

      IF(PT_15isOK) THEN
       E_OPT = 1.9d0 * EXP( 1.25d-1 * ( PT_15 - 3.01d2 ) )
      ELSE
       E_OPT = 1.9d0
      END IF

      ! T_OPT: temperature at which E_OPT occurs

      IF(PT_15isOK) THEN
       T_OPT = 3.16d2 + ( 5.0d-1 * ( PT_15 - 3.01d2 ) )
      ELSE
       T_OPT = 3.16d2
      END IF

      ! Energy of activation 
      CT1   = 76d0

      ! Energy of deactivation 
      CT2   = 160d0

      ! Variable related to temperature 
      X     = ( 1.d0/T_OPT - 1.d0/T ) / R

      ! C_T: Effect of temperature on leaf BVOC emission, including 
      ! effect of average temperature over previous 15 days, based on 
      ! Eq 5a, 5b, 5c from Guenther et al, 1999.
      C_T   = E_OPT * CT2 * EXP( CT1 * X ) /  &
     &        ( CT2 - CT1 * ( 1.d0 - EXP( CT2 * X ) ) )

      ! Hourly emission activity = C_T
      ! Prevent negative values
      HEA_T = MAX( C_T, 0d0 )

      ! Return to calling program
      END FUNCTION GET_HEA_T
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GET_MEA
!
! !INTERFACE:
!
      FUNCTION GET_MEA( CMLAI, PMLAI, T ) RESULT( MEA )

      implicit none
!
! !INPUT PARAMETERS:
      REAL*8, INTENT(IN) :: T       ! Number of days between current and previous LAI
      REAL*8, INTENT(IN) :: CMLAI   ! Current month leaf area index at gridbox
      REAL*8, INTENT(IN) :: PMLAI   ! Last month leaf area index at gridbox
!
! !RETURN VALUE:
      REAL*8             :: MEA        
!
! !DESCRIPTION:
!  Computes the monthly exchange activity factor (MEA). (tmf, 10/24/05)
!
!     $$ MEA = M\_LAI \times M\_AGE \times M\_H $$
!  
! \begin{verbatim}
!  References (see above for full citations):
!  ============================================================================
!  (3) Guenther et al, 1999
!  (5) Guenther et al, 2004
! \end{verbatim}
!
!  !REMARKS:
!  (1) Original code by Dorian Abbot (9/2003).  Modified for the standard 
!      code by May Fu (11/2004)
!  (2) Update to publically released (as of 11/2004) MEGAN algorithm and 
!      modified for the standard code by May Fu (11/2004).
!  (3) Algorithm is based on the latest User's Guide (tmf, 11/19/04)
!
! !LOCAL VARIABLES:
      ! M_LAI: leaf area factor
      REAL*8             :: M_LAI 

      ! M_AGE: leaf age factor
      REAL*8             :: M_AGE

      ! M_H: sensible heat flux factor      
      REAL*8             :: M_H

      ! FNEW, ANEW: new foliage that emits negligible amounts of isoprene
      REAL*8             :: FNEW
      REAL*8,  PARAMETER :: ANEW = 1.d-2

      ! FGRO, AGRO: growing foliage that emits isoprene at low rates
      REAL*8             :: FGRO
      REAL*8,  PARAMETER :: AGRO = 5.d-1

      ! FMAT, AMAT: mature foliage that emits isoprene at peak rates
      REAL*8             :: FMAT
      REAL*8,  PARAMETER :: AMAT = 1.d0

      ! FSEN, ASEN: senescing foliage that emits isoprene at reduced rates
      REAL*8             :: FSEN
      REAL*8,  PARAMETER :: ASEN = 3.3d-1

      ! TI: number of days after budbreak required to induce iso. em.
      REAL*8,  PARAMETER :: TI   = 12.d0

      ! TM: number of days after budbreak required to reach peak iso. em. rates
      REAL*8,  PARAMETER :: TM   = 28.d0

      ! Variable for storing T or TM
      REAL*8             :: TG
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC
      !-----------------------
      ! Compute M_LAI
      !-----------------------
      M_LAI = 0.49d0 * CMLAI / SQRT( 1.d0 + 0.2d0 * CMLAI*CMLAI )
      
      !-----------------------
      ! Compute M_AGE
      !-----------------------
      IF ( T > TM ) THEN 
         TG = TM 
      ELSE 
         TG = T
      ENDIF

      IF ( CMLAI == PMLAI ) THEN

         FMAT = 1.d0  
         FNEW = 0.d0
         FGRO = 0.d0
         FSEN = 0.d0

      ELSE IF ( CMLAI > PMLAI ) THEN

         FSEN = 0.d0

         IF ( T > TI ) THEN
            FNEW = ( TI / T ) * ( 1.d0 -  PMLAI / CMLAI )
            FGRO = ( ( TG - TI ) / T ) * ( 1.d0 -  PMLAI / CMLAI )
         ELSE
            FNEW = 1.d0 - ( PMLAI / CMLAI )
            FGRO = 0.d0
         ENDIF

         IF ( T > TM ) THEN
            FMAT = ( PMLAI / CMLAI ) + &
     &             ( ( T - TM ) / T ) * ( 1.d0 -  PMLAI / CMLAI )
         ELSE 
            FMAT = ( PMLAI / CMLAI )
         ENDIF

      ELSE

         FSEN = ( PMLAI - CMLAI ) / PMLAI
         FMAT = 1.d0 - FSEN
         FGRO = 0.d0
         FNEW = 0.d0

      ENDIF

      ! Age factor
      M_AGE = FNEW*ANEW + FGRO*AGRO + FMAT*AMAT + FSEN*ASEN 

      !--------------------------------------------------------------
      ! Compute M_H = 1.d0 + 0.03d0 * (H_MONTH - H_SEASON)
      ! 
      ! where H_MONTH is the daytime sensible heat flux of the past 
      ! month and H_SEASON is the daytime sensible heat flux of the 
      ! entire growing season
      !--------------------------------------------------------------

      ! Disable M_H for now (tmf, 10/24/05)
      M_H = 1.d0     

      !-------------------------------------
      ! Compute MEA = M_LAI * M_AGE * M_H
      !-------------------------------------

      ! Prevent negative values
      MEA = MAX( M_LAI * M_AGE * M_H, 0d0 )

      END FUNCTION GET_MEA
!EOC
!------------------------------------------------------------------------------

      end module GmiEmissionMEGAN_mod
