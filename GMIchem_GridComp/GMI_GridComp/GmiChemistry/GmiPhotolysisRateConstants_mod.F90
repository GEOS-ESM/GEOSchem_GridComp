!------------------------------------------------------------------------------
! NASA GSFC - SIVO Code 610.3
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: GmiPhotRateConst_mod
!
      module GmiPhotRateConst_mod
!
! !USES:
      use fastJX65_mod             , only : controlFastJX65
      use fastJX65_mod             , only : getQAA_RAAinFastJX65
      use CloudJ_mod               , only : controlFastJX74
      use GmiAerDustODSA_mod       , only : Aero_OptDep_SurfArea
      use GmiAerDustODSA_mod       , only : Dust_OptDep_SurfArea
      use GmiTimeControl_mod       , only : GmiSplitDateTime, GetDaysFromJanuary1
      use GmiTimeControl_mod       , only : ConvertTimeToSeconds
      use GmiArrayBundlePointer_mod, only : t_GmiArrayBundle
      use GmiFastJX_includeMod     , only : t_fastJXbundle
      use MAPL_BaseMod             , only : MAPL_UNDEF
!
      implicit none
!
      private
!
! !PUBLIC MEMBER FUNCTIONS:
      public  :: calcPhotolysisRateConstants
!
! !AUTHOR:
!   John Tannahill, LLNL, jrt@llnl.gov
!   Jules Kouatchou, NASA/GSFC, Jules.Kouatchou-1@nasa.gov
!
! !REVISION HISTORY:
!   - 24January2005 * Jules Kouatchou
!     Modified the routine Update_Qj to compute aerosol/dust
!     surface area and optical depth. The routines used for calculations
!     are Aero_OptDep_SurfArea and Dust_OptDep_SurfArea.
!     The changes are valid only for the troposphere chemical mechanism.
!     The argument "gridBoxHeight" was also added in the subroutine
!     Update_Qj for the above calculations.
!     Pass on the values of optical depht (ODAER_ij, ODMDUST_ij)
!     aerosols/dust to the fastj control routine.
!   - 9March2005 * Jules Kouatchou
!     The call of Control_Fastj can now be done for any chemical mechanism.
!     Two arguments (ODAER_ij, ODMDUST_ij) were added to Control_FastJX and
!     the variable "ozone_ij" is now an optional argument of Control_FastJX.
!     For this reason Control_FastJX is declared here as an interface.
!     "ozone_ij" is not used for the tropospheric chemical mechanism because
!     the code uses climatology provide by FastJX instead of the model
!     climatology. For the combo mechanism, there is a namelist variable
!     "do_ozone_inFastJX" for selecting a particular climatology.
!   - 14March2005 * Jules Kouatchou
!     Computations of optical depth and surface area for different
!     aerosols/dust are done only if "do_AerDust_Calc" is set to TRUE.
!   - 26May2016 * Luke Oman
!     Added option for FastJX 6.5
!     And use Pressures on Edges instead of on Centers (fixed long-standing bug)
!
!EOP
!------------------------------------------------------------------------------
      CONTAINS
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: calcPhotolysisRateConstants
!
! !INTERFACE:
!

      subroutine calcPhotolysisRateConstants( JXbundle, chem_mecha, tropp,     &
     &               pr_qj_o3_o1d,  pr_qj_opt_depth,                           &
     &               pctm2, mass3, pres3e, pres3c, temp3, concentration,       &
     &               solarZenithAngle, mcor, surf_alb_uv,                      &
     &               fracCloudCover, tau_cloud, tau_clw, tau_cli,              &
     &               totalCloudFraction, qi, ql, ri, rl,                       &
!    &               cnv_frc, frland,                                          &
     &               overheadO3col, qjgmi, gridBoxHeight, OptDepth,            &
     &               Eradius, Tarea, Odaer, relativeHumidity, Odmdust, Dust,   &
     &               Waersl, Daersl, humidity, num_AerDust, phot_opt,          &
     &               fastj_opt, fastj_offset_sec, do_clear_sky,                &
     &               do_AerDust_Calc, do_ozone_inFastJX, do_synoz, qj_timpyr,  &
     &               io3_num, ih2o_num, isynoz_num, chem_mask_khi, nymd, nhms, &
     &               pr_diag, loc_proc, synoz_threshold, AerDust_Effect_opt,   &
     &               num_species, num_qjs, num_qjo, ilo, ihi, julo, jhi, &
     &               i1, i2, ju1, j2, k1, k2, jNOindex, jNOamp, cldflag)
!
      implicit none
!
#     include "gmi_phys_constants.h"
#     include "gmi_time_constants.h"
#     include "setkin_par.h"
#     include "gmi_AerDust_const.h"
!
! !INPUT PARAMETERS:
      CHARACTER(LEN=*), INTENT(IN) :: chem_mecha
      logical, intent(in) :: pr_diag
      logical, intent(in) :: do_AerDust_calc, do_clear_sky
      logical, intent(in) :: do_ozone_inFastJX, do_synoz
                             ! should special reaction O3->O1d be saved?
      logical, intent(in) :: pr_qj_o3_o1d
                             ! should optical depth be saved at the end of qj file?
      logical, intent(in) :: pr_qj_opt_depth
      integer, intent(in) :: loc_proc
      integer, intent(in) :: num_species, num_qjs, num_qjo
      integer, intent(in) :: ilo, ihi, julo, jhi
      integer, intent(in) :: i1, i2, ju1, j2, k1, k2
      integer, intent(in) :: isynoz_num, io3_num, ih2o_num
      integer, intent(in) :: nymd, nhms
      integer, intent(in) :: chem_mask_khi, qj_timpyr
      integer, intent(in) :: AerDust_Effect_opt
      integer, intent(in) :: phot_opt, fastj_opt
      integer, intent(in) :: num_AerDust
      integer, intent(in) :: cldflag
      INTEGER, INTENT(IN) :: jNOindex
      REAL   , INTENT(IN) :: jNOamp
      real*8 , intent(in) :: synoz_threshold
      real*8 , intent(in) :: fastj_offset_sec
      real*8 , intent(in) :: mcor  (i1:i2, ju1:j2)                   ! area of grid box (m^2)
      REAL*8 , INTENT(IN) :: tropp(i1:i2, ju1:j2)                    ! Tropopause pressure (hPa)
      real*8 , intent(in) :: solarZenithAngle(i1:i2, ju1:j2)
      real*8 , intent(in) :: fracCloudCover(i1:i2, ju1:j2)
      real*8 , intent(in) :: gridBoxHeight (i1:i2, ju1:j2, k1:k2)    ! height of each grid box (m)
      real*8 , intent(in) :: humidity (i1:i2, ju1:j2, k1:k2)         ! specific humidity
      REAL*8 , intent(in) :: relativeHumidity (i1:i2, ju1:j2, k1:k2) ! relative humidity
      real*8 , intent(in) :: Dust  (i1:i2,ju1:j2,k1:k2,NSADdust)
      real*8 , intent(in) :: Waersl(i1:i2,ju1:j2,k1:k2,NSADaer)
      real*8 , intent(in) :: Daersl(i1:i2,ju1:j2,k1:k2,2      )
      real*8 , intent(in) :: pctm2 (ilo:ihi, julo:jhi)          ! CTM surface pressure at t1+tdt (mb)
      real*8 , intent(in) :: mass3 (i1:i2, ju1:j2, k1:k2)       ! total mass of the atmosphere within each grid box (kg)
      real*8 , intent(in) :: pres3c(ilo:ihi, julo:jhi, k1:k2)   ! atmospheric pressure at the center of each grid box (mb)
      real*8 , intent(in) :: pres3e(ilo:ihi, julo:jhi, k1-1:k2) ! atmospheric pressure at the edge of each grid box (mb)
      real*8 , intent(in) :: temp3 (ilo:ihi, julo:jhi, k1:k2)   ! temperature (degK)
      real*8 , intent(in) :: surf_alb_uv        (i1:i2, ju1:j2)        ! bulk surface albedo (fraction 0-1)
      real*8 , intent(in) :: totalCloudFraction (i1:i2, ju1:j2, k1:k2) ! total Cloud Fraction
      real*8 , intent(in) :: tau_cloud          (i1:i2, ju1:j2, k1:k2) ! optical depth (dimensionless)
      real*8 , intent(in) :: tau_clw            (i1:i2, ju1:j2, k1:k2) ! optical thickness for liquid cloud
      real*8 , intent(in) :: tau_cli            (i1:i2, ju1:j2, k1:k2) ! optical thickness for ice cloud
!     real*8 , intent(in) :: cnv_frc            (i1:i2, ju1:j2)
!     real*8 , intent(in) :: frland             (i1:i2, ju1:j2)
      real*8 , intent(in) :: qi                 (i1:i2, ju1:j2, k1:k2) ! in-cloud ice content
      real*8 , intent(in) :: ql                 (i1:i2, ju1:j2, k1:k2) ! in-cloud liquid content
      real*8 , intent(in) :: ri                 (i1:i2, ju1:j2, k1:k2) ! effective radius for ice
      real*8 , intent(in) :: rl                 (i1:i2, ju1:j2, k1:k2) ! effective radius for liquid
! !INPUT/OUTPUT PARAMETERS:
      real*8 , intent(inOut) :: OptDepth(i1:i2, ju1:j2, k1:k2, num_AerDust)
      real*8 , intent(inOut) :: ERADIUS (i1:i2, ju1:j2, k1:k2, NSADdust+NSADaer)
      real*8 , intent(inOut) :: TAREA   (i1:i2, ju1:j2, k1:k2, NSADdust+NSADaer)
      REAL*8 , intent(inOut) :: ODAER   (i1:i2, ju1:j2, k1:k2, NSADaer*NRH_b)
      real*8 , intent(inOut) :: ODmdust (i1:i2, ju1:j2, k1:k2, NSADdust)
! !OUTPUT PARAMTERS
      real*8 , intent(out) :: overheadO3col (i1:i2, ju1:j2, k1:k2)

      type (t_GmiArrayBundle), intent(inOut) :: qjgmi(num_qjo) ! photolysis rate constants (s^-1)
      type (t_GmiArrayBundle), intent(inOut) :: concentration(num_species) ! species concentration, at zone centers (mixing ratio)
      type (t_fastJXbundle)  , intent(inOut) :: JXbundle
!
! !DESCRIPTION:
!  This routine updates the photolysis rate constants (i.e., qjs).
!
! !LOCAL VARIABLES:
      integer :: idumday, idumyear
      integer :: il, ij, it, ic, n
      integer :: jday, ich4_num
      integer :: month_gmi
      integer, parameter :: four = 4
      real*8  :: time_sec, sza_ij, lat_ij
      real*8  :: overheadO3col_ij (k1:k2)
      real*8  ::           kel_ij (k1:k2)
      real*8  ::         cldOD_ij (k1:k2)
      real*8  :: gridBoxHeight_ij (k1:k2)
      real*8  ::         ozone_ij (k1:k2)
      real*8  ::           ch4_ij (k1:k2)
      real*8  ::           h2o_ij (k1:k2)
                 ! Column optical depth for aerosol
      real*8  :: ODcAER_ij  (k1:k2,2)
      real*8  :: ODAER_ij  (k1:k2,NSADaer*NRH_b)
      real*8  :: HYGRO_ij  (k1:k2,NSADaer)
                 ! Column optical depth for mineral dust
      real*8  :: ODMDUST_ij(k1:k2,NSADdust) 
      real*8  :: qjgmi_ij   (k1:chem_mask_khi, num_qjs)
      real*8  :: n2adj(i1:i2, ju1:j2, k1:k2)
      real*8  :: o2adj(i1:i2, ju1:j2, k1:k2)
      real*8  :: RAA_b(4, NP_b), QAA_b(4, NP_b)
      real*8  :: ERADIUS_ij (k1:k2, NSADdust+NSADaer)
      real*8  :: TAREA_ij   (k1:k2, NSADdust+NSADaer)
!
!EOP
!------------------------------------------------------------------------------
!BOC
      if (pr_diag) then
        Write (6,*) 'Update_Qj called by ', loc_proc
      end if

!      if (first) then
!         first = .false.
         if (phot_opt == 3) then
            if ((TRIM(chem_mecha) ==        'troposphere') .OR. &
     &          (TRIM(chem_mecha) ==         'strat_trop') .OR. &
     &          (TRIM(chem_mecha) == 'strat_trop_aerosol')) THEN
               if (do_AerDust_Calc) then
                  if (fastj_opt == 4) then
                     call  GetQAA_RAAinFastJX65 (RAA_b, QAA_b, four, NP_b)
                  end if
               end if
            end if
         end if

!      end if
!
!     ==================
      if (phot_opt == 3) then
!     ==================
!
!       --------------------------------------------------------------
!       First form some non-grid dependent arguments needed for Fastj.
!       --------------------------------------------------------------

        call GmiSplitDateTime (nymd, idumyear, month_gmi, idumday)

        call GetDaysFromJanuary1 (jday, nymd)

        time_sec = ConvertTimeToSeconds (nhms)

        ! For aerosol/dust optical depth/surface area calculations
        if ((TRIM(chem_mecha) ==        'troposphere') .OR. &
     &      (TRIM(chem_mecha) ==         'strat_trop') .OR. &
     &      (TRIM(chem_mecha) == 'strat_trop_aerosol')) THEN
           if (do_AerDust_Calc .and. fastj_opt /= 5) then
              call Aero_OptDep_SurfArea(gridBoxHeight, concentration, tropp,   &
     &                 pres3c, OptDepth, Eradius, Tarea, Odaer,                &
     &                 relativeHumidity, Daersl, Waersl, RAA_b, QAA_b,         &
     &                 do_synoz, isynoz_num, synoz_threshold,                  &
     &                 AerDust_Effect_opt, i1, i2, ju1, j2, k1, k2, ilo, ihi,  &
     &                 julo, jhi, num_species, num_AerDust)

              call Dust_OptDep_SurfArea(gridBoxHeight, concentration, tropp,   &
     &                 pres3c, OptDepth, Eradius, Tarea, Odmdust, Dust, RAA_b, &
     &                 QAA_b, do_synoz, isynoz_num, synoz_threshold,           &
     &                 AerDust_Effect_opt, i1, i2, ju1, j2, k1, k2,            &
     &                 num_species, num_AerDust)
           end if
        end if

        if (pr_qj_opt_depth) qjgmi(num_qjo)%pArray3D(:,:,:) = tau_cloud(:,:,:)
!
!... get methane index for CloudJ longwave calc, NOT TESTED (or needed?) yet
        if (fastj_opt .eq. 5) ich4_num = ICH4
!
!       ------------------------------------------------------------------
!       Now loop over all latitudes and longitudes for this processor
!       because Fastj is set up to be a column calculation.
!       When doing OpenMP this would be a natural place to split the work.
!       ------------------------------------------------------------------

        do ij = ju1, j2
          do il = i1, i2

            sza_ij    = solarZenithAngle(il,ij)
            kel_ij(:) = temp3(il,ij,:)
            gridBoxHeight_ij(:) = gridBoxHeight(il,ij,:)
!... NEED TO SET PROPERLY
            lat_ij    = 0.
!
            if ((TRIM(chem_mecha) ==        'troposphere') .OR. &
     &          (TRIM(chem_mecha) ==         'strat_trop') .OR. &
     &          (TRIM(chem_mecha) == 'strat_trop_aerosol')) THEN
!
!... for CloudJ aerosol OD calc input
               if (do_AerDust_calc) then
                  if (fastj_opt .eq. 5) then 
!... preprocess for CloudJ AOD calc
!... hydrophobic BC aerosols
                    ODcAER_ij(:,1) = daersl(il,ij,:,1)
!... hydrophobic OC aerosols 
                    ODcAER_ij(:,2) = daersl(il,ij,:,2)
!... dust aerosols
                    do N=1,NSADdust
                      ODMDUST_ij(:,N) = DUST(il,ij,:,N)
                    enddo
!... hydrophilic aerosols
                    do N=1,NSADaer
                      ODAER_ij(:,N) = waersl(il,ij,:,N)
                    enddo
!... other FastJ calls need AOD calc'd from Aero_OptDep_SurfArea and Dust_OptDep_SurfArea
                  else
                    ODAER_ij  (:,:) = ODAER  (il,ij,:,:)
                    ODMDUST_ij(:,:) = ODMDUST(il,ij,:,:)
                  endif
               else
                 ODAER_ij  (:,:) = 0.0d0
                 ODMDUST_ij(:,:) = 0.0d0
               endif
!
            endif
!
            ozone_ij(:) = concentration(io3_num)%pArray3D(il,ij,:)
!... needed for longwave phot calc in CloudJ (SolarJ)
            if (fastj_opt .eq. 5) then
              ch4_ij(:) = concentration(ich4_num)%pArray3D(il,ij,:)
              h2o_ij(:) = concentration(ih2o_num)%pArray3D(il,ij,:)
            endif
!
!... set up cloud OD for fastj_opt = 0 to 4
            if(do_clear_sky) then 
              cldOD_ij(:) = 0.0d0
            else
              cldOD_ij(:) = tau_cloud(il,ij,:)
            endif
!
            if (fastj_opt == 4) then
               if (.not. do_ozone_inFastJX) then
                  call controlFastJX65 (k1, k2, chem_mask_khi, num_qjs, month_gmi,          &
     &                        jday, time_sec, sza_ij, do_clear_sky,                         &
     &                        tau_clw(il,ij,k1:k2), tau_cli(il,ij,k1:k2),                   &
     &                        pres3e(il,ij,k1:k2), pctm2(il,ij), kel_ij,                    &
     &                        surf_alb_uv(il,ij), qjgmi_ij, relativeHumidity(il,ij,k1:k2),  &
     &                        overheadO3col_ij, ODAER_ij, ODMDUST_ij, cldOD_ij,             &
     &                        JXbundle%fjx_solar_cycle_param, ozone_ij)
               else
                  call controlFastJX65 (k1, k2, chem_mask_khi, num_qjs, month_gmi,          &
     &                        jday, time_sec, sza_ij, do_clear_sky,                         &
     &                        tau_clw(il,ij,k1:k2), tau_cli(il,ij,k1:k2),                   &
     &                        pres3e(il,ij,k1:k2), pctm2(il,ij), kel_ij,                    &
     &                        surf_alb_uv(il,ij), qjgmi_ij, relativeHumidity(il,ij,k1:k2),  &
     &                        overheadO3col_ij, ODAER_ij, ODMDUST_ij, cldOD_ij,             &
     &                        JXbundle%fjx_solar_cycle_param)
               endif
!
!... use CloudJ (aka FastJX7.4)
            elseif (fastj_opt == 5) then
               cldOD_ij(:)     = MAPL_UNDEF
               eradius_ij(:,:) = 0.0d0
               tArea_ij(:,:)   = 0.0d0
!
               if (do_ozone_inFastJX) then
                  call controlFastJX74 (k1, k2, chem_mask_khi, lat_ij, num_qjs, month_gmi,          &
     &                        jday, time_sec, do_clear_sky, cldflag, gridBoxHeight_ij(k1:k2),       &
     &                        sza_ij, totalCloudFraction(il,ij,k1:k2),                              &
     &                        qi(il,ij,k1:k2), ql(il,ij,k1:k2),                                     &
     &                        ri(il,ij,k1:k2), rl(il,ij,k1:k2),                                     &
!    &                        cnv_frc(il,ij), frland(il,ij),                                        &
     &                        pres3e(il,ij,k1-1:k2), pctm2(il,ij), kel_ij,                          &
     &                        surf_alb_uv(il,ij), qjgmi_ij, relativeHumidity(il,ij,k1:k2),          &
     &                        overheadO3col_ij, ODAER_ij, ODMDUST_ij, ODcAER_ij, HYGRO_ij,          &
     &                        do_AerDust_Calc, AerDust_Effect_opt, cldOD_ij, eradius_ij, tArea_ij,  &
     &                        JXbundle%fjx_solar_cycle_param, CH4_ij, H2O_ij)
               else
                  call controlFastJX74 (k1, k2, chem_mask_khi, lat_ij, num_qjs, month_gmi,          &
     &                        jday, time_sec, do_clear_sky, cldflag, gridBoxHeight_ij(k1:k2),       &
     &                        sza_ij, totalCloudFraction(il,ij,k1:k2),                              &
     &                        qi(il,ij,k1:k2), ql(il,ij,k1:k2),                                     &
     &                        ri(il,ij,k1:k2), rl(il,ij,k1:k2),                                     &
!    &                        cnv_frc(il,ij), frland(il,ij),                                        &
     &                        pres3e(il,ij,k1-1:k2), pctm2(il,ij), kel_ij,                          &
     &                        surf_alb_uv(il,ij), qjgmi_ij, relativeHumidity(il,ij,k1:k2),          &
     &                        overheadO3col_ij, ODAER_ij, ODMDUST_ij, ODcAER_ij, HYGRO_ij,          &
     &                        do_AerDust_Calc, AerDust_Effect_opt, cldOD_ij, eradius_ij, tArea_ij,  &
     &                        JXbundle%fjx_solar_cycle_param, CH4_ij, H2O_ij, ozone_ij)
               endif
!
               eradius(il,ij,:,:) = eradius_ij(:,:)
               tArea(il,ij,:,:)   = tArea_ij(:,:)
!... Dust surface areas - optDepth(5)
               do N = 1, NSADdust
                 optDepth(il,ij,:,5) = optDepth(il,ij,:,5) + tArea(il,ij,:,N)
               enddo
!... wAersl and dAersl aerosol diagnostics
               do N = 1, NSADaer
!... hydrophilic hygroscopic growth diagnostic  - optDepth(7,10,13,16,19)
                 optDepth(il,ij,:,4+3*N) = optDepth(il,ij,:,4+3*N) + HYGRO_ij(:,N)
               enddo
!... hydrophilic surface areas diagnostic - optDepth(8,11,14,17,20)
               do N = 1, NSADaer
		 optDepth(il,ij,:,5+3*N) = optDepth(il,ij,:,5+3*N) + tArea(il,ij,:,NSADdust+N)
               enddo
!
!... Dust optical depths - optDepth(4)
               do N = 1, NSADdust
                 optDepth(il,ij,:,4) = optDepth(il,ij,:,4) + ODMDUST_ij(:,N)
               enddo
!
!... hydrophilic aerosol optical depths diagnostic - optDepth(6,9,12,15,18)
               do N = 1, NSADaer
                 optDepth(il,ij,:,3+3*N) = optDepth(il,ij,:,3+3*N) + ODAER_ij(:,N)
               enddo
!
!...  add in hydrophobic carbon (BC/OC) optical depths - optDepth(9,12)
               do N = 2,3
                 optDepth(il,ij,:,3+3*N) = optDepth(il,ij,:,3+3*N) + ODcAER_ij(:,N-1)
               enddo
!... convert cloud optical depth from in-cloud to gridbox average
               where(cldOD_ij(:).le.5.0d5) cldOD_ij(:) = cldOD_ij(:)*totalCloudFraction(il,ij,:)
!... end CloudJ stuff
            endif
!
!... FastJX used cloud OD diagnostic - optDepth(1)
            optDepth(il,ij,:,1) = cldOD_ij(:)
!... FastJX used cloud fraction diagnostic - optDepth(2)
            optDepth(il,ij,:,2) = totalCloudFraction(il,ij,:)
!
            do ic = 1, num_qjs
               qjgmi(ic)%pArray3D(il,ij,k1:chem_mask_khi) = qjgmi_ij(k1:chem_mask_khi,ic)
            end do

!           ----------------------------------------------------------------------------------------------------------
!                                      Michael Prather recommends a run-time adjustment to j(NO)
!           ----------------------------------------------------------------------------------------------------------
            qjgmi(jNOindex)%pArray3D(il,ij,k1:chem_mask_khi) = qjgmi(jNOindex)%pArray3D(il,ij,k1:chem_mask_khi)*jNOamp
!           ----------------------------------------------------------------------------------------------------------
 
            overheadO3col(il,ij,:) = overheadO3col_ij(:)

          end do
        end do

      end if

!!     ----------------------------------------------------------------------
!!     Check to see if the mechanism has a photolysis reaction O3 + hv = 2OH.
!!     If so, save the rate in the last entry of qjgmi and then the real rate
!!     needs to be adjusted.
!!     ----------------------------------------------------------------------
!
!!      if (n_qj_O3_2OH > 0) then
!
!        n2adj(:,:,:) =  &
!!    &   1.8d-11 * Exp (110.0d0 / temp3(i1:i2,ju1:j2,:)) * MXRN2
!!    Bryan's Experiment
!     &   2.1d-11 * Exp (115.0d0 / temp3(i1:i2,ju1:j2,:)) * MXRN2
!        o2adj(:,:,:) =  &
!     &   3.2d-11 * Exp ( 70.0d0 / temp3(i1:i2,ju1:j2,:)) * MXRO2
!
!        if (pr_qj_o3_o1d) then
!          qjgmi(:,:,:,num_qjs+1) = qjgmi(:,:,:,n_qj_O3_2OH)
!        end if
!
!        qjgmi(:,:,:,n_qj_O3_2OH) =  &
!     &    qjgmi(:,:,:,n_qj_O3_2OH) /  &
!     &    (1.0d0 +  &
!     &     ((n2adj(:,:,:) + o2adj(:,:,:)) /  &
!     &      (1.62d-10 * Exp(65.0d0/temp3(i1:i2,ju1:j2,:))  &
!     &       * concentration(ih2o_num)%pArray3D(:,:,:))))
!
!!    &      (2.2d-10 * concentration(ih2o_num)%pArray3D(:,:,:))))
!
!       end if


      return

      end subroutine calcPhotolysisRateConstants
!
!EOC
!------------------------------------------------------------------------------
      end module GmiPhotRateConst_mod
