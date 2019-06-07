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
      use fastj                    , only : Control_Fastj, GetQAA_RAAinFastj
      use fast_JX                  , only : Control_Fast_JX, GetQAA_RAAinFastJX
      use fast_JX53b               , only : Control_Fast_JX53b
      use fast_JX53b               , only : GetQAA_RAAinFastJX53b
      use fast_JX53c               , only : Control_Fast_JX53c
      use fast_JX53c               , only : GetQAA_RAAinFastJX53c
      use fastJX65_mod             , only : controlFastJX65
      use fastJX65_mod             , only : getQAA_RAAinFastJX65
      use GmiAerDustODSA_mod       , only : Aero_OptDep_SurfArea
      use GmiAerDustODSA_mod       , only : Dust_OptDep_SurfArea
      use GmiTimeControl_mod       , only : GmiSplitDateTime, GetDaysFromJanuary1
      use GmiTimeControl_mod       , only : ConvertTimeToSeconds
      use GmiArrayBundlePointer_mod, only : t_GmiArrayBundle
      use GmiFastJX_includeMod,      ONLY : t_fastJXbundle
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
     &               pr_qj_o3_o1d,  pr_qj_opt_depth, photintv, rsec_jan1,      &
     &               pctm2, mass3, pres3e, pres3c, temp3, concentration,       &
     &               solarZenithAngle, mcor, surf_alb_uv,              &
     &               fracCloudCover, tau_cloud, tau_clw, tau_cli,              &
     &               overheadO3col, qjgmi, gridBoxHeight, OptDepth,            &
     &               Eradius, Tarea, Odaer, relativeHumidity, Odmdust, Dust,   &
     &               Waersl, Daersl, humidity, num_AerDust, phot_opt,          &
     &               fastj_opt, fastj_offset_sec, do_clear_sky,                &
     &               do_AerDust_Calc, do_ozone_inFastJX, do_synoz, qj_timpyr,  &
     &               io3_num, ih2o_num, isynoz_num, chem_mask_khi, nymd, nhms, &
     &               pr_diag, loc_proc, synoz_threshold, AerDust_Effect_opt,   &
     &               num_species, num_qjs, num_qjo, ilo, ihi, julo, jhi, &
     &               i1, i2, ju1, j2, k1, k2, jNOindex, jNOamp)
!
      implicit none
!
#     include "gmi_phys_constants.h"
#     include "gmi_time_constants.h"
#     include "phot_lookup.h"
#     include "phot_monthly.h"
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
      INTEGER, INTENT(IN) :: jNOindex
      REAL   , INTENT(IN) :: jNOamp
      real*8 , intent(in) :: synoz_threshold
                             ! photolysis time step  (s)
      real*8 , intent(in) :: photintv
                             ! seconds from Jan. 1st (s)
      real*8 , intent(in) :: rsec_jan1
      real*8 , intent(in) :: fastj_offset_sec
                             ! area of grid box (m^2)
      real*8 , intent(in) :: mcor  (i1:i2, ju1:j2)
                             ! Tropopause pressure (hPa)
      REAL*8 , INTENT(IN) :: tropp(i1:i2, ju1:j2)
      real*8 , intent(in) :: solarZenithAngle(i1:i2, ju1:j2)
      real*8 , intent(in) :: fracCloudCover(i1:i2, ju1:j2)
                             ! height of each grid box (m)
      real*8 , intent(in) :: gridBoxHeight (i1:i2, ju1:j2, k1:k2)
                             ! specific humidity
      real*8 , intent(in) :: humidity (i1:i2, ju1:j2, k1:k2)
                             ! relative humidity
      REAL*8 , intent(in) :: relativeHumidity (i1:i2, ju1:j2, k1:k2)
      real*8 , intent(in) :: Dust  (i1:i2,ju1:j2,k1:k2,NSADdust)
      real*8 , intent(in) :: Waersl(i1:i2,ju1:j2,k1:k2,NSADaer)
      real*8 , intent(in) :: Daersl(i1:i2,ju1:j2,k1:k2,2      )
                             ! CTM surface pressure at t1+tdt (mb)
      real*8 , intent(in) :: pctm2 (ilo:ihi, julo:jhi)
                             ! total mass of the atmosphere within each grid box (kg)
      real*8 , intent(in) :: mass3 (i1:i2, ju1:j2, k1:k2)
                             ! atmospheric pressure at the center of each 
                             ! grid box (mb)
      real*8 , intent(in) :: pres3c(ilo:ihi, julo:jhi, k1:k2)
                             ! atmospheric pressure at the edge of each 
                             ! grid box (mb)
      real*8 , intent(in) :: pres3e(ilo:ihi, julo:jhi, k1-1:k2)
                             ! emperature (degK)
      real*8 , intent(in) :: temp3 (ilo:ihi, julo:jhi, k1:k2)
                             ! bulk surface albedo (fraction 0-1)
      real*8 , intent(in) :: surf_alb_uv(i1:i2, ju1:j2)
                             ! optical depth (dimensionless)
      real*8 , intent(in) :: tau_cloud  (i1:i2, ju1:j2, k1:k2)
                             ! optical thickness for liquid cloud
      real*8 , intent(in) :: tau_clw    (i1:i2, ju1:j2, k1:k2)
                             ! optical thickness for ice cloud
      real*8 , intent(in) :: tau_cli    (i1:i2, ju1:j2, k1:k2)
      real*8  :: overheadO3col (i1:i2, ju1:j2, k1:k2)
!
! !INPUT/OUTPUT PARAMETERS:
      real*8 , intent(inOut) :: OptDepth(i1:i2, ju1:j2, k1:k2, num_AerDust)
      real*8 , intent(inOut) :: ERADIUS (i1:i2, ju1:j2, k1:k2, NSADdust+NSADaer)
      real*8 , intent(inOut) :: TAREA   (i1:i2, ju1:j2, k1:k2, NSADdust+NSADaer)
      REAL*8 , intent(inOut) :: ODAER   (i1:i2, ju1:j2, k1:k2, NSADaer*NRH_b)
      real*8 , intent(inOut) :: ODmdust (i1:i2, ju1:j2, k1:k2, NSADdust)
                             ! photolysis rate constants (s^-1)
      type (t_GmiArrayBundle), intent(inOut) :: qjgmi(num_qjo)
                             ! species concentration, at zone centers (mixing ratio)
      type (t_GmiArrayBundle), intent(inOut) :: concentration(num_species)
      type (t_fastJXbundle)  , intent(inOut) :: JXbundle
!
! !DESCRIPTION:
!  This routine updates the photolysis rate constants (i.e., qjs).
!
! !LOCAL VARIABLES:
      integer :: idumday, idumyear
      integer :: il, ij, it, ic
      integer :: jday
      integer :: month_gmi
      real*8  :: time_sec, sza_ij
      real*8  :: overheadO3col_ij     (k1:k2)
      real*8  :: kel_ij     (k1:k2)
      real*8  :: optdepth_ij(k1:k2)
      real*8  :: ozone_ij   (k1:k2)
                 ! Column optical depth for aerosol
      real*8  :: ODAER_ij  (k1:k2,NSADaer*nrh_b)
                 ! Column optical depth for mineral dust
      real*8  :: ODMDUST_ij(k1:k2,NSADdust) 
      real*8  :: qjgmi_ij   (k1:chem_mask_khi, num_qjs)
      real*8  :: n2adj(i1:i2, ju1:j2, k1:k2)
      real*8  :: o2adj(i1:i2, ju1:j2, k1:k2)
      real*8 :: RAA_b(4, NP_b), QAA_b(4, NP_b)
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
                  if (fastj_opt == 0) then
                     call  GetQAA_RAAinFastj (RAA_b, QAA_b)
                  end if
                  if (fastj_opt == 1) then
                     call  GetQAA_RAAinFastJX (JXbundle, RAA_b, QAA_b)
                  end if
                  if (fastj_opt == 2) then
                     call  GetQAA_RAAinFastJX53b (RAA_b, QAA_b)
                  end if
                  if (fastj_opt == 3) then
                     call  GetQAA_RAAinFastJX53c (RAA_b, QAA_b)
                  end if         
                  if (fastj_opt == 4) then
                     call  GetQAA_RAAinFastJX65 (RAA_b, QAA_b)
                  end if
               end if
            end if
         end if

!      end if

!     ==================
      if (phot_opt == 2) then
!     ==================

        if (qj_timpyr == MONTHS_PER_YEAR) then
          call GmiSplitDateTime (nymd, idumyear, month_gmi, idumday)
          it = month_gmi
        else
          it = 1
        endif 

          do ic = 1, num_qjs
             qjgmi(ic)%pArray3D(i1:i2,ju1:j2,k1:k2) =  &
     &                          qjmon(i1:i2,ju1:j2,k1:k2,ic,it)
          end do

!     ==================
      elseif (phot_opt == 3) then
!     ==================

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
           if (do_AerDust_Calc) then
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

!       ------------------------------------------------------------------
!       Now loop over all latitudes and longitudes for this processor
!       because Fastj is set up to be a column calculation.
!       When doing OpenMP this would be a natural place to split the work.
!       ------------------------------------------------------------------

        do ij = ju1, j2
          do il = i1, i2

            sza_ij    = solarZenithAngle(il,ij)
            kel_ij(:) = temp3(il,ij,:)

!           set cloud OD for FastJ, FastJX, etc
            if ( do_clear_sky ) then
              optdepth_ij(:) = 0.0d0
            else
              optdepth_ij(:) = tau_cloud(il,ij,:)
            end if

            ODAER_ij   = 0.0d0
            ODMDUST_ij = 0.0d0

            if ((TRIM(chem_mecha) ==        'troposphere') .OR. &
     &          (TRIM(chem_mecha) ==         'strat_trop') .OR. &
     &          (TRIM(chem_mecha) == 'strat_trop_aerosol')) THEN
               if (do_AerDust_calc) then
                  ODAER_ij  (:,:) = ODAER  (il,ij,:,:)
                  ODMDUST_ij(:,:) = ODMDUST(il,ij,:,:)
               end if
            end if

            ozone_ij(:) = concentration(io3_num)%pArray3D(il,ij,:)

            if (fastj_opt == 0) then

!                 ==================
                  call Control_Fastj  &
!                 ==================
     &              (k1, k2, chem_mask_khi,  &
     &               num_qjs, month_gmi, jday, time_sec, fastj_offset_sec,  &
     &               sza_ij, pres3e(il,ij,k1-1:k2),  &
     &               kel_ij, optdepth_ij, surf_alb_uv(il,ij), qjgmi_ij, overheadO3col_ij, &
     &               ODAER_ij, ODMDUST_ij)

            elseif (fastj_opt == 1) then
               if (.not. do_ozone_inFastJX) then
!                 ==================
                  call Control_Fast_JX  &
!                 ==================
     &              (JXbundle, k1, k2, chem_mask_khi,  &
     &               num_qjs, month_gmi, jday, time_sec, fastj_offset_sec,  &
     &               sza_ij, pres3e(il,ij,k1:k2), pctm2(il,ij),  &
     &               kel_ij, optdepth_ij, surf_alb_uv(il,ij), qjgmi_ij, overheadO3col_ij, &
     &               ODAER_ij, ODMDUST_ij, ozone_ij)
               else
!                 ==================
                  call Control_Fast_JX  &
!                 ==================
     &              (JXbundle, k1, k2, chem_mask_khi,  &
     &               num_qjs, month_gmi, jday, time_sec, fastj_offset_sec,  &
     &               sza_ij, pres3e(il,ij,k1:k2), pctm2(il,ij),  &
     &               kel_ij, optdepth_ij, surf_alb_uv(il,ij), qjgmi_ij, overheadO3col_ij, &
     &               ODAER_ij, ODMDUST_ij)
               endif
            elseif (fastj_opt == 2) then
!                    ==================
                     call Control_Fast_JX53b  &
!                    ==================
     &                 (k1, k2, chem_mask_khi,  &
     &                  num_qjs, month_gmi, jday, time_sec,  &
     &                  sza_ij, pres3e(il,ij,k1:k2), pctm2(il,ij),  &
     &                  kel_ij, optdepth_ij, surf_alb_uv(il,ij), qjgmi_ij, overheadO3col_ij, &
     &                  ozone_ij)
            elseif (fastj_opt == 3) then
!                    ==================
                     call Control_Fast_JX53c  &
!                    ==================
     &                 (k1, k2, chem_mask_khi, num_qjs, month_gmi, jday, time_sec,  &
     &                  sza_ij, pres3e(il,ij,k1:k2), pctm2(il,ij),  &
     &                  kel_ij, optdepth_ij, surf_alb_uv(il,ij), qjgmi_ij, overheadO3col_ij, &
     &                   ODAER_ij, ODMDUST_ij, ozone_ij)
!!                    ==================
!                     call RunFastJX53c  &
!!                    ==================
!     &                 (k1, k2, jday, time_sec, month_gmi,  &
!     &                  sza_ij, pres3e(il,ij,k1:k2), pctm2(il,ij),  &
!     &                  kel_ij, optdepth_ij, surf_alb_uv(il,ij), qjgmi_ij,  &
!     &                  ozone_ij)
            elseif (fastj_opt == 4) then
               if (.not. do_ozone_inFastJX) then
                  call controlFastJX65 (k1, k2, chem_mask_khi, num_qjs, month_gmi,          &
     &                        jday, time_sec, sza_ij, do_clear_sky,                         &
     &                        tau_clw(il,ij,k1:k2), tau_cli(il,ij,k1:k2),                   &
     &                        pres3e(il,ij,k1:k2), pctm2(il,ij), kel_ij,                    &
     &                        surf_alb_uv(il,ij), qjgmi_ij, relativeHumidity(il,ij,k1:k2),  &
     &                        overheadO3col_ij, ODAER_ij, ODMDUST_ij, JXbundle%fjx_solar_cycle_param, ozone_ij)
               else
                  call controlFastJX65 (k1, k2, chem_mask_khi, num_qjs, month_gmi,          &
     &                        jday, time_sec, sza_ij, do_clear_sky,                         &
     &                        tau_clw(il,ij,k1:k2), tau_cli(il,ij,k1:k2),                   &
     &                        pres3e(il,ij,k1:k2), pctm2(il,ij), kel_ij,                    &
     &                        surf_alb_uv(il,ij), qjgmi_ij, relativeHumidity(il,ij,k1:k2),  &
     &                        overheadO3col_ij, ODAER_ij, ODMDUST_ij, JXbundle%fjx_solar_cycle_param)
               endif
            endif

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

!     ==============================================
      else if ((phot_opt == 4) .or. (phot_opt == 5)) then
!     ==============================================

!       ==============
        call Lookup_Qj  &
!       ==============
     &    (chem_mecha, do_clear_sky, phot_opt, io3_num, nymd, photintv,  &
     &     rsec_jan1, pres3c, temp3, concentration, solarZenithAngle,    &
     &     mcor, mass3,  fracCloudCover, qjgmi, &
     &     pr_diag, loc_proc, num_species, num_qjs, num_qjo, &
     &     ilo, ihi, julo, jhi, i1, i2, ju1, j2, k1, k2)

      end if

      if (TRIM(chem_mecha) == 'troposphere') then
!        ----------------------------------------------------------------------
!        Check to see if the mechanism has a photolysis reaction O3 + hv = 2OH.
!        If so, save the rate in the last entry of qjgmi and then the real rate
!        needs to be adjusted.  Updated to JPL 06-2 (Bryan Duncan 10/2006).
!        ----------------------------------------------------------------------
                   
         if (n_qj_O3_2OH > 0) then
                   
            n2adj(:,:,:) = 2.15d-11 * Exp (110.0d0 / temp3(i1:i2,ju1:j2,:)) * MXRN2
            o2adj(:,:,:) = 3.30d-11 * Exp ( 55.0d0 / temp3(i1:i2,ju1:j2,:)) * MXRO2
                   
            if (pr_qj_o3_o1d) then
               qjgmi(num_qjs+1)%pArray3D(:,:,:) = &
                     qjgmi(n_qj_O3_2OH)%pArray3D(:,:,:)
            end if

            qjgmi(n_qj_O3_2OH)%pArray3D(:,:,:) = &
                     qjgmi(n_qj_O3_2OH)%pArray3D(:,:,:) / &
     &                         (1.0d0 + ((n2adj(:,:,:) + o2adj(:,:,:)) /     &
     &                         (1.63d-10 * Exp(60.0d0/temp3(i1:i2,ju1:j2,:)) &
     &                        * concentration(ih2o_num)%pArray3D(:,:,:))))
                   
         end if
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
