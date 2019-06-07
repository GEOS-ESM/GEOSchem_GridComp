!------------------------------------------------------------------------------
! NASA GSFC - SIVO Code 610.3
!------------------------------------------------------------------------------
!BOP
! !MODULE: GmiThermalRateConstants_mod
!
      module GmiThermalRateConstants_mod
!
! !USES:
      use GmiTimeControl_mod, only : GmiSplitDateTime
      use GmiArrayBundlePointer_mod, only : t_GmiArrayBundle
      use GmiFlush_mod, only : GmiFlush
!
      implicit none
!
      private
!
! !PUBLIC MEMBER FUNCTIONS:
      public  :: calcThermalRateConstants, Accum_Qqjk

#     include "GmiParameters.h"
#     include "gmi_phys_constants.h"
!
! !AUTHOR:
! John Tannahill, LLNL, jrt@llnl.gov
!
! !REVISION HISTORY:
!   24January2005 - Jules Kouatchou
!     Modified the routine Update_Qk.
!     Added the variables "sadcol2", "radA" (surface area)
!     and "rhcol" (relative humidity column) that are arguments of the
!     routine Kcalc.
!   15September2005 - Jules Kouatchou
!     Modified the routine Update_Qk.
!     Added  "phot_opt" as argument and the condition
!     "if ((phot_opt == 3) .or. (phot_opt == 8))" before
!     "if (do_AerDust_Calc)". This was done to make sure that "sadcol2"
!     and "radA" are initialized only if fastJ or fastJX is employed.
!   10January 2011 - Jules Kouatchou
!     Used the provided cloud liquid water content instead of 
!     calculating it internally. THis is done when do_wetchem is set to true.
!     
!EOP
!------------------------------------------------------------------------------
      CONTAINS
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: calcThermalRateConstants
!
! !INTERFACE:
!

      subroutine calcThermalRateConstants (do_wetchem, chem_mecha, rootProc,    &
     &               num_time_steps, ih2o_num, imgas_num, nymd, rxnr_adjust_map,&
     &               pres3c, tropp, temp3, clwc, cmf, sadgmi, qkgmi,            &
     &               concentration, rxnr_adjust, Eradius, Tarea,                &
     &               relativeHumidity, conPBLFlag, do_AerDust_Calc, phot_opt,   &
     &               pr_diag, loc_proc, num_rxnr_adjust, rxnr_adjust_timpyr,    &
     &               ivert, num_sad, num_qks, num_molefrac, num_species,        &
     &               ilo, ihi, julo, jhi, i1, i2, ju1, j2, k1, k2)


      implicit none

#     include "gmi_time_constants.h"
#     include "setkin_par.h"
#     include "gmi_AerDust_const.h"
!
! !INPUT PARAMETERS:
      logical, intent(in) :: pr_diag
      logical, intent(in) :: do_AerDust_Calc, rootProc
           ! do wet chemistry?
      logical, intent(in) :: do_wetchem
      integer, intent(in) :: loc_proc
      integer, intent(in) :: num_species
      integer, intent(in) :: num_sad, num_qks, num_molefrac
      integer, intent(in) :: ivert
      integer, intent(in) :: num_rxnr_adjust, rxnr_adjust_timpyr
      integer, intent(in) :: i1, i2, ju1, j2, k1, k2
      integer, intent(in) :: ilo, ihi, julo, jhi
                             ! number of time steps
      integer, intent(in) :: num_time_steps
                             ! species concentration array index for water vapor
      integer, intent(in) :: ih2o_num
                             ! species concentration array index for air density
      integer, intent(in) :: imgas_num
                             ! current year/month/day (YYYYMMDD)
      integer, intent(in) :: nymd
      integer, intent(in) :: phot_opt
                             ! mapping of reaction rate adjustment number to 
                             ! reaction rate number
      integer, intent(in) :: rxnr_adjust_map(num_rxnr_adjust)
      integer, intent(in) :: conPBLFlag(i1:i2,   ju1:j2,   k1:k2)
                             ! atmospheric pressure at the center of each grid 
                             ! box (mb)
      real*8 , intent(in) :: pres3c(ilo:ihi, julo:jhi, k1:k2)
                             ! temperature (degK)
      real*8 , intent(in) :: temp3 (ilo:ihi, julo:jhi, k1:k2)
                             ! tropopause pressure (mb)
      real*8 , intent(in) :: tropp (i1:i2,   ju1:j2)
                             ! convective mass flux   (kg/m^2*s)
      real*8 , intent(in) :: cmf   (i1:i2,   ju1:j2,   k1:k2)
                             ! cloud liquid water content, grid box average (g/m^3)
      real*8 , intent(in) :: clwc  (i1:i2,   ju1:j2,   k1:k2)
      REAL*8 , intent(in) :: relativeHumidity (i1:i2, ju1:j2, k1:k2)
      REAL*8 , intent(in) :: ERADIUS (i1:i2, ju1:j2, k1:k2, NSADdust+NSADaer)
      REAL*8 , intent(in) :: TAREA   (i1:i2, ju1:j2, k1:k2, NSADdust+NSADaer)
                             ! Chemical mechanism
      character(len=*), intent(in) :: chem_mecha
                             ! array of reaction rate adjustment factors
      real*8 , intent(in) :: rxnr_adjust(i1:i2, ju1:j2, k1:k2, num_rxnr_adjust, &
     &                       rxnr_adjust_timpyr)
                             ! surface area densities (cm^2/cm^3)
      type (t_GmiArrayBundle), intent(in) :: sadgmi(num_sad)
                             ! species concentration, known at zone centers 
                             ! (mixing ratio)
      type (t_GmiArrayBundle), intent(in) :: concentration(num_species)

!
! !INPUT/OUTPUT PARAMETERS:
                             ! thermal rate constants (units vary)
      type (t_GmiArrayBundle), intent(inOut) :: qkgmi(num_qks)
!
! !DESCRIPTION:
!   This routine updates the thermal rate constants (i.e., qk's).
!
! !DEFINED PARAMETERS:
      logical, parameter :: DO_QK_STATS = .false.
!
! !LOCAL VARIABLES:
      integer :: idumyear, idumday
      integer :: il, ij, ic
      integer :: im, iq
      real*8  :: adcol   (k1:k2)
      real*8  :: prescol (k1:k2)
      real*8  :: tempcol (k1:k2)
      real*8  :: lwccol  (k1:k2)
      integer :: cPBLcol (k1:k2)
      real*8  :: constcol(num_molefrac, k1:k2)
      real*8  :: qkcol   (num_qks,      k1:k2)
      real*8  :: sadcol  (num_sad,      k1:k2)
      real*8  :: sadcol2 (NSADdust+NSADaer,   k1:k2)
      real*8  :: radA    (NSADdust+NSADaer,   k1:k2)
      real*8  :: rhcol   (k1:k2)
!
!EOP
!------------------------------------------------------------------------------
!BOC
      if (pr_diag) then
        Write (6,*) 'calcThermalRateConstants called by ', loc_proc
      end if

      lwccol (:) = 0.0d0

      qkcol(:,:) = 0.0d0

      do ij = ju1, j2
        do il = i1, i2

          prescol(:) = pres3c(il,ij,:)

          tempcol(:) = temp3 (il,ij,:)

          if (do_wetchem) then

              lwccol(:) = clwc(il,ij,:)

          end if

          !==================================================='
          ! Variables needed for gas/heterogeneous chemistry.'
          ! adcol   = air density (molec/cm3)'
          ! FRH     = relative humidity fraction (0-1)'
          ! radA    = effective radius of aerosol (cm)'
          ! sadcol2 = surface area of aerosol/volume of air (cm2/cm3)'
          !==================================================='

          adcol(:) = concentration(imgas_num)%pArray3D(il,ij,:)

          do ic = 1, num_molefrac
            constcol(ic,:) = concentration(ic)%pArray3D(il,ij,:) * adcol(:)
          end do

          do ic = 1, num_sad
             sadcol(ic,:) = sadgmi(ic)%pArray3D(il,ij,:)
          end do

          do ic = 1, NSADdust+NSADaer
             sadcol2(ic,:) = 0.0d0
             radA   (ic,:) = 0.0d0
          end do
          rhcol(:)     = 0.0d0

          if ((TRIM(chem_mecha) ==        'troposphere') .OR. &
     &        (TRIM(chem_mecha) ==         'strat_trop') .OR. &
     &        (TRIM(chem_mecha) == 'strat_trop_aerosol')) THEN
             if ((phot_opt == 3) .or. (phot_opt == 8)) then
                if (do_AerDust_Calc) then
                   do ic = 1, NSADdust+NSADaer
                      sadcol2(ic,:) = TAREA  (il,ij,:,ic)
                      radA   (ic,:) = ERADIUS(il,ij,:,ic)
                   end do
                   rhcol(:)     = relativeHumidity(il,ij,:)
                end if
             end if
          end if
	  
	  cPBLcol(:) = conPBLFlag(il,ij,:)

!         ==========
          call Kcalc  &
!         ==========
     &        (ivert, sadcol, sadcol2, prescol, tropp(il,ij),  cPBLcol, &
     &         tempcol, lwccol, adcol, constcol, qkcol, radA, rhcol)

          do iq = 1, num_qks
            qkgmi(iq)%pArray3D(il,ij,:) = Max (0.0d0, qkcol(iq,:))
          end do

        end do
      end do

      if (rxnr_adjust_timpyr == MONTHS_PER_YEAR) then
        call GmiSplitDateTime  &
     &    (nymd, idumyear, im, idumday)
      else
        im = 1
      end if

      IF(num_rxnr_adjust > 0) THEN
       do iq = 1, num_rxnr_adjust
        qkgmi(rxnr_adjust_map(iq))%pArray3D(:,:,:) =  &
     &    qkgmi(rxnr_adjust_map(iq))%pArray3D(:,:,:) * rxnr_adjust(:,:,:,iq,im)
       end do
      END IF

      if (DO_QK_STATS) then
!       ==================
        call Calc_Qk_Stats  &
!       ==================
     &    (qkgmi, num_qks, num_time_steps, &
     &     pr_diag, loc_proc, i1, i2, ju1, j2, k1, k2)
      end if

      return

      end subroutine calcThermalRateConstants
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Calc_Qk_Stats
!
! !INTERFACE:
!
      subroutine Calc_Qk_Stats  &
     &  (rate_array, dim4, num_time_steps, &
     &   pr_diag, loc_proc, i1, i2, ju1, j2, k1, k2)
!
      implicit none
!
! !INPUT PARAMETERS:
      logical, intent(in) :: pr_diag
      integer, intent(in) :: loc_proc
      integer, intent(in) :: i1, i2, ju1, j2, k1, k2
                             ! 4th dimension of rate_array
      integer, intent(in) :: dim4
                             ! number of time steps
      integer, intent(in) :: num_time_steps
                             ! thermal rate constants (i.e., qk's)
      type (t_GmiArrayBundle), intent(in) :: rate_array(dim4)
!
! !DESCRIPTION:
!   Writes out some qk statistics (min/maxes for each dim4).
!
! !LOCAL VARIABLES:
      integer :: iq
      real*8  :: minrate
      real*8  :: maxrate
!EOP
!------------------------------------------------------------------------------
!BOC
      if (pr_diag) then
        Write (6,*) 'Calc_Qk_Stats called by ', loc_proc
      end if

      Write (77,*) 'QK mins/maxes follow =>'
      Write (77,*) '  Step #:  ', num_time_steps+1
      Write (77,*) ' '

      do iq = 1, dim4
        minrate = Minval (rate_array(iq)%pArray3D(:,:,:))
        maxrate = Maxval (rate_array(iq)%pArray3D(:,:,:))

        Write (77,900) iq, minrate, maxrate
      end do

      Write (77,*) ' '

      call GmiFlush (77)

 900  format (i4, ', ', e20.10, ', ', e20.10)

      return

      end subroutine Calc_Qk_Stats
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Accum_Qqjk
!
! !INTERFACE:
!
      subroutine Accum_Qqjk (do_qqjk_reset, imgas_num, concentration, qjgmi,   &
     &                qkgmi, qqjgmi, qqkgmi, num_molefrac, num_species,        &
     &                num_qks, num_qjs, num_qjo, pr_diag, loc_proc, ilong, i1, &
     &                i2, ju1, j2, k1, k2)
!
      implicit none
!
! !INPUT PARAMETERS:
      logical, intent(in) :: pr_diag
      integer, intent(in) :: num_molefrac, num_species, num_qks, num_qjs, num_qjo
                             ! const array index for air density
      integer, intent(in) :: imgas_num
      integer, intent(in) :: loc_proc
      integer, intent(in) :: ilong
      integer, intent(in) :: i1, i2, ju1, j2, k1, k2
                             ! species concentration, known at zone centers
                             ! (mixing ratio)
      type (t_GmiArrayBundle), intent(in) :: concentration(num_species)
                             ! photolysis rate constants (s^-1)
      type (t_GmiArrayBundle), intent(in) :: qjgmi(num_qjo)
                             ! thermal    rate constants (units vary)
      type (t_GmiArrayBundle), intent(in) :: qkgmi(num_qks)
                             ! rates of photolytic processes (moles/m3/s)
      type (t_GmiArrayBundle), intent(inout) :: qqjgmi(num_qjo)
                             ! rates of thermal processes (moles/m3/s)
      type (t_GmiArrayBundle), intent(inout) :: qqkgmi(num_qks)
!
! !INPUT/OUTPUT PARAMETERS:
                                ! reset qqjk accumulators (i.e., qqjk output 
                                ! occurred on previous time step)?
      logical, intent(inOut) :: do_qqjk_reset
!
! !DESCRIPTION:
!   This routine accumulates the rates of the photolysis and thermal
!   processes.
!
! !LOCAL VARIABLES:
      integer :: ic, iq
      integer :: il, ij, ik
      real*8  :: const2 (i1:i2, imgas_num-1)
      real*8  :: qjgmi2 (i1:i2, num_qjs)
      real*8  :: qqjgmi2(i1:i2, num_qjs)
      real*8  :: qkgmi2 (i1:i2, num_qks)
      real*8  :: qqkgmi2(i1:i2, num_qks)

      real*8, PARAMETER ::   Avogadro = 6.0221415d+23

!
! !REVISION HISTORY:
!   Initial code.
!
!EOP
!------------------------------------------------------------------------------
!EOC
      if (pr_diag) then
        Write (6,*) 'Accum_Qqjk called by ', loc_proc
      end if

      do_qqjk_reset = .TRUE.

      const2 (:,:) = 0.0d0
      qjgmi2 (:,:) = 0.0d0
      qqjgmi2(:,:) = 0.0d0
      qkgmi2 (:,:) = 0.0d0
      qqkgmi2(:,:) = 0.0d0

      do ik = k1, k2
        do ij = ju1, j2

          do ic = 1, num_molefrac
            do il = i1, i2
              const2(il,ic) = concentration(ic)%pArray3D(il,ij,ik) * &
     &                        concentration(imgas_num)%pArray3D(il,ij,ik)
            end do
          end do

          if ((imgas_num - 1) > num_molefrac) then
            do ic = num_molefrac+1, imgas_num-1
              do il = i1, i2
                const2(il,ic) = concentration(ic)%pArray3D(il,ij,ik)
              end do
            end do
          end if

          do iq = 1, num_qjs
            do il = i1, i2
              qjgmi2(il,iq) = qjgmi(iq)%pArray3D(il,ij,ik)
            end do
          end do

          do iq = 1, num_qks
            do il = i1, i2
              qkgmi2(il,iq) = qkgmi(iq)%pArray3D(il,ij,ik)
            end do
          end do

!         =====================
          call Calc_Rate_Setkin  &
!         =====================
     &      (ilong, num_qjs, num_qks, imgas_num-1,  &
     &       qkgmi2, qjgmi2, const2, qqkgmi2, qqjgmi2)

!         Convert from molec/cm3/sec  to  moles/m3/sec

          do iq = 1, num_qjs
            do il = i1, i2
              qqjgmi(iq)%pArray3D(il,ij,ik) =  &
     &          qqjgmi2(il,iq) * (1.0e6 / Avogadro)
            end do
          end do

          do iq = 1, num_qks
            do il = i1, i2
              qqkgmi(iq)%pArray3D(il,ij,ik) =  &
     &          qqkgmi2(il,iq) * (1.0e6 / Avogadro)
            end do
          end do

        end do
      end do

      return 
      
      end subroutine Accum_Qqjk
!EOC
!------------------------------------------------------------------------------
      end module GmiThermalRateConstants_mod
