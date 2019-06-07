!------------------------------------------------------------------------------
! NASA GSFC - SIVO Code 610.3
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: GmiUpdateSAD_mod
!
! !INTERFACE:
!
      module GmiUpdateSAD_mod
!
! !USES:
      use GmiArrayBundlePointer_mod, only : t_GmiArrayBundle
      use GmiTimeControl_mod       , only : GmiSplitDateTime
      use GmiASCIIoperations_mod   , only : AsciiOpenWrite
      use GmiFlush_mod             , only : GmiFlush
!
      implicit none
!
#     include "GmiParameters.h"
#     include "gmi_phys_constants.h"
#     include "gmi_sad_constants.h"
!
      private
!
! !PUBLIC MEMBER FUNCTIONS:
      public :: updateSurfaceAreaDensities
!
! !AUTHOR:
! David B. Considine, NASA GSFC
!
! !REVISION HISTORY:
!  Original code written by David B. Considine
!  LLNL modifications by John Tannahill, jrt@llnl.gov
!  14Dec2007, included the routines in a module - Jules Kouatchou
!  14Dec2007, passed tropp and nameOfModel to Update_Sad2 - Jules Kouatchou
!  06Jan2011, removed nameOfModel from Update_Sad2 - Jules Kouatchou
!
!EOP
!------------------------------------------------------------------------------
      CONTAINS
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: updateSurfaceAreaDensities
!
! !INTERFACE:
!
      subroutine updateSurfaceAreaDensities (rateintv, tropp, press3c, press3e,  &
     &                 kel, concentration, ch4clim, h2oclim, hno3cond, hno3gas,  &
     &                 lbssad, sadgmi, h2oback, h2ocond, reffice, reffsts, vfall,&
     &                 dehydmin, dehyd_opt, h2oclim_opt, lbssad_opt, sad_opt,    &
     &                 ihno3cond_num, idehyd_num, ich4_num,       &
     &                 ihno3_num, ih2o_num, nymd, pr_diag,         &
     &                 loc_proc, num_species, num_sad,  &
     &                 lbssad_timpyr, h2oclim_timpyr, &
     &                 ilo, ihi, julo, jhi, i1, i2, ju1, j2, k1, k2,      &
     &                 londeg, latdeg, NoPSCZone, PSCMaxP, chem_mecha)

      implicit none

! !INPUT PARAMETERS:
      integer           , intent(in) :: ilo, ihi, julo, jhi
      integer           , intent(in) :: i1, i2, ju1, j2, k1, k2
      integer           , intent(in) :: num_species, num_sad
      integer           , intent(in) :: lbssad_timpyr, h2oclim_timpyr
      integer           , intent(in) :: ih2o_num, ihno3_num
      integer           , intent(in) :: ihno3cond_num, nymd
      integer           , intent(in) :: idehyd_num, ich4_num
      integer           , intent(in) :: loc_proc, dehyd_opt
      integer           , intent(in) :: sad_opt, h2oclim_opt, lbssad_opt
      logical           , intent(in) :: pr_diag
      real*8            , intent(in) :: rateintv
      real*8            , intent(in) :: kel    (ilo:ihi, julo:jhi, k1:k2)
      real*8            , intent(in) :: tropp      (i1:i2, ju1:j2)
      real*8            , intent(in) :: press3c(ilo:ihi, julo:jhi, k1:k2)
      real*8            , intent(in) :: press3e(ilo:ihi, julo:jhi, k1-1:k2)
      character (len=*) , intent(in) :: chem_mecha
      real*8            , intent(in) :: h2ocond(i1:i2,   ju1:j2,   k1:k2)

      integer           , intent(in) :: NoPSCZone
      integer           , intent(in) :: PSCMaxP
      real*8            , intent(in) :: londeg(i1:i2, ju1:j2), latdeg(i1:i2, ju1:j2)
!
! !OUTPUT PARAMETERS:
      real*8 , intent(out) :: dehydmin
      real*8  :: ch4clim (i1:i2,   ju1:j2,   k1:k2, h2oclim_timpyr)
      real*8  :: h2oclim (i1:i2,   ju1:j2,   k1:k2, h2oclim_timpyr)
      real*8  :: hno3gas (i1:i2,   ju1:j2,   k1:k2)
      real*8  :: lbssad  (i1:i2,   ju1:j2,   k1:k2, lbssad_timpyr)
      real*8  :: h2oback (i1:i2,   ju1:j2,   k1:k2)
      real*8  :: hno3cond (i1:i2,   ju1:j2,   k1:k2)
      real*8  :: reffice (i1:i2,   ju1:j2,   k1:k2)
      real*8  :: reffsts (i1:i2,   ju1:j2,   k1:k2)
      real*8  :: vfall   (i1:i2,   ju1:j2,   k1:k2)
!
! !INPUT/OUTPUT PARAMETERS:
      type (t_GmiArrayBundle), intent(inOut) :: sadgmi(num_sad)
      type (t_GmiArrayBundle), intent(inOut) :: concentration(num_species)

!
! !DESCRIPTION:
! Interface for the Surface Area Densities calculations.
!
! !LOCAL VARIABLES:
      integer :: idumday, idumyear, month, ic
      real*8, allocatable  :: h2ocombo(:, :, :, :)
      real*8, allocatable  :: sadcombo(:, :, :, :)
!
!EOP
!-----------------------------------------------------------------------------
!BOC
      if (pr_diag) Write (6,*) 'updateSurfaceAreaDensities called by ', loc_proc

      if (sad_opt == 1) then
        call Update_Sad1 (ihno3_num, concentration, hno3cond, hno3gas, sadgmi, &
     &          pr_diag, loc_proc, i1, i2, ju1, j2, k1, k2, num_sad, num_species)
      else if (sad_opt == 2) then
         IF (TRIM(chem_mecha) ==         'strat_trop' .OR. &
	     TRIM(chem_mecha) == 'strat_trop_aerosol') THEN
            allocate(h2ocombo(i1:i2,  ju1:j2,  k1:k2,  h2oclim_timpyr))
            allocate(sadcombo(i1:i2,  ju1:j2,  k1:k2,  num_sad       ))

            ! Note: h2ocombo is not used when dehyd_opt = 0.
            ! ---------------------------------------------

            h2ocombo(:,:,:,:) = h2oclim(:,:,:,:)

            call GmiSplitDateTime (nymd, idumyear, month, idumday)

            where (press3c(i1:i2,ju1:j2,:) > Spread (tropp(:,:), 3, k2))
               h2ocombo(:,:,:,month) = concentration(ih2o_num)%pArray3D(:,:,:)
            end where

            call Update_Sad2 (rateintv, tropp, press3c, press3e, kel,          &
     &                 concentration, ch4clim, h2ocombo, hno3cond, hno3gas,    &
     &                 lbssad, sadgmi, h2oback, h2ocond, reffice, reffsts,     &
     &                 vfall, dehydmin, dehyd_opt, h2oclim_opt, lbssad_opt,    &
     &                 pr_diag, loc_proc, londeg, latdeg, NoPSCZone, PSCMaxP, nymd, &
     &                 ihno3_num, ihno3cond_num, idehyd_num, ich4_num,         &
     &                 ih2o_num, ilo, ihi, julo, jhi, i1, i2, ju1, j2, k1,     &
     &                 k2, lbssad_timpyr, h2oclim_timpyr, num_sad, num_species)

            do ic = 1, num_sad
               sadcombo(:,:,:,ic) = sadgmi(ic)%pArray3D(:,:,:)
            end do

            call Update_Sad3 (press3c, kel, concentration, lbssad, lbssad_opt, &
     &                 sadcombo, pr_diag, loc_proc, nymd, ih2o_num, i1, i2,    &
     &                 ju1, j2, k1, k2, ilo, ihi, julo, jhi, lbssad_timpyr,    &
     &                 num_sad, num_species)

            where (press3c(i1:i2,ju1:j2,:) > Spread (tropp(:,:), 3, k2))
               sadgmi(ILBSSAD)%pArray3D(:,:,:) = sadcombo(:,:,:,ILBSSAD)
            end where

            deallocate(h2ocombo)
            deallocate(sadcombo)
         ELSE

            call Update_Sad2 (rateintv, tropp, press3c, press3e, kel,          &
     &                 concentration, ch4clim, h2oclim, hno3cond, hno3gas,     &
     &                 lbssad, sadgmi, h2oback, h2ocond, reffice, reffsts,     &
     &                 vfall, dehydmin, dehyd_opt, h2oclim_opt, lbssad_opt,    &
     &                 pr_diag, loc_proc, londeg, latdeg, NoPSCZone, PSCMaxP, nymd, &
     &                 ihno3_num, ihno3cond_num, idehyd_num, ich4_num,         &
     &                 ih2o_num, ilo, ihi, julo, jhi, i1, i2, ju1, j2, k1,     &
     &                 k2, lbssad_timpyr, h2oclim_timpyr, num_sad, num_species)

         ENDIF
      else if (sad_opt == 3) then
        allocate(sadcombo(i1:i2,  ju1:j2,  k1:k2,  num_sad       ))
        do ic = 1, num_sad
           sadcombo(:,:,:,ic) = sadgmi(ic)%pArray3D(:,:,:)
        end do

        call Update_Sad3 (press3c, kel, concentration, lbssad, lbssad_opt,     &
     &             sadcombo, pr_diag, loc_proc, nymd, ih2o_num, i1, i2, ju1, j2, &
     &             k1, k2, ilo, ihi, julo, jhi, lbssad_timpyr, num_sad,        &
     &             num_species)

        do ic = 1, num_sad
           sadgmi(ic)%pArray3D(:,:,:) = sadcombo(:,:,:,ic)
        end do
        deallocate(sadcombo)
      end if

      end subroutine updateSurfaceAreaDensities
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Update_Sad1
!
! !INTERFACE:
!
      subroutine Update_Sad1 (ihno3_num, concentration, hno3cond, hno3gas,     &
     &                 sadgmi, pr_diag, loc_proc, i1, i2, ju1, j2, k1, k2,     &
     &                 num_sad, num_species)
!
      implicit none
!
! !INPUT PARAMETERS:
      logical, intent(in) :: pr_diag
      integer, intent(in) :: loc_proc, i1, i2, ju1, j2, k1, k2
      integer, intent(in) :: num_sad, num_species
           ! index for HNO3 in the species concentration pointer
      integer, intent(in) :: ihno3_num
           ! species concentration, known at zone centers (mixing ratio)
      type (t_GmiArrayBundle), intent(in) :: concentration(num_species)
!
! !OUTPUT PARAMETERS:
           ! condensed phase hno3 array (mixing ratio)
      real*8 , intent(out) :: hno3cond(i1:i2, ju1:j2, k1:k2)
           ! gas       phase hno3 array (mixing ratio)
      real*8 , intent(out) :: hno3gas (i1:i2, ju1:j2, k1:k2)
           ! surface area densities (cm^2/cm^3)
      type (t_GmiArrayBundle), intent(out) :: sadgmi(num_sad)
!
! !DESCRIPTION:
!  Sets the stratospheric sulfuric acid aerosol surface area densities.
!
! !LOCAL VARIABLES:
      integer :: ic
!
!EOP
!------------------------------------------------------------------------------
!BOC
      if (pr_diag) then
        Write (6,*) 'Update_Sad1 called by ', loc_proc
      end if

      do ic = 1, num_sad
         sadgmi(ic)%pArray3D(:,:,:) = 0.0d0
      end do

      hno3cond(:,:,:) = 0.0d0

      hno3gas (:,:,:) = concentration(ihno3_num)%pArray3D(:,:,:)

      return

      end subroutine Update_Sad1
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Update_Sad2
!
! !INTERFACE:
!
      subroutine Update_Sad2 (sadintv, tropp, pres3c, pres3e, temp3,           &
     &                 concentration, ch4clim, h2oclim,  hno3cond, hno3gas,    &
     &                 lbssad, sadgmi, h2oback, h2ocond, reffice, reffsts,     &
     &                 vfall, dehydmin, dehyd_opt, h2oclim_opt, lbssad_opt,    &
     &                 pr_diag, loc_proc, londeg, latdeg, NoPSCZone, PSCMaxP, nymd, &
     &                 ihno3_num, ihno3cond_num, idehyd_num, ich4_num,         &
     &                 ih2o_num, ilo, ihi, julo, jhi, i1, i2, ju1, j2, k1,     &
     &                 k2, lbssad_timpyr, h2oclim_timpyr, num_sad, num_species)
!
      implicit none
!
! !INPUT PARAMETERS
      logical, intent(in) :: pr_diag
      integer, intent(in) :: loc_proc
      integer, intent(in) :: ilo, ihi, julo, jhi
      integer, intent(in) :: i1, i2, ju1, j2, k1, k2
           ! 0=disabled (AGCM mode) WARNING: Must be disabled.
      integer, intent(in) :: dehyd_opt
      INTEGER, INTENT(IN) :: h2oclim_opt, lbssad_opt
      integer, intent(in) :: lbssad_timpyr, h2oclim_timpyr, num_sad, num_species
      integer, intent(in) :: nymd
      integer, intent(in) :: ihno3_num, idehyd_num
      integer, intent(in) :: ich4_num, ih2o_num
      integer, intent(in) :: ihno3cond_num
           ! Tropopause pressure (hPa)
      REAL*8 , INTENT(IN) :: tropp(i1:i2, ju1:j2)
           ! surface area density time step (s)
      real*8 , intent(in) :: sadintv
           ! atmospheric pressure at the center of each grid box (mb)
      real*8 , intent(in) :: pres3c  (ilo:ihi, julo:jhi, k1:k2)
           ! atmospheric pressure at the edge   of each grid box (mb)
      real*8 , intent(in) :: pres3e  (ilo:ihi, julo:jhi, k1-1:k2)
           ! temperature (degK)
      real*8 , intent(in) :: temp3   (ilo:ihi, julo:jhi, k1:k2)
           ! condensed phase h2o array (mixing ratio)
      real*8 , intent(in) :: h2ocond (i1:i2,   ju1:j2,   k1:k2)

      integer, intent(in) :: NoPSCZone
      integer, intent(in) :: PSCMaxP
      real*8 , intent(in) :: londeg (i1:i2, ju1:j2), latdeg (i1:i2, ju1:j2)
!
! !OUTPUT PARAMETERS
      real*8 , intent(out):: dehydmin
           ! condensed phase hno3 array (mixing ratio)
      real*8 , intent(out) :: hno3cond(i1:i2,   ju1:j2,   k1:k2)
           ! gas       phase hno3 array (mixing ratio)
      real*8 , intent(out) :: hno3gas (i1:i2,   ju1:j2,   k1:k2)
           ! background      h2o array (mixing ratio)
      real*8 , intent(out)  :: h2oback (i1:i2,   ju1:j2,   k1:k2)

           ! array of ch4 climatology
      real*8  :: ch4clim (i1:i2,   ju1:j2,   k1:k2, h2oclim_timpyr)
           ! array of h2o climatology
      real*8  :: h2oclim (i1:i2,   ju1:j2,   k1:k2, h2oclim_timpyr)
           ! liquid binary sulfate background surface area density (cm^-1)
      real*8  :: lbssad  (i1:i2,   ju1:j2,   k1:k2, lbssad_timpyr)
           ! effective radius of ICE aerosols  (cm)
      real*8  :: reffice (i1:i2,   ju1:j2,   k1:k2)
           ! effective radius of STS aerosols  (cm)
      real*8  :: reffsts (i1:i2,   ju1:j2,   k1:k2)
           ! effective aerosol fall velocities (cm/s)
      real*8  :: vfall   (i1:i2,   ju1:j2,   k1:k2)
!
! !INPUT/OUTPUT PARAMETERS
           ! surface area densities (cm^2/cm^3)
      type (t_GmiArrayBundle), intent(inOut) :: sadgmi(num_sad)
           ! species concentration, known at zone centers (mixing ratio)
      type (t_GmiArrayBundle), intent(inout) :: concentration(num_species)
!
! !DESCRIPTION:
!  Sets the stratospheric sulfuric acid aerosol surface area densities.
!
! !DEFINED PARAMETERS:
      logical, parameter :: WR_ARGS = .false.
      integer, parameter :: WR_PROC =  1
      integer, parameter :: X1 = 36
      integer, parameter :: X2 = 36
      integer, parameter :: Y1 =  4
      integer, parameter :: Y2 =  4
      integer, parameter :: Z1 =  1
      integer, parameter :: Z2 =  1
!
! !LOCAL VARIABLES:
      logical :: is_before
      integer :: idumday
      integer :: idumyear
      integer :: ik
      integer :: month
      real*8  :: fac
      real*8  :: dehyd   (i1:i2, ju1:j2, k1:k2)
      real*8  :: denssts (i1:i2, ju1:j2, k1:k2)
      real*8  :: dz      (i1:i2, ju1:j2, k1:k2)
      real*8  :: h2so4gas(i1:i2, ju1:j2, k1:k2)
      real*8  :: reffnat (i1:i2, ju1:j2, k1:k2)
      real*8  :: sadgmi_iice(i1:i2, ju1:j2, k1:k2)
      real*8  :: sadgmi_ilbs(i1:i2, ju1:j2, k1:k2)
      real*8  :: sadgmi_inat(i1:i2, ju1:j2, k1:k2)
      real*8  :: sadgmi_ists(i1:i2, ju1:j2, k1:k2)
!
!EOP
!------------------------------------------------------------------------------
!BOC

      if (pr_diag) then
        Write (6,*) 'Update_Sad2 called by ', loc_proc
      end if

      dehydmin = 0.00D+00
      denssts (:,:,:) = 0.0d0

      h2so4gas(:,:,:) = 0.0d0
      reffnat(:,:,:) = 0.0d0

      call GmiSplitDateTime (nymd, idumyear, month, idumday)

      fac = ((GAS_CONST_J * GPKG) / (GMI_G * MWTAIR)) * CMPM

      hno3cond(:,:,:) = 0.0d0

      hno3gas(:,:,:) = concentration(ihno3_num)%pArray3D(:,:,:)
      IF(h2oclim_opt == 3) THEN
       hno3gas(:,:,:) = hno3gas(:,:,:)+concentration(ihno3cond_num)%pArray3D(:,:,:)
      END IF
      
      IF(dehyd_opt == 0) THEN
        dehyd(:,:,:) = 0.0d0
      ELSE
        PRINT *,"Update_Sad2: dehyd_opt must be zero."
	STOP 
      END IF

      IF(h2oclim_opt == 3) THEN
       h2oback(:,:,:) =  concentration(ih2o_num)%pArray3D(:,:,:) +  &
     &                   h2ocond(:,:,:)
      ELSE
       h2oback(:,:,:) = h2oclim(:,:,:,month) + 2.0d0 *  &
     &  (concentration(ich4_num)%pArray3D(:,:,:) - ch4clim(:,:,:,month))
      END IF

      IF(lbssad_opt == 4) THEN
        sadgmi(ILBSSAD)%pArray3D(:,:,:) = lbssad(:,:,:,1)
      ELSE
        sadgmi(ILBSSAD)%pArray3D(:,:,:) = lbssad(:,:,:,month)
      END IF

      do ik = k1, k2

        dz(i1:i2,ju1:j2,ik) =  &
     &    fac * temp3(i1:i2,ju1:j2,ik) *  &
     &    Log (pres3e(i1:i2,ju1:j2,ik-1) / pres3e(i1:i2,ju1:j2,ik))

      end do

      sadgmi_iice(:,:,:) = sadgmi(IICESAD)%pArray3D(:,:,:)
      sadgmi_ilbs(:,:,:) = sadgmi(ILBSSAD)%pArray3D(:,:,:)
      sadgmi_inat(:,:,:) = sadgmi(INATSAD)%pArray3D(:,:,:)
      sadgmi_ists(:,:,:) = sadgmi(ISTSSAD)%pArray3D(:,:,:)

!     =============
      call Condense  &
!     =============
     &  (sadintv, tropp, pres3c, temp3, dehyd, dehyd_opt, dz,  &
     &   h2oback, hno3cond, hno3gas, sadgmi_ilbs,  &
     &   denssts, h2ocond, h2so4gas, sadgmi_iice,  &
     &   sadgmi_inat, sadgmi_ists, reffice, reffnat, reffsts, vfall, &
     &   pr_diag, loc_proc, londeg, latdeg, NoPSCZone, PSCMaxP, &
     &   ilo, ihi, julo, jhi, i1, i2, ju1, j2, k1, k2)

      sadgmi(IICESAD)%pArray3D(:,:,:) = sadgmi_iice(:,:,:)
      sadgmi(ILBSSAD)%pArray3D(:,:,:) = sadgmi_ilbs(:,:,:)
      sadgmi(INATSAD)%pArray3D(:,:,:) = sadgmi_inat(:,:,:)
      sadgmi(ISTSSAD)%pArray3D(:,:,:) = sadgmi_ists(:,:,:)

      call GmiSplitDateTime (nymd, idumyear, month, idumday)

!     For now, just zero out sootsad.
      sadgmi(ISOOTSAD)%pArray3D(:,:,:) = 0.0d0

      IF(dehyd_opt /= 0) THEN
        PRINT *,"Update_Sad2: dehyd_opt must be zero."
	STOP 
      END IF

      concentration(ihno3_num)%pArray3D(:,:,:)  = hno3gas(:,:,:)
      IF(h2oclim_opt == 3) THEN
        concentration(ihno3cond_num)%pArray3D(:,:,:)  = hno3cond(:,:,:)
      END IF

! In GEOS-5 apply only this bounds check on water
! -----------------------------------------------

      concentration(ih2o_num)%pArray3D(:,:,:) =  &
     &  Max (concentration(ih2o_num)%pArray3D(:,:,:), 0.01d-06)

      return

      end subroutine Update_Sad2
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Update_Sad3
!
! !INTERFACE:
!
      subroutine Update_Sad3 (pres3c, temp3, concentration, lbssad, lbssad_opt,&
     &                 loc_sadgmi, pr_diag, loc_proc, nymd, ih2o_num, i1, i2, ju1, &
     &                 j2, k1, k2, ilo, ihi, julo, jhi, lbssad_timpyr, num_sad,&
     &                 num_species)
!
      implicit none
!
! !INPUT PARAMETERS:
      logical, intent(in) :: pr_diag
      integer, intent(in) :: loc_proc
      integer, intent(in) :: nymd
      integer, intent(in) :: ih2o_num
      integer, intent(in) :: i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi
      integer, intent(in) :: lbssad_timpyr, lbssad_opt, num_sad, num_species
      real*8 , intent(in) :: pres3c(ilo:ihi, julo:jhi, k1:k2)
      real*8 , intent(in) :: temp3 (ilo:ihi, julo:jhi, k1:k2)
      real*8 , intent(in) :: lbssad(i1:i2, ju1:j2, k1:k2, lbssad_timpyr)
           ! species concentration, known at zone centers (mixing ratio)
      type (t_GmiArrayBundle), intent(in) :: concentration(num_species)
!
! !INPUT/OUTPUT PARAMETERS:
           ! surface area densities (cm^2/cm^3)
      real*8, intent(inOut) :: loc_sadgmi(i1:i2, ju1:j2, k1:k2,num_sad)
!
! !DESCRIPTION:
!  Sets the tropospheric sulfuric aerosol surface area densities.
!
! !LOCAL VARIABLES:
      integer :: idumday, idumyear
      integer :: month
      real*8  :: rh(i1:i2, ju1:j2, k1:k2)
!EOP
!------------------------------------------------------------------------------
!BOC
      if (pr_diag) then
        Write (6,*) 'Update_Sad3 called by ', loc_proc
      end if

!     -----------------------------------------------------------
!     Calculate relative humidity from Seinfeld (1986) p. 181.
!     The first rh is the temperature dependent parameter a.
!     The second rh is the saturation vapor pressure of water.
!     The third rh is the actual relative humidity as a fraction.
!     Then make sure rh is between 0 and 1.
!     -----------------------------------------------------------

      rh(:,:,:) = 1.0d0 - (373.15d0 / temp3(i1:i2,ju1:j2,:))

      rh(:,:,:) =  &
     &  1013.25d0 * Exp (13.3185d0 * rh(:,:,:)    -  &
     &                    1.9760d0 * rh(:,:,:)**2 -  &
     &                    0.6445d0 * rh(:,:,:)**3 -  &
     &                    0.1299d0 * rh(:,:,:)**4)

      rh(:,:,:) =  &
     &  concentration(ih2o_num)%pArray3D(:,:,:) *  &
     &  pres3c(i1:i2,ju1:j2,:) / rh(:,:,:)

      rh(:,:,:) = Max (Min (rh(:,:,:), 1.0d0), 0.0d0)

!     -----------------------------------------------------------------
!     Choose month only when all 12 months of lbssad are available.
!     -----------------------------------------------------------------

      IF(lbssad_opt == 4) THEN
       month = 1
      ELSE
       CALL GmiSplitDateTime (nymd, idumyear, month, idumday)
      END IF

!     -----------------------------------------------------------------
!     Now calculate the change in surface area density due to humidity.
!     This is a fit to data from Chuang in the form suggested by Grant.
!     -----------------------------------------------------------------

      loc_sadgmi(:,:,:, 1) =  &
     &  lbssad(:,:,:,month) *  &
     &  (1.0d0 +  &
     &   1.25824d0 * (Exp (2.71811d0 * rh(:,:,:)**11.3214) - 1.0d0) +  &
     &   24.94300d0 * (1.0d0 - Exp (-0.0329662 * rh(:,:,:))))**2

      return

      end subroutine Update_Sad3
!EOC
!------------------------------------------------------------------------------
      end module GmiUpdateSAD_mod
