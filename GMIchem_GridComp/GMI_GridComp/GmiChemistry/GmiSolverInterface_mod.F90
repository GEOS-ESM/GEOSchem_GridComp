!------------------------------------------------------------------------------
! NASA GSFC - SIVO Code 610.3
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: GmiSolverInterface_mod
!
      module GmiSolverInterface_mod
!
! !USES:
      use GmiArrayBundlePointer_mod, only : t_GmiArrayBundle
      use GmiSavedVariables_mod,     only : t_ChemistrySaved
!
      implicit none
!
      private
!
! !PUBLIC MEMBER FUNCTIONS:
      public  :: Update_Smv2chem
!
#     include "GmiParameters.h"
#     include "gmi_phys_constants.h"
#     include "gmi_time_constants.h"
!
!EOP
!------------------------------------------------------------------------------
      contains
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Update_Smv2chem
!
! !INTERFACE:
!
      subroutine Update_Smv2chem (savedVars, chemintv, chem_mecha, surfEmissForChem,      &
     &                 humidity, qjgmi, qkgmi, press3e, pctm2, kel,  &
     &                 concentration, pr_diag, pr_qqjk, pr_smv2, &
     &                 do_smv_reord, do_synoz, do_qqjk_inchem,    &
     &                 do_semiss_inchem, imgas_num, initrogen_num, ioxygen_num,   &
     &                 isynoz_num, yda, qqkda, qqjda, pr_nc_period, tdt,       &
     &                 chem_mask_klo, chem_mask_khi, loc_proc, synoz_threshold,&
     &                 ilong, ilat, ivert, itloop, i1, i2, ju1, j2, k1, k2,    &
     &                 ilo, ihi, julo, jhi, num_molefrac, num_qjo,  &
     &                 num_qks, num_qjs, num_active, num_species, rootProc)
!
      implicit none
!
! !INPUT PARAMETERS:
      logical, intent(in) :: pr_diag, rootProc
      logical, intent(in) :: pr_qqjk
      logical, intent(in) :: pr_smv2
      logical, intent(in) :: do_smv_reord
      logical, intent(in) :: do_synoz
      logical, intent(in) :: do_qqjk_inchem
      logical, intent(in) :: do_semiss_inchem
      integer, intent(in) :: imgas_num, initrogen_num, ioxygen_num, isynoz_num
      integer, intent(in) :: ilong, ilat, ivert, itloop
      integer, intent(in) :: i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi
      integer, intent(in) :: chem_mask_klo, chem_mask_khi
      integer, intent(in) :: num_molefrac
      integer, intent(in) :: num_qjo, num_qks, num_qjs
      integer, intent(in) :: num_active, num_species
      integer, intent(in) :: loc_proc
      real*8 , intent(in) :: pr_nc_period
      real*8 , intent(in) :: synoz_threshold
                             ! chemistry time step (s)
      real*8 , intent(in) :: chemintv
      real*8 , intent(in) :: tdt
                             ! specific humidity (g/kg)
      real*8 , intent(in) :: humidity(i1:i2, ju1:j2, k1:k2)
                             ! CTM surface pressure at t1+tdt (mb)
      real*8 , intent(in) :: pctm2(ilo:ihi, julo:jhi)
      real*8 , intent(in) :: press3e(ilo:ihi, julo:jhi, k1-1:k2)
                             ! temperature (degK)
      real*8 , intent(in) :: kel  (ilo:ihi, julo:jhi, k1:k2)
      CHARACTER(LEN=*), INTENT(IN) :: chem_mecha
                             ! photolysis rate constants (s^-1)
      type (t_GmiArrayBundle), intent(in) :: qjgmi(num_qjo)
                             ! thermal rate constants (units vary)
      type (t_GmiArrayBundle), intent(in) :: qkgmi(num_qks)

      real*8  :: surfEmissForChem(i1:i2, ju1:j2, num_species)
!
! !INPUT/OUTPUT PARAMETERS:
      type(t_ChemistrySaved), intent(inOut) :: savedVars
      real*8, intent(inOut) :: qqjda(i1:i2, ju1:j2, k1:k2, num_qjs)
      real*8, intent(inOut) :: qqkda(i1:i2, ju1:j2, k1:k2, num_qks)
      real*8, intent(inOut) :: yda  (i1:i2, ju1:j2, k1:k2, num_active)
                             ! species concentration, known at zone centers (mixing ratio)
      type (t_GmiArrayBundle), intent(inOut) :: concentration(num_species)
!
! !DESCRIPTION:
!   This routine updates the Smvgear2 chemistry.
!
! !LOCAL VARIABLES:
      character (len=MAX_LENGTH_FILE_NAME) :: smv_filnam
      logical :: do_cell_chem(i1:i2, ju1:j2, k1:k2)
      integer :: ic
      real*8, allocatable  :: speciesConc(:,:,:,:)
      real*8, allocatable  :: loc_qjgmi(:,:,:,:)
      real*8, allocatable  :: loc_qkgmi(:,:,:,:)
!
! !REVISION HISTORY:
!   Initial code.
!
!EOP
!------------------------------------------------------------------------------
!BOC

      if (pr_diag) then
        Write (6,*) 'Update_Smv2chem called by ', loc_proc
      end if


!!     ==========
!      if (first) then
!!     ==========
!
!        first = .false.
!
!        if (pr_smv2) then
!           smv_filnam = 'smv2chem.asc'
!        end if
!
!!       =================
!        call Do_Smv2_Init  &
!!       =================
!     &    (smv_filnam, do_smv_reord, pr_diag, pr_smv2, loc_proc,  &
!     &     itloop, chemintv)
!
!      end if

!     ------------------------------------------------
!     Change units from mixing ratio to concentration.
!     ------------------------------------------------
      do ic = 1, num_molefrac

        concentration(ic)%pArray3D(:,:,:) =  &
     &    concentration(ic)%pArray3D(:,:,:) * concentration(imgas_num)%pArray3D(:,:,:)

      end do


      if (initrogen_num /= 0) then
        concentration(initrogen_num)%pArray3D(:,:,:) =  &
     &    MXRN2 * concentration(imgas_num)%pArray3D(:,:,:)
      end if

      if (ioxygen_num   /= 0) then
        concentration(ioxygen_num)%pArray3D(:,:,:)   =  &
     &    MXRO2 * concentration(imgas_num)%pArray3D(:,:,:)
      end if

      do_cell_chem(:,:,:) = .false.

      do_cell_chem(:,:,chem_mask_klo:chem_mask_khi) = .true.

      ic = INDEX(TRIM(chem_mecha), 'strat_trop')

      IF (ic == 0 .AND. do_synoz) THEN

        where (concentration(isynoz_num)%pArray3D(:,:,:) > synoz_threshold)
          do_cell_chem(:,:,:) = .false.
        end where

      END IF

       allocate(speciesConc(i1:i2, ju1:j2, k1:k2, num_species))
       do ic = 1, num_species
          speciesConc(:,:,:,ic) = concentration(ic)%pArray3D(:,:,:)
       end do

       allocate(loc_qjgmi(i1:i2, ju1:j2, k1:k2, num_qjo))
       do ic = 1, num_qjo
          loc_qjgmi(:,:,:,ic) = qjgmi(ic)%pArray3D(:,:,:)
       end do

       allocate(loc_qkgmi(i1:i2, ju1:j2, k1:k2, num_qks))
       do ic = 1, num_qks
          loc_qkgmi(:,:,:,ic) = qkgmi(ic)%pArray3D(:,:,:)
       end do

!     IF(rootProc) PRINT *,"Calling Do_Smv2_Solver with do_semiss_inchem=",do_semiss_inchem
!     ===================
      call Do_Smv2_Solver  &
!     ===================
     &  (savedVars, do_qqjk_inchem, do_semiss_inchem, pr_diag, pr_qqjk, pr_smv2,  &
     &   loc_proc, ilat, ilong, ivert, itloop, pr_nc_period, tdt,  &
     &   do_cell_chem, loc_qkgmi, loc_qjgmi, surfEmissForChem, speciesConc, &
     &   yda, qqkda, qqjda, loc_qkgmi, loc_qjgmi, &
     &   i1, i2, ju1, j2, k1, k2, &
     &   num_qjo, num_qks, num_qjs, num_active)

       do ic = 1, num_species
          concentration(ic)%pArray3D(:,:,:) = speciesConc(:,:,:,ic)
       end do
       deallocate(speciesConc)

!     -----------------------------------------------------
!     Change units from concentration back to mixing ratio.
!     -----------------------------------------------------

      do ic = 1, num_molefrac

        concentration(ic)%pArray3D(:,:,:) =  &
     &    concentration(ic)%pArray3D(:,:,:) / concentration(imgas_num)%pArray3D(:,:,:)

      end do

      deallocate(loc_qkgmi)
      deallocate(loc_qjgmi)

      return

      end subroutine Update_Smv2chem
!EOC
!-------------------------------------------------------------------------
      end module GmiSolverInterface_mod
