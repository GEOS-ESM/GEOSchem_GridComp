!------------------------------------------------------------------------------
! NASA GSFC - SIVO Code 610.3
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: GmiSurfaceEmissionInChemistry_mod
!
      module GmiSurfaceEmissionInChemistry_mod
!
! !USES:
      use GmiTimeControl_mod, only : GmiSplitDateTime
      use GmiArrayBundlePointer_mod, only : t_GmiArrayBundle

      implicit none
!
      private
!
! !PUBLIC MEMBER FUNCTIONS:
      public  :: updateSurfEmissionInChemistry
!
#     include "gmi_phys_constants.h"
#     include "gmi_time_constants.h"
#     include "gmi_emiss_constants.h"
!
! !DESCRIPTION:
! Module for computing surface emissions for Chemistry.
!
! !AUTHOR:
! Jules Kouatchou, NASA/GSFC, Jules.Kouatchou-1@nasa.gov
!
!EOP
!------------------------------------------------------------------------------
      CONTAINS
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: updateSurfEmissionInChemistry
!
! !INTERFACE:
!
      subroutine updateSurfEmissionInChemistry (pr_surf_emiss, pr_emiss_3d,    &
     &                 emiss_isop, emiss_monot, emiss_nox, do_ShipEmission,    &
     &                 emiss_hno3, emiss_o3, ihno3_num, io3_num, mcor,         &
     &                 surf_emiss_out, surf_emiss_out2, emiss_3d_out,          &
     &                 emissionArray, surfEmissForChem,     &
     &                 gridBoxHeight, emiss_timpyr, num_emiss, emiss_opt,      &
     &                 emiss_map, tdt, nymd, ico_num, ino_num, ipropene_num,   &
     &                 iisoprene_num, mw, pr_diag, loc_proc, i1, i2, ju1, j2,  &
     &                 k1, k2, ilo, ihi, julo, jhi, num_species)
!
      implicit none
!
! !INPUT PARAMETERS:
      logical, intent(in) :: pr_surf_emiss, do_ShipEmission
      logical, intent(in) :: pr_diag, pr_emiss_3d
      integer, intent(in) :: loc_proc, i1, i2, ju1, j2, k1, k2
      integer, intent(in) :: ilo, ihi, julo, jhi, num_species
      integer, intent(in) :: num_emiss, emiss_timpyr, emiss_opt
      integer, intent(in) :: emiss_map(num_species)
                             ! model time step (s)
      real*8 , intent(in) :: tdt
                             !  array of species' molecular weights (g/mol)
      real*8 , intent(in) :: mw(num_species)
      integer, intent(in) :: nymd
                             ! indices for CO, NO, propene, isoprene, O3, HNO3 
                             ! in the species concentration array
      integer, intent(in) :: ico_num, ino_num, ipropene_num, iisoprene_num
      integer, intent(in) :: io3_num, ihno3_num
                             ! area of grid box  (m^2)
      real*8 , intent(in) :: mcor           (i1:i2, ju1:j2)
                             ! isoprene emissions (kg/g)
      real*8 , intent(in) :: emiss_isop     (i1:i2, ju1:j2)
                             ! monoterpene emissions (kg/g)
      real*8 , intent(in) :: emiss_monot    (i1:i2, ju1:j2)
                             ! NOx emissions (kg/g)
      real*8 , intent(in) :: emiss_nox      (i1:i2, ju1:j2)
                             ! HNO3 emissions (kg/g)
      real*8 , intent(in) :: emiss_hno3     (i1:i2, ju1:j2)
                             ! O3 emissions (kg/g)
      real*8 , intent(in) :: emiss_o3    (i1:i2, ju1:j2)
                             ! height of each grid box (m)
      real*8 , intent(in) :: gridBoxHeight(i1:i2, ju1:j2, k1:k2)
! !OUTPUT PARAMETERS:
      ! Surface emission
      real*8 , intent(out) :: surfEmissForChem (i1:i2, ju1:j2, num_species)
!
! !INPUT/OUTPUT PARAMETERS:
      ! Surface emission diagnistics
      real*8 , intent(inOut) :: surf_emiss_out (i1:i2, ju1:j2, num_species)
      real*8 , intent(inOut) :: surf_emiss_out2 (i1:i2, ju1:j2, 6)
      real*8 , intent(inOut) :: emiss_3d_out (i1:i2, ju1:j2, k1:k2,num_species)
      ! array of emissions (kg/s)
      type(t_GmiArrayBundle), intent(inOut) :: emissionArray(num_emiss)
!
! !DEFINED PARAMETERS:
      real*8, parameter :: CMPM3 = CMPM * CMPM * CMPM
!
! !DESCRIPTION:
! Updates species concentration based on surface emissions, and does it in 
! the chemistry solver.
!
! !LOCAL VARIABLES:
      integer :: ic, icx, ij, il, ik
      integer :: idumday, idumyear, imon, inum, month
      real*8  :: box_height(i1:i2, ju1:j2)
      real*8  :: conv_emiss(i1:i2, ju1:j2)
      real*8  :: surf_emiss2 (i1:i2, ju1:j2, 6)
      real*8  :: emiss_3d(i1:i2, ju1:j2, k1:k2, num_species)
!EOP
!------------------------------------------------------------------------------
!BOC
      if (pr_diag) Write (6,*) 'Update_Semiss_Inchem called by ', loc_proc

      if (emiss_timpyr == MONTHS_PER_YEAR) then
         call GmiSplitDateTime (nymd, idumyear, month, idumday)
         imon = month
      else
         imon = 1
      end if

      !---------------------------------------------------------------
      ! Calculate conversion factor to go from kg/box/s to mol/cm^3/s,
      ! but leave out the molecular weight term for each species.
      !---------------------------------------------------------------

      conv_emiss(:,:) = (AVOGAD*GPKG) / mcor(:,:) / CMPM3 / gridBoxheight(:,:,k1)

      !--------------------------------------------------------------
      ! Set surf_emiss; convert units from kg/s to ? for each species
      ! that is emitted.
      !--------------------------------------------------------------

      inum = 0

      surfEmissForChem (:,:,:)   = 0.0d0
      surf_emiss2(:,:,:)   = 0.0d0
      emiss_3d   (:,:,:,:) = 0.0d0

      !===============================
      SPCLOOP: do icx = 1, num_species
      !===============================

         ic = emiss_map(icx)
         if (ic > 0) then
            inum = inum + 1

            surfEmissForChem(:,:,ic) = surfEmissForChem(:,:,ic) +  &
     &          emissionArray(inum)%pArray3D(:,:,k1) * conv_emiss(:,:) / mw(ic)

            emiss_3d(:,:,k1,ic) = emiss_3d(:,:,k1,ic) +  &
     &           emissionArray(inum)%pArray3D(:,:,k1) * conv_emiss(:,:) / mw(ic)
                                   !============
            if (inum == num_emiss) exit SPCLOOP
                                   !============
         end if
      !==============
      end do SPCLOOP
      !==============


      if (emiss_opt == 2) then
         call addHarvardSurfaceEmission (iisoprene_num, ico_num, ipropene_num, &
     &                 ino_num, mw, emiss_isop, emiss_monot, emiss_nox,        &
     &                 do_ShipEmission, emiss_hno3, emiss_o3, ihno3_num,       &
     &                 io3_num, conv_emiss, surfEmissForChem, surf_emiss2, emiss_3d, &
     &                 pr_diag, loc_proc, i1, i2, ju1, j2, k1, k2, num_species)
      end if

      if (pr_surf_emiss) then
         do ic = 1, num_species
            surf_emiss_out(:,:,ic) = surf_emiss_out(:,:,ic) +  &
     &                surfEmissForChem(:,:,ic) / conv_emiss(:,:) * mw(ic) *  &
     &                tdt / mcor(:,:)
         end do

         surf_emiss_out2(:,:,1) = surf_emiss_out2(:,:,1) +  &
     &             surf_emiss2(:,:,1) / conv_emiss(:,:) * mw(ico_num) *  &
     &             tdt / mcor(:,:)

         surf_emiss_out2(:,:,2) = surf_emiss_out2(:,:,2) +  &
     &             surf_emiss2(:,:,2) / conv_emiss(:,:) * mw(ico_num) *  &
     &             tdt / mcor(:,:)

         surf_emiss_out2(:,:,3) = surf_emiss_out2(:,:,3) +  &
     &             surf_emiss2(:,:,3) / conv_emiss(:,:) * mw(ipropene_num) *  &
     &             tdt / mcor(:,:)

         surf_emiss_out2(:,:,4) = surf_emiss_out2(:,:,4) +  &
     &             surf_emiss2(:,:,4) / conv_emiss(:,:) * mw(ino_num) *  &
     &             tdt / mcor(:,:)

         if (do_ShipEmission) then
            surf_emiss_out2(:,:,5) = surf_emiss_out2(:,:,5) +  &
     &                surf_emiss2(:,:,5) / conv_emiss(:,:) * mw(ihno3_num) *  &
     &                tdt / mcor(:,:)

            surf_emiss_out2(:,:,6) = surf_emiss_out2(:,:,6) +  &
     &                surf_emiss2(:,:,6) / conv_emiss(:,:) * mw(io3_num) *  &
     &                tdt / mcor(:,:)
          end if
      end if

      if (pr_emiss_3d) then
         do ic = 1, num_species
            emiss_3d_out(:,:,k1,ic) =  emiss_3d_out(:,:,k1,ic) +  &
     &           emiss_3d(:,:,k1,ic) / conv_emiss(:,:) * mw(ic) *  &
     &            tdt / mcor(:,:)
         end do
      end if

      return

      end subroutine updateSurfEmissionInChemistry
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: addHarvardSurfaceEmission
! 
! !INTERFACE:
!
      subroutine addHarvardSurfaceEmission (iisoprene_num, ico_num,            &
     &                     ipropene_num, ino_num, mw, emiss_isop, emiss_monot, &
     &                     emiss_nox, do_ShipEmission, emiss_hno3, emiss_o3,   &
     &                     ihno3_num, io3_num, conv_emiss, surfEmissForChem,         &
     &                     surf_emiss2, emiss_3d, pr_diag, loc_proc, i1, i2,   &
     &                     ju1, j2, k1, k2, num_species)
!
      implicit none
!
! !INPUT PARAMETERS:
      logical, intent(in) :: pr_diag, do_ShipEmission
      integer, intent(in) :: loc_proc, i1, i2, ju1, j2, k1, k2, num_species
      integer, intent(in) :: iisoprene_num, ico_num, ipropene_num
      integer, intent(in) :: ino_num, io3_num, ihno3_num
      real*8 , intent(in) :: mw(num_species)
      real*8 , intent(in) :: emiss_isop (i1:i2, ju1:j2)
      real*8 , intent(in) :: emiss_monot(i1:i2, ju1:j2)
      real*8 , intent(in) :: emiss_nox  (i1:i2, ju1:j2)
      real*8 , intent(in) :: emiss_hno3 (i1:i2, ju1:j2)
      real*8 , intent(in) :: emiss_o3(i1:i2, ju1:j2)
      real*8 , intent(in) :: conv_emiss (i1:i2, ju1:j2)
!
! !INPUT/OUTPUT PARAMETERS:
      real*8 , intent(inOut) :: surfEmissForChem (i1:i2, ju1:j2, num_species)
      real*8 , intent(inOut) :: surf_emiss2 (i1:i2, ju1:j2, 6)
      real*8 , intent(inOut) :: emiss_3d (i1:i2, ju1:j2, k1:k2, num_species)
!
! !DESCRIPTION: 
! Adds in the Harvard emissions, and is used just prior tp calling the
! Chemistry solver.
!
!EOP
!------------------------------------------------------------------------------
!BOC
      if (pr_diag) Write (6,*) 'Add_Semiss_Harvard_Inchem called by ', loc_proc
      
      if (iisoprene_num > 0) then
         !-----------------------------
         ! Biogenic source of isoprene.
         !-----------------------------
         surfEmissForChem(:,:,iisoprene_num) =  surfEmissForChem(:,:,iisoprene_num) +  &
     &                (emiss_isop(:,:) * conv_emiss(:,:) / mw(iisoprene_num))

         emiss_3d(:,:,1,iisoprene_num) = emiss_3d(:,:,1,iisoprene_num) +  &
     &           emiss_isop(:,:) * conv_emiss(:,:) / mw(iisoprene_num)
      end if

      if (ico_num > 0) then
         !-----------------------------------------------
         ! Biogenic source of CO from methanol oxidation,
         ! scaled from isoprene.
         !-----------------------------------------------

         surfEmissForChem(:,:,ico_num) = surfEmissForChem(:,:,ico_num) +  &
     &             (emiss_isop(:,:) * conv_emiss(:,:) / mw(ico_num)) *  &
     &             ICO_FAC_ISOP

        surf_emiss2(:,:,1) = surf_emiss2(:,:,1) +  &
     &             (emiss_isop(:,:) * conv_emiss(:,:) / mw(ico_num)) *  &
     &             ICO_FAC_ISOP

        emiss_3d(:,:,1,ico_num) = emiss_3d(:,:,1,ico_num) +  &
     &           emiss_isop(:,:) * conv_emiss(:,:) / mw(ico_num) *  &
     &           ICO_FAC_ISOP

        !--------------------------------------------------
        ! Biogenic source of CO from monoterpene oxidation.
        !--------------------------------------------------

        surfEmissForChem(:,:,ico_num) = surfEmissForChem(:,:,ico_num) +  &
     &            (emiss_monot(:,:) *  conv_emiss(:,:) / mw(ico_num)) *  &
     &            ICO_FAC_MONOT

        surf_emiss2(:,:,2) = surf_emiss2(:,:,2) +  &
     &            (emiss_monot(:,:) *  conv_emiss(:,:) / mw(ico_num)) *  &
     &            ICO_FAC_MONOT

        emiss_3d(:,:,1,ico_num) = emiss_3d(:,:,1,ico_num) +  &
     &           emiss_monot(:,:) * conv_emiss(:,:) / mw(ico_num) *  &
     &           ICO_FAC_MONOT

      end if

      if (ipropene_num > 0) then
         !--------------------------------------------------
         ! Biogenic source of propene, scaled from isoprene.
         !--------------------------------------------------

         surfEmissForChem(:,:,ipropene_num) = surfEmissForChem(:,:,ipropene_num)  +  &
     &           (emiss_isop(:,:) * conv_emiss(:,:) / mw(iisoprene_num)) *  &
     &           BIOSCAL * (ATOMSC_PER_MOLECISOP / ATOMSC_PER_MOLECPRPE)

         surf_emiss2(:,:,3) =  surf_emiss2(:,:,3)  +  &
     &           (emiss_isop(:,:) * conv_emiss(:,:) / mw(iisoprene_num)) *  &
     &           BIOSCAL * (ATOMSC_PER_MOLECISOP / ATOMSC_PER_MOLECPRPE)

         emiss_3d(:,:,1,ipropene_num) = emiss_3d(:,:,1,ipropene_num) +  &
     &           emiss_isop(:,:) * conv_emiss(:,:) / mw(iisoprene_num) *  &
     &           BIOSCAL * (ATOMSC_PER_MOLECISOP / ATOMSC_PER_MOLECPRPE)

      endif

      if (ino_num > 0) then
         !--------------------
         ! Soil source of NOx.
         !--------------------
         surfEmissForChem(:,:,ino_num) = surfEmissForChem(:,:,ino_num) +  &
     &          (emiss_nox(:,:) * conv_emiss(:,:) / mw(ino_num))

         surf_emiss2(:,:,4) = surf_emiss2(:,:,4) +  &
     &          (emiss_nox(:,:) * conv_emiss(:,:) / mw(ino_num))

         emiss_3d(:,:,1,ino_num) = emiss_3d(:,:,1,ino_num) +  &
     &           emiss_nox(:,:) * conv_emiss(:,:) / mw(ino_num)
      end if

      if (do_ShipEmission) then

         if (ihno3_num > 0) then
            ! --------------
            ! Source of HNO3
            ! --------------
            surf_emiss2(:,:,5) =  surf_emiss2(:,:,5) +  &
     &          (emiss_hno3(:,:) * conv_emiss(:,:) / mw(ihno3_num))
         end if

         if (io3_num > 0) then
            ! ------------
            ! Source of O3
            ! ------------
            surf_emiss2(:,:,6) =  surf_emiss2(:,:,6) +  &
     &          (emiss_o3(:,:) * conv_emiss(:,:) / mw(io3_num))
         end if

      end if

      return

      end subroutine addHarvardSurfaceEmission
!EOC
!------------------------------------------------------------------------------
      end module GmiSurfaceEmissionInChemistry_mod
