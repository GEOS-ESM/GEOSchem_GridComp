!=======================================================================
!
! $Id$
!
! CODE DEVELOPER
!   Peter Connell, LLNL
!   connell2@llnl.gov
!
! FILE
!   setkin_group_specifics.h
!
! DESCRIPTION
!   This include file contains information about grouping species for
!   transport purposes, that is specific to the mechanism. In general,
!   this file can not be generated automatically, but must be manually
!   produced.
!
!  Chemistry input file:    GMIS2 4:59 PM 7/30/2002
!  Reaction dictionary:     JPL00_reactions.db
!  Setkin files generated:  Mon Jan 13 22:23:01 2003
!
!=======================================================================

      qq2(:,:,:)                   = 0.0d0

      max_bry_adjust               = 0.0d0
      max_cly_adjust               = 0.0d0
      max_nox_adjust               = 0.0d0

      select case (ig)
        case (1)

!.... Bry

          do im                    = 1, MAXGRP_ELEM
            imsgrp                 = sgrp_elem_map(im,ig)
            if (imsgrp > 0)  &
     &        qq2(:,:,:)           = qq2(:,:,:) +  &
     &                               const(:,:,:,imsgrp) *  &
     &                               sgrp_fac(im,ig)
          end do

          group_factor(:,:,:)      = qq1(i1:i2,ju1:j2,:) / qq2(:,:,:)

          do im                    = 1, MAXGRP_ELEM
            imsgrp                 = sgrp_elem_map(im,ig)
            if (imsgrp > 0)  &
     &        const(:,:,:,imsgrp)  = const(:,:,:,imsgrp) *  &
     &                              group_factor(:,:,:)
          end do

        case (2)

!.... Cly

        do im                      = 1, MAXGRP_ELEM
          imsgrp                   = sgrp_elem_map(im,ig)
          if (imsgrp > 0)  &
     &      qq2(:,:,:)             = qq2(:,:,:) +  &
     &                               const(:,:,:,imsgrp) *  &
     &                               sgrp_fac(im,ig)
        end do

        group_factor(:,:,:)        = qq1(i1:i2,ju1:j2,:) / qq2(:,:,:)

        do im                      = 1, MAXGRP_ELEM
          imsgrp                   = sgrp_elem_map(im,ig)

!.... Exclude BrCl from adjustment

          if ((imsgrp > 0) .and. (imsgrp /= IBRCL))  &
     &      const(:,:,:,imsgrp)    = const(:,:,:,imsgrp) *  &
     &                               group_factor(:,:,:)
        end do

!.... Adjust Cl2 in place of BrCl

        group_adjust(:,:,:)        = (group_factor(:,:,:) - 1.0d0) *  &
     &                               0.5d0 * const(:,:,:,IBRCL)

        const(:,:,:,ICL2)          = const(:,:,:,ICL2) +  &
     &                               group_adjust(:,:,:)

        if (Any(const(:,:,:,ICL2) < 0.0d0)) then

!.... In locations where Cl2 is too small to absorb the entire BrCl
!.... decrement, extend it to OClO

          group_adjust(:,:,:)      = 0.0d0
          Where (const(:,:,:,ICL2) < 0.0d0)
            group_adjust(:,:,:)    = 2.0d0 * const(:,:,:,ICL2)
            const(:,:,:,ICL2)      = 0.0d0
            const(:,:,:,IOCLO)     = const(:,:,:,IOCLO) +  &
     &                               group_adjust(:,:,:)
          end where

        end if

        if (Any(const(:,:,:,IOCLO) < 0.0d0)) then

!.... In locations where OClO is too small to absorb the BrCl
!.... residual decrement,extend it to Cl2O2

          group_adjust(:,:,:)      = 0.0d0
          Where (const(:,:,:,IOCLO) < 0.0d0)
            group_adjust(:,:,:)    = 0.5d0 * const(:,:,:,IOCLO)
            const(:,:,:,IOCLO)     = 0.0d0
            const(:,:,:,ICL2O2)    = const(:,:,:,ICL2O2) +  &
     &                               group_adjust(:,:,:)
          end where

        end if

        if (Any(const(:,:,:,ICL2O2) < 0.0d0)) then

!.... In locations where Cl2O2 is too small to absorb the BrCl
!.... residual decrement,extend it to ClO

          group_adjust(:,:,:)      = 0.0d0
          Where (const(:,:,:,ICL2O2) < 0.0d0)
            group_adjust(:,:,:)    = 2.0d0 * const(:,:,:,ICL2O2)
            const(:,:,:,ICL2O2)    = 0.0d0
            const(:,:,:,ICLO)      = const(:,:,:,ICLO) +  &
     &                               group_adjust(:,:,:)
          end where

        end if

        if (Any(const(:,:,:,ICLO) < 0.0d0)) then

!.... In locations where ClO is too small to absorb the BrCl
!.... residual decrement,extend it to HCl

          group_adjust(:,:,:)      = 0.0d0
          Where (const(:,:,:,ICLO) < 0.0d0)
            group_adjust(:,:,:)    = const(:,:,:,ICLO)
            const(:,:,:,ICLO)      = 0.0d0
            const(:,:,:,IHCL)      = const(:,:,:,IHCL) +  &
     &                               group_adjust(:,:,:)
          end where

          if (Any(const(:,:,:,IHCL) < 0.0d0)) then
            max_cly_adjust         = min(max_cly_adjust,  &
     &                                   Minval(group_factor - 1.0d0,  &
     &                                   MASK=((group_factor - 1.0d0)  &
     &                                         < 0.0d0)))
            Write (6,*) 'BrCl adjustment drove HCl negative',  &
     &        max_cly_adjust, loc_proc
            Write (6,*) Minloc(const(:,:,:,IHCL))
          end if

        end if

!.... Account for changes in NOy reservoir ClONO2

        group_adjust(:,:,:)       = const(:,:,:,ICLONO2) *  &
     &                              (1.0d0 - group_factor(:,:,:)) /  &
     &                              group_factor(:,:,:)

!.... Park difference in N2O5

        const(:,:,:,IN2O5)        = const(:,:,:,IN2O5) +  &
     &                              0.5d0 * group_adjust(:,:,:)

        if (Any(const(:,:,:,IN2O5) < 0.0d0)) then

!.... Park remainder in NOx and then HNO3 if necessary

          group_adjust(:,:,:)     = 0.0d0
          Where (const(:,:,:,IN2O5) < 0.0d0)
            group_adjust(:,:,:)   = 2.0d0 * const(:,:,:,IN2O5)
            const(:,:,:,IN2O5)        = 0.0d0
            const(:,:,:,INO)      = const(:,:,:,INO) +  &
     &                              group_adjust(:,:,:) *  &
     &                              (const(:,:,:,INO) /  &
     &                              (const(:,:,:,INO) + const(:,:,:,INO2)))
            const(:,:,:,INO2)      = const(:,:,:,INO2) +  &
     &                              group_adjust(:,:,:) *  &
     &                              (const(:,:,:,INO2) /  &
     &                              (const(:,:,:,INO) + const(:,:,:,INO2)))
          end where

          if (Any(const(:,:,:,INO) < 0.0d0)) then

              group_adjust(:,:,:) = 0.0d0
            Where (const(:,:,:,INO) < 0.0d0)
              group_adjust(:,:,:) = const(:,:,:,INO)
              const(:,:,:,INO)    = 0.0d0
              const(:,:,:,IHNO3)  = const(:,:,:,IHNO3) +  &
     &                              group_adjust(:,:,:)
            end where

          endif

          if (Any(const(:,:,:,INO2) < 0.0d0)) then

              group_adjust(:,:,:) = 0.0d0
            Where (const(:,:,:,INO2) < 0.0d0)
              group_adjust(:,:,:) = const(:,:,:,INO2)
              const(:,:,:,INO2)    = 0.0d0
              const(:,:,:,IHNO3)  = const(:,:,:,IHNO3) +  &
     &                              group_adjust(:,:,:)
            end where

          endif

          if (Any(const(:,:,:,IHNO3) < 0.0d0)) then
            max_nox_adjust        = min(max_nox_adjust,  &
     &                                  minval(group_adjust(:,:,:)))
            Write (6,*) 'ClONO2 adjustment drove NOy negative',  &
     &        max_nox_adjust, loc_proc
          end if

          end if

      end select

