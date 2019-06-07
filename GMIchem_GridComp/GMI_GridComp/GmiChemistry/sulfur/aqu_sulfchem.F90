!=============================================================================
!
! $Id$
!
! CODE DEVELOPER
!   originated from GRANTOUR programmed by John Walton
!   modified by Xiaohong Liu
!
! FILE
!   aqu_sulfchem.F
!
! ROUTINES
!   Do_Aqu_Sulfchem
!
! HISTORY
!   - December 8, 2005 * Bigyani Das
!     Added variables sfso4, snso4 for primary so4 emission
!     as 2.5% of SO2 emission and also in the SO4 chemical source
!     and sink terms add sfso4, snso4.  used dms_no3 array as the 
!     index for primary so4 emission budget
!     Added controlling variables do_aerocom and do_emiss_dust

!=============================================================================

!-----------------------------------------------------------------------------
!
! ROUTINE
!   Do_Aqu_Sulfchem (chem30)
!
! DESCRIPTION
!
!     ******************************************************************
!     * Uses a simple centered difference scheme to solve sulfur
!     * chemistry.
!     ******************************************************************
!
!     *     ndt          ==  the number of time steps
!     *     ndtmin       ==  the initial number of time steps to
!                            be used in this routine.
!     *     dtchem       ==  model time step in seconds
!     *     dtclc        ==  the time step (s) for the coupling of
!                            the clear & cloudy sky chemistry
!     *     dtaq         ==  working time step (sec)
!     *     acc          ==  max difference permitted between
!                            succesive solutions for the chemisty
!                            to be considered solved.
!     *     ndone        ==  the number of grid points completed
!     *     emin         ==  the minimum concentration (cm-3) to
!                            which a species is permitted to go.
!     *     numcl        ==  no. of grid points in cloud
!     *     numwk        ==  no. of grid points currently being solved
!     *     nfin         ==  -1,  grid points finished
!                             1,  grid points still integrating
!     *     idxcl        ==  the numbers of the numcl points in cloud
!     *     idxwk        ==  the numbers of the numwk points still
!                            being solved
!
!     *     The following species are categorized by the type
!           of source driving them.
!
!     *     fso4a      ==  H2SO4,  fossil - aqueous
!     *     fso4n      ==  H2SO4,  fossil - clear sky
!     *     fso2       ==  SO2  ,    "
!
!     *     nso4a      ==  H2SO4,  natural - aqueous
!     *     nso4n      ==  H2SO4,     "    - clear sky
!     *     nso2       ==  SO2  ,     "
!     *     ndms       ==  DMS  ,     "
!
!     *     h2o2       ==  H2O2
!
!     ********** output:
!
!     *     dfso4a      ==  H2SO4,  fossil - aqueous
!     *     dfso2       ==  SO2  ,    "
!
!     *     dnso4a      ==  H2SO4,  natural- aqueous
!     *     dnso2       ==  SO2  ,     "
!     *     dndms       ==  DMS  ,     "
!
!     *     dh2o2       ==  H2O2
!
!
!     ********** notes:
!
!     *     The calculation is performed with a centered difference,
!           with and initial predictor step followed by a corrector.
!
!     *     This routine is designed so that the total number of
!           grid points being integrated, numwk at each time, are
!           contiguous.  That is, there are no gaps between the slots
!           being integrated.
!
!     *     The routine is entered with the contribution of species
!           change due to clear sky reactions, stored in dcp.  It
!           computes the aqueous contribution and adds it to that for
!           clear sky.
!
!----------------------------------------------------------------------

      subroutine Do_Aqu_Sulfchem  &
     &  (acc, dtaq, dtchem, dtcl, dtclc,  &
     &   itloop, numcl, idxcl, semiss,  &
     &   massc, pr_sulf_src, do_aerocom, do_dust_emiss,  &
     &   dms_oh, dms_no3, so2_oh, so2_h2o2, so2_o3,  &
     &   cp, dcp, cfp, zmp, a, loc_proc)

      implicit none

#     include "setkin_par.h"
#     include "sulfchem.h"

!     ----------------------
!     GMIMOD-related values.
!     ----------------------

      integer :: itloop, numcl
      integer :: idxcl(itloop)

      real*8  &
     &  acc,  &
     &  dtaq,    & ! time step used in integration of the aqueous equations in seconds
     &  dtchem,  & ! model time step in seconds
     &  dtcl,    & ! time step between mixing calculations in seconds
     &  dtclc  ! time step for clear air chemistry

      real*8  :: cp    (itloop, NSP)
      real*8  :: dcp   (itloop, NSP)
      real*8  :: semiss(itloop, NSP)
      real*8  :: a     (itloop, 9)
      real*8  :: cfp   (itloop)
      real*8  :: zmp   (itloop)

      real*8  :: massc   (itloop)
      logical :: pr_sulf_src
      logical :: do_aerocom
      logical :: do_dust_emiss
      real*8  :: dms_oh  (itloop)
      real*8  :: dms_no3 (itloop)
      real*8  :: so2_oh  (itloop)
      real*8  :: so2_h2o2(itloop)
      real*8  :: so2_o3  (itloop)

      integer :: loc_proc

! -------------------------------------------------------------
!  *** LOCAL  VARIABLES WERE DECLARED HERE.
! -------------------------------------------------------------

      real*8  ndtmin, emin

      data ndtmin / 24.0d0 /
      data emin / 100.0d0 /

      integer :: nfin(numcl)
      integer :: idxwk(numcl)

      real*8 :: h2o2(itloop)
      real*8 :: ndms(itloop)
      real*8 :: nso2(itloop)
      real*8 :: fso2(itloop)

      real*8 :: dh2o2 (itloop)
      real*8 :: dndms (itloop)
      real*8 :: dnso2 (itloop)
      real*8 :: dfso2 (itloop)
      real*8 :: dnso4a(itloop)
      real*8 :: dfso4a(itloop)

      real*8 :: sfso2(itloop)
      real*8 :: snso2(itloop)
      real*8 :: sndms(itloop)

      real*8 :: sfso4(itloop)    ! primary so4
      real*8 :: snso4(itloop)    ! primary so4
      real*8 :: zsfso4, zsnso4

      real*8 :: h2o2w1(numcl)
      real*8 :: fso2w1(numcl)
      real*8 :: nso2w1(numcl)
      real*8 :: ndmsw1(numcl)
      real*8 :: h2o2w2(numcl)
      real*8 :: fso2w2(numcl)
      real*8 :: nso2w2(numcl)
      real*8 :: ndmsw2(numcl)

      real*8 :: dh2o21(numcl)
      real*8 :: dfso21(numcl)
      real*8 :: dnso21(numcl)
      real*8 :: dndms1(numcl)
      real*8 :: dh2o22(numcl)
      real*8 :: dfso22(numcl)
      real*8 :: dnso22(numcl)
      real*8 :: dndms2(numcl)

      real*8 :: rh2o2(numcl)
      real*8 :: rfso2(numcl)
      real*8 :: rnso2(numcl)

      real*8 :: h2o2o(numcl)
      real*8 :: fso2o(numcl)
      real*8 :: nso2o(numcl)
      real*8 :: ndmso(numcl)

      real*8 :: zdms_oh  (itloop)
      real*8 :: zdms_no3 (itloop)
      real*8 :: zso2_oh  (itloop)
      real*8 :: zso2_h2o2(itloop)
      real*8 :: zso2_o3  (itloop)

      nfin = 0.0d0
      idxwk = 0.0d0

      h2o2 = 0.0d0
      ndms = 0.0d0
      nso2 = 0.0d0
      fso2 = 0.0d0

      dh2o2  = 0.0d0
      dndms  = 0.0d0
      dnso2  = 0.0d0
      dfso2  = 0.0d0
      dnso4a = 0.0d0
      dfso4a = 0.0d0

      sfso2  = 0.0d0
      snso2  = 0.0d0
      sndms  = 0.0d0

      sfso4  = 0.0d0
      snso4  = 0.0d0

      h2o2w1 = 0.0d0
      fso2w1 = 0.0d0
      nso2w1 = 0.0d0
      ndmsw1 = 0.0d0
      h2o2w2 = 0.0d0
      fso2w2 = 0.0d0
      nso2w2 = 0.0d0
      ndmsw2 = 0.0d0

      dh2o21 = 0.0d0
      dfso21 = 0.0d0
      dnso21 = 0.0d0
      dndms1 = 0.0d0
      dh2o22 = 0.0d0
      dfso22 = 0.0d0
      dnso22 = 0.0d0
      dndms2 = 0.0d0

      rh2o2 = 0.0d0
      rfso2 = 0.0d0
      rnso2 = 0.0d0

      h2o2o = 0.0d0
      fso2o = 0.0d0
      nso2o = 0.0d0
      ndmso = 0.0d0

      zdms_oh   = 0.0d0
      zdms_no3  = 0.0d0
      zso2_oh   = 0.0d0
      zso2_h2o2 = 0.0d0
      zso2_o3   = 0.0d0


      h2o2(:) = cp(:, IH2O2)
      ndms(:) = cp(:, INDMS)
      nso2(:) = cp(:, INSO2)
      fso2(:) = cp(:, IFSO2)

      if (do_aerocom) then
         sfso2(:) = semiss(:, IFSO2) * 0.975d0
         snso2(:) = semiss(:, INSO2) * 0.975d0
         sndms(:) = semiss(:, INDMS)

         sfso4(:) = semiss(:, IFSO2) * 0.025d0
         snso4(:) = semiss(:, INSO2) * 0.025d0
      else
         sfso2(:) = semiss(:, IFSO2)
         snso2(:) = semiss(:, INSO2)
         sndms(:) = semiss(:, INDMS)   ! no scaling
!     sndms(:) = semiss(:, INDMS) /2.0d0  ! scale by 2 to reduce DMS emission
      endif

      dh2o2(:) = dcp(:, IH2O2)
      dndms(:) = dcp(:, INDMS)
      dfso2(:) = dcp(:, IFSO2)
      dnso2(:) = dcp(:, INSO2)


!     *--------------------------------------------------
!     * Initialize parameters for this set of grid-points
!     *--------------------------------------------------

!.... ndone, the number of grid points finished.  When this = numcl
!            the routine is finished

      ndone = 0

!.... nfin =  1,  if the grid-point is still being integrated
!          = -1,  if the grid-point is done

!.... idxwk, a working index containing the numbers of grid-points to be
!            integrated. This will be initially set = idxcl, the numbers
!            of the grid-points within clouds.
!            idxwk and its length, numwk, will change as the integration
!            advances.  idxcl and numcl, on the other hand, will remain
!            unchanged.

      do i = 1, numcl
        nfin(i)  = 1
        idxwk(i) = idxcl(i)
      enddo

!.... initially the grid-point index and number is the same as
!     you come in with

      numwk = numcl

!     The initial dtaq will be dtclc / ndtmin.

      ndt = ndtmin


!     *-----------------------------------------------------------------
!     * This is the point to which one returns when a new time step,
!       dtaq, integration is to be started.
!     *-----------------------------------------------------------------

!     ========
  100 continue
!     ========

!.... set dtaq

      dtaq = dtclc / float( ndt )


!     ******************************************************************
!       If dtaq < 0.1 seconds the method has failed to converge.
!       An implicit calculation with a 1 second time step will be used.
!       The relation solved, for species f, is
!
!         [ f(t+dt) - f(t) ] / dt
!         = Source + Production( f' ) - L( f' ) [ f(t+dt) + f(t) ] / 2
!
!       where f' refers to a species list that does not include f.
!     ******************************************************************

      if( dtaq .lt. 0.1 ) then

        write(6, 10)  numwk, 2. * dtaq, loc_proc
   10   format( / 'chem30: ', i6, ' grid cells failed to converge with',  &
     &            ' dtaq = ', e12.4, ' sec.',  &
     &            'This loc_proc is ', i4,  &
     &          / 'Set dtaq = 1 sec. and ',  &
     &            'do an implicit integration for these grid cells.' )

        dtaq = 1.
        ndt  = ( dtclc + 0.1 ) / dtaq


        do 920 i = 1, numwk

          iw = idxwk(i)

!     Species
          zh2o2 = h2o2( iw )
          zfso2 = fso2( iw )
          znso2 = nso2( iw )
          zndms = ndms( iw )
!         zbso2 = bso2( iw )

!     Sources
          zsfso2 = sfso2( iw )
          zsndms = sndms( iw )
          zsnso2 = snso2( iw )
!         zsbso2 = sbso2( iw )

          if (do_aerocom) then
            zsfso4 = sfso4( iw )
            zsnso4 = snso4( iw )
          end if
!     a's
          za1 = a( iw, 1 )
          za2 = a( iw, 2 )
          za3 = a( iw, 3 )
          za4 = a( iw, 4 )
          za5 = a( iw, 5 )
          za6 = a( iw, 6 )
          za7 = a( iw, 7 )
          za8 = a( iw, 8 )
          za9 = a( iw, 9 )

!     Sulf budget
          zdms_oh  (iw) = 0.0d0
          zdms_no3 (iw) = 0.0d0
          zso2_oh  (iw) = 0.0d0
          zso2_h2o2(iw) = 0.0d0
          zso2_o3  (iw) = 0.0d0


          do 900 m = 1, ndt

!     Calculate the production and lose terms.

!.... H2O2
            prh2o2 = za3 + za8
            dsh2o2 = ( za4 + za5  &
     &               + za6 * ( zfso2 + znso2 ) ) * 0.5 * dtaq

!.... fossil
            prfso2 = zsfso2
            dsfso2 = ( za2 + za7 + za6 * zh2o2 ) * 0.5 * dtaq

!.... natural
            prndms = zsndms
            dsndms = ( za1 + za9 ) * 0.5 * dtaq

            prnso2 = zsnso2 + ( za1 + za9 ) * zndms
            dsnso2 = ( za2 + za7 + za6 * zh2o2 ) * 0.5 * dtaq

!.... biomass
!           prbso2 = zsbso2
!           dsbso2 = ( za2 + za7 + za6 * zh2o2 ) * 0.5 * dtaq
!           dbso21 = sbso2
!    #               - ( a2 + a7
!    #               + a6 * h2o2 ) * bso2

!     Advance species concentrations.

            zh2o2 = ( prh2o2 * dtaq + ( 1. - dsh2o2  ) * zh2o2 )  &
     &            /                   ( 1. + dsh2o2  )
            zfso2 = ( prfso2 * dtaq + ( 1. - dsfso2  ) * zfso2 )  &
     &            /                   ( 1. + dsfso2  )
            zndms = ( prndms * dtaq + ( 1. - dsndms  ) * zndms )  &
     &            /                   ( 1. + dsndms  )
            znso2 = ( prnso2 * dtaq + ( 1. - dsnso2  ) * znso2 )  &
     &            /                   ( 1. + dsnso2  )
!           zbso2 = ( prbso2 * dtaq + ( 1. - dsbso2  ) * zbso2 )
!    &            /                   ( 1. + dsbso2  )

            if(pr_sulf_src) then

              zdms_oh(iw)   = zdms_oh(iw) +  &
     &                        za1 * zndms * dtaq
              if (do_aerocom) then
                 zdms_no3(iw)  = zdms_no3(iw) +  &
     &                        (zsfso4 + zsnso4) * dtaq
              else
                 zdms_no3(iw)  = zdms_no3(iw) +  &
     &                        za9 * zndms * dtaq
              end if
              zso2_oh(iw)   = zso2_oh(iw) +  &
     &                        za2 * (zfso2 + znso2) * dtaq
              zso2_h2o2(iw) = zso2_h2o2(iw) +  &
     &                        za6 * zh2o2 * (zfso2 + znso2) * dtaq
              zso2_o3(iw)   = zso2_o3(iw) +  &
     &                        za7 * (zfso2 + znso2) * dtaq

            end if

  900     continue


!     Include the species change due to aqueous reaction in the
!     total change due to chemistry.

!.... H2O2
          dh2o2( iw )  = dh2o2( iw )  &
     &                 + cfp( iw ) * ( zh2o2 - h2o2( iw ) )

!.... fossil
          dfso2( iw )  = dfso2( iw )  &
     &                 + cfp( iw ) * ( zfso2 - fso2( iw ) )
          if (do_aerocom) then
             dfso4a( iw ) = cfp( iw )  &
     &                 * (   sfso2(iw) * dtclc + sfso4(iw) * dtclc      & ! dtchem
     &                   - ( zfso2 - fso2( iw ) ) )
          else
             dfso4a( iw ) = cfp( iw )  &
     &                 * (   sfso2(iw) * dtclc       & ! dtchem
     &                   - ( zfso2 - fso2( iw ) ) )
          end if

!.... natural
          dnso2( iw )  = dnso2( iw )  &
     &                 + cfp( iw ) * ( znso2 - nso2( iw ) )
          dndms( iw )  = dndms( iw )  &
     &                 + cfp( iw ) * ( zndms - ndms( iw ) )
          dnso4a( iw ) = cfp( iw )  &
     &                 * ( ( snso2(iw) + sndms(iw) ) * dtclc    & ! dtchem
     &                   - ( znso2 - nso2( iw ) )  &
     &                   - ( zndms - ndms( iw ) ) )

!.... biomass
!         dbso2( iw )  = dbso2( iw )
!     &                + cfp( iw ) * ( zbso2 - bso2( iw ) )
!         dbso4a( iw ) = cfp( iw )
!     &                * (   sbso2(iw) * dtclc     ! dtchem
!     &                  - ( zbso2 - bso2( iw ) ) )

  920   continue


!     This completes the calculation using the implicit scheme

!       =========
        go to 787
!       =========

!.... this is the end of the dtaq if statement
!     =====
      endif
!     =====


!     Note, since the only change to H2SO4 is through chemical
!     sources, it can be set = 0 here and the routine will
!     calculate the H2SO4 change due to aqueous reactions.

      do i = 1, numwk

        iw = idxwk(i)

!.... H2O2
        h2o2w1(i) = h2o2( iw )

!.... fossil
        fso2w1(i) = fso2( iw )

!.... natural
        nso2w1(i) = nso2( iw )
        ndmsw1(i) = ndms( iw )

!.... biomass
!       bso2w1(i) = bso2( iw )

      enddo


!     ******************************************************************
!     * Advance the chemistry for each active grid-points.
!     * In what follows, it is necessary to point to the correct a's
!       and sources for each grid point being run.  The idxwk does this.
!     ******************************************************************

!     *--------------------
!     Begin time step loop.
!     *--------------------

      do i = 1, numwk

        iw = idxwk(i)

        zdms_oh  (iw) = 0.0d0
        zdms_no3 (iw) = 0.0d0
        zso2_oh  (iw) = 0.0d0
        zso2_h2o2(iw) = 0.0d0
        zso2_o3  (iw) = 0.0d0

      enddo

      do 200 m = 1, ndt


!     *-----------------------------------------------------------------
!     * Tentative step.
!     * Compute species derivatives based on initial values.
!     *-----------------------------------------------------------------

        do i = 1, numwk

          iw = idxwk(i)

!.... H2O2
          dh2o21(i) = a(iw, 3) + a(iw, 8) - ( a(iw, 4) + a(iw, 5)  &
     &              + a(iw, 6) * ( fso2w1(i) + nso2w1(i) ) )  &
     &              * h2o2w1(i)

!.... fossil
          dfso21(i) = sfso2(iw)  &
     &              - ( a(iw, 2) + a(iw, 7)  &
     &              + a(iw, 6) * h2o2w1(i) ) * fso2w1(i)

!.... natural
          dndms1(i) = sndms(iw) - ( a(iw, 1) + a(iw, 9) ) * ndmsw1(i)

          dnso21(i) = snso2(iw)  &
     &              + ( a(iw, 1) + a(iw, 9) ) * ndmsw1(i)  &
     &              - ( a(iw, 2) + a(iw, 7)  &
     &              + a(iw, 6) * h2o2w1(i) ) * nso2w1(i)

!.... biomass
!         dbso21(i) = sbso2(iw)
!    &              - ( a2(iw) + a7(iw)
!    &              + a6(iw) * h2o2w1(i) ) * bso2w1(i)

        enddo


!     *-----------------------------------------------------------------
!     * Evaluate tentative the species based on initial values.
!     *-----------------------------------------------------------------

        do i = 1, numwk

!.... H2O2
          h2o2w2(i) = max( h2o2w1(i) + dh2o21(i) * dtaq, 0.0d0)

!.... fossil
          fso2w2(i) = max( fso2w1(i) + dfso21(i) * dtaq, 0.0d0)

!.... natural
          ndmsw2(i) = max( ndmsw1(i) + dndms1(i) * dtaq, 0.0d0)
          nso2w2(i) = max( nso2w1(i) + dnso21(i) * dtaq, 0.0d0)

!.... biomass
!         bso2w2(i) = max( bso2w1(i) + dbso21(i) * dtaq, 0.0d0)

        enddo


!     *-----------------------------------------------------------------
!     * Final step.
!     * Compute species derivatives based on the tentative
!       concentrations.
!     *-----------------------------------------------------------------

        do i = 1, numwk

          iw = idxwk(i)

!.... H2O2
          dh2o22(i) = a(iw, 3) + a(iw, 8) - ( a(iw, 4) + a(iw, 5)  &
     &              + a(iw, 6) * ( fso2w2(i) + nso2w2(i) ) )  &
     &              * h2o2w2(i)

!.... fossil
          dfso22(i) = sfso2(iw)  &
     &              - ( a(iw, 2) + a(iw, 7)  &
     &              + a(iw, 6) * h2o2w2(i) ) * fso2w2(i)

!.... natural
          dndms2(i) = sndms(iw) - ( a(iw, 1) + a(iw, 9) ) * ndmsw2(i)

          dnso22(i) = snso2(iw)  &
     &              + ( a(iw, 1) + a(iw, 9) ) * ndmsw2(i)  &
     &              - ( a(iw, 2) + a(iw, 7)  &
     &              + a(iw, 6) * h2o2w2(i) ) * nso2w2(i)

!.... biomass
!         dbso22(i) = sbso2(iw)
!    &              - ( a2(iw) + a7(iw)
!    &                + a6(iw) * h2o2w2(i) ) * bso2w2(i)

        enddo

!c sulf budget
        if(pr_sulf_src) then

          do i = 1, numwk

            iw = idxwk(i)

              zdms_oh(iw) = zdms_oh(iw) +  &
     &          a(iw, 1) * 0.5 *  &
     &          (ndmsw1(i) + ndmsw2(i)) * dtaq
              if (do_aerocom) then
                zdms_no3(iw) = zdms_no3(iw) +  &
!    &           a(iw, 9) * 0.5 *
!    &           (ndmsw1(i) + ndmsw2(i)) * dtaq
     &           (sfso4(iw) + snso4(iw)) * dtaq
              else
                zdms_no3(iw) = zdms_no3(iw) +  &
     &           a(iw, 9) * 0.5 *  &
     &           (ndmsw1(i) + ndmsw2(i)) * dtaq
              end if
              zso2_oh(iw) = zso2_oh(iw) +  &
     &          a(iw, 2) * 0.5 *  &
     &          (fso2w1(i) + nso2w1(i) + fso2w2(i) + nso2w2(i)) * dtaq
              zso2_h2o2(iw) = zso2_h2o2(iw) +  &
     &          a(iw, 6) * 0.5 *  &
     &          (h2o2w1(i) * (fso2w1(i) + nso2w1(i)) +  &
     &           h2o2w2(i) * (fso2w2(i) + nso2w2(i))) * dtaq
              zso2_o3(iw) = zso2_o3(iw) +  &
     &          a(iw, 7) * 0.5 *  &
     &          (fso2w1(i) + nso2w1(i) + fso2w2(i) + nso2w2(i)) * dtaq

          enddo

      end if

!     *-----------------------------------------------------------------
!     * Compute a final species concentration based on the mean of
!       the two derivatives.
!     *-----------------------------------------------------------------

        do i = 1, numwk

!.... H2O2
          h2o2w1(i) = max( h2o2w1(i)  &
     &              + 0.5 * ( dh2o21(i) + dh2o22(i) ) * dtaq, 0.0d0)

!.... fossil
          fso2w1(i) = max( fso2w1(i)  &
     &              + 0.5 * ( dfso21(i) + dfso22(i) ) * dtaq, 0.0d0)

!.... natural
          ndmsw1(i) = max( ndmsw1(i)  &
     &              + 0.5 * ( dndms1(i) + dndms2(i) ) * dtaq, 0.0d0)
          nso2w1(i) = max( nso2w1(i)  &
     &              + 0.5 * ( dnso21(i) + dnso22(i) ) * dtaq, 0.0d0)

!.... biomass
!         bso2w1(i) = max( bso2w1(i)
!    &              + 0.5 * ( dbso21(i) + dbso22(i) ) * dtaq, 0.0d0)

        enddo


!     *--------------------
!     * End time step loop.
!     *--------------------

  200 continue


!     ******************************************************************
!     * Perform time step analysis.
!     ******************************************************************

      if( ndt .eq. ndtmin ) then


!     *-----------------------------------------------------------------
!     * If ndt = ndtmin halve the time step.
!     * No grid-point solution can exist at this point.
!     *-----------------------------------------------------------------

        ndt = 2 * ndtmin

        do i = 1, numwk

!.... H2O2
          h2o2o(i) = h2o2w1(i)

!.... fossil
          fso2o(i) = fso2w1(i)

!.... natural
          ndmso(i) = ndmsw1(i)
          nso2o(i) = nso2w1(i)

!.... biomass
!         bso2o(i) = bso2w1(i)

        enddo

!       =========
        go to 100
!       =========

      else


!     *-----------------------------------------------------------------
!     * If ndt > ndtmin check the solution.
!     * If the solution satisfies the accuracy condition store it.
!     * Note, aqueous H2SO4, ()so4a, has no clear sky contribution.
!     *-----------------------------------------------------------------

!.... compute the fractional changes from the previous cycle

        do i = 1, numwk

          rh2o2(i) = ( h2o2o(i)  - h2o2w1(i) )  &
     &             / ( h2o2w1(i) + emin )

          rfso2(i) = ( fso2o(i)  - fso2w1(i) )  &
     &             / ( fso2w1(i) + emin )

!         rndms(i) = ( ndmso(i)  - ndmsw1(i) )
!    &             / ( ndmsw1(i) + emin )

          rnso2(i) = ( nso2o(i)  - nso2w1(i) )  &
     &             / ( nso2w1(i) + emin )

!         rbso2(i) = ( bso2o(i)  - bso2w1(i) )
!    &             / ( bso2w1(i) + emin )

        enddo


!.... if the solution passes the accuray test store the solution
!     from the working slot "i" back into grid-point idxwk(i)

!     The change in H2SO4 will be computed using the DMS and SO2
!     sources and the changes in DMS and SO2.

!     The changes in H2O2, SO2 and DMS due to clear sky chemistry have
!     already been computed, and weighted by ( 1 - cfp ), in chem11,
!     so here it remains to add their aqueous changes, weighted by
!     cfp, to these clear air changes.

        do i = 1, numwk

          if( abs( rh2o2(i) ) .lt. acc .and.  &
     &        abs( rfso2(i) ) .lt. acc .and.  &
     &        abs( rnso2(i) ) .lt. acc       ) then
!    &        abs( rndms(i) ) .lt. acc .and.

            ndone   = ndone + 1
            nfin(i) = -1

            iw = idxwk(i)

!.... H2O2
            dh2o2( iw )  = dh2o2( iw )  &
     &                   + cfp( iw ) * ( h2o2w1(i) - h2o2( iw ) )

!.... fossil
            dfso2( iw )  = dfso2( iw )  &
     &                   + cfp( iw ) * ( fso2w1(i) - fso2( iw ) )

            if (do_aerocom) then
            dfso4a( iw ) = cfp( iw )  &
     &                   * ( sfso2(iw) * dtclc + sfso4(iw) * dtclc    & ! dtchem
     &                   - ( fso2w1(i) - fso2( iw ) ) )

            else
            dfso4a( iw ) = cfp( iw )  &
     &                   * ( sfso2(iw) * dtclc     & ! dtchem
     &                   - ( fso2w1(i) - fso2( iw ) ) )
            end if

!.... natural
            dnso2( iw )  = dnso2( iw )  &
     &                   + cfp( iw ) * ( nso2w1(i) - nso2( iw ) )

            dndms( iw )  = dndms( iw )  &
     &                   + cfp( iw ) * ( ndmsw1(i) - ndms( iw ) )

            dnso4a( iw ) = cfp( iw )  &
     &                   * ( ( snso2(iw) + sndms(iw) ) * dtclc     & ! dtchem
     &                     - ( nso2w1(i) - nso2( iw ) )  &
     &                     - ( ndmsw1(i) - ndms( iw ) ) )

!.... biomass
!           dbso2( iw )  = dbso2( iw )
!     &                  + cfp( iw ) * ( bso2w1(i) - bso2( iw ) )
!           dbso4a( iw ) = cfp( iw )
!     &                  * (   sbso2(iw) * dtclc   ! dtchem
!     &                    - ( bso2w1(i) - bso2( iw ) ) )

          endif

        enddo

      endif  ! if( ndt .eq. ndtmin ) then


!     *****************************************************************
!     * If ndone = numcl all grid points have a solution, return.
!     * If ndone < numcl halve the time step of those that have
!       not converged and try again.
!     *****************************************************************


!.... if all grid points have been solved, return

      if(ndone .eq. numcl) then
!       =========
        go to 787
!       =========
      endif

!.... some grid points are yet to be solved,
!     compress out solved grid points and store the latest
!     computed values in h2o2o, etc.

      ndt = 2 * ndt

      num1 = 0

      do i = 1, numwk

        if( nfin(i) .eq. 1 ) then

          num1 = num1 + 1

          nfin( num1 )  = 1

          idxwk( num1 ) = idxwk(i)

          h2o2o( num1 ) = h2o2w1(i)
          fso2o( num1 ) = fso2w1(i)
          ndmso( num1 ) = ndmsw1(i)
          nso2o( num1 ) = nso2w1(i)
!         bso2o( num1 ) = bso2w1(i)

        endif

      enddo

      numwk = num1

!     =========
      go to 100
!     =========

 787  continue


!     conv: convert SO4 concentration from # cm-3 to gm gm-1
      conv = wtmair / wtmso4

!.... H2O2

        dcp(:, IH2O2)  = dh2o2(:)

!.... fossil

        dcp(:, IFSO2)  = dfso2(:)
        dcp(:, IFSO4A) = dfso4a(:)  &
     &                     / (conv * zmp(:))   ! unit: g/g

!.... natural

        dcp(:, INSO2)  = dnso2(:)
        dcp(:, INDMS)  = dndms(:)
        dcp(:, INSO4A) = dnso4a(:)  &
     &                     / (conv * zmp(:))   ! unit: g/g

        if(pr_sulf_src) then

          dms_oh(:) = dms_oh(:)  &
     &                + cfp(:) * zdms_oh(:) * massc(:)  &
     &                / zmp(:) * wtmso2 / wtmair
          dms_no3(:) = dms_no3(:)  &
     &                + cfp(:) * zdms_no3(:) * massc(:)  &
     &                / zmp(:) * wtmso2 / wtmair
          so2_oh(:) = so2_oh(:)  &
     &                + cfp(:) * zso2_oh(:) * massc(:)  &
     &                / zmp(:) * wtmso4 / wtmair
          so2_h2o2(:) = so2_h2o2(:)  &
     &                + cfp(:) * zso2_h2o2(:) * massc(:)  &
     &                / zmp(:) * wtmso4 / wtmair
          so2_o3(:) = so2_o3(:)  &
     &                + cfp(:) * zso2_o3(:) * massc(:)  &
     &                / zmp(:) * wtmso4 / wtmair

        end if


      return

      end

