!=============================================================================
!
! $Id$
!
! CODE DEVELOPER
!   originated from GRANTOUR programmed by John Walton
!   modified by Xiaohong Liu
!
! FILE
!   const_so2topar.F
!
! ROUTINES
!   Do_Const_So2toPar
!
!=============================================================================

!-----------------------------------------------------------------------------
!
! ROUTINE
!   Do_Const_So2toPar (chem60)
!
! DESCRIPTION
!
!     *     nocyc     :   no. of sub-cycles of "dtchem"
!     *     dtcyc     :   cycle time step
!
!     *     frdust    :   the fraction of dust molecules
!                         available to dissolve SO4
!     *     frsslt    :   the fraction of sea salt molecules
!                         available to dissolve SO4
!
!     *     fmd       :   fossil SO4 on dust 1 - 4           [#/cm3]
!           fms       :   fossil SO4 on sslt 1 - 4              "
!
!     *     nmd       :   natural SO4 on dust 1 - 4          [#/cm3]
!           nms       :   natural SO4 on sslt 1 - 4             "
!
!     *     rtcd      :   rate coef. for SO2 -> dust transfer  [s-1]
!     *     rtcs      :   rate coef. for SO2 -> sslt transfer    "
!
!     *     rtctot    :   this is "s" in the notes below       [s-1]
!
!     *     dust      :   dust concentration [1 - 4]         [#/cm3]
!     *     sslt      :   sslt concentration [1 - 4]         [#/cm3]
!
!       The limits on the number of SO4 molecules that can go
!       into solution on dust and sea salt particles.  They
!       will be a prescribed fraction, "frdust" and "frsslt",
!       of the species concentrations.  These are the "M" in
!       the notes below.  Once these have been computed the
!       dust and sea salt concentrations will not be used.
!
!     *     ldso4 [1 - 4]                                    [#/cm3]
!     *     lsso4 [1 - 4]                                    [#/cm3]
!
!.... the calculated changes in concentration of SO2 vapor
!     and the dissolved SO4, "dZ" & "dm" in the notes
!
!     *     dfz, dnz, dbz                                    [#/cm3]
!
!     *     dfmd [1 - 4]                                     [#/cm3]
!           dfms [1 - 4]                                        "
!
!     *     dnmd [1 - 4]                                     [#/cm3]
!           dnms [1 - 4]                                        "
!
!.... the sum of the 3 categories of sulfate
!     on dust and sea salt
!
!     *     smd [1 - 4]                                      [#/cm3]
!     *     sms [1 - 4]                                      [#/cm3]
!
!.... total sulfate change
!
!     *     dzeta                                            [#/cm3]
!
!.... total SO4 change on dust and sea salt
!     when dissolved SO4 increasing
!
!     *     dmud [1 - 4]                                     [#/cm3]
!     *     dmus [1 - 4]                                        "
!
!     ********** notes:
!
!  * The solution for the transfer of SO2 vapor to the surface
!    of a particle and its subsequent conversion to SO4.
!
!  n       sulfur source index
!          1 - fossil, 2 - natural, 3 - biomass
!
!  j       particle type index
!          1 - 4  for dust 1 - 4,
!          5 - 8  for sea salt 1 - 4
!
!  mj,n    [SO4] dissolved on a particle
!          (dust or sea salt)                               [#/cm3]
!
!  mtj     tot., over source type, of
!          [SO4] dissolved on particle type j
!
!          = mj,1 + mj,2 + mj,3                                "
!
!  Mj      max [SO4] dissolved on particle type j           [#/cm3]
!
!  Zn      [SO2] vapor, from source type n                  [#/cm3]
!
!  Z       total [SO2] vapor
!
!          = Z1 + Z2 + Z3                                      "
!
!  dt      local time step                                      [s]
!
!  kj      transfer rate coefficients                         [s-1]
!
!  approximations
!  --------------
!
!      (1)  The variation of mtj is small relative to Mj
!
!      (2)  Z is allowed to vary while integrating m
!
!  [SO2] equation, source type n
!  -----------------------------
!
!  (1)    dZn/dt = - [ k1 ( 1 - mt1 / M1 ) + k2 ( 1 - mt2 / M2 )
!                    + k3 ( 1 - mt3 / M3 ) + k4 ( 1 - mt4 / M4 ) ] Zn
!
!                = - dm1,n/dt - dm2,n/dt - dm3,n/dt - dm4,n/dt
!
!                = - s Zn
!  where
!
!  (2)     s = k1 ( 1 - mt1 / M1 ) + k2 ( 1 - mt2 / M2 )
!            + k3 ( 1 - mt3 / M3 ) + k4 ( 1 - mt4 / M4 )
!
!  By the first approximation I can write
!
!  (3)     Zn(t) = Zn(0) exp( - s t )
!
!  and
!
!  (4)     Z(t) = Z(0) exp( - s t ).
!
!  The change in Zn, or Z, over dt is
!
!          dZn = Zn(t+dt) - Zn(t)
!
!              = Zn(0) ( exp[ - s (t+dt) ] - exp[ - s t ] )
!
!  (5)     dZn = Zn(t) ( exp[ - s dt ] - 1 )
!
!  and
!
!  (6)     dZ = Z(t) ( exp[ - s dt ] - 1 )
!
!  j = 1-4 and j = 5-8 - [SO4] equation (dust and sea salt)
!  ---------------------------------------------------------
!
!  (7)     dmj,n/dt = kj ( 1 - mtj / Mj ) Zn
!
!                   = kj ( 1 - mtj / Mj ) * Zn(0) exp( - s t )
!
!  Summing over n in (7) and using the definition of mtj
!
!  (8)     dmtj/dt = kj ( 1 - mtj / Mj ) Z
!
!                  = kj ( 1 - mtj / Mj ) Z(0) exp( - s t )
!
!  Remember,
!
!          Z * dt = - dZ / s
!  so
!          dmtj = - (kj/s) * ( 1 - mtj / Mj ) * dZ
!
!  Integrating (8) over t -> t+dt gives
!
!  (9)     mtj(t+dt) = Mj - [ Mj - mtj(t) ] exp[ (kj/s) * (dZ/Mj) ].
!
!  The change in mtj over dt is then
!
!  (10)    dmtj = [ Mj - mtj ] ( 1 - exp[ (kj/s) * (dZ/Mj) ] ).
!
!  Finally, dmtj is apportioned among the 3 groups of species,
!  fossil, natural and biomass, by
!
!  (11)    dmj,n = dmtj ( dZn / dZ ).
!
!  Given the dmj,n, dZn can be re-evaluated
!
!  (12)    dZn = - dm1,n - dm2,n - dm3,n - dm4,n
!
!  which gurantees number concentration.
!
!     *----------------------------------------------------------------
!
!     This is a rewriting of the change calculation.
!
!  M   (ldso4 and lsso4) limit on the total dissolved SO4 on
!      dust and sea salt particles
!
!  f   (fms) the current amount of dissolved fossil SO4
!  n   (nms) the current amount of dissolved natural SO4
!  b   (bms) the current amount of dissolved biomass SO4
!
!  mt  = f + n + b, (smd) total dissolved SO4 on dust
!  mt  = f + n + b, (sms) total dissolved SO4 on sea salt
!
!  df  (dfms) change in dissolved fossil SO4
!  dn  (dnms) change in dissolved natural SO4
!  db  (dbms) change in dissolved biomass SO4
!
!  dm  (dmus) change in total dissolved SO4
!
!  q   fractional change in dissolved SO4 when there is a loss
!
!  e   the exponential term in (10)
!
!  if M > mt
!  ---------
!
!  The total dissolved SO4 will approach the limit M, (10)
!  and will be partitioned by (11).
!
!  if M = mt
!  ---------
!
!  No change.
!
!  if M < mt
!  ---------
!
!  No change.   (as per Michael Herzog)
!
!     *----------------------------------------------------------------
!       Dust will now be treated as sea salt is, so limits etc.
!       will have to be calculated.
!
!     *----------------------------------------------------------------
!
!-----------------------------------------------------------------------------

      subroutine Do_Const_So2toPar  &
     &  (dtchem, itloop, zmp,  &
     &   rtcd, rtcs, cp)

      implicit none

#     include "setkin_par.h"
#     include "sulfchem.h"

!     ----------------------
!     GMIMOD-related values.
!     ----------------------

      integer :: itloop

      real*8  :: dtchem  ! model time step in seconds

      real*8  :: zmp(itloop)

      real*8  :: rtcd(itloop, 4)
      real*8  :: rtcs(itloop, 4)

      real*8  :: cp(itloop, NSP)

! -------------------------------------------------------------
!  *** LOCAL  VARIABLES WERE DECLARED HERE.
! -------------------------------------------------------------

      real*8 :: frso2(itloop)

      real*8 :: dust(itloop, 4)
      real*8 :: sslt(itloop, 4)

      real*8 :: ldso4(itloop, 4)
      real*8 :: lsso4(itloop, 4)

      real*8 :: smd(itloop, 4)
      real*8 :: sms(itloop, 4)

      real*8 :: rtctot(itloop)
      real*8 :: dfz   (itloop)
      real*8 :: dnz   (itloop)
      real*8 :: dzeta (itloop)

      real*8 :: dmud(itloop, 4)
      real*8 :: dmus(itloop, 4)
      real*8 :: dfmd(itloop, 4)
      real*8 :: dfms(itloop, 4)
      real*8 :: dnmd(itloop, 4)
      real*8 :: dnms(itloop, 4)

      dust = 0.0d0
      sslt = 0.0d0

      ldso4 = 0.0d0
      lsso4 = 0.0d0

      smd = 0.0d0
      sms = 0.0d0

      rtctot = 0.0d0
      dfz    = 0.0d0
      dnz    = 0.0d0
      dzeta  = 0.0d0

      dmud = 0.0d0
      dmus = 0.0d0
      dfmd = 0.0d0
      dfms = 0.0d0
      dnmd = 0.0d0
      dnms = 0.0d0


!.... the "dtchem" sub-time step

      dtcyc = dtchem / float( nocyc )

!     *****************************************************************
!     * Calculate the concentrations of dust and sea salt, and using
!       these, obtain for each particle size range the maximum
!       concentration of SO4 that can be reached.
!     * Note, the units conversion is to molecules / cm3 NOT
!       particles / cm3.
!     * The concentration limits depend only upon the size distribution
!       of the particles and its mass mixing ratios.  It does not
!       depend on the source from which the sulfur comes.
!     * The sea salt concentrations are not needed after the limits
!       have been computed.
!     *****************************************************************

!.... dust and sea salt concentrations [molecules/cm3]

!     The mass of Ca is 0.05 * the mass of dust.
!     The molecular weight of Ca is 40.08.

      conv     = wtmair / wtmso4
      convdust = wtmair / wtmdust
      convsslt = wtmair / wtmsslt

      do ijk = 1, itloop

        dust(ijk,1) = cp(ijk, IDUST1)
        dust(ijk,2) = cp(ijk, IDUST2)
        dust(ijk,3) = cp(ijk, IDUST3)
        dust(ijk,4) = cp(ijk, IDUST4)

        sslt(ijk,1) = cp(ijk, ISSLT1)
        sslt(ijk,2) = cp(ijk, ISSLT2)
        sslt(ijk,3) = cp(ijk, ISSLT3)
        sslt(ijk,4) = cp(ijk, ISSLT4)

      enddo

      do i = 1, 4
        dust(:,i) = dust(:,i) * convdust * zmp(:)
        sslt(:,i) = sslt(:,i) * convsslt * zmp(:)
      enddo


!.... SO4 concentration limits [molecules/cm3]

      do i = 1, 4
        ldso4(:,i) = frdust * dust(:,i)
        lsso4(:,i) = frsslt * sslt(:,i)
      enddo

!.... neglect limits < 0.1 molecule/cm3

      do i = 1, 4
        where ( ldso4(:,i) .lt. 0.1 ) ldso4(:,i) = 0.0d0
        where ( lsso4(:,i) .lt. 0.1 ) lsso4(:,i) = 0.0d0
      enddo


!     *****************************************************************
!     * Advance the SO2 vapor and SO4 dissolved on particles.
!     *****************************************************************

      do 100 nc = 1, nocyc       ! time step, dtcyc, loop begin

!.... the sum of dissolved sulfate on the five size ranges
!     of dust and sea salt

!     These are the mtj

        do ijk = 1, itloop

          smd(ijk,1) = cp(ijk, IFSO4D1) + cp(ijk, INSO4D1)
          smd(ijk,2) = cp(ijk, IFSO4D2) + cp(ijk, INSO4D2)
          smd(ijk,3) = cp(ijk, IFSO4D3) + cp(ijk, INSO4D3)
          smd(ijk,4) = cp(ijk, IFSO4D4) + cp(ijk, INSO4D4)

          sms(ijk,1) = cp(ijk, IFSO4S1) + cp(ijk, INSO4S1)
          sms(ijk,2) = cp(ijk, IFSO4S2) + cp(ijk, INSO4S2)
          sms(ijk,3) = cp(ijk, IFSO4S3) + cp(ijk, INSO4S3)
          sms(ijk,4) = cp(ijk, IFSO4S4) + cp(ijk, INSO4S4)

        enddo

        do i = 1, 4
          smd(:,i) = smd(:,i) * conv * zmp(:)
          sms(:,i) = sms(:,i) * conv * zmp(:)
        enddo


!.... neglect totals < 0.1 molecule/cm3

        do i = 1, 4
          where ( smd(:,i) .lt. 0.1 ) smd(:,i) = 0.0d0
          where ( sms(:,i) .lt. 0.1 ) sms(:,i) = 0.0d0
        enddo

!     rtctot is the [SO2] rate coefficient "s" in (2).
!     add 1.d-20 to eliminate 0 divides

        rtctot(:) = 1.d-20

!.... dust part of the rate coefficient, "s" in (2)

        do i = 1, 4

          where ( ldso4(:,i) .gt. 1.d-20 .and.  &
     &        ldso4(:,i) .gt. smd(:,i) )  &
     &      rtctot(:) = rtctot(:)  &
     &        + rtcd(:,i) * ( 1. - smd(:,i) / ldso4(:,i) )

        enddo

!.... add in sea salt particles

        do i = 1, 4

          where ( lsso4(:,i) .gt. 1.d-20 .and.  &
     &        lsso4(:,i) .gt. sms(:,i) )  &
     &      rtctot(:) = rtctot(:)  &
     &        + rtcs(:,i) * ( 1. - sms(:,i) / lsso4(:,i) )

        enddo


!     *----------------------------------------------------------------
!     * Change in SO2, "dZn" and "dZ" in (5) and (6)
!       the 1.d-20 in "dzeta" prevents a zero divide.
!         dbz(ijk) = bso2(ijk) * ( exp( - rtctot(ijk) * dtcyc ) - 1. )
!         dzeta(ijk) = dfz(ijk) + dnz(ijk) + dbz(ijk) - 1.d-20
!     *----------------------------------------------------------------

        do ijk = 1, itloop

          dfz(ijk) = cp(ijk, IFSO2) *  &
     &               ( exp( - rtctot(ijk) * dtcyc ) - 1.0d0 )
          dnz(ijk) = cp(ijk, INSO2) *  &
     &               ( exp( - rtctot(ijk) * dtcyc ) - 1.0d0 )
          dzeta(ijk) = dfz(ijk) + dnz(ijk) - 1.d-20

        enddo


!     *----------------------------------------------------------------
!     * Change in SO4 on dust, "dmj,n", see (10) and (11)
!     *----------------------------------------------------------------

!.... total change

          do i = 1, 4

            where ( ldso4(:,i) .gt. 1.d-20 .and.  &
     &          ldso4(:,i) .gt. smd(:,i) )
              dmud(:,i) = ( ldso4(:,i) - smd(:,i) )  &
     &                    * ( 1.0d0 - exp( ( rtcd(:,i) / rtctot(:) )  &
     &                                * dzeta(:) / ldso4(:,i) ) )
            elsewhere
              dmud(:,i) = 0.
            end where

          enddo

!.... apportion total change

          do i = 1, 4
            dfmd(:,i) = dmud(:,i) * dfz(:) / dzeta(:)
            dnmd(:,i) = dmud(:,i) * dnz(:) / dzeta(:)
!           dbmd(:,i) = dmud(:,i) * dbz(:) / dzeta(:)
          enddo


!     *----------------------------------------------------------------
!     * Change in SO4 on sea salt, "dmj,n", see (10) and (11).
!     *----------------------------------------------------------------

!.... total change

          do i = 1, 4

            where ( lsso4(:,i) .gt. 1.d-20 .and.  &
     &          lsso4(:,i) .gt. sms(:,i) )
              dmus(:,i) = ( lsso4(:,i) - sms(:,i) )  &
     &                    * ( 1.0d0 - exp( ( rtcs(:,i) / rtctot(:) )  &
     &                                * dzeta(:) / lsso4(:,i) ) )
            elsewhere
              dmus(:,i) = 0.
            end where

          enddo

!.... apportion total change

          do i = 1, 4
            dfms(:,i) = dmus(:,i) * dfz(:) / dzeta(:)
            dnms(:,i) = dmus(:,i) * dnz(:) / dzeta(:)
!           dbms(:,i) = dmus(:,i) * dbz(:) / dzeta(:)
          enddo


!     *----------------------------------------------------------------
!     * Update dso2 using the changes in dissolved SO4 concentrations
!       just computed, (12) this guarantees # conservation.
!     *----------------------------------------------------------------

        dfz(:) = 0.0d0
        dnz(:) = 0.0d0
!       dbz(:) = 0.0d0

        do i = 1, 4
          dfz(:) = dfz(:) - dfmd(:,i) - dfms(:,i)
          dnz(:) = dnz(:) - dnmd(:,i) - dnms(:,i)
!         dbz(:) = dbz(:) - dbmd(:,i) - dbms(:,i)
        enddo


!     The following loops guarantee that SO2 will not go negative
!     while at the same time conserving #.

!.... fossil

          do i = 1, 4

            where ( dfz(:) .lt. - cp(:, IFSO2) )

              frso2    = - cp(:, IFSO2) / dfz(:)
              dfz(:) = - cp(:, IFSO2)

              dfmd(:,i) = frso2 * dfmd(:,i)
              dfms(:,i) = frso2 * dfms(:,i)

            end where

          end do

!.... natural

          do i = 1, 4

            where ( dnz(:) .lt. - cp(:, INSO2) )

              frso2    = - cp(:, INSO2) / dnz(:)
              dnz(:) = - cp(:, INSO2)

              dnmd(:,i) = frso2 * dnmd(:,i)
              dnms(:,i) = frso2 * dnms(:,i)

            end where

          end do


!     *----------------------------------------------------------------
!     * update the SO2 vapor and dissolved SO4
!     *----------------------------------------------------------------

!.... debug begin
!           make sure there is no source of SO2

        do ijk = 1, itloop
          if( (dfz(ijk) .gt. 0.0d0) .or.  &
     &        (dnz(ijk) .gt. 0.0d0) ) then
            write(*  ,10)
            stop
          endif
        enddo
   10   format( / 'chem60 error exit, there is a source of SO2' /)

!.... debug end


        do i = 1, 4
          dfmd(:,i) = dfmd(:,i) / (conv * zmp(:))
          dfms(:,i) = dfms(:,i) / (conv * zmp(:))
          dnmd(:,i) = dnmd(:,i) / (conv * zmp(:))
          dnms(:,i) = dnms(:,i) / (conv * zmp(:))
        enddo


        cp(:, IFSO2) = cp(:, IFSO2) + dfz(:)
        cp(:, INSO2) = cp(:, INSO2) + dnz(:)

        cp(:, IFSO4D1) = cp(:, IFSO4D1) + dfmd(:,1)
        cp(:, IFSO4D2) = cp(:, IFSO4D2) + dfmd(:,2)
        cp(:, IFSO4D3) = cp(:, IFSO4D3) + dfmd(:,3)
        cp(:, IFSO4D4) = cp(:, IFSO4D4) + dfmd(:,4)

        cp(:, IFSO4S1) = cp(:, IFSO4S1) + dfms(:,1)
        cp(:, IFSO4S2) = cp(:, IFSO4S2) + dfms(:,2)
        cp(:, IFSO4S3) = cp(:, IFSO4S3) + dfms(:,3)
        cp(:, IFSO4S4) = cp(:, IFSO4S4) + dfms(:,4)

        cp(:, INSO4D1) = cp(:, INSO4D1) + dnmd(:,1)
        cp(:, INSO4D2) = cp(:, INSO4D2) + dnmd(:,2)
        cp(:, INSO4D3) = cp(:, INSO4D3) + dnmd(:,3)
        cp(:, INSO4D4) = cp(:, INSO4D4) + dnmd(:,4)

        cp(:, INSO4S1) = cp(:, INSO4S1) + dnms(:,1)
        cp(:, INSO4S2) = cp(:, INSO4S2) + dnms(:,2)
        cp(:, INSO4S3) = cp(:, INSO4S3) + dnms(:,3)
        cp(:, INSO4S4) = cp(:, INSO4S4) + dnms(:,4)


  100 continue

      write (6,*) " I am here"

      return

      end

