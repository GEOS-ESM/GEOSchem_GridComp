      SUBROUTINE sediment(levels,dt,cellDepth,m,kel,hno3c,h2oc,
     &  	          rmedsts,rmednat,rmedice,denssts,vfall,gcSC)

! Condensed nitric acid and water have units of mole fraction.

      USE SC_GridCompMod
      IMPLICIT NONE

      INTEGER, PARAMETER :: DBL = KIND(0.00D+00)
 
      TYPE(SC_GridComp), INTENT(INOUT) :: gcSC  ! SC Grid Component

      INTEGER, INTENT(IN) :: levels

      REAL, INTENT(IN) :: dt,cellDepth(levels)
      REAL, INTENT(IN) :: m(levels),kel(levels)
      REAL, INTENT(IN) :: denssts(levels),rmednat(levels)
      REAL, INTENT(IN) :: rmedice(levels),rmedsts(levels)

      REAL, INTENT(INOUT) :: hno3c(levels),h2oc(levels)

      REAL, INTENT(OUT) :: vfall(levels)

      LOGICAL :: go
      INTEGER :: maxht, minht
      INTEGER :: k, idtsed, itime
      REAL :: dtsed, fluxcorr, vfallmax, dens, fallvel, dzmin
      REAL(KIND=DBL) :: prod, loss, cna, cw
      REAL, PARAMETER :: MIN_VALUE = 1.00E-20

!    create 3d pressure array press & 3d array of vertical gridbox heights dz:
          do k=1,levels
            vfall(k) = 0.
          enddo

!    at each lat and lon, first decide if sedimentation needs
!    to be calculated.  This is done by checking to see if hno3cond
!    or h2ocond is nonzero somewhere in the vertical column.

!    Sedimentation of H2O ice is done in AGCM.  Test only for HNO3.
          go = .false.
          do k = 1,levels
            IF(hno3c(k) > MIN_VALUE) go=.TRUE.
          enddo

          if(go) then

!    determine fall velocity field and max fall velocity,
!    and upper alt of condensed phase.  Here, only one fall
!    velocity field needs to be calculated, because we never
!    have a coexisting population of Type 1 and Type 2 particles.
            maxht = 1
            minht = levels
            vfallmax = 0.0

            do k=1,levels
              IF(hno3c(k) > MIN_VALUE .OR. h2oc(k) > MIN_VALUE) THEN
                if(k.gt.maxht) maxht = k
                if(k.lt.minht) minht = k
               endif
             enddo

            if(maxht.eq.levels) maxht = levels-1
            if(minht.eq.1) minht = 2

!    initialize dzmin with an arbitrarily large number
            dzmin = 1.e7
            do k=minht,maxht
!    find the minimum gridbox height (depth is in cm)
              dzmin = min(dzmin,cellDepth(k))

!    If there is condensed phase h2o at the gridpoint, then
!    code assumes all condensed hno3 is contained in h2o particles.
!    Therefore, sediment hno3 using h2o particle size distribution.
              IF(h2oc(k) > MIN_VALUE) THEN
                fluxcorr = exp(8.0*LOG(gcSC%sigice)*LOG(gcSC%sigice))
!    calculate fall velocity using ice particle size distribution when there is
!    any condensed h2o. If constant particle number option is chosen, use
!    the particle median radius calculated in saice (rmedice).
                if(gcSC%constantnice) then
                   vfall(k) = fluxcorr*fallvel(rmedice(k)*1.e4
     &               ,gcSC%densice,m(k),kel(k))
                 else
                   vfall(k) = fluxcorr*
     &               fallvel(gcSC%rice,gcSC%densice,m(k),kel(k))
                 endif
               else

                 if(gcSC%calcsts)then
!    particles are sts - use sts size distribution width
                   fluxcorr = exp(8.0*LOG(gcSC%sigsts)*LOG(gcSC%sigsts))
                   if(gcSC%constantnsts)then
                     vfall(k) = fluxcorr*
     &                 fallvel(rmedsts(k)*1.e4,denssts(k),
     &                 m(k),kel(k))
                    else
                      vfall(k) = fluxcorr*fallvel(gcSC%rsts,
     &                  denssts(k),m(k),kel(k))
                    endif

                  else   
!    particles are NAT - use nat size distribution width
                    fluxcorr = exp(8.0*LOG(gcSC%signat)*LOG(gcSC%signat))
                    if(gcSC%constantnnat)then
                       vfall(k) = fluxcorr*
     &                   fallvel(rmednat(k)*1.e4,gcSC%densnat,
     &                   m(k),kel(k))
                     else
                       vfall(k) = fluxcorr*
     &                   fallvel(gcSC%rnat,gcSC%densnat,m(k),kel(k))
                     endif
                  endif
               endif

!    vfallmax is recalculated at each longitude
              vfallmax = max(vfall(k),vfallmax)
             enddo
! JEN ======================================================================
            IF(vfallmax <= 0.00) THEN
	     PRINT *,"StratChem::SEDIMENT: *WARNING* Found vanishing fall speed."
	    ELSE
! JEN ======================================================================
!    calculate time step for sedimentation.  The time step is calculated
!    so that no particle can fall through the gridbox below it in one
!    sedimentation time step.
            dtsed = dzmin/vfallmax
            if(dt/dtsed.le.1.0)then
              idtsed = 1
              dtsed = dt             ! no need to change timestep
             else
               idtsed = int(dt/dtsed)+1 ! number of iterations necessary
               idtsed = min(idtsed,10)  ! limit number of iterations
               dtsed = dt/(idtsed*1.0)
               vfallmax = dzmin/dtsed 
             endif

!    Do sedimentation. Loop idtsedal times, use time step of dtsedal:
            do itime=1,idtsed
              do k=1,maxht

!    calculate production and loss terms for sedimentation of hno3
!    these terms should be calculated using the fall velocity interpolated
!    to the gridbox edges.  I am not doing this simply for ease of coding. 
                prod = dtsed/cellDepth(k+1)
     &            *min(vfall(k+1),vfallmax)
     &            *hno3c(k+1)*m(k+1)

                loss = dtsed/cellDepth(k)
     &            *min(vfall(k),vfallmax)
     &            *hno3c(k)*m(k)

!    actually sediment the hno3.  The code uses a very simple explicit
!    scheme that also has the advantage of conserving column hno3cond.
!    The commented out implicit update does not conserve. Note that the
!    formula used is not a strict discretization of a differential 
!    equation.  This is because vfall changes in each gridbox 
!    depending on whether the PSC type in the gridbox is a type 1
!    or a type 2 PSC.  A discretization of a differential equation
!    would involve vertical gradients of the vall velocity.  But
!    what does that mean when one of the gridboxes has type 1
!    pscs in it, and the other has type 2 (resulting in a large
!    gradient in vfall)?  If I want to be more "correct" in future
!    incarnations, I need to define fall velocity fields for type 1
!    and type 2 PSCs, etc.
		cna = hno3c(k)
                cna = cna +(prod-loss)/m(k)
		hno3c(k) = cna

!    This shouldn't be needed because we have an explicit update scheme
!    and we are changing the time step, but do it anyway.
                IF(hno3c(k) < MIN_VALUE) hno3c(k) = 0.0
                IF(h2oc(k) > MIN_VALUE .OR. h2oc(k+1) > MIN_VALUE) THEN

!    sediment h2o.  Reversal of prod and loss is legacy from
!    version using "dehyd"
                  prod = dtsed/cellDepth(k)*min(vfall(k),vfallmax)
     &              *h2oc(k)
                  loss = dtsed/cellDepth(k+1)*min(vfall(k+1),vfallmax)
     &              *h2oc(k+1)*m(k+1)

!    production of dehydration is a loss for h2ocond.
                  cw = h2oc(k)
                  cw = cw +loss/m(k)-prod
                  h2oc(k) = cw

!    This shouldn't be needed because we have an explicit update scheme
!    and we are changing the time step, but do it anyway.
                  h2oc(k) = max(h2oc(k),0.0)

                 endif
               enddo      ! vertical loop
             enddo        ! time loop
! JEN ======================================================================
	    END IF
! JEN ======================================================================
           endif          ! matches if(go)
 
      RETURN
      END SUBROUTINE sediment
