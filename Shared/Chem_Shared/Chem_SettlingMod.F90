#include "unused_dummy.H"
!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  Chem_SettlingMod --- Gravitional Sedimentation & Settling Speed
!
! !INTERFACE:
!

   module  Chem_SettlingMod

! !USES:

   use Chem_Mod
   use Chem_ConstMod, only: grav        ! Constants !
   use Chem_UtilMod

   use m_mpout

   implicit none

! !PUBLIC TYPES:
!
   PRIVATE

!
! !PUBLIC MEMBER FUNCTIONS:
!

   PUBLIC  Chem_Settling
   PUBLIC  Chem_SettlingSimple
   PUBLIC  Chem_CalcVsettle

!
! !DESCRIPTION:
!
!  This module implements the (pre-ESMF) DU Grid Component. 
!
! !REVISION HISTORY:
!
!  16Sep2003 da Silva  First crack.
!
!EOP
!-------------------------------------------------------------------------
CONTAINS
!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Chem_Settling - 
!
! !INTERFACE:
!

   subroutine Chem_Settling ( i1, i2, j1, j2, km, nbeg, nend, nbins, flag, &
                              radiusInp, rhopInp, cdt, w_c, tmpu, rhoa, &
                              hsurf, hghte, fluxout, rc, &
                              vsettleOut, correctionMaring )

! !USES:

  implicit NONE

! !INPUT PARAMETERS:

   integer, intent(in) :: i1, i2, j1, j2, km, nbeg, nend, nbins
   integer, intent(in) :: flag     ! flag to control particle swelling (see note)
   real, intent(in)    :: cdt
   real, pointer, dimension(:)     :: radiusInp, rhopInp
   real, pointer, dimension(:,:)   :: hsurf
   real, pointer, dimension(:,:,:) :: tmpu, rhoa, hghte

! !OUTPUT PARAMETERS:

   type(Chem_Bundle), intent(inout) :: w_c        ! Chemical tracer fields
   type(Chem_Array), pointer        :: fluxout(:) ! Mass lost by settling
                                                  ! to surface, kg/m2/s
   integer, intent(out)             :: rc         ! Error return code:
                                                  !  0 - all is well
                                                  !  1 - 
!  Optionally output the settling velocity calculated
   type(Chem_Array), pointer, optional, dimension(:)  :: vsettleOut

!  Optionally correct the settling velocity following Maring et al, 2003
   logical, optional, intent(in)    :: correctionMaring

   character(len=*), parameter :: myname = 'Settling'

! !DESCRIPTION: Gravitational settling of aerosol between vertical
!               layers.  Assumes input radius in [m] and density (rhop) 
!               in [kg m-3]. If flag is set, use the Fitzgerald 1975 (flag = 1)
!               or Gerber 1985 (flag = 2) parameterization to update the 
!               particle radius for the calculation (local variables radius
!               and rhop).
!
! !REVISION HISTORY:
!
!  15May2019  Darmenov  Refactor and speed up code
!  17Sep2004  Colarco   Strengthen sedimentation flux out at surface
!                       by setting removal to be valid from middle of
!                       surface layer
!  06Nov2003  Colarco   Based on Ginoux
!  23Jan2003  da Silva  Standardization
!
!EOP
!-------------------------------------------------------------------------

! !Local Variables
   integer  ::  i, j, k, n, dk
   
   real, parameter ::  rhow = 1000.  ! Density of water [kg m-3]

   real :: vsettle(i1:i2,j1:j2,km)               ! fall speed [m s-1]
   real :: dz(i1:i2,j1:j2,km)                    ! layer thickness [m]

   real, pointer, dimension(:,:,:) :: qa

   real*8 :: cmass_before(i1:i2,j1:j2)
   real*8 :: cmass_after(i1:i2,j1:j2)

   real, dimension(i1:i2,j1:j2,km) :: radius
   real, dimension(i1:i2,j1:j2,km) :: rhop

!  parameters from Maring et al, 2003
   real, parameter :: v_upwardMaring = 0.33e-2   ! upward velocity, [m s-1]
   real, parameter :: diameterMaring = 7.30e-6   ! particle diameter, [m]

   real, parameter :: one_over_g = 1.0/grav

   _UNUSED_DUMMY(nend)
   rc = 0


!  Handle the fact that hghte may be in the range [1,km+1] or [0,km]
!  -----------------------------------------------------------------
   dk = lbound(hghte,3) - 1  ! This is either 0 or 1
  
!  Layer thickness from hydrostatic equation
   k = km
   dz(:,:,k) = hghte(:,:,k+dk)-hsurf(:,:)
   do k = km-1, 1, -1
    dz(:,:,k) = hghte(:,:,k+dk) - hghte(:,:,k+dk+1)
   enddo


!  Loop over the number of dust bins
   do n = 1, nbins

    ! alias
    qa => w_c%qa(nbeg+n-1)%data3d


    if( associated(fluxout(n)%data2d) ) fluxout(n)%data2d(i1:i2,j1:j2) = 0.0

    cmass_before(:,:) = 0.d0
    cmass_after(:,:)  = 0.d0

!   If radius le 0 then get out of loop
    if(radiusInp(n) .le. 0.) cycle

!   Find the column dry mass before sedimentation
    do k = 1, km
     do j = j1, j2
      do i = i1, i2
       cmass_before(i,j) = cmass_before(i,j) + qa(i,j,k) * w_c%delp(i,j,k) * one_over_g
      enddo
     enddo
    enddo

!   Particle swelling
    call ParticleSwelling(i1, i2, j1, j2, km, w_c%rh, radiusInp(n), rhopInp(n), radius, rhop, flag)

!   Settling velocity of the wet particle
    do k = 1, km
     do j = j1, j2
      do i = i1, i2
       call Chem_CalcVsettle(radius(i,j,k), rhop(i,j,k), rhoa(i,j,k), &
                             tmpu(i,j,k), vsettle(i,j,k))
      end do
     end do
    end do


    if(present(correctionMaring)) then
     if (correctionMaring) then
       vsettle = max(1.0e-9, vsettle - v_upwardMaring)
     endif
    endif

    if(present(vsettleOut)) then
     vsettleOut(n)%data3d = vsettle
    endif

!   Time integration
    call SettlingSolver(i1, i2, j1, j2, km, cdt, w_c%delp, dz, vsettle, qa)

!   Find the column dry mass after sedimentation and thus the loss flux
    do k = 1, km
     do j = j1, j2
      do i = i1, i2
       cmass_after(i,j) = cmass_after(i,j) + qa(i,j,k) * w_c%delp(i,j,k) * one_over_g
      enddo
     enddo
    enddo

    if( associated(fluxout(n)%data2d) ) then
       fluxout(n)%data2d(i1:i2,j1:j2) &
        = (cmass_before(i1:i2,j1:j2) - cmass_after(i1:i2,j1:j2))/cdt
    endif

   end do   ! n

 end  subroutine Chem_Settling


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Chem_SettlingSimple - support for single bin settling call
!
! !INTERFACE:
!

   subroutine Chem_SettlingSimple ( i1, i2, j1, j2, km, ibin, flag, &
                              radiusInp, rhopInp, cdt, w_c, tmpu, rhoa, &
                              hsurf, hghte, fluxout, rc, &
                              vsettleOut, correctionMaring )

! !USES:

  implicit NONE

! !INPUT PARAMETERS:

   integer, intent(in) :: i1, i2, j1, j2, km, ibin
   integer, intent(in) :: flag     ! flag to control particle swelling (see note)
   real, intent(in)    :: cdt, radiusInp, rhopInp
   real, pointer, dimension(:,:)   :: hsurf
   real, pointer, dimension(:,:,:) :: tmpu, rhoa, hghte

! !OUTPUT PARAMETERS:

   type(Chem_Bundle), intent(inout) :: w_c        ! Chemical tracer fields
   type(Chem_Array), pointer        :: fluxout    ! Mass lost by settling
                                                  ! to surface, kg/m2/s
   integer, intent(out)             :: rc         ! Error return code:
                                                  !  0 - all is well
                                                  !  1 - 
!  Optionally output the settling velocity calculated
   type(Chem_Array), pointer, optional  :: vsettleOut

!  Optionally correct the settling velocity following Maring et al, 2003
   logical, optional, intent(in)    :: correctionMaring

   character(len=*), parameter :: myname = 'Settling'

! !DESCRIPTION: Gravitational settling of aerosol between vertical
!               layers.  Assumes input radius in [m] and density (rhop) 
!               in [kg m-3]. If flag is set, use the Fitzgerald 1975 (flag = 1)
!               or Gerber 1985 (flag = 2) parameterization to update the 
!               particle radius for the calculation (local variables radius
!               and rhop).
!
! !REVISION HISTORY:
!
!  17Sep2004  Colarco   Strengthen sedimentation flux out at surface
!                       by setting removal to be valid from middle of
!                       surface layer
!  06Nov2003  Colarco   Based on Ginoux
!  23Jan2003  da Silva  Standardization
!
!EOP
!-------------------------------------------------------------------------

! !Local Variables
   integer  ::  i, j, k, dk

   real, parameter ::  rhow = 1000.  ! Density of water [kg m-3]

   real :: vsettle(i1:i2,j1:j2,km)               ! fall speed [m s-1]
   real :: dz(i1:i2,j1:j2,km)                    ! layer thickness [m]

   real, pointer, dimension(:,:,:) :: qa

   real*8 :: cmass_before(i1:i2,j1:j2)
   real*8 :: cmass_after(i1:i2,j1:j2)

   real, dimension(i1:i2,j1:j2,km) :: radius
   real, dimension(i1:i2,j1:j2,km) :: rhop

!  parameters from Maring et al, 2003
   real, parameter :: v_upwardMaring = 0.33e-2   ! upward velocity, [m s-1]
   real, parameter :: diameterMaring = 7.30e-6   ! particle diameter, [m]


   real, parameter :: one_over_g = 1.0/grav

   rc = 0


!  Handle the fact that hghte may be in the range [1,km+1] or [0,km]
!  -----------------------------------------------------------------
   dk = lbound(hghte,3) - 1  ! This is either 0 or 1
  
!  Layer thickness from hydrostatic equation
   k = km
   dz(:,:,k) = hghte(:,:,k+dk)-hsurf(:,:)
   do k = km-1, 1, -1
    dz(:,:,k) = hghte(:,:,k+dk) - hghte(:,:,k+dk+1)
   enddo

   ! alias
   qa => w_c%qa(ibin)%data3d

   if( associated(fluxout%data2d) ) fluxout%data2d(i1:i2,j1:j2) = 0.0

   cmass_before(:,:) = 0.d0
   cmass_after(:,:)  = 0.d0

!  If radius le 0 then get out of loop
   if(radiusInp .le. 0.) return

!  Find the column dry mass before sedimentation
   do k = 1, km
    do j = j1, j2
     do i = i1, i2
      cmass_before(i,j) = cmass_before(i,j) + qa(i,j,k) * w_c%delp(i,j,k)*one_over_g
     enddo
    enddo
   enddo

!  Particle swelling
   call ParticleSwelling(i1, i2, j1, j2, km, w_c%rh, radiusInp, rhopInp, radius, rhop, flag)

!  Settling velocity of the wet particle
   do k = 1, km
    do j = j1, j2
     do i = i1, i2
      call Chem_CalcVsettle(radius(i,j,k), rhop(i,j,k), rhoa(i,j,k), &
                            tmpu(i,j,k), vsettle(i,j,k))
     end do
    end do
   end do


   if(present(correctionMaring)) then
    if (correctionMaring) then
      vsettle = max(1.0e-9, vsettle - v_upwardMaring)
    endif
   endif

   if(present(vsettleOut)) then
    if(associated(vsettleOut%data3d)) vsettleOut%data3d = vsettle
   endif

!  Time integration
   call SettlingSolver(i1, i2, j1, j2, km, cdt, w_c%delp, dz, vsettle, qa)

!  Find the column dry mass after sedimentation and thus the loss flux
   do k = 1, km
    do j = j1, j2
     do i = i1, i2
      cmass_after(i,j) = cmass_after(i,j) + qa(i,j,k) * w_c%delp(i,j,k) * one_over_g
     enddo
    enddo
   enddo

   if( associated(fluxout%data2d) ) then
      fluxout%data2d(i1:i2,j1:j2) &
       = (cmass_before(i1:i2,j1:j2) - cmass_after(i1:i2,j1:j2))/cdt
   endif

 end subroutine Chem_SettlingSimple


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Chem_CalcVsettle - Calculate the aerosol settling velocity
!
! !INTERFACE:
!

   subroutine Chem_CalcVsettle ( radius, rhop, rhoa, tmpu, &
                                 vsettle )

! !USES:

   implicit NONE

! !INPUT PARAMETERS:

   real, intent(in)    :: radius              ! Particle radius [m]
   real, intent(in)    :: rhop                ! Particle density [kg m-3]
   real, intent(in)    :: rhoa                ! Layer air density [kg m-3]
   real, intent(in)    :: tmpu                ! Layer temperature [K]

! !OUTPUT PARAMETERS:

   real, intent(out)   :: vsettle                 ! Layer fall speed [m s-1]

   character(len=*), parameter :: myname = 'Vsettle'

! !DESCRIPTION: Calculates the aerosol settling velocity and Brownian diffusion
!               coefficient
!               Follows discussions in Seinfeld and Pandis, Pruppacher and
!               Klett, and the coding in CARMA (Toon et al., 1988)
!               Should work satisfactorily for al reasonable sized aerosols
!               (up to Reynolds number 300)
!
! !REVISION HISTORY:
!
!  06Nov2003  Colarco   Initial version.
!  23Jan2003  da Silva  Standardization
!
!EOP
!-------------------------------------------------------------------------

! !Local Variables
   real :: rmu                       ! Dynamic viscosity [kg m-1 s-1]
   real :: vt                        ! Thermal velocity of air molecule [m s-1]
   real :: rmfp                      ! Air molecule mean free path [m]
   real :: bpm                       ! Cunningham slip correction factor
   real :: rkn                       ! Knudsen number
   real :: re, x, y                  ! reynold's number and parameters
   real, parameter :: kb = 1.3807e-23 ! Boltzmann constant [kg m2 s-1 K-1 mol-1]
   real, parameter :: m_air = 4.8096e-26 ! Mass of <avg> air molecule [kg]
   real, parameter :: pi = 3.141529265
   
   real, parameter :: f_vt = 8*kb/pi/m_air
   real, parameter :: f_diff_coef = kb/(6*pi)
   real, parameter :: two_over_nine = 2./9.

   real, parameter :: a0 = -3.18657
   real, parameter :: a1 =  0.992696
   real, parameter :: a2 = -1.53193e-3
   real, parameter :: a3 = -9.870593e-4
   real, parameter :: a4 = -5.78878e-4
   real, parameter :: a5 =  8.55176e-5 
   real, parameter :: a6 = -3.27815e-6


!  Dynamic viscosity from corrected Sutherland's Equation
   rmu = 1.8325e-5*(416.16/(tmpu+120.))*(tmpu/296.16)**1.5

!  Thermal velocity of air molecule
   vt = sqrt(tmpu * f_vt)

!  Air molecule mean free path
   rmfp = 2*rmu/(rhoa*vt)

!  Knudsen number
   rkn = rmfp/radius

!  Cunningham slip correction factor
!  bpm = 1. + 1.246*rkn + 0.42*rkn*exp(-0.87/rkn)
!  linearized form, Binkowski and Shankar (equation A27, 1995) 
   bpm = 1 + 1.246*rkn

!  Brownian diffusion coefficient
!  diff_coef = tmpu*bpm/(rmu*radius) * f_diff_coef

!  Fall speed (assumes Reynolds # < 0.01)
   vsettle = two_over_nine*rhop*radius*radius*grav*bpm/rmu

!  Check the Reynold's number to see if we need a drag correction
!  First guess at Reynold's number using Stoke's calculation
   re = 2.*rhoa*radius*vsettle/rmu

!  If Re > 0.01 then apply drag correction following Pruppacher and
!  Klett regime 2 (eq. 10-142).  Assuming reasonable aerosols we
!  do not consider that particle Re may exceed 300.
   if(re .gt. 0.01) then
    x  = log(24.*re/bpm)
    y  = a0 + x*(a1 + x*(a2 + x*(a3 + x*(a4 + x*(a5 + a6*x)))))
    re = exp(y)*bpm
    vsettle = 0.5*rmu*re/(rhoa*radius)
   endif

   end subroutine Chem_CalcVsettle

   
   subroutine SettlingSolver(i1, i2, j1, j2, km, cdt, delp, dz, vs, qa)

    implicit none

    integer, intent(in) :: i1, i2
    integer, intent(in) :: j1, j2
    integer, intent(in) :: km

    real,    intent(in) :: cdt

    real, dimension(i1:i2,j1:j2,km), intent(in) :: delp
    real, dimension(i1:i2,j1:j2,km), intent(in) :: dz
    real, dimension(i1:i2,j1:j2,km), intent(in) :: vs

    real, dimension(i1:i2,j1:j2,km), intent(inout) :: qa


    ! local
    integer :: i, j, iit
    integer :: nSubSteps

    real, dimension(i1:i2, j1:j2, km) :: tau

    real, dimension(km) :: dp_
    real, dimension(km) :: tau_

    real :: dt, dt_cfl
    

    tau = vs/dz

    do j = j1, j2
      do i = i1, i2

          dp_  = delp(i,j,:)
          tau_ = tau(i,j,:)
          
          dt_cfl  = 1 / maxval(tau_)

          if (dt_cfl > cdt) then
              ! no need for time sub-splitting
              nSubSteps = 1
              dt = cdt
          else
              nSubSteps = ceiling(cdt / dt_cfl)
              dt = cdt/nSubSteps
          end if

          do iit = 1, nSubSteps
              qa(i,j,   1) = qa(i,j,   1) * (1 - dt*tau_(1))
              qa(i,j,2:km) = qa(i,j,2:km) + ( (dp_(1:km-1)/dp_(2:km))*(dt*tau_(1:km-1))*qa(i,j,1:km-1) ) &
                                          - dt*tau_(2:km)*qa(i,j,2:km)
          end do

      enddo
    enddo

   end subroutine SettlingSolver



   subroutine ParticleSwelling(i1, i2, j1, j2, km, rh, radius_dry, rhop_dry, radius, rhop, flag)

    implicit none

    integer, intent(in) :: i1, i2
    integer, intent(in) :: j1, j2
    integer, intent(in) :: km

    real, dimension(i1:i2,j1:j2,km), intent(in)  :: rh

    integer, intent(in) :: flag

    real, intent(in) :: radius_dry
    real, intent(in) :: rhop_dry

    real, dimension(i1:i2,j1:j2,km), intent(out) :: radius  ! radius  of the wet particle
    real, dimension(i1:i2,j1:j2,km), intent(out) :: rhop    ! density of the wet particle


    select case (flag)
        case (0)
            radius = radius_dry
            rhop   = rhop_dry

        case (1)
            call ParticleSwelling_Fitzgerald(i1, i2, j1, j2, km, rh, radius_dry, rhop_dry, radius, rhop)

        case (2)
            call ParticleSwelling_Gerber(i1, i2, j1, j2, km, rh, radius_dry, rhop_dry, radius, rhop)

        case (3)
            call ParticleSwelling_Gerber_AmmoniumSulfate(i1, i2, j1, j2, km, rh, radius_dry, rhop_dry, radius, rhop)  

        case (4)
            call ParticleSwelling_PK2007(i1, i2, j1, j2, km, rh, radius_dry, rhop_dry, radius, rhop)

        case default
            radius = radius_dry
            rhop   = rhop_dry
    end select 

   end subroutine ParticleSwelling



   subroutine ParticleSwelling_Fitzgerald(i1, i2, j1, j2, km, rh, radius_dry, rhop_dry, radius, rhop)

    implicit none

    integer, intent(in) :: i1, i2
    integer, intent(in) :: j1, j2
    integer, intent(in) :: km

    real, dimension(i1:i2,j1:j2,km), intent(in)  :: rh

    real, intent(in) :: radius_dry 
    real, intent(in) :: rhop_dry

    real, dimension(i1:i2,j1:j2,km), intent(out) :: radius  ! radius  of the wet particle
    real, dimension(i1:i2,j1:j2,km), intent(out) :: rhop    ! density of the wet particle


    ! local
    real, parameter :: rhow = 1000.  ! density of water [kg m-3]

    ! the following parameters relate to the swelling of seasalt like particles
    ! following Fitzgerald, Journal of Applied Meteorology, 1975.
    real, parameter :: epsilon = 1.   ! soluble fraction of deliqeuscing particle
    real, parameter :: alphaNaCl = 1.35

    real :: alpha, alpha1, alpharat, beta, theta, f1, f2
    real :: sat, rrat

    integer :: i, j, k


!   Adjust particle size for relative humidity effects,
!   based on Fitzgerald, Journal of Applied Meteorology, 1975

    do k = 1, km
     do j = j1, j2
      do i = i1, i2

       radius(i,j,k) = radius_dry
       rhop(i,j,k) = rhop_dry

       sat = rh(i,j,k)

       if (sat > 0.80) then
!       parameterization blows up for RH > 0.995, so set that as max
!       rh needs to be scaled 0 - 1
        sat = min(0.995, sat)

!       Calculate the alpha and beta parameters for the wet particle
!       relative to amonium sulfate
        beta = exp( (0.00077*sat) / (1.009-sat) )

        if(sat .le. 0.97) then
         theta = 1.058
        else
         theta = 1.058 - (0.0155*(sat-0.97)) /(1.02-sat**1.4)
        endif

        alpha1 = 1.2*exp( (0.066*sat) / (theta-sat) )

! no need of this calculations, because epsilon == 1
!       f1 = 10.2 - 23.7*sat + 14.5*sat*sat
!       f2 = -6.7 + 15.5*sat -  9.2*sat*sat
!       alpharat = 1. - f1*(1.-epsilon) - f2*(1.-epsilon**2.)
!       alpha = alphaNaCl * (alpha1*alpharat)
! instead, it is faster to do 
        alpha = alphaNaCl * alpha1

        radius(i,j,k) = alpha * radius_dry**beta

        rrat = radius_dry/radius(i,j,k)
        rrat = rrat*rrat*rrat

        rhop(i,j,k) = rrat*rhop_dry + (1 - rrat)*rhow
       endif

     end do
    end do
   end do

   end subroutine ParticleSwelling_Fitzgerald


   subroutine ParticleSwelling_Gerber(i1, i2, j1, j2, km, rh, radius_dry, rhop_dry, radius, rhop)

    implicit none

    integer, intent(in) :: i1, i2
    integer, intent(in) :: j1, j2
    integer, intent(in) :: km

    real, dimension(i1:i2,j1:j2,km), intent(in)  :: rh

    real, intent(in)  :: radius_dry 
    real, intent(in)  :: rhop_dry

    real, dimension(i1:i2,j1:j2,km), intent(out) :: radius  ! radius  of the wet particle
    real, dimension(i1:i2,j1:j2,km), intent(out) :: rhop    ! density of the wet particle


    ! local
    real, parameter :: rhow = 1000.  ! density of water [kg m-3]

    ! parameters from Gerber 1985 (units require radius in cm, see the variable rcm)
    real, parameter :: c1=0.7674, c2=3.079, c3=2.573e-11, c4=-1.424

    real :: sat, rrat, rcm

    integer :: i, j, k


    ! Adjust the particle size for relative humidity effects,
    ! based on Gerber 1985
    do k = 1, km
     do j = j1, j2
      do i = i1, i2
 
       sat = max(rh(i,j,k), tiny(1.0)) ! to avoid zero FPE
       sat = min(0.995, sat)

       rcm = radius_dry*100. ! radius in 'cm'

       radius(i,j,k) = 0.01 * ( c1*rcm**c2 / (c3*rcm**c4-alog10(sat)) &
                                + rcm*rcm*rcm )**(1./3.)

       rrat = radius_dry/radius(i,j,k)
       rrat = rrat*rrat*rrat

       rhop(i,j,k) = rrat*rhop_dry + (1 - rrat)*rhow
 
      end do
     end do
    end do

   end subroutine ParticleSwelling_Gerber


   subroutine ParticleSwelling_Gerber_AmmoniumSulfate(i1, i2, j1, j2, km, rh, radius_dry, rhop_dry, radius, rhop)

    implicit none

    integer, intent(in) :: i1, i2
    integer, intent(in) :: j1, j2
    integer, intent(in) :: km

    real, dimension(i1:i2,j1:j2,km), intent(in)  :: rh

    real, intent(in) :: radius_dry 
    real, intent(in) :: rhop_dry

    real, dimension(i1:i2,j1:j2,km), intent(out) :: radius  ! radius  of the wet particle
    real, dimension(i1:i2,j1:j2,km), intent(out) :: rhop    ! density of the wet particle


    ! local
    real, parameter :: rhow = 1000.  ! density of water [kg m-3]

    ! parameters for ammonium sulfate from Gerber 1985 (units require radius in cm, see the variable rcm)
    real, parameter :: SU_c1=0.4809, SU_c2=3.082, SU_c3=3.110e-11, SU_c4=-1.428

    real :: sat, rrat, rcm

    integer :: i, j, k


    ! Adjust the particle size for relative humidity effects,
    ! based on Gerber parameterization for Ammonium Sulfate

    do k = 1, km
     do j = j1, j2
      do i = i1, i2
 
       sat = max(rh(i,j,k), tiny(1.0)) ! to avoid zero FPE
       sat = min(0.995, sat)

       rcm = radius_dry*100. ! radius in 'cm'
       radius(i,j,k) = 0.01 * ( SU_c1*rcm**SU_c2 / (SU_c3*rcm**SU_c4-alog10(sat)) &
                                + rcm*rcm*rcm )**(1./3.)

       rrat = radius_dry/radius(i,j,k)
       rrat = rrat*rrat*rrat

       rhop(i,j,k) = rrat*rhop_dry + (1 - rrat)*rhow
 
      end do
     end do
    end do

   end subroutine ParticleSwelling_Gerber_AmmoniumSulfate


   subroutine ParticleSwelling_PK2007(i1, i2, j1, j2, km, rh, radius_dry, rhop_dry, radius, rhop)

    implicit none

    integer, intent(in) :: i1, i2
    integer, intent(in) :: j1, j2
    integer, intent(in) :: km

    real, dimension(i1:i2,j1:j2,km), intent(in)  :: rh

    real, intent(in) :: radius_dry 
    real, intent(in) :: rhop_dry

    real, dimension(i1:i2,j1:j2,km), intent(out) :: radius  ! radius  of the wet particle
    real, dimension(i1:i2,j1:j2,km), intent(out) :: rhop    ! density of the wet particle


    ! local
    real, parameter :: rhow = 1000.  ! density of water [kg m-3]

    real :: sat, rrat

    integer :: i, j, k

 
    ! Adjust the particle size for relative humidity effects,
    ! based on Petters and Kreidenweis (ACP2007) parameterization

    do k = 1, km
     do j = j1, j2
      do i = i1, i2

       sat = rh(i,j,k)
       sat = min(0.99, sat)

       radius(i,j,k) = radius_dry * (1+1.19*sat/(1-sat))**(1./3.)

       rrat = radius_dry/radius(i,j,k)
       rrat = rrat*rrat*rrat

       rhop(i,j,k) = rrat*rhop_dry + (1 - rrat)*rhow

     end do
    end do
   end do

   end subroutine ParticleSwelling_PK2007


   end module Chem_SettlingMod
