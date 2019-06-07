#include "MAPL_Exceptions.h"
#include "MAPL_Generic.h"


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: MAML_WetRemovalMod - Gravitational sedimentation/settling and dry
!                               deposition of aerosol particles and gases
!
! !INTERFACE:
!
   module MAML_WetRemovalMod
!
! !USES:
!
   use MAPL_Mod

   use MAPL_ConstantsMod, only: MAPL_GRAV

   use WetRemovalMod, only: WetRemovalGOCART


   implicit None
   private

!
! !PUBLIC MEMBER FUNCTIONS:

   public MAML_WetRemoval

!
! !PUBLIC PARAMETERS:

!
! !PRIVATE PARAMETERS:

   real, private, parameter :: g_E = MAPL_GRAV     ! standard gravity,        'm s-2'

!
! !DESCRIPTION: 
!
!  {\tt MAML\_SettlingMod} provides a collection of methods for 
!  modeling graviational sedimentation/settling of aerosol particles
!
!
! !REVISION HISTORY:
!
!  29Nov2011  A. Darmenov  Initial version
!
!EOP
!-------------------------------------------------------------------------


   contains


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: MAML_WetRemovalAerosol --- 
!
! !INTERFACE:

   subroutine MAML_WetRemoval(qa, f_wet, aero_type, kin, ple, tmpu, rhoa, &
                              pfllsan, pfilsan, precc, precl, cdt, flux, rc)

! !USES:

   implicit None


! !INPUT/OUTPUT PARAMETERS:
   real, pointer, dimension(:,:,:), intent(inout) :: qa       ! mixing ratio,     'kg kg-1' or '# kg-1'
   real, pointer, dimension(:,:),   intent(inout) :: flux     ! wet removal flux, 'kg m-2 s-1' or '# m-2 s-1'

! !INPUT PARAMETERS:
   real, intent(in)                            :: f_wet       ! large scale wet removal efficiency
   character(len=*)                            :: aero_type   ! aerosol type
   logical, intent(in)                         :: kin         ! true for aerosol

   real, pointer, dimension(:,:,:), intent(in) :: ple         ! pressure at model layer edges, 'Pa'
   real, pointer, dimension(:,:,:), intent(in) :: tmpu        ! temperature
   real, pointer, dimension(:,:,:), intent(in) :: rhoa
   real, pointer, dimension(:,:,:), intent(in) :: pfllsan, pfilsan
   real, pointer, dimension(:,:),   intent(in) :: precc
   real, pointer, dimension(:,:),   intent(in) :: precl
   real, intent(in)                            :: cdt         ! time step, 's'

! !OUTPUT PARAMETERS:
   integer, intent(out)                        :: rc          ! return code


! !DESCRIPTION: Based on WetRemovalGOCART() subroutine.
!
! !REVISION HISTORY:
!
!  10Dec2012  A. Darmenov
!
!EOP
!-------------------------------------------------------------------------
   
                  __Iam__('MAML_WetRemoval')


   integer             :: i1, i2, j1, j2, k1, km    ! indexes
   integer             :: i, j, k, kk, n            ! loop counters
   
   integer             :: LH
 
   real                :: F, B, BT                  ! temporary variables on cloud, freq.
   real                :: qmx, qd, A                ! temporary variables on moisture

   real                :: pls, pcv, pac             ! ls, cv, tot precip [mm day-1]
   real, dimension(:), allocatable :: qls, qcv      ! ls, cv portion dqcond [kg m-3 s-1]

   real, dimension(:,:,:), allocatable :: pdog      ! air mass factor dp/g [kg m-2]
   real, dimension(:,:,:), allocatable :: delz      ! box height  dp/g/rhoa [m]

   real                :: WASHFRAC, WASHFRAC_F_14 
   real, allocatable   :: fd(:,:)                   ! flux across layers [kg m-2]
   real, allocatable   :: dpfli(:,:,:)              ! vertical gradient of LS ice+rain precip flux 
   real, allocatable   :: DC(:)                     ! scavenge change in mass mixing ratio
   real, allocatable   :: c_h2o(:,:,:), cldliq(:,:,:), cldice(:,:,:)


   real                :: effRemoval

   integer, parameter  :: nbins = 1                  ! work with a single bin tracer

   real, parameter     :: grav = g_E 

   !  Rain parameters from Liu et al.
   real, parameter     :: B0_ls = 1.0e-4
   real, parameter     :: F0_ls = 1.0
   real, parameter     :: XL_ls = 5.0e-4
   real, parameter     :: B0_cv = 1.5e-3
   real, parameter     :: F0_cv = 0.3
   real, parameter     :: XL_cv = 2.0e-3
   real, parameter     :: k_wash = 1.d0  ! first order washout rate, constant, [cm^-1]

!  Duration of rain: ls = model timestep, cv = 1800 s (<= cdt)
   real                :: Td_ls
   real, parameter     :: Td_cv = 1800.0
   real(kind=8), parameter :: R = 8.2057d-2  ! universal gas constant [L*atm/moles/K]
   real(kind=8), parameter :: INV_T0 = 1d0 / 298d0
   real(kind=8), parameter :: conv_NH3 = 5.69209978831d-1 ! 0.6*SQRT(0.9) for ice to gas ratio
   real(kind=8)        :: k_rain, Kstar298, H298_R, I2G, L2G, C_TOT, F_L, F_I
   real(kind=8)        :: PP, LP

   logical :: snow_scavenging

   i1 = lbound(qa, 1); i2 = ubound(qa, 1)
   j1 = lbound(qa, 2); j2 = ubound(qa, 2)
   k1 = lbound(qa, 3); km = ubound(qa, 3)
   
   effRemoval = f_wet

   rc = 0

!  Initialize local variables
!  --------------------------

!  c_h2o, cldliq, and cldice are respectively intended to be the 
!  water mixing ratio (liquid or vapor?, in or out of cloud?)
!  cloud liquid water mixing ratio
!  cloud ice water mixing ratio

   allocate( c_h2o(i1:i2,j1:j2,km), __STAT__)
   allocate(cldliq(i1:i2,j1:j2,km), __STAT__)
   allocate(cldice(i1:i2,j1:j2,km), __STAT__)


   c_h2o  = (10d0**(-2663.5d0/tmpu(:,:,:) + 12.537d0 ) ) /  &
                   (ple(:,:,0:km-1)+ple(:,:,1:km)) /2d0   
   cldliq = 0.d0
   where(tmpu >  248.) cldliq = 1.d-6 * ( ( tmpu - 248.d0) / 20.d0 )
   where(tmpu >= 268.) cldliq = 1.d-6
   cldice = 1.d-6 - cldliq

   Td_ls = cdt

   if (associated(flux)) flux = 0.0

!  Allocate the dynamic arrays
   allocate(fd(km,nbins), __STAT__)
   allocate(dc(nbins),    __STAT__)

   allocate(qls(km),      __STAT__)         
   allocate(qcv(km),      __STAT__)
   
   allocate(dpfli(i1:i2,j1:j2,km), __STAT__)
   allocate(pdog(i1:i2,j1:j2,km), __STAT__)
   allocate(delz(i1:i2,j1:j2,km), __STAT__)


!  Accumulate the 3-dimensional arrays of rhoa and pdog
   pdog = (ple(:,:,1:km)-ple(:,:,0:km-1)) / grav
   delz = pdog / rhoa
   
   dpfli = pfllsan(:,:,1:km)-pfllsan(:,:,0:km-1)+pfilsan(:,:,1:km)-pfilsan(:,:,0:km-1)

   if (.not. KIN) then              ! Gases
    if (aero_type == 'NH3') then  ! Only for NH3 at present
    ! values adopted in Umich/IMPACT and GMI, effective Henry's law coefficient at pH=5
      Kstar298 = 1.05d6
      H298_R = -4.2d3
    else
      if (MAPL_AM_I_ROOT()) print *, 'stop in WetRemoval, need Kstar298 and H298_R'
      stop
    endif
   endif

!  Snow scavenging flag
   snow_scavenging = .true.

   if (aero_type == 'BC' .and. n == 2) then
       snow_scavenging = .false.
   end if

   if ( (aero_type == 'OC'      ) .or. &
        (aero_type == 'sea_salt') .or. &
        (aero_type == 'sulfur'  ) .or. &
        (aero_type == 'seasalt' ) .or. &
        (aero_type == 'sulfate' ) .or. &
        (aero_type == 'NH3'     ) .or. &
        (aero_type == 'NH4a'    ) .or. &
        (aero_type == 'nitrate' ) .or. &
        (aero_type == 'dust'    ) ) then
      snow_scavenging = .false.
   end if

!  Loop over spatial indices
   do j = j1, j2
    do i = i1, i2

!    Check for total precipitation amount
!    Assume no precip in column if precl+precc = 0
     pac = precl(i,j) + precc(i,j)
     if(pac .le. 0.) goto 100
     pls = precl(i,j)
     pcv = precc(i,j)

!    Initialize the precipitation fields
     qls(:)  = 0.
     qcv(:)  = 0.
     fd(:,:) = 0.

!    Find the highest model layer experiencing rainout.  Assumes no
!    scavenging if T < 258 K
     LH = 0
     do k = 1, km
      if(dpfli(i,j,k) .gt. 0. ) then
       LH = k
       goto 15
      endif
     end do
 15  continue
     if(LH .lt. 1) goto 100

     do k = LH, km
      qls(k) = dpfli(i,j,k)/pdog(i,j,k)*rhoa(i,j,k)
     end do

!    Loop over vertical to do the scavenging!
     do k = LH, km

!-----------------------------------------------------------------------------
!   (1) LARGE-SCALE RAINOUT:             
!       Tracer loss by rainout = TC0 * F * exp(-B*dt)
!         where B = precipitation frequency,
!               F = fraction of grid box covered by precipitating clouds.
!       We assume that tracer scavenged by rain is falling down to the
!       next level, where a fraction could be re-evaporated to gas phase
!       if Qls is less then 0 in that level.
!-----------------------------------------------------------------------------
      if (qls(k) .gt. 0.) then
       F  = F0_ls / (1. + F0_ls*B0_ls*XL_ls/(qls(k)*cdt/Td_ls))
       k_rain  = B0_ls/F0_ls +1./(F0_ls*XL_ls/qls(k)) 
       if ( kin ) then     ! Aerosols
          B = k_rain
       else                ! Gases
        ! ice to gas ratio
          if ( c_h2o(i,j,k) > 0.d0) then
             I2G = (cldice(i,j,k) / c_h2o(i,j,k)) * conv_NH3
          else
             I2G = 0.d0
          endif
          L2G = cldliq(i,j,k) * R * tmpu(i,j,k) * &
                  Kstar298 * EXP( -H298_R * ( ( 1d0 / tmpu(i,j,k) ) - INV_T0 ) )
        ! fraction of NH3 in liquid & ice phases
          C_TOT = 1d0 + L2G + I2G
          F_L = L2G / C_TOT
          F_I = I2G / C_TOT
        ! compute kg, the retention factor for liquid NH3 is 0 at T < 248K and
        ! 0.05 at 248K < T < 268K
          if (tmpu(i,j,k) >=268d0) then
             B = k_rain * ( F_L+F_I )
          elseif ( (248d0 < tmpu(i,j,k)) .and. (tmpu(i,j,k) < 268d0) ) then
             B = k_rain * ( (0.05*F_L)+F_I )
          else
             B = k_rain * F_I
          endif
       endif ! kin
       BT = B * Td_ls
       if (BT.gt.10.) BT = 10.               !< Avoid overflow >
!      Adjust du level:
       do n = 1, nbins
! implementing QQ Wang's change for snow scavening by Huisheng Bian (4/24/2015)
! supress scavenging at cold T except BC1 (hydrophobic), dust, and hno3
        if (tmpu(i,j,k) < 258d0 .and. .not.snow_scavenging) then
           F = 0.d0
        endif

        effRemoval = f_wet
        DC(n) = qa(i,j,k) * F * effRemoval *(1.-exp(-BT))
        if (DC(n).lt.0.) DC(n) = 0.
        qa(i,j,k) = qa(i,j,k)-DC(n)
        if (qa(i,j,k) .lt. 1.0E-32) qa(i,j,k) = 1.0E-32
       end do
!      Flux down:  unit is kg m-2
!      Formulated in terms of production in the layer.  In the revaporation step
!      we consider possibly adding flux from above...
       do n = 1, nbins
        Fd(k,n) = DC(n)*pdog(i,j,k)
       end do

      end if                                    ! if Qls > 0  >>>

!-----------------------------------------------------------------------------
! * (2) LARGE-SCALE WASHOUT:
! *     Occurs when rain at this level is less than above.
!-----------------------------------------------------------------------------
      if(k .gt. LH .and. qls(k) .ge. 0.) then
       if(qls(k) .lt. qls(k-1)) then
!       Find a maximum F overhead until the level where Qls<0.
        Qmx   = 0.
        do kk = k-1,LH,-1
         if (Qls(kk).gt.0.) then
          Qmx = max(Qmx,Qls(kk))
         else
          goto 333
         end if
        end do

 333    continue
        F = F0_ls / (1. + F0_ls*B0_ls*XL_ls/(Qmx*cdt/Td_ls))
        if (F.lt.0.01) F = 0.01
!-----------------------------------------------------------------------------
!  The following is to convert Q(k) from kgH2O/m3/sec to mm/sec in order
!  to use the Harvard formula.  Convert back to mixing ratio by multiplying
!  by rhoa.  Multiply by pdog gives kg/m2/s of precip.  Divide by density
!  of water (=1000 kg/m3) gives m/s of precip and multiply by 1000 gives
!  units of mm/s (omit the multiply and divide by 1000).
!-----------------------------------------------------------------------------

!       Aerosols
        Qd = Qmx /rhoa(i,j,k)*pdog(i,j,k)
        if (Qd.ge.50.) then
         B = 0.
        else
         B = Qd * 0.1
        end if
        BT = B * cdt
        if (BT.gt.10.) BT = 10.

!       Gases
        if ( .not. KIN ) then
           IF ( tmpu(i,j,k) >= 268d0 ) THEN
            !------------------------
            ! T >= 268K: Do washout
            !------------------------
            ! Rainwater content in the grid box (Eq. 17, Jacob et al, 2000)
            PP = (PFLLSAN(i,j,k)/1000d0 + PFILSAN(i,j,k)/917d0 )*100d0 ! from kg H2O/m2/s to cm3 H2O/cm2/s
            LP = ( PP * cdt ) / ( F * delz(i,j,k)*100.d0 )  ! DZ*100.d0 in cm
            ! Compute liquid to gas ratio for H2O2, using the appropriate 
            ! parameters for Henry's law -- also use rainwater content Lp
            ! (Eqs. 7, 8, and Table 1, Jacob et al, 2000)
            !CALL COMPUTE_L2G( Kstar298, H298_R, tmpu(i,j,k), LP, L2G )
            L2G = Kstar298 * EXP( -H298_R*((1d0/tmpu(i,j,k))-INV_T0) ) &
                  * LP * R * tmpu(i,j,k)
            ! Washout fraction from Henry's law (Eq. 16, Jacob et al, 2000)
            WASHFRAC = L2G / ( 1d0 + L2G )
            ! Washout fraction / F from Eq. 14, Jacob et al, 2000
            ! Note: WASHFRAC_F_14 should match what's used for HNO3 (hma, 13aug2011)
            WASHFRAC_F_14 = 1d0 - EXP( -K_WASH * ( PP / F ) * cdt )
            ! Do not let the Henry's law washout fraction exceed
            IF ( WASHFRAC > WASHFRAC_F_14 ) THEN
              WASHFRAC = WASHFRAC_F_14
            ENDIF
           ELSE
            !------------------------
            ! T < 268K: No washout
            !------------------------
            WASHFRAC = 0d0
           ENDIF
        endif

!       Adjust du level:
        do n = 1, nbins
         if ( KIN ) then
            DC(n) = qa(i,j,k) * F * (1.-exp(-BT))
         else
            DC(n) = qa(i,j,k) * F * WASHFRAC
         endif
         if (DC(n).lt.0.) DC(n) = 0.
         qa(i,j,k) = qa(i,j,k)-DC(n)
         if (qa(i,j,k) .lt. 1.0E-32) & 
          qa(i,j,k) = 1.0E-32
         if( associated(flux) ) then
          flux(i,j) = flux(i,j)+DC(n)*pdog(i,j,k)/cdt
         endif
        end do

       end if
      end if                                    ! if ls washout  >>>

!-----------------------------------------------------------------------------
!  (3) CONVECTIVE RAINOUT:
!      Tracer loss by rainout = dd0 * F * exp(-B*dt)
!        where B = precipitation frequency,
!              F = fraction of grid box covered by precipitating clouds.
!-----------------------------------------------------------------------------

      if (qcv(k) .gt. 0.) then
       F  = F0_cv / (1. + F0_cv*B0_cv*XL_cv/(Qcv(k)*cdt/Td_cv))
       B  = B0_cv
       BT = B * Td_cv
       if (BT.gt.10.) BT = 10.               !< Avoid overflow >

!      Adjust du level: 
       do n = 1, nbins
        effRemoval = f_wet
        DC(n) = qa(i,j,k) * F * effRemoval * (1.-exp(-BT))
        if (DC(n).lt.0.) DC(n) = 0.
        qa(i,j,k) = qa(i,j,k)-DC(n)
        if (qa(i,j,k) .lt. 1.0E-32) qa(i,j,k) = 1.0E-32
       end do

!------  Flux down:  unit is kg. Including both ls and cv.
       do n = 1, nbins
        Fd(k,n) = Fd(k,n) + DC(n)*pdog(i,j,k)
       end do

      end if                                  ! if Qcv > 0   >>>

!-----------------------------------------------------------------------------
!  (4) CONVECTIVE WASHOUT:
!      Occurs when rain at this level is less than above.
!-----------------------------------------------------------------------------

      if (k.gt.LH .and. Qcv(k).ge.0.) then
       if (Qcv(k).lt.Qcv(k-1)) then
!-----  Find a maximum F overhead until the level where Qls<0.
        Qmx   = 0.
        do kk = k-1, LH, -1
         if (Qcv(kk).gt.0.) then
          Qmx = max(Qmx,Qcv(kk))
         else
          goto 444
         end if
        end do

 444    continue
        F = F0_cv / (1. + F0_cv*B0_cv*XL_cv/(Qmx*cdt/Td_cv))
        if (F.lt.0.01) F = 0.01
!-----------------------------------------------------------------------------
!  The following is to convert Q(k) from kgH2O/m3/sec to mm/sec in order
!  to use the Harvard formula.  Convert back to mixing ratio by multiplying
!  by rhoa.  Multiply by pdog gives kg/m2/s of precip.  Divide by density
!  of water (=1000 kg/m3) gives m/s of precip and multiply by 1000 gives
!  units of mm/s (omit the multiply and divide by 1000).
!-----------------------------------------------------------------------------

        Qd = Qmx / rhoa(i,j,k)*pdog(i,j,k)
        if (Qd.ge.50.) then
         B = 0.
        else
         B = Qd * 0.1
        end if
        BT = B * cdt
        if (BT.gt.10.) BT = 10.

!       Adjust du level:
        do n = 1, nbins
         DC(n) = qa(i,j,k) * F * (1.-exp(-BT))
         if (DC(n).lt.0.) DC(n) = 0.
         qa(i,j,k) = qa(i,j,k)-DC(n)
         if (qa(i,j,k) .lt. 1.0E-32) & 
          qa(i,j,k) = 1.0E-32
         if( associated(flux) ) then
          flux(i,j) = flux(i,j)+DC(n)*pdog(i,j,k)/cdt
         endif
        end do

       end if
      end if                                    ! if cv washout  >>>

!-----------------------------------------------------------------------------
!  (5) RE-EVAPORATION.  Assume that SO2 is re-evaporated as SO4 since it
!      has been oxidized by H2O2 at the level above. 
!-----------------------------------------------------------------------------
!     Add in the flux from above, which will be subtracted if reevaporation occurs
      if(k .gt. LH) then
       do n = 1, nbins
        Fd(k,n) = Fd(k,n) + Fd(k-1,n)
       end do

!      Is there evaporation in the currect layer?
       if (dpfli(i,j,k) .lt. 0.) then
!       Fraction evaporated = H2O(k)evap / H2O(next condensation level).
        if (dpfli(i,j,k-1) .gt. 0.) then

          A =  abs(  dpfli(i,j,k) / dpfli(i,j,k-1)   )
          if (A .gt. 1.) A = 1.

!         Adjust tracer in the level
          do n = 1, nbins
           DC(n) =  Fd(k-1,n) / pdog(i,j,k) * A
           qa(i,j,k) = qa(i,j,k) + DC(n)
           qa(i,j,k) = max(qa(i,j,k),1.e-32)
!          Adjust the flux out of the bottom of the layer
           Fd(k,n)  = Fd(k,n) - DC(n)*pdog(i,j,k)
          end do

        endif
       endif                                   ! if -moistq < 0
      endif
     end do  ! k

     do n = 1, nbins
      if( associated(flux) ) then
       flux = flux + Fd(km,n)/cdt
      endif
     end do

 100 continue
    end do   ! i
   end do    ! j

   deallocate(fd,     __STAT__)
   deallocate(dpfli,  __STAT__)
   deallocate(DC,     __STAT__)
   deallocate(c_h2o,  __STAT__)
   deallocate(cldliq, __STAT__)
   deallocate(cldice, __STAT__)
   deallocate(qls,    __STAT__)         
   deallocate(qcv,    __STAT__)
   deallocate(pdog,   __STAT__)
   deallocate(delz,   __STAT__)

   end subroutine MAML_WetRemoval


   end module MAML_WetRemovalMod

