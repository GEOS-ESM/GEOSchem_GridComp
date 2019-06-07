#include "MAPL_Generic.h"

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 610.1, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  MAM_SeasaltMod --- MAM Sea-salt processes and diagnostics
!
! !INTERFACE:
!

 module  MAM_SeasaltMod

! !USES:

   use ESMF

   use MAPL_Mod

   use MAPL_SimpleBundleMod

   use Chem_ConstMod,      only: grav    

   use SeasaltEmissionMod, only: SeasaltEmission

   use MAM_BaseMod
   use MAM3_DataMod
   use MAM7_DataMod

   implicit none

   private

! !PUBLIC TYPES:
!


!
! !PUBLIIC MEMBER FUNCTIONS:
!

   public MAM_SS_Emission
   public MAM_SS_Diagnostics

!
! !DESCRIPTION:
!
!  This module implements MAM Sea-salt (emission, etc) and 
!  diagnostic fields
!
! !REVISION HISTORY:
!
!  07 Sep 2011    A. Darmenov    Initial implementation.
!
!EOP
!-------------------------------------------------------------------------


 contains

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  MAM_SS_Emission --- The Emission Driver
!
! !INTERFACE:
!

   subroutine MAM_SS_Emission (self, import, export, qa, f_emiss, cdt, rc)

! !USES:

   implicit NONE

! !INPUT/OUTPUT PARAMETERS:
   type(MAPL_SimpleBundle), intent(inout) :: qa         ! interstitial aerosol tracer fields
   type(ESMF_State), intent(inout)        :: export     ! export fields

! !INPUT PARAMETERS:
   type(MAM_Scheme), intent(in)           :: self       ! MAM scheme/configuration
   real, intent(in)                       :: f_emiss    ! tuning parameter for the seasalt emissions

   type(ESMF_State), intent(inout)        :: import     ! import fields
   real,    intent(in)                    :: cdt        ! chemical timestep (secs)

! !OUTPUT PARAMETERS:
   integer, intent(out)                   :: rc         ! error return code:
                                                        !    0 - all is well
                                                        !    1 -
 
! !DESCRIPTION: This routine implements the Seasalt Emissions Driver. That
!               is, adds tendencies due to emission.
!
! !REVISION HISTORY:
!
!  07 Sep 2011    A. Darmenov    Initial implementation.
!
!EOP
!-------------------------------------------------------------------------

   character(len=*), parameter :: Iam = 'MAM_SS_Emission'

   integer :: STATUS
   integer :: i1, i2, j1, j2, k1, km, n, i
   integer :: ijl, ijkl
   real    :: qmin, qmax
   real    :: rUp, rLow
   real, pointer,     dimension(:,:) :: emission_mass, emission_num
   real, allocatable, dimension(:,:) :: f_grid_efficiency, f_sst_emis, tskin_c, w10m 
   real, allocatable, dimension(:,:) :: dqa_mass, dqa_num

   integer, parameter :: method = 3       ! seasalt emission scheme (hardwired)
   integer, parameter :: sstemisFlag = 2  ! SST correction parameterization (hardwired)

!  Mode parameters
!  ------------------------
   integer :: nmodes  ! number of modes with seasalt emission
   real, dimension(MAM_MAX_NUMBER_MODES) :: d_cutoff_low
   real, dimension(MAM_MAX_NUMBER_MODES) :: d_cutoff_up

   integer, dimension(MAM_MAX_NUMBER_MODES) :: mode
   character(len=ESMF_MAXSTR) :: mode_name(MAM_MAX_NUMBER_MODES)
   character(len=ESMF_MAXSTR) :: mmr_name, nmr_name
   character(len=ESMF_MAXSTR) :: emiss_name

!  Input fields from fvGCM
!  -----------------------
   real, pointer, dimension(:,:)   :: u10m, v10m, ustar, tskin, frocean, frseaice
   real, pointer, dimension(:,:,:) :: rhoa, ple, delp

!  Exports
!  -----------------------
   real, pointer, dimension(:,:) :: emission


!  Initialize local variables
!  --------------------------
   rc = 0

   ASSERT_(self%id == MAM7_SCHEME .or. self%id == MAM3_SCHEME) 
 
   if (self%id == MAM7_SCHEME) then
       nmodes = size(MAM7_SS_EMISSION_MODE_ID)

       mode(1:nmodes)         = MAM7_SS_EMISSION_MODE_ID
       d_cutoff_low(1:nmodes) = MAM7_SS_EMISSION_D_CUTOFF_LOW
       d_cutoff_up(1:nmodes)  = MAM7_SS_EMISSION_D_CUTOFF_UP

       mode_name(1:nmodes)    = MAM7_MODE_NAME(mode(1:nmodes))
   else
       nmodes = size(MAM3_SS_EMISSION_MODE_ID)

       mode(1:nmodes)         = MAM3_SS_EMISSION_MODE_ID
       d_cutoff_low(1:nmodes) = MAM3_SS_EMISSION_D_CUTOFF_LOW
       d_cutoff_up(1:nmodes)  = MAM3_SS_EMISSION_D_CUTOFF_UP

       mode_name(1:nmodes)    = MAM3_MODE_NAME(mode(1:nmodes))
   end if


!  Get Imports
!  --------------
   call MAPL_GetPointer(import, frocean,  'FROCEAN', __RC__)
   call MAPL_GetPointer(import, frseaice, 'FRACI',   __RC__)

   call MAPL_GetPointer(import, u10m,     'U10M',    __RC__)
   call MAPL_GetPointer(import, v10m,     'V10M',    __RC__)
   call MAPL_GetPointer(import, ustar,    'USTAR',   __RC__)
   call MAPL_GetPointer(import, tskin,    'TS',      __RC__)
   
   call MAPL_GetPointer(import, rhoa,     'AIRDENS', __RC__)
   call MAPL_GetPointer(import, ple,      'PLE',     __RC__)
   call MAPL_GetPointer(import, delp,     'DELP',    __RC__)


!  Local dimensions
!  ----------------
   i1 = lbound(rhoa, 1)
   i2 = ubound(rhoa, 1)
   j1 = lbound(rhoa, 2)
   j2 = ubound(rhoa, 2)
   k1 = lbound(rhoa, 3)
   km = ubound(rhoa, 3)

   ijl  = (i2 - i1 + 1) * (j2 - j1 + 1)
   ijkl = ijl * km

#ifdef DEBUG
   call write_parallel(trim(Iam) // '::DEBUG::Indexes:')
   call write_parallel((/i1, i2/), format='(("i1, i2 = ", (X2I3)))')
   call write_parallel((/j1, j2/), format='(("j1, j2 = ", (X2I3)))')
   call write_parallel((/k1, km/), format='(("k1, k2 = ", (X2I3)))')

   call write_parallel(trim(Iam) // '::DEBUG::Inputs:')
   call write_parallel(self%id,    format='(("model        = ", (I5)))')
   call write_parallel(f_emiss,    format='(("emiss factor = ", (F5.3)))')
#endif

!  Seasalt Emissions
!  -----------------
   allocate(emission_mass(i1:i2,j1:j2),     __STAT__)
   allocate(emission_num(i1:i2,j1:j2),      __STAT__)
   allocate(dqa_mass(i1:i2,j1:j2),          __STAT__)
   allocate(dqa_num(i1:i2,j1:j2),           __STAT__)
   allocate(f_grid_efficiency(i1:i2,j1:j2), __STAT__) 
   allocate(w10m(i1:i2,j1:j2),              __STAT__)
   allocate(f_sst_emis(i1:i2,j1:j2),        __STAT__ )

!  Define 10-m wind speed
   w10m = sqrt(u10m*u10m + v10m*v10m)

!  Define grid emission efficiency
   f_grid_efficiency = min(max(0.0, frocean - frseaice), 1.0)

!  Apply SST correction to emissions
   f_sst_emis = 1.0

   if (sstemisFlag == 1) then          ! SST correction folowing Jaegle et al. 2011
       f_sst_emis = 0.0

       allocate(tskin_c(i1:i2,j1:j2), __STAT__)
       tskin_c  = tskin - 273.15
       f_sst_emis = (0.3 + 0.1*tskin_c - 0.0076*tskin_c**2 + 0.00021*tskin_c**3)

       where(f_sst_emis < 0.0) f_sst_emis = 0.0
       deallocate(tskin_c, __STAT__)
   else if (sstemisFlag == 2) then     ! GEOS5 SST correction
       f_sst_emis = 0.0

       allocate(tskin_c(i1:i2,j1:j2), __STAT__)
       tskin_c  = tskin - 273.15
    
       where(tskin_c < -0.1) tskin_c = -0.1    ! temperature range (0, 36) C 
       where(tskin_c > 36.0) tskin_c = 36.0    !

       f_sst_emis = (-1.107211 -0.010681*tskin_c -0.002276*tskin_c**2 + 60.288927*1.0/(40.0 - tskin_c))
       where(f_sst_emis < 0.0) f_sst_emis = 0.0
       where(f_sst_emis > 7.0) f_sst_emis = 7.0

       deallocate(tskin_c, __STAT__)
   endif




   do n = 1, nmodes
       mmr_name   = 'SS_A_'  // trim(mode_name(n))  ! name of the mass mixing ratio
       nmr_name   = 'NUM_A_' // trim(mode_name(n))  ! name of the number mixing ratio
       emiss_name = 'SSEM'   // trim(mode_name(n))  ! name of the emission export

       rLow = 1.0e6 * d_cutoff_low(n) / 2           ! convert from [m] to [um]
       rUp  = 1.0e6 * d_cutoff_up(n)  / 2           ! convert from [m] to [um]
       
       emission_num  = 0.0
       emission_mass = 0.0

       call SeasaltEmission(rLow, rUp, method, w10m, ustar, &
                            emission_mass, emission_num, rc)

       emission_mass = f_emiss * f_grid_efficiency * f_sst_emis * emission_mass
       emission_num  = f_emiss * f_grid_efficiency * f_sst_emis * emission_num

#ifdef DEBUG
       call write_parallel('SS: mode ' // trim(mode_name(n)))
       
       call pmaxmin('SS: emission_mass  ', emission_mass, qmin, qmax, ijl, 1, 1.)
       call pmaxmin('SS: emission_number', emission_num,  qmin, qmax, ijl, 1, 1.)
#endif

       dqa_mass = emission_mass * cdt * grav / delp(:,:,km)
       dqa_num  = emission_num  * cdt * grav / delp(:,:,km)

       ! update the mass and number mixing ratios due to emission
       i = MAPL_SimpleBundleGetIndex(qa, mmr_name, 3, __RC__)
       qa%r3(i)%q(:,:,km) = qa%r3(i)%q(:,:,km) + dqa_mass

       i = MAPL_SimpleBundleGetIndex(qa, nmr_name, 3, __RC__)
       qa%r3(i)%q(:,:,km) = qa%r3(i)%q(:,:,km) + dqa_num

       call MAPL_GetPointer(export, emission, emiss_name, __RC__)
       if (associated(emission)) then
           emission = emission_mass
       endif
   end do


!  Clean up
!  --------
   deallocate(emission_mass,     __STAT__)
   deallocate(emission_num,      __STAT__)
   deallocate(dqa_mass,          __STAT__)
   deallocate(dqa_num,           __STAT__)
   deallocate(f_grid_efficiency, __STAT__)
   deallocate(f_sst_emis,        __STAT__)
   deallocate(w10m,              __STAT__)


   RETURN_(ESMF_SUCCESS)

 end subroutine MAM_SS_Emission



!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  MAM_SS_Diagnostics --- The Diagnostics Driver
!
! !INTERFACE:
!

   subroutine MAM_SS_Diagnostics (self, import, export, qa, cdt, rc)

! !USES:

   implicit NONE

! !INPUT/OUTPUT PARAMETERS:
   type(MAPL_SimpleBundle), intent(inout) :: qa         ! interstitial aerosol tracer fields
   type(ESMF_State), intent(inout)        :: export     ! export fields

! !INPUT PARAMETERS:
   type(MAM_Scheme), intent(in)           :: self      ! MAM scheme/configuration

   type(ESMF_State), intent(inout)        :: import     ! import fields
   real,    intent(in)                    :: cdt        ! chemical timestep (secs)

! !OUTPUT PARAMETERS:
   integer, intent(out)                   :: rc         ! error return code:
                                                        !    0 - all is well
                                                        !    1 -
 
! !DESCRIPTION: This routine calculates a number of diagnostic fields.
!
! !REVISION HISTORY:
!
!  13 Oct 2011    A. Darmenov    Initial implementation.
!
!EOP
!-------------------------------------------------------------------------

   character(len=*), parameter :: Iam = 'MAM_SS_Diagnostics'

   integer :: STATUS
   integer :: i1, i2, j1, j2, k1, km, n, i, k
   integer :: ijl, ijkl
   real    :: qmin, qmax

!  Mode parameters
!  ------------------------
   integer :: nmodes  ! number of modes with dust emission

   integer, dimension(MAM_MAX_NUMBER_MODES) :: mode
   character(len=ESMF_MAXSTR) :: mode_name(MAM_MAX_NUMBER_MODES)
   character(len=ESMF_MAXSTR) :: mmr_name, nmr_name

!  Input fields from fvGCM
!  -----------------------
   real, pointer, dimension(:,:,:) :: rhoa, ple, delp
   real, pointer, dimension(:,:,:) :: u, v

!  Exports - diagnostic fields
!  ---------------------------
   real, pointer, dimension(:,:)   :: sfcmass      ! surface mass concentration, kg m-3
   real, pointer, dimension(:,:)   :: sfcmass25    ! surface PM2.5 mass concentration, kg m-3
   real, pointer, dimension(:,:)   :: colmass      ! column integrated mass density, kg m-2
   real, pointer, dimension(:,:)   :: colmass25    ! column integrated PM2.5 mass density, kg m-2

   real, pointer, dimension(:,:)   :: fluxu        ! Column mass flux in x direction
   real, pointer, dimension(:,:)   :: fluxv        ! Column mass flux in y direction

   real, pointer, dimension(:,:,:) :: conc         ! mass concentration, kg m-3
   real, pointer, dimension(:,:,:) :: mass         ! mass mixing ratio, kg kg-1
   real, pointer, dimension(:,:,:) :: mass25       ! PM2.5 mass mixing ratio, kg kg-1


!  Initialize local variables
!  --------------------------
   rc = 0

   ASSERT_(self%id == MAM7_SCHEME .or. self%id == MAM3_SCHEME) 
 
   if (self%id == MAM7_SCHEME) then
       nmodes = size(MAM7_SS_EMISSION_MODE_ID)
       mode(1:nmodes) = MAM7_SS_EMISSION_MODE_ID

       mode_name(1:nmodes) = MAM7_MODE_NAME(mode(1:nmodes))
   else
       nmodes = size(MAM3_SS_EMISSION_MODE_ID)
       mode(1:nmodes) = MAM3_SS_EMISSION_MODE_ID

       mode_name(1:nmodes) = MAM3_MODE_NAME(mode(1:nmodes))
   end if


!  Get Imports
!  --------------
   call MAPL_GetPointer(import, rhoa,      'AIRDENS',   __RC__)
   call MAPL_GetPointer(import, ple,       'PLE',       __RC__)
   call MAPL_GetPointer(import, delp,      'DELP',      __RC__)
   call MAPL_GetPointer(import, u,         'U',         __RC__)
   call MAPL_GetPointer(import, v,         'V',         __RC__)

!  Get Exports
!  --------------
   call MAPL_GetPointer(export, sfcmass,   'SSSMASS',   __RC__)
   call MAPL_GetPointer(export, sfcmass25, 'SSSMASS25', __RC__)
   call MAPL_GetPointer(export, colmass,   'SSCMASS',   __RC__)
   call MAPL_GetPointer(export, colmass25, 'SSCMASS25', __RC__)
   
   call MAPL_GetPointer(export, fluxu,     'SSFLUXU',   __RC__)
   call MAPL_GetPointer(export, fluxv,     'SSFLUXV',   __RC__)

   call MAPL_GetPointer(export, conc,      'SSCONC',    __RC__)
   call MAPL_GetPointer(export, mass,      'SSMASS',    __RC__)
   call MAPL_GetPointer(export, mass25,    'SSMASS25',  __RC__)


!  Local dimensions
!  ----------------
   i1 = lbound(rhoa, 1)
   i2 = ubound(rhoa, 1)
   j1 = lbound(rhoa, 2)
   j2 = ubound(rhoa, 2)
   k1 = lbound(rhoa, 3)
   km = ubound(rhoa, 3)

   ijl  = (i2 - i1 + 1) * (j2 - j1 + 1)
   ijkl = ijl * km

#ifdef DEBUG
   call write_parallel(trim(Iam) // '::DEBUG::Indexes:')
   call write_parallel((/i1, i2/), format='(("i1, i2 = ", (X2I3)))')
   call write_parallel((/j1, j2/), format='(("j1, j2 = ", (XI3)))')
   call write_parallel((/k1, km/), format='(("k1, k2 = ", (XI3)))')

   call write_parallel(trim(Iam) // '::DEBUG::Inputs:')
   call write_parallel(self%id,    format='(("model = ", (I5)))')
#endif


!  Initialize diagnostic fields
!  ----------------------------
   if (associated(sfcmass))    sfcmass   = 0.0
   if (associated(sfcmass25))  sfcmass25 = 0.0
   if (associated(colmass))    colmass   = 0.0
   if (associated(colmass25))  colmass25 = 0.0

   if (associated(fluxu))      fluxu     = 0.0
   if (associated(fluxv))      fluxv     = 0.0

   if (associated(conc))       conc      = 0.0
   if (associated(mass))       mass      = 0.0
   if (associated(mass25))     mass25    = 0.0


!  Calculate diagnostic fields
!  ---------------------------
   do n = 1, nmodes
       mmr_name   = 'SS_A_'  // trim(mode_name(n))  ! name of the mass mixing ratio
       nmr_name   = 'NUM_A_' // trim(mode_name(n))  ! name of the number mixing ratio

       i = MAPL_SimpleBundleGetIndex(qa, mmr_name, 3, __RC__)


       if (associated(sfcmass)) then 
           sfcmass(:,:) = sfcmass(:,:) + qa%r3(i)%q(:,:,km) * rhoa(:,:,km)
       end if

       if (associated(sfcmass25)) then ! placeholder for now
           sfcmass25(:,:) = sfcmass25(:,:) + 0.0 * qa%r3(i)%q(:,:,km) * rhoa(:,:,km)
       end if
  
       if (associated(colmass)) then
           do k = 1, km
               colmass(:,:) = colmass(:,:) + qa%r3(i)%q(:,:,k) * delp(:,:,k)/grav
           end do
       end if

       if (associated(colmass25)) then   ! placeholder for now
           do k = 1, km
               colmass25(:,:) = colmass25(:,:) + 0.0 * qa%r3(i)%q(:,:,k) * delp(:,:,k)/grav
           end do
       end if

       if (associated(fluxu)) then
           do k = 1, km
               fluxu(:,:) = fluxu(:,:) + u(:,:,k) * qa%r3(i)%q(:,:,k) * delp(:,:,k)/grav
           end do
       end if

       if (associated(fluxv)) then
           do k = 1, km
               fluxv(:,:) = fluxv(:,:) + v(:,:,k) * qa%r3(i)%q(:,:,k) * delp(:,:,k)/grav
           end do
       end if

       if (associated(conc)) then
           conc = conc + qa%r3(i)%q * rhoa
       end if

       if (associated(mass)) then  
           mass = mass + qa%r3(i)%q
       end if
 
       if (associated(mass25)) then     ! placeholder for now
           mass25 = mass25 + 0.0 * qa%r3(i)%q
       end if

   end do


!  Clean up
!  --------

   RETURN_(ESMF_SUCCESS)

 end subroutine MAM_SS_Diagnostics



 end module  MAM_SeasaltMod
