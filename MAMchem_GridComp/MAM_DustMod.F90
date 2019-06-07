#include "MAPL_Generic.h"

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 610.1, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  MAM_DustMod --- MAM dust processes and diagnostics
!
! !INTERFACE:
!

   module  MAM_DustMod

! !USES:

   USE ESMF
   USE MAPL_Mod

   USE MAPL_SimpleBundleMod

   use Chem_ConstMod,   only: grav

   use DustEmissionMod, only: MAM_DustEmissionGOCART, MAM_DustEmission

   use MAM_BaseMod
   use MAM3_DataMod
   use MAM7_DataMod

   implicit none

! !PUBLIC TYPES:
!
   PRIVATE

!
! !PUBLIIC MEMBER FUNCTIONS:
!

   PUBLIC  MAM_DU_Emission
   PUBLIC  MAM_DU_Diagnostics

!
! !DESCRIPTION:
!
!  This module implements MAM dust processes (emission, etc) and 
!  diagnostic fields
!
! !REVISION HISTORY:
!
!  07 Sep 2011    A. Darmenov    Initial implementation.
!
!EOP
!-------------------------------------------------------------------------


CONTAINS

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  MAM_DU_Emission --- The Emission Driver
!
! !INTERFACE:
!

   subroutine MAM_DU_Emission (self, import, export, qa, f_emiss, cdt, rc)

! !USES:

   implicit NONE

! !INPUT/OUTPUT PARAMETERS:
   type(MAPL_SimpleBundle), intent(inout) :: qa         ! interstitial aerosol tracer fields
   type(ESMF_State), intent(inout)        :: export     ! export fields

! !INPUT PARAMETERS:
   type(MAM_Scheme), intent(in)           :: self       ! MAM scheme/configuration
   real, intent(in)                       :: f_emiss    ! tuning parameter for the dust emissions

   type(ESMF_State), intent(inout)        :: import     ! import fields
   real,    intent(in)                    :: cdt        ! chemical timestep (secs)

! !OUTPUT PARAMETERS:
   integer, intent(out)                   :: rc         ! error return code:
                                                        !    0 - all is well
                                                        !    1 -
 
! !DESCRIPTION: This routine implements the Dust Emissions Driver. That
!               is, adds tendencies due to emission.
!
! !REVISION HISTORY:
!
!  07 Sep 2011    A. Darmenov    Initial implementation.
!
!EOP
!-------------------------------------------------------------------------

   character(len=*), parameter :: Iam = 'MAM_DU_Emission'

   integer :: STATUS
   integer :: i1, i2, j1, j2, k1, km, n, i
   integer :: ijl, ijkl
   real    :: qmin, qmax
   real    :: rUp, rLow
   real, pointer, dimension(:,:) :: emission_total
   real, pointer, dimension(:,:) :: emission_mass, emission_num
   real, pointer, dimension(:,:) :: dqa_mass, dqa_num

!  Mode parameters
!  ------------------------
   integer :: nmodes  ! number of modes with dust emission
   real, dimension(MAM_MAX_NUMBER_MODES) :: d_cutoff_low
   real, dimension(MAM_MAX_NUMBER_MODES) :: d_cutoff_up

   integer, dimension(MAM_MAX_NUMBER_MODES) :: mode
   character(len=ESMF_MAXSTR) :: mode_name(MAM_MAX_NUMBER_MODES)
   character(len=ESMF_MAXSTR) :: mmr_name, nmr_name
   character(len=ESMF_MAXSTR) :: emiss_name

!  Input fields from fvGCM
!  -----------------------
   real, pointer, dimension(:,:)   :: fsrc, oro, u10m, v10m
   real, pointer, dimension(:,:)   :: fraclake, gwettop
   real, pointer, dimension(:,:,:) :: rhoa, ple, delp

!  Exports
!  -----------------------
   real, pointer, dimension(:,:) :: emission


!  Initialize local variables
!  --------------------------
   rc = 0

   ASSERT_(self%id == MAM7_SCHEME .or. self%id == MAM3_SCHEME) 
 
   if (self%id == MAM7_SCHEME) then
       nmodes = size(MAM7_DU_EMISSION_MODE_ID)

       mode(1:nmodes)         = MAM7_DU_EMISSION_MODE_ID
       d_cutoff_low(1:nmodes) = MAM7_DU_EMISSION_D_CUTOFF_LOW
       d_cutoff_up(1:nmodes)  = MAM7_DU_EMISSION_D_CUTOFF_UP

       mode_name(1:nmodes)    = MAM7_MODE_NAME(mode(1:nmodes))
   else
       nmodes = size(MAM3_DU_EMISSION_MODE_ID)

       mode(1:nmodes)         = MAM3_DU_EMISSION_MODE_ID
       d_cutoff_low(1:nmodes) = MAM3_DU_EMISSION_D_CUTOFF_LOW
       d_cutoff_up(1:nmodes)  = MAM3_DU_EMISSION_D_CUTOFF_UP

       mode_name(1:nmodes)    = MAM3_MODE_NAME(mode(1:nmodes))
   end if


!  Get Imports
!  --------------
   call MAPL_GetPointer(import, fsrc,     'GINOUX_DU', __RC__)
   call MAPL_GetPointer(import, oro,      'LWI',       __RC__)
   call MAPL_GetPointer(import, u10m,     'U10M',      __RC__)
   call MAPL_GetPointer(import, v10m,     'V10M',      __RC__)
   call MAPL_GetPointer(import, fraclake, 'FRLAKE',    __RC__)
   call MAPL_GetPointer(import, gwettop,  'WET1',      __RC__)
   call MAPL_GetPointer(import, rhoa,     'AIRDENS',   __RC__)
   call MAPL_GetPointer(import, ple,      'PLE',       __RC__)
   call MAPL_GetPointer(import, delp,     'DELP',      __RC__)


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
   call write_parallel(trim(Iam) // '::DEBUG(before)::Indexes:')
   call write_parallel((/i1, i2/), format='(("i1, i2 = ", (X2I3)))')
   call write_parallel((/j1, j2/), format='(("j1, j2 = ", (X2I3)))')
   call write_parallel((/k1, km/), format='(("k1, k2 = ", (X2I3)))')
   
   do i = 1, qa%n3d
       call pmaxmin('CAM:qa:'//trim(qa%r3(i)%name)//'   : ', &
                     qa%r3(i)%q(i1:i2,j1:j2,k1:km), qmin, qmax, ijl, km, 1.)
   end do
#endif


#ifdef DEBUG
   call write_parallel(trim(Iam) // '::DEBUG::Indexes:')
   call write_parallel((/i1, i2/), format='(("i1, i2 = ", (X2I3)))')
   call write_parallel((/j1, j2/), format='(("j1, j2 = ", (X2I3)))')
   call write_parallel((/k1, km/), format='(("k1, k2 = ", (X2I3)))')

   call write_parallel(trim(Iam) // '::DEBUG::Inputs:')
   call write_parallel(self%id,    format='(("model        = ", (I5)))')
   call write_parallel(f_emiss,    format='(("emiss factor = ", (F5.3)))')

   call pmaxmin('DU: fsrc     ', fsrc,     qmin, qmax, ijl, k1, 1.)
   call pmaxmin('DU: oro      ', oro,      qmin, qmax, ijl, k1, 1.)
   call pmaxmin('DU: u10m     ', u10m,     qmin, qmax, ijl, k1, 1.)
   call pmaxmin('DU: v10m     ', v10m,     qmin, qmax, ijl, k1, 1.)
   call pmaxmin('DU: fraclake ', fraclake, qmin, qmax, ijl, k1, 1.)
   call pmaxmin('DU: gwettop  ', gwettop,  qmin, qmax, ijl, k1, 1.)

   call pmaxmin('DU: rhoa     ', rhoa,     qmin, qmax, ijl, km, 1.)
   call pmaxmin('DU: ple      ', ple,      qmin, qmax, ijl, km, 1.)
   call pmaxmin('DU: delp     ', delp,     qmin, qmax, ijl, km, 1.)
#endif

!  Dust Emissions
!  --------------
   allocate(emission_total(i1:i2,j1:j2), __STAT__)
   allocate(emission_mass(i1:i2,j1:j2),  __STAT__)
   allocate(emission_num(i1:i2,j1:j2),   __STAT__)
   allocate(dqa_mass(i1:i2,j1:j2),       __STAT__)
   allocate(dqa_num(i1:i2,j1:j2),        __STAT__)


   emission_total = 0.0
   call MAM_DustEmissionGOCART(i1, i2, j1, j2, km, &
                               fraclake, gwettop, oro, u10m, v10m, &
                               emission_total, rc)

   ! apply the dust emission tuning coefficient [kg s2 m-5] and Ginoux dust source function
   emission_total = (f_emiss * 1e-9) * fsrc * emission_total

   do n = 1, nmodes
       mmr_name   = 'DU_A_'  // trim(mode_name(n))  ! name of the mass mixing ratio
       nmr_name   = 'NUM_A_' // trim(mode_name(n))  ! name of the number mixing ratio
       emiss_name = 'DUEM'   // trim(mode_name(n))  ! name of the emission export

       rLow = d_cutoff_low(n) / 2  ! [m]
       rUp  = d_cutoff_up(n)  / 2  ! [m]
       
       emission_num  = 0.0
       emission_mass = 0.0

       call MAM_DustEmission(i1, i2, j1, j2, km, &
                             rLow, rUp, &
                             emission_total, &
                             emission_mass, emission_num, rc)

       emission_mass = emission_mass
       emission_num  = emission_num

#ifdef DEBUG
       call pmaxmin('DU: emission_total ', emission_total, qmin, qmax, ijl, 1, 1.)

       call write_parallel('DU: mode ' // trim(mode_name(n)))
       call pmaxmin('DU: emission_mass  ', emission_mass,  qmin, qmax, ijl, 1, 1.)
       call pmaxmin('DU: emission_number', emission_num,   qmin, qmax, ijl, 1, 1.)
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


#ifdef DEBUG
   call write_parallel(trim(Iam) // '::DEBUG(after)::Indexes:')
   call write_parallel((/i1, i2/), format='(("i1, i2 = ", (X2I3)))')
   call write_parallel((/j1, j2/), format='(("j1, j2 = ", (X2I3)))')
   call write_parallel((/k1, km/), format='(("k1, k2 = ", (X2I3)))')
   
   do i = 1, qa%n3d
       call pmaxmin('CAM:qa:'//trim(qa%r3(i)%name)//'   : ', &
                     qa%r3(i)%q(i1:i2,j1:j2,k1:km), qmin, qmax, ijl, km, 1.)
   end do
#endif


!  Clean up
!  --------
   deallocate(emission_total, __STAT__)
   deallocate(emission_mass,  __STAT__)
   deallocate(emission_num,   __STAT__)
   deallocate(dqa_mass,       __STAT__)
   deallocate(dqa_num,        __STAT__)


   RETURN_(ESMF_SUCCESS)

 end subroutine MAM_DU_Emission


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  MAM_DU_Diagnostics --- The Diagnostics Driver
!
! !INTERFACE:
!

   subroutine MAM_DU_Diagnostics (self, import, export, qa, cdt, rc)

! !USES:

   implicit NONE

! !INPUT/OUTPUT PARAMETERS:
   type(MAPL_SimpleBundle), intent(inout) :: qa         ! interstitial aerosol tracer fields
   type(ESMF_State), intent(inout)        :: export     ! export fields

! !INPUT PARAMETERS:
   type(MAM_Scheme), intent(in)           :: self       ! MAM scheme/configuration

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

   character(len=*), parameter :: Iam = 'MAM_DU_Diagnostics'

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
       nmodes = size(MAM7_DU_EMISSION_MODE_ID)
       mode(1:nmodes) = MAM7_DU_EMISSION_MODE_ID

       mode_name(1:nmodes) = MAM7_MODE_NAME(mode(1:nmodes))
   else
       nmodes = size(MAM3_DU_EMISSION_MODE_ID)
       mode(1:nmodes) = MAM3_DU_EMISSION_MODE_ID

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
   call MAPL_GetPointer(export, sfcmass,   'DUSMASS',   __RC__)
   call MAPL_GetPointer(export, sfcmass25, 'DUSMASS25', __RC__)
   call MAPL_GetPointer(export, colmass,   'DUCMASS',   __RC__)
   call MAPL_GetPointer(export, colmass25, 'DUCMASS25', __RC__)
   
   call MAPL_GetPointer(export, fluxu,     'DUFLUXU',   __RC__)
   call MAPL_GetPointer(export, fluxv,     'DUFLUXV',   __RC__)

   call MAPL_GetPointer(export, conc,      'DUCONC',    __RC__)
   call MAPL_GetPointer(export, mass,      'DUMASS',    __RC__)
   call MAPL_GetPointer(export, mass25,    'DUMASS25',  __RC__)


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
       mmr_name   = 'DU_A_'  // trim(mode_name(n))  ! name of the mass mixing ratio
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

 end subroutine MAM_DU_Diagnostics


 end module  MAM_DustMod
