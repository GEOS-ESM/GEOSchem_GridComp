#include "MAPL_Generic.h"

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 610.1, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  MAM_SulfateMod --- MAM Sulfate (SO4) processes and diagnostics
!
! !INTERFACE:
!

   module  MAM_SulfateMod

! !USES:

   USE ESMF
   USE MAPL_Mod

   USE MAPL_SimpleBundleMod

   use Chem_ConstMod, only: grav, undef
   use Chem_UtilMod,  only: Chem_BiomassDiurnal

   use MAM_BaseMod
   use MAM3_DataMod
   use MAM7_DataMod

   implicit none

! !PUBLIC TYPES:
!
   PRIVATE
   real, private, parameter :: pi = MAPL_PI

!
! !PUBLIIC MEMBER FUNCTIONS:
!

   PUBLIC  MAM_SO4_Emission
   PUBLIC  MAM_SO4_Diagnostics

!
! !DESCRIPTION:
!
!  This module implements MAM Sulfate (SO4) processes (emission, etc) and 
!  diagnostic fields
!
! !REVISION HISTORY:
!
!  01 Aug 2016    A. Darmenov    Initial implementation.
!
!EOP
!-------------------------------------------------------------------------


CONTAINS

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  MAM_SO4_Emission --- The Emission Driver
!
! !INTERFACE:
!

   subroutine MAM_SO4_Emission (self, import, export, qa, cdt, rc)

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
 
! !DESCRIPTION: This routine implements the Sulfate (SO4) Emissions Driver. That
!               is, adds tendencies due to emission.
!
! !REVISION HISTORY:
!
!  11 Jul 2012    A. Darmenov    Initial implementation.
!
!EOP
!-------------------------------------------------------------------------

   character(len=*), parameter :: Iam = 'MAM_SO4_Emission'

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

   integer, dimension(MAM_MAX_NUMBER_MODES) :: mode
   character(len=ESMF_MAXSTR) :: mode_name
   character(len=ESMF_MAXSTR) :: mmr_name, nmr_name
   character(len=ESMF_MAXSTR) :: emiss_name

!  Input fields from fvGCM
!  -----------------------
   real, pointer, dimension(:,:)   :: gwettop
   real, pointer, dimension(:,:,:) :: rhoa, ple, delp

!  Input fields from ExtData
!  -------------------------
   real, pointer, dimension(:,:) :: emiss_sh

!  Exports
!  -----------------------
   real, pointer, dimension(:,:) :: emission

!
!  Parameters of primary aerosol emissions 
!  ---------------------------------------
   type PAE
       integer :: mode_id        ! mode ID

       real    :: weight         ! weight by mass
       real    :: sigma          ! geometric standard deviation
       real    :: diameter       ! geometric mean diameter of number size distribution

       real, pointer, dimension(:,:)   :: emission  => null()
       real, pointer, dimension(:,:,:) :: injection => null()
   end type


   ! note that the SO4 from ship emissions is emitted with D_mass_emiss = 0.261 micrometers
   type(PAE), parameter :: pae_sh = PAE(1, 1.0, 0.0, 0.261, null(), null())

   real :: D_emiss_sh, f_sh

!  Initialize local variables
!  --------------------------
   rc = 0

   D_emiss_sh = pae_sh%diameter

   f_sh = 1 / (pi/6 * D_emiss_sh**3)


   ASSERT_(self%id == MAM7_SCHEME .or. self%id == MAM3_SCHEME) 
 
   if (self%id == MAM7_SCHEME) then
       mode_name = MAM7_ACCUMULATION_MODE_NAME
   else
       mode_name = MAM3_ACCUMULATION_MODE_NAME
   end if


!  Get Imports
!  --------------
   call MAPL_GetPointer(import, gwettop,  'WET1',                __RC__)
   call MAPL_GetPointer(import, rhoa,     'AIRDENS',             __RC__)
   call MAPL_GetPointer(import, ple,      'PLE',                 __RC__)
   call MAPL_GetPointer(import, delp,     'DELP',                __RC__)

   call MAPL_GetPointer(import, emiss_sh, 'SO4_EMIS_SHIP',       __RC__)



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
   call write_parallel(self%id, format='(("model     = ", (I5)))')

   call pmaxmin('SO4: gwettop  ', gwettop,  qmin, qmax, ijl, k1, 1.)

   call pmaxmin('SO4: rhoa     ', rhoa,     qmin, qmax, ijl, km, 1.)
   call pmaxmin('SO4: ple      ', ple,      qmin, qmax, ijl, km, 1.)
   call pmaxmin('SO4: delp     ', delp,     qmin, qmax, ijl, km, 1.)
#endif

!  Sulfate (SO4) Emissions
!  -----------------------
   allocate(emission_total(i1:i2,j1:j2), __STAT__)
   allocate(emission_mass(i1:i2,j1:j2),  __STAT__)
   allocate(emission_num(i1:i2,j1:j2),   __STAT__)
   allocate(dqa_mass(i1:i2,j1:j2),       __STAT__)
   allocate(dqa_num(i1:i2,j1:j2),        __STAT__)

   emission_total = 0.0

   mmr_name   = 'SU_A_'  // trim(mode_name)  ! name of the mass mixing ratio
   nmr_name   = 'NUM_A_' // trim(mode_name)  ! name of the number mixing ratio
   emiss_name = 'SUEM'   // trim(mode_name)  ! name of the emission export

   emission_num  = 0.0
   emission_mass = 0.0


   ! set undefined to 0
   where ((emiss_sh - undef) > abs(emiss_sh)*epsilon(emiss_sh)) emiss_sh = 0.0


   emission_mass = emiss_sh

   emission_num  = f_sh * emiss_sh
                       


#ifdef DEBUG
   call pmaxmin('SO4: emission_total ', emission_total, qmin, qmax, ijl, 1, 1.)

   call write_parallel('SO4: mode ' // trim(mode_name))
   call pmaxmin('SO4: emission_mass  ', emission_mass,  qmin, qmax, ijl, 1, 1.)
   call pmaxmin('SO4: emission_number', emission_num,   qmin, qmax, ijl, 1, 1.)
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

 end subroutine MAM_SO4_Emission


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  MAM_SO4_Diagnostics --- The Diagnostics Driver
!
! !INTERFACE:
!

   subroutine MAM_SO4_Diagnostics (self, import, export, qa, cdt, rc)

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
!  11 Jul 2012    A. Darmenov    Initial implementation.
!
!EOP
!-------------------------------------------------------------------------

   character(len=*), parameter :: Iam = 'MAM_SO4_Diagnostics'

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
#if (0) 
   if (self%id == MAM7_SCHEME) then
       nmodes = size(MAM7_SO4_EMISSION_MODE_ID)
       mode(1:nmodes) = MAM7_SO4_EMISSION_MODE_ID

       mode_name(1:nmodes) = MAM7_MODE_NAME(mode(1:nmodes))
   else
       nmodes = size(MAM3_SO4_EMISSION_MODE_ID)
       mode(1:nmodes) = MAM3_SO4_EMISSION_MODE_ID

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
   call MAPL_GetPointer(export, sfcmass,   'SUSMASS',   __RC__)
   call MAPL_GetPointer(export, sfcmass25, 'SUSMASS25', __RC__)
   call MAPL_GetPointer(export, colmass,   'SUCMASS',   __RC__)
   call MAPL_GetPointer(export, colmass25, 'SUCMASS25', __RC__)
   
   call MAPL_GetPointer(export, fluxu,     'SUFLUXU',   __RC__)
   call MAPL_GetPointer(export, fluxv,     'SUFLUXV',   __RC__)

   call MAPL_GetPointer(export, conc,      'SUCONC',    __RC__)
   call MAPL_GetPointer(export, mass,      'SUMASS',    __RC__)
   call MAPL_GetPointer(export, mass25,    'SUMASS25',  __RC__)


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
   call write_parallel(self%id,    format='(("model  = ", (I5)))')
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
       mmr_name   = 'SU_A_'  // trim(mode_name(n))  ! name of the mass mixing ratio
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
#endif

!  Clean up
!  --------

   RETURN_(ESMF_SUCCESS)

 end subroutine MAM_SO4_Diagnostics


 end module  MAM_SulfateMod
