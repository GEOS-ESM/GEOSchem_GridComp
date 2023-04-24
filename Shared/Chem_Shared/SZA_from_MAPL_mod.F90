module SZA_from_MAPL_mod

#include "MAPL_Generic.h"

!----------------------------------------------------------------------
!BOP

! !MODULE: 

!    SZA_from_MAPL_mod -- Call MAPL routine for Solar Zenith Angle

! !USES:

use ESMF
use MAPL

!use GEOS_UtilsMod
!
!use, intrinsic :: iso_fortran_env, only: REAL64


   IMPLICIT NONE
   INTEGER, PARAMETER :: DBL = KIND(0.00D+00)

   PRIVATE

  
! PUBLIC MEMBER FUNCTIONS:
!

   PUBLIC :: compute_SZA

   interface compute_SZA
      module procedure compute_sza_main         ! args include LONS, LATS, ORBIT
      module procedure compute_sza_wrapper_1    ! args include GC
   end interface compute_SZA



! !PARAMETERS:


! !DESCRIPTION:

!  This module provides a routine to determine Solar Zenith Angle
!  using the MAPL routine MAPL_SunGetInsolation.
!
! !REVISION HISTORY:
!
! March 15 2023 - Mike Manyin - Adapted from GMI photolysis code
!
!EOP
!-------------------------------------------------------------------------


contains


!-----------------------------------------------------------------------
!BOP
! !IROUTINE: compute_sza_main
!
! !DESCRIPTION:
!  
!  ORIGIN AND CONTACT
!
! !REVISION HISTORY:
! 15 Mar 2023 Manyin     First crack
!EOP
!-----------------------------------------------------------------------

subroutine compute_sza_main ( LONS, LATS, ORBIT, CLOCK, tdt, label, SZA, RC )

   REAL,                 intent(in)  :: LONS (:,:) ! radians
   REAL,                 intent(in)  :: LATS (:,:) ! radians
   TYPE(MAPL_SunOrbit),  intent(in)  :: ORBIT
   TYPE(ESMF_Clock),     intent(in)  :: CLOCK
   REAL,                 intent(in)  :: tdt        ! caller's timestep (sec)
   CHARACTER(LEN=*),     intent(in)  :: label      ! name of caller
   REAL(KIND=DBL),       intent(out) :: SZA (:,:)  ! degrees
   INTEGER, OPTIONAL,    intent(out) :: RC

! LOCAL
   REAL, ALLOCATABLE ::   ZTH(:,:)
   REAL, ALLOCATABLE ::   SLR(:,:)
   REAL, ALLOCATABLE ::  ZTHP(:,:)   ! if SZA becomes Single Precision, we can use it as ZTHP

   INTEGER :: IM, JM   ! array dimensions

   REAL :: pi,radToDeg

   TYPE (ESMF_TimeInterval)    :: CALLER_timestep
   TYPE (ESMF_TimeInterval)    :: MAPL_timestep
   TYPE (ESMF_Time)            :: CURRENTTIME
   TYPE (ESMF_Time)            :: SZA_start_time   ! compute average SZA starting at this time
   TYPE (ESMF_Time)            :: SZA_midpoint

   LOGICAL :: verbose_time   ! To see details on SZA time averaging

   INTEGER :: STATUS
   CHARACTER(LEN=*), PARAMETER :: IAm = 'compute_sza_main'

   verbose_time = .TRUE.


   IM = SIZE(LONS,1)
   JM = SIZE(LONS,2)

   ALLOCATE( ZTH(IM,JM), SLR(IM,JM), ZTHP(IM,JM),  __STAT__ )

   call ESMF_ClockGet(CLOCK, TIMESTEP=MAPL_timestep, currTIME=CURRENTTIME, __RC__ )

   IF( verbose_time .AND. MAPL_AM_I_ROOT() ) THEN
     call ESMF_TimePrint(CURRENTTIME, preString="CURRENTTIME = ", __RC__ )
     print *, "MAPL_timestep = "
     call ESMF_TimeIntervalPrint(MAPL_timestep, options="string", __RC__ )
   ENDIF

   call ESMF_TimeIntervalSet(CALLER_timestep, s=INT(tdt+0.1), __RC__ )

   IF( verbose_time .AND. MAPL_AM_I_ROOT() ) THEN
     print *, TRIM(label)//" timestep = "
     call ESMF_TimeIntervalPrint(CALLER_timestep, options="string", __RC__ )
     print *, "computing SZA w/ tdt = ", tdt
   ENDIF

   ! We want a starting time = MAPL time + 1/2 MAPL timestep - 1/2 CALLER timestep
   ! We want a time interval == CALLER timestep

   ! Position SZA_midpoint to be midpoint of MAPL_timestep
   SZA_midpoint = CURRENTTIME + (MAPL_timestep/2)

   ! Position SZA_start_time to be half of a CALLER timestep earlier
   SZA_start_time = SZA_midpoint - (CALLER_timestep/2)

   IF( verbose_time .AND. MAPL_AM_I_ROOT() ) THEN
     call ESMF_TimePrint(SZA_start_time,                 preString=TRIM(label)//" SZA averaging start_time = ", __RC__ )
     call ESMF_TimePrint(SZA_start_time+CALLER_timestep, preString=TRIM(label)//" SZA averaging   end_time = ", __RC__ )
   ENDIF


   ! For more complete calling sequence, see SOLAR GridComp

   CALL MAPL_SunGetInsolation(           &
           LONS     = LONS,              &
           LATS     = LATS,              &
           ORBIT    = ORBIT,             &
           ZTH      = ZTH,               &
           SLR      = SLR,               &
           CURRTIME = SZA_start_time,    &
           INTV     = CALLER_timestep,   &
           ZTHP     = ZTHP,              &
           __RC__ )


   pi = 4.00*ATAN(1.00)   !  LATER replace with the MAPL version
   radToDeg = 180.00/pi

   SZA = ACOS( ZTHP ) * radToDeg

   DEALLOCATE( ZTH, SLR, ZTHP, __STAT__ )

   RETURN_(ESMF_SUCCESS)

end subroutine compute_sza_main

      
!-----------------------------------------------------------------------
!BOP
! !IROUTINE: compute_sza_wrapper_1
!
! !DESCRIPTION:
!  
!  ORIGIN AND CONTACT
!
! !REVISION HISTORY:
! 16 Mar 2023 Manyin     First crack
!EOP
!-----------------------------------------------------------------------

subroutine compute_sza_wrapper_1 ( GC, CLOCK, tdt, label, SZA, RC )

   TYPE(ESMF_GridComp),  intent(inout) :: GC   
   TYPE(ESMF_Clock),     intent(in   ) :: CLOCK
   REAL,                 intent(in   ) :: tdt        ! caller's timestep (sec)
   CHARACTER(LEN=*),     intent(in   ) :: label      ! name of caller
   REAL(KIND=DBL),       intent(  out) :: SZA (:,:)  ! degrees
   INTEGER, OPTIONAL,    intent(  out) :: RC

! LOCAL
   INTEGER :: STATUS
   CHARACTER(LEN=*), PARAMETER :: IAm = 'compute_sza_wrapper_1'

   TYPE(MAPL_MetaComp), POINTER    :: MAPLobj   ! GEOS Generic State    
   REAL, POINTER, DIMENSION(:,:)   :: LONS      ! radians
   REAL, POINTER, DIMENSION(:,:)   :: LATS      ! radians
   TYPE(MAPL_SunOrbit)             :: ORBIT  

   CALL MAPL_GetObjectFromGC ( GC, MAPLobj, __RC__ )
   CALL MAPL_Get ( MAPLobj, LONS=LONS, LATS=LATS, ORBIT=ORBIT, __RC__ )

   CALL compute_sza_main ( LONS=LONS, LATS=LATS, ORBIT=ORBIT, CLOCK=CLOCK, tdt=tdt, label=label, SZA=SZA, __RC__ )

   RETURN_(ESMF_SUCCESS)

end subroutine compute_sza_wrapper_1

      
  end module SZA_from_MAPL_mod
