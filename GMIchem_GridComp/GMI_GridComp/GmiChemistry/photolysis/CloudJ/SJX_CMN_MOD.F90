!---------------------------------------------------------------------------
!BOC
!------------------------------------------------------------------------------
!    'sjx_cmn_mod.f90'  for Cloud-J to move onto Solar-J  v 7.5+ (prather 12/16)
!------------------------------------------------------------------------------
!
!INTERFACE:
!
      MODULE SJX_CMN_MOD
!
      USE FJX_CMN_MOD
!
      implicit none
      public
!
! This is a cmn module to allow Cloud-J to operate but without the full
!   RRTMG bins.
!-----------------------------------------------------------------------
! Turn on/off RRTM-G absorption cross sections
! For this version Cloud-J 7.5 .true. means that bins 19-27 are calculated
!   using the clouds and aerosol data to get heating rates
!   (w/o the RRTM gaseous absorbers)
!
!.sds.. no LW heating needed so turn off LRRTMG...
!      logical, parameter::  LRRTMG= .true.
      logical, parameter::  LRRTMG= .false.
!
!-----------------------------------------------------------------------
!
! W_   = dim = no. of Fast-J Wavelength bins:  currenly only 18, TROP-ONLY is done by zeroing FL fluxes
!      integer, parameter ::  W_=18
! S_   = dim = number of wavelength bins INLCUDING the Solar-J extensions (RRTMG value = 27)
!      integer, parameter ::  S_=27
! L1_  = dim = number of CTM layers + 1
!
!
! Solar J parameters - for Cloud-J v7.5 just do one sub-bin per RRTMg superbin
!-----------------------------------------------------------------------
      integer, parameter:: W_r = S_-W_
      real*8, dimension(L1_,W_r)  ::  TAUG_RRTMG
      real*8, dimension(L1_) :: COLH2O,COLCO2,COLO3,COLCH4,COLDRY,COLO2
!
      integer,  parameter, dimension(S_) ::   NGC =        &
        [  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,          &
           1,  1,  1,  1,  1,  1,  1,  1,  1,  1,          &
           1,  1,  1,  1,  1,  1,  1  ]
!
! actual RRMTMg values
!      integer,  parameter, dimension(S_) ::   NGC =       &
!       [   1,  1,  1,  1,  1,  1,  1,  1,  1,  1,         &
!           1,  1,  1,  1,  1,  1,  1,  5, 10,  2,         &
!          10, 10,  8,  8, 12,  6, 12]
!
!
      END MODULE SJX_CMN_MOD
!EOP
