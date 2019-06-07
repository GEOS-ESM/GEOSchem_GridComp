
!=============================================================================
!
! $Id$
!
! CODE DEVELOPER
!   John Tannahill, LLNL (Original code from Mark Jacobson, UCLA/Stanford)
!   jrt@llnl.gov
!
! FILE
!   smv2chem_par.h
!
! DESCRIPTION
!   This include file contains SmvgearII parameters.
!
!=============================================================================


!     ===========================
#     include "gem_sys_options.h"
#     include "setkin_smv2par.h"
!     ===========================


!     -------------------
!     Integer parameters.
!     -------------------

      integer, parameter ::  &
     &  IGAS   = SK_IGAS,      & ! # of SmvgearII spc
     &  IPHOT  = SK_IPHOT,     & ! # of SmvgearII qjs
     &  ITHERM = SK_ITHERM,    & ! # of SmvgearII qks
     &  NNACT  = SK_NACT     ! # of SmvgearII active spc


!     ---------------------------------------------------------------------
!     Coordinate-System parameters.
!
!     KBLOOP   : max # of grid-points in a vectorized block; should range
!                from  512 (below which vectorization decreases)
!                to   1024 (above which array space is limited)
!     KULOOPIN : kuloop = Min (KULOOPIN, KBLOOP, itloop),
!                  where KULOOP = intended # of grid-cells in a grid-block,
!                        KBLOOP is defined above, and
!                        itloop = ilat * ilong * ivert
!     MXBLOCK  : max # of grid-point blocks
!     ---------------------------------------------------------------------

!     -------------------------------------------------------------
!     Provide a mechanism for setting KBLOOP/KULOOPIN to one set of
!     numbers on parallel machines and another set on non-parallel
!     machines.
!
!     Note that MXBLOCK must be greater than (itloop / KBLOOP) + 1,
!     and that itloop changes with the number of processors and the
!     lon/lat/alt resolution.
!
!     For convenience, and at not much in the way of cost, we just
!     set MXBLOCK to be a number big enough to handle most cases.
!     If it isn't big enough an error message will be generated.
!     If we wanted to be more precise, we could use numbers like
!     the following for parallel cases:
!       For  3 processors, 72x46x25:
!       parameter  (MXBLOCK  = 280)
!       For 30 processors, 72x46x25:
!       parameter  (MXBLOCK  =  30)
!     -------------------------------------------------------------

        integer, parameter :: KBLOOP   =   32 ! Westmere nodes
        integer, parameter :: KULOOPIN =   32
        integer, parameter :: MXBLOCK  = 7500

      integer, parameter ::  &
     &  MAXORD   =  5,    & ! max allowable order of integration method
     &  MBETWEEN = 50   ! max allowable # steps between calls to Pderiv


!     ---------------------
!     Gas-Phase parameters.
!     ---------------------

      integer, parameter ::  &
     &  NCSGAS  =   1,     & ! # of gas     chemistry sets solved for (max=3)
     &  NCSAQ   =   0    ! # of aqueous chamistry sets solved for (max=1)

      integer, parameter ::  &
!c   &  ICS     =   1,   ! # of smvgear eqn sets;
!c                       ! 3 gaschem + 1 aqchem + 1 growth
     &  ICS     = NCSGAS + NCSAQ,    & ! total # of chemistry sets solved for
     &  ICP     = ICS*2  ! # of smvgear chemistry sets x 2 + 1 (for growth)

      integer, parameter ::  &
!c   &  MAXGL   = 140,              ! max # of gains/losses for given array
!c   &  MAXGL   = 350,              ! max # of gains/losses for given array
     &  MAXGL   = 400,                & ! max # of gains/losses for given array
                                    ! MAXGL increased for GMI-T3 {PJC}.
     &  MAXGL2  = (MAXGL + 1) / 2,    & ! 1/2 MAXGL
     &  MAXGL3  = (MAXGL + 4) / 5   ! 1/5 MAXGL

      integer, parameter ::  &
     &  MORDER  =   7,                & ! max order for gear params. for
                                    ! dimension purposes
     &  NMRPROD =  16,                & ! max # of spc in a rxn rate
     &  NMTRATE = ITHERM + IPHOT    ! max # of kinetic + phot rxns


!     -----------------------------------
!     Parameters to minimize array space.
!     -----------------------------------

      integer, parameter ::  &
     &  MXARRAY   = 2100,    & ! max length of matrix put in 1D array
     &  MXGSAER   = IGAS   ! used to dimension and loop on some arrays

      integer, parameter ::  &
     &  MXCOUNT2  = MXGSAER * 100,     & ! array size used to minimize matrix space
     &  MXCOUNT3  = MXGSAER *  32,     & ! array size used to minimize matrix space
!c   &  MXCOUNT4  = MXGSAER *   8    ! array size used to minimize matrix space
     &  MXCOUNT4  = MXGSAER *   9    ! array size used to minimize matrix space
                                     ! MXCOUNT4 increased for GMI-T3 {PJC}.

!     ----------------
!     Real parameters.
!     ----------------

      real*8, parameter ::  &
     &  HMIN  = 1.0d-15   ! min allowable time step (s)

      real*8, parameter ::  &
     &  SMAL1 = 1.0d-06,    & ! small   number
     &  SMAL2 = 1.0d-100  ! smaller number


!     --------------------------------
!     Error tolerance parameters, etc.
!     --------------------------------

      real*8, parameter ::  &
     &  ERRMAXIN  = 1.0d-03,     & ! relative error tolerance
     &  YLOWIN    = 1.0d+04,     & ! absolute error tolerance
     &  YHIIN     = 1.0d+06

      real*8, parameter ::  &
     &  FRACDECIN = 0.25d+00,    & ! fractional cut in smvgear time step
     &  HMAXDAYIN = 8.64d+04,    & ! max time step for day   (s)
     &  HMAXNITIN = 8.64d+04   ! max time step for night (s)

