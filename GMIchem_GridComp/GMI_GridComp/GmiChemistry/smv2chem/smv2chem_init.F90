

!=============================================================================
!
! $Id$
!
! CODE DEVELOPER
!   Original code from Mark Z. Jacobson ((C) COPYRIGHT, 1993).
!   LLNL modifications:  John Tannahill
!                        jrt@llnl.gov
!
! FILE
!   smv2chem_init.F
!
! ROUTINES
!   Do_Smv2_Init
!
!=============================================================================


!-----------------------------------------------------------------------------
!
! ROUTINE
!   Do_Smv2_Init
!
! DESCRIPTION
!   This is the initialization routine for the ordinary differential
!   equation solver, "Smvgear II" (Sparse Matrix Vectorized Gear code).
!
! ARGUMENTS
!   smv_filnam   : filename for smv2 diagnostics
!   do_smv_reord : reorder the grid-boxes in order of stiffness?
!   pr_diag      : print some diagnostic output to screen?
!   pr_smv2      : should the SmvgearII output file be written
!                  (non-parallel mode only)?
!   loc_proc     : local processor #
!   itloop       : # of zones (ilong * ilat * ivert)
!   chemintv     : chemistry time step (s)
!
!-----------------------------------------------------------------------------

      subroutine Do_Smv2_Init (savedVars, smv_filnam, do_smv_reord, &
                           pr_diag, pr_smv2, loc_proc, itloop, chemintv)

      use GmiPrintError_mod, only : GmiPrintError
      use GmiASCIIoperations_mod, only : AsciiOpenWrite
      use GmiSavedVariables_mod, only : t_ChemistrySaved

      IMPLICIT NONE

#     include "smv2chem_par.h"

!#     include "smv2chem1.h"
!#     include "smv2chem2.h"

!     ----------------------
!     Argument declarations.
!     ----------------------

      character (len=*), intent(in) :: smv_filnam
      logical, intent(in) :: do_smv_reord
      logical, intent(in) :: pr_diag
      logical, intent(in) :: pr_smv2
      integer, intent(in) :: loc_proc
      integer, intent(in) :: itloop
      real*8,  intent(in) :: chemintv

      type(t_ChemistrySaved), intent(inOut) :: savedVars

!     ----------------------
!     Variable declarations.
!     ----------------------

      character (len=75) :: err_msg
      integer :: ii, nn
      real*8  :: abhi
      real*8  :: ablo
      real*8  :: mim6


!     ----------------
!     Begin execution.
!     ----------------

      if (pr_diag) then
        Write (6,*) 'Do_Smv2_Init called by ', loc_proc
      end if

!     --------------------------------------
!     Check to see if MXBLOCK is big enough.
!     --------------------------------------

      if (MXBLOCK < ((itloop / KBLOOP) + 1)) then
        err_msg = 'MXBLOCK too small in Do_Smv2_Init.'
        call GmiPrintError  &
     &    (err_msg, .false., 1, MXBLOCK, 0, 0, 0.0d0, 0.0d0)
        Write (6,*) 'Needs to be at least:  ', (itloop / KBLOOP) + 1
        call GmiPrintError  &
     &    (err_msg, .true., 2, itloop, KBLOOP, 0, 0.0d0, 0.0d0)
      end if

      if (pr_smv2) then
        call AsciiOpenWrite (savedVars%lunsmv, smv_filnam)
      end if

!     ---------------------------------
!     Determine appropriate block size.
!     ---------------------------------

      savedVars%kuloop = Min (KULOOPIN, KBLOOP, itloop)

!     --------------------------------
!     Chemistry solver-related values.
!     --------------------------------

      if (do_smv_reord) then
        savedVars%ifreord = 1
      else
        savedVars%ifreord = 0
      end if


      call initializeSavedVars(savedVars)


      savedVars%fracdec = FRACDECIN
      savedVars%hmaxnit = HMAXNITIN
      savedVars%ncs     = NCSGAS

      savedVars%errmax   (savedVars%ncs) = ERRMAXIN
      savedVars%hmaxday  (savedVars%ncs) = HMAXDAYIN
!c    ylow     (ncs) = YLOWIN

      savedVars%ntloopncs(savedVars%ncs) = itloop
      savedVars%timeintv (savedVars%ncs) = chemintv

      savedVars%abst2    (savedVars%ncs) = 1.0d0 / (chemintv * chemintv)

      savedVars%abtol(1,savedVars%ncs) = YHIIN
      savedVars%abtol(6,savedVars%ncs) = YLOWIN

      do nn = 1, NCSGAS
         abhi = Log10 (savedVars%abtol(1,nn))
         ablo = Log10 (savedVars%abtol(6,nn))

         do ii = 2, 5
            mim6 = 6 - ii
            savedVars%abtol(ii,nn) = &
                     10.0d0**(ablo + ((abhi - ablo) * mim6 / 5.0d0))
         end do
      end do

      call initJsparse(savedVars)
      call initKsparse(savedVars)
      call initSMV2solver(savedVars, itloop)
      call initPhysproc(savedVars)

      call JsparseGMI(savedVars, pr_smv2, savedVars%lunsmv, savedVars%ncs,&
            savedVars%jphotrat, savedVars%nrates,   &
            savedVars%ntspec,   savedVars%inewold,  savedVars%npphotrat,&
            savedVars%ncsp,     savedVars%lzero,    savedVars%jarraypt, &
            savedVars%iccount,  savedVars%jccount,  savedVars%kccount,  &
            savedVars%mccount,  savedVars%icnt,     savedVars%jcnt,     &
            savedVars%kcnt,     savedVars%mcnt,     savedVars%idecomp,  &
            savedVars%ijtot,    savedVars%kztot,    savedVars%mztot,    &
            savedVars%kbsub,    savedVars%mbsub, savedVars%fkoef, &
            savedVars%irm, savedVars%namesp2)

      return

      end subroutine Do_Smv2_Init
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
      subroutine initJsparse (savedVars)

      use GmiSavedVariables_mod, only : t_ChemistrySaved

      IMPLICIT NONE

#     include "smv2chem_par.h"

      type(t_ChemistrySaved), intent(inOut) :: savedVars
      INTEGER :: STATUS
      !character(len=ESMF_MAXSTR) :: IAm = "initJsparse"
!     ----------------------------------------------------------------------
!     jarraypt : identifies 1d array point for each 2d point i,j
!     lzero    : = 1 if an array spot is filled with a non-zero value; it is
!                updated as we simulate the order of calculations during a
!                practice l-u decomposition (?)
!     ----------------------------------------------------------------------
       allocate(savedVars%jarraypt(MXGSAER,MXGSAER), STAT=STATUS)
       savedVars%jarraypt = 0
       allocate(savedVars%lzero(MXGSAER,MXGSAER), STAT=STATUS)
       savedVars%lzero    = 0
!      ncsp  ! NCSGAS       => for daytime   gas chemistry
!            ! NCSGAS + ICS => for nighttime gas chemistry
       savedVars%ncsp = NCSGAS + ICS
      end subroutine initJsparse
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
    subroutine initKsparse(savedVars)

      use GmiSavedVariables_mod, only : t_ChemistrySaved

      IMPLICIT NONE

#     include "smv2chem_par.h"

!# include "smv2chem2.h" ! aset,enqq1,enqq2,enqq3,conp15,conpst,pertst2 
      type(t_ChemistrySaved), intent(inOut) :: savedVars

      integer :: i, STATUS
      integer :: nqq ! order of integration method; varies between 1 and MAXORD
      real*8  :: real_nqq
      !character(len=ESMF_MAXSTR) :: IAm = "initKsparse"
!     ---------------------------------------------------------------
!     pertst : coefficients used to select step-size and order; thus,
!              only about 1% accuracy needed; see Gear (1971) or
!              Hindmarsh (1973, UCID-30059).
!     ---------------------------------------------------------------
      real*8 :: pertst(MORDER,3)
      real*8, dimension(MORDER*3), parameter :: &
        pertstData = (/  &
        2.0d0,       4.5d0,       7.333d0,  10.42d0,  &
       13.7d0,      17.15d0,      1.0d0,  &
        3.0d0,       6.0d0,       9.167d0,  12.5d0,  &
       15.98d0,      1.0d0,       1.0d0,  &
        1.0d0,       1.0d0,       0.5d0,     0.1667d0,  &
        0.04133d0,   0.008267d0,  1.0d0 /)

      pertst = reshape( pertstData, (/ MORDER, 3 /) )

      allocate( savedVars%enqq1 (MORDER), savedVars%enqq2 (MORDER), savedVars%enqq3 (MORDER), STAT=STATUS)
      savedVars%enqq1 = 0.0d0
      savedVars%enqq2 = 0.0d0
      savedVars%enqq3 = 0.0d0

      allocate( savedVars%conp15(MORDER), savedVars%conpst(MORDER), STAT=STATUS)
      savedVars%conp15 = 0.0d0
      savedVars%conpst = 0.0d0

      allocate( savedVars%pertst2(MORDER, 3), STAT=STATUS)
      savedVars%pertst2 = 0.0d0

      allocate( savedVars%aset(10, 8), STAT=STATUS)
      savedVars%aset = 0.0d0

      savedVars%icnt    = 0
      savedVars%jcnt    = 0
      savedVars%kcnt    = 0
      savedVars%mcnt    = 0

      savedVars%iccount = 0
      savedVars%jccount = 0
      savedVars%kccount = 0
      savedVars%mccount = 0

      savedVars%idecomp = 0

      savedVars%ijtot   = 0
      savedVars%kztot   = 0
      savedVars%mztot   = 0

      savedVars%kbsub   = 0
      savedVars%mbsub   = 0

      do i = 1, 6
        savedVars%aset(i,2) = 1.0d0
        savedVars%aset(i,8) = 0.0d0
      end do

      savedVars%aset(1,1) = 1.0d0

      savedVars%aset(2,1) = 2.0d0    /    3.0d0
      savedVars%aset(2,3) = 1.0d0    /    3.0d0

      savedVars%aset(3,1) = 6.0d0    /   11.0d0
      savedVars%aset(3,3) = 6.0d0    /   11.0d0
      savedVars%aset(3,4) = 1.0d0    /   11.0d0

      savedVars%aset(4,1) = 12.0d0   /   25.0d0
      savedVars%aset(4,3) =   .70d0
      savedVars%aset(4,4) =   .20d0
      savedVars%aset(4,5) =   .020d0

      savedVars%aset(5,1) =   60.0d0 /  137.0d0
      savedVars%aset(5,3) =  225.0d0 /  274.0d0
      savedVars%aset(5,4) =   85.0d0 /  274.0d0
      savedVars%aset(5,5) =   15.0d0 /  274.0d0
      savedVars%aset(5,6) =    1.0d0 /  274.0d0

      savedVars%aset(6,1) =  180.0d0 /  441.0d0
      savedVars%aset(6,3) =  406.0d0 /  441.0d0
      savedVars%aset(6,4) =  735.0d0 / 1764.0d0
      savedVars%aset(6,5) =  175.0d0 / 1764.0d0
      savedVars%aset(6,6) =   21.0d0 / 1764.0d0
      savedVars%aset(6,7) =    1.0d0 / 1764.0d0

      do nqq = 1, 7
        real_nqq       = nqq
        savedVars%enqq1(nqq)     = 0.5d0 / real_nqq
        real_nqq       = nqq + 1
        savedVars%enqq2(nqq)     = 0.5d0 / real_nqq
        real_nqq       = nqq + 2
        savedVars%enqq3(nqq)     = 0.5d0 / real_nqq
        savedVars%conpst(nqq)    = 1.0d0 / (pertst(nqq,1) * savedVars%enqq3(nqq))
        savedVars%conp15(nqq)    = 1.5d0 * savedVars%conpst(nqq)
        savedVars%pertst2(nqq,1) = pertst(nqq,1) * pertst(nqq,1)
        savedVars%pertst2(nqq,2) = pertst(nqq,2) * pertst(nqq,2)
        savedVars%pertst2(nqq,3) = pertst(nqq,3) * pertst(nqq,3)
      end do

    end subroutine initKsparse
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
    subroutine initSMV2solver(savedVars, itloop)

      use GmiSavedVariables_mod, only : t_ChemistrySaved

      IMPLICIT NONE

#     include "smv2chem_par.h"

      type(t_ChemistrySaved), intent(inOut) :: savedVars

      integer :: itloop
      INTEGER :: STATUS
      real*8  :: fkoef(NMRPROD, NMTRATE, ICS)
      integer :: irm      (NMRPROD, NMTRATE, ICS)
      character (len=14) :: namesp2(0:MXGSAER, ICS)
      integer :: nm3bod (ICS)
      integer :: nmair  (ICS)
      integer :: nmn2   (ICS)
      integer :: nmo2   (ICS)

      integer :: lgas3bod(MAXGL3,  ICS)
      integer :: nreac3b (MAXGL3,  ICS)
      integer :: nreacair(MAXGL3,  ICS)
      integer :: nreacn2 (MAXGL3,  ICS)
      integer :: nreaco2 (MAXGL3,  ICS)

      integer :: IH2O, IMGAS, IOXYGEN, INITROGEN
      integer :: i, j, jgas

      Allocate (savedVars%csuma(itloop), STAT=STATUS)
      Allocate (savedVars%csumc(itloop), STAT=STATUS)
      savedVars%csuma = 0.0d0; savedVars%csumc = 0.0d0

#     include "setkin_smv2.h"
      ! Coming from setkin_smv2.h
      savedVars%lgas3bod(:,:)    = LGAS3BOD(:,:)
      savedVars%nreac3b (:,:)    = NREAC3B (:,:)
      savedVars%nreacair(:,:)    = NREACAIR(:,:)
      savedVars%nreaco2 (1:4, 1) = NREACO2 (1:4, 1)
      savedVars%nreacn2 (1,1)    = NREACN2(1,1)

      savedVars%nmn2  (1) = NMN2(1)
      savedVars%nmo2  (1) = NMO2(1)
      savedVars%nmair (1) = NMAIR(1)
      savedVars%nm3bod(1) = NM3BOD(1)

      savedVars%namesp2 = namesp2
      savedVars%fkoef   = fkoef
      savedVars%irm     = irm  

      savedVars%IH2O      = IH2O
      savedVars%IMGAS     = IMGAS
      savedVars%IOXYGEN   = IOXYGEN
      savedVars%INITROGEN = INITROGEN

    end subroutine initSMV2solver
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
    subroutine initPhysproc(savedVars)

      use GmiSavedVariables_mod, only : t_ChemistrySaved

      IMPLICIT NONE

      type(t_ChemistrySaved), intent(inOut) :: savedVars
      savedVars%nblockuse_max = -1
    end subroutine initPhysproc
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
    subroutine initSMV2diag(savedVars, tdt, i1, i2, ju1, j2, k1, k2, &
                            ilong, ilat, ivert, itloop, &
                            num_qjs, num_qks, num_active, pr_nc_period)

      use GmiSavedVariables_mod, only : t_ChemistrySaved

      IMPLICIT NONE

#     include "smv2chem_par.h"

      type(t_ChemistrySaved), intent(inOut) :: savedVars
      integer, intent(in) :: ilong, ilat, ivert, itloop
      integer, intent(in) :: i1, i2, ju1, j2, k1, k2
      integer, intent(in) :: num_qjs, num_qks, num_active
      real*8 , intent(in) :: tdt 
      integer, intent(in) :: pr_nc_period !! need to properly initialize
      INTEGER :: STATUS

      integer :: ik, ij, il, ind
      integer ::  ic
      !character(len=ESMF_MAXSTR) :: IAm = "initSMV2diag"
!
      Allocate (savedVars%lonloop(itloop), STAT=STATUS)
      !VERIFY_(STATUS)
      Allocate (savedVars%latloop(itloop), STAT=STATUS)
      !VERIFY_(STATUS)
      Allocate (savedVars%altloop(itloop), STAT=STATUS)
      !VERIFY_(STATUS)
      savedVars%lonloop = 0
      savedVars%latloop = 0
      savedVars%altloop = 0

      Allocate (savedVars%isteps(ilong, ilat, ivert), STAT=STATUS)
      !VERIFY_(STATUS)
      savedVars%isteps = 0

      Allocate (savedVars%taccum(ilong, ilat, ivert), STAT=STATUS)
      !VERIFY_(STATUS)
      savedVars%taccum = 0.0d0

      Allocate (savedVars%qqjts(ilong, ilat, ivert, num_qjs), STAT=STATUS)
      !VERIFY_(STATUS)
      Allocate (savedVars%qqkts(ilong, ilat, ivert, num_qks), STAT=STATUS)
      !VERIFY_(STATUS)
      Allocate (savedVars%yts  (ilong, ilat, ivert, num_active), STAT=STATUS)
      !VERIFY_(STATUS)
      savedVars%qqjts = 0.0d0
      savedVars%qqkts = 0.0d0
      savedVars%yts = 0.0d0

      ind = 0

      do ik = k1, k2
         do ij = ju1, j2
            do il = i1, i2
               ind = ind + 1
               savedVars%lonloop(ind) = il
               savedVars%latloop(ind) = ij
               savedVars%altloop(ind) = ik
            end do
         end do
      end do

      savedVars%nsteps = 0
      savedVars%nsteps_per_period = Nint (pr_nc_period / tdt)

      savedVars%end_period = .false.

      savedVars%isteps (:,:,:) = 0
      savedVars%taccum (:,:,:) = 0.0d0
      savedVars%qqjts(:,:,:,:) = 0.0d0
      savedVars%qqkts(:,:,:,:) = 0.0d0
      savedVars%yts  (:,:,:,:) = 0.0d0

    end subroutine initSMV2diag
!------------------------------------------------------------------------------
!-------------------------------------------------------------------------
    subroutine finalizeJsparse (savedVars)

      use GmiSavedVariables_mod, only : t_ChemistrySaved

      IMPLICIT NONE

      type(t_ChemistrySaved), intent(inOut) :: savedVars
       deallocate(savedVars%jarraypt)
       deallocate(savedVars%lzero)
    end subroutine finalizeJsparse
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
    subroutine finalizeKsparse (savedVars)

      use GmiSavedVariables_mod, only : t_ChemistrySaved

      IMPLICIT NONE

      type(t_ChemistrySaved), intent(inOut) :: savedVars
      deallocate(savedVars%enqq1, savedVars%enqq2, savedVars%enqq3)
      deallocate(savedVars%conp15, savedVars%conpst)
      deallocate(savedVars%pertst2)
      deallocate(savedVars%aset)
    end subroutine finalizeKsparse
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
    subroutine finalizeSMV2solver (savedVars)

      use GmiSavedVariables_mod, only : t_ChemistrySaved

      IMPLICIT NONE

      type(t_ChemistrySaved), intent(inOut) :: savedVars
      deallocate(savedVars%csuma, savedVars%csumc)
    end subroutine finalizeSMV2solver
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
    subroutine finalizePhysproc(savedVars)

      use GmiSavedVariables_mod, only : t_ChemistrySaved

      IMPLICIT NONE

      type(t_ChemistrySaved), intent(inOut) :: savedVars
      return
    end subroutine finalizePhysproc
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
    subroutine finalizeSMV2diag (savedVars)

      use GmiSavedVariables_mod, only : t_ChemistrySaved

      IMPLICIT NONE

      type(t_ChemistrySaved), intent(inOut) :: savedVars
      deallocate(savedVars%lonloop)
      deallocate(savedVars%latloop)
      deallocate(savedVars%altloop)
      deallocate(savedVars%isteps)
      deallocate(savedVars%taccum)
      deallocate(savedVars%qqjts)
      deallocate(savedVars%qqkts)
      deallocate(savedVars%yts)
    end subroutine finalizeSMV2diag
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
      subroutine initializeSavedVars(savedVars)

      use GmiSavedVariables_mod, only : t_ChemistrySaved

      IMPLICIT NONE

      type(t_ChemistrySaved), intent(inOut) :: savedVars

      integer, parameter :: initVal = 0

      ! Begin initialize integer arrays to a default value
      savedVars%inewold   = initVal
      savedVars%npphotrat = initVal
      savedVars%mappl     = initVal
      savedVars%lgasbino  = initVal
      savedVars%nreacoth  = initVal
      savedVars%lgas3bod  = initVal
      savedVars%losinacp  = initVal
      savedVars%nreac3b   = initVal
      savedVars%nreacair  = initVal
      savedVars%nreacn2   = initVal
      savedVars%nreaco2   = initVal
      savedVars%jphotnk   = initVal
      savedVars%noldfnew  = initVal
      savedVars%irm2      = initVal

      savedVars%ikztot    = initVal

      savedVars%kbh1      = initVal
      savedVars%kbh2      = initVal
      savedVars%kbh3      = initVal
      savedVars%kbh4      = initVal
      savedVars%kbh5      = initVal
      savedVars%kbl1      = initVal
      savedVars%kbl2      = initVal
      savedVars%kbl3      = initVal
      savedVars%kbl4      = initVal
      savedVars%kbl5      = initVal

      savedVars%mbh1      = initVal
      savedVars%mbh2      = initVal
      savedVars%mbh3      = initVal
      savedVars%mbh4      = initVal
      savedVars%mbh5      = initVal
      savedVars%mbl1      = initVal
      savedVars%mbl2      = initVal
      savedVars%mbl3      = initVal
      savedVars%mbl4      = initVal
      savedVars%mbl5      = initVal

      savedVars%idh1      = initVal
      savedVars%idh2      = initVal
      savedVars%idh3      = initVal
      savedVars%idh4      = initVal
      savedVars%idh5      = initVal
      savedVars%idl1      = initVal
      savedVars%idl2      = initVal
      savedVars%idl3      = initVal
      savedVars%idl4      = initVal
      savedVars%idl5      = initVal

      savedVars%nph1      = initVal
      savedVars%nph2      = initVal
      savedVars%nph3      = initVal
      savedVars%nph4      = initVal
      savedVars%nph5      = initVal
      savedVars%npl1      = initVal
      savedVars%npl2      = initVal
      savedVars%npl3      = initVal
      savedVars%npl4      = initVal
      savedVars%npl5      = initVal

      savedVars%kzeroa    = initVal
      savedVars%kzerob    = initVal
      savedVars%kzeroc    = initVal
      savedVars%kzerod    = initVal
      savedVars%kzeroe    = initVal

      savedVars%mzeroa    = initVal
      savedVars%mzerob    = initVal
      savedVars%mzeroc    = initVal
      savedVars%mzerod    = initVal
      savedVars%mzeroe    = initVal

      savedVars%imztot    = initVal

      savedVars%ijval     = initVal
      savedVars%jzeroa    = initVal

      savedVars%ikdeca = initVal
      savedVars%ikdecb = initVal
      savedVars%ikdecc = initVal
      savedVars%ikdecd = initVal
      savedVars%ikdece = initVal

      savedVars%kjdeca = initVal
      savedVars%kjdecb = initVal
      savedVars%kjdecc = initVal
      savedVars%kjdecd = initVal
      savedVars%kjdece = initVal

      savedVars%ijthi    = initVal
      savedVars%ijtlo    = initVal
      savedVars%jarrdiag = initVal
      savedVars%jhiz1    = initVal
      savedVars%jloz1    = initVal
      savedVars%iialpd   = initVal
      savedVars%ipospd   = initVal
      savedVars%nkpdterm = initVal
      savedVars%jspcnfr  = initVal

      savedVars%jspnpl = initVal
      savedVars%nknfr  = initVal

      savedVars%lossra = initVal
      savedVars%lossrb = initVal
      savedVars%lossrc = initVal
      savedVars%lossrd = initVal
      savedVars%lossre = initVal

      savedVars%newfold   = initVal
      savedVars%nknlosp   = initVal
      savedVars%nknphotrt = initVal

      savedVars%abtol   = 0.0d0
      savedVars%fracpl  = 0.0d0
      savedVars%fracnfr = 0.0d0
      ! End initialize integer arrays to a default value

      RETURN

      end subroutine initializeSavedVars

!------------------------------------------------------------------------------
