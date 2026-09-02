!-------------------------------------------------------------------------
! NASA GSFC - SSSO Code 610.3
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: GmiSavedVariables_mod
!
! !INTERFACE:
!
      module GmiSavedVariables_mod

      private
      public  :: t_ChemistrySaved

#     include "smv2chem_par.h"

      type t_ChemistrySaved
          ! Variables in setkin_smv2.h that are initialized
          ! with data statements.
          real*8  :: fkoef(NMRPROD, NMTRATE, ICS)
          character (len=14) :: namesp2(0:MXGSAER, ICS)
          integer :: irm      (NMRPROD, NMTRATE, ICS)

          ! Jsparse
          integer, allocatable :: jarraypt(:,:)
          integer, allocatable :: lzero   (:,:)
          integer :: ncsp
    
          ! Ksparse
          integer :: iccount, jccount, kccount, mccount
          integer :: icnt, jcnt, kcnt, mcnt
          integer :: idecomp
          integer :: ijtot, kztot, mztot
          integer :: kbsub, mbsub
          integer :: nqq
          real*8  :: real_nqq
          real*8, allocatable  :: enqq1 (:), enqq2 (:), enqq3 (:)
          real*8, allocatable  :: conp15(:), conpst(:)
          !   coefficients used in selecting step and order (see Ksparse);
          !   pertst2 = original pertst^2
          real*8, allocatable  :: pertst2(:,:)
          !   coefficients for determining order of integration method and
          !   for calculating matrix P
          real*8, allocatable  :: aset(:,:)
    
          ! Do_Smv2_Solver
          real*8, allocatable :: csuma(:)
          real*8, allocatable :: csumc(:)
    
          ! Do_Smv2_Diag
          integer              :: nsteps
          integer              :: nsteps_per_period
          integer, allocatable :: lonloop(:)
          integer, allocatable :: latloop(:)
          integer, allocatable :: altloop(:)
          logical              :: end_period
          integer, allocatable :: isteps (:,:,:)
          real*8,  allocatable :: taccum(:,:,:)
          real*8,  allocatable :: qqjts (:,:,:,:)
          real*8,  allocatable :: qqkts (:,:,:,:)
          real*8,  allocatable :: yts   (:,:,:,:)

          ! Physproc
          integer :: nblockuse_max

          !   smv2chem1.h
          ! ----------------------------------------------------------------
          ! ifreord   : if 1, then reorder grid-cells by stiffness
          ! ih2o      : identifies spc # of water vapor
          ! imgas     : tbd
          ! initrogen : identifies spc # of nitrogen gas
          ! ioxygen   : identifies spc # of oxygen   gas
          ! kuloop    : intended # of grid-cells in a grid-block
          ! lunsmv    : logical unit number to write to when pr_smv2 is true
          ! ncs       : identifies gas chemistry type (1..NCSGAS)
          ! ----------------------------------------------------------------

          integer :: ifreord
          integer :: ih2o
          integer :: imgas
          integer :: initrogen, ioxygen
          integer :: kuloop
          integer :: lunsmv
          integer :: ncs

          ! -----------------------------------------
          ! jphotrat  : tbd
          ! nrates    : # of kinetic rxns (non-photo)
          ! ntloopncs : tbd
          ! ntspec    : # of active + inactive gases
          ! -----------------------------------------

          integer :: jphotrat (ICS)
          integer :: nrates   (ICS)
          integer :: ntloopncs(ICS)
          integer :: ntspec   (ICS)

          ! -----------------------------------------------
          ! inewold   : original spc # of each new jnew spc
          ! npphotrat : tbd
          ! -----------------------------------------------

          integer :: inewold  (MXGSAER, ICS)
          integer :: npphotrat(IPHOT,   ICS)

          ! -------------------------------------------------------------------
          ! fracdec : fraction time step is decreased in Smvgear if convergence
          !           test fails
          ! hmaxnit : max time step for night all chem (s)
          ! -------------------------------------------------------------------

          real*8  :: fracdec
          real*8  :: hmaxnit

          !=============
          !  smv2chem2.h
          !=============
          ! ----------------------------------------------------------------------
          ! ioner    : # of rxns with one active reactant
          !
          ! inorep   : last reordered rxn # prior to sets of two rxns with two
          !            reactants
          ! ithrr    : # of rxns with three active reactants
          ! itwor    : # of rxns with two active reactants
          ! nmair    : # rxns where spc in third position is M = O2 + N2
          ! nmn2     : # rxns where spc in third position is N2
          ! nmo2     : # rxns where spc in third position is O2
          ! nmoth    : # occurences of  spc in third position that are not O2, N2,
          !            or M, or of spc in any position that are inactive (?);
          !            # of occurrences where inactive spc appears in rate eqn (?)
          ! ntrates  : # of kinetic + photo rxns
          !
          ! mappl    : maps original spc #s to spc #s reordered for chemistry
          !
          ! lgasbino : = jold
          !
          ! noldfnew : old rxn rate # corresponding to each reordered rxn
          !
          ! irm2     : new spc # of each active product in each rxn
          ! ----------------------------------------------------------------------

          integer :: ioner  (ICP)
          integer :: nallrat(ICP)

          integer :: inorep (ICS)
          integer :: ithrr  (ICS)
          integer :: itwor  (ICS)
          integer :: nm3bod (ICS)
          integer :: nmair  (ICS)
          integer :: nmn2   (ICS), nmo2 (ICS)
          integer :: nmoth  (ICS)
          integer :: ntrates(ICS)

          integer :: mappl   (MXGSAER, ICS)

          integer :: lgasbino(MAXGL2,  ICS)
          integer :: nreacoth(MAXGL2,  ICS)

          integer :: lgas3bod(MAXGL3,  ICS)
          integer :: losinacp(MAXGL3,  ICS)
          integer :: nreac3b (MAXGL3,  ICS)
          integer :: nreacair(MAXGL3,  ICS)
          integer :: nreacn2 (MAXGL3,  ICS)
          integer :: nreaco2 (MAXGL3,  ICS)

          integer :: jphotnk (NMTRATE, ICS)
          integer :: noldfnew(NMTRATE, ICS)

          integer :: irm2(NMRPROD, NMTRATE, ICS)

          ! ---------------------------------------------------
          ! ischang : # of original nspec spc with >= 1 pd term
          ! ---------------------------------------------------

          integer :: ischang(ICS)

          ! ---------------------------------------------------
          ! kzeroa..e : arrays identifying terms in gloss array
          ! ---------------------------------------------------

          integer :: kzthi(ICP), kztlo(ICP)

          integer :: ikztot(MXCOUNT4)

          integer :: kbh1(MXCOUNT4), kbh2(MXCOUNT4)
          integer :: kbh3(MXCOUNT4), kbh4(MXCOUNT4), kbh5(MXCOUNT4)
          integer :: kbl1(MXCOUNT4), kbl2(MXCOUNT4)
          integer :: kbl3(MXCOUNT4), kbl4(MXCOUNT4), kbl5(MXCOUNT4)

          integer :: mbh1(MXCOUNT4), mbh2(MXCOUNT4)
          integer :: mbh3(MXCOUNT4), mbh4(MXCOUNT4), mbh5(MXCOUNT4)
          integer :: mbl1(MXCOUNT4), mbl2(MXCOUNT4)
          integer :: mbl3(MXCOUNT4), mbl4(MXCOUNT4), mbl5(MXCOUNT4)

          integer :: kzeroa(MXCOUNT4)
          integer :: kzerob(MXCOUNT4), kzeroc(MXCOUNT4)
          integer :: kzerod(MXCOUNT4), kzeroe(MXCOUNT4)

          integer :: mzeroa(MXCOUNT4)
          integer :: mzerob(MXCOUNT4), mzeroc(MXCOUNT4)
          integer :: mzerod(MXCOUNT4), mzeroe(MXCOUNT4)

          integer :: imztot(MXGSAER, ICP)

          ! ------------------------------------------------------------------
          ! jzeroa   : identifies the array position of each jloz1..jhiz1 term
          !
          ! jarrdiag : diagonal term of decompostion
          ! ------------------------------------------------------------------
          
          integer :: ijval (MXCOUNT3)
          integer :: jzeroa(MXCOUNT3)

          integer :: idh1  (MXCOUNT3), idh2  (MXCOUNT3)
          integer :: idh3  (MXCOUNT3), idh4  (MXCOUNT3), idh5(MXCOUNT3)
          integer :: idl1  (MXCOUNT3), idl2  (MXCOUNT3)
          integer :: idl3  (MXCOUNT3), idl4  (MXCOUNT3), idl5(MXCOUNT3)

          integer :: ikdeca(MXCOUNT3), ikdecb(MXCOUNT3)
          integer :: ikdecc(MXCOUNT3), ikdecd(MXCOUNT3)
          integer :: ikdece(MXCOUNT3)

          integer :: kjdeca(MXCOUNT3)
          integer :: kjdecb(MXCOUNT3), kjdecc(MXCOUNT3)
          integer :: kjdecd(MXCOUNT3), kjdece(MXCOUNT3)

          integer :: ijthi   (MXGSAER, ICP), ijtlo(MXGSAER, ICP)
          integer :: jarrdiag(MXGSAER, ICP)
          integer :: jhiz1   (MXGSAER, ICP), jloz1(MXGSAER, ICP)

          ! --------------------------------------------------------------
          ! iarray : length of 1d array holding all sparse matrix points =
          !          sparse matrix dimension (?); = total # of matrix
          !          positions filled after matrix processes (?)
          ! --------------------------------------------------------------

          integer :: iarray(ICP)
          integer :: npdhi (ICP), npdlo(ICP)

          integer :: iialpd  (MXCOUNT2)
          integer :: ipospd  (MXCOUNT2)
          integer :: nkpdterm(MXCOUNT2)

          ! ----------------------------------------------------------------
          ! lossra..e : reaordered rxn rate #s for each loss (and prod) term
          ! ----------------------------------------------------------------

          integer :: nfrhi(ICP), nfrlo(ICP)
          integer :: nplhi(ICP), npllo(ICP)

          integer :: jspcnfr(MXCOUNT4), jspnpl(MXCOUNT4)
          integer :: nknfr  (MXCOUNT4)

          integer :: lossra (MXCOUNT4)
          integer :: lossrb (MXCOUNT4), lossrc(MXCOUNT4)
          integer :: lossrd (MXCOUNT4), lossre(MXCOUNT4)

          integer :: nph1(MXCOUNT4), nph2(MXCOUNT4)
          integer :: nph3(MXCOUNT4), nph4(MXCOUNT4), nph5(MXCOUNT4)
          integer :: npl1(MXCOUNT4), npl2(MXCOUNT4)
          integer :: npl3(MXCOUNT4), npl4(MXCOUNT4), npl5(MXCOUNT4)

          ! ---------------------------------------------------------------
          ! newfold  : new rxn rate # corresponding to each original rate #
          ! ---------------------------------------------------------------

          integer :: nolosp(ICP)
          integer :: newfold(NMTRATE*2, ICS)
          integer :: nknlosp(MAXGL3, ICS)
          integer :: nknphotrt(IPHOT,ICS)

          ! -----------------------------------------------------------------------
          ! abst2    : 1 / (chemintv * chemintv) (s^-2)
          ! errmax   : relative error tolerance; eps should be < 1 for speedy and
          !            reliable results, 10^-3 is reasonable; for many decimal
          !            places of accuracy, decrease eps
          !
          ! abtol    : pre-defined absolute error tolerances; if it is too small,
          !            then integration will take too long; if it is too large,
          !            convergence will be too easy and errors will accumulate, the
          !            time step may be cut too small, and the integration may stop
          !            (delt < hmin or floating point exception in Decomp);
          !            typical gas-phase values of abstol are 10^3 cm^-3 (?),
          !            typical aq -phase values of abstol are
          !            10^-13 to 10^-15 m l^-1 (?)
          ! -----------------------------------------------------------------------

          real*8  :: abst2   (ICS)
          real*8  :: errmax  (ICS)
          real*8  :: hmaxday (ICS)
          real*8  :: timeintv(ICS)
    
          real*8  :: abtol(6, ICS)

          ! -----------------------------------------------------------------------
          ! pertst2  : coefficients used in selecting step and order (see Ksparse);
          !            pertst2 = original pertst^2
          !
          ! aset     : coefficients for determining order of integration method and
          !            for calculating matrix P
          ! -----------------------------------------------------------------------

          !real*8  :: enqq1 (MORDER), enqq2 (MORDER), enqq3 (MORDER)
          !real*8  :: conp15(MORDER), conpst(MORDER)
          !real*8  :: pertst2(MORDER, 3)
          !real*8  :: aset(10, 8)

          ! --------------------------------------------
          ! fracpl  : = -1 for all reactants;
          !           = +1 or +fraction for all products
          ! --------------------------------------------

          real*8  :: fracpl (MXCOUNT2)
          real*8  :: fracnfr(MXCOUNT4)

      end type t_ChemistrySaved

      end module GmiSavedVariables_mod
