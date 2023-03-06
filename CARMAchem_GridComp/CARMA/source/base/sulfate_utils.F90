! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

module sulfate_utils
  use carma_precision_mod
  use carma_enums_mod
  use carma_constants_mod
  use carma_types_mod
  use carmastate_mod
  use carma_mod
 
  implicit none
  
  ! Declare the public methods.
  public wtpct_tabaz
  public wtpct_sts
  public sulfate_density
  public sulfate_surf_tens
  
  real(kind=f), public:: dnwtp(46), dnc0(46), dnc1(46)
  real(kind=f), public:: at_hno3(30), at_h2so4(20), amt_hno3(30), amt_h2so4(30)
   
  data dnwtp / 0._f, 1._f, 5._f, 10._f, 20._f, 25._f, 30._f, 35._f, 40._f, &
     41._f, 45._f, 50._f, 53._f, 55._f, 56._f, 60._f, 65._f, 66._f, 70._f, &
     72._f, 73._f, 74._f, 75._f, 76._f, 78._f, 79._f, 80._f, 81._f, 82._f, &
     83._f, 84._f, 85._f, 86._f, 87._f, 88._f, 89._f, 90._f, 91._f, 92._f, &
     93._f, 94._f, 95._f, 96._f, 97._f, 98._f, 100._f /
     
   data dnc0 / 1._f, 1.13185_f, 1.17171_f, 1.22164_f, 1.3219_f, 1.37209_f,       &
     1.42185_f, 1.4705_f, 1.51767_f, 1.52731_f, 1.56584_f, 1.61834_f, 1.65191_f, &
     1.6752_f, 1.68708_f, 1.7356_f, 1.7997_f, 1.81271_f, 1.86696_f, 1.89491_f,   &
     1.9092_f, 1.92395_f, 1.93904_f, 1.95438_f, 1.98574_f, 2.00151_f, 2.01703_f, &
     2.03234_f, 2.04716_f, 2.06082_f, 2.07363_f, 2.08461_f, 2.09386_f, 2.10143_f,&
     2.10764_f, 2.11283_f, 2.11671_f, 2.11938_f, 2.12125_f, 2.1219_f, 2.12723_f, &
     2.12654_f, 2.12621_f, 2.12561_f, 2.12494_f, 2.12093_f /
     
   data dnc1 / 0._f,  -0.000435022_f, -0.000479481_f, -0.000531558_f, -0.000622448_f,&
     -0.000660866_f, -0.000693492_f, -0.000718251_f, -0.000732869_f, -0.000735755_f, &
     -0.000744294_f, -0.000761493_f, -0.000774238_f, -0.00078392_f, -0.000788939_f,  &
     -0.00080946_f, -0.000839848_f, -0.000845825_f, -0.000874337_f, -0.000890074_f,  &
     -0.00089873_f, -0.000908778_f, -0.000920012_f, -0.000932184_f, -0.000959514_f,  &
     -0.000974043_f, -0.000988264_f, -0.00100258_f, -0.00101634_f, -0.00102762_f,    &
     -0.00103757_f, -0.00104337_f, -0.00104563_f, -0.00104458_f, -0.00104144_f,      &
     -0.00103719_f, -0.00103089_f, -0.00102262_f, -0.00101355_f, -0.00100249_f,      &
     -0.00100934_f, -0.000998299_f, -0.000990961_f, -0.000985845_f, -0.000984529_f,  &
     -0.000989315_f /  

  data at_hno3 / &
          -9.3085785070e-1,1.1200784716e-2,-8.7594232370e-5,2.9290261722e-7,-3.6297845637e-10, &
          -9.9276927926e+0,1.3861173987e-1,-7.5302447966e-4,1.9053537417e-6,-1.8847180104e-9,  &
          8.9976745608e-1,-1.1682398549e-2,6.1056862242e-5,-1.5087523503e-7,1.4643716979e-10,  &
          -3.8389447725e-2,4.8922229154e-4,-2.5494288719e-6,6.3306350216e-9,-6.1901374001e-12, &
          7.1911444008e-4,-8.7957856299e-6,4.4035804399e-8,-1.0509519536e-10,9.8591778862e-14, &
          -5.2736179784e-6,6.3762209490e-8,-3.1557072358e-10,7.4508569217e-13,-6.9083781268e-16 /

  data at_h2so4 / &
          -9.8727713620e+1,1.5892180900e+0,-1.0611069051e-2,3.1437317659e-5,-3.5694366687e-8,  &
          2.6972534510e+1,-4.1774114259e-1,2.7534704937e-3,-8.0885350553e-6,9.0919984894e-9,   &
          -3.1506575361e+0,5.1477027299e-2,-3.4697470359e-4,1.0511865215e-6,-1.2167638793e-9,  &
          8.9194643751e-2,-1.4398498884e-3,9.5874823381e-6,-2.8832930837e-8,3.3199717594e-11 /

  data amt_h2so4 / &
          4.9306007769e+0,-2.8124576227e+2,3.6171943540e+4,-7.3921080947e+5,-1.1640936469e+8,  &
          -1.6902946223e+1,5.7843291724e+3,-1.2462848248e+5,3.1325022591e+7,-2.2068275308e+9,  &
          -3.9722280419e+1,1.2350607474e+4,-3.4299494505e+5,6.2642389672e+7,-3.9709694493e+9,  &
          -5.5968384906e+1,1.2922351288e+4,1.3504086346e+6,-1.7890533860e+8,8.8498119334e+9,   &
          -8.2938840352e+1,2.2079294414e+4,2.9469683691e+5,-3.1424855089e+7,1.0884875646e+9,   &
          -1.0647596744e+2,2.7525067463e+4,4.2061852240e+5,-5.1877378665e+7,2.2849838182e+9 /

  data amt_hno3 / &
          2.5757237579e-1,3.4615149493e+3,-1.1460419802e+6,1.6003066569e+8,-8.2005020704e+9,   &
          -2.3081801501e+1,9.7545732474e+3,-1.0751476647e+6,1.2845681641e+8,-5.6387338050e+9,  &
          -1.1454916074e+2,6.7557746435e+4,-1.5833469853e+7,2.0068038322e+9,-9.5789893230e+10, &
          -3.7614906671e+2,2.5721043666e+5,-6.9093891724e+7,8.8886258262e+9,-4.3013954256e+11, &
          -1.1566205559e+3,8.5657451707e+5,-2.4386088798e+8,3.1789555772e+10,-1.5566652191e+12,&
          -2.9858872606e+3,2.2912842052e+6,-6.6825883931e+8,8.7835101682e+10,-4.3330365673e+12 /
contains

  !!  This function calculates the weight % H2SO4 composition of 
  !!  sulfate aerosol, using Tabazadeh et. al. (GRL, 1931, 1997).
  !!  Rated for T=185-260K, activity=0.01-1.0
  !!
  !!  Argument list input:   
  !!    temp = temperature (K)
  !!    h2o_mass = water vapor mass concentration (g/cm3)
  !!    h2o_vp = water eq. vaper pressure (dynes/cm2)
  !!
  !!  Output:
  !!    wtpct_tabaz = weight % H2SO4 in H2O/H2SO4 particle (0-100)
  !!
  !!  Include global constants and variables (BK=Boltzman constant,
  !!   AVG=Avogadro's constant)
  !!
  !! @author Jason English
  !! @ version Apr-2010
  function wtpct_tabaz(carma, temp, h2o_mass, h2o_vp, rc)
   
    real(kind=f)                         :: wtpct_tabaz
    type(carma_type), intent(in)         :: carma     !! the carma object
    real(kind=f), intent(in)             :: temp      !! temperature [K]
    real(kind=f), intent(in)             :: h2o_mass  !! water vapor mass concentration (g/cm3)
    real(kind=f), intent(in)             :: h2o_vp    !! water eq. vaper pressure (dynes/cm2) 
    integer, intent(inout)               :: rc        !! return code, negative indicates failure
      
    !  Declare variables for this routine only
    real(kind=f)     :: atab1,btab1,ctab1,dtab1,atab2,btab2,ctab2,dtab2
    real(kind=f)     :: h2o_num, p_h2o, vp_h2o
    real(kind=f)     :: contl, conth, contt, conwtp
    real(kind=f)     :: activ 
         
    ! Get number density of water (/cm3) from mass concentration (g/cm3)
    h2o_num=h2o_mass*AVG/gwtmol(1)

    !  Get partial pressure of water (dynes/cm2) from concentration (/cm3)
    ! Ideal gas law: P=nkT
    p_h2o=h2o_num*bk*temp

    !  Convert from dynes/cm2 to mb (hPa)
    p_h2o=p_h2o/1000.0_f     ! partial pressure
    vp_h2o=h2o_vp/1000.0_f   ! eq. vp

    !  Activity = water pp in mb / water eq. vp over pure water in mb
    activ = p_h2o/vp_h2o
 
    if (activ.lt.0.05_f) then
      activ = max(activ,1.e-32_f)    ! restrict minimum activity
      atab1 	= 12.37208932_f	
      btab1 	= -0.16125516114_f
      ctab1 	= -30.490657554_f
      dtab1 	= -2.1133114241_f
      atab2 	= 13.455394705_f	
      btab2 	= -0.1921312255_f
      ctab2 	= -34.285174607_f
      dtab2 	= -1.7620073078_f
    elseif (activ.ge.0.05_f.and.activ.le.0.85_f) then
      atab1 	= 11.820654354_f
      btab1 	= -0.20786404244_f
      ctab1 	= -4.807306373_f
      dtab1 	= -5.1727540348_f
      atab2 	= 12.891938068_f	
      btab2 	= -0.23233847708_f
      ctab2 	= -6.4261237757_f
      dtab2 	= -4.9005471319_f
    elseif (activ.gt.0.85_f) then
      activ = min(activ,1._f)      ! restrict maximum activity
      atab1 	= -180.06541028_f
      btab1 	= -0.38601102592_f
      ctab1 	= -93.317846778_f
      dtab1 	= 273.88132245_f
      atab2 	= -176.95814097_f
      btab2 	= -0.36257048154_f
      ctab2 	= -90.469744201_f
      dtab2 	= 267.45509988_f
    else
      if (do_print) write(LUNOPRT,*) 'invalid activity: activity,pp,vp=',activ, p_h2o
      rc = RC_ERROR
      return
    endif

    contl = atab1*(activ**btab1)+ctab1*activ+dtab1
    conth = atab2*(activ**btab2)+ctab2*activ+dtab2
      
    contt = contl + (conth-contl) * ((temp -190._f)/70._f)
    conwtp = (contt*98._f) + 1000._f

    wtpct_tabaz = (100._f*contt*98._f)/conwtp
    wtpct_tabaz = min(max(wtpct_tabaz,1._f),100._f) ! restrict between 1 and 100 %
      
    return
  end function wtpct_tabaz

  !!  This function calculates the weight percent of nitric acid and sulfuric
  !!    acid in STS liquid PSCs based on sulfate mass and mass concentrations
  !!    of water vapor and nitric acid vapor. This routine is based on
  !!    Jacobson et al., 1994 JGR and Tabazadeh et al., 1994 JGR.
  !!  Adapted from WACCM STS equilibrium code by Tabazadeh.
  !!
  !! @ author Parker Case
  !! @ version Apr-2019
  subroutine wtpct_sts( temper, h2so4m, hno3_avail, h2o_avail, press, &
                          wts, wtn, rc)
  !----------------------------------------------------------------------
  !
  ! temper     = Temperature (Kelvin)
  !	h2so4m     = Total mass of H2SO4 (g cm-3)
  !	hno3_avail = HNO3 volume mixing ratio (mol mol-1)
  !	h2o_avail  = H2O volume mixing ratio (mol mol-1)
  !	press      = Total atmospheric pressure (mb)
  ! wts        = Weight percent of sulfuric acid in ternary solution (%)
  ! wtn        = Weight percent of nitric acid in ternary solution (%)
  ! 
  !----------------------------------------------------------------------

    implicit none

  !----------------------------------------------------------------------
  !	... dummy arguments
  !----------------------------------------------------------------------
    real(kind=f), intent(in)  :: h2so4m    
    real(kind=f), intent(in)  :: hno3_avail    
    real(kind=f), intent(in)  :: h2o_avail    
    real(kind=f), intent(in)  :: press
    real(kind=f), intent(in)  :: temper
    real(kind=f), intent(out) :: wts
    real(kind=f), intent(out) :: wtn
    integer, intent(inout)    :: rc        !! return code, negative indicates failure
  !----------------------------------------------------------------------
  !	... local variables
  !----------------------------------------------------------------------
    integer, parameter  :: itermax = 100
    real(kind=f), parameter :: con_lim  = .00005_f
    real(kind=f), parameter :: t0       = 298.15_f
    real(kind=f), parameter :: ks0      = 2.45e6_f
    real(kind=f), parameter :: lower_delx = 1.e-10_f
    real(kind=f), parameter :: upper_delx = .98_f
    real(kind=f), parameter :: con_crit_chem = 5.e-5_f

    logical      :: is_chem
    logical      :: converged
    integer      :: iter, l
    real(kind=f) :: ch2so4
    real(kind=f) :: molh2so4
    real(kind=f) :: molhno3
    real(kind=f) :: wts0
    real(kind=f) :: reduction_factor
    real(kind=f) :: p_h2o
    real(kind=f) :: tr
    real(kind=f) :: wtn0
    real(kind=f) :: pures
    real(kind=f) :: puren
    real(kind=f) :: chno3
    real(kind=f) :: chplus
    real(kind=f) :: cno3
    real(kind=f) :: wrk
    real(kind=f) :: z, num, den
    real(kind=f) :: deltax
    real(kind=f) :: chplusnew
    real(kind=f) :: cno3new
    real(kind=f) :: stren
    real(kind=f) :: sm
    real(kind=f) :: actn
    real(kind=f) :: acts
    real(kind=f) :: nm
    real(kind=f) :: ks
    real(kind=f) :: lnks
    real(kind=f) :: lnks0
    real(kind=f) :: mixyln
    real(kind=f) :: wrk_h2so4
    real(kind=f) :: cphno3new
    real(kind=f) :: con_val
    real(kind=f) :: ti, t1, t2, fi, f1, f2, ymix, hplus, wtotal, ratio 
    real(kind=f) :: con_crit
    real(kind=f) :: h2o_cond
    real(kind=f) :: fratio(0:itermax)
    real(kind=f) :: delx(0:itermax)
    real(kind=f) :: delz(0:itermax)
    real(kind=f) :: a_hno3(5,6)
    real(kind=f) :: a_h2so4(5,4)
    real(kind=f) :: am_hno3(5,6)
    real(kind=f) :: am_h2so4(5,6)
    real(kind=f) :: c_hno3(6)
    real(kind=f) :: c_h2so4(6)
    real(kind=f) :: d_hno3(6)
    real(kind=f) :: d_h2so4(4)
    logical  :: interval_set
    logical  :: positive

    a_hno3 = reshape(at_hno3, shape(a_hno3))
    a_h2so4 = reshape(at_h2so4, shape(a_h2so4))
    am_hno3 = reshape(amt_hno3, shape(am_hno3))
    am_h2so4 = reshape(amt_h2so4, shape(am_h2so4))
    
    converged = .false.

    lnks0 = log( ks0 )
    if( is_chem ) then
       con_crit = con_crit_chem
    else
       con_crit = con_crit_chem
    endif
    p_h2o = h2o_avail * press * .7501_f ! mbar -> torr
  !----------------------------------------------------------------------
  !	Calculating the molality for pure binary systems of H2SO4/H2O
  !	and HNO3/H2O at a given temperature and water vapor pressure
  !	profile (relative humiditiy). Water activities were used to
  !	calculate the molalities as described in Tabazadeh et al. (1994).
  !----------------------------------------------------------------------
    ti  = max( 180._f,temper )
    tr = 1._f/ti
    do l = 1,6
       c_hno3(l) = exp( am_hno3(1,l) + tr*(am_hno3(2,l) + tr*(am_hno3(3,l) + tr*(am_hno3(4,l) + tr*am_hno3(5,l)))) )
       c_h2so4(l) = exp( am_h2so4(1,l) + tr*(am_h2so4(2,l) + tr*(am_h2so4(3,l) + tr*(am_h2so4(4,l) + tr*am_h2so4(5,l)))) )
       ! TODO
    end do
!----------------------------------------------------------------------
!	... H2SO4/H2O pure weight percent and molality
!----------------------------------------------------------------------
    wts0  = max( 0.01_f,c_h2so4(1) + p_h2o*(-1._f*c_h2so4(2) + p_h2o*(c_h2so4(3) + p_h2o*(-1._f*c_h2so4(4) + p_h2o*(c_h2so4(5) - p_h2o*c_h2so4(6))))))
    pures = (wts0 * 1000._f)/(100._f - wts0)
    pures = pures / 98._f
!----------------------------------------------------------------------
!	... HNO3/H2O pure weight percent and molality
!----------------------------------------------------------------------
    puren = max( 0._f,c_hno3(1) + p_h2o*(-1._f*c_hno3(2) + p_h2o*(c_hno3(3) + p_h2o*(-1._f*c_hno3(4) + p_h2o*(c_hno3(5) - p_h2o*c_hno3(6))))) )
!----------------------------------------------------------------------
!	The solving scheme is described both in Jacobson et al. and Tabazadeh
!	et al.. Assumptions:
!	(1) H2SO4 is present only in the aqueous-phase
!	(2) H2SO4 and HNO3 in solution are fully dissocated into H+
!	    SO42- and NO3-
!	(3) PHNO3 + NO3- = constant
!----------------------------------------------------------------------
    ch2so4 = (h2so4m) / 98._f
    if( pures > 0._f ) then
      wrk_h2so4 = (1000._f*ch2so4)/(pures*18._f)
    else
      wrk_h2so4 = 0._f
    endif
    chno3 = 1.2029e-5_f * press * tr * hno3_avail
    do l = 1,6
      d_hno3(l) = a_hno3(1,l) + ti*(a_hno3(2,l) + ti*(a_hno3(3,l) + ti*(a_hno3(4,l) + ti*a_hno3(5,l))))
    end do
    do l = 1,4
      d_h2so4(l) = a_h2so4(1,l) + ti*(a_h2so4(2,l) + ti*(a_h2so4(3,l) + ti*(a_h2so4(4,l) + ti*a_h2so4(5,l))))
    end do
!----------------------------------------------------------------------
!	Note that KS depends only on the temperature
!----------------------------------------------------------------------
    t1       = (ti - t0)/(ti*t0)
    t2       = t0/ti - 1._f - log( t0/ti )
    lnks     = lnks0 - 8792.3984_f * t1  - 16.8439_f * t2
    ks       = exp( lnks )

    converged = .false.
!----------------------------------------------------------------------
!	Setting up initial guesses for the equations above.  Note that
!	for the initial choices the mass and the charge must be conserved.
!----------------------------------------------------------------------
    delx(0)  = .5_f
    z        = .5_f
    delz(0)  = .5_f
    fratio(0) = 0._f
    reduction_factor = .1_f
    interval_set = .false.
Iter_loop : do iter = 1,itermax
!----------------------------------------------------------------------
!	Cwater is the water equation as described in Tabazadeh et
!	al. and Jacobson et al.
!----------------------------------------------------------------------
      cno3new   = chno3 * delx(iter-1)
      cphno3new = chno3 * (1._f - delx(iter-1))
      if( puren > 0._f ) then
        t1 = (1000._f*cno3new)/(puren*18._f)
      else
        t1 = 0._f
      endif
      h2o_cond = t1 + wrk_h2so4
      if( h2o_cond > 0._f ) then
        wrk      = 1.e3_f / (18._f * h2o_cond)
        molhno3  = cno3new * wrk
        molh2so4 = ch2so4 * wrk
      else
        molhno3  = 0._f
        molh2so4 = 0._f
      endif
      stren = molhno3 + 3._f * molh2so4
!----------------------------------------------------------------------
!	(1) Calculate the activity of H2SO4 at a given STREN
!----------------------------------------------------------------------
      sm = stren/3._f
      acts = d_h2so4(1) + sm*(d_h2so4(2) + sm*(d_h2so4(3) + sm*d_h2so4(4)))
!----------------------------------------------------------------------
!	(2) Calculate the activity for HNO3 at a given STREN
!----------------------------------------------------------------------
      nm = stren
      actn = d_hno3(1) + nm*(d_hno3(2) + nm*(d_hno3(3) + nm*(d_hno3(4) + nm*(d_hno3(5) + nm*d_hno3(6)))))
!----------------------------------------------------------------------
!	(3) Calculate the mixed activity coefficient for HNO3 at STREN
!	    as described by Tabazadeh et al.
!----------------------------------------------------------------------
      f1 = 2._f * (molh2so4 + molhno3) * actn
      f2 = 2.25_f * molh2so4 * acts
      if (stren > 0._f) then
        mixyln = (f1 + f2) / (2._f * stren)
      else
        mixyln = 0._f
      endif
      ymix = exp( mixyln )
      hplus = 2._f * molh2so4 + molhno3
      num = ymix**2 * hplus * molhno3
      den = 1000._f * cphno3new * .0820578_f * ti * ks
      if (den .eq. 0._f) PRINT *, 'delx, chno3, cphno3new, ti, ks', delx(iter-1), chno3, cphno3new, ti, ks
      if( chno3 == 0._f ) then
        converged = .true.
        exit Iter_loop
      endif
!----------------------------------------------------------------------
!       the denominator
!       Calculate the ratio F, check convergence
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!       Calculate the ratio F and reset the deltaX (see Jacobson et al.)
!----------------------------------------------------------------------
!       When the numerator is zero, it can drive the denominator
!       to 0, which resulted in a NaN for f and also the fraction
!       ratio. Assume that in this case, the limit of f would
!       really approach 1, not infinity and thus converge the
!       solution.
      if ((num .eq. 0._f) .and. (den .eq. 0._f)) then
        fi = 1._f
      else
        fi = num / den
      endif
      fratio(iter) = abs( fi ) - 1._f
      con_val      = abs( fi - 1._f )
      if( con_val <= con_lim ) then
        converged  = .true.
        exit Iter_loop
      endif
!----------------------------------------------------------------------
!       non-convergence; setup next iterate
!----------------------------------------------------------------------
      if( interval_set ) then
        z = reduction_factor * z
        delz(iter) = z
        if( fi > 1._f ) then
          deltax = -z
        else
          deltax = z
        endif
          delx(iter) = delx(iter-1) + deltax
      else
        if( iter == 1 ) then
          if( fratio(iter) >= 1._f ) then
            positive = .false.
          else
            positive = .true.
          endif
        endif
        if( fratio(iter)*fratio(iter-1) < 0._f ) then
          interval_set = .true.
          reduction_factor = .5_f
          delx(iter) = .5_f*(delx(iter-1) + delx(iter-2))
          z = .5_f*abs( delx(iter-1) - delx(iter-2) )
        else
          if( .not. positive ) then
            delx(iter) = reduction_factor * delx(iter-1)
          else
            delx(iter) = reduction_factor + delx(iter-1)
            if( delx(iter) > upper_delx ) then
              delx(iter) = .5_f
              interval_set = .true.
              reduction_factor = .5_f
            endif
          endif
        endif
      endif
    end do Iter_loop

    wtotal   = molhno3 * 63._f + molh2so4 * 98._f + 1000._f
    wts = (molh2so4 * 9800._f) / wtotal
    wtn = (molhno3 *6300._f)/ wtotal
    if( cno3new /= 0._f .or. cphno3new /= 0._f ) then
      ratio = max( 0._f,min( 1._f,cno3new/(cphno3new + cno3new) ) )
    endif

  end subroutine wtpct_sts
           
  !! Calculates specific gravity (g/cm3) of sulfate of 
  !! different compositions as a linear function of temperature,
  !! based of measurements of H2SO4/H2O solution densities made 
  !! at 0 to 100C tabulated in the International Critical Tables 
  !! (Washburn, ed., NRC, 1928). Measurements have confirmed that 
  !! this data may be linearly extrapolated to stratospheric 
  !! temperatures (180-380K) with excellent accuracy 
  !! (Beyer, Ravishankara, & Lovejoy, JGR, 1996).
  !!
  !! Argument list input:
  !!    wtp = aerosol composition in weight % H2SO4 (0-100)
  !!    temp = temperature in Kelvin
  !!
  !! Output:
  !!    sulfate_density (g/cm3) [function name]
  !!
  !! This function requires setup_sulfate_density to be run
  !! first to read in the density coefficients DNC0 and DNC1
  !! and the tabulated weight percents DNWTP.
  !!
  !! @author Mike Mills
  !! @version Mar-2013
  function sulfate_density(carma, wtp, temp, rc)

  !! Include global constants and variables

    real(kind=f)                         :: sulfate_density
    type(carma_type), intent(in)         :: carma   !! the carma object
    real(kind=f), intent(in)             :: wtp     !! weight percent
    real(kind=f), intent(in)             :: temp    !! temperature 
    integer, intent(inout)               :: rc      !! return code, negative indicates failure
    
    ! Local declarations
    integer           :: i
    real(kind=f)      :: den1, den2
    real(kind=f)      :: frac, temp_loc

    if (wtp .lt. 0.0_f .or. wtp .gt. 100.0_f) then
      if (do_print) write(LUNOPRT,*)'sulfate_density: Illegal value for wtp:',wtp
      rc = RC_ERROR
      return
    endif

    ! limit temperature to bounds of extrapolation
    temp_loc=min(temp, 380.0_f)
    temp_loc=max(temp_loc, 180.0_f)

    i=1

    do while (wtp .gt. dnwtp(i))
     i=i+1
    end do

    den2=dnc0(i)+dnc1(i)*temp_loc

    if (i.eq.1 .or. wtp.eq.dnwtp(i)) then
      sulfate_density=den2
      return
    endif

    den1=dnc0(i-1)+dnc1(i-1)*temp_loc
    frac=(dnwtp(i)-wtp)/(dnwtp(i)-dnwtp(i-1))
    sulfate_density=den1*frac+den2*(1.0_f-frac)

    return
  end function sulfate_density

  !! Calculates the density of ternary HNO3/H2SO4/H2O solutions based on the
  !!   method described by Carslaw et al., 1995, GRL. Valid for temperatures
  !!   between 185K and 240K.
  !! Note that weight percent in this calculation is assumed to be 0-1 not 
  !!   /100.
  !! PAC: check for temperature bounds?
  !!
  !! Argument list input:
  !!    wts  = weight percent of H2SO4
  !!    wtn  = wieght percent of HNO3
  !!    temp = temperature in Kelvin
  !!
  !! Output:
  !!    STS density (g/cm3)
  !!
  !! @author Parker Case
  !! @version Apr-2019
  function sts_density(carma, wts, wtn, temp, rc)

  !! Include global constants and variables

    real(kind=f)                         :: sts_density
    type(carma_type), intent(in)         :: carma   !! the carma object
    real(kind=f), intent(in)             :: wts     !! weight percent H2SO4
    real(kind=f), intent(in)             :: wtn     !! weight percent HNO3
    real(kind=f), intent(in)             :: temp    !! temperature 
    integer, intent(inout)               :: rc      !! return code, negative indicates failure
    
    ! Local declarations
    real(kind=f) :: ms, mn !! molality of sulfuric and nitric acid
    real(kind=f) :: rhos, rhon !! density of binary solutions

    ms = -25000._f * wts / (2453 * (wtn + wts - 1))
    mn = -100000._f * wtn / (6303 * (wtn + wts - 1))
    rhos = 1000+(123.64-5.6e-4*temp**2)*ms-(29.54-1.81E-4*temp**2)*ms**1.5 + &
      (2.343-1.487E-3*temp-1.324E-5*temp**2)*ms**2 
    rhon = 1000+(85.11-5.04e-4*temp**2)*mn-(18.96-1.427e-4*temp**2)*mn**1.5 + &
      (1.458-1.198e-3*temp-9.703e-6*temp**2)*mn**2 
    sts_density = 1/ ((1/rhos) *ms/ (mn+ms)+ (1/rhon) *mn/(mn+ms) ) * 1.e-3

    return
  end function sts_density

  !!  Calculates surface tension (erg/cm2 = dyne/cm) of sulfate of 
  !!  different compositions as a linear function of temperature,
  !!  as described in Mills (Ph.D. Thesis, 1996), derived from
  !!  the measurements of Sabinina and Terpugow (1935).
  !!
  !!  Argument list input:
  !!     WTP = aerosol composition in weight % H2SO4 (0-100)
  !!     TEMP = temperature in Kelvin
  !!
  !!  Output:
  !!     sulfate_surf_tens (erg/cm2) [function name]
  !!
  !!  This function requires setup_sulfate_density to be run
  !!  first to read in the density coefficients DNC0 and DNC1
  !!  and the tabulated weight percents DNWTP.
  !!
  !! @author Mike Mills
  !! @version Mar-2013
  function sulfate_surf_tens(carma, wtp, temp, rc)  
  
    real(kind=f)                         :: sulfate_surf_tens
    type(carma_type), intent(in)         :: carma   !! the carma object
    real(kind=f), intent(in)             :: wtp     !! weight percent
    real(kind=f), intent(in)             :: temp    !! temperature 
    integer, intent(inout)               :: rc      !! return code, negative indicates failure
    
    ! Local declarations
    integer           :: i  
    real(kind=f)      :: sig1, sig2
    real(kind=f)      :: frac, temp_loc
    real(kind=f)      :: stwtp(15), stc0(15), stc1(15)
    
    data stwtp/0._f, 23.8141_f, 38.0279_f, 40.6856_f, 45.335_f, 52.9305_f, 56.2735_f, &
       & 59.8557_f, 66.2364_f, 73.103_f, 79.432_f, 85.9195_f, 91.7444_f, 97.6687_f, 100._f/
    
    data stc0/117.564_f, 103.303_f, 101.796_f, 100.42_f, 98.4993_f, 91.8866_f,     &
       & 88.3033_f, 86.5546_f, 84.471_f, 81.2939_f, 79.3556_f, 75.608_f, 70.0777_f,  &
       & 63.7412_f, 61.4591_f /
  
    data stc1/-0.153641_f, -0.0982007_f, -0.0872379_f, -0.0818509_f,           &
       & -0.0746702_f, -0.0522399_f, -0.0407773_f, -0.0357946_f, -0.0317062_f,   &
       & -0.025825_f, -0.0267212_f, -0.0269204_f, -0.0276187_f, -0.0302094_f,    &
       & -0.0303081_f /

    ! limit temperature to reasonable bounds of extrapolation
    temp_loc=min(temp, 380.0_f)
    temp_loc=max(temp_loc, 180.0_f)
       
    if (wtp .lt. 0.0_f .OR. wtp .gt. 100.0_f) then
      if (do_print) write(LUNOPRT,*)'sulfate_surf_tens: Illegal value for wtp:',wtp
      if (do_print) write(LUNOPRT,*)'sulfate_surf_tens: temp=',temp
      rc = RC_ERROR
      return
    endif
  
    i=1
  
    do while (wtp.gt.stwtp(i))
      i=i+1
    end do
  
    sig2=stc0(i)+stc1(i)*temp_loc
  
    if (i.eq.1 .or. wtp.eq.stwtp(i)) then
      sulfate_surf_tens=sig2
      return
    end if
  
    sig1=stc0(i-1)+stc1(i-1)*temp_loc
    frac=(stwtp(i)-wtp)/(stwtp(i)-stwtp(i-1))
    sulfate_surf_tens=sig1*frac+sig2*(1.0_f-frac)
  
    return
  end function sulfate_surf_tens
            
end module sulfate_utils
