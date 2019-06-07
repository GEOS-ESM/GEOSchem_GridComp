!=============================================================================
!
! $Id$
!
! CODE DEVELOPER
!   originated from GRANTOUR programmed by John Walton
!   modified by Xiaohong Liu
!
! FILE
!   so2topar_rates.F
!
! ROUTINES
!   Do_So2toPar_Rates
!
!=============================================================================

!-----------------------------------------------------------------------------
!
! ROUTINE
!   Do_So2toPar_Rates (chemst60)
!
! DESCRIPTION
!
!     ******************************************************************
!     * This set of routines will compute the vapor removal rate
!     * coefficient for the transfer of vapor to the surface of dust
!     * and sea salt particles.
!     ******************************************************************
!
! **********************************************************************
! CALCULATION OF THE VAPOR REMOVAL RATE COEFFICIENT.             5/14/97
!
!
! *----------
! DEFINITIONS
! *----------
!
! alph     The fraction of the molecules colliding with an
!          aerosol particle that are absorbed by it.                  ()
!          (= 0.1, or Joyce's simple funtion)
!
! ctr(Kn)  The Fuchs and Sutugin interpolation term accounting
!          for the transition between particle size ranges
!          described by continuum equations and gas-kinetic
!          equations.                                                 ()
!
!          ctr = [ 0.71 + (4/3) Kn ] / [ 1 + Kn ]
!
! dmol     The molecular diffusion coefficient for molecules
!          of vapor (v) in air(a).                               (cm2/s)
!
!          dmol = 0.1 * (1.0/p) * (T/298) ** 3/2
!
! gc       The gas constant.
!          (= 8.3169e7 [(cm3 dynes/cm2) / (mole K)])
!
! G(r)     Fuchs & Sutugin particle size dependent mass transfer
!          rate coefficient.                                     (cm3/s)
!
!          G = (4 pi r dmol)
!            / ( 1 + Kn [ ctr + (4/3) ( 1/alph - 1 ) ] ).
!
! Kn       The Knudsen number.                                        ()
!
!          Kn = lambda / r
!             = (3 dmol) / (r cbar)
!
! kv       The rate coefficient for vapor -> particle removal      (s-1)
!          for each type of particle, i.e. 2 dust particle size
!          ranges and 2 size ranges for sea salt.
!
! lambda   The mean free path of vapor molecules in air.
!          This is just one of the possible definitions.  However,
!          it is required for the function G to have the proper
!          limit for large Knudsen number.
!
!          lambda = (3 dmol) / cbar
!
! n        normalized number conentration, i.e.
!
!          dndrt = dn/dlogr
!          integral[ (dn/dlogr) dlogr ] = 1
!
! Mt       particle mass concentration                          (gm/cm3)
!
! N        number concentration of the aerosol, as a function
!          of r, described by the number of particles with        (cm-3)
!          radii between r and r+dr
!
!          (dN/dlogr)     or
!          Nt (dn/dlogr)
!
! Nt       the total number of particles                          (cm-3)
!
!          = integral[ (dN/dlogr) dlogr ]
!
! p        Air pressure.                                          (bars)
!
! r        Aerosol particle radius.
!
! rhop     the density of an aerosol particle                   (gm/cm3)
!
! T        Air temperature.                                          (K)
!
! cbar     The mean velocity of vapor molecules in air.           (cm/s)
!
!          cbar = [ (8 gc T) / ( pi wtm(3) ) ] ** 1/2
!
!
! *------------------------------------
! THE PARTICULATE DISTRIBUTION FUNCTION
! *------------------------------------
!
! For each type of particulate there is a distribution function
! that describes the fraction of the total number of particles that
! is of a given size.  Regardless of the total number of particles,
! this distribution function will be the same.  I will designate it
! as n(r) (dimensionless), and it is defined by its logarithmic
! (base 10) derivative
!
!        dn/dlogr.                                                   (1)
!
! By its definition, the integral of dn/dlogr over the full particle
! size range is unity
!
!        integral[ (dn/dlogr) dlogr ] = 1.                           (2)
!
! For a specific observation of total particle number, Nt (cm-3),
! N(r) can be defined, with
!
!        dN/dlogr = Nt (dn/dlogr).                                   (3)
!
! Since in the final relations a ratio of integrals involving
! dn/dlogr will appear, the normalization (2) is not required.
! I have included it to preserve some mathematical neatness and
! to avoid some confusion for me.
!
!
! Given the function n(r), the particle mass concentration,
! Mt (gm/cm3), (at GRANTOUR nodes) and the particle density,
! rhop (gm/cm3), the total number of particles, Nt (cm-3), can
! be calculated.
!
! The mass of all particles with radii between r and r+dr is
!
!        dM(r) = rhop [ (4/3) pi r**3 ] (dN/dlogr) dlogr
!              = Nt rhop [ (4/3) pi r**3 ] (dn/dlogr) dlogr.         (4)
!
! Integrating (4) over the range of r gives the total mass
!
!        Mt = Nt rhop integral[ (4/3) pi r**3 (dn/dlogr) dlogr ].    (5)
!
! All quantities, except Nt, being given, this expression
! can be solved for total number
!
!    Nt = Mt / ( rhop integral[ (4/3) pi r**3 (dn/dlogr) dlogr ] ).  (6)
!
! NOTE:  Nt is always obtained using the "dry radius" distribution
!        function, dn/dlogr.  However, if the aerosol includes
!        water, i.e. wet & dry diameters different, the integration
!        in (6) must use the "wet radius" distribution function.
!        The same number of particles are present, but their size
!        distribution will be different.
!
!
! *---------------------------------
! THE VAPOR REMOVAL RATE COEFFICIENT
! *---------------------------------
!
! The rate coefficient, kv (s-1), for the removal of the gas
! phase of a species to an aerosol surface is
!
!        kv = integral[ G(r,T,p) (dN/dlogr) dlogr ]
!           = Nt integral[ G(r,T,p) (dn/dlogr) dlogr ],              (7)
!
! where G (cm3/s) is the Fuchs & Sutugin expression for the mass
! transfer rate coefficient
!
!        G = (4 pi r dmol)
!          / ( 1 + Kn [ ctr + (4/3) ( 1/alph - 1 ) ] ),              (8)
!
! and ctr is an interpolation factor that accounts for the
! transition between particle size ranges in which continuum
! and gas-kinetic equations are applicable
!
!        ctr = [ 0.71 + (4/3) Kn ] / [ 1 + Kn ].                     (9)
!
! The dependence of G upon temperature and pressure is through its
! relation to dmol and cbar which in turn depend upon T and p.
!
!
! Substituting (6) for Nt in (7) we get the removal rate
! coefficient
!
!    kv = Mt * integral[ G(r,T,p) (dn/dlogr) dlogr ]
!       / ( rhop * integral[ (4/3) pi r**3 (dn/dlogr) dlogr ] ).    (10)
!
!
! ---------------------------------------------------------------
! This routine will call those that evaluate the integrals in the
! numerator and denominator of (10).
!
! The denominator need be evaluated only once for each type
! of particle since it depends only upon the prescribed size
! distribution function.  In the numerator however, since the
! function G does depend upon T and p as well as r, this
! integral must be evaluated at all nodes.
! ---------------------------------------------------------------
!
! *------------------
! DATA FILES REQUIRED
! *------------------
!
!     seasalt.75RH                   sea salt
!
!     ******************************************************************
!     * Arrays used only in the routine.
!       1) pgcm       pressure at nodes                  [bars]
!       2) cbar       mean molecular speed               [cm/s]
!       3) dmol       molecular diffusion coefficient   [cm2/s]
!       4) lamb       mean free path                       [cm]
!       5) alph       sticking coefficient                  [ ]
!     ******************************************************************

      subroutine Do_So2toPar_Rates  &
     &  (itloop, pgcm, tgcm,  &
     &   relh, aqua_infile_name,  &
     &   gmair, cp, rtcd, rtcs)

      implicit none

#     include "setkin_par.h"
#     include "sulfchem.h"

!     ----------------------
!     GMIMOD-related values.
!     ----------------------

      character*128 :: aqua_infile_name

      integer :: itloop

      real*8  :: pgcm (itloop)
      real*8  :: tgcm (itloop)
      real*8  :: relh (itloop)
      real*8  :: gmair(itloop)
      real*8  :: cp   (itloop, NSP)

!     ----------------------
!     Argument declarations.
!     ----------------------

      integer :: npi, nr(4)

      real*8  :: denom(8)
      real*8  :: rtcd(itloop, 4)
      real*8  :: rtcs(itloop, 4)


! -------------------------------------------------------------
!  *** LOCAL  VARIABLES WERE DECLARED HERE.
! -------------------------------------------------------------

!     the molecular diffusion coefficient (cm2/s) for SO2.
      real*8 :: dmol(itloop)

!     the molecular mean free path (cm) for SO2.
      real*8 :: lamb(itloop)

!     the fraction of vapor molecules colliding
!     with an aerosol particle that "stick" to it.
      real*8 :: alph(itloop)

!     the integral in the numerator of (10), units [cm3/s]
      real*8 :: nmrtr(itloop, 12)

      real*8 :: r    (2001, 4)
      real*8 :: rtwk (2001)
      real*8 :: dndrt(2001, 4)
      real*8 :: rdry (2001, 4)

      dmol  = 0.0d0
      lamb  = 0.0d0
      alph  = 0.0d0
      nmrtr = 0.0d0

      r     = 0.0d0
      rtwk  = 0.0d0
      dndrt = 0.0d0
      rdry  = 0.0d0


!.... calculate dmol and lamb
!           ========
       call chemst70  &
!           ========
     &   (itloop, tgcm, pgcm,  &
     &    dmol, lamb)


!     ******************************************************************
!     * The integrals appearing in the numerator and denominator of the
!       expression used to evaluate the SO2 -> particle mass transfer
!       rate coefficients do depend upon the particle size distribution
!       function.
!     * 990511 - There are now two sets of numerators for dust,
!                corresponding to the two values of the sticking coef.,
!                "alph", for relative humidity > 0.5 and > 0.05.
!     ******************************************************************

! ----
! dust
! ----

!.... the size distribution function and the denominator integral
!          ========
      call chemst61  &
!          ========
     &  (denom, nr, r, dndrt)

!.... dust sticking coefficient for relative humidity < 0.5

      npi = 1

!          ========
      call chemst64  &
!          ========
     &  (npi, tgcm, itloop, alph)


!.... numerator integral
!          ========
      call chemst63  &
!          ========
     &  (nr, npi, itloop, alph,  &
     &   dmol, lamb, r, rtwk, dndrt, nmrtr)


!.... dust sticking coefficient for relative humidity > 0.5

      npi = 5

!          ========
      call chemst64  &
!          ========
     &  (npi, tgcm, itloop, alph)


!.... numerator integral
!          ========
      call chemst63  &
!          ========
     &  (nr, npi, itloop, alph,  &
     &   dmol, lamb, r, rtwk, dndrt, nmrtr)


! --------
! sea salt
! --------

!.... the size distribution function and the denominator integral
!          ========
      call chemst62  &
!          ========
     &  (denom, nr, r, dndrt, rdry, aqua_infile_name)

!.... sea salt sticking coefficient

      npi = 9

!          ========
      call chemst64  &
!          ========
     &  (npi, tgcm, itloop, alph)


!.... numerator integral
!          ========
      call chemst63  &
!          ========
     &  (nr, npi, itloop, alph,  &
     &   dmol, lamb, r, rtwk, dndrt, nmrtr)


!     ----------------------------------------------------------
!     calculate the SO2 mass transfer rate coefficients for dust
!     and sea salt
!     ----------------------------------------------------------
!
!          ======
      call chem06  &
!          ======
     &  (itloop, denom, nmrtr, relh,  &
     &   gmair, cp, rtcd, rtcs)


      return

      end


!-----------------------------------------------------------------------------
!
! ROUTINE
!   chem06
!
! DESCRIPTION
!
!     *****************************************************************
!     * Calculate the SO2 mass transfer rate coefficients.
!     *****************************************************************
!
!     *   denom :   particle density times the integral in the
!                   denominator in (10)    [gm]
!     *   nmrtr :   the integral in the numerator of (10)
!                   (see chemst60)
!     *   relh  :   relative humidity [0, 1]
!     *   rtcd  :   (s-1) rate coef. for SO2 -> dust transfer
!     *   rtcs  :   (s-1) rate coef. for SO2 -> sslt transfer
!
!     *****************************************************************
!     * denom has units [gm]
!       nmrtr has units [cm3/s]
!       cp    has units [gm spec./gm air for aerosols]
!       gmair has units [gm air/cm3]
!       "rtc" has units [s-1]
!     *****************************************************************
!
!     *----------------------------------------------------------------

      subroutine chem06  &
     &  (itloop, denom, nmrtr, relh,  &
     &   gmair, cp, rtcd, rtcs)

      implicit none

#     include "setkin_par.h"
#     include "sulfchem.h"

!     ----------------------
!     GMIMOD-related values.
!     ----------------------

      integer :: itloop

      real*8  :: relh (itloop)
      real*8  :: gmair(itloop)
      real*8  :: cp   (itloop, NSP)

!     ----------------------
!     Argument declarations.
!     ----------------------

      real*8  :: denom(8)
      real*8  :: nmrtr(itloop, 12)
      real*8  :: rtcd (itloop, 4)
      real*8  :: rtcs (itloop, 4)


!.... dust, two possibilities


        where ( relh(:) .lt. 0.5d0 )

          rtcd(:, 1) = nmrtr(:, 1) * gmair(:)  &
     &                 * cp(:, IDUST1) / denom(1)
          rtcd(:, 2) = nmrtr(:, 2) * gmair(:)  &
     &                 * cp(:, IDUST2) / denom(2)
          rtcd(:, 3) = nmrtr(:, 3) * gmair(:)  &
     &                 * cp(:, IDUST3) / denom(3)
          rtcd(:, 4) = nmrtr(:, 4) * gmair(:)  &
     &                 * cp(:, IDUST4) / denom(4)

        elsewhere

          rtcd(:, 1) = nmrtr(:, 5) * gmair(:)  &
     &                 * cp(:, IDUST1) / denom(1)
          rtcd(:, 2) = nmrtr(:, 6) * gmair(:)  &
     &                 * cp(:, IDUST2) / denom(2)
          rtcd(:, 3) = nmrtr(:, 7) * gmair(:)  &
     &                 * cp(:, IDUST3) / denom(3)
          rtcd(:, 4) = nmrtr(:, 8) * gmair(:)  &
     &                 * cp(:, IDUST4) / denom(4)

        end where


!.... sea salt


        rtcs(:, 1) = nmrtr(:, 9) * gmair(:)  &
     &               * cp(:, ISSLT1) / denom(5)
        rtcs(:, 2) = nmrtr(:,10) * gmair(:)  &
     &               * cp(:, ISSLT2) / denom(6)
        rtcs(:, 3) = nmrtr(:,11) * gmair(:)  &
     &               * cp(:, ISSLT3) / denom(7)
        rtcs(:, 4) = nmrtr(:,12) * gmair(:)  &
     &               * cp(:, ISSLT4) / denom(8)



      return

      end


!-----------------------------------------------------------------------------
!
! ROUTINE
!   chemst61
!
! DESCRIPTION
!
!     *****************************************************************
!     * This is a tri-modal distribution for desert dust.
!     *
!     * For desert dust the parameters for the 3 distributions are:
!     * [de Reus et al., 2000]
!     *    nd (cm-3) =   33.,  200.,    20
!     *    rg (cm)   = 1.e-6, 4.5 e-6, 2.75 e-5
!     *    log(s)    = 0.146, 0.204,   0.398 (in order)
!     ******************************************************************
!
!     ******************************************************************
!      COMMENT ON THE NORMALIZATION OVER DIFFERENT SIZE INTERVALS,
!      FOR DUST AND SEA SALT.
!
!      Particles of some sort, say dust, are described over an interval
!      r1 --> r2 by the function dN/dlogr, which presumably represents
!      some sort of observational data.  The total number of particles
!      is
!
!             Nt = integral[ (dN/dlogr) dlogr ].
!
!      If I divide the interval at some arbitrary point, ri, I can
!      calculate the number of particles in each sub-interval
!
!             Nt1 = integral[ (dN/dlogr) dlogr ]
!                   r1-->ri
!      and
!             Nt2 = integral[ (dN/dlogr) dlogr ] .
!                   ri-->r2
!
!      I could write similar expressions for total mass, Mt, and mass
!      in the sub-intervals Mt1 and Mt2.
!
!     ******************************************************************
!
!     ******************************************************************
!     * 1) The distribution functions will be evaluated for r over
!     *    the five size ranges.
!     * 2) Each distribution function will then be normalized.
!     *****************************************************************

       subroutine chemst61  &
     &   (denom, nr, r, dndrt)

      implicit none

#     include "sulfchem.h"

!     ----------------------
!     Argument declarations.
!     ----------------------

      integer :: nr(4)

      real*8  :: denom(8)
      real*8  :: r(2001, 4)
      real*8  :: dndrt(2001, 4)

! -------------------------------------------------------------
!  *** LOCAL  VARIABLES WERE DECLARED HERE.
! -------------------------------------------------------------

      real*8  ::  &
     &  rhodust,  &
     &  nd1, nd2, nd3,  &
     &  rg1, rg2, rg3,  &
     &  logs1, logs2, logs3,  &
     &  ri, r1, r2

      real*8  :: ri1(4), ri2(4)
      real*8  :: totvol(4)

      data rhodust / 2.5d0 /
      data ri1 / 5.d-6, 6.3d-5, 1.3d-4, 2.5d-4 /
      data ri2 / 6.3d-5,1.3d-4, 2.5d-4, 1.d-3  /

!.... dust density 2.5 gm cm-3

      rhop = rhodust

!.... normalization factors

      nd1 = 33.d0
      nd2 = 200.d0
      nd3 = 20.d0

!.... radius (cm)

      rg1 = 1.0d-6
      rg2 = 4.5d-6
      rg3 = 2.75d-5

!.... log (base 10) of the standard deviation

      logs1 = 0.146d0
      logs2 = 0.204d0
      logs3 = 0.398d0


!     *****************************************************************
!     * Compute dn/dlogr over the range 0.05 um < r < 10 um.
!     * dN/dlogr, dndrt, is computed for nr radii.
!          nr(1), r1 = 5.e-6  - 6.3e-5,  particle bin #1
!          nr(2), r2 = 6.3e-5 - 1.3e-4,  particle bin #2
!          nr(3), r3 = 1.3e-4 - 2.5e-4,  particle bin #3
!          nr(4), r4 = 2.5e-4 - 1.e-3,   particle bin #4
!     *****************************************************************

      do i = 1, 4

        r1 = ri1(i)
        r2 = ri2(i)

!.... compute r(,i) for particle bin #i
!            ========
        call chemst74  &
!            ========
     &    (r1, r2, r(1,i), nr(i))


        do m = 1, nr(i)

          ri = r(m,i)

          dndrt(m,i) = ( nd1 / ( logs1 * sqrt( 2.d0 * pi ) ) )  &
     &               * exp( - 0.5d0 * ( log10( ri/rg1 ) / logs1 ) ** 2 )  &
     &               + ( nd2 / ( logs2 * sqrt( 2. * pi ) ) )  &
     &               * exp( - 0.5d0 * ( log10( ri/rg2 ) / logs2 ) ** 2 )  &
     &               + ( nd3 / ( logs3 * sqrt( 2. * pi ) ) )  &
     &               * exp( - 0.5d0 * ( log10( ri/rg3 ) / logs3 ) ** 2 )

        enddo

      enddo


!     *****************************************************************
!     * Calculate the total number of particles in each size range, and
!       use these to normalize the dN/dlogr.
!     * dndrt (=dn/dlogr) will then be dimensionless & normalized.
!     * Compute total volume to get denom.
!     *****************************************************************

      do i = 1, 4

!            ========
        call chemst71  &
!            ========
     &    (nr(i), r(1,i), dndrt(1,i))

!            ========
        call chemst72  &
!            ========
     &    (nr(i), r(1,i), dndrt(1,i), totvol(i))

        denom(i) = rhop * totvol(i)

      enddo


      return

      end


!-----------------------------------------------------------------------------
!
! ROUTINE
!   chemst62
!
! DESCRIPTION
!
!     *****************************************************************
!     * Sea salt
!
!     * The distribution will be read from the file "seasalt.75RH",
!     * which was produced from a RH = 0% file provided by Wenwei using
!     * the wet radius code "wetrad.970624".
!
!     * Column 1 in the file is particle DRY radius (um), column 2 the
!     * particle WET radius and column 3 is dN/dlogr (cm-3).
!
!     * The radii and distribution function will be in the ,1) locations
!     * of the r(,) and dndrt(,) arrays.
!
!     * See "dfnso4" for the conversion from dN/dlogrd to dN/dlogrw.
!     ******************************************************************
!
!-----------------------------------------------------------------------------

      subroutine chemst62  &
     &  (denom, nr, r, dndrt, rdry, aqua_infile_name)

      use GmiASCIIoperations_mod, only : AsciiOpenRead

      implicit none

#     include "sulfchem.h"

!     ----------------------
!     Argument declarations.
!     ----------------------

      character*(*) :: aqua_infile_name

      integer :: nr(4)

      real*8  :: denom(8)
      real*8  :: r    (2001, 4)
      real*8  :: dndrt(2001, 4)
      real*8  :: rdry (2001, 4)

! -------------------------------------------------------------
!  *** LOCAL  VARIABLES WERE DECLARED HERE.
! -------------------------------------------------------------

      integer :: lun
      logical, save :: first = .true.

      real*8  :: fwk(3)

      real*8  :: rhosslt
      real*8  :: r1(4)
      real*8  :: totvol(4)

      data rhosslt / 2.165 /
      data r1 / 0.63d0, 1.26d0, 2.5d0, 10.d0/   ! um


!.... sea salt density (gm cm-3)

      rhop = rhosslt


!     ==========
      if (first) then
!     ==========

        first = .false.


!       Open and read the sea salt distribution file.
!            ===========
        call AsciiOpenRead  &
!            ===========
     &    (lun, aqua_infile_name)

!....   find the first entry

        do 100 i = 1, 9
  100   read(lun,*)

        read(lun,10) (fwk(j), j=1, 3)

!....   convert radius from (um) to (cm) - store the 3 variables

!       ===========
        do i = 1, 4
!       ===========

          nr1 = 1

          rdryss(nr1,i)  = 1.d-4 * fwk(1)
          rss(nr1,i)     = 1.d-4 * fwk(2)
          dndrtss(nr1,i) = fwk(3)

!....   read the remainder of the entries

  102     read(lun,10,end=990) (fwk(j), j=1,3)

          nr1 = nr1 + 1

!....   convert radius from (um) to (cm) - store the 3 variables

          rdryss(nr1,i)  = 1.d-4 * fwk(1)
          rss(nr1,i)     = 1.d-4 * fwk(2)
          dndrtss(nr1,i) = fwk(3)

!....   this list ends when rdry = 1 um, 1.e-6 is a little fuzz factor

          if( fwk(1) .lt. r1(i) - 1.d-6 ) go to 102

 990      nrss(i) = nr1

!       ======
        end do
!       ======

 10     format(3(1pe12.4))

        Close (lun)

      end if


!     ===========
      do i = 1, 4
!     ===========

      nr(i)      = nrss(i)
      rdry (:,i) = rdryss (:,i)
      r    (:,i) = rss    (:,i)
      dndrt(:,i) = dndrtss(:,i)

!     *---------------------------------------------------------------*
!     * Compute the unscaled number, dndrint, and use it to normalize
!       the distribution.
!     * dndrt (=dn/dlogr) will then be dimensionless & normalized.
!     *---------------------------------------------------------------*

!          ========
      call chemst71  &
!          ========
     &  (nr(i), rdry(1,i), dndrt(1,i))

!          ========
      call chemst72  &
!          ========
     &  (nr(i), rdry(1,i), dndrt(1,i), totvol(i))

      denom(4+i) = rhop * totvol(i)


!     *---------------------------------------------------------------*
!     * Use wet and dry radii, r(,1) and rdry, and the dry radius
!       distribution function, dndrt, to generate the wet radius
!       distribution function.
!     * Normalization will be checked again.
!     *---------------------------------------------------------------*

!.... change of variable
!          ========
      call chemst73  &
!          ========
     &  (nr(i), rdry(1,i), r(1,i), dndrt(1,i))

!.... check normalization
!          ========
      call chemst71  &
!          ========
     &  (nr(i), r(1,i), dndrt(1,i))

!     ======
      end do
!     ======


      return

      end


!-----------------------------------------------------------------------------
!
! ROUTINE
!   chemst63
!
! DESCRIPTION
!
!     **************************************************************
!     *... dust (tri-modal, dry particles) and
!          sea salt (tri-modal, soluble particles)
!     *    Evaluate the integral in the numerator of (10)
!
!            integral[ G(r,T,p) (dn/dlogr) dlogr ]
!
!     **************************************************************
!
!     alph        The fraction of vapor molecules colliding with an
!                 aerosol particle that "stick" to it.
!
!     dmol        The molecular diffusion coefficient (cm2/s) for SO2
!
!     dndrt       The normalized number distribution  (dimensionless)
!
!     rhop        The mass density of a particle      (gm/cm3)
!
!     r           aerosol particle radius             (cm)
!
!     nmrtr       the numerator of (10), excluding the
!                 mass mixing ratio.
!
!     lamb        (cm) the molecular mean free path for SO2.  This
!                 relation is somewhat arbitray, but is required for
!                 the size dependent mass transfer coefficient to
!                 have the correct limit for large Knudsen number.
!                 = 3 * dmol / cbar
!
!     npi         the location of the small particle values indexed
!                 to 1 for dust1.
!
!     xnud        the (dimensionless) Knudsen number is defined as
!                 = lamb / r
!
!     **************************************************************
!     * AT EACH NODE:
!     * 1) Compute the size dependent mass transfer rate coef.
!     * 2) Calculate the integral of the product of this coef.
!          with the distribution function, dndrt.
!          This is the numerator of (10), excluding the
!          mass mixng ratio
!     **************************************************************
!
!-----------------------------------------------------------------------------

      subroutine chemst63  &
     &  (nr, npi, itloop, alph,  &
     &   dmol, lamb, r, rtwk, dndrt, nmrtr)

      implicit none

#     include "sulfchem.h"

!     ----------------------
!     GMIMOD-related values.
!     ----------------------

      integer :: itloop

!     ----------------------
!     Argument declarations.
!     ----------------------

      integer :: npi
      integer :: nr(4)

      real*8  :: dmol (itloop)
      real*8  :: lamb (itloop)
      real*8  :: alph (itloop)
      real*8  :: nmrtr(itloop, 12)

      real*8  :: r(2001, 4)
      real*8  :: rtwk(2001)
      real*8  :: dndrt(2001, 4)

! -------------------------------------------------------------
!  *** LOCAL  VARIABLES WERE DECLARED HERE.
! -------------------------------------------------------------

      real*8  alpha, xnud, ctr

      ftd = 4.d0 / 3.d0
      fpi = 4.d0 * pi

!
! 990521 - To better accomodate small alphas, "alph", I will
!     multiply the numerator and denominator of the size
!     dependent mass transfer coefficient by alpha.  This
!     will change the expressions for "dmoltrm", "alphtrm",
!     and "rtwk".

      do ijk = 1, itloop

        alpha   = alph(ijk)
        dmoltrm = fpi * dmol(ijk) * alpha
        alphtrm = ftd * ( 1.d0 - alpha )

        do i = 1, 4

!....   G, the size dependent mass trans. coef.

          do m = 1, nr(i)

            xnud = lamb(ijk) / r(m,i)
            ctr  = ( 0.71d0 + ftd * xnud ) / ( 1.d0 + xnud )

            rtwk(m) = r(m,i) * dmoltrm  &
     &              / ( alpha + xnud * ( ctr * alpha + alphtrm ) )

          enddo

!....   G (dn/dlogr)

          do m = 1, nr(i)
            rtwk(m) = rtwk(m) * dndrt(m,i)
          enddo

!....   integral[ G (dn/dlogr) dlogr ]
!              ========
          call chemst75  &
!              ========
     &      (nr(i), r(1,i), rtwk(1), nmrtr(ijk, npi+i-1))

        enddo

      enddo


      return

      end


!-----------------------------------------------------------------------------
!
! ROUTINE
!   chemst64
!
! DESCRIPTION
!
!     **************************************************************
!     * Specify the sticking coefficient.
!     **************************************************************
!
!     ***** input
!
!     *    npi  :   particle type indicator
!                   = 1  - dust, RH < 0.5
!                   = 5  - dust, RH > 0.5
!                   = 9  - sea salt
!
!     *    tgcm :   air temperature (K)
!
!     *    alph :   the fraction of vapor molecules colliding
!                   with an aerosol particle that "stick" to it.
!
!     *-------------------------------------------------------------
!
!     * The sticking coefficient for dust now depends upon the
!       relative humidity, RH, see Dentener et al, "Role of mineral
!       aerosol as a reactive surface in the global troposphere",
!       JGR V 101, No. D17, October 20, 1996.
!
!     * For RH < 0.5, alph = 3.e-4
!           RH > 0.5, alph = 0.1
!
!     * Two values for the dust mass transfer rate coefficient will
!       be calculated using these alph's.
!
!
!     *-------------------------------------------------------------

      subroutine chemst64  &
     &  (npi, tgcm, itloop, alph)

      implicit none

#     include "sulfchem.h"

!     ----------------------
!     GMIMOD-related values.
!     ----------------------

      integer :: itloop

      real*8  :: tgcm(itloop)
      real*8  :: alph(itloop)

!     ----------------------
!     Argument declarations.
!     ----------------------

      integer :: npi

! -------------------------------------------------------------
!  *** LOCAL  VARIABLES WERE DECLARED HERE.
! -------------------------------------------------------------

      real*8, parameter ::  &
     &  stckcoef = 1.d-4


      if( npi .eq. 1 ) then

!.... dust, RH < 0.5

!     990521 - no relative humidity dependence

        do ijk = 1, itloop
          alph(ijk) = stckcoef
        end do

      elseif( npi .eq. 5 ) then

!.... dust, RH > 0.5

!     990521 - no relative humidity dependence

        do ijk = 1, itloop
          alph(ijk) = stckcoef
        end do

      else

!.... sea salt

        do ijk = 1, itloop
          if( tgcm(ijk) .le. 273.d0 ) then
            alph(ijk) = 0.1d0
          elseif( tgcm(ijk) .ge. 288.d0 ) then
            alph(ijk) = 0.03d0
          else
            alph(ijk) = 0.03d0 + 0.07d0 * ( 288.d0 - tgcm(ijk) ) / 15.d0
          endif
        end do

      endif


      return

      end


!-----------------------------------------------------------------------------
!
! ROUTINE
!   chemst70
!
! DESCRIPTION
!
!     **************************************************************
!     * Note:  This routine is independent of the type of
!     *        particles being treated.
!     *
!     * 1) Compute SO2 molecular speed based on temperatures.
!     * 2) Compute the molecular duffusion coefficient for SO2,
!     *    as a function of temperatures and pressures
!     * 3) Compute the molecular mean free path number at each
!     *    model node.
!     **************************************************************
!
!     tgcm    :   (K) air temperature.
!     pgcm    :   (bars) grid point pressures
!     gc      :   gas constant, 8.317e7 (ergs/K/mole or cm3 dynes/cm3),
!     wtmso2  :   = 64.06d0, molecular weight of SO2
!     cbar    :   the mean speed of SO2 molecules (cm/sec)
!     dmol    :   the molecular diffusion coefficient (cm2/s) for SO2.
!     lamb    :   (cm) the molecular mean free path for SO2.  This
!                 relation is somewhat arbitray, but is required for
!                 the size dependent mass transfer coefficient to
!                 have the correct limit for large Knudsen number.
!                 = 3 * dmol / cbar
!
!     The (dimensionless) Knudsen number is defined as
!
!         mean free path / particle radius = lamb / r
!
!-----------------------------------------------------------------------------

       subroutine chemst70  &
     &   (itloop, tgcm, pgcm,  &
     &    dmol, lamb)

      implicit none

#     include "sulfchem.h"

!     ----------------------
!     GMIMOD-related values.
!     ----------------------

      integer :: itloop

      real*8  :: tgcm(itloop)
      real*8  :: pgcm(itloop)

!     ----------------------
!     Argument declarations.
!     ----------------------

      real*8  :: dmol(itloop)
      real*8  :: lamb(itloop)

! -------------------------------------------------------------
!  *** LOCAL  VARIABLES WERE DECLARED HERE.
! -------------------------------------------------------------

!     the mean speed of SO2 molecules (cm/sec)
      real*8 :: cbar(itloop)

!.... the parameters for calculating the molecular diffusion coefficient

      data dmol0, temp0, pres0 / 0.1d0, 298.d0, 1.d0/

      cbar = 0.0d0

!.... compute the mean molecular speed (cm/s) of SO2
!     wtmso2 is the moleculare weight of SO2

      do ijk = 1, itloop
        cbar(ijk) = sqrt( 8.d0 * gc * tgcm(ijk)  &
     &                    / (pi * wtmso2) )
      end do

!.... compute the molecular diffusion coefficient (cm2/s) for SO2
!     This is from the Chapman-Enskog theory, eq. 8.8 in Seinfeld.
!     I have used the dependence upon temperature and pressure to
!     relate this to its value at some refererence level.

      do ijk = 1, itloop
        dmol(ijk) = dmol0 * ( pres0 / pgcm(ijk) )  &
     &                    * ( tgcm(ijk) / temp0 ) ** 1.5
      end do

!.... compute the mean free path, lamb (cm)

      do ijk = 1, itloop
        lamb(ijk) = 3. * dmol(ijk) / cbar(ijk)
      end do


      return

      end


!-----------------------------------------------------------------------------
!
! ROUTINE
!   chemst71
!
! DESCRIPTION
!
!      subroutine dfn01( nr, r, dndlogr )
!
!     *****************************************************************
!     * Compute the unscaled number, dndrint, and use it to
!       normalize the distribution.
!       dndlogr (=dn/dlogr) will then be dimensionless & normalized.
!     *****************************************************************
!-----------------------------------------------------------------------------

      subroutine chemst71  &
     &  (nr, r, dndlogr)

      implicit none

#     include "sulfchem.h"

!     ----------------------
!     Argument declarations.
!     ----------------------

      integer :: nr       ! the number of r's computed

      real*8  :: r(2001), dndlogr(2001)

! -------------------------------------------------------------
!  *** LOCAL  VARIABLES WERE DECLARED HERE.
! -------------------------------------------------------------

      real*8  :: dndlrint


!     *****************************************************************
!     The normalization of dndlogr.
!     dn/dlogr can then be thought of as dimensionless.
!     *****************************************************************

!.... unscaled number
!          ========
      call chemst75  &
!          ========
     &  (nr, r(1), dndlogr(1), dndlrint)

!.... normalize - eliminate dimensions

      do m = 1, nr
        dndlogr(m) = dndlogr(m) / dndlrint
      enddo


      return

      end


!-----------------------------------------------------------------------------
!
! ROUTINE
!   chemst72
!
! DESCRIPTION
!
!      subroutine dfn02( nr, r, dndlogr, totvol )
!
!     *****************************************************************
!     * With the normalized dndlogr, compute the total volume,
!       totvol.
!       Multiplying totvol by a particle density gives denom.
!     *****************************************************************
!-----------------------------------------------------------------------------

      subroutine chemst72  &
     &  (nr, r, dndlogr, totvol)

      implicit none

#     include "sulfchem.h"

!     ----------------------
!     Argument declarations.
!     ----------------------

      integer :: nr       ! the number of r's computed

      real*8  :: totvol
      real*8  :: r(2001), dndlogr(2001)

! -------------------------------------------------------------
!  *** LOCAL  VARIABLES WERE DECLARED HERE.
! -------------------------------------------------------------

      real*8 :: rtwk(2001)

      rtwk = 0.0d0


      ftdpi = (4.d0 / 3.d0) * pi


!     *****************************************************************
!     * Compute the total normalized particle volume (cm3), totvol.
!       This is the denominator of (10)
!     *****************************************************************

      do m = 1, nr
        rtwk(m) = (ftdpi * r(m) ** 3) * dndlogr(m)
      enddo

!          ========
      call chemst75  &
!          ========
     &  (nr, r(1), rtwk(1), totvol)


      return

      end


!-----------------------------------------------------------------------------
!
! ROUTINE
!   chemst73
!
! DESCRIPTION
!
!      subroutine dfn03( nr, rd, rw, dndlogr )
!
!     *****************************************************************
!     * Use wet and dry radii, rw & rd and the dry radius distribution
!       function, dndlogr, to generate the wet radius distribution
!       function.
!     * The wet distribution function is returned as dndlogr.
!     *****************************************************************
!
!-----------------------------------------------------------------------------

      subroutine chemst73  &
     &  (nr, rd, rw, dndlogr)

      implicit none

#     include "sulfchem.h"

!     ----------------------
!     Argument declarations.
!     ----------------------

      integer :: nr       ! the number of r's computed

      real*8  :: rd(2001), rw(2001), dndlogr(2001)


! -------------------------------------------------------------
!  *** LOCAL  VARIABLES WERE DECLARED HERE.
! -------------------------------------------------------------

      real*8 :: dlrddlrw(2001)

      dlrddlrw = 0.0d0

!.... compute dlog(rd)/dlog(rw)

      dlrddlrw( 1) = log10( rd(2) / rd(1) )  &
     &             / log10( rw(2) / rw(1) )
      dlrddlrw(nr) = log10( rd(nr) / rd(nr-1) )  &
     &             / log10( rw(nr) / rw(nr-1) )

      do i = 2, nr-1

        dlrddlrw(i) = 0.5d0 * log10( rd(i+1) / rd(i) )  &
     &                    / log10( rw(i+1) / rw(i) )  &
     &              + 0.5d0 * log10( rd(i) / rd(i-1) )  &
     &                    / log10( rw(i) / rw(i-1) )

      enddo

!.... perform the change of variable

      do i = 1, nr
        dndlogr(i) = dndlogr(i) * dlrddlrw(i)
      enddo


      return

      end


!-----------------------------------------------------------------------------
!
! ROUTINE
!   chemst74
!
! DESCRIPTION
!
!     subroutine rgen( r1, rnr, r, nr )
!
!     *****************************************************************
!     * Generates a set of nr radii, r.
!     *****************************************************************
!
!     nd :   the number of values in a decade of r
!     nr :   the number of r's computed
!     r1 :   minimum value of r                                  (cm)
!     rnr:   maximum value of r                                  (cm)
!     r  :   the nr values computed from r1 through rnr          (cm)
!
!     The relation defining r is
!
!         r(i) = r1 * 10 ** [ beta * ( i - 1 ) ]
!
!-----------------------------------------------------------------------------

      subroutine chemst74  &
     &  (r1, rnr, r, nr)

      implicit none

#     include "sulfchem.h"

!     ----------------------
!     Argument declarations.
!     ----------------------

      integer :: nr

      real*8  :: r(2001), r1, rnr

! -------------------------------------------------------------
!  *** LOCAL  VARIABLES WERE DECLARED HERE.
! -------------------------------------------------------------

      integer :: nd

      real*8  :: beta

!     Defining the parameters.

!.... set nd

      nd = 1001

!.... determine beta

      beta = 1.d0 / float( nd - 1 )

!.... determine nr

      nr   = 1.d0 + ( nd - 1.d0 ) * log10( rnr / r1 )

!     Generate the nr values of r.

      do i = 1, nr

        r(i) = r1 * 10 ** ( beta * float( i - 1 ) )

      enddo


      return

      end


!-----------------------------------------------------------------------------
!
! ROUTINE
!   chemst75
!
! DESCRIPTION
!
!      subroutine dlogrint( nr, r, dndlogr, dlri )
!
!     *****************************************************************
!     * nr : no. of variables passed
!     * r  : radius (cm)
!     * dndlogr
!     *    : dN/dlogr (cm-3), This can actually be any function
!     *                       of r to be integrated wrt dlogr.
!     * Look at the integral, dlri, of dndlogr * dlogr
!     * over the range of r.
!     *****************************************************************
!
!-----------------------------------------------------------------------------

      subroutine chemst75  &
     &  (nr, r, dndlogr, dlri)

      implicit none

#     include "sulfchem.h"

!     ----------------------
!     Argument declarations.
!     ----------------------

      integer :: nr       ! the number of r's computed

      real*8  :: dlogr, dndlrb
      real*8  :: r(2001), dndlogr(2001), dlri


      dlri = 0.d0

      do i = 1, nr-1

        dlogr = log10( r(i+1) ) - log10( r(i) )

        dndlrb = 0.5d0 * ( dndlogr(i+1) + dndlogr(i) )

        if( dndlrb .gt. 1.d-20 ) then
          dlri = dlri + dndlrb * dlogr
        endif

      enddo


      return

      end


