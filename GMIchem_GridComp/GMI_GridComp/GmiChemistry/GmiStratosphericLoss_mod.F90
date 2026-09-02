!------------------------------------------------------------------------------
!BOP
!
! !MODULE: GmiStratosphericLoss_mod
!
      module GmiStratosphericLoss_mod
!
! !USES:
      use GmiArrayBundlePointer_mod, only : t_GmiArrayBundle
      use GmiTimeControl_mod, only : GmiSplitDateTime
!
      implicit none
!
! !PUBLIC FUNCTIONS:
      public  :: updateStratosphericLoss
!
#     include "gmi_phys_constants.h"
#     include "gmi_time_constants.h"
!
!EOP
!------------------------------------------------------------------------------
      contains
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: updateStratosphericLoss
!
! !INTERFACE:
! 
      subroutine updateStratosphericLoss (chemintv, dlatr, loss_freq,          &
     &                 lossData, concentration, averagePressEdge,              &
     &                 loss_freq_opt, nymd, kmin_loss, kmax_loss,              &
     &                 const_init_val, num_species, pr_diag, loc_proc, i1, i2, &
     &                 ju1, j2, k1, k2, ivert, ju1_gl, j2_gl, KDIM, JDIM, MDIM)
!
      implicit none
!
! !INPUT PARAMETERS:
      logical, intent(in   ) :: pr_diag
      integer, intent(in   ) :: loc_proc
      integer, intent(in   ) :: num_species
      integer, intent(in   ) :: i1, i2, ju1, j2, k1, k2, ivert, ju1_gl, j2_gl
      integer, intent(in   ) :: loss_freq_opt, KDIM, JDIM, MDIM
      integer, intent(in   ) :: nymd
      integer, intent(in   ) :: kmin_loss, kmax_loss
      real*8 , intent(in   ) :: chemintv          ! chemistry time step (s)
      real*8 , intent(in   ) :: averagePressEdge (k1-1:k2)
      real*8 , intent(in   ) :: dlatr    (ju1_gl:j2_gl)
!                             ! latitude of zone center in latitude direction (rad)
      real*8 , intent(in   ) :: const_init_val(num_species)
      real*8 , intent(in   ) :: lossData(KDIM, JDIM, MDIM, 1)
!
! !OUTPUT PARAMETERS:
      real*8 , intent(  out) :: loss_freq(i1:i2, ju1:j2, k1:k2, num_species) ! loss_rates (s^-1)
!
! !INPUT/OUTPUT PARAMETERS:
      type (t_GmiArrayBundle), intent(inOut) :: concentration(num_species)
!
! !DESCRIPTION:
!   Updates the simple loss chemistry.
!
! !LOCAL VARIABLES:
      integer :: idumday, idumyear, month
      integer :: ic, ik, numj
      integer, save :: last_month = 0
!
! !REVISION HISTORY:
!  Original code.
!
!EOP
!-------------------------------------------------------------------------
!BOC
      !======================
      if (loss_freq_opt == 2) then
      !======================
         call GmiSplitDateTime (nymd, idumyear, month, idumday)

         if (month /= last_month) then
            last_month = month
            if (pr_diag) Write (6,*) 'Stratl_Gmi called by ', loc_proc

            numj = j2_gl - ju1_gl + 1

            call Stratl_Gmi (month, dlatr, loss_freq, lossData, numj,          &
     &                  num_species, averagePressEdge, i1, i2, ju1, j2, k1, k2,&
     &                  ivert, ju1_gl, j2_gl, KDIM, JDIM, MDIM)
         end if
      end if

      do ik = k1, k2
         if ((ik >= kmin_loss) .and. (ik <= kmax_loss)) then
            do ic = 1, num_species
               concentration(ic)%pArray3D(:,:,ik) =  &
     &                           concentration(ic)%pArray3D(:,:,ik) -  &
     &                          (concentration(ic)%pArray3D(:,:,ik) *  &
     &                           loss_freq(:,:,ik,ic) * chemintv)
            end do

         else if (ik < kmin_loss) then
            do ic = 1, num_species
               concentration(ic)%pArray3D(:,:,ik) = const_init_val(ic)
            end do
         end if
      end do

      return

      end subroutine updateStratosphericLoss
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Stratl_Gmi
!
! !INTERFACE:
!
      subroutine Stratl_Gmi (month, dlatr, loss_freq, lossData,                &
     &                numj, num_species, plevl, i1, i2, ju1, j2, k1, k2, ivert,&
     &                ju1_gl, j2_gl, KDIM, JDIM, MDIM)

      implicit none

!
! !INPUT DESCRIPTION:
      integer, intent(in) :: i1, i2, ju1, j2, k1, k2, ivert, ju1_gl, j2_gl
           ! current month of simulation
      integer, intent(in) :: month
           ! number of species
      integer, intent(in) :: num_species
           ! number of global latitudes
      integer, intent(in) :: numj
      integer, intent(in) :: KDIM, JDIM, MDIM
           ! atitude of zone center in latitude direction  (deg)
      real*8 , intent(in) :: dlatr(ju1_gl:j2_gl)
           ! pressure at edges of vertical grid (bottom up) (mb)
      real*8 , intent(in) :: plevl(ivert+1)
           ! 4D array of loss frequencies known on preset grid
           ! defined for 18 lats (85S, 75S, ...85N) & 12 months
      real*8 , intent(in) :: lossData(KDIM, JDIM, MDIM, 1)
!
! !OUTPUT PARAMETERS:
           ! loss rates (s^-1)
      real*8 , intent(out) :: loss_freq(i1:i2,ju1:j2,k1:k2,num_species)
!
! !DESCRIPTION:
!   Monthly fixup of chemistry param's.
!
! !LOCAL VARIABLES:
      integer :: i, j, jj, jx
      integer :: k, lr, n
      integer :: nstrtc
      integer :: numk
           ! INDEX of the preset lossData array closest to dlatr
      integer :: jlatmdx(numj)
      real*8  :: rjx
           ! column of loss rates at closest latitude location
      real*8  :: strtx (KDIM)
           ! pressure at edges of altitude grid (top down)
      real*8  :: p0l   (ivert+1)
      real*8  :: strt0l(ivert)
           ! loss rates at needed locations
      real*8  :: tlstt (numj, ivert, num_species)
!
!EOP
!------------------------------------------------------------------------------
!BOC

      jlatmdx(:) = 0

      strt0l (:) = 0.0d0

      numk   = ivert
      nstrtc = numk

!     ---------------------------------------------------------------------
!     The loss rates are known only at locations of 85S, 75S, by 10 degrees
!     to 85N.  The following loop is finding out which of these locations
!     is closest to the simulation latitude zone.  No interpolation is
!     done, just pick closest.
!     ---------------------------------------------------------------------

      do j = ju1_gl, j2_gl

        rjx = (dlatr(j) * DEGPRAD / 10.0d0) + 10.0d0

!       -----------------------------------------------------
!       Add a small delta to take care of precision problems.
!       -----------------------------------------------------
        rjx = rjx + 1.0d-06

        jx = rjx
        jj = j + 1 - ju1_gl

        jlatmdx(jj) = Min (JDIM, Max (1, jx))

      end do

!     -------------------------------------------
!     The following loop is flipping the presure,
!       i.e., flipping up/down to down/up,
!       e.g., (not actual values)
!             plevl(1)      = 1000
!             plevl(2)      =  800
!             . . .
!             plevl(numk+1) = ptop
!                       &
!             p0l(1)        = ptop
!             . . .
!             p0l(numk)     =  800
!             p0l(numk+1)   = 1000
!     -------------------------------------------

      do lr = 1, numk + 1
        p0l(lr) = plevl(numk + 2 - lr)
      end do

!     ----------------------------------------------------------------------
!     The following loops are calculating the loss rates.  One loops through
!     the number of tracers and latitudes for a given month and gets a
!     column of the preset loss rates at preset vertical locations.  Strt2m
!     does the interpolation to simulation vertical levels.
!     ----------------------------------------------------------------------

!     ---------------------------------------------------------
!     lossData defined for 18 lats (85S, 75S, ...85N) & 12 months.
!     Skip interpolating, pick nearest latitude.
!     ---------------------------------------------------------

      do n = 1, num_species
        do j = 1, numj

          jj = jlatmdx(j)

          do k = 1, KDIM
            strtx(k) = lossData(k,jj,month,n)
          end do

!         ===============
          call Strt2m_Gmi (strtx, KDIM, strt0l, p0l, nstrtc)
!         ===============

!         -------------------------------------------------
!         Store loss freq for exact CTM layers nstrtc down.
!         -------------------------------------------------

          do lr = 1, numk

            if (lr <= nstrtc) then
              tlstt(j,lr,n) = strt0l(lr)
            else
              tlstt(j,lr,n) = 0.0d0
            end if

          end do

        end do
      end do

      do n = 1, num_species
        do k = k1, k2
          do j = ju1, j2

            jj = j + 1 - ju1_gl

            do i = i1, i2

              loss_freq(i,j,k,n) = tlstt(jj,k2+1-k,n)

            end do

          end do
        end do
      end do

      return

      end subroutine Stratl_Gmi
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Strt2m_Gmi
!
! !INTERFACE:
!
      subroutine Strt2m_Gmi (strtx, nk, strt0l, p0l, nstrt)
!
      implicit none
!
! !INPUT PARAMETERS:
           ! number of vertical levels in stratospheric loss data
      integer, intent(in) :: nk
           ! number of vertical levels in the CTM
      integer, intent(in) :: nstrt
           ! pressure at edges of vertical grid (top down) (mb)
      real*8 , intent(in) :: p0l   (nstrt+1)
           ! column of loss rates at closest latitude location (s^-1)
      real*8 , intent(in) :: strtx (nk)
!
! !OUTPUT PARAMETERS:
           ! loss frequency at each level (s^-1)
      real*8 , intent(out) :: strt0l(nstrt)
!
! !DEFINED PARAMETERS:
      integer, parameter :: NPARAM = 27
!
! !DESCRIPTION:
!   This routine sets up std z* atmosphere:  p = 1000 * 10^(-z*/16 km).
!   Assume that stratospheric chemical parameters always start at
!   52 km (n=27).  Scan downward from 52 km to 14 km (nk=20) by 2 km
!   intervals, constant >52km.  N.B. F(@30km) assumed to be constant from
!   29-31 km (by mass).
!
! !DECLARED VARIABLES:
      integer :: i, k
      real*8  :: f0, f1, f2
      real*8  :: p1, p2
      real*8  :: xps0
      real*8  :: f  (NPARAM)
      real*8  :: ps (NPARAM+1)
      real*8  :: xps(NPARAM+1)
!
!EOP
!------------------------------------------------------------------------------
!BOC
!     ---------------------------------------------------------
!     These loops simply set up the grid for the known data for
!     interpolation.
!     ---------------------------------------------------------

      xps0   = 10.0d0**(-0.125d0)

      xps(1) = 1000.0d0

      do i = 2, NPARAM + 1
       xps(i) = xps(i-1) * xps0
      end do

      ps(1) = 1000.0d0

      do i = 2, NPARAM
        ps(i) = 0.5d0 * (xps(i-1) + xps(i))
      end do

      ps(NPARAM+1) = 0.0d0

      do i = 1, NPARAM - nk
        f(i) = 0.0d0
      end do

!     -----------------------------------
!     k=1 should be at top of atmosphere.
!     -----------------------------------

      do k = 1, nk
        f(NPARAM+1-k) = strtx(k)
      end do

!     ----------------------------------
!     Calling the interpolation routine.
!     ----------------------------------

      do k = 1, nstrt

        p1 = p0l(k+1)
        p2 = p0l(k)

!       ===============
        call Somlfq_Gmi (p1, p2, f0, f1, f2, ps, f, NPARAM)
!       ===============

        strt0l(k) = f0

      end do

      return

      end subroutine Strt2m_Gmi
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Somlfq_Gmi
!
! !INTERFACE:
!
      subroutine Somlfq_Gmi (p1, p2, f0, f1, f2, ps, f, n)
!
      implicit none
!
! !INPUT PARAMETERS:
           ! number of loss rates in column
      integer, intent(in) :: n       
           ! pressure at interface at bottom of grid box (mb)
      real*8 , intent(in) :: p1
           ! pressure at interface at top of grid box (mb)
      real*8 , intent(in) :: p2 
           ! column of loss rates (s^-1)
      real*8 , intent(in) :: f (n) 
           ! ! vertical grid of loss data (mb)
      real*8 , intent(in) :: ps(n+1) 
!
! !OUTPUT PARAMETERS:
           ! calculated loss frequency (s^-1)
      real*8 , intent(out) :: f0 
           ! higher order moment of loss frequency (unused)
      real*8 , intent(out) :: f1
           ! ext higher order moment of loss frequency (unused)
      real*8 , intent(out) :: f2 
!
! !DESCRIPTION:
!   This routine calculates the loss freq moments from a set of loss freq's
!   at std z*.
!
!   The MOMENTS for a square-wave or 'bar':  
! \begin{verbatim}
!        F(x) = f0  b<=x<=c,  = 0.0 else
!        S0 =   f0 (x)                      (from x=b to x=c)
!        S1 = 3 f0 (x^2 - x)                (from x=b to x=c)
!        S2 = 5 f0 (2x^3 - 3x^2 + x)        (from x=b to x=c)
! \end{verbatim}
!
! !LOCAL VARIABLES:
      integer :: i
      real*8  :: pc, pb
      real*8  :: xb, xc
!EOP
!------------------------------------------------------------------------------
!BOC
      f0 = 0.0d0
      f1 = 0.0d0
      f2 = 0.0d0

!     --------------------------------------------------------------------
!     Given a CTM model interval pressure range:  p1 > p2 (decreasing up),
!     the pressure levels between z* values are:
!                        ps(I) > ps(i+1) bounds z*(i)
!     n:  z* levels, ==> ps(n+1) = 0  (extrapolate chemical loss to top)
!       z1 = 16.0d0 * Log10 (1000.0d0 / p1)
!       z2 = 16.0d0 * Log10 (1000.0d0 / p2)
!     --------------------------------------------------------------------

      do i = 1, n

        pc = Min (p1, ps(i))
        pb = Max (p2, ps(i+1))

        if (pc > pb) then

!         -------------------------------------
!         Have condition:  p1 >= pc > pb >= p2.
!         -------------------------------------

          xc = (pc - p2) / (p1 - p2)
          xb = (pb - p2) / (p1 - p2)

!         ----------------------------------------------------------------
!         Have condition:  0 <= xb < xc <= 1.
!         Assume that the loss freq, f, is constant over interval (xb,xc),
!         f0: (c-b), f1: 6((c2-c)-(b2-b)), f2: 5((2c3-3c2+c)-(2b3-3b2+b)).
!         Calculate its contribution to the moments in the interval (0,1).
!         ----------------------------------------------------------------

          f0 = f0 + f(i) * (xc - xb)

          f1 = f1 + f(i) * 3.0d0 * ((xc * xc - xc) - (xb * xb - xb))

          f2 = f2 + f(i) * 5.0d0 *  &
     &         ((xc + xc - 1.0d0) * (xc * xc - xc) -  &
     &          (xb + xb - 1.0d0) * (xb * xb - xb))

        end if

      end do

!     ---------------------------------------------------------------
!     Restrain moments:  force monotonicity & positive at min end pt.
!     ---------------------------------------------------------------

      f0 = Max (f0, 0.0d0)

      if (f2 > 0.0d0) then

!       --------------------------------------------
!       Do not allow reversal of curvature:  f2 > 0.
!       --------------------------------------------

        f2 = Min (f2, Abs (f1) / 3.0d0, 0.5d0 * f0)

        if (f1 < 0.0d0) then
          f1 = Max (-(f0 + f2), f1)
        else
          f1 = Min (+(f0 + f2), f1)
        end if

      else

!       -----------------------------------------------
!       f2 < 0 = curved down at ends, allow if f1 < f0.
!       -----------------------------------------------

        f1 = Min (+f0, Max(-f0, f1))
        f2 = Max (f2, -f0 + Abs (f1), -Abs (f1) / 3.0d0)

      end if

      return

      end subroutine Somlfq_Gmi

!EOC
!------------------------------------------------------------------------------
      end module GmiStratosphericLoss_mod
