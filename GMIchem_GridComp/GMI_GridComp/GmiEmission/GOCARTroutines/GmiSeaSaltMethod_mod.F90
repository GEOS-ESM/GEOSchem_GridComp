!-------------------------------------------------------------------------
!  NASA/GFSC, SIVO, Code 610.3
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: GmiSeaSaltMethod_mod
!
! !INTERFACE:
!
module GmiSeaSaltMethod_mod
!
  implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
!
  private 
  public  :: SourceSeaSalt
  public  :: InitializationSeaSalt
  public  :: Allocate_srcEmissSeaSalt
!
! !PUBLIC DATA MEMBERS:
!
  public  :: nSeaSaltBin
  public  :: srcEmissSeaSalt
!
  INTEGER, PARAMETER :: nSeaSaltBin = 4           ! number of sea salt bins
  INTEGER, PARAMETER :: msub        = 4
  REAL*8             :: ssaltden (nSeaSaltBin)    ! Density of sea salt (kg/m3)
  REAL*8             :: ssaltreff(nSeaSaltBin)    ! Main effective radius (m)
  REAL*8             :: ra       (nSeaSaltBin)
  REAL*8             :: rb       (nSeaSaltBin)
  REAL*8             :: ch_ss    (nSeaSaltBin,12)
  REAL*8 , POINTER   :: srcEmissSeaSalt(:,:,:) => null() ! Sea Salt emission computed from source
!
! !DESCRIPTION:
!  Provides routines to manipulate Sea Salt related variables.
!
! !AUTHOR:
!  Jules Kouatchou
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------------
CONTAINS

!-------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Allocate_srcEmissSeaSalt
!
! !INTERFACE:
!
  subroutine Allocate_srcEmissSeaSalt(i1, i2, ju1, j2)
!
  IMPLICIT NONE
!
! !INPUT PARAMETERS:
  integer, intent(in) :: i1, i2, ju1, j2
!
! !DESCRIPTION:
! Allocate the variable "srcEmissSeaSalt".
!
! !AUTHOR:
!  Jules Kouatchou
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC
  allocate (srcEmissSeaSalt (i1:i2, ju1:j2, nSeaSaltBin))
  srcEmissSeaSalt(:,:,:) = 0.0d0
  return
  end subroutine Allocate_srcEmissSeaSalt
!EOC
!-------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: InitializationSeaSalt
!
! !INTERFACE:
!
  subroutine InitializationSeaSalt(  )
!
  IMPLICIT NONE
!
!
! !DESCRIPTION:
! Initialize sea salt related variables: 
! ra, rb, ssaltden, ssaltreff, and ch_ss.
!
! !AUTHOR:
!  Jules Kouatchou
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC

  ! Density of sea salt (kg/m3)

  ssaltden(1:nSeaSaltBin) = 2200.0
     
  ! Main effective radius (m)

  ssaltreff(1) = 0.26E-6
  ssaltreff(2) = 1.19E-6
  ssaltreff(3) = 2.43E-6
  ssaltreff(4) = 7.57E-6

  ! Bin edges in um:

  ra(1) = 0.1
  ra(2) = 0.5
  ra(3) = 1.5
  ra(4) = 5.0

  rb(1) = 0.5
  rb(2) = 1.5
  rb(3) = 5.0
  rb(4) = 10.0
     
  ch_ss(:,:) =  1.0
 !         RESHAPE(SOURCE=ch_ss_aerocom_1(:),  SHAPE=(/4,12/)) / &
 !         RESHAPE(SOURCE=ch_ss_geos4_mod_1(:),SHAPE=(/4,12/)) 
     
  return

  end subroutine InitializationSeaSalt
!EOC
!-------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:
!
! !INTERFACE:
!
  SUBROUTINE SourceSeaSalt (mcor, tdt, nymd, i1, i2, ju1, j2, k1, k2)
!
  use GocartDerivedVariables_mod, only : w10m, airmas
  use GocartDerivedVariables_mod, only : lwi_flags => lwi_flags_1
  use GmiWaterMethod_mod        , only : water
  use GmiTimeControl_mod        , only : GmiSplitDateTime
!
  IMPLICIT NONE
!
# include "gmi_phys_constants.h"
!
! !INPUT PARAMETERS:
  INTEGER, INTENT(IN)    :: i1, i2, ju1, j2, k1, k2
  INTEGER, INTENT(IN)    :: nymd
  REAL*8 , INTENT(IN)    :: tdt
  REAL*8 , INTENT(IN)    :: mcor  (i1:i2, ju1:j2)         ! area of grid box (m2)
!
! !OUTPUT PARAMETERS:
!  REAL*8 , INTENT(OUT)   :: bems      (i1:i2, ju1:j2, nSeaSaltBin)
  REAL*8  bems      (i1:i2, ju1:j2, nSeaSaltBin)
                                                          ! dust tracer mixing ratio at the
                                                          ! first model layer
!
! !DESCRIPTION:
! Evaluate the source of each seasalt particles size classes  (kg/m3)
! by soil emission.
!
! Number flux density: Original formula by Monahan et al. (1986) adapted
! by Sunling Gong (JGR 1997 (old) and GBC 2003 (new)).  The new version is
! to better represent emission of sub-micron sea salt particles.
!
! \begin{eqnarray}
! \frac{dFn}{dr} = \frac{c1\times u10^{c2}}{r^A} \times (1+c3\times r^{c4})\times 10^{c5 \times \exp{-B^2}}
! \end{eqnarray}
! where
! \begin{eqnarray}
!  B = \frac{b1 - \log{(r)}}{b2}
! \end{eqnarray}
!
! see c_old, c_new, b_old, b_new below for the constants.
! number fluxes are at $80\%$ RH.
! 
! To calculate the flux:
! \begin{enumerate}
! \item Calculate dFn based on Monahan et al. (1986) and Gong (2003)
! \item Assume that wet radius r at 80% RH = dry radius r_d *frh
! \item Convert particles flux to mass flux :
! \end{enumerate}
! \begin{eqnarray}
! \begin{array}{lcl}
! \frac{dFM}{dr\_d} &=& \frac{4}{3}\times \pi\times rho\_d\times r\_d^3 \times \frac{dr}{dr\_d}\times \frac{dFn}{dr}\\
!           & = & \frac{4}{3}\times \pi\times rho_d\times r\_d^3 \times frh \times \frac{dFn}{dr}
! \end{array}
! \end{eqnarray}
! where $rho\_p$ is the particle density (kg/m3).
! The factor $10^{-18}$ is to convert in micro-meter $r\_d^3$.
!
! !LOCAL VARIABLES:
!  INTEGER :: lwi_flags (i1:i2, ju1:j2)     ! land, water, ice indicator
  REAL*8            :: c0(5), b0(2)
  REAL*8, PARAMETER :: c_old(5)=(/1.373, 3.41, 0.057, 1.05, 1.190/)
  REAL*8, PARAMETER :: c_new(5)=(/1.373, 3.41, 0.057, 3.45, 1.607/)
  REAL*8, PARAMETER :: b_old(2)=(/0.380, 0.650/)
  REAL*8, PARAMETER :: b_new(2)=(/0.433, 0.433/)
  REAL*8, PARAMETER :: dr=5.0E-2 ! um
  REAL*8, PARAMETER :: theta=30.0
  ! Swelling coefficient frh (d rwet / d rd)
!!!  REAL,    PARAMETER :: frh = 1.65
  REAL*8,    PARAMETER :: frh = 2.0
!  LOGICAL, PARAMETER :: old=.FALSE., new=.TRUE.
  LOGICAL, PARAMETER :: old=.TRUE., new=.FALSE.
  REAL*8 :: rho_d, r0, r1, r, r_w, a, b, dfn, r_d, dfm, src
  INTEGER :: i, j, n, nr, ir
  integer :: idumday, idumyear
  integer :: month_gmi
! 
! !AUTHOR:
!  Jules Kouatchou
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC

  srcEmissSeaSalt(:,:,:) = 0.0d0
  call GmiSplitDateTime (nymd, idumyear, month_gmi, idumday)

  DO n = 1,nSeaSaltBin
     bems(:,:,n) = 0.0
     rho_d = ssaltden(n)
     r0 = ra(n)*frh
     r1 = rb(n)*frh
     r = r0
     nr = INT((r1-r0)/dr)
     DO ir = 1,nr
        r_w = r + dr*0.5
        r = r + dr
        IF (new) THEN
           a = 4.7*(1.0 + theta*r_w)**(-0.017*r_w**(-1.44))
           c0 = c_new
           b0 = b_new
        ELSE
           a = 3.0
           c0 = c_old
           b0 = b_old
        END IF
        !
        b = (b0(1) - LOG10(r_w))/b0(2)
        dfn = (c0(1)/r_w**a)*(1.0 + c0(3)*r_w**c0(4))* &
             10**(c0(5)*EXP(-(b**2)))

        r_d = r_w/frh*1.0E-6  ! um -> m
        dfm = 4.0/3.0*GMI_PI*r_d**3*rho_d*frh*dfn*dr
!        dfm = 4.0/3.0*GMI_PI*r_d**3*rho_d*frh*dfn*dr*tdt
        DO i = i1,i2
           DO j = ju1,j2
              IF (water(i,j) > 0.0) THEN
!              IF (lwi_flags(i,j) == 0) THEN
                 src = dfm*mcor(i,j)*water(i,j)*w10m(i,j)**c0(2)
!                 src = dfm*mcor(i,j)*w10m(i,j)**c0(2)
!                 src = ch_ss(n,month_gmi)*dfm*mcor(i,j)*w10m(i,j)**c0(2)
                 srcEmissSeaSalt(i,j,n) = srcEmissSeaSalt(i,j,n) + src
!                 srcEmissSeaSalt(i,j,n) = srcEmissSeaSalt(i,j,n) + src/airmas(i,j,1)
              ELSE
                 src = 0.0
              END IF
              bems(i,j,n) = bems(i,j,n) + src
           END DO  ! i
        END DO ! j
     END DO ! ir
  END DO ! n

  return

  END SUBROUTINE SourceSeaSalt
!EOC
!-------------------------------------------------------------------------------

end module GmiSeaSaltMethod_mod
