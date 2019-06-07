      FUNCTION fallvel(rad,density,ndens,tem)

        USE MAPL_ConstantsMod,  only: MAPL_GRAV

!    Function fallvel ---------- David B. Considine, 11/11/96

!    Code to calculate the fall velocity for a particle of radius r,
!    which is passed to the function.  The fall velocity is calculated
!    according to Kasten [1968]. (Kasten, F., Falling Speed of Aerosol
!    Particles, J. Appl. Met., Oct, 1968

!    output

!    fallvel - speed of falling particle in cm/s

!    variables
!    
!    rad - input radius of falling particle (microns)
!    density - input density of particle (grams/cm**3)
!    ndens - atmospheric number density (molecules/cm**3)
!    dynvis - dynamic viscosity
!    tem - temperature (Kelvin)
!    mfp - mean free path
!    
!    parameters
!    
!    sigsq - square of effective collision diameter
!    bet -
!    s - sutherland's constant
!    a,b,cc - constants used in Kasten's formula
!    g - force of gravity (meters/sec**2)

!    declare variables
      implicit none
      real radius,rho,ndens,tem,mfp,dynvis
     &     ,fallvel,rad,density

!    set internal parameters
      real, parameter :: sigsq=1.3323e-19
      real, parameter :: bet=1.458e-6
      real, parameter :: s=110.4
      real, parameter :: a=1.249
      real, parameter :: b=0.42
      real, parameter :: cc=0.87

!     REAL, PARAMETER :: g = 9.8

!    don't cause problems if subroutine called with rad or density
!    equal to 0:
      IF(rad == 0.0 .OR. density == 0.0)THEN
         fallvel = 0.0
      ELSE

!    code assumes radius is in microns and rho is in grams/cm^3.
!    converts radius to meters and rho to kg/m^3:
         radius = rad*1.e-6
         rho = density*1.e3 

!    Calculate mean free path. According to the U.S. Standard Atmosphere,
!    1976, the mean free path is given by the relationship: mfp =
!    sqrt(2)/(2*pi*sig**2*N), where sig is the effective collision diameter
!    (3.65e-10 meters) and N is the number density (units of #/meter**3)
!    mfp units are 1/meters
         mfp = .22508/(sigsq*ndens*1.e6)

!    Calculate dynamic viscosity. The formula is:
!    dynvis=bet*T**(3/2)/(T+S), where bet is ..., s is Sutherland's
!    constant (110.4 Kelvin).  The units are Newtons/sec/meter^2.
         dynvis = bet*tem**1.5/(tem+s)

!    Calculate fall velocity.  The factor of 100 results in cm/s
!    output.
         fallvel = 0.2222*rho*radius*radius*MAPL_GRAV/dynvis*100.0*
     &                    (1.+ mfp/radius*(a+b*exp(-cc*radius/mfp)))
      END IF

      RETURN
      END FUNCTION fallvel
