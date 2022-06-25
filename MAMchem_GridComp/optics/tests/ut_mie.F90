program testMie
    
    use mie, only : mie_scattering_lognormal => scattering_lognormal, &
                    integration_method_midpoint, &
                    integration_method_simpson,  &
                    psd_number,                  &
                    psd_surface

    implicit none

    complex :: m      = (1.75, -4.5e-1) ! OPAC/soot refractive index
    real    :: sigma  = 2.00            ! standard deviation
    real    :: median = 0.0118d-6       ! radius, 1e-6 m
    
    real    :: r_min      = 0.005d-6    ! 'm'
    real    :: r_max      = 20.00d-6    ! 'm'
    real    :: wavelength = 0.500d-6    ! 'm'

    integer :: N          = 10000       ! number of intervals for integration
    logical :: mass_specific = .false.  ! if true return mass specific optical properties

    real    :: ext, sca, g

    ! reported values
    print *, "GADS/OPAC results for 'soot': (reported with three significant digits)"
    print *, 'extinction   = ', 6.385d-7,       'm-1, [N] = 1 particle cm-3'
    print *, 'scattering   = ', 1.441d-7,       'm-1, [N] = 1 particle cm-3'
    print *, 'ssa          = ', 2.257d-1
    print *, 'asimm. fact. = ', 0.353d0
    print *
    
    ! base-line (midpoint method with very large number of intervals): number distribution
    call mie_scattering_lognormal(psd_number, 2*median, sigma, m, wavelength, ext, sca, g, 2*r_min, 2*r_max, 100*N, mass_specific, integration_method_midpoint)
    print *, "Results using 'midpoint' method: N = ", 100*N
    print *, 'extinction   = ', ext      * 1d9, 'm-1, [N] = 1 particle cm-3'
    print *, 'scattering   = ', sca      * 1d9, 'm-1, [N] = 1 particle cm-3'
    print *, 'ssa          = ', sca/ext
    print *, 'asimm. fact. = ', g
    print *
    
    ! midpoint method
    call mie_scattering_lognormal(psd_number, 2*median, sigma, m, wavelength, ext, sca, g, 2*r_min, 2*r_max, N, mass_specific, integration_method_midpoint)
    print *, "Results using 'midpoint' method: N = ", N
    print *, 'extinction   = ', ext      * 1d9, 'm-1, [N] = 1 particle cm-3'
    print *, 'scattering   = ', sca      * 1d9, 'm-1, [N] = 1 particle cm-3'
    print *, 'ssa          = ', sca/ext  
    print *, 'asimm. fact. = ', g
    print *

    ! Simpson's  composite method
    call mie_scattering_lognormal(psd_number, 2*median, sigma, m, wavelength, ext, sca, g, 2*r_min, 2*r_max, N, mass_specific, integration_method_simpson)
    print *, "Results using Simpson's composite method: N = ", N
    print *, 'extinction   = ', ext      * 1d9, 'm-1, [N] = 1 particle cm-3'
    print *, 'scattering   = ', sca      * 1d9, 'm-1, [N] = 1 particle cm-3'
    print *, 'ssa          = ', sca/ext
    print *, 'asimm. fact. = ', g
    print *


    print *, 'TEST SURFACE SIZE DISTRIBUTION: MASS_SPECIFIC = False'
    mass_specific = .false.

    ! base-line (midpoint method with very large number of intervals): number distribution
    call mie_scattering_lognormal(psd_number, 2*median, sigma, m, wavelength, ext, sca, g, 2*r_min, 2*r_max, 100*N, mass_specific, integration_method_midpoint)
    print *, "Results using 'midpoint' method: N = ", 100*N
    print *, 'extinction   = ', ext      * 1d9, 'm-1, [N] = 1 particle cm-3'
    print *, 'scattering   = ', sca      * 1d9, 'm-1, [N] = 1 particle cm-3'
    print *, 'ssa          = ', sca/ext
    print *, 'asimm. fact. = ', g
    print *

    ! base-line (midpoint method with very large number of intervals): number distribution
    call mie_scattering_lognormal(psd_surface, 2*median*exp(2*log(sigma)**2), sigma, m, wavelength, ext, sca, g, 2*r_min, 2*r_max, 100*N, mass_specific, integration_method_midpoint)
    ext = (3.14159265358979d0 * (2 * median)**2 * exp(2 * log(sigma)**2)) * ext    ! normalize N = 1
    sca = (3.14159265358979d0 * (2 * median)**2 * exp(2 * log(sigma)**2)) * sca    ! ...

    print *, "Results using 'midpoint' method: N = ", 100*N
    print *, 'extinction   = ', ext      * 1d9, 'm-1, [N] = 1 particle cm-3'
    print *, 'scattering   = ', sca      * 1d9, 'm-1, [N] = 1 particle cm-3'
    print *, 'ssa          = ', sca/ext
    print *, 'asimm. fact. = ', g
    print *

    print *, 'TEST SURFACE SIZE DISTRIBUTION: MASS_SPECIFIC = True'
    mass_specific = .true.

    ! base-line (midpoint method with very large number of intervals): number distribution
    call mie_scattering_lognormal(psd_number, 2*median, sigma, m, wavelength, ext, sca, g, 2*r_min, 2*r_max, 100*N, mass_specific, integration_method_midpoint)
    print *, "Results using 'midpoint' method: N = ", 100*N
    print *, 'extinction   = ', ext      * 1d9, 'm-1, [N] = 1 particle cm-3'
    print *, 'scattering   = ', sca      * 1d9, 'm-1, [N] = 1 particle cm-3'
    print *, 'ssa          = ', sca/ext
    print *, 'asimm. fact. = ', g
    print *

    ! base-line (midpoint method with very large number of intervals): number distribution
    call mie_scattering_lognormal(psd_surface, 2*median*exp(2*log(sigma)**2), sigma, m, wavelength, ext, sca, g, 2*r_min, 2*r_max, 100*N, mass_specific, integration_method_midpoint)
    print *, "Results using 'midpoint' method: N = ", 100*N
    print *, 'extinction   = ', ext      * 1d9, 'm-1, [N] = 1 particle cm-3'
    print *, 'scattering   = ', sca      * 1d9, 'm-1, [N] = 1 particle cm-3'
    print *, 'ssa          = ', sca/ext
    print *, 'asimm. fact. = ', g
    print *



end program testMie
