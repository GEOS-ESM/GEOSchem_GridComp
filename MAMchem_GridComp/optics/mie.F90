module mie 

    use omp_lib

    implicit none

    private

    real, parameter :: pi = 3.141592653589793d0

    public scattering_lognormal

    integer, public, parameter :: integration_method_midpoint = 1  ! midpoint  method for numerical integration
    integer, public, parameter :: integration_method_simpson  = 2  ! Simpson's method for numerical integration
    
    integer, public, parameter :: psd_number  = 1                  ! particle size distrbution: number
    integer, public, parameter :: psd_surface = 2                  ! particle size distrbution: surface area
    integer, public, parameter :: psd_volume  = 3                  ! particle size distrbution: volume

    integer, public, parameter :: ERROR_UNSUPORTED_PSD = 100


contains

    subroutine scattering_lognormal(psd, Dg, sigma, refractive_index, wavelength, ext, sca, g, size_min, size_max, intervals, specific, method)
        
        !
        ! Returns (optionally specific) extinction and scattering, and asimmetry parameter of 
        ! a population of homogeneous spherical particles with lognormal size distribution
        ! and density equal to density of water.
        !
        ! Note: PDF of the size distribution is assumed to be normalized, therefore extinction and
        ! scattering computed with surface distribution parameters need to be multiplied by a
        ! factor (pi * (2 * Dgn)**2 * exp(2 * log(sigma)**2)) in order to compare them with the 
        ! results obtained with number distribution parameters. However there is no need to 
        ! apply this factor if mass specific extinction and scattering are computed, they can be
        ! directly compared regardless of the types of size distribution.
        !

        implicit none

        ! input
        integer, intent(in)  :: psd               ! particle syze distribution: number of surface
        real,    intent(in)  :: Dg                ! geometric mean (median) diameter of number or surface distribution
        real,    intent(in)  :: sigma             ! standard deviation of the lognormal size distribution
        complex, intent(in)  :: refractive_index  ! particle refractive index
        real,    intent(in)  :: wavelength        ! radiation wavelength
        real,    intent(in)  :: size_min          ! lower bound for integration over the size distribution
        real,    intent(in)  :: size_max          ! upper bound for integration over the size distribution
        integer, intent(in)  :: intervals         ! number of intervals/bins for unerical integration
        logical, intent(in)  :: specific          ! mass specific optical properties
        integer, intent(in)  :: method            ! integration method
!       real,    intent(in)  :: spike_threshold   ! MIEV0 specific, values below about 0.3 signify a ripple spike

        ! output
        real,    intent(out) :: ext               ! extinction
        real,    intent(out) :: sca               ! scattering
        real,    intent(out) :: g                 ! asimmetry factor

        ! local
        real, parameter :: density_water = 1000d0 ! density of water, 'kg m-3'

        select case (method)
            case (integration_method_simpson)
                call integrate_simpson_(psd, Dg, sigma, refractive_index, wavelength, ext, sca, g, size_min, size_max, intervals, specific)
            case default
                call integrate_midpoint_(psd, Dg, sigma, refractive_index, wavelength, ext, sca, g, size_min, size_max, intervals, specific)
        end select

    contains 
        subroutine integrate_midpoint_(psd, Dg, sigma, refractive_index, wavelength, ext, sca, g, size_min, size_max, intervals, specific)

            !
            ! Use midpoint rule to do integration of optical properties
            ! over the number size distribution.
            !
 
            implicit none 

            ! input
            integer, intent(in)  :: psd               ! particle syze distribution: number of surface
            real,    intent(in)  :: Dg                ! geometric mean (median) diameter of number of surface distribution
            real,    intent(in)  :: sigma             ! standard deviation of the lognormal size distribution
            complex, intent(in)  :: refractive_index  ! particle refractive index
            real,    intent(in)  :: wavelength        ! radiation wavelength
            real,    intent(in)  :: size_min          ! lower bound for integration over the size distribution
            real,    intent(in)  :: size_max          ! upper bound for integration over the size distribution
            integer, intent(in)  :: intervals         ! number of intervals/bins for unerical integration
            logical, intent(in)  :: specific          ! mass specific optical properties

            ! output
            real,    intent(out) :: ext               ! extinction
            real,    intent(out) :: sca               ! scattering
            real,    intent(out) :: g                 ! asimmetry factor

            ! local
            real    :: delta                          ! integration step
            real    :: d, s                           ! diameter, crossection
            real    :: n, m                           ! number, mass~volume
            real    :: dN, dS                         ! number of particles, volume of particles in a size bin
            integer :: i                              ! index
            real    :: q_ext, q_sca                   ! extinction amd scattering efficiencies
            real    :: gq_sca                         ! asimmetry factor times scattering efficiency
            real    :: spike                          ! MIEV0 spike indicator
            real    :: q_ext_, q_sca_, gq_sca_        ! last good values without resonance spike


            if ((psd /= psd_number) .and. (psd /= psd_surface)) then
                print *, "ERROR:: 'Unsuported PSD'       RC=", ERROR_UNSUPORTED_PSD
                stop ERROR_UNSUPORTED_PSD
            end if

            delta = (size_max - size_min)/intervals

            ext = 0.0d0
            sca = 0.0d0
            g   = 0.0d0

            m   = 0.0d0
            n   = 0.0d0

            q_ext_  = 0.0
            q_sca_  = 0.0
            gq_sca_ = 0.0

            NUMBER_DISTRIBUTION: if (psd == psd_number) then
                !$omp parallel default( none ) shared ( intervals, size_min, delta, refractive_index, wavelength, Dg, sigma, ext, sca, g, m ) private ( i, d, q_ext, q_sca, gq_sca, spike, dN, s ) 
                !$omp do reduction ( + : ext, sca, g, m ) &
                !$omp& schedule(guided)
                do i = 0, intervals - 1

                    d  = size_min + (i + 0.5d0)*delta

                    call scattering_sphere(d, refractive_index, wavelength, q_ext, q_sca, gq_sca, spike)

                    !if (spike > 0.3) then
                    !    q_ext  = q_ext_
                    !    q_sca  = q_sca_
                    !    gq_sca = gq_sca_
                    !else
                    !    print *, 'mie: spike resonance detected in size integration... [excluded]'
                    !end if

                    dN  = psd_lognormal(d, Dg, sigma) * delta   ! number of particles in (d, d+delta) size interval
                    s   = (pi / 4d0) * d**2                     ! cross-section 

                    ext = ext + (s * q_ext)  * dN               ! integrated size dependent optical properties
                    sca = sca + (s * q_sca)  * dN               ! ...
                    g   = g   + (s * gq_sca) * dN               ! ...

                    m   = m + (pi / 6d0) * d**3 * dN            ! population volume 
                end do
                !$omp end do
                !$omp end parallel
            end if NUMBER_DISTRIBUTION

            SURFACE_DISTRIBUTION: if (psd == psd_surface) then
                !$omp parallel default( none ) shared ( intervals, size_min, delta, refractive_index, wavelength, Dg, sigma, ext, sca, g, m) private ( i, d, q_ext, q_sca, gq_sca, spike, dS ) 
                !$omp do reduction ( + : ext, sca, g, m ) &
                !$omp& schedule(guided)
                do i = 0, intervals - 1
                    d  = size_min + (i + 0.5d0)*delta

                    call scattering_sphere(d, refractive_index, wavelength, q_ext, q_sca, gq_sca, spike)

                    !if (spike > 0.3) then
                    !    q_ext  = q_ext_
                    !    q_sca  = q_sca_
                    !    gq_sca = gq_sca_
                    !else
                    !    print *, 'mie: spike resonance detected in size integration... [excluded]'
                    !end if

                    dS  = psd_lognormal(d, Dg, sigma) * delta   ! surface of particles in (d, d+delta) size interval

                    ext = ext + q_ext  * dS                     ! integrated size dependent optical properties
                    sca = sca + q_sca  * dS                     ! ...
                    g   = g   + gq_sca * dS                     ! ...

                    m   = m + (1d0 / 6d0) * d * dS              ! population volume
                end do
                !$omp end do
                !$omp end parallel

                ext = 0.25d0 * ext                              ! surface factor = 1/4
                sca = 0.25d0 * sca                              ! ...
                g   = 0.25d0 * g                                ! ...

            end if SURFACE_DISTRIBUTION

            if (specific) then
                m = density_water * m                           ! population mass
            else
                m = 1.0d0
            end if

            g   = g   / sca
            ext = ext / m
            sca = sca / m

        end subroutine integrate_midpoint_


        subroutine integrate_simpson_(psd, Dg, sigma, refractive_index, wavelength, ext, sca, g, size_min, size_max, intervals, specific)

            !
            ! Use midpoint rule to do integration of optical properties
            ! over the number size distribution.
            !
            
            implicit none 

            ! input
            integer, intent(in)  :: psd               ! particle syze distribution: number of surface
            real,    intent(in)  :: Dg                ! geometric mean (median) diameter of number size distribution
            real,    intent(in)  :: sigma             ! standard deviation of the lognormal size distribution
            complex, intent(in)  :: refractive_index  ! particle refractive index
            real,    intent(in)  :: wavelength        ! radiation wavelength
            real,    intent(in)  :: size_min          ! lower bound for integration over the size distribution
            real,    intent(in)  :: size_max          ! upper bound for integration over the size distribution
            integer, intent(in)  :: intervals         ! number of intervals/bins for unerical integration
            logical, intent(in)  :: specific          ! mass specific optical properties

            ! output
            real,    intent(out) :: ext               ! extinction
            real,    intent(out) :: sca               ! scattering
            real,    intent(out) :: g                 ! asimmetry factor


            ! local
            real    :: delta                          ! integration step
            real    :: w                              ! weight factor
            real    :: d                              ! diameter
            real    :: s, m, dN                       ! crossection, mass|volume, number of particles
            integer :: i, n                           ! index
            real    :: q_ext, q_sca                   ! extinction amd scattering efficiencies
            real    :: gq_sca                         ! asimmetry factor times scattering efficiency
            real    :: spike                          ! MIEV0 spike indicator


            if ((psd /= psd_number)) then
                stop ERROR_UNSUPORTED_PSD
            end if

            ! number if intervals has to be even
            if ((intervals/2)*2 == intervals) then
              n = intervals
            else
              n = intervals + 1
            end if

            delta = (size_max - size_min)/n

            ext = 0.0d0
            sca = 0.0d0
            g   = 0.0d0

            m   = 0.0d0

            ! step 1
            d = size_min
            call scattering_sphere(d, refractive_index, wavelength, q_ext, q_sca, gq_sca, spike)

            dN  = psd_lognormal(d, Dg, sigma) * delta   ! number of particles in (d, d+delta) size interval
            s   = (pi / 4d0) * d**2                     ! cross-section 

            ext = ext + (s * q_ext)  * dN               ! integrated size dependent optical properties
            sca = sca + (s * q_sca)  * dN               !
            g   = g   + (s * gq_sca) * dN               !

            m  = m + (pi / 6d0) * d**3 * dN             ! population volume

            
            ! step 2
            d = size_max
            call scattering_sphere(d, refractive_index, wavelength, q_ext, q_sca, gq_sca, spike)

            dN  = psd_lognormal(d, Dg, sigma) * delta   ! number of particles in (d, d+delta) size interval
            s   = (pi / 4d0) * d**2                     ! cross-section 

            ext = ext + (s * q_ext)  * dN               ! integrated size dependent optical properties
            sca = sca + (s * q_sca)  * dN               !
            g   = g   + (s * gq_sca) * dN               !

            m  = m + (pi / 6d0) * d**3 * dN             ! population volume


            ! step 3
            w = 4.0d0
            do i = 1, (n - 1), 2
                d  = size_min + i*delta

                call scattering_sphere(d, refractive_index, wavelength, q_ext, q_sca, gq_sca, spike)

                dN  = psd_lognormal(d, Dg, sigma) * delta   ! number of particles in (d, d+delta) size interval
                s   = (pi / 4d0) * d**2                     ! cross-section 

                ext = ext + w * (s * q_ext)  * dN           ! integrated size dependent optical properties
                sca = sca + w * (s * q_sca)  * dN           !
                g   = g   + w * (s * gq_sca) * dN           !

                m   = m   + w * (pi / 6d0) * d**3 * dN      ! population volume
            end do

            
            ! step 4
            w = 2.0d0
            do i = 2, (n - 2), 2
                d  = size_min + i*delta

                call scattering_sphere(d, refractive_index, wavelength, q_ext, q_sca, gq_sca, spike)

                dN  = psd_lognormal(d, Dg, sigma) * delta   ! number of particles in (d, d+delta) size interval
                s   = (pi / 4d0) * d**2                     ! cross-section 

                ext = ext + w * (s * q_ext)  * dN           ! integrated size dependent optical properties
                sca = sca + w * (s * q_sca)  * dN           !
                g   = g   + w * (s * gq_sca) * dN           !

                m   = m   + w * (pi / 6d0) * d**3 * dN      ! population volume
            end do


            ! step 5
            ext = ext / 3.0d0
            sca = sca / 3.0d0
            g   = g   / 3.0d0
            m   = m   / 3.0d0

           
            ! return optics
            if (specific) then
                m = density_water * m                       ! population mass
            else
                m = 1.0d0
            end if

            g   = g   / sca
            ext = ext / m
            sca = sca / m

        end subroutine integrate_simpson_

    end subroutine scattering_lognormal

#if (0)
    subroutine scattering_sphere(diameter, refractive_index, wavelength, q_ext, q_sca, gq_sca, spike)

        !
        ! Returns extinction and scattering efficiencies, and asimmetry parameter of 
        ! a homogeneous spherical particle at some wavelegth.
        !
        ! This subroutine provides simplified interface to the bhmie code.
        ! 

        implicit none

        ! input
        real,    intent(in)  :: diameter          ! sphere diameter
        complex, intent(in)  :: refractive_index  ! sphere refractive index 
        real,    intent(in)  :: wavelength        ! radiation wavelength 

        ! output
        real,    intent(out) :: q_ext             ! extinction efficiency
        real,    intent(out) :: q_sca             ! scattering efficiency
        real,    intent(out) :: gq_sca            ! asymmetry parameter (g) times scattering efficiency (q_sca)
        real,    intent(out) :: spike             ! defaults to 1.0 when not on a spike.
                                                  ! values below about 0.3 signify a ripple spike.
    
        ! local
        real    :: x          ! size paramter 

        integer, parameter :: MXNANG=1000         ! see BHMIEV/README
        complex            :: S1(2*MXNANG-1)
        complex            :: S2(2*MXNANG-1)

        integer            :: n_ang
        real               :: q_back

        ! initialize some of the BHMIE inputs
        
        n_ang  = 2

        q_ext  = 0.0
        q_sca  = 0.0
        q_back = 0.0
        gq_sca = 0.0

        x = size_parameter(diameter, wavelength)

        call bhmie(x, refractive_index, n_ang, S1, S2, q_ext, q_sca, q_back, gq_sca)
        gq_sca = gq_sca*q_sca

        spike  = 1.0

    end subroutine scattering_sphere
#else
    subroutine scattering_sphere(diameter, refractive_index, wavelength, q_ext, q_sca, gq_sca, spike)

        !
        ! Returns extinction and scattering efficiencies, and asimmetry parameter of 
        ! a homogeneous spherical particle at some wavelegth.
        !
        ! This subroutine provides simplified interface to the MIEV0 code.
        ! 

        implicit none

        ! input
        real,    intent(in)  :: diameter          ! sphere diameter
        complex, intent(in)  :: refractive_index  ! sphere refractive index 
        real,    intent(in)  :: wavelength        ! radiation wavelength 

        ! output
        real,    intent(out) :: q_ext             ! extinction efficiency
        real,    intent(out) :: q_sca             ! scattering efficiency
        real,    intent(out) :: gq_sca            ! asymmetry parameter (g) times scattering efficiency (q_sca)
        real,    intent(out) :: spike             ! defaults to 1.0 when not on a spike.
                                                  ! values below about 0.3 signify a ripple spike.
    
        ! local
        real    :: x          ! size paramter 

        logical :: perfct     ! see MIEV0.doc   
        real    :: mimcut     !
        logical :: anyang     !
        integer :: numang     !
        real    :: xmu(2)     !
        integer :: nmom       !
        integer :: ipolzn     !
        integer :: momdim     ! determines first dimension of PMOM
        logical :: prnt(2)    !
        real    :: p_mom(0:1,4) !
        complex :: sforw      !
        complex :: sback      !
        complex :: s1(2)      !
        complex :: s2(2)      !
        complex :: t_forw(2)  !
        complex :: t_back(2)  !


        ! initialize some of the MIEV0 inputs
        perfct = .false.
        mimcut = 1.0d-6
        anyang = .false.
        numang = 0            ! skip calculation of S1 and S2
        xmu    = 0
        nmom   = 0            ! prevent calculation of PMOM
        ipolzn = 0
        momdim = 1
        prnt   = .false.
        p_mom  = 0.0
        s1     = (0.0, 0.0)
        s2     = (0.0, 0.0)
        t_forw = (0.0, 0.0)
        t_back = (0.0, 0.0)

        q_ext  = 0.0
        q_sca  = 0.0
        gq_sca = 0.0
        spike  = 0.0
        
        x = size_parameter(diameter, wavelength)

        call miev0(x, refractive_index, perfct, mimcut, anyang,      &
                   numang, xmu, nmom, ipolzn, momdim, prnt, q_ext, q_sca, gq_sca, &
                   p_mom, sforw, sback, s1, s2, t_forw, t_back, spike)
        
    end subroutine scattering_sphere
#endif

    function size_parameter(diameter, wavelength)

        !
        ! Returns the size parameter of a spherical
        ! particle at some wavelegth.
        !

        implicit none

        real :: size_parameter            ! result
        
        ! input
        real, intent(in) :: diameter      ! sphere size/diameter
        real, intent(in) :: wavelength    ! radiation wavelength

        size_parameter = pi * diameter / wavelength

    end function size_parameter


    function psd_lognormal(x, median, sigma)

        !
        ! Computes lognormal particle size distribution.
        !

        implicit none

        real :: psd_lognormal             ! result
        
        ! input
        real, intent(in) :: x             ! parameter
        real, intent(in) :: median        ! geometric median
        real, intent(in) :: sigma         ! standard deviation

        ! local
        real :: log_sigma

        log_sigma = log(sigma)

        psd_lognormal = 1/(x * log_sigma * sqrt(2*pi)) * exp(-(log(x) - log(median))**2 / (2*log_sigma*log_sigma))

    end function psd_lognormal

end module mie
