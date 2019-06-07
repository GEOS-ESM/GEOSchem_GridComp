!------------------------------------------------------------------
!  Array needed for relative time variation in solar flux at top of atmosphere
!   in NUMLAM wavelength bins, read cheecks that NUMLAM = NUM_SOLAR_CYC_LAM
!   iscyr is index in array based on current model year and month
!
!------------------------------------------------------------------
      module solar_cycle_mod

        integer, parameter ::  &
! number of wavelengths bins for the solar cycle dataset
     &   NUM_SOLAR_CYC_LAM = 79,  &
! number of months for the solar cycle dataset
     &   NUM_SOLAR_CYC_MON = 658

!... variable for relative solar flux in wavelength bins used in look up
!...  phot tables monthly data beginning Jan 1947
        real (kind=8), save ::  &
     &   s_cycle(NUM_SOLAR_CYC_LAM, NUM_SOLAR_CYC_MON)

      end module solar_cycle_mod

