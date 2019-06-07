      module gcr_mod
!------------------------------------------------------------------
!  Arrays and dimensions needed for calculation of the NOx emission
!   from Galactic Cosmic Rays from the sun. Depends on number of
!   sunspots and a linear fit to the emission.
!------------------------------------------------------------------

        integer, parameter ::  &
!... parameters in input file
     &    NUM_GCR_YEARS = 100  &
!     &  , NUM_GCR_LATS  = 91
!     &  , NUM_GCR_LEVS  = 33
     &  , NUM_GCR_LEVS = 100  &
     &  , NUM_GCR_LATS = 181

!... variables for calculation of the Galactic Cosmic Ray production
!...  of NOx (read in)
        real (kind=8), save ::  &
     &    gcr_date (NUM_GCR_YEARS)  &
     &  , gcr_sunspot (NUM_GCR_YEARS)

!... variables used in model
        real (kind=8), save ::  &
     &    gcr_slope (NUM_GCR_LATS, NUM_GCR_LEVS)  &
     &  , gcr_aintcp (NUM_GCR_LATS, NUM_GCR_LEVS)
!     &    gcr_slope(ju1:j2, k1:k2)
!     &  , gcr_aintcp(ju1:j2, k1:k2)

      end module gcr_mod
