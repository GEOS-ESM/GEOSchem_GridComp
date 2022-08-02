      module module_data_mosaic_main

      use module_data_mosaic_kind, only:  r8
      use module_data_mosaic_constants, only:  &
          avogad, deg2rad, pi, piover4, piover6, third


      implicit none

      integer, parameter ::   &
                ngas_com = 40,   &
                ngas_urb = 19,   &
                ngas_bio =  7,   &
                ngas_mar = 11
      !BSINGH - 05/28/2013(RCE updates)
      integer, parameter ::   &
                  naer_tot = 24 		      ! total num of 3-D variables per bin

      integer, save ::   &
                naerbin  = -999888777         ! number of bins (set at run time)
      !BSINGH - 05/28/2013(RCE updates ENDS)
!               naerbin  = 41760  	      ! ( 48 size)*(29 wbc)*(30 kappa)
!               naerbin  = 3240  	      ! ( 24 size)*(15 wbc)*( 9 kappa)
!               naerbin  = 90000 	      ! (100 size)*(30 wbc)*(30 kappa)

      integer, parameter ::   &
                ncld_tot = 13,		   &  ! + 8 = total num of 3-D variables/bin
                ncldbin  =  4,		   &  ! num of cloud bins
                ncld     = 22		! num of dynamic cloud species/bin

      integer, parameter :: ngas_max = ngas_com + ngas_urb + ngas_bio + ngas_mar
      
      integer, parameter :: ncld_max = ncld_tot*ncldbin
      !BSINGH - 05/28/2013(RCE updates)
      integer, save :: naer_max = -999888777  ! set at run time to naer_tot*naerbin

      integer, save :: ntot_max = -999888777  ! set at run time to (ngas_max + naer_max + ncld_max)
      !BSINGH - 05/28/2013(RCE updates ENDS)

      integer, save ::   &
      		naerbin_used=0,   &   ! num of aerosol bins being used
      		ncldbin_used=0,   &   ! num of  cloud  bins being used
      		ntot_used=ngas_max    ! portion of cnn array being used

      integer, save ::   &
      		ipmcmos = 0        ! if > 0, do emissions, dilution, air density,
      		                   ! and relative humidity as in partmc_mosaic 

      real(r8), parameter :: press0_pa = 1.01325d5  ! pressure of 1 atm [Pa]
      real(r8), parameter :: mw_air = 28.966d0      ! dry-air mean molecular weight [g/mol]

!------------------------------------------------------------------------
! Global Species Indices
!
      integer, save ::   &
       kh2so4,      khno3,       khcl,        knh3,        kno,   &
       kno2,        kno3,        kn2o5,       khono,       khno4,   &
       ko3,         ko1d,        ko3p,        koh,         kho2,   &
       kh2o2,       kco,         kso2,        kch4,        kc2h6,   &
       kch3o2,      kethp,       khcho,       kch3oh,      kanol,   &
       kch3ooh,     kethooh,     kald2,       khcooh,      krcooh,   &
       kc2o3,       kpan,   &
       karo1,       karo2,       kalk1,       kole1,       kapi1,   &
       kapi2,       klim1,       klim2,   &
       kpar,        kaone,       kmgly,       keth,        kolet,   &
       kolei,       ktol,        kxyl,        kcres,       kto2,   &
       kcro,        kopen,       konit,       krooh,       kro2,   &
       kano2,       knap,        kxo2,        kxpar,   &
       kisop,       kisoprd,     kisopp,      kisopn,      kisopo2,   &
       kapi,        klim,   &
       kdms,        kmsa,        kdmso,       kdmso2,      kch3so2h,   &
       kch3sch2oo,  kch3so2,     kch3so3,     kch3so2ch2oo,kch3so2oo,   &
       ksulfhox

      integer, save ::   &
       knum_a,      kdpdry_a,    ksigmag_a,  kjhyst_a,   &
       kwater_a,    kso4_a,      kno3_a,     kcl_a,       knh4_a,   &
       koc_a,       kmsa_a,      kco3_a,     kna_a,       kca_a,   &
       kbc_a,       koin_a,      karo1_a,    karo2_a,     kalk1_a,   &
       kole1_a,     kapi1_a,     kapi2_a,    klim1_a,     klim2_a

      integer, save ::   &
       knum_c,      kwater_c,    kso4_c,     kno3_c,     kcl_c,   &
       kmsa_c,      kco3_c,      knh4_c,     kna_c,      kca_c,   &
       koc_c,       kbc_c,       koin_c,   &
       karo1_c,     karo2_c,     kalk1_c,    kole1_c,    kapi1_c,   &
       kapi2_c,     klim1_c,     klim2_c



!-------------------------------------------------------------------------

      integer, save ::	m_partmc_mosaic  ! >0 for partmc_mosaic, <=0 for mosaic box model!BSINGH - 05/28/2013(RCE updates)

      integer, save ::	mgas, maer, mcld

      integer, save ::	maeroptic, mshellcore

      integer, save ::	msolar, mphoto


!------------------------------------------------------------------------

      end module module_data_mosaic_main
