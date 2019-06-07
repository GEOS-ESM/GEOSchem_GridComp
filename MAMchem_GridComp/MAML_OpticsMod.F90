#include "MAPL_Exceptions.h"

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1    !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  MAML_OpticsMod--- aerosol optics.
!
! !INTERFACE:
!

   module  MAML_OpticsMod

! !USES:

   use ESMF
   use MAPL_Mod

   use MAPL_SimpleBundleMod

   use MAM_BaseMod
   use MAM_ConstituentsDataMod
   use MAM_ComponentsDataMod

   use MAML_OpticsTableMod, only: MAML_OpticsTable

   implicit none


! !PUBLIC TYPES:

   private

!
! !PUBLIC MEMBER FUNCTIONS:
!
   public  MAML_OpticsInterpolate     ! Computes aerosol optical properties

!
! !GLOBAL PARAMETERS
!

   integer,      parameter :: f = MAPL_R8

!
! !DESCRIPTION:
!
!  This module computes aerosol optical properties.
!
! !REVISION HISTORY:
!
!  09May2013 Darmenov - Initial code.
!
!EOP
!-------------------------------------------------------------------------

contains

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  MAML_OpticsCalculator --- Compute aerosol ext, sca, asy...
!
! !INTERFACE:
!

   subroutine MAML_OpticsInterpolate(mie_table, band, q, density, dgn_wet, delp, ext, sca, asy, nc, i1, i2, j1, j2, k1, k2, rc)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

! !INPUT PARAMETERS:
   integer, intent(in)                               :: nc             ! number of aerosol components

   integer, intent(in)                               :: i1, i2         ! upper and lower bounds of spatial dims
   integer, intent(in)                               :: j1, j2         ! ...
   integer, intent(in)                               :: k1, k2         ! ...

   type(MAML_OpticsTable), intent(in)                :: mie_table
   integer, intent(in)                               :: band           ! band index
   real, dimension(nc,i1:i2,j1:j2,k1:k2), intent(in) :: q              ! mass mixing ratios of aerosol components (including absorbed water), 'kg kg-1'
   real, dimension(nc), intent(in)                   :: density        ! density of aerosol components, 'kg m-3'
   real, dimension(i1:i2,j1:j2,k1:k2), intent(in)    :: dgn_wet        ! wet size of number size distribution

   real, dimension(i1:i2,j1:j2,k1:k2), intent(in)    :: delp           ! 

! !OUTPUT PARAMETERS:
   real, dimension(i1:i2,j1:j2,k1:k2), intent(out)   :: ext, sca, asy  ! extinction, scattering and asymmetry parameter
   integer, optional, intent(out)                    :: rc             ! Error return code:
                                                                       !  0 - all is well
                                                                       !  1 -

! !DESCRIPTION:
!  Interpolates Mie lookup table and returns extinction, scattering and asymmetry parameter.
!  NOTE: It is important that the q and density arrays follow the same order of components 
!        as the mie_table.
!
! !REVISION HISTORY:
!
!  09May2013 Darmenov    API.
!
!EOP
!-------------------------------------------------------------------------

                           __Iam__("MAML_OpticsInterpolate")

   real(kind=f) :: Dgs_min, Dgs_max
   real(kind=f) :: logDgs_min, logDgs_max, logDgs
   real(kind=f) :: sigma, log_sigma
   real(kind=f) :: n_re_min, n_re_max
   real(kind=f) :: n_im_min, n_im_max
   real(kind=f) :: q_max

   real(kind=f) :: vol, vol_
   real(kind=f) :: n_re_mix, n_im_mix, n_re_mix_, n_im_mix_

   real(kind=f) :: x                   ! chebyshev polynomial parameter
   real(kind=f) :: tt, uu              ! bilinear interpolation and weights
   real(kind=f) :: w00, w01, w10, w11  ! ...

   real(kind=f), allocatable, dimension(:) :: T
   real(kind=f), allocatable, dimension(:) :: c_ext
   real(kind=f), allocatable, dimension(:) :: c_sca
   real(kind=f), allocatable, dimension(:) :: c_asy

   real(kind=f), allocatable, dimension(:) :: n_re
   real(kind=f), allocatable, dimension(:) :: n_im

   real(kind=f), allocatable, dimension(:) :: component_n_re
   real(kind=f), allocatable, dimension(:) :: component_n_im
   real(kind=f), allocatable, dimension(:) :: q_

   real(kind=f), allocatable, dimension(:,:,:) :: mie_c_ext
   real(kind=f), allocatable, dimension(:,:,:) :: mie_c_sca
   real(kind=f), allocatable, dimension(:,:,:) :: mie_c_asy

   integer :: n_ac            ! number of aerosol components

   integer :: i_re, i_im      ! indexes

   integer :: i, j, k, n      ! loop counters

   integer :: size_im, size_re, size_n

   real(kind=f)                            :: factor                     ! multiplication factor
   real(kind=f), parameter                 :: density_water = 1000.0_f   ! kg m-3


   if (present(rc)) rc = MAM_NOT_IMPLEMENTED_ERROR


   ! basic dim-size check
   ASSERT_(mie_table%n_aerosol_components == nc)


   ! short hand and likely to improve data storage
   n_ac    = mie_table%n_aerosol_components
   Dgs_min = mie_table%Dgs_min
   Dgs_max = mie_table%Dgs_max
   sigma   = mie_table%sigma

   allocate(component_n_re(n_ac), __STAT__)
   allocate(component_n_im(n_ac), __STAT__)
   allocate(q_(n_ac), __STAT__)

   allocate(n_re(mie_table%n_refractive_index_re), __STAT__)
   allocate(n_im(mie_table%n_refractive_index_im), __STAT__)

   component_n_re = mie_table%component_refractive_index_re(band,:)
   component_n_im = mie_table%component_refractive_index_im(band,:)

   n_re     = mie_table%refractive_index_re(band, :)
   n_im     = mie_table%refractive_index_im(band, :)

   n_re_min = n_re(1) ! values are monotonically increasing
   n_re_max = n_re(mie_table%n_refractive_index_re)

   n_im_min = n_im(1) ! values are monotonically increasing
   n_im_max = n_im(mie_table%n_refractive_index_im)

   allocate(    T(mie_table%n_cheb), __STAT__)
   allocate(c_ext(mie_table%n_cheb), __STAT__)
   allocate(c_sca(mie_table%n_cheb), __STAT__)
   allocate(c_asy(mie_table%n_cheb), __STAT__)

   size_im = size(mie_table%c_ext,dim=2)
   size_re = size(mie_table%c_ext,dim=3)
   size_n  = size(mie_table%c_ext,dim=4)

   allocate(mie_c_ext(size_n,size_im,size_re), __STAT__)
   allocate(mie_c_sca(size_n,size_im,size_re), __STAT__)
   allocate(mie_c_asy(size_n,size_im,size_re), __STAT__)

   do k = 1, size_re
      do j = 1, size_im
         do i = 1, size_n
            mie_c_ext(i,j,k) = mie_table%c_ext(band,j,k,i)
            mie_c_sca(i,j,k) = mie_table%c_sca(band,j,k,i)
            mie_c_asy(i,j,k) = mie_table%c_asy(band,j,k,i)
         end do
      end do
   end do

   logDgs_min = log(Dgs_min)
   logDgs_max = log(Dgs_max)
   log_sigma  = log(sigma)

   ext = 0.0
   sca = 0.0
   asy = 0.0

   do k = k1, k2
       do j = j1, j2
           do i = i1, i2

               if (any(q(:,i,j,k) > 0.0)) then
                   ! normalize mass mixing ratios to avoid numerical issues 
                   q_max = maxval(q(:,i,j,k))
                   q_    = q(:,i,j,k) / q_max
                   vol   = sum(q_/density)

                   ! compute effective refractive index using volume mixing rule
                   n_re_mix = sum(component_n_re * q_/density) / vol
                   n_im_mix = sum(component_n_im * q_/density) / vol

                   vol = q_max * vol
               else
                   q_max = 0.0_f
                   vol   = 0.0_f
                   n_re_mix = component_n_re(1)
                   n_im_mix = component_n_im(1)
               end if   


               if (n_re_mix < n_re_min) then
                   n_re_mix = n_re_min + tiny(0.0)
               end if

               if (n_re_mix > n_re_max) then
                   n_re_mix = n_re_max - tiny(0.0)
               end if    
 
               if (n_im_mix < n_im_min) then
                   n_im_mix = n_im_min + tiny(0.0)
               end if

               if (n_im_mix > n_im_max) then
                   n_im_mix = n_im_max - tiny(0.0)
               end if

               i_re = locate(n_re_mix, n_re, mie_table%n_refractive_index_re)
               i_im = locate(n_im_mix, n_im, mie_table%n_refractive_index_im)

               if (i_re < 1                                 ) i_re = 1
               if (i_re > mie_table%n_refractive_index_re -1) i_re = mie_table%n_refractive_index_re -1
               if (i_im < 1                                 ) i_im = 1
               if (i_im > mie_table%n_refractive_index_im -1) i_im = mie_table%n_refractive_index_im -1


! bilinear interpolation: notes
!              t = (x - x1) / (x2 - x1)
!              u = (y - y1) / (y2 - y1)
!
!              w00 = (1.0_f - t) * (1.0_f - u), f00 = f1
!              w10 = (        t) * (1.0_f - u), f10 = f2
!              w01 = (1.0_f - t) * (        u), f01 = f4
!              w11 = (        t) * (        u), f11 = f3
!
!              bilinear_interpolation = (1.0_f - t) * (1.0_f - u) * f1 + &
!                                       (        t) * (1.0_f - u) * f2 + &
!                                       (        t) * (        u) * f3 + &
!                                       (1.0_f - t) * (        u) * f4
!
!              bilinear_interpolation = w00*f00 + w10*f10 + w01*f01 + w11*f11
!
!
!              mapping: 
!                       x -> re
!                       y -> im 

               tt = (    n_re_mix - n_re(i_re)) / &
                    (n_re(i_re+1) - n_re(i_re))

               uu = (    n_im_mix - n_im(i_im)) / &
                    (n_im(i_im+1) - n_im(i_im))


               ASSERT_ (tt >= 0 .and. tt <= 1.0)
               ASSERT_ (uu >= 0 .and. uu <= 1.0)

               w00 = (1.0_f - tt) * (1.0_f - uu)     ! weight for f(i_re+0, i_im+0)
               w10 = (        tt) * (1.0_f - uu)     ! weight for f(i_re+1, i_im+0)
               w01 = (1.0_f - tt) * (        uu)     ! weight for f(i_re+0, i_im+1)
               w11 = (        tt) * (        uu)     ! weight for f(i_re+1, i_im+1)

               T     = 0.0_f
               c_ext = 0.0_f
               c_sca = 0.0_f
               c_asy = 0.0_f

               do n = 1, mie_table%n_cheb
                   c_ext(n) = w00 * mie_c_ext(n, i_im  , i_re  ) + &
                              w10 * mie_c_ext(n, i_im  , i_re+1) + &
                              w01 * mie_c_ext(n, i_im+1, i_re  ) + &
                              w11 * mie_c_ext(n, i_im+1, i_re+1)
                   
                   c_sca(n) = w00 * mie_c_sca(n, i_im  , i_re  ) + &
                              w10 * mie_c_sca(n, i_im  , i_re+1) + &
                              w01 * mie_c_sca(n, i_im+1, i_re  ) + &
                              w11 * mie_c_sca(n, i_im+1, i_re+1)

                   c_asy(n) = w00 * mie_c_asy(n, i_im  , i_re  ) + &
                              w10 * mie_c_asy(n, i_im  , i_re+1) + &
                              w01 * mie_c_asy(n, i_im+1, i_re  ) + &
                              w11 * mie_c_asy(n, i_im+1, i_re+1)
               end do

               logDgs = log(dgn_wet(i,j,k)) + 2.0_f*log_sigma*log_sigma
               x = (2.0_f*logDgs - logDgs_max - logDgs_min) / (logDgs_max - logDgs_min)

               if (x > 1.0_f) then
                   x = 1.0_f - tiny(0.0_f)
               end if

               if (x < -1.0_f) then
                   x = -1.0_f + tiny(0.0_f)
               end if


               ! compute Chebyshev polynomials
               T = cheb_poly(mie_table%n_cheb, x)

               ! extinction and scattering are unitless, i.e., AOT(i,j) = sum(ext(i,j,1:km)
               factor = vol * density_water * delp(i,j,k) / MAPL_GRAV

               ext(i,j,k) = sum(c_ext * T) * factor
               sca(i,j,k) = sum(c_sca * T) * factor
               asy(i,j,k) = sum(c_asy * T)

               if (ext(i,j,k) < 0.0_f) then
                   ext(i,j,k) = tiny(0.0_f)
               end if

               if (sca(i,j,k) > ext(i,j,k)) then
                   sca(i,j,k) = ext(i,j,k)
               end if

           end do
       end do
   end do
 
   deallocate(component_n_re, component_n_im, __STAT__)
   deallocate(n_re, n_im, q_,                 __STAT__)
   deallocate(T, c_ext, c_sca, c_asy,         __STAT__)
   deallocate(mie_c_ext, mie_c_sca, mie_c_asy,__STAT__)


   RETURN_(ESMF_SUCCESS)

  end subroutine MAML_OpticsInterpolate



!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  MAML_RefractiveIndex --- Computes refractive index of internaly mixed components.
!
! !INTERFACE:
!

   subroutine MAML_RefractiveIndex(q_mass,   &
                                   density,  &
                                   n_re,     &
                                   n_im,     &
                                   n_re_mix, &
                                   n_im_mix, &
                                   rc)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

! !INPUT PARAMETERS:
   real, dimension(:), intent(in) :: q_mass      ! mass mixing ratios of components, including absorbed water
   real, dimension(:), intent(in) :: density     ! densities of components
   real, dimension(:), intent(in) :: n_re        ! real      part of refractive index 
   real, dimension(:), intent(in) :: n_im        ! imaginary part of refractive index


! !OUTPUT PARAMETERS:
   real, intent(out)              :: n_re_mix    ! real      part of effective refractive index
   real, intent(out)              :: n_im_mix    ! imaginary part of effective refractive index
   integer, optional, intent(out) :: rc          ! Error return code:
                                                 !  0 - all is well
                                                 !  1 - 

! !DESCRIPTION: Computes refractive index of internally mixed components using 
!               volume mixing rule.
!
!
! !REVISION HISTORY:
!
!  08May2014 Darmenov    API.
!
!EOP
!-------------------------------------------------------------------------

                           __Iam__("MAML_RefractiveIndex")

   ! local variables
   real :: vol         ! volume


   if (present(rc)) rc = MAM_NOT_IMPLEMENTED_ERROR


   vol = sum(q_mass / density)

   n_re_mix = sum(n_re_mix * q_mass / density) / vol
   n_im_mix = sum(n_im_mix * q_mass / density) / vol


   RETURN_(ESMF_SUCCESS)

  end subroutine MAML_RefractiveIndex




function locate(v, x, nx)

    !
    ! Returns the index of the largest element from the array x
    ! that is smaller than v.
    !
    ! Assumes that the array x is sorted in increasing order,
    ! and that min(x) <= v <= max(x)
    !

    implicit none

    integer :: locate

    real(kind=f), intent(in)                :: v
    real(kind=f), dimension(nx), intent(in) :: x
    integer, intent(in)                     :: nx

    ! local
    integer i
   
    do i = 1, nx
        if (x(i) > v) exit         
    end do
    
    locate = (i - 1)

    return
end function locate


function linear_interpolation(x1, x2, f1, f2, x)

    !
    ! Linear interpolation.
    !

    implicit none

    real(kind=f) :: linear_interpolation

    real(kind=f), intent(in) :: x1
    real(kind=f), intent(in) :: x2
    real(kind=f), intent(in) :: f1
    real(kind=f), intent(in) :: f2
    real(kind=f), intent(in) :: x

    linear_interpolation = f1 + ((x - x1) / (x2 - x1)) * (f2 - f1)
end function linear_interpolation


function bilinear_interpolation(x1, x2, y1, y2, f1, f2, f3, f4, x, y)
    
    ! 
    ! Bilinear interpolation.
    ! 

    implicit none

    real(kind=f) :: bilinear_interpolation

    real(kind=f), intent(in) :: x1
    real(kind=f), intent(in) :: x2
    real(kind=f), intent(in) :: y1
    real(kind=f), intent(in) :: y2
    real(kind=f), intent(in) :: f1
    real(kind=f), intent(in) :: f2
    real(kind=f), intent(in) :: f3
    real(kind=f), intent(in) :: f4
    real(kind=f), intent(in) :: x
    real(kind=f), intent(in) :: y

    ! local
    real(kind=f) :: t, u

    t = (x - x1) / (x2 - x1)
    u = (y - y1) / (y2 - y1)

    bilinear_interpolation = (1.0_f - t) * (1.0_f - u) * f1 + &
                             (        t) * (1.0_f - u) * f2 + &
                             (1.0_f - t) * (        u) * f4 + &
                             (        t) * (        u) * f3 

end function bilinear_interpolation


pure function cheb_poly(n, x) result(T)

    !
    ! Returns the value of the first n (where, n > 2) Chebyshev 
    ! polynomials. Note that the code is optimized for
    ! n > 2, but no safety check are done to ensure this 
    ! that this condition is met!
    ! 
    ! Uses requrrence equation:
    !     T_n(x) = 2xT_(n-1)(x) - T_(n-2)(x)
    !     T_(n=1)(x) = 1.0
    !     T_(n=2)(x) = x
    !

    implicit none

    real(kind=f), dimension(n) :: T
    integer, intent(in)        :: n
    real(kind=f), intent(in)   :: x

    ! local 
    integer      :: k
    real(kind=f) :: x2
    
    x2 = 2.0_f * x

    T(1) = 1.0_f
    T(2) = x

    do k = 3, n
        T(k) = x2 * T(k-1) - T(k-2)
    end do

    return
end function cheb_poly


subroutine lookup(lut_refractive_index_re, & 
                  lut_refractive_index_im, &
                  lut_c_ext,               &
                  lut_c_sca,               &
                  lut_c_g,                 &
                  refractive_index_re,     & 
                  refractive_index_im,     &
                  x,                       &
                  ext,                     &
                  sca,                     &
                  g,                       &
                  n_lut_re,                &
                  n_lut_im,                &
                  n_lut_cheb,              &
                  n_lev,                   &
                  n_clm)

    !
    ! Returns 
    !

    implicit none

    integer, intent(in) :: n_lut_re
    integer, intent(in) :: n_lut_im
    integer, intent(in) :: n_lut_cheb

    integer, intent(in) :: n_lev
    integer, intent(in) :: n_clm

    real(kind=f), dimension(n_lut_re), intent(in) :: lut_refractive_index_re
    real(kind=f), dimension(n_lut_im), intent(in) :: lut_refractive_index_im

    real(kind=f), dimension(n_lut_cheb, n_lut_re, n_lut_im), intent(in) :: lut_c_ext
    real(kind=f), dimension(n_lut_cheb, n_lut_re, n_lut_im), intent(in) :: lut_c_sca
    real(kind=f), dimension(n_lut_cheb, n_lut_re, n_lut_im), intent(in) :: lut_c_g

    real(kind=f), dimension(n_clm, n_lev), intent(in) :: refractive_index_re
    real(kind=f), dimension(n_clm, n_lev), intent(in) :: refractive_index_im 
    real(kind=f), dimension(n_clm, n_lev), intent(in) :: x

    real(kind=f), dimension(n_clm, n_lev), intent(out) :: ext
    real(kind=f), dimension(n_clm, n_lev), intent(out) :: sca
    real(kind=f), dimension(n_clm, n_lev), intent(out) :: g


    ! local
    integer :: c, k, l
    integer :: i_re
    integer :: i_im
    real(kind=f), dimension(n_lut_cheb) :: T
    real(kind=f), dimension(n_lut_cheb) :: c_ext
    real(kind=f), dimension(n_lut_cheb) :: c_sca
    real(kind=f), dimension(n_lut_cheb) :: c_g


loop_levels:  do l = 1, n_lev
loop_columns: do c = 1, n_clm

    ! initialize 
    T     = 0.0_f
    c_ext = 0.0_f
    c_sca = 0.0_f
    c_g   = 0.0_f

    ! locate the indexes in refractive index plane
    i_re = locate(refractive_index_re(c,l), lut_refractive_index_re, n_lut_re)
    i_im = locate(refractive_index_im(c,l), lut_refractive_index_im, n_lut_im)

    ! interpolate in the refractive index plane
    do k = 1, n_lut_cheb
        c_ext(k) = bilinear_interpolation(lut_refractive_index_re(i_re),   &
                                          lut_refractive_index_re(i_re+1), &
                                          lut_refractive_index_im(i_im),   &
                                          lut_refractive_index_im(i_im+1), &
                                          lut_c_ext(k,i_re  ,i_im  ),      &
                                          lut_c_ext(k,i_re+1,i_im  ),      &
                                          lut_c_ext(k,i_re+1,i_im+1),      &
                                          lut_c_ext(k,i_re  ,i_im+1),      &
                                          refractive_index_re(c,l),        &
                                          refractive_index_im(c,l))

        c_sca(k) = bilinear_interpolation(lut_refractive_index_re(i_re),   &
                                          lut_refractive_index_re(i_re+1), &
                                          lut_refractive_index_im(i_im),   &
                                          lut_refractive_index_im(i_im+1), &
                                          lut_c_sca(k,i_re  ,i_im  ),      &
                                          lut_c_sca(k,i_re+1,i_im  ),      &
                                          lut_c_sca(k,i_re+1,i_im+1),      &
                                          lut_c_sca(k,i_re  ,i_im+1),      &
                                          refractive_index_re(c,l),        &
                                          refractive_index_im(c,l))

        c_g(k)   = bilinear_interpolation(lut_refractive_index_re(i_re),   &
                                          lut_refractive_index_re(i_re+1), &
                                          lut_refractive_index_im(i_im),   &
                                          lut_refractive_index_im(i_im+1), &
                                          lut_c_g(k,i_re  ,i_im  ),        &
                                          lut_c_g(k,i_re+1,i_im  ),        &
                                          lut_c_g(k,i_re+1,i_im+1),        &
                                          lut_c_g(k,i_re  ,i_im+1),        &
                                          refractive_index_re(c,l),        &
                                          refractive_index_im(c,l))
    end do


    ! compute Chebyshev polynomials
    T = cheb_poly(n_lut_cheb, x(c,l))

    ext(c,l) = sum(c_ext * T)
    sca(c,l) = sum(c_sca * T)
    g(c,l)   = sum(c_g   * T)

end do loop_columns
end do loop_levels

end subroutine lookup

end module MAML_OpticsMod
