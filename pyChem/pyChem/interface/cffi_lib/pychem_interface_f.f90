module pychem_interface_mod

   use iso_c_binding, only: c_int, c_float, c_double, c_bool

   implicit none

   !> Single precision real numbers, 6 digits, range 10⁻³⁷ to 10³⁷-1; 32 bits
   integer, parameter :: sp = selected_real_kind(6, 37)
   !> Double precision real numbers, 15 digits, range 10⁻³⁰⁷ to 10³⁰⁷-1; 64 bits
   integer, parameter :: dp = selected_real_kind(15, 307)

   private
   public :: pychem_f_init
   public :: pychem_f_pchem

   interface

      subroutine pychem_f_init( &
         !inputs
         clim_years, &
         tau, &
         USE_H2O_ProdLoss &
         !inouts
         !outputs
         ) bind(c, name='pychem_c_init')
         import c_int, c_float, c_double
         implicit none

         integer(kind=c_int), value, intent(in) :: clim_years

         real(kind=c_float), value, intent(in) :: tau

         integer(kind=c_int), value, intent(in) :: USE_H2O_ProdLoss

      end subroutine pychem_f_init

      subroutine pychem_f_pchem( &
         !inputs
         NN_CH4, &
         NN_N20, &
         NN_CFC11, &
         NN_CFC12, &
         NN_HCFC22, &
         NN_OX, &
         NN_H2O, &
         ple, ple_dim_sizes, ple_rank, &
         delp, &
         fac, &
         mncv, mncv_dim_sizes, mncv_rank, &
         pcrit, &
         lats, lats_dim_sizes, lats_rank, &
         pchem_levs, pchem_levs_dim_sizes, pchem_levs_rank, &
         pchem_lats, pchem_lats_dim_sizes, pchem_lats_rank, &
         dt, &
         tropp, tropp_dim_sizes, tropp_rank, &
         XX_CH4_in, XX_CH4_in_dim_sizes, XX_CH4_in_rank, &
         XX_N2O_in, XX_N2O_in_dim_sizes, XX_N2O_in_rank, &
         XX_CFC11_in, XX_CFC11_in_dim_sizes, XX_CFC11_in_rank, &
         XX_CFC12_in, XX_CFC12_in_dim_sizes, XX_CFC12_in_rank, &
         XX_HCFC22_in, XX_HCFC22_in_dim_sizes, XX_HCFC22_in_rank, &
         XX_OX_in, XX_OX_in_dim_sizes, XX_OX_in_rank, &
         XX_H2O_in, XX_H2O_in_dim_sizes, XX_H2O_in_rank, &
         AOA_in, AOA_in_dim_sizes, AOA_in_rank, &
         O3_pointer, &
         O3PPMV_pointer, &
         TO3_pointer, &
         TTO3_pointer, &
         O3VMR, O3VMR_dim_sizes, O3VMR_rank, &
         OX, OX_dim_sizes, OX_rank, &
         ZTH, ZTH_dim_sizes, ZTH_rank, &
         !inouts
         !outputs
         CH4, CH4_dim_sizes, CH4_rank, &
         N2O, N2O_dim_sizes, N2O_rank, &
         CFC11, CFC11_dim_sizes, CFC11_rank, &
         CFC12, CFC12_dim_sizes, CFC12_rank, &
         HCFC22, HCFC22_dim_sizes, HCFC22_rank, &
         H2O, H2O_dim_sizes, H2O_rank, &
         O3, O3_dim_sizes, O3_rank, &
         O3PPMV, O3PPMV_dim_sizes, O3PPMV_rank, &
         AOA, AOA_dim_sizes, AOA_rank &
         ) bind(c, name='pychem_c_pchem')
         import c_int, c_float, c_double
         implicit none

         integer(kind=c_int), value, intent(in) :: NN_CH4

         integer(kind=c_int), value, intent(in) :: NN_N20

         integer(kind=c_int), value, intent(in) :: NN_CFC11

         integer(kind=c_int), value, intent(in) :: NN_CFC12

         integer(kind=c_int), value, intent(in) :: NN_HCFC22

         integer(kind=c_int), value, intent(in) :: NN_OX

         integer(kind=c_int), value, intent(in) :: NN_H2O

         real(kind=c_float), dimension(*), intent(in) :: ple
         integer(kind=c_int), dimension(*), intent(in) :: ple_dim_sizes
         integer(kind=c_int), value, intent(in) :: ple_rank

         real(kind=c_float), value, intent(in) :: delp

         real(kind=c_float), value, intent(in) :: fac

         real(kind=c_float), dimension(*), intent(in) :: mncv
         integer(kind=c_int), dimension(*), intent(in) :: mncv_dim_sizes
         integer(kind=c_int), value, intent(in) :: mncv_rank

         real(kind=c_float), value, intent(in) :: pcrit

         real(kind=c_float), dimension(*), intent(in) :: lats
         integer(kind=c_int), dimension(*), intent(in) :: lats_dim_sizes
         integer(kind=c_int), value, intent(in) :: lats_rank

         real(kind=c_float), dimension(*), intent(in) :: pchem_levs
         integer(kind=c_int), dimension(*), intent(in) :: pchem_levs_dim_sizes
         integer(kind=c_int), value, intent(in) :: pchem_levs_rank

         real(kind=c_float), dimension(*), intent(in) :: pchem_lats
         integer(kind=c_int), dimension(*), intent(in) :: pchem_lats_dim_sizes
         integer(kind=c_int), value, intent(in) :: pchem_lats_rank

         real(kind=c_float), value, intent(in) :: dt

         real(kind=c_float), dimension(*), intent(in) :: tropp
         integer(kind=c_int), dimension(*), intent(in) :: tropp_dim_sizes
         integer(kind=c_int), value, intent(in) :: tropp_rank

         real(kind=c_float), dimension(*), intent(in) :: XX_CH4_in
         integer(kind=c_int), dimension(*), intent(in) :: XX_CH4_in_dim_sizes
         integer(kind=c_int), value, intent(in) :: XX_CH4_in_rank

         real(kind=c_float), dimension(*), intent(in) :: XX_N2O_in
         integer(kind=c_int), dimension(*), intent(in) :: XX_N2O_in_dim_sizes
         integer(kind=c_int), value, intent(in) :: XX_N2O_in_rank

         real(kind=c_float), dimension(*), intent(in) :: XX_CFC11_in
         integer(kind=c_int), dimension(*), intent(in) :: XX_CFC11_in_dim_sizes
         integer(kind=c_int), value, intent(in) :: XX_CFC11_in_rank

         real(kind=c_float), dimension(*), intent(in) :: XX_CFC12_in
         integer(kind=c_int), dimension(*), intent(in) :: XX_CFC12_in_dim_sizes
         integer(kind=c_int), value, intent(in) :: XX_CFC12_in_rank

         real(kind=c_float), dimension(*), intent(in) :: XX_HCFC22_in
         integer(kind=c_int), dimension(*), intent(in) :: XX_HCFC22_in_dim_sizes
         integer(kind=c_int), value, intent(in) :: XX_HCFC22_in_rank

         real(kind=c_float), dimension(*), intent(in) :: XX_OX_in
         integer(kind=c_int), dimension(*), intent(in) :: XX_OX_in_dim_sizes
         integer(kind=c_int), value, intent(in) :: XX_OX_in_rank

         real(kind=c_float), dimension(*), intent(in) :: XX_H2O_in
         integer(kind=c_int), dimension(*), intent(in) :: XX_H2O_in_dim_sizes
         integer(kind=c_int), value, intent(in) :: XX_H2O_in_rank

         real(kind=c_float), dimension(*), intent(in) :: AOA_in
         integer(kind=c_int), dimension(*), intent(in) :: AOA_in_dim_sizes
         integer(kind=c_int), value, intent(in) :: AOA_in_rank

         integer(kind=c_int), value, intent(in) :: O3_pointer

         integer(kind=c_int), value, intent(in) :: O3PPMV_pointer

         integer(kind=c_int), value, intent(in) :: TO3_pointer

         integer(kind=c_int), value, intent(in) :: TTO3_pointer

         real(kind=c_float), dimension(*), intent(in) :: O3VMR
         integer(kind=c_int), dimension(*), intent(in) :: O3VMR_dim_sizes
         integer(kind=c_int), value, intent(in) :: O3VMR_rank

         real(kind=c_float), dimension(*), intent(in) :: OX
         integer(kind=c_int), dimension(*), intent(in) :: OX_dim_sizes
         integer(kind=c_int), value, intent(in) :: OX_rank

         real(kind=c_float), dimension(*), intent(in) :: ZTH
         integer(kind=c_int), dimension(*), intent(in) :: ZTH_dim_sizes
         integer(kind=c_int), value, intent(in) :: ZTH_rank

         real(kind=c_float), dimension(*), intent(out) :: CH4
         integer(kind=c_int), dimension(*), intent(in) :: CH4_dim_sizes
         integer(kind=c_int), value, intent(in) :: CH4_rank

         real(kind=c_float), dimension(*), intent(out) :: N2O
         integer(kind=c_int), dimension(*), intent(in) :: N2O_dim_sizes
         integer(kind=c_int), value, intent(in) :: N2O_rank

         real(kind=c_float), dimension(*), intent(out) :: CFC11
         integer(kind=c_int), dimension(*), intent(in) :: CFC11_dim_sizes
         integer(kind=c_int), value, intent(in) :: CFC11_rank

         real(kind=c_float), dimension(*), intent(out) :: CFC12
         integer(kind=c_int), dimension(*), intent(in) :: CFC12_dim_sizes
         integer(kind=c_int), value, intent(in) :: CFC12_rank

         real(kind=c_float), dimension(*), intent(out) :: HCFC22
         integer(kind=c_int), dimension(*), intent(in) :: HCFC22_dim_sizes
         integer(kind=c_int), value, intent(in) :: HCFC22_rank

         real(kind=c_float), dimension(*), intent(out) :: H2O
         integer(kind=c_int), dimension(*), intent(in) :: H2O_dim_sizes
         integer(kind=c_int), value, intent(in) :: H2O_rank

         real(kind=c_float), dimension(*), intent(out) :: O3
         integer(kind=c_int), dimension(*), intent(in) :: O3_dim_sizes
         integer(kind=c_int), value, intent(in) :: O3_rank

         real(kind=c_float), dimension(*), intent(out) :: O3PPMV
         integer(kind=c_int), dimension(*), intent(in) :: O3PPMV_dim_sizes
         integer(kind=c_int), value, intent(in) :: O3PPMV_rank

         real(kind=c_float), dimension(*), intent(out) :: AOA
         integer(kind=c_int), dimension(*), intent(in) :: AOA_dim_sizes
         integer(kind=c_int), value, intent(in) :: AOA_rank

      end subroutine pychem_f_pchem

   end interface

contains

end module pychem_interface_mod
