#include "mpi.h"
#include "pychem_interface.h"
#include <stdio.h>
#include <time.h>

void pychem_c_init(
    // inputs
    int clim_years, float tau, int USE_H2O_ProdLoss
    // inputs-outputs

    // outputs

) {
  pychem_py_init(clim_years, tau, USE_H2O_ProdLoss

  );

  // if (return_code < 0)
  //{
  //     exit(return_code);
  // }
}

void pychem_c_pchem(
    // inputs
    int NN_CH4, int NN_N20, int NN_CFC11, int NN_CFC12, int NN_HCFC22,
    int NN_OX, int NN_H2O, const float *ple, int *ple_dim_sizes, int ple_rank,
    float delp, float fac, const float *mncv, int *mncv_dim_sizes,
    int mncv_rank, float pcrit, const float *lats, int *lats_dim_sizes,
    int lats_rank, const float *pchem_levs, int *pchem_levs_dim_sizes,
    int pchem_levs_rank, const float *pchem_lats, int *pchem_lats_dim_sizes,
    int pchem_lats_rank, float dt, const float *tropp, int *tropp_dim_sizes,
    int tropp_rank, const float *XX_CH4_in, int *XX_CH4_in_dim_sizes,
    int XX_CH4_in_rank, const float *XX_N2O_in, int *XX_N2O_in_dim_sizes,
    int XX_N2O_in_rank, const float *XX_CFC11_in, int *XX_CFC11_in_dim_sizes,
    int XX_CFC11_in_rank, const float *XX_CFC12_in, int *XX_CFC12_in_dim_sizes,
    int XX_CFC12_in_rank, const float *XX_HCFC22_in,
    int *XX_HCFC22_in_dim_sizes, int XX_HCFC22_in_rank, const float *XX_OX_in,
    int *XX_OX_in_dim_sizes, int XX_OX_in_rank, const float *XX_H2O_in,
    int *XX_H2O_in_dim_sizes, int XX_H2O_in_rank, const float *AOA_in,
    int *AOA_in_dim_sizes, int AOA_in_rank, int O3_pointer, int O3PPMV_pointer,
    int TO3_pointer, int TTO3_pointer, const float *O3VMR, int *O3VMR_dim_sizes,
    int O3VMR_rank, const float *OX, int *OX_dim_sizes, int OX_rank,
    const float *ZTH, int *ZTH_dim_sizes, int ZTH_rank,
    // inputs-outputs

    // outputs
    float *CH4, int *CH4_dim_sizes, int CH4_rank, float *N2O,
    int *N2O_dim_sizes, int N2O_rank, float *CFC11, int *CFC11_dim_sizes,
    int CFC11_rank, float *CFC12, int *CFC12_dim_sizes, int CFC12_rank,
    float *HCFC22, int *HCFC22_dim_sizes, int HCFC22_rank, float *H2O,
    int *H2O_dim_sizes, int H2O_rank, float *O3, int *O3_dim_sizes, int O3_rank,
    float *O3PPMV, int *O3PPMV_dim_sizes, int O3PPMV_rank, float *AOA,
    int *AOA_dim_sizes, int AOA_rank) {
  pychem_py_pchem(
      NN_CH4, NN_N20, NN_CFC11, NN_CFC12, NN_HCFC22, NN_OX, NN_H2O, ple,
      ple_dim_sizes, ple_rank, delp, fac, mncv, mncv_dim_sizes, mncv_rank,
      pcrit, lats, lats_dim_sizes, lats_rank, pchem_levs, pchem_levs_dim_sizes,
      pchem_levs_rank, pchem_lats, pchem_lats_dim_sizes, pchem_lats_rank, dt,
      tropp, tropp_dim_sizes, tropp_rank, XX_CH4_in, XX_CH4_in_dim_sizes,
      XX_CH4_in_rank, XX_N2O_in, XX_N2O_in_dim_sizes, XX_N2O_in_rank,
      XX_CFC11_in, XX_CFC11_in_dim_sizes, XX_CFC11_in_rank, XX_CFC12_in,
      XX_CFC12_in_dim_sizes, XX_CFC12_in_rank, XX_HCFC22_in,
      XX_HCFC22_in_dim_sizes, XX_HCFC22_in_rank, XX_OX_in, XX_OX_in_dim_sizes,
      XX_OX_in_rank, XX_H2O_in, XX_H2O_in_dim_sizes, XX_H2O_in_rank, AOA_in,
      AOA_in_dim_sizes, AOA_in_rank, O3_pointer, O3PPMV_pointer, TO3_pointer,
      TTO3_pointer, O3VMR, O3VMR_dim_sizes, O3VMR_rank, OX, OX_dim_sizes,
      OX_rank, ZTH, ZTH_dim_sizes, ZTH_rank, CH4, CH4_dim_sizes, CH4_rank, N2O,
      N2O_dim_sizes, N2O_rank, CFC11, CFC11_dim_sizes, CFC11_rank, CFC12,
      CFC12_dim_sizes, CFC12_rank, HCFC22, HCFC22_dim_sizes, HCFC22_rank, H2O,
      H2O_dim_sizes, H2O_rank, O3, O3_dim_sizes, O3_rank, O3PPMV,
      O3PPMV_dim_sizes, O3PPMV_rank, AOA, AOA_dim_sizes, AOA_rank

  );

  // if (return_code < 0)
  //{
  //     exit(return_code);
  // }
}
