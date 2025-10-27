#pragma once

#include <stdbool.h>
#include <stdlib.h>

typedef union {
  int comm_int;
  void *comm_ptr;
} MPI_Comm_t;

extern void pychem_py_init(
    // inputs
    int clim_years, float tau, int USE_H2O_ProdLoss
    // inputs-outputs

    // outputs

);

extern void pychem_py_pchem(
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
    int *AOA_dim_sizes, int AOA_rank);
