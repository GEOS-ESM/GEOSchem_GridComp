import cffi
from mpi4py import MPI


TMPFILEBASE = "pychem_interface"

ffi = cffi.FFI()

# MPI_Comm can be int or void*
if MPI._sizeof(MPI.Comm) == ffi.sizeof("int"):
    _mpi_comm_t = "int"
else:
    _mpi_comm_t = "void*"


source = f"""
from {TMPFILEBASE} import ffi
from datetime import datetime
from mpi4py import MPI
import traceback

from pychem_hook import pychem



@ffi.def_extern()
def pychem_py_init(
      # inputs
      clim_years:int,
      tau:float,
      USE_H2O_ProdLoss:int,

      # inputs-outputs

      # outputs

    ) -> int:

    # Transform init code


    pychem.init(

      clim_years=clim_years,
      tau=tau,
      USE_H2O_ProdLoss=USE_H2O_ProdLoss,



    )

@ffi.def_extern()
def pychem_py_pchem(
      # inputs
      NN_CH4:int,
      NN_N20:int,
      NN_CFC11:int,
      NN_CFC12:int,
      NN_HCFC22:int,
      NN_OX:int,
      NN_H2O:int,
      ple:'cffi.FFI.CData',
      ple_dim_sizes:'cffi.FFI.CData',
      ple_rank:int,
      delp:float,
      fac:float,
      mncv:'cffi.FFI.CData',
      mncv_dim_sizes:'cffi.FFI.CData',
      mncv_rank:int,
      pcrit:float,
      lats:'cffi.FFI.CData',
      lats_dim_sizes:'cffi.FFI.CData',
      lats_rank:int,
      pchem_levs:'cffi.FFI.CData',
      pchem_levs_dim_sizes:'cffi.FFI.CData',
      pchem_levs_rank:int,
      pchem_lats:'cffi.FFI.CData',
      pchem_lats_dim_sizes:'cffi.FFI.CData',
      pchem_lats_rank:int,
      dt:float,
      tropp:'cffi.FFI.CData',
      tropp_dim_sizes:'cffi.FFI.CData',
      tropp_rank:int,
      XX_CH4_in:'cffi.FFI.CData',
      XX_CH4_in_dim_sizes:'cffi.FFI.CData',
      XX_CH4_in_rank:int,
      XX_N2O_in:'cffi.FFI.CData',
      XX_N2O_in_dim_sizes:'cffi.FFI.CData',
      XX_N2O_in_rank:int,
      XX_CFC11_in:'cffi.FFI.CData',
      XX_CFC11_in_dim_sizes:'cffi.FFI.CData',
      XX_CFC11_in_rank:int,
      XX_CFC12_in:'cffi.FFI.CData',
      XX_CFC12_in_dim_sizes:'cffi.FFI.CData',
      XX_CFC12_in_rank:int,
      XX_HCFC22_in:'cffi.FFI.CData',
      XX_HCFC22_in_dim_sizes:'cffi.FFI.CData',
      XX_HCFC22_in_rank:int,
      XX_OX_in:'cffi.FFI.CData',
      XX_OX_in_dim_sizes:'cffi.FFI.CData',
      XX_OX_in_rank:int,
      XX_H2O_in:'cffi.FFI.CData',
      XX_H2O_in_dim_sizes:'cffi.FFI.CData',
      XX_H2O_in_rank:int,
      AOA_in:'cffi.FFI.CData',
      AOA_in_dim_sizes:'cffi.FFI.CData',
      AOA_in_rank:int,
      O3_pointer:int,
      O3PPMV_pointer:int,
      TO3_pointer:int,
      TTO3_pointer:int,
      O3VMR:'cffi.FFI.CData',
      O3VMR_dim_sizes:'cffi.FFI.CData',
      O3VMR_rank:int,
      OX:'cffi.FFI.CData',
      OX_dim_sizes:'cffi.FFI.CData',
      OX_rank:int,
      ZTH:'cffi.FFI.CData',
      ZTH_dim_sizes:'cffi.FFI.CData',
      ZTH_rank:int,

      # inputs-outputs

      # outputs
      CH4:'cffi.FFI.CData',
      CH4_dim_sizes:'cffi.FFI.CData',
      CH4_rank:int,
      N2O:'cffi.FFI.CData',
      N2O_dim_sizes:'cffi.FFI.CData',
      N2O_rank:int,
      CFC11:'cffi.FFI.CData',
      CFC11_dim_sizes:'cffi.FFI.CData',
      CFC11_rank:int,
      CFC12:'cffi.FFI.CData',
      CFC12_dim_sizes:'cffi.FFI.CData',
      CFC12_rank:int,
      HCFC22:'cffi.FFI.CData',
      HCFC22_dim_sizes:'cffi.FFI.CData',
      HCFC22_rank:int,
      H2O:'cffi.FFI.CData',
      H2O_dim_sizes:'cffi.FFI.CData',
      H2O_rank:int,
      O3:'cffi.FFI.CData',
      O3_dim_sizes:'cffi.FFI.CData',
      O3_rank:int,
      O3PPMV:'cffi.FFI.CData',
      O3PPMV_dim_sizes:'cffi.FFI.CData',
      O3PPMV_rank:int,
      AOA:'cffi.FFI.CData',
      AOA_dim_sizes:'cffi.FFI.CData',
      AOA_rank:int,

    ) -> int:

    # Transform init code


    pychem.pchem(

      NN_CH4=NN_CH4,
      NN_N20=NN_N20,
      NN_CFC11=NN_CFC11,
      NN_CFC12=NN_CFC12,
      NN_HCFC22=NN_HCFC22,
      NN_OX=NN_OX,
      NN_H2O=NN_H2O,
      ple=ple,
      ple_dim_sizes=ple_dim_sizes,
      ple_rank=ple_rank,
      delp=delp,
      fac=fac,
      mncv=mncv,
      mncv_dim_sizes=mncv_dim_sizes,
      mncv_rank=mncv_rank,
      pcrit=pcrit,
      lats=lats,
      lats_dim_sizes=lats_dim_sizes,
      lats_rank=lats_rank,
      pchem_levs=pchem_levs,
      pchem_levs_dim_sizes=pchem_levs_dim_sizes,
      pchem_levs_rank=pchem_levs_rank,
      pchem_lats=pchem_lats,
      pchem_lats_dim_sizes=pchem_lats_dim_sizes,
      pchem_lats_rank=pchem_lats_rank,
      dt=dt,
      tropp=tropp,
      tropp_dim_sizes=tropp_dim_sizes,
      tropp_rank=tropp_rank,
      XX_CH4_in=XX_CH4_in,
      XX_CH4_in_dim_sizes=XX_CH4_in_dim_sizes,
      XX_CH4_in_rank=XX_CH4_in_rank,
      XX_N2O_in=XX_N2O_in,
      XX_N2O_in_dim_sizes=XX_N2O_in_dim_sizes,
      XX_N2O_in_rank=XX_N2O_in_rank,
      XX_CFC11_in=XX_CFC11_in,
      XX_CFC11_in_dim_sizes=XX_CFC11_in_dim_sizes,
      XX_CFC11_in_rank=XX_CFC11_in_rank,
      XX_CFC12_in=XX_CFC12_in,
      XX_CFC12_in_dim_sizes=XX_CFC12_in_dim_sizes,
      XX_CFC12_in_rank=XX_CFC12_in_rank,
      XX_HCFC22_in=XX_HCFC22_in,
      XX_HCFC22_in_dim_sizes=XX_HCFC22_in_dim_sizes,
      XX_HCFC22_in_rank=XX_HCFC22_in_rank,
      XX_OX_in=XX_OX_in,
      XX_OX_in_dim_sizes=XX_OX_in_dim_sizes,
      XX_OX_in_rank=XX_OX_in_rank,
      XX_H2O_in=XX_H2O_in,
      XX_H2O_in_dim_sizes=XX_H2O_in_dim_sizes,
      XX_H2O_in_rank=XX_H2O_in_rank,
      AOA_in=AOA_in,
      AOA_in_dim_sizes=AOA_in_dim_sizes,
      AOA_in_rank=AOA_in_rank,
      O3_pointer=O3_pointer,
      O3PPMV_pointer=O3PPMV_pointer,
      TO3_pointer=TO3_pointer,
      TTO3_pointer=TTO3_pointer,
      O3VMR=O3VMR,
      O3VMR_dim_sizes=O3VMR_dim_sizes,
      O3VMR_rank=O3VMR_rank,
      OX=OX,
      OX_dim_sizes=OX_dim_sizes,
      OX_rank=OX_rank,
      ZTH=ZTH,
      ZTH_dim_sizes=ZTH_dim_sizes,
      ZTH_rank=ZTH_rank,


      CH4=CH4,
      CH4_dim_sizes=CH4_dim_sizes,
      CH4_rank=CH4_rank,
      N2O=N2O,
      N2O_dim_sizes=N2O_dim_sizes,
      N2O_rank=N2O_rank,
      CFC11=CFC11,
      CFC11_dim_sizes=CFC11_dim_sizes,
      CFC11_rank=CFC11_rank,
      CFC12=CFC12,
      CFC12_dim_sizes=CFC12_dim_sizes,
      CFC12_rank=CFC12_rank,
      HCFC22=HCFC22,
      HCFC22_dim_sizes=HCFC22_dim_sizes,
      HCFC22_rank=HCFC22_rank,
      H2O=H2O,
      H2O_dim_sizes=H2O_dim_sizes,
      H2O_rank=H2O_rank,
      O3=O3,
      O3_dim_sizes=O3_dim_sizes,
      O3_rank=O3_rank,
      O3PPMV=O3PPMV,
      O3PPMV_dim_sizes=O3PPMV_dim_sizes,
      O3PPMV_rank=O3PPMV_rank,
      AOA=AOA,
      AOA_dim_sizes=AOA_dim_sizes,
      AOA_rank=AOA_rank,

    )



"""


with open(f"{TMPFILEBASE}.h") as f:
    data = "".join([line for line in f if not line.startswith("#")])
    data = data.replace("CFFI_DLLEXPORT", "")
    ffi.embedding_api(data)

ffi.set_source(TMPFILEBASE, f'#include "{TMPFILEBASE}.h"')

ffi.embedding_init_code(source)
ffi.compile(target="lib" + TMPFILEBASE + ".so", verbose=True)
