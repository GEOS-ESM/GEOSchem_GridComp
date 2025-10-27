import numpy as np
from mpi4py import MPI
from typing import TYPE_CHECKING
from ndsl.dsl.typing import Float, Int

from tcn.py_ftn_interface.templates.data_conversion import FortranPythonConversion

if TYPE_CHECKING:
    import cffi


class PYCHEM:
    def __init__(self):
        # Default converter on Host memory. Use `cupy`
        # insteady of Numpy to get a Host<>Device converter
        self._f2py = FortranPythonConversion(target_numpy_module=np)

    def init(
        self,
        # input
        clim_years: int,
        tau: float,
        USE_H2O_ProdLoss: int,
        # input-output
        # output
    ):
        python_state_data = {}
        # Note: Will need to add in dimensions for call to .fortran_to_python

        print("My code for pychem_init goes here.")

        # Note: Will need to add in potential offset for call to .python_to_fortran

    def pchem(
        self,
        # input
        NN_CH4: int,
        NN_N20: int,
        NN_CFC11: int,
        NN_CFC12: int,
        NN_HCFC22: int,
        NN_OX: int,
        NN_H2O: int,
        ple: "cffi.FFI.CData",
        ple_dim_sizes: "cffi.FFI.CData",
        ple_rank: int,
        delp: float,
        fac: float,
        mncv: "cffi.FFI.CData",
        mncv_dim_sizes: "cffi.FFI.CData",
        mncv_rank: int,
        pcrit: float,
        lats: "cffi.FFI.CData",
        lats_dim_sizes: "cffi.FFI.CData",
        lats_rank: int,
        pchem_levs: "cffi.FFI.CData",
        pchem_levs_dim_sizes: "cffi.FFI.CData",
        pchem_levs_rank: int,
        pchem_lats: "cffi.FFI.CData",
        pchem_lats_dim_sizes: "cffi.FFI.CData",
        pchem_lats_rank: int,
        dt: float,
        tropp: "cffi.FFI.CData",
        tropp_dim_sizes: "cffi.FFI.CData",
        tropp_rank: int,
        XX_CH4_in: "cffi.FFI.CData",
        XX_CH4_in_dim_sizes: "cffi.FFI.CData",
        XX_CH4_in_rank: int,
        XX_N2O_in: "cffi.FFI.CData",
        XX_N2O_in_dim_sizes: "cffi.FFI.CData",
        XX_N2O_in_rank: int,
        XX_CFC11_in: "cffi.FFI.CData",
        XX_CFC11_in_dim_sizes: "cffi.FFI.CData",
        XX_CFC11_in_rank: int,
        XX_CFC12_in: "cffi.FFI.CData",
        XX_CFC12_in_dim_sizes: "cffi.FFI.CData",
        XX_CFC12_in_rank: int,
        XX_HCFC22_in: "cffi.FFI.CData",
        XX_HCFC22_in_dim_sizes: "cffi.FFI.CData",
        XX_HCFC22_in_rank: int,
        XX_OX_in: "cffi.FFI.CData",
        XX_OX_in_dim_sizes: "cffi.FFI.CData",
        XX_OX_in_rank: int,
        XX_H2O_in: "cffi.FFI.CData",
        XX_H2O_in_dim_sizes: "cffi.FFI.CData",
        XX_H2O_in_rank: int,
        AOA_in: "cffi.FFI.CData",
        AOA_in_dim_sizes: "cffi.FFI.CData",
        AOA_in_rank: int,
        O3_pointer: int,
        O3PPMV_pointer: int,
        TO3_pointer: int,
        TTO3_pointer: int,
        O3VMR: "cffi.FFI.CData",
        O3VMR_dim_sizes: "cffi.FFI.CData",
        O3VMR_rank: int,
        OX: "cffi.FFI.CData",
        OX_dim_sizes: "cffi.FFI.CData",
        OX_rank: int,
        ZTH: "cffi.FFI.CData",
        ZTH_dim_sizes: "cffi.FFI.CData",
        ZTH_rank: int,
        # input-output
        # output
        CH4: "cffi.FFI.CData",
        CH4_dim_sizes: "cffi.FFI.CData",
        CH4_rank: int,
        N2O: "cffi.FFI.CData",
        N2O_dim_sizes: "cffi.FFI.CData",
        N2O_rank: int,
        CFC11: "cffi.FFI.CData",
        CFC11_dim_sizes: "cffi.FFI.CData",
        CFC11_rank: int,
        CFC12: "cffi.FFI.CData",
        CFC12_dim_sizes: "cffi.FFI.CData",
        CFC12_rank: int,
        HCFC22: "cffi.FFI.CData",
        HCFC22_dim_sizes: "cffi.FFI.CData",
        HCFC22_rank: int,
        H2O: "cffi.FFI.CData",
        H2O_dim_sizes: "cffi.FFI.CData",
        H2O_rank: int,
        O3: "cffi.FFI.CData",
        O3_dim_sizes: "cffi.FFI.CData",
        O3_rank: int,
        O3PPMV: "cffi.FFI.CData",
        O3PPMV_dim_sizes: "cffi.FFI.CData",
        O3PPMV_rank: int,
        AOA: "cffi.FFI.CData",
        AOA_dim_sizes: "cffi.FFI.CData",
        AOA_rank: int,
    ):
        python_state_data = {}
        # Note: Will need to add in dimensions for call to .fortran_to_python

        python_state_data["ple"] = self._f2py.fortran_to_python(
            ple, dim=ple_dim_sizes, rank=ple_rank
        )

        python_state_data["mncv"] = self._f2py.fortran_to_python(
            mncv, dim=mncv_dim_sizes, rank=mncv_rank
        )

        python_state_data["lats"] = self._f2py.fortran_to_python(
            lats, dim=lats_dim_sizes, rank=lats_rank
        )
        python_state_data["pchem_levs"] = self._f2py.fortran_to_python(
            pchem_levs, dim=pchem_levs_dim_sizes, rank=pchem_levs_rank
        )
        python_state_data["pchem_lats"] = self._f2py.fortran_to_python(
            pchem_lats, dim=pchem_lats_dim_sizes, rank=pchem_lats_rank
        )

        python_state_data["tropp"] = self._f2py.fortran_to_python(
            tropp, dim=tropp_dim_sizes, rank=tropp_rank
        )
        python_state_data["XX_CH4_in"] = self._f2py.fortran_to_python(
            XX_CH4_in, dim=XX_CH4_in_dim_sizes, rank=XX_CH4_in_rank
        )
        python_state_data["XX_N2O_in"] = self._f2py.fortran_to_python(
            XX_N2O_in, dim=XX_N2O_in_dim_sizes, rank=XX_N2O_in_rank
        )
        python_state_data["XX_CFC11_in"] = self._f2py.fortran_to_python(
            XX_CFC11_in, dim=XX_CFC11_in_dim_sizes, rank=XX_CFC11_in_rank
        )
        python_state_data["XX_CFC12_in"] = self._f2py.fortran_to_python(
            XX_CFC12_in, dim=XX_CFC12_in_dim_sizes, rank=XX_CFC12_in_rank
        )
        python_state_data["XX_HCFC22_in"] = self._f2py.fortran_to_python(
            XX_HCFC22_in, dim=XX_HCFC22_in_dim_sizes, rank=XX_HCFC22_in_rank
        )
        python_state_data["XX_OX_in"] = self._f2py.fortran_to_python(
            XX_OX_in, dim=XX_OX_in_dim_sizes, rank=XX_OX_in_rank
        )
        python_state_data["XX_H2O_in"] = self._f2py.fortran_to_python(
            XX_H2O_in, dim=XX_H2O_in_dim_sizes, rank=XX_H2O_in_rank
        )
        python_state_data["AOA_in"] = self._f2py.fortran_to_python(
            AOA_in, dim=AOA_in_dim_sizes, rank=AOA_in_rank
        )

        python_state_data["O3VMR"] = self._f2py.fortran_to_python(
            O3VMR, dim=O3VMR_dim_sizes, rank=O3VMR_rank
        )
        python_state_data["OX"] = self._f2py.fortran_to_python(
            OX, dim=OX_dim_sizes, rank=OX_rank
        )
        python_state_data["ZTH"] = self._f2py.fortran_to_python(
            ZTH, dim=ZTH_dim_sizes, rank=ZTH_rank
        )

        python_state_data["CH4"] = self._f2py.fortran_to_python(
            CH4, dim=CH4_dim_sizes, rank=CH4_rank
        )
        python_state_data["N2O"] = self._f2py.fortran_to_python(
            N2O, dim=N2O_dim_sizes, rank=N2O_rank
        )
        python_state_data["CFC11"] = self._f2py.fortran_to_python(
            CFC11, dim=CFC11_dim_sizes, rank=CFC11_rank
        )
        python_state_data["CFC12"] = self._f2py.fortran_to_python(
            CFC12, dim=CFC12_dim_sizes, rank=CFC12_rank
        )
        python_state_data["HCFC22"] = self._f2py.fortran_to_python(
            HCFC22, dim=HCFC22_dim_sizes, rank=HCFC22_rank
        )
        python_state_data["H2O"] = self._f2py.fortran_to_python(
            H2O, dim=H2O_dim_sizes, rank=H2O_rank
        )
        python_state_data["O3"] = self._f2py.fortran_to_python(
            O3, dim=O3_dim_sizes, rank=O3_rank
        )
        python_state_data["O3PPMV"] = self._f2py.fortran_to_python(
            O3PPMV, dim=O3PPMV_dim_sizes, rank=O3PPMV_rank
        )
        python_state_data["AOA"] = self._f2py.fortran_to_python(
            AOA, dim=AOA_dim_sizes, rank=AOA_rank
        )

        print("My code for pychem_pchem goes here.")

        # Note: Will need to add in potential offset for call to .python_to_fortran

        self._f2py.python_to_fortran(python_state_data["CH4"], CH4)
        self._f2py.python_to_fortran(python_state_data["N2O"], N2O)
        self._f2py.python_to_fortran(python_state_data["CFC11"], CFC11)
        self._f2py.python_to_fortran(python_state_data["CFC12"], CFC12)
        self._f2py.python_to_fortran(python_state_data["HCFC22"], HCFC22)
        self._f2py.python_to_fortran(python_state_data["H2O"], H2O)
        self._f2py.python_to_fortran(python_state_data["O3"], O3)
        self._f2py.python_to_fortran(python_state_data["O3PPMV"], O3PPMV)
        self._f2py.python_to_fortran(python_state_data["AOA"], AOA)


pychem = PYCHEM()
