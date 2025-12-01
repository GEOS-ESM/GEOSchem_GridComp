from f90nml import Namelist

from ndsl import StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM, Z_INTERFACE_DIM
from ndsl.dsl.typing import Float, Int
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from ndsl.utils import safe_assign_array
from pyChem.PChem.config import PChemConfiguration
from pyChem.PChem.pchem import PChem
import ndsl.xumpy as xp


class TranslatePChem(TranslateFortranData2Py):
    def __init__(
        self,
        grid,
        namelist: Namelist,
        stencil_factory: StencilFactory,
    ):
        super().__init__(grid, stencil_factory)
        self.stencil_factory = stencil_factory
        self.quantity_factory = grid.quantity_factory

        # FloatField Inputs
        self.in_vars["data_vars"] = {
            "LATS": {},
            "PCHEM_LATS": {},
            "PCHEM_LEVS": {},
            "PLE": {},
            "XX_CH4_in": {},
            "XX_N2O_in": {},
            "XX_CFC11_in": {},
            "XX_CFC12_in": {},
            "XX_HCFC22_in": {},
            "XX_OX_in": {},
            "XX_H2O_in": {},
            "mncv": {},
            "TROPP": {},
            "OX": {},
            "O3VMR": {},
            "ZTH": {},
            "AOA": {},
        }

        # Float/Int Inputs
        self.in_vars["parameters"] = [
            "NN_CH4",
            "NN_N2O",
            "NN_CFC11",
            "NN_CFC12",
            "NN_HCFC22",
            "NN_OX",
            "NN_H2O",
            "clim_years",
            "delp",
            "fac",
            "pcrit",
            "tau",
            "dt",
            "USE_H2O_ProdLoss",
            "O3_pointer",
            "O3PPMV_pointer",
            "TO3_pointer",
            "TTO3_pointer",
        ]

        # FloatField Outputs
        self.out_vars = {
            "XX_CH4": self.grid.compute_dict(),
            "XX_N2O": self.grid.compute_dict(),
            "XX_CFC11": self.grid.compute_dict(),
            "XX_CFC12": self.grid.compute_dict(),
            "XX_HCFC22": self.grid.compute_dict(),
            "XX_OX": self.grid.compute_dict(),
            "XX_H2O": self.grid.compute_dict(),
            "O3": self.grid.compute_dict(),
            "O3PPMV": self.grid.compute_dict(),
            "AOA": self.grid.compute_dict(),
        }

    def compute(self, inputs):
        self.pchem_config = PChemConfiguration(
            Int(inputs["clim_years"]),
            Float(inputs["tau"]),
            Int(inputs["USE_H2O_ProdLoss"]),
        )

        pchem = PChem(
            self.stencil_factory,
            self.grid.quantity_factory,
            self.pchem_config,
        )

        # Float/Int Inputs
        NN_CH4 = Int(inputs["NN_CH4"])
        NN_N2O = Int(inputs["NN_N2O"])
        NN_CFC11 = Int(inputs["NN_CFC11"])
        NN_CFC12 = Int(inputs["NN_CFC12"])
        NN_HCFC22 = Int(inputs["NN_HCFC22"])
        NN_OX = Int(inputs["NN_OX"])
        NN_H2O = Int(inputs["NN_H2O"])
        delp = Float(inputs["delp"])
        fac = Float(inputs["fac"])
        pcrit = Float(inputs["pcrit"])
        dt = Float(inputs["dt"])
        O3_pointer = Int(inputs["O3_pointer"])
        O3PPMV_pointer = Int(inputs["O3PPMV_pointer"])
        TO3_pointer = Int(inputs["TO3_pointer"])
        TTO3_pointer = Int(inputs["TTO3_pointer"])

        # Field inputs
        ple = self.quantity_factory.zeros(
            dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a"
        )
        safe_assign_array(ple.view[:, :, :], inputs["PLE"])

        XX_CH4_in = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(XX_CH4_in.view[:, :, :], inputs["XX_CH4_in"])

        XX_N2O_in = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(XX_N2O_in.view[:, :, :], inputs["XX_N2O_in"])

        XX_CFC11_in = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(XX_CFC11_in.view[:, :, :], inputs["XX_CFC11_in"])

        XX_CFC12_in = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(XX_CFC12_in.view[:, :, :], inputs["XX_CFC12_in"])

        XX_HCFC22_in = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(XX_HCFC22_in.view[:, :, :], inputs["XX_HCFC22_in"])

        XX_OX_in = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(XX_OX_in.view[:, :, :], inputs["XX_OX_in"])

        XX_H2O_in = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(XX_H2O_in.view[:, :, :], inputs["XX_H2O_in"])

        PCHEM_LEVS = self.quantity_factory.zeros([Z_DIM], units="n/a")
        safe_assign_array(PCHEM_LEVS.view[:], inputs["PCHEM_LEVS"])

        PCHEM_LATS = xp.zeros(shape=[91], backend=self.stencil_factory.backend)
        safe_assign_array(PCHEM_LATS[:], inputs["PCHEM_LATS"])

        LATS = self.quantity_factory.zeros([X_DIM, Y_DIM], units="n/a")
        safe_assign_array(LATS.view[:, :], inputs["LATS"])

        mncv = xp.zeros(shape=(91, 72, 7, 2), backend=self.stencil_factory.backend)
        safe_assign_array(mncv[:, :, :, :], inputs["mncv"])

        tropp = self.quantity_factory.zeros([X_DIM, Y_DIM], units="n/a")
        safe_assign_array(tropp.view[:, :], inputs["TROPP"])

        O3VMR = self.quantity_factory.zeros([X_DIM, Y_DIM], units="n/a")
        safe_assign_array(O3VMR.view[:, :], inputs["O3VMR"])

        OX = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(OX.view[:, :, :], inputs["OX"])

        ZTH = self.quantity_factory.zeros([X_DIM, Y_DIM], units="n/a")
        safe_assign_array(ZTH.view[:, :], inputs["ZTH"])

        AOA_in = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(AOA_in.view[:, :, :], inputs["AOA"])

        CH4 = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], units="n/a")

        N2O = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], units="n/a")

        CFC11 = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], units="n/a")

        CFC12 = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], units="n/a")

        HCFC22 = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], units="n/a")

        OX = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], units="n/a")

        H2O = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], units="n/a")
        O3 = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], units="n/a")

        O3PPMV = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], units="n/a")
        AOA = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], units="n/a")

        pchem(
            # In
            NN_CH4=NN_CH4,
            NN_N2O=NN_N2O,
            NN_CFC11=NN_CFC11,
            NN_CFC12=NN_CFC12,
            NN_HCFC22=NN_HCFC22,
            NN_OX=NN_OX,
            NN_H2O=NN_H2O,
            ple=ple,
            delp=delp,
            fac=fac,
            mncv=mncv,
            pcrit=pcrit,
            dt=dt,
            lats=LATS,
            pchem_lats=PCHEM_LATS,
            pchem_levs=PCHEM_LEVS,
            tropp=tropp,
            XX_CH4_in=XX_CH4_in,
            XX_N2O_in=XX_N2O_in,
            XX_CFC11_in=XX_CFC11_in,
            XX_CFC12_in=XX_CFC12_in,
            XX_HCFC22_in=XX_HCFC22_in,
            XX_OX_in=XX_OX_in,
            XX_H2O_in=XX_H2O_in,
            AOA_in=AOA_in,
            O3_pointer=O3_pointer,
            O3PPMV_pointer=O3PPMV_pointer,
            TO3_pointer=TO3_pointer,
            TTO3_pointer=TTO3_pointer,
            O3VMR=O3VMR,
            OX=OX,
            ZTH=ZTH,
            # Out
            CH4=CH4,
            N2O=N2O,
            CFC11=CFC11,
            CFC12=CFC12,
            HCFC22=HCFC22,
            H2O=H2O,
            O3=O3,
            O3PPMV=O3PPMV,
            AOA=AOA,
        )

        return {
            "XX_CH4": CH4.view[:],
            "XX_N2O": N2O.view[:],
            "XX_CFC11": CFC11.view[:],
            "XX_CFC12": CFC12.view[:],
            "XX_HCFC22": HCFC22.view[:],
            "XX_OX": OX.view[:],
            "XX_H2O": H2O.view[:],
            "O3": O3.view[:],
            "O3PPMV": O3PPMV.view[:],
            "AOA": AOA.view[:],
        }
