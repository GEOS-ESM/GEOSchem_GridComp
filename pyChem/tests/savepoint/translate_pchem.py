from ndsl import Namelist, Quantity, QuantityFactory, StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM, Z_INTERFACE_DIM
from ndsl.dsl.typing import Float, FloatField, Int
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from ndsl.utils import safe_assign_array
from pyChem.PChem.pchem import PChem
from pyChem.PChem.config import PChemConfiguration
import numpy as np


class TranslateUpdate(TranslateFortranData2Py):
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
            "XX_in": {},
            "mncv": {},
            "TROPP": {},
        }

        # Float/Int Inputs
        self.in_vars["parameters"] = [
            "NN_CH4",
            "clim_years",
            "delp",
            "fac",
            "pcrit",
            "tau",
            "dt",
        ]

        # FloatField Outputs
        self.out_vars = {
            "XX_CH4": self.grid.compute_dict(),
        }

    def make_ijk_field(self, data, dtype=FloatField) -> Quantity:
        qty = self.quantity_factory.empty([X_DIM, Y_DIM, Z_DIM], "n/a", dtype=dtype)
        qty.view[:, :, :] = qty.np.asarray(data[:, :, :])
        return qty

    def compute(self, inputs):
        self.pchem_config = PChemConfiguration(Int(inputs["clim_years"]))

        pchem = PChem(
            self.stencil_factory,
            self.grid.quantity_factory,
            self.pchem_config,
        )

        # Float/Int Inputs
        NN_CH4 = Int(inputs["NN_CH4"])
        delp = Float(inputs["delp"])
        fac = Float(inputs["fac"])
        pcrit = Float(inputs["pcrit"])
        dt = Float(inputs["dt"])
        tau = Float(inputs["tau"])

        # Field inputs
        ple = QuantityFactory.zeros(
            self.quantity_factory, dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a"
        )
        safe_assign_array(ple.view[:, :, :], inputs["PLE"])
        xx_in = QuantityFactory.zeros(
            self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a"
        )
        safe_assign_array(xx_in.view[:, :, :], inputs["XX_in"])

        PCHEM_LEVS = QuantityFactory.zeros(
            self.quantity_factory, dims=[Z_DIM], units="n/a"
        )
        safe_assign_array(PCHEM_LEVS.view[:], inputs["PCHEM_LEVS"])

        PCHEM_LATS = np.zeros(shape=[91])
        safe_assign_array(PCHEM_LATS[:], inputs["PCHEM_LATS"])

        LATS = QuantityFactory.zeros(
            self.quantity_factory, dims=[X_DIM, Y_DIM], units="n/a"
        )
        safe_assign_array(LATS.view[:, :], inputs["LATS"])

        mncv = np.zeros(shape=(91, 72, 7, 2))
        safe_assign_array(mncv[:, :, :, :], inputs["mncv"])

        tropp = QuantityFactory.zeros(
            self.quantity_factory, dims=[X_DIM, Y_DIM], units="n/a"
        )
        safe_assign_array(tropp.view[:, :], inputs["TROPP"])

        # FloatFields
        xx_CH4 = QuantityFactory.zeros(
            self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a"
        )

        PROD = QuantityFactory.zeros(
            self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a"
        )

        PROD_INT = QuantityFactory.zeros(
            self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a"
        )
        PL = QuantityFactory.zeros(
            self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a"
        )

        pchem(
            # Field inputs
            NN=NN_CH4,
            ple=ple,
            XX_in=xx_in,
            delp=delp,
            fac=fac,
            mncv=mncv,
            pcrit=pcrit,
            dt=dt,
            tau=tau,
            prod=PROD,
            prod_int=PROD_INT,
            pl=PL,
            lats=LATS,
            pchem_lats=PCHEM_LATS,
            pchem_levs=PCHEM_LEVS,
            tropp=tropp,
            # Outputs
            XX=xx_CH4,
        )

        return {
            "XX_CH4": xx_CH4.view[:],
        }
