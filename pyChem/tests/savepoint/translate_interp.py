from ndsl import Namelist, Quantity, QuantityFactory, StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import FloatField
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from ndsl.utils import safe_assign_array

# from pyChem.MAPL.interp_no_extrap import Interp
import numpy as np


class TranslateInterp(TranslateFortranData2Py):
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
            "PL": {},
            "PCHEM_LEVS": {},
            "Prod1": {},
        }

        # FloatField Outputs
        self.out_vars = {
            "PROD_INT": self.grid.compute_dict(),
        }

    def compute(self, inputs):

        interp = Interp(
            self.stencil_factory,
            self.grid.quantity_factory,
        )

        # Field inputs
        PL = QuantityFactory.zeros(
            self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a"
        )
        safe_assign_array(PL.view[:, :, :], inputs["PL"])

        Prod1 = np.zeros(shape=(91, 72))
        safe_assign_array(Prod1[:, :], inputs["Prod1"])

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

        # FloatFields
        PROD = QuantityFactory.zeros(
            self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a"
        )

        PROD_INT = QuantityFactory.zeros(
            self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a"
        )

        # interp(
        #     # Field inputs
        #     PL=PL,
        #     PROD=PROD,
        #     PCHEM_LEVS=PCHEM_LEVS,
        #     # Outputs
        #     PROD_INT=PROD_INT,
        # )

        im, jm, n_levs = 24, 24, 72

        for j in range(jm):
            for k in range(n_levs):
                PROD.field[:, j, k] = interp_no_extrap(
                    OX_list=LATS.field[:, j],
                    IY=Prod1[:, k],
                    IX=PCHEM_LATS[:],
                )
            for i in range(im):
                PROD_INT.field[i, j, :] = interp_no_extrap(
                    OX_list=PL.field[i, j, :],
                    IY=PROD.field[i, j, :],
                    IX=PCHEM_LEVS.field[:],
                )

        print("Done")

        return {
            "PROD_INT": PROD_INT.field[:],
        }
