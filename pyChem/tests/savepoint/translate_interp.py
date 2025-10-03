from ndsl import Namelist, Quantity, QuantityFactory, StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import FloatField
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from ndsl.utils import safe_assign_array
from pyChem.PChem.interp_no_extrap import Interp
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
            "PL": {},
            "PCHEM_LEVS": {},
            "PROD1": {},
            "PROD2": {},
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

        PROD1 = QuantityFactory.zeros(
            self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a"
        )
        safe_assign_array(PROD1.view[:, :, :], inputs["PROD1"])

        PROD2 = QuantityFactory.zeros(
            self.quantity_factory, dims=[X_DIM, Z_DIM], units="n/a"
        )
        safe_assign_array(PROD2.view[:, :], inputs["PROD2"])

        PCHEM_LEVS = QuantityFactory.zeros(
            self.quantity_factory, dims=[Z_DIM], units="n/a"
        )
        safe_assign_array(PCHEM_LEVS.view[:], inputs["PCHEM_LEVS"])

        # FloatFields
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

        def mapl_interp(OX, IY, IX):
            """
            Python version of INTERP_LIN_0011_1.
            """
            # Find the interval index J such that IX[J] <= OX <= IX[J+1]
            J = min(max(np.count_nonzero(IX <= OX), 1), IX.size - 1) - 1
            # print(J)

            # Linear interpolation
            if IX[J + 1] != IX[J]:
                OY = IY[J] + ((OX - IX[J]) / (IX[J + 1] - IX[J])) * (IY[J + 1] - IY[J])

            else:
                OY = IY[J]

            return OY

        def mapl_interp_1111_1(OX_list, IY, IX):
            """
            Python version of INTERP_LIN_1111_1.
            """
            OY = []
            for ox in OX_list:
                oy = mapl_interp(ox, IY, IX)
                OY.append(oy)

            return OY

        def interp_no_extrap(OX_list, IY, IX):
            """
            Python version of INTERP_NO_EXTRAP.
            """
            max_index = len(IX)

            OY = mapl_interp_1111_1(OX_list, IY, IX)

            # If OX is below the first IX, clamp to the first value
            for i in range(0, max_index):
                if OX_list[i] <= IX[0]:
                    OY[i] = IY[0]

                if OX_list[i] >= IX[max_index - 1]:
                    OY[i] = IY[max_index - 1]

            return OY

        n_lat, n_lon = 24, 24

        for i in range(n_lat):
            for j in range(n_lon):
                PROD_INT.field[i, j, :] = interp_no_extrap(
                    OX_list=PL.field[i, j, :],
                    IY=PROD1.field[i, j, :],
                    IX=PCHEM_LEVS.field[:],
                )

        print("Done")

        return {
            "PROD_INT": PROD_INT.field[:],
        }
