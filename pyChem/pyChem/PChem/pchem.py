from ndsl import QuantityFactory, StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import (
    Float,
    FloatField,
    Int,
    FloatFieldIJ,
)

from pyChem.PChem.config import PChemConfiguration
import pyChem.PChem.stencils as stencils
import pyChem.Shared.interp_no_extrap as interp
import numpy as np


class PChem:
    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
        pchem_config: PChemConfiguration,
    ) -> None:

        self.pchem_config = pchem_config
        # self.temporaries = Temporaries.make(quantity_factory)
        self.stencil_factory = stencil_factory
        self.quantity_factory = quantity_factory

        # Raise an error if clim_years > 1 and tau < 0
        # if self.pchem_config.clim_years > 1 and self.pchem_config.tau < 0:
        #     raise NotImplementedError(
        #         f"Cannot run PCHEM in P & L mode with climYears > 1."
        #     )

        self._update1 = self.stencil_factory.from_dims_halo(
            func=stencils.update1,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self._update2 = self.stencil_factory.from_dims_halo(
            func=stencils.update2,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

    def __call__(
        self,
        # Inputs
        NN: Int,
        ple: FloatField,
        XX_in: FloatField,
        delp: Float,
        fac: Float,
        mncv: FloatField,
        # mnpl: FloatField,
        pcrit: Float,
        lats: FloatFieldIJ,
        pchem_levs: FloatField,
        pchem_lats: FloatField,
        dt: Float,
        tau: Float,
        prod: FloatField,
        prod_int: FloatField,
        pl: FloatField,
        tropp: FloatFieldIJ,
        # loss: FloatField,
        # loss_int: FloatField,
        # Outputs
        XX: FloatField,
    ):
        lm = 73
        pl.view[:, :, :] = 0.5 * (ple.view[:, :, 0 : lm - 1] + ple.view[:, :, 1:lm])

        self._update1(
            # Inputs
            XX_in=XX_in,
            # Outputs
            XX=XX,
        )

        im, jm, n_levs = 24, 24, 72

        # if tau <= 0.0:
        #     prod1 = mnpl[:, :, NN, 1, 1] * fac + mnpl[:, :, NN, 1, 2] * (1.0 - fac)
        #     loss1 = mnpl[:, :, NN, 2, 1] * fac + mnpl[:, :, NN, 2, 2] * (1.0 - fac)

        #     for j in range(jm):
        #         for k in range(n_levs):
        #             prod[:, j, k] = interp.interp_no_extrap(
        #                 OX_list=lats[:, j],
        #                 IY=prod1[:, k],
        #                 IX=pchem_lats[:],
        #             )
        #             loss[:, j, k] = interp.interp_no_extrap(
        #                 OX_list=lats[:, j],
        #                 IY=loss1[:, k],
        #                 IX=pchem_lats[:],
        #             )
        #         for i in range(im):
        #             prod_int[i, j, :] = interp.interp_no_extrap(
        #                 OX_list=pl[i, j, :],
        #                 IY=prod[i, j, :],
        #                 IX=pchem_levs[:],
        #             )
        #             loss_int[i, j, :] = interp.interp_no_extrap(
        #                 OX_list=pl[i, j, :],
        #                 IY=loss[i, j, :],
        #                 IX=pchem_levs[:],
        #             )

        #     XX = (XX + dt * prod_int) / (1.0 + dt * loss_int)

        if tau > 0.0:
            prod1 = mncv[:, :, NN - 1, 0] * fac + mncv[:, :, NN - 1, 1] * (1.0 - fac)

            for j in range(jm):
                for k in range(n_levs):
                    prod.view[:, j, k] = interp.interp_no_extrap(
                        OX_list=lats.view[:, j],
                        IY=prod1[:, k],
                        IX=pchem_lats[:],
                    )
                for i in range(im):
                    prod_int.view[i, j, :] = interp.interp_no_extrap(
                        OX_list=pl.view[i, j, :],
                        IY=prod.view[i, j, :],
                        IX=pchem_levs.view[:],
                    )

        self._update2(
            # Inputs
            NN=NN,
            tau=tau,
            pl=pl,
            delp=delp,
            pcrit=pcrit,
            prod_int=prod_int,
            dt=dt,
            tropp=tropp,
            # In/Outs
            XX=XX,
        )
