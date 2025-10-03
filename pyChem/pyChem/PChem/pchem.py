from ndsl import QuantityFactory, StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import (
    Float,
    FloatField,
    Int,
)

from pyChem.PChem.config import PChemConfiguration
import pyChem.PChem.stencils as stencils


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
        if self.pchem_config.clim_years > 1 and self.pchem_config.tau < 0:
            raise NotImplementedError(
                f"Cannot run PCHEM in P & L mode with climYears > 1."
            )

        self._update = self.stencil_factory.from_dims_halo(
            func=stencils.update,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

    def __call__(
        self,
        # Inputs
        NN: Int,
        pl: FloatField,
        XX_in: FloatField,
        delp: Float,
        fac: Float,
        # mncv: FloatField,
        pcrit: Float,
        prod: FloatField,
        prod_int: FloatField,
        dt: Float,
        # Outputs
        XX: FloatField,
    ):

        self._update(
            # Inputs
            NN=NN,
            pl=pl,
            XX_in=XX_in,
            delp=delp,
            fac=fac,
            # mncv=mncv,
            pcrit=pcrit,
            prod=prod,
            prod_int=prod_int,
            dt=dt,
            tau=self.pchem_config.tau,
            # Outputs
            XX=XX,
        )
