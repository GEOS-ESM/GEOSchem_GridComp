from ndsl import Namelist, Quantity, QuantityFactory, StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import Float, FloatField, Int
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from ndsl.utils import safe_assign_array
from pyChem.PChem.pchem import PChem
from pyChem.PChem.config import PChemConfiguration


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
            "PL": {},
            "XX_in": {},
            "mncv": {},
            "prod": {},
            "prod_int": {},
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
        self.pchem_config = PChemConfiguration(
            Int(inputs["clim_years"]), Float(inputs["tau"])
        )

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

        # Field inputs
        pl = QuantityFactory.zeros(
            self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a"
        )
        safe_assign_array(pl.view[:, :, :], inputs["PL"])
        xx_in = QuantityFactory.zeros(
            self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a"
        )
        safe_assign_array(xx_in.view[:, :, :], inputs["XX_in"])
        prod_int = QuantityFactory.zeros(
            self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a"
        )
        safe_assign_array(prod_int.view[:, :, :], inputs["prod_int"])
        prod = QuantityFactory.zeros(
            self.quantity_factory, dims=[X_DIM, Z_DIM], units="n/a"
        )
        safe_assign_array(prod.view[:, :], inputs["prod"])

        # FloatFields
        xx_CH4 = QuantityFactory.zeros(
            self.quantity_factory, dims=[X_DIM, Y_DIM, Z_DIM], units="n/a"
        )

        pchem(
            # Field inputs
            NN=NN_CH4,
            pl=pl,
            XX_in=xx_in,
            delp=delp,
            fac=fac,
            # mncv=mncv, Need to revisit this 4D field
            pcrit=pcrit,
            prod=prod,
            prod_int=prod_int,
            dt=dt,
            # Outputs
            XX=xx_CH4,
        )

        return {
            "XX_CH4": xx_CH4.view[:],
        }
