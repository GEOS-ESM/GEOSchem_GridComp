from dataclasses import dataclass

from ndsl.dsl.typing import Float, Int


@dataclass
class PChemConfiguration:
    clim_years: Int
    tau: Float
    USE_H2O_ProdLoss: Int
