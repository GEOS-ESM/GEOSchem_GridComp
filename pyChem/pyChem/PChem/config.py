from dataclasses import dataclass

from ndsl.dsl.typing import Int, Float


@dataclass
class PChemConfiguration:
    clim_years: Int
    tau: Float
