from dataclasses import dataclass

from ndsl import Quantity, QuantityFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM


@dataclass
class Temporaries:
    XX: Quantity
    PL: Quantity
    PROD_INT: Quantity
    PROD: Quantity

    @classmethod
    def make(cls, quantity_factory: QuantityFactory):
        # FloatFields
        XX = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        PL = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        PROD_INT = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        PROD = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")

        return cls(
            XX,
            PL,
            PROD_INT,
            PROD,
        )
