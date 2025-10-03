from gt4py.cartesian.gtscript import (
    PARALLEL,
    computation,
    interval,
    FORWARD,
    int32,
    K,
)
import gt4py.cartesian.gtscript as gtscript
from ndsl.dsl.typing import (
    Float,
    FloatField,
    FloatFieldK,
)
from ndsl import QuantityFactory, StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM


@gtscript.function
def interp(
    OX: Float,
    IY: FloatFieldK,  # FloatFieldIK
    IX: FloatFieldK,
) -> Float:

    sum = 0
    lev = 0
    while lev <= 71:
        if IX.at(K=lev) <= OX:
            sum += 1  # count how many IX <= OX
        lev += 1

    J = max(min(sum, 72 - 1), 1) - 1

    if IX.at(K=J + 1) != IX.at(K=J):
        OY = IY.at(K=J) + ((OX - IX.at(K=J)) / (IX.at(K=J + 1) - IX.at(K=J))) * (
            IY.at(K=J + 1) - IY.at(K=J)
        )
    else:
        OY = IY.at(K=J)

    return OY


@gtscript.function
def interp_no_extrap(
    OX: FloatFieldK,
    IY: FloatFieldK,  # FloatFieldIK
    IX: FloatFieldK,
) -> Float:
    # max_index = size(IX)
    max_index = 72

    OY = interp(OX, IY, IX)

    if OX < IX.at(K=0):
        OY = IY.at(K=0)

    if OX > IX.at(K=max_index - 1):
        OY = IY.at(K=max_index - 1)

    return OY


def test_interp_no_extrap(
    OY: FloatField,
    OX: FloatFieldK,
    IY: FloatFieldK,  # FloatFieldIK
    IX: FloatFieldK,
):
    with computation(PARALLEL), interval(...):
        max_index = 72

        sum = 0
        lev = 0
        while lev <= 71:
            if IX.at(K=lev) <= OX:
                sum += 1  # count how many IX <= OX
            lev += 1

        J = max(min(sum, 72 - 1), 1) - 1

        if IX.at(K=J + 1) != IX.at(K=J):
            OY = IY.at(K=J) + ((OX - IX.at(K=J)) / (IX.at(K=J + 1) - IX.at(K=J))) * (
                IY.at(K=J + 1) - IY.at(K=J)
            )
        else:
            OY = IY.at(K=J)

        if OX < IX.at(K=0):
            OY = IY.at(K=0)

        if OX > IX.at(K=max_index - 1):
            OY = IY.at(K=max_index - 1)


class Interp:
    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
    ) -> None:

        self.stencil_factory = stencil_factory
        self.quantity_factory = quantity_factory

        self._test_interp_no_extrap = self.stencil_factory.from_dims_halo(
            func=test_interp_no_extrap,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

    def __call__(
        self,
        PROD_INT: FloatField,
        PL: FloatFieldK,
        PROD: FloatFieldK,  # FloatFieldIK
        PCHEM_LEVS: FloatFieldK,
    ):

        self._test_interp_no_extrap(
            # Output
            OY=PROD_INT,
            # Inputs
            OX=PL,
            IY=PROD,
            IX=PCHEM_LEVS,
        )
