from gt4py.cartesian.gtscript import FORWARD, PARALLEL, K, computation, exp, float32, interval, log10

import pyChem.constants as constants
from ndsl.dsl.typing import Float, FloatField, FloatFieldIJ, Int


def update1(
    # In
    XX_in: FloatField,
    # Out
    XX: FloatField,
):
    with computation(PARALLEL), interval(...):
        # Load in concentration levels
        XX = XX_in


def update2(
    # In
    NN: Int,
    tau: Float,
    pl: FloatField,
    delp: Float,
    pcrit: Float,
    prod_int: FloatField,
    dt: Float,
    tropp: FloatFieldIJ,
    # In/Out
    XX: FloatField,
):
    with computation(FORWARD), interval(...):
        if NN == 7:
            pcrit_H2O = float32(20000)

            if tropp == constants.MAPL_UNDEF:
                wrk_H2O = pcrit_H2O
            else:
                wrk_H2O = tropp

            wrk_H2O = min(wrk_H2O, pcrit_H2O)
            loss_int = (1.0 / tau) * max(min((wrk_H2O - pl) / delp, 1.0), 0.0)

    with computation(PARALLEL), interval(...):
        if NN != 7:
            loss_int = 1.0 / tau * max(min((pcrit - pl) / delp, 1.0), 0.0)

        prod_int = loss_int * prod_int

        XX = (XX + dt * prod_int) / (1.0 + dt * loss_int)


def update_ozone(
    # In
    O3_pointer: Int,
    O3PPMV_pointer: Int,
    O3VMR: FloatFieldIJ,
    OX: FloatField,
    PL: FloatField,
    ZTH: FloatFieldIJ,
    # Out
    O3: FloatField,
    O3PPMV: FloatField,
):
    with computation(FORWARD), interval(...):
        if PL < 100.0 and ZTH > 0.0:
            O3VMR = OX * exp(-1.5 * (log10(PL) - 2.0) ** 2)
        else:
            O3VMR = OX

        if O3_pointer == 1:
            O3 = O3VMR * (constants.MAPL_O3MW / constants.MAPL_AIRMW)
        if O3PPMV_pointer == 1:
            O3PPMV = O3VMR * 1.0e6


def age_of_air(
    AOA_in: FloatField,
    dt: Float,
    AOA: FloatField,
) -> None:
    with computation(PARALLEL), interval(0, -1):
        AOA = AOA_in + (dt / 86400.0)
    with computation(PARALLEL), interval(-1, None):
        AOA = 0.0


def init_pl(
    pl: FloatField,
    ple: FloatField,
) -> None:
    with computation(PARALLEL), interval(...):
        pl = 0.5 * (ple + ple[K + 1])
