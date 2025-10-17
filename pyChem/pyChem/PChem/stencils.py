from gt4py.cartesian.gtscript import (
    PARALLEL,
    FORWARD,
    computation,
    interval,
    int32,
    float32,
    log10,
    exp,
    K,
)

from ndsl.dsl.typing import (
    Float,
    FloatField,
    Int,
    FloatFieldIJ,
)

import pyChem.constants as constants


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
):
    from __externals__ import k_end

    with computation(PARALLEL), interval(...):
        AOA = AOA_in + (dt / 86400.0)

        if K == k_end:
            AOA = 0.0
