import numpy as np


def mapl_interp(OX, IY, IX):
    """
    Python version of INTERP_LIN_0011_1.
    """
    # Find the interval index J such that IX[J] <= OX <= IX[J+1]
    J = min(max(np.count_nonzero(IX <= OX), 1), IX.size - 1) - 1

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
    for i in range(0, len(OX_list)):
        if OX_list[i] <= IX[0]:
            OY[i] = IY[0]

        if OX_list[i] >= IX[max_index - 1]:
            OY[i] = IY[max_index - 1]

    return OY
