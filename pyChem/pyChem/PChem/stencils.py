from gt4py.cartesian.gtscript import (
    PARALLEL,
    FORWARD,
    computation,
    interval,
    int32,
    float32,
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
    """
    Stencil that updates species (e.g., CH4, N20, CFC 11, etc.) accoridng to tabulated production
    and loss rates.

    TODO: This was a first attemp at porting the 'Update' portion of PChem. There is still more
    to do - MAPL calls need to be revisited.
    """

    with computation(PARALLEL), interval(...):
        # if (trim(NAME) == "H2O") then
        #    call MAPL_GetPointer ( IMPORT,   XX,  'Q', RC=STATUS )
        #    VERIFY_(STATUS)
        #    _ASSERT(associated(XX),'needs informative message')
        # else
        #    call MAPL_GetPointer ( INTERNAL, XX, NAME, RC=STATUS )
        #    VERIFY_(STATUS)
        #    _ASSERT(associated(XX),'needs informative message')
        #    call MAPL_GetResource(MAPL, VALUE,LABEL=trim(NAME)//"_FIXED_VALUE:", default=-1., RC=STATUS)
        #    VERIFY_(STATUS)
        #    if(VALUE>=0.0) then
        #       XX = VALUE
        #       return
        #    end if
        # endif

        XX = XX_in

        # call MAPL_GetResource(MAPL,   TAU,LABEL=trim(NAME)//"_RELAXTIME:", DEFAULT=0.0 ,RC=STATUS)
        # VERIFY_(STATUS)
        # call MAPL_GetPointer ( EXPORT, XX_PROD, trim(NAME)//'_PROD', RC=STATUS )
        # VERIFY_(STATUS)
        # call MAPL_GetPointer ( EXPORT, XX_LOSS, trim(NAME)//'_LOSS', RC=STATUS )
        # VERIFY_(STATUS)


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
        # call MAPL_GetResource(MAPL, DELP,  LABEL=trim(NAME)//"_DELP:" , DEFAULT=5000. ,RC=STATUS)
        # VERIFY_(STATUS)
        # delp = max(delp, 1.0e-16)

        if NN == 7:
            # call MAPL_GetResource(MAPL, PCRIT, LABEL=trim(NAME)//"_PCRIT:", DEFAULT=20000. ,RC=STATUS)
            # VERIFY_(STATUS)
            # allocate(WRK(IM,JM),stat=STATUS)
            # VERIFY_(STATUS)
            pcrit_H2O = float32(20000)

            if tropp == constants.MAPL_UNDEF:
                wrk_H2O = pcrit_H2O
            else:
                wrk_H2O = tropp

            wrk_H2O = min(wrk_H2O, pcrit_H2O)
            loss_int = (1.0 / tau) * max(min((wrk_H2O - pl) / delp, 1.0), 0.0)

    with computation(PARALLEL), interval(...):

        if NN != 7:
            # call MAPL_GetResource(MAPL, PCRIT, LABEL=trim(NAME)//"_PCRIT:", DEFAULT=1.e+16 ,RC=STATUS)
            # VERIFY_(STATUS)
            loss_int = 1.0 / tau * max(min((pcrit - pl) / delp, 1.0), 0.0)

        prod_int = loss_int * prod_int

        XX = (XX + dt * prod_int) / (1.0 + dt * loss_int)

        # if(associated(XX_PROD)) XX_PROD =  PROD_INT
        # if(associated(XX_LOSS)) XX_LOSS = -LOSS_INT*XX

        # if(trim(NAME)=='OX') then
        #    call MAPL_GetPointer ( EXPORT, OX_TEND, 'OX_TEND', RC=STATUS )
        #    VERIFY_(STATUS)
        #    if(associated(OX_TEND)) OX_TEND = (PROD_INT - LOSS_INT*XX)
        # end if

        # if(trim(NAME)=='H2O') then
        #    call MAPL_GetPointer ( EXPORT, H2O_TEND, 'H2O_TEND', RC=STATUS )
        #    VERIFY_(STATUS)
        #    if(associated(H2O_TEND)) H2O_TEND = (PROD_INT - LOSS_INT*XX)
        # end if
