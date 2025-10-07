from gt4py.cartesian.gtscript import (
    PARALLEL,
    computation,
    interval,
)
import gt4py.cartesian.gtscript as gtscript
from ndsl.dsl.typing import (
    Float,
    FloatField,
    Int,
    FloatFieldIJ,
)

import pyChem.constants as constants


def update1(
    # Inputs
    XX_in: FloatField,
    # Outputs
    XX: FloatField,
):
    """
    Function that updates species (e.g., CH4, N20, CFC 11, etc.) accoridng to tabulated production
    and loss rates.

    TODO: This was a first attemp at porting the 'update' subroutine of PChem. There is still much
    to do. Anything commented out needs to still be ported. MNCV (shape=[91, 72, 7, 2]) needs to be
    revisited. There is also alot of MAPL calls that need to be revisited. interp_no_extrap needs
    to be ported as well.
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
    # Inputs
    NN: Int,
    tau: Float,
    pl: FloatField,
    delp: Float,
    pcrit: Float,
    prod_int: FloatField,
    dt: Float,
    tropp: FloatFieldIJ,
    # In/Outs
    XX: FloatField,
):
    with computation(PARALLEL), interval(...):
        # call MAPL_GetResource(MAPL, DELP,  LABEL=trim(NAME)//"_DELP:" , DEFAULT=5000. ,RC=STATUS)
        # VERIFY_(STATUS)
        # delp = max(delp, 1.0e-16)

        if NN == 7:
            # call MAPL_GetResource(MAPL, PCRIT, LABEL=trim(NAME)//"_PCRIT:", DEFAULT=20000. ,RC=STATUS)
            # VERIFY_(STATUS)
            # allocate(WRK(IM,JM),stat=STATUS)
            # VERIFY_(STATUS)
            if tropp == constants.MAPL_UNDEF:
                wrk = pcrit
            else:
                wrk = tropp

            wrk = 0.9  # min(wrk, pcrit)

            loss_int = (1.0 / tau) * max(min(wrk - pl / delp, 1.0), 0.0)

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
