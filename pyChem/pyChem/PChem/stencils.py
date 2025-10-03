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
    FloatFieldK,
)


def update(
    # Inputs
    NN: Int,
    pl: FloatField,
    XX_in: FloatField,
    delp: Float,
    fac: Float,
    # mncv: FloatField,  # Fix this (shape=[91,72,7,2])
    pcrit: Float,
    prod: FloatField,
    prod_int: FloatField,
    tau: Float,
    dt: Float,
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

    with computation(PARALLEL), interval(...):
        # if tau <= 0.0:  # By convention this is the prod(index 1) and loss(index 2) case
        #     prod1 = mnpl(:,:,NN,1,1)*fac + mnpl(:,:,NN,1,2)*(1.-fac)
        #     loss1 = mnpl(:,:,NN,2,1)*fac + mnpl(:,:,NN,2,2)*(1.-fac)

        #     prod = interp_no_extrap( prod(:,L), lat(:,J), prod1(:,L), pchem_lats)
        #     loss = interp_no_extrap( loss(:,L), lats(:,J), loss1(:,L), pchem_lats)

        #     prod_int = interp_no_extrap( prod_int(i,j,:), pl(i,j,:), prod(i,:), pchem_levs)
        #     loss_int = interp_no_extrap( loss_int(i,j,:), pl(i,j,:), loss(i,:), pchem_levs)

        #     XX = (XX + dt*prod_int) / (1.0 + dt*loss_int)

        if tau > 0.0:

            # prod1 = mncv.at(K=NN - 1, ddim=0) * fac + mncv.at(K=NN - 1, ddim=1) * (
            #     1.0 - fac
            # )

            # prod = interp_no_extrap( prod(:,L), lats(:,J), prod1(:,L), pchem_lats)
            # prod_int = interp_no_extrap( prod_int(i,j,:), pl(i,j,:), prod(i,:), pchem_levs)

            # call MAPL_GetResource(MAPL, DELP,  LABEL=trim(NAME)//"_DELP:" , DEFAULT=5000. ,RC=STATUS)
            # VERIFY_(STATUS)

            # delp = max(delp, 1.0e-16)

            # if NN == 7:
            # call MAPL_GetResource(MAPL, PCRIT, LABEL=trim(NAME)//"_PCRIT:", DEFAULT=20000. ,RC=STATUS)
            # VERIFY_(STATUS)
            # allocate(WRK(IM,JM),stat=STATUS)
            # VERIFY_(STATUS)
            # where (TROPP==MAPL_UNDEF)
            #     WRK = PCRIT
            # elsewhere
            #     WRK = TROPP

            # WRK = min(WRK, PCRIT)

            # LOSS_INT(:,:,L) = (1./TAU) * max( min( (WRK-PL(:,:,L))/DELP, 1.0), 0.0)

            # deallocate(WRK)
            if NN != 7:
                # call MAPL_GetResource(MAPL, PCRIT, LABEL=trim(NAME)//"_PCRIT:", DEFAULT=1.e+16 ,RC=STATUS)
                # VERIFY_(STATUS)
                loss_int = (1.0 / tau) * max(min((pcrit - pl) / delp, 1.0), 0.0)

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
