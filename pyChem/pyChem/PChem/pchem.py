from gt4py.cartesian.gtscript import (
    PARALLEL,
    computation,
    interval,
)

from ndsl import QuantityFactory, StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import (
    Float,
    FloatField,
    Int,
)


def update(
    # Inputs
    NN: Int,
    pl: FloatField,
    XX_in: FloatField,
    delp: Float,
    fac: Float,
    mncv: FloatField,  # Fix this (shape=[91,72,7,2])
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

    TODO: This was a first attemp at porting the 'update' subroutine of PChem. There is still much to
    do. Anything commented out needs to still be ported. MNCV (shape=[91, 72, 7, 2]) needs to be
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
        # if (TAU<=0.0) then  ! By convention this is the prod(index 1) and loss(index 2) case

        #     _ASSERT(trim(NAME)/="H2O",'needs informative message')

        #     PROD1 = PCHEM_STATE%MNPL(:,:,NN,1,1)*FAC + PCHEM_STATE%MNPL(:,:,NN,1,2)*(1.-FAC)
        #     LOSS1 = PCHEM_STATE%MNPL(:,:,NN,2,1)*FAC + PCHEM_STATE%MNPL(:,:,NN,2,2)*(1.-FAC)

        #     do j=1,jm
        #         do l=1,nlevs
        #             call INTERP_NO_EXTRAP( PROD(:,L), LATS(:,J), Prod1(:,L), PCHEM_STATE%LATS)
        #             call INTERP_NO_EXTRAP( LOSS(:,L), LATS(:,J), Loss1(:,L), PCHEM_STATE%LATS)
        #         enddo
        #         do i=1,im
        #             call INTERP_NO_EXTRAP( PROD_INT(i,j,:), PL(i,j,:), PROD(i,:), PCHEM_STATE%LEVS)
        #             call INTERP_NO_EXTRAP( LOSS_INT(i,j,:), PL(i,j,:), LOSS(i,:), PCHEM_STATE%LEVS)
        #         enddo
        #     end do

        #     XX = (XX + DT*PROD_INT) / (1.0 + DT*LOSS_INT)

        if tau > 0.0:

            prod1 = mncv.at(K=NN - 1, ddim=0) * fac + mncv.at(K=NN - 1, ddim=1) * (
                1.0 - fac
            )

            # INTERP_NO_EXTRAP( PROD(:,L), LATS(:,J), Prod1(:,L), PCHEM_STATE%LATS)
            # INTERP_NO_EXTRAP( PROD_INT(i,j,:), PL(i,j,:), PROD(i,:), PCHEM_STATE%LEVS)

            # call MAPL_GetResource(MAPL, DELP,  LABEL=trim(NAME)//"_DELP:" , DEFAULT=5000. ,RC=STATUS)
            # VERIFY_(STATUS)

            delp = max(delp, 1.0e-16)

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


class PChem:
    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
        # pchem_config: PChemConfiguration,
    ) -> None:

        # self.pchem_config = pchem_config
        # self.temporaries = Temporaries.make(quantity_factory)
        self.stencil_factory = stencil_factory
        self.quantity_factory = quantity_factory

        # Raise an error if clim_years > 1 and tau < 0
        # if self.pchem_config.clim_years > 1 and self.pchem_config.tau < 0:
        #     raise NotImplementedError(
        #         f"Cannot run PCHEM in P & L mode with climYears > 1."
        #     )

        self._update = self.stencil_factory.from_dims_halo(
            func=update,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

    def __call__(
        self,
        # Inputs
        NN: Int,
        clim_years: Int,
        pl: FloatField,
        XX_in: FloatField,
        delp: Float,
        fac: Float,
        mncv: FloatField,
        pcrit: Float,
        prod: FloatField,
        prod_int: FloatField,
        tau: Float,
        dt: Float,
        # Outputs
        XX: FloatField,
    ):

        self._update(
            # Inputs
            NN=NN,
            pl=pl,
            XX_in=XX_in,
            delp=delp,
            fac=fac,
            mncv=mncv,
            pcrit=pcrit,
            prod=prod,
            prod_int=prod_int,
            tau=tau,
            dt=dt,
            # Outputs
            XX=XX,
        )
