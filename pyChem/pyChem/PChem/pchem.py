from ndsl import QuantityFactory, StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import (
    Float,
    FloatField,
    Int,
    FloatFieldIJ,
)

from pyChem.PChem.config import PChemConfiguration
from pyChem.PChem.temporaries import Temporaries
import pyChem.PChem.stencils as stencils
import pyChem.Shared.interp_no_extrap as interp


class PChem:
    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
        pchem_config: PChemConfiguration,
    ) -> None:

        self.pchem_config = pchem_config
        self.temporaries = Temporaries.make(quantity_factory)
        self.stencil_factory = stencil_factory
        self.quantity_factory = quantity_factory

        # Raise an error if clim_years > 1 and tau < 0
        if self.pchem_config.clim_years > 1 and self.pchem_config.tau < 0:
            raise NotImplementedError(
                f"Cannot run PCHEM in P & L mode with climYears > 1."
            )
        # Raise error if UPDATE_H2O_PL is triggered
        if self.pchem_config.USE_H2O_ProdLoss != 0:
            raise NotImplementedError(f"Warning: This code has not been ported!!")

        if self.pchem_config.tau < 0.0:
            raise NotImplementedError(f"Warning: This code has not been ported!!")

        self._update1 = self.stencil_factory.from_dims_halo(
            func=stencils.update1,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self._update2 = self.stencil_factory.from_dims_halo(
            func=stencils.update2,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self._update_ozone = self.stencil_factory.from_dims_halo(
            func=stencils.update_ozone,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self._age_of_air = self.stencil_factory.from_dims_halo(
            func=stencils.age_of_air,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

    def __call__(
        self,
        # In
        NN_CH4: Int,
        NN_N2O: Int,
        NN_CFC11: Int,
        NN_CFC12: Int,
        NN_HCFC22: Int,
        NN_OX: Int,
        NN_H2O: Int,
        ple: FloatField,
        delp: Float,
        fac: Float,
        mncv: FloatField,
        pcrit: Float,
        lats: FloatFieldIJ,
        pchem_levs: FloatField,
        pchem_lats: FloatField,
        dt: Float,
        tropp: FloatFieldIJ,
        XX_CH4_in: FloatField,
        XX_N2O_in: FloatField,
        XX_CFC11_in: FloatField,
        XX_CFC12_in: FloatField,
        XX_HCFC22_in: FloatField,
        XX_OX_in: FloatField,
        XX_H2O_in: FloatField,
        AOA_in: FloatField,
        O3_pointer: Int,
        O3PPMV_pointer: Int,
        TO3_pointer: Int,
        TTO3_pointer: Int,
        O3VMR: FloatFieldIJ,
        OX: FloatField,
        ZTH: FloatFieldIJ,
        # Out
        CH4: FloatField,
        N2O: FloatField,
        CFC11: FloatField,
        CFC12: FloatField,
        HCFC22: FloatField,
        H2O: FloatField,
        O3: FloatField,
        O3PPMV: FloatField,
        AOA: FloatField,
    ):
        """
        PChem Driver

        Updates 7 chemical species (CH4, N2O, CFC11, CFC12, HCFC22, OX, H2O) and O3 based on
        production and loss rates.
        """
        lm = 73
        self.temporaries.PL.view[:, :, :] = 0.5 * (
            ple.view[:, :, 0 : lm - 1] + ple.view[:, :, 1:lm]
        )

        state = [NN_CH4, NN_N2O, NN_CFC11, NN_CFC12, NN_HCFC22, NN_OX, NN_H2O]
        species = [
            XX_CH4_in,
            XX_N2O_in,
            XX_CFC11_in,
            XX_CFC12_in,
            XX_HCFC22_in,
            XX_OX_in,
            XX_H2O_in,
        ]
        XX_list = []

        for NN, XX_in in zip(state, species):

            self._update1(
                # In
                XX_in=XX_in,
                # Out
                XX=self.temporaries.XX,
            )

            im, jm, n_levs = 24, 24, 72

            if self.pchem_config.tau > 0.0:
                prod1 = mncv[:, :, NN - 1, 0] * fac + mncv[:, :, NN - 1, 1] * (
                    1.0 - fac
                )

                for j in range(jm):
                    for k in range(n_levs):
                        self.temporaries.PROD.view[:, j, k] = interp.interp_no_extrap(
                            OX_list=lats.view[:, j],
                            IY=prod1[:, k],
                            IX=pchem_lats[:],
                        )
                    for i in range(im):
                        self.temporaries.PROD_INT.view[i, j, :] = (
                            interp.interp_no_extrap(
                                OX_list=self.temporaries.PL.view[i, j, :],
                                IY=self.temporaries.PROD.view[i, j, :],
                                IX=pchem_levs.view[:],
                            )
                        )

            self._update2(
                # In
                NN=NN,
                tau=self.pchem_config.tau,
                pl=self.temporaries.PL,
                delp=delp,
                pcrit=pcrit,
                prod_int=self.temporaries.PROD_INT,
                dt=dt,
                tropp=tropp,
                # In/Out
                XX=self.temporaries.XX,
            )

            XX_list.append(self.temporaries.XX.view[:].copy())

        species_list = [CH4, N2O, CFC11, CFC12, HCFC22, OX, H2O]

        for species, XX in zip(species_list, XX_list):
            species.view[:] = XX

        if TO3_pointer or TTO3_pointer == 1:
            raise NotImplementedError(f"Warning: This code has not been ported!!")

        self._update_ozone(
            # In
            O3_pointer=O3_pointer,
            O3PPMV_pointer=O3PPMV_pointer,
            O3VMR=O3VMR,
            OX=OX,
            PL=self.temporaries.PL,
            ZTH=ZTH,
            # Out
            O3=O3,
            O3PPMV=O3PPMV,
        )

        self._age_of_air(
            # In
            AOA_in=AOA_in,
            dt=dt,
            # Out
            AOA=AOA,
        )
