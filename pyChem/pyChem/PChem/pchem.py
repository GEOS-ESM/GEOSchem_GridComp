import pyChem.PChem.stencils as stencils
from ndsl import Quantity, QuantityFactory, StencilFactory, NDSLRuntime
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import Float, FloatField, FloatFieldIJ, Int
from pyChem.PChem.config import PChemConfiguration
from pyChem.PChem.temporaries import Temporaries
from pyChem.PChem.species_interpolate import PChemSpeciesInterpolate
from ndsl.stencils.basic_operations import copy_defn


class PChem(NDSLRuntime):
    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
        pchem_config: PChemConfiguration,
    ) -> None:
        super().__init__(stencil_factory)

        self.pchem_config = pchem_config
        self.temporaries = Temporaries.make(quantity_factory)
        self.stencil_factory = stencil_factory
        self.quantity_factory = quantity_factory

        # Raise an error if clim_years > 1 and tau < 0
        if self.pchem_config.clim_years > 1 and self.pchem_config.tau < 0:
            raise NotImplementedError(
                "Cannot run PCHEM in P & L mode with climYears > 1."
            )
        # Raise error if UPDATE_H2O_PL is triggered
        if self.pchem_config.USE_H2O_ProdLoss != 0:
            raise NotImplementedError("Warning: This code has not been ported!!")

        if self.pchem_config.tau < 0.0:
            raise NotImplementedError("Warning: This code has not been ported!!")

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

        self._init_pl = self.stencil_factory.from_dims_halo(
            func=stencils.init_pl, compute_dims=[X_DIM, Y_DIM, Z_DIM]
        )

        self._interpolate = PChemSpeciesInterpolate(
            stencil_factory, [X_DIM, Y_DIM, Z_DIM]
        )

        self._copy = self.stencil_factory.from_dims_halo(
            copy_defn, [X_DIM, Y_DIM, Z_DIM]
        )

    def _update_species(
        self,
        species_in,
        species_index,
        species_out,
        mncv,
        fac,
        lats,
        pchem_lats,
        pchem_levs,
        delp,
        pcrit,
        dt,
        tropp,
    ):
        self._update1(
            # In
            XX_in=species_in,
            # Out
            XX=self.temporaries.XX,
        )

        if self.pchem_config.tau > 0.0:
            self._interpolate(
                mncv_offgrid=mncv,
                species_index=species_index,
                fac=fac,
                lats_IJ=lats,
                pchem_lats_K_offgrid=pchem_lats,
                pchem_levs_K=pchem_levs,
                prod=self.temporaries.PROD,
                pl=self.temporaries.PL,
                prod_int=self.temporaries.PROD_INT,
            )

        self._update2(
            # In
            NN=species_index,
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

        self._copy(self.temporaries.XX, species_out)

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
        XX_CH4_in: Quantity,
        XX_N2O_in: Quantity,
        XX_CFC11_in: Quantity,
        XX_CFC12_in: Quantity,
        XX_HCFC22_in: Quantity,
        XX_OX_in: Quantity,
        XX_H2O_in: Quantity,
        AOA_in: FloatField,
        O3_pointer: Int,
        O3PPMV_pointer: Int,
        TO3_pointer: Int,
        TTO3_pointer: Int,
        O3VMR: FloatFieldIJ,
        OX: Quantity,
        ZTH: FloatFieldIJ,
        # Out
        CH4: Quantity,
        N2O: Quantity,
        CFC11: Quantity,
        CFC12: Quantity,
        HCFC22: Quantity,
        H2O: Quantity,
        O3: FloatField,
        O3PPMV: FloatField,
        AOA: FloatField,
    ) -> None:
        """
        PChem Driver

        Updates 7 chemical species (CH4, N2O, CFC11, CFC12, HCFC22, OX, H2O) and O3 based on
        production and loss rates.
        """
        self._init_pl(self.temporaries.PL, ple)

        self._update_species(
            species_index=NN_CH4,
            species_in=XX_CH4_in,
            species_out=CH4,
            mncv=mncv,
            fac=fac,
            lats=lats,
            pchem_lats=pchem_lats,
            pchem_levs=pchem_levs,
            delp=delp,
            pcrit=pcrit,
            dt=dt,
            tropp=tropp,
        )
        self._update_species(
            species_index=NN_N2O,
            species_in=XX_N2O_in,
            species_out=N2O,
            mncv=mncv,
            fac=fac,
            lats=lats,
            pchem_lats=pchem_lats,
            pchem_levs=pchem_levs,
            delp=delp,
            pcrit=pcrit,
            dt=dt,
            tropp=tropp,
        )
        self._update_species(
            species_index=NN_CFC11,
            species_in=XX_CFC11_in,
            species_out=CFC11,
            mncv=mncv,
            fac=fac,
            lats=lats,
            pchem_lats=pchem_lats,
            pchem_levs=pchem_levs,
            delp=delp,
            pcrit=pcrit,
            dt=dt,
            tropp=tropp,
        )
        self._update_species(
            species_index=NN_CFC12,
            species_in=XX_CFC12_in,
            species_out=CFC12,
            mncv=mncv,
            fac=fac,
            lats=lats,
            pchem_lats=pchem_lats,
            pchem_levs=pchem_levs,
            delp=delp,
            pcrit=pcrit,
            dt=dt,
            tropp=tropp,
        )
        self._update_species(
            species_index=NN_HCFC22,
            species_in=XX_HCFC22_in,
            species_out=HCFC22,
            mncv=mncv,
            fac=fac,
            lats=lats,
            pchem_lats=pchem_lats,
            pchem_levs=pchem_levs,
            delp=delp,
            pcrit=pcrit,
            dt=dt,
            tropp=tropp,
        )
        self._update_species(
            species_index=NN_OX,
            species_in=XX_OX_in,
            species_out=OX,
            mncv=mncv,
            fac=fac,
            lats=lats,
            pchem_lats=pchem_lats,
            pchem_levs=pchem_levs,
            delp=delp,
            pcrit=pcrit,
            dt=dt,
            tropp=tropp,
        )
        self._update_species(
            species_index=NN_H2O,
            species_in=XX_H2O_in,
            species_out=H2O,
            mncv=mncv,
            fac=fac,
            lats=lats,
            pchem_lats=pchem_lats,
            pchem_levs=pchem_levs,
            delp=delp,
            pcrit=pcrit,
            dt=dt,
            tropp=tropp,
        )

        if TO3_pointer or TTO3_pointer == 1:
            raise NotImplementedError("Warning: This code has not been ported!!")

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
