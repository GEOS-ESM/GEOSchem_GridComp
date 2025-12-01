from ndsl import Quantity, StencilFactory
from ndsl.dsl.typing import Float, FloatFieldK
from typing import Sequence
import ndsl.xumpy as xp
import dace

from dace.sdfg.analysis.schedule_tree import treenodes as tn
import gt4py.cartesian.gtc.dace.treeir as tir


class PChemSpeciesInterpolate:
    def __init__(
        self,
        stencil_factory: StencilFactory,
        compute_dims: Sequence[str],
    ) -> None:
        self._grid_indexing = stencil_factory.grid_indexing
        self._is_cpu = not stencil_factory.config.is_gpu_backend
        self._is_orch = stencil_factory.config.dace_config.is_dace_orchestrated()
        self._dims = compute_dims

    def _mapl_interp(self, OX, IY, IX):
        """
        Python version of INTERP_LIN_0011_1.
        """
        # Find the interval index J such that IX[J] <= OX <= IX[J+1]
        J = (
            min(
                max(xp.count_nonzero(IX <= OX), 1),
                IX.size - 1,
            )
            - 1
        )

        # Linear interpolation
        if IX[J + 1] != IX[J]:
            OY = IY[J] + ((OX - IX[J]) / (IX[J + 1] - IX[J])) * (IY[J + 1] - IY[J])

        else:
            OY = IY[J]

        return OY

    def _mapl_interp_1111_1(self, OX, OY, IY, IX, dim_size):
        """
        Python version of INTERP_LIN_1111_1.
        """
        for i in range(dim_size):
            OY[i] = self._mapl_interp(OX[i], IY, IX)

    def _interpolate_no_extrapolate(
        self,
        OX: FloatFieldK,
        OY: FloatFieldK,
        Odim: int,
        IY: FloatFieldK,
        IX: FloatFieldK,
        Idim: int,
    ):
        self._mapl_interp_1111_1(OX, OY, IY, IX, Odim)

        # If OX is below the first IX, clamp to the first value
        for i in range(Odim):
            if OX[i] <= IX[0]:
                OY[i] = IY[0]

            if OX[i] >= IX[Idim - 1]:
                OY[i] = IY[Idim - 1]

    def __call__(
        self,
        mncv_offgrid: Quantity,
        species_index: int,
        fac: Float,
        lats_IJ: Quantity,
        pchem_lats_K_offgrid: Quantity,
        pchem_levs_K: Quantity,
        prod: Quantity,
        pl: Quantity,
        prod_int: Quantity,
    ) -> None:
        species_offgrid = mncv_offgrid[:, :, species_index - 1, 0] * fac + mncv_offgrid[
            :, :, species_index - 1, 1
        ] * (1.0 - fac)

        for j in range(self._grid_indexing.domain[1]):
            for k in range(self._grid_indexing.domain[2]):
                self._interpolate_no_extrapolate(
                    OX=lats_IJ.field[:, j],
                    OY=prod.field[:, j, k],
                    Odim=self._grid_indexing.domain[0],
                    IY=species_offgrid[:, k],
                    IX=pchem_lats_K_offgrid[:],
                    Idim=pchem_lats_K_offgrid.shape[0],
                )
            for i in range(self._grid_indexing.domain[0]):
                self._interpolate_no_extrapolate(
                    OX=pl.field[i, j, :],
                    OY=prod_int.field[i, j, :],
                    Odim=self._grid_indexing.domain[2],
                    IY=prod.field[i, j, :],
                    IX=pchem_levs_K.field[:],
                    Idim=self._grid_indexing.domain[2],
                )

    def __sdfg__(
        self,
        mncv_offgrid: dace.data.Array,
        species_index: int,
        fac: Float,
        lats_IJ: dace.data.Array,
        pchem_lats_K_offgrid: dace.data.Array,
        pchem_levs_K: dace.data.Array,
        prod: dace.data.Array,
        pl: dace.data.Array,
        prod_int: dace.data.Array,
    ):
        # TODO: In progress work
        return
        containers = {
            "mncv_offgrid": mncv_offgrid,
            "lats_IJ": lats_IJ,
            "pchem_lats_K_offgrid": pchem_lats_K_offgrid,
            "pchem_levs_K": pchem_levs_K,
            "prod": prod,
            "pl": pl,
            "prod_int": prod_int,
        }
        symbols = {
            axis.domain_symbol(): dace.dtypes.int32 for axis in tir.Axis.dims_3d()
        } | {axis.iteration_symbol(): dace.dtypes.int32 for axis in tir.Axis.dims_3d()}

        origin = self._grid_indexing.origin
        domain = self._grid_indexing.domain

        # Tasklet (_mapl_interp)
        code = (
            "J = (min(max(xp.count_nonzero(IX <= OX), 1), IX.size - 1) - 1)"
            "if IX[J + 1] != IX[J]:"
            "    OY = IY[J] + ((OX - IX[J]) / (IX[J + 1] - IX[J])) * (IY[J + 1] - IY[J])"
            "else:"
            "    OY = IY[J]"
        )

        tasklet_input = {
            "OX": dace.Memlet.from_array("in_buffer_3D", in_buffer_3D),
            "IX": "",
            "IY": "",
        }

        tasklet_output = {
            "OY": "",
        }

        dace_tasklet = dace.nodes.Tasklet(
            label=f"_mapl_interp__execution",
            code=code,
            inputs=list(tasklet_input.keys()),
            outputs=list(tasklet_output.keys()),
        )
        tasklet = tn.TaskletNode(
            node=dace_tasklet,
            in_memlets=tasklet_input,
            out_memlets=tasklet_output,
        )

        # Deep Interp map
        dace_map = dace.nodes.Map(
            label=f"interp_loop_{id(self)}",
            params=["idx_interp"],
            ndrange=dace.subsets.Range(
                [
                    # -1 because range bounds are inclusive
                    (f"{origin[2]}", f"{domain[2]} - 1", 1),
                ]
            ),
            schedule=dace.dtypes.ScheduleType.CPU_Multicore
            if self._is_cpu
            else dace.dtypes.ScheduleType.GPU_Device,
        )
        vertical_map = tn.MapScope(
            node=dace.nodes.MapEntry(dace_map), children=[tasklet]
        )

        # IJ sequential map
        dace_map = dace.nodes.Map(
            label=f"horizontal_loop_{id(self)}",
            params=[axis.iteration_symbol() for axis in tir.Axis.dims_horizontal()],
            ndrange=dace.subsets.Range(
                [
                    # -1 because range bounds are inclusive
                    (f"{origin[0]}", f"{domain[0]} - 1", 1),
                    (f"{origin[1]}", f"{domain[1]} - 1", 1),
                ]
            ),
            schedule=dace.dtypes.ScheduleType.Sequential,
        )
        horizontal_map = tn.MapScope(
            node=dace.nodes.MapEntry(dace_map), children=[vertical_map]
        )

        tree = tn.ScheduleTreeRoot(
            name=self.__class__.__name__,
            containers=containers,
            arg_names=list(containers.keys()),
            symbols=symbols,
            constants={},
            children=[horizontal_map],
        )

        sdfg = tree.as_sdfg(validate=True)
        return sdfg
