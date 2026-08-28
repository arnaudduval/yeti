from .postprocessing import Postprocessing
from typing import Union, List, Tuple, Optional, Literal, Any, Dict
from matplotlib import pyplot as plt
from copy import deepcopy
import numpy as np
import logging

logger = logging.getLogger("SRC.IGA.POSTPROCESSING")

try:
    import pyvista as pv
except:
    logger.debug(
        "pyvista module is not installed. Some functions in post-processing may be disabled"
    )


class IgaPostprocessing:

    _ndim_default = 3

    @staticmethod
    def _set_points_to_evaluate(
        patch: Any, sample_size: Optional[Union[int, np.ndarray]]
    ) -> Tuple[np.ndarray, List[np.ndarray]]:
        if sample_size is None:
            knots_list = [quadrule.quadpts for quadrule in patch.quadrule_list]
            sample_size = np.copy(patch.nbqp)
        elif np.isscalar(sample_size):
            sample_size = int(np.real(sample_size))
            knots_list = [np.linspace(0.0, 1.0, sample_size)] * patch.ndim
            sample_size = sample_size * np.ones(patch.ndim, dtype=int)
        elif isinstance(sample_size, np.ndarray):
            knots_list = [
                np.linspace(0.0, 1.0, int(sample_size[i])) for i in range(patch.ndim)
            ]
        else:
            raise NotImplementedError()

        list_size = np.ones(IgaPostprocessing._ndim_default, dtype=int)
        for i in range(min(len(sample_size), len(list_size))):
            list_size[i] = sample_size[i]
        return list_size, knots_list

    @staticmethod
    def _extract_position(
        patch: Any,
        sample_size: np.ndarray,
        knots_list: Optional[List[np.ndarray]],
    ) -> Tuple[List[np.ndarray], np.ndarray, np.ndarray]:

        # Interpolate field
        _, knots_phy, det_jac, _ = patch.interpolate_field(
            knots_list=knots_list, eval_geometry=True
        )

        # Extract x, y, and z coordinates from physical points
        ndim = knots_phy.shape[0]
        x_data = knots_phy[0]
        y_data = np.zeros_like(x_data) if ndim < 2 else knots_phy[1]
        z_data = np.zeros_like(x_data) if ndim < 3 else knots_phy[2]

        # Reshape data into the required format for visualization
        position_data = [
            np.reshape(data, (sample_size[0], sample_size[1], -1), order="F")
            for data in [x_data, y_data, z_data]
        ]
        return position_data, det_jac, knots_phy

    @staticmethod
    def _update_data(
        field_name: str,
        field_values: np.ndarray,
        sample_size: np.ndarray,
        point_data: dict,
    ):
        nrows = field_values.shape[0]
        for idx_row in range(nrows):
            newfieldname = f"{field_name}{'_'}{idx_row+1}" if nrows > 1 else field_name
            newfielvalue = np.reshape(
                field_values[idx_row],
                (sample_size[0], sample_size[1], -1),
            )
            point_data.update({newfieldname: newfielvalue})

    @staticmethod
    def _save_vtk(
        path: str, position_data: List[np.ndarray], point_data: Dict[str, np.ndarray]
    ):
        assert isinstance(position_data, list) and isinstance(point_data, dict)
        grid = pv.StructuredGrid(position_data[0], position_data[1], position_data[2])
        logger.info(repr(grid))
        for name, array in point_data.items():
            grid.point_data[name] = array.flatten()
        grid.save(f"{path}.vts")

    @staticmethod
    def _export_primal_data(
        patch: Any,
        filename: str = "output",
        folder: Optional[str] = None,
        fields: dict = {},
        write_file: bool = True,
        warping: Optional[np.ndarray] = None,
        **extra_args,
    ):
        """Export solution in VTK format.
        It is possible to use Paraview to visualize data
        """
        if write_file:
            folder = Postprocessing.verify_folder_to_save(folder)

        # Get sample size and knots to evaluate
        sample_size = extra_args.get("sample_size")
        sample_size, knots_list = IgaPostprocessing._set_points_to_evaluate(
            patch, sample_size
        )

        # Get position data
        patch_copy = deepcopy(patch)
        if isinstance(warping, np.ndarray):
            patch_copy.ctrlpts += warping
        position_data, det_jac, pts_phy = IgaPostprocessing._extract_position(
            patch_copy,
            sample_size,
            knots_list,
        )

        # Create point data
        point_data = {}
        for fieldname, fieldvalue in fields.items():
            if isinstance(fieldvalue, np.ndarray):
                fieldinterp = patch_copy.interpolate_field(
                    knots_list=knots_list,
                    u_ctrlpts=np.atleast_2d(fieldvalue),
                    eval_geometry=False,
                )[0]
            elif callable(fieldvalue):
                extra_args.update({"position": pts_phy})
                fieldinterp = np.atleast_2d(np.asarray(fieldvalue(extra_args)))
            else:
                continue
            IgaPostprocessing._update_data(
                fieldname,
                fieldinterp,
                sample_size,
                point_data,
            )
        if fields.get("det_jac") is None:
            IgaPostprocessing._update_data(
                "det_jac",
                np.atleast_2d(det_jac),
                sample_size,
                point_data,
            )

        if write_file:
            IgaPostprocessing._save_vtk(
                f"{folder}/{filename}", position_data, point_data
            )

        return position_data, point_data

    @staticmethod
    def _export_dual_data(
        patch: Any,
        filename: str = "output",
        folder: Optional[str] = None,
        fields: dict = {},
        write_file: bool = True,
        warping: Optional[np.ndarray] = None,
        **extra_args,
    ):

        if write_file:
            folder = Postprocessing.verify_folder_to_save(folder)

        # Get sample size and knots to evaluate
        sample_size = IgaPostprocessing._set_points_to_evaluate(patch, None)[0]

        # Get position data
        patch_copy = deepcopy(patch)
        if isinstance(warping, np.ndarray):
            patch_copy.ctrlpts += warping
        position_data, det_jac, _ = IgaPostprocessing._extract_position(
            patch_copy,
            sample_size,
            None,
        )

        # Create point data
        point_data = {}
        for fieldname, fieldvalue in fields.items():
            if fieldvalue is None:
                continue
            fieldinterp: np.ndarray = np.atleast_2d(fieldvalue)
            IgaPostprocessing._update_data(
                fieldname,
                fieldinterp,
                sample_size,
                point_data,
            )
        if fields.get("det_jac") is None:
            IgaPostprocessing._update_data(
                "det_jac",
                np.atleast_2d(det_jac),
                sample_size,
                point_data,
            )

        if write_file:
            IgaPostprocessing._save_vtk(
                f"{folder}{filename}", position_data, point_data
            )

        return position_data, point_data

    @staticmethod
    def export_patch(
        patch: Any,  # iga-singlepatch object
        field_type: Literal["primal", "dual"],
        filename: str = "output",
        folder: Optional[str] = None,
        fields: dict = {},
        write_file: bool = True,
        warping: Optional[np.ndarray] = None,
        **extra_args,
    ):
        assert field_type in ["primal", "dual"]
        func = {
            "primal": IgaPostprocessing._export_primal_data,
            "dual": IgaPostprocessing._export_dual_data,
        }[field_type]
        return func(
            patch,
            filename=filename,
            folder=folder,
            fields=fields,
            write_file=write_file,
            warping=warping,
            **extra_args,
        )

    @staticmethod
    def plot_patch(
        patch_list: List[Any],
        filename: str = "plotpatch",
        folder: Optional[str] = None,
        primal_list: Optional[List[np.ndarray]] = None,
        add_determinant: bool = False,
        add_ctrlpts_net: bool = False,
        add_colorbar: bool = False,
        clim: Optional[tuple] = None,
        figsize: tuple = (5, 5),
    ):
        folder = Postprocessing.verify_folder_to_save(folder)
        # Verify input data
        if not isinstance(patch_list, List):
            patch_list = [patch_list]

        if any(patch.ndim != 2 for patch in patch_list):
            logger.warning("Plot cannot be created. Exit code.")
            return

        if primal_list is not None:
            if not isinstance(primal_list, list):
                logger.warning("Plot cannot be created. Exit code.")
                return
            if not all(
                isinstance(primal, np.ndarray) and primal.ndim == 1
                for primal in primal_list
            ):
                logger.warning("Plot cannot be created. Exit code.")
                return

        fig, ax = plt.subplots(figsize=figsize)
        ax.grid(None)
        ax.set_axis_off()

        # Preprocessing
        SAMPLESIZE = 101
        SAMPLEKNOTS = np.linspace(0, 1, SAMPLESIZE)
        X_list, Y_list, Z_list = [], [], []
        min_Z_list, max_Z_list = [], []
        for ii, patch in enumerate(patch_list):
            u_ctrlpts = None if primal_list is None else np.atleast_2d(primal_list[ii])
            u_interp, evalpts, det_jac, _ = patch.interpolate_field(
                knots_list=[SAMPLEKNOTS for _ in range(patch.ndim)],
                u_ctrlpts=u_ctrlpts,
                eval_geometry=True,
            )

            det_jac = det_jac.reshape((SAMPLESIZE, SAMPLESIZE)) / np.max(det_jac)
            if isinstance(u_interp, np.ndarray):
                u_interp = u_interp.reshape((SAMPLESIZE, SAMPLESIZE))
            else:
                u_interp = np.zeros_like(det_jac)
            X_list.append(evalpts[0].reshape((SAMPLESIZE, SAMPLESIZE)))
            Y_list.append(evalpts[1].reshape((SAMPLESIZE, SAMPLESIZE)))
            Z_list.append(det_jac if add_determinant else u_interp)
            min_Z_list.append(Z_list[-1].min())
            max_Z_list.append(Z_list[-1].max())

        if isinstance(clim, (tuple, list)):
            min_Z = clim[0]
            max_Z = clim[1]
        else:
            min_Z = np.floor(np.min(min_Z_list))
            max_Z = np.ceil(np.max(max_Z_list))

        # Add colormaps
        alpha_selector = lambda x: (
            0.25 if not (np.isclose(x, 0.0) or np.isclose(x, 1.0)) else 1.0
        )
        c = None
        for X, Y, Z, patch in zip(X_list, Y_list, Z_list, patch_list):
            c = ax.contourf(
                X,
                Y,
                Z,
                cmap="coolwarm",
                levels=np.linspace(min_Z, max_Z, 25),
            )

            # Add knots from knot-vector in physical space
            for knot in np.unique(patch.knotvector[0]):
                u_knots = patch.interpolate_field(
                    knots_list=[np.array([knot]), SAMPLEKNOTS]
                )[1]
                alpha = alpha_selector(knot)
                ax.plot(u_knots[0], u_knots[1], color="k", linestyle="-", alpha=alpha)

            for knot in np.unique(patch.knotvector[1]):
                u_knots = patch.interpolate_field(
                    knots_list=[SAMPLEKNOTS, np.array([knot])]
                )[1]
                alpha = alpha_selector(knot)
                ax.plot(u_knots[0], u_knots[1], color="k", linestyle="-", alpha=alpha)

            # Add control points
            if add_ctrlpts_net:
                ctrlptsX = np.reshape(patch.ctrlpts[0], patch.nbctrlpts, order="F")
                ctrlptsY = np.reshape(patch.ctrlpts[1], patch.nbctrlpts, order="F")
                ax.plot(
                    patch.ctrlpts[0],
                    patch.ctrlpts[1],
                    color="tab:orange",
                    marker="o",
                    linestyle="",
                    label="Control points",
                )
                for i in range(ctrlptsX.shape[0]):
                    ax.plot(
                        ctrlptsX[i],
                        ctrlptsY[i],
                        color="tab:orange",
                        linestyle="--",
                        alpha=0.3,
                    )

                for i in range(ctrlptsX.shape[1]):
                    ax.plot(
                        ctrlptsX[:, i],
                        ctrlptsY[:, i],
                        color="tab:orange",
                        linestyle="--",
                        alpha=0.3,
                    )

        if c is not None and add_determinant:
            fig.colorbar(c, ax=ax, shrink=0.7, label="Det. of Jacobien")

        if c is not None and add_colorbar:
            fig.colorbar(c, ax=ax, shrink=0.7, label="Scalar field")

        ax.autoscale()
        if add_ctrlpts_net:
            ax.legend(bbox_to_anchor=(0.5, 1.25), loc="upper center")
        ax.set_aspect("equal", adjustable="box")
        fig.tight_layout()
        fig.savefig(f"{folder}/{filename}.pdf", bbox_inches=None, pad_inches=0)
