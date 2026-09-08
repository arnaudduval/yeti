from yeti_iga.pymfiga.common.base.cls import BasePatch
from yeti_iga.pymfiga.common.numerics.operations import BsplineOperations, NurbsOperations
from typing import Union, List, Tuple
from geomdl import BSpline, NURBS
import numpy as np


class BasePatchIGA(BasePatch):

    _ndim_default = 3

    def __init__(
        self,
        obj: Union[
            BSpline.Curve,
            BSpline.Surface,
            BSpline.Volume,
            NURBS.Curve,
            NURBS.Surface,
            NURBS.Volume,
        ],
    ):
        # Extract information
        self.ndim: int = self._read_dimensionality(obj)
        self.degree: np.ndarray = self._read_degree(obj)
        self.knotvector: List[np.ndarray] = self._read_knotvector(obj)
        self.nbctrlpts: np.ndarray = self._compute_nbctrlpts()
        self.ctrlpts, self.nurbs_weights = self._read_control_points(obj)
        self._box_corners = None

    @property
    def nbctrlpts_total(self) -> int:
        "Total number of control points (product over all parametric directions)."
        return int(np.prod(self.nbctrlpts[: self.ndim]))

    @property
    def operator_engine(self):
        return NurbsOperations if self.nurbs_weights.size > 0 else BsplineOperations

    def _read_dimensionality(self, obj) -> int:
        is_curve = (
            hasattr(obj, "degree")
            and hasattr(obj, "knotvector")
            and not hasattr(obj, "knotvector_v")
        )

        is_surface = (
            hasattr(obj, "degree_u")
            and hasattr(obj, "degree_v")
            and not hasattr(obj, "degree_w")
        )

        is_volume = (
            hasattr(obj, "degree_u")
            and hasattr(obj, "degree_v")
            and hasattr(obj, "degree_w")
        )

        if is_volume:
            ndim = 3
        elif is_surface:
            ndim = 2
        elif is_curve:
            ndim = 1
        else:
            raise TypeError("Geometry is not a supported geomdl Curve/Surface/Volume.")
        return ndim

    def _read_degree(self, obj) -> np.ndarray:
        degree = np.ones(self._ndim_default, dtype=int)
        if self.ndim == 1:
            degree[0] = int(obj.degree)
        else:
            degree[0] = int(obj.degree_u)
            degree[1] = int(obj.degree_v)
            if self.ndim == 3:
                degree[2] = int(obj.degree_w)
        return degree

    def _read_knotvector(self, obj) -> List[np.ndarray]:
        kv: List[np.ndarray] = []
        if self.ndim == 1:
            kv.append(np.asarray(obj.knotvector, dtype=float))
        else:
            kv.append(np.asarray(obj.knotvector_u, dtype=float))
            kv.append(np.asarray(obj.knotvector_v, dtype=float))
            if self.ndim == 3:
                kv.append(np.asarray(obj.knotvector_w, dtype=float))
        return kv

    def _compute_nbctrlpts(self) -> np.ndarray:
        nb = np.ones(self._ndim_default, dtype=int)
        for i in range(self.ndim):
            nb[i] = len(self.knotvector[i]) - self.degree[i] - 1
        return nb

    def _read_control_points(self, obj) -> Tuple[np.ndarray, np.ndarray]:
        nbctrlpts = self.nbctrlpts
        ctrlpts = np.array(obj.ctrlpts)
        weights = np.array(obj.weights or np.array([]))

        # NOTE: by default Geomdl object has 3D coordinates
        ctrlpts = np.reshape(
            ctrlpts, (nbctrlpts[0], nbctrlpts[1], nbctrlpts[2], -1)
        )  # i, j, k, l
        ctrlpts = np.moveaxis(ctrlpts, 1, 0)  # j, i, k, l
        ctrlpts = np.moveaxis(ctrlpts, -1, 0)  # l, j, i, k

        # Get the correct dimension
        ndim = ctrlpts.shape[0]
        ctrlpts = np.reshape(ctrlpts, (ndim, -1))
        new_ctrlpts = np.zeros((self._ndim_default, self.nbctrlpts_total))
        new_ctrlpts[:ndim] = ctrlpts

        if weights.size > 0:
            weights = np.reshape(weights, (nbctrlpts[0], nbctrlpts[1], nbctrlpts[2]))
            weights = np.moveaxis(weights, 1, 0)
            weights = np.ravel(weights)

        return new_ctrlpts, weights

    @property
    def box_corners(self) -> np.ndarray:
        "Return the coordinates of the corners of the bounding box of the patch."
        if self._box_corners is not None:
            return self._box_corners

        min_coords = np.min(self.ctrlpts[: self.ndim, :], axis=1)
        max_coords = np.max(self.ctrlpts[: self.ndim, :], axis=1)
        if self.ndim == 1:
            corners = np.array([[min_coords[0]], [max_coords[0]]])
        elif self.ndim == 2:
            corners = np.array(
                [
                    [min_coords[0], min_coords[1]],
                    [max_coords[0], min_coords[1]],
                    [min_coords[0], max_coords[1]],
                    [max_coords[0], max_coords[1]],
                ]
            ).T
        elif self.ndim == 3:
            corners = np.array(
                [
                    [min_coords[0], min_coords[1], min_coords[2]],
                    [max_coords[0], min_coords[1], min_coords[2]],
                    [min_coords[0], max_coords[1], min_coords[2]],
                    [max_coords[0], max_coords[1], min_coords[2]],
                    [min_coords[0], min_coords[1], max_coords[2]],
                    [max_coords[0], min_coords[1], max_coords[2]],
                    [min_coords[0], max_coords[1], max_coords[2]],
                    [max_coords[0], max_coords[1], max_coords[2]],
                ]
            ).T
        else:
            raise ValueError("Invalid number of dimensions.")
        self._box_corners = corners
        return self._box_corners

    def characteristic_length(self) -> float:
        """
        Return a characteristic length of the patch,
        defined as the maximum distance in the box corners.
        """
        corners = self.box_corners
        max_distance = 0.0
        for i in range(corners.shape[1]):
            for j in range(i + 1, corners.shape[1]):
                distance = np.linalg.norm(corners[:, i] - corners[:, j])
                if distance > max_distance:
                    max_distance = float(distance)
        return max_distance
