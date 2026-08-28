from .helpers import (
    make_uniform_knotvector,
    make_BSPLINE_line,
    make_BSPLINE_quarter_circonference,
    create_NURBS_arc,
)
from geomdl import BSpline, NURBS, operations
from typing import Tuple, List, Union, Dict, Callable, Sequence, Literal
import numpy as np


class GeomdlGenerator:
    def __init__(
        self,
        filename: Literal[
            "line",
            "circle",
            "quarter_annulus",
            "nurbs_quarter_annulus",
            "nurbs_plate_hole",
            "square",
            "trapezium",
            "cube",
            "thick_ring",
            "rotated_quarter_annulus",
            "prism",
        ],
        geo_args: dict,
    ):
        assert isinstance(filename, str), "Insert name of geometry"
        self.filename = filename
        self._ndim: int = 0
        degree: Union[int, float, np.floating, np.ndarray] = geo_args.get("degree", 0)
        if np.isscalar(degree):
            degree = np.array([degree] * 3, dtype=int)
        self._degree = np.asarray(degree)
        nbel: Union[int, float, np.floating, np.ndarray] = geo_args.get("nbel", 0)
        if np.isscalar(nbel):
            nbel = np.array([nbel] * 3, dtype=int)
        self._nbel = np.asarray(nbel)
        self._extra_args: dict = geo_args.get("parameters", {})
        assert isinstance(
            self._extra_args, dict
        ), "Extra arguments should be dictionary"

    def export_geometry(self) -> Union[BSpline.Curve, BSpline.Surface, BSpline.Volume]:

        geometry_map: Dict[str, Tuple[int, dict, Callable]] = {
            "line": (1, {"L": 1.0}, self._create_line),
            "circle": (1, {"R": 1.0}, self._create_quarter_circle),
            "quarter_annulus": (
                2,
                {"Rin": 0.25, "Rex": 1.0},
                self._create_quarter_annulus,
            ),
            "nurbs_quarter_annulus": (
                2,
                {"Rin": 0.25, "Rex": 1.0},
                self._create_nurbs_quarter_annulus,
            ),
            "nurbs_plate_hole": (
                2,
                {"R": 1, "L": 4},
                self._create_nurbs_plate_with_hole,
            ),
            "square": (
                2,
                {"XY": np.array([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]])},
                self._create_quadrilateral,
            ),
            "trapezium": (
                2,
                {"XY": np.array([[0.0, -0.5], [0.4, -0.2], [0.4, 0.2], [0.0, 0.5]])},
                self._create_quadrilateral,
            ),
            "cube": (3, {"Lx": 1.0, "Ly": 1.0, "Lz": 1.0}, self._create_parallelepiped),
            "thick_ring": (
                3,
                {"Rin": 1.0, "Rex": 2.0, "height": 1.0},
                self._create_thick_ring,
            ),
            "rotated_quarter_annulus": (
                3,
                {"Rin": 1.0, "Rex": 2.0, "exc": 1.0},
                self._create_rotated_quarter_annulus,
            ),
            "prism": (
                3,
                {
                    "XY": np.array([[0.0, -7.5], [6.0, -2.5], [6.0, 2.5], [0.0, 7.5]]),
                    "height": 1.0,
                },
                self._create_prism,
            ),
        }

        if self.filename not in geometry_map:
            raise ValueError("Not developed in this library")

        ndim, default_extra_args, func = geometry_map[self.filename]
        self._ndim = ndim
        self._extra_args = default_extra_args | self._extra_args
        self._degree = self._degree[: self._ndim]
        self._nbel = self._nbel[: self._ndim]
        return func(*self._degree, *self._nbel, **self._extra_args)

    # ----------------
    # CREATE GEOMETRY
    # ----------------

    # 1D
    def _create_line(self, degree_u: int, nbel_u: int, **kwargs) -> BSpline.Curve:
        "Creates a line segment"

        length = kwargs.get("L")

        # Get uniform control points
        knotvector, ctrlpts = make_BSPLINE_line(degree_u, nbel_u)
        ctrlpts = length * ctrlpts

        # Create curve
        obj = BSpline.Curve()
        obj.degree = degree_u
        obj.ctrlpts = [[x, 0.0] for x in ctrlpts]
        obj.knotvector = knotvector
        return obj

    def _create_quarter_circle(
        self, degree_u: int, nbel_u: int, **kwargs
    ) -> BSpline.Curve:
        "Creates a quater of a circle"

        Radius = kwargs.get("R")

        # Construction of the arc
        knotvector, ctrlpts = make_BSPLINE_quarter_circonference(degree_u, nbel_u)
        ctrlpts = Radius * ctrlpts
        x_list = ctrlpts[:, 0]
        y_list = ctrlpts[:, 1]

        # Create curve
        obj = BSpline.Curve()
        obj.degree = degree_u
        obj.ctrlpts = [[x, y] for x, y in zip(x_list, y_list)]
        obj.knotvector = knotvector
        return obj

    # 2D

    def _create_nurbs_plate_with_hole(
        self, degree_u: int, degree_v: int, nbel_u: int, nbel_v: int, **kwargs
    ) -> NURBS.Surface:

        radius, length = kwargs.get("R", 0.0), kwargs.get("L", 0.0)
        assert degree_u > 1 and degree_v > 0 and nbel_u % 2 == 0
        assert radius < length

        def create_l_shape(degree):
            ones = np.ones(degree)
            ctrlpts_1d = np.hstack((ones, make_BSPLINE_line(degree, 1)[-1][::-1]))
            ctrlpts_3d = [
                [-x, y, 0.0, 1.0] for x, y in zip(ctrlpts_1d, ctrlpts_1d[::-1])
            ]
            knotvector = make_uniform_knotvector(degree, 2, multiplicity=degree)
            obj = NURBS.Curve()
            obj.degree = degree
            obj.ctrlptsw = ctrlpts_3d
            obj.knotvector = knotvector
            return obj

        def make_nurbs_contour(degree: int) -> Sequence[np.ndarray]:

            # Create line
            obj_line = create_l_shape(degree)

            # Create circle
            obj_circle = create_NURBS_arc(degree, 1, np.pi, np.pi / 2)
            operations.insert_knot(obj_circle, [0.5], [degree])

            # Get output
            info = []
            for obj in [obj_line, obj_circle]:
                info.append(obj.ctrlpts)
                info.append(obj.weights)
            info.append(obj_circle.knotvector)
            return info

        def interpolate(p1: List[list], p2: List[list], where: np.ndarray):
            assert len(p1) == len(p2)
            return [
                [(p1[i] * (1 - w) + p2[i] * w) for i in range(len(p1))] for w in where
            ]

        # Construction of the arc
        (
            ctrlpts_l_shape,
            _,
            ctrlpts_arc,
            weights_arc,
            knotvector_u,
        ) = make_nurbs_contour(
            degree_u
        )  # 2 elements

        # Construction of line
        knotvector_v, greville_hline = make_BSPLINE_line(degree_v, 1)  # 1 element

        ctrlpts = []
        for ii in range(len(ctrlpts_l_shape)):
            # First point
            w_arc = weights_arc[ii]
            c_arc = ctrlpts_arc[ii]
            p1 = [radius * c_arc[0] * w_arc, radius * c_arc[1] * w_arc, 0, w_arc]

            # Second point
            c_vline = ctrlpts_l_shape[ii]
            p2 = [length * c_vline[0], length * c_vline[1], 0, 1]

            # Add new point
            newpoint = interpolate(p1, p2, greville_hline)
            ctrlpts.extend(newpoint)

        # Create surface
        obj = NURBS.Surface()
        obj.degree_u = degree_u
        obj.degree_v = degree_v
        obj.set_ctrlpts(ctrlpts, len(ctrlpts_arc), len(greville_hline))
        obj.knotvector_u = knotvector_u
        obj.knotvector_v = knotvector_v

        # Add knot refinement
        for k, (nbel, is_forward) in enumerate(zip([nbel_u, nbel_v], [True, False])):
            for knot in np.linspace(0.0, 1.0, nbel + 1)[1:-1]:
                kpdir = [knot, 0] if is_forward else [0, knot]
                opdir = [1, 0] if is_forward else [0, 1]
                if k == 0 and np.isclose(knot, 0.5):
                    continue
                operations.insert_knot(obj, kpdir, opdir)

        return obj

    def _create_nurbs_quarter_annulus(
        self, degree_u: int, degree_v: int, nbel_u: int, nbel_v: int, **kwargs
    ) -> NURBS.Surface:
        "Creates a quarter of a ring (or annulus)"

        Rin, Rex = kwargs.get("Rin", 0.0), kwargs.get("Rex", 0.0)

        # Construction of the arc
        obj_arc = create_NURBS_arc(degree_v, nbel_v, np.pi, np.pi / 2)
        knotvector_v = obj_arc.knotvector
        ctrlpts_arc = obj_arc.ctrlpts
        weights_arc = obj_arc.weights

        # Construction of line
        knotvector_u, ctrlpts_line = make_BSPLINE_line(degree_u, nbel_u)
        ctrlpts_line = Rin + ctrlpts_line * (Rex - Rin)

        # Construction of annulus sector
        ctrlpts = [
            [x_line * x_arc * w_arc, x_line * y_arc * w_arc, 0.0, w_arc]
            for x_line in ctrlpts_line
            for (x_arc, y_arc, _), w_arc in zip(ctrlpts_arc, weights_arc)
        ]

        # Create surface
        obj = NURBS.Surface()
        obj.degree_u = degree_u
        obj.degree_v = degree_v
        obj.set_ctrlpts(ctrlpts, len(ctrlpts_line), len(ctrlpts_arc))
        obj.knotvector_u = knotvector_u
        obj.knotvector_v = knotvector_v

        return obj

    def _create_quarter_annulus(
        self, degree_u: int, degree_v: int, nbel_u: int, nbel_v: int, **kwargs
    ) -> BSpline.Surface:
        "Creates a quarter of a ring (or annulus)"

        Rin, Rex = kwargs.get("Rin", 0.0), kwargs.get("Rex", 0.0)

        # Construction of the arc
        nb_ctrlpts_v = degree_v + nbel_v
        knotvector_v, ctrlpts_arc = make_BSPLINE_quarter_circonference(degree_v, nbel_v)

        # Construction of line
        nb_ctrlpts_u = degree_u + nbel_u
        knotvector_u, ctrlpts_line = make_BSPLINE_line(degree_u, nbel_u)
        ctrlpts_line = Rin + ctrlpts_line * (Rex - Rin)

        # Construction of annulus sector
        ctrlpts = [
            [x_line * x_arc, x_line * y_arc, 0.0]
            for x_line in ctrlpts_line
            for x_arc, y_arc in ctrlpts_arc
        ]

        # Create surface
        obj = BSpline.Surface()
        obj.degree_u = degree_u
        obj.degree_v = degree_v
        obj.ctrlpts_size_u, obj.ctrlpts_size_v = int(nb_ctrlpts_u), int(nb_ctrlpts_v)
        obj.set_ctrlpts(ctrlpts, nb_ctrlpts_u, nb_ctrlpts_v)
        obj.knotvector_u = knotvector_u
        obj.knotvector_v = knotvector_v

        return obj

    def _create_quadrilateral(
        self, degree_u: int, degree_v: int, nbel_u: int, nbel_v: int, **kwargs
    ) -> BSpline.Surface:
        "Creates a quadrilateral given coordinates in counterclockwise direction"

        xy = kwargs.get("XY", np.array([[]]))

        # Set reference position and real position
        x0 = [0.0, 1.0, 1.0, 0.0]
        y0 = [0.0, 0.0, 1.0, 1.0]
        x1 = xy[:, 0]
        y1 = xy[:, 1]

        # Transformation of control points
        # x1 = ax1 x0 + ax2 y0 + ax3 x0 y0 + ax4
        # y1 = ay1 x0 + ay2 y0 + ay3 x0 y0 + ay4
        T = [[x0[i], y0[i], x0[i] * y0[i], 1] for i in range(4)]
        ax = np.linalg.solve(T, x1)
        ay = np.linalg.solve(T, y1)

        # Set control points
        nb_ctrlpts_u = degree_u + nbel_u
        nb_ctrlpts_v = degree_v + nbel_v
        knotvector_u, ctrlpts_u = make_BSPLINE_line(degree_u, nbel_u)
        knotvector_v, ctrlpts_v = make_BSPLINE_line(degree_v, nbel_v)

        ctrlpts = [
            [
                ax[0] * xt + ax[1] * yt + ax[2] * xt * yt + ax[3],
                ay[0] * xt + ay[1] * yt + ay[2] * xt * yt + ay[3],
                0.0,
            ]
            for xt in ctrlpts_u
            for yt in ctrlpts_v
        ]

        # Create surface
        obj = BSpline.Surface()
        obj.degree_u = degree_u
        obj.degree_v = degree_v
        obj.ctrlpts_size_u, obj.ctrlpts_size_v = int(nb_ctrlpts_u), int(nb_ctrlpts_v)
        obj.set_ctrlpts(ctrlpts, nb_ctrlpts_u, nb_ctrlpts_v)
        obj.knotvector_u = knotvector_u
        obj.knotvector_v = knotvector_v

        return obj

    # 3D
    def _create_parallelepiped(
        self,
        degree_u: int,
        degree_v: int,
        degree_w: int,
        nbel_u: int,
        nbel_v: int,
        nbel_w: int,
        **kwargs,
    ) -> BSpline.Volume:
        "Creates a brick (or parallelepiped)"

        Lx, Ly, Lz = kwargs.get("Lx"), kwargs.get("Ly"), kwargs.get("Lz")

        # Set number of control points
        nb_ctrlpts_u = degree_u + nbel_u
        nb_ctrlpts_v = degree_v + nbel_v
        nb_ctrlpts_w = degree_w + nbel_w

        # Get uniform control points
        knotvector_u, ctrlpts_u = make_BSPLINE_line(degree_u, nbel_u)
        knotvector_v, ctrlpts_v = make_BSPLINE_line(degree_v, nbel_v)
        knotvector_w, ctrlpts_w = make_BSPLINE_line(degree_w, nbel_w)

        # Create control points of the volume
        ctrlpts = [
            [cptu * Lx, cptv * Ly, cptw * Lz]
            for cptw in ctrlpts_w
            for cptu in ctrlpts_u
            for cptv in ctrlpts_v
        ]

        # Create a B-spline volume
        obj = BSpline.Volume()
        obj.degree_u, obj.degree_v, obj.degree_w = degree_u, degree_v, degree_w
        obj.ctrlpts_size_u, obj.ctrlpts_size_v, obj.ctrlpts_size_w = (
            int(nb_ctrlpts_u),
            int(nb_ctrlpts_v),
            int(nb_ctrlpts_w),
        )
        obj.set_ctrlpts(ctrlpts, nb_ctrlpts_u, nb_ctrlpts_v, nb_ctrlpts_w)
        obj.knotvector_u = knotvector_u
        obj.knotvector_v = knotvector_v
        obj.knotvector_w = knotvector_w

        return obj

    def _create_thick_ring(
        self,
        degree_u: int,
        degree_v: int,
        degree_w: int,
        nbel_u: int,
        nbel_v: int,
        nbel_w: int,
        **kwargs,
    ) -> BSpline.Volume:
        "Creates a thick ring (quarter of annulus extruded)"

        Rin, Rex, height = (
            kwargs.get("Rin", 0.0),
            kwargs.get("Rex", 0.0),
            kwargs.get("height", 0.0),
        )

        # construction of the arc
        nb_ctrlpts_v = degree_v + nbel_v
        knotvector_v, ctrlpts_arc = make_BSPLINE_quarter_circonference(degree_v, nbel_v)

        # construction of line
        nb_ctrlpts_u = degree_u + nbel_u
        knotvector_u, ctrlpts_line = make_BSPLINE_line(degree_u, nbel_u)
        ctrlpts_line = Rin + ctrlpts_line * (Rex - Rin)

        nb_ctrlpts_w = degree_w + nbel_w
        knotvector_w, ctrlpts_height = make_BSPLINE_line(degree_w, nbel_w)
        ctrlpts_height = height * ctrlpts_height

        # construction of annulus sector
        ctrlpts = [
            [x_line * x_arc, x_line * y_arc, z]
            for z in ctrlpts_height
            for x_line in ctrlpts_line
            for x_arc, y_arc in ctrlpts_arc
        ]

        # Create volume
        obj = BSpline.Volume()
        obj.degree_u, obj.degree_v, obj.degree_w = degree_u, degree_v, degree_w
        obj.ctrlpts_size_u, obj.ctrlpts_size_v, obj.ctrlpts_size_w = (
            int(nb_ctrlpts_u),
            int(nb_ctrlpts_v),
            int(nb_ctrlpts_w),
        )
        obj.set_ctrlpts(ctrlpts, nb_ctrlpts_u, nb_ctrlpts_v, nb_ctrlpts_w)
        obj.knotvector_u = knotvector_u
        obj.knotvector_v = knotvector_v
        obj.knotvector_w = knotvector_w

        return obj

    def _create_rotated_quarter_annulus(
        self,
        degree_u: int,
        degree_v: int,
        degree_w: int,
        nbel_u: int,
        nbel_v: int,
        nbel_w: int,
        **kwargs,
    ) -> BSpline.Volume:
        "Creates a quarter of a ring rotated (or revolted)"

        Rin, Rex, exc = (
            kwargs.get("Rin", 0.0),
            kwargs.get("Rex", 0.0),
            kwargs.get("exc", 0.0),
        )

        # construction of the arc 1
        nb_ctrlpts_v = degree_v + nbel_v
        knotvector_v, ctrlpts_arc_1 = make_BSPLINE_quarter_circonference(
            degree_v, nbel_v
        )

        # construction of line
        nb_ctrlpts_u = degree_u + nbel_u
        knotvector_u, ctrlpts_line = make_BSPLINE_line(degree_u, nbel_u)
        ctrlpts_line = Rin + ctrlpts_line * (Rex - Rin)

        # construction of the arc 2
        nb_ctrlpts_w = degree_w + nbel_w
        knotvector_w, ctrlpts_arc_2 = make_BSPLINE_quarter_circonference(
            degree_w, nbel_w
        )

        # Get control points
        ctrlpts = [
            [
                x_line * x_arc_1,
                (x_line * y_arc_1 + exc) * y_arc_2,
                (x_line * y_arc_1 + exc) * z_arc_2,
            ]
            for y_arc_2, z_arc_2 in ctrlpts_arc_2
            for x_line in ctrlpts_line
            for x_arc_1, y_arc_1 in ctrlpts_arc_1
        ]

        # Create volume
        obj = BSpline.Volume()
        obj.degree_u, obj.degree_v, obj.degree_w = degree_u, degree_v, degree_w
        obj.ctrlpts_size_u, obj.ctrlpts_size_v, obj.ctrlpts_size_w = (
            int(nb_ctrlpts_u),
            int(nb_ctrlpts_v),
            int(nb_ctrlpts_w),
        )
        obj.set_ctrlpts(ctrlpts, nb_ctrlpts_u, nb_ctrlpts_v, nb_ctrlpts_w)
        obj.knotvector_u = knotvector_u
        obj.knotvector_v = knotvector_v
        obj.knotvector_w = knotvector_w

        return obj

    def _create_prism(
        self,
        degree_u: int,
        degree_v: int,
        degree_w: int,
        nbel_u: int,
        nbel_v: int,
        nbel_w: int,
        **kwargs,
    ) -> BSpline.Volume:
        """Creates a prism using a quadrilateral as a base.
        The quadrilateral coordinates are given in counterclockwise direction"""

        xy, height = kwargs.get("XY", np.array([[]])), kwargs.get("height", 0.0)

        # Set reference position and real position
        x0 = [0.0, 1.0, 1.0, 0.0]
        y0 = [0.0, 0.0, 1.0, 1.0]
        x1 = xy[:, 0]
        y1 = xy[:, 1]

        # Transformation of control points
        # x1 = ax1 x0 + ax2 y0 + ax3 x0 y0 + ax4
        # y1 = ay1 x0 + ay2 y0 + ay3 x0 y0 + ay4
        T = [[x0[i], y0[i], x0[i] * y0[i], 1] for i in range(4)]
        ax = np.linalg.solve(T, x1)
        ay = np.linalg.solve(T, y1)

        # Set control points
        nb_ctrlpts_u = degree_u + nbel_u
        nb_ctrlpts_v = degree_v + nbel_v
        nb_ctrlpts_w = degree_w + nbel_w
        knotvector_u, ctrlpts_u = make_BSPLINE_line(degree_u, nbel_u)
        knotvector_v, ctrlpts_v = make_BSPLINE_line(degree_v, nbel_v)
        knotvector_w, ctrlpts_w = make_BSPLINE_line(degree_w, nbel_w)

        ctrlpts = [
            [
                ax[0] * xt + ax[1] * yt + ax[2] * xt * yt + ax[3],
                ay[0] * xt + ay[1] * yt + ay[2] * xt * yt + ay[3],
                zt * height,
            ]
            for zt in ctrlpts_w
            for xt in ctrlpts_u
            for yt in ctrlpts_v
        ]

        # Create volume
        obj = BSpline.Volume()
        obj.degree_u, obj.degree_v, obj.degree_w = degree_u, degree_v, degree_w
        obj.ctrlpts_size_u, obj.ctrlpts_size_v, obj.ctrlpts_size_w = (
            int(nb_ctrlpts_u),
            int(nb_ctrlpts_v),
            int(nb_ctrlpts_w),
        )
        obj.set_ctrlpts(ctrlpts, nb_ctrlpts_u, nb_ctrlpts_v, nb_ctrlpts_w)
        obj.knotvector_u = knotvector_u
        obj.knotvector_v = knotvector_v
        obj.knotvector_w = knotvector_w

        return obj
