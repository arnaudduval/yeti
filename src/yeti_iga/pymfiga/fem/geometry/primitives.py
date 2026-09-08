from typing import List, Union, Literal
try:
    from meshpy import triangle
except ImportError:
    triangle = None
import numpy as np


def make_arc_circle(
    center: tuple = (0, 0),
    radius: float = 0.2,
    start_theta: float = 0.0,
    end_theta: float = np.pi / 2,
    n_points: int = 10,
) -> List:
    x0, y0 = center
    return [
        (x0 + radius * np.cos(theta), y0 + radius * np.sin(theta))
        for theta in np.linspace(start_theta, end_theta, int(n_points))
    ]


def make_line(p1: Union[list, tuple], p2: Union[list, tuple], n_points: int) -> List:
    return [
        (p1[0] + i / n_points * (p2[0] - p1[0]), p1[1] + i / n_points * (p2[1] - p1[1]))
        for i in range(int(n_points) + 1)
    ]


class MeshpyGenerator:
    def __init__(
        self,
        filename: Literal[
            "triangle",
            "square",
            "square_modified",
            "plate_hole",
            "semicircle",
            "quarter_annulus",
        ],
        max_volume: float = 1e-2,
    ):
        assert isinstance(filename, str), "Insert name of geometry"
        self.filename = filename
        self._max_volume = max_volume
        return

    def mesh_geometry(self):
        geometry_map = {
            "triangle": (dict(a=1.0, b=1.0), self._create_triangle),
            "square": (dict(sqlen=1.0), self._create_square),
            "square_modified": (
                dict(sqlen=1.0, alpha=0.05),
                self._create_square_modified,
            ),
            "plate_hole": (dict(sqlen=1.0, radius=0.25), self._create_plate_hole),
            "semicircle": (
                dict(radius=0.2, offset=(0, 1.2), alpha=0.15, angle=3 * np.pi / 8),
                self._create_quart_circle,
            ),
            "quarter_annulus": (
                dict(Rint=0.2, Rext=1.0),
                self._create_quart_annulus,
            ),
        }

        assert self.filename in geometry_map, "Unknown geometry"
        default_extra_args, func = geometry_map[self.filename]
        extra_args = {**default_extra_args, "max_volume": self._max_volume}
        mesh = func(**extra_args)
        return mesh

    def _create_square(self, sqlen: float, max_volume: float):
        # Define geometry
        points = [(0, 0), (sqlen, 0), (sqlen, sqlen), (0, sqlen)]
        facets = [(0, 1), (1, 2), (2, 3), (3, 0)]

        # Create mesh
        meshinfo = triangle.MeshInfo()
        meshinfo.set_points(points)
        meshinfo.set_facets(facets)

        # Generate mesh
        mesh = triangle.build(meshinfo, max_volume=max_volume)
        return mesh

    def _create_triangle(self, a: float, b: float, max_volume: float):
        # Define geometry
        points = [(0, 0), (a, 0), (0, b)]
        facets = [(0, 1), (1, 2), (2, 0)]

        # Create mesh
        meshinfo = triangle.MeshInfo()
        meshinfo.set_points(points)
        meshinfo.set_facets(facets)

        # Generate mesh
        mesh = triangle.build(meshinfo, max_volume=max_volume)
        return mesh

    def _create_square_modified(self, sqlen: float, alpha: float, max_volume: float):

        # Define geometry
        points = [
            (0.0, 0.0),
            (sqlen, 0.0),
            (sqlen, sqlen),
            (alpha * sqlen, sqlen),
            *make_line(
                (alpha * sqlen, sqlen), (0.0, sqlen), np.ceil(1 / np.sqrt(max_volume))
            )[1:],
            *make_line(
                (0.0, sqlen),
                (0.0, sqlen * (1 - alpha)),
                np.ceil(1 / np.sqrt(max_volume)),
            )[1:],
        ]
        facets = [*[(i, i + 1) for i in range(len(points) - 1)], (len(points) - 1, 0)]

        # Create mesh
        meshinfo = triangle.MeshInfo()
        meshinfo.set_points(points)
        meshinfo.set_facets(facets)

        # Generate mesh
        mesh = triangle.build(meshinfo, max_volume=max_volume)
        return mesh

    def _create_quart_circle(
        self,
        radius: float,
        offset: tuple,
        alpha: float,
        angle: float,
        max_volume: float,
    ):

        # Define geometry
        arc_circle = make_arc_circle(
            offset,
            radius=radius,
            start_theta=-np.pi / 2,
            end_theta=-angle,
            n_points=np.ceil(np.sqrt(1 / max_volume)),
        )[:-1] + make_arc_circle(
            offset,
            radius=radius,
            start_theta=-angle,
            end_theta=0,
            n_points=np.ceil(np.sqrt(1 / max_volume)),
        )

        line = make_line(
            (offset[0], offset[1] - radius * (1 - alpha)),
            (offset[0], offset[1] - radius),
            n_points=np.ceil(np.sqrt(1 / max_volume)),
        )[1:-1]
        polygon = [offset, (offset[0], offset[1] - radius * (1 - alpha)), *line]

        # Create mesh
        info = triangle.MeshInfo()
        points = arc_circle + polygon
        info.set_points(points)
        facets = [*[(i, i + 1) for i in range(len(points) - 1)], (len(points) - 1, 0)]
        info.set_facets(facets)

        # Generate mesh
        mesh = triangle.build(info, max_volume=max_volume)
        return mesh

    def _create_quart_annulus(
        self,
        Rint: float,
        Rext: float,
        max_volume: float,
    ):

        # Define geometry
        nb_pts = np.sqrt(1 / max_volume)
        arc_circle_rext = make_arc_circle(
            (0, 0),
            radius=Rext,
            start_theta=0,
            end_theta=np.pi / 2,
            n_points=np.ceil(nb_pts / 2),
        )[:-1]

        arc_circle_rint = make_arc_circle(
            (0, 0),
            radius=Rint,
            start_theta=np.pi / 2,
            end_theta=0,
            n_points=np.ceil(2 * nb_pts),
        )[:-1]

        vertical_line = make_line(
            (0, Rext),
            (0, Rint),
            n_points=np.ceil(nb_pts),
        )[:-1]

        horizontal_line = make_line(
            (Rint, 0),
            (Rext, 0),
            n_points=np.ceil(nb_pts),
        )[:-1]

        # Create mesh
        info = triangle.MeshInfo()
        points = arc_circle_rext + vertical_line + arc_circle_rint + horizontal_line
        info.set_points(points)
        facets = [*[(i, i + 1) for i in range(len(points) - 1)], (len(points) - 1, 0)]
        info.set_facets(facets)

        # Generate mesh
        mesh = triangle.build(info, max_volume=max_volume)
        return mesh

    def _create_plate_hole(self, sqlen: float, radius: float, max_volume: float):

        # Define geometry
        arc_circle = make_arc_circle(
            (0, 0), radius=radius, n_points=int(np.ceil(1 / np.sqrt(max_volume)))
        )
        polygon = [(0, sqlen), (sqlen, sqlen), (sqlen, 0)]

        # Create mesh
        info = triangle.MeshInfo()
        points = arc_circle + polygon
        info.set_points(points)
        facets = [*[(i, i + 1) for i in range(len(points) - 1)], (len(points) - 1, 0)]
        info.set_facets(facets)

        # Generate mesh
        mesh = triangle.build(info, max_volume=max_volume)
        return mesh
