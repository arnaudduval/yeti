from .postprocessing import Postprocessing
from typing import Dict, Optional, Literal, Any
import numpy as np
import logging

logger = logging.getLogger("SRC.FEA.POSTPROCESSING")

try:
    import meshio
except:
    logger.debug(
        "meshio module is not installed. Some functions in post-processing may be disabled"
    )


class FeaPostprocessing:
    @staticmethod
    def export_meshpy(
        mesh: Any,  # meshpy object
        filename: str = "meshpy",
        folder: Optional[str] = None,
        point_data: Dict = {},
        cell_data: Dict = {},
        warping: Optional[np.ndarray] = None,
        file_format: Literal["vtk", "vtu", "xdmf"] = "vtk",
        mesh_order: Optional[Literal["linear", "quadratic"]] = None,
    ):
        folder = Postprocessing.verify_folder_to_save(folder)

        if mesh_order is None:
            mesh_order = mesh.mesh_order
        if mesh_order == "linear":
            cell_type = "triangle"
        elif mesh_order == "quadratic":
            cell_type = "triangle6"
        else:
            raise NotImplementedError()

        # Create mesh to export
        order = 1 if mesh_order == "linear" else 2
        elements = mesh.recover_elements("lagrange", order)
        points = mesh.recover_points("lagrange", order)
        modpoints = (points + warping) if warping is not None else points
        cells = [(cell_type, elements)]
        mesh2export = meshio.Mesh(points=modpoints, cells=cells)

        # Add point data to the mesh if provided
        if point_data:
            for key, val in point_data.items():
                mesh2export.point_data[key] = val

        # Add cell data to the mesh if provided
        if cell_data:
            for key, val in cell_data.items():
                mesh2export.cell_data[key] = val

        # Save the mesh to a file
        supported_formats = ["vtk", "vtu", "xdmf"]
        if file_format not in supported_formats:
            text_1 = f"Unsupported file format: {str(file_format)}."
            text_2 = f"Supported formats are: {supported_formats}"
            logger.debug(f"{text_1}\n{text_2}")
            raise NotImplementedError()

        meshio.write(
            f"{folder}/{filename}.{file_format}", mesh2export, file_format=file_format
        )

    @staticmethod
    def plot_meshpy(
        mesh: Any,
        filename: str = "meshpy",
        folder: Optional[str] = None,
        mesh_order: Optional[Literal["linear", "quadratic"]] = None,
    ):
        from matplotlib import pyplot as plt
        from matplotlib.collections import PolyCollection
        from collections import Counter

        if mesh_order is None:
            mesh_order = mesh.mesh_order
            assert mesh_order == "linear"

        fig, ax = plt.subplots(figsize=(5, 5))
        ax.grid(None)
        ax.set_axis_off()

        order = 1 if mesh_order == "linear" else 2
        elements = mesh.recover_elements("lagrange", order)
        points = np.array(mesh.recover_points("lagrange", order))
        polys = [points[elem] for elem in elements]
        cmap = plt.get_cmap("GnBu")
        pc = PolyCollection(polys, facecolors=cmap(0), edgecolor=(0, 0, 0, 0.25))
        ax.add_collection(pc)
        edges = []
        for tri in elements:
            edges.extend(
                [
                    tuple(sorted((tri[0], tri[1]))),
                    tuple(sorted((tri[1], tri[2]))),
                    tuple(sorted((tri[2], tri[0]))),
                ]
            )
        counts = Counter(edges)
        boundary_edges = [edge for edge, c in counts.items() if c == 1]
        for i, j in boundary_edges:
            x = [points[i][0], points[j][0]]
            y = [points[i][1], points[j][1]]
            ax.plot(x, y, color="k")

        ax.autoscale()
        ax.set_aspect("equal", adjustable="box")
        fig.tight_layout()
        folder = Postprocessing.verify_folder_to_save(folder)
        fig.savefig(f"{folder}/{filename}.pdf", bbox_inches="tight", pad_inches=0)
