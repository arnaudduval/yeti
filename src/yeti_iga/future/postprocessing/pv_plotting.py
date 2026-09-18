"""
pyvista-based 3D visualization helpers for B-spline/NURBS patches.

Generalizes the ad-hoc VTK/pyvista code from 13_kirchhoff_love_shell.ipynb's
Part 5 into reusable functions, in the same spirit as plotting.py's 2D
matplotlib helpers: operate on a caller-supplied pyvista.Plotter (or a mesh),
return the objects added, so the caller stays in control of the scene.

Requires pyvista (``pip install -e ".[viz]"``), imported lazily by callers --
this module is kept separate from vtu.py (which has zero extra dependencies
beyond numpy) so that writing a .vtu file never requires pyvista to be
installed.
"""

from __future__ import annotations

import numpy as np
import pyvista as pv


def element_boundary_lines(mesh, points=None):
    """
    True element boundaries of a write_bezier_patch_vtu()-exported mesh, as a
    PolyData of closed 4-point line loops -- one per original B-spline/NURBS
    element.

    This is *not* the same as toggling "edges" on the mesh itself: VTK must
    tessellate each curved Bezier cell into flat triangles to render it, and
    that tessellation's edges (what "show edges" draws) are far denser than,
    and unrelated to, the actual element grid. This function reads the
    original cell connectivity instead, so it is exact and independent of
    tessellation resolution.

    Relies on write_bezier_patch_vtu()'s point ordering, where each cell's
    first 4 points are always its corners (vtu.py's ``_perm_2d`` docstring:
    "Corners (CCW): (0,0),(pu,0),(pu,pv),(0,pv)"), regardless of degree.

    Parameters
    ----------
    mesh : pyvista.UnstructuredGrid
        A mesh read from a write_bezier_patch_vtu(...) export (or any
        UnstructuredGrid following the same per-cell corner-first point
        ordering).
    points : array_like, shape (n_points, 3), optional
        Point positions to build the lines from, instead of ``mesh.points``
        -- e.g. ``mesh.warp_by_vector(...).points``, to draw the boundary in
        a deformed configuration while keeping the *original* (undeformed)
        cell connectivity. Must have the same length/ordering as
        ``mesh.points``.

    Returns
    -------
    pyvista.PolyData
        One closed 5-point line loop (4 corners + repeated first corner) per
        cell of ``mesh``.
    """
    if points is None:
        points = mesh.points

    lines = []
    for i in range(mesh.n_cells):
        pids = mesh.get_cell(i).point_ids
        a, b, c, d = pids[0], pids[1], pids[2], pids[3]
        lines.append([5, a, b, c, d, a])

    return pv.PolyData(points, lines=np.hstack(lines))


def add_deformed_shell(
    pl,
    mesh,
    field_name='displacement',
    factor=1.0,
    component=None,
    cmap='coolwarm',
    show_element_edges=True,
    edge_color='black',
    edge_line_width=2,
    scalar_bar_title=None,
    smooth_shading=True,
):
    """
    Add a shell's deformed shape, colored by a point field, to a pyvista
    Plotter -- warp, color, and (optionally) true element-boundary overlay
    in one call.

    Parameters
    ----------
    pl : pyvista.Plotter
        Plotter to add the mesh(es) to.
    mesh : pyvista.UnstructuredGrid
        A mesh read from a write_bezier_patch_vtu(...) export, with
        ``field_name`` as one of its point data arrays (as written by
        ``write_bezier_patch_vtu(..., field=..., field_name=...)``).
    field_name : str, default 'displacement'
        Point data array to warp by and to color by.
    factor : float, default 1.0
        Warp scale factor -- e.g. 100.0 to exaggerate millimeter-scale
        deflections on a meters-scale structure for visibility.
    component : int, optional
        Vector component to color by (e.g. 2 for the Z displacement). If
        None, colors by the field's magnitude.
    cmap : default 'coolwarm'
        Matplotlib colormap name.
    show_element_edges : bool, default True
        If True, overlay the true B-spline/NURBS element boundaries (via
        element_boundary_lines()) on the deformed shape. This is what you
        want instead of the mesh's own "show edges" toggle, which draws the
        much finer rendering tessellation instead (see element_boundary_lines'
        docstring).
    edge_color : default 'black'
        Element boundary line color.
    edge_line_width : float, default 2
        Element boundary line width.
    scalar_bar_title : str, optional
        Scalar bar title. Defaults to ``field_name`` (plus the component
        index, if given).
    smooth_shading : bool, default True
        Gouraud/Phong-interpolated shading across the surface.

    Returns
    -------
    warped : pyvista.PolyData
        The warped, colored surface (already added to `pl`).
    edges : pyvista.PolyData or None
        The element-boundary overlay (already added to `pl`), or None if
        show_element_edges is False.
    """
    warped = mesh.warp_by_vector(field_name, factor=factor)

    title = scalar_bar_title
    if title is None:
        title = field_name if component is None else f'{field_name}[{component}]'

    pl.add_mesh(warped, scalars=field_name, component=component, cmap=cmap,
                smooth_shading=smooth_shading, scalar_bar_args={'title': title})

    edges = None
    if show_element_edges:
        edges = element_boundary_lines(mesh, points=warped.points)
        pl.add_mesh(edges, color=edge_color, line_width=edge_line_width)

    return warped, edges
