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
    Straight-line element boundaries of a write_bezier_patch_vtu()-exported
    mesh, as a PolyData of closed 4-point line loops -- one per original
    B-spline/NURBS element, connecting each element's 4 corners directly.

    This is *not* the same as toggling "edges" on the mesh itself: VTK must
    tessellate each curved Bezier cell into flat triangles to render it, and
    that tessellation's edges (what "show edges" draws) are far denser than,
    and unrelated to, the actual element grid. This function reads the
    original cell connectivity instead, giving the exact element grid rather
    than the rendering tessellation.

    Straight corner-to-corner segments are *exact* whenever an element's
    edges genuinely are straight lines (e.g. a flat plate with control points
    on a regular grid -- true by the affine/linear precision of the B-spline
    basis, regardless of degree). For an element with a genuinely curved edge
    (a NURBS arc, a curved shell), this draws a straight-line approximation
    instead of the true curve -- use element_boundary_curves() for those.

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


def _iso_line_walk(patch, n_samples):
    """
    Yield (spans, params) for each knot-break line (element edge) of a patch,
    in each parametric direction -- the shared iso-parametric-line walk
    behind element_boundary_curves() and element_boundary_curves_deformed(),
    mirroring plotting.py's plot_patches_2d().
    """
    su = patch.tensor.components[0]
    sv = patch.tensor.components[1]
    u_breaks = np.unique(su.knot_vector)
    v_breaks = np.unique(sv.knot_vector)

    for u_val in u_breaks:
        vs = np.linspace(v_breaks[0], v_breaks[-1], n_samples)
        spans = np.array(
            [[su.find_span(u_val), sv.find_span(v)] for v in vs], dtype=np.int32)
        params = np.column_stack([np.full(n_samples, u_val), vs])
        yield spans, params

    for v_val in v_breaks:
        us = np.linspace(u_breaks[0], u_breaks[-1], n_samples)
        spans = np.array(
            [[su.find_span(u), sv.find_span(v_val)] for u in us], dtype=np.int32)
        params = np.column_stack([us, np.full(n_samples, v_val)])
        yield spans, params


class _PolylineBuilder:
    """Accumulates polylines (each a list of points) into one pyvista.PolyData."""

    def __init__(self):
        self._points = []
        self._lines = []
        self._offset = 0

    def add(self, pts):
        if pts.shape[1] < 3:  # pv.PolyData needs 3D points even for a 2D patch
            pts = np.column_stack([pts, np.zeros((len(pts), 3 - pts.shape[1]))])
        n = len(pts)
        self._points.append(pts)
        self._lines.append(np.hstack([[n], np.arange(self._offset, self._offset + n)]))
        self._offset += n

    def build(self):
        return pv.PolyData(np.vstack(self._points), lines=np.hstack(self._lines))


def control_net_lines(patch, points=None):
    """
    Control net (control polygon) of a Patch: straight segments connecting
    each control point to its immediate neighbors in u and v, as a PolyData.

    This is the control *cage*, not the surface -- always drawn as straight
    lines between control points, by definition, whether or not the surface
    itself is curved (unlike element_boundary_curves(), which traces the
    true curved surface boundary).

    Parameters
    ----------
    patch : Patch
        2D patch (any physical dimension).
    points : array_like, shape (n_cp, 3), optional
        Control point positions to use instead of each
        ``patch.control_point(i)`` -- e.g. to show a *deformed* control net
        by adding a (scaled) displacement to each control point. Must be in
        the same order (u-fastest) as ``patch.control_point(i)``.

    Returns
    -------
    pyvista.PolyData
    """
    nu, nv = patch.local_shape
    if points is None:
        points = np.array([patch.control_point(i) for i in range(patch.n_cp)])
    points = np.asarray(points)
    if points.shape[1] < 3:  # pv.PolyData needs 3D points even for a 2D patch
        points = np.column_stack([points, np.zeros((len(points), 3 - points.shape[1]))])

    def idx(iu, iv):
        return iu + nu * iv

    lines = []
    for iv in range(nv):
        for iu in range(nu - 1):
            lines.append([2, idx(iu, iv), idx(iu + 1, iv)])
    for iu in range(nu):
        for iv in range(nv - 1):
            lines.append([2, idx(iu, iv), idx(iu, iv + 1)])

    return pv.PolyData(points, lines=np.hstack(lines))


def element_boundary_curves(patch, n_samples=20):
    """
    True (curved) element boundaries of a Patch, evaluated directly from its
    B-spline/NURBS basis -- one polyline per knot-break line in each
    parametric direction (i.e. every element edge, including interior ones
    for a multi-element patch), sampled at ``n_samples`` points each.

    Unlike element_boundary_lines() (straight corner-to-corner segments, read
    from an already-exported/tessellated file), this calls
    Patch.evaluate_patch_nd_omp() directly -- the same NURBS-aware evaluator
    used everywhere else in ``future`` -- so a genuinely curved element edge
    (a NURBS arc, a curved shell) comes out curved here, not as a chord
    approximation. Mirrors plotting.py's plot_patches_2d() iso-parametric-line
    walk, generalized to any physical dimension (not just 2D) and returning a
    pyvista.PolyData instead of drawing on a matplotlib Axes.

    This traces the *reference* (undeformed) geometry. For the boundary of a
    deformed FE solution, use element_boundary_curves_deformed() instead --
    a straight chord between two warped corner points is only exact if the
    displacement field happens to vary affinely along that edge, which is not
    guaranteed in general.

    Parameters
    ----------
    patch : Patch
        2D patch (any physical dimension) to trace.
    n_samples : int, default 20
        Number of points sampled along each element edge -- raise this for a
        very sharply curved edge if the polyline still looks faceted.

    Returns
    -------
    pyvista.PolyData
        One polyline per knot-break line (per direction), each with
        n_samples points.
    """
    builder = _PolylineBuilder()
    for spans, params in _iso_line_walk(patch, n_samples):
        builder.add(patch.evaluate_patch_nd_omp(spans, params))
    return builder.build()


def element_boundary_curves_deformed(patch, u_global, factor=1.0, n_samples=20):
    """
    Like element_boundary_curves(), but for a *deformed* FE solution: adds
    the solution's displacement (scaled by ``factor``) to the reference
    geometry at the same fine parametric samples along each element edge,
    via PatchEvaluator.

    A straight chord between two warped corner points (as
    element_boundary_lines() draws) is only exact when the displacement
    field varies affinely along that edge -- not guaranteed for a general FE
    solution, and visibly wrong wherever it doesn't (e.g. near a sharp local
    feature such as a point support). This evaluates the true deformed edge
    at ``n_samples`` points instead, exactly like the geometry-only case.

    Parameters
    ----------
    patch : Patch
        Patch with a PatchDOFManager attached (required by PatchEvaluator).
    u_global : ndarray, shape (n_dof,)
        Global solution vector matching `patch`'s dof numbering (e.g. from
        solving a PatchIntegrator.integrate_shell_stiffness() system). Only
        the first 3 dof components per control point (translations) are
        used.
    factor : float, default 1.0
        Displacement scale factor -- pass the same value given to
        warp_by_vector()/add_deformed_shell() so this overlay matches the
        colored surface.
    n_samples : int, default 20
        Number of points sampled along each element edge.

    Returns
    -------
    pyvista.PolyData
    """
    from ..bspline import PatchEvaluator  # local import: only needed here, avoids a hard bspline dep at module load
    evaluator = PatchEvaluator(patch)

    builder = _PolylineBuilder()
    for spans, params in _iso_line_walk(patch, n_samples):
        pts_ref = patch.evaluate_patch_nd_omp(spans, params)
        u_field = evaluator.evaluate_solution(params, u_global)
        builder.add(pts_ref + factor * u_field[:, :pts_ref.shape[1]])
    return builder.build()


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
    nonlinear_subdivision=1,
    patch=None,
    u_global=None,
    n_edge_samples=20,
    show_control_net=False,
    control_net_color='gray',
    control_net_line_width=1,
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
        If True, overlay the element boundaries on the deformed shape --
        curved (via element_boundary_curves_deformed()) if `patch` and
        `u_global` are given, straight corner-to-corner chords (via
        element_boundary_lines(), only exact when the *deformed* edge
        happens to stay straight) otherwise. This is what you want instead
        of the mesh's own "show edges" toggle, which draws the much finer
        rendering tessellation instead.
    edge_color : default 'black'
        Element boundary line color.
    edge_line_width : float, default 2
        Element boundary line width.
    scalar_bar_title : str, optional
        Scalar bar title. Defaults to ``field_name`` (plus the component
        index, if given).
    smooth_shading : bool, default True
        Gouraud/Phong-interpolated shading across the surface.
    nonlinear_subdivision : int, default 1
        Tessellation level for `mesh`'s curved Bezier cells before warping
        (see ``UnstructuredGrid.extract_surface``'s parameter of the same
        name). The default (1) is coarse -- fine for elements with little
        curvature and little variation in `field_name` (e.g. a barely-warped
        flat plate), but raise it (e.g. 4-5) whenever the *deformed* shape
        has real curvature the coarse tessellation would facet.
    patch : Patch, optional
        The original Patch `mesh` was exported from (with a PatchDOFManager
        attached). Needed, together with `u_global`, for a curved element-
        boundary overlay; without it, the overlay falls back to straight
        corner-to-corner chords.
    u_global : ndarray, shape (n_dof,), optional
        Global solution vector matching `patch`'s dof numbering. See `patch`.
    n_edge_samples : int, default 20
        Points sampled along each element edge when `patch`/`u_global` are
        given (passed to element_boundary_curves_deformed()).
    show_control_net : bool, default False
        If True, also overlay the control net (control_net_lines()) -- the
        straight-line control polygon, not the surface itself. Requires
        `patch`; deformed (by `u_global`, scaled by `factor`) if `u_global`
        is also given, otherwise the reference (undeformed) control net.
    control_net_color : default 'gray'
        Control net line color.
    control_net_line_width : float, default 1
        Control net line width.

    Returns
    -------
    warped : pyvista.PolyData
        The warped, colored surface (already added to `pl`).
    edges : pyvista.PolyData or None
        The element-boundary overlay (already added to `pl`), or None if
        show_element_edges is False.
    control_net : pyvista.PolyData or None
        The control net overlay (already added to `pl`), or None if
        show_control_net is False.
    """
    surf = mesh if nonlinear_subdivision <= 1 else mesh.extract_surface(
        nonlinear_subdivision=nonlinear_subdivision)
    warped = surf.warp_by_vector(field_name, factor=factor)

    title = scalar_bar_title
    if title is None:
        title = field_name if component is None else f'{field_name}[{component}]'

    pl.add_mesh(warped, scalars=field_name, component=component, cmap=cmap,
                smooth_shading=smooth_shading, scalar_bar_args={'title': title})

    edges = None
    if show_element_edges:
        if patch is not None and u_global is not None:
            edges = element_boundary_curves_deformed(
                patch, u_global, factor=factor, n_samples=n_edge_samples)
        else:
            edges = element_boundary_lines(
                mesh, points=mesh.warp_by_vector(field_name, factor=factor).points)
        pl.add_mesh(edges, color=edge_color, line_width=edge_line_width)

    control_net = None
    if show_control_net:
        if patch is None:
            raise ValueError(
                "add_deformed_shell: show_control_net=True requires `patch`.")
        cp_ref = np.array([patch.control_point(i) for i in range(patch.n_cp)])
        cp_points = cp_ref
        if u_global is not None:
            cp_field = np.array([
                u_global[patch.dof_manager.get_global_dof_indices(i)]
                for i in range(patch.n_cp)])
            cp_points = cp_ref + factor * cp_field[:, :cp_ref.shape[1]]
        control_net = control_net_lines(patch, points=cp_points)
        pl.add_mesh(control_net, color=control_net_color, line_width=control_net_line_width)

    return warped, edges, control_net
