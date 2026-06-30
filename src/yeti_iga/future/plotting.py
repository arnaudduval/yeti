"""
2D visualization helpers for B-spline/NURBS patches.

Generalizes the ad-hoc plotting code from the 02_patch_basics notebook into
a single reusable function, usable on one patch, a list of patches, or a
PatchAssembly (e.g. to visualize a multi-patch assembly on one set of axes).

plot_dirichlet_bc() and plot_distributed_load() add boundary-condition
markers on top of an existing plot_patches_2d() axes -- one for Dirichlet
(fixed control points), one for Neumann (distributed edge loads).
"""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt

from .bspline import Patch, PatchAssembly


def _as_patch_list(patches):
    """Normalize a Patch, an iterable of Patch, or a PatchAssembly into a list of Patch."""
    if isinstance(patches, PatchAssembly):
        return list(patches.get_patchs())
    if isinstance(patches, Patch):
        return [patches]
    return list(patches)


def plot_patches_2d(
    patches,
    n_samples=100,
    show_control_points=False,
    show_control_point_indices=False,
    ax=None,
    title=None,
):
    """
    Plot the iso-parametric element borders of one or several 2D patches.

    Parameters
    ----------
    patches : Patch, iterable of Patch, or PatchAssembly
        Patch(es) to plot. Every patch is drawn on the same axes, each in
        its own color -- handy to visualize a multi-patch assembly and
        check that shared edges line up.
    n_samples : int, default 100
        Number of points sampled along each iso-parametric line.
    show_control_points : bool, default False
        If True, draw each patch's control points as square markers.
    show_control_point_indices : bool, default False
        If True, annotate each control point with its global id in the
        shared ControlPointManager pool.
    ax : matplotlib.axes.Axes, optional
        Axes to draw on. A new figure is created and shown if not given.
    title : str, optional
        Plot title.

    Returns
    -------
    ax : matplotlib.axes.Axes
    """
    patch_list = _as_patch_list(patches)
    if not patch_list:
        raise ValueError("plot_patches_2d: no patch to plot.")

    created_fig = ax is None
    if created_fig:
        _, ax = plt.subplots()

    colors = plt.cm.tab10(np.arange(len(patch_list)) % 10)

    for patch, color in zip(patch_list, colors):
        if len(patch.tensor.components) != 2:
            raise ValueError("plot_patches_2d only supports 2D patches.")

        su = patch.tensor.components[0]
        sv = patch.tensor.components[1]
        u_breaks = np.unique(su.knot_vector)
        v_breaks = np.unique(sv.knot_vector)

        # iso-u lines
        for u_val in u_breaks:
            vs = np.linspace(v_breaks[0], v_breaks[-1], n_samples)
            spans = np.array(
                [[su.find_span(u_val), sv.find_span(v)] for v in vs], dtype=np.int32
            )
            params = np.column_stack([np.full(n_samples, u_val), vs])
            pts = patch.evaluate_patch_nd_omp(spans, params)
            ax.plot(pts[:, 0], pts[:, 1], '-', lw=1, color=color)

        # iso-v lines
        for v_val in v_breaks:
            us = np.linspace(u_breaks[0], u_breaks[-1], n_samples)
            spans = np.array(
                [[su.find_span(u), sv.find_span(v_val)] for u in us], dtype=np.int32
            )
            params = np.column_stack([us, np.full(n_samples, v_val)])
            pts = patch.evaluate_patch_nd_omp(spans, params)
            ax.plot(pts[:, 0], pts[:, 1], '-', lw=1, color=color)

        if show_control_points or show_control_point_indices:
            cps = np.array([patch.control_point(i) for i in range(patch.n_cp)])

            if show_control_points:
                ax.scatter(cps[:, 0], cps[:, 1], marker='s', s=30,
                           color=color, zorder=3)

            if show_control_point_indices:
                for local_i, (x, y) in zip(patch.global_indices, cps):
                    ax.annotate(str(local_i), (x, y), textcoords="offset points",
                                xytext=(6, 6), color=color, fontsize=8)

    ax.set_aspect('equal')
    ax.grid(True)
    if title:
        ax.set_title(title)

    if created_fig:
        plt.tight_layout()
        plt.show()

    return ax


def plot_dirichlet_bc(ax, patch, local_cp_indices, dofs=(0, 1), size=80, color='red'):
    """
    Draw simple-support markers at given control points for a Dirichlet
    (displacement) boundary condition.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to draw on (e.g. the one returned by plot_patches_2d()).
    patch : Patch
        Patch the control points belong to.
    local_cp_indices : iterable of int
        Patch-LOCAL control point positions to mark, e.g. from
        Patch.boundary_control_points() (edge/span-range granularity) or any
        other selection (control-point granularity).
    dofs : tuple of int, default (0, 1)
        Which dof components are fixed at every point in local_cp_indices:
        (0, 1) for both x and y (drawn as a filled square), (0,) for x only
        (a left-pointing triangle), (1,) for y only (a downward triangle).
        Call this function once per dof subset if different points in a
        selection are fixed differently.
    size : float, default 80
        Marker size (as in matplotlib's `s` scatter parameter).
    color : default 'red'
        Marker fill color.

    Returns
    -------
    ax : matplotlib.axes.Axes
    """
    dofs = tuple(sorted(set(dofs)))
    marker = {(0, 1): 's', (0,): '<', (1,): 'v'}.get(dofs)
    if marker is None:
        raise ValueError(
            f"plot_dirichlet_bc: unsupported dofs {dofs} -- use (0,), (1,), or (0, 1).")

    local_cp_indices = list(local_cp_indices)
    if not local_cp_indices:
        return ax

    pts = np.array([patch.control_point(i) for i in local_cp_indices])
    ax.scatter(pts[:, 0], pts[:, 1], marker=marker, s=size,
               facecolor=color, edgecolor='black', zorder=5)
    return ax


def plot_distributed_load(
    ax,
    patch,
    direction,
    side,
    traction,
    span_min=-1,
    span_max=-1,
    n_arrows=7,
    scale=None,
    color='C3',
):
    """
    Draw a distributed (Neumann) boundary load as a row of arrows sampled
    along a patch edge.

    Parameters mirror PatchIntegrator.integrate_boundary_load(): direction
    and side select the edge (fixing `direction` at its first (side=0) or
    last (side=1) parameter), and span_min/span_max optionally restrict the
    arrows to a sub-range of the edge (raw knot-span indices of the OTHER
    direction, like Patch.boundary_control_points()).

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to draw on (e.g. the one returned by plot_patches_2d()).
    patch : Patch
        2D patch the loaded edge belongs to.
    direction, side : int
        Edge selection, see above.
    traction : Traction
        Evaluated at each sampled point (traction.evaluate(point)) -- works
        for ConstantTraction as well as any future spatially-varying kernel.
    span_min, span_max : int, default -1, -1
        Raw knot-span range of the OTHER direction to restrict the arrows
        to; -1/-1 (default) samples the whole edge.
    n_arrows : int, default 7
        Number of arrows sampled along the edge.
    scale : float, optional
        Arrow length per unit of traction magnitude. If not given, picked
        automatically so the largest arrow spans ~15% of the patch's
        bounding-box diagonal.
    color : default 'C3'
        Arrow color.

    Returns
    -------
    ax : matplotlib.axes.Axes
    """
    if len(patch.tensor.components) != 2:
        raise ValueError("plot_distributed_load only supports 2D patches.")

    varying = 1 - direction
    s_fixed = patch.tensor.components[direction]
    s_varying = patch.tensor.components[varying]

    u_fixed = s_fixed.knot_vector[0] if side == 0 else s_fixed.knot_vector[-1]

    kv_varying = s_varying.knot_vector
    u_lo, u_hi = kv_varying[0], kv_varying[-1]
    if span_min >= 0:
        u_lo, u_hi = kv_varying[span_min], kv_varying[span_max + 1]

    params_varying = np.linspace(u_lo, u_hi, n_arrows)

    points = np.empty((n_arrows, 2))
    vectors = np.empty((n_arrows, 2))
    for i, u_var in enumerate(params_varying):
        span_uv = [0, 0]
        param_uv = [0.0, 0.0]
        span_uv[direction] = s_fixed.find_span(u_fixed)
        span_uv[varying] = s_varying.find_span(u_var)
        param_uv[direction] = u_fixed
        param_uv[varying] = u_var

        pt = patch.evaluate_patch_nd_omp(
            np.array([span_uv], dtype=np.int32), np.array([param_uv]))[0]
        points[i] = pt
        vectors[i] = traction.evaluate(pt)

    if scale is None:
        max_mag = np.linalg.norm(vectors, axis=1).max()
        cps = np.array([patch.control_point(i) for i in range(patch.n_cp)])
        diag = np.linalg.norm(cps.max(axis=0) - cps.min(axis=0))
        scale = 0.15 * diag / max(max_mag, 1.e-12)

    ax.quiver(points[:, 0], points[:, 1], vectors[:, 0] * scale, vectors[:, 1] * scale,
              angles='xy', scale_units='xy', scale=1, color=color, width=0.005, zorder=4)
    return ax
