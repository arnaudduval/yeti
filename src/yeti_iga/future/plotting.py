"""
2D visualization helpers for B-spline/NURBS patches.

Generalizes the ad-hoc plotting code from the 02_patch_basics notebook into
a single reusable function, usable on one patch, a list of patches, or a
PatchAssembly (e.g. to visualize a multi-patch assembly on one set of axes).
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
