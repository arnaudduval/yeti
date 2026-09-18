"""
Post-processing and visualization helpers for the ``future`` module.

Three independent, separately-imported submodules, kept apart deliberately so
that each one's dependency footprint stays isolated -- nothing here is
imported eagerly by this package's __init__:

- ``vtu``: VTU export as high-order VTK Bezier cells. No dependency beyond
  numpy.
- ``plotting``: 2D matplotlib visualization of B-spline/NURBS patches.
  Requires matplotlib (``pip install -e ".[viz]"``).
- ``pv_plotting``: 3D pyvista-based rendering (deformed shapes, interactive
  Jupyter widgets). Requires pyvista (``pip install -e ".[viz]"``).

Import what you need directly, e.g.::

    from yeti_iga.future.postprocessing.vtu import write_bezier_patch_vtu
    from yeti_iga.future.postprocessing.plotting import plot_patches_2d
    from yeti_iga.future.postprocessing.pv_plotting import add_deformed_shell
"""
