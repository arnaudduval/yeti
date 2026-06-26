API reference
=============

.. contents:: Table of contents
    :depth: 1
    :local:
    :backlinks: none

This page is generated from the docstrings of the compiled :mod:`yeti_iga.future.bspline`
pybind11 extension module. It documents *what* each class and function does and its
signature. For *why* the module is structured this way (shared control point pools,
protected vs. private control points during refinement, the two-level DOF manager
split, the multipatch assembly call-order contract...), see :doc:`design`.

.. note::

    A few methods are not yet documented at the binding level (missing docstrings in
    ``bindings.cpp``): :class:`~yeti_iga.future.bspline.BSpline`'s ``find_span`` and
    ``basis_funs``, and :class:`~yeti_iga.future.bspline.HRefiner`'s ``refine_1d`` and
    ``get_type``. This is a documentation gap to close in the bindings themselves, not
    something this page can paper over.

B-spline basis
---------------

.. autoclass:: yeti_iga.future.bspline.BSpline
    :special-members: __init__
    :members:

.. autoclass:: yeti_iga.future.bspline.BSplineTensor
    :special-members: __init__
    :members:

.. autoclass:: yeti_iga.future.bspline.BSplineSurface
    :special-members: __init__
    :members:

.. autoclass:: yeti_iga.future.bspline.BSplineVolume
    :special-members: __init__
    :members:

Control points and patches
---------------------------

.. autoclass:: yeti_iga.future.bspline.ControlPointManager
    :special-members: __init__
    :members:

.. autoclass:: yeti_iga.future.bspline.Patch
    :special-members: __init__
    :members:

DOF management
----------------

.. autoclass:: yeti_iga.future.bspline.GlobalDOFManager
    :special-members: __init__
    :members:

.. autoclass:: yeti_iga.future.bspline.PatchDOFManager
    :special-members: __init__
    :members:

Refinement operators
----------------------

.. autoclass:: yeti_iga.future.bspline.RefinementOperator
    :special-members: __init__
    :members:

.. autoclass:: yeti_iga.future.bspline.HRefiner
    :special-members: __init__
    :members:

.. autoclass:: yeti_iga.future.bspline.SubdivisionRefiner
    :special-members: __init__
    :members:

.. autoclass:: yeti_iga.future.bspline.PRefiner
    :special-members: __init__
    :members:

.. autofunction:: yeti_iga.future.bspline.nd_transition_from_1d

Multipatch assembly
----------------------

.. autoclass:: yeti_iga.future.bspline.PatchAssembly
    :special-members: __init__
    :members:

Integration and assembly
---------------------------

.. autoclass:: yeti_iga.future.bspline.SpanIterator
    :special-members: __init__
    :members:

.. autoclass:: yeti_iga.future.bspline.SpanGauss1D
    :special-members: __init__
    :members:

.. autoclass:: yeti_iga.future.bspline.IGABasis1D
    :special-members: __init__
    :members:

.. autoclass:: yeti_iga.future.bspline.MaterialProperties
    :special-members: __init__
    :members:

.. autoclass:: yeti_iga.future.bspline.PatchIntegrator
    :special-members: __init__
    :members:

.. autoclass:: yeti_iga.future.bspline.LocalOperator
    :special-members: __init__
    :members:

Bézier extraction
--------------------

.. autoclass:: yeti_iga.future.bspline.BezierElementND
    :special-members: __init__
    :members:

.. autoclass:: yeti_iga.future.bspline.BezierExtractor
    :special-members: __init__
    :members:
