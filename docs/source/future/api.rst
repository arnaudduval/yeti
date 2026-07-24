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

.. autoclass:: yeti_iga.future.bspline.Material
    :special-members: __init__
    :members:

.. autoclass:: yeti_iga.future.bspline.ConstitutiveLaw
    :special-members: __init__
    :members:

.. autoclass:: yeti_iga.future.bspline.PlaneStress
    :special-members: __init__
    :members:

.. autoclass:: yeti_iga.future.bspline.PlaneStrain
    :special-members: __init__
    :members:

.. autoclass:: yeti_iga.future.bspline.PatchIntegrator
    :special-members: __init__
    :members:

.. autoclass:: yeti_iga.future.bspline.LocalOperator
    :special-members: __init__
    :members:

.. autoclass:: yeti_iga.future.bspline.ScalarLocalOperator
    :special-members: __init__
    :members:

Solution evaluation
----------------------

.. autoclass:: yeti_iga.future.bspline.PatchEvaluator
    :special-members: __init__
    :members:

Boundary conditions and loads
---------------------------------

.. autoclass:: yeti_iga.future.bspline.Traction
    :special-members: __init__
    :members:

.. autoclass:: yeti_iga.future.bspline.ConstantTraction
    :special-members: __init__
    :members:

.. autoclass:: yeti_iga.future.bspline.BoundaryLoadSpec
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
