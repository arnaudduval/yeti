"""
Regression test for a lifetime bug discovered while writing the 3D-solid
tests: PatchIntegrator (and PatchEvaluator) store C++ REFERENCES to their
`patch`/`basis_u`/`basis_v`/`basis_w` constructor arguments (not copies), so
without py::keep_alive on the pybind11 constructor bindings, Python is free
to garbage-collect those objects as soon as no Python variable references
them -- even while the object built from them is still alive and in use.
Before the fix this crashed with an "unordered_map::at" exception (reading
basis_u_.span_indices.at(...) through a dangling reference) the first time
integrate_stiffness()/integrate_mass()/etc. ran after the arguments were
collected. It reproduced deterministically (not just occasionally) whenever
patch/basis construction was factored into a helper function that returned
only the integrator -- exactly the pattern test_stiffness_3d.py's
_build_integrator() and friends use.

These tests build the integrator/evaluator inside a nested function that
returns only that object (discarding every local reference to patch/basis_u/
v/w), then force garbage collection before using it -- the worst case for
the bug.
"""

import gc

import numpy as np

# pylint: disable=no-name-in-module
from yeti_iga.future.bspline import (
    BSpline, BSplineSurface, BSplineVolume, ControlPointManager, Patch,
    GlobalDOFManager, PatchDOFManager, IGABasis1D, Material, PlaneStress,
    Elasticity3D, PatchIntegrator, PatchEvaluator,
)


def _make_2d_integrator():
    """Patch and basis_u/basis_v are local to this function; only the
    integrator is returned."""
    mgr = ControlPointManager(dim=2)
    mgr.add_point([0.0, 0.0])
    mgr.add_point([1.0, 0.0])
    mgr.add_point([0.0, 1.0])
    mgr.add_point([1.0, 1.0])

    dof_manager = GlobalDOFManager([2] * mgr.n_points)
    su = BSpline(1, np.array([0., 0., 1., 1.]))
    sv = BSpline(1, np.array([0., 0., 1., 1.]))
    surf = BSplineSurface(su, sv)
    mapping = [0, 1, 2, 3]
    dof_manager_patch = PatchDOFManager(2, mapping, dof_manager)
    patch = Patch(surf, mgr, mapping, [2, 2], dof_manager_patch)

    basis_u = IGABasis1D.build(su, 2)
    basis_v = IGABasis1D.build(sv, 2)

    return PatchIntegrator(patch, basis_u, basis_v, PlaneStress(Material(E=210000, nu=0.3)))


def _make_3d_integrator():
    """Patch and basis_u/basis_v/basis_w are local to this function; only
    the integrator is returned."""
    kv = np.array([0., 0., 0., 1., 1., 1.])
    splines = [BSpline(2, kv) for _ in range(3)]
    tensor = BSplineVolume(*splines)

    mgr = ControlPointManager(dim=3)
    n1 = 3
    for kw in range(n1):
        for jv in range(n1):
            for iu in range(n1):
                mgr.add_point([float(iu), float(jv), float(kw)])

    n_cp = n1 ** 3
    gdm = GlobalDOFManager([3] * n_cp)
    pdm = PatchDOFManager(3, list(range(n_cp)), gdm)
    patch = Patch(tensor, mgr, list(range(n_cp)), [n1, n1, n1], pdm)

    basis_u = IGABasis1D.build(patch.tensor.components[0], 3)
    basis_v = IGABasis1D.build(patch.tensor.components[1], 3)
    basis_w = IGABasis1D.build(patch.tensor.components[2], 3)

    return PatchIntegrator(patch, basis_u, basis_v, basis_w, Elasticity3D(Material(E=210000, nu=0.3)))


def test_2d_constructor_keeps_patch_and_bases_alive():
    integrator = _make_2d_integrator()
    gc.collect()  # force collection of the now-unreferenced patch/basis_u/basis_v
    K = integrator.integrate_stiffness()
    assert K.shape == (8, 8)


def test_3d_constructor_keeps_patch_and_bases_alive():
    integrator = _make_3d_integrator()
    gc.collect()  # force collection of the now-unreferenced patch/basis_u/v/w
    K = integrator.integrate_stiffness()
    assert K.shape == (81, 81)


def _make_evaluator():
    """Patch is local to this function; only the evaluator is returned."""
    mgr = ControlPointManager(dim=2)
    mgr.add_point([0.0, 0.0])
    mgr.add_point([1.0, 0.0])
    mgr.add_point([0.0, 1.0])
    mgr.add_point([1.0, 1.0])

    dof_manager = GlobalDOFManager([2] * mgr.n_points)
    su = BSpline(1, np.array([0., 0., 1., 1.]))
    sv = BSpline(1, np.array([0., 0., 1., 1.]))
    surf = BSplineSurface(su, sv)
    mapping = [0, 1, 2, 3]
    dof_manager_patch = PatchDOFManager(2, mapping, dof_manager)
    patch = Patch(surf, mgr, mapping, [2, 2], dof_manager_patch)

    return PatchEvaluator(patch)


def test_evaluator_constructor_keeps_patch_alive():
    evaluator = _make_evaluator()
    gc.collect()  # force collection of the now-unreferenced patch
    params = np.array([[0.5, 0.5]])
    u_global = np.zeros(8)
    result = evaluator.evaluate_solution(params, u_global)
    assert result.shape == (1, 2)


if __name__ == '__main__':
    test_2d_constructor_keeps_patch_and_bases_alive()
    test_3d_constructor_keeps_patch_and_bases_alive()
    test_evaluator_constructor_keeps_patch_alive()
