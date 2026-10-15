"""
Tests for a previously-unexercised code path: refining a Patch that already
has a PatchDOFManager attached (as opposed to the two-stage "refine a
geometry-only Patch, then attach DOFs afterward" pattern every existing
notebook/test otherwise uses).

Patch.dof_manager is read-only from Python (no setter is bound), so there is
no way to attach a PatchDOFManager to an EXISTING Patch object from Python --
only the C++ refiners can rewrite it in place (direct field access, not
through the binding). PRefiner::refine() and HRefiner::refine() (and
SubdivisionRefiner, which delegates to HRefiner) both contain:

    if (patch.dof_manager) {
        patch.dof_manager = std::make_shared<PatchDOFManager>(
            *patch.dof_manager, patch.global_indices);
    }

i.e. they keep an already-attached DOF manager in sync with the refined
control net automatically. This was discovered (and confirmed to work) while
reviewing examples/future/18_solid_3d_pure_bending.ipynb, which attaches the
DOF manager before k-refining instead of rebuilding one afterward -- a
simplification only possible for a single, isolated patch (the two-stage
pattern is still needed for multi-patch assemblies, where global DOF
numbering must wait until every patch has reached its final refined size).
"""

import numpy as np

# pylint: disable=no-name-in-module
from yeti_iga.future.bspline import (
    BSpline, BSplineSurface, BSplineVolume, ControlPointManager, Patch,
    GlobalDOFManager, PatchDOFManager, PRefiner, HRefiner, SubdivisionRefiner,
    IGABasis1D, Material, PlaneStress, Elasticity3D, PatchIntegrator,
)


def _identity_dof_indices_match(patch, dofs_per_cp):
    """For a single, un-shared patch, global_indices is the identity mapping,
    so get_global_dof_indices(i) must simply be
    [dofs_per_cp*i, dofs_per_cp*i+1, ..., dofs_per_cp*i+dofs_per_cp-1]
    for every control point, both before and after refinement."""
    for i in range(patch.n_cp):
        expected = list(range(dofs_per_cp * i, dofs_per_cp * i + dofs_per_cp))
        assert patch.dof_manager.get_global_dof_indices(i) == expected


def _build_2d_patch_with_dof_manager(p, kv, dofs_per_cp=2):
    n = len(kv) - p - 1
    mgr = ControlPointManager(dim=2)
    for iv in range(n):
        for iu in range(n):
            mgr.add_point([float(iu), float(iv)])
    surf = BSplineSurface(BSpline(p, kv), BSpline(p, kv))

    n_cp = n * n
    gdm = GlobalDOFManager([dofs_per_cp] * n_cp)
    pdm = PatchDOFManager(dofs_per_cp, list(range(n_cp)), gdm)
    return Patch(surf, mgr, list(range(n_cp)), [n, n], pdm)


def test_prefiner_keeps_attached_dof_manager_in_sync():
    kv = np.array([0., 0., 1., 1.])
    patch = _build_2d_patch_with_dof_manager(1, kv)
    _identity_dof_indices_match(patch, 2)

    PRefiner(direction=0, n_elevations=1).refine(patch)
    PRefiner(direction=1, n_elevations=1).refine(patch)

    assert patch.dof_manager is not None
    _identity_dof_indices_match(patch, 2)


def test_hrefiner_keeps_attached_dof_manager_in_sync():
    kv = np.array([0., 0., 0., 1., 1., 1.])
    patch = _build_2d_patch_with_dof_manager(2, kv)
    _identity_dof_indices_match(patch, 2)

    HRefiner(direction=0, knot=0.5).refine(patch)

    assert patch.dof_manager is not None
    _identity_dof_indices_match(patch, 2)


def test_subdivision_refiner_keeps_attached_dof_manager_in_sync():
    kv = np.array([0., 0., 0., 1., 1., 1.])
    patch = _build_2d_patch_with_dof_manager(2, kv)
    _identity_dof_indices_match(patch, 2)

    SubdivisionRefiner(direction=0, n_levels=2).refine(patch)
    SubdivisionRefiner(direction=1, n_levels=2).refine(patch)

    assert patch.dof_manager is not None
    _identity_dof_indices_match(patch, 2)


def test_refine_3d_patch_with_attached_dof_manager_then_assemble():
    """
    End-to-end check mirroring notebook 18's simplified workflow: attach the
    DOF manager to a 3D patch before k-refining it (PRefiner + Subdivision-
    Refiner), then assemble the stiffness matrix directly -- no intermediate
    "rebuild the Patch with a fresh DOF manager after refining" step.
    """
    kv = np.array([0., 0., 1., 1.])
    splines = [BSpline(1, kv) for _ in range(3)]
    tensor = BSplineVolume(*splines)
    mgr = ControlPointManager(dim=3)
    for iw in range(2):
        for iv in range(2):
            for iu in range(2):
                mgr.add_point([float(iu), float(iv), float(iw)])

    gdm = GlobalDOFManager([3] * 8)
    pdm = PatchDOFManager(3, list(range(8)), gdm)
    patch = Patch(tensor, mgr, list(range(8)), [2, 2, 2], pdm)
    _identity_dof_indices_match(patch, 3)

    for d in range(3):
        PRefiner(direction=d, n_elevations=1).refine(patch)
    for d in range(3):
        SubdivisionRefiner(direction=d, n_levels=1).refine(patch)

    assert patch.n_cp == 64
    _identity_dof_indices_match(patch, 3)

    su, sv, sw = patch.tensor.components
    basis_u = IGABasis1D.build(su, su.degree + 1)
    basis_v = IGABasis1D.build(sv, sv.degree + 1)
    basis_w = IGABasis1D.build(sw, sw.degree + 1)
    law = Elasticity3D(Material(E=210000.0, nu=0.3))
    integrator = PatchIntegrator(patch, basis_u, basis_v, basis_w, law)

    K = integrator.integrate_stiffness()
    assert K.shape == (192, 192)
    np.testing.assert_allclose(K.toarray(), K.toarray().T, atol=1.e-6)


def test_refined_dof_manager_matches_rebuilt_from_scratch():
    """
    The DOF manager kept in sync by refining an attached one must be
    EQUIVALENT (same DOF-index assignment) to the standard two-stage
    pattern -- refine a geometry-only patch, then build a fresh DOF manager
    from scratch afterward -- confirming this is not merely "didn't crash"
    but produces the exact same DOF numbering either way.
    """
    kv = np.array([0., 0., 0., 1., 1., 1.])

    # Path A: DOF manager attached before refining.
    patch_a = _build_2d_patch_with_dof_manager(2, kv, dofs_per_cp=2)
    PRefiner(direction=0, n_elevations=1).refine(patch_a)
    SubdivisionRefiner(direction=1, n_levels=1).refine(patch_a)

    # Path B: geometry-only patch, DOF manager rebuilt from scratch after refining.
    mgr_b = ControlPointManager(dim=2)
    n0 = len(kv) - 2 - 1
    for iv in range(n0):
        for iu in range(n0):
            mgr_b.add_point([float(iu), float(iv)])
    patch_b = Patch(BSplineSurface(BSpline(2, kv), BSpline(2, kv)), mgr_b,
                     list(range(n0 * n0)), [n0, n0])
    PRefiner(direction=0, n_elevations=1).refine(patch_b)
    SubdivisionRefiner(direction=1, n_levels=1).refine(patch_b)

    gdm_b = GlobalDOFManager([2] * patch_b.n_cp)
    pdm_b = PatchDOFManager(2, list(range(patch_b.n_cp)), gdm_b)
    patch_b = Patch(patch_b.tensor, patch_b.cp_manager,
                     list(patch_b.global_indices), list(patch_b.local_shape), pdm_b)

    assert patch_a.n_cp == patch_b.n_cp
    for i in range(patch_a.n_cp):
        assert (patch_a.dof_manager.get_global_dof_indices(i)
                == patch_b.dof_manager.get_global_dof_indices(i))


if __name__ == '__main__':
    test_prefiner_keeps_attached_dof_manager_in_sync()
    test_hrefiner_keeps_attached_dof_manager_in_sync()
    test_subdivision_refiner_keeps_attached_dof_manager_in_sync()
    test_refine_3d_patch_with_attached_dof_manager_then_assemble()
    test_refined_dof_manager_matches_rebuilt_from_scratch()
