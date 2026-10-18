"""
Test the NURBS (rational) path of the 3D solid (BSplineVolume) kernels:
PatchIntegrator's stiffness/mass/boundary-load kernels, dispatched at
runtime on patch.cp_manager.is_rational() -- see
evaluateGaussPointGeometry3DNURBS and friends in PatchIntegrator.cpp.

No legacy Fortran 3D reference exists (same situation as the B-spline-only
3D solid work), and unlike the Saint-Venant pure-bending validation case
(test_stiffness_3d.py), a NURBS-specific check needs a genuinely CURVED
geometry with non-trivial weights -- a uniformly-weighted NURBS patch
reduces identically to the B-spline case (the quotient rule's W(xi) is
then constant, so dW/dxi == 0 and nothing new is exercised).

Validation geometry: an OCTANT OF A HOLLOW SPHERE (a spherical shell
wedge), built via the standard "surface of revolution" NURBS technique
(Piegl & Tiller): a quarter-circle PROFILE curve (polar angle u, weights
[1, 1/sqrt(2), 1]) swept through a quarter-circle ROTATION (azimuthal
angle v, SAME weight pattern) about the z-axis, with the radial
(thickness, w) direction a plain affine B-spline (no weights needed --
the natural choice for an offset/thickness direction, exactly as in the
cylinder case this geometry replaces).

Crucially, W(u,v) = W_u(u) * W_v(v) with BOTH factors genuinely
non-constant, so dW/du depends on v (through the W_v(v) factor) and dW/dv
depends on u -- the two quotient-rule terms are *coupled*, unlike a
geometry rational in only one direction (where the other direction's dW
is identically zero and the cross-coupling is never exercised). Confirmed
directly: W(0.1, 0.1) != W(0.1, 0.3) * W_u(0.1)/W_u(0.1) trivially, but the
real check is dW/du at v=0.1 vs v=0.3 differing (verified in the research
that produced this file).
"""

import numpy as np

# pylint: disable=no-name-in-module
from yeti_iga.future.bspline import (
    BSpline, BSplineVolume, Patch, ControlPointManager,
    GlobalDOFManager, PatchDOFManager, IGABasis1D, Material, Elasticity3D,
    PatchIntegrator, ConstantTraction3D,
)

RIN, REX = 1.0, 2.0
W0 = 1.0 / np.sqrt(2.0)
RHO = 7800.0

# Quarter-circle PROFILE (polar angle u: 0=pole on the z-axis, 1=equator on
# the x-axis), radius R: control points (x, z, weight).
def _profile(R):
    return [(0.0, R, 1.0), (R, R, W0), (R, 0.0, 1.0)]


# Quarter-circle ROTATION (azimuthal angle v: 0=x-axis plane, 1=y-axis
# plane): control points (Rx, Ry, weight) -- same weight pattern as the
# profile, applied to the in-plane direction instead of a radius.
_ROT = [(1.0, 0.0, 1.0), (1.0, 1.0, W0), (0.0, 1.0, 1.0)]


def _build_shell_octant_patch():
    """
    Surface-of-revolution construction: P(iu, iv) = (x_p(iu)*Rx(iv),
    x_p(iu)*Ry(iv), z_p(iu)), weight = w_p(iu)*w_rot(iv) -- the standard
    technique for sweeping a profile curve through a circular arc exactly
    (Piegl & Tiller, "The NURBS Book", Sec. 8.4). One copy at R=Rin
    (w=0 layer), one at R=Rex (w=1 layer); the w direction is a plain
    affine (non-rational) B-spline. local_shape = [3, 3, 2], u-fastest.
    """
    mgr = ControlPointManager(dim=3)
    n_cp = 3 * 3 * 2
    for radius in (RIN, REX):
        profile = _profile(radius)
        for iv in range(3):
            rx, ry, wv = _ROT[iv]
            for iu in range(3):
                x_p, z_p, wu = profile[iu]
                mgr.add_point([x_p * rx, x_p * ry, z_p], w=wu * wv)

    su = BSpline(2, np.array([0., 0., 0., 1., 1., 1.]))   # polar angle, degree 2
    sv = BSpline(2, np.array([0., 0., 0., 1., 1., 1.]))   # azimuthal angle, degree 2
    sw = BSpline(1, np.array([0., 0., 1., 1.]))            # radial thickness, degree 1
    tensor = BSplineVolume(su, sv, sw)

    gdm = GlobalDOFManager([3] * n_cp)
    pdm = PatchDOFManager(3, list(range(n_cp)), gdm)
    patch = Patch(tensor, mgr, list(range(n_cp)), [3, 3, 2], pdm)
    return patch, n_cp


def _build_integrator(patch, law, ngauss_u=None, ngauss_v=None, ngauss_w=None):
    su, sv, sw = patch.tensor.components
    basis_u = IGABasis1D.build(su, ngauss_u or su.degree + 1)
    basis_v = IGABasis1D.build(sv, ngauss_v or sv.degree + 1)
    basis_w = IGABasis1D.build(sw, ngauss_w or sw.degree + 1)
    return PatchIntegrator(patch, basis_u, basis_v, basis_w, law)


def test_shell_octant_is_rational():
    patch, _ = _build_shell_octant_patch()
    assert patch.cp_manager.is_rational


def test_shell_octant_weight_couples_both_angular_directions():
    """
    Confirms (independently of PatchIntegrator) that this geometry's
    W(u,v) genuinely couples u and v -- dW/du depends on v and vice versa
    -- unlike a geometry rational in only one direction. This is the
    property that makes this test geometry exercise the 3D quotient rule's
    cross terms, which a single-direction-rational geometry cannot.
    """
    patch, _ = _build_shell_octant_patch()
    su, sv, _ = patch.tensor.components
    w = np.array([1.0, W0, 1.0])

    def W_u(u):
        return w @ np.array([su.one_basis_fun(u, a) for a in range(3)])

    def W_v(v):
        return w @ np.array([sv.one_basis_fun(v, a) for a in range(3)])

    # W_u (and W_v, identical pattern) must be non-constant...
    assert not np.isclose(W_u(0.1), W_u(0.5))
    # ...so that dW/du (proportional to W_v(v)) genuinely depends on v.
    assert not np.isclose(W_v(0.1), W_v(0.3))


def test_nurbs3d_stiffness_symmetric_rigid_body_modes():
    patch, n_cp = _build_shell_octant_patch()
    law = Elasticity3D(Material(E=210000.0, nu=0.3))
    integrator = _build_integrator(patch, law)

    K = integrator.integrate_stiffness()
    np.testing.assert_allclose(K.toarray(), K.toarray().T, atol=1.e-6)

    coords = np.array([patch.control_point(i) for i in range(n_cp)])
    center = coords.mean(axis=0)

    for t in ([1., 0., 0.], [0., 1., 0.], [0., 0., 1.]):
        u = np.tile(t, n_cp)
        assert abs(u @ (K @ u)) < 1.e-6

    for axis in ([1., 0., 0.], [0., 1., 0.], [0., 0., 1.]):
        u = np.zeros(3 * n_cp)
        for i in range(n_cp):
            u[3*i:3*i+3] = np.cross(axis, coords[i] - center)
        energy = u @ (K @ u)
        assert abs(energy) < 1.e-6 * max(np.linalg.norm(u) ** 2, 1.0)


def test_nurbs3d_mass_matches_exact_volume():
    """
    u^T M u == rho * Volume for a rigid unit-translation field, with
    Volume = (pi/6)*(Rex^3-Rin^3) (closed-form octant-of-hollow-sphere
    volume). NURBS basis functions are rational, not polynomial, so this
    is a CONVERGENCE check (more Gauss points -> smaller error), not an
    exact equality -- verified by checking the error shrinks monotonically
    as ngauss increases, then requiring a tight relative error at a
    generously oversampled quadrature order.
    """
    patch, n_cp = _build_shell_octant_patch()
    law = Elasticity3D(Material(E=210000.0, nu=0.3, rho=RHO))
    exact_volume = (np.pi / 6.0) * (REX**3 - RIN**3)
    exact_mass_x = RHO * exact_volume

    u = np.tile([1.0, 0.0, 0.0], n_cp)

    rel_errors = []
    for ngauss in (3, 6, 10):
        integrator = _build_integrator(patch, law, ngauss_u=ngauss, ngauss_v=ngauss)
        M = integrator.integrate_mass()
        mass_x = u @ (M @ u)
        rel_errors.append(abs(mass_x - exact_mass_x) / exact_mass_x)

    assert rel_errors[0] > rel_errors[1] > rel_errors[2]
    assert rel_errors[-1] < 1.e-9


def test_nurbs3d_boundary_load_resultant_converges_on_coupled_face():
    """
    The v=0 face (direction=1, side=0) lies in the FLAT y=0 plane (the
    rotation's v=0 control points all have Ry=0), but is still CURVED
    within that plane (a quarter-circle arc swept by the polar angle u) --
    so, unlike the cylinder geometry's flat-and-non-rational face, this
    resultant is also a convergence check, not an exact equality.
    Area = (pi/4)*(Rex^2-Rin^2) (quarter-annulus cross-section).
    """
    patch, _ = _build_shell_octant_patch()
    law = Elasticity3D(Material(E=210000.0, nu=0.3))

    p = 1000.0
    traction = ConstantTraction3D(np.array([0.0, 0.0, -p]))
    exact_area = (np.pi / 4.0) * (REX**2 - RIN**2)
    exact_force = -p * exact_area

    rel_errors = []
    for ngauss in (3, 6, 10):
        integrator = _build_integrator(patch, law, ngauss_u=ngauss, ngauss_v=ngauss)
        F = integrator.integrate_boundary_load_3d(1, 0, traction)
        np.testing.assert_allclose(F[0::3].sum(), 0.0, atol=1.e-8)
        np.testing.assert_allclose(F[1::3].sum(), 0.0, atol=1.e-8)
        rel_errors.append(abs(F[2::3].sum() - exact_force) / abs(exact_force))

    assert rel_errors[0] > rel_errors[1] > rel_errors[2]
    assert rel_errors[-1] < 1.e-9


if __name__ == '__main__':
    test_shell_octant_is_rational()
    test_shell_octant_weight_couples_both_angular_directions()
    test_nurbs3d_stiffness_symmetric_rigid_body_modes()
    test_nurbs3d_mass_matches_exact_volume()
    test_nurbs3d_boundary_load_resultant_converges_on_coupled_face()
