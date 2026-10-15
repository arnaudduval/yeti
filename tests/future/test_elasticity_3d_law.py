"""
Test the Elasticity3D constitutive law: n_dofs_per_cp()==3, and the B-free
stiffness_density() formula (Planas, Romero & Sancho 2012) directly for
3-vectors -- the same formula PlaneStress/PlaneStrain already use for n=2,
confirmed here unmodified at n=3 (stiffness_density() reads n=v.size() at
runtime, see ConstitutiveLaw.hpp).
"""

import numpy as np

# pylint: disable=no-name-in-module
from yeti_iga.future.bspline import Material, Elasticity3D


def test_elasticity3d_n_dofs_per_cp():
    law = Elasticity3D(Material(E=210000, nu=0.3))
    assert law.n_dofs_per_cp() == 3


def test_elasticity3d_matches_bfree_formula_directly():
    """
    C{v,w} = lambda*(v w^T) + mu*(w v^T) + mu*(v.w)*I, with
    lambda = lambda_3d() = nu*E/((1+nu)*(1-2*nu)), mu = E/(2*(1+nu)).
    """
    E, nu = 210000.0, 0.3
    law = Elasticity3D(Material(E=E, nu=nu))
    mu = E / (2.0 * (1.0 + nu))
    lam = nu * E / ((1.0 + nu) * (1.0 - 2.0 * nu))

    rng = np.random.default_rng(2)
    for _ in range(5):
        v = rng.normal(size=3)
        w = rng.normal(size=3)
        x = rng.normal(size=3)
        expected = (lam * np.outer(v, w) + mu * np.outer(w, v)
                    + mu * np.dot(v, w) * np.eye(3))
        actual = law.stiffness_density(v, w, x)
        np.testing.assert_allclose(actual, expected, atol=1e-8)


def test_elasticity3d_stiffness_density_symmetry_invariant():
    """
    stiffness_density(a, b, x).T == stiffness_density(b, a, x) -- the
    invariant computeLocalStiffnessContribution3D relies on to only compute
    the a<=b pairs and mirror the rest (see PatchIntegrator.cpp).
    """
    law = Elasticity3D(Material(E=210000, nu=0.3))
    rng = np.random.default_rng(1)

    for _ in range(5):
        a = rng.normal(size=3)
        b = rng.normal(size=3)
        x = rng.normal(size=3)
        Cab = law.stiffness_density(a, b, x)
        Cba = law.stiffness_density(b, a, x)
        np.testing.assert_allclose(Cab.T, Cba, atol=1e-10)


def test_elasticity3d_material_accessors_match_lambda_3d():
    E, nu = 210000.0, 0.3
    material = Material(E=E, nu=nu)
    law = Elasticity3D(material)

    assert np.isclose(law.material.mu(), E / (2.0 * (1.0 + nu)))
    assert np.isclose(law.material.lambda_3d(),
                       nu * E / ((1.0 + nu) * (1.0 - 2.0 * nu)))


if __name__ == '__main__':
    test_elasticity3d_n_dofs_per_cp()
    test_elasticity3d_matches_bfree_formula_directly()
    test_elasticity3d_stiffness_density_symmetry_invariant()
    test_elasticity3d_material_accessors_match_lambda_3d()
