#pragma once
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <string>
#include "Patch.hpp"
#include "ConstitutiveLaw.hpp"
#include "WeightedQuadrature1D.hpp"

// Matrix-free elasticity stiffness action K @ v, for a single (standalone)
// 2D solid Patch with an isotropic PlaneStress/PlaneStrain law, using
// weighted quadrature (WeightedQuadrature1D, Phase A) and the 2D Kronecker
// sandwich product (matrix_free_apply_2d, Phase B) instead of ever
// assembling K. Direct C++ port of pymfiga's
// iga/single_model/mechanical.py: MechanicalModel.compute_stiffness_property
// + compute_mf_stiffness -- with the material-tensor convention corrected to
// match `future`'s OWN validated ConstitutiveLaw::stiffness_density (see the
// .cpp for the derivation: `future`'s B-free identity
// C{v,w} = lambda*(v x w^T) + mu*(w x v^T) + mu*(v.w)*I corresponds to the
// 4th-order tensor D_IJlm = lambda*d_Il*d_Jm + mu*d_Im*d_Jl + mu*d_IJ*d_lm,
// NOT the textbook C_ijkl = lambda*d_ij*d_kl + mu*(d_ik*d_jl+d_il*d_jk) that
// pymfiga's own material classes build -- the two differ termwise (though
// both are mathematically valid elasticity tensors under full strain
// symmetrization), and only the former reproduces
// PatchIntegrator::integrateStiffness() bit-for-bit on this codebase's
// existing, legacy-Fortran-validated kernel.
//
// Scope: 2D only (n_dofs_per_cp() == 2, i.e. PlaneStress/PlaneStrain),
// non-periodic tensor-product patch, standalone (not part of a multi-patch
// PatchAssembly -- global dof indices are read via patch.dof_manager, so a
// shared-dof patch would still index correctly, but this class does not
// merge contributions across patches the way PatchIntegrator::assembleStiffness does).
//
// KNOWN CAVEAT, discovered while building an iterative matrix-free solver on
// top of apply() (see yeti_iga.future.matrix_free_solver): on CURVED
// geometry, apply() is only approximately symmetric. WeightedQuadrature1D's
// gather matrices (B0/B1) and scatter matrices (W00..W11) are fit
// independently (a least-squares problem each, see its own docstring) and do
// not form an exact adjoint pair; for an affine patch this is harmless (the
// pulled-back stiffness_property tensor is constant, apply() matches
// integrateStiffness() to machine precision -- see test_wq_stiffness.py),
// but for a curved patch stiffness_property varies per WQ point and sits
// between that independently-fit gather/scatter pair, breaking exact
// symmetry by a margin that shrinks under mesh refinement (same convergence
// as apply()'s own K-approximation error). Plain Conjugate Gradient can fail
// to converge on the coarsest/most-curved configurations as a result --
// matrix_free_solver.solve(..., method="gmres"/"bicgstab") tolerates it.
//
// NURBS RATIONALITY: for a rational (weighted) patch, apply() applies the
// same quotient-rule correction to the DISPLACEMENT FIELD interpolation that
// the constructor already applies to the GEOMETRY map (x_ctrl/y_ctrl -> J).
// A rational basis function R_a = w_a*N_a/W(xi) has
// dR_a/dxi = (w_a/W)*dN_a/dxi - R_a*(dW/dxi)/W -- both the trial field
// (interpolated from v) and the implicit test function (index a, realized by
// the W00..W11 scatter matrices) are rational, so the bilinear form expands
// into 4 terms (direct port of pymfiga's
// common/numerics/operations/Nurbs.py: NurbsOperations
// .compute_mf_scalar_gradu_gradv, which this codebase's WQMatrixFreeStiffness
// was found to disagree with on curved NURBS patches -- ~17% relative before
// this fix on a quarter-ring benchmark, down to matching pymfiga's own
// matrix-free result to machine precision. Confirmed via a control
// experiment: forcing all control-point weights to 1 made the two
// codebases' WQ actions agree to ~1e-16 even before this fix, isolating the
// gap to exactly this missing rational correction).
class WQMatrixFreeStiffness {
public:
    // Precomputes both directions' WeightedQuadrature1D rules and, once, the
    // pulled-back material+geometry tensor at every (tensor-product) WQ
    // point -- mirrors compute_stiffness_property(). quadtype: "1" or "2"
    // (default "2", matching pymfiga's own validated default).
    WQMatrixFreeStiffness(const Patch& patch, const ConstitutiveLaw& law,
                           const std::string& quadtype = "2");

    // K @ v without ever assembling K. v/output: length nbctrlpts_u *
    // nbctrlpts_v * 2, indexed via patch.dof_manager->get_global_dof() (the
    // same interleaved convention -- component fastest per control point --
    // as PatchIntegrator::integrateStiffness()'s assembled matrix).
    Eigen::VectorXd apply(const Eigen::VectorXd& v) const;

    int nq_u() const { return nq_u_; }
    int nq_v() const { return nq_v_; }

private:
    const Patch& patch_;
    WeightedQuadrature1D wq_u_, wq_v_;
    int nbctrlpts_u_, nbctrlpts_v_;
    int nq_u_, nq_v_;

    // stiffness_property_[I][J][alpha][beta]: length (nq_u_*nq_v_) vector
    // (u-fastest, matching matrix_free_apply_2d's convention), one per
    // (displacement component I, displacement component J, parametric test
    // direction alpha, parametric trial direction beta) combination.
    Eigen::VectorXd stiffness_property_[2][2][2][2];

    // NURBS rationality (see the class docstring's "NURBS RATIONALITY" note).
    // Only populated/used when rational_ is true.
    bool rational_ = false;
    Eigen::VectorXd w_ctrl_;             // per-control-point NURBS weight w_a
    Eigen::VectorXd w_proj_;             // 1 / W(xi), at every WQ point
    Eigen::VectorXd z_proj_u_, z_proj_v_; // (dW/dxi_dir) / W(xi)^2, per direction
};
