#pragma once
#include <Eigen/Dense>
#include <Eigen/Sparse>

// Matrix-free application of a 2D tensor-product (Kronecker) operator to a
// vector, without ever forming the Kronecker product matrix. Direct 2D port
// of pymfiga's MatrixFree.apply
// (src/yeti_iga/pymfiga/common/numerics/operations/matrix_free.py), which is
// generic over an arbitrary number of modes via tensor reshaping; `future`
// is 2D parametric, so this specializes to exactly 2 modes, which reduces to
// a plain "sandwich" of two matrix products -- no ND tensor machinery
// needed.
//
// Convention: `v_in` is interpreted as a (nu_in x nv_in) matrix with u
// varying fastest -- i.e. v_in[iv * nu_in + iu] -- matching `future`'s own
// control-point/DOF ordering convention used throughout (BSplineSurface,
// ControlPointManager, PatchDOFManager, ...). Mu acts along the u direction,
// Mv along v. Mathematically:
//   is_transpose = false: result = ravel(Mu @ V @ Mv^T), i.e. (Mv (x) Mu) @ v_in
//   is_transpose = true:  result = ravel(Mu^T @ V @ Mv), i.e. (Mv (x) Mu)^T @ v_in
// where V is v_in reshaped to (nu_in x nv_in) and (x) is the Kronecker
// product -- never actually formed here. `is_transpose` lets the same pair
// of matrices serve both as a "gather" (parametric-space interpolation,
// e.g. WeightedQuadrature1D's B0/B1 mapping control points to quadrature
// points) and, via the flag, a "scatter" back (e.g. its W matrices mapping
// quadrature-point contributions back to control points), matching
// pymfiga's own usage pattern.
Eigen::VectorXd matrix_free_apply_2d(
    const Eigen::SparseMatrix<double>& Mu,
    const Eigen::SparseMatrix<double>& Mv,
    const Eigen::VectorXd& v_in,
    bool is_transpose = false);
