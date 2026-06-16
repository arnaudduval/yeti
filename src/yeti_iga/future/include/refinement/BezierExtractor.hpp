#pragma once
#include <vector>
#include <Eigen/Dense>
#include "BSpline.hpp"
#include "Patch.hpp"

// ND Bézier element: result of extract_nd() for one tensor-product element.
//
// The extraction operator C satisfies (u-fastest Bernstein ordering):
//   N_active(xi) = C * B_nd(xi_hat)
//
// where:
//   N_active  — the prod_d(p_d+1) active B-spline basis functions on this element
//   B_nd      — the tensor-product Bernstein polynomials, u-fastest on [0,1]^ndim
//   xi_hat_d  = (xi_d - xi_d_a) / (xi_d_b - xi_d_a)   (local coord, direction d)
//
// active: u-fastest flat indices into the patch's CP array (used to
//   extract the locally active control points before applying C).
//
// elem_index: 0-based multi-index (e_0, e_1, ...) identifying the element.
//
struct BezierElementND {
    Eigen::MatrixXd C;           // prod(p_d+1) × prod(p_d+1)
    std::vector<size_t> active;  // u-fastest flat CP indices, length = prod(p_d+1)
    std::vector<int> elem_index; // element multi-index (e_0, e_1, ...)
};


// Bézier extraction operator (Borden et al. 2011 - https://doi.org/10.1002/nme.2968).
//
// 1D API  — BezierExtractor(direction).extract(patch)
//   Returns ne matrices of shape (p+1, p+1), one per element in `direction`.
//
// ND API  — BezierExtractor::extract_nd(patch)
//   Returns one BezierElementND per tensor-product element (u-fastest ordering).
//   C = kron(C_{n-1}, kron(..., kron(C_1, C_0)...))  →  u-fastest Bézier DOFs.
//
class BezierExtractor {
public:
    explicit BezierExtractor(int direction) : direction_(direction) {}

    // 1D extraction per direction (convenience wrapper).
    std::vector<Eigen::MatrixXd> extract(const Patch& patch) const;

    // Low-level: 1D extraction operators for a single BSpline.
    static std::vector<Eigen::MatrixXd> extract_1d(const BSpline& spline);

    // Span index k (knot vector) for each element: U[k] <= xi < U[k+1].
    static std::vector<int> element_spans(const BSpline& spline);

    // ND extraction: one BezierElementND per tensor-product element.
    static std::vector<BezierElementND> extract_nd(const Patch& patch);

private:
    int direction_;

    // Kronecker product: kron(A, B).
    static Eigen::MatrixXd kron(const Eigen::MatrixXd& A,
                                const Eigen::MatrixXd& B);
};
