#pragma once
#include "refinement/RefinementOperator.hpp"
#include <vector>


class HRefiner : public RefinementOperator {
public:
    ~HRefiner() override = default;
    HRefiner(int direction, double knot) : direction_(direction), knot_(knot) {}

    // Refine patch in-place by inserting one knot (h-refinement / Boehm's algorithm).
    // Returns the full nD transition matrix (nb_new_cp × nb_old_cp).
    void refine(
        Patch& patch,
        Eigen::MatrixXd& transition_matrix
    ) const override;

    std::string getType() const override { return "h"; }

    // -----------------------------------------------------------------------
    // Low-level 1D helpers (used by SubdivisionRefiner / PRefiner fast path)
    // -----------------------------------------------------------------------

    // Build the 1D knot-insertion transition matrix for inserting `knot` into
    // `spline`.  T_1d has shape (n+1, n) where n = number of old CPs.
    // Does NOT touch the patch.
    static void compute_1d_transition(
        const BSpline& spline, double knot,
        Eigen::MatrixXd& T_1d,
        std::vector<double>& new_kv);

    // Apply a precomputed 1D transition matrix T_1d in `direction` to the
    // patch control points.  Updates cp_manager, local_shape, and
    // global_indices in-place; does NOT build the full nD matrix.
    // The caller is responsible for updating patch.tensor afterwards.
    static void apply_1d_cp_update(
        Patch& patch, int direction,
        const Eigen::MatrixXd& T_1d);

private:
    int direction_;
    double knot_;
};

// ---------------------------------------------------------------------------
// Utility: build full nD transition matrix from a 1D one (Kronecker product).
// u-fastest ordering (direction 0 fastest).
//   T_1d           — 1D matrix, shape (n_new_d, n_old_d)
//   direction      — parametric direction that was refined
//   shape_before   — local_shape of the patch before refinement
//   shape_after    — local_shape of the patch after refinement
// ---------------------------------------------------------------------------
Eigen::MatrixXd nd_transition_from_1d(
    const Eigen::MatrixXd& T_1d,
    int direction,
    const std::vector<size_t>& shape_before,
    const std::vector<size_t>& shape_after);
