#pragma once
#include "refinement/RefinementOperator.hpp"


// Subdivide all knot spans in one parametric direction by bisection.
// Each level halves every existing span; n_levels applications double
// the number of elements per level.
class SubdivisionRefiner : public RefinementOperator {
public:
    SubdivisionRefiner(int direction, int n_levels = 1)
        : direction_(direction), n_levels_(n_levels) {}

    // Full nD transition matrix (nb_final_cp × nb_initial_cp).
    void refine(Patch& patch, Eigen::MatrixXd& transition_matrix) const override;

    // Fast path: compose only the 1D transition matrix for `direction_`.
    // T_1d has shape (n_final_1d × n_initial_1d) — much smaller than the full nD.
    // Build the full nD matrix afterwards with nd_transition_from_1d() if needed.
    void refine_1d(Patch& patch, Eigen::MatrixXd& T_1d) const;

    std::string getType() const override { return "subdivision"; }

private:
    int direction_;
    int n_levels_;
};
