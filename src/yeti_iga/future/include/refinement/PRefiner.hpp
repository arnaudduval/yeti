#pragma once
#include "refinement/RefinementOperator.hpp"


// Degree elevation along one parametric direction (Piegl & Tiller A5.9).
// Each interior knot multiplicity increases by 1 per elevation → continuity preserved.
// For maximum-continuity k-refinement: apply PRefiner BEFORE HRefiner/SubdivisionRefiner.
class PRefiner : public RefinementOperator {
public:
    explicit PRefiner(int direction, int n_elevations = 1)
        : direction_(direction), n_elevations_(n_elevations) {}

    void refine(Patch& patch, Eigen::MatrixXd& transition_matrix) const override;

    std::string getType() const override { return "p"; }

private:
    int direction_;
    int n_elevations_;

    // Apply a single elevation step in-place and return its T matrix.
    void apply_one_elevation(Patch& patch, Eigen::MatrixXd& T_step) const;

    // Build the 1-D transition matrix and elevated knot vector for one B-spline.
    static void elevate_1d(int p,
                           const std::vector<double>& U,
                           Eigen::MatrixXd& T1d,
                           std::vector<double>& new_kv);
};
