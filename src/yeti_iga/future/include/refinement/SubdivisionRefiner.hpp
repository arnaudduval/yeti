#pragma once
#include "refinement/RefinementOperator.hpp"


// Subdivide all knot spans in one parametric direction by bisection.
// Each level halves every existing span; n_levels applications double
// the number of elements per level.
class SubdivisionRefiner : public RefinementOperator {
public:
    SubdivisionRefiner(int direction, int n_levels = 1)
        : direction_(direction), n_levels_(n_levels) {}

    // Refine patch in-place. Returns the composed transition matrix
    // T = T_last @ ... @ T_first  (shape: nb_final_cp x nb_initial_cp).
    void refine(Patch& patch, Eigen::MatrixXd& transition_matrix) const override;

    std::string getType() const override { return "subdivision"; }

private:
    int direction_;
    int n_levels_;
};
