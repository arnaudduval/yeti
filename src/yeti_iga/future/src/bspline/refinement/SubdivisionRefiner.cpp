#include "refinement/SubdivisionRefiner.hpp"
#include "refinement/HRefiner.hpp"


void SubdivisionRefiner::refine(
        Patch& patch,
        Eigen::MatrixXd& transition_matrix
) const {
    if (direction_ < 0 || direction_ >= static_cast<int>(patch.tensor.components.size())) {
        throw std::invalid_argument("Invalid direction for subdivision.");
    }

    size_t nb_cp = patch.global_indices.size();
    transition_matrix = Eigen::MatrixXd::Identity(nb_cp, nb_cp);

    for (int level = 0; level < n_levels_; ++level) {
        // Snapshot midpoints of all non-zero spans BEFORE any insertion this level.
        const auto& knots = patch.tensor.components[direction_].getKnotVector();
        std::vector<double> midpoints;
        for (size_t i = 0; i + 1 < knots.size(); ++i) {
            if (knots[i + 1] - knots[i] > 1e-14)
                midpoints.push_back(0.5 * (knots[i] + knots[i + 1]));
        }

        // Insert each midpoint and compose transition matrices.
        // T_local maps (pos in previous state) → (pos in new state).
        // Composition: T_total = T_last @ ... @ T_first.
        for (double mid : midpoints) {
            Eigen::MatrixXd T_local;
            HRefiner(direction_, mid).refine(patch, T_local);
            transition_matrix = T_local * transition_matrix;
        }
    }
}
