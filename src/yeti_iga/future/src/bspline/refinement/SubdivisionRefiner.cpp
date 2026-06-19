#include "refinement/SubdivisionRefiner.hpp"
#include "refinement/HRefiner.hpp"
#include "BSplineTensor.hpp"
#include <pybind11/numpy.h>


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


void SubdivisionRefiner::refine_1d(Patch& patch, Eigen::MatrixXd& T_1d,
                                    const std::unordered_set<size_t>& protected_global_ids) const {
    if (direction_ < 0 || direction_ >= static_cast<int>(patch.tensor.components.size()))
        throw std::invalid_argument("Invalid direction for subdivision.");

    // Track the evolving knot vector locally — no CP update per step.
    BSpline spline_1d = patch.tensor.components[direction_];
    int n = static_cast<int>(patch.local_shape[direction_]);
    T_1d = Eigen::MatrixXd::Identity(n, n);

    for (int level = 0; level < n_levels_; ++level) {
        const auto& kv = spline_1d.getKnotVector();
        std::vector<double> midpoints;
        for (size_t i = 0; i + 1 < kv.size(); ++i)
            if (kv[i + 1] - kv[i] > 1e-14)
                midpoints.push_back(0.5 * (kv[i] + kv[i + 1]));

        for (double mid : midpoints) {
            Eigen::MatrixXd T_step;
            std::vector<double> new_kv;
            HRefiner::compute_1d_transition(spline_1d, mid, T_step, new_kv);
            spline_1d = BSpline(spline_1d.getDegree(), py::cast(new_kv));
            T_1d = T_step * T_1d;
        }
    }

    // Apply composed T_1d to patch CPs once (instead of once per insertion).
    HRefiner::apply_1d_cp_update(patch, direction_, T_1d, protected_global_ids);

    std::vector<BSpline> new_components = patch.tensor.components;
    new_components[direction_] = spline_1d;
    patch.tensor = BSplineTensor(new_components);
}
