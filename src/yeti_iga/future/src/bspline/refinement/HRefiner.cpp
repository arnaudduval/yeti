#include <memory>
#include <pybind11/numpy.h>
#include "refinement/HRefiner.hpp"
#include "BSplineTensor.hpp"


std::shared_ptr<Patch> HRefiner::refine(
        const Patch& patch,
        Eigen::MatrixXd& transition_matrix
) const {
    if (direction_ < 0 || direction_ >= static_cast<int>(patch.tensor.components.size())) {
        throw std::invalid_argument("Invalid direction for knot insertion.");
    }

    const BSpline& spline = patch.tensor.components[direction_];
    const std::vector<double>& knots = spline.getKnotVector();
    int p = spline.getDegree();
    size_t dim_phys = patch.cp_manager->dim_phys;
    size_t ndim = patch.tensor.components.size();

    int k = spline.FindSpan(knot_);

    // Build new knot vector
    std::vector<double> new_knots = knots;
    auto it = std::lower_bound(new_knots.begin(), new_knots.end(), knot_);
    new_knots.insert(it, knot_);

    // Old and new CP counts
    size_t nb_old_cp = patch.global_indices.size();
    size_t nb_new_cp = 1;
    for (size_t d = 0; d < ndim; ++d) {
        nb_new_cp *= (d == static_cast<size_t>(direction_))
            ? patch.local_shape[d] + 1
            : patch.local_shape[d];
    }

    // Strides for u-fastest flat indexing:
    //   stride[0] = 1, stride[d] = product(local_shape[0..d-1])
    std::vector<size_t> old_stride(ndim), new_local_shape(patch.local_shape);
    old_stride[0] = 1;
    for (size_t d = 1; d < ndim; ++d)
        old_stride[d] = old_stride[d - 1] * patch.local_shape[d - 1];

    new_local_shape[direction_] += 1;
    std::vector<size_t> new_stride(ndim);
    new_stride[0] = 1;
    for (size_t d = 1; d < ndim; ++d)
        new_stride[d] = new_stride[d - 1] * new_local_shape[d - 1];

    // Allocate output
    std::vector<std::vector<double>> new_control_points(nb_new_cp, std::vector<double>(dim_phys, 0.0));
    transition_matrix = Eigen::MatrixXd::Zero(nb_new_cp, nb_old_cp);

    // Iterate over all lines parallel to direction_
    // (one line per combination of the other-direction indices)
    size_t nb_lines = nb_old_cp / patch.local_shape[direction_];
    std::vector<size_t> other_idx(ndim, 0);

    for (size_t line = 0; line < nb_lines; ++line) {
        // Starting flat index for this line in old and new arrays
        size_t ls_old = 0, ls_new = 0;
        for (size_t d = 0; d < ndim; ++d) {
            if (d != static_cast<size_t>(direction_)) {
                ls_old += other_idx[d] * old_stride[d];
                ls_new += other_idx[d] * new_stride[d];
            }
        }

        refine1DLine(patch, k,
                     ls_old, old_stride[direction_],
                     ls_new, new_stride[direction_],
                     new_control_points, transition_matrix);

        // Increment multi-index (carry, skipping direction_)
        for (int d = static_cast<int>(ndim) - 1; d >= 0; --d) {
            if (d == direction_) continue;
            if (++other_idx[d] < patch.local_shape[d]) break;
            other_idx[d] = 0;
        }
    }

    // Append new CPs to the shared ControlPointManager
    std::vector<size_t> new_global_indices(nb_new_cp);
    {
        std::lock_guard<std::mutex> lock(patch.cp_manager->mtx);
        size_t start_idx = patch.cp_manager->n_points();
        patch.cp_manager->coords.reserve((start_idx + nb_new_cp) * dim_phys);
        for (size_t i = 0; i < nb_new_cp; ++i) {
            for (size_t d = 0; d < dim_phys; ++d)
                patch.cp_manager->coords.push_back(new_control_points[i][d]);
            new_global_indices[i] = start_idx + i;
        }
    }

    // Build the refined patch
    std::vector<BSpline> new_components = patch.tensor.components;
    py::array_t<double> py_new_knots = py::cast(new_knots);
    new_components[direction_] = BSpline(p, py_new_knots);
    BSplineTensor new_tensor(new_components);

    auto new_patch = std::make_shared<Patch>(
        new_tensor, patch.cp_manager, new_global_indices, new_local_shape
    );

    if (patch.dof_manager) {
        new_patch->dof_manager = std::make_shared<PatchDOFManager>(
            *patch.dof_manager, new_global_indices
        );
    }

    return new_patch;
}


void HRefiner::refine1DLine(
    const Patch& patch,
    int k,
    size_t line_start_old,
    size_t line_stride_old,
    size_t line_start_new,
    size_t line_stride_new,
    std::vector<std::vector<double>>& new_control_points,
    Eigen::MatrixXd& T
) const {
    const BSpline& spline = patch.tensor.components[direction_];
    const auto& knots = spline.getKnotVector();
    int p = spline.getDegree();
    size_t dim_phys = patch.cp_manager->dim_phys;
    int n = static_cast<int>(patch.local_shape[direction_]);  // number of old CPs in this line

    // Accessors for old and new CPs via flat patch index
    auto old_cp = [&](int i) -> const double* {
        return patch.local_cp_ptr(line_start_old + static_cast<size_t>(i) * line_stride_old);
    };
    auto set_new_cp = [&](int i, const double* src) {
        size_t flat = line_start_new + static_cast<size_t>(i) * line_stride_new;
        new_control_points[flat].assign(src, src + dim_phys);
    };

    // Boehm's knot insertion (Piegl & Tiller A5.1, single insertion s=1):
    //   Q[i] = P[i]                              for 0 <= i <= k-p
    //   Q[i] = alpha[i]*P[i] + (1-alpha[i])*P[i-1]  for k-p+1 <= i <= k
    //     alpha[i] = (knot_ - U[i]) / (U[i+p] - U[i])
    //   Q[i] = P[i-1]                            for k+1 <= i <= n

    // Unchanged before insertion point: 0..k-p
    for (int i = 0; i <= k - p; ++i) {
        set_new_cp(i, old_cp(i));
        size_t nf = line_start_new + static_cast<size_t>(i) * line_stride_new;
        size_t of = line_start_old + static_cast<size_t>(i) * line_stride_old;
        T(nf, of) = 1.0;
    }

    // Blended CPs: k-p+1..k
    for (int i = k - p + 1; i <= k; ++i) {
        double denom = knots[i + p] - knots[i];
        double alpha = (denom > 1e-14) ? (knot_ - knots[i]) / denom : 0.0;

        std::vector<double> pt(dim_phys);
        const double* pi   = old_cp(i);
        const double* pi_1 = old_cp(i - 1);
        for (size_t d = 0; d < dim_phys; ++d)
            pt[d] = alpha * pi[d] + (1.0 - alpha) * pi_1[d];

        size_t nf  = line_start_new + static_cast<size_t>(i) * line_stride_new;
        size_t of_i  = line_start_old + static_cast<size_t>(i)     * line_stride_old;
        size_t of_i1 = line_start_old + static_cast<size_t>(i - 1) * line_stride_old;
        new_control_points[nf] = std::move(pt);
        T(nf, of_i)  = alpha;
        T(nf, of_i1) = 1.0 - alpha;
    }

    // Unchanged after insertion point: k+1..n  (old index i-1, new index i)
    for (int i = k + 1; i <= n; ++i) {
        set_new_cp(i, old_cp(i - 1));
        size_t nf = line_start_new + static_cast<size_t>(i)     * line_stride_new;
        size_t of = line_start_old + static_cast<size_t>(i - 1) * line_stride_old;
        T(nf, of) = 1.0;
    }
}
