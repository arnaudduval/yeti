#include <memory>
#include <algorithm>
#include <pybind11/numpy.h>
#include "refinement/HRefiner.hpp"
#include "BSplineTensor.hpp"


void HRefiner::refine(
        Patch& patch,
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

    size_t nb_old_cp = patch.global_indices.size();
    size_t nb_new_cp = 1;
    for (size_t d = 0; d < ndim; ++d) {
        nb_new_cp *= (d == static_cast<size_t>(direction_))
            ? patch.local_shape[d] + 1
            : patch.local_shape[d];
    }

    // u-fastest strides: direction 0 fastest (stride[0]=1)
    std::vector<size_t> old_stride(ndim), new_local_shape(patch.local_shape);
    new_local_shape[direction_] += 1;

    old_stride[0] = 1;
    for (size_t d = 1; d < ndim; ++d)
        old_stride[d] = old_stride[d-1] * patch.local_shape[d-1];

    std::vector<size_t> new_stride(ndim);
    new_stride[0] = 1;
    for (size_t d = 1; d < ndim; ++d)
        new_stride[d] = new_stride[d-1] * new_local_shape[d-1];

    size_t nb_lines = nb_old_cp / patch.local_shape[direction_];
    int n = static_cast<int>(patch.local_shape[direction_]);

    transition_matrix = Eigen::MatrixXd::Zero(nb_new_cp, nb_old_cp);
    std::vector<size_t> new_global_indices(nb_new_cp);

    // Blended CP coordinates (only p * nb_lines new points, not a full copy)
    size_t nb_blended = static_cast<size_t>(p) * nb_lines;
    std::vector<std::vector<double>> blended_coords;
    blended_coords.reserve(nb_blended);

    // --- First pass: T matrix + reuse/collect indices ---
    std::vector<size_t> other_idx(ndim, 0);

    for (size_t line = 0; line < nb_lines; ++line) {
        size_t ls_old = 0, ls_new = 0;
        for (size_t d = 0; d < ndim; ++d) {
            if (d != static_cast<size_t>(direction_)) {
                ls_old += other_idx[d] * old_stride[d];
                ls_new += other_idx[d] * new_stride[d];
            }
        }

        // Unchanged before insertion point: 0..k-p
        for (int i = 0; i <= k - p; ++i) {
            size_t nf = ls_new + static_cast<size_t>(i) * new_stride[direction_];
            size_t of = ls_old + static_cast<size_t>(i) * old_stride[direction_];
            new_global_indices[nf] = patch.global_indices[of];
            transition_matrix(nf, of) = 1.0;
        }

        // Blended CPs: k-p+1..k  (truly new — computed from neighbours)
        for (int i = k - p + 1; i <= k; ++i) {
            double denom = knots[i + p] - knots[i];
            double alpha = (denom > 1e-14) ? (knot_ - knots[i]) / denom : 0.0;

            const double* pi   = patch.local_cp_ptr(ls_old + static_cast<size_t>(i)     * old_stride[direction_]);
            const double* pi_1 = patch.local_cp_ptr(ls_old + static_cast<size_t>(i - 1) * old_stride[direction_]);

            std::vector<double> pt(dim_phys);
            for (size_t d = 0; d < dim_phys; ++d)
                pt[d] = alpha * pi[d] + (1.0 - alpha) * pi_1[d];
            blended_coords.push_back(std::move(pt));

            size_t nf    = ls_new + static_cast<size_t>(i)     * new_stride[direction_];
            size_t of_i  = ls_old + static_cast<size_t>(i)     * old_stride[direction_];
            size_t of_i1 = ls_old + static_cast<size_t>(i - 1) * old_stride[direction_];
            transition_matrix(nf, of_i)  = alpha;
            transition_matrix(nf, of_i1) = 1.0 - alpha;
            // new_global_indices[nf] filled in second pass
        }

        // Unchanged after: k+1..n  (old index i-1, new index i)
        for (int i = k + 1; i <= n; ++i) {
            size_t nf = ls_new + static_cast<size_t>(i)     * new_stride[direction_];
            size_t of = ls_old + static_cast<size_t>(i - 1) * old_stride[direction_];
            new_global_indices[nf] = patch.global_indices[of];
            transition_matrix(nf, of) = 1.0;
        }

        // Carry increment over non-direction_ indices
        for (int d = static_cast<int>(ndim) - 1; d >= 0; --d) {
            if (d == direction_) continue;
            if (++other_idx[d] < patch.local_shape[d]) break;
            other_idx[d] = 0;
        }
    }

    // --- Second pass: append blended CPs to cp_manager, fill their global indices ---
    {
        std::lock_guard<std::mutex> lock(patch.cp_manager->mtx);
        size_t start_idx = patch.cp_manager->n_points();
        patch.cp_manager->coords.reserve((start_idx + nb_blended) * dim_phys);

        size_t blend_counter = 0;
        std::fill(other_idx.begin(), other_idx.end(), 0);

        for (size_t line = 0; line < nb_lines; ++line) {
            size_t ls_new = 0;
            for (size_t d = 0; d < ndim; ++d) {
                if (d != static_cast<size_t>(direction_))
                    ls_new += other_idx[d] * new_stride[d];
            }

            for (int i = k - p + 1; i <= k; ++i) {
                size_t nf = ls_new + static_cast<size_t>(i) * new_stride[direction_];
                for (size_t d = 0; d < dim_phys; ++d)
                    patch.cp_manager->coords.push_back(blended_coords[blend_counter][d]);
                new_global_indices[nf] = start_idx + blend_counter;
                ++blend_counter;
            }

            for (int d = static_cast<int>(ndim) - 1; d >= 0; --d) {
                if (d == direction_) continue;
                if (++other_idx[d] < patch.local_shape[d]) break;
                other_idx[d] = 0;
            }
        }
    }

    // --- Update patch in-place ---
    patch.global_indices = std::move(new_global_indices);
    patch.local_shape    = std::move(new_local_shape);

    std::vector<BSpline> new_components = patch.tensor.components;
    py::array_t<double> py_new_knots = py::cast(new_knots);
    new_components[direction_] = BSpline(p, py_new_knots);
    patch.tensor = BSplineTensor(new_components);

    if (patch.dof_manager) {
        patch.dof_manager = std::make_shared<PatchDOFManager>(
            *patch.dof_manager, patch.global_indices
        );
    }

    // --- Compact cp_manager: rebuild coords in u-fastest flat order ---
    // This ensures mgr.coords_view() always reflects the u-fastest ordering
    // of the patch (mgr[i] == patch.control_point(i) after compaction).
    {
        std::lock_guard<std::mutex> lock(patch.cp_manager->mtx);
        size_t dim = patch.cp_manager->dim_phys;

        // Build compact coords in flat (u-fastest) order
        std::vector<double> new_coords;
        new_coords.reserve(nb_new_cp * dim);
        for (size_t flat = 0; flat < nb_new_cp; ++flat) {
            size_t gid = patch.global_indices[flat];
            const double* src = patch.cp_manager->coords.data() + gid * dim;
            new_coords.insert(new_coords.end(), src, src + dim);
        }

        // Sequential mapping: mgr[i] == flat i
        std::iota(patch.global_indices.begin(), patch.global_indices.end(), 0);

        patch.cp_manager->coords = std::move(new_coords);
    }
}
