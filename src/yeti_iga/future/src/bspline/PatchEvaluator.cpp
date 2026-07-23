#include "PatchEvaluator.hpp"
#include <stdexcept>
#include <cstring>

PatchEvaluator::PatchEvaluator(const Patch& patch) : patch_(patch) {
    if (!patch_.dof_manager)
        throw std::invalid_argument(
            "PatchEvaluator: patch must have a PatchDOFManager attached.");
}

py::array_t<double> PatchEvaluator::evaluateSolutionOMP(
    const py::array_t<double>& params,
    const Eigen::VectorXd& u_global) const
{
    auto par = params.unchecked<2>();
    const ssize_t n_pts  = params.shape(0);
    const ssize_t n_dims = params.shape(1);
    const int     n_comp = patch_.dof_manager->dofs_per_control_point;

    py::array_t<double> result({n_pts, (ssize_t)n_comp});
    double* res_ptr = result.mutable_data();
    std::fill(res_ptr, res_ptr + n_pts * n_comp, 0.0);

    const bool rational = patch_.cp_manager->is_rational();

    #pragma omp parallel
    {
        // Thread-local buffers to avoid repeated allocations in the hot loop
        std::vector<int>     span_buf((size_t)n_dims);
        std::vector<ssize_t> sizes;
        std::vector<double>  basis_vals;
        std::vector<ssize_t> idx((size_t)n_dims);
        std::vector<double>  val((size_t)n_comp);

        #pragma omp for schedule(static)
        for (ssize_t k = 0; k < n_pts; ++k) {
            const double* param_ptr = &par(k, 0);

            // 1. Find the knot span in each parametric direction
            patch_.tensor.FindSpanND(param_ptr, span_buf.data());
            const int* span_ptr = span_buf.data();

            // 2. Evaluate tensor-product B-spline basis (raw: u-fastest)
            patch_.tensor.BasisFunsND_raw(span_ptr, param_ptr, sizes, basis_vals);
            const ssize_t total_size = (ssize_t)basis_vals.size();

            idx.assign((size_t)n_dims, 0);
            std::fill(val.begin(), val.end(), 0.0);
            double W = 0.0;

            // 3. Accumulate weighted DOF contributions
            for (ssize_t n = 0; n < total_size; ++n) {
                // Recover the u-fastest local flat index of this active basis function
                ssize_t lin_idx = 0, stride = 1;
                for (ssize_t d = 0; d < n_dims; ++d) {
                    int p = patch_.tensor.components[d].getDegree();
                    lin_idx += ((ssize_t)(span_ptr[d] - p) + idx[d]) * stride;
                    stride  *= (ssize_t)patch_.local_shape[d];
                }

                // Basis weight: w_a * N_a for NURBS, N_a for B-spline
                double R;
                if (rational) {
                    R  = patch_.cp_manager->get_weight(patch_.global_indices[lin_idx])
                         * basis_vals[n];
                    W += R;
                } else {
                    R = basis_vals[n];
                }

                // u_h += R_a * d_a  (one entry per DOF component)
                for (int c = 0; c < n_comp; ++c)
                    val[c] += R * u_global[
                        (ssize_t)patch_.dof_manager->get_global_dof((size_t)lin_idx, c)];

                // Increment the N-D multi-index (u-fastest carry)
                for (ssize_t d = n_dims - 1; d >= 0; --d) {
                    if (++idx[d] < sizes[d]) break;
                    idx[d] = 0;
                }
            }

            // 4. Write result; NURBS: divide accumulated (w_a N_a d_a) by W
            double* dst = res_ptr + k * n_comp;
            if (rational && W > 0.0) {
                const double inv_W = 1.0 / W;
                for (int c = 0; c < n_comp; ++c)
                    dst[c] = val[c] * inv_W;
            } else {
                for (int c = 0; c < n_comp; ++c)
                    dst[c] = val[c];
            }
        }
    }

    return result;
}
