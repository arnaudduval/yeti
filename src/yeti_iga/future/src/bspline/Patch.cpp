#include "Patch.hpp"
#include "SpanNDIterator.hpp"
#include <algorithm>
#include <iostream> // temp for debug

double* Patch::local_cp_ptr(size_t i_local) {
    size_t gid = global_indices[i_local];
    return cp_manager->coords.data() + gid * cp_manager->dim_phys;
}

const double* Patch::local_cp_ptr(size_t i_local) const {
    size_t gid = global_indices[i_local];
    return cp_manager->coords.data() + gid * cp_manager->dim_phys;
}

py::array_t<double> Patch::local_cp_view(size_t i_local) {
    double* ptr = local_cp_ptr(i_local);
    auto capsule = py::capsule(ptr); // do not deallocate
    return py::array_t<double>(
        {cp_manager->dim_phys},       // shape 1D
        {sizeof(double)},             // stride
        ptr,
        capsule
    );
}

py::array_t<double> Patch::local_control_point_view() const {
    std::vector<ssize_t> shape;
    for (auto s : local_shape) shape.push_back((ssize_t)s);
    shape.push_back((ssize_t)cp_manager->dim_phys);

    // return a view on the *global* coords with shape (n_points_global, dim_phys).
    std::vector<ssize_t> gshape = {(ssize_t)cp_manager->n_points(), (ssize_t)cp_manager->dim_phys};
    auto arr = py::array_t<double>(gshape, { (ssize_t)sizeof(double)*cp_manager->dim_phys, (ssize_t)sizeof(double) }, cp_manager->coords.data());
    // User can index arr[global_indices] in Python to build local array (zero-copy slicing not trivial).
    return arr;
}

py::array_t<double> Patch::EvaluatePatchND(const py::array_t<int> spans,
                                           const py::array_t<double>& u) const {
    auto params = u.unchecked<2>();   // shape(n_points, dim_params)
    auto sp = spans.unchecked<2>();
    const ssize_t n_points = u.shape(0);
    const ssize_t n_dims   = u.shape(1);
    const ssize_t dim_phys = cp_manager->dim_phys;

    // output array Python (contigu)
    py::array_t<double> result({n_points, dim_phys});
    double* res_ptr = result.mutable_data();
    // initialization
    std::memset(res_ptr, 0, n_points * dim_phys * sizeof(double));

    // pointeurs vers points de contrôle locaux
    std::vector<const double*> local_pts_ptrs(global_indices.size());
    for (size_t i = 0; i < global_indices.size(); ++i)
        local_pts_ptrs[i] = cp_manager->coords.data() + global_indices[i] * dim_phys;

    // buffer temporaire pour BasisFunsND_raw
    std::vector<ssize_t> sizes;
    std::vector<double> basis_vals;

    // évaluation pour chaque point
    for (ssize_t k = 0; k < n_points; ++k) {
        // pointeurs vers spans et params pour ce point
        const int* span_ptr   = &sp(k,0);
        const double* param_ptr = &params(k,0);

        // appel à la version raw pour remplir basis_vals
        tensor.BasisFunsND_raw(span_ptr, param_ptr, sizes, basis_vals);

        // calcul de l’indice linéaire pour chaque combinaison des dimensions
        // pour N-dimensions, on fait un produit tensoriel
        std::vector<ssize_t> idx(n_dims, 0);
        ssize_t total_size = 1;
        for (auto s : sizes) total_size *= s;

        for (ssize_t n = 0; n < total_size; ++n) {
            // linear index = j*n_u + i
            ssize_t lin_idx = 0;
            ssize_t stride = 1;
            for (ssize_t d = 0; d < n_dims; ++d) {
                // lin_idx += idx[d] * stride;
                // stride *= sizes[d];
                int p = tensor.components[d].getDegree();
                int global_idx = (span_ptr[d] - p) + idx[d];

                lin_idx += global_idx * stride;
                // stride *= sizes[d];
                stride *= local_shape[d];
            }


            // récupérer le point de contrôle correspondant
            const double* pt = local_pts_ptrs[lin_idx];

            // accumulation dans res_ptr pour ce point

            for (ssize_t d = 0; d < dim_phys; ++d)
                res_ptr[k*dim_phys + d] += basis_vals[n] * pt[d];

            // incrément de l’indice multi-dim
            for (ssize_t d = n_dims - 1; d >= 0; --d) {
                if (++idx[d] < sizes[d]) break;
                idx[d] = 0;
            }
        }
    }

    return result;
}

py::array_t<double> Patch::EvaluatePatchNDOMP(const py::array_t<int> spans,
                                              const py::array_t<double>& u) const
{
    auto params = u.unchecked<2>();   // shape(n_points, dim_params)
    auto sp = spans.unchecked<2>();
    const ssize_t n_points = u.shape(0);
    const ssize_t n_dims   = u.shape(1);
    const ssize_t dim_phys = cp_manager->dim_phys;

    // output array Python (contigu)
    py::array_t<double> result({n_points, dim_phys});
    double* res_ptr = result.mutable_data();
    // initialization
    std::memset(res_ptr, 0, (size_t)n_points * (size_t)dim_phys * sizeof(double));

    // pointeurs vers points de contrôle locaux
    std::vector<const double*> local_pts_ptrs(global_indices.size());
    for (size_t i = 0; i < global_indices.size(); ++i)
        local_pts_ptrs[i] = cp_manager->coords.data() + global_indices[i] * dim_phys;


    // -------------------------
    // Parallel loop over points
    // -------------------------
    #pragma omp parallel
    {
        // buffers thread-local (éviter réallocations fréquentes)
        std::vector<ssize_t> sizes;
        std::vector<double> basis_vals;
        std::vector<ssize_t> idx;
        std::vector<double> val;
        #pragma omp for schedule(static)
        for (ssize_t kk = 0; kk < n_points; ++kk) {
            // pointeurs vers spans et params pour ce point
            const int* span_ptr   = &sp(kk,0);
            const double* param_ptr = &params(kk,0);

            // compute tensor-product basis (RAW version)
            tensor.BasisFunsND_raw(span_ptr, param_ptr, sizes, basis_vals);
            const ssize_t total_size = basis_vals.size();


            idx.assign((size_t)n_dims, 0);
            val.assign((size_t)dim_phys, 0.0);



            // loop over tensor-product basis
            for (ssize_t n = 0; n < total_size; ++n)
            {
                // compute linear index (row-major: [d] least significant)
                ssize_t lin_idx = 0;
                ssize_t stride = 1;

                for (ssize_t d = 0; d < n_dims; ++d)
                {
                    // lin_idx += idx[d] * stride;
                    // stride *= sizes[d];
                    int p = tensor.components[d].getDegree();
                    int global_idx = (span_ptr[d] - p) + idx[d];

                    lin_idx += global_idx * stride;
                    stride *= local_shape[d];

                }

                const double b = basis_vals[n];
                const double* pt = local_pts_ptrs[lin_idx];

                for (ssize_t d = 0; d < dim_phys; ++d)
                    val[d] += b * pt[d];

                // increment N-D index
                for (ssize_t d = n_dims - 1; d >= 0; --d)
                {
                    if (++idx[d] < sizes[d])
                        break;
                    idx[d] = 0;
                }
            }

            // write result
            double* dst = res_ptr + kk * dim_phys;
            for (ssize_t d = 0; d < dim_phys; ++d)
                dst[d] = val[d];
        }

    }
    return result;
}

std::vector<const double*> Patch::control_points_for_span(const std::vector<int>& span) const {
    const ssize_t n_dims = tensor.components.size();
    const ssize_t dim_phys = cp_manager->dim_phys;

    // 1D lists of active local indices in each parametric direction
    std::vector<std::vector<size_t>> active_idx(n_dims);

    for (ssize_t d = 0; d < n_dims; ++d) {
                int deg = tensor.components[d].getDegree();
        int s = span[d];

        ssize_t start = std::max(0, s - deg);
        ssize_t end   = std::min(int(local_shape[d]-1), s);

        active_idx[d].resize(end - start + 1);
        for (ssize_t i = 0; i <= end - start; ++i)
            active_idx[d][i] = start + i;
    }

    // Build tensor-product of local indices
    std::vector<size_t> idx(n_dims, 0);
    std::vector<const double*> pts;
    ssize_t total_size = 1;
    for (auto& v : active_idx) total_size *= v.size();
    pts.reserve(total_size);

    for (ssize_t n = 0; n < total_size; ++n) {
        // compute linear index in local ordering
        ssize_t lin_idx = 0;
        ssize_t stride = 1;
    for (ssize_t d = 0; d < n_dims; ++d) {
        lin_idx += active_idx[d][idx[d]] * stride;
        stride *= local_shape[d];
    }

        pts.push_back(local_cp_ptr(lin_idx));

        // increment multi-dim index
        for (ssize_t d = n_dims - 1; d >= 0; --d) {
            if (++idx[d] < active_idx[d].size())
                break;
            idx[d] = 0;
        }
    }

    return pts;
}

void Patch::Test()
{
    SpanNDIterator it(tensor);
    // while(!it.is_done()) {
    //     auto span = it.current();
    //     std::cout << "boucle finale : " << span[0] << "\t" << span[1] << "\n";
    //     it.next();
    // }
    std::cout << it.size() << "\n";

    for (auto span : it) {
        auto pts = this->control_points_for_span(span);
        std::cout << "Span: ";
        for (auto s : span) std::cout << s << " ";
        std::cout << " --> " << pts.size() << " active points\n";
        for (auto pt : pts)
            std::cout << pt[0] << "\t" << pt[1] << "\n";

    }

}


