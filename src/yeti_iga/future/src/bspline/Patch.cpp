#include "Patch.hpp"
#include "SpanNDIterator.hpp"
#include "IGABasis1D.hpp"
#include "PatchIntegrator.hpp"
#include <algorithm>
#include <stdexcept>
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

        const bool rational = cp_manager->is_rational();
        double W = 0.0;

        for (ssize_t n = 0; n < total_size; ++n) {
            // u-fastest linear index: direction 0 fastest (stride = 1)
            ssize_t lin_idx = 0;
            ssize_t stride = 1;
            for (ssize_t d = 0; d < n_dims; ++d) {
                int p = tensor.components[d].getDegree();
                int global_idx = (span_ptr[d] - p) + idx[d];
                lin_idx += global_idx * stride;
                stride *= local_shape[d];
            }

            const double* pt = local_pts_ptrs[lin_idx];

            if (rational) {
                double wn = cp_manager->get_weight(global_indices[lin_idx]);
                double wN = wn * basis_vals[n];
                for (ssize_t d = 0; d < dim_phys; ++d)
                    res_ptr[k*dim_phys + d] += wN * pt[d];
                W += wN;
            } else {
                for (ssize_t d = 0; d < dim_phys; ++d)
                    res_ptr[k*dim_phys + d] += basis_vals[n] * pt[d];
            }

            // incrément de l’indice multi-dim
            for (ssize_t d = n_dims - 1; d >= 0; --d) {
                if (++idx[d] < sizes[d]) break;
                idx[d] = 0;
            }
        }

        // NURBS: divide accumulated weighted sum by W
        if (rational && W > 0.0) {
            for (ssize_t d = 0; d < dim_phys; ++d)
                res_ptr[k*dim_phys + d] /= W;
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

            const bool rational = cp_manager->is_rational();
            double W = 0.0;

            // loop over tensor-product basis
            for (ssize_t n = 0; n < total_size; ++n)
            {
                // u-fastest linear index: direction 0 fastest (stride = 1)
                ssize_t lin_idx = 0;
                ssize_t stride = 1;
                for (ssize_t d = 0; d < n_dims; ++d)
                {
                    int p = tensor.components[d].getDegree();
                    int global_idx = (span_ptr[d] - p) + idx[d];
                    lin_idx += global_idx * stride;
                    stride *= local_shape[d];
                }

                const double* pt = local_pts_ptrs[lin_idx];

                if (rational) {
                    double wN = cp_manager->get_weight(global_indices[lin_idx]) * basis_vals[n];
                    for (ssize_t d = 0; d < dim_phys; ++d)
                        val[d] += wN * pt[d];
                    W += wN;
                } else {
                    const double b = basis_vals[n];
                    for (ssize_t d = 0; d < dim_phys; ++d)
                        val[d] += b * pt[d];
                }

                // increment N-D index
                for (ssize_t d = n_dims - 1; d >= 0; --d)
                {
                    if (++idx[d] < sizes[d])
                        break;
                    idx[d] = 0;
                }
            }

            // write result (NURBS: divide by W)
            double* dst = res_ptr + kk * dim_phys;
            if (rational && W > 0.0) {
                for (ssize_t d = 0; d < dim_phys; ++d)
                    dst[d] = val[d] / W;
            } else {
                for (ssize_t d = 0; d < dim_phys; ++d)
                    dst[d] = val[d];
            }
        }

    }
    return result;
}

std::vector<const double*> Patch::control_points_for_span(const std::vector<int>& span) const {
    std::vector<const double*> pts;
    int p_u = tensor.components[0].getDegree();
    int p_v = tensor.components[1].getDegree();

    int start_u = span[0] - p_u;
    int start_v = span[1] - p_v;

    ssize_t n_u = local_shape[0];
    ssize_t n_v = local_shape[1];

    // Réserver de l'espace pour les pointeurs
    pts.reserve((p_u + 1) * (p_v + 1));

    // u-fastest: direction 0 (u) fastest
    for (int jv = 0; jv <= p_v; ++jv) {
        int lv = start_v + jv;
        for (int iu = 0; iu <= p_u; ++iu) {
            int lu = start_u + iu;
            size_t local_linear = static_cast<size_t>(lv * n_u + lu);
            pts.push_back(local_cp_ptr(local_linear));
        }
    }

    return pts;
}


std::vector<double> Patch::weights_for_span(const std::vector<int>& span) const {
    if (!cp_manager->is_rational()) return {};

    int p_u = tensor.components[0].getDegree();
    int p_v = tensor.components[1].getDegree();

    int start_u = span[0] - p_u;
    int start_v = span[1] - p_v;

    ssize_t n_u = local_shape[0];

    std::vector<double> w;
    w.reserve((p_u + 1) * (p_v + 1));
    for (int jv = 0; jv <= p_v; ++jv) {
        int lv = start_v + jv;
        for (int iu = 0; iu <= p_u; ++iu) {
            int lu = start_u + iu;
            size_t local_linear = static_cast<size_t>(lv * n_u + lu);
            w.push_back(cp_manager->get_weight(global_indices[local_linear]));
        }
    }
    return w;
}

std::vector<size_t> Patch::boundary_control_points(int direction, int side,
                                                    int span_min, int span_max) const {
    size_t ndim = tensor.components.size();
    if (ndim != 2)
        throw std::invalid_argument(
            "Patch::boundary_control_points: only 2D patches are supported (Phase 1 scope).");
    if (direction != 0 && direction != 1)
        throw std::invalid_argument("Patch::boundary_control_points: direction must be 0 or 1.");
    if (side != 0 && side != 1)
        throw std::invalid_argument("Patch::boundary_control_points: side must be 0 (min) or 1 (max).");
    if ((span_min < 0) != (span_max < 0))
        throw std::invalid_argument(
            "Patch::boundary_control_points: span_min and span_max must be given together "
            "(both >= 0), or both left at -1 for the whole edge.");

    int varying = 1 - direction;
    size_t fixed_index = (side == 0) ? 0 : local_shape[direction] - 1;

    size_t lo = 0;
    size_t hi = local_shape[varying] - 1;
    if (span_min >= 0) {
        int p_varying = tensor.components[varying].getDegree();
        int local_lo = span_min - p_varying;
        int local_hi = span_max;
        lo = static_cast<size_t>(std::max(local_lo, 0));
        hi = static_cast<size_t>(std::min(local_hi, static_cast<int>(local_shape[varying]) - 1));
    }

    std::vector<size_t> stride(ndim);
    stride[0] = 1;
    for (size_t d = 1; d < ndim; ++d)
        stride[d] = stride[d-1] * local_shape[d-1];

    std::vector<size_t> result;
    if (lo > hi) return result;
    result.reserve(hi - lo + 1);
    for (size_t i = lo; i <= hi; ++i) {
        size_t flat = fixed_index * stride[direction] + i * stride[varying];
        result.push_back(flat);
    }
    return result;
}

void Patch::Test()
{
    int ngauss_u = this->tensor.components[0].getDegree() + 1;
    int ngauss_v = this->tensor.components[1].getDegree() + 1;

    std::cout << ngauss_u << "\t" << ngauss_v << "\n";

    // build IGA basis per direction
    IGABasis1D basis_u = IGABasis1D::build(this->tensor.components[0], ngauss_u);
    IGABasis1D basis_v = IGABasis1D::build(this->tensor.components[1], ngauss_v);

    // // Assembler
    // IGAAssembler2D assembler(*this, basis_u, basis_v);
    // auto elems = assembler.assemble_stiffness();
    // std::cout << elems.size() << "\n";

    // for (auto e : elems) {
    //     std::cout << e.nb_loc << "\n";
    //     for (auto i : e.global_indices)
    //         std::cout << i << "\t";
    //     std::cout << std::endl;
    //     std::cout << e.K << std::endl;


    // }

    // return;


    // SpanNDIterator it(tensor);
    // // while(!it.is_done()) {
    // //     auto span = it.current();
    // //     std::cout << "boucle finale : " << span[0] << "\t" << span[1] << "\n";
    // //     it.next();
    // // }
    // std::cout << it.size() << "\n";

    // for (auto span : it) {
    //     auto pts = this->control_points_for_span(span);
    //     std::cout << "Span: ";
    //     for (auto s : span) std::cout << s << " ";
    //     std::cout << " --> " << pts.size() << " active points\n";
    //     for (auto pt : pts)
    //         std::cout << pt[0] << "\t" << pt[1] << "\n";

    // }

}


