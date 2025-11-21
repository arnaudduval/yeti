#include "Patch.hpp"

double* Patch::local_cp_ptr(size_t i_local) {
    size_t gid = global_indices.at(i_local);
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
    auto params = u.unchecked<2>();     // shape(n_points, dim_params)
    auto sp = spans.unchecked<2>();
    ssize_t n_points = u.shape(0);
    ssize_t n_dims = u.shape(1);
    ssize_t dim_phys = cp_manager->dim_phys;

    // output array
    py::array_t<double> result({n_points, dim_phys});
    auto res = result.mutable_unchecked<2>();

    // Get local points
    std::vector<double*> local_pts_ptrs(global_indices.size());
    for (size_t i = 0; i < global_indices.size(); ++i)
        local_pts_ptrs[i] = cp_manager->coords.data() + global_indices[i] * dim_phys;

    // Compute for each point
    for (ssize_t k = 0; k < n_points; ++k) {
        // Temporary 1D array for current point
        // TODO not optimized : there is a copy
        py::array_t<double> u_point({n_dims});
        py::array_t<int> span_point({n_dims});
        auto u_point_mut = u_point.mutable_unchecked<1>();
        auto span_point_mut = span_point.mutable_unchecked<1>();
        for (ssize_t d = 0; d < n_dims; ++d) {
            u_point_mut(d) = params(k,d);
            span_point_mut(d) = sp(k,d);
        }
        py::array_t<double> tensor_basis = tensor.BasisFunsND(span_point, u_point);

        // intialize value
        std::vector<double> val(dim_phys, 0.0);

        auto basis = tensor_basis.unchecked<2>(); // pour 2D, adaptables pour N-d
        ssize_t nu = basis.shape(0);
        ssize_t nv = basis.shape(1);

        for (ssize_t i = 0; i < nu; ++i) {
            for (ssize_t j = 0; j < nv; ++j) {
                const double* pt = local_pts_ptrs[i*nv + j];
                for (ssize_t d = 0; d < dim_phys; ++d)
                    val[d] += basis(i,j) * pt[d];
            }
        }

        for (ssize_t d = 0; d < dim_phys; ++d)
            res(k,d) = val[d];
    }

    return result;
}
