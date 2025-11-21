#include "BSplineTensor.hpp"

py::array_t<int> BSplineTensor::FindSpanND(const py::array_t<double>& u) const {
    auto params = u.unchecked<1>();
    const ssize_t n_dims = params.shape(0);

    // Remove throw to improve computation time
    // if (n_dims != components.size())
    //     throw std::runtime_error("Wrong parameters dimensions, (" + std::to_string(components.size()) + " expected");

    py::array_t<int> result({n_dims});
    auto res = result.mutable_unchecked<1>();

    for (ssize_t d = 0; d < n_dims; ++d)
        res(d) = components[d].FindSpan(params(d));

    return result;
}

/*
Computes basis functions for N dimensions
Each dimension has degree+1 non zero functions
Tensor product is computed recursively for N dimensions
result is a numpy array of shape (P0+1, p1+1, ..., pN-1 +1) containing all non zero functions
*/
py::array_t<double> BSplineTensor::BasisFunsND(const py::array_t<int>& span, const py::array_t<double>& u) const {
    auto params = u.unchecked<1>();
    auto sp = span.unchecked<1>();
    const ssize_t n_dims = components.size();

    std::vector<std::vector<double>> basis_vals(n_dims);
    std::vector<ssize_t> sizes(n_dims);

    for (ssize_t d = 0; d < n_dims; ++d) {
        py::array_t<double> N = components[d].BasisFuns(sp(d), params(d));
        auto Nv = N.unchecked<1>();
        basis_vals[d].resize(Nv.shape(0));
        sizes[d] = Nv.shape(0);
        for (ssize_t i = 0; i < Nv.shape(0); ++i)
            basis_vals[d][i] = Nv(i);
    }

    py::array_t<double> result(sizes);
    double* out = result.mutable_data();

    // Precompute cumulative strides
    std::vector<ssize_t> strides(n_dims, 1);
    for (ssize_t d = n_dims - 2; d >= 0; --d)
        strides[d] = strides[d+1] * sizes[d+1];

    std::vector<ssize_t> idx(n_dims, 0);

    // Recursive filling
    std::function<void(ssize_t, double)> fill;

    fill = [&](ssize_t dim, double accum)
    {
        if (dim == n_dims) {
            // Compute linear index
            ssize_t linear = 0;
            for (ssize_t d = 0; d < n_dims; ++d)
                linear += idx[d] * strides[d];

            out[linear] = accum;
            return;
        }

        for (ssize_t i = 0; i < sizes[dim]; ++i) {
            idx[dim] = i;
            fill(dim + 1, accum * basis_vals[dim][i]);
        }
    };


    fill(0, 1.0);
    return result;
}

