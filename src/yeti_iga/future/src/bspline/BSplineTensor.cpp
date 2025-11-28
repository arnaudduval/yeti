#include "BSplineTensor.hpp"

py::array_t<int> BSplineTensor::FindSpanND(const py::array_t<double>& u) const {
    auto pr = u.unchecked<1>();
    py::array_t<int> spans({(ssize_t)components.size()});
    auto sp = spans.mutable_unchecked<1>();

    FindSpanND(&pr(0), &sp(0));  // pointeur → zéro copie
    return spans;
}

void BSplineTensor::FindSpanND(const double* u, int* spans) const
{
    int dim = components.size();
    for (int d = 0; d < dim; ++d)
        spans[d] = components[d].FindSpan(u[d]);
}

void BSplineTensor::BasisFunsND_raw(const int* spans,
                                    const double* params,
                                    std::vector<ssize_t>& sizes,
                                    std::vector<double>& out) const
{
    const ssize_t n_dims = components.size();
    sizes.resize(n_dims);

    // 1) compute 1D basis functions
    // basis_vals[d] : vector containing N_i(u[d]) for dimension d
    std::vector<std::vector<double>> basis_vals(n_dims);

    for (ssize_t d = 0; d < n_dims; ++d) {
        const BSpline& B = components[d];
        int span = spans[d];
        int p = B.getDegree();

        // p-1 non zero functions
        size_t nb = p + 1;
        sizes[d] = nb;

        basis_vals[d].resize(nb);
        B.BasisFuns_raw(span, params[d], basis_vals[d].data());
    }

    // 2) build flattened ND array
    ssize_t total_size = 1;
    for (ssize_t d = 0; d < n_dims; ++d)
        total_size *= sizes[d];

    out.resize(total_size);

    std::vector<ssize_t> strides(n_dims, 1);
    for (ssize_t d = n_dims - 2; d >= 0; --d)
        strides[d] = strides[d+1] * sizes[d+1];

    // 3) tensor filling
    std::vector<ssize_t> idx(n_dims, 0);

    for (;;)
    {
        // compute linear index
        ssize_t lin = 0;
        for (ssize_t d = 0; d < n_dims; ++d)
            lin += idx[d] * strides[d];

        // compute product of basis
        double v = 1.0;
        for (ssize_t d = 0; d < n_dims; ++d)
            v *= basis_vals[d][ idx[d] ];

        out[lin] = v;

        // Increment multi-index
        ssize_t d = n_dims - 1;
        for (;;)
        {
            idx[d]++;
            if (idx[d] < sizes[d])
                break;
            idx[d] = 0;
            if (d == 0)
                return;
            d--;
        }
    }
}

/*
Computes basis functions for N dimensions
Each dimension has degree+1 non zero functions
Tensor product is computed recursively for N dimensions
result is a numpy array of shape (P0+1, p1+1, ..., pN-1 +1) containing all non zero functions
*/
py::array_t<double> BSplineTensor::BasisFunsND(const py::array_t<int>& span, const py::array_t<double>& u) const {
    ssize_t n_dims = components.size();

    if (span.ndim() != 1 || span.shape(0) != n_dims)
        throw std::runtime_error("span must be a 1D array of size n_dims");

    if (u.ndim() != 1 || u.shape(0) != n_dims)
        throw std::runtime_error("u must be a 1D array of size n_dims");

    auto sp = span.unchecked<1>();
    auto pr = u.unchecked<1>();

    // Préparation des arguments RAW
    std::vector<ssize_t> sizes;
    std::vector<double> raw_output;

    // Appel à la version optimisée
    BasisFunsND_raw(&sp(0), &pr(0), sizes, raw_output);

    // Création du tableau numpy multidimensionnel
    py::array_t<double> result(sizes);
    double* out = result.mutable_data();

    // Copie unique du buffer
    std::memcpy(out, raw_output.data(),
                raw_output.size() * sizeof(double));

    return result;
}

