#include "refinement/BezierExtractor.hpp"
#include <cmath>
#include <stdexcept>

// Algorithm 1 from Borden et al. 2011 — "Isogeometric finite element data structures
// based on Bézier extraction of NURBS".
std::vector<Eigen::MatrixXd> BezierExtractor::extract_1d(const BSpline& spline) {
    int p = spline.getDegree();
    const auto& U = spline.getKnotVector();
    int m = static_cast<int>(U.size()) - 1;  // last knot index

    // Count elements (non-zero knot spans within the parametric domain)
    int ne = 0;
    for (int k = p; k < m - p; ++k)
        if (U[k+1] - U[k] > 1e-14) ++ne;

    if (ne == 0)
        return {};

    // Each element starts as the identity
    std::vector<Eigen::MatrixXd> Ce(ne, Eigen::MatrixXd::Identity(p+1, p+1));

    int a  = p;
    int b  = p + 1;
    int el = 0;  // current element index (0-based)

    while (b < m) {
        int i = b;
        // Find multiplicity of the knot at U[b]
        while (b < m && std::abs(U[b+1] - U[b]) < 1e-14) ++b;
        int mult = b - i + 1;

        if (mult < p) {
            double numer = U[b] - U[a];
            int r = p - mult;  // number of extra insertions needed

            // alpha[j-mult-1] for j = mult+1 .. p
            std::vector<double> alphas(r);
            for (int j = mult + 1; j <= p; ++j)
                alphas[j - mult - 1] = numer / (U[a + j] - U[a]);

            for (int j = 1; j <= r; ++j) {
                int save = r - j;
                int s    = mult + j;

                // Update columns of the current element operator
                for (int k = p; k >= s; --k) {
                    double alpha = alphas[k - s];
                    Ce[el].col(k) = alpha * Ce[el].col(k) + (1.0 - alpha) * Ce[el].col(k-1);
                }

                // Transfer boundary entries to the next element's operator
                if (b < m && el + 1 < ne) {
                    for (int row = save; row <= save + j; ++row)
                        Ce[el+1](row, save) = Ce[el](p - j + (row - save), p);
                }
            }
        }

        ++el;
        if (el >= ne) break;
        a = b;
        ++b;
    }

    return Ce;
}

std::vector<int> BezierExtractor::element_spans(const BSpline& spline) {
    int p = spline.getDegree();
    const auto& U = spline.getKnotVector();
    int m = static_cast<int>(U.size()) - 1;
    std::vector<int> spans;
    for (int k = p; k < m - p; ++k)
        if (U[k+1] - U[k] > 1e-14)
            spans.push_back(k);
    return spans;
}

std::vector<Eigen::MatrixXd> BezierExtractor::extract(const Patch& patch) const {
    if (direction_ < 0 ||
        direction_ >= static_cast<int>(patch.tensor.components.size()))
        throw std::invalid_argument("BezierExtractor: invalid direction");
    return extract_1d(patch.tensor.components[direction_]);
}

Eigen::MatrixXd BezierExtractor::kron(const Eigen::MatrixXd& A,
                                      const Eigen::MatrixXd& B) {
    Eigen::MatrixXd result(A.rows() * B.rows(), A.cols() * B.cols());
    for (int i = 0; i < A.rows(); ++i)
        for (int j = 0; j < A.cols(); ++j)
            result.block(i * B.rows(), j * B.cols(), B.rows(), B.cols())
                = A(i, j) * B;
    return result;
}

std::vector<BezierElementND> BezierExtractor::extract_nd(const Patch& patch) {
    size_t ndim = patch.tensor.components.size();

    // 1D extraction operators and span indices per direction
    std::vector<std::vector<Eigen::MatrixXd>> Ce_d(ndim);
    std::vector<std::vector<int>>             sp_d(ndim);
    std::vector<int>                          ne_d(ndim);
    for (size_t d = 0; d < ndim; ++d) {
        Ce_d[d] = extract_1d(patch.tensor.components[d]);
        sp_d[d] = element_spans(patch.tensor.components[d]);
        ne_d[d] = static_cast<int>(Ce_d[d].size());
    }

    // u-fastest strides for global CP flat index
    std::vector<size_t> cp_stride(ndim);
    cp_stride[0] = 1;
    for (size_t d = 1; d < ndim; ++d)
        cp_stride[d] = cp_stride[d-1] * patch.local_shape[d-1];

    // Total number of ND elements
    size_t ne_total = 1;
    for (size_t d = 0; d < ndim; ++d) ne_total *= static_cast<size_t>(ne_d[d]);

    std::vector<BezierElementND> elements;
    elements.reserve(ne_total);

    // Iterate over all ND elements in u-fastest order (direction 0 fastest)
    std::vector<int> em(ndim, 0);  // element multi-index

    for (size_t e = 0; e < ne_total; ++e) {

        // --- Kronecker product: C = kron(C_{n-1}, kron(..., kron(C_1, C_0)...)) ---
        // Gives u-fastest Bézier DOF ordering: flat = j_0 + j_1*(p_0+1) + ...
        Eigen::MatrixXd C = Ce_d[0][em[0]];
        for (size_t d = 1; d < ndim; ++d)
            C = kron(Ce_d[d][em[d]], C);

        // --- Active B-spline CP flat indices (u-fastest), length = prod(p_d+1) ---
        std::vector<int> p_d(ndim);
        std::vector<int> n_loc(ndim);   // p_d + 1 per direction
        size_t n_local = 1;
        for (size_t d = 0; d < ndim; ++d) {
            p_d[d]  = patch.tensor.components[d].getDegree();
            n_loc[d] = p_d[d] + 1;
            n_local *= static_cast<size_t>(n_loc[d]);
        }

        std::vector<size_t> active;
        active.reserve(n_local);
        std::vector<int> li(ndim, 0);  // local multi-index within the element

        for (size_t l = 0; l < n_local; ++l) {
            size_t flat = 0;
            for (size_t d = 0; d < ndim; ++d) {
                int global_i = (sp_d[d][em[d]] - p_d[d]) + li[d];
                flat += static_cast<size_t>(global_i) * cp_stride[d];
            }
            active.push_back(flat);

            // Increment li in u-fastest order
            for (size_t d = 0; d < ndim; ++d) {
                if (++li[d] < n_loc[d]) break;
                li[d] = 0;
            }
        }

        elements.push_back({std::move(C), std::move(active), em});

        // Increment element multi-index in u-fastest order
        for (size_t d = 0; d < ndim; ++d) {
            if (++em[d] < ne_d[d]) break;
            em[d] = 0;
        }
    }

    return elements;
}
