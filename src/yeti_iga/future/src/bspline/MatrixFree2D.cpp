#include "MatrixFree2D.hpp"
#include <stdexcept>

Eigen::VectorXd matrix_free_apply_2d(
    const Eigen::SparseMatrix<double>& Mu,
    const Eigen::SparseMatrix<double>& Mv,
    const Eigen::VectorXd& v_in,
    bool is_transpose)
{
    int nu_in = static_cast<int>(is_transpose ? Mu.rows() : Mu.cols());
    int nv_in = static_cast<int>(is_transpose ? Mv.rows() : Mv.cols());
    if (static_cast<int>(v_in.size()) != nu_in * nv_in)
        throw std::invalid_argument(
            "matrix_free_apply_2d: v_in size does not match Mu/Mv dimensions "
            "(expected nu_in * nv_in, consistent with is_transpose).");

    // V(iu, iv) = v_in(iv * nu_in + iu): Eigen's default column-major layout
    // makes this reshape a zero-copy reinterpretation, no transpose needed.
    Eigen::Map<const Eigen::MatrixXd> V(v_in.data(), nu_in, nv_in);

    Eigen::MatrixXd Result = is_transpose
        ? Eigen::MatrixXd(Mu.transpose() * V * Mv)
        : Eigen::MatrixXd(Mu * V * Mv.transpose());

    // Result is already column-major (u fastest) -- ravel is zero-copy too.
    return Eigen::Map<Eigen::VectorXd>(Result.data(), Result.size());
}
