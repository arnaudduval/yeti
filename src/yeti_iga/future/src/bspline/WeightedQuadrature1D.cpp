#include "WeightedQuadrature1D.hpp"
#include "IGABasis1D.hpp"
#include <Eigen/QR>
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace {

std::vector<double> unique_knots_01(const std::vector<double>& kv) {
    std::vector<double> out;
    for (double k : kv) {
        if (k < -1e-12 || k > 1.0 + 1e-12) continue;
        if (out.empty() || k > out.back() + 1e-12) out.push_back(k);
    }
    return out;
}

// Global (dense) basis value/derivative matrices for `bsp` at arbitrary
// points `pts` (not tied to any quadrature rule): B0(q,i)=N_i(pts[q]),
// B1(q,i)=dN_i/du(pts[q]). Shape (npts x nbctrlpts).
void eval_global_basis(const BSpline& bsp, const std::vector<double>& pts,
                        Eigen::MatrixXd& B0, Eigen::MatrixXd& B1)
{
    int p = bsp.getDegree();
    int nbctrlpts = static_cast<int>(bsp.getKnotVector().size()) - p - 1;
    int npts = static_cast<int>(pts.size());
    B0 = Eigen::MatrixXd::Zero(npts, nbctrlpts);
    B1 = Eigen::MatrixXd::Zero(npts, nbctrlpts);
    std::vector<double> ders(2 * (p + 1));
    for (int q = 0; q < npts; ++q) {
        int span = bsp.FindSpan(pts[q]);
        bsp.BasisFunsDerivatives(span, pts[q], 1, ders.data());
        for (int j = 0; j <= p; ++j) {
            int col = span - p + j;
            B0(q, col) = ders[j];
            B1(q, col) = ders[(p + 1) + j];
        }
    }
}

Eigen::MatrixXd eval_global_basis_values(const BSpline& bsp, const std::vector<double>& pts) {
    Eigen::MatrixXd B0, B1;
    eval_global_basis(bsp, pts, B0, B1);
    return B0;
}

// Direct port of pymfiga's Operations.increase_multiplicity_to_knotvector:
// every INTERIOR knot's multiplicity goes up by `repeat` (capped at degree+1).
std::vector<double> increase_multiplicity_to_knotvector(
    int repeat, int degree, const std::vector<double>& kv)
{
    if (static_cast<int>(kv.size()) < 2 * (degree + 1))
        throw std::runtime_error("WeightedQuadrature1D: knot vector too short.");

    std::vector<double> out(kv.begin(), kv.begin() + degree + 1);
    out.insert(out.end(), kv.end() - (degree + 1), kv.end());

    std::vector<double> interior_unique;
    for (int i = degree + 1; i < static_cast<int>(kv.size()) - (degree + 1); ++i) {
        double k = kv[i];
        if (interior_unique.empty() || k > interior_unique.back() + 1e-12)
            interior_unique.push_back(k);
    }
    for (double knot : interior_unique) {
        int mult = 0;
        for (double k : kv) if (std::abs(k - knot) <= 1e-12) ++mult;
        mult += repeat;
        if (mult > degree + 1)
            throw std::runtime_error(
                "WeightedQuadrature1D: knot multiplicity would exceed degree+1.");
        for (int r = 0; r < mult; ++r) out.push_back(knot);
    }
    std::sort(out.begin(), out.end());
    return out;
}

// Midpoint-rule WQ quadrature points (include_boundaries=true always, the
// only mode `future` needs). `s`/`r` are method-dependent (see build()).
// Direct port of WeightedQuadrature._midpoint_rule / _generate_quadrature_points.
std::vector<double> midpoint_rule_points(int s, int r, const std::vector<double>& unique_kv) {
    int nbelem = static_cast<int>(unique_kv.size()) - 1;

    auto linspace = [](double a, double b, int n) {
        std::vector<double> out(n);
        if (n == 1) { out[0] = a; return out; }
        for (int i = 0; i < n; ++i) out[i] = a + (b - a) * i / (n - 1);
        return out;
    };

    std::vector<double> pts;
    auto append = [&](const std::vector<double>& seg) {
        pts.insert(pts.end(), seg.begin(), seg.end());
    };
    append(linspace(unique_kv[0], unique_kv[1], r));
    append(linspace(unique_kv[nbelem - 1], unique_kv[nbelem], r));
    for (int i = 1; i < nbelem - 1; ++i)
        append(linspace(unique_kv[i], unique_kv[i + 1], 2 + s));

    std::sort(pts.begin(), pts.end());
    pts.erase(std::unique(pts.begin(), pts.end(),
                          [](double a, double b) { return std::abs(a - b) < 1e-12; }),
              pts.end());
    return pts;
}

// Per-point "knot support" width, used to regularize the least-squares
// weight computation. Direct port of
// WeightedQuadrature._compute_knot_support.
std::vector<double> compute_knot_support(const std::vector<double>& quadpts) {
    int n = static_cast<int>(quadpts.size());
    std::vector<double> extended(n + 2);
    extended[0] = -quadpts[0];
    for (int i = 0; i < n; ++i) extended[i + 1] = quadpts[i];
    extended[n + 1] = 2.0 - quadpts[n - 1];

    std::vector<double> mean(n + 1);
    for (int i = 0; i <= n; ++i) mean[i] = 0.5 * (extended[i] + extended[i + 1]);

    std::vector<double> support(n);
    for (int i = 0; i < n; ++i) support[i] = mean[i + 1] - mean[i];
    return support;
}

// Solve, for ONE basis function row: minimize ||diag(Z)^{-1} w||^2 subject to
// A @ w = b (Z mostly zero, nonzero only at quadrature points in this basis
// function's support). Reformulated (substituting w = diag(Z) @ y) as the
// minimum-norm least-squares problem min||y|| s.t. (A @ diag(Z)) @ y = b,
// solved via CompleteOrthogonalDecomposition -- like numpy.linalg.lstsq
// (rcond=None), this returns the minimum-norm solution for a rank-deficient
// system, matching pymfiga's solve_optimization_problem exactly.
Eigen::RowVectorXd solve_wq_row(const Eigen::MatrixXd& A, const Eigen::RowVectorXd& Z,
                                 const Eigen::VectorXd& b)
{
    Eigen::MatrixXd A_scaled = A * Z.asDiagonal();
    Eigen::VectorXd y = A_scaled.completeOrthogonalDecomposition().solve(b);
    return (Z.asDiagonal() * y).transpose();
}

Eigen::SparseMatrix<double> to_sparse(const Eigen::MatrixXd& dense) {
    Eigen::SparseMatrix<double> out = dense.sparseView();
    out.makeCompressed();
    return out;
}

} // namespace


WeightedQuadrature1D WeightedQuadrature1D::build(const BSpline& bspline, const std::string& quadtype) {
    if (quadtype != "1" && quadtype != "2")
        throw std::invalid_argument("WeightedQuadrature1D::build: quadtype must be \"1\" or \"2\".");

    WeightedQuadrature1D out;
    const auto& kv = bspline.getKnotVector();
    int p = bspline.getDegree();
    out.degree = p;
    out.quadtype = quadtype;
    out.nbctrlpts = static_cast<int>(kv.size()) - p - 1;

    std::vector<double> unique_kv = unique_knots_01(kv);

    // --- Test space: standard Gauss-Legendre quadrature (degree+1 points per
    // span), reusing future's existing IGABasis1D exactly as-is. Shared by
    // both methods. ---
    IGABasis1D gauss = IGABasis1D::build(bspline, p + 1, /*deriv_order=*/1);
    int nq_gauss = 0;
    for (auto& sg : gauss.gauss_spans) nq_gauss += static_cast<int>(sg.u_param.size());

    Eigen::MatrixXd B0_gauss = Eigen::MatrixXd::Zero(nq_gauss, out.nbctrlpts);
    Eigen::MatrixXd B1_gauss = Eigen::MatrixXd::Zero(nq_gauss, out.nbctrlpts);
    Eigen::VectorXd gauss_w(nq_gauss);
    std::vector<double> gauss_pts(nq_gauss);
    {
        int row = 0;
        for (const auto& kv_pair : gauss.span_indices) {
            int span = kv_pair.first;
            const auto& sg = gauss.gauss_spans[kv_pair.second];
            for (size_t q = 0; q < sg.u_param.size(); ++q, ++row) {
                for (int j = 0; j <= p; ++j) {
                    int col = span - p + j;
                    B0_gauss(row, col) = sg.N[q](j);
                    B1_gauss(row, col) = sg.dN[q](j);
                }
                gauss_w(row) = sg.weight[q];
                gauss_pts[row] = sg.u_param[q];
            }
        }
    }
    // W0cgg_test(i,q) = B0_gauss(q,i)*gauss_w(q); W1cgg_test likewise with B1.
    // Shape (nbctrlpts x nq_gauss).
    Eigen::MatrixXd W0_gauss_test =
        (B0_gauss.array().colwise() * gauss_w.array()).matrix().transpose();
    Eigen::MatrixXd W1_gauss_test =
        (B1_gauss.array().colwise() * gauss_w.array()).matrix().transpose();

    // --- WQ quadrature points: method-dependent (s, r) ---
    // Method "1": s=1, r=degree+2. Method "2": s=2, r=degree+3.
    int s = (quadtype == "1") ? 1 : 2;
    int r = (quadtype == "1") ? p + 2 : p + 3;
    out.quadpts = midpoint_rule_points(s, r, unique_kv);
    int nq_wq = static_cast<int>(out.quadpts.size());

    // B0_wq/B1_wq: the ORIGINAL (test) space basis at the WQ points -- shape
    // (nq_wq x nbctrlpts). Transposed, this is pymfiga's B0wq_test/B1wq_test.
    Eigen::MatrixXd B0_wq, B1_wq;
    eval_global_basis(bspline, out.quadpts, B0_wq, B1_wq);

    std::vector<double> support = compute_knot_support(out.quadpts);
    Eigen::RowVectorXd support_vec = Eigen::Map<const Eigen::RowVectorXd>(support.data(), nq_wq);

    // --- Target space: method-dependent ---
    // Method "2": same degree, +1 interior knot multiplicity (space S^p_{r-1}).
    // Method "1": degree-1, knot vector stripped of its first/last entry
    // (space S^{p-1}_{r-1}).
    int degree_target;
    std::vector<double> kv_target;
    if (quadtype == "2") {
        degree_target = p;
        kv_target = increase_multiplicity_to_knotvector(1, p, kv);
    } else {
        degree_target = p - 1;
        if (degree_target < 0)
            throw std::runtime_error("WeightedQuadrature1D: method \"1\" requires degree >= 1.");
        kv_target.assign(kv.begin() + 1, kv.end() - 1);
    }
    BSpline bsp_target(degree_target, kv_target);

    Eigen::MatrixXd B0_gauss_target = eval_global_basis_values(bsp_target, gauss_pts);
    // A_target (fixed across every row): shape (nbctrlpts_target x nq_wq).
    Eigen::MatrixXd A_target = eval_global_basis_values(bsp_target, out.quadpts).transpose();

    Eigen::MatrixXd W0 = Eigen::MatrixXd::Zero(out.nbctrlpts, nq_wq);
    Eigen::MatrixXd W1 = Eigen::MatrixXd::Zero(out.nbctrlpts, nq_wq);

    if (quadtype == "2") {
        // Only 2 distinct least-squares solves per basis function, both
        // constrained against the target space; reused for W00==W01, W10==W11.
        Eigen::MatrixXd integral0 = W0_gauss_test * B0_gauss_target;
        Eigen::MatrixXd integral1 = W1_gauss_test * B0_gauss_target;
        for (int i = 0; i < out.nbctrlpts; ++i) {
            Eigen::RowVectorXd Z0 = B0_wq.col(i).transpose().cwiseProduct(support_vec);
            Eigen::RowVectorXd Z1 = B1_wq.col(i).transpose().cwiseProduct(support_vec);
            W0.row(i) = solve_wq_row(A_target, Z0, integral0.row(i).transpose());
            W1.row(i) = solve_wq_row(A_target, Z1, integral1.row(i).transpose());
        }
        out.B0 = to_sparse(B0_wq);
        out.B1 = to_sparse(B1_wq);
        out.W00 = to_sparse(W0);
        out.W01 = to_sparse(W0);
        out.W10 = to_sparse(W1);
        out.W11 = to_sparse(W1);
        return out;
    }

    // Method "1": 4 distinct least-squares solves per basis function. W00/W10
    // are constrained against the ORIGINAL (test) space basis at the WQ
    // points; W01/W11 against the reduced target space.
    Eigen::MatrixXd A_test = B0_wq.transpose();  // (nbctrlpts x nq_wq)
    Eigen::MatrixXd integral_00 = W0_gauss_test * B0_gauss;         // -> W00
    Eigen::MatrixXd integral_01 = W0_gauss_test * B0_gauss_target;  // -> W01
    Eigen::MatrixXd integral_10 = W1_gauss_test * B0_gauss;         // -> W10
    Eigen::MatrixXd integral_11 = W1_gauss_test * B0_gauss_target;  // -> W11

    Eigen::MatrixXd W00 = Eigen::MatrixXd::Zero(out.nbctrlpts, nq_wq);
    Eigen::MatrixXd W01 = Eigen::MatrixXd::Zero(out.nbctrlpts, nq_wq);
    Eigen::MatrixXd W10 = Eigen::MatrixXd::Zero(out.nbctrlpts, nq_wq);
    Eigen::MatrixXd W11 = Eigen::MatrixXd::Zero(out.nbctrlpts, nq_wq);
    for (int i = 0; i < out.nbctrlpts; ++i) {
        Eigen::RowVectorXd Z0 = B0_wq.col(i).transpose().cwiseProduct(support_vec);
        Eigen::RowVectorXd Z1 = B1_wq.col(i).transpose().cwiseProduct(support_vec);
        W00.row(i) = solve_wq_row(A_test, Z0, integral_00.row(i).transpose());
        W01.row(i) = solve_wq_row(A_target, Z0, integral_01.row(i).transpose());
        W10.row(i) = solve_wq_row(A_test, Z1, integral_10.row(i).transpose());
        W11.row(i) = solve_wq_row(A_target, Z1, integral_11.row(i).transpose());
    }

    out.B0 = to_sparse(B0_wq);
    out.B1 = to_sparse(B1_wq);
    out.W00 = to_sparse(W00);
    out.W01 = to_sparse(W01);
    out.W10 = to_sparse(W10);
    out.W11 = to_sparse(W11);

    return out;
}
