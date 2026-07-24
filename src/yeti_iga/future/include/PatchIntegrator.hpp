#pragma once
#include <algorithm>
#include <functional>
#include <memory>
#include <stdexcept>
#include <vector>
#include <cstddef>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "IGABasis1D.hpp"
#include "Patch.hpp"
#include "PatchAssembly.hpp"
#include "LocalOperator.hpp"
#include "ScalarLocalOperator.hpp"
#include "Traction.hpp"
#include "ConstitutiveLaw.hpp"



// One distributed boundary load to assemble over a whole PatchAssembly (see
// PatchIntegrator::assembleBoundaryLoad()). Unlike `laws`/`operators`
// (one entry per patch), a boundary load only applies to specific
// patches/edges, so specs are given as a sparse list instead.
struct BoundaryLoadSpec {
    size_t patch_index;
    int direction, side;
    std::shared_ptr<Traction> traction;
    int span_min = -1, span_max = -1;
};

class PatchIntegrator {
private:
    const Patch& patch_;
    // TODO: generalize it for N dimensions
    const IGABasis1D& basis_u_;
    const IGABasis1D& basis_v_;
    std::shared_ptr<const ConstitutiveLaw> law_;

    // No-law constructor for assembleOperator() / assembleBoundaryLoad() which
    // never call the stiffness/mass kernels. Static methods can access it.
    PatchIntegrator(const Patch& patch, const IGABasis1D& basis_u, const IGABasis1D& basis_v)
        : patch_(patch), basis_u_(basis_u), basis_v_(basis_v), law_(nullptr) {}

    // Geometry shared by the stiffness and mass kernels at one Gauss point:
    // basis values/derivatives in physical-mapping order (u-fastest), the
    // Jacobian components and its determinant. Mass only needs R and detJ;
    // stiffness derives invJ/grads from J11..J22 and dRdu/dRdv on top of
    // this. detJ == 0.0 signals a degenerate point to be skipped by the caller.
    struct GaussPointGeometry {
        std::vector<double> R, dRdu, dRdv;
        double J11, J12, J21, J22;
        double detJ;
    };
    GaussPointGeometry evaluateGaussPointGeometry(
        const std::vector<const double*>& pts,
        const Eigen::VectorXd& Nu, const Eigen::VectorXd& dNu,
        const Eigen::VectorXd& Nv, const Eigen::VectorXd& dNv) const;

    // NURBS rational variant: same as evaluateGaussPointGeometry but applies
    // the quotient rule (R_a = w_a N_a / W) AFTER computing the raw tensor
    // product and BEFORE computing the Jacobian, so J/detJ are already those
    // of the NURBS surface. Weights are the active CPs' weights for this
    // span, in the same u-fastest order as pts.
    GaussPointGeometry evaluateGaussPointGeometryNURBS(
        const std::vector<const double*>& pts,
        const Eigen::VectorXd& Nu, const Eigen::VectorXd& dNu,
        const Eigen::VectorXd& Nv, const Eigen::VectorXd& dNv,
        const std::vector<double>& weights) const;

    // Compute local stiffness contribution for a given span (B-free kernel).
    Eigen::MatrixXd computeLocalStiffnessContribution(const Patch& patch, const SpanGauss1D& sg_u, const SpanGauss1D& sg_v, const std::vector<int>& span);
    // NURBS variant -- same but uses evaluateGaussPointGeometryNURBS
    Eigen::MatrixXd computeLocalStiffnessContributionNURBS(const Patch& patch, const SpanGauss1D& sg_u, const SpanGauss1D& sg_v, const std::vector<int>& span, const std::vector<double>& weights);

    // Compute local mass contribution for a given span
    Eigen::MatrixXd computeLocalMassContribution(const Patch& patch, const SpanGauss1D& sg_u, const SpanGauss1D& sg_v, const std::vector<int>& span);
    // NURBS variant
    Eigen::MatrixXd computeLocalMassContributionNURBS(const Patch& patch, const SpanGauss1D& sg_u, const SpanGauss1D& sg_v, const std::vector<int>& span, const std::vector<double>& weights);

    // Assemble local contribution into global matrix (shared by stiffness and mass)
    void assembleLocalContribution(const Eigen::MatrixXd& local_contribution, const std::vector<int>& span, std::vector<Eigen::Triplet<double>>& tripletList);

    // Compute local boundary-load contribution for a given boundary span (a
    // regular ND span whose `direction` component is fixed at the edge).
    // sg_varying is the SpanGauss1D of the OTHER ("varying") direction for
    // this span; N_fixed_boundary is the fixed direction's basis row
    // evaluated once at the boundary parameter (no Gauss loop needed there).
    Eigen::VectorXd computeLocalBoundaryLoadContribution(
        const Patch& patch, const std::vector<int>& span, int direction,
        const SpanGauss1D& sg_varying, const Eigen::VectorXd& N_fixed_boundary,
        const Traction& traction);
    // NURBS variant
    Eigen::VectorXd computeLocalBoundaryLoadContributionNURBS(
        const Patch& patch, const std::vector<int>& span, int direction,
        const SpanGauss1D& sg_varying, const Eigen::VectorXd& N_fixed_boundary,
        const Traction& traction, const std::vector<double>& weights);

    // if constexpr dispatch helpers: IsRational=false generates code identical
    // to the current B-spline-only path; IsRational=true calls the NURBS
    // variants. The dispatch is done once per collect*() call (per patch per
    // assembly), so there is zero overhead inside the Gauss loop for B-splines.
    template<bool IsRational>
    void collectTripletsImpl(std::vector<Eigen::Triplet<double>>& tripletList);
    template<bool IsRational>
    void collectMassTripletsImpl(std::vector<Eigen::Triplet<double>>& tripletList);
    template<bool IsRational>
    void collectOperatorTripletsImpl(LocalOperator& op, std::vector<Eigen::Triplet<double>>& tripletList);

    // Same as assembleLocalContribution(), but scatters a full local load
    // vector (not a sparse matrix contribution) directly into a dense global
    // vector using global dof indices.
    void assembleLocalLoadContribution(
        const Eigen::VectorXd& local_load, const std::vector<int>& span, Eigen::VectorXd& global_load);

    // Patch-LOCAL flat positions (u-fastest) of the span's active control points.
    std::vector<size_t> buildSpanLocalIndices(const Patch& patch, const std::vector<int>& span) const;

    // Shared orchestration for assembleStiffness()/assembleMass(): validates
    // laws, builds one IGABasis1D pair + PatchIntegrator per patch, invokes
    // `collect` into one shared triplet list, sizes the result to the
    // assembly-wide total dof count.
    static Eigen::SparseMatrix<double> assembleGeneric(
        const PatchAssembly& assembly,
        const std::vector<std::shared_ptr<const ConstitutiveLaw>>& laws,
        int gauss_n,
        const std::function<void(PatchIntegrator&, std::vector<Eigen::Triplet<double>>&)>& collect);

public:
    PatchIntegrator(const Patch& patch, const IGABasis1D& basis_u, const IGABasis1D& basis_v,
                    std::shared_ptr<const ConstitutiveLaw> law)
        : patch_(patch), basis_u_(basis_u), basis_v_(basis_v), law_(std::move(law)) {}

    // Append this patch's contribution to a (possibly shared) triplet list,
    // using GLOBAL dof indices. Does not size or build any matrix -- this is
    // what lets assembleStiffness() merge several patches' contributions
    // before a single setFromTriplets() call. Dispatches to
    // collectTripletsImpl<IsRational> once (outside the span loop), so the
    // Gauss loop is zero-overhead for B-splines.
    void collectTriplets(std::vector<Eigen::Triplet<double>>& tripletList) {
        if (!law_)
            throw std::logic_error(
                "PatchIntegrator::collectTriplets: no ConstitutiveLaw set.");
        if (patch_.cp_manager->is_rational())
            collectTripletsImpl<true>(tripletList);
        else
            collectTripletsImpl<false>(tripletList);
    }

    // Same as collectTriplets(), but for the consistent mass matrix
    // contribution (requires law->material().rho > 0).
    void collectMassTriplets(std::vector<Eigen::Triplet<double>>& tripletList) {
        if (!law_)
            throw std::logic_error(
                "PatchIntegrator::collectMassTriplets: no ConstitutiveLaw set.");
        if (law_->material().rho <= 0.0)
            throw std::invalid_argument(
                "PatchIntegrator::collectMassTriplets: material.rho "
                "must be set (> 0) to assemble a mass matrix.");
        if (patch_.cp_manager->is_rational())
            collectMassTripletsImpl<true>(tripletList);
        else
            collectMassTripletsImpl<false>(tripletList);
    }

    // This patch's own highest global dof + 1 (scanning every local control
    // point, not just the last one -- global dof order need not follow
    // local position order once global_indices isn't the identity mapping,
    // e.g. for any patch that isn't the first one added to a PatchAssembly).
    size_t localTotalDofs() const {
        const auto& local_to_global = patch_.dof_manager->local_to_global_dofs;
        return *std::max_element(local_to_global.begin(), local_to_global.end()) + 1;
    }

    Eigen::SparseMatrix<double> integrateStiffness() {
        std::vector<Eigen::Triplet<double>> tripletList;
        collectTriplets(tripletList);

        size_t total_dofs = localTotalDofs();

        Eigen::SparseMatrix<double> stiffness_matrix(total_dofs, total_dofs);
        stiffness_matrix.setFromTriplets(tripletList.begin(), tripletList.end());

        return stiffness_matrix;
    }

    // Same as integrateStiffness(), but for the consistent mass matrix
    // (requires law->material().rho > 0).
    Eigen::SparseMatrix<double> integrateMass() {
        std::vector<Eigen::Triplet<double>> tripletList;
        collectMassTriplets(tripletList);

        size_t total_dofs = localTotalDofs();

        Eigen::SparseMatrix<double> mass_matrix(total_dofs, total_dofs);
        mass_matrix.setFromTriplets(tripletList.begin(), tripletList.end());

        return mass_matrix;
    }

    // Assemble the global stiffness matrix of a whole PatchAssembly.
    // laws must have one entry per patch, in assembly.get_patchs() (add_patch()) order.
    // gauss_n: Gauss points per span per direction; 0 (default) means degree+1 per direction.
    static Eigen::SparseMatrix<double> assembleStiffness(
        const PatchAssembly& assembly,
        const std::vector<std::shared_ptr<const ConstitutiveLaw>>& laws,
        int gauss_n = 0);

    // Same as assembleStiffness(), but for the consistent mass matrix.
    // Every entry in laws must have material().rho > 0.
    static Eigen::SparseMatrix<double> assembleMass(
        const PatchAssembly& assembly,
        const std::vector<std::shared_ptr<const ConstitutiveLaw>>& laws,
        int gauss_n = 0);

    // Same as collectTriplets()/collectMassTriplets(), but the per-Gauss-
    // point term comes from a caller-supplied LocalOperator instead of a
    // built-in kernel. Dispatches to collectOperatorTripletsImpl<IsRational>
    // so the Gauss loop uses the correct (rational or B-spline) geometry.
    void collectOperatorTriplets(LocalOperator& op, std::vector<Eigen::Triplet<double>>& tripletList) {
        if (patch_.cp_manager->is_rational())
            collectOperatorTripletsImpl<true>(op, tripletList);
        else
            collectOperatorTripletsImpl<false>(op, tripletList);
    }

    // Same as integrateStiffness()/integrateMass(), but for a custom
    // LocalOperator over this single patch.
    Eigen::SparseMatrix<double> integrateOperator(LocalOperator& op) {
        std::vector<Eigen::Triplet<double>> tripletList;
        collectOperatorTriplets(op, tripletList);

        size_t total_dofs = localTotalDofs();

        Eigen::SparseMatrix<double> result(total_dofs, total_dofs);
        result.setFromTriplets(tripletList.begin(), tripletList.end());

        return result;
    }

    // Assemble a custom LocalOperator over a whole PatchAssembly. operators
    // must have one entry per patch, in assembly.get_patchs() (add_patch())
    // order.
    static Eigen::SparseMatrix<double> assembleOperator(
        const PatchAssembly& assembly,
        const std::vector<std::shared_ptr<LocalOperator>>& operators,
        int gauss_n = 0);

    // Integrate a distributed boundary load (force per unit length) over
    // this single patch's edge.
    Eigen::VectorXd integrateBoundaryLoad(
        int direction, int side, const Traction& traction,
        int span_min = -1, int span_max = -1);

    // Assemble several distributed boundary loads over a whole PatchAssembly.
    static Eigen::VectorXd assembleBoundaryLoad(
        const PatchAssembly& assembly,
        const std::vector<BoundaryLoadSpec>& specs,
        int gauss_n = 0);

    // Integrate a ScalarLocalOperator over this single patch.
    double integrateScalarOperator(ScalarLocalOperator& op,
                                   const Eigen::VectorXd& u_global);

};

// ─────────────────────────────────────────────────────────────────────────────
// Template bodies (must be visible at instantiation; kept here to avoid an
// explicit-instantiation .cpp).  Each Impl<false> generates code identical to
// the old B-spline-only methods; Impl<true> calls the NURBS variants.
// ─────────────────────────────────────────────────────────────────────────────

template<bool IsRational>
void PatchIntegrator::collectTripletsImpl(std::vector<Eigen::Triplet<double>>& tripletList) {
    SpanNDIterator it = patch_.spans();
    for (auto span : it) {
        int idx_u = basis_u_.span_indices.at(span[0]);
        int idx_v = basis_v_.span_indices.at(span[1]);
        const SpanGauss1D& sg_u = basis_u_.gauss_spans[idx_u];
        const SpanGauss1D& sg_v = basis_v_.gauss_spans[idx_v];

        Eigen::MatrixXd lc;
        if constexpr (IsRational) {
            auto w = patch_.weights_for_span(span);
            lc = computeLocalStiffnessContributionNURBS(patch_, sg_u, sg_v, span, w);
        } else {
            lc = computeLocalStiffnessContribution(patch_, sg_u, sg_v, span);
        }
        assembleLocalContribution(lc, span, tripletList);
    }
}

template<bool IsRational>
void PatchIntegrator::collectMassTripletsImpl(std::vector<Eigen::Triplet<double>>& tripletList) {
    SpanNDIterator it = patch_.spans();
    for (auto span : it) {
        int idx_u = basis_u_.span_indices.at(span[0]);
        int idx_v = basis_v_.span_indices.at(span[1]);
        const SpanGauss1D& sg_u = basis_u_.gauss_spans[idx_u];
        const SpanGauss1D& sg_v = basis_v_.gauss_spans[idx_v];

        Eigen::MatrixXd lc;
        if constexpr (IsRational) {
            auto w = patch_.weights_for_span(span);
            lc = computeLocalMassContributionNURBS(patch_, sg_u, sg_v, span, w);
        } else {
            lc = computeLocalMassContribution(patch_, sg_u, sg_v, span);
        }
        assembleLocalContribution(lc, span, tripletList);
    }
}

template<bool IsRational>
void PatchIntegrator::collectOperatorTripletsImpl(
    LocalOperator& op, std::vector<Eigen::Triplet<double>>& tripletList)
{
    SpanNDIterator it = patch_.spans();
    for (auto span : it) {
        int idx_u = basis_u_.span_indices.at(span[0]);
        int idx_v = basis_v_.span_indices.at(span[1]);
        const SpanGauss1D& sg_u = basis_u_.gauss_spans[idx_u];
        const SpanGauss1D& sg_v = basis_v_.gauss_spans[idx_v];

        std::vector<const double*> pts = patch_.control_points_for_span(span);
        size_t nb_loc = pts.size();

        [[maybe_unused]] std::vector<double> w;
        if constexpr (IsRational) w = patch_.weights_for_span(span);

        int ngauss_u = static_cast<int>(sg_u.u_param.size());
        int ngauss_v = static_cast<int>(sg_v.u_param.size());

        Eigen::MatrixXd local_contribution;

        for (int gu = 0; gu < ngauss_u; ++gu) {
            for (int gv = 0; gv < ngauss_v; ++gv) {
                double wg = sg_u.weight[gu] * sg_v.weight[gv];

                const Eigen::VectorXd& Nu  = sg_u.N[gu];
                const Eigen::VectorXd& dNu = sg_u.dN[gu];
                const Eigen::VectorXd& Nv  = sg_v.N[gv];
                const Eigen::VectorXd& dNv = sg_v.dN[gv];

                GaussPointGeometry g;
                if constexpr (IsRational)
                    g = evaluateGaussPointGeometryNURBS(pts, Nu, dNu, Nv, dNv, w);
                else
                    g = evaluateGaussPointGeometry(pts, Nu, dNu, Nv, dNv);

                if (g.detJ == 0.0) continue;

                double invJ11 = g.J22 / g.detJ;
                double invJ12 = -g.J12 / g.detJ;
                double invJ21 = -g.J21 / g.detJ;
                double invJ22 = g.J11 / g.detJ;

                std::vector<double> dRdx(nb_loc), dRdy(nb_loc);
                for (size_t a = 0; a < nb_loc; ++a) {
                    dRdx[a] = invJ11 * g.dRdu[a] + invJ21 * g.dRdv[a];
                    dRdy[a] = invJ12 * g.dRdu[a] + invJ22 * g.dRdv[a];
                }

                Eigen::MatrixXd term = op.computeIntegrand(g.R, dRdx, dRdy);
                if (local_contribution.size() == 0)
                    local_contribution = Eigen::MatrixXd::Zero(term.rows(), term.cols());
                local_contribution += term * wg * std::abs(g.detJ);
            }
        }
        assembleLocalContribution(local_contribution, span, tripletList);
    }
}
