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



struct MaterialProperties {
    // TODO Generalize it
    double E;
    double nu;
    double thickness = 1.0;   // Thickness (for plane problems)
    double rho = 0.0;         // Mass density (required by integrateMass()/assembleMass())
};

class PatchIntegrator {
private:
    const Patch& patch_;
    // TODO: generalize it for N dimensions
    const IGABasis1D& basis_u_;
    const IGABasis1D& basis_v_;
    MaterialProperties material_properties_;

    // Geometry shared by the stiffness and mass kernels at one Gauss point:
    // basis values/derivatives in physical-mapping order (u-fastest), the
    // Jacobian components and its determinant. Mass only needs R and detJ;
    // stiffness derives invJ/grads/B/D from J11..J22 and dRdu/dRdv on top of
    // this (storing J11..J22 here, instead of just detJ, avoids recomputing
    // the Jacobian a second time in computeLocalStiffnessContribution). detJ == 0.0
    // signals a degenerate point to be skipped by the caller (mirrors the
    // previous inline `if (std::abs(detJ) < 1e-14) continue;`).
    struct GaussPointGeometry {
        std::vector<double> R, dRdu, dRdv;
        double J11, J12, J21, J22;
        double detJ;
    };
    GaussPointGeometry evaluateGaussPointGeometry(
        const std::vector<const double*>& pts,
        const Eigen::VectorXd& Nu, const Eigen::VectorXd& dNu,
        const Eigen::VectorXd& Nv, const Eigen::VectorXd& dNv) const;

    // Compute local stiffness contribution for a given span
    Eigen::MatrixXd computeLocalStiffnessContribution(const Patch& patch, const SpanGauss1D& sg_u, const SpanGauss1D& sg_v, const std::vector<int>& span);
    // Compute local mass contribution for a given span
    Eigen::MatrixXd computeLocalMassContribution(const Patch& patch, const SpanGauss1D& sg_u, const SpanGauss1D& sg_v, const std::vector<int>& span);
    // Assemble local contriubution into global matrix (shared by stiffness and mass)
    void assembleLocalContribution(const Eigen::MatrixXd& local_contribution, const std::vector<int>& span, std::vector<Eigen::Triplet<double>>& tripletList);

    // Patch-LOCAL flat positions (u-fastest) of the span's active control points.
    std::vector<size_t> buildSpanLocalIndices(const Patch& patch, const std::vector<int>& span) const;

    // Shared orchestration for assemble()/assembleMass(): validates
    // materials, builds one IGABasis1D pair + PatchIntegrator per patch,
    // invokes `collect` (collectTriplets or collectMassTriplets) into one
    // shared triplet list, sizes the result to the assembly-wide total dof
    // count. Called once per patch, not per Gauss point -- the std::function
    // indirection has no measurable cost here.
    static Eigen::SparseMatrix<double> assembleGeneric(
        const PatchAssembly& assembly,
        const std::vector<MaterialProperties>& materials,
        int gauss_n,
        const std::function<void(PatchIntegrator&, std::vector<Eigen::Triplet<double>>&)>& collect);

public:
    PatchIntegrator(const Patch& patch, const IGABasis1D& basis_u, const IGABasis1D& basis_v, const MaterialProperties& material_properties)
        : patch_(patch), basis_u_(basis_u), basis_v_(basis_v), material_properties_(material_properties) {}

    // Append this patch's contribution to a (possibly shared) triplet list,
    // using GLOBAL dof indices. Does not size or build any matrix -- this is
    // what lets assembleStiffness() merge several patches' contributions
    // before a single setFromTriplets() call (see assembleStiffness()).
    void collectTriplets(std::vector<Eigen::Triplet<double>>& tripletList) {
        SpanNDIterator it = patch_.spans();

        for (auto span : it) {
            int span_u = span[0];
            int span_v = span[1];

            // Find spans indices in pre-computed basis
            int idx_u = basis_u_.span_indices.at(span_u);
            int idx_v = basis_v_.span_indices.at(span_v);

            const SpanGauss1D& sg_u = basis_u_.gauss_spans[idx_u];
            const SpanGauss1D& sg_v = basis_v_.gauss_spans[idx_v];

            // Compute local contribution
            Eigen::MatrixXd local_contribution = computeLocalStiffnessContribution(patch_, sg_u, sg_v, span);

            // Assemble local contribution into global matrix
            assembleLocalContribution(local_contribution, span, tripletList);
        }
    }

    // Same as collectTriplets(), but for the consistent mass matrix
    // contribution (requires material_properties_.rho > 0).
    void collectMassTriplets(std::vector<Eigen::Triplet<double>>& tripletList) {
        if (material_properties_.rho <= 0.0)
            throw std::invalid_argument(
                "PatchIntegrator::collectMassTriplets: material_properties_.rho "
                "must be set (> 0) to assemble a mass matrix.");

        SpanNDIterator it = patch_.spans();

        for (auto span : it) {
            int span_u = span[0];
            int span_v = span[1];

            int idx_u = basis_u_.span_indices.at(span_u);
            int idx_v = basis_v_.span_indices.at(span_v);

            const SpanGauss1D& sg_u = basis_u_.gauss_spans[idx_u];
            const SpanGauss1D& sg_v = basis_v_.gauss_spans[idx_v];

            Eigen::MatrixXd local_contribution = computeLocalMassContribution(patch_, sg_u, sg_v, span);

            assembleLocalContribution(local_contribution, span, tripletList);
        }
    }

    // This patch's own highest global dof + 1 (scanning every local control
    // point, not just the last one -- global dof order need not follow
    // local position order once global_indices isn't the identity mapping,
    // e.g. for any patch that isn't the first one added to a PatchAssembly).
    // For a single, exclusively-owned patch this is the full matrix size
    // (see integrateStiffness()); for a patch that is part of a
    // PatchAssembly, the assembly-wide total is generally larger -- see
    // assembleStiffness().
    size_t localTotalDofs() const {
        const auto& local_to_global = patch_.dof_manager->local_to_global_dofs;
        return *std::max_element(local_to_global.begin(), local_to_global.end()) + 1;
    }

    Eigen::SparseMatrix<double> integrateStiffness() {
        std::vector<Eigen::Triplet<double>> tripletList;
        collectTriplets(tripletList);

        size_t total_dofs = localTotalDofs();

        // build sparse matrix from triplets
        Eigen::SparseMatrix<double> stiffness_matrix(total_dofs, total_dofs);
        stiffness_matrix.setFromTriplets(tripletList.begin(), tripletList.end());

        return stiffness_matrix;
    }

    // Same as integrateStiffness(), but for the consistent mass matrix
    // (requires material_properties_.rho > 0).
    Eigen::SparseMatrix<double> integrateMass() {
        std::vector<Eigen::Triplet<double>> tripletList;
        collectMassTriplets(tripletList);

        size_t total_dofs = localTotalDofs();

        Eigen::SparseMatrix<double> mass_matrix(total_dofs, total_dofs);
        mass_matrix.setFromTriplets(tripletList.begin(), tripletList.end());

        return mass_matrix;
    }

    // Assemble the global stiffness matrix of a whole PatchAssembly: builds
    // one PatchIntegrator per patch (with its own IGABasis1D in each
    // direction and its own MaterialProperties), merges every patch's
    // triplets into ONE list, and sizes the result to the assembly-wide
    // total dof count. Control points shared between patches already
    // resolve to the same global dof (see GlobalDOFManager/PatchAssembly
    // documentation) and Eigen::setFromTriplets() sums duplicate (row, col)
    // entries automatically, so shared-boundary coupling falls out of this
    // merge with no special-casing.
    //
    // materials must have one entry per patch, in assembly.get_patchs()
    // (add_patch()) order.
    //
    // gauss_n: number of Gauss points per span, per direction, used for
    // every patch. If 0 (default), each direction of each patch uses
    // degree + 1 (sufficient to integrate the stiffness bilinear form
    // exactly, see IGABasis1D::build()).
    static Eigen::SparseMatrix<double> assembleStiffness(
        const PatchAssembly& assembly,
        const std::vector<MaterialProperties>& materials,
        int gauss_n = 0);

    // Same as assembleStiffness(), but for the consistent mass matrix of
    // the whole PatchAssembly. Every entry in `materials` must have rho > 0.
    static Eigen::SparseMatrix<double> assembleMass(
        const PatchAssembly& assembly,
        const std::vector<MaterialProperties>& materials,
        int gauss_n = 0);

    // Same as collectTriplets()/collectMassTriplets(), but the per-Gauss-
    // point term comes from a caller-supplied LocalOperator instead of a
    // built-in kernel: PatchIntegrator still does the Gauss-point loop and
    // the Jacobian/gradient computation (via the shared
    // evaluateGaussPointGeometry() helper -- same one used by
    // computeLocalStiffnessContribution()/computeLocalMassContribution()), the
    // operator only supplies the term to weight and sum. Development/
    // testing convenience -- does not touch, and is not used by,
    // collectTriplets()/collectMassTriplets().
    void collectOperatorTriplets(LocalOperator& op, std::vector<Eigen::Triplet<double>>& tripletList) {
        SpanNDIterator it = patch_.spans();

        for (auto span : it) {
            int span_u = span[0];
            int span_v = span[1];

            int idx_u = basis_u_.span_indices.at(span_u);
            int idx_v = basis_v_.span_indices.at(span_v);

            const SpanGauss1D& sg_u = basis_u_.gauss_spans[idx_u];
            const SpanGauss1D& sg_v = basis_v_.gauss_spans[idx_v];

            std::vector<const double*> pts = patch_.control_points_for_span(span);
            size_t nb_loc = pts.size();

            int ngauss_u = static_cast<int>(sg_u.u_param.size());
            int ngauss_v = static_cast<int>(sg_v.u_param.size());

            Eigen::MatrixXd local_contribution;  // sized from the operator's first returned term

            for (int gu = 0; gu < ngauss_u; ++gu) {
                for (int gv = 0; gv < ngauss_v; ++gv) {
                    double w = sg_u.weight[gu] * sg_v.weight[gv];

                    const Eigen::VectorXd& Nu = sg_u.N[gu];
                    const Eigen::VectorXd& dNu = sg_u.dN[gu];
                    const Eigen::VectorXd& Nv = sg_v.N[gv];
                    const Eigen::VectorXd& dNv = sg_v.dN[gv];

                    GaussPointGeometry g = evaluateGaussPointGeometry(pts, Nu, dNu, Nv, dNv);
                    if (g.detJ == 0.0) {
                        continue;
                    }

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
                    if (local_contribution.size() == 0) {
                        local_contribution = Eigen::MatrixXd::Zero(term.rows(), term.cols());
                    }
                    local_contribution += term * w * std::abs(g.detJ);
                }
            }

            assembleLocalContribution(local_contribution, span, tripletList);
        }
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
    // order -- the same convention as `materials` in assembleStiffness()/
    // assembleMass(). Deliberately a standalone implementation (does not
    // call assembleGeneric()), so that this development/testing path can
    // never affect assembleStiffness()/assembleMass() behavior.
    static Eigen::SparseMatrix<double> assembleOperator(
        const PatchAssembly& assembly,
        const std::vector<std::shared_ptr<LocalOperator>>& operators,
        int gauss_n = 0);

};

