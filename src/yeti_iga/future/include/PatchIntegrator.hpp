#pragma once
#include <algorithm>
#include <vector>
#include <cstddef>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "IGABasis1D.hpp"
#include "Patch.hpp"
#include "PatchAssembly.hpp"



struct MaterialProperties {
    // TODO Generalize it
    double E;
    double nu;
    double thickness;   // Thickness (for plane problems)
};

class PatchIntegrator {
private:
    const Patch& patch_;
    // TODO: generalize it for N dimensions
    const IGABasis1D& basis_u_;
    const IGABasis1D& basis_v_;
    MaterialProperties material_properties_;

    // Compute local contrinution for a given span
    Eigen::MatrixXd computeLocalContribution(const Patch& patch, const SpanGauss1D& sg_u, const SpanGauss1D& sg_v, const std::vector<int>& span);
    // Assemble local contriubution into global matrix
    void assembleLocalContribution(const Eigen::MatrixXd& local_contribution, const std::vector<int>& span, std::vector<Eigen::Triplet<double>>& tripletList);

    // Patch-LOCAL flat positions (u-fastest) of the span's active control points.
    std::vector<size_t> buildSpanLocalIndices(const Patch& patch, const std::vector<int>& span) const;

public:
    PatchIntegrator(const Patch& patch, const IGABasis1D& basis_u, const IGABasis1D& basis_v, const MaterialProperties& material_properties)
        : patch_(patch), basis_u_(basis_u), basis_v_(basis_v), material_properties_(material_properties) {}

    // Append this patch's contribution to a (possibly shared) triplet list,
    // using GLOBAL dof indices. Does not size or build any matrix -- this is
    // what lets assemble() merge several patches' contributions before a
    // single setFromTriplets() call (see assemble()).
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
            Eigen::MatrixXd local_contribution = computeLocalContribution(patch_, sg_u, sg_v, span);

            // Assemble local contribution into global matrix
            assembleLocalContribution(local_contribution, span, tripletList);
        }
    }

    // This patch's own highest global dof + 1 (scanning every local control
    // point, not just the last one -- global dof order need not follow
    // local position order once global_indices isn't the identity mapping,
    // e.g. for any patch that isn't the first one added to a PatchAssembly).
    // For a single, exclusively-owned patch this is the full matrix size
    // (see integrate()); for a patch that is part of a PatchAssembly, the
    // assembly-wide total is generally larger -- see assemble().
    size_t localTotalDofs() const {
        const auto& local_to_global = patch_.dof_manager->local_to_global_dofs;
        return *std::max_element(local_to_global.begin(), local_to_global.end()) + 1;
    }

    Eigen::SparseMatrix<double> integrate() {
        std::vector<Eigen::Triplet<double>> tripletList;
        collectTriplets(tripletList);

        size_t total_dofs = localTotalDofs();

        // build sparse matrix from triplets
        Eigen::SparseMatrix<double> stiffness_matrix(total_dofs, total_dofs);
        stiffness_matrix.setFromTriplets(tripletList.begin(), tripletList.end());

        return stiffness_matrix;
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
    static Eigen::SparseMatrix<double> assemble(
        const PatchAssembly& assembly,
        const std::vector<MaterialProperties>& materials,
        int gauss_n = 0);

};

