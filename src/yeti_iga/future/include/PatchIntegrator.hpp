#pragma once
#include <vector>
#include <cstddef>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "IGABasis1D.hpp"
#include "Patch.hpp"



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

    // Compute local contrinution for a given span
    Eigen::MatrixXd computeLocalContribution(const Patch& patch, const SpanGauss1D& sg_u, const SpanGauss1D& sg_v, const std::vector<int>& span);
    // Assemble local contriubution into global matrix
    void assembleLocalContribution(const Eigen::MatrixXd& local_contribution, const std::vector<int>& span, std::vector<Eigen::Triplet<double>>& tripletList);

    std::vector<size_t> buildLocalToGlobalMapping(const Patch& patch, const std::vector<int>& span) const;

public:
    PatchIntegrator(const Patch& patch, const IGABasis1D& basis_u, const IGABasis1D& basis_v)
        : patch_(patch), basis_u_(basis_u), basis_v_(basis_v) {}

    Eigen::SparseMatrix<double> integrate() {
        std::vector<Eigen::Triplet<double>> tripletList;
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

        // Set size of global matrix
        size_t total_dofs = patch_.dof_manager->get_global_dof_indices(patch_.global_indices.back()).back() + 1;

        // build sparse matrix from triplets
        Eigen::SparseMatrix<double> stiffness_matrix(total_dofs, total_dofs);
        stiffness_matrix.setFromTriplets(tripletList.begin(), tripletList.end());

        return stiffness_matrix;
    }

};

