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
    Eigen::MatrixXd computeLocalContribution(const Patch& patch, const SpanGauss1D& sg_u, const SpanGauss1D& sg_v, const std::vector<int>& span) {
        int ngauss_u = static_cast<int>(sg_u.u_param.size());
        int ngauss_v = static_cast<int>(sg_v.u_param.size());

        // get pointers to active CP of current span
        std::vector<const double*> pts = patch.control_points_for_span(span);
        size_t nb_loc = pts.size();

        // Initialize elementary stiffness matrix
        Eigen::MatrixXd K_loc = Eigen::MatrixXd::Zero(2*nb_loc, 2*nb_loc);

        // Loop over Gauss points
        for (int gu = 0; gu < ngauss_u; ++gu) {
            for (int gv = 0; gv < ngauss_v; ++gv) {
                double w = sg_u.weight[gu] * sg_v.weight[gv];

                // get basis functions and derivatives
                const Eigen::VectorXd& Nu = sg_u.N[gu];
                const Eigen::VectorXd& dNu = sg_u.dN[gu];
                const Eigen::VectorXd& Nv = sg_v.N[gv];
                const Eigen::VectorXd& dNv = sg_v.dN[gv];

                // Compute R, dR/du, dR/dv
                std::vector<double> R(nb_loc);
                std::vector<double> dRdu(nb_loc);
                std::vector<double> dRdv(nb_loc);

                size_t idx = 0;
                for (size_t jv = 0; jv < Nv.size(); ++jv) {
                    for (size_t iu = 0; iu < Nu.size(); ++iu) {
                        R[idx] = Nu(iu) * Nv(jv);
                        dRdu[idx] = dNu(iu) * Nv(jv);
                        dRdv[idx] = Nu(iu) * dNv(jv);
                        ++idx;
                    }
                }

                // Compute mapping
                double J11 = 0.0, J12 = 0.0, J21 = 0.0, J22 = 0.0;
                for (size_t a = 0; a < nb_loc; ++a) {
                    const double* P = pts[a];
                    const double px = P[0];
                    const double py = P[1];

                    J11 += dRdu[a] * px;
                    J21 += dRdu[a] * py;
                    J12 += dRdv[a] * px;
                    J22 += dRdv[a] * py;

                }

                double detJ = J11*J22 - J12*J21;
                if (std::abs(detJ) < 1.e-14) {
                    continue;
                }

                double invJ11 = J22 / detJ;         // du/dx
                double invJ12 = - J12 / detJ;       // du/dy
                double invJ21 = - J21 / detJ;       // dv/dx
                double invJ22 = J11 / detJ;         // dv/dy

                // Compute gradients
                std::vector<std::array<double, 2>> grads(nb_loc);
                for (size_t a = 0; a < nb_loc; ++a) {
                    grads[a][0] = invJ11 * dRdu[a] + invJ21 * dRdv[a];  // dRdx
                    grads[a][1] = invJ12 * dRdu[a] + invJ22 * dRdv[a];  // dRdx
                }

                // Build matrix B
                // TODO implement B free ...
                Eigen::MatrixXd B = Eigen::MatrixXd::Zero(3, 2 * nb_loc);
                for (size_t a = 0; a < nb_loc; ++a) {
                    // TODO : those affectation must be verified (Voigt convention ?)
                    B(0, 2*a) = grads[a][0];            // dN/dx for u_x
                    B(1, 2*a + 1) = grads[a][1];        // dN/dy for u_y
                    B(2, 2*a) = grads[a][1];            // dN/dy for shear (u_x)
                    B(2, 2*a + 1) = grads[a][0];        // dN/dx for shear (u_y)
                }

                // Constitutive matrix D
                // TODO : verify if Voigt convention is enforced
                Eigen::Matrix3d D;
                // TODO : handle material properties with proper dedicated object
                double E = 210000.0;
                double nu = 0.3;
                double factor = E / (1.0 - nu*nu);
                D <<
                factor, factor * nu, 0.0,
                factor * nu, factor, 0.0,
                0.0, 0.0, factor * (1.0 - nu) / 2.0;

                // Local contribution : B^T * D * B * w * detJ
                K_loc += B.transpose() * D * B * w * std::abs(detJ);

            }
        }
        return K_loc;
    }

    // Assemble local contriubution into global matrix
    void assembleLocalContribution(const Eigen::MatrixXd& local_contribution, const std::vector<int>& span, std::vector<Eigen::Triplet<double>>& tripletList) {
        // Get indices of control points for given span
        std::vector<const double*> pts = patch_.control_points_for_span(span);
        size_t nb_loc = pts.size();

        // Get indices of globel DOFs
        std::vector<size_t> local_to_global = buildLocalToGlobalMapping(patch_, span);

        // Assemble local contribution into global matrix
        for (size_t i = 0; i < local_contribution.rows(); ++i) {
            for (size_t j = 0; j < local_contribution.cols(); ++j) {
                size_t local_control_point_i = i / patch_.dof_manager->dofs_per_control_point;
                size_t local_dof_i = i % patch_.dof_manager->dofs_per_control_point;
                size_t local_control_point_j = j / patch_.dof_manager->dofs_per_control_point;
                size_t local_dof_j = j % patch_.dof_manager->dofs_per_control_point;

                size_t global_i = patch_.dof_manager->get_global_dof_indices(local_to_global[local_control_point_i])[local_dof_i];
                size_t global_j = patch_.dof_manager->get_global_dof_indices(local_to_global[local_control_point_j])[local_dof_j];

                tripletList.emplace_back(global_i, global_j, local_contribution(i, j));
            }
        }
    }

    std::vector<size_t> buildLocalToGlobalMapping(const Patch& patch, const std::vector<int>& span) const {
        // TODO This function could be improved AND VERIFIED

        std::vector<size_t> local_to_global;

        // Get degrees
        int p_u = patch.tensor.components[0].getDegree();
        int p_v = patch.tensor.components[1].getDegree();

        // Get start index for current span
        int start_u = span[0] - p_u;
        int start_v = span[1] - p_v;

        // Get local dimensions of patch
        ssize_t n_u = patch.local_shape[0];
        ssize_t n_v = patch.local_shape[1];

        // Compute number of control points for curent patch
        size_t nb_loc = (p_u + 1) * (p_v + 1);

        // Fill local to global mapping patch.global_indices
        for (int jv = 0; jv <= p_v; ++jv) {
            int lv = start_v + jv;
            for (int iu = 0; iu <= p_u; ++iu) {
                int lu = start_u + iu;
                // Compute local linear index
                size_t local_linear = static_cast<size_t>(lv * n_u + lu);
                // Get global index
                size_t global_index = patch.global_indices[local_linear];
                local_to_global.push_back(global_index);
            }
        }

        return local_to_global;
    }

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

