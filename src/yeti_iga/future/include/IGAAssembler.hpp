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
                const Eigen::VectorXd& dNu = sg_u.N[gu];
                const Eigen::VectorXd& Nv = sg_v.N[gv];
                const Eigen::VectorXd& dNv = sg_v.dN[gv];

                // Compute R, dR/du, dR/dv
                std::vector<double> R(nb_loc);
                std::vector<double> dRdu(nb_loc);
                std::vector<double> dRdv(nb_loc);

                for (size_t i=0; i < nb_loc; ++i) {
                    R[i] = Nu(i) * Nv(i);
                    dRdu[i] = dNu(i) * Nv(i);
                    dRdv[i] = Nu(i) * dNv(i);
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

                double invJ11 = J22 / detJ;
                double invJ12 = - J12 / detJ;
                double invJ21 = - J21 / detJ;
                double invJ22 = J11 / detJ;

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
        for (size_t i = 0; i < 2 * nb_loc; ++i) {
            for (size_t j = 0; j < 2 * nb_loc; ++j) {
                // TODO check the utility of [i / 2] and [i % 2]. Is the number of DOF per CP is taken into account according to data structure ?
                size_t global_i = patch_.dof_manager->get_global_dof_indices(local_to_global[i / 2])[i % 2];
                size_t global_j = patch_.dof_manager->get_global_dof_indices(local_to_global[j / 2])[j % 2];
                tripletList.emplace_back(global_i, global_j, local_contribution(i, j));
            }
        }
    }

    // Build local to global mapping
    std::vector<size_t> buildLocalToGlobalMapping(const Patch& patch, const std::vector<int>& span) const {
        std::vector<size_t> local_to_global;
        for (size_t i = 0; i < span.size(); ++i) {
            local_to_global.push_back(span[i]);
        }
        return local_to_global;
    }

public:
    PatchIntegrator(const Patch& patch, const IGABasis1D& basis_u, const IGABasis1D& basis_v)
        : patch_(patch), basis_u_(basis_u), basis_v_(basis_v) {}

};





// struct ElementMatrix {
//     // global indices for the local basis
//     std::vector<size_t> global_indices;
//     // dense local stiffness
//     Eigen::MatrixXd K;
//     size_t nb_loc = 0;

//     // Constructor (initialize K with required shape)
//     ElementMatrix(size_t size) : nb_loc(size), K(2*size, 2*size) {
//         K.setZero();
//     }
// };

// class IGAAssembler2D {
// public:
//     // Construct with the patch and pre-computed per-direction IGABasis1D
//     IGAAssembler2D(const Patch& patch,
//                    const IGABasis1D& basis_u,
//                    const IGABasis1D& basis_v)
//         : patch_(patch), basis_u_(basis_u), basis_v_(basis_v)
//     {
//         // sanity checks
//         assert(patch_.local_shape.size() == 2);
//         p_u_ = patch_.tensor.components[0].getDegree();
//         p_v_ = patch_.tensor.components[1].getDegree();
//         nb_loc_ = (p_u_+1) * (p_v_+1);

//         // Hard coding of material properties, TODO : should be set properly
//         material_.E = 210000.;
//         material_.nu = 0.3;
//     }

//     // Assemble per-element stiffness blocks and return them
//     // (it DOES NOT assemble into a global sparse matrix)
//     std::vector<ElementMatrix> assemble_stiffness() const;

//     Eigen::Matrix3d computeConstitutiveMatrix() const;

// private:
//     const Patch& patch_;
//     const IGABasis1D& basis_u_;
//     const IGABasis1D& basis_v_;
//     int p_u_, p_v_;
//     size_t nb_loc_;
//     MaterialProperties material_;       // Should be defined as a member of patch ?
// };
