#include "PatchIntegrator.hpp"
#include <algorithm>
#include <stdexcept>

// Geometry shared by the stiffness and mass kernels at one Gauss point.
// detJ == 0.0 signals a degenerate point that the caller must skip.
PatchIntegrator::GaussPointGeometry PatchIntegrator::evaluateGaussPointGeometry(
    const std::vector<const double*>& pts,
    const Eigen::VectorXd& Nu, const Eigen::VectorXd& dNu,
    const Eigen::VectorXd& Nv, const Eigen::VectorXd& dNv) const
{
    size_t nb_loc = pts.size();
    GaussPointGeometry g;
    g.R.resize(nb_loc);
    g.dRdu.resize(nb_loc);
    g.dRdv.resize(nb_loc);

    size_t idx = 0;
    for (size_t jv = 0; jv < static_cast<size_t>(Nv.size()); ++jv) {
        for (size_t iu = 0; iu < static_cast<size_t>(Nu.size()); ++iu) {
            g.R[idx] = Nu(iu) * Nv(jv);
            g.dRdu[idx] = dNu(iu) * Nv(jv);
            g.dRdv[idx] = Nu(iu) * dNv(jv);
            ++idx;
        }
    }

    // Compute mapping
    g.J11 = 0.0; g.J12 = 0.0; g.J21 = 0.0; g.J22 = 0.0;
    for (size_t a = 0; a < nb_loc; ++a) {
        const double* P = pts[a];
        const double px = P[0];
        const double py = P[1];

        g.J11 += g.dRdu[a] * px;
        g.J21 += g.dRdu[a] * py;
        g.J12 += g.dRdv[a] * px;
        g.J22 += g.dRdv[a] * py;
    }

    double detJ = g.J11*g.J22 - g.J12*g.J21;
    g.detJ = (std::abs(detJ) < 1.e-14) ? 0.0 : detJ;

    return g;
}

// Compute local stiffness contribution for a given span
Eigen::MatrixXd PatchIntegrator::computeLocalContribution(const Patch& patch, const SpanGauss1D& sg_u, const SpanGauss1D& sg_v, const std::vector<int>& span) {
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

            GaussPointGeometry g = evaluateGaussPointGeometry(pts, Nu, dNu, Nv, dNv);
            if (g.detJ == 0.0) {
                continue;
            }
            const std::vector<double>& dRdu = g.dRdu;
            const std::vector<double>& dRdv = g.dRdv;
            double detJ = g.detJ;

            double invJ11 = g.J22 / detJ;        // du/dx
            double invJ12 = - g.J12 / detJ;       // du/dy
            double invJ21 = - g.J21 / detJ;       // dv/dx
            double invJ22 = g.J11 / detJ;         // dv/dy

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
                B(0, 2*a) = grads[a][0];            // dN/dx for u_x
                B(1, 2*a + 1) = grads[a][1];        // dN/dy for u_y
                B(2, 2*a) = grads[a][1];            // dN/dy for shear (u_x)
                B(2, 2*a + 1) = grads[a][0];        // dN/dx for shear (u_y)
            }

            // Constitutive matrix D
            Eigen::Matrix3d D;
            // TODO : handle material properties with proper dedicated object
            double E = material_properties_.E;
            double nu = material_properties_.nu;
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

// Compute local mass contribution for a given span
Eigen::MatrixXd PatchIntegrator::computeLocalMassContribution(const Patch& patch, const SpanGauss1D& sg_u, const SpanGauss1D& sg_v, const std::vector<int>& span) {
    int ngauss_u = static_cast<int>(sg_u.u_param.size());
    int ngauss_v = static_cast<int>(sg_v.u_param.size());

    std::vector<const double*> pts = patch.control_points_for_span(span);
    size_t nb_loc = pts.size();

    Eigen::MatrixXd M_loc = Eigen::MatrixXd::Zero(2*nb_loc, 2*nb_loc);
    double rho = material_properties_.rho;

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

            // N matrix: same block placement as B, but with shape values
            // (no gradients/Jacobian inverse needed beyond detJ).
            Eigen::MatrixXd N = Eigen::MatrixXd::Zero(2, 2 * nb_loc);
            for (size_t a = 0; a < nb_loc; ++a) {
                N(0, 2*a) = g.R[a];
                N(1, 2*a + 1) = g.R[a];
            }

            M_loc += rho * N.transpose() * N * w * std::abs(g.detJ);
        }
    }
    return M_loc;
}

// Assemble local contriubution into global matrix
void PatchIntegrator::assembleLocalContribution(const Eigen::MatrixXd& local_contribution, const std::vector<int>& span, std::vector<Eigen::Triplet<double>>& tripletList) {
    // Get indices of control points for given span
    std::vector<const double*> pts = patch_.control_points_for_span(span);
    size_t nb_loc = pts.size();

    // Patch-LOCAL flat position (u-fastest) of each of this span's active
    // control points -- this is what PatchDOFManager::get_global_dof()
    // expects (it maps patch-local position -> global dof; it already
    // accounts for control points shared with other patches, via the global
    // ids the PatchDOFManager was built from).
    std::vector<size_t> span_local_positions = buildSpanLocalIndices(patch_, span);

    // Assemble local contribution into global matrix
    for (size_t i = 0; i < local_contribution.rows(); ++i) {
        for (size_t j = 0; j < local_contribution.cols(); ++j) {
            size_t local_control_point_i = i / patch_.dof_manager->dofs_per_control_point;
            size_t local_dof_i = i % patch_.dof_manager->dofs_per_control_point;
            size_t local_control_point_j = j / patch_.dof_manager->dofs_per_control_point;
            size_t local_dof_j = j % patch_.dof_manager->dofs_per_control_point;

            size_t global_i = patch_.dof_manager->get_global_dof(span_local_positions[local_control_point_i], local_dof_i);
            size_t global_j = patch_.dof_manager->get_global_dof(span_local_positions[local_control_point_j], local_dof_j);

            tripletList.emplace_back(global_i, global_j, local_contribution(i, j));
        }
    }
}

std::vector<size_t> PatchIntegrator::buildSpanLocalIndices(const Patch& patch, const std::vector<int>& span) const {
    std::vector<size_t> span_local_positions;

    // Get degrees
    int p_u = patch.tensor.components[0].getDegree();
    int p_v = patch.tensor.components[1].getDegree();

    // Get start index for current span
    int start_u = span[0] - p_u;
    int start_v = span[1] - p_v;

    // Get local dimensions of patch
    ssize_t n_u = patch.local_shape[0];
    ssize_t n_v = patch.local_shape[1];

    // u-fastest: direction 0 (u) fastest
    for (int jv = 0; jv <= p_v; ++jv) {
        int lv = start_v + jv;
        for (int iu = 0; iu <= p_u; ++iu) {
            int lu = start_u + iu;
            size_t local_linear = static_cast<size_t>(lv * n_u + lu);
            span_local_positions.push_back(local_linear);
        }
    }

    return span_local_positions;
}

Eigen::SparseMatrix<double> PatchIntegrator::assembleGeneric(
    const PatchAssembly& assembly,
    const std::vector<MaterialProperties>& materials,
    int gauss_n,
    const std::function<void(PatchIntegrator&, std::vector<Eigen::Triplet<double>>&)>& collect)
{
    const auto& patches = assembly.getPatchs();
    if (materials.size() != patches.size())
        throw std::invalid_argument(
            "PatchIntegrator::assembleGeneric: materials.size() must equal the "
            "number of patches in the assembly (one entry per patch, in "
            "add_patch() order).");

    std::vector<Eigen::Triplet<double>> tripletList;
    size_t total_dofs = 0;

    // Keep one IGABasis1D per patch alive until collect() has run for every
    // patch -- PatchIntegrator only stores references to its bases.
    std::vector<IGABasis1D> bases_u, bases_v;
    bases_u.reserve(patches.size());
    bases_v.reserve(patches.size());

    for (size_t p = 0; p < patches.size(); ++p) {
        const Patch& patch = *patches[p];

        int p_u = patch.tensor.components[0].getDegree();
        int p_v = patch.tensor.components[1].getDegree();
        int n_u = (gauss_n > 0) ? gauss_n : p_u + 1;
        int n_v = (gauss_n > 0) ? gauss_n : p_v + 1;

        bases_u.push_back(IGABasis1D::build(patch.tensor.components[0], n_u));
        bases_v.push_back(IGABasis1D::build(patch.tensor.components[1], n_v));

        PatchIntegrator integrator(patch, bases_u.back(), bases_v.back(), materials[p]);
        collect(integrator, tripletList);
        total_dofs = std::max(total_dofs, integrator.localTotalDofs());
    }

    Eigen::SparseMatrix<double> result(total_dofs, total_dofs);
    result.setFromTriplets(tripletList.begin(), tripletList.end());

    return result;
}

Eigen::SparseMatrix<double> PatchIntegrator::assembleStiffness(
    const PatchAssembly& assembly,
    const std::vector<MaterialProperties>& materials,
    int gauss_n)
{
    return assembleGeneric(assembly, materials, gauss_n,
        [](PatchIntegrator& pi, std::vector<Eigen::Triplet<double>>& t) { pi.collectTriplets(t); });
}

Eigen::SparseMatrix<double> PatchIntegrator::assembleMass(
    const PatchAssembly& assembly,
    const std::vector<MaterialProperties>& materials,
    int gauss_n)
{
    return assembleGeneric(assembly, materials, gauss_n,
        [](PatchIntegrator& pi, std::vector<Eigen::Triplet<double>>& t) { pi.collectMassTriplets(t); });
}