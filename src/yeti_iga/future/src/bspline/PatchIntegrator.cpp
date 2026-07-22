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
Eigen::MatrixXd PatchIntegrator::computeLocalStiffnessContribution(const Patch& patch, const SpanGauss1D& sg_u, const SpanGauss1D& sg_v, const std::vector<int>& span) {
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

// Standalone orchestration for assembleOperator() -- intentionally does not
// call assembleGeneric() (which is tied to MaterialProperties), so this
// development/testing path cannot affect assembleStiffness()/assembleMass().
Eigen::SparseMatrix<double> PatchIntegrator::assembleOperator(
    const PatchAssembly& assembly,
    const std::vector<std::shared_ptr<LocalOperator>>& operators,
    int gauss_n)
{
    const auto& patches = assembly.getPatchs();
    if (operators.size() != patches.size())
        throw std::invalid_argument(
            "PatchIntegrator::assembleOperator: operators.size() must equal the "
            "number of patches in the assembly (one entry per patch, in "
            "add_patch() order).");

    std::vector<Eigen::Triplet<double>> tripletList;
    size_t total_dofs = 0;

    std::vector<IGABasis1D> bases_u, bases_v;
    bases_u.reserve(patches.size());
    bases_v.reserve(patches.size());

    MaterialProperties unused_material{0.0, 0.0};

    for (size_t p = 0; p < patches.size(); ++p) {
        const Patch& patch = *patches[p];

        int p_u = patch.tensor.components[0].getDegree();
        int p_v = patch.tensor.components[1].getDegree();
        int n_u = (gauss_n > 0) ? gauss_n : p_u + 1;
        int n_v = (gauss_n > 0) ? gauss_n : p_v + 1;

        bases_u.push_back(IGABasis1D::build(patch.tensor.components[0], n_u));
        bases_v.push_back(IGABasis1D::build(patch.tensor.components[1], n_v));

        PatchIntegrator integrator(patch, bases_u.back(), bases_v.back(), unused_material);
        integrator.collectOperatorTriplets(*operators[p], tripletList);
        total_dofs = std::max(total_dofs, integrator.localTotalDofs());
    }

    Eigen::SparseMatrix<double> result(total_dofs, total_dofs);
    result.setFromTriplets(tripletList.begin(), tripletList.end());

    return result;
}

// Local contribution for one boundary span. R/dR-w.r.t.-the-varying-
// parameter use the same u-fastest tensor product as
// evaluateGaussPointGeometry(), except one of the two factors (the fixed
// direction's) is the boundary-evaluated row passed in, with a zero
// derivative -- the chain rule below collapses to the right term in either
// case (direction == 0 or 1) without a separate branch for each.
Eigen::VectorXd PatchIntegrator::computeLocalBoundaryLoadContribution(
    const Patch& patch, const std::vector<int>& span, int direction,
    const SpanGauss1D& sg_varying, const Eigen::VectorXd& N_fixed_boundary,
    const Traction& traction)
{
    std::vector<const double*> pts = patch.control_points_for_span(span);
    size_t nb_loc = pts.size();
    int ngauss = static_cast<int>(sg_varying.u_param.size());

    Eigen::VectorXd F_loc = Eigen::VectorXd::Zero(2 * nb_loc);
    Eigen::VectorXd zero_deriv = Eigen::VectorXd::Zero(N_fixed_boundary.size());

    for (int g = 0; g < ngauss; ++g) {
        double w = sg_varying.weight[g];

        const Eigen::VectorXd& Nu  = (direction == 0) ? N_fixed_boundary : sg_varying.N[g];
        const Eigen::VectorXd& dNu = (direction == 0) ? zero_deriv       : sg_varying.dN[g];
        const Eigen::VectorXd& Nv  = (direction == 1) ? N_fixed_boundary : sg_varying.N[g];
        const Eigen::VectorXd& dNv = (direction == 1) ? zero_deriv       : sg_varying.dN[g];

        std::vector<double> R(nb_loc), dRdvar(nb_loc);
        size_t idx = 0;
        for (size_t jv = 0; jv < static_cast<size_t>(Nv.size()); ++jv) {
            for (size_t iu = 0; iu < static_cast<size_t>(Nu.size()); ++iu) {
                R[idx] = Nu(iu) * Nv(jv);
                // Exactly one of the two terms below is nonzero, depending
                // on which direction is fixed (its dN is the zero vector).
                dRdvar[idx] = dNu(iu) * Nv(jv) + Nu(iu) * dNv(jv);
                ++idx;
            }
        }

        Eigen::Vector2d T = Eigen::Vector2d::Zero();
        Eigen::Vector2d X = Eigen::Vector2d::Zero();
        for (size_t a = 0; a < nb_loc; ++a) {
            const double* P = pts[a];
            T[0] += dRdvar[a] * P[0];
            T[1] += dRdvar[a] * P[1];
            X[0] += R[a] * P[0];
            X[1] += R[a] * P[1];
        }
        double ds = T.norm();
        if (ds < 1.e-14) continue;

        Eigen::Vector2d t = traction.evaluate(X);

        for (size_t a = 0; a < nb_loc; ++a) {
            F_loc[2*a]     += R[a] * t[0] * w * ds;
            F_loc[2*a + 1] += R[a] * t[1] * w * ds;
        }
    }

    return F_loc;
}

void PatchIntegrator::assembleLocalLoadContribution(
    const Eigen::VectorXd& local_load, const std::vector<int>& span, Eigen::VectorXd& global_load)
{
    std::vector<size_t> span_local_positions = buildSpanLocalIndices(patch_, span);

    for (size_t i = 0; i < static_cast<size_t>(local_load.size()); ++i) {
        size_t local_control_point = i / patch_.dof_manager->dofs_per_control_point;
        size_t local_dof = i % patch_.dof_manager->dofs_per_control_point;
        size_t global_i = patch_.dof_manager->get_global_dof(span_local_positions[local_control_point], local_dof);
        global_load[global_i] += local_load[i];
    }
}

Eigen::VectorXd PatchIntegrator::integrateBoundaryLoad(
    int direction, int side, const Traction& traction, int span_min, int span_max)
{
    if (direction != 0 && direction != 1)
        throw std::invalid_argument("PatchIntegrator::integrateBoundaryLoad: direction must be 0 or 1.");
    if (side != 0 && side != 1)
        throw std::invalid_argument("PatchIntegrator::integrateBoundaryLoad: side must be 0 (min) or 1 (max).");

    int varying = 1 - direction;
    const IGABasis1D& basis_varying = (direction == 0) ? basis_v_ : basis_u_;
    const BSpline& bspline_fixed = patch_.tensor.components[direction];

    // First/last valid span of the fixed direction (same "valid span" scan
    // as SpanNDIterator/IGABasis1D::build).
    const auto& kv_fixed = bspline_fixed.getKnotVector();
    int p_fixed = bspline_fixed.getDegree();
    int m_fixed = static_cast<int>(kv_fixed.size()) - 1;
    int span_fixed = -1;
    if (side == 0) {
        for (int i = p_fixed; i <= m_fixed - p_fixed - 1; ++i)
            if (kv_fixed[i+1] > kv_fixed[i]) { span_fixed = i; break; }
    } else {
        for (int i = m_fixed - p_fixed - 1; i >= p_fixed; --i)
            if (kv_fixed[i+1] > kv_fixed[i]) { span_fixed = i; break; }
    }
    if (span_fixed < 0)
        throw std::runtime_error(
            "PatchIntegrator::integrateBoundaryLoad: no valid span found for the fixed direction.");

    double u_boundary = (side == 0) ? kv_fixed.front() : kv_fixed.back();
    std::vector<double> N_fixed_vals(p_fixed + 1);
    bspline_fixed.BasisFuns_raw(span_fixed, u_boundary, N_fixed_vals.data());
    Eigen::VectorXd N_fixed_boundary(p_fixed + 1);
    for (int i = 0; i <= p_fixed; ++i) N_fixed_boundary[i] = N_fixed_vals[i];

    Eigen::VectorXd global_load = Eigen::VectorXd::Zero(localTotalDofs());

    const bool rational = patch_.cp_manager->is_rational();
    SpanNDIterator it = patch_.spans();
    for (auto span : it) {
        if (span[direction] != span_fixed) continue;
        if (span_min >= 0 && (span[varying] < span_min || span[varying] > span_max)) continue;

        int idx_varying = basis_varying.span_indices.at(span[varying]);
        const SpanGauss1D& sg_varying = basis_varying.gauss_spans[idx_varying];

        Eigen::VectorXd local_load;
        if (rational) {
            auto w = patch_.weights_for_span(span);
            local_load = computeLocalBoundaryLoadContributionNURBS(
                patch_, span, direction, sg_varying, N_fixed_boundary, traction, w);
        } else {
            local_load = computeLocalBoundaryLoadContribution(
                patch_, span, direction, sg_varying, N_fixed_boundary, traction);
        }

        assembleLocalLoadContribution(local_load, span, global_load);
    }

    return global_load;
}

// Standalone orchestration -- does not call assembleGeneric() (one entry
// per spec, not one per patch, unlike assembleStiffness()/assembleMass()).
Eigen::VectorXd PatchIntegrator::assembleBoundaryLoad(
    const PatchAssembly& assembly,
    const std::vector<BoundaryLoadSpec>& specs,
    int gauss_n)
{
    const auto& patches = assembly.getPatchs();

    // Size to the assembly-wide total dof count, independent of which
    // patches the specs actually touch, so the result composes directly
    // with assembleStiffness()/assembleMass() (same total size).
    size_t total_dofs = 0;
    for (const auto& patch : patches) {
        if (!patch->dof_manager) continue;
        const auto& l2g = patch->dof_manager->local_to_global_dofs;
        if (!l2g.empty())
            total_dofs = std::max(total_dofs, *std::max_element(l2g.begin(), l2g.end()) + 1);
    }
    Eigen::VectorXd global_load = Eigen::VectorXd::Zero(total_dofs);

    MaterialProperties unused_material{0.0, 0.0};
    std::vector<IGABasis1D> bases_u, bases_v;
    bases_u.reserve(specs.size());
    bases_v.reserve(specs.size());

    for (const auto& spec : specs) {
        if (spec.patch_index >= patches.size())
            throw std::out_of_range("PatchIntegrator::assembleBoundaryLoad: patch_index out of range.");

        const Patch& patch = *patches[spec.patch_index];

        int p_u = patch.tensor.components[0].getDegree();
        int p_v = patch.tensor.components[1].getDegree();
        int n_u = (gauss_n > 0) ? gauss_n : p_u + 1;
        int n_v = (gauss_n > 0) ? gauss_n : p_v + 1;

        bases_u.push_back(IGABasis1D::build(patch.tensor.components[0], n_u));
        bases_v.push_back(IGABasis1D::build(patch.tensor.components[1], n_v));

        PatchIntegrator integrator(patch, bases_u.back(), bases_v.back(), unused_material);
        Eigen::VectorXd local_result = integrator.integrateBoundaryLoad(
            spec.direction, spec.side, *spec.traction, spec.span_min, spec.span_max);

        global_load.head(local_result.size()) += local_result;
    }

    return global_load;
}

// ─────────────────────────────────────────────────────────────────────────────
// NURBS implementations
// ─────────────────────────────────────────────────────────────────────────────

PatchIntegrator::GaussPointGeometry PatchIntegrator::evaluateGaussPointGeometryNURBS(
    const std::vector<const double*>& pts,
    const Eigen::VectorXd& Nu, const Eigen::VectorXd& dNu,
    const Eigen::VectorXd& Nv, const Eigen::VectorXd& dNv,
    const std::vector<double>& weights) const
{
    size_t nb_loc = pts.size();
    GaussPointGeometry g;
    g.R.resize(nb_loc);
    g.dRdu.resize(nb_loc);
    g.dRdv.resize(nb_loc);

    // Raw tensor-product B-spline basis (same as evaluateGaussPointGeometry)
    size_t idx = 0;
    for (size_t jv = 0; jv < static_cast<size_t>(Nv.size()); ++jv) {
        for (size_t iu = 0; iu < static_cast<size_t>(Nu.size()); ++iu) {
            g.R[idx] = Nu(iu) * Nv(jv);
            g.dRdu[idx] = dNu(iu) * Nv(jv);
            g.dRdv[idx] = Nu(iu) * dNv(jv);
            ++idx;
        }
    }

    // NURBS rationalization (quotient rule: R_a = w_a N_a / W)
    double W = 0.0, dWdu = 0.0, dWdv = 0.0;
    for (size_t a = 0; a < nb_loc; ++a) {
        W    += weights[a] * g.R[a];
        dWdu += weights[a] * g.dRdu[a];
        dWdv += weights[a] * g.dRdv[a];
    }
    const double inv_W = 1.0 / W;
    for (size_t a = 0; a < nb_loc; ++a) {
        const double Ra = weights[a] * g.R[a] * inv_W;
        g.dRdu[a] = (weights[a] * g.dRdu[a] - Ra * dWdu) * inv_W;
        g.dRdv[a] = (weights[a] * g.dRdv[a] - Ra * dWdv) * inv_W;
        g.R[a] = Ra;  // written last: derivative above uses old g.R[a]
    }

    // Jacobian from rationalized R/dRdu/dRdv (same formula as B-spline)
    g.J11 = 0.0; g.J12 = 0.0; g.J21 = 0.0; g.J22 = 0.0;
    for (size_t a = 0; a < nb_loc; ++a) {
        const double* P = pts[a];
        g.J11 += g.dRdu[a] * P[0];
        g.J21 += g.dRdu[a] * P[1];
        g.J12 += g.dRdv[a] * P[0];
        g.J22 += g.dRdv[a] * P[1];
    }
    double detJ = g.J11*g.J22 - g.J12*g.J21;
    g.detJ = (std::abs(detJ) < 1.e-14) ? 0.0 : detJ;

    return g;
}

Eigen::MatrixXd PatchIntegrator::computeLocalStiffnessContributionNURBS(
    const Patch& patch, const SpanGauss1D& sg_u, const SpanGauss1D& sg_v,
    const std::vector<int>& span, const std::vector<double>& weights)
{
    int ngauss_u = static_cast<int>(sg_u.u_param.size());
    int ngauss_v = static_cast<int>(sg_v.u_param.size());

    std::vector<const double*> pts = patch.control_points_for_span(span);
    size_t nb_loc = pts.size();

    Eigen::MatrixXd K_loc = Eigen::MatrixXd::Zero(2*nb_loc, 2*nb_loc);

    for (int gu = 0; gu < ngauss_u; ++gu) {
        for (int gv = 0; gv < ngauss_v; ++gv) {
            double w = sg_u.weight[gu] * sg_v.weight[gv];

            const Eigen::VectorXd& Nu = sg_u.N[gu];
            const Eigen::VectorXd& dNu = sg_u.dN[gu];
            const Eigen::VectorXd& Nv = sg_v.N[gv];
            const Eigen::VectorXd& dNv = sg_v.dN[gv];

            GaussPointGeometry g = evaluateGaussPointGeometryNURBS(pts, Nu, dNu, Nv, dNv, weights);
            if (g.detJ == 0.0) continue;

            const std::vector<double>& dRdu = g.dRdu;
            const std::vector<double>& dRdv = g.dRdv;
            double detJ = g.detJ;

            double invJ11 = g.J22 / detJ;
            double invJ12 = -g.J12 / detJ;
            double invJ21 = -g.J21 / detJ;
            double invJ22 = g.J11 / detJ;

            std::vector<std::array<double, 2>> grads(nb_loc);
            for (size_t a = 0; a < nb_loc; ++a) {
                grads[a][0] = invJ11 * dRdu[a] + invJ21 * dRdv[a];
                grads[a][1] = invJ12 * dRdu[a] + invJ22 * dRdv[a];
            }

            Eigen::MatrixXd B = Eigen::MatrixXd::Zero(3, 2 * nb_loc);
            for (size_t a = 0; a < nb_loc; ++a) {
                B(0, 2*a)     = grads[a][0];
                B(1, 2*a + 1) = grads[a][1];
                B(2, 2*a)     = grads[a][1];
                B(2, 2*a + 1) = grads[a][0];
            }

            Eigen::Matrix3d D;
            double E = material_properties_.E;
            double nu = material_properties_.nu;
            double factor = E / (1.0 - nu*nu);
            D << factor, factor*nu, 0.0,
                 factor*nu, factor, 0.0,
                 0.0, 0.0, factor*(1.0-nu)/2.0;

            K_loc += B.transpose() * D * B * w * std::abs(detJ);
        }
    }
    return K_loc;
}

Eigen::MatrixXd PatchIntegrator::computeLocalMassContributionNURBS(
    const Patch& patch, const SpanGauss1D& sg_u, const SpanGauss1D& sg_v,
    const std::vector<int>& span, const std::vector<double>& weights)
{
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

            GaussPointGeometry g = evaluateGaussPointGeometryNURBS(pts, Nu, dNu, Nv, dNv, weights);
            if (g.detJ == 0.0) continue;

            Eigen::MatrixXd N = Eigen::MatrixXd::Zero(2, 2 * nb_loc);
            for (size_t a = 0; a < nb_loc; ++a) {
                N(0, 2*a)     = g.R[a];
                N(1, 2*a + 1) = g.R[a];
            }

            M_loc += rho * N.transpose() * N * w * std::abs(g.detJ);
        }
    }
    return M_loc;
}

Eigen::VectorXd PatchIntegrator::computeLocalBoundaryLoadContributionNURBS(
    const Patch& patch, const std::vector<int>& span, int direction,
    const SpanGauss1D& sg_varying, const Eigen::VectorXd& N_fixed_boundary,
    const Traction& traction, const std::vector<double>& weights)
{
    std::vector<const double*> pts = patch.control_points_for_span(span);
    size_t nb_loc = pts.size();
    int ngauss = static_cast<int>(sg_varying.u_param.size());

    Eigen::VectorXd F_loc = Eigen::VectorXd::Zero(2 * nb_loc);
    Eigen::VectorXd zero_deriv = Eigen::VectorXd::Zero(N_fixed_boundary.size());

    for (int g = 0; g < ngauss; ++g) {
        double w = sg_varying.weight[g];

        const Eigen::VectorXd& Nu  = (direction == 0) ? N_fixed_boundary : sg_varying.N[g];
        const Eigen::VectorXd& dNu = (direction == 0) ? zero_deriv       : sg_varying.dN[g];
        const Eigen::VectorXd& Nv  = (direction == 1) ? N_fixed_boundary : sg_varying.N[g];
        const Eigen::VectorXd& dNv = (direction == 1) ? zero_deriv       : sg_varying.dN[g];

        std::vector<double> R(nb_loc), dRdvar(nb_loc);
        size_t idx = 0;
        for (size_t jv = 0; jv < static_cast<size_t>(Nv.size()); ++jv) {
            for (size_t iu = 0; iu < static_cast<size_t>(Nu.size()); ++iu) {
                R[idx]      = Nu(iu) * Nv(jv);
                dRdvar[idx] = dNu(iu) * Nv(jv) + Nu(iu) * dNv(jv);
                ++idx;
            }
        }

        // NURBS rationalization along the boundary curve
        double W = 0.0, dWdvar = 0.0;
        for (size_t a = 0; a < nb_loc; ++a) {
            W      += weights[a] * R[a];
            dWdvar += weights[a] * dRdvar[a];
        }
        const double inv_W = 1.0 / W;
        for (size_t a = 0; a < nb_loc; ++a) {
            const double Ra = weights[a] * R[a] * inv_W;
            dRdvar[a] = (weights[a] * dRdvar[a] - Ra * dWdvar) * inv_W;
            R[a] = Ra;
        }

        Eigen::Vector2d T = Eigen::Vector2d::Zero();
        Eigen::Vector2d X = Eigen::Vector2d::Zero();
        for (size_t a = 0; a < nb_loc; ++a) {
            const double* P = pts[a];
            T[0] += dRdvar[a] * P[0];
            T[1] += dRdvar[a] * P[1];
            X[0] += R[a] * P[0];
            X[1] += R[a] * P[1];
        }
        double ds = T.norm();
        if (ds < 1.e-14) continue;

        Eigen::Vector2d t = traction.evaluate(X);

        for (size_t a = 0; a < nb_loc; ++a) {
            F_loc[2*a]     += R[a] * t[0] * w * ds;
            F_loc[2*a + 1] += R[a] * t[1] * w * ds;
        }
    }

    return F_loc;
}
double PatchIntegrator::integrateScalarOperator(ScalarLocalOperator& op,
                                                const Eigen::VectorXd& u_global)
{
    const bool rational = patch_.cp_manager->is_rational();
    const size_t dpc = patch_.dof_manager->dofs_per_control_point;
    double result = 0.0;

    SpanNDIterator it = patch_.spans();
    for (auto span : it) {
        int idx_u = basis_u_.span_indices.at(span[0]);
        int idx_v = basis_v_.span_indices.at(span[1]);
        const SpanGauss1D& sg_u = basis_u_.gauss_spans[idx_u];
        const SpanGauss1D& sg_v = basis_v_.gauss_spans[idx_v];

        std::vector<const double*> pts = patch_.control_points_for_span(span);
        size_t nb_loc = pts.size();

        std::vector<double> weights_span;
        if (rational) weights_span = patch_.weights_for_span(span);

        std::vector<size_t> span_local = buildSpanLocalIndices(patch_, span);

        int ngauss_u = static_cast<int>(sg_u.u_param.size());
        int ngauss_v = static_cast<int>(sg_v.u_param.size());

        for (int gu = 0; gu < ngauss_u; ++gu) {
            for (int gv = 0; gv < ngauss_v; ++gv) {
                double wg = sg_u.weight[gu] * sg_v.weight[gv];

                const Eigen::VectorXd& Nu  = sg_u.N[gu];
                const Eigen::VectorXd& dNu = sg_u.dN[gu];
                const Eigen::VectorXd& Nv  = sg_v.N[gv];
                const Eigen::VectorXd& dNv = sg_v.dN[gv];

                GaussPointGeometry g = rational
                    ? evaluateGaussPointGeometryNURBS(pts, Nu, dNu, Nv, dNv, weights_span)
                    : evaluateGaussPointGeometry(pts, Nu, dNu, Nv, dNv);

                if (g.detJ == 0.0) continue;

                double invJ11 =  g.J22 / g.detJ;
                double invJ12 = -g.J12 / g.detJ;
                double invJ21 = -g.J21 / g.detJ;
                double invJ22 =  g.J11 / g.detJ;

                std::vector<double> dRdx(nb_loc), dRdy(nb_loc);
                Eigen::Vector2d x_gp = Eigen::Vector2d::Zero();
                for (size_t a = 0; a < nb_loc; ++a) {
                    dRdx[a] = invJ11 * g.dRdu[a] + invJ21 * g.dRdv[a];
                    dRdy[a] = invJ12 * g.dRdu[a] + invJ22 * g.dRdv[a];
                    x_gp[0] += g.R[a] * pts[a][0];
                    x_gp[1] += g.R[a] * pts[a][1];
                }

                // Gather local DOF values: layout [u_x0, u_y0, u_x1, u_y1, ...]
                Eigen::VectorXd u_local(static_cast<Eigen::Index>(nb_loc * dpc));
                for (size_t a = 0; a < nb_loc; ++a) {
                    for (size_t d = 0; d < dpc; ++d) {
                        size_t global_dof = patch_.dof_manager->get_global_dof(span_local[a], d);
                        u_local[static_cast<Eigen::Index>(a * dpc + d)] = u_global[static_cast<Eigen::Index>(global_dof)];
                    }
                }

                double integrand = op.computeScalarIntegrand(g.R, dRdx, dRdy, x_gp, u_local);
                result += integrand * wg * std::abs(g.detJ);
            }
        }
    }

    return result;
}
