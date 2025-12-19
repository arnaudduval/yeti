#include "IGAAssembler.hpp"
#include <cstring>
#include <cmath>
#include <iostream>

Eigen::Matrix3d IGAAssembler2D::computeConstitutiveMatrix() const {
    double E = material_.E;
    double nu = material_.nu;
    double factor = E / (1.0 - nu*nu);
    Eigen::Matrix3d D;
    D <<
        factor, factor * nu, 0.0,
        factor * nu, factor, 0.0,
        0.0, 0.0, factor * (1.0 - nu) / 2.0;
    return D;
}



std::vector<ElementMatrix> IGAAssembler2D::assemble_stiffness() const {
    std::vector<ElementMatrix> elems;
    SpanNDIterator it(patch_.tensor);

    for (auto span : it) {
        int span_u = span[0];
        int span_v = span[1];

        // Find corresponding span in basis pre-computation
        // To find index: need to know the list of valid spans
        // Build vector of valid spans in each direction, by comparing knots locations
        const auto& kv_u = patch_.tensor.components[0].getKnotVector();
        const auto& kv_v = patch_.tensor.components[1].getKnotVector();

        // Determine span index in precomputed arrays
        int idx_u = -1;
        {
            int p = p_u_;
            int m = static_cast<int>(kv_u.size()) - 1;
            int count = 0;
            for (int i = p; i<= m - p - 1; ++i) {
                if (!(kv_u[i+1] > kv_u[i])) continue;
                if (i == span_u) {idx_u = count; break; }
                ++count;
            }
        }
        int idx_v = -1;
        {
            int p = p_v_;
            int m = static_cast<int>(kv_v.size()) - 1;
            int count = 0;
            for (int i = p; i <= m - p - 1; ++i) {
                if (!(kv_v[i+1] > kv_v[i])) continue;
                if (i == span_v) { idx_v = count; break; }
                ++count;
            }
        }
        if (idx_u < 0 || idx_v < 0) {
            // shoud not happend => raise an exception ?
            continue;
        }

        const SpanGauss1D& sg_u = basis_u_.spans[idx_u];
        const SpanGauss1D& sg_v = basis_v_.spans[idx_v];

        int ngauss_u = static_cast<int>(sg_u.u_param.size());
        int ngauss_v = static_cast<int>(sg_v.u_param.size());

        // Elementary matrix
        ElementMatrix E(nb_loc_);
        Eigen::Matrix<double, 12, 12> K_loc = Eigen::Matrix<double, 12, 12>::Zero();
        Eigen::Matrix3d D = computeConstitutiveMatrix();

        // Get active control points pointers for current span
        std::vector<const double*> pts = patch_.control_points_for_span(span);
        if (pts.size() != (size_t)nb_loc_) {
            std::cerr << "WARNING: pts.size() != nb_loc_ for span (" << span_u << ", " << span_v << "\n";
            // continue to next span
            continue;
        }

        // We have the global mapping for patch
        // We must build local mapping for current span
        // start_u = span_u + p_u ...
        int start_u = span_u - p_u_;
        int start_v = span_v - p_v_;
        ssize_t n_u = patch_.local_shape[0];
        ssize_t n_v = patch_.local_shape[1];

        std::vector<size_t> local_to_global(nb_loc_);
        // ordering : i_u + i_v * n_u
        int idx = 0;
        for (int jv = 0; jv <= p_v_; ++jv) {
            int lv = start_v + jv;
            for (int iu = 0; iu <= p_u_; ++iu) {
                int lu = start_u + iu;
                size_t local_linear = static_cast<size_t>(lv * n_u + lu);
                size_t gid = patch_.global_indices.at(local_linear);
                local_to_global[idx++] = gid;
            }
        }
        E.global_indices = local_to_global;

        // Loop over Gauss points
        for (int gu = 0; gu < ngauss_u; ++gu) {
            for (int gv = 0; gv < ngauss_v; ++gv) {
                double w = sg_u.weight[gu] * sg_v.weight[gv];

                // Get pre-computed basis functions and their derivatives
                const std::vector<double>& Nu = sg_u.N[gu];
                const std::vector<double>& dNu = sg_u.dN[gu];
                const std::vector<double>& Nv = sg_v.N[gv];
                const std::vector<double>& dNv = sg_v.dN[gv];

                // Compute R, dR/du, dR/dv
                std::vector<double> R(nb_loc_);
                std::vector<double> dRdu(nb_loc_);
                std::vector<double> dRdv(nb_loc_);

                idx = 0;
                for (int jv = 0; jv <= p_v_; ++jv) {
                    for (int iu = 0; iu <= p_u_; ++iu) {
                        R[idx] = Nu[iu] * Nv[jv];
                        dRdu[idx] = dNu[iu] * Nv[jv];
                        dRdv[idx] = Nu[iu] * dNv[jv];
                        ++idx;
                    }
                }

                // Compute geometry mapping
                double J11 = 0.0, J12 = 0.0, J21 = 0.0, J22 = 0.0;
                double Xx = 0.0, Yy = 0.0;

                for (int a = 0; a < nb_loc_; ++a) {
                    const double* P = pts[a];
                    const double px = P[0];
                    const double py = P[1];

                    Xx += R[a] * px;
                    Yy += R[a] * py;

                    J11 += dRdu[a] * px;
                    J21 += dRdu[a] * py;
                    J12 += dRdv[a] * px;
                    J22 += dRdv[a] * py;
                }

                double detJ = J11 * J22 - J12 * J21;
                if (std::abs(detJ) < 1.e-14) {
                    // Degenerate mapping : handle exception, skip or warning ?
                    continue;
                }

                double invJ11 = J22 / detJ;
                double invJ12 = -J12 / detJ;
                double invJ21 = -J21 / detJ;
                double invJ22 = J11 / detJ;


                // Local gradient of basis function wrto physical coords
                // grad R_a = [dR/dx, dR/dy] = invJ^T * [dR/du, dR/dv]
                std::vector<std::array<double, 2>> grads(nb_loc_);
                for (int a = 0; a < nb_loc_; ++a) {
                    grads[a][0] = invJ11 * dRdu[a] + invJ21 * dRdv[a]; // dR/dx
                    grads[a][1] = invJ12 * dRdu[a] + invJ22 * dRdv[a]; // dR/dy
                }
                for (int a = 0; a < nb_loc_; ++a) {
                    Eigen::Matrix<double, 3, 12> B;
                    B.setZero();
                    for (int a = 0; a < nb_loc_; ++a) {
                        B(0, 2*a)   = grads[a][0];  // dN/dx for u_x
                        B(1, 2*a+1) = grads[a][1];  // dN/dy for u_y
                        B(2, 2*a)   = grads[a][1];  // dN/dy for shear (u_x)
                        B(2, 2*a+1) = grads[a][0];  // dN/dx for shear (u_y)
                    }

                    // Local contribution : B^T * D * B * w * detJ
                    K_loc += B.transpose() * D * B * w * std::abs(detJ);
                }
            }
        }
        E.K = K_loc;
        elems.push_back(std::move(E));
    }
    return elems;
}
