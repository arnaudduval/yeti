#include "IGAAssembler.hpp"
#include <cstring>
#include <cmath>
#include <iostream>

std::vector<ElementMatrix> IGAAssembler2D::assemble_stiffness() const {
    std::vector<ElementMatrix> elems;

    // Get iterator over spans
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

        int ngauss_u = static_cast<int>(sg_u.u_phys.size());
        int ngauss_v = static_cast<int>(sg_v.u_phys.size());

        // Elementary matrix
        ElementMatrix E;
        E.nb_loc = nb_loc_;
        E.K.assign(nb_loc_ * nb_loc_, 0.0);

        // Get active control points pointers for current span
        std::vector<const double*> pts = patch_.control_points_for_span(span);
        if (pts.size() != (size_t)nb_loc_) {
            std::cerr << "WARNING: pts.size() != nb_loc_ for span (" << span_u << ", " << span_v << "\n";
            // continue to next span
            continue;
        }

        // Global indices for local DOFs
        E.global_indices = patch_.global_indices;
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
                size_t local_linear = static_cast<size_t>(lv * n_u + iu);
                size_t gid = patch_.global_indices.at(local_linear);
                local_to_global[idx++] = gid;
            }
        }
        E.global_indices = local_to_global;

        // Loop over Gauss points
        for (int gu = 0; gu < ngauss_u; ++gu) {
            for(int gv = 0; gv < ngauss_v; ++gv) {
                double wu = sg_u.w_phys[gu];
                double wv = sg_v.w_phys[gv];
                double w = wu * wv;

                // Build arrays of 1D N and dN pointers
                const std::vector<double>& Nu = sg_u.N[gu];
                const std::vector<double>& dNu = sg_u.dN[gu];
                const std::vector<double>& Nv = sg_v.N[gv];
                const std::vector<double>& dNv = sg_v.dN[gv];

                // Compute geometry mapping
                double J11 = 0.0, J12 = 0.0, J21 = 0.0, J22 = 0.0;
                double Xx = 0.0, Yy = 0.0;

                idx = 0;
                for (int jv = 0; jv <= p_v_; ++jv) {
                    for (int iu = 0; iu <= p_u_; ++iu) {
                        double R = Nu[iu] * Nv[jv];
                        double dRdu = dNu[iu] * Nv[jv];
                        double dRdv = Nu[iu] * dNv[jv];
                        const double* P = pts[idx]; // Pointer to coordinates x, y
                        double px = P[0], py = P[1];
                        Xx += R * px;
                        Yy += R * py;
                        J11 += dRdu * px;   // dX/du
                        J21 += dRdu * py;   // dY/du
                        J12 += dRdv * px;   // dX/dv
                        J22 += dRdv * py;   // dY/dv
                        ++idx;
                    }
                }
                double detJ = J11 * J22 - J12 * J21;
                if (std::abs(detJ) < 1.e-14) {
                    // Degenerate mapping : handle exception, skip or warning ?
                    continue;
                }

                double invJ11 =  J22 / detJ;
                double invJ12 = -J12 / detJ;
                double invJ21 = -J21 / detJ;
                double invJ22 =  J11 / detJ;

                // Local gradient of basis function in physical coords
                // grad R_a = [dR/dx, dR/dy] = invJ^T * [dR/du, dR/dv]
                std::vector<std::array<double, 2>> grads(nb_loc_);
                idx = 0;
                for (int jv = 0; jv <= p_v_; ++jv) {
                    for (int iu = 0; iu <= p_u_; ++iu) {
                        // TODO : on a déjà calculé ces quantités
                        double dRdu = dNu[iu] * Nv[jv];
                        double dRdv = Nu[iu] * dNv[jv];
                        // [dR/dx] = invJ11 * dRdu + invJ21 * dRdv
                        // [dR/dy] = invJ12 * dRdu + invJ22 * dRdv
                        double dRx = invJ11 * dRdu + invJ21 * dRdv;
                        double dRy = invJ12 * dRdu + invJ22 * dRdv;
                        grads[idx][0] = dRx;
                        grads[idx][1] = dRy;
                        ++idx;
                    }
                }

                // Fill K_loc: K_ab += (grad_a . grad_b) * w * detJ
                double weight = w * std::abs(detJ);
                for (size_t a = 0; a < (size_t)nb_loc_; ++a) {
                    for (size_t b = 0; b < (size_t)nb_loc_; ++b) {
                        double dot = grads[a][0]*grads[b][0] + grads[a][1]*grads[b][1];
                        E.K[a*nb_loc_ + b] += dot * weight;
                    }
                }
            }
        }

        elems.push_back(std::move(E));

    }
    return elems;
}
