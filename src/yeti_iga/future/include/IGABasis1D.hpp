#pragma once
#include <vector>
#include <cassert>
#include <Eigen/Dense>
#include "GaussQuadrature.hpp"
#include "BSpline.hpp"

// Precomputed basis values for one span at gauss points
struct SpanGauss1D {
    std::vector<double> u_param;    // Gauss point in BSpline parametric space
    std::vector<double> weight;     // Gauss integration wieght

    std::vector<Eigen::VectorXd> N;  // functions values at Gauss points
    std::vector<Eigen::VectorXd> dN; // functions derivatives at Gauss points
    std::vector<Eigen::VectorXd> d2N;// 2nd derivatives at Gauss points; empty unless
                                      // IGABasis1D::build() was called with deriv_order>=2
                                      // (e.g. Kirchhoff-Love shell curvature terms)
};


// Object handling precomputed values for all spans of a 1D BSpline parametric space
struct IGABasis1D {
    std::vector<SpanGauss1D> gauss_spans;
    std::unordered_map<int, int> span_indices;


    // Build from a BSpline object
    // - gauss_n : number of Gauss points per span
    // - deriv_order : highest basis-function derivative order to precompute (1 by
    //   default, matching every existing solid/plane-stress caller byte-for-byte).
    //   Pass 2 to additionally fill SpanGauss1D::d2N (needed by shell curvature terms).
    static IGABasis1D build(const BSpline& bspline, int ngauss, int deriv_order = 1) {
        IGABasis1D out;
        const auto& kv = bspline.getKnotVector();
        int p = bspline.getDegree();
        int m = static_cast<int>(kv.size()) - 1;

        std::vector<double> gauss_points, gauss_weights;
        gauss_legendre_table(ngauss, gauss_points, gauss_weights);

        // find valid spans and create indices map
        int span_count = 0;
        for (int span = p; span <= m - p - 1; ++span) {
            if (!(kv[span+1] > kv[span])) continue;
            out.span_indices[span] = span_count;
            span_count++;

            SpanGauss1D sg;
            sg.u_param.reserve(ngauss);
            sg.weight.reserve(ngauss);
            sg.N.reserve(ngauss);
            sg.dN.reserve(ngauss);
            if (deriv_order >= 2) sg.d2N.reserve(ngauss);

            double a = kv[span];
            double b = kv[span+1];
            double half = 0.5 * (b - a);
            double mid = 0.5 * (b + a);

            for (int i = 0; i < ngauss; ++i) {
                double xi = gauss_points[i];
                double up = mid + half*xi;
                sg.u_param.push_back(up);
                sg.weight.push_back(gauss_weights[i] * half);    // weight * jacobian of mapping

                // compute functions + derivatives up to deriv_order (output in ders,
                // one row of (p+1) values per derivative order 0..deriv_order).
                // Zero-initialized: if deriv_order > p, BasisFunsDerivatives caps its
                // internal order at p and leaves higher rows untouched, which is
                // mathematically correct (a degree-p basis has zero derivatives beyond
                // order p) as long as the buffer starts at zero.
                std::vector<double> ders((deriv_order+1)*(p+1), 0.0);
                bspline.BasisFunsDerivatives(span, up, deriv_order, ders.data());

                // Extract functions and derivatives
                Eigen::VectorXd N(p+1);
                Eigen::VectorXd dN(p+1);
                for (int j = 0; j <= p; ++j) {
                    N(j) = ders[j];
                    dN(j) = ders[(p+1) + j];
                }

                sg.N.push_back(std::move(N));
                sg.dN.push_back(std::move(dN));

                if (deriv_order >= 2) {
                    Eigen::VectorXd d2N(p+1);
                    for (int j = 0; j <= p; ++j)
                        d2N(j) = ders[2*(p+1) + j];
                    sg.d2N.push_back(std::move(d2N));
                }
            }
            out.gauss_spans.push_back(std::move(sg));
        }
        return out;
    }
};