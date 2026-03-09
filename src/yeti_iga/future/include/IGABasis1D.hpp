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

    std::vector<Eigen::VectorXd> N; // functions values at Gauss points
    std::vector<Eigen::VectorXd> dN;// functions derivatives at Gauss points
};


// Object handling precomputed values for all spans of a 1D BSpline parametric space
struct IGABasis1D {
    std::vector<SpanGauss1D> gauss_spans;


    // Build from a BSpline object
    // - gauss_n : number of Gauss points per span
    static IGABasis1D build(const BSpline& bspline, int ngauss) {
        IGABasis1D out;
        const auto& kv = bspline.getKnotVector();
        int p = bspline.getDegree();
        int m = static_cast<int>(kv.size()) - 1;

        std::vector<double> gauss_points, gauss_weights;
        gauss_legendre_table(ngauss, gauss_points, gauss_weights);

        // find valid spans: span in [p .. m-p-1] where kv[span] < kv[span+1]
        for (int span = p; span <= m - p - 1; ++span) {
            if (!(kv[span+1] > kv[span])) continue;

            SpanGauss1D sg;
            sg.u_param.reserve(ngauss);
            sg.weight.reserve(ngauss);
            sg.N.reserve(ngauss);
            sg.dN.reserve(ngauss);

            double a = kv[span];
            double b = kv[span+1];
            double half = 0.5 * (b - a);
            double mid = 0.5 * (b + a);

            for (int i = 0; i < ngauss; ++i) {
                double xi = gauss_points[i];
                double up = mid + half*xi;
                sg.u_param.push_back(up);
                sg.weight.push_back(gauss_weights[i] * half);    // weight * jacobian of mapping

                // compute functions + 1st derivative (output in ders)
                std::vector<double> ders(2*(p+1));  // 2 lines : 1 for function, 1 for derivatives
                bspline.BasisFunsDerivatives(span, up, 1, ders.data());

                // Extract functions and derivatives
                Eigen::VectorXd N(p+1);
                Eigen::VectorXd dN(p+1);
                for (int j = 0; j <= p; ++j) {
                    N(j) = ders[j];
                    dN(j) = ders[(p+1) + j];
                }

                sg.N.push_back(std::move(N));
                sg.dN.push_back(std::move(dN));
            }
            out.gauss_spans.push_back(std::move(sg));
        }
        return out;
    }
};