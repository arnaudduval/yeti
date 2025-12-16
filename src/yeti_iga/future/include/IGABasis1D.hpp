#pragma once
#include <vector>
#include <cassert>
#include "GaussQuadrature.hpp"
#include "BSpline.hpp"

// Precomputed basis values for one span at gauss points
struct SpanGauss1D {
    std::vector<double> u_param;
    std::vector<double> weight;
    // Functions et derivatives values N, dN
    // N[g][a] : g in 0..ngauss-1, a in 0..p
    std::vector<std::vector<double>> N;
    std::vector<std::vector<double>> dN;
};

struct IGABasis1D {
    std::vector<SpanGauss1D> spans;

    // Build fram a BSpline object
    // - gauss_n : number of Gauss points per span
    // - derivative_order = 1 (1st derivative only)
    static IGABasis1D build(const BSpline& b, int gauss_n) {
        IGABasis1D out;
        const auto& kv = b.getKnotVector();
        int p = b.getDegree();
        int m = static_cast<int>(kv.size()) - 1;

        // find valid spans: i in [p .. m-p-1] where kv[i] < kv[i+1]
        for (int i = p; i <= m - p - 1; ++i) {
            if (!(kv[i+1] > kv[i])) continue;
            SpanGauss1D sg;
            std::vector<double> gx, gw;
            gauss_legendre_table(gauss_n, gx, gw);

            sg.u_param.resize(gauss_n);
            sg.weight.resize(gauss_n);
            sg.N.assign(gauss_n, std::vector<double>(p+1));
            sg.dN.assign(gauss_n, std::vector<double>(p+1));

            double a = kv[i];
            double bval = kv[i+1];
            double half = 0.5 * (bval - a);
            double mid = 0.5 * (bval + a);

            std::vector<double> ders((i+1)*(p+1));

            for (int g = 0; g < gauss_n; ++g) {
                double xi = gx[g];
                double up = mid + half*xi;
                sg.u_param[g] = up;
                sg.weight[g] = gw[g] * half;    // weight * jacobian of mapping

                // compute functions + 1st derivative (output in ders)
                b.BasisFunsDerivatives(i, up, 1, ders.data());
                // ders layout: ders[k*(p+1) + a], k=0..1
                for (int aidx = 0; aidx <= p; ++aidx) {

                    sg.N[g][aidx] = ders[aidx];   // == ders[0*(p+1) + aidx];
                    sg.dN[g][aidx] = ders[1*(p+1) + aidx];
                }
            }
            out.spans.push_back(std::move(sg));
        }
        return out;
    }
};