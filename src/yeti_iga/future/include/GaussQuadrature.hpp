#pragma once
#include <vector>
#include <cassert>
#include <stdexcept>

// To be completed for higher orders
inline void gauss_legendre_table(int ngauss, std::vector<double>& points, std::vector<double>& weights) {
    points.resize(ngauss);
    weights.resize(ngauss);

    if (ngauss == 1) {
        points[0] = 0.0;
        weights[0] = 2.0;
    } else if (ngauss == 2) {
        points[0] = -0.5773502691896257;
        points[1] = -points[0];
        weights[0] = weights[1] = 1.0;
    } else if (ngauss == 3) {
        points[0] = -0.7745966692414834;
        points[1] = 0.0;
        points[2] = -points[0];
        weights[0] = weights[2] = 0.5555555555555556;
        weights[1] = 0.8888888888888888;
    } else if (ngauss == 4) {
        points[0] = -0.8611363115940526;
        points[1] = -0.3399810435848563;
        points[2] = -points[1];
        points[3] = -points[0];
        weights[0] = weights[3] = 0.3478548451374538;
        weights[1] = weights[2] = 0.6521451548625461;
    } else if (ngauss == 5) {
        points[0] = -0.9061798459386640;
        points[1] = -0.5384693101056831;
        points[2] = 0.0;
        points[3] = -points[1];
        points[4] = -points[0];
        weights[0] = weights[4] = 0.2369268850561891;
        weights[1] = weights[3] = 0.4786286704993665;
        weights[2] = 0.5688888888888889;
    } else if (ngauss == 6) {
        points[0] = -0.9324695142031521;
        points[1] = -0.6612093864662645;
        points[2] = -0.2386191860831969;
        points[3] = -points[2];
        points[4] = -points[1];
        points[5] = -points[0];
        weights[0] = weights[5] = 0.1713244923791704;
        weights[1] = weights[4] = 0.3607615730481386;
        weights[2] = weights[3] = 0.4679139345726910;
    }
    else {
        throw std::runtime_error("Number of gauss points not supported");
    }

}
