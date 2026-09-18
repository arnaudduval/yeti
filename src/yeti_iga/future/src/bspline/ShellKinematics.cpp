#include "ShellKinematics.hpp"

// Direct translation of curvilinearCoordinates.f's `curvilinear` subroutine.
ShellGeometry computeShellGeometry(
    const std::vector<const double*>& pts,
    const std::vector<double>& dRdu, const std::vector<double>& dRdv,
    const std::vector<double>& d2Rdu2, const std::vector<double>& d2Rdv2,
    const std::vector<double>& d2Rdudv)
{
    ShellGeometry g;
    g.a1.setZero(); g.a2.setZero();
    g.a11.setZero(); g.a12.setZero(); g.a22.setZero();

    const size_t nb_loc = pts.size();
    for (size_t a = 0; a < nb_loc; ++a) {
        const Eigen::Map<const Eigen::Vector3d> P(pts[a]);
        g.a1  += dRdu[a]    * P;
        g.a2  += dRdv[a]    * P;
        g.a11 += d2Rdu2[a]  * P;
        g.a22 += d2Rdv2[a]  * P;
        g.a12 += d2Rdudv[a] * P;
    }

    Eigen::Vector3d n = g.a1.cross(g.a2);
    g.Area = n.norm();
    // Area == 0 signals a degenerate Gauss point (a1 parallel to a2, e.g. a
    // collapsed/singular map); callers should skip it, mirroring the solid
    // kernels' `if (detJ == 0.0) continue;` convention. Guard the divide here
    // so a3 doesn't silently become NaN/Inf.
    g.a3 = (g.Area > 1.e-14) ? Eigen::Vector3d(n / g.Area) : Eigen::Vector3d::Zero();

    g.AAI(0,0) = g.a1.dot(g.a1);
    g.AAI(0,1) = g.AAI(1,0) = g.a1.dot(g.a2);
    g.AAI(1,1) = g.a2.dot(g.a2);
    g.AAE = g.AAI.inverse();

    return g;
}

// Direct translation of USFMEM_shell.f's BoJ construction.
Eigen::Matrix3d membraneB(const ShellGeometry& g, double dRdu_a, double dRdv_a) {
    Eigen::Matrix3d B;
    B.row(0) = dRdu_a * g.a1;
    B.row(1) = dRdv_a * g.a2;
    B.row(2) = dRdv_a * g.a1 + dRdu_a * g.a2;
    return B;
}

// Direct translation of USFBND_shell.f's BoJ construction. The cross/dot
// products of a1,a2,a3,a11,a12,a22 are Gauss-point-level (node-independent)
// quantities; recomputing them per node here is a small, deliberate
// simplification versus the Fortran (which precomputes them once per Gauss
// point outside the node loop) -- fine given bendingB's cost is dominated by
// the cross products either way, and PatchIntegrator's caller can hoist this
// if profiling ever shows it matters.
Eigen::Matrix3d bendingB(const ShellGeometry& g, double dRdu_a, double dRdv_a,
                         double d2Rdu2_a, double d2Rdv2_a, double d2Rdudv_a)
{
    const Eigen::Vector3d dA1d1_A2 = g.a11.cross(g.a2);
    const Eigen::Vector3d dA2d2_A2 = g.a22.cross(g.a2);
    const Eigen::Vector3d dA1d2_A2 = g.a12.cross(g.a2);
    const Eigen::Vector3d A1_dA1d1 = g.a1.cross(g.a11);
    const Eigen::Vector3d A1_dA2d2 = g.a1.cross(g.a22);
    const Eigen::Vector3d A1_dA1d2 = g.a1.cross(g.a12);
    const Eigen::Vector3d A2_A3    = g.a2.cross(g.a3);
    const Eigen::Vector3d A3_A1    = g.a3.cross(g.a1);
    const double A3dA1d1 = g.a3.dot(g.a11);
    const double A3dA2d2 = g.a3.dot(g.a22);
    const double A3dA1d2 = g.a3.dot(g.a12);

    const Eigen::Vector3d common = dRdu_a * A2_A3 + dRdv_a * A3_A1;

    const Eigen::Vector3d B1 = -d2Rdu2_a   * g.a3
        + (dRdu_a * dA1d1_A2 + dRdv_a * A1_dA1d1 + A3dA1d1 * common) / g.Area;
    const Eigen::Vector3d B2 = -d2Rdv2_a   * g.a3
        + (dRdu_a * dA2d2_A2 + dRdv_a * A1_dA2d2 + A3dA2d2 * common) / g.Area;
    const Eigen::Vector3d B3 = -d2Rdudv_a  * g.a3
        + (dRdu_a * dA1d2_A2 + dRdv_a * A1_dA1d2 + A3dA1d2 * common) / g.Area;

    Eigen::Matrix3d B;
    B.row(0) = B1;
    B.row(1) = B2;
    B.row(2) = 2.0 * B3;
    return B;
}
