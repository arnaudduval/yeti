#pragma once
#include <vector>
#include <Eigen/Dense>

// Base class for a custom integration term, meant to be subclassed from
// Python via the LocalOperator pybind11 trampoline (see bindings.cpp) for
// development/testing convenience. PatchIntegrator does the Gauss-point
// loop, the Jacobian/geometry computation and the weighting/summation --
// the operator only has to return the term to integrate at one Gauss point
// (e.g. B^T*D*B for stiffness, rho*N^T*N for mass), built from the basis
// values R and physical gradients dRdx, dRdy already computed by
// PatchIntegrator. No C++ rebuild needed to try a new physical kernel.
class LocalOperator {
public:
    virtual ~LocalOperator() = default;

    // R, dRdx, dRdy: basis values and physical-space gradients of this
    // span's active basis functions at one Gauss point (u-fastest order,
    // matching patch.control_points_for_span()'s row order). The returned
    // matrix is multiplied by the Gauss weight and |detJ| and summed over
    // every Gauss point of every span by PatchIntegrator -- do not apply
    // weight/detJ here.
    virtual Eigen::MatrixXd computeIntegrand(
        const std::vector<double>& R,
        const std::vector<double>& dRdx,
        const std::vector<double>& dRdy) const = 0;
};
