#pragma once
#include <vector>
#include <Eigen/Dense>

// Base class for a scalar-valued integration term (error norms, energy
// functionals, etc.), designed to be subclassed from Python via the
// ScalarLocalOperator pybind11 trampoline (see bindings.cpp).
//
// PatchIntegrator::integrateScalarOperator() provides the Gauss loop,
// Jacobian, physical-space gradients, physical coordinates, and the
// LOCAL DOF values extracted from the global solution vector.
// The subclass only implements the scalar integrand at one Gauss point.
//
// Arguments passed to computeScalarIntegrand:
//   R, dRdx, dRdy  -- basis values and physical-space gradients of the
//                      span's active functions (u-fastest order, same as
//                      LocalOperator::computeIntegrand).
//   physical_point  -- (x, y) physical coordinates of the Gauss point.
//   u_local         -- local DOF values in dof_manager order:
//                      [comp0_cp0, comp1_cp0, comp0_cp1, comp1_cp1, ...]
//                      (size = n_active_CPs * dofs_per_CP).
//                      For a 2-DOF 2-D problem: u_local[2a] = u_x of CP a,
//                      u_local[2a+1] = u_y of CP a.
//
// PatchIntegrator multiplies the returned scalar by the Gauss weight and
// |detJ|, then sums it over every Gauss point of every span.
class ScalarLocalOperator {
public:
    virtual ~ScalarLocalOperator() = default;

    virtual double computeScalarIntegrand(
        const std::vector<double>& R,
        const std::vector<double>& dRdx,
        const std::vector<double>& dRdy,
        const Eigen::Vector2d& physical_point,
        const Eigen::VectorXd& u_local) const = 0;
};
