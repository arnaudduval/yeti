#pragma once
#include <pybind11/numpy.h>
#include <Eigen/Dense>
#include "Patch.hpp"

namespace py = pybind11;

// Evaluates a FE solution field at a set of parametric-space points.
//
// Mirrors the EvaluatePatchNDOMP pattern (OpenMP-parallel, one thread-local
// buffer per thread), but instead of summing control-point coordinates it
// sums DOF values from the global solution vector weighted by the NURBS
// rational basis functions.  B-spline patches (is_rational() == false) use
// the same zero-overhead fast path as EvaluatePatchNDOMP.
//
// The patch must have a PatchDOFManager attached at construction.
class PatchEvaluator {
public:
    explicit PatchEvaluator(const Patch& patch);

    // Evaluate the FE solution at parametric-space points.
    //   params   : (n_pts, n_param_dims) array of parameter values
    //   u_global : global solution vector (total DOF count)
    //   returns  : (n_pts, n_dofs_per_cp) array of field values
    py::array_t<double> evaluateSolutionOMP(const py::array_t<double>& params,
                                            const Eigen::VectorXd& u_global) const;

private:
    const Patch& patch_;
};
