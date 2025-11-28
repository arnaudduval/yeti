#pragma once
#include <vector>
#include <pybind11/numpy.h>
#include "BSpline.hpp"
#include "BSplineTensor.hpp"
#include "ControlPointManager.hpp"

namespace py = pybind11;

struct Patch {
    BSplineTensor tensor;

    // mapping local -> global
    std::vector<size_t> global_indices;

    // dimensions nu, nv, nw (for indexing)
    std::vector<size_t> local_shape;

    ControlPointManager* cp_manager = nullptr;

    Patch(const BSplineTensor& t,
          ControlPointManager* mgr,
          const std::vector<size_t>& mapping,
          const std::vector<size_t>& local_shape_)
        : tensor(t), global_indices(mapping), local_shape(local_shape_), cp_manager(mgr) {}

    // get view (zero copy) to coordinates of local control point i_local
    double* local_cp_ptr(size_t i_local);

    // Numpy 1D view on a local point
    py::array_t<double> local_cp_view(size_t i_local);

    // get a numpy view for all local control points
    // shape = local_shape + (dim_phys,)
    py::array_t<double> local_control_point_view() const;

    py::array_t<double> EvaluatePatchND(const py::array_t<int> spans,
                                        const py::array_t<double>& u) const;

    py::array_t<double> EvaluatePatchNDOMP(const py::array_t<int> spans,
                                        const py::array_t<double>& u) const;
};