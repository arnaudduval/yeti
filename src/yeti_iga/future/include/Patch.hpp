#pragma once
#include <vector>
#include <pybind11/numpy.h>
#include "BSpline.hpp"

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
    double* local_cp_ptr(size_t i_local) {
        size_t gid = global_indices.at(i_local);
        return cp_manager->coords.data() + gid * cp_manager->dim_phys;
    }

    // Numpy 1D view on a local point
    py::array_t<double> local_cp_view(size_t i_local) {
        double* ptr = local_cp_ptr(i_local);
        auto capsule = py::capsule(ptr); // do not deallocate
        return py::array_t<double>(
            {cp_manager->dim_phys},       // shape 1D
            {sizeof(double)},             // stride
            ptr,
            capsule
        );
    }

    // get a numpy view for all local control points
    // shape = local_shape + (dim_phys,)
    py::array_t<double> local_control_point_view() const {
        std::vector<ssize_t> shape;
        for (auto s : local_shape) shape.push_back((ssize_t)s);
        shape.push_back((ssize_t)cp_manager->dim_phys);

        // return a view on the *global* coords with shape (n_points_global, dim_phys).
        std::vector<ssize_t> gshape = {(ssize_t)cp_manager->n_points(), (ssize_t)cp_manager->dim_phys};
        auto arr = py::array_t<double>(gshape, { (ssize_t)sizeof(double)*cp_manager->dim_phys, (ssize_t)sizeof(double) }, cp_manager->coords.data());
        // User can index arr[global_indices] in Python to build local array (zero-copy slicing not trivial).
        return arr;
    }
};