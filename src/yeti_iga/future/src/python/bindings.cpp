#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include "BSpline.hpp"
#include "BSplineTensor.hpp"
#include "ControlPointManager.hpp"
#include "Patch.hpp"


namespace py = pybind11;

PYBIND11_MODULE(bspline, m)
{
    py::class_<BSpline>(m, "BSpline")
        .def(py::init<int, py::array_t<double>>())
        .def("find_span", &BSpline::FindSpan)
        .def("basis_funs", &BSpline::BasisFuns)
        .def("one_basis_fun", &BSpline::OneBasisFun)
        .def_property_readonly("degree", &BSpline::getDegree)
        .def_property_readonly("knot_vector", &BSpline::kvView);

    py::class_<BSplineTensor>(m, "BSplineTensor")
        .def_property_readonly("components", [](const BSplineTensor& t) {return t.components; })
        .def("basis_funs_nd", &BSplineTensor::BasisFunsND)
        .def("find_span_nd",
             static_cast<py::array_t<int>(BSplineTensor::*)(const py::array_t<double>&) const>
                (&BSplineTensor::FindSpanND),
             py::arg("u"),
             "Compute span indices in all parameter dimensions.");

    py::class_<BSplineSurface, BSplineTensor>(m, "BSplineSurface")
        .def(py::init<const BSpline&, const BSpline&>(),
             py::arg("su"), py::arg("sv"));

    py::class_<BSplineVolume, BSplineTensor>(m, "BSplineVolume")
        .def(py::init<const BSpline&, const BSpline&, const BSpline&>(),
             py::arg("su"), py::arg("sv"), py::arg("sw"));

    py::class_<ControlPointManager>(m, "ControlPointManager")
        .def(py::init<int>(), py::arg("dim")=3)
        .def_property_readonly("dim_phys", [](const ControlPointManager& mgr) {return mgr.dim_phys; })
        .def("add_point", &ControlPointManager::add_point)
        .def_property_readonly("n_points", &ControlPointManager::n_points)
        .def("coords_view", [](ControlPointManager& self){
            auto capsule = py::capsule(&self);
            std::vector<ssize_t> shape = {(ssize_t)self.n_points(), (ssize_t)self.dim_phys};
            std::vector<ssize_t> strides = {
                static_cast<std::ptrdiff_t>(self.dim_phys * sizeof(double)),
                static_cast<std::ptrdiff_t>(sizeof(double))
            };
            return py::array_t<double>(shape, strides, self.coords.data(), capsule);
        });

    py::class_<Patch>(m, "Patch")
        .def(py::init<const BSplineTensor&, ControlPointManager*, const std::vector<size_t>&, const std::vector<size_t>&>(),
             py::arg("tensor"), py::arg("cp_manager"), py::arg("global_indices"), py::arg("local_shape"))
        .def("local_cp_ptr", &Patch::local_cp_ptr,
             py::arg("i_local"),
             "Return pointer to local control point (as int or PyCapsule for Python?)")
        .def("local_cp_view", &Patch::local_cp_view,
            py::arg("i_local"),
            "Return Numpy 1D view on a local point")
        .def("local_control_point_view", &Patch::local_control_point_view,
             "Return a zero-copy NumPy array view of global control points (indexing from Python needed)")
        .def("evaluate_patch_nd", &Patch::EvaluatePatchND)
        .def("evaluate_patch_nd_omp", &Patch::EvaluatePatchNDOMP);
}
