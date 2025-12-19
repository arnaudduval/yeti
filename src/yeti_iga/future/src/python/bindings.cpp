#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include "BSpline.hpp"
#include "BSplineTensor.hpp"
#include "ControlPointManager.hpp"
#include "Patch.hpp"
#include "SpanNDIterator.hpp"
#include "IGAAssembler.hpp"


namespace py = pybind11;

PYBIND11_MODULE(bspline, m)
{
    py::class_<BSpline>(m, "BSpline")
        .def(py::init<int, py::array_t<double>>())
        .def("find_span", &BSpline::FindSpan)
        .def("basis_funs", &BSpline::BasisFuns)
        .def("basis_funs_derivatives", &bspline_basis_funs_derivatives,
                                        py::arg("span"),
                                        py::arg("u"),
                                        py::arg("d"),
                                        R"pbdoc(
                                            Compute B-spline basis functions and their derivatives.

                                            Parameters
                                            ----------
                                            span : int
                                                Valid span index.
                                            u : float
                                                Parameter value in knot support.
                                            d : int
                                                Order of derivatives to compute.

                                            Returns
                                            -------
                                            numpy.ndarray
                                                Shape (d+1, degree+1):
                                                - row 0 = basis values
                                                - row k = kth derivative
                                        )pbdoc")
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
        .def("local_cp_ptr", static_cast<double*(Patch::*)(size_t)>(&Patch::local_cp_ptr),
             py::arg("i_local"),
             "Return pointer to local control point (as int or PyCapsule for Python?)")
        .def("control_point",
             [](const Patch& p, size_t i_local) {
                auto ptr = p.local_cp_ptr(i_local);
                ssize_t dim = p.cp_manager->dim_phys;
                return py::array_t<double>({dim}, ptr);
            },
            "Return Numpy 1D view on a local point")
        .def("local_control_point_view", &Patch::local_control_point_view,
             "Return a zero-copy NumPy array view of global control points (indexing from Python needed)")
        .def("evaluate_patch_nd", &Patch::EvaluatePatchND)
        .def("evaluate_patch_nd_omp", &Patch::EvaluatePatchNDOMP)
        .def("spans", &Patch::spans)
        .def("control_points_for_span",
            [](const Patch &self, const std::vector<int>& span) {
                auto pts = self.control_points_for_span(span);

                ssize_t dim = self.cp_manager->dim_phys;
                ssize_t n   = pts.size();

                // Build a (n, dim) numpy array
                py::array_t<double> arr({n, dim});
                double* data = arr.mutable_data();

                for (ssize_t i = 0; i < n; ++i)
                    std::memcpy(data + i*dim, pts[i], dim * sizeof(double));

                return arr;
            })
        .def("test", &Patch::Test);


    py::class_<SpanNDIterator>(m, "SpanIterator")
        .def("__iter__", [](SpanNDIterator &self) -> SpanNDIterator& {
            return self;
        })
        .def("__next__", [](SpanNDIterator &self) {
            if (self.is_done())
                throw py::stop_iteration();
            auto result = self.current();
            self.next();
            return py::cast(result);
        })
        .def("size", &SpanNDIterator::size)
        .def("current", &SpanNDIterator::current)
        .def("next", &SpanNDIterator::next);


    py::class_<ElementMatrix>(m, "ElementMatrix")
        .def_property_readonly("nb_loc", [](const ElementMatrix& self) {return self.nb_loc;})
        .def_property_readonly("global_indices", [](const ElementMatrix& self) {return self.global_indices;})
        .def("get_K_as_numpy", [](const ElementMatrix& self) {
            py::array_t<double> arr({self.nb_loc * 2, self.nb_loc * 2});
            double* data = arr.mutable_data();
            Eigen::Map<Eigen::MatrixXd> K_map(data, self.nb_loc * 2, self.nb_loc * 2);
            K_map = self.K;
            return arr;
        });
        // .def_property_readonly("K", [](const ElementMatrix& self) {return self.K;});

    py::class_<SpanGauss1D>(m, "SpanGauss1D")
        .def_property_readonly("u_param", [](const SpanGauss1D& self) { return self.u_param;})
        .def_property_readonly("weight", [](const SpanGauss1D& self) { return self.weight;})
        .def_property_readonly("N", [](const SpanGauss1D& self) {return self.N;})
        .def_property_readonly("dN", [](const SpanGauss1D& self) {return self.dN;});

    py::class_<IGABasis1D>(m, "IGABasis1D")
        .def_property_readonly("spans", [](const IGABasis1D& self) {return self.spans;})
        .def_static("build", &IGABasis1D::build, py::arg("b"), py::arg("gauss_n"));


    py::class_<IGAAssembler2D>(m, "IGAAssembler2D")
        .def(py::init<const Patch&, const IGABasis1D&, const IGABasis1D&>())
        .def("assemble_stiffness", &IGAAssembler2D::assemble_stiffness);
}
