#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>
#include "BSpline.hpp"
#include "BSplineTensor.hpp"
#include "ControlPointManager.hpp"
#include "DOFManager.hpp"
#include "Patch.hpp"
#include "SpanNDIterator.hpp"
#include "PatchIntegrator.hpp"
#include "PatchAssembly.hpp"


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

    py::class_<GlobalDOFManager>(m, "GlobalDOFManager")
        .def(py::init<const std::vector<int>&>(), py::arg("dofs_per_control_point"))
        .def("get_dof_indices", &GlobalDOFManager::get_dof_indices, py::arg("control_point_idx"));

    py::class_<PatchDOFManager, std::shared_ptr<PatchDOFManager>>(m, "PatchDOFManager")
        .def(py::init<int, const std::vector<size_t>&, const GlobalDOFManager&>(),
             py::arg("dofs_per_control_point"), py::arg("control_points"), py::arg("global_dof_manager"))
        .def("get_global_dof_indices", &PatchDOFManager::get_global_dof_indices, py::arg("local_control_point_idx"));

    py::class_<Patch, std::shared_ptr<Patch>>(m, "Patch")
        .def(py::init<const BSplineTensor&, ControlPointManager*, const std::vector<size_t>&, const std::vector<size_t>&>(),
             py::arg("tensor"), py::arg("cp_manager"), py::arg("global_indices"), py::arg("local_shape"))
        .def(py::init<const BSplineTensor&, ControlPointManager*, const std::vector<size_t>&, const std::vector<size_t>&, std::shared_ptr<PatchDOFManager>>(),
             py::arg("tensor"), py::arg("cp_manager"), py::arg("global_indices"), py::arg("local_shape"), py::arg("dof_manager"))
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
            },
            py::arg("span"),
            "Return array containing control points coordinates for a given span. Warning : data is return as stored in memory and need to be reshaped/transposed for proper use"
        )
        .def("test", &Patch::Test)
        .def_property_readonly("dof_manager", [](const Patch& self) -> std::shared_ptr<PatchDOFManager> { return self.dof_manager; });


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

    py::class_<SpanGauss1D>(m, "SpanGauss1D")
        .def_property_readonly("u_param", [](const SpanGauss1D& self) { return self.u_param;})
        .def_property_readonly("weight", [](const SpanGauss1D& self) { return self.weight;})
        .def_property_readonly("N", [](const SpanGauss1D& self) {
            std::vector<py::array_t<double>> arrays;
            for (const auto& vec : self.N) {
                py::array_t<double> arr(vec.size());
                Eigen::Map<Eigen::VectorXd>(arr.mutable_data(), vec.size()) = vec;
                arrays.push_back(arr);
            }
            return arrays;
        })
        .def_property_readonly("dN", [](const SpanGauss1D& self) {
            std::vector<py::array_t<double>> arrays;
            for (const auto& vec : self.dN) {
                py::array_t<double> arr(vec.size());
                Eigen::Map<Eigen::VectorXd>(arr.mutable_data(), vec.size()) = vec;
                arrays.push_back(arr);
            }
            return arrays;
        });
        // .def_property_readonly("N", [](const SpanGauss1D& self) {return self.N;})
        // .def_property_readonly("dN", [](const SpanGauss1D& self) {return self.dN;});

    py::class_<IGABasis1D>(m, "IGABasis1D")
        .def_property_readonly("gauss_spans", [](const IGABasis1D& self) {return self.gauss_spans;})
        .def_static("build", &IGABasis1D::build, py::arg("b"), py::arg("gauss_n"));

    py::class_<MaterialProperties>(m, "MaterialProperties")
        .def(py::init<double, double, double>(),
             py::arg("E"), py::arg("nu"), py::arg("thickness") = 1.0)
        .def_readwrite("E", &MaterialProperties::E)
        .def_readwrite("nu", &MaterialProperties::nu)
        .def_readwrite("thickness", &MaterialProperties::thickness);

    py::class_<PatchIntegrator>(m, "PatchIntegrator")
        .def(py::init<const Patch&, const IGABasis1D&, const IGABasis1D&, const MaterialProperties&>(),
             py::arg("patch"), py::arg("basis_u"), py::arg("basis_v"), py::arg("material_properties"))
        .def("integrate", &PatchIntegrator::integrate);

    py::class_<PatchAssembly>(m, "PatchAssembly")
        .def(py::init<>())
        .def("add_patch", &PatchAssembly::addPatch, py::arg("patch"))
        .def("detect_shared_control_points", &PatchAssembly::detectSharedControlPoints)
        .def("get_patchs_sharing_control_points", &PatchAssembly::getPatchsSharingControlPoints, py::arg("global_cp_index"))
        .def("get_shared_control_poins", [](const PatchAssembly& self, const Patch& patch1, const Patch& patch2) {
                return self.getSharedControlPoints(patch1, patch2);
            }, py::arg("patch1"), py::arg("patch2"))
        .def("get_control_points_for_patch", &PatchAssembly::getControlPointsForPatch, py::arg("patch_index"))
        .def("apply_transformation_to_control_points", [](const PatchAssembly& self, const std::vector<Eigen::VectorXd>& old_control_points, size_t patch_index) {
                return self.applyTransformationToControlPoints(old_control_points, patch_index);
            }, py::arg("old_control_points"), py::arg("patch_index")
        )
        .def("apply_transformation_to_dofs", [](const PatchAssembly& self, const Eigen::VectorXd& old_dofs, size_t patch_index, size_t dim_phys) {
                return self.applyTransformationToDOFs(old_dofs, patch_index, dim_phys);
            }, py::arg("old_dofs"), py::arg("patch_index"), py::arg("dim_phys")
        )
        .def("set_transformation_matrix", [](PatchAssembly& self, size_t patch_index, const Eigen::MatrixXd matrix) {
                    return self.setTransformationMatrix(patch_index, matrix);
            }, py::arg("patch_index"), py::arg("matrix")
        )
        .def("get_patchs", &PatchAssembly::getPatchs)
        .def("get_shared_control_points_map", [](const PatchAssembly& self) {
                const auto& shared_map = self.getSharedControlPoints();
                py::dict result;
                for (const auto& [cp_idx, patch_indices] : shared_map) {
                    py::list patch_list;
                    for (size_t idx : patch_indices) {
                        patch_list.append(py::cast(idx));
                    }
                    result[py::cast(cp_idx)] = patch_list;
                }
                return result;
            })
        .def("get_transformation_matrices", &PatchAssembly::getTransformationMatrices);
}
