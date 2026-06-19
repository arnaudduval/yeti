#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
#include "BSpline.hpp"
#include "BSplineTensor.hpp"
#include "ControlPointManager.hpp"
#include "DOFManager.hpp"
#include "Patch.hpp"
#include "SpanNDIterator.hpp"
#include "PatchIntegrator.hpp"
#include "PatchAssembly.hpp"
#include "refinement/RefinementOperator.hpp"
#include "refinement/HRefiner.hpp"
#include "refinement/SubdivisionRefiner.hpp"
#include "refinement/PRefiner.hpp"
#include "refinement/BezierExtractor.hpp"


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

    py::class_<ControlPointManager, std::shared_ptr<ControlPointManager>>(m, "ControlPointManager")
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

    py::class_<GlobalDOFManager>(m, "GlobalDOFManager",
        "Pool-wide map from a control-point id (in the shared "
        "ControlPointManager) to its global dof indices. Since lookup is a "
        "pure function of cp id, two patches that reference the SAME cp id "
        "(e.g. a boundary control point merged across a shared interface) "
        "automatically resolve to the SAME global dofs -- this is what "
        "makes PatchAssembly.update_dof_managers() work without any "
        "explicit tracking of which control points were merged.")
        .def(py::init<const std::vector<int>&>(), py::arg("dofs_per_control_point"))
        .def("get_dof_indices", &GlobalDOFManager::get_dof_indices, py::arg("control_point_idx"))
        .def("n_control_points", &GlobalDOFManager::n_control_points,
            "Number of control points currently covered (ids 0..n-1 are "
            "valid arguments to get_dof_indices()).")
        .def("grow", &GlobalDOFManager::grow, py::arg("new_n_cp"), py::arg("dofs_per_cp"),
            "Grow to cover new_n_cp control points (no-op if already >= that). "
            "Newly covered control points get dofs_per_cp fresh dofs each.");

    py::class_<PatchDOFManager, std::shared_ptr<PatchDOFManager>>(m, "PatchDOFManager",
        "Maps a patch's LOCAL control point position (u-fastest index into "
        "Patch.global_indices) to global dof indices. Independent of "
        "control point ids once built, so it survives PatchAssembly.compact() "
        "(which renumbers ids) unchanged.")
        .def(py::init<int, const std::vector<size_t>&, const GlobalDOFManager&>(),
             py::arg("dofs_per_control_point"), py::arg("control_points"), py::arg("global_dof_manager"))
        .def("get_global_dof_indices", &PatchDOFManager::get_global_dof_indices, py::arg("local_control_point_idx"));

    py::class_<Patch, std::shared_ptr<Patch>>(m, "Patch")
        .def(py::init<const BSplineTensor&, std::shared_ptr<ControlPointManager>, const std::vector<size_t>&, const std::vector<size_t>&>(),
             py::arg("tensor"), py::arg("cp_manager"), py::arg("global_indices"), py::arg("local_shape"))
        .def(py::init<const BSplineTensor&, std::shared_ptr<ControlPointManager>, const std::vector<size_t>&, const std::vector<size_t>&, std::shared_ptr<PatchDOFManager>>(),
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
        .def_readonly("tensor", &Patch::tensor)
        .def_property_readonly("n_cp", [](const Patch& self) { return self.global_indices.size(); },
            "Number of control points in the mapping (= product of local_shape).")
        .def_property_readonly("global_indices", [](const Patch& self) { return self.global_indices; },
            "u-fastest list mapping each local control point index to its global "
            "id in the shared ControlPointManager pool.")
        .def_property_readonly("local_shape", [](const Patch& self) { return self.local_shape; },
            "Number of basis functions per parametric direction [n_u, n_v, ...].")
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

    py::class_<IGABasis1D>(m, "IGABasis1D",
        R"doc(
Precomputed Gauss quadrature data for a 1D B-spline parametric direction.

For each non-empty knot span, stores the physical Gauss-point coordinates,
the integration weights (already multiplied by the span Jacobian), and the
values of the (p+1) active B-spline basis functions and their first
parametric derivatives at every Gauss point.

Building this object once and reusing it across assembly calls avoids
redundant evaluations of the basis functions.

Attributes
----------
gauss_spans : list[SpanGauss1D]
    One entry per non-empty knot span, in knot-vector order.  Each entry
    contains ``u_param``, ``weight``, ``N`` and ``dN`` arrays of length
    ``gauss_n``.

See Also
--------
PatchIntegrator : uses two IGABasis1D objects (one per parametric direction)
    to assemble stiffness / mass matrices over a 2-D patch.
        )doc")
        .def_property_readonly("gauss_spans", [](const IGABasis1D& self) {return self.gauss_spans;})
        .def_static("build", &IGABasis1D::build, py::arg("b"), py::arg("gauss_n"),
            R"doc(
Build an IGABasis1D from a BSpline and a Gauss-point count.

Iterates over every non-empty knot span ``[U[i], U[i+1]]``, maps the
``gauss_n`` Gauss-Legendre reference points from ``[-1, 1]`` into the span,
and evaluates both the B-spline basis functions and their first derivatives
there via de Boor's algorithm.

Parameters
----------
b : BSpline
    The 1-D B-spline whose spans are to be precomputed.
gauss_n : int
    Number of Gauss-Legendre quadrature points per span.  For an
    order-``p`` B-spline, ``gauss_n = p + 1`` integrates polynomials of
    degree ``2p`` exactly (sufficient for the stiffness matrix of a
    Laplacian-like bilinear form).

Returns
-------
IGABasis1D
    Object whose ``gauss_spans[k]`` holds the precomputed data for the
    k-th non-empty span.
            )doc");

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

    py::class_<PatchAssembly>(m, "PatchAssembly",
        R"doc(
Orchestrates multi-patch operations: detecting control points and edges
shared between patches, propagating refinement across shared edges, and
keeping DOF numbering consistent across the assembly.

Typical workflow, in order:

1. ``add_patch(patch)`` for every patch in the assembly.
2. ``detect_shared_control_points()`` -- find control point ids common to
   several patches (requires patches to already share ids in their
   ``global_indices``, e.g. by construction).
3. ``detect_interfaces()`` -- identify which (direction, side) facet of
   each patch carries a shared edge, and the orientation between them.
4. ``refine_with_propagation(patch_index, direction, refine_1d_fn)`` --
   refine one patch and automatically propagate to its neighbors across
   detected interfaces, merging the new boundary control points.
5. ``update_dof_managers(global_dof_manager, dofs_per_cp)`` -- grow the
   global DOF pool and rebuild every patch's PatchDOFManager so merged
   control points get the same global dofs on every patch.
6. ``compact()`` -- once ALL needed refinements are done, reclaim memory
   orphaned by repeated refinements of shared patches. Re-run
   ``detect_shared_control_points()`` / ``detect_interfaces()`` afterwards
   if their results are still needed (compact() invalidates both).
        )doc")
        .def(py::init<>())
        .def("add_patch", &PatchAssembly::addPatch, py::arg("patch"),
            "Add a patch to the assembly.")
        .def("detect_shared_control_points", &PatchAssembly::detectSharedControlPoints,
            "Detect control point ids shared by two or more patches (a control "
            "point is shared if the same id appears in more than one patch's "
            "global_indices). Populates the map returned by "
            "get_shared_control_points_map().")
        .def("compact", &PatchAssembly::compact,
            "Reclaim CPs orphaned by repeated refinements of shared patches. "
            "Call once, after all needed refinements are done -- rewrites the "
            "whole shared CP pool and every patch's global_indices. Clears "
            "any previously detected shared-CP map (call "
            "detect_shared_control_points() again afterwards if needed).")
        .def("detect_interfaces", &PatchAssembly::detectInterfaces,
            "Detect compatible interfaces (shared edges) between all pairs of 2D "
            "patches. For each match, records which (direction, side) facet of "
            "each patch carries the edge, the direction to refine on each side "
            "to propagate refinement (varying_direction_a/b), and whether "
            "traversal order is reversed between the two patches. Throws if "
            "two facets share control points but neither direct nor reversed "
            "order matches.")
        .def("refine_with_propagation", &PatchAssembly::refineWithPropagation,
            py::arg("patch_index"), py::arg("direction"), py::arg("refine_1d_fn"),
            "Refine patches[patch_index] along `direction`, automatically "
            "propagating the same refinement to every neighbor connected "
            "through a detected interface (see detect_interfaces()) whose "
            "varying direction matches, then merging the boundary control "
            "points created independently on both sides. "
            "refine_1d_fn(patch, direction, protected_global_ids) must refine "
            "`patch` in-place along `direction` while leaving "
            "protected_global_ids untouched. `direction` is passed explicitly "
            "because on a crossed interface (e.g. u of one patch glued to v "
            "of the other) the neighbor must be refined along its OWN "
            "matching direction, not necessarily the one given to "
            "refine_with_propagation(). Example: "
            "lambda patch, d, ids: SubdivisionRefiner(d, n_levels).refine_1d(patch, ids). "
            "Requires detect_interfaces() (and detect_shared_control_points()) "
            "to have been called first. Does not update any PatchDOFManager "
            "(see update_dof_managers()).")
        .def("update_dof_managers", &PatchAssembly::updateDOFManagers,
            py::arg("global_dof_manager"), py::arg("dofs_per_cp"),
            "Grow `global_dof_manager` to cover every control point currently "
            "referenced by the assembly's patches, then rebuild each patch's "
            "PatchDOFManager from its CURRENT global_indices. Because "
            "GlobalDOFManager maps dofs by control-point id, two patches "
            "referencing the same (merged) boundary CP automatically get the "
            "same global dofs. `dofs_per_cp` only applies to control points "
            "newly covered by the growth (pre-existing CPs keep their dofs; "
            "each patch keeps its own already-set dofs_per_control_point). "
            "Call this AFTER refine_with_propagation() and BEFORE compact() "
            "-- compact() renumbers control point ids, which would "
            "desynchronize this cp-id-indexed mapping if dofs were assigned "
            "first against soon-to-be-stale ids.")
        .def("get_interfaces", [](const PatchAssembly& self) {
                py::list result;
                for (const auto& itf : self.getInterfaces()) {
                    py::dict d;
                    d["patch_a"] = itf.patch_a;
                    d["direction_a"] = itf.direction_a;
                    d["side_a"] = itf.side_a;
                    d["varying_direction_a"] = itf.varying_direction_a;
                    d["patch_b"] = itf.patch_b;
                    d["direction_b"] = itf.direction_b;
                    d["side_b"] = itf.side_b;
                    d["varying_direction_b"] = itf.varying_direction_b;
                    d["reversed"] = itf.reversed;
                    result.append(d);
                }
                return result;
            })
        .def("get_patchs_sharing_control_points", &PatchAssembly::getPatchsSharingControlPoints, py::arg("global_cp_index"),
            "Return the list of patch indices sharing the given global "
            "control point id (empty if detect_shared_control_points() "
            "hasn't found it, or hasn't been called).")
        .def("get_shared_control_points_between", [](const PatchAssembly& self, const Patch& patch1, const Patch& patch2) {
                return self.getSharedControlPoints(patch1, patch2);
            }, py::arg("patch1"), py::arg("patch2"),
            "Return the global control point ids common to exactly these two "
            "patches' global_indices. Computed directly from the two patches "
            "-- does not require detect_shared_control_points() to have "
            "been called first.")
        .def("get_control_points_for_patch", &PatchAssembly::getControlPointsForPatch, py::arg("patch_index"),
            "Return the number of control points of patches[patch_index].")
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
        .def("get_patchs", &PatchAssembly::getPatchs,
            "Return the list of patches added to this assembly, in add_patch() order.")
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
            },
            "Return the map {global_cp_id: [patch indices sharing it]} built "
            "by detect_shared_control_points() (empty if not yet called).")
        .def("get_transformation_matrices", &PatchAssembly::getTransformationMatrices,
            "Return the per-patch refinement transformation matrices "
            "(identity until set_transformation_matrix() is called).");

    py::class_<RefinementOperator, std::shared_ptr<RefinementOperator>>(m, "RefinementOperator")
        .def("get_type", &RefinementOperator::getType)
        .def("refine", [](const RefinementOperator& self, Patch& patch) {
            Eigen::MatrixXd T;
            self.refine(patch, T);
            return T;
        }, py::arg("patch"),
           "Refine patch in-place. Returns transition_matrix (nb_new_cp x nb_old_cp).");

    py::class_<HRefiner, RefinementOperator, std::shared_ptr<HRefiner>>(m, "HRefiner")
        .def(py::init<int, double>(), py::arg("direction"), py::arg("knot"))
        .def("refine", [](const HRefiner& self, Patch& patch) {
            Eigen::MatrixXd T;
            self.refine(patch, T);
            return T;
        }, py::arg("patch"),
           "Insert knot in-place (h-refinement). Returns transition_matrix (nb_new_cp x nb_old_cp).");

    py::class_<SubdivisionRefiner, RefinementOperator, std::shared_ptr<SubdivisionRefiner>>(m, "SubdivisionRefiner")
        .def(py::init<int, int>(), py::arg("direction"), py::arg("n_levels") = 1)
        .def("refine", [](const SubdivisionRefiner& self, Patch& patch) {
            Eigen::MatrixXd T;
            self.refine(patch, T);
            return T;
        }, py::arg("patch"),
           "Bisect all knot spans in-place (subdivision). Returns composed transition_matrix (nb_final_cp x nb_initial_cp).")
        .def("refine_1d", [](const SubdivisionRefiner& self, Patch& patch,
                              const std::unordered_set<size_t>& protected_global_ids) {
            Eigen::MatrixXd T;
            self.refine_1d(patch, T, protected_global_ids);
            return T;
        }, py::arg("patch"), py::arg("protected_global_ids") = std::unordered_set<size_t>{},
           "Fast path: refine in-place, return only the 1D transition matrix (n_new_1d x n_old_1d). "
           "Use nd_transition_from_1d() to build the full nD matrix if needed. "
           "protected_global_ids: ids borrowed from another patch (shared interface) that must "
           "never be recomputed/renumbered; if empty (default), this patch is assumed to "
           "exclusively own its cp_manager range.");

    py::class_<PRefiner, RefinementOperator, std::shared_ptr<PRefiner>>(m, "PRefiner")
        .def(py::init<int, int>(), py::arg("direction"), py::arg("n_elevations") = 1)
        .def("refine", [](const PRefiner& self, Patch& patch) {
            Eigen::MatrixXd T;
            self.refine(patch, T);
            return T;
        }, py::arg("patch"),
           "Elevate degree by 1 in-place (P&T A5.9). Returns transition_matrix (nb_new_cp x nb_old_cp).")
        .def("refine_1d", [](const PRefiner& self, Patch& patch,
                              const std::unordered_set<size_t>& protected_global_ids) {
            Eigen::MatrixXd T;
            self.refine_1d(patch, T, protected_global_ids);
            return T;
        }, py::arg("patch"), py::arg("protected_global_ids") = std::unordered_set<size_t>{},
           "Fast path: elevate in-place, return only the 1D transition matrix (n_new_1d x n_old_1d). "
           "Use nd_transition_from_1d() to build the full nD matrix if needed. "
           "protected_global_ids: ids borrowed from another patch (shared interface) that must "
           "never be recomputed/renumbered; if empty (default), this patch is assumed to "
           "exclusively own its cp_manager range.");

    m.def("nd_transition_from_1d",
        [](const Eigen::MatrixXd& T_1d, int direction,
           const std::vector<size_t>& shape_before,
           const std::vector<size_t>& shape_after) {
            return nd_transition_from_1d(T_1d, direction, shape_before, shape_after);
        },
        py::arg("T_1d"), py::arg("direction"),
        py::arg("shape_before"), py::arg("shape_after"),
        "Build the full nD transition matrix from a 1D one via Kronecker products.\n"
        "shape_before/after are the local_shape of the patch before/after refinement.");

    py::class_<BezierElementND>(m, "BezierElementND",
        R"pbdoc(
        Result of BezierExtractor.extract_nd() for one tensor-product element.

        Attributes
        ----------
        C : ndarray, shape (n_local, n_local)
            Extraction matrix: N_active(xi) = C @ B_nd(xi_hat)
            where n_local = prod_d(p_d + 1).
        active_indices : list[int]
            u-fastest flat CP indices of the active B-spline functions.
            Use patch.control_point(i) for i in active_indices.
        elem_index : list[int]
            0-based element multi-index (e_0, e_1, ...).
        )pbdoc")
        .def_readonly("C", &BezierElementND::C)
        .def_property_readonly("active_indices",
            [](const BezierElementND& self) { return self.active; })
        .def_property_readonly("elem_index",
            [](const BezierElementND& self) { return self.elem_index; });

    py::class_<BezierExtractor>(m, "BezierExtractor",
        R"pbdoc(
        Bézier extraction operator (Borden et al. 2011).

        For a B-spline of degree p with ne elements, returns ne matrices C^e of
        shape (p+1, p+1) such that:

            N_active_e(xi) = C^e @ B^p(xi_hat)

        where N_active_e are the p+1 active B-spline basis functions on element e,
        B^p are the Bernstein polynomials on [0,1], and xi_hat is the local
        coordinate mapped to [0,1].

        For a 2D patch, the element operator is the Kronecker product of 1D operators:
            C_element = np.kron(Ce_v[ev], Ce_u[eu])

        Parameters
        ----------
        direction : int
            Parametric direction (0 = u, 1 = v, 2 = w).

        Example
        -------
        Ce_u = BezierExtractor(direction=0).extract(patch)  # list of (p_u+1, p_u+1) matrices
        Ce_v = BezierExtractor(direction=1).extract(patch)  # list of (p_v+1, p_v+1) matrices
        spans_u = BezierExtractor.spans(patch, direction=0) # span index per element
        )pbdoc")
        .def(py::init<int>(), py::arg("direction"))
        .def("extract", [](const BezierExtractor& self, const Patch& patch) -> py::list {
            auto Ce = self.extract(patch);
            py::list result;
            for (const auto& mat : Ce)
                result.append(mat);
            return result;
        }, py::arg("patch"),
           "Return list of (p+1, p+1) extraction matrices, one per element.")
        .def_static("extract_1d", [](const BSpline& spline) -> py::list {
            auto Ce = BezierExtractor::extract_1d(spline);
            py::list result;
            for (const auto& mat : Ce)
                result.append(mat);
            return result;
        }, py::arg("spline"),
           "1D extraction: list of (p+1, p+1) matrices for a single BSpline.")
        .def_static("element_spans", [](const BSpline& spline) {
            return BezierExtractor::element_spans(spline);
        }, py::arg("spline"),
           "Return span index k (in the knot vector) for each element.")
        .def_static("extract_nd", &BezierExtractor::extract_nd,
           py::arg("patch"),
           R"pbdoc(
           ND Bézier extraction: one BezierElementND per tensor-product element.

           Elements are returned in u-fastest order (direction 0 varies fastest).
           The extraction matrix C satisfies N_active(xi) = C @ B_nd(xi_hat)
           where B_nd is the tensor-product Bernstein basis (u-fastest).

           Parameters
           ----------
           patch : Patch

           Returns
           -------
           list[BezierElementND]
           )pbdoc");
}
