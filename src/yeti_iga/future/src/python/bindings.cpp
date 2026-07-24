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
#include "LocalOperator.hpp"
#include "ScalarLocalOperator.hpp"
#include "Traction.hpp"
#include "refinement/RefinementOperator.hpp"
#include "refinement/HRefiner.hpp"
#include "refinement/SubdivisionRefiner.hpp"
#include "refinement/PRefiner.hpp"
#include "refinement/BezierExtractor.hpp"
#include "PatchEvaluator.hpp"


namespace py = pybind11;

// Trampoline letting LocalOperator be subclassed from Python (development/
// testing convenience -- see LocalOperator.hpp).
class PyLocalOperator : public LocalOperator {
public:
    using LocalOperator::LocalOperator;
    Eigen::MatrixXd computeIntegrand(
        const std::vector<double>& R, const std::vector<double>& dRdx,
        const std::vector<double>& dRdy) const override {
        PYBIND11_OVERRIDE_PURE_NAME(
            Eigen::MatrixXd, LocalOperator, "compute_integrand", computeIntegrand,
            R, dRdx, dRdy);
    }
};

// Trampoline letting ScalarLocalOperator be subclassed from Python (scalar-
// valued integration terms such as error norms -- see ScalarLocalOperator.hpp).
class PyScalarLocalOperator : public ScalarLocalOperator {
public:
    using ScalarLocalOperator::ScalarLocalOperator;
    double computeScalarIntegrand(
        const std::vector<double>& R, const std::vector<double>& dRdx,
        const std::vector<double>& dRdy,
        const Eigen::Vector2d& physical_point,
        const Eigen::VectorXd& u_local) const override {
        PYBIND11_OVERRIDE_PURE_NAME(
            double, ScalarLocalOperator, "compute_scalar_integrand",
            computeScalarIntegrand,
            R, dRdx, dRdy, physical_point, u_local);
    }
};

// Trampoline letting Traction be subclassed from Python (position-dependent
// boundary loads -- see Traction.hpp).
class PyTraction : public Traction {
public:
    using Traction::Traction;
    Eigen::Vector2d evaluate(const Eigen::Vector2d& physical_point) const override {
        PYBIND11_OVERRIDE_PURE_NAME(
            Eigen::Vector2d, Traction, "evaluate", evaluate,
            physical_point);
    }
};

// Trampoline letting ConstitutiveLaw be subclassed from Python.
// Explicit public constructor calls the protected base constructor.
class PyConstitutiveLaw : public ConstitutiveLaw {
public:
    explicit PyConstitutiveLaw(const Material& m) : ConstitutiveLaw(m) {}
    int n_dofs_per_cp() const override {
        PYBIND11_OVERRIDE_PURE_NAME(
            int, ConstitutiveLaw, "n_dofs_per_cp", n_dofs_per_cp);
    }
    Eigen::MatrixXd stiffness_density(
        const Eigen::VectorXd& grad_a,
        const Eigen::VectorXd& grad_b,
        const Eigen::VectorXd& x_phys) const override
    {
        PYBIND11_OVERRIDE_PURE_NAME(
            Eigen::MatrixXd, ConstitutiveLaw, "stiffness_density", stiffness_density,
            grad_a, grad_b, x_phys);
    }
};

PYBIND11_MODULE(bspline, m)
{
    py::class_<BSpline>(m, "BSpline",
        "1D B-spline basis of degree p over a knot vector. "
        "find_span(u) locates the knot span enclosing the parameter value u; "
        "basis_funs(span, u) evaluates the p+1 non-zero basis functions there; "
        "basis_funs_derivatives(span, u, d) adds derivatives up to order d.")
        .def(py::init<int, py::array_t<double>>(),
             py::arg("degree"), py::arg("knot_vector"),
             "Construct a B-spline of the given polynomial degree over knot_vector. "
             "knot_vector must be non-decreasing with length >= 2*(degree+1).")
        .def("find_span", &BSpline::FindSpan,
             py::arg("u"),
             "Return the index i of the knot span [U[i], U[i+1]) containing u.\n\n"
             "Uses binary search (Cox-de Boor convention). u must lie in [U[p], U[n+1]].\n"
             "The result is the span to pass to basis_funs(), basis_funs_derivatives(), "
             "and one_basis_fun().")
        .def("basis_funs", &BSpline::BasisFuns,
             py::arg("span"), py::arg("u"),
             "Evaluate the p+1 non-zero B-spline basis functions N_{span-p,p}..N_{span,p} at u.\n\n"
             "span must be a valid knot-span index returned by find_span(u). "
             "Returns a 1D NumPy array of length degree+1. "
             "For all non-zero functions AND their derivatives use basis_funs_derivatives().")
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
        .def("one_basis_fun", &BSpline::OneBasisFun,
             py::arg("u"), py::arg("i"),
             "Evaluate the single B-spline basis function N_{i,p}(u). "
             "Useful when only one function value is needed; for all p+1 "
             "non-zero functions at a span use basis_funs().")
        .def_property_readonly("degree", &BSpline::getDegree,
            "Polynomial degree p.")
        .def_property_readonly("knot_vector", &BSpline::kvView,
            "Knot vector as a 1D NumPy read-only view.");

    py::class_<BSplineTensor>(m, "BSplineTensor",
        "Tensor-product extension of several 1D B-splines. Evaluates the "
        "multivariate basis over a Cartesian product of parametric directions. "
        "BSplineSurface (2D) and BSplineVolume (3D) are the concrete subclasses.")
        .def_property_readonly("components", [](const BSplineTensor& t) {return t.components; },
            "list[BSpline] — one 1D B-spline per parametric direction (u, v[, w]).")
        .def("basis_funs_nd", &BSplineTensor::BasisFunsND,
             py::arg("span"), py::arg("u"),
             "Evaluate all non-zero tensor-product basis functions at a parametric point.\n"
             "span : int array, shape (n_param_dims,) -- knot-span indices (from find_span_nd).\n"
             "u    : float array, shape (n_param_dims,) -- parametric coordinates.\n"
             "Returns: 1D float array of length prod_d(degree_d + 1), u-fastest order.")
        .def("find_span_nd",
             static_cast<py::array_t<int>(BSplineTensor::*)(const py::array_t<double>&) const>
                (&BSplineTensor::FindSpanND),
             py::arg("u"),
             "Compute span indices in all parameter dimensions.");

    py::class_<BSplineSurface, BSplineTensor>(m, "BSplineSurface",
        "Tensor-product B-spline surface basis built from two 1D BSplines (u and v). "
        "Pass to Patch as the tensor to define the parametric-to-physical mapping of a 2D patch.")
        .def(py::init<const BSpline&, const BSpline&>(),
             py::arg("su"), py::arg("sv"),
             "Construct from a BSpline in the u direction and one in the v direction.");

    py::class_<BSplineVolume, BSplineTensor>(m, "BSplineVolume",
        "Tensor-product B-spline volume basis built from three 1D BSplines (u, v, w). "
        "Pass to Patch as the tensor to define the parametric-to-physical mapping of a 3D patch.")
        .def(py::init<const BSpline&, const BSpline&, const BSpline&>(),
             py::arg("su"), py::arg("sv"), py::arg("sw"),
             "Construct from a BSpline in each of the u, v, and w directions.");

    py::class_<ControlPointManager, std::shared_ptr<ControlPointManager>>(m, "ControlPointManager",
        "Shared pool of physical control points (and optional NURBS weights) "
        "referenced by one or more patches. Add points via add_point(); patches "
        "then reference them by global index. Setting any weight != 1.0 activates "
        "NURBS rational mode (is_rational becomes True); otherwise the pure "
        "B-spline code path runs with no overhead.")
        .def(py::init<int>(), py::arg("dim")=3,
             "Construct an empty pool for control points in dim-dimensional "
             "physical space (default 3).")
        .def_property_readonly("dim_phys", [](const ControlPointManager& mgr) {return mgr.dim_phys; },
            "Physical space dimension (2 for 2D, 3 for 3D).")
        .def("add_point", &ControlPointManager::add_point,
             py::arg("coords"), py::arg("w") = 1.0,
             "Add a control point with optional NURBS weight (default 1.0 = pure "
             "B-spline, no overhead). When w != 1.0 or is_rational is already true, "
             "activates rational mode and stores w alongside the point.")
        .def("set_weight", &ControlPointManager::set_weight,
             py::arg("id"), py::arg("w"),
             "Set the NURBS weight of an existing control point. Activates rational "
             "mode (all previous points get weight 1.0 if not yet set). Reverts to "
             "B-spline mode automatically if all weights are 1.0 afterwards.")
        .def_property_readonly("is_rational", &ControlPointManager::is_rational,
             "True if any control point has a weight != 1.0 (NURBS mode). When false "
             "all integration and evaluation uses the B-spline fast path with no "
             "extra overhead.")
        .def_property_readonly("n_points", &ControlPointManager::n_points,
            "Number of control points currently in the pool.")
        .def("coords_view", [](ControlPointManager& self){
            auto capsule = py::capsule(&self);
            std::vector<ssize_t> shape = {(ssize_t)self.n_points(), (ssize_t)self.dim_phys};
            std::vector<ssize_t> strides = {
                static_cast<std::ptrdiff_t>(self.dim_phys * sizeof(double)),
                static_cast<std::ptrdiff_t>(sizeof(double))
            };
            return py::array_t<double>(shape, strides, self.coords.data(), capsule);
        }, "Zero-copy NumPy view of all control point coordinates, shape (n_points, dim_phys). "
           "Row i is the physical coordinates of the control point with global id i.")
        .def("weights_view", [](ControlPointManager& self) -> py::array_t<double> {
            if (!self.is_rational())
                throw std::runtime_error(
                    "ControlPointManager.weights_view(): manager is not rational "
                    "(all weights are 1.0). Call add_point(..., w) or set_weight() first.");
            auto capsule = py::capsule(&self);
            return py::array_t<double>(
                {(ssize_t)self.n_points()},
                {(ssize_t)sizeof(double)},
                self.weights.data(), capsule);
        }, "Zero-copy view of NURBS weights as a 1D numpy array (shape n_points). "
           "Only available when is_rational is True.");

    py::class_<GlobalDOFManager>(m, "GlobalDOFManager",
        "Pool-wide map from a control-point id (in the shared "
        "ControlPointManager) to its global dof indices. Since lookup is a "
        "pure function of cp id, two patches that reference the SAME cp id "
        "(e.g. a boundary control point merged across a shared interface) "
        "automatically resolve to the SAME global dofs -- this is what "
        "makes PatchAssembly.update_dof_managers() work without any "
        "explicit tracking of which control points were merged.")
        .def(py::init<const std::vector<int>&>(), py::arg("dofs_per_control_point"),
            "Construct with a list of dof counts (one per existing control point). "
            "Prefer GlobalDOFManager([dofs_per_cp] * n_cp) for uniform problems.")
        .def("get_dof_indices", &GlobalDOFManager::get_dof_indices, py::arg("control_point_idx"),
            "Return the global DOF indices for control point control_point_idx. "
            "The returned list has length dofs_per_cp (typically 2 for 2-D, 3 for 3-D). "
            "Raises if control_point_idx >= n_control_points().")
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
             py::arg("dofs_per_control_point"), py::arg("control_points"), py::arg("global_dof_manager"),
             "Build from the patch's dofs_per_cp, its global_indices list, and the "
             "assembly-wide GlobalDOFManager. Typically called by "
             "PatchAssembly.update_dof_managers() rather than directly.")
        .def("get_global_dof_indices", &PatchDOFManager::get_global_dof_indices, py::arg("local_control_point_idx"),
            "Return the assembly-wide DOF indices for local control point position "
            "local_control_point_idx (u-fastest index into Patch.global_indices). "
            "Combine with Patch.boundary_control_points() to identify constrained "
            "DOFs for Dirichlet boundary conditions.");

    py::class_<Patch, std::shared_ptr<Patch>>(m, "Patch",
        "Links a BSplineTensor (parametric basis), a shared ControlPointManager "
        "(physical control points and optional NURBS weights), and an optional "
        "PatchDOFManager (DOF mapping). One patch of a multi-patch IGA model. "
        "The DOF manager is required for integration and solution evaluation; "
        "attach it at construction or via PatchAssembly.update_dof_managers().")
        .def(py::init<const BSplineTensor&, std::shared_ptr<ControlPointManager>, const std::vector<size_t>&, const std::vector<size_t>&>(),
             py::arg("tensor"), py::arg("cp_manager"), py::arg("global_indices"), py::arg("local_shape"))
        .def(py::init<const BSplineTensor&, std::shared_ptr<ControlPointManager>, const std::vector<size_t>&, const std::vector<size_t>&, std::shared_ptr<PatchDOFManager>>(),
             py::arg("tensor"), py::arg("cp_manager"), py::arg("global_indices"), py::arg("local_shape"), py::arg("dof_manager"))
        .def("control_point",
             [](const Patch& p, size_t i_local) {
                auto ptr = p.local_cp_ptr(i_local);
                ssize_t dim = p.cp_manager->dim_phys;
                return py::array_t<double>({dim}, ptr);
            },
            py::arg("i_local"),
            "Zero-copy 1D NumPy view (length dim_phys) of the coordinates of "
            "local control point i_local (u-fastest index).")
        .def("local_control_point_view", &Patch::local_control_point_view,
             "Zero-copy NumPy view of all local control points. "
             "Shape: local_shape + (dim_phys,), row-major. "
             "Last axis is the coordinate (x, y[, z]).")
        .def("evaluate_patch_nd", &Patch::EvaluatePatchND,
             py::arg("spans"), py::arg("params"),
             "Evaluate the geometric mapping (parametric -> physical) at n_pts points.\n"
             "spans  : int array, shape (n_pts, n_param_dims) -- knot-span index per direction.\n"
             "params : float array, shape (n_pts, n_param_dims) -- parametric coordinates.\n"
             "Returns: float array, shape (n_pts, dim_phys) -- physical coordinates.\n"
             "Single-threaded. For the parallel variant use evaluate_patch_nd_omp().\n"
             "Note: evaluates the geometry, not an FE solution field. "
             "For solution evaluation use PatchEvaluator.evaluate_solution().")
        .def("evaluate_patch_nd_omp", &Patch::EvaluatePatchNDOMP,
             py::arg("spans"), py::arg("params"),
             "OpenMP-parallel version of evaluate_patch_nd(). Same signature and "
             "return value; use for large point sets where parallelism pays off.")
        .def("spans", &Patch::spans,
             "Iterate over all non-degenerate knot spans of this patch. "
             "Returns a SpanIterator whose items are lists [i_u, i_v, ...] "
             "(one span index per parametric direction, u-fastest order). "
             "Zero-length spans from knot repetitions are skipped automatically.")
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
            "Return the coordinates of the active control points for a given span "
            "(u-fastest order, same order as PatchIntegrator uses internally). "
            "Shape: (n_active, dim_phys) where n_active = product of (degree+1) "
            "per direction. Low-level; used internally by PatchIntegrator.")
        .def("boundary_control_points", &Patch::boundary_control_points,
            py::arg("direction"), py::arg("side"), py::arg("span_min") = -1, py::arg("span_max") = -1,
            "Patch-LOCAL flat positions (u-fastest) of the control points on "
            "the edge obtained by fixing `direction` at its first (side=0) "
            "or last (side=1) local index. If span_min/span_max are given "
            "(both >= 0, raw knot-span indices like the ones spans() "
            "yields), restricts the selection to control points active "
            "over those spans of the OTHER direction; -1/-1 (default) "
            "selects the whole edge. Building block for Dirichlet boundary "
            "conditions at edge/span-range granularity -- combine with "
            "dof_manager.get_global_dof_indices(local_pos) to get global "
            "dofs (control-point granularity needs no helper: any local "
            "position works directly with the DOF manager).")
        .def_readonly("tensor", &Patch::tensor,
            "The BSplineTensor defining the parametric-to-physical mapping of this patch.")
        .def_property_readonly("n_cp", [](const Patch& self) { return self.global_indices.size(); },
            "Number of control points in the mapping (= product of local_shape).")
        .def_property_readonly("global_indices", [](const Patch& self) { return self.global_indices; },
            "u-fastest list mapping each local control point index to its global "
            "id in the shared ControlPointManager pool.")
        .def_property_readonly("local_shape", [](const Patch& self) { return self.local_shape; },
            "Number of basis functions per parametric direction [n_u, n_v, ...].")
        .def_property_readonly("dof_manager", [](const Patch& self) -> std::shared_ptr<PatchDOFManager> { return self.dof_manager; },
            "PatchDOFManager attached to this patch, or None if the patch was "
            "constructed without one (geometry-only phase, before DOF assignment). "
            "Set by PatchAssembly.update_dof_managers() or by passing dof_manager "
            "to the constructor. Required by PatchIntegrator and PatchEvaluator.")
        .def_property_readonly("cp_manager", [](const Patch& self) -> std::shared_ptr<ControlPointManager> { return self.cp_manager; },
            "Shared ControlPointManager holding the coordinates (and optional NURBS weights) "
            "of every control point referenced by this patch.");


    py::class_<SpanNDIterator>(m, "SpanIterator",
        "Iterator over the non-degenerate knot spans of a tensor-product B-spline. "
        "Each item is a list of knot-span indices [i_u, i_v, ...] (one per "
        "parametric direction, u-fastest). Zero-length spans produced by knot "
        "repetitions (e.g. a C0 interface knot) are skipped automatically. "
        "Obtain via Patch.spans().")
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
        .def("size", &SpanNDIterator::size,
             "Total number of non-degenerate spans in this iterator.")
        .def("current", &SpanNDIterator::current,
             "Current span as a list of knot-span indices [i_u, i_v, ...].")
        .def("next", &SpanNDIterator::next,
             "Advance to the next non-degenerate span.");

    py::class_<SpanGauss1D>(m, "SpanGauss1D",
        "Precomputed Gauss quadrature data for one (span, Gauss-point) pair in "
        "one parametric direction. Stored inside IGABasis1D.gauss_spans. "
        "Fields: u_param (parametric coordinate), weight (Gauss weight), "
        "N (list of basis function values, one entry per active function), "
        "dN (corresponding derivatives).")
        .def_property_readonly("u_param", [](const SpanGauss1D& self) { return self.u_param;},
            "Parametric coordinate of this Gauss point in the knot-vector space.")
        .def_property_readonly("weight", [](const SpanGauss1D& self) { return self.weight;},
            "Gauss quadrature weight (already includes the span half-length Jacobian factor).")
        .def_property_readonly("N", [](const SpanGauss1D& self) {
            std::vector<py::array_t<double>> arrays;
            for (const auto& vec : self.N) {
                py::array_t<double> arr(vec.size());
                Eigen::Map<Eigen::VectorXd>(arr.mutable_data(), vec.size()) = vec;
                arrays.push_back(arr);
            }
            return arrays;
        }, "list of arrays, one per Gauss point: basis function values for the p+1 "
           "active basis functions in this span (u-fastest order).")
        .def_property_readonly("dN", [](const SpanGauss1D& self) {
            std::vector<py::array_t<double>> arrays;
            for (const auto& vec : self.dN) {
                py::array_t<double> arr(vec.size());
                Eigen::Map<Eigen::VectorXd>(arr.mutable_data(), vec.size()) = vec;
                arrays.push_back(arr);
            }
            return arrays;
        }, "list of arrays, one per Gauss point: first parametric derivatives of the "
           "p+1 active basis functions (same order as N).");
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
        .def_property_readonly("gauss_spans", [](const IGABasis1D& self) {return self.gauss_spans;},
            "list[SpanGauss1D] — one entry per non-empty knot span, in knot-vector order. "
            "Each SpanGauss1D holds the precomputed u_param, weight, N and dN arrays "
            "for gauss_n quadrature points in that span.")
        .def_property_readonly("span_indices",
            [](const IGABasis1D& self) {
                // Return as a dict {knot_span_index: gauss_span_index} (sorted for determinism)
                return std::map<int,int>(self.span_indices.begin(), self.span_indices.end());
            },
            "Dict mapping each non-empty knot span index to its position in gauss_spans. "
            "Mirrors the internal span_indices used by PatchIntegrator.")
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

    py::class_<Material>(m, "Material",
        "Physical material constants (E, nu, rho, thickness) with derived elastic "
        "moduli (mu, lambda_3d, bulk_modulus). Pass to a ConstitutiveLaw to set "
        "the mechanical behaviour for integration.")
        .def(py::init([](double E, double nu, double rho, double thickness) {
                return Material{E, nu, rho, thickness};
            }),
            py::arg("E"), py::arg("nu"), py::arg("rho") = 0.0, py::arg("thickness") = 1.0,
            "Construct a material. E: Young's modulus; nu: Poisson's ratio; "
            "rho: mass density (>0 for mass integration, default 0); "
            "thickness: out-of-plane thickness for 2D plane problems (default 1).")
        .def_readwrite("E", &Material::E, "Young's modulus.")
        .def_readwrite("nu", &Material::nu, "Poisson's ratio.")
        .def_readwrite("rho", &Material::rho,
            "Mass density. Must be > 0 to use integrate_mass()/assemble_mass().")
        .def_readwrite("thickness", &Material::thickness,
            "Out-of-plane thickness for 2D plane problems (default 1.0).")
        .def("mu", &Material::mu, "Shear modulus E / (2*(1+nu)).")
        .def("lambda_3d", &Material::lambda_3d,
            "3D first Lame parameter nu*E / ((1+nu)*(1-2*nu)).")
        .def("bulk_modulus", &Material::bulk_modulus, "Bulk modulus E / (3*(1-2*nu)).")
        .def_static("steel",     &Material::steel,     "Structural steel (E=210000, nu=0.30, rho=7850).")
        .def_static("aluminium", &Material::aluminium, "Aluminium alloy (E=70000, nu=0.33, rho=2700).")
        .def_static("concrete",  &Material::concrete,  "Concrete (E=30000, nu=0.20, rho=2400).");

    py::class_<ConstitutiveLaw, PyConstitutiveLaw, std::shared_ptr<ConstitutiveLaw>>(m, "ConstitutiveLaw",
        "Abstract B-free constitutive law (Planas, Romero & Sancho 2012). "
        "Subclass from Python and override n_dofs_per_cp() -> int and "
        "stiffness_density(grad_a, grad_b, x_phys) -> ndarray (n x n). "
        "PlaneStress and PlaneStrain are the built-in concrete subclasses.")
        .def(py::init<const Material&>(), py::arg("material"),
            "Base-class constructor for Python subclasses. Prefer PlaneStress(mat) "
            "or PlaneStrain(mat) unless implementing a custom law.")
        .def("n_dofs_per_cp", &ConstitutiveLaw::n_dofs_per_cp,
            "Number of displacement DOFs per control point (2 for plane problems, 3 for 3D solid).")
        .def("stiffness_density", &ConstitutiveLaw::stiffness_density,
            py::arg("grad_a"), py::arg("grad_b"), py::arg("x_phys"),
            "Elementary stiffness block K^{ab} for one Gauss point.\n"
            "grad_a / grad_b: physical-space gradients of basis functions a and b.\n"
            "x_phys: physical coordinates of the Gauss point.\n"
            "Returns an (n x n) matrix where n = n_dofs_per_cp().")
        .def_property_readonly("material", &ConstitutiveLaw::material,
            "The Material this law was constructed with.");

    py::class_<PlaneStress, ConstitutiveLaw, std::shared_ptr<PlaneStress>>(m, "PlaneStress",
        "2D plane-stress constitutive law (sigma_33 = 0). "
        "lambda_eff = nu*E / (1 - nu^2), mu = E / (2*(1+nu)).")
        .def(py::init<const Material&>(), py::arg("material"),
            "Construct from a Material. "
            "Use for thin structures where the out-of-plane stress is zero.");

    py::class_<PlaneStrain, ConstitutiveLaw, std::shared_ptr<PlaneStrain>>(m, "PlaneStrain",
        "2D plane-strain constitutive law (eps_33 = 0). "
        "lambda_eff = nu*E / ((1+nu)*(1-2*nu)), mu = E / (2*(1+nu)).")
        .def(py::init<const Material&>(), py::arg("material"),
            "Construct from a Material. "
            "Use for thick cross-sections where the out-of-plane strain is zero.");

    py::class_<Traction, PyTraction, std::shared_ptr<Traction>>(m, "Traction",
        "Base class for a distributed boundary load (force per unit length), "
        "evaluated at a physical point. Subclass from Python and override "
        "evaluate(physical_point) -> np.ndarray to define position-dependent "
        "tractions. ConstantTraction is the built-in concrete subclass for "
        "uniform loads.")
        .def(py::init<>(), "Construct a Traction base object; subclass and override evaluate().")
        .def("evaluate", &Traction::evaluate, py::arg("physical_point"),
            "Override in a Python subclass: return the traction vector (force per unit "
            "length) at physical_point as a 1D NumPy array of length dim_phys. "
            "PatchIntegrator.integrate_boundary_load() calls this once per Gauss point "
            "on the loaded edge.");

    py::class_<ConstantTraction, Traction, std::shared_ptr<ConstantTraction>>(m, "ConstantTraction",
        "A uniform traction vector (tx, ty), constant over the whole "
        "loaded edge/span-range.")
        .def(py::init<const Eigen::Vector2d&>(), py::arg("value"),
            "Construct with a constant traction vector value (shape (2,) for 2D).");

    py::class_<BoundaryLoadSpec>(m, "BoundaryLoadSpec",
        "One distributed boundary load to assemble over a whole "
        "PatchAssembly (see PatchIntegrator.assemble_boundary_load()). "
        "Unlike `materials`/`operators` (one entry per patch), a boundary "
        "load only applies to a specific patch/edge, so specs are given as "
        "a sparse list instead.")
        .def(py::init([](size_t patch_index, int direction, int side,
                          std::shared_ptr<Traction> traction, int span_min, int span_max) {
                BoundaryLoadSpec spec;
                spec.patch_index = patch_index;
                spec.direction = direction;
                spec.side = side;
                spec.traction = traction;
                spec.span_min = span_min;
                spec.span_max = span_max;
                return spec;
            }),
            py::arg("patch_index"), py::arg("direction"), py::arg("side"),
            py::arg("traction"), py::arg("span_min") = -1, py::arg("span_max") = -1,
            "Construct a boundary load spec. patch_index: 0-based index in the assembly. "
            "direction: fixed parametric direction (0=u, 1=v). side: 0=first, 1=last. "
            "traction: Traction object (e.g. ConstantTraction). "
            "span_min/span_max: restrict to a sub-range of the edge (-1/-1 = whole edge).")
        .def_readwrite("patch_index", &BoundaryLoadSpec::patch_index,
            "0-based index of the patch this load applies to, in PatchAssembly.get_patchs() order.")
        .def_readwrite("direction", &BoundaryLoadSpec::direction,
            "Fixed parametric direction of the loaded edge (0=u, 1=v).")
        .def_readwrite("side", &BoundaryLoadSpec::side,
            "0 for the first boundary (u=0 or v=0), 1 for the last (u=1 or v=1).")
        .def_readwrite("traction", &BoundaryLoadSpec::traction,
            "Traction object evaluated at each Gauss point along the loaded edge.")
        .def_readwrite("span_min", &BoundaryLoadSpec::span_min,
            "First knot-span index (in the non-fixed direction) of the loaded sub-range. "
            "-1 means start from the first non-degenerate span (default: load whole edge).")
        .def_readwrite("span_max", &BoundaryLoadSpec::span_max,
            "Last knot-span index (inclusive) of the loaded sub-range. "
            "-1 means end at the last non-degenerate span (default: load whole edge).");

    py::class_<PatchIntegrator>(m, "PatchIntegrator",
        "Gauss-quadrature integrator over a single 2D patch. Assembles stiffness "
        "and mass matrices, boundary load vectors, and custom scalar integrals via "
        "a ConstitutiveLaw. The integration dispatches at construction time between "
        "a B-spline and a NURBS code path (if constexpr), so pure B-spline patches "
        "run with no rational-basis overhead. "
        "For multi-patch assemblies use the static assemble_*() class methods.")
        .def(py::init<const Patch&, const IGABasis1D&, const IGABasis1D&,
                      std::shared_ptr<const ConstitutiveLaw>>(),
             py::arg("patch"), py::arg("basis_u"), py::arg("basis_v"), py::arg("law"),
             "Construct for the given patch and precomputed Gauss data (one IGABasis1D "
             "per parametric direction). law defines the mechanical behaviour and the "
             "number of DOFs per control point.")
        .def("integrate_stiffness", &PatchIntegrator::integrateStiffness,
            "Assemble the stiffness matrix of this single patch. "
            "Returns a scipy.sparse.csc_matrix of size (n_dof, n_dof) where "
            "n_dof = n_cp * law.n_dofs_per_cp(). "
            "For multi-patch assemblies use assemble_stiffness() instead.")
        .def("integrate_mass", &PatchIntegrator::integrateMass,
            "Same as integrate_stiffness(), but for the consistent mass matrix of "
            "this single patch. Requires law.material().rho > 0.")
        .def_static("assemble_stiffness", &PatchIntegrator::assembleStiffness,
            py::arg("assembly"), py::arg("laws"), py::arg("gauss_n") = 0,
            "Assemble the global stiffness matrix of a whole PatchAssembly.\n"
            "laws must have one entry per patch, in assembly.get_patchs() "
            "(add_patch()) order. Control points shared between patches "
            "(see PatchAssembly) automatically resolve to the same global dof, "
            "so contributions from every patch touching a shared boundary are "
            "summed into the same matrix entry -- no special-casing needed. "
            "gauss_n: Gauss points per span per direction for every patch; if "
            "0 (default), each direction of each patch uses its own degree + 1.")
        .def_static("assemble_mass", &PatchIntegrator::assembleMass,
            py::arg("assembly"), py::arg("laws"), py::arg("gauss_n") = 0,
            "Same as assemble_stiffness(), but for the consistent mass matrix of "
            "the whole PatchAssembly. Every entry in laws must have material().rho > 0.")
        .def("integrate_operator", &PatchIntegrator::integrateOperator, py::arg("op"),
            "Integrate a custom LocalOperator over this single patch. For "
            "development/testing: subclass LocalOperator in Python and override "
            "compute_integrand(). Not intended for performance-critical "
            "assembly -- prefer integrate_stiffness()/integrate_mass() for that.")
        .def_static("assemble_operator", &PatchIntegrator::assembleOperator,
            py::arg("assembly"), py::arg("operators"), py::arg("gauss_n") = 0,
            "Same as assemble_stiffness()/assemble_mass(), but for a custom "
            "LocalOperator. operators must have one entry per patch, in "
            "add_patch() order.")
        .def("integrate_boundary_load", &PatchIntegrator::integrateBoundaryLoad,
            py::arg("direction"), py::arg("side"), py::arg("traction"),
            py::arg("span_min") = -1, py::arg("span_max") = -1,
            "Integrate a distributed boundary load (force per unit length) "
            "over this single patch's edge obtained by fixing `direction` "
            "at its first (side=0) or last (side=1) span. span_min/span_max "
            "(raw knot-span indices of the OTHER direction, like "
            "Patch.boundary_control_points()) restrict integration to a "
            "sub-range of the edge; -1/-1 (default) integrates the whole "
            "edge. Returns a vector sized like integrate_stiffness()'s "
            "matrix, with nonzero entries only at dofs of control points on "
            "the loaded edge/span-range.")
        .def_static("assemble_boundary_load", &PatchIntegrator::assembleBoundaryLoad,
            py::arg("assembly"), py::arg("specs"), py::arg("gauss_n") = 0,
            "Assemble several distributed boundary loads (each on a "
            "specific patch/edge, via a list of BoundaryLoadSpec) over a "
            "whole PatchAssembly. Sized to the assembly-wide total dof "
            "count -- same size as assemble_stiffness()/assemble_mass(), "
            "regardless of which patches the specs actually touch -- so the "
            "result can be added directly to a stiffness/mass right-hand "
            "side.")
        .def("integrate_scalar_operator", &PatchIntegrator::integrateScalarOperator,
            py::arg("op"), py::arg("u_global"),
            "Integrate a ScalarLocalOperator over this single patch and return "
            "the accumulated scalar. At each Gauss point the operator receives: "
            "R, dRdx, dRdy (rationalized NURBS basis values and physical-space "
            "gradients), physical_point (x,y), and u_local "
            "(local DOF values in [u_x0,u_y0,u_x1,u_y1,...] order). Typical "
            "use: error norms and energy functionals that depend on the current "
            "FE solution. Subclass ScalarLocalOperator in Python and override "
            "compute_scalar_integrand().");

    py::class_<ScalarLocalOperator, PyScalarLocalOperator,
               std::shared_ptr<ScalarLocalOperator>>(m, "ScalarLocalOperator",
        "Base class for a scalar-valued integration term (error norms, energy "
        "functionals, etc.) to be subclassed FROM PYTHON. "
        "PatchIntegrator.integrate_scalar_operator() provides the Gauss loop, "
        "Jacobian, physical-space gradients, physical coordinates, and local "
        "DOF values; the subclass only needs to implement the scalar integrand "
        "at one Gauss point. "
        "Override compute_scalar_integrand(R, dRdx, dRdy, physical_point, "
        "u_local) -> float, where u_local has layout "
        "[u_x_cp0, u_y_cp0, u_x_cp1, u_y_cp1, ...] (2 DOFs per CP for 2D "
        "problems). PatchIntegrator multiplies the returned value by the Gauss "
        "weight and abs(detJ) before summing.")
        .def(py::init<>(), "Construct; subclass and override compute_scalar_integrand().")
        .def("compute_scalar_integrand", &ScalarLocalOperator::computeScalarIntegrand,
             py::arg("R"), py::arg("dRdx"), py::arg("dRdy"),
             py::arg("physical_point"), py::arg("u_local"),
             "Override in a Python subclass: return the scalar integrand at one Gauss point. "
             "R, dRdx, dRdy: lists (length n_active) of basis values and physical-space gradients. "
             "physical_point: (x, y) physical coordinates (shape (2,)). "
             "u_local: flat DOF values [u_x0, u_y0, u_x1, u_y1, ...] for the active control points. "
             "The returned float is multiplied by |detJ| * gauss_weight and accumulated over the patch.");

    py::class_<LocalOperator, PyLocalOperator, std::shared_ptr<LocalOperator>>(m, "LocalOperator",
        "Base class to subclass FROM PYTHON for a custom integration term. "
        "PatchIntegrator does the Gauss-point loop and the Jacobian/gradient "
        "computation -- override compute_integrand(R, dRdx, dRdy) to return "
        "just the term to integrate at one Gauss point (e.g. B^T*D*B for "
        "stiffness, rho*N^T*N for mass), built from the basis values R and "
        "physical gradients dRdx, dRdy of this span's active basis functions. "
        "PatchIntegrator multiplies the returned matrix by the Gauss weight "
        "and abs(detJ) and sums it over every Gauss point of every span. "
        "Intended for development/testing convenience, not performance: every "
        "Gauss point triggers one Python call.")
        .def(py::init<>(), "Construct; subclass and override compute_integrand().")
        .def("compute_integrand", &LocalOperator::computeIntegrand,
             py::arg("R"), py::arg("dRdx"), py::arg("dRdy"),
             "Override in a Python subclass: return the local matrix contribution at one "
             "Gauss point. R, dRdx, dRdy are lists (length n_active) of basis values and "
             "their physical-space gradients for the active basis functions in this span "
             "(u-fastest order). The returned (n x n) matrix is multiplied by |detJ| * "
             "gauss_weight and summed over every Gauss point of every span.");

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
        .def(py::init<>(), "Construct an empty assembly; call add_patch() for each patch.")
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
            },
            "Return the list of detected interfaces (one dict per compatible shared edge). "
            "Each dict has keys: "
            "patch_a/patch_b (patch indices), "
            "direction_a/direction_b (fixed parametric direction of the shared edge on each patch), "
            "side_a/side_b (0=first, 1=last), "
            "varying_direction_a/varying_direction_b (direction to refine on each patch "
            "to propagate refinement across this interface), "
            "reversed (bool: True if the shared edge is traversed in opposite senses on "
            "the two patches, e.g. u of patch_a glued to v of patch_b). "
            "Empty until detect_interfaces() has been called.")
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
            "Return the number of control points of patches[patch_index] "
            "(= product of local_shape entries).")
        .def("apply_transformation_to_control_points", [](const PatchAssembly& self, const std::vector<Eigen::VectorXd>& old_control_points, size_t patch_index) {
                return self.applyTransformationToControlPoints(old_control_points, patch_index);
            }, py::arg("old_control_points"), py::arg("patch_index"),
            "Apply patches[patch_index]'s stored transition matrix to a list of "
            "pre-refinement control points and return the post-refinement list. "
            "Each element of old_control_points is a 1D array of physical coordinates "
            "for one control point. Requires set_transformation_matrix() to have been "
            "called first (or a refiner's transition matrix to be stored)."
        )
        .def("apply_transformation_to_dofs", [](const PatchAssembly& self, const Eigen::VectorXd& old_dofs, size_t patch_index, size_t dim_phys) {
                return self.applyTransformationToDOFs(old_dofs, patch_index, dim_phys);
            }, py::arg("old_dofs"), py::arg("patch_index"), py::arg("dim_phys"),
            "Apply patches[patch_index]'s stored transition matrix to a flat "
            "pre-refinement DOF vector (length n_old_cp * dim_phys) and return "
            "the post-refinement DOF vector (length n_new_cp * dim_phys). "
            "Useful for warm-starting a refined problem from an existing solution. "
            "Requires set_transformation_matrix() to have been called first."
        )
        .def("set_transformation_matrix", [](PatchAssembly& self, size_t patch_index, const Eigen::MatrixXd matrix) {
                    return self.setTransformationMatrix(patch_index, matrix);
            }, py::arg("patch_index"), py::arg("matrix"),
            "Store a refinement transition matrix for patches[patch_index]. "
            "Typically the matrix returned by a refiner (e.g. HRefiner.refine(), "
            "SubdivisionRefiner.refine()), shape (n_new_cp, n_old_cp). "
            "Enables apply_transformation_to_control_points() and "
            "apply_transformation_to_dofs() for this patch. "
            "Defaults to identity (no transformation stored) until set."
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

    py::class_<RefinementOperator, std::shared_ptr<RefinementOperator>>(m, "RefinementOperator",
        "Abstract base class for in-place patch refinement. Concrete subclasses: "
        "HRefiner (single knot insertion), SubdivisionRefiner (uniform bisection), "
        "PRefiner (degree elevation). Each refine() call updates the patch's "
        "BSplineTensor and ControlPointManager in-place and returns the transition "
        "matrix (n_new_cp x n_old_cp).")
        .def("get_type", &RefinementOperator::getType,
            "Return a short string identifying the refiner type ('h', 'subdivision', 'p').")
        .def("refine", [](const RefinementOperator& self, Patch& patch) {
            Eigen::MatrixXd T;
            self.refine(patch, T);
            return T;
        }, py::arg("patch"),
           "Refine patch in-place. Returns transition_matrix (nb_new_cp x nb_old_cp).");

    py::class_<HRefiner, RefinementOperator, std::shared_ptr<HRefiner>>(m, "HRefiner",
        "Single knot insertion refiner (h-refinement). Inserts one knot value into "
        "one parametric direction via Boehm's algorithm (Piegl & Tiller A5.1). "
        "Chain multiple calls to insert several knots; for uniform refinement "
        "SubdivisionRefiner is more efficient.")
        .def(py::init<int, double>(), py::arg("direction"), py::arg("knot"),
            "direction: parametric direction (0=u, 1=v). "
            "knot: the value to insert (must lie within the existing knot span).")
        .def("refine", [](const HRefiner& self, Patch& patch) {
            Eigen::MatrixXd T;
            self.refine(patch, T);
            return T;
        }, py::arg("patch"),
           "Insert knot in-place (h-refinement). Returns transition_matrix (nb_new_cp x nb_old_cp).")
        .def("refine_1d", [](const HRefiner& self, Patch& patch,
                              const std::unordered_set<size_t>& protected_global_ids) {
            Eigen::MatrixXd T;
            self.refine_1d(patch, T, protected_global_ids);
            return T;
        }, py::arg("patch"), py::arg("protected_global_ids") = std::unordered_set<size_t>{},
           "Fast path: insert one knot in-place, return only the 1D transition matrix "
           "(n_new_1d x n_old_1d). Use nd_transition_from_1d() to reconstruct the full "
           "nD matrix if needed. "
           "protected_global_ids: ids borrowed from another patch (shared interface) that "
           "must never be recomputed/renumbered; if empty (default), this patch is assumed "
           "to exclusively own its cp_manager range. "
           "Use as refine_1d_fn in PatchAssembly.refine_with_propagation(): "
           "lambda patch, d, ids: HRefiner(d, knot).refine_1d(patch, ids).");

    py::class_<SubdivisionRefiner, RefinementOperator, std::shared_ptr<SubdivisionRefiner>>(m, "SubdivisionRefiner",
        "Uniform knot insertion refiner: bisects every existing knot span n_levels times "
        "in one parametric direction (element count multiplied by 2^n_levels). "
        "refine_1d() exposes a fast path returning only the 1D transition matrix, "
        "suitable for PatchAssembly.refine_with_propagation().")
        .def(py::init<int, int>(), py::arg("direction"), py::arg("n_levels") = 1,
            "direction: parametric direction (0=u, 1=v). "
            "n_levels: number of bisection levels (default 1, i.e. double the element count).")
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

    py::class_<PRefiner, RefinementOperator, std::shared_ptr<PRefiner>>(m, "PRefiner",
        "Degree elevation refiner (p-refinement). Elevates the polynomial degree by "
        "n_elevations steps in one parametric direction (Piegl & Tiller A5.9). "
        "Typically combined with knot insertion for k-refinement: elevate degree first, "
        "then insert knots to recover continuity control.")
        .def(py::init<int, int>(), py::arg("direction"), py::arg("n_elevations") = 1,
            "direction: parametric direction (0=u, 1=v). "
            "n_elevations: number of degree elevation steps (default 1).")
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
        .def_readonly("C", &BezierElementND::C,
            "Extraction matrix, shape (n_local, n_local): N_active(xi) = C @ B_nd(xi_hat) "
            "where n_local = prod_d(p_d + 1).")
        .def_property_readonly("active_indices",
            [](const BezierElementND& self) { return self.active; },
            "u-fastest flat CP indices of the active B-spline basis functions on this element.")
        .def_property_readonly("elem_index",
            [](const BezierElementND& self) { return self.elem_index; },
            "0-based element multi-index (e_0, e_1, ...) of this element in the patch grid.");

    py::class_<PatchEvaluator>(m, "PatchEvaluator",
        R"doc(
Evaluates a FE solution field at arbitrary parametric-space points.

The patch must carry a PatchDOFManager (i.e. it was constructed with one, or
one was attached before building this evaluator).  NURBS patches are handled
exactly — the rational basis R_a = w_a N_a / W is applied before accumulating
DOF contributions, using the same zero-overhead dispatch as PatchIntegrator.
The evaluation loop is OpenMP-parallel over evaluation points.

Parameters
----------
patch : Patch
    Patch with a PatchDOFManager attached.
        )doc")
        .def(py::init<const Patch&>(), py::arg("patch"),
            "Construct for the given patch. The patch must already carry a "
            "PatchDOFManager (attached at construction or via "
            "PatchAssembly.update_dof_managers()).")
        .def("evaluate_solution", &PatchEvaluator::evaluateSolutionOMP,
             py::arg("params"), py::arg("u_global"),
             R"doc(
Evaluate the FE solution at parametric-space points.

Parameters
----------
params : ndarray, shape (n_pts, n_param_dims)
    Parameter values at which to evaluate.
u_global : ndarray, shape (n_dof,)
    Global solution vector (same size as the stiffness matrix).

Returns
-------
ndarray, shape (n_pts, n_dofs_per_cp)
    Field values at each evaluation point.
    For a 2-D solid with 2 DOFs per CP: columns are [u_x, u_y].
             )doc");

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
        .def(py::init<int>(), py::arg("direction"),
            "Construct for the given parametric direction (0=u, 1=v, 2=w).")
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
