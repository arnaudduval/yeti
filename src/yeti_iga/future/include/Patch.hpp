#pragma once
#include <vector>
#include <memory>
#include <pybind11/numpy.h>
#include "BSpline.hpp"
#include "BSplineTensor.hpp"
#include "ControlPointManager.hpp"
#include "DOFManager.hpp"
#include "SpanNDIterator.hpp"

namespace py = pybind11;

struct Patch {
    BSplineTensor tensor;

    // mapping local -> global
    // Row-major indexing : [i_v * n_u + i_u]
    std::vector<size_t> global_indices;

    // dimensions nu, nv, nw (for indexing)
    std::vector<size_t> local_shape;

    std::shared_ptr<ControlPointManager> cp_manager;
    std::shared_ptr<PatchDOFManager> dof_manager;

    // Constructor without DOF manager (geometry only)
    Patch(const BSplineTensor& t,
          std::shared_ptr<ControlPointManager> mgr,
          const std::vector<size_t>& mapping,
          const std::vector<size_t>& local_shape_)
        : tensor(t), global_indices(mapping), local_shape(local_shape_), cp_manager(std::move(mgr)), dof_manager(nullptr) {}

    Patch(const BSplineTensor& t,
          std::shared_ptr<ControlPointManager> mgr,
          const std::vector<size_t>& mapping,
          const std::vector<size_t>& local_shape_,
          std::shared_ptr<PatchDOFManager> dof_mgr)
        : tensor(t), global_indices(mapping), local_shape(local_shape_), cp_manager(std::move(mgr)), dof_manager(dof_mgr) {}

    // get view (zero copy) to coordinates of local control point i_local
    double* local_cp_ptr(size_t i_local);
    const double* local_cp_ptr(size_t i_local) const;

    // Numpy 1D view on a local point
    py::array_t<double> local_cp_view(size_t i_local);

    // get a numpy view for all local control points
    // shape = local_shape + (dim_phys,)
    py::array_t<double> local_control_point_view() const;

    py::array_t<double> EvaluatePatchND(const py::array_t<int> spans,
                                        const py::array_t<double>& u) const;

    py::array_t<double> EvaluatePatchNDOMP(const py::array_t<int> spans,
                                           const py::array_t<double>& u) const;

    std::vector<const double*> control_points_for_span(const std::vector<int>& span) const;

    // Patch-LOCAL flat positions (u-fastest) of the control points on the
    // edge obtained by fixing `direction` at its first (side=0) or last
    // (side=1) local index. If span_min/span_max are given (both >= 0,
    // raw knot-span indices like the ones Patch::spans() yields), restricts
    // the selection to control points active over those spans of the
    // OTHER ("varying") direction; -1/-1 (default) selects the whole edge.
    // 2D patches only (Phase 1 scope, like PatchAssembly's interfaces).
    //
    // Returned positions are suitable for
    // PatchDOFManager::get_global_dof_indices()/get_global_dof() -- this is
    // the building block for specifying Dirichlet boundary conditions at
    // edge or boundary-span granularity (control-point granularity needs no
    // helper: any local position works directly with the DOF manager).
    std::vector<size_t> boundary_control_points(int direction, int side,
                                                 int span_min = -1, int span_max = -1) const;

    SpanNDIterator spans() const { return SpanNDIterator(tensor); }

    void Test();
};

