#pragma once
#include <vector>
#include <memory>


// Structure handling degrees of freedom of a patch
// TODO use sepcific directory to separate IGA and pure geometry ?
struct PatchDOFManager {
    int dofs_per_control_point;
    std::vector<size_t> local_dofs;
    std::vector<size_t> global_dofs;

    PatchDOFManager(int dofs_per_control_point, size_t n_control_points, size_t global_dof_offset) : dofs_per_control_point(dofs_per_control_point) {
        local_dofs.resize(n_control_points * dofs_per_control_point);
        global_dofs.resize(n_control_points * dofs_per_control_point);
        for (size_t i = 0; i < local_dofs.size(); ++i) {
            local_dofs[i] = i;
            global_dofs[i] = global_dof_offset + i;
        }
    }

    // Get local dofs indices of a given control point
    std::vector<size_t> get_local_dof_indices(size_t control_point_idx) const {
        std::vector<size_t> dof_indices(dofs_per_control_point);
        for (int i = 0; i < dofs_per_control_point; ++i) {
            dof_indices[i] = local_dofs[control_point_idx * dofs_per_control_point + i];
        }
        return dof_indices;
    }

    // Get global dofs indices of a given control point
    std::vector<size_t> get_global_dof_indices(size_t control_point_idx) const {
        std::vector<size_t> dof_indices(dofs_per_control_point);
        for (int i = 0; i < dofs_per_control_point; ++i) {
            dof_indices[i] = global_dofs[control_point_idx * dofs_per_control_point + i];
        }
        return dof_indices;
    }
};