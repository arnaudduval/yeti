#pragma once
#include <algorithm>
#include <vector>
#include <memory>

// TODO : make a full restructuration of this file
// TODO use specific directory to separate IGA and pure geometry ?

// Handle global degrees of freedom for all control points
// Patchs can share degrees of freedom
// Different patchs can have different numbers of DOF by CP
struct GlobalDOFManager {
    std::vector<int> dofs_per_control_point;
    std::vector<size_t> global_dofs;
    std::vector<size_t> offsets;

    GlobalDOFManager(const std::vector<int>& dofs_per_control_point)
        : dofs_per_control_point(dofs_per_control_point) {
        offsets.resize(dofs_per_control_point.size() + 1);
        offsets[0] = 0;
        for (size_t i =0; i < dofs_per_control_point.size(); ++i) {
            offsets[i+1] = offsets[i] + dofs_per_control_point[i];
        }
        global_dofs.resize(offsets.back());
        for (size_t i = 0; i < global_dofs.size(); ++i) {
            global_dofs[i] = i;
        }
    }

    // Get DOF indices for a given control point
    std::vector<size_t> get_dof_indices(size_t control_point_idx) const {
        std::vector<size_t> dof_indices;
        size_t start = offsets[control_point_idx];
        size_t end = offsets[control_point_idx + 1];
        for (size_t i = start; i < end; ++i) {
            dof_indices.push_back(global_dofs[i]);
        }
        return dof_indices;
    }
};



// Structure handling degrees of freedom of a patch
struct PatchDOFManager {
    int dofs_per_control_point;
    std::vector<size_t> local_to_global_dofs;   // Mapping of local DOFs to global DOFs

    PatchDOFManager(int dofs_per_control_point, const std::vector<size_t>& control_points, const GlobalDOFManager& global_dof_manager)
    : dofs_per_control_point(dofs_per_control_point) {
        local_to_global_dofs.resize(control_points.size() * dofs_per_control_point);
        for (size_t i = 0; i < control_points.size(); ++i) {
            std::vector<size_t> global_dofs = global_dof_manager.get_dof_indices(control_points[i]);
            for (int j = 0 ; j < dofs_per_control_point; ++j) {
                local_to_global_dofs[i*dofs_per_control_point + j] = global_dofs[j];
            }
        }
    }

    // Build a refined DOF manager: copy old mapping, assign new consecutive DOFs for extra CPs
    PatchDOFManager(const PatchDOFManager& old_manager, const std::vector<size_t>& new_global_indices)
        : dofs_per_control_point(old_manager.dofs_per_control_point) {
        local_to_global_dofs.resize(new_global_indices.size() * dofs_per_control_point);
        size_t nb_old = old_manager.local_to_global_dofs.size();
        for (size_t i = 0; i < nb_old; ++i)
            local_to_global_dofs[i] = old_manager.local_to_global_dofs[i];
        size_t next_dof = nb_old == 0 ? 0 :
            *std::max_element(old_manager.local_to_global_dofs.begin(),
                              old_manager.local_to_global_dofs.end()) + 1;
        for (size_t i = nb_old; i < local_to_global_dofs.size(); ++i)
            local_to_global_dofs[i] = next_dof++;
    }

    // Get global dofs indices of a given local control point
    std::vector<size_t> get_global_dof_indices(size_t local_control_point_idx) const {
        std::vector<size_t> dof_indices(dofs_per_control_point);
        for (int i = 0; i < dofs_per_control_point; ++i) {
            dof_indices[i] = local_to_global_dofs[local_control_point_idx * dofs_per_control_point + i];
        }
        return dof_indices;
    }
};