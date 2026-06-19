#pragma once
#include <algorithm>
#include <stdexcept>
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
        if (control_point_idx + 1 >= offsets.size())
            throw std::out_of_range(
                "GlobalDOFManager::get_dof_indices: control_point_idx out of range "
                "(this manager was sized for " + std::to_string(offsets.size() - 1) +
                " control points; call grow() after refining a shared assembly).");
        std::vector<size_t> dof_indices;
        size_t start = offsets[control_point_idx];
        size_t end = offsets[control_point_idx + 1];
        for (size_t i = start; i < end; ++i) {
            dof_indices.push_back(global_dofs[i]);
        }
        return dof_indices;
    }

    // Number of control points currently covered (indices 0..n-1 are valid
    // for get_dof_indices()).
    size_t n_control_points() const { return dofs_per_control_point.size(); }

    // Grow to cover `new_n_cp` control points (no-op if already >= that).
    // Every newly covered control point gets `dofs_per_cp` fresh, globally
    // unique dofs, appended after the highest dof currently in use.
    //
    // Call this BEFORE PatchAssembly::compact(): compact() renumbers control
    // point ids, which would silently desynchronize this manager's
    // cp-id-indexed mapping if dofs were assigned first against stale ids.
    void grow(size_t new_n_cp, int dofs_per_cp) {
        size_t old_n_cp = dofs_per_control_point.size();
        if (new_n_cp <= old_n_cp) return;

        dofs_per_control_point.resize(new_n_cp, dofs_per_cp);
        offsets.resize(new_n_cp + 1);
        for (size_t i = old_n_cp; i < new_n_cp; ++i)
            offsets[i + 1] = offsets[i] + dofs_per_control_point[i];

        size_t next_dof = global_dofs.empty() ? 0 :
            *std::max_element(global_dofs.begin(), global_dofs.end()) + 1;
        global_dofs.resize(offsets.back());
        for (size_t i = offsets[old_n_cp]; i < global_dofs.size(); ++i)
            global_dofs[i] = next_dof++;
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

    // Get global dofs indices of a given local control point.
    // Checked (throws on out-of-range) — use in high-level / Python-facing code.
    std::vector<size_t> get_global_dof_indices(size_t local_control_point_idx) const {
        if ((local_control_point_idx + 1) * dofs_per_control_point > local_to_global_dofs.size())
            throw std::out_of_range(
                "PatchDOFManager::get_global_dof_indices: local_control_point_idx out of range.");
        std::vector<size_t> dof_indices(dofs_per_control_point);
        for (int i = 0; i < dofs_per_control_point; ++i) {
            dof_indices[i] = local_to_global_dofs[local_control_point_idx * dofs_per_control_point + i];
        }
        return dof_indices;
    }

    // Single-DOF accessor, no heap allocation — use in hot assembly loops
    // (e.g. PatchIntegrator::assembleLocalContribution). No bounds check:
    // caller (internal C++ code with already-validated indices) is trusted.
    size_t get_global_dof(size_t local_control_point_idx, int local_dof_idx) const {
        return local_to_global_dofs[local_control_point_idx * dofs_per_control_point + local_dof_idx];
    }
};