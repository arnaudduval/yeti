#pragma once

#include <vector>
#include <unordered_map>
#include <unordered_set>
#include "Patch.hpp"

class PatchAssembly {
private:
    // List of assembly patchs
    std::vector<std::shared_ptr<Patch>> patches_;

    // Structure storing connexions between patches
    // Key: global index of shared control point
    // Value: list of indices of patchs sharing this point
    std::unordered_map<size_t, std::vector<size_t>> shared_control_points_;

    // transformation matrices for each patch
    std::vector<Eigen::MatrixXd> transformation_matrices_;

public:
    PatchAssembly() = default;
    
    // Add a patch to assembly
    void addPatch(const std::shared_ptr<Patch>& patch) {
        patches_.push_back(patch);

        // Initialize identity matrix for this patch
        size_t nb_cp = patch->global_indices.size();
        transformation_matrices_.emplace_back(Eigen::MatrixXd::Identity(nb_cp, nb_cp));
    }

    // Detect shared control points between all patchs
    void detectSharedControlPoints() {
        shared_control_points_.clear();
        for (size_t i = 0; i < patches_.size(); ++i) {
            for (size_t j = i+1; j < patches_.size(); ++j) {
                const auto& patch_i = patches_[i];
                const auto& patch_j = patches_[j]; 
                std::unordered_set<size_t> global_indices_i(
                    patch_i->global_indices.begin(),
                    patch_i->global_indices.end()
                );
                for (size_t idx : patch_j->global_indices) {
                    if (global_indices_i.count(idx)) {
                        shared_control_points_[idx].push_back(i);
                        shared_control_points_[idx].push_back(j);
                    }
                }
            }
        }
    }

    // Return patchs sharing a given control point
    const std::vector<size_t>& getPatchsSharingControlPoints(size_t global_cp_idx) const {
        static const std::vector<size_t> empty;
        auto it = shared_control_points_.find(global_cp_idx);
        return (it != shared_control_points_.end()) ? it->second : empty;
    }

    // Return control points shared between 2 patchs
    std::vector<size_t> getSharedControlPoints(
        const Patch& patch1,
        const Patch& patch2
    ) const {
        std::vector<size_t> shared;
        std::unordered_set<size_t> indices_patch1(
            patch1.global_indices.begin(),
            patch2.global_indices.end()
        );
        for (size_t idx : patch2.global_indices) {
            if (indices_patch1.count(idx)) {
                shared.push_back(idx);
            }
        }
        return shared;
    }

    // Return number of control points for a given patch
    size_t getControlPointsForPatch(size_t patch_index) const {
        if (patch_index >= patches_.size()) {
            throw std::out_of_range("Invalid patch index.");
        }
        return patches_[patch_index]->global_indices.size();
    }

    // Apply transformation matrix to control points of a patch
    // Return new control points after refinement
    std::vector<Eigen::VectorXd> applyTransformationToControlPoints(
        const std::vector<Eigen::VectorXd>& old_control_points,
        size_t patch_index
    ) const {
        if (patch_index >= transformation_matrices_[patch_index].size()) {
            throw std::out_of_range("Invalid atch index for transformation matrix.");
        }
        const auto& matrix = transformation_matrices_[patch_index];

        size_t nb_new_cp = matrix.rows();
        size_t nb_old_cp = old_control_points.size();

        // Verify data coherence
        if (nb_old_cp != static_cast<size_t>(matrix.cols())) {
            throw std::runtime_error("Incompatible matrix size for control points transformation.");
        }

        std::vector<Eigen::VectorXd> new_control_points(nb_new_cp);
        for (size_t i = 0; i < nb_new_cp; ++i) {
            new_control_points[i] = Eigen::VectorXd::Zero(old_control_points[0].size());
            for (size_t j = 0; j < nb_old_cp; ++j) {
                new_control_points[i] += matrix(i, j) * old_control_points[j];
            }
        }
        return new_control_points;
    }

    // Apply transformation matrix to DOFs of a patch
    // Return new DOFs after refinement
    Eigen::VectorXd applyTransformationToDOFs(
        const Eigen::VectorXd& old_dofs,
        size_t patch_index,
        size_t dim_phys
    ) const {
        if (patch_index >= transformation_matrices_.size()) {
            throw std::out_of_range("Invalid patch index for transformation matrix.");
        }
        const auto matrix = transformation_matrices_[patch_index];

        size_t nb_old_cp = matrix.cols();
        size_t nb_new_cp = matrix.rows();

        // Verify size
        if (old_dofs.size() != nb_old_cp * dim_phys) {
            throw std::runtime_error("Incompatible DIFs size for transformation.");
        }

        Eigen::VectorXd new_dofs(nb_new_cp * dim_phys);
        new_dofs.setZero();

        for (size_t d = 0; d < dim_phys; ++d) {
            for (size_t i = 0; i < nb_new_cp; ++i) {
                for (size_t j = 0; j < nb_old_cp; ++j) {
                    new_dofs[i * dim_phys + d] += matrix(i, j) * old_dofs[j * dim_phys + d];
                }
            }
        }

        return new_dofs;
    }
    
    // Update transformation matrix for a given patch
    void setTransformationMatrix(size_t patch_index, const Eigen::MatrixXd& matrix) {
        if (patch_index >= transformation_matrices_.size()) {
            throw std::out_of_range("Invalid patch index for transformation matrix.");
        }
        transformation_matrices_[patch_index] = matrix;
    }



    // Getters
    const std::vector<std::shared_ptr<Patch>>& getPatchs() const { return patches_; }
    const std::unordered_map<size_t, std::vector<size_t>> getSharedControlPoints() const {
        return shared_control_points_;
    }
    const std::vector<Eigen::MatrixXd>& getTransformationMatrices() const {
        return transformation_matrices_;
    }

};