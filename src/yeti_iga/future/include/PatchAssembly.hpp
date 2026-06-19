#pragma once

#include <algorithm>
#include <functional>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <stdexcept>
#include "Patch.hpp"

// Describes a compatible interface between two 2D patches: a shared edge,
// identified by which (direction, side) facet of each patch carries it.
//
// "direction_a"/"side_a": the facet of patch_a is obtained by fixing
// direction_a's local index at 0 (side_a == 0) or at its last index
// (side_a == 1); the facet's control points vary over the OTHER direction
// (varying_direction_a) -- this is the direction to refine on patch_a to
// grow the shared edge.
//
// "reversed": true if the shared control points, walked in increasing
// varying_direction_a order on patch_a, appear in DECREASING
// varying_direction_b order on patch_b (e.g. u of patch_a glued to v of
// patch_b but parametrized in opposite senses).
//
// Phase 1 scope: only 2D patches (tensor.components.size() == 2) are
// detected -- a 3D face interface has two varying directions and isn't
// handled yet.
struct PatchInterface {
    size_t patch_a, patch_b;
    int direction_a, side_a;
    int direction_b, side_b;
    int varying_direction_a;
    int varying_direction_b;
    bool reversed;
};

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

    // Compatible interfaces (shared edges) detected between patches.
    std::vector<PatchInterface> interfaces_;

    // Ordered list of global ids on the facet of `patch` obtained by fixing
    // `direction` at index 0 (side == 0) or at its last index (side == 1).
    // Order follows increasing index of the other (varying) direction.
    // 2D patches only (see PatchInterface).
    static std::vector<size_t> extractFacet(const Patch& patch, int direction, int side) {
        size_t ndim = patch.tensor.components.size();
        if (ndim != 2)
            throw std::invalid_argument("extractFacet: only 2D patches are supported (Phase 1 scope).");

        int varying = 1 - direction;
        size_t fixed_index = (side == 0) ? 0 : patch.local_shape[direction] - 1;

        std::vector<size_t> stride(ndim);
        stride[0] = 1;
        for (size_t d = 1; d < ndim; ++d)
            stride[d] = stride[d-1] * patch.local_shape[d-1];

        std::vector<size_t> facet;
        facet.reserve(patch.local_shape[varying]);
        for (size_t i = 0; i < patch.local_shape[varying]; ++i) {
            size_t flat = fixed_index * stride[direction] + i * stride[varying];
            facet.push_back(patch.global_indices[flat]);
        }
        return facet;
    }

public:
    PatchAssembly() = default;

    // Add a patch to assembly
    void addPatch(const std::shared_ptr<Patch>& patch) {
        patches_.push_back(patch);

        // Initialize identity matrix for this patch
        size_t nb_cp = patch->global_indices.size();
        transformation_matrices_.emplace_back(Eigen::MatrixXd::Identity(nb_cp, nb_cp));
    }

    // Reclaim memory orphaned by repeated refinements (see
    // HRefiner::apply_1d_cp_update / refine_1d's protected_global_ids: each
    // refinement of a shared patch relocates its private CPs to a fresh
    // dense block, leaving their old slots unused but still present in
    // cp_manager).
    //
    // Call this explicitly once all the refinements you need are done — not
    // after every refinement step, since it rewrites the whole shared CP
    // pool and every patch's global_indices.
    //
    // Walks all patches (in the order they were added), keeps each CP's
    // FIRST-encountered global id, assigns it a fresh dense id, and reuses
    // that same fresh id for every later occurrence (i.e. for CPs shared
    // with an already-processed patch). DOF managers are untouched: they
    // map local position -> global dof, independent of CP global ids.
    void compact() {
        if (patches_.empty()) return;
        auto cp_manager = patches_[0]->cp_manager;
        size_t dim = cp_manager->dim_phys;

        std::unordered_map<size_t, size_t> old_to_new;
        std::vector<double> new_coords;
        new_coords.reserve(cp_manager->coords.size());

        for (auto& patch : patches_) {
            for (size_t& gid : patch->global_indices) {
                size_t old_gid = gid;
                auto it = old_to_new.find(old_gid);
                size_t new_id;
                if (it == old_to_new.end()) {
                    new_id = new_coords.size() / dim;
                    const double* src = cp_manager->coords.data() + old_gid * dim;
                    new_coords.insert(new_coords.end(), src, src + dim);
                    old_to_new.emplace(old_gid, new_id);
                } else {
                    new_id = it->second;
                }
                gid = new_id;
            }
        }

        {
            std::lock_guard<std::mutex> lock(cp_manager->mtx);
            cp_manager->coords = std::move(new_coords);
        }

        // Old global ids are no longer valid -- caller must call
        // detectSharedControlPoints() / detectInterfaces() again if needed.
        shared_control_points_.clear();
        interfaces_.clear();
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

    // Detect compatible interfaces (shared edges) between all pairs of
    // patches: for every pair of facets (one per patch) whose control point
    // SETS match exactly, record a PatchInterface describing which
    // direction/side carries the edge on each side, and whether traversal
    // order is reversed between the two.
    //
    // Throws if two facets share the same set of control points but neither
    // direct nor reversed traversal order matches -- this violates the
    // "compatible boundary" assumption (same knot vector/degree/CPs).
    //
    // 2D patches only (see PatchInterface); patches with more than 2
    // parametric directions are silently skipped.
    void detectInterfaces() {
        interfaces_.clear();
        for (size_t i = 0; i < patches_.size(); ++i) {
            const Patch& A = *patches_[i];
            if (A.tensor.components.size() != 2) continue;

            for (size_t j = i + 1; j < patches_.size(); ++j) {
                const Patch& B = *patches_[j];
                if (B.tensor.components.size() != 2) continue;

                std::unordered_set<size_t> a_ids(A.global_indices.begin(), A.global_indices.end());
                bool any_shared = false;
                for (size_t idx : B.global_indices) {
                    if (a_ids.count(idx)) { any_shared = true; break; }
                }
                if (!any_shared) continue;

                for (int da = 0; da < 2; ++da) {
                    for (int sa = 0; sa < 2; ++sa) {
                        std::vector<size_t> facet_a = extractFacet(A, da, sa);
                        std::unordered_set<size_t> set_a(facet_a.begin(), facet_a.end());

                        for (int db = 0; db < 2; ++db) {
                            for (int sb = 0; sb < 2; ++sb) {
                                std::vector<size_t> facet_b = extractFacet(B, db, sb);
                                if (facet_b.size() != facet_a.size()) continue;

                                std::unordered_set<size_t> set_b(facet_b.begin(), facet_b.end());
                                if (set_a != set_b) continue;

                                bool same_order = (facet_a == facet_b);
                                std::vector<size_t> facet_b_rev(facet_b.rbegin(), facet_b.rend());
                                bool reversed_order = (facet_a == facet_b_rev);

                                if (!same_order && !reversed_order) {
                                    throw std::runtime_error(
                                        "Incompatible boundary ordering between patches "
                                        + std::to_string(i) + " and " + std::to_string(j) +
                                        ": shared control points don't match in either "
                                        "direct or reversed order.");
                                }

                                PatchInterface interface;
                                interface.patch_a = i;
                                interface.direction_a = da;
                                interface.side_a = sa;
                                interface.varying_direction_a = 1 - da;
                                interface.patch_b = j;
                                interface.direction_b = db;
                                interface.side_b = sb;
                                interface.varying_direction_b = 1 - db;
                                interface.reversed = reversed_order;
                                interfaces_.push_back(interface);
                            }
                        }
                    }
                }
            }
        }
    }

    // Signature of a 1D refinement step: refine `patch` in-place along the
    // given direction, never touching/renumbering the given protected ids.
    // The direction argument matters: on a crossed interface (e.g. u of
    // patch A glued to v of patch B), the neighbor must be refined along
    // its OWN matching direction, which generally differs from the
    // direction passed to refineWithPropagation() for the initial patch.
    // Matches a thin wrapper around SubdivisionRefiner::refine_1d /
    // PRefiner::refine_1d (modulo the discarded transition matrix, which
    // propagation doesn't need).
    using Refine1DFn = std::function<void(Patch&, int, const std::unordered_set<size_t>&)>;

    // Refine patches_[patch_index] along `direction`, automatically
    // propagating the SAME refinement to every neighbor connected through a
    // detected PatchInterface whose varying direction matches, then merging
    // the boundary control points created independently on both sides
    // (they are guaranteed to be coordinate-identical: same algorithm, same
    // input boundary CPs, deterministic).
    //
    // Requires detectInterfaces() (and detectSharedControlPoints()) to have
    // been called first. Does NOT touch any PatchDOFManager: the merged
    // boundary CPs on each side currently keep their own, unrelated dofs --
    // that is Phase 3's job.
    //
    // Limitation: only propagates one hop (to direct neighbors of
    // patch_index) -- it does not chain transitively through the neighbor
    // to ITS other neighbors.
    void refineWithPropagation(size_t patch_index, int direction, const Refine1DFn& refine_1d_fn) {
        if (patch_index >= patches_.size())
            throw std::out_of_range("Invalid patch index.");

        auto shared_ids_of = [&](size_t idx) {
            std::unordered_set<size_t> ids;
            for (const auto& kv : shared_control_points_)
                for (size_t p : kv.second)
                    if (p == idx) { ids.insert(kv.first); break; }
            return ids;
        };

        refine_1d_fn(*patches_[patch_index], direction, shared_ids_of(patch_index));

        for (const auto& itf : interfaces_) {
            size_t other_index;
            int other_direction;
            if (itf.patch_a == patch_index && itf.varying_direction_a == direction) {
                other_index = itf.patch_b;
                other_direction = itf.varying_direction_b;
            } else if (itf.patch_b == patch_index && itf.varying_direction_b == direction) {
                other_index = itf.patch_a;
                other_direction = itf.varying_direction_a;
            } else {
                continue;
            }

            refine_1d_fn(*patches_[other_index], other_direction, shared_ids_of(other_index));

            // Merge the boundary CPs grown independently on both sides.
            Patch& PA = *patches_[itf.patch_a];
            Patch& PB = *patches_[itf.patch_b];
            std::vector<size_t> facet_a = extractFacet(PA, itf.direction_a, itf.side_a);
            std::vector<size_t> facet_b = extractFacet(PB, itf.direction_b, itf.side_b);
            if (itf.reversed) std::reverse(facet_b.begin(), facet_b.end());

            if (facet_a.size() != facet_b.size())
                throw std::runtime_error(
                    "Interface mismatch after propagated refinement: edge lengths differ "
                    "between patches " + std::to_string(itf.patch_a) + " and "
                    + std::to_string(itf.patch_b) + ".");

            for (size_t k = 0; k < facet_a.size(); ++k) {
                size_t keep = facet_a[k];
                if (facet_a[k] != facet_b[k]) {
                    size_t drop = facet_b[k];
                    for (size_t& gid : PB.global_indices)
                        if (gid == drop) gid = keep;
                }
                auto& plist = shared_control_points_[keep];
                if (std::find(plist.begin(), plist.end(), itf.patch_a) == plist.end())
                    plist.push_back(itf.patch_a);
                if (std::find(plist.begin(), plist.end(), itf.patch_b) == plist.end())
                    plist.push_back(itf.patch_b);
            }
        }
    }

    // Rebuild every patch's PatchDOFManager from its CURRENT global_indices,
    // growing `global_dof_manager` first so it covers every control point
    // id now in use. Because GlobalDOFManager::get_dof_indices(cp_id) is a
    // pure function of cp_id, two patches that reference the SAME cp_id
    // (e.g. a boundary CP merged by refineWithPropagation()) automatically
    // get assigned the SAME global dof -- no special-casing needed.
    //
    // Call this after refineWithPropagation() (so newly merged boundary CPs
    // get a shared dof) and BEFORE compact() (which renumbers control point
    // ids and would desynchronize global_dof_manager's cp-id-indexed
    // mapping if dofs were assigned against ids that compact() later
    // changes).
    //
    // dofs_per_cp is used only for control points not yet covered by
    // global_dof_manager; existing control points keep whatever dof count
    // they already had.
    void updateDOFManagers(GlobalDOFManager& global_dof_manager, int dofs_per_cp) {
        size_t max_cp_plus_one = global_dof_manager.n_control_points();
        for (const auto& patch : patches_)
            for (size_t gid : patch->global_indices)
                max_cp_plus_one = std::max(max_cp_plus_one, gid + 1);

        global_dof_manager.grow(max_cp_plus_one, dofs_per_cp);

        for (auto& patch : patches_) {
            if (!patch->dof_manager) continue;
            int dpc = patch->dof_manager->dofs_per_control_point;
            patch->dof_manager = std::make_shared<PatchDOFManager>(
                dpc, patch->global_indices, global_dof_manager);
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
            patch1.global_indices.end()
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
        if (patch_index >= transformation_matrices_.size()) {
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
    const std::vector<PatchInterface>& getInterfaces() const {
        return interfaces_;
    }

};