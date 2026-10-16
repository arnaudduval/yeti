#pragma once

#include <algorithm>
#include <array>
#include <functional>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <stdexcept>
#include "Patch.hpp"

// Describes a compatible interface between two patches: a shared edge (2D)
// or a shared face (3D), identified by which (direction, side) facet of
// each patch carries it.
//
// "direction_a"/"side_a": the facet of patch_a is obtained by fixing
// direction_a's local index at 0 (side_a == 0) or at its last index
// (side_a == 1).
//
// "ndim" tags which field group below is meaningful -- added as genuinely
// parallel fields for the 3D case rather than generalizing the 2D ones in
// place, so 2D interfaces (and the 3 existing notebooks that read
// get_interfaces() dicts) are completely unaffected.
//
// ndim == 2 (original, unchanged): the facet is an EDGE with a single
// varying direction (varying_direction_a) -- the direction to refine on
// patch_a to grow the shared edge. "reversed": true if the shared control
// points, walked in increasing varying_direction_a order on patch_a, appear
// in DECREASING varying_direction_b order on patch_b (e.g. u of patch_a
// glued to v of patch_b but parametrized in opposite senses). Only 2
// possible correspondences exist for a 1D edge (same order / reversed).
//
// ndim == 3: the facet is a FACE with TWO varying directions
// (varying_directions_a/b, each (v0, v1) in increasing index order, same
// convention as Patch::boundary_control_points). Two matching faces can
// then correspond in up to 8 ways -- the symmetry group of a rectangle:
// which of patch_b's two varying directions lines up with patch_a's v0 vs
// v1 ("transpose"), and independently whether each matched axis is
// traversed in the same or opposite sense on patch_b (reverse_0/reverse_1,
// relative to patch_a's canonical (v0, v1) order). NOTE: "transpose" here
// describes a correspondence of LOCAL INDICES between the two patches'
// facets -- it does not imply patch_b is physically rotated in space.
struct PatchInterface {
    size_t patch_a, patch_b;
    int direction_a, side_a;
    int direction_b, side_b;

    // ndim == 2 only (meaningless/untouched otherwise -- see ndim below).
    int varying_direction_a = -1;
    int varying_direction_b = -1;
    bool reversed = false;

    // 2 => the fields above are meaningful; 3 => the fields below are.
    int ndim = 2;

    // ndim == 3 only.
    std::array<int, 2> varying_directions_a{-1, -1};
    std::array<int, 2> varying_directions_b{-1, -1};
    bool transpose = false;
    bool reverse_0 = false;
    bool reverse_1 = false;
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
    // 2D: a single varying direction, increasing order (unchanged
    // behavior). 3D: a whole face, v0-fastest over the two varying
    // directions (v0 < v1) -- delegates to Patch::boundary_control_points,
    // already proven N-D generic and equivalent to this function's own
    // former 2D-only arithmetic. boundary_control_points() returns LOCAL
    // flat positions, not global ids -- mapped through global_indices here.
    static std::vector<size_t> extractFacet(const Patch& patch, int direction, int side) {
        std::vector<size_t> local_positions = patch.boundary_control_points(direction, side);
        std::vector<size_t> facet;
        facet.reserve(local_positions.size());
        for (size_t pos : local_positions)
            facet.push_back(patch.global_indices[pos]);
        return facet;
    }

    // The two parametric directions other than `direction`, in increasing
    // index order (v0, v1) -- same convention as
    // Patch::boundary_control_points's internal "varying_dirs". 3D patches
    // only (ndim == 3); meaningless for 2D (only 1 varying direction).
    static std::array<int, 2> facetVaryingDirections(const Patch& patch, int direction) {
        std::array<int, 2> result{-1, -1};
        int idx = 0;
        for (size_t d = 0; d < patch.tensor.components.size(); ++d)
            if (static_cast<int>(d) != direction) result[idx++] = static_cast<int>(d);
        return result;
    }

    // Try to match a 3D face grid `facet_b` (shape n0_b x n1_b, v0-fastest)
    // onto `facet_a` (shape n0_a x n1_a) under one of the 8 symmetries of a
    // rectangle: swap the two axes ("transpose") and/or reverse each axis
    // independently. On success, fills out_transpose/out_reverse_0/
    // out_reverse_1 (relative to facet_a's canonical (v0, v1) axes) and
    // returns true. Returns false if no such correspondence makes the grids
    // coincide position-by-position (not just as sets).
    static bool matchFacetGrids(
        const std::vector<size_t>& facet_a, size_t n0_a, size_t n1_a,
        const std::vector<size_t>& facet_b, size_t n0_b, size_t n1_b,
        bool& out_transpose, bool& out_reverse_0, bool& out_reverse_1)
    {
        for (bool transpose : {false, true}) {
            size_t expect_n0_b = transpose ? n1_a : n0_a;
            size_t expect_n1_b = transpose ? n0_a : n1_a;
            if (n0_b != expect_n0_b || n1_b != expect_n1_b) continue;

            for (bool rev0 : {false, true}) {
                for (bool rev1 : {false, true}) {
                    bool ok = true;
                    for (size_t j = 0; j < n1_a && ok; ++j) {
                        for (size_t i = 0; i < n0_a; ++i) {
                            size_t i_b, j_b;
                            if (!transpose) {
                                i_b = rev0 ? (n0_b - 1 - i) : i;
                                j_b = rev1 ? (n1_b - 1 - j) : j;
                            } else {
                                i_b = rev1 ? (n0_b - 1 - j) : j;
                                j_b = rev0 ? (n1_b - 1 - i) : i;
                            }
                            if (facet_a[i + j * n0_a] != facet_b[i_b + j_b * n0_b]) {
                                ok = false;
                                break;
                            }
                        }
                    }
                    if (ok) {
                        out_transpose = transpose;
                        out_reverse_0 = rev0;
                        out_reverse_1 = rev1;
                        return true;
                    }
                }
            }
        }
        return false;
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
        bool rational = cp_manager->is_rational();

        std::unordered_map<size_t, size_t> old_to_new;
        std::vector<double> new_coords;
        std::vector<double> new_weights;
        new_coords.reserve(cp_manager->coords.size());
        if (rational) new_weights.reserve(cp_manager->weights.size());

        for (auto& patch : patches_) {
            for (size_t& gid : patch->global_indices) {
                size_t old_gid = gid;
                auto it = old_to_new.find(old_gid);
                size_t new_id;
                if (it == old_to_new.end()) {
                    new_id = new_coords.size() / dim;
                    const double* src = cp_manager->coords.data() + old_gid * dim;
                    new_coords.insert(new_coords.end(), src, src + dim);
                    if (rational) new_weights.push_back(cp_manager->weights[old_gid]);
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
            if (rational) cp_manager->weights = std::move(new_weights);
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

        // 3D: a parallel loop, not a generalization of the one above -- a
        // face has two varying directions, so matching needs the 8-way
        // (transpose x 2 independent flips) check in matchFacetGrids()
        // instead of the 2D loop's simple same/reversed-order comparison.
        for (size_t i = 0; i < patches_.size(); ++i) {
            const Patch& A = *patches_[i];
            if (A.tensor.components.size() != 3) continue;

            for (size_t j = i + 1; j < patches_.size(); ++j) {
                const Patch& B = *patches_[j];
                if (B.tensor.components.size() != 3) continue;

                std::unordered_set<size_t> a_ids(A.global_indices.begin(), A.global_indices.end());
                bool any_shared = false;
                for (size_t idx : B.global_indices) {
                    if (a_ids.count(idx)) { any_shared = true; break; }
                }
                if (!any_shared) continue;

                for (int da = 0; da < 3; ++da) {
                    for (int sa = 0; sa < 2; ++sa) {
                        std::vector<size_t> facet_a = extractFacet(A, da, sa);
                        std::array<int, 2> vdirs_a = facetVaryingDirections(A, da);
                        size_t n0_a = A.local_shape[vdirs_a[0]];
                        size_t n1_a = A.local_shape[vdirs_a[1]];
                        std::unordered_set<size_t> set_a(facet_a.begin(), facet_a.end());

                        for (int db = 0; db < 3; ++db) {
                            for (int sb = 0; sb < 2; ++sb) {
                                std::vector<size_t> facet_b = extractFacet(B, db, sb);
                                if (facet_b.size() != facet_a.size()) continue;

                                std::unordered_set<size_t> set_b(facet_b.begin(), facet_b.end());
                                if (set_a != set_b) continue;

                                std::array<int, 2> vdirs_b = facetVaryingDirections(B, db);
                                size_t n0_b = B.local_shape[vdirs_b[0]];
                                size_t n1_b = B.local_shape[vdirs_b[1]];

                                bool transpose, rev0, rev1;
                                if (!matchFacetGrids(facet_a, n0_a, n1_a, facet_b, n0_b, n1_b,
                                                      transpose, rev0, rev1)) {
                                    throw std::runtime_error(
                                        "Incompatible 3D boundary ordering between patches "
                                        + std::to_string(i) + " and " + std::to_string(j) +
                                        ": shared control points' face grids don't match under "
                                        "any axis permutation/reversal.");
                                }

                                PatchInterface interface;
                                interface.patch_a = i;
                                interface.direction_a = da;
                                interface.side_a = sa;
                                interface.patch_b = j;
                                interface.direction_b = db;
                                interface.side_b = sb;
                                interface.ndim = 3;
                                interface.varying_directions_a = vdirs_a;
                                interface.varying_directions_b = vdirs_b;
                                interface.transpose = transpose;
                                interface.reverse_0 = rev0;
                                interface.reverse_1 = rev1;
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
            if (itf.ndim == 2) {
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
            } else if (itf.ndim == 3) {
                // A 3D face has TWO varying directions: find which one (if
                // any) matches `direction` on whichever side patch_index is
                // on, then pick the corresponding direction on the other
                // side via the stored transpose flag (see PatchInterface).
                bool is_a_side;
                int k = -1;
                if (itf.patch_a == patch_index) {
                    is_a_side = true;
                    for (int kk = 0; kk < 2; ++kk)
                        if (itf.varying_directions_a[kk] == direction) { k = kk; break; }
                } else if (itf.patch_b == patch_index) {
                    is_a_side = false;
                    for (int kk = 0; kk < 2; ++kk)
                        if (itf.varying_directions_b[kk] == direction) { k = kk; break; }
                } else {
                    continue;
                }
                if (k < 0) continue;  // direction isn't one of this face's varying directions

                size_t other_index;
                int other_direction;
                if (is_a_side) {
                    other_index = itf.patch_b;
                    other_direction = itf.varying_directions_b[itf.transpose ? (1 - k) : k];
                } else {
                    other_index = itf.patch_a;
                    other_direction = itf.varying_directions_a[itf.transpose ? (1 - k) : k];
                }

                refine_1d_fn(*patches_[other_index], other_direction, shared_ids_of(other_index));

                // Merge the two independently-grown 2D grids, position by
                // position, using the transform detectInterfaces() already
                // resolved (transpose/reverse_0/reverse_1) to align B's grid
                // onto A's -- same precondition as the 2D merge above
                // (deterministic, coordinate-identical new CPs on both
                // sides), just over 2 indices instead of 1.
                Patch& PA = *patches_[itf.patch_a];
                Patch& PB = *patches_[itf.patch_b];
                std::vector<size_t> facet_a = extractFacet(PA, itf.direction_a, itf.side_a);
                std::vector<size_t> facet_b = extractFacet(PB, itf.direction_b, itf.side_b);
                size_t n0_a = PA.local_shape[itf.varying_directions_a[0]];
                size_t n1_a = PA.local_shape[itf.varying_directions_a[1]];
                size_t n0_b = PB.local_shape[itf.varying_directions_b[0]];
                size_t n1_b = PB.local_shape[itf.varying_directions_b[1]];

                if (facet_a.size() != facet_b.size())
                    throw std::runtime_error(
                        "Interface mismatch after propagated refinement: face sizes differ "
                        "between patches " + std::to_string(itf.patch_a) + " and "
                        + std::to_string(itf.patch_b) + ".");

                for (size_t j = 0; j < n1_a; ++j) {
                    for (size_t i = 0; i < n0_a; ++i) {
                        size_t i_b, j_b;
                        if (!itf.transpose) {
                            i_b = itf.reverse_0 ? (n0_b - 1 - i) : i;
                            j_b = itf.reverse_1 ? (n1_b - 1 - j) : j;
                        } else {
                            i_b = itf.reverse_1 ? (n0_b - 1 - j) : j;
                            j_b = itf.reverse_0 ? (n1_b - 1 - i) : i;
                        }
                        size_t keep = facet_a[i + j * n0_a];
                        size_t drop = facet_b[i_b + j_b * n0_b];
                        if (keep != drop) {
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