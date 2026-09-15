"""
Multi-hop refinement propagation across a PatchAssembly.

`PatchAssembly.refine_with_propagation(patch_index, direction, refine_1d_fn)` refines one
patch along `direction` and propagates to its DIRECT neighbors only (one hop): any interface
whose *varying direction* (the direction the shared edge itself runs along) matches
`direction` gets its neighbor refined too, using whichever direction the C++ side determines
is the matching one on that neighbor (crossed interfaces are handled automatically).

`refine_from()` below repeatedly re-triggers propagation from each newly-touched patch --
using the direction that was actually applied to it, captured via a guarded callback -- until
no new patch is reached. This covers topologies where the refinement must travel further than
one hop, e.g. all the way around a loop of quarter-patches.

The kind of refinement applied at each hop -- uniform h-refinement (`SubdivisionRefiner`),
degree elevation (`PRefiner`), or a single targeted knot insertion (`HRefiner`) -- is
pluggable via the `refiner_factory` argument, since all three expose the same
`refine_1d(patch, protected_global_ids)` interface.
"""

from __future__ import annotations

from .bspline import SubdivisionRefiner


def refine_from(assembly, patches, start_index, start_direction, refiner_factory=None):
    """Refine one patch along `start_direction`, propagating through every
    interface reachable via matching varying directions (crossed interfaces
    included), no matter how many hops away.

    Parameters
    ----------
    assembly : PatchAssembly
        Assembly the patches belong to. `detect_shared_control_points()` and
        `detect_interfaces()` are called on it before refining.
    patches : list of Patch
        Patches in `assembly.add_patch()` order -- used to map a `Patch`
        object back to its integer index for the returned dict.
    start_index : int
        Index (into `patches`) of the patch to refine first.
    start_direction : int
        Parametric direction (0 or 1) to refine `start_index` along.
    refiner_factory : callable(direction) -> refiner, optional
        Builds the refiner to apply at each patch touched during
        propagation, given the direction determined for that particular
        patch (which may differ from `start_direction` across a crossed
        interface). The returned object must expose
        `refine_1d(patch, protected_global_ids)` -- e.g.
        `SubdivisionRefiner(direction=direction, n_levels=...)` (uniform
        h-refinement), `PRefiner(direction=direction, n_elevations=...)`
        (degree elevation), or `HRefiner(direction=direction, knot=...)` (a
        single targeted knot insertion). Defaults to
        `SubdivisionRefiner(direction=direction, n_levels=1)`.

    Returns
    -------
    dict[int, int]
        Maps every patch index that was refined to the direction (0 or 1)
        actually applied to it.
    """
    if refiner_factory is None:
        refiner_factory = lambda direction: SubdivisionRefiner(direction=direction, n_levels=1)

    assembly.detect_shared_control_points()
    assembly.detect_interfaces()

    # pybind11 keeps one Python wrapper per shared_ptr<Patch>, so id() is stable.
    patch_id_to_index = {id(p): i for i, p in enumerate(patches)}
    refined_direction = {}

    def refine_1d_fn(patch, direction, protected_global_ids):
        i = patch_id_to_index[id(patch)]
        if i in refined_direction:
            return  # already refined -- but let the merge step proceed
        refined_direction[i] = direction
        refiner_factory(direction).refine_1d(patch, protected_global_ids)

    assembly.refine_with_propagation(
        patch_index=start_index, direction=start_direction, refine_1d_fn=refine_1d_fn)

    expanded = {start_index}
    frontier = [i for i in refined_direction if i != start_index]
    while frontier:
        next_frontier = []
        for idx in frontier:
            if idx in expanded:
                continue
            expanded.add(idx)
            before = set(refined_direction)
            assembly.refine_with_propagation(
                patch_index=idx, direction=refined_direction[idx], refine_1d_fn=refine_1d_fn)
            next_frontier.extend(set(refined_direction) - before)
        frontier = next_frontier

    assembly.detect_shared_control_points()
    return refined_direction
