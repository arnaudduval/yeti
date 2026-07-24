Design choices
================

.. contents:: Table of Contents
    :depth: 2
    :local:
    :backlinks: none

This page explains *why* the :mod:`yeti_iga.future.bspline` module is built the way it
is, with a focus on the multipatch assembly machinery (:class:`~yeti_iga.future.bspline.PatchAssembly`
and the refinement operators). For the precise signature of every class and method, see
:doc:`api`.

Overview
----------

``future`` (:file:`src/yeti_iga/future/`) is a C++17 reimplementation of YETI's core
B-spline/NURBS machinery, exposed to Python through pybind11 as
:mod:`yeti_iga.future.bspline`. It is developed alongside, not as a replacement (yet)
for, the legacy Fortran/f2py layer documented under *API reference*. It uses Eigen for
linear algebra and OpenMP for parallelism, and is built independently of the legacy
layer (see :file:`src/yeti_iga/future/CMakeLists.txt`).

The rest of this page focuses on a problem that doesn't exist in a single-patch
B-spline library: several :class:`~yeti_iga.future.bspline.Patch` objects need to share
some of their control points (a glued/conforming multipatch interface), and every
operation on one patch — refinement, DOF numbering, memory layout — has to remain
correct for its neighbors without either patch knowing about the other's internals.
Everything below is a consequence of that one constraint.

The diagram below summarizes the core classes involved and how they relate; the
sections that follow explain the *why* behind each relationship.

.. mermaid::

    classDiagram
        class BSplineTensor {
            +components : BSpline[]
            find_span_nd()
            basis_funs_nd()
        }
        class BSplineSurface
        class BSplineVolume
        BSplineTensor <|-- BSplineSurface
        BSplineTensor <|-- BSplineVolume

        class ControlPointManager {
            +coords
            +weights
            +is_rational : bool
            add_point(w=1.0)
            coords_view()
            weights_view()
        }

        class Patch {
            +global_indices
            +local_shape
            spans()
        }
        Patch *-- "1" BSplineTensor : tensor
        Patch --> "1" ControlPointManager : cp_manager (shared pool)
        Patch --> "0..1" PatchDOFManager : dof_manager

        class GlobalDOFManager {
            get_dof_indices(cp_id)
            grow()
        }
        class PatchDOFManager {
            get_global_dof_indices(local_pos)
        }
        PatchDOFManager ..> GlobalDOFManager : built from

        class RefinementOperator {
            <<abstract>>
            refine()
            get_type()
        }
        class HRefiner
        class SubdivisionRefiner
        class PRefiner
        RefinementOperator <|-- HRefiner
        RefinementOperator <|-- SubdivisionRefiner
        RefinementOperator <|-- PRefiner
        RefinementOperator ..> Patch : refines in-place

        class PatchAssembly {
            add_patch()
            detect_shared_control_points()
            detect_interfaces()
            refine_with_propagation()
            update_dof_managers()
            compact()
        }
        PatchAssembly o-- "*" Patch : patches_
        PatchAssembly ..> GlobalDOFManager : update_dof_managers()

Three relationships are worth noting up front, since they are easy to miss on a first
read of the diagram: several ``Patch`` instances point ``-->`` (reference, not own) the
*same* ``ControlPointManager`` — that shared pointer is the entire mechanism behind
control point sharing; ``PatchDOFManager`` is built ``from`` a ``GlobalDOFManager`` but
keeps no further link to it, which is precisely what lets it survive
``PatchAssembly.compact()``; and ``RefinementOperator`` depends on ``Patch`` (it mutates
one) but a ``Patch`` has no reverse dependency on any refiner.

Shared control points: one pool, integer ids
-----------------------------------------------

A :class:`~yeti_iga.future.bspline.Patch` does not own its control point coordinates.
Coordinates live in a single :class:`~yeti_iga.future.bspline.ControlPointManager`
instance — a contiguous ``[x0, y0, z0, x1, y1, z1, ...]`` buffer — and a patch only
holds a ``global_indices`` array mapping its *local* control point positions (in
u-fastest order) to *ids* in that shared pool.

Two patches that are meant to share a boundary are built referencing the **same**
``ControlPointManager`` instance and the **same** ids for their common control points.
Sharing is therefore established once, at construction time, by id equality — there is
no separate "constraint" or "glue" object to keep in sync afterwards.
:meth:`PatchAssembly.detect_shared_control_points() <yeti_iga.future.bspline.PatchAssembly.detect_shared_control_points>`
simply looks for ids that appear in more than one patch's ``global_indices``.

This is also why :meth:`ControlPointManager.coords_view() <yeti_iga.future.bspline.ControlPointManager.coords_view>`
returns a zero-copy NumPy view rather than a fresh array: the buffer is the single
source of truth for every patch's geometry, and a small ``mutex`` on the manager keeps
concurrent appends safe.

NURBS: rational bases with zero B-spline overhead
----------------------------------------------------

B-splines are the default — all basis functions are evaluated as plain polynomials and no
weight-related computation occurs. NURBS are activated by passing ``w != 1.0`` to
:meth:`ControlPointManager.add_point() <yeti_iga.future.bspline.ControlPointManager.add_point>`.
Once activated, ``is_rational`` returns ``True`` and the rational basis
:math:`R_a = w_a N_a / W` (where :math:`W = \sum_b w_b N_b`) replaces the polynomial
basis in every integration and evaluation routine.

The key constraint is that **B-spline patches pay zero additional cost** — no extra
memory, no division, no branch inside the Gauss loop. This is enforced by a
``if constexpr`` dispatch at the level of the span-collection methods
(``collectTripletsImpl``, ``collectMassTripletsImpl``, ``collectOperatorTripletsImpl``,
``computeLocalBoundaryLoadContribution``). Each checks
``patch.cp_manager->is_rational()`` **once**, outside every loop, then instantiates one
of two fully separate template paths:

- ``..Impl<false>`` (B-spline) — identical to the pre-NURBS code, with every NURBS
  branch compiled away by the optimizer.
- ``..Impl<true>`` (NURBS) — fetches per-span weights via ``Patch::weights_for_span()``
  and applies the quotient-rule rationalisation inside the Gauss loop.

:class:`~yeti_iga.future.bspline.PatchEvaluator`'s evaluation loop follows the same
split.

``ControlPointManager`` activates rational mode lazily: the first ``add_point(...,
w=...)`` call with ``w != 1.0`` allocates the weights vector and backfills 1.0 for
every previously added point. Patches that never use a non-unit weight keep ``weights``
empty, ``is_rational()`` returning ``False``, and the B-spline path stays active — zero
overhead, by construction.

Refinement of shared patches: protected vs. private control points
-----------------------------------------------------------------------

The corruption problem
~~~~~~~~~~~~~~~~~~~~~~~~~

Refining a single, unshared patch is simple: recompute every control point from the
knot-insertion or degree-elevation formula and overwrite the pool with the new, dense
set of ids ``0..n_new-1``. That is unsafe the moment the pool is shared: a neighboring
patch's ``global_indices`` still point at the *old* ids, which the naive algorithm has
just renumbered and overwritten out from under it.

Protected vs. private control points
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Every refinement operator (:class:`~yeti_iga.future.bspline.HRefiner`,
:class:`~yeti_iga.future.bspline.SubdivisionRefiner`,
:class:`~yeti_iga.future.bspline.PRefiner`) therefore accepts a ``protected_global_ids``
set (see :meth:`~yeti_iga.future.bspline.HRefiner.refine_1d`,
:meth:`~yeti_iga.future.bspline.SubdivisionRefiner.refine_1d`,
:meth:`~yeti_iga.future.bspline.PRefiner.refine_1d`) and classifies every control point of
the refined patch into one of two categories:

- **Protected** — an id *borrowed* from another patch (a shared interface control
  point). It is never recomputed, renumbered, or written to: it keeps exactly the id
  and coordinates it already has in the pool.
- **Private** — everything else, i.e. control points this patch owns exclusively,
  whether their position is geometrically unchanged by the refinement or genuinely new
  (blended). Private control points are always (re)assigned to a **fresh, dense block**
  of new ids via ``cp_manager.add_point()``, even when their coordinates did not
  actually change.

Concretely, in ``HRefiner::apply_1d_cp_update`` (the shared low-level C++ routine behind
all three operators' ``refine_1d`` fast path), a row of the 1D transition matrix that is
a pure copy of one old control point is checked against ``protected_global_ids``: if the
old id is protected, the new position simply re-uses it unchanged; otherwise it is
queued for (re)numbering as a private control point.

Why refinement never deletes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Always moving private control points to a fresh block — instead of reusing their old
slot in place — means repeated refinement of a shared patch leaves old private control
point slots behind in the pool, unused but still present ("orphaned"). This is a
deliberate trade-off: at the moment one patch is refined, another patch (or another
thread) might still be reading the pool through its own, still-valid ``global_indices``.
Refinement therefore never deletes or reuses a slot in place when the pool might be
shared; cleanup is a separate, explicit step — see
:meth:`PatchAssembly.compact() <yeti_iga.future.bspline.PatchAssembly.compact>` below.

The exclusive-ownership fast path
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When ``protected_global_ids`` is empty (the default), the patch is assumed to
exclusively own its range of the pool: there is nothing to protect, so the
implementation falls back to the cheaper, original single-patch behavior — overwrite
the pool in place with sequential ids ``0..nb_new_cp-1``, with zero extra memory
overhead. Multipatch awareness therefore costs nothing for patches that don't actually
share anything.

Two-level DOF management
---------------------------

Degrees of freedom are split across two cooperating classes
(:file:`src/yeti_iga/future/include/DOFManager.hpp`), each indexed differently on
purpose:

- :class:`~yeti_iga.future.bspline.GlobalDOFManager` maps a **control point id** (in the
  shared pool) to its global dof indices. Because this lookup is a pure function of the
  id, two patches that reference the *same* id (e.g. a boundary control point merged
  across a shared interface) automatically resolve to the *same* global dofs — with no
  explicit bookkeeping of which control points were merged. This is what makes
  :meth:`PatchAssembly.update_dof_managers() <yeti_iga.future.bspline.PatchAssembly.update_dof_managers>`
  correct without ever inspecting the merge history.
- :class:`~yeti_iga.future.bspline.PatchDOFManager` maps a patch's **local control
  point position** (not its id) to global dof indices. Being position-indexed rather
  than id-indexed is exactly what lets a ``PatchDOFManager`` survive
  :meth:`PatchAssembly.compact() <yeti_iga.future.bspline.PatchAssembly.compact>`
  unchanged, even though ``compact()`` renumbers every control point id: the patch's
  *positions* don't move, only the ids they point at.

This split is the enabling mechanism behind the assembly's call-order contract below: if
``PatchDOFManager`` were id-indexed like ``GlobalDOFManager``, ``compact()`` would
silently desynchronize every patch's dof mapping the moment it renumbers ids.

Stale DOF manager pitfall
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Every refinement operator updates ``patch.dof_manager`` **in-place** after refining a
patch, and :meth:`~yeti_iga.future.bspline.PatchAssembly.update_dof_managers` rebuilds
every patch's ``PatchDOFManager`` in-place. Any Python variable that captured the old
``PatchDOFManager`` before either of these calls becomes **stale**: it still holds the
pre-refinement DOF mapping and will silently return wrong indices if used afterwards.
Always read ``patch.dof_manager`` **after** refining or calling ``update_dof_managers``,
never cache it beforehand.

.. code-block:: python

   # Wrong: pdm is captured before refinement, becomes stale.
   pdm = patch.dof_manager
   HRefiner(0, 0.5).refine(patch, ...)
   pdm.get_global_dof_indices(...)     # stale — returns pre-refinement indices

   # Correct: read dof_manager from the patch after refining.
   HRefiner(0, 0.5).refine(patch, ...)
   patch.dof_manager.get_global_dof_indices(...)  # always up to date

PatchAssembly's call-order contract
---------------------------------------

:class:`~yeti_iga.future.bspline.PatchAssembly` orchestrates every multipatch operation,
and its methods are designed to be called in a specific order, each step depending on
the result or side effect of the one before it:

1. :meth:`~yeti_iga.future.bspline.PatchAssembly.add_patch` for every patch in the
   assembly.
2. :meth:`~yeti_iga.future.bspline.PatchAssembly.detect_shared_control_points` — finds
   ids common to several patches. Requires patches to already share ids in their
   ``global_indices`` by construction (see above) — it does not infer sharing from
   coordinates.
3. :meth:`~yeti_iga.future.bspline.PatchAssembly.detect_interfaces` — uses the shared-id
   information to identify, for each pair of compatible 2D patches, which
   ``(direction, side)`` facet carries the shared edge on each side and the orientation
   between them (see :ref:`design-crossed-interface` below). It needs step 2's result
   to know which patches are even candidates for an interface.
4. :meth:`~yeti_iga.future.bspline.PatchAssembly.refine_with_propagation` — refines one
   patch and propagates to its neighbors across the interfaces found in step 3, merging
   the new boundary control points. It needs step 3's interfaces to know *which*
   neighbors to propagate to and along *which* direction.
5. :meth:`~yeti_iga.future.bspline.PatchAssembly.update_dof_managers` — grows the global
   DOF pool and rebuilds every patch's ``PatchDOFManager``. It must run **after** step 4
   so that the newly merged boundary control points get assigned a *shared* dof, and
   **before** ``compact()`` (step 6): since ``GlobalDOFManager`` is indexed by control
   point id, assigning dofs against ids that ``compact()`` is about to renumber would
   silently desynchronize the mapping.
6. :meth:`~yeti_iga.future.bspline.PatchAssembly.compact` — once *all* the refinements
   you need are done (not after every single one), reclaims the control points
   orphaned by repeated refinements of shared patches (see above): it walks the patches
   in ``add_patch()`` order, keeps each control point's first-encountered id, and
   reassigns it a fresh, dense id, reusing that same new id for every later occurrence
   of the same old id (i.e. for control points shared with an already-processed patch).
   It leaves every ``PatchDOFManager`` untouched (position-indexed, see above) but
   invalidates the shared-control-point map and the detected interfaces — call steps 2
   and 3 again afterwards if you still need them.

.. _design-crossed-interface:

Interface detection and the crossed-interface case
-------------------------------------------------------

:meth:`~yeti_iga.future.bspline.PatchAssembly.detect_interfaces` only handles 2D
patches (a 3D volume face interface has two varying directions and isn't supported
yet — this is an explicit "Phase 1" scope limitation, not an oversight). For each pair
of patches, it tries every ``(direction, side)`` facet combination on both sides and
matches the ones whose control point *sets* are equal; it then checks whether the two
facets, walked in order, agree (direct order) or are reversed (e.g. ``u`` of one patch
glued to ``v`` of the other, parametrized in opposite senses). It raises if the sets
match but neither ordering does, since that would violate the "compatible boundary"
assumption (same knot vector, degree, and control points along the edge).

This *crossed interface* case (the varying direction differs between the two sides) is
why :meth:`~yeti_iga.future.bspline.PatchAssembly.refine_with_propagation` takes a
``refine_1d_fn`` callback that receives the refinement ``direction`` as an **explicit
argument**, rather than the direction being baked into a closure created once for the
triggering patch: the neighbor generally has to be refined along its *own* matching
direction, which is not necessarily the same value as the direction given for the patch
that triggered the propagation.

One-hop propagation limitation
----------------------------------

``refine_with_propagation`` only propagates to the *direct* neighbors of the patch being
refined — it does not chain transitively through a neighbor to that neighbor's other
neighbors. For an assembly with chains of more than two patches sharing edges, this
means propagation currently has to be triggered explicitly, patch by patch, along the
chain. This is a documented scope limitation rather than a bug, left for a later phase.

Carrying a prior solution across refinement
-------------------------------------------------

When progressively refining a problem (e.g. adaptive h-refinement), it is often useful to
**warm-start** the refined solve from the solution on the coarser mesh. The transfer
formula is :math:`u_\text{new} = T \, u_\text{old}`, where :math:`T` is the nD
control-point transition matrix produced by the refinement step.

:class:`~yeti_iga.future.bspline.PatchAssembly` stores one nD transition matrix per
patch, initialized to the identity when the patch is added. The intended workflow is:

1. Instead of ``refine_with_propagation``, call ``refine_1d()`` directly — it returns
   the 1D transition matrix ``T_1d``.
2. Expand it to nD:
   ``T_nd = nd_transition_from_1d(T_1d, direction, shape_before, shape_after)``.
3. Store it in the assembly:
   :meth:`assembly.set_transformation_matrix(patch_index, T_nd) <yeti_iga.future.bspline.PatchAssembly.set_transformation_matrix>`.
4. After solving on the fine mesh — or to construct an initial guess — apply it to the
   flat DOF vector:
   :meth:`assembly.apply_transformation_to_dofs(u_old, patch_index, dim_phys) <yeti_iga.future.bspline.PatchAssembly.apply_transformation_to_dofs>`.
   Use :meth:`~yeti_iga.future.bspline.PatchAssembly.apply_transformation_to_control_points`
   instead when working directly with per-control-point geometry vectors.

The reason ``refine_with_propagation`` does not expose the transition matrices directly
is that its ``refine_1d_fn`` callback has no return value by design (see
:ref:`design-crossed-interface`): the direction it passes to the callback is
patch-specific and may differ from the direction that triggered the propagation. Having
the callback return a matrix would force the user to track which matrix belongs to which
patch, which is exactly what ``PatchAssembly``'s per-patch matrix store already does.
Use ``refine_1d`` directly whenever warm-starting is needed.

Performance rationale: the 1D fast path
-------------------------------------------

A naive ND refinement recomputes the full ``(nb_new_cp × nb_old_cp)`` transition matrix
directly. Instead, every refinement operator exposes a ``refine_1d`` fast path that only
builds the **1D** transition matrix for the refined direction — much smaller than the
full ND one — and the free function
:func:`~yeti_iga.future.bspline.nd_transition_from_1d` reconstructs the full ND matrix
from it afterwards (via a Kronecker product against identities for the unaffected
directions), only if and when that full matrix is actually needed.

This "precompute the small/cheap thing once, derive the large/expensive thing from it
only on demand" pattern shows up again on the integration side:
:class:`~yeti_iga.future.bspline.IGABasis1D` precomputes Gauss point coordinates,
weights, and basis function values/derivatives once per 1D parametric direction, and
:class:`~yeti_iga.future.bspline.PatchIntegrator` reuses that precomputed data across
every span of a 2D patch rather than re-evaluating basis functions per assembly call.

Skipping degenerate spans
~~~~~~~~~~~~~~~~~~~~~~~~~~~

:class:`~yeti_iga.future.bspline.SpanIterator` filters out zero-measure knot spans
(those where ``kv[i+1] - kv[i] == 0``) at construction time — before the integration
loop ever runs. This matters for NURBS patches with C^0 knot repetitions and after
degree elevation, both of which introduce internal zero-length spans that contribute
nothing to any integral but can cause division by zero when computing the parametric
Jacobian. By handling this once, in the iterator, every caller — stiffness, mass,
boundary load, ``LocalOperator``, ``ScalarLocalOperator`` — benefits automatically
with no per-integration check.

A generic integration term: ``LocalOperator``
--------------------------------------------------

Two built-in kernels, one shared shape
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Before ``LocalOperator``, :class:`~yeti_iga.future.bspline.PatchIntegrator` already had
two ways to integrate, and they are structurally near-identical:
``computeLocalStiffnessContribution()`` (``B^T*D*B``, plane-stress stiffness) and
``computeLocalMassContribution()`` (``rho*N^T*N``, consistent mass) are both private,
per-span C++ methods sharing the exact same signature (``Patch``, the span's
``SpanGauss1D`` in each direction, the span itself). Both loop the same Gauss points;
both call the shared ``evaluateGaussPointGeometry()`` helper to get the basis values,
gradients, and Jacobian/``detJ`` at each point; both hand their resulting per-span
matrix to the same ``assembleLocalContribution()`` to scatter into the global triplet
list. The *only* thing that differs between them is the algebra applied to that shared
geometry — stiffness contracts the inverse-Jacobian-transformed gradients through a
``B`` matrix and the plane-stress ``D`` matrix, mass just outer-products the basis
values weighted by ``rho``.

The orchestration above them mirrors this symmetry: ``collectTriplets()`` and
``collectMassTriplets()`` are identical apart from which ``computeLocal*Contribution``
they call, and
:meth:`assemble_stiffness() <yeti_iga.future.bspline.PatchIntegrator.assemble_stiffness>`/
:meth:`assemble_mass() <yeti_iga.future.bspline.PatchIntegrator.assemble_mass>` are both
thin wrappers around the same private ``assembleGeneric()`` helper, parametrized by
which ``collect*`` method to invoke per patch.

Why a third way to integrate
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Both of the kernels above are fixed C++ implementations: adding a third physical term,
or trying a variant of an existing one, means writing a new method with the same shape
and recompiling the extension.
:class:`~yeti_iga.future.bspline.LocalOperator` is a third way in: subclass it
**from Python** and override
:meth:`compute_integrand() <yeti_iga.future.bspline.LocalOperator.compute_integrand>`
to plug in a different physical term, without touching C++ or rebuilding the
extension.

It exists purely for development/testing convenience — prototyping a new term, or
reproducing a textbook formula to sanity-check the built-in kernels — not for
performance. Every Gauss point the operator integrates triggers one Python call, so
:meth:`integrate_stiffness() <yeti_iga.future.bspline.PatchIntegrator.integrate_stiffness>`/
:meth:`integrate_mass() <yeti_iga.future.bspline.PatchIntegrator.integrate_mass>`
remain the methods to use once a kernel is settled, or to bake a kernel validated this
way into C++.

Per-Gauss-point granularity, not per-span
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Two granularities were considered for the boundary between code that stays fixed (in
C++) and code that is user-supplied (in Python):

- **Per-span** — the operator receives the precomputed basis/Gauss data for a whole
  span and loops its own Gauss points, the same shape
  ``computeLocalStiffnessContribution()``/``computeLocalMassContribution()`` above
  already have.
- **Per-Gauss-point** — ``PatchIntegrator`` does the Gauss-point loop *and* the
  Jacobian/gradient computation, and the operator only supplies the algebraic term to
  integrate at one point (e.g. ``B^T*D*B`` for stiffness).

Per-Gauss-point was chosen, even though it triggers more Python calls: the point of
``LocalOperator`` is letting a user write down a textbook formula, not re-derive the
Jacobian inversion or the basis-function tensor product every time. The Jacobian/
gradient code is exactly the part most worth *not* duplicating in Python — it is also
the part the built-in stiffness kernel got wrong twice during its own development (two
separate GLOBAL-vs-LOCAL index bugs, caught by the dedicated multipatch test suite
before they shipped), a good argument for keeping that logic in one, already-validated
place rather than asking every operator author to get it right again.

Sharing C++ machinery without risking the existing kernels
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``collectOperatorTriplets()`` reuses the same ``evaluateGaussPointGeometry()`` helper,
and the same ``assembleLocalContribution()`` triplet-scattering routine, that
``collectTriplets()``/``collectMassTriplets()`` use — so a ``LocalOperator``-based
assembly goes through exactly the same Jacobian and dof-mapping code the built-in
kernels already exercise.

The new methods (``collectOperatorTriplets``, ``integrate_operator``,
``assemble_operator``) were added purely additively: not a single line of
``computeLocalStiffnessContribution``, ``computeLocalMassContribution``, or the private
``assembleGeneric`` helper was modified.
:meth:`assemble_operator() <yeti_iga.future.bspline.PatchIntegrator.assemble_operator>`
in particular is a deliberately standalone implementation — it does **not** route
through ``assembleGeneric()``, the helper
:meth:`assemble_stiffness() <yeti_iga.future.bspline.PatchIntegrator.assemble_stiffness>`/
:meth:`assemble_mass() <yeti_iga.future.bspline.PatchIntegrator.assemble_mass>` share,
even though that duplicates a small amount of per-patch orchestration. The duplication
buys a guarantee: this newer, less battle-tested path can never affect the existing
ones, by construction rather than by careful review.

Subclassing a C++ class from Python: the pybind11 trampoline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Making ``LocalOperator`` overridable from a pure-Python subclass needed a pybind11
*trampoline* (``PyLocalOperator`` in ``bindings.cpp``) — a small C++ shim that forwards
the virtual call to whatever Python method is found on the instance, falling back to a
pure-virtual error if none is. This is the first use of this pattern in ``future``: the
existing :class:`~yeti_iga.future.bspline.RefinementOperator` hierarchy
(:class:`~yeti_iga.future.bspline.HRefiner`,
:class:`~yeti_iga.future.bspline.SubdivisionRefiner`,
:class:`~yeti_iga.future.bspline.PRefiner`) only goes the other way: C++ subclasses
exposed to Python, never a Python subclass of a C++ base.

One pitfall worth recording: ``PYBIND11_OVERRIDE_PURE`` looks up the Python override by
the *literal C++ method name* (``computeIntegrand``), not by whatever name the binding
exposes it under (``compute_integrand``, to match this module's snake_case Python
convention). Since the two differ here, the macro needs its explicit-name variant,
``PYBIND11_OVERRIDE_PURE_NAME(..., "compute_integrand", computeIntegrand, ...)`` — the
plain macro silently looks up the wrong attribute and always falls through to the
pure-virtual error, even though registration, inheritance, and the binding itself are
all otherwise correct.

What a per-Gauss-point term can and cannot express
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Because the operator only ever sees one Gauss point at a time, it cannot directly
express an operation that depends on a span's *whole* integral, such as literal
row-sum ("HRZ") mass lumping, which sums a span's full consistent mass matrix over its
columns *after* integrating it. It turns out not to matter for this example: row-summing
``outer(R, R)`` at a single Gauss point already gives ``R`` itself, since the span's
active basis functions sum to 1 there (partition of unity) — so lumping reduces to a
per-Gauss-point diagonal term after all, with no need to wait for the full span integral
before lumping. See :file:`examples/future/07_local_operator.ipynb` for the worked
example, with both this lumped-mass term and a from-Python reproduction of the built-in
stiffness kernel.

Boundary conditions and loads: bricks, not a solver
--------------------------------------------------------

A full linear elasticity computation needs three more things on top of ``K``:
Dirichlet (displacement) boundary conditions, a Neumann (distributed load)
right-hand side, and a way to solve the reduced system. ``future`` provides the first
two — :meth:`Patch.boundary_control_points() <yeti_iga.future.bspline.Patch.boundary_control_points>`
selects control points on an edge or a span sub-range of it (for Dirichlet), and
:meth:`PatchIntegrator.integrate_boundary_load() <yeti_iga.future.bspline.PatchIntegrator.integrate_boundary_load>`/
:meth:`assemble_boundary_load() <yeti_iga.future.bspline.PatchIntegrator.assemble_boundary_load>`
integrate a :class:`~yeti_iga.future.bspline.Traction` over an edge into a load vector
the same size as ``K`` — but deliberately stops there: eliminating the fixed dofs and
solving the reduced system is a handful of ``scipy`` lines once the dof list and ``F``
exist (see :file:`examples/future/08_boundary_conditions.ipynb`), not new code here.
This is the same "matrices in, matrices out" scope ``future`` has kept since the first
stiffness/mass kernels — a solver is a separate, later concern.

:class:`~yeti_iga.future.bspline.BoundaryLoadSpec` is a plain data struct that packages
a single load specification — patch index, direction, side, traction object, and optional
``span_min``/``span_max`` for sub-range application — into a value the batch API
:meth:`~yeti_iga.future.bspline.PatchIntegrator.assemble_boundary_load` can iterate
over. It carries no logic; all design decisions sit in ``Traction`` and the integration
loop.

``Traction`` mirrors ``LocalOperator``'s extensibility approach for the same reason:
:meth:`Traction.evaluate() <yeti_iga.future.bspline.Traction.evaluate>` already takes
the physical point, not just internal parameters, so a future Python-callback-based
load (for a non-uniform traction) can be added as a new subclass — with a pybind11
trampoline, exactly like :class:`~yeti_iga.future.bspline.LocalOperator`'s — without
changing a single line of the boundary-load integration loop.
:class:`~yeti_iga.future.bspline.ConstantTraction` is the only kernel implemented so
far; the abstraction exists ahead of that need, the Python-callback subclass does not
(yet).

Constitutive law architecture: B-free formulation
----------------------------------------------------

Material data vs. mechanical behaviour
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The integration infrastructure needs two independent pieces of information: the physical
constants that characterise a material (Young's modulus, Poisson's ratio, density…) and
the mechanical formulation that converts those constants into a stiffness contribution at
a Gauss point. These are separated into two distinct objects:

- :class:`~yeti_iga.future.bspline.Material` — a plain data struct. It stores ``E``,
  ``nu``, ``rho`` (default 0), and ``thickness`` (default 1), and computes derived
  elastic moduli (``mu()``, ``lambda_3d()``, ``bulk_modulus()``). It is *not*
  constructible through the abstract law hierarchy and carries no formulation knowledge.
- :class:`~yeti_iga.future.bspline.ConstitutiveLaw` — an abstract strategy class. It
  holds a reference to a ``Material`` (accessible via ``material()``) and declares two
  pure-virtual methods: ``n_dofs_per_cp() -> int`` and
  ``stiffness_density(grad_a, grad_b, x_phys) -> matrix``. Concrete subclasses
  implement those two methods for a specific mechanical formulation.

This separation lets the same ``Material`` instance be reused across different
formulations without coupling them: ``PlaneStress(mat)`` and ``PlaneStrain(mat)`` both
accept the same ``mat``; future additions (``Solid3D``, ``Axisymmetric``, J2 plastic)
will do the same. It also opens the door to Python subclassing — see below.

The B-free stiffness kernel
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The built-in formulations (``PlaneStress``, ``PlaneStrain``) implement a **B-free**
stiffness kernel instead of the classical
:math:`\mathbf{B}^T \mathbf{D} \mathbf{B}` matrix product. For a pair of basis-function
physical-space gradients :math:`\nabla R_a` and :math:`\nabla R_b`, the elementary
stiffness block at a Gauss point is:

.. math::

   \mathbf{K}^{ab} = \lambda\,(\nabla R_a \otimes \nabla R_b^T)
                   + \mu\,(\nabla R_b \otimes \nabla R_a^T)
                   + \mu\,(\nabla R_a \cdot \nabla R_b)\,\mathbf{I}

where :math:`\lambda` and :math:`\mu` are the Lamé constants of the chosen formulation
(:math:`\lambda_\text{PS} = \nu E/(1-\nu^2)` for plane stress,
:math:`\lambda_\text{PE} = \nu E/((1+\nu)(1-2\nu))` for plane strain, :math:`\mu`
identical). For linear isotropic elasticity this is mathematically equivalent to
:math:`\mathbf{B}^T \mathbf{D} \mathbf{B}` — the two formulations produce the same
global stiffness matrix.

The practical advantages over the B-matrix approach are:

- **Dimension-agnostic**: the formula holds for 2D, 3D, shell, and axisymmetric problems
  without changing the integration loop — only the Lamé constants differ.
- **No Voigt encoding**: the outer-product form avoids assembling the ``B`` matrix and
  the Voigt-packed ``D`` matrix entirely; for an isotropic law the number of operations
  per Gauss point is :math:`O(n^2)` instead of :math:`O(n^2 d^2)` where :math:`n` is
  ``n_dofs_per_cp`` and :math:`d` the spatial dimension.
- **Natural extension point**: anisotropic or path-dependent materials (J2 plasticity,
  thermoelasticity) implement ``stiffness_density`` with their own algebra, with no
  changes to the integration loop and no Voigt convention to maintain.

Python subclassing via pybind11 trampoline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:class:`~yeti_iga.future.bspline.ConstitutiveLaw` is the third class in ``future``
exposed through a pybind11 trampoline (after
:class:`~yeti_iga.future.bspline.LocalOperator` and
:class:`~yeti_iga.future.bspline.Traction`). Subclass it from Python and override
``n_dofs_per_cp()`` and ``stiffness_density(grad_a, grad_b, x_phys)`` to define a
custom material model — the integration loop in ``PatchIntegrator`` calls these methods
via the virtual dispatch and weights the result by the Gauss weight and ``|detJ|``
before accumulating it into the local stiffness matrix.

One implementation detail worth noting: ``ConstitutiveLaw``'s constructor is
``protected`` (it is abstract — direct instantiation makes no sense). pybind11's
standard ``using Base::Base;`` inheritance keeps the same access level, making the
constructor inaccessible from the binding layer. The trampoline therefore defines an
explicit ``public`` constructor that delegates to the protected base, which is the
canonical fix for this pattern.

Solution evaluation and scalar integrals
------------------------------------------

Evaluating a FE solution at arbitrary points
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Once a linear system is solved, the displacement field
:math:`u_h(\xi) = \sum_a R_a(\xi)\,d_a` must be evaluated at a set of parametric
points — typically to generate a deformed-mesh plot or to compute error norms against an
analytical solution. :class:`~yeti_iga.future.bspline.PatchEvaluator` handles this:
given a ``Patch`` with a ``PatchDOFManager`` and a global solution vector ``u_global``,
``evaluate_solution(params, u_global)`` returns the field values at every parametric
point in ``params`` (shape ``(n_pts, n_dofs_per_cp)``).

The evaluation loop is OpenMP-parallel over points and NURBS-aware: the same
``if constexpr (IsRational)`` dispatch used by ``PatchIntegrator`` ensures that B-spline
patches never pay any rationalisation cost.

Post-processing scalar quantities: ScalarLocalOperator
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:class:`~yeti_iga.future.bspline.ScalarLocalOperator` addresses a different
post-processing need: computing a single scalar over the patch (e.g. an
:math:`L^2` or :math:`H^1` error norm, an energy functional) that depends on both the
geometry and the current FE solution :math:`u_h`. It is the scalar counterpart of
:class:`~yeti_iga.future.bspline.LocalOperator`: subclass it from Python and override
``compute_scalar_integrand(R, dRdx, dRdy, physical_point, u_local) -> float``; the
integration loop in
:meth:`PatchIntegrator.integrate_scalar_operator() <yeti_iga.future.bspline.PatchIntegrator.integrate_scalar_operator>`
supplies the rationalized basis values, physical-space gradients, physical coordinates,
and local DOF values at each Gauss point, then weights the returned scalar by the Gauss
weight and ``|detJ|`` before accumulating.

Compared to ``LocalOperator`` (which returns a matrix and is used during assembly),
``ScalarLocalOperator`` returns a float and is used during post-processing. Both rely on
the same per-Gauss-point granularity rationale: the costly Jacobian inversion and
NURBS rationalisation remain in validated C++; only the problem-specific algebra moves
to Python.

Bézier extraction's role
----------------------------

:class:`~yeti_iga.future.bspline.BezierExtractor` implements the Bézier extraction
operator of :cite:`borden_isogeometric_2011`: for each tensor-product element, it builds
the matrix ``C`` such that the active B-spline basis functions on that element are a
linear combination (via ``C``) of the standard Bernstein polynomials on a local
``[0, 1]`` (or ``[0,1]^ndim``, via Kronecker products of the 1D operators) reference
element. Results are returned as :class:`~yeti_iga.future.bspline.BezierElementND`
instances — plain data containers each holding the extraction matrix ``C``, the indices
of the active B-spline functions, and the span (element) index — iterated over all
non-degenerate elements of the patch. It operates on the same
:class:`~yeti_iga.future.bspline.Patch`/:class:`~yeti_iga.future.bspline.BSpline` data
structures as the rest of the module, but it is not load-bearing for the multipatch
assembly workflow described above — it exists to provide an alternate, per-element basis
representation (e.g. for export/interop with finite-element-style tools or alternate
assembly strategies), independent of patch sharing or refinement propagation.

Worked example
------------------

A complete, narrated run through every mechanism described on this page — building two
patches that share control points, safely refining one of them with
``protected_global_ids``, propagating a refinement across a crossed interface, rebuilding
the DOF managers, and finally reclaiming orphaned control points with ``compact()`` — is
available as a Jupyter notebook at :file:`examples/future/05_multipatch.ipynb` in the
repository.
