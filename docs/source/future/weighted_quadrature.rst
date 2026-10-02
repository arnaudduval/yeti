Weighted quadrature
====================

.. contents:: Table of Contents
    :depth: 2
    :local:
    :backlinks: none

This page explains, from first principles, *why* weighted quadrature (WQ) exists, *what*
problem it solves relative to standard Gauss quadrature, and *how*
:class:`~yeti_iga.future.bspline.WeightedQuadrature1D` builds it. It is a direct C++ port
of ``pymfiga``'s own
:file:`src/yeti_iga/pymfiga/common/numerics/quadrature_rules/weighted_quadrature.py`; for
the precise class/method signatures, see :doc:`api`.

The problem: Gauss quadrature doesn't know the basis is smooth
------------------------------------------------------------------

Assembling a stiffness or mass matrix in isogeometric analysis means evaluating
integrals such as

.. math::

    K_{ij} = \int_0^1 \frac{dN_i}{du}\,\frac{dN_j}{du}\; du,
    \qquad
    M_{ij} = \int_0^1 N_i\, N_j\; du

for every pair of degree-:math:`p` B-spline basis functions :math:`N_i`, :math:`N_j`.
Standard Gauss-Legendre quadrature does this element by element: on each knot span,
:math:`p+1` points are enough to integrate a degree-:math:`2p` polynomial exactly (the
degree of a product of two degree-:math:`p` functions), so the whole mesh needs

.. math::

    n_\text{Gauss} = n_\text{el} \,(p+1)

points, where :math:`n_\text{el}` is the number of elements.

This is exact, and it is also wasteful, for a reason that has nothing to do with the
integrand's polynomial degree: a B-spline basis function of degree :math:`p` is
:math:`C^{p-1}`-continuous across every non-repeated interior knot. Gauss quadrature
treats each element as if its basis functions were completely unrelated to their
neighbors' — it re-derives, independently per element, information that is actually
shared across many elements by the basis's own smoothness. The number of quadrature
points it uses scales with the number of *elements*, when the number of *independent*
degrees of freedom in a smooth spline space grows much more slowly.

Weighted quadrature, introduced by :cite:`calabro_fast_2017`, is built to exploit
exactly this: it uses a **global**, **reduced** set of points — not tied to individual
elements — and recovers the lost per-point generality by giving **each basis function
its own weight** at each point, rather than one weight shared by every basis function
(as Gauss does).

The core idea: decouple "where" from "how much"
---------------------------------------------------

In Gauss quadrature, a single quantity — the point's location :math:`u_q` and its weight
:math:`w_q` — is shared by every basis function evaluated there:
:math:`\int f\,du \approx \sum_q w_q\, f(u_q)` for *any* :math:`f`. Weighted quadrature
breaks that sharing. It keeps a small set of points :math:`u_q`, but replaces the single
scalar :math:`w_q` with a **matrix** of weights: a different number
:math:`W_i(u_q)` for every basis function :math:`N_i`. Concretely, instead of one rule
that integrates *anything*, WQ builds, **for every basis function separately**, a
custom one-row quadrature rule

.. math::

    \int_0^1 \left(\frac{d^a N_i}{du^a}\right) M\; du \;\approx\; \sum_q W^{(a)}_i(u_q)\, M(u_q)

that is built to be **exact** — not approximate — whenever :math:`M` belongs to a
chosen reference space (the *target space*, see below). Because this per-basis-function
rule only has to work for :math:`N_i`'s own (local) support, and because the exactness
requirement is against a deliberately smaller target space rather than the full product
space Gauss has to handle, far fewer points are needed overall — while every basis
function's row can reuse the same small, shared, mesh-global point set.

Building the rule, step by step
-----------------------------------

:meth:`WeightedQuadrature1D::build() <yeti_iga.future.bspline.WeightedQuadrature1D.build>`
follows the same four stages as ``pymfiga``'s own implementation:

1. Point placement
~~~~~~~~~~~~~~~~~~~~~

The WQ points are placed with a simple **midpoint rule**: a handful of evenly-spaced
points inside each knot span, with the two boundary spans treated slightly
differently from interior ones (method-dependent point counts :math:`r` and
:math:`s`). This placement is purely geometric and does not depend on which basis
function is being considered — the same point set is reused by every row's
least-squares solve in the next stage.

2. The target space: relaxing what "exact" means
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This is the crux of the method. The weights are **not** fit to reproduce integrals
exactly against the full, smooth degree-:math:`p` B-spline space the basis functions
themselves live in — the paper's result is that it suffices to be exact against a
smaller, specifically-chosen *target space*, and that doing so does not cost any
approximation order in the resulting Galerkin matrices. Two variants are implemented,
matching ``pymfiga``'s own ``"1"``/``"2"`` naming (selected via ``quadtype``):

.. list-table:: Target space by method
    :widths: 15 35 50
    :header-rows: 1

    * - Method
      - Target space
      - In Calabro-Sangalli-Tani's notation
    * - ``"2"``
      - Same degree :math:`p`, every **interior knot's multiplicity increased by
        one** relative to the original space
      - :math:`S^p_{r-1}` (one continuity class less than the original
        :math:`S^p_r`)
    * - ``"1"``
      - Degree :math:`p-1`, original knot vector with its first and last entries
        stripped
      - :math:`S^{p-1}_{r-1}`

Both target spaces are, in a precise sense, "one notch less smooth" than the original
space — method ``"2"`` by repeating interior knots (same degree, lower continuity),
method ``"1"`` by dropping the degree by one instead. This is implemented by
``increase_multiplicity_to_knotvector()`` (a direct port of ``pymfiga``'s
``Operations.increase_multiplicity_to_knotvector``) for method ``"2"``, and by simply
slicing the knot vector for method ``"1"``.

3. One small least-squares problem per basis function
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For each basis function :math:`N_i` and each derivative order :math:`a \in \{0, 1\}`,
the unknowns are the weights :math:`W^{(a)}_i(u_q)` at the (few) WQ points inside
:math:`N_i`'s support. They are determined by requiring exact reproduction of the
Gauss-quadrature value (itself exact, computed with the existing
:class:`~yeti_iga.future.bspline.IGABasis1D` machinery) of

.. math::

    \int_0^1 \left(\frac{d^a N_i}{du^a}\right) M_k\; du

for **every** basis function :math:`M_k` of the target space — one linear equation per
target-space basis function, one unknown per WQ point in :math:`N_i`'s support. This
system is typically under-determined (more candidate points than constraints), so it is
solved as the **minimum-norm** least-squares solution (via
``Eigen::CompleteOrthogonalDecomposition``, the same minimum-norm solution
``numpy.linalg.lstsq(rcond=None)`` would return) — among every weight vector that
reproduces the target integrals exactly, the one with the smallest norm is kept. This
local, per-row solve is what ``solve_wq_row()`` does, and it is repeated independently
for every basis function and for :math:`a \in \{0, 1\}`.

Isn't computing an *exact* Gauss integral, just to build a *cheaper* rule, defeating the
purpose? No — this full-Gauss evaluation
(:meth:`IGABasis1D::build() <yeti_iga.future.bspline.IGABasis1D.build>` with
:math:`p+1` points per element, line 165 of
:file:`WeightedQuadrature1D.cpp`) costs exactly as much as one ordinary Gauss-based
assembly pass — nothing more, since it reuses the same class any plain Gauss assembly in
``future`` would already build. The difference is *when* and *how often* that cost is
paid: :meth:`WeightedQuadrature1D::build() <yeti_iga.future.bspline.WeightedQuadrature1D.build>`
runs it **once per parametric direction**, while constructing the rule — not once per
``apply()`` call. Every iterative solve that follows
(:func:`~yeti_iga.future.matrix_free_solver.solve`, calling
:meth:`WQMatrixFreeStiffness::apply() <yeti_iga.future.bspline.WQMatrixFreeStiffness.apply>`
dozens of times) reuses the resulting few-point rule at its reduced cost every single
time. The one-time price of knowing the "correct" target integrals to fit against is
what buys every *later* evaluation its savings — it would be wasted effort only if the
rule, once built, were used just once.

4. Four weight matrices, one per (test, trial) derivative pair
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The result is stored as four sparse matrices
:attr:`W00, W01, W10, W11 <yeti_iga.future.bspline.WeightedQuadrature1D.W00>`
(shape :math:`n_\text{ctrlpts} \times n_\text{WQ}`), alongside the ordinary basis value
matrices
:attr:`B0, B1 <yeti_iga.future.bspline.WeightedQuadrature1D.B0>`
(shape :math:`n_\text{WQ} \times n_\text{ctrlpts}`, the original space's values/
derivatives at the WQ points — nothing special, the same role as
:class:`~yeti_iga.future.bspline.IGABasis1D`'s per-span ``N``/``dN``, just evaluated
globally). :math:`W_{\alpha\beta}` holds the weights for a **test** function
differentiated :math:`\alpha` times against a **trial** function differentiated
:math:`\beta` times:

.. list-table::
    :widths: 15 40 45
    :header-rows: 1

    * - Matrix
      - :math:`(\alpha, \beta)`
      - Role
    * - ``W00``
      - :math:`(0, 0)`
      - value · value — mass-like
    * - ``W01``
      - :math:`(0, 1)`
      - mixed
    * - ``W10``
      - :math:`(1, 0)`
      - mixed
    * - ``W11``
      - :math:`(1, 1)`
      - derivative · derivative — stiffness-like

The actual matrix-vector product
:math:`(K_{\alpha\beta})_{ij} = \int (d^\alpha N_i)(d^\beta N_j)\,du` applied to a DOF
vector :math:`d` is computed in two steps, **gather** then **scatter**, and this is the
contract every caller (:class:`~yeti_iga.future.bspline.WQMatrixFreeStiffness`,
``pymfiga``'s own ``compute_mf_scalar_gradu_gradv``) actually uses:

.. math::

    K_{\alpha\beta}\, d \;=\; W_{\alpha\beta} \,\big(B_\beta\, d\big)

— gather the trial field at the WQ points with :math:`B_\beta` (the trial side's own
derivative order :math:`\beta` selects which of ``B0``/``B1`` to use), then scatter it
back into test-function space with :math:`W_{\alpha\beta}` (whose *first* index,
:math:`\alpha`, is the test function's derivative order). For method ``"2"`` specifically
only two independent least-squares solves are needed per basis function (not four):
because its target space does not distinguish value and derivative constraints the way
method ``"1"``'s does, ``W00 == W01`` and ``W10 == W11`` by construction.

A fully worked example: degree 3, 4 elements
--------------------------------------------------

The three building blocks above are easiest to see on one small, concrete case: a
degree-3, open/clamped B-spline over 4 elements (interior knots of multiplicity 1 at
1/4, 1/2, 3/4 — boundary knots at multiplicity :math:`p+1=4`, so 7 control points):

.. code-block:: python

    kv = [0, 0, 0, 0, 0.25, 0.5, 0.75, 1, 1, 1, 1]   # 11 entries, 7 control points

**Step 1 — unique breakpoints.** ``unique_knots_01()`` collapses the repeated boundary
knots down to one occurrence each, leaving the element boundaries:

.. code-block:: python

    unique_knots_01(kv) == [0, 0.25, 0.5, 0.75, 1]   # 5 breakpoints, 4 elements

**Step 2 — the target space (method "2").**
``increase_multiplicity_to_knotvector(repeat=1, degree=3, kv)`` copies the two boundary
clusters unchanged and raises each of the 3 *interior* knots' multiplicity by one (from
1 to 2 — one continuity class less, :math:`C^2 \to C^1`, same degree):

.. code-block:: python

    kv_target == [0, 0, 0, 0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1, 1, 1, 1]
    # 14 entries, 10 control points (3 more than the original 7)

**Step 3 — WQ point placement.** ``midpoint_rule_points(s, r, unique_kv)`` places
:math:`r` points across the first and last elements, :math:`2+s` across every interior
one, then merges shared endpoints. For method ``"2"`` (:math:`s=2`, :math:`r=p+3=6`
here):

.. code-block:: python

    midpoint_rule_points(s=2, r=6, unique_kv=[0, 0.25, 0.5, 0.75, 1]) == [
        0.0, 0.05, 0.1, 0.15, 0.2, 0.25,   # 1st element: r=6 points
        0.3333, 0.4167,                   # interior element [0.25, 0.5]: 2+s=4 points,
                                           # 2 new (0.25, 0.5 already listed)
        0.5833, 0.6667,                   # interior element [0.5, 0.75]: 2 new likewise
        0.75, 0.8, 0.85, 0.9, 0.95, 1.0,   # last element: r=6 points
    ]                                      # 17 points total, matching the n_el=4 row below

Method ``"1"`` (:math:`s=1`, :math:`r=p+2=5`) places 13 points on the same mesh by the
same construction — both counts match the :math:`n_\text{el}=4` row of the table just
below, which is exactly how that table was produced: this same three-step construction,
repeated for growing element counts.

How many points does it actually save?
-------------------------------------------

For a degree-3 B-spline, growing the element count :math:`n_\text{el}`:

.. list-table::
    :widths: 15 15 20 20 15 15
    :header-rows: 1

    * - :math:`n_\text{el}`
      - Gauss
      - WQ method "1"
      - WQ method "2"
      - savings "1"
      - savings "2"
    * - 4
      - 16
      - 13
      - 17
      - 18.8 %
      - −6.2 %
    * - 8
      - 32
      - 21
      - 29
      - 34.4 %
      - 9.4 %
    * - 16
      - 64
      - 37
      - 53
      - 42.2 %
      - 17.2 %
    * - 32
      - 128
      - 69
      - 101
      - 46.1 %
      - 21.1 %
    * - 64
      - 256
      - 133
      - 197
      - 48.0 %
      - 23.0 %

(reproducible with :func:`~yeti_iga.future.bspline.WeightedQuadrature1D.build`, counting
``len(wq.quadpts)`` against ``n_el * (degree + 1)`` — see also
:file:`examples/future/14_weighted_quadrature.ipynb`). Two things stand out:

- At very coarse meshes the savings can be negative — WQ's point placement has a
  fixed-size overhead near the two domain boundaries that a handful of interior
  elements cannot yet amortize. The benefit only appears, and then grows, as the mesh is
  refined.
- Method ``"1"`` consistently places fewer points than method ``"2"``. The trade-off
  is its lower-degree (:math:`p-1`) target space, which is a strictly weaker
  constraint — this codebase's own regression suite
  (:file:`tests/future/test_weighted_quadrature.py`,
  :file:`tests/future/test_wq_stiffness.py`) and ``pymfiga``'s own validated benchmarks
  only exercise method ``"2"``, which is accordingly the default
  (``quadtype="2"``) everywhere in ``future``.

For a tensor-product 2D or 3D problem, these per-direction savings compound: a 23%
reduction per direction is already a ~40% reduction in total 2D point count
(:math:`1 - 0.77^2`), ~54% in 3D (:math:`1 - 0.77^3`).

Is the result merely approximate?
--------------------------------------

Not for the case that matters for correctness: on an **affine** patch (straight edges,
constant Jacobian — including the plain 1D case above), the WQ-assembled mass and
stiffness matrices match the Gauss-assembled ones to machine precision, for both
methods, at every degree and mesh size — not an approximation that merely converges,
an exact reproduction with far fewer points. This is because the target-space
exactness constraint, although formally weaker than exactness against the full smooth
space, turns out (Calabro-Sangalli-Tani's actual result) to already be strong enough
to reproduce these particular bilinear forms exactly when the integrand is polynomial.

The one place this stops being exact is **curved (rational/NURBS) geometry**, discussed
next.

Payoff: assembling without ever assembling
------------------------------------------------

The reason this codebase implements WQ at all is
:class:`~yeti_iga.future.bspline.WQMatrixFreeStiffness`: it applies the elasticity
stiffness operator, :math:`K\,v`, **without ever forming** the matrix :math:`K`. Per
direction, it precomputes one :class:`~yeti_iga.future.bspline.WeightedQuadrature1D`
rule and the pulled-back material/geometry tensor at every (tensor-product) WQ point;
``apply(v)`` then runs the gather/scatter product above through
``matrix_free_apply_2d`` — a 2D Kronecker-product ("sum-factorization") contraction that
never materializes the full tensor-product operator either. Reducing the WQ point count
directly reduces the cost of *every* ``apply()`` call, which is why it matters for an
iterative solver (:func:`yeti_iga.future.matrix_free_solver.solve`) that calls it dozens
of times per solve. The "material/geometry tensor" pulled back at each WQ point follows
the same B-free convention as :meth:`ConstitutiveLaw::stiffness_density()
<yeti_iga.future.bspline.ConstitutiveLaw.stiffness_density>` — see :doc:`design`'s
*Constitutive law architecture* section for why that specific tensor convention, and not
the more common Voigt one, is what has to be reproduced here.

A known caveat: only approximately symmetric on curved geometry
--------------------------------------------------------------------

:class:`~yeti_iga.future.bspline.WeightedQuadrature1D`'s gather matrices (``B0``/``B1``)
and scatter matrices (``W00``..``W11``) are fit **independently** — two separate
least-squares problems, not one that enforces them to be exact adjoints of one another.
On an affine patch this is harmless: the pulled-back ``stiffness_property`` tensor is
constant across the patch, and ``apply()`` matches
:meth:`PatchIntegrator::integrate_stiffness() <yeti_iga.future.bspline.PatchIntegrator.integrate_stiffness>`
to machine precision regardless (see
:file:`tests/future/test_wq_stiffness.py`). On a **curved** patch, that tensor varies
per WQ point and sits in between the independently-fit gather/scatter pair, which
breaks exact symmetry of ``apply()`` by a margin that shrinks under mesh refinement —
the same convergence rate as the method's own approximation error. In practice, plain
Conjugate Gradient can fail to converge on the coarsest/most-curved configurations as a
result, which is why
:func:`~yeti_iga.future.matrix_free_solver.solve` defaults to GMRES (and also accepts
BiCGSTAB) rather than CG — both tolerate the loss of exact symmetry.

For a **rational** (NURBS-weighted) patch specifically, ``apply()`` additionally has to
apply the same quotient-rule correction to the displacement field that the geometry map
itself needs (a rational basis function :math:`R_a = w_a N_a / W(\xi)` differentiates to
:math:`dR_a/d\xi = (w_a/W)\,dN_a/d\xi - R_a\,(dW/d\xi)/W`): both the trial field and the
implicit test function are rational, so the bilinear form expands into four terms
instead of one. This is a direct port of ``pymfiga``'s own
``NurbsOperations.compute_mf_scalar_gradu_gradv``, and was itself found, while building
this class, to disagree with it by about 17% on a curved quarter-ring benchmark before
this correction was added — isolated, via a control experiment forcing every
control-point weight to 1 (making the patch effectively non-rational), to exactly this
missing term.

Seeing it for yourself
----------------------------

:file:`examples/future/14_weighted_quadrature.ipynb` builds both methods' quadrature
rules side by side for a degree-3, 10-element B-spline, plots the point placement
against Gauss's, reproduces the point-count table above over a range of mesh sizes, and
cross-checks every quantity (points, basis values, all four weight matrices) against
``pymfiga``'s own Python implementation directly.
:file:`tests/future/test_weighted_quadrature.py` and
:file:`tests/future/test_wq_stiffness.py` are this codebase's regression tests for,
respectively, the quadrature rule itself and the matrix-free stiffness action built on
top of it.

Reference: :cite:`calabro_fast_2017`.
