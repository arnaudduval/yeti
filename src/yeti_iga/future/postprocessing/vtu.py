"""
VTU export of B-spline patches as Bézier elements (VTK high-order cells).

Each B-spline element is exported as:
  - 2D patch  →  VTK_BEZIER_QUADRILATERAL  (cell type 77)
  - 3D patch  →  VTK_BEZIER_HEXAHEDRON     (cell type 79)

The Bézier control points are computed element-by-element via the extraction
operator:  P_bezier = (C^e)^T @ P_active

Requires Paraview 5.9+ / VTK 9.0+ for correct high-order Bézier rendering.
"""

from __future__ import annotations

import numpy as np
from pathlib import Path

from ..bspline import BezierExtractor


# ---------------------------------------------------------------------------
# Point-ordering permutations: u-fastest Bézier CPs → VTK cell ordering
# ---------------------------------------------------------------------------
# For both VTK_BEZIER_QUADRILATERAL and VTK_BEZIER_HEXAHEDRON, VTK expects:
#   1. Corner vertices
#   2. Edge interior points  (one edge at a time, in VTK edge order)
#   3. Face interior points  (one face at a time)        ← 3D only
#   4. Volume interior points                            ← 3D only
#
# Our u-fastest flat index:  flat = iu + iv*(pu+1) [+ iw*(pu+1)*(pv+1)]
#
# Both the exact point order below AND the file-level metadata it depends on
# were confirmed directly against VTK itself (VTK 9.7.0) rather than assumed
# from documentation, after a rendered 3D solid came out visibly twisted
# (self-intersecting "pinwheel" folds) even for a trivial, non-rational,
# isotropic-degree box -- proving the bug lived in this permutation / the
# file it's written into, not in the geometry, NURBS rationalization, or
# BezierExtractor (each independently verified correct first). Two distinct,
# compounding issues were found, both fixed now:
#
# 1. The 6 face-interior points were grouped/ordered wrong (an unverified
#    guess at VTK's convention). The *actual* order --
#    iu=0, iu=pu, iv=0, iv=pv, iw=0, iw=pw, with the lower-lettered free
#    direction fastest within each face -- was read off directly from VTK's
#    own ground truth: calling the static
#    `vtkBezierHexahedron.PointIndexFromIJK(i, j, k, order)` for every grid
#    point and inverting it, cross-checked at several isotropic and
#    anisotropic degrees (not just the one case that first exposed the bug).
#
# 2. Independently, `write_bezier_patch_vtu`'s hand-rolled XML declared
#    `<VTKFile ... version="2.0">`. VTK's XML reader silently mis-parses
#    VTK_BEZIER_HEXAHEDRON connectivity under that (too old) format version
#    -- confirmed by writing the identical point/connectivity data with
#    vtkXMLUnstructuredGridWriter (which stamps version 2.3) and seeing no
#    corruption, then bisecting the version string down to "2.2" as the
#    boundary. This is also what caused the previously-diagnosed
#    "anisotropic HigherOrderDegrees misread" crash (vtkBezierHexahedron
#    segfaulting under warp_by_vector()/clip() for e.g. degree (2,2,1)) --
#    that was the SAME root cause wearing a different symptom, not a
#    separate bug; the degree-elevation clone that used to work around it
#    (_equalize_hexahedron_degrees) is no longer needed and has been
#    removed. The `HigherOrderDegrees` CellData array name is also now
#    spelled correctly (it was "HighOrderDegrees", missing "er") for the
#    same reason: self-consistent but non-canonical naming happened not to
#    be the trigger here, but there's no reason to keep a typo that
#    diverges from what VTK's own writer emits.


def _perm_2d(pu: int, pv: int) -> list[int]:
    """
    u-fastest → VTK_BEZIER_QUADRILATERAL ordering.
    perm[vtk_rank] = ufastest_flat_idx
    """
    def f(iu, iv):
        return iu + iv * (pu + 1)

    p = []
    # Corners (CCW): (0,0), (pu,0), (pu,pv), (0,pv)
    p += [f(0, 0), f(pu, 0), f(pu, pv), f(0, pv)]
    # Bottom edge interior  (iv=0,  iu=1..pu-1)
    for iu in range(1, pu):
        p.append(f(iu, 0))
    # Right edge interior   (iu=pu, iv=1..pv-1)
    for iv in range(1, pv):
        p.append(f(pu, iv))
    # Top edge interior     (iv=pv, iu=1..pu-1)  ← VTK uses ascending i
    for iu in range(1, pu):
        p.append(f(iu, pv))
    # Left edge interior    (iu=0,  iv=1..pv-1)  ← VTK uses ascending j
    for iv in range(1, pv):
        p.append(f(0, iv))
    # Face interior         (iv=1..pv-1, iu=1..pu-1)
    for iv in range(1, pv):
        for iu in range(1, pu):
            p.append(f(iu, iv))
    return p


def _perm_3d(pu: int, pv: int, pw: int) -> list[int]:
    """
    u-fastest → VTK_BEZIER_HEXAHEDRON ordering.
    perm[vtk_rank] = ufastest_flat_idx

    VTK hex vertex layout (iw=0 bottom, iw=pw top):
        3----2   7----6
        |    |   |    |
        0----1   4----5
    (iu,iv,iw): 0=(0,0,0), 1=(pu,0,0), 2=(pu,pv,0), 3=(0,pv,0),
                4=(0,0,pw), 5=(pu,0,pw), 6=(pu,pv,pw), 7=(0,pv,pw)
    """
    def f(iu, iv, iw):
        return iu + iv * (pu + 1) + iw * (pu + 1) * (pv + 1)

    p = []
    # 8 corners
    p += [f(0, 0, 0), f(pu, 0, 0), f(pu, pv, 0), f(0, pv, 0),
          f(0, 0, pw), f(pu, 0, pw), f(pu, pv, pw), f(0, pv, pw)]

    # 12 edge interiors (in VTK edge order)
    # Bottom face (iw=0): edges 0→1, 1→2, 2→3, 3→0 (all ascending per VTK impl)
    for iu in range(1, pu):     p.append(f(iu, 0, 0))
    for iv in range(1, pv):     p.append(f(pu, iv, 0))
    for iu in range(1, pu):     p.append(f(iu, pv, 0))
    for iv in range(1, pv):     p.append(f(0, iv, 0))
    # Top face (iw=pw): edges 4→5, 5→6, 6→7, 7→4 (all ascending per VTK impl)
    for iu in range(1, pu):     p.append(f(iu, 0, pw))
    for iv in range(1, pv):     p.append(f(pu, iv, pw))
    for iu in range(1, pu):     p.append(f(iu, pv, pw))
    for iv in range(1, pv):     p.append(f(0, iv, pw))
    # Vertical edges: 0→4, 1→5, 2→6, 3→7
    for iw in range(1, pw): p.append(f(0,  0,  iw))
    for iw in range(1, pw): p.append(f(pu, 0,  iw))
    for iw in range(1, pw): p.append(f(pu, pv, iw))
    for iw in range(1, pw): p.append(f(0,  pv, iw))

    # 6 face interiors, in VTK order: iu=0, iu=pu, iv=0, iv=pv, iw=0, iw=pw --
    # confirmed directly against vtkBezierHexahedron.PointIndexFromIJK (not
    # the order this code originally guessed; see design.rst for how this
    # was diagnosed). Within each face, the lower-lettered free direction
    # varies fastest.
    for iw in range(1, pw):
        for iv in range(1, pv): p.append(f(0, iv, iw))    # face iu=0
    for iw in range(1, pw):
        for iv in range(1, pv): p.append(f(pu, iv, iw))   # face iu=pu
    for iw in range(1, pw):
        for iu in range(1, pu): p.append(f(iu, 0, iw))    # face iv=0
    for iw in range(1, pw):
        for iu in range(1, pu): p.append(f(iu, pv, iw))   # face iv=pv
    for iv in range(1, pv):
        for iu in range(1, pu): p.append(f(iu, iv, 0))    # face iw=0
    for iv in range(1, pv):
        for iu in range(1, pu): p.append(f(iu, iv, pw))   # face iw=pw

    # Volume interior
    for iw in range(1, pw):
        for iv in range(1, pv):
            for iu in range(1, pu):
                p.append(f(iu, iv, iw))
    return p


# ---------------------------------------------------------------------------
# Main export function
# ---------------------------------------------------------------------------

def _bezier_cps_and_weights(elem, patch, rational, all_weights):
    """
    Compute physical Bézier CPs and (for NURBS) Bézier weights for one element.

    B-spline: P_bz = C^T @ P_active
    NURBS   : blend in homogeneous coordinates then divide:
              w_bz  = C^T @ w_active
              P_bz  = (C^T @ (w * P)_active) / w_bz
    """
    active   = list(elem.active_indices)
    P_active = np.array([patch.control_point(j) for j in active], dtype=np.float64)
    if rational:
        w_active = all_weights[active]
        w_bz     = elem.C.T @ w_active
        wP_bz    = elem.C.T @ (w_active[:, None] * P_active)
        return wP_bz / w_bz[:, None], w_bz
    else:
        return elem.C.T @ P_active, None


def _local_flat_to_global(elem_idx, local_flat, degrees, n_pts_per_dir):
    """
    Map a local u-fastest flat index to the global shared-connectivity point ID.
    """
    ndim = len(degrees)
    strides = [1] * ndim
    for d in range(1, ndim):
        strides[d] = strides[d - 1] * n_pts_per_dir[d - 1]
    global_flat = 0
    rem = local_flat
    for d in range(ndim):
        i_d = rem % (degrees[d] + 1)
        rem //= (degrees[d] + 1)
        global_flat += (elem_idx[d] * degrees[d] + i_d) * strides[d]
    return global_flat


def write_bezier_patch_vtu(
    patch,
    filename: str | Path,
    field=None,
    field_name: str = "field",
) -> None:
    """
    Write a B-spline or NURBS patch as Bézier elements to a VTU file.

    Uses **VTK_BEZIER_QUADRILATERAL** (type 77) / **VTK_BEZIER_HEXAHEDRON**
    (type 79).  NURBS patches write a ``RationalWeights`` PointData array so
    VTK applies the rational Bézier formula exactly.

    Adjacent elements share boundary point IDs (shared connectivity), which
    prevents rendering cracks between elements.

    Parameters
    ----------
    patch : Patch
        The patch to export (B-spline or NURBS).
    filename : str or Path
        Output .vtu file path.
    field : array_like, shape (n_cp,) or (n_cp, k), optional
        Scalar or vector field at the B-spline control points (u-fastest flat
        order).  Transformed to Bézier CPs via the extraction operator.
    field_name : str
        Name of the field in the VTU file.

    Notes
    -----
    Requires Paraview 5.9+ / VTK 9.0+. The VTKFile header must declare
    format version 2.2 or later -- version 2.0 (an older default some code
    uses) makes VTK's XML reader silently mis-parse higher-order
    (Bezier/Lagrange) cell connectivity, see the module-level comment above
    _perm_3d() for how this was diagnosed.
    """
    ndim = len(patch.tensor.components)
    dim_phys = len(patch.control_point(0))

    if ndim == 2:
        cell_type = 77   # VTK_BEZIER_QUADRILATERAL
    elif ndim == 3:
        cell_type = 79   # VTK_BEZIER_HEXAHEDRON
    else:
        raise NotImplementedError(f"ndim={ndim} not supported (only 2D and 3D)")

    degrees = [int(s.degree) for s in patch.tensor.components]
    n_local = 1
    for p in degrees:
        n_local *= (p + 1)

    perm  = _perm_2d(*degrees) if ndim == 2 else _perm_3d(*degrees)
    elems = BezierExtractor.extract_nd(patch)

    rational    = patch.cp_manager.is_rational
    all_weights = patch.cp_manager.weights_view() if rational else None

    n_cells = len(elems)

    # Shared connectivity: unique Bézier CP grid
    # Total unique CPs per dir d = n_elems_d * degree_d + 1
    max_eidx      = [max(e.elem_index[d] for e in elems) for d in range(ndim)]
    n_pts_per_dir = [max_eidx[d] * degrees[d] + degrees[d] + 1 for d in range(ndim)]
    n_pts_tot = 1
    for n in n_pts_per_dir:
        n_pts_tot *= n

    coords      = np.zeros((n_pts_tot, 3), dtype=np.float64)
    rat_weights = np.ones(n_pts_tot, dtype=np.float64) if rational else None

    if field is not None:
        farr   = np.asarray(field, dtype=np.float64)
        scalar = farr.ndim == 1
        if scalar:
            farr = farr[:, np.newaxis]
        field_bz = np.zeros((n_pts_tot, farr.shape[1]), dtype=np.float64)
    else:
        field_bz = None
        scalar   = True

    connectivity = np.zeros(n_cells * n_local, dtype=np.int64)

    for i, elem in enumerate(elems):
        P_bezier, w_bezier = _bezier_cps_and_weights(
            elem, patch, rational, all_weights)
        elem_idx = list(elem.elem_index)

        if field_bz is not None:
            f_act    = farr[list(elem.active_indices)]
            f_bezier = elem.C.T @ f_act

        for vtk_rank, local_flat in enumerate(perm):
            gid = _local_flat_to_global(elem_idx, local_flat, degrees, n_pts_per_dir)
            connectivity[i * n_local + vtk_rank] = gid
            coords[gid, :dim_phys] = P_bezier[local_flat]
            if rational:
                rat_weights[gid] = w_bezier[local_flat]
            if field_bz is not None:
                field_bz[gid] = f_bezier[local_flat]

    offsets = np.arange(n_local, n_cells * n_local + n_local, n_local, dtype=np.int64)
    types   = np.full(n_cells, cell_type, dtype=np.uint8)

    hod_tuple = np.zeros(3, dtype=np.int32)
    for d, deg in enumerate(degrees):
        hod_tuple[d] = deg
    hod = np.tile(hod_tuple, (n_cells, 1))

    _write_vtu_xml(
        Path(filename), n_pts_tot, n_cells,
        coords, connectivity, offsets, types, hod,
        field_bz if field is not None else None,
        field_name, scalar,
        rat_weights,
    )


# ---------------------------------------------------------------------------
# VTU XML writer (ASCII, no external dependency)
# ---------------------------------------------------------------------------

def _fmt(arr: np.ndarray) -> str:
    """Flatten array to space-separated string."""
    return " ".join(repr(float(x)) if arr.dtype.kind == 'f' else str(x)
                    for x in arr.ravel())


def _write_vtu_xml(
    path: Path,
    n_pts: int,
    n_cells: int,
    coords: np.ndarray,
    connectivity: np.ndarray,
    offsets: np.ndarray,
    types: np.ndarray,
    hod: np.ndarray,
    field_bz: np.ndarray | None,
    field_name: str,
    scalar: bool,
    rat_weights: np.ndarray | None = None,
) -> None:
    n_comp = 1 if scalar else (field_bz.shape[1] if field_bz is not None else 1)

    # Flatten hod: "pu pv pw  pu pv pw  ..." — one triple per cell
    hod_vals = " ".join(str(v) for v in hod.ravel())

    L = [
        '<?xml version="1.0"?>',
        '<VTKFile type="UnstructuredGrid" version="2.2" byte_order="LittleEndian">',
        '  <UnstructuredGrid>',
        f'    <Piece NumberOfPoints="{n_pts}" NumberOfCells="{n_cells}">',
        # HigherOrderDegrees attribute tells VTK which CellData array holds degrees.
        '      <CellData HigherOrderDegrees="HigherOrderDegrees">',
        '        <DataArray type="Int32" Name="HigherOrderDegrees"'
        ' NumberOfComponents="3" format="ascii">',
        f'          {hod_vals}',
        '        </DataArray>',
        '      </CellData>',
    ]

    # PointData: optional RationalWeights (not set as active Scalars so the user
    # field remains the default coloring array) + optional user field.
    has_point_data = (rat_weights is not None) or (field_bz is not None)
    if has_point_data:
        # Active Scalars = user field (if any), otherwise leave unset
        scalars_attr = f' Scalars="{field_name}"' if field_bz is not None else ''
        # RationalWeights attribute tells VTK which PointData array holds the
        # per-control-point weights for rational Bezier cells -- exactly like
        # HigherOrderDegrees does for CellData above. Without it, VTK silently
        # falls back to treating the cell as non-rational (plain polynomial),
        # which for a weighted control net (e.g. a circular arc) distorts the
        # geometry rather than erroring -- confirmed against the VTK/Kitware
        # "Implementation of rational Bezier cells into VTK" report.
        rational_attr = ' RationalWeights="RationalWeights"' if rat_weights is not None else ''
        L.append(f'      <PointData{scalars_attr}{rational_attr}>')
        if rat_weights is not None:
            L += [
                '        <DataArray type="Float64" Name="RationalWeights"'
                ' NumberOfComponents="1" format="ascii">',
                f'          {_fmt(rat_weights)}',
                '        </DataArray>',
            ]
        if field_bz is not None:
            arr_to_write = field_bz.ravel() if scalar else field_bz
            L += [
                f'        <DataArray type="Float64" Name="{field_name}"'
                f' NumberOfComponents="{n_comp}" format="ascii">',
                f'          {_fmt(arr_to_write)}',
                '        </DataArray>',
            ]
        L.append('      </PointData>')

    # Points
    L += [
        '      <Points>',
        '        <DataArray type="Float64" NumberOfComponents="3" format="ascii">',
        f'          {_fmt(coords)}',
        '        </DataArray>',
        '      </Points>',
    ]

    # Cells
    L += [
        '      <Cells>',
        '        <DataArray type="Int64" Name="connectivity" format="ascii">',
        f'          {_fmt(connectivity)}',
        '        </DataArray>',
        '        <DataArray type="Int64" Name="offsets" format="ascii">',
        f'          {_fmt(offsets)}',
        '        </DataArray>',
        '        <DataArray type="UInt8" Name="types" format="ascii">',
        f'          {_fmt(types)}',
        '        </DataArray>',
        '      </Cells>',
        '    </Piece>',
        '  </UnstructuredGrid>',
        '</VTKFile>',
    ]

    path.write_text('\n'.join(L))
