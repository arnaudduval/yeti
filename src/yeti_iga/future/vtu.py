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

from .bspline import BezierExtractor


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
    # Top edge interior     (iv=pv, iu=pu-1..1)  ← reversed
    for iu in range(pu - 1, 0, -1):
        p.append(f(iu, pv))
    # Left edge interior    (iu=0,  iv=pv-1..1)  ← reversed
    for iv in range(pv - 1, 0, -1):
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
    # Bottom face (iw=0): edges 0→1, 1→2, 2→3 (rev), 3→0 (rev)
    for iu in range(1, pu):     p.append(f(iu, 0, 0))
    for iv in range(1, pv):     p.append(f(pu, iv, 0))
    for iu in range(pu-1, 0, -1): p.append(f(iu, pv, 0))
    for iv in range(pv-1, 0, -1): p.append(f(0, iv, 0))
    # Top face (iw=pw): edges 4→5, 5→6, 6→7 (rev), 7→4 (rev)
    for iu in range(1, pu):     p.append(f(iu, 0, pw))
    for iv in range(1, pv):     p.append(f(pu, iv, pw))
    for iu in range(pu-1, 0, -1): p.append(f(iu, pv, pw))
    for iv in range(pv-1, 0, -1): p.append(f(0, iv, pw))
    # Vertical edges: 0→4, 1→5, 2→6, 3→7
    for iw in range(1, pw): p.append(f(0,  0,  iw))
    for iw in range(1, pw): p.append(f(pu, 0,  iw))
    for iw in range(1, pw): p.append(f(pu, pv, iw))
    for iw in range(1, pw): p.append(f(0,  pv, iw))

    # 6 face interiors (one per face of the hex)
    for iv in range(1, pv):
        for iu in range(1, pu): p.append(f(iu, iv, 0))   # face iw=0
    for iv in range(1, pv):
        for iu in range(1, pu): p.append(f(iu, iv, pw))  # face iw=pw
    for iw in range(1, pw):
        for iu in range(1, pu): p.append(f(iu, 0, iw))   # face iv=0
    for iw in range(1, pw):
        for iv in range(1, pv): p.append(f(pu, iv, iw))  # face iu=pu
    for iw in range(1, pw):
        for iu in range(pu-1, 0, -1): p.append(f(iu, pv, iw))  # face iv=pv (rev)
    for iw in range(1, pw):
        for iv in range(pv-1, 0, -1): p.append(f(0, iv, iw))   # face iu=0 (rev)

    # Volume interior
    for iw in range(1, pw):
        for iv in range(1, pv):
            for iu in range(1, pu):
                p.append(f(iu, iv, iw))
    return p


# ---------------------------------------------------------------------------
# Main export function
# ---------------------------------------------------------------------------

def write_bezier_patch_vtu(
    patch,
    filename: str | Path,
    field=None,
    field_name: str = "field",
) -> None:
    """
    Write a B-spline patch as Bézier elements to a VTU file.

    Parameters
    ----------
    patch : Patch
        The B-spline patch to export.
    filename : str or Path
        Output .vtu file path.
    field : array_like, shape (n_cp,) or (n_cp, k), optional
        Scalar or vector field at the B-spline control points (u-fastest flat
        order).  Transformed to Bézier CPs via the same extraction operator
        so Paraview can interpolate it correctly inside each element.
    field_name : str
        Name of the field in the VTU file (default ``"field"``).

    Notes
    -----
    Requires Paraview 5.9+ / VTK 9.0+ for VTK_BEZIER_QUADRILATERAL /
    VTK_BEZIER_HEXAHEDRON cells.
    """
    ndim     = len(patch.tensor.components)
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

    n_cells   = len(elems)
    n_pts_tot = n_cells * n_local

    # Physical coordinates of Bézier CPs (VTK always 3-component)
    coords = np.zeros((n_pts_tot, 3), dtype=np.float64)

    # Optional field transformed to Bézier CPs
    if field is not None:
        farr   = np.asarray(field, dtype=np.float64)
        scalar = farr.ndim == 1
        if scalar:
            farr = farr[:, np.newaxis]
        field_bz = np.zeros((n_pts_tot, farr.shape[1]), dtype=np.float64)
    else:
        field_bz = None
        scalar   = True

    for i, elem in enumerate(elems):
        P_active = np.array(
            [patch.control_point(j) for j in elem.active_indices],
            dtype=np.float64,
        )                                          # (n_local, dim_phys)
        P_bezier = elem.C.T @ P_active             # (n_local, dim_phys)
        P_vtk    = P_bezier[perm]                  # reordered for VTK

        off = i * n_local
        coords[off:off + n_local, :dim_phys] = P_vtk

        if field_bz is not None:
            f_act = farr[list(elem.active_indices)]  # (n_local, k)
            field_bz[off:off + n_local] = (elem.C.T @ f_act)[perm]

    # VTK arrays
    connectivity = np.arange(n_pts_tot, dtype=np.int64)
    offsets      = np.arange(n_local, n_pts_tot + n_local, n_local, dtype=np.int64)
    types        = np.full(n_cells, cell_type, dtype=np.uint8)

    # HighOrderDegrees: one Int32 tuple [pu, pv, pw] per cell, in <CellData>.
    hod_tuple = np.zeros(3, dtype=np.int32)
    for d, deg in enumerate(degrees):
        hod_tuple[d] = deg
    hod = np.tile(hod_tuple, (n_cells, 1))   # (n_cells, 3), Int32

    _write_vtu_xml(
        Path(filename), n_pts_tot, n_cells,
        coords, connectivity, offsets, types, hod,
        field_bz if field is not None else None,
        field_name, scalar,
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
) -> None:
    n_comp = 1 if scalar else (field_bz.shape[1] if field_bz is not None else 1)

    # L2-norm range for HighOrderDegrees (matches VTK's own writer metadata)
    norms = np.linalg.norm(hod.astype(float), axis=1)
    hod_range_min = float(norms.min())
    hod_range_max = float(norms.max())

    # Flatten hod to one row of values per cell: "pu pv pw pu pv pw ..."
    hod_vals = " ".join(str(v) for v in hod.ravel())

    L = [
        '<?xml version="1.0"?>',
        '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">',
        '  <UnstructuredGrid>',
        f'    <Piece NumberOfPoints="{n_pts}" NumberOfCells="{n_cells}">',
        # HighOrderDegrees in <CellData> — required by VTK 9.0+ / Paraview 5.9+.
        # The HigherOrderDegrees="..." attribute on <CellData> marks it as the
        # *active* attribute (like Scalars=/Vectors=); without it VTK ignores
        # the array even though it is present, and silently assumes a wrong
        # uniform degree for direction-dependent Bezier cells.
        '      <CellData HigherOrderDegrees="HighOrderDegrees">',
        f'        <DataArray type="Int32" Name="HighOrderDegrees"'
        f' NumberOfComponents="3" format="ascii"'
        f' RangeMin="{hod_range_min}" RangeMax="{hod_range_max}">',
        f'          {hod_vals}',
        '          <InformationKey name="L2_NORM_RANGE" location="vtkDataArray" length="2">',
        f'            <Value index="0">{hod_range_min}</Value>',
        f'            <Value index="1">{hod_range_max}</Value>',
        '          </InformationKey>',
        '        </DataArray>',
        '      </CellData>',
    ]

    # PointData: optional field
    if field_bz is not None:
        arr_to_write = field_bz.ravel() if scalar else field_bz
        L += [
            '      <PointData>',
            f'        <DataArray type="Float64" Name="{field_name}"'
            f' NumberOfComponents="{n_comp}" format="ascii">',
            f'          {_fmt(arr_to_write)}',
            '        </DataArray>',
            '      </PointData>',
        ]

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
