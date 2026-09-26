import os

os.environ["OMP_NUM_THREADS"] = "1"  # OpenMP
os.environ["OPENBLAS_NUM_THREADS"] = "1"  # OpenBLAS
os.environ["MKL_NUM_THREADS"] = "1"  # Intel MKL
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"  # Accelerate (macOS)
os.environ["NUMEXPR_NUM_THREADS"] = "1"  # NumExpr
from yeti_iga.pymfiga.common.material import ThermalMaterial
from yeti_iga.pymfiga.common.io import IgaPostprocessing
from yeti_iga.pymfiga.iga.geometry import GeomdlGenerator, SinglePatch
from yeti_iga.pymfiga.iga.boundary import BoundaryCondition, DOF
from yeti_iga.pymfiga.iga.single_model import SpaceTimeThermalModel
from typing import Literal
from time import time
import pandas as pd
import numpy as np

FILEPATH = os.path.dirname(os.path.realpath(__file__))
FOLDER2DATA = f"{FILEPATH}/data"
if not os.path.isdir(FOLDER2DATA):
    os.mkdir(FOLDER2DATA)


def simulation(
    degree: int,
    nbel: int,
    quadclass: Literal["gs", "wq"],
    ndim: int,
    export_patch=False,
):
    geo_args = {
        "degree": degree,
        "nbel": nbel,
    }
    if ndim == 2:
        filename = "square"
    elif ndim == 3:
        filename = "cube"
    else:
        raise ValueError()
    geometry = GeomdlGenerator(filename=filename, geo_args=geo_args).export_geometry()
    patch = SinglePatch(geometry, quadclass=quadclass)
    patch.generate()
    if export_patch:
        IgaPostprocessing.export_patch(
            patch=patch, field_type="primal", filename=f"performance_{filename}"
        )

    geometry = GeomdlGenerator(filename="line", geo_args=geo_args).export_geometry()
    time_patch = SinglePatch(geometry, quadclass=quadclass)
    time_patch.generate()

    material = ThermalMaterial()
    material.add_capacity(1, is_uniform=True)
    material.add_conductivity(2, is_uniform=True, ndim=ndim)

    boundary = BoundaryCondition(nbctrlpts=patch.nbctrlpts, dofs=(DOF.T,))

    model = SpaceTimeThermalModel(material, patch, time_patch, boundary)

    nnz = patch.nbctrlpts_total * time_patch.nbctrlpts_total
    v = np.random.random(nnz)
    model.compute_mf_sptm_mass(v)
    model.compute_mf_sptm_stiffness(v)

    start = time()
    model.compute_mf_sptm_mass(v)
    model.compute_mf_sptm_stiffness(v)
    finish = time()
    return finish - start


results = []
degree_list = np.arange(1, 7, dtype=int)
nbel_list = np.array([2**c for c in range(3, 7)]).astype(int)
for n in nbel_list:
    for p in degree_list:
        t = simulation(degree=p, nbel=n, ndim=2, quadclass="wq", export_patch=False)
        results.append({"degree": p, "nbel": n, "time": t})

        filename = f"{FOLDER2DATA}/cpu_time_sptm_matvec.csv"
        df = pd.DataFrame(results)
        df.to_csv(filename, index=False)
