import os

os.environ["OMP_NUM_THREADS"] = "1"  # OpenMP
os.environ["OPENBLAS_NUM_THREADS"] = "1"  # OpenBLAS
os.environ["MKL_NUM_THREADS"] = "1"  # Intel MKL
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"  # Accelerate (macOS)
os.environ["NUMEXPR_NUM_THREADS"] = "1"  # NumExpr
from yeti_iga.pymfiga.common.material import LinearElasticity
from yeti_iga.pymfiga.iga.geometry import GeomdlGenerator, SinglePatch
from yeti_iga.pymfiga.iga.boundary import BoundaryCondition, DOF
from yeti_iga.pymfiga.iga.single_model import MechanicalModel
from typing import Literal
from time import time
import pandas as pd
import numpy as np

FILEPATH = os.path.dirname(os.path.realpath(__file__))
FOLDER2DATA = f"{FILEPATH}/data"
if not os.path.isdir(FOLDER2DATA):
    os.mkdir(FOLDER2DATA)


def simulation(degree: int, nbel: int, quadclass: Literal["gs", "wq"]):
    geo_args = {
        "degree": degree,
        "nbel": nbel,
    }
    geometry = GeomdlGenerator(filename="cube", geo_args=geo_args).export_geometry()
    patch = SinglePatch(geometry, quadclass=quadclass)
    patch.generate()

    material = LinearElasticity({"elastic_modulus": 1.0, "poisson_ratio": 0.3})

    boundary = BoundaryCondition(
        nbctrlpts=patch.nbctrlpts, dofs=(DOF.UX, DOF.UY, DOF.UZ)
    )

    model = MechanicalModel(material, patch, boundary)

    v = np.random.random(model.get_size_of_arrays())
    model.compute_mf_stiffness(v)

    start = time()
    model.compute_mf_stiffness(v)
    finish = time()
    return finish - start


results = []
degree_list = np.arange(1, 7, dtype=int)
nbel_list = np.array([2**c for c in range(3, 7)]).astype(int)
for n in nbel_list:
    for p in degree_list:
        t = simulation(degree=p, nbel=n, quadclass="wq")
        results.append({"degree": p, "nbel": n, "time": t})

        filename = f"{FOLDER2DATA}/cpu_time_elliptic_matvec.csv"
        df = pd.DataFrame(results)
        df.to_csv(filename, index=False)
