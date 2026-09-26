import os

os.environ["OMP_NUM_THREADS"] = "1"  # OpenMP
os.environ["OPENBLAS_NUM_THREADS"] = "1"  # OpenBLAS
os.environ["MKL_NUM_THREADS"] = "1"  # Intel MKL
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"  # Accelerate (macOS)
os.environ["NUMEXPR_NUM_THREADS"] = "1"  # NumExpr
from yeti_iga.pymfiga.iga.fastdiagonalization import SingleFD
from yeti_iga.pymfiga.iga.geometry.primitives import make_uniform_knotvector
from yeti_iga.pymfiga.common.numerics.quadrature_rules import StandardGauss
from time import time
import pandas as pd
import numpy as np

FILEPATH = os.path.dirname(os.path.realpath(__file__))
FOLDER2DATA = f"{FILEPATH}/data"
if not os.path.isdir(FOLDER2DATA):
    os.mkdir(FOLDER2DATA)


def simulation(degree: int, nbel: int, ndim=2):
    ndof = 1
    knotvector = make_uniform_knotvector(degree, nbel)
    quadrule = StandardGauss(
        degree, knotvector, quadtype="legendre"
    ).export_quadrature_rules()
    quadrule_list = [quadrule] * ndim

    fd = SingleFD()
    table_dirichlet = np.ones((ndof, ndim, 2), dtype=bool)
    fd.compute_space_eigendecomposition(quadrule_list, table_dirichlet)
    fd.compute_time_schurdecomposition(quadrule)

    v = np.random.random(
        ndof * fd.space_preconditioner.sp_nnz * fd.sptm_preconditioner.tm_nnz
    )
    fd.apply_spacetime_preconditioner(v)

    start = time()
    fd.apply_spacetime_preconditioner(v)
    finish = time()
    return finish - start


results = []
degree_list = np.arange(1, 5, dtype=int)
nbel_list = np.array([2**c for c in range(3, 7)]).astype(int)
for p in degree_list:
    for n in nbel_list:
        t = simulation(degree=p, nbel=n, ndim=2)
        results.append({"degree": p, "nbel": n, "time": t})

filename = f"{FOLDER2DATA}/cpu_time_sptm_fastdiag.csv"
df = pd.DataFrame(results)
df.to_csv(filename, index=False)
