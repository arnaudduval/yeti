"""
.. Test of mecanical displacement 1D
.. Author: Fabio MADIE
.. Joaquin Cornejo added some corrections 28 nov. 2024
"""

from yeti_iga.pymfiga.common.material import LinearElasticity
from yeti_iga.pymfiga.common.physics import EigenProblem
from yeti_iga.pymfiga.iga.geometry import GeomdlGenerator, SinglePatch
from yeti_iga.pymfiga.iga.boundary import BoundaryCondition, DOF
from yeti_iga.pymfiga.iga.single_model import MechanicalModel
from matplotlib import pyplot as plt
from time import time
import numpy as np
import os

RESULT_FOLDER = os.path.join(os.getcwd(), "results/")
if not os.path.isdir(RESULT_FOLDER):
    os.mkdir(RESULT_FOLDER)


# Global constants
YOUNG, RHO, LENGTH = 210e9, 7800, 1.0


def simulate(degree):
    nbel = 128
    geometry = GeomdlGenerator(
        filename="line",
        geo_args={
            "degree": degree,
            "nbel": nbel,
            "parameters": {"L": LENGTH},
        },
    ).export_geometry()
    patch = SinglePatch(geometry, quadclass="gs", quadtype="legendre")
    patch.generate()
    boundary = BoundaryCondition(nbctrlpts=patch.nbctrlpts, dofs=(DOF.UX,))
    material = LinearElasticity({"elastic_modulus": YOUNG}, is_unidimensional=True)
    material.add_density(RHO, is_uniform=True)
    model = MechanicalModel(material, patch, boundary)
    frequency = np.sqrt(EigenProblem().solve(model, which="SM", k=nbel - 2)[0][1:])
    return frequency


degree_labels = ["Linear", "Quadratic", "Cubic", "Quartic", "Quintic"]
fig, ax = plt.subplots()
for label, degree in zip(degree_labels, range(1, 6)):
    start = time()
    approx_freq = simulate(degree)
    print(f"Computations in {time() - start} seconds")
    freq_indices = np.arange(1, len(approx_freq) + 1)
    exact_freq = (freq_indices * np.pi / LENGTH) * np.sqrt(YOUNG / RHO)
    ax.plot(freq_indices / len(approx_freq), approx_freq / exact_freq, label=label)

ax.set_xlabel(r"$n/N$")
ax.set_ylabel(r"$\omega^{app}/\omega^{exact}$")
ax.legend()
ax.grid(False)
ax.set_xlim((0, 1))
ax.set_ylim((0.9, 1.5))
fig.tight_layout()
fig.savefig(f"{RESULT_FOLDER}/iga_eigenspectrum.pdf")

reference = [
    16300.923854,
    32601.847708,
    48902.771562,
    65203.695416,
    81504.61927,
    97805.543123,
    114106.466977,
    130407.390831,
    146708.314685,
    163009.238539,
]
approx_freq = simulate(5)
np.testing.assert_allclose(approx_freq[:10], reference)
