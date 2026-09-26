"""
.. Test of mecanical dynamics 1D
.. Author: Fabio MADIE
.. Joaquin Cornejo added some corrections 28 nov. 2024
"""

from yeti_iga.pymfiga.common.material import LinearElasticity
from yeti_iga.pymfiga.common.physics import ExplicitLinearDynamics, EigenProblem
from yeti_iga.pymfiga.common.physics.iga_l2projection import L2projection
from yeti_iga.pymfiga.iga.geometry import GeomdlGenerator, SinglePatch
from yeti_iga.pymfiga.iga.boundary import BoundaryCondition, DOF
from yeti_iga.pymfiga.iga.single_model import MechanicalModel, ExplicitDynamicsModel
from matplotlib import pyplot as plt
import numpy as np
import os

RESULT_FOLDER = os.path.join(os.getcwd(), "results/")
if not os.path.isdir(RESULT_FOLDER):
    os.mkdir(RESULT_FOLDER)

# Set global variables
YOUNG, RHO, LENGTH = 2.0e3, 7.8e-6, 1000
WAVEVEL = np.sqrt(YOUNG / RHO)


def compute_initial_displacement(args: dict):
    x = args["position"]
    u = np.exp(-1.0e-3 / 2 * (x - LENGTH / 2) ** 2)
    return u


# Create geometry
degree, nbel = 3, 256
geometry = GeomdlGenerator(
    filename="line",
    geo_args={"degree": degree, "nbel": nbel, "parameters": {"L": LENGTH}},
).export_geometry()
patch = SinglePatch(geometry, quadclass="gs", quadtype="legendre")
patch.generate()

# Create boundary condition
boundary = BoundaryCondition(nbctrlpts=patch.nbctrlpts, dofs=(DOF.UX,))

# Set material
material = LinearElasticity({"elastic_modulus": YOUNG}, is_unidimensional=True)
material.add_density(RHO, is_uniform=True)

# Set mechanical model
model = MechanicalModel(material, patch, boundary)
timespan = LENGTH / WAVEVEL
freqmax = np.sqrt(np.max(EigenProblem().solve(model, which="LM", k=2)[0]))
nbsteps_min = int(1.005 * np.ceil(timespan * freqmax / 2))

# Create external force
time_list = np.linspace(0, timespan, nbsteps_min)
external_force = np.zeros((len(time_list), patch.nbctrlpts_total))
displacement = np.zeros_like(external_force)

# Compute initial values
displacement[0] = L2projection().solve(
    model, compute_initial_displacement({"position": model.part.qp_phy})
)

# Solve linear dynamics with newmark scheme
model = ExplicitDynamicsModel(material, patch, boundary)
ExplicitLinearDynamics().solve(model, displacement, external_force, time_list)
initial = displacement[0]
final = displacement[-1]
err = np.linalg.norm(final - initial) / np.linalg.norm(initial)
# TODO: investigate if the error is normal
# Maybe it is due to numerical dissipation (common in explicit codes)
np.testing.assert_array_less(err, 4e-2)

# Post processing
from yeti_iga.pymfiga.common.numerics.operations import BsplineOperations

fig, ax = plt.subplots()
knots_interp = np.linspace(0, 1, 501)
for quadrule in patch.quadrule_list:
    quadrule.knots_to_sample = knots_interp

disp_interp_initial = np.ravel(
    BsplineOperations.interpolate_meshgrid(
        quadrule_list=patch.quadrule_list,
        u_ctrlpts=np.atleast_2d(displacement[0]),
    )
)
ax.plot(knots_interp * LENGTH, np.ravel(disp_interp_initial), "--", label="Initial")

disp_interp_final = np.ravel(
    BsplineOperations.interpolate_meshgrid(
        quadrule_list=patch.quadrule_list,
        u_ctrlpts=np.atleast_2d(displacement[-1]),
    )
)
ax.plot(knots_interp * LENGTH, np.ravel(disp_interp_final), ".", label="Final")

disp_exact = compute_initial_displacement({"position": knots_interp * LENGTH})
ax.plot(knots_interp * LENGTH, disp_exact, label="Exact", alpha=0.5)

ax.set_xlabel("Position")
ax.set_ylabel("Displacement")
ax.legend()
fig.savefig(f"{RESULT_FOLDER}/iga_dynamics_1d.pdf")
