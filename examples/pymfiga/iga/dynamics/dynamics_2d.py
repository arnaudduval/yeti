from yeti_iga.pymfiga.common.material import J2General
from yeti_iga.pymfiga.common.physics import ExplicitLinearDynamics, EigenProblem
from yeti_iga.pymfiga.common.io import IgaPostprocessing
from yeti_iga.pymfiga.iga.geometry import GeomdlGenerator, SinglePatch
from yeti_iga.pymfiga.iga.boundary import BoundaryCondition, DOF, ParametricDirection, BoundarySide
from yeti_iga.pymfiga.iga.single_model import MechanicalModel, ExplicitDynamicsModel
import numpy as np
from time import time
import os

# Create subfolder
RESULT_FOLDER = os.path.join(os.getcwd(), "results/")
if not os.path.isdir(RESULT_FOLDER):
    os.mkdir(RESULT_FOLDER)

SUBFOLDER = f"{RESULT_FOLDER}/iga_dynamics_2d/"
if not os.path.isdir(SUBFOLDER):
    os.mkdir(SUBFOLDER)

# Set global variables
DEGREE, NBEL = 4, 32
TRACTION, RINT, REXT = 1.0, 0.2, 1.0
YOUNG, POISSON, RHO = 2.0e3, 0.25, 7.8e-6
WAVEVEL = np.sqrt(YOUNG / RHO * (1 - POISSON) / ((1 + POISSON) * (1 - 2 * POISSON)))


def surface_force(args: dict):
    position = args["position"]
    x = position[0]
    y = position[1]

    r_square = x**2 + y**2
    b = RINT**2 / r_square
    b2 = b**2

    r = np.sqrt(r_square)
    cos_theta = x / r
    sin_theta = y / r

    cos_3theta = 4 * cos_theta**3 - 3 * cos_theta
    sin_3theta = 3 * sin_theta - 4 * sin_theta**3

    force = np.zeros_like(position)
    force[0] = (
        TRACTION
        / 2
        * (2 * cos_theta - b * (2 * cos_theta + 3 * cos_3theta) + 3 * b2 * cos_3theta)
    )
    force[1] = TRACTION / 2 * 3 * sin_3theta * (b2 - b)
    return force


# Create geometry
geo_args = {
    "degree": DEGREE,
    "nbel": NBEL,
    "parameters": {"Rin": RINT, "Rex": REXT},
}
geometry = GeomdlGenerator(
    filename="quarter_annulus", geo_args=geo_args
).export_geometry()
patch = SinglePatch(geometry, quadclass="gs")
patch.generate()

# Set Dirichlet boundaries
boundary = BoundaryCondition(nbctrlpts=patch.nbctrlpts, dofs=(DOF.UX, DOF.UY))
boundary.add_constraint(
    constraint_info=[
        {
            "direction": ParametricDirection.ETA,
            "face": BoundarySide.MAX,
            "dofs": (DOF.UX,),
        },
        {
            "direction": ParametricDirection.ETA,
            "face": BoundarySide.MIN,
            "dofs": (DOF.UY,),
        },
    ],
    constraint_type="dirichlet",
)

# Set material
material = J2General({"elastic_modulus": YOUNG, "poisson_ratio": POISSON})
material.add_density(RHO, is_uniform=True)

# Set mechanical model
model = MechanicalModel(material, patch, boundary)
timespan = (REXT - RINT + np.pi / 2 * REXT) / WAVEVEL
freqmax = np.sqrt(np.max(EigenProblem().solve(model, which="LM", k=2)[0]))
nbsteps_min = int(1.01 * np.ceil(timespan * freqmax / 2))

# Create external force
time_list = np.linspace(0, timespan, nbsteps_min)
force_ref = model.assemble_surface_force(
    {
        DOF.ALL: (
            {
                "direction": ParametricDirection.XI,
                "face": BoundarySide.MAX,
            },
            surface_force,
        )
    }
)
external_force = np.tensordot(np.ones_like(time_list), force_ref, axes=0)
displacement = np.zeros_like(external_force)
acceleration = np.zeros_like(external_force)

# Solve linear dynamics with newmark scheme
model = ExplicitDynamicsModel(material, patch, boundary)
start = time()
ExplicitLinearDynamics().solve(model, displacement, external_force, time_list)
finish = time()
print("Solve problem in %.2f seconds" % (finish - start))


# Post processing
def export_vtk(disp, name):
    strain_3d = model.interpolate_strain(disp, convert_to_3d=True)
    stress_3d = material.eval_elastic_stress(strain_3d)
    vmstress = material.eval_von_mises_stress(stress_3d)
    IgaPostprocessing.export_patch(
        model.part, "dual", folder=SUBFOLDER, fields={"vms": vmstress}, filename=name
    )
    return


export_vtk(displacement[1], f"dyn2d_{0}")
counter = 0
for counter, idx in enumerate(range(2, len(time_list) - 1, 15)):
    export_vtk(displacement[idx], f"dyn2d_{counter+1}")
export_vtk(displacement[-1], f"dyn2d_{counter+1}")
