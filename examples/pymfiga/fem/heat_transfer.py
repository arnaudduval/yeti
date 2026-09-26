from yeti_iga.pymfiga.common.material import ThermalMaterial
from yeti_iga.pymfiga.common.numerics.quadrature_rules import FiniteElementQuadrature
from yeti_iga.pymfiga.fem.geometry import MeshpyGenerator, MeshConverter
from yeti_iga.pymfiga.fem.boundary import BoundaryCondition, DOF
from yeti_iga.pymfiga.fem.model import ThermalModel
from yeti_iga.pymfiga.common.io import FeaPostprocessing
import numpy as np

# Material definition
material = ThermalMaterial()
material.add_conductivity(inpt=1, is_uniform=True, ndim=2)
material.add_capacity(inpt=1, is_uniform=True)

# Geometry definition
mesh = MeshConverter(
    MeshpyGenerator(filename="square", max_volume=0.01).mesh_geometry()
)

# Boundary conditions
boundary = BoundaryCondition(mesh, dofs=(DOF.T,))
boundary.add_constraint({DOF.T: True}, constraint_type="dirichlet")

# Model definition
model = ThermalModel(
    mesh, material, FiniteElementQuadrature(interior_order=1), boundary
)


# Compute external force
def volumetric_source(pts):
    x, y = pts[0], pts[1]
    return 2 * np.pi**2 * np.sin(np.pi * x) * np.sin(np.pi * y)


heat_source = model.assemble_volumetric_force({DOF.T: volumetric_source})
heat_source_extended = np.hstack((heat_source, np.zeros_like(heat_source)))

# Solve
temperature = model.solve_linearized_system(
    heat_source_extended,
    scalar_coefs=(0, 1),
    flux_factor=0.0,
)
np.testing.assert_almost_equal(temperature.max(), 0.9818103375097156)

# Postprocessing
FeaPostprocessing.export_meshpy(
    mesh, filename="fem_heat", point_data={"temperature": temperature}
)
