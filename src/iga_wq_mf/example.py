
"""
.. TEST OF ELASTICITY 2D.
.. Infinite plate with a hole under uniaxial traction. 
.. The analytical solution of this problem is given by Timoshenko.
.. Joaquin Cornejo 
"""

from pysrc.lib.__init__ import *
from pysrc.lib.lib_base import vtk2png
from pysrc.lib.lib_geomdl import Geomdl
from pysrc.lib.lib_part import part
from pysrc.lib.lib_material import mechamat
from pysrc.lib.lib_boundary import boundaryCondition
from pysrc.lib.lib_job3d import mechaproblem

# Set global variables
TRACTION, RINT, REXT = 1.0, 1.0, 2.0
YOUNG, POISSON, ELASTICLIM = 1e3, 0.3, 1e10
MATARGS = {'elastic_modulus':YOUNG, 
			'elastic_limit':ELASTICLIM, 
			'poisson_ratio':POISSON,
			'isoHardLaw':{'name':'none'}}

# Define external force
def forceSurf_infPlate(position, traction=TRACTION, radius_int=RINT):
	assert isinstance(position, np.ndarray), 'Variable needs to be an array'
	nnz = np.size(position, axis=1)
	x = position[0, :]; y = position[1, :]
	r_square = x**2 + y**2
	b = radius_int**2/r_square
	theta = np.arcsin(y/np.sqrt(r_square))

	F = np.zeros((2, nnz))
	F[0, :] = traction/2*(2*np.cos(theta) - b*(2*np.cos(theta) 
						+ 3*np.cos(3*theta)) + 3*b**2*np.cos(3*theta))
	F[1, :] = traction/2*3*np.sin(3*theta)*(b**2 - b)
	return F

# Define YETI class
degree, cuts = 6, 6
quadrature_args = {'quadrule': 'iga', 'type': 'leg'}
geometry_args = {'name': 'QA', 'degree': degree*np.ones(3, dtype=int), 
				'nb_refinementByDirection': cuts*np.ones(3, dtype=int), 
				'extra':{'Rin':RINT, 'Rex':REXT}
				}
mechanical_material = mechamat(MATARGS)
yetimodel = Geomdl(geometry_args).getIGAParametrization()
yetipart = part(yetimodel, quadArgs=quadrature_args)

# Set Dirichlet boundaries
boundary_condition = boundaryCondition(yetipart.nbctrlpts)
dirichlet_table = np.zeros((2, 2, 2), dtype=int)
dirichlet_table[1, 1, 0] = 1; dirichlet_table[1, 0, 1] = 1
boundary_condition.add_DirichletDisplacement(table=dirichlet_table)

# Solve elastic problem
static_problem = mechaproblem(mechanical_material, yetipart, boundary_condition)
force_surf = static_problem.compute_surfForce(forceSurf_infPlate, nbFacePosition=1)[0]
displacement = static_problem._solveLinearizedElasticityProblem(force_surf)[0]

# Postprocessing
strain_voigt = static_problem.interpolate_strain(displacement)
strain_tensor = np.zeros((6, static_problem.part.nbqp_total))
strain_tensor[:2,:] = strain_voigt[:2,:]; strain_tensor[3,:] = strain_voigt[-1,:]
stress_vonmises = static_problem.mechamaterial.evalElasticStress(strain_tensor)
static_problem.part.postProcessingDual(fields={'vms':stress_vonmises}, name='stress')
vtk2png(filename='stress', fieldname='vms', cmap='coolwarm', title='Von Mises stress')
