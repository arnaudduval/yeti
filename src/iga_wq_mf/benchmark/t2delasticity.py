"""
.. Test of elastoplasticity 2D
.. We test how elasticity module works
.. SI (Steel) : 
..      - Stress : MPa (200e3)
..      - Length : mm
.. Joaquin Cornejo 
"""

from pysrc.lib.__init__ import *
from pysrc.lib.lib_geomdl import Geomdl
from pysrc.lib.lib_part import part
from pysrc.lib.lib_material import mechamat, block_dot_product
from pysrc.lib.lib_boundary import boundaryCondition
from pysrc.lib.lib_job import mechaproblem

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/d2delasticity/'
if not os.path.isdir(folder): os.mkdir(folder)

def forceSurf_infPlate(P:list):
	Tx, rin = 1.0, 1.0
	x = P[0, :]; y = P[1, :]; nnz = np.size(P, axis=1)
	r_square = x**2 + y**2
	theta_x2 = 2.0*np.arcsin(y/np.sqrt(r_square))
	div = rin**2/r_square # Already squared

	CartSt = np.zeros((2, 2, nnz))
	CartSt[0, 0, :] = Tx*(1.0 - div*(1.5*np.cos(theta_x2) + np.cos(2.0*theta_x2)) + 1.5*div**2*np.cos(2.0*theta_x2))
	CartSt[1, 1, :] = Tx*(-div*(0.5*np.cos(theta_x2) - np.cos(2.0*theta_x2)) - 1.5*div**2*np.cos(2.0*theta_x2))
	CartSt[0, 1, :] = CartSt[1, 0, :] = Tx*(-div*(0.5*np.sin(theta_x2) + np.sin(2.0*theta_x2)) + 1.5*div**2*np.sin(2.0*theta_x2))	
	
	F = np.zeros((2, nnz))
	F[0, :] = CartSt[0, 0, :]*np.cos(theta_x2/2.0) + CartSt[0, 1, :]*np.sin(theta_x2/2.0)
	F[1, :] = CartSt[1, 0, :]*np.cos(theta_x2/2.0) + CartSt[1, 1, :]*np.sin(theta_x2/2.0)
	return F

def exactDisplacement_infPlate(P:list):
	Tx, rin, E, nu = 1.0, 1.0, 1.e3, 0.3
	x = P[0, :]; y = P[1, :]; nnz = np.size(P, axis=1)
	r_square = x**2 + y**2
	theta_x2 = 2*np.arcsin(y/np.sqrt(r_square))
	div   = rin**2/r_square # Already squared
	a     = rin*Tx*(1.0 + nu)/(2*E)

	Polardisp = np.zeros((2, nnz))
	Polardisp[0, :] = a*(div + 1 -nu + np.cos(theta_x2)*(-div**2 + 4.0*(1.0 - nu)*div + 1))
	Polardisp[1, :] = -a*np.sin(theta_x2)*(div**2 + 2.0*(1.0 - 2.0*nu)*div + 1)

	disp = np.zeros((2, nnz))
	disp[0, :] = Polardisp[0, :]*np.cos(theta_x2/2.0) - Polardisp[1, :]*np.sin(theta_x2/2.0)
	disp[1, :] = Polardisp[0, :]*np.sin(theta_x2/2.0) + Polardisp[1, :]*np.cos(theta_x2/2.0)
	return disp

# Set global variables
E, nu = 1e3, 0.3
trueEnergy = -135/32768*np.pi/E*(1024*nu**2 + 5*nu - 1019)

quadArgs = {'quadrule': 'iga', 'type': 'leg'}
# quadArgs = {'quadrule': 'wq', 'type': 2}
matArgs  = {'elastic_modulus':E, 'elastic_limit':1e10, 'poisson_ratio': nu}
solverArgs = {'nbIterationsPCG':150, 'PCGThreshold':1e-15}

# Create model 
degree_list = np.array([2, 3, 4, 6, 8])
cuts_list   = np.arange(2, 9)
error_list  = np.ones(len(cuts_list))

fig, ax  = plt.subplots(figsize=(8, 4))
for i, degree in enumerate(degree_list):
	for j, cuts in enumerate(cuts_list):
		geoArgs = {'name': 'QA', 'degree': degree*np.ones(3, dtype=int), 
					'nb_refinementByDirection': cuts*np.ones(3, dtype=int), 
					'extra':{'Rin':1.0, 'Rex':4.0}
		}
		blockPrint()
		material = mechamat(matArgs)
		modelGeo = Geomdl(geoArgs)
		modelIGA = modelGeo.getIGAParametrization()
		modelPhy = part(modelIGA, quadArgs=quadArgs)

		# Set Dirichlet boundaries
		boundary = boundaryCondition(modelPhy.nbctrlpts)
		table = np.zeros((2, 2, 2), dtype=int)
		table[1, 1, 0] = 1
		table[1, 0, 1] = 1
		boundary.add_DirichletDisplacement(table=table)
		enablePrint()

		# Solve elastic problem
		problem = mechaproblem(material, modelPhy, boundary)
		problem.addSolverConstraints(solverArgs=solverArgs)
		Fext = problem.eval_surfForce(forceSurf_infPlate, nbFacePosition=1)
		displacement, _, stress_qp = problem.solveElasticityProblemFT(Fext=Fext)
		# error_list[j] = abs(trueEnergy -  block_dot_product(2, Fext, displacement))/trueEnergy*100
		error_list[j] = problem.L2NormOfError(exactDisplacement_infPlate, displacement)

	nbctrlpts_list = (2**cuts_list+degree)**2
	ax.loglog(nbctrlpts_list, error_list, marker=markerSet[i], label='degree p='+str(degree))

	if str(quadArgs['quadrule']) == 'iga':
		slope = np.polyfit(np.log10(nbctrlpts_list[:3]),np.log10(error_list[:3]), 1)[0]
		slope = round(slope, 1)
		annotation.slope_marker((nbctrlpts_list[2], error_list[2]), slope, 
								poly_kwargs={'facecolor': (0.73, 0.8, 1)})
	
	ax.set_ylabel(r'$\displaystyle\sqrt{\frac{\int_\Omega \left(u-u^h\right)^2 d\Omega}{\int_\Omega u^2 d\Omega}}$')
	# ax.set_ylabel('Relative error of energy ' + r'$\frac{|U-U^{h}|}{|U|}$')
	ax.set_xlabel('Total number of DOF')
	ax.set_ylim(top=1e1, bottom=1e-16)
	ax.set_xlim(left=10, right=1e5)

	ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
	fig.tight_layout()
	fig.savefig(folder + 'FigInfinitePlate_' + str(quadArgs['quadrule']) +'.png')
