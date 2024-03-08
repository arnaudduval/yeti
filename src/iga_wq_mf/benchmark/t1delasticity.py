"""
.. Test of elasticity 1D
.. Joaquin Cornejo 
"""

import pickle
from pysrc.lib.__init__ import *
from pysrc.lib.lib_base import createUniformCurve
from pysrc.lib.lib_part import part1D
from pysrc.lib.lib_1d import mechaproblem1D
from pysrc.lib.lib_material import mechamat
from pysrc.lib.lib_boundary import boundaryCondition

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/d1elastoplasticity/'
if not os.path.isdir(folder): os.mkdir(folder)

# Global variables
YOUNG, CST, LENGTH = 2e11, 4.e7, 1
MATARGS = {'elastic_modulus':YOUNG, 'elastic_limit':1e10, 'poisson_ratio':0.3,
		'isoHardLaw': {'name':'none'}}
MECHAMATERIAL = mechamat(MATARGS)
isReference = True

def forceVol(P:list):
	force = CST*(P - 1/10*P**2)
	return force

def simulate(degree, nbel, kwargs):
	geometry = createUniformCurve(degree, nbel, LENGTH)
	modelPhy = part1D(geometry, kwargs)
	boundary = boundaryCondition(modelPhy.nbctrlpts)
	boundary.add_DirichletConstTemperature(table=np.array([[1, 1]]))
	problem = mechaproblem1D(mechanical_material=MECHAMATERIAL, part=modelPhy, boundary=boundary)
	model2return = deepcopy(problem)
	Fref = np.atleast_2d(problem.compute_volForce(forceVol)).transpose()
	Fext_list = np.kron(Fref, [0, 1])
	displacement = np.zeros(np.shape(Fext_list))
	problem.solvePlasticityProblem(displacement, Fext_list)
	return model2return, displacement

if isReference:

	degree, nbel = 2, 1024
	args = {'quadArgs': {'quadrule': 'iga', 'type': 'leg'}}
	modelPhy, displacement = simulate(degree, nbel, args)
	np.save(folder + 'dispel', displacement)
	with open(folder + 'refpartel.pkl', 'wb') as outp:
		pickle.dump(modelPhy, outp, pickle.HIGHEST_PROTOCOL)
	
else: 

	def exactDisplacement(P:list):
		return CST/YOUNG*(-P**3/6 + P**4/120 + (LENGTH**2/6 - LENGTH**3/120)*P)
	
	def exactDisplacementDers(P:list):
		return CST/YOUNG*(-P**2/2 + P**3/30 + (LENGTH**2/6 - LENGTH**3/120))

	disp_ref = np.load(folder + 'dispel.npy')
	with open(folder + 'refpartel.pkl', 'rb') as inp:
		part_ref = pickle.load(inp)

	degree_list = np.arange(1, 4)
	cuts_list   = np.arange(1, 9)
	error_list  = np.zeros(len(cuts_list))

	fig, ax = plt.subplots(figsize=(9, 6))
	for degree in degree_list:
		for j, cuts in enumerate(cuts_list):
			nbel = 2**cuts
			args = {'quadArgs': {'quadrule': 'iga', 'type': 'leg'}}
			modelPhy, displacement = simulate(degree, nbel, args)
			# error_list[j], _ = modelPhy.normOfError(displacement[:, -1], normArgs={'type':'H1', 
			# 														'exactFunction': exactDisplacement, 
			# 														'exactFunctionDers': exactDisplacementDers})
			error_list[j], _ = modelPhy.normOfError(displacement[:, -1], normArgs={'type':'H1', 
																	'part_ref': part_ref, 
																	'u_ref': disp_ref[:, -1]})		

		ax.loglog(2**cuts_list, error_list, marker='s', markerfacecolor='w',
					markersize=10, linestyle='-', label='IGA-GL deg. ' + str(degree))
		ax.set_ylabel(r'$||u-u^h||_{H^1(\Omega)}$')
		ax.set_xlabel('Number of elements')
		ax.set_ylim(bottom=1e-14, top=1e-4)
		ax.set_xlim(left=1, right=10**3)

		ax.legend()
		fig.tight_layout()
		fig.savefig(folder + 'FigElasticityH1app' +'.pdf')
