"""
.. Test of elastoplasticity 1D
.. We test how elasticity module works
.. SI (Steel) : 
..      - Stress : MPa
..      - Length : mm
..      - Force  : N
..      - Mass   : metric ton 
.. Joaquin Cornejo 
"""

from lib.__init__ import *
from lib.lib_geomdl import Geomdl
from lib.lib_part import part
from lib.lib_material import mechamat
from lib.lib_boundary import boundaryCondition
from lib.lib_job import mechaproblem

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/t2delastoplasticity/'
if not os.path.isdir(folder): os.mkdir(folder)

def run_simulation(degree, cuts, matArgs, nbSteps, quadrule='iga'):

	# Create geometry
	geoName = 'QA'
	if   quadrule == 'iga': quadType = 'leg'
	elif quadrule == 'wq' : quadType = 1
	else: raise Warning('Not possible')
	geoArgs = {'name': geoName, 'degree': degree*np.ones(3, dtype=int), 
			'nb_refinementByDirection': cuts*np.ones(3, dtype=int), 
			'extra':{'Rin':5.e2, 'Rex':1.e3}}
	quadArgs  = {'quadrule': quadrule, 'type': quadType}

	modelGeo = Geomdl(geoArgs)
	modelIGA = modelGeo.getIGAParametrization()
	model    = part(modelIGA, quadArgs=quadArgs)

	# Add material
	material = mechamat(matArgs)

	# Set Dirichlet boundaries
	boundary = boundaryCondition(model.nbctrlpts)
	table = np.zeros((2, 2, 2), dtype=int)
	table[1, 1, 0] = 1
	table[1, 0, 1] = 1
	boundary.add_DirichletDisplacement(table=table)

	# Elasticity problem
	problem = mechaproblem(material, model, boundary)

	def forceSurfFun(P:list):
		x = P[0, :]
		y = P[1, :]
		ref  = np.array([2.e1, 0.0])
		prop = np.zeros((2, len(x)))
		for i in range(2): prop[i, :] = ref[i] 
		return prop
	
	Fextref = problem.eval_surfForce(forceSurfFun, nbFacePosition=1)
	Fext   = np.zeros((*np.shape(Fextref), 2*nbSteps + 1))
	for i in range(0, nbSteps+1): Fext[:, :, i] = i/nbSteps*Fextref
	for i in range(nbSteps+1, 2*nbSteps+1): Fext[:, :, i] = (2*nbSteps - i)/nbSteps*Fextref

	# Solve
	disp_cp, _, VMstress_qp = problem.solvePlasticityProblemPy(Fext=Fext)
	# plastic_cp  = model.interpolate_CntrlPtsField(plastic_qp)
	VMstress_cp = model.interpolate_CntrlPtsField(VMstress_qp) #!!!!!!!!!!!!!!!!!
	# plot_results(model.quadRule, geoArgs['length'], disp_cp, plastic_cp, stress_cp, folder=folder, method=quadrule)

	return problem, disp_cp, VMstress_cp

isReference = True
quadrule = 'wq'

# Set global variables
samplesize = 2500
nbSteps    = 50
matArgs    = {'elastic_modulus':2e5, 'elastic_limit':100, 'poisson_ratio': 0.3,
			'plasticLaw': {'name': 'swift', 'K':2e4, 'exp':0.5}}

if isReference:
	degree, cuts = 4, 7
	problem, disp_cp, VMstress_cp = run_simulation(degree, cuts, matArgs, nbSteps, quadrule='iga')
	stress_interp = problem.part.interpolateField(VMstress_cp, nbDOF=1, sampleSize=samplesize)[-1]
	# disp_interp   = model.interpolate_sampleField(disp_cp,   sampleSize=samplesize)[0]
	# np.save(folder+'disp_interp_refcont.npy', disp_interp)
	np.save(folder+'stress_interp_refcont.npy', stress_interp)

