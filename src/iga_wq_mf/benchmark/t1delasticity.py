"""
.. Test of elasticity 1D
.. Joaquin Cornejo 
"""

import pickle
from pysrc.lib.__init__ import *
from pysrc.lib.lib_base import createUniformCurve
from pysrc.lib.lib_1d import mechamat1D

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/d1elastoplasticity/'
if not os.path.isdir(folder): os.mkdir(folder)

# Global variables
YOUNG, CST, LENGTH = 2e11, 3.5e8, 1
MATARGS = {'elastic_modulus':YOUNG, 'elastic_limit':1e16, 'plasticLaw': {'Isoname': 'none'}}
isReference = False

def forceVol(P:list):
	force = CST*(P - 1/10*P**2)
	return force

def simulate(degree, nbel, args):
	crv = createUniformCurve(degree, nbel, LENGTH)
	modelPhy = mechamat1D(crv, args)
	modelPhy.activate_mechanical(MATARGS)
	modelPhy.add_DirichletCondition(table=[1, 1])
	Fref = np.atleast_2d(modelPhy.compute_volForce(forceVol)).transpose()
	Fext = np.kron(Fref, [0, 1])
	displacement = modelPhy.solve(Fext_list=Fext)[0]
	return modelPhy, displacement

if isReference:

	degree, nbel = 2, 4096
	args = {'quadArgs': {'quadrule': 'iga', 'type': 'leg'}}
	modelPhy, displacement = simulate(degree, nbel, args)
	np.save(folder + 'dispel', displacement)
	with open(folder + 'refpartpl.pkl', 'wb') as outp:
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
	cuts_list   = np.arange(2, 9)
	error_list  = np.zeros(len(cuts_list))

	fig, ax = plt.subplots()
	for degree in degree_list:
		for j, cuts in enumerate(cuts_list):
			nbel = 2**cuts
			args = {'quadArgs': {'quadrule': 'iga', 'type': 'leg'}}
			modelPhy, displacement = simulate(degree, nbel, args)
			# error_list[j] = modelPhy.normOfError(displacement[:, -1], normArgs={'type':'H1', 'exactFunction': exactDisplacement, 
			# 														'exactFunctionDers': exactDisplacementdiff})
			error_list[j] = modelPhy.normOfError(displacement[:, -1], normArgs={'type':'H1', 'part_ref': part_ref, 
																	'u_ref': disp_ref[:, -1]})		

		ax.loglog(2**cuts_list, error_list, label='degree '+str(degree), marker='o')
		ax.set_ylabel(r'$H^1$'+ ' Relative error (\%)')
		ax.set_xlabel('Number of elements')
		ax.set_ylim(bottom=1e-10, top=1e0)
		ax.set_xlim(left=1, right=10**3)

		ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
		fig.tight_layout()
		fig.savefig(folder + 'FigElasticityH1app' +'.png')
