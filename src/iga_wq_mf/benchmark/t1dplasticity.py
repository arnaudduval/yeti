"""
.. Test of elastoplasticity 1D
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
YOUNG, CST, LENGTH  = 2e11, 3.5e8, 1
NBSTEPS = 251
TIME_LIST = np.linspace(0, np.pi, NBSTEPS)
MATARGS   = {'elastic_modulus':YOUNG, 'elastic_limit':1e8, 'plasticLaw': {'Isoname': 'linear', 'Eiso':YOUNG/10}}
isReference = False

def forceVol(P:list):
	force = CST*(P - 1/10*P**2)
	return force

def simulate(degree, nbel, args, step=-2):
	crv = createUniformCurve(degree, nbel, LENGTH)
	modelPhy = mechamat1D(crv, args)
	modelPhy.activate_mechanical(MATARGS)
	modelPhy.add_DirichletCondition(table=[1, 1])
	Fref = np.atleast_2d(modelPhy.compute_volForce(forceVol)).transpose()
	Fext = np.kron(Fref, np.sin(TIME_LIST))
	blockPrint()
	displacement = modelPhy.solve(Fext=Fext[:, :step+1])[0]
	enablePrint()
	return modelPhy, displacement


if isReference:

	degree, nbel = 2, 4096
	args = {'quadArgs': {'quadrule': 'iga', 'type': 'leg'}}
	modelPhy, displacement = simulate(degree, nbel, args)
	np.save(folder + 'dispel', displacement)
	with open(folder + 'refpartpl.pkl', 'wb') as outp:
		pickle.dump(modelPhy, outp, pickle.HIGHEST_PROTOCOL)
	
else: 

	disp_ref = np.load(folder + 'disppl.npy')
	with open(folder + 'refpartpl.pkl', 'rb') as inp:
		part_ref = pickle.load(inp)

	degree_list = np.arange(1, 4)
	cuts_list   = np.arange(2, 9)
	error_list  = np.zeros(len(cuts_list))

	fig, ax = plt.subplots()
	for degree in degree_list:
		for j, cuts in enumerate(cuts_list):
			nbel = 2**cuts
			args = {'quadArgs': {'quadrule': 'iga', 'type': 'leg'}}
			step = 96
			modelPhy, displacement = simulate(degree, nbel, args, step=step)
			error_list[j] = modelPhy.normOfError(displacement[:, -1], normArgs={'type':'H1', 'part_ref':part_ref, 
																			'u_ref': disp_ref[:, step]})		

		ax.loglog(2**cuts_list, error_list, label='degree '+str(degree), marker='o')
		ax.set_ylabel(r'$H^1$'+ ' Relative error (\%)')
		ax.set_xlabel('Number of elements')
		ax.set_ylim(bottom=1e-5, top=1e0)
		ax.set_xlim(left=1, right=10**3)

		ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
		fig.tight_layout()
		fig.savefig(folder + 'FigPlasticity' + str(step) +'.pdf')