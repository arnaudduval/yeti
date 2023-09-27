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
E, cst  = 2e11, 3.5e8
length  = 1
matArgs   = {'elastic_modulus':E, 'elastic_limit':1e16, 'plasticLaw': {'Isoname': 'linear', 'Eiso':2e10}}
isReference = False

def forceVol(P:list):
	force = cst*(P - 1/10*P**2)
	return force

if isReference:

	degree, nbel = 2, 4096
	crv = createUniformCurve(degree, nbel, length)
	args = {'quadArgs': {'quadrule': 'iga', 'type': 'leg'}}
	modelPhy = mechamat1D(crv, args)

	with open(folder + 'refpartel.pkl', 'wb') as outp:
		pickle.dump(modelPhy, outp, pickle.HIGHEST_PROTOCOL)

	modelPhy.activate_mechanical(matArgs)
	modelPhy.add_DirichletCondition(table=[1, 1])
	Fend = np.atleast_2d(modelPhy.compute_volForce(forceVol)).transpose()
	Fext = np.kron(Fend, [0, 1])
	disp_cp, strain_qp, stress_qp, plastic_qp, Cep_qp = modelPhy.solve(Fext=Fext)
	np.save(folder+'dispel', disp_cp)

else: 
	def exactDisplacement(P:list):
		u = cst/E*(-P**3/6 + P**4/120 + (length**2/6 - length**3/120)*P)
		return u
	
	def exactDisplacementdiff(P:list):
		u = cst/E*(-P**2/2 + P**3/30 + (length**2/6 - length**3/120))
		return u

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
			crv = createUniformCurve(degree, nbel, length)

			args = {'quadArgs': {'quadrule': 'iga', 'type': 'leg'}}
			modelPhy = mechamat1D(crv, args)

			modelPhy.activate_mechanical(matArgs)
			modelPhy.add_DirichletCondition(table=[1, 1])
			Fend = np.atleast_2d(modelPhy.compute_volForce(forceVol)).transpose()
			Fext = np.kron(Fend, [0, 1])

			blockPrint()
			disp_cp = modelPhy.solve(Fext=Fext)[0]
			enablePrint()

			# error_list[j] = modelPhy.H1NormOfError(disp_cp[:, -1], H1NormArgs={'exactFunction0': exactDisplacement, 
			# 														'exactFunction1': exactDisplacementdiff})
			error_list[j] = modelPhy.H1NormOfError(disp_cp[:, -1], H1NormArgs={'part_ref': part_ref, 
																	'u_ref': disp_ref[:, -1]})		

		ax.loglog(2**cuts_list, error_list, label='degree '+str(degree), marker='o')
		ax.set_ylabel(r'$H^1$'+ ' Relative error (\%)')
		ax.set_xlabel('Number of elements')
		ax.set_ylim(bottom=1e-10, top=1e0)
		ax.set_xlim(left=2, right=10**3)

		ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
		fig.tight_layout()
		fig.savefig(folder + 'FigElasticityH1app' +'.png')
