"""
.. Test of basis and weights 
.. We test if functions done in python and fortran for WQ approach gives the expected results.
.. Joaquin Cornejo 
"""
from pysrc.lib.__init__ import *
from pysrc.lib.lib_base import createUniformKnotvector_Rmultiplicity, evalDersBasisCSRPy
from pysrc.lib.lib_quadrules import GaussQuadrature, WeightedQuadrature

def relativeError(array_interp, array_th, relType='inf'):
	error = array_th - array_interp
	if   relType == 'inf': arg = np.inf
	elif relType == 'fro': arg = None
	try:    relError = np.linalg.norm(error, ord=arg)/np.linalg.norm(array_th, ord=arg)
	except: relError = sp.linalg.norm(error, arg)/sp.linalg.norm(array_th, arg)
	return relError

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path)
if not os.path.isdir(folder): os.mkdir(folder)

# evalDersBasisCSRPy(0, [0, 0.2, 0.4, 0.6, 0.8, 1.0], [0, 0.1, 0.2, 1.0], isfortran=True)

fig, axs  = plt.subplots(nrows=2, ncols=2, figsize=(12, 12))
nbel_list = [2**i for i in np.arange(1, 6)]

for ax, varName in zip(np.ravel(axs), ['I00', 'I01', 'I10', 'I11']):
	for degree in range(1, 5):

		error_list = []

		for nbel in nbel_list: 
			knotvector = createUniformKnotvector_Rmultiplicity(degree, nbel)
			nb_ctrlpts = len(knotvector) - degree - 1

			# WQ
			weightedQuad = WeightedQuadrature(degree, knotvector, {'type': 1})
			weightedQuad.getQuadratureRulesInfo()
			basis, weights = weightedQuad.getDenseQuadRules()
			quadPos = weightedQuad.quadPtsPos
			[B0f, B1f] = basis; [W00f, W01f, W10f, W11f] = weights
			I00f = W00f @ B0f.T; I01f = W01f @ B1f.T
			I10f = W10f @ B0f.T; I11f = W11f @ B1f.T

			fig, ax = plt.subplots(figsize=(8, 5))
			ax.plot(weightedQuad._uniqueKV, np.zeros(len(weightedQuad._uniqueKV)), marker='s', color='k')
			weightsmatrix = weights[0].todense()
			for i in range(weightedQuad.nbctrlpts):
				ax.plot(quadPos, np.ravel(weightsmatrix[i, :]))
			fig.savefig(folder+'/weights'+'W00'+'.png')

			fig, ax = plt.subplots(figsize=(8, 5))
			ax.plot(weightedQuad._uniqueKV, np.zeros(len(weightedQuad._uniqueKV)), marker='s', color='k')
			weightsmatrix = weights[-1].todense()
			for i in range(weightedQuad.nbctrlpts):
				ax.plot(quadPos, np.ravel(weightsmatrix[i, :]))
			fig.savefig(folder+'/weights'+'W11'+'.png')

			fig, ax = plt.subplots(figsize=(8, 5))
			ax.plot(weightedQuad._uniqueKV, np.zeros(len(weightedQuad._uniqueKV)), marker='s', color='k')
			basismatrix = basis[0].todense()
			for i in range(weightedQuad.nbctrlpts):
				ax.plot(quadPos, np.ravel(basismatrix[i, :]))
			fig.savefig(folder+'/basis'+'B0'+'.png')
			
			fig, ax = plt.subplots(figsize=(8, 5))
			ax.plot(weightedQuad._uniqueKV, np.zeros(len(weightedQuad._uniqueKV)), marker='s', color='k')
			basismatrix = basis[-1].todense()
			for i in range(weightedQuad.nbctrlpts):
				ax.plot(quadPos, np.ravel(basismatrix[i, :]))
			fig.savefig(folder+'/basis'+'B1'+'.png')

			# IGA
			gaussQuad = GaussQuadrature(degree, knotvector, {})
			gaussQuad.getQuadratureRulesInfo()
			basis, weights = gaussQuad.getDenseQuadRules()
			[B0, B1] = basis
			[W00, W01, W10, W11] = weights

			I00 = W00 @ B0.T; I01 = W01 @ B1.T
			I10 = W10 @ B0.T; I11 = W11 @ B1.T

			# Compare results 
			if varName   == 'I00': var1 = I00; var2 = I00f
			elif varName == 'I01': var1 = I01; var2 = I01f 
			elif varName == 'I10': var1 = I10; var2 = I10f
			elif varName == 'I11': var1 = I11; var2 = I11f

			error = relativeError(var2, var1)
			if error > 1e-5: raise Warning('Something happend. Fortran basis are wrong')
			error_list.append(error)

		label = 'Degree $p = $ ' + str(degree)
		ax.loglog(nbel_list, error_list, '-o', label=label)

	ax.set_xlabel('Discretization level ' + r'$h^{-1}$')
	ax.set_ylabel('Relative error')

axs[-1, -1].legend(bbox_to_anchor= (1.05, 1.0), loc='upper left')
fig.tight_layout()
fig.savefig(folder + '/Error_basisweights' + '.png')