"""
.. Test of basis and weights 
.. We test if functions done in python and fortran for WQ approach 
.. gives the expected results.
.. Joaquin Cornejo 
"""
from lib.__init__ import *
from lib.lib_base import createKnotVector, relativeError
from lib.lib_quadrules import GaussQuadrature, WeightedQuadrature

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/test/'
if not os.path.isdir(folder): os.mkdir(folder)

nbel_list = [2**i for i in np.arange(3, 6)]

for varName in ['I00', 'I01', 'I10', 'I11']:
	
	fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 4))

	for degree in range(2, 8):

		norm_fortran = []; norm_python = []
		color = next(ax._get_lines.prop_cycler)['color']

		for nbel in nbel_list: 
			knotvector = createKnotVector(degree, nbel)
			nb_ctrlpts = len(knotvector) - degree - 1

			# --------
			# FORTRAN
			# --------
			weightedQuad = WeightedQuadrature(degree, knotvector, {'type': 1, 'extra':{'r': 3, 's': 2}})
			weightedQuad.getQuadratureRulesInfo()
			basis, weights = weightedQuad.getDenseQuadRules()
			[B0f, B1f] = basis
			[W00f, W01f, W10f, W11f] = weights

			# Calculate I
			I00f = W00f @ B0f.T; I01f = W01f @ B1f.T
			I10f = W10f @ B0f.T; I11f = W11f @ B1f.T

			# ----------
			# REFERENCE
			# ----------
			gaussQuad = GaussQuadrature(degree, knotvector, {})
			gaussQuad.getQuadratureRulesInfo()
			basis, weights = gaussQuad.getDenseQuadRules()
			[B0, B1] = basis
			[W00, W01, W10, W11] = weights

			# Calculate I
			I00 = W00 @ B0.T; I01 = W01 @ B1.T
			I10 = W10 @ B0.T; I11 = W11 @ B1.T

			# ---------------
			# Compare results 
			# ---------------
			if varName   == 'I00': var1 = I00; var2 = I00f
			elif varName == 'I01': var1 = I01; var2 = I01f 
			elif varName == 'I10': var1 = I10; var2 = I10f
			elif varName == 'I11': var1 = I11; var2 = I11f

			norm_temp = relativeError(var2, var1)
			dif = var2-var1
			if norm_temp > 1e-5: 
				raise Warning('Something happend. Fortran basis are wrong')
			norm_fortran.append(norm_temp)

		label = 'Degree $p = $ ' + str(degree)
		ax.loglog(nbel_list, norm_fortran, '-o', label=label, color=color)

	ax.set_xlabel('Discretization level ' + r'$h^{-1}$')
	ax.set_ylabel('Relative error')
	ax.legend(bbox_to_anchor= (1.05, 1.0), loc='upper left')
	fig.tight_layout()
	fig.savefig(folder + 'Error_basisweights_' + varName + '.png')