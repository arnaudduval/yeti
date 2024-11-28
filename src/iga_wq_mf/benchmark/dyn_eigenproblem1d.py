"""
.. Test of mecanical displacement 1D
.. Author: Fabio MADIE
.. Joaquin Cornejo added some corrections 28 nov. 2024
"""

from pysrc.lib.__init__ import *
from pysrc.lib.lib_base import createUniformOpenCurve
from pysrc.lib.lib_part import part1D
from pysrc.lib.lib_job1d import mechaproblem1D
from pysrc.lib.lib_boundary import boundaryCondition
from pysrc.lib.lib_material import mechamat

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/'
if not os.path.isdir(folder): os.mkdir(folder)

YOUNG, RHO, LENGTH = 210e9, 7800, 1
MECHAMATERIAL = mechamat({'elastic_modulus':YOUNG, 
						'elastic_limit':1e6, 
						'poisson_ratio':0.3,
						'isoHardLaw': {'name':'None'}})
MECHAMATERIAL.addDensity(RHO, isIsotropic=True)

def plot_eigenValues(problem:mechaproblem1D):
	fig, ax = plt.subplots()
	D, W = problem.compute_EigenvalueProblem()
	D = np.sqrt(D[1:])
	n = np.arange(1,len(D)+1)
	w_h = (n*np.pi/LENGTH)*np.sqrt(YOUNG/RHO)
	ax.plot(n/len(D), D/w_h)
	ax.xlabel(r'$n/N$')
	ax.ylabel(r'$\omega^{h}/\omega^{th}$')
	ax.grid(False)
	ax.xlim([0, 1]); ax.ylim([0.9, 3])
	fig.tight_layout()
	fig.savefig(folder+'eigspectrum')
	return 

def simulate():
	label_list = ['Linear', 'Quadratic', 'Cubic', 'Quartic', 'Quintic']
	fig, ax = plt.subplots()
	for i, d in enumerate(np.arange(1, 6)):
		degree, nbel = d, 512
		crv = createUniformOpenCurve(degree, nbel, LENGTH)
		modelIGA = part1D(crv, kwargs={'quadArgs': {'quadrule': 'iga', 'type': 'leg'}}) 
		boundary = boundaryCondition(modelIGA.nbctrlpts)
		boundary.add_DirichletConstTemperature(table=np.array([[1 , 0]]))
		problem = mechaproblem1D(MECHAMATERIAL, modelIGA, boundary)
		
		eigs_app = np.sqrt(problem.compute_EigenvalueProblem()[0][1:])
		neigsvals = np.arange(1, len(eigs_app) + 1)
		eigs_exact = (neigsvals*np.pi/LENGTH)*np.sqrt(YOUNG/RHO)
		ax.plot(neigsvals/len(eigs_app), eigs_app/eigs_exact, label = label_list[i])

	ax.xlabel(r'$n/N$')
	ax.ylabel(r'$\omega^{app}/\omega^{exact}$')
	ax.legend()
	ax.grid(False)
	ax.xlim([0, 1]); plt.ylim([0.9, 3])
	fig.tight_layout()
	fig.savefig(folder+'eigspectrum_all')
	return 

degree, nbel = 4, 128
curveRef = createUniformOpenCurve(degree, nbel, LENGTH)
geometry = part1D(curveRef, kwargs={'quadArgs': {'quadrule': 'iga', 'type': 'leg'}})
boundary = boundaryCondition(geometry.nbctrlpts)
problem = mechaproblem1D(MECHAMATERIAL, geometry, boundary)
plot_eigenValues(problem)
simulate()
