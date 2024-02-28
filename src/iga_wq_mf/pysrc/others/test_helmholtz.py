"""
	Helmholtz equation:
	-div(u) = lambda u, over Omega
	u = 0, overall boundary Omega
	The lowest eigenvalue of this problem is related to the Poincaré constant C
	from the inequality:
	||u||_{L2} <= C ||grad(u)||_{L2}, since the better approximation C = 1/lambda
"""

from pysrc.lib.__init__ import *
from pysrc.lib.lib_geomdl import Geomdl
from pysrc.lib.lib_part import part
from pysrc.lib.lib_boundary import boundaryCondition
from pysrc.lib.lib_material import heatmat, mechamat
from pysrc.lib.lib_job import heatproblem, mechaproblem
from pysrc.lib.lib_base import solver

full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/'
if not os.path.isdir(folder): os.mkdir(folder)

# Set global variables
MATARGS = {'elastic_modulus':1e0, 'elastic_limit':1e10, 'poisson_ratio':0.3,
				'isoHardLaw': {'name':'none'}}
DEGREE, CUTS, N = 1, 3, 40
solv = solver()

def simulation(degree, cuts, quadArgs): 
	geoArgs = {'name': 'QA', 'degree': degree*np.ones(3, dtype=int), 
						'nb_refinementByDirection': cuts*np.ones(3, dtype=int), 
						'extra':{'Rin':1.0, 'Rex':2.0}
			}
	blockPrint()			
	modelGeo = Geomdl(geoArgs)
	modelIGA = modelGeo.getIGAParametrization()
	modelPhy = part(modelIGA, quadArgs=quadArgs)
	modelPhy_ref = part(modelIGA, quadArgs={'quadrule':'iga', 'type':'leg'})

	heatmaterial = heatmat()
	heatmaterial.addCapacity(inpt=1.0, isIsotropic=True)
	heatmaterial.addConductivity(inpt=1.0, isIsotropic=True, shape=2)
	elasticmaterial = mechamat(matArgs=MATARGS)

	# Set Dirichlet boundaries
	boundary = boundaryCondition(modelPhy.nbctrlpts)
	boundary.add_DirichletConstTemperature(table=np.ones((2, 2), dtype=int))
	boundary.add_DirichletDisplacement(table=np.ones((2, 2, 2), dtype=int))
	enablePrint()

	# Solve elastic problem
	heatprob_ref = heatproblem(heatmaterial, modelPhy_ref, boundary)
	mecaprob_ref = mechaproblem(elasticmaterial, modelPhy_ref, boundary)

	heatprob = heatproblem(heatmaterial, modelPhy, boundary)
	mecaprob = mechaproblem(elasticmaterial, modelPhy, boundary)

	def Capacity(x):
		x_in = np.zeros(boundary._nbctrlpts_total)
		x_in[boundary.thdof] = x
		y = heatprob_ref.compute_mfCapacity(x_in, args=modelPhy_ref.qpPhy)
		x_out = y[boundary.thdof]
		return x_out
	
	def Mass(x):
		x_in = np.zeros((boundary._dim, boundary._nbctrlpts_total))
		c = 0
		for i, dof in enumerate(boundary.mchdof):
			x_in[i, dof] = x[c:c+len(dof)]
			c += len(dof)
		y = mecaprob_ref.compute_mfMass(x_in, args=modelPhy_ref.qpPhy)
		x_out = np.array([])
		for i, dof in enumerate(boundary.mchdof):
			x_out = np.append(x_out, y[i, dof])
		return x_out

	def Conductivity(x):
		x_in = np.zeros(boundary._nbctrlpts_total)
		x_in[boundary.thdof] = x
		y = heatprob.compute_mfConductivity(x_in, args=modelPhy.qpPhy)
		x_out = y[boundary.thdof]
		return x_out
	
	def Stiffness(x):
		x_in = np.zeros((boundary._dim, boundary._nbctrlpts_total))
		c = 0
		for i, dof in enumerate(boundary.mchdof):
			x_in[i, dof] = x[c:c+len(dof)]
			c += len(dof)
		y = mecaprob.compute_mfStiffness(x_in)
		x_out = np.array([])
		for i, dof in enumerate(boundary.mchdof):
			x_out = np.append(x_out, y[i, dof])
		return x_out

	theigvals, _ = solv.eigs(N=boundary._thndof, Afun=Conductivity, Bfun=Capacity, neigvals=N, sigma=0)
	mcheigvals, _ = solv.eigs(N=boundary._mchndof, Afun=Stiffness, Bfun=Mass, neigvals=N, sigma=0)
	return theigvals, mcheigvals

# for quadrule, quadtype in zip(['iga', 'wq', 'wq'], ['leg', 1, 2]):
# 	for degree in range(1, 4):
# 		for cuts in range(2, 5):
# 			quadArgs = {'quadrule': quadrule, 'type': quadtype}
# 			theigvals, mcheigvals = simulation(degree=degree, cuts=cuts, quadArgs=quadArgs)
# 			print('degree:%d, cuts:%d, thmin:%.5e, mchmin:%.5e' %(degree, cuts, np.min(theigvals), np.min(mcheigvals)))
# 	print('*********')

theigvals_ref, mcheigvals_ref = simulation(degree=DEGREE, cuts=CUTS, quadArgs={'quadrule': 'iga', 'type': 'leg'})

fig, [ax0, ax1] = plt.subplots(nrows=1, ncols=2, figsize=(11, 5))
for quadrule, quadtype in zip(['wq', 'wq'], [1, 2]):
	quadArgs = {'quadrule': quadrule, 'type': quadtype}
	theigvals, mcheigvals = simulation(degree=DEGREE, cuts=CUTS, quadArgs=quadArgs)
	therror = np.abs(theigvals/theigvals_ref)
	mcherror = np.abs(mcheigvals/mcheigvals_ref)
	if quadtype==1: quadnumber = 4
	if quadtype==2: quadnumber = 2

	ax0.semilogy(np.linspace(0, 1, N), therror, label='IGA-WQ '+str(quadnumber))
	ax1.semilogy(np.linspace(0, 1, N), mcherror, label='IGA-WQ '+str(quadnumber))

for ax in [ax0, ax1]:
	ax.set_xlabel('Normalized number of eigenvalues')
	ax.set_ylabel(r'$\lambda/\lambda_{Gauss}$')
	ax.set_ylim(top=1.01, bottom=0.99)
	ax.legend(loc='upper right')
	
ax0.set_title('Helmholtz problem')
ax1.set_title('Elasticity problem')
fig.tight_layout()
fig.savefig(folder+'eigenvalueproblem'+str(DEGREE)+str(CUTS)+'.pdf')