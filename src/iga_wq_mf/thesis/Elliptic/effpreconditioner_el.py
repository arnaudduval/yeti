from thesis.Elliptic.__init__ import *
from pysrc.lib.lib_geomdl import Geomdl
from pysrc.lib.lib_part import part
from pysrc.lib.lib_material import mechamat
from pysrc.lib.lib_boundary import boundaryCondition
from pysrc.lib.lib_job3d import mechaproblem
from pysrc.lib.lib_base import solver

def simulate_el(degree, cuts, preconditioner='JMC'):
	geoArgs = {'name': GEONAME, 'degree': degree*np.ones(3, dtype=int), 
				'nb_refinementByDirection': cuts*np.ones(3, dtype=int), 
				'extra':{'Rin':RINT, 'Rex':REXT}
				}
	material = mechamat(MATARGS)

	blockPrint()
	modelGeo = Geomdl(geoArgs)
	modelIGA = modelGeo.getIGAParametrization()
	modelPhy = part(modelIGA, quadArgs={'quadrule':'wq', 'type':1})

	# Set Dirichlet boundaries
	boundary = boundaryCondition(modelPhy.nbctrlpts)
	table = np.zeros((2, 2, 2), dtype=int)
	table[1, 1, 0] = 1; table[1, 0, 1] = 1
	boundary.add_DirichletDisplacement(table=table)
	enablePrint()

	# Solve elastic problem
	problem = mechaproblem(material, modelPhy, boundary)
	Fext = problem.compute_surfForce(forceSurf_infPlate, nbFacePosition=1)[0]

	if preconditioner == 'ilu':
		stiffnessmatrix = buildmatrix_el(problem)
		residue = solvesystem_el(problem, stiffnessmatrix, np.ravel(Fext))

	else:
		problem._thresLin = 1.e-12; problem._linPreCond = preconditioner
		_, residue = problem._solveLinearizedElasticityProblem(Fext=Fext)
		
	return problem, residue

degree, cuts = 6, 6
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5.5, 5.5))
for j, preconditioner in enumerate(['WP', 'ilu', 'C', 'JMC']):
	start = time.process_time()
	problem, residue = simulate_el(degree, cuts, preconditioner=preconditioner)
	stop = time.process_time()
	print('time:%.2e'%(stop-start))

	if preconditioner == 'WP': labelfig = 'w.o. preconditioner'
	elif preconditioner == 'ilu': labelfig = 'Incomplete LU'
	elif preconditioner == 'C' : labelfig = 'Classic FD'
	elif preconditioner == 'JMC' : labelfig = 'This work'
	ax.semilogy(residue, marker=MARKERLIST[j], label=labelfig)

ax.legend(ncol=2, bbox_to_anchor=(0.5, 1.2), loc='upper center')
ax.set_ylabel('Relative residue')
ax.set_xlabel('Number of iterations (GMRES)')
ax.set_ylim([1e-12, 1e1])
ax.set_xlim([0, 100])
fig.savefig(FOLDER2SAVE+'preconditioner_el'+'.pdf')

for degree in range(4, 7):
	for cuts in range(6, 9):
		start = time.process_time()
		problem, residue = simulate_el(degree, cuts, preconditioner='JMC')
		stop = time.process_time()
		print('%d, %d, %.2f, %d' %(degree, cuts, stop-start, len(residue[residue>0.0])))