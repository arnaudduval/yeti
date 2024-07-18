from thesis.IncrementalHeat.__init__ import *
from pysrc.lib.lib_geomdl import Geomdl
from pysrc.lib.lib_part import part
from pysrc.lib.lib_material import heatmat
from pysrc.lib.lib_boundary import boundaryCondition
from pysrc.lib.lib_job3d import heatproblem

# Set global variables
RUNSIMU = False
degList = range(1, 11)
cuts    = 6

if RUNSIMU:

	timeMF_Mass_matrix = np.zeros((len(degList), 2))
	timeMF_Mass_matrix[:, 0] = degList
	timeMF_Stiff_matrix = np.zeros((len(degList), 2))
	timeMF_Stiff_matrix[:, 0] = degList
	timeMF_SM_matrix = np.zeros((len(degList), 2))
	timeMF_SM_matrix[:, 0] = degList
	timePython = np.zeros((len(degList), 2))
	timePython[:, 0] = degList

	for i, degree in enumerate(degList):
		blockPrint()
		geoArgs = {
					'name': 'CB', 'degree': degree*np.ones(3, dtype=int), 
					'nb_refinementByDirection': cuts*np.ones(3, dtype=int),
				}

		modelGeo = Geomdl(geoArgs)
		modelIGA = modelGeo.getIGAParametrization()
		modelPhy = part(modelIGA, quadArgs={'quadrule':'wq'})

		material = heatmat()
		material.addCapacity(inpt=1.0, isIsotropic=True)
		material.addConductivity(inpt=1.0, isIsotropic=True, shape=3)
	
		# Set Dirichlet boundaries
		boundary = boundaryCondition(modelPhy.nbctrlpts)
		boundary.add_DirichletConstTemperature(table=np.ones((3, 2), dtype=int))

		# Solve elastic problem
		problem = heatproblem(material, modelPhy, boundary)
		v_in = np.random.random(boundary._nbctrlpts_total)
		enablePrint()
		print('******')

		# ------------------
		# Compute MF product
		# ------------------
		start = time.process_time()
		problem.compute_mfCapacity(v_in)
		stop = time.process_time()
		print('Time Capacity:%.2e' %(stop-start))
		timeMF_Mass_matrix[i, 1] = stop - start
		# np.savetxt(FOLDER2SAVE+'matvecMF_Mass'+'.dat', timeMF_Mass_matrix)

		start = time.process_time()
		problem.compute_mfConductivity(v_in)
		stop = time.process_time()
		print('Time Conductivity:%.2e' %(stop-start))
		timeMF_Stiff_matrix[i, 1] = stop - start
		# np.savetxt(FOLDER2SAVE+'matvecMF_Stiff'+'.dat', timeMF_Stiff_matrix)

		# -------------------
		# Compute Python
		# -------------------
		"We improvise a matrix"
		univariateList = []
		for dim in range(problem.part.dim):
			quadrule = problem.part._quadraturerules[dim]
			B0 = quadrule._denseBasis[0]
			W0 = quadrule._denseWeights[0]
			univariateList.append(W0 @ B0.T)
		matrix = sp.kron(sp.kron(univariateList[0],univariateList[1]),univariateList[2])

		start = time.process_time()
		v_out = matrix @ v_in
		stop = time.process_time()
		print('Time Python:%.2e' %(stop-start))
		timePython[i, 1] = stop - start
		# np.savetxt(FOLDER2SAVE+'matvecPython'+'.dat', timePython)

else:

	fig, ax = plt.subplots()

	# Load data
	file_P = pd.read_table(FOLDER2SAVE + 'matvecPython'   + '.dat', sep=' ', names=['degree', 'P1', 'P2']) 
	file_M = pd.read_table(FOLDER2SAVE + 'matvecMF_Mass'  + '.dat', sep=' ', names=['degree', 'P1', 'P2']) 
	file_K = pd.read_table(FOLDER2SAVE + 'matvecMF_Stiff' + '.dat', sep=' ', names=['degree', 'P1', 'P2']) 

	degree = file_M.degree
	values = [file_M.P1, file_K.P1]
	labels = ['MF Capacity', 'MF Conductivity']

	ax.semilogy(file_P.degree, file_P.P1, '-', label='Built-in Python', marker=MARKERLIST[0])
	for i, [value, label] in enumerate(zip(values, labels)): 
		ax.semilogy(degree, value, '--', label=label, marker=MARKERLIST[i+1])

	ax.legend(loc='best')
	ax.set_xlabel('Degree ' + r'$p$')
	ax.set_ylabel('CPU time (s)')
	ax.set_xlim([0, 10])
	ax.set_ylim([0.01, 100])
	fig.tight_layout()
	fig.savefig(FOLDER2SAVE + 'ProductMF' + '.pdf')