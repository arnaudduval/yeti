"""
.. Test of fast diagonalization
.. We verify speed of fast diagonalization as direct method to solve Sylvester system
.. TO DO:
	rebuild the functions to assemble the capacity matrix
.. Joaquin Cornejo 
"""

from pysrc.lib.__init__ import *
from pysrc.lib.lib_geomdl import Geomdl
from pysrc.lib.lib_part import part
from pysrc.lib.lib_material import thermomat
from pysrc.lib.lib_boundary import boundaryCondition
from pysrc.lib.lib_job import heatproblem

# Select folder
full_path = os.path.realpath(__file__)
folder2save = os.path.dirname(full_path) + '/results/biblio/'
folder2find = os.path.dirname(full_path) + '/data/'

# Set global variables
dataExist     = True
withReference = True
degree_list   = range(2, 3)
cuts          = 6

if not dataExist:

	timeMF_Mass_matrix = np.zeros((len(degree_list), 2))
	timeMF_Mass_matrix[:, 0] = degree_list
	timeMF_Stiff_matrix = np.zeros((len(degree_list), 2))
	timeMF_Stiff_matrix[:, 0] = degree_list
	timeMF_SM_matrix = np.zeros((len(degree_list), 2))
	timeMF_SM_matrix[:, 0] = degree_list
	timePython = np.zeros((len(degree_list), 2))
	timePython[:, 0] = degree_list

	for i, degree in enumerate(degree_list):
		
		blockPrint()

		# Create model
		geoArgs = {'name': 'CB', 'degree': degree*np.ones(3, dtype=int), 
			'nb_refinementByDirection': cuts*np.ones(3, dtype=int)}
		quadArgs  = {'quadrule': 'iga', 'type': 'leg'} 

		modelGeo = Geomdl(geoArgs)
		modelIGA = modelGeo.getIGAParametrization()
		modelPhy = part(modelIGA, quadArgs=quadArgs)

		# Add material 
		material = thermomat()
		material.addConductivity(np.array([[1, 0.5, 0.1],[0.5, 2, 0.25], [0.1, 0.25, 3]]), isIsotropic=True) 
		material.addCapacity(1.0, isIsotropic=True) 

		# Block boundaries
		boundary = boundaryCondition(modelPhy.nbctrlpts)
		boundary.add_DirichletConstTemperature(table=np.array([[1, 1], [1, 1], [1, 1]]))
		problem = heatproblem(material, modelPhy, boundary)

		nrows = modelPhy.nbctrlpts[0] - 2
		V = np.random.random(nrows**3)

		# # --------------
		# # Compute matrix !!! to add
		# # --------------
		# dof = problem.boundary.thdof 
		# CC  = problem.eval_capacity_matrix()[:, dof][dof, :]
		
		# start = time.process_time()
		# R = CC @ V
		# stop = time.process_time()
		# timePython[i, 1] = stop - start
		# np.savetxt(folder_data+'matvec_Python_'+'.dat', timePython)

		# ------------------
		# Compute MF product
		# ------------------
		start = time.process_time()
		problem.compute_mfCapacity(V, args=problem.part.qpPhy)
		stop = time.process_time()
		timeMF_Mass_matrix[i, 1] = stop - start
		# np.savetxt(folder_data+'matvec_MF_Mass_'+'.dat', timeMF_Mass_matrix)

		start = time.process_time()
		problem.compute_mfConductivity(V, args=problem.part.qpPhy)
		stop = time.process_time()
		timeMF_Stiff_matrix[i, 1] = stop - start
		# np.savetxt(folder_data+'matvec_MF_Stiff_'+'.dat', timeMF_Stiff_matrix)

		enablePrint()


else:

	fig, ax = plt.subplots()

	# Load data
	file_P = pd.read_table(folder2find + 'matvec_Python'   + '.dat', sep=' ', names=['degree', 'P1', 'P2']) 
	file_M = pd.read_table(folder2find + 'matvec_MF_Mass'  + '.dat', sep=' ', names=['degree', 'P1', 'P2']) 
	file_K = pd.read_table(folder2find + 'matvec_MF_Stiff' + '.dat', sep=' ', names=['degree', 'P1', 'P2']) 
	file_A = pd.read_table(folder2find + 'matvec_MF_SM'    + '.dat', sep=' ', names=['degree', 'P1', 'P2']) 
	degree = file_M.degree
	values = [file_M.P1, file_K.P1, file_A.P1]
	labels = [  r'$\mathsf{M}x$' + ', ' + r'$h^{-1}=$ ' + str(2**cuts), 
				r'$\mathsf{K}x$' + ', ' + r'$h^{-1}=$ ' + str(2**cuts), 
				r'$\mathsf{A}x$' + ', ' + r'$h^{-1}=$ ' + str(2**cuts)]

	ax.semilogy(file_P.degree, file_P.P1, '--', label='Python'+', '+r'$h^{-1}=$ ' + str(2**cuts), marker=markerList[0])
	for i, [value, label] in enumerate(zip(values, labels)): ax.semilogy(degree, value, '--', label=label, marker=markerList[i+1])

	ax.legend(loc='best')
	ax.set_xlabel('Degree ' + r'$p$')
	ax.set_ylabel('CPU time (s)')
	ax.set_xlim([1, 11])
	ax.set_ylim([0.01, 100])
	fig.tight_layout()
	fig.savefig(folder2save + 'ProductMF' + '.pdf')