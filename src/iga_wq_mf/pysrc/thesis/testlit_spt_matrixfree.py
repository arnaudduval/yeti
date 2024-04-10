from pysrc.lib.__init__ import *
from pysrc.lib.lib_base import createUniformCurve
from pysrc.lib.lib_geomdl import Geomdl
from pysrc.lib.lib_part import part, part1D
from pysrc.lib.lib_material import heatmat
from pysrc.lib.lib_boundary import boundaryCondition
from pysrc.lib.lib_stjob import stheatproblem

# Select folder
full_path = os.path.realpath(__file__)
folder2save = os.path.dirname(full_path) + '/results/biblio/'
folder2find = os.path.dirname(full_path) + '/data/'

# Set global variables
dataExist     = False
degree_list   = range(1, 11)
cuts          = 5

if not dataExist:

	timeMF_Mass_matrix = np.zeros((len(degree_list), 2))
	timeMF_Mass_matrix[:, 0] = degree_list
	timeMF_Stiff_matrix = np.zeros((len(degree_list), 2))
	timeMF_Stiff_matrix[:, 0] = degree_list
	timeMF_dersMass_matrix = np.zeros((len(degree_list), 2))
	timeMF_dersMass_matrix[:, 0] = degree_list
	timeMF_dersStiff_matrix = np.zeros((len(degree_list), 2))
	timeMF_dersStiff_matrix[:, 0] = degree_list


	for i, degree in enumerate(degree_list):
		
		geoArgs = {'name': 'CB', 'degree': degree*np.ones(3, dtype=int), 
					'nb_refinementByDirection': cuts*np.ones(3, dtype=int)
		}
		blockPrint()			
		modelGeo = Geomdl(geoArgs)
		modelIGA = modelGeo.getIGAParametrization()
		modelPhy = part(modelIGA, quadArgs={'quadrule':'iga'})
		time_spt = part1D(createUniformCurve(degree, int(2**cuts), 1.0), {'quadArgs':{'quadrule':'iga'}})

		heatmaterial = heatmat()
		heatmaterial.addCapacity(inpt=1.0, isIsotropic=True)
		heatmaterial.addConductivity(inpt=1.0, isIsotropic=True, shape=3)
		heatmaterial.addCapacityDers(inpt=1.0, isIsotropic=True)
		heatmaterial.addConductivityDers(inpt=1.0, isIsotropic=True, shape=3)
	
		# Set Dirichlet boundaries	
		sptnbctrlpts = np.array([*modelPhy.nbctrlpts[:modelPhy.dim], time_spt.nbctrlpts_total])
		boundary = boundaryCondition(sptnbctrlpts)
		dirichlet_table = np.ones((4, 2)); dirichlet_table[-1, -1] = 0
		boundary.add_DirichletConstTemperature(table=dirichlet_table)
		enablePrint()

		# Solve elastic problem
		problem = stheatproblem(heatmaterial, modelPhy, time_spt, boundary)
		v_in = np.random.random(boundary._nbctrlpts_total)

		# ------------------
		# Compute MF product
		# ------------------
		enablePrint()
		print('******')

		# start = time.time()
		# problem.compute_mfSTCapacity(v_in)
		# stop = time.time()
		# print('Time Capacity:%.2e' %(stop-start))
		# timeMF_Mass_matrix[i, 1] = stop - start
		# # np.savetxt(folder2save+'matvec_stMF_Mass_'+'.dat', timeMF_Mass_matrix)

		start = time.time()
		problem.compute_mfSTConductivity(v_in)
		stop = time.time()
		print('Time Conductivity:%.2e' %(stop-start))
		timeMF_Stiff_matrix[i, 1] = stop - start
		# np.savetxt(folder2save+'matvec_stMF_Stiff_'+'.dat', timeMF_Stiff_matrix)

		# start = time.time()
		# problem.compute_mfSTCapacityDers(v_in)
		# stop = time.time()
		# print('Time Ders Capacity:%.2e' %(stop-start))
		# timeMF_dersMass_matrix[i, 1] = stop - start
		# # np.savetxt(folder2save+'matvec_stMF_MassDers_'+'.dat', timeMF_dersMass_matrix)

		# start = time.time()
		# problem.compute_mfSTConductivityDers(v_in)
		# stop = time.time()
		# print('Time Ders Conductivity:%.2e' %(stop-start))
		# timeMF_dersStiff_matrix[i, 1] = stop - start
		# # np.savetxt(folder2save+'matvec_stMF_StiffDers_'+'.dat', timeMF_dersStiff_matrix)

else:

	fig, ax = plt.subplots()

	# Load data
	file_M = pd.read_table(folder2find + 'matvec_stMF_Mass_'  + '.dat', sep=' ', names=['degree', 'P1', 'P2']) 
	file_K = pd.read_table(folder2find + 'matvec_stMF_Stiff_' + '.dat', sep=' ', names=['degree', 'P1', 'P2']) 
	file_Mders = pd.read_table(folder2find + 'matvec_stMF_MassDers_'    + '.dat', sep=' ', names=['degree', 'P1', 'P2']) 
	file_Kders = pd.read_table(folder2find + 'matvec_stMF_StiffDers_'    + '.dat', sep=' ', names=['degree', 'P1', 'P2']) 

	degree = file_M.degree
	values = [file_M.P1+file_K.P1, file_Mders.P1+file_Kders.P1]
	labels = [  r'$\mathsf{A}x$' + ', for ' + r'$32^4$' + ' elements', 
				r'$\mathsf{B}x$' + ', for ' + r'$32^4$' + ' elements']

	for i, [value, label] in enumerate(zip(values, labels)): ax.semilogy(degree, value, '--', label=label, marker=MARKERLIST[i+1])
	ax.legend(loc='lower left')
	ax.set_xlabel('Degree ' + r'$p$')
	ax.set_ylabel('Wall time (s)')
	ax.set_xlim([1, 11])
	ax.set_ylim([1, 100])
	fig.tight_layout()
	fig.savefig(folder2save + 'sptProductMF' + '.pdf')