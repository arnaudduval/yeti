"""
.. Test of transient heat solver
.. ATTENTION: IT ONLY WORKS WITH 'ISOTROPIC' MATERIALS
.. Joaquin Cornejo 
"""

from lib.__init__ import *
from lib.physics import *
from lib.create_geomdl import geomdlModel
from lib.fortran_mf_wq import fortran_mf_wq
from lib.base_functions import create_table_properties, sigmoid

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/test7/'
if not os.path.isdir(folder): os.mkdir(folder)

# Set global variables
dataExist   = True
geolist     = ['VB', 'CB']
method_list = ['WP', 'C', 'JMC']
# method_list = ['JMC']

if not dataExist:

	degree, cuts = 6, 5
	conductivity, capacity = 0.1, 1.0
	theta = 1.0
	time_list = np.linspace(0, 10, 41)  
	table_Kprop = create_table_properties(setKprop, prop=conductivity)
	table_Cprop = create_table_properties(setCprop, prop=capacity)     

	for geoname in geolist:
		for PCGmethod in method_list:
			filename = folder + 'ResPCG_' + geoname + '_' + PCGmethod + '.dat'        

			# Create model 
			geometry = {'degree':[degree, degree, degree]}
			modelGeo = geomdlModel(geoname, **geometry)
			modelIGA = modelGeo.export_IGAparametrization(nb_refinementByDirection=
														np.array([cuts, cuts, cuts]))
			modelPhy = fortran_mf_wq(modelIGA)

			# Add material 
			material = {'capacity':capacity, 'conductivity':conductivity*np.eye(3)}
			modelPhy._set_material(material)

			# Block boundaries
			Dirichlet = {'thermal':np.array([[1, 1], [0, 0], [0, 0]])}
			modelPhy._set_dirichlet_boundaries(Dirichlet)

			# Add constant temperature
			modelPhy._add_thermal_IBC(np.array([[0, 1], [0, 0], [0, 0]]), temperature=1.0)

			# ---------------------
			# Transient model
			# ---------------------
			# Interpolate temperature on boundaries over time 
			GBound = np.zeros((len(modelPhy._thermal_dod), len(time_list)))
			for i in range(len(time_list)): GBound[:, i] = modelPhy._get_thermal_IBC()

			# Add external force (transient)
			Fend  = modelPhy.eval_source_vector(powden)
			Fendt = np.atleast_2d(Fend).reshape(-1, 1)
			Fext  = np.kron(Fendt, sigmoid(time_list))

			# Solve
			N = 21
			Tsol, resPCG, properties = modelPhy.MFtransientHeatNL(Fext=Fext[:, :N], G=GBound[:, :N], time_list=time_list[:N],
											table_Kprop=table_Kprop, table_Cprop=table_Cprop, 
											methodPCG=PCGmethod, theta=theta)
											
			modelPhy.export_results(u_ctrlpts=properties, nbDOF=2, folder=folder+geoname+'_'+PCGmethod+'/')
			print('Finish')
			np.savetxt(filename, resPCG)
			# modelPhy.export_results(u_ctrlpts=Tsol[:, -1], folder=folder, nbDOF=1)

else:

	for geoname in geolist:
		fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(6, 4))

		for PCGmethod in method_list:
			filename = folder + 'ResPCG_' + geoname + '_' + PCGmethod + '.dat'
			resPCG = np.loadtxt(filename)

			if PCGmethod   == "WP" : labelmethod = 'w.o. preconditioner'
			elif PCGmethod == "C"  : labelmethod = 'Classic FD method'
			elif PCGmethod == "JMC": labelmethod = 'This work'

			# Print the first or last
			step = resPCG[0, 0]; iterNL = resPCG[1, 0]
			newresidue = resPCG[2:, 0]; newresidue = newresidue[newresidue>0]
			ax.semilogy(np.arange(len(newresidue)), newresidue, 'o-', linewidth=2.5,
						label=labelmethod)

			ax.legend(loc=0)
			ax.set_xlabel('Number of iterations of BiCGSTAB solver')
			ax.set_ylabel('Relative residue ' + r'$\displaystyle\frac{||r||_\infty}{||b||_\infty}$')
			ax.set_ybound(lower=1e-12, upper=10)

		filename = folder + 'TransientNL_' + geoname + '.pdf'
		fig.tight_layout()
		fig.savefig(filename)