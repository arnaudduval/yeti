"""
.. Test of transient heat solver
.. ATTENTION: IT ONLY WORKS WITH 'ISOTROPIC' MATERIALS
.. Joaquin Cornejo 
"""

from lib.__init__ import *
from lib.physics import *
from lib.create_geomdl import geomdlModel
from lib.fortran_mf_wq import fortran_mf_wq
from lib.base_functions import create_table_properties, sigmoid, cropImage

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/test7/'
if not os.path.isdir(folder): os.mkdir(folder)

# Set global variables
dataExist   = True
plotVTK	    = False
# geolist     = ['VB', 'CB']
geolist     = ['VB']
method_list = ['WP', 'C', 'JMC']
# method_list = ['JMC']

if not dataExist:

	degree, cuts = 6, 5
	conductivity, capacity = 0.1, 1.0
	theta = 1.0
	time_list = np.linspace(0, 20, 41)  
	table_Kprop = create_table_properties(setKprop, prop=conductivity)
	table_Cprop = create_table_properties(setCprop, prop=capacity)     

	for geoName in geolist:
		for PCGmethod in method_list:
			filename = folder + 'ResPCG_' + geoName + '_' + PCGmethod + '.dat'        

			# Create model 
			geometry = {'degree':[degree, degree, degree]}
			modelGeo = geomdlModel(geoName, **geometry)
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
			N = 40
			Tsol, resPCG, properties_first, properties_last = modelPhy.MFtransientHeatNL(Fext=Fext[:, :N], G=GBound[:, :N], time_list=time_list[:N],
											table_Kprop=table_Kprop, table_Cprop=table_Cprop, 
											methodPCG=PCGmethod, theta=theta)

			# folderVTK = folder + geoName + '_' + PCGmethod + '/'					
			# modelPhy.export_results(u_ctrlpts=properties_first, nbDOF=2, folder=folderVTK, name='Iga00')
			# properties_first[[0, 1], :] = properties_first[[1, 0], :]
			# modelPhy.export_results(u_ctrlpts=properties_first, nbDOF=2, folder=folderVTK, name='Iga01')

			# modelPhy.export_results(u_ctrlpts=properties_last, nbDOF=2, folder=folderVTK, name='Iga10')
			# properties_last[[0, 1], :] = properties_last[[1, 0], :]
			# modelPhy.export_results(u_ctrlpts=properties_last, nbDOF=2, folder=folderVTK, name='Iga11')
			
			print('Finish')
			np.savetxt(filename, resPCG)
			# modelPhy.export_results(u_ctrlpts=Tsol[:, -1], folder=folder, nbDOF=1)

else:

	for geoName in geolist:
		fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(6, 4.7))

		for i, PCGmethod in enumerate(method_list):
			filename = folder + 'ResPCG_' + geoName + '_' + PCGmethod + '.dat'
			resPCG = np.loadtxt(filename)

			if PCGmethod   == "WP" : labelmethod = 'w.o. \npreconditioner'
			elif PCGmethod == "C"  : labelmethod = 'Classic FD \nmethod'
			elif PCGmethod == "JMC": labelmethod = 'This work'

			# Print the first 
			step = resPCG[0, 0]; iterNL = resPCG[1, 0]
			newresidue = resPCG[2:, 0]; newresidue = newresidue[newresidue>0]
			ax.semilogy(np.arange(len(newresidue)), newresidue, '-', linewidth=2.5, marker=markerSet[i],
						label=labelmethod)
			
			# # Print the last
			# step = resPCG[0, -1]; iterNL = resPCG[1, -1]
			# newresidue = resPCG[2:, -1]; newresidue = newresidue[newresidue>0]
			# ax.semilogy(np.arange(len(newresidue)), newresidue, '-', linewidth=2.5, marker=markerSet[i],
			# 			label=labelmethod)


		ax.legend(bbox_to_anchor=(-0.25, 1.02, 1.25, 0.2), loc="lower left",
                mode="expand", ncol=3)
		ax.set_xlabel('Number of iterations of BiCGSTAB solver')
		ax.set_ylabel('Relative residue ' + r'$\displaystyle\frac{||r||_\infty}{||b||_\infty}$')
		ax.set_ybound(lower=1e-12, upper=10)

		filename = folder + 'TransientNL_' + geoName + '.pdf'
		# fig.tight_layout()
		fig.savefig(filename)

if plotVTK:

	for geoName in geolist:
		PCGmethod = 'JMC'
		folderVTK = folder + geoName + '_' + PCGmethod + '/'

		for i in [0, 1]: # first or last
			for j in [0, 1]:
				if j == 0: propName = 'Capacity'
				if j == 1: propName = 'Conductivity'
				gridName = 'Iga' + str(i) + str(j) + '.vts'
				grid = pv.read(folderVTK + gridName)

				sargs = dict(
					title = propName,
					title_font_size=50,
					label_font_size=40,
					shadow=True,
					n_labels=2,
					fmt="%.1f",
					position_x=0.2, 
					position_y=0.1,
				)
				pv.start_xvfb()
				plotter = pv.Plotter(off_screen=True)
				plotter.add_mesh(grid, cmap='viridis', scalar_bar_args=sargs)
				if geoName == 'CB': plotter.camera.zoom(0.6)				
				if geoName == 'VB': 
					plotter.camera_position  = 'yz'
					plotter.camera.elevation = 45

				plotter.background_color = 'white'
				plotter.window_size = [1600, 1600]
				filename = folderVTK + 'VTK_Results' + str(i) + str(j) + '.png'
				plotter.screenshot(filename)

				cropImage(filename)