"""
.. Isogeometric analysis figures
.. Select case: 
.. CASE 0: B-spline curve
.. CASE 1: Univariate B-spline functions in parametric space
.. CASE 2: Bivariate B-spline functions in parametric space
.. CASE 3: Quadrature points in IGA-Galerkin approach
.. CASE 4: Quadrature points in IGA-WQ approach
.. CASE 5: B-spline surface
.. CASE 6: FEM basis 
.. CASE 7: Convergence curve
.. CASE 8: Quadrature rules W00 or W11
.. CASE 9: Plot 3D geometries
.. CASE 10: Plot example of results
"""

from lib.__init__ import *
from lib.lib_base import (createKnotVector, evalDersBasisPy, cropImage)
from lib.lib_quadrules import *
from lib.lib_geomdl import Geomdl
from lib.lib_model import part
# from lib.fortran_mf_wq import fortran_mf_wq

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/phd/'
if not os.path.isdir(folder): os.mkdir(folder)

def plot2DGeo(model:part):
	"Plots a 2D geometry "

	def plot_mesh(pts, shape, ax):
		"Plots mesh of control points"

		if pts.shape[0] == 3: pts = pts[:2, :]

		pts2D = []
		for j in range(shape[1]):
			pts2D_temp = []
			for i in range(shape[0]):
				pos = i + j*shape[0]
				pts2D_temp.append(pts[:, pos].tolist())
			pts2D.append(pts2D_temp)
		pts2D = np.asarray(pts2D)

		# In the first direction
		for _ in range(shape[1]): 
			x = pts2D[_, :, 0]; y = pts2D[_, :, 1]
			ax.plot(x, y, 'k--')

		# In the second direction
		for _ in range(shape[0]):
			x = pts2D[:, _, 0]; y = pts2D[:, _, 1]
			ax.plot(x, y, 'k--')

		return

	samplesize = model._sampleSize
	ctrlpts = model._ctrlpts
	evalpts = model.interpolateField()[1]

	X = np.asarray(evalpts[0, :].reshape((samplesize, samplesize)).tolist())
	Y = np.asarray(evalpts[1, :].reshape((samplesize, samplesize)).tolist())
	Z = np.zeros(X.shape)

	fig, ax = plt.subplots()
	ax.grid(None)
	ax.pcolormesh(X, Y, Z, cmap=plt.cm.Pastel1, shading='gouraud')
	ax.plot(ctrlpts[0, :], ctrlpts[1, :], 'o', label='Control points')
	plot_mesh(ctrlpts, model._nbctrlpts, ax)
	ax.set_xticks(np.arange(0, max(evalpts[:,0])+2, 1.0))
	ax.set_yticks(np.arange(0, max(evalpts[:,1])+1, 1.0))

	ax.axis('equal')
	ax.legend()
	ax.set_xlabel(r'$X_1$')
	ax.set_ylabel(r'$X_2$')
	fig.tight_layout()
	
	return fig

def plotVerticalLine(x, y, ax=None, color='k'):
	for xi, yi in zip(x, y):
		ax.plot([xi, xi], [0, yi], color=color)
	return

# Set global variables
CASE      = 8
extension = '.png'

if CASE == 0: # B-spline curve

	def case0(folder, extension):
		# Set filename
		filename = folder + 'BsplineCurve' + extension

		# Create the curve 
		degree, nbel = 2, 4
		knotvector   = createKnotVector(degree, nbel)
		crv            = BSpline.Curve()
		crv.degree     = degree
		crv.ctrlpts    = [[-1, 1, 0], [-0.5, 0.25, 0], [0, 2, 0], 
							[0.75, -0.5, 0], [1.5, 1, 0], [2, 0, 0]]
		crv.knotvector = knotvector
		crv.delta      = 0.01

		# Get data
		evalpts = np.asarray(crv.evalpts)
		ctrlpts = np.asarray(crv.ctrlpts)
		x = evalpts[:, 0]; y = evalpts[:, 1]
		
		fig, ax = plt.subplots(nrows=1, ncols=1)
		ax.plot(x, y, label='B-spline curve')
		ax.plot(ctrlpts[:, 0], ctrlpts[:, 1], 'o--', markersize=10, label='Control points')

		ax.legend()
		ax.set_xlabel(r'$X_1$')
		ax.set_ylabel(r'$X_2$')
		ax.axis('equal')
		fig.tight_layout()
		fig.savefig(filename)
		return

	case0(folder, extension)
	
elif CASE == 1: # Univariate functions

	def case1(folder, extension): 
		# Set filename
		filename = folder + 'UnivariateFunctions' + extension

		# B-spline properties 
		degree, nbel = 3, 4
		multiplicity = 1
		knotvector   = createKnotVector(degree, nbel, multiplicity=multiplicity)
		quadRule     = QuadratureRules(degree, knotvector)
		basis, knots = quadRule.getGeneralizedBasis()
		B0 = basis[0].toarray(); B1 = basis[0].toarray()

		fig, [ax1, ax2] = plt.subplots(nrows=1, ncols=2, figsize=(10, 4))
		for i in range(np.shape(B0)[0]): 
			ax1.plot(knots, B0[i, :], linewidth=2)
			ax2.plot(knots, B1[i, :], linewidth=2)

		for ax in [ax1, ax2]:
			ax.set_xlabel(r'$\xi$')
			ax.set_xticks(np.linspace(0, 1, nbel+1))
		ax1.set_ylabel(r'$B_{i,p}(\xi)$')
		ax2.set_ylabel(r"${B'}_{i,p}(\xi)$")
		fig.tight_layout()
		fig.savefig(filename)
		return
	
	case1(folder, extension)

elif CASE == 2: # Bivariate functions

	def case2(folder, extension, is2D=True):
		# Set filename
		filename = folder + 'BivariateFunctions' + extension

		# B-Spline properties
		degree, nbel = 3, 4
		knotvector   = createKnotVector(degree, nbel)
		quadRule     = QuadratureRules(degree, knotvector)
		basis, knots = quadRule.getGeneralizedBasis()
		B0 = basis[0].toarray()
		B02plot = B0[1, :]

		# B-Spline 2D
		X, Y = np.meshgrid(knots, knots)
		Z = np.kron(B02plot, B02plot).reshape((len(knots), len(knots)))

		if is2D:
			fig, axs = plt.subplots(2, 2, sharex="col", sharey="row", 
									gridspec_kw=dict(height_ratios=[1,3],
													width_ratios=[3,1]))
			axs[0,1].set_visible(False)
			axs[0,0].set_box_aspect(1/3)
			axs[1,0].set_box_aspect(1)
			axs[1,1].set_box_aspect(3/1)
			axs[1,0].grid(None)
			axs[1,0].pcolormesh(X, Y, Z, cmap='GnBu', shading='gouraud', rasterized=True)
			axs[1,0].set_yticks([0, 0.5, 1])

			for i in range(degree+nbel): 
				axs[0, 0].plot(knots, B0[i, :], color="0.8")
				axs[1, 1].plot(B0[i, :], knots, color="0.8")

			axs[0,0].plot(knots, B02plot); axs[0, 0].axis(ymin=0, ymax=1)
			axs[0,0].set_xlabel(r'$\xi_1$')
			axs[0,0].set_ylabel(r'$B_{i_1,p_1}(\xi_1)$')
			axs[1,1].plot(B02plot, knots); axs[1, 1].axis(xmin=0, xmax=1)
			axs[1,1].set_ylabel(r'$\xi_2$')
			axs[1,1].set_xlabel(r'$B_{i_2,p_2}(\xi_2)$')
			axs[1,1].set_xticks([0, 1])
			fig.tight_layout()
			fig.savefig(filename, dpi=300) 

		else:
			fig = plt.figure()
			ax = plt.axes(projection='3d')
			ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap='GnBu', edgecolor='none')
			ax.grid(False)
			ax.set_xticks([0, 0.5, 1])
			ax.set_yticks([0, 0.5, 1])
			ax.set_zticks([0, 0.1, 0.2])
			ax.set_xlabel(r'$\xi_1$')
			ax.set_ylabel(r'$\xi_2$')
			ax.set_zlabel(r'$B_{i,p}(\xi)$')
			fig.tight_layout()
			fig.savefig(filename) 
		return
	
	case2(folder, extension)

elif CASE == 3: # Quadrature points in IGA

	def case3(folder, extension):
		# Set filename
		filename = folder + 'QuadPtsIGA' + extension

		fig, [ax1, ax2] = plt.subplots(nrows=1, ncols=2, figsize=(8, 4))
		degree = 4
		for ax, nbel in zip([ax1, ax2], [8, 32]):
			knotvector  = createKnotVector(degree, nbel)
			quadRule    = GaussQuadrature(degree, knotvector)
			quadRule.getQuadratureRulesInfo()
			XX, YY      = np.meshgrid(quadRule._quadPtsPos, quadRule._quadPtsPos)
			ax.plot(XX, YY, 'ko', markersize=1.2)

			grid = np.linspace(0, 1, nbel+1)
			for i in grid:
				ax.plot([i, i], [0, 1], 'grey', linewidth=0.5, alpha=0.8)
				ax.plot([0, 1], [i, i], 'grey', linewidth=0.5, alpha=0.8)

			ax.set_xticks([0, 0.5, 1])
			ax.set_yticks([0, 0.5, 1])
			ax.axis('equal')
			ax.set_ylabel(r'$\xi_2$')
			ax.set_xlabel(r'$\xi_1$')   

		fig.tight_layout()
		fig.savefig(filename) 
		return

	case3(folder, extension)

elif CASE == 4: # Quadrature points in WQ

	def case4(folder, extension):
		# Set filename
		filename = folder + 'QuadPtsWQ' + extension

		fig, [ax1, ax2] = plt.subplots(nrows=1, ncols=2, figsize=(8, 4))
		degree = 4
		for ax, nbel in zip([ax1, ax2], [8, 32]):
			knotvector  = createKnotVector(degree, nbel)
			quadRule    = WeightedQuadrature(degree, knotvector)
			quadRule.getQuadratureRulesInfo()
			XX, YY      = np.meshgrid(quadRule._quadPtsPos, quadRule._quadPtsPos)
			ax.plot(XX, YY, 'ko', markersize=1.2)

			grid = np.linspace(0.,1,nbel+1)
			for i in grid:
				ax.plot([i, i], [0, 1], 'grey', linewidth=0.5, alpha=0.8)
				ax.plot([0, 1], [i, i], 'grey', linewidth=0.5, alpha=0.8)

			ax.set_xticks([0, 0.5, 1])
			ax.set_yticks([0, 0.5, 1])
			ax.axis('equal')
			ax.set_ylabel(r'$\xi_2$')
			ax.set_xlabel(r'$\xi_1$')

		fig.tight_layout()
		fig.savefig(filename) 
		return
	
	case4(folder, extension)

elif CASE == 5: # B-spline surface

	# Set filename
	filename = folder + 'BsplineSurface' + extension

	# Surface properties
	modelGeo = Geomdl(**{'name':'quarter_annulus'})
	modelIGA = modelGeo.getIGAParametrization()
	model    = part(modelIGA)
	fig = plot2DGeo(model)
	fig.savefig(filename) 

elif CASE == 6: # FEM functions 

	def case6(folder, extension):
		# Set filename
		filename = folder + 'FEM_Functions' + extension

		# Functions in isoparametric space 
		x  = np.linspace(-1, 1, 90)
		f1 = x*(x-1)/2
		f2 = (1-x**2)
		f3 = x*(x+1)/2

		# Parametric space properties
		nbel  = 4
		knots = np.linspace(0, 1, nbel+1)

		fig, ax = plt.subplots()
		for _ in range(nbel):
			x0 = knots[_]; x1 = knots[_+1]
			xtemp = x0*(1-x)/2 + x1*(1+x)/2
			ax.plot(xtemp, f1, color=colorSet[_])
			ax.plot(xtemp, f2, color=colorSet[_])
			ax.plot(xtemp, f3, color=colorSet[_])
		
		ax.set_xlabel(r'$\xi$')
		ax.set_xticks([0, 0.5, 1])
		fig.tight_layout()
		fig.savefig(filename)
		return
	
	case6(folder, extension)

# elif CASE == 7: # Convergence curve

# 	# Set filename
# 	filename = folder + 'ConvergenceIGA'+ extension

# 	def powden(P:list):
# 		x, y = P
# 		f = np.sin(np.pi*x)*np.sin(np.pi*y)
# 		return f

# 	def solution(P:list): 
# 		x, y = P
# 		t = 1/(2*np.pi**2)*np.sin(np.pi*x)*np.sin(np.pi*y)
# 		return t

# 	fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 4))
# 	for degree in range(2, 8):
# 		norm = []; nbel_list =[]
# 		for cuts in range(1, 7):
# 			print([degree, 2**cuts])

# 			blockPrint()
# 			# Create model
# 			name     = 'quadrilateral'
# 			modelGeo = geomdlModel(name=name, **{'degree':[degree, degree, degree]})
# 			modelGeo.knot_refinement(np.array([cuts, cuts, cuts]))
# 			modelPhy = fortran_mf_wq(modelGeo)

# 			# Add material 
# 			material = {'capacity':1, 'conductivity':np.eye(2)}
# 			modelPhy._set_material(material)

# 			# Block boundaries
# 			Dirichlet = {'thermal':np.array([[1, 1], [1, 1], [1, 1]])}
# 			modelPhy._set_dirichlet_boundaries(Dirichlet)
# 			dof = modelPhy._thermal_dof
# 			dod = modelPhy._thermal_dod 

# 			# Solve
# 			Kdd = modelPhy.eval_conductivity_matrix(indi=dof, indj=dof)
# 			Fd  = modelPhy.eval_source_vector(powden, indi=dof)
# 			Td  = sp.linalg.spsolve(Kdd, Fd)  
# 			T   = np.zeros(modelPhy._nb_ctrlpts_total); T[dof] = Td
# 			enablePrint()

# 			# Interpolate
# 			output  = modelPhy.interpolate_field(u_ctrlpts=T, nbDOF=1)
# 			qp_interp, u_interp = output[1], output[-1]
# 			u_exact = [solution(qp_interp[:, i]) for i in range(len(u_interp))]
# 			u_exact = np.array(u_exact)

# 			# Relative error
# 			error = relativeError(u_interp, u_exact)
# 			norm.append(error)
# 			nbel_list.append(2**cuts)
			
# 		ax.loglog(nbel_list, norm, label='Degree $p =$ ' + str(degree))
# 		slope = np.polyfit(np.log10(nbel_list[1:5]),np.log10(norm[1:5]), 1)[0]
# 		slope = round(slope)
# 		annotation.slope_marker((nbel_list[3], norm[3]), slope, 
# 								poly_kwargs={'facecolor': (0.73, 0.8, 1)})

# 	ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
# 	ax.set_xlabel('Discretization level ' + r'$h^{-1}$')
# 	ax.set_ylabel('Relative error ' + r'$\displaystyle\frac{||u-u^h||_\infty}{||u||_\infty}$')   
# 	fig.tight_layout()
# 	fig.savefig(filename)

elif CASE == 8: # Weights W00 and W11

	WeightName = 1 # or 0
	if WeightName not in [0, 1]: raise Warning('Not possible')
	if WeightName == 0: ylim1 = [0, 1];  ylim2 = [0, 0.25]
	else:               ylim1 = [-3, 3]; ylim2 = [-0.6, 0.6]

	# Set filename
	filename = folder + 'WeightsW' + str(WeightName) + extension

	# B-spline properties 
	WeightPos    = 2
	degree, nbel = 2, 3
	knotvector   = createKnotVector(degree, nbel)
	quadRule     = WeightedQuadrature(degree, knotvector)
	quadRule.getQuadratureRulesInfo()
	basis, knots = quadRule.getGeneralizedBasis()
	B = basis[WeightName].toarray()

	# Get weights
	quadRule.getDenseQuadRules()
	W00 = quadRule._denseWeights[-WeightName]
	weights = W00.toarray()[WeightPos, :]

	if WeightName == 0:
		Bref = B
	else:
		kvref = createKnotVector(degree-1, nbel)
		Bref = evalDersBasisPy(degree-1, kvref, knots)[0]
		Bref = Bref.toarray()

	fig, ax1 = plt.subplots(nrows=1, ncols=1, figsize=(8,4))
	for i in range(np.shape(Bref)[0]): 
		ax1.plot(knots, Bref[i, :], linewidth=2)

	# Fill basis chosen
	ax1.fill_between(x=knots, y1=B[WeightPos, :], color="g", alpha=0.2)
	ax1.set_ylim(ylim1)

	ax2 = ax1.twinx()
	ax2.plot(quadRule._quadPtsPos, weights, 'ko')
	plotVerticalLine(quadRule._quadPtsPos, weights, ax2)
	ax2.set_ylim(ylim2)
	ax2.grid(None)

	ax1.set_xlabel(r'$\xi$')
	ax1.set_ylabel('Basis')
	ax2.set_ylabel('Weights')
	ax1.set_xticks(np.linspace(0, 1, nbel+1), ['0', '1/3', '2/3', '1'])
	fig.tight_layout()
	fig.savefig(filename)

# elif CASE == 9: # 3D Geometries

# 	geoName   = 'TR'
# 	dataExist = True

# 	# Set filename
# 	filename = folder + 'VTK_' + geoName + '.png'

# 	if not dataExist:

# 		# Create model
# 		modelGeo = geomdlModel(geoName)
# 		modelGeo.knot_refinement(nb_refinementByDirection=np.array([1, 1, 1]))
# 		modelPhy = fortran_mf_wq(modelGeo)
# 		modelPhy.export_results()    

# 		# Read data
# 		full_path = os.path.realpath(__file__)
# 		fileVTK   = os.path.dirname(full_path) + '/results/' + geoName
# 		grid      = pv.read(fileVTK + '.vts')
		
# 		sargs = dict(
# 				title = 'det J',
# 				title_font_size=50,
# 				label_font_size=40,
# 				shadow=True,
# 				n_labels=2,
# 				fmt="%.1f",
# 				position_x=0.2, 
# 				position_y=0.1,
# 		)
# 		pv.start_xvfb()
# 		plotter = pv.Plotter(off_screen=True)
# 		if geoName == 'CB': 
# 			boring_cmap = plt.cm.get_cmap("GnBu", 1)
# 			plotter.add_mesh(grid, cmap=boring_cmap, scalar_bar_args=sargs)
# 			plotter.camera.zoom(0.6)
# 		else: 
# 			plotter.add_mesh(grid, cmap='viridis', scalar_bar_args=sargs)

# 		if geoName == 'VB': 
# 			plotter.camera_position  = 'yz'
# 			plotter.camera.elevation = 45
# 		elif geoName == 'RQA': 
# 			plotter.camera_position  = 'xz'
# 			plotter.camera.azimuth   = -45
# 			plotter.camera.zoom(0.8)

# 		plotter.background_color = 'white'
# 		plotter.window_size = [1600, 1600]
# 		plotter.screenshot(filename)
# 		os.remove(fileVTK + '.vts')

# 	cropImage(filename)

# elif CASE == 10: # Results

# 	# Read data
# 	full_path = os.path.realpath(__file__)
# 	folderVTK = os.path.dirname(full_path) + '/results/test8/'
# 	grid = pv.read(folderVTK + 'IGAparametrization.vts')
# 	filename = folder + 'VTK_Results' + '.png'

# 	sargs = dict(
# 			title = 'Temperature (K)',
# 			title_font_size=50,
# 			label_font_size=40,
# 			shadow=True,
# 			n_labels=2,
# 			fmt="%.1f",
# 			position_x=0.2, 
# 			position_y=0.1,
# 	)

# 	pv.start_xvfb()
# 	plotter = pv.Plotter(off_screen=True)
# 	plotter.add_mesh(grid, cmap='viridis', scalar_bar_args=sargs)
# 	plotter.camera.zoom(0.6)
# 	plotter.background_color = 'white'
# 	plotter.window_size = [1600, 1600]
# 	plotter.screenshot(filename)

# 	cropImage(filename)

else: raise Warning('Case unkwnon')