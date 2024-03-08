"""
In this file, we analyze the stability of the time-integration scheme
using IsoGeometric Analysis in a 1D problem. 
The Laplace problem is:
	sigma*dT/dt -div(k grad(T)) = f, 
	with T(x, 0) = 0,
	T(0, t) = 0 and T(1, t) = 1
"""

from pysrc.lib.__init__ import *
from pysrc.lib.lib_base import createUniformKnotvector_Rmultiplicity, eraseRowsCSR, array2csr_matrix
from pysrc.lib.lib_quadrules import GaussQuadrature
from scipy import interpolate

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/stability/'
if not os.path.isdir(folder): os.mkdir(folder)

def build_sparse_matrix(basis, indi_in, indj_in):
	B0  = array2csr_matrix(basis[:, 0], indi_in, indj_in, isfortran=True)
	B1  = array2csr_matrix(basis[:, -1], indi_in, indj_in, isfortran=True)
	return B0.toarray(), B1.toarray()

def scheme_analysis(prop, degree, cuts=None, nbel=None):
	""" Here, prop = k/rho*Cp (isotropic material)
		First, we analyse time scheme stability : time step <= time_stab. 
		We use eigenvalues in such a way that :
			time_stab = 2/lambda_max (in an Euler explicit scheme). 
		We solve the generalized eigenvalue problem => prop K U = lambda M U
		Then, we analyse space oscillations : time step  >= time_osc.
		We build A = M + theta*(time step)*prop*K,
		All non-diagonal terms of A must be less or equal than 0
	"""

	if cuts is not None: nbel = int(2**cuts)
	knotvector = createUniformKnotvector_Rmultiplicity(degree, nbel, multiplicity=degree)

	gaussQuad = GaussQuadrature(degree, knotvector, {})
	info = gaussQuad.getQuadratureRulesInfo()
	qp, [indi_in, indj_in], basis_in, weights_in = info

	# Get basis and weights in IgA 
	indi, indj, [basis, weights] = eraseRowsCSR([0, -1], indi_in, indj_in, [basis_in, weights_in])

	# Time scheme stability
	mcoefs = np.ones(len(qp)); kcoefs = prop*np.ones(len(qp))
	eigenvalues, eigenvectors = geophy.eigen_decomposition_sp(mcoefs, kcoefs, indi, indj, basis, weights)
	lambda_max  = max(eigenvalues)
	t_stab      = 2/lambda_max

	# Space oscillations
	dof = np.arange(1, len(indi_in)-2)
	B0, B1 = build_sparse_matrix(basis_in, indi_in, indj_in)
	W0, W1 = build_sparse_matrix(weights_in, indi_in, indj_in)
	Mass  = W0 @ B0.transpose()
	Stiff = -prop * W1 @ B1.transpose()

	Mass  = Mass[dof, :][:, :]; Stiff = Stiff[dof, :][:, :]
	newmark_matrix = np.divide(Mass, Stiff, out=np.zeros_like(Mass), where=np.abs(Stiff)>1.e-12)
	t_osc = newmark_matrix.max()

	return t_stab, t_osc

def plot_analysis(degree_list, cuts_list=None, nbel_list=None, filenumber=0, folder=None, extension='.png', option=1):
	fig, ax = plt.subplots()

	if cuts_list is not None: nbel_list = [2**i for i in cuts_list]
	if option == 1: name = 'time_stab_' + str(filenumber); title = 'Maximum time step'
	if option == 2: name = 'time_osc_' + str(filenumber);  title = 'Minimum time step'
	time_plot = np.loadtxt(folder + name + '.dat')
	nbel_new, degree_new = np.meshgrid(nbel_list, degree_list)
	plot = ax.pcolormesh(nbel_new, degree_new, time_plot,
						norm=mpl.colors.LogNorm(vmin=time_plot.min(), vmax=time_plot.max()),
						cmap='viridis', shading='auto')
	cbar = fig.colorbar(plot, format='%.1e')
	cbar.ax.set_ylabel(title)
	ax.set_xscale('log')
	ax.set_ylabel('Polynomial degree '+ r'$p$')
	ax.set_xlabel('Number of elements')
	fig.tight_layout()
	fig.savefig(folder + name + extension)
	return

def plot_variable_degree(nbel, degree_list, cuts_list=None, nbel_list=None, filenumber=0, folder=None, extension='.png'):
	
	fig, ax = plt.subplots(nrows=1, ncols=1)

	if cuts_list is not None: nbel_list = [2**i for i in cuts_list]
	name = 'time_osc_' + str(filenumber);  title = 'Minimum time step'
	time_plot = np.loadtxt(folder + name + '.dat')

	f = interpolate.interp2d(degree_list, nbel_list, time_plot.T, kind='linear')
	xnew = degree_list
	ynew = np.array([nbel])
	znew = f(xnew, ynew)

	ax.semilogy(xnew, znew)
	ax.set_xlabel('Polynomial degree '+ r'$p$')
	ax.set_ylabel(title)
	ax.set_ylim(bottom=1e-3, top=1e-1)
	fig.tight_layout()
	fig.savefig(folder + name + '_' + str(nbel) + '_' + extension)
	return

dataExist   = True
degree_list = range(1, 7)
nbel_list   =    [i for i in range(4, 8)]
nbel_list.extend([i for i in range(8, 32, 2)])
nbel_list.extend([i for i in range(32, 128, 16)])
prop_list   = [1., 1e-2, 1e-4, 1e-6]

if not dataExist:

	for k, prop in enumerate(prop_list):
		save_data1 = np.zeros((len(degree_list), len(nbel_list)))
		save_data2 = np.zeros((len(degree_list), len(nbel_list)))
		for i, degree in enumerate(degree_list):
			for j, nbel in enumerate(nbel_list):
				t_stab, t_osc = scheme_analysis(prop, degree, nbel=nbel)
				save_data1[i, j] = t_stab
				save_data2[i, j] = t_osc
		name = 'time_stab_' + str(k) + '.dat'
		np.savetxt(folder + name, save_data1)
		name = 'time_osc_'  + str(k) + '.dat'
		np.savetxt(folder + name, save_data2)
else:
	for filenumber in range(len(prop_list)):
		plot_analysis(degree_list, nbel_list=nbel_list, filenumber=filenumber, folder=folder, option=1)
		plot_analysis(degree_list, nbel_list=nbel_list, filenumber=filenumber, folder=folder, option=2)

	plot_variable_degree(10, degree_list, nbel_list=nbel_list, filenumber=0, folder=folder,)
