"""
In this file, we analyze the stability of the time-integration scheme
using IsoGeometric Analysis in a 1D problem. 
The Laplace problem is:
	sigma*dT/dt -div(k grad(T)) = f, 
	with T(x, 0) = 0,
	T(0, t) = 0 and T(1, t) = 1
"""

from lib.__init__ import *
from lib.base_functions import (create_knotvector,
								wq_find_basis_weights_fortran,
								eigen_decomposition,
								erase_rows_csr
)

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/stability/'
if not os.path.isdir(folder): os.mkdir(folder)

def build_sparse_matrix(basis, weights, indi_in, indj_in):
	indi = np.copy(indi_in); indj = np.copy(indj_in)
	indi -= 1; indj -= 1
	nrows = len(indi) - 1; ncols = np.max(indj)
	B0  = sp.csr_matrix((basis[:,0], indj, indi), shape=(nrows, ncols)).toarray()
	B1  = sp.csr_matrix((basis[:,1], indj, indi), shape=(nrows, ncols)).toarray()
	W00 = sp.csr_matrix((weights[:,0], indj, indi), shape=(nrows, ncols)).toarray()
	W11 = sp.csr_matrix((weights[:,3], indj, indi), shape=(nrows, ncols)).toarray()
	return B0, B1, W00, W11

def scheme_analysis(degree, cuts, prop):
	""" Here, prop = k/rho*Cp (isotropic material)
		First, we analyse time scheme stability : time step <= time_stab. 
		We use eigenvalues in such a way that :
			time_stab = 2/lambda_max (in an Euler explicit scheme). 
		We solve the generalized eigenvalue problem => prop K U = lambda M U
		Then, we analyse space oscillations : time step  >= time_osc.
		We build A = M + theta*(time step)*prop*K,
		All non-diagonal terms of A must be less or equal than 0
	"""
	nbel 	   = int(2**cuts)
	knotvector = create_knotvector(degree, nbel)

	# Get basis and weights in IgA 
	qp, basis_in, weights_in, indi_in, indj_in = wq_find_basis_weights_fortran(degree, knotvector)[1:]
	indi, indj, data = erase_rows_csr([0], indi_in, indj_in, [basis_in, weights_in])
	basis, weights   = data 
	data_B0   = basis[:, 0]  ; data_B1  = basis[:, 1]
	data_W00  = weights[:, 0]; data_W11 = weights[:, -1]
	data = [data_B0, data_B1, data_W00, data_W11]

	# Time scheme stability
	mcoefs = np.ones(len(qp)); kcoefs = prop*np.ones(len(qp)); coefs = [mcoefs, kcoefs]
	eigenvalues = eigen_decomposition(indi, indj, data, coefs=coefs)[0]
	lambda_max  = max(eigenvalues)
	t_stab      = 2/lambda_max

	# Space oscillations
	B0, B1, W00, W11 = build_sparse_matrix(basis, weights, indi, indj)
	Mass  = W00 @ B0.transpose()
	Stiff = -prop * W11 @ B1.transpose()
	newmark_matrix = np.divide(Mass, Stiff, out=np.zeros_like(Mass), where=Stiff!=0)
	t_osc = newmark_matrix.max()

	return t_stab, t_osc

def plot_analysis(degree_list, cuts_list, filenumber=0, folder=None, extension='.png'):
	fig, ax   = plt.subplots(nrows=1, ncols=1)
	nbel_list = [2**i for i in cuts_list]
	name      = 'time_stab_' + str(filenumber)
	# name      = 'time_osc_' + str(filenumber)
	time_plot = np.loadtxt(folder + name + '.dat')
	plot = ax.contourf(nbel_list, degree_list, time_plot, 
						norm=mpl.colors.LogNorm(vmin=time_plot.min(), vmax=time_plot.max()))
	plt.xscale('log')
	ax.grid(None)
	ax.set_ylabel('Polynomial degree '+ r'$p$')
	ax.set_xlabel('Mesh discretization level ' + r'$h$')
	fig.colorbar(plot)
	fig.tight_layout()
	fig.savefig(folder + name + extension)
	return

data_exist  = True
degree_list = range(2, 11)
cuts_list   = range(3, 10)
prop_list   = [1e-5, 1e-3, 0.1, 1, 10]

if not data_exist:

	for k, prop in enumerate(prop_list):
		save_data1 = np.zeros((len(degree_list), len(cuts_list)))
		save_data2 = np.zeros((len(degree_list), len(cuts_list)))
		for i, degree in enumerate(degree_list):
			for j, cuts in enumerate(cuts_list):
				t_stab, t_osc = scheme_analysis(degree, cuts, prop)
				save_data1[i, j] = t_stab
				save_data2[i, j] = t_osc
		name = 'time_stab_' + str(k) + '.dat'
		np.savetxt(folder + name, save_data1)
		name = 'time_osc_'  + str(k) + '.dat'
		np.savetxt(folder + name, save_data2)
else:
	for filenumber in range(len(prop_list)):
		plot_analysis(degree_list, cuts_list, filenumber=filenumber, folder=folder)

