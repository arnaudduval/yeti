"""
.. Test of elastoplasticity 1D
.. Joaquin Cornejo 
"""

import pickle
from pysrc.lib.__init__ import *
from pysrc.lib.lib_base import createUniformKnotvector_Rmultiplicity
from pysrc.lib.lib_quadrules import QuadratureRules
from pysrc.lib.lib_1d import mechamat1D

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/d1elastoplasticity/'
if not os.path.isdir(folder): os.mkdir(folder)

def plot_results(quadRule:QuadratureRules, JJ, disp_cp, plastic_cp, stress_cp, folder=None, method='iga', extension='.png'):
	from mpl_toolkits.axes_grid1 import make_axes_locatable
	basis, knots = quadRule.getSampleBasis(sampleSize=101)
	displacement   = basis[0].T @ disp_cp
	strain_interp  = basis[1].T @ disp_cp / JJ
	plastic_interp = basis[0].T @ plastic_cp
	stress_interp  = basis[0].T @ stress_cp

	# Plot fields
	nbsteps   = np.shape(disp_cp)[1]
	XX, STEPS = np.meshgrid(knots*JJ, np.arange(nbsteps))
	names = ['Displacement field', 'Plastic strain field', 'Stress field']
	units = ['mm', '', 'MPa']
	fig, [ax1, ax2, ax3] = plt.subplots(nrows=1, ncols=3, figsize=(16, 4))
	for ax, variable, name, unit in zip([ax1, ax2, ax3], [displacement, plastic_interp, stress_interp], names, units):
		im = ax.pcolormesh(XX, STEPS, variable.T, cmap='PuBu_r', shading='linear')
		ax.set_title(name)
		ax.set_ylabel('Step')
		ax.set_xlabel('Position (mm)')
		ax.grid(False)
		divider = make_axes_locatable(ax)
		cax = divider.append_axes('right', size='5%', pad=0.05)
		cbar = fig.colorbar(im, cax=cax)
		cbar.ax.set_title(unit)

	fig.tight_layout()
	fig.savefig(folder + 'ElastoPlasticity' + method + extension)

	# Plot stress-strain of single point
	fig, [ax1, ax2, ax3] = plt.subplots(nrows=1, ncols=3, figsize=(14,4))
	for ax, pos in zip([ax1, ax2, ax3], [25, 50, 75]):
		ax.plot(strain_interp[pos, :], stress_interp[pos, :])
		ax.set_ylabel('Stress (MPa)')
		ax.set_xlabel('Total strain (-)')
		ax.set_ylim(bottom=0.0, top=200)
		ax.set_xlim(left=0.0, right=strain_interp.max())

	fig.tight_layout()
	fig.savefig(folder + 'TractionCurve' + method + extension)
	return


def forceVol(P:list):
	force = 0.4*np.sin(P/1e3)
	return force

# Set global variables
nbsteps = 50
geoArgs = {'length': 1.e3}
matArgs = {'elastic_modulus':2e5, 'elastic_limit':100,
			'plasticLaw': {'Isoname': 'swift', 'e0':2e4, 'n':0.5}} #!!!!!!!!!!!!!!!!!!

degree, nbel = 8, 1024
knotvector   = createUniformKnotvector_Rmultiplicity(degree, nbel, multiplicity=degree)

# Create geometry
quadArgs  = {'degree': degree, 'knotvector': knotvector, 'quadrule': 'iga', 'type': 'leg'}
args  = {'quadArgs': quadArgs, 'geoArgs': geoArgs}
modelPhy = mechamat1D(args)

with open(folder + 'refpart.pkl', 'wb') as outp:
    pickle.dump(modelPhy, outp, pickle.HIGHEST_PROTOCOL)

# Add material
modelPhy.activate_mechanical(matArgs)

# Add boundary condition
modelPhy.add_DirichletCondition(table=[1, 0])

# Define boundaries conditions
Fext    = np.zeros((modelPhy.nbctrlpts, 2*nbsteps + 1))
Fextref = modelPhy.compute_volForce(forceVol(modelPhy.qpPhy))
for i in range(0, nbsteps+1): Fext[:, i] = i/nbsteps*Fextref
for i in range(nbsteps+1, 2*nbsteps+1): Fext[:, i] = (2*nbsteps - i)/nbsteps*Fextref

# Solve
disp_cp, strain_qp, stress_qp, plastic_qp, Cep_qp = modelPhy.solve(Fext=Fext)
np.save(folder+'disp', disp_cp)
plastic_cp  = modelPhy.L2projectionCtrlptsVol(plastic_qp)
stress_cp 	= modelPhy.L2projectionCtrlptsVol(stress_qp)
plot_results(modelPhy.quadRule, geoArgs['length'], disp_cp, plastic_cp, stress_cp, folder=folder, method=quadArgs['quadrule'])