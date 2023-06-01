from lib.__init__ import *
from lib.lib_base import (createKnotVector)
from lib.thermomecha1D import mechamat1D, plot_results
from lib.lib_load import *

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/d1plasticity/'
if not os.path.isdir(folder): os.mkdir(folder)

# Set global variables
matArgs   = {'elastic_modulus':200e3, 'elastic_limit':506, 'plasticLaw': {'name': 'swift', 'K':2e4, 'exp':0.5}}
length    = 1.0
nbSteps   = 101
geoArgs   = {'length': length}
sampleSize = 2500
lastStep  = 53

# ---------------
# IGA
# ---------------
degree, nbel = 7, 16
knotvector   = createKnotVector(degree, nbel)
quadArgs  = {'degree': degree, 'knotvector': knotvector, 'quadrule': 'iga', 'type': 'leg'}
args  = {'quadArgs': quadArgs, 'geoArgs': geoArgs}
model = mechamat1D(args)

# Add material
model.activate_mechanical(matArgs)

# Add boundary condition
model.add_DirichletCondition(table=[1, 0])
nbqpiga = model.nbqp
print('nb quad points: %s' %str(nbqpiga))

# Define boundaries conditions
Fext        = np.zeros((model.nbctrlpts, nbSteps))
Fext[:, -1] = model.compute_volForce(forceVol(model.qpPhy))
for i in range(1, nbSteps-1): Fext[:, i] = i/(nbSteps-1)*Fext[:, -1]

# blockPrint()
# Solve 
disp_cp, strain_iga, stress_iga, plastic_iga, Cep_iga = model.solve(Fext=Fext[:, :lastStep])
strain_cp   = model.interpolate_CntrlPtsField(strain_iga)
plastic_cp  = model.interpolate_CntrlPtsField(plastic_iga)
stress_cp 	= model.interpolate_CntrlPtsField(stress_iga)
# plot_results(model.quadRule, length, disp_cp, plastic_cp, stress_cp, folder=folder, method='IGA')

interp_iga  = []
interp_iga.append(model.interpolate_sampleField(disp_cp, sampleSize=sampleSize)[0])
interp_iga.append(model.interpolate_sampleField(strain_cp, sampleSize=sampleSize)[0])
interp_iga.append(model.interpolate_sampleField(stress_cp, sampleSize=sampleSize)[0])
# enablePrint()

# -------------------
# Weighted Quadrature
# -------------------
# degree, nbel = 6, 442
knotvector   = createKnotVector(degree, nbel)

quadArgs  = {'degree': degree, 'knotvector': knotvector, 'quadrule': 'wq', 'type': 1}
# quadArgs  = {'degree': degree, 'knotvector': knotvector, 'quadrule': 'iga', 'type': 'lob'}
args  = {'quadArgs': quadArgs, 'geoArgs': geoArgs}
model = mechamat1D(args)

# Add material
model.activate_mechanical(matArgs)

# Add boundary condition
model.add_DirichletCondition(table=[1, 0])
nbqpwq = model.nbqp
print('nb quad points: %s' %str(nbqpwq))

# Define boundaries conditions
Fext        = np.zeros((model.nbctrlpts, nbSteps))
Fext[:, -1] = model.compute_volForce(forceVol(model.qpPhy))
for i in range(1, nbSteps-1): Fext[:, i] = i/(nbSteps-1)*Fext[:, -1]

# blockPrint()
# Solve
disp_cp, strain_wq, stress_wq, plastic_wq, Cep_wq = model.solve(Fext=Fext[:, :lastStep])
strain_cp   = model.interpolate_CntrlPtsField(strain_wq)
plastic_cp  = model.interpolate_CntrlPtsField(plastic_wq)
stress_cp 	= model.interpolate_CntrlPtsField(stress_wq)
# plot_results(model.quadRule, length, disp_cp, plastic_cp, stress_cp, folder=folder, method='WQ')

interp_wq   = []
interp_wq.append(model.interpolate_sampleField(disp_cp, sampleSize=sampleSize)[0])
interp_wq.append(model.interpolate_sampleField(strain_cp, sampleSize=sampleSize)[0])
interp_wq.append(model.interpolate_sampleField(stress_cp, sampleSize=sampleSize)[0])
# enablePrint()

tmp = interp_iga[-1][:, -1] - interp_wq[-1][:, -1]
relerror = np.linalg.norm(tmp)/np.linalg.norm(interp_iga[-1][:, -1])*100
print('Relative error: %.3e' %relerror)

# # ------------------
# # Post-treatement
# # ------------------
# # fig, [ax0, ax1] = plt.subplots(nrows=1, ncols=2, figsize=(10, 4))
# # ax0.plot(model.qpPhy, stress_wq[:, 51])
# # ax0.set_ylabel('Stress (MPa)')
# # ax1.plot(model.qpPhy, Cep_wq[:, 51])
# # ax1.set_ylabel('Tangent modulus (MPa)')
# # for ax in [ax0, ax1]:
# # 	ax.set_ylim(bottom=0.0)
# # 	ax.set_xlim(left=0.0, right=1.0)
# # 	ax.set_xlabel('Quadrature point position')
# # fig.tight_layout()
# # fig.savefig(folder+'data_step51')

# ref = []
# ref.append(np.load(folder + 'disp_interp_ref.npy'))
# ref.append(np.load(folder + 'strain_interp_ref.npy'))
# ref.append(np.load(folder + 'stress_interp_ref.npy'))
# error_iga = np.zeros((sampleSize, nbSteps, 3))
# error_wq  = np.zeros((sampleSize, nbSteps, 3))

# for i in range(3):
# 	error_iga[:, :, i] = ref[i] - interp_iga[i]
# 	error_wq[:, :, i]  = ref[i] - interp_wq[i]

# fig, axs  = plt.subplots(nrows=1, ncols=3, figsize=(14, 4))
# for [i, ax], title in zip(enumerate(axs), ['Displacement', 'Strain', 'Stress']):
	
# 	norm_ref     = np.linalg.norm(ref[i], axis=0)
	
# 	relerror_iga = np.linalg.norm(error_iga[:, :, i], axis=0)
# 	relerror_iga = np.divide(relerror_iga, norm_ref, out=np.zeros_like(relerror_iga), where=np.abs(norm_ref)>1.e-12)*100
# 	ax.semilogy(relerror_iga[1:],label='iga-leg w.r.t. ref.')

# 	relerror_wq  = np.linalg.norm(error_wq[:, :, i], axis=0)
# 	relerror_wq  = np.divide(relerror_wq, norm_ref, out=np.zeros_like(relerror_wq), where=np.abs(norm_ref)>1.e-12)*100
# 	# ax.semilogy(relerror_wq[1:], label='iga-lob w.r.t. ref.')
# 	ax.semilogy(relerror_wq[1:], label='WQ w.r.t. ref.')

# 	ax.set_title(title)
# 	ax.set_ylim(bottom=1e-12, top=1e0)
# 	ax.set_ylabel('L2 Relative error (\%)')
# 	ax.set_xlabel('Step')
# 	ax.legend()

# fig.tight_layout()
# fig.savefig(folder + 'Relative_error.png')
