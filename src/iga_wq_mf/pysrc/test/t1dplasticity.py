from pysrc.lib.__init__ import *
from pysrc.lib.lib_base import createUniformMaxregularKnotvector
from pysrc.lib.thermomecha1D import mechamat1D
import pickle

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/d1elastoplasticity/'
if not os.path.isdir(folder): os.mkdir(folder)

def forceVol(P:list):
	force = 0.4*np.sin(P/1e3)
	return force

# Set global variables
sampleSize = 2500
nbSteps   = 50
geoArgs   = {'length': 1.e3}
matArgs   = {'elastic_modulus':2e5, 'elastic_limit':100, 
			'plasticLaw': {'name': 'swift', 'K':2e4, 'exp':0.5}}

FigPlot = 0
disp_ref = np.load(folder + 'disp.npy')
with open(folder + 'refpart.pkl', 'rb') as inp:
    refPart = pickle.load(inp)

if FigPlot == 0:
	degree_list = np.arange(1, 5)
	cuts_list   = np.arange(3, 9)
	error = np.zeros(len(cuts_list))

	# First curve
	fig, ax = plt.subplots()
	for degree in degree_list:
		for j, cuts in enumerate(cuts_list):
			nbel = 2**cuts
			knotvector = createUniformMaxregularKnotvector(degree, nbel)
			# quadArgs   = {'degree': degree, 'knotvector': knotvector, 'quadrule': 'wq', 'type': 1}
			quadArgs   = {'degree': degree, 'knotvector': knotvector, 'quadrule': 'iga', 'type': 'leg'}
			args  = {'quadArgs': quadArgs, 'geoArgs': geoArgs}
			model = mechamat1D(args)

			model.activate_mechanical(matArgs)
			model.add_DirichletCondition(table=[1, 0])
			Fext    = np.zeros((model.nbctrlpts, 2*nbSteps + 1))
			Fextref = model.compute_volForce(forceVol(model.qpPhy))
			for i in range(0, nbSteps+1): Fext[:, i] = i/nbSteps*Fextref
			for i in range(nbSteps+1, 2*nbSteps+1): Fext[:, i] = (2*nbSteps - i)/nbSteps*Fextref
			
			blockPrint()
			# Solve 
			lastStep = nbSteps+1 # -2 or nbSteps+1
			disp_cp  = model.solve(Fext=Fext[:, :lastStep+1])[0]
			enablePrint()

			error[j] = model.L2NormOfError(disp_cp[:, -1], L2NormArgs={'referencePart':refPart, 'u_ref': disp_ref[:, lastStep]})		

		ax.semilogy(2**cuts_list, error, label='deg. '+str(degree))
		ax.set_ylabel('L2 Relative error (\%)')
		ax.set_xlabel('Discretization level ' + r'$h^{-1}$')
		ax.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
		ax.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())
		ax.set_xticks([8, 32, 128, 256])
		ax.set_ylim(bottom=1e-16, top=1e-1)

		ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
		fig.tight_layout()
		fig.savefig(folder + 'Fig1ErrorNBelP_' + str(quadArgs['quadrule']) +'.png')
