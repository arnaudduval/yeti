"""
.. Test of elastoplasticity 2D
.. We test how elastoplasticity module works
.. Joaquin Cornejo 
"""

from thesis.Incremental.__init__ import *

# Set global variables
degree, cuts = 3, 6
RUNSIMU = False

if RUNSIMU:

	for quadrule, quadtype in zip(['iga', 'wq', 'wq'], ['leg', 1, 2]):
		quadArgs = {'quadrule':quadrule, 'type': quadtype}
		problem, displacement, _, intvar = simulate_2d(degree, cuts, quadArgs)
		NLres = intvar['NonLinear'][-1]
		np.savetxt(FOLDER2DATA+'nonlinearres'+quadrule+str(quadtype)+'.dat',NLres)

fig, ax = plt.subplots()
ax.set_xlabel('Number of nonlinear iterations')
ax.set_ylabel('Relative nonlinear residue')

for quadrule, quadtype in zip(['iga', 'wq', 'wq'], ['leg', 1, 2]):
	NLres = np.loadtxt(FOLDER2DATA+'nonlinearres'+quadrule+str(quadtype)+'.dat')
	if quadtype == 'leg': label = 'IGA-GL'
	elif quadtype == 1: label = 'IGA-WQ 1'
	else: label = 'IGA-WQ 2'
	ax.semilogy(np.array(NLres)/NLres[0], label=label)
ax.set_xlim([0, 8])
ax.set_ylim([1e-10, 1])
ax.legend()
fig.savefig(FOLDER2SAVE+'NonlinearResidual'+'.pdf')