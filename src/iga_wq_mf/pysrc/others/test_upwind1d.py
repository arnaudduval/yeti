"""
	In this file we develop some ideas about high order nonlinear upwind in 1D.
	The following equation, du/dt = f, with u(0) = 0, it could present spurious oscillations
	when using FEM/IGA approaches. 
"""

from pysrc.lib.__init__ import *
from pysrc.lib.lib_base import createUniformKnotvector_Rmultiplicity, evalDersBasisCSRPy
from pysrc.lib.lib_quadrules import GaussQuadrature

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/'
if not os.path.isdir(folder): os.mkdir(folder)

def externalforce(qpPhy):
	f = 50*np.cos(50*qpPhy)
	return np.atleast_2d(f).T

# Time discretization
degree_tm, nbel_tm = 2, 10
knotvector_tm = createUniformKnotvector_Rmultiplicity(degree_tm, nbel_tm)
tmDscrt       = GaussQuadrature(degree_tm, knotvector_tm, {})
tmDscrt.getQuadratureRulesInfo()
tmDenseBasis, tmDenseWeights = tmDscrt.getDenseQuadRules()

# Construct of matrices and vectors
Adv  = tmDenseWeights[1] @ np.diag(np.ones(tmDscrt.nbqp)) @ tmDenseBasis[1].T
Fext = tmDenseWeights[0] @ np.diag(np.ones(tmDscrt.nbqp)) @ externalforce(tmDscrt.quadPtsPos)
usol = np.zeros(tmDscrt.nbctrlpts)

# Solve linear system
dof = np.arange(1, tmDscrt.nbctrlpts, dtype=int)
usol[dof] = np.ravel(np.linalg.solve(Adv[np.ix_(dof, dof)], Fext[dof, :]))

# Post processing
basis_tm, qp_tm = tmDscrt.getSampleBasis(sampleSize=1001)
uinterp = basis_tm[0].T @ usol

fig, ax = plt.subplots(nrows=1, ncols=1)
ax.plot(qp_tm, uinterp)
ax.grid(False)
ax.set_ylabel(r'$u(t)$')
ax.set_xlabel(r'$t$')
fig.tight_layout()
fig.savefig(folder + 'upwind.png')