from pysrc.lib.__init__ import *
from pysrc.lib.lib_part import part
from pysrc.lib.lib_boundary import boundaryCondition
from pysrc.lib.lib_job3d import problem

FOLDER = os.path.dirname(os.path.realpath(__file__)) + '/'

DEGREE, CUTS = 0, 1
QUADARGS = {'quadrule':'iga', 'type':'leg'}

def exactfunction(P:list):
	x = P[0, :]; y = P[1, :]
	F = np.sin(np.pi*x)*np.sin(np.pi*y)*(x+4)*(y-4)*(x**2+y**2-1)
	return F

# Creation of the IGA object
# --------------------------
modelIGA = IGAparametrization(filename=FOLDER+'platewithhole')
modelIGA.refine(nb_refinementByDirection=CUTS*np.array([1, 1, 1]), 
				nb_degreeElevationByDirection=DEGREE*np.array([1, 1, 1]),)
modelPhy = part(modelIGA, quadArgs=QUADARGS)
NURBSweights = modelIGA._vectWeight

boundary = boundaryCondition()
Mproblem = problem(modelPhy, boundary, {})
NURBSfactor = Mproblem._interpolate_meshgrid(NURBSweights)

u_atQuadpts = ...
u_atCtrlpts = Mproblem.L2projectionCtrlpts(u_atQuadpts, table=np.ones((2, 2)), prop=...)
