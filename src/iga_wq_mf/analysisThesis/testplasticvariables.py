from pysrc.lib.__init__ import *
from pysrc.lib.lib_geomdl import Geomdl
from pysrc.lib.lib_part import part
from pysrc.lib.lib_material import mechamat
from pysrc.lib.lib_boundary import boundaryCondition
from pysrc.lib.lib_job import mechaproblem

# Set global variables
YOUNG, POISSON = 2400, 0.2
H, beta = 100, 0.3
MATARGS = {'elastic_modulus':YOUNG, 'elastic_limit':300, 'poisson_ratio': POISSON, 
			'isoHardLaw': {'name':'linear', 'Eiso':(1-beta)*H}, 
			'kineHardLaw':{'parameters':np.array([[2/3*beta*H, 0]])}
			}

material = mechamat(MATARGS)
strain = np.zeros((6, 1))
strain[0, 0] = 0.125 + 0.1; strain[1, 0] = -0.025 - 0.02; strain[2, 0] = -0.025 - 0.02
pls_n0 = np.zeros((6, 1))
a_n0   = np.zeros((1, 1))
b_n0   = np.zeros((1, 6, 1))
out, _ = material.J2returnMappingAlgorithm3D(strain, pls_n0, a_n0, b_n0)
lame_mu = material.lame_mu
dgamma, norm_eta_trial = 0.0948, 180*np.sqrt(6)
c1 = 4*lame_mu**2/(2*lame_mu+2/3*H)
c2 = 4*lame_mu**2*dgamma/norm_eta_trial
print(c1/(2*lame_mu), c2/(2*lame_mu))