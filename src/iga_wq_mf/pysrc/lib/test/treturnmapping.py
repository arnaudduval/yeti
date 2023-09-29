from pysrc.lib.__init__ import *
from pysrc.lib.lib_base import createUniformCurve
from pysrc.lib.lib_1d import mechanics1D
from pysrc.lib.lib_material import mechamat, symtensor2array4All

# # 1 DIMENSION
# length = 1
# degree, nbel = 1, 1
# crv = createUniformCurve(degree, nbel, length)
# args = {'quadArgs': {'quadrule': 'iga', 'type': 'leg'}}
# modelPhy = mechamat1D(crv, args)

# matArgs = {'elastic_modulus':2e11, 'elastic_limit':250e6, 'plasticLaw': {'Isoname': 'linear', 'Eiso':25e9}}
# modelPhy.activate_mechanical(matArgs)
# a_n0   = 0.0001*np.ones(2)
# pls_n0 = 0.0001*np.ones(2)
# b_n0   = np.zeros(2)
# strain = 0.00075 + a_n0 + 0.002
# output = modelPhy.returnMappingAlgorithm(strain, pls_n0, a_n0, b_n0)
# stress, pls_n1, a_n1, b_n1, Cep = output[0, :], output[1, :], output[2, :], output[3, :], output[4, :]
# stress_ref = 285.6e6
# a_n1_ref  = 1.422e-3

# if np.any(np.abs(stress-stress_ref)/stress_ref >= 1e-3) : raise Warning('Problem')
# if np.any(np.abs(a_n1-a_n1_ref)/a_n1_ref >= 1e-3) : raise Warning('Problem')

# 3 DIMENSIONS
E  = 2.4e9
nu = 0.2
sigma_Y0 = 300e6
Eiso = 100e6
matArgs = {'elastic_modulus':E, 'elastic_limit':sigma_Y0, 'poisson_ratio': nu, 
			'plasticLaw': {'Isoname':'linear', 'Eiso':Eiso}}
s = 300e6
Tstrain = np.zeros((3, 3, 1))
Tstrain[:, :, 0] = np.diag([0.1, -0.02, -0.02]) + np.diag([s/E, -nu*s/E, -nu*s/E])
strain = symtensor2array4All(Tstrain, 3)
a_n0 = np.array([0.0])
pls_n0 = np.zeros((6, 1))
b_n0 = np.zeros((6, 1))
material = mechamat(matArgs)
material.returnMappingAlgorithm(strain, pls_n0, a_n0, b_n0)