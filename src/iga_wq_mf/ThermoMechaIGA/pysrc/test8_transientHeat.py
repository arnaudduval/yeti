"""
.. Test of transient heat solver
.. Joaquin Cornejo 
"""

from lib.__init__ import *
from lib.physics import *
from lib.create_geomdl import geomdlModel
from lib.fortran_mf_wq import fortran_mf_wq

# Set global variables
degree, cuts = 4, 3

# Create geometry 
geometry = {'degree':[degree, degree, degree]}
modelGeo = geomdlModel('VB', **geometry)
modelIGA = modelGeo.export_IGAparametrization(nb_refinementByDirection=
                                            np.array([cuts, cuts, cuts]))

# Create model
modelPhy = fortran_mf_wq(modelIGA)

# Add material 
material = {'capacity':1, 'conductivity':np.eye(3)}
modelPhy._set_material(material)

# Block boundaries
Dirichlet = {'thermal':np.array([[1, 1], [1, 1], [1, 1]])}
modelPhy._set_dirichlet_boundaries(Dirichlet)

# Add external force 
Fvol = modelPhy.eval_source_vector(powden_prism)
