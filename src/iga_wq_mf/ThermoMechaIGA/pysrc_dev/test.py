"""
.. Joaquin Cornejo 
"""
import numpy as np
from lib.create_geomdl import geomdlModel
from lib.create_model import thermoMechaModel

# Ceate geometry
name = 'VB'
Geomodel = geomdlModel(name=name)
IGAmodel = Geomodel.export_IGAparametrization(nb_refinementByDirection=np.array([2, 2, 2]))

# Create model
phymodel = thermoMechaModel(IGAmodel)
nnz = phymodel._nb_ctrlpts_total
u_ctrlpts = np.arange(nnz)
phymodel.export_results(u_ctrlpts= u_ctrlpts)

