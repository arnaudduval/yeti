"""
.. This module helps to read information from a given geometry (YETI format)
.. and to create thermo-mechanical model (material, boundary conditions, etc.)
.. Joaquin Cornejo
"""

from .__init__ import *

# My libraries
from .create_material import *

# class thermoMechaModel(): 

# 	def __init__(self):
# 		self._geometry = None
# 		self._meshDiscretization = None
# 		return
	
# 	class material():
# 		def __init__(self):
# 			return
	
# 	def addGeometry(self, name=None, filename=None, **kwargs):
# 		self._geometry = geomdlModel(name=name, filename=filename, **kwargs)
# 		return
	
# 	def addMeshDiscretization(self, mesh):
# 		return