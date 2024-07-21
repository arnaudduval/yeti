from pysrc.lib.__init__ import *
from pysrc.lib.lib_base import createUniformOpenCurve
from pysrc.lib.lib_geomdl import Geomdl
from pysrc.lib.lib_part import part1D, part
from pysrc.lib.lib_job1d import mechaproblem1D, heatproblem1D, problem1D
from pysrc.lib.lib_job3d import mechaproblem, heatproblem
from pysrc.lib.lib_material import mechamat
from pysrc.lib.lib_boundary import boundaryCondition
import pickle
from pyevtk.vtk import VtkGroup

FOLDER2SAVE = os.path.dirname(os.path.realpath(__file__)) + '/results/'
if not os.path.isdir(FOLDER2SAVE): os.mkdir(FOLDER2SAVE)

FOLDER2DATA = FOLDER2SAVE + '/datafromsimu/'
if not os.path.isdir(FOLDER2DATA): os.mkdir(FOLDER2DATA)

def run(filename=None, folder=None, nbFiles=1):
	assert folder is not None, 'Folder unknown'
	if filename is None: filename = 'out'
	print("Running group...")
	g = VtkGroup(folder)
	for i in range(nbFiles):
		g.addFile(filepath = folder + filename + str(i) + '.vts', sim_time = i)
	g.save()
	return
