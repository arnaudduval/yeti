from lib.__init__ import *
from lib.lib_load import *
from lib.lib_simulation import *
from lib.lib_base import relativeError

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/t3dsteady/'
if not os.path.isdir(folder): os.mkdir(folder)

dataExist   = True
degree_list = np.arange(6, 7)
cuts_list   = np.arange(4, 5)
# name_list   = ['cb', 'vb', 'tr', 'rqa']
name_list   = ['vb', 'rqa']
IterMethods = ["WP", "C", "JMC", "TDC"]

def setKprop(P:list):
	x = P[0, :]
	y = P[1, :]
	z = P[2, :]
	Kref  = np.array([[1, 0.5, 0.1],[0.5, 2, 0.25], [0.1, 0.25, 3]])
	Kprop = np.zeros((3, 3, len(x)))
	for i in range(3): 
		for j in range(3):
			Kprop[i, j, :] = Kref[i, j] 
	Kprop[0, 0, :] += 0.75*np.cos(np.pi*y)
	Kprop[1, 1, :] += 2*np.exp(-(z-0.5)**2)
	Kprop[2, 2, :] += 2.5*np.cos(np.pi*x)**2
	# for i in range(len(x)):
	# 	if not np.all(np.linalg.eigvals(Kprop[:, :, i]) > 0):
	# 		raise Warning('It is not a positive matrix')
	return Kprop 

for cuts in cuts_list:
	for degree in degree_list:
		for name in name_list: 

			if name   == 'cb' : funpow, funtemp = powden_cube, None 
			elif name == 'vb' : funpow, funtemp = powden_prism, None 
			elif name == 'tr' : funpow, funtemp = powden_thickring, None 
			elif name == 'rqa': funpow, funtemp = powden_rotring, temperature_rotring 

			inputs = {'degree': degree, 'cuts': cuts, 'name': name, 'isGauss': False, 
					'funPowerDensity': funpow, 'funTemperature': funtemp, 'IterMethods': IterMethods,
					'folder': folder}
			simu   = encoder(**inputs)  
			
			if not dataExist:
				mat = thermomat()
				mat.addConductivity(setKprop, isIsotropic=False)				
				table    = np.ones((3, 2), dtype=bool)
				boundary = boundaryCondition(simu._nbctrlpts*np.ones(3, dtype=int))
				boundary.add_DirichletTemperature(table=table)
				un = simu.simulate(material=mat, boundary=boundary, overwrite=True)
				np.save(folder+'solution', un)

			else :
				filename   = simu._filename
				simuOutput = decoder(filename)
				simuOutput.plot_results(extension='.png', plotLegend=True)
				
