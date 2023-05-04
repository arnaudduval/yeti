from lib.__init__ import *
from lib.lib_load import *
from lib.lib_simulation import *

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/test6/'
if not os.path.isdir(folder): os.mkdir(folder)

dataExist   = False
degree_list = np.arange(6, 7)
cuts_list   = np.arange(5, 6)
name_list   = ['cb', 'vb', 'tr', 'rqa']
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
	Kprop[0, 0, :] += 0.5*np.abs(y)**1.5 
	Kprop[1, 1, :] += np.abs(z)
	Kprop[2, 2, :] += np.abs(x)**0.5
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
				simu.create_model()
				mat = thermomat()
				# mat.addConductivity(np.array([[1, 0.5, 0.1],[0.5, 2, 0.25], [0.1, 0.25, 3]]),
				# 						isIsotropic=True, shape=(3, 3))	
				mat.addConductivity(setKprop, isIsotropic=False)				
				table    = np.ones((3, 2), dtype=bool)
				boundary = step(simu._part._nbctrlpts)
				boundary.add_DirichletTemperature(table=table)
				simu.simulate(material=mat, boundary=boundary, overwrite=True)

			else :
				filename   = simu._filename
				simuOutput = decoder(filename)
				simuOutput.plot_results(extension='.png', plotLegend=False)
