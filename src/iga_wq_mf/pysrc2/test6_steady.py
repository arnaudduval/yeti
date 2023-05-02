from lib.__init__ import *
from lib.lib_load import *
from lib.lib_simulation import *

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/test6/'
if not os.path.isdir(folder): os.mkdir(folder)

dataExist   = True
degree_list = np.arange(6, 7)
cuts_list   = np.arange(5, 7)
name_list   = ['cb', 'vb', 'tr', 'rqa']
IterMethods = ["WP", "C", "JMC", "TDC"]

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
				heatmat = thermomat()
				heatmat.addConductivity(np.array([[1, 0.5, 0.1],[0.5, 2, 0.25], [0.1, 0.25, 3]]),
										isIsotropic=True, shape=(3, 3))				
				table    = np.ones((3, 2), dtype=bool)
				boundary = step(simu._part._nbctrlpts)
				boundary.add_DirichletTemperature(table=table)
				simu.simulate(material=heatmat, boundary=boundary, overwrite=True)

			else :
				filename   = simu._filename
				simuOutput = decoder(filename)
				simuOutput.plot_results(extension='.png')
