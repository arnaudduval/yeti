
from pysrc.lib.__init__ import *
from pysrc.lib.lib_base import createUniformOpenCurve, createUniformKnotvector_Rmultiplicity
from pysrc.lib.lib_geomdl import Geomdl
from pysrc.lib.lib_part import part, part1D
from pysrc.lib.lib_boundary import boundaryCondition
from pysrc.lib.lib_material import heatmat
from pysrc.lib.lib_job3d import heatproblem, stheatproblem
from pysrc.lib.lib_job1d import heatproblem1D, stheatproblem1D
from numpy import pi, sin, cos, abs, exp, sign, tanh

GEONAME = 'QA'
IS1DIM = False
ISLINEAR = False
CUTS_TIME = 6

if GEONAME == 'SQ' or IS1DIM: ISISOTROPIC = True
else: ISISOTROPIC = False

if GEONAME == 'SQ' or IS1DIM: CST = 100
else: CST = 50

def nonlinearfunc(args:dict):
	temperature = args.get('temperature')
	if ISLINEAR: Kprop1d = 2*np.ones(shape=np.shape(temperature))
	else: Kprop1d = 3.0 + 2.0*tanh(temperature/50)
	return np.atleast_2d(Kprop1d)

def nonlineardersfunc(args:dict):
	temperature = args.get('temperature')
	if ISLINEAR: Kprop = np.zeros(shape=np.shape(temperature))
	else: Kprop = (2.0/50)/(np.cosh(temperature/50))**2
	return Kprop

def createAsymmetricalKnotvector(level, xasym=0.1):

	def discretize(array, level):
		n = len(array) 
		newarray = np.copy(array)
		if level == 2:
			for i in range(n-1): 
				newarray = np.append(newarray, (array[i]+array[i+1])/2)
			newarray = np.sort(newarray)
		else:
			for _ in range(1, level): 
				tmp = discretize(newarray, 2)
				newarray = np.copy(tmp)
		return newarray
	
	assert level > 0, 'Not possible. Try higher level'
	assert xasym < 0.25, 'Try lower value'
	knotvector = np.array([0.0, xasym, 0.5-xasym/2, 0.5+xasym/2, 1.0-xasym, 1.0]) # This is level 1
	if level > 1: 
		tmp = discretize(knotvector, level)
		knotvector = np.copy(tmp)

	return knotvector

def createAsymmetricalCurve(degree, level, length, xasym=0.05):
	crv = BSpline.Curve()
	crv.degree  = degree
	crv.ctrlpts = [[i*length/degree, 0.0] for i in range(degree+1)]
	crv.knotvector = createUniformKnotvector_Rmultiplicity(degree, 1)
	knotList = createAsymmetricalKnotvector(level, xasym=xasym)
	for knot in knotList[1:-1]:
		operations.insert_knot(crv, [knot], [1])
	return crv

def capacityProperty(args:dict):
	temperature = args.get('temperature')
	capacity = np.ones(shape=np.shape(temperature))
	return capacity

def capacityDersProperty(args:dict):
	temperature = args.get('temperature')
	capacity = np.zeros(shape=np.shape(temperature))
	return capacity

def conductivityProperty(args:dict, isIsotropic=ISISOTROPIC):
	temperature = args.get('temperature')
	Kprop1d = np.ravel(nonlinearfunc(args), order='F')
	if IS1DIM: 
		return Kprop1d
	else:
		reference = np.array([[1., 0.0],[0.0, 1.0]])
		if not isIsotropic: reference = np.array([[1., 0.5],[0.5, 2.0]])
		Kprop2d = np.zeros((2, 2, len(temperature)))
		for i in range(2): 
			for j in range(2):
				Kprop2d[i, j, :] = reference[i, j]*Kprop1d
		return Kprop2d
	
def conductivityDersProperty(args:dict, isIsotropic=ISISOTROPIC):
	temperature = args.get('temperature')
	Kprop1d = np.ravel(nonlineardersfunc(args), order='F')
	if IS1DIM: 
		return Kprop1d
	else:
		reference = np.array([[1., 0.0],[0.0, 1.0]])
		if not isIsotropic: reference = np.array([[1., 0.5],[0.5, 2.0]])
		Kprop2d = np.zeros((2, 2, len(temperature)))
		for i in range(2): 
			for j in range(2):
				Kprop2d[i, j, :] = reference[i, j]*Kprop1d
		return Kprop2d
	
# Square shape

def exactTemperatureSquare_inc(args:dict):
	t = args['time']
	if IS1DIM: x = args['position']
	else: x = args['position'][0, :]
	u = CST*sin(2*pi*x)*sin(pi/2*t)*(1+0.75*cos(3*pi/2*t)) 
	return u

def exactTemperatureSquare_spt(qpPhy):
	x = qpPhy[0, :]; t = qpPhy[-1, :]
	u = CST*sin(2*pi*x)*sin(pi/2*t)*(1+0.75*cos(3*pi/2*t)) 
	return u

def powerDensitySquare_inc(args:dict):
	t = args['time']
	if IS1DIM: x = args['position']
	else: x = args['position'][0, :]

	if not ISISOTROPIC: raise Warning('Not possible')
	u = sin((pi*t)/2)*sin(2*pi*x)*((3*cos((3*pi*t)/2))/4 + 1) 
	if ISLINEAR:
		f = (
			(CST*pi*cos((pi*t)/2)*sin(2*pi*x)*((3*cos((3*pi*t)/2))/4 + 1))/2 
			- (9*CST*pi*sin((pi*t)/2)*sin((3*pi*t)/2)*sin(2*pi*x))/8 
			+ 8*CST*pi**2*u
		)
	else: 
		f = ((CST*pi*cos((pi*t)/2)*sin(2*pi*x)*((3*cos((3*pi*t)/2))/4 + 1))/2 
			- (9*CST*pi*sin((pi*t)/2)*sin((3*pi*t)/2)*sin(2*pi*x))/8 
			+ 4*CST*pi**2*u*(2*tanh((CST*u)/50) + 3) 
			+ (4*CST**2*pi**2*cos(2*pi*x)**2*sin((pi*t)/2)**2*((3*cos((3*pi*t)/2))/4 + 1)**2*(tanh((CST*u)/50)**2 - 1))/25
		)	

	return f

def powerDensitySquare_spt(args:dict):
	time = args['time']
	position = args['position']
	if IS1DIM: nc_sp = len(position)
	else: nc_sp = np.size(position, axis=1)
	nc_tm = np.size(time); f = np.zeros((nc_sp, nc_tm))
	for i in range(nc_tm):
		t = time[i]
		f[:, i] = powerDensitySquare_inc(args={'time':t, 'position':position})
	return np.ravel(f, order='F')

# Trapezoidal shape

def exactTemperatureTrap_inc(args:dict):
	t = args['time']
	if IS1DIM: raise Warning('Try higher dimension')
	x = args['position'][0, :]
	y = args['position'][1, :]
	u = CST*sin(pi*y)*sin(pi*(y+0.75*x-0.5)*(-y+0.75*x-0.5))*sin(5*pi*x)*sin(pi/2*t)*(1+0.75*cos(3*pi/2*t)) 
	return u

def exactTemperatureTrap_spt(qpPhy):
	x = qpPhy[0, :]; y=qpPhy[1, :]; t = qpPhy[-1, :]
	u = CST*sin(pi*y)*sin(pi*(y+0.75*x-0.5)*(-y+0.75*x-0.5))*sin(5*pi*x)*sin(pi/2*t)*(1+0.75*cos(3*pi/2*t)) 
	return u

def powerDensityTrap_inc(args:dict):
	t = args['time']
	x = args['position'][0, :]
	y = args['position'][1, :]
	if ISISOTROPIC: raise Warning('Not possible')
	
	u1 = pi*(y - (3*x)/4 + 1/2)*((3*x)/4 + y - 1/2); u2 = sin((pi*t)/2); u3 = sin(5*pi*x)
	u4 = ((3*cos((3*pi*t)/2))/4 + 1); u5 = pi*(y - (3*x)/4 + 1/2) + pi*((3*x)/4 + y - 1/2)
	u6 = (3*pi*(y - (3*x)/4 + 1/2))/4 - (3*pi*((3*x)/4 + y - 1/2))/4
	if ISLINEAR:
		f = ((23*CST*pi*cos(u1)*u2*u3*sin(pi*y)*u4)/4 
			- (CST*pi*sin(u1)*cos((pi*t)/2)*u3*sin(pi*y)*u4)/2 
			- 4*CST*sin(u1)*u2*u3*sin(pi*y)*(u5)**2*u4 
			- 2*CST*sin(u1)*u2*u3*sin(pi*y)*(u6)**2*u4 
			+ 10*CST*pi**2*sin(u1)*cos(5*pi*x)*cos(pi*y)*u2*u4 
			+ (9*CST*pi*sin(u1)*u2*sin((3*pi*t)/2)*u3*sin(pi*y))/8 
			- 54*CST*pi**2*sin(u1)*u2*u3*sin(pi*y)*u4 
			- 2*CST*sin(u1)*u2*u3*sin(pi*y)*(u5)*(u6)*u4 
			+ 10*CST*pi*cos(u1)*cos(5*pi*x)*u2*sin(pi*y)*(u5)*u4 
			+ 8*CST*pi*cos(u1)*cos(pi*y)*u2*u3*(u5)*u4 
			+ 20*CST*pi*cos(u1)*cos(5*pi*x)*u2*sin(pi*y)*(u6)*u4 
			+ 2*CST*pi*cos(u1)*cos(pi*y)*u2*u3*(u6)*u4
		)
	else: 
		f = ((2*tanh((CST*sin(u1)*u2*u3*sin(pi*y)*u4)/50) - 3)*((9*CST*pi*cos(u1)*u2*u3*sin(pi*y)*u4)/8 
					+ CST*sin(u1)*u2*u3*sin(pi*y)*(u6)**2*u4 
					+ 25*CST*pi**2*sin(u1)*u2*u3*sin(pi*y)*u4 
					- 10*CST*pi*cos(u1)*cos(5*pi*x)*u2*sin(pi*y)*(u6)*u4) 
			- (4*tanh((CST*sin(u1)*u2*u3*sin(pi*y)*u4)/50) - 6)*(2*CST*pi*cos(u1)*u2*u3*sin(pi*y)*u4 
					- CST*sin(u1)*u2*u3*sin(pi*y)*(u5)**2*u4 
					- CST*pi**2*sin(u1)*u2*u3*sin(pi*y)*u4 
					+ 2*CST*pi*cos(u1)*cos(pi*y)*u2*u3*(u5)*u4) 
			- 2*(tanh((CST*sin(u1)*u2*u3*sin(pi*y)*u4)/50) - 3/2)*(5*CST*pi**2*sin(u1)*cos(5*pi*x)*cos(pi*y)*u2*u4 
					- CST*sin(u1)*u2*u3*sin(pi*y)*(u5)*(u6)*u4 
					+ 5*CST*pi*cos(u1)*cos(5*pi*x)*u2*sin(pi*y)*(u5)*u4 
					+ CST*pi*cos(u1)*cos(pi*y)*u2*u3*(u6)*u4) 
			+ (CST*cos(u1)*u2*u3*sin(pi*y)*(u5)*u4 + CST*pi*sin(u1)*cos(pi*y)*u2*u3*u4)*((CST*cos(u1)*u2*u3*sin(pi*y)*(u6)*u4)/50 
					+ (CST*pi*sin(u1)*cos(5*pi*x)*u2*sin(pi*y)*u4)/10)*(tanh((CST*sin(u1)*u2*u3*sin(pi*y)*u4)/50)**2 - 1) 
			+ 2*(CST*cos(u1)*u2*u3*sin(pi*y)*(u6)*u4 + 5*CST*pi*sin(u1)*cos(5*pi*x)*u2*sin(pi*y)*u4)*((CST*cos(u1)*u2*u3*sin(pi*y)*(u6)*u4)/50 
					+ (CST*pi*sin(u1)*cos(5*pi*x)*u2*sin(pi*y)*u4)/10)*(tanh((CST*sin(u1)*u2*u3*sin(pi*y)*u4)/50)**2 - 1) 
			+ 4*(CST*cos(u1)*u2*u3*sin(pi*y)*(u5)*u4 + CST*pi*sin(u1)*cos(pi*y)*u2*u3*u4)*((CST*cos(u1)*u2*u3*sin(pi*y)*(u5)*u4)/50 
					+ (CST*pi*sin(u1)*cos(pi*y)*u2*u3*u4)/50)*(tanh((CST*sin(u1)*u2*u3*sin(pi*y)*u4)/50)**2 - 1) 
					+ (CST*cos(u1)*u2*u3*sin(pi*y)*(u6)*u4 + 5*CST*pi*sin(u1)*cos(5*pi*x)*u2*sin(pi*y)*u4)*((CST*cos(u1)*u2*u3*sin(pi*y)*(u5)*u4)/50 
							+ (CST*pi*sin(u1)*cos(pi*y)*u2*u3*u4)/50)*(tanh((CST*sin(u1)*u2*u3*sin(pi*y)*u4)/50)**2 - 1) 
							- (CST*pi*sin(u1)*cos((pi*t)/2)*u3*sin(pi*y)*u4)/2 + (9*CST*pi*sin(u1)*u2*sin((3*pi*t)/2)*u3*sin(pi*y))/8
		)

	return f

def powerDensityTrap_spt(args:dict):
	time = args['time']
	position = args['position']
	nc_sp = np.size(position, axis=1); nc_tm = np.size(time); f = np.zeros((nc_sp, nc_tm))
	for i in range(nc_tm):
		t = time[i]
		f[:, i] = powerDensityTrap_inc(args={'time':t, 'position':position})
	return np.ravel(f, order='F')

# Ring shape

def exactTemperatureRing_inc(args:dict):
	t = args['time']
	if IS1DIM: raise Warning('Try higher dimension')
	x = args['position'][0, :]
	y = args['position'][1, :]
	u = -CST*tanh(x**2+y**2-1.0)*sin(pi*(x**2+y**2-0.25**2))*sin(pi*x*y)*sin(pi/2*t)*(1+0.75*cos(3*pi/2*t))
	return u

def exactTemperatureRing_spt(qpPhy):
	x = qpPhy[0, :]; y=qpPhy[1, :]; t = qpPhy[-1, :]
	u = -CST*tanh(x**2+y**2-1.0)*sin(pi*(x**2+y**2-0.25**2))*sin(pi*x*y)*sin(pi/2*t)*(1+0.75*cos(3*pi/2*t))
	return u

def powerDensityRing_inc(args:dict):
	t = args['time']
	x = args['position'][0, :]
	y = args['position'][1, :]

	if ISISOTROPIC: raise Warning('Not possible')
	u1 = pi*(x**2 + y**2 - 1/16); u2 = x**2 + y**2 - 1; u3 = (3*cos((3*pi*t)/2))/4 + 1; u4 = sin((pi*t)/2)
	if ISLINEAR:
		f = (2*CST*pi*cos(pi*x*y)*sin(u1)*tanh(u2)*u4*(u3) 
		- 12*CST*sin(pi*x*y)*sin(u1)*u4*(tanh(u2)**2 - 1)*(u3) 
		+ 12*CST*pi*sin(pi*x*y)*cos(u1)*tanh(u2)*u4*(u3) 
		- (CST*pi*sin(pi*x*y)*sin(u1)*tanh(u2)*cos((pi*t)/2)*(u3))/2 
		+ (9*CST*pi*sin(pi*x*y)*sin(u1)*tanh(u2)*u4*sin((3*pi*t)/2))/8 
		- 4*CST*x**2*pi*cos(pi*x*y)*sin(u1)*u4*(tanh(u2)**2 - 1)*(u3) 
		- 16*CST*x**2*pi*sin(pi*x*y)*cos(u1)*u4*(tanh(u2)**2 - 1)*(u3) 
		- 4*CST*y**2*pi*cos(pi*x*y)*sin(u1)*u4*(tanh(u2)**2 - 1)*(u3) 
		- 32*CST*y**2*pi*sin(pi*x*y)*cos(u1)*u4*(tanh(u2)**2 - 1)*(u3) 
		+ 16*CST*x**2*sin(pi*x*y)*sin(u1)*tanh(u2)*u4*(tanh(u2)**2 - 1)*(u3) 
		+ 32*CST*y**2*sin(pi*x*y)*sin(u1)*tanh(u2)*u4*(tanh(u2)**2 - 1)*(u3) 
		+ 4*CST*x**2*pi**2*cos(pi*x*y)*cos(u1)*tanh(u2)*u4*(u3) 
		+ 4*CST*y**2*pi**2*cos(pi*x*y)*cos(u1)*tanh(u2)*u4*(u3) 
		- 12*CST*x**2*pi**2*sin(pi*x*y)*sin(u1)*tanh(u2)*u4*(u3) 
		- 18*CST*y**2*pi**2*sin(pi*x*y)*sin(u1)*tanh(u2)*u4*(u3) 
		- 24*CST*x*y*pi*cos(pi*x*y)*sin(u1)*u4*(tanh(u2)**2 - 1)*(u3) 
		- 16*CST*x*y*pi*sin(pi*x*y)*cos(u1)*u4*(tanh(u2)**2 - 1)*(u3) 
		+ 16*CST*x*y*sin(pi*x*y)*sin(u1)*tanh(u2)*u4*(tanh(u2)**2 - 1)*(u3) 
		+ 24*CST*x*y*pi**2*cos(pi*x*y)*cos(u1)*tanh(u2)*u4*(u3) 
		- 10*CST*x*y*pi**2*sin(pi*x*y)*sin(u1)*tanh(u2)*u4*(u3)
	)
	else: 

		f = ((2*tanh((CST*sin(pi*x*y)*sin(u1)*tanh(u2)*u4*(u3))/50) - 3)*(2*CST*sin(pi*x*y)*sin(u1)*u4*(tanh(u2)**2 - 1)*(u3) 
					- 2*CST*pi*sin(pi*x*y)*cos(u1)*tanh(u2)*u4*(u3) 
					+ 8*CST*x**2*pi*sin(pi*x*y)*cos(u1)*u4*(tanh(u2)**2 - 1)*(u3) 
					- 8*CST*x**2*sin(pi*x*y)*sin(u1)*tanh(u2)*u4*(tanh(u2)**2 - 1)*(u3) 
					+ 4*CST*x**2*pi**2*sin(pi*x*y)*sin(u1)*tanh(u2)*u4*(u3) 
					+ CST*y**2*pi**2*sin(pi*x*y)*sin(u1)*tanh(u2)*u4*(u3) 
					+ 4*CST*x*y*pi*cos(pi*x*y)*sin(u1)*u4*(tanh(u2)**2 - 1)*(u3) 
					- 4*CST*x*y*pi**2*cos(pi*x*y)*cos(u1)*tanh(u2)*u4*(u3)) 
			- 2*(tanh((CST*sin(pi*x*y)*sin(u1)*tanh(u2)*u4*(u3))/50) - 3/2)*(CST*pi*cos(pi*x*y)*sin(u1)*tanh(u2)*u4*(u3) 
					- 2*CST*x**2*pi*cos(pi*x*y)*sin(u1)*u4*(tanh(u2)**2 - 1)*(u3) 
					- 2*CST*y**2*pi*cos(pi*x*y)*sin(u1)*u4*(tanh(u2)**2 - 1)*(u3) 
					+ 2*CST*x**2*pi**2*cos(pi*x*y)*cos(u1)*tanh(u2)*u4*(u3) 
					+ 2*CST*y**2*pi**2*cos(pi*x*y)*cos(u1)*tanh(u2)*u4*(u3) 
					- 8*CST*x*y*pi*sin(pi*x*y)*cos(u1)*u4*(tanh(u2)**2 - 1)*(u3) 
					+ 8*CST*x*y*sin(pi*x*y)*sin(u1)*tanh(u2)*u4*(tanh(u2)**2 - 1)*(u3) 
					- 5*CST*x*y*pi**2*sin(pi*x*y)*sin(u1)*tanh(u2)*u4*(u3)) 
			+ (4*tanh((CST*sin(pi*x*y)*sin(u1)*tanh(u2)*u4*(u3))/50) - 6)*(2*CST*sin(pi*x*y)*sin(u1)*u4*(tanh(u2)**2 - 1)*(u3) 
							- 2*CST*pi*sin(pi*x*y)*cos(u1)*tanh(u2)*u4*(u3) 
							+ 8*CST*y**2*pi*sin(pi*x*y)*cos(u1)*u4*(tanh(u2)**2 - 1)*(u3) 
							- 8*CST*y**2*sin(pi*x*y)*sin(u1)*tanh(u2)*u4*(tanh(u2)**2 - 1)*(u3) 
							+ CST*x**2*pi**2*sin(pi*x*y)*sin(u1)*tanh(u2)*u4*(u3) 
							+ 4*CST*y**2*pi**2*sin(pi*x*y)*sin(u1)*tanh(u2)*u4*(u3) 
							+ 4*CST*x*y*pi*cos(pi*x*y)*sin(u1)*u4*(tanh(u2)**2 - 1)*(u3) 
							- 4*CST*x*y*pi**2*cos(pi*x*y)*cos(u1)*tanh(u2)*u4*(u3)) 
			+ 2*(tanh((CST*sin(pi*x*y)*sin(u1)*tanh(u2)*u4*(u3))/50)**2 - 1)*(2*CST*x*pi*sin(pi*x*y)*cos(u1)*tanh(u2)*u4*(u3) 
							- 2*CST*x*sin(pi*x*y)*sin(u1)*u4*(tanh(u2)**2 - 1)*(u3) 
							+ CST*y*pi*cos(pi*x*y)*sin(u1)*tanh(u2)*u4*(u3))*((CST*x*pi*sin(pi*x*y)*cos(u1)*tanh(u2)*u4*(u3))/25 
									- (CST*x*sin(pi*x*y)*sin(u1)*u4*(tanh(u2)**2 - 1)*(u3))/25 
									+ (CST*y*pi*cos(pi*x*y)*sin(u1)*tanh(u2)*u4*(u3))/50) 
			+ (tanh((CST*sin(pi*x*y)*sin(u1)*tanh(u2)*u4*(u3))/50)**2 - 1)*(2*CST*x*pi*sin(pi*x*y)*cos(u1)*tanh(u2)*u4*(u3) 
					- 2*CST*x*sin(pi*x*y)*sin(u1)*u4*(tanh(u2)**2 - 1)*(u3) 
					+ CST*y*pi*cos(pi*x*y)*sin(u1)*tanh(u2)*u4*(u3))*((CST*x*pi*cos(pi*x*y)*sin(u1)*tanh(u2)*u4*(u3))/50 
							- (CST*y*sin(pi*x*y)*sin(u1)*u4*(tanh(u2)**2 - 1)*(u3))/25 
							+ (CST*y*pi*sin(pi*x*y)*cos(u1)*tanh(u2)*u4*(u3))/25) 
			+ (tanh((CST*sin(pi*x*y)*sin(u1)*tanh(u2)*u4*(u3))/50)**2 - 1)*((CST*x*pi*sin(pi*x*y)*cos(u1)*tanh(u2)*u4*(u3))/25
					- (CST*x*sin(pi*x*y)*sin(u1)*u4*(tanh(u2)**2 - 1)*(u3))/25 
					+ (CST*y*pi*cos(pi*x*y)*sin(u1)*tanh(u2)*u4*(u3))/50)*(CST*x*pi*cos(pi*x*y)*sin(u1)*tanh(u2)*u4*(u3) 
							- 2*CST*y*sin(pi*x*y)*sin(u1)*u4*(tanh(u2)**2 - 1)*(u3) 
							+ 2*CST*y*pi*sin(pi*x*y)*cos(u1)*tanh(u2)*u4*(u3)) 
			+ 4*(tanh((CST*sin(pi*x*y)*sin(u1)*tanh(u2)*u4*(u3))/50)**2 - 1)*(CST*x*pi*cos(pi*x*y)*sin(u1)*tanh(u2)*u4*(u3) 
					- 2*CST*y*sin(pi*x*y)*sin(u1)*u4*(tanh(u2)**2 - 1)*(u3) 
					+ 2*CST*y*pi*sin(pi*x*y)*cos(u1)*tanh(u2)*u4*(u3))*((CST*x*pi*cos(pi*x*y)*sin(u1)*tanh(u2)*u4*(u3))/50 
							- (CST*y*sin(pi*x*y)*sin(u1)*u4*(tanh(u2)**2 - 1)*(u3))/25 
							+ (CST*y*pi*sin(pi*x*y)*cos(u1)*tanh(u2)*u4*(u3))/25)
			- (CST*pi*sin(pi*x*y)*sin(u1)*tanh(u2)*cos((pi*t)/2)*(u3))/2 
			+ (9*CST*pi*sin(pi*x*y)*sin(u1)*tanh(u2)*u4*sin((3*pi*t)/2))/8
	)
	
	return f

def powerDensityRing_spt(args:dict):
	time = args['time']
	position = args['position']
	nc_sp = np.size(position, axis=1); nc_tm = np.size(time); f = np.zeros((nc_sp, nc_tm))
	for i in range(nc_tm):
		t = time[i]
		f[:, i] = powerDensityRing_inc(args={'time':t, 'position':position})
	return np.ravel(f, order='F')

def simulate_incremental(degree, cuts, powerdensity, geoArgs=None, nbel_time=None, 
						quadArgs=None, solveSystem=True, is1dim=False, alpha=0.5):

	# Create geometry
	if quadArgs is None: quadArgs = {'quadrule':'iga', 'type':'leg'}
	if nbel_time is None: nbel_time = 2**CUTS_TIME
	dirichlet_table = np.ones((2, 2))

	if is1dim:
		geometry = createUniformOpenCurve(degree, int(2**cuts), 1.0)
		modelPhy = part1D(geometry, kwargs={'quadArgs':quadArgs})
	else:
		if geoArgs is None: geoArgs = {'name': 'SQ', 'degree': degree*np.ones(3, dtype=int), 
										'nb_refinementByDirection': np.array([cuts, 1, 1])}
		modelGeo = Geomdl(geoArgs)
		modelIGA = modelGeo.getIGAParametrization()
		modelPhy = part(modelIGA, quadArgs=quadArgs)

	time_inc = np.linspace(0, 1.0, nbel_time+1)

	# Add material 
	material = heatmat()
	material.addConductivity(conductivityProperty, isIsotropic=False) 
	material.addCapacity(capacityProperty, isIsotropic=False) 

	# Block boundaries
	boundary_inc = boundaryCondition(modelPhy.nbctrlpts)
	boundary_inc.add_DirichletConstTemperature(table=dirichlet_table)

	# Transient model
	if is1dim: problem_inc = heatproblem1D(material, modelPhy, boundary_inc)
	else: problem_inc = heatproblem(material, modelPhy, boundary_inc)
	Tinout = np.zeros((modelPhy.nbctrlpts_total, len(time_inc)))

	if not solveSystem: return problem_inc, time_inc, Tinout

	# Add external force 
	Fext_list = np.zeros((problem_inc.part.nbctrlpts_total, len(time_inc)))
	for i, t in enumerate(time_inc):
		Fext_list[:, i] = problem_inc.compute_volForce(powerdensity, 
														args={'position':problem_inc.part.qpPhy, 'time':t})

	# Solve
	problem_inc._thresNL = 1e-8; problem_inc._itersNL = 50
	problem_inc.solveFourierTransientProblem(Tinout=Tinout, Fext_list=Fext_list, 
											time_list=time_inc, alpha=alpha)
	return problem_inc, time_inc, Tinout

def simulate_spacetime(degree, cuts, powerdensity, geoArgs=None,  
					degree_time=None, nbel_time=None, quadArgs=None, solveSystem=True,
					isfull=False, isadaptive=True, getOthers=False, is1dim=False):
	
	if quadArgs is None: quadArgs = {'quadrule':'iga', 'type':'leg'}
	if nbel_time is None: nbel_time=2**CUTS_TIME
	if degree_time is None: degree_time = 2
	dirichlet_table = np.ones((3, 2)); dirichlet_table[-1, 1] = 0

	if is1dim:
		geometry = createUniformOpenCurve(degree, int(2**cuts), 1.0)
		modelPhy = part1D(geometry, kwargs={'quadArgs':quadArgs})
	else:
		if geoArgs is None: geoArgs = {'name': 'SQ', 'degree': degree*np.ones(3, dtype=int), 
										'nb_refinementByDirection': np.array([cuts, 1, 1])}
		modelGeo = Geomdl(geoArgs)
		modelIGA = modelGeo.getIGAParametrization()
		modelPhy = part(modelIGA, quadArgs=quadArgs)
	
	time_spt = part1D(createUniformOpenCurve(degree_time, nbel_time, 1.0), {'quadArgs':quadArgs})

	# Add material 
	material = heatmat()
	material.addConductivity(conductivityProperty, isIsotropic=False) 
	material.addCapacity(capacityProperty, isIsotropic=False)
	if isfull: material.addCapacityDers(capacityDersProperty, isIsotropic=False)
	if isfull: material.addConductivityDers(conductivityDersProperty, isIsotropic=False)

	# Block boundaries
	if is1dim: sptnbctrlpts = np.array([modelPhy.nbctrlpts_total, time_spt.nbctrlpts_total, 1])
	else: sptnbctrlpts = np.array([*modelPhy.nbctrlpts[:modelPhy.dim], time_spt.nbctrlpts_total])
	boundary_spt = boundaryCondition(sptnbctrlpts)
	boundary_spt.add_DirichletConstTemperature(table=dirichlet_table)

	# Transient model
	if is1dim: problem_spt = stheatproblem1D(material, modelPhy, time_spt, boundary_spt)
	else: problem_spt = stheatproblem(material, modelPhy, time_spt, boundary_spt)
	
	# Add external force
	Fext = problem_spt.compute_volForce(powerdensity, 
									{'position':problem_spt.part.qpPhy, 
									'time':problem_spt.time.qpPhy})
	
	Tinout = np.zeros(np.prod(sptnbctrlpts))
	if not solveSystem: return problem_spt, time_spt, Tinout

	problem_spt._thresNL = 1e-8; problem_spt._itersNL = 50
	output=problem_spt.solveFourierSTHeatProblem(Tinout=Tinout, Fext=Fext, isfull=isfull, isadaptive=isadaptive)
	if getOthers: return problem_spt, time_spt, Tinout, output
	return problem_spt, time_spt, Tinout