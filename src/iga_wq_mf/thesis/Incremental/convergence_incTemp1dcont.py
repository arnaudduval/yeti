from thesis.Incremental.__init__ import *
from pysrc.lib.lib_job1d import stheatproblem1D
from scipy.integrate import solve_ivp
from pysrc.lib.lib_base import bdf
from numpy import sin, cos, pi, tanh

CST = 100
CUTS_TIME = 6
ISLINEAR = False

def capacityProperty(args:dict):
	temperature = args.get('temperature')
	capacity = np.ones(shape=np.shape(temperature))
	return capacity

def conductivityProperty(args:dict):
	temperature = args.get('temperature')
	if ISLINEAR: 
		conductivity = 2*np.ones(shape=np.shape(temperature))
	else:
		conductivity = 3.0 + 2.0*tanh(temperature/50)
	return conductivity

def exactTemperature_inc(args:dict):
	t = args['time']
	x = args['position']
	u = CST*sin(2*pi*x)*sin(pi/2*t)*(1+0.75*cos(3*pi/2*t)) 
	return u

def powerDensity_inc(args:dict):
	t = args['time']
	x = args['position']

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

def simulate_incremental(degree, cuts, powerdensity, dirichlet_table=None, nbel_time=None, 
						quadArgs=None, IVPmethod='alpha', alpha=0.5):

	# Create geometry
	if quadArgs is None: quadArgs = {'quadrule':'iga', 'type':'leg'}
	if nbel_time is None: nbel_time = 2**CUTS_TIME
	if dirichlet_table is None: dirichlet_table = np.zeros((2, 2)); dirichlet_table[0, :] = 1

	geometry = createUniformOpenCurve(degree, int(2**cuts), 1.0)
	modelPhy = part1D(geometry, kwargs={'quadArgs':quadArgs})

	time_inc = np.linspace(0, 1.0, nbel_time+1)

	# Add material 
	material = heatmat()
	material.addConductivity(conductivityProperty, isIsotropic=False) 
	material.addCapacity(capacityProperty, isIsotropic=False) 

	# Block boundaries
	boundary_inc = boundaryCondition(modelPhy.nbctrlpts)
	boundary_inc.add_DirichletConstTemperature(table=dirichlet_table)

	# Transient model
	problem_inc = heatproblem1D(material, modelPhy, boundary_inc)
	Tinout = np.zeros((modelPhy.nbctrlpts_total, len(time_inc)))

	if IVPmethod == 'alpha':
		# Add external force 
		Fext_list = np.zeros((problem_inc.part.nbctrlpts_total, len(time_inc)))
		for i, t in enumerate(time_inc):
			Fext_list[:, i] = problem_inc.compute_volForce(powerdensity, 
								args={'position':problem_inc.part.qpPhy, 'time':t})

		# Solve
		problem_inc._itersNL = 50; problem_inc._thresNL = 1e-8
		problem_inc.solveFourierTransientProblem(Tinout=Tinout, Fext_list=Fext_list, 
												time_list=time_inc, alpha=alpha)
	
	else:
		args = {'temperature':np.ones(problem_inc.part.nbqp), 'position':problem_inc.part.qpPhy}
		dof = problem_inc.boundary.thdof
		cprop = problem_inc.heatmaterial.capacity(args)*problem_inc.heatmaterial.density(args)*problem_inc.part.detJ
		capacity = problem_inc.part._denseweights[0] @ np.diag(cprop) @ problem_inc.part._densebasis[0].T
		invcapacity = np.linalg.pinv(capacity[np.ix_(dof, dof)])

		def fun(t, y):
			"y is of size ndof"
			tmp = np.zeros(problem_inc.part.nbctrlpts_total); tmp[dof] = y
			tmpinterp = problem_inc.interpolate_temperature(tmp)
			args['time'] = t; args['temperature'] = tmpinterp
			force = problem_inc.compute_volForce(powerdensity, args=args)[dof]
			MFKu = problem_inc.compute_mfConductivity(tmp, args=args)[dof]
			vel = invcapacity @ (force - MFKu)
			return vel

		t_span = (time_inc[0], time_inc[-1])
		y0 = np.zeros(problem_inc.boundary._thndof)

		if IVPmethod in ['BDF1', 'BDF2', 'BDF3', 'BDF4']:
			y = np.transpose(bdf(fun, t_span, y0, nbel_time, norder=int(IVPmethod[-1]))[-1])
			Tinout[dof, :] = np.copy(y)

		else:
			print("SCIPY METHODS: WE DO NOT CONTROL TIME STEP")
			output = solve_ivp(fun, t_span, y0, 
							method=IVPmethod, 
							t_eval=time_inc, 
							first_step=time_inc[1]-time_inc[0],
							max_step=(time_inc[1]-time_inc[0]),
							rtol=1e-8, atol=1e-10)
			Tinout[dof, :] = np.copy(output['y'])

	return problem_inc, time_inc, Tinout

# Set global variables
SUFIX = 'lin_1d' if ISLINEAR else 'nonlin_1d'
PLOTRELATIVE = True
RUNSIMU = True
EXTENSION = '.dat'

filenameA1 = FOLDER2DATA + '1incheatAbs'+SUFIX
filenameR1 = FOLDER2DATA + '1incheatRel'+SUFIX
filenameT1 = FOLDER2DATA + '1incheatTim'+SUFIX

degList = np.array([1, 2, 3, 4])
cutList = np.arange(1, 8)
cutstime = 5

if RUNSIMU:

	A1errorList = np.ones((len(degList), len(cutList)))
	R1errorList = np.ones((len(degList), len(cutList)))
	T1timeList = np.ones((len(degList), len(cutList)))

	A2errorList = np.ones((len(degList), len(cutList)))
	R2errorList = np.ones((len(degList), len(cutList)))
	T2timeList = np.ones((len(degList), len(cutList)))
	
	for j, cuts in enumerate(cutList):
		for i, degree in enumerate(degList):

			start = time.process_time()
			dirichlet_table = np.ones((2, 2))
			problem_inc, time_inc, temp_inc = simulate_incremental(degree, cuts, powerDensity_inc, 
												dirichlet_table=dirichlet_table, nbel_time=2**cutstime)
			finish = time.process_time()
			T1timeList[i, j] = finish - start
			
			# Error of last step
			A1errorList[i, j], R1errorList[i, j] = problem_inc.normOfError(temp_inc[:, -1], 
													normArgs={'type':'L2',
																'exactFunction':exactTemperature_inc,
																'exactExtraArgs':{'time':time_inc[-1]}})

		np.savetxt(filenameA1+EXTENSION, A1errorList)
		np.savetxt(filenameR1+EXTENSION, R1errorList)
		np.savetxt(filenameT1+EXTENSION, T1timeList)

errorList1 = np.loadtxt(filenameR1+EXTENSION)
fig, ax = plt.subplots()
for i, degree in enumerate(degList):
	color = COLORLIST[i]
	ax.loglog(2**cutList, errorList1[i, :], color=color, marker=CONFIGLINE4['marker'], markerfacecolor='w',
				markersize=CONFIGLINE4['markersize'], linestyle=CONFIGLINE4['linestyle'], label='IGA-GL deg. '+str(degree))

ax.set_ylabel('Relative ' + r'$L^2$' + ' error at last time-step')
ax.set_xlabel('Number of elements')
ax.set_xlim(left=1, right=200)
ax.set_ylim(top=1e1, bottom=1e-6)
ax.legend(loc='upper right')
fig.tight_layout()
fig.savefig(FOLDER2SAVE + 'INCL1Convergence' + str(cutstime)  + SUFIX +  '.pdf')
plt.close(fig)
