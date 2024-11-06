"""
... My original code does not include NURBS
... In this file, I implement some ideas to eventually include them into YETI
"""

from pysrc.lib.__init__ import *
from pysrc.lib.lib_part import part
from pysrc.lib.lib_boundary import boundaryCondition
from pysrc.lib.lib_job3d import problem
from pysrc.lib.lib_quadrules import GaussQuadrature
from pysrc.lib.lib_base import vtk2png

FOLDER = os.path.dirname(os.path.realpath(__file__)) + '/'

def exactfunction(P:list):
	x = P[0, :]; y = P[1, :]
	F = np.sin(np.pi*x)*np.sin(np.pi*y)*np.sin(np.pi*(x+4)*(y-4))*(x**2+y**2-1)
	return F

def normOfError(l2problem:problem, u_ctrlpts, normArgs:dict):
	""" Computes the norm L2 or H1 of the error. The exactfun is the function of the exact solution. 
		and u_ctrlpts is the field at the control points. We compute the integral using Gauss Quadrature
		whether the default quadrature is weighted quadrature. 
	"""	
	typeNorm = normArgs.get('type', 'l2').lower()
	if all(norm != typeNorm  for norm in ['l2']): raise Warning('Unknown norm')

	# Compute u interpolation
	nbqp, quadPts, indices, basis, parametricWeights = [], [], [], [], []
	for i in range(l2problem.part.dim):
		quadRule = GaussQuadrature(l2problem.part.degree[i], l2problem.part.knotvector[i], quadArgs={'type':'leg'})
		quadPtsByDir, indicesByDir, basisByDir, _ = quadRule.getQuadratureRulesInfo()
		indi, indj = indicesByDir; parweightsByDir = quadRule._parametricWeights
		
		nbqp.append(quadRule.nbqp); quadPts.append(quadPtsByDir); indices.append(indi); indices.append(indj)
		basis.append(basisByDir); parametricWeights.append(parweightsByDir)

	inpts = [*nbqp, *indices, *basis, np.atleast_2d(l2problem.part.NURBSwgt)]
	if l2problem.part.dim == 2:
		NURBScorrection = geophy.interpolate_meshgrid_2d(*inpts)[0, :]
		dersNURBScorrection = geophy.eval_jacobien_2d(*inpts)[0, :, :]

	elif l2problem.part.dim == 3: 
		NURBScorrection = geophy.interpolate_meshgrid_3d(*inpts)
		dersNURBScorrection = geophy.eval_jacobien_3d(*inpts)

	inpts = [*nbqp, *indices, *basis, l2problem.part.ctrlpts]
	if l2problem.part.dim == 2:
		qpPhy = geophy.interpolate_meshgrid_2d(*inpts)
		pseudojac = geophy.eval_jacobien_2d(*inpts)
	if l2problem.part.dim == 3:
		qpPhy = geophy.interpolate_meshgrid_3d(*inpts)
		pseudojac = geophy.eval_jacobien_3d(*inpts)
	for i in range(l2problem.part.dim):
		qpPhy[i, :] /= NURBScorrection
		
	Jqp = np.zeros((l2problem.part.dim, l2problem.part.dim, len(NURBScorrection)))
	for i in range(l2problem.part.dim):
		for j in range(l2problem.part.dim):
			Jqp[i, j, :] = (pseudojac[i, j, :] - qpPhy[i, :]*dersNURBScorrection[j, :])/NURBScorrection

	detJ = geophy.eval_inverse_det(Jqp)[0]

	inpts = [*nbqp, *indices, *basis, np.atleast_2d(u_ctrlpts)]
	if l2problem.part.dim == 2:   u_interp = geophy.interpolate_meshgrid_2d(*inpts)    
	elif l2problem.part.dim == 3: u_interp = geophy.interpolate_meshgrid_3d(*inpts)
	for i in range(np.size(u_interp, axis=0)):
		u_interp[i, :] /= NURBScorrection

	# --------------------------------
	# Compute u exact
	u_exact = None
	exactfun = normArgs.get('exactFunction', None)
	exactextraArgs = normArgs.get('exactExtraArgs', None)
	if exactextraArgs is not None:
		assert isinstance(exactextraArgs, dict), 'Error type of extra args'
		if not 'position' in exactextraArgs.keys(): exactextraArgs['position'] = qpPhy
		if callable(exactfun): u_exact = np.atleast_2d(exactfun(exactextraArgs))
	else:
		if callable(exactfun): u_exact = np.atleast_2d(exactfun(qpPhy))

	# Compute error
	uedfuf2_l2 = np.einsum('il->l', (u_exact - u_interp)**2)
	ue2_l2     = np.einsum('il->l', u_exact**2)
	norm1 = uedfuf2_l2*detJ
	norm2 = ue2_l2*detJ

	norm1 = np.reshape(norm1, tuple(nbqp), order='F')
	norm2 = np.reshape(norm2, tuple(nbqp), order='F')
	if l2problem.part.dim == 2: 
		tmp1 = np.einsum('i,j,ij->', parametricWeights[0], parametricWeights[1], norm1)
		tmp2 = np.einsum('i,j,ij->', parametricWeights[0], parametricWeights[1], norm2)
	if l2problem.part.dim == 3: 
		tmp1 = np.einsum('i,j,k,ijk->', parametricWeights[0], parametricWeights[1], parametricWeights[2], norm1)
		tmp2 = np.einsum('i,j,k,ijk->', parametricWeights[0], parametricWeights[1], parametricWeights[2], norm2)
		
	abserror = np.sqrt(tmp1)
	relerror = np.sqrt(tmp1/tmp2) if tmp2!=0 else np.sqrt(tmp1)
	if tmp2 == 0: print('Warning: Dividing by zero')

	return abserror, relerror

def simulate_l2problem(degree, cuts, quadArgs, plot=False):
	blockPrint()
	# Creation of the IGA object
	modelIGA = IGAparametrization(filename=FOLDER+'platewithhole')
	modelIGA.refine(nb_refinementByDirection=np.array([cuts, cuts+1, 1]), 
					nb_degreeElevationByDirection=degree*np.array([1, 1, 1]),)

	# Definition of the L2 projection problem
	modelPhy = part(modelIGA, quadArgs=quadArgs)
	boundary = boundaryCondition()
	l2problem = problem(modelPhy, boundary, {})
	l2problem._thresLin = 1e-14

	u_atQuadpts = exactfunction(modelPhy.qpPhy)/modelPhy.NURBScorrection
	additional_property = 1/modelPhy.NURBScorrection**2
	u_atCtrlpts = l2problem.L2projectionCtrlpts(u_atQuadpts, table=np.ones((2, 2)), prop=additional_property)
	abserror, relerror = normOfError(l2problem, u_atCtrlpts, normArgs={'exactFunction':exactfunction})
	if plot: modelPhy.postProcessingPrimal(fields={'approxfield':u_atCtrlpts}, folder=FOLDER, name='l2problem', sampleSize=201)
	enablePrint()
	return abserror, relerror

degList = np.arange(0, 5)
cutList = np.arange(0, 8)
RUNSIMU = True

if RUNSIMU:

	AbserrorList = np.ones((len(degList), len(cutList)))
	RelerrorList = np.copy(AbserrorList)

	for quadrule, quadtype in zip(['iga', 'wq'], ['leg', 3]):
		quadArgs = {'quadrule': quadrule, 'type': quadtype}
		
		for i, degree in enumerate(degList):
			for j, cuts in enumerate(cutList):
				AbserrorList[i, j], RelerrorList[i, j] = simulate_l2problem(degree, cuts, quadArgs)
			np.savetxt(FOLDER+'AbsError_l2_'+quadrule+'_'+str(quadtype)+'.dat', AbserrorList)
			np.savetxt(FOLDER+'RelError_l2_'+quadrule+'_'+str(quadtype)+'.dat', RelerrorList)

simulate_l2problem(6, 6, {'quadrule':'iga', 'type':'leg'}, plot=True)
vtk2png(folder=FOLDER, filename='l2problem', fieldname='approxfield', cmap='coolwarm', title='Exact function', position_y=0.1)

fig, ax = plt.subplots()
figname = FOLDER + 'ConvergenceL2_projection'
for quadrule, quadtype, plotpars in zip(['iga', 'wq'], ['leg', 3], [CONFIGLINE0, CONFIGLINE2]):
	quadArgs = {'quadrule': quadrule, 'type': quadtype}
	errorList = np.loadtxt(FOLDER+'RelError_l2_'+quadrule+'_'+str(quadtype)+'.dat')

	for i, degree in enumerate(degList):
		color = COLORLIST[i+1]
		nbelList = 2*2**cutList

		if quadrule == 'iga': 
			ax.loglog(nbelList, errorList[i, :], label='IGA-GL deg. '+str(degree+2), color=color, marker=plotpars['marker'], 
					markerfacecolor='w', markersize=plotpars['markersize'], linestyle=plotpars['linestyle'])
			# slope = round(np.polyfit(np.log(nbelList), np.log(errorList[i, :]), 1)[0], 1)
			# annotation.slope_marker((nbelList[-2],  errorList[i, -2]), slope, 
			# 				poly_kwargs={'facecolor': (0.73, 0.8, 1)}, ax=ax)
		else: 
			ax.loglog(nbelList, errorList[i, :], color=color, marker=plotpars['marker'], markerfacecolor='w',
				markersize=plotpars['markersize'], linestyle=plotpars['linestyle'])

ax.loglog([], [], color='k', marker=CONFIGLINE2['marker'], markerfacecolor='w',
		markersize=CONFIGLINE2['markersize'], linestyle=CONFIGLINE2['linestyle'], label='IGA-WQ 2')

ax.set_ylabel('Relative '+r'$L^2$'+' error')
ax.set_xlabel('Number of elements by dimension')
ax.set_ylim(top=1e2, bottom=1e-14)
ax.set_xlim(left=1, right=400)
ax.legend(loc='lower left')
fig.tight_layout()
fig.savefig(figname+'.png')
