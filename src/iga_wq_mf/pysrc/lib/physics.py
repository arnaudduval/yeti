"""
.. module :: Physics
.. author :: Joaquin Cornejo
.. This module provides functions used in heat transfer equation
"""
import numpy as np

def power_density(P: list):
	" Compute power density at point P in physical space"
	f = 1 
	return f

def powden_cube(P: list):
	""" u = sin(pi*x)*sin(pi*y)*sin(pi*z)
		f = -div(lambda * grad(u))
	"""
	x = P[0, :]
	y = P[1, :]
	z = P[2, :]

	# # Isotropy
	# f = 3*np.pi**2*np.sin(np.pi*x)*np.sin(np.pi*y)*np.sin(np.pi*z) 

	# Anisotropy
	f = (6*np.pi**2*np.sin(np.pi*x)*np.sin(np.pi*y)*np.sin(np.pi*z) 
	- (np.pi**2*np.cos(np.pi*x)*np.cos(np.pi*z)*np.sin(np.pi*y))/5 
	- (np.pi**2*np.cos(np.pi*y)*np.cos(np.pi*z)*np.sin(np.pi*x))/2 
	- np.pi**2*np.cos(np.pi*x)*np.cos(np.pi*y)*np.sin(np.pi*z)
	)

	return f

def powden_prism(P: list):
	""" u = (-5*x+6*y+45)*(5*x+6*y-45)*x*(x-6)*sin(pi*z)
		f = -div(lambda * grad(u))
	"""
	x = P[0, :]
	y = P[1, :]
	z = P[2, :]

	# # Isotropy
	# f = (10*x*np.sin(np.pi*z)*(5*x + 6*y - 45) - 22*x*np.sin(np.pi*z)*(x - 6) 
	#     - 10*x*np.sin(np.pi*z)*(6*y - 5*x + 45) - 2*np.sin(np.pi*z)*(6*y - 5*x + 45)*(5*x + 6*y - 45) 
	#     - 10*np.sin(np.pi*z)*(x - 6)*(6*y - 5*x + 45) + 10*np.sin(np.pi*z)*(x - 6)*(5*x + 6*y - 45) 
	#     + x*np.pi**2*np.sin(np.pi*z)*(x - 6)*(6*y - 5*x + 45)*(5*x + 6*y - 45)
	# )/1000

	# Anisotropy
	f = (4*x*np.sin(np.pi*z)*(5*x + 6*y - 45) 
	- 94*x*np.sin(np.pi*z)*(x - 6) 
	- 16*x*np.sin(np.pi*z)*(6*y - 5*x + 45) 
	- 2*np.sin(np.pi*z)*(6*y - 5*x + 45)*(5*x + 6*y - 45) 
	- 16*np.sin(np.pi*z)*(x - 6)*(6*y - 5*x + 45) 
	+ 4*np.sin(np.pi*z)*(x - 6)*(5*x + 6*y - 45) 
	- (np.pi*np.cos(np.pi*z)*(x - 6)*(6*y - 5*x + 45)*(5*x + 6*y - 45))/5 
	- 4*x*np.pi*np.cos(np.pi*z)*(x - 6)*(6*y - 5*x + 45) 
	- 2*x*np.pi*np.cos(np.pi*z)*(x - 6)*(5*x + 6*y - 45) 
	- (x*np.pi*np.cos(np.pi*z)*(6*y - 5*x + 45)*(5*x + 6*y - 45))/5 
	+ 3*x*np.pi**2*np.sin(np.pi*z)*(x - 6)*(6*y - 5*x + 45)*(5*x + 6*y - 45)
	)/1000
	
	return f

def powden_thickring(P: list):
	""" u = sin(5*pi*x)*sin(5*pi*y)*sin(5*pi*z)*(x**2+y**2-1)*(x**2+y**2-4)
		f = -div(lambda * grad(u))
	"""
	x = P[0, :]
	y = P[1, :]
	z = P[2, :] 

	# # Isotropy
	# f = (75*np.pi**2*np.sin(5*np.pi*x)*np.sin(5*np.pi*y)*np.sin(5*np.pi*z)*(x**2 + y**2 - 1)*(x**2 + y**2 - 4) 
	#     - 8*y**2*np.sin(5*np.pi*x)*np.sin(5*np.pi*y)*np.sin(5*np.pi*z) 
	#     - 4*np.sin(5*np.pi*x)*np.sin(5*np.pi*y)*np.sin(5*np.pi*z)*(x**2 + y**2 - 1) 
	#     - 4*np.sin(5*np.pi*x)*np.sin(5*np.pi*y)*np.sin(5*np.pi*z)*(x**2 + y**2 - 4) 
	#     - 8*x**2*np.sin(5*np.pi*x)*np.sin(5*np.pi*y)*np.sin(5*np.pi*z) 
	#     - 20*x*np.pi*np.cos(5*np.pi*x)*np.sin(5*np.pi*y)*np.sin(5*np.pi*z)*(x**2 + y**2 - 1) 
	#     - 20*x*np.pi*np.cos(5*np.pi*x)*np.sin(5*np.pi*y)*np.sin(5*np.pi*z)*(x**2 + y**2 - 4) 
	#     - 20*y*np.pi*np.cos(5*np.pi*y)*np.sin(5*np.pi*x)*np.sin(5*np.pi*z)*(x**2 + y**2 - 1) 
	#     - 20*y*np.pi*np.cos(5*np.pi*y)*np.sin(5*np.pi*x)*np.sin(5*np.pi*z)*(x**2 + y**2 - 4)
	# )

	# Anisotropy
	f = (150*np.pi**2*np.sin(5*np.pi*x)*np.sin(5*np.pi*y)*np.sin(5*np.pi*z)*(x**2 + y**2 - 1)*(x**2 + y**2 - 4)
	- 16*y**2*np.sin(5*np.pi*x)*np.sin(5*np.pi*y)*np.sin(5*np.pi*z) 
	- 6*np.sin(5*np.pi*x)*np.sin(5*np.pi*y)*np.sin(5*np.pi*z)*(x**2 + y**2 - 1) 
	- 6*np.sin(5*np.pi*x)*np.sin(5*np.pi*y)*np.sin(5*np.pi*z)*(x**2 + y**2 - 4) 
	- 8*x*y*np.sin(5*np.pi*x)*np.sin(5*np.pi*y)*np.sin(5*np.pi*z) 
	- 8*x**2*np.sin(5*np.pi*x)*np.sin(5*np.pi*y)*np.sin(5*np.pi*z) 
	- 25*np.pi**2*np.cos(5*np.pi*x)*np.cos(5*np.pi*y)*np.sin(5*np.pi*z)*(x**2 + y**2 - 1)*(x**2 + y**2 - 4) 
	- 5*np.pi**2*np.cos(5*np.pi*x)*np.cos(5*np.pi*z)*np.sin(5*np.pi*y)*(x**2 + y**2 - 1)*(x**2 + y**2 - 4) 
	- (25*np.pi**2*np.cos(5*np.pi*y)*np.cos(5*np.pi*z)*np.sin(5*np.pi*x)*(x**2 + y**2 - 1)*(x**2 + y**2 - 4))/2 
	- 20*x*np.pi*np.cos(5*np.pi*x)*np.sin(5*np.pi*y)*np.sin(5*np.pi*z)*(x**2 + y**2 - 1) 
	- 10*x*np.pi*np.cos(5*np.pi*y)*np.sin(5*np.pi*x)*np.sin(5*np.pi*z)*(x**2 + y**2 - 1) 
	- 2*x*np.pi*np.cos(5*np.pi*z)*np.sin(5*np.pi*x)*np.sin(5*np.pi*y)*(x**2 + y**2 - 1) 
	- 20*x*np.pi*np.cos(5*np.pi*x)*np.sin(5*np.pi*y)*np.sin(5*np.pi*z)*(x**2 + y**2 - 4) 
	- 10*x*np.pi*np.cos(5*np.pi*y)*np.sin(5*np.pi*x)*np.sin(5*np.pi*z)*(x**2 + y**2 - 4) 
	- 2*x*np.pi*np.cos(5*np.pi*z)*np.sin(5*np.pi*x)*np.sin(5*np.pi*y)*(x**2 + y**2 - 4) 
	- 10*y*np.pi*np.cos(5*np.pi*x)*np.sin(5*np.pi*y)*np.sin(5*np.pi*z)*(x**2 + y**2 - 1) 
	- 40*y*np.pi*np.cos(5*np.pi*y)*np.sin(5*np.pi*x)*np.sin(5*np.pi*z)*(x**2 + y**2 - 1) 
	- 5*y*np.pi*np.cos(5*np.pi*z)*np.sin(5*np.pi*x)*np.sin(5*np.pi*y)*(x**2 + y**2 - 1) 
	- 10*y*np.pi*np.cos(5*np.pi*x)*np.sin(5*np.pi*y)*np.sin(5*np.pi*z)*(x**2 + y**2 - 4) 
	- 40*y*np.pi*np.cos(5*np.pi*y)*np.sin(5*np.pi*x)*np.sin(5*np.pi*z)*(x**2 + y**2 - 4) 
	- 5*y*np.pi*np.cos(5*np.pi*z)*np.sin(5*np.pi*x)*np.sin(5*np.pi*y)*(x**2 + y**2 - 4)
	)

	return f

def powden_rotring(P: list):
	""" u = -(x**2 + y**2 - 1)*(x**2 + y**2 - 4)*x*(y**2)*sin(pi*z)
		f = -div(lambda * grad(u))
	"""
	x = P[0, :]
	y = P[1, :]
	z = P[2, :]

	# # Isotropy
	# f = (8*x*y**4*np.sin(np.pi*z) + 8*x**3*y**2*np.sin(np.pi*z) 
	#     + 16*x*y**2*np.sin(np.pi*z)*(x**2 + y**2 - 1) + 16*x*y**2*np.sin(np.pi*z)*(x**2 + y**2 - 4) 
	#     + 2*x*np.sin(np.pi*z)*(x**2 + y**2 - 1)*(x**2 + y**2 - 4) 
	#     - x*y**2*np.pi**2*np.sin(np.pi*z)*(x**2 + y**2 - 1)*(x**2 + y**2 - 4)
	# )

	# Anisotropy
	f = (16*x*y**4*np.sin(np.pi*z) 
	+ 8*x**2*y**3*np.sin(np.pi*z) 
	+ 8*x**3*y**2*np.sin(np.pi*z) 
	+ 2*y**3*np.sin(np.pi*z)*(x**2 + y**2 - 1) 
	+ 2*y**3*np.sin(np.pi*z)*(x**2 + y**2 - 4) 
	+ 26*x*y**2*np.sin(np.pi*z)*(x**2 + y**2 - 1) 
	+ 4*x**2*y*np.sin(np.pi*z)*(x**2 + y**2 - 1) 
	+ 26*x*y**2*np.sin(np.pi*z)*(x**2 + y**2 - 4) 
	+ 4*x**2*y*np.sin(np.pi*z)*(x**2 + y**2 - 4) 
	+ 4*x*np.sin(np.pi*z)*(x**2 + y**2 - 1)*(x**2 + y**2 - 4) 
	+ 2*y*np.sin(np.pi*z)*(x**2 + y**2 - 1)*(x**2 + y**2 - 4) 
	+ (y**2*np.pi*np.cos(np.pi*z)*(x**2 + y**2 - 1)*(x**2 + y**2 - 4))/5 
	+ x*y**3*np.pi*np.cos(np.pi*z)*(x**2 + y**2 - 1) 
	+ x*y**3*np.pi*np.cos(np.pi*z)*(x**2 + y**2 - 4) 
	+ (2*x**2*y**2*np.pi*np.cos(np.pi*z)*(x**2 + y**2 - 1))/5 
	+ (2*x**2*y**2*np.pi*np.cos(np.pi*z)*(x**2 + y**2 - 4))/5 
	- 3*x*y**2*np.pi**2*np.sin(np.pi*z)*(x**2 + y**2 - 1)*(x**2 + y**2 - 4) 
	+ x*y*np.pi*np.cos(np.pi*z)*(x**2 + y**2 - 1)*(x**2 + y**2 - 4)
	)

	return f

def temperature_rotring(P: list):
	" T = -(x**2 + y**2 - 1)*(x**2 + y**2 - 4)*x*(y**2)*sin(pi*z) "
	x = P[0, :]
	y = P[1, :]
	z = P[2, :]
	f = -(x**2 + y**2 - 1)*(x**2 + y**2 - 4)*x*(y**2)*np.sin(np.pi*z)

	return f

def powden_annulus(P: list):
	""" u = (x**2 + y**2 - 1)*(x**2 + y**2 - 4)*sin(pi*x)*sin(pi*y)
		f = -div(lambda * grad(u))
	"""
	x = P[0, :]
	y = P[1, :]

	# # Isotropy
	# f = (2*np.pi**2*np.sin(np.pi*x)*np.sin(np.pi*y)*(x**2 + y**2 - 1)*(x**2 + y**2 - 4) 
	#     - 8*y**2*np.sin(np.pi*x)*np.sin(np.pi*y) 
	#     - 4*np.sin(np.pi*x)*np.sin(np.pi*y)*(x**2 + y**2 - 1) 
	#     - 4*np.sin(np.pi*x)*np.sin(np.pi*y)*(x**2 + y**2 - 4) 
	#     - 4*x*np.pi*np.cos(np.pi*x)*np.sin(np.pi*y)*(x**2 + y**2 - 1) 
	#     - 4*x*np.pi*np.cos(np.pi*x)*np.sin(np.pi*y)*(x**2 + y**2 - 4) 
	#     - 4*y*np.pi*np.cos(np.pi*y)*np.sin(np.pi*x)*(x**2 + y**2 - 1) 
	#     - 4*y*np.pi*np.cos(np.pi*y)*np.sin(np.np.pi*x)*(x**2 + y**2 - 4) 
	#     - 8*x**2*np.sin(np.pi*x)*np.sin(np.pi*y)
	# )

	# Anisotropy lambda = [1, 0; 0, 0.1]
	f = ((11*np.pi**2*np.sin(np.pi*x)*np.sin(np.pi*y)*(x**2 + y**2 - 1)*(x**2 + y**2 - 4))/10
		- (4*y**2*np.sin(np.pi*x)*np.sin(np.pi*y))/5
		- (11*np.sin(np.pi*x)*np.sin(np.pi*y)*(x**2 + y**2 - 1))/5
		- (11*np.sin(np.pi*x)*np.sin(np.pi*y)*(x**2 + y**2 - 4))/5
		- 4*x*np.pi*np.cos(np.pi*x)*np.sin(np.pi*y)*(x**2 + y**2 - 1)
		- 4*x*np.pi*np.cos(np.pi*x)*np.sin(np.pi*y)*(x**2 + y**2 - 4)
		- (2*y*np.pi*np.cos(np.pi*y)*np.sin(np.pi*x)*(x**2 + y**2 - 1))/5
		- (2*y*np.pi*np.cos(np.pi*y)*np.sin(np.pi*x)*(x**2 + y**2 - 4))/5
		- 8*x**2*np.sin(np.pi*x)*np.sin(np.pi*y)
	)

	return f

# ------------------

def setKprop(T, prop=0.1):
	# y = prop + prop*np.exp(-0.1*abs(T))
	# y = prop + prop*2.0/(1.0 + np.exp(-5*(T-1.0)))
	y   = np.ones(len(T))
	return y

def setCprop(T, prop=1.0):
	# y = prop + prop*np.exp(-2.0*abs(T))
	y   = np.ones(len(T))
	return y

def powden(P:list, dim=1):
	if dim == 1  : x = P[:]
	elif dim == 2: x = P[0, :]
	elif dim == 3: x = P[0, :]
	f = 0.0*np.sin(np.pi*x)
	return f

# ------------------

def forceVol(P:list):
	# force = 5e3*np.sin(np.pi*P)
	force = 2e3*P
	# force = 1e2*np.ones(len(P))
	return force
