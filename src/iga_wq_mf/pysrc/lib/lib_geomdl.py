"""
.. This module contains functions to construct specific geometries using only B-splines
.. Joaquin Cornejo
"""

from lib.__init__ import *
from lib.lib_base import createUniformMaxregularKnotvector

class Geomdl():
	def __init__(self, geoArgs:dict):
		self._dim     = None
		self._name    = geoArgs.get('name', '').lower()
		self._degree  = geoArgs.get('degree', None)
		self._nbcuts  = geoArgs.get('nb_refinementByDirection', np.array([1, 1, 1]))
		self._extraArgs = geoArgs.get('extra', {})
		return
	
	def __getInfo(self, obj): 
		" Updates and saves important properties of the geometry created "

		info  = {}
		dimen = self._dim
		if obj is None: raise Warning('Geometry unknown')

		# Set degree
		degree = np.ones(3, dtype= int)
		degree[:dimen] = np.array(obj._degree)
		if any(p == 1 for p in degree[:dimen]): 
			raise Warning('Model must have at least degree 2')
		info['degree'] = degree

		# Set knot vector
		knotvector = [np.array(obj._knot_vector[dim]) for dim in range(dimen)]
		info['knotvector'] = knotvector

		# Set size knot-vector in each dimension
		size_kv = np.ones(3, dtype= int)  
		for i in range(dimen):
			size_kv[i] = np.size(knotvector[i])

		# Set number of quadrature points in each dimension
		nbctrlpts = np.ones(3, dtype= int)  
		for i in range(dimen):
			nbctrlpts[i] = size_kv[i] - degree[i] - 1
		nbctrlpts_total = np.product(nbctrlpts)
		info['nbctrlpts'] = nbctrlpts

		if dimen == 2: 
			c = 0
			ctrlpts_old = obj._control_points
			ctrlpts_new = np.zeros((3, nbctrlpts_total))
			for j in range(nbctrlpts[1]):
				for i in range(nbctrlpts[0]):
					pos = j + i*nbctrlpts[1]
					ctrlpts_temp = ctrlpts_old[pos]
					ctrlpts_new[:, c] = ctrlpts_temp
					c += 1 

		elif dimen == 3: 
			c =  0
			ctrlpts_old = obj._control_points
			ctrlpts_new = np.zeros((3, nbctrlpts_total))
			for k in range(nbctrlpts[2]):
				for j in range(nbctrlpts[1]):
					for i in range(nbctrlpts[0]):
						pos = j + i*nbctrlpts[1] + k*nbctrlpts[1]*nbctrlpts[0]
						ctrlpts_temp = ctrlpts_old[pos]
						ctrlpts_new[:, c] = ctrlpts_temp
						c += 1 
		info['ctrlpts'] = ctrlpts_new
		return info

	def __writeAbaqusFile(self, info):
		" Writes an inp and NB file. By the moment, it only works with one patch"

		def array2txt(array: np.array, format= '%.2f'):
			return ','.join([format %(i) for i in array])
		
		dimen 	   = self._dim
		name       = info['name']
		degree     = info['degree']
		knotvector = info['knotvector']
		nbctrlpts  = info['nbctrlpts']
		ctrlpts    = info['ctrlpts']
		nbctrlpts_total = np.product(nbctrlpts)

		# .inp file
		inpfile = name + '.inp'
		introduction =  [
			'** Copyright 2020 Thibaut Hirschler',
			'** Copyright 2020 Arnaud Duval',
			'** This file is part of Yeti.',
			'**',
			'** Yeti is free software: you can redistribute it and/or modify it under the terms',
			'** of the GNU Lesser General Public License as published by the Free Software',
			'** Foundation, either version 3 of the License, or (at your option) any later version.',
			'**',
			'** Yeti is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;',
			'** without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR',
			'** PURPOSE. See the GNU Lesser General Public License for more details.',
			'**',
			'** You should have received a copy of the GNU Lesser General Public License along',
			'** with Yeti. If not, see <https://www.gnu.org/licenses/>',
			'**', 
			'*HEADING',
			'**NurbsABQ - Laboratoire de Mecanique des Contacts et des Solides - INSA-Lyon'
			]
		with open(inpfile, 'w') as f:
			f.write('\n'.join(introduction))
			f.write('\n')
			f.write('*Part, name=%s\n' %name)
			f.write('*USER ELEMENT, NODES=%d, TYPE=U1, COORDINATES=%d, INTEGRATION=%d\n' 
					%(nbctrlpts_total, dimen, nbctrlpts_total))
			f.write(array2txt(np.arange(dimen)+1, format='%d'))
			f.write('\n*Node,nset=AllNode\n')
			for i in range(nbctrlpts_total):
				CP = ctrlpts[:, i]
				f.write('%d, %.15f, %.15f, %.15f\n' %(i+1, CP[0], CP[1], CP[2]))
			f.write('*Element,type=U1,elset=AllEls\n1,\t')
			f.write(array2txt(np.arange(nbctrlpts_total, 0, -1), format='%d'))
			f.write('\n')
			f.write('*ELSET,ELSET=EltPatch1,generate\n1,1,1\n')
			f.write('*UEL PROPERTY, ELSET=EltPatch1, MATERIAL=Mat\n1\n')
			f.write('*End Part\n')
			f.write('**ASSEMBLY\n*Assembly, name=Assembly\n')
			f.write('*Instance, name=I1, part=%s\n' %name)
			f.write('*End Instance\n*End Assembly\n')
			f.write('**MATERIAL\n*MATERIAL,NAME=Mat\n*Elastic\n')
			f.write('%f, %f\n' %(3e3, 0.3)) 
			f.write('*STEP,extrapolation=NO,NLGEOM=NO\n*Static\n')
			f.write('** OUTPUT REQUESTS\n*node file,frequency=1\nU,RF,CF\n*el file,frequency=1\nSDV\n*End Step')

		# .NB file
		NBfile = name + '.NB'
		introduction =  [
			'** Copyright 2020 Thibaut Hirschler',
			'** Copyright 2020 Arnaud Duval',
			'** This file is part of Yeti.',
			'**',
			'** Yeti is free software: you can redistribute it and/or modify it under the terms',
			'** of the GNU Lesser General Public License as published by the Free Software',
			'** Foundation, either version 3 of the License, or (at your option) any later version.',
			'**',
			'** Yeti is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;',
			'** without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR',
			'** PURPOSE. See the GNU Lesser General Public License for more details.',
			'**',
			'** You should have received a copy of the GNU Lesser General Public License along',
			'** with Yeti. If not, see <https://www.gnu.org/licenses/>'
			]
		with open(NBfile, 'w') as f:
			f.write('\n'.join(introduction))
			f.write('\n\n')
			f.write('*Dimension\n%d\n' %(dimen))
			f.write('*Number of CP by element\n%d\n' %(nbctrlpts_total))
			f.write('*Number of patch\n%d\n' %(1))
			f.write('*Total number of element \n%d\n' %(1))
			f.write('*Number of element by patch\n%d\n' %(1))
			f.write('*Patch(1)\n')
			for i in range(dimen):
				kv = knotvector[i]
				f.write('%d\n' %(len(kv)))
				f.write(array2txt(kv))
				f.write('\n')
			f.write('*Jpqr\n')
			f.write(array2txt(degree[:dimen], format='%d'))
			f.write('\n')
			f.write('*Nijk\n1,\t')
			f.write(array2txt(nbctrlpts[:dimen], format='%d'))
			f.write('\n')
			f.write('*Weight\n1,\t')
			f.write(array2txt(np.ones(nbctrlpts_total)))
			
		return

	def getIGAParametrization(self):
		name = self._name
		print('\nCreating geometry: ' + name + '...')
		if name == 'quarter_annulus' or name == 'qa':
			dimen = 2
			Rin = self._extraArgs.get('Rin', 1.0)
			Rex = self._extraArgs.get('Rex', 2.0)
			if self._degree is None: self._degree = [2, 3, 1]
			part = self.__create_quarterAnnulus(Rin, Rex, *self._degree[:dimen]) 

		elif name == 'quadrilateral' or name == 'sq':
			dimen = 2
			XY = self._extraArgs.get('XY', np.array([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]))
			if self._degree is None: self._degree = [2, 2, 1]
			part = self.__create_quadrilateral(XY, *self._degree[:dimen]) 

		elif name == 'cube' or name == 'cb': 
			dimen = 3
			Lx = self._extraArgs.get('Lx', 1.0)
			Ly = self._extraArgs.get('Ly', 1.0)
			Lz = self._extraArgs.get('Lz', 1.0)
			if self._degree is None: self._degree = [2, 2, 2]
			part = self.__create_parallelepiped(Lx, Ly, Lz, *self._degree[:dimen]) 

		elif name == 'thick_ring' or name == 'tr':
			dimen = 3
			Rin = self._extraArgs.get('Rin', 1.0)
			Rex = self._extraArgs.get('Rex', 2.0)
			height = self._extraArgs.get('height', 1.0)
			if self._degree is None: self._degree = [4, 4, 4]
			part = self.__create_thickRing(Rin, Rex, height, *self._degree[:dimen]) 

		elif name == 'rotated_quarter_annulus' or name == 'rqa':
			dimen = 3
			Rin = self._extraArgs.get('Rin', 1.0)
			Rex = self._extraArgs.get('Rex', 2.0)
			exc = self._extraArgs.get('exc', 1.0) 
			if self._degree is None: self._degree = [4, 4, 4]
			part = self.__create_rotatedQuarterAnnulus(Rin, Rex, exc, *self._degree[:dimen]) 

		elif name == 'prism' or name == 'vb':
			dimen = 3
			XY     = self._extraArgs.get('xy', np.array([[0.0, -7.5], [6.0, -2.5], [6.0, 2.5], [0.0, 7.5]]))
			height = self._extraArgs.get('height', 1)
			if self._degree is None: self._degree = [2, 2, 2]
			part = self.__create_prism(XY, height, *self._degree[:dimen]) 

		else: raise Warning("Not developped in this library")

		self._dim = dimen
		info = self.__getInfo(part); info['name'] = name
		self.__writeAbaqusFile(info)
		modelIGA = IGAparametrization(filename=name)
		modelIGA.refine(nb_refinementByDirection=self._nbcuts)

		# Clean files created
		os.remove(name + '.inp')
		os.remove(name + '.NB')
		os.remove(name + '.save')

		return modelIGA

	# ----------------
	# CREATE GEOMETRY
	# ----------------
	# 2D
	def __create_quarterAnnulus(self, Rin, Rex, degree_u, degree_v):
		" Creates a quarter of a ring (or annulus) "

		# -------------------------------------
		# First part : construction of the arc
		# -------------------------------------
		nb_ctrlpts_v  = degree_v + 1 
		knot_vector_v = createUniformMaxregularKnotvector(degree_v, 1)

		# Construct points to be interpolated
		theta = np.linspace(0.0, np.pi/2.0, nb_ctrlpts_v)
		radius = 1.0
		pts_interp = [[radius, 0.0]]
		for angle in theta[1:-1]: 
			x = radius * np.cos(angle)
			y = radius * np.sin(angle)
			pts_interp.append([x, y])
		pts_interp.append([0.0, radius])

		curve = fitting.interpolate_curve(pts_interp, degree_v)
		ctrlpts_arc = np.asarray(curve.ctrlpts)

		# -------------------------------------
		# Second part : construction of line
		# -------------------------------------
		nb_ctrlpts_u = degree_u + 1 
		ctrlpts_line = np.linspace(Rin, Rex, nb_ctrlpts_u)
		knot_vector_u = createUniformMaxregularKnotvector(degree_u, 1)

		# -------------------------------------------
		# Third part : construction of annulus sector
		# -------------------------------------------
		# Set control points
		ctrlpts = []
		for x_line in ctrlpts_line: 
			for x_arc, y_arc in ctrlpts_arc:
				x = x_line * x_arc
				y = x_line * y_arc
				ctrlpts.append([x, y, 0])

		# Create surface
		obj = BSpline.Surface()
		obj.degree_u = degree_u
		obj.degree_v = degree_v
		obj.ctrlpts_size_u, obj.ctrlpts_size_v = int(nb_ctrlpts_u), int(nb_ctrlpts_v)
		obj.ctrlpts = ctrlpts
		obj.knotvector_u = knot_vector_u
		obj.knotvector_v = knot_vector_v

		return obj

	def __create_quadrilateral(self, XY, degree_u, degree_v):
		" Creates a quadrilateral given coordinates in counterclockwise direction "

		# Set reference position and real position
		x0 = [0, 1, 1, 0]; y0 = [0, 0, 1, 1]
		x1 = XY[:, 0];     y1 = XY[:, 1]

		# Transformation of control points
		# x1 = ax1 x0 + ax2 y0 + ax3 x0 y0 + ax4
		# y1 = ay1 x0 + ay2 y0 + ay3 x0 y0 + ay4
		T = [[x0[i], y0[i], x0[i]*y0[i], 1] for i in range(4)]
		ax = np.linalg.solve(T, x1)
		ay = np.linalg.solve(T, y1)
		
		# Set control points
		nb_ctrlpts_u = degree_u + 1 
		ctrlpts_u = np.linspace(0, 1, nb_ctrlpts_u)
		nb_ctrlpts_v = degree_v + 1 
		ctrlpts_v = np.linspace(0, 1, nb_ctrlpts_v)
		
		ctrlpts = []
		for xt in ctrlpts_u: 
			for yt in ctrlpts_v:   

				x = ax[0]*xt + ax[1]*yt + ax[2]*xt*yt + ax[3]
				y = ay[0]*xt + ay[1]*yt + ay[2]*xt*yt + ay[3]

				ctrlpts.append([x, y, 0])

		# Create surface
		obj = BSpline.Surface() 
		obj.degree_u = degree_u
		obj.degree_v = degree_v
		obj.ctrlpts_size_u, obj.ctrlpts_size_v = int(nb_ctrlpts_u), int(nb_ctrlpts_v)
		obj.ctrlpts = ctrlpts
		obj.knotvector_u = createUniformMaxregularKnotvector(degree_u, 1)
		obj.knotvector_v = createUniformMaxregularKnotvector(degree_v, 1)

		return obj

	# 3D
	def __create_parallelepiped(self, Lx, Ly, Lz, degree_u, degree_v, degree_w):
		" Creates a brick (or parallelepiped) "

		# Set number of control points
		nb_ctrlpts_u = degree_u + 1 
		nb_ctrlpts_v = degree_v + 1 
		nb_ctrlpts_w = degree_w + 1 

		# Get uniform control points
		ctrlpts_u = np.linspace(0.0, Lx, nb_ctrlpts_u)
		ctrlpts_v = np.linspace(0.0, Ly, nb_ctrlpts_v)
		ctrlpts_w = np.linspace(0.0, Lz, nb_ctrlpts_w)

		# Create control points of the volume
		control_points = []
		for cptw in ctrlpts_w:
			for cptu in ctrlpts_u: 
				for cptv in ctrlpts_v: 
					control_points.append([cptu, cptv, cptw])

		# Create a B-spline volume
		obj = BSpline.Volume()
		obj.degree_u, obj.degree_v, obj.degree_w = degree_u, degree_v, degree_w
		obj.ctrlpts_size_u, obj.ctrlpts_size_v, obj.ctrlpts_size_w = int(nb_ctrlpts_u), int(nb_ctrlpts_v), int(nb_ctrlpts_w)
		obj.ctrlpts = control_points
		obj.knotvector_u = createUniformMaxregularKnotvector(degree_u, 1)
		obj.knotvector_v = createUniformMaxregularKnotvector(degree_v, 1)
		obj.knotvector_w = createUniformMaxregularKnotvector(degree_w, 1)

		return obj

	def __create_thickRing(self, Rin, Rex, height, degree_u, degree_v, degree_w):
		" Creates a thick ring (quarter of annulus extruded) "

		# -------------------------------------
		# First part : construction of the arc
		# -------------------------------------
		nb_ctrlpts_u = degree_v + 1 
		knot_vector_u = createUniformMaxregularKnotvector(degree_v, 1)

		# Construct points to be interpolated
		theta = np.linspace(0.0, np.pi/2.0, nb_ctrlpts_u)
		radius = 1.0
		pts_interp = [[radius, 0.0]]
		for angle in theta[1:-1]: 
			x = radius*np.cos(angle)
			y = radius*np.sin(angle)
			pts_interp.append([x, y])
		pts_interp.append([0.0, radius])

		curve = fitting.interpolate_curve(pts_interp, degree_v)
		ctrlpts_arc = np.asarray(curve.ctrlpts)

		# -------------------------------------
		# Second part : construction of line
		# -------------------------------------
		# Define degree in u direction
		nb_ctrlpts_u = degree_u + 1 
		ctrlpts_line = np.linspace(Rin, Rex, nb_ctrlpts_u)
		knot_vector_u = createUniformMaxregularKnotvector(degree_u, 1)

		# Define degree in w direction
		nb_ctrlpts_w = degree_w + 1 
		ctrlpts_height = np.linspace(0, height, nb_ctrlpts_w)
		knot_vector_w = createUniformMaxregularKnotvector(degree_w, 1)

		# -------------------------------------------
		# Third part : construction of annulus sector
		# -------------------------------------------
		# Get control points
		ctrlpts = []
		for z in ctrlpts_height:
			for x_line in ctrlpts_line: 
				for x_arc, y_arc in ctrlpts_arc:
					x = x_line*x_arc
					y = x_line*y_arc
					ctrlpts.append([x, y, z])

		# Create volume
		obj = BSpline.Volume()
		obj.degree_u, obj.degree_v, obj.degree_w = degree_u, degree_v, degree_w
		obj.ctrlpts_size_u, obj.ctrlpts_size_v, obj.ctrlpts_size_w = int(nb_ctrlpts_u), int(nb_ctrlpts_u), int(nb_ctrlpts_w)
		obj.ctrlpts = ctrlpts
		obj.knotvector_u = knot_vector_u
		obj.knotvector_v = knot_vector_u
		obj.knotvector_w = knot_vector_w

		return obj

	def __create_rotatedQuarterAnnulus(self, Rin, Rex, exc, degree_u, degree_v, degree_w):
		" Creates a quarter of a ring rotated (or revolted) "

		# -------------------------------------
		# First part : construction of the arc 1
		# -------------------------------------
		nb_ctrlpts_v = degree_v + 1 
		knot_vector_v = createUniformMaxregularKnotvector(degree_v, 1)

		# Construct points to be interpolated
		theta = np.linspace(0.0, np.pi/2.0, nb_ctrlpts_v)
		radius = 1.0
		pts_interp = [[radius, 0.0]]
		for angle in theta[1:-1]: 
			x = radius*np.cos(angle)
			y = radius*np.sin(angle)
			pts_interp.append([x, y])
		pts_interp.append([0.0, radius])

		curve = fitting.interpolate_curve(pts_interp, degree_v)
		ctrlpts_arc_1 = np.asarray(curve.ctrlpts)

		# -------------------------------------
		# Second part : construction of line
		# -------------------------------------
		nb_ctrlpts_u = degree_u + 1 
		ctrlpts_line = np.linspace(Rin, Rex, nb_ctrlpts_u)
		knot_vector_u = createUniformMaxregularKnotvector(degree_u, 1)

		# -------------------------------------
		# Third part : construction of the arc 2
		# -------------------------------------
		nb_ctrlpts_w = degree_w + 1 
		knot_vector_w = createUniformMaxregularKnotvector(degree_w, 1)

		# Construct points to be interpolated
		theta = np.linspace(0.0, np.pi/2.0, nb_ctrlpts_w)
		radius = 1.0
		pts_interp = [[radius, 0.0]]
		for angle in theta[1:-1]: 
			x = radius * np.cos(angle)
			y = radius * np.sin(angle)
			pts_interp.append([x, y])
		pts_interp.append([0.0, radius])

		curve = fitting.interpolate_curve(pts_interp, degree_w)
		ctrlpts_arc_2 = np.asarray(curve.ctrlpts)

		# -------------------------------------------
		# Fourth part : construction of annulus sector
		# -------------------------------------------
		# Get control points
		ctrlpts = []
		for y_arc_2, z_arc_2 in ctrlpts_arc_2:
			for x_line in ctrlpts_line: 
				for x_arc_1, y_arc_1 in ctrlpts_arc_1:
					# ctrlpts_x[i] = radius_1
					x = x_line*x_arc_1
					y = (x_line*y_arc_1 + exc)*y_arc_2 
					z = (x_line*y_arc_1 + exc)*z_arc_2
					ctrlpts.append([x, y, z])

		# Create volume
		obj = BSpline.Volume()
		obj.degree_u, obj.degree_v, obj.degree_w = degree_u, degree_v, degree_w
		obj.ctrlpts_size_u, obj.ctrlpts_size_v, obj.ctrlpts_size_w = int(nb_ctrlpts_u), int(nb_ctrlpts_v), int(nb_ctrlpts_w)
		obj.ctrlpts = ctrlpts
		obj.knotvector_u = knot_vector_u
		obj.knotvector_v = knot_vector_v
		obj.knotvector_w = knot_vector_w

		return obj

	def __create_prism(self, XY, height, degree_u, degree_v, degree_w):
		""" Creates a prism using a quadrilateral as a base.
		The quadrilateral coordinates are given in counterclockwise direction """

		# Set reference position and real position
		x0 = [0, 1, 1, 0]; y0 = [0, 0, 1, 1]
		x1 = XY[:, 0];     y1 = XY[:, 1]

		# Transformation of control points
		# x1 = ax1 x0 + ax2 y0 + ax3 x0 y0 + ax4
		# y1 = ay1 x0 + ay2 y0 + ay3 x0 y0 + ay4
		T = [[x0[i], y0[i], x0[i]*y0[i], 1] for i in range(4)]
		ax = np.linalg.solve(T, x1)
		ay = np.linalg.solve(T, y1)

		# Set control points
		nb_ctrlpts_u = degree_u + 1 
		ctrlpts_u = np.linspace(0, 1, nb_ctrlpts_u)
		nb_ctrlpts_v = degree_v + 1 
		ctrlpts_v = np.linspace(0, 1, nb_ctrlpts_v)
		nb_ctrlpts_w = degree_w + 1 
		ctrlpts_w = np.linspace(0, 1, nb_ctrlpts_w)
		
		ctrlpts = []
		for zt in ctrlpts_w:
			for xt in ctrlpts_u: 
				for yt in ctrlpts_v:      
					z = zt*height
					x = ax[0]*xt + ax[1]*yt + ax[2]*xt*yt + ax[3]
					y = ay[0]*xt + ay[1]*yt + ay[2]*xt*yt + ay[3]
					ctrlpts.append([x, y, z])

		# Create volume
		obj = BSpline.Volume()
		obj.degree_u, obj.degree_v, obj.degree_w = degree_u, degree_v, degree_w
		obj.ctrlpts_size_u, obj.ctrlpts_size_v, obj.ctrlpts_size_w = int(nb_ctrlpts_u), int(nb_ctrlpts_v), int(nb_ctrlpts_w)
		obj.ctrlpts = ctrlpts
		obj.knotvector_u = createUniformMaxregularKnotvector(degree_u, 1)
		obj.knotvector_v = createUniformMaxregularKnotvector(degree_v, 1)
		obj.knotvector_w = createUniformMaxregularKnotvector(degree_w, 1)

		return obj

