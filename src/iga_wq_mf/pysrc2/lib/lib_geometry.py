"""
.. This module contains functions to construct specific geometries using only B-splines
.. Joaquin Cornejo
"""

from .__init__ import *
from .lib_base import createKnotVector

class geomdlModel(): 

	def __init__(self, name=None, **geometry): 

		if name is None: raise Warning('Insert the name of the part')
		self._name = name
		self._sample_size = 101
		self._geometry = None

		print('\nCreating geometry: ' + name + '...')
		start = time.process_time()
		if name == 'quarter_annulus' or name == 'QA':
			self._dim = 2
			Rin = geometry.get('Rin', 1.0)
			Rex = geometry.get('Rex', 2.0)
			degree_u, degree_v, _ = geometry.get('degree', [2, 3, 2])
			self._geometry = self.create_quarterAnnulus(Rin, Rex, degree_u, degree_v) 

		elif name == 'quadrilateral' or name == 'SQ':
			self._dim = 2
			XY = geometry.get('XY', np.array([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]))
			degree_u, degree_v, _ = geometry.get('degree', [2, 2, 2])
			self._geometry = self.create_quadrilateral(XY, degree_u, degree_v) 

		elif name == 'parallelepiped' or name == 'CB': 
			self._dim = 3
			Lx = geometry.get('Lx', 1.0)
			Ly = geometry.get('Ly', 1.0)
			Lz = geometry.get('Lz', 1.0)
			degree_u, degree_v, degree_w = geometry.get('degree', [2, 2, 2])
			self._geometry = self.create_parallelepiped(Lx, Ly, Lz, degree_u, degree_v, degree_w) 

		elif name == 'thick_ring' or name == 'TR':
			self._dim = 3
			Rin = geometry.get('Rin', 1.0)
			Rex = geometry.get('Rex', 2.0)
			height = geometry.get('Height', 1.0)
			degree_u, degree_v, degree_w = geometry.get('degree', [4, 4, 4])
			self._geometry = self.create_thickRing(Rin, Rex, height, degree_u, degree_v, degree_w) 

		elif name == 'rotated_quarter_annulus' or name == 'RQA':
			self._dim = 3
			Rin = geometry.get('Rin', 1.0)
			Rex = geometry.get('Rex', 2.0)
			exc = geometry.get('exc', 1.0) 
			degree_u, degree_v, degree_w = geometry.get('degree', [4, 4, 4])
			self._geometry = self.create_rotatedQuarterAnnulus(Rin, Rex, exc, degree_u, degree_v, degree_w) 

		elif name == 'prism' or name == 'VB':
			self._dim = 3
			XY = geometry.get('XY', np.array([[0.0, -7.5], [6.0, -2.5], [6.0, 2.5], [0.0, 7.5]]))
			height = geometry.get('Height', 1)
			degree_u, degree_v, degree_w = geometry.get('degree', [2, 2, 2])
			self._geometry = self.create_prism(XY, height, degree_u, degree_v, degree_w) 

		else: raise Warning("Not developped in this library")
		
		stop = time.process_time()
		print('\tBasic geometry created in: %.3e s' %(stop-start))

		self.updateGeometry()

		return

	def updateGeometry(self): 
		" Updates and saves important properties of the geometry created "

		start = time.process_time()

		obj = self._geometry
		if obj is None: raise Warning('Geometry unknown')

		# Set degree
		self._degree = np.ones(3, dtype= int)
		self._degree[:self._dim] = np.array(obj._degree)
		if any(p == 1 for p in self._degree[:self._dim]): 
			raise Warning('Model must have at least degree 2')

		# Set knot vector
		self._knotvector = [np.array(obj._knot_vector[dim]) for dim in range(self._dim)]

		# Set size knot-vector in each dimension
		self._size_kv = np.ones(3, dtype= int)  
		for dim in range(self._dim):
			self._size_kv[dim] = np.size(self._knotvector[dim])

		# Set number of elements in each dimension
		self._nbel = np.ones(3, dtype= int)  
		for dim in range(self._dim):
			self._nbel[dim] = len(np.unique(self._knotvector[dim])) - 1
		self._nb_el_total = np.product(self._nbel)

		# Set number of quadrature points in each dimension
		self._nb_ctrlpts = np.ones(3, dtype= int)  
		for dim in range(self._dim):
			self._nb_ctrlpts[dim] = self._size_kv[dim] - self._degree[dim] - 1
		self._nb_ctrlpts_total = np.product(self._nb_ctrlpts)

		if self._dim == 2: 
			c = 0
			ctrlpts_old = obj._control_points
			ctrlpts_new = np.zeros((3, self._nb_ctrlpts_total))
			for j in range(self._nb_ctrlpts[1]):
				for i in range(self._nb_ctrlpts[0]):
					pos = j + i*self._nb_ctrlpts[1]
					ctrlpts_new[:, c] = ctrlpts_old[pos]
					c += 1 

		elif self._dim == 3: 
			c =  0
			ctrlpts_old = obj._control_points
			ctrlpts_new = np.zeros((3, self._nb_ctrlpts_total))
			for k in range(self._nb_ctrlpts[2]):
				for j in range(self._nb_ctrlpts[1]):
					for i in range(self._nb_ctrlpts[0]):
						pos = j + i*self._nb_ctrlpts[1] + k*self._nb_ctrlpts[1]*self._nb_ctrlpts[0]
						ctrlpts_new[:, c] = ctrlpts_old[pos]
						c += 1 
		self._ctrlpts = ctrlpts_new

		stop = time.process_time()
		print('\tGeometry properties updated in: %.3e s\n' %(stop-start))

		return

	def writeAbaqusFile(self, filename):
		" Writes an inp and NB file. By the moment, it only works with one patch"

		def array2txt(array: np.array, format= '%.2f'):
			return ','.join([format %(i) for i in array])

		# .inp file
		inpfile = filename + '.inp'
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
			f.write('*Part, name=%s\n' %self._name)
			f.write('*USER ELEMENT, NODES=%d, TYPE=U1, COORDINATES=%d, INTEGRATION=%d\n' 
					%(self._nb_ctrlpts_total, self._dim, self._nb_ctrlpts_total))
			f.write(array2txt(np.arange(self._dim)+1, format='%d'))
			f.write('\n*Node,nset=AllNode\n')
			for i in range(self._nb_ctrlpts_total):
				CP = self._ctrlpts[:, i]
				f.write('%d, %.15f, %.15f, %.15f\n' %(i+1, CP[0], CP[1], CP[2]))
			f.write('*Element,type=U1,elset=AllEls\n1,\t')
			f.write(array2txt(np.arange(self._nb_ctrlpts_total, 0, -1), format='%d'))
			f.write('\n')
			f.write('*ELSET,ELSET=EltPatch1,generate\n1,1,1\n')
			f.write('*UEL PROPERTY, ELSET=EltPatch1, MATERIAL=Mat\n1\n')
			f.write('*End Part\n')
			f.write('**ASSEMBLY\n*Assembly, name=Assembly\n')
			f.write('*Instance, name=I1, part=%s\n' %self._name)
			f.write('*End Instance\n*End Assembly\n')
			f.write('**MATERIAL\n*MATERIAL,NAME=Mat\n*Elastic\n')
			f.write('%f, %f\n' %(3e3, 0.3)) 
			f.write('*STEP,extrapolation=NO,NLGEOM=NO\n*Static\n')
			f.write('** OUTPUT REQUESTS\n*node file,frequency=1\nU,RF,CF\n*el file,frequency=1\nSDV\n*End Step')

		# .NB file
		NBfile = filename + '.NB'
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
			f.write('*Dimension\n%d\n' %(self._dim))
			f.write('*Number of CP by element\n%d\n' %(self._nb_ctrlpts_total))
			f.write('*Number of patch\n%d\n' %(1))
			f.write('*Total number of element \n%d\n' %(1))
			f.write('*Number of element by patch\n%d\n' %(1))
			f.write('*Patch(1)\n')
			for i in range(self._dim):
				kv = self._knotvector[i]
				f.write('%d\n' %(len(kv)))
				f.write(array2txt(kv))
				f.write('\n')
			f.write('*Jpqr\n')
			f.write(array2txt(self._degree[:self._dim], format='%d'))
			f.write('\n')
			f.write('*Nijk\n1,\t')
			f.write(array2txt(self._nb_ctrlpts[:self._dim], format='%d'))
			f.write('\n')
			f.write('*Weight\n1,\t')
			f.write(array2txt(np.ones(self._nb_ctrlpts_total)))
			
		return

	def knotRefinement(self, nb_refinementByDirection=np.array([0,0,0])):
		""" Refine geometry following each dimension. 
			It is slow because it uses python methods. 
			This functions is deprecated. Instead use YETI functions. 
		"""

		start = time.process_time()

		geometry = deepcopy(self._geometry)
		cuts = nb_refinementByDirection
		nbel = [2**cuts[dim]*self._nbel[dim] for dim in range(self._dim)]
		knotvector_insert = [np.linspace(0.0, 1.0, i+1)[1:-1] for i in nbel]

		for i in range(self._dim):
			multiplicity = np.zeros(self._dim, dtype= int)
			multiplicity[i] = 1
			for knot in knotvector_insert[i]: 
				knot_insert = np.zeros(self._dim)
				knot_insert[i] = knot
				operations.insert_knot(geometry, knot_insert, multiplicity.tolist())
		self._geometry = geometry

		stop = time.process_time()
		print('Knot refinement in: %.3e s' %(stop-start))

		self.updateGeometry()

		return

	def fastKnotRefinement(self, nb_refinementByDirection=np.array([0,0,0])):
		""" Refine geometry following each dimension. 
			It has a better performance than knot-refinement function
		"""

		self.writeAbaqusFile(filename=self._name)
		modelIGA = IGAparametrization(filename=self._name)
		modelIGA.refine(nb_refinementByDirection=nb_refinementByDirection)

		# Clean files created
		os.remove(self._name + '.inp')
		os.remove(self._name + '.NB')
		os.remove(self._name + '.save')

		return modelIGA

	# ----------------
	# CREATE GEOMETRY
	# ----------------
	# 2D
	def create_quarterAnnulus(self, Rin, Rex, degree_u, degree_v):
		" Creates a quarter of a ring (or annulus) "

		# -------------------------------------
		# First part : construction of the arc
		# -------------------------------------
		nb_ctrlpts_v  = degree_v + 1 
		knot_vector_v = createKnotVector(degree_v, 1)

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
		knot_vector_u = createKnotVector(degree_u, 1)

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
		obj.ctrlpts_size_u, obj.ctrlpts_size_v = nb_ctrlpts_u, nb_ctrlpts_v
		obj.ctrlpts = ctrlpts
		obj.knotvector_u = knot_vector_u
		obj.knotvector_v = knot_vector_v
		obj.sample_size  = self._sample_size

		return obj

	def create_quadrilateral(self, XY, degree_u, degree_v):
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
		obj.ctrlpts_size_u, obj.ctrlpts_size_v = nb_ctrlpts_u, nb_ctrlpts_v
		obj.ctrlpts = ctrlpts
		obj.knotvector_u = createKnotVector(degree_u, 1)
		obj.knotvector_v = createKnotVector(degree_v, 1)
		obj.sample_size = self._sample_size

		return obj

	# 3D
	def create_parallelepiped(self, Lx, Ly, Lz, degree_u, degree_v, degree_w):
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
			for cptv in ctrlpts_v: 
				for cptu in ctrlpts_u: 
					control_points.append([cptu, cptv, cptw])

		# Create a B-spline volume
		obj = BSpline.Volume()
		obj.degree_u, obj.degree_v, obj.degree_w = degree_u, degree_v, degree_w
		obj.ctrlpts_size_u, obj.ctrlpts_size_v, obj.ctrlpts_size_w = int(nb_ctrlpts_u), int(nb_ctrlpts_v), int(nb_ctrlpts_w)
		obj.ctrlpts = control_points
		obj.knotvector_u = createKnotVector(degree_u, 1)
		obj.knotvector_v = createKnotVector(degree_v, 1)
		obj.knotvector_w = createKnotVector(degree_w, 1)
		obj.sample_size = self._sample_size

		return obj

	def create_thickRing(self, Rin, Rex, height, degree_u, degree_v, degree_w):
		" Creates a thick ring (quarter of annulus extruded) "

		# -------------------------------------
		# First part : construction of the arc
		# -------------------------------------
		nb_ctrlpts_u = degree_v + 1 
		knot_vector_u = createKnotVector(degree_v, 1)

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
		knot_vector_u = createKnotVector(degree_u, 1)

		# Define degree in w direction
		nb_ctrlpts_w = degree_w + 1 
		ctrlpts_height = np.linspace(0, height, nb_ctrlpts_w)
		knot_vector_w = createKnotVector(degree_w, 1)

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
		obj.sample_size = self._sample_size

		return obj

	def create_rotatedQuarterAnnulus(self, Rin, Rex, exc, degree_u, degree_v, degree_w):
		" Creates a quarter of a ring rotated (or revolted) "

		# -------------------------------------
		# First part : construction of the arc 1
		# -------------------------------------
		nb_ctrlpts_v = degree_v + 1 
		knot_vector_v = createKnotVector(degree_v, 1)

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
		knot_vector_u = createKnotVector(degree_u, 1)

		# -------------------------------------
		# Third part : construction of the arc 2
		# -------------------------------------
		nb_ctrlpts_w = degree_w + 1 
		knot_vector_w = createKnotVector(degree_w, 1)

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
		obj.sample_size = self._sample_size

		return obj

	def create_prism(self, XY, height, degree_u, degree_v, degree_w):
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
		obj.knotvector_u = createKnotVector(degree_u, 1)
		obj.knotvector_v = createKnotVector(degree_v, 1)
		obj.knotvector_w = createKnotVector(degree_w, 1)
		obj.sample_size = self._sample_size

		return obj
