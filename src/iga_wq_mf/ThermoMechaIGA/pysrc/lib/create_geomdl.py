"""
.. This module helps to construct some geometries using only B-splines
.. Joaquin Cornejo
"""

# Python libraries
from geomdl import (
    fitting, 
    BSpline, 
    operations,
)
import matplotlib.pyplot as plt
import math, numpy as np
import copy, time

# My libraries
from .base_functions import create_knotvector

class geomdlModel(): 

    def __init__(self, filename= None, **geometry: None): 

        if filename is None:
            raise Warning('Insert filename')
        else : 
            # Set name of object
            self._name = filename

        # Set number of samples
        self._sample_size = 61

        print('\nCreating geometry: ' + filename + '...')
        start = time.time()
        if filename == 'quarter_annulus':
            # Set number of dimension 
            self._dim = [2]

            # Create quarter annulus
            Rin = geometry.get('Rin', 1.0)
            Rout = geometry.get('Rout', 2.0)
            degree_xi, degree_nu = geometry.get('degree', [2, 3])
            self._geometry = self.create_quarter_annulus(Rin, Rout, degree_xi, degree_nu) 

        elif filename == 'quadrilateral':
            # Set number of dimension
            self._dim = [2]
            
            # Create parallelepiped
            XY = geometry.get('XY', np.array([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]))
            degree_xi, degree_nu = geometry.get('degree', [2, 2])
            self._geometry = self.create_quadrilateral(XY, degree_xi, degree_nu) 

        elif filename == 'parallelepiped': 
            # Set number of dimension
            self._dim = [3]
            
            # Create parallelepiped
            Lx = geometry.get('Lx', 1.0)
            Ly = geometry.get('Ly', 1.0)
            Lz = geometry.get('Lz', 1.0)
            degree_xi, degree_nu, degree_eta = geometry.get('degree', [2, 2, 2])
            self._geometry = self.create_parallelepiped(Lx, Ly, Lz, degree_xi, degree_nu, degree_eta) 

        elif filename == 'thick_ring':
            # Set number of dimension 
            self._dim = [3]

            # Create quarter annulus
            Rin = geometry.get('Rin', 1.0)
            Rout = geometry.get('Rout', 2.0)
            Height = geometry.get('Height', 1.0)
            degree_xi, degree_nu, degree_eta = geometry.get('degree', [4, 4, 4])
            self._geometry = self.create_thick_ring(Rin, Rout, Height, degree_xi, degree_nu, degree_eta) 

        elif filename == 'rotated_quarter_annulus':
            # Set number of dimension 
            self._dim = [3]

            # Create quarter annulus
            Rin = geometry.get('Rin', 1.0)
            Rout = geometry.get('Rout', 2.0)
            exc = geometry.get('Rout', 1.0) # excentricity
            degree_xi, degree_nu, degree_eta = geometry.get('degree', [4, 4, 4])
            self._geometry = self.create_rotated_quarter_annulus(Rin, Rout, exc, degree_xi, degree_nu, degree_eta) 

        elif filename == 'prism':
            # Set number of dimension
            self._dim = [3]
            
            # Create parallelepiped
            XY = geometry.get('XY', np.array([[0.0, -7.5], [6.0, -2.5], [6.0, 2.5], [0.0, 7.5]]))
            Height = geometry.get('Height', 2)
            degree_xi, degree_nu, degree_eta = geometry.get('degree', [2, 2, 2])
            self._geometry = self.create_prisme(XY, Height, degree_xi, degree_nu, degree_eta) 

        else: 
            raise Warning("Not a shape in this library")
        
        stop = time.time()
        print('\tBasic geometry created in: %.3e s' %(stop-start))

        # Update data
        self.update_geometry()

        return

    def knot_refinement(self, nb_refinementByDirection= np.array([0,0,0])):
        
        start = time.time()
        # Copy geometry
        geometry = copy.deepcopy(self._geometry)

        # Set new number of elements
        cuts = nb_refinementByDirection
        nbel = [2**cuts[dim]*self._nb_el[dim][0] for dim in range(self._dim[0])]

        # Create knots to be inserted
        knotvector_insert = [np.linspace(0.0, 1.0, i+1)[1:-1] for i in nbel]

        # Insert knots
        for dim in range(self._dim[0]):
            multiplicity = np.zeros(self._dim[0], dtype= int)
            multiplicity[dim] = 1
            for knot in knotvector_insert[dim]: 
                knot_insert = np.zeros(self._dim[0])
                knot_insert[dim] = knot
                operations.insert_knot(geometry, knot_insert, multiplicity.tolist())
        self._geometry = geometry
        stop = time.time()
        print('Knot refinement in: %.3e s' %(stop-start))

        # Update values
        self.update_geometry()

        return

    def update_geometry(self): 
        start = time.time()
        # Set geometry object
        obj = self._geometry
        if obj is None: raise Warning('Geometry unknown')

        # Set degree
        self._degree = np.zeros((3, 1), dtype= int)
        self._degree[:self._dim[0], 0] = np.asarray(obj._degree)
        if any(p == 1 for p in self._degree[:self._dim[0]]): 
           raise Warning('Model must have at least degree p = 2')

        # Set knot vector
        self._knotvector = [[np.asarray(obj._knot_vector[dim]) for dim in range(self._dim[0])]]

        # Set size knot-vector in each dimension
        self._size_kv = np.zeros((3, 1), dtype= int)  
        for dim in range(self._dim[0]):
            self._size_kv[dim] = np.size(self._knotvector[0][dim])

        # Set number of elements in each dimension
        self._nb_el = np.zeros((3, 1), dtype= int)  
        self._nb_el[:self._dim[0], 0] = self._size_kv[:self._dim[0], 0]\
                                    - (2*self._degree[:self._dim[0], 0] + 1)

        # Set total number of elements
        self._nb_el_total = np.product(self._nb_el[:self._dim[0]][:])

        # Set number of quadrature points/functions in each dimension
        self._nb_ctrlpts = self._degree + self._nb_el
        self._nb_ctrlpts_total = np.product(self._nb_ctrlpts[:self._dim[0]][:])

        if self._dim[0] == 2: 
            ctrlpts_old = np.asarray(obj._control_points)
            ctrlpts_new = []
            for j in range(self._nb_ctrlpts[1][0]):
                for i in range(self._nb_ctrlpts[0][0]):
                    pos = j + i*self._nb_ctrlpts[1][0]
                    ctrlpts_temp = ctrlpts_old[pos]
                    ctrlpts_new.append(ctrlpts_temp)

        elif self._dim[0] == 3: 
            ctrlpts_old = obj._control_points
            ctrlpts_new = []
            for k in range(self._nb_ctrlpts[2][0]):
                for j in range(self._nb_ctrlpts[1][0]):
                    for i in range(self._nb_ctrlpts[0][0]):
                        pos = j + i*self._nb_ctrlpts[1][0] + k*self._nb_ctrlpts[1][0]*self._nb_ctrlpts[0][0]
                        ctrlpts_temp = ctrlpts_old[pos]
                        ctrlpts_new.append(ctrlpts_temp)

        # Update control points
        self._ctrlpts = ctrlpts_new
        stop = time.time()
        print('\tGeometry properties updated in: %.3e s\n' %(stop-start))

        return

    # -----------------------------------
    # CREATE GEOMETRY
    # -----------------------------------
    # 2D
    def create_quarter_annulus(self, Rin, Rout, degree_xi, degree_nu):
        " Creates a quarter of a ring "

        # -------------------------------------
        # First part : construction of the arc
        # -------------------------------------
        # Set knot-vector characteristics
        nb_ctrlpts_nu = degree_nu + 1 
        knot_vector_nu = create_knotvector(degree_nu, 1)

        # Construct points to be interpolated
        theta = np.linspace(0.0, math.pi/2.0, nb_ctrlpts_nu)
        radius = 1.0
        pts_interp = [[radius, 0.0]]
        for angle in theta[1:-1]: 
            x = radius * np.cos(angle)
            y = radius * np.sin(angle)
            pts_interp.append([x, y])

        pts_interp.append([0.0, radius])

        # Do interpolation
        curve = fitting.interpolate_curve(pts_interp, degree_nu)

        # Get control points and knot-vector
        ctrlpts_arc = np.asarray(curve.ctrlpts)

        # -------------------------------------
        # Second part : construction of line
        # -------------------------------------
        # Define degree in this direction
        nb_ctrlpts_xi = degree_xi + 1 
        ctrlpts_x = np.linspace(Rin, Rout, nb_ctrlpts_xi)
        knot_vector_xi = create_knotvector(degree_xi, 1)

        # -------------------------------------------
        # Third part : construction of annulus sector
        # -------------------------------------------
        # Get control points
        ctrlpts = []

        # Over x direction 
        for i in range(nb_ctrlpts_xi): 
            for x_arc, y_arc in ctrlpts_arc:
                x = ctrlpts_x[i] * x_arc
                y = ctrlpts_x[i] * y_arc
                ctrlpts.append([x, y, 0])

        # Create surface
        srf = BSpline.Surface()

        # Set degree
        srf.degree_u = degree_xi
        srf.degree_v = degree_nu

        # Set number of control points 
        srf.ctrlpts_size_u, srf.ctrlpts_size_v = nb_ctrlpts_xi, nb_ctrlpts_nu

        # Set control points
        srf.ctrlpts = ctrlpts

        # Set knot-vector
        srf.knotvector_u = knot_vector_xi
        srf.knotvector_v = knot_vector_nu

        srf.sample_size = self._sample_size

        return srf

    def create_quadrilateral(self, XY, degree_xi, degree_nu):
        " Creates a quadrilateral given coordinates in counterclockwise direction "

        # Position of reference
        x0 = [0, 1, 1, 0]
        y0 = [0, 0, 1, 1]

        # Position of quadrilateral corners
        x1 = XY[:, 0]
        y1 = XY[:, 1]

        # Transformation of control points
        # x1 = ax1 x0 + ax2 y0 + ax3 x0 y0 + ax4
        # y1 = ay1 x0 + ay2 y0 + ay3 x0 y0 + ay4
        T = []
        for i in range(4): 
            T.append([x0[i], y0[i], x0[i] * y0[i], 1])

        ax = np.linalg.solve(T, x1)
        ay = np.linalg.solve(T, y1)

        # Control points
        nb_ctrlpts_xi = degree_xi + 1 
        nb_ctrlpts_nu = degree_nu + 1 
        
        ctrlpts_xi = np.linspace(0, 1, nb_ctrlpts_xi)
        ctrlpts_nu = np.linspace(0, 1, nb_ctrlpts_nu)
        
        ctrlpts = []
        for i in range(nb_ctrlpts_xi): 
            for j in range(nb_ctrlpts_nu):    
                xt = ctrlpts_xi[i]
                yt = ctrlpts_nu[j]

                x = ax[0]*xt + ax[1]*yt + ax[2]*xt*yt + ax[3]
                y = ay[0]*xt + ay[1]*yt + ay[2]*xt*yt + ay[3]

                ctrlpts.append([x, y, 0])

        # Get knot-vector
        knot_vector_xi = create_knotvector(degree_xi, 1)
        knot_vector_nu = create_knotvector(degree_nu, 1)

        # Creation of Quadralateral surface
        srf = BSpline.Surface() 

        # Set degree
        srf.degree_u = degree_xi
        srf.degree_v = degree_nu
        
        # Set number of control points 
        srf.ctrlpts_size_u, srf.ctrlpts_size_v = nb_ctrlpts_xi, nb_ctrlpts_nu

        # Set control points
        srf.ctrlpts = ctrlpts

        # Set knot-vector
        srf.knotvector_u = knot_vector_xi
        srf.knotvector_v = knot_vector_nu

        srf.sample_size = self._sample_size

        return srf

    # 3D
    def create_parallelepiped(self, Lx, Ly, Lz, 
                                degree_xi, degree_nu, degree_eta):
        " Creates a brick "

        # Get control points
        nb_ctrlpts_xi = degree_xi + 1 
        nb_ctrlpts_nu = degree_nu + 1 
        nb_ctrlpts_eta = degree_eta + 1 

        # Get control points
        ctrlpts_u = np.linspace(0.0, Lx, nb_ctrlpts_xi)
        ctrlpts_v = np.linspace(0.0, Ly, nb_ctrlpts_nu)
        ctrlpts_w = np.linspace(0.0, Lz, nb_ctrlpts_eta)

        # Get knot-vector
        knot_vector_xi = create_knotvector(degree_xi, 1)
        knot_vector_nu = create_knotvector(degree_nu, 1)
        knot_vector_eta = create_knotvector(degree_eta, 1)

        # Create control points of the volume
        control_points = []
        for k in range(nb_ctrlpts_eta):
            for i in range(nb_ctrlpts_xi): 
                for j in range(nb_ctrlpts_nu): 
                    control_points.append([ctrlpts_u[i], ctrlpts_v[j], ctrlpts_w[k]])

        # Create a B-spline surface
        vol = BSpline.Volume()

        # Set degree
        vol.degree_u = degree_xi
        vol.degree_v = degree_nu
        vol.degree_w = degree_eta

        # Set number of control points 
        vol.ctrlpts_size_u = nb_ctrlpts_xi
        vol.ctrlpts_size_v = nb_ctrlpts_nu
        vol.ctrlpts_size_w = nb_ctrlpts_eta

        # Set control points
        vol.ctrlpts = control_points

        # Set knot-vector
        vol.knotvector_u = knot_vector_xi
        vol.knotvector_v = knot_vector_nu
        vol.knotvector_w = knot_vector_eta

        vol.sample_size = self._sample_size

        return vol

    def create_thick_ring(self, Rin, Rout, Height, degree_xi, degree_nu, degree_eta):
        # -------------------------------------
        # First part : construction of the arc
        # -------------------------------------
        # Set knot-vector characteristics
        nb_ctrlpts_nu = degree_nu + 1 
        knot_vector_nu = create_knotvector(degree_nu, 1)

        # Construct points to be interpolated
        theta = np.linspace(0.0, math.pi/2.0, nb_ctrlpts_nu)
        radius = 1.0
        pts_interp = [[radius, 0.0]]
        for angle in theta[1:-1]: 
            x = radius * np.cos(angle)
            y = radius * np.sin(angle)
            pts_interp.append([x, y])

        pts_interp.append([0.0, radius])

        # Do interpolation
        curve = fitting.interpolate_curve(pts_interp, degree_nu)

        # Get control points and knot-vector
        ctrlpts_arc = np.asarray(curve.ctrlpts)

        # -------------------------------------
        # Second part : construction of line
        # -------------------------------------
        # Define degree in xi direction
        nb_ctrlpts_xi = degree_xi + 1 
        ctrlpts_x = np.linspace(Rin, Rout, nb_ctrlpts_xi)
        knot_vector_xi = create_knotvector(degree_xi, 1)

        # Define degree in eta direction
        nb_ctrlpts_eta = degree_eta + 1 
        ctrlpts_z = np.linspace(0, Height, nb_ctrlpts_eta)
        knot_vector_eta = create_knotvector(degree_eta, 1)

        # -------------------------------------------
        # Third part : construction of annulus sector
        # -------------------------------------------
        # Get control points
        ctrlpts = []

        # Over x and z direction 
        for k in range(nb_ctrlpts_eta):
            for i in range(nb_ctrlpts_xi): 
                for x_arc, y_arc in ctrlpts_arc:
                    x = ctrlpts_x[i] * x_arc
                    y = ctrlpts_x[i] * y_arc
                    z = ctrlpts_z[k]
                    ctrlpts.append([x, y, z])

        # Create Volume
        vol = BSpline.Volume()

        # Set degree
        vol.degree_u = degree_xi
        vol.degree_v = degree_nu
        vol.degree_w = degree_eta

        # Set number of control points 
        vol.ctrlpts_size_u, vol.ctrlpts_size_v, vol.ctrlpts_size_w = \
            nb_ctrlpts_xi, nb_ctrlpts_nu, nb_ctrlpts_eta

        # Set control points
        vol.ctrlpts = ctrlpts

        # Set knot-vector
        vol.knotvector_u = knot_vector_xi
        vol.knotvector_v = knot_vector_nu
        vol.knotvector_w = knot_vector_eta

        vol.sample_size = self._sample_size

        return vol

    def create_rotated_quarter_annulus(self, Rin, Rout, exc, degree_xi, degree_nu, degree_eta):
        " Creates a quarter of a ring "

        # -------------------------------------
        # First part : construction of the arc 1
        # -------------------------------------
        # Set knot-vector characteristics
        nb_ctrlpts_nu = degree_nu + 1 
        knot_vector_nu = create_knotvector(degree_nu, 1)

        # Construct points to be interpolated
        theta = np.linspace(0.0, math.pi/2.0, nb_ctrlpts_nu)
        radius = 1.0
        pts_interp = [[radius, 0.0]]
        for angle in theta[1:-1]: 
            x = radius * np.cos(angle)
            y = radius * np.sin(angle)
            pts_interp.append([x, y])

        pts_interp.append([0.0, radius])

        # Do interpolation
        curve = fitting.interpolate_curve(pts_interp, degree_nu)

        # Get control points and knot-vector
        ctrlpts_arc_1 = np.asarray(curve.ctrlpts)

        # -------------------------------------
        # Second part : construction of line
        # -------------------------------------
        # Define degree in xi direction
        nb_ctrlpts_xi = degree_xi + 1 
        ctrlpts_x = np.linspace(Rin, Rout, nb_ctrlpts_xi)
        knot_vector_xi = create_knotvector(degree_xi, 1)

        # -------------------------------------
        # Third part : construction of the arc 2
        # -------------------------------------
        # Set knot-vector characteristics
        nb_ctrlpts_eta = degree_eta + 1 
        knot_vector_eta = create_knotvector(degree_eta, 1)

        # Construct points to be interpolated
        theta = np.linspace(0.0, math.pi/2.0, nb_ctrlpts_eta)
        radius = 1.0
        pts_interp = [[radius, 0.0]]
        for angle in theta[1:-1]: 
            x = radius * np.cos(angle)
            y = radius * np.sin(angle)
            pts_interp.append([x, y])

        pts_interp.append([0.0, radius])

        # Do interpolation
        curve = fitting.interpolate_curve(pts_interp, degree_eta)

        # Get control points and knot-vector
        ctrlpts_arc_2 = np.asarray(curve.ctrlpts)

        # -------------------------------------------
        # Fourth part : construction of annulus sector
        # -------------------------------------------
        # Get control points
        ctrlpts = []

        # Over x and z direction 
        for y_arc_2, z_arc_2 in ctrlpts_arc_2:
            for i in range(nb_ctrlpts_xi): 
                for x_arc_1, y_arc_1 in ctrlpts_arc_1:
                    # ctrlpts_x[i] = radius_1
                    x = ctrlpts_x[i] * x_arc_1
                    y = (ctrlpts_x[i] * y_arc_1 + exc)*y_arc_2 - exc
                    z = (ctrlpts_x[i] * y_arc_1 + exc)*z_arc_2
                    ctrlpts.append([x, y, z])

        # Create Volume
        vol = BSpline.Volume()

        # Set degree
        vol.degree_u = degree_xi
        vol.degree_v = degree_nu
        vol.degree_w = degree_eta

        # Set number of control points 
        vol.ctrlpts_size_u, vol.ctrlpts_size_v, vol.ctrlpts_size_w = \
            nb_ctrlpts_xi, nb_ctrlpts_nu, nb_ctrlpts_eta

        # Set control points
        vol.ctrlpts = ctrlpts

        # Set knot-vector
        vol.knotvector_u = knot_vector_xi
        vol.knotvector_v = knot_vector_nu
        vol.knotvector_w = knot_vector_eta

        vol.sample_size = self._sample_size

        return vol

    def create_prisme(self, XY, Height, degree_xi, degree_nu, degree_eta):
        """ Creates a prisme using a quadrilateral as a base.
        The quadrilateral coordinates are given in counterclockwise direction """

        # Position of reference
        x0 = [0, 1, 1, 0]
        y0 = [0, 0, 1, 1]

        # Position of quadrilateral corners
        x1 = XY[:, 0]
        y1 = XY[:, 1]

        # Transformation of control points
        # x1 = ax1 x0 + ax2 y0 + ax3 x0 y0 + ax4
        # y1 = ay1 x0 + ay2 y0 + ay3 x0 y0 + ay4
        T = []
        for i in range(4): 
            T.append([x0[i], y0[i], x0[i] * y0[i], 1])

        ax = np.linalg.solve(T, x1)
        ay = np.linalg.solve(T, y1)

        # Control points
        nb_ctrlpts_xi = degree_xi + 1 
        nb_ctrlpts_nu = degree_nu + 1 
        nb_ctrlpts_eta = degree_eta + 1 
        
        ctrlpts_xi = np.linspace(0, 1, nb_ctrlpts_xi)
        ctrlpts_nu = np.linspace(0, 1, nb_ctrlpts_nu)
        ctrlpts_eta = np.linspace(0, 1, nb_ctrlpts_eta)
        
        ctrlpts = []
        for k in range(nb_ctrlpts_eta):
            for i in range(nb_ctrlpts_xi): 
                for j in range(nb_ctrlpts_nu):    
                    xt = ctrlpts_xi[i]
                    yt = ctrlpts_nu[j]
                    z = ctrlpts_eta[k]*Height

                    x = ax[0]*xt + ax[1]*yt + ax[2]*xt*yt + ax[3]
                    y = ay[0]*xt + ay[1]*yt + ay[2]*xt*yt + ay[3]

                    ctrlpts.append([x, y, z])

        # Get knot-vector
        knot_vector_xi = create_knotvector(degree_xi, 1)
        knot_vector_nu = create_knotvector(degree_nu, 1)
        knot_vector_eta = create_knotvector(degree_eta, 1)

        # Creation of Quadralateral surface
        vol = BSpline.Volume() 

        # Set degree
        vol.degree_u = degree_xi
        vol.degree_v = degree_nu
        vol.degree_w = degree_eta
        
        # Set number of control points 
        vol.ctrlpts_size_u, vol.ctrlpts_size_v, vol.ctrlpts_size_w = nb_ctrlpts_xi, nb_ctrlpts_nu, nb_ctrlpts_eta

        # Set control points
        vol.ctrlpts = ctrlpts

        # Set knot-vector
        vol.knotvector_u = knot_vector_xi
        vol.knotvector_v = knot_vector_nu
        vol.knotvector_w = knot_vector_eta

        vol.sample_size = self._sample_size

        return vol

    # -----------------------------------
    # PLOT 2D GEOMETRY
    # -----------------------------------
    def plot_2D_geometry(self):

        def plot_mesh(pts, shape, ax):

            if pts.shape[1] == 3:
                pts = pts[:, :2]

            pts2D = []
            for j in range(shape[1]):
                pts_temp = []
                for i in range(shape[0]):
                    pos = i + j*shape[0]
                    pts_temp.append(pts[pos, :].tolist())
                pts2D.append(pts_temp)
            pts2D = np.asarray(pts2D)

            # In the first direction
            for _ in range(shape[1]): 
                x = pts2D[_, :, 0]
                y = pts2D[_, :, 1]
                ax.plot(x, y, 'k--')

            # In the second direction
            for _ in range(shape[0]):
                x = pts2D[:, _, 0]
                y = pts2D[:, _, 1]
                ax.plot(x, y, 'k--')

            return

        if self._dim[0] != 2:
            raise Warning('Only can plot 2D parts')

        # Get geometry
        geometry = self._geometry
        N = geometry.sample_size[0]

        # Control points
        ctrlpts = np.asarray(self._ctrlpts)

        # Control points in 2d
        shape_ctrlpts = [1, 1]
        for _ in range(self._dim[0]):
            shape_ctrlpts[_] = self._nb_ctrlpts[_][0]
        shape_ctrlpts = tuple(shape_ctrlpts)

        # Set shape
        shape = [1, 1]
        for _ in range(self._dim[0]):
            shape[_] = N
        shape = tuple(shape)

        # Get eval points
        evalpts_old = np.array(geometry.evalpts)
        evalpts_new = []
        for j in range(N):
            for i in range(N):
                pos = j + i*N
                evalpts_temp = evalpts_old[pos, :]
                evalpts_new.append(evalpts_temp.tolist())

        # Set values
        evalpts_new = np.asarray(evalpts_new) 
        X = np.asarray(evalpts_new[:,0].reshape(shape).tolist())
        Y = np.asarray(evalpts_new[:,1].reshape(shape).tolist())
        Z = np.zeros(X.shape)

        # Plot
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(ctrlpts[:, 0], ctrlpts[:, 1], 'o', markersize=10, label= 'Control points')
        ax.pcolormesh(X, Y, Z, cmap=plt.cm.Pastel1, shading = 'gouraud', label='B-spline surface')
        plot_mesh(ctrlpts, shape_ctrlpts, ax)
        ax.set_xticks(np.arange(0, max(evalpts_new[:,0])+1, 1.0))
        ax.set_yticks(np.arange(0, max(evalpts_new[:,1])+1, 1.0))

        # Set parameters
        ax.axis('equal')
        ax.legend(prop={'size': 14})
        ax.tick_params(axis='both', which='major', labelsize=16)
        ax.set_xlabel('X', fontsize=16)
        ax.set_ylabel('Y', fontsize=16)
        fig.tight_layout()
        
        return fig

def create_geometry(degree, cuts, geometry_case):

    if geometry_case == 'CB': filename = 'parallelepiped'
    elif geometry_case == 'VB': filename = 'prism'
    elif geometry_case == 'TR': filename = 'thick_ring'
    elif geometry_case == 'RQA': filename = 'rotated_quarter_annulus'
    else: raise Warning('Geometry does not exist')

    # Create and refine model
    geometry = {'degree': [degree, degree, degree]}
    modelGeo = geomdlModel(filename=filename, **geometry)
    modelGeo.knot_refinement(nb_refinementByDirection= cuts*np.array([1, 1, 1]))

    return modelGeo
