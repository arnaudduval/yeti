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
import os, copy, time, math, numpy as np

# My libraries
from .base_functions import create_knotvector

# YETI libraries 
from preprocessing.igaparametrization import IGAparametrization

class geomdlModel(): 

    def __init__(self, name= None, **geometry: None): 

        if name is None:
            raise Warning('Insert the name of the part')
        else : 
            # Set name of object
            self._name = name

        # Set number of samples
        self._sample_size = 101

        print('\nCreating geometry: ' + name + '...')
        start = time.time()
        if name == 'quarter_annulus' or name == 'QA':
            # Set number of dimension 
            self._dim = 2

            # Create quarter annulus
            Rin = geometry.get('Rin', 1.0)
            Rout = geometry.get('Rout', 2.0)
            degree_xi, degree_nu, _ = geometry.get('degree', [2, 3, 2])
            self._geometry = self.create_quarter_annulus(Rin, Rout, degree_xi, degree_nu) 

        elif name == 'quadrilateral' or name == 'SQ':
            # Set number of dimension
            self._dim = 2
            
            # Create parallelepiped
            XY = geometry.get('XY', np.array([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]))
            degree_xi, degree_nu, _ = geometry.get('degree', [2, 2, 2])
            self._geometry = self.create_quadrilateral(XY, degree_xi, degree_nu) 

        elif name == 'parallelepiped' or name == 'CB': 
            # Set number of dimension
            self._dim = 3
            
            # Create parallelepiped
            Lx = geometry.get('Lx', 1.0)
            Ly = geometry.get('Ly', 1.0)
            Lz = geometry.get('Lz', 1.0)
            degree_xi, degree_nu, degree_eta = geometry.get('degree', [2, 2, 2])
            self._geometry = self.create_parallelepiped(Lx, Ly, Lz, degree_xi, degree_nu, degree_eta) 

        elif name == 'thick_ring' or name == 'TR':
            # Set number of dimension 
            self._dim = 3

            # Create quarter annulus
            Rin = geometry.get('Rin', 1.0)
            Rout = geometry.get('Rout', 2.0)
            Height = geometry.get('Height', 1.0)
            degree_xi, degree_nu, degree_eta = geometry.get('degree', [4, 4, 4])
            self._geometry = self.create_thick_ring(Rin, Rout, Height, degree_xi, degree_nu, degree_eta) 

        elif name == 'rotated_quarter_annulus' or name == 'RQA':
            # Set number of dimension 
            self._dim = 3

            # Create quarter annulus
            Rin = geometry.get('Rin', 1.0)
            Rout = geometry.get('Rout', 2.0)
            exc = geometry.get('Rout', 1.0) # excentricity
            degree_xi, degree_nu, degree_eta = geometry.get('degree', [4, 4, 4])
            self._geometry = self.create_rotated_quarter_annulus(Rin, Rout, exc, degree_xi, degree_nu, degree_eta) 

        elif name == 'prism' or name == 'VB':
            # Set number of dimension
            self._dim = 3
            
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

    def update_geometry(self): 
        start = time.time()
        # Set geometry object
        obj = self._geometry
        if obj is None: raise Warning('Geometry unknown')

        # Set degree
        self._degree = np.ones(3, dtype= int)
        self._degree[:self._dim] = np.array(obj._degree)
        if any(p == 1 for p in self._degree[:self._dim]): 
           raise Warning('Model must have at least degree p = 2')

        # Set knot vector
        self._knotvector = [np.array(obj._knot_vector[dim]) for dim in range(self._dim)]

        # Set size knot-vector in each dimension
        self._size_kv = np.ones(3, dtype= int)  
        for dim in range(self._dim):
            self._size_kv[dim] = np.size(self._knotvector[dim])

        # Set number of elements in each dimension
        self._nb_el = np.ones(3, dtype= int)  
        for dim in range(self._dim):
            self._nb_el[dim] = len(np.unique(self._knotvector[dim])) - 1

        # Set total number of elements
        self._nb_el_total = np.product(self._nb_el)

        # Set number of quadrature points/functions in each dimension
        self._nb_ctrlpts = np.ones(3, dtype= int)  
        for dim in range(self._dim):
            self._nb_ctrlpts[dim] = self._size_kv[dim] - self._degree[dim] - 1
        self._nb_ctrlpts_total = np.product(self._nb_ctrlpts)

        if self._dim == 2: 
            c = 0
            ctrlpts_old = obj._control_points
            ctrlpts_new = np.zeros((self._dim, self._nb_ctrlpts_total))
            for j in range(self._nb_ctrlpts[1]):
                for i in range(self._nb_ctrlpts[0]):
                    pos = j + i*self._nb_ctrlpts[1]
                    ctrlpts_temp = ctrlpts_old[pos]
                    ctrlpts_new[:, c] = ctrlpts_temp
                    c += 1 

        elif self._dim == 3: 
            c =  0
            ctrlpts_old = obj._control_points
            ctrlpts_new = np.zeros((self._dim, self._nb_ctrlpts_total))
            for k in range(self._nb_ctrlpts[2]):
                for j in range(self._nb_ctrlpts[1]):
                    for i in range(self._nb_ctrlpts[0]):
                        pos = j + i*self._nb_ctrlpts[1] + k*self._nb_ctrlpts[1]*self._nb_ctrlpts[0]
                        ctrlpts_temp = ctrlpts_old[pos]
                        ctrlpts_new[:, c] = ctrlpts_temp
                        c += 1 

        # Update control points
        self._ctrlpts = ctrlpts_new
        stop = time.time()
        print('\tGeometry properties updated in: %.3e s\n' %(stop-start))

        return

    def write_abaqus_file(self, filename):
        "Returns the inp and NB file. It only works with one patch"

        def array2txt(array: np.array, format= '%.2f'):
            return ','.join([format %(i) for i in array])

        # ------------
        # With inp file
        # ------------
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

        # ------------
        # With NB file
        # ------------
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

    def knot_refinement(self, nb_refinementByDirection= np.array([0,0,0])):
        "Refine geometry following each direction. It is slow because it uses python methods"

        start = time.time()
        # Copy geometry
        geometry = copy.deepcopy(self._geometry)

        # Set new number of elements
        cuts = nb_refinementByDirection
        nbel = [2**cuts[dim]*self._nb_el[dim] for dim in range(self._dim)]

        # Create knots to be inserted
        knotvector_insert = [np.linspace(0.0, 1.0, i+1)[1:-1] for i in nbel]

        # Insert knots
        for dim in range(self._dim):
            multiplicity = np.zeros(self._dim, dtype= int)
            multiplicity[dim] = 1
            for knot in knotvector_insert[dim]: 
                knot_insert = np.zeros(self._dim)
                knot_insert[dim] = knot
                operations.insert_knot(geometry, knot_insert, multiplicity.tolist())
        self._geometry = geometry
        stop = time.time()
        print('Knot refinement in: %.3e s' %(stop-start))

        # Update values
        self.update_geometry()

        return

    def export_IGAparametrization(self, nb_refinementByDirection= np.array([0,0,0])):
        "Refine geometry following each direction. It is very fast"

        # Write file abaqus
        self.write_abaqus_file(filename=self._name)

        # Discretisize the part
        modelIGA = IGAparametrization(filename= self._name)
        modelIGA.refine(nb_refinementByDirection=nb_refinementByDirection)

        # Clean files created
        os.remove(self._name+'.inp')
        os.remove(self._name+'.NB')
        os.remove(self._name+'.save')

        return modelIGA

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

        # Create a B-spline volume
        vol = BSpline.Volume()

        # Set degree
        vol.degree_u = degree_xi
        vol.degree_v = degree_nu
        vol.degree_w = degree_eta

        # Set number of control points 
        vol.ctrlpts_size_u = int(nb_ctrlpts_xi)
        vol.ctrlpts_size_v = int(nb_ctrlpts_nu)
        vol.ctrlpts_size_w = int(nb_ctrlpts_eta)

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

        # Create volume
        vol = BSpline.Volume()

        # Set degree
        vol.degree_u = degree_xi
        vol.degree_v = degree_nu
        vol.degree_w = degree_eta

        # Set number of control points 
        vol.ctrlpts_size_u, vol.ctrlpts_size_v, vol.ctrlpts_size_w = \
            int(nb_ctrlpts_xi), int(nb_ctrlpts_nu), int(nb_ctrlpts_eta)

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
                    y = (ctrlpts_x[i] * y_arc_1 + exc)*y_arc_2 
                    z = (ctrlpts_x[i] * y_arc_1 + exc)*z_arc_2
                    ctrlpts.append([x, y, z])

        # Create volume
        vol = BSpline.Volume()

        # Set degree
        vol.degree_u = degree_xi
        vol.degree_v = degree_nu
        vol.degree_w = degree_eta

        # Set number of control points 
        vol.ctrlpts_size_u, vol.ctrlpts_size_v, vol.ctrlpts_size_w = \
            int(nb_ctrlpts_xi), int(nb_ctrlpts_nu), int(nb_ctrlpts_eta)

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
        vol.ctrlpts_size_u, vol.ctrlpts_size_v, vol.ctrlpts_size_w = \
        int(nb_ctrlpts_xi), int(nb_ctrlpts_nu), int(nb_ctrlpts_eta)

        # Set control points
        vol.ctrlpts = ctrlpts

        # Set knot-vector
        vol.knotvector_u = knot_vector_xi
        vol.knotvector_v = knot_vector_nu
        vol.knotvector_w = knot_vector_eta

        vol.sample_size = self._sample_size

        return vol

