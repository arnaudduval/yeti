"""
.. This module helps to read information from geometries (YETI format)
.. Joaquin Cornejo
"""

# Python libraries
import statistics, time
import numpy as np
from scipy import sparse as sp
from pyevtk.hl import gridToVTK 

# Yeti libraries
from preprocessing.igaparametrization import IGAparametrization

# My libraries
from .base_functions import eval_basis_fortran
from .create_geomdl import geomdlModel
from iga_wq_mf import assembly

class thermoMechaModel(): 

    def __init__(self, modelIGA: IGAparametrization, material=None, Dirichlet=None):
        
        print('\nInitializing thermo-mechanical model')

        # Initialize B-Spline properties
        print('Setting B-spline properties')
        self._sample_size = 101
        self._r_ = 2
        self._dim = self.read_dimensions(modelIGA)
        self._degree = self.read_degree(modelIGA)
        self._knotvector, self._size_kv = self.read_knotvector(modelIGA)
        self._ctrlpts = self.read_ControlPoints(modelIGA)
        self._set_parametric_properties()

        # Initialize thermo-mechanical properties
        self._set_material(material)

        # Initialize Dirichlet boundaries
        self._thermal_dof, self._thermal_dod,\
        self._mechanical_dof, self._mechanical_dod = self._set_blocked_boundaries(Dirichlet)

        return

    def _set_parametric_properties(self):
        " Sets some B-spline properties "

        # Number of elements in each dimension
        self._nb_el = np.ones(3, dtype= int)  
        for dim in range(self._dim):
            self._nb_el[dim] = len(np.unique(self._knotvector[dim])) - 1
        
        # Whether or not Bspline respect IGA-WQ conditions 
        if any(self._nb_el[:self._dim]<2): 
            raise Warning('Model need at least 2 elements in each dimension')

        # Total number of elements
        self._nb_el_total = np.product(self._nb_el)

        # Set number of control points in each dimension
        self._nb_ctrlpts = np.ones(3, dtype= int)  
        for dim in range(self._dim):
            self._nb_ctrlpts[dim] = self._size_kv[dim] - self._degree[dim] - 1

        # Total number of control points
        self._nb_ctrlpts_total = np.product(self._nb_ctrlpts)

        # Set number of quadrature points
        # In WQ approach
        self._nb_qp_wq = np.ones(3, dtype= int) 
        self._nb_qp_wq[:self._dim] = 2*(self._degree[:self._dim] 
                                        + self._nb_el[:self._dim] + self._r_) - 5 
        # In IGA approach
        self._nb_qp_cgg = np.ones(3, dtype= int)
        self._nb_qp_cgg[:self._dim] = (self._degree[:self._dim] + 1)*self._nb_el[:self._dim]

        # Set total number of quadrature points
        self._nb_qp_cgg_total = np.prod(self._nb_qp_cgg)
        self._nb_qp_wq_total = np.prod(self._nb_qp_wq)

        print('B-spline properties :\n\
                -Dimensions: %d\n\
                -Degree: %s\n\
                -Number of elements: %s\n\
                -Number of control points : %s\n\
                -Number of WQ quadrature points: %s\n\
                -Number of Gauss quadrature points: %s\n\
                '
                %(
                self._dim, 
                self._degree,
                self._nb_el,
                self._nb_ctrlpts, 
                self._nb_qp_wq, 
                self._nb_qp_cgg
                )
            )

        return

    def _set_material(self, material): 
        " Define material properties (thermal and mechanical) "

        try: self._conductivity = material["conductivity"]
        except: self._conductivity = None

        try: self._capacity = material["capacity"]
        except: self._capacity = None

        try: self._poissonCoef = material["poisson"]
        except: self._poissonCoef = None

        try: self._youngModule = material["young"]
        except: self._youngModule = None

        return

    def _set_blocked_boundaries(self, Dirichlet):
        " Gets free and blocked control points "
        
        # Thermal 
        try: 
            blockedboundaries = Dirichlet['thermal']
            Tdof, Tdod = self.block_boundaries(table_dirichlet= blockedboundaries)
        except: 
            Tdof, Tdod = None, None

        # Mechanical 
        try: 
            blockedboundaries = Dirichlet['mecanical']
            Mdof, Mdod = self.block_boundaries(table_dirichlet= blockedboundaries)
        except: 
            Mdof, Mdod = None, None

        return Tdof, Tdod, Mdof, Mdod

    # =======================
    # READ FILE
    # =======================
    
    def read_name(self, modelIGA): 
        " Reads name from model "
        if isinstance(modelIGA, IGAparametrization): 
            try: name = modelIGA._name
            except: name = 'IGAparametrization'
        elif isinstance(modelIGA, geomdlModel): name = modelIGA._name
        return name

    def read_degree(self, modelIGA): 
        " Reads degree from model "
        if isinstance(modelIGA, IGAparametrization): degree = modelIGA._Jpqr.flatten()
        elif isinstance(modelIGA, geomdlModel): degree = modelIGA._degree
        if any(p == 1 for p in degree[:self._dim]): 
            raise Warning('Model must have at least degree p = 2')
        return degree

    def read_dimensions(self, modelIGA):
        " Reads dimensions from model"
        if isinstance(modelIGA, IGAparametrization): dim = modelIGA._dim[0]
        elif isinstance(modelIGA, geomdlModel): dim = modelIGA._dim
        if dim != 3: raise Warning('Model must be 3D')
        return dim

    def read_knotvector(self, modelIGA):
        " Reads knot-vector from model"
        if isinstance(modelIGA, IGAparametrization): 
            knotvector = modelIGA._Ukv[0]
            size_kv = modelIGA._Nkv.flatten()
        elif isinstance(modelIGA, geomdlModel): 
            knotvector = modelIGA._knotvector
            size_kv = modelIGA._size_kv
        return knotvector, size_kv

    def read_ControlPoints(self, modelIGA): 
        " Reads control points from model"
        if isinstance(modelIGA, IGAparametrization): 
            ctrlpts = modelIGA._COORDS[:self._dim, :]
        elif isinstance(modelIGA, geomdlModel): 
            ctrlpts = modelIGA._ctrlpts
        return ctrlpts

    # ========================
    # LOCAL FUNCTIONS
    # ========================

    def block_boundaries(self, table_dirichlet): 

        def get_NURBScoordinates(): 
            """ Sets topology table: 
            INC: NURBS coordinates
            """

            # Get number of control points in each dimension
            nb_ctrlpts = self._nb_ctrlpts

            # Find total number of control points 
            nb_ctrlpts_total = self._nb_ctrlpts_total

            # ----------------------
            # INC: NURBS coordinates
            # ----------------------
            INC = np.zeros((nb_ctrlpts_total, 3), dtype= int)

            for i3 in range(nb_ctrlpts[2]): 
               for i2 in range(nb_ctrlpts[1]): 
                   for i1 in range(nb_ctrlpts[0]):
                       genPos = i1 + i2*nb_ctrlpts[0] + i3*nb_ctrlpts[0]*nb_ctrlpts[1]
                       INC[genPos, :] = [i1, i2, i3]

            return INC

        # The table of dirichlet boundaries must be at least 3D
        table_dirichlet = np.atleast_3d(table_dirichlet)

        # Get number of degree of freedom (DOF) per node
        nbDOF = np.shape(table_dirichlet)[2]
        if np.shape(table_dirichlet)[0] < self._dim or np.shape(table_dirichlet)[1] != 2:
            raise Warning('Table is not well defined')

        # Get total number of control points
        nb_ctrlpts_total = self._nb_ctrlpts_total

        # Get topology 
        INC = get_NURBScoordinates()

        # Find Dirichlet nodes
        nodes_dir_total = []
        for dof in range(nbDOF):
            nodes_dir = []
            for dim in range(self._dim):
                block_bound_dim = table_dirichlet[dim, :, dof]
                
                if block_bound_dim[0]: 
                    nodes_dir.extend(np.where(INC[:, dim] == 0)[0])

                if block_bound_dim[1]:
                    nodes_dir.extend(np.where(INC[:, dim]  == self._nb_ctrlpts[dim]-1)[0])

            # Rearrange
            nodes_dir = np.unique(nodes_dir)
            nodes_dir += dof*nb_ctrlpts_total*np.ones(len(nodes_dir), dtype= int)
            nodes_dir_total.extend(list(nodes_dir))

        # Find blocked nodes
        dod = nodes_dir_total

        # Find equations nodes
        dof = set(np.arange(nb_ctrlpts_total, dtype= int)) - set(dod)
        dof = list(dof)

        return dof, dod

    def array2coo_matrix(self, data, indi, indj):
        " Computes coo sparse matrix "

        # Computes number of rows and cols
        nb_rows = len(indi) - 1
        nb_cols = max(indj) + 1

        # Set sparse coo matrix
        sparse_matrix = sp.coo_matrix((data, (indi, indj)), 
                                        shape=(nb_rows, nb_cols))
        return sparse_matrix

    def array2csr_matrix(self, data, indi, indj):
        " Computes csr sparse matrix "

        # Computes number of rows and cols
        nb_rows = len(indi) - 1
        nb_cols = max(indj) + 1

        # Set sparse coo matrix
        sparse_matrix = sp.csr_matrix((data, indj, indi), 
                                        shape=(nb_rows, nb_cols))
                                        
        return sparse_matrix

    # ===========================
    # SPECIAL  FUNCTIONS 
    # ===========================

    def eval_jacobien_physicalPosition(self, dim, nnz, ctrlpts, DB): 
        """ Computes Jacobien matrix and find the position in physical space 
        of P (in pts list) in parametric space
        """
        print('Evaluating jacobien and physical position')
        start = time.time()
        
        # Initialize 
        J = np.zeros((dim, dim, nnz))
        PPS = np.zeros((dim, nnz))
        detJ = np.zeros(nnz)

        # Evaluate jacobien
        for j in range(dim):
            alpha = np.zeros(dim, dtype = int); alpha[j] = 1

            B = 1
            for dim in range(dim):
                at = alpha[dim] 
                B = sp.kron(DB[dim][at], B)

            for i in range(dim):
                J[i, j, :] = sp.coo_matrix.dot(B.T, ctrlpts[i, :])
            
        # Evaluate position in physical space
        for i in range(dim):
            B = 1
            for dim in range(dim):
                B = sp.kron(DB[dim][0], B)

            PPS[i, :] = sp.coo_matrix.dot(B.T, ctrlpts[i, :])

        # Evaluate determinant
        for i in range(nnz):
            detJ[i] = np.linalg.det(J[:, :, i])
        
        stop = time.time()
        print('\tJacobian in : %.5f s' %(stop-start))
        
        return J, PPS, detJ

    def eval_thermal_coefficient(self, nnz, JJ, KK, CC): 
        " Computes coefficients at points P in parametric space "

        print('Getting conductivity and capacity coefficients')
        start = time.time()

        Kcoef = np.zeros(np.shape(JJ))
        Ccoef = np.zeros(nnz)

        # Transform Kprop and Cprop
        KK = np.atleast_3d(KK)
        CC = np.atleast_1d(CC)

        if np.shape(KK)[2] == 1:
            for i in range(nnz): 
                # Find determinant of Jacobien 
                det_J = np.linalg.det(JJ[:, :, i])

                # Find inverse of Jacobien 
                inv_J = np.linalg.inv(JJ[:, :, i])

                # Find coefficient of conductivity matrix
                Kcoef[:, :, i] = inv_J @ KK[:, :, 0] @ inv_J.T * det_J

        elif np.shape(KK)[2] == nnz:
            for i in range(nnz): 
                # Find determinant of Jacobien 
                det_J = np.linalg.det(JJ[:, :, i])

                # Find inverse of Jacobien 
                inv_J = np.linalg.inv(JJ[:, :, i])

                # Find coefficient of conductivity matrix
                Kcoef[:, :, i] = inv_J @ KK[:, :, i] @ inv_J.T * det_J

        else: 
            raise Warning('Something happen, it is not possible to compute coefficients')

        if len(CC) == 1:
            for i in range(nnz): 
                # Find determinant of Jacobien 
                det_J = np.linalg.det(JJ[:, :, i])

                # Find coefficient of capacity matrix or heat vector
                Ccoef[i] = CC[0] * det_J

        elif len(CC) == nnz:
            for i in range(nnz): 
                # Find determinant of Jacobien 
                det_J = np.linalg.det(JJ[:, :, i])

                # Find coefficient of capacity matrix or heat vector
                Ccoef[i] = CC[i] * det_J

        else: 
            raise Warning('Something happen, it is not possible to compute coefficients')

        stop = time.time()
        print('\tConductivity and capacity coefficients in : %.5f s' %(stop-start))

        return Kcoef, Ccoef

    def eval_source_coefficient(self, fun, qp, det): 
        " Computes coefficients at points P in parametric space "

        print('Getting source coefficients')
        start = time.time()
        # Get source coefficient
        qp = np.atleast_2d(qp)
        source_coef = fun(qp) * det
        stop = time.time()
        print('\tSource coefficients in : %.5f s' %(stop-start))

        return source_coef

    def solver_scipy(self, A, b, nbIterations=100, epsilon=1e-10, PreCond='ilu', isCG=True):
        "Solve system using conjugate gradient or Bi-conjugate gradient. Matrix A must be assembled"

        # Find preconditionner
        if PreCond == 'ilu': 
            B = sp.linalg.spilu(A)
            Mx = lambda x: B.solve(x)
            M = sp.linalg.LinearOperator(A.shape, Mx)
        # ADD NEW METHODS ...

        # Solve with iterative method
        if isCG: x, info = sp.linalg.cg(A, b, tol=epsilon, maxiter=nbIterations, M=M)
        else: x, info = sp.linalg.bicgstab(A, b, tol=epsilon, maxiter=nbIterations, M=M)

        return x

    # ===========================
    # POST-PROCESSING 
    # ===========================

    def interpolate_field(self, nnz= None, u_ctrlpts= None):
        "Interpolates the input field. In all cases, it returns jacobien."

        if nnz == None: nnz = self._sample_size

        # =====================
        # Get Basis
        # =====================
        # Define knots
        knots = np.linspace(0, 1, nnz)

        # Set basis and indices
        DB, ind = [], []
        for dim in range(self._dim):  
            B0, B1, indi, indj = eval_basis_fortran(self._degree[dim], self._knotvector[dim], knots)
            DB.append([B0, B1])
            ind.append([indi, indj])

        # ==============================
        # Get position and determinant
        # ==============================
        indices, data, ctrlpts = [], [], []
        for dim in range(self._dim):
            indices.append(ind[dim][0])
            indices.append(ind[dim][1])
            data.append(DB[dim][0])
            data.append(DB[dim][1])
            ctrlpts.append(self._ctrlpts[dim, :])
        inputs = [self._dim*[nnz], *indices, *data, *ctrlpts]

        if self._dim == 2: jacobien_PS, qp_PS, detJ = assembly.jacobien_physicalposition_2d(*inputs)    
        elif self._dim == 3: jacobien_PS, qp_PS, detJ = assembly.jacobien_physicalposition_3d(*inputs)

        # ==============================
        # Get interpolation
        # ==============================
        if u_ctrlpts is not None:
            indices, data = [], []
            for dim in range(self._dim):
                indices.append(ind[dim][0])
                indices.append(ind[dim][1])
                data.append(DB[dim][0])
            inputs = [self._dim*[nnz], *indices, *data, u_ctrlpts]

            if self._dim == 2: u_interp = assembly.interpolation_2d(*inputs)    
            elif self._dim == 3: u_interp = assembly.interpolation_3d(*inputs)
        else: 
            u_interp = None

        return jacobien_PS, qp_PS, detJ, u_interp
    
    def export_results(self, u_ctrlpts= None, folder=None): 
        " Returns solution using geometry basis "

        if folder == None: 
            import os
            full_path = os.path.realpath(__file__)
            dirname = os.path.dirname
            folder = dirname(dirname(full_path)) + '/results/'
            if not os.path.isdir(folder): os.mkdir(folder)

        if u_ctrlpts == None: pass
        elif isinstance(u_ctrlpts, np.ndarray): 
            if len(u_ctrlpts) == self._nb_ctrlpts_total: pass
            else: raise Warning('Not enough control points')
        else: raise Warning('Solution must be ndarray type')

        # Set shape
        shape_pts = [1, 1, 1]
        for dim in range(self._dim): shape_pts[dim] = self._sample_size
        shape_pts = tuple(shape_pts)

        # ==============================
        # Get interpolation
        # ==============================
        # Interpolate 
        _, qp_PS, detJ, u_interp = self.interpolate_field(u_ctrlpts=u_ctrlpts)
        mean_detJ = statistics.mean(detJ)
        detJ /= mean_detJ

        # ==============================
        # Export results
        # ==============================
        X1 = np.zeros(shape_pts)
        X2 = np.zeros(shape_pts)
        X3 = np.zeros(shape_pts)
        U = np.zeros(shape_pts)
        DET = np.zeros(shape_pts)

        for k in range(shape_pts[2]):
            for j in range(shape_pts[1]):
                for i in range(shape_pts[0]):
                    pos = i + j * self._sample_size + k * self._sample_size**2
                    X1[i,j,k] = qp_PS[0, pos]
                    X2[i,j,k] = qp_PS[1, pos]
                    DET[i,j,k] = detJ[pos]
                    if self._dim == 3: X3[i,j,k] = qp_PS[2, pos]
                    if u_interp != None: U[i,j,k] = u_interp[pos]
                    
        # Export geometry
        name = folder + self._name
        gridToVTK(name, X1, X2, X3, pointData= {"U1" : U, 
                                                "detJ" : DET,})
        return
