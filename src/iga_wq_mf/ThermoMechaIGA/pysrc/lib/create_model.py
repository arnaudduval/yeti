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
from .D3viscoplasticity import create_der2sym, fourth_order_identity, one_kron_one 
from .base_functions import eval_basis_fortran
from .create_geomdl import geomdlModel
from iga_wq_mf import assembly

class thermoMechaModel(): 

    def __init__(self, modelIGA: None, material=None, Dirichlet=None, Neumann=None):
        
        print('\nInitializing thermo-mechanical model')

        # Initialize B-Spline properties
        print('Setting B-spline properties')
        self._sample_size = 101
        self._r_ = 2
        self._name = self.read_name(modelIGA)
        self._dim = self.read_dimensions(modelIGA)
        self._degree = self.read_degree(modelIGA)
        self._knotvector, self._size_kv = self.read_knotvector(modelIGA)
        self._ctrlpts = self.read_ControlPoints(modelIGA)
        self._set_parametric_properties()

        # Initialize thermo-mechanical properties
        self._set_material(material)

        # Initialize Dirichlet boundaries
        self._set_dirichlet_boundaries(Dirichlet)

        # Initialize Neumman boundaries
        self._set_neumann_boundaries(Neumann)

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

        try: self._density = material["density"]
        except: self._density = None

        try: self._poissonCoef = material["poisson"]
        except: self._poissonCoef = None

        try: self._youngModule = material["young"]
        except: self._youngModule = None

        try: self._Hardening = material["hardening"]
        except: self._Hardening = None

        try: self._betaHard = material["betahard"]
        except: self._betaHard = None

        try: self._sigmaY = material["sigmaY"]
        except: self._sigmaY = None

        # Initialize others properties
        self._Ctensor = None

        return
    
    def _set_extra_mechanical_properties(self):
        " Compute mechanical properties from E and nu"

        # Initialize
        d = self._dim

        if self._youngModule is None or self._poissonCoef is None:
            pass
        else:
            
            # Create tensor
            identity = fourth_order_identity(d)
            onekronone = one_kron_one(d)

            # Get material properties
            E = self._youngModule
            nu = self._poissonCoef
            lamb = nu*E/((1+nu)*(1-2*nu))
            mu = E/(2*(1+nu))
            bulk = lamb + 2.0/3.0*mu

            # Create material tensor
            Idev = identity - 1.0/3.0*onekronone
            C = lamb*onekronone + 2*mu*identity
            S = 1.0/(9.0*bulk)*onekronone + 1.0/(2.0*mu)*Idev

            # Update
            self._Ctensor = C
            self._Stensor = S
            self._Idev = Idev
            self._lame_lambda = lamb
            self._lame_mu = mu
            self._lame_bulk = bulk

        return

    def _set_dirichlet_boundaries(self, Dirichlet):
        " Gets free and blocked control points "
        
        # Thermal 
        try: 
            TTable = Dirichlet['thermal']
            Tdof, Tdod = self.Dirichlet_boundaries(table= TTable)
            Tdof, Tdod = Tdof[0], Tdod[0]
        except: 
            TTable, Tdof, Tdod = None, None, None
        self._thermalDirichlet = TTable
        self._thermal_dof = Tdof
        self._thermal_dod = Tdod

        # Mechanical 
        try: 
            MTable = Dirichlet['mechanical']
            Mdof, Mdod = self.Dirichlet_boundaries(table= MTable)
        except: 
            MTable, Mdof, Mdod = None, None, None
        self._MechanicalDirichlet = MTable
        self._mechanical_dof = Mdof
        self._mechanical_dod = Mdod

        return TTable, Tdof, Tdod, MTable, Mdof, Mdod

    def _set_neumann_boundaries(self, Neumann):
        " Gets Neumann control points and quadrature points"
        
        # Thermal 
        try: 
            TTable = Neumann['thermal']
        except: 
            TTable = None
        self._thermalNeumann = TTable

        # Mechanical 
        try: 
            MTable = Neumann['mechanical']
        except: 
            MTable = None
        self._MechanicalNeumann = MTable

        return TTable, MTable

    def _clear_material(self): 
        " Clears material"
        self._conductivity = None
        self._capacity = None
        self._poissonCoef = None
        self._youngModule = None
        self._Hardening = None
        self._betaHard = None
        self._sigmaY = None

        return

    def _clear_Dirichlet(self):
        " Clears blocked boundaries "
        self._thermalDirichlet = None
        self._thermal_dof = None
        self._thermal_dod = None
        self._MechanicalDirichlet = None
        self._mechanical_dof = None
        self._mechanical_dod = None
        return

    def _clear_Neumann(self):
        " Clears Neumann boundaries "
        self._thermalNeumann = None
        self._MechanicalNeumann = None
        return

    def _verify_mechanics(self): 
        " Verifies if mechanical properties exits"
        prop = [self._youngModule, self._Hardening, self._betaHard, self._poissonCoef, self._sigmaY]
        if any([var is None for var in prop]): raise Warning('Mechanics not well defined')
        if self._Ctensor is None:
            self._set_extra_mechanical_properties()
        return
    
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
    def get_NURBScoordinates(self, nnz_dim: list): 
        """ Sets topology table: 
        INC: NURBS coordinates
        """

        # Find total number of nnz
        nnz_total = np.prod(nnz_dim)

        # ----------------------
        # INC: NURBS coordinates
        # ----------------------
        INC = np.zeros((nnz_total, 3), dtype= int)

        for i3 in range(nnz_dim[2]): 
            for i2 in range(nnz_dim[1]): 
                for i1 in range(nnz_dim[0]):
                    genPos = i1 + i2*nnz_dim[0] + i3*nnz_dim[0]*nnz_dim[1]
                    INC[genPos, :] = [i1, i2, i3]

        return INC

    def Dirichlet_boundaries(self, table): 

        # The table of dirichlet boundaries must be at least 3D
        table = np.atleast_3d(table)

        # Get number of degree of freedom (DOF) per node
        nbDOF = np.shape(table)[2]
        if np.shape(table)[0] < self._dim or np.shape(table)[1] != 2:
            raise Warning('Table is not well defined')

        # Get total number of control points
        nb_ctrlpts = self._nb_ctrlpts
        nb_ctrlpts_total = self._nb_ctrlpts_total

        # Get topology 
        INC = self.get_NURBScoordinates(nb_ctrlpts)

        # Find nodes
        dod_total = []; dof_total = []
        for i in range(nbDOF):
            dod = []
            for dim in range(self._dim):
                block_bound_dim = table[dim, :, i]
                
                if block_bound_dim[0]: 
                    dod.extend(np.where(INC[:, dim] == 0)[0])

                if block_bound_dim[1]:
                    dod.extend(np.where(INC[:, dim]  == nb_ctrlpts[dim]-1)[0])

            # Rearrange
            dod = np.unique(dod)
            dof = set(np.arange(nb_ctrlpts_total, dtype= int)) - set(dod)
            dod_total.append(list(dod))
            dof_total.append(list(dof))

        return dof_total, dod_total

    def array2coo_matrix(self, data, indi, indj):
        " Computes csr sparse matrix "

        # Computes number of rows and cols
        nb_rows = max(indi) + 1
        nb_cols = max(indj) + 1

        # Set sparse coo matrix
        sparse_matrix = sp.coo_matrix((data, (indi, indj)), shape=(nb_rows, nb_cols))
                                        
        return sparse_matrix

    def array2csr_matrix(self, data, indi, indj):
        " Computes csr sparse matrix "

        # Computes number of rows and cols
        nb_rows = len(indi) - 1
        nb_cols = max(indj) + 1

        # Set sparse csr matrix
        sparse_matrix = sp.csr_matrix((data, indj, indi), shape=(nb_rows, nb_cols))
                                        
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
            for d in range(dim):
                at = alpha[d] 
                B = sp.kron(DB[d][at], B)

            for i in range(dim):
                J[i, j, :] = sp.coo_matrix.dot(B.T, ctrlpts[i, :])
            
        # Evaluate position in physical space
        for i in range(dim):
            B = 1
            for d in range(dim):
                B = sp.kron(DB[d][0], B)

            PPS[i, :] = sp.coo_matrix.dot(B.T, ctrlpts[i, :])

        # Evaluate determinant
        for i in range(nnz):
            detJ[i] = np.linalg.det(J[:, :, i])
        
        stop = time.time()
        print('\tJacobian in : %.5f s' %(stop-start))
        
        return J, PPS, detJ

    def eval_conductivity_coefficient(self, JJ, KK): 
        " Computes coefficients at points P in parametric space "

        print('Getting conductivity and capacity coefficients')
        start = time.time()
        
        # Transform Kprop and Cprop
        KK = np.atleast_3d(KK)
        nnz = np.shape(JJ)[2]
        Kcoef = np.zeros(np.shape(JJ))

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

        stop = time.time()
        print('\tConductivity coefficients in : %.5f s' %(stop-start))

        return Kcoef

    def eval_capacity_coefficient(self, JJ, CC): 
        " Computes coefficients at points P in parametric space "

        print('Getting conductivity and capacity coefficients')
        start = time.time()

        # Transform Kprop and Cprop
        CC = np.atleast_1d(CC)
        nnz = np.shape(JJ)[2]
        Ccoef = np.zeros(nnz)
        
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
        print('\tCapacity coefficients in : %.5f s' %(stop-start))

        return Ccoef

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

    def compute_elastic_coefficient(self, JJ, isnoised=False):
        " Computes coefficients at points P in parametric space. This function only consider linear isotropic case"

        # Set shape
        d = self._dim
        ddl = int(d*(d+1)/2)
        nnz = np.shape(JJ)[2]

        # Create tensors
        MM = create_der2sym(d)
        coefs = np.zeros((d*d, d*d, nnz))
        Gamma = np.zeros((d*d, d*d))
        
        # Create material tensor
        CC = self._Ctensor
        if CC is None: raise Warning('C tensor not define')

        # Initialize
        C = np.zeros((ddl, ddl, nnz))
        for i in range(nnz): 
            C[:, :, i] = CC

        if isnoised:
            # Inset noise by hand
            noise = np.random.normal(loc=-100, scale=13, size=(nnz))
            C[0, 0, :] += noise

            noise = np.random.normal(loc=-50, scale=8, size=(nnz))
            C[1, 1, :] += noise
            C[2, 2, :] += noise

            noise = np.random.normal(loc=-60, scale=12, size=(nnz))
            C[3, 3, :] += noise
            C[4, 4, :] += noise
            C[5, 5, :] += noise

            noise = np.random.normal(loc=50, scale=7, size=(nnz))
            C[0, 1, :] += noise
            C[1, 0, :] += noise
            C[0, 2, :] += noise
            C[2, 0, :] += noise

            noise = np.random.normal(loc=5, scale=5, size=(nnz))
            C[1, 2, :] += noise
            C[2, 1, :] += noise

            # --------------
            noise = np.random.normal(loc=0, scale=8, size=(nnz))
            C[0, 3, :] += noise
            C[0, 5, :] += noise
            C[3, 0, :] += noise
            C[5, 0, :] += noise

            noise = np.random.normal(loc=0, scale=4, size=(nnz))
            C[1, 3, :] += noise
            C[3, 1, :] += noise
            C[2, 5, :] += noise
            C[5, 2, :] += noise

            noise = np.random.normal(loc=0, scale=4, size=(nnz))
            C[2, 3, :] += noise
            C[3, 2, :] += noise
            C[1, 5, :] += noise
            C[5, 1, :] += noise

            noise = np.random.normal(loc=0, scale=3, size=(nnz))
            C[5, 3, :] += noise
            C[3, 5, :] += noise

        for i in range(nnz):

            # Compute M.T D M
            MDM = MM.T @ C[:, :, i] @ MM

            # Compute inverse of JJ
            invJJ = np.linalg.inv(JJ[:, :, i])

            # Create extended version
            for j in range(d):
                Gamma[j*d:(j+1)*d, j*d:(j+1)*d] = invJJ

            # Create coefs
            coefs[:, :, i] = Gamma.T @ MDM @ Gamma 
            coefs[:, :, i] *= np.linalg.det(JJ[:, :, i])

        return coefs

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
            B, indi, indj = eval_basis_fortran(self._degree[dim], self._knotvector[dim], knots)
            DB.append(B)
            ind.append([indi, indj])

        # ==============================
        # Get position and determinant
        # ==============================
        indices, data = [], []
        for dim in range(self._dim):
            indices.append(ind[dim][0])
            indices.append(ind[dim][1])
            data.append(DB[dim])
        inputs = [*self._dim*[nnz], *indices, *data, self._ctrlpts]
        
        if self._dim == 2:
            jacobien_PS, detJ, _ = assembly.eval_jacobien_2d(*inputs)
            qp_PS = assembly.interpolate_fieldphy_2d(*inputs)
        elif self._dim == 3: 
            jacobien_PS, detJ, _ = assembly.eval_jacobien_3d(*inputs)
            qp_PS = assembly.interpolate_fieldphy_3d(*inputs)

        # ==============================
        # Get interpolation
        # ==============================
        if u_ctrlpts is not None:
            indices, data = [], []
            for dim in range(self._dim):
                indices.append(ind[dim][0])
                indices.append(ind[dim][1])
                data.append(DB[dim])
            inputs = [*self._dim*[nnz], *indices, *data, u_ctrlpts]

            if self._dim == 2: u_interp = assembly.interpolate_fieldphy_2d(*inputs)    
            elif self._dim == 3: u_interp = assembly.interpolate_fieldphy_3d(*inputs)
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

        if u_ctrlpts is None: pass
        elif isinstance(u_ctrlpts, np.ndarray): 
            u_ctrlpts = np.atleast_2d(u_ctrlpts)
            if np.size(u_ctrlpts)%self._nb_ctrlpts_total == 0: 
                ndof =  len(u_ctrlpts)
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
        U = np.zeros((ndof, *shape_pts))
        DET = np.zeros(shape_pts)

        for k in range(shape_pts[2]):
            for j in range(shape_pts[1]):
                for i in range(shape_pts[0]):
                    pos = i + j * self._sample_size + k * self._sample_size**2
                    X1[i,j,k] = qp_PS[0, pos]
                    X2[i,j,k] = qp_PS[1, pos]
                    DET[i,j,k] = detJ[pos]
                    if self._dim == 3: X3[i,j,k] = qp_PS[2, pos]
                    if u_interp is not None: 
                        for l in range(ndof):
                            U[l,i,j,k] = u_interp[l, pos]
        
        # Create point data 
        pointData= {"detJ" : DET}

        for l in range(ndof):
            varname = 'U' + str(l+1)
            pointData[varname] = U[l,:,:,:]

        # Export geometry
        name = folder + self._name
        gridToVTK(name, X1, X2, X3, pointData= pointData)
        
        return
