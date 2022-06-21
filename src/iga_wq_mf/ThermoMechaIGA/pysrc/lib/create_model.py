"""
.. This module helps to read information from geometries (YETI format and Geomdl)
.. Joaquin Cornejo
"""

# Python libraries
import statistics
import numpy as np
from scipy import sparse as sp
import time
from pyevtk.hl import gridToVTK 

# Yeti libraries
from preprocessing.igaparametrization import IGAparametrization

# My libraries
from .base_functions import eval_basis_fortran
from .create_geomdl import geomdlModel
from iga_wq_mf import assembly

class thermoMechaModel(): 

    def __init__(self, modelIGA= None, isThermal= True, isMechanical=False, 
                thermalblockedboundaries= None, mechablockedboundaries= None):

        print('\nInitializing thermo-mechanical model')
        if isinstance(modelIGA, IGAparametrization):
            self._geometry_type = 'yeti'
            self._name = 'IGAparametrization'
        elif isinstance(modelIGA, geomdlModel):
            self._geometry_type = 'geomdl'
            self._name = modelIGA._name
        else: raise Warning("Type unknown")

        # Set number of samples
        try: self._sample_size = modelIGA._sample_size
        except: self._sample_size = 101

        # Set type of model
        self._isThermal = isThermal
        self._isMechanical = isMechanical

        # Initialize b-spline properties
        print('Setting B-spline properties')
        self.__set_bspline_properties(modelIGA)
        
        if self._isThermal:
            # Initialize thermal properties
            print('Settig thermal properties')
            self.__set_thermal_properties(conductivity=np.eye(self._dim), capacity= 1.)

            # Set free and bloqued nodes
            print('Settig free and blocked control points')
            if thermalblockedboundaries == None:
                thermalblockedboundaries = [[1, 1], [1, 1], [1, 1]]
                print('Dirichlet not defined. Default: all blocked')
            self._thermalblockedboundaries = np.asarray(thermalblockedboundaries)

        if self._isMechanical:
            # Initialize elastic properties
            print('Settig mechanical properties')
            self.__set_elastic_properties(modelIGA, E=210e3, nu=0.3)

            # Set free and bloqued nodes
            print('Settig free and blocked control points')
            if mechablockedboundaries == None:
                mechablockedboundaries = [[1, 1], [1, 1], [1, 1]]
                print('Dirichlet not defined. Default: all blocked')
            self._mechablockedboundaries = np.asarray(mechablockedboundaries)

        if self._isMechanical and self._isThermal:
            # Initialize thermal exapansion
            print('Setting thermo-mechanical properties')
            self._thermalexpansion = 9.0e-6

        # Get nodes
        self.__get_free_blocked_nodes()

        return

    def __set_thermal_properties(self, conductivity=None, capacity=None):
        " Set thermal properties " 

        # (Later must be set in .inp file)
        # Assuming isotropic behavior :
        # To be changed or modified: the way we get material properties
        self._conductivity = conductivity
        self._capacity = capacity
        
        print('Thermal properties (isotropic behavior) :\n\
                -Conductivity: %s\n\
                -Capacity: %s' 
                %(
                self._conductivity, 
                self._capacity
                )
            )

        return

    def __set_elastic_properties(self, E=None, nu=None, isPlaneStress= True):
        " Set elasto-mechanical properties " 

        # Assuming isotropic behavior :
        self._youngModule = E # MPa
        self._poissonCoef = nu

        print(
                'Mechanical properties (isotropic behavior) :\n\
                -Young module: %s\n\
                -Poisson coefficient: %f' 
                %(
                self._youngModule, 
                self._poissonCoef
                )
            )

        if self._dim == 2:
            if not isPlaneStress: 
                # Plane strain
                ElasticMatrix = [[1-nu, nu, 0],
                                [nu, 1-nu, 0], 
                                [0, 0, 0.5-nu]]
            elif isPlaneStress:
                # Plane stress
                ElasticMatrix = [[1, nu, 0],
                                [nu, 1, 0], 
                                [0, 0, 1-nu]]

        if self._dim == 3:
            ElasticMatrix = [[1-nu, nu, nu, 0, 0, 0],
                            [nu, 1-nu, nu, 0, 0, 0], 
                            [nu, nu, 1-nu, 0, 0, 0], 
                            [0, 0, 0, 0.5-nu, 0, 0], 
                            [0, 0, 0, 0, 0.5-nu, 0], 
                            [0, 0, 0, 0, 0, 0.5-nu]]

        self._ElasticMatrix = np.asarray(ElasticMatrix)

        return

    def __set_bspline_properties(self, modelIGA): 
        " Sets some B-spline properties "

        # Set some variables we will use :
        # --------------------------------- 
        # Degree in each dimension
        if isinstance(modelIGA, IGAparametrization): self._degree = modelIGA._Jpqr
        elif isinstance(modelIGA, geomdlModel): self._degree = modelIGA._degree

        # Dimensions
        self._dim = modelIGA._dim[0]
        if any(p == 1 for p in self._degree[:self._dim]): 
            raise Warning('Model must have at least degree p = 2')

        # knot-vector in each dimension
        if isinstance(modelIGA, IGAparametrization): self._knotvector = modelIGA._Ukv
        elif isinstance(modelIGA, geomdlModel): self._knotvector = modelIGA._knotvector

        # Size knot-vector in each dimension
        if isinstance(modelIGA, IGAparametrization): self._size_kv = modelIGA._Nkv
        elif isinstance(modelIGA, geomdlModel): self._size_kv = modelIGA._size_kv

        # Number of elements in each dimension
        self._nb_el = np.zeros((3, 1), dtype= int)  
        if isinstance(modelIGA, IGAparametrization): 
            self._nb_el[:self._dim, 0] = self._size_kv[:self._dim, 0]\
                                    - (2*self._degree[:self._dim, 0] + 1)
        elif isinstance(modelIGA, geomdlModel):
             self._nb_el[:self._dim, 0] = modelIGA._nb_el[:self._dim, 0]

        # Total number of elements
        self._nb_el_total = np.product(self._nb_el[:self._dim][:])

        # Control points
        if isinstance(modelIGA, IGAparametrization): 
            self._ctrlpts = modelIGA._COORDS[:self._dim, :].T
        elif isinstance(modelIGA, geomdlModel): 
            self._ctrlpts = np.asarray(modelIGA._ctrlpts)

        if isinstance(modelIGA, IGAparametrization): 
            self._nb_ctrlpts_total = modelIGA._nb_cp
        elif isinstance(modelIGA, geomdlModel):
            self._nb_ctrlpts_total = modelIGA._nb_ctrlpts_total

        # Set number of quadrature points/functions in each dimension
        self._nb_ctrlpts = self._degree + self._nb_el

        print('B-spline properties :\n\
                -Dimensions: %d\n\
                -Degree: %s\n\
                -Number of elements: %s\n\
                -Number of control points : %s'
                %(
                self._dim, 
                self._degree.reshape(1, -1),
                self._nb_el.reshape(1, -1),
                self._nb_ctrlpts.reshape(1, -1)
                )
            )

        # Whether or not Bspline respect IGA-WQ conditions 
        if any(self._nb_el[:self._dim]<2): 
            raise Warning('Model need at least 2 elements in each dimension')

        # Set quadrature points
        self.__set_qp_properties()

        return

    def __set_qp_properties(self):
        " Sets some WQ properties "

        # Set nulber of additional points at the end
        self._r_ = 2

        # Set number of quadrature points
        # In WQ approach
        self._nb_qp_wq = np.zeros((3, 1), dtype= int) 
        self._nb_qp_wq[:self._dim, 0] = 2*(self._degree[:self._dim, 0] 
                                            + self._nb_el[:self._dim, 0] + self._r_) - 5 
        # In IGA approach
        self._nb_qp_cgg = np.zeros((3, 1), dtype= int)
        self._nb_qp_cgg[:self._dim, 0] = (self._degree[:self._dim, 0] + 1)\
                                            *self._nb_el[:self._dim, 0]

        # Set total number of quadrature points
        self._nb_qp_cgg_total = np.prod(self._nb_qp_cgg[:self._dim, 0])
        self._nb_qp_wq_total = np.prod(self._nb_qp_wq[:self._dim, 0])

        print(' -Number of WQ quadrature points: %s\n\
                -Number of Gauss quadrature points: %s'
                %(
                    self._nb_qp_wq.reshape(1, -1), 
                    self._nb_qp_cgg.reshape(1, -1)
                )
            )

        return
    
    def __get_free_blocked_nodes(self):
        " Gets free and blocked control points "

        # ----------------------------
        # Thermal points
        # ----------------------------
        if self._isThermal:
            # Dirichlet 
            dof, dod = self.block_boundaries(blockedboundaries= self._thermalblockedboundaries, typeEl='T')

            # Set 
            self._thermal_dof = dof
            self._thermal_dod = dod
            
        # ----------------------------
        # Mechanical points
        # ----------------------------
        if self._isMechanical:
            # Dirichlet 
            dof, dod = self.block_boundaries(blockedboundaries= self._mechablockedboundaries, typeEl='D')

            # Set 
            self._mechanical_dof = dof
            self._mechanical_dod = dod
            
        return

    # ===========================
    # LOCAL FUNCTIONS 
    # ===========================

    def array2coo_matrix(self, nb_rows, nb_cols, data, indi, indj):
        " Computes coo sparse matrix "

        # Set sparse coo matrix
        sparse_matrix = sp.coo_matrix((data, (indi, indj)), 
                                        shape=(nb_rows, nb_cols))
        return sparse_matrix

    def array2csr_matrix(self, nb_rows, nb_cols, data, indi, indj):
        " Computes csr sparse matrix "

        # Set sparse coo matrix
        sparse_matrix = sp.csr_matrix((data, indj, indi), 
                                        shape=(nb_rows, nb_cols))
                                        
        return sparse_matrix

    def block_boundaries(self, blockedboundaries= None, typeEl= None): 
        # This function it will not be implemented in future code
        # For now, we consider Dirichlet condition 

        def get_NURBScoordinates(): 
            """ Sets topology table: 
            INC: NURBS coordinates
            """

            # Get number of dimensions
            dimensions = self._dim

            # Get number of control points in each dimension
            nb_ctrlpts = np.ones(3, dtype= int)
            for dim in range(dimensions):
                nb_ctrlpts[dim] = self._nb_ctrlpts[dim] 

            # Find total number of control points 
            nb_ctrlpts_total = self._nb_ctrlpts_total

            # ----------------------
            # INC: NURBS coordinates
            # ----------------------
            INC = np.zeros((nb_ctrlpts_total, 3), dtype= int)

            for i3 in range(nb_ctrlpts[2]): 
               for i2 in range(nb_ctrlpts[1]): 
                   for i1 in range(nb_ctrlpts[0]):
                       genPos = i1 + i2*(nb_ctrlpts[0]) + i3*(nb_ctrlpts[0]*nb_ctrlpts[1])
                       INC[genPos, :] = [i1, i2, i3]

            return INC

        # Get number of dimensions
        dimensions = self._dim

        # Get number of control points in each dimension
        nb_ctrlpts = self._nb_ctrlpts

        # Find total number of control points 
        nb_ctrlpts_total = self._nb_ctrlpts_total

        # Verify if list of bloeck bounds has at least the same dimension of geometry
        if len(blockedboundaries) < dimensions:
            raise Warning('Not enough data')

        # Get topology 
        INC = get_NURBScoordinates()

        # Find Dirichlet nodes
        nd = []
        for dim in range(dimensions):
            block_bound_dim = blockedboundaries[dim]
            
            if block_bound_dim[0]: 
                nd.extend(np.where(INC[:, dim] == 0)[0])

            if block_bound_dim[1]:
                nd.extend(np.where(INC[:, dim]  == nb_ctrlpts[dim][0]-1)[0])

        nd = np.unique(nd).tolist()

        if typeEl == 'T': 
            # Find blocked nodes
            dod = nd

            # Find equations nodes
            dof = set(np.arange(nb_ctrlpts_total, dtype= int)) - set(dod)
            dof = list(dof)

            return dof, dod

        elif typeEl == 'D':
            # Find blocked nodes
            dod = nd
            if dimensions == 2:
                new_nd = nd + self._nb_ctrlpts_total
                dod.extend(new_nd)
            elif dimensions == 3:
                new_nd = np.array(nd) + self._nb_ctrlpts_total
                dod.extend(new_nd.tolist())
                new_nd = np.array(nd) + 2*self._nb_ctrlpts_total
                dod.extend(new_nd.tolist())

            # Find equations nodes
            dof = set(np.arange(dimensions*nb_ctrlpts_total, dtype= int)) - set(dod)
            dof = list(dof)

            return dof, dod

        else: 
            return None

    # ===========================
    # ASSEMBLY  FUNCTIONS 
    # ===========================

    def eval_jacobien_physicalPosition(self, dimensions, ctrlpts, data_DB, nb_pts): 
        """ Computes Jacobien matrix and find the position in physical space 
        of P (in pts list) in parametric space
        """
        print('Evaluating jacobien and physical position')
        start = time.time()

        # Get control points
        CP = ctrlpts
        
        # Initialize 
        J = np.zeros((dimensions, dimensions, nb_pts))
        PPS = np.zeros((dimensions, nb_pts))

        # Evaluate jacobien
        for j in range(dimensions):
            alpha = np.zeros(dimensions, dtype = int); alpha[j] = 1

            B = 1
            for dim in range(dimensions):
                at = alpha[dim] 
                B = sp.kron(data_DB[dim][at], B)

            for i in range(dimensions):
                J[i, j, :] = sp.coo_matrix.dot(B.transpose(), CP[:, i])
            
        # Evaluate position in physical space
        for i in range(dimensions):
            B = 1
            for dim in range(dimensions):
                B = sp.kron(data_DB[dim][0], B)

            PPS[i, :] = sp.coo_matrix.dot(B.transpose(), CP[:, i])
        
        stop = time.time()
        print('\tJacobian in : %.5f s' %(stop-start))
        
        return J, PPS

    def eval_thermal_coefficient(self, nb_pts, J_pts, Kprop, Cprop): 
        " Computes coefficients at points P in parametric space "

        print('Getting conductivity and capacity coefficients')
        start = time.time()

        Kcoef = np.zeros(np.shape(J_pts))
        Ccoef = np.zeros(nb_pts)
        detJ = np.zeros(nb_pts)

        for _ in range(nb_pts): 
            # Find determinant of Jacobien 
            det_J = np.linalg.det(J_pts[:, :, _])

            # Find inverse of Jacobien 
            inv_J = np.linalg.inv(J_pts[:, :, _])

            # Find coefficient of conductivity matrix
            Kcoef_temp = inv_J @ inv_J.T @ Kprop * det_J
            Kcoef[:, :, _] = Kcoef_temp

            # Find coefficient of capacity matrix or heat vector
            Ccoef_temp = Cprop * det_J
            Ccoef[_] = Ccoef_temp

            # Find determinant of Jacobien
            detJ[_] = det_J

        stop = time.time()
        print('\tConductivity and capacity coefficients in : %.5f s' %(stop-start))

        return Kcoef, Ccoef, detJ

    def eval_source_coefficient(self, fun, qp, det): 
        " Computes coefficients at points P in parametric space "

        print('Getting source coefficients')
        start = time.time()
        # Get source coefficient
        source_coef = [fun(qp[:, i]) * det[i] 
                    for i in range(len(det))]
        source_coef = np.array(source_coef)
        stop = time.time()
        print('\tSource coefficients in : %.5f s' %(stop-start))

        return source_coef

    def eval_bodyforce_coefficient(self, fun, qp_PS, detJ): 
        " Computes coefficients at points P in parametric space "

        # Set material properties
        E = self._youngModule
        nu = self._poissonCoef

        if self._dim == 2: 
            # coef = (1 + nu)*(1 - 2*nu)/E #!!!!!!!!!!!!!!!!!!!!!! To modify in plane stress 
            coef = (1 - nu*nu)/E
        elif self._dim == 3: 
            coef = (1 + nu)*(1 - 2*nu)/E

        print('Getting body force coefficients')   
        start = time.time()
        # Get source coefficient
        bodyforce_coef = [fun(qp_PS[:, _]) * detJ[_] * coef for _ in range(len(detJ))]
        stop = time.time()
        print('\tBody force coefficients in : %.5f s' %(stop-start))

        return np.asarray(bodyforce_coef)

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

    def interpolate_field(self, u_ctrlpts= None):
        "Interpolates the input field. In all cases, it returns jacobien."

        # =====================
        # Get Basis
        # =====================
        # Set shape
        shape = [1, 1, 1]
        for _ in range(self._dim): shape[_] = self._sample_size
        shape = tuple(shape)

        # Define knots
        knots = np.linspace(0, 1, self._sample_size)

        # Set basis and indices
        DB, ind = [], []
        for _ in range(self._dim):  
            B0, B1, indi, indj = eval_basis_fortran(self._degree[_][0], 
                                                self._nb_el[_][0], knots)
            DB.append([B0, B1])
            ind.append([indi, indj])

        # ==============================
        # Get position and determinant
        # ==============================
        shape_matrices, indices, data, ctrlpts = [], [], [], []
        for dim in range(self._dim):
            shape_matrices.append(self._sample_size)
            indices.append(ind[dim][0])
            indices.append(ind[dim][1])
            data.append(DB[dim][0])
            data.append(DB[dim][1])
            ctrlpts.append(self._ctrlpts[:, dim])
        inputs = [*shape_matrices, *ctrlpts, *indices, *data]

        if self._dim == 2: jacobien_PS, qp_PS, detJ = assembly.jacobien_physicalposition_2d(*inputs)    
        elif self._dim == 3: jacobien_PS, qp_PS, detJ = assembly.jacobien_physicalposition_3d(*inputs)

        # ==============================
        # Get interpolation
        # ==============================
        if u_ctrlpts is not None:
            shape_matrices, indices, data = [], [], []
            for dim in range(self._dim):
                shape_matrices.append(self._sample_size)
                indices.append(ind[dim][0])
                indices.append(ind[dim][1])
                data.append(DB[dim][0])
            inputs = [*shape_matrices, u_ctrlpts, *indices, *data]

            if self._dim == 2: u_interp = assembly.interpolation_2d(*inputs)    
            elif self._dim == 3: u_interp = assembly.interpolation_3d(*inputs)
        else: 
            u_interp = None

        return jacobien_PS, qp_PS, detJ, u_interp
    
    def export_results(self, u_ctrlpts= None, filename= None, folder=None): 
        " Returns solution using geometry basis "

        if folder == None: 
            import os
            full_path = os.path.realpath(__file__)
            dirname = os.path.dirname
            folder = dirname(dirname(full_path)) + '/results/'
            if not os.path.isdir(folder): os.mkdir(folder)

        # Warnings
        if self._geometry_type == None: raise Warning('Try another method to export results')
        if self._isMechanical == True: raise Warning('Not coded, try another method')

        if u_ctrlpts is None:
            UctrlptsExist = False
            u_ctrlpts_new = None
        elif isinstance(u_ctrlpts, np.ndarray): 
            UctrlptsExist = True
            u_ctrlpts_new = u_ctrlpts
        else: 
            raise Warning('Solution must be ndarray type')

        # Set shape
        shape = [1, 1, 1]
        for _ in range(self._dim): shape[_] = self._sample_size
        shape = tuple(shape)

        # ==============================
        # Get interpolation
        # ==============================
        # Interpolate 
        _, qp_PS, detJ, u_interp = self.interpolate_field(u_ctrlpts_new)

        # Find statistics
        mean_detJ = statistics.mean(detJ)
        detJ /= mean_detJ
        variance_detJ = statistics.variance(detJ)
        print("The variance of jacobien normalized is: %.5f" %(variance_detJ))

        # ==============================
        # Export results
        # ==============================
        X1 = np.zeros(shape)
        X2 = np.zeros(shape)
        X3 = np.zeros(shape)
        U = np.zeros(shape)
        DET = np.zeros(shape)

        for k in range(shape[2]):
            for j in range(shape[1]):
                for i in range(shape[0]):
                    pos = i + j * self._sample_size + k * self._sample_size**2
                    X1[i,j,k] = qp_PS[0, pos]
                    X2[i,j,k] = qp_PS[1, pos]
                    if self._dim == 3: X3[i,j,k] = qp_PS[2, pos]

                    if UctrlptsExist: U[i,j,k] = u_interp[pos]
                    DET[i,j,k] = detJ[pos]

        # Export geometry
        if filename == None: 
            try: name = self._name 
            except: name = "VTKResults"
        else: name = filename 
        name =  folder + name
        gridToVTK(name, X1, X2, X3, pointData= {"U1" : U, 
                                                "detJ" : DET,})
        return

