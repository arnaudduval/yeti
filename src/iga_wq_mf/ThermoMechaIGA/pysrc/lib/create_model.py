"""
.. This module helps to read information from geometries (YETI format and Geomdl)
.. Joaquin Cornejo
"""

# Python libraries
import statistics
import numpy as np
from scipy import sparse as sp
import time
import matplotlib.pyplot as plt
from pyevtk.hl import gridToVTK 

# Yeti libraries
from preprocessing.igaparametrization import IGAparametrization

# My libraries
from .base_functions import eval_basis_fortran
from .create_geomdl import geomdlModel
from iga_wq_mf import assembly

def write_text_file(filename, method_list, inputs): 

    # Define inputs
    time_assembly = inputs["TimeAssembly"]
    time_direct = inputs["TimeDirect"]
    memory_direct = inputs["MemDirect"]
    time_iter = inputs["TimeIter"]
    time_noiter = inputs["TimeNoIter"]
    residue = inputs["Res"]
    error = inputs["Error"]
    memory_iter = inputs["MemIter"]
    memory_noiter = inputs["MemNoIter"]

    with open(filename, 'w') as outputfile:
        outputfile.write('** RESULTS **\n')
        outputfile.write('* DIRECT SOLVER\n')
        outputfile.write("{:E}".format(time_assembly) + "\t"
                        +"{:E}".format(time_direct) + "\t" 
                        +"{:E}".format(memory_direct) +"\n"
        )

        for i, pcgmethod in enumerate(method_list):
            outputfile.write('* ITERATIVE SOLVER : ' + pcgmethod + '\n')
            outputfile.write("{:E}".format(time_noiter[i]) + '\t' 
                            +"{:E}".format(memory_noiter[i])+ '\t' 
                            +"{:E}".format(time_iter[i]) + '\t' 
                            +"{:E}".format(memory_iter[i]) +'\n'
            )
            outputfile.writelines(["{:E}".format(res) + "\t"+ "{:E}".format(err)+"\n"
                            for res, err in zip(residue[i], error[i])]) 
    return

def read_text_file(filename):
    # Read file 
    residue_list = []
    error_list = []
    time_iter = []
    time_noiter = []
    memory_iter = []
    memory_noiter = []

    with open(filename, 'r') as inputfile:
        lines = inputfile.readlines()

        text_position = []
        for _ in range(len(lines)):
            if lines[_][0] == '*': 
                text_position.append(_)
        text_position.append(len(lines)+1)

        # We know that our first numerical line is direct method: 
        lines_data_split = lines[2].split('\t')
        time_assembly = float(lines_data_split[0])
        time_direct = float(lines_data_split[1])
        memory_direct = float(lines_data_split[2])
        del text_position[0:2]

        # For the different type of solvers
        for _ in range(len(text_position)-1):
            ind = [text_position[_]+1, text_position[_+1]]
            lines_data = lines[ind[0]:ind[1]]
            lines_data_split = lines_data[0].split('\t')
            time_noiter.append(float(lines_data_split[0]))
            memory_noiter.append(float(lines_data_split[1]))
            time_iter.append(float(lines_data_split[2]))
            memory_iter.append(float(lines_data_split[3]))

            residue_tmp = []
            error_tmp = []
            for i  in range(1, len(lines_data)):
                lines_data_split = lines_data[i].split('\t')
                residue_tmp.append(float(lines_data_split[0]))
                error_tmp.append(float(lines_data_split[1]))
                
            residue_list.append(residue_tmp)
            error_list.append(error_tmp)

        output = {"TimeAssembly": time_assembly, "TimeDirect": time_direct, "MemDirect": memory_direct, 
                "TimeNoIter":time_noiter, "TimeIter": time_iter, "Res": residue_list, 
                "Error": error_list, "MemNoIter": memory_noiter, "MemIter": memory_iter}

    return output

def plot_iterative_solver(filename, inputs, method_list, extension ='.png'):
    # Define color
    CB_color_cycle = ['#377eb8', '#ff7f00', '#4daf4a',
                    '#f781bf', '#a65628', '#984ea3',
                    '#999999', '#e41a1c', '#dede00']

    # Define inputs
    residue = inputs["Res"]
    error = inputs["Error"]

    # Get new names
    new_method_list = []
    for i, pcgmethod in enumerate(method_list):
        if pcgmethod == "WP": new_method_list.append('Without preconditioner')
        elif pcgmethod == "C": new_method_list.append('Fast Diagonalisation (FD)')
        elif pcgmethod == "TDS": new_method_list.append('FD + tensor decomp. + scaling') 
        elif pcgmethod == "JM": new_method_list.append('FD + jacobien mean')  
        elif pcgmethod == "TD": new_method_list.append('FD + tensor decomp.') 
        elif pcgmethod == "JMS": new_method_list.append('FD + jacobien mean + scaling')

    # Select important values
    tol = 1.e-17
    
    # Set figure parameters
    fig, [ax1, ax2] = plt.subplots(nrows=1, ncols=2, figsize=(12,6))

    for i, pcgmethod in enumerate(new_method_list):
        error_method = np.asarray(error[i])
        residue_method = np.asarray(residue[i])

        error_method = error_method[residue_method>tol]
        residue_method = residue_method[residue_method>tol]
        
        ax1.plot(np.arange(len(error_method)), error_method*100, label=pcgmethod, color=CB_color_cycle[i])
        ax2.plot(np.arange(len(residue_method)), residue_method*100, label=pcgmethod, color=CB_color_cycle[i])

    # Set properties
    ax1.set_yscale("log")
    ax1.set_xlabel('Number of iterations', fontsize=16)
    ax1.set_ylabel('Relative error (%)', fontsize=16)
    ax1.tick_params(axis='x', labelsize=16)
    ax1.tick_params(axis='y', labelsize=16)
    ax1.legend(loc='best',fontsize=12)
    ax1.grid()

    ax2.set_yscale("log")
    ax2.set_xlabel('Number of iterations', fontsize=16)
    ax2.set_ylabel('Relative residue (%)', fontsize=16)
    ax2.tick_params(axis='x', labelsize=16)
    ax2.tick_params(axis='y', labelsize=16)
    ax2.legend(loc='best',fontsize=12)
    ax2.grid()

    # Save figure
    plt.tight_layout()
    plt.savefig(filename + extension)
    plt.close(fig)

    return

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
        except: self._sample_size = 51

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
            self._thermalblockedboundaries = np.asarray(thermalblockedboundaries)
            self._mechablockedboundaries = np.asarray(mechablockedboundaries)

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

    def eval_jacobien_pps(self, dimensions, ctrlpts, data_DB, nb_pts): 
        """ Computes Jacobien matrix and find the position in physical space 
        of P (in pts list) in parametric space
        """
        print('Evaluating jacobien and physical position')
        start = time.time()

        # Get control points
        CP = ctrlpts
        
        # Initialize 
        J = np.zeros((dimensions, dimensions, nb_pts))
        PPS = np.zeros((1, dimensions, nb_pts))

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

            PPS[0, i, :] = sp.coo_matrix.dot(B.transpose(), CP[:, i])
        
        stop = time.time()
        print('\tJacobian in : %.5f s' %(stop-start))
        
        return J, PPS

    def eval_K_C_coefficient(self, nb_pts, J_pts, Kprop, Cprop): 
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

    def eval_F_coefficient(self, fun, qp, det): 
        " Computes coefficients at points P in parametric space "

        print('Getting source coefficients')
        start = time.time()
        # Get source coefficient
        source_coef = [fun(qp[:, :, i][0]) * det[i] 
                    for i in range(len(det))]
        source_coef = np.array(source_coef)
        stop = time.time()
        print('\tSource coefficients in : %.5f s' %(stop-start))

        return source_coef

    def eval_BodyForce_coefficient(self, fun, qp_PS, detJ): 
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
        bodyforce_coef = [fun(qp_PS[:, :, _][0]) * detJ[_] * coef for _ in range(len(detJ))]
        stop = time.time()
        print('\tBody force coefficients in : %.5f s' %(stop-start))

        return np.asarray(bodyforce_coef)

    def conjugate_gradient_python(self, fun_Au, bi, dof, nbIterations=100, epsilon=1e-10):   
        " Evaluate K u at choosen equations "

        # ------------------
        # Conjugate Gradient algorithm
        # ------------------
        x = np.zeros(len(bi))
        r = bi
        p = r
        rsold = np.dot(r, r)
        RelRes = np.zeros(nbIterations+1)
        RelRes[0] = 1.0

        for k in range(nbIterations):
            Ap = fun_Au(p, dof)
            alpha = rsold/np.dot(p, Ap)
            x = x + alpha*p
            r = r - alpha*Ap
            RelRes[k+1] = (np.linalg.norm(r, np.inf)/np.linalg.norm(bi, np.inf))

            if RelRes[k+1]<epsilon:
                break
            rsnew = np.dot(r, r)
            p = r + rsnew/rsold * p
            rsold = rsnew

        return x, RelRes

    def conjugate_gradient_scipy(self, A, b, nbIterations=100, epsilon=1e-10, PreCond='ilu', isCG=True):

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

    def reconstruct_solution(self, x, dof):
        " Returns the data reconstructed "

        # Degrees of freedom
        MCRD = 3

        # Number of control points
        nb_ctrlpts_total = self._nb_ctrlpts_total
        
        # Initialize solution
        SOL_vector = np.zeros(nb_ctrlpts_total * MCRD)
        SOL_matrix = np.zeros((nb_ctrlpts_total, MCRD))
        
        # Set solution
        SOL_vector[dof] = x

        # Array to matrix 
        for _ in range(nb_ctrlpts_total): 
            SOL_matrix[_, :] = [SOL_vector[_], SOL_vector[_+nb_ctrlpts_total], 
                                SOL_vector[_+2*nb_ctrlpts_total]]

        return SOL_matrix

    def interpolate_results(self, u_ctrlpts= None):
        # =====================
        # Get Basis
        # =====================
        # Set shape
        shape = [1, 1, 1]
        for _ in range(self._dim):
            shape[_] = self._sample_size
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
            return jacobien_PS, qp_PS, detJ, u_interp
        else: 
            return jacobien_PS, qp_PS, detJ
    
    def export_results(self, u_ctrlpts= None, filename= None): 
        " Returns solution using geometry basis "

        if self._geometry_type == None:
            raise Warning('Try another method to export results')
        
        if self._isMechanical == True:
            raise Warning('Not coded, try another method')

        # For now, only consider scalar field as input !!!!!!!!!!!!!!!!!!!
        if  u_ctrlpts is None:
            UctrlptsExist = False
        elif isinstance(u_ctrlpts, np.ndarray): 
            UctrlptsExist = True
            try: u_ctrlpts_new = u_ctrlpts[:, 0] 
            except: u_ctrlpts_new = u_ctrlpts
        else: 
            raise Warning('Solution must be ndarray type')

        # Set shape
        shape = [1, 1, 1]
        for _ in range(self._dim):
            shape[_] = self._sample_size
        shape = tuple(shape)

        # ==============================
        # Get interpolation
        # ==============================
        if UctrlptsExist:
            jacobien_PS, qp_PS, detJ, u_interp = self.interpolate_results(u_ctrlpts_new)
        else :
            jacobien_PS, qp_PS, detJ = self.interpolate_results()

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
                    X1[i,j,k] = qp_PS[0, 0, pos]
                    X2[i,j,k] = qp_PS[0, 1, pos]
                    if self._dim == 3: X3[i,j,k] = qp_PS[0, 2, pos]

                    if UctrlptsExist: U[i,j,k] = u_interp[pos]
                    DET[i,j,k] = detJ[pos]

        # Export geometry
        if filename == None: 
            try: name = self._name 
            except: name = "Results_VTK"
        else: name = filename 
        gridToVTK(name, X1, X2, X3, pointData= {"Temp" : U, 
                                                "detJ" : DET,})
        return

    def plot_pojection_jacobien(self, normalvector=0):

        # =====================
        # Get Basis
        # =====================
        # Define knots
        knots = np.linspace(0, 1, self._sample_size)

        # Set basis 
        DB = []

        # Set indices 
        ind = []

        for dim in range(self._dim):  
            B0, B1, indi, indj = eval_basis_fortran(self._degree[dim][0], self._nb_el[dim][0], knots)
            DB.append([B0, B1])
            ind.append([indi, indj])

        # =====================
        # Get inputs
        # =====================
        shape_matrices = []
        indices = []
        data = []
        ctrlpts = []

        for dim in range(self._dim):
            shape_matrices.append(self._sample_size)
            indices.append(ind[dim][0])
            indices.append(ind[dim][1])
            data.append(DB[dim][0])
            data.append(DB[dim][1])
            ctrlpts.append(self._ctrlpts[:, dim])
        
        inputs = [self._nb_ctrlpts_total, *shape_matrices, *ctrlpts, *indices, *data]

        # =========================
        # Get position and jacobien
        # =========================
        if self._dim == 2: raise Warning('Until now not done')    
                
        if self._dim == 3:
            _, qp_PS, detJ = assembly.jacobien_physicalposition_3d(*inputs)

        # Find statistics
        mean_detJ = statistics.mean(detJ)

        # Normalize data
        detJ /= mean_detJ
        variance_detJ = statistics.variance(detJ)
        print("The variance of jacobien normalized is: %.5f" %(variance_detJ))

        # ----------------------------------------------------------------------
        # FIRST TRY: det J
        # ----------------------------------------------------------------------
        if self._dim == 2:
            X1 = np.zeros((self._sample_size, self._sample_size))
            X2 = np.zeros((self._sample_size, self._sample_size))
            U = np.zeros((self._sample_size, self._sample_size))
            for j in range(self._sample_size):
                for i in range(self._sample_size):
                    pos = i + j * self._sample_size
                    X1[i, j] = qp_PS[0, 0, pos]
                    X2[i, j] = qp_PS[0, 1, pos]
                    U[i, j] = detJ[pos]
            X1label = 'X (m)'
            X2label = 'Y (m)'

        elif self._dim ==3: 
            X1 = np.zeros((self._sample_size, self._sample_size))
            X2 = np.zeros((self._sample_size, self._sample_size))
            U = np.zeros((self._sample_size, self._sample_size))
            indpos = round((self._sample_size+1)/2)

            if normalvector == 0:
                # Normal following x 
                for k in range(self._sample_size):
                    for j in range(self._sample_size):
                        for i in [0]:
                            pos = i + j * self._sample_size + k * self._sample_size**2
                            X1[j,k] = qp_PS[0, 1, pos]
                            X2[j,k] = qp_PS[0, 2, pos]
                            U[j,k] = detJ[pos]
                X1label = 'Y (m)'
                X2label = 'Z (m)'

            elif normalvector == 1:
                # Normal following y
                for k in range(self._sample_size):
                    for j in [0]:
                        for i in range(self._sample_size):
                            pos = i + j * self._sample_size + k * self._sample_size**2
                            X1[i,k] = qp_PS[0, 0, pos]
                            X2[i,k] = qp_PS[0, 2, pos]
                            U[i,k] = detJ[pos]
                X1label = 'X (m)'
                X2label = 'Z (m)'
            
            elif normalvector == 2:
                # Normal following z
                for k in [indpos]:
                    for j in range(self._sample_size):
                        for i in range(self._sample_size):
                            pos = i + j * self._sample_size + k * self._sample_size**2
                            X1[i,j] = qp_PS[0, 0, pos]
                            X2[i,j] = qp_PS[0, 1, pos]
                            U[i,j] = detJ[pos]
                X1label = 'X (m)'
                X2label = 'Y (m)'
                            
            else:
                raise Warning('That direction is not coded yet. Try another method')

        # Plot the surface
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.pcolormesh(X1, X2, U, cmap='inferno', shading = 'gouraud')
        ax.set_xlabel(X1label, fontsize=14); 
        ax.set_ylabel(X2label, fontsize=14)

       # Set parameters
        ax.axis('equal')
        ax.tick_params(axis='both', which='major', labelsize=16)
        ax.set_xlabel('X', fontsize=16)
        ax.set_ylabel('Y', fontsize=16)
        fig.tight_layout()

        return fig

        