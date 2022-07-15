"""
.. This module is an adaptation of Matrix free - Weighted quadrature method
.. for Isogeometric analysis (MF-WQ-IGA)
.. Joaquin Cornejo 
"""

# Python libraries
import time, numpy as np

# My libraries
from .base_functions import erase_rows_csr, wq_find_basis_weights_fortran
from .create_model import thermoMechaModel
from iga_wq_mf import assembly, solver

class fortran_mf_wq(thermoMechaModel):

    def __init__(self, modelIGA: None, material=None, Dirichlet=None):

        super().__init__(modelIGA, material= material, Dirichlet= Dirichlet)

        # Set basis and weights
        self.eval_basis_weigths()
        
        # Get jacobian and physical position 
        self.eval_jacobien_physicalPosition()

        # Initialize thermal properties
        self._conductivity_coef = None
        self._capacity_coef = None

        return

    def eval_basis_weigths(self): 
        " Computes basis and weights "

        print('Evaluating basis and weights')
        start = time.time()
        # Set number of non-zeros values of I
        self._nnz_I = []

        # Set basis and weights 
        self._qp_wq, self._DB, self._DW, self._indices = [], [], [], []

        for dim in range(self._dim):  
            nnz_I, qp_pos, B, W, indi, indj = wq_find_basis_weights_fortran(self._degree[dim], 
                                                                            self._knotvector[dim])
            self._nnz_I.append(nnz_I)
            self._qp_wq.append(qp_pos)
            self._DB.append(B)
            self._DW.append(W)
            self._indices.append(indi)
            self._indices.append(indj)

        stop = time.time()
        print('\tBasis and weights in : %.5f s' %(stop-start))

        return

    def eval_jacobien_physicalPosition(self): 
        " Computes jacobien and physical position "

        print('Evaluating jacobien and physical position')
        start = time.time()
        # Get inputs
        inputs = [*self._nb_qp_wq, *self._indices, *self._DB, self._ctrlpts]
        
        if self._dim == 2:
            self._Jqp, self._qp_PS, self._detJ = assembly.jacobien_physicalposition_2d(*inputs)
        if self._dim == 3:
            self._Jqp, self._qp_PS, self._detJ = assembly.jacobien_physicalposition_3d(*inputs)
        stop = time.time()
        print('\t Time jacobien: %.5f s' %(stop-start))

        return

    def _verify_conductivity(self): 
        " Verifies if conductivity exits"

        if self._conductivity is None: raise Warning('Conductivity not defined')
        if self._conductivity_coef is None:
            print('Getting conductivity coefficients')
            start = time.time()
            coef, info = assembly.eval_conductivity_coefficient(self._Jqp, self._conductivity)
            if info == 0: raise Warning('Something happen computing coefficients')
            else: self._conductivity_coef = coef
            stop = time.time()
            print('\tConductivity coefficients in : %.5f s' %(stop-start))
        return

    def _verify_capacity(self): 
        " Verifies if capacity exits"

        if self._capacity is None: raise Warning('Capacity not defined')
        if self._capacity_coef is None:
            print('Getting capacity coefficients')
            start = time.time()
            coef, info = assembly.eval_capacity_coefficient(self._Jqp, self._capacity)
            if info == 0: raise Warning('Something happen computing coefficients')
            else: self._capacity_coef = coef
            stop = time.time()
            print('\tCapacity coefficients in : %.5f s' %(stop-start))
        return

    # ===========================
    # WEIGHTED QUADRATURE
    # ===========================

    def eval_source_coefficient(self, fun): 
        " Computes source coefficients "
        source_coef = super().eval_source_coefficient(fun, self._qp_PS, self._detJ)
        return source_coef

    def eval_capacity_matrix(self, indi= None, indj= None): 
        " Computes capacity matrix "

        if indi is None: indi = np.arange(self._nb_ctrlpts_total, dtype=int)
        if indj is None: indj = np.arange(self._nb_ctrlpts_total, dtype=int)
        
        # Get inputs
        self._verify_capacity()
        inputs = [self._capacity_coef, *self._nb_qp_wq, *self._indices, *self._DB, *self._DW, *self._nnz_I]

        start = time.time()
        if self._dim == 2: val_, indi_, indj_ = assembly.wq_get_capacity_2d(*inputs)
        if self._dim == 3: val_, indi_, indj_ = assembly.wq_get_capacity_3d(*inputs)

        # Convert results in csr sparse matrix
        C = super().array2csr_matrix(val_, indi_, indj_).tocsc()[indi, :][:, indj]

        stop = time.time()
        print('Capacity matrix assembled in : %.5f s' %(stop-start))
        
        return  C

    def eval_conductivity_matrix(self, indi= None, indj= None): 
        " Computes conductivity matrix "

        if indi is None: indi = np.arange(self._nb_ctrlpts_total, dtype=int)
        if indj is None: indj = np.arange(self._nb_ctrlpts_total, dtype=int)
        
        # Get inputs
        self._verify_conductivity()
        inputs = [self._conductivity_coef, *self._nb_qp_wq, *self._indices, *self._DB, *self._DW, *self._nnz_I]

        start = time.time()
        if self._dim == 2: val_, indi_, indj_ = assembly.wq_get_conductivity_2d(*inputs)
        if self._dim == 3: val_, indi_, indj_ = assembly.wq_get_conductivity_3d(*inputs)
                
        # Convert results in csr sparse matrix
        K = super().array2csr_matrix(val_, indi_, indj_).tocsc()[indi, :][:, indj]

        stop = time.time()
        print('Conductivity matrix assembled in : %5f s' %(stop-start))
        
        return K

    def eval_Ku(self, u, table= None): 
        " Computes K u "

        # Get inputs
        self._verify_conductivity()
        inputs = self.get_input4MatrixFree(table)

        if self._dim == 2: raise Warning('Until now not done')
        if self._dim == 3: result = solver.mf_wq_get_ku_3d_csr(self._conductivity_coef, *inputs, u)

        return result

    def eval_Cu(self, u, table= None): 
        " Computes K u "

        # Get inputs
        self._verify_capacity()
        inputs = self.get_input4MatrixFree(table)

        if self._dim == 2: raise Warning('Until now not done')
        if self._dim == 3: result = solver.mf_wq_get_cu_3d_csr(self._capacity_coef, *inputs, u)

        return result

    def eval_source_vector(self, fun, indi= None): 
        " Computes source vector Fn - Knd Td "

        if indi is None: indi = np.arange(self._nb_ctrlpts_total, dtype=int)

        # Get source coefficients
        self._source_coef = self.eval_source_coefficient(fun)
        inputs = [self._source_coef, *self._nb_qp_wq, *self._indices, *self._DW]

        start = time.time()
        if self._dim == 2: F = assembly.wq_get_source_2d(*inputs)[indi]
        if self._dim == 3: F = assembly.wq_get_source_3d(*inputs)[indi]
        stop = time.time()
        print('Source vector assembled in : %.5f s' %(stop-start))

        return F

    def eval_diag_K(self): 
        " Computes the diagonal of conductivity matrix "

        # Get inputs
        self._verify_conductivity()
        inputs = [self._conductivity_coef, *self._nb_qp_wq, *self._indices, *self._DB, *self._DW]
        
        start = time.time()
        if self._dim == 2: raise Warning('Until now not done')
        if self._dim == 3: kdiag = assembly.wq_find_conductivity_diagonal_3d(*inputs)

        stop = time.time()
        print('Conductivity matrix assembled in : %5f s' %(stop-start))
        
        return kdiag

    # =============================
    # MATRIX FREE SOLUTION
    # =============================   
     
    def get_input4MatrixFree(self, table= None):
        " Returns necessary inputs to compute the product between a matrix and a vector"

        # Initialize
        indices, data_B, data_W = [], [], []
        if table is None: table = np.asarray([[0, 0], [0, 0], [0, 0]])
        for dim in range(self._dim):
            # Select data
            if np.array_equal(table[dim, :], [0, 0]): rows2erase = []
            if np.array_equal(table[dim, :], [0, 1]): rows2erase = [-1]
            if np.array_equal(table[dim, :], [1, 0]): rows2erase = [0]
            if np.array_equal(table[dim, :], [1, 1]): rows2erase = [0, -1]
            indi_t, indj_t, data_t = erase_rows_csr(rows2erase, 
                                    self._indices[2*dim], self._indices[2*dim+1],  
                                    [self._DB[dim], self._DW[dim]])
            
            # Extract data and append to list
            [dB, dW] = data_t
            indices.append(indi_t); indices.append(indj_t) 
            data_B.append(dB); data_W.append(dW)

        inputs = [*self._nb_qp_wq, *indices, *data_B, *data_W]

        return inputs

    def MFsolver(self, u, nbIterations, epsilon, method, directsol): 

        # Get inputs 
        self._verify_conductivity()
        if self._thermalDirichlet is None: raise Warning('Ill conditionned. It needs Dirichlet conditions')
        inputs_tmp = self.get_input4MatrixFree(self._thermalDirichlet)
        inputs = [self._conductivity_coef, *inputs_tmp, u, nbIterations, epsilon, method, 
                self._conductivity, self._Jqp, directsol]

        if self._dim == 2: raise Warning('Until now not done')
        if self._dim == 3: sol, residue, error = solver.wq_mf_bicgstab_3d(*inputs)

        return sol, residue, error

    def interpolate_ControlPoints(self, fun, nbIter=100, eps=1e-14):
        
        # Get temperature coeficients 
        coef_F = fun(self._qp_PS) * self._detJ

        # Get inputs
        inputs = [coef_F, *self._nb_qp_wq, *self._indices, *self._DW]

        # Calculate temperature vector
        if self._dim == 2: raise Warning('Until now not done')
        if self._dim == 3: F = assembly.wq_get_source_3d(*inputs)

        # Solve linear system with fortran
        inputs = [self._detJ, *self._nb_qp_wq, *self._indices, *self._DB, *self._DW, F, nbIter, eps]
        start = time.time()
        u_interp, relres = solver.wq_mf_interp_3d(*inputs)
        stop = time.time()
        res_end = relres[np.nonzero(relres)][-1]
        print('Interpolation in: %.3e s with relative residue %.3e' %(stop-start, res_end))

        return u_interp
        