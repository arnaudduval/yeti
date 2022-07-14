"""
.. This module is an adaptation of Matrix free method  
.. for Isogeometric analysis (MF-IGA)
.. Joaquin Cornejo 
"""

# Python libraries
import numpy as np
import time

# My libraries
from .base_functions import erase_rows_csr, iga_find_basis_weights_fortran
from .create_model import thermoMechaModel
from iga_wq_mf import assembly, solver
    
class fortran_mf_iga(thermoMechaModel):
    
    def __init__(self, modelIGA: None, isThermal= True, thermalblockedboundaries= None, **properties):
        super().__init__(modelIGA, isThermal= isThermal, thermalblockedboundaries= thermalblockedboundaries, **properties)

        # Set basis and weights
        self.eval_positions_basis_weights()

        # Get jacobian and physical position 
        self.eval_jacobien_physicalPosition()

        # Initialize thermal properties
        self._conductivity_coef = None
        self._capacity_coef = None

        return

    def eval_positions_basis_weights(self):
        " Computes basis and weights "

        print('Evaluating basis and weights')
        start = time.time()

        # Set number of non-zeros values of I
        self._nnz_I = []

        # Set basis and weights 
        self._qp_cgg, self._DB, self._DW, self._indices = [], [], [], []

        for dim in range(self._dim):  
            nnz_I, qp_pos, W, B, indi, indj  = iga_find_basis_weights_fortran(self._degree[dim], self._knotvector[dim])
            
            self._nnz_I.append(nnz_I)
            self._qp_cgg.append(qp_pos)
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
        inputs = [*self._nb_qp_cgg, *self._indices, *self._DB, self._ctrlpts]

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
    # ISOGEOMETRIC ANALYSIS
    # ===========================

    def eval_source_coefficient(self, fun): 
        " Computes source coefficients "
        source_coef = super().eval_source_coefficient(fun, self._qp_PS, self._detJ)
        return source_coef

    def eval_capacity_matrix(self, indi=None, indj=None): 
        " Computes capacity matrix "

        if indi is None: indi = np.arange(self._nb_ctrlpts_total, dtype=int)
        if indj is None: indj = np.arange(self._nb_ctrlpts_total, dtype=int)
        
        # Get inputs
        self._verify_capacity()
        inputs = [self._capacity_coef, *self._indices, *self._DB, *self._DW, *self._nnz_I]

        start = time.time()
        if self._dim == 2: val, indi, indj = assembly.iga_get_capacity_2d(*inputs)
        if self._dim == 3: val, indi, indj = assembly.iga_get_capacity_3d(*inputs)
                
        # Convert results in csr sparse matrix
        C = super().array2csr_matrix(val, indi, indj).tocsc()[indi, :][:, indj]

        stop = time.time()
        print('Capacity matrix assembled in : %.5f s' %(stop-start))
        
        return  C

    def eval_conductivity_matrix(self, indi=None, indj=None): 
        " Computes conductivity matrix "

        if indi is None: indi = np.arange(self._nb_ctrlpts_total, dtype=int)
        if indj is None: indj = np.arange(self._nb_ctrlpts_total, dtype=int)

        # Get inputs
        self._verify_conductivity()
        inputs = [self._conductivity_coef, *self._indices, *self._DB, *self._DW, *self._nnz_I]

        start = time.time()
        if self._dim == 2: val, indi, indj = assembly.iga_get_conductivity_2d(*inputs)
        if self._dim == 3: val, indi, indj = assembly.iga_get_conductivity_3d(*inputs)
                
        # Convert results in csr sparse matrix
        K = super().array2csr_matrix(val, indi, indj).tocsc()[indi, :][:, indj]

        stop = time.time()
        print('Conductivity matrix assembled in : %5f s' %(stop-start))
        
        return K

    def eval_Ku(self, u): 
        " Computes K u "

        # Get inputs
        self._verify_conductivity()
        inputs = self.get_input4MatrixFree(u)

        if self._dim == 2: raise Warning('Until now not done')
        if self._dim == 3: result = solver.mf_iga_get_ku_3d_csr(*inputs)

        return result
    
    def eval_source_vector(self, fun, indi= None, indj=None, Td=None): 
        " Computes source vector "

        if indi is None: indi = np.arange(self._nb_ctrlpts_total, dtype=int)

        # Get source coefficients
        self._source_coef = self.eval_source_coefficient(fun)
        inputs = [self._source_coef, *self._indices, *self._DB, *self._DW]

        start = time.time()
        if self._dim == 2: F = assembly.iga_get_source_2d(*inputs)
        if self._dim == 3: F = assembly.iga_get_source_3d(*inputs)
        stop = time.time()
        print('Source vector assembled in : %.5f s' %(stop-start))

        return F[indi]

    # =============================
    # MATRIX FREE SOLUTION
    # =============================

    def get_input4MatrixFree(self, table= None):
        " Returns necessary inputs to compute the product between a matrix and a vector"

        # Initialize
        indices, data = [], []
        if table is None: table = np.asarray([[0, 0], [0, 0], [0, 0]])
        for dim in range(self._dim):
            # Select data
            if np.array_equal(table[dim, :], [0, 0]): rows2erase = []
            if np.array_equal(table[dim, :], [0, 1]): rows2erase = [-1]
            if np.array_equal(table[dim, :], [1, 0]): rows2erase = [0]
            if np.array_equal(table[dim, :], [1, 1]): rows2erase = [0, -1]
            indi_t, indj_t, data_t = erase_rows_csr(rows2erase, *self._indices[dim], 
                                    [self._DB[dim], self._DW[dim]])
            
            # Extract data and append to list
            [dB, dW] = data_t
            indices.append(indi_t); indices.append(indj_t) 
            data.append(dB); data.append(dW)

        inputs = [*indices, *data]

        return inputs

    def MFsolver(self, u, nbIterations, epsilon, method, directsol): 

        # Get inputs
        self._verify_conductivity()
        if self._thermalDirichlet is None: raise Warning('Ill conditionned. It needs Dirichlet conditions')
        inputs_tmp = self.get_input4MatrixFree(self._thermalDirichlet)
        inputs = [self._conductivity_coef, *inputs_tmp, u, nbIterations, epsilon, method, 
                self._conductivity, self._Jqp, directsol]

        if self._dim == 2: raise Warning('Until now not done')
        if self._dim == 3: 
            sol, residue, error = solver.iga_mf_cg_3d(*inputs, self._conductivity, self._Jqp, directsol)

        return sol, residue, error
        
    def interpolate_ControlPoints(self, fun, nbIter=100, eps=1e-14):
        
        # Get temperature coeficients 
        coef_F = fun(self._qp_PS)  * self._detJ

        # Define inputs 
        inputs = [coef_F *self._indices, *self._DB, *self._DW]

        # Calculate capacity matrix and temperature vector
        if self._dim == 2: raise Warning('Until now not done')
        if self._dim == 3: F = assembly.iga_get_source_3d(*inputs)

        # Solve linear system with fortran
        inputs = [self._detJ, *self._indices, *self._DB, *self._DW, F, nbIter, eps]
        start = time.time()
        u_interp, relres = solver.wq_mf_interp_3d(*inputs)
        stop = time.time()
        res_end = relres[np.nonzero(relres)][-1]
        print('Interpolation in: %.3e s with relative residue %.3e' %(stop-start, res_end))

        return u_interp