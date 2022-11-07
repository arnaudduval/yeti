"""
.. This module is an adaptation of Matrix free method  
.. for Isogeometric analysis (MF-IGA)
.. Joaquin Cornejo 
"""

from .__init__ import *

# My libraries
from .base_functions import erase_rows_csr, iga_find_basis_weights_fortran
from .create_model import thermoMechaModel
    
class fortran_mf_iga(thermoMechaModel):
    
    def __init__(self, modelIGA: None, material={}, Dirichlet={}, Neumann={}):
        super().__init__(modelIGA, material=material, Dirichlet=Dirichlet, Neumann=Neumann)

        self._nb_qp, self._nb_qp_total = np.ones(self._dim, dtype=int), None
        self.eval_basis_weights()
        self.eval_jacobien_physicalPosition()

        return

    def eval_basis_weights(self):
        " Computes basis and weights "

        print('Evaluating basis and weights')
        start = time.process_time()

        self._nnz_I, self._qp_dim, self._DB, self._DW, self._indices = [], [], [], [], []
        for dim in range(self._dim):  
            nnz_I, qp_position, \
            weights, basis, indi, indj = iga_find_basis_weights_fortran(self._degree[dim], self._knotvector[dim])
            self._nb_qp[dim] = len(qp_position)
            self._nnz_I.append(nnz_I); self._qp_dim.append(qp_position)
            self._DB.append(basis); self._DW.append(weights)
            self._indices.append(indi); self._indices.append(indj)

        # Update number of quadrature points
        self._nb_qp_total = np.prod(self._nb_qp)

        stop = time.process_time()
        print('\tBasis and weights in : %.5f s' %(stop-start))

        return

    def eval_jacobien_physicalPosition(self): 
        " Computes jacobien and physical position "

        print('Evaluating jacobien and physical position')
        start = time.process_time()

        inputs = [*self._nb_qp, *self._indices, *self._DB, self._ctrlpts]
        if self._dim == 2:
            self._Jqp, self._detJ, self._invJ = assembly.eval_jacobien_2d(*inputs)
            self._qp_PS = assembly.interpolate_fieldphy_2d(*inputs)
        if self._dim == 3:
            self._Jqp, self._detJ, self._invJ = assembly.eval_jacobien_3d(*inputs)
            self._qp_PS = assembly.interpolate_fieldphy_3d(*inputs)

        stop = time.process_time()
        print('\t Time jacobien: %.5f s' %(stop-start))

        return

    # ----------------------
    # ISOGEOMETRIC ANALYSIS
    # ----------------------

    def eval_source_coefficient(self, fun): 
        " Computes source coefficients "
        coefs = super().eval_source_coefficient(fun, self._qp_PS, self._detJ)
        return coefs

    def eval_capacity_matrix(self, indi=None, indj=None): 
        " Computes capacity matrix "

        if indi is None: indi = np.arange(self._nb_ctrlpts_total, dtype=int)
        if indj is None: indj = np.arange(self._nb_ctrlpts_total, dtype=int)
        
        super()._verify_thermal()
        coefs = super().eval_capacity_coefficient(self._detJ, self._capacity)
        inputs = [coefs, *self._indices, *self._DB, *self._DW, *self._nnz_I]

        start = time.process_time()
        if self._dim == 2: val_, indi_, indj_ = assembly.iga_get_capacity_2d(*inputs)
        if self._dim == 3: val_, indi_, indj_ = assembly.iga_get_capacity_3d(*inputs)
        matrix = super().array2csr_matrix(val_, indi_, indj_).tocsc()[indi, :][:, indj]
        stop = time.process_time()
        
        print('Capacity matrix assembled in : %.5f s' %(stop-start))

        return  matrix

    def eval_conductivity_matrix(self, indi=None, indj=None): 
        " Computes conductivity matrix "

        if indi is None: indi = np.arange(self._nb_ctrlpts_total, dtype=int)
        if indj is None: indj = np.arange(self._nb_ctrlpts_total, dtype=int)

        super()._verify_thermal()
        coefs = super().eval_conductivity_coefficient(self._invJ, self._detJ, self._conductivity)
        inputs = [coefs, *self._indices, *self._DB, *self._DW, *self._nnz_I]

        start = time.process_time()
        if self._dim == 2: val_, indi_, indj_ = assembly.iga_get_conductivity_2d(*inputs)
        if self._dim == 3: val_, indi_, indj_ = assembly.iga_get_conductivity_3d(*inputs)
        matrix = super().array2csr_matrix(val_, indi_, indj_).tocsc()[indi, :][:, indj]
        stop = time.process_time()

        print('Conductivity matrix assembled in : %5f s' %(stop-start))
        
        return matrix

    def eval_Ku(self, u, table=None): 
        " Computes K u where K is conductivity matrix "

        super()._verify_thermal()
        coefs = super().eval_conductivity_coefficient(self._invJ, self._detJ, self._conductivity)
        inputs = self.get_input4MatrixFree(table=table)

        start = time.process_time()
        if self._dim == 2: raise Warning('Until now not done')
        if self._dim == 3: result = solver.mf_iga_get_ku_3d_py(coefs, *inputs, u)  
        stop = time.process_time()
        timeCPU = stop - start

        return result, timeCPU
    
    def eval_source_vector(self, fun, indi= None): 
        " Computes source vector "

        if indi is None: indi = np.arange(self._nb_ctrlpts_total, dtype=int)
        coefs = self.eval_source_coefficient(fun)
        inputs = [coefs, *self._indices, *self._DB, *self._DW]

        start = time.process_time()
        if self._dim == 2: vector = assembly.iga_get_source_2d(*inputs)[indi]
        if self._dim == 3: vector = assembly.iga_get_source_3d(*inputs)[indi]
        stop = time.process_time()

        print('Source vector assembled in : %.5f s' %(stop-start))

        return vector

    # ----------------------
    # MATRIX FREE SOLUTION
    # ----------------------

    def get_input4MatrixFree(self, table=None):
        " Returns necessary inputs to compute the product between a matrix and a vector"
        
        if table is None: table = np.asarray([[0, 0], [0, 0], [0, 0]])

        indices, data = [], []
        for dim in range(self._dim):
            # Select data
            if np.array_equal(table[dim, :], [0, 0]): rows2erase = []
            if np.array_equal(table[dim, :], [0, 1]): rows2erase = [-1]
            if np.array_equal(table[dim, :], [1, 0]): rows2erase = [0]
            if np.array_equal(table[dim, :], [1, 1]): rows2erase = [0, -1]
            indi_t, indj_t, data_t = erase_rows_csr(rows2erase, 
                                    self._indices[2*dim], self._indices[2*dim+1],  
                                    [self._DB[dim]])
            # Extract data and append to list
            indices.append(indi_t); indices.append(indj_t) 
            data.append(*data_t)

        inputs = [*indices, *data, *self._DW]

        return inputs

    def MFsteadyHeat(self, b, nbIterPCG=100, threshold=1e-12, methodPCG='FDC'): 
        " Solves steady heat problems using directly substitution method "

        if self._thermalDirichlet is None: raise Warning('Ill conditionned. It needs Dirichlet conditions')

        super()._verify_thermal()
        coefs = super().eval_conductivity_coefficient(self._invJ, self._detJ, self._conductivity)
        inputs_tmp = self.get_input4MatrixFree(table=self._thermalDirichlet)
        inputs = [coefs, *inputs_tmp, b, nbIterPCG, threshold, methodPCG]

        if self._dim == 2: raise Warning('Until now not done')
        if self._dim == 3: sol, residue = solver.mf_iga_steady_heat_3d(*inputs)

        return sol, residue
        
    def interpolate_ControlPoints(self, funfield=None, datafield=None, nbIterPCG=100, threshold=1e-14):
        " Interpolation from parametric space to physical space "
        
        coefs = None
        if datafield is not None: coefs = datafield * self._detJ
        if funfield is not None: coefs = funfield(self._qp_PS) * self._detJ
        if coefs is None: raise Warning('Missing data')

        # Calculate vector
        inputs = [coefs, *self._indices, *self._DB, *self._DW]
        if self._dim == 2: raise Warning('Until now not done')
        if self._dim == 3: F = assembly.iga_get_source_3d(*inputs)

        # Solve linear system with fortran
        inputs = [self._detJ, *self._indices, *self._DB, *self._DW, F, nbIterPCG, threshold]
        start = time.process_time()
        u_interp, relres = solver.mf_iga_interpolate_cp_3d(*inputs)
        stop = time.process_time()
        res_end = relres[np.nonzero(relres)][-1]
        print('Interpolation in: %.3e s with relative residue %.3e' %(stop-start, res_end))

        return u_interp