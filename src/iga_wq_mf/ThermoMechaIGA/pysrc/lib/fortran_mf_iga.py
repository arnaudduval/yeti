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
    
    def __init__(self, modelIGA: None, isThermal= True, thermalblockedboundaries= None):
        super().__init__(modelIGA, isThermal= isThermal, thermalblockedboundaries= thermalblockedboundaries)

        # Set basis and weights
        self.eval_positions_basis_weights()

        # Get jacobian and physical position 
        if self._isThermal or self._isMechanical:
            self.eval_jacobien_physicalPosition()

        if self._isThermal:
            # Get conductivity and capacity coefficients
            print('Getting conductivity and capacity coefficients')
            start = time.time()
            self._conductivity_coef, self._capacity_coef = \
                assembly.eval_thermal_coefficient(self._Jqp, self._conductivity, self._capacity)
            stop = time.time()
            print('\tConductivity and capacity coefficients in : %.5f s' %(stop-start))

        return

    # ===========================
    # INPUTS 
    # ===========================

    def get_input4jacobien(self): 
        " Returns necessary inputs to compute jacobien "

        # Initialize
        shape_matrices, indices, data, ctrlpts = [], [], [], []

        for dim in range(self._dim):
            shape_matrices.append(self._nb_qp_cgg[dim][0])
            indices.append(self._indices[dim][0])
            indices.append(self._indices[dim][1])
            data.append(self._DB[dim][0])
            data.append(self._DB[dim][1])
            ctrlpts.append(self._ctrlpts[:, dim])
        
        inputs = [*shape_matrices, *ctrlpts, *indices, *data]

        return inputs

    # Assemble matrices

    def get_input4Assembly_capacity(self):
        " Returns necessary inputs to compute capacity matrix "

        # Initialize
        indices, data, size_I = [], [], []

        for dim in range(self._dim):
            indices.append(self._indices[dim][0])
            indices.append(self._indices[dim][1])
            data.append(self._DB[dim][0])
            data.append(self._DW[dim])
            size_I.append(self._nnz_I_dim[dim])

        inputs = [self._capacity_coef, *indices, *data, *size_I]
        
        return inputs

    def get_input4Assembly_conductivity(self):
        " Returns necessary inputs to compute conducivity matrix "

        # Initialize
        indices, data, size_I = [], [], []

        for dim in range(self._dim):
            indices.append(self._indices[dim][0])
            indices.append(self._indices[dim][1])
            data.append(self._DB[dim][0])
            data.append(self._DB[dim][1])
            data.append(self._DW[dim])
            size_I.append(self._nnz_I_dim[dim])

        inputs = [self._conductivity_coef, *indices, *data, *size_I]
        
        return inputs

    # Assemble vectors

    def get_input4Assembly_source(self):
        " Returns necessary inputs to compute source vector "

        # Initialize
        indices, data = [], []

        for dim in range(self._dim):
            indices.append(self._indices[dim][0])
            indices.append(self._indices[dim][1])
            data.append(self._DB[dim][0])
            data.append(self._DW[dim])

        inputs = [self._source_coef, *indices, *data]
        
        return inputs

    # Matrix-Free
    def get_input4MatrixFree_Ku(self, u):
        " Returns necessary inputs to compute K u "

        # Initialize
        indices, data = [], []

        for dim in range(self._dim):
            indices.append(self._indices[dim][0])
            indices.append(self._indices[dim][1])
            data.append(self._DB[dim][0])
            data.append(self._DB[dim][1])
            data.append(self._DW[dim])

        inputs = [self._conductivity_coef, *indices, *data, u]

        return inputs

    def get_input4ConjugateGradient(self, bi, nbIterations, epsilon, method):
        " Returns necessary inputs to compute K u "

        # Initialize
        indices, data = [], []
        table = self._thermalblockedboundaries

        for dim in range(self._dim):

            # Select data
            if np.array_equal(table[dim, :], [0, 0]): rows2erase = []
            if np.array_equal(table[dim, :], [0, 1]): rows2erase = [-1]
            if np.array_equal(table[dim, :], [1, 0]): rows2erase = [0]
            if np.array_equal(table[dim, :], [1, 1]): rows2erase = [0, -1]
            indi_t, indj_t, data_t = erase_rows_csr(rows2erase, *self._indices[dim], [*self._DB[dim]])
            
            # Extract data and append to list
            [dB0, dB1] = data_t
            indices.append(indi_t); indices.append(indj_t) 
            data.append(dB0); data.append(dB1); data.append(self._DW[dim])

        inputs = [self._conductivity_coef, *indices, 
                    *data, bi, nbIterations, epsilon, method]
        return inputs

    # ===========================
    # ISOGEOMETRIC ANALYSIS
    # ===========================

    def eval_positions_basis_weights(self):
        " Computes basis and weights "

        print('Evaluating basis and weights')
        start = time.time()

        # Set number of non-zeros values of I = B . B.transpose 
        self._nnz_I_dim = []

        # Set quadrature points 
        self._qp_wq_dim = []

        # Set basis 
        self._DB = []

        # Set weights 
        self._DW = []

        # Set indices 
        self._indices = []

        for dim in range(self._dim):  
            nnz_I, qp_pos, qp_weights, \
            B0, B1, indi, indj  = iga_find_basis_weights_fortran(self._degree[dim][0], self._nb_el[dim][0])
            
            self._nnz_I_dim.append(nnz_I)
            self._qp_wq_dim.append(qp_pos)
            self._DB.append([B0, B1])
            self._DW.append(qp_weights)
            self._indices.append([indi, indj])

        stop = time.time()
        print('\tBasis and weights in : %.5f s' %(stop-start))

        return

    def eval_jacobien_physicalPosition(self): 
        " Computes jacobien and physical position "

        print('Evaluating jacobien and physical position')
        start = time.time()

        # Get inputs
        inputs = self.get_input4jacobien()

        if self._dim < 2 and self._dim > 3:
            raise Warning('Until now not done')

        if self._dim == 2:
            self._Jqp, self._qp_PS, self._detJ = assembly.jacobien_physicalposition_2d(*inputs)

                
        if self._dim == 3:
            self._Jqp, self._qp_PS, self._detJ = assembly.jacobien_physicalposition_3d(*inputs)
            
        stop = time.time()
        print('\t Time jacobien: %.5f s' %(stop-start))

        return

    def eval_source_coefficient(self, fun): 
        " Computes source coefficients "
        source_coef = super().eval_source_coefficient(fun, self._qp_PS, self._detJ)
        return source_coef

    # Assemble matrices

    def eval_capacity_matrix(self, indi=None, indj=None): 
        " Computes capacity matrix "

        if indi is None: 
            indi = np.arange(self._nb_ctrlpts_total, dtype=int)
        if indj is None:
            indj = np.arange(self._nb_ctrlpts_total, dtype=int)
        
        # Get inputs
        inputs = self.get_input4Assembly_capacity()

        start = time.time()
        if self._dim < 2 and self._dim > 3:
            raise Warning('Until now not done')

        if self._dim == 2:
            val_C, indi_C, indj_C = assembly.iga_get_capacity_2d(*inputs)
                
        if self._dim == 3:
            val_C, indi_C, indj_C = assembly.iga_get_capacity_3d(*inputs)
                
        # Convert results in coo sparse matrix
        C_coo = super().array2csr_matrix(self._nb_ctrlpts_total, self._nb_ctrlpts_total,  
                                            val_C, indi_C, indj_C).tocsc()[indi, :][:, indj]

        stop = time.time()
        print('Capacity matrix assembled in : %.5f s' %(stop-start))
        
        return  C_coo

    def eval_conductivity_matrix(self, indi=None, indj=None): 
        " Computes conductivity matrix "

        if indi is None: 
            indi = np.arange(self._nb_ctrlpts_total, dtype=int)
        if indj is None:
            indj = np.arange(self._nb_ctrlpts_total, dtype=int)

        # Get inputs
        inputs = self.get_input4Assembly_conductivity()

        start = time.time()
        if self._dim < 2 and self._dim > 3:
            raise Warning('Until now not done')

        if self._dim == 2:
            val_K, indi_K, indj_K = assembly.iga_get_conductivity_2d(*inputs)
                
        if self._dim == 3:
            val_K, indi_K, indj_K = assembly.iga_get_conductivity_3d(*inputs)
                
        # Convert results in coo sparse matrix
        K_coo = super().array2csr_matrix(self._nb_ctrlpts_total, self._nb_ctrlpts_total,  
                                            val_K, indi_K, indj_K).tocsc()[indi, :][:, indj]

        stop = time.time()
        print('Conductivity matrix assembled in : %5f s' %(stop-start))
        
        return K_coo

    def eval_Ku(self, u): 
        " Computes K u "

        # Get inputs
        inputs = self.get_input4MatrixFree_Ku(u)

        if self._dim < 2 and self._dim > 3:
            raise Warning('Until now not done')

        if self._dim == 2:
            raise Warning('Until now not done')

        if self._dim == 3:
            result = solver.mf_iga_get_ku_3d_csr(*inputs)

        return result
    
    def eval_source_vector(self, fun, indi= None, indj=None, Td=None): 
        " Computes source vector "

        if indi is None: 
            indi = np.arange(self._nb_ctrlpts_total, dtype=int)

        # Get source coefficients
        self._source_coef = self.eval_source_coefficient(fun)
        inputs = self.get_input4Assembly_source()   

        start = time.time()
        if self._dim < 2 and self._dim > 3:
            raise Warning('Until now not done')

        if self._dim == 2:
            F = assembly.iga_get_source_2d(*inputs)

        if self._dim == 3:
            F = assembly.iga_get_source_3d(*inputs)
            
        stop = time.time()
        print('Source vector assembled in : %.5f s' %(stop-start))

        # Evaluate K T*, where T* = [0, Td]
        if (indj is not None) and (Td is not None):
            if len(Td) != len(indj):
                raise Warning('Different dimensions')

            # Initialize Ttilde
            Ttilde = np.zeros(self._nb_ctrlpts_total)
            Ttilde[indj] = Td
            
            # Eval K@Ttilde
            if self._dim == 2:
                raise Warning('Try another method')
                
            if self._dim == 3:
                KTtilde = self.eval_Ku(Ttilde)

            # Recalculate F
            F -= KTtilde 

        return F[indi]

    # =============================
    # MATRIX FREE SOLUTION
    # =============================

    def MFsolver(self, bi, nbIterations, epsilon, method, directsol, isCG): 

        # Get inputs
        inputs = self.get_input4ConjugateGradient(bi, nbIterations, epsilon, method)
        
        if self._dim < 2 and self._dim > 3:
            raise Warning('Until now not done')

        if self._dim == 2:
            raise Warning('Until now not done')

        if self._dim == 3:
            if isCG :
                sol, residue, error = solver.iga_mf_cg_3d(*inputs, self._Jqp, directsol)
            else: 
                raise Warning('It can only be Conjugate Gradient')

        return sol, residue, error
        
    def interpolate_ControlPoints(self, fun, nbIter=100, eps=1e-14):
        
        # Get temperature coeficients 
        coef_F = [fun(self._qp_PS[:, _])*self._detJ[_] 
                    for _ in range(self._nb_qp_cgg_total)]

        # Define inputs for C and F
        indices, data_F, data_interp = [], [], []

        for dim in range(self._dim):
            indices.append(self._indices[dim][0])
            indices.append(self._indices[dim][1])
            data_F.append(self._DB[dim][0])
            data_F.append(self._DW[dim])
            data_interp.append(self._DB[dim][0])
            data_interp.append(self._DB[dim][1])
            data_interp.append(self._DW[dim])

        inputs_F = [coef_F, *indices, *data_F]

        # Calculate capacity matrix and temperature vector
        if self._dim < 2 and self._dim > 3:
            raise Warning('Until now not done')
        if self._dim == 2:
            raise Warning('Until now not done')
        if self._dim == 3:
            F = assembly.iga_get_source_3d(*inputs_F)

        # Solve linear system with fortran
        start = time.time()
        inputs_interp = [self._detJ, *indices, *data_interp, F, nbIter, eps]
        Tf, relres = solver.iga_mf_interp_3d(*inputs_interp)
        lastres = relres[np.nonzero(relres)][-1]
        Tdir = Tf[self._thermal_dod]
        stop = time.time()
        print('Interpolation in: %.3e s with relative residue %.3e' %(stop-start, lastres))
        
        return Tf, Tdir