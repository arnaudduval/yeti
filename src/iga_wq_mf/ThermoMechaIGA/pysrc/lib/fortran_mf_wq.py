"""
.. This module is an adaptation of Matrix free - Weighted quadrature method
.. for Isogeometric analysis (MF-WQ-IGA)
.. Joaquin Cornejo 
"""

# Python libraries
import numpy as np
import time

# My libraries
from .base_functions import erase_rows_csr, wq_find_basis_weights_fortran
from .create_model import thermoMechaModel
from iga_wq_mf import assembly, solver

class fortran_mf_wq(thermoMechaModel):

    def __init__(self, modelIGA: None, isThermal= True, isMechanical= False, 
                thermalblockedboundaries= None, mechablockedboundaries= None, **properties):

        super().__init__(modelIGA= modelIGA, isThermal= isThermal, isMechanical= isMechanical,
                        thermalblockedboundaries= thermalblockedboundaries,
                        mechablockedboundaries= mechablockedboundaries, **properties)

        # Set basis and weights
        self.eval_basis_weigths()
        
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

        if self._isMechanical:
            # Get stiffness coefficients
            print('Getting elastic coefficients')
            start = time.time()
            self._stiff_coef = \
                assembly.eval_mech_coefficient(self._Jqp, self._ElasticMatrix)
            stop = time.time()
            print('\tElastic coefficients in : %.5f s' %(stop-start))

        if self._isMechanical and self._isThermal:
            # Get thermal stiffness coefficients
            print('Getting thermo-elastic coefficients')
            start = time.time()
            self._thermalstiff_coef = \
                assembly.eval_thermomech_coefficient(self._Jqp, self._ElasticMatrix, self._thermalexpansion)
            stop = time.time()
            print('\tThermo-elastic coefficients in : %.5f s' %(stop-start))

        return
        
    # ===========================
    # INPUTS 
    # ===========================

    def get_input4jacobien(self): 
        " Returns necessary inputs to compute jacobien matrix "

        # Initialize
        shape_matrices, indices, data, ctrlpts = [], [], [], []

        for dim in range(self._dim):
            shape_matrices.append(self._nb_qp_wq[dim][0])
            indices.append(self._indices[dim][0])
            indices.append(self._indices[dim][1])
            data.append(self._DB[dim][0])
            data.append(self._DB[dim][1])
            ctrlpts.append(self._ctrlpts[:, dim])
        
        inputs = [self._nb_qp_wq_total, *shape_matrices, *ctrlpts, *indices, *data]

        return inputs

    # Assemble matrices

    def get_input4Assembly_capacity(self): 
        " Returns necessary inputs to compute capacity matrix "

        # Initialize
        shape_matrices, indices, data, size_I = [], [], [], []

        for dim in range(self._dim):
            shape_matrices.append(self._nb_qp_wq[dim][0])
            indices.append(self._indices[dim][0])
            indices.append(self._indices[dim][1])
            data.append(self._DB[dim][0])
            data.append(self._DW[dim][0][0])
            size_I.append(self._nnz_I_dim[dim])

        inputs = [self._capacity_coef, *shape_matrices, 
                *indices, *data, *size_I]

        return inputs

    def get_input4Assembly_conductivity(self):
        " Returns necessary inputs to compute conductivity matrix "
        
        # Initialize
        shape_matrices, indices, data, size_I = [], [], [], []

        for dim in range(self._dim):
            shape_matrices.append(self._nb_qp_wq[dim][0])
            indices.append(self._indices[dim][0])
            indices.append(self._indices[dim][1])
            data.append(self._DB[dim][0])
            data.append(self._DB[dim][1])
            data.append(self._DW[dim][0][0])
            data.append(self._DW[dim][0][1])
            data.append(self._DW[dim][1][0])
            data.append(self._DW[dim][1][1])
            size_I.append(self._nnz_I_dim[dim])

        inputs = [self._conductivity_coef, *shape_matrices, 
                    *indices, *data, *size_I]
        
        return inputs

    def get_input4Assembly_stiffness(self):
        " Returns necessary inputs to compute stiffness matrix "

        # Initialize
        shape_matrices, indices, data, size_I = [], [], [], []

        for dim in range(self._dim):
            shape_matrices.append(self._nb_qp_wq[dim][0])
            indices.append(self._indices[dim][0])
            indices.append(self._indices[dim][1])
            data.append(self._DB[dim][0])
            data.append(self._DB[dim][1])
            data.append(self._DW[dim][0][0])
            data.append(self._DW[dim][0][1])
            data.append(self._DW[dim][1][0])
            data.append(self._DW[dim][1][1])
            size_I.append(self._nnz_I_dim[dim])

        inputs = [self._stiff_coef, *shape_matrices,
                *indices, *data, *size_I]

        return inputs

    def get_input4Assembly_thermalstiffness(self):

        " Returns necessary inputs to compute thermal stiffness matrix "

        # Initialize
        shape_matrices, indices, data, size_I = [], [], [], []

        for dim in range(self._dim):
            shape_matrices.append(self._nb_qp_wq[dim][0])
            indices.append(self._indices[dim][0])
            indices.append(self._indices[dim][1])
            data.append(self._DB[dim][0])
            data.append(self._DW[dim][0][0])
            data.append(self._DW[dim][1][0])
            size_I.append(self._nnz_I_dim[dim])

        inputs = [self._thermalstiff_coef, *shape_matrices, 
                    *indices, *data, *size_I]

        return inputs

    # Assemble vectors

    def get_input4Assembly_source(self):
        " Returns necessary inputs to compute source vector "

        # Initialize
        shape_matrices, indices, data = [], [], []

        for dim in range(self._dim):
            shape_matrices.append(self._nb_qp_wq[dim][0])
            indices.append(self._indices[dim][0])
            indices.append(self._indices[dim][1])
            data.append(self._DW[dim][0][0])

        inputs = [self._source_coef, *shape_matrices, *indices, *data]
        
        return inputs

    def get_input4Assembly_force(self):
        " Returns necessary inputs to compute force vector "

        # Initialize
        shape_matrices, indices, data = [], [], []

        for dim in range(self._dim):
            shape_matrices.append(self._nb_qp_wq[dim][0])
            indices.append(self._indices[dim][0])
            indices.append(self._indices[dim][1])
            data.append(self._DW[dim][0][0])

        inputs = [self._bodyforce_coef, *shape_matrices, *indices, *data]
        
        return inputs

    def get_input4diagonal_K(self): 
        " Returns necessary inputs to compute the diagonal of conductivity matrix "
        
        # Initialize
        shape_matrices, indices, data = [], [], [], []

        for dim in range(self._dim):
            shape_matrices.append(self._nb_qp_wq[dim][0])
            indices.append(self._indices[dim][0])
            indices.append(self._indices[dim][1])
            data.append(self._DB[dim][0])
            data.append(self._DB[dim][1])
            data.append(self._DW[dim][0][0])
            data.append(self._DW[dim][0][1])
            data.append(self._DW[dim][1][0])
            data.append(self._DW[dim][1][1])

        inputs = [self._conductivity_coef, *shape_matrices, *indices, *data]
        
        return inputs
    
    # Matrix-Free

    def get_input4FastDiagonalization(self, b): 
        " Returns necessary inputs to solve P r = s, where P is preconditioner of Fast diagoalization "

        # Initialize
        shape_matrices, indices, data = [], [], []
        table = self._thermalblockedboundaries

        for dim in range(self._dim):
            shape_matrices.append(self._nb_qp_wq[dim][0])

            # Select data
            if np.array_equal(table[dim, :], [0, 0]): rows2erase = []
            if np.array_equal(table[dim, :], [0, 1]): rows2erase = [-1]
            if np.array_equal(table[dim, :], [1, 0]): rows2erase = [0]
            if np.array_equal(table[dim, :], [1, 1]): rows2erase = [0, -1]
            indi_t, indj_t, data_t = erase_rows_csr(rows2erase, *self._indices[dim], 
                                    [*self._DB[dim], *self._DW[dim][0], *self._DW[dim][1]])
            
            # Extract data and append to list
            [dB0, dB1, dW00, _, _, dW11] = data_t
            indices.append(indi_t); indices.append(indj_t) 
            data.append(dB0); data.append(dB1)
            data.append(dW00); data.append(dW11)

        inputs = [*shape_matrices, *indices, *data, b]

        return inputs

    def get_input4MatrixFree_Ku(self, u, table=None):
        " Returns necessary inputs to compute K u "

        # Initialize
        shape_matrices, indices, data = [], [], []
        if table is None: table = np.asarray([[0, 0], [0, 0], [0, 0]])

        for dim in range(self._dim):
            shape_matrices.append(self._nb_qp_wq[dim][0])

            # Select data
            if np.array_equal(table[dim, :], [0, 0]): rows2erase = []
            if np.array_equal(table[dim, :], [0, 1]): rows2erase = [-1]
            if np.array_equal(table[dim, :], [1, 0]): rows2erase = [0]
            if np.array_equal(table[dim, :], [1, 1]): rows2erase = [0, -1]
            indi_t, indj_t, data_t = erase_rows_csr(rows2erase, *self._indices[dim], 
                                    [*self._DB[dim], *self._DW[dim][0], *self._DW[dim][1]])
            
            # Extract data and append to list
            [dB0, dB1, dW00, dW01, dW10, dW11] = data_t
            indices.append(indi_t); indices.append(indj_t) 
            data.append(dB0); data.append(dB1)
            data.append(dW00); data.append(dW01); data.append(dW10); data.append(dW11)

        inputs = [self._conductivity_coef, *shape_matrices, *indices, *data, u]

        return inputs

    def get_input4MatrixFree_Cu(self, u, table=None):
        " Returns necessary inputs to compute C u "

        # Initialize
        shape_matrices, indices, data = [], [], []
        if table is None: table = np.asarray([[0, 0], [0, 0], [0, 0]])

        for dim in range(self._dim):
            shape_matrices.append(self._nb_qp_wq[dim][0])

            # Select data
            if np.array_equal(table[dim, :], [0, 0]): rows2erase = []
            if np.array_equal(table[dim, :], [0, 1]): rows2erase = [-1]
            if np.array_equal(table[dim, :], [1, 0]): rows2erase = [0]
            if np.array_equal(table[dim, :], [1, 1]): rows2erase = [0, -1]
            indi_t, indj_t, data_t = erase_rows_csr(rows2erase, *self._indices[dim], 
                                    [*self._DB[dim], *self._DW[dim][0], *self._DW[dim][1]])
            
            # Extract data and append to list
            [dB0, dB1, dW00, dW01, dW10, dW11] = data_t
            indices.append(indi_t); indices.append(indj_t) 
            data.append(dB0); data.append(dW00)

        inputs = [self._capacity_coef, *shape_matrices, *indices, *data, u]

        return inputs

    def get_input4ConjugateGradient(self, bi, nbIterations, epsilon, method):
        " Returns necessary inputs to solve matrix system with matrix free approach "

        # Initialize
        shape_matrices, indices, data = [], [], []
        table = self._thermalblockedboundaries

        for dim in range(self._dim):
            shape_matrices.append(self._nb_qp_wq[dim][0])

            # Select data
            if np.array_equal(table[dim, :], [0, 0]): rows2erase = []
            if np.array_equal(table[dim, :], [0, 1]): rows2erase = [-1]
            if np.array_equal(table[dim, :], [1, 0]): rows2erase = [0]
            if np.array_equal(table[dim, :], [1, 1]): rows2erase = [0, -1]
            indi_t, indj_t, data_t = erase_rows_csr(rows2erase, *self._indices[dim], 
                                    [*self._DB[dim], *self._DW[dim][0], *self._DW[dim][1]])
            
            # Extract data and append to list
            [dB0, dB1, dW00, dW01, dW10, dW11] = data_t
            indices.append(indi_t); indices.append(indj_t) 
            data.append(dB0); data.append(dB1)
            data.append(dW00); data.append(dW01); data.append(dW10); data.append(dW11)

        inputs = [self._conductivity_coef, *shape_matrices, *indices, 
                *data, bi, nbIterations, epsilon, method]

        return inputs

    # ===========================
    # WEIGHTED QUADRATURE
    # ===========================

    def eval_basis_weigths(self): 
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
            nnz_I, qp_pos, B0, B1, W00, W01, \
            W10, W11, indi, indj = wq_find_basis_weights_fortran(self._degree[dim][0], self._nb_el[dim][0])
            
            self._nnz_I_dim.append(nnz_I)
            self._qp_wq_dim.append(qp_pos)
            self._DB.append([B0, B1])
            self._DW.append([[W00, W01], [W10, W11]])
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

    def eval_bodyforce_coefficient(self, fun): 
        " Computes source coefficients "
        bodyforce = super().eval_bodyforce_coefficient(fun, self._qp_PS, self._detJ)
        return bodyforce

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
            val_C, indi_C, indj_C = assembly.wq_get_capacity_2d(*inputs)
                
        if self._dim == 3:
            val_C, indi_C, indj_C = assembly.wq_get_capacity_3d(*inputs)

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
            val_K, indi_K, indj_K = assembly.wq_get_conductivity_2d(*inputs)
                
        if self._dim == 3:
            val_K, indi_K, indj_K = assembly.wq_get_conductivity_3d(*inputs)
                
        # Convert results in coo sparse matrix
        K_coo = super().array2csr_matrix(self._nb_ctrlpts_total, self._nb_ctrlpts_total,  
                                            val_K, indi_K, indj_K).tocsc()[indi, :][:, indj]

        stop = time.time()
        print('Conductivity matrix assembled in : %5f s' %(stop-start))
        
        return K_coo

    def eval_Ku(self, u, table=None): 
        " Computes K u "

        # Get inputs
        inputs = self.get_input4MatrixFree_Ku(u, table)

        if self._dim < 2 and self._dim > 3:
            raise Warning('Until now not done')

        if self._dim == 2:
            raise Warning('Until now not done')

        if self._dim == 3:
            result = solver.mf_wq_get_ku_3d_csr(*inputs)

        return result

    def eval_Cu(self, u, table=None): 
        " Computes K u "

        # Get inputs
        inputs = self.get_input4MatrixFree_Cu(u, table)

        if self._dim < 2 and self._dim > 3:
            raise Warning('Until now not done')

        if self._dim == 2:
            raise Warning('Until now not done')

        if self._dim == 3:
            result = solver.mf_wq_get_cu_3d_csr(*inputs)

        return result

    def eval_source_vector(self, fun, indi= None, indj=None, Td=None): 
        " Computes source vector Fn - Knd Td "

        if indi is None: 
            indi = np.arange(self._nb_ctrlpts_total, dtype=int)

        # Get source coefficients
        self._source_coef = self.eval_source_coefficient(fun)
        inputs = self.get_input4Assembly_source()   

        start = time.time()
        if self._dim < 2 and self._dim > 3:
            raise Warning('Until now not done')

        if self._dim == 2:
            F = assembly.wq_get_source_2d(*inputs)
                
        if self._dim == 3:
            F = assembly.wq_get_source_3d(*inputs)
            
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

    # Thermo-mechanical matrices ! Until not used or well-defined

    def eval_stiffness_matrix(self):
        " Computes stiffness matrix "

        # Get inputs
        inputs = self.get_input4Assembly_stiffness()

        start = time.time()
        if self._dim < 2 and self._dim > 3:
            raise Warning('Until now not done')

        if self._dim == 2:
            val_K, indi_K, indj_K = assembly.wq_get_stiffness_2d(*inputs)

        if self._dim == 3:
            val_K, indi_K, indj_K = assembly.wq_get_stiffness_3d(*inputs)

        # Convert results in coo sparse matrix
        S_coo = super().array2coo_matrix(self._dim*self._nb_ctrlpts_total, 
                                         self._dim*self._nb_ctrlpts_total,  
                                         val_K, indi_K, indj_K)

        stop = time.time()
        print('Stiffness matrix assembled in : %.5f s' %(stop-start))
        
        return S_coo

    def eval_bodyforce_vector(self, fun): 
        " Computes source vector "

        # Get source coefficients
        self._bodyforce_coef = self.eval_bodyforce_coefficient(fun)
        inputs = self.get_input4Assembly_force()   

        start = time.time()
        if self._dim < 2 and self._dim > 3:
            raise Warning('Until now not done')

        if self._dim == 2:
            F = assembly.wq_get_force_2d(*inputs)
                
        if self._dim == 3:
            F = assembly.wq_get_force_3d(*inputs)
            
        stop = time.time()
        print('Body force vector assembled in : %.5f s' %(stop-start))

        return F

    def eval_thermalstiffness_matrix(self):
        " Computes thermal stiffness matrix "

        # Get inputs
        inputs = self.get_input4Assembly_thermalstiffness()

        start = time.time()
        if self._dim < 2 and self._dim > 3:
            raise Warning('Until now not done')

        if self._dim == 2:
            raise Warning('Until now not done')

        if self._dim == 3:
            val_K, indi_K, indj_K = assembly.wq_get_thermalstiffness_3d(*inputs)

        # Convert results in coo sparse matrix
        S_coo = super().array2coo_matrix( self._dim*self._nb_ctrlpts_total, 
                                         self._nb_ctrlpts_total,  
                                         val_K, indi_K, indj_K)

        stop = time.time()
        print('Thermal stiffness matrix assembled in : %.5f s' %(stop-start))
        
        return S_coo

    def eval_diag_K(self): 
        " Computes conductivity matrix "

        # Get inputs
        inputs = self.get_input4diagonal_K()

        start = time.time()
        if self._dim < 2 and self._dim > 3:
            raise Warning('Until now not done')

        if self._dim == 2:
            raise Warning('Until now not done')
                
        if self._dim == 3:
            kdiag = assembly.wq_find_conductivity_diagonal_3d(*inputs)

        stop = time.time()
        print('Conductivity matrix assembled in : %5f s' %(stop-start))
        
        return kdiag

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
            if isCG:
                sol, residue, error = solver.wq_mf_cg_3d(*inputs, self._Jqp, directsol)
            else:
                sol, residue, error = solver.wq_mf_bicgstab_3d(*inputs, self._Jqp, directsol)

        return sol, residue, error

    def interpolate_ControlPoints(self, fun, nbIter=100, eps=1e-14):
        
        # Get temperature coeficients 
        fun_qp = fun(self._qp_PS)  
        coef_F = fun_qp * self._detJ

        # Define inputs for C and F
        shape_matrices, indices, data_interp, data_F = [], [], [], []

        for dim in range(self._dim):
            shape_matrices.append(self._nb_qp_wq[dim][0])
            indices.append(self._indices[dim][0])
            indices.append(self._indices[dim][1])
            data_interp.append(self._DB[dim][0])
            data_interp.append(self._DB[dim][1])
            data_interp.append(self._DW[dim][0][0])
            data_interp.append(self._DW[dim][1][1])
            data_F.append(self._DW[dim][0][0])

        inputs_F = [coef_F, *shape_matrices, *indices, *data_F]

        # Calculate capacity matrix and temperature vector
        if self._dim < 2 and self._dim > 3:
            raise Warning('Until now not done')
        if self._dim == 2:
            raise Warning('Until now not done')
        if self._dim == 3:
            F = assembly.wq_get_source_3d(*inputs_F)

        # Solve linear system with fortran
        start = time.time()
        inputs_interp = [self._detJ, *shape_matrices, *indices, *data_interp, F, nbIter, eps]
        Tf, relres = solver.wq_mf_interp_3d(*inputs_interp)
        lastres = relres[np.nonzero(relres)][-1]
        Tdir = Tf[self._thermal_dod]
        stop = time.time()
        print('Interpolation in: %.3e s with relative residue %.3e' %(stop-start, lastres))

        return Tf, Tdir
        