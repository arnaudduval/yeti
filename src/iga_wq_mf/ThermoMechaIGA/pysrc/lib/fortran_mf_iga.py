"""
.. This module is an adaptation of Matrix free method  
.. for Isogeometric analysis (MF-IGA)
.. Joaquin Cornejo 
"""

# Python libraries
import numpy as np
import time

# My libraries
from .create_model import thermoMechaModel
from iga_wq_mf import basis_weights, assembly, solver

def iga_get_basis_weights(degree, nb_el): 
    " Computes basis and weights "

    # Set number of quadrature points
    nb_qp = (degree + 1) * nb_el

    # Set size guessed of data arrays
    nnz_B = (degree + 1) * nb_qp

    # Get basis and weights from fortran
    qp_pos, qp_weights, data_b0, data_b1, \
    data_ind, nnz_I = basis_weights.iga_get_data(degree, nb_el, nnz_B)
 
    # Resize outputs
    B0 = data_b0
    B1 = data_b1
    ind = data_ind

    return nnz_I, qp_pos, qp_weights, B0, B1, ind
    
class fortran_mf_iga(thermoMechaModel):
    
    def __init__(self, modelIGA: None, thermalblockedboundaries= None):
        super().__init__(modelIGA, thermalblockedboundaries= thermalblockedboundaries)
        if self._dim != 3:
            raise Warning('This class is only avaliable for 3D geometries')

        # Set basis and weights
        self.eval_positions_basis_weights()
        
        # Get jacobian and physical position 
        self.eval_jacobien_physicalPosition()

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
        shape_matrices, indexes, data, ctrlpts = [], [], [], []

        for dim in range(self._dim):
            shape_matrices.append(self._nb_ctrlpts[dim][0])
            shape_matrices.append(self._nb_qp_cgg[dim][0])
            indexes.append(self._indexes[dim][:, 0])
            indexes.append(self._indexes[dim][:, 1])
            data.append(self._DB[dim][0])
            data.append(self._DB[dim][1])
            ctrlpts.append(self._ctrlpts[:, dim])
        
        inputs = [*shape_matrices, *ctrlpts, *indexes, *data]

        return inputs

    # Assemble matrices

    def get_input4Assembly_capacity(self):
        " Returns necessary inputs to compute capacity matrix "

        # Initialize
        shape_matrices, indexes, data, size_I = [], [], [], []

        for dim in range(self._dim):
            shape_matrices.append(self._nb_ctrlpts[dim][0])
            indexes.append(self._indexes[dim][:, 0])
            indexes.append(self._indexes[dim][:, 1])
            data.append(self._DB[dim][0])
            data.append(self._DW[dim])
            size_I.append(self._nnz_I_dim[dim])

        inputs = [self._capacity_coef, *shape_matrices, *indexes, *data, *size_I]
        
        return inputs

    def get_input4Assembly_conductivity(self):
        " Returns necessary inputs to compute conducivity matrix "

        # Initialize
        shape_matrices, indexes, data, size_I = [], [], [], []

        for dim in range(self._dim):
            shape_matrices.append(self._nb_ctrlpts[dim][0])
            indexes.append(self._indexes[dim][:, 0])
            indexes.append(self._indexes[dim][:, 1])
            data.append(self._DB[dim][0])
            data.append(self._DB[dim][1])
            data.append(self._DW[dim])
            size_I.append(self._nnz_I_dim[dim])

        inputs = [self._conductivity_coef, *shape_matrices, *indexes, *data, *size_I]
        
        return inputs

    # Assemble vectors

    def get_input4Assembly_source(self):
        " Returns necessary inputs to compute source vector "

        # Initialize
        shape_matrices, indexes, data = [], [], []

        for dim in range(self._dim):
            shape_matrices.append(self._nb_ctrlpts[dim][0])
            indexes.append(self._indexes[dim][:, 0])
            indexes.append(self._indexes[dim][:, 1])
            data.append(self._DB[dim][0])
            data.append(self._DW[dim])

        inputs = [self._source_coef, *shape_matrices, *indexes, *data]
        
        return inputs

    # Matrix-Free

    def get_input4ConjugateGradient(self, bi, dof, nbIterations, epsilon):
        " Returns necessary inputs to compute K u "

        # Initialize
        shape_matrices, indexes, data = [], [], []
        
        # Remember: Indexes in fortran starts at 1
        dof = np.asarray(dof) + 1

        for dim in range(self._dim):
            shape_matrices.append(self._nb_ctrlpts[dim][0])
            indexes.append(self._indexes[dim][:, 0])
            indexes.append(self._indexes[dim][:, 1])
            data.append(self._DB[dim][0])
            data.append(self._DB[dim][1])
            data.append(self._DW[dim])

        inputs = [self._conductivity_coef, self._thermalblockedboundaries, *shape_matrices, *indexes, 
                    *data, bi, dof, nbIterations, epsilon]
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

        # Set indexes 
        self._indexes = []

        for dim in range(self._dim):  
            nnz_I, qp_pos, qp_weights, \
            B0, B1, ind  = iga_get_basis_weights(self._degree[dim][0], self._nb_el[dim][0])
            
            self._nnz_I_dim.append(nnz_I)
            self._qp_wq_dim.append(qp_pos)
            self._DB.append([B0, B1])
            self._DW.append(qp_weights)
            self._indexes.append(ind)

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
        source_coef = super().eval_F_coefficient(fun, self._dim, self._qp_PS, self._detJ)
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
            raise Warning('Until now not done')
                
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
            raise Warning('Until now not done')
                
        if self._dim == 3:
            val_K, indi_K, indj_K = assembly.iga_get_conductivity_3d(*inputs)
                
        # Convert results in coo sparse matrix
        K_coo = super().array2csr_matrix(self._nb_ctrlpts_total, self._nb_ctrlpts_total,  
                                            val_K, indi_K, indj_K).tocsc()[indi, :][:, indj]

        stop = time.time()
        print('Conductivity matrix assembled in : %5f s' %(stop-start))
        
        return K_coo

    def eval_source_vector(self, fun, indi= None): 
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
            raise Warning('Until now not done')   

        if self._dim == 3:
            F = assembly.iga_get_source_3d(*inputs)
            
        stop = time.time()
        print('Source vector assembled in : %.5f s' %(stop-start))

        return F[indi]

    # =============================
    # MATRIX FREE SOLUTION
    # =============================

    def mf_conj_grad(self, bi, dof, nbIterations, epsilon, method, directsol, isCG): 

        # Get inputs
        inputs = self.get_input4ConjugateGradient(bi, dof, nbIterations, epsilon)
        
        if self._dim < 2 and self._dim > 3:
            raise Warning('Until now not done')

        if self._dim == 2:
            raise Warning('Until now not done')

        if self._dim == 3:
            if isCG :
                sol, residue, error = solver.iga_mf_cg_3d(*inputs, method, self._Jqp, directsol)
            else: 
                raise Warning('It can only be Conjugate Gradient')

        return sol, residue, error
        
    def MSE_ControlPoints(self, fun):
        
        # Get temperature coeficients 
        coef_F = [fun(self._dim, self._qp_PS[:, :, _][0])*self._detJ[_] 
                    for _ in range(self._nb_qp_cgg_total)]

        # Define inputs for C and F
        shape_matrices, indexes, data, size_I = [], [], [], []

        for dim in range(self._dim):
            shape_matrices.append(self._nb_ctrlpts[dim][0])
            indexes.append(self._indexes[dim][:, 0])
            indexes.append(self._indexes[dim][:, 1])
            data.append(self._DB[dim][0])
            data.append(self._DW[dim])
            size_I.append(self._nnz_I_dim[dim])

        inputs_C = [self._detJ, *shape_matrices, *indexes, *data, *size_I]
        inputs_F = [coef_F, *shape_matrices, *indexes, *data]

        # Calculate capacity matrix and temperature vector
        if self._dim < 2 and self._dim > 3:
            raise Warning('Until now not done')
        if self._dim == 2:
            raise Warning('Until now not done')
        if self._dim == 3:
            val_C, indi_C, indj_C = assembly.iga_get_capacity_3d(*inputs_C)
            F = assembly.iga_get_source_3d(*inputs_F)
        C = super().array2csr_matrix(self._nb_ctrlpts_total, self._nb_ctrlpts_total,  
                                            val_C, indi_C, indj_C).tocsc()

        # Solve linear system
        T = self.conjugate_gradient_scipy(C, F)
        Tdir = T[self._thermal_dod]

        return T, Tdir