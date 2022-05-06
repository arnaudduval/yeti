"""
.. This module is an adaptation of Matrix free - Weighted quadrature method
.. for Isogeometric analysis (MF-WQ-IGA)
.. Joaquin Cornejo 
"""

# Python libraries
import numpy as np
import time

# My libraries
from .create_model import thermoMechaModel
from iga_wq_mf import basis_weights, assembly, solver

def wq_get_basis_weights(degree, nb_el): 
    " Computes basis and weights "

    # Set size guessed of data arrays
    nnz_B, size_qp = basis_weights.wq_get_size_data(degree, nb_el)

    # Get basis and weights from fortran
    qp_pos, B0, B1, \
    W00, W01, W10, W11, \
    ind, nnz_I = basis_weights.wq_get_data(degree, nb_el, nnz_B, size_qp)

    return nnz_I, qp_pos, B0, B1, W00, W01, W10, W11, ind

def wq_get_basis_weights_csr(degree, nb_el): 
    " Computes basis and weights "

    # Set size guessed of data arrays
    nnz_B, size_qp = basis_weights.wq_get_size_data(degree, nb_el)

    # Get basis and weights from fortran
    qp_pos, B0, B1, \
    W00, W01, W10, W11, \
    indi, indj, nnz_I = basis_weights.wq_get_data_csr(degree, nb_el, nnz_B, size_qp)

    return nnz_I, qp_pos, B0, B1, W00, W01, W10, W11, indi, indj

class fortran_mf_wq(thermoMechaModel):

    # ===========================
    # INITIALIZE 
    # ===========================

    def __init__(self, modelIGA: None, isThermal= True, isMechanical= False, 
                thermalblockedboundaries= None, 
                mechablockedboundaries= None):
        super().__init__(modelIGA= modelIGA, isThermal= isThermal, isMechanical= isMechanical,
                        thermalblockedboundaries= thermalblockedboundaries,
                        mechablockedboundaries= mechablockedboundaries)

        # Set basis and weights
        self.eval_basis_weigths()
        
        # Get jacobian and physical position 
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

    def get_input4basisweights(self): 
        " Returns necessary inputs to compute basis and weights "

        # Initialize
        inputs = []

        for dim in range(self._dim):
            inputs_t = [self._degree[dim][0], 
                        self._nb_el[dim][0]]
            inputs.append(inputs_t)

        return inputs

    def get_input4jacobien(self): 
        " Returns necessary inputs to compute jacobien "

        # Set shape of matrices B0 and B1
        shape_matrices = []

        # Set indexes 
        indexes = []

        # Set data 
        data = []

        # Set control points
        ctrlpts = []

        for dim in range(self._dim):
            shape_matrices.append(self._nb_ctrlpts[dim][0])
            shape_matrices.append(self._nb_qp_wq[dim][0])
            indexes.append(self._indexes[dim][:, 0])
            indexes.append(self._indexes[dim][:, 1])
            data.append(self._DB[dim][0])
            data.append(self._DB[dim][1])
            ctrlpts.append(self._ctrlpts[:, dim])
        
        inputs = [*shape_matrices, *ctrlpts, *indexes, *data]

        return inputs

    def get_input4Assembly_capacity(self): 
        " Returns necessary inputs to compute capacity matrix "

        # Set shape of matrices B0 and B1
        shape_matrices = []

        # Set indexes 
        indexes = []

        # Set data 
        data = []

        # Set size of I
        size_I = []

        for dim in range(self._dim):
            shape_matrices.append(self._nb_ctrlpts[dim][0])
            shape_matrices.append(self._nb_qp_wq[dim][0])
            indexes.append(self._indexes[dim][:, 0])
            indexes.append(self._indexes[dim][:, 1])
            data.append(self._DB[dim][0])
            data.append(self._DW[dim][0][0])
            size_I.append(self._nnz_I_dim[dim])

        inputs = [self._capacity_coef, *shape_matrices, 
                *indexes, *data, *size_I]

        return inputs

    def get_input4Assembly_conductivity(self):

        # Set shape of matrices B0 and B1
        shape_matrices = []

        # Set indexes 
        indexes = []

        # Set data 
        data = []

        # Set size of I
        size_I = []

        for dim in range(self._dim):
            shape_matrices.append(self._nb_ctrlpts[dim][0])
            shape_matrices.append(self._nb_qp_wq[dim][0])
            indexes.append(self._indexes[dim][:, 0])
            indexes.append(self._indexes[dim][:, 1])
            data.append(self._DB[dim][0])
            data.append(self._DB[dim][1])
            data.append(self._DW[dim][0][0])
            data.append(self._DW[dim][0][1])
            data.append(self._DW[dim][1][0])
            data.append(self._DW[dim][1][1])
            size_I.append(self._nnz_I_dim[dim])

        inputs = [self._conductivity_coef, *shape_matrices, 
                    *indexes, *data, *size_I]
        
        return inputs

    def get_input4Assembly_source(self):
        " Returns necessary inputs to compute source vector "

        # Set shape of matrices B0 and B1
        shape_matrices = []

        # Set indexes 
        indexes = []

        # Set data 
        data = []

        for dim in range(self._dim):
            shape_matrices.append(self._nb_ctrlpts[dim][0])
            shape_matrices.append(self._nb_qp_wq[dim][0])
            indexes.append(self._indexes[dim][:, 0])
            indexes.append(self._indexes[dim][:, 1])
            data.append(self._DW[dim][0][0])

        inputs = [self._source_coef, *shape_matrices, *indexes, *data]
        
        return inputs

    def get_input4ConjugateGradient(self, bi, dof, nbIterations, epsilon):
        " Returns necessary inputs to compute K u "

        # Set shape of matrices B0 and B1
        shape_matrices = []

        # Set indexes 
        indexes = []

        # Set data 
        data = []

        # !!! Remember: Indexes in fortran starts at 1
        dof = np.asarray(dof) + 1

        for dim in range(self._dim):
            shape_matrices.append(self._nb_ctrlpts[dim][0])
            shape_matrices.append(self._nb_qp_wq[dim][0])
            indexes.append(self._indexes[dim][:, 0])
            indexes.append(self._indexes[dim][:, 1])
            data.append(self._DB[dim][0])
            data.append(self._DB[dim][1])
            data.append(self._DW[dim][0][0])
            data.append(self._DW[dim][0][1])
            data.append(self._DW[dim][1][0])
            data.append(self._DW[dim][1][1])

        inputs = [self._conductivity_coef, self._thermalblockedboundaries, *shape_matrices, *indexes, 
                *data, bi, dof, nbIterations, epsilon]

        return inputs

    def get_input4Assembly_stiffness(self):

        # Set shape of matrices B0 and B1
        shape_matrices = []

        # Set indexes 
        indexes = []

        # Set data 
        data = []

        # Set size of I
        size_I = []

        for dim in range(self._dim):
            shape_matrices.append(self._nb_ctrlpts[dim][0])
            shape_matrices.append(self._nb_qp_wq[dim][0])
            indexes.append(self._indexes[dim][:, 0])
            indexes.append(self._indexes[dim][:, 1])
            data.append(self._DB[dim][0])
            data.append(self._DB[dim][1])
            data.append(self._DW[dim][0][0])
            data.append(self._DW[dim][0][1])
            data.append(self._DW[dim][1][0])
            data.append(self._DW[dim][1][1])
            size_I.append(self._nnz_I_dim[dim])

        inputs = [self._stiff_coef, *shape_matrices,
                *indexes, *data, *size_I]

        return inputs

    def get_input4Assembly_force(self):
        " Returns necessary inputs to compute source vector "

        # Set shape of matrices B0 and B1
        shape_matrices = []

        # Set indexes 
        indexes = []

        # Set data 
        data = []

        for dim in range(self._dim):
            shape_matrices.append(self._nb_ctrlpts[dim][0])
            shape_matrices.append(self._nb_qp_wq[dim][0])
            indexes.append(self._indexes[dim][:, 0])
            indexes.append(self._indexes[dim][:, 1])
            data.append(self._DW[dim][0][0])

        inputs = [self._bodyforce_coef, *shape_matrices, *indexes, *data]
        
        return inputs

    def get_input4Assembly_thermalstiffness(self):

        # Set shape of matrices B0 and B1
        shape_matrices = []

        # Set indexes 
        indexes = []

        # Set data 
        data = []

        # Set size of I
        size_I = []

        for dim in range(self._dim):
            shape_matrices.append(self._nb_ctrlpts[dim][0])
            shape_matrices.append(self._nb_qp_wq[dim][0])
            indexes.append(self._indexes[dim][:, 0])
            indexes.append(self._indexes[dim][:, 1])
            data.append(self._DB[dim][0])
            data.append(self._DW[dim][0][0])
            data.append(self._DW[dim][1][0])
            size_I.append(self._nnz_I_dim[dim])

        inputs = [self._thermalstiff_coef, *shape_matrices, 
                    *indexes, *data, *size_I]

        return inputs

    def get_input4Find_conductivity_diagonal(self, indi_csr):
        # Set shape of matrices B0 and B1
        shape_matrices = []

        # Set indexes 
        indexes = []

        # Set data 
        data = []

        for dim in range(self._dim):
            indexes.append(indi_csr[dim])
            indexes.append(self._indexes[dim][:, 1])
            data.append(self._DB[dim][0])
            data.append(self._DB[dim][1])
            data.append(self._DW[dim][0][0])
            data.append(self._DW[dim][0][1])
            data.append(self._DW[dim][1][0])
            data.append(self._DW[dim][1][1])
            shape_matrices.append(self._nb_qp_wq[dim][0])

        inputs = [self._conductivity_coef, *shape_matrices, *indexes, *data]

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

        # Set indexes 
        self._indexes = []

        # Get inputs
        inputs = self.get_input4basisweights()

        for dim in range(self._dim):  
            nnz_I, qp_pos, B0, B1, W00, W01, \
            W10, W11, indexes = wq_get_basis_weights(*inputs[dim])
            
            self._nnz_I_dim.append(nnz_I)
            self._qp_wq_dim.append(qp_pos)
            self._DB.append([B0, B1])
            self._DW.append([[W00, W01], [W10, W11]])
            self._indexes.append(indexes)

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

    def eval_bodyforce_coefficient(self, fun): 
        " Computes source coefficients "
        bodyforce = super().eval_BodyForce_coefficient(fun, self._dim, self._qp_PS, self._detJ)
        return bodyforce

    def eval_capacity_matrix(self): 
        " Computes capacity matrix "
        
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
                                            val_C, indi_C, indj_C)

        stop = time.time()
        print('Capacity matrix assembled in : %.5f s' %(stop-start))
        
        return  C_coo

    def eval_conductivity_matrix(self): 
        " Computes conductivity matrix "

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
                                            val_K, indi_K, indj_K)

        stop = time.time()
        print('Conductivity matrix assembled in : %5f s' %(stop-start))
        
        return K_coo

    def eval_source_vector(self, fun): 
        " Computes source vector "

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

        return F

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
        S_coo = super().array2coo_matrix( self._dim*self._nb_ctrlpts_total, 
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

    # =============================
    # MATRIX FREE SOLUTION
    # =============================
    def mf_find_K_diagonal(self): 

        def get_basis_weigths_csr(): 
            indi_csr = []

            for dim in range(self._dim):  
                _, _, _, _, _, _, _, _, indi, _ = wq_get_basis_weights_csr(self._degree[dim][0], 
                                                self._nb_el[dim][0])
                indi_csr.append(indi)

            return indi_csr

        # Get inputs
        indi_csr = get_basis_weigths_csr()
        inputs = self.get_input4Find_conductivity_diagonal(indi_csr)
        
        if self._dim < 2 and self._dim > 3:
            raise Warning('Until now not done')

        if self._dim == 2:
            raise Warning('Until now not done')

        if self._dim == 3:
            Kdiag = solver.wq_find_conductivity_diagonal_3d(*inputs)

        return Kdiag
    
    def mf_conj_grad(self, bi, dof, nbIterations, epsilon, method, directsol, isCG): 

        # Get inputs
        inputs = self.get_input4ConjugateGradient(bi, dof, nbIterations, epsilon)
        
        if self._dim < 2 and self._dim > 3:
            raise Warning('Until now not done')

        if self._dim == 2:
            raise Warning('Until now not done')

        if self._dim == 3:
            if isCG:
                sol, residue, error = solver.wq_mf_cg_3d(*inputs, method, self._Jqp, directsol)
            else:
                sol, residue, error = solver.wq_mf_bicgstab_3d(*inputs, method, self._Jqp, directsol)

        return sol, residue, error
        