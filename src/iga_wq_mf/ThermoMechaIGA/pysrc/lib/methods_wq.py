"""
.. WQ methods
.. by Joaquin Cornejo
.. Disclaimer :: This module will hardly ever be used and 
..               and maybe it does not work as expected
..               Fortran functions to calculate basis are more efficient
"""

# Python libraries
import numpy as np
from scipy import sparse as sp
import time

# My libraries
from .base_functions import wq_find_basis_weights_opt, get_matrix_3D, get_indexes_3D
from .create_model import thermoMechaModel

# Yeti libraries
from preprocessing.igaparametrization import IGAparametrization

class WQ(thermoMechaModel): 

    def __init__(self, modelIGA: IGAparametrization):
        super().__init__(modelIGA)

        # Evaluate basis and weights
        self.eval_basis_weights()

        # Evaluate jacobian 
        self._Jqp, self._qp_PS = super().eval_jacobien_pps(self._dim, self._ctrlpts, 
                                                        self._DB, self._nb_qp_wq_total)

        # Evaluate conductivity and capacity coefficients
        self._conductivity_coef, self._capacity_coef, self._detJ = \
            super().eval_K_C_coefficient(self._nb_qp_wq_total, self._Jqp, 
                                        self._conductivity, self._capacity)
        
        return 

    def eval_basis_weights(self): 
        " Computes Basis and weights in WQ approach "
        
        print('Evaluating basis and weights')
        start = time.time()

        # Initalize storage vectors
        self._qp_wq_dim = []
        self._DB = [] # Stocks B-spline functions
        self._DW = [] # Stocks weights

        for dim in range(self._dim): 
            # Find basis and weights 
            qp_pos, B0, B1, W00, W01, W10, W11 = \
            wq_find_basis_weights_opt(self._degree[dim][0], self._knotvector[0][dim], self._r_)  

            # Save quadrature points
            self._qp_wq_dim.append(qp_pos)

            # Save DB
            self._DB.append([B0, B1])

            # Save DW
            self._DW.append([[W00, W01], [W10, W11]])

        stop = time.time()
        print('\tBasis and weights in : %.5f s' %(stop-start))

    def eval_source_coefficient(self, fun): 
        " Computes source coefficients "
        source_coef = super().eval_F_coefficient(fun, self._qp_PS, self._detJ)
        return source_coef
    
    # ===========================
    # ASSEMBLE 
    # ===========================
    def eval_conductivity_matrix(self):
        " Assemble conductivity matrix K "

        start = time.time()
        # Initialize conductivity matrix 
        K = sp.csr_matrix((self._nb_ctrlpts_total, self._nb_ctrlpts_total))
        for j in range(self._dim):
            beta = np.zeros(self._dim, dtype = int); beta[j] = 1
            Bt = 1 
            for dim in range(self._dim): 
                bt = beta[dim]
                Bt = sp.kron(self._DB[dim][bt], Bt)

            for i in range(self._dim):
                Cij = self._conductivity_coef[i, j, :]
                alpha = np.zeros(self._dim, dtype = int); alpha[i] = 1
                
                Wt = 1
                for dim in range(self._dim): 
                    at = alpha[dim]
                    bt = beta[dim]
                    Wt = sp.kron(self._DW[dim][at][bt], Wt)  

                # Evaluates Cij * W in each point
                Wt = sp.csr_matrix.dot(Wt, sp.diags(Cij))

                # Find K = W C B
                K += sp.csr_matrix.dot(Wt.tocsr()[:,:], Bt.tocsr()[:,:].T)
        
        stop = time.time()
        print('Conductivity matrix assembled in : %.5f s' %(stop-start))

        return K

    def eval_capacity_matrix(self):
        " Assemble capacity matrix C "

        start = time.time()
        # Initialize basis and weights
        B = 1; W = 1

        for dim in range(self._dim): 
            # Find basis and weights
            B0 = self._DB[dim][0]
            W00 = self._DW[dim][0][0]

            # Find W and B
            B = sp.kron(B0, B)
            W = sp.kron(W00, W)

        # Assemble C
        W = sp.csr_matrix.dot(W, sp.diags(self._capacity_coef))
        C = sp.csr_matrix.dot(W.tocsr()[:,:], B.tocsr()[:,:].T)

        stop = time.time()
        print('Capacity matrix assembled in : %.5f s' %(stop-start))

        return C

    def eval_source_vector(self, fun):
        " Assemble power density vector F "

        # Get source coefficients
        self._source_coef = self.eval_source_coefficient(fun)       

        start = time.time()
        # Initialize weights
        W = 1; 
        for dim in range(self._dim): 
            # Find weights
            W00 = self._DW[dim][0][0]

            # Find W
            W = sp.kron(W00, W)

        # Assemble F
        F = sp.csr.csr_matrix.dot(W.tocsr()[:,:], self._source_coef)
        stop = time.time()
        print('Source vector assembled in : %.5f s' %(stop-start))

        return F

    def eval_C_TenProd(self): 
        " Computes capacity matrix with tensor product algorithm "
        " !!!!!!!!!!!!!!!!!!!!!!!!!!! TO BE TESTED BEFORE "

        # These algorithms only have been tested with 3D geometries
        if self._dim != 3:
            raise Warning("Use other method. Only 3D geometries")
        
        start = time.time()
        # For C matrix we only need B0 and W0, we will define 
        DB = []
        DW = []
        DI = []
        for _ in range(self._dim):
            DB.append(self._DB[_][0])
            DW.append(self._DW[_][0][0])
            DI.append(self._DB[_][0]@self._DB[_][0].T)

        # Compute capacity matrix
        indi_C, indj_C = get_indexes_3D(DI)
        val_C = get_matrix_3D(self._capacity_coef, DB, DW, DI)

        C = self.array2csr_matrix(self._nb_ctrlpts_total, self._nb_ctrlpts_total,  
                                            val_C, indi_C, indj_C).tocsc()
        
        stop = time.time()
        print('Capacity matrix assembled in : %.5f s' %(stop-start))

        return 
