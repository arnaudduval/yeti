"""
.. IGA methods
.. by Joaquin Cornejo
.. Disclaimer :: This module will hardly ever be used and 
..               and maybe it does not work as expected
"""

# Python libraries
import numpy as np
from scipy import sparse as sp
import time

# My libraries
from .base_functions import iga_find_basis_weights_opt
from .create_model import thermoMechaModel

class IGA(thermoMechaModel):

    def __init__(self, modelIGA: None, material=None, Dirichlet=None):

        super().__init__(modelIGA, material= material, Dirichlet= Dirichlet)

        # Evaluate basis and weights
        self.eval_basis_weights()

        # Get jacobian and physical position
        self._Jqp, self._qp_PS, self._detJ = super().eval_jacobien_physicalPosition(self._dim, 
                                                        self._nb_qp_cgg_total, self._ctrlpts, self._DB)
        
        # Initialize thermal properties
        self._conductivity_coef = None
        self._capacity_coef = None
        
        return

    def eval_basis_weights(self): 
        " Computes Basis and weights in WQ approach "
        
        print('Evaluating basis and weights')
        start = time.time()

        # Initalize storage vectors
        self._qp_cgg, self._DB, self._DW = [], [], []

        for dim in range(self._dim): 
            # Find basis and weights 
            qp_pos, B0, B1, W = iga_find_basis_weights_opt(self._degree[dim], self._knotvector[dim]) 

            # Save 
            self._qp_cgg.append(qp_pos)
            self._DB.append([B0, B1])
            self._DW.append(W)

        stop = time.time()
        print('\tBasis and weights in : %.5f s' %(stop-start))

    def _verify_conductivity(self): 
        " Verifies if conductivity exits"

        if self._conductivity is None: raise Warning('Conductivity not defined')
        if self._conductivity_coef is None:
            print('Getting conductivity coefficients')
            start = time.time()
            coef, info = super().eval_conductivity_coefficient(self._Jqp, self._conductivity)
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
            coef, info = super().eval_capacity_coefficient(self._Jqp, self._capacity)
            if info == 0: raise Warning('Something happen computing coefficients')
            else: self._capacity_coef = coef
            stop = time.time()
            print('\tCapacity coefficients in : %.5f s' %(stop-start))
        return

    # ===========================
    # ASSEMBLE 
    # ===========================

    def eval_source_coefficient(self, fun): 
        " Computes source coefficients "
        source_coef = super().eval_source_coefficient(fun, self._qp_PS, self._detJ)
        return source_coef

    def eval_conductivity_matrix(self):
        " Assemble conductivity matrix K "

        start = time.time()
        
        # Initialize conductivity matrix
        self._verify_conductivity()
        K = sp.csr_matrix((self._nb_ctrlpts_total, self._nb_ctrlpts_total))

        Wt = 1
        for dim in range(self._dim):
            Wt = sp.kron(sp.diags(self._DW[dim]), Wt)

        for j in range(self._dim):
            beta = np.zeros(self._dim, dtype = int); beta[j] = 1
            Btr = 1 
            for dim in range(self._dim):
                bt = beta[dim]
                Btr = sp.kron(self._DB[dim][bt], Btr)
            
            for i in range(self._dim):
                Cij = self._conductivity_coef[i, j, :]
                alpha = np.zeros(self._dim, dtype = int); alpha[i] = 1
                
                Btl = 1
                for dim in range(self._dim):
                    at = alpha[dim] 
                    Btl = sp.kron(self._DB[dim][at], Btl)
                    
                # Evaluates Cij * B in each point
                Btr_Cij = sp.csr_matrix.dot(Btr, sp.diags(Cij))

                # Find K
                BW = sp.csr_matrix.dot(Btl.tocsr()[:,:], Wt.tocsr()[:,:])
                K += sp.csr_matrix.dot(BW.tocsr()[:,:], Btr_Cij.tocsr()[:,:].T)

        stop = time.time()
        print('Conductivity matrix assembled in : %.5f s' %(stop-start))

        return K

    def eval_capacity_matrix(self):
        " Assemble capacity matrix C "

        start = time.time()
        
        # Initialize weights
        self._verify_capacity()
        B = 1; W = 1
        for dim in range(self._dim): 
            # Find basis
            B0 = self._DB[dim][0]

            # Find weights
            W00 = self._DW[dim]

            # Find W and B
            B = sp.kron(B0, B)
            W = np.kron(np.array(W00), W)
            
        # Assemble C
        W = sp.csr_matrix.dot(W, sp.diags(self._capacity_coef))
        C = sp.csr_matrix.dot(sp.csr_matrix.dot(B.tocsr()[:,:], sp.diags(W)), B.tocsr()[:,:].T)

        stop = time.time()
        print('Capacity matrix assembled in : %.5f s' %(stop-start))

        return C

    def eval_source_vector(self, fun):
        " Assemble power density vector F "

        start = time.time()
        
        # Get source coefficients
        self._source_coef = self.eval_source_coefficient(fun)

        # Initialize weights
        B = 1; W = 1
        for dim in range(self._dim): 
            # Find basis
            B0 = self._DB[dim][0]

            # Find weights
            W00 = self._DW[dim]

            # Find W and B
            W = sp.kron(sp.diags(W00), W)
            B = sp.kron(B0, B)

        # Assemble F
        F = sp.csr_matrix.dot(sp.csr_matrix.dot(B.tocsr()[:,:], W.tocsr()[:,:]), self._source_coef)

        stop = time.time()
        print('Source vector assembled in : %.5f s' %(stop-start))

        return F
