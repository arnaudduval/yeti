"""
.. IGA methods
.. by Joaquin Cornejo
.. Disclaimer :: This module will hardly ever be used and 
..               and maybe it does not work as expected
"""

from .__init__ import *

# My libraries
from .base_functions import iga_find_basis_weights_opt
from .create_model import thermoMechaModel

class IGA(thermoMechaModel):

    def __init__(self, modelIGA, material={}, Dirichlet={}):

        super().__init__(modelIGA, material=material, Dirichlet=Dirichlet)

        # Evaluate basis and weights
        self._nb_qp, self._nb_qp_total = np.ones(3, dtype=int), None
        self.eval_basis_weights()

        # Get jacobian and physical position
        self._Jqp, self._detJ, self._invJ, self._qp_PS = super().eval_jacobien_physicalPosition(self._dim, 
                                                        self._nb_qp_cgg_total, self._ctrlpts, self._DB)
        
        return

    def eval_basis_weights(self): 
        " Computes Basis and weights in WQ approach "
        
        print('Evaluating basis and weights')
        start = time.time()

        # Initalize 
        self._qp_dim, self._DB, self._DW = [], [], []

        for dim in range(self._dim): 
            qp_position, B0, B1, W = iga_find_basis_weights_opt(self._degree[dim], self._knotvector[dim]) 
            self._nb_qp[dim] = len(qp_position)
            self._qp_dim.append(qp_position)
            self._DB.append([B0, B1]); self._DW.append(W)
        
        # Update number of quadrature points
        self._nb_qp_total = np.prod(self._nb_qp)

        stop = time.time()
        print('\tBasis and weights in : %.5f s' %(stop-start))

    # ---------------
    # ASSEMBLE 
    # ---------------

    def eval_source_coefficient(self, fun): 
        " Computes source coefficients "
        coefs = super().eval_source_coefficient(fun, self._qp_PS, self._detJ)
        return coefs

    def eval_conductivity_matrix(self):
        " Assemble conductivity matrix K "

        start = time.time()
        
        # Initialize 
        super()._verify_thermal()
        coefs = super().eval_conductivity_coefficient(self._invJ, self._detJ, self._conductivity)
        matrix = sp.csr_matrix((self._nb_ctrlpts_total, self._nb_ctrlpts_total))

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
                Cij = coefs[i, j, :]
                alpha = np.zeros(self._dim, dtype = int); alpha[i] = 1
                
                Btl = 1
                for dim in range(self._dim):
                    at = alpha[dim] 
                    Btl = sp.kron(self._DB[dim][at], Btl)
                    
                # Evaluates Cij * B in each point
                Btr_Cij = sp.csr_matrix.dot(Btr, sp.diags(Cij))

                # Find K
                BW = sp.csr_matrix.dot(Btl.tocsr()[:,:], Wt.tocsr()[:,:])
                matrix += sp.csr_matrix.dot(BW.tocsr()[:,:], Btr_Cij.tocsr()[:,:].T)

        stop = time.time()
        print('Conductivity matrix assembled in : %.5f s' %(stop-start))

        return matrix

    def eval_capacity_matrix(self):
        " Assemble capacity matrix C "

        start = time.time()
        
        # Initialize 
        super()._verify_thermal()
        coefs = super().eval_capacity_coefficient(self._detJ, self._capacity)
        
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
        W = sp.csr_matrix.dot(W, sp.diags(coefs))
        matrix = sp.csr_matrix.dot(sp.csr_matrix.dot(B.tocsr()[:,:], sp.diags(W)), B.tocsr()[:,:].T)

        stop = time.time()
        print('Capacity matrix assembled in : %.5f s' %(stop-start))

        return matrix

    def eval_source_vector(self, fun):
        " Assemble power density vector F "

        start = time.time()
        
        # Get source coefficients
        coefs = self.eval_source_coefficient(fun)

        B = 1; W = 1
        for dim in range(self._dim): 
            # Find basis
            B0 = self._DB[dim][0]

            # Find weights
            W00 = self._DW[dim]

            # Find W and B
            W = sp.kron(sp.diags(W00), W)
            B = sp.kron(B0, B)

        # Assemble vector
        vector = sp.csr_matrix.dot(sp.csr_matrix.dot(B.tocsr()[:,:], W.tocsr()[:,:]), coefs)
        stop = time.time()
        print('Source vector assembled in : %.5f s' %(stop-start))

        return vector
