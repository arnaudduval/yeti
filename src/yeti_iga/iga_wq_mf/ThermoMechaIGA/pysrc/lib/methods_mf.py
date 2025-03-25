"""
.. WQ-MF methods
.. by Joaquin Cornejo
.. Disclaimer :: This module will hardly ever be used and 
..               and maybe it does not work as expected
"""

# Python libraries
import numpy as np
from scipy import sparse as sp

# My libraries
from .methods_wq import WQ

class MF(WQ): 
    
    def __init__(self, obj: None):
        super().__init__(obj)

    def eval_Ku(self, u, indi, indj):
        " Computes K u "

        # Initialize conductivity matrix 
        Ku = np.zeros(len(indi))

        for j in range(self._dim):
            beta = np.zeros(self._dim, dtype = int); beta[j] = 1
            Bt = 1 
            for dim in range(self._dim): 
                    bt = beta[dim]
                    Bt = sp.kron(self._DB[dim][bt], Bt)
            But = sp.csr.csr_matrix.dot(Bt.tocsr()[indj,:].T, u)

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
                Ku += sp.csr_matrix.dot(Wt.tocsr()[indi,:], But)

        return Ku

    def eval_Cu(self, u, indi, indj):
        " Computes C u "

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
        W = W * sp.diags(self._capacity_coef)
        Bu = sp.csr_matrix.dot(B.tocsr()[indj,:].T, u)
        Cu = sp.csr_matrix.dot(W.tocsr()[indi,:], Bu)

        return Cu

    def eval_F(self, fun, indi):
        " Computes F "

        # Get source coefficients
        self._source_coef = super().eval_source_coefficient(fun)

        # Initialize weights
        W = 1; 

        for dim in range(self._dim): 
            # Find weights
            W00 = self._DW[dim][0][0]

            # Find W
            W = sp.kron(W, W00)

        # Assemble F
        F = sp.csr_matrix.dot(W.tocsr()[indi,:], self._source_coef)

        return F
    
    def mf_wq_evaluate_Au(self, ui, dof):
        "Return (A x)i"
        # Here x is all the vector 
        # In this fist case A = K
        Axi = self.eval_Ku(ui, dof, dof)
        return Axi

    def mf_wq_conj_grad(self, bi, dof, nbIterations, epsilon): 

        fun = self.mf_wq_evaluate_Au
        x, residue = super().conjugate_gradient_python(fun, bi, dof, nbIterations, epsilon)

        return x, residue
        