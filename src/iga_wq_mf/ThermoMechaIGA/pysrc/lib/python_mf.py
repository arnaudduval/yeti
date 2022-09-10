"""
.. WQ-MF methods
.. by Joaquin Cornejo
.. Disclaimer :: This module will hardly ever be used and 
..               and maybe it does not work as expected
"""

from .__init__ import *

# My libraries
from .python_wq import WQ

class MF(WQ): 
    
    def __init__(self, modelIGA, material=None, Dirichlet=None):
        super().__init__(modelIGA, material=material, Dirichlet=Dirichlet)
        return

    def eval_Ku(self, u, indi=None, indj=None):
        " Computes K u "

        if indi is None: indi = np.arange(self._nb_ctrlpts_total, dtype=int)
        if indj is None: indj = np.arange(self._nb_ctrlpts_total, dtype=int)
        super()._verify_thermal()
        coefs = super().eval_conductivity_coefficient(self._invJ, self._detJ, self._conductivity)

        # Initialize
        Ku = np.zeros(len(indi))
        for j in range(self._dim):
            beta = np.zeros(self._dim, dtype = int); beta[j] = 1
            Bt = 1 
            for dim in range(self._dim): 
                    bt = beta[dim]
                    Bt = sp.kron(self._DB[dim][bt], Bt)
            But = sp.csr.csr_matrix.dot(Bt.tocsr()[indj,:].T, u)

            for i in range(self._dim):
                Cij = coefs[i, j, :]
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

    def eval_Cu(self, u, indi=None, indj=None):
        " Computes C u "

        if indi is None: indi = np.arange(self._nb_ctrlpts_total, dtype=int)
        if indj is None: indj = np.arange(self._nb_ctrlpts_total, dtype=int)
        super()._verify_thermal()
        coefs = super().eval_capacity_coefficient(self._detJ, self._capacity)
        
        B = 1; W = 1
        for dim in range(self._dim): 
            # Find basis and weights
            B0 = self._DB[dim][0]
            W00 = self._DW[dim][0][0]

            # Find W and B
            B = sp.kron(B0, B)
            W = sp.kron(W00, W)

        # Assemble C
        W = W * sp.diags(coefs)
        Bu = sp.csr_matrix.dot(B.tocsr()[indj,:].T, u)
        Cu = sp.csr_matrix.dot(W.tocsr()[indi,:], Bu)

        return Cu

    def eval_F(self, fun, indi=None):
        " Computes F "

        if indi is None: indi = np.arange(self._nb_ctrlpts_total, dtype=int)

        # Get source coefficients
        coefs = super().eval_source_coefficient(fun)

        W = 1; 
        for dim in range(self._dim): 
            # Find weights
            W00 = self._DW[dim][0][0]

            # Find W
            W = sp.kron(W, W00)

        # Assemble F
        vector = sp.csr_matrix.dot(W.tocsr()[indi,:], coefs)

        return vector
    
    def mf_wq_evaluate_Au(self, u, dof):
        "Return (A x)i"
        # Here x is all the vector 
        # In this fist case A = K
        Au = self.eval_Ku(u, dof, dof)
        return Au

    def conjugate_gradient(self, fun_Au, bi, dof, nbIterations=100, epsilon=1e-10):   
        " Evaluate K u at choosen equations "

        # -----------------------------
        # Conjugate Gradient algorithm
        # -----------------------------
        x = np.zeros(len(bi))
        r = bi
        p = r
        rsold = np.dot(r, r)
        RelRes = np.zeros(nbIterations+1)
        RelRes[0] = 1.0

        for k in range(nbIterations):
            Ap = fun_Au(p, dof)
            alpha = rsold/np.dot(p, Ap)
            x = x + alpha*p
            r = r - alpha*Ap
            RelRes[k+1] = (np.linalg.norm(r, np.inf)/np.linalg.norm(bi, np.inf))
            if RelRes[k+1]<epsilon: break
            
            rsnew = np.dot(r, r)
            p = r + rsnew/rsold * p
            rsold = rsnew

        return x, RelRes

    def MFsolver(self, bi, dof, nbIterations, epsilon): 

        fun = self.mf_wq_evaluate_Au
        x, residue = self.conjugate_gradient(fun, bi, dof, nbIterations, epsilon)

        return x, residue
        