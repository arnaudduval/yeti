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

def tensor_decomposition_3D(n_list, coefs_matrix: np.ndarray):

    # We consider that coefs is a 3 x 3 x nbpts table   
    # Number of dimensions
    dim = 3

    # Set shape of coefs
    n1, n2, n3 = n_list

    try: n1 = n1[0]; n2 = n2[0]; n3 = n3[0]
    except: pass

    # Get diagonal
    coefs = np.zeros((dim, n1*n2*n3))
    for _  in range(dim): 
        coefs[_, :] = coefs_matrix[_, _, :]

    # Initialize
    u1 = np.ones(n1)
    u2 = np.ones(n2)
    u3 = np.ones(n3)
    w1 = np.ones(n1)
    w2 = np.ones(n2)
    w3 = np.ones(n3)

    Vscript = np.zeros((n1, n2, n3))
    Wscript = np.zeros((2, n1, n2, n3))
    Nscript = np.zeros((n1, n2, n3))
    Mscript = np.zeros((n1, n2, n3))

    # Transform coefs to tensor of coefs
    coefstens = np.zeros((dim, n1, n2, n3))
    for k in range(dim):
        for i3 in range(n3):
            for i2 in range(n2):
                for i1 in range(n1):
                    pos = i1 + i2*n1 + i3*n1*n2
                    coefstens[k, i1, i2, i3] = coefs[k, pos]

    for iter in range(2):
        for k in range(dim):
            # Set Dscript
            for i3 in range(n3):
                for i2 in range(n2):
                    for i1 in range(n1):
                        U = [u1[i1], u2[i2], u3[i3]]
                        Vscript[i1, i2, i3] = coefstens[k, i1, i2, i3]*U[k]\
                                                /(U[0]*U[1]*U[2])

            # Update w
            if k == 0:
                for j in range(n1):
                    m = Vscript[j, :, :].min()
                    M = Vscript[j, :, :].max()
                    w1[j] = np.sqrt(m*M)
            elif k == 1: 
                for j in range(n2):
                    m = Vscript[:, j, :].min()
                    M = Vscript[:, j, :].max()
                    w2[j] = np.sqrt(m*M)        
            elif k == 2: 
                for j in range(n2):
                    m = Vscript[:, :, j].min()
                    M = Vscript[:, :, j].max()
                    w3[j] = np.sqrt(m*M)

        for k in range(dim):
            cont = -1
            # Compute Wscript temporary
            for l in range(dim):
                if k != l:
                    cont += 1
                    for i3 in range(n3):
                        for i2 in range(n2):
                            for i1 in range(n1):
                                U = [u1[i1], u2[i2], u3[i3]]
                                W = [w1[i1], w2[i2], w3[i3]]
                                Wscript[cont, i1, i2, i3] = coefstens[k, i1, i2, i3]*U[l]*U[k]\
                                                            /(U[0]*U[1]*U[2]*W[k])

            # Compute Nscript and Mscript
            for i3 in range(n3):
                for i2 in range(n2):
                    for i1 in range(n1): 
                        WWlk = [Wscript[_, i1, i2, i3] for _ in range(2)]
                        Nscript[i1, i2, i3] = min(WWlk)
                        Mscript[i1, i2, i3] = max(WWlk)

            # Update u
            if k == 0:
                for j in range(n1):
                    m = Nscript[j, :, :].min()
                    M = Mscript[j, :, :].max()
                    u1[j] = np.sqrt(m*M)
            elif k == 1: 
                for j in range(n2):
                    m = Nscript[:, j, :].min()
                    M = Mscript[:, j, :].max()
                    u2[j] = np.sqrt(m*M)        
            elif k == 2: 
                for j in range(n2):
                    m = Nscript[:, :, j].min()
                    M = Mscript[:, :, j].max()
                    u3[j] = np.sqrt(m*M) 

    return u1, u2, u3, w1, w2, w3

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
        x, residue = super().conjugate_gradient(fun, bi, dof, nbIterations, epsilon)

        return x, residue
        