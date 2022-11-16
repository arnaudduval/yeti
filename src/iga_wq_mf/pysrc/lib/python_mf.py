"""
.. WQ-MF methods
.. by Joaquin Cornejo
.. Disclaimer :: This module will hardly ever be used and 
..               and maybe it does not work as expected
"""

from .__init__ import *

# My libraries
from .python_wq import WQ
from .base_functions import CG

class MF(WQ): 
	
	def __init__(self, modelIGA, material={}, Dirichlet={}):
		super().__init__(modelIGA, material=material, Dirichlet=Dirichlet)
		return

	def eval_Ku(self, u, indi=None, indj=None):
		" Computes K u "

		if indi is None: indi = np.arange(self._nb_ctrlpts_total, dtype=int)
		if indj is None: indj = np.arange(self._nb_ctrlpts_total, dtype=int)
		super()._verify_thermal()
		coefs = super().eval_conductivity_coefficient(self._invJ, self._detJ, self._conductivity)

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
			B0 = self._DB[dim][0]
			W00 = self._DW[dim][0][0]
			B = sp.kron(B0, B)
			W = sp.kron(W00, W)

		W = W * sp.diags(coefs)
		Bu = sp.csr_matrix.dot(B.tocsr()[indj,:].T, u)
		Cu = sp.csr_matrix.dot(W.tocsr()[indi,:], Bu)

		return Cu

	def eval_F(self, fun, indi=None):
		" Computes F "

		if indi is None: indi = np.arange(self._nb_ctrlpts_total, dtype=int)
		coefs = super().eval_source_coefficient(fun)

		W = 1
		for dim in range(self._dim): 
			W00 = self._DW[dim][0][0]
			W = sp.kron(W, W00)

		vector = sp.csr_matrix.dot(W.tocsr()[indi,:], coefs)

		return vector
	
	def mf_wq_evaluate_Au(self, u, dof):
		"Return (A x)i"
		# Here x is all the vector 
		# In this fist case A = K
		Au = self.eval_Ku(u, dof, dof)
		return Au
