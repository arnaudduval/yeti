"""
.. This module is an adaptation of Matrix Free (MF) 
.. Weighted Quadrature (WQ) for IsoGeometric Analysis (MF-WQ-IGA)
.. Joaquin Cornejo 
"""

# Python libraries
from .__init__ import *

# My libraries
from .base_functions import (compute_eig_diag, 
                            eigen_decomposition, 
                            erase_rows_csr, 
                            fast_diagonalization, 
                            wq_find_basis_weights_fortran
)
from .create_model import thermoMechaModel
from .D3viscoplasticity import *

class fortran_mf_wq(thermoMechaModel):

    def __init__(self, modelIGA, material=None, Dirichlet=None, Neumann=None):

        super().__init__(modelIGA, material=material, Dirichlet=Dirichlet, Neumann=Neumann)

        # Set basis and weights
        self._nb_qp, self._nb_qp_total = np.ones(3, dtype=int), None
        self.eval_basis_weigths()
        
        # Get jacobian and physical position 
        self.eval_jacobien_physicalPosition()

        return

    def eval_basis_weigths(self): 
        " Computes basis and weights "

        print('Evaluating basis and weights')
        start = time.time()

        # Initialize 
        self._nnz_I, self._qp_dim, self._DB, self._DW, self._indices = [], [], [], [], []

        for dim in range(self._dim):  
            nnz_I, qp_position, basis, \
            weights, indi, indj = wq_find_basis_weights_fortran(self._degree[dim], self._knotvector[dim])
            self._nb_qp[dim] = len(qp_position)
            self._nnz_I.append(nnz_I); self._qp_dim.append(qp_position)
            self._DB.append(basis); self._DW.append(weights)
            self._indices.append(indi); self._indices.append(indj)

        # Update number of quadrature points
        self._nb_qp_total = np.prod(self._nb_qp)

        stop = time.time()
        print('\tBasis and weights in : %.5f s' %(stop-start))

        return

    def eval_jacobien_physicalPosition(self): 
        " Computes jacobien and physical position "

        print('Evaluating jacobien and physical position')
        start = time.time()
        
        # Get inputs
        inputs = [*self._nb_qp, *self._indices, *self._DB, self._ctrlpts]
        if self._dim == 2:
            self._Jqp, self._detJ, self._invJ = assembly.eval_jacobien_2d(*inputs)
            self._qp_PS = assembly.interpolate_fieldphy_2d(*inputs)
        if self._dim == 3:
            self._Jqp, self._detJ, self._invJ = assembly.eval_jacobien_3d(*inputs)
            self._qp_PS = assembly.interpolate_fieldphy_3d(*inputs)
        stop = time.time()
        print('\t Time jacobien: %.5f s' %(stop-start))

        return

    # ----------------------
    # WEIGHTED QUADRATURE
    # ----------------------

    def eval_source_coefficient(self, fun): 
        " Computes source coefficients "
        coefs = super().eval_source_coefficient(fun, self._qp_PS, self._detJ)
        return coefs

    def eval_bodyforce_coefficient(self, fun):
        " Computes body force coefficients "
        factor = self._detJ*self._density
        coefs = super().eval_source_coefficient(fun, self._qp_PS, factor)
        return coefs

    def eval_capacity_matrix(self, indi=None, indj=None): 
        " Computes capacity matrix "

        if indi is None: indi = np.arange(self._nb_ctrlpts_total, dtype=int)
        if indj is None: indj = np.arange(self._nb_ctrlpts_total, dtype=int)
        
        # Get inputs
        super()._verify_thermal()
        coefs = super().eval_capacity_coefficient(self._detJ, self._capacity)
        inputs = [coefs, *self._nb_qp, *self._indices, *self._DB, *self._DW, *self._nnz_I]

        start = time.time()
        if self._dim == 2: val_, indi_, indj_ = assembly.wq_get_capacity_2d(*inputs)
        if self._dim == 3: val_, indi_, indj_ = assembly.wq_get_capacity_3d(*inputs)
        matrix = super().array2csr_matrix(val_, indi_, indj_).tocsc()[indi, :][:, indj]
        stop = time.time()

        print('Capacity matrix assembled in : %.5f s' %(stop-start))
        
        return  matrix

    def eval_conductivity_matrix(self, indi= None, indj= None): 
        " Computes conductivity matrix "

        if indi is None: indi = np.arange(self._nb_ctrlpts_total, dtype=int)
        if indj is None: indj = np.arange(self._nb_ctrlpts_total, dtype=int)
        
        # Get inputs
        super()._verify_thermal()
        coefs = super().eval_conductivity_coefficient(self._invJ, self._detJ, self._conductivity)
        inputs = [coefs, *self._nb_qp, *self._indices, *self._DB, *self._DW, *self._nnz_I]

        start = time.time()
        if self._dim == 2: val_, indi_, indj_ = assembly.wq_get_conductivity_2d(*inputs)
        if self._dim == 3: val_, indi_, indj_ = assembly.wq_get_conductivity_3d(*inputs)
        matrix = super().array2csr_matrix(val_, indi_, indj_).tocsc()[indi, :][:, indj]
        stop = time.time()

        print('Conductivity matrix assembled in : %5f s' %(stop-start))
        
        return matrix

    def eval_stiffness_matrix(self, coefs=None): 
        " Computes conductivity matrix "
        
        if self._dim != 3: raise Warning('Not yet')
        super()._verify_mechanics()
        
        # Get inputs
        if coefs is None: coefs = super().eval_elastic_coefficient(self._Jqp)
        inputs = [coefs, *self._nb_qp, *self._indices, *self._DB, *self._DW, *self._nnz_I]
        start = time.time()
        val_, indi_, indj_ = assembly.wq_get_stiffness_3d(*inputs)
        matrix = super().array2coo_matrix(val_, indi_, indj_).tocsc()
        stop = time.time()

        print('Stiffness matrix assembled in : %5f s' %(stop-start))
        
        return matrix

    def eval_Ku(self, u, table=None): 
        " Computes K u where K is conductivity matrix "

        # Get inputs
        super()._verify_thermal()
        coefs = super().eval_conductivity_coefficient(self._invJ, self._detJ, self._conductivity)
        inputs = self.get_input4MatrixFree(table=table)
        if self._dim == 2: raise Warning('Until now not done')
        if self._dim == 3: result = solver.mf_wq_get_ku_3d_csr(coefs, *inputs, u)

        return result

    def eval_Cu(self, u, table= None): 
        " Computes C u where C is capacity matrix "

        # Get inputs
        super()._verify_thermal()
        coefs = super().eval_capacity_coefficient(self._detJ, self._capacity)
        inputs = self.get_input4MatrixFree(table=table)
        if self._dim == 2: raise Warning('Until now not done')
        if self._dim == 3: result = solver.mf_wq_get_cu_3d_csr(coefs, *inputs, u)

        return result

    def eval_Su(self, u, coefs=None):
        " Computes S u where S is stiffness matrix "

        # Get inputs
        if self._dim != 3: raise Warning('Until now not done')
        super()._verify_mechanics()
        if coefs is None: coefs = super().eval_elastic_coefficient(self._Jqp)
        inputs = [coefs, *self._nb_qp, *self._indices, *self._DB, *self._DW]
        result = solver.mf_wq_get_su_3d_csr(*inputs, u)

        return result

    def eval_source_vector(self, fun, indi=None): 
        " Computes source vector "

        if indi is None: indi = np.arange(self._nb_ctrlpts_total, dtype=int)

        # Get source coefficients
        coefs = self.eval_source_coefficient(fun)
        inputs = [coefs, *self._nb_qp, *self._indices, *self._DW]
        start = time.time()
        if self._dim == 2: vector = assembly.wq_get_source_2d(*inputs)[indi]
        if self._dim == 3: vector = assembly.wq_get_source_3d(*inputs)[indi]
        stop = time.time()

        print('Source vector assembled in : %.5f s' %(stop-start))

        return vector

    def eval_force_surf(self):
        " Returns force vector at the surface. In 3D: surface integrals. "
        ## For now, we only consider constants forces at the boundaries

        if self._dim != 3: raise Warning('Method only for 3D geometries')
        if self._mechanicalNeumann is None: raise Warning('Define Neumann boundaries')

        def get_info(nb):
            direction = int(np.floor(nb/2))
            if nb%2 == 1: side = 1
            else: side = 0
            return direction, side

        # Initialize 
        Ftemp = np.zeros((self._dim+1, self._nb_ctrlpts_total))

        # Get INC of control points and INC of quadrature points
        INC_CP = super().get_NURBScoordinates(self._nb_ctrlpts)
        INC_QP = super().get_NURBScoordinates(self._nb_qp)

        for _ in range(self._dim*2):
            # Get direction 
            direction, side = get_info(_)

            # Get force
            force = self._mechanicalNeumann[_]
            if np.array_equal(force, np.zeros(self._dim)): continue
            else:
                # Get control points and quadrature points list
                if side == 0: 
                    CPList = np.where(INC_CP[:, direction] == 0)[0]
                    QPList = np.where(INC_QP[:, direction] == 0)[0]

                elif side == 1: 
                    CPList = np.where(INC_CP[:, direction] == self._nb_ctrlpts[direction]-1)[0]
                    QPList = np.where(INC_QP[:, direction] == self._nb_qp[direction]-1)[0]
                
                # Order my list
                CPList = list(np.sort(CPList))
                QPList = list(np.sort(QPList))

                # Get modified Jacobien matrix
                valrange = [i for i in range(self._dim)]
                valrange.pop(direction)
                JJ = self._Jqp[:, :, QPList]
                JJ = JJ[:, valrange, :]

                # Get inputs to compute vector
                nnz, indices, data_W = [], [], []
                for _ in valrange:
                    nnz.append(self._nb_qp[_]); data_W.append(self._DW[_])
                    indices.append(self._indices[2*_]); indices.append(self._indices[2*_+1]) 
                
                # Compute surface force
                FSurf = solver.wq_get_forcesurf_3d(force, JJ, *nnz, *indices, *data_W)
                Ftemp[:-1, CPList] += FSurf
                Ftemp[-1, CPList] += 1

        # Final update of the vector (average)
        F = Ftemp[:-1, :]
        FList = np.where(Ftemp[-1, :] >= 2)[0]
        F[:, FList] /= Ftemp[-1, FList]

        return F

    def eval_force_body(self, fun):
        " Computes body (volumetric) force "

        # Get source coefficients
        if self._dim != 3: raise Warning('Method only for 3D geometries')
        coefs = self.eval_bodyforce_coefficient(fun)
        inputs = [coefs, *self._nb_qp, *self._indices, *self._DW]
        start = time.time()
        vector = solver.wq_get_forcevol_3d(*inputs)
        stop = time.time()

        print('Body force vector assembled in : %.5f s' %(stop-start))

        return vector

    def eval_diag_K(self): 
        " Computes the diagonal of conductivity matrix "

        # Get inputs
        super()._verify_thermal()
        coefs = super().eval_conductivity_coefficient(self._invJ, self._detJ, self._conductivity)
        inputs = [coefs, *self._nb_qp, *self._indices, *self._DB, *self._DW]
        start = time.time()
        if self._dim == 2: raise Warning('Until now not done')
        if self._dim == 3: vector = assembly.wq_find_conductivity_diagonal_3d(*inputs)
        stop = time.time()

        print('Conductivity matrix assembled in : %5f s' %(stop-start))
        
        return vector

    # ----------------------------------
    # MATRIX FREE SOLUTION (IN FORTRAN)
    # ----------------------------------   

    def get_input4MatrixFree(self, table=None):
        " Returns necessary inputs to compute the product between a matrix and a vector"

        # Initialize
        indices, data_B, data_W = [], [], []
        
        # Compute inputs
        if table is None: table = np.asarray([[0, 0], [0, 0], [0, 0]])
        for dim in range(self._dim):
            # Select data
            if np.array_equal(table[dim, :], [0, 0]): rows2erase = []
            if np.array_equal(table[dim, :], [0, 1]): rows2erase = [-1]
            if np.array_equal(table[dim, :], [1, 0]): rows2erase = [0]
            if np.array_equal(table[dim, :], [1, 1]): rows2erase = [0, -1]
            indi_t, indj_t, data_t = erase_rows_csr(rows2erase, 
                                    self._indices[2*dim], self._indices[2*dim+1],  
                                    [self._DB[dim], self._DW[dim]])
            
            # Extract data and append to list
            [dB, dW] = data_t
            indices.append(indi_t); indices.append(indj_t) 
            data_B.append(dB); data_W.append(dW)

        inputs = [*self._nb_qp, *indices, *data_B, *data_W]

        return inputs

    def MFsteadyHeat(self, f, nbIterations, epsilon, method, directsol=None): 
        " Solves steady heat problems using directly substitution method "

        if self._thermalDirichlet is None: raise Warning('Ill conditionned. It needs Dirichlet conditions')
        if directsol is None: directsol = np.ones(len(f))

        # Get inputs 
        super()._verify_thermal()
        coefs = super().eval_conductivity_coefficient(self._invJ, self._detJ, self._conductivity)
        inputs_tmp = self.get_input4MatrixFree(table=self._thermalDirichlet)
        inputs = [coefs, *inputs_tmp, f, nbIterations, epsilon, method, 
                    self._conductivity, self._Jqp, directsol]

        if self._dim == 2: raise Warning('Until now not done')
        if self._dim == 3: sol, residue, error = solver.mf_wq_steady_heat_3d(*inputs)

        return sol, residue, error

    def MFsteadyHeat_PLS(self, f, nbIterations, epsilon, ud=None, method_pls='S', methodPCG='TDC'): 
        " Solves steady heat with penalty, Lagrange or substitution method "

        # Get inputs 
        super()._verify_thermal()
        coefs = super().eval_conductivity_coefficient(self._invJ, self._detJ, self._conductivity)
        if self._thermalDirichlet is None: raise Warning('Ill conditionned. It needs Dirichlet conditions')
        if ud is None: ud = np.zeros(len(self._thermal_dod)) 

        # Convert Python to Fortran
        dod = np.copy(self._thermal_dod); dod += 1
        inputs = [coefs, *self._nb_qp, *self._indices, *self._DB, *self._DW, 
                    self._thermalDirichlet, dod, f, ud, nbIterations, epsilon]

        if self._dim == 2: raise Warning('Until now not done')
        if self._dim == 3: 
            if method_pls == 'S': 
                sol, residue = solver. mf_wq_solve_shs_3d(*inputs, 
                                methodPCG, self._conductivity, self._Jqp)   
            else: raise Warning('Method is not well implemented')

        return sol, residue

    def interpolate_ControlPoints(self, fun, nbIter=100, eps=1e-14):
        " Interpolation from parametric space to physical space "

        # Get coeficients 
        coefs = fun(self._qp_PS) * self._detJ

        # Get inputs
        inputs = [coefs, *self._nb_qp, *self._indices, *self._DW]

        # Calculate vector
        if self._dim == 2: raise Warning('Until now not done')
        if self._dim == 3: vector = assembly.wq_get_source_3d(*inputs)

        # Solve linear system with fortran
        inputs = [self._detJ, *self._nb_qp, *self._indices, *self._DB, *self._DW, vector, nbIter, eps]
        start = time.time()
        u_interp, relres = solver.mf_wq_interpolate_cp_3d(*inputs)
        stop = time.time()
        res_end = relres[np.nonzero(relres)][-1]
        print('Interpolation in: %.3e s with relative residue %.3e' %(stop-start, res_end))

        return u_interp

    def MFplasticity_fortran(self, Fext=None, indi=None):
        " Solves a plasticity problem "

        # Get inputs 
        if self._dim != 3: raise Warning('Not yet')
        super()._verify_mechanics()
        if self._mechanicalDirichlet is None: raise Warning('Ill conditionned. It needs Dirichlet conditions')
        if indi is None or Fext is None: raise Warning('Impossible')

        dod = deepcopy(indi)
        for _ in range(len(dod)):
            newdod = np.array(dod[_])
            newdod += 1
            dod[_] = list(newdod)

        prop = np.array([self._youngModule, self._hardening, self._betaHard, self._poissonCoef, self._sigmaY])       
        inputs = [*self._nb_qp, *self._indices, *self._DB, *self._DW, prop, self._mechanicalDirichlet, 
                    *dod, self._invJ, self._detJ]
        result = solver.mf_wq_plasticity_3d(*inputs, Fext)

        return result
    
    def MFelasticity_fortran(self, coefs=None, Fext=None, indi=None, nbIterations=300, isPrecond=True, isnoised=False):
        " Solves a elasticity problem "
        
        # Get inputs 
        if self._mechanicalDirichlet is None: raise Warning('Ill conditionned. It needs Dirichlet conditions')
        if indi is None or Fext is None: raise Warning('Impossible')
        super()._verify_mechanics()
        if coefs is None: coefs = super().eval_elastic_coefficient(self._Jqp, isnoised=isnoised)

        dod = deepcopy(indi)
        for _ in range(len(dod)):
            newdod = np.array(dod[_])
            newdod += 1
            dod[_] = list(newdod)

        inputs = [coefs, *self._nb_qp, *self._indices, *self._DB, *self._DW, 
                isPrecond, nbIterations, self._mechanicalDirichlet, *dod]
        result = solver.mf_wq_elasticity_3d_py(*inputs, Fext)

        return result

    # ----------------------------------
    # MATRIX FREE SOLUTION (IN PYTHON)
    # ---------------------------------- 

    def compute_eigen_all(self, table=None, ndof=3):
        """ Computes the eigen values and vectors considering Robin condition
            If Dirichlet condition is needed, then is better to try other function
        """

        # Initialize
        Deig = np.zeros((3, self._nb_ctrlpts_total))
        eigvec_U = np.zeros((self._nb_ctrlpts[0], self._nb_ctrlpts[0], 3))
        eigvec_V = np.zeros((self._nb_ctrlpts[1], self._nb_ctrlpts[1], 3))
        eigvec_W = np.zeros((self._nb_ctrlpts[2], self._nb_ctrlpts[2], 3))

        for idof in range(ndof):
            list_eig, list_vectors = [], []
            for dim in range(self._dim):
                indit = self._indices[2*dim]
                indjt =  self._indices[2*dim+1]
                B = self._DB[dim]
                W = self._DW[dim]
                datat = [B[:, 0], B[:, 1], W[:, 0], W[:, -1]]
                robin = table[dim, :, idof]
                eigt, Ut = eigen_decomposition(indit, indjt, datat, robin_condition=robin)
                list_eig.append(eigt)
                list_vectors.append(Ut)

            Deig[idof, :] = compute_eig_diag(*list_eig)
            eigvec_U[:, :, idof] = list_vectors[0]
            eigvec_V[:, :, idof] = list_vectors[1]
            eigvec_W[:, :, idof] = list_vectors[2]
        
        # Save data
        DU = [eigvec_U, eigvec_V, eigvec_W, Deig]

        return DU

    def compute_strain(self, u=None):
        " Compute strain field from displacement field "

        if self._dim != 3: raise Warning('Not yet')
        if u is None: raise Warning('Insert displacement')
        
        inputs = [*self._nb_qp, *self._indices, *self._DB, self._invJ, u]
        eps = solver.interpolate_strain_3d(*inputs)

        return eps

    def compute_internal_force(self, coefs=None):
        "Compute internal force using sigma coefficients "

        if self._dim != 3: raise Warning('Not yet')
        if coefs is None: raise Warning('Insert coefficients')

        inputs = [coefs, *self._nb_qp, *self._indices, *self._DW]
        Fint = solver.wq_get_forceint_3d(*inputs)

        return Fint

    def MFelasticity_py(self, coefs=None, DU=None, Fext=None, indi=None, 
                            nbIterations=200, tol=1e-7, isPrecond=True, isnoised=False):
        " Solve linear system using Bi-CG Stab algorithm for elasticity problems "

        # Initialize
        if self._dim != 3: raise Warning('Not yet')
        super()._verify_mechanics()
        if coefs is None: coefs = super().eval_elastic_coefficient(self._Jqp, isnoised=isnoised)
        
        x = np.zeros(np.shape(Fext)); r = Fext
        clean_dirichlet_3d(r, indi)
        rhat, p = r, r
        rsold = block_dot_product(self._dim, r, rhat)
        norm2b = np.linalg.norm(r)

        if not isPrecond: # Without preconditioner 
            
            for i in range(nbIterations):
                Ap = self.eval_Su(p, coefs=coefs)
                clean_dirichlet_3d(Ap, indi)

                alpha = rsold/block_dot_product(self._dim, Ap, rhat)
                s = r -alpha*Ap

                As = self.eval_Su(s, coefs=coefs)
                clean_dirichlet_3d(As, indi)

                omega = block_dot_product(self._dim, As, s)/block_dot_product(self._dim, As, As)
                x += alpha*p + omega*s
                r = s - omega*As

                relerror = np.sqrt(np.linalg.norm(r)/norm2b)
                print(relerror)
                if relerror <= tol: break

                rsnew = block_dot_product(self._dim, r, rhat)
                beta = (alpha/omega)*(rsnew/rsold)
                p = r + beta*(p - omega*Ap)
                rsold = rsnew

        else: # With preconditioner
            
            # Set values
            if DU is None: DU = self.compute_eigen_all(table=self._mechanicalDirichlet)
            U, V, W, D = DU[0], DU[1], DU[2], DU[3]
            
            for i in range(nbIterations):
                ptilde = fast_diagonalization(U, V, W, D, p, fdtype='elastic')
                clean_dirichlet_3d(ptilde, indi)

                Aptilde = self.eval_Su(ptilde, coefs=coefs)
                clean_dirichlet_3d(Aptilde, indi)

                alpha = rsold/block_dot_product(self._dim, Aptilde, rhat)
                s = r - alpha*Aptilde

                stilde = fast_diagonalization(U, V, W, D, s, fdtype='elastic')
                clean_dirichlet_3d(stilde, indi)

                Astilde = self.eval_Su(stilde, coefs=coefs)
                clean_dirichlet_3d(Astilde, indi)

                omega = block_dot_product(self._dim, Astilde, s)/block_dot_product(self._dim, Astilde, Astilde)
                x += alpha*ptilde + omega*stilde
                r = s - omega*Astilde

                relerror = np.sqrt(np.linalg.norm(r)/norm2b)
                print(relerror)
                if relerror <= tol: break

                rsnew = block_dot_product(self._dim, r, rhat)
                beta = (alpha/omega)*(rsnew/rsold)
                p = r + beta*(p - omega*Aptilde)
                rsold = rsnew

        print('After %d iteration, the relative residue is %f' %(i, relerror))

        return x

    def MFplasticity_py(self, Fext=None, indi=None, tol=1e-6, d=3):
        " Solves plasticity problem "

        if self._dim != 3: raise Warning('Only for 3D')
        super()._verify_mechanics()

        # Initialize
        ddl = int(d*(d+1)/2)
        ep_n0 = np.zeros(self._nb_qp_total)
        sigma_n0 = np.zeros((ddl, self._nb_qp_total))
        ep_n1 = np.zeros(self._nb_qp_total)
        sigma_n1 = np.zeros((ddl, self._nb_qp_total))
        Dalg = np.zeros((ddl, ddl, self._nb_qp_total))
        disp = np.zeros(np.shape(Fext))
        inputs = [self._Ctensor, self._sigmaY, self._lame_mu, self._betaHard, self._hardening, self._Idev]

        # Compute eigen values and vectors
        DU = self.compute_eigen_all(table=self._mechanicalDirichlet)

        # Solver
        for i in range(1, np.shape(Fext)[2]):

            # Initialize
            ddisp = np.zeros(np.shape(disp[:, :, i-1]))
            Fext_t = Fext[:, :, i]
            prod2 = block_dot_product(d, Fext_t, Fext_t)

            # Solver Newton-Raphson
            for j in range(2):
                print('Pas %d, iteration %d' %(i+1, j))

                # Compute strain as function of displacement
                deps = self.compute_strain(ddisp)
    
                # Closest point projection in perfect plasticity
                for k in range(self._nb_qp_total):
                    sigma_n1t, ep_n1t, Dalgt = cpp_combined_hardening(inputs, deps[:, k], sigma_n0[:, k], ep_n0[k])
                    sigma_n1[:, k], ep_n1[k], Dalg[:, :, k] = sigma_n1t, ep_n1t, Dalgt

                # Compute coefficients to compute Fint and Stiffness
                coef_Fint, coef_Stiff = compute_plasticity_coef(sigma_n1, Dalg, self._invJ, self._detJ, d=d)

                # Compute Fint 
                Fint = self.compute_internal_force(coef_Fint)
                dF = Fext_t - Fint
                clean_dirichlet_3d(dF, indi) 
                print("WARNING: FIND A NEW STOP CRITERION")
                prod1 = block_dot_product(d, dF, dF)
                relerror = np.sqrt(prod1/prod2)
                print('Relative error: %.5f' %relerror)

                if relerror <= self._sigmaY*tol:
                    break
                else:
                    ddisp = self.MFelasticity_py(coefs=coef_Stiff, DU=DU, indi=indi, Fext=dF)
                    # ddisp = self.MFelasticity(coefs=coef_Stiff, indi=indi, Fext=dF)

            # Set values
            disp[:, :, i] = disp[:, :, i-1] + ddisp
            ep_n0 = ep_n1
            sigma_n0 = sigma_n1

        return
