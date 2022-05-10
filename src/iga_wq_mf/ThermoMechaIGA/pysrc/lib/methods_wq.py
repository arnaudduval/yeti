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
from .base_functions import create_knotvector, eval_basis_python
from .create_model import thermoMechaModel
from .methods_iga import iga_find_positions_weights

# Yeti libraries
from preprocessing.igaparametrization import IGAparametrization

def wq_get_shape_B(degree, nbel, r): 
    " Return the shape of basis in WQ approach "

    # Set number of functions 
    nbfunct = degree + nbel

    # Set number of quadrature points WQ:
    nbx = 2 * (nbel + degree + r) - 5

    # Set table of positions of quadrature points 
    tableOfPointsOnSpan = np.zeros((nbel, 2), dtype= int)
    tableOfPointsOnSpan[0, 0] = 0; tableOfPointsOnSpan[0, 1] = degree + r -1 
    tableOfPointsOnSpan[-1, 0] = nbx - (degree + r) ; tableOfPointsOnSpan[-1,1] = nbx - 1
     
    for i in range(1, nbel - 1):
        tableOfPointsOnSpan[i, 0] = tableOfPointsOnSpan[i - 1, 1]
        tableOfPointsOnSpan[i, 1] = tableOfPointsOnSpan[i, 0] + 2
    
    # Set table of functions on every element
    tableOfFunctionsOnSpan = []
    for i in range(nbel): 
        tableOfFunctionsOnSpan.append(np.arange(i, i + degree + 1, dtype= int))
    tableOfFunctionsOnElement = np.asarray(tableOfFunctionsOnSpan)

    # Set table of knot-span for every function
    tableOfSpansForFunction = []
    for i in range(nbfunct):
        indi, _ = np.asarray(tableOfFunctionsOnElement == i).nonzero()
        indi = np.sort(indi)
        tableOfSpansForFunction.append(indi)

    # Find B
    indi_B0, indj_B0 = [], []
    indi_B1, indj_B1 = [], []
    for i in range(nbfunct) : 
        spanOfFunction = tableOfSpansForFunction[i]
        spanLeft = np.min(spanOfFunction)
        spanRight = np.max(spanOfFunction)

        support = np.arange(tableOfPointsOnSpan[spanLeft,0], \
                            tableOfPointsOnSpan[spanRight,1] + 1, dtype=int)
        
        support_B0 = support
        if i == 0: 
            support_B0 = np.delete(support, [-1])
        elif i == nbfunct - 1:
            support_B0 = np.delete(support, [0])
        else : 
            support_B0 = np.delete(support, [0 , -1])
        
        indi_B0.extend(i*np.ones(len(support_B0), dtype= int))
        indj_B0.extend(support_B0)

        support_B1 = support
        if i == 0 or i == 1: 
            support_B1 = np.delete(support, [-1])
        elif i == nbfunct - 1 or i == nbfunct - 2:
            support_B1 = np.delete(support, [0])
        else : 
            support_B1 = np.delete(support, [0 , -1])
        
        indi_B1.extend(i*np.ones(len(support_B1), dtype= int))
        indj_B1.extend(support_B1)
    
    # Set shape of B0 and B1
    data_B0 = np.ones(len(indi_B0))
    B0shape = sp.csr_matrix((data_B0, (indi_B0, indj_B0)), shape=(nbfunct, nbx))
    
    data_B1 = np.ones(len(indi_B1))
    B1shape = sp.csr_matrix((data_B1, (indi_B1, indj_B1)), shape=(nbfunct, nbx))
    
    return B0shape, B1shape

def wq_solve_equation_system(B, I):
    " Return solution of matrix system B w = I "

    # Convert type of B to array
    B = B.toarray()
    I = I.toarray()

    # Solve system 
    # It does not matter if B is under, well or over determined
    sol1, _, _, _ = np.linalg.lstsq(B, I, rcond=None)

    # Set solution
    w = sol1.reshape((1, -1)).tolist()
    
    return w[0]

def wq_find_positions(degree, knotvector, r):
    " Return position of quadrature points in WQ approach "
    
    # Set knots
    kvUnique = np.unique(knotvector)

    # Set number of elements
    nbel = len(kvUnique) - 1

    # Initialize quadrature points positions
    xwq = []
    
    # Find position of quadrature points
    for i in range(nbel):
        # Select quadrature points
        if i == 0 or i == nbel - 1:
            # Case: boundary knot-spans
            xt = np.linspace(kvUnique[i], kvUnique[i+1], degree+r)[1:-1]
        else:
            # case: interior knot-spans
            xt = np.linspace(kvUnique[i], kvUnique[i+1], 3)[1:-1]
        
        # Assemble vector
        xwq = np.append(xwq, xt)
    
    #  Include knots of knot-vector
    xwq = np.append(xwq, kvUnique)

    # Sort vector
    xwq = np.sort(xwq).tolist()

    return xwq

def wq_find_weights(degree, knotvector, r):
    " Returns weights at quadrature points in WQ approach using Space p-1 technique "
    
    # Set number of elements
    nbel = len(np.unique(knotvector))-1

    # ------------------------------------
    # Degree p
    # ------------------------------------
    # Find positions and weights in IGA approach 
    xcgg, wcgg = iga_find_positions_weights(degree, knotvector)

    # Find basis function values at Gauss points 
    B0cgg_p0, B1cgg_p0 = eval_basis_python(degree, knotvector, xcgg)

    # Find positions in WQ approach
    xwq = wq_find_positions(degree, knotvector, r)

    # Find basis function values at WQ points
    B0wq_p0 = eval_basis_python(degree, knotvector, xwq)[0]

    # ------------------------------------
    # Degree p - 1
    # ------------------------------------
    # Create curve
    knotvector_p1 = knotvector[1:-1]

    # Find basis function values at Gauss points
    B0cgg_p1 = eval_basis_python(degree-1, knotvector_p1, xcgg)[0]

    # Find basis function values at WQ points
    B0wq_p1 = eval_basis_python(degree-1, knotvector_p1, xwq)[0]  

    # ------------------------------------
    # Integrals
    # ------------------------------------
    # Compute exact integral I00 = int Bi(p) Bj(p) dx to calculate W00
    I00 = B0cgg_p0 @ sp.diags(wcgg) @ B0cgg_p0.T
    
    # Compute exact integral I01 = int Bi(p) Bj'(p) dx to calculate W01
    I01 = B0cgg_p1 @ sp.diags(wcgg) @ B0cgg_p0.T 

    # Compute exact integral I01 = int Bi'(p) Bj(p) dx to calculate W01
    I10 = B0cgg_p0 @ sp.diags(wcgg) @ B1cgg_p0.T
    
    # Compute exact integral I11 = int Bi'(p) Bj'(p) dx to calculate W11
    I11 = B0cgg_p1 @ sp.diags(wcgg) @ B1cgg_p0.T

    # Set number of functions
    nbfunct = degree + nbel

    # Set number of WQ quadrature points 
    nbx = len(xwq)

    # Initialize data for W00, W01, W10 and W11
    data_W00, data_W01, data_W10, data_W11 = [], [], [], []
    
    # Find shape of B and I
    shape_B_p0 = wq_get_shape_B(degree, nbel, r)[0]
    shape_B_p1 = wq_get_shape_B(degree-1, nbel, r+1)[0]

    # -----------------------------------------------
    # Computation of W00
    # -----------------------------------------------
    Ishape = shape_B_p0 @ shape_B_p0.T 
    indi = []; indj = []
    for i in range(nbfunct):
        F = np.nonzero(Ishape[:, i])[0]
        P = np.nonzero(shape_B_p0[i, :])[1]

        Bt = B0wq_p0[np.ix_(F,P)]
        It = I00[np.ix_(F, [i])]
      
        data_W00.extend(wq_solve_equation_system(Bt, It))
        indi.extend(i*np.ones(len(P), dtype= int))
        indj.extend(P)
    W00 = sp.csr_matrix((data_W00, (indi, indj)), shape=(nbfunct, nbx))
    
    # -----------------------------------------------
    # Computation of W01
    # -----------------------------------------------
    Ishape = shape_B_p1 @ shape_B_p0.T 
    indi = []; indj = []
    for i in range(nbfunct):
        F = np.nonzero(Ishape[:,i])[0]
        P = np.nonzero(shape_B_p0[i, :])[1]

        Bt = B0wq_p1[np.ix_(F,P)]
        It = I01[np.ix_(F, [i])]
        
        data_W01.extend(wq_solve_equation_system(Bt, It))
        indi.extend(i*np.ones(len(P), dtype= int))
        indj.extend(P)
    W01 = sp.csr_matrix((data_W01, (indi, indj)), shape=(nbfunct, nbx))
    
    # -----------------------------------------------
    # Computation of W10
    # -----------------------------------------------
    Ishape = shape_B_p0 @ shape_B_p0.T 
    indi = []; indj = []
    for i in range(nbfunct):
        F = np.nonzero(Ishape[:, i])[0]
        P = np.nonzero(shape_B_p0[i, :])[1]

        Bt = B0wq_p0[np.ix_(F,P)]
        It = I10[np.ix_(F, [i])]

        data_W10.extend(wq_solve_equation_system(Bt, It))
        indi.extend(i*np.ones(len(P), dtype= int))
        indj.extend(P)

    W10 = sp.csr_matrix((data_W10, (indi, indj)), shape=(nbfunct, nbx))
    
    # -----------------------------------------------
    # Computation of W11
    # -----------------------------------------------
    Ishape = shape_B_p1 @ shape_B_p0.T
    indi = []; indj = []
    for i in range(nbfunct):
        F = np.nonzero(Ishape[:,i])[0]
        P = np.nonzero(shape_B_p0[i, :])[1]

        Bt = B0wq_p1[np.ix_(F,P)]
        It = I11[np.ix_(F, [i])]

        data_W11.extend(wq_solve_equation_system(Bt, It))
        indi.extend(i*np.ones(len(P), dtype= int))
        indj.extend(P)
    W11 = sp.csr_matrix((data_W11, (indi, indj)), shape=(nbfunct, nbx))
    
    return W00, W01, W10, W11 

def wq_find_basis_weights_opt(degree, knotvector, r): 
    " Return Basis and Weights computed in an efficient way "
    
    # Set number of elements
    nbel = len(np.unique(knotvector))-1

    # Set quadrature points
    xwq = wq_find_positions(degree, knotvector, r)

    if nbel <= degree + 3 : 
        B0, B1 = eval_basis_python(degree, knotvector, xwq)
        W00, W01, W10,  W11 = wq_find_weights(degree, knotvector, r)
    else : 
        # Create model
        nbel_model = degree + 3
        knotvector_model = create_knotvector(degree, nbel_model)
        xwq_model = wq_find_positions(degree, knotvector_model, r)
        B0_model, B1_model = eval_basis_python(degree, knotvector_model, xwq_model)
        W00_model, W01_model, W10_model,  W11_model = wq_find_weights(degree, knotvector_model, r)
        shape_B0_model, shape_B1_model = wq_get_shape_B(degree, nbel, r)

        # Modify the results
        W00_model *= nbel_model / nbel
        W01_model *= nbel_model / nbel
        B1_model *= nbel / nbel_model

        #------------------------------
        # Algorithm to get the results
        #------------------------------
        # Set number of points
        nbx = 2*(nbel + degree + r) - 5

        # Set number of functions
        nbfunct = degree + nbel

        # Initialize weight matrices 
        data_B0, data_B1, data_W00, data_W01, data_W10, data_W11 = [], [], [], [], [], []

        # Find shape of B
        shape_B0, shape_B1 = wq_get_shape_B(degree, nbel, r)

        function_B0 = []; support_B0 = []
        function_B1 = []; support_B1 = []

        for i in range(degree+1):
            # For B0 and normal
            indj = shape_B0.getrow(i).nonzero()[1]
            function_B0.extend(i * np.ones(len(indj), dtype=int))
            support_B0.extend(indj)
            data_B0.extend(B0_model[np.ix_([i],indj)].toarray()[0].tolist())
            data_W00.extend(W00_model[np.ix_([i],indj)].toarray()[0].tolist())
            data_W01.extend(W01_model[np.ix_([i],indj)].toarray()[0].tolist())

            # For B0 and flip
            indj_flip = shape_B0.getrow(nbfunct-i-1).nonzero()[1]
            function_B0.extend((nbfunct-i-1) * np.ones(len(indj_flip), dtype=int))
            support_B0.extend(indj_flip)
            data_B0.extend(np.flip(B0_model[np.ix_([i],indj)]).toarray()[0].tolist())         
            data_W00.extend(np.flip(W00_model[np.ix_([i],indj)]).toarray()[0].tolist())   
            data_W01.extend(np.flip(W01_model[np.ix_([i],indj)]).toarray()[0].tolist())   
            # -------------------------
            
            # For B1 and normal
            indj = shape_B1.getrow(i).nonzero()[1]
            function_B1.extend(i * np.ones(len(indj), dtype=int))
            support_B1.extend(indj)
            data_B1.extend(B1_model[np.ix_([i],indj)].toarray()[0].tolist())
            data_W10.extend(W10_model[np.ix_([i],indj)].toarray()[0].tolist())
            data_W11.extend(W11_model[np.ix_([i],indj)].toarray()[0].tolist())

            # For B1 and flip
            indj_flip = shape_B1.getrow(nbfunct-i-1).nonzero()[1]
            function_B1.extend((nbfunct-i-1) * np.ones(len(indj_flip), dtype=int))
            support_B1.extend(indj_flip)
            data_B1.extend(np.flip(-B1_model[np.ix_([i],indj)]).toarray()[0].tolist())         
            data_W10.extend(np.flip(-W10_model[np.ix_([i],indj)]).toarray()[0].tolist())   
            data_W11.extend(np.flip(-W11_model[np.ix_([i],indj)]).toarray()[0].tolist())   

        for i in range(degree+1,nbfunct-degree-1): 
            # For B0
            indj_model = shape_B0_model.getrow(degree+1).nonzero()[1]
            indj = shape_B0.getrow(i).nonzero()[1]
            function_B0.extend(i * np.ones(len(indj), dtype=int))
            support_B0.extend(indj)
            data_B0.extend(B0_model[np.ix_([degree+1],indj_model)].toarray()[0].tolist())
            data_W00.extend(W00_model[np.ix_([degree+1],indj_model)].toarray()[0].tolist())
            data_W01.extend(W01_model[np.ix_([degree+1],indj_model)].toarray()[0].tolist())

            # For B1
            indj_model = shape_B1_model.getrow(degree+1).nonzero()[1]
            indj = shape_B1.getrow(i).nonzero()[1]
            function_B1.extend(i * np.ones(len(indj), dtype=int))
            support_B1.extend(indj)
            data_B1.extend(B1_model[np.ix_([degree+1],indj_model)].toarray()[0].tolist())
            data_W10.extend(W10_model[np.ix_([degree+1],indj_model)].toarray()[0].tolist())
            data_W11.extend(W11_model[np.ix_([degree+1],indj_model)].toarray()[0].tolist())
    
        B0 = sp.csr_matrix((data_B0, (function_B0, support_B0)), shape=(nbfunct, nbx))
        W00 = sp.csr_matrix((data_W00, (function_B0, support_B0)), shape=(nbfunct, nbx))    
        W01 = sp.csr_matrix((data_W01, (function_B0, support_B0)), shape=(nbfunct, nbx))    

        B1 = sp.csr_matrix((data_B1, (function_B1, support_B1)), shape=(nbfunct, nbx))    
        W10 = sp.csr_matrix((data_W10, (function_B1, support_B1)), shape=(nbfunct, nbx))    
        W11 = sp.csr_matrix((data_W11, (function_B1, support_B1)), shape=(nbfunct, nbx))  

    return xwq, B0, B1, W00, W01, W10, W11

class WQ(thermoMechaModel): 

    # ===========================
    # INITIALIZE 
    # ===========================

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
        source_coef = super().eval_F_coefficient(fun, self._dim, self._qp_PS, self._detJ)
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
