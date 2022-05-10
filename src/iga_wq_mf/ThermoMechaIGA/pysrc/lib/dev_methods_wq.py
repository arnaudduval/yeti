"""
.. New WQ methods
.. All functions were tested but not implemented in principal code
.. by Joaquin Cornejo
"""

# Python libraries
import os
import numpy as np
from scipy import sparse as sp
import matplotlib.pyplot as plt

# My libraries
from .base_functions import (eval_basis_python, 
                                create_knotvector,
                                wq_find_positions,
                                wq_get_shape_B,
                                wq_solve_equation_system,
                                iga_find_positions_weights
)

# Choose folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/'
if not os.path.isdir(folder):
    os.mkdir(folder)

def wq_find_weights(degree, knotvector, r, maxrule= 1): 

    # Set number of elements
    nbel = len(np.unique(knotvector))-1
    # ------------------------------------
    # Degree p and continuity p-1
    # ------------------------------------
    # Find positions and weights in IGA approach 
    xcgg, wcgg = iga_find_positions_weights(degree, knotvector)

    # Find basis function values at Gauss points 
    B0cgg_p0, B1cgg_p0 = eval_basis_python(degree, knotvector, xcgg)

    # Find positions in WQ approach
    xwq = wq_find_positions(degree, knotvector, r, maxrule= maxrule)

    # ------------------------------------
    # Degree p and continuity 2
    # ------------------------------------
    # Create curve
    knotvector_p1 = create_knotvector(degree, nbel, multiplicity= 2)

    # Find basis function values at Gauss points
    B0cgg_p1 = eval_basis_python(degree, knotvector_p1, xcgg, multiplicity= 2)[0]

    # Find basis function values at WQ points
    B0wq_p1 = eval_basis_python(degree, knotvector_p1, xwq, multiplicity= 2)[0]  

    # ------------------------------------
    # Integrals
    # ------------------------------------
    # Compute exact integral I01 = int Bi(p) Bj'(p) dx to calculate W01
    I0 = B0cgg_p1 @ sp.diags(wcgg) @ B0cgg_p0.T 

    # Compute exact integral I11 = int Bi'(p) Bj'(p) dx to calculate W11
    I1 = B0cgg_p1 @ sp.diags(wcgg) @ B1cgg_p0.T

    # Set number of functions
    nbfunct = degree + nbel

    # Set number of WQ quadrature points 
    nbx = len(xwq)

    # Initialize data for W00, W01, W10 and W11
    data_W0, data_W1 = [], []
    
    # Find shape of B and I
    shape_B_p0 = wq_get_shape_B(degree, nbel, r, maxrule= maxrule)[0]
    shape_B_p1 = wq_get_shape_B(degree, nbel, r, multiplicity= 2, maxrule= maxrule)[0]
    Ishape = shape_B_p1 @ shape_B_p0.T

    # -----------------------------------------------
    # Computation of W00
    # -----------------------------------------------
    indi = []; indj = []
    for i in range(nbfunct):
        F = np.nonzero(Ishape[:, i])[0]
        P = np.nonzero(shape_B_p0[i, :])[1]

        Bt = B0wq_p1[np.ix_(F,P)]
        It = I0[np.ix_(F, [i])]
        Wt = wq_solve_equation_system(Bt, It)

        data_W0.extend(Wt)
        indi.extend(i*np.ones(len(P), dtype= int))
        indj.extend(P)
    W0 = sp.csr_matrix((data_W0, (indi, indj)), shape=(nbfunct, nbx))
    
    # -----------------------------------------------
    # Computation of W11
    # -----------------------------------------------
    indi = []; indj = []
    for i in range(nbfunct):
        F = np.nonzero(Ishape[:,i])[0]
        P = np.nonzero(shape_B_p0[i, :])[1]

        Bt = B0wq_p1[np.ix_(F,P)]
        It = I1[np.ix_(F, [i])]

        data_W1.extend(wq_solve_equation_system(Bt, It))
        indi.extend(i*np.ones(len(P), dtype= int))
        indj.extend(P)
    W1 = sp.csr_matrix((data_W1, (indi, indj)), shape=(nbfunct, nbx))
    
    return W0, W1

def wq_find_basis_weights_opt(degree, knotvector, r, maxrule= 1): 
    " Return Basis and Weights computed in an efficient way "
    
    # Set number of elements
    nbel = len(np.unique(knotvector))-1

    # Set quadrature points
    xwq = wq_find_positions(degree, knotvector, r, maxrule= maxrule)

    if nbel <= degree + 3 : 
        B0, B1 = eval_basis_python(degree, knotvector, xwq)
        W0, W1 = wq_find_weights(degree, knotvector, r, maxrule= maxrule)
    else : 
        # Create model
        nbel_model = degree + 3
        knotvector_model = create_knotvector(degree, nbel_model)
        xwq_model = wq_find_positions(degree, knotvector_model, r, maxrule= maxrule)
        B0_model, B1_model = eval_basis_python(degree, knotvector_model, xwq_model)
        W0_model, W1_model = wq_find_weights(degree, knotvector_model, r, maxrule= maxrule)
        shape_B0_model, shape_B1_model = wq_get_shape_B(degree, nbel, r, maxrule= maxrule)

        # Modify the results
        W0_model *= nbel_model / nbel
        B1_model *= nbel / nbel_model

        #------------------------------
        # Algorithm to get the results
        #------------------------------
        # Set number of points
        nbx = 2*(degree + r) + nbel*(maxrule + 1) - 2*maxrule - 3

        # Set number of functions
        nbfunct = degree + nbel

        # Initialize weight matrices 
        data_B0, data_B1, data_W0, data_W1 = [], [], [], []

        # Find shape of B
        shape_B0, shape_B1 = wq_get_shape_B(degree, nbel, r, maxrule= maxrule)

        function_B0 = []; support_B0 = []
        function_B1 = []; support_B1 = []

        for i in range(degree+1):
            # For B0 and normal
            indj = shape_B0.getrow(i).nonzero()[1]
            function_B0.extend(i * np.ones(len(indj), dtype=int))
            support_B0.extend(indj)
            data_B0.extend(B0_model[np.ix_([i],indj)].toarray()[0].tolist())
            data_W0.extend(W0_model[np.ix_([i],indj)].toarray()[0].tolist())

            # For B0 and flip
            indj_flip = shape_B0.getrow(nbfunct-i-1).nonzero()[1]
            function_B0.extend((nbfunct-i-1) * np.ones(len(indj_flip), dtype=int))
            support_B0.extend(indj_flip)
            data_B0.extend(np.flip(B0_model[np.ix_([i],indj)]).toarray()[0].tolist())         
            data_W0.extend(np.flip(W0_model[np.ix_([i],indj)]).toarray()[0].tolist())   
            # -------------------------
            
            # For B1 and normal
            indj = shape_B1.getrow(i).nonzero()[1]
            function_B1.extend(i * np.ones(len(indj), dtype=int))
            support_B1.extend(indj)
            data_B1.extend(B1_model[np.ix_([i],indj)].toarray()[0].tolist())
            data_W1.extend(W1_model[np.ix_([i],indj)].toarray()[0].tolist())

            # For B1 and flip
            indj_flip = shape_B1.getrow(nbfunct-i-1).nonzero()[1]
            function_B1.extend((nbfunct-i-1) * np.ones(len(indj_flip), dtype=int))
            support_B1.extend(indj_flip)
            data_B1.extend(np.flip(-B1_model[np.ix_([i],indj)]).toarray()[0].tolist())         
            data_W1.extend(np.flip(-W1_model[np.ix_([i],indj)]).toarray()[0].tolist())   

        for i in range(degree+1,nbfunct-degree-1): 
            # For B0
            indj_model = shape_B0_model.getrow(degree+1).nonzero()[1]
            indj = shape_B0.getrow(i).nonzero()[1]
            function_B0.extend(i * np.ones(len(indj), dtype=int))
            support_B0.extend(indj)
            data_B0.extend(B0_model[np.ix_([degree+1],indj_model)].toarray()[0].tolist())
            data_W0.extend(W0_model[np.ix_([degree+1],indj_model)].toarray()[0].tolist())

            # For B1
            indj_model = shape_B1_model.getrow(degree+1).nonzero()[1]
            indj = shape_B1.getrow(i).nonzero()[1]
            function_B1.extend(i * np.ones(len(indj), dtype=int))
            support_B1.extend(indj)
            data_B1.extend(B1_model[np.ix_([degree+1],indj_model)].toarray()[0].tolist())
            data_W1.extend(W1_model[np.ix_([degree+1],indj_model)].toarray()[0].tolist())
    
        B0 = sp.csr_matrix((data_B0, (function_B0, support_B0)), shape=(nbfunct, nbx))
        W0 = sp.csr_matrix((data_W0, (function_B0, support_B0)), shape=(nbfunct, nbx))    

        B1 = sp.csr_matrix((data_B1, (function_B1, support_B1)), shape=(nbfunct, nbx))    
        W1 = sp.csr_matrix((data_W1, (function_B1, support_B1)), shape=(nbfunct, nbx))  

    return xwq, B0, B1, W0, W1
