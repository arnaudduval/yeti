"""
.. New WQ methods 
.. by Joaquin Cornejo
"""

# Python libraries
import os
import numpy as np
from scipy import sparse as sp
from geomdl import helpers
import matplotlib.pyplot as plt

# My libraries
from lib.base_functions import (eval_basis_python, 
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

# # ================
# # CODE 
# # ================
# degree = 2
# nbel = 4
# knotvector = create_knotvector(degree, nbel)
# _, B0p, B1p, W0p, W1p = wq_find_basis_weights_opt(degree, knotvector, 2, maxrule= 2)
# print(W0p.toarray())

# CUTS = np.arange(2, 5)
# NBEL = [2**_ for _ in CUTS]

# for varName in ['I00', 'I01', 'I10', 'I11']:
#     plt.figure(1)
#     ax = plt.gca()
#     for degree in range(2, 6):
#         norm_python = []; ddl =[]
#         color = next(ax._get_lines.prop_cycler)['color']
#         for nbel in NBEL: 
#             # ========================================
#             # PYTHON
#             # ========================================
#             knotvector = create_knotvector(degree, nbel)
#             _, B0p, B1p, W0p, W1p = wq_find_basis_weights_opt(degree, knotvector, 2, maxrule= 2)

#             print(B0p.toarray())

#             # Calculate I
#             I00p = W0p @ B0p.transpose()
#             I01p = W0p @ B1p.transpose()
#             I10p = W1p @ B0p.transpose()
#             I11p = W1p @ B1p.transpose()

#             # ========================================
#             # REFERENCE
#             # ========================================
#             qp_cgg, Wcgg = iga_find_positions_weights(degree, knotvector)
#             B0, B1 = eval_basis(degree, knotvector, qp_cgg)

#             # Calculate I
#             I00 = B0 @ np.diag(Wcgg) @ B0.transpose()
#             I01 = B0 @ np.diag(Wcgg) @ B1.transpose()
#             I10 = B1 @ np.diag(Wcgg) @ B0.transpose()
#             I11 = B1 @ np.diag(Wcgg) @ B1.transpose()

#             # To choose variables
#             if varName == 'I00': var1 = I00; var2 = I00p
#             elif varName == 'I01': var1 = I01; var2 = I01p
#             elif varName == 'I10': var1 = I10; var2 = I10p
#             elif varName == 'I11': var1 = I11; var2 = I11p

#             # Compare results 
#             error_python = var1 - var2
#             norm_temp = np.linalg.norm(error_python, np.inf)/np.linalg.norm(var1, np.inf)
#             norm_python.append(norm_temp)

#         # Change type 
#         norm_python = np.asarray(norm_python)
#         ddl = np.asarray(NBEL)

#         # Figure 
#         plt.figure(1)
#         plt.plot(ddl, norm_python*100, '--P', color=color)

#     # Plot configurations
#     plt.grid()
#     plt.xscale("log")
#     plt.xlabel("Number of elements $nb_{el}$", fontsize= 16)
#     plt.yscale("log")
#     plt.ylabel("Relative error (%)", fontsize= 16)
#     plt.xlim([1, 100])
#     plt.xticks(fontsize=16)
#     plt.yticks(fontsize=16)
#     plt.tight_layout()

#     # Maybe we can save images in svg format (vectorized)
#     plt.savefig(folder + 'Error_basisweights2_' + varName +'.png')
#     plt.figure(1).clear()
