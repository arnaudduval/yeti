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
from lib.methods_iga import iga_find_positions_weights

# Choose folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/'
if not os.path.isdir(folder):
    os.mkdir(folder)

def create_knotvector(p, nbel, multiplicity= 1):
    " Creates an uniform and open knot-vector "

    # Set knot-vector to be inserted
    knotvector_Unique = np.linspace(0., 1., nbel + 1)[1 : -1]

    # Create knot-vector 
    knotvector = []
    for _ in range(p+1): 
        knotvector.append(0.0)

    for knot in knotvector_Unique: 
        for _ in range(multiplicity): 
            knotvector.append(knot)

    for _ in range(p+1): 
        knotvector.append(1.0)
    
    return knotvector

def eval_basis(degree, knotvector, knots, multiplicity= 1): 
    " Evaluates B-spline functions at given knots "

    # Find number of points x
    nbx = len(knots)

    # Find number of elements 
    nbel = len(np.unique(knotvector)) - 1

    # Find number of functions 
    nbfunct = degree + multiplicity*(nbel - 1) + 1

    # Set table of functions per element 
    table_functions_element = np.zeros((nbel, degree + 2), dtype= int); 
    table_functions_element[0, 0] = degree; table_functions_element[0, 1:] = np.arange(degree + 1) 

    for _ in range(1, nbel): 
        # Set values of the table
        table_functions_element[_, :2] = table_functions_element[_-1, :2] + multiplicity
        table_functions_element[_, 2:] = table_functions_element[_, 1] + np.arange(1, degree + 1) 

    # Evaluate B0 and B1
    B0 = sp.lil_matrix((nbfunct, nbx))
    B1 = sp.lil_matrix((nbfunct, nbx))

    for _ in range(len(knots)):
        # Get knot
        knot = knots[_]    
    
        # Find knot-span
        knot_span = helpers.find_span_linear(degree, knotvector, nbfunct, knot)
        
        # Find element
        element = np.where(table_functions_element[:, 0] == knot_span)[0].tolist()
        
        # Find functions at the element
        functions_element = table_functions_element[element, 1:][0]

        # Evaluate B0 and B1 at the knot
        B0t, B1t = helpers.basis_function_ders(degree, knotvector, knot_span, knot, 1)

        # Set procedure if knot is in the knot-vector
        if knot in np.unique(knotvector)[1:-1]:             
            # Erase zeros
            B0t = B0t[:-multiplicity] 
            B1t = B1t[:-multiplicity] 

            # Erase zeros functions
            functions_element = functions_element[:-multiplicity]

        # Replace values
        B0[np.ix_(functions_element, [_])] = np.asarray(B0t).reshape((-1,1))
        B1[np.ix_(functions_element, [_])] = np.asarray(B1t).reshape((-1,1))

    return B0, B1

def wq_get_shape_B(degree, nbel, r, multiplicity= 1, maxrule= 1): 
    " Return the shape of basis in WQ approach "

    # Set number of functions 
    nbfunct = degree + multiplicity*(nbel - 1) + 1

    # Set number of quadrature points WQ:
    nbx = 2*(degree + r) + nbel*(maxrule + 1) - 2*maxrule - 3

    # Set table of positions of quadrature points 
    tableOfPointsOnSpan = np.zeros((nbel, 2), dtype= int)
    tableOfPointsOnSpan[0, 0] = 0; tableOfPointsOnSpan[0, 1] = degree + r -1 
    tableOfPointsOnSpan[-1, 0] = nbx - (degree + r) ; tableOfPointsOnSpan[-1,1] = nbx - 1
     
    for i in range(1, nbel - 1):
        tableOfPointsOnSpan[i, 0] = tableOfPointsOnSpan[i - 1, 1]
        tableOfPointsOnSpan[i, 1] = tableOfPointsOnSpan[i, 0] + 1 + maxrule
    
    # Set table of functions on every element 
    tableOfFunctionsOnElement = np.zeros((nbel, degree + 1), dtype= int); 
    tableOfFunctionsOnElement[0, :] = np.arange(degree + 1) 

    for _ in range(1, nbel): 
        # Set values of the table
        tableOfFunctionsOnElement[_, 0] = tableOfFunctionsOnElement[_-1, 0] + multiplicity
        tableOfFunctionsOnElement[_, 1:] = tableOfFunctionsOnElement[_, 0] + np.arange(1, degree + 1) 

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

def wq_find_positions(degree, knotvector, r, maxrule= 1):
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
            xt = np.linspace(kvUnique[i], kvUnique[i+1], 2 + maxrule)[1:-1]
        
        # Assemble vector
        xwq = np.append(xwq, xt)
    
    #  Include knots of knot-vector
    xwq = np.append(xwq, kvUnique)

    # Sort vector
    xwq = np.sort(xwq).tolist()

    return xwq

def wq_find_weights(degree, knotvector, r, maxrule= 1): 

    # Set number of elements
    nbel = len(np.unique(knotvector))-1
    # ------------------------------------
    # Degree p and continuity p-1
    # ------------------------------------
    # Find positions and weights in IGA approach 
    xcgg, wcgg = iga_find_positions_weights(degree, knotvector)

    # Find basis function values at Gauss points 
    B0cgg_p0, B1cgg_p0 = eval_basis(degree, knotvector, xcgg)

    # Find positions in WQ approach
    xwq = wq_find_positions(degree, knotvector, r, maxrule= maxrule)

    # ------------------------------------
    # Degree p and continuity 2
    # ------------------------------------
    # Create curve
    knotvector_p1 = create_knotvector(degree, nbel, multiplicity= 2)

    # Find basis function values at Gauss points
    B0cgg_p1 = eval_basis(degree, knotvector_p1, xcgg, multiplicity= 2)[0]

    # Find basis function values at WQ points
    B0wq_p1 = eval_basis(degree, knotvector_p1, xwq, multiplicity= 2)[0]  

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
        B0, B1 = eval_basis(degree, knotvector, xwq)
        W0, W1 = wq_find_weights(degree, knotvector, r, maxrule= maxrule)
    else : 
        # Create model
        nbel_model = degree + 3
        knotvector_model = create_knotvector(degree, nbel_model)
        xwq_model = wq_find_positions(degree, knotvector_model, r, maxrule= maxrule)
        B0_model, B1_model = eval_basis(degree, knotvector_model, xwq_model)
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
