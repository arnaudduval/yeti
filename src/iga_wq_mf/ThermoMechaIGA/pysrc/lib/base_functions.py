"""
.. This module contains basis functions which use python or fortran
.. Joaquin Cornejo
"""

# Python libraries
import numpy as np
from geomdl import helpers
from scipy import sparse as sp

# My libraries
from iga_wq_mf import basis_weights

# ==========================
# GENERAL FUNCTIONS
# ==========================

def erase_rows_csr(rows2er, indi_in, indj_in, data_in, isfortran=True):
    "Returns new data after erasing some rows in CSR format"
    
    # Initialize outputs
    indi_int, indj_int = np.copy(indi_in), np.copy(indj_in) 
    if isfortran: 
        indi_int -= 1
        indj_int -= 1 

    indi_outt = np.copy(indi_int)
    indj_out = []
    data_out = []

    # Delete indices
    indi_outt = np.delete(indi_outt, rows2er)
    
    # Copy some indices
    for _ in range(len(indi_outt)-1): 
        indj_out.extend(indj_int[indi_outt[_]:indi_outt[_+1]])
    indj_out = np.array(indj_out, dtype=int)

    # Copy some values
    for a_in in data_in: 
        a_out = []
        for _ in range(len(indi_outt)-1): 
            a_out.extend(a_in[indi_outt[_]:indi_outt[_+1]])
        a_out = np.array(a_out)
        data_out.append(a_out)

    # Update indi output
    indi_out = [0]
    for _ in range(len(indi_outt)-1): 
        nnz = indi_outt[_+1] - indi_outt[_]
        newvalue = indi_out[-1] + nnz
        indi_out.append(newvalue)
    indi_out = np.array(indi_out, dtype=int)

    if isfortran: 
        indi_out += 1
        indj_out += 1 

    return indi_out, indj_out, data_out

# ==========================
# B-SPLINE FUNCTIONS
# ==========================

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

def eval_basis_python(degree, knotvector, knots, multiplicity= 1): 
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

def eval_basis_fortran(degree, nbel, knots, multiplicity= 1):
    " Evaluates B-spline functions at given knots using fortran libraries "

    B0, B1, indi, indj = basis_weights.get_basis_generalized_csr(
                            degree, nbel, len(knots), knots, multiplicity)

    return B0, B1, indi, indj

# ==========================
# IGA FUNCTIONS
# ==========================

def gaussTable(pgaus):
    " Computes Gauss weights and positions in isogeometric space for a known degree "

    if pgaus == 1:
        x = [0.0]
        w = [2.0]
    elif pgaus == 2:
        x = [-0.577350269189625764509148780502, \
              0.577350269189625764509148780502]
        # ------------
        w = [1.0,\
             1.0]
    elif pgaus == 3:
        x = [-0.774596669241483377035853079956, \
              0.0, \
              0.774596669241483377035853079956]
        # ------------
        w = [5.0 / 9.0, \
             8.0 / 9.0, \
             5.0 / 9.0]
    elif pgaus == 4:
        x = [- 0.861136311594052575223946488893, \
             - 0.339981043584856264802665759103, \
               0.339981043584856264802665759103, \
               0.861136311594052575223946488893]
        # ------------
        w = [0.347854845137453857373063949222, \
             0.652145154862546142626936050778, \
             0.652145154862546142626936050778, \
             0.347854845137453857373063949222]
    elif pgaus == 5:
        x = [- 0.906179845938663992797626878299, \
             - 0.538469310105683091036314420700, \
               0.0, \
               0.538469310105683091036314420700, \
               0.906179845938663992797626878299]
        # -----------
        w = [0.236926885056189087514264040720, \
             0.478628670499366468041291514836, \
             0.568888888888888888888888888889, \
             0.478628670499366468041291514836, \
             0.236926885056189087514264040720]
    elif pgaus == 6:
        x = [- 0.932469514203152027812301554494, \
             - 0.661209386466264513661399595020, \
             - 0.238619186083196908630501721681, \
               0.238619186083196908630501721681, \
               0.661209386466264513661399595020, \
               0.932469514203152027812301554494]
        # -----------
        w = [0.171324492379170345040296142173, \
             0.360761573048138607569833513838, \
             0.467913934572691047389870343990, \
             0.467913934572691047389870343990, \
             0.360761573048138607569833513838, \
             0.171324492379170345040296142173]
    elif pgaus == 7:
        x = [- 0.949107912342758524526189684048, \
             - 0.741531185599394439863864773281, \
             - 0.405845151377397166906606412077, \
               0.0, \
               0.405845151377397166906606412077, \
               0.741531185599394439863864773281, \
               0.949107912342758524526189684048]
        # -----------
        w = [0.129484966168869693270611432679, \
             0.279705391489276667901467771424, \
             0.381830050505118944950369775489, \
             0.417959183673469387755102040816, \
             0.381830050505118944950369775489, \
             0.279705391489276667901467771424, \
             0.129484966168869693270611432679]
    elif pgaus == 8:
        x = [- 0.960289856497536231683560868569, \
             - 0.796666477413626739591553936476, \
             - 0.525532409916328985817739049189, \
             - 0.183434642495649804939476142360, \
               0.183434642495649804939476142360, \
               0.525532409916328985817739049189, \
               0.796666477413626739591553936476, \
               0.960289856497536231683560868569]
        # -----------
        w = [0.101228536290376259152531354310, \
             0.222381034453374470544355994426, \
             0.313706645877887287337962201987, \
             0.362683783378361982965150449277, \
             0.362683783378361982965150449277, \
             0.313706645877887287337962201987, \
             0.222381034453374470544355994426, \
             0.101228536290376259152531354310]
    elif pgaus == 9:
        x = [- 0.968160239507626089835576202904, \
             - 0.836031107326635794299429788070, \
             - 0.613371432700590397308702039341, \
             - 0.324253423403808929038538014643, \
               0.0, \
               0.324253423403808929038538014643, \
               0.613371432700590397308702039341, \
               0.836031107326635794299429788070, \
               0.968160239507626089835576202904]
        # -----------
        w = [0.812743883615744119718921581105E-01, \
             0.180648160694857404058472031243, \
             0.260610696402935462318742869419, \
             0.312347077040002840068630406584, \
             0.330239355001259763164525069287, \
             0.312347077040002840068630406584, \
             0.260610696402935462318742869419, \
             0.180648160694857404058472031243, \
             0.812743883615744119718921581105E-01]
    
    # Change type of arrays
    x = np.asarray(x); w = np.asarray(w)

    return x, w

def iga_find_positions_weights_element(degree, knotvector, el):
    " Computes Gauss weights and positions in isogeometric space for a known element "

    # Find knots of the knot-vector
    knots = np.unique(knotvector)

    # Find position and weight of Gauss points in isogeometric space
    xg, wg = gaussTable(degree + 1)
    
    # Find position of quadrature points at the element (scaling)
    xg = 0.5*((knots[el+1] - knots[el])*xg + knots[el] + knots[el+1])

    # Find weight of quadrature points at the element (scaling)
    wg = 0.5*(knots[el+1] - knots[el])*wg

    return xg, wg

def iga_find_positions_weights(degree, knotvector):
    " Computes Gauss weights and positions in parametric space for a known curve "

    # Find number of elements
    nbel = len(np.unique(knotvector)) - 1

    # Initialize vectors
    xg = []; wg = []

    for _ in range(nbel): 
        # Find quadrature pints for all elements
        xg_el, wg_el = iga_find_positions_weights_element(degree, knotvector, _)

        # Assmbly vectors
        xg.extend(xg_el); wg.extend(wg_el)

    return xg, wg

# =========================
# WQ FUNCTIONS
# =========================

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
    " Return solution of matrix system B w = I. B and I are in CSR format "

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

# =========================
# MF FUNCTIONS
# =========================

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

# =========================
# USING TENSOR ALGEBRA
# =========================

def get_indices_3D(DI):
    "Returns the indices of I1 x I2 x I3, where x is kronecker product"

    def get_indices_kron_product(indi_A, indj_A, indi_B, indj_B):

        # Defines some values
        nb_rows_A = len(indi_A) - 1
        nb_rows_B = len(indi_B) - 1
        nb_rows_C = nb_rows_A * nb_rows_B
        nb_cols_B = max(indj_B)

        # Initialize
        indi_C = np.zeros(nb_rows_C+1, dtype=int)
        indi_C[0] = 1
        
        # Set indices i in CSR
        for i in range(nb_rows_A):
            for j in range(nb_rows_B):
                k = i*nb_rows_B + j
                nnz_A = indi_A[i+1] - indi_A[i]
                nnz_B = indi_B[j+1] - indi_B[j]
                nnz_C = nnz_A*nnz_B
                indi_C[k+1] = indi_C[k] + nnz_C

        # Set indices j in CSR
        indj_C = []
        for i in range(nb_rows_A):
            for j in range(nb_rows_B):
                k = i*nb_rows_B + j
                
                indj_C_temp = []
                for m in range(indi_A[i], indi_A[i+1]):
                    for n in range(indi_B[j], indi_B[j+1]):
                        indj_C_temp.append(indj_A[m]*nb_cols_B + indj_B[n])

                # Update values
                indj_C.extend(indj_C_temp)
        
        indj_C = np.array(indj_C, dtype=int)

        return indi_C, indj_C

    # We assume DI has 3 dimensions 
    I1, I2, I3 = DI[0], DI[1], DI[2]
    indi_I1, indi_I2, indi_I3 = I1.indptr, I2.indptr,I3.indptr
    indj_I1, indj_I2, indj_I3 = I1.indices, I2.indices, I3.indices
    
    # Find indices of I1 x I2
    indi_temp, indj_temp = get_indices_kron_product(indi_I1, indj_I1, indi_I2, indj_I2)

    # Find indices of I1 x I2 x I3
    indi, indj = get_indices_kron_product(indi_temp, indj_temp, indi_I3, indj_I3)

    return indi, indj

def tensor_n_mode_product(X, U, n): 
    " Evaluates tensor n-mode product with a matrix R = X xn U. We assume 3D geometries or inferior"

    # Find shape of the tensor and the matrix
    I1, I2, I3 = np.shape(X)
    J, In = np.shape(U)

    if n == 0: # In = I1
        R = np.zeros((J, I2, I3))
        for ii3 in range(I3):
            for ii2 in range(I2): 
                for jj in range(J):
                    s = 0.0
                    for ii1 in range(I1): 
                        s += X[ii1, ii2, ii3]*U[jj, ii1] 
                    R[jj, ii2, ii3] = s

    elif n == 1: # In = I2
        R = np.zeros((I1, J, I3))
        for ii3 in range(I3):
            for jj in range(J):
                for ii1 in range(I1): 
                    s = 0.0
                    for ii2 in range(I2): 
                        s += X[ii1, ii2, ii3]*U[jj, ii2] 
                    R[ii1, jj, ii3] = s

    elif n == 2: # In = I3
        R = np.zeros((I1, I2, J))
        for jj in range(J):
            for ii2 in range(I2):
                for ii1 in range(I1): 
                    s = 0.0
                    for ii3 in range(I3): 
                        s += X[ii1, ii2, ii3]*U[jj, ii3] 
                    R[ii1, ii2, jj] = s

    return R

def get_matrix_3D(coef, DB, DW, DI): 
    """Returns non-zero values of matrix given by W.diag(coef).B, where W and B are n-rank tensors.
    DB, DW, DI contains B, W and I matrices for each dimension. 
    The functions is built to 3D cases but it may be generalized for other dimensions
    """

    # Initialize
    nonzerovalues = []
    nbrows = np.ones(3, dtype=int)
    nbcols = np.ones(3, dtype=int)

    for dim in range(3): # 3D figures
        nrow, ncol = np.shape(DB[dim])
        nbrows[dim] = nrow
        nbcols[dim] = ncol

    for i3 in range(nbrows[2]):
        for i2 in range(nbrows[1]): 
            for i1 in range(nbrows[0]):
                ii = [i1, i2, i3]
                nnz_indi = []
                nnz_indj = []

                for dim in range(3):
                    F = np.nonzero(DI[dim][:, ii[dim]])[0]
                    P = np.nonzero(DW[dim][ii[dim], :])[1]
                    nnz_indi.append(F)
                    nnz_indj.append(P)

                # Create tensor from coefficients
                X = np.zeros((len(nnz_indj[0]), len(nnz_indj[1]), len(nnz_indj[2])))
                for j3 in range(len(nnz_indj[2])):
                    for j2 in range(len(nnz_indj[1])):
                        for j1 in range(len(nnz_indj[0])):
                            genPosCoef = (nnz_indj[0][j1] + nnz_indj[1][j2]*nbcols[0]
                                        + nnz_indj[2][j3]*nbcols[0]*nbcols[1])
                            X[j1, j2, j3] = coef[genPosCoef]

                Xnew = X
                for dim in range(2, -1, -1):
                    F = nnz_indi[dim]
                    P = nnz_indj[dim]

                    Bt = DB[dim][np.ix_(F,P)].toarray()
                    Wt = DW[dim][np.ix_([ii[dim]],P)].toarray()
                    Ut = Bt @ np.diag(Wt[0])
                    R = tensor_n_mode_product(Xnew, Ut, dim) 
                    Xnew = R

                for j3 in range(len(nnz_indi[2])):
                    for j2 in range(len(nnz_indi[1])):
                        for j1 in range(len(nnz_indi[0])):
                            nonzerovalues.append(R[j1, j2, j3])

    return nonzerovalues