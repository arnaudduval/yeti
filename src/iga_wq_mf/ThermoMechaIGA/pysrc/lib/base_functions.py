"""
.. This module contains functions that will be used in other modules. 
.. Some of them use fortran functions to accelerate computing time
.. Joaquin Cornejo
"""

from .__init__ import *

# ==========================
# GENERAL FUNCTIONS
# ==========================

def erase_rows_csr(rows2er, indi_in, indj_in, data_in, isfortran=True):
    " Returns new data after erasing rows in CSR format "
    
    # Initialize 
    indi_int, indj_int = np.copy(indi_in), np.copy(indj_in) 
    if isfortran: indi_int -= 1; indj_int -= 1 
    indi_outt = np.copy(indi_int)
    indj_out, data_out = [], []

    # Delete indices
    indi_outt = np.delete(indi_outt, rows2er)
    
    # Copy column indices
    for _ in range(len(indi_outt)-1): 
        indj_out.extend(indj_int[indi_outt[_]:indi_outt[_+1]])
    indj_out = np.array(indj_out, dtype=int)

    # Copy data 
    for a_in in data_in: 
        a_in = np.atleast_2d(a_in); a_out = []
        for _ in range(len(indi_outt)-1): 
            a_out.extend(a_in[indi_outt[_]:indi_outt[_+1], :])
        a_out = np.array(a_out)
        data_out.append(a_out)

    # Update row indices
    indi_out = [0]
    for _ in range(len(indi_outt)-1): 
        nnz = indi_outt[_+1] - indi_outt[_]
        newvalue = indi_out[-1] + nnz
        indi_out.append(newvalue)
    indi_out = np.array(indi_out, dtype=int)

    if isfortran: indi_out += 1; indj_out += 1 

    return indi_out, indj_out, data_out

def generate_rand_positive_matrix(dim, nnz):
    " Return nnz random symetric positive definite matrix of size (dim x dim) "
    
    # Generate random matrix 
    A = np.random.random((dim, dim, nnz))
    I = np.eye(dim)
    
    for i in range(nnz):
        # Construct symmetric matrix 
        B = A[:, :, i]
        B = np.matmul(B, B.T) + I
        A[:, :, i] = B

    return A 

def compute_jacobien_mean(J):
    """ Returns the average of the diagonal of H. 
        where H is the stretch tensor in the polar decomposition of J.
    """

    # Get the size of J
    J = np.atleast_3d(J)
    nnz = np.shape(J)[2]
    dim = np.shape(J)[1] # or [0]

    # Initialize diagonal
    diag = np.zeros(dim)

    for i in range(nnz):
        _, P = sclin.polar(J[:, :, i])
        diag += np.diagonal(P)/nnz

    return diag

# ==========================
# B-SPLINE FUNCTIONS
# ==========================

def create_knotvector(p, nbel, multiplicity=1):
    " Creates an uniform, open, with maximum regularity, knot-vector "

    # Set knot-vector to be inserted
    kv_unique = np.linspace(0., 1., nbel + 1)[1 : -1]

    # Create knot-vector 
    knotvector = []
    for _ in range(p+1): 
        knotvector.append(0.0)

    for knot in kv_unique: 
        for _ in range(multiplicity): 
            knotvector.append(knot)

    for _ in range(p+1): 
        knotvector.append(1.0)

    # Change type
    knotvector = np.array(knotvector)
    
    return knotvector

def eval_basis_python(degree, knotvector, knots, multiplicity=1): 
    """ Evaluates B-spline functions at given knots. 
        Knot-vector needs to be regular
    """

    # Set number of points x
    nbx = len(knots)

    # Set number of elements 
    nbel = len(np.unique(knotvector)) - 1

    # Find number of control points 
    nb_ctrlpts = degree + multiplicity*(nbel - 1) + 1

    # Set table of functions per element 
    table_functions_element = np.zeros((nbel, degree + 2), dtype=int); 
    table_functions_element[0, 0] = degree; table_functions_element[0, 1:] = np.arange(degree + 1) 

    for _ in range(1, nbel): 
        # Set values of the table
        table_functions_element[_, :2] = table_functions_element[_-1, :2] + multiplicity
        table_functions_element[_, 2:] = table_functions_element[_, 1] + np.arange(1, degree + 1) 

    # Evaluate B0 and B1
    B0 = sp.lil_matrix((nb_ctrlpts, nbx))
    B1 = sp.lil_matrix((nb_ctrlpts, nbx))

    for i, knot in enumerate(knots):
    
        # Find knot-span
        knot_span = helpers.find_span_linear(degree, knotvector, nb_ctrlpts, knot)
        
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
        B0[np.ix_(functions_element, [i])] = np.asarray(B0t).reshape((-1,1))
        B1[np.ix_(functions_element, [i])] = np.asarray(B1t).reshape((-1,1))

    # Convert COO to CSR format 
    B0, B1 = B0.tocsr(), B1.tocsr()

    return B0, B1

def eval_basis_fortran(degree, knotvector, knots):
    " Evaluates B-spline functions at given knots using fortran libraries "

    B, indi, indj = basis_weights.get_basis_generalized_csr(
                                    degree, knotvector, knots)
                            
    return B, indi, indj

# ==========================
# IGA FUNCTIONS
# ==========================

def gaussTable(order):
    " Computes Gauss weights and positions in isoparametric space for a given degree "

    if order == 1:
        pos = [0.0]
        wgt = [2.0]
    elif order == 2:
        pos = [-0.577350269189625764509148780502,
                0.577350269189625764509148780502]

        wgt = [1.0,
                1.0]
    elif order == 3:
        pos = [-0.774596669241483377035853079956,
                0.0,
                0.774596669241483377035853079956]

        wgt = [5.0 / 9.0,
                8.0 / 9.0,
                5.0 / 9.0]
    elif order == 4:
        pos = [- 0.861136311594052575223946488893,
                - 0.339981043584856264802665759103,
                0.339981043584856264802665759103, 
                0.861136311594052575223946488893]

        wgt = [0.347854845137453857373063949222,
                0.652145154862546142626936050778,
                0.652145154862546142626936050778,
                0.347854845137453857373063949222]
    elif order == 5:
        pos = [- 0.906179845938663992797626878299,
                - 0.538469310105683091036314420700,
                0.0,
                0.538469310105683091036314420700,
                0.906179845938663992797626878299]

        wgt = [0.236926885056189087514264040720,
                0.478628670499366468041291514836,
                0.568888888888888888888888888889,
                0.478628670499366468041291514836,
                0.236926885056189087514264040720]
    elif order == 6:
        pos = [- 0.932469514203152027812301554494,
                - 0.661209386466264513661399595020,
                - 0.238619186083196908630501721681,
                0.238619186083196908630501721681,
                0.661209386466264513661399595020,
                0.932469514203152027812301554494]

        wgt = [0.171324492379170345040296142173,
                0.360761573048138607569833513838, 
                0.467913934572691047389870343990,
                0.467913934572691047389870343990,
                0.360761573048138607569833513838,
                0.171324492379170345040296142173]
    elif order == 7:
        pos = [- 0.949107912342758524526189684048,
                - 0.741531185599394439863864773281,
                - 0.405845151377397166906606412077,
                0.0,
                0.405845151377397166906606412077,
                0.741531185599394439863864773281,
                0.949107912342758524526189684048]

        wgt = [0.129484966168869693270611432679,
                0.279705391489276667901467771424,
                0.381830050505118944950369775489,
                0.417959183673469387755102040816,
                0.381830050505118944950369775489,
                0.279705391489276667901467771424,
                0.129484966168869693270611432679]
    elif order == 8:
        pos = [- 0.960289856497536231683560868569,
                - 0.796666477413626739591553936476,
                - 0.525532409916328985817739049189,
                - 0.183434642495649804939476142360,
                0.183434642495649804939476142360,
                0.525532409916328985817739049189,
                0.796666477413626739591553936476,
                0.960289856497536231683560868569]

        wgt = [0.101228536290376259152531354310,
                0.222381034453374470544355994426,
                0.313706645877887287337962201987,
                0.362683783378361982965150449277,
                0.362683783378361982965150449277,
                0.313706645877887287337962201987,
                0.222381034453374470544355994426,
                0.101228536290376259152531354310]
    elif order == 9:
        pos = [- 0.968160239507626089835576202904,
                - 0.836031107326635794299429788070,
                - 0.613371432700590397308702039341,
                - 0.324253423403808929038538014643,
                0.0,
                0.324253423403808929038538014643,
                0.613371432700590397308702039341,
                0.836031107326635794299429788070,
                0.968160239507626089835576202904]

        wgt = [0.812743883615744119718921581105E-01,
                0.180648160694857404058472031243,
                0.260610696402935462318742869419,
                0.312347077040002840068630406584,
                0.330239355001259763164525069287,
                0.312347077040002840068630406584,
                0.260610696402935462318742869419,
                0.180648160694857404058472031243,
                0.812743883615744119718921581105E-01]
    elif order == 10:
        pos = [-0.9739065285171717,
                -0.8650633666889845,
                -0.6794095682990244,
                -0.4333953941292472,
                -0.1488743389816312,
                0.1488743389816312,
                0.4333953941292472,
                0.6794095682990244,
                0.8650633666889845,
                0.9739065285171717]

        wgt = [0.0666713443086881,
                0.1494513491505806,
                0.2190863625159820,
                0.2692667193099963,
                0.2955242247147529,
                0.2955242247147529,
                0.2692667193099963,
                0.2190863625159820,
                0.1494513491505806,
                0.0666713443086881]                
    elif order == 11:
        pos = [-0.9782286581460570,
                -0.8870625997680953,
                -0.7301520055740494,
                -0.5190961292068118,
                -0.2695431559523450,
                0.0,
                0.2695431559523450,
                0.5190961292068118,
                0.7301520055740494,
                0.8870625997680953,
                0.9782286581460570]

        wgt = [0.0556685671161737,
                0.1255803694649046,
                0.1862902109277343,
                0.2331937645919905,
                0.2628045445102467,
                0.2729250867779006,
                0.2628045445102467,
                0.2331937645919905,
                0.1862902109277343,
                0.1255803694649046,
                0.0556685671161737]
    
    # Change type of arrays
    pos, wgt = np.array(pos), np.array(wgt)

    return pos, wgt

def iga_find_positions_weights(degree, knotvector):
    " Computes Gauss positions and weights in parametric space for a given knotvector "

    def iga_find_positions_weights_element(degree, knotvector, el):
        " Computes Gauss weights and positions in parametric space for a given element "

        # Find knots of the knot-vector
        knots = np.unique(knotvector)

        # Find position and weight of Gauss points
        xg, wg = gaussTable(degree + 1)
        
        # Find position of quadrature points at the element (scaling)
        xg = 0.5*((knots[el+1] - knots[el])*xg + knots[el] + knots[el+1])

        # Find weight of quadrature points at the element (scaling)
        wg = 0.5*(knots[el+1] - knots[el])*wg

        return xg, wg

    # Find number of elements
    nbel = len(np.unique(knotvector)) - 1

    # Initialize
    xg = []; wg = []

    for _ in range(nbel): 
        # Find quadrature points for all elements
        xg_el, wg_el = iga_find_positions_weights_element(degree, knotvector, _)

        # Save data
        xg.extend(xg_el); wg.extend(wg_el)

    # Change type
    xg, wg = np.array(xg), np.array(wg)

    return xg, wg

def iga_find_basis_weights_opt(degree, knotvector):
    " Computes basis and weights in IGA approach "

    # Find positions and weights
    qp_position, qp_weight = iga_find_positions_weights(degree, knotvector)

    # Find basis
    B0, B1 = eval_basis_python(degree, knotvector, qp_position)

    return qp_position, B0, B1, qp_weight

def iga_find_basis_weights_fortran(degree, knotvector): 
    " Computes basis and weights in IGA approach using fortran "

    # Set number of elements
    nbel = len(np.unique(knotvector)) - 1

    # Set number of quadrature points
    nb_qp = (degree + 1) * nbel

    # Set guessed size of data 
    nnz_B = (degree + 1) * nb_qp

    # Get basis and weights from fortran
    qp_position, qp_weight, basis, \
    indi, indj, nnz_I = basis_weights.iga_get_data_csr(degree, knotvector, nb_qp, nnz_B)
 
    return nnz_I, qp_position, qp_weight, basis, indi, indj

# =========================
# WQ FUNCTIONS
# =========================

def wq_get_shape_B(degree, nbel, r, multiplicity=1, maxrule=1): 
    """ Return the shape of basis in WQ approach 
        It considers an open and uniform knot-vector
    """

    # Set number of control points 
    nb_ctrlpts = degree + multiplicity*(nbel - 1) + 1

    # Set number of quadrature points
    nb_qp = 2*(degree + r) + nbel*(maxrule + 1) - 2*maxrule - 3

    # Set table of positions of quadrature points 
    tableOfPointsOnSpan = np.zeros((nbel, 2), dtype=int)
    tableOfPointsOnSpan[0, 0] = 0; tableOfPointsOnSpan[0, 1] = degree + r -1 
    tableOfPointsOnSpan[-1, 0] = nb_qp - (degree + r) ; tableOfPointsOnSpan[-1,1] = nb_qp - 1

    for i in range(1, nbel-1):
        tableOfPointsOnSpan[i, 0] = tableOfPointsOnSpan[i - 1, 1]
        tableOfPointsOnSpan[i, 1] = tableOfPointsOnSpan[i, 0] + 1 + maxrule
    
    # Set table of functions on every element 
    tableOfFunctionsOnElement = np.zeros((nbel, degree + 1), dtype=int); 
    tableOfFunctionsOnElement[0, :] = np.arange(degree + 1) 

    for _ in range(1, nbel): 
        # Set values of the table
        tableOfFunctionsOnElement[_, 0] = tableOfFunctionsOnElement[_-1, 0] + multiplicity
        tableOfFunctionsOnElement[_, 1:] = tableOfFunctionsOnElement[_, 0] + np.arange(1, degree + 1) 

    # Set table of knot-span for every function
    tableOfSpansForFunction = []
    for i in range(nb_ctrlpts):
        indi, _ = np.asarray(tableOfFunctionsOnElement == i).nonzero()
        indi = np.sort(indi)
        tableOfSpansForFunction.append(indi)

    # Find basis shape
    indi_B0, indj_B0 = [], []
    indi_B1, indj_B1 = [], []
    for i in range(nb_ctrlpts) : 
        spanOfFunction = tableOfSpansForFunction[i]
        spanLeft = np.min(spanOfFunction)
        spanRight = np.max(spanOfFunction)

        support = np.arange(tableOfPointsOnSpan[spanLeft,0], \
                            tableOfPointsOnSpan[spanRight,1] + 1, dtype=int)
        
        support_B0 = support
        if i == 0: 
            support_B0 = np.delete(support, [-1])
        elif i == nb_ctrlpts - 1:
            support_B0 = np.delete(support, [0])
        else : 
            support_B0 = np.delete(support, [0 , -1])
        
        indi_B0.extend(i*np.ones(len(support_B0), dtype= int))
        indj_B0.extend(support_B0)

        support_B1 = support
        if i == 0 or i == 1: 
            support_B1 = np.delete(support, [-1])
        elif i == nb_ctrlpts - 1 or i == nb_ctrlpts - 2:
            support_B1 = np.delete(support, [0])
        else : 
            support_B1 = np.delete(support, [0 , -1])
        
        indi_B1.extend(i*np.ones(len(support_B1), dtype= int))
        indj_B1.extend(support_B1)
    
    # Set shape of B0 and B1
    data_B0 = np.ones(len(indi_B0))
    B0shape = sp.csr_matrix((data_B0, (indi_B0, indj_B0)), shape=(nb_ctrlpts, nb_qp))
    
    data_B1 = np.ones(len(indi_B1))
    B1shape = sp.csr_matrix((data_B1, (indi_B1, indj_B1)), shape=(nb_ctrlpts, nb_qp))
    
    return B0shape, B1shape

def wq_solve_equation_system(B, I):
    """ Solves linear system B w = I where B and I are written in CSR format. 
        System can be under, well, over determined
    """

    # Convert type of B and I to array
    B = B.toarray()
    I = I.toarray()

    # Solve system 
    sol1, _, _, _ = np.linalg.lstsq(B, I, rcond=None)

    # Save solution
    w = sol1.reshape((1, -1)).tolist()[0]
    
    return w

def wq_find_positions(degree, knotvector, r, maxrule= 1):
    " Return position of quadrature points in WQ approach "
    
    # Set unique knots
    kv_unique = np.unique(knotvector)

    # Set number of elements
    nbel = len(kv_unique) - 1

    # Initialize 
    qp_position = []
    
    for i in range(nbel):

        if i == 0 or i == nbel - 1:
            # Case: boundary knot-spans
            xt = np.linspace(kv_unique[i], kv_unique[i+1], degree+r)[1:-1]
        else:
            # Case: interior knot-spans
            xt = np.linspace(kv_unique[i], kv_unique[i+1], 2 + maxrule)[1:-1]
        
        # Save
        qp_position.extend(xt)
    
    # Include knots of knot-vector
    qp_position.extend(list(kv_unique))

    # Change type and sort
    qp_position = np.array(qp_position)
    qp_position.sort()
    
    return qp_position

def wq_find_weights(degree, knotvector, r):
    " Returns weights at quadrature points in WQ approach using space S^[p-1] method "
    
    # Set number of elements
    nbel = len(np.unique(knotvector))-1

    # ------------
    # Space S^p_r
    # ------------
    # Find positions and weights in IGA approach 
    xcgg, wcgg = iga_find_positions_weights(degree, knotvector)

    # Find basis function values at Gauss points 
    B0cgg_p0, B1cgg_p0 = eval_basis_python(degree, knotvector, xcgg)

    # Find positions in WQ approach
    xwq = wq_find_positions(degree, knotvector, r)

    # Find basis function values at WQ points
    B0wq_p0 = eval_basis_python(degree, knotvector, xwq)[0]

    # ---------------------
    # Sapce S^[p-1]_[r-1]
    # ---------------------
    # Create knot-vector
    knotvector_p1 = knotvector[1:-1]

    # Find basis function values at Gauss points
    B0cgg_p1 = eval_basis_python(degree-1, knotvector_p1, xcgg)[0]

    # Find basis function values at WQ points
    B0wq_p1 = eval_basis_python(degree-1, knotvector_p1, xwq)[0]  

    # ------------------
    # Compute Integrals
    # ------------------
    # Compute exact integral I00 = int Bi(p) Bj(p) dx to calculate W00
    I00 = B0cgg_p0 @ sp.diags(wcgg) @ B0cgg_p0.T
    
    # Compute exact integral I01 = int Bi(p) Bj'(p) dx to calculate W01
    I01 = B0cgg_p1 @ sp.diags(wcgg) @ B0cgg_p0.T 

    # Compute exact integral I01 = int Bi'(p) Bj(p) dx to calculate W01
    I10 = B0cgg_p0 @ sp.diags(wcgg) @ B1cgg_p0.T
    
    # Compute exact integral I11 = int Bi'(p) Bj'(p) dx to calculate W11
    I11 = B0cgg_p1 @ sp.diags(wcgg) @ B1cgg_p0.T

    # Set number of control points
    nb_ctrlpts = degree + nbel

    # Set number of WQ quadrature points 
    nb_qp_wq = len(xwq)

    # Initialize data for W00, W01, W10 and W11
    data_W00, data_W01, data_W10, data_W11 = [], [], [], []
    
    # Find shape of B and I
    shape_B0_p0, shape_B1_p0 = wq_get_shape_B(degree, nbel, r)
    shape_B_p1 = wq_get_shape_B(degree-1, nbel, r+1)[0]

    # -------------------
    # Computation of W00
    # -------------------
    Ishape = shape_B0_p0 @ shape_B0_p0.T 
    indi = []; indj = []
    for i in range(nb_ctrlpts):
        F = np.nonzero(Ishape[:, i])[0]
        P = np.nonzero(shape_B0_p0[i, :])[1]

        Bt = B0wq_p0[np.ix_(F,P)]
        It = I00[np.ix_(F, [i])]

        data_W00.extend(wq_solve_equation_system(Bt, It))
        indi.extend(i*np.ones(len(P), dtype= int))
        indj.extend(P)
    W00 = sp.csr_matrix((data_W00, (indi, indj)), shape=(nb_ctrlpts, nb_qp_wq))

    # --------------------
    # Computation of W01
    # --------------------
    Ishape = shape_B_p1 @ shape_B0_p0.T 
    indi = []; indj = []
    for i in range(nb_ctrlpts):
        F = np.nonzero(Ishape[:,i])[0]
        P = np.nonzero(shape_B0_p0[i, :])[1]

        Bt = B0wq_p1[np.ix_(F,P)]
        It = I01[np.ix_(F, [i])]
        
        data_W01.extend(wq_solve_equation_system(Bt, It))
        indi.extend(i*np.ones(len(P), dtype= int))
        indj.extend(P)
    W01 = sp.csr_matrix((data_W01, (indi, indj)), shape=(nb_ctrlpts, nb_qp_wq))
    
    # --------------------
    # Computation of W10
    # --------------------
    Ishape = shape_B1_p0 @ shape_B0_p0.T 
    indi = []; indj = []
    for i in range(nb_ctrlpts):
        F = np.nonzero(Ishape[:, i])[0]
        P = np.nonzero(shape_B1_p0[i, :])[1]

        Bt = B0wq_p0[np.ix_(F,P)]
        It = I10[np.ix_(F, [i])]

        data_W10.extend(wq_solve_equation_system(Bt, It))
        indi.extend(i*np.ones(len(P), dtype= int))
        indj.extend(P)

    W10 = sp.csr_matrix((data_W10, (indi, indj)), shape=(nb_ctrlpts, nb_qp_wq))

    # -------------------
    # Computation of W11
    # -------------------
    Ishape = shape_B_p1 @ shape_B1_p0.T
    indi = []; indj = []
    for i in range(nb_ctrlpts):
        F = np.nonzero(Ishape[:,i])[0]
        P = np.nonzero(shape_B1_p0[i, :])[1]

        Bt = B0wq_p1[np.ix_(F,P)]
        It = I11[np.ix_(F, [i])]

        data_W11.extend(wq_solve_equation_system(Bt, It))
        indi.extend(i*np.ones(len(P), dtype= int))
        indj.extend(P)
    W11 = sp.csr_matrix((data_W11, (indi, indj)), shape=(nb_ctrlpts, nb_qp_wq))
    
    return W00, W01, W10, W11 

def wq_find_basis_weights_opt(degree, knotvector, r): 
    " Return basis and weights in WQ approach computed in an efficient way "
    
    # Set number of elements
    nbel = len(np.unique(knotvector)) - 1

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
        W00_model *= nbel_model/nbel
        W01_model *= nbel_model/nbel
        B1_model *= nbel/nbel_model

        #------------------------------
        # Algorithm to get the results
        #------------------------------
        # Set number of points
        nbx = 2*(nbel + degree + r) - 5

        # Set number of functions
        nbfunct = degree + nbel

        # Initialize 
        data_B0, data_B1, data_W00, data_W01, data_W10, data_W11 = [], [], [], [], [], []
        function_B0, support_B0, function_B1, support_B1 = [], [], [], []

        # Find shape of B
        shape_B0, shape_B1 = wq_get_shape_B(degree, nbel, r)

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

def wq_find_basis_weights_fortran(degree, knotvector): 
    " Computes basis and weights in WQ approach using Fortran "

    # Set guessed size of data 
    nnz_B, nb_qp = basis_weights.wq_get_size_data(degree, knotvector)

    # Get basis and weights from fortran
    qp_position, basis, weights, indi, indj, nnz_I = basis_weights.wq_get_data_csr(
                                                    degree, knotvector, nnz_B, nb_qp)

    return nnz_I, qp_position, basis, weights, indi, indj
    
# =========================
# MF FUNCTIONS
# =========================

def solver_scipy(A, b, nbIterations=100, epsilon=1e-10, PreCond='ilu', isCG=True):
    """ Solve system using iterative method : conjugate gradient or bi-conjugate gradient. 
        It can use ILU preconditioner to accelerate convergence. 
    """

    # Select preconditionner (by the moment only ILU)
    if PreCond == 'ilu': 
        B = sclin.spilu(A)
        Mx = lambda x: B.solve(x)
        M = sclin.LinearOperator(A.shape, Mx)

    # Solve with iterative method
    if isCG: x, info = sclin.cg(A, b, tol=epsilon, maxiter=nbIterations, M=M)
    else: x, info = sclin.bicgstab(A, b, tol=epsilon, maxiter=nbIterations, M=M)

    return x

def eigen_decomposition(indi, indj, data, robin_condition=[0, 0], coefs=None):
    """ 
        Eigen decomposition generalized KU = MUD
        K: stiffness matrix, K = int B1 B1 dx = W11 * B1
        M: mass matrix, M = int B0 B0 dx = W00 * B0
        U: eigenvectors matrix
        D: diagonal of eigenvalues
    """

    # Get number of quadrature points
    nb_qp = np.max(indj)

    # Unpack B0, B1 and W00, W11
    [B0, B1, W0, W1] = data

    # Create mcoef and kcoef is are none type
    if coefs is None: mcoefs = np.ones(nb_qp); kcoefs = np.ones(nb_qp)
    else: [mcoefs, kcoefs] = coefs

    # Compute eigen values and vectors
    eigenvalues, eigenvectors = solver.eigen_decomposition_py(indi, indj, B0, W0, 
                                        B1, W1, mcoefs, kcoefs, robin_condition)

    return eigenvalues, eigenvectors

def tensor_decomposition_3D(n_list, CC_matrix, dim=3):
    """ Tensor decomposition of CC matrix to improve Fast diagonalization precontionner
        Based on "Preconditioners for Isogemetric Analysis" by M. Montardini
    """
    if dim != 3: raise Warning('Method not implemented')

    # Unpack shape
    [n1, n2, n3] = n_list
    try: n1 = n1[0]; n2 = n2[0]; n3 = n3[0]
    except: pass

    # Get diagonal of CC matrix
    coefs = np.zeros((dim, n1*n2*n3))
    for i in range(dim): coefs[i, :] = CC_matrix[i, i, :]

    # Initialize
    u1, u2, u3 = np.ones(n1), np.ones(n2), np.ones(n3)
    w1, w2, w3 = np.ones(n1), np.ones(n2), np.ones(n3)
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
                for j in range(n3):
                    m = Vscript[:, :, j].min()
                    M = Vscript[:, :, j].max()
                    w3[j] = np.sqrt(m*M)

        for k in range(dim):
            c = -1
            # Compute Wscript temporary
            for l in range(dim):
                if k != l:
                    c += 1
                    for i3 in range(n3):
                        for i2 in range(n2):
                            for i1 in range(n1):
                                U = [u1[i1], u2[i2], u3[i3]]
                                W = [w1[i1], w2[i2], w3[i3]]
                                Wscript[c, i1, i2, i3] = coefstens[k, i1, i2, i3]*U[l]*U[k]\
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
                for j in range(n3):
                    m = Nscript[:, :, j].min()
                    M = Mscript[:, :, j].max()
                    u3[j] = np.sqrt(m*M) 

    return u1, u2, u3, w1, w2, w3

def compute_eig_diag(eig_u, eig_v, eig_w, coefs=[1, 1, 1]):
    " Computes diagonal of eigen values "

    def kron3vec(A, B, C):
        " Computes kron product of 3 vectors: A x B x C"
        result = np.kron(np.kron(A, B), C)
        return result

    # Initialize
    one_u = np.ones(len(eig_u))
    one_v = np.ones(len(eig_v))
    one_w = np.ones(len(eig_w))

    # Compute eigen diagonal
    eig_diag = coefs[0]*kron3vec(one_w, one_v, eig_u)
    eig_diag += coefs[1]*kron3vec(one_w, eig_v, one_u)
    eig_diag += coefs[2]*kron3vec(eig_w, one_v, one_u)
    
    return eig_diag

def fast_diagonalization(U, V, W, D, array_in, fdtype='steady'):
    " Compute fast diagonalization using Fortran"
    if fdtype == 'steady':
        array_out = solver.fd_steady_heat_3d(U, V, W, D, array_in)
    elif fdtype == 'elastic':
        array_out = solver.fd_elasticity_3d(U, V, W, D, array_in)
    return array_out

# =========================
# USING TENSOR ALGEBRA
# =========================

def tensor_n_mode_product(x, cols, Ucsr=None, mode=1, isdense=False, istransp=False): 
    """ Computes tensor n-mode product with a matrix R = X xn U,
        where X = vec(x). It uses fortran functions. 
        WARNING : It has not been tested.
    """
    # Verify if input are well defined
    nc_total = np.prod(cols)
    if nc_total != len(x): raise Warning('Shape of X is not well defined') 
    if mode not in [1, 2, 3]: raise Warning('n-mode product can only be 1, 2 or 3')

    # Unpack data of matrix
    [indi, indj, data] = Ucsr

    # Get number of rows of U
    nrU = len(indi) - 1 

    # Get size of R
    rows = cols; rows[mode-1] = nrU
    nrR = np.prod(rows)

    # Compute tensor n-mode product
    R = assembly.tensor_n_mode_product_py(cols, x, data, indi, indj, mode, isdense, istransp, nrR)
    
    return R
