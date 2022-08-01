"""
.. Integrating elastoplasticity in python
"""

import numpy as np

def block_dot_product(d, A, B):
    """ Computes dot product of A and B. 
    Both are actually vectors arranged following each dimension
    A = [Au, Av, Aw] and B = [Bu, Bv, Bw]. Then A.B = Au.Bu + Av.Bv + Aw.Bw
    """

    # Initialize
    result = 0.0

    for i in range(d):
        result += A[i, :] @ B[i, :]

    return result

def array_maps_tensor(d, p):
    "Converts an array index into a symmetric second-order tensor index "

    if d == 3:
        if p < d:
            i = p; j = p
        else:
            if p == 3:
                i = 0; j = 1
            elif p == 4:
                i = 0; j = 2
            elif p == 5:
                i = 1; j = 2
    elif d == 2:
        if p < d:
            i = p, j = p
        elif p == 2:
            i = 0; j = 1
    else: 
        raise Warning('Not possible')

    return i, j

def array2st(d, array): 
    """Returns second-order tensor from array
    i.e. from [t11, t22, t12], one gets [[t11, t12], [t21, t22]] (with t21 = t12)
    """

    # Initialize
    ddl = int(d*(d+1)/2)
    tensor = np.zeros((d, d))

    for p in range(ddl):
        i, j = array_maps_tensor(d, p)
        tensor[i, j] = array[p]

    # Set diagonal of tensor
    diag = np.diag(np.diag(tensor))

    # Send back tensor
    tensor += tensor.T - diag

    return tensor

def stdcst(d, A, B):
    """
    Returns A double contracted with B
    A and B are second order tensor in Voigt representation
    With inditial notation : result = A_ij * B_ij
    """

    # Convert to tensor
    ddl = int(d*(d+1)/2)
    At = array2st(d, ddl, A)
    Bt = array2st(d, ddl, B)

    # Double contraction 
    result = np.einsum('ij, ij', At, Bt)

    return result 

def stkronst(A, B):
    """
    Returns kron product of tensors A and B
    A and B are second order tensor in Voigt representation
    """

    # Set kron product
    At = np.atleast_2d(A)
    Bt = np.atleast_2d(B)
    result = np.kron(At.T, Bt)

    return result

def compute_deviatoric(d, tensor):
    "Returns deviatoric of a second-order tensor "

    # Initialize
    ddl = int(d*(d+1)/2)
    one = np.zeros(ddl)

    # Compute trace of tensor and one 
    trace = 0.0
    for i in range(d):
        trace += tensor[i]
        one[i] = 1.0
        
    # Compute deviatoric
    dev = tensor - 1.0/3.0*trace*one

    return dev

def clean_dirichlet_3d(A, dod):

    """
    Set to 0 (Dirichlet condition) the values of an array using the indices in each dimension
    A is actually a vector arranged following each dimension [Au, Av, Aw]
    """

    for i in range(3):
        A[i, dod[i]] = 0.0

    return

def fourth_order_identity(d=3):
    " Creates a fourth-order identity (Voigt representation)"
    
    # Initialize 
    ddl = int(d*(d+1)/2)
    I = np.eye(ddl)

    return I

def one_kron_one(d=3): 
    "Creates a one kron one tensor (Voigt representation)"

    # Initialize
    ddl = int(d*(d+1)/2)
    one = np.zeros(ddl)
    one[:d] = 1.0

    # Define 
    onekronone = stkronst(one, one)

    return onekronone

def create_der2sym(d=3):

    """
    Creates M matrix. M is the passage matrix from derivative to actual symetric values. 
    If we multiply a vector of values u_(i, j) with M matrix, one obtains the vector: 
    us_ij = 0.5*(u_(i,j) + u_(j,i))  
    """

    # Create MM
    ddl = int(d*(d+1)/2)
    MM = np.zeros((ddl, d*d))

    if d==3: 
        MM[0, 0] = 1; MM[1, 4] = 1; MM[2, 8] = 1
        MM[3, 1] = 0.5; MM[3, 3] = 0.5
        MM[4, 5] = 0.5; MM[4, 7] = 0.5
        MM[5, 2] = 0.5; MM[5, 6] = 0.5
    elif d == 2: 
        MM[0, 0] = 1; MM[1, 3] = 1
        MM[2, 1] = 0.5; MM[2, 2] = 0.5

    return MM

def consistency_perfect_plasticity(sigma_Y, sigma, d=3):
    "Computes the value of f (consistency condition) and its derivatives in perfect plasticity criteria"

    # Compute deviatoric
    dev = compute_deviatoric(d, sigma)

    # Compute the norm of dev
    norm = np.sqrt(stdcst(d, dev, dev))

    # Compute f
    f = norm - np.sqrt(2.0/3.0)*sigma_Y

    # Compute gradient of f
    grad = 1/norm * dev

    # Compute gradient of the gradient
    devdev = stkronst(dev, dev)
    Iden = fourth_order_identity(d)

    gradgrad = 1.0/norm*Iden - 1.0/(norm**3)*devdev

    return f, grad, gradgrad

def cpp_perfect_plasticity(CC, SS, sigma_Y, e_n1, ep_n0, nbIter=20, d=3):
    "Return closest point proyection (cpp) in perfect plasticity criteria"

    # Compute elastic predictor
    diff_e_ep = e_n1 - ep_n0
    sigma_n1 = CC @ diff_e_ep

    # Review condition
    f, grad, grad2 = consistency_perfect_plasticity(sigma_Y, sigma_n1, d=d)

    if f <= 0.0: # Elastic point
        ep_n1 = ep_n0
        dSdE = CC
    else:
        # Initialize
        dgamma = 0.0
        ep_n1 = ep_n0

        # Return mapping iterative algorithm
        for i in range(nbIter):
            # Compute stress
            diff_e_ep = e_n1 - ep_n0
            sigma_n1 = CC @ diff_e_ep

            # Review condition
            f, grad, grad2 = consistency_perfect_plasticity(sigma_Y, sigma_n1)

            # Compute residual
            r_n1 = ep_n1 - ep_n0 - dgamma*grad

            # Check convergence
            norm = np.sqrt(stdcst(d, r_n1, r_n1))   
            if (norm < 1e-10) & (f < 1e-10): break

            # Compute tangent moduli
            xi_inv = SS + dgamma*grad2

            # Compute delta gamma
            xi_r = np.linalg.solve(xi_inv, r_n1)
            xi_grad = np.linalg.solve(xi_inv, grad)
            d2gamma = (f + stdcst(d, grad, xi_r))/stdcst(d, grad, xi_grad)

            # Compute increment
            diff_grad = d2gamma*grad - r_n1
            xi_diff = np.linalg.solve(xi_inv, diff_grad)
            d_ep = SS @ xi_diff

            # Update strains 
            dgamma += d2gamma
            ep_n1 += d_ep

        # Compute dSdE
        dSdE = np.linalg.inv(SS + dgamma*grad2)

        # Compute N
        xi_grad = dSdE @ grad
        N = xi_grad/np.sqrt(stdcst(d, grad, xi_grad))
        NNT = stkronst(N, N)

        # Update 
        dSdE -= NNT

    return ep_n1, sigma_n1, dSdE

def compute_plasticity_coef(sigma, dSdE, invJ, detJ, d=3):
    "Computes the coefficients to use in internal force vector and stiffness matrix"

    # Computes passage matrix
    MM = create_der2sym(d)

    # Initialize
    invJext = np.zeros((d*d, d*d))

    for i, det in enumerate(detJ):

        # Computes inverse of jacobian extended
        for j in range(d):
            invJext[j*d:(j+1)*d, j*d:(j+1)*d] = invJ[:, :, i]

        # Computes the coefficients to use in Fint vector
        coef_Fint = invJext.T @ (MM.T @ sigma) * det

        # Computes the coefficients to use in Stiffness matrix
        MDM = MM.T @ dSdE[:, :, i] @ MM
        coef_Stiff = invJext.T @ MDM @ invJext * det

    return coef_Fint, coef_Stiff

