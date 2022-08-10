"""
.. Integrating elastoplasticity in python
"""

import numpy as np

def clean_dirichlet_3d(A, dod):

    """
    Set to 0 (Dirichlet condition) the values of an array using the indices in each dimension
    A is actually a vector arranged following each dimension [Au, Av, Aw]
    """

    for i in range(3):
        A[i, dod[i]] = 0.0

    return

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

def compute_stress_deviatoric(d, tensor):
    "Returns deviatoric of a second-order stress-like tensor "

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

def compute_stress_norm(d, tensor): 
    "Returns frobenius norm of a second-order stress-like tensor "

    # Initialize
    ddl = int(d*(d+1)/2)
    norm = 0.0

    for i in range(d):
        norm += tensor(i)**2

    for i in range(d, ddl):
        norm += 2.0*tensor(i)**2
    
    norm = np.sqrt(norm)

    return norm

def fourth_order_identity(d=3):
    " Creates a fourth-order identity (Voigt representation)"
    
    # Initialize 
    ddl = int(d*(d+1)/2)
    I = np.diag(np.array([1.0, 1.0, 1.0, 0.5, 0.5, 0.5]))

    return I

def one_kron_one(d=3): 
    "Creates a one kron one tensor (Voigt representation)"

    # Initialize
    ddl = int(d*(d+1)/2)
    onekronone = np.zeros((ddl, ddl))

    for i in range(d):
        for j in range(d):
            one_kron_one[i,j] = 1.0

    return onekronone

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

def create_incidence(d=3):

    """
    Creates M matrix. M is the passage matrix from derivative to actual symetric values. 
    If we multiply a vector of values u_(i, j) with M matrix, one obtains the vector: 
    us_ij = 0.5*(u_(i,j) + u_(j,i))  
    """

    # Create MM
    ddl = int(d*(d+1)/2)
    EE = np.zeros((ddl, d, d))

    if d==3: 
        EE[0, 0, 0] = 1.0; EE[4, 2, 0] = 1.0; EE[5, 1, 0] = 1.0
        EE[1, 1, 1] = 1.0; EE[3, 2, 1] = 1.0; EE[5, 0, 1] = 1.0
        EE[2, 2, 2] = 1.0; EE[3, 1, 0] = 1.0; EE[4, 0, 0] = 1.0
    elif d == 2: 
        EE[0, 0, 0] = 1.0; EE[2, 1, 0] = 1.0
        EE[1, 1, 1] = 1.0; EE[2, 0, 1] = 1.0

    return EE

def cpp_combined_hardening(inputs, deps, sigma_n, ep_n, alpha_n, d=3):
    "Return closest point proyection (cpp) in combined hardening plasticity criteria"

    def yield_function(sigma_Y, beta, H, sigma, alpha, ep_n, d=3):
        "Computes the value of f (consistency condition) in perfect plasticity criteria"

        # Compute deviatoric
        eta = compute_stress_deviatoric(d, sigma) - alpha

        # Compute the norm of dev
        norm = compute_stress_norm(d, eta)

        # Compute f
        f = norm - np.sqrt(2.0/3.0)*(sigma_Y + (1-beta)*H*ep_n)

        # Compute unit deviatoric tensor
        if norm > 0.0: N = 1.0/norm*eta
        else: N = np.zeros(np.shape(eta))

        return f, N, norm

    # Initialize
    CC = inputs[0]
    sigma_Y = inputs[1]
    mu = inputs[2]
    beta = inputs[3]
    H = inputs[4]
    Idev = inputs[5]

    # Compute trial sigma 
    sigma_trial = sigma_n + CC @ deps

    # Compute yield function and unit deviatoric tensor
    f, N, norm = yield_function(sigma_Y, beta, H, sigma_trial, alpha_n, ep_n, d=d)

    # Check status
    if f < 0:
        # Copy old variables
        sigma_n1 = sigma_trial
        ep_n1 = ep_n
        alpha_n1 = alpha_n
        Dalg = CC
    else:
        # Consistency parameter
        dgamma = f/(2.0*mu+2.0/3.0*H)

        # Update stress
        sigma_n1 = sigma_trial - 2*mu*dgamma*N

        # Update back stress
        alpha_n1 = alpha_n + 2.0/3.0*beta*H*dgamma*N

        # Update effective plastic strain
        ep_n1 = ep_n + np.sqrt(2.0/3.0)*dgamma

        # Compute consistent tangent matrix
        c1 = 4.0*mu**2.0/(2.0*mu+2.0/3.0*H)
        c2 = 4.0*mu**2.0*dgamma/norm
        NNT = stkronst(N, N)
        Dalg = CC - (c1 - c2)*NNT - c2*Idev

    return sigma_n1, ep_n1, alpha_n1, Dalg

def compute_plasticity_coef(sigma, Dalg, invJ, detJ, d=3):
    "Computes the coefficients to use in internal force vector and stiffness matrix"

    # Computes passage matrix
    EE = create_incidence(d)

    # Initialize
    coef_Fint = np.zeros((d*d, len(detJ)))
    coef_Stiff = np.zeros((d*d, d*d, len(detJ)))

    for k, det in enumerate(detJ):

        for i in range(d):
            for j in range(d):
                Dij = invJ[:,:,k].T @ EE[:,:,i].T @ Dalg[:,:,k] @ EE[:,:,j] @ invJ[:,:,k]
                coef_Stiff[i*d:(i+1)*d, j*d:(j+1)*d] = Dij*det

            Si = invJ[:,:,k].T @ EE[:,:,i].T @ sigma[:,k]
            coef_Fint[i*d:(i+1)*d] = Si*det

    return coef_Fint, coef_Stiff

