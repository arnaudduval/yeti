import numpy as np

def return_mapping_point(E, K, H, sigma_Y, eps_n1, eps_pn0=None, alpha_n0=None, q_n0=None):
    """Return mapping algorithm for one-dimensional, rate-independent plasticity. 
    Combined Isotropic/Kinematic hardening"""
    
    # Number of points 
    nnz = len(eps_n1)

    if eps_pn0 is None: eps_pn0 = np.zeros(nnz)
    if alpha_n0 is None: alpha_n0 = np.zeros(nnz)
    if q_n0 is None: q_n0 = np.zeros(nnz)

    # Initialize output
    coef, sigma_n1, eps_pn1, alpha_n1, q_n1 = np.zeros(nnz), np.zeros(nnz), eps_pn0, alpha_n0, q_n0

    # Compute elastic trial stress and test for plastic loading
    sigma_trial = E*(eps_n1 - eps_pn0)
    xi_trial = sigma_trial - q_n0
    f_trial = np.abs(xi_trial) - (sigma_Y + K*alpha_n0)

    for i in range(nnz):
        if f_trial[i] <= 0:
            # Elastic step
            sigma_n1[i] = sigma_trial[i]
        else: 
            # Plastic step
            sign = np.sign(xi_trial[i])
            delta_gamma = f_trial[i]/(E + K + H)
            sigma_n1[i] = sigma_trial[i] - delta_gamma*E*sign
            eps_pn1[i] = eps_pn0[i] + delta_gamma*sign
            alpha_n1[i] = alpha_n0[i] + delta_gamma
            q_n1[i] = q_n0[i] + delta_gamma*H*sign

    for i in range(nnz):
        if f_trial[i] <= 0: coef[i] = E
        else: coef[i] = E*(K + H)/(E + K + H)
    
    return [coef, sigma_n1, eps_pn1, alpha_n1, q_n1]

def compute_Fint(DB, W, sigma_n1):
    "Returns vector F int "

    # Compute F int
    Fint = DB[1] @ np.diag(W) @ sigma_n1.T

    return Fint

def compute_Fext(DB, W, b, sigma_ext): 
    " Returns vector F ext"

    Fext = DB[0] @ np.diag(W) @ b.T
    Fext += sigma_ext

    return Fext

def compute_K(L, DB, W, coefs):
    " Returns matrix K "

    # Compute WC
    WC = W * coefs
    K = DB[1] @ np.diag(WC) @ DB[1].T / L

    return K

def newton_raphson(L, E, K, H, sigma_Y, DB, W, d_n0, b, sigma_ext, dof, 
                    eps_pn0, alpha_n0, q_n0, nbIter=100, threshold=1e-6):
    " Computes displacement at tn+1"

    # Define the first candidate to be dn1
    d_n1 = d_n0

    for _ in range(nbIter):
        # Compute strain field at quadrature points
        eps_n1 = DB[1].T @ d_n1/L

        # Compute variables using return maping for each quadrature point
        coefs, sigma_n1, eps_pn1, alpha_n1, q_n1 = return_mapping_point(E, K, H, sigma_Y, eps_n1, eps_pn0, alpha_n0, q_n0)

        # Update values
        eps_pn0, alpha_n0, q_n0 = eps_pn1, alpha_n1, q_n1

        # Compute Fint and Fext
        Fint = compute_Fint(DB, W, sigma_n1)[dof]
        Fext = compute_Fext(DB, W, b, sigma_ext)[dof]
        error = np.linalg.norm(Fint - Fext, np.inf)/np.linalg.norm(Fext, np.inf)

        if error < threshold: 
            d_n1 = d_n0
            break
        else:
            # Compute coefs
            KK = compute_K(L, DB, W, coefs)[np.ix_(dof, dof)]
            diffF = Fint - Fext
            delta_d_dof = -np.linalg.solve(KK, diffF)
            d_n1[dof] = d_n0[dof] + delta_d_dof
        
    return d_n1, sigma_n1, eps_pn1, alpha_n1, q_n1

def interpolate_CP(DB, W, u_ref):
    " Returns the values computed at control points "

    M = DB[0] @ np.diag(W) @ DB[0].T
    WU = W*u_ref
    R = DB[0] @ WU

    sol = np.linalg.solve(M, R)

    return sol