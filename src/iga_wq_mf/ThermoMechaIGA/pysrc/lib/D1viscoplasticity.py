"""
.. Integrating elastoplasticity 1D in python
"""

import numpy as np

def return_mapping_point_1D(E, H, beta, sigma_Y0, deps, sigma_n0, alpha_n0, ep_n0):
    """Return mapping algorithm for one-dimensional rate-independent plasticity. 
    Combined Isotropic/Kinematic hardening. For one quadrature point """
    
    # Elastic predictor
    sigma_trial = sigma_n0 + E*deps 
    eta_trial = sigma_trial - alpha_n0

    # Check yield status
    ftol = sigma_Y0*1e-6
    f_trial = abs(eta_trial) - (sigma_Y0 +(1-beta)*H*ep_n0)

    if f_trial <= ftol: # Elastic
        sigma_n1 = sigma_trial
        alpha_n1 = alpha_n0
        ep_n1 = ep_n0
        Dalg = E
    else:
        dep = f_trial/(E + H)
        sigma_n1 = sigma_trial - np.sign(eta_trial)*E*dep
        alpha_n1 = alpha_n0 + np.sign(eta_trial)*beta*H*dep
        ep_n1 = ep_n0 + dep
        Dalg = E*H/(E+H)

    return [Dalg, sigma_n1, alpha_n1, ep_n1]

def compute_static_Fint_1D(DB, W, sigma):
    """Returns vector F int. 
    Fint = int_Omega dB/dx sigma dx = int_[0, 1] J^-1 dB/dxi sigma detJ dxi.
    But in 1D: detJ times J^-1 get cancelled.
    """

    # Compute F int
    Fint = DB[1] @ np.diag(W) @ sigma.T

    return Fint

def compute_Fvol_1D(DB, W, b): 
    " Returns volumetric force"
    F = DB[0] @ np.diag(W) @ b.T
    return F

def compute_tangent_static_matrix_1D(JJ, DB, W, Scoefs):
    """Returns stiffness matrix 
    S = int_Omega dB/dx Dalg dB/dx dx = int_[0, 1] J^-1 dB/dxi Dalg J^-1 dB/dxi detJ dxi.
    But in 1D: detJ times J^-1 get cancelled.
    """

    # Compute WC
    coefs = W * Scoefs * 1.0/JJ
    S = DB[1] @ np.diag(coefs) @ DB[1].T 

    return S

def compute_strain_1D(JJ, DB, disp):
    " Computes strain given displacement field "

    # Compute strain
    eps = DB[1].T @ disp / JJ
    return eps

def solve_plasticity_1D(properties, DB=None, W =None, Fext=None, dof=None, tol=1e-6, nbIter=200):
    " Solves plasticity problem in 1D "

    # Initialize
    JJ, E, H, beta, sigma_Y0, nnz = properties
    sigma_n0, alpha_n0, ep_n0 = np.zeros(nnz), np.zeros(nnz), np.zeros(nnz)
    sigma_n1, alpha_n1, ep_n1 = np.zeros(nnz), np.zeros(nnz), np.zeros(nnz)
    Dalg = np.zeros(nnz)
    disp = np.zeros(np.shape(Fext))

    # Set some variables to return
    ep = np.zeros((nnz, np.shape(Fext)[1]))
    sigma = np.zeros((nnz, np.shape(Fext)[1]))

    for i in range(1, np.shape(Fext)[1]):

        # Initialize
        ddisp = np.zeros(np.shape(disp[:, i-1]))
        F = Fext[:, i]
        prod2 = np.dot(F, F)

        # Newton Raphson
        for j in range(nbIter):

            # Compute strain as function of displacement
            deps = compute_strain_1D(JJ, DB, ddisp)

            # Closest point projection in perfect plasticity
            for k in range(nnz):
                result = return_mapping_point_1D(E, H, beta, sigma_Y0, 
                            deps[k], sigma_n0[k], alpha_n0[k], ep_n0[k])
                Dalg[k], sigma_n1[k], alpha_n1[k], ep_n1[k] = result

            # Compute Fint
            Fint = compute_static_Fint_1D(DB, W, sigma_n1)
            dF = F[dof] - Fint[dof]
            prod1 = np.dot(dF, dF)
            relerror = np.sqrt(prod1/prod2)

            if relerror <= sigma_Y0*tol:
                break
            else:
                # Compute stiffness
                S = compute_tangent_static_matrix_1D(JJ, DB, W, Dalg)[np.ix_(dof, dof)]
                ddisp[dof] = np.linalg.solve(S, dF)
                
        print(i, j)
        # Set values
        disp[:, i] = disp[:, i-1] + ddisp
        ep[:, i] = ep_n1
        sigma[:, i] = sigma_n1
        ep_n0 = ep_n1
        sigma_n0 = sigma_n1
 
    return disp, ep, sigma

def interpolate_controlPoints_1D(DB, W, u_ref):
    " Returns the values computed at control points "

    M = DB[0] @ np.diag(W) @ DB[0].T
    WU = np.diag(W) @ u_ref
    R = DB[0] @ WU

    u_ctrlpts = np.linalg.solve(M, R)

    return u_ctrlpts