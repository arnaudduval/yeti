"""
.. This module contains functions to treat elastoplasticity problems in 1D. 
.. It is adapted to IGA problems (using Gauss quadrature)
.. Joaquin Cornejo
"""

import numpy as np

def cpp_combined_hardening_1D(E, H, beta, sigma_Y0, deps, sigma_n0, alpha_n0, ep_n0):
    """ Return mapping algorithm for one-dimensional rate-independent plasticity. 
        It uses combined isotropic/kinematic hardening theory.  
    """

    # Elastic predictor
    sigma_trial = sigma_n0 + E*deps 
    eta_trial = sigma_trial - alpha_n0

    # Check yield status
    fthreshold = sigma_Y0*1e-6
    f_trial = abs(eta_trial) - (sigma_Y0 + (1 - beta)*H*ep_n0)

    if f_trial <= fthreshold: # Elastic
        sigma_n1 = sigma_trial
        alpha_n1 = alpha_n0
        ep_n1 = ep_n0
        Dalg = E
    else: # Plastic
        dep = f_trial/(E + H)
        sigma_n1 = sigma_trial - np.sign(eta_trial)*E*dep
        alpha_n1 = alpha_n0 + np.sign(eta_trial)*beta*H*dep
        ep_n1 = ep_n0 + dep
        Dalg = E*H/(E + H)

    return [Dalg, sigma_n1, alpha_n1, ep_n1]

def compute_static_Fint_1D(DB, W, sigma):
    """ Computes internal force Fint. 
        Fint = int_Omega dB/dx sigma dx = int_[0, 1] J^-1 dB/dxi sigma detJ dxi.
        But in 1D: detJ times J^-1 get cancelled.
    """
    Fint = DB[1] @ np.diag(W) @ sigma.T
    return Fint

def compute_volForce_1D(DB, W, b): 
    " Computes volumetric force in 1D. Usualy b = rho*g*L (Weight) "
    F = DB[0] @ np.diag(W) @ b.T
    return F

def compute_tangent_static_matrix_1D(JJ, DB, W, Dalg):
    """ Computes stiffness matrix in elastoplasticity
        S = int_Omega dB/dx Dalg dB/dx dx = int_[0, 1] J^-1 dB/dxi Dalg J^-1 dB/dxi detJ dxi.
        But in 1D: detJ times J^-1 get cancelled.
    """
    Scoefs = W*Dalg*1.0/JJ
    S = DB[1] @ np.diag(Scoefs) @ DB[1].T 
    return S

def interpolate_strain_1D(JJ, DB, disp):
    " Computes strain field from a given displacement field "
    eps = DB[1].T @ disp / JJ
    return eps

def solve_plasticity_1D(properties, DB=None, W=None, Fext=None, dof=None, tol=1e-6, nbIterNL=200):
    " Solves elasto-plasticity problem in 1D. It considers Dirichlet boundaries equal to 0 "

    # Initialize
    JJ, E, H, beta, sigma_Y0, nb_qp = properties
    sigma_n0, alpha_n0, ep_n0 = np.zeros(nb_qp), np.zeros(nb_qp), np.zeros(nb_qp)
    sigma_n1, alpha_n1, ep_n1 = np.zeros(nb_qp), np.zeros(nb_qp), np.zeros(nb_qp)
    Dalg = np.zeros(nb_qp)
    disp = np.zeros(np.shape(Fext))

    # Set some variables to return
    ep = np.zeros((nb_qp, np.shape(Fext)[1]))
    sigma = np.zeros((nb_qp, np.shape(Fext)[1]))

    for i in range(1, np.shape(Fext)[1]):

        # Initialize
        ddisp = np.zeros(np.shape(disp[:, i-1]))
        Fstep = Fext[:, i]
        prod2 = np.dot(Fstep, Fstep)

        # Solver Newton-Raphson
        for j in range(nbIterNL):

            # Compute strain as function of displacement
            deps = interpolate_strain_1D(JJ, DB, ddisp)

            # Find closest point projection 
            for k in range(nb_qp):
                result = cpp_combined_hardening_1D(E, H, beta, sigma_Y0, 
                            deps[k], sigma_n0[k], alpha_n0[k], ep_n0[k])
                Dalg[k], sigma_n1[k], alpha_n1[k], ep_n1[k] = result

            # Compute Fint
            Fint = compute_static_Fint_1D(DB, W, sigma_n1)
            dF = Fstep[dof] - Fint[dof]
            prod1 = np.dot(dF, dF)
            relerror = np.sqrt(prod1/prod2)
            if relerror <= sigma_Y0*tol: break

            # Compute stiffness
            S = compute_tangent_static_matrix_1D(JJ, DB, W, Dalg)[np.ix_(dof, dof)]
            ddisp[dof] = np.linalg.solve(S, dF)
                
        print(j+1, relerror)

        # Update values in output
        disp[:, i] = disp[:, i-1] + ddisp
        ep[:, i] = ep_n1
        sigma[:, i] = sigma_n1
        ep_n0 = ep_n1
        sigma_n0 = sigma_n1

    return disp, ep, sigma

def interpolate_controlPoints_1D(DB, W, u_ref):
    " Interpolate control point field (from parametric to physical space) "
    masse = DB[0] @ np.diag(W) @ DB[0].T
    force = DB[0] @ np.diag(W) @ u_ref
    u_ctrlpts = np.linalg.solve(masse, force)
    return u_ctrlpts