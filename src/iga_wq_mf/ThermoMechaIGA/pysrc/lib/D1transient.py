"""
.. Integrating elastoplasticity 1D in python
"""

import numpy as np

def compute_source_1D(JJ, DB, W, Fcoefs):
    "Computes source vector in 1D"
    f_coefs = W * Fcoefs * JJ
    F = DB[0] @ f_coefs 
    return F

def compute_temperature_1D(DB, T_ctrlpts):
    " Computes temperature "

    # Compute temperature
    T_qp = DB[0].T @ T_ctrlpts
    return T_qp

def compute_transient_residue_1D(JJ, DB, W, Kcoefs, Ccoefs, T, dT, F):
    "Returns the residue in transient heat"

    # Compute conductivity matrix
    k_coefs = W * Kcoefs * 1.0/JJ
    K = DB[1] @ np.diag(k_coefs) @ DB[1].T 

    # Compute capacitiy matrix 
    c_coefs = W * Ccoefs * JJ
    C = DB[0] @ np.diag(c_coefs) @ DB[0].T

    # Compute tangent matrix 
    dF = F - C @ dT - K @ T

    return dF

def compute_tangent_transient_matrix_1D(JJ, DB, W, Kcoefs, Ccoefs, alpha=0.5, dt=0.1):
    """Returns tangent matrix in transient heat
    S = C + alpha dt K
    K = int_Omega dB/dx Kcoefs dB/dx dx = int_[0, 1] J^-1 dB/dxi Kcoefs J^-1 dB/dxi detJ dxi.
    But in 1D: detJ times J^-1 get cancelled.
    C = int_Omega B Ccoefs C dx = int [0, 1] B Ccoefs det J B dxi
    """

    # Compute conductivity matrix
    k_coefs = W * Kcoefs * 1.0/JJ
    K = DB[1] @ np.diag(k_coefs) @ DB[1].T # In fact we have to add radiation matrix and convection matrix

    # Compute capacitiy matrix 
    c_coefs = W * Ccoefs * JJ
    C = DB[0] @ np.diag(c_coefs) @ DB[0].T

    # Compute tangent matrix 
    S = C + alpha*dt*K

    return S

def solve_transient_1D(properties, DB=None, W =None, Fext=None, time_list=None, dof=None, tol=1e-6, nbIter=200):
    " Solves transient heat problem in 1D. It only solves Dirichlet (g=0) boundary."

    # Initialize
    JJ, conductivity, capacity, alpha = properties
    TT = np.zeros(np.shape(Fext)) # Temperature
    ddTT = np.zeros(np.shape(Fext)) # d Temperature/ d time

    for i in range(1, np.shape(Fext)[1]):
        # Get delta time
        delta_t = time_list[i] - time_list[i-1]

        # Get "force" of new step
        F = Fext[:, i]
        prod2 = np.dot(F, F)

        # Get values of last step
        TTn0 = TT[:, i-1] 
        ddTTn0 = ddTT[:, i-1]
        d_dT = np.zeros(np.shape(TTn0))

        # Predict values of new step
        TTn1 = TTn0 + delta_t*(1-alpha)*ddTTn0
        ddTTn1 = np.zeros(np.shape(TTn1))

        # Corrector (Newton-Raphson)
        for i in range(nbIter):

            # Interpolate temperature
            T_nq = compute_temperature_1D(DB, TTn1)

            # Get capacity and conductivity coefficients
            Kcoefs = conductivity(T_nq)
            Ccoefs = capacity(T_nq)

            # Compute tangent matrix
            SS = compute_tangent_transient_matrix_1D(JJ, DB, W, 
                                    Kcoefs, Ccoefs, alpha=alpha, dt=delta_t)[np.ix_(dof, dof)]

            # Compute residue
            dF = compute_transient_residue_1D(JJ, DB, W, Kcoefs, Ccoefs, TTn1, ddTTn1, F)[dof]
            prod1 = np.dot(dF, dF)
            relerror = np.sqrt(prod1/prod2)

            if relerror <= tol:
                break
            else:
                # Compute delta dT 
                d_dT[dof] = np.linalg.solve(SS, dF)

                # Update values
                ddTTn1 += d_dT
                TTn1 += alpha*delta_t*ddTTn1

        # Update values in output
        TT[:, i] = TTn1
        ddTT[:, i] = ddTTn1

    return TT



