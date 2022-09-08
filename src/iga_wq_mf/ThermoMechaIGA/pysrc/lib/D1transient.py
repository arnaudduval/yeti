"""
.. Integrating elastoplasticity 1D in python
"""

from copy import deepcopy
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

def compute_transient_Fint_1D(JJ, DB, W, Kcoefs, Ccoefs, T, dT):
    "Returns the residue in transient heat"

    # Compute conductivity matrix
    k_coefs = W * Kcoefs * 1.0/JJ
    K = DB[1] @ np.diag(k_coefs) @ DB[1].T 

    # Compute capacitiy matrix 
    c_coefs = W * Ccoefs * JJ
    C = DB[0] @ np.diag(c_coefs) @ DB[0].T

    # Compute tangent matrix 
    Fint = C @ dT + K @ T

    return Fint

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
    M = C + alpha*dt*K

    return M

def solve_transient_1D(properties, DB=None, W=None, Fext=None, time_list=None, dof=None, dod=None, Tin=None, tol=1e-4, nbIter=100):
    " Solves transient heat problem in 1D. It only solves Dirichlet (g=0) boundary."

    # Initialize
    JJ, conductivity, capacity, alpha = properties
    ddGG = np.zeros(len(dod)) # d Temperature/ d time
    VVn0 = np.zeros(len(dof))

    # Compute initial velocity from boundry conditions (for i = 0)
    if np.shape(Tin)[1] == 2:
        delta_t = time_list[1] - time_list[0]
        ddGG = (Tin[dod, 1] - Tin[dod, 0])/delta_t
    elif np.shape(Tin)[1] >= 2:
        delta_t1 = time_list[1] - time_list[0]
        delta_t2 = time_list[2] - time_list[0]
        factork = delta_t2/delta_t1
        ddGG = (Tin[dod, 2] - (factork**2)*Tin[dod, 1] - (1-factork**2)*Tin[dod, 0])/(delta_t1*(factork-factork**2))
    else:
        raise Warning('We need more than 2 steps')

    for i in range(1, np.shape(Tin)[1]):
        # Get delta time
        delta_t = time_list[i] - time_list[i-1]

        # Get values of last step
        TTn0 = Tin[:, i-1]; TTn1 = TTn0

        # Predict values of new step
        TTn1[dof] = TTn0[dof] + delta_t*(1-alpha)*VVn0
        TTn1[dod] = Tin[dod, i]
        TTn1i0 = deepcopy(TTn1)
        ddTT = np.zeros(len(TTn1))
        ddTT[dod] = 1.0/alpha*(1.0/delta_t*(Tin[dod, i] - Tin[dod, i-1]) - (1-alpha)*ddGG)
        VVn1 = np.zeros(len(dof))

        # Get "force" of new step
        F = Fext[:, i]
        prod2 = np.dot(F, F)

        # Corrector (Newton-Raphson)
        for j in range(nbIter):

            # Interpolate temperature
            T_qp = compute_temperature_1D(DB, TTn1)

            # Get capacity and conductivity coefficients
            Kcoefs = conductivity(T_qp)
            Ccoefs = capacity(T_qp)

            # Compute residue
            Fint = compute_transient_Fint_1D(JJ, DB, W, Kcoefs, Ccoefs, TTn1, ddTT)
            dF = F[dof] - Fint[dof]
            prod1 = np.dot(dF, dF)
            relerror = np.sqrt(prod1/prod2)

            if relerror <= tol:
                break
            else:
                # Compute tangent matrix
                MM = compute_tangent_transient_matrix_1D(JJ, DB, W, 
                                        Kcoefs, Ccoefs, alpha=alpha, dt=delta_t)[np.ix_(dof, dof)]

                # Compute delta dT 
                ddVV = np.linalg.solve(MM, dF)

                # Update values
                VVn1 += ddVV
                TTn1[dof] = TTn1i0[dof] + alpha*delta_t*VVn1
                ddTT[dof] += VVn1 

        print(j+1, relerror)

        # Update values in output
        Tin[:, i] = TTn1
        VVn0 = VVn1
        ddGG = ddTT[dod]
        
    return 
