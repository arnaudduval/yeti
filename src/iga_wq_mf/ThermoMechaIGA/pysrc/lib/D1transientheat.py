"""
.. This module contains functions to treat transient heat problems in 1D. 
.. It is adapted to IGA problems (using Gauss quadrature)
.. Joaquin Cornejo
"""

import numpy as np

def compute_volsource_1D(JJ, DB, W, Fprop):
    " Computes 'volumetric' source vector in 1D "
    Fcoefs = W * Fprop * JJ
    F = DB[0] @ Fcoefs 
    return F

def interpolate_temperature_1D(DB, T_ctrlpts):
    " Interpolate temperature in 1D "
    T_interp = DB[0].T @ T_ctrlpts
    return T_interp

def compute_thermal_Fint_1D(JJ, DB, W, Kprop, Cprop, T, dT):
    "Returns the internal heat force in transient heat"

    # Compute conductivity matrix
    Kcoefs = W * Kprop * 1.0/JJ
    K = DB[1] @ np.diag(Kcoefs) @ DB[1].T 

    # Compute capacitiy matrix 
    Ccoefs = W * Cprop * JJ
    C = DB[0] @ np.diag(Ccoefs) @ DB[0].T

    # Compute internal heat force 
    Fint = C @ dT + K @ T

    return Fint

def compute_tangent_thermal_matrix_1D(JJ, DB, W, Kprop, Cprop, newmark=0.5, dt=0.1):
    """ Computes tangent matrix in transient heat
        S = C + newmark dt K
        K = int_Omega dB/dx Kprop dB/dx dx = int_[0, 1] J^-1 dB/dxi Kprop J^-1 dB/dxi detJ dxi.
        But in 1D: detJ times J^-1 get cancelled.
        C = int_Omega B Cprop B dx = int [0, 1] B Cprop det J B dxi
    """

    # Compute conductivity matrix
    Kcoefs = W * Kprop * 1.0/JJ
    K = DB[1] @ np.diag(Kcoefs) @ DB[1].T # In fact we have to add radiation matrix and convection matrix

    # Compute capacitiy matrix 
    Ccoefs = W * Cprop * JJ
    C = DB[0] @ np.diag(Ccoefs) @ DB[0].T

    # Compute tangent matrix 
    M = C + newmark*dt*K

    return M

def solve_transient_heat_1D(properties, DB=None, W=None, Fext=None, time_list=None, 
                            dof=None, dod=None, Tinout=None, threshold=1e-12, nbIterNL=20):
    " Solves transient heat problem in 1D. "

    # Initialize
    JJ, conductivity, capacity, newmark = properties
    ddGG = np.zeros(len(dod)) # d Temperature/ d time
    VVn0 = np.zeros(len(dof))

    # Compute initial velocity from boundry conditions (for i = 0)
    if np.shape(Tinout)[1] == 2:
        delta_t = time_list[1] - time_list[0]
        ddGG = (Tinout[dod, 1] - Tinout[dod, 0])/delta_t
    elif np.shape(Tinout)[1] >= 2:
        delta_t1 = time_list[1] - time_list[0]
        delta_t2 = time_list[2] - time_list[0]
        factor = delta_t2/delta_t1
        ddGG = (Tinout[dod, 2] - (factor**2)*Tinout[dod, 1] - (1-factor**2)*Tinout[dod, 0])/(delta_t1*(factor-factor**2))
    else: raise Warning('We need more than 2 steps')

    for i in range(1, np.shape(Tinout)[1]):
        # Get delta time
        delta_t = time_list[i] - time_list[i-1]

        # Get values of last step
        TTn0 = Tinout[:, i-1]; TTn1 = np.copy(TTn0)

        # Predict values of new step
        TTn1[dof] = TTn0[dof] + delta_t*(1-newmark)*VVn0; TTn1[dod] = Tinout[dod, i]
        TTn1i0 = np.copy(TTn1); ddTT = np.zeros(len(TTn1))
        ddTT[dod] = 1.0/newmark*(1.0/delta_t*(Tinout[dod, i] - Tinout[dod, i-1]) - (1-newmark)*ddGG)
        VVn1 = np.zeros(len(dof))

        # Get "force" of new step
        Fstep = Fext[:, i]

        # Corrector (Newton-Raphson)
        for j in range(nbIterNL):

            # Interpolate temperature
            T_interp = interpolate_temperature_1D(DB, TTn1)

            # Get capacity and conductivity properties
            Kprop = conductivity(T_interp)
            Cprop = capacity(T_interp)

            # Compute residue
            Fint = compute_thermal_Fint_1D(JJ, DB, W, Kprop, Cprop, TTn1, ddTT)
            dF = Fstep[dof] - Fint[dof]
            prod1 = np.dot(dF, dF)
            relerror = np.sqrt(prod1)
            if relerror <= threshold: break

            # Compute tangent matrix
            MM = compute_tangent_thermal_matrix_1D(JJ, DB, W, 
                    Kprop, Cprop, newmark=newmark, dt=delta_t)[np.ix_(dof, dof)]

            # Compute delta dT 
            ddVV = np.linalg.solve(MM, dF)

            # Update values
            VVn1 += ddVV
            TTn1[dof] = TTn1i0[dof] + newmark*delta_t*VVn1
            ddTT[dof] += VVn1 

        print(j+1, relerror)

        # Update values in output
        Tinout[:, i] = TTn1
        VVn0 = np.copy(VVn1)
        ddGG = np.copy(ddTT[dod])
        
    return 
