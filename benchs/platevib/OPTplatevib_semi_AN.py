# Copyright 2020 Arnaud Duval

# This file is part of Yeti.
#
# Yeti is free software: you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version.
#
# Yeti is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with Yeti. If not, see <https://www.gnu.org/licenses/>
#
# #!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Shape of a 3D solid plate is optimized versus 2 target eigenfrequencies
Semi-analytical gradients are used
validation is performed to verify is target frequencies are obtained
"""

# Python module
import sys
import numpy as np

# IGA module
import yeti_iga.reconstructionSOL as rsol
from yeti_iga.preprocessing.igaparametrization import IGAparametrization,    \
                                             IGAmanip as manip,     \
                                             OPTmodelling
import yeti_iga.postprocessing.postproc as pp

if __name__ == "__main__":
    # Base path for .inp and .NB files
    FILENAME = 'plateVolume'

    # Creation of the IGA object
    modeleIGA = IGAparametrization(filename=FILENAME)

    # Refine IGA object along 1st and 2nd directions
    nb_deg = np.array([1, 1, 0])
    nb_ref = np.array([2, 2, 0])
    modeleIGA.refine(nb_ref, nb_deg)

    # Parametrization
    # Get index of CP at top and bottom
    botcps = manip.get_directionCP(modeleIGA, 5, 0, 0) - 1
    topcps = manip.get_directionCP(modeleIGA, 6 ,0, 0) - 1

    nb_var = botcps.size

    # Define thickness map depending on design variables
    H_MAX = 0.5
    H_MIN = 0.05

    def thickness(coords0, igapara, var):
        """
        Shape parametrization.
        Move control points as a function of esign variables corresponding to
        thickness
        """
        var = (H_MAX - H_MIN)*var + H_MIN

        igapara.coords[:, :] = coords0[:, :]
        igapara.coords[2, botcps] = -0.5 * ((H_MAX - H_MIN)*var + H_MIN)
        igapara.coords[2, topcps] = 0.5 * ((H_MAX-H_MIN)*var + H_MIN)

    # Refinement from optimization model to analysis model
    nb_deg = np.array([0, 0, 1])
    nb_ref = np.array([2, 2, 0])

    # Define optim problem
    # modeleIGA -> initial parametrization
    # nb_var    -> number of design variables
    # thickness -> function updating coordinates of the coarse (optim) model
    # nb_degreeElevationByDirection -> degree elevation for each diretion
    # nb_refinementByDirection -> refinement in each direction
    optPB = OPTmodelling(modeleIGA, nb_var, thickness,
                         nb_degreeElevationByDirection=nb_deg,
                         nb_refinementByDirection=nb_ref)

    # OPTIMIZATION

    from scipy.optimize import minimize

    # Initial values of design variables
    x0 = 0.6*np.ones(nb_var)
    # Compute initial 1st eigenvalue (m0) and eigenvector (v0)
    m0, v0 = optPB.compute_vibrationMode(x0)

    # Define target eigenvalue(s)
    NB_FRQCIBLE = 2
    xcible = 0.2*np.ones(nb_var)
    wcible, vcible = optPB.compute_vibrationMode(xcible, nb_frq=NB_FRQCIBLE)

    # Define a function ans its gradient w.r.t design variables
    # Gap with target eigenvalues
    def frecIGA(x_k):
        """
        Cost function
        Return relative difference between eigenfrequencies and reference
        values
        """
        valsi, _ = optPB.compute_vibrationMode(x_k, nb_frq=NB_FRQCIBLE)
        m_i = np.sum((valsi - wcible)**2.)
        return m_i / np.sum(wcible**2.)


    def gradFrecIGA(xk):
        """
        Gradient of cost function
        """
        gradF = np.zeros_like(xk)
        valsi, _ = optPB.compute_vibrationMode(xk, nb_frq=NB_FRQCIBLE)
        gradFreq = optPB.compute_gradVibration_semiAN(xk, eps=1.e-6,
                                                      nb_frq=NB_FRQCIBLE)
        for i in np.arange(NB_FRQCIBLE):
            gradF[:] += 2.*gradFreq[:, i]*(valsi[i] - wcible[i])
        return gradF / np.sum(wcible**2.)

    # Define a function and its gradient w.r.t design variables
    # Volume relative to initial volume
    vol0 = optPB.compute_volume(x0)

    def volIGA(x_k):
        """
        Constraint function : volume relative to initial volume
        """
        global vol0
        return optPB.compute_volume(x_k)/vol0


    def gradVolIGA(x_k):
        """
        Gradient of constraint function
        """
        global vol0
        return optPB.compute_gradVolume_AN(x_k)/vol0

    tab = []
    iopt = 0

    # Define callback function to be run at each optimization iteration (iopt)
    def saveXk(x_k):
        """
        Callback function.
        Save results at each optimization iteration
        """
        global iopt
        print(f'\nIteration{iopt:02}')
        valsi, vecti = optPB.compute_vibrationMode(x_k)
        tab.append([x_k[0], valsi[0]])

        # - Plot Thickness
        ticknessfield = np.zeros_like(optPB._coarseParametrization.coords)
        ticknessfield[2, :] = np.tile(
            optPB._coarseParametrization.coords[2, topcps]
            - optPB._coarseParametrization.coords[2, botcps],
            (2, 1)).flatten('C')
        pp.generatevtu(*optPB._coarseParametrization.get_inputs4postprocVTU(
            f'vibOpt{iopt:02}', ticknessfield,
            nb_ref=np.array([4, 4, 1]), Flag=np.array([True, False, False])))
        np.savetxt(f'results/cps{iopt:02}.txt',
                   optPB._coarseParametrization.coords.T,
                   delimiter=',')

        # - Plot Analysis

        sol, _ = rsol.reconstruction(
            **optPB._fineParametrization.get_inputs4solution(vecti))
        pp.generatevtu(*optPB._fineParametrization.get_inputs4postprocVTU(
            f'vibAn%0.2d{iopt:02}',  sol.transpose(),
            nb_ref=np.array([2, 2, 2]), Flag=np.array([True, False, False])))

        iopt += 1

    # Bounds for design variables : [0 ; 1]
    bds = ((0., 1.),)*nb_var

    # Run optimization

    # 1 - get given 1st eigenvalue
    res = minimize(frecIGA, x0, method='SLSQP', jac=gradFrecIGA,
                   bounds=bds, callback=saveXk, tol=1.e-9)

    if not res['success'] or res['fun'] > 0.02:
        print("res['success'] : ", res['success'])
        print("res['fun'] : ", res['fun'])
        sys.exit(-1)
    else:
        sys.exit(0)
