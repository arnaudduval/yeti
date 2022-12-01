# Copyright 2020 Thibaut Hirschler
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

# !/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This cas is described in the following publication :
Hirschler, T., Bouclier, R., Duval, A. et al.
A New Lighting on Analytical Discrete Sensitivities in the Context of
IsoGeometric Shape Optimization. Arch Computat Methods Eng (2020).
https://doi.org/10.1007/s11831-020-09458-6

The shape of a Kirchhoff-Love shell roof is optimized versus its maximal
displacement.
Volume is kept constant
Resuting shape is compared to reference numerical results

"""


# Python module
import numpy as np
import nlopt

# IGA module
from preprocessing.igaparametrization import IGAparametrization, \
                                             IGAmanip as manip
from preprocessing.igaparametrization import OPTmodelling
import postprocessing.postproc as pp
import reconstructionSOL as rsol


# Selection of .INP and .NB file
# ------------------------------
modeleIGA = IGAparametrization(filename='squareShellPlate')

nb_degDV = np.array([1, 1, 0])
nb_refDV = np.array([2, 2, 0])

modeleIGA.refine(nb_refDV, nb_degDV)


# --
# Shape Parametrization

vertex = manip.get_vertexCPindice(modeleIGA._Nkv, modeleIGA._Jpqr,
                                  modeleIGA._dim)[:4]
freecp = np.setxor1d(np.arange(0, modeleIGA._nb_cp), vertex)
nb_var = freecp.size


def altitude(coords0, igapara, var):
    """
    Change z coordinate of control points as a function of design variables
    """
    igapara._COORDS[:, :] = coords0[:, :]
    igapara._COORDS[2, freecp] = coords0[2, freecp] + 5.*var[:]


# --
# Build the optimization pb

nb_degAN = np.maximum(np.array([1, 1, 0])-nb_degDV, 0)
nb_refAN = np.maximum(np.array([5, 5, 0])-nb_refDV, 0)

optPB = OPTmodelling(modeleIGA, nb_var, altitude,
                     nb_degreeElevationByDirection=nb_degAN,
                     nb_refinementByDirection=nb_refAN)


# --
# Initialization and Definition of the objective and
# constraints (using nlopt)


P_NORM = 20
x0 = np.ones(nb_var)*0.1
n0 = optPB.compute_displacementAggreg(x0, pnorm=P_NORM)
V0 = 1.10*optPB.compute_volume(np.zeros(nb_var))

i = 0
ref_plot = np.array([2, 2, 1])
Output = np.array([True, False, False])


def maxdisp(x_d, grad_d):
    """
    Compute maximum displacement and its gradient
    """
    n_i = optPB.compute_displacementAggreg(x_d, pnorm=P_NORM)
    if grad_d.size > 0:
        global i
        i += 1
        print('\n--')
        print(f"Iter {i:03}")
        grad_d[:] = \
            optPB.compute_gradDisplacementAggreg_AN(x_d, pnorm=P_NORM)/n0

        # postprocessing
        pp.generatevtu(*optPB._coarseParametrization.get_inputs4postprocVTU(
            f"OPT5-coarse{i:02}",
            optPB._coarseParametrization._COORDS - optPB._initialCOORDS,
            nb_ref=2*ref_plot, Flag=Output))
        optPB._coarseParametrization.generate_vtk4controlMeshVisu(
            f"OPT5-coarse{i:02}", 0)

        sol, _ = rsol.reconstruction(
            **optPB._fineParametrization.get_inputs4solution(
                optPB._save_sol_fine))
        pp.generatevtu(*optPB._fineParametrization.get_inputs4postprocVTU(
            f"OPT5-fine{i:02}",  sol.transpose(),
            nb_ref=1*ref_plot, Flag=Output))
    return n_i/n0


def vol(x_v, grad_v):
    """
    Return volume and compute volume gradient as a function of design variables
    """
    if grad_v.size > 0:
        grad_v[:] = optPB.compute_gradVolume_AN(x_v)/V0
    return optPB.compute_volume(x_v)/V0 - 1.


minimize = nlopt.opt(nlopt.LD_SLSQP, nb_var)

minimize.set_min_objective(maxdisp)
minimize.add_inequality_constraint(vol, 1e-5)

minimize.set_ftol_rel(1.0e-06)
minimize.set_xtol_rel(1.0e-06)
minimize.set_maxeval(200)

minimize.set_lower_bounds(0.*np.ones(nb_var))
minimize.set_upper_bounds(1.*np.ones(nb_var))


# --
# Run optimization
x = minimize.optimize(x0)

# Verify results
# Numerical reference result
x_ref = np.array([0.1472491,  0.53377687, 0.53385731, 0.14725017, 0.14735826,
                  0.24791125, 0.50500495, 0.50501984, 0.24804221, 0.14727647,
                  0.53389746, 0.50526637, 0.66379548, 0.66400121, 0.5049124,
                  0.53418678, 0.53375444, 0.50543068, 0.66376945, 0.66413801,
                  0.504866,   0.53431571, 0.1474783,  0.2477457, 0.50501584,
                  0.50456767, 0.24810994, 0.14724689, 0.1472503,  0.53384915,
                  0.53436835, 0.14711689])
error = np.sqrt(sum((x-x_ref)**2.))

print(error)

assert error < 5.e-4
