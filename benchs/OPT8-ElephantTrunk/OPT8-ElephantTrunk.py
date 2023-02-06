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
# FOR A PARTICULAR
# PURPOSE. See the GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with Yeti. If not, see <https://www.gnu.org/licenses/>

# !/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This case is described in the following publication :
Hirschler, T., Bouclier, R., Duval, A. et al.
A New Lighting on Analytical Discrete Sensitivities in the Context of
IsoGeometric Shape Optimization.
Arch Computat Methods Eng (2020). https://doi.org/10.1007/s11831-020-09458-6

The shape of a solid 3D elephant trunk is optimized versus its eigenfrequencies
Volume is kept constant
Resuting shape is compared to reference numerical results

"""


# Python module
import numpy as np
import sys
import time
import nlopt

# IGA module
from preprocessing.igaparametrization import IGAparametrization, IGAmanip as manip
import postprocessing.postproc as pp
import reconstructionSOL as rsol

from preprocessing.igaparametrization import OPTmodelling


# Selection of .INP and .NB file
# ------------------------------
modeleIGA = IGAparametrization(filename='elephantTrunk')

nb_degDV = np.array([0, 0, 1])
nb_refDV = np.array([0, 0, 2])

modeleIGA.refine(nb_refDV, nb_degDV)

# --
# Shape Parametrization

nb_var = np.int64((modeleIGA._Nkv[2, 0] - modeleIGA._Jpqr[2, 0] - 1)*2)


def dilatation(coords0, igapara, var):
    igapara._COORDS[:, : ] = coords0[:, :]
    for i in np.arange(0, var.size/2, dtype=np.intp):
        icps = np.intp(manip.get_directionCP(igapara, 5, 0, i) - 1)
        igapara._COORDS[0, icps] *= var[i]
        igapara._COORDS[1, icps] *= var[int(i+var.size/2)]
    return None


# --
# Build the optinization pb

nb_degAN = np.maximum(np.array([0, 0, 1]) - nb_degDV, 0)
nb_refAN = np.maximum(np.array([2, 2, 4]) - nb_refDV, 0)

optPB = OPTmodelling(modeleIGA, nb_var, dilatation,
                     nb_degreeElevationByDirection=nb_degAN,
                     nb_refinementByDirection=nb_refAN)

# --
# Initialization and Definition of the objective and constraints (using nlopt)

Pnorm = 20
nb_frqAggreg = 2
nb_frqPlot = 5
x0 = np.block([np.ones(int(nb_var/2))*1.25, np.ones(int(nb_var/2))])
vals0, vect0 = optPB.compute_vibrationMode(x0, nb_frqPlot)
m0 = np.sum((1./vals0[:nb_frqAggreg])**Pnorm)**(1./Pnorm)

i = 0
freqHistory = np.array([vals0])
ref_plot = np.array([2, 2, 2])
Output = np.array([True, False, False])


def MAC(vectA, vectB):
    mat = np.dot(vectA[:, :].T, vectB[:, :])**2.
    mat /= np.tile(np.diag(vectA[:, :].T.dot(vectA[:, :])),
                   (vectA.shape[1], 1)).T
    mat /= np.tile(np.diag(vectB[:, :].T.dot(vectB[:, :])), (vectB.shape[1], 1))
    return mat


def vibration(xV, gradV):
    global i, m0, freqHistory, vect0
    valsi, vecti = optPB.compute_vibrationMode(xV, nb_frqPlot)
    test_switch = np.argmax(MAC(vecti, vect0), axis=1)
    valsi[:] = valsi[test_switch[:]]
    vecti[:, :] = vecti[:, test_switch[:]]
    mi = np.sum((1./valsi[:nb_frqAggreg])**Pnorm)**(1./Pnorm)
    if gradV.size > 0:
        i += 1
        print('\n--')
        print('Iter %3i' % i)

        gradFq = optPB.compute_gradVibration_AN(xV, nb_frq=nb_frqPlot)
        gradV[:] = 0.
        for n in np.arange(nb_frqAggreg):
            gradV[:] -= gradFq[:, test_switch[n]]*(1./valsi[n])**(Pnorm+1)
        gradV[:] *= np.sum((1./valsi[:nb_frqAggreg])**Pnorm)**(1./Pnorm-1)
        gradV /= m0

        pp.generatevtu(*optPB._coarseParametrization.get_inputs4postprocVTU(
            'OPT8-coarse%0.2d' % i, optPB._coarseParametrization._COORDS,
            nb_ref=2*ref_plot, Flag=Output))
        optPB._coarseParametrization.generate_vtk4controlMeshVisu(
            'OPT8-coarse%0.2d' % i, 0)

        for n in np.arange(nb_frqPlot):
            ampl = np.sign(vect0[:, n].dot(vecti[:, n]))
            SOL, u = rsol.reconstruction(
                **optPB._fineParametrization.get_inputs4solution(
                    ampl*vecti[:, n]))
            pp.generatevtu(*optPB._fineParametrization.get_inputs4postprocVTU(
                'OPT8-fine-m%0.1d-%0.2d' % (n, i),  SOL.transpose(),
                nb_ref=ref_plot, Flag=Output))
        freqHistory = np.block([[freqHistory], [valsi]])
    return mi/m0


V0 = optPB.compute_volume(x0)


def vol(xV, gradV):
    if gradV.size > 0:
        gradV[:] = optPB.compute_gradVolume_AN(xV)/V0
    return optPB.compute_volume(xV)/V0 - 1.


minimize = nlopt.opt(nlopt.LD_SLSQP, nb_var)

minimize.set_min_objective(vibration)
minimize.add_equality_constraint(vol, 1e-5)

minimize.set_ftol_rel(1.0e-06)
minimize.set_xtol_rel(1.0e-06)
minimize.set_maxeval(200)

minimize.set_lower_bounds(0.25*np.ones(nb_var))
minimize.set_upper_bounds(3.00*np.ones(nb_var))


# --
# Run optimization
x = minimize.optimize(x0)

# Verify results
# Numerical reference result
x_ref = np.array([1.83742868, 1.65492556, 1.31414974, 0.66352119, 0.25, 0.25,
                  1.83109299, 1.6576408,  1.31434919, 0.66363835, 0.25, 0.25])

error = np.sqrt(sum((x-x_ref)**2.))

print(error)

if error > 1.e-5:
    sys.exit(-1)
else:
    sys.exit(0)
