# Copyright 2022 Thibaut Hirschler

# This file is part of Yeti.
#
# Yeti is free software: you can redistribute it and/or modify it under the terms 
# of the GNU Lesser General Public License as published by the Free Software 
# Foundation, either version 3 of the License, or (at your option) any later version.
#
# Yeti is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
# PURPOSE. See the GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License along 
# with Yeti. If not, see <https://www.gnu.org/licenses/>

#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Use a parameterized model.
Create a simple Reduced Order Model.
"""

import numpy as np
import scipy.sparse as sp

from preprocessing.igaparametrization import IGAparametrization
from stiffmtrx_elemstorage import sys_linmat_lindef_static as build_stiffmatrix
import reconstructionSOL as rsol
import postprocessing.postproc as pp

# Read data and create IGAparametrization object
modelIGA = IGAparametrization(filename='parametericPlateWithHole')

# Refine model
nb_deg = np.zeros((3,modelIGA._nb_patch),dtype=np.intp)
nb_ref = np.zeros((3,modelIGA._nb_patch),dtype=np.intp)
nb_deg[0,:] = 0
nb_deg[1,:] = 1
nb_ref[:2,:] = 2
modelIGA.refine(nb_ref, nb_deg)

# --
# Example of shape update.
# The design parameters are initially defined in the input file 'parametericPlateWithHole.inp'.
# Once the file read, these parameters are contained in the dictionary: see `modelIGA._design_parameters`.
designvar_name = list(modelIGA._design_parameters.keys())[0] # should be '<DESIGNVAR>'
designvar_value= modelIGA._design_parameters[designvar_name]

# a. Initial geometry
print('\nInitial Geometry')
print(' - design variable: %s = %f'%(designvar_name,designvar_value))
nb_ref_visu = np.array([2, 2, 0])         # Refinement for visu: [xi, eta, zeta]
output = np.array([False, False, False])  # Output type: [disp, stress, VM]
pp.generatevtu(*modelIGA.get_inputs4postprocVTU(
    'initialgeom',np.zeros((2,modelIGA._nb_cp),order='F'), nb_ref=nb_ref_visu, Flag=output))

# b. Shape modified geometry
# The shape modification is done in two steps:
# - step1: assign the new value of the design variable.
# - step2: call `modelIGA.shapeupdate()` to update the control points.
print('\nModified Geometry')
designvar_value = 0.25
modelIGA._design_parameters[designvar_name] = designvar_value # step1
modelIGA.shapeupdate() # step2
print(' - design variable: %s = %f'%(designvar_name,designvar_value))
nb_ref_visu = np.array([2, 2, 0])         # Refinement for visu: [xi, eta, zeta]
output = np.array([False, False, False])  # Output type: [disp, stress, VM]
pp.generatevtu(*modelIGA.get_inputs4postprocVTU(
    'modifiedgeom',np.zeros((2,modelIGA._nb_cp),order='F'), nb_ref=nb_ref_visu, Flag=output))



# --
# Let us try to build a simple Reduced Order Model by using this parameterized model.
print('\n\nBuild a simple Reduced Order Model.')

# Preliminaries: define a function to automatically build the FE systems.
ndof = modelIGA._nb_dof_free
idof = modelIGA._ind_dof_free[:ndof]-1
def buildsystem(designvar):
    '''Build the finite element operators for a given geometrical configuration.

    Parameters
    ----------
    designvar : float
        The value of the design variable, i.e. the radius of the hole.
        Should be between 0.05 and 0.95.
    
    Returns
    -------
    matK : sparse matrix
        Stiffness matrix.
    vecF : array of floats
        Load vector.
    '''
    if designvar<0.05 or designvar>0.95:
        raise ValueError('Input parameter should be a float in interval [0.05, 0.95].')

    modelIGA._design_parameters[designvar_name] = designvar
    modelIGA.shapeupdate()

    data, row, col, Fb = build_stiffmatrix(*modelIGA.get_inputs4system_elemStorage())
    Kside = sp.coo_matrix((data, (row, col)),
                      shape=(modelIGA._nb_dof_tot, modelIGA._nb_dof_tot),
                      dtype='float64').tocsc()
    Ktot = Kside + Kside.transpose()
    matK = Ktot[idof, :][:, idof]
    matK.sort_indices()
    vecF = Fb[idof]
    return matK,vecF

def refsolution(designvar):
    '''Perform the IGA for a given geometrical configuration.
    
    Parameters
    ----------
    designvar : float
        The value of the design variable, i.e. the radius of the hole.
        Should be between 0.05 and 0.95.
    
    Returns
    -------
    vecU : array of floats
        Solution vector.
    '''
    matK,vecF = buildsystem(designvar)
    vecU = sp.linalg.spsolve(matK,vecF)
    return vecU


# Step1: build snapshots of K and F.
print('- Generate the snapshots for K and F.')
from scipy.stats import qmc
rmin = 0.05
rmax = 0.95
sample = np.ravel((rmax-rmin)*qmc.LatinHypercube(1).random(n=100)+rmin)
sample = np.append(sample,[0.05,0.95])

snapshotsK = []
snapshotsF = []
for i in range(sample.size):
    matK,vecF = buildsystem(sample[i])
    snapshotsK.append([matK.data.copy()])
    snapshotsF.append([vecF.copy()])
snapshotsK = np.block(snapshotsK).T
snapshotsF = np.block(snapshotsF).T
print('  Done.')


# Step2: build a reduced basis for K and F.
# We do it by using a truncated SVD.
rtol = 1e-6 # tolerance for truncating the basis.
print('- Build reduced basis for K and F.')
print('  Tolerance %1.1e'%rtol)
from scipy.linalg import svd as scipySVD, solve as scipy_solve
from scipy import interpolate

U,s,Vh = scipySVD(snapshotsK)
nBasisK= np.count_nonzero(s/s.max()>rtol)
basisK = U[:,:nBasisK]*s[:nBasisK]
coefK  = [interpolate.interp1d(sample,Vh[I,:],'cubic') for I in range(nBasisK)]
print('  Size of the reduced basis for K:',nBasisK)

U,s,Vh = scipySVD(snapshotsF)
nBasisF = np.count_nonzero(s/s.max()>rtol)
basisF = U[:,:nBasisF]*s[:nBasisF]
coefF  = [interpolate.interp1d(sample,Vh[I,:],'cubic') for I in range(nBasisF)]
print('  Size of the reduced basis for F:',nBasisF)
print('  Done.')

# Check approximation
print('- Check accuracy')
# Define a function to automatically build the FE systems as given by the ROM.
indices = matK.indices.copy()
indptr = matK.indptr.copy()
def buildapproxsystem(designvar):
    '''Build the finite element operators as approximated by the ROM for a given geometrical configuration.
    
    Parameters
    ----------
    designvar : float
        The value of the design variable, i.e. the radius of the hole.
        Should be between 0.05 and 0.95.
    
    Returns
    -------
    matK : sparse matrix
        Stiffness matrix.
    vecF : array of floats
        Load vector.
    '''
    if designvar<0.05 or designvar>0.95:
        raise ValueError('Input parameter should be a float in interval [0.05, 0.95].')
    
    modelIGA._design_parameters[designvar_name] = designvar
    modelIGA.shapeupdate()
    
    alphaK = np.array([c(designvar) for c in coefK])
    dataK = basisK.dot(alphaK)
    matK = sp.csc_matrix((dataK,indices,indptr))
    
    alphaF = np.array([c(designvar) for c in coefF])
    vecF = basisF.dot(alphaF)
    return matK,vecF

designvar_value = 0.45
Kexact,Fexact = buildsystem(designvar_value)
Krom,From = buildapproxsystem(designvar_value)
print('  Parameter value used for the test: %f'%designvar_value)
print('  Error on K: %.2e' % (sp.linalg.norm(Kexact - Krom)/sp.linalg.norm(Kexact)))
print('  Error on F: %.2e' % (np.linalg.norm(Fexact - From)/np.linalg.norm(Fexact)))



# Step3: build a reduced basis for U.
print('- Generate the snapshots for U.')
from scipy.stats import qmc
rmin = 0.05
rmax = 0.95
sample = np.ravel((rmax-rmin)*qmc.LatinHypercube(1).random(n=100)+rmin)
sample = np.append(sample,[0.05,0.95])

snapshotsU = []
for i in range(sample.size):
    matK,vecF = buildapproxsystem(sample[i])
    vecU = sp.linalg.spsolve(matK,vecF)
    snapshotsU.append([vecU.copy()])
snapshotsU = np.block(snapshotsU).T
print('  Done.')

# Step4: build a reduced basis for U.
# We do it by using a truncated SVD.
rtolsol = 1e-4 # tolerance for truncating the basis.
print('- Build reduced basis for K and F.')
print('  Tolerance %1.1e'%rtolsol)

U,s,Vh = scipySVD(snapshotsU)
nBasisU= np.count_nonzero(s/s.max()>=rtolsol)
basisU = U[:,:nBasisU]
print('  Size of the reduced basis for U:',nBasisU)
print('  Done.')

def romsolution(designvar):
    '''Get the solution using the Reduced Order Model.
    
    Parameters
    ----------
    designvar : float
        The value of the design variable, i.e. the radius of the hole.
        Should be between 0.05 and 0.95.
    
    Returns
    -------
    vecU : array of floats
        Solution vector.
    '''
    if designvar<0.05 or designvar>0.95:
        raise ValueError('Input parameter should be a float in interval [0.05, 0.95].')
    
    matK,vecF = buildapproxsystem(designvar)
    
    matUKU = basisU.T.dot(matK.dot(basisU))
    vecUF = basisU.T.dot(vecF)
    alphaU = scipy_solve(matUKU,vecUF,assume_a='pos')
    
    vecU = basisU.dot(alphaU)
    return vecU

# Check approximation
print('\nCheck accuracy of ROM solution for different values of the shape parameter.')
sample = np.linspace(rmin,rmax,15)
errorROM = []
speedup = []
i = 0
nb_ref_visu = np.array([3, 3, 0])         # Refinement for visu: [xi, eta, zeta]
output = np.array([True, True, False])  # Output type: [disp, stress, VM]
for designvar_value in sample:
    
    Uiga = refsolution(designvar_value)
    SOL, u = rsol.reconstruction(**modelIGA.get_inputs4solution(Uiga))
    pp.generatevtu(*modelIGA.get_inputs4postprocVTU(
        'igasolution%02d'%i,SOL.T, nb_ref=nb_ref_visu, Flag=output))
    
    Urom = romsolution(designvar_value)
    SOL, u = rsol.reconstruction(**modelIGA.get_inputs4solution(Urom))
    pp.generatevtu(*modelIGA.get_inputs4postprocVTU(
        'romsolution%02d'%i,SOL.T, nb_ref=nb_ref_visu, Flag=output))
    SOL, u = rsol.reconstruction(**modelIGA.get_inputs4solution(Urom-Uiga))
    pp.generatevtu(*modelIGA.get_inputs4postprocVTU(
        'romerror%02d'%i,SOL.T, nb_ref=nb_ref_visu, Flag=np.array([True,False,False])))
    
    errorROM.append(np.linalg.norm(Uiga - Urom)/np.linalg.norm(Uiga))
    i += 1

print('\nShape var | ROM error')
for i in range(sample.size):
    print('  %f %.3e'%(sample[i],errorROM[i]))
