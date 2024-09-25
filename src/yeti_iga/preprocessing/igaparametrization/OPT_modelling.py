# Copyright 2018-2020 Thibaut Hirschler
# Copyright 2020 Marie Guerder
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

# -*- coding: utf-8 -*-

import os
import sys
import time
from copy import deepcopy
import logging

import numpy as np
from scipy.sparse import csc_matrix, linalg as sla
import scipy.sparse as sp

from ...stiffmtrx_elemstorage import sys_linmat_lindef_static as build_stiffmatrix
from ...massmtrx import build_cmassmatrix
from ...coupling.cplgmatrix import cplg_matrix
from ...coupling.cplgmatrix import cplg_matrixu5 as cplg_matrixU5
from ... import reconstructionSOL as rsol
from  ...postprocessing import postproc as pp

from ...optim import volume as vol
#from optim import gradvolume as gvol
from ...optim.gradadjointwork import gradlinelastwork_an
from ...optim import disp, stress, vibration, volume
from ...optim.gradcoupling import gradcplg_an


class OPTmodelling:
    """A class defining shape optimization problems.

    Parameters
    ----------
    useSPLU : bool
        Boolean indicating if optimization problem should use splu linear solver (default=True)
        WARNING : this parameter is not taken into account in all cases

    """

    def __init__(self, initialParametrization, nb_DesignVar,
                 fct_updateCoarseCoords,
                 nb_degreeElevationByDirection=np.array([0, 0, 0]),
                 nb_refinementByDirection=np.array([0, 0, 0]),
                 fct_dervShapeParam=None,
                 useSPLU=True):

        self._coarseParametrization = deepcopy(initialParametrization)
        self._fineParametrization = deepcopy(initialParametrization)

        self._transformationMatrices_coarse2fine = \
            self._fineParametrization.refine_and_getTransformationMatrices(
                nb_refinementByDirection, nb_degreeElevationByDirection)

        self._nb_DesignVar = int(nb_DesignVar)
        self._special_updateCoarseCoords = fct_updateCoarseCoords
        self._fct_dervShapeParametrization = fct_dervShapeParam

        # save
        self._initialCOORDS = np.copy(initialParametrization._COORDS)
        self._current_vectX = np.random.random(nb_DesignVar)
        self._save_sol_fine = np.zeros(self._fineParametrization._nb_dof_free)
        self._save_secondmembre = \
            np.zeros(self._fineParametrization._nb_dof_free)
        self._cp_by_var = self.get_indCPbyDesignVar()
        self._dir_by_var = self.get_dirCPbyDesignVar()
        self._elem_by_var = self.get_ElembyDesignVar()
        self._movable_patch = self.get_movablePatchesInd()

        self._useSPLU=useSPLU

        return None

    def get_indCPbyDesignVar(self):
        '''Get the indices of the C.P. linked with each Design Variable.'''
        self._updateCoarseCoords(np.zeros(self._nb_DesignVar))
        coords0 = np.copy(self._coarseParametrization._COORDS.transpose())
        cp_by_var = np.zeros(self._nb_DesignVar, dtype=object)
        for num_var in np.arange(0, self._nb_DesignVar):
            vect_test = np.zeros(self._nb_DesignVar)
            vect_test[num_var] = 1.
            self._updateCoarseCoords(vect_test)
            coords_test = self._coarseParametrization._COORDS.transpose()
            test = np.logical_not(np.isclose(coords_test, coords0))
            cp_by_var[num_var] = np.unique(np.where(test)[0])

        self._updateCoarseCoords(np.zeros(self._nb_DesignVar))

        return cp_by_var

    def get_dirCPbyDesignVar(self):
        '''Get moving direction of C.P. linked with each Design Variable.'''
        self._updateCoarseCoords(np.zeros(self._nb_DesignVar))
        coords0 = np.copy(self._coarseParametrization._COORDS)
        dir_by_var = np.zeros(self._nb_DesignVar, dtype=object)
        for num_var in np.arange(0, self._nb_DesignVar):
            vect_test = np.zeros(self._nb_DesignVar)
            vect_test[num_var] = 1.
            self._updateCoarseCoords(vect_test)
            coords_test = self._coarseParametrization._COORDS
            dir_by_var[num_var] = \
                (coords_test-coords0)[:, self._cp_by_var[num_var]]
        self._updateCoarseCoords(np.zeros(self._nb_DesignVar))

        return dir_by_var

    def get_ElembyDesignVar(self):
        '''Get the list of Elements linked with each Design Variable.'''
        self._updateCoarseCoords(np.zeros(self._nb_DesignVar))
        self._updateFineCoords()
        coords0 = np.copy(self._fineParametrization._COORDS.transpose())
        offset = np.cumsum(
            np.concatenate(([0], self._fineParametrization._elementsByPatch)))
        elem_by_var = np.zeros(self._nb_DesignVar, dtype=object)
        for num_var in np.arange(0, self._nb_DesignVar):
            vect_test = np.zeros(self._nb_DesignVar)
            vect_test[num_var] = 1.
            self._updateCoarseCoords(vect_test)
            self._updateFineCoords()
            coords_test = self._fineParametrization._COORDS.transpose()
            test = np.logical_not(np.isclose(coords_test, coords0))
            list_cps = np.unique(np.where(test)[0]) + 1
            activeElem = np.zeros(self._fineParametrization._nb_elem,
                                  dtype=np.intp)
            for patch in np.arange(0, self._fineParametrization._nb_patch):
                for cp in list_cps:
                    index = np.where(cp ==
                                     self._fineParametrization._IEN[patch])[0]
                    activeElem[index + offset[patch]] = 1
            elem_by_var[num_var] = activeElem
        self._updateCoarseCoords(np.zeros(self._nb_DesignVar))
        self._updateFineCoords()

        return elem_by_var

    def get_movablePatchesInd(self):
        if not hasattr(self, '_cp_by_var'):
            self._cp_by_var = self.get_indCPbyDesignVar()
        cpall = np.concatenate(self._cp_by_var)
        ipatch = []
        for p in np.arange(self._coarseParametrization._nb_patch):
            if np.any(np.isin(self._coarseParametrization._indCPbyPatch[p] - 1,
                              cpall)):
                ipatch.append(p)

        return np.array(ipatch, dtype=np.intp)

    def differentiateShapeParametrization(self, vectX, numVar):
        # size (n_cp x 3)
        dShape_coarse = np.zeros_like(self._coarseParametrization._COORDS)
        if self._fct_dervShapeParametrization is None:
            # linear shape param
            dShape_coarse[:, self._cp_by_var[numVar]] = \
                self._dir_by_var[numVar]
        else:
            dShape_coarse = self._fct_dervShapeParametrization(
                    self._initialCOORDS, self._coarseParametrization, vectX,
                    numVar)
        return dShape_coarse

    def _updateCoarseCoords(self, vectX):
        self._special_updateCoarseCoords(self._initialCOORDS,
                                         self._coarseParametrization, vectX)
        return None

    def _updateFineCoords(self):
        self._fineParametrization._COORDS[:,:] = self._mapCoarse2Fine(
            self._coarseParametrization._COORDS.T).T
        return None

    def _mapCoarse2Fine(self, coarseQ):
        vectWeight_coarse = self._coarseParametrization._vectWeight
        vectWeight_fine = self._fineParametrization._vectWeight
        fineQ = np.zeros(self._fineParametrization._COORDS.shape[::-1])
        indCPbyPatch_coarse = self._coarseParametrization._indCPbyPatch
        indCPbyPatch_fine = self._fineParametrization._indCPbyPatch

        for num_patch in range(0, self._coarseParametrization._nb_patch):

            nb_cp_thisPatch = np.size(indCPbyPatch_coarse[num_patch], 0)
            Pwi_thisPatch = np.zeros((nb_cp_thisPatch, 3))
            Pwi_thisPatch[:, 0] = \
                coarseQ[indCPbyPatch_coarse[num_patch] - 1, 0] * \
                vectWeight_coarse[indCPbyPatch_coarse[num_patch] - 1]
            Pwi_thisPatch[:, 1] = \
                coarseQ[indCPbyPatch_coarse[num_patch] - 1, 1] * \
                vectWeight_coarse[indCPbyPatch_coarse[num_patch] - 1]
            Pwi_thisPatch[:, 2] = \
                coarseQ[indCPbyPatch_coarse[num_patch] - 1, 2] * \
                vectWeight_coarse[indCPbyPatch_coarse[num_patch] - 1]
            #Pwi_thisPatch[:, 3] = \
            #    vectWeight_coarse[indCPbyPatch_coarse[num_patch] - 1]

            Qwi_xi_thisPatch = \
                self._transformationMatrices_coarse2fine[num_patch][0] * \
                Pwi_thisPatch
            Qwi_xi_eta_thisPatch = \
                self._transformationMatrices_coarse2fine[num_patch][1] * \
                Qwi_xi_thisPatch
            Qwi_thisPatch = \
                self._transformationMatrices_coarse2fine[num_patch][2] * \
                Qwi_xi_eta_thisPatch

            nb_cp_fine_thisPatch = np.size(Qwi_thisPatch, 0)
            Qi_thisPatch = np.zeros((nb_cp_fine_thisPatch, 3))
            wi_thisPatch = vectWeight_fine[indCPbyPatch_fine[num_patch]-1]
            Qi_thisPatch[:, 0] = Qwi_thisPatch[:, 0] / wi_thisPatch[:]
            Qi_thisPatch[:, 1] = Qwi_thisPatch[:, 1] / wi_thisPatch[:]
            Qi_thisPatch[:, 2] = Qwi_thisPatch[:, 2] / wi_thisPatch[:]

            fineQ[indCPbyPatch_fine[num_patch]-1, :] = Qi_thisPatch[:, :]

        return fineQ

    def _transposeFine2Coarse(self, fineQ):
        vectWeight_coarse = self._coarseParametrization._vectWeight
        vectWeight_fine = self._fineParametrization._vectWeight
        coarseQ = np.zeros_like(self._coarseParametrization._COORDS)
        indCPbyPatch_coarse = self._coarseParametrization._indCPbyPatch
        indCPbyPatch_fine = self._fineParametrization._indCPbyPatch

        for num_patch in range(0, self._fineParametrization._nb_patch):
            ind_cp_fine_thisPatch = indCPbyPatch_fine[num_patch]-1
            vectWeight_fine_thisPatch = vectWeight_fine[ind_cp_fine_thisPatch]
            nb_cp_thisPatch = np.size(ind_cp_fine_thisPatch, 0)
            Qwi_thisPatch = np.zeros((3,nb_cp_thisPatch))
            Qwi_thisPatch[0,:] = fineQ[0,ind_cp_fine_thisPatch] / vectWeight_fine_thisPatch
            Qwi_thisPatch[1,:] = fineQ[1,ind_cp_fine_thisPatch] / vectWeight_fine_thisPatch
            Qwi_thisPatch[2,:] = fineQ[2,ind_cp_fine_thisPatch] / vectWeight_fine_thisPatch

            Pwi_zeta_thisPatch = \
                Qwi_thisPatch*self._transformationMatrices_coarse2fine[num_patch][2]
            Pwi_eta_zeta_thisPatch = \
                Pwi_zeta_thisPatch*self._transformationMatrices_coarse2fine[num_patch][1]
            Pwi_thisPatch = \
                Pwi_eta_zeta_thisPatch*self._transformationMatrices_coarse2fine[num_patch][0]

            nb_cp_coarse_thisPatch = np.size(Pwi_thisPatch, 1)
            ind_cp_coarse_thisPatch = indCPbyPatch_coarse[num_patch]-1
            vectWeight_coarse_thisPatch = vectWeight_coarse[ind_cp_coarse_thisPatch]
            Pi_thisPatch = np.zeros((3,nb_cp_coarse_thisPatch))
            Pi_thisPatch[0,:] = Pwi_thisPatch[0,:] * vectWeight_coarse_thisPatch
            Pi_thisPatch[1,:] = Pwi_thisPatch[1,:] * vectWeight_coarse_thisPatch
            Pi_thisPatch[2,:] = Pwi_thisPatch[2,:] * vectWeight_coarse_thisPatch

            coarseQ[:,ind_cp_coarse_thisPatch] = Pi_thisPatch[:, :]
        return coarseQ


    # --
    # ANALYSIS
    def get_LinElastSystem(self, vectX, activeelem=None):
        if np.all(vectX == self._current_vectX):
            return self._save_raideur,self._save_secondmembre
        self._updateCoarseCoords(vectX)
        self._updateFineCoords()
        ndof = self._fineParametrization._nb_dof_free
        idof = self._fineParametrization._ind_dof_free[:ndof] - 1
        data, row, col, Fb = build_stiffmatrix(
            *self._fineParametrization.get_inputs4system_elemStorage(
                activeElem=activeelem))
        Kside = sp.coo_matrix((data, (row, col)),
                              shape=(self._fineParametrization._nb_dof_tot,
                                     self._fineParametrization._nb_dof_tot),
                              dtype='float64').tocsc()
        Ktot = Kside + Kside.transpose()
        K2solve = Ktot[idof, :][:, idof]
        del Kside, data, row, col, Ktot
        if np.any(self._fineParametrization._ELT_TYPE == 'U00'):
            # Coupling
            Cdata, Crow, Ccol = \
                cplg_matrix(*self._fineParametrization.get_inputs4cplgmatrix())
            Cside = sp.coo_matrix(
                (Cdata, (Crow, Ccol)),
                shape=(self._fineParametrization._nb_dof_tot,
                       self._fineParametrization._nb_dof_tot),
                dtype='float64').tocsc()
            Ctot = Cside + Cside.transpose()
            C2solve = Ctot[idof, :][:, idof]
            del Cdata, Crow, Ccol, Cside, Ctot
            K2solve = K2solve + C2solve
        if np.any(self._fineParametrization._ELT_TYPE == 'U5'):
            # Coupling U5
            Cdata, Crow, Ccol = \
                cplg_matrixU5(*self._fineParametrization.get_inputs4cplgmatrixU5())
            Cside = sp.coo_matrix(
                (Cdata, (Crow, Ccol)),
                shape=(self._fineParametrization._nb_dof_tot,
                       self._fineParametrization._nb_dof_tot),
                dtype='float64').tocsc()
            Ctot = Cside + Cside.transpose()
            C2solve = Ctot[idof, :][:, idof]
            del Cdata, Crow, Ccol, Cside, Ctot
            K2solve = K2solve + C2solve

        return K2solve, Fb[idof]

    def compute_analysisSol(self, vectX):
        if np.all(vectX == self._current_vectX):
            return self._save_sol_fine
        else:
            print(' Linear elasticity analysis...')
            K, F = self.get_LinElastSystem(vectX)

            if self._useSPLU:
                LU = sp.linalg.splu(K)
                sol_fine = LU.solve(F)
            else:
                sol_fine = sp.linalg.spsolve(K, F)

            self._save_sol_fine = sol_fine
            self._save_secondmembre = F
            self._save_raideur = K
            if self._useSPLU:
                self._save_LU = LU
            self._current_vectX[:] = vectX[:]

            return sol_fine

    def compute_LinElastAdjWork(self,vectX,stateSOL,adjointSOL,activeelem=None,
                                return_internalWork=True,return_externalWork=True):
        if not (return_internalWork or return_externalWork):
            return None
        K,F = self.get_LinElastSystem(vectX, activeelem=activeelem)
        self._updateCoarseCoords(self._current_vectX)
        self._updateFineCoords()
        if return_internalWork:
            Wint = -np.dot(adjointSOL.T,K.dot(stateSOL))
            if not return_externalWork:
                return Wint
        if return_externalWork:
            Wext = np.dot(adjointSOL.T,F)
            if not return_internalWork:
                return Wext
        return Wint,Wext

    def get_NaturalFrqSystem(self, vectX):
        self._updateCoarseCoords(vectX)
        self._updateFineCoords()
        ndof = self._fineParametrization._nb_dof_free
        idof = self._fineParametrization._ind_dof_free[: ndof] - 1

        data, row, col, Fb = build_stiffmatrix(
            *self._fineParametrization.get_inputs4system_elemStorage())
        Kside = sp.coo_matrix((data, (row, col)),
                              shape=(self._fineParametrization._nb_dof_tot,
                                     self._fineParametrization._nb_dof_tot),
                              dtype='float64').tocsc()
        Ktot = Kside + Kside.transpose()
        K2solve = Ktot[idof, :][:, idof]
        del Kside, data, row, col, Ktot

        data, row, col = build_cmassmatrix(
            *self._fineParametrization.get_inputs4massmat())
        Mside = sp.coo_matrix((data, (row, col)),
                              shape=(self._fineParametrization._nb_dof_tot,
                                     self._fineParametrization._nb_dof_tot),
                              dtype='float64').tocsc()
        Mtot = Mside + Mside.transpose()
        M2solve = Mtot[idof, :][:, idof]
        del Mside, data, row, col, Mtot

        return K2solve, M2solve

    def compute_vibrationMode(self, vectX, nb_frq=1):
        if np.all(vectX == self._current_vectX) and nb_frq<= self._save_vals.size:
            return self._save_vals[:nb_frq], self._save_vect[:,:nb_frq]
        else:
            print(' Natural frequency analysis...')
            K2solve, M2solve = self.get_NaturalFrqSystem(vectX)
            vals, vecs = sp.linalg.eigsh(K2solve, k=nb_frq, M=M2solve,
                                         sigma=0.)
            self._save_vals = vals.copy()
            self._save_vect = vecs.copy()
            self._save_raideur = K2solve
            self._save_masse = M2solve
            self._current_vectX[:] = vectX[:]
        return vals, vecs

    def compute_NaturalFrqWork(self,vectX,vals,vecs):
        K,M = self.get_NaturalFrqSystem(vectX)
        self._updateCoarseCoords(self._current_vectX)
        self._updateFineCoords()
        Wint = []; Wcin = []
        for ifrq in np.arange(vals.size):
            veci = vecs[:,ifrq]
            Wint.append(-np.dot(veci,K.dot(veci)))
            Wcin.append(-np.dot(veci,M.dot(veci))*vals[ifrq])
        return np.array(Wcin)-np.array(Wint)


    # --
    # SENSITIVITY
    #   Baseline function for response functions and sensitivity analysis
    def _sensitivityFD(self,fctF,vectX,eps=1.e-6,centerFD=False):
        '''Baseline function for Finite Difference Approximation of the Sensitivities'''
        f0 = fctF(vectX)
        nbVar = self._nb_DesignVar
        sizeF = np.size(f0)
        if sizeF==1:
            gradF = np.zeros(nbVar,dtype=np.float64)
        else:
            gradF = np.zeros((nbVar,sizeF),dtype=np.float64)
        vectX_thisVar = np.zeros(nbVar)
        vectX_save = vectX.copy()
        for numVar in np.arange(nbVar):
            sys.stdout.write(' numVar:%3i/%i\r' % (numVar+1, self._nb_DesignVar))
            sys.stdout.flush()
            vectX_thisVar[:] = vectX_save[:]
            vectX_thisVar[numVar] += eps
            fL = fctF(vectX_thisVar)
            if centerFD is False:
                gradF[numVar] = (fL-f0)/eps
            else:
                vectX_thisVar[:] = vectX_save[:]
                vectX_thisVar[numVar] -= eps
                fR = fctF(vectX_thisVar)
                gradF[numVar] = (fL-fR)/2./eps
        print(' Done.          ')

        return gradF

    def _propagateSensitivityAnalysis(self,vectX,gradFine):
        gradDV = np.zeros(self._nb_DesignVar)
        gradCoarse = self._transposeFine2Coarse(gradFine)
        for numVar in np.arange(self._nb_DesignVar):
            dShape = self.differentiateShapeParametrization(vectX, numVar)
            gradDV[numVar] = np.tensordot(dShape,gradCoarse,axes=2)
        return gradDV

    def compute_dervLinElastAdjWork(self,vectX,U0,UA,computeWint=True,computeWext=True):
        self._updateCoarseCoords(vectX)
        self._updateFineCoords()
        nb_cp = self._fineParametrization._nb_cp
        mcrd  = self._fineParametrization._mcrd
        if UA.ndim == 1:
            UA = np.vstack(UA)
        nbAdj = np.size(UA,1)

        SOL, u = rsol.reconstruction(
            **self._fineParametrization.get_inputs4solution(U0))
        ADJ = np.zeros((nbAdj, mcrd, nb_cp), dtype=np.float64)
        for iA in np.arange(nbAdj):
            adjiA, u = rsol.reconstruction(
                **self._fineParametrization.get_inputs4solution(UA[:, iA]))
            ADJ[iA, :, :] = adjiA[:, :].T

        activeElem = np.zeros(self._fineParametrization._nb_elem,
                              dtype=np.intp)
        activeElem[np.where(np.sum(self._elem_by_var) > 0)[0]] = 1
        activeDir = np.zeros(3, dtype=np.intp)
        savemax = np.zeros(3)
        for var in self._dir_by_var:
            savemax = np.maximum(savemax, np.max(np.abs(var), axis=1))
        activeDir[np.where(savemax > 0)[0]] = 1


        dervWint,dervWext = gradlinelastwork_an(
            *self._fineParametrization._get_inputs4gradTotalwork(
                SOL.T, ADJ, activeElem=activeElem, activeDir=activeDir,computeWint=computeWint,
                computeWext=computeWext))

        self._updateCoarseCoords(self._current_vectX)
        self._updateFineCoords()
        return dervWint,dervWext

    def _set_greville(self):
        from ...fitting.interpolate import grevilleabscissae
        pts = []
        for ipatch in np.arange(self._fineParametrization._nb_patch):
            u = self._fineParametrization._Ukv[ipatch]
            p = self._fineParametrization._Jpqr[:,ipatch]
            n = self._fineParametrization._Nkv[:,ipatch] - p-1
            ptspatch = [np.zeros(1)]*3
            for i in np.arange(self._fineParametrization._dim[ipatch]):
                ptspatch[i] = grevilleabscissae(uknot=u[i],p=p[i],n=n[i],nknots=p[i]+n[i]+1)
            pts.append(
                np.block([[np.tile(ptspatch[0],ptspatch[1].size*ptspatch[2].size)],
                          [np.tile(np.repeat(ptspatch[1],ptspatch[0].size),ptspatch[2].size)],
                          [np.repeat(ptspatch[2],ptspatch[0].size*ptspatch[1].size)]]) )
        self._grevillepts = pts
        return None

    def _get_greville(self):
        if not hasattr(self,'_grevillepts'):
            self._set_greville()
        return self._grevillepts

    # --
    # RESPONSE FUNCTIONS
    # 1. Volume/area
    def compute_volume(self, vectX, listpatch=None):
        logging.info("Compute Volume")
        self._updateCoarseCoords(vectX)
        self._updateFineCoords()
        V = volume.computevolume(
            *self._fineParametrization.get_inputs4area(
                activepatch=listpatch))
        return V

    def compute_gradVolume_DF(self, vectX, eps=1.e-6, centerFD=False, listpatch=None):
        print("Compute FD gradient of the volume")
        fctF = lambda x: self.compute_volume(x,listpatch=listpatch)
        gradV = self._sensitivityFD(fctF,vectX,eps=eps,centerFD=centerFD)
        return gradV

    def compute_gradVolume_AN(self, vectX, listpatch=None):
        print("Compute AN gradient of the volume")
        self._updateCoarseCoords(vectX)
        self._updateFineCoords()
        dVol_fine = volume.computegradvolume(
                *self._fineParametrization.get_inputs4area(
                    activepatch=listpatch))
        gradV = self._propagateSensitivityAnalysis(vectX,dVol_fine)
        return gradV

    # 2. Compliance
    def compute_compliance_discrete(self, vectX):
        print("Compute compliance")
        U0 = self.compute_analysisSol(vectX)
        F0 = self._save_secondmembre
        comp = 0.5 * np.dot(U0, F0)
        self._saveComp = comp
        return comp

    def compute_gradCompliance_FD(self, vectX, eps=1.e-6, centerFD=False):
        print("Compute FD gradient of the compliance")
        fctF = lambda x: self.compute_compliance_discrete(x)
        gradF = self._sensitivityFD(fctF,vectX,eps=eps,centerFD=centerFD)
        return gradF

    def compute_gradCompliance_semiAN(self, vectX, eps=1.e-6, centerFD=False):
        print("Compute semiAN gradient of the compliance")
        U0 = self.compute_analysisSol(vectX)
        UA = 0.5*U0
        def workW(x):
            Wint,Wext = self.compute_LinElastAdjWork(x,U0,UA)
            return Wint + 2.0*Wext
        gradW = self._sensitivityFD(workW,vectX,eps=eps,centerFD=centerFD)
        return gradW

    def compute_gradCompliance_AN(self, vectX):
        print("Compute AN gradient of the compliance")
        U0 = self.compute_analysisSol(vectX)
        UA = 0.5*U0
        dervWint,dervWext = self.compute_dervLinElastAdjWork(vectX,U0,UA)
        dervW = dervWint[0]+2*dervWext[0]
        gradC = self._propagateSensitivityAnalysis(vectX,dervW)
        return gradC

    # 3. Displacement
    def compute_displacement(self, vectX, xi, numpatch=1):
        print("Compute displacement")
        U0 = self.compute_analysisSol(vectX)
        SOL, U = rsol.reconstruction(
            **self._fineParametrization.get_inputs4solution(U0))
        disp = pp.evaldisp(
            *self._fineParametrization.get_inputs4evaldisp(SOL.T, xi,
                                                           numpatch=numpatch))
        return disp

    def compute_gradDisplacement_FD(self, vectX, xi, eps=1.e-6, centerFD=False):
        print("Compute FD gradient of the displacement")
        fctF = lambda x: self.compute_displacement(x,xi)
        gradF = self._sensitivityFD(fctF,vectX,eps=eps,centerFD=centerFD)
        return gradF

    def compute_gradDisplacement_semiAN(self, vectX, xi, eps=1.e-6, centerFD=False):
        print("Compute semiAN gradient of the displacement")
        U0 = self.compute_analysisSol(vectX) # state solution
        ndof = self._fineParametrization._nb_dof_free
        idof = self._fineParametrization._ind_dof_free[:ndof]-1
        FA = disp.adjointrhsdisp( *self._fineParametrization.get_inputs4adjointdisp(xi) ).T
        UA = self._save_LU.solve(FA[idof]) # adjoint solution
        workW = lambda x: np.sum(self.compute_LinElastAdjWork(x,U0,UA),axis=0)
        gradW = self._sensitivityFD(workW,vectX,eps=eps,centerFD=centerFD)
        return gradW

    def compute_gradDisplacement_AN(self, vectX, xi):
        print("Compute AN gradient of the displacement")
        U0 = self.compute_analysisSol(vectX) # state solution
        ndof = self._fineParametrization._nb_dof_free
        idof = self._fineParametrization._ind_dof_free[:ndof]-1
        FA = disp.adjointrhsdisp( *self._fineParametrization.get_inputs4adjointdisp(xi) ).T
        UA = self._save_LU.solve(FA[idof]) # adjoint solution
        dervWint,derWext = self.compute_dervLinElastAdjWork(vectX,U0,UA)
        dervW = dervWint+derWext
        nDV = self._nb_DesignVar;nA = np.size(UA,1)
        gradD = np.zeros((nDV,nA),dtype=np.float64)
        for i in np.arange(nA):
            gradD[:,i] = self._propagateSensitivityAnalysis(vectX,dervW[i])
        return gradD

    def compute_displacementAggreg(self,vectX,pnorm=10,ptseval=None):
        print("Compute maximal displacement (approximated through P-Norm)")
        if ptseval is None:
            ptseval = self._get_greville()
        U0 = self.compute_analysisSol(vectX)
        SOL,U = rsol.reconstruction(
            **self._fineParametrization.get_inputs4solution(U0))
        dispAggreg = disp.pnormdispmagn(
            *self._fineParametrization.get_inputs4dispaggregDiscrete(SOL.T,ptseval,pnorm=pnorm) )
        return dispAggreg

    def compute_gradDisplacementAggreg_FD(self, vectX, pnorm=10, ptseval=None,
                                          eps=1.e-6, centerFD=False):
        print("Compute FD gradient of the maximal displacement")
        fctF = lambda x: self.compute_displacementAggreg(x,pnorm=pnorm,ptseval=ptseval)
        gradF = self._sensitivityFD(fctF,vectX,eps=eps,centerFD=centerFD)
        return gradF

    def compute_gradDisplacementAggreg_semiAN(self, vectX, pnorm=10, ptseval=None,
                                              eps=1.e-6, centerFD=False):
        print("Compute semiAN gradient of the maximal displacement")
        if ptseval is None:
            ptseval = self._get_greville()
        U0 = self.compute_analysisSol(vectX) # state solution
        SOL,U = rsol.reconstruction(
            **self._fineParametrization.get_inputs4solution(U0))
        ndof = self._fineParametrization._nb_dof_free
        idof = self._fineParametrization._ind_dof_free[:ndof]-1
        FA = disp.adjointrhspnormdispmagn(
            *self._fineParametrization.get_inputs4dispaggregDiscrete(SOL.T,ptseval,pnorm=pnorm) )
        UA = self._save_LU.solve(FA[idof]) # adjoint solution
        workW = lambda x: np.sum(self.compute_LinElastAdjWork(x,U0,UA),axis=0)
        gradW = self._sensitivityFD(workW,vectX,eps=eps,centerFD=centerFD)
        return gradW

    def compute_gradDisplacementAggreg_AN(self, vectX, pnorm=10, ptseval=None):
        print("Compute AN gradient of the maximal displacement")
        if ptseval is None:
            ptseval = self._get_greville()
        U0 = self.compute_analysisSol(vectX) # state solution
        SOL,U = rsol.reconstruction(
            **self._fineParametrization.get_inputs4solution(U0))
        ndof = self._fineParametrization._nb_dof_free
        idof = self._fineParametrization._ind_dof_free[:ndof]-1
        FA = disp.adjointrhspnormdispmagn(
            *self._fineParametrization.get_inputs4dispaggregDiscrete(SOL.T,ptseval,pnorm=pnorm) )
        UA = self._save_LU.solve(FA[idof]) # adjoint solution
        dervWint,derWext = self.compute_dervLinElastAdjWork(vectX,U0,UA)
        dervW = dervWint+derWext
        gradD = self._propagateSensitivityAnalysis(vectX,dervW[0])
        return gradD

    # 4. Stress
    def compute_stressAggreg(self, vectX, pnorm=10, ptseval=None):
        print("Compute maximal stress (approximated through P-Norm)")
        if ptseval is None:
            ptseval = self._get_greville()
        U0 = self.compute_analysisSol(vectX)
        SOL, U = rsol.reconstruction(
            **self._fineParametrization.get_inputs4solution(U0))
        stressAggreg = stress.pnormstress(
            *self._fineParametrization.get_inputs4stressaggregDiscrete(
                SOL.T,ptseval,pnorm=pnorm) )
        return stressAggreg

    def compute_gradStressAggreg_FD(self, vectX, pnorm=10, ptseval=None,
                                     eps=1.e-6, centerFD=False):
        print("Compute FD gradient of the maximal stress")
        fctF = lambda x: self.compute_stressAggreg(x,pnorm=pnorm,ptseval=ptseval)
        gradF = self._sensitivityFD(fctF,vectX,eps=eps,centerFD=centerFD)
        return gradF

    def compute_gradStressAggreg_semiAN(self, vectX, pnorm=10, ptseval=None,
                                           eps=1.e-6, centerFD=False):
        print("Compute semiAN gradient of the maximal stress")
        if ptseval is None:
            ptseval = self._get_greville()
        U0 = self.compute_analysisSol(vectX) # state solution
        SOL,U = rsol.reconstruction(
            **self._fineParametrization.get_inputs4solution(U0))
        ndof = self._fineParametrization._nb_dof_free
        idof = self._fineParametrization._ind_dof_free[:ndof]-1
        FA = stress.adjointrhspnormstress(
            *self._fineParametrization.get_inputs4stressaggregDiscrete(
                SOL.T,ptseval,pnorm=pnorm) ).T
        UA = self._save_LU.solve(FA[idof]) # adjoint solution

        workW = lambda x: np.sum(self.compute_LinElastAdjWork(x,U0,UA),axis=0)
        def partW(x):
            self._updateCoarseCoords(x)
            self._updateFineCoords()
            vecstress = stress.pnormstress(
                *self._fineParametrization.get_inputs4stressaggregDiscrete(
                    SOL.T,ptseval,pnorm=pnorm) )
            self._updateCoarseCoords(self._current_vectX)
            self._updateFineCoords()
            return vecstress
        gradW = self._sensitivityFD(workW,vectX,eps=eps,centerFD=centerFD)
        gradW+= self._sensitivityFD(partW,vectX,eps=eps,centerFD=centerFD)
        return gradW

    def compute_gradStressAggreg_AN(self, vectX, pnorm=10, ptseval=None):
        print("Compute semiAN gradient of the maximal stress")
        if ptseval is None:
            ptseval = self._get_greville()
        U0 = self.compute_analysisSol(vectX) # state solution
        SOL,U = rsol.reconstruction(
            **self._fineParametrization.get_inputs4solution(U0))
        ndof = self._fineParametrization._nb_dof_free
        idof = self._fineParametrization._ind_dof_free[:ndof]-1
        FA = stress.adjointrhspnormstress(
            *self._fineParametrization.get_inputs4stressaggregDiscrete(
                SOL.T,ptseval,pnorm=pnorm) ).T
        UA = self._save_LU.solve(FA[idof]) # adjoint solution

        dervWint,derWext = self.compute_dervLinElastAdjWork(vectX,U0,UA)
        dervW = dervWint+derWext
        dervW+= stress.partialdervcppnormstress(
            *self._fineParametrization.get_inputs4stressaggregDiscrete(
                SOL.T,ptseval,pnorm=pnorm) )

        nDV = self._nb_DesignVar;nA = np.size(UA,1)
        gradW = np.zeros((nDV,nA))
        for i in np.arange(nA):
            gradW[:,i] = self._propagateSensitivityAnalysis(vectX,dervW[i])
        return gradW

    def compute_vonmisesAggreg(self, vectX, pnorm=10, ptseval=None):
        print("Compute maximal von-mises stress (approximated through P-Norm)")
        if ptseval is None:
            ptseval = self._get_greville()
        U0 = self.compute_analysisSol(vectX)
        SOL, U = rsol.reconstruction(
            **self._fineParametrization.get_inputs4solution(U0))
        vmAggreg = stress.pnormvm(
            *self._fineParametrization.get_inputs4stressaggregDiscrete(
                SOL.T,ptseval,pnorm=pnorm) )
        return vmAggreg

    def compute_gradVonMisesAggreg_FD(self, vectX, pnorm=10, ptseval=None,
                                       eps=1.e-6, centerFD=False):
        print("Compute FD gradient of the maximal von-mises stress")
        fctF = lambda x: self.compute_vonmisesAggreg(x,pnorm=pnorm,ptseval=ptseval)
        gradF = self._sensitivityFD(fctF,vectX,eps=eps,centerFD=centerFD)
        return gradF

    def compute_gradVonMisesAggreg_semiAN(self, vectX, pnorm=10, ptseval=None,
                                           eps=1.e-6, centerFD=False):
        print("Compute semiAN gradient of the maximal von-mises stress")
        if ptseval is None:
            ptseval = self._get_greville()
        U0 = self.compute_analysisSol(vectX) # state solution
        SOL,U = rsol.reconstruction(
            **self._fineParametrization.get_inputs4solution(U0))
        ndof = self._fineParametrization._nb_dof_free
        idof = self._fineParametrization._ind_dof_free[:ndof]-1
        FA = stress.adjointrhspnormvm(
            *self._fineParametrization.get_inputs4stressaggregDiscrete(
                SOL.T,ptseval,pnorm=pnorm) )
        UA = self._save_LU.solve(FA[idof]) # adjoint solution

        workW = lambda x: np.sum(self.compute_LinElastAdjWork(x,U0,UA),axis=0)
        def partW(x):
            self._updateCoarseCoords(x)
            self._updateFineCoords()
            vm = stress.pnormvm(
                *self._fineParametrization.get_inputs4stressaggregDiscrete(
                    SOL.T,ptseval,pnorm=pnorm) )
            self._updateCoarseCoords(self._current_vectX)
            self._updateFineCoords()
            return vm
        gradW = self._sensitivityFD(workW,vectX,eps=eps,centerFD=centerFD)
        gradW+= self._sensitivityFD(partW,vectX,eps=eps,centerFD=centerFD)
        return gradW


    def compute_gradVonMisesAggreg_AN(self, vectX, pnorm=10, ptseval=None):
        print("Compute AN gradient of the maximal von-mises stress")
        if ptseval is None:
            ptseval = self._get_greville()
        U0 = self.compute_analysisSol(vectX) # state solution
        SOL,U = rsol.reconstruction(
            **self._fineParametrization.get_inputs4solution(U0))
        ndof = self._fineParametrization._nb_dof_free
        idof = self._fineParametrization._ind_dof_free[:ndof]-1
        FA = stress.adjointrhspnormvm(
            *self._fineParametrization.get_inputs4stressaggregDiscrete(
                SOL.T,ptseval,pnorm=pnorm) )
        UA = self._save_LU.solve(FA[idof]) # adjoint solution

        dervWint,derWext = self.compute_dervLinElastAdjWork(vectX,U0,UA)
        dervW = dervWint+derWext
        dervW+= stress.partialdervcppnormvm(
            *self._fineParametrization.get_inputs4stressaggregDiscrete(
                SOL.T,ptseval,pnorm=pnorm) )
        gradW = self._propagateSensitivityAnalysis(vectX,dervW[0])
        return gradW


    # 5. Natural frequencies
    def compute_gradVibration_FD(self, vectX, nb_frq=1, eps=1.e-6, centerFD=False):
        print("Compute FD gradient of the natural frequencies")
        fctF = lambda x: self.compute_vibrationMode(x,nb_frq=nb_frq)[0]
        gradF = self._sensitivityFD(fctF,vectX,eps=eps,centerFD=centerFD)
        return gradF

    def compute_gradVibration_semiAN(self,vectX, nb_frq=1, eps=1.e-6, centerFD=False):
        print("Compute semiAN gradient of the natural frequencies")
        vals0, vecs0 = self.compute_vibrationMode(vectX, nb_frq=nb_frq)
        workW = lambda x: self.compute_NaturalFrqWork(x,vals0,vecs0)
        gradW = self._sensitivityFD(workW,vectX,eps=eps,centerFD=centerFD)
        return gradW

    def compute_gradVibration_AN(self, vectX, nb_frq=1):
        print("Compute AN gradient of the natural frequencies")
        nbVar = self._nb_DesignVar
        mcrd = self._fineParametrization._mcrd
        nb_cp = self._fineParametrization._nb_cp
        w0, V0 = self.compute_vibrationMode(vectX, nb_frq=nb_frq)
        VECT = np.zeros((nb_frq, mcrd, nb_cp), dtype=np.float64)
        for iV in np.arange(nb_frq):
            SOL, u = rsol.reconstruction(
                **self._fineParametrization.get_inputs4solution(V0[:, iV]))
            VECT[iV] = SOL[:, :].T
        ndof = self._fineParametrization._nb_dof_free
        idof = self._fineParametrization._ind_dof_free[: ndof] - 1
        activeElem = np.zeros(self._fineParametrization._nb_elem,
                              dtype=np.intp)
        activeElem[np.where(np.sum(self._elem_by_var) > 0)[0]] = 1
        activeDir = np.zeros(3, dtype=np.intp)
        savemax = np.zeros(3)
        for var in self._dir_by_var:
            savemax = np.maximum(savemax, np.max(np.abs(var), axis=1))
        activeDir[np.where(savemax > 0)[0]] = 1
        dMode_fine = vibration.gradvibration_an(
            *self._fineParametrization.get_inputs4gradVibration(
                VECT, w0, activeElem=activeElem, activeDir=activeDir))

        dMode = np.zeros((nbVar,nb_frq),np.float64)
        for i in np.arange(nb_frq):
            dMode[:,i] = self._propagateSensitivityAnalysis(vectX,dMode_fine[i])

        return dMode




    # ---------------------------------------------------------------------
    # Old stuff

    def compute_gradCompliance_cplgOnly_AN(self, vectX, listpatch=None):

        print("Calcul Gradient analytique compliance "
              "(terme couplage uniquement)...")
        if listpatch is None:
            listpatch = np.zeros(self._coarseParametrization._nb_patch,
                                 dtype=np.intp)
            listpatch[self._movable_patch] = 1
        U0 = self.compute_analysisSol(vectX)
        nbVar = self._nb_DesignVar
        dComp = np.zeros(nbVar)
        SOL, u = rsol.reconstruction(
            **self._fineParametrization.get_inputs4solution(U0))

        # dComp_fine =-gcplgAN.gradcplg_an(
        dComp_fine = -gradcplg_an(
            *self._fineParametrization.get_inputs4gradCoupling(
                SOL.transpose(), activepatch=listpatch))
        dShape_coarse = np.zeros(
            self._coarseParametrization._COORDS.shape[:: -1])
        for numVar in np.arange(0, nbVar):
            dShape_coarse[self._cp_by_var[numVar], :] = \
                self._dir_by_var[numVar].transpose()
            dShape_fine = self._mapCoarse2Fine(dShape_coarse)
            dComp[numVar] = np.sum(np.multiply(
                dShape_fine, dComp_fine.transpose()))
            dShape_coarse[self._cp_by_var[numVar], :] = 0.

        return dComp


    def _compute_gradCoupling_semiAN(self, vectX, eps=1.e-6, centerFD=False):
        U0 = self.compute_analysisSol(vectX) # state solution
        stateSOL   = U0.copy()
        adjointSOL = 0.5*U0.copy()

        ndof = self._fineParametrization._nb_dof_free
        idof = self._fineParametrization._ind_dof_free[:ndof] - 1

        def workWcplg(x):
            self._updateCoarseCoords(x)
            self._updateFineCoords()

            Cdata, Crow, Ccol = \
                cplg_matrix(*self._fineParametrization.get_inputs4cplgmatrix())
            Cside = sp.coo_matrix(
                (Cdata, (Crow, Ccol)),
                shape=(self._fineParametrization._nb_dof_tot,
                       self._fineParametrization._nb_dof_tot),
                dtype='float64').tocsc()
            Ctot = Cside + Cside.transpose()
            C2solve = Ctot[idof, :][:, idof]

            Wcplg = -np.dot(adjointSOL.T,C2solve.dot(stateSOL))

            self._updateCoarseCoords(self._current_vectX)
            self._updateFineCoords()

            return Wcplg

        gradW = self._sensitivityFD(workWcplg,vectX,eps=eps,centerFD=centerFD)
        return gradW

    @property
    def fine_parametrization(self):
        """
        IGAparametrization object containing analysis model
        """
        return self._fineParametrization

    @property
    def coarse_parametrization(self):
        """
        IGA parametrization object containing optimization model
        """
        return self._coarseParametrization

    @property
    def save_sol_fine(self):
        """
        Saved solution computed on analysis model
        """
        return self._save_sol_fine
