# Copyright 2020 Thibaut Hirschler
# Copyright 2020 Arnaud Duval

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

# Python modules
# --------------
import os
import sys
import time

import numpy as np
from scipy.sparse import csc_matrix, linalg as sla
import scipy.sparse as sp
from copy import deepcopy


# IGA Modules
# -----------
from ...stiffmtrx_elemstorage import sys_linmat_lindef_static as build_stiffmatrix
from ...massmtrx import build_cmassmatrix
from ...coupling.cplgmatrix import cplg_matrix
from ... import reconstructionSOL as rsol
from ...postprocessing import postproc as pp

from ...optim import volume as vol
#from optim.gradcompliance import gradcomp_semian, gradcomp_an, gradcplg_an




class DDOPTmodelling:
    'A class defining shape optimization problems.'

    def __init__(self,initialSubDomain, nb_DesingVar, fct_updateCoarseCoords,
                 nb_degreeElevationByDirection = np.array([0,0,0]),
                 nb_refinementByDirection      = np.array([0,0,0])):

        self._coarseSubDomain = deepcopy(initialSubDomain)
        self._fineSubDomain   = deepcopy(initialSubDomain)

        self._transformationMatrices_coarse2fine \
            = self._fineSubDomain._modeleIGA.refine_and_getTransformationMatrices(
                nb_refinementByDirection, nb_degreeElevationByDirection)
        self._fineSubDomain.set_dofInfos()

        self._nb_DesingVar = int(nb_DesingVar)
        self._special_updateCoarseCoords = fct_updateCoarseCoords

        # save
        self._initialCOORDS     = np.copy(initialSubDomain._modeleIGA._COORDS)
        self._current_vectX     = np.random.random(nb_DesingVar)
        self._save_sol_fine     = np.zeros(self._fineSubDomain._modeleIGA._nb_dof_free)
        self._save_secondmembre = np.zeros(self._fineSubDomain._modeleIGA._nb_dof_free)
        self._cp_by_var         = self.get_indCPbyDesignVar()
        self._dir_by_var        = self.get_dirCPbyDesignVar()
        self._elem_by_var       = self.get_ElembyDesignVar()
        self._movable_patch     = self.get_movablePatchesInd()

        # only for optimisation history plots
        self._OPT_FILENAME   = 'Opt'
        self._optimIteration = 0
        self._maxiter  = 10
        self._saveComp = 1.0
        self._comp0    = 1.0

        return None



    def get_indCPbyDesignVar(self):
        ''' Get the indices of the Control Points linked with each Design Variable '''
        self._updateCoarseCoords(np.zeros(self._nb_DesingVar))
        coords0 = np.copy(self._coarseSubDomain._modeleIGA._COORDS.transpose())

        cp_by_var = np.zeros(self._nb_DesingVar, dtype=object)
        for num_var in np.arange(0, self._nb_DesingVar):
            vect_test = np.zeros(self._nb_DesingVar)
            vect_test[num_var] = 1.
            self._updateCoarseCoords(vect_test)

            coords_test = self._coarseSubDomain._modeleIGA._COORDS.transpose()
            test        = np.logical_not(np.isclose(coords_test, coords0))
            cp_by_var[num_var] = np.unique(np.where(test)[0])

        self._updateCoarseCoords(np.zeros(self._nb_DesingVar))
        return cp_by_var

    def get_dirCPbyDesignVar(self):
        ''' Get the moving direction of the Control Points linked with each Design Variable '''
        self._updateCoarseCoords(np.zeros(self._nb_DesingVar))
        coords0 = np.copy(self._coarseSubDomain._modeleIGA._COORDS)

        dir_by_var = np.zeros(self._nb_DesingVar, dtype=object)
        for num_var in np.arange(0, self._nb_DesingVar):
            vect_test = np.zeros(self._nb_DesingVar)
            vect_test[num_var] = 1.
            self._updateCoarseCoords(vect_test)

            coords_test = self._coarseSubDomain._modeleIGA._COORDS
            dir_by_var[num_var] = (coords_test-coords0)[:,self._cp_by_var[num_var]]

        self._updateCoarseCoords(np.zeros(self._nb_DesingVar))
        return dir_by_var

    def get_ElembyDesignVar(self):
        ''' Get the list of Elements linked with each Design Variable '''
        self._updateCoarseCoords(np.zeros(self._nb_DesingVar))
        self._updateFineCoords()
        coords0 = np.copy(self._fineSubDomain._modeleIGA._COORDS.transpose())

        offset = np.cumsum(np.concatenate(([0],self._fineSubDomain._modeleIGA._elementsByPatch)))
        elem_by_var = np.zeros(self._nb_DesingVar, dtype=object)
        for num_var in np.arange(0,self._nb_DesingVar):
            vect_test = np.zeros(self._nb_DesingVar)
            vect_test[num_var] = 1.
            self._updateCoarseCoords(vect_test)
            self._updateFineCoords()

            coords_test = self._fineSubDomain._modeleIGA._COORDS.transpose()
            test        = np.logical_not(np.isclose(coords_test, coords0))
            list_cps    = np.unique(np.where(test)[0]) + 1
            activeElem  = np.zeros(self._fineSubDomain._modeleIGA._nb_elem,dtype=np.intp)
            for patch in np.arange(0,self._fineSubDomain._modeleIGA._nb_patch):
                for cp in list_cps:
                    activeElem[np.where(cp == self._fineSubDomain._modeleIGA._IEN[patch])[0]
                               + offset[patch]] = 1
            elem_by_var[num_var] = activeElem

        self._updateCoarseCoords(np.zeros(self._nb_DesingVar))
        self._updateFineCoords()
        return elem_by_var

    def get_movablePatchesInd(self):
        if not hasattr(self,'_cp_by_var'):
            self._cp_by_var = self.get_indCPbyDesignVar()
        cpall = np.concatenate(self._cp_by_var)
        ipatch= []
        for p in np.arange(self._coarseSubDomain._modeleIGA._nb_patch):
            if np.any(np.isin(self._coarseSubDomain._modeleIGA._indCPbyPatch[p]-1,cpall)):
                ipatch.append(p)
        return np.array(ipatch,dtype=np.intp)

    def _updateCoarseCoords(self, vectX):
        self._special_updateCoarseCoords(
            self._initialCOORDS, self._coarseSubDomain._modeleIGA, vectX)
        return None



    def _updateFineCoords(self):
        COORDS_coarse       = self._coarseSubDomain._modeleIGA._COORDS
        vectWeight_coarse   = self._coarseSubDomain._modeleIGA._vectWeight
        vectWeight_fine     = self._fineSubDomain._modeleIGA._vectWeight
        indCPbyPatch_coarse = self._coarseSubDomain._modeleIGA._indCPbyPatch
        indCPbyPatch_fine   = self._fineSubDomain._modeleIGA._indCPbyPatch

        for num_patch in range(0,self._coarseSubDomain._modeleIGA._nb_patch):

            nb_cp_thisPatch = np.size(indCPbyPatch_coarse[num_patch],0)
            Pwi_thisPatch = np.zeros((nb_cp_thisPatch, 4))
            Pwi_thisPatch[:,0] = COORDS_coarse[0, indCPbyPatch_coarse[num_patch]-1] \
                                 * vectWeight_coarse[indCPbyPatch_coarse[num_patch]-1]
            Pwi_thisPatch[:,1] = COORDS_coarse[1, indCPbyPatch_coarse[num_patch]-1] \
                                 * vectWeight_coarse[indCPbyPatch_coarse[num_patch]-1]
            Pwi_thisPatch[:,2] = COORDS_coarse[2, indCPbyPatch_coarse[num_patch]-1] \
                                 * vectWeight_coarse[indCPbyPatch_coarse[num_patch]-1]
            Pwi_thisPatch[:,3] = vectWeight_coarse[indCPbyPatch_coarse[num_patch]-1]

            Qwi_xi_thisPatch     = self._transformationMatrices_coarse2fine[num_patch][0] \
                                   * Pwi_thisPatch
            Qwi_xi_eta_thisPatch = self._transformationMatrices_coarse2fine[num_patch][1] \
                                   * Qwi_xi_thisPatch
            Qwi_thisPatch        = self._transformationMatrices_coarse2fine[num_patch][2] \
                                   * Qwi_xi_eta_thisPatch

            nb_cp_fine_thisPatch = np.size(Qwi_thisPatch,0)
            Qi_thisPatch = np.zeros((nb_cp_fine_thisPatch,3))
            Qi_thisPatch[:,0] = Qwi_thisPatch[:,0] / Qwi_thisPatch[:,3]
            Qi_thisPatch[:,1] = Qwi_thisPatch[:,1] / Qwi_thisPatch[:,3]
            Qi_thisPatch[:,2] = Qwi_thisPatch[:,2] / Qwi_thisPatch[:,3]

            self._fineSubDomain._modeleIGA._COORDS[:,indCPbyPatch_fine[num_patch]-1] \
                = Qi_thisPatch.transpose()
        return None

    def _mapCoarse2Fine(self, coarseQ):
        vectWeight_coarse   = self._coarseSubDomain._modeleIGA._vectWeight
        vectWeight_fine     = self._fineSubDomain._modeleIGA._vectWeight
        fineQ               = np.zeros(self._fineSubDomain._modeleIGA._COORDS.shape[::-1])
        indCPbyPatch_coarse = self._coarseSubDomain._modeleIGA._indCPbyPatch
        indCPbyPatch_fine   = self._fineSubDomain._modeleIGA._indCPbyPatch

        for num_patch in range(0,self._coarseSubDomain._modeleIGA._nb_patch):

            nb_cp_thisPatch = np.size(indCPbyPatch_coarse[num_patch],0)
            Pwi_thisPatch = np.zeros((nb_cp_thisPatch, 4))
            Pwi_thisPatch[:,0] = coarseQ[indCPbyPatch_coarse[num_patch]-1, 0] \
                                 * vectWeight_coarse[indCPbyPatch_coarse[num_patch]-1]
            Pwi_thisPatch[:,1] = coarseQ[indCPbyPatch_coarse[num_patch]-1, 1] \
                                 * vectWeight_coarse[indCPbyPatch_coarse[num_patch]-1]
            Pwi_thisPatch[:,2] = coarseQ[indCPbyPatch_coarse[num_patch]-1, 2] \
                                 * vectWeight_coarse[indCPbyPatch_coarse[num_patch]-1]
            Pwi_thisPatch[:,3] = vectWeight_coarse[indCPbyPatch_coarse[num_patch]-1]

            Qwi_xi_thisPatch     = self._transformationMatrices_coarse2fine[num_patch][0] \
                                   * Pwi_thisPatch
            Qwi_xi_eta_thisPatch = self._transformationMatrices_coarse2fine[num_patch][1] \
                                   * Qwi_xi_thisPatch
            Qwi_thisPatch        = self._transformationMatrices_coarse2fine[num_patch][2] \
                                   * Qwi_xi_eta_thisPatch

            nb_cp_fine_thisPatch = np.size(Qwi_thisPatch,0)
            Qi_thisPatch = np.zeros((nb_cp_fine_thisPatch,3))
            Qi_thisPatch[:,0] = Qwi_thisPatch[:,0] / Qwi_thisPatch[:,3]
            Qi_thisPatch[:,1] = Qwi_thisPatch[:,1] / Qwi_thisPatch[:,3]
            Qi_thisPatch[:,2] = Qwi_thisPatch[:,2] / Qwi_thisPatch[:,3]

            fineQ[indCPbyPatch_fine[num_patch]-1,:] = Qi_thisPatch[:,:]
        return fineQ

    def _transposeFine2Coarse(self, fineQ):
        vectWeight_coarse   = self._coarseSubDomain._modeleIGA._vectWeight
        vectWeight_fine     = self._fineSubDomain._modeleIGA._vectWeight
        coarseQ             = np.zeros(self._coarseSubDomain._modeleIGA._COORDS.shape[::-1])
        indCPbyPatch_coarse = self._coarseSubDomain._modeleIGA._indCPbyPatch
        indCPbyPatch_fine   = self._fineSubDomain._modeleIGA._indCPbyPatch

        for num_patch in range(0,self._fineSubDomain._modeleIGA._nb_patch):

            nb_cp_thisPatch = np.size(indCPbyPatch_fine[num_patch],0)
            Qwi_thisPatch = np.zeros((nb_cp_thisPatch, 4))
            Qwi_thisPatch[:,0] = fineQ[indCPbyPatch_fine[num_patch]-1, 0] \
                                 / vectWeight_fine[indCPbyPatch_fine[num_patch]-1]
            Qwi_thisPatch[:,1] = fineQ[indCPbyPatch_fine[num_patch]-1, 1] \
                                 / vectWeight_fine[indCPbyPatch_fine[num_patch]-1]
            Qwi_thisPatch[:,2] = fineQ[indCPbyPatch_fine[num_patch]-1, 2] \
                                 / vectWeight_fine[indCPbyPatch_fine[num_patch]-1]


            Pwi_zeta_thisPatch     = self._transformationMatrices_coarse2fine[num_patch][2].T \
                                     * Qwi_thisPatch
            Pwi_eta_zeta_thisPatch = self._transformationMatrices_coarse2fine[num_patch][1].T \
                                     * Pwi_zeta_thisPatch
            Pwi_thisPatch          = self._transformationMatrices_coarse2fine[num_patch][0].T \
                                     * Pwi_eta_zeta_thisPatch

            nb_cp_coarse_thisPatch = np.size(Pwi_thisPatch,0)
            vectWeight_coarse_thisPatch = vectWeight_coarse[indCPbyPatch_coarse[num_patch]-1]
            Pi_thisPatch = np.zeros((nb_cp_coarse_thisPatch,3))
            Pi_thisPatch[:,0] = Pwi_thisPatch[:,0] * vectWeight_coarse_thisPatch
            Pi_thisPatch[:,1] = Pwi_thisPatch[:,1] * vectWeight_coarse_thisPatch
            Pi_thisPatch[:,2] = Pwi_thisPatch[:,2] * vectWeight_coarse_thisPatch

            coarseQ[indCPbyPatch_coarse[num_patch]-1,:] = Pi_thisPatch[:,:]
        return coarseQ




    def _mapFine2Coarse(self,fineQ):
        coarseQ             = np.zeros(self._coarseSubDomain._modeleIGA._COORDS.shape[::-1])
        vectWeight_coarse   = self._coarseSubDomain._modeleIGA._vectWeight
        vectWeight_fine     = self._fineSubDomain._modeleIGA._vectWeight
        indCPbyPatch_coarse = self._coarseSubDomain._modeleIGA._indCPbyPatch
        indCPbyPatch_fine   = self._fineSubDomain._modeleIGA._indCPbyPatch

        for num_patch in range(0,self._coarseSubDomain._modeleIGA._nb_patch):

            nb_cp_thisPatch = np.size(indCPbyPatch_fine[num_patch],0)
            Pwi_thisPatch = np.zeros((nb_cp_thisPatch, 4))
            Pwi_thisPatch[:,0] = fineQ[indCPbyPatch_fine[num_patch]-1, 0] \
                                 * vectWeight_fine[indCPbyPatch_fine[num_patch]-1]
            Pwi_thisPatch[:,1] = fineQ[indCPbyPatch_fine[num_patch]-1, 1] \
                                 * vectWeight_fine[indCPbyPatch_fine[num_patch]-1]
            Pwi_thisPatch[:,2] = fineQ[indCPbyPatch_fine[num_patch]-1, 2] \
                                 * vectWeight_fine[indCPbyPatch_fine[num_patch]-1]
            Pwi_thisPatch[:,3] = vectWeight_fine[indCPbyPatch_fine[num_patch]-1]


            Qwi_zeta_thisPatch = self._transformationMatrices_fine2coarse[num_patch][2] \
                                 * Pwi_thisPatch
            Qwi_eta_thisPatch  = self._transformationMatrices_fine2coarse[num_patch][1] \
                                 * Qwi_zeta_thisPatch
            Qwi_thisPatch      = self._transformationMatrices_fine2coarse[num_patch][0] \
                                 * Qwi_eta_thisPatch

            nb_cp_coarse_thisPatch = np.size(Qwi_thisPatch,0)
            Qi_thisPatch = np.zeros((nb_cp_coarse_thisPatch,3))
            Qi_thisPatch[:,0] = Qwi_thisPatch[:,0] / Qwi_thisPatch[:,3]
            Qi_thisPatch[:,1] = Qwi_thisPatch[:,1] / Qwi_thisPatch[:,3]
            Qi_thisPatch[:,2] = Qwi_thisPatch[:,2] / Qwi_thisPatch[:,3]

            coarseQ[indCPbyPatch_coarse[num_patch]-1,:] = Qi_thisPatch[:,:]
        return coarseQ



    def _build_transformationMatrices_fine2coarse(self):

        saveR = []
        for num_patch in range(0, self._coarseSubDomain._modeleIGA._nb_patch):

            S_xi_thisPatch   = self._transformationMatrices_coarse2fine[num_patch][0]
            R_xi_thisPatch   = sla.inv( (S_xi_thisPatch.transpose() * S_xi_thisPatch).tocsc() ) \
                             * S_xi_thisPatch.transpose()

            S_eta_thisPatch  = self._transformationMatrices_coarse2fine[num_patch][1]
            R_eta_thisPatch  = sla.inv( (S_eta_thisPatch.transpose() * S_eta_thisPatch).tocsc() ) \
                              * S_eta_thisPatch.transpose()

            S_zeta_thisPatch = self._transformationMatrices_coarse2fine[num_patch][2]
            R_zeta_thisPatch = sla.inv( (S_zeta_thisPatch.transpose() * S_zeta_thisPatch).tocsc())\
                              * S_zeta_thisPatch.transpose()

            saveR.append([R_xi_thisPatch, R_eta_thisPatch, R_zeta_thisPatch])

        return saveR



    def set_DDanalysis(self, vectX, activeelem=None, tol=1e-6, pseudo=True):
        if np.all(vectX == self._current_vectX):
            localG = self._fineSubDomain._localG
            localt = self._fineSubDomain.compute_condensedRHSvect()
            locale = self._fineSubDomain.compute_rigidbodyRHSvect()
            localC = self._fineSubDomain._C2solve
            localinvdK = self._fineSubDomain._invDiagK.diagonal()
        else:
            print(" Set sub-domain operators...")
            self._updateCoarseCoords(vectX)
            self._updateFineCoords()

            self._fineSubDomain.set_stiffnessMATRIX()
            self._fineSubDomain.set_couplingMATRIX()
            if pseudo is True:
                self._fineSubDomain.set_factorizationMATRIX(tol=tol)
            self._fineSubDomain.set_admissibleconstMATRIX()
            self._fineSubDomain.set_invdiaKMATRIX()
            self._fineSubDomain.set_factorizationInternalMATRIX()

            localG = self._fineSubDomain._localG
            localt = self._fineSubDomain.compute_condensedRHSvect()
            locale = self._fineSubDomain.compute_rigidbodyRHSvect()
            localC = self._fineSubDomain._C2solve
            localinvdK = self._fineSubDomain._invDiagK.diagonal()

            self._current_vectX[:] = vectX[:]

        return localG,localt,locale,localC,localinvdK


    def set_analysisSol(self, lmbdas,alphas):
        self._save_lmbda_fine   = lmbdas.copy()
        self._save_sol_fine     = self._fineSubDomain.evaluate_displacement(lmbdas,alphas)
        self._save_secondmembre = self._fineSubDomain._f2solve
        self._save_raideur      = self._fineSubDomain._K2solve
        self._save_coupling     = self._fineSubDomain._C2solve

        return None


    def compute_compliance_discrete(self):
        U0 = self._save_sol_fine
        F0 = self._save_secondmembre
        comp = 0.5*np.dot(U0,F0)
        self._saveComp = comp
        return comp


    def compute_volume(self,vectX,listpatch=None):
        print("Calcul Volume...")
        self._updateCoarseCoords(vectX)
        self._updateFineCoords()
        V = vol.computevolume( *self._fineSubDomain._modeleIGA.get_inputs4area(
            activepatch=listpatch))
        return V


    def compute_gradVolume_AN(self,vectX,listpatch=None):
        print("Calcul Gradient Analytique volume...")

        self._updateCoarseCoords(vectX)
        self._updateFineCoords()

        nbVar = self._nb_DesingVar
        dVol  = np.zeros(nbVar)

        dVol_fine = vol.computegradvolume_wembded(
                *self._fineSubDomain._modeleIGA.get_inputs4area(activepatch=listpatch) )

        dShape_coarse = np.zeros( self._coarseSubDomain._modeleIGA._COORDS.shape[::-1] )
        for numVar in np.arange(0,nbVar):
            dShape_coarse[self._cp_by_var[numVar],:] = self._dir_by_var[numVar].transpose()
            dShape_fine = self._mapCoarse2Fine(dShape_coarse)

            dVol[numVar] = np.sum(np.multiply( dShape_fine,dVol_fine.transpose() ))

            dShape_coarse[self._cp_by_var[numVar],:] = 0.

        return dVol


    def compute_gradVolume_DF(self,vectX,eps=1.e-4,listpatch=None):
        print("Calcul Gradient differences finies volume...")
        self._updateCoarseCoords(vectX)
        self._updateFineCoords()
        V0 = vol.computevolume(
            *self._fineSubDomain._modeleIGA.get_inputs4area(activepatch=listpatch) )

        save  = vectX
        nbVar = self._nb_DesingVar
        gradV_df      = np.zeros(nbVar)
        vectX_thisVar = np.zeros(nbVar)
        for numVar in np.arange(0,nbVar):

            sys.stdout.write(' numVar:%3i/%i\r' % (numVar, self._nb_DesingVar) )
            sys.stdout.flush()

            vectX_thisVar[:]       = save[:]
            vectX_thisVar[numVar] += eps

            self._updateCoarseCoords(vectX_thisVar)
            self._updateFineCoords()
            V = vol.computevolume(
                *self._fineSubDomain._modeleIGA.get_inputs4area(activepatch=listpatch) )

            gradV_df[numVar] = (V-V0)/eps

        self._updateCoarseCoords(save)
        self._updateFineCoords()

        print(' Done.'+' '*30)

        return gradV_df


    def compute_gradCompliance_semiAN(self,vectX,eps = 1.e-6):

        print("Calcul Gradient semi-analytique compliance...")

        U0 = self._save_sol_fine
        l0 = self._save_lmbda_fine
        F0 = self._save_secondmembre
        K0 = self._save_raideur
        C0 = self._save_coupling

        #epsilon = 1.e-5
        save  = vectX
        nbVar = self._nb_DesingVar
        dComp = np.zeros(nbVar)
        vectX_thisVar = np.zeros(nbVar)

        for numVar in np.arange(0,nbVar):

            sys.stdout.write(' numVar:%3i/%i\r' % (numVar, self._nb_DesingVar) )
            sys.stdout.flush()

            vectX_thisVar[:]       = save[:]
            vectX_thisVar[numVar] += eps

            self._updateCoarseCoords(vectX_thisVar)
            self._updateFineCoords()

            self._fineSubDomain.set_stiffnessMATRIX()
            self._fineSubDomain.set_couplingMATRIX()

            K_thisVar = self._fineSubDomain._K2solve
            F_thisVar = self._fineSubDomain._f2solve
            C_thisVar = self._fineSubDomain._C2solve

            Re_thisVar = (F_thisVar-F0)/eps - 0.5*(K_thisVar-K0)*U0/eps
            Re_thisVar-= (C_thisVar - C0).T*l0/eps
            dComp[numVar] = np.dot(U0, Re_thisVar)

            del K_thisVar
            del F_thisVar
            del C_thisVar

        self._updateCoarseCoords(save)
        self._updateFineCoords()

        self._fineSubDomain._f2solve = F0.copy()
        self._fineSubDomain._K2solve = K0.copy()
        self._fineSubDomain._C2solve = C0.copy()

        print(' Done.          ')

        return dComp



    def compute_gradCompliance_AN(self,vectX):

        print("Calcul Gradient analytique compliance...")

        U0 = self.compute_analysisSol(vectX)

        nbVar = self._nb_DesingVar
        dComp = np.zeros(nbVar)

        SOL,u = rsol.reconstruction(
            *self._fineParametrization.get_inputs4solution(U0))
        if self._fineParametrization._mcrd == 2:
            SOL3D = np.c_[SOL,np.zeros(SOL.shape[0])]
        else:
            SOL3D = SOL
        activeElem = np.zeros(self._fineParametrization._nb_elem, dtype=np.intp)
        activeElem[np.where(np.sum(self._elem_by_var)>0)[0]] = 1

        activeDir = np.zeros(3, dtype=np.intp)
        savemax   = np.zeros(3)
        for var in self._dir_by_var:
            savemax = np.maximum(savemax, np.max(np.abs(var), axis=1))
        activeDir[np.where(savemax>0)[0]] = 1

        #dComp_fine = gcompAN.gradcomp_an(
        dComp_fine = gradcomp_an(
            *self._fineParametrization.get_inputs4gradCompliance(
                SOL3D.transpose(),
                activeElem= activeElem,
                activeDir = activeDir) )

        dComp_coarse = np.transpose( self._transposeFine2Coarse(dComp_fine.T) )

        dShape_coarse = np.zeros( self._coarseParametrization._COORDS.shape[::-1] )
        for numVar in np.arange(0,nbVar):
            dShape_coarse[self._cp_by_var[numVar],:] = self._dir_by_var[numVar].transpose()
            dShape_fine = self._mapCoarse2Fine(dShape_coarse)
            dComp[numVar] = np.sum(np.multiply( dShape_fine,dComp_fine.transpose() ))
            dShape_coarse[self._cp_by_var[numVar],:] = 0.

        #dShape_coarse = np.zeros_like( self._coarseParametrization._COORDS )
        #for numVar in np.arange(0,nbVar):
        #    dShape_coarse[:,self._cp_by_var[numVar]] = self._dir_by_var[numVar]
        #    dComp[numVar] = np.tensordot(dShape_coarse, dComp_coarse, axes=2)
        #    dShape_coarse[:,self._cp_by_var[numVar]] = 0.

        return dComp


    def compute_gradCompliance_cplgOnly_AN(self,vectX,listpatch=None):

        print("Calcul Gradient analytique compliance (terme couplage uniquement)...")

        if listpatch is None:
            listpatch = np.zeros(self._coarseParametrization._nb_patch,dtype=np.intp)
            listpatch[self._movable_patch] = 1

        U0 = self.compute_analysisSol(vectX)

        nbVar = self._nb_DesingVar
        dComp = np.zeros(nbVar)

        SOL,u = rsol.reconstruction(
            *self._fineParametrization.get_inputs4solution(U0))

        #dComp_fine =-gcplgAN.gradcplg_an(
        dComp_fine =-gradcplg_an(
            *self._fineParametrization.get_inputs4gradCoupling(
                SOL.transpose(),activepatch=listpatch) )

        dShape_coarse = np.zeros( self._coarseParametrization._COORDS.shape[::-1] )
        for numVar in np.arange(0,nbVar):
            dShape_coarse[self._cp_by_var[numVar],:] = self._dir_by_var[numVar].transpose()
            dShape_fine = self._mapCoarse2Fine(dShape_coarse)
            dComp[numVar] = np.sum(np.multiply( dShape_fine,dComp_fine.transpose() ))
            dShape_coarse[self._cp_by_var[numVar],:] = 0.

        return dComp




    def compute_gradCompliance_cplgOnly_semiAN(self,vectX,eps = 1.e-6):

        print("Calcul Gradient semi-analytique compliance (terme couplage uniquement)...")

        U0 = self.compute_analysisSol(vectX)

        ndof = self._fineParametrization._nb_dof_free
        idof = self._fineParametrization._ind_dof_free[:ndof]-1

        Cdata,Crow,Ccol = cplg_matrix( *self._fineParametrization.get_inputs4cplgmatrix() )
        Cside = sp.coo_matrix((Cdata,(Crow,Ccol)),
                              shape=(self._fineParametrization._nb_dof_tot,
                                     self._fineParametrization._nb_dof_tot),
                              dtype='float64').tocsc()
        Ctot  = Cside + Cside.transpose()
        C0 = Ctot[idof,:][:,idof]
        del Cdata,Crow,Ccol,Cside,Ctot

        #epsilon = 1.e-5
        save  = vectX
        nbVar = self._nb_DesingVar
        dComp = np.zeros(nbVar)
        vectX_thisVar = np.zeros(nbVar)

        for numVar in np.arange(0,nbVar):

            sys.stdout.write(' numVar:%3i/%i\r' % (numVar, self._nb_DesingVar) )
            sys.stdout.flush()

            vectX_thisVar[:]       = save[:]
            vectX_thisVar[numVar] += eps

            self._updateCoarseCoords(vectX_thisVar)
            self._updateFineCoords()

            Cdata,Crow,Ccol = cplg_matrix( *self._fineParametrization.get_inputs4cplgmatrix() )
            Cside = sp.coo_matrix((Cdata,(Crow,Ccol)),
                                  shape=(self._fineParametrization._nb_dof_tot,
                                         self._fineParametrization._nb_dof_tot),
                                  dtype='float64').tocsc()
            Ctot  = Cside + Cside.transpose()
            C_thisVar = Ctot[idof,:][:,idof]
            del Cdata,Crow,Ccol,Cside,Ctot

            Re_thisVar = - 0.5*(C_thisVar-C0)*U0/eps
            dComp[numVar] = np.dot(U0, Re_thisVar)

            del C_thisVar

        self._updateCoarseCoords(save)
        self._updateFineCoords()

        print(' Done.          ')

        return dComp
