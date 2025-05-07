# Copyright 2018-2021 Thibaut Hirschler
# Copyright 2020-2023 Arnaud Duval
# Copyright 2021 Marie Guerder

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

"""
CONTAINS:
IGAparametrization class to generate the datasetting for an IsoGeometric
Analysis.
"""

import numpy as np
from scipy import sparse as sp

from ...DOF import getinddof
from ...postprocessing.postproc import generate_vtk, generate_vtk_wsol
from ..geometricmodel import NBfile
from ..mechanicalmodel import MechanicalModel as INPfile
from .IGA_refinementFcts import iga_refinement
from .IGA_manipulation import bezier_decomposition_patch



class IGAparametrization:
    """A class containing the necessary parameters for IGA."""

    def __init__(self, mechanicalSettings=None, geometricSettings=None,
                 filename=None):
        """Initialise IGA parametrization class object.

        :param str filename: short name (without exetension) of input files to read,
            defaults to ``None``
        :param list mechanicalSettings: list containing mechanical settings, defaults to ``None``
        :param list geometricSettings: list containing geometric settings, defaults to ``None``

        Raises
        ------
        Exception
            DESCRIPTION.

        """
        if filename is not None:
            # *read files
            # - NB file
            nbfile = NBfile(filename.strip()+'.NB')
            geometricSettings = nbfile.getGeoInfos()
            # - INP file
            inpfile = INPfile()
            inpfile.read(filename=filename+'.inp', format='ABAQUS')
            mechanicalSettings = []
            mechanicalSettings.append(inpfile.get_parameters())             # [0]
            mechanicalSettings.append(inpfile.get_bc())                     # [1]
            mechanicalSettings.append(inpfile.get_load())                   # [2]
            mechanicalSettings.append(inpfile.get_nodes())                  # [3]
            mechanicalSettings.append(inpfile.get_ien())                    # [4]
            mechanicalSettings.append(inpfile.get_material_properties())    # [5]
            mechanicalSettings.append(inpfile.get_properties())             # [6]
            mechanicalSettings.append(inpfile.get_tables())                 # [7]
            mechanicalSettings.append(inpfile._get_shapeparametrization())
        elif not mechanicalSettings:
            return
        elif not geometricSettings:
            return

        # Mechanical properties:
        # ---------------------
        # Parameters
        self._ELT_TYPE = mechanicalSettings[0][0]
        self._NBPINT = mechanicalSettings[0][1]
        self._TENSOR = mechanicalSettings[0][2]
        self._mcrd = mechanicalSettings[0][3]

        # Boundary conditions
        self._bc_target = mechanicalSettings[1][0]
        self._bc_values = mechanicalSettings[1][1]
        self._bc_target_nbelem = mechanicalSettings[1][2]
        self._nb_bc = mechanicalSettings[1][3]

        # Loads
        self._indDLoad = mechanicalSettings[2][0]
        self._JDLType = mechanicalSettings[2][1]
        self._ADLMAG = mechanicalSettings[2][2]
        self._load_target_nbelem = mechanicalSettings[2][3]
        self._nb_load = mechanicalSettings[2][4]

        # Case where additional loading information are defined
        # (centrifugal force axis, distribution,...)
        if len(mechanicalSettings[2]) > 5:
            self._additionalLoadInfos = mechanicalSettings[2][5]
        else:
            self._additionalLoadInfos = []
        self._nb_additionalLoadInfos = np.zeros(self._nb_load, dtype=int)

        for i_load, load_type in enumerate(self._JDLType):
            if load_type == 101:
                # centrifugal body force, additional infos = rotation axis
                self._nb_additionalLoadInfos[i_load] = 6
            elif load_type % 10 == 4:
                # distributed pressure, additional infos = distribution index
                self._nb_additionalLoadInfos[i_load] = 1

        # self._nb_cload= mechanicalSettings[2][5]

        # Nodes info
        self._COORDS = mechanicalSettings[3][0]
        self._nb_cp = mechanicalSettings[3][1]

        # Elements info
        self._IEN = mechanicalSettings[4]  # [0]
        # self._nb_elem = mechanicalSettings[4][1]

        # Material
        # TODO Handling of material properties could be improved
        self._MATERIAL_PROPERTIES = mechanicalSettings[5][0]
        self._N_MATERIAL_PROPERTIES = mechanicalSettings[5][1]

        # Properties
        self._PROPS = mechanicalSettings[6][0]
        self._JPROPS = mechanicalSettings[6][1]

        # Distributions
        # Assuming nodal distributions only
        self._nodal_distributions = {}
        # Keep a copy of the nodal distribution defined on the initial geometry
        self._nodal_distributions_init = {}

        if len(mechanicalSettings) > 7:
            for distrib in mechanicalSettings[7]:
                if distrib.opts['location'].strip() == 'NODE':
                    a = np.array(distrib)
                    n = np.size(a, axis=1) - 1
                    d = np.zeros((self._nb_cp, n), order='F')
                    d[:] = distrib.get_defaultvalues()
                    idx = a[:, 0].astype(int)
                    d[idx - 1] = a[:, 1:]
                    headers = distrib.opts['dtype']
                    self._nodal_distributions[distrib.get_name()] = \
                        [d.copy(), headers]
                    self._nodal_distributions_init[distrib.get_name()] = \
                        [d.copy(), headers]
                else:
                    raise Exception('Only nodal distributions are handled')

        # Shape parametrization
        if len(mechanicalSettings) > 8:
            self._design_parameters = mechanicalSettings[8][0]
            self._shapeparametrization_def = mechanicalSettings[8][1]

        # Geometric properties:
        # --------------------
        self._dim = geometricSettings[0]
        self._Nkv = geometricSettings[1]
        self._Ukv = geometricSettings[2]
        self._Jpqr = geometricSettings[3]
        self._Nijk = geometricSettings[4]
        self._weight = geometricSettings[5]
        self._elementsByPatch = geometricSettings[6]
        self._nb_patch = geometricSettings[7]
        self._nnode = geometricSettings[8]
        self._nb_elem = geometricSettings[9]

        # Other attributes:
        # ---------------
        self._flatten_data()
        self._indCPbyPatch = self._autoset_indCPbyPatch()
        self._compute_vectWeight()
        self._update_dof_info()
        self._initRefinementMatHistory()

    def get_mechanicalSettings(self):
        """Get the object's mechanical settings."""
        mechSet = []
        mechSet.append([self._ELT_TYPE, self._NBPINT, self._TENSOR,
                        self._mcrd])
        mechSet.append([self._bc_target, self._bc_values,
                        self._bc_target_nbelem, self._nb_bc])
        mechSet.append([self._indDLoad, self._JDLType, self._ADLMAG,
                        self._load_target_nbelem, self._nb_load])
        mechSet.append([self._COORDS, self._nb_cp])
        mechSet.append(self._IEN)
        mechSet.append([self._MATERIAL_PROPERTIES])
        mechSet.append([self._PROPS, self._JPROPS])
        return mechSet

    def get_mechanicalSettings_somePatch(self, listpatch, with_bc=False,
                                         updatePROPS=False):
        """Get the mechanical settings for a given patch.

        Parameters
        ----------
        listpatch : list of ints
            List of indices corresponding to the patche(s) of interest.
        with_bc : boolean, optional
            Adds boundary conditions to the returned list if set to True. The
            default is False.
        updatePROPS : boolean, optional
            Updates the object's properties if set to True. The default is
            False.

        Returns
        -------
        mechSet : list
            List of ordered mechanical settings.
        """
        mechSet = []
        mechSet.append([self._ELT_TYPE[listpatch], self._NBPINT[listpatch],
                        self._TENSOR[listpatch], self._mcrd])
        indCP = np.unique(np.concatenate(
            np.array([self._indCPbyPatch[i] for i in listpatch], dtype=int))) - 1
        tabCP = np.zeros(self._nb_cp, dtype=np.intp)
        tabCP[indCP] = np.arange(1, indCP.size+1)
        mechSet.append([self._COORDS[:, indCP], indCP.size])
        IEN = [tabCP[self._IEN[i] - 1] for i in listpatch]
        mechSet.append(IEN)
        mechSet.append((self._MATERIAL_PROPERTIES[:, listpatch],
                        self._N_MATERIAL_PROPERTIES[listpatch]))
        PROPS = [self._PROPS[i].copy() for i in listpatch]
        if updatePROPS:
            for i in range(0, len(listpatch)):
                patch = listpatch[i]
                PROPS[i][0] = i+1
                if (self._ELT_TYPE[patch] == 'U10' or
                        self._ELT_TYPE[patch] == 'U30'):
                    test = np.where(listpatch+1 == PROPS[i][1])
                    if np.size(test) > 0:
                        PROPS[i][1] = test[0]+1
                    else:
                        print('Warning: props of patch%i cannot be updated' +
                              ' (no corresponding volume).' % (i+1))
                if self._ELT_TYPE[patch] == 'U00':
                    test1 = np.where(listpatch+1 == PROPS[i][1])
                    test2 = np.where(listpatch+1 == PROPS[i][2])
                    if np.size(test1) > 0:
                        PROPS[i][1] = test1[0] + 1
                    else:
                        PROPS[i][1] = 0
                        print('Warning: props of curve%i cannot be updated' +
                              ' (no corresponding patch).' % (i+1))
                    if np.size(test2) > 0:
                        PROPS[i][2] = test2[0] + 1
                    else:
                        PROPS[i][2] = 0
                        print('Warning: props of curve%i cannot be updated' +
                              ' (no corresponding lgrge).' % (i+1))
        mechSet.append([PROPS, self._JPROPS[listpatch]])

        if with_bc:
            # dirichlet
            bc_target = []
            bc_target_nbelem = []
            indbcfind = np.zeros(self._nb_bc, dtype=np.intp)
            i = 0
            for bc in self._bc_target:
                thisbc = tabCP[np.array(bc, dtype=np.intp) - 1]
                nbelem = np.count_nonzero(thisbc)
                if nbelem > 0:
                    bc_target.append(list(thisbc[np.nonzero(thisbc)]))
                    bc_target_nbelem.append(nbelem)
                    indbcfind[i] = 1
                i += 1
            nb_bc = len(bc_target)
            bcData = [bc_target, self._bc_values[:, np.nonzero(indbcfind)[0]],
                      bc_target_nbelem, nb_bc]
            mechSet.insert(1, bcData)

            # neumann
            cumelem = np.insert(np.cumsum(self._elementsByPatch), 0, 0)
            indEL = np.concatenate(
                [np.arange(cumelem[i], cumelem[i+1]) for i in listpatch])
            tabEL = np.zeros(self._nb_elem, dtype=np.intp)
            tabEL[indEL] = np.arange(1, indEL.size+1)
            indDLoad = []
            load_target_nbelem = []
            indloadfind = np.zeros(self._nb_load, dtype=np.intp)
            i = 0
            for load in self._indDLoad:
                if self._JDLType[i] > 9:
                    thisload = tabEL[np.array(load, dtype=np.intp) - 1]
                else:
                    thisload = tabCP[np.array(load, dtype=np.intp) - 1]
                nbelem = np.count_nonzero(thisload)
                if nbelem > 0:
                    indDLoad.append(list(thisload[np.nonzero(thisload)]))
                    load_target_nbelem.append(nbelem)
                    indloadfind[i] = 1
                i += 1
            nb_load = len(indDLoad)
            index = np.nonzero(indloadfind)[0]
            loadData = [np.array(indDLoad), self._JDLType[index],
                        self._ADLMAG[index], np.array(load_target_nbelem),
                        nb_load, np.array([])]
            mechSet.insert(2, loadData)

        return mechSet

    def set_mechanicalSettings(self, newMechSet):
        """Set the mechanical settings."""
        [self._ELT_TYPE, self._NBPINT, self._TENSOR, self._mcrd] = \
            newMechSet[0]
        [self._bc_target, self._bc_values, self._bc_target_nbelem,
         self._nb_bc] = newMechSet[1]
        [self._indDLoad, self._JDLType, self._ADLMAG, self._load_target_nbelem,
         self._nb_load] = newMechSet[2]
        [self._COORDS, self._nb_cp] = newMechSet[3]
        self._IEN = newMechSet[4]
        [self._MATERIAL_PROPERTIES] = newMechSet[5]
        [self._PROPS, self._JPROPS] = newMechSet[6]

        return None

    def get_geometricSettings(self):
        """Get the geometrical settings."""
        geoSet = [self._dim, self._Nkv, self._Ukv, self._Jpqr, self._Nijk,
                  self._weight, self._elementsByPatch, self._nb_patch,
                  self._nnode, self._nb_elem]

        return geoSet

    def get_geometricSettings_somePatch(self, listpatch):
        """Get the geometrical settings of a given patch."""
        cumelem = np.insert(np.cumsum(self._elementsByPatch), 0, 0)
        listelem = np.concatenate(
            [np.arange(cumelem[i], cumelem[i+1]) for i in listpatch])
        Ukv = [self._Ukv[i] for i in listpatch]
        weight = [self._weight[e] for e in listelem]
        geoSet = [self._dim[listpatch], self._Nkv[:, listpatch], Ukv,
                  self._Jpqr[:, listpatch], self._Nijk[:, listelem], weight,
                  self._elementsByPatch[listpatch], listpatch.size,
                  self._nnode[listpatch], listelem.size]

        return geoSet

    def set_geomectricSettings(self, newGeoSet):
        """Set the geometrical settings."""
        [self._dim, self._Nkv, self._Ukv, self._Jpqr, self._Nijk, self._weight,
         self._elementsByPatch, self._nb_patch, self._nnode, self._nb_elem] = \
            newGeoSet

        return None

    def add_patch(self, patchMechSet, patchGeoSet, mode=0):
        """Add a patch to current IGAparametrization object.

        Parameters
        ----------
        patchMechSet : TYPE
            DESCRIPTION.
        patchGeoSet : TYPE
            DESCRIPTION.
        mode : TYPE, optional
            DESCRIPTION. The default is 0.

        Returns
        -------
        None.
        """
        if mode == 0:
            # --> independant new patch
            # Update mechanical props
            self._ELT_TYPE = np.append(self._ELT_TYPE, patchMechSet[0][0])
            self._NBPINT = np.append(self._NBPINT, patchMechSet[0][1])
            self._TENSOR = np.append(self._TENSOR, patchMechSet[0][2])
            self._COORDS = np.concatenate((self._COORDS, patchMechSet[1][0]),
                                          axis=1)
            self._nb_cp += patchMechSet[1][1]
            self._IEN.append(patchMechSet[2][0] +
                             self._nb_cp-patchMechSet[1][1])
            self._MATERIAL_PROPERTIES = \
                np.concatenate((self._MATERIAL_PROPERTIES,
                                np.vstack(patchMechSet[3][0])), axis=1)
            self._PROPS.append(patchMechSet[4][0])
            self._JPROPS = np.append(self._JPROPS, patchMechSet[4][1])
            # Update geometric props
            self._dim = np.append(self._dim, patchGeoSet[0])
            self._Nkv = np.concatenate((self._Nkv, np.vstack(patchGeoSet[1])),
                                       axis=1)
            self._Ukv.append(patchGeoSet[2])
            self._Jpqr = np.concatenate((self._Jpqr,
                                         np.vstack(patchGeoSet[3])), axis=1)
            self._Nijk = np.concatenate((self._Nijk,
                                         np.vstack(patchGeoSet[4])), axis=1)
            for w in patchGeoSet[5]:
                self._weight.append(w)
            self._elementsByPatch = np.append(self._elementsByPatch,
                                              patchGeoSet[7])
            self._nb_patch += 1
            self._nnode = np.append(self._nnode, patchGeoSet[6])
            self._nb_elem += patchGeoSet[7]

        self._flatten_data()
        self._indCPbyPatch = self._autoset_indCPbyPatch()
        self._compute_vectWeight()
        self._update_dof_info()

        return None

    def add_multiplePatch(self, patchMechSet, patchGeoSet):
        """

        Parameters
        ----------
        patchMechSet : TYPE
            DESCRIPTION.
        patchGeoSet : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        # --> independant patch
        # Update mechanical props
        self._ELT_TYPE = np.append(self._ELT_TYPE, patchMechSet[0][0])
        self._NBPINT = np.append(self._NBPINT, patchMechSet[0][1])
        self._TENSOR = np.append(self._TENSOR, patchMechSet[0][2])
        self._COORDS = np.concatenate((self._COORDS, patchMechSet[1][0]),
                                      axis=1)
        for e in patchMechSet[2]:
            self._IEN.append(e + self._nb_cp)
        self._nb_cp += patchMechSet[1][1]
        self._MATERIAL_PROPERTIES = np.concatenate(
            (self._MATERIAL_PROPERTIES, np.vstack(patchMechSet[3])),
            axis=1)
        for p in patchMechSet[4][0]:
            self._PROPS.append(p)
        self._JPROPS = np.append(self._JPROPS, patchMechSet[4][1])
        # Update geometric props
        self._dim = np.append(self._dim, patchGeoSet[0])
        self._Nkv = np.concatenate((self._Nkv, np.vstack(patchGeoSet[1])),
                                   axis=1)
        for u in patchGeoSet[2]:
            self._Ukv.append(u)
        self._Jpqr = np.concatenate((self._Jpqr, np.vstack(patchGeoSet[3])),
                                    axis=1)
        self._Nijk = np.concatenate((self._Nijk, np.vstack(patchGeoSet[4])),
                                    axis=1)
        for w in patchGeoSet[5]:
            self._weight.append(w)
        self._elementsByPatch = np.append(self._elementsByPatch,
                                          patchGeoSet[6])
        self._nb_patch += patchGeoSet[7]
        self._nnode = np.append(self._nnode, patchGeoSet[8])
        self._nb_elem += patchGeoSet[9]

        self._flatten_data()
        self._indCPbyPatch = self._autoset_indCPbyPatch()
        self._compute_vectWeight()
        self._update_dof_info()

        return None

    def _get_indCPbyPatch(self):
        """Get control points indices for each patch."""
        return self._indCPbyPatch

    def _set_indCPbyPatch(self, indCPbyPatch):
        """Set control points indices for each patch."""
        self._indCPbyPatch = indCPbyPatch

        return None

    def _get_vectWeight(self):
        """Get NURBS weight vector."""
        return self._vectWeight

    def _compute_vectWeight(self):
        """Build ordered NURBS weight vector."""
        self._vectWeight = np.ones(self._nb_cp)
        Jelem = 0
        for num_patch in range(0, self._nb_patch):
            for num_elem in range(0, self._elementsByPatch[num_patch]):
                self._vectWeight[self._IEN[num_patch][num_elem, :] - 1] = \
                    self._weight[Jelem][:]
                Jelem += 1

        return None

    def _update_dof_info(self):
        """Set degrees of freedom information."""
        self._flatten_data()
        self._COUPLG_flag = False
        self._nb_dof_tot = int(self._mcrd * self._nb_cp)
        self._nb_dof_bloq = int(np.sum(self._bc_target_nbelem))
        self._nb_dof_free = int(self._nb_dof_tot - self._nb_dof_bloq)
        if self._nb_dof_free == 0:
            self._ind_dof_bloq = np.arange(1, self._nb_dof_tot+1)
            self._ind_dof_free = np.array([])
        elif self._nb_dof_bloq == 0:
            self._ind_dof_free = np.arange(1, self._nb_dof_tot+1)
            self._ind_dof_bloq = np.array([])
        else:
            self._nb_dof_bloq, self._nb_dof_free, self._ind_dof_bloq, \
                self._ind_dof_free, self._COUPLG_flag, self._BNDSTRIP_flag = \
                getinddof(*self.get_inputs4indDOF())

        return None

    def _autoset_indCPbyPatch(self):
        """Automatic computation of control point indices by patch."""
        indCPbyPatch = np.zeros(self._nb_patch, dtype=object)
        JELEM = -1
        for num_patch in range(0, self._nb_patch):
            nnodeD = np.zeros(3, int)
            nnodeD[:] = self._Jpqr[:, num_patch] + 1
            nb_cpD = np.zeros(3, int)
            nb_cpD[:] = np.maximum(
                self._Nkv[:, num_patch] - (self._Jpqr[:, num_patch] + 1),
                np.ones(3))
            indCP_thisPatch = np.zeros(np.prod(nb_cpD), int)
            for num_elem in range(0, self._elementsByPatch[num_patch]):
                JELEM += 1
                indFirstCPbyDir = np.zeros(3, int)
                indFirstCPbyDir[:] = np.maximum(
                    self._Nijk[:, JELEM] - (self._Jpqr[:, num_patch] + 1),
                    np.zeros(3))
                list_xi = np.arange(0, nnodeD[0]) + indFirstCPbyDir[0]
                list_eta = np.array([])
                list_thisElem = np.array([])
                for eta in range(0, nnodeD[1]):
                    list_eta = np.append(
                        list_eta,
                        list_xi + indFirstCPbyDir[1]*nb_cpD[0] + eta*nb_cpD[0])
                for zeta in range(0, nnodeD[2]):
                    list_thisElem = np.append(
                        list_thisElem,
                        list_eta + indFirstCPbyDir[2]*np.prod(nb_cpD[0:2]) +
                        zeta*np.prod(nb_cpD[0: 2]))
                indCP_thisPatch[list_thisElem.astype(int)] = \
                    self._IEN[num_patch][num_elem, ::-1]

            indCPbyPatch[num_patch] = indCP_thisPatch

        return indCPbyPatch

    def _flatten_data(self):
        """Flatten all geometric or mechanical attributes arrays.

        This method is used for passing arguments to the FORTRAN part of the
        code, which is taking only flattened arrays.
        """
        self._weight_flat = np.concatenate(self._weight)

        self._Ukv_flat = np.array([], dtype=np.float64)
        for u in self._Ukv:
            self._Ukv_flat = np.hstack((self._Ukv_flat, np.concatenate(u)))

        self._IEN_flat = np.array([], dtype=np.intp)
        for ien in self._IEN:
            self._IEN_flat = np.hstack((self._IEN_flat, np.concatenate(ien)))

        self._ELT_TYPE_flat = ''
        for et in self._ELT_TYPE:
            self._ELT_TYPE_flat += et

        self._TENSOR_flat = ''
        for tn in self._TENSOR:
            self._TENSOR_flat += '/'+tn

        self._PROPS_flat = np.array([], dtype=np.float64)
        for p in self._PROPS:
            self._PROPS_flat = np.hstack((self._PROPS_flat, p))

        self._indDLoad_flat = np.array([], dtype=np.intp)
        for load in self._indDLoad:
            self._indDLoad_flat = np.hstack((self._indDLoad_flat, load))

        try:
            self._additionalLoadInfos_flat =  \
                np.concatenate(self._additionalLoadInfos)
        except:
            self._additionalLoadInfos_flat = np.array([])

        self._bc_target_flat = np.array([], dtype=np.intp)
        for b in self._bc_target:
            self._bc_target_flat = np.hstack((self._bc_target_flat, b))

        return None

    def _get_load_info(self):
        """Get attributes related to loads."""
        if self._nb_load > 0:
            loadinfos = [self._indDLoad_flat, self._JDLType, self._ADLMAG,
                         self._load_target_nbelem, self._nb_load,
                         self._additionalLoadInfos_flat,
                         self._nb_additionalLoadInfos]
            distribnodalload = []
            icount = 0
            for nodaldistrib in self._nodal_distributions.values():
                if nodaldistrib[1] == 'distributedload':
                    distribnodalload.append([nodaldistrib[0].T])
                    icount += 1
            distribnodalload = np.block(distribnodalload) \
                if icount > 0 else np.zeros((1, self._nb_cp))
            loadinfos.append(distribnodalload.copy('F'))
            return loadinfos
        else:
            i0 = np.zeros(1, dtype=np.intp)
            v0 = np.zeros(1, dtype=np.float64)
        return [i0, i0, v0, i0, 1, np.array([]), np.array([]), np.zeros((1, self._nb_cp))]

    def _get_bcs_info(self):
        """Get attributes related to boundary conditions."""
        if self._nb_bc > 0:
            return [self._bc_values, self._bc_target_flat,
                    self._bc_target_nbelem, self._nb_bc]
        i0 = np.zeros(1, dtype=np.intp)
        return [np.zeros((2, 1), dtype=np.float64), i0, i0, 1]

    """
    ---------------------------------------------------------------------------
    Les methodes suivantes renvoient les donnees necessaires pour les
    differentes analyses et les calculs des grandeurs d'interet.

    Principales methodes :
     * get_inputs4analysis() ........................ Calcul lineaire elastique
     * get_inputs4indDOF() ............ Determination des ddl bloques et libres
     * get_inputs4solution() ........... Reconstruction des Deplacements Nodaux
     * get_inputs4postproc() .. Postprocessing : ecriture de fichiers de sortie
     * get_inputs4vtu() .................... Fichier VTU pour le postprocessing
     * get_inputs4controlMesh() ............ Ficher VTK du polygone de controle
     * get_inputs4area() ..................... Calcul de l'aire de la structure
     * get_inputs4system() ............... Construction du systeme lineaire IGA
     * get_inputs4dlmmat() ........... Calcul Matrice Diagonale de Masse Lumped
     * get_inputs4strainExtrmm() ......... Calcul des extremum des deformations

    ---------------------------------------------------------------------------
    """

    #### Not used anymore ?
    # def get_inputs4analysis(self):
    #     """Return the necessary data for linear eleastic analysis."""
    #     inputs = [self._COORDS, self._IEN, self._elementsByPatch, self._Nkv,
    #               self._Ukv, self._Nijk, self._weight, self._Jpqr,
    #               self._ELT_TYPE, self._PROPS, self._JPROPS,
    #               self._MATERIAL_PROPERTIES, self._TENSOR, self._bc_target,
    #               self._indDLoad, self._JDLType, self._ADLMAG,
    #               self._bc_target_nbelem, self._load_target_nbelem,
    #               self._bc_values, self._mcrd, self._NBPINT, self._nb_bc,
    #               self._nb_load, self._nb_cload, self._nb_patch, self._nb_elem,
    #               self._nnode, self._nb_cp]

    #     return inputs

    def get_inputs4indDOF(self):
        """Get the settings concerning fixed and free degrees of freedom."""
        inputs = [self._mcrd, self._bc_target_flat, self._bc_target_nbelem,
                  self._bc_values, self._nb_dof_tot, self._nb_bc,
                  np.sum(self._bc_target_nbelem)]

        return inputs

    def get_inputs4solution(self, U_inv):
        """Return the necessary data for reconstructing nodal displacement."""
        COUPLG_flag = False
        bcs_infos = self._get_bcs_info()

        inputs = {'u_inv': U_inv,
                  'coords': self._COORDS,
                  'ind_dof_free': self._ind_dof_free,
                  'bc_target': bcs_infos[1].flatten(),
                  'bc_target_nbelem': bcs_infos[2],
                  'bc_values': bcs_infos[0],
                  'couplg_flag': COUPLG_flag,
                  'nb_bc': bcs_infos[-1],
                  'mcrd': self._mcrd,
                  'nb_cp': self._nb_cp,
                  'nb_dof_free': self._nb_dof_free}

        return inputs

    def get_inputs4eval_coupling_gpts(self, filename, npts_u, npts_v):
        """
        Return the necessary data to evaluate coupling Gauss points coordinates

        :param str filename: base name for resulting files
        :param int npts_u: number of points to compute along u parametric direction
        :param int npts_v: number of points to compute along v parametric direction
        :return: the necessary input parameters for Fortran subroutine ``eval_coupling_gpts``
        :rtype: dict
        """
        inputs = {'filename': filename,
                  'npts_u': npts_u,
                  'npts_v': npts_v,
                  'coords3d': self._COORDS,
                  'ien': self._IEN_flat,
                  'nb_elem_patch': self._elementsByPatch,
                  'nkv': self._Nkv,
                  'ukv': self._Ukv_flat,
                  'nijk': self._Nijk,
                  'weight': self._weight_flat,
                  'jpqr': self._Jpqr,
                  'elt_type': self._ELT_TYPE_flat,
                  'props': self._PROPS_flat,
                  'jprops': self._JPROPS,
                  'tensor': self._TENSOR_flat,
                  'mcrd': self._mcrd,
                  'nb_patch': self._nb_patch,
                  'nb_elem': self._nb_elem,
                  'nnode': self._nnode,
                  'nb_cp': self._nb_cp}

        return inputs

    def get_inputs4postprocVTU(self, filename, sol, nb_ref=np.ones(3, dtype=int),
                               Flag=np.array([True, True, True]), output_path='results'):
        """Get the necessary inputs for .vtu file generation.

        Parameters
        ----------
        filename : str
            Name of the output .vtu file.
        sol : numpy array
            Displacement field to plot. The function `get_inputs4solution` can
            be used for this purpose.
        nb_ref : numpy array, optional
            Mesh refinement for visualisation. It takes the form:
            np.array([xi-dir., eta-dir., zeta-dir.]. The default is np.ones(3).
        Flag : numpy array, optional
            Field ouput to save. Three are available: Nodal displacement,
            stresses, and Von Mises stresses.
            The default is np.array([True, True, True]).
        output_path: string
            Path to the directory where the output file should be written

        Returns
        -------
        inputs : list
            Necessary inputs for generating the .vtu file.
        """
        nb_ref = np.maximum(nb_ref, np.ones(3, dtype=int))
        inputs = [filename, output_path, Flag, nb_ref, sol,
                  self._COORDS, self._IEN_flat, self._elementsByPatch,
                  self._Nkv, self._Ukv_flat, self._Nijk, self._weight_flat,
                  self._Jpqr, self._ELT_TYPE_flat,
                  self._MATERIAL_PROPERTIES[:2, :], self._TENSOR_flat,
                  self._PROPS_flat, self._JPROPS, self._nnode, self._nb_patch,
                  self._nb_elem, self._nb_cp, self._mcrd]

        return inputs

    def get_inputs4post_curve_2D(self, i_patch, i_face, n_sample, sol):
        """
        Get the necessary inputs for function postproc.XXXXXXXXXXXXX in order to compute
        mechanical quantities (displacement and displacement gradient) along
        a curve for a 2D problem

        Parameters
        ----------
        i_patch : int
            Index of patch to process (starts at 1)
        i_face : int
            Index of side curve to process (1, 2, 3 or 4)
        n_sample : int
            Number of sample points to generate
        sol : numpy array
            Problem solution

        Returns
        -------
        inputs : dict
            Necessary inputs for function postproc.postproc_curve_2D
        """

        assert i_face > 0
        assert i_face < 5
        assert i_patch > 0

        if self._ELT_TYPE[i_patch-1] != 'U1':
            raise Exception('Element type ' +
                            self._ELT_TYPE[i_patch-1] +
                            'is not handled')

        inputs = {'sol': sol,
                  'n_sample': n_sample,
                  'i_patch': i_patch,
                  'i_face': i_face,
                  'ien': self._IEN_flat,
                  'props': self._PROPS_flat,
                  'jprops': self._JPROPS,
                  'nnode': self._nnode,
                  'nb_elem_patch': self._elementsByPatch,
                  'elt_type': self._ELT_TYPE_flat,
                  'tensor': self._TENSOR_flat,
                  'nkv': self._Nkv,
                  'jpqr': self._Jpqr,
                  'nijk': self._Nijk,
                  'ukv': self._Ukv_flat,
                  'weight': self._weight_flat,
                  'coords': self._COORDS
                  }

        return inputs

    def get_inputs4postproc_faces_vtu(self, filename, nb_ref=np.ones(3), output_path='results'):
        """
        Get the necessary inputs for function postproc.generate_faces_vtu

        Parameters
        ----------
        filename : str
            Name of the output .vtu file
        nb_ref : numpy array, optional
            Mesh refinement for the visualisation. It takes the form:
            np.array([xi-dir., eta-dir., zeta-dir.]). The default is np.ones(3)
        output_path: string
            Path to the directory where the output file should be written.
            The default is 'results'

        Returns
        -------
        inputs : dict
            Necessary inputs for function postproc.generate_faces_vtu
        """
        nb_ref = np.maximum(nb_ref, np.ones(3))
        inputs = { 'filename': filename,
                   'output_path': output_path,
                   'nb_refinement': nb_ref,
                   'coords3d': self._COORDS,
                   'ien': self._IEN_flat,
                   'nb_elem_patch': self._elementsByPatch,
                   'nkv': self._Nkv,
                   'ukv': self._Ukv_flat,
                   'nijk': self._Nijk,
                   'weight': self._weight_flat,
                   'jpqr': self._Jpqr,
                   'elt_type': self._ELT_TYPE_flat,
                   'tensor': self._TENSOR_flat,
                   'props': self._PROPS_flat,
                   'jprops': self._JPROPS,
                   'nnode': self._nnode,
                   'nb_patch': self._nb_patch,
                   'nb_elem': self._nb_elem,
                   'nb_cp': self._nb_cp,
                   'mcrd': self._mcrd }

        return inputs

    def get_inputs4postproc_bezier(self, i_patch, filename, sol, output_path='results'):
        """
        Get the necessary inputs for VTU postprocessing using Bezier cells
        with function postprocessing.postproc.generate_vtu_bezier
        Ouptut is made for a given patch of type U1

        Parameters
        ----------
        i_patch : int
            Index of patch to process (starts at 1)
        filename : str
            Name of output .vtu file
        sol : numpy array
            Displacement field to plot

        Returns
        -------
        inputs : dict
            Necessary inputs for VTU postprocessing using Bezier cells with
            function postprocessing.postproc.generate_vtu_bezier
        """

        if self._ELT_TYPE[i_patch-1] != 'U1':
            raise Exception('Element type ' +
                            self._ELT_TYPE[i_patch-1] +
                            'is not handled')


        # Make Bezier extraction
        CPs, Jpqr, ien, M = bezier_decomposition_patch(self, numpatch=i_patch-1,
                                                 return_ien=True,
                                                 return_mat=True)

        # Create an array containing local coordinates and weights of Bezier decomposition
        coords = np.zeros((4, M.shape[0]))
        for i_elem in range(ien.shape[0]):
            for i_cp in range(ien.shape[1]):
                coords[:, ien[i_elem, i_cp]] = CPs[i_elem, i_cp, :]

        # Convert solution for Bezier extraction of patch
        sol_patch_bezier = np.zeros((3, M.shape[0]))
        for i in range(3):
            sol_patch_bezier[i, :] = M @ sol[i, self._indCPbyPatch[i_patch-1]-1]

        # New inputs with built-in Bezier extraction
        inputs = {'filename': filename,
                  'output_path': output_path,
                  'i_patch': i_patch,
                  'sol': sol_patch_bezier,
                  'coords': coords[:3, :],
                  'weights': coords[3, :],
                  'ien': ien,
                  'jpqr': Jpqr
                  }

        return inputs

    def get_inputs4evaldisp(self, sol, xi, numpatch=1):
        """

        Parameters
        ----------
        sol : TYPE
            DESCRIPTION.
        xi : TYPE
            DESCRIPTION.
        numpatch : TYPE, optional
            DESCRIPTION. The default is 1.

        Returns
        -------
        inputs : TYPE
            DESCRIPTION.

        """
        if not (numpatch > 0 and numpatch < self._nb_patch + 1):
            print('Error: numpatch should be between 1'
                  'and {}'.format(self._nb_patch))
            return None

        inputs = [numpatch, xi, sol,
                  self._COORDS, self._IEN_flat, self._elementsByPatch,
                  self._Nkv, self._Ukv_flat, self._Nijk, self._weight_flat,
                  self._Jpqr, self._ELT_TYPE_flat,
                  self._MATERIAL_PROPERTIES[:2, :], self._TENSOR_flat,
                  self._PROPS_flat, self._JPROPS, self._nnode, self._nb_patch,
                  self._nb_elem, self._nb_cp, self._mcrd]

        return inputs

    def get_inputs4evaldispmulti(self, sol, n_xi, n_eta, n_zeta, numpatch=1):
        """

        Parameters
        ----------
        sol : TYPE
            DESCRIPTION.
        n_xi : TYPE
            DESCRIPTION.
        n_eta : TYPE
            DESCRIPTION.
        n_zeta : TYPE
            DESCRIPTION.
        numpatch : TYPE, optional
            DESCRIPTION. The default is 1.

        Returns
        -------
        inputs : TYPE
            DESCRIPTION.

        """
        if numpatch < 1 or numpatch > self._nb_patch:
            print('Error: numpatch should be between 1'
                  'and {}'.format(self._nb_patch))
            return None

        inputs = [numpatch, n_xi, n_eta, n_zeta, sol,
                  self._COORDS, self._IEN_flat, self._elementsByPatch,
                  self._Nkv, self._Ukv_flat, self._Nijk, self._weight_flat,
                  self._Jpqr, self._ELT_TYPE_flat,
                  self._MATERIAL_PROPERTIES[:2, :], self._TENSOR_flat,
                  self._PROPS_flat, self._JPROPS, self._nnode, self._nb_patch,
                  self._nb_elem, self._nb_cp, self._mcrd]

        return inputs

    def get_inputs4evalstress(self, sol, xi, numpatch=1):
        """

        Parameters
        ----------
        sol : TYPE
            DESCRIPTION.
        xi : TYPE
            DESCRIPTION.
        numpatch : TYPE, optional
            DESCRIPTION. The default is 1.

        Returns
        -------
        inputs : TYPE
            DESCRIPTION.

        """
        if not (numpatch > 0 and numpatch < self._nb_patch + 1):
            print('Error: numpatch should be between 1'
                  'and {}'.format(self._nb_patch))
            return None

        ntens = 2*self._mcrd
        inputs = [numpatch, xi, sol,
                  self._COORDS, self._IEN_flat, self._elementsByPatch,
                  self._Nkv, self._Ukv_flat, self._Nijk, self._weight_flat,
                  self._Jpqr, self._ELT_TYPE_flat,
                  self._MATERIAL_PROPERTIES[:2, :], self._TENSOR_flat,
                  self._PROPS_flat, self._JPROPS, ntens, self._nnode,
                  self._nb_patch, self._nb_elem, self._nb_cp, self._mcrd]

        return inputs

    def get_inputs4stressaggreg(self, sol, pnorm=10):
        """

        Parameters
        ----------
        sol : TYPE
            DESCRIPTION.
        pnorm : TYPE, optional
            DESCRIPTION. The default is 10.

        Returns
        -------
        inputs : TYPE
            DESCRIPTION.

        """
        inputs = [pnorm, sol,
                  self._COORDS, self._IEN_flat, self._elementsByPatch,
                  self._Nkv, self._Ukv_flat, self._Nijk, self._weight_flat,
                  self._Jpqr, self._ELT_TYPE_flat,
                  self._MATERIAL_PROPERTIES[:2, :], self._TENSOR_flat,
                  self._PROPS_flat, self._JPROPS, self._nnode, self._NBPINT,
                  self._nb_patch, self._nb_elem, self._nb_cp, self._mcrd]

        return inputs

    def get_inputs4stressaggregDiscrete(self, sol, ptseval, pnorm=10):
        """

        Parameters
        ----------
        sol : TYPE
            DESCRIPTION.
        ptseval : TYPE
            DESCRIPTION.
        pnorm : TYPE, optional
            DESCRIPTION. The default is 10.

        Returns
        -------
        inputs : TYPE
            DESCRIPTION.

        """
        if len(ptseval) != self._nb_patch:
            print('Error: wrong input ptseval. Should be a list with length'
                  'equals to nb_patch.')
            return None

        nb_pts_eval = np.zeros(self._nb_patch, dtype=np.intp)
        for i in np.arange(self._nb_patch):
            nb_pts_eval[i] = np.shape(ptseval[i])[1]
        nb_pts = np.sum(nb_pts_eval)
        ptsevalconc = np.zeros((3, nb_pts))
        j = 0
        for i in np.arange(self._nb_patch):
            if nb_pts_eval[i] > 0:
                ptsevalconc[:, j: j+nb_pts_eval[i]] = ptseval[i][:, :]
                j += nb_pts_eval[i]
        inputs = [pnorm, ptsevalconc, nb_pts_eval, sol,
                  self._COORDS, self._IEN_flat, self._elementsByPatch,
                  self._Nkv, self._Ukv_flat, self._Nijk, self._weight_flat,
                  self._Jpqr, self._ELT_TYPE_flat,
                  self._MATERIAL_PROPERTIES[:2, :], self._TENSOR_flat,
                  self._PROPS_flat, self._JPROPS, self._nnode, nb_pts,
                  self._nb_patch, self._nb_elem, self._nb_cp, self._mcrd]

        return inputs

    def get_inputs4adjointdisp(self, xi, numpatch=1):
        """

        Parameters
        ----------
        xi : TYPE
            DESCRIPTION.
        numpatch : TYPE, optional
            DESCRIPTION. The default is 1.

        Returns
        -------
        inputs : TYPE
            DESCRIPTION.

        """
        if not (numpatch > 0 and numpatch < self._nb_patch+1):
            print('Error: numpatch should be between 1 '
                  'and {}'.format(self._nb_patch))
            return None

        inputs = [numpatch, xi,
                  self._COORDS, self._IEN_flat, self._elementsByPatch,
                  self._Nkv, self._Ukv_flat, self._Nijk, self._weight_flat,
                  self._Jpqr, self._ELT_TYPE_flat,
                  self._MATERIAL_PROPERTIES[:2, :], self._TENSOR_flat,
                  self._PROPS_flat, self._JPROPS, self._nnode, self._mcrd,
                  self._nb_patch, self._nb_elem, self._nb_cp]

        return inputs

    def get_inputs4dispaggregDiscrete(self, sol, ptseval, pnorm=10):

        return self.get_inputs4stressaggregDiscrete(sol, ptseval, pnorm=pnorm)

    def get_inputs4regularizationaggreg(self, ptseval, pnorm=10):
        """

        Parameters
        ----------
        ptseval : TYPE
            DESCRIPTION.
        pnorm : TYPE, optional
            DESCRIPTION. The default is 10.

        Returns
        -------
        inputs : TYPE
            DESCRIPTION.

        """
        if len(ptseval) != self._nb_patch:
            print('Error: wrong input ptseval. Should be a list with length'
                  ' equals to nb_patch.')
            return None

        nb_pts_eval = np.zeros(self._nb_patch, dtype=np.intp)
        for i in np.arange(self._nb_patch):
            nb_pts_eval[i] = np.shape(ptseval[i])[1]
        nb_pts = np.sum(nb_pts_eval)
        ptsevalconc = np.zeros((3, nb_pts))
        j = 0
        for i in np.arange(self._nb_patch):
            if nb_pts_eval[i] > 0:
                ptsevalconc[:, j:j+nb_pts_eval[i]] = ptseval[i][:, :]
                j += nb_pts_eval[i]
        inputs = [pnorm, ptsevalconc, nb_pts_eval,
                  self._COORDS, self._IEN_flat, self._elementsByPatch,
                  self._Nkv, self._Ukv_flat, self._Nijk, self._weight_flat,
                  self._Jpqr, self._ELT_TYPE_flat,
                  self._MATERIAL_PROPERTIES[:2, :], self._TENSOR_flat,
                  self._PROPS_flat, self._JPROPS, self._nnode, nb_pts,
                  self._nb_patch, self._nb_elem, self._nb_cp]

        return inputs

    def get_inputs4postprocCPLG(self, filename, sol, nb_ref=5,
                                Flag=np.array([True, False, False])):
        """

        Parameters
        ----------
        filename : TYPE
            DESCRIPTION.
        sol : TYPE
            DESCRIPTION.
        nb_ref : TYPE, optional
            DESCRIPTION. The default is 5.
        Flag : TYPE, optional
            DESCRIPTION. The default is np.array([True, False, False]).

        Returns
        -------
        inputs : TYPE
            DESCRIPTION.

    """
        nb_ref = np.maximum(nb_ref, 2)
        inputs = [filename, Flag, nb_ref, sol,
                  self._COORDS, self._IEN_flat, self._elementsByPatch,
                  self._Nkv, self._Ukv_flat, self._Nijk, self._weight_flat,
                  self._Jpqr, self._ELT_TYPE_flat, self._PROPS_flat,
                  self._JPROPS, self._MATERIAL_PROPERTIES[:2, :],
                  self._TENSOR_flat, self._mcrd, self._NBPINT, self._nnode,
                  self._nb_patch, self._nb_elem, self._nb_cp]

        return inputs

    def get_inputs4postprocBasisFcts(self, filename, activepatch=None,
                                     dervorder=0, nb_pts=5):
        """

        Parameters
        ----------
        filename : TYPE
            DESCRIPTION.
        activepatch : TYPE, optional
            DESCRIPTION. The default is None.
        dervorder : TYPE, optional
            DESCRIPTION. The default is 0.
        nb_pts : TYPE, optional
            DESCRIPTION. The default is 5.

        Returns
        -------
        inputs : TYPE
            DESCRIPTION.

        """
        if activepatch is None or np.size(activepatch) != self._nb_patch:
            activepatch = np.ones(self._nb_patch, dtype=np.intp)
        inputs = [filename, activepatch, dervorder, nb_pts,
                  self._COORDS, self._IEN_flat, self._elementsByPatch,
                  self._Nkv, self._Ukv_flat, self._Nijk, self._weight_flat,
                  self._Jpqr, self._ELT_TYPE_flat, self._PROPS_flat,
                  self._JPROPS, self._MATERIAL_PROPERTIES[:2, :],
                  self._TENSOR_flat, self._mcrd, self._NBPINT, self._nnode,
                  self._nb_patch, self._nb_elem, self._nb_cp]

        return inputs

    def get_inputs4postprocPtsCloud(self, filename, activepatch=None,
                                    nb_pts=5):
        """

        Parameters
        ----------
        filename : TYPE
            DESCRIPTION.
        activepatch : TYPE, optional
            DESCRIPTION. The default is None.
        nb_pts : TYPE, optional
            DESCRIPTION. The default is 5.

        Returns
        -------
        inputs : TYPE
            DESCRIPTION.

        """
        if activepatch is None or np.size(activepatch) != self._nb_patch:
            activepatch = np.ones(self._nb_patch, dtype=np.intp)
        inputs = [filename, activepatch, np.maximum(2, nb_pts),
                  self._COORDS, self._IEN_flat, self._elementsByPatch,
                  self._Nkv, self._Ukv_flat, self._Nijk, self._weight_flat,
                  self._Jpqr, self._ELT_TYPE_flat, self._PROPS_flat,
                  self._JPROPS, self._MATERIAL_PROPERTIES[:2, :],
                  self._TENSOR_flat, self._mcrd, self._NBPINT, self._nnode,
                  self._nb_patch, self._nb_elem, self._nb_cp]

        return inputs

    def get_inputs4controlMesh(self, filename, num_patch, output_path='results', sol=None):
        """

        Parameters
        ----------
        filename : TYPE
            DESCRIPTION.
        num_patch : TYPE
            DESCRIPTION.
        sol : TYPE, optional
            DESCRIPTION. The default is None.

        Returns
        -------
        inputs : TYPE
            DESCRIPTION.

        """
        nb_cp_dir_thisPatch = np.maximum(
            self._Nkv[0:3, num_patch] - (self._Jpqr[0:3, num_patch]+1),
            np.ones(3))
        nb_cp_thisPatch = np.prod(nb_cp_dir_thisPatch)
        COORDS_thisPatch = self._COORDS[0:3, self._indCPbyPatch[num_patch] - 1]
        if sol is None:
            inputs = [filename, output_path, COORDS_thisPatch, nb_cp_dir_thisPatch,
                      nb_cp_thisPatch]
        else:
            sol_thisPatch = sol[:, self._indCPbyPatch[num_patch] - 1]
            inputs = [filename, output_path, COORDS_thisPatch, sol_thisPatch,
                      nb_cp_dir_thisPatch, nb_cp_thisPatch]

        return inputs

    def get_inputs4area(self, COORDS=None, activepatch=None):
        """

        Parameters
        ----------
        COORDS : TYPE, optional
            DESCRIPTION. The default is None.
        activepatch : TYPE, optional
            DESCRIPTION. The default is None.

        Returns
        -------
        inputs : TYPE
            DESCRIPTION.

        """
        if np.shape(COORDS) != np.shape(self._COORDS):
            COORDS = self._COORDS
        if activepatch is None or np.size(activepatch) != self._nb_patch:
            activepatch = np.ones(self._nb_patch, dtype=np.intp)
        inputs = [activepatch, COORDS, self._IEN_flat, self._elementsByPatch,
                  self._Nkv, self._Ukv_flat, self._Nijk, self._weight_flat,
                  self._Jpqr, self._ELT_TYPE_flat, self._PROPS_flat,
                  self._JPROPS, self._TENSOR_flat, self._mcrd, self._NBPINT,
                  self._nnode, self._nb_patch, self._nb_elem, self._nb_cp]

        return inputs

    def get_inputs4gradCompliance(self, sol, epsilon=1.e-4, COORDS=None,
                                  activeElem=None,
                                  activeDir=np.ones(3, np.intp)):
        """

        Parameters
        ----------
        sol : TYPE
            DESCRIPTION.
        epsilon : TYPE, optional
            DESCRIPTION. The default is 1.e-4.
        COORDS : TYPE, optional
            DESCRIPTION. The default is None.
        activeElem : TYPE, optional
            DESCRIPTION. The default is None.
        activeDir : TYPE, optional
            DESCRIPTION. The default is np.ones(3, np.intp).

        Returns
        -------
        inputs : TYPE
            DESCRIPTION.

        """
        if np.shape(COORDS) != np.shape(self._COORDS):
            COORDS = self._COORDS
        if np.size(activeElem) != self._nb_elem:
            activeElem = np.ones(self._nb_elem, dtype=np.intp)
        load_infos = self._get_load_info()
        inputs = [sol, epsilon, activeElem, activeDir, COORDS, self._IEN_flat,
                  self._elementsByPatch, self._Nkv, self._Ukv_flat, self._Nijk,
                  self._weight_flat, self._Jpqr, self._ELT_TYPE_flat,
                  self._PROPS_flat, self._JPROPS,
                  self._MATERIAL_PROPERTIES[:2, :], self._TENSOR_flat,
                  load_infos[0], load_infos[1], load_infos[2], load_infos[3],
                  self._mcrd, self._NBPINT, self._nnode, load_infos[4],
                  self._nb_patch, self._nb_elem, self._nb_cp]

        return inputs

    def get_inputs4gradVibration(self, eigenvect, eigenval, COORDS=None,
                                 activeElem=None,
                                 activeDir=np.ones(3, np.intp)):
        """

        Parameters
        ----------
        eigenvect : TYPE
            DESCRIPTION.
        eigenval : TYPE
            DESCRIPTION.
        COORDS : TYPE, optional
            DESCRIPTION. The default is None.
        activeElem : TYPE, optional
            DESCRIPTION. The default is None.
        activeDir : TYPE, optional
            DESCRIPTION. The default is np.ones(3, np.intp).

        Returns
        -------
        inputs : TYPE
            DESCRIPTION.

        """
        if np.shape(COORDS) != np.shape(self._COORDS):
            COORDS = self._COORDS
        if np.size(activeElem) != self._nb_elem:
            activeElem = np.ones(self._nb_elem, dtype=np.intp)
        nb_frq = np.size(eigenval)
        inputs = [eigenvect, eigenval, activeElem, activeDir, COORDS,
                  self._IEN_flat, self._elementsByPatch, self._Nkv,
                  self._Ukv_flat, self._Nijk, self._weight_flat, self._Jpqr,
                  self._ELT_TYPE_flat, self._PROPS_flat, self._JPROPS,
                  self._MATERIAL_PROPERTIES[:2, :],
                  self._MATERIAL_PROPERTIES[2, :], self._TENSOR_flat,
                  self._mcrd, self._NBPINT, self._nnode, nb_frq,
                  self._nb_patch, self._nb_elem, self._nb_cp]

        return inputs

    def get_inputs4gradCoupling(self, sol, COORDS=None, activepatch=None):
        """

        Parameters
        ----------
        sol : TYPE
            DESCRIPTION.
        COORDS : TYPE, optional
            DESCRIPTION. The default is None.
        activepatch : TYPE, optional
            DESCRIPTION. The default is None.

        Returns
        -------
        inputs : TYPE
            DESCRIPTION.

        """
        if np.shape(COORDS) != np.shape(self._COORDS):
            COORDS = self._COORDS
        if activepatch is None or np.size(activepatch) != self._nb_patch:
            activepatch = np.ones(self._nb_patch, dtype=np.intp)
        inputs = [activepatch, sol, COORDS, self._IEN_flat,
                  self._elementsByPatch, self._Nkv, self._Ukv_flat, self._Nijk,
                  self._weight_flat, self._Jpqr, self._ELT_TYPE_flat,
                  self._PROPS_flat, self._JPROPS, self._TENSOR_flat,
                  self._mcrd, self._NBPINT, self._nnode, self._nb_patch,
                  self._nb_elem, self._nb_cp]

        return inputs

    def get_inputs4gradDisplacement(self, sol, ADJ, COORDS=None,
                                    activeElem=None,
                                    activeDir=np.ones(3, np.intp)):
        """

        Parameters
        ----------
        sol : TYPE
            DESCRIPTION.
        ADJ : TYPE
            DESCRIPTION.
        COORDS : TYPE, optional
            DESCRIPTION. The default is None.
        activeElem : TYPE, optional
            DESCRIPTION. The default is None.
        activeDir : TYPE, optional
            DESCRIPTION. The default is np.ones(3, np.intp).

        Returns
        -------
        inputs : TYPE
            DESCRIPTION.

        """
        if np.shape(COORDS) != np.shape(self._COORDS):
            COORDS = self._COORDS
        if np.size(activeElem) != self._nb_elem:
            activeElem = np.ones(self._nb_elem, dtype=np.intp)
        nadj = np.size(ADJ, 0)
        load_infos = self._get_load_info()
        inputs = [sol, ADJ, activeElem, activeDir, COORDS,
                  self._IEN_flat, self._elementsByPatch, self._Nkv,
                  self._Ukv_flat, self._Nijk, self._weight_flat, self._Jpqr,
                  self._ELT_TYPE_flat, self._PROPS_flat, self._JPROPS,
                  self._MATERIAL_PROPERTIES[:2, :], self._TENSOR_flat,
                  load_infos[0], load_infos[1], load_infos[2], load_infos[3],
                  self._NBPINT, self._nnode, nadj, self._mcrd, load_infos[4],
                  self._nb_patch, self._nb_elem, self._nb_cp]

        return inputs

    def get_inputs4gradTotalwork(self, sol, ADJ, COORDS=None, activeElem=None,
                                 activeDir=np.ones(3, np.intp)):
        """

        Parameters
        ----------
        sol : TYPE
            DESCRIPTION.
        ADJ : TYPE
            DESCRIPTION.
        COORDS : TYPE, optional
            DESCRIPTION. The default is None.
        activeElem : TYPE, optional
            DESCRIPTION. The default is None.
        activeDir : TYPE, optional
            DESCRIPTION. The default is np.ones(3, np.intp).

        Returns
        -------
        inputs : TYPE
            DESCRIPTION.

        """
        if np.shape(COORDS) != np.shape(self._COORDS):
            COORDS = self._COORDS
        if np.size(activeElem) != self._nb_elem:
            activeElem = np.ones(self._nb_elem, dtype=np.intp)
        nadj = np.size(ADJ, 0)
        load_infos = self._get_load_info()
        inputs = [sol, ADJ, activeElem, activeDir, COORDS,
                  self._IEN_flat, self._elementsByPatch, self._Nkv,
                  self._Ukv_flat, self._Nijk, self._weight_flat, self._Jpqr,
                  self._ELT_TYPE_flat, self._PROPS_flat, self._JPROPS,
                  self._MATERIAL_PROPERTIES[:2, :],
                  self._MATERIAL_PROPERTIES[2, :], self._TENSOR_flat,
                  load_infos[0], load_infos[1], load_infos[2], load_infos[3],
                  load_infos[5], self._NBPINT, self._nnode, nadj, self._mcrd,
                  load_infos[4], self._nb_patch, self._nb_elem, self._nb_cp]

        return inputs

    def _get_inputs4gradTotalwork(self, sol, ADJ, COORDS=None, activeElem=None,
                                  activeDir=np.ones(3, np.intp),
                                  computeWint=True, computeWext=True):
        if np.shape(COORDS) != np.shape(self._COORDS):
            COORDS = self._COORDS
        if np.size(activeElem) != self._nb_elem:
            activeElem = np.ones(self._nb_elem, dtype=np.intp)
        nadj = np.size(ADJ, 0)
        load_infos = self._get_load_info()
        inputs = [computeWint, computeWext, sol, ADJ, activeElem, activeDir,
                  COORDS,
                  self._IEN_flat, self._elementsByPatch, self._Nkv,
                  self._Ukv_flat, self._Nijk, self._weight_flat, self._Jpqr,
                  self._ELT_TYPE_flat, self._PROPS_flat, self._JPROPS,
                  self._MATERIAL_PROPERTIES[:2, :],
                  self._MATERIAL_PROPERTIES[2, :], self._TENSOR_flat,
                  load_infos[0], load_infos[1], load_infos[2], load_infos[3],
                  load_infos[5], self._NBPINT, self._nnode, nadj, self._mcrd,
                  load_infos[4], self._nb_patch, self._nb_elem, self._nb_cp]

        return inputs

    def get_inputs4bandInfos(self):
        """SHORT DESCRIPTION."""
        inputs = [self._IEN_flat, self._elementsByPatch, self._Nkv,
                  self._Ukv_flat, self._Nijk, self._weight_flat, self._Jpqr,
                  self._ELT_TYPE_flat, self._PROPS_flat, self._JPROPS,
                  self._TENSOR_flat, self._ind_dof_free, self._nb_dof_free,
                  self._mcrd,  self._NBPINT, self._nnode, self._nb_patch,
                  self._nb_elem, self._nb_dof_tot]

        return inputs


    def get_inputs4system(self, COORDS=None):
        """SHORT DESCRIPTION."""
        load_infos = self._get_load_info()
        bcs_infos = self._get_bcs_info()

        if np.shape(COORDS) != np.shape(self._COORDS):
            COORDS = self._COORDS

        inputs = [COORDS, self._IEN_flat, self._elementsByPatch, self._Nkv,
                  self._Ukv_flat, self._Nijk, self._weight_flat, self._Jpqr,
                  self._ELT_TYPE_flat, self._PROPS_flat, self._JPROPS,
                  self._MATERIAL_PROPERTIES[:2, :], self._TENSOR_flat,
                  load_infos[0], load_infos[1], load_infos[2], load_infos[3],
                  bcs_infos[0], bcs_infos[1], bcs_infos[2], self._ind_dof_free,
                  self._nb_dof_free, self._mcrd,  self._NBPINT, self._nnode,
                  bcs_infos[-1], load_infos[4], self._nb_patch,self._nb_elem,
                  self._nb_cp, self._nb_dof_tot]

        return inputs

    def get_inputs4system_elemStorage(self, COORDS=None, activeElem=None):
        """Return parameters to compute the sparse stiffness matrix.

        Data is structured to comply with the function input parameters. The
        functions is `stiffmtrx_elemstorage.sys_linmat_lindef_static`

        Parameters
        ----------
        COORDS : numpy array, optional
            Control points coordinates. The default is None. Default
            correponds to the coordinates of the current geometry.
        activeElem : numpy array (dtype=intp), optional
            Array indicating elements considered to build the stiffness matrix.
            Value `1` stands for active elements. The default is None (all
            elements are activated).

        Returns
        -------
        inputs : list
            All necessary inputs to build stiffness matrix.

        """
        if activeElem is None:
            activeElem = np.ones(self._nb_elem, dtype=np.intp)
            indpatch2rm = np.where(np.all(np.array([self._ELT_TYPE != 'U1',
                                                    self._ELT_TYPE != 'U99',
                                                    self._ELT_TYPE != 'U98',
                                                    self._ELT_TYPE != 'U2',
                                                    self._ELT_TYPE != 'U3',
                                                    self._ELT_TYPE != 'U10',
                                                    self._ELT_TYPE != 'U30']),
                                          axis=0))[0]
            e = np.cumsum(np.concatenate(([0], self._elementsByPatch)),
                          dtype=np.intp)
            for patch in indpatch2rm:
                activeElem[e[patch]: e[patch+1]] = 0

        ndofel = self._nnode*self._mcrd
        ndofel_list = np.repeat(ndofel*(ndofel + 1) / 2, self._elementsByPatch)
        nb_data = np.sum(activeElem*ndofel_list)

        load_infos = self._get_load_info()
        bcs_infos = self._get_bcs_info()
        if np.shape(COORDS) != np.shape(self._COORDS):
            COORDS = self._COORDS

        inputs = [activeElem, nb_data, COORDS, self._IEN_flat,
                  self._elementsByPatch, self._Nkv, self._Ukv_flat,
                  self._Nijk, self._weight_flat, self._Jpqr,
                  self._ELT_TYPE_flat, self._PROPS_flat, self._JPROPS,
                  self._MATERIAL_PROPERTIES[:2, :], self._N_MATERIAL_PROPERTIES[:],
                  self._MATERIAL_PROPERTIES[2, :], self._TENSOR_flat,
                  load_infos[0], load_infos[1], load_infos[2], load_infos[3],
                  load_infos[5], load_infos[6], bcs_infos[0], bcs_infos[1], bcs_infos[2],
                  self._ind_dof_free, self._nb_dof_free, self._mcrd,
                  self._NBPINT, self._nnode, load_infos[-1], bcs_infos[-1],
                  load_infos[4], self._nb_patch, self._nb_elem, self._nb_cp,
                  self._nb_dof_tot]

        return inputs

    def get_inputs4massmat(self, COORDS=None, activeElem=None):
        """

        Parameters
        ----------
        COORDS : TYPE, optional
            DESCRIPTION. The default is None.
        activeElem : TYPE, optional
            DESCRIPTION. The default is None.

        Returns
        -------
        inputs : TYPE
            DESCRIPTION.

        """
        if activeElem == None:
            activeElem = np.ones(self._nb_elem, dtype=np.intp)
            indpatch2rm = np.where(np.all(np.array([self._ELT_TYPE != 'U1',
                                                    self._ELT_TYPE != 'U2',
                                                    self._ELT_TYPE != 'U3',
                                                    self._ELT_TYPE != 'U10',
                                                    self._ELT_TYPE != 'U30']),
                                          axis=0))[0]
            e = np.cumsum(np.concatenate(([0], self._elementsByPatch)),
                          dtype=np.intp)
            for patch in indpatch2rm:
                activeElem[e[patch]: e[patch+1]] = 0

        ndofel = self._nnode*self._mcrd
        ndofel_list = np.repeat(ndofel*(ndofel + 1) / 2, self._elementsByPatch)
        nb_data = int(np.sum(activeElem*ndofel_list))

        if np.shape(COORDS) != np.shape(self._COORDS):
            COORDS = self._COORDS

        inputs = [activeElem, nb_data, COORDS, self._IEN_flat,
                  self._elementsByPatch, self._Nkv, self._Ukv_flat, self._Nijk,
                  self._weight_flat, self._Jpqr, self._ELT_TYPE_flat,
                  self._PROPS_flat, self._JPROPS,
                  self._MATERIAL_PROPERTIES[2, :], self._TENSOR_flat,
                  self._ind_dof_free, self._nb_dof_free, self._mcrd,
                  self._NBPINT, self._nnode, self._nb_patch, self._nb_elem,
                  self._nb_cp, self._nb_dof_tot]

        return inputs

    def get_inputs4grammat(self):
        """
        Generate input for the routine computing Gram matrix
        """
        activeElem = np.ones(self._nb_elem, dtype=np.intp)
        indpatch2rm = np.where(np.all(np.array([self._ELT_TYPE != 'U1',
                                                self._ELT_TYPE != 'U2',
                                                self._ELT_TYPE != 'U3',
                                                self._ELT_TYPE != 'U10',
                                                self._ELT_TYPE != 'U30']),
                                      axis=0))[0]
        e = np.cumsum(np.concatenate(([0], self._elementsByPatch)),
                      dtype=np.intp)
        for patch in indpatch2rm:
            activeElem[e[patch]: e[patch+1]] = 0

        nnode_list = np.repeat(
            self._nnode*(self._nnode + 1) / 2, self._elementsByPatch)
        nb_data = np.sum(activeElem*nnode_list)

        inputs = [activeElem, nb_data, self._COORDS, self._IEN_flat,
                  self._elementsByPatch, self._Nkv, self._Ukv_flat,
                  self._Nijk, self._weight_flat, self._Jpqr,
                  self._ELT_TYPE_flat, self._PROPS_flat, self._JPROPS,
                  self._TENSOR_flat,self._ind_dof_free, self._mcrd,
                  self._NBPINT, self._nnode, self._nb_dof_tot]

        return inputs

    def get_inputs4svarsrhs(self, sol):
        """
        Generate input for the routine computing RHS vectors of the least square
        projection of values at Gauss points

        Parameters
        ----------
        sol : np.array
            displacement solution at control points
        """
        activeElem = np.ones(self._nb_elem, dtype=np.intp)
        indpatch2rm = np.where(np.all(np.array([self._ELT_TYPE != 'U1',
                                                self._ELT_TYPE != 'U2',
                                                self._ELT_TYPE != 'U3',
                                                self._ELT_TYPE != 'U10',
                                                self._ELT_TYPE != 'U30']),
                                      axis=0))[0]
        e = np.cumsum(np.concatenate(([0], self._elementsByPatch)),
                      dtype=np.intp)
        for patch in indpatch2rm:
            activeElem[e[patch]: e[patch+1]] = 0

        inputs = [sol, activeElem, self._COORDS, self._IEN_flat,
                  self._elementsByPatch, self._Nkv, self._Ukv_flat, self._Nijk,
                  self._weight_flat, self._Jpqr, self._ELT_TYPE_flat,
                  self._TENSOR_flat, self._PROPS_flat, self._JPROPS,
                  self._MATERIAL_PROPERTIES[:2, :], self._NBPINT, self._nnode]

        return inputs

    def get_inputs4proj_vtu(self, filename, sol, svars, nb_ref=np.ones(3), output_path='results'):
        """Return the inputs for least squer projection in a VTU file.

        Parameters
        ----------
        filename: string
            Name of the file to write without extension
        sol: numpy.array
            Displacement field solution
        svars: numpy.array
            variables originally defined at integration points and mean-square
            projected on CP
        nb_ref: numpy.array
            Mesh refinement level in each direction. Default is [1, 1, 1]
        outout_path: string
            path to output directory, default is 'results'

        Returns
        -------
        inputs : list
            list of input parameters for pp.generate_pro_vtu
        """

        nb_ref = np.maximum(nb_ref, np.ones(3))
        inputs = [filename, output_path, nb_ref, sol, svars, self._COORDS, self._IEN_flat,
                  self._elementsByPatch, self._Nkv, self._Ukv_flat, self._Nijk,
                  self._weight_flat, self._Jpqr, self._ELT_TYPE_flat,
                  self._TENSOR, self._PROPS_flat, self._JPROPS, self._nnode]

        return inputs

    def get_inputs4postproc_cplg_vtu(self, filename, lgrge_patch_number, sol,
                                     nb_ref=np.ones(2), output_path='results'):
        """
        Returns the inputs for the generation of VTU file for results
        at coupling interface

        Parameters
        ----------
        filename : string
            Name of the file to write without extension
        lgrge_patch_number : int
            Index of Lagrange patch number of the interface to process
        sol : numpy.array
            Displacement field solution
        nb_ref : numpy.array
            mesh refinement level in each direction. default is [1, 1]

        Returns
        -------
        inputs : list
            List of input parameters for pp.generate_coupling_vtu
        """

        nb_ref = np.maximum(nb_ref, np.ones(2))
        inputs = [filename, output_path, lgrge_patch_number, nb_ref, sol, self._COORDS,
                  self._IEN_flat, self._elementsByPatch, self._Nkv,
                  self._Ukv_flat, self._Nijk, self._weight_flat, self._Jpqr,
                  self._ELT_TYPE_flat, self._TENSOR_flat, self._PROPS_flat,
                  self._JPROPS, self._nnode, self._mcrd
                  ]

        return inputs

    def get_inputs4cplgmatrix(self, COORDS=None):
        """SHORT DESCRIPTION."""
        ndofel = self._nnode*self._mcrd
        ind = np.where(self._ELT_TYPE == 'U00')[0]
        nb_data = 0
        for c in ind:
            ph = np.intp(self._PROPS[c][1]-1)
            lg = np.intp(self._PROPS[c][2]-1)
            nb_data += self._elementsByPatch[c] * self._NBPINT[c] * \
                ndofel[ph] * ndofel[lg]

        if np.shape(COORDS) != np.shape(self._COORDS):
            COORDS = self._COORDS

        inputs = [nb_data, COORDS, self._IEN_flat, self._elementsByPatch,
                  self._Nkv, self._Ukv_flat, self._Nijk, self._weight_flat,
                  self._Jpqr, self._ELT_TYPE_flat, self._PROPS_flat,
                  self._JPROPS, self._MATERIAL_PROPERTIES[:2, :],
                  self._TENSOR_flat, self._ind_dof_free, self._nb_dof_free,
                  self._mcrd,self._NBPINT, self._nnode, self._nb_patch,
                  self._nb_elem, self._nb_cp, self._nb_dof_tot]

        return inputs

    def get_inputs4cplgmatrixU5(self, COORDS=None, integrationOrder=0,
                                output_path='results'):
        """
        Generate input arguments for ``coupling.cplgmatrix.cplg_matrixu5``.

        Parameters
        ----------
        COORDS : list, optional
            Control points coordinates. The default is None.
        integrationOrder : int, optional
            Integration order for coupling equation (used only if greater than
            auto-evaluated one). The default is 0.
        output_path : string
            Path to directory where debug files (for projection) should be
            written.

        Returns
        -------
        inputs : list
            List of parameters for ``coupling.cplgmatrix.cplg_matrixu5``.

        """
        ind = np.where(self._ELT_TYPE == 'U5')[0]
        nb_data = 0
        for c in ind:
            lg = np.intp(self._PROPS[c][0] - 1)   # patch lagrange
            pm = np.intp(self._PROPS[c][1] - 1)   # patch master
            fm = np.intp(self._PROPS[c][2])       # face master
            ps = np.intp(self._PROPS[c][3] - 1)   # patch slave
            fs = np.intp(self._PROPS[c][4])       # face slave

            # TODO : differenciate 2D and 3D cases in a more stylish way.

            # Number of elements/face - master
            if self._dim[c] == 2:       # 3D case : interface has dimension 2
                if fm in (1, 2):
                    nb_el_m = (self._Nkv[1, pm] - 2 * self._Jpqr[1, pm] - 1) * \
                        (self._Nkv[2, pm] - 2 * self._Jpqr[2, pm] - 1)
                elif fm in (3, 4):
                    nb_el_m = (self._Nkv[0, pm] - 2 * self._Jpqr[0, pm] - 1) * \
                        (self._Nkv[2, pm] - 2 * self._Jpqr[2, pm] - 1)
                elif fm in (5, 6):
                    nb_el_m = (self._Nkv[0, pm] - 2 * self._Jpqr[0, pm] - 1) * \
                        (self._Nkv[1, pm] - 2 * self._Jpqr[1, pm] - 1)
            elif self._dim[c] == 1:     # 2D case : interface has dimension 1
                if fm == 1 or fm == 2:
                    nb_el_m = (self._Nkv[1, pm] - 2 * self._Jpqr[1, pm] - 1)
                elif fm == 3 or fm == 4:
                    nb_el_m = (self._Nkv[0, pm] - 2 * self._Jpqr[0, pm] - 1)

            # Number of elements/face - slave
            if self._dim[c] == 2:       # 3D case : interface has dimension 2
                if fs in (1, 2):
                    nb_el_s = (self._Nkv[1, ps] - 2 * self._Jpqr[1, ps] - 1) * \
                        (self._Nkv[2, ps] - 2 * self._Jpqr[2, ps] - 1)
                elif fs in (3, 4):
                    nb_el_s = (self._Nkv[0, ps] - 2 * self._Jpqr[0, ps] - 1) * \
                        (self._Nkv[2, ps] - 2 * self._Jpqr[2, ps] - 1)
                elif fs in (5, 6):
                    nb_el_s = (self._Nkv[0, ps] - 2 * self._Jpqr[0, ps] - 1) * \
                        (self._Nkv[1, ps] - 2 * self._Jpqr[1, ps] - 1)
            elif self._dim[c] == 1:     # 2D case : interface has dimension 1
                if fs in (1, 2):
                    nb_el_s = (self._Nkv[1, ps] - 2 * self._Jpqr[1, ps] - 1)
                elif fs in (3, 4):
                    nb_el_s = (self._Nkv[0, ps] - 2 * self._Jpqr[0, ps] - 1)

            # Number of elements/face - Lagrange
            if self._dim[c] == 2:       # 3D case : interface has dimension 2
                nb_el_l = (self._Nkv[0, lg] - 2 * self._Jpqr[0, lg] - 1) * \
                    (self._Nkv[1, lg] - 2 * self._Jpqr[1, lg] - 1)
            elif self._dim[c] == 1:     # 2D case : interface has dimension 1
                nb_el_l = (self._Nkv[0, lg] - 2 * self._Jpqr[0, lg] - 1)

            # Compute integration order
            order = max(self._Jpqr[:, lg]) + \
                max(max(self._Jpqr[:, pm]), max(self._Jpqr[:, ps]))
            if order < integrationOrder:
                order = integrationOrder
            # nb_data += (order**2) * \
            #     (nb_el_m * self._nnode[pm] + nb_el_s * self._nnode[ps]) * \
            #     nb_el_l * self._nnode[lg] * self._mcrd

            # TODO : verify if number of elements taken into account is the good
            # one in the following line (use nb_el_l instead ???)
            nb_gps = (order**self._dim[c]) * nb_el_m
            nb_data += nb_gps * self._mcrd * self._nnode[lg] * \
                (self._nnode[ps]+self._nnode[pm])

        if np.shape(COORDS) != np.shape(self._COORDS):
            COORDS = self._COORDS

        inputs = [output_path, nb_data, COORDS, self._IEN_flat, self._elementsByPatch,
                  self._Nkv, self._Ukv_flat, self._Nijk, self._weight_flat,
                  self._Jpqr, self._ELT_TYPE_flat, self._PROPS_flat,
                  self._JPROPS, self._MATERIAL_PROPERTIES[:2, :],
                  self._TENSOR_flat, self._ind_dof_free, self._nb_dof_free,
                  self._mcrd, self._NBPINT, self._nnode, integrationOrder]

        return inputs

    def get_inputs4dirichletmatrix(self):
        """SHORT DESCRIPTION."""
        nb_data = self._nb_dof_bloq
        inputs = [nb_data, self._mcrd, self._bc_target_flat,
                  self._bc_target_nbelem, self._bc_values, self._nb_dof_tot,
                  self._nb_bc, np.sum(self._bc_target_nbelem)]

        return inputs

    def get_inputs4interpolation(self, activePatch=None, listParam=None):
        """

        Parameters
        ----------
        activePatch : TYPE, optional
            DESCRIPTION. The default is None.
        listParam : TYPE, optional
            DESCRIPTION. The default is None.

        Returns
        -------
        inputs : TYPE
            DESCRIPTION.

        """
        if activePatch is None:
            activePatch = np.ones(self._nb_patch, dtype=np.intp)

        cpDir = np.maximum(self._Nkv - (self._Jpqr + 1), 0).astype(np.intp)
        if listParam is None:
            # greville
            nbdataDir = np.maximum(
                (self._Jpqr + 1) * (cpDir - 2) + 2, 0).astype(np.intp)
            nbdata = np.sum(nbdataDir[:, activePatch == 1])
            nbrow = np.sum(cpDir[:, activePatch == 1]) + 1
            inputs = [nbdata, nbrow, activePatch, self._IEN_flat,
                      self._elementsByPatch, self._Nkv, self._Ukv_flat,
                      self._Nijk, self._weight_flat, self._Jpqr,
                      self._ELT_TYPE_flat, self._TENSOR_flat, self._PROPS_flat,
                      self._JPROPS, self._nnode, self._nb_patch, self._nb_elem]
        else:
            # user parameters
            xipara = np.array([], dtype=np.float64)
            nbactivedir = len(listParam)
            nb_rowByDir = np.zeros(nbactivedir, dtype=np.intp)
            i = 0
            for xi in listParam:
                xipara = np.hstack((xipara, xi))
                nb_rowByDir[i] = xi.size
                i += 1
            temp = np.ravel(self._Jpqr[:, activePatch == 1], order='F')
            nbdata = np.sum(
                np.multiply(nb_rowByDir, temp[np.nonzero(temp)] + 1))
            nbrow = xipara.size
            inputs = [nbdata, activePatch, xipara, nb_rowByDir, self._IEN_flat,
                      self._elementsByPatch, self._Nkv, self._Ukv_flat,
                      self._Nijk, self._weight_flat, self._Jpqr,
                      self._ELT_TYPE_flat, self._TENSOR_flat, self._PROPS_flat,
                      self._JPROPS, self._nnode, nbrow, nbactivedir,
                      self._nb_patch, self._nb_elem]

        return inputs

    def get_inputs4grevphyscoords(self, activePatch=None):
        """SHORT DESCRIPTION."""
        if activePatch is None:
            activePatch = np.ones(self._nb_patch, dtype=np.intp)
        cpDir = np.maximum(self._Nkv - (self._Jpqr + 1), 1).astype(np.intp)
        nb_grev = np.sum(np.prod(cpDir, axis=0)[activePatch == 1])
        inputs = [nb_grev, activePatch, self._COORDS, self._IEN_flat,
                  self._elementsByPatch, self._Nkv, self._Ukv_flat, self._Nijk,
                  self._weight_flat, self._Jpqr, self._ELT_TYPE_flat,
                  self._TENSOR_flat, self._PROPS_flat, self._JPROPS,
                  self._nnode, self._nb_patch, self._nb_elem, self._nb_cp]

        return inputs

    def get_inputs4cplginfos(self, nb_pts=3):
        """SHORT DESCRIPTION."""
        ind = np.where(self._ELT_TYPE == 'U00')[0]
        nb_interface = ind.size
        dim_interface = np.max(self._dim[ind])
        nb_data = nb_pts**dim_interface*3 + 4
        inputs = [nb_pts, nb_data, nb_interface, self._COORDS, self._IEN_flat,
                  self._elementsByPatch, self._Nkv, self._Ukv_flat, self._Nijk,
                  self._weight_flat, self._Jpqr, self._ELT_TYPE_flat,
                  self._PROPS_flat, self._JPROPS,
                  self._MATERIAL_PROPERTIES[:2, :], self._TENSOR_flat,
                  self._mcrd, self._NBPINT, self._nnode, self._nb_patch,
                  self._nb_elem, self._nb_cp]

        return inputs

    def get_inputs4system_bandStorage(self, nb_diag=None, NUMDof2Diag=None,
                                      COORDS=None, activeElem=None):
        """

        Parameters
        ----------
        nb_diag : TYPE, optional
            DESCRIPTION. The default is None.
        NUMDof2Diag : TYPE, optional
            DESCRIPTION. The default is None.
        COORDS : TYPE, optional
            DESCRIPTION. The default is None.
        activeElem : TYPE, optional
            DESCRIPTION. The default is None.

        Returns
        -------
        inputs : TYPE
            DESCRIPTION.

        """
        if not nb_diag:
            nb_diag, NUMDiag2Dof, NUMDof2Diag = \
                bs.bandstorageinfo(*self.get_inputs4bandInfos())
        if activeElem is None:
            activeElem = np.ones(self._nb_elem, dtype=np.intp)

        load_infos = self._get_load_info()
        bcs_infos = self._get_bcs_info()
        if np.shape(COORDS) != np.shape(self._COORDS):
            COORDS = self._COORDS
        inputs = [nb_diag, NUMDof2Diag, activeElem,
                  COORDS, self._IEN_flat, self._elementsByPatch, self._Nkv,
                  self._Ukv_flat, self._Nijk, self._weight_flat, self._Jpqr,
                  self._ELT_TYPE_flat, self._PROPS_flat, self._JPROPS,
                  self._MATERIAL_PROPERTIES[:2,:], self._TENSOR_flat, load_infos[0],
                  load_infos[1], load_infos[2], load_infos[3], bcs_infos[0],
                  bcs_infos[1], bcs_infos[2], self._ind_dof_free, self._nb_dof_free,
                  self._mcrd,  self._NBPINT, self._nnode, bcs_infos[-1], load_infos[4],
                  self._nb_patch,self._nb_elem, self._nb_cp, self._nb_dof_tot]

        return inputs

    def get_inputs4entities2param(self):
        """SHORT DESCRIPTION."""
        numEntities = np.where(self._ELT_TYPE == 'U31')
        nb_points = np.sum(
            self._elementsByPatch[numEntities]*self._NBPINT[numEntities])
        inputs = [nb_points, self._COORDS, self._IEN_flat,
                  self._elementsByPatch, self._Nkv, self._Ukv_flat,self._Nijk,
                  self._weight_flat, self._Jpqr, self._ELT_TYPE_flat,
                  self._MATERIAL_PROPERTIES[:2, :], self._TENSOR_flat,
                  self._PROPS_flat, self._JPROPS,self._nnode, self._NBPINT,
                  self._mcrd, self._nb_patch, self._nb_elem, self._nb_cp]

        return inputs

    def get_inputs4cplgEmbddEntites(self, nb_diag=None, NUMDof2Diag=None):
        """

        Parameters
        ----------
        nb_diag : TYPE, optional
            DESCRIPTION. The default is None.
        NUMDof2Diag : TYPE, optional
            DESCRIPTION. The default is None.

        Returns
        -------
        inputs : TYPE
            DESCRIPTION.

        """
        if not nb_diag:
            nb_diag,NUMDiag2Dof,NUMDof2Diag = bs.bandstorageinfo( *self.get_inputs4bandInfos() )
        inputs = [nb_diag,NUMDof2Diag,
                  self._COORDS, self._IEN_flat, self._elementsByPatch, self._Nkv, self._Ukv_flat,
                  self._Nijk,self._weight_flat, self._Jpqr, self._ELT_TYPE_flat, self._PROPS_flat,
                  self._JPROPS, self._MATERIAL_PROPERTIES[:2,:], self._TENSOR_flat, self._mcrd,
                  self._NBPINT, self._nnode, self._nb_patch,self._nb_elem, self._nb_cp,
                  self._nb_dof_tot]

        return inputs



    def get_inputs4dlmmat(self):
        """SHORT DESCRIPTION."""
        inputs = [self._COORDS, self._IEN, self._elementsByPatch, self._Nkv, self._Ukv,
                  self._Nijk, self._weight, self._Jpqr, self._ELT_TYPE, self._mcrd, self._NBPINT,
                  self._nb_patch, self._nb_elem, self._nnode, self._nb_cp]

        return inputs

    def get_inputs4geomat(self, sol, case=1):
        """

        Parameters
        ----------
        sol : TYPE
            DESCRIPTION.
        case : TYPE, optional
            DESCRIPTION. The default is 1.

        Returns
        -------
        inputs : TYPE
            DESCRIPTION.

        """
        if case==1:
            inputs = [self._COORDS, sol, self._IEN, self._elementsByPatch, self._ind_dof_free,
                      self._Nkv, self._Ukv, self._Nijk, self._weight, self._Jpqr,
                      self._ELT_TYPE, self._PROPS, self._JPROPS, self._MATERIAL_PROPERTIES,
                      self._TENSOR, self._nb_dof_free, self._NBPINT, self._mcrd,
                      self._nb_patch, self._nb_elem, self._nnode, self._nb_cp, self._nb_dof_tot]
        else:
            inputs = [self._COORDS, sol, self._IEN, self._elementsByPatch, self._ind_dof_free,
                      self._Nkv, self._Ukv, self._Nijk, self._weight, self._Jpqr,
                      self._ELT_TYPE, self._PROPS, self._JPROPS, self._MATERIAL_PROPERTIES,
                      self._TENSOR, self._bc_target, self._bc_target_nbelem, self._bc_values,
                      self._nb_dof_free, self._NBPINT, self._nb_bc, self._mcrd, self._nb_patch,
                      self._nb_elem, self._nnode, self._nb_cp, self._nb_dof_tot]

        return inputs


    def get_inputs4strainExtrmm(self,sol):
        """SHORT DESCRIPTION."""
        inputs = [self._COORDS, self._IEN, self._elementsByPatch, sol, self._Nkv, self._Ukv,
                  self._Nijk, self._weight, self._Jpqr, self._ELT_TYPE, self._NBPINT, self._nnode,
                  self._nb_patch, self._nb_elem, self._nb_cp, self._mcrd]

        return inputs


    def get_inputs4gauss3D(self,KNumFace=0):
        """SHORT DESCRIPTION."""
        nb_gauss_tot = 0
        j=0 if KNumFace==0 else 1
        for i in np.arange(0,self._nb_patch):
            NbPtInt = np.int(np.int(self._NBPINT[i])**(1./self._dim[i]))
            if NbPtInt**self._dim[i]<np.int(self._NBPINT[i]):
                NbPtInt += 1
            nb_gauss_tot += NbPtInt**(self._dim[i]-j) * self._elementsByPatch[i]
        inputs = [self._COORDS, self._IEN_flat, self._elementsByPatch, self._Nkv,
                  self._Ukv_flat, self._Nijk, self._weight_flat, self._Jpqr,
                  self._ELT_TYPE_flat, self._PROPS_flat, self._JPROPS,
                  self._TENSOR_flat, KNumFace, self._NBPINT, nb_gauss_tot,
                  self._mcrd, self._nnode, self._nb_patch, self._nb_elem,
                  self._nb_cp]

        return inputs


    def get_inputs4lagrangeCplgMatrx(self, Mult, IGJELEM, GaussCoordsGlo):
        """

        Parameters
        ----------
        Mult : TYPE
            DESCRIPTION.
        IGJELEM : TYPE
            DESCRIPTION.
        GaussCoordsGlo : TYPE
            DESCRIPTION.

        Returns
        -------
        inputs : TYPE
            DESCRIPTION.

        """
        nb_gauss_tot = np.size(GaussCoordsGlo,1)
        inputs = [self._ind_dof_free, self._nb_dof_free, IGJELEM, GaussCoordsGlo,
                  self._COORDS, self._IEN_flat, self._elementsByPatch, self._Nkv, self._Ukv_flat,
                  self._Nijk, self._weight_flat, self._Jpqr, self._ELT_TYPE_flat, self._PROPS_flat,
                  self._JPROPS, self._TENSOR_flat, self._nnode, self._mcrd,
                  Mult._COORDS, Mult._IEN_flat, Mult._elementsByPatch, Mult._Nkv, Mult._Ukv_flat,
                  Mult._Nijk, Mult._weight_flat, Mult._Jpqr, Mult._ELT_TYPE_flat, Mult._PROPS_flat,
                  Mult._JPROPS, Mult._TENSOR_flat, Mult._nnode, Mult._NBPINT, Mult._dim[0],
                  self._nb_dof_tot,nb_gauss_tot,self._nb_cp,self._nb_patch,self._nb_elem,
                  Mult._nb_cp,Mult._nb_patch,Mult._nb_elem]

        return inputs

    def refine_and_getTransformationMatrices(
            self, nb_refinementByDirection,
            nb_degreeElevationByDirection=[0, 0, 0], additional_knots=None):
        """Refinement by degree elevation and knot insertion.

        For knot insertion, a level of refinement consists in inserting all of
        the average knots of non-zero intervals. The number of elements is
        hence doubled.

        Parameters
        ----------
        nb_refinementByDirection : list of ints
            Level of refinement by direction.
                It takes the form: [xi-dir., eta-dir., zeta-dir.].
        nb_degreeElevationByDirection : list of ints, optional
            Number of degree elevation by direction. The default is [0, 0, 0].
                It takes the form:[xi-dir., eta-dir., zeta-dir.].

        additional_knots : dict, optional
            Dictionnary for specific knot insertion. The default is None.
              The dictionay takes the following keywords: `patches`, `1`, `2`,
              `3`.
            Example:
                Refine the patch indexed as `2` at parameters u=0.5 and w=0.8.
                Hence we have: `additional_knots = {'patches': np.array([2]),
                '1': np.array([0.5]), '2': np.array([]),
                '3': np.array([0.8])}`.

        Returns
        -------
        transformationMatrices : numpy array
            Resulting transformation matrices for refinement.
        """
        # Initialization refinement inputs
        # - Level of knot insertions
        if np.size(np.shape(nb_refinementByDirection)) == 1:
            nb_refinementByDirection = np.append(
                nb_refinementByDirection,
                np.zeros(3 - nb_refinementByDirection.size, dtype=np.intp))
            nb_refinementByDirection_byPatch = np.tile(
                np.vstack(nb_refinementByDirection),
                (1, self._nb_patch))
        elif np.size(nb_refinementByDirection, 1) < self._nb_patch:
            nb_refinementByDirection_byPatch = np.append(
                nb_refinementByDirection,
                np.zeros((3,
                          self._nb_patch - nb_refinementByDirection.shape[1]),
                         dtype=np.intp), axis=1)
        else:
            nb_refinementByDirection_byPatch = nb_refinementByDirection

        # - Degree elevation
        if np.size(np.shape(nb_degreeElevationByDirection)) == 1:
            nb_degreeElevationByDirection = np.append(
                nb_degreeElevationByDirection,
                np.zeros(3 - nb_degreeElevationByDirection.size,
                         dtype=np.intp))
            nb_degreeElevationByDirection_byPatch = np.tile(
                np.vstack(nb_degreeElevationByDirection),
                (1, self._nb_patch))
        elif np.size(nb_degreeElevationByDirection, 1) < self._nb_patch:
            nb_degreeElevationByDirection_byPatch = np.append(
                nb_degreeElevationByDirection,
                np.zeros(
                    (3,
                     self._nb_patch - nb_degreeElevationByDirection.shape[1]),
                    dtype=np.intp), axis=1)
        else:
            nb_degreeElevationByDirection_byPatch = \
                nb_degreeElevationByDirection

        # - Additional knot insertions
        if additional_knots is None:
            additional_knots = {"patches": np.array([]),
                                "1": np.array([]),
                                "2": np.array([]),
                                "3": np.array([])}

        # NURBS Refinement
        old_mechSet = self.get_mechanicalSettings()
        old_geoSet = self.get_geometricSettings()
        old_indCPbyPatch = self._get_indCPbyPatch()
        [new_mechSet, new_geoSet, new_indCPbyPatch, transformationMatrices] = \
            iga_refinement(nb_degreeElevationByDirection_byPatch,
                           nb_refinementByDirection_byPatch,
                           old_mechSet, old_geoSet, old_indCPbyPatch,
                           additional_knots)

        self.set_geomectricSettings(new_geoSet)
        self.set_mechanicalSettings(new_mechSet)
        self._set_indCPbyPatch(new_indCPbyPatch)
        self._flatten_data()
        self._compute_vectWeight()
        self._update_dof_info()
        self._updateRefinementMatHistory(transformationMatrices)

        # # # Refine nodal distributions
        # # # !!! Assuming only one patch
        # num_patch = 0
        # new_nodal_distributions = {}
        # for key, distrib in self._nodal_distributions.items():
        #     if self._nb_patch > 1:
        #         raise Exception('Nodal distribution refinement is only' +
        #                         ' avaibale for single patch.')
        #     new_dist = transformationMatrices[num_patch][2] * (
        #         transformationMatrices[num_patch][1] * (
        #             transformationMatrices[num_patch][0] * distrib[0]))

        #     new_nodal_distributions[key] = [new_dist.copy(), distrib[1]]

        # self._nodal_distributions = new_nodal_distributions

        # Refine nodal distributions
        # Nodal distribution is refined from the discretisation state of the initial geometry
        new_nodal_distributions = {}
        for key, distrib in self._nodal_distributions_init.items():
            # Get pressure magnitude at nodes
            init_field = distrib[0]
            # Reshape to a two-dimensionnal array
            init_field = np.vstack(init_field)
            # Compute new field
            new_field = self._updateNodalField(init_field)
            # Get initial type of load (does not change)
            l_name = distrib[1]
            # Fill new dictionnary with refined pressure field
            new_nodal_distributions[key] = [new_field.copy(), l_name]
        # Update internal property
        self._nodal_distributions = new_nodal_distributions

        return transformationMatrices

    def refine(self, nb_refinementByDirection,
               nb_degreeElevationByDirection=np.array([0, 0, 0]),
               additional_knots=None):
        """Refinement by degree elevation and knot insertion.

        For more details on the function's arguments see the function
        ``refine_and_getTransformationMatrices``. The present function does not
        return the transformation matrices for refinement.

        """
        self.refine_and_getTransformationMatrices(
          nb_refinementByDirection, nb_degreeElevationByDirection,
          additional_knots)

    def _initRefinementMatHistory(self):
        """Initialise `_refinementMatHistory` to track refinement steps."""
        vectweight = self._get_vectWeight()
        indCPbyPatch = self._indCPbyPatch
        self._refinementMatHistory = []
        for ipatch in range(self._nb_patch):
            cps = indCPbyPatch[ipatch] - 1
            data = vectweight[cps]
            row = np.arange(cps.size)
            col = cps
            M = sp.coo_matrix((data,
                               (row, col)),
                              shape=(row.size, self._nb_cp)).tocsr()
            self._refinementMatHistory.append(M.copy())


    def _updateRefinementMatHistory(self, transformationMatrices):
        """Update `_refinementMatHistory` according to a refinement step.

        Parameters
        ----------
        transformationMatrices : list of sparse matrices
            The matrices that transform the control points from a coarser to
            a finer mesh. Those matrices are considered to be as given by the
            function `refine_and_getTransformationMatrices`.
        """
        for ipatch in range(self._nb_patch):
            M0 = transformationMatrices[ipatch][0] * \
                self._refinementMatHistory[ipatch]
            M1 = transformationMatrices[ipatch][1] * M0
            M2 = transformationMatrices[ipatch][2]*M1
            self._refinementMatHistory[ipatch] = M2

        return None

    def _updateNodalField(self, fieldinit):
        """Update a nodal field according to the previous refinement steps.

        Parameters
        ----------
        fieldinit : array of floats
            The initial field which discretization level corresponds to the
            initial geometry.

        Returns
        -------
        field : array of floats
            The same field but discretized as the current geometry.
        """
        nf = fieldinit.shape[1]
        field = np.zeros((self._nb_cp, nf))
        indCPbyPatch = self._indCPbyPatch
        for ipatch in range(self._nb_patch):
            fieldpatch = self._refinementMatHistory[ipatch].dot(fieldinit)
            field[indCPbyPatch[ipatch]-1, :] = fieldpatch[:, :]
        vectweight = self._get_vectWeight()
        field = sp.diags(1/vectweight).dot(field)

        return field

    def shapeupdate(self):
        """Apply the shape update.

        The shape update is defined by defined by the attribute
        `_shapeparametrization_def`.

        The design parameters (names and values) are set in the attribute
        `_design_parameters`, and are taken as the inputs.

        The function updates the control point coordinates `_COORDS`.
        """
        if not hasattr(self, '_design_parameters'):
            return None
        if len(self._design_parameters) > 0:
            for key, value in self._design_parameters.items():
                paramname = key.translate({ord('<'): None, ord('>'): None})
                exec("%s=%s" % (paramname, value))
            exec(self._shapeparametrization_def, locals())
            coords_coarse = eval('coords').T
            coords_fine = self._updateNodalField(coords_coarse)
            self._COORDS[:, :] = coords_fine.T

    def writeCOORDS(self, filename):
        """Write the control points coordinates to file, ordered by patch."""
        for num_patch in range(0, self._nb_patch):
            np.savetxt(f'{filename}_{num_patch}.txt',
                       self._COORDS[self._indCPbyPatch[num_patch] - 1],
                       delimiter=', ')


    def writeCOORDS_tot(self, filename):
        """Write all control points coordinates."""
        np.savetxt('{}.txt'.format(filename),
                   self._COORDS,
                   delimiter=', ')


    def generate_vtk4controlMeshVisu(self, filename, num_patch, sol=None, output_path='results'):
        """Generate control mesh visualisation for a given patch number."""
        inputs = self.get_inputs4controlMesh(filename, num_patch, output_path=output_path, sol=sol)
        # print(f'Generate {filename}.vtk ...')

        if sol is None:
            generate_vtk(*inputs)
        else:
            generate_vtk_wsol(*inputs)

    @property
    def coords(self):
        """
        Coordinates of control points
        """
        return self._COORDS

    @coords.setter
    def coords(self, coords):
        assert coords.shape == self._COORDS.shape
        self._COORDS = coords

    @property
    def nb_patch(self):
        """"
        Total number of patchs in the model
        """
        return self._nb_patch

    @nb_patch.setter
    def nb_patch(self, nb_patch):
        raise ValueError("Cannot set number of patchs.")

    @property
    def nb_dof_free(self):
        """
        Number of unblocked degrees of freedom
        """
        return self._nb_dof_free

    @property
    def ind_dof_free(self):
        """
        Indices of unblocked degrees of freedom
        """
        return self._ind_dof_free

    @property
    def nb_dof_tot(self):
        """
        Total number of degrees of freedom
        """
        return self._nb_dof_tot

    @property
    def n_kv(self):
        """
        knot vectors size for all patchs
        """
        return self._Nkv

    @property
    def j_pqr(self):
        """
        Degrees for all patchs
        """
        return self._Jpqr

    @property
    def nb_cp(self):
        """
        TODO : document this property
        """
        return self._nb_cp

    @property
    def dim(self):
        """
        TODO: document this property
        """
        return self._dim

    @property
    def ind_cp_by_patch(self):
        """
        TODO: document this property
        """
        return self._indCPbyPatch

    @property
    def ien(self):
        """
        TODO: document this property
        """
        return self._IEN

    @property
    def elt_type(self):
        """
        TODO: document this property
        """
        return self._ELT_TYPE

    @property
    def num_integration_points(self):
        """
        TODO: document this property
        """
        return self._NBPINT

    @property
    def design_parameters(self):
        """
        TODO Document this property
        """
        return self._design_parameters




