# !/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import scipy.sparse as sp

from ..preprocessing.igaparametrization.IGA_parametrization import IGAparametrization
from ..stiffmtrx_elemstorage import sys_linmat_lindef_static

from .material import ElasticMaterial

class IgaModel:
    def __init__(self, model_type):
        if model_type not in ['2D solid', '3D solid', '3D shell']:
            raise ValueError(f'{model_type} is not a valid model type')
        self.model_type = model_type

        self.patchs = []
        self.local_global_cp = []
        self.local_global_el = []
        self.boundary_conditions = []

        self.iga_param = IGAparametrization()

        # Global
        if model_type == '3D solid':
            # TODO add an assertaion for de dimension of coordinates array
            self.iga_param._mcrd = 3
        self.iga_param._nb_patch = 0
        self.iga_param._nb_elem = 0
        self.iga_param._nb_cp = 0


        # Patchs
        self.iga_param._ELT_TYPE = np.array([])
        self.iga_param._TENSOR = np.array([])
        self.iga_param._NBPINT = np.array([], dtype=int)
        self.iga_param._Jpqr = np.array([[],[],[]], dtype=int)
        self.iga_param._Ukv = []
        self.iga_param._Nkv = np.array([[],[],[]], dtype=int)
        self.iga_param._dim = np.array([], dtype=int)
        self.iga_param._COORDS = np.array([[],[],[]], dtype=float)
        self.iga_param._IEN = []
        self.iga_param._Nijk = np.array([[],[],[]], dtype=int)
        self.iga_param._weight = []
        self.iga_param._elementsByPatch = np.array([], dtype=int)
        self.iga_param._MATERIAL_PROPERTIES = np.array([[],[],[]], dtype=float)
        self.iga_param._N_MATERIAL_PROPERTIES = np.array([], dtype=int)
        self.iga_param._elementsByPatch = np.array([], dtype=int)
        self.iga_param._nnode = np.array([], dtype=int)
        # ERROR : props and jprops applies on elements, not patch
        self.iga_param._PROPS = []
        self.iga_param._JPROPS = np.array([], dtype=int)
        # Pas nécessaire : il existe un autoset
        # self.iga_param._indCPbyPatch = np.array([], dtype=object)

        # Boundary conditions
        self.iga_param._nb_bc = 0
        self.iga_param._bc_target = []
        self.iga_param._bc_values = np.array([[],[]])
        self.iga_param._bc_target_nbelem = np.array([], dtype=int)

        # Loads
        self.iga_param._nb_load = 0
        self.iga_param._indDLoad = np.array([[]], dtype=int)
        self.iga_param._JDLType = np.array([], dtype=int)
        self.iga_param._ADLMAG = np.array([], dtype=float)
        self.iga_param._load_target_nbelem = np.array([], dtype=int)


        self.iga_param._additionalLoadInfos = []
        self.iga_param._nb_additionalLoadInfos = np.array([], dtype=int)

        # Nodal distributions
        self.iga_param._nodal_distributions = {}
        self.iga_param._nodal_distributions_init = {}




        return

    def add_patch(self, patch):
        """
        Add a patch to model

        Parameters
        ----------
        patch : Patch
            patch to add
        """

        self.patchs.append(patch)

        self.iga_param._nb_patch += 1
        self.iga_param._nb_elem += patch.connectivity.shape[0]
        self.iga_param._nb_cp += patch.control_points.shape[0]

        self.iga_param._ELT_TYPE = np.append(
            self.iga_param._ELT_TYPE, patch.element_type
        )
        if patch.element_type == 'U1':
            self.iga_param._TENSOR = np.append(
                self.iga_param._TENSOR, 'THREED',
            )

        self.iga_param._NBPINT = np.append(
            self.iga_param._NBPINT, np.prod(patch.degrees + 1)
        )
        self.iga_param._Jpqr = np.append(
            self.iga_param._Jpqr, np.array([patch.degrees]).T, axis=1
        )
        self.iga_param._Ukv.append(patch.knot_vectors)
        self.iga_param._COORDS = np.append(
            self.iga_param._COORDS, patch.control_points.T, axis=1
        )

        # Set global cp numbering
        if len(self.patchs) == 1:
            max_cp_idx = -1
        else:
            max_cp_idx = np.max([arr.max() for arr in self.local_global_cp])
        self.local_global_cp.append(np.arange(patch.control_points.shape[0])+max_cp_idx+1)

        Nkv = np.array([0, 0, 0], dtype=int)
        for i, kv in enumerate(patch.knot_vectors):
            Nkv[i] = kv.size
        self.iga_param._Nkv = np.append(
            self.iga_param._Nkv, np.array([Nkv]).T, axis=1
        )

        self.iga_param._dim = np.append(self.iga_param._dim, len(patch.knot_vectors))

        self.iga_param._IEN.append(patch.connectivity[:, self.local_global_cp[-1]]+1)
        # TODO Possible de faire plus simple sans ajouter de None.
        # On veut en sortie un array de array,
        # ex : '_indCPbyPatch': array([array([1, 2, 3, 4, 5, 6, 7, 8])], dtype=object)
        # self.iga_param._indCPbyPatch = np.append(
        #     self.iga_param._indCPbyPatch, None)
        # self.iga_param._indCPbyPatch[-1] = self.local_global_cp[-1]+1
        # ==> Fait par une routine déja existante dans IGA_parametrization

        # Set global elements numbering
        if len(self.patchs) == 1:
            max_el_index = -1
        else:
            max_el_index = np.max([arr.max() for arr in self.local_global_el])

        self.local_global_el.append(np.arange(patch.connectivity.shape[0])+max_el_index+1)

        self.iga_param._elementsByPatch = np.append(
            self.iga_param._elementsByPatch, patch.connectivity.shape[0]
        )
        self.iga_param._nnode = np.append(
            self.iga_param._nnode, patch.connectivity.shape[1]
        )
        self.iga_param._Nijk = np.append(
            self.iga_param._Nijk, (patch.spans+1).T, axis=1
        )
        # Reorder weights to assign them per element
        for elt in patch.connectivity:
            self.iga_param._weight.append(patch.weights[elt])

        if isinstance(patch.material, ElasticMaterial):
            self.iga_param._MATERIAL_PROPERTIES = np.append(
                self.iga_param._MATERIAL_PROPERTIES,
                np.array([
                    [patch.material.young_modulus],
                    [patch.material.poisson_ratio],
                    [patch.material.density]
                ]),
                axis=1
            )
            self.iga_param._N_MATERIAL_PROPERTIES = np.append(
                self.iga_param._N_MATERIAL_PROPERTIES, 3
            )
        else:
            raise TypeError("Only ElasticMaterial type is supported")

        # Patch properties : first property is the index of patch, starting at 1
        self.iga_param._PROPS.append(np.array([float(len(self.patchs))], dtype=object))
        # WARNING: JPROPS is not a list of arrays ???
        self.iga_param._JPROPS = np.append(
            self.iga_param._JPROPS, len(self.patchs)
        )

        self.iga_param._flatten_data()
        self.iga_param._indCPbyPatch = self.iga_param._autoset_indCPbyPatch()
        self.iga_param._compute_vectWeight()
        self.iga_param._update_dof_info()
        self.iga_param._initRefinementMatHistory()

    def add_boundary_condition(self, ipatch, bc):
        """
        Add a boundary condition to the model

        Parameters
        ----------
        ipatch : int
            Patch index on which boundary condition is applied
        bc : BoundaryCondition
            Boundary condition to add
        """

        self.boundary_conditions.append((ipatch, bc))

        self.iga_param._nb_bc += 1

        self.iga_param._bc_target.append(self.local_global_cp[ipatch][bc.cp_index]+1)
        # Warning : not clear if is a number of CP or not ???
        self.iga_param._bc_target_nbelem = np.append(
            self.iga_param._bc_target_nbelem, bc.cp_index.size
        )
        self.iga_param._bc_values = np.append(
            self.iga_param._bc_values, np.array([[bc.dof+1],[bc.value]]), axis=1
        )

        self.iga_param._flatten_data()
        self.iga_param._update_dof_info()


    def add_distributed_load(self, ipatch, dload):
        """
        Add a distributed load

        Paramatees
        ----------
        ipatch: int
            Patch index on which distributed load is applied
        dload : DistributedLoad
            Distributed load to add
        """

        self.iga_param._nb_load += 1

        self.iga_param._indDLoad = np.append(
            self.iga_param._indDLoad, [self.local_global_el[ipatch][dload.el_index]+1], axis=1
        )
        self.iga_param._JDLType = np.append(
            self.iga_param._JDLType, dload.dl_type
        )
        self.iga_param._ADLMAG = np.append(
            self.iga_param._ADLMAG, dload.magnitude
        )
        self.iga_param._load_target_nbelem = np.append(
            self.iga_param._load_target_nbelem, dload.el_index.size
        )

        # TODO faire quelque chose avec
        # _additionalLoadInfos
        # _nb_additionalLoadInfos
        # ==> A tester avec un cas de chargement centrifuge
        self.iga_param._additionalLoadInfos.append(np.array([], dtype=float))
        self.iga_param._nb_additionalLoadInfos = np.append(
            self.iga_param._nb_additionalLoadInfos, 0
        )

        self.iga_param._flatten_data()

    def refine_patch(self, ipatch,
                     nb_degree_elevation=np.zeros(3),
                     nb_refinements=np.zeros(3),
                     additional_knots=None):
        """
        Refine a patch of the model with k-refinement:
         - degree elevation
         - subdivision
         - knot insertion

        Parameters
        ----------
        ipatch: int
            Index of patch to refine
        nb_degree_elevation : numpy.array(shape=(3,), dtype=numpy.intp)
            Number of degree elevation for each parametric direction
            default = numpy.zeros(3)
        nb_refinement : numpy.array(shape=(3,), dtype=numpy.intp)
            Number of knot vector subdivision for each parametric direction
            default = numpy.zeros(3)
        additional_knots : list of np.array(dtype=float)
            knots to insert for each parametric direction
            default = None
        """

        if ipatch >= len(self.patchs):
            raise ValueError(f'id_patch must be <= {len(self.patchs) - 1}')

        nb_deg = np.zeros((3, self.iga_param.nb_patch),dtype=np.intp)
        nb_ref = np.zeros((3, self.iga_param.nb_patch),dtype=np.intp)
        nb_deg[:, ipatch] = nb_degree_elevation
        nb_ref[:, ipatch] = nb_refinements

        if additional_knots is not None:
            additional_knots_legacy = {"patches": np.array([ipatch]),
                                    "1": additional_knots[0],
                                    "2": additional_knots[1],
                                    "3": additional_knots[2]
                                    }
            self.iga_param.refine(nb_ref, nb_deg, additional_knots=additional_knots_legacy)
        else:
            self.iga_param.refine(nb_ref, nb_deg)


    def build_stiffness_matrix(self):
        """
        Build stiffness matrix and right hand member

        Retuns
        ------
        stiff : scipy.sparse.csc_matrix
            Assembled stifness matrix
        rhs : numpy.array
            Rught hans side vector
        """

        params = self.iga_param.get_inputs4system_elemStorage()

        data, row, col, rhs = sys_linmat_lindef_static(*params)

        stiff_side = sp.coo_matrix((data, (row, col)),
                               shape=(self.iga_param.nb_dof_tot,
                                      self.iga_param.nb_dof_tot),
                               dtype='float64').tocsc()
        return stiff_side + stiff_side.transpose(), rhs

    @property
    def idof_free(self):
        """
        Return indices of free degrees of freedom
        """

        return self.iga_param.ind_dof_free[:self.iga_param.nb_dof_free]-1
