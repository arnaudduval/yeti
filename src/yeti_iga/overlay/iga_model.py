# !/usr/bin/env python
# -*- coding: utf-8 -*-

"""
A module handling an IGA model for structural analysis
"""

import os

import numpy as np
import scipy.sparse as sp

from ..preprocessing.igaparametrization.IGA_parametrization import IGAparametrization
from ..stiffmtrx_elemstorage import sys_linmat_lindef_static
from .. import reconstructionSOL as rsol
from ..postprocessing import postproc as pp
from ..preprocessing.igaparametrization import IGAmanip as manip

from .material import ElasticMaterial

class IgaModel:
    def __init__(self, model_type):
        """
        Parameters
        ----------
        model_type : string
            Model type : '2D solid', '3D solid' or '3D shell'
        """
        if model_type not in ['2D solid', '3D solid', '3D shell']:
            raise ValueError(f'{model_type} is not a valid model type')
        self.model_type = model_type

        self.patchs = []
        self.local_global_cp = []       # Maps local CP numbering to global one
        self.local_global_el = []       # Maps local element numbering to global one
        self.boundary_conditions = []

        self._refined_flag = False           # Flag set to true when a refinement has been made
                                             # To prevent new patch or BC addition

        self.iga_param = IGAparametrization()

        # Global
        if model_type in('3D solid', '3D shell'):
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
        # TODO ERROR : props and jprops applies on elements, not patch
        self.iga_param._PROPS = []
        self.iga_param._JPROPS = np.array([], dtype=int)

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


    def add_patch(self, patch):
        """
        Add a patch to model

        Parameters
        ----------
        patch : Patch
            patch to add

        Return
        ------
        index : int
            index of added patch
        """

        if self._refined_flag:
            raise Exception("Can not add a new patch if model has already been refined")

        self.patchs.append(patch)

        self.iga_param._nb_patch += 1
        self.iga_param._nb_elem += patch.connectivity.shape[0]
        self.iga_param._nb_cp += patch.control_points.shape[0]

        self.iga_param._ELT_TYPE = np.append(
            self.iga_param._ELT_TYPE, patch.element_type
        )
        if patch.element_type == 'U1':
            self.iga_param._TENSOR = np.append(
                self.iga_param._TENSOR, 'THREED'
            )
        elif patch.element_type == 'U3':
            self.iga_param._TENSOR = np.append(
                self.iga_param._TENSOR, 'PSTRESS'
            )

        self.iga_param._dim = np.append(self.iga_param._dim, len(patch.knot_vectors))

        self.iga_param._NBPINT = np.append(
            self.iga_param._NBPINT, np.prod(patch.degrees + 1)
        )
        if self.iga_param._dim == 3:
            self.iga_param._Jpqr = np.append(
                self.iga_param._Jpqr, np.array([patch.degrees]).T, axis=1
                )
            # TODO Verify if a [] is needed, like in the 2 dim case
            self.iga_param._Nijk = np.append(
                self.iga_param._Nijk, (patch.spans+1).T, axis=1
                )
        elif self.iga_param._dim == 2:
            self.iga_param._Jpqr = np.append(
                self.iga_param._Jpqr, np.array([np.append(patch.degrees, 0)]).T, axis=1
                )
            self.iga_param._Nijk = np.append(
                self.iga_param._Nijk, np.array([np.append(patch.spans+1, 0)]).T, axis=1
                )
        else:
            raise ValueError("Patch dimension must be 2 or 3")



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
        # Add other properties
        self.iga_param._PROPS[-1 ] = np.append(self.iga_param._PROPS[-1 ],
                                               patch.properties)

        self.iga_param._JPROPS = np.append(
            self.iga_param._JPROPS, patch.properties.size + 1
            )

        self.iga_param._flatten_data()
        self.iga_param._indCPbyPatch = self.iga_param._autoset_indCPbyPatch()
        self.iga_param._compute_vectWeight()
        self.iga_param._update_dof_info()
        self.iga_param._initRefinementMatHistory()

        return self.nb_patch - 1

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

        if self._refined_flag:
            raise Exception("Can not add a new boundary condition if model " \
            "has already been refined")

        self.boundary_conditions.append((ipatch, bc))



        for dof in bc.dof:
            self.iga_param._nb_bc += 1
            self.iga_param._bc_target.append(self.local_global_cp[ipatch][bc.cp_index]+1)

            # Warning : not clear if is a number of CP or not ???
            # See /src/yeti_iga/preprocessing/mechanicalmodel/model.py
            self.iga_param._bc_target_nbelem = np.append(
                self.iga_param._bc_target_nbelem, bc.cp_index.size
                )
            self.iga_param._bc_values = np.append(
                self.iga_param._bc_values, np.array([[dof+1],[bc.value]]), axis=1
                )

        self.iga_param._flatten_data()
        self.iga_param._update_dof_info()


    def add_distributed_load(self, ipatch, dload):
        """
        Add a distributed load

        Parameters
        ----------
        ipatch: int
            Patch index on which distributed load is applied
        dload : DistributedLoad
            Distributed load to add
        """

        if self._refined_flag:
            raise Exception("Can not add a new load if model has already been refined")

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

        # TODO  Test with a centrifugal force case
        self.iga_param._additionalLoadInfos.append(np.array([], dtype=float))
        self.iga_param._nb_additionalLoadInfos = np.append(
            self.iga_param._nb_additionalLoadInfos, 0
        )

        self.iga_param._flatten_data()

    def refine_patch(self, ipatch,
                     nb_degree_elevation=np.zeros(3),
                     nb_subdivision=np.zeros(3),
                     additional_knots=None):
        """
        Refine a patch of the model with k-refinement:
         1 - degree elevation
         2 - knot insertion
         3 - subdivision


        Parameters
        ----------
        ipatch: int
            Index of patch to refine
        nb_degree_elevation : numpy.array(shape=(3,), dtype=numpy.intp)
            Number of degree elevation for each parametric direction
            default = numpy.zeros(3)
        nb_subdivision : numpy.array(shape=(3,), dtype=numpy.intp)
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

        dim = self.iga_param._dim[ipatch]

        nb_deg[:dim, ipatch] = nb_degree_elevation[:dim]
        nb_ref[:dim, ipatch] = nb_subdivision[:dim]

        if additional_knots is not None:
            additional_knots_legacy = {"patches": np.array([ipatch]),
                                    "1": additional_knots[0],
                                    "2": additional_knots[1],
                                    "3": additional_knots[2]
                                    }
            self.iga_param.refine(nb_ref, nb_deg, additional_knots=additional_knots_legacy)
        else:
            self.iga_param.refine(nb_ref, nb_deg)

        self._refined_flag = True



    def build_stiffness_matrix(self):
        """
        Build stiffness matrix and right hand member

        Returns
        -------
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

    @property
    def nb_patch(self):
        """
        Return the number of patchs of the model
        """

        return len(self.patchs)

    @property
    def cp_indices(self):
        """
        Returns an array containing all control points indices of the model

        Returns
        -------
        cp_indices : numpy.array(dtype=int)
            All control points indices
        """
        return np.arange(self.iga_param.nb_cp)

    @property
    def cp_coordinates(self):
        """
        Returns all control points coordinates of the model

        Returns
        ------
        cp_coordinates : numpy.array(shape=(., 3), dtype=float)
            Coordinates of control points
        """

        return self.iga_param._COORDS.T

    @cp_coordinates.setter
    def cp_coordinates(self, coords):
        """
        Set the control points coordinates of the model

        Parameters
        ----------
        coords : numpy.array(shape=(., 3), dtype=float)
            new coordinates
        """

        if coords.T.shape != self.iga_param._COORDS:
            raise ValueError("coords parameters must have {self.iga_param._COORDS.T.shape} shape")


        self.iga_param._COORDS = coords.T







    def write_solution_vtu(self, x, filename,
                           per_patch=False,
                           refinement=np.array([3, 3, 3]),
                           data_flag=np.array([True, False, False])):
        """
        Write the solution of an IGA analysis to a VTU file
        TODO: make BSpline output when possible

        Parameters
        ----------
        x : numpy.array(dtype=float)
            Array containing the result of an IGA analysis
        filename : string
            File name
        per_patch : bool (default=False)
            If set to true, a separate file is generated for each patch.
            Non-embedded patchs are generated using Bezier extraction feature
            of VTK format. Produced files will have a _pxxx suffix.
        refinement : numpy.array(shape=(3,), dtype=int) (default=[3, 3, 3])
            Refinement applied when files are egenerated using a classical FE
            data structure
        data_flag : numpy.array(shape=(3,), dtype=bool)
            (default=[True, False, False])
            Boolean array indicating generated fields : [0] : displacement,
            [1] : stress, [2] : strain

        """

        if per_patch:
            raise NotImplementedError(
                "Per patch feature not yet implemented"
                )

        if len(os.path.splitext(os.path.basename(filename))[1]) == 0:
            raise ValueError(
                "File name has no extension"
            )

        if os.path.splitext(os.path.basename(filename))[-1].lower() != '.vtu':
            raise ValueError(
                "File name must have .vtu extension"
            )

        sol, _ = rsol.reconstruction(**self.iga_param.get_inputs4solution(x))

        output_path = os.path.dirname(os.path.realpath(filename))

        filename = os.path.splitext(os.path.basename(filename))[0]

        pp.generatevtu(*self.iga_param.get_inputs4postprocVTU(
            filename, sol.T, nb_ref=refinement,
            Flag=data_flag,
            output_path=output_path))

    def write_control_mesh_vtk(self, filename, ipatch=0):
        """
        Write control mesh of the model in a VTU file
        A separate file is generated for each patch

        Parameters
        ----------
        filename : string
            File name
        ipatch : int
            Index of patch to write, default=0
        """

        if len(os.path.splitext(os.path.basename(filename))[1]) == 0:
            raise ValueError(
                "File name has no extension"
            )

        if os.path.splitext(os.path.basename(filename))[-1].lower() != '.vtk':
            raise ValueError(
                "File name must have .vtk extension"
            )

        output_path = os.path.dirname(os.path.realpath(filename))
        filename = os.path.splitext(os.path.basename(filename))[0]


        self.iga_param.generate_vtk4controlMeshVisu(filename, ipatch, output_path=output_path)


    def get_corners_cp_indices(self, ipatch=0):
        """
        Return global indices of control points corresponding to a given patch
        corners

        Parameters
        ----------
        ipatch : int (default = 0)
            index of patch

        Returns
        -------
        indices : numpy.array
            indices of corners control points
        """

        indices = manip.get_vertexCPindice(self.iga_param._Nkv,
                                           self.iga_param._Jpqr,
                                           self.iga_param._dim,
                                           num_patch=ipatch)
        return indices

    def get_face_cp_indices(self, id_face, ipatch=0):
        """
        Return global indices of control points corresponding to a face of a
        given patch

        Parameters
        ----------
        id_face : int
            identifier of face : 1, 2, 3, 4, 5 or 6
        ipatch : int (default = 0)
            index of patch

        Returns
        -------
        indices : numpy.array(dtype=int)
            indices of face control points
        """

        if id_face < 1 or id_face > 6:
            raise ValueError('id_face parameter must be in [1 6]')

        indices = manip.get_boundCPindice(self.iga_param.n_kv,
                                          self.iga_param.j_pqr,
                                          id_face,
                                          num_patch=ipatch)

        return indices
