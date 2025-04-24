"""
Shape optimization af an IGA model
"""

# !/usr/bin/env python
# -*- coding: utf-8 -*-

import os

import numpy as np

# pylint: disable=no-name-in-module
from ..preprocessing.igaparametrization import OPTmodelling
from ..postprocessing import postproc as pp
from .. import reconstructionSOL as rsol


class IgaOptimization:
    """
    An object to handle shape opyimization of an IGA model
    """
    def __init__(self, iga_model, nb_var, update_function, refinement,
                 d_update_function=None):
        """
        Parameters
        ----------
        iga_model : IgaModel
            IGA model on which optimization is made
        nb_var : int
            Number of design variables
        update_function : function
            Function updating IGA model, based on initial CP coordinates and
            design variables
        refinement : Refinement
            Refinement from design model to analysis model
        d_update_function : function
            Function comuting the derivative of the shape update function with
            respect to the design variables (default=None)
        """

        self.update_function = update_function
        self.d_update_function = d_update_function
        self.iga_model = iga_model

        def legacy_update_function(coords_0, iga_param, x):
            self.iga_model.iga_param = iga_param
            self.update_function(coords_0.T, self.iga_model, x)

        if d_update_function is not None:
            def legacy_d_update_function(coords_0, iga_param, x, ivar):
                self.iga_model.iga_param = iga_param
                return self.d_update_function(coords_0.T, self.iga_model, x, ivar).T
        else:
            legacy_d_update_function = None

        self._opt_pb = OPTmodelling(self.iga_model.iga_param,
                                    nb_var,
                                    legacy_update_function,
                                    nb_degreeElevationByDirection=refinement.degrees_legacy,
                                    nb_refinementByDirection=refinement.subdivision_legacy,
                                    fct_dervShapeParam=legacy_d_update_function
                                    )
    # Response function and thier gradients :
    #  - volume
    #  - compliance
    #  - displacement


    def volume(self, x, listpatch=None):
        """
        Compute volume for a given set of design variables

        Parameters
        ----------
        x : np.array(dtype=float)
            Array containing the design variables
        listpatch : numpy.array(dtype=int)
            An array of 0 ro 1 indocating which patch must be taken into
            account for volume computation
            Default = None

        Returns
        -------
        volume : float
            Volume of the model
        """

        if listpatch:
            if listpatch.size != self.iga_model.nb_patch:
                raise ValueError("Length of listpatch parameter must be equal " \
                "to {self.iga_model.nb_patch}")

            if not np.all(np.isin(listpatch, [0, 1])):
                raise ValueError("listpatch parameters must exclusively " \
                "contain 0 or 1 values.")

        return self._opt_pb.compute_volume(x, listpatch=listpatch)


    def grad_volume_analytic(self, x, listpatch=None):
        """
        Compute gradient of the volume with respect to a set of deign variables
        using analytic method

        Parameters
        ----------
        x : numpy.array(dtype=float)
            Array containing the design variables
        listpatch : numpy.array(dtype=int)
            An array of 0 ro 1 indicating which patch must be taken into
            account for volume computation
            Default = None

        Returns
        -------
        grad_volume : numpy.array(dtype=float)
            gradient of the volume
        """
        if listpatch:
            if listpatch.size != self.iga_model.nb_patch:
                raise ValueError("Length of listpatch parameter must be equal " \
                "to {self.iga_model.nb_patch}")

            if not np.all(np.isin(listpatch, [0, 1])):
                raise ValueError("listpatch parameters must exclusively " \
                "contain 0 or 1 values.")

        return self._opt_pb.compute_gradVolume_AN(x, listpatch)


    def grad_volume_finite_differences(self, x, eps=1.e-6, centered=False, listpatch=None):
        """
        Compute gradient of the volume with respect to a set of design variables
        using finite differences

        Parameters
        ----------
        x : numpy.array(dtype=float)
            Array containing the design variables
        eps : float
            epsilon value for finite difference computation
            Default = 1.e-6
        centered : bool
            Boolean indicating if centered finite differences must be used
            Default = False
        listpatch : numpy.array(dtype=int)
            An array of 0 ro 1 indicating which patch must be taken into
            account for volume computation
            Default = None

        Returns
        -------
        grad_volume : numpy.array(dtype=float)
            gradient of the volume
        """

        return self._opt_pb.compute_gradVolume_DF(x, eps=eps,
                                                  centerFD=centered,
                                                  listpatch=listpatch)


    def compliance(self, x):
        """
        Compute compliance for a given set of design variables

        Parameters
        ----------
        x : np.array(dtype=float)
            Array containing the design variables

        Returns
        -------
        compliance : float
            Discrete compliance of the model
        """
        return self._opt_pb.compute_compliance_discrete(x)


    def grad_compliance_analytic(self, x):
        """
        Compute gradient of the compliance with respect to a set of design
        variables using analytic method

        Parameters
        ----------
        x : numpy.array(dtype=float)
            Array containing the design variables

        Returns
        -------
        grad_compliance : numpy.array(dtype=float)
            Gradient of the compliance
        """

        return self._opt_pb.compute_gradCompliance_AN(x)


    def grad_compliance_finite_differences(self, x, eps=1.e-6, centered=False):
        """
        Compute gradient of the compliance with respect to a set of design
        variables using finite differences

        Parameters
        ----------
        x : numpy.array(dtype=float)
            Array containing the design variables
        eps : float
            epsilon value for finite difference computation
            Default = 1.e-6
        centered : bool
            Boolean indicating if centered finite differences must be used
            Default = False

        Returns
        -------
        grad_compliance : numpy.array(dtype=float)
            Gradient of the compliance
        """

        return self._opt_pb.compute_gradCompliance_FD(x, eps=eps, centerFD=centered)


    def displacement(self, x, xi, ipatch=0):
        """
        Compute displacement at a point of a patch for a given set of design
        variables

        Parameters
        ----------
        x : np.array(dtype=float)
            Array continaing the design variables
        xi : np.array(dtype=float)
            Array containing the paralmetric coordinates of the point
        ipatch : int
            Index of the considered patch (default=0)

        Returns
        -------
        displacement : np.array(dtype=float)
            Displacement at point
        """

        return self._opt_pb.compute_displacement(x, xi, ipatch+1)

    def grad_displacement_analytic(self, x, xi, ipatch=0):
        """
        Compute gradient of the displacement at a points with respect to a set
        of design variables using analytic method

        Parameters
        ----------
        x : np.array(dtype=float)
            Array continaing the design variables
        xi : np.array(dtype=float)
            Array containing the paralmetric coordinates of the point
        ipatch : int
            Index of the considered patch (default = 0)
            Not taken into account, computation is for patch 0

        Returns
        -------
        grad_displacement : np.array(dtype=float)
            Gradient of the displacement
        """

        return self._opt_pb.compute_gradDisplacement_AN(x, xi)

    def grad_displacement_finite_differences(self, x, xi, ipatch=0,
                                             eps=1.e-6, centered=False):
        """
        Compute gradient of the displacement at a points with respect to a set
        of design variables using finite differences

        Parameters
        ----------
        x : np.array(dtype=float)
            Array continaing the design variables
        xi : np.array(dtype=float)
            Array containing the paralmetric coordinates of the point
        ipatch : int
            Index of the considered patch (default=0)
            Not taken into account, computation is for patch 0
        eps : float
            epsilon value for finite difference computation
            Default = 1.e-6
        centered : bool
            Boolean indicating if centered finite differences must be used
            Default = False

        Returns
        -------
        grad_displacement : np.array(dtype=float)
            Gradient of the displacement
        """

        return self._opt_pb.compute_gradDisplacement_FD(x, xi, eps=eps, centerFD=centered)



    def write_analysis_solution_vtu(self, filename,
                                    x=None,
                                    per_patch=False,
                                    refinement=np.array([3, 3, 3]),
                                    data_flag=np.array([True, False, False])):
        """
        Write analysis result of the model for a set of design variables x
        TODO: make BSpline output when possible

        Parameters
        ----------
        filename : string
            Path to the file to write
        x : numpy.array(dtype=float)
            Array containing the design variables
            default=None
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

        output_path = os.path.dirname(os.path.realpath(filename))
        filename = os.path.splitext(os.path.basename(filename))[0]

        sol, _ = rsol.reconstruction(
            **self._opt_pb.fine_parametrization.get_inputs4solution(
                self._opt_pb.save_sol_fine))
        pp.generatevtu(*self._opt_pb.fine_parametrization.get_inputs4postprocVTU(
            filename,  sol.transpose(),
            nb_ref=refinement,
            Flag=data_flag,
            output_path=output_path))


    def write_design_model_control_mesh(self, filename,
                                        ipatch=0,
                                        x=None):
        """
        Write control mesh of the design model for a set of design
        variables x

        Parameters
        ----------
        filename : string
            Path to the file to write
        ipatch : int
            Index of patch top write, default=0
        x : numpy.array(dtype=float)
            Array containing the design variables
            default=None
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


        self._opt_pb.coarse_parametrization.generate_vtk4controlMeshVisu(
            filename, ipatch, output_path=output_path)
