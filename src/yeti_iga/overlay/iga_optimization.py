"""
Shape optimization af an IGA model
"""

# !/usr/bin/env python
# -*- coding: utf-8 -*-

import os

import numpy as np

from ..preprocessing.igaparametrization import OPTmodelling
from ..postprocessing import postproc as pp
from .. import reconstructionSOL as rsol

class IgaOptimization:
    """
    An object to handle shape opyimization of an IGA model
    """
    def __init__(self, iga_model, nb_var, update_function, refinement):
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
        """
        self._opt_pb = OPTmodelling(iga_model.iga_param,
                                    nb_var,
                                    update_function,
                                    nb_degreeElevationByDirection=refinement.degrees_legacy,
                                    nb_refinementByDirection=refinement.subdivision_legacy
                                    )

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

    def volume(self, x):
        """
        Compute volume for a given set of design variables

        Parameters
        ----------
        x : np.array(dtype=float)
            Array containing the design variables

        Returns
        -------
        volume : float
            Volume of the model
        """
        # TODO : add mask for patchs taken into account

        return self._opt_pb.compute_volume(x)

    def grad_volume_analytic(self, x):
        """
        Compute gradient of the volume with respect to a set of deign variables
        using analytic method

        Parameters
        ----------
        x : numpy.array(dtype=float)
            Array containing the design variables

        Returns
        -------
        grad_volume : numpy.array(dtype=float)
            gradient of the volume
        """
        # TODO : add mask for patchs taken into account

        return self._opt_pb.compute_gradVolume_AN(x)

    def grad_volume_finite_differences(self, x):
        """
        Compute gradient of the volume with respect to a set of design variables
        using finite differences

        Parameters
        ----------
        x : numpy.array(dtype=float)
            Array containing the design variables

        Returns
        -------
        grad_volume : numpy.array(dtype=float)
            gradient of the volume
        """
        # TODO : add mask for patchs taken into account
        # TODO add specific parameters for finite differences

        return self._opt_pb.compute_gradVolume_DF(x)

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

    def grad_compliance_finite_differences(self, x):
        """
        Compute gradient of the compliance with respect to a set of design
        variables using finite differences

        Parameters
        ----------
        x : numpy.array(dtype=float)
            Array containing the design variables

        Returns
        -------
        grad_compliance : numpy.array(dtype=float)
            Gradient of the compliance
        """
        # TODO add specific parameters for finite differences

        return self._opt_pb.compute_gradCompliance_FD(x)

    def write_analysis_solution_vtu(self, x, filename,
                                    per_patch=False,
                                    refinement=np.array([3, 3, 3]),
                                    data_flag=np.array([True, False, False])):
        """
        Write analysis result of the model for a set of design variables x
        TODO: make BSpline output when possible

        Parameters
        ----------
        x : numpy.array(dtype=float)
            Array containing the design variables
        filename : string
            Path to the file to write
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

        sol, _ = rsol.reconstruction(**self._opt_pb.fine_parametrization.get_inputs4solution(self._opt_pb.save_sol_fine))
        pp.generatevtu(*self._opt_pb.fine_parametrization.get_inputs4postprocVTU(
            filename,  sol.transpose(),
            nb_ref=refinement,
            Flag=data_flag,
            output_path=output_path))



