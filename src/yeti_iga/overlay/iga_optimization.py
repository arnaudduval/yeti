"""
Shape optimization af an IGA model
"""

# !/usr/bin/env python
# -*- coding: utf-8 -*-

from yeti_iga.preprocessing.igaparametrization import OPTmodelling


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

