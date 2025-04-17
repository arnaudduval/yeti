# !/usr/bin/env python
# -*- coding: utf-8 -*-

class IgaOptimization:
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

