# !/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Module handling distributed load on an IGA model
"""

class DistributedLoad:
    """
    Object handling distributed load on an IGA model
    A distributed load applies on a set of elements, has a type and a magnitude
    """
    def __init__(self, el_index, dl_type, magnitude):
        """
        Parameters
        ----------
        el_index : np.array([], dtype=int)
            local indices of elements on which distributed load is applied
        dl_type : int
            Type of applied dload (U10, U11, etc.)
        magnitude : float
            Magnitude of applied distributed load

        """

        self.el_index = el_index
        self.dl_type = dl_type
        self.magnitude = magnitude
