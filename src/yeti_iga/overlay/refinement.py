# !/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Module handling refinement properties for an IGA model
"""

import numpy as np

class Refinement:
    """
    An object defining an IGA refinement (degree elevation and subdivision) on
    a set of patchs
    """
    def __init__(self, nb_patch):
        """
        Parameters
        ----------
        nb_patch : int
            Number of patch on which refinement will operate
        """

        self.perpatch = []
        for _ in range(nb_patch):
            self.perpatch.append({'degree_elevation': np.array([0, 0, 0], dtype=int),
                                  'subdivision': np.array([0, 0, 0], dtype=int)
                                  })

    def set_refinement(self, ipatch, degree_elevation, subdivision):
        """
        Set refinement (degree levation and subdivision for a given patch)

        Parameters
        ----------
        ipatch : int
            Index of patch to refine
        degree_elevation : np.array(dtype=int)
            Number of degree elevation for each parametric direction
        subdivision : np.array(dtype=int)
            Number subdivision for each parametric direction
        """

        if ipatch > len(self.perpatch):
            raise ValueError(f"Maximum index for patch refinementb is {len(self.perpatch)}")

        self.perpatch [ipatch]['degree_elevation'][:degree_elevation.size] = degree_elevation[:]
        self.perpatch [ipatch]['subdivision'][:subdivision.size] = subdivision[:]

    @property
    def degrees_legacy(self):
        """
        Degree elevation for each patch of an IGA paramatrization in legacy format
        """
        degrees = np.zeros((3, len(self.perpatch)),dtype=np.intp)
        for i, _ in enumerate(self.perpatch):
            degrees[:, i] = self.perpatch[i]['degree_elevation']

        return degrees


    @property
    def subdivision_legacy(self):
        """
        Subdivision number for each patch of an IGA paramatrization in legacy format"
        """
        subdivisions = np.zeros((3, len(self.perpatch)),dtype=np.intp)
        for i, _ in enumerate(self.perpatch):
            subdivisions[:, i] = self.perpatch[i]['subdivision']

        return subdivisions
