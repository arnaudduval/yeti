# !/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Module handling a NURBS or BSpline patch
"""

import numpy as np

class Patch:
    """
    Object handling a patch
    """
    def __init__(self, element_type, degrees, knot_vectors, control_points,
                 weights, connectivity, spans, material,
                 properties=np.array([], dtype=float)):
        """
        Parameters
        ----------
        element_type : string
            type of element
        degrees : numpy.array
            degree in each direction of parametric space
        knot_vectors : list of numpy.array
            knot vectors for each direction of parametric space
        control_points : numpy.array
            control points coordinates : control_points[.,.]
        weights : numpy.array
            control points weights
        connectivity : numpy.array
            connectivity matrix, indexing starts at 0, numbering is local to patch
        spans : numpy.array
            spans values for each element, numbering starts at 0
        material : Material
            material properties
        properties : numpy.array(dtype=float)
            extra patch properties (thickness, hull index, ...)
        """

        if element_type == 'U3' and properties.size == 0:
            raise ValueError('Shell patch should have at least 1 property defining its thuckness')


        self.element_type = element_type
        self.degrees = degrees
        self.knot_vectors = knot_vectors
        self.control_points = control_points
        assert weights.shape[0] == control_points.shape[0]
        self.weights = weights
        self.connectivity = connectivity
        self.spans = spans
        self.material = material
        self.properties = properties
