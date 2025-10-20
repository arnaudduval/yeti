"""
Patch geometry templates
"""

import numpy as np
from math import sqrt
from . import Patch


class QuarterPipe(Patch):
    """
    A patch with a quarter pipe geometry
    """

    def __init__(self, length, inner, outer, material):
        """
        Parameters
        ----------
        length: float
            quarter pipe length
        inner: float
            quarter pipe inner radius
        outer: float
            quarter pipe outer radius
        material: Material
            material of the quarter pipe
        """

        super().__init__(element_type='U1',
                         degrees=np.array([1, 2, 1]),
                         knot_vectors=[np.array([0., 0., 1., 1.]),
                                       np.array([0., 0., 0., 1., 1., 1.]),
                                       np.array([0., 0., 1., 1.])],
                         control_points=np.array([[inner, 0., 0.],
                                                  [outer, 0., 0.],
                                                  [inner, inner, 0.],
                                                  [outer, outer, 0.],
                                                  [0., inner, 0.],
                                                  [0., outer, 0.],
                                                  [inner, 0., length],
                                                  [outer, 0., length],
                                                  [inner, inner, length],
                                                  [outer, outer, length],
                                                  [0., inner, length],
                                                  [0., outer, length]]),
                         weights=np.array([1., 1., 1./sqrt(2.), 1./sqrt(2.), 1., 1., 1., 1., 1./sqrt(2.), 1./sqrt(2.), 1., 1.]),
                         connectivity=np.array([[11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0]]),
                         spans=np.array([[1, 2, 1]]),
                         material=material)

