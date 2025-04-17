# !/usr/bin/env python
# -*- coding: utf-8 -*-

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
        for i in range(nb_patch):
            self.perpatch.append({'degree_elevation': np.array([0, 0, 0], dtype=int),
                                  'subdivision': np.array([0, 0, 0], dtype=int)
                                  })