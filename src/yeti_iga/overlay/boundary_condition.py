# !/usr/bin/env python
# -*- coding: utf-8 -*-

class BoundaryCondition:
    def __init__(self, cp_index, dof, value):
        """
        Parameters
        ----------
        cp_index : numpy.array([], dtype=int)
            local indices of control points on which boundary condition is applied
        dof : numpy.array([], dtype=int)
            indices of degrees of freedom on which the boundary condition is applied
        value : float
            prescribed value of degree of freedom
        """
        self.cp_index = cp_index
        self.dof = dof
        self.value = value