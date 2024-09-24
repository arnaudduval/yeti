# Copyright 2016-2018 Thibaut Hirschler

# This file is part of Yeti.
#
# Yeti is free software: you can redistribute it and/or modify it under the terms
# of the GNU Lesser General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later version.
#
# Yeti is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
# PURPOSE. See the GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License along
# with Yeti. If not, see <https://www.gnu.org/licenses/>

from ..common import Item
import numpy as np

class _BoundaryCondition(Item):

    def __init__(self, label, name, target=None):
        Item.__init__(self, label, name)
        self._target = target

    def get_target(self):
        return self._target

    def set_target(self, value):
        self._target = value

class _MechanicalBoundaryCondition(_BoundaryCondition):

    def __init__(self, label, name, target=None, x=None, y=None, z=None, rx=None, ry=None, rz=None):
        _BoundaryCondition.__init__(self, label, name, target)
        self._x = x
        self._y = y
        self._z = z
        self._rx = rx
        self._ry = ry
        self._rz = rz

    def get_x(self):
        return self._x

    def get_y(self):
        return self._y

    def get_z(self):
        return self._z

    def get_rx(self):
        return self._rx

    def get_ry(self):
        return self._ry

    def get_rz(self):
        return self._rz

    def get_values(self):
        return [self._x, self._y, self._z, self._rx, self._ry, self._rz]

    def get_use_value(self):
        A=np.zeros(2) # first value : direction of imposed displacement #2nd value : value of the displacement
        if (self._x!=None):
            A=[1,self._x]
        elif self._y!=None :
            A=[2,self._y]
        elif self._z!=None :
            A=[3,self._z]

        return A

    def set_x(self, value):
        self._x = value

    def set_y(self, value):
        self._y = value

    def set_z(self, value):
        self._z = value

    def set_rx(self, value):
        self._rx = value

    def set_ry(self, value):
        self._ry = value

    def set_rz(self, value):
        self._rz = value

    def _values_to_str(self):
        s = '('
        for v in self.get_values():
            if v == None:
                s += '-,'
            else:
                s += '%1.1e,' % v
        s = s[:-1] + ')'
        return s

class DisplacementBoundaryCondition(_MechanicalBoundaryCondition):

    def __init__(self, label, name, target=None, x=None, y=None, z=None, rx=None, ry=None, rz=None):
        _MechanicalBoundaryCondition.__init__(self, label, name, target, x, y, z, rx, ry, rz)

    def __str__(self):
        return 'Displacement boundary condition : label=%i, name=%s, target=%s, values=%s' % (self.get_label(), self.get_name(), self.get_target().get_name(), self._values_to_str())

class VelocityBoundaryCondition(_MechanicalBoundaryCondition):

    def __init__(self, label, name, target=None, x=None, y=None, z=None, rx=None, ry=None, rz=None):
        _MechanicalBoundaryCondition.__init__(self, label, name, target, x, y, z, rx, ry, rz)

    def __str__(self):
        return 'Velocity boundary condition : label=%i, name=%s, target=%s, values=%s' % (self.get_label(), self.get_name(), self.get_target().get_name(), self._values_to_str())

class AccelerationBoundaryCondition(_MechanicalBoundaryCondition):

    def __init__(self, label, name, target=None, x=None, y=None, z=None, rx=None, ry=None, rz=None):
        _MechanicalBoundaryCondition.__init__(self, label, name, target, x, y, z, rx, ry, rz)

    def __str__(self):
        return 'Acceleration boundary condition : label=%i, name=%s, target=%s, values=%s' % (self.get_label(), self.get_name(), self.get_target().get_name(), self._values_to_str())

