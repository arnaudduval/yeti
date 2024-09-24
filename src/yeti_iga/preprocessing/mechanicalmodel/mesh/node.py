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

class Node(Item):

    def __init__(self, label, coordinates, misc={'assembly':False}):
        Item.__init__(self, label, misc=misc)
        self._x = coordinates[0]
        self._y = coordinates[1]
        self._z = coordinates[2]

    def get_label(self):
        return self._label

    def set_label(self,label):
        self._label = label

    def get_x(self):
        return self._x

    def get_y(self):
        return self._y

    def get_z(self):
        return self._z

    def get_coordinates(self):
        return np.array([self._x,self._y,self._z])

    def set_coordinates(self,coords):
        self.set_x(coords[0])
        self.set_y(coords[0])
        self.set_z(coords[0])

    def dist(self, other):
        return np.linalg.norm(self.get_coordinates()-other.get_coordinates())

    def set_x(self, x):
        self._x = x

    def set_y(self, y):
        self._y = y

    def set_z(self, z):
        self._z = z

    def __str__(self):
        return 'Node %i : coordinates=(%1.2f,%1.2f,%1.2f)' % (self.get_label(), self.get_x(), self.get_y(), self.get_z())

    def __repr__(self):
        return 'Node %i : coordinates=(%1.2f,%1.2f,%1.2f)' % (self.get_label(), self.get_x(), self.get_y(), self.get_z())


class NodeParam(Node):

    def __init__(self, label, coordinates, parametrization,misc={'assembly':False}):
        Node.__init__(self, label, coordinates, misc=misc)

        if ('<' in parametrization[0]) and ('>' in parametrization[0]):
            self._isxparam = True
            self._xstr = parametrization[0].translate({ord('<'): None,ord('>'): None})
        else:
            self._isxparam = False
        if ('<' in parametrization[1]) and ('>' in parametrization[1]):
            self._isyparam = True
            self._ystr = parametrization[1].translate({ord('<'): None,ord('>'): None})
        else:
            self._isyparam = False
        if ('<' in parametrization[2]) and ('>' in parametrization[2]):
            self._iszparam = True
            self._zstr = parametrization[2].translate({ord('<'): None,ord('>'): None})
        else:
            self._iszparam = False

    def _update_x(self,parameters):
        if not type(parameters)==type({}):
            raise ValueError('Input should be a dictionary')
        if self._isxparam:
            self.set_x(float(parameters[self._xstr]))

    def _update_y(self,parameters):
        if not type(parameters)==type({}):
            raise ValueError('Input should be a dictionary')
        if self._isyparam:
            self.set_y(float(parameters[self._ystr]))

    def _update_z(self,parameters):
        if not type(parameters)==type({}):
            raise ValueError('Input should be a dictionary')
        if self._iszparam:
            self.set_z(float(parameters[self._zstr]))

    def update_coords(self,parameters):
        if not type(parameters)==type({}):
            raise ValueError('Input should be a dictionary')
        self._update_x(parameters)
        self._update_y(parameters)
        self._update_z(parameters)


    def get_parametrization(self):
        parametrization = []
        if self._isxparam:
            parametrization.append("<"+self._xstr+">")
        else:
            parametrization.append(str(self.get_x()))
        if self._isyparam:
            parametrization.append("<"+self._ystr+">")
        else:
            parametrization.append(str(self.get_y()))
        if self._iszparam:
            parametrization.append("<"+self._zstr+">")
        else:
            parametrization.append(str(self.get_z()))
        return parametrization

