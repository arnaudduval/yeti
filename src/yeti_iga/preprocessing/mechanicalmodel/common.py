# Copyright 2016-2021 Thibaut Hirschler

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

import math
import numpy as np

class Item(object):

    def __init__(self, label, name=None, misc=None):
        self._label = label
        self._name = name
        self._misc = misc

    def get_label(self):
        return self._label

    def get_name(self):
        return self._name

    def get_misc(self):
        return self._misc

    def set_label(self, label):
        self._label = label

    def set_name(self, name):
        self._name = name

    def set_misc(self, name):
        self._misc = name

class Container(dict):

    def __init__(self):
        dict.__init__(self)

    def __getitem__(self, key):
        if isinstance(key, int):
            # Case of label
            return dict.__getitem__(self, key)
        elif isinstance(key, str):
            # Case of name
            for item in list(self.values()):
                if item.get_name() == key:
                    return dict.__getitem__(self, item.get_label())

    def __contains__(self, key):
        if isinstance(key, int):
            dict.__contains__(self, key)
        elif isinstance(key, str):
            for item in list(self.values()):
                if item.get_name() == key:
                    return True
            return False

    def __iter__(self):
        for item in dict.__iter__(self):
            yield dict.__getitem__(self, item)

    def add(self, item):
        if hasattr(item, '__iter__'):
            # Case of iterable
            for oneitem in item:
                dict.__setitem__(self, oneitem.get_label(), oneitem)
        else:
            # Case of one item
            dict.__setitem__(self, item.get_label(), item)

    def split(self):
        groups = {}
        for item in self:
            if type(item) not in groups:
                groups[type(item)] = [item]
            else:
                groups[type(item)].append(item)
        return groups

    def get_labels(self):
        return [item.get_label() for item in self]

    def get_names(self):
        return [item.get_name() for item in self]

    def get_one(self):
        key, item = self.popitem()
        self.add(item)
        return item


class TransformationMatrix(object):

    def __init__(self):
        self._matrix = np.identity(4)

    def get(self):
        return self._matrix

    def translate(self, v):
        x, y, z = v
        mat = np.array([[1, 0, 0, x],
                        [0, 1, 0, y],
                        [0, 0, 1, z],
                        [0, 0, 0, 1]], dtype=np.float64)
        self._matrix = np.dot(self._matrix, mat)

    def translate_x(self, x):
        v = np.array([x,0,0])
        self.translate(v)

    def translate_y(self, y):
        v = np.array([0,y,0])
        self.translate(v)

    def translate_z(self, z):
        v = np.array([0,0,z])
        self.translate(v)

    def rotate_x(self, theta):
        c = math.cos(theta)
        s = math.sin(theta)
        mat = N.array([[1, 0, 0, 0],
                       [0, c, -s, 0],
                       [0, s, c, 0],
                       [0, 0, 0, 1]],dtype=np.float64)
        self._matrix = N.dot(self._matrix, mat)

    def rotate_y(self, theta):
        c = math.cos(theta)
        s = math.sin(theta)
        mat = N.array([[c, 0, s, 0],
                       [0, 1, 0, 0],
                       [-s, 0, c, 0],
                       [0, 0, 0, 1]],dtype=np.float64)
        self._matrix = N.dot(self._matrix, mat)

    def rotate_z(self, theta):
        c = math.cos(theta)
        s = math.sin(theta)
        mat = N.array([[c, -s, 0, 0],
                       [s, c, 0, 0],
                       [0, 0, 1, 0],
                       [0, 0, 0, 1]],dtype=np.float64)
        self._matrix = N.dot(self._matrix, mat)

    def __str__(self):
        return str(self._matrix)


class Table(Item):

    def __init__(self, label, name, data=np.empty(0), defaultvalues=0., opts={}):
        Item.__init__(self, label, name)
        self._data = np.array(data)
        self.opts = opts
        self.set_defaultvalues(defaultvalues)

    def append(self, values, axis=0):
        if len(self._data) == 0:
            self._data = np.append(self._data, values).reshape((1, np.array(values).size))
        else:
            self._data = np.append(self._data, values, axis)

    def set_defaultvalues(self, values):
        self._defaultdata = np.ascontiguousarray(values,dtype=float)

    def get_defaultvalues(self):
        return np.copy(self._defaultdata)

    def insert(self, pos, values, axis=0):
        self._data = np.insert(self._data, pos, values, axis)

    def delete(self, pos, axis=0):
        self._data = np.delete(self._data, pos, axis)

    def shape(self):
        return self._data.shape

    def values(self):
        for x in self._data:
            yield x

    def __getitem__(self, pos):
        if isinstance(pos, str):
            return self.opts[pos]
        else:
            return self._data.__getitem__(pos)

    def __len__(self):
        return len(self._data)

    def __repr__(self):
        return str(self._data)

    def __str__(self):
        return str(self._data)

