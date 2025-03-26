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

from numpy import size
from ..common import Item

class ThinShell(Item):

    def __init__(self, label, name, thickness=None, nint=5, offset=0, elset=None):
        Item.__init__(self, label, name)
        self._thickness = thickness
        self._nint = nint
        self._offset = offset
        self._elset = elset

    def get_thickness(self):
        return self._thickness

    def get_nint(self):
        return self._nint

    def get_offset(self):
        return self._offset

    def get_elset(self):
        return self._elset

    def set_thickness(self, value):
        self._thickness = value

    def set_nint(self, value):
        self._nint = nint

    def set_offset(self, value):
        self._offset = value

    def set_elset(self, value):
        self._elset = value

    def __str__(self):
        if isinstance(self._offset, float):
            return 'Property Thin shell %s : label=%i, thickness=%1.2f, offset=%1.1f' % (self.get_name(), self.get_label(), self.get_thickness(), self.get_offset())
        else:
            return 'Property Thin shell %s : label=%i, thickness=%1.2f, offset=%s' % (self.get_name(), self.get_label(), self.get_thickness(), self.get_offset().get_name())

class CompositeThinShell(ThinShell):

    def __init__(self, label, name, layup=None, elset=None):
        ThinShell.__init__(self, label=label, name=name, elset=elset)
        self.set_layup(layup)

    def get_layup(self):
        return self._layup

    def set_layup(self, layup):
        self._layup = layup
        t = 0
        for ply in self._layup:
            t += ply.get_thickness()
        self.set_thickness(t)

    def set_thickness(self, thickness):
        self._thickness = thickness

    def __str__(self):
        return 'Composite Thin shell %s : label=%i, layup=%s, thickness=%1.2f' % (self.get_name(), self.get_label(), self.get_layup().get_name(), self.get_thickness())

class Ply(object):
    def __init__(self, name=None, thickness=None, orientation=None, material=None, nint=3):
        self._name = name
        self._thickness = thickness
        self._orientation = orientation
        self._material = material
        self._nint = nint

    def get_name(self):
        return self._name

    def get_thickness(self):
        return self._thickness

    def get_orientation(self):
        return self._orientation

    def get_material(self):
        return self._material

    def get_nint(self):
        return self._nint

    def set_name(self, value):
        self._name = value

    def set_thickness(self, value):
        self._thickness = value

    def set_orientation(self, value):
        self._orientation = value

    def set_material(self, value):
        self._material = value

    def set_nint(self, value):
        self._nint = value

class CompositeLayup(list):

    def __init__(self, label=None, name=None):
        list.__init__(self)
        self._label = label
        self._name = name

    def get_label(self):
        return self._label

    def get_name(self):
        return self._name

    def set_label(self, value):
        self._label = value

    def set_name(self, value):
        self._name = value

class Solid(Item):

    def __init__(self, label=None, name=None, elset=None):
        Item.__init__(self, label, name)
        self._elset = elset

    def get_elset(self):
        return self._elset

    def set_elset(self, value):
        self._elset = value

    def __str__(self):
        return 'Solid property %s : label=%i' % (self.get_name(), self.get_label())




class userElementProps(Item):

    def __init__(self, label=None, name=None, elset=None, props=None, material=None):
        Item.__init__(self, label, name)
        self._props    = props
        self._elset    = elset
        self._material = material

    def get_props(self):
        return self._props

    def get_jprops(self):
        return size(self.get_props())

    def get_elset(self):
        return self._elset

    def get_material(self):
        return self._material

    def set_props(self, value):
        self._props = value

    def set_elset(self, value):
        self._elset = value

    def set_material(self, value):
        self._material = value

