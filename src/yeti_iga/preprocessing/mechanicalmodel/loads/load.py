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

class _Load(Item):
    def __init__(self, label, name, target=None):
        Item.__init__(self, label, name)
        self._target = target

    def get_target(self):
        return self._target

    def set_target(self, value):
        self._target = value


class Load(_Load):
    def __init__(self, label, name, target=None, type='unknown', ADLMAG=0, JDLTYPE=0):
        _Load.__init__(self, label, name, target)
        self._type    = type
        self._ADLMAG  = 0
        self._JDLTYPE = 0
        self._additionalLoadInfos = np.array([])

    def get_type(self):
        return self._type

    def get_ADLMAG(self):
        return self._ADLMAG

    def get_JDLTYPE(self):
        return self._JDLTYPE

    def get_values(self):
        return [self._ADLMAG,self._JDLTYPE]

    def get_additionalLoadInfos(self):
        return self._additionalLoadInfos

    def set_ADLMAG(self,value):
        self._ADLMAG  = value

    def set_JDLTYPE(self,value):
        self._JDLTYPE = value

    def set_additionalInfos(self,array):
        self._additionalLoadInfos = np.asarray(array)
