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

class _OutputRequest(Item):

    def __init__(self, label, name, nodal_variables, elemental_variables):
        Item.__init__(self, label, name)
        self._nodal_variables = nodal_variables
        self._elemental_variables = elemental_variables

    def get_group(self):
        return self._group

    def get_nodal_variables(self):
        return self._nodal_variables

    def get_elemental_variables(self):
        return self._elemental_variables

    def set_group(self, value):
        self._group = value

    def set_nodal_variables(self, value):
        self._nodal_variables = value

    def set_elemental_variables(self, value):
        self._elemental_variables = value


class FieldOutputRequest(_OutputRequest):

    def __init__(self, label, name, nodal_variables=[], elemental_variables=[]):
        _OutputRequest.__init__(self, label, name, nodal_variables, elemental_variables)

class HistoryOutputRequest(_OutputRequest):

    def __init__(self, label, name, nodal_variables=[], elemental_variables=[], group=None):
        _OutputRequest.__init__(self, label, name, nodal_variables, elemental_variables)
        self._group = group

    def get_group(self):
        return self._group

    def set_group(self, value):
        self._group = value


