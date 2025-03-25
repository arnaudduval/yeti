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

from ..common import Item, Container
import numpy as np

class _Step(Item):

    def __init__(self, label, name, family=None, type=None):
        Item.__init__(self, label, name)
        self._family = family
        self._type = type
        self.bcs = Container()
        self.loads = Container()
        self.outputs = Container()

    def get_family(self):
        return self._family

    def get_type(self):
        return self._type

    def set_family(self, value):
        self._family = value

    def set_type(self, value):
        self._type = value


class _GeneralStep(_Step):

    def __init__(self, label, name, type=None):
        _Step.__init__(self, label, name, family='general', type=type)

class InitialStep(_GeneralStep):
    def __init__(self, label, name):
        _GeneralStep.__init__(self, label, name, type='initial')

    def __str__(self):
        return 'Initial step : label=%i' % self.get_label()

class GeneralStaticStep(_GeneralStep):

    def __init__(self, label, name, period=1.0, non_linear=False, nincmax=None, increment_init=1.0,
                 increment_min=1e-5, increment_max=1.0):
        _GeneralStep.__init__(self, label, name, type='general_static')
        self._period         = period
        self._non_linear     = non_linear
        self._nincmax        = nincmax
        self._increment_init = increment_init
        self._increment_min  = increment_min
        self._increment_max  = increment_max

    def get_period(self):
        return self._period

    def get_non_linear(self):
        return self._non_linear

    def get_nincmax(self):
        return self._nincmax

    def get_increment_init(self):
        return self._increment_init

    def get_increment_min(self):
        return self._increment_min

    def get_increment_max(self):
        return self._increment_max

    def set_period(self, value):
        self._period = value

    def set_non_linear(self, value):
        self._non_linear = value

    def set_nincmax(self, value):
        self._nincmax = value

    def set_increment_init(self, value):
        self._increment_init = value

    def set_increment_min(self, value):
        self._increment_min = value

    def set_increment_max(self, value):
        self._increment_max = value

    def __str__(self):
        return 'General static step : label=%i, name=%s, number_bcs=%i, number_loads=%i, period=%1.1f, non_linear=%s, inc_init=%1.1e, inc_min=%1.1e, inc_max=%1.1e' % (self.get_label(), self.get_name(), len(self.bcs), len(self.loads), self.get_period(), self.get_non_linear(), self.get_increment_init(), self.get_increment_min(), self.get_increment_max())

class RiksStaticStep(_GeneralStep):

    def __init__(self, label, name):
        _GeneralStep.__init__(self, label, name, type='riks_static')

class _LinearPerturbationStep(_Step):

    def __init__(self, label, name):
        _Step.__init__(self, label, name, family='linear_perburbation')


class BuckleStep(_LinearPerturbationStep):

    def __init__(self, label, name, neigenvalues=1, maxeigenvalues=None, nvectors=2, maxiter=30):
        _LinearPerturbationStep.__init__(self, label, name, type='buckle')
        self._neigenvalues = neigenvalues
        self._maxeigenvalues = maxeigenvalues
        self._nvectors = nvectors
        self._maxiter = maxiter

    def get_neigenvalues(self):
        return self._neigenvalues

    def get_maxeigenvalues(self):
        return self._maxeigenvalues

    def get_nvectors(self):
        return self._nvectors

    def get_maxiter(self):
        return self._maxiter

    def set_neigenvalues(self, value):
        self._neigenvalues = value

    def set_maxeigenvalues(self, value):
        self._maxeigenvalues = value

    def set_nvectors(self, value):
        self._nvectors = value

    def set_maxiter(self, value):
        self._maxiter = value

    def __str__(self):
        return 'Perturbation buckle step : label=%i, name=%s, number_bcs=%i, number_loads=%i, number_eigenvalues=%i, number_vectors=%i, max_iter=%i' % (self.get_label(), self.get_name(), len(self.bcs), len(self.loads), self.get_neigenvalues(), self.get_nvectors(), self.get_maxiter())


