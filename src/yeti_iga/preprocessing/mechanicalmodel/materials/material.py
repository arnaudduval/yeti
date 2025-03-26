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

from ..common import Item

class _Material(Item):
    def __init__(self, label, name, density=0., isconstant=True, type=None):
        Item.__init__(self, label, name)
        self._density = density
        self._type = type
        self._isconstant = isconstant

    def get_density(self):
        return self._density

    def set_density(self, density):
        self._density = density

    def get_label(self):
        return self._label

    def get_type(self):
        return self._type

    def set_type(self, value):
        self._type = value

    def link2distribution(self,distributionname):
        self._isconstant = False
        self._distributionname = distributionname

    def get_distributionname(self):
        if self._isconstant==True:
            return None
        else:
            return self._distributionname


class Isotropic(_Material):

    def __init__(self, label, name, density=0., e=0., nu=0.):
        _Material.__init__(self, label, name, density, type='Isotropic')
        self._e = e
        self._nu = nu

    def get_e(self):
        return self._e

    def get_nu(self):
        return self._nu

    def set_e(self, value):
        self._e = value

    def set_nu(self, value):
        self._nu = value

    def __str__(self):
        return 'Material Isotropic %s : label=%i, E=%1.2f, Nu=%1.1f' \
            % (self.get_name(), self.get_label(), self.get_e(), self.get_nu())

class Orthotropic(_Material):

    def __init__(self, label, name, density=0., ex=0., ey=0., nuxy=0., gxy=0., gxz=0.,
                 gyz=0., theta=None, type='Orthotropic'):
        _Material.__init__(self, label, name, density)
        self._ex = ex
        self._ey = ey
        self._nuxy = nuxy
        self._gxy = gxy
        self._gxz = gxz
        self._gyz = gyz
        self._theta = theta

    def get_ex(self):
        return self._ex

    def get_ey(self):
        return self._ey

    def get_nuxy(self):
        return self._nuxy

    def get_gxy(self):
        return self._gxy

    def get_gxz(self):
        return self._gxz

    def get_gyz(self):
        return self._gyz

    def get_theta(self):
        return self._theta

    def set_ex(self, value):
        self._ex = value

    def set_ey(self, value):
        self._ey = value

    def set_nuxy(self, value):
        self._nuxy = value

    def set_gxy(self, value):
        self._gxy = value

    def set_gxz(self, value):
        self._gxz = value

    def set_gyz(self, value):
        self._gyz = value

    def set_theta(self, value):
        self._theta = value

    def __str__(self):
        return 'Material Orthotropic %s : label=%i, Ex=%1.2f, Ey=%1.2f, Nuxy=%1.2f, Gxy=%1.2f, Gxz=%1.2f, Gyz=%1.2f, theta=%1.2f' % (self.get_name(), self.get_label(), self.get_ex(), self.get_ey(), self.get_nuxy(), self.get_gxy(), self.get_gxz(), self.get_gyz(), self.get_theta())

class HighOrderElastic(_Material):
    """
    Class representing a material with high order derivatives of strain
    dependant elatic bahaviour law.
    See Mindlin 1964 publication
    Derives from _Material class

    Attributes
    ----------
    _lambda : float
        Lamé coefficient _lambda
    _mu : float
        Lamé coefficient mu
    _a : [float]
        coefficients a_1 to a_5 soterd in a list (size = 5)
    _b : [float]
        coefficients b_0 to b_7 stored in a list (size = 8)
    _c : [float]
        coefficients c_1 to c_3 stored in a list (size = 3)

    Methods
    -------
    TODO
    """
    def __init__(self, label, name, density=None, lambd=None, mu=None, a=None, b=None, c=None):
        _Material.__init__(self, label, name, density, type='HighOrderElastic')
        self._lambda = lambd
        self._mu = mu
        self._a = a
        self._b = b
        self._c = c

    def get_lambda(self):
        return self._lambda

    def get_mu(self):
        return self._mu

    def get_a(self):
        return self._a

    def get_b(self):
        return self._b

    def get_c(self):
        return self._c

    def set_lambda(self, value):
        self._lambda = value

    def set_mu(self, value):
        self._mu = value

    def set_a(self, values):
        assert len(values) == 5
        self._a = values

    def set_b(self, values):
        assert len(values) == 8
        self._b = values

    def set_c(self, values):
        assert len(values) == 3
        self._c = values

    # TODO definir __str__
