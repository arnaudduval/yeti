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

import numpy as np
from math import sqrt

from ..common import Item

class _Element(Item):

    def __init__(self, label, nodes, property=None, material=None, type=None):
        Item.__init__(self, label)
        self._nodes    = nodes
        self._property = property
        self._material = material
        self._dimension= None
        self._order    = None
        self._type     = type

    def get_nodes(self):
        return self._nodes

    def get_type(self):
        return self._type

    def get_node_labels(self):
        return [node.get_label() for node in self._nodes]

    def get_dimension(self):
        return self._dimension

    def get_order(self):
        return self._order

    def get_property(self):
        return self._property

    def get_material(self):
        return self._material

    def get_volume(self):
        pass

    def get_mass(self):
        return self.get_volume()*self.get_material().get_density()

    def get_centre(self):
        n = len(self.get_nodes())
        xc = yc = zc = 0
        for node in self.get_nodes():
            x,y,z = node.get_coordinates()
            xc += x
            yc += y
            zc += z
        xc /= n
        yc /= n
        zc /= n
        return xc,yc,zc

    def set_nodes(self, nodes):
        self._nodes = nodes

    def set_type(self,elt_type):
        self._type = elt_type

    def set_property(self, property):
        self._property = property

    def set_material(self, material):
        self._material = material


class _OneDimensionElement(_Element):

    def __init__(self, label, nodes, property=None, material=None):
        _Element.__init__(self, label=label, nodes=nodes, property=property, material=material)
        self._dimension = 1

    def get_length(self):
        n1, n2 = self._nodes
        return n1.dist(n2)

    def get_volume(self):
        return self.get_length()*self.get_property().get_area()

class _TwoDimensionElement(_Element):

    def __init__(self, label, nodes, property=None, material=None):
        _Element.__init__(self, label=label, nodes=nodes, property=property, material=material)
        self._dimension = 2

    def get_volume(self):
        return self.get_area()*self.get_property().get_thickness()

    def get_normal(self):
        # Let n1,n2,n3 the 3 first nodes of element definition
        # Normal vector n = vectorial product between n1n2 and n2n3
        n1 = self._nodes[0].get_coordinates()
        n2 = self._nodes[1].get_coordinates()
        n3 = self._nodes[2].get_coordinates()
        n1n2 = n2-n1
        n2n3 = n3-n2
        normal = np.cross(n1n2, n2n3)
        return normal/np.linalg.norm(normal)

class _ThreeDimensionElement(_Element):

    def __init__(self, label, nodes, property=None, material=None):
        _Element.__init__(self, label=label, nodes=nodes, property=property, material=material)
        self._dimension = 3

class Beam(_OneDimensionElement):

    def __init__(self, label, nodes, property=None, material=None):
        _OneDimensionElement.__init__(self, label=label, nodes=nodes, property=property, material=material)

    def __str__(self):
        return 'Element 1D Beam : label=%i, connectivity=(%i,%i)' % (self.get_label(), self.get_nodes()[0].get_label(), self.get_nodes()[1].get_label())

class Tria(_TwoDimensionElement):

    def __init__(self, label, nodes, property=None, material=None):
        _TwoDimensionElement.__init__(self, label=label, nodes=nodes, property=property, material=material)

    def get_area(self):
        n1, n2, n3 = self._nodes
        a = n1.dist(n2)
        b = n2.dist(n3)
        c = n3.dist(n1)
        s = 0.5*(a + b + c)
        return sqrt(s*(s-a)*(s-b)*(s-c))

    def __str__(self):
        return 'Element 2D Tria : label=%i, connectivity=(%i,%i,%i)' % (self.get_label(), self.get_nodes()[0].get_label(), self.get_nodes()[1].get_label(), self.get_nodes()[2].get_label())

class Quad(_TwoDimensionElement):

    def __init__(self, label, nodes, property=None, material=None):
        _TwoDimensionElement.__init__(self, label=label, nodes=nodes, property=property, material=material)

    def get_area(self):
        n1, n2, n3, n4 = self._nodes
        t1 = Tria(label=1, nodes=[n1,n2,n4])
        t2 = Tria(label=2, nodes=[n2,n3,n4])
        return t1.get_area() + t2.get_area()

    def __str__(self):
        return 'Element 2D Quad : label=%i, connectivity=(%i,%i,%i,%i)' % (self.get_label(), self.get_nodes()[0].get_label(), self.get_nodes()[1].get_label(), self.get_nodes()[2].get_label(), self.get_nodes()[3].get_label())

class Penta(_ThreeDimensionElement):

    def __init__(self, label, nodes, property=None, material=None):
        _ThreeDimensionElement.__init__(self, label=label, nodes=nodes, property=property, material=material)

    def __str__(self):
        return 'Element 3D Penta : label=%i, connectivity=(%i,%i,%i,%i,%i,%i,%i)' % (self.get_label(), self.get_nodes()[0].get_label(), self.get_nodes()[1].get_label(), self.get_nodes()[2].get_label(), self.get_nodes()[3].get_label(), self.get_nodes()[4].get_label(), self.get_nodes()[5].get_label())

class Tetra(_ThreeDimensionElement):

    def __init__(self, label, nodes, property=None, material=None):
        _ThreeDimensionElement.__init__(self, label=label, nodes=nodes, property=property, material=material)

    def get_volume(self):
        n1, n2, n3, n4  = self._nodes
        return 1/6.*(abs(np.det(np.array([
                n1.get_coordinates()- n4.get_coordinates(),
                n2.get_coordinates()- n4.get_coordinates(),
                n3.get_coordinates()- n4.get_coordinates()]))))

    def __str__(self):
        return 'Element 3D Tetra : label=%i, connectivity=(%i,%i,%i,%i,%i)' % (self.get_label(), self.get_nodes()[0].get_label(), self.get_nodes()[1].get_label(), self.get_nodes()[2].get_label(), self.get_nodes()[3].get_label(), self.get_nodes()[4].get_label())

class Hexa(_ThreeDimensionElement):

    def __init__(self, label, nodes, property=None, material=None):
        _ThreeDimensionElement.__init__(self, label=label, nodes=nodes, property=property, material=material)

    def get_volume(self):
        n1, n2, n3, n4, n5, n6, n7, n8 = self._nodes

        n1n7 = n7.get_coordinates() - n1.get_coordinates()
        return 1/6.*(
        abs(np.linalg.det(np.array([n1n7, n2.get_coordinates() - n1.get_coordinates(),
                                          n3.get_coordinates() - n6.get_coordinates()])))
        + \
        abs(np.linalg.det(np.array([n1n7, n5.get_coordinates() - n1.get_coordinates(),
                                          n6.get_coordinates() - n8.get_coordinates()])))
        + \
        abs(np.linalg.det(np.array([n1n7, n4.get_coordinates() - n1.get_coordinates(),
                                          n6.get_coordinates() - n3.get_coordinates()]))))

    def __str__(self):
        return 'Element 3D Hexa : label=%i, connectivity=(%i,%i,%i,%i,%i,%i,%i,%i)' % (self.get_label(), self.get_nodes()[0].get_label(), self.get_nodes()[1].get_label(), self.get_nodes()[2].get_label(), self.get_nodes()[3].get_label(), self.get_nodes()[4].get_label(), self.get_nodes()[5].get_label(), self.get_nodes()[6].get_label(), self.get_nodes()[7].get_label())



class UserElement(_Element):
    def __init__(self, label, nodes, property=None, material=None, type='U1'):
        _Element.__init__(self, label=label, nodes=nodes, property=property, material=material,
                          type=type)

class UserElementParam(Item):
    def __init__(self, label, name, nnode, coords=1,var=0,iprop=0,integration=1,tensor='THREED'):
        Item.__init__(self, label, name)
        self._nnode       = nnode
        self._coordinates = coords
        self._variables   = var
        self._iproperties = iprop
        self._integration = integration
        self._tensor      = tensor

    def get_nnode(self):
        return self._nnode
    def set_nnode(self,value):
        self._nnode=value

    def get_coordinates(self):
        return self._coordinates
    def set_coordinates(self,value):
        self._coordinates=value

    def get_variables(self):
        return self._variables
    def set_variables(self,value):
        self._variables=value

    def get_iproperties(self):
        return self._iproperties
    def set_iproperties(self,value):
        self._iproperties=value

    def get_integration(self):
        return self._integration
    def set_integration(self,value):
        self._integration=value

    def get_tensor(self):
        return self._tensor
    def set_tensor(self,value):
        self._tensor=value

