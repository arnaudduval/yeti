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

from .common import Item, Container

class Group(Item):

    def __init__(self, label, name, misc={'assembly':False, 'internal':False, 'from_instance':False}, type=None):

        Item.__init__(self, label, name, misc=misc)
        self.nodes = Container()
        self.elements = Container()
        self._type = type

    def get_type(self):
        return self._type

    def set_type(self, type):
        self._type = type

    def __str__(self):
        return 'Group %i : name = %s, %i nodes and %i elements, infos=%s' % (self.get_label(), self.get_name(), len(self.nodes), len(self.elements), self.get_infos())

    def get_nodes_labels(self):
        return self.nodes.get_labels()


    def get_elts_labels(self):
        return self.elements.get_labels()



#     def extend_nodes(self):
#
#         _nodes = [_node for _element in self.elements for _node in _element.nodes]
#         self.nodes.add_multiple(_nodes)
#
#     def extend_elements(self,steps=-1):
#
#         n = 0
#         l1 = len(self.elements)
#         #for element in self.elements:
#         #    element.set_neighbours()
#         while True:
#             _elements = set()
#             for element in self.elements:
#                 for _neighbour in element.neighbours:
#                     if _neighbour not in self.elements:
#                         _elements.add(_neighbour)
#             #for element in _elements:
#             #    element.set_neighbours()
#             self.elements.add_multiple(_elements)
#             l2 = len(self.elements)
#             n += 1
#             if (steps > 0 and n >= steps) or (l1 == l2):
#                 break
#             else:
#                 l1 = l2
#
#     def calculate_area(self,force=False):
#
#         area = 0.
#         if force:
#             for element in self.elements:
#                 area += element.calculate_area()
#         else:
#             for element in self.elements:
#                 area += element.get_area()
#         self._area = area
#         return self._area
#
#     def calculate_surface(self,force=False):
#
#         return self.calculate_area(force)
#
#     def get_area(self):
#
#         try:
#             return self._area
#         except AttributeError:
#             return self.calculate_area()
#
#     def get_surface(self):
#
#         return get_area()
#
#     def get_volume(self):
#         volume = 0.
#         for element in self.elements:
#             volume += element.get_volume()
#         return volume
#
#     def get_mass(self):
#         mass = 0.
#         for element in self.elements:
#             mass += element.get_mass()
#         return mass
#
#     def transformation(self,grid,matrix):
#
#         filter = [ node.grids[grid] for node in self.nodes ]
#         grid.transformation(matrix,filter)
#
#     def _get_info(self):
#
#         info = db.Item._get_info(self)
#         try:
#             info.append(('area',self._area))
#         except AttributeError:
#             pass
#         return info
#
#     def __or__(self,other):
#
#         _group = Group()
#         _group.elements.add_multiple(self.elements | other.elements)
#         _group.nodes.add_multiple(self.nodes | other.nodes)
#         return _group
#
#     def __ior__(self,other):
#
#         self.elements.__ior__(other.elements)
#         self.nodes.__ior__(other.nodes)
#         return self
#
#     def __and__(self,other):
#
#         _group = Group()
#         _group.elements.add_multiple(self.elements & other.elements)
#         _group.nodes.add_multiple(self.nodes & other.nodes)
#         return _group
#
#     def __iand__(self,other):
#
#         self.elements.__iand__(other.elements)
#         self.nodes.__iand__(other.nodes)
#         return self
#
#     def __sub__(self,other):
#
#         _group = Group()
#         _group.elements.add_multiple(self.elements - other.elements)
#         _group.nodes.add_multiple(self.nodes - other.nodes)
#         return _group
#
#     def __isub__(self,other):
#
#         self.elements.__isub__(other.elements)
#         self.nodes.__isub__(other.nodes)
#         return self
#
#     def __xor__(self,other):
#
#         _group = Group()
#         _group.elements.add_multiple(self.elements ^ other.elements)
#         _group.nodes.add_multiple(self.nodes ^ other.nodes)
#         return _group
#
#     def __ixor__(self,other):
#
#         self.elements.__ixor__(other.elements)
#         self.nodes.__ixor__(other.nodes)
#         return self
#
#     def __eq__(self,other):
#
#         return self.elements == other.elements and self.nodes == other.nodes
#
#     def __ne__(self,other):
#
#         return self.elements != other.elements and self.nodes != other.nodes
#
#     def __le__(self,other):
#
#         return self.elements <= other.elements and self.nodes <= other.nodes
#
#     def __lt__(self,other):
#
#         return self.elements < other.elements and self.nodes < other.nodes
#
#     def __ge__(self,other):
#
#         return self.elements >= other.elements and self.nodes >= other.nodes
#
#     def __gt__(self,other):
#
#         return self.elements > other.elements and self.nodes > other.nodes
#
#     def isdisjoint(self,other):
#         return self.elements.isdisjoint(other.elements)
#
#     def _child_to_str(self, child='elements'):
#
#         labels = [c.get_label() for c in getattr(self,child)]
#         labels.sort()
#         n = len(labels)
#         start = None
#         last = None
#         ranges = []
#         for i, label in enumerate(labels):
#             if start is None:
#                 start = label
#                 last = label
#             else:
#                 if i == n-1:
#                     if label == last + 1:
#                         ranges.append("%i:%i" % (start, label))
#                     else:
#                         if last == start:
#                             ranges.append(str(start))
#                         else:
#                             ranges.append("%i:%i" % (start, last))
#                         ranges.append(str(label))
#                 else:
#                     if (label == last + 1):
#                         last = label
#                     else:
#                         if last == start:
#                             ranges.append(str(start))
#                         else:
#                             ranges.append("%i:%i" % (start, last))
#                         start = label
#                         last = label
#         return ' '.join(ranges)
#
#     def elements_to_str(self):
#
#         return self._child_to_str(child='elements')
#
#     def nodes_to_str(self):
#
#         return self._child_to_str(child='nodes')
