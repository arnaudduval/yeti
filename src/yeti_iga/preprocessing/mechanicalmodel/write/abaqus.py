# Copyright 2016-2019 Thibaut Hirschler

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

from .common import Writer
from ..mesh.node import Node
from ..mesh.element import Beam, Tria, Quad, Penta, Tetra, Hexa
from ..properties.property import ThinShell, CompositeThinShell, Solid
from ..materials.material import Isotropic, Orthotropic
from ..boundary_conditions.bc import DisplacementBoundaryCondition, VelocityBoundaryCondition, AccelerationBoundaryCondition
from ..steps import InitialStep, GeneralStaticStep, BuckleStep, FieldOutputRequest, HistoryOutputRequest

class AbaqusWriter(Writer):
    """
    Abaqus inp writer
    Limitations:
      - creates only one part and one instance for the whole model
      - part nodes coordinates are actually assembly nodes coordinates
      - loads not supported
    """

    def __init__(self, model, filename):
        Writer.__init__(self, model, filename)

    def write(self):
        print('\nWriting Abaqus inp file %s...\n' % self._filename)
        self._write_parts()
        self._write_assembly()
        self._write_materials()
        self._write_distribution_tables()
        self._write_steps()
        self._file.close()
        print('Successfully done.')

    def _write_parts(self):
        print('Writing parts...')
        #### SECTION TO BE CHANGED TO SUPPORT ALL PARTS AND INSTANCES ##########
        #### TEMPORARY SOLUTION ###
        self._partname = self._model.info['name']
        self._instancename = self._partname + '-1'
        #self._instancename = self._model.info['instances'][self._partname][0]
        #### TEMPORARY SOLUTION ###

        self._file.write('**\n** PARTS\n**\n')
        self._file.write('*Part, name=%s\n' % self._partname)
        self._write_nodes()
        self._write_elements()
        self._write_groups()
        self._write_properties()
        self._file.write('*End Part\n**\n')
        print('Parts successfully written.\n')

    def _write_nodes(self):
        print('\tWriting nodes...')
        self._file.write('*Node\n')
        for node in self._model.nodes:
            if not node.get_misc()['assembly']:
                self._file.write('%7i,%13f,%13f,%13f\n' % (node.get_label(), node.get_x(), node.get_y(), node.get_z()))
        print('%i nodes written.' % len(self._model.nodes))

    def _write_elements(self):
        print('\tWriting elements...')
        element_groups = self._model.elements.split()
        for type in element_groups:
            if type == Beam:
                self._file.write('*Element, type=T3D2\n')
                format = str('%6i,'*3)[:-1] + '\n'
            elif type == Tria:
                self._file.write('*Element, type=S3R\n')
                format = str('%6i,'*4)[:-1] + '\n'
            elif type == Quad:
                self._file.write('*Element, type=S4R\n')
                format = str('%6i,'*5)[:-1] + '\n'
            elif type == Hexa:
                self._file.write('*Element, type=C3D8R\n')
                format = str('%6i,'*9)[:-1] + '\n'
            for element in element_groups[type]:
                vars = (element.get_label(),) + tuple(element.get_node_labels())
                self._file.write(format % vars)
        print('%i elements written.' % len(self._model.elements))

    def _write_assembly(self):
        print('Writing assembly...')
        self._file.write('**\n** ASSEMBLY\n**\n*Assembly, name=Assembly\n')
        self._write_instances()
        self._write_assembly_nodes()
        self._write_assembly_groups()
        self._file.write('**\n*End Assembly\n')
        print('Assembly successfully written.\n')

    def _write_instances(self):
        print('\tWriting instances...')
        self._file.write('*Instance, name=%s, part=%s\n' % (self._instancename, self._partname))
        #self._write_groups()
        #self._write_properties()
        self._file.write('*End Instance\n')
        print(' 1 instances successfully written.')

    def _write_groups(self):
        print('\tWriting sets...')
        cpt = 0
        for group in self._model.groups:
            if not group.get_misc()['assembly']:
                if len(group.nodes) > 0:
                    #self._file.write('*Nset, nset=%s, instance=%s\n' % (groupname, instance)
                    self._file.write('*Nset, nset=%s\n' % group.get_name())
                    self._write_group(group.nodes)
                    cpt += 1
                if len(group.elements) > 0:
                    #self._file.write('*Elset, elset=%s, instance=%s\n' % (groupname, instance))
                    self._file.write('*Elset, elset=%s\n' % group.get_name())
                    self._write_group(group.elements)
                    cpt += 1
        print('%i sets successfully written.' % cpt)

    def _write_assembly_nodes(self):
        print('\tWriting nodes...')
        self._file.write('*Node\n')
        for node in self._model.nodes:
            if node.get_misc()['assembly']:
                self._file.write('%7i,%13f,%13f,%13f\n' % (node.get_label(), node.get_x(), node.get_y(), node.get_z()))
        print('%i nodes written.' % len(self._model.nodes))

    def _write_assembly_groups(self):
        print('\tWriting assembly groups...')
        cpt = 0
        for group in self._model.groups:
            if group.get_misc()['assembly']:
                if len(group.nodes) > 0:
                    if group.get_misc()['from_instance']:
                        self._file.write('*Nset, nset=%s, instance=%s\n' % (group.get_name(),self._instancename))
                    else:
                        self._file.write('*Nset, nset=%s\n' % (group.get_name()))
                    self._write_group(group.nodes)
                    cpt += 1
                if len(group.elements) > 0:
                    if group.get_misc()['from_instance']:
                        self._file.write('*Elset, elset=%s, instance=%s\n' % (group.get_name(),self._instancename))
                    else:
                        self._file.write('*Elset, elset=%s\n' % (group.get_name()))
                    self._write_group(group.elements)
                    cpt += 1
        print(' %i sets successfully written.' % cpt)

    def _write_group(self, group):
        labels = group.get_labels()
        subs = [labels[i:i+11] for i in range(0,len(labels),11)]
        format = str('%6i,'*11)[:-1] + '\n'
        for sub in subs:
            if len(sub) < 11:
                format = str('%6i,'*len(sub))[:-1] + '\n'
                self._file.write(format % tuple(sub))
            else:
                self._file.write(format % tuple(sub))

    def _write_properties(self):
        print('\tWriting properties...')
        cpt = 0
        for property in self._model.properties:
            if isinstance(property, CompositeThinShell):
                elset = property.get_elset()
                layup = property.get_layup()
                self._file.write('** Section: %s\n' % property.get_name())
                self._file.write('*Shell Section, elset=%s, composite, layup=%s\n' % (elset, layup.get_name()))
                for ply in layup:
                    self._file.write('%f, %i, %s, %1.1f, %s\n' % (ply.get_thickness(), ply.get_nint(), ply.get_material(), ply.get_orientation(), ply.get_name()))
                cpt += 1
            elif isinstance(property, ThinShell):
                elset = property.get_elset()
                offset = property.get_offset()
                t = property.get_thickness()
                nint = property.get_nint()
                #for group in self._model.groups:
                #    if elset in group.get_name():
                #        elsetname = group.get_name()
                #        break
                mat = self._model.groups[elset].elements.get_one().get_material().get_name()
                self._file.write('** Section: %s\n' % property.get_name())
                if offset:
                    if isinstance(offset, float):
                        self._file.write('*Shell Section, elset=%s, material=%s, offset=%1.1f\n' % (elset, mat, offset))
                        self._file.write('%1.1f, %i\n' % (t, nint))
                    else:
                        self._file.write('*Shell Section, elset=%s, material=%s, offset=%s\n' % (elset, mat, offset.get_name()))
                        self._file.write('%1.1f, %i\n' % (t, nint))
                        self._file.write('*Distribution, name=%s, location=%s, table=%s\n' % (offset.get_name(), offset['location'], offset['dtable']))
                        for label, value in list(offset.values()):
                            if int(label) == 0:
                                self._file.write(', %f\n' % value)
                            else:
                                self._file.write('%i, %f\n' % (label, value))
                else:
                    self._file.write('*Shell Section, elset=%s, material=%s\n' % (elset, mat))
                    self._file.write('%1.1f, %i\n' % (t, nint))
                cpt += 1
            elif isinstance(property, Solid):
                elset = property.get_elset()
                for group in self._model.groups:
                    if elset in group.get_name():
                        elsetname = group.get_name()
                        break
                mat = self._model.groups[elsetname].elements.get_one().get_material().get_name()
                self._file.write('** Section: %s\n' % property.get_name())
                self._file.write('*Solid Section, elset=%s, material=%s\n' % (property.get_elset(), mat))
                self._file.write(',\n')
                cpt += 1

        print(' %i properties successfully written.' % cpt)

    def _write_materials(self):
        print('Writing materials...')
        self._file.write('**\n** MATERIALS\n**\n')
        cpt = 0
        for material in self._model.materials:
            self._file.write('*Material, name=%s\n' % material.get_name())
            if material.get_density():
                self._file.write('*Density\n%1.1e,\n' % material.get_density())
            if isinstance(material, Isotropic):
                if isinstance(material.get_e(), float):
                    self._file.write('*Elastic\n%1.1f, %1.1f\n' % (material.get_e(), material.get_nu()))
                else:
                    table = material.get_e()
                    nr, nc = table.shape()
                    self._file.write('*Elastic, dependencies=%i\n' % (nc-2))
                    for xi in table:
                        e, nu = xi[:2]
                        f = ()
                        for fi in xi[2:]:
                            f += (fi,)
                        format = '%f, %f, %s' + ', %f'*(nc-2) + '\n'
                        values = (e,nu, '') + f
                        self._file.write(format % values)
                    self._file.write('*User Defined Field\n')
                cpt += 1
            elif isinstance(material, Orthotropic):
                self._file.write('*Elastic, type=LAMINA\n%1.1f, %1.1f, %1.1f, %1.1f, %1.1f, %1.1f, %1.1f\n' % (material.get_ex(), material.get_ey(), material.get_nuxy(), material.get_gxy(), material.get_gxz(), material.get_gyz(), material.get_theta()))
                cpt += 1
        print('%i materials successfully written.\n' % cpt)

    def _write_distribution_tables(self):
        if len(self._model.tables) > 0:
            print('Writing distribution tables...')
            self._file.write('**\n** DISTRIBUTION TABLES\n**\n')
            cpt = 0
            for table in self._model.tables:
                if len(table.opts) > 0:
                    self._file.write('*Distribution table, name=%s\n' % table['dtable'])
                    self._file.write('%s\n' % table['dtype'])
                    cpt += 1
            print('%i distribution tables successfully written.\n' % cpt)

    def _write_steps(self):
        print('Writing steps...')
        cpt = 0
        for step in self._model.steps:
            if isinstance(step, InitialStep):
                self._write_initial_step(step)
                cpt += 1
            elif isinstance(step, GeneralStaticStep):
                self._write_general_static_step(step)
                cpt += 1
            elif isinstance(step, BuckleStep):
                self._write_buckle_step(step)
                cpt += 1
        print('%i steps successfully written.\n' % cpt)

    def _write_initial_step(self, step):
        print('\tWriting initial boundary conditions...')
        self._write_bcs_and_loads(step)
        print('\tInitial boundary conditions successfully written.')

    def _write_general_static_step(self, step):
        print('\tWriting general static step %s...' % step.get_name())
        self._file.write('** ----------------------------------------------------------------\n')
        self._file.write('**\n** STEP: %s\n**\n' % step.get_name())
        line = '*Step, name=%s, nlgeom=%s'
        if step.get_non_linear():
            nl = 'YES'
        else:
            nl = 'NO'
        if step.get_nincmax():
            line += ', inc=%i\n'
            self._file.write(line % (step.get_name(), nl, step.get_nincmax()))
        else:
            line += '\n'
            self._file.write(line % (step.get_name(), nl))
        self._file.write('*Static\n')
        self._file.write('%f, %f, %e, %f\n' % (step.get_increment_init(), step.get_period(), step.get_increment_min(), step.get_increment_max()))
        self._write_bcs_and_loads(step)
        self._write_output_requests(step.outputs)
        self._file.write('*End Step\n')
        print('\tGeneral static step %s successfully written.' % step.get_name())

    def _write_bcs_and_loads(self, step):
        self._file.write('**\n** BOUNDARY CONDITIONS\n**\n')
        for bc in step.bcs:
            self._write_bc(bc)
        for load in step.loads:
            self._write_load(load)

    def _write_bc(self, bc):
        self._file.write('** Name: %s Type: Displacement/Rotation\n' % bc.get_name())
        if isinstance(bc, DisplacementBoundaryCondition):
            self._file.write('*Boundary, type=DISPLACEMENT\n')
        elif isinstance(bc, VelocityBoundaryCondition):
            self._file.write('*Boundary, type=VELOCITY\n')
        elif isinstance(bc, AccelerationBoundaryCondition):
            self._file.write('*Boundary, type=ACCELERATION\n')
        values = bc.get_values()
        target = bc.get_target()
        for i in range(len(values)):
            if values[i] != None:
                if isinstance(target, Node):
                    self._file.write('%i, %i, %i, %1.1f\n' % (target.get_label(), i+1, i+1, values[i]))
                else:
                    self._file.write('%s, %i, %i, %1.1f\n' % (target.get_name(), i+1, i+1, values[i]))

    def _write_loads(self, load):
        pass

    def _write_output_requests(self, output_requests):
        self._file.write(
'''**
** OUTPUT REQUESTS
**
*Restart, write, frequency=0
**
** FIELD OUTPUT: F-Output-1
**
*Output, field
*Node Output
U,
*Element Output, directions=YES
E, S
**
** HISTORY OUTPUT: H-Output-1
**
*Output, history, variable=PRESELECT
'''
)

        for output in output_requests:
            if isinstance(output, FieldOutputRequest):
                self._write_fo_request(output)
            elif isinstance(output, HistoryOutputRequest):
                self._write_no_request(output)

    def _write_fo_request(self, output):
        self._file.write('**\n** FIELD OUTPUT: %s\n**\n' % output.get_name())
        self._file.write('*Output, field\n')
        nvariables = output.get_nodal_variables()
        evariables = output.get_elemental_variables()
        if nvariables != None:
            self._file.write('*Node Output\n')
            s = ''
            for variable in nvariables:
                s += variable + ','
            s = s[:-1] + '\n'
            self._file.write(s)
        if evariables != None:
            self._file.write('*Element Output\n')
            s = ''
            for variable in evariables:
                s += variable + ','
            s = s[:-1] + '\n'
            self._file.write(s)

    def _write_no_request(self, output):
        self._file.write('**\n** HISTORY OUTPUT: %s\n**\n' % output.get_name())
        self._file.write('*Output, history\n')
        self._file.write('*Node Output, nset=%s\n' % output.get_group().get_name())
        s = ''
        for variable in output.get_nodal_variables():
            s += variable + ','
        s = s[:-1] + '\n'
        self._file.write(s)


