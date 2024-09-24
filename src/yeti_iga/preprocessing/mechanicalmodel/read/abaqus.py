# Copyright 2016-2019 Thibaut Hirschler
# Copyright 2020 Arnaud Duval

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

from .common import Parser

from ..mesh import Node, NodeParam, Beam, Tria, Quad, Penta, Hexa, UserElement, UserElementParam
from ..properties import ThinShell, Ply, CompositeLayup, CompositeThinShell, Solid,\
    userElementProps
from ..materials import Isotropic, Orthotropic, HighOrderElastic
from ..group import Group
from ..boundary_conditions import DisplacementBoundaryCondition,\
    VelocityBoundaryCondition,AccelerationBoundaryCondition
from ..loads import Load
#from preprocessing.mechanicalmodel.loads import Force
from ..steps import InitialStep, GeneralStaticStep, RiksStaticStep, BuckleStep
from ..common import TransformationMatrix, Table, Container

import numpy as np
from math import *
import pickle as pickle
import copy


class Part(object):
    def __init__(self, name=None):
        self.name  = name
        self.nodes      = Container()
        self.elements   = Container()
        self.properties = Container()
        self.groups     = Container()



class AbaqusParser(Parser):
    """
    Abaqus inp parser
    Limitations:
      - CS : only global coordinate system is supported for all entities
      - elements : only T3D2(1D bar), S3R (2D tria), S4R (2D quad), C3D8R (3D hexa) supported
      - properties : only shell, composite layups and solid sections supported
      - materials : only isotropic and orthotropic supported (density+elastic [type=LAMINA] cards)
      - BCs : only displacement, velocity and acceleration bcs
      - loads : distributed and concentrated loads
      - steps : only general static and perturbation buckle steps
    Supported cards:
      - Distributions :
          *DISTRIBUTION, name=, location=
          *DISTRIBUTION TABLE, name=
      - Materials :
          *MATERIAL, name=
          *ELASTIC, type=ISOTROPIC/LAMINA
          *DENSITY
      - Parts :
          *PART, name=
          *NODE
          *ELEMENT, type=S3R/S4R/C3D8R/T3D2
          *NSET, nset=, [internal], [generate]
          *ELSET, elset=, [internal], [generate]
          *SHELL SECTION, elset=, material=, offset=value/table
          *SHELL SECTION, elset=, composite, layup=
          *SOLID SECTION, elset=, material=
      - Assembly :
          *ASSEMBLY
          *INSTANCE, name=, part=
          *NODE
    To do : *PLASTIC
            *DAMAGE
            Interactions : tie, contact
            *SURFACE
            *ORIENTATION
            *CONTROLS
            *AMPLITUDE
            *SECTION CONTROLS
            *OUTPUT, HISTORY/FIELD
    """

    def __init__(self, model, filename):
        Parser.__init__(self, model=model, filename=filename)
        self._parts       = {}
        self._uel_params  = {}
        self._ndist       = 0
        self._nmats       = 0
        self._nparts      = 0
        self._ngroupspart = 0
        self._ngroups     = 0
        self._npropspart  = 0
        self._nprops      = 0
        self._ninstances  = 0
        self._nsteps      = 0
        self._nbcs        = 0
        self._nloads      = 0
        self._nb_nodes    = 0
        # Assembly counters
        self._nnodes      = 0
        self._nnodes_ass  = 0
        self._nelements   = 0
        self._nuelements  = 0
        self._ngroups_ass = 0
        self._node_labels          = {}
        self._assembly_node_labels = {}
        self._element_labels       = {}
        self._property_labels      = {}
        # errors
        self._error = []

    def parse(self):
        print("\nReading Abaqus inp file from %s..." % self._filename)
        self._read_nb_nodes()
        self._read_uel_parameters()
        self._read_model_name()
        self._read_parameter()
        self._read_distributions()
        self._read_materials()
        self._read_parts()
        self._read_assembly()
        self._read_initial_bcs()
        self._read_steps()
        self._file.close()

        if len(self._error)==0:
            print("\nAbaqus inp file %s successfully read." % self._filename)
        else:
            print("\nAbaqus inp file %s read with the following errors" % self._filename)
            for err in self._error:
                print(' - '+ err)
        return None

    def _read_nb_nodes(self):
        # read the number of nodes by element in case of user elements
        print("\nReading number of nodes by element...")
        while 1:
            line = self._file.readline()
            if line.upper().startswith('*USER ELEMENT'):
                self._nb_nodes = int(line.upper().split('NODES=')[1].split(',')[0])
                #print('Number of nodes by element : ', self._nb_nodes)
                break
            elif len(line)==0:
                print("Number of nodes not available")
                self._file.seek(0)
                break

    def _read_uel_parameters(self):
        print("\nReading parameters defining the user elements...")
        self._file.seek(0)
        while 1:
            line = self._file.readline()
            if line.upper().startswith('*USER ELEMENT'):
                self._nuelements += 1
                try:
                    name =line.upper().split('TYPE=')[1].split(',')[0].strip()
                except:
                    self._error.append('error: type of user element is not defined')
                    print(self._error[-1])
                    break
                try:
                    nnode=int(line.upper().split('NODES=')[1].split(',')[0])
                except:
                    self._error.append('error: node of user element is not defined')
                    print(self._error[-1])
                    break

                uel_param = UserElementParam(label=name,name=name,nnode=nnode)

                options = line.rstrip().upper().split(',')
                for opt in options:
                    if 'COORDINATES' in opt:
                        uel_param.set_coordinates(int(opt.split('=')[1]))
                    elif 'VARIABLES' in opt:
                        uel_param.set_variables(int(opt.split('=')[1]))
                    elif 'I PROPERTIES' in opt:
                        uel_param.set_iproperties(int(opt.split('=')[1]))
                    elif 'INTEGRATION' in opt:
                        uel_param.set_integration(int(opt.split('=')[1]))
                    elif 'TENSOR' in opt:
                        uel_param.set_tensor(opt.split('=')[1].strip())
                self._uel_params[name] = uel_param
                self._model.elements.add(uel_param)

            if not line:
                break


    def _read_model_name(self):
        print("\nReading model name...")
        while 1:
           line = self._file.readline()
           if 'Model name' in line:
                self._model.info['name'] = line.split(':')[2].strip()
                print("Model name : %s..." % self._model.info['name'])
                break
           elif len(line)==0:
                print("No model name.")
                self._file.seek(0)
                break

    def _updateparameters(self):
        if not hasattr(self,'_design_parameters'):
            return None
        if len(self._design_parameters)>0:
            for key,value in self._design_parameters.items():
                paramname = key.translate({ord('<'): None,ord('>'): None})
                exec("%s=%s"%(paramname,value))
            exec(self._parameters_def)
            for key in self._parameters.keys():
                self._parameters[key] = eval(key)

    def _read_parameter(self):
        print('\nReading parameters...')
        self._read_designparameter()
        self._parameters = {}
        #if len(self._design_parameters)>0:
        self._parameters_def = ''
        while 1:
            line = self._file.readline()
            if line.startswith('**'):
                continue
            elif line.upper().startswith('*PARAMETER'):
                self._assignValueToParameters()
            elif not line:
                break
        self._file.seek(0)
        if len(self._parameters)>0:
            print('\tRead %i parameters'%len(self._parameters))
            for key, value in self._parameters.items():
                print('\t',key, ' : ', value)

    def _read_designparameter(self):
        self._design_parameters = {}
        while 1:
            line = self._file.readline()
            if line.upper().startswith('*DESIGN PARAMETER'):
                line = self._file.readline()
                while not line.startswith('*'):
                    dsgparam = line.split(',')
                    for d in dsgparam:
                        paramname = d.strip()
                        if len(paramname)>0:
                            self._design_parameters[paramname] = None
                    line = self._file.readline()
            elif not line:
                break
        self._file.seek(0)
        if len(self._design_parameters)>0:
            print('\tFound %i design parameter(s)'%len(self._design_parameters))
            for key in self._design_parameters.keys():
                print('\t',key)

    def _assignValueToParameters(self):
        line = self._file.readline()
        savedef = True if len(self._design_parameters) else False
        while not line.startswith('*'):
            exec(line)
            paramName = line.split('=')[0].strip()
            if savedef:
                if not ('<'+paramName+'>' in self._design_parameters.keys()):
                    self._parameters_def += line
                else:
                    self._design_parameters['<'+paramName+'>'] = eval(paramName)
            self._parameters[paramName] = eval(paramName)
            line = self._file.readline()

    def _read_distributions(self):
        print('\nReading distributions...')
        while 1:
            line = self._file.readline()
            if line.startswith('**'):
                continue
            elif line.upper().startswith('*DISTRIBUTION,'):
                self._read_distribution(line)
            elif not line:
                break
        self._file.seek(0)
        while 1:
            line = self._file.readline()
            if line.startswith('**'):
                continue
            elif line.upper().startswith('*DISTRIBUTION TABLE'):
                self._read_distribution_table(line)
            elif not line:
                break
        print('%i distributions read.' % len(self._model.tables))
        self._file.seek(0)

    def _read_distribution(self, line):
        self._ndist += 1
        distribtable = None
        for option in line.split(',')[1:]:
            option_name, value = option.split('=')
            if option_name.strip().upper() == 'NAME':
                name = value.strip()
            elif option_name.strip().upper() == 'LOCATION':
                loc = value.strip()
            elif option_name.strip().upper() == 'TABLE':
                distribtable = value.strip()
        print('\tReading distribution %s...' % name)
        table = Table(name=name, label=self._ndist,
                      opts={'location': loc,
                            'dtable': distribtable,
                            'dtype': 'unused'})
        while 1:
            pos = self._file.tell()
            line = self._file.readline()
            if line.upper().startswith('**'):
                continue
            elif line.startswith('*'):
                # End of reading table
                print('\t%i values read.' % len(table))
                print('\tDistribution %s successfully read.' % name)
                self._model.tables.add(table)
                self._file.seek(pos)
                return
            else:
                datastr = line.strip().split(',')
                if datastr[0]=='':
                    table.set_defaultvalues(datastr[1:])
                else:
                    data = [int(datastr[0])]
                    for v in datastr[1:]:
                        data.append(float(v))
                    table.append([data])

    def _read_distribution_table(self, line):
        name = line.split('=')[1].strip()
        linetype = self._file.readline().split(',')
        tabletype = [t.strip() for t in linetype]
        for table in self._model.tables:
            if table.opts['dtable'].upper() == name.upper():
                table.opts['dtype'] = tabletype

    def _read_materials(self):
        print('\nReading materials...')
        while 1:
            line = self._file.readline()
            if line.startswith('**'):
                continue
            elif line.upper().startswith('*MATERIAL'):
                self._read_material(line)
            elif not line:
                break
        print('%i materials read.' % len(self._model.materials))
        self._file.seek(0)

    def _read_material(self, line):
        self._nmats += 1
        name = line.split(',')[1].split('=')[1].strip()
        print('\tReading material %s...' % name)
        rho  = 0.
        while 1:
            pos  = self._file.tell()
            line = self._file.readline()
            if line.upper().startswith('*ELASTIC'):
                if 'TYPE' in line.upper():
                    mat_type = line[line.index('TYPE'):].split('=')[1].strip().upper()
                    if mat_type == 'LAMINA':
                        mat = Orthotropic(label=self._nmats, name=name, density=rho)
                        ex, ey, nuxy, gxy, gxz, gyz, theta = self._file.readline().split(',')
                        mat.set_ex(float(ex))
                        mat.set_ey(float(ex))
                        mat.set_nuxy(float(nuxy))
                        mat.set_gxy(float(gxy))
                        mat.set_gxz(float(gxz))
                        mat.set_gyz(float(gyz))
                        mat.set_theta(float(theta))
                        print('\t\tType orthotropic.')
                        print('\t\tEx = %1.1f' % float(ex))
                        print('\t\tEy = %1.1f' % float(ey))
                        print('\t\tNuxy = %1.1f' % float(nuxy))
                        print('\t\tGxy = %1.1f' % float(gxy))
                        print('\t\tGxz = %1.1f' % float(gxz))
                        print('\t\tGyz = %1.1f' % float(gyz))
                        print('\t\tTheta = %1.1f' % float(theta))
                else:
                    mat = Isotropic(label=self._nmats, name=name, density=rho)
                    linedata = self._file.readline().strip().split(',')
                    if len(linedata)==2:
                        e, nu = linedata
                        mat.set_e(float(e))
                        mat.set_nu(float(nu))
                    else:
                        distributionname = linedata[0]
                        if distributionname in self._model.tables.get_names():
                            distribution = self._model.tables[distributionname]
                            if (distribution.opts['dtype'][0].upper() == 'MODULUS' \
                                and distribution.opts['dtype'][1].upper() == 'RATIO'):
                                values = distribution.get_defaultvalues()
                                e  = values[0]
                                nu = values[1]
                                mat.set_e(float(e))
                                mat.set_nu(float(nu))
                                mat.link2distribution(distributionname)
                    print('\t\tType isotropic.')
                    print('\t\tElastic modulus = %1.1f' % mat.get_e())
                    print('\t\tPoisson = %1.1f' % mat.get_nu())
            elif line.upper().startswith('*HIGHORDERELASTIC'):
                mat = HighOrderElastic(label=self._nmats, name=name)
                params = self._file.readline().split(',')
                fparams = [float(ele) for ele in params]
                mat.set_lambda(fparams[0])
                mat.set_mu(fparams[1])
                mat.set_a(fparams[2:7])
                mat.set_b(fparams[7:15])
                mat.set_c(fparams[15:])
                print('\t\tType isotropic 2nd order.')
                print('\t\tlambda = %1.1f' % fparams[0])
                print('\t\tmu = %1.1f' % fparams[1])
                print('\t\ta = ', fparams[2:7])
                print('\t\tb = ', fparams[7:15])
                print('\t\tc = ', fparams[15:])

            elif line.upper().startswith('*DENSITY'):
                rho = float(self._file.readline().split(',')[0])
                try:
                    mat.set_density(float(rho))
                    print('\t\tDensity = %1.1f' % float(rho))
                except:
                    pass
            elif line.startswith('**'):
                continue
            elif line.upper().startswith('*MATERIAL'):
                break
            elif line == '':
                break
        print('\tMaterial %s successfully read.' % name)
        self._file.seek(pos)
        self._model.materials.add(mat)

    def _read_parts(self):
        print('\nReading parts...')
        self._model.info['parts'] = []
        while 1:
            line = self._file.readline()
            if line.startswith('**'):
                continue
            elif line.upper().startswith('*PART'):
                self._read_part(line)
                self._nparts += 1
            elif not line:
                break
            else:
                continue
        print('%i parts read.' % self._nparts)
        print('\nSaving part objects...')
        for part in self._parts:
            pickle.dump(self._parts[part], open( "%s.save" % part, "wb" ) )
        print('Parts successfully saved.')
        del self._parts
        self._nprops = 0
        self._file.seek(0)

    def _read_part(self, line):
        self._ngroupspart = 0
        self._npropspart = 0
        self._partname = line.split(',')[1].split('=')[1].strip()
        print('\tReading part %s...' % self._partname)
        self._model.info['parts'].append(self._partname)
        #self._parts[self._partname] = mm.Part(name=self._partname)
        self._parts[self._partname] = Part(name=self._partname)
        while 1:
            line = self._file.readline()
            if line.upper().startswith('** SECTION'):
                self._sectionname = line.split(':')[1].strip()
            elif line.startswith('**'):
                continue
            elif '*NODE' in line.upper():
                self._read_nodes(mode='part')
            elif '*ELEMENT' in line.upper():
                self._read_elements(line)
            elif '*NSET' in line.upper():
                self._read_set(line, type='node', mode='part')
            elif '*ELSET' in line.upper():
                self._read_set(line, type='element', mode='part')
            elif '*UEL' in line.upper():
                self._read_uel_property(line)
            elif '*SHELL SECTION' in line.upper() and not 'COMPOSITE' in line.upper():
                self._read_shell_section(line)
            elif '*SHELL SECTION' in line.upper() and 'COMPOSITE' in line.upper():
                self._read_composite_shell_section(line)
            elif '*SOLID SECTION' in line.upper():
                self._read_solid_section(line)
            elif '*END PART' in line.upper():
                print('\t\t%i groups read.' % self._ngroupspart)
                print('\t\t%i properties read.' % self._npropspart)
                print('\tPart %s successfully read.' % self._partname)
                break

    def _read_nodes(self, mode):
        while 1:
            line = self._file.readline()
            if line.startswith('**'):
                continue
            elif '*' in line:
                self._file.seek(pos)
                break
            else:
                pos = self._file.tell()
                if len(line.split(','))==4:
                    label, x, y, z = line.split(',')
                elif len(line.split(','))==3:
                    label, x, y = line.split(',')
                    z = '0'

                coordsstr = [x.strip(),y.strip(),z.strip()]
                coordsflt = []
                flagparam = False
                for cstr in coordsstr:
                    try:
                        cflt = float(cstr)
                    except:
                        param = cstr.translate({ord('<'): None,ord('>'): None})
                        cflt = float(self._parameters[param])
                        flagparam = True
                    coordsflt.append(cflt)
                if not flagparam:
                    node = Node(label=int(label), coordinates=tuple(coordsflt))
                else:
                    node = NodeParam(label=int(label), coordinates=tuple(coordsflt), parametrization=coordsstr)
                #node = Node(label=int(label), coordinates=(float(x),float(y),float(z)))

                if mode == 'part':
                    self._parts[self._partname].nodes.add(node)

                elif mode == 'assembly':
                    self._nnodes += 1
                    self._nnodes_ass += 1
                    node.set_misc({'assembly':True})
                    node.set_label(self._nnodes)
                    self._assembly_node_labels[int(label)] = self._nnodes
                    self._model.nodes.add(node)
        if mode == 'part':
            print('\t\t%i nodes read.' % len(self._parts[self._partname].nodes))


    def _read_elements(self, line, mode='part'):
        type = line.split(',')[1].split('=')[1].strip().upper()
        while 1:
            line = self._file.readline()
            if '*' in line:
                self._file.seek(pos)
                break
            elif type == 'T3D2':
                label, n1, n2 = line.split(',')
                element = Beam(label=int(label),
                               nodes=[self._parts[self._partname].nodes[int(n1)],
                                      self._parts[self._partname].nodes[int(n2)]])
            elif type == 'S3R':
                label, n1, n2, n3 = line.split(',')
                element = Tria(label=int(label),
                               nodes=[self._parts[self._partname].nodes[int(n1)],
                                      self._parts[self._partname].nodes[int(n2)],
                                      self._parts[self._partname].nodes[int(n3)]])
            elif type == 'S4R':
                label, n1, n2, n3, n4 = line.split(',')
                element = Quad(label=int(label),
                               nodes=[self._parts[self._partname].nodes[int(n1)],
                                      self._parts[self._partname].nodes[int(n2)],
                                      self._parts[self._partname].nodes[int(n3)],
                                      self._parts[self._partname].nodes[int(n4)]])
                nodes=[self._parts[self._partname].nodes[int(n1)],
                                      self._parts[self._partname].nodes[int(n2)],
                                      self._parts[self._partname].nodes[int(n3)],
                                      self._parts[self._partname].nodes[int(n4)]]
            elif type == 'C3D8R':
                label, n1, n2, n3, n4, n5, n6, n7, n8 = line.split(',')
                element = Hexa(label=int(label),
                               nodes=[self._parts[self._partname].nodes[int(n1)],
                                      self._parts[self._partname].nodes[int(n2)],
                                      self._parts[self._partname].nodes[int(n3)],
                                      self._parts[self._partname].nodes[int(n4)],
                                      self._parts[self._partname].nodes[int(n5)],
                                      self._parts[self._partname].nodes[int(n6)],
                                      self._parts[self._partname].nodes[int(n7)],
                                      self._parts[self._partname].nodes[int(n8)]])
            elif type[0]=='U':
                label = line.split(',')[0]
                line2 = line
                while len(line2.rstrip().split(',')[-1].split())==0:
                    line2 = self._file.readline()
                    line  = line + line2
                nodes=[]
                for i_str in line.split(',')[1:]:
                    nodes.append(self._parts[self._partname].nodes[int(i_str)])
                element = UserElement(label=int(label),nodes=nodes,type=type.split()[0])
            else:
                print('Element type %s not yet supported.' % type)
                break

            if mode == 'part':
                self._parts[self._partname].elements.add(element)

            elif mode == 'assembly':
                self._error.append('error: element cannot be defined in assembly')
                print(self._error[-1])

            pos = self._file.tell()
            #self._model.elements.add(element)


        print('\t\t%i elements read.' % len(self._parts[self._partname].elements))


    def _read_set(self, line, type, mode):
        line_split = line.split(',')
        groupname = line_split[1].split('=')[1].strip()
        generate = False
        misc = {'assembly':False, 'internal': False, 'from_instance':False}
        if len(line_split) > 2:
            for k in range(2,len(line_split)):
                if line_split[k].strip().upper() == 'GENERATE':
                    generate = True
                if line_split[k].strip().upper() == 'INTERNAL':
                    misc['internal'] = True
                try:
                    keyword, value = line_split[k].split('=')
                    if keyword.strip().upper() == 'INSTANCE':
                        misc['from_instance'] = True
                except:
                    pass
        if mode == 'part':
            groups = self._parts[self._partname].groups
            if groupname in groups:
                group = groups[groupname]
            else:
                self._ngroupspart += 1
                group = Group(label=self._ngroupspart, name=groupname, type=type, misc=misc)
        elif mode == 'assembly':
            groups = self._model.groups
            if groupname in groups:
                group = groups[groupname]
            else:
                self._ngroups += 1
                group = Group(label=self._ngroups, name=groupname, type=type, misc=misc)
        if generate:
            line = self._file.readline()
            start, end, step = line.split(',')
            for k in range(int(start),int(end)+1,int(step)):

                if type == 'node':
                    if mode == 'part':
                        group.nodes.add(self._parts[self._partname].nodes[k])
                    elif mode == 'assembly':
                        if group.get_misc()['from_instance']:
                            group.nodes.add(self._model.nodes[k])
                        else:
                            group.nodes.add(self._model.nodes[self._assembly_node_labels[k]])
                elif type == 'element':
                    if mode == 'part':
                        group.elements.add(self._parts[self._partname].elements[k])
                    elif mode == 'assembly':
                        group.elements.add(self._model.elements[k])
        else:
            while 1:

                line = self._file.readline()
                if line.startswith('**'):
                    continue
                elif '*' in line:
                    self._file.seek(pos)
                    break
                for k in line.split(','):

                    pos = self._file.tell()
                    if k.strip():
                        if type == 'node':
                            if mode == 'part':
                                group.nodes.add(self._parts[self._partname].nodes[int(k)])

                            elif mode == 'assembly':
                                if group.get_misc()['from_instance']:
                                    group.nodes.add(self._model.nodes[int(k)])
                                else:
                                    group.nodes.add(self._model.nodes[
                                        self._assembly_node_labels[int(k)]])

                        elif type == 'element':
                            if mode == 'part':
                                group.elements.add(self._parts[self._partname].elements[int(k)])
                            elif mode == 'assembly':
                                group.elements.add(self._model.elements[int(k)])


        if mode == 'part':
            self._parts[self._partname].groups.add(group)
        elif mode == 'assembly':
            group.get_misc()['assembly']=True
            self._model.groups.add(group)

    def _read_shell_section(self, line):
        self._npropspart += 1
        self._nprops += 1
        offset = 0.
        for option in line.split(',')[1:]:
            option_name, value = option.split('=')
            if option_name.strip().upper() == 'ELSET':
                elset = value.strip()
            elif option_name.strip().upper() == 'MATERIAL':
                material = value.strip()
            elif option_name.strip().upper() == 'OFFSET':
                offset = value.strip()
        thickness, nint = self._file.readline().split(',')
        try:
            offset = float(offset)
        except ValueError:
            offset = self._model.tables[offset]
        property = ThinShell(label=self._nprops, name=self._partname+'.'+self._sectionname,
                             elset=self._partname+'.'+elset, thickness=float(thickness),
                             nint=int(nint), offset=offset)
        #self._model.properties.add(property)
        self._parts[self._partname].properties.add(property)
        for element in self._parts[self._partname].groups[elset].elements:
            element.set_material(self._model.materials[material])
            element.set_property(self._parts[self._partname].properties[property.get_name()])

    def _read_uel_property(self,line):
        self._npropspart += 1
        self._nprops += 1
        for option in line.split(',')[1:]:
            option_name, value = option.split('=')
            if option_name.strip().upper() == 'ELSET':
                elset = value.strip()
            elif option_name.strip().upper() == 'MATERIAL':
                material = value.strip()
        #
        line = self._file.readline()
        props= np.array(line.split(',')[:], dtype=float)
        property = userElementProps(
            label=self._nprops,
            name='UEL-props'+self._partname+'.'+str(self._nprops),
            elset=self._partname+'.'+elset, props=props[:],
            material=material)
        self._parts[self._partname].properties.add(property)
        #
        for element in self._parts[self._partname].groups[elset].elements:
            element.set_material(self._model.materials[material])
            element.set_property(self._parts[self._partname].properties[property.get_name()])


    def _read_composite_shell_section(self, line):
        self._npropspart += 1
        self._nprops += 1
        for option in line.split(',')[1:]:
            try:
                option_name, value = option.split('=')
            except ValueError:
                # Composite option : no from like name = value
                continue
            if option_name.strip().upper()   == 'ELSET':
                elset = value.strip()
            elif option_name.strip().upper() == 'LAYUP':
                layupname = value.strip()
        layup = CompositeLayup(name=layupname)
        while 1:
            pos = self._file.tell()
            line = self._file.readline()
            if '*' in line:
                break
            t, nint, mat, theta, name = line.split(',')
            layup.append(Ply(name=name.strip(), thickness=float(t), orientation=float(theta),
                             material=mat.strip(), nint=int(nint)))

        self._file.seek(pos)
        property = CompositeThinShell(label=self._nprops, name='Section-composite-'+self._partname+'.'+str(self._nprops), elset=self._partname+'.'+elset, layup=layup)
        #self._model.properties.add(property)
        self._parts[self._partname].properties.add(property)
        for element in self._parts[self._partname].groups[elset].elements:
            element.set_material(self._model.materials[mat.strip()])
            element.set_property(self._parts[self._partname].properties[property.get_name()])

    def _read_solid_section(self, line):
        self._npropspart += 1
        self._nprops += 1
        line_split = line.split(',')
        for option in line_split[1:]:
            option_name, value = option.split('=')
            if option_name.strip().upper() == 'ELSET':
                elset = value.strip()
            elif option_name.strip().upper() == 'MATERIAL':
                material = value.strip()
        property = Solid(label=self._nprops, name=self._partname+'.'+self._sectionname,
                         elset=self._partname+'.'+elset)
        #self._model.properties.add(property)
        self._parts[self._partname].properties.add(property)
        for element in self._parts[self._partname].groups[elset].elements:
            element.set_material(self._model.materials[material])
            element.set_property(self._parts[self._partname].properties[property.get_name()])


    def _read_assembly(self):
        print('\nReading assembly...')
        while 1:
            line = self._file.readline()
            if line.upper().startswith('*ASSEMBLY'):
                break
        self._model.info['instances'] = {}
        while 1:
            line = self._file.readline()
            if line.startswith('**'):
                continue
            elif line.upper().startswith('*NODE'):
                self._read_nodes(mode='assembly')
            #elif line.upper().startswith('*ELEMENT'):
                #self._read_elements(line,mode='assembly')
            elif line.upper().startswith('*INSTANCE'):
                self._read_instance(line)
            elif '*NSET' in line.upper():
                self._read_set(line, type='node', mode='assembly')
            elif '*ELSET' in line.upper():
                self._read_set(line, type='element', mode='assembly')
            elif line.upper().startswith('*END ASSEMBLY'):
                print('%i nodes read.' % self._nnodes_ass)
                print('%i instances created.' % self._ninstances)
                break
            else:
                continue
        self._file.seek(0)

    def _read_instance(self, line):
        self._ninstances += 1
        for option in line.split(',')[1:]:
            option_name, value = option.split('=')
            if   option_name.strip().upper() == 'NAME':
                instancename = value.strip()
            elif option_name.strip().upper() == 'PART':
                part = value.strip()
        if part not in self._model.info['instances']:
            self._model.info['instances'][part] = [instancename]
        else:
            self._model.info['instances'][part].append(instancename)
        pos = self._file.tell()
        line = self._file.readline()
        if line.upper().startswith('*END INSTANCE'):
            m = TransformationMatrix()
            self._create_instance(part, instancename, m)
            return
        else:
            # First line (if present) : translation matrix
            m = TransformationMatrix()
            tx, ty, tz = line.split(',')
            m.translate([float(tx),float(ty),float(tz)])
        pos = self._file.tell()
        line = self._file.readline()
        if line.upper().startswith('*END INSTANCE'):
            self._create_instance(part, instancename, m)
            return
        else:
            # Second line (if present) : rotation matrix
            ax, ay, az, bx, by, bz, theta = line.split(',')
            ab = np.array([float(bx)-float(ax), float(by)-float(ay), float(bz)-float(az)])
            abx, aby, abz = ab/np.linalg.norm(ab)
            theta = np.radians(float(theta))
            c = np.cos(theta)
            s = np.sin(theta)
            m.get()[0][0] = abx**2 + (1-abx**2)*c
            m.get()[0][1] = abx*aby*(1-c) - abz*s
            m.get()[0][2] = abx*abz*(1-c) + aby*s
            m.get()[1][0] = abx*aby*(1-c) + abz*s
            m.get()[1][1] = aby**2 + (1-aby**2)*c
            m.get()[1][2] = aby*abz*(1-c) - abx*s
            m.get()[2][0] = abx*abz*(1-c) - aby*s
            m.get()[2][1] = aby*abz*(1-c) + abx*s
            m.get()[2][2] = abz**2 + (1-abz**2)*c
            self._create_instance(part, instancename, m)

    def _create_instance(self, partname, instancename, matrix):
        print('Creating instance %s of part %s...' % (instancename, partname))
        # Load part
        part = pickle.load( open("%s.save" % partname, "rb" ))
        # Create nodes
        flagshapeparam = False
        for node in part.nodes:
            self._nnodes += 1
            self._node_labels[(instancename,node.get_label())] = self._nnodes
            coords = np.dot(matrix.get(), np.append(node.get_coordinates(),1))[:3]
            if node.__class__ == NodeParam:
                parametrization = node.get_parametrization()
                self._model.nodes.add(NodeParam(label=self._nnodes, coordinates=coords, parametrization=parametrization))
                flagshapeparam = True
                #self._model.nodes.update(self._parameters)
            else:
                self._model.nodes.add(Node(label=self._nnodes, coordinates=coords))
        if flagshapeparam:
            self._model._design_parameters = self._design_parameters
            self._model._parameters_def = self._parameters_def
            self._model._parameters = self._parameters

        # Create elements
        for element in part.elements:
            self._nelements += 1
            self._element_labels[(instancename,element.get_label())] = self._nelements
            label = self._nelements
            if isinstance(element, Beam):
                n1,n2 = element.get_nodes()
                self._model.elements.add(
                    Beam(label=label,
                         nodes=[
                             self._model.nodes[self._node_labels[(instancename,n1.get_label())]],
                             self._model.nodes[self._node_labels[(instancename,n2.get_label())]]]))
            elif isinstance(element, Tria):
                n1,n2,n3 = element.get_nodes()
                self._model.elements.add(
                    Tria(label=label,
                         nodes=[
                             self._model.nodes[self._node_labels[(instancename,n1.get_label())]],
                             self._model.nodes[self._node_labels[(instancename,n2.get_label())]],
                             self._model.nodes[self._node_labels[(instancename,n3.get_label())]]]))
            elif isinstance(element, Quad):
                n1,n2,n3,n4 = element.get_nodes()
                self._model.elements.add(
                    Quad(label=label,
                         nodes=[
                             self._model.nodes[self._node_labels[(instancename,n1.get_label())]],
                             self._model.nodes[self._node_labels[(instancename,n2.get_label())]],
                             self._model.nodes[self._node_labels[(instancename,n3.get_label())]],
                             self._model.nodes[self._node_labels[(instancename,n4.get_label())]]]))
            elif isinstance(element, Hexa):
                n1,n2,n3,n4,n5,n6,n7,n8 = element.get_nodes()
                self._model.elements.add(
                    Hexa(label=label,
                         nodes=[
                             self._model.nodes[self._node_labels[(instancename,n1.get_label())]],
                             self._model.nodes[self._node_labels[(instancename,n2.get_label())]],
                             self._model.nodes[self._node_labels[(instancename,n3.get_label())]],
                             self._model.nodes[self._node_labels[(instancename,n4.get_label())]],
                             self._model.nodes[self._node_labels[(instancename,n5.get_label())]],
                             self._model.nodes[self._node_labels[(instancename,n6.get_label())]],
                             self._model.nodes[self._node_labels[(instancename,n7.get_label())]],
                             self._model.nodes[self._node_labels[(instancename,n8.get_label())]]]))
            elif isinstance(element,UserElement):
                elt_type = element.get_type()
                tmp = []
                queue = element.get_nodes()
                for ni in queue:
                    tmp.append(
                        self._model.nodes[self._node_labels[(instancename,ni.get_label())]])
                self._model.elements.add(
                    UserElement(label=label,nodes=tmp,type=elt_type))

            # Set element property and material
            #self._model.elements[label].set_property(
            #    self._model.properties[element.get_property().get_label()])
            self._model.elements[label].set_material(
                self._model.materials[element.get_material().get_label()])

        # Create sets
        # Group name = "instancename-groupname"
        for group in part.groups:
            self._ngroups += 1
            groupname = instancename+'-'+group.get_name()
            g = Group(label=self._ngroups, name=groupname, type=group.get_type(),
                      misc=group.get_misc())
            for node in group.nodes:
                g.nodes.add(self._model.nodes[self._node_labels[(instancename,node.get_label())]])
            for element in group.elements:
                g.elements.add(
                    self._model.elements[self._element_labels[(instancename,element.get_label())]])
            self._model.groups.add(g)

        # Properties : redefine elset name with new group names
        new_properties = []
        for property in part.properties:
            if partname in property.get_name():
                self._nprops += 1
                new_property = copy.deepcopy(property)
                new_property.set_label(self._nprops)
                new_property.set_name(
                    property.get_name().replace(partname,'',1).replace('.','')+'-'+instancename)
                new_property.set_elset(instancename+'-'+property.get_elset().split('.')[1])
                new_properties.append(new_property)
                self._property_labels[property.get_label()] = new_property.get_label()

        for new_property in new_properties:
            self._model.properties.add(new_property)

        # Set element property
        for element in part.elements:
            element.set_property(
                self._model.properties[self._property_labels[element.get_property().get_label()]])

        print('Successfully created.')

    def _read_initial_bcs(self):
        print('\nReading initial boundary conditions...')
        self._nsteps += 1
        self._step = InitialStep(label=self._nsteps, name='Initial')
        self._model.steps.add(self._step)
        print('\tInitial step created.')
        while 1:
            line = self._file.readline()
            if  line.upper().startswith('*BOUNDARY'):
                self._read_bc(line)
            elif line.upper().startswith('*STEP'):
                break
        print('Initial boundary conditions read.')
        self._file.seek(0)

    def _read_steps(self):
        print('\nReading steps...')
        while 1:
            line = self._file.readline()
            if line.startswith('**'):
                continue
            elif line.upper().startswith('*STEP'):
                self._read_step(line)
            elif not line:
                break
        print('%i steps read.' % len(self._model.steps))
        self._file.seek(0)

    def _read_step(self, line):
        self._nsteps += 1
        nlgeom  = 'NO'
        nincmax = None
        name    ='No name'
        for option in line.split(',')[1:]:
            option_name, value = option.split('=')
            if option_name.strip().upper()   == 'NAME':
                name    = value.strip()
            elif option_name.strip().upper() == 'NLGEOM':
                nlgeom  = value.strip()
            elif option_name.strip().upper() == 'INC':
                nincmax = float(value.strip())

        print('\tReading step ...' , name)
        while 1:
            line = self._file.readline()
            if line.startswith('**'):
                continue
            elif line.upper().startswith('*STATIC'):
                line = self._file.readline()
                if (len(line.split(',')) != 4 or line.startswith('**')):
                    inc_init, period, inc_min, inc_max = 1., 1., 1e-05, 1 #default values
                else :
                    inc_init, period, inc_min, inc_max = line.split(',')
                self._step = GeneralStaticStep(
                    label=self._nsteps, name=name, period=float(period), non_linear=nlgeom,
                    nincmax=nincmax, increment_init=float(inc_init), increment_min=float(inc_min),
                    increment_max=float(inc_max))
                self._model.steps.add(self._step)
            elif line.upper().startswith('*BUCKLE'):
                line = self._file.readline()
                neigs, maxeigs, nvectors, maxiter = line.split(',')
                if not maxeigs:
                    maxeigs = None
                else:
                    maxeigs = int(maxeigs)
                self._step = BuckleStep(label=self._nsteps, name=name, neigenvalues=int(neigs),
                                        maxeigenvalues=maxeigs, nvectors=int(nvectors),
                                        maxiter=int(maxiter))
                self._model.steps.add(self._step)
            elif line.upper().startswith('*BOUNDARY'):
                self._read_bc(line)
            elif line.upper().startswith('*DLOAD'):
                self._read_dload(line)
            elif line.upper().startswith('*CLOAD'):
                self._read_cload(line)
            if line.upper().startswith('*END STEP'):
                print('\tStep %s successfully read with %i bcs and %i loads.' \
                    % (name, len(self._model.steps[name].bcs), len(self._model.steps[name].loads)))
                break


    def _read_bc(self, line):
        name = 'bc with no name'
        print('\t\tReading boundary condition...')
        if line.upper().startswith('*BOUNDARY'):
            line_split = line.split(',')
            try:
                bc_type = line.rstrip().upper().split('TYPE=')[1]
            except:
                bc_type = 'DISPLACEMENT'

            while 1:
                pos  = self._file.tell()
                line = self._file.readline()
                if line.startswith('*'):
                    self._file.seek(pos)
                    break

                self._nbcs += 1
                bc_name = name+str(self._nbcs)
                if   bc_type == 'DISPLACEMENT':
                    bc = DisplacementBoundaryCondition(label=self._nbcs, name=bc_name)
                elif bc_type == 'VELOCITY':
                    bc = VelocityBoundaryCondition(    label=self._nbcs, name=bc_name)
                elif bc_type == 'ACCELERATION':
                    bc = AccelerationBoundaryCondition(label=self._nbcs, name=bc_name)
                else:
                    self._error.append('error: unknown type of bc: ' + bc_type)
                    print(self._error[-1])
                    break
                line_split = line.split(',')
                if len(line_split) == 2:
                    target, keyword = line_split
                    keyword = keyword.strip().upper()
                    value = 0
                    if   keyword == 'XSYMM':
                        ddls = [1,5,6]
                    elif keyword == 'YSYMM':
                        ddls = [2,4,6]
                    elif keyword == 'ZSYMM':
                        ddls = [3,4,5]
                    if keyword == 'XASYMM':
                        ddls = [2,3,4]
                    elif keyword == 'YASYMM':
                        ddls = [1,3,5]
                    elif keyword == 'ZASYMM':
                        ddls = [1,2,6]
                    elif keyword == 'PINNED':
                        ddls = [1,2,3]
                    elif keyword == 'ENCASTRE':
                        ddls = [1,2,3,4,5,6]
                else:
                    if len(line_split) == 3:
                        target, first_ddl, last_ddl = line_split
                        value = 0
                    elif len(line_split) == 4:
                        target, first_ddl, last_ddl, value = line_split

                    first_ddl = int(first_ddl)
                    last_ddl  = int(last_ddl)
                    value     = float(value)
                    ddls      = list(range(first_ddl, last_ddl+1))
                    target = target.split('.')[0] + '-' + target.split('.')[1]

                try:
                    # Case of node label
                    target = self._model.nodes[int(target)]
                except:
                    # Case of node set
                    target = self._model.groups[target.strip()]
                bc.set_target(target)

                for ddl in ddls:
                    if ddl == 1:
                        bc.set_x(value)
                    elif ddl == 2:
                        bc.set_y(value)
                    elif ddl == 3:
                        bc.set_z(value)
                    elif ddl == 4:
                        bc.set_rx(value)
                    elif ddl == 5:
                        bc.set_ry(value)
                    elif ddl == 6:
                        bc.set_rz(value)

                self._model.steps[self._step.get_name()].bcs.add(bc)
                print('\t\t - Boundary condition (%s) successfully added.' % bc.get_name())



    def _read_dload(self,line):
        print('\t\tReading Dload ...')
        while 1:
            pos  = self._file.tell()
            line = self._file.readline()
            if line.startswith('*'):
                self._file.seek(pos)
                break

            self._nloads += 1
            DLOAD = Load(label=self._nloads,name='no name', type='dload')

            line_split = line.rstrip().split(',')
            target=line_split[0]
            try:
                # Case of element label
                target = self._model.elements[int(target)]
            except:
                # Case of element set
                target = target.split('.')[0] + '-' + target.split('.')[1]
                target = self._model.groups[target.strip()]

            #print target.get_name()

            JDLTYPE = int(line_split[1].split('U')[1])
            #JDLTYPE = int(line_split[1].strip().upper())
            DLOAD.set_target(target)
            DLOAD.set_JDLTYPE(JDLTYPE)
            try:
                # Case of scalar load magnitude
                ADLMAG  = float(line_split[2]) #np.array(line_split[2:],dtype=np.float64)
                DLOAD.set_ADLMAG(ADLMAG)
                if len(line_split)>3:
                    DLOAD.set_additionalInfos(np.array(line_split[3:],dtype=np.float64))
            except:
                # Case of magnitude defined by a nodal distribution
                DLOAD.set_ADLMAG(0.)
                if self._model.tables[line_split[2].strip()].opts['location'].strip() == "NODE":
                    DLOAD.set_additionalInfos(np.array([self._model.tables[line_split[2].strip()].get_label()],dtype=np.float64))
                    self._model.tables[line_split[2].strip()].opts['dtype'] = 'distributedload'
                else:
                    raise Exception("%s is not a field defined at nodes" % line_split[2].strip())


            self._model.steps[self._step.get_name()].loads.add(DLOAD)
            print('\t\t - Distributed load (%s) successfully added.' % DLOAD.get_name())

        return None


    def _read_cload(self,line):
        print('\t\tReading Cload ...')
        while 1:
            pos  = self._file.tell()
            line = self._file.readline()
            if line.startswith('*'):
                self._file.seek(pos)
                break

            self._nloads += 1
            CLOAD = Load(label=self._nloads, name='no name', type='cload')

            line_split = line.rstrip().split(',')
            target = line_split[0]
            try:
                # Case of node label
                target = self._model.nodes[int(target)]
            except:
                # Case of node set
                target = target.split('.')[0] + '-' + target.split('.')[1]
                target = self._model.groups[target.strip()]

            JDLTYPE = int(line_split[1])
            ADLMAG  = float(line_split[2])
            CLOAD.set_target(target)
            CLOAD.set_ADLMAG(ADLMAG)
            CLOAD.set_JDLTYPE(JDLTYPE)
            self._model.steps[self._step.get_name()].loads.add(CLOAD)
            print('\t\t - Concentrated load (%s) successfully added.' % CLOAD.get_name())

        return None
