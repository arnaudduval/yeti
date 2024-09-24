# Copyright 2016-2020 Thibaut Hirschler
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

from .read   import *
from .write  import *
from .common import Container
import numpy as np

class PART(object):
    def __init__(self, name=None):
        self.name  = name
        self.nodes      = Container()
        self.elements   = Container()
        self.properties = Container()
        self.groups     = Container()

class MechanicalModel(PART):
    def __init__(self, name=None):
        PART.__init__(self, name)
        self.coords     = Container()
        self.tables     = Container()
        self.materials  = Container()
        #self.properties = Container()
        self.steps      = Container()
        self.info       = {}

    def read(self, format, filename, config=default):
        if format == "ABAQUS":
            parser = AbaqusParser(model=self, filename=filename)
        else:
            print('Parser %s not supported, only ABAQUS is available' % format)
            return
        parser.configure(**config)
        parser.parse()

    def write(self, format, filename, config=default):
        if   format == "NASTRAN":
            writer = NastranWriter(model=self, filename=filename)
        elif format == "ABAQUS":
            writer = AbaqusWriter(model=self, filename=filename)
        elif format == "IDEAS":
            writer = UnvWriter(model=self, filename=filename)
        else:
            print('Writer %s not supported, only NASTRAN, ABAQUS and IDEAS are available' % format)
            return

        writer.configure(**config)
        writer.write()

    def get_nodes(self):
        nb_cp = 0
        for node in self.nodes:
            x,y,z = node.get_coordinates()
            if nb_cp==0:
                COORDS = np.array([x,y,z])
            else:
                COORDS = np.vstack((COORDS,np.array([x,y,z])))
            nb_cp += 1
        return COORDS.transpose(),nb_cp

    def get_ien_patch(self,num_patch):
        for grp in self.groups:
            if 'ELTPATCH%i'%num_patch in grp.get_name().upper():
                break

        #print grp.get_elts_labels()
        for num_elem in np.sort(grp.get_elts_labels()):
            eleme = self.elements.get(num_elem)

            try:
                IEN_patch = np.vstack((IEN_patch,np.array(eleme.get_node_labels())))
            except:
                IEN_patch = np.array([eleme.get_node_labels()])
        return IEN_patch

    def get_ien(self):
        IEN = []
        nb_patch = 0
        for grp in self.groups:
            grp_name = grp.get_name()
            if 'ELTPATCH' in grp_name.upper():
                nb_patch += 1
        for num_patch in range(1,nb_patch+1):
            IEN.append(self.get_ien_patch(num_patch))
        return IEN

    def get_bc(self,num_step=2):
        nb_bc = 0
        bc_target = []
        bc_target_nbelem = []
        bc_values = np.vstack([[],[]])
        step  = self.steps.get(num_step)
        for bc in step.bcs:
            thisbc_target = bc.get_target().get_nodes_labels()
            thisbc_values = bc.get_values()
            for ddl in range(0,6):
                if not thisbc_values[ddl]==None:
                    nb_bc += 1
                    bc_target.append(np.array(thisbc_target,dtype=np.intp))
                    bc_target_nbelem.append(len(thisbc_target))
                    disp = thisbc_values[ddl]
                    if nb_bc>1:
                        bc_values = np.insert(bc_values, [nb_bc-1], [[ddl+1],[disp]], axis=1)
                    else:
                        bc_values = np.vstack([[ddl+1],[disp]])
        return bc_target, bc_values, np.array(bc_target_nbelem,dtype=np.intp), nb_bc

    def get_load(self,num_step=2):
        nb_load = 0
        inddload= []
        jdltype = []
        adlmag  = []
        load_target_nbelem = []
        additionalLoadInfos= []
        step  = self.steps.get(num_step)
        for load in step.loads:
            nb_load += 1
            load_grp = load.get_target()
            if load.get_type().lower()=='dload':
                thisload_ind = load_grp.get_elts_labels()
            elif load.get_type().lower()=='cload':
                try:
                    thisload_ind = load_grp.get_nodes_labels()
                except:
                    thisload_ind = [load_grp.get_label()]
            inddload.append(thisload_ind)
            adlmag.append(load.get_values()[0])
            jdltype.append(load.get_values()[1])
            load_target_nbelem.append(len(thisload_ind))
            additionalLoadInfos.append(load.get_additionalLoadInfos())
        return np.array(inddload),np.array(jdltype),np.array(adlmag),np.array(load_target_nbelem),\
            nb_load,additionalLoadInfos


    def get_parameters(self):
        idx      = []
        nnode    = []
        tensor   = []
        mcrd     = []
        elt_type = []
        nbint    = []
        for grp in self.groups:
            grp_name = grp.get_name()
            i = grp_name.upper().find('ELTPATCH')
            if i>=0:
                num_patch = int(grp_name[i+len('ELTPATCH'):])
                idx.append(num_patch)
                elt_type_patch = grp.elements.get_one().get_type()
                for elem in grp.elements:
                    if not elem.get_type()==elt_type_patch:
                        print('error: elements of patch %i have not the same elt_type.' % num_patch)
                        return None
                nnode.append (self.elements[elt_type_patch].get_nnode())
                tensor.append(self.elements[elt_type_patch].get_tensor())
                mcrd.append  (self.elements[elt_type_patch].get_coordinates())
                nbint.append (self.elements[elt_type_patch].get_integration())
                elt_type.append(elt_type_patch)
        mcrd = max(mcrd)
        idx  = np.array(idx).argsort()
        nnode    = np.array(nnode   )[idx]
        tensor   = np.array(tensor  )[idx]
        elt_type = np.array(elt_type)[idx]
        nbint    = np.array(nbint   )[idx]
        return elt_type,nbint,tensor,mcrd,nnode


    def get_material_properties(self):
        """
        Return material properties for each patche read from inp file

        Parameters
        ----------


        Returns
        -------
        material_properties : numpy.ndarray(float)
            Bidimensional array containing material properties for each patch.
            First index is property index (size = size of the material with max value of material properties)
            Second index is patch index

        n_mat_props : numpy.array(int)
            Array containing the number of effective material properties for each patch
        """
        from operator import itemgetter
        MAT_all = []
        for props in self.properties:
            props_name = props.get_name()
            if props_name.upper().startswith('UEL'):
                elset_name   = props.get_elset()
                test_ispatch = elset_name.lower().find('eltpatch')
                if test_ispatch>=0:
                    num_patch = int(elset_name[test_ispatch+len('eltpatch'):])
                    MAT = self.materials[props.get_material()]
                    if   MAT.get_type()=='Isotropic':
                        E  = MAT.get_e()
                        nu = MAT.get_nu()
                        rho= MAT.get_density()
                    elif MAT.get_type()=='Orthotropic':
                        print('warning: material type %s not available yet.' % MAT.get_type())
                        E  = MAT.get_ex()
                        nu = MAT.get_nuxy()
                        rho= MAT.get_density()
                        print(' --> replaced by Isotropic material with E=%1.2f and nu=%1.2f' % (E,nu))
                    elif MAT.get_type()=='HighOrderElastic':
                        lambd = MAT.get_lambda()
                        mu = MAT.get_mu()
                        a = MAT.get_a()
                        b = MAT.get_b()
                        c = MAT.get_c()
                        rho = MAT.get_density()
                    else:
                        print('error: material type %s is unknown' % MAT.get_type())
                        return None
                    if rho==None:
                        rho= 0.
                    if MAT.get_type()=='HighOrderElastic':
                        MAT_all.append([num_patch, lambd, mu,
                                        a[0], a[1], a[2], a[3], a[4],
                                        b[0], b[1], b[2], b[3], b[4],
                                        b[5], b[6], b[7],
                                        c[0], c[1], c[2], rho])
                    else:
                        MAT_all.append([num_patch, E, nu, rho])

        # If there are several materials with different number of properties
        # dimension array as the largest number of props and save number of
        # used properties in array n_mat_props

        # sort by patych number
        MAT_all = sorted(MAT_all, key=itemgetter(0))
        n_mat_props = np.array([len(mat)-1 for mat in MAT_all], dtype=int)
        material_properties = np.zeros((max(n_mat_props), len(MAT_all)), dtype=float)
        for ipatch in range(len(MAT_all)):
            material_properties[:n_mat_props[ipatch],ipatch] = MAT_all[ipatch][1:]


        # OLD VERSION FOR ONLY ELASTIC MATERIALS WITH ALL THE SAME NUMBER OF PROPS
        #MAT_all = np.array(MAT_all)
        #MAT_all = MAT_all[MAT_all[:,0].argsort()] # --sort by patch number
        #material_properties = MAT_all[:,1:].transpose()

        return material_properties, n_mat_props


    def get_properties(self):
        properties  = []
        jproperties = []
        idx = []
        for props in self.properties:
            props_name = props.get_name()
            if props_name.upper().startswith('UEL'):
                elset_name   = props.get_elset()
                test_ispatch = elset_name.lower().find('eltpatch')
                if test_ispatch>=0:
                    num_patch = int(elset_name[test_ispatch+len('eltpatch'):])
                    idx.append(num_patch)
                    jproperties.append(props.get_jprops())
                    properties.append( props.get_props())

        jproperties = np.array(jproperties)
        properties  = np.array(properties, dtype=object)
        idx = np.array(idx).argsort()
        return list(properties[idx]),jproperties[idx]

    def get_tables(self):
        return list(self.tables)

    def _updateshape(self):
        if not hasattr(self,'_design_parameters'):
            return None
        if len(self._design_parameters)>0:
            for key,value in self._design_parameters.items():
                paramname = key.translate({ord('<'): None,ord('>'): None})
                exec("%s=%s"%(paramname,value))
            exec(self._parameters_def)
            for key in self._parameters.keys():
                self._parameters[key] = eval(key)
            for node in self.nodes:
                if hasattr(node,'update_coords'):
                    node.update_coords(self._parameters)

    def _get_shapeparametrization(self):
        if not hasattr(self,'_design_parameters'):
            designvar = {}
        else:
            designvar = self._design_parameters
        coords_str = 'coords = np.array([\n'
        for node in self.nodes:
            if hasattr(node,'get_parametrization'):
                s = '\t[%s,%s,%s]'%tuple(node.get_parametrization())
                coords_str += s.translate({ord('<'): None,ord('>'): None})
            else:
                coords_str += '\t'+str(node.get_coordinates().tolist())
            coords_str += ',\n'
        coords_str+= '\t]).T'
        shapeparam_str = 'from math import *\nimport numpy as np\n'
        if hasattr(self,'_parameters_def'):
            shapeparam_str += self._parameters_def

        return designvar,shapeparam_str+coords_str

    def __str__(self):
        return '\nModel is composed of :\n - %i nodes\n - %i elements\n - %i properties\n' + \
            ' - %i materials\n - %i groups\n - %i tables\n - %i steps\n' % \
            (len(self.nodes), len(self.elements), len(self.properties),
             len(self.materials), len(self.groups), len(self.tables), len(self.steps))



