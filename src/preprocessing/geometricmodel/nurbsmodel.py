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

import numpy as np
import os

class NBfile(object):
    
    def __init__(self,filename):
        self._filename = filename
        
    def getGeoInfos(self):
        lines = self._get_cleanLines()
        nb_patch = self._read_nb_patch(lines)
        dim      = self._read_dimension(lines,nb_patch)
        nnode    = self._read_NNODE(lines,nb_patch)
        nb_elem  = self._read_nb_elem(lines)
        ukv,nkv  = self._read_knots(lines,dim,nb_patch)
        jpqr     = self._read_jpqr(lines,dim,nb_patch)
        nijk     = self._read_nijk(lines,nb_elem)
        nb_elem_patch = self._read_nbelembypatch(lines)
        weight   = self._read_weight(lines,nb_patch,nb_elem_patch,nnode)
        geoSet   = [dim,nkv,ukv,jpqr,nijk,weight,nb_elem_patch,nb_patch,nnode,nb_elem]
        return geoSet        
    
    def _get_cleanLines(self,filename=None):
        if not filename:
            inputFile  = open(self._filename, 'r')
        else:
            inputFile  = open(filename, 'r')
        lines      = inputFile.readlines()
        linesClean = []
        theLine    = ''
        # Joining lines ending with ',' in a new list of lines
        # This allow to merge definitions written on more than a single line
        for i in range(0, len(lines)):
            words = lines[i].rstrip().split(',')
            # Removing trailing spaces
            lastWord = words[-1].split()
            if(len(lastWord) == 0):
                theLine = theLine + lines[i]
                # Removing '\n' character
                theLine = theLine.rstrip()
            else:
                theLine = theLine + lines[i]
                linesClean.append(theLine.rstrip())
                theLine = ''
        inputFile.close()
        return linesClean

    def _get_num_line(self,lines,str2find):
        i=0
        for line in lines:
            if line.lower().startswith(str2find):
                break
            i += 1
        if i==len(lines):
            print("Error: keyword " + str2find + " is missing in NB file.")
        return i
        
    def _read_dimension(self,lines,nb_patch):
        i = self._get_num_line(lines,'*dimension')
        dimension = np.fromstring(lines[i+1],dtype=np.intp,sep=',')
        if dimension.size==1 and nb_patch>1:
            dimension = dimension*np.ones(nb_patch, dtype=np.intp)
        return dimension
    
    def _read_nb_patch(self, lines):
        i = self._get_num_line(lines,'*number of patch')
        nb_patch = np.intp(lines[i+1])
        return nb_patch
        
    def _read_nb_elem(self, lines):
        i = self._get_num_line(lines,'*total number of element')
        nb_elem = np.intp(lines[i+1])
        return nb_elem
        
    def _read_nbelembypatch(self, lines):
        i = self._get_num_line(lines,'*number of element by patch')
        nb_elem_patch = np.fromstring(lines[i+1],dtype=np.intp,sep=',')
        return nb_elem_patch

    def _read_NNODE(self,lines,nb_patch):
        i = self._get_num_line(lines,'*number of cp by element')
        NNODE = np.fromstring(lines[i+1],dtype=np.intp,sep=',')
        if NNODE.size==1 and nb_patch>1:
            NNODE = NNODE*np.ones(nb_patch, dtype=np.intp)
        return NNODE
    
    def _read_knots(self,lines,dimension,nb_patch):
        Ukv = []
        Nkv = np.zeros((3,nb_patch), dtype=np.int64)
        for num_patch in np.arange(0,nb_patch):
            ukv_patch,nkv_patch = self._read_knots_patch(lines,dimension,num_patch+1)
            Ukv.append(ukv_patch)
            Nkv[:,num_patch] = nkv_patch[:]
            
        #Ukv = np.array(Ukv)
        return Ukv,Nkv
        
    def _read_knots_patch(self,lines,dimension,num_patch):
        ukv_patch = []
        nkv_patch = np.zeros(3, dtype=np.intp)
        
        i = self._get_num_line(lines,'*patch(%i)'%num_patch)
        i += 1
        for n in np.arange(0,dimension[num_patch-1]):
            nkv_patch[n] = np.float64(lines[i].rstrip())
            ukv_patch.append(np.fromstring(lines[i+1],dtype=np.float64,sep=','))
            i += 2
        return ukv_patch,nkv_patch
    
    def _read_jpqr(self,lines,dimension,nb_patch):
        jpqr = np.zeros((3,nb_patch), dtype=np.intp)
        i = self._get_num_line(lines,'*jpqr')
        jpqr[:dimension[0],0] = np.fromstring(lines[i+1],dtype=np.intp,sep=',')
        if nb_patch>1:
            i += 1
            try:
                lines[i+1].startswith('*')
            except IndexError:
                jpqr[:,1:] = np.repeat(np.vstack(jpqr[:,0]), nb_patch-1, axis=1)
            else:
                if lines[i+1].startswith('*'):
                    jpqr[:,1:] = np.repeat(np.vstack(jpqr[:,0]), nb_patch-1, axis=1)
                else:
                    for j in np.arange(1,nb_patch):
                        jpqr[:dimension[j],j] = np.fromstring(lines[i+1],dtype=np.intp,sep=',')
                        i += 1
        return jpqr
    
    def _read_nijk(self,lines,nb_elem):
        i = self._get_num_line(lines,'*nijk')
        nijk = np.zeros((3,nb_elem), dtype=np.intp)
        for num_elem in np.arange(0,nb_elem):
            nijk_elem = np.fromstring(lines[i+1],dtype=np.intp,sep=',')
            nijk[:nijk_elem.size-1,num_elem] = nijk_elem[1:]
            i += 1
        return nijk
    
    def _read_weight(self,lines,nb_patch,nb_elem_patch,NNODE):
        i = self._get_num_line(lines,'*weight')
        weight = []
        for num_patch in np.arange(0,nb_patch):
            #weight_patch = np.zeros((NNODE[num_patch],nb_elem_patch[num_patch]), dtype=np.float64)
            for num_elem in np.arange(0,nb_elem_patch[num_patch]):
                weight_elem = np.fromstring(lines[i+1],dtype=np.float64,sep=',')
                #weight_patch[:,num_elem] = weight_elem[1:]
                i += 1
                weight.append(weight_elem[1:])
            #weight.append(weight_patch)
        #weight = np.array(weight)
        return weight
