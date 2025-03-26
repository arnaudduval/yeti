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

default = {"nodes"              : True,  \
           "elements"           : True,  \
           "properties"         : True,  \
           "materials"          : True,  \
           "groups"             : True,  \
           "boundary_conditions": False, \
           "loads"              : False, \
           "steps"              : False, \
           "outputs"            : False, \
          }
          
          
class Parser(object):
    
    def __init__(self, model, filename):
        self._model = model
        self._filename = filename
        #self._file = open(self._filename, 'r', encoding='utf-8')
        self._file = open(self._filename, 'r')
        self._config = default
    
    def configure(self, **kwargs):
        for (key, value) in kwargs.items():
            self._config[key] = value

    def parse(self):
        pass
        
