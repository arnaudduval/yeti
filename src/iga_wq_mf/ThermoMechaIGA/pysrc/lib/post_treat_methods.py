"""
.. module :: Post - treatment
    :synopsis: Provides functions used to plot our results
.. author :: Joaquin Cornejo
"""
# Python libraries
import numpy as np
import matplotlib.pyplot as plt
import tracemalloc
import scipy, numpy as np
import os, time

# My libraries
from . import enablePrint, blockPrint
from .create_geomdl import create_geometry
from .fortran_mf_iga import fortran_mf_iga
from .fortran_mf_wq import fortran_mf_wq

def plot_iterative_solver(filename, inputs, extension ='.png'):
    
    # Get new name
    savename = filename.split('.')[0] + extension
    
    # Define color
    CB_color_cycle = ['#377eb8', '#ff7f00', '#4daf4a',
                    '#f781bf', '#a65628', '#984ea3',
                    '#999999', '#e41a1c', '#dede00']

    # Define inputs
    residue = inputs["Res"]
    error = inputs["Error"]
    method_list = inputs["Methods"]

    # Get new names
    new_method_list = []
    for i, pcgmethod in enumerate(method_list):
        if pcgmethod == "WP": new_method_list.append('Without preconditioner')
        elif pcgmethod == "C": new_method_list.append('Fast Diagonalisation (FD)')
        elif pcgmethod == "TDS": new_method_list.append('FD + tensor decomp. + scaling') 
        elif pcgmethod == "JM": new_method_list.append('FD + jacobien mean')  
        elif pcgmethod == "TD": new_method_list.append('FD + tensor decomp.') 
        elif pcgmethod == "JMS": new_method_list.append('FD + jacobien mean + scaling')

    # Select important values
    tol = 1.e-16
    
    # Set figure parameters
    fig, [ax1, ax2] = plt.subplots(nrows=1, ncols=2, figsize=(12,6))

    for i, pcgmethod in enumerate(new_method_list):
        error_method = np.asarray(error[i])
        residue_method = np.asarray(residue[i])

        error_method = error_method[residue_method>tol]
        residue_method = residue_method[residue_method>tol]
        
        ax1.semilogy(np.arange(len(error_method)), error_method*100, label=pcgmethod, color=CB_color_cycle[i])
        ax2.semilogy(np.arange(len(residue_method)), residue_method*100, label=pcgmethod, color=CB_color_cycle[i])

    # Set properties
    ax1.set_xlabel('Number of iterations', fontsize=16)
    ax1.set_ylabel('Relative error (%)', fontsize=16)
    ax1.tick_params(axis='x', labelsize=16)
    ax1.tick_params(axis='y', labelsize=16)
    ax1.legend(loc='best', fontsize=12)
    ax1.grid()

    ax2.set_xlabel('Number of iterations', fontsize=16)
    ax2.set_ylabel('Relative residue (%)', fontsize=16)
    ax2.tick_params(axis='x', labelsize=16)
    ax2.tick_params(axis='y', labelsize=16)
    ax2.legend(loc='best', fontsize=12)
    ax2.grid()

    # Save figure
    plt.tight_layout()
    plt.savefig(savename)
    plt.close(fig)

    return

class ThermalSimulation():

    def __init__(self, inputs: dict, folder):

        # Get essential data to run simulation
        self._degree = inputs.get('degree', 2)
        self._cuts = inputs.get('cuts', 2)
        self._geoCase = inputs.get('case', 'TR')
        self._isIGA = inputs.get('isIGA', False)
        self._isCG = inputs.get('isCG', True)
        self._funPowDen = inputs.get('funPowDen', None)
        self._funTemp = inputs.get('funTemp', None)   

        # Get extra data to run simulation  
        self._iterMethods = inputs.get('IterMethods', ['WP'])
        self._isOnlyDirect = inputs.get('isOnlyDirect', False)
        self._isOnlyIter = inputs.get('isOnlyIter', True)

        # Get filename
        filename = self._get_filename()
        self._filename = folder + filename + '.txt'
        
        # Create thermal mode
        self._thermalModel = None

        # Create source vector 
        self._Fn = None

        # Set number of iterations
        self._nbIter = 100
        
        return

    def _get_filename(self):

        # Get text file name
        filename = self._geoCase + '_p_' + str(self._degree) + '_nbel_' + str(2**self._cuts)
        if self._isIGA: filename += '_IGAG'
        else: filename += '_IGAWQ'
        if self._isCG: filename += '_CG'
        else: filename += '_BiCG'

        return filename

    def write_text_file(self, inputs: dict): 

        # Define inputs
        time_assembly = inputs['TimeAssembly']
        time_direct = inputs['TimeDirect']
        memory_direct = inputs['MemDirect']
        
        iterMethods = inputs['Methods']
        time_noiter = inputs['TimeNoIter']
        memory_noiter = inputs['MemNoIter']
        time_iter = inputs['TimeIter']
        memory_iter = inputs['MemIter']
        residue = inputs['Res']
        error = inputs['Error']
        
        # Write file
        with open(self._filename, 'w') as f:
            f.write('** RESULTS **\n')
            f.write('** Direct solver\n')
            f.write('*Time assembly\n')
            f.write('{:E}\n'.format(time_assembly))
            f.write('*Time direct\n')
            f.write('{:E}\n'.format(time_direct))
            f.write('*Memory direct\n')
            f.write('{:E}\n'.format(memory_direct))
            f.write('** Iterative solver ' + ','.join([item for item in iterMethods]) + '\n')
            f.write('** Number of iterations ' + '{:d}\n'.format(self._nbIter))

            for i, method in enumerate(iterMethods):
                f.write('**' + method + '\n')
                f.write('*Time prepare ' + method +'\n')
                f.write('{:E}\n'.format(time_noiter[i]))
                f.write('*Memory prepare ' + method +'\n')
                f.write('{:E}\n'.format(memory_noiter[i]))
                f.write('*Time iter ' + method +'\n')
                f.write('{:E}\n'.format(time_iter[i]))
                f.write('*Memory iter ' + method +'\n')
                f.write('{:E}\n'.format(memory_iter[i]))
                f.write('*Residue ' + method + '\n')
                f.writelines(['{:E}'.format(res) + ',' + '{:E}'.format(err) + '\n'
                                for res, err in zip(residue[i], error[i])]) 
        return

    # Run simulations

    def create_geometryModel(self):
        geoModel = create_geometry(self._degree, self._cuts, self._geoCase)
        return geoModel

    def create_thermalModel(self, **properties):
        geoModel = self.create_geometryModel()
        if self._isIGA: thermalModel = fortran_mf_iga(geoModel, **properties)
        else: thermalModel = fortran_mf_wq(geoModel, **properties)
        del geoModel
        return thermalModel

    def compute_source(self, model, funPowDen, funTemp=None):
        # Get data
        dof = model._thermal_dof
        dod = model._thermal_dod

        # Assemble source vector F
        if funTemp is not None:  
            Td = model.interpolate_ControlPoints(funTemp)[1]
            Fn = model.eval_source_vector(funPowDen, indi=dof, indj=dod, Td=Td)
        else:
            Td = None
            Fn = model.eval_source_vector(funPowDen, indi=dof)
            
        return Fn, Td
    
    def run_iterative_solver(self, model, Fn, nbIter=100, eps=1e-15, method='WP', solDir=None):

        if solDir == None: 
            solDir = np.ones(len(Fn))
            print("Direct solution unknown. Default: ones chosen. Be aware of residue results")

        # Run matrix free solver
        start = time.process_time()
        solIter, residue, error = model.MFsolver(Fn, nbIter, eps, method, solDir, self._isCG)
        stop = time.process_time()
        time_MF = stop - start 

        return solIter, residue, error, time_MF

    def run_direct_solver(self, model, Fn):
       
        # Assemble conductivity matrix K
        start = time.process_time()
        Knn = model.eval_conductivity_matrix(indi=model._thermal_dof, indj=model._thermal_dof)
        stop = time.process_time()
        time_assembly = stop - start 

        # Solve system
        start = time.process_time()
        solDir = scipy.sparse.linalg.spsolve(Knn, Fn)
        stop = time.process_time()
        time_direct = stop - start 

        return solDir, time_assembly, time_direct

    def run_simulation(self, **properties):

        # Initialize
        time_assembly, time_direct, memory_direct = 0.0, 0.0, 0.0
        time_noiter, memory_noiter = 0.0, 0.0
        time_iter, memory_iter = 0.0, 0.0
        residue, error = np.zeros(self._nbIter+1), np.zeros(self._nbIter+1)

        # Create thermal model
        solDir = None
        self._thermalModel = self.create_thermalModel(**properties)

        # Create source vector 
        self._Fn, Td = self.compute_source(self._thermalModel, self._funPowDen, self._funTemp)

        # Define actions 
        doDirect, doIterative = True, True
        if self._isOnlyDirect: doIterative = False
        if self._isOnlyIter: doDirect = False

        if doDirect :
            solDir, time_assembly, time_direct = self.run_direct_solver(self._thermalModel, self._Fn)

        if doIterative:
            # Only compute time to prepare method before iterations
            time_noiter, memory_noiter = [], []
            for method in self._iterMethods:
                time_noiter_t = self.run_iterative_solver(self._thermalModel, self._Fn, nbIter=0, method=method, solDir=solDir)[3]
                memory_noiter_t = 0.0
                
                # Save data
                time_noiter.append(time_noiter_t)
                memory_noiter.append(memory_noiter_t)

            # With and without preconditioner
            time_iter, memory_iter, residue, error = [], [], [], []
            for method in self._iterMethods:
                Tn_t, residue_t, error_t, time_iter_t = self.run_iterative_solver(self._thermalModel, self._Fn, method=method, solDir=solDir)
                memory_iter_t = 0.0

                # Save data
                time_iter.append(time_iter_t)
                residue.append(residue_t)
                error.append(error_t)
                memory_iter.append(memory_iter_t)

        # Export results
        solution = np.zeros(self._thermalModel._nb_ctrlpts_total)
        dof = self._thermalModel._thermal_dof

        if Td == None:
            solution[dof] = Tn_t
        else: 
            dod = self._thermalModel._thermal_dod
            solution[dod] = Td
            solution[dof] = Tn_t
        
        self._thermalModel.export_results(u_ctrlpts=solution)
        
        # Write file
        output = {"Methods": self._iterMethods, "TimeAssembly": time_assembly, 
                "TimeDirect": time_direct, "MemDirect": memory_direct, 
                "TimeNoIter":time_noiter, "TimeIter": time_iter, "Res": residue, 
                "Error": error, "MemNoIter": memory_noiter, "MemIter": memory_iter}

        self.write_text_file(output)

        return 

class SimulationData(): 

    def __init__(self, filename): 

        # Set filename
        self._filename = filename
        
        # Get important data from title
        full_path = os.path.realpath(filename)
        name = os.path.basename(full_path).split('.')[0]
        self._interpret_title(name)

        # Get data from file
        self._dataDir = self.getInfosDirect()
        self._dataIter = self.getInfosIter()
    
        return

    # Read data

    def _interpret_title(self, filename):

        # Split string
        ls = filename.split('_')

        # Get data
        self._geoCase = ls[0]
        self._degree = int(ls[2])
        self._cuts = int(np.log2(int(ls[4])))

        if ls[5] == 'IGAG': self._isIGA = True
        else: self._isIGA = False

        if ls[6] == 'CG': self._isCG = True
        else: self._isCG = False

        return

    def getInfosDirect(self):
        lines = self._get_cleanLines()
        time_assembly = self._read_time_assembly(lines)
        time_direct = self._read_time_solve(lines)
        memory_direct = self._read_memory_direct(lines)

        output = {"TimeAssembly": time_assembly, "TimeDirect": time_direct, "MemDirect": memory_direct}
        return output

    def getInfosIter(self):
        lines = self._get_cleanLines()
        iterMethods, nbIter = self._read_iter_methods(lines)

        time_noiter, memory_noiter = [], []
        time_iter, memory_iter, residue, error = [], [], [], []
        for method in iterMethods:
            tni = self._read_time_preparation(lines, method)
            ti = self._read_time_iterations(lines, method)
            mni = self._read_memory_preparation(lines, method)
            mi = self._read_memory_iterations(lines, method)
            res, err = self._read_residue_error(lines, method, nbIter)

            time_noiter.append(tni)
            memory_noiter.append(mni)
            time_iter.append(ti)
            memory_iter.append(mi)
            residue.append(res)
            error.append(err)

        output = {"Methods": iterMethods, "TimeNoIter":time_noiter, "TimeIter": time_iter, "Res": residue, 
                "Error": error, "MemNoIter": memory_noiter, "MemIter": memory_iter}
        return output

    def _get_cleanLines(self):
        inputFile  = open(self._filename, 'r')
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

    def _get_num_line(self, lines, str2find):
        i=0
        for line in lines:
            if line.startswith(str2find):
                break
            i += 1
        if i==len(lines):
            print("Error: keyword " + str2find + " is missing in data file.")
        return i

    def _read_time_assembly(self, lines):
        i = self._get_num_line(lines,'*Time assembly')
        time = float(lines[i+1])
        return time

    def _read_time_solve(self, lines): 
        i = self._get_num_line(lines,'*Time direct')
        time = float(lines[i+1])
        return time

    def _read_memory_direct(self, lines): 
        i = self._get_num_line(lines,'*Memory direct')
        memory = float(lines[i+1])
        return memory

    def _read_iter_methods(self, lines):
        i = self._get_num_line(lines,'** Iterative solver')
        ls = lines[i].split(' ')[-1]
        methods = ls.split(',')
        ls = lines[i+1].split(' ')[-1]
        nbIter = int(ls)
        return methods, nbIter

    def _read_time_preparation(self, lines, method):
        i = self._get_num_line(lines,'*Time prepare ' + method)
        time = float(lines[i+1])
        return time

    def _read_time_iterations(self, lines, method):
        i = self._get_num_line(lines,'*Time iter ' + method)
        time = float(lines[i+1])
        return time

    def _read_memory_preparation(self, lines, method):
        i = self._get_num_line(lines,'*Memory prepare ' + method)
        memory = float(lines[i+1])
        return memory

    def _read_memory_iterations(self, lines, method):
        i = self._get_num_line(lines,'*Memory iter ' + method)
        memory = float(lines[i+1])
        return memory
    
    def _read_residue_error(self, lines, method, length=100):
        i = self._get_num_line(lines,'*Residue ' + method)
        i += 1
        lines_data = lines[i:i+length+1]

        residue, error = [], []
        for line in lines_data:
            ls = line.split(',')
            residue.append(float(ls[0]))
            error.append(float(ls[1]))

        return residue, error
    