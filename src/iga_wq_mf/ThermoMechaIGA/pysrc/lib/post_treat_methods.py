"""
.. This module contains functions used to run a simulation and plot the results
.. Joaquin Cornejo
"""

from .__init__ import *

# My libraries
from .create_geomdl import geomdlModel
from .fortran_mf_iga import fortran_mf_iga
from .fortran_mf_wq import fortran_mf_wq

def plot_iterative_solver(filename, inputs:dict, extension='.png', threshold=1.e-16):
    
    # Get new name
    savename = filename.split('.')[0] + extension
    
    # Define inputs
    residue = inputs["Res"]
    error = inputs["Error"]
    method_list = inputs["Methods"]

    # Get new names
    new_method_list = []
    for pcgmethod in method_list:
        if pcgmethod == "WP": new_method_list.append('Without preconditioner')
        elif pcgmethod == "C": new_method_list.append('Fast Diagonalisation (FD)')
        elif pcgmethod == "TDS": new_method_list.append('FD + tensor decomp. + scaling') 
        elif pcgmethod == "JMC": new_method_list.append('FD + jacobien mean')  
        elif pcgmethod == "TDC": new_method_list.append('FD + tensor decomp.') 
        elif pcgmethod == "JMS": new_method_list.append('FD + jacobien mean + scaling')
    
    # Set figure parameters
    fig, [ax1, ax2] = plt.subplots(nrows=1, ncols=2, figsize=(15,6))
    makers = ['o', 'v', 'X', 's', '+', 'p']
    for i, pcgmethod in enumerate(new_method_list):
        error_method = np.asarray(error[i])
        residue_method = np.asarray(residue[i])

        error_method = error_method[residue_method>threshold]
        residue_method = residue_method[residue_method>threshold]
        
        ax1.semilogy(np.arange(len(error_method)), error_method*100, '--',
                    label=pcgmethod, marker=makers[i])
        ax2.semilogy(np.arange(len(residue_method)), residue_method*100, '--',
                    label=pcgmethod, marker=makers[i])

    # Set properties
    for ax in [ax1, ax2]: ax.set_xlabel('Number of iterations')
    ax1.set_ylabel('Relative error (%)')
    ax2.set_ylabel('Relative residue (%)')
    ax2.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    fig.tight_layout()
    fig.savefig(savename)

    return

class ThermalSimulation():

    def __init__(self, inputs:dict, folder=None):

        # Get data to run simulation
        self._degree = inputs.get('degree', 2)
        self._cuts = inputs.get('cuts', 2)
        self._geoCase = inputs.get('case', 'TR')
        self._isIGA = inputs.get('isIGA', False)
        self._funPowDen = inputs.get('funPowDen', None)
        self._funTemp = inputs.get('funTemp', None)   
        self._nbIter = 100
        self._iterMethods = inputs.get('IterMethods', ['WP'])
        self._isOnlyIter = inputs.get('isOnlyIter', True)
        self._thermalModel = None

        # Get filename
        filename = self._get_filename()
        self._filename = folder + filename

        return

    def _get_filename(self):
        " Set filename using input information "
        # Get text file name
        filename = (self._geoCase 
                    + '_p_' + str(self._degree) 
                    + '_nbel_' + str(2**self._cuts)
        )
        if self._isIGA: filename += '_IGAG'
        else: filename += '_IGAWQ'
        filename += '.txt'

        return filename

    def create_geometryModel(self):
        " Create model using YETI refinement "
        # Set properties
        degree = self._degree
        cuts = self._cuts
        geometry = {'degree':degree*np.ones(3, dtype=int)}

        # Create geometry using geomdl
        modelGeo = geomdlModel(self._geoCase, **geometry)

        # Refine mesh using YETI function
        modelIGA = modelGeo.export_IGAparametrization(nb_refinementByDirection=cuts*np.ones(3, dtype=int))

        return modelIGA

    def create_thermalModel(self, material=None, Dirichlet=None):
        " Creates thermal model in IGA-G or IGA-WQ approach "

        # Creates geometry
        modelIGA = self.create_geometryModel()

        # Create thermal model
        if self._isIGA: thermalModel = fortran_mf_iga(modelIGA, material=material, Dirichlet=Dirichlet)
        else: thermalModel = fortran_mf_wq(modelIGA, material=material, Dirichlet=Dirichlet)

        return thermalModel

    def compute_Fn_SM(self, model:fortran_mf_wq, funPowDen=None, funTemp=None):
        """ Compute vector b = Fn - Knd.ud where F is source vector and K is conductivity matrix. 
            This equation is used in substitution method (SM)
        """
        
        # Initialize
        dof = model._thermal_dof
        dod = model._thermal_dod

        if funTemp is not None:  
            ud = model.interpolate_ControlPoints(funfield=funTemp)[dod]
            u = np.zeros(model._nb_ctrlpts_total); u[dod] = ud
            Fn = model.eval_source_vector(funPowDen, indi=dof) 
            Knd_ud = model.eval_Ku(u)
            Fn -= Knd_ud
        else:
            ud = None
            Fn = model.eval_source_vector(funPowDen, indi=dof)
            
        return Fn, ud
    
    def run_iterative_solver(self, model:fortran_mf_wq, Fn=None, nbIterPCG=100, 
                            threshold=1e-15, method='WP', solDir=None):
        " Solve steady heat problems using iterative solver "
        
        if solDir is None: 
            solDir = np.ones(len(Fn))
            print("Direct solution unknown. Default: ones chosen. Be aware of error results")

        # Run solver
        start = time.process_time()
        if model._thermalDirichlet is None: raise Warning('Ill conditionned. It needs Dirichlet conditions')
        un, residue, error = model.MFsteadyHeat(Fn, nbIterPCG, threshold, method, solDir)
        stop = time.process_time()
        time_itersolver = stop - start 

        return un, residue, error, time_itersolver

    def run_direct_solver(self, model:fortran_mf_wq, Fn=None):
        " Solve steady heat problems using direct solver "

        # Assemble conductivity matrix 
        start = time.process_time()
        if model._thermalDirichlet is None: raise Warning('Ill conditionned. It needs Dirichlet conditions')
        Knn = model.eval_conductivity_matrix(indi=model._thermal_dof, indj=model._thermal_dof)
        stop = time.process_time()
        time_assembly = stop - start 

        # Solve system
        start = time.process_time()
        un = sp.linalg.spsolve(Knn, Fn)
        stop = time.process_time()
        time_dirsolver = stop - start 

        return un, time_assembly, time_dirsolver

    def run_simulation(self, material=None, Dirichlet=None):
        " Runs simulation using given input information "

        # Initialize
        doDirect = True
        if self._isOnlyIter: doDirect = False
        time_assembly, time_dirsolver = 0.0, 0.0
        time_noiter, time_itersolver = 0.0, 0.0
        residue, error = np.zeros(self._nbIter+1), np.zeros(self._nbIter+1)

        # Create thermal model
        self._thermalModel = self.create_thermalModel(material=material, Dirichlet=Dirichlet)

        # Create force vector F = Fn - Knd ud 
        Fn = self.compute_Fn_SM(self._thermalModel, self._funPowDen, self._funTemp)[0]

        # Run direct solver 
        if doDirect: solDir, time_assembly, time_dirsolver = self.run_direct_solver(self._thermalModel, Fn)
        else: solDir = None

        # Compute time to prepare method before iterations
        time_noiter = []
        for itmethod in self._iterMethods:
            time_temp = self.run_iterative_solver(self._thermalModel, Fn=Fn, nbIterPCG=0, method=itmethod, solDir=solDir)[3]
            time_noiter.append(time_temp)

        # Run iterative methods
        time_itersolver, residue, error = [], [], []
        for itmethod in self._iterMethods:
            _, residue_t, error_t, time_temp = self.run_iterative_solver(self._thermalModel, Fn=Fn, 
                                                                        method=itmethod, solDir=solDir)
            time_itersolver.append(time_temp)
            residue.append(residue_t); error.append(error_t)
                
        # Write file
        output = {"Methods": self._iterMethods, "TimeAssembly": time_assembly, "TimeDirect": time_dirsolver, 
                "TimeNoIter":time_noiter, "TimeIter": time_itersolver, "Res": residue, "Error": error}
        self.write_text_file(output)

        return 

    def write_text_file(self, inputs:dict): 
        " Writes and exports simulation data in txt file "
        
        # Define inputs
        iterMethods = inputs['Methods']
        time_assembly = inputs['TimeAssembly']
        time_dirsolver = inputs['TimeDirect']
        time_itersolver = inputs['TimeIter']
        time_noiter = inputs['TimeNoIter']
        residue = inputs['Res']
        error = inputs['Error']
        
        # Write file
        with open(self._filename, 'w') as f:
            f.write('** RESULTS **\n')
            f.write('** Direct solver\n')
            f.write('*Time assembly\n')
            f.write('{:E}\n'.format(time_assembly))
            f.write('*Time direct\n')
            f.write('{:E}\n'.format(time_dirsolver))
            f.write('** Iterative solver ' + ','.join([item for item in iterMethods]) + '\n')
            f.write('** Number of iterations ' + '{:d}\n'.format(self._nbIter))

            for i, method in enumerate(iterMethods):
                f.write('**' + method + '\n')
                f.write('*Time prepare ' + method +'\n')
                f.write('{:E}\n'.format(time_noiter[i]))
                f.write('*Time iter ' + method +'\n')
                f.write('{:E}\n'.format(time_itersolver[i]))
                f.write('*Residue ' + method + '\n')
                f.writelines(['{:E}'.format(res) + ',' + '{:E}'.format(err) + '\n'
                                for res, err in zip(residue[i], error[i])]) 
        return

class SimulationData(): 

    def __init__(self, filename=None): 

        # Set filename
        self._filename = filename
        
        # Get important data from title
        full_path = os.path.realpath(self._filename)
        name = os.path.basename(full_path).split('.')[0]
        self._interpret_title(name)

        # Get data from file
        self._dataSimulation = self.get_infos_simulation()
    
        return

    # Read data

    def _interpret_title(self, name):
        " Gets important information from title "

        # Split string
        ls = name.split('_')

        # Get data
        self._geoCase = ls[0]
        self._degree = int(ls[2])
        self._cuts = int(np.log2(int(ls[4])))

        if ls[5] == 'IGAG': self._isIGA = True
        else: self._isIGA = False

        return

    def get_infos_simulation(self):
        " Gets the information of the simulation "

        lines = self._get_cleanLines()
        iterMethods, nbIter = self._read_iter_methods(lines)

        time_noiter, time_itersolver, residue, error = [], [], [], []
        for itmethod in iterMethods:
            tp = self._read_time_preparation(lines, itmethod)
            ti = self._read_time_iterations(lines, itmethod)
            res, err = self._read_residue_error(lines, itmethod, nbIter)

            time_noiter.append(tp); time_itersolver.append(ti)
            residue.append(res); error.append(err)

        output = {"Methods": iterMethods, "TimeNoIter":time_noiter, "TimeIter": time_itersolver, 
                    "Res": residue, "Error": error}
        return output

    def _get_cleanLines(self):
        " Gets the lines from the file that could have important information "

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
    
    def _read_residue_error(self, lines, method, length=100):
        i = self._get_num_line(lines,'*Residue ' + method) + 1
        lines_data = lines[i:i+length+1]

        residue, error = [], []
        for line in lines_data:
            ls = line.split(',')
            residue.append(float(ls[0]))
            error.append(float(ls[1]))

        return residue, error
    