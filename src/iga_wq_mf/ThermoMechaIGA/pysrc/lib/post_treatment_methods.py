# Python libraries
import numpy as np
import matplotlib.pyplot as plt


def write_text_file(filename, method_list, inputs): 

    # Define inputs
    time_assembly = inputs["TimeAssembly"]
    time_direct = inputs["TimeDirect"]
    memory_direct = inputs["MemDirect"]
    time_iter = inputs["TimeIter"]
    time_noiter = inputs["TimeNoIter"]
    residue = inputs["Res"]
    error = inputs["Error"]
    memory_iter = inputs["MemIter"]
    memory_noiter = inputs["MemNoIter"]

    with open(filename, 'w') as outputfile:
        outputfile.write('** RESULTS **\n')
        outputfile.write('* DIRECT SOLVER\n')
        outputfile.write("{:E}".format(time_assembly) + "\t"
                        +"{:E}".format(time_direct) + "\t" 
                        +"{:E}".format(memory_direct) +"\n"
        )

        for i, pcgmethod in enumerate(method_list):
            outputfile.write('* ITERATIVE SOLVER : ' + pcgmethod + '\n')
            outputfile.write("{:E}".format(time_noiter[i]) + '\t' 
                            +"{:E}".format(memory_noiter[i])+ '\t' 
                            +"{:E}".format(time_iter[i]) + '\t' 
                            +"{:E}".format(memory_iter[i]) +'\n'
            )
            outputfile.writelines(["{:E}".format(res) + "\t"+ "{:E}".format(err)+"\n"
                            for res, err in zip(residue[i], error[i])]) 
    return

def read_text_file(filename):
    # Read file 
    residue_list = []
    error_list = []
    time_iter = []
    time_noiter = []
    memory_iter = []
    memory_noiter = []

    with open(filename, 'r') as inputfile:
        lines = inputfile.readlines()

        text_position = []
        for _ in range(len(lines)):
            if lines[_][0] == '*': 
                text_position.append(_)
        text_position.append(len(lines)+1)

        # We know that our first numerical line is direct method: 
        lines_data_split = lines[2].split('\t')
        time_assembly = float(lines_data_split[0])
        time_direct = float(lines_data_split[1])
        memory_direct = float(lines_data_split[2])
        del text_position[0:2]

        # For the different type of solvers
        for _ in range(len(text_position)-1):
            ind = [text_position[_]+1, text_position[_+1]]
            lines_data = lines[ind[0]:ind[1]]
            lines_data_split = lines_data[0].split('\t')
            time_noiter.append(float(lines_data_split[0]))
            memory_noiter.append(float(lines_data_split[1]))
            time_iter.append(float(lines_data_split[2]))
            memory_iter.append(float(lines_data_split[3]))

            residue_tmp = []
            error_tmp = []
            for i  in range(1, len(lines_data)):
                lines_data_split = lines_data[i].split('\t')
                residue_tmp.append(float(lines_data_split[0]))
                error_tmp.append(float(lines_data_split[1]))
                
            residue_list.append(residue_tmp)
            error_list.append(error_tmp)

        output = {"TimeAssembly": time_assembly, "TimeDirect": time_direct, "MemDirect": memory_direct, 
                "TimeNoIter":time_noiter, "TimeIter": time_iter, "Res": residue_list, 
                "Error": error_list, "MemNoIter": memory_noiter, "MemIter": memory_iter}

    return output

def plot_iterative_solver(filename, inputs, method_list, extension ='.png'):
    # Define color
    CB_color_cycle = ['#377eb8', '#ff7f00', '#4daf4a',
                    '#f781bf', '#a65628', '#984ea3',
                    '#999999', '#e41a1c', '#dede00']

    # Define inputs
    residue = inputs["Res"]
    error = inputs["Error"]

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
    tol = 1.e-17
    
    # Set figure parameters
    fig, [ax1, ax2] = plt.subplots(nrows=1, ncols=2, figsize=(12,6))

    for i, pcgmethod in enumerate(new_method_list):
        error_method = np.asarray(error[i])
        residue_method = np.asarray(residue[i])

        error_method = error_method[residue_method>tol]
        residue_method = residue_method[residue_method>tol]
        
        ax1.plot(np.arange(len(error_method)), error_method*100, label=pcgmethod, color=CB_color_cycle[i])
        ax2.plot(np.arange(len(residue_method)), residue_method*100, label=pcgmethod, color=CB_color_cycle[i])

    # Set properties
    ax1.set_yscale("log")
    ax1.set_xlabel('Number of iterations', fontsize=16)
    ax1.set_ylabel('Relative error (%)', fontsize=16)
    ax1.tick_params(axis='x', labelsize=16)
    ax1.tick_params(axis='y', labelsize=16)
    ax1.legend(loc='best',fontsize=12)
    ax1.grid()

    ax2.set_yscale("log")
    ax2.set_xlabel('Number of iterations', fontsize=16)
    ax2.set_ylabel('Relative residue (%)', fontsize=16)
    ax2.tick_params(axis='x', labelsize=16)
    ax2.tick_params(axis='y', labelsize=16)
    ax2.legend(loc='best',fontsize=12)
    ax2.grid()

    # Save figure
    plt.tight_layout()
    plt.savefig(filename + extension)
    plt.close(fig)

    return
