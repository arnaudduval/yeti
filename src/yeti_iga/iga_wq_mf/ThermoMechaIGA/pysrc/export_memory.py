"""
.. Test memory
.. We compare memory usage of diferent methods using mprofiler library
.. Joaquin Cornejo 
"""

# Python libraries
import pandas as pd 
import matplotlib.pyplot as plt
import os 

# Chosee folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/'

# Select files
fileList = ['mprofile_I.dat']
# labels = ['Direct Method', 'Iterative method']

plt.figure(1)
ax = plt.gca()
for i, file in enumerate(fileList): 
    color = next(ax._get_lines.prop_cycler)['color']

    # Import data
    data = pd.read_table(folder + file, sep=' ', names=['MEM', 'memory', 'time'])
    memory =  data.memory
    time = data.time - data.time[0]

    if i == 1:
        time_zoomed = time[time<80]
        memory_zoommed = memory[time<80]

    plt.plot(time, memory, color=color)
    # plt.plot(time, memory, label=labels[i], color=color)

# Set properties
plt.grid()
plt.xlabel("Time (s)", fontsize= 16)
plt.ylabel("Memory used (in MiB)", fontsize= 16)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.legend(loc='best')
plt.tight_layout()

# # Plot the zoomed portion
# # Location for the zoomed portion 
# sub_axes = plt.axes([.38, .3, .25, .25]) 
# sub_axes.plot(time_zoomed, memory_zoommed, color=color) 
# plt.xticks(fontsize=16)
# plt.yticks(fontsize=16)
# plt.grid()

# Maybe we can save images in svg format (vectorized)
plt.savefig(folder + 'memoryusage' + '.png')
