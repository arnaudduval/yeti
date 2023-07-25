"""
Plot error as a function of refinement for different degrees
for the clamped square case
"""


import numpy as np
import matplotlib.pyplot as plt

error = np.loadtxt('convergence.txt')

refinements = 2**(np.arange(error.shape[1]-1)+1)
print(refinements)


for degree in range(error.shape[0]):
    plt.plot(refinements, error[degree, 1:], 'o-',
             label="degree "+str(int(error[degree, 0])))

plt.legend()
plt.xscale("log")
plt.yscale("log")
plt.show()
