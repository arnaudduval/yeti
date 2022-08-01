# Python libraries
import numpy as np
from matplotlib import pyplot as plt

# My libraries
from lib.D3viscoplasticity import stdcst, stkronst, one_kron_one

# Define constants
d = 3
ddl = int(d*(d+1)/2)

A = np.array([1, 2, 3, 4, 5, 6])
B = np.array([1, 2, 3, 4, 5, 6])

onekronone = one_kron_one(d)
print(onekronone)