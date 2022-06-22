"""
.. Test of memory
.. We test which library is better to get memory used
.. Joaquin Cornejo 
"""

# Python libraries
import memory_profiler as mp
import tracemalloc
import os, psutil
import numpy as np
from resource import getrusage, RUSAGE_SELF

# My libraries
from iga_wq_mf import basis_weights

def fun_f2py():
    basis_weights.test_memory('Y')
    return 

def fun_py():
    # It is the same function in fortran
    N = 500
    A = np.zeros((N, N), dtype=float)
    for i in range(N):
        for j in range(N):
            A[i, j] = i + 0.5*j
    B = A @ A
    return

# Define some cases
case = 1
if case == 0: fun = fun_f2py
else: fun = fun_py

tracemalloc.start()
fun()
_, memory_direct = tracemalloc.get_traced_memory()
memory_direct /= 1024*1024
print('used mem in tracemalloc', memory_direct)
tracemalloc.clear_traces()
tracemalloc.stop()

process = psutil.Process(os.getpid())
memi = process.memory_info().rss/(1024*1024)
fun()
memf = process.memory_info().rss/(1024*1024)
print('used mem in psutil', memf-memi)

start_mem = mp.memory_usage(max_usage=True)
res = mp.memory_usage(proc=(fun,), max_usage=True, interval=.0001) 
print('used mem in mprofiler', res-start_mem)

print("Peak memory (MiB):",
      getrusage(RUSAGE_SELF).ru_maxrss / 1024)