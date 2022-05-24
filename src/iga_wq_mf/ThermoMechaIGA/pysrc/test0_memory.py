"""
.. Test of memory
.. We test which library is better to get memory used
.. Joaquin Cornejo 
"""

# Python libraries
import tracemalloc
import os, psutil
import numpy as np

# My libraries
from iga_wq_mf import basis_weights

tracemalloc.start()
basis_weights.test_memory('Y')
_, memory_direct = tracemalloc.get_traced_memory()
memory_direct /= 1024*1024
print(memory_direct)
tracemalloc.clear_traces()
tracemalloc.stop()

process = psutil.Process(os.getpid())
memi = process.memory_info().rss/(1024*1024)
vec = np.ones(100000)
basis_weights.test_memory('Y')
memf = process.memory_info().rss/(1024*1024)
print(memf-memi)

# ===========================
import memory_profiler as mp
def fun():
    basis_weights.test_memory('Y')
    return "XXXXX"
start_mem = mp.memory_usage(max_usage=True)
res = mp.memory_usage(proc=(fun,), max_usage=True, interval=.001) 
print('start mem', start_mem)
print('max mem', res)
print('used mem', res-start_mem)
