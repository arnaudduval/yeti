"""
Author: Joaquin Cornejo
This files contains the libraries used in other files.
"""

# Python libraries
import os, sys, time
from copy import deepcopy
import numpy as np, math, statistics
from scipy import sparse as sp, linalg as sclin
from geomdl import (
    helpers,
    fitting, 
    BSpline, 
    operations,
)
from pyevtk.hl import gridToVTK
import matplotlib.pyplot as plt

# YETI libraries
from iga_wq_mf import basis_weights, assembly, solver
from preprocessing.igaparametrization import IGAparametrization

# Define global functions
def blockPrint():
    sys.stdout = open(os.devnull, 'w')

def enablePrint():
    sys.stdout = sys.__stdout__