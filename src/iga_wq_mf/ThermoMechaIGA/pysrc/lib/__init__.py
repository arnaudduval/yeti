"""
Author: Joaquin Cornejo
This files contains the libraries used in other files.
"""

# Python libraries
import os, sys, time
from copy import deepcopy
import numpy as np, math, statistics, pandas as pd
from scipy import sparse as sp, linalg as sclin
from geomdl import (
    helpers,
    fitting, 
    BSpline, 
    operations,
)
from pyevtk.hl import gridToVTK
import matplotlib as mpl
from matplotlib import pyplot as plt
from mpltools import annotation

# Default properties 
mpl.rcParams['figure.figsize'] = (5.0, 4.0)    
mpl.rcParams['axes.unicode_minus'] = False
mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['font.family'] = 'STIXGeneral'
mpl.rcParams["axes.grid"] =True

# Define size
SMALL_SIZE, MEDIUM_SIZE, BIGGER_SIZE = 14, 16, 18
plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

# YETI libraries
from iga_wq_mf import basis_weights, assembly, solver
from preprocessing.igaparametrization import IGAparametrization

# Define global functions
def blockPrint():
    sys.stdout = open(os.devnull, 'w')

def enablePrint():
    sys.stdout = sys.__stdout__