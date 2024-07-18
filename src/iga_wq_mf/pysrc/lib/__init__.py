"""
	Author: Joaquin Cornejo
	This files contains the libraries used in other files.
"""

# YETI libraries
from iga_wq_mf import basisweights, geophy, heatsolver, plasticitysolver, eigensolver, stheatsolver
from preprocessing.igaparametrization import IGAparametrization

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

import matplotlib as mpl
from matplotlib import pyplot as plt
from mpltools import annotation
from cycler import cycler

# Default properties 
MARKERLIST = ['o', 'v', 's', 'X', '+', 'p', '*']
# COLORLIST  = ['#377eb8', '#ff7f00', '#4daf4a',
# 			'#f781bf', '#a65628', '#984ea3',
# 			'#999999', '#e41a1c', '#dede00']
COLORLIST = [
			'#1F77B4', '#AEC7E8', '#FF7F0E', '#FFBB78',
			'#2CA02C', '#98DF8A', '#D62728', '#FF9896', 
			'#9467BD', '#C5B0D5', '#8C564B', '#C49C94',
			'#E377C2', '#F7B6D2', '#7F7F7F', '#C7C7C7',
			'#BCBD22', '#DBDB8D', '#17BECF', '#9EDAE5',
			]

mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}'] 
mpl.rcParams.update({'figure.autolayout': True})
mpl.rcParams['figure.figsize'] = (6.0, 4.0) 
mpl.rcParams['figure.dpi'] = 300 
mpl.rcParams['axes.unicode_minus'] = True
mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['font.family'] = 'STIXGeneral'
mpl.rcParams["axes.grid"] = True
mpl.rcParams['axes.prop_cycle'] = cycler('color', COLORLIST)

# Define size
SMALL_SIZE, MEDIUM_SIZE, BIGGER_SIZE = 14, 16, 18
plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

import pyvista as pv
pv.rcParams['transparent_background'] = True
pv.global_theme.font.family = 'times'
pv.global_theme.font.color  = 'black'

# Define global functions
def blockPrint(): sys.stdout = open(os.devnull, 'w')
def enablePrint(): sys.stdout = sys.__stdout__