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
			'#1F77B4', '#FF7F0E', 
			'#2CA02C', '#D62728', 
			'#9467BD', '#8C564B', 
			'#E377C2', '#7F7F7F', 
			'#BCBD22', '#17BECF',
			'#AEC7E8', '#FFBB78',
			'#98DF8A', '#FF9896',
			'#C5B0D5', '#C49C94',
			'#F7B6D2', '#C7C7C7',
			'#DBDB8D', '#9EDAE5',
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

CONFIGLINE0  = {'marker': 's', 'linestyle': '-', 'markersize': 10}
CONFIGLINE1 = {'marker': 'o', 'linestyle': '--', 'markersize': 6}
CONFIGLINE2 = {'marker': 'x', 'linestyle': ':', 'markersize': 6}
CONFIGLINE4 = {'marker': 'd', 'linestyle': '-.', 'markersize': 6}

import pyvista as pv
pv.rcParams['transparent_background'] = True
pv.global_theme.font.family = 'times'
pv.global_theme.font.color  = 'black'

# Define global functions
def blockPrint(): sys.stdout = open(os.devnull, 'w')
def enablePrint(): sys.stdout = sys.__stdout__