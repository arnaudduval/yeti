"""
.. Isogeometric analysis figures
.. Select case: 
.. CASE 0: B-spline curve
.. CASE 1: Univariate B-spline functions in parametric space
.. CASE 2: Bivariate B-spline functions in parametric space
.. CASE 3: Quadrature points in IGA-Galerkin approach
.. CASE 4: Quadrature points in IGA-WQ approach
.. CASE 5: B-spline surface
.. CASE 6: FEM basis 
"""

# Python libraries
import os
from geomdl import BSpline
import numpy as np
import matplotlib.pyplot as plt

# My libraries
from lib.base_functions import (create_knotvector, 
                                wq_find_positions, 
                                eval_basis_python, 
                                iga_find_positions_weights
)
from lib.create_geomdl import geomdlModel

# Choose folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/phd/'
if not os.path.isdir(folder):
    os.mkdir(folder)

CB_color_cycle = ['#377eb8', '#ff7f00', '#4daf4a',
                '#f781bf', '#a65628', '#984ea3',
                '#999999', '#e41a1c', '#dede00']

# Select extension
extension = '.png'

# Select case
CASE = 4

if CASE == 0: # B-spline curve

    # B-spline properties
    degree = 2
    nbel = 4
    knotvector = create_knotvector(degree, nbel)

    # Create the curve 
    crv = BSpline.Curve()
    crv.degree = degree
    crv.ctrlpts = [[-1, 1, 0], [-0.5, 0.25, 0], [0, 2, 0], 
                    [0.75, -0.5, 0], [1.5, 1, 0], [2, 0, 0]]
    crv.knotvector = knotvector
    crv.delta = 0.01

    # Get data points
    pts = np.asarray(crv.evalpts)

    # Get x and y
    x = pts[:, 0]
    y = pts[:, 1]
    ctrlpts = np.asarray(crv.ctrlpts)

    # Plot control points
    plt.figure()
    for _ in range(6):
        plt.plot(ctrlpts[_, 0], ctrlpts[_, 1], 'o', markersize=15)
    plt.plot(x, y, label='B-spline curve')
    plt.plot(ctrlpts[:, 0], ctrlpts[:, 1], '--')

    # Set filename
    filename = folder + 'BsplineCurve'+ extension

    # Properties
    plt.grid()
    plt.xlabel('X', fontsize=16)
    plt.ylabel('Y', fontsize=16)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.axis('equal')
    plt.legend(prop={'size': 14})
    plt.tight_layout()
    plt.savefig(filename)

elif CASE == 1: # Univariate functions

    # B-spline properties 
    degree = 2
    nbel = 4
    knots = np.linspace(0, 1, 181)
    knotvector = create_knotvector(degree, nbel)
    qp_wq = wq_find_positions(degree, knotvector, 2)
    B0, B1 = eval_basis_python(degree, knotvector, knots)
    B0 = B0.toarray()
    B1 = B1.toarray()

    fig, [ax1, ax2] = plt.subplots(nrows=1, ncols=2, figsize=(10,4))
    for i in range(degree+nbel): 
        ax1.plot(knots, B0[i, :], linewidth=1, color= CB_color_cycle[i])
        ax2.plot(knots, B1[i, :], linewidth=1, color= CB_color_cycle[i])

    # Set filename
    filename = folder + 'basis_p' + str(degree) + '_nbel' + str(nbel) + extension

    # Properties
    for ax in [ax1, ax2]:
        ax.grid()
        ax.set_xlabel(r'$\xi$', fontsize=16)
        ax.set_xticks(np.linspace(0, 1, nbel+1))
        ax.tick_params(axis='x', labelsize=16)
        ax.tick_params(axis='y', labelsize=16)

    plt.tight_layout()
    plt.savefig(filename)

elif CASE == 2: # Bivariate functions

    # B-Spline properties
    degree = 2
    nbel = 1
    knots = np.linspace(0, 1, 181)
    knotvector = create_knotvector(degree, nbel)
    B0, _ = eval_basis_python(degree, knotvector, knots)
    B0 = B0.toarray()
    B0_plot = B0[1, :]

    # B-Spline 2D
    X, Y = np.meshgrid(knots, knots)
    Z = np.kron(B0_plot, B0_plot).reshape((len(knots), len(knots)))

    # Plot
    fig, axs = plt.subplots(2, 2, sharex="col", sharey="row",
                             gridspec_kw=dict(height_ratios=[1,3],
                                              width_ratios=[3,1]))
    axs[0,1].set_visible(False)
    axs[0,0].set_box_aspect(1/3)
    axs[1,0].set_box_aspect(1)
    axs[1,1].set_box_aspect(3/1)
    axs[1,0].pcolormesh(X, Y, Z, cmap='inferno', shading = 'gouraud')
    axs[1,0].tick_params(axis='x', labelsize=14)
    axs[1,0].tick_params(axis='y', labelsize=14)

    for i in range(degree+nbel): 
        axs[0, 0].plot(knots, B0[i, :], color="0.8")
        axs[1, 1].plot(B0[i, :], knots, color="0.8")

    axs[0,0].plot(knots, B0_plot); axs[0, 0].axis(ymin=0,ymax=1)
    axs[0,0].set_xlabel(r'$\xi$', fontsize=14)
    axs[0,0].tick_params(axis='y', labelsize=14)
    axs[1,1].plot(B0_plot, knots); axs[1, 1].axis(xmin=0,xmax=1)
    axs[1,1].set_ylabel(r'$\eta$', fontsize=14)
    axs[1,1].tick_params(axis='x', labelsize=14)

    # Set filename
    filename = folder + 'Basis_2D' + extension

    plt.tight_layout()
    plt.savefig(filename) 

elif CASE == 3: # Quadrature points in IGA

    fig, [ax1, ax2, ax3] = plt.subplots(nrows=1, ncols=3, figsize=(12,4))

    # B-spline properties
    degree = 4
    for ax, nbel in zip([ax1, ax2, ax3], [8, 16, 32]):
        knotvector = create_knotvector(degree, nbel)
        wq, _ = iga_find_positions_weights(degree, knotvector)
        Xwq, Ywq = np.meshgrid(wq, wq)
        ax.plot(Xwq, Ywq, 'ko', markersize=1)

        # Properties
        ax.grid()
        ax.set_xticks([0, 0.5, 1])
        ax.set_yticks([0, 0.5, 1])
        ax.axis('equal')
        ax.set_ylabel(r'$\eta$', fontsize=14)
        ax.set_xlabel(r'$\xi$', fontsize=14)
        ax.tick_params(axis='x', labelsize=14)
        ax.tick_params(axis='y', labelsize=14)

    # Set filename
    filename = folder + 'QuadPtsIGA'+ str(nbel) + extension

    plt.tight_layout()
    plt.savefig(filename) 

elif CASE == 4: # Quadrature points in WQ
    # B-spline properties
    r = 2 
    degree = 4

    fig, [ax1, ax2, ax3] = plt.subplots(nrows=1, ncols=3, figsize=(12,4))

    # B-spline properties
    degree = 4
    for ax, nbel in zip([ax1, ax2, ax3], [8, 16, 32]):
        knotvector = create_knotvector(degree, nbel)
        wq = wq_find_positions(degree, knotvector, r)
        Xwq, Ywq = np.meshgrid(wq, wq)
        ax.plot(Xwq, Ywq, 'ko', markersize=1)

        # Properties
        ax.grid()
        ax.set_xticks([0, 0.5, 1])
        ax.set_yticks([0, 0.5, 1])
        ax.axis('equal')
        ax.set_ylabel(r'$\eta$', fontsize=14)
        ax.set_xlabel(r'$\xi$', fontsize=14)
        ax.tick_params(axis='x', labelsize=14)
        ax.tick_params(axis='y', labelsize=14)

    # Set filename
    filename = folder + 'QuadPtsWQ'+str(nbel) + extension

    plt.tight_layout()
    plt.savefig(filename) 

elif CASE == 5: # B-spline surface

    # Surface properties
    name = 'quarter_annulus'
    geometry = {'degree': [3, 3]}
    geometry = geomdlModel(filename= name, **geometry)
    fig = geometry.plot_2D_geometry()

    # Set filename
    filename = folder + 'BsplineSurface' + extension
    fig.savefig(filename) 

elif CASE == 6: # FEM functions $

    # Functions in isoparametric space 
    x = np.linspace(-1, 1, 90)
    f1 = x*(x-1)/2
    f2 = (1-x**2)
    f3 = x*(x+1)/2

    # Parametric space properties
    nbel = 4
    knots = np.linspace(0, 1, nbel+1)

    ax = plt.gca()
    color = []
    for _ in range(3): 
        color.append(next(ax._get_lines.prop_cycler)['color'])

    for _ in range(nbel):
        x0 = knots[_]
        x1 = knots[_+1]
        xtemp = x0*(1-x)/2 + x1*(1+x)/2
        plt.figure(1)
        plt.plot(xtemp, f1, color=color[(_)%3])
        plt.plot(xtemp, f2, color=color[(_)%3])
        plt.plot(xtemp, f3, color=color[(_)%3])
    
    # Properties
    plt.grid()
    plt.gca().set_xlabel(r'$\xi$', fontsize=16)
    plt.xticks([0, 0.5, 1])
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    
    # Set filename
    filename = folder + 'FEM' + '_nbel' + str(nbel) + extension

    plt.tight_layout()
    plt.savefig(filename)

else: 
    raise Warning('Case unkwnon')
