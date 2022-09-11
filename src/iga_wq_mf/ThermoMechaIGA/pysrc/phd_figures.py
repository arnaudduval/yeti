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
.. CASE 7: Convergence curve
"""

from lib.__init__ import *
from lib.base_functions import (create_knotvector, 
                                wq_find_positions, 
                                eval_basis_python, 
                                iga_find_positions_weights
)
from lib.create_geomdl import geomdlModel
from lib.fortran_mf_wq import fortran_mf_wq

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/phd/'
if not os.path.isdir(folder): os.mkdir(folder)

def plot_geometry2D(geo:geomdlModel):
    "Plots a 2D geometry "

    def plot_mesh(pts, shape, ax):
        "Plots mesh of control points"

        if pts.shape[0] == 3: pts = pts[:2, :]

        pts2D = []
        for j in range(shape[1]):
            pts2D_temp = []
            for i in range(shape[0]):
                pos = i + j*shape[0]
                pts2D_temp.append(pts[:, pos].tolist())
            pts2D.append(pts2D_temp)
        pts2D = np.asarray(pts2D)

        # In the first direction
        for _ in range(shape[1]): 
            x = pts2D[_, :, 0]; y = pts2D[_, :, 1]
            ax.plot(x, y, 'k--')

        # In the second direction
        for _ in range(shape[0]):
            x = pts2D[:, _, 0]; y = pts2D[:, _, 1]
            ax.plot(x, y, 'k--')

        return

    # Get geometry
    geometry = geo._geometry
    
    # Get control points
    ctrlpts = np.asarray(geo._ctrlpts)

    # Set shape
    samplesize = geo._sample_size
    shape_sample = (samplesize, samplesize)

    # Get eval points
    evalpts_old = np.array(geometry.evalpts)
    evalpts_new = []
    for j in range(samplesize):
        for i in range(samplesize):
            pos = j + i*samplesize
            evalpts_temp = evalpts_old[pos, :]
            evalpts_new.append(evalpts_temp.tolist())

    # Set values
    evalpts_new = np.asarray(evalpts_new) 
    X = np.asarray(evalpts_new[:,0].reshape(shape_sample).tolist())
    Y = np.asarray(evalpts_new[:,1].reshape(shape_sample).tolist())
    Z = np.zeros(X.shape)

    # Plot
    fig, ax = plt.subplots(nrows=1, ncols=1)
    ax.grid(None)
    ax.pcolormesh(X, Y, Z, cmap=plt.cm.Pastel1, shading='gouraud')
    ax.plot(ctrlpts[0, :], ctrlpts[1, :], 'o', label='Control points')
    plot_mesh(ctrlpts, geo._nb_ctrlpts, ax)
    ax.set_xticks(np.arange(0, max(evalpts_new[:,0])+1, 1.0))
    ax.set_yticks(np.arange(0, max(evalpts_new[:,1])+1, 1.0))

    # Set properties
    ax.axis('equal')
    ax.legend()
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    fig.tight_layout()
    
    return fig

# Set global variables
CASE = 5
extension = '.png'

if CASE == 0: # B-spline curve

    # Set filename
    filename = folder + 'BsplineCurve'+ extension

    # B-spline properties
    degree, nbel = 2, 4
    knotvector = create_knotvector(degree, nbel)

    # Create the curve 
    crv = BSpline.Curve()
    crv.degree = degree
    crv.ctrlpts = [[-1, 1, 0], [-0.5, 0.25, 0], [0, 2, 0], 
                    [0.75, -0.5, 0], [1.5, 1, 0], [2, 0, 0]]
    crv.knotvector = knotvector
    crv.delta = 0.01

    # Extract data
    pts = np.asarray(crv.evalpts)
    ctrlpts = np.asarray(crv.ctrlpts)
    x = pts[:, 0]; y = pts[:, 1]
    
    # Plot 
    fig, ax = plt.subplots(nrows=1, ncols=1)
    ax.plot(x, y, label='B-spline curve')
    ax.plot(ctrlpts[:, 0], ctrlpts[:, 1], 'o--', markersize=10, label='Control points')

    # Set properties
    ax.legend(prop={'size': 14})
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.axis('equal')
    fig.tight_layout()
    fig.savefig(filename)

elif CASE == 1: # Univariate functions

    # Set filename
    filename = folder + 'UnivariateFunctions' + extension

    # B-spline properties 
    degree, nbel = 2, 4
    knots = np.linspace(0, 1, 201)
    knotvector = create_knotvector(degree, nbel)
    qp_position = wq_find_positions(degree, knotvector, 2)
    B0, B1 = eval_basis_python(degree, knotvector, knots)
    B0 = B0.toarray(); B1 = B1.toarray()

    # Plot
    fig, [ax1, ax2] = plt.subplots(nrows=1, ncols=2, figsize=(10,4))
    for i in range(degree+nbel): 
        ax1.plot(knots, B0[i, :], linewidth=2)
        ax2.plot(knots, B1[i, :], linewidth=2)

    # Set properties
    for ax in [ax1, ax2]:
        ax.set_xlabel(r'$\xi$')
        ax.set_xticks(np.linspace(0, 1, nbel+1))
    fig.tight_layout()
    fig.savefig(filename)

elif CASE == 2: # Bivariate functions

    # Set filename
    filename = folder + 'BivariateFunctions' + extension

    # B-Spline properties
    degree, nbel = 3, 1
    knots = np.linspace(0, 1, 201)
    knotvector = create_knotvector(degree, nbel)
    B0 = eval_basis_python(degree, knotvector, knots)[0]
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
    axs[1,0].grid(None)
    axs[1,0].pcolormesh(X, Y, Z, cmap='inferno', shading = 'gouraud')
    axs[1,0].set_yticks([0, 0.5, 1])

    for i in range(degree+nbel): 
        axs[0, 0].plot(knots, B0[i, :], color="0.8")
        axs[1, 1].plot(B0[i, :], knots, color="0.8")

    axs[0,0].plot(knots, B0_plot); axs[0, 0].axis(ymin=0,ymax=1)
    axs[0,0].set_xlabel(r'$\xi$', fontsize=14)
    axs[1,1].plot(B0_plot, knots); axs[1, 1].axis(xmin=0,xmax=1)
    axs[1,1].set_ylabel(r'$\eta$', fontsize=14)
    fig.tight_layout()
    fig.savefig(filename) 

elif CASE == 3: # Quadrature points in IGA

    # Set filename
    filename = folder + 'QuadPtsIGA' + extension

    # Plot
    fig, [ax1, ax3] = plt.subplots(nrows=1, ncols=2, figsize=(8,4))

    # B-spline properties
    degree = 4
    for ax, nbel in zip([ax1, ax3], [8, 32]):
        knotvector = create_knotvector(degree, nbel)
        qp_position = iga_find_positions_weights(degree, knotvector)[0]
        XX, YY = np.meshgrid(qp_position, qp_position)
        ax.plot(XX, YY, 'ko', markersize=1.25)

        xy = np.linspace(0, 1, nbel+1)
        for i in xy:
            ax.plot([i, i], [0, 1], 'grey', linewidth=0.5, alpha=0.8)
            ax.plot([0, 1], [i, i], 'grey', linewidth=0.5, alpha=0.8)

        # Set properties
        ax.set_xticks([0, 0.5, 1])
        ax.set_yticks([0, 0.5, 1])
        ax.axis('equal')
        ax.set_ylabel(r'$\eta$')
        ax.set_xlabel(r'$\xi$')   
    fig.tight_layout()
    fig.savefig(filename) 

elif CASE == 4: # Quadrature points in WQ

    # Set filename
    filename = folder + 'QuadPtsWQ' + extension

    # Plot
    fig, [ax1, ax3] = plt.subplots(nrows=1, ncols=2, figsize=(8,4))

    # B-spline properties
    degree, r = 4, 2
    for ax, nbel in zip([ax1, ax3], [8, 32]):
        knotvector = create_knotvector(degree, nbel)
        qp_position = wq_find_positions(degree, knotvector, r)
        XX, YY = np.meshgrid(qp_position, qp_position)
        ax.plot(XX, YY, 'ko', markersize=1.25)

        xy = np.linspace(0.,1,nbel+1)
        for i in xy:
            ax.plot([i, i], [0, 1], 'grey', linewidth=0.5, alpha=0.8)
            ax.plot([0, 1], [i, i], 'grey', linewidth=0.5, alpha=0.8)

        # Set properties
        ax.set_xticks([0, 0.5, 1])
        ax.set_yticks([0, 0.5, 1])
        ax.axis('equal')
        ax.set_ylabel(r'$\eta$')
        ax.set_xlabel(r'$\xi$')
    fig.tight_layout()
    fig.savefig(filename) 

elif CASE == 5: # B-spline surface

    # Surface properties
    name = 'quarter_annulus'
    geometry = geomdlModel(name=name, **{'degree':[3, 3, 3]})
    fig = plot_geometry2D(geometry)

    # Set filename
    filename = folder + 'BsplineSurface' + extension
    fig.savefig(filename) 

elif CASE == 6: # FEM functions 

    # Set filename
    filename = folder + 'FEM_Functions' + extension

    # Functions in isoparametric space 
    x  = np.linspace(-1, 1, 90)
    f1 = x*(x-1)/2
    f2 = (1-x**2)
    f3 = x*(x+1)/2

    # Parametric space properties
    nbel  = 4
    knots = np.linspace(0, 1, nbel+1)

    # Plot
    fig, ax = plt.subplots(nrows=1, ncols=1)
    color = [next(ax._get_lines.prop_cycler)['color'] for _ in range(3)]

    for _ in range(nbel):
        x0 = knots[_]; x1 = knots[_+1]
        xtemp = x0*(1-x)/2 + x1*(1+x)/2
        ax.plot(xtemp, f1, color=color[(_)%3])
        ax.plot(xtemp, f2, color=color[(_)%3])
        ax.plot(xtemp, f3, color=color[(_)%3])
    
    # Set properties
    ax.set_xlabel(r'$\xi$', fontsize=16)
    ax.set_xticks([0, 0.5, 1])
    fig.tight_layout()
    fig.savefig(filename)

elif CASE == 7: # Convergence curve

    def power_density(P:list):
        x = P[0]
        y = P[1]
        f = np.sin(np.pi*x)*np.sin(np.pi*y)
        # f = (2*np.pi**2*np.sin(np.pi*x)*np.sin(np.pi*y)*(x**2 + y**2 - 1)*(x**2 + y**2 - 4) 
        # - 8*y**2*np.sin(np.pi*x)*np.sin(np.pi*y) 
        # - 4*np.sin(np.pi*x)*np.sin(np.pi*y)*(x**2 + y**2 - 1) 
        # - 4*np.sin(np.pi*x)*np.sin(np.pi*y)*(x**2 + y**2 - 4) 
        # - 4*x*np.pi*np.cos(np.pi*x)*np.sin(np.pi*y)*(x**2 + y**2 - 1) 
        # - 4*x*np.pi*np.cos(np.pi*x)*np.sin(np.pi*y)*(x**2 + y**2 - 4) 
        # - 4*y*np.pi*np.cos(np.pi*y)*np.sin(np.pi*x)*(x**2 + y**2 - 1) 
        # - 4*y*np.pi*np.cos(np.pi*y)*np.sin(np.pi*x)*(x**2 + y**2 - 4) 
        # - 8*x**2*np.sin(np.pi*x)*np.sin(np.pi*y)
        # )

        return f

    def solution(P:list): 
        x = P[0]
        y = P[1]
        t = 1/(2*np.pi**2)*np.sin(np.pi*x)*np.sin(np.pi*y)
        # t = np.sin(np.pi*x)*np.sin(np.pi*y)*(x**2+y**2-1)*(x**2+y**2-4)
        return t

    fig, ax = plt.subplots(nrows=1, ncols=1)
    for degree in range(3, 8):
        norm = []; ddl =[]
        for cuts in range(1, 7):
            print([degree, 2**cuts])

            blockPrint()
            # Number of elements
            nbel = 2**cuts

            # Geometry
            name = 'quarter_annulus'
            geometry = {'degree': [degree, degree, degree]}
            geometry = geomdlModel(filename= name, **geometry)
            geometry.knot_refinement(np.array([cuts, cuts, cuts]))

            # Model
            model = fortran_mf_wq(geometry)
            dof = model._thermal_dof

            # Thermal equation
            Kdd = model.eval_conductivity_matrix(indi=dof, indj=dof)
            Fd = model.eval_source_vector(power_density, indi=dof)
            Td = sp.linalg.spsolve(Kdd, Fd)  
            T = np.zeros(model._nb_ctrlpts_total)
            T[dof] = Td

            # Interpolation
            _, qp_PS, _, u_interp = model.interpolate_field(T)
            u_exact = [solution(qp_PS[:, i][0]) for i in range(len(u_interp))]
            u_exact = np.array(u_exact)

            # Error
            error = np.linalg.norm(u_exact - u_interp, np.inf)/np.linalg.norm(u_exact, np.inf)

            norm.append(error*100)
            ddl.append(nbel)
            enablePrint()

        norm = np.asarray(norm)
        ddl = np.asarray(ddl)

        # Figure 
        ax.loglog(ddl, norm, label='p = ' + str(degree))
        
        # Get slope
        slope, _ = np.polyfit(np.log10(ddl[1:5]),np.log10(norm[1:5]), 1)
        slope = round(slope)
        annotation.slope_marker((ddl[3], norm[3]), slope, 
                                text_kwargs={'fontsize': 14},
                                poly_kwargs={'facecolor': (0.73, 0.8, 1)})

    # Set filename
    filename = folder + 'ConvergenceIGA_Annulus2'+ extension

    # Properties
    plt.grid()
    plt.xlabel("Number of elements $nb_{el}$", fontsize= 16)
    plt.ylabel("Relative error (%)", fontsize= 16)   
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.legend(prop={'size': 14})
    plt.tight_layout()
    plt.savefig(filename)

else: raise Warning('Case unkwnon')

