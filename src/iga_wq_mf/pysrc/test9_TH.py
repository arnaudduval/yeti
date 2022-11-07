"""
.. Test of transient heat solver
.. ATTENTION: IT ONLY WORKS WITH 'ISOTROPIC' MATERIALS
.. Joaquin Cornejo 
"""

from lib.__init__ import *
from lib.physics import *
from lib.create_geomdl import geomdlModel
from lib.fortran_mf_wq import fortran_mf_wq
from lib.base_functions import create_table_properties, sigmoid

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/test8/'
if not os.path.isdir(folder): os.mkdir(folder)

# Set global variables
dataExist   = True
geolist     = ['VB', 'CB']
method_list = ['WP', 'C', 'JMC']

if not dataExist:

    degree, cuts = 6, 5
    conductivity, capacity = 0.1, 1.0
    theta = 1.0
    time_list = np.linspace(0, 25, 81)  
    table_Kprop = create_table_properties(setKprop, prop=conductivity)
    table_Cprop = create_table_properties(setCprop, prop=capacity)     

    for geoname in geolist:
        for PCGmethod in method_list:
            filename = folder + 'ResPCG_' + geoname + '_' + PCGmethod + '.dat'        

            # Create model 
            geometry = {'degree':[degree, degree, degree]}
            modelGeo = geomdlModel(geoname, **geometry)
            modelIGA = modelGeo.export_IGAparametrization(nb_refinementByDirection=
                                                        np.array([cuts, cuts, cuts]))
            modelPhy = fortran_mf_wq(modelIGA)

            # Add material 
            material = {'capacity':capacity, 'conductivity':conductivity*np.eye(3)}
            modelPhy._set_material(material)

            # Block boundaries
            Dirichlet = {'thermal':np.array([[1, 1], [0, 0], [0, 0]])}
            modelPhy._set_dirichlet_boundaries(Dirichlet)

            # Add constant temperature
            modelPhy._add_thermal_IBC(np.array([[0, 1], [0, 0], [0, 0]]), temperature=1.0)

            # ---------------------
            # Transient model
            # ---------------------
            # Interpolate temperature on boundaries over time 
            GBound = np.zeros((len(modelPhy._thermal_dod), len(time_list)))
            for i in range(len(time_list)): GBound[:, i] = modelPhy._get_thermal_IBC()

            # Add external force (transient)
            Fend  = modelPhy.eval_source_vector(powden)
            Fendt = np.atleast_2d(Fend).reshape(-1, 1)
            Fext  = np.kron(Fendt, sigmoid(time_list))

            # Solve
            Tsol, resPCG = modelPhy.MFtransientHeatNL(F=Fext, G=GBound, time_list=time_list,
                                            table_Kprop=table_Kprop, table_Cprop=table_Cprop, 
                                            methodPCG=PCGmethod, theta=theta)
            print('Finish')
            np.savetxt(filename, resPCG)
            # modelPhy.export_results(u_ctrlpts=Tsol[:, -1], folder=folder, nbDOF=1)

else:
    # for geoname in geolist:
    #     for PCGmethod in mlist:
    #         filename = folder + 'ResPCG_' + geoname + '_' + PCGmethod + '.dat'
            
    #         # --------------
    #         # Post-treatment
    #         # --------------
    #         resPCG = np.loadtxt(filename)
    #         resPCG = resPCG[:, resPCG[0, :]>0]

    #         # Colors
    #         colorset = ['#377eb8', '#ff7f00', '#4daf4a',
    #                     '#f781bf', '#a65628', '#984ea3',
    #                     '#999999', '#e41a1c', '#dede00']

    #         fig, [ax1, ax2] = plt.subplots(nrows=1, ncols=2, figsize=(12, 4), 
    #                                     gridspec_kw=dict(width_ratios=[1,2]))
    #         # Initialize new data
    #         step_list, iterNl_list, niter_list = [], [], []
    #         for _ in range(np.shape(resPCG)[1]):
    #             step = resPCG[0, _]; iterNL = resPCG[1, _]
    #             newresidue = resPCG[2:, _]; newresidue = newresidue[newresidue>0]*100
    #             step_list.append(int(step)) # Maybe it will no be used
    #             iterNl_list.append(int(iterNL))
    #             niter_list.append(len(newresidue))
    #             ax1.semilogy(np.arange(len(newresidue)), newresidue, 
    #                         color=colorset[int(step%len(colorset))], alpha=1.0/iterNL)
    #         # Print the first
    #         step = resPCG[0, 0]; iterNL = resPCG[1, 0]
    #         newresidue = resPCG[2:, 0]; newresidue = newresidue[newresidue>0]*100
    #         ax1.semilogy(np.arange(len(newresidue)), newresidue, 'o-', linewidth=2.5,
    #                     color='black', label='First step NR1')
    #         ax1.legend(loc=1)
    #         ax1.set_xlabel('Number of iterations')
    #         ax1.set_ylabel('Relative residue ' + r'$\displaystyle\frac{||r||_\infty}{||b||_\infty}$' + ' (\%)')

    #         # ----------------------
    #         # Get labels of bar-diagram
    #         labels_NL   = np.unique(iterNl_list)
    #         labels_iter = np.unique(niter_list)
            
    #         # Count data
    #         Data2plot = np.zeros((len(labels_iter), len(labels_NL)+1), dtype=int)
    #         Data2plot[:, 0] = np.array(labels_iter) - 1
    #         for i, iter in enumerate(labels_iter):
    #             for j, NL in enumerate(labels_NL):
    #                 c = 0
    #                 for k in range(len(step_list)):
    #                     if iterNl_list[k] == NL and niter_list[k] == iter:
    #                         c += 1
    #                 Data2plot[i, j+1] = c
            
    #         sumData = np.sum(Data2plot, axis=0)
    #         for i in range(1, np.shape(Data2plot)[1]):
    #             Data2plot[:, i] = Data2plot[:, i]/sumData[i]*100

    #         # Set columns
    #         columns = ['Iterations']
    #         columns.extend(['NR' + str(NL) + ': ' + str(nb) for NL, nb in zip(labels_NL, sumData[1:])])

    #         # Set dataframe
    #         df = pd.DataFrame(Data2plot, columns=columns)

    #         # Plot
    #         df.plot(x='Iterations', kind='bar', ax=ax2)
    #         for container in ax2.containers:
    #             ax2.bar_label(container, fontsize=10)
    #         ax2.grid(None)
    #         ax2.set_ylim(top=np.ceil(ax2.get_ylim()[-1]/10+1)*10)
    #         ax2.set_ylabel('Frequency (\%)')
    #         ax2.legend(loc=2)

    #         # Set properties
    #         filename = folder + 'TransientNL_' + geoname + '_' + PCGmethod + '.png'
    #         fig.tight_layout()
    #         fig.savefig(filename)

    for geoname in geolist:
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(6, 4))

        for PCGmethod in method_list:
            filename = folder + 'ResPCG_' + geoname + '_' + PCGmethod + '.dat'
            resPCG = np.loadtxt(filename)
            resPCG = resPCG[:, resPCG[0, :]>0]

            if PCGmethod == "WP": labelmethod='w.o. preconditioner'
            elif PCGmethod == "C": labelmethod='Classic FD method'
            elif PCGmethod == "JMC": labelmethod='This work'

            # # Initialize new data
            # step_list, iterNl_list, niter_list = [], [], []
            # for _ in range(np.shape(resPCG)[1]):
            #     step = resPCG[0, _]; iterNL = resPCG[1, _]
            #     newresidue = resPCG[2:, _]; newresidue = newresidue[newresidue>0]
            #     step_list.append(int(step)) # Maybe it will no be used
            #     iterNl_list.append(int(iterNL))
            #     niter_list.append(len(newresidue))
            #     ax.semilogy(np.arange(len(newresidue)), newresidue, 
            #                 color=colorSet[int(step%len(colorSet))], alpha=1.0/iterNL)

            # Print the first
            step = resPCG[0, 0]; iterNL = resPCG[1, 0]
            newresidue = resPCG[2:, 0]; newresidue = newresidue[newresidue>0]
            ax.semilogy(np.arange(len(newresidue)), newresidue, 'o-', linewidth=2.5,
                        label=labelmethod)
            # ax.legend(loc=0)
            ax.set_xlabel('Number of iterations of BiCGSTAB solver')
            ax.set_ylabel('Relative residue ' + r'$\displaystyle\frac{||r||_\infty}{||b||_\infty}$')
            ax.set_ybound(lower=1e-12, upper=10)

        filename = folder + 'TransientNL_' + geoname + '.pdf'
        fig.tight_layout()
        fig.savefig(filename)