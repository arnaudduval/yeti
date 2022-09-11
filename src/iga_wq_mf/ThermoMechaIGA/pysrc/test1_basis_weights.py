"""
.. Test of basis and weights 
.. We test if functions done in python and fortran for WQ approach 
.. gives the expected results.
.. Joaquin Cornejo 
"""
from lib.__init__ import *
from lib.base_functions import (create_knotvector, 
                                eval_basis_python,
                                iga_find_positions_weights,
                                wq_find_basis_weights_opt, 
                                wq_find_basis_weights_fortran
)

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/test1/'
if not os.path.isdir(folder): os.mkdir(folder)

# Set number of elements
nbel_list = [2**i for i in np.arange(1, 6)]

for varName in ['I00', 'I01', 'I10', 'I11']:
    
    # Create plot
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8,4))

    for degree in range(3, 8):

        # Initialize
        norm_fortran = []; norm_python = []
        color = next(ax._get_lines.prop_cycler)['color']

        for nbel in nbel_list: 

            # Initialize
            knotvector = create_knotvector(degree, nbel)
            nb_ctrlpts = degree + nbel

            # ----------
            # FORTRAN
            # ---------- 
            B, W, indi, indj = wq_find_basis_weights_fortran(degree, knotvector)[2:]
            nb_qp_wq = np.max(indj); indi -= 1; indj -= 1

            # Create basis and weights from fortran
            B0f  = sp.csr_matrix((B[:,0], indj, indi), shape=(nb_ctrlpts, nb_qp_wq))
            B1f  = sp.csr_matrix((B[:,1], indj, indi), shape=(nb_ctrlpts, nb_qp_wq))
            W00f = sp.csr_matrix((W[:,0], indj, indi), shape=(nb_ctrlpts, nb_qp_wq))
            W01f = sp.csr_matrix((W[:,1], indj, indi), shape=(nb_ctrlpts, nb_qp_wq))
            W10f = sp.csr_matrix((W[:,2], indj, indi), shape=(nb_ctrlpts, nb_qp_wq))
            W11f = sp.csr_matrix((W[:,3], indj, indi), shape=(nb_ctrlpts, nb_qp_wq))

            # Calculate I
            I00f = W00f @ B0f.T; I01f = W01f @ B1f.T
            I10f = W10f @ B0f.T; I11f = W11f @ B1f.T

            # ---------
            # PYTHON
            # ---------
            B0p, B1p, W00p, W01p, W10p, W11p = wq_find_basis_weights_opt(degree, knotvector, 2)[1:]

            # Calculate I
            I00p = W00p @ B0p.T; I01p = W01p @ B1p.T
            I10p = W10p @ B0p.T; I11p = W11p @ B1p.T

            # -------------
            # REFERENCE
            # -------------
            qp_cgg, Wcgg = iga_find_positions_weights(degree, knotvector)
            B0, B1 = eval_basis_python(degree, knotvector, qp_cgg)

            # Calculate I
            I00 = B0 @ np.diag(Wcgg) @ B0.T; I01 = B0 @ np.diag(Wcgg) @ B1.T
            I10 = B1 @ np.diag(Wcgg) @ B0.T; I11 = B1 @ np.diag(Wcgg) @ B1.T

            # To choose variables
            if varName   == 'I00': var1 = I00; var2 = I00f; var3 = I00p
            elif varName == 'I01': var1 = I01; var2 = I01f; var3 = I01p
            elif varName == 'I10': var1 = I10; var2 = I10f; var3 = I10p
            elif varName == 'I11': var1 = I11; var2 = I11f; var3 = I11p

            # Compare results 
            error_fortran = var1 - var2
            norm_temp = np.linalg.norm(error_fortran, np.inf)/np.linalg.norm(var1, np.inf)*100
            if norm_temp > 1e-5: raise Warning("Something happend. Fortran basis are wrong")
            norm_fortran.append(norm_temp)

            error_python = var1 - var3
            norm_temp = np.linalg.norm(error_python, np.inf)/np.linalg.norm(var1, np.inf)*100
            if norm_temp > 1e-5: raise Warning("Something happend. Python basis are wrong")
            norm_python.append(norm_temp)

        # Figure 
        strlabel = 'Degree p = ' + str(degree)
        ax.loglog(nbel_list, norm_fortran, '-o', label=strlabel, color=color)
        ax.loglog(nbel_list, norm_python, '--P', color=color)

    # Plot configurations
    ax.grid()
    ax.set_xlabel("Number of elements $nb_{el}$")
    ax.set_ylabel("Relative error (%)")
    ax.legend(bbox_to_anchor= (1.05, 1.0), loc='upper left')
    fig.tight_layout()
    fig.savefig(folder + 'Error_basisweights_' + varName +'.png')
