from thesis.StateOfArt.__init__ import *

FIGCASE = 1

if FIGCASE == 0:

	# Set filename
	filenameDat = FOLDER2FIND + 'AssemblyWQ' + '.dat'
	filenameFig = FOLDER2SAVE + 'WeightedQuadrature' + '.pdf' 

	fig, ax = plt.subplots(figsize=(6, 4))

	# Load data
	tabLiterature = pd.read_table(filenameDat, sep='\t', names=['degree', 'wq', 'iga'])
	ax.semilogy(tabLiterature.degree, tabLiterature.iga, 's-', label='Gauss quadrature')
	ax.semilogy(tabLiterature.degree, tabLiterature.wq, 'o--', label='Weighted quadrature')

	ax.legend()
	ax.set_xlim([2, 10])
	ax.set_ylim([1e0, 1e6])
	ax.set_ylabel('Setup time (s)')
	ax.set_xlabel('Polynomial degree')
	fig.tight_layout()
	fig.savefig(filenameFig)

elif FIGCASE == 1:

	IgaPlot = {'marker': 's', 'linestyle': '-', 'markersize': 10}
	degList = np.arange(1, 11)
	filenameDat_error = FOLDER2FIND + 'PoissonError' + '.dat'
	filenameDat_time  = FOLDER2FIND + 'PoissonTime' + '.dat'

	fig, ax = plt.subplots(figsize=(6,4))
	cmap = plt.get_cmap('coolwarm', 10)

	Elist = np.loadtxt(filenameDat_error)
	Tlist = np.loadtxt(filenameDat_time)
	for pos in range(np.size(Elist, axis=1)):
		im = ax.scatter(Tlist[:len(degList), pos], Elist[:len(degList), pos], c=degList,
						cmap=cmap, marker=IgaPlot['marker'], s=10*IgaPlot['markersize'])
			
		ax.loglog(Tlist[:len(degList), pos], Elist[:len(degList), pos], 
				color='k', marker='', linestyle=IgaPlot['linestyle'], alpha=0.5)
		
		if pos>0:
			ax.text(Tlist[-1, pos]/2, Elist[-1, pos]/20, str(int(2**(pos+4)))+r'$^3$'+' el.')
		else:
			ax.text(Tlist[int(len(degList)/2), pos]/2, Elist[-1, pos]/20, str(int(2**(pos+4)))+r'$^3$'+' el.')


	cbar = plt.colorbar(im)
	cbar.set_label('Degree')
	tick_locs = 1+(np.arange(len(degList)) + 0.5)*(len(degList)-1)/len(degList)
	cbar.set_ticks(tick_locs)
	cbar.set_ticklabels(degList)

	ax.grid(False)
	ax.set_ylim(top=1e0, bottom=1e-15)
	ax.set_xlim(left=1e-1, right=1e4)
	ax.set_ylabel('Relative '+r'$H^1$'+' error')
	ax.set_xlabel('Computation time (s)')
	fig.tight_layout()
	fig.savefig(FOLDER2SAVE + 'MatrixFree' +  '.pdf')
