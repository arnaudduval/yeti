from pysrc.lib.__init__ import *

# Select folder
full_path = os.path.realpath(__file__)
folder2save = os.path.dirname(full_path) + '/results/biblio/'
folder2find = os.path.dirname(full_path) + '/data/'

FIGCASE = 0

if FIGCASE == 0:

	# Set filename
	filename_data   = folder2find + 'PoissonAssembly' + '.dat'
	filename_figure = folder2save + 'WeightedQuadrature' + '.pdf' 

	fig, ax = plt.subplots(figsize=(6,4))

	# Load data
	litterature = pd.read_table(filename_data, sep='\t', names=['degree', 'wq', 'iga'])
	ax.semilogy(litterature.degree, litterature.wq, 'o--', label='Weighted quadrature')
	ax.semilogy(litterature.degree, litterature.iga, 'o--', label='Gauss quadrature')

	ax.legend()
	ax.set_xlim([2, 10])
	ax.set_ylim([1e0, 1e6])
	ax.set_ylabel('Setup time (s)')
	ax.set_xlabel('Polynomial degree')
	fig.tight_layout()
	fig.savefig(filename_figure)

elif FIGCASE == 1:

	IgaPlot = {'marker': 's', 'linestyle': '-', 'markersize': 10}
	degList = np.arange(1, 11)
	filenameError = folder2find + 'PoissonError' + '.dat'
	filenameTime  = folder2find + 'PoissonTime' + '.dat'

	fig, ax = plt.subplots(figsize=(6,4))
	cmap = plt.get_cmap('coolwarm', 10)

	Elist = np.loadtxt(filenameError)
	Tlist = np.loadtxt(filenameTime)
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
	fig.savefig(folder2save + 'MatrixFree' +  '.pdf')
