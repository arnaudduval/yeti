from thesis.Elliptic.__init__ import *

degree, cuts = 6, 6
fig, ax = plt.subplots(figsize=(5.5, 5.5))
colors = [COLORLIST[0], COLORLIST[1], COLORLIST[4], COLORLIST[3]]
for j, [preconditioner, color] in enumerate(zip(['WP', 'C', 'ilu', 'JMC'], colors)):
	start = time.process_time()
	problem, _, residue = simulate_el(degree, cuts, 
									preconditioner=preconditioner,
									linsolver='GMRES',
									)
	stop = time.process_time()
	print('time:%.2e, residue:%.2e'%(stop-start,residue[-1]))

	if preconditioner == 'WP': labelfig = 'w.o. preconditioner'
	elif preconditioner == 'ilu': labelfig = 'Incomplete LU'
	elif preconditioner == 'C' : labelfig = 'Standard Fast Diag.'
	elif preconditioner == 'JMC' : labelfig = 'Our contribution'
	ax.semilogy(residue, marker=MARKERLIST[j], label=labelfig, color=color)

ax.legend(ncol=2, bbox_to_anchor=(0.5, 1.2), loc='upper center')
ax.set_ylabel('Relative residue')
ax.set_xlabel('Number of iterations (GMRES)')
# ax.set_xlabel('Number of iterations (CG)')
ax.set_ylim([1e-12, 1e1])
ax.set_xlim([0, 100])
fig.tight_layout()
fig.savefig(FOLDER2SAVE+'preconditioner_el'+'.png')

# for degree in range(4, 7):
# 	for cuts in range(6, 9):
# 		start = time.process_time()
# 		residue = simulate_el(degree, cuts, preconditioner='JMC')[-1]
# 		stop = time.process_time()
# 		print('%d, %d, %.2f, %d' %(degree, cuts, stop-start, len(residue[residue>0.0])))