from pysrc.lib.__init__ import *

XI = np.linspace(0, 1, 1000)
ALPHA = np.logspace(-2, 1, 1000)

KAPPA1, KAPPA2 = [], []
for alpha in ALPHA:
	A = np.zeros((2, len(XI)))
	B = np.zeros((2, len(XI)))
	A[0, :] = np.pi*(1/alpha + XI)/2
	A[1, :] = 2/(np.pi*(1/alpha + XI))
	B[0, :] = 2*(1/alpha + XI)/(1 + 2/alpha)
	B[1, :] = 1/((1/alpha + XI)*np.log(1+alpha))
	KAPPA1.append(A.max()/A.min())
	KAPPA2.append(B.max()/B.min())

KAPPA1 = np.array(KAPPA1); KAPPA2 = np.array(KAPPA2)

plt.loglog(1+ALPHA, KAPPA1, label='Classic FD')
plt.loglog(1+ALPHA, KAPPA2, label='Modified FD')
plt.ylabel('Condition number '+r'$\kappa$')

# plt.loglog(1+ALPHA, KAPPA1/KAPPA2 - 1)
# plt.ylabel('Gain ' + r'$\displaystyle\frac{\kappa_{cl.}-\kappa_{mod.}}{\kappa_{mod.}}$')

plt.xlabel(r'$R_{ext}/R_{int}$')
plt.xticks([i for i in range(0, 11, 2)], [i for i in range(0, 11, 2)])
plt.minorticks_off()
plt.xlim([1, 10])
plt.ylim([1e-2, 1e3])
plt.savefig('spectral_number_QA')

N = 100
XI1, XI2 = np.meshgrid(np.linspace(0, 1, N), np.linspace(0, 1, N))
alpha = np.linspace(0, 0.5, 2*N+1)[1:-1]
beta  = np.logspace(-2, 1, N)
ALPHA, BETA = np.meshgrid(alpha, beta)

KAPPA1, KAPPA2 = [], []
for a in alpha:
	for b in beta:
		A = np.zeros((2, N, N))
		A[0, :, :] = (b**2 + a**2*(1 - 2*XI1)**2)/(b*(1 - 2*a*XI2)) 
		A[1, :, :] = (1 - 2*a*XI2)/b

		B = np.zeros((2, N, N))
		tmp = (1 - 2*a*XI2)*(b**2 + 1/3*a**2)*np.log(1 - 2*a)
		B[0, :, :] = -2*a*(b**2 + a**2*(1 - 2*XI1)**2)/((1 - 2*a*XI2)*(b**2 + 1/3*a**2)*np.log(1 - 2*a))
		B[1, :, :] = (1 - 2*a*XI2)/(1 - a)
		KAPPA1.append(A.max()/A.min())
		KAPPA2.append(B.max()/B.min())

KAPPA1 = np.reshape(np.array(KAPPA1), newshape=(len(alpha), len(beta)))
KAPPA2 = np.reshape(np.array(KAPPA2), newshape=(len(alpha), len(beta)))
DIFF   = KAPPA1/KAPPA2 - 1

fig, ax = plt.subplots(nrows=1, ncols=1)

bounds = np.linspace(0, 10, 11)
norm = mpl.colors.BoundaryNorm(boundaries=bounds, ncolors=256, extend='both')
im = ax.pcolormesh(1-2*ALPHA, BETA, DIFF.T, norm=norm, cmap='RdBu_r')
cbar = fig.colorbar(im, extend='both')
ax.grid(False)
ax.set_xlabel(r'$a/b$')
ax.set_ylabel(r'$h/b$')
# ax.set_yscale('log')
ax.set_xlim(left=0, right=1)
ax.set_ylim(bottom=0, top=5)
fig.tight_layout()
fig.savefig('spectral_number_VB')

