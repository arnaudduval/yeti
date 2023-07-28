import numpy as np
from matplotlib import pyplot as plt

Tx = 10
R  = 10

Np = int(1e4)
theta, radius = np.meshgrid(np.linspace(0.0, np.pi/2, Np), np.linspace(R, 30*R, Np))
xx = radius * np.cos(theta)
yy = radius * np.sin(theta)

r_square = xx**2 + yy**2
theta_times2 = 2.0*theta
div = R**2/r_square
PolarSt = np.zeros((2, 2, Np, Np))
PolarSt[0, 0, :, :] = Tx/2.0*(1.0 - div + np.cos(theta_times2)*(1.0 - 4.0*div + 3.0*div**2))
PolarSt[1, 1, :, :] = Tx/2.0*(1.0 + div - np.cos(theta_times2)*(1.0 + 3.0*div**2))
PolarSt[0, 1, :, :] = PolarSt[1, 0, :, :] = -Tx/2.0*(1.0 + 2*div - 3*div**2)*np.sin(theta_times2)
Rot = np.zeros((2, 2, Np, Np))
Rot[0, 0, :, :] = Rot[1, 1, :, :] = np.cos(theta_times2/2.0)
Rot[0, 1, :, :] = np.sin(theta_times2/2.0); Rot[1, 0, :, :] = -Rot[0, 1, :, :]
tmp = np.einsum('iklm,kjlm->ijlm', PolarSt, Rot)
CartSt = np.einsum('kilm,kjlm->ijlm', Rot, tmp)
VMSt = np.sqrt(np.einsum('ijkl,ijkl->kl', CartSt, CartSt))

fig, ax = plt.subplots()
plot = ax.pcolormesh(xx, yy, VMSt, cmap='viridis', shading='auto')
cbar = fig.colorbar(plot)
ax.set_xlim([0, 20])
ax.set_ylim([0, 20])
ax.set_aspect('equal')
fig.savefig('Plot2D.png')