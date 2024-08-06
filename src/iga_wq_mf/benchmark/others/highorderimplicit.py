import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

def dynamics(u, A_inv, F, K):
    return A_inv @ (F - K @ u)

def gauss_legendre_4th_order_step(u, dynamics, dt, A_inv, F, K):
    sq3 = np.sqrt(3)
    c1 = 1/2 - sq3/6
    c2 = 1/2 + sq3/6

    def implicit_system(stages):
        k1, k2 = stages[:len(u)], stages[len(u):]
        eq1 = k1 - dynamics(u + dt * (c1 * k1 + (1/4 - sq3/6) * k2), A_inv, F, K)
        eq2 = k2 - dynamics(u + dt * ((1/4 + sq3/6) * k1 + c2 * k2), A_inv, F, K)
        return np.concatenate([eq1, eq2])

    initial_guess = np.tile(dynamics(u, A_inv, F, K), 2)
    k1, k2 = fsolve(implicit_system, initial_guess).reshape(2, -1)
    
    u_next = u + dt / 2 * (k1 + k2)
    return u_next

# Example usage
# Define A, F, K matrices and initial condition u
A = np.eye(3) # Example A matrix
F = np.array([0, 0, 0]) # Example F vector
K = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]]) # Example K matrix
A_inv = np.linalg.inv(A)

u = np.array([10.5440, 4.1124, 35.8233])
dt = 0.01
N = 10000
u_series = [u]

for i in range(N):
    u = gauss_legendre_4th_order_step(u, dynamics, dt, A_inv, F, K)
    u_series.append(u)

u_series = np.array(u_series).T

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(u_series[0], u_series[1], u_series[2])
ax.set_xticks([])
ax.set_yticks([])
ax.set_zticks([])
ax.set_title('System Dynamics')
plt.savefig('HighOrder')