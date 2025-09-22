"""
utility tools for stiffened structure optimization
"""

import numpy as np

from stiffened_cylinder import RADIUS, HEIGHT, NB_STIFFENERS_BY_HALF

def smoothmin(x, p=10):
    """
    Compute the min value of an array with a p-norm in order to be differenciable
    Gradient is also computed
    """

    x = np.array(x)
    assert np.all(x > 0)
    S = np.sum(x ** (-p))
    value = S ** (-1 / p)
    grad = -x ** (-p - 1) / (S ** (1 + 1/p))
    return value, grad

def len_stiffener(x):
    """
    Compute the length of a stiffener given its coordinatesdesign variables
    Design variables are given as a size 4 array of design variables that contains
    parametric coordinates u1, u2, v1, v2
    """
    circ = RADIUS*np.pi
    return np.sqrt((circ*(x[1]-x[0]))**2. + (HEIGHT*(x[3]-x[2]))**2.)

def grad_len_stiffener(x):
    """
    Gradient of the length of one stiffener
    """
    circ = RADIUS*np.pi
    A = circ * (x[1] - x[0])
    B = HEIGHT * (x[3] - x[2])
    denom = np.sqrt(A**2 + B**2)
    grad = np.array([
        -circ * A / denom,
        +circ * A / denom,
        -HEIGHT * B / denom,
        +HEIGHT * B / denom
    ])
    return grad

def len_stiffeners(xk):
    """
    Compute length of stiffeners
    Parameters ar a full set of design variables, grouped by 4 for a given stiffener
    """
    lengths = np.empty(NB_STIFFENERS_BY_HALF)
    for istiff in range(NB_STIFFENERS_BY_HALF):
        lengths[istiff] = len_stiffener(xk[4*istiff:4*istiff+4])
    return lengths

def grad_len_stiffeners(xk):
    """
    Compute gradient of length of stiffeners
    """
    grad = np.empty((NB_STIFFENERS_BY_HALF, len(xk)))

    for istiff in range(NB_STIFFENERS_BY_HALF):
        i0 = 4*istiff
        x = xk[i0:i0+4]
        grad[istiff, i0:i0+4] = grad_len_stiffener(x)

    return grad