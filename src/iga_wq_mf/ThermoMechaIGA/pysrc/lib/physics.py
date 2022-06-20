"""
.. module :: Physics - Heat transfer
    :synopsis: Provides functions used in heat transfer equation
.. author :: Joaquin Cornejo
"""
import numpy as np
from numpy import sin, cos, pi

def power_density(P: list):
    " Compute power density at point P in physical space"

    f = 1 

    return f

def powden_cube(P: list):
    """ u = sin(pi*x)*sin(pi*y)*sin(pi*z)
        f = -(d2u/dx2 + d2u/dy2 + d2u/dz2)
    """
    x = P[0]
    y = P[1]
    z = P[2]

    f = 3*pi**2*sin(pi*x)*sin(pi*y)*sin(pi*z)

    return f

def powden_prism(P: list):
    """ u = (-5*x+6*y+45)*(5*x+6*y-45)*x*(x-6)*sin(pi*z)
        f = -(d2u/dx2 + d2u/dy2 + d2u/dz2)
    """
    x = P[0]
    y = P[1]
    z = P[2]

    f = (10*x*sin(pi*z)*(5*x + 6*y - 45) - 22*x*sin(pi*z)*(x - 6) 
        - 10*x*sin(pi*z)*(6*y - 5*x + 45) - 2*sin(pi*z)*(6*y - 5*x + 45)*(5*x + 6*y - 45) 
        - 10*sin(pi*z)*(x - 6)*(6*y - 5*x + 45) + 10*sin(pi*z)*(x - 6)*(5*x + 6*y - 45) 
        + x*pi**2*sin(pi*z)*(x - 6)*(6*y - 5*x + 45)*(5*x + 6*y - 45)
    )
    
    return f

def powden_thickring(P: list):
    """ u = sin(5*pi*x)*sin(5*pi*y)*sin(5*pi*z)*(x**2+y**2-1)*(x**2+y**2-4)
        f = -(d2u/dx2 + d2u/dy2 + d2u/dz2)
    """
    x = P[0]
    y = P[1]
    z = P[2] 

    f = (75*pi**2*sin(5*pi*x)*sin(5*pi*y)*sin(5*pi*z)*(x**2 + y**2 - 1)*(x**2 + y**2 - 4) 
        - 8*y**2*sin(5*pi*x)*sin(5*pi*y)*sin(5*pi*z) - 4*sin(5*pi*x)*sin(5*pi*y)*sin(5*pi*z)*(x**2 + y**2 - 1) 
        - 4*sin(5*pi*x)*sin(5*pi*y)*sin(5*pi*z)*(x**2 + y**2 - 4) - 8*x**2*sin(5*pi*x)*sin(5*pi*y)*sin(5*pi*z) 
        - 20*x*pi*cos(5*pi*x)*sin(5*pi*y)*sin(5*pi*z)*(x**2 + y**2 - 1) 
        - 20*x*pi*cos(5*pi*x)*sin(5*pi*y)*sin(5*pi*z)*(x**2 + y**2 - 4) 
        - 20*y*pi*cos(5*pi*y)*sin(5*pi*x)*sin(5*pi*z)*(x**2 + y**2 - 1) 
        - 20*y*pi*cos(5*pi*y)*sin(5*pi*x)*sin(5*pi*z)*(x**2 + y**2 - 4)
    )
    
    return f

def powden_rotring(P: list):
    """
    u = -(x**2 + y**2 - 1)*(x**2 + y**2 - 4)*x*(y**2)*sin(pi*z)
    f = -(d2u/dx2 + d2u/dy2 + d2u/dz2)
    """
    x = P[0]
    y = P[1]
    z = P[2]
    f = (8*x*y**4*np.sin(np.pi*z) + 8*x**3*y**2*np.sin(np.pi*z) 
        + 16*x*y**2*np.sin(np.pi*z)*(x**2 + y**2 - 1) + 16*x*y**2*np.sin(np.pi*z)*(x**2 + y**2 - 4) 
        + 2*x*np.sin(np.pi*z)*(x**2 + y**2 - 1)*(x**2 + y**2 - 4) 
        - x*y**2*np.pi**2*np.sin(np.pi*z)*(x**2 + y**2 - 1)*(x**2 + y**2 - 4)
    )
    
    return f

def temperature_rotring(P):
    "T = -(x**2 + y**2 - 1)*(x**2 + y**2 - 4)*x*(y**2)*sin(pi*z)"
    x = P[0]
    y = P[1]
    z = P[2]
    f = -(x**2 + y**2 - 1)*(x**2 + y**2 - 4)*x*(y**2)*np.sin(np.pi*z)

    return f