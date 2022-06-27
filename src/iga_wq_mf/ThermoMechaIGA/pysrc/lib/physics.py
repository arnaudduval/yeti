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
        f = -div(lambda * grad(u))
    """
    x = P[0, :]
    y = P[1, :]
    z = P[2, :]

    # # Isotropy
    # f = 3*pi**2*sin(pi*x)*sin(pi*y)*sin(pi*z) 

    # Anisotropy
    f = (pi**2*cos(pi*x)*cos(pi*y)*sin(pi*z) 
    + (pi**2*cos(pi*x)*cos(pi*z)*sin(pi*y))/5 
    + (pi**2*cos(pi*y)*cos(pi*z)*sin(pi*x))/2 
    + 6*pi**2*sin(pi*x)*sin(pi*y)*sin(pi*z)
    )

    return f

def powden_prism(P: list):
    """ u = (-5*x+6*y+45)*(5*x+6*y-45)*x*(x-6)*sin(pi*z)
        f = -div(lambda * grad(u))
    """
    x = P[0, :]
    y = P[1, :]
    z = P[2, :]

    # # Isotropy
    # f = (10*x*sin(pi*z)*(5*x + 6*y - 45) - 22*x*sin(pi*z)*(x - 6) 
    #     - 10*x*sin(pi*z)*(6*y - 5*x + 45) - 2*sin(pi*z)*(6*y - 5*x + 45)*(5*x + 6*y - 45) 
    #     - 10*sin(pi*z)*(x - 6)*(6*y - 5*x + 45) + 10*sin(pi*z)*(x - 6)*(5*x + 6*y - 45) 
    #     + x*pi**2*sin(pi*z)*(x - 6)*(6*y - 5*x + 45)*(5*x + 6*y - 45)
    # )

    # Anisotropy
    f = (16*x*sin(pi*z)*(5*x + 6*y - 45) 
    - 94*x*sin(pi*z)*(x - 6) 
    - 4*x*sin(pi*z)*(6*y - 5*x + 45) 
    - 2*sin(pi*z)*(6*y - 5*x + 45)*(5*x + 6*y - 45) 
    - 4*sin(pi*z)*(x - 6)*(6*y - 5*x + 45) 
    + 16*sin(pi*z)*(x - 6)*(5*x + 6*y - 45) 
    + (pi*cos(pi*z)*(x - 6)*(6*y - 5*x + 45)*(5*x + 6*y - 45))/5 
    + 4*x*pi*cos(pi*z)*(x - 6)*(6*y - 5*x + 45) 
    + 2*x*pi*cos(pi*z)*(x - 6)*(5*x + 6*y - 45) 
    + (x*pi*cos(pi*z)*(6*y - 5*x + 45)*(5*x + 6*y - 45))/5 
    + 3*x*pi**2*sin(pi*z)*(x - 6)*(6*y - 5*x + 45)*(5*x + 6*y - 45)
    )    

    return f

def powden_thickring(P: list):
    """ u = sin(5*pi*x)*sin(5*pi*y)*sin(5*pi*z)*(x**2+y**2-1)*(x**2+y**2-4)
        f = -div(lambda * grad(u))
    """
    x = P[0, :]
    y = P[1, :]
    z = P[2, :] 

    # # Isotropy
    # f = (75*pi**2*sin(5*pi*x)*sin(5*pi*y)*sin(5*pi*z)*(x**2 + y**2 - 1)*(x**2 + y**2 - 4) 
    #     - 8*y**2*sin(5*pi*x)*sin(5*pi*y)*sin(5*pi*z) 
    #     - 4*sin(5*pi*x)*sin(5*pi*y)*sin(5*pi*z)*(x**2 + y**2 - 1) 
    #     - 4*sin(5*pi*x)*sin(5*pi*y)*sin(5*pi*z)*(x**2 + y**2 - 4) 
    #     - 8*x**2*sin(5*pi*x)*sin(5*pi*y)*sin(5*pi*z) 
    #     - 20*x*pi*cos(5*pi*x)*sin(5*pi*y)*sin(5*pi*z)*(x**2 + y**2 - 1) 
    #     - 20*x*pi*cos(5*pi*x)*sin(5*pi*y)*sin(5*pi*z)*(x**2 + y**2 - 4) 
    #     - 20*y*pi*cos(5*pi*y)*sin(5*pi*x)*sin(5*pi*z)*(x**2 + y**2 - 1) 
    #     - 20*y*pi*cos(5*pi*y)*sin(5*pi*x)*sin(5*pi*z)*(x**2 + y**2 - 4)
    # )

    # Anisotropy
    f = (8*x*y*sin(5*pi*x)*sin(5*pi*y)*sin(5*pi*z) 
    - 16*y**2*sin(5*pi*x)*sin(5*pi*y)*sin(5*pi*z) 
    - 6*sin(5*pi*x)*sin(5*pi*y)*sin(5*pi*z)*(x**2 + y**2 - 1) 
    - 6*sin(5*pi*x)*sin(5*pi*y)*sin(5*pi*z)*(x**2 + y**2 - 4) 
    - 8*x**2*sin(5*pi*x)*sin(5*pi*y)*sin(5*pi*z) 
    + 150*pi**2*sin(5*pi*x)*sin(5*pi*y)*sin(5*pi*z)*(x**2 + y**2 - 1)*(x**2 + y**2 - 4) 
    + 25*pi**2*cos(5*pi*x)*cos(5*pi*y)*sin(5*pi*z)*(x**2 + y**2 - 1)*(x**2 + y**2 - 4) 
    + 5*pi**2*cos(5*pi*x)*cos(5*pi*z)*sin(5*pi*y)*(x**2 + y**2 - 1)*(x**2 + y**2 - 4) 
    + (25*pi**2*cos(5*pi*y)*cos(5*pi*z)*sin(5*pi*x)*(x**2 + y**2 - 1)*(x**2 + y**2 - 4))/2 
    - 20*x*pi*cos(5*pi*x)*sin(5*pi*y)*sin(5*pi*z)*(x**2 + y**2 - 1) 
    + 10*x*pi*cos(5*pi*y)*sin(5*pi*x)*sin(5*pi*z)*(x**2 + y**2 - 1) 
    + 2*x*pi*cos(5*pi*z)*sin(5*pi*x)*sin(5*pi*y)*(x**2 + y**2 - 1) 
    - 20*x*pi*cos(5*pi*x)*sin(5*pi*y)*sin(5*pi*z)*(x**2 + y**2 - 4) 
    + 10*x*pi*cos(5*pi*y)*sin(5*pi*x)*sin(5*pi*z)*(x**2 + y**2 - 4) 
    + 2*x*pi*cos(5*pi*z)*sin(5*pi*x)*sin(5*pi*y)*(x**2 + y**2 - 4) 
    + 10*y*pi*cos(5*pi*x)*sin(5*pi*y)*sin(5*pi*z)*(x**2 + y**2 - 1) 
    - 40*y*pi*cos(5*pi*y)*sin(5*pi*x)*sin(5*pi*z)*(x**2 + y**2 - 1) 
    + 5*y*pi*cos(5*pi*z)*sin(5*pi*x)*sin(5*pi*y)*(x**2 + y**2 - 1) 
    + 10*y*pi*cos(5*pi*x)*sin(5*pi*y)*sin(5*pi*z)*(x**2 + y**2 - 4) 
    - 40*y*pi*cos(5*pi*y)*sin(5*pi*x)*sin(5*pi*z)*(x**2 + y**2 - 4) 
    + 5*y*pi*cos(5*pi*z)*sin(5*pi*x)*sin(5*pi*y)*(x**2 + y**2 - 4)
    )

    return f

def powden_rotring(P: list):
    """ u = -(x**2 + y**2 - 1)*(x**2 + y**2 - 4)*x*(y**2)*sin(pi*z)
        f = -div(lambda * grad(u))
    """
    x = P[0, :]
    y = P[1, :]
    z = P[2, :]

    # Isotropy
    f = (8*x*y**4*np.sin(np.pi*z) + 8*x**3*y**2*np.sin(np.pi*z) 
        + 16*x*y**2*np.sin(np.pi*z)*(x**2 + y**2 - 1) + 16*x*y**2*np.sin(np.pi*z)*(x**2 + y**2 - 4) 
        + 2*x*np.sin(np.pi*z)*(x**2 + y**2 - 1)*(x**2 + y**2 - 4) 
        - x*y**2*np.pi**2*np.sin(np.pi*z)*(x**2 + y**2 - 1)*(x**2 + y**2 - 4)
    )

    return f

def temperature_rotring(P: list):
    " T = -(x**2 + y**2 - 1)*(x**2 + y**2 - 4)*x*(y**2)*sin(pi*z) "
    x = P[0, :]
    y = P[1, :]
    z = P[2, :]
    f = -(x**2 + y**2 - 1)*(x**2 + y**2 - 4)*x*(y**2)*np.sin(np.pi*z)

    return f

def powden_annulus(P: list):
    """ u = (x**2 + y**2 - 1)*(x**2 + y**2 - 4)*sin(pi*x)*sin(pi*y)
        f = -div(lambda * grad(u))
    """
    x = P[0, :]
    y = P[1, :]

    # # Isotropy
    # f = (2*pi**2*sin(pi*x)*sin(pi*y)*(x**2 + y**2 - 1)*(x**2 + y**2 - 4) 
    #     - 8*y**2*sin(pi*x)*sin(pi*y) 
    #     - 4*sin(pi*x)*sin(pi*y)*(x**2 + y**2 - 1) 
    #     - 4*sin(pi*x)*sin(pi*y)*(x**2 + y**2 - 4) 
    #     - 4*x*pi*cos(pi*x)*sin(pi*y)*(x**2 + y**2 - 1) 
    #     - 4*x*pi*cos(pi*x)*sin(pi*y)*(x**2 + y**2 - 4) 
    #     - 4*y*pi*cos(pi*y)*sin(pi*x)*(x**2 + y**2 - 1) 
    #     - 4*y*pi*cos(pi*y)*sin(pi*x)*(x**2 + y**2 - 4) 
    #     - 8*x**2*sin(pi*x)*sin(pi*y)
    # )

    # Anisotropy lambda = [1, 0; 0, 0.1]
    f = ((11*pi**2*sin(pi*x)*sin(pi*y)*(x**2 + y**2 - 1)*(x**2 + y**2 - 4))/10 
        - (4*y**2*sin(pi*x)*sin(pi*y))/5 
        - (11*sin(pi*x)*sin(pi*y)*(x**2 + y**2 - 1))/5 
        - (11*sin(pi*x)*sin(pi*y)*(x**2 + y**2 - 4))/5 
        - 4*x*pi*cos(pi*x)*sin(pi*y)*(x**2 + y**2 - 1) 
        - 4*x*pi*cos(pi*x)*sin(pi*y)*(x**2 + y**2 - 4) 
        - (2*y*pi*cos(pi*y)*sin(pi*x)*(x**2 + y**2 - 1))/5 
        - (2*y*pi*cos(pi*y)*sin(pi*x)*(x**2 + y**2 - 4))/5 
        - 8*x**2*sin(pi*x)*sin(pi*y)
    )
    
    return f