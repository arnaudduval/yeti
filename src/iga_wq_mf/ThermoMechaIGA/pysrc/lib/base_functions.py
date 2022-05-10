"""
.. This module contains basis functions which use python or fortran
.. Joaquin Cornejo
"""

# Python libraries
import numpy as np
from geomdl import helpers
from scipy import sparse as sp

# My libraries
from iga_wq_mf import basis_weights

# ==========================
# B-spline functions 
# ==========================

def create_knotvector(p, nbel, multiplicity= 1):
    " Creates an uniform and open knot-vector "

    # Set knot-vector to be inserted
    knotvector_Unique = np.linspace(0., 1., nbel + 1)[1 : -1]

    # Create knot-vector 
    knotvector = []
    for _ in range(p+1): 
        knotvector.append(0.0)

    for knot in knotvector_Unique: 
        for _ in range(multiplicity): 
            knotvector.append(knot)

    for _ in range(p+1): 
        knotvector.append(1.0)
    
    return knotvector

def eval_basis_python(degree, knotvector, knots, multiplicity= 1): 
    " Evaluates B-spline functions at given knots "

    # Find number of points x
    nbx = len(knots)

    # Find number of elements 
    nbel = len(np.unique(knotvector)) - 1

    # Find number of functions 
    nbfunct = degree + multiplicity*(nbel - 1) + 1

    # Set table of functions per element 
    table_functions_element = np.zeros((nbel, degree + 2), dtype= int); 
    table_functions_element[0, 0] = degree; table_functions_element[0, 1:] = np.arange(degree + 1) 

    for _ in range(1, nbel): 
        # Set values of the table
        table_functions_element[_, :2] = table_functions_element[_-1, :2] + multiplicity
        table_functions_element[_, 2:] = table_functions_element[_, 1] + np.arange(1, degree + 1) 

    # Evaluate B0 and B1
    B0 = sp.lil_matrix((nbfunct, nbx))
    B1 = sp.lil_matrix((nbfunct, nbx))

    for _ in range(len(knots)):
        # Get knot
        knot = knots[_]    
    
        # Find knot-span
        knot_span = helpers.find_span_linear(degree, knotvector, nbfunct, knot)
        
        # Find element
        element = np.where(table_functions_element[:, 0] == knot_span)[0].tolist()
        
        # Find functions at the element
        functions_element = table_functions_element[element, 1:][0]

        # Evaluate B0 and B1 at the knot
        B0t, B1t = helpers.basis_function_ders(degree, knotvector, knot_span, knot, 1)

        # Set procedure if knot is in the knot-vector
        if knot in np.unique(knotvector)[1:-1]:             
            # Erase zeros
            B0t = B0t[:-multiplicity] 
            B1t = B1t[:-multiplicity] 

            # Erase zeros functions
            functions_element = functions_element[:-multiplicity]

        # Replace values
        B0[np.ix_(functions_element, [_])] = np.asarray(B0t).reshape((-1,1))
        B1[np.ix_(functions_element, [_])] = np.asarray(B1t).reshape((-1,1))

    return B0, B1

def eval_basis_fortran(degree, nbel, knots, multiplicity=1):
    " Evaluates B-spline functions at given knots using fortran libraries "

    B0, B1, indi, indj = basis_weights.get_basis_generalized(
                            degree, nbel, len(knots), knots, multiplicity)

    return B0, B1, indi, indj

# ==========================
# IGA - WQ FUNCTIONS
# ==========================
