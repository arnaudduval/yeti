# "In this file we will try to implement our code to a multipatch case"
# from pysrc.lib.__init__ import *
# from pysrc.lib.lib_part import part

# #IGA module
# from preprocessing.igaparametrization import IGAparametrization

# # Selection of .INP and .NB file
# modelIGA = IGAparametrization(filename='tensileSpecimen')
# nb_degDV = np.ones((3, modelIGA._nb_patch), dtype=np.intp)
# nb_refDV = np.ones((3, modelIGA._nb_patch), dtype=np.intp)
# modelIGA.refine(nb_refDV, nb_degDV)

# # Create model 
# quadArgs = {'quadrule': 'iga', 'type': 'leg'}
# model    = part(modelIGA, quadArgs=quadArgs)
